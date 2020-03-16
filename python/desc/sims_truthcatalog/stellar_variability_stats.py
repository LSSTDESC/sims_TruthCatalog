"""
Module to compute statistics of stellar light curves using lsst_sims
VariabilityMixin code.
"""
import os
import shutil
import json
import sqlite3
import numpy as np
import pandas as pd
from lsst.sims.utils import defaultSpecMap, angularSeparation
import lsst.sims.photUtils as sims_photUtils
import lsst.sims.catUtils.mixins.VariabilityMixin as variability


__all__ = ['write_star_variability_stats', 'merge_sqlite3_dbs',
           'StellarLightCurveFactory', 'write_star_variability_truth']


class VariabilityGenerator(variability.StellarVariabilityModels,
                           variability.MLTflaringMixin,
                           variability.ParametrizedLightCurveMixin):
    """
    This is a class that mimics the API of an InstanceCatalog. This makes it
    easier for us to generate light curves from sources with differing
    variability models.

    This class was copied from
    https://github.com/LSSTDESC/sims_GCRCatSimInterface/blob/master/workspace/truth/generate_stellar_light_curves.py#L25
    """

    def __init__(self, chunk):
        """
        chunk = (simobjid, hpid, sedFilename, magNorm, ebv, varParamStr,
                 parallax, ra, dec)
        """
        bands = 'ugrizy'
        # initialize some member variables expected by the
        # stellar variability models
        self.photParams = sims_photUtils.PhotometricParameters(nexp=1,
                                                               exptime=30)
        self.lsstBandpassDict \
            = sims_photUtils.BandpassDict.loadTotalBandpassesFromFiles()
        self._actually_calculated_columns = []
        for bp in bands:
            self._actually_calculated_columns.append('lsst_%s' % bp)

        sed_dir = os.environ['SIMS_SED_LIBRARY_DIR']
        self.quiescent_mags = {}
        for bp in bands:
            self.quiescent_mags[bp] = np.zeros(len(chunk), dtype=float)

        self.parallax = np.zeros(len(chunk), dtype=float)
        self.ebv = np.zeros(len(chunk), dtype=float)
        self.varParamStr = np.empty(len(chunk), dtype=(str, 300))
        self.simobjid = np.empty(len(chunk), dtype=int)
        ccm_wav = None

        # find the quiescent magnitudes of the stars
        for i_star, star in enumerate(chunk):
            self.simobjid[i_star] = int(star[0])
            full_file_name = os.path.join(sed_dir, defaultSpecMap[star[2]])
            spec = sims_photUtils.Sed()
            spec.readSED_flambda(full_file_name)
            fnorm = sims_photUtils.getImsimFluxNorm(spec, float(star[3]))
            spec.multiplyFluxNorm(fnorm)
            if ccm_wav is None or not np.array_equal(spec.wavelen, ccm_wav):
                ccm_wav = np.copy(spec.wavelen)
                a_x, b_x = spec.setupCCM_ab()
            ebv_val = float(star[4])
            spec.addDust(a_x, b_x, ebv=ebv_val, R_v=3.1)
            mag_list = self.lsstBandpassDict.magListForSed(spec)
            for i_bp, bp in enumerate(bands):
                self.quiescent_mags[bp][i_star] = mag_list[i_bp]
            self.ebv[i_star] = ebv_val
            self.varParamStr[i_star] = star[5]
            self.parallax[i_star] = float(star[6])

        self.parallax *= (np.pi/648000000.0) # convert from mas to radians

    def column_by_name(self, name):
        """Columns expected by lsst_sims stellar variability code."""
        if name.startswith('quiescent'):
            bp = name[-1]
            return np.copy(self.quiescent_mags[bp])
        if name == 'ebv':
            return np.copy(self.ebv)
        if name == 'parallax':
            return np.copy(self.parallax)
        if name == 'simobjid':
            return np.copy(self.simobjid)
        raise RuntimeError("\n\nCannot get column %s\n\n" % name)


def get_variability_model(varParamStr):
    """
    Extract the variability model from the varParamStr column.
    """
    if varParamStr == 'None':
        return 'None'
    params = json.loads(varParamStr)
    try:
        return params['m']
    except KeyError:
        return params['varMethodName']


def write_star_variability_stats(stars_db_file, outfile, row_min, row_max,
                                 chunk_size=10000, mjds=None):
    """
    Write a sqlite3 db file with the mean and stdev of the delta
    magnitudes for each star in the `ugrizy` bands.

    Parameters
    ----------
    stars_db_file: str
        Name of the sqlite3 file containing the `stars` db table
        containing the stellar model parameters.
    outfile: str
        Output filename.
    row_min: int
        First row to process in the stars table in stars_db_file.
    row_max: int
        Last row to process in the stars table in stars_db_file.
    chunk_size: int [10000]
        Number of rows to fetch in each chunk read from stars_db_file.
    mjds: list-like [None]
        Sequence of MJDs at which to evalute the delta mag values for
        computing the statistics.  If None, then generate 1000 random
        values between MJD 59580 and MJD 59580 + 5*365.
    """
    # Make sure that the Kepler light curves are loaded into memory
    variability.ParametrizedLightCurveMixin()\
               .load_parametrized_light_curves()

    if mjds is None:
        t0 = 59580.
        mjds = np.sort(np.random.uniform(t0, t0 + 365*5, size=1000))

    # Set up the sqlite3 cursor for the star db.
    star_db = sqlite3.connect(stars_db_file)
    query = f"""SELECT simobjid, -1, sedFilename, magNorm, ebv,
                varParamStr, parallax, ra, decl FROM stars
                limit {row_min}, {row_max - row_min}"""
    star_curs = star_db.execute(query)

    # Create the stellar_variability table.
    table_name = 'stellar_variability'
    cmd = f'''create table if not exists {table_name}
              (id TEXT, model TEXT, mean_u, mean_g, mean_r,
               mean_i, mean_z, mean_y, stdev_u, stdev_g,
               stdev_r, stdev_i, stdev_z, stdev_y)'''
    with sqlite3.connect(outfile) as conn:
        cursor = conn.cursor()
        cursor.execute(cmd)
        conn.commit()

        # Loop over chunks and write each processed chunk to the output table.
        chunk = star_curs.fetchmany(chunk_size)
        num_rows = 0
        while chunk:
            print(num_rows)
            var_gen = VariabilityGenerator(chunk)
            dmags = var_gen.applyVariability(var_gen.varParamStr, expmjd=mjds)
            # Transpose the columns so that the ordering is star, filter, mjd
            dmags = dmags.transpose([1, 0, 2])
            if len(chunk) != len(dmags):
                raise RuntimeError('# of rows in chunk do not match # of lcs')
            values = []
            for istar, lcs in enumerate(dmags):
                model = get_variability_model(var_gen.varParamStr[istar])
                row = [str(chunk[istar][0]), model]
                nbands = len(lcs)
                row.extend([np.mean(lcs[iband]) for iband in range(nbands)])
                row.extend([np.std(lcs[iband]) for iband in range(nbands)])
                values.append(row)
            cursor.executemany(f'''insert into {table_name} values
                                   (?, ?, ?, ?, ?,
                                    ?, ?, ?, ?, ?,
                                    ?, ?, ?, ?)''',
                               values)
            conn.commit()
            num_rows += len(chunk)
            chunk = star_curs.fetchmany(chunk_size)
    star_db.close()


def merge_sqlite3_dbs(infiles, outfile):
    """
    Function to merge the db tables in a list of input sqlite3 files.
    The corresponding tables in each file are assumed to have the same
    schemas.

    Parameters
    ----------
    infiles: list-like
        List of sqlite3 filenames.
    outfile: str
        Filename of output sqlite3 file to contain the merged tables.
    """
    shutil.copy(infiles[0], outfile)
    with sqlite3.connect(outfile) as conn:
        for infile in infiles[1:]:
            conn.execute(f"attach '{infile}' as in_db")
            for row in conn.execute("""select * from in_db.sqlite_master
                                       where type='table'"""):
                combine = f"insert into {row[1]} select * from in_db.{row[1]}"
                conn.execute(combine)
            conn.commit()
            conn.execute("detach database in_db")


class StellarLightCurveFactory:
    """
    Factory class to generate stellar light curves using model parameters
    from a stars db file.
    """
    def __init__(self, stars_db_file, mjds=None):
        """
        Parameters
        ----------
        stars_db_file: str
            Name of the sqlite3 file containing the `stars` db table
            containing the stellar model parameters.
        mjds: list-like [None]
            Sequence of MJDs at which to evaluate the delta mag values.
            If None, then np.arange(0, 365*5, 5 + 59580 will be used.
            This sequence can be overridden in the `.create` method.
        """
        self.conn = sqlite3.connect(stars_db_file)
        # Make sure that the Kepler light curves are loaded into memory
        variability.ParametrizedLightCurveMixin()\
                   .load_parametrized_light_curves()
        t0 = 59580.
        self.mjds = np.arange(0, 365*5, 5) + t0 if mjds is None else mjds

    def __del__(self):
        self.conn.close()

    def create(self, simobjid, mjds=None):
        """
        Parameters
        ----------
        simobjid: int
            ID number of the star in the stars db table.
        mjds: list-like [None]
            Sequence of MJDs at which to evaluate the delta mag values.

        Results
        -------
        dict of np.arrays of delta mag values keyed by `ugrizy` band.
        """
        if mjds is None:
            mjds = self.mjds
        bands = 'ugrizy'
        query = f"""SELECT simobjid, -1, sedFilename, magNorm, ebv,
                    varParamStr, parallax, ra, decl FROM stars where
                    simobjid={simobjid} limit 1"""
        chunk = list(self.conn.execute(query))
        var_gen = VariabilityGenerator(chunk)
        dmags = var_gen.applyVariability(var_gen.varParamStr, expmjd=mjds)
        # Transpose the columns so that the ordering is star, filter, mjd
        dmags = dmags.transpose([1, 0, 2])

        return dict(zip(bands, dmags[0])), var_gen.quiescent_mags


def write_star_variability_truth(outfile, stars_db_file, opsim_db_file,
                                 star_lc_stats_db_file,
                                 fp_radius=2.05, max_num_stars=None,
                                 dmag_threshold=1e-3):
    """
    Write the variability truth table for stars that show standard
    deviations in delta magnitudes > dmag_threshold.

    Parameters
    ----------
    outfile: str
        Output sqlite3 file to contain the variability truth table.
    star_db_file: str
        Sqlite3 db file containing the information on the properties
        of each star.
    opsim_db_file: str
        The sqlite3 file containing the OpSim Summary table which
        has the pointing information for each visit.
    star_lc_stats_db_file: str
        Sqlite3 db file containing summary statistics for each star
        in the star_db_file.
    fp_radius: float [2.05]
        Effective radius of the focal plane in degrees.  This defines
        the acceptance cone centered on the pointing direction for
        determining if an object is being observed by LSST for the
        purpose of computing a flux entry for the visit to be entered
        in the Variability Truth Table.
    max_num_stars: int [None]
        Maximum number of stars to process from the stars_db_file.
    dmag_threshold: float [1e-3]
        Threshold in magnitudes for computing light curve data for
        a given star.  If the stdev of the delta magnitude values in
        any band is above threshold, then the light curvve data are
        computed.
    """
    bands = 'ugrizy'

    # Read in the needed columns from the opsim_db file and add
    # ra, dec columns in degrees to enable per-object visit selections.
    with sqlite3.connect(opsim_db_file) as conn:
        opsim_db = pd.read_sql('''select obsHistID, descDitheredRA,
                                  descDitheredDec, filter, expMJD from Summary
                                  order by expMJD asc''',
                               conn)
    opsim_db['ra'] = np.degrees(opsim_db['descDitheredRA'])
    opsim_db['dec'] = np.degrees(opsim_db['descDitheredDec'])

    # Find the stars with stdev(delta_mag) > dmag_threshold in any band.
    query = ('select id, model from stellar_variability where '
             + ' or '.join([f'stdev_{band} > {dmag_threshold}'
                            for band in bands]))
    with sqlite3.connect(star_lc_stats_db_file) as conn:
        variable_stars = dict(conn.execute(query))

    # DC2 Run 2 boundaries
    ra_bounds = (49.92, 73.79)      # bounds at dec=-44.33
    dec_bounds = (-44.33, -27.25)

    # Create variabilty truth table in the output sqlite3 file.
    table_name = 'stellar_variability_truth'
    cmd = f'''create table if not exists {table_name}
              (id TEXT, obsHistID INTEGER, MJD FLOAT, bandpass TEXT,
               delta_flux FLOAT)'''
    output = sqlite3.connect(outfile)
    cursor = output.cursor()
    cursor.execute(cmd)
    output.commit()

    # Loop over stars in the stars_db_file and compute flux points
    # for visits in which they are observed.
    slc_factory = StellarLightCurveFactory(stars_db_file=stars_db_file)
    num_stars = 0
    with sqlite3.connect(stars_db_file) as conn:
        for star_id in variable_stars:
            # Get the RA, Dec of the current star.
            query = f'''select ra, decl from stars where
                        simobjid={star_id} limit 1'''
            ra, dec = list(conn.execute(query))[0]

            # Rough cut to skip stars not within the Run 2 boundaries.
            if not (ra_bounds[0] <= ra <= ra_bounds[1] and
                    dec_bounds[0] <= dec <= dec_bounds[1]):
                continue
            num_stars += 1
            if num_stars > max_num_stars:
                break

            # Find visits with pointing directions within fp_radius of
            # the star, selecting in ra and dec ranges first, then
            # calculating and selecting on the angular separations.
            dec_min, dec_max = dec - fp_radius, dec + fp_radius
            cos_dec = np.cos(np.radians(dec))
            ra_min, ra_max = (ra - fp_radius/cos_dec, ra + fp_radius/cos_dec)
            query = (f'{dec_min} <= dec <= {dec_max} and '
                     + f'{ra_min} <= ra <= {ra_max}')
            df = pd.DataFrame(opsim_db.query(query))
            df['ang_sep'] = angularSeparation(
                df['ra'].to_numpy(), df['dec'].to_numpy(), ra, dec)
            visits = df.query(f'ang_sep < {fp_radius}')
            if len(visits) == 0:
                # No visits found for this star, so proceed to next one.
                continue

            # Visit mjds and bandpasses:
            obsHistIDs = visits['obsHistID'].to_numpy()
            mjds = visits['expMJD'].to_numpy()
            bps = visits['filter'].to_numpy()

            # Compute the delta flux values in nano-Janskys.
            dm, m0 = slc_factory.create(star_id, mjds)
            dfluxes = dict()
            for band in bands:
                dfluxes[band] = (10.**((8.9 - (m0[band] + dm[band]))/2.5)
                                 - 10.**((8.9 - m0[band])/2.5))*1e9

            # Write the rows to the output table.
            values = []
            for i, band in enumerate(bps):
                values.append((star_id, int(obsHistIDs[i]), mjds[i], band,
                               dfluxes[band][i]))
            cursor.executemany(f'''INSERT INTO {table_name} VALUES
                                   (?, ?, ?, ?, ?)''', values)
            output.commit()
        output.close()
