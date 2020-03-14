import os
import sys
import json
import sqlite3
import multiprocessing
import numpy as np
import pandas as pd
from lsst.sims.utils import defaultSpecMap
import lsst.sims.photUtils as sims_photUtils
import lsst.sims.catUtils.mixins.VariabilityMixin as variability


__all__ = ['write_star_variability_stats', 'merge_sqlite3_dbs']


class VariabilityGenerator(variability.StellarVariabilityModels,
                           variability.MLTflaringMixin,
                           variability.ParametrizedLightCurveMixin):
    """
    This is a class that mimics the API of an InstanceCatalog. This makes it
    easier for us to generate light curves from sources with differing
    variability models.
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
        self.varParamStr = np.empty(len(chunk), dtype=(str,300))
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
        if name.startswith('quiescent'):
            bp = name[-1]
            return np.copy(self.quiescent_mags[bp])
        elif name == 'ebv':
            return np.copy(self.ebv)
        elif name=='parallax':
            return np.copy(self.parallax)
        elif name=='simobjid':
            return np.copy(self.simobjid)
        raise RuntimeError("\n\nCannot get column %s\n\n" % name)


def get_variability_model(varParamStr):
    if varParamStr == 'None':
        return 'None'
    params = json.loads(varParamStr)
    try:
        return params['m']
    except KeyError:
        return params['varMethodName']


def write_star_variability_stats(stars_db_file, outfile, row_min, row_max,
                                 chunk_size=10000):
    # Make sure that the Kepler light curves are loaded into memory
    variability.ParametrizedLightCurveMixin()\
               .load_parametrized_light_curves()

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
            for istar in range(len(dmags)):
                nbands = len(dmags[istar])
                model = get_variability_model(var_gen.varParamStr[istar])
                row = [str(chunk[istar][0]), model]
                row.extend([np.mean(dmags[istar][iband]) for iband
                            in range(nbands)])
                row.extend([np.std(dmags[istar][iband]) for iband
                            in range(nbands)])
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
