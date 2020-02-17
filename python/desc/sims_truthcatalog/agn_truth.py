"""
Module to write truth catalogs for AGNs using the AGNs parameters db as input.
"""
import os
import math
from collections import defaultdict
import json
import sqlite3
import numpy as np
import pandas as pd
from lsst.sims.utils import angularSeparation
from .synthetic_photometry import SyntheticPhotometry, find_sed_file
from .write_sqlite import write_sqlite


__all__ = ['AGNTruthWriter', 'agn_mag_norms']


def agn_mag_norms(mjds, redshift, tau, sf, seed):
    """
    Return the delta mag_norm values wrt the infinite-time average
    mag_norm for the provided AGN light curve parameters.  mag_norm is
    the object's un-reddened monochromatic magnitude at 500nm.

    Parameters
    ----------
    mjds: np.array
        Times at which to evaluate the light curve delta flux values in MJD.
        Observer frame.
    redshift: float
        Redshift of the AGN, used to account for time-dilation between
        rest frame and observer frame.
    tau: float
        Variability time scale in days.
    sf: float
        Structure function parameter, i.e., asymptotic rms variability on
        long time scales.
    seed: int
        Random number seed.

    Returns
    -------
    np.array of delta mag_norm values.

    Notes
    -----
    This code is stolen from
    https://github.com/astroML/astroML/blob/master/astroML/time_series/generate.py
    """
    rng = np.random.RandomState(seed)
    t_rest = mjds/(1 + redshift)/tau
    nbins = len(t_rest)
    steps = rng.normal(0, 1, nbins)
    delta_mag_norm = np.zeros(nbins)
    delta_mag_norm[0] = steps[0]*sf
    for i in range(1, nbins):
        dt = t_rest[i] - t_rest[i - 1]
        delta_mag_norm[i] \
            = delta_mag_norm[i - 1]*(1. - dt) + np.sqrt(2*dt)*sf*steps[i]
    return delta_mag_norm


class AGNTruthWriter:
    '''
    Write Summary and Variable truth tables for unlensed AGNs.
    '''
    agn_type_id = 117
    def __init__(self, outfile, agn_db_file):
        '''
        Parameters
        ----------
        outfile: str
            Name of the sqlite3 file to contain the truth tables.
        agn_db_file: str
            The sqlite3 file containing the AGN model parameters.
        '''
        self.outfile = outfile
        if os.path.isfile(outfile):
            raise OSError(f'{outfile} already exists.')
        if not os.path.isfile(agn_db_file):
            raise FileNotFoundError(f'{agn_db_file} not found.')
        self.conn = sqlite3.connect(agn_db_file)
        self.query = '''select galaxy_id, magNorm, redshift, M_i, ra, dec,
                     varParamStr from agn_params'''
        curs = self.conn.execute(self.query)
        self.icol = {_[0]: icol for icol, _ in enumerate(curs.description)}

    @staticmethod
    def object_id(galaxy_id):
        """Return the AGN object ID based on the host galaxy ID"""
        return str(galaxy_id*1024 + AGNTruthWriter.agn_type_id)

    def write(self, chunk_size=10000, verbose=False):
        '''
        Extract the column data from the agn db file and write
        the summary truth table to the sqlite file.

        Parameters
        ----------
        chunk_size: int [10000]
            Number of records to read in at a time from the star db
            file and write to the output file.
        verbose: bool [False]
            Flag to write the number of records that have been processed.
        '''
        bands = 'ugrizy'
        curs = self.conn.execute(self.query)
        irec = 0
        while True:
            ids, galaxy_ids, ra, dec, redshift = [], [], [], [], []
            is_variable, is_pointsource, good_ixes = [], [], []
            flux_by_band_MW = {_: [] for _ in bands}
            flux_by_band_noMW = {_: [] for _ in bands}
            chunk = curs.fetchmany(chunk_size)
            if not chunk:
                break
            for row in chunk:
                if verbose:
                    print(irec)
                irec += 1
                # All AGNs are variable point sources:
                is_pointsource.append(1)
                is_variable.append(1)
                # AGN-dependent entries:
                ra.append(row[self.icol['ra']])
                dec.append(row[self.icol['dec']])
                redshift.append(row[self.icol['redshift']])
                ids.append(self.object_id(row[self.icol['galaxy_id']]))
                galaxy_ids.append(row[self.icol['galaxy_id']])
                sed_file = find_sed_file('agnSED/agn.spec.gz')

                # Create SyntheticPhotometry object initially without
                # Milky Way dust parameters.
                synth_phot = SyntheticPhotometry(sed_file,
                                                 row[self.icol['magNorm']],
                                                 redshift[-1])
                for band in bands:
                    flux_by_band_noMW[band].append(synth_phot.calcFlux(band))

                # Set Milky Way dust parameters and compute ugrizy fluxes.
                synth_phot.add_MW_dust(ra[-1], dec[-1], Rv=3.1)
                for band in bands:
                    flux_by_band_MW[band].append(synth_phot.calcFlux(band))
            write_sqlite(self.outfile,
                         ids=ids,
                         galaxy_ids=galaxy_ids,
                         ra=ra,
                         dec=dec,
                         redshift=redshift,
                         is_variable=is_variable,
                         is_pointsource=is_pointsource,
                         flux_by_band_MW=flux_by_band_MW,
                         flux_by_band_noMW=flux_by_band_noMW,
                         good_ixes=range(len(ids)))

    def write_auxiliary_truth(self, chunk_size=10000, verbose=False):
        """
        Write the AGN auxiliary truth table from the AGN db file.

        Parameters
        ----------
        chunk_size: int [10000]
            Number of records to read in at a time from the star db
            file and write to the output file.
        verbose: bool [False]
            Flag to write the number of records that have been processed.
        """
        bands = 'ugrizy'
        curs = self.conn.execute(self.query)

        table_name = 'agn_auxiliary_info'
        cmd = f'''CREATE TABLE IF NOT EXISTS {table_name}
                  (id TEXT, host_galaxy BIGINT, M_i DOUBLE, seed BIGINT,
                  tau_u DOUBLE, tau_g DOUBLE, tau_r DOUBLE,
                  tau_i DOUBLE, tau_z DOUBLE, tau_y DOUBLE,
                  sf_u DOUBLE, sf_g DOUBLE, sf_r DOUBLE,
                  sf_i DOUBLE, sf_z DOUBLE, sf_y DOUBLE)'''
        with sqlite3.connect(self.outfile) as conn:
            cursor = conn.cursor()
            cursor.execute(cmd)
            conn.commit()

            irec = 0
            while True:
                ids, host_galaxies, M_is, seeds = [], [], [], []
                taus = {_: [] for _ in bands}
                sfs = {_: [] for _ in bands}
                chunk = curs.fetchmany(chunk_size)
                if not chunk:
                    break
                values = []
                for row in chunk:
                    if verbose:
                        print(irec)
                    irec += 1
                    pars = json.loads(row[self.icol['varParamStr']])['p']
                    my_row = [self.object_id(row[self.icol['galaxy_id']]),
                              row[self.icol['galaxy_id']],
                              row[self.icol['M_i']],
                              pars['seed']]
                    my_row.extend([pars[f'agn_tau_{band}'] for band in bands])
                    my_row.extend([pars[f'agn_sf_{band}'] for band in bands])
                    values.append(my_row)
                cursor.executemany(f'''INSERT INTO {table_name} VALUES
                                       (?, ?, ?, ?,
                                        ?, ?, ?, ?, ?, ?,
                                        ?, ?, ?, ?, ?, ?)''', values)
                conn.commit()

    def write_variability_truth(self, opsim_db_file, fp_radius=2.05,
                                max_rows=None, verbose=False):
        """
        Write the AGN fluxes for each visit.

        Parameters
        ----------
        opsim_db_file: str
            The sqlite3 file containing the OpSim Summary table which
            has the pointing information for each visit.
        fp_radius: float [2.05]
            Effective radius of the focal plane in degrees.  This defines
            the acceptance cone centered on the pointing direction for
            determining if an object is being observed by LSST for the
            purpose of computing a flux entry for the visit to be entered
            in the Variability Truth Table.
        max_rows: int [None]
            Threshold number of rows to write to the table.  This is useful
            for testing.  If None, then write all entries for all SNe in
            the sne_db_file.
        """
        bands = 'ugrizy'

        # Retrieve the pointing information for each visit from the opsim db.
        with sqlite3.connect(opsim_db_file) as conn:
            opsim_df = pd.read_sql(
                '''select obsHistID, descDitheredRA, descDitheredDec, filter,
                   expMJD from Summary''', conn)
        opsim_df['ra'] = np.degrees(opsim_df['descDitheredRA'])
        opsim_df['dec'] = np.degrees(opsim_df['descDitheredDec'])

        # Create the Variability Truth table.
        table_name = 'agn_variability_truth'
        cmd = f'''CREATE TABLE IF NOT EXISTS {table_name}
                  (id TEXT, obsHistID INTEGER, MJD FLOAT, bandpass TEXT,
                  delta_flux FLOAT, magNorm FLOAT)'''
        with sqlite3.connect(self.outfile) as conn:
            cursor = conn.cursor()
            cursor.execute(cmd)
            conn.commit()

            # Loop over rows in AGN db and add the flux for each
            # observation where the AGN is observed by LSST.
            sed_file = find_sed_file('agnSED/agn.spec.gz')
            agn_db_curs = self.conn.execute(self.query)
            row = agn_db_curs.fetchone()
            num_rows = 0
            while row is not None:
                # Extract the AGN info and model parameters.
                object_id = self.object_id(row[self.icol['galaxy_id']])
                ra = row[self.icol['ra']]
                dec = row[self.icol['dec']]
                magNorm = row[self.icol['magNorm']]
                redshift = row[self.icol['redshift']]
                params = json.loads(row[self.icol['varParamStr']])['p']
                seed = params['seed']

                # Compute baseline fluxes in each band.
                synth_phot = SyntheticPhotometry(sed_file, magNorm, redshift)
                gAv, gRv = synth_phot.add_MW_dust(ra, dec)
                flux0 = {band: synth_phot.calcFlux(band) for band in bands}

                # Select the visits from the opsim db in which the AGN
                # is observed by applying cuts on the sky coordinates.
                dec_cut = f'{dec - fp_radius} <= dec <= {dec + fp_radius}'
                df = pd.DataFrame(opsim_df.query(dec_cut))
                df['ang_sep'] = angularSeparation(df['ra'].to_numpy(),
                                                  df['dec'].to_numpy(), ra, dec)
                df = df.query(f'ang_sep <= {fp_radius}')
                # Compute delta fluxes for each band.
                for band in bands:
                    tau = params[f'agn_tau_{band}']
                    sf = params[f'agn_sf_{band}']
                    my_df = df.query(f'filter == "{band}"')
                    if len(my_df) == 0:
                        continue
                    obsHistIDs = my_df['obsHistID'].to_list()
                    mjds = my_df['expMJD'].to_numpy()
                    mag_norms = (agn_mag_norms(mjds, redshift, tau, sf, seed)
                                 + magNorm)
                    values = []
                    for obsHistID, mjd, mag_norm in zip(obsHistIDs, mjds,
                                                        mag_norms):
                        synth_phot = SyntheticPhotometry(sed_file, mag_norm,
                                                         redshift=redshift,
                                                         gAv=gAv, gRv=gRv)
                        delta_flux = synth_phot.calcFlux(band) - flux0[band]
                        values.append((object_id, obsHistID, mjd, band,
                                       delta_flux, mag_norm))
                    cursor.executemany(f'''INSERT INTO {table_name} VALUES
                                           (?, ?, ?, ?, ?, ?)''', values)
                    conn.commit()
                    num_rows += len(values)
                    if verbose:
                        print(num_rows)
                    if max_rows is not None and num_rows > max_rows:
                        return
                row = agn_db_curs.fetchone()
