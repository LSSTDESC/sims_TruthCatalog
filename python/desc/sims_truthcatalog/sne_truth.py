"""
Module to write truth catalogs for SNe using the SNe parameters db as input.
"""
import os
import sqlite3
import numpy as np
import pandas as pd
from lsst.sims.catUtils.supernovae import SNObject
from lsst.sims.photUtils import BandpassDict
from lsst.sims.utils import angularSeparation
from .synthetic_photometry import SyntheticPhotometry
from .write_sqlite import write_sqlite


__all__ = ['SNeTruthWriter']


class SNSynthPhotFactory:
    """
    Factory class to return the SyntheticPhotometry objects for a SN Ia
    as a function of time.
    """
    lsst_bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
    def __init__(self, row, bp_dict=None):
        """
        Parameters
        ----------
        row: pandas.Series
             Row of a pandas DataFrame that contains the sne_params db table.
        bp_dict: lsst.sims.photUtils.BandpassDict
             Dictionary containing the bandpasses. If None, the standard
             LSST total bandpasses will be used.
        """
        self.sn_obj = SNObject(row.snra_in, row.sndec_in)
        self.sn_obj.set(z=row.z_in, t0=row.t0_in, x0=row.x0_in,
                        x1=row.x1_in, c=row.c_in, hostebv=0, hostr_v=3.1,
                        mwebv=0, mwr_v=3.1)
        self.bp_dict = self.lsst_bp_dict if bp_dict is None else bp_dict

    def create(self, mjd):
        """
        Return a SyntheticPhotometry object for the specified time in MJD.

        Parameters
        ----------
        mjd: float
            Observation time in MJD.

        Returns
        -------
        desc.sims_truthcatalog.SyntheticPhotometry
        """
        sed = self.sn_obj.SNObjectSED(mjd, bandpass=self.bp_dict,
                                      applyExtinction=False)
        synth_phot = SyntheticPhotometry.create_from_sed(sed, redshift=0)
        synth_phot.add_MW_dust(*np.degrees(self.sn_obj.skycoord).flatten())
        return synth_phot

    def __getattr__(self, attr):
        return getattr(self.sn_obj, attr)


class SNeTruthWriter:
    '''
    Write Summary and Variable truth tables for SNe.
    '''
    def __init__(self, outfile, sne_db_file):
        """
        Parameters
        ----------
        outfile: str
            Name of the sqlite3 file to contain the truth tables.
        sne_db_file: str
            The sqlite3 file containing the SNe model parameters.
        """
        self.outfile = outfile
        if os.path.isfile(outfile):
            raise OSError(f'{outfile} already exists.')
        if not os.path.isfile(sne_db_file):
            raise FileNotFoundError(f'{sne_db_file} not found.')
        with sqlite3.connect(sne_db_file) as conn:
            self.sne_df = pd.read_sql('select * from sne_params', conn)

    def write(self):
        '''
        Extract the column data from the SNe db file and write the
        sqlite file.
        '''
        zeros = np.zeros(len(self.sne_df))
        ones = np.ones(len(self.sne_df))
        write_sqlite(self.outfile,
                     ids=self.sne_df['snid_in'],
                     galaxy_ids=self.sne_df['galaxy_id'],
                     ra=self.sne_df['snra_in'],
                     dec=self.sne_df['sndec_in'],
                     redshift=self.sne_df['z_in'],
                     is_variable=ones,
                     is_pointsource=ones,
                     flux_by_band_MW={_: zeros for _ in 'ugrizy'},
                     flux_by_band_noMW={_: zeros for _ in 'ugrizy'},
                     good_ixes=range(len(self.sne_df)))

    def write_auxiliary_truth(self):
        """
        Write the SN Auxiliary truth table from the SNe db file.
        This is almost a direct transcription of the columns.
        """
        table_name = 'sn_auxiliary_info'
        cmd = f'''CREATE TABLE IF NOT EXISTS {table_name}
        (id TEXT, host_galaxy_id BIGINT, ra DOUBLE, dec DOUBLE,
        t0 FLOAT, x0 FLOAT, x1 FLOAT, c FLOAT, redshift FLOAT)'''
        with sqlite3.connect(self.outfile) as conn:
            cursor = conn.cursor()
            cursor.execute(cmd)
            conn.commit()
            values = ((str(self.sne_df['snid_in'][i_obj]),
                       int(self.sne_df['galaxy_id'][i_obj]),
                       self.sne_df['snra_in'][i_obj],
                       self.sne_df['sndec_in'][i_obj],
                       self.sne_df['t0_in'][i_obj],
                       self.sne_df['x0_in'][i_obj],
                       self.sne_df['x1_in'][i_obj],
                       self.sne_df['c_in'][i_obj],
                       self.sne_df['z_in'][i_obj])
                      for i_obj in range(len(self.sne_df)))
            cursor.executemany(f'''INSERT INTO {table_name}
                                   VALUES (?,?,?,?,?,?,?,?,?)''', values)
            conn.commit()

    def write_variability_truth(self, opsim_db_file, fp_radius=2.1,
                                max_rows=None):
        """
        Write the Variability Truth Table. This will contain light curve
        points for visits in which the SN was being observed by LSST.

        Parameters
        ----------
        opsim_db_file: str
            The sqlite3 file containing the OpSim Summary table which
            has the pointing information for each visit.
        fp_radius: float [2.1]
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
        # Read in pointing and filter info from the OpSim db file.
        with sqlite3.connect(opsim_db_file) as conn:
            opsim_df = pd.read_sql(
                '''select obsHistID, descDitheredRA, descDitheredDec, filter,
                   expMJD from Summary''', conn)
        opsim_df['ra'] = np.degrees(opsim_df['descDitheredRA'])
        opsim_df['dec'] = np.degrees(opsim_df['descDitheredDec'])

        # Create the Variability Truth table.
        table_name = 'sn_variability_truth'
        cmd = f'''CREATE TABLE IF NOT EXISTS {table_name}
                  (id TEXT, obsHistID INT, MJD FLOAT, bandpass TEXT,
                  delta_flux FLOAT)'''
        with sqlite3.connect(self.outfile) as conn:
            cursor = conn.cursor()
            cursor.execute(cmd)
            conn.commit()

            # Loop over rows in the SN database and add the flux for
            # each observation where the SN is observed by LSST.
            num_rows = 0
            for iloc in range(len(self.sne_df)):
                row = self.sne_df.iloc[iloc]
                sp_factory = SNSynthPhotFactory(row)

                # Make cuts on time based on sncosmo valid model range
                # and on allowed range of Declination values.
                tmin, tmax = sp_factory.mintime(), sp_factory.maxtime()
                dmin, dmax = row.sndec_in - fp_radius, row.sndec_in + fp_radius
                df = pd.DataFrame(
                    opsim_df.query(f'{tmin} <= expMJD <= {tmax} and '
                                   f'{dmin} <= dec <= {dmax}'))

                # Compute angular separation from SN position to use
                # for applying acceptance cone cut.
                df['ang_sep'] = angularSeparation(df['ra'].to_numpy(),
                                                  df['dec'].to_numpy(),
                                                  row.snra_in, row.sndec_in)
                df = df.query(f'ang_sep <= {fp_radius}')

                # Insert the rows into the variability truth table.
                values = []
                for visit, band, mjd in zip(df['obsHistID'], df['filter'],
                                            df['expMJD']):
                    synth_phot = sp_factory.create(mjd)
                    values.append((row.snid_in, visit, mjd, band,
                                   synth_phot.calcFlux(band)))
                cursor.executemany(f'''INSERT INTO {table_name}
                                       VALUES (?,?,?,?,?)''', values)
                conn.commit()
                num_rows += len(values)
                if max_rows is not None and num_rows > max_rows:
                    break
