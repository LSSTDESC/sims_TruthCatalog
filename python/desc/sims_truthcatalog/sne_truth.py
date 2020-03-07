"""
Module to write truth catalogs for SNe using the SNe parameters db as input.
"""
import os
import sqlite3
import numpy as np
import pandas as pd
from multiprocessing import Process, Lock
from lsst.sims.catUtils.supernovae import SNObject
from lsst.sims.photUtils import BandpassDict
from lsst.sims.utils import angularSeparation
from .synthetic_photometry import SyntheticPhotometry
from .write_sqlite import write_sqlite


__all__ = ['SNeTruthWriter', 'SNSynthPhotFactory']


# For use in parallelizing variability table
_opsim_df = None

class SNSynthPhotFactory:
    """
    Factory class to return the SyntheticPhotometry objects for a SN Ia
    as a function of time.  This uses the 'salt2-extended' model in
    sncosmo as implemented in sims_catUtils.
    """
    lsst_bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
    def __init__(self, z=0, t0=0, x0=1, x1=0, c=0, snra=0, sndec=0):
        """
        Parameters
        ----------
        z: float [0]
             Redshift of the SNIa object to pass to sncosmo.
        t0: float [0]
             Time in mjd of phase=0, corresponding to the B-band
             maximum.
        x0: float [1]
             Normalization factor for the lightcurves.
        x1: float [0]
             Empirical parameter controlling the stretch in time of the
             light curves
        c: float [0]
             Empirical parameter controlling the colors.
        snra: float [0]
             RA in degrees of the SNIa.
        sndec: float [0]
             Dec in degrees of the SNIa.
        """
        self.sn_obj = SNObject(snra, sndec)
        self.sn_obj.set(z=z, t0=t0, x0=x0, x1=x1, c=c,
                        hostebv=0, hostr_v=3.1, mwebv=0, mwr_v=3.1)
        self.bp_dict = self.lsst_bp_dict

    def set_bp_dict(self, bp_dict):
        """
        Set the bandpass dictionary.

        Parameters
        ----------
        bp_dict: lsst.sims.photUtils.BandpassDict
             Dictionary containing the bandpasses. If None, the standard
             LSST total bandpasses will be used.
        """
        self.bp_dict = bp_dict

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
        # Create the Sed object.  Milky Way extinction will be
        # below, so set applyExtinction=False.
        sed = self.sn_obj.SNObjectSED(mjd, bandpass=self.bp_dict,
                                      applyExtinction=False)
        # The redshift was applied to the model SED computed by
        # sncosmo in the __init__ function, so set redshift=0 here.
        synth_phot = SyntheticPhotometry.create_from_sed(sed, redshift=0)
        synth_phot.add_MW_dust(*np.degrees(self.sn_obj.skycoord).flatten())
        return synth_phot

    def __getattr__(self, attr):
        # Pass any attribute access requests to the SNObject instance.
        return getattr(self.sn_obj, attr)

    
def _process_chunk(db_lock, outfile, sne_db_file, sne_db_query, opsim_df,
                   fp_radius, process_num, max_rows, dry_run, verbose):
    '''
    Compute and write out variability for a manageable chunk of SNe
    '''

    if verbose or dry_run:
        print("_process_chunk called for process #{} ".format(process_num))
        print("query_string: \n{}".format(sne_db_query))

    # table_name should perhaps be yet another argument
    table_name = 'sn_variability_truth'
    
    if dry_run:
        return

    # Get our batch from sne_db_file
    with sqlite3.connect(sne_db_file) as conn:
        sne_df = pd.read_sql(sne_db_query, conn)
    if (verbose):
        print("Process {} has  {} sne".format(process_num, len(sne_df)))

    # Loop over rows in the SN database and add the flux for
    # each observation where the SN is observed by LSST.
    num_rows = 0

    for iloc in range(len(sne_df)):
        row = sne_df.iloc[iloc]
        params = {_: row[f'{_}_in'] for _ in
                  'z t0 x0 x1 c snra sndec'.split()}
        sp_factory = SNSynthPhotFactory(**params)

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
        if len(values) == 0:
            if (verbose):
                print("Process {}: no rows found for id {}".format(process_num, row.snid_in))                
            continue
        
        if not db_lock.acquire(timeout=60.0):
            print('Process {} failed to acquire outpout lock'.format(i_p))
            exit(1)

        with sqlite3.connect(outfile) as out_conn:
            cursor = out_conn.cursor()
        

            cursor.executemany(f'''INSERT INTO {table_name} VALUES (?,?,?,?,?)''', values)
            out_conn.commit()
            
        db_lock.release()
        
        num_rows += len(values)
        if verbose:
            print('Process {} inserting {} rows for id {} '.format(process_num,
                                                                   len(values),
                                                                   row.snid_in))
        if max_rows is not None and num_rows > max_rows:
            break

class SNeTruthWriter:
    '''
    Write Summary and Variable truth tables for SNe.
    '''
    def __init__(self, outfile, sne_db_file, max_parallel=1, dry_run=False):
        """
        Parameters
        ----------
        outfile: str
            Name of the sqlite3 file to contain the truth tables.
        sne_db_file: str
            The sqlite3 file containing the SNe model parameters.
        sne_limit: int 
            Maximum number of SNe to process. Default: no limit
        dry_run: boolean
            No actual db write.  Default: False
        """
        self.outfile = outfile
        if os.path.isfile(outfile):
            raise OSError(f'{outfile} already exists.')
        if not os.path.isfile(sne_db_file):
            raise FileNotFoundError(f'{sne_db_file} not found.')
        self.dry_run = dry_run
        self.sne_db_file = sne_db_file

        q = 'select galaxy_id, c_in, t0_in, x0_in, x1_in, z_in, snid_in, snra_in, sndec_in from sne_params order by rowid '
        print('The query: ', q)
        with sqlite3.connect(sne_db_file) as conn:
            #self.sne_df = pd.read_sql('select * from sne_params', conn)
            self.sne_df = pd.read_sql(q, conn)

            curs = conn.execute('select max(rowid) from sne_params')
            self.rowid_max = curs.fetchone()[0]

        self.per_process = int((self.rowid_max + max_parallel - 1)/max_parallel)
        
    def write(self):
        '''
        Extract the column data from the SNe db file and write the
        sqlite file.
        '''
        if (self.dry_run):
            return
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
        if (self.dry_run):
            return
        table_name = 'sn_auxiliary_info'
        cmd = f'''CREATE TABLE IF NOT EXISTS {table_name}
        (id TEXT, host_galaxy BIGINT, ra DOUBLE, dec DOUBLE,
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

    def write_variability_truth(self, opsim_db_file, fp_radius=2.05,
                                max_rows=None, max_parallel=1, verbose=False):
        """
        Write the Variability Truth Table. This will contain light curve
        points for visits in which the SN was being observed by LSST.

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
            Threshold number of rows per process to write to the table.  This is useful
            for testing.  If None, then write all entries for all SNe in
            the sne_db_file.
        max_parallel: int [1]
            split work among this number of processes
        verbose: boolean [False]
            if True write debug output
        """
        # Read in pointing and filter info from the OpSim db file.
        with sqlite3.connect(opsim_db_file) as conn:
            opsim_df = pd.read_sql(
                '''select obsHistID, descDitheredRA, descDitheredDec, filter,
                   expMJD from Summary''', conn)
        opsim_df['ra'] = np.degrees(opsim_df['descDitheredRA'])
        opsim_df['dec'] = np.degrees(opsim_df['descDitheredDec'])

        _opsim_df = opsim_df

        # Create the Variability Truth table.
        table_name = 'sn_variability_truth'
        cmd = f'''CREATE TABLE IF NOT EXISTS {table_name}
                  (id TEXT, obsHistID INT, MJD FLOAT, bandpass TEXT,
                  delta_flux FLOAT)'''
        if not self.dry_run:
            with sqlite3.connect(self.outfile) as conn:
                cursor = conn.cursor()
                cursor.execute(cmd)
                conn.commit()

        db_lock = Lock()
        p_list = []

        # Same query as before except don't need galaxy_id and include conditions on rowid
        q = 'select c_in, t0_in, x0_in, x1_in, z_in, snid_in, snra_in, sndec_in from sne_params '
        
        for i_p in range(max_parallel):
            row_cut = 'where rowid > {} and rowid <= {} order by rowid'.format(i_p * self.per_process, (i_p + 1) * self.per_process)
            i_query = q + row_cut

            p = Process(target=_process_chunk, name='proc_{}'.format(i_p),
                        args=(db_lock, self.outfile, self.sne_db_file, i_query,
                              _opsim_df, fp_radius, i_p, max_rows, self.dry_run,
                              verbose))
            p.start()
            p_list.append(p)

        for p in p_list:
            p.join()

        # Finally, index to make look-up of light curves faster
        if not self.dry_run:
            with sqlite3.connect(self.outfile) as conn:
                conn.cursor().execute('create index snid_ix on sn_variability_truth(id)')
                conn.commit()
