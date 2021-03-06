"""
Module to write truth catalogs for SNe using the SNe parameters db as input.
"""
import os
import sys
import re
import sqlite3
import numpy as np
import pandas as pd
import warnings
from multiprocessing import Process, Lock
from yaml import load as yload
from yaml import FullLoader

from lsst.sims.catUtils.supernovae import SNObject
from lsst.sims.photUtils import BandpassDict
from lsst.sims.utils import angularSeparation
from .synthetic_photometry import SyntheticPhotometry
from .write_sqlite import write_sqlite
from .sqlite_utils import merge_sqlite3_dbs
from .script_utils import print_date as print_date


__all__ = ['SNeTruthWriter', 'SNSynthPhotFactory', 'get_chunk_intervals',
           'sne_create_var_index', 'sne_merge_var']


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

    
def _process_chunk(outfile, logfile, sne_db_file, sne_db_query, opsim_df,
                   fp_radius, process_num, max_rows, dry_run, verbose):
    '''
    Compute and write out variability for a manageable chunk of SNe
    '''
    if logfile is not None:
        lg = open(logfile, 'w')
    else:
        lg = sys.stdout
    print(f"_process_chunk called for process #{process_num} ", file=lg)
    print(f"query_string: \n{sne_db_query}", file=lg)

    # table_name should perhaps be yet another argument
    table_name = 'sn_variability_truth'
    
    if dry_run:           
        print(f'Dry run: "finish" chunk {process_num}', file=lg)
        return

    # Create the Variability Truth table.
    table_name = 'sn_variability_truth'
    cmd = f'''CREATE TABLE IF NOT EXISTS {table_name}
    (id TEXT, obsHistID INT, MJD FLOAT, bandpass TEXT,
    delta_flux FLOAT)'''
    with sqlite3.connect(outfile) as conn:
        cursor = conn.cursor()
        cursor.execute(cmd)
        conn.commit()
    
    # Get our batch from sne_db_file
    with sqlite3.connect(sne_db_file) as read_conn:
        sne_df = pd.read_sql(sne_db_query, read_conn)
    if verbose:
        print(f"Process {process_num} has  {len(sen_df)} sne", file=lg)

    # Loop over rows in the SN database and add the flux for
    # each observation where the SN is observed by LSST.
    num_rows = 0

    for iloc in range(len(sne_df)):
        row = sne_df.iloc[iloc]
        if iloc % 10000 == 0:
            print(f'Process {process_num} is on its {iloc} object having id {row.snid_in}',
                  file=lg)
            print_date(file=lg)
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
                print(f"Process {process_num}: no rows found for id {row.snid_in}", file=lg)                
            continue

        with sqlite3.connect(outfile) as out_conn:
            cursor = out_conn.cursor()
            cursor.executemany(f'''INSERT INTO {table_name} VALUES (?,?,?,?,?)''', values)
            out_conn.commit()
        
        num_rows += len(values)
        if verbose:
            print(f'Process {process_num} inserting {len(values)} rows for id {row.snid_in} ',
                  file=lg)
        if max_rows is not None and num_rows > max_rows:
            break
    print(f'Done with chunk {process_num}', file=lg)

def _get_chunk_intervals(fpath, max_chunks):
    '''
    Parameters:
          fpath (string) Path to properly-formatted yaml file containing
             intervals and (optional) max # SNe per chunk
          max_chunks (int) max # of returned intervals allowed
    Returns:  
          intervals - list of (start, end) values for chunks
    '''
    y = yload(open(fpath), Loader=FullLoader)

    if not isinstance(y, dict):
        raise ValueError('variability yaml file is not a dict!')
    if not 'intervals' in y:
        raise ValueError('variability yaml file missing intervals keyword')
    intervals_yaml = y['intervals']
    if not isinstance(intervals_yaml, list):
        raise ValueError('variability yaml file has improper intervals list')
    if len(intervals_yaml) > max_chunks:
        raise ValueError('variability yaml file has too many intervals')
    intervals = []
    lens = []
    ix = 0
    for iy in intervals_yaml:
        start = int(iy['start'])
        end = int(iy['end'])
        if (start < 0) or (start > end):
            raise ValueError('variability yaml file has bad interval')
        intervals.append((start, end))
        lens.append(end + 1 - start)
        ix += 1
                                                                         
    if not 'max_chunk_size' in y:
        return intervals
    mcs = y['max_chunk_size']
    if mcs >= max(lens): return intervals

    # Split intervals as needed up to max_chunks
    ret_intervals = []
    remaining = max_chunks - len(intervals)
    if remaining == 0:
        return intervals
    for i in range(len(intervals)):
        if remaining == 0 or lens[i] <= mcs:
            ret_intervals.append(intervals[i])
            continue
        # split
        to_split = int(min(remaining + 1, np.ceil(lens[i]/mcs)))
        add_for_end = int(np.ceil(lens[i]/to_split)) - 1
        start = int(intervals[i][0])
        for k in range(to_split):
            end = min(start+add_for_end, intervals[i][1])
            ret_intervals.append((start, end))
            start = end + 1
        remaining -= (to_split - 1)

    return ret_intervals

# Public version for debugging
def get_chunk_intervals(fpath, max_chunks):
    return _get_chunk_intervals(fpath, max_chunks)


class SNeTruthWriter:
    '''
    Write Summary and Variable truth tables for SNe.
    '''
    def __init__(self, outfile, sne_db_file, max_parallel=1, dry_run=False,
                 no_summary=False):
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
        no_summary: boolean
            Do not write summary table
        """
        self.outfile = outfile
        if os.path.isfile(outfile) and not no_summary:
            raise OSError(f'{outfile} already exists.')
        if not os.path.isfile(sne_db_file):
            raise FileNotFoundError(f'{sne_db_file} not found.')
        self.dry_run = dry_run
        self.sne_db_file = sne_db_file

        q = 'select galaxy_id, c_in, t0_in, x0_in, x1_in, z_in, snid_in, snra_in, sndec_in from sne_params order by rowid '
        print('The query: ', q)
        with sqlite3.connect(sne_db_file) as conn:
            self.sne_df = pd.read_sql(q, conn)

            curs = conn.execute('select max(rowid) from sne_params')
            self.rowid_max = curs.fetchone()[0]

        # Used for normal case (no interval file) only
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

    def write_variability_truth(self, opsim_db_file, chunk_log=None,
                                fp_radius=2.05,
                                max_rows=None, max_parallel=1, verbose=False,
                                interval_file=None):
        """
        Write the Variability Truth Table. This will contain light curve
        points for visits in which the SN was being observed by LSST.

        Parameters
        ----------
        opsim_db_file: str
            The sqlite3 file containing the OpSim Summary table which
            has the pointing information for each visit.
        chunk_log: str
            File(s) where per-chunk output is written.  If '+' is not in the
            string a '+' will be appended.  Then, for each process handling
            a chunk, replace '+' with the process number and log to that path.
            If None write to standard output
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
        interval_file: string [None]
            if not None contains list of intervals (rows within SNe input
            db) to be processed.  Else process all.  May also specify
            max chunk size, in which case intervals may be split
        """
        # Read in pointing and filter info from the OpSim db file.
        with sqlite3.connect(opsim_db_file) as conn:
            opsim_df = pd.read_sql(
                '''select obsHistID, descDitheredRA, descDitheredDec, filter,
                   expMJD from Summary''', conn)
        opsim_df['ra'] = np.degrees(opsim_df['descDitheredRA'])
        opsim_df['dec'] = np.degrees(opsim_df['descDitheredDec'])
        _opsim_df = opsim_df
        p_list = []

        if chunk_log is not None:
            if not '+' in chunk_log: chunk_log += '+'

        # Same query as before except don't need galaxy_id and include conditions on rowid
        q = 'select c_in, t0_in, x0_in, x1_in, z_in, snid_in, snra_in, sndec_in from sne_params '

        n_proc = max_parallel
        if interval_file is not None:
            intervals = _get_chunk_intervals(interval_file,max_parallel)
            n_proc = len(intervals)
            
        for i_p in range(n_proc):
            if interval_file is not None:
                # compose row_cut using interval information
                row_cut = f'where rowid >= {intervals[i_p][0]} and rowid <= {intervals[i_p][1]} order by rowid'
            else:
                row_cut = f'where rowid > {i_p*self.per_process} and rowid <= {(i_p+1)*self.per_process} order by rowid'
            i_query = q + row_cut

            if chunk_log is None:
                i_out = self.outfile
                i_log = None
            else:
                i_out = f'{self.outfile}_{i_p}'
                i_log = chunk_log.replace('+', str(i_p))
                
            p = Process(target=_process_chunk, name=f'proc{i_p}',
                        args=(i_out, i_log, self.sne_db_file, i_query, _opsim_df,
                              fp_radius, i_p, max_rows, self.dry_run,
                              verbose))
            p.start()
            p_list.append(p)

        for p in p_list:
            p.join()

def sne_create_var_index(db_file):
    with sqlite3.connect(db_file) as conn:
        conn.cursor().execute('create index snid_ix on sn_variability_truth(id)')
        conn.commit()

def sne_merge_var(outfile, in_dir, in_pattern):
    '''
    Merge collection into outfile
    Parameters
    ----------
    outfile:  May or may not already exist.  May or may not already have
              variability table.  If so, schema had better match.
              Also assuming all infiles have the same schema
    in_dir:   dbs to be merged are in this directory
    in_pattern:  regular expresssion. Matches files we need from the directory
                 Note:  it has to be a full match

    Return
    ------
    True if merge was done; else False
    '''

    # Get collection of input files
    dir_contents = os.listdir(in_dir)
    infiles = []
    for f in dir_contents:
        if re.fullmatch(in_pattern, f):
            infiles.append(os.path.join(in_dir,f))

    if len(infiles) == 0:
        warnings.warn('No matching files to merge')
        return False

    if not os.path.exists(outfile):
        merge_sqlite3_dbs(infiles, outfile)
        return  True

    out_dir = os.path.dirname(outfile)
    tmp_merge = os.path.join(out_dir,'merge_sqlite3_dbs.db_tmp' )
    while os.path.exists(tmp_merge):
        tmp_merge += '_tmp'
    merge_sqlite3_dbs(infiles, tmp_merge)

    # get table definition
    tbl = 'sn_variability_truth'
    with sqlite3.connect(tmp_merge) as conn:
        c = conn.cursor()
        c.execute(f"select sql from sqlite_master where name='{tbl}'")
        sch = c.fetchone()[0][len('CREATE TABLE'):]

    create_string = 'CREATE TABLE IF NOT EXISTS ' + sch
    with sqlite3.connect(outfile) as conn:
        conn.execute(create_string)
        conn.execute(f"attach '{tmp_merge}' as in_db")
        conn.execute(f"insert into {tbl} select * from in_db.{tbl}")
        conn.commit()
        conn.execute("detach database in_db")
    os.remove(tmp_merge)
    
    return True
        
        

           
