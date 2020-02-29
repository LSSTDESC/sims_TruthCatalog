"""
Module to write truth catalogs for stars using the star parameters db as input.
"""
import os
import sqlite3
from multiprocessing import Process, Lock

import numpy as np
import pandas as pd
from .write_sqlite import write_sqlite
from .synthetic_photometry import SyntheticPhotometry, find_sed_file


__all__ = ['StarTruthWriter']

def _process_chunk(db_lock, star_db_file, query_string, outfile,
                   chunk_size, verbose, process_num, dry_run):

    if (verbose or dry_run):
        print("_process_chunk called for process #{} ".format(process_num))
        print("query_string: \n{}".format(query_string))
    if (dry_run):
        return
    
    with sqlite3.connect(star_db_file) as conn:
        cur = conn.execute(query_string)

        icol = {_[0]: icol for icol, _ in enumerate(cur.description)}
        irec = 0
        while True:
            ids, galaxy_ids, ra, dec, redshift = [], [], [], [], []
            is_variable, is_pointsource, good_ixes = [], [], []
            flux_by_band_MW = {_: [] for _ in 'ugrizy'}
            flux_by_band_noMW = {_: [] for _ in 'ugrizy'}
            chunk = cur.fetchmany(chunk_size)
            if not chunk:
                # No more rows to retrieve so exit the while loop.
                break
            for row in chunk:
                if verbose:
                    print('Process# {}: irec={}'.format(process_num, irec))
                irec += 1
                redshift.append(0)  # All stars are at redshift = 0.
                is_pointsource.append(1)  # All stars are point sources.
                ids.append(str(row[icol['simobjid']]))
                galaxy_ids.append(-1)
                ra.append(row[icol['ra']])
                dec.append(row[icol['decl']])
                if row[icol['varParamStr']] == 'None':
                    is_variable.append(0)
                else:
                    is_variable.append(1)
                sed_file = find_sed_file(row[icol['sedFilename']])

                # Create SyntheticPhotometry object initially without
                # Milky Way dust parameters.
                synth_phot = SyntheticPhotometry(sed_file,
                                                 row[icol['magNorm']],
                                                 redshift[-1])
                for band in 'ugrizy':
                    flux_by_band_noMW[band].append(synth_phot.calcFlux(band))

                # Set Milky Way dust parameters and compute ugrizy fluxes.
                synth_phot.add_MW_dust(ra[-1], dec[-1], Rv=3.1)
                for band in 'ugrizy':
                    flux_by_band_MW[band].append(synth_phot.calcFlux(band))

            # Get the output lock
            if not db_lock.acquire(timeout=120.0):
                print('Process {} failed to acquire output lock'.format(process_num))
                exit(1)
            write_sqlite(outfile,
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
            db_lock.release()


###
class StarTruthWriter:
    '''
    Write Summary and Variable truth tables for stars.
    '''
    
    def __init__(self, outfile, star_db_file, radec_bounds=None,
                 row_limit=None, max_parallel=20, dry_run=False):
        """
        Parameters
        ----------
        outfile: str
            Filename of output sqlite3 file to contain the truth_summary
            table.
        star_db_file: str
            Sqlite3 db file containing the information on the properties
            of each star.
        radec_bounds: (float, float, float, float) [None]
            Selection region in degrees as (ra_min, ra_max, dec_min, dec_max).
            If None, then no selection on ra, dec will be made.
        row_limit: int [None]
            Limit on number of rows to return from the query to the star
            database.  If None, then no limit will be applied.
        max_parallel: int [20]
            Limit on number of spawned processes running in parallel
        dry_run:    boolean [False]. If True just go through some motions
        """
        self.outfile = outfile
        self.dry_run = dry_run
        self.max_parallel = max_parallel
        self.star_db_file = star_db_file
        self.row_limit = row_limit
        if os.path.isfile(outfile):
            raise OSError(f'{outfile} already exists.')
        if not os.path.isfile(star_db_file):
            raise FileNotFoundError(f'{star_db_file} not found.')
        # Normal sqlite tables have special 64-bit integer column rowid which
        # which is a primary key.  Searches using conditions on it are
        # faster than using an ordinary user-declared index
        with sqlite3.connect(star_db_file) as conn:
           curs = conn.execute('select max(rowid) from stars')
           self.rowid_max = curs.fetchone()[0]
        self.per_process = int((self.rowid_max + max_parallel - 1)/max_parallel)
        
        if row_limit is not None:
            self.per_process = min(self.per_process, row_limit)

        ###self.conn = sqlite3.connect(star_db_file)
        self.query = '''select simobjid, ra, decl, varParamStr, sedFilename,
                     magNorm from stars where '''
        if radec_bounds is not None:
            self.query += (f' {radec_bounds[0]} <= ra and ' +
                           f'ra <= {radec_bounds[1]} and ' +
                           f'{radec_bounds[2]} <= decl and ' +
                           f'decl <= {radec_bounds[3]} and ')

        # self.curs = self.conn.execute(query)
        # self.icol = {_[0]: icol for icol, _ in enumerate(self.curs.description)}
        # Also will have to add a condition on range of rowid
        # and need to order by rowid

    def write(self, chunk_size=10000, verbose=False):
        '''
        Spawn max_parallel processes, each of which handles a fraction
        of the data; i.e., it makes a quety to extract its part of the 
        data from the star db file and writes to the output sqlite file.

        Parameters
        ----------
        chunk_size: int [10000]
            Number of records to read in at a time from the star db
            file and write to the output file.
        verbose: bool [False]
            Flag to write the number of records that have been processed.
        '''

        db_lock = Lock()
        p_list = []
        # For each process to be spawned, assemble query string
        for i_p in range(self.max_parallel):
            row_cut = 'rowid > {} and rowid <= {} order by rowid'.format(i_p * self.per_process, (i_p + 1) * self.per_process)
            i_query = self.query + row_cut

            if self.row_limit is not None:
                chunk_size = min(chunk_size, self.row_limit)
            p = Process(target=_process_chunk, name='proc_{}'.format(i_p),
                        args=(db_lock, self.star_db_file, i_query, self.outfile,
                              chunk_size, verbose, i_p,
                              self.dry_run))
            p.start()
            p_list.append(p)

        for p in p_list:
            p.join()

        print('All done with star truth')
