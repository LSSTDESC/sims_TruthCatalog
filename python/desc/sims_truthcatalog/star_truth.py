"""
Module to write truth catalogs for stars using the star parameters db as input.
"""
import os
import sqlite3
import numpy as np
import pandas as pd
from .write_sqlite import write_sqlite
from .synthetic_photometry import SyntheticPhotometry, find_sed_file


__all__ = ['StarTruthWriter']


class StarTruthWriter:
    '''
    Write Summary and Variable truth tables for stars.
    '''
    def __init__(self, outfile, star_db_file, radec_bounds=None,
                 row_limit=None):
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
        """
        self.outfile = outfile
        if os.path.isfile(outfile):
            raise OSError(f'{outfile} already exists.')
        if not os.path.isfile(star_db_file):
            raise FileNotFoundError(f'{star_db_file} not found.')
        self.conn = sqlite3.connect(star_db_file)
        query = '''select simobjid, ra, decl, varParamStr, sedFilename,
                magNorm from stars'''
        if radec_bounds is not None:
            query += (f' where {radec_bounds[0]} <= ra and ' +
                      f'ra <= {radec_bounds[1]} and ' +
                      f'{radec_bounds[2]} <= decl and ' +
                      f'decl <= {radec_bounds[3]}')
        if row_limit is not None:
            query += f' limit {row_limit}'
        self.curs = self.conn.execute(query)
        self.icol = {_[0]: icol for icol, _ in enumerate(self.curs.description)}

    def write(self, chunk_size=10000, verbose=False):
        '''
        Extract the column data from the star db file and write the
        sqlite file.

        Parameters
        ----------
        chunk_size: int [10000]
            Number of records to read in at a time from the star db
            file and write to the output file.
        verbose: bool [False]
            Flag to write the number of records that have been processed.
        '''
        irec = 0
        while True:
            ids, galaxy_ids, ra, dec, redshift = [], [], [], [], []
            is_variable, is_pointsource, good_ixes = [], [], []
            flux_by_band_MW = {_: [] for _ in 'ugrizy'}
            flux_by_band_noMW = {_: [] for _ in 'ugrizy'}
            chunk = self.curs.fetchmany(chunk_size)
            if not chunk:
                # No more rows to retrieve so exit the while loop.
                break
            for row in chunk:
                if verbose:
                    print(irec)
                irec += 1
                redshift.append(0)  # All stars are at redshift = 0.
                is_pointsource.append(1)  # All stars are point sources.
                ids.append(str(row[self.icol['simobjid']]))
                galaxy_ids.append(-1)
                ra.append(row[self.icol['ra']])
                dec.append(row[self.icol['decl']])
                if row[self.icol['varParamStr']] == 'None':
                    is_variable.append(0)
                else:
                    is_variable.append(1)
                sed_file = find_sed_file(row[self.icol['sedFilename']])

                # Create SyntheticPhotometry object initially without
                # Milky Way dust parameters.
                synth_phot = SyntheticPhotometry(sed_file,
                                                 row[self.icol['magNorm']],
                                                 redshift[-1])
                for band in 'ugrizy':
                    flux_by_band_noMW[band].append(synth_phot.calcFlux(band))

                # Set Milky Way dust parameters and compute ugrizy fluxes.
                synth_phot.add_MW_dust(ra[-1], dec[-1], Rv=3.1)
                for band in 'ugrizy':
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
