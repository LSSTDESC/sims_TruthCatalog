"""
Module to write truth catalogs for SNe using the SNe parameters db as input.
"""
import os
import sqlite3
import numpy as np
import pandas as pd
from .write_sqlite import write_sqlite


__all__ = ['SNeTruthWriter']


class SNeTruthWriter:
    '''
    Write Summary and Variable truth tables for SNe.
    '''
    def __init__(self, outfile, sne_db_file):
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
