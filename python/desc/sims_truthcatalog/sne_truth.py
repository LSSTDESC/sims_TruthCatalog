"""
Module to write truth catalogs for SNe using the SNe parameters db as input.
"""
import os
import sqlite3
import numpy as np
import pandas as pd
from .galaxy_truth import GalaxyTruthWriter, _write_sqlite as write_sqlite


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
            raise OSError(f'{sne_db_file} not found.')
        with sqlite3.connect(sne_db_file) as conn:
            self.sne_df = pd.read_sql('select * from sne_params', conn)

    def write(self):
        '''
        Extract the column data from the SNe db file and write the
        sqlite file.
        '''
        flux_by_band = {_: np.zeros(len(self.sne_df)) for _ in 'ugrizy'}
        write_sqlite(self.outfile,
                     self.sne_df['galaxy_id'],
                     self.sne_df['snra_in'],
                     self.sne_df['sndec_in'],
                     self.sne_df['z_in'],
                     flux_by_band,
                     flux_by_band,
                     range(len(self.sne_df)),
                     is_variable=1, is_pointsource=1)
