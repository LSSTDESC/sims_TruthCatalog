"""
Unit tests for star truth catalog code.
"""
import os
import unittest
import sqlite3
import numpy as np
import pandas as pd
from desc.sims_truthcatalog import StarTruthWriter


class StarTruthWriterTestCase(unittest.TestCase):
    """
    Test case class for StarTruthWriter.
    """
    def setUp(self):
        self.outfile = 'star_truth_cat.db'
        self.star_db_file = os.path.join(os.environ['SIMS_TRUTHCATALOG_DIR'],
                                         'data', 'dc2_stars_small.db')

    def tearDown(self):
        if os.path.isfile(self.outfile):
            os.remove(self.outfile)

    def test_write(self):
        """
        Test the StarTruthWriter.
        """
        # Check exceptions in __init__.
        with self.assertRaises(FileNotFoundError):
            StarTruthWriter(self.outfile, self.star_db_file[:-1])

        star_truth_writer = StarTruthWriter(self.outfile, self.star_db_file)
        star_truth_writer.write()

        with self.assertRaises(OSError):
            StarTruthWriter(self.outfile, self.star_db_file)

        # Read in the output catalog and compute values for
        # down-selection options.
        query = 'select ra, dec from truth_summary'
        with sqlite3.connect(self.outfile) as conn:
            df = pd.read_sql(query, conn)
        nrows = len(df)
        ra_min, ra_max = min(df['ra']), max(df['ra'])
        dec_min, dec_max = min(df['dec']), max(df['dec'])

        os.remove(self.outfile)

        # Check application of row_limit.
        row_limit = nrows//2
        star_truth_writer = StarTruthWriter(self.outfile, self.star_db_file,
                                            row_limit=row_limit)
        star_truth_writer.write()
        with sqlite3.connect(self.outfile) as conn:
            df = pd.read_sql(query, conn)
        self.assertEqual(len(df), row_limit)

        os.remove(self.outfile)

        # Check application of limits on ra, dec values.
        dra = (ra_max - ra_min)/4
        ddec = (dec_max - dec_min)/4
        radec_bounds = (ra_min + dra, ra_max - dra,
                        dec_min + ddec, dec_max - ddec)
        star_truth_writer = StarTruthWriter(self.outfile, self.star_db_file,
                                            radec_bounds=radec_bounds)
        star_truth_writer.write()
        with sqlite3.connect(self.outfile) as conn:
            df = pd.read_sql(query, conn)
        self.assertLessEqual(radec_bounds[0], min(df['ra']))
        self.assertLessEqual(max(df['ra']), radec_bounds[1])
        self.assertLessEqual(radec_bounds[2], min(df['dec']))
        self.assertLessEqual(max(df['dec']), radec_bounds[3])


if __name__ == '__main__':
    unittest.main()
