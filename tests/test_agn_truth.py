"""
Unit test code for agn_truth.py module.
"""
import os
import unittest
import sqlite3
import pandas as pd
from desc.sims_truthcatalog import AGNTruthWriter


class AGNTruthWriterTestCase(unittest.TestCase):
    '''
    Test case class for AGNTruthWriter.
    '''
    def setUp(self):
        self.outfile = 'agn_truth_cat.db'
        self.agn_db_file = os.path.join(os.environ['SIMS_TRUTHCATALOG_DIR'],
                                        'data', 'agn_cosmoDC2_v1.1.4_small.db')

    def tearDown(self):
        if os.path.isfile(self.outfile):
            os.remove(self.outfile)

    def test_write(self):
        '''
        Test the AGNTruthWriter class.
        '''
        # Check exceptions in __init__.
        with self.assertRaises(FileNotFoundError):
            AGNTruthWriter(self.outfile, self.agn_db_file[:-1])

        agn_truth_writer = AGNTruthWriter(self.outfile, self.agn_db_file)
        agn_truth_writer.write()

        with self.assertRaises(OSError):
            AGNTruthWriter(self.outfile, self.agn_db_file)

        with sqlite3.connect(self.outfile) as conn:
            df = pd.read_sql('select * from truth_summary', conn)

        row = df.iloc[0]
        self.assertEqual(row['host_galaxy'], 1250442285)
        self.assertAlmostEqual(row['redshift'], 0.517463326454163)

if __name__ == '__main__':
    unittest.main()
