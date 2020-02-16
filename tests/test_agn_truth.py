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
        self.data_dir \
            = os.path.join(os.environ['SIMS_TRUTHCATALOG_DIR'], 'data')
        self.agn_db_file \
            = os.path.join(self.data_dir, 'agn_cosmoDC2_v1.1.4_small.db')

    def tearDown(self):
        if os.path.isfile(self.outfile):
            os.remove(self.outfile)

    def test_truth_summary(self):
        '''
        Test the AGNTruthWriter.write method.
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

    def test_auxiliary_truth(self):
        '''
        Test the AGNTruthWriter.write_auxiliary_truth method.
        '''
        agn_truth_writer = AGNTruthWriter(self.outfile, self.agn_db_file)
        agn_truth_writer.write_auxiliary_truth()
        with sqlite3.connect(self.outfile) as conn:
            df = pd.read_sql('select * from agn_auxiliary_info', conn)

        row = df.iloc[0]
        galaxy_id = 1250442285
        self.assertEqual(row['id'], str(1024*galaxy_id + 117))
        self.assertEqual(row['host_galaxy'], galaxy_id)
        self.assertAlmostEqual(row['M_i'], -25.2830250010622)
        self.assertEqual(row['seed'], 2463655)

        row = df.iloc[len(df) - 1]
        galaxy_id = 1253710043
        self.assertEqual(row['id'], str(1024*galaxy_id + 117))
        self.assertEqual(row['host_galaxy'], galaxy_id)
        self.assertAlmostEqual(row['M_i'], -22.5594939688854)
        self.assertEqual(row['seed'], 7649944)

    def test_variability_truth(self):
        """
        Test some light curve values for AGNs from the test AGN catalog.
        """
        opsim_db_file = os.path.join(self.data_dir,
                                     'minion_1016_desc_dithered_v4_small.db')
        agn_truth_writer = AGNTruthWriter(self.outfile, self.agn_db_file)
        agn_truth_writer.write_variability_truth(opsim_db_file, max_rows=100)


if __name__ == '__main__':
    unittest.main()
