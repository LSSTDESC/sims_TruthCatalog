"""
Unit tests for SNIa truth catalog code.
"""
import os
import unittest
import sqlite3
import numpy as np
import pandas as pd
from desc.sims_truthcatalog import SNeTruthWriter, SNSynthPhotFactory


class SNSynthPhotFactoryTestCase(unittest.TestCase):
    """
    Test case class for SNIa synthetic photometry factory class.
    """
    def test_SNSythPhotFactory(self):
        """
        Test some flux calculations using the underlying SNObject
        and SyntheticPhotometry classes.
        """
        sp_factory = SNSynthPhotFactory(z=0.6322702169418335,
                                        t0=61719.9950436545,
                                        x0=4.2832710977804034e-06,
                                        x1=-1.207738485943195,
                                        c=-0.0069750402968899936,
                                        snra=55.26407314527358,
                                        sndec=-40.81575605788344)
        mjds = (61689.150791, 61697.354470, 61712.258685)
        bands = ('z', 'i', 'r')
        fluxes = (2.6401569864737633, 71.18561504923377, 1048.0327802379868)
        for mjd, band, flux in zip(mjds, bands, fluxes):
            sp = sp_factory.create(mjd)
            self.assertAlmostEqual(sp.calcFlux(band), flux)


class SNeTruthWriterTestCase(unittest.TestCase):
    """
    Test case class for SNIa truth catalog generation class.
    """
    def setUp(self):
        self.outfile = 'test_sne_truth_cat.db'
        self.data_dir = os.path.join(os.environ['SIMS_TRUTHCATALOG_DIR'],
                                     'data')
        sn_db_file = os.path.join(self.data_dir,
                                  'sne_cosmoDC2_v1.1.4_MS_DDF_small.db')
        self.sne_truth_writer = SNeTruthWriter(self.outfile, sn_db_file)

    def tearDown(self):
        if os.path.isfile(self.outfile):
            os.remove(self.outfile)

    def test_truth_summary(self):
        """Test that the truth_summary columns are filled out as expected."""
        self.sne_truth_writer.write()
        with sqlite3.connect(self.outfile) as conn:
            df = pd.read_sql('select * from truth_summary', conn)
        zeros = np.zeros(len(df))
        ones = np.ones(len(df))
        np.testing.assert_equal(df['is_variable'], ones)
        np.testing.assert_equal(df['is_pointsource'], ones)
        for band in 'ugrizy':
            flux_col = f'flux_{band}'
            np.testing.assert_equal(df[flux_col], zeros)
            flux_col += '_noMW'
            np.testing.assert_equal(df[flux_col], zeros)

    def test_auxiliary_truth(self):
        """
        Test that the columns from the sne_params table are transcribed
        correctly.
        """
        self.sne_truth_writer.write_auxiliary_truth()
        with sqlite3.connect(self.outfile) as conn:
            df = pd.read_sql('select * from sn_auxiliary_info', conn)
        np.testing.assert_equal(self.sne_truth_writer.sne_df['snid_in'],
                                df['id'].to_numpy())
        np.testing.assert_equal(self.sne_truth_writer.sne_df['galaxy_id'],
                                df['host_galaxy'].to_numpy())
        np.testing.assert_equal(self.sne_truth_writer.sne_df['snra_in'],
                                df['ra'].to_numpy())
        np.testing.assert_equal(self.sne_truth_writer.sne_df['t0_in'],
                                df['t0'].to_numpy())
        np.testing.assert_equal(self.sne_truth_writer.sne_df['z_in'],
                                df['redshift'].to_numpy())

    def test_variability_truth(self):
        """
        Test some expected values for a SNIa in the test SNe catalog
        using a small opsim db table.
        """
        opsim_db_file = os.path.join(self.data_dir,
                                     'minion_1016_desc_dithered_v4_small.db')
        self.sne_truth_writer.write_variability_truth(opsim_db_file,
                                                      max_rows=60)
        with sqlite3.connect(self.outfile) as conn:
            df = pd.read_sql('select * from sn_variability_truth', conn)
        my_object = 'MS_10195_1375'
        self.assertIn(my_object, df['id'].to_list())
        my_df = df.query(f'id == "{my_object}"')
        for visit in (1425850, 1433860, 1495410):
            self.assertIn(visit, my_df['obsHistID'].to_list())


if __name__ == '__main__':
    unittest.main()
