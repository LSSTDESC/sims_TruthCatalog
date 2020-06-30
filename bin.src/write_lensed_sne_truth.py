#!/usr/bin/env python
import argparse
import desc.sims_truthcatalog as stc

parser = argparse.ArgumentParser(
    description='Write truth tables for lensed SNe for Run3.0i')
parser.add_argument('--opsim_db_file', type=str, help='OpSim db file for DC2',
                    default=('/global/cfs/cdirs/descssim/DC2/',
                             'minion_1016_desc_dithered_v4_trimmed.db'))
parser.add_argument('--lensed_sne_truth_cat', type=str,
                    help='sqlite3 db file containing lensed SNe parameters',
                    default=('/global/cfs/cdirs/descssim/DC2/Run3.0i/'
                             'truth_tables/updated_lensed_sne_truth.db'))
parser.add_argument('--outfile', type=str,
                    help='Filename of output sqlite3 file',
                    default='lensed_sne_truth_cat.db')
args = parser.parse_args()

stc.write_lensed_sn_variability(args.opsim_db_file, args.lensed_sne_truth_cat,
                                args.outfile)



