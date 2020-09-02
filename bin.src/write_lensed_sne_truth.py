#!/usr/bin/env python
import os
import argparse
import desc.sims_truthcatalog as stc

parser = argparse.ArgumentParser(
    description='Write truth tables for lensed SNe for Run3.1i')
parser.add_argument('--opsim_db_file', type=str, help='OpSim db file for DC2',
                    default=('/global/cfs/cdirs/descssim/DC2/'
                             'minion_1016_desc_dithered_v4_trimmed.db'))
parser.add_argument('--lensed_sne_truth_cat', type=str,
                    help='sqlite3 db file containing lensed SNe parameters',
                    default=('/global/cfs/cdirs/descssim/DC2/Run3.0i/'
                             'truth_tables/updated_lensed_sne_truth.db'))
parser.add_argument('--outfile', type=str,
                    help='Filename of output sqlite3 file',
                    default='lensed_sne_truth_cat.db')
parser.add_argument('--start_mjd', type=float, default=59580,
                    help=('Starting MJD of variability data. The default is'
                          'the start of the minion_1016 candence.'))
parser.add_argument('--end_mjd', type=float, default=61405,
                    help=('Ending MJD of variability data.  The default is '
                          'the end of Y05 of the minion_1016 cadence.'))
parser.add_argument('--verbose', default=False, action='store_true',
                    help='Verbosity flag.')
parser.add_argument('--num_objects', type=int, default=None,
                    help=('Number of objects to process. '
                          'If None, then process all objects.'))
parser.add_argument('--overwrite', default=False, action='store_true',
                    help='Flag to overwrite an existing output file.')
args = parser.parse_args()

if os.path.isfile(args.outfile):
    if args.overwrite:
        os.remove(args.outfile)
    else:
        raise RuntimeError(f'{args.outfile} already exists')

stc.write_lensed_sn_truth_summary(args.lensed_sne_truth_cat, args.outfile,
                                  verbose=args.verbose)

stc.write_lensed_sn_variability_truth(args.opsim_db_file,
                                      args.lensed_sne_truth_cat,
                                      args.outfile, verbose=args.verbose,
                                      num_objects=args.num_objects)
