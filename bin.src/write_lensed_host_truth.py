#!/usr/bin/env python
import os
import argparse
import desc.sims_truthcatalog as stc

parser = argparse.ArgumentParser(
    description='Write lensed host truth catalog for Run3.0i')
parser.add_argument('--host_truth_db_file', type=str,
                    help=('db file containing model parameters for lensed '
                          'AGN and SNe hosts.'),
                    default=('/global/cfs/cdirs/descssim/DC2/Run3.0i/'
                             'truth_tables/host_truth.db'))
parser.add_argument('--image_dir', type=str,
                    help=('directory containing the FITS stamps for the '
                          'lensed hosts'),
                    default=('/global/cfs/cdirs/descssim/DC2/Run3.0i/'
                             'FITS_stamps'))
parser.add_argument('--outfile', type=str, help='output sqlite3 filename',
                    default='lensed_host_truth.db')
parser.add_argument('--verbose', default=False, action='store_true',
                    help='Verbosity flag.')
parser.add_argument('--overwrite', default=False, action='store_true',
                    help='Overwrite existing output file')
args = parser.parse_args()

if os.path.isfile(args.outfile):
    if args.overwrite:
        os.remove(args.outfile)
    else:
        raise RuntimeError(f'{args.outfile} already exists')

stc.write_lensed_host_truth(args.host_truth_db_file, args.image_dir,
                            args.outfile, verbose=args.verbose)
