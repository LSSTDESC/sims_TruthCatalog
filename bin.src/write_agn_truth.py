#!/usr/bin/env python
import os
import argparse
import desc.sims_truthcatalog as stc

parser = argparse.ArgumentParser(
    description='Write AGN truth catalogs for Run3.0i')
parser.add_argument('--agn_db_file', type=str, help='AGN db file',
                    default=('/global/cfs/cdirs/descssim/DC2/Run3.0i/'
                             'agn_cosmoDC2_v1.1.4_ddf.db'))
parser.add_argument('--opsim_db_file', type=str, help='OpSim db file',
                    default=('/global/cfs/cdirs/descssim/DC2/'
                             'minion_1016_desc_dithered_v4_trimmed.db'))
parser.add_argument('--outfile', type=str, help='output sqlite3 filename',
                    default='agn_truth_cat.db')
parser.add_argument('--start_mjd', type=float, default=59580,
                    help=('Starting MJD of variability data. The default is'
                          'the start of the minion_1016 candence.'))
parser.add_argument('--end_mjd', type=float, default=61395,
                    help=('Ending MJD of variability data.  The default is '
                          'the end of Y05 of the minion_1016 cadence.'))
parser.add_argument('--verbose', action='store_true',
                    default=False, help='enable verbose output')
parser.add_argument('--num_objects', type=int, default=None,
                    help=('Number of objects to process. '
                          'If None, then process all objects.'))
parser.add_argument('--overwrite', default=False, action='store_true',
                    help='Overwrite existing output file.')
args = parser.parse_args()

if os.path.isfile(args.outfile):
    if args.overwrite:
        os.remove(args.outfile)
    else:
        raise OSError(f'{args.outfile} already exists.')

agn_truth_writer = stc.AGNTruthWriter(args.outfile, args.agn_db_file)

agn_truth_writer.write(verbose=args.verbose)

agn_truth_writer.write_auxiliary_truth(verbose=args.verbose)

agn_truth_writer.write_variability_truth(args.opsim_db_file,
                                         start_mjd=args.start_mjd,
                                         end_mjd=args.end_mjd,
                                         num_objects=args.num_objects,
                                         verbose=args.verbose)
