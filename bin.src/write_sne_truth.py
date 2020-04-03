#!/usr/bin/env python
"""
Script to write the SNe truth catalog.
"""
import argparse
import os
import sys
import desc.sims_truthcatalog
from datetime import datetime as dt

parser = argparse.ArgumentParser(description="Write a truth catalog (sqlite) for SNe given a SNe parameters db file")
parser.add_argument('outfile', type=str,
                    help='Filename for output sqlite file')
parser.add_argument('sne_db_file', type=str,
                    help='Filename for input SNe parameters db file')
parser.add_argument('--dry-run', action='store_true',
                    help='no computation or database write')
parser.add_argument('--variable-table', action='store_true',
                    help='create and fill variable truth table')
parser.add_argument('--aux-table', action='store_true',
                    help='create and fill auxilliary truth table')
parser.add_argument('--max-parallel', default=1, dest='max_parallel', type=int,
                    help='max # of processes to be run in parallel; applies only to variability table (default: %(default)s)')
parser.add_argument('--row-limit', default=None, type=int,
                    help='write at most this number per process (applies only to variability table (default: no limit)')
parser.add_argument('--no-summary', action='store_true', help='speed up debugging other tables')
parser.add_argument('--verbose', action='store_true',
                    help='write debug output while making aux table')
parser.add_argument('--interval-file', default=None,
                    help='to process subset of SNe specify yaml file of intervals')
parser.add_argument('--chunk-log', default='chunk_+.log',
                    help='Output file per chunk for variable table. +, if present, is replace with chunk number')

args = parser.parse_args()

desc.sims_truthcatalog.print_callinfo(sys.argv[0], args)
time_fmt = desc.sims_truthcatalog.TIME_TO_SECOND_FMT
writer = desc.sims_truthcatalog.SNeTruthWriter(args.outfile, args.sne_db_file,
                                               max_parallel=args.max_parallel,
                                               dry_run=args.dry_run,
                                               no_summary=args.no_summary)
print('{}   Created SNeTruthWriter instance'.format(dt.now().strftime(time_fmt)))
if not args.no_summary:
    writer.write()
    print('{}   Summary step table done'.format(dt.now().strftime(time_fmt)))
    sys.stdout.flush()
    
if args.aux_table:
    writer.write_auxiliary_truth()
    print('{}   Aux table step done'.format(dt.now().strftime(time_fmt)))
    sys.stdout.flush()

if args.variable_table:
    opsim_db_dir = '/global/projecta/projectdirs/lsst/groups/SSim/DC2'
    writer.write_variability_truth(os.path.join(opsim_db_dir, 
                                                'minion_1016_desc_dithered_v4_sfd.db'),
                                   chunk_log=args.chunk_log,
                                   max_rows=args.row_limit,
                                   max_parallel=args.max_parallel,
                                   interval_file = args.interval_file,
                                   verbose=args.verbose)
    print('{}   Variability step done'.format(dt.now().strftime(time_fmt)))

print('{}   write_sne_truth done'.format(dt.now().strftime(time_fmt)))
