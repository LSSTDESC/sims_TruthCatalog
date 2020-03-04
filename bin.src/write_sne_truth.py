#!/usr/bin/env python
"""
Script to write the SNe truth catalog.
"""
import argparse
import os
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
parser.add_argument('--max-parallel', default=1, dest='parallel', type=int,
                    help='max # of processes to be run in parallel (default: %(default)s)')
parser.add_argument('--sne-limit', default=None, type=int,
                    help='process at most this number (default: no limit)')

args = parser.parse_args()

print('write_sne_truth invoked with arguments')
for e in dir(args):
    if not e.startswith('_'):
        nm = 'args.' + e
        print('{}: {}'.format(e, eval(nm)))

print(dt.now())
writer = desc.sims_truthcatalog.SNeTruthWriter(args.outfile, args.sne_db_file,
                                               sne_limit=args.sne_limit,
                                               dry_run=args.dry_run)
print('Created instance at ', dt.now()) 
writer.write()
print('Wrote summary table at ', dt.now()) 
if args.aux_table:
    writer.write_auxiliary_truth()
    print('Wrote aux table at ', dt.now()) 

if args.variable_table:
    opsim_db_dir = '/global/projecta/projectdirs/lsst/groups/SSim/DC2'
    writer.write_variability_truth(os.path.join(opsim_db_dir,
                                                'minion_1016_desc_dithered_v4_sfd.db'),
                                   sne_limit=args.sne_limit, max_parallel=args.max_parallel)
    print('Wrote variability table at ', dt.now())
