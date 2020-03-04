#!/usr/bin/env python
'''
Interface to module which writes a truth catalog or part of a truth
catalog for some or all stars coming from a star db file.
For usage type
  $ python write_star_truth.py --help

'''
import argparse
import os
from desc.sims_truthcatalog.star_truth import StarTruthWriter


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Write a truth catalog (sqlite) for stars from the specified healpixel')
    parser.add_argument('outfile', type=str,
                        help='Destination path for output sqlite file')
    parser.add_argument('--star-db', type=str,
                        default='/global/projecta/projectdirs/lsst/groups/SSim/DC2/dc2_stellar_healpixel.db',
                        help='where to find star db (default: %(default)s)')
    parser.add_argument('--row-limit', dest='row_limit',
                        type=int, default=None,
                        help='max stars to handle per spawned process (default: no limit)')
    parser.add_argument('--max-parallel', dest='parallel', type=int, default=20,
                        help='Max number of chunks allowed to run concurrently (default: %(default)s)')
    parser.add_argument('--dry-run', action='store_true',
                        help='if specified, no computation or database write (default: False)')
    args = parser.parse_args()
    print('write_star_truth invoked with arguments')
    for e in dir(args):
        if not e.startswith('__'):
            nm = 'args.' + e
            print('{}: {}'.format(e, eval(nm)))

    writer = StarTruthWriter(args.outfile, args.star_db, 
                             row_limit=args.row_limit,
                             max_parallel=args.parallel,
                             dry_run=args.dry_run)
    writer.write()
