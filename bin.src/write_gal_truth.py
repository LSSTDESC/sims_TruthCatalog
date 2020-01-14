import argparse
import os
from desc.sims_truthcatalog.galaxy_truth import GalaxyTruthWriter

'''
Interface to module which writes a truth catalog or part of a truth
catalog for specified healpixel to an sqlite3 file.
For usage type
  $ python write_gal_truth.py --help

'''

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Write a truth catalog (sqlite) for galaxies from the specified healpixel')
    parser.add_argument('output_dir', type=str,
                        help='Destination directory for output sqlite file')
    parser.add_argument('--healpixel', '--hpid', dest='hpid',
                        type=int, default=None,
                        help='healpixel for which entries will be written')
    parser.add_argument('--sed-fit-dir', type=str,
                        default='/global/projecta/projectdirs/lsst/groups/SSim/DC2/cosmoDC2_v1.1.4/sedLookup',
                        help='where to find sed fit files')
    parser.add_argument('--magnitude-cut', type=float, default=29.0,
                        help='only galaxies at least this bright will be include')
    parser.add_argument('--chunk-size', type=int, default=100000,
                        help='#galaxies to process in one go')
    parser.add_argument('--start', type=int, default=0,
                        help='galaxy to start with, as ordered in sed fit file')
    parser.add_argument('--nchunk', type=int, default=None,
                        help='if not None, stop after NCHUNK chunks')
    parser.add_argument('--dry-run', action='store_true',
                        help='no computation or database write')
    parser.add_argument('--max-parallel', dest='parallel', type=int, default=10,
                        help='Max number of chunks allowed to run concurrently')
    parser.add_argument('--call', action='store_true',
                        help='call _process_chunk rather than forking')
    parser.add_argument('--knl', action='store_true',
                        help='Use this flag for knl node; default is haswell')
    args = parser.parse_args()

    sed_fit_dir = args.sed_fit_dir
    assert os.path.isdir(sed_fit_dir)
    
    writer = GalaxyTruthWriter(args.output_dir, args.hpid, sed_fit_dir,
                               args.magnitude_cut, args.chunk_size,
                               args.start,args.nchunk, args.parallel,
                               args.dry_run, args.call, args.knl)
    writer.write()
