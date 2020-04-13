#!/usr/bin/env python
"""
Script to generate stellar_variability table using multiprocessing
on a single 32 core Cori-Haswell node.
"""
import numpy as np
import multiprocessing
from desc.sims_truthcatalog import write_star_variability_stats


stars_db_file = ('/global/projecta/projectdirs/lsst/'
                 'groups/SSim/DC2/dc2_stellar_healpixel.db')
row_offset = 0
processes = 32
num_rows = 20269973   # total number of rows in stars_db_file
chunk_size = 10000
row_bounds = list(range(0, num_rows, num_rows//processes))
row_bounds.append(num_rows)
row_bounds = np.array(row_bounds) + row_offset
workers = []
with multiprocessing.Pool(processes=processes) as pool:
    for row_min, row_max in zip(row_bounds[:-1], row_bounds[1:]):
        outfile = f'star_lc_stats_{row_min:08d}_{row_max:08d}.db'
        workers.append(pool.apply_async(write_star_variability_stats,
                                        (stars_db_file, outfile, row_min,
                                         row_max), dict(chunk_size=chunk_size)))
    pool.close()
    pool.join()
    [_.get() for _ in workers]
