import os
import multiprocessing
import numpy as np
import desc.sims_truthcatalog as stc

dc2_info = '/global/cfs/cdirs/descssim/DC2'
lensed_agn_truth_cat = os.path.join(dc2_info, 'Run3.0i', 'truth_tables',
                                    'updated_lensed_agn_truth.db')
opsim_db_file = os.path.join(dc2_info, 'minion_1016_desc_dithered_v4_trimmed.db')
verbose = True

outfile = 'lensed_agn_truth_cat.db'
stc.write_lensed_agn_truth_summary(lensed_agn_truth_cat, outfile, verbose=verbose)

num_agn = 2437
processes = 32
#num_agn = 21
#processes = 3

object_ranges = [int(_) for _ in np.linspace(0, num_agn, processes + 1)]
with multiprocessing.Pool(processes=processes) as pool:
    workers = []
    for xmin, xmax in zip(object_ranges[:-1], object_ranges[1:]):
        outfile = f'lensed_agn_variability_truth_{xmin:04d}_{xmax:04d}.db'
        func = stc.write_lensed_agn_variability_truth
        args = (opsim_db_file, lensed_agn_truth_cat, outfile)
        kwds = dict(verbose=verbose, object_range=(xmin, xmax))
        workers.append(pool.apply_async(func, args, kwds))
    pool.close()
    pool.join()
    _ = [worker.get() for worker in workers]
