import os
import multiprocessing
import numpy as np
import desc.sims_truthcatalog as stc

dc2_info = '/global/cfs/cdirs/descssim/DC2'
agn_db_file = os.path.join(dc2_info, 'Run3.0i', 'agn_cosmoDC2_v1.1.4_ddf.db')
opsim_db_file = os.path.join(dc2_info, 'minion_1016_desc_dithered_v4_trimmed.db')
outfile = 'agn_truth_cat.db'
verbose = True

agn_truth_writer = stc.AGNTruthWriter(outfile, agn_db_file)
#agn_truth_writer.write(verbose=verbose)
#agn_truth_writer.write_auxiliary_truth(verbose=verbose)

num_agn = 11441
processes = 32
#num_agn = 20
#processes = 5

object_ranges = [int(_) for _ in np.linspace(0, num_agn, processes + 1)]
with multiprocessing.Pool(processes=processes) as pool:
    workers = []
    for xmin, xmax in zip(object_ranges[:-1], object_ranges[1:]):
        outfile = f'agn_truth_cat_{xmin:05d}_{xmax:05d}.db'
        func = stc.write_agn_variability_truth
        args = (agn_db_file, agn_truth_writer.query, opsim_db_file,)
        kwds = dict(object_range=(xmin, xmax), verbose=verbose, outfile=outfile)
        workers.append(pool.apply_async(func, args, kwds))
    pool.close()
    pool.join()
    _ = [worker.get() for worker in workers]
