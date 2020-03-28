#!/usr/bin/env python
import os
import desc.sims_truthcatalog as stc
import multiprocessing

opsim_db_file = '/global/cscratch1/sd/jchiang8/desc/Run2.2i/minion_1016_desc_dithered_v4_trimmed.db'
assert os.path.isfile(opsim_db_file)

star_lc_stats_db_file = '/global/cscratch1/sd/jchiang8/desc/Run2.2i/stellar_variability/merged_star_db/star_lc_stats_merged_trimmed.db'
assert os.path.isfile(star_lc_stats_db_file)

stars_db_file = '/global/projecta/projectdirs/lsst/groups/SSim/DC2/dc2_stellar_healpixel.db'
assert os.path.isfile(stars_db_file)

num_rows = None
processes = 32

outfile_prefix = 'svt'

svt = stc.StellarVariabilityTruth(stars_db_file, opsim_db_file,
                                  star_lc_stats_db_file)
svt.write_tables(outfile_prefix, processes, num_rows=num_rows)
