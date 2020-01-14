#!/bin/bash

out_dir="${@:1:1}"
max_parallel="${@:2:1}"
chunk_size="${@:3:1}"
hp_id="${@:4:1}"
sims_TruthCatalog_root="${@:5:1}"
###knl_option="${@:6:1}"

export OMP_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

source /opt/lsst/software/stack/loadLSST.bash
setup  lsst_distrib                 # is this needed?
setup  lsst_sims
setup -j -r ${sims_TruthCatalog_root}

# Assume access to GCR has been set up, e.g. by pip install --user

## setup -j -r -t DC2production throughputs   #already in the image

python -W'ignore' \      # ignore all warnings
$SIMS_TRUTHCATALOG_DIR/bin.src/write_gal_truth.py \
${out_dir} \
--healpixel ${hp_id} \
--max-parallel ${max_parallel} \
--chunk-size ${chunk_size} \


