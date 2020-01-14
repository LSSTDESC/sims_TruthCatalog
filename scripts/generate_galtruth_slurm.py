import os
import sys
import numpy as np
import argparse
import os

if __name__ == "__main__":
    '''
    Write a script which may be submitted to slurm which runs jobs on
    potentially many nodes.  Each node is responsible for at most
    one healpixel
    '''

    our_image='docker:lsstdesc/stack-sims:w_2019_42-sims_w_2019_42-v2\n'
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_script', type=str,
                        default='local/slurm_scripts/batch_gal_truth.sl')
    parser.add_argument('--out_dir', type=str, default=None,
                        help='Dest dir for jobs run by out_script')
    parser.add_argument('--jobname', type=str,default='testjob',
                        help='used to identify batch job')
    parser.add_argument('--healpixel_file', type=str,
                        default='data/healpixels.txt',
                        help='File listing the healpixels to simulate')
    parser.add_argument('--chunk_size', type=int, default=100000)
    parser.add_argument('--max_parallel', type=int, default=28,
                        help='max number of processes to start per node (default 20 for knl, 28 for haswell)')
    parser.add_argument('--knl',  action='store_true',
                        help='by default use haswell')
    parser.add_argument('--debug', action='store_true',
                        help='Use debug queue and small chunk size')
                       
    args = parser.parse_args()

    print('generate_galtruth_slurm will use arguments ')
    print('out_script={}'.format(args.out_script))
    print('out_dir={}'.format(args.out_dir))
    print('jobname={}'.format(args.jobname))
    print('healpixl_file={}'.format(args.healpixel_file))
    print('chunk_size={}'.format(args.chunk_size))
    print('max_parallel={}'.format(args.max_parallel))

    pkg_root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    print('Package root is ', pkg_root)

    local_out = os.path.join(os.path.join(pkg_root, 'local'),'slurm_logs')
    slurm_out = os.path.join(local_out, '{}_out.txt'.format(args.jobname))
    slurm_err = os.path.join(local_out, '{}_err.txt'.format(args.jobname))
    out_dir = args.out_dir
    per_node_script = os.path.join(os.path.join(pkg_root, 'scripts'),
                                   'runshift_galtruth.sh')
    out_script_dir = os.path.dirname(args.out_script)
    os.makedirs(out_script_dir, exist_ok=True)

    # Read in the healpixels
    hp_path = args.healpixel_file
    hp = []
    with open(hp_path, 'r') as hp_file:
        for line in hp_file:
            hp.append(int(line))

    n_hp = len(hp)
    if not args.knl:
        n_hrs = 6          # should be conservative
        node_type = 'haswell'
    else:
        node_type = 'knl'
        print('--knl not supported')
        exit(1)

    # If directory to contain output script doesn't exist, make it
        out_file.write("if [ ! -d ${local_out} ]; then\n")
        out_file.write("    mkdir -p ${local_out}\n")
        out_file.write("fi\n")
        out_file.write('\n')
    
    ##exit(0)
    #Now create a script which starts up a job on one node for each healpixel
    #in the list.    Start with a bunch of instructions to SBATCH

    with open(args.out_script, 'w') as out_file:
        out_file.write('#!/bin/bash -l\n')
        out_file.write('#SBATCH ')
        out_file.write('--image={}'.format(our_image))
        out_file.write('#SBATCH -N {}\n'.format(n_hp))
        if (args.debug):
            out_file.write('#SBATCH -t 0:30:00\n'.format(n_hrs))
        else:
            out_file.write('#SBATCH -t {}:00:00\n'.format(n_hrs))
        out_file.write('#SBATCH -o {}\n'.format(slurm_out))
        out_file.write('#SBATCH -e {}\n'.format(slurm_err))
        if (args.debug):
            out_file.write('#SBATCH -q debug\n')
        else:
            out_file.write('#SBATCH -q regular\n')
        out_file.write('#SBATCH -A m1727\n')
        out_file.write('#SBATCH --tasks-per-node=1\n')
        out_file.write('#SBATCH --cpus-per-task=64\n')
        out_file.write('#SBATCH -C {}\n'.format(node_type))
        out_file.write('#SBATCH --job-name {}\n'.format(args.jobname))

        # Define some symbols
        out_file.write('per_node_script={}\n'.format(per_node_script))
        out_file.write('out_dir={}\n'.format(args.out_dir))
        ##out_file.write('local_out={}\n'.format(local_out))

        # Make output directory if absent
        out_file.write("if [ ! -d ${out_dir} ]; then\n")
        out_file.write("    mkdir -p ${out_dir}\n")
        out_file.write("fi\n")
        out_file.write('\n')


        # For each healpixel start up a shell script on a node which
        # does various kinds of setup, then invokes write_gal_truth.py
        # for that healpixel
        for i_hp in range(0, n_hp):
            out_file.write('\n')
            out_file.write('srun -N 1 -n 1 -c 64 --exclusive \\\n')
            out_file.write('shifter ${per_node_script} \\\n')
            out_file.write('  ${out_dir} ')
            out_file.write('{} '.format(args.max_parallel))
            out_file.write('{} '.format(args.chunk_size))
            out_file.write('{} '.format(hp[i_hp]))
            out_file.write('{} &\n\n'.format(pkg_root))

        out_file.write('\nwait\n')
        out_file.write("\necho 'master all done for {} '\n".format(args.out_dir))
        out_file.write('date\n')

        
