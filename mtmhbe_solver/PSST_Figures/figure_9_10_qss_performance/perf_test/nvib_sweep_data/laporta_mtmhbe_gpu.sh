#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=00:20:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=8   # 36 processor core(s) per node 
#SBATCH --mem=8G   # maximum memory per node
#SBATCH --gres=gpu:a100:1
#SBATCH --partition=gpu
#SBATCH --job-name="boltz_test_gpu"
#SBATCH --output="laporta_run_mtmhbe_ET0_log_gpu.out"

module load matlab/R2022a
export omp_num_threads=8
cd /home/lynch/boltzmann_solvers/mtmhbe_solver/performance/laporta/perf_test/nvib_sweep_data/
matlab -nodisplay -nosplash -nodesktop -r "laporta_mtmhbe_ET2_log_gpu;exit;"