#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=04:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=18   # 36 processor core(s) per node 
#SBATCH --mem=8G   # maximum memory per node
#SBATCH --job-name="multibolt_laporta"
#SBATCH --output="laporta_run_multibolt_linear_800.out"

export omp_num_threads=8
/home/lynch/boltzmann_solvers/multibolt_timings/applications/laporta_run
./laporta_run_800_linear > '/home/lynch/boltzmann_solvers/mtmhbe_solver/performance/laporta/perf_test/nvib_sweep_data/laporta_multibolt_linear_800.txt'
