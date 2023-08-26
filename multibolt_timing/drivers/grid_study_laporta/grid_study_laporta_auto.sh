#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=06:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=18   # 36 processor core(s) per node 
#SBATCH --mem=8G   # maximum memory per node
#SBATCH --job-name="multibolt_laporta"
#SBATCH --output="laporta_run_multibolt_linear_100.out"

export omp_num_threads=8
cd /home/lynch/boltzmann_solvers/multibolt_timings/applications/grid_study_laporta
./grid_study_laporta_auto > '/home/lynch/boltzmann_solvers/mtmhbe_solver/performance/laporta/grid_study/grid_sweep_data/laporta_multibolt_auto.txt'
