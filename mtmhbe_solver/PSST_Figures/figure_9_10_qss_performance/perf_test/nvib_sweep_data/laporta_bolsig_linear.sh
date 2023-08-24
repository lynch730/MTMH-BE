#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=01:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=18   # 36 processor core(s) per node 
#SBATCH --mem=8G   # maximum memory per node
#SBATCH --job-name="lap_bol_lin"
#SBATCH --output="laporta_run_bolsig_linear.out"

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load matlab/R2022a
module load intel
export omp_num_threads=8
cd /home/lynch/boltzmann_solvers/mtmhbe_solver/performance/laporta/perf_test/nvib_sweep_data/
matlab -nodisplay -nosplash -nodesktop -r "laporta_bolsig_linear;exit;"