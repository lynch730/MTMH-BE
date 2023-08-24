#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=04:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=18   # 36 processor core(s) per node 
#SBATCH --mem=64G   # maximum memory per node
#SBATCH --job-name="gst_mat_0log"
#SBATCH --output="laporta_run_mtmhbe_ET0_log_grid_study.out"

module load matlab/R2022a
cd /home/lynch/boltzmann_solvers/mtmhbe_solver/performance/laporta/grid_study/grid_sweep_data/
matlab -nodisplay -nosplash -nodesktop -r "laporta_mtmhbe_ET0_log;exit;"

