#!/bin/bash

#SBATCH --time=06:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 36 processor core(s) per node 
#SBATCH --mem=64G   # maximum memory per node
#SBATCH --partition=compute    # community node(s)
#SBATCH --job-name="temporal"
#SBATCH --output="temporal.out"

module load matlab/R2022a
cd /home/lynch/boltzmann_solvers/mtmhbe_solver/performance/pulse_loki_b/
matlab -nodisplay -nosplash -nodesktop -r "time_convergence;exit;"
