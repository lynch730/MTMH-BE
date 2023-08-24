#!/bin/bash

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 36 processor core(s) per node 
#SBATCH --mem=64G   # maximum memory per node
#SBATCH --partition=community    # community node(s)
#SBATCH --job-name="sweep_test"
#SBATCH --output="NU_sweep.out"

module load matlab/R2022a
cd /home/lynch/boltzmann_solvers/mtmhbe_solver/performance/pulse_loki_b/
matlab -nodisplay -nosplash -nodesktop -r "loki_b_pulse_Neps;exit;"
