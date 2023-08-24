#!/bin/bash
 
MATLAB_VER=R2021a
#FFLAGS="-O2 -fpp -qopenmp -qmkl=parallel -xHost -parallel -qopt-matmul -mtune=skylake-avx512  -warn all,noexternal"
#FFLAGS="-O2 -fpp -qopenmp -qmkl=parallel -warn all,noexternal"
FFLAGS="-O3 -fpp -qopenmp -qmkl=parallel -xHost"

cd /home/lynch/boltzmann_solvers/mtmhbe_solver/solvers_fortran/laporta_perf_test
rm *.o *.mod
module purge
module load intel
module load intel-mkl
module load matlab/$MATLAB_VER
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/shared/hpc/matlab/$MATLAB_VER/bin/glnxa64:/shared/hpc/matlab/$MATLAB_VER/sys/os/glnxa64
export PATH=$PATH:/shared/hpc/matlab/$MATLAB_VER/extern/include

ifort $FFLAGS -I/shared/hpc/matlab/$MATLAB_VER/extern/include -c ../source/timers.f90
ifort $FFLAGS -I/shared/hpc/matlab/$MATLAB_VER/extern/include -c ../source/mdat.f90
ifort $FFLAGS -I/shared/hpc/matlab/$MATLAB_VER/extern/include -c ../source/lsolve.f90
ifort $FFLAGS -I/shared/hpc/matlab/$MATLAB_VER/extern/include -c ../source/bz_solution.f90
ifort $FFLAGS -I/shared/hpc/matlab/$MATLAB_VER/extern/include -c ../source/mbsol.f90
ifort $FFLAGS -I/shared/hpc/matlab/$MATLAB_VER/extern/include -c ../source/solver_common.f90
ifort $FFLAGS -I/shared/hpc/matlab/$MATLAB_VER/extern/include -c ../source/solver_qss.f90
ifort $FFLAGS -I/shared/hpc/matlab/$MATLAB_VER/extern/include -c laporta_perf_test.f90
ifort *.o  $FFLAGS -L/shared/hpc/matlab/$MATLAB_VER/bin/glnxa64 -L/shared/hpc/matlab/$MATLAB_VER/sys/os/glnxa64 -o laporta_exe -lmx -lmat -lstdc++
