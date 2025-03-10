#!/bin/bash
##SBATCH --mail-user=julianhille@zedat.fu-berlin.de
#SBATCH --qos=standard                         
#SBATCH --job-name=MRMD
#SBATCH --ntasks=8
#SBATCH --time=01:00:00
#SBATCH --mem=10GB    
#SBATCH --nodes=1-1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

# load modules into shell environment
module purge
module load foss/2023a CMake/3.26.3-GCCcore-12.3.0  make/4.4.1-GCCcore-12.3.0 Anaconda3/2022.05 HDF5/1.14.0-gompi-2023a CUDA/12.4.0 # toolchain for MRMD with OpenMP, MPI and CUDA support
conda activate MRMDvenv

# OpenMP thread-affinity
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

# configure
cmake -S .. -B ./  -DCMAKE_BUILD_TYPE=Debug -DMRMD_ENABLE_COVERAGE=ON -DMRMD_ENABLE_TESTING=ON -DMRMD_ENABLE_HDF5=ON -DMRMD_ENABLE_MPI=ON -DMRMD_ENABLE_PYTHON=ON -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON -DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_PASCAL61=ON -DKokkos_ENABLE_CUDA_CONSTEXPR=ON &> cmake.txt

# build
cmake --build ./ --target all -j $SLURM_NTASKS &> cmake_all.txt  

# test
cmake --build ./tests/  --target test -j $SLURM_NTASKS &> cmake_test.txt
make test -j $SLURM_NTASKS &> make_test.txt 
examples/Argon/Argon --nsteps 10 