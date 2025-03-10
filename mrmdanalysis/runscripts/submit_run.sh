#!/bin/bash
##SBATCH --mail-user=julianhille@zedat.fu-berlin.de
#SBATCH --qos=standard                         
#SBATCH --job-name=MRMD
#SBATCH --ntasks=8
#SBATCH --time=00:05:00
#SBATCH --mem=10GB    
#SBATCH --nodes=1-1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

# load modules into shell environment
module purge
module load foss/2023a CMake/3.26.3-GCCcore-12.3.0  make/4.4.1-GCCcore-12.3.0 Anaconda3/2022.05 HDF5/1.14.0-gompi-2023a CUDA/12.4.0 # toolchain for MRMD with OpenMP, MPI and CUDA support

# OpenMP thread-affinity
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

# non-critical metadata
module list 2> modules_list.txt
hostname
date
pwd

# run 
./Argon --nsteps 10000