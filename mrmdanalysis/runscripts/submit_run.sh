#!/bin/bash
##SBATCH --mail-user=julianhille@zedat.fu-berlin.de
#SBATCH --qos=standard                         
#SBATCH --job-name=MRMD
#SBATCH --ntasks=4
#SBATCH --time=12:00:00
#SBATCH --mem=8GB    
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

# manage input
inpfile=thermoForce_2025_03_19_22443935 
inpdir=/home/julianhille/project/mrmd_simulations/$inpfile
# manage output
outfile=tracerProduction_2025_03_20_$SLURM_JOB_ID
outdir=/home/julianhille/project/mrmd_simulations/$outfile
mkdir $outdir


# equilibrate 
#./eqNVTBerendsen
#./eqNVTLangevin

# thermodynamic force
#./ThermoForceIterations --nsteps 40000001 --outint 1000000 --sampling 200 --update 1000000 --outfile $outfile --binwidth 0.15 --damping 1 --neighbors 0 --forcemod 1 > $outfile.out

# production run
./TracerProductionNVT --nsteps 40000001 --outint 100000 --outfile $outfile --inpfile $inpdir/$inpfile > $outfile.out

# move output to directory
scontrol write batch_script $SLURM_JOB_ID slurm-$SLURM_JOB_ID.sh
mv slurm-$SLURM_JOB_ID.sh $outdir
mv $outfile* $outdir
mv slurm-$SLURM_JOB_ID.out $outdir