#!/bin/bash
#SBATCH -N 8
#SBATCH -n 384
#SBATCH -c 1 # specify 6 threads per process
#SBATCH -t 01:00:00
#SBATCH -p workq
#SBATCH -A hpc_proteus02o
#SBATCH -o o.out # optional, name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e e.err # optional, name of the stderr, using job and first node values
#SBATCH -J emi2023_3d_floating_tld

date

module purge
module load proteus/fct
module load intel/2021.5.0
module load mvapich2/2.3.7/intel-2021.5.0
module load gcc/11.2.0
export LD_LIBRARY_PATH=/home/packages/compilers/intel/compiler/2022.0.2/linux/compiler/lib/intel64_lin:${LD_LIBRARY_PATH}
export MV2_HOMOGENEOUS_CLUSTER=1

mkdir -p $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cd $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cp $SLURM_SUBMIT_DIR/*.stl .
cp $SLURM_SUBMIT_DIR/*.py .
cp $SLURM_SUBMIT_DIR/*.sh .


#srun parun beji_battjes_so.py -F -l 5 -C "he=0.001 T=300.0" -O petsc.options.superlu_dist
srun parun TN_with_box_so.py -F -l 5 -C "he=1. T=3."  #-O petsc.options.asm

date

exit 0

