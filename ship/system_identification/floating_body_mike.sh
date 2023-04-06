#!/bin/bash
##SBATCH -N 4
##SBATCH -n 256
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 00:10:00
#SBATCH -p workq
#SBATCH -A hpc_proteus02o
#SBATCH -o o.out
#SBATCH -e e.err
#SBATCH -J emi2023_free_decay_ship
#load proteus module and ensure proteus's python is in path

date

module purge
module load proteus/fct
export LD_LIBRARY_PATH=/home/packages/compilers/intel/compiler/2022.0.2/linux/compiler/lib/intel64_lin:${LD_LIBRARY_PATH}
export MV2_HOMOGENEOUS_CLUSTER=1

mkdir -p $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cd $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cp $SLURM_SUBMIT_DIR/*.py .
cp $SLURM_SUBMIT_DIR/*.sh .

parun --TwoPhaseFlow pmtld.py -l 5 -C "he=0.5"

exit 0