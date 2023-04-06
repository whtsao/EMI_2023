#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1 # specify 6 threads per process
#SBATCH -t 12:0:00
#SBATCH -p workq
#SBATCH -A loni_proteus01s
#SBATCH -o o.out # optional, name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e e.err
#SBATCH -J emi2023_free_decay_ship

date

module purge
module load proteus/1.8.1

mkdir -p $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cd $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cp $SLURM_SUBMIT_DIR/*.py .
#cp $SLURM_SUBMIT_DIR/petsc.options.superlu_dist .
#cp $SLURM_SUBMIT_DIR/petsc.options.asm .
cp $SLURM_SUBMIT_DIR/*.sh .

parun --TwoPhaseFlow pmtld.py -l 5 -C "he=0.1"
#srun parun --TwoPhaseFlow pmtld.py -F -l 5 -C "he=0.002 fr=1.0" -O petsc.options.asm

date

exit 0

