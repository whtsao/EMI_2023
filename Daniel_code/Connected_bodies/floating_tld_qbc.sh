#!/bin/bash
#SBATCH -N 8
#SBATCH -n 384
#SBATCH -c 1 # specify 6 threads per process
#SBATCH -t 01:00:00
#SBATCH -p workq
#SBATCH -A loni_proteus01s
#SBATCH -o o.out # optional, name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e e.err # optional, name of the stderr, using job and first node values
#SBATCH -J emi2023_3d_floating_tld

date

module purge
#pip install pycatenary
#pip install py2gmsh
module load proteus/1.8.1
#pip install pycatenary
#pip install py2gmsh

mkdir -p $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cd $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cp $SLURM_SUBMIT_DIR/*.stl .
cp $SLURM_SUBMIT_DIR/*.py .
cp $SLURM_SUBMIT_DIR/*.sh .


#srun parun beji_battjes_so.py -F -l 5 -C "he=0.001 T=300.0" -O petsc.options.superlu_dist
srun parun TN_with_box_so.py -F -l 5 -C "he=1. T=1."  #-O petsc.options.asm

date

exit 0

