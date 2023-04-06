#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -c 1 # specify 6 threads per process
#SBATCH -t 00:10:00
#SBATCH -p workq
#SBATCH -A loni_proteus01s
#SBATCH -o o.out # optional, name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e e.err # optional, name of the stderr, using job and first node values
#SBATCH -J emi2023_2d_floating_tld_fr100

date

module purge
#pip install pycatenary
#pip install py2gmsh
module load proteus/1.8.1
#pip install pycatenary
#pip install py2gmsh

mkdir -p $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cd $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
#cp $SLURM_SUBMIT_DIR/*.stl .
cp $SLURM_SUBMIT_DIR/petsc.options.superlu_dist .
cp $SLURM_SUBMIT_DIR/*.py .
cp $SLURM_SUBMIT_DIR/*.sh .


#srun parun beji_battjes_so.py -F -l 5 -C "he=0.001 T=300.0" -O petsc.options.superlu_dist
parun --TwoPhaseFlow two_floating_bodies.py -F -l 5 -C "he=0.5 fr=1." #-O petsc.options.superlu_dist #-O petsc.options.asm

date

exit 0

