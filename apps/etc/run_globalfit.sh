#!/bin/bash
#SBATCH --job-name=LDC2a_Prod
#SBATCH --output=slurm_12.out
#SBATCH --error=slurm_12.err
#SBATCH --partition=hpc-ondemand
#SBATCH --ntasks=624
#SBATCH --tasks-per-node=8
#SBATCH --cpus-per-task=12

globalfit="${GLASS_PREFIX}/global_fit"
data="~/ldc/sangria/data/LDC2_sangria_training_v2.h5"
vgb="${GLASS_REPO}/apps/etc/hmcnc.dat"
outdir="${GLASS_TEST}/globalfit"

#path/to/search_sources.dat (contains starting point for MBH sampler)
mbh="~/ldasoft/apps/etc/"

fmin=0.0003
samples=128
samples_max=128


Tobs=7864320
padding=32

#Tobs=31457280
#padding=128
#outdir="globalfit"

Tstart=0
sources=40


cmd="${globalfit} --h5-data ${data} --sangria --fmin ${fmin} --chains 12 --threads 12 --start-time ${Tstart} --duration ${Tobs} --samples ${samples} --padding ${padding} --sources ${sources} --rundir ${outdir} --known-sources ${vgb} "


#--mbh-search-path ${mbh}"


echo $cmd
export OMP_NUM_THREADS=12

mpirun -np 4 --oversubscribe $cmd
