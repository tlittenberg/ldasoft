#!/bin/bash
#SBATCH --job-name=LDC2a_Prod
#SBATCH --output=slurm_12.out
#SBATCH --error=slurm_12.err
#SBATCH --partition=hpc-ondemand
#SBATCH --ntasks=624
#SBATCH --tasks-per-node=8
#SBATCH --cpus-per-task=12

globalfit=/shared/opt/bin/global_fit
data=/benchmarks/ldc/sangria/LDC2_sangria_training_v2.h5
vgb=/benchmarks/ldc/sangria/ldc_sangria_vgb_list.dat

#path/to/search_sources.dat (contains starting point for MBH sampler)
mbh=/benchmarks/ldc/sangria

fmin=0.0003
samples=128
samples_max=128

#Tobs=3932160
#padding=16
#outdir=/shared/ldc/sangria/prod/sangria_training_01mo

Tobs=7864320
padding=32
#ucb=/benchmarks/ldc/sangria/sangria_training_01mo/gb_catalog.cache
outdir=/shared/ldc/sangria/prod/sangria_training_03mo

#Tobs=15728640
#padding=64
#ucb=/benchmarks/ldc/sangria/sangria_training_03mo/gb_catalog.cache
#outdir=/shared/ldc/sangria/prod/sangria_training_06mo

#Tobs=31457280
#padding=128
#ucb=/benchmarks/ldc/sangria/sangria_training_06mo/gb_catalog.cache
#outdir=/benchmarks/ldc/sangria/sangria_training_12mo

Tstart=0
sources=40

#Set up whatever package we need to run with
module load gsl-2.7.1-gcc-9.4.0-ylwugg3
module load hdf5-1.12.2-gcc-9.4.0-f5p3zoy

cmd="${globalfit} \
--h5-data ${data} \
--sangria \
--fmin ${fmin} \
--chains 24 \
--start-time ${Tstart} \
--duration ${Tobs} \
--samples ${samples} \
--padding ${padding} \
--sources ${sources} \
--rundir ${outdir} \
--mbh-search-path ${mbh} \
--known-sources ${vgb} \
"
#--catalog ${ucb} \


echo $cmd
export OMP_NUM_THREADS=18

mpirun -np $SLURM_NTASKS $cmd
#(time mpirun -np $SLURM_NTASKS $cmd) 1> run.out 2>&1
