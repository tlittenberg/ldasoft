#!/bin/bash
# All Sbatch stuff needs to be a continuous comment
#
#SBATCH --job-name=global_fit_LDC2a_dev
#SBATCH --mail-type=NONE
#SBATCH --mail-user=mtauraso
#SBATCH --account=gwastro
#SBATCH --partition=compute
#
# Note: new hyak nodes have 40 cores.
# Tyson's config has 8 tasks per node * 12 cpus per task = 96 cores
# Then something like 600 tasks are run... however 5 tasks is the minimum
#
# If we want to start 5 tasks * 12 cpus per task = 60 cores
# I'm going to reduce this to 4 cpus per task so we can get 
# 5*4 = 20 cores and run on half of the single node gwastro has
#
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=5
#SBATCH --mem-per-cpu=4G
#
#Don't move any environment variables into slurm-land
#SBATCH --export=ALL
#
#Kill us after a time
#SBATCH --time=10:00:00


# SLURM says we get them for free... with a job ID, which is better actually
# stdout and stderr logging (Slurm docs says we'll get these files anyway)
# SBATCH --output=slurm_12.out
# SBATCH --error=slurm_12.err

set -o xtrace


DATA_DIR="/gscratch/gwastro/mtauraso"
LDASOFT_DIR="/mmfs1/home/mtauraso/ldasoft/build-conda/ldasoft-install"

# TODO: Would like to have the build script build the environment in the build_conda dir
# That way we can build this all off of a single directory tree
#
# It may then be possible in the preamble to unpack miniconda and the env into /tmp
# Which if successful could give us the same perf as apptainer one-file strategy
CONDA_DIR="${LDASOFT_DIR}/../conda_env"
OMPI_VERSION="4.1.4"

# input data and verification galactic binary files
data="${DATA_DIR}/LDC2_sangria_training_v2.h5"
vgb="${LDASOFT_DIR}/data/ldc_sangria_vgb_list.dat"

# path/to/search_sources.dat (contains starting point for MBH sampler)
mbh="${LDASOFT_DIR}/data/"

fmin=0.0003
samples=128
samples_max=128

#Tobs=3932160
#padding=16
#outdir=/shared/ldc/sangria/prod/sangria_training_01mo

Tobs=7864320
padding=32
#ucb=/benchmarks/ldc/sangria/sangria_training_01mo/gb_catalog.cache
outdir="${DATA_DIR}/global_fit_LDC2a_dev"

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
#module load gsl-2.7.1-gcc-9.4.0-ylwugg3
#module load hdf5-1.12.2-gcc-9.4.0-f5p3zoy

global_fit_cmd="global_fit \
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



# Activate conda first, and then load modules
# This ensures that any hyak module stuff appears first in $PATH/$LD_LIBRARY_PATH
eval $(conda shell.bash activate ${CONDA_DIR})
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"

module load ompi/${OMPI_VERSION}


srun -n ${SLURM_NTASKS} --chdir ${mbh} ${LDASOFT_DIR}/bin/${global_fit_cmd}


