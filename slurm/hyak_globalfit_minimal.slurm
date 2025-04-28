#!/bin/bash
# All Sbatch stuff needs to be a continuous comment
#
#SBATCH --job-name=globalfit
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
#Kill us after 10h for now
#SBATCH --time=10:00:00

set -o xtrace

DATA_DIR="/gscratch/gwastro/mtauraso"

# TODO: Pack PACKAGE_DIR
# Depending on things like inode limits and f/s speed, it may be desirable to
# zip this directory and then unpack it on a node in /tmp somewhere.
# As written, hyak's backplane is responsible for moving many small files which
# are mostly our conda environment
PACKAGE_DIR="/mmfs1/home/mtauraso/ldasoft/build_conda"

LDASOFT_DIR="${PACKAGE_DIR}/ldasoft-install"


# input data from LDC
data="${DATA_DIR}/LDC2_sangria_training_v2.h5"

# verification binaries data file
vgb="${LDASOFT_DIR}/data/ldc_sangria_vgb_list.dat"
# path/to/search_sources.dat (contains starting point for MBH sampler)
mbh="${LDASOFT_DIR}/data/"
# working directory is where we look for ucb_frequency_spacing.dat
workingdir="${LDASOFT_DIR}/data"

fmin=0.0003
samples=128
samples_max=128

#Tobs=3932160
#padding=16
#outdir=/shared/ldc/sangria/prod/sangria_training_01mo

Tobs=7864320
padding=32
#ucb=/benchmarks/ldc/sangria/sangria_training_01mo/gb_catalog.cache
outdir="${DATA_DIR}/global_fit_LDC2a_training_01mo"

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
CONDA_ENV_DIR="${PACKAGE_DIR}/conda-env"
eval $(conda shell.bash activate ${CONDA_ENV_DIR})
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"

module load ompi

#export UCX_LOG_LEVEL="info"

cleanup() 
{
	echo "In cleanup function looking for core dumps"
	# List core dumps available
	coredumpctl list ${LDASOFT_DIR}/bin/global_fit

	# Dump the last coredump to ${outdir}/core-<epoch seconds>
	DATE=$(date +%s)
	coredumpctl -1 --output=${outdir}/core-$DATE.coredump dump

	# Now dump dmesg to the log
	dmesg | grep -i oom
}

trap cleanup SIGINT SIGABRT SIGQUIT SIGTERM

# Run the global fit command ensuring we have set the working directory properly
srun -n ${SLURM_NTASKS} --chdir ${workingdir} ${LDASOFT_DIR}/bin/${global_fit_cmd}


