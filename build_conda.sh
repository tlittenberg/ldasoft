#!/bin/bash
set -o errexit
set -o pipefail
set -o nounset

# Build LDASOFT and MBH from scratch in a conda environment.
# Pack the environment and exes into a single folder structure

# Get the directory of the current script
SOURCE_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# We use build-conda as the staging area
BUILD_DIR="${SOURCE_DIR}/build_conda"

rm -rf ${BUILD_DIR}
mkdir -p ${BUILD_DIR}
pushd ${BUILD_DIR}

	# Build a brand-new conda env.
	CONDA_ENV_YML="${SOURCE_DIR}/environment.yml"
	CONDA_ENV_DIR="${BUILD_DIR}/conda-env"
	conda env create -f ${CONDA_ENV_YML} -p ${CONDA_ENV_DIR}

	# Activate our conda env
	eval $(conda shell.bash activate ${CONDA_ENV_DIR})

	# Tell cmake where to find conda environment
	export CMAKE_PREFIX_PATH="${CONDA_PREFIX}:${CMAKE_PREFIX_PATH:-}"

	# Module load our stuff from hyak
	module load ompi # Note this also loads gcc and ucx
	module load cmake

	# Hyaks load of gcc/11.2.x sets $CC to mpiCC which is technically a c++ compiler
	# cmake gets mad about this so we set CC to a c compiler (mpicc) in order to keep
	# cmake happy.
	export CC="mpicc"

	# Build and install MBH from source
	#
	# TODO: This is pointed mtauraso's fork of MBH because some patches are not yet merged.
	# build-fixup branch includes these two PRs:
	#
	# https://github.com/eXtremeGravityInstitute/LISA-Massive-Black-Hole/pull/5
	# https://github.com/eXtremeGravityInstitute/LISA-Massive-Black-Hole/pull/4
	#
	MBH_GIT="https://github.com/mtauraso/LISA-Massive-Black-Hole.git"
	MBH_BRANCH="build-fixup"
	MBH_SOURCE_DIR="${BUILD_DIR}/mbh-src"
	MBH_INSTALL_DIR="${BUILD_DIR}/mbh-install"

	git clone --branch ${MBH_BRANCH} ${MBH_GIT} ${MBH_SOURCE_DIR}
	pushd ${MBH_SOURCE_DIR}
		./install.sh ${MBH_INSTALL_DIR}
	popd
	
	# Tell cmake where to find MBH
	CMAKE_PREFIX_PATH="${MBH_INSTALL_DIR}:${CMAKE_PREFIX_PATH}"
	LDASOFT_INSTALL_DIR="${BUILD_DIR}/ldasoft-install"
	
	# Build and install ldasoft
	pushd ${SOURCE_DIR}
		./install.sh ${LDASOFT_INSTALL_DIR}
	popd
popd
