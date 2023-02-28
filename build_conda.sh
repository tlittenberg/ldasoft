#!/bin/bash
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace

# Build LDASOFT and MBH from scratch w/ Conda

# Assume we run from ldasoft root (TODO fix this so we don't need to assume)
# We use build-conda as a staging area.
BUILD_DIR="$(pwd)/build-conda"

SOURCE_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

rm -rf ${BUILD_DIR}
mkdir -p ${BUILD_DIR}
pushd ${BUILD_DIR}

	# Build a brand-new conda env in BUILD_DIR
	conda env create -f $SOURCE_DIR/environment.yml -p ${BUILD_DIR}/conda_env

	# Activate our conda env
	eval $(conda shell.bash activate ${BUILD_DIR}/conda_env)

	# Tell cmake where to find conda stuff
	export CMAKE_PREFIX_PATH="${CONDA_PREFIX}:${CMAKE_PREFIX_PATH:-}"

	# Module load our stuff from hyak
	module load ompi/4.1.4 # Note this also loads gcc and ucx
	module load cmake

	# Hyaks load of gcc/11.2.x sets $CC to mpiCC which is technically a c++ compiler
	# cmake gets mad about this so we set CC in order to point it right
	export CC="mpicc"

	# Build and install MBH from source
	#
	# TODO Ergonomics: What is the right place to get MBH from? is git okay, or
	# will this typically be a local checkout that we should build the ability
	# to specify

	export MBH_GIT="https://github.com/mtauraso/LISA-Massive-Black-Hole.git"
	export MBH_BRANCH="build-fixup"

	git clone --branch ${MBH_BRANCH} ${MBH_GIT} ${BUILD_DIR}/mbh

	pushd ${BUILD_DIR}/mbh
		./install.sh ${BUILD_DIR}/mbh-install
	popd
	
	# Tell cmake where mbh is
	export CMAKE_PREFIX_PATH="${BUILD_DIR}/mbh-install/:${CMAKE_PREFIX_PATH}"
	
	# Build and install ldasoft
	pushd ${SOURCE_DIR}
		./install.sh ${BUILD_DIR}/ldasoft-install
	popd
popd



