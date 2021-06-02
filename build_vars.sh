#!/bin/sh

# Directories
SWMM5PLUS_DIR=$PWD
PARALLEL_DIR='parallel'
INIT_DIR='initialization'
API_DIR='initialization/interface/flibswmm'
API_TEST_DIR='initialization/interface/tests'
JSON_DIR='json-fortran'
SIM_DIR='simulation'
TEST_DIR='test_cases'
UTIL_DIR='utilities'
VARS_DIR='vars'
ARRAY_DIR='arrays'

# Compiler vars
FC='gfortran'
OUPTFLAGS=-g
FFLAGS=-O3
PROGRAM=SWMM

# Dependencies source code
MPICH_SOURCE="$SWMM5PLUS_DIR/mpich"
CMAKE_SOURCE="$SWMM5PLUS_DIR/cmake"
COARRAY_SOURCE="$SWMM5PLUS_DIR/opencoarray"

# Dependencies install
MPICH_INSTALL="$MPICH_SOURCE/mpich-install"
CMAKE_INSTALL="$CMAKE_SOURCE/cmake-install"
COARRAY_INSTALL="$COARRAY_SOURCE/opencoarray-install"
COARRAY_FC="${COARRAY_INSTALL}/bin/caf"

# Debugging
DEBUG_API=false
DEBUG_SOURCES=""

CAF_VERSION="0.0.0"
CMAKE_VERSION="0.0.0" #initialize by using 0.0.0, not version can lower than 0.0.0
MPICH_VERSION="0.0.0"

#GCC_REQUIRE_VERSION="10.1.0"
CMAKE_REQUIRE_VERSION="3.10.0"
MPICH_REQUIRE_VERSION="3.2.0"

INSTALLATION_LOG="$SWMM5PLUS_DIR/installation_log.txt"

package_executable_array=(
    "cmake:cmake"
    "mpich:mpifort"
)

# Check OS
unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=linux;;
    Darwin*)    machine=mac;;
    CYGWIN*)    machine=cygwin;;
    MINGW*)     machine=mingw;;
    *)          machine="UNKNOWN:${unameOut}"
esac

# Version Comparison
vercomp () {
    if [[ $1 == $2 ]]
    then
        return 0
    fi
    local IFS=.
    local i ver1=($1) ver2=($2)
    # fill empty fields in ver1 with zeros
    for ((i=${#ver1[@]}; i<${#ver2[@]}; i++))
    do
        ver1[i]=0
    done
    for ((i=0; i<${#ver1[@]}; i++))
    do
        if [[ -z ${ver2[i]} ]]
        then
            # fill empty fields in ver2 with zeros
            ver2[i]=0
        fi
        if ((10#${ver1[i]} > 10#${ver2[i]}))
        then
            return 1
        fi
        if ((10#${ver1[i]} < 10#${ver2[i]}))
        then
            return 2
        fi
    done
    return 0
}