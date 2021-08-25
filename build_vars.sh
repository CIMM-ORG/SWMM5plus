#!/bin/sh

# Check OS
unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=linux;;
    Darwin*)    machine=mac;;
    CYGWIN*)    machine=cygwin;;
    MINGW*)     machine=mingw;;
    *)          machine="UNKNOWN:${unameOut}"
esac

# Directories
SWMM5PLUS_DIR=$PWD
INSTALLATION_LOG="$SWMM5PLUS_DIR/installation_log.txt"

# Installation log
if [[ -f $INSTALLATION_LOG ]]
then
    rm $INSTALLATION_LOG
fi

INIT_DIR='initialization'
API_DIR='interface'
UTIL_DIR='utility'
DEF_DIR='definitions'
DDIR="$SWMM5PLUS_DIR/dependencies"
GEO_DIR='geometry'
TL_DIR='timeloop'
FIN_DIR='finalization'
TEST_DIR='tests'
OUT_DIR='output'

if [[ -f $TEST_DIR/main.f08 ]]
then
    MAIN_DIR=$TEST_DIR
else
    MAIN_DIR='main'
fi

if [[ ! -d $TEST_DIR ]]
then
    mkdir $TEST_DIR
fi

# Compiler vars
FC='gfortran'
OUPTFLAGS=-g
FFLAGS=-O3
PROGRAM=SWMM

# Dependencies source code
JSON_DIR="$DDIR/json-fortran"
MPICH_SOURCE="$DDIR/mpich"
CMAKE_SOURCE="$DDIR/cmake"
COARRAY_SOURCE="$DDIR/opencoarray"

# Dependencies install
MPICH_INSTALL="$MPICH_SOURCE/mpich-install"
CMAKE_INSTALL="$CMAKE_SOURCE/cmake-install"
COARRAY_INSTALL="$COARRAY_SOURCE/opencoarray-install"

if [[ $machine = "linux" ]]
then
    CAF="$COARRAY_INSTALL/bin/caf"
    echo "export CAF_RUN=$COARRAY_INSTALL/bin/cafrun" >> $INSTALLATION_LOG
elif [[ $machine = "mac" ]]
then
    CAF= "$COARRAY_INSTALL/bin/caf" #"caf"
    echo "export CAF_RUN=$COARRAY_INSTALL/bin/cafrun" >> $INSTALLATION_LOG
fi

CAF_VERSION="0.0.0"
CMAKE_VERSION="0.0.0" #initialize by using 0.0.0, not version can lower than 0.0.0
MPICH_VERSION="0.0.0"

#GCC_REQUIRE_VERSION="10.1.0"
CMAKE_REQUIRE_VERSION="3.10.0"
MPICH_REQUIRE_VERSION="3.2.0"

package_executable_array=(
    "cmake:cmake"
    "mpich:mpifort"
)

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
