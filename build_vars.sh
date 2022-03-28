#!/bin/sh

num_processors=1
compile_fortran="true"
compile_swmmc="true"
download_swmmc="true" # brh20211209 allows skip of download (dangerous!)

# Check flags
while [ ! $# -eq 0 ]
do
	case "$1" in
		--swmmc | -s) # Compile SWMM C only
			compile_fortran="false" ;;
		--fortran | -f) # Compile SWMM5+ only
			compile_swmmc="false" ;; 
        --processors | -n)
            num_processors="$2" ;;
        --clean | -c)
            clean="true" ;;
        --existing_swmm | -z) # compile existing SWMM C (no download)  brh2021209 
            download_swmmc="false" ;; 
	esac
	shift
done

# Check OS
unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=linux;;
    Darwin*)    machine=mac;;
    CYGWIN*)    machine=cygwin;;
    MINGW*)     machine=mingw;;
    *)          machine="UNKNOWN:${unameOut}"
esac


# --------------------------------------------------------------------------------------
# Set number of processors

if grep "export FOR_COARRAY_NUM_IMAGES" ~/.bashrc
then
    sed -i "s/export FOR_COARRAY_NUM_IMAGES=[0-9]*/export FOR_COARRAY_NUM_IMAGES=$num_processors/g" ~/.bashrc
else
    echo "export FOR_COARRAY_NUM_IMAGES=$num_processors" >> ~/.bashrc
fi
source ~/.bashrc
# --------------------------------------------------------------------------------------


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
PROGRAM=SWMM

# Dependencies source code
JSON_DIR="$DDIR/json-fortran"
