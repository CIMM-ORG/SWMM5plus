#!/bin/bash

# --------------------------------------------------------------------------------------
#  Alpha Release:
#  Shell script for installing dependencies
#
#  [README]
#  This shell script is for downloading the dependencies of SWMM 5+. The dependencies
#  include:
#  (1) json-fortran, from: https://github.com/jacobwilliams/json-fortran.git
#  (2) SWMM 5 libraries, from: https://github.com/USEPA/Stormwater-Management-Model/archive/v5.1.13.tar.gz
#  (3) OpenCoarray and its dependencies, from: https://github.com/sourceryinstitute/OpenCoarrays
#
#  OpenCoarray is a Fortran Coarray libaray designed for GNC compiler (e.g. gfortran series).
#  Before we start the installation/compile of SWMM5+, we strongly recommand users to check the
#  of compiler version. If user is using Intel Compiler, then the OpenCoarray dependencies
#  installation can be skipped by simply adding a "-t" flag after the command "./build.sh". If
#  user is using gcc compiler, please make sure the compiler version is higher than 6.1.0.
#  The dependencies of OpenCoarray include:
#
#  opencoarrays
#  ├── cmake-3.4.0
#  └── mpich-3.2
#      └── gcc-6.1.0
#          ├── flex-2.6.0
#          │   └── bison-3.0.4
#          │       └── m4-1.4.17
#          ├── gmp
#          ├── mpc
#          └── mpfr
#
#  This script will first check if the above-mentioned pakcages exist in the local directory,
#  If not found -> then download and install.
#  For the OpenCoarray dependencies, script will check the existencce of local copy of cmake and mpich,
#  For the user first time download/compile SWMM5+, all these packages will be downloaded and installed.
#  This is an one-time installation work if user keeps working in the same SWMM5+ directory. However,
#  the dependencies installation process is required everytime after cloning a new version of SWMM 5+.
#
#  More information for the packages can refer to :
#  json fortran : https://github.com/jacobwilliams/json-fortran#brief-description
#  SWMM 5 library: https://www.epa.gov/water-research/storm-water-management-model-swmm
#  OpenCoarray: http://www.opencoarrays.org/
# --------------------------------------------------------------------------------------

shopt -s extglob

# Create directory for dependencies
if [[ ! -d $DDIR ]]
then
    mkdir $DDIR
fi

# --------------------------------------------------------------------------------------
# Download json-fortran
# --------------------------------------------------------------------------------------

if [ ! -d "$JSON_DIR" ]
then
    echo
    echo "Downloading json-fortran"
    echo
    cd $DDIR
    git clone 'https://github.com/jacobwilliams/json-fortran.git'
    cd json-fortran
    rm -rf -v !('src'|'LICENSE')
    rm -rf .*
    mv src/*.* .
    rm -rf src
    cd $SWMM5PLUS_DIR
fi

# --------------------------------------------------------------------------------------
# Download EPA-SWMM
# --------------------------------------------------------------------------------------

if [[ $compile_swmmc = "true" ]]
then
    if [[ $compile_fortran = "true" ]]
    then
        if [[ -d interface/src ]]
        then
            rm -rf interface/src
        fi
    fi
fi

if [ ! -d "$API_DIR/src" ]
then
    echo
    echo "Downloading SWMM 5.1.13"
    echo
    if [[ $machine = "linux" ]]
    then
        wget https://github.com/USEPA/Stormwater-Management-Model/archive/v5.1.13.tar.gz
    elif [[ $machine = "mac" ]]
    then
        curl -L "https://github.com/USEPA/Stormwater-Management-Model/archive/v5.1.13.tar.gz" > v5.1.13.tar.gz
    else
        echo
        echo "OS is not supported (only mac, linux)"
        echo
        exit
    fi
    tar -xvf *.tar.gz
    rm *.tar.gz
    cp -r Stormwater*/src "$API_DIR/src"
    rm -rf Stormwater*
fi
 
# Compile EPA-SWMM

echo
echo Compiling SWMM 5.13 ...
echo

cp -f "$API_DIR/api.h" "$API_DIR/src/"
cp -f "$API_DIR/api.c" "$API_DIR/src/"
cp -f "$API_DIR/api_error.h" "$API_DIR/src/"
cp -f "$API_DIR/api_error.c" "$API_DIR/src/"

# Insert new files in EPA-SWMM Makefile

SCRIPTS='api_error.o api.o'
OBJECTS='api_error.o   : headers.h api_error.h\
api.o   : api.h api_error.h headers.h\
'

API_TEST_FILES=""
for fname in $(find $TEST_DIR -name '*.c')
do
    F=$(basename -- "$fname")
    F="${F%.*}"
    SCRIPTS="$SCRIPTS $F.o"
    if [[ -f "${fname%.*}.h" ]]
    then
        OBJECTS="${OBJECTS}$F.o       : headers.h $F.h\
        "
        cp -f "$TEST_DIR/$F.h" "$API_DIR/src/"
    else
        OBJECTS="${OBJECTS}$F.o       : headers.h\
        "
    fi
    cp -f "$TEST_DIR/$F.c" "$API_DIR/src/"
done

sed "s#{{SCRIPTS}}#$SCRIPTS#" "$API_DIR/Makefile" > "$API_DIR/src/Makefile"
sed -i'.original' -e "s#{{OBJECTS}}#$OBJECTS#" "$API_DIR/src/Makefile"
rm "$API_DIR/src/Makefile.original"

cd "$API_DIR/src" && make
cd $SWMM5PLUS_DIR
cp $API_DIR/src/libswmm5.so $SWMM5PLUS_DIR/libswmm5.so

# --------------------------------------------------------------------------------------

echo
echo Completed Installation of Dependencies!
echo
