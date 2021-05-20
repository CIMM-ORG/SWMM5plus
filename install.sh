#!/bin/bash
shopt -s extglob

# shell script for installing dependencies

# Constants
SWMM5PLUS_DIR=$PWD

# Dependencies source code
MPICH_SOURCE="$SWMM5PLUS_DIR/mpich"
CMAKE_SOURCE="$SWMM5PLUS_DIR/cmake"
COARRAY_SOURCE="$SWMM5PLUS_DIR/opencoarray"

# Dependencies install
MPICH_INSTALL="$MPICH_SOURCE/mpich-install"
CMAKE_INSTALL="$CMAKE_SOURCE/cmake-install"
COARRAY_INSTALL="$COARRAY_SOURCE/opencoarray-install"




# Download dependencies

if [ ! -d 'json-fortran' ]
then
    git clone 'https://github.com/jacobwilliams/json-fortran.git'
    cd json-fortran
    rm -r -v !('src'|'LICENSE')
    sudo rm -r .*
    mv src/*.* .
    rm -r src
    cd ..
fi

# Download SWMM C

if [ ! -d "$API_DIR/src" ]
then
    echo Downloading SWMM from US EPA repository ...
    wget "https://github.com/USEPA/Stormwater-Management-Model/archive/v5.1.13.tar.gz"
    tar -xvf *.tar.gz
    rm *.tar.gz
    mv Stormwater*/src "$API_DIR/src"
    rm -r Stormwater*
fi

if [ ! -d "$MPICH_SOURCE" ]  
then
    echo "Installing the prerequisite (mpich) for opencoarray fortran ..."
    sleep 3.0
    mkdir $MPICH_SOURCE
    cd $MPICH_SOURCE
    mkdir $MPICH_INSTALL
    wget "https://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz"
    tar -xvf *.tar.gz
    rm *.tar.gz
    if [ -d "/tmp/mpich-build" ]
    then
        rm -r -f /tmp/mpich-build
    fi
    mkdir /tmp/mpich-build
    cd /tmp/mpich-build
    $MPICH_SOURCE/mpich-3.2/configure -prefix=$MPICH_INSTALL
    make
    make install
    cd $SWMM5PLUS_DIR
fi

if [ ! -d $CMAKE_SOURCE ] 
then 
    echo "cmake is not found in current directory."
    echo "Installing cmake - the prerequisite for opencoarray fortran ..."
    sleep 3.0
    mkdir $CMAKE_SOURCE
    cd $CMAKE_SOURCE
    mkdir $CMAKE_INSTALL
    wget "https://cmake.org/files/v3.3/cmake-3.3.2.tar.gz"  # use 3.3.2 for now. Versions: https://cmake.org/files/
    tar -xvf *.tar.gz
    rm *.tar.gz
    cd cmake-3.3.2
    ./configure --prefix=$CMAKE_INSTALL
    make
    make install
    cd ../../
fi


# Download Opencoarray
if [ ! -d $COARRAY_SOURCE ]   
then
    echo "opencoarray is not found in current directory."
    sleep 3.0
    mkdir $COARRAY_SOURCE
    cd $COARRAY_SOURCE
    mkdir $COARRAY_INSTALL
    echo Installing Opencoarray from https://github.com/sourceryinstitute/OpenCoarrays
    sleep 3.0
    git clone --branch 1.9.3 https://github.com/sourceryinstitute/OpenCoarrays
    cd OpenCoarrays
    mkdir opencoarrays-build
    cd opencoarrays-build
    CC=gcc FC=gfortran ${CMAKE_INSTALL}/bin/cmake .. -DCMAKE_INSTALL_PREFIX=$COARRAY_INSTALL -DMPI_HOME=$MPICH_INSTALL
    make
    sudo make install
    cd ../../../
fi

echo
echo Dependencies Installation Complete!
echo
