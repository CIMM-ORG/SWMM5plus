#!/bin/bash
# --------------------------------------------------------------------------------------
# Shell script for installing dependencies
# --------------------------------------------------------------------------------------

shopt -s extglob

# Download json-fortran
# --------------------------------------------------------------------------------------
if [ ! -d 'json-fortran' ]
then
    echo
    echo "Downloading json-fortran"
    echo
    git clone 'https://github.com/jacobwilliams/json-fortran.git'
    cd json-fortran
    rm -r -v !('src'|'LICENSE')
    sudo rm -r .*
    mv src/*.* .
    rm -r src
    cd ..
fi
# --------------------------------------------------------------------------------------


# Download SWMM C
# --------------------------------------------------------------------------------------

if [ ! -d "$API_DIR/src" ]
then
    echo
    echo "Downloading SWMM 5.1.13"
    echo
    wget "https://github.com/USEPA/Stormwater-Management-Model/archive/v5.1.13.tar.gz"
    tar -xvf *.tar.gz
    rm *.tar.gz
    mv Stormwater*/src "$API_DIR/src"
    rm -r Stormwater*
fi
# --------------------------------------------------------------------------------------


# Download MPICH (required for OpenCAF)
# --------------------------------------------------------------------------------------
if [ ! -d "$MPICH_SOURCE" ]
then
    echo
    echo "Installing MPICH locally"
    echo
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
# --------------------------------------------------------------------------------------


# Download CMAKE (required for OpenCAF)
# --------------------------------------------------------------------------------------
if [ ! -d $CMAKE_SOURCE ]
then
    echo
    echo "Installing CMAKE locally"
    echo
    sleep 3.0
    mkdir $CMAKE_SOURCE
    cd $CMAKE_SOURCE
    mkdir $CMAKE_INSTALL
    wget "https://cmake.org/files/v3.11/cmake-3.11.0.tar.gz"  # use 3.3.2 for now. Versions: https://cmake.org/files/
    tar -xvf *.tar.gz
    rm *.tar.gz
    cd cmake-3.11.0
    ./configure --prefix=$CMAKE_INSTALL
    make
    make install
    cd ../../
fi
# --------------------------------------------------------------------------------------


# Download OpenCAF
# --------------------------------------------------------------------------------------
if [ ! -d $COARRAY_SOURCE ]
then
    echo
    echo "Installing OpenCAF locally"
    echo
    mkdir $COARRAY_SOURCE
    cd $COARRAY_SOURCE
    mkdir $COARRAY_INSTALL
    echo Installing Opencoarray from https://github.com/sourceryinstitute/OpenCoarrays
    git clone https://github.com/sourceryinstitute/OpenCoarrays
    cd OpenCoarrays
    mkdir opencoarrays-build
    cd opencoarrays-build
    CC=gcc FC=gfortran ${CMAKE_INSTALL}/bin/cmake .. -DCMAKE_INSTALL_PREFIX=$COARRAY_INSTALL -DMPI_HOME=$MPICH_INSTALL
    make
    sudo make install
    cd ../../../
fi
# --------------------------------------------------------------------------------------

echo
echo Completed - Installation of dependencies
echo
