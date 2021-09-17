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
    rm -r -v !('src'|'LICENSE')
    sudo rm -r .*
    mv src/*.* .
    rm -r src
    cd $SWMM5PLUS_DIR
fi

# --------------------------------------------------------------------------------------
# Download EPA-SWMM
# --------------------------------------------------------------------------------------

if [[ ! $skip_swmm = "true" ]]
then
    if [[ ! $skip_fortran = "true" ]]
    then
        if [[ -d interface/src ]]
        then
            rm -r interface/src
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
        wget "https://github.com/USEPA/Stormwater-Management-Model/archive/v5.1.13.tar.gz"
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
    mv Stormwater*/src "$API_DIR/src"
    rm -r Stormwater*
fi

# --------------------------------------------------------------------------------------
# Install MPICH (required for OpenCAF)
# --------------------------------------------------------------------------------------

install_mpich()
{
    if [ ! -e $MPICH_INSTALL/bin/mpifort ] # local mpich not found
    then
        if [ -d $MPICH_SOURCE ]
        then
            rm -r -f $MPICH_SOURCE # clean and remake
        fi
        echo "Local mpich installation does not exist, creating folder: $MPICH_SOURCE "
        echo "Installing the prerequisite (mpich) for opencoarray fortran ..."
        mkdir -p $MPICH_SOURCE
        cd $MPICH_SOURCE
        mkdir -p $MPICH_INSTALL
        if [[ $machine = "linux" ]]
        then
            wget "https://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz"
        elif [[ $machine = "mac" ]]
        then
            curl -L "https://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz" > mpich-3.2.tar.gz
        else
            echo
            echo "OS is not supported (only mac, linux)"
            echo
            exit
        fi
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
        if [ ! -e $MPICH_INSTALL/bin/mpifort ] # if the installation fail
        then
            echo "[ERROR] MPICH installation failed."
            exit
        else
            echo "MPICH installed in $MPICH_INSTALL"
            export MPICH_PATH=$MPICH_INSTALL
            echo "mpich path: ${MPICH_PATH}" >> $INSTALLATION_LOG
        fi
    else
        echo "MPICH found in local directory: $MPICH_INSTALL/bin/mpifort"
        export MPICH_PATH=$MPICH_INSTALL
    fi

}

# --------------------------------------------------------------------------------------
# Install CMAKE (required for OpenCAF)
# --------------------------------------------------------------------------------------

install_cmake()
{
    if [ ! -e $CMAKE_INSTALL/bin/cmake ] # if local cmake not found
    then
        cd $SWMM5PLUS_DIR # make sure we are at the right level
        if [ -d $CMAKE_SOURCE ]
        then
            rm -r -f $CMAKE_SOURCE # clean and remake
        fi
        mkdir -p $CMAKE_SOURCE
        cd $CMAKE_SOURCE
        mkdir -p $CMAKE_INSTALL
        if [[ $machine = "linux" ]]
        then
            wget "https://cmake.org/files/v3.11/cmake-3.11.0.tar.gz"  # Versions: https://cmake.org/files/
        elif [[ $machine = "mac" ]]
        then
            curl -L "https://cmake.org/files/v3.11/cmake-3.11.0.tar.gz" > cmake-3.11.0.tar.gz
        else
            echo
            echo "OS is not supported (only mac, linux)"
            echo
            exit
        fi
        tar -xvf *.tar.gz
        rm *.tar.gz
        cd cmake-3.11.0
        ./configure --prefix=$CMAKE_INSTALL
        make
        make install
        if [ ! -e $CMAKE_INSTALL/bin/cmake ] # local cmake not found -> install fail
        then
            echo "[ERROR] CMAKE installation failed"
            exit
        else
            echo "CMAKE installed in $CMAKE_INSTALL"
            export CMAKE_EXEC=$CMAKE_INSTALL/bin/cmake
            echo "cmake path: $CMAKE_EXEC" >> $INSTALLATION_LOG
        fi
        cd $SWMM5PLUS_DIR # back to the SWMM directory
    else
        echo "CMAKE found in local directory: $CMAKE_INSTALL/bin/cmake"
        export CMAKE_EXEC=$CMAKE_INSTALL/bin/cmake
    fi
}

# --------------------------------------------------------------------------------------
# Install OpenCAF
# --------------------------------------------------------------------------------------

install_opencoarray_linux()
{
    # Download Opencoarray
    if [ ! -e $COARRAY_INSTALL/bin/caf ] || [ ! -e $COARRAY_INSTALL/bin/cafrun ]
    then
        echo "opencoarray not found in local directory ... "
        rm -r -f $COARRAY_SOURCE # clean it and remake
        mkdir -p $COARRAY_INSTALL
        cd $COARRAY_SOURCE
        echo Installing Opencoarray from https://github.com/sourceryinstitute/OpenCoarrays
        git clone https://github.com/sourceryinstitute/OpenCoarrays
        cd OpenCoarrays
        mkdir opencoarrays-build
        cd opencoarrays-build
        echo ${CMAKE_EXEC}
        CC=gcc FC=gfortran ${CMAKE_EXEC} .. -DCMAKE_INSTALL_PREFIX=$COARRAY_INSTALL -DMPI_HOME=$MPICH_PATH
        make
        echo "Installing Opencoarrays ... "
        sudo make install
        cd $SWMM5PLUS_DIR
    fi
    if [ -e $COARRAY_INSTALL/bin/caf ]
    then
        echo "opencoarray installed in $COARRAY_INSTALL/bin/caf"
        export CAFRUN=$COARRAY_INSTALL/bin/cafrun # export variable
        echo "opencoarray path: $COARRAY_INSTALL/bin/cafrun" >> $INSTALLATION_LOG
    else
        echo "[ERROR] opencoarray local installation failed."
        echo "Now trying default installation process in OpenCoarray."
        cd $COARRAY_SOURCE
        ./install.sh --install-prefix $COARRAY_INSTALL
        if [ ! -e $COARRAY_INSTALL/bin/caf ]
        then
            echo "OpenCoarray Installation Failed."
            echo "Please check https://github.com/sourceryinstitute/OpenCoarrays/blob/master/INSTALL.md#linux for library compatibility."
            exit
        fi
    fi
}

install_opencoarray_mac()
{
    if [ ! -e $COARRAY_INSTALL/bin/caf ] || [ ! -e $COARRAY_INSTALL/bin/cafrun ]
    then
        echo "opencoarray not found in local directory ... "
        rm -r -f $COARRAY_SOURCE # clean it and remake
        mkdir -p $COARRAY_INSTALL
        cd $COARRAY_SOURCE
        echo Installing Opencoarray from https://github.com/sourceryinstitute/OpenCoarrays
        git clone https://github.com/sourceryinstitute/OpenCoarrays
        cd OpenCoarrays
        ./install.sh --install-prefix $COARRAY_INSTALL
        cd $SWMM5PLUS_DIR
    fi
}

opencoarray_prerequisite()
{   # For simplicity, install everything in local directory.
    echo "Installing cmake in $CMAKE_INSTALL ..."
    install_cmake
    echo "Cmake installation complete. "
    echo "Installing mpich in $MPICH_INSTALL ..."
    install_mpich
    echo "MPICH installation complete. "
}


# --------------------------------------------------------------------------------------

if [[ ! -f $CAF ]]
then
    if [[ $machine = "linux" ]]
    then
        if [[ $tacc = true ]]
        then
            CAF="ifort -coarray=distributed"
        else
            opencoarray_prerequisite
            install_opencoarray_linux
        fi
    elif [[ $machine = "mac" ]]
    then
        install_opencoarray_mac # If user want to use Homebrew to install OpenCoarrays, please comment out this line
        CAF="$COARRAY_INSTALL/bin/caf" #If user want to use Homebrew to install OpenCoarrays, please change the CAF path
    fi
fi
# --------------------------------------------------------------------------------------

# Compile EPA-SWMM

echo
echo Compiling SWMM 5.13 ...
echo

cp -f "$API_DIR/api.h" "$API_DIR/src/"
cp -f "$API_DIR/api.c" "$API_DIR/src/"

# Insert new files in EPA-SWMM Makefile

SCRIPTS='api.o'
OBJECTS='api.o   : headers.h api.h\
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
