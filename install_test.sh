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

CAF_VERSION="0.0.0"
CMAKE_VERSION="0.0.0" #initialize by using 0.0.0, not version can lower than 0.0.0
MPICH_VERSION="0.0.0"

#GCC_REQUIRE_VERSION="10.1.0"
CMAKE_REQUIRE_VERSION="3.10.0"
MPICH_REQUIRE_VERSION="3.2.0"

INSTALLATION_LOG="$SWMM5PLUS_DIR/dependencies_installation_log.txt"

# Default version for packages
# cmake (default version 3.10.0)        https://www.cmake.org/files/v3.10/cmake-3.10.0-Linux-x86_64.sh
# gcc (default version 10.1.0)          https://ftp.gnu.org/gnu/gcc/gcc-10.1.0/gcc-10.1.0.tar.gz
# mpich (default version 3.2)           https://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz
# wget (default version 1.16.3)         https://ftpmirror.gnu.org/gnu/wget/wget-1.16.3.tar.gz
# flex (default version 2.6.0)          https://sourceforge.net/projects/flex/files/flex-2.6.0.tar.bz2
# bison (default version 3.0.4)         https://ftpmirror.gnu.org/gnu/bison/bison-3.0.4.tar.gz
# pkg-config (default version 0.28)     https://pkgconfig.freedesktop.org/releases/pkg-config-0.28.tar.gz
# make (default version 4.1)            https://ftpmirror.gnu.org/gnu/make/make-4.1.tar.bz2
# m4 (default version 1.4.17)           https://ftpmirror.gnu.org/gnu/m4/m4-1.4.17.tar.bz2
# subversion (default version 1.9.4)    https://www.eu.apache.org/dist/subversion/subversion-1.9.4.tar.gz
# ofp (default version sdf) - this is only needed for OSX
# opencoarrays (version 2.9.2)          

if [ -f $INSTALLATION_LOG ]
then
    rm $INSTALLATION_LOG # make sure we create a new one everytime
fi


package_executable_array=(
    "cmake:cmake"
    "mpich:mpifort"
    "opencoarray:caf"
  )

  package_default_version=(
    "cmake:3.4.0"
    "mpich:3.2"
    "opencoarray:2.9.2"
  )

if [[ -x "/usr/bin/caf" && -x "usr/bin/cafrun" ]]
then
    echo "Opencoarray is already installed in /usr/bin ..."
    # set the compiler path =caf/cafrun
    export COARRAY_FC='caf'
    echo "coarray path: $COARRAY_FC" >> $INSTALLATION_LOG
    exit 0
else
    

for element in "${package_executable_array[@]}" ; do
    KEY="${element%%:*}"
    VALUE="${element##*:}"
    
    if ! command -v $VALUE &> /dev/null; # if cannot find the command in root (usually usr/bin or usr/local/bin)
    then
        case $KEY in 
            "cmake")
                echo "${KEY} is not installed in root. Will install in under ${CMAKE_SOURCE}"

                ;;
            "mpich")
                echo "${KEY} is not installed in root. Will install in under ${MPICH_SOURCE}"
                ;;

        esac

    

    elif command -v $VALUE &> /dev/null; # command exists
    then
        case $KEY in
            "cmake")
                CMAKE_VERSION=$(echo $(cmake --version) | head -1 | sed 's/[^0-9.]*\([0-9.]*\).*/\1/')
                if [ "$(printf '%s\n' "$CMAKE_REQUIRE_VERSION" "$CMAKE_VERSION" | sort -V | head -n1)" = "$CMAKE_REQUIRE_VERSION" ]
                then
                    echo "Current cmake version is ${CMAKE_VERSION}, higher than required version (${CMAKE_REQUIRE_VERSION})."
                    export CMAKE_EXEC=$(which cmake)
                    echo "cmake path: $CMAKE_EXEC" >> $INSTALLATION_LOG
                else
                    echo "Current cmake version is: ${CMAKE_VERSION}, lower than the required version: ${CMAKE_REQUIRE_VERSION}."
                    echo "Install local cmake under ${CMAKE_SOURCE} ...."
                    if [ -d $CMAKE_SOURCE ] && [ -f $CMAKE_INSTALL/bin/cmake ] # if cmake already in local
                    then 
                        CMAKE_VERSION=$(echo $(${CMAKE_INSTALL}/bin/cmake --version) | head -1 | sed 's/[^0-9.]*\([0-9.]*\).*/\1/')
                        if [ "$(printf '%s\n' "$CMAKE_REQUIRE_VERSION" "$CMAKE_VERSION" | sort -V | head -n1)" = "$CMAKE_REQUIRE_VERSION" ]
                        then
                            export CMAKE_EXEC=$CMAKE_INSTALL/bin/cmake
                        else #the local cmake is also outdated, make install
                            echo "Local cmake is outdated, install the updated version (cmake 3.11.0) ..."
                            read -p "Delete the local ${KEY} directory and start the re-installation? [Y/N]:" -n 1 -r
                            echo    # (optional) move to a new line
                            if [[ $REPLY =~ ^[Yy]$ ]]
                            then
                                sudo rm -rf $CMAKE_SOURCE
                            else
                                echo "Did not delete the existing ${CMAKE_SOURCE}, ${KEY} reinstallation suspended."
                                echo "[WARNING] cmake is not installed [WARNING]"
                                exit 0
                            fi
                            ### cmake installation body ### [Should make this become a function later]
                            cd $SWMM5PLUS_DIR # make sure we are at the right level
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
                            cd $SWMM5PLUS_DIR # back to the SWMM directory
                            export CMAKE_EXEC=$CMAKE_INSTALL/bin/cmake
                            echo "cmake path: $CMAKE_EXEC" >> $INSTALLATION_LOG
                        fi
                    else # no cmake in local -> make the install
                        if [ -d $CMAKE_SOURCE ]; then
                            sudo rm -rf $CMAKE_SOURCE
                        fi
                        cd $SWMM5PLUS_DIR # make sure we are at the right level
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
                        cd $SWMM5PLUS_DIR # back to the SWMM directory
                    fi
                fi

            ;;
            "mpich")
                MPICH_VERSION=$(echo $(mpiexec --version) | head -1 | sed 's/[^0-9.]*\([0-9.]*\).*/\1/')
                if [ "$(printf '%s\n' "$MPICH_REQUIRE_VERSION" "$MPICH_VERSION" | sort -V | head -n1)" = "$MPICH_REQUIRE_VERSION" ]
                then
                    echo "Current cmake version is ${MPICH_VERSION}, higher than required version (${MPICH_REQUIRE_VERSION})."
                    export MPICH_EXEC=$(which mpiexec)
                else
                    echo "Current mpich version is: ${MPICH_VERSION}, lower than the required version: ${MPI_REQUIRE_VERSION}."
                    echo "Install local mpich under ${MPICH_SOURCE} ...."
                fi

            ;;

            "opencoarray")

            ;;
        esac
    fi
    
done


# Version Check - cmake
if ! command -v cmake &> /dev/null
then
    CMAKE_VERSION=$(echo $(cmake --version) | head -1 | sed 's/[^0-9.]*\([0-9.]*\).*/\1/')
fi

# Version Check - MPICH
if command -v mpiexec &> /dev/null 
then 
    MPICH_VERSION=$(echo $(mpiexec --version) | head -1 | sed 's/[^0-9.]*\([0-9.]*\).*/\1/')
fi


if [ "$(printf '%s\n' "$MPICH_REQUIRE_VERSION" "$MPICH_VERSION" | sort -V | head -n1)" = "$MPICH_REQUIRE_VERSION" ]; then
        echo "Greater than or equal to ${MPICH_REQUIRE_VERSION}"
 else
        echo "Less than ${MPICH_REQUIRE_VERSION}"
 fi


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


# Download Opencoarray
if [ ! -d $COARRAY_SOURCE ]   
then
    echo "opencoarray is not found in current directory."
    mkdir $COARRAY_SOURCE
    cd $COARRAY_SOURCE
    mkdir $COARRAY_INSTALL
    echo Installing Opencoarray from https://github.com/sourceryinstitute/OpenCoarrays
    sleep 3.0
    git clone https://github.com/sourceryinstitute/OpenCoarrays 
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
