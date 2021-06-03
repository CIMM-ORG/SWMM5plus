#!/bin/bash

# --------------------------------------------------------------------------------------
# Shell script for installing dependencies
# --------------------------------------------------------------------------------------

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

shopt -s extglob

# Create directory for dependencies
if [[ ! -d $DDIR ]]
then
    mkdir $DDIR
fi

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

# Download SWMM C
# --------------------------------------------------------------------------------------
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
    export MPICH_PATH=$MPICH_INSTALL
    echo "mpich path: ${MPICH_PATH}" >> $INSTALLATION_LOG
}
# --------------------------------------------------------------------------------------

# Install CMAKE (required for OpenCAF)
# --------------------------------------------------------------------------------------
install_cmake()
{
    cd $SWMM5PLUS_DIR # make sure we are at the right level
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
    cd $SWMM5PLUS_DIR # back to the SWMM directory
    export CMAKE_EXEC=$CMAKE_INSTALL/bin/cmake
    echo "cmake path: $CMAKE_EXEC" >> $INSTALLATION_LOG
}
# --------------------------------------------------------------------------------------

# Install OpenCAF
# --------------------------------------------------------------------------------------
opencoarray_prerequisite()
{
    for element in "${package_executable_array[@]}" ; do
        KEY="${element%%:*}"
        VALUE="${element##*:}"

        if ! command -v $VALUE &> /dev/null; # if cannot find the command in root (usually usr/bin or usr/local/bin)
        then
            case $KEY in
                "cmake")
                    echo "${KEY} is not installed in root. Searching in ${CMAKE_SOURCE} ..."
                    if [ -d $CMAKE_SOURCE ] && [ -f $CMAKE_INSTALL/bin/cmake ] # if cmake already in local
                    then
                        echo "Found cmake in $CMAKE_INSTALL ..."

                        CMAKE_VERSION=$(echo $(${CMAKE_INSTALL}/bin/cmake --version) | head -1 | sed 's/[^0-9.]*\([0-9.]*\).*/\1/')
                        vercomp $CMAKE_REQUIRE_VERSION $CMAKE_VERSION
                        result=$?
                        if [[ $result = 0 ]] || [[ $result = 2 ]]
                        then
                            echo "Current cmake version is ${CMAKE_VERSION}, higher than required version (${CMAKE_REQUIRE_VERSION})."
                            export CMAKE_EXEC=$(which cmake)
                            echo "cmake path: $CMAKE_EXEC" >> $INSTALLATION_LOG
                        else # version is outdated
                            echo "Local cmake is outdated. Installing a newer version of cmake ..."
                            sudo rm -rf $CMAKE_SOURCE
                            install_cmake
                        fi

                    else
                        echo "No local cmake found ... Installing cmake in ${CMAKE_SOURCE}"
                        if [ -d $CMAKE_SOURCE ]; then
                            sudo rm -rf $CMAKE_SOURCE
                        fi
                        install_cmake
                    fi
                    ;;

                "mpich")
                    echo "${KEY} is not installed in root. Searching in ${MPICH_SOURCE} ..."
                    if [ -d $MPICH_SOURCE ] && [ -x $MPICH_INSTALL/bin/mpifort ]
                    then
                        echo "Found mpich in $MPICH_INSTALL ..."
                        MPICH_VERSION=$(echo $(${MPICH_INSTALL}/bin/mpichversion --version) | head -1 | sed 's/[^0-9.]*\([0-9.]*\).*/\1/')
                        vercomp $MPICH_REQUIRE_VERSION $MPICH_VERSION
                        result=$?
                        if [[ $result = 0 ]] || [[ $result = 2 ]]
                        then
                            echo "Current ${KEY} version is ${MPICH_VERSION}, higher than required version (${MPICH_REQUIRE_VERSION})."
                            export MPICH_PATH=$MPICH_INSTALL # use the mpich in root
                            echo "mpich path: $MPICH_PATH" >> $INSTALLATION_LOG
                        else
                            echo "Local ${KEY} is outdated, install the updated version (mpich 3.2.0) ..."
                            read -p "Delete the local ${KEY} directory and start the re-installation? [Y/N]:" -n 1 -r
                            echo    # (optional) move to a new line
                            if [[ $REPLY =~ ^[Yy]$ ]]
                            then
                                sudo rm -rf $MPICH_SOURCE
                            else
                                echo "Did not delete the existing ${MPICH_SOURCE}, ${KEY} reinstallation suspended."
                                echo "[WARNING] mpich is not installed [WARNING]"
                                exit 0
                            fi
                            install_mpich
                        fi
                    else
                        echo "No local ${KEY} found .... Installing ${KEY} in ${MPICH_SOURCE}"
                        if [ -d $MPICH_SOURCE ]
                        then
                            sudo rm -rf $MPICH_SOURCE
                        fi
                        install_mpich
                    fi
                    ;;
            esac

        elif command -v $VALUE &> /dev/null; # command exists
        then
            case $KEY in
                "cmake")
                    CMAKE_VERSION=$(echo $($VALUE --version) | head -1 | sed 's/[^0-9.]*\([0-9.]*\).*/\1/')
                    vercomp $CMAKE_REQUIRE_VERSION $CMAKE_VERSION
                    result=$?
                    if [[ $result = 0 ]] || [[ $result = 2 ]]
                    then
                        echo "Current ${KEY} version is ${CMAKE_VERSION}, higher than required version (${CMAKE_REQUIRE_VERSION})."
                        export CMAKE_EXEC=$(which $VALUE)
                        echo "${KEY} path: $CMAKE_EXEC" >> $INSTALLATION_LOG
                    else
                        echo "Current ${KEY} version is: ${CMAKE_VERSION}, lower than the required version: ${CMAKE_REQUIRE_VERSION}."
                        echo "Install local ${KEY} under ${CMAKE_SOURCE} ...."
                        if [ -d $CMAKE_SOURCE ] && [ -f $CMAKE_INSTALL/bin/cmake ] # if cmake already in local
                        then
                            CMAKE_VERSION=$(echo $(${CMAKE_INSTALL}/bin/cmake --version) | head -1 | sed 's/[^0-9.]*\([0-9.]*\).*/\1/')
                            vercomp $CMAKE_REQUIRE_VERSION $CMAKE_VERSION
                            result=$?
                            if [[ $result = 0 ]] || [[ $result = 2 ]]
                            then
                                export CMAKE_EXEC=$CMAKE_INSTALL/bin/cmake
                            else #the local cmake is also outdated, make install
                                echo "Local ${KEY} is outdated, install the updated version (cmake 3.11.0) ..."
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
                                ### cmake installation body
                                install_cmake
                            fi
                        else # no cmake in local -> make the install
                            if [ -d $CMAKE_SOURCE ]
                            then
                                sudo rm -rf $CMAKE_SOURCE
                            fi
                            install_cmake
                        fi
                    fi
                ;;

                "mpich")
                    MPICH_VERSION=$(echo $(mpichversion --version) | head -1 | sed 's/[^0-9.]*\([0-9.]*\).*/\1/')
                    vercomp $MPICH_REQUIRE_VERSION $MPICH_VERSION
                    result=$?
                    if [[ $result = 0 ]] || [[ $result = 2 ]]
                    then
                        echo "Current ${KEY} version is ${MPICH_VERSION}, higher than required version (${MPICH_REQUIRE_VERSION})."
                        export MPICH_PATH="/usr/" # return the mpich path in root
                        echo "${KEY} path : $MPICH_PATH" >> $INSTALLATION_LOG
                    else
                        echo "Current ${KEY} version in root is: ${MPICH_VERSION}, lower than the required version: ${MPICH_REQUIRE_VERSION}."
                        echo "Install local ${KEY} under ${MPICH_SOURCE} ...."
                        if [ -d $MPICH_SOURCE ] && [ -x $MPICH_INSTALL/bin/mpifort ] # if we have mpich in local
                        then
                            echo "Found existing ${KEY} in local directory."
                            MPICH_VERSION=$(echo $(${MPICH_INSTALL}/bin/mpichversion --version) | head -1 | sed 's/[^0-9.]*\([0-9.]*\).*/\1/')
                            vercomp $MPICH_REQUIRE_VERSION $MPICH_VERSION
                            result=$?
                            if [[ $result = 0 ]] || [[ $result = 2 ]]
                            then
                                echo "Current ${KEY} version is ${MPICH_VERSION}, higher than required version (${MPICH_REQUIRE_VERSION})."
                                export MPICH_PATH=$MPICH_INSTALL # use the mpich in root
                                echo "${KEY} path: $MPICH_PATH" >> $INSTALLATION_LOG
                            else
                                echo "Local ${KEY} is outdated, install the updated version (mpich 3.2.0) ..."
                                read -p "Delete the local ${KEY} directory and start the re-installation? [Y/N]:" -n 1 -r
                                echo    # (optional) move to a new line
                                if [[ $REPLY =~ ^[Yy]$ ]]
                                then
                                    sudo rm -rf $MPICH_SOURCE
                                else
                                    echo "Did not delete the existing ${MPICH_SOURCE}, ${KEY} reinstallation suspended."
                                    echo "[WARNING] mpich is not installed [WARNING]"
                                    exit 0
                                fi
                                install_mpich
                            fi
                        else # no mpich in root and local -> make install
                            if [ -d $MPICH_SOURCE ]
                            then
                                sudo rm -rf $MPICH_SOURCE
                            fi
                            install_mpich
                        fi
                    fi
                ;;

            esac
        fi

    done
}

install_opencoarray()
{
    # Download Opencoarray
    if [ ! -d $COARRAY_SOURCE ]
    then
        echo "opencoarray is not found in $COARRAY_SOURCE"
        mkdir -p $COARRAY_INSTALL
        cd $COARRAY_SOURCE
        echo Installing Opencoarray from https://github.com/sourceryinstitute/OpenCoarrays
        git clone https://github.com/sourceryinstitute/OpenCoarrays
        cd OpenCoarrays
        mkdir opencoarrays-build
        cd opencoarrays-build
        CC=gcc FC=gfortran ${CMAKE_EXEC} .. -DCMAKE_INSTALL_PREFIX=$COARRAY_INSTALL -DMPI_HOME=$MPICH_PATH
            make
        echo "Installing Opencoarrays ... "
        sudo make install
        cd $SWMM5PLUS_DIR
    fi
    echo "opencoarray path: $COARRAY_INSTALL/bin/cafrun" >> $INSTALLATION_LOG
}

# --------------------------------------------------------------------------------------

if [[ ! -f $CAF ]]
then
    opencoarray_prerequisite
    install_opencoarray
fi
# --------------------------------------------------------------------------------------

# Compile SWMM C

echo
echo Compiling SWMM 5.13 ...
echo

cp -f "$API_DIR/interface.h" "$API_DIR/src/"
cp -f "$API_DIR/interface.c" "$API_DIR/src/"

# Insert new files in SWMM C Makefile

SCRIPTS="interface.o"
OBJECTS="interface.o   : headers.h interface.h\n"

API_TEST_FILES=""
for fname in $(find $TEST_DIR -name '*.c')
do
    F=$(basename -- "$fname")
    F="${F%.*}"
    SCRIPTS="$SCRIPTS $F.o"
    if [[ -f "${fname%.*}.h" ]]
    then
        OBJECTS="${OBJECTS}$F.o       : headers.h $F.h\n"
        cp -f "$TEST_DIR/$F.h" "$API_DIR/src/"
    else
        OBJECTS="${OBJECTS}$F.o       : headers.h\n"
    fi
    cp -f "$TEST_DIR/$F.c" "$API_DIR/src/"
done

sed "s#{{SCRIPTS}}#$SCRIPTS#" "$API_DIR/Makefile" > "$API_DIR/src/Makefile"
sed -i "s#{{OBJECTS}}#$OBJECTS#" "$API_DIR/src/Makefile"

cd "$API_DIR/src" && make
cd $SWMM5PLUS_DIR
cp $API_DIR/src/libswmm5.so $SWMM5PLUS_DIR/libswmm5.so

# --------------------------------------------------------------------------------------

echo
echo Completed Installation of Dependencies!
echo
