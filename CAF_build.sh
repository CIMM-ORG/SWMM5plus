#!/bin/bash
shopt -s extglob

# Constants

ARRAY_DIR='arrays'
INIT_DIR='initialization'
API_DIR='initialization/interface'
JSON_DIR='json-fortran'
SIM_DIR='simulation'
TEST_DIR='test_cases'
UTIL_DIR='utilities'
VARS_DIR='vars'

FC='gfortran'
OUPTFLAGS=-g
FFLAGS=-O3
PROGRAM=CAF_SWMM
DEBUG=true

if [ $DEBUG = true ]
then
    rm -r debug
    mkdir debug
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

# Compile SWMM C

echo Compiling SWMM 5.13 ...
cp "$API_DIR/interface.h" "$API_DIR/src/"
cp "$API_DIR/interface.c" "$API_DIR/src/"
cp "$API_DIR/tests.h" "$API_DIR/src/"
cp "$API_DIR/tests.c" "$API_DIR/src/"
cp "$API_DIR/Makefile" "$API_DIR/src/"
cd "$API_DIR/src"
make
cd ..
cp src/libswmm5.so ../../libswmm5.so
cd ../..

# Compile

SOURCESF="$JSON_DIR/json_kinds.F90\
          $JSON_DIR/json_parameters.F90\
          $JSON_DIR/json_string_utilities.F90\
          $JSON_DIR/json_value_module.F90\
          $JSON_DIR/json_file_module.F90\
          $JSON_DIR/json_module.F90\
          $VARS_DIR/data_keys.f08\
          $VARS_DIR/type_definitions.f08\
          $INIT_DIR/setting_definition.f08\
          $VARS_DIR/globals.f08\
          $API_DIR/dll.f08\
          $UTIL_DIR/datetime.f08\
          $ARRAY_DIR/dynamic_array.f08\
          $ARRAY_DIR/tables.f08\
          $API_DIR/interface.f08\
          $VARS_DIR/assign_index.f08\
          $UTIL_DIR/utility.f08\
          $INIT_DIR/allocate_storage.f08\
          $INIT_DIR/initialization.f08\
          $UTIL_DIR/BIPquick.f08
          main_caf.f08"

echo Compiling ...
$FC $SOURCESF -fcoarray=lib -lcaf_mpi -ldl -o $PROGRAM

$clean:
    echo Clean Object files ...
    rm -rf *.o *.mod *.out

echo Complete!
