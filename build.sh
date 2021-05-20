#!/bin/bash
shopt -s extglob

# Flags

DEBUG_API=false
DEBUG_SOURCES=""

while getopts 'd' flag; do
  case "${flag}" in
    d) DEBUG_API=true;;
    *) error "Unexpected option ${flag}"
       exit 1;;
  esac
done

# Constants
SWMM5PLUS_DIR=$PWD

ARRAY_DIR='arrays'
INIT_DIR='initialization'
API_DIR='initialization/interface/flibswmm'
API_TEST_DIR='initialization/interface/tests'
JSON_DIR='json-fortran'
SIM_DIR='simulation'
TEST_DIR='test_cases'
UTIL_DIR='utilities'
VARS_DIR='vars'

# Dependencies source code
MPICH_SOURCE="$SWMM5PLUS_DIR/mpich"
CMAKE_SOURCE="$SWMM5PLUS_DIR/cmake"
COARRAY_SOURCE="$SWMM5PLUS_DIR/opencoarray"

# Dependencies install
MPICH_INSTALL="$MPICH_SOURCE/mpich-install"
CMAKE_INSTALL="$CMAKE_SOURCE/cmake-install"
COARRAY_INSTALL="$COARRAY_SOURCE/opencoarray-install"



COARRAY_FC="${COARRAY_INSTALL}/bin/caf"
FC='gfortran'
OUPTFLAGS=-g
FFLAGS=-O3
PROGRAM=SWMM

# Update settings
TOKEN='"DebugAPI" : '
NEW_VAL="${TOKEN}$DEBUG_API"
sed -i "s#$TOKEN.*#$NEW_VAL,#" $INIT_DIR/settings.json

if [ $DEBUG_API = true ]
then
    rm -r debug
    mkdir debug

    # Compile tests
    cp "$API_TEST_DIR/tests.c" "$API_DIR/src/"
    cp "$API_TEST_DIR/tests.h" "$API_DIR/src/"

    # DEBUG source files
    DEBUG_SOURCES="$API_TEST_DIR/interface_tests.f08"

fi

# Download dependencies
./install.sh


# Compile SWMM C

echo
echo Compiling SWMM 5.13 ...
echo

rm "$API_DIR/src/interface.h"
rm "$API_DIR/src/interface.c"
rm "$API_DIR/src/Makefile"

if [ -f "$API_DIR/src/tests.c" ]
then
    rm "$API_DIR/src/tests.c"
fi

if [ -f "$API_DIR/src/tests.h" ]
then
    rm "$API_DIR/src/tests.h"
fi

cp "$API_DIR/interface.h" "$API_DIR/src/"
cp "$API_DIR/interface.c" "$API_DIR/src/"

# Insert new files in SWMM C Makefile

SCRIPTS="interface.o"
OBJECTS="interface.o   : headers.h interface.h\n"

if [ $DEBUG_API = true ]
then
    SCRIPTS="${SCRIPTS} tests.o"
    OBJECTS="${OBJECTS}tests.o       : headers.h tests.h\n"
    cp "$API_TEST_DIR/tests.h" "$API_DIR/src/"
    cp "$API_TEST_DIR/tests.c" "$API_DIR/src/"
fi

sed "s#{{SCRIPTS}}#$SCRIPTS#" "$API_DIR/../Makefile" > "$API_DIR/src/Makefile"
sed -i "s#{{OBJECTS}}#$OBJECTS#" "$API_DIR/src/Makefile"

cd "$API_DIR/src" && make
cd .. && cp src/libswmm5.so ../../../libswmm5.so
cd ../../..

# Compile SWMM5+

SOURCESF="$JSON_DIR/json_kinds.F90\
          $JSON_DIR/json_parameters.F90\
          $JSON_DIR/json_string_utilities.F90\
          $JSON_DIR/json_value_module.F90\
          $JSON_DIR/json_file_module.F90\
          $JSON_DIR/json_module.F90\
          $VARS_DIR/data_keys.f08\
          $VARS_DIR/type_definitions.f08\
          $UTIL_DIR/string_utils.f08\
          $INIT_DIR/setting_definition.f08\
          $VARS_DIR/globals.f08\
          $API_DIR/dll.f08\
          $UTIL_DIR/datetime.f08\
          $ARRAY_DIR/dynamic_array.f08\
          $ARRAY_DIR/tables.f08\
          $API_DIR/interface.f08\
          $VARS_DIR/array_index.f08\
          $UTIL_DIR/utility.f08\
          $INIT_DIR/allocate_storage.f08\
          $UTIL_DIR/BIPquick.f08\
          $ARRAY_DIR/coarray_partition.f08\
          $INIT_DIR/initialization.f08\
          $UTIL_DIR/partitioning.f08"

echo
echo Compiling SWMM5+ ...
echo

#caf $SOURCESF $DEBUG_SOURCES main.f08 -ldl -o $PROGRAM
$COARRAY_FC $SOURCESF $DEBUG_SOURCES main.f08 -ldl -o $PROGRAM


$clean:
    echo
    echo Clean Object files ...
    echo
    rm -rf *.o *.mod *.out

echo
echo Complete!
echo
