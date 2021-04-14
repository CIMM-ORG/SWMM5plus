#!/bin/bash
# Download dependencies
shopt -s extglob
if [ ! -d 'json-fortran' ]
then
    git clone 'https://github.com/jacobwilliams/json-fortran.git'
    cd json-fortran
    rm -r -v !('src'|'LICENSE')
    sudo rm -r .*
    mv src/*.* .
    rm -r src
fi

INIT_DIR='initialization'
JSON_DIR='json-fortran'
VARS_DIR='vars'

FC='gfortran-9'
OUPTFLAGS=-g
FFLAGS=-O3
PROGRAM=SWMM
PROGRAM_OBJ=$PROGRAM.o

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
          main.f08"

echo Compiling ...
$FC $SOURCESF -ldl -o $PROGRAM

$clean:
    echo Clean Object files...
    rm -rf *.o *.mod

echo Complete!