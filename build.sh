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

SOURCESF="$INIT_DIR/dll.f08\
            $INIT_DIR/setting_definition.f08\
            $INIT_DIR/interface.f08"

echo Compiling ...
$FC $SOURCESF -ldl -o $PROGRAM

$clean:
    echo Clean Object files...
    rm -rf *.o *.mod

echo Complete!