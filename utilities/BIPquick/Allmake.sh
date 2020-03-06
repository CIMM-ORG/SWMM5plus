#!/bin/bash

echo Compiling BIPquick ...

# Compiler/Linker settings
FC=gfortran
OPTFLAGS=-g
FFLAGS=-02
PROGRAM=BIPquick
PRG_OBJ=$PROGRAM.o

# Find all source files, create a list of corresponding object files
SOURCESF="  BIPquick.f08
    "

# Linker
echo Compiling ...
$FC -o $PROGRAM $SOURCESF

$clean:
    echo Clean Object files...
    rm -rf *.o *.mod
    
echo Complete!