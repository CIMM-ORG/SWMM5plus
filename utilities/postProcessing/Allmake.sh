#!/bin/bash

echo Compiling Plot ...

# Compiler/Linker settings
FC=gfortran
OPTFLAGS=-g
FFLAGS=-02
PROGRAM=Plot
PRG_OBJ=$PROGRAM.o

# Find all source files, create a list of corresponding object files
SOURCESF="  postProcessing.f08
			main.f08
    	"

# Linker
echo Compiling ...
$FC -o $PROGRAM $SOURCESF

$clean:
    echo Clean Object files...
    rm -rf *.o *.mod
    
echo Complete!