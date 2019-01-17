#!/bin/bash

echo Cleaning SWMM ...
DIR=debugoutputA

PROGRAM=SWMM
PRG_OBJ=$PROGRAM.o

$clean:
    echo Clean Object files...
    rm -rf *.o *.mod
    echo Clean debug folder...
    rm -rf "$DIR"
    echo Clean "$PROGRAM" ...
    rm -rf "$PROGRAM"
    
echo Complete!