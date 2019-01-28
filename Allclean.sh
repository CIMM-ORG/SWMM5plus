#!/bin/bash

echo Cleaning SWMM ...
DIRDebug=debugoutputA
DIRThreaded=OutputThreaded

PROGRAM=SWMM
PRG_OBJ=$PROGRAM.o

$clean:
    echo Clean Object files...
    rm -rf *.o *.mod
    echo Clean debug folder...
    rm -rf "$DIRDebug"
    echo Clean threaded folder...
    rm -rf "$DIRThreaded"
    echo Clean "$PROGRAM" ...
    rm -rf "$PROGRAM"
    
echo Complete!