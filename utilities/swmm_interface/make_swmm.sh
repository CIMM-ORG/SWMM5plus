#!/bin/bash

echo downloading the Stormwater Management Model from US EPA repository ...
wget "https://github.com/USEPA/Stormwater-Management-Model/archive/v5.1.13.tar.gz"
tar -xvf *.tar.gz
rm *.tar.gz

mv Stormwater*/src .
rm -r Stormwater*

cp interface.h src/
cp interface.c src/
cp Makefile src/
cd src

echo Compiling SWMM C code
make
cd ..

echo Generating DLL
cp src/libswmm5.so ../../
rm -r src