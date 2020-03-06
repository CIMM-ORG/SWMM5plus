#!/bin/bash

echo downloading the Stormwater Management Model from US EPA repository ...
git clone "https://github.com/USEPA/Stormwater-Management-Model.git"

echo Copying the source files and combine it with interface code ...
cp -R Stormwater-Management-Model/src/. src/

echo removing the unnecessary folders and files ...
rm -rf "Stormwater-Management-Model"

#removing the copied c files
mkdir tmp
mv src/interface.h /tmp/
mv src/interface.c /tmp/
mv src/Makefile /tmp/

rm -rf src

mkdir src
mv /tmp/interface.h src/
mv /tmp/interface.c src/
mv /tmp/Makefile src/
rm -rf tmp