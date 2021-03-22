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