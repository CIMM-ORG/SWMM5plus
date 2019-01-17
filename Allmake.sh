#!/bin/bash

echo Making the debug directory ...
DIR=debugoutputA
mkdir "$DIR"

echo Compiling the SWMMengine ...

# Compiler/Linker settings
FC=gfortran
OPTFLAGS=-g
FFLAGS=-02
PROGRAM=SWMM
PRG_OBJ=$PROGRAM.o

echo Choosing $FC as a compiler ...

# Find all source files, create a list of corresponding object files
echo Find all source files, create a list of corresponding object files ...

# Find all source files, create a list of corresponding object files
SOURCESF="  setting_definition.f08 \
        globals.f08 \
        array_index.f08\
        data_keys.f08\
        utility.f08\
        adjustments.f08\
        allocate_storage.f08\
        bc.f08\
        case_simple_channel.f08\
        checking.f08\
        debug.f08\
        diagnostic.f08\
        element_dynamics.f08\
        element_geometry.f08\
        face_values.f08\
        explicit_euler.f08\
        friction_model.f08\
        junction.f08\
        initial_condition.f08\
        initialization.f08\
        link_node.f08\
        network_define.f08\
        test_cases.f08\
        runge_kutta.f08\
        time_loop.f08\
        main.f08\
        stub.f08\
    "

# Linker
echo Compiling ...
$FC -o $PROGRAM $SOURCESF

$clean:
    echo Clean Object files...
    rm -rf *.o *.mod
    
echo Complete!