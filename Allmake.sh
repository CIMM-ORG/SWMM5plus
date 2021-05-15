#!/bin/bash

# echo Download SWMM C interface
# git clone git@github.austin.utexas.edu:gr24269/SWMMwrapper.git
cd SWMMwrapper
./Allmake.sh
cp libswmm5.so ../
cd ..
# sudo rm -r SWMMwrapper

echo Making the debug directory ...
DIRDebug=debugoutputA
DIRThreaded=OutputThreaded
mkdir "$DIRDebug"
mkdir "$DIRThreaded"

echo Compiling the SWMMengine ...

# Compiler/Linker settings
FC=gfortran-9
OPTFLAGS=-g
FFLAGS=-O3
PROGRAM=SWMM
PRG_OBJ=$PROGRAM.o

echo Compiling with $FC ...

# Find all source files, create a list of corresponding object files
echo Find all source files, create a list of corresponding object files ...

# Find all source files, create a list of corresponding object files
SOURCESF="  errors.f08  datetime.f08  dll_mod.f08  type_definitions.f08  dynamic_array.f08\
        setting_definition.f08  globals.f08  selectors.f08  array_index.f08 data_keys.f08\
        objects.f08  tables.f08  interface.f08  utility.f08 xsect_tables.f08 bc.f08  control.f08\
        inflow.f08 network_graph.f08  adjustments.f08  allocate_storage.f08  read_width_depth.f08\
        checking.f08  debug.f08  diagnostic.f08  junction.f08  face_values.f08  storage.f08\
        weir.f08  orifice.f08  element_dynamics.f08  element_geometry.f08  explicit_euler.f08\
        friction_model.f08  initial_condition.f08  initialization.f08  project.f08 link_node.f08\
        network_define.f08  output.f08  artificial_compressibility.f08 runge_kutta.f08  ac_convergence_loop.f08\
        time_loop.f08 stub.f08  case_simple_channel.f08  case_simple_orifice.f08  case_simple_weir.f08\
        case_waller_creek.f08  case_y_channel.f08  case_y_storage_channel.f08  case_simple_pipe.f08 case_swashes.f08\
        trajkovic_cases.f08 weir_steady_state_test.f08 test_cases.f08 main.f08"

# Linker
echo Compiling ...
$FC $SOURCESF -ldl -o $PROGRAM

$clean:
    echo Clean Object files...
    rm -rf *.o *.mod

echo Complete!