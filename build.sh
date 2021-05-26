#!/bin/bash
shopt -s extglob

. ./build_vars.sh

# Flags
while getopts 'd' flag; do
  case "${flag}" in
    d) DEBUG_API=true;;
    *) error "Unexpected option ${flag}"
       exit 1;;
  esac
done

. ./build_dependencies.sh
. ./build_debug.sh

# Compile SWMM5+
. ./build_swmm5.sh

SOURCESF="$JSON_DIR/json_kinds.F90\
          $JSON_DIR/json_parameters.F90\
          $JSON_DIR/json_string_utilities.F90\
          $JSON_DIR/json_value_module.F90\
          $JSON_DIR/json_file_module.F90\
          $JSON_DIR/json_module.F90\
          $VARS_DIR/data_keys.f08\
          $VARS_DIR/type_definitions.f08\
          $UTIL_DIR/string_utils.f08\
          $INIT_DIR/setting_definition.f08\
          $VARS_DIR/globals.f08\
          $API_DIR/dll.f08\
          $UTIL_DIR/datetime.f08\
          $ARRAY_DIR/dynamic_array.f08\
          $ARRAY_DIR/tables.f08\
          $API_DIR/interface.f08\
          $VARS_DIR/array_index.f08\
          $UTIL_DIR/utility.f08\
          $INIT_DIR/allocate_storage.f08\
          $UTIL_DIR/BIPquick.f08\
          $ARRAY_DIR/coarray_partition.f08\
          $INIT_DIR/initialization.f08\
          $UTIL_DIR/partitioning.f08"

echo
echo Compiling SWMM5+ ...
echo

#caf $SOURCESF $DEBUG_SOURCES main.f08 -ldl -o $PROGRAM
$COARRAY_FC $SOURCESF $DEBUG_SOURCES main.f08 -ldl -o $PROGRAM

$clean:
    echo
    echo Clean Object files ...
    echo
    rm -rf *.o *.mod *.out

echo
echo Complete!
echo
