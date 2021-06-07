#!/bin/bash
shopt -s extglob

# --------------------------------------------------------------------------------------

. ./build_vars.sh
. ./build_dependencies.sh

# --------------------------------------------------------------------------------------

SOURCE_FILES="$JSON_DIR/json_kinds.F90\
              $JSON_DIR/json_parameters.F90\
              $JSON_DIR/json_string_utilities.F90\
              $JSON_DIR/json_value_module.F90\
              $JSON_DIR/json_file_module.F90\
              $JSON_DIR/json_module.F90\
              $UTIL_DIR/utility_string.f08\
              $UTIL_DIR/utility_datetime.f08\
              $DEF_DIR/define_keys.f08\
              $DEF_DIR/define_types.f08\
              $DEF_DIR/define_settings.f08\
              $DEF_DIR/define_globals.f08\
              $DEF_DIR/define_indexes.f08\
              $UTIL_DIR/utility.f08\
              $UTIL_DIR/utility_allocate.f08\
              $UTIL_DIR/utility_array.f08\
              $API_DIR/dll.f08\
              $API_DIR/interface.f08\
              $INIT_DIR/discretization.f08\
              $INIT_DIR/partitioning.f08\
              $INIT_DIR/network_define.f08\
              $INIT_DIR/initialization.f08\
              $FIN_DIR/finalization.f08"

# --------------------------------------------------------------------------------------

TEST_FILES=""
for i in $(find $TEST_DIR -name '*.f08')
do
    fname=$(basename -- "$i")
    fname="${fname%.*}"
    if [[ $fname != _* ]]
    then
        TEST_FILES="$TEST_FILES $i"
    fi
done

echo
echo Compiling SWMM5+ ...
echo

$CAF $SOURCE_FILES $TEST_FILES $MAIN_DIR/main.f08 -ldl -o $PROGRAM

# --------------------------------------------------------------------------------------

$clean:
    echo
    echo Clean Object files ...
    echo
    rm -rf *.o *.mod *.out
    if [[ -d debug ]]; then rm -r debug; fi
    mkdir debug

echo
echo Complete!
echo