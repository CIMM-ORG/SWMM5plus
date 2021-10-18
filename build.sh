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
              $UTIL_DIR/utility_string.f90\
              $DEF_DIR/define_types.f90\
              $DEF_DIR/define_api_keys.f90\
              $DEF_DIR/define_globals.f90\
              $DEF_DIR/define_keys.f90\
              $DEF_DIR/define_settings.f90\
              $DEF_DIR/define_indexes.f90\
              $UTIL_DIR/utility.f90\
              $UTIL_DIR/utility_datetime.f90\
              $API_DIR/c_library.f90\
              $API_DIR/interface.f90\
              $UTIL_DIR/utility_allocate.f90\
              $UTIL_DIR/utility_deallocate.f90\
              $UTIL_DIR/utility_profiler.f90\
              $UTIL_DIR/utility_array.f90\
              $UTIL_DIR/utility_debug.f90\
	          $OUT_DIR/output.f90\
              $UTIL_DIR/utility_output.f90\
              $UTIL_DIR/utility_interpolate.f90\
              $UTIL_DIR/utility_files.f90\
              $INIT_DIR/pack_mask_arrays.f90\
              $INIT_DIR/discretization.f90\
              $INIT_DIR/BIPquick.f90\
              $INIT_DIR/partitioning.f90\
              $INIT_DIR/network_define.f90\
              $UTIL_DIR/utility_unit_testing.f90\
              $TL_DIR/adjust.f90\
              $TL_DIR/jump.f90\
              $GEO_DIR/xsect_tables.f90\
              $GEO_DIR/rectangular_channel.f90\
              $GEO_DIR/trapezoidal_channel.f90\
              $GEO_DIR/circular_conduit.f90\
              $GEO_DIR/geometry.f90\
              $TL_DIR/common_elements.f90\
              $TL_DIR/weir_elements.f90\
              $TL_DIR/orifice_elements.f90\
              $TL_DIR/pump_elements.f90\
              $TL_DIR/lowlevel_rk2.f90\
              $TL_DIR/update.f90\
              $TL_DIR/face.f90\
              $TL_DIR/boundary_conditions.f90\
              $TL_DIR/diagnostic_elements.f90\
              $TL_DIR/runge_kutta2.f90\
              $TL_DIR/timeloop.f90\
              $INIT_DIR/initial_condition.f90\
              $INIT_DIR/initialization.f90\
              $FIN_DIR/finalization.f90"

# --------------------------------------------------------------------------------------

TEST_FILES=""
for i in $(find $TEST_DIR -name '*.f90')
do
    fname=$(basename -- "$i")
    fname="${fname%.*}"
    if [[ $fname != _* ]] && [[ $fname != "main" ]]
    then
        TEST_FILES="$TEST_FILES $i"
    fi
done

echo
echo Compiling SWMM5+ ...
echo

if [[ ! $skip_fortran = "true" ]]
then
    $CAF $SOURCE_FILES $TEST_FILES $MAIN_DIR/main.f90 -ldl -o $PROGRAM
fi

# --------------------------------------------------------------------------------------

$clean:
echo
echo Clean Object files ...
echo
rm -rf *.o *.mod *.out
if [[ -d debug_input ]]; then rm -r debug_input; fi
if [[ -d debug_output ]]; then rm -r debug_output; fi
if [[ -d swmm5_output ]]; then rm -r swmm5_output; fi

echo
echo Complete!
echo
