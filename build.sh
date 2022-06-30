#!/bin/bash
shopt -s extglob

# --------------------------------------------------------------------------------------

# source /opt/intel/oneapi/setvars.sh 
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
              $UTIL_DIR/utility_crash.f90\
              $DEF_DIR/define_xsect_tables.f90\
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
              $GEO_DIR/triangular_channel.f90\
              $GEO_DIR/circular_conduit.f90\
              $GEO_DIR/irregular_channel.f90\
              $GEO_DIR/storage_geometry.f90
              $GEO_DIR/geometry.f90\
              $TL_DIR/common_elements.f90\
              $TL_DIR/weir_elements.f90\
              $TL_DIR/orifice_elements.f90\
              $TL_DIR/outlet_elements.f90\
              $TL_DIR/pump_elements.f90\
              $TL_DIR/control_hydraulics.f90\
              $TL_DIR/update.f90\
              $TL_DIR/face.f90\
              $TL_DIR/lowlevel_rk2.f90\
              $TL_DIR/boundary_conditions.f90\
              $TL_DIR/diagnostic_elements.f90\
              $TL_DIR/runge_kutta2.f90\
              $TL_DIR/hydrology.f90\
              $TL_DIR/timeloop.f90\
              $INIT_DIR/initial_condition.f90\
              $INIT_DIR/initialization.f90\
              $FIN_DIR/finalization.f90"

# --------------------------------------------------------------------------------------
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


if [[ $compile_fortran = "true" ]]
then
    ifort -coarray=distributed $SOURCE_FILES $TEST_FILES $MAIN_DIR/main.f90 -I/usr/local/hdf5/include \
    -L/usr/local/hdf5/lib /usr/local/hdf5/lib/libhdf5hl_fortran.a /usr/local/hdf5/lib/libhdf5_hl.a \
    /usr/local/hdf5/lib/libhdf5_fortran.a /usr/local/hdf5/lib/libhdf5.a -lm -Wl,-rpath -Wl,/usr/local/hdf5/lib -ldl -o SWMM
fi

echo
echo Clean Object files ...
echo
rm -rf *.o *.mod *.out
if [[ -d debug_input ]]; then rm -rf debug_input; fi
if [[ -d debug_output ]]; then rm -rf debug_output; fi
if [[ -d swmm5_output ]]; then rm -rf swmm5_output; fi

echo
echo Complete!
echo "To update number of processors in the system:"
echo "Please execute >>> source ~/.bashrc"

# source ~/.bashrc

