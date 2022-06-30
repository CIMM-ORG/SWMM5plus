include(FetchContent)
include(ExternalProject)

cmake_minimum_required(VERSION 3.20.0)
project(SWMM5Plus VERSION 0.0.1)
enable_language (Fortran)

if (NOT CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    message(FATAL_ERROR "Intel Fortran compiler not found. Exiting.")
endif()

set(HDF5_PREFER_PARALLEL true)
find_package(HDF5 "1.12.2" COMPONENTS Fortran)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -lm -lpthread -shared")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -coarray=distributed -ldl")
set(FETCHCONTENT_QUIET "OFF")

set(BUILD_COARRAY_NUM_IMAGES $ENV{FOR_COARRAY_NUM_IMAGES} CACHE STRING "Set up the number of processors to use for parallel SWMM5+.")



FetchContent_Declare(
    epaswmm5
    URL https://github.com/USEPA/Stormwater-Management-Model/archive/v5.1.13.tar.gz
)

FetchContent_Declare(
    json_fortran
    GIT_REPOSITORY https://github.com/jacobwilliams/json-fortran.git
    GIT_TAG 3ab8f98209871875325c6985dd0e50085d1c82c2 #Release 8.3.0
)
FetchContent_Declare(
    hdf5
    GIT_REPOSITORY https://github.com/HDFGroup/hdf5.git
    GIT_TAG 3e847e003632bdd5fdc189ccbffe25ad2661e16f #Release 1.12.2
)
Fetchcontent_MakeAvailable(epaswmm5)
FetchContent_GetProperties(json_fortran hdf5)


if(NOT json_fortran_POPULATED)
  FetchContent_Populate(json_fortran)
endif()

if(NOT HDF5_FOUND)
    message(WARNING "HDF5 not found. Will download and install locally.")
    if(NOT hdf5_POPULATED)
        FetchContent_Populate(hdf5)
        ExternalProject_Add(hdf5
        SOURCE_DIR ${hdf5_SOURCE_DIR}
        CONFIGURE_COMMAND ${hdf5_SOURCE_DIR}/configure --enable-fortran
        BINARY_DIR ${hdf5_BINARY_DIR}
        BUILD_COMMAND ${MAKE})
    endif()
endif()

#Following sh script
file(COPY ${epaswmm5_SOURCE_DIR}/src/ DESTINATION ${CMAKE_SOURCE_DIR}/interface/src)
file(COPY ${CMAKE_SOURCE_DIR}/interface/api.h DESTINATION ${CMAKE_SOURCE_DIR}/interface/src)
file(COPY ${CMAKE_SOURCE_DIR}/interface/api.c DESTINATION ${CMAKE_SOURCE_DIR}/interface/src)
file(COPY ${CMAKE_SOURCE_DIR}/interface/api_error.h DESTINATION ${CMAKE_SOURCE_DIR}/interface/src)
file(COPY ${CMAKE_SOURCE_DIR}/interface/api_error.c DESTINATION ${CMAKE_SOURCE_DIR}/interface/src)

file(READ ${CMAKE_SOURCE_DIR}/interface/add_to_controls.c CONTROLS)
file(READ ${CMAKE_SOURCE_DIR}/interface/add_to_funcs.h FUNCS)
file(APPEND ${CMAKE_SOURCE_DIR}/interface/src/controls.c "${CONTROLS}")
file(APPEND ${CMAKE_SOURCE_DIR}/interface/src/funcs.h "${FUNCS}")

file(GLOB
    SWMM_SOURCES
        ${CMAKE_SOURCE_DIR}/interface/src/*.c 
        ${CMAKE_SOURCE_DIR}/*.h
)
file(GLOB
    SWMM5X_SOURCES
        ${json_fortran_SOURCE_DIR}/src/json_kinds.F90
        ${json_fortran_SOURCE_DIR}/src/json_parameters.F90
        ${json_fortran_SOURCE_DIR}/src/json_string_utilities.F90
        ${json_fortran_SOURCE_DIR}/src/json_value_module.F90
        ${json_fortran_SOURCE_DIR}/src/json_file_module.F90
        ${json_fortran_SOURCE_DIR}/src/json_module.F90

        ${CMAKE_SOURCE_DIR}/utility/*.f90 
        ${CMAKE_SOURCE_DIR}/definitions/*.f90
        ${CMAKE_SOURCE_DIR}/interface/*.f90
        ${CMAKE_SOURCE_DIR}/initialization/*.f90
        ${CMAKE_SOURCE_DIR}/geometry/*.f90
        ${CMAKE_SOURCE_DIR}/timeloop/*.f90
        ${CMAKE_SOURCE_DIR}/output/*.f90
        ${CMAKE_SOURCE_DIR}/finalization/*.f90
        ${CMAKE_SOURCE_DIR}/main/*.f90

        ${hdf5_BINARY_DIR}/hdf5/include/
        ${hdf5_BINARY_DIR}/hdf5/lib/
        ${hdf5_BINARY_DIR}/hdf5/lib/libhdf5hl_fortran.a
        ${hdf5_BINARY_DIR}/hdf5/lib/libhdf5_hl.a
        ${hdf5_BINARY_DIR}/hdf5/lib/libhdf5_fortran.a
        ${hdf5_BINARY_DIR}/hdf5/lib/libhdf5.a
)

#Targets

add_library(swmm5 SHARED 
    ${SWMM_SOURCES}
)
set_target_properties(swmm5 PROPERTIES LINKER_LANGUAGE C)


add_executable(SWMM
    ${SWMM5X_SOURCES}
)
target_compile_options(SWMM PUBLIC "$<$<COMPILE_LANGUAGE:Fortran>:-coarray-num-images=${BUILD_COARRAY_NUM_IMAGES}>")
add_dependencies(SWMM hdf5)

if (DEFINED ENV{FOR_COARRAY_NUM_IMAGES})
    message(WARNING "Environment variable FOR_COARRAY_NUM_IMAGES is set. Compiler will ignore CMake settings for -coarray-num-images. If unsure, run the following before using CMake \n$ unset FOR_COARRAY_NUM_IMAGES")
endif()


#Testing
enable_testing()
add_test(A_Single_Link SWMM
${CMAKE_SOURCE_DIR}/test_cases/A_Single_Link/SL_sub_IN\=con_OUT\=fix.inp
-s ${CMAKE_SOURCE_DIR}/test_cases/A_Single_Link/settings_SL.json
)
