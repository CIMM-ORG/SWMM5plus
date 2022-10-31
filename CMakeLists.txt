include(FetchContent)
include(ExternalProject)

cmake_minimum_required(VERSION 3.16.0)
project(SWMM5Plus VERSION 0.0.1)
enable_language (Fortran)


find_package (Python3 COMPONENTS Interpreter)
set(testing_available ON)

function(check_python_dependencies module)
message(STATUS "Looking for Python module '${module}'")
execute_process(COMMAND pip show ${module}
    OUTPUT_QUIET   
    ERROR_QUIET 
    RESULT_VARIABLE EXIT_CODE_PIPSHOW
    )
if (NOT ${EXIT_CODE_PIPSHOW} EQUAL 0)
    set(testing_available OFF PARENT_SCOPE)
    message(WARNING "\nMissing Python module '${module}'\n**Testing unavailable**\nPlease run python_dependencies.sh\n")
else()
    message(STATUS "'${module}' found")
endif()
endfunction()






# if(NOT _swmmtoolbox_status)
#     set(PY_${module_upper} ${_swmmtoolbox_location} CACHE STRING 
#         "Location of Python module ${module}")
# endif(NOT _swmmtoolbox_status)

if (NOT CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    message(FATAL_ERROR "Intel Fortran compiler not found. Found ${CMAKE_Fortran_COMPILER_ID} instead. Run:
    \n$ source path/to/intel/oneapi/setvars.sh
    \n$ export FC=ifort
    \nRemove all files in build directory: $ rm -rf *
    \nRun cmake again. Exiting.")
endif()

function(commentout IN_FILE line_num)
  file (STRINGS ${IN_FILE} LINES)
  file(WRITE ${IN_FILE} "")
  set(line_counter 1)
  foreach(LINE IN LISTS LINES)
    if (${line_counter} STREQUAL line_num)
        string(PREPEND LINE "//")
    endif()
    MATH(EXPR line_counter "${line_counter}+1")
    file(APPEND ${IN_FILE} "${LINE}\n")
  endforeach()
endfunction()

get_filename_component(build_configuration ${CMAKE_BINARY_DIR} NAME)

if (build_configuration STREQUAL "build")
    set(default_build_type "Release")
    message(STATUS "Building inside '${build_configuration}' named folder, setting configuration as '${default_build_type}'")
elseif (build_configuration STREQUAL "dwflow")
    set(default_build_type "Release")
    message(STATUS "Building inside '${build_configuration}' named folder, setting configuration as '${default_build_type}'")
    message(WARNING "Lines [249-252] will be commented out from 'interface/src/dwflow.c'")
elseif(build_configuration STREQUAL "debug")
    set(default_build_type "Debug")
    message(STATUS "Building inside '${build_configuration}' named folder, setting configuration as '${default_build_type}'")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback -debug extended")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -lm -lpthread -shared -fcommon")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -coarray=distributed -ldl -lz")

set(FETCHCONTENT_QUIET "ON")

set(BUILD_COARRAY_NUM_IMAGES 1 CACHE STRING "Set up the number of processors to use for parallel SWMM5+.")
set(BUILD_HDF5 "/hdf5/hdf5-1.12.1/hdf5")
if (DEFINED ENV{HDF5_DIR})
    set(BUILD_HDF5 $ENV{HDF5_DIR})
endif()




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

if (EXISTS ${BUILD_HDF5}/lib 
    AND EXISTS ${BUILD_HDF5}/include
    AND EXISTS ${BUILD_HDF5}/lib/libhdf5hl_fortran.a
    AND EXISTS ${BUILD_HDF5}/lib/libhdf5_hl.a
    AND EXISTS ${BUILD_HDF5}/lib/libhdf5_fortran.a
    AND EXISTS ${BUILD_HDF5}/lib/libhdf5.a
)

    ExternalProject_Add(hdf5
    SOURCE_DIR ${BUILD_HDF5}
    CONFIGURE_COMMAND ""
    BINARY_DIR ${BUILD_HDF5}
    BUILD_COMMAND ""
    INSTALL_COMMAND "")
    set(HDF5_LIB_DIR ${BUILD_HDF5}/lib)
    set(HDF5_INCLUDE_DIR ${BUILD_HDF5}/include)
    set(HDF5_LIBRARIES ${BUILD_HDF5}/lib/libhdf5hl_fortran.a ${BUILD_HDF5}/lib/libhdf5_hl.a 
    ${BUILD_HDF5}/lib/libhdf5_fortran.a ${BUILD_HDF5}/lib/libhdf5.a)
    set(HDF5_FOUND ON)
endif()

if (HDF5_FOUND)
    message(WARNING "HDF5 found in ${BUILD_HDF5}")
else ()
    message(WARNING "HDF5 not found, will download and build with required options.\n If HDF5 is already installed, point to its root directory by running:\n
$ export HDF5_DIR=path/to/hdf5/")
    FetchContent_Populate(hdf5)
    ExternalProject_Add(hdf5
    SOURCE_DIR ${hdf5_SOURCE_DIR}
    CONFIGURE_COMMAND FC=ifort ${hdf5_SOURCE_DIR}/configure --enable-fortran
    BINARY_DIR ${hdf5_BINARY_DIR}
    BUILD_COMMAND ${MAKE}
    COMMAND make install)
    

    set(HDF5_LIB_DIR ${hdf5_BINARY_DIR}/hdf5/lib/)
    set(HDF5_INCLUDE_DIR ${hdf5_BINARY_DIR}/hdf5/include/)
    set(HDF5_LIBRARIES ${hdf5_BINARY_DIR}/hdf5/lib/libhdf5hl_fortran.a ${hdf5_BINARY_DIR}/hdf5/lib/libhdf5_hl.a 
    ${hdf5_BINARY_DIR}/hdf5/lib/libhdf5_fortran.a ${hdf5_BINARY_DIR}/hdf5/lib/libhdf5.a)

    set(BUILD_HDF5 ${hdf5_BINARY_DIR})
    set(HDF5_FOUND ON)
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
message(STATUS "Appended to controls.c and funcs.h")

if (build_configuration STREQUAL "dwflow")
    commentout(${CMAKE_SOURCE_DIR}/interface/src/dwflow.c 248)
    commentout(${CMAKE_SOURCE_DIR}/interface/src/dwflow.c 249)
    commentout(${CMAKE_SOURCE_DIR}/interface/src/dwflow.c 250)
    commentout(${CMAKE_SOURCE_DIR}/interface/src/dwflow.c 251)
    commentout(${CMAKE_SOURCE_DIR}/interface/src/dwflow.c 252)
    message(STATUS "Commented out lines [248-252] from 'interface/src/dwflow.c'")
endif()


file(GLOB
    SWMM_SOURCES
        ${CMAKE_SOURCE_DIR}/interface/src/*.c 
        ${CMAKE_SOURCE_DIR}/interface/src/*.h
        # ${CMAKE_SOURCE_DIR}/*.h
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
        # ${CMAKE_SOURCE_DIR}/interface/*.f90
        ${CMAKE_SOURCE_DIR}/interface/c_library.f90
        ${CMAKE_SOURCE_DIR}/interface/interface.f90
        ${CMAKE_SOURCE_DIR}/initialization/*.f90
        ${CMAKE_SOURCE_DIR}/geometry/*.f90
        ${CMAKE_SOURCE_DIR}/timeloop/*.f90
        ${CMAKE_SOURCE_DIR}/output/*.f90
        ${CMAKE_SOURCE_DIR}/finalization/*.f90
        ${CMAKE_SOURCE_DIR}/main/*.f90
)

#Targets

add_library(swmm5 SHARED 
    ${SWMM_SOURCES}
)
set_target_properties(swmm5 PROPERTIES LINKER_LANGUAGE C
CXX_STANDARD 14)


add_executable(SWMM
    ${SWMM5X_SOURCES}
)
target_include_directories(SWMM PUBLIC
    ${HDF5_INCLUDE_DIR}
)
target_link_directories(SWMM PUBLIC
    ${HDF5_LIB_DIR}
)
target_link_libraries(SWMM PUBLIC
    ${HDF5_LIBRARIES}
)
add_dependencies(SWMM hdf5)

# if (DEFINED ENV{FOR_COARRAY_NUM_IMAGES})
#     message(WARNING "Environment variable FOR_COARRAY_NUM_IMAGES is set. Compiler will ignore CMake settings for -coarray-num-images. If unsure, run the following before using CMake \n$ unset FOR_COARRAY_NUM_IMAGES")
# endif()


#Testing
check_python_dependencies(tornado)
check_python_dependencies(tabulate)
check_python_dependencies(h5py)
check_python_dependencies(swmmtoolbox)

#Tests defined in separate file: ctest/CMakeLists.txt
if (testing_available)
    message(STATUS "Making tests available")
    add_subdirectory(ctest)    
    enable_testing()
endif()
