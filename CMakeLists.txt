include(FetchContent)
include(ExternalProject)

cmake_minimum_required(VERSION 3.16.0)
project(SWMM5Plus VERSION 1.0.3)
enable_language (Fortran)


find_package (Python3 COMPONENTS Interpreter)

set(testing_available ON)
set(CPACK_SET_DESTDIR ON)

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
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -traceback -coarray=distributed -ldl")

set(FETCHCONTENT_QUIET OFF)

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
    GIT_PROGRESS ON
)
FetchContent_Declare(
    hdf5
    GIT_REPOSITORY https://github.com/HDFGroup/hdf5.git
    GIT_TAG 3e847e003632bdd5fdc189ccbffe25ad2661e16f #Release 1.12.2
    GIT_PROGRESS ON
)
Fetchcontent_MakeAvailable(epaswmm5)
FetchContent_GetProperties(json_fortran hdf5)


if(NOT json_fortran_POPULATED)
  FetchContent_Populate(json_fortran)
endif()

FetchContent_Populate(hdf5)
ExternalProject_Add(hdf5
SOURCE_DIR ${hdf5_SOURCE_DIR}
CONFIGURE_COMMAND FC=ifort ${hdf5_SOURCE_DIR}/configure --enable-fortran --without-zlib
BINARY_DIR ${hdf5_BINARY_DIR}
BUILD_COMMAND ${MAKE}
COMMAND make install)

set(HDF5_LIB_DIR ${hdf5_BINARY_DIR}/hdf5/lib/)
set(HDF5_INCLUDE_DIR ${hdf5_BINARY_DIR}/hdf5/include/)
set(HDF5_LIBRARIES ${hdf5_BINARY_DIR}/hdf5/lib/libhdf5hl_fortran.a ${hdf5_BINARY_DIR}/hdf5/lib/libhdf5_hl.a 
${hdf5_BINARY_DIR}/hdf5/lib/libhdf5_fortran.a ${hdf5_BINARY_DIR}/hdf5/lib/libhdf5.a)

set(BUILD_HDF5 ${hdf5_BINARY_DIR})
set(HDF5_FOUND ON)

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

file(READ ${CMAKE_SOURCE_DIR}/interface/add_to_lid.c LIDC)
file(READ ${CMAKE_SOURCE_DIR}/interface/add_to_lid.h LIDH)
file(APPEND ${CMAKE_SOURCE_DIR}/interface/src/lid.c "${LIDC}")
file(APPEND ${CMAKE_SOURCE_DIR}/interface/src/lid.h "${LIDH}")
message(STATUS "Appended to lid.c and lid.h")

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
        ${CMAKE_SOURCE_DIR}/special_elements/*.f90 
        ${CMAKE_SOURCE_DIR}/definitions/*.f90
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


add_executable(SWMM5p
    ${SWMM5X_SOURCES}
)
target_include_directories(SWMM5p PUBLIC
    ${HDF5_INCLUDE_DIR}
)
target_link_directories(SWMM5p PUBLIC
    ${HDF5_LIB_DIR}
)
target_link_libraries(SWMM5p PUBLIC
    ${HDF5_LIBRARIES}
)
add_dependencies(SWMM5p hdf5)

#Testing
check_python_dependencies(tornado)
check_python_dependencies(tabulate)
check_python_dependencies(h5py)
check_python_dependencies(swmmtoolbox)

# Tests defined in separate file: ctest/CMakeLists.txt
if (testing_available)
    enable_testing()
    message(STATUS "Making tests available")
    add_subdirectory(ctest)    
endif()


#INSTALLATION
install(TARGETS SWMM5p DESTINATION bin)
install(TARGETS swmm5 DESTINATION lib)

#cpack
set(CPACK_PACKAGE_CONTACT "NCIMM")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "NCIMM")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "intel-hpckit")
include(CPack)
