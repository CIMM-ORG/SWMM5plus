enable_testing()

add_test(NAME circular_01
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/circular_01.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

