enable_testing()


add_test(NAME Force_main 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_sub_G=FM_IN=con_OUT=fix.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/temp
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
