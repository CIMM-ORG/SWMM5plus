enable_testing()

add_test(NAME A_Single_Link 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./test_cases/A_Single_Link/SL_sub_IN\=con_OUT\=fix.inp -s ./test_cases/A_Single_Link/settings_SL.json -o ./ctest
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)