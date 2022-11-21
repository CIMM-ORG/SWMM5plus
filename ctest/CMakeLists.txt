enable_testing()


add_test(NAME Arch_geometry 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_sub_G=arch_IN=con_OUT=fix.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/temp
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_test(NAME Basket_handle 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_sub_G=basket_handle_IN=con_OUT=fix.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/temp
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_test(NAME Circular 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_sub_G=circ_IN=con_OUT=fix.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/temp
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_test(NAME Egg_shaped 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_sub_G=egg_IN=con_OUT=fix.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/temp
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_test(NAME Filled_circular 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_sub_G=filled_circular_IN=con_OUT=fix.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/temp
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_test(NAME Force_main 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_sub_G=FM_IN=con_OUT=fix.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/temp
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_test(NAME Gothic 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_sub_G=gothic_IN=con_OUT=fix.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/temp
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_test(NAME Horiz_ellipse 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_sub_G=horiz_ellipse_IN=con_OUT=fix.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/temp
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)


add_test(NAME Horse_shoe 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_sub_G=horse_shoe_IN=con_OUT=fix.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/temp
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_test(NAME Mod_basket 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_sub_G=mod_basket_IN=con_OUT=fix.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/temp
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_test(NAME Rect_open
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_sub_G=rect_IN=con_OUT=fix.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/temp
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_test(NAME Rect_round 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_sub_G=rect_round_IN=con_OUT=fix.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/temp
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_test(NAME Semi_circular 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_sub_G=semi_circular_IN=con_OUT=fix.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/temp
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_test(NAME Semi_elliptical 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_sub_G=semi_elliptical_IN=con_OUT=fix.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/temp
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_test(NAME Vert_ellipse 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_sub_G=vert_ellipse_IN=con_OUT=fix.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/temp
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
