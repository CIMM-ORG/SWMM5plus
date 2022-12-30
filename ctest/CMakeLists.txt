enable_testing()

# add_test(NAME extran1_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME extran1_arch_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_arch_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME extran1_baskethandle_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_baskethandle_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME extran1_catenary_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_catenary_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME extran1_circular_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_circular_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # # FAILS custom NOT WORKING AS OF 20221201
# # add_test(NAME extran1_custom_mod
# # COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_custom_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# # WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# # )

# add_test(NAME extran1_egg_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_egg_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# ## FAILS (SWMM5+ ERROR) AS OF 20221201
# # add_test(NAME extran1_filled_circular_mod
# # COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_filled_circular_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# # WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# # )

# ## FAILS (SWMM5+ has nan) AS OF 20221201
# # add_test(NAME extran1_forcemain_dw_mod
# # COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_forcemain_dw_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# # WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# # )

# ## FAILS (SWMM5+ gets an error in reading EPA SWMM input)
# # add_test(NAME extran1_forcemain_hw_mod
# # COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_forcemain_hw_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# # WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# # )

# add_test(NAME extran1_gothic_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_gothic_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME extran1_horizontal_ellipse_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_horizontal_ellipse_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME extran1_horseshoe_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_horseshoe_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# FAILED (SWMM5+ SEGMENTATION) 20221201
# OK based on figures 20221219
# add_test(NAME extran1_irregular_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_irregular_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME extran1_modified_baskethandle_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_modified_baskethandle_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME extran1_parabolic_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_parabolic_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# ## FAILED (SWMM5+ not implemented) 20222101
# # add_test(NAME extran1_power_mod
# # COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_power_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# # WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# # )

# add_test(NAME extran1_rect_closed_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_rect_closed_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME extran1_rect_open_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_rect_open_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME extran1_rect_round_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_rect_round_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME extran1_rect_triangular_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_rect_triangular_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME extran1_semi_circular_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_semi_circular_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME extran1_semi_elliptical_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_semi_elliptical_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME extran1_vertical_ellipse_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_vertical_ellipse_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# ## FAILED (SWMM5+)  20221201
# # add_test(NAME force_main_example_mod
# # COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/force_main_example_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# # WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# # )

# # # FAILED (SWMM5+)  20221201
# # add_test(NAME force_main_gravity_example_mod
# # COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/force_main_gravity_example_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/temp/
# # WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# # )

# # FAILED (splitting downstream mismatch)  20221221
# add_test(NAME extran6_type2_pump_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_type2_pump_mod.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/temp
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED (splitting downstream mismatch) 20221222
# add_test(NAME extran6_type3_pump_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_type3_pump_mod.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/temp
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED (splitting downstream mismatch)  20221222
# add_test(NAME extran6_type4_pump_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_type4_pump_mod.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/temp
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED (Functional storage mismatch) 20221221
# add_test(NAME extran6_storage_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_storage_mod.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/temp
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # # NOT WORKING 20221227
# add_test(NAME extran6_control_type2_pump_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_control_type2_pump_mod.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/temp
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Failed, but should be fixed with python head fix at junction 20221230
# add_test(NAME extran6_storage_functional_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_storage_functional_final.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/temp
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # # 
# add_test(NAME extran6_pump_type1_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_pump_type1_mod.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/temp
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # 
add_test(NAME extran6_loop_mod
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_loop_mod.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/temp
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

