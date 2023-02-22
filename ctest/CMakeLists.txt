enable_testing()

###################################### GEOMETRY #############################

# # Passed 20230220
# add_test(NAME Circular_FixedBC 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_sub_G=circ_IN=con_OUT=fix.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230220
# add_test(NAME Circular_FreeBC 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_circ_01.inp -s ./ctest/internal_testing/geometry_test/settings_Co_circ_01.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230220
# add_test(NAME Circular_FreeBC 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_circ_02.inp -s ./ctest/internal_testing/geometry_test/settings_Co_circ_02.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#)

###################################### LEVEL 1 PASSED #############################
#  passed 20230222
# add_test(NAME extran1_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_final.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230222
# add_test(NAME extran1_arch_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_arch_final.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230222
# add_test(NAME extran1_baskethandle_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_baskethandle_final.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230222
# add_test(NAME extran1_catenary_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_catenary_final.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230222
# add_test(NAME extran1_circular_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_circular_final.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230222
# add_test(NAME extran1_egg_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_egg_final.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230222
# add_test(NAME extran1_filled_circular_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_filled_circular_final.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230223
# add_test(NAME extran1_forcemain_hw_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_forcemain_hw_final.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230222
# add_test(NAME extran1_gothic_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_gothic_final.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# FAILED 20230223
add_test(NAME extran1_horizontal_ellipse_final
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_horizontal_ellipse_final.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

# # FAILED 20230223
# add_test(NAME extran1_horseshoe_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_horseshoe_final.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# FAILED 20230223
# add_test(NAME extran1_parabolic_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_parabolic_final.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230222
# add_test(NAME extran1_rect_open_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_rect_open_final.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230222
# add_test(NAME extran1_rect_open_001
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_rect_open_001.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230222
# add_test(NAME extran1_rect_open_002
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_rect_open_002.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230222
# add_test(NAME extran1_semi_circular_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_semi_circular_final.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230223
# add_test(NAME extran1_semi_elliptical_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_semi_elliptical_final.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230223
# add_test(NAME extran1_vertical_ellipse_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_vertical_ellipse_final.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

###############################################  LEVEL 2 PASSED ###############################

# # Passed 20230104  Failed 20230123
# add_test(NAME extran3_bottom_orifice_pid
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_bottom_orifice_pid.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230104   Failed 20230123
# add_test(NAME extran3_side_open_close_orifice
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_side_open_close_orifice.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230104   Failed 20230123
# add_test(NAME extran3_side_open_close_rect_offset_orifice
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_side_open_close_rect_offset_orifice.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230104   Failed 20230123
# add_test(NAME extran3_side_open_close_rect_orifice
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_side_open_close_rect_orifice.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230104   Failed 20230123
# add_test(NAME extran3_side_orifice
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_side_orifice.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230104   Failed 20230123
# add_test(NAME extran3_side_rect_orifice
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_side_rect_orifice.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230104 Failed 20230123
# add_test(NAME extran4_roadway
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran4_roadway.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230104 Failed 20230123
# add_test(NAME extran4_sideflow
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran4_sideflow.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230104 Failed 20230123
# add_test(NAME extran4_trapezoidal
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran4_trapezoidal.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230104  Failed 20230123
# add_test(NAME extran4_vnotch
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran4_vnotch.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230104  Failed 20230123
# add_test(NAME extran4
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran4.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230102  Failed 20230123
# add_test(NAME extran6_storage_functional_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_storage_functional_final.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # # Passed 20230102 Passed 20230123
# add_test(NAME extran6_pump_type2_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_pump_type2_final.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230102  Failed 20230123
# add_test(NAME extran6_pump_type3_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_pump_type3_final.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230102  Passed 20230123
# add_test(NAME extran6_pump_type4_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_pump_type4_final.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

###################### LEVEL 3 --- WORKING

# add_test(NAME x
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/x.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# NOTHING WORKED! 20230104


###################### FAILURES  -- LEVEL 1 #####################################################

# # FAILS custom NOT WORKING AS OF 20221201
# add_test(NAME extran1_custom_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_custom_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILS (SWMM5+ has nan) AS OF 20221201
# add_test(NAME extran1_forcemain_dw_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_forcemain_dw_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # PASSES BUT NEEDS INVESTIGATION
# add_test(NAME extran1_irregular_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_irregular_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILS with a file error
# add_test(NAME extran1_mod_baskethandle_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_mod_baskethandle_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# ## FAILED (SWMM5+ not implemented) 20222101
# # add_test(NAME extran1_power_mod
# # COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_power_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# # WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# # )

# # FAILED 20230104
# add_test(NAME extran1_rect_closed_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_rect_closed_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230104
# add_test(NAME extran1_rect_round_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_rect_round_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230104
# add_test(NAME extran1_rect_triangular_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_rect_triangular_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED (SWMM5+)  20230104
# add_test(NAME force_main_example_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/force_main_example_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED (SWMM5+)  20230104
# add_test(NAME force_main_gravity_example_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/force_main_gravity_example_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

###################### FAILURES  -- LEVEL 2 #####################################################

# # FAILED 20230104 --- looks like problem in BC -- only link_culvert1_Q has the correct inflow Q
# # while the others appear to immediately increase the flowrate. Probably a problem associated with
# # re-using a BC for multiple inflows
# add_test(NAME all_culverts
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/all_culverts.inp -s ./ctest/internal_testing/level02_test/settings_all_culverts.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )



# # FAILED 20230104  -- something wrong in converting units?
# add_test(NAME extran3_bottom_orifice_si_units
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_bottom_orifice_si_units.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )
# # FAILED 20230104 -- something wrong in converting units?
# add_test(NAME extran3_bottom_rect_orifice_si_units
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_bottom_rect_orifice_si_units.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )
# # FAILED 20230104 -- something wrong in converting units?
# add_test(NAME extran3_side_open_close_orifice_si_units
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_side_open_close_orifice_si_units.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )
# # FAILED 20230104 -- something wrong in converting units?
# add_test(NAME extran3_side_open_close_rect_offset_orifice_si_units
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_side_open_close_rect_offset_orifice_si_units.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )
# # FAILED 20230104 -- something wrong in converting units?
# add_test(NAME extran3_side_open_close_rect_orifice_si_units
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_side_open_close_rect_orifice_si_units.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )
# # FAILED 20230104 -- something wrong in converting units?
# add_test(NAME extran3_side_orifice_si_units
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_side_orifice_si_units.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )
# # FAILED 20230104 -- something wrong in converting units?
# add_test(NAME extran3_side_rect_orifice_si_units
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_side_rect_orifice_si_units.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )




# # FAILED -- depth at orifice is problem -- EPA SWMM wrong?
# # NOTE THAT extran3_side_open_close_orifice depth is OK, so what's going on in EPA SWMM?
# # OR is the problem in compare.py?
# add_test(NAME extran3_bottom_orifice
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_bottom_orifice.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )
# # FAILED -- depth at orifice is problem -- EPA SWMM wrong?
# add_test(NAME extran3_bottom_rect_orifice
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_bottom_rect_orifice.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )



# FAILED -- something strange
# add_test(NAME extran4_no_surcharge
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran4_no_surcharge.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )




# # FAILED 20230102 NEED to set pump depths to zero
# add_test(NAME extran6_pump_type1_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_pump_type1_mod.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Failed -- but looks OK 20230103
# add_test(NAME extran6_control_pump_type2_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_control_pump_type2_final.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/temp_output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

###################### FAILURES  -- LEVEL 3

# # FAILED 20230104 -- crash onconservation issue with JM
# add_test(NAME catchment_as_outfall
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/catchment_as_outfall.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230104 -- crash in compare.py, line 313 with excluded surcharge linke
# add_test(NAME Example1
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/Example1.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230104 -- segmentation fault, problem in BC
# add_test(NAME Example3
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/Example3.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAIL 20230104 -- actually passes, but flowrates look wrong
# # cut simulation time down to 30 minutes and examine more closely
# add_test(NAME gate_control_2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/gate_control_2.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAIL 20230104 -- requires more than 5 upstream connections
# add_test(NAME Greenville_EPA
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/Greenville_EPA.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAIL 20230104 -- crash in compare.py at line 313 dealing with surcharge link
# add_test(NAME user1
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user1.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAIL 20230104 -- crash in network define
# add_test(NAME user2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user2.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAIL 20230104 -- initial conditions problem
# add_test(NAME user3
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user3.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAIL 20230104 -- crash on BC type
# add_test(NAME user4
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user4.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAIL 20230104 -- 
# add_test(NAME user5
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_F08.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# add_test(NAME extran1_rect_open_test_A05
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/extran1_rect_open_test_A05.inp -s ./ctest/internal_testing/brh/settings_extran1_test_A.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )
