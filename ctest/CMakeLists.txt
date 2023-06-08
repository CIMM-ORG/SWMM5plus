enable_testing()

##########################################################################################
## Final tests for code cleanup 20230601
##########################################################################################

add_test(NAME Geometry_test 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/release_tests/all_geometry_test.inp -s ./ctest/internal_testing/release_tests/settings_dx10.json -o ./ctest/outbrh/
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_test(NAME Y_Junction_test 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/release_tests/two_Y_fixed_inflow.inp -s ./ctest/internal_testing/release_tests/settings_dx10.json -o ./ctest/outbrh/
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_test(NAME Weir_test 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/release_tests/all_weir_test.inp -s ./ctest/internal_testing/release_tests/settings_dx10.json -o ./ctest/outbrh/
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_test(NAME Orifice_test 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/release_tests/all_orifice_test.inp -s ./ctest/internal_testing/release_tests/settings_dx10.json -o ./ctest/outbrh/
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_test(NAME Storage_Unit_test 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/release_tests/all_storage_emptying_test.inp -s ./ctest/internal_testing/release_tests/settings_dx10.json -o ./ctest/outbrh/
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)



# PUMP 1 FAILS
# add_test(NAME Pump test 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/finalized_tests/all_pump_test.inp -s ./ctest/internal_testing/finalized_tests/settings_dx5.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )




###############3333
## PROBLEM 20230326
## T002_RC_Fix_dx0100 is smooth
## T002_CC_Fix_dx0100 has periodic oscillations, WHY?

##########################################################################################
## T001 downstream boundary condition tests
##########################################################################################
## single link of 1000 m 
## no internal junctions, 
## initial zero water levels
## ramp up flow
## end at steady-state flow value
## no surcharge
## CC = circular conduit; RC = Rectangular closed, RO = Rectangular open
## Fix/Normal/Free = BC
## dxNNN = subelement length

# # passed 20230506
# add_test(NAME T001_CC_Fix_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_CC_Fix_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T001_CC_Free_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_CC_Free_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T001_CC_Normal_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_CC_Normal_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # passed 20230506
# add_test(NAME T001_RC_Free_dx0500 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_RC_Free_dx0500.inp -s ./ctest/internal_testing/brh/settings_500.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T001_RC_Free_dx0333 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_RC_Free_dx0333.inp -s ./ctest/internal_testing/brh/settings_333.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T001_RC_Free_dx0250 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_RC_Free_dx0250.inp -s ./ctest/internal_testing/brh/settings_250.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T001_RC_Free_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_RC_Free_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T001_RO_Free_dx0500 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_RO_Free_dx0500.inp -s ./ctest/internal_testing/brh/settings_500.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T001_RO_Free_dx0250 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_RO_Free_dx0250.inp -s ./ctest/internal_testing/brh/settings_250.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T001_RO_Free_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_RO_Free_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T001_RO_Free_dx050 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_RO_Free_dx0050.inp -s ./ctest/internal_testing/brh/settings_50.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )



##########################################################################################
## T002 junction tests
##########################################################################################
## Single flow path
## two links of 500 m 
## one internal junction, 
## initial zero water levels
## ramp up flow
## end at steady-state flow value
## no surcharge
## CC = circular conduit; RC = Rectangular closed, RO = Rectangular open
## Fix/Normal/Free = BC
## dxNNN = subelement length

    ###### FIXED DOWNSTREAM BOUNDARY #####

# # Passed 20230506
# add_test(NAME T002_CC_Fix_dx0010_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Fix_dx0010_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_10_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T002_CC_Fix_dx0100_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Fix_dx0100_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_100_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T002_CC_Fix_dx0250_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Fix_dx0250_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_250_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

    ###### FREE DOWNSTREAM BOUNDARY #####

# Passed 20230506
# add_test(NAME T002_CC_Free_dx0010_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Free_dx0010_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_10_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T002_CC_Free_dx0100_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Free_dx0100_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_100_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T002_CC_Free_dx0250_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Free_dx0250_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_250_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # Passed 20230506
# add_test(NAME T002_RC_Fix_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_RC_Fix_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T002_RO_Free_dx0250 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_RO_Free_dx0250.inp -s ./ctest/internal_testing/brh/settings_250.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T002_RO_Free_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_RO_Free_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )



# # Passed 20230506
# add_test(NAME T002_RC_Free_dx0250 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_RC_Free_dx0250.inp -s ./ctest/internal_testing/brh/settings_250.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T002_RC_Free_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_RC_Free_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T002_RC_Free_dx0010 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_RC_Free_dx0010.inp -s ./ctest/internal_testing/brh/settings_10.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # Passed 20230506
# add_test(NAME T002_CC_Fix_dx0250 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Fix_dx0250.inp -s ./ctest/internal_testing/brh/settings_250.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T002_CC_Fix_dx0167 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Fix_dx0167.inp -s ./ctest/internal_testing/brh/settings_167.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T002_CC_Fix_dx0125 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Fix_dx0125.inp -s ./ctest/internal_testing/brh/settings_125.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T002_CC_Fix_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Fix_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T002_CC_Fix_dx0010 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Fix_dx0010.inp -s ./ctest/internal_testing/brh/settings_10.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T002_CC_Free_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Free_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T002_CC_Free_dx0010 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Free_dx0010.inp -s ./ctest/internal_testing/brh/settings_10.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T002_CC_Normal_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Normal_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

##########################################################################################
## T003 junction tests
##########################################################################################
## Single Flow paths two junctions
## two links of 500 m 
## two internal junction, 
## initial zero water levels
## ramp up flow
## end at steady-state flow value
## no surcharge
## CC = circular conduit; RC = Rectangular closed, RO = Rectangular open
## Fix/Normal/Free = BC
## dxNNN = subelement length

    ###### FIXED DOWNSTREAM BOUNDARY #####


# # Passed 20230506
# add_test(NAME T003_CC_Fix_dx0010_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_CC_Fix_dx0010_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_10_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T003_CC_Fix_dx0111_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_CC_Fix_dx0111_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_111_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )



    ##### FREE DOWNSTREAM BOUNDARY #####


# #  Passed 20230506
# add_test(NAME T003_CC_Free_dx0010_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_CC_Free_dx0010_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_10_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# #  Passed 20230506
# add_test(NAME T003_CC_Free_dx0111_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_CC_Free_dx0111_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_111_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # passed 20230506
# add_test(NAME T003_RC_Free_dx0167 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_RC_Free_dx0167.inp -s ./ctest/internal_testing/brh/settings_167.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 202304506
# add_test(NAME T003_RC_Free_dx0111 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_RC_Free_dx0111.inp -s ./ctest/internal_testing/brh/settings_111.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T003_RC_Free_dx0048 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_RC_Free_dx0048.inp -s ./ctest/internal_testing/brh/settings_48.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# ## passed 20230506
# add_test(NAME T003_CC_Free_dx0048 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_CC_Free_dx0048.inp -s ./ctest/internal_testing/brh/settings_48.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T003_CC_Free_dx0010 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_CC_Free_dx0010.inp -s ./ctest/internal_testing/brh/settings_10.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )



##########################################################################################
## T004 junction tests
##########################################################################################
## Single flow path
## Similar to T001 to T003, but with 9 junctions
## ten links of 100 m 
## nine internal junction, 
## initial zero water levels
## ramp up flow
## end at steady-state flow value
## no surcharge
## CC = circular conduit; RC = Rectangular closed, RO = Rectangular open
## Fix/Normal/Free = BC
## dxNNN = subelement length

    ###### FREE DOWNSTREAM BOUNDARY #####

# # Passed 20230506
# add_test(NAME T004_CC_Free_dx0010_NJ2 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T004_CC_Free_dx0010_NJ2.inp -s ./ctest/internal_testing/brh/settings_10_faceNJ2.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T004_CC_Free_dx0010_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T004_CC_Free_dx0010_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_10_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T004_CC_Free_dx0050_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T004_CC_Free_dx0050_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_50_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


    ###### FIXED DOWNSTREAM BOUNDARY #####         

# # Passed 20230506
# add_test(NAME T004_CC_Fix_dx0010_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T004_CC_Fix_dx0010_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_10_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T004_CC_Free_dx0010 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T004_CC_Free_dx0010.inp -s ./ctest/internal_testing/brh/settings_10.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T004_CC_Free_dx0050 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T004_CC_Free_dx0050.inp -s ./ctest/internal_testing/brh/settings_50.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )



##########################################################################################
## T005 junction tests
##########################################################################################
## Similar to T004, but with flow ramp down
## ten links of 100 m 
## nine internal junction, 
## initial zero water levels
## ramp up flow
## ramp down flow after 2 hour steady
## no surcharge
## CC = circular conduit; RC = Rectangular closed, RO = Rectangular open
## Fix/Normal/Free = BC
## dxNNN = subelement length

# # passed 20230506
# add_test(NAME T005_CC_Free_dx0010_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T005_CC_Free_dx0010_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_10_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T005_CC_Free_dx0050_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T005_CC_Free_dx0050_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_50_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T005_CC_Free_dx0010_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T005_CC_Free_dx0010_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_10_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T005_CC_Free_dx0050_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T005_CC_Free_dx0050_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_50_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )



##########################################################################################
## T005W weir and junction tests
##########################################################################################
## Similar to T005, but with weir in the middle
## using ForceJM to test junctions with weir
## ten links of 100 m 
## nine internal junction, 
## initial zero water levels
## ramp up flow
## ramp down flow after 2 hour steady
## no surcharge
## CC = circular conduit; RC = Rectangular closed, RO = Rectangular open
## Fix/Normal/Free = BC
## dxNNN = subelement length


# # passed 20230506
# add_test(NAME T005Wr_RO_Free_dx0010_ForceJM 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T005Wr_RO_Free_dx0010_ForceJM.inp -s ./ctest/internal_testing/brh/settings_10_ForceJM.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T005Wr_RO_Free_dx0010_weirJM 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T005Wr_RO_Free_dx0010_weirJM.inp -s ./ctest/internal_testing/brh/settings_10_weirJM.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T005Wr_RO_Free_dx0010_NJ2 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T005Wr_RO_Free_dx0010_NJ2.inp -s ./ctest/internal_testing/brh/settings_10_faceNJ2.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


##########################################################################################
## T006 junction tests
##########################################################################################
## Y Junction with 2 upstream branches (5 links each) and 1 downstream branch (5 links)
## slight differences in upstream branch slope
## each link 100 m
## initial zero water levels
## ramp up flow
## ramp down flow after 2 hour steady
## no surcharge
## CC = circular conduit; RC = Rectangular closed, RO = Rectangular open
## Fix/Normal/Free = BC
## dxNNN = subelement length

# # passed 20230506
# add_test(NAME T006_CC_Free_dx0010_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T006_CC_Free_dx0010_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_10_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T006_CC_Free_dx0010_NJ2 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T006_CC_Free_dx0010_NJ2.inp -s ./ctest/internal_testing/brh/settings_10_faceNJ2.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )



###################################### GEOMETRY #############################

# # Passed 20230506
# add_test(NAME Circular_FixedBC 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_sub_G=circ_IN=con_OUT=fix.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # # Passed 20230506
# add_test(NAME Circular_FreeBC_01 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_circ_01.inp -s ./ctest/internal_testing/geometry_test/settings_Co_circ.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME Circular_FreeBC_02 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_circ_02.inp -s ./ctest/internal_testing/geometry_test/settings_Co_circ.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME Circular_FreeBC_03 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_circ_03.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME RectClosed_FreeBC_03 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_rect_closed_03.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME RectOpen_FreeBC_03 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_rect_open_03.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # Passed 20230506
# add_test(NAME YJ_Geom
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/YJ_Geom_Test.inp -s ./ctest/internal_testing/geometry_test/settings_YJ.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME YJ_Geom_mod01
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/YJ_Geom_Test_mod01.inp -s ./ctest/internal_testing/geometry_test/settings_YJ.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # ###################################### LEVEL 1 PASSED #############################

# # passed 20230506
# add_test(NAME extran1_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME extran1_final_mod2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_final_mod2.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # passed 20230506
# add_test(NAME extran1_arch_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_arch_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME extran1_baskethandle_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_baskethandle_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230506
# add_test(NAME extran1_catenary_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_catenary_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME extran1_circular_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_circular_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME extran1_egg_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_egg_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME extran1_filled_circular_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_filled_circular_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME extran1_forcemain_hw_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_forcemain_hw_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME extran1_gothic_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_gothic_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME extran1_horseshoe_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_horseshoe_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME extran1_irregular_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_irregular_mod.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME extran1_parabolic_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_parabolic_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # # Passed 20230506
# add_test(NAME extran1_rect_closed_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_rect_closed_mod.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME extran1_rect_open_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_rect_open_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # # Passed 20230506
# add_test(NAME extran1_semi_circular_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_semi_circular_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME extran1_semi_elliptical_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_semi_elliptical_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME extran1_vertical_ellipse_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_vertical_ellipse_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

###############################################  LEVEL 2 PASSED ###############################

# # FAILED 20230424
# add_test(NAME extran3_bottom_orifice_pid
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_bottom_orifice_pid.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME extran3_side_open_close_orifice
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_side_open_close_orifice.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME extran3_side_open_close_rect_offset_orifice
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_side_open_close_rect_offset_orifice.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME extran3_side_open_close_rect_orifice
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_side_open_close_rect_orifice.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME extran3_side_orifice
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_side_orifice.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # # FAILED 20230424
# add_test(NAME extran3_side_rect_orifice
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran3_side_rect_orifice.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # # FAILED 20230424
# add_test(NAME extran4_roadway
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran4_roadway.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # # FAILED 20230424
# add_test(NAME extran4_sideflow
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran4_sideflow.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # # FAILED 20230424
# add_test(NAME extran4_trapezoidal
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran4_trapezoidal.inp -s ./ctest/internal_testing/level02_test/settings_10.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # # FAILED 20230424
# add_test(NAME extran4_vnotch
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran4_vnotch.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME extran4
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran4.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME extran6_storage_functional_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_storage_functional_final.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# add_test(NAME extran6_pump_type1_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_pump_type1_mod.inp -s ./ctest/internal_testing/level02_test/settings_10.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME extran6_pump_type2_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_pump_type2_final.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME extran6_pump_type3_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_pump_type3_final.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
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

# # FAILS custom geometry NOT WORKING AS OF 20221201
# add_test(NAME extran1_custom_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_custom_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# ## FAILED power geometry NOT WORKING AS OF 20222101
# # add_test(NAME extran1_power_mod
# # COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_power_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# # WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# # )

# # FAILED 20230506
# add_test(NAME extran1_forcemain_dw_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_forcemain_dw_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230506
# add_test(NAME extran1_horizontal_ellipse_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_horizontal_ellipse_final.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230506
# add_test(NAME extran1_mod_baskethandle_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_mod_baskethandle_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230506
# add_test(NAME extran1_rect_open_test_A05
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/extran1_rect_open_test_A05.inp -s ./ctest/internal_testing/brh/settings_extran1_test_A.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230506
# add_test(NAME extran1_rect_round_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_rect_round_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230506
# add_test(NAME extran1_rect_triangular_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_rect_triangular_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230506
# add_test(NAME force_main_example_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/force_main_example_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# FAILED 20230506
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



# # FAILED 20230424
# add_test(NAME extran4_no_surcharge
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran4_no_surcharge.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )



# # FAILED 20230424
# add_test(NAME extran6_pump_type1_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_pump_type1_mod.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME extran6_control_pump_type2_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran6_control_pump_type2_final.inp -s ./ctest/internal_testing/level02_test/settings_default.json -o ./ctest/temp_output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

###################### FAILURES  -- LEVEL 3

# # FAILED 20230424
# add_test(NAME catchment_as_outfall
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/catchment_as_outfall.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME Example1_mod01
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/Example1_mod01.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME Example1
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/Example1.inp -s ./ctest/internal_testing/level03_test/settings_10.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME Example3
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/Example3.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME gate_control_2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/gate_control_2.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAIL 20230104 -- requires more than 5 upstream connections
# add_test(NAME Greenville_EPA
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/Greenville_EPA.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME user1
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user1.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME user2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user2.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME user3
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user3.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME user4
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user4.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME user5_F08
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_F08.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # FAILED 20230424
# add_test(NAME user5_C
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_C.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME user5_C02
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_C02.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME user5_E
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_E.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424
# add_test(NAME user5_F
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_F.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424 -- controlled oscillations, very similar to EPA SWMM in .13_NAT, .21NAT, 0.06_NAT, 700pipe, 702pipe, 717, 740pipe, 801,
# add_test(NAME user5_G
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_G.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424 -- controlled oscillations, very similar to EPA SWMM -- see pipes in user5_G case
# add_test(NAME user5_H
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_H.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424 -- controlled oscillations, very similar to EPA SWMM -- limited pipe output
# add_test(NAME user5_I
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_I.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230424, but strange results - plus/minus oscillations of small flowrate
# add_test(NAME user5_JN
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_JN.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230424 (zero flow case)
# add_test(NAME user5_KN
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_KN.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230424 -- zero flows, looks fine, EPA SWMM oscillating
# add_test(NAME user5_L
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_L.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# #  FAILED 20230424 (NO OUTPUT)
# add_test(NAME user5_LN
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_LN.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230424 -- zero flows, EPA SWMM is wildly oscillating
# add_test(NAME user5_M
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_M.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230424 -- zero flows (1e-12)
# add_test(NAME user5_MN
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_MN.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230424 -- zero flows
# add_test(NAME user5_O
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_O.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# FAILED 20230503
# add_test(NAME user5_A01
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh05/user5_A01.inp -s ./ctest/internal_testing/brh05/settings_dx10.json -o ./ctest/outbrh5/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # oscillations due to inconsistent initial conditions -- not a good test case
# add_test(NAME user5_K
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh05/user5_K.inp -s ./ctest/internal_testing/brh05/settings_dx10.json -o ./ctest/outbrh5/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230503 zero flows
# add_test(NAME user5_KN
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh05/user5_KN.inp -s ./ctest/internal_testing/brh05/settings_dx10.json -o ./ctest/outbrh5/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230503 zero flows
# add_test(NAME user5_KNa6
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh05/user5_KNa6.inp -s ./ctest/internal_testing/brh05/settings_dx10.json -o ./ctest/outbrh5/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230503 zero flows
# add_test(NAME user5_KNa6alt
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh05/user5_KNa6alt.inp -s ./ctest/internal_testing/brh05/settings_dx10.json -o ./ctest/outbrh5/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Pass?  EVERYTHING IS OSCILLATING
# add_test(NAME user5_L
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh05/user5_L.inp -s ./ctest/internal_testing/brh05/settings_dx10.json -o ./ctest/outbrh5/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230503 zero flows
# add_test(NAME user5_LN
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh05/user5_LN.inp -s ./ctest/internal_testing/brh05/settings_dx10.json -o ./ctest/outbrh5/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # Passed -- EPA SWMM has strong oscillations
# add_test(NAME user5_M
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh05/user5_M.inp -s ./ctest/internal_testing/brh05/settings_dx10.json -o ./ctest/outbrh5/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230503 zero flows
# add_test(NAME user5_MN
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh05/user5_MN.inp -s ./ctest/internal_testing/brh05/settings_dx10.json -o ./ctest/outbrh5/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230503 zero flows
# add_test(NAME user5_N
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh05/user5_N.inp -s ./ctest/internal_testing/brh05/settings_dx10.json -o ./ctest/outbrh5/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230503 zero flows
# add_test(NAME user5_O
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh05/user5_O.inp -s ./ctest/internal_testing/brh05/settings_dx10.json -o ./ctest/outbrh5/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230503
# add_test(NAME user5_reduced01
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh05/user5_reduced01.inp -s ./ctest/internal_testing/brh05/settings_dx10.json -o ./ctest/outbrh5/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # Passed 20230503  zero flow
# add_test(NAME user5_reduced02
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh05/user5_reduced02.inp -s ./ctest/internal_testing/brh05/settings_dx50.json -o ./ctest/outbrh5/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )






# # # Passed 20230506
# add_test(NAME T007_RO_Free_dx0010_NJ2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T007_RO_Free_dx0010_NJ2.inp -s ./ctest/internal_testing/brh/settings_10_faceNJ2.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T007_RO_Free_dx0100_NJ2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T007_RO_Free_dx0100_NJ2.inp -s ./ctest/internal_testing/brh/settings_100_faceNJ2.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T007_RO_Free_dx0333_NJ2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T007_RO_Free_dx0333_NJ2.inp -s ./ctest/internal_testing/brh/settings_333_faceNJ2.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T007_RO_Free_dx0010_NJ2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007_RO_Free_dx0010_NJ2.inp -s ./ctest/internal_testing/brh03/settings_10_NJ2.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T007_RO_Free_dx0010
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007_RO_Free_dx0010.inp -s ./ctest/internal_testing/brh03/settings_10.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T007_RO_Free_dx0010_T20
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007_RO_Free_dx0010_T20.inp -s ./ctest/internal_testing/brh03/settings_10_T20.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T007_RO_Free_dx0100_NJ2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007_RO_Free_dx0100_NJ2.inp -s ./ctest/internal_testing/brh03/settings_100_NJ2.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T007_RO_Free_dx0100
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007_RO_Free_dx0100.inp -s ./ctest/internal_testing/brh03/settings_100.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T007_RO_Free_dx0100_T20
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007_RO_Free_dx0100_T20.inp -s ./ctest/internal_testing/brh03/settings_100_T20.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T007_RO_Free_dx0250_NJ2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007_RO_Free_dx0250_NJ2.inp -s ./ctest/internal_testing/brh03/settings_250_NJ2.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T007_RO_Free_dx0250
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007_RO_Free_dx0250.inp -s ./ctest/internal_testing/brh03/settings_250.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T007_RO_Free_dx0250_T20
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007_RO_Free_dx0250_T20.inp -s ./ctest/internal_testing/brh03/settings_250_T20.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T007_RO_Free_dx0333
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007_RO_Free_dx0333.inp -s ./ctest/internal_testing/brh03/settings_333.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T007_RO_Free_dx0333_T20
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007_RO_Free_dx0333_T20.inp -s ./ctest/internal_testing/brh03/settings_333_T20.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T007a_RO_Free_dx0010_NJ2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007a_RO_Free_dx0010_NJ2.inp -s ./ctest/internal_testing/brh03/settings_10_NJ2.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T007a_RO_Free_dx0010
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007a_RO_Free_dx0010.inp -s ./ctest/internal_testing/brh03/settings_10.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T007a_RO_Free_dx0010_T20
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007a_RO_Free_dx0010_T20.inp -s ./ctest/internal_testing/brh03/settings_10_T20.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T007a_RO_Free_dx0100
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007a_RO_Free_dx0100.inp -s ./ctest/internal_testing/brh03/settings_100.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T007a_RO_Free_dx0100_T20
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007a_RO_Free_dx0100_T20.inp -s ./ctest/internal_testing/brh03/settings_100_T20.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T007a_RO_Free_dx0250_NJ2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007a_RO_Free_dx0250_NJ2.inp -s ./ctest/internal_testing/brh03/settings_250_NJ2.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T007a_RO_Free_dx0250
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007a_RO_Free_dx0250.inp -s ./ctest/internal_testing/brh03/settings_250.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230506
# add_test(NAME T007c_RO_Free_dx0333
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007c_RO_Free_dx0333.inp -s ./ctest/internal_testing/brh03/settings_333.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T007a_RO_Free_dx0250_T20
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007a_RO_Free_dx0250_T20.inp -s ./ctest/internal_testing/brh03/settings_250_T20.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# #  Passed 20230506
# add_test(NAME T007a_RO_Free_dx0333
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007a_RO_Free_dx0333.inp -s ./ctest/internal_testing/brh03/settings_333.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T007a_RO_Free_dx0333_T20
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007a_RO_Free_dx0333_T20.inp -s ./ctest/internal_testing/brh03/settings_333_T20.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T007_RO_Free_dx0010_NJ2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007_RO_Free_dx0010_NJ2.inp -s ./ctest/internal_testing/brh03/settings_10_NJ2.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T007_RO_Free_dx0100
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007_RO_Free_dx0100.inp -s ./ctest/internal_testing/brh03/settings_100.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T007_RO_Free_dx0100_NJ2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007_RO_Free_dx0100_NJ2.inp -s ./ctest/internal_testing/brh03/settings_100_NJ2.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T007_RO_Free_dx0250
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007_RO_Free_dx0250.inp -s ./ctest/internal_testing/brh03/settings_250.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T007_RO_Free_dx0250_T20
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007_RO_Free_dx0250_T20.inp -s ./ctest/internal_testing/brh03/settings_250_T20.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# Passed 20230506
# add_test(NAME T007_RO_Free_dx0333
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007_RO_Free_dx0333.inp -s ./ctest/internal_testing/brh03/settings_333.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# Passed 20230506
# add_test(NAME T007_RO_Free_dx0333_NJ2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T007_RO_Free_dx0333_NJ2.inp -s ./ctest/internal_testing/brh03/settings_333_NJ2.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T008_RO_Free_dx0010
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T008_RO_Free_dx0010.inp -s ./ctest/internal_testing/brh03/settings_10.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T008_RO_Free_dx0010_NJ2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T008_RO_Free_dx0010_NJ2.inp -s ./ctest/internal_testing/brh03/settings_10_NJ2.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T008_RO_Free_dx0100_NJ2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T008_RO_Free_dx0100_NJ2.inp -s ./ctest/internal_testing/brh03/settings_100_NJ2.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T008_RO_Free_dx0250
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T008_RO_Free_dx0250.inp -s ./ctest/internal_testing/brh03/settings_250.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# Passed 20230506
# add_test(NAME T008_RO_Free_dx0100
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T008_RO_Free_dx0100.inp -s ./ctest/internal_testing/brh03/settings_100.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230506
# add_test(NAME T008_RO_Free_dx0333_NJ2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T008_RO_Free_dx0333_NJ2.inp -s ./ctest/internal_testing/brh03/settings_333_NJ2.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# passed 20230506
# add_test(NAME T008_RO_Free_dx0333
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T008_RO_Free_dx0333.inp -s ./ctest/internal_testing/brh03/settings_333.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


#  Passed 20230506
# add_test(NAME T009Wr_RO_Free_dx0010_NJ2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T009Wr_RO_Free_dx0010_NJ2.inp -s ./ctest/internal_testing/brh03/settings_10_NJ2.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

#  Passed 20230506
# add_test(NAME T009Wr_RO_Free_dx0010
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T009Wr_RO_Free_dx0010.inp -s ./ctest/internal_testing/brh03/settings_10.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# Passed 20230506
# add_test(NAME T009Wr_RO_Free_dx0100
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T009Wr_RO_Free_dx0100.inp -s ./ctest/internal_testing/brh03/settings_100.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# Passed 20230506
# add_test(NAME T009Wr_RO_Free_dx0100_NJ2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh03/T009Wr_RO_Free_dx0100_NJ2.inp -s ./ctest/internal_testing/brh03/settings_100_NJ2.json -o ./ctest/outbrh3/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME T010b_CC_Dmax0_sur100_dx0010
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh05/T010b_CC_Dmax0_sur100_dx0010.inp -s ./ctest/internal_testing/brh05/settings_dx10.json -o ./ctest/outbrh5/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME T010b_CC_Dmax0_sur100_dx0100
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh05/T010b_CC_Dmax0_sur100_dx0100.inp -s ./ctest/internal_testing/brh05/settings_dx100.json -o ./ctest/outbrh5/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME T010b_CC_Dmax0_sur100_dx0333
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh05/T010b_CC_Dmax0_sur100_dx0333.inp -s ./ctest/internal_testing/brh05/settings_dx333.json -o ./ctest/outbrh5/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# add_test(NAME T012_CC_Dmax0_sur1a_dx0100
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh05/T012_CC_Dmax0_sur1a_dx0100.inp -s ./ctest/internal_testing/brh05/settings_dx100.json -o ./ctest/outbrh5/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# FAILED 20230510 -- OK functioning, but timing of changes is different
# add_test(NAME T013_pumptype1.inp
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh06/T013_pumptype1.inp -s ./ctest/internal_testing/brh06/settings_dx10.json -o ./ctest/outbrh6/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230510
# add_test(NAME T013_pumptype2.inp
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh06/T013_pumptype2.inp -s ./ctest/internal_testing/brh06/settings_dx10.json -o ./ctest/outbrh6/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230510
# add_test(NAME T013_pumptype3.inp
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh06/T013_pumptype3.inp -s ./ctest/internal_testing/brh06/settings_dx10.json -o ./ctest/outbrh6/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# FAILED 20230510 -- functionality looks OK, but slightly differen behavior
# add_test(NAME T013_pumptype4.inp
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh06/T013_pumptype4.inp -s ./ctest/internal_testing/brh06/settings_dx10.json -o ./ctest/outbrh6/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # OK Failure 20230511
# add_test(NAME user5_reduced04
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh06/user5_reduced04.inp -s ./ctest/internal_testing/brh06/settings_dx100.json -o ./ctest/outbrh6/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # OK Failure 20230511
# add_test(NAME user5_reduced05
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh06/user5_reduced05.inp -s ./ctest/internal_testing/brh06/settings_dx150.json -o ./ctest/outbrh6/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # OK Failure 20230511
# add_test(NAME user5_reduced06
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh06/user5_reduced06.inp -s ./ctest/internal_testing/brh06/settings_dx150.json -o ./ctest/outbrh6/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# Passed 20230512 --- oscillations and differences are explainable based on profile views.
#  For example, 740pipe has upstream supercritical with jump to downstream subcritical and surcharged
#  which leads to waves and is not well-represented by simply averaging the flowrates and
#  depths. Also, pipe .21NAT is a deep pipe where 1e-3 ft (0.01 in) free surface oscillations lead to large  (+/-1 cfs)
#  flowrate oscillations
# add_test(NAME user5_reduced01a
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh06/user5_reduced01a.inp -s ./ctest/internal_testing/brh06/settings_dx5.json -o ./ctest/outbrh6/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230512
# add_test(NAME user5_reduced01
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh06/user5_reduced01.inp -s ./ctest/internal_testing/brh06/settings_dx5.json -o ./ctest/outbrh6/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# add_test(NAME user5
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh06/user5.inp -s ./ctest/internal_testing/brh06/settings_dx5.json -o ./ctest/outbrh6/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME user5_flowpath102_only
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh06/user5_flowpath102_only.inp -s ./ctest/internal_testing/brh06/settings_dx5.json -o ./ctest/outbrh6/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

#=================================================================================================
#  V1 TEST SET 
#=================================================================================================
# add_test(NAME user5a
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/V1test_prelim/user5a.inp -s ./ctest/V1test_prelim/user5a_settings.json -o ./ctest/outV1prelim/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME T007_RO_Free_dx100_forceJunction
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/V1test_prelim/T007_RO_Free_dx100_forceJunction.inp -s ./ctest/V1test_prelim/T007_RO_Free_dx100_forceJunction_settings.json -o ./ctest/outV1prelim/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME T007_RO_Free_dx100
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/V1test_prelim/T007_RO_Free_dx100.inp -s ./ctest/V1test_prelim/T007_RO_Free_dx100_settings.json -o ./ctest/outV1prelim/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME T007a_RO_Free_dx100_forceJunction
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/V1test_prelim/T007a_RO_Free_dx100_forceJunction.inp -s ./ctest/V1test_prelim/T007a_RO_Free_dx100_forceJunction_settings.json -o ./ctest/outV1prelim/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME T007a_RO_Free_dx100
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/V1test_prelim/T007a_RO_Free_dx100.inp -s ./ctest/V1test_prelim/T007a_RO_Free_dx100_settings.json -o ./ctest/outV1prelim/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME T007b_RO_Free_dx100_forceJunction
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/V1test_prelim/T007b_RO_Free_dx100_forceJunction.inp -s ./ctest/V1test_prelim/T007b_RO_Free_dx100_forceJunction_settings.json -o ./ctest/outV1prelim/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME T007b_RO_Free_dx100
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/V1test_prelim/T007b_RO_Free_dx100.inp -s ./ctest/V1test_prelim/T007b_RO_Free_dx100_settings.json -o ./ctest/outV1prelim/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME T007c_RO_Free_dx100_forceJunction
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/V1test_prelim/T007c_RO_Free_dx100_forceJunction.inp -s ./ctest/V1test_prelim/T007c_RO_Free_dx100_forceJunction_settings.json -o ./ctest/outV1prelim/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# add_test(NAME T007c_RO_Free_dx100
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/V1test_prelim/T007c_RO_Free_dx100.inp -s ./ctest/V1test_prelim/T007c_RO_Free_dx100_settings.json -o ./ctest/outV1prelim/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )