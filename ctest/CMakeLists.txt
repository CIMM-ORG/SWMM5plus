enable_testing()


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

## passed 20230403
# add_test(NAME T001_CC_Fix_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_CC_Fix_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230404
# add_test(NAME T001_CC_Free_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_CC_Free_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230404
# add_test(NAME T001_CC_Normal_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_CC_Normal_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## passed 20230404
# add_test(NAME T001_RC_Free_dx1000 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_RC_Free_dx1000.inp -s ./ctest/internal_testing/brh/settings_1000.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## passed 20230404
# add_test(NAME T001_RC_Free_dx0500 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_RC_Free_dx0500.inp -s ./ctest/internal_testing/brh/settings_500.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## passed 20230404
# add_test(NAME T001_RC_Free_dx0333 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_RC_Free_dx0333.inp -s ./ctest/internal_testing/brh/settings_333.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## passed 20230404
# add_test(NAME T001_RC_Free_dx0250 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_RC_Free_dx0250.inp -s ./ctest/internal_testing/brh/settings_250.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


## passed 20230404
# add_test(NAME T001_RC_Free_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_RC_Free_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## passed 20230404
# add_test(NAME T001_RO_Free_dx1000 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_RC_Free_dx1000.inp -s ./ctest/internal_testing/brh/settings_1000.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILS in init_IC_from_linkdata from segmentation fault
# add_test(NAME T001_RO_Free_dx0500 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_RO_Free_dx0500.inp -s ./ctest/internal_testing/brh/settings_500.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILS in init_IC_from_linkdata from segmentation fault
# add_test(NAME T001_RO_Free_dx0250 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_RO_Free_dx0250.inp -s ./ctest/internal_testing/brh/settings_250.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILS in init_IC_from_linkdata from segmentation fault
# add_test(NAME T001_RO_Free_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T001_RO_Free_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILS in init_IC_from_linkdata from segmentation fault
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

## Passed 20230407
# add_test(NAME T002_CC_Fix_dx0010_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Fix_dx0010_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_10_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230405
# add_test(NAME T002_CC_Fix_dx0100_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Fix_dx0100_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_100_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230405
# add_test(NAME T002_CC_Fix_dx0250_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Fix_dx0250_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_250_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230407 -- but oscillations in flowrate of upstream link
# add_test(NAME T002_CC_Fix_dx0500_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Fix_dx0500_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_500_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230407 -- increasing upstream conduit length removed problem in T002_CC_Fix_dx0500
# add_test(NAME T002a_CC_Fix_dx0500_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002a_CC_Fix_dx0500_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_500_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230407 -- but oscillations in flowrate of upstream link as in T002_CC_Fix_dx500_ForceStorage
# add_test(NAME T002b_CC_Fix_dx0500_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002b_CC_Fix_dx0500_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_500_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

############################################################################################
##### LESSON: MUST BE AT LEAST TWO ELEMENTS BETWEEN AN UPSTREAM INFLOW AND A JUNCTION
############################################################################################


    ###### FREE DOWNSTREAM BOUNDARY #####

## Passed 20230407
# add_test(NAME T002_CC_Free_dx0010_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Free_dx0010_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_10_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230407
# add_test(NAME T002_CC_Free_dx0100_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Free_dx0100_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_100_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230407
# add_test(NAME T002_CC_Free_dx0250_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Free_dx0250_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_250_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passesd 20230407 --- BUT OSCILLATIONS WITH UPSTREAM LINK DUE TO ONLY 1 ELEMENT BETWEEN INFLOW AND NODE
# add_test(NAME T002_CC_Free_dx0500_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Free_dx0500_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_500_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILS at stop 6698734 with error message "small Rquad "
# add_test(NAME T002_RC_Fix_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_RC_Fix_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILS in init_IC_from_linkdata from segmentation fault
# add_test(NAME T002_RO_Free_dx0500 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_RO_Free_dx0500.inp -s ./ctest/internal_testing/brh/settings_500.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILS in init_IC_from_linkdata from segmentation fault
# add_test(NAME T002_RO_Free_dx0250 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_RO_Free_dx0250.inp -s ./ctest/internal_testing/brh/settings_250.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILS in init_IC_from_linkdata from segmentation fault
# add_test(NAME T002_RO_Free_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_RO_Free_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILS at stop 239874 with error message "small bquad "
# add_test(NAME T002_RC_Free_dx0500 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_RC_Free_dx0500.inp -s ./ctest/internal_testing/brh/settings_500.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILS at stop 6698734 with error message "small Rquad "
# add_test(NAME T002_RC_Free_dx0250 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_RC_Free_dx0250.inp -s ./ctest/internal_testing/brh/settings_250.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILS at stop 6698734 with error message "small Rquad "
# add_test(NAME T002_RC_Free_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_RC_Free_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## passed 20230404 (terminal shows "skipping junction mass conservation" message)
# add_test(NAME T002_RC_Free_dx0010 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_RC_Free_dx0010.inp -s ./ctest/internal_testing/brh/settings_10.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILS at stop 239874 with error message "small bquad "
# add_test(NAME T002_CC_Fix_dx0250 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Fix_dx0250.inp -s ./ctest/internal_testing/brh/settings_250.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILS at stop 6698734 with error message "small Rquad "
# add_test(NAME T002_CC_Fix_dx0167 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Fix_dx0167.inp -s ./ctest/internal_testing/brh/settings_167.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILS at stop 6698734 with error message "small Rquad "
# add_test(NAME T002_CC_Fix_dx0125 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Fix_dx0125.inp -s ./ctest/internal_testing/brh/settings_125.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILS at stop 6698734 with error message "small Rquad "
# add_test(NAME T002_CC_Fix_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Fix_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILS at stop 6698734 with error message "small Rquad "
# add_test(NAME T002_CC_Fix_dx0010 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Fix_dx0010.inp -s ./ctest/internal_testing/brh/settings_10.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILS at stop 6698734 with error message "small Rquad "
# add_test(NAME T002_CC_Free_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Free_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILS at stop 6698734 with error message "small Rquad "
# add_test(NAME T002_CC_Free_dx0010 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T002_CC_Free_dx0010.inp -s ./ctest/internal_testing/brh/settings_10.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILS at stop 6698734 with error message "small Rquad "
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


# # Passed 20230410
# add_test(NAME T003_CC_Fix_dx0010_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_CC_Fix_dx0010_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_10_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230410
# add_test(NAME T003_CC_Fix_dx0111_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_CC_Fix_dx0111_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_111_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # Failed 20230410 -- but looks OK
# add_test(NAME T003_CC_Fix_dx0333_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_CC_Fix_dx0333_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_333_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


    ##### FREE DOWNSTREAM BOUNDARY #####


# # Passed 20230410
# add_test(NAME T003_CC_Free_dx0010_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_CC_Free_dx0010_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_10_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230410 
# add_test(NAME T003_CC_Free_dx0111_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_CC_Free_dx0111_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_111_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230410 -- due do depths but looks similar to others with dx333 -- time lage is aproblem
# add_test(NAME T003_CC_Free_dx0333_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_CC_Free_dx0333_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_333_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILED 20230410 -- time lag
# add_test(NAME T003_RC_Free_dx0333 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_RC_Free_dx0333.inp -s ./ctest/internal_testing/brh/settings_333.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## passed 20230410
# add_test(NAME T003_RC_Free_dx0167 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_RC_Free_dx0167.inp -s ./ctest/internal_testing/brh/settings_167.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## passed 20230410
# add_test(NAME T003_RC_Free_dx0111 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_RC_Free_dx0111.inp -s ./ctest/internal_testing/brh/settings_111.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## passed 20230410
#add_test(NAME T003_RC_Free_dx0048 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_RC_Free_dx0048.inp -s ./ctest/internal_testing/brh/settings_48.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# ## passed 20230410
# add_test(NAME T003_CC_Free_dx0048 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T003_CC_Free_dx0048.inp -s ./ctest/internal_testing/brh/settings_48.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # # passed 20230410
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

# Passed 20230410
# add_test(NAME T004_CC_Free_dx0010_NJ2 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T004_CC_Free_dx0010_NJ2.inp -s ./ctest/internal_testing/brh/settings_10_faceNJ2.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230410
# add_test(NAME T004_CC_Free_dx0010_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T004_CC_Free_dx0010_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_10_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230410
# add_test(NAME T004_CC_Free_dx0050_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T004_CC_Free_dx0050_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_50_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230410 -- Link C7 flow rate too slow to increase
# add_test(NAME T004_CC_Free_dx0100_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T004_CC_Free_dx0100_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_100_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230410  -- link flowrates C8, C9, C10 (slow increase)
# add_test(NAME T004a_CC_Free_dx0100_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T004a_CC_Free_dx0100_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_100_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


    ###### FIXED DOWNSTREAM BOUNDARY #####         

# # Passed 20230410
# add_test(NAME T004_CC_Fix_dx0010_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T004_CC_Fix_dx0010_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_10_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230410
# add_test(NAME T004_CC_Free_dx0010 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T004_CC_Free_dx0010.inp -s ./ctest/internal_testing/brh/settings_10.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230410
# add_test(NAME T004_CC_Free_dx0050 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T004_CC_Free_dx0050.inp -s ./ctest/internal_testing/brh/settings_50.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# FAILED 20230410 -- time lag
# add_test(NAME T004_CC_Free_dx0100 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T004_CC_Free_dx0100.inp -s ./ctest/internal_testing/brh/settings_100.json -o ./ctest/outbrh/
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

# # passed 20230410
# add_test(NAME T005_CC_Free_dx0010_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T005_CC_Free_dx0010_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_10_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230410 -
# add_test(NAME T005_CC_Free_dx0050_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T005_CC_Free_dx0050_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_50_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# ## FAILED 20230410 -- needs more elements in each link
# add_test(NAME T005_CC_Free_dx0100_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T005_CC_Free_dx0100_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_100_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230410
# add_test(NAME T005_CC_Free_dx0010_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T005_CC_Free_dx0010_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_10_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230410
# add_test(NAME T005_CC_Free_dx0050_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T005_CC_Free_dx0050_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_50_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Failed 20230410 -- time lag
# add_test(NAME T005_CC_Free_dx0100_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T005_CC_Free_dx0100_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_100_ForceStorage.json -o ./ctest/outbrh/
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


## passed 20230410
# add_test(NAME T005Wr_RO_Free_dx0010_ForceJM 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T005Wr_RO_Free_dx0010_ForceJM.inp -s ./ctest/internal_testing/brh/settings_10_ForceJM.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## passed 20230410
# add_test(NAME T005Wr_RO_Free_dx0010_weirJM 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T005Wr_RO_Free_dx0010_weirJM.inp -s ./ctest/internal_testing/brh/settings_10_weirJM.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# passed 20230410
add_test(NAME T005Wr_RO_Free_dx0010_NJ2 
COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T005Wr_RO_Free_dx0010_NJ2.inp -s ./ctest/internal_testing/brh/settings_10_faceNJ2.json -o ./ctest/outbrh/
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)


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

# # passed 20230410
# add_test(NAME T006_CC_Free_dx0010_ForceStorage 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T006_CC_Free_dx0010_ForceStorage.inp -s ./ctest/internal_testing/brh/settings_10_ForceStorage.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230410
# add_test(NAME T006_CC_Free_dx0010_NJ2 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/T006_CC_Free_dx0010_NJ2.inp -s ./ctest/internal_testing/brh/settings_10_faceNJ2.json -o ./ctest/outbrh/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )



###################################### GEOMETRY #############################

# # Passed 20230324
# add_test(NAME Circular_FixedBC 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_sub_G=circ_IN=con_OUT=fix.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # # Passed 20230324
# add_test(NAME Circular_FreeBC_01 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_circ_01.inp -s ./ctest/internal_testing/geometry_test/settings_Co_circ.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230324
# add_test(NAME Circular_FreeBC_02 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_circ_02.inp -s ./ctest/internal_testing/geometry_test/settings_Co_circ.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # P
# add_test(NAME Circular_FreeBC_03 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_circ_03.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# #
# add_test(NAME RectClosed_FreeBC_03 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_rect_closed_03.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# #
# add_test(NAME RectOpen_FreeBC_03 
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/Co_rect_open_03.inp -s ./ctest/internal_testing/geometry_test/settings_Co.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # Passed 2023034
# add_test(NAME YJ_Geom
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/YJ_Geom_Test.inp -s ./ctest/internal_testing/geometry_test/settings_YJ.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230324
# add_test(NAME YJ_Geom_mod01
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/geometry_test/YJ_Geom_Test_mod01.inp -s ./ctest/internal_testing/geometry_test/settings_YJ.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # ###################################### LEVEL 1 PASSED #############################

#  # passed 202304010
# add_test(NAME extran1_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # passed 20230410
# add_test(NAME extran1_final_mod2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_final_mod2.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # passed 20230316
# add_test(NAME extran1_arch_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_arch_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230316
# add_test(NAME extran1_baskethandle_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_baskethandle_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230316
# add_test(NAME extran1_catenary_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_catenary_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230316
# add_test(NAME extran1_circular_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_circular_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230316
# add_test(NAME extran1_egg_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_egg_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230316
# add_test(NAME extran1_filled_circular_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_filled_circular_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230316
# add_test(NAME extran1_forcemain_hw_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_forcemain_hw_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230316
# add_test(NAME extran1_gothic_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_gothic_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # passed 20230316
# add_test(NAME extran1_horseshoe_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_horseshoe_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230316
# add_test(NAME extran1_irregular_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_irregular_mod.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230316
# add_test(NAME extran1_parabolic_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_parabolic_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # # Passed 20230316
# add_test(NAME extran1_rect_closed_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_rect_closed_mod.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230316
# add_test(NAME extran1_rect_open_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_rect_open_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # Passed 20230316
# add_test(NAME extran1_semi_circular_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_semi_circular_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230316
# add_test(NAME extran1_semi_elliptical_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_semi_elliptical_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230316
# add_test(NAME extran1_vertical_ellipse_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_vertical_ellipse_final.inp -s ./ctest/internal_testing/level01_test/settings_extran1.json -o ./ctest/output/
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

# # FAILED 20230410
# add_test(NAME extran4_trapezoidal
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level02_test/extran4_trapezoidal.inp -s ./ctest/internal_testing/level02_test/settings_10.json -o ./ctest/output/
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

# # FAILED No flows 20230316
# add_test(NAME extran1_forcemain_dw_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_forcemain_dw_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # failed 20230316 -- twonodes with RSME head > 3.3 ft
# add_test(NAME extran1_horizontal_ellipse_final
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_horizontal_ellipse_final.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230316 but strange results for link 8060 needs investigating
# add_test(NAME extran1_mod_baskethandle_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_mod_baskethandle_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230316
# add_test(NAME extran1_rect_open_test_A05
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/brh/extran1_rect_open_test_A05.inp -s ./ctest/internal_testing/brh/settings_extran1_test_A.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # failed 20230316 -- flowrate links 1600 1602
# add_test(NAME extran1_rect_round_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_rect_round_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # failed 20230316 --- BLOWUP
# add_test(NAME extran1_rect_triangular_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/extran1_rect_triangular_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230316 (configuration problem)  
# add_test(NAME force_main_example_mod
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level01_test/force_main_example_mod.inp -s ./ctest/internal_testing/level01_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # # FAILED (configuration problem)  20230316
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



# FAILED 20230410 due to inflow condition.
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

# # Passed 20230410
# add_test(NAME catchment_as_outfall
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/catchment_as_outfall.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

## FAILED 20230410 -- Q looks OK, depths are off
# add_test(NAME Example1_mod01
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/Example1_mod01.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230410 -- flow rate residual
# add_test(NAME Example1
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/Example1.inp -s ./ctest/internal_testing/level03_test/settings_10.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230224 -- requires PUMPS
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

# # FAIL 20230316 -- no plots
# add_test(NAME user1
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user1.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAIL 20230224 -- crash due to diagnostic next to junction
# add_test(NAME user2
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user2.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAIL 20230224 -- requires PUMPS
# add_test(NAME user3
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user3.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAIL 20230224 -- crash on BC type? no plots
# add_test(NAME user4
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user4.inp -s ./ctest/internal_testing/level03_test/settings_default.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAIL 20230316 -- very different flowrates
# add_test(NAME user5_F08
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_F08.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # FAIL 20230316 - no plot
# add_test(NAME user5_C
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_C.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# FAIL 20230316 -- node head at 8SB3_OF, but flows look good.
# add_test(NAME user5_C02
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_C02.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # FAILS 20230316 -- very different flows
# add_test(NAME user5_E
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_E.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILS 20230316 -- no plots
# add_test(NAME user5_F
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_F.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILS 20230316 -- controlled oscillations, very similar to EPA SWMM in .13_NAT, .21NAT, 0.06_NAT, 700pipe, 702pipe, 717, 740pipe, 801,
# add_test(NAME user5_G
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_G.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILS 20230316 -- controlled oscillations, very similar to EPA SWMM -- see pipes in user5_G case
# add_test(NAME user5_H
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_H.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILS 20230316 -- controlled oscillations, very similar to EPA SWMM -- limited pipe output
# add_test(NAME user5_I
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_I.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # Passed 20230316, but strange results - plus/minus oscillations of small flowrate
# add_test(NAME user5_JN
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_JN.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230316 -- zero flows, looks fine
# add_test(NAME user5_KN
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_KN.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230316 -- zero flows, looks fine, EPA SWMM oscillating
# add_test(NAME user5_L
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_L.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230316 -- zero flows, looks fine
# add_test(NAME user5_LN
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_LN.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230316 -- zero flows, EPA SWMM is wildly oscillating
# add_test(NAME user5_M
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_M.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230316 -- zero flows (1e-12)
# add_test(NAME user5_MN
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_MN.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # passed 20230316 -- zero flows
# add_test(NAME user5_O
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_O.inp -s ./ctest/internal_testing/level03_test/settings_user5.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )



# # FAILED 20230316 -- no plots
# add_test(NAME user5_reduced04
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_reduced04.inp -s ./ctest/internal_testing/level03_test/settings_user5reduced04.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


# # FAILED 20230316 -- due to head, but flowrate looks good
# add_test(NAME user5_reduced05
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_reduced05.inp -s ./ctest/internal_testing/level03_test/settings_user5reduced04.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )

# # FAILED 20230315 -- no plot -- conservation issue
# add_test(NAME user5_reduced06
# COMMAND python ${CMAKE_SOURCE_DIR}/compare.py -i ./ctest/internal_testing/level03_test/user5_reduced06.inp -s ./ctest/internal_testing/level03_test/settings_user5reduced04.json -o ./ctest/output/
# WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# )


