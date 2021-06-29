! module define_keys
!
! Provides relationship between integers and keys used for different
! data types.
!
! For example, the elem2I(:,:) array has column ei_elem_type that provides
! the index to the element type for each element. The possible values are
! provided below as eChannel, ePipe, etc.
!
!==========================================================================
module define_keys

    implicit none

    !% BRH 20210608 New data keyes needed for identifiying time loops
    enum, bind(c)
        enumerator :: hydrology = 1         !% indicates hydrology loop
        enumerator :: hydraulics            !% indicates hydraulics loop
        enumerator :: ALLtm                 !% indictes all time marching types
        enumerator :: ETM                   !% Explicit time march
        enumerator :: ETM_AC                !% Explicit time march and AC
        enumerator :: AC                    !% AC
        enumerator :: CCJM                  !% channel, conduit or junction main
        enumerator :: ALL                   !% all elements
        !% different link and their geometry types (HACK: probably could be consolidated to element types)
        enumerator :: lchannel              !% channel link
        enumerator :: lpipe                 !% pipe link
        enumerator :: lweir                 !% weir link
        enumerator :: lTransverseWeir       !% transverse weir link
        enumerator :: lSideFlowWeir         !% sideflow weir link
        enumerator :: lRoadWayWeir          !% roadway weir link
        enumerator :: lVnotchWeir           !% vnotch weir link
        enumerator :: lTrapezoidalWeir      !% trapezoidal weir link
        enumerator :: lOrifice              !% orifice link
        enumerator :: lBottomOrifice        !% bottom orifice link
        enumerator :: lSideOrifice          !% side orifice link
        enumerator :: lPump                 !% pump link
        enumerator :: lType1Pump            !% type 1 pump link
        enumerator :: lType2Pump            !% type 2 pump link
        enumerator :: lType3Pump            !% type 3 pump link
        enumerator :: lType4Pump            !% type 4 pump link
        enumerator :: lRectangular          !% rectangular link geometry
        enumerator :: lParabolic            !% parabolic link geometry
        enumerator :: lTrapezoidal          !% trapezoidal link geometry
        enumerator :: lTriangular           !% triangular link geometry
        enumerator :: lWidthDepth           !% widthdepth link geometry
        enumerator :: lCircular             !% circular link geometry
        !% different link roughness types
        enumerator :: lManningsN            !% ManningsN roughness type
        enumerator :: lCD                   !% drag coefficient roughness type
        !% different node types (HACK: probably could be consolidated to element types)
        enumerator :: nJ2                   !% junction node with 2 links
        enumerator :: nJm                   !% junction node with multiple links
        enumerator :: nStorage              !% stroage node
        enumerator :: nBCdn                 !% downstream BC node
        enumerator :: nBCup                 !% upstream BC node
        !% SWMM5+ elements types
        enumerator :: CC                    !% conduit or channel element
        enumerator :: weir                  !% weir element
        enumerator :: orifice               !% orifice element
        enumerator :: pump                  !% pump element
        enumerator :: JM                    !% junction main element
        enumerator :: JB                    !% junction branch element
        enumerator :: storage               !% storage element
        enumerator :: manhole               !% manhole elemen (HACK: not sure if we need this)
        enumerator :: dummy                 !% dummy element type
        !% SWMM5+ CC geometry types
        enumerator :: rectangular           !% rectangular channel or conduit
        enumerator :: circular              !% circular channel or conduit
        enumerator :: triangular            !% triangular channel or conduit
        enumerator :: trapezoidal           !% trapezoidal channel or conduit
        enumerator :: parabolic             !% parabolic channel or conduit
        !% SWMM5+ CC roughness type
        enumerator :: ManningsN             !% ID for mannings n for roughness_type
        enumerator :: CD                    !% ID for using drag coefficient for roughness_type
        !% SWMM5+ element types based on time marching
        enumerator :: diagnostic            !% diagnostic element
        enumerator :: time_march            !% indicates a time marched
        !% SWMM5+ special element types
        enumerator :: transverse_weir       !% transverse weir type
        enumerator :: side_flow             !% sideflow weir type
        enumerator :: roadway_weir          !% roadway weir type
        enumerator :: vnotch_weir           !% vnotch weir type
        enumerator :: trapezoidal_weir      !% trapezoidal weir type
        enumerator :: bottom_orifice        !% bottom orifice type
        enumerator :: side_orifice          !% side orifice type
        enumerator :: type1_Pump            !% type 1 pump type
        enumerator :: type2_Pump            !% type 2 pump type
        enumerator :: type3_Pump            !% type 3 pump type
        enumerator :: type4_Pump            !% type 4 pump type
        !% BC face types
        enumerator :: BCup                  !% Up BC face
        enumerator :: BCdn                  !% Dn BC face
        enumerator :: critical              !% critical depth BC designator
        enumerator :: depth                 !% depth (subcritical) BC designator
        enumerator :: flowrate              !% flowrate BC designato
        !% face jump types
        enumerator :: jump_none             !% type of hydraulic jump
        enumerator :: jump_from_upstream    !% type of hydraulic jump
        enumerator :: jump_from_downstream  !% type of hydraulic jump
        !% type of momentum source
        enumerator :: T00                   !% type of momentum source
        enumerator :: T10                   !% type of momentum source
        enumerator :: T20                   !% type of momentum source
        !% MISC keys
        enumerator :: doesnotexist          !% used in various places to indicate something doesn't exist
        enumerator :: vshape                !% type of adjustment
        enumerator :: vshape_surcharge_only !% type of adjustment
        enumerator :: FroudeNumber          !% data types for limiter BC approach
    end enum

    ! data types for bcdata
    integer, parameter :: bc_updn_downstream = 1
    integer, parameter :: bc_updn_upstream = 0
    integer, parameter :: bc_category_elevation = 0
    integer, parameter :: bc_category_inflowrate = 1

    ! data types for Partitioing Algorithm type (setting%Partitioning%PartitioningMethod)
    enum, bind(c)
        enumerator :: Default = 1
        enumerator :: BQuick
        enumerator :: Random
        enumerator :: BLink
    end enum

    !% rm 20210610 brh because of issues with ALLtm key
    ! ! data types for solver (setting%Solver%SolverSelect)
    ! enum, bind(c)
    !     enumerator :: ETM = 1
    !     enumerator :: ETM_AC
    !     enumerator :: AC
    ! end enum

end module define_keys
