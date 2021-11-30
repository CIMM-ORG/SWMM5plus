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
        enumerator :: lChannel              !% channel link
        enumerator :: lPipe                 !% pipe link
        enumerator :: lWeir                 !% weir link
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
        enumerator :: lTypeIdealPump        !% ideal pump link
        enumerator :: lOutlet               !% outlet link
        enumerator :: lFunctionalDepth      !% outlet having functional depth relationship  
        enumerator :: lFunctionalHead       !% outlet having functional head relationship 
        enumerator :: lTabularDepth         !% outlet having curve relationship to depth
        enumerator :: lTabularHead          !% outlet having curve relationship to head
        !% open channel cross-sectional geometry types
        enumerator :: lRectangular           !% rectangular open channel
        enumerator :: lTrapezoidal           !% trapezoidal open channel
        enumerator :: lTriangular            !% triangular open channel
        enumerator :: lParabolic             !% parabolic open channel
        enumerator :: lPower_function        !% power function open channel
        enumerator :: lRect_triang           !% rectangular-triangular open channel
        enumerator :: lRect_round            !% rectangular-round open channel
        enumerator :: lMod_basket            !% modified basket handle open channel
        enumerator :: lIrregular             !% irregular open channel
        !% closed conduit cross-sectional geometry types
        enumerator :: lCircular              !% circular closed conduit
        enumerator :: lFilled_circular       !% filled circular closed conduit
        enumerator :: lRectangular_closed    !% rectangular closed conduit
        enumerator :: lHoriz_ellipse         !% horizontal ellipse closed conduit
        enumerator :: lVert_ellipse          !% vertical ellipse closed conduit
        enumerator :: lArch                  !% arch closed conduit
        enumerator :: lEggshaped             !% eggshaped closed conduit
        enumerator :: lHorseshoe             !% horseshoe closed conduit
        enumerator :: lGothic                !% gothic closed conduit
        enumerator :: lCatenary              !% catenary closed conduit
        enumerator :: lSemi_elliptical       !% semi-elliptical closed conduit
        enumerator :: lBasket_handle         !% basket handle closed conduit
        enumerator :: lSemi_circular         !% semi-circular closed conduit
        enumerator :: lCustom                !% custom closed conduit
        enumerator :: lForce_main            !% force main closed conduit
        !% different link roughness types
        enumerator :: lManningsN            !% ManningsN roughness type
        enumerator :: lCD                   !% drag coefficient roughness type
        !% different node types (HACK: probably could be consolidated to element types)
        enumerator :: nJ2                   !% junction node with 2 links
        enumerator :: nJm                   !% junction node with multiple links
        enumerator :: nStorage              !% storage node
        enumerator :: nBCdn                 !% downstream BC node
        enumerator :: nBCup                 !% upstream BC node
        enumerator :: nBClat                !% lateral BC node
        !% SWMM5+ elements types
        enumerator :: CC                    !% conduit or channel element
        enumerator :: weir                  !% weir element
        enumerator :: orifice               !% orifice element
        enumerator :: pump                  !% pump element
        enumerator :: outlet                !% outlet element
        enumerator :: JM                    !% junction main element
        enumerator :: JB                    !% junction branch element
        enumerator :: ArtificalStorage      !% junction main storage is artificially created
        enumerator :: FunctionalStorage     !% junction main storage is cauculated from a user provided function
        enumerator :: TabularStorage        !% junction main storage is cauculated from a user provided table
        enumerator :: storage               !% storage element
        enumerator :: manhole               !% manhole elemen (HACK: not sure if we need this)
        enumerator :: dummy                 !% dummy element type
        !% SWMM5+ CC geometry types
        !% open channel cross-sectional geometry types
        enumerator :: rectangular           !% rectangular open channel
        enumerator :: trapezoidal           !% trapezoidal open channel
        enumerator :: triangular            !% triangular open channel
        enumerator :: parabolic             !% parabolic open channel
        enumerator :: power_function        !% power function open channel
        enumerator :: rect_triang           !% rectangular-triangular open channel
        enumerator :: rect_round            !% rectangular-round open channel
        enumerator :: mod_basket            !% modified basket handle open channel
        enumerator :: irregular             !% irregular open channel
        !% closed conduit cross-sectional geometry types
        enumerator :: circular              !% circular closed conduit
        enumerator :: filled_circular       !% filled circular closed conduit
        enumerator :: rectangular_closed    !% rectangular closed conduit
        enumerator :: horiz_ellipse         !% horizontal ellipse closed conduit
        enumerator :: vert_ellipse          !% vertical ellipse closed conduit
        enumerator :: arch                  !% arch closed conduit
        enumerator :: eggshaped             !% eggshaped closed conduit
        enumerator :: horseshoe             !% horseshoe closed conduit
        enumerator :: gothic                !% gothic closed conduit
        enumerator :: catenary              !% catenary closed conduit
        enumerator :: semi_elliptical       !% semi-elliptical closed conduit
        enumerator :: basket_handle         !% basket handle closed conduit
        enumerator :: semi_circular         !% semi-circular closed conduit
        enumerator :: custom                !% custom closed conduit
        enumerator :: force_main            !% force main closed conduit
        !% SWMM5+ CC roughness type
        enumerator :: ManningsN             !% ID for mannings n for roughness_type
        enumerator :: CD                    !% ID for using drag coefficient for roughness_type
        !% SWMM5+ element types based on time marching
        enumerator :: diagnostic            !% diagnostic element
        enumerator :: time_march            !% indicates a time marched
        enumerator :: notused               !% where no Q method is used (junctions)
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
        enumerator :: functional_depth      !% outlet having functional depth relationship  
        enumerator :: functional_head       !% outlet having functional head relationship 
        enumerator :: tabular_depth         !% outlet having curve relationship to depth
        enumerator :: tabular_head          !% outlet having curve relationship to head
        !% BC types
        enumerator :: BCFlow
        enumerator :: BCHead
        !% BC category
        enumerator :: BCdn                  !% downstream BC
        enumerator :: BCup                  !% upstream BC
        enumerator :: BClat                 !% lateral BC
        !% BC subcategory (Q - flowrate)
        enumerator :: BCQ_fixed
        enumerator :: BCQ_tseries
        !% BC subcategory (H - head)
        enumerator :: BCH_free             !% minimum of critical flow depth an normal flow depth
        enumerator :: BCH_normal           !% outfall stage based on normal flow depth in connecting conduit
        enumerator :: BCH_fixed            !% outfall stage set to a fixed value
        enumerator :: BCH_tidal            !% outfall stage given by a table of tide elevation versus time of day
        enumerator :: BCH_tseries          !% outfall stage supplied from a time series of elevations
        !% face jump types
        enumerator :: jump_none             !% type of hydraulic jump
        enumerator :: jump_from_upstream    !% type of hydraulic jump
        enumerator :: jump_from_downstream  !% type of hydraulic jump
        !% type of momentum source
        enumerator :: T00                   !% type of momentum source
        enumerator :: T10                   !% type of momentum source
        enumerator :: T20                   !% type of momentum source
        !% type of face interpolation for downstream JB
        enumerator :: static
        enumerator :: dynamic
        !% MISC keys
        enumerator :: doesnotexist          !% used in various places to indicate something doesn't exist
        enumerator :: vshape                !% type of adjustment
        enumerator :: vshape_surcharge_only !% type of adjustment
        enumerator :: FroudeNumber          !% data types for limiter BC approachs
        ! data types for Partitioing Algorithm type (setting%Partitioning%PartitioningMethod)
        enumerator :: Default
        enumerator :: BQuick
        enumerator :: Random
        enumerator :: BLink
        enumerator :: StaticSlot
        enumerator :: VariableSlot
        !% keys for report time processing
        enumerator :: InSeconds
        enumerator :: InMinutes
        enumerator :: InHours  
        enumerator :: InDays    
        !% keys for output link processing
        enumerator :: AverageElements
        enumerator :: SumElements
        enumerator :: MaximumValue
        enumerator :: SingleValue
        !% keys for output FeatureType
        enumerator :: LinkOut
        enumerator :: NodeElemOut
        enumerator :: NodeFaceOut
        !% keys for curve types
        enumerator :: StorageCurve             !% surf. area v. depth for storage node
        enumerator :: DiversionCurve           !% diverted flow v. inflow for divider node
        enumerator :: TidalCurve               !% water elev. v. hour of day for outfall
        enumerator :: RatingCurve              !% flow rate v. head for outlet link
        enumerator :: ControlCurve             !% control setting v. controller variable
        enumerator :: ShapeCurve               !% width v. depth for custom x-section
        enumerator :: WeirCurve                !% discharge coeff. v. head for weir
        enumerator :: Pump1Curve               !% flow v. wet well volume for pump
        enumerator :: Pump2Curve               !% flow v. depth for pump (discrete)
        enumerator :: Pump3Curve               !% flow v. head for pump (continuous)
        enumerator :: Pump4Curve               !% flow v. depth for pump (continuous)
        !% keys for minimum element size adjustment
        enumerator :: ElemLengthAdjust         !% only adjust the element length and set it to user defined min
        enumerator :: RawElemLength            !% keep the raw data and do not adjust
    end enum

    !% rm 20210610 brh because of issues with ALLtm key
    ! ! data types for solver (setting%Solver%SolverSelect)
    ! enum, bind(c)
    !     enumerator :: ETM = 1
    !     enumerator :: ETM_AC
    !     enumerator :: AC
    ! end enum



end module define_keys
