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

    use define_globals, only: reverseKey
    implicit none
    public

    !% Make sure that subroutine define_keys_reverse() is updated when a
    !% new key is added
    enum, bind(c)
        !rm brh20211211 enumerator :: hydrology = 1         !% indicates hydrology loop
        !rm brh20211211 enumerator :: hydraulics            !% indicates hydraulics loop
        enumerator :: ALLtm = 1                 !% indicates elements using any time marching type
        enumerator :: ETM                   !% Explicit time march elements or solver
        enumerator :: ETM_AC                !% Explicit time march and AC solver
        enumerator :: AC                    !% AC elements or solver
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
        enumerator :: lNodeDepth            !% outlet having functional\curve flow vs depth relationship  
        enumerator :: lNodeHead             !% outlet having functional\curve flow vs head relationship 
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
        !!% different link roughness types (OBSOLETE)
        !enumerator :: lManningsN            !% ManningsN roughness type
        !enumerator :: lCD                   !% drag coefficient roughness type
        !% different node types 
        enumerator :: nJ1                   !% a node without an inflow connecting to only 1 link    
        enumerator :: nJ2                   !% junction node with 2 links
        enumerator :: nJm                   !% junction node with multiple links
        enumerator :: nStorage              !% storage node
        enumerator :: nBCdn                 !% downstream BC node
        enumerator :: nBCup                 !% upstream BC node
        !enumerator :: nBClat                !% lateral BC node
        !% SWMM5+ elements types
        enumerator :: CC                    !% conduit or channel element
        enumerator :: weir                  !% weir element
        enumerator :: orifice               !% orifice element
        enumerator :: pump                  !% pump element
        enumerator :: outlet                !% outlet element
        enumerator :: JM                    !% junction main element
        enumerator :: JB                    !% junction branch element
        enumerator :: ImpliedStorage        !% junction main storage is artificially created
        enumerator :: FunctionalStorage     !% junction main storage is cauculated from a user provided function
        enumerator :: TabularStorage        !% junction main storage is cauculated from a user provided table
        !enumerator :: storage               !% storage element NOT USED AS OF 20220626
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
        !enumerator :: ManningsN             !% ID for mannings n for roughness_type
        !enumerator :: CD                    !% ID for using drag coefficient for roughness_type
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
        enumerator :: type_IdealPump        !% type ideal pump
        enumerator :: func_depth_outlet     !% outlet having functional depth relationship  
        enumerator :: func_head_outlet      !% outlet having functional head relationship 
        enumerator :: tabl_depth_outlet     !% outlet having curve relationship to depth
        enumerator :: tabl_head_outlet      !% outlet having curve relationship to head
        !% BC types
        enumerator :: BCFlow
        enumerator :: BCHead
        !% BC category
        enumerator :: BCdn                  !% downstream BC
        enumerator :: BCup                  !% upstream BC
        enumerator :: BClat                 !% lateral BC
        enumerator :: BCnone                !% end face without a BC
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
        enumerator :: T10s2                 !% type of momentum source
        enumerator :: TA1                   !% type of momentum source
        enumerator :: TA2                   !% type of momentum source
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
        enumerator :: DynamicSlot
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
        !% data types used for ZeroValues
        enumerator :: DepthValue
        enumerator :: VolumeValue
        enumerator :: AreaValue
        !% keys for initial depth distribution
        enumerator :: LinearlyVaryingDepth
        enumerator :: UniformDepth
        enumerator :: ExponentialDepth
        enumerator :: FixedHead
        !% Keys for uniformTable data types
        enumerator :: DepthData
        enumerator :: AreaData
        enumerator :: PerimeterData
        enumerator :: VolumeData
        enumerator :: SectionFactorData
        enumerator :: QcriticalData
        !% Keys for Force Main
        enumerator :: NotForceMain
        enumerator :: HazenWilliams
        enumerator :: DarcyWeisbach
        !% last items for bookkeeping
        enumerator :: undefinedKey
        enumerator :: keys_lastplusone
    end enum
    

 contains   
!%==========================================================================
!%
    subroutine define_keys_reverse ()
        !%------------------------------------------------------------------
        !% Description:
        !% creates the reverseKey global that provide the string 
        !% name for the keys defined in define_keys
        !%------------------------------------------------------------------
        !% Declarations
        !%------------------------------------------------------------------
        !% Preliminaries
        !%------------------------------------------------------------------
        !% allocate the global space for the reverse keys
        allocate(reverseKey(keys_lastplusone))

        !% define the reverse keys
        reverseKey(ALLtm) = 'ALLtm'
        reverseKey(ETM) = 'ETM'
        reverseKey(ETM_AC) = 'ETM_AC'
        reverseKey(AC) = 'AC'
        reverseKey(CCJM) = 'CCJM'
        reverseKey(ALL) = 'ALL'
        reverseKey(lChannel) = 'lChannel'
        reverseKey(lPipe) = 'lPipe'
        reverseKey(lWeir) = 'lWeir'
        reverseKey(lTransverseWeir) = 'lTransverseWeir'
        reverseKey(lSideFlowWeir) = 'lSideFlowWeir'
        reverseKey(lRoadWayWeir) = 'lRoadWayWeir'
        reverseKey(lVnotchWeir) = 'lVnotchWeir'
        reverseKey(lTrapezoidalWeir) = 'lTrapezoidalWeir'
        reverseKey(lOrifice) = 'lOrifice'
        reverseKey(lBottomOrifice) = 'lBottomOrifice'
        reverseKey(lSideOrifice) = 'lSideOrifice'
        reverseKey(lPump) = 'lPump'
        reverseKey(lType1Pump) = 'lType1Pump'
        reverseKey(lType2Pump) = 'lType2Pump'
        reverseKey(lType3Pump) = 'lType3Pump'
        reverseKey(lType4Pump) = 'lType4Pump'
        reverseKey(lTypeIdealPump) = 'lTypeIdealPump'
        reverseKey(lOutlet) = 'lOutlet'
        reverseKey(lNodeDepth) = 'lNodeDepth'
        reverseKey(lNodeHead) = 'lNodeHead'
        reverseKey(lRectangular) = 'lRectangular'
        reverseKey(lTrapezoidal) = 'lTrapezoidal'
        reverseKey(lTriangular) = 'lTriangular'
        reverseKey(lParabolic) = 'lParabolic'
        reverseKey(lPower_function) = 'lPower_function'
        reverseKey(lRect_triang) = 'lRect_triang '
        reverseKey(lRect_round) = 'lRect_round'
        reverseKey(lMod_basket) = 'lMod_basket'
        reverseKey(lIrregular) = 'lIrregular'
        reverseKey(lCircular) = 'lCircular'
        reverseKey(lFilled_circular) = 'lFilled_circular'
        reverseKey(lRectangular_closed) = 'lRectangular_closed'
        reverseKey(lHoriz_ellipse) = 'lHoriz_ellipse'
        reverseKey(lVert_ellipse) = 'lVert_ellipse'
        reverseKey(lArch) = 'lArch'
        reverseKey(lEggshaped) = 'lEggshaped'
        reverseKey(lHorseshoe) = 'lHorseshoe'
        reverseKey(lGothic) = 'lGothic'
        reverseKey(lCatenary) = 'lCatenary'
        reverseKey(lSemi_elliptical) = 'lSemi_elliptical'
        reverseKey(lBasket_handle) = 'lBasket_handle'
        reverseKey(lSemi_circular) = 'lSemi_circular'
        reverseKey(lCustom) = 'lCustom'
        reverseKey(lForce_main) = 'lForce_main'
        !reverseKey(lManningsN) = 'lManningsN'
        !reverseKey(lCD) = 'lCD'
        reverseKey(nJ1) = 'nJ1'
        reverseKey(nJ2) = 'nJ2'
        reverseKey(nJm) = 'nJm'
        reverseKey(nStorage) = 'nStorage'
        reverseKey(nBCdn) = 'nBCdn'
        reverseKey(nBCup) = 'nBCup'
        !reverseKey(nBClat) = 'nBClat'
        reverseKey(CC) = 'CC'
        reverseKey(weir) = 'weir'
        reverseKey(orifice) = 'orifice'
        reverseKey(pump) = 'pump'
        reverseKey(outlet) = 'outlet'
        reverseKey(JM) = 'JM'
        reverseKey(JB) = 'JB'
        reverseKey(ImpliedStorage) = 'ImpliedStorage'
        reverseKey(FunctionalStorage) = 'FunctionalStorage'
        reverseKey(TabularStorage) = 'TabularStorage'
        !reverseKey(storage) = 'storage'
        reverseKey(manhole) = 'manhole'
        reverseKey(dummy) = 'dummy'
        reverseKey(rectangular) = 'rectangular'
        reverseKey(trapezoidal) = 'trapezoidal'
        reverseKey(triangular) = 'triangular'
        reverseKey(parabolic) = 'parabolic'
        reverseKey(power_function) = 'power_function'
        reverseKey(rect_triang) = 'rect_triang'
        reverseKey(rect_round) = 'rect_round'
        reverseKey(mod_basket) = 'mod_basket'
        reverseKey(irregular) = 'irregular'
        reverseKey(circular) = 'circular'
        reverseKey(filled_circular) = 'filled_circular'
        reverseKey(rectangular_closed) = 'rectangular_closed'
        reverseKey(horiz_ellipse) = 'horiz_ellipse'
        reverseKey(vert_ellipse) = 'vert_ellipse'
        reverseKey(arch) = 'arch'
        reverseKey(eggshaped) = 'eggshaped'
        reverseKey(horseshoe) = 'horseshoe'
        reverseKey(gothic) = 'gothic'
        reverseKey(catenary) = 'catenary'
        reverseKey(semi_elliptical) = 'semi_elliptical'
        reverseKey(basket_handle) = 'basket_handle'
        reverseKey(semi_circular) = 'semi_circular'
        reverseKey(custom) = 'custom'
        reverseKey(force_main) = 'force_main'
        !reverseKey(ManningsN) = 'ManningsN'
        !reverseKey(CD) = 'CD'
        reverseKey(diagnostic) = 'diagnostic'
        reverseKey(time_march) = 'time_march'
        reverseKey(notused) = 'notused'
        reverseKey(transverse_weir) = 'transverse_weir'
        reverseKey(side_flow) = 'side_flow'
        reverseKey(roadway_weir) = 'roadway_weir'
        reverseKey(vnotch_weir) = 'vnotch_weir'
        reverseKey(trapezoidal_weir) = 'trapezoidal_weir'
        reverseKey(bottom_orifice) = 'bottom_orifice'
        reverseKey(side_orifice) = 'side_orifice'
        reverseKey(type1_Pump) = 'type1_Pump'
        reverseKey(type2_Pump) = 'type2_Pump'
        reverseKey(type3_Pump) = 'type3_Pump'
        reverseKey(type4_Pump) = 'type4_Pump'
        reverseKey(type_IdealPump) = 'type_IdealPump'
        reverseKey(func_depth_outlet) = 'func_depth_outlet'
        reverseKey(func_head_outlet) = 'func_head_outlet'
        reverseKey(tabl_depth_outlet) = 'tabl_depth_outlet'
        reverseKey(tabl_head_outlet)  = 'tabl_head_outlet'
        reverseKey(BCFlow) = 'BCFlow'
        reverseKey(BCHead) = 'BCHead'
        reverseKey(BCdn) = 'BCdn'
        reverseKey(BCup) = 'BCup'
        reverseKey(BClat) = 'BClat'
        reverseKey(BCnone) = 'BCnone'
        reverseKey(BCQ_fixed) = 'BCQ_fixed'
        reverseKey(BCQ_tseries) = 'BCQ_tseries'
        reverseKey(BCH_free) = 'BCH_free'
        reverseKey(BCH_normal) = 'BCH_normal'
        reverseKey(BCH_fixed) = 'BCH_fixed'
        reverseKey(BCH_tidal) = 'BCH_tidal'
        reverseKey(BCH_tseries) = 'BCH_tseries'
        reverseKey(jump_none) = 'jump_none'
        reverseKey(jump_from_upstream) = 'jump_from_upstream'
        reverseKey(jump_from_downstream) = 'jump_from_downstream'
        reverseKey(T00) = 'T00'
        reverseKey(T10) = 'T10'
        reverseKey(T20) = 'T20'
        reverseKey(T10s2) = 'T10s2'
        reverseKey(TA1) = 'TA1'
        reverseKey(TA2) = 'TA2'
        reverseKey(static) = 'static'
        reverseKey(dynamic) = 'dynamic'
        reverseKey(doesnotexist) = 'doesnotexist'
        reverseKey(vshape) = 'vshape'
        reverseKey(vshape_surcharge_only) = 'vshape_surcharge_only'
        reverseKey(FroudeNumber) = 'FroudeNumber'
        reverseKey(Default) = 'Default'
        reverseKey(BQuick) = 'BQuick'
        reverseKey(Random) = 'Random'
        reverseKey(BLink) = 'BLink'
        reverseKey(StaticSlot) = 'StaticSlot'
        reverseKey(DynamicSlot) = 'DynamicSlot'
        reverseKey(InSeconds) = 'InSeconds'
        reverseKey(InMinutes) = 'InMinutes'
        reverseKey(InHours) = 'InHours'
        reverseKey(InDays) = 'InDays'
        reverseKey(AverageElements) = 'AverageElements'
        reverseKey(SumElements) = 'SumElements'
        reverseKey(MaximumValue) = 'MaximumValue'
        reverseKey(SingleValue) = 'SingleValue'
        reverseKey(LinkOut) = 'LinkOut'
        reverseKey(NodeElemOut) = 'NodeElemOut'
        reverseKey(NodeFaceOut) = 'NodeFaceOut'
        reverseKey(StorageCurve) = 'StorageCurve'
        reverseKey(DiversionCurve) = 'DiversionCurve'
        reverseKey(TidalCurve) = 'TidalCurve'
        reverseKey(RatingCurve) = 'RatingCurve'
        reverseKey(ControlCurve) = 'ControlCurve'
        reverseKey(ShapeCurve) = 'ShapeCurve'
        reverseKey(WeirCurve) = 'WeirCurve'
        reverseKey(Pump1Curve) = 'Pump1Curve'
        reverseKey(Pump2Curve) = 'Pump2Curve'
        reverseKey(Pump3Curve) = 'Pump3Curve'
        reverseKey(Pump4Curve) = 'Pump4Curve'
        reverseKey(ElemLengthAdjust) = 'ElemLengthAdjust'
        reverseKey(RawElemLength) = 'RawElemLength'
        reverseKey(DepthValue)  = 'DepthValue'
        reverseKey(VolumeValue) = 'VolumeValue'
        reverseKey(AreaValue)   = 'AreaValue'
        reverseKey(LinearlyVaryingDepth)   = 'LinearlyVaryingDepth'
        reverseKey(UniformDepth)           = 'UniformDepth'
        reverseKey(ExponentialDepth)       = 'ExponentialDepth'
        reverseKey(FixedHead)              = 'FixedHead'
        reverseKey(DepthData) = 'DepthData'
        reverseKey(AreaData) = 'AreaData'
        reverseKey(PerimeterData) = 'PerimeterData'
        reverseKey(VolumeData) = 'VolumeData'
        reverseKey(SectionFactorData) = 'SectionFactorData'
        reverseKey(QcriticalData) = 'QcriticalData'
        reverseKey(NotForceMain) = 'NotForceMain'
        reverseKey(HazenWilliams) = 'HazenWilliams'
        reverseKey(DarcyWeisbach) = 'DarcyWeisbach'
        reverseKey(undefinedKey)     = 'undefinedKey'
        reverseKey(keys_lastplusone) = 'keys_lastplusone'

        !%------------------------------------------------------------------
        !% Closing
    end subroutine define_keys_reverse
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine define_keys_printByNumber ()
        !%------------------------------------------------------------------
        !% Description:
        !% prints all the keys from define_keys (for debugging)
        !%------------------------------------------------------------------
        !% Declarations
            integer :: ii
        !%------------------------------------------------------------------
        !% Preliminaries
        !%------------------------------------------------------------------
        do ii = 1,keys_lastplusone
            write(*,'(i4," = ",A)') ii, trim(reverseKey(ii))
        end do
        !%------------------------------------------------------------------
        !% Closing
    end subroutine define_keys_printByNumber
!%    
!%==========================================================================
!%==========================================================================
!%    
    subroutine define_keys_printByName ()
        !%------------------------------------------------------------------
        !% Description:
        !% prints all the keys from define_keys in alphabetical order (for debugging)
        !%------------------------------------------------------------------
        !% Declarations
        !%------------------------------------------------------------------
        !% Preliminaries
        !%------------------------------------------------------------------

        write(*,'(A," = ",i4)') trim(reverseKey(AC)) , AC
        write(*,'(A," = ",i4)') trim(reverseKey(ALL)) , ALL
        write(*,'(A," = ",i4)') trim(reverseKey(ALLtm)) , ALLtm
        write(*,'(A," = ",i4)') trim(reverseKey(arch)) , arch
        write(*,'(A," = ",i4)') trim(reverseKey(AreaValue)), AreaValue
        write(*,'(A," = ",i4)') trim(reverseKey(ImpliedStorage)) , ImpliedStorage
        write(*,'(A," = ",i4)') trim(reverseKey(AverageElements)) , AverageElements
        write(*,'(A," = ",i4)') trim(reverseKey(basket_handle)) , basket_handle
        write(*,'(A," = ",i4)') trim(reverseKey(BCFlow)) , BCFlow
        write(*,'(A," = ",i4)') trim(reverseKey(BCHead)) , BCHead
        write(*,'(A," = ",i4)') trim(reverseKey(BCH_fixed)) , BCH_fixed
        write(*,'(A," = ",i4)') trim(reverseKey(BCH_free)) , BCH_free
        write(*,'(A," = ",i4)') trim(reverseKey(BCH_normal)) , BCH_normal
        write(*,'(A," = ",i4)') trim(reverseKey(BCH_tidal)) , BCH_tidal
        write(*,'(A," = ",i4)') trim(reverseKey(BCH_tseries)) , BCH_tseries
        write(*,'(A," = ",i4)') trim(reverseKey(BCdn)) , BCdn
        write(*,'(A," = ",i4)') trim(reverseKey(BClat)) , BClat
        write(*,'(A," = ",i4)') trim(reverseKey(BCup)) , BCup
        write(*,'(A," = ",i4)') trim(reverseKey(BCQ_fixed)) , BCQ_fixed
        write(*,'(A," = ",i4)') trim(reverseKey(BCQ_tseries)) , BCQ_tseries
        write(*,'(A," = ",i4)') trim(reverseKey(BCnone)) , BCnone
        write(*,'(A," = ",i4)') trim(reverseKey(BLink)) , BLink
        write(*,'(A," = ",i4)') trim(reverseKey(bottom_orifice)) , bottom_orifice
        write(*,'(A," = ",i4)') trim(reverseKey(BQuick)) , BQuick
        write(*,'(A," = ",i4)') trim(reverseKey(catenary)) , catenary
        write(*,'(A," = ",i4)') trim(reverseKey(CC)) , CC
        write(*,'(A," = ",i4)') trim(reverseKey(CCJM)) , CCJM
        !write(*,'(A," = ",i4)') trim(reverseKey(CD)) , CD
        write(*,'(A," = ",i4)') trim(reverseKey(circular)) , circular
        write(*,'(A," = ",i4)') trim(reverseKey(ControlCurve)) , ControlCurve
        write(*,'(A," = ",i4)') trim(reverseKey(custom)) , custom
        write(*,'(A," = ",i4)') trim(reverseKey(DarcyWeisbach)), DarcyWeisbach
        write(*,'(A," = ",i4)') trim(reverseKey(Default)) , Default
        write(*,'(A," = ",i4)') trim(reverseKey(DepthValue)), DepthValue
        write(*,'(A," = ",i4)') trim(reverseKey(diagnostic)) , diagnostic
        write(*,'(A," = ",i4)') trim(reverseKey(DiversionCurve)) , DiversionCurve
        write(*,'(A," = ",i4)') trim(reverseKey(doesnotexist)) , doesnotexist
        write(*,'(A," = ",i4)') trim(reverseKey(dummy)) , dummy
        write(*,'(A," = ",i4)') trim(reverseKey(dynamic)) , dynamic
        write(*,'(A," = ",i4)') trim(reverseKey(eggshaped)) , eggshaped
        write(*,'(A," = ",i4)') trim(reverseKey(ElemLengthAdjust)) , ElemLengthAdjust
        write(*,'(A," = ",i4)') trim(reverseKey(ETM)) , ETM
        write(*,'(A," = ",i4)') trim(reverseKey(ETM_AC)) , ETM_AC
        write(*,'(A," = ",i4)') trim(reverseKey(filled_circular)) , filled_circular
        write(*,'(A," = ",i4)') trim(reverseKey(force_main)) , force_main
        write(*,'(A," = ",i4)') trim(reverseKey(FroudeNumber)) , FroudeNumber
        write(*,'(A," = ",i4)') trim(reverseKey(func_depth_outlet)) , func_depth_outlet
        write(*,'(A," = ",i4)') trim(reverseKey(func_head_outlet)) , func_head_outlet
        write(*,'(A," = ",i4)') trim(reverseKey(FunctionalStorage)) , FunctionalStorage
        write(*,'(A," = ",i4)') trim(reverseKey(gothic)) , gothic
        write(*,'(A," = ",i4)') trim(reverseKey(HazenWilliams)), HazenWilliams
        write(*,'(A," = ",i4)') trim(reverseKey(horseshoe)) , horseshoe
        write(*,'(A," = ",i4)') trim(reverseKey(horiz_ellipse)) , horiz_ellipse
        write(*,'(A," = ",i4)') trim(reverseKey(InDays)) , InDays
        write(*,'(A," = ",i4)') trim(reverseKey(InHours)) , InHours
        write(*,'(A," = ",i4)') trim(reverseKey(InMinutes)) , InMinutes
        write(*,'(A," = ",i4)') trim(reverseKey(InSeconds)) , InSeconds
        write(*,'(A," = ",i4)') trim(reverseKey(irregular)) , irregular
        write(*,'(A," = ",i4)') trim(reverseKey(JB)) , JB
        write(*,'(A," = ",i4)') trim(reverseKey(JM)) , JM
        write(*,'(A," = ",i4)') trim(reverseKey(jump_from_upstream)) , jump_from_upstream
        write(*,'(A," = ",i4)') trim(reverseKey(jump_from_downstream)) , jump_from_downstream
        write(*,'(A," = ",i4)') trim(reverseKey(jump_none)) , jump_none
        write(*,'(A," = ",i4)') trim(reverseKey(keys_lastplusone)) , keys_lastplusone
        write(*,'(A," = ",i4)') trim(reverseKey(lArch)) , lArch
        write(*,'(A," = ",i4)') trim(reverseKey(lBasket_handle)) , lBasket_handle
        write(*,'(A," = ",i4)') trim(reverseKey(lBottomOrifice)) , lBottomOrifice
        write(*,'(A," = ",i4)') trim(reverseKey(lCatenary)) , lCatenary
        !write(*,'(A," = ",i4)') trim(reverseKey(lCD)) , lCD
        write(*,'(A," = ",i4)') trim(reverseKey(lCircular)) , lCircular
        write(*,'(A," = ",i4)') trim(reverseKey(lChannel)) , lChannel
        write(*,'(A," = ",i4)') trim(reverseKey(lEggshaped)) , lEggshaped
        write(*,'(A," = ",i4)') trim(reverseKey(lFilled_circular)) , lFilled_circular
        write(*,'(A," = ",i4)') trim(reverseKey(lForce_main)) , lForce_main
        write(*,'(A," = ",i4)') trim(reverseKey(lGothic)) , lGothic
        write(*,'(A," = ",i4)') trim(reverseKey(lHoriz_ellipse)) , lHoriz_ellipse
        write(*,'(A," = ",i4)') trim(reverseKey(lHorseshoe)) , lHorseshoe
        write(*,'(A," = ",i4)') trim(reverseKey(lIrregular)) , lIrregular
        !write(*,'(A," = ",i4)') trim(reverseKey(lManningsN)) , lManningsN
        write(*,'(A," = ",i4)') trim(reverseKey(lMod_basket)) , lMod_basket
        write(*,'(A," = ",i4)') trim(reverseKey(lNodeDepth)) , lNodeDepth
        write(*,'(A," = ",i4)') trim(reverseKey(lNodeHead)) , lNodeHead
        write(*,'(A," = ",i4)') trim(reverseKey(lOrifice)) , lOrifice
        write(*,'(A," = ",i4)') trim(reverseKey(lOutlet)) , lOutlet
        write(*,'(A," = ",i4)') trim(reverseKey(lParabolic)) , lParabolic
        write(*,'(A," = ",i4)') trim(reverseKey(lPipe)) , lPipe
        write(*,'(A," = ",i4)') trim(reverseKey(lPower_function)) , lPower_function
        write(*,'(A," = ",i4)') trim(reverseKey(lPump)) , lPump
        write(*,'(A," = ",i4)') trim(reverseKey(lRectangular)) , lRectangular
        write(*,'(A," = ",i4)') trim(reverseKey(lRectangular_closed)) , lRectangular_closed
        write(*,'(A," = ",i4)') trim(reverseKey(lRect_round)) , lRect_round
        write(*,'(A," = ",i4)') trim(reverseKey(lRect_triang)) , lRect_triang 
        write(*,'(A," = ",i4)') trim(reverseKey(lRoadWayWeir)) , lRoadWayWeir
        write(*,'(A," = ",i4)') trim(reverseKey(lSemi_circular)) , lSemi_circular
        write(*,'(A," = ",i4)') trim(reverseKey(lSemi_elliptical)) , lSemi_elliptical
        write(*,'(A," = ",i4)') trim(reverseKey(lSideFlowWeir)) , lSideFlowWeir
        write(*,'(A," = ",i4)') trim(reverseKey(lSideOrifice)) , lSideOrifice
        write(*,'(A," = ",i4)') trim(reverseKey(lTransverseWeir)) , lTransverseWeir
        write(*,'(A," = ",i4)') trim(reverseKey(lTrapezoidalWeir)) , lTrapezoidalWeir
        write(*,'(A," = ",i4)') trim(reverseKey(lTrapezoidal)) , lTrapezoidal
        write(*,'(A," = ",i4)') trim(reverseKey(lTriangular)) , lTriangular
        write(*,'(A," = ",i4)') trim(reverseKey(lType1Pump)) , lType1Pump
        write(*,'(A," = ",i4)') trim(reverseKey(lType2Pump)) , lType2Pump
        write(*,'(A," = ",i4)') trim(reverseKey(lType3Pump)) , lType3Pump
        write(*,'(A," = ",i4)') trim(reverseKey(lType4Pump)) , lType4Pump
        write(*,'(A," = ",i4)') trim(reverseKey(lTypeIdealPump)) , lTypeIdealPump
        write(*,'(A," = ",i4)') trim(reverseKey(lVert_ellipse)) , lVert_ellipse
        write(*,'(A," = ",i4)') trim(reverseKey(lVnotchWeir)) , lVnotchWeir
        write(*,'(A," = ",i4)') trim(reverseKey(lWeir)) , lWeir
        write(*,'(A," = ",i4)') trim(reverseKey(LinkOut)) , LinkOut
        write(*,'(A," = ",i4)') trim(reverseKey(manhole)) , manhole
        !write(*,'(A," = ",i4)') trim(reverseKey(ManningsN)) , ManningsN
        write(*,'(A," = ",i4)') trim(reverseKey(MaximumValue)) , MaximumValue
        write(*,'(A," = ",i4)') trim(reverseKey(mod_basket)) , mod_basket
        write(*,'(A," = ",i4)') trim(reverseKey(nBCdn)) , nBCdn
       !write(*,'(A," = ",i4)') trim(reverseKey(nBClat)) , nBClat
        write(*,'(A," = ",i4)') trim(reverseKey(nBCup)) , nBCup
        write(*,'(A," = ",i4)') trim(reverseKey(nJ1)) , nJ1
        write(*,'(A," = ",i4)') trim(reverseKey(nJ2)) , nJ2
        write(*,'(A," = ",i4)') trim(reverseKey(nJm)) , nJm
        write(*,'(A," = ",i4)') trim(reverseKey(NodeElemOut)) , NodeElemOut
        write(*,'(A," = ",i4)') trim(reverseKey(NodeFaceOut)) , NodeFaceOut
        write(*,'(A," = ",i4)') trim(reverseKey(notused)) , notused
        write(*,'(A," = ",i4)') trim(reverseKey(nStorage)) , nStorage
        write(*,'(A," = ",i4)') trim(reverseKey(orifice)) , orifice
        write(*,'(A," = ",i4)') trim(reverseKey(outlet)) , outlet
        write(*,'(A," = ",i4)') trim(reverseKey(parabolic)) , parabolic
        write(*,'(A," = ",i4)') trim(reverseKey(power_function)) , power_function
        write(*,'(A," = ",i4)') trim(reverseKey(pump)) , pump
        write(*,'(A," = ",i4)') trim(reverseKey(Pump1Curve)) , Pump1Curve
        write(*,'(A," = ",i4)') trim(reverseKey(Pump2Curve)) , Pump2Curve
        write(*,'(A," = ",i4)') trim(reverseKey(Pump3Curve)) , Pump3Curve
        write(*,'(A," = ",i4)') trim(reverseKey(Pump4Curve)) , Pump4Curve
        write(*,'(A," = ",i4)') trim(reverseKey(Random)) , Random
        write(*,'(A," = ",i4)') trim(reverseKey(RatingCurve)) , RatingCurve
        write(*,'(A," = ",i4)') trim(reverseKey(RawElemLength)) , RawElemLength
        write(*,'(A," = ",i4)') trim(reverseKey(rectangular)) , rectangular
        write(*,'(A," = ",i4)') trim(reverseKey(rectangular_closed)) , rectangular_closed
        write(*,'(A," = ",i4)') trim(reverseKey(rect_round)) , rect_round
        write(*,'(A," = ",i4)') trim(reverseKey(rect_triang)) , rect_triang
        write(*,'(A," = ",i4)') trim(reverseKey(roadway_weir)) , roadway_weir
        write(*,'(A," = ",i4)') trim(reverseKey(semi_circular)) , semi_circular
        write(*,'(A," = ",i4)') trim(reverseKey(semi_elliptical)) , semi_elliptical
        write(*,'(A," = ",i4)') trim(reverseKey(ShapeCurve)) , ShapeCurve
        write(*,'(A," = ",i4)') trim(reverseKey(side_flow)) , side_flow
        write(*,'(A," = ",i4)') trim(reverseKey(side_orifice)) , side_orifice
        write(*,'(A," = ",i4)') trim(reverseKey(SingleValue)) , SingleValue
        write(*,'(A," = ",i4)') trim(reverseKey(static)) , static
        write(*,'(A," = ",i4)') trim(reverseKey(StaticSlot)) , StaticSlot
        !write(*,'(A," = ",i4)') trim(reverseKey(storage)) , storage
        write(*,'(A," = ",i4)') trim(reverseKey(StorageCurve)) , StorageCurve
        write(*,'(A," = ",i4)') trim(reverseKey(SumElements)) , SumElements
        write(*,'(A," = ",i4)') trim(reverseKey(T00)) , T00
        write(*,'(A," = ",i4)') trim(reverseKey(T10)) , T10
        write(*,'(A," = ",i4)') trim(reverseKey(T20)) , T20
        write(*,'(A," = ",i4)') trim(reverseKey(T10s2)) , T10s2
        write(*,'(A," = ",i4)') trim(reverseKey(TA1)) , TA1
        write(*,'(A," = ",i4)') trim(reverseKey(TA2)) , TA2
        write(*,'(A," = ",i4)') trim(reverseKey(tabl_depth_outlet)) , tabl_depth_outlet
        write(*,'(A," = ",i4)') trim(reverseKey(tabl_head_outlet)) , tabl_head_outlet
        write(*,'(A," = ",i4)') trim(reverseKey(TabularStorage)) , TabularStorage
        write(*,'(A," = ",i4)') trim(reverseKey(TidalCurve)) , TidalCurve
        write(*,'(A," = ",i4)') trim(reverseKey(time_march)) , time_march
        write(*,'(A," = ",i4)') trim(reverseKey(transverse_weir)) , transverse_weir
        write(*,'(A," = ",i4)') trim(reverseKey(trapezoidal_weir)) , trapezoidal_weir
        write(*,'(A," = ",i4)') trim(reverseKey(trapezoidal)) , trapezoidal
        write(*,'(A," = ",i4)') trim(reverseKey(triangular)) , triangular
        write(*,'(A," = ",i4)') trim(reverseKey(type1_Pump)) , type1_Pump
        write(*,'(A," = ",i4)') trim(reverseKey(type2_Pump)) , type2_Pump
        write(*,'(A," = ",i4)') trim(reverseKey(type3_Pump)) , type3_Pump
        write(*,'(A," = ",i4)') trim(reverseKey(type4_Pump)) , type4_Pump
        write(*,'(A," = ",i4)') trim(reverseKey(undefinedKey)) , undefinedKey
        write(*,'(A," = ",i4)') trim(reverseKey(StaticSlot)) , StaticSlot
        write(*,'(A," = ",i4)') trim(reverseKey(DynamicSlot)) , DynamicSlot
        write(*,'(A," = ",i4)') trim(reverseKey(vert_ellipse)) , vert_ellipse
        write(*,'(A," = ",i4)') trim(reverseKey(vnotch_weir)) , vnotch_weir
        write(*,'(A," = ",i4)') trim(reverseKey(VolumeValue)) , VolumeValue
        write(*,'(A," = ",i4)') trim(reverseKey(vshape)) , vshape
        write(*,'(A," = ",i4)') trim(reverseKey(vshape_surcharge_only)) , vshape_surcharge_only
        write(*,'(A," = ",i4)') trim(reverseKey(weir)) , weir
        write(*,'(A," = ",i4)') trim(reverseKey(WeirCurve)) , WeirCurve
        !%------------------------------------------------------------------
        !% Closing
    end subroutine define_keys_printByName 
!%
!%==========================================================================
! END OF MODULE
!%==========================================================================
!%
end module define_keys
