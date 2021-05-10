! module array_index
!
! The following is based on the ontology outlined in Google Sheet
! https://docs.google.com/spreadsheets/d/126yXPLGS9F-jz6I9kuMcW_up87TurdYftQiwQTvJWSU/edit#gid=0
! This sets up all the columns that are in the elem* and face* arrays.
! Each array is associated with an integer value Ncol_xxxx,
! e.g., Ncol_elemI that provides the number of columns in the array.
! Keep in mind that we cannot use pointer at parameters and enumerated
! types (in Fortran, a pointer must be a variable). Thus, it will be
! useful to also define a set of integer vectors that simply store
! the number of columns of each array.
!
! Ben R. Hodges
! 20200411
!

module array_index

    use globals
    use iso_c_binding

    implicit none
    !%  FIRST INDEXES (DO NOT CHANGE) --------------------------------------------
    ! In theory, we can use different first index values for the arrays.
    ! This was initially used in debugging, but it seemed the gfortran compiler
    ! had some behaviors with non-unity starting points that I couldn't figure out.
    integer, parameter :: first_face_index  = 1
    integer, parameter :: first_elemM_index = 1
    integer, parameter :: first_elem2_index = 1

    !%  SIZES FOR MULTI_FACE ELEMENTS --------------------------------------------
    ! design of array depends on faces per element and elements per face
    ! CODE NOTE: face_per_elemM must equal sum of upstream and downstream
    integer, parameter :: face_per_elemM          = 6  !maximum number of faces per elemM 6,3 and 3
    integer, parameter :: upstream_face_per_elemM = 3  !maximum number of upstream faces
    integer, parameter :: dnstream_face_per_elemM = 3

    !%  linkI COLUMN INDEXES FOR INTEGER DATA OF LINKS IN LINK/NODE SYSTEM -------
    ! column index for integer data in linkI array
    enum, bind(c)
        enumerator :: li_idx = 1
        enumerator :: li_link_type
        enumerator :: li_weir_type  ! type of weir link
        enumerator :: li_orif_type  ! type of orifice link
        enumerator :: li_pump_type  ! type of pump link
        enumerator :: li_geometry
        enumerator :: li_roughness_type
        enumerator :: li_N_element  ! Number of elements in this link
        enumerator :: li_Mnode_u  ! map to upstream node connecting to link
        enumerator :: li_Mnode_d ! map to downstram node connecting to link
        enumerator :: li_Melem_u ! element ID of upstream element of link
        enumerator :: li_Melem_d ! element ID of downstream element of link
        enumerator :: li_Mface_u ! face ID of upstream face of link
        enumerator :: li_Mface_d ! face ID of downstream face of link
        enumerator :: li_assigned ! given 1 when link is assigned
        enumerator :: li_InitialDepthType ! 1=uniform, 2= lineary change, 3=exponential decay
        enumerator :: li_BQ_image ! image number assigned from BIPquick
    end enum
    integer, parameter :: li_idx_max = li_InitialDepthType

    !%  nodeI COLUMN INDEXES FOR INTEGER DATA OF NODES IN LINK/NODE SYSTEM -------
    ! column index for integer data in nodeI array
    enum, bind(c)
        enumerator :: ni_idx = 1
        enumerator :: ni_node_type
        enumerator :: ni_N_link_u ! number of upstream links at this node
        enumerator :: ni_N_link_d ! number of downstram links at this node
        enumerator :: ni_curve_type ! ID for nodal storage surface area curve type. 1 for functional and 2 for tabular
        enumerator :: ni_assigned ! given 1 when node has been assigned to face/elem,
        enumerator :: ni_total_inflow ! index to total_inflow (-1 if not total_inflow)
        enumerator :: ni_BQ_image ! image number assigned from BIPquick
        enumerator :: ni_BQ_edge  ! 0=this node has nothing to do with image communication; 1=this node is an edge
    end enum
    integer, parameter :: ni_idx_base1 = ni_total_inflow

    ! column indexes for multi-branch nodes
    integer, parameter :: ni_Mlink_u1   = ni_idx_base1+1 ! map to link of upstream branch 1
    integer, parameter :: ni_Mlink_u2   = ni_idx_base1+2 ! map to link up dowstream branch 1
    integer, parameter :: ni_Mlink_u3   = ni_idx_base1+3
    integer, parameter :: ni_idx_base2  = ni_idx_base1 + upstream_face_per_elemM

    integer, parameter :: ni_Mlink_d1   = ni_idx_base2+1
    integer, parameter :: ni_Mlink_d2   = ni_idx_base2+2
    integer, parameter :: ni_Mlink_d3   = ni_idx_base2+3
    integer, parameter :: ni_idx_max = ni_idx_base2 + dnstream_face_per_elemM

    ! storage for link index for upstream and downstream links
    integer, dimension(upstream_face_per_elemM) :: ni_MlinkUp = nullvalueI
    integer, dimension(dnstream_face_per_elemM) :: ni_MlinkDn = nullvalueI

    !%  linkR COLUMN INDEXES FOR REAL DATA OF LINKS IN LINK/NODE SYSTEM  --------
    ! column index for real data in the linkR array
    enum, bind(c)
        enumerator :: lr_Length = 1
        enumerator :: lr_InletOffset ! Every links should have a inlet and oulet offset
        enumerator :: lr_OutletOffset ! to make it consistent with SWMM.
        enumerator :: lr_BreadthScale
        enumerator :: lr_TopWidth
        enumerator :: lr_ElementLength
        enumerator :: lr_Slope
        enumerator :: lr_LeftSlope
        enumerator :: lr_RightSlope
        enumerator :: lr_Roughness
        enumerator :: lr_InitialFlowrate
        enumerator :: lr_InitialDepth
        enumerator :: lr_InitialUpstreamDepth
        enumerator :: lr_InitialDnstreamDepth
        enumerator :: lr_ParabolaValue
        enumerator :: lr_SideSlope ! for weirs only
        enumerator :: lr_DischargeCoeff1 ! discharge coefficient for triangular weir part or orifice element
        enumerator :: lr_DischargeCoeff2 ! discharge coefficient for rectangular weir part
        enumerator :: lr_FullDepth ! vertical opening of pipe, weir, orifice
        enumerator :: lr_EndContractions
        enumerator :: lr_Flowrate
        enumerator :: lr_Depth
        enumerator :: lr_DepthUp
        enumerator :: lr_DepthDn
        enumerator :: lr_Volume
        enumerator :: lr_Velocity
        enumerator :: lr_Capacity
    end enum
    integer, parameter :: lr_idx_max = lr_Capacity

    !%  nodeR COLUMN INDEXES FOR REAL DATA OF NODES IN LINK/NODE SYSTEM ----------
    ! column index for real data in the nodeR array
    enum, bind(c)
        enumerator :: nr_Zbottom = 1
        enumerator :: nr_InitialDepth
        enumerator :: nr_FullDepth
        enumerator :: nr_StorageConstant
        enumerator :: nr_StorageCoeff
        enumerator :: nr_StorageExponent
        enumerator :: nr_PondedArea
        enumerator :: nr_SurchargeDepth
        enumerator :: nr_MaxInflow
        enumerator :: nr_Eta
        enumerator :: nr_Depth
        enumerator :: nr_Volume
        enumerator :: nr_LateralInflow
        enumerator :: nr_TotalInflow
        enumerator :: nr_Flooding
    end enum
    integer, parameter :: nr_idx_base1 = nr_Flooding

    ! column index for real data on multiple branches of a node
    integer, parameter :: nr_ElementLength_u1 = nr_idx_base1 + 1 ! used for subdividing junctions
    integer, parameter :: nr_ElementLength_u2 = nr_idx_base1 + 2 ! used for subdividing junctions
    integer, parameter :: nr_ElementLength_u3 = nr_idx_base1 + 3 ! used for subdividing junctions
    integer, parameter :: nr_idx_base2 = nr_idx_base1 + upstream_face_per_elemM
    integer, parameter :: nr_ElementLength_d1 = nr_idx_base2 + 1 ! used for subdividing junctions
    integer, parameter :: nr_ElementLength_d2 = nr_idx_base2 + 2 ! used for subdividing junctions
    integer, parameter :: nr_ElementLength_d3 = nr_idx_base2 + 3 ! used for subdividing junctions
    integer, parameter :: nr_idx_max = nr_idx_base2 + dnstream_face_per_elemM

    ! storage of node indexes for multi-branch data
    integer, dimension(upstream_face_per_elemM) :: nr_ElementLengthUp = nullvalueI
    integer, dimension(dnstream_face_per_elemM) :: nr_ElementLengthDn = nullvalueI

    !%  nodeYN COLUMN INDEXES FOR LOGICAL DATA ON NODES --------------------------
    ! column index for logical data in nodeYN array
    integer, parameter :: nYN_temp1 = 1
    integer, parameter :: nYN_idx_max = 1

    !%  linkYN COLUMN INDEXES FOR LOGICAL DATA ON LINKS ---------------------------
    ! column index for logical data in linkYN array
    integer, parameter :: lYN_temp1 = 1
    integer, parameter :: lYN_idx_max = 1

    !% define the column indexes for elemI(:,:) array
    enum, bind(c)
        enumerator :: ei_Lidx = 1 !% local element index
        enumerator :: ei_Gidx !% global element index
        enumerator :: ei_HeqType !% head equation type
        enumerator :: ei_QeqType !% flowrate equation type
        enumerator :: ei_elementType !% general element type
        enumerator :: ei_specificType !% specific element type
        enumerator :: ei_geometryType !% cross-sectional geometry type
        enumerator :: ei_frictionType !% friction model
        enumerator :: ei_Mface_uG !% map to upstream face global index
        enumerator :: ei_Mface_dG !% map to downstream face global index
        enumerator :: ei_Mface_uL !% map to upstream face local index
        enumerator :: ei_Mface_dL !% map to downstream face local index
        enumerator :: ei_link_Gidx_SWMM !% link index from global SWMM network
        enumerator :: ei_node_Gidx_SWMM !% node index from global SWMM network
        enumerator :: ei_link_Gidx_BIPquick !% link index from global BIPquick network
        enumerator :: ei_node_Gidx_BIPquick !% node index from global BIPquick network
        enumerator :: ei_tmType !% time march type
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemI = ei_tmType

    !% define the column indexes for elemR(:,:)
    enum, bind(c)
        enumerator :: er_Volume = 1 !% volume (latest)
        enumerator :: er_Volume_N0 !% volume (time N)
        enumerator :: er_Volume_N1 !% volume (time N-1)
        enumerator :: er_SmallVolume !% the value of a "small volume" for this element (static)
        !% SmallVolumeRatio is ratio used in blending solved velocity with ad hoc velocity for small volume.
        enumerator :: er_SmallVolumeRatio
        enumerator :: er_FullVolume !% ?
        enumerator :: er_Flowrate !% flowrate (latest)
        enumerator :: er_Flowrate_N0 !% flowrate (time N)
        enumerator :: er_Flowrate_N1 !% flowrate (time N-1)
        enumerator :: er_Velocity !% velocity (latest)
        enumerator :: er_InterpWeight_uQ !% Interpolation Weight, Upstream, for flowrate
        enumerator :: er_InterpWeight_dQ
        enumerator :: er_InterpWeight_uG !% Interpolation Weight, Upstream for geometry
        enumerator :: er_InterpWeight_dG
        enumerator :: er_Friction !% friction term (meaning depends on friction model)
        enumerator :: er_Head !% Piezometric head (latest) -- water surface elevation in open channel
        enumerator :: er_Head_N0 !% Piezometric head (time N)
        enumerator :: er_Area !% Cross-sectional flow area (latest)
        enumerator :: er_Area_N0 !% Cross-sectional flow area (time N)
        enumerator :: er_FullArea !% full-flow area for closed conduit (static)
        enumerator :: er_Topwidth !% Topwidth of flow at free surface
        enumerator :: er_Perimeter !% Wetted perimeter of flow
        enumerator :: er_Depth !% Actual maximum depth of open-channel flow
        enumerator :: er_FullDepth !% Maximum possible flow depth in closed conduit
        enumerator :: er_HydDepth !% Hydraulic depth of flow
        enumerator :: er_HydRadius !% Hydraulic radius of flow
        enumerator :: er_X !% linear X location in system (static)
        enumerator :: er_Length !% length of element (static)
        enumerator :: er_Zbottom !% bottom elevation of element (static)
        enumerator :: er_Zcrown !% inside crown elevation of closed conduit (static)
        enumerator :: er_Roughness !% baseline roughness value for friction model (static)
        enumerator :: er_VolumeConservation !% net volume conservation for this element
        enumerator :: er_FroudeNumber !% Froude number of flow in this element
        enumerator :: er_ell !% the ell (lower case L) length scale in AC solver
        enumerator :: er_dHdA !% geometric change in elevation with area
        enumerator :: er_CtestH0 !% Convergence test for AC: d/dtau of H * A
        enumerator :: er_CtestQ0 !% Convergence test for AC: d/dtau of Q
        enumerator :: er_CtestH1 !% Convergence test for AC: d/dtau of H * A
        enumerator :: er_CtestQ1 !% Convergence test for AC: d/dtau of Q
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemR = er_CtestQ1

    enum, bind(c)
        !% columns for packed array of indexes
        enumerator :: em_CCJM = 1 !% channels, conduits, junction main
        enumerator :: em_CCJM_Htm !% channels, conduits, junction main with H time march
        enumerator :: em_CC_Qtm !% channels, conduits with Q time march
        enumerator :: em_CC_Htm !% channels, conduits with H time march
        enumerator :: em_CC_AC !% channels, conduits with AC method
        enumerator :: em_CC_Q_AC !% channels, conduits with AC method for Q time march
        enumerator :: em_CC_Q_ETM !% channels, conduits with ETM method for Q time march
        enumerator :: em_CC_H_AC !% channels, conduits with AC method for H time march
        enumerator :: em_CCJM_H_AC !% channels, conduits, junction mains with AC method for H time march
        enumerator :: em_CCJM_H_AC_open !% open channels, conduits, junction mains with AC method for H time march
        enumerator :: em_CCJM_H_AC_surcharged !% surcharged channels, conduits, junction mains with AC method for
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemM = em_CCJM_H_AC_surcharged

    enum, bind(c)
        enumerator :: ep_ALLtm = 1 !% all ETM, AC elements
        enumerator :: ep_ETM = 1 !% all ETM elements
        enumerator :: ep_AC !% all AC elements
        enumerator :: ep_Diag !% diagnostic elements
        enumerator :: ep_Surcharged_ALLtm !% all time march surcharged
        enumerator :: ep_Surcharged_ETM!% all surcharged with ETM
        enumerator :: ep_Surcharged_AC!% all surcharged with AC
        enumerator :: ep_NonSurcharged_ALLtm = 1 !% all time march nonsurcharged
        enumerator :: ep_NonSurcharged_ETM!% all surcharged with ETM
        enumerator :: ep_NonSurcharged_AC!% all surcharged with AC enumerator :: ep_AC !% all AC elements
        enumerator :: ep_JM_ALLtm !% Junction mains with all time march
        enumerator :: ep_JM_AC !% junction mains using AC method
        enumerator :: ep_JM_ETM !% junction mains using ETM method
        enumerator :: ep_smallvolume_ALLtm !% small volume with any time march
        enumerator :: ep_smallvolume_ETM !% small volume cells with ETM
        enumerator :: ep_smallvolume_AC !% small volume cells with AC
        enumerator :: ep_CC_ALLTm !% all CC elements that are ETM or AC
        enumerator :: ep_CC_ETM !% all CC elements that are ETM
        enumerator :: ep_CC_AC !% all CC elements that are AC
        enumerator :: ep_CCJB_eETM_i_fAC !% CC and junction branches where element ETM has face AC
        enumerator :: ep_CCJB_eAC_i_fETM !% CC and junction branches where element AC has face ETM
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemP = ep_CCJB_eAC_i_fETM

    !% define the column indexes for elemYN(:,:)
    enum, bind(c)
        enumerator :: eYN_isClosedTop = 1 !% TRUE is closed conduit, FALSE is open channel
        enumerator :: eYN_isSmallVolume !% TRUE is use small volume algorithm
        enumerator :: eYN_isAdhocFlowrate !% TRUE is use ad hoc flowrate algorithm
        enumerator :: eYN_isSurcharged !% TRUE is a surcharged conduit, FALSE is open channel flow
        enumerator :: eYN_canSurcharge !% TRUE for element that can surcharge, FALSE where it cannot
        enumerator :: eYN_isMultiface !% TRUE for elements that use the multiface junction framework
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemYN = eYN_isMultiface

    !% define the column indexes for elemSI(:,:) junctions only
    enum, bind(c)
        enumerator :: eSi_nfaces_u = 1 !% Junction
        enumerator :: eSi_nfaces_d !% Junction
        enumerator :: eSI_storage_curveType !% Junction
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSI = eSI_storage_curveType

    !% define the column indexes for elemSR(:,:) for circular geometry
    enum, bind(c)
        enumerator :: eSr_Circular_Radius = 1 !% radius for circular geometry (Circular Cross-Section)
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSR_circular = eSr_Circular_Radius

    !% define the column indexes for elemSR(:,:) for rectangular pipe or channel
    enum, bind(c)
        enumerator :: eSr_Rectangular_Breadth = 1 !% breadth for rectangular geometry
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSR_rectangular = eSr_Rectangular_Breadth

    !% define the column indexes for elemSR(:,:) for trapezodial geometry
    enum, bind(c)
        enumerator :: eSr_Trapezoidal_Breadth = 1 !% bottom breadth for trapezoidal geometry
        enumerator :: eSr_LeftSlope !% left inverse slope of trapezoid (looking downstream)
        enumerator :: eSr_RightSlope !% right inverse slope of trapezoid (looking downstream)
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSR_trapezoidal = eSr_RightSlope

    !% define the column indexes for elemSR(:,:) for parabolic geometry
    enum, bind(c)
        enumerator :: eSr_ParabolaValue = 1 !% ?
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSR_parabolic = eSr_ParabolaValue

    !% define the column indexes for elemSR(:,:) for geometry that has not yet been confirmed and assigned:
    enum, bind(c)
        enumerator :: eSr_DischargeCoeff1 = 1 !% discharge coefficient for triangular weir part or orifice element
        enumerator :: eSr_DischargeCoeff2 !% discharge coefficient for rectangular weir part
        enumerator :: eSr_EndContractions !% End contractions for rectengular and trapezoidal weir
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSR_weir_orifice = eSr_EndContractions

    !% find the max array size for any type of special element data stored in elemSR
    integer, parameter :: Ncol_elemSR = &
        max(Ncol_elemSR_circular, &
            Ncol_elemSR_rectangular, &
            Ncol_elemSR_trapezoidal, &
            Ncol_elemSR_parabolic, &
            Ncol_elemSR_weir_orifice)

    !% define the column indexes for faceI(:,:)
    enum, bind(c)
        enumerator :: fi_Lidx = 1 !% local array index (row)
        enumerator :: fi_Gidx !% global (unique) index
        enumerator :: fi_Melem_uG !% Map to element upstream (global index)
        enumerator :: fi_Melem_dG !% Map to element upstream (local index)
        enumerator :: fi_Melem_uL !% Map to element upstream (local index)
        enumerator :: fi_Melem_dL !% Map to element downstream (local index)
        enumerator :: fi_eHeqType_u !% Type of H solution on element upstream
        enumerator :: fi_eHeqType_d !% Type of H solution on element downstream
        enumerator :: fi_eQeqType_u !% Type of Q solution on element upstream
        enumerator :: fi_eQeqType_d !% Type of Q solution on element downstream
        enumerator :: fi_jump_type !% Type of hydraulic jump
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_faceI = fi_jump_type

    !% define the column indexes for faceR(:,:)
    enum, bind(c)
        enumerator :: fr_Area_u = 1 !% cross-sectional area on upstream side of face
        enumerator :: fr_Area_d !% cross-sectional area on downstream side of face
        enumerator :: fr_Head_u !% Piezometric head on upstream side of face
        enumerator :: fr_Head_d !% Piezometric head on downstream side of face
        enumerator :: fr_Flowrate !% Flowrate through face (latest)
        enumerator :: fr_Flowrate_N0 !% Flowrate through face (time N)
        enumerator :: fr_Flowrate_net !% ?
        enumerator :: fr_HydDepth_u !% Hydraulic Depth on upstream side of face
        enumerator :: fr_HydDepth_d !% Hydraulic Depth on downstream side of face
        enumerator :: fr_Topwidth_u !% Topwidth on upstream side of face
        enumerator :: fr_Topwidth_d !% Topwidth on downstream side of face
        enumerator :: fr_Velocity_u !% Velocity on upstream side of face
        enumerator :: fr_Velocity_d !% Velocity on downstream side of face
        enumerator :: fr_Zbottom_u !% Bottom elevation on upstream side of face
        enumerator :: fr_Zbottom_d !% Bottom elevation on downstream side of face
        enumerator :: fr_X !% Linear X location in system
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_faceR = fr_X

    enum, bind(c)
        enumerator :: fm_all = 1
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_faceM = fm_all

    enum, bind(c)
        enumerator :: fp_Diag = 1 !% face with adjacent diagnostic element
        enumerator :: fp_JumpUp !% face with hydraulic jump from nominal upstream to downstream
        enumerator :: fp_JumpDn !% face with hydraulic jump from nominal downstream to upstream
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_faceP = fp_JumpDn

    enum, bind(c)
        enumerator :: fYN_isDiag_adjacent = 1
        enumerator :: fYN_isBCface
        enumerator :: fYN_isETM_adjacent
        enumerator :: fYN_isAC_adjacent
        enumerator :: fYN_isnull
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_faceYN = fYN_isnull
end module array_index