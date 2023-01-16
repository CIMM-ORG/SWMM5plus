! module define_indexes
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
! Note that indexes marked with KEY must have values in the range 1...undefinedKey
! as enuerated in the define_keys module

module define_indexes

    use define_globals
    !use iso_c_binding

    implicit none
!%
!%==========================================================================
!% LINKS
!%==========================================================================
!%
    !%-------------------------------------------------------------------------
    !% Define the column indexes for link%I(:,:) arrays
    !% These are the for the full arrays of integer data
    !%-------------------------------------------------------------------------
    enum, bind(c)
        enumerator :: li_idx = 1
        enumerator :: li_link_type           ! KEY type of links (i.e. conduit, orifice, weir, etc.)   
        enumerator :: li_link_sub_type       ! KEY link subtype (i.e. vnotch weir, side orifice, etc.)
        enumerator :: li_link_direction      ! link direction
        enumerator :: li_geometry            ! KEY link geometry type
        enumerator :: li_barrels             ! KEY link # of barrels
        enumerator :: li_culvertCode         ! KEY culvert code for conduit
        !enumerator :: li_roughness_type obsolete 20220708brh
        enumerator :: li_N_element           ! Number of elements in this link
        enumerator :: li_Mnode_u             ! map to upstream node connecting to link
        enumerator :: li_Mnode_d             ! map to downstram node connecting to link
        enumerator :: li_assigned            ! given 1 when link is assigned
        enumerator :: li_InitialDepthType    ! KEY UniformDepth, LinearlyVaryingDepth, ExponentialDepth, FixedHead
        enumerator :: li_length_adjusted     ! 1 = length was not adjusted, 2 = one side was adjusted, 3 = both side was adjusted
        enumerator :: li_P_image             ! image number assigned from BIPquick
        enumerator :: li_parent_link         ! A map to the corresponding SWMM link after a BIPquick link-split
        !enumerator :: li_num_phantom_links   ! Number of phantom links associated 
        enumerator :: li_weir_EndContractions ! (0,1) to indicate contraction
        enumerator :: li_RoadSurface           ! roadsurface type for roadway weir
        enumerator :: li_curve_id            ! curve id if the link is associated with any curve
        enumerator :: li_first_elem_idx
        enumerator :: li_last_elem_idx
        enumerator :: li_transect_idx         ! transect index if the link is associated with an irregular geometry transect
        enumerator :: li_lastplusone !% must be last enum item
    end enum
    integer, target :: Ncol_linkI = li_lastplusone-1

    !%-------------------------------------------------------------------------
    !% Define the column indexes for link%R(:,:) arrays
    !% These are the for the full arrays of real data
    !%-------------------------------------------------------------------------
    enum, bind(c)
        enumerator :: lr_Length = 1
        enumerator :: lr_AdjustedLength ! length adjustment if multi-link junction is present
        enumerator :: lr_InletOffset    ! Every links should have a inlet and oulet offset
        enumerator :: lr_OutletOffset   ! to make it consistent with SWMM.
        enumerator :: lr_FullArea
        enumerator :: lr_FullHydRadius
        enumerator :: lr_BottomDepth 
        enumerator :: lr_BottomRadius  
        enumerator :: lr_BreadthScale
        enumerator :: lr_TopWidth
        enumerator :: lr_ElementLength
        enumerator :: lr_Slope
        enumerator :: lr_LeftSlope
        enumerator :: lr_RightSlope
        enumerator :: lr_Roughness
        enumerator :: lr_Kentry_MinorLoss           !% K factor for entry minor loss
        enumerator :: lr_Kexit_MinorLoss            !% K factor for exit minor loss
        enumerator :: lr_Kconduit_MinorLoss         !% K factor over the body of the conduit
        enumerator :: lr_SeepRate                   !% seepage rate (converted to m/s)
        enumerator :: lr_FlowrateInitial
        enumerator :: lr_FlowrateLimit           ! user.inp file Qmax (0 is does not apply)
        enumerator :: lr_ForceMain_Coef
        !enumerator :: lr_InitialDepth
        enumerator :: lr_InitialUpstreamDepth
        enumerator :: lr_InitialDnstreamDepth
        enumerator :: lr_ParabolaValue
        enumerator :: lr_SideSlope             ! for weirs only
        enumerator :: lr_DischargeCoeff1       ! discharge coefficient for triangular weir part or orifice element
        enumerator :: lr_DischargeCoeff2       ! discharge coefficient for rectangular weir part
        enumerator :: lr_RoadWidth             ! road width for roadway weir
        enumerator :: lr_initSetting           ! initial pump speed setting 
        enumerator :: lr_yOn                   ! startup depth for pumps   
        enumerator :: lr_yOff                  ! shutoff depth for pumps   
        enumerator :: lr_FullDepth             ! vertical opening of pipe, weir, orifice
        enumerator :: lr_Setting               !% the 0 to 1 open/close setting of EPA-SWMM
        enumerator :: lr_TargetSetting         !% target setting of a control action
        enumerator :: lr_TimeLastSet           !% the time (in seconds) the link setting was last changed
        !enumerator :: lr_Flowrate
        !enumerator :: lr_Depth
        !enumerator :: lr_DepthUp
        !enumerator :: lr_DepthDn
        !enumerator :: lr_Volume
        !enumerator :: lr_Velocity
        !enumerator :: lr_Capacity
        enumerator :: lr_ZbottomUp             ! Z bottom of upstream node
        enumerator :: lr_ZbottomDn             ! Z bottom of downstream node
        enumerator :: lr_lastplusone !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is
    integer, target :: Ncol_linkR = lr_lastplusone-1

    !%-------------------------------------------------------------------------
    !% Define the column indexes for link%YN(:,:) arrays
    !% These are the for the full arrays of logical
    !%-------------------------------------------------------------------------
    enum, bind(c)
        enumerator :: lYN_weir_CanSurcharge = 1
        enumerator :: lYN_isOutput
        enumerator :: lYN_isPhantomLink
        enumerator :: lYN_hasFlapGate
        enumerator :: lYN_temp1
        enumerator :: lYN_lastplusone !% must be last enum item
    end enum
    integer, target :: Ncol_linkYN  = lYN_lastplusone-1

!%
!%==========================================================================
!% NODES
!%==========================================================================
!%    
    !%-------------------------------------------------------------------------
    !% Define the column indexes for node%I(:,:) arrays
    !% These are the for the full arrays of integer data
    !%-------------------------------------------------------------------------
    enum, bind(c)
        enumerator :: ni_idx = 1
        enumerator :: ni_node_type     ! KEY
        enumerator :: ni_N_link_u      ! number of upstream links at this node
        enumerator :: ni_N_link_d      ! number of downstram links at this node
        enumerator :: ni_curve_ID      ! ID for nodal storage surface area curve
        enumerator :: ni_assigned      ! given 1 when node has been assigned to face/elem,
        enumerator :: ni_P_image       ! image number assigned from BIPquick
        enumerator :: ni_P_is_boundary ! 0=this node has nothing to do with image communication; >0=this node is a partition boundary
        ! if node is BCup or BCdn, ni_elemface_idx is the index of its associated BC face
        ! if node is nJm or nJ2, ni_elemface_idx is the index of the associated element
        !enumerator :: ni_elemface_idx  ! OBSOLETE
        enumerator :: ni_elem_idx      !% this is the element of an nJM node, upstream element of BCdn, downstream element of BCup
        enumerator :: ni_face_idx      !% for nJ2, BCup, BCdn, nJ1, this is the face associated with the node, not defined for nJM
        enumerator :: ni_pattern_resolution ! minimum resolution of patterns associated with node BC
        enumerator :: ni_routeTo       !% subcatchment indx that node outfall is routed to
        enumerator :: ni_routeFrom     !% subcatchment indx with outlet to this node
        enumerator :: ni_SWMMoutfallIdx !% Outfall index in EPA SWMM for an outfall node
        enumerator :: ni_lastplusone !% must be last enum item
    end enum
    integer, parameter :: ni_idx_base1 = ni_lastplusone-1

    !% column indexes for multi-branch nodes
    integer, parameter :: ni_Mlink_u1   = ni_idx_base1+1 ! map to link of upstream branch 1
    integer, parameter :: ni_Mlink_u2   = ni_idx_base1+2 
    integer, parameter :: ni_Mlink_u3   = ni_idx_base1+3
    integer, parameter :: ni_Mlink_u4   = ni_idx_base1+4               !% ADDBRANCH
    integer, parameter :: ni_Mlink_u5   = ni_idx_base1+5               !% ADDBRANCH

    integer, parameter :: ni_idx_base2  = ni_idx_base1 + max_branch_per_node/2

    integer, parameter :: ni_Mlink_d1   = ni_idx_base2+1
    integer, parameter :: ni_Mlink_d2   = ni_idx_base2+2
    integer, parameter :: ni_Mlink_d3   = ni_idx_base2+3
    integer, parameter :: ni_Mlink_d4   = ni_idx_base2+4               !% ADDBRANCH
    integer, parameter :: ni_Mlink_d5   = ni_idx_base2+5               !% ADDBRANCH

    !% start of the multi-branch connections, needed in partitioning
    integer, parameter :: ni_MlinkStart = ni_Mlink_u1

    !% end of the multi-branch connections -- update for change in # of branches
    !integer, parameter :: ni_MlinkEnd   = ni_Mlink_d3                 !% ADDBRANCH
    !integer, parameter :: ni_MlinkEnd   = ni_Mlink_d4                  !% ADDBRANCH
    integer, parameter :: ni_MlinkEnd   = ni_Mlink_d5                  !% ADDBRANCH

    !% storage for link index for upstream and downstream links
    integer, dimension(max_branch_per_node/2) :: ni_MlinkUp = nullvalueI
    integer, dimension(max_branch_per_node/2) :: ni_MlinkDn = nullvalueI

    integer, target :: Ncol_nodeI = ni_idx_base2 + max_branch_per_node/2

    !%-------------------------------------------------------------------------
    !% Define the column indexes for node%R(:,:) arrays
    !% These are the for the full arrays of real data
    !%-------------------------------------------------------------------------
    enum, bind(c)
        enumerator :: nr_Zbottom = 1
        enumerator :: nr_InitialDepth
        enumerator :: nr_FullDepth
        enumerator :: nr_StorageConstant
        enumerator :: nr_StorageCoeff
        enumerator :: nr_StorageExponent
        enumerator :: nr_StorageFevap
        enumerator :: nr_PondedArea
        enumerator :: nr_SurchargeExtraDepth
        enumerator :: nr_MaxInflow
        enumerator :: nr_Eta
        enumerator :: nr_Depth
        enumerator :: nr_head
        enumerator :: nr_Volume
        enumerator :: nr_Flooding
        enumerator :: nr_JunctionBranch_Kfactor
        enumerator :: nr_lastplusone !% must be last enum item
    end enum
    integer, parameter :: nr_idx_base1 = nr_lastplusone-1

    integer, target :: Ncol_nodeR = nr_idx_base1

    !%-------------------------------------------------------------------------
    !% Define the column indexes for node%YN(:,:) arrays
    !% These are the for the full arrays of logical
    !%-------------------------------------------------------------------------
    enum, bind(c)
        enumerator :: nYN_has_inflow = 1
        enumerator :: nYN_has_extInflow
        enumerator :: nYN_has_dwfInflow
        enumerator :: nYN_has_storage
        enumerator :: nYN_isOutput
        enumerator :: nYN_is_phantom_node
        enumerator :: nYN_hasFlapGate
        enumerator :: nYN_lastplusone !% must be last enum item
    end enum
    integer, target :: Ncol_nodeYN  = nYN_lastplusone-1

!%
!%==========================================================================
!% BOUNDARY CONDITIONS
!%==========================================================================
!%  
    !% Column indexes for BC%xR(:,:) where x is head or flow
    enum, bind(c)
        enumerator :: br_value = 1    !% interpolated value for BC at this time step
        enumerator :: br_timeInterval !% time interval for latest forcing data
        enumerator :: br_Temp01       !% temporary array
        enumerator :: br_lastplusone  !% must be last enum item
    end enum

    !% Column indexes for BC%xI(:,:) where x is head or flow
    enum, bind(c)
        enumerator :: bi_idx = 1
        enumerator :: bi_node_idx
        enumerator :: bi_face_idx    ! Index of face nBCup/dn nodes
        enumerator :: bi_elem_idx    ! Index of element associated with either nJ2 or nJm node with lateral inflow
        enumerator :: bi_category    ! KEY
        enumerator :: bi_subcategory ! KEY
        enumerator :: bi_fetch       ! 1 if BC%xR_timeseries needs to be fetched, 0 otherwise
        enumerator :: bi_TS_upper_idx  !% index of the current level in the timeseries storage
        enumerator :: bi_UTidx       !% index in uniform table array associated with this BC.
        enumerator :: bi_lastplusone !% must be last enum item
    end enum

    !% Column indexes for BC%xYN(:,:)
    enum, bind(c)
        enumerator :: bYN_read_input_series = 1
        enumerator :: bYN_hasFlapGate
        enumerator :: bYN_lastplusone !% must be last enum item
    end enum

    !% Column indexes (3rd index) for BC%xTimeseries(:,:,:) where X is flow or head
    enum, bind(c)
        enumerator :: brts_time = 1
        enumerator :: brts_value
        enumerator :: brts_lastplusone !% must be last enum item
    end enum


    !% HACK - we will probably want to create a different set of indexes for BC%flowI and BC%headI tables
    !% For instance, BC%flowI tables will probably need addititonal information to distribute flowrates
    !% over link elements.
    integer, parameter :: N_flowI = bi_lastplusone-1
    integer, parameter :: N_headI = bi_lastplusone-1
    integer, parameter :: N_flowYN = bYN_lastplusone-1
    integer, parameter :: N_headYN = bYN_lastplusone-1
    integer, parameter :: N_flowR  = br_lastplusone-1
    integer, parameter :: N_headR  = br_lastplusone-1
    integer, parameter :: N_flowR_TS = brts_lastplusone - 1
    integer, parameter :: N_headR_TS = brts_lastplusone - 1 

!%
!%==========================================================================
!% ELEMENTS
!%==========================================================================
!%  
    !%-------------------------------------------------------------------------
    !% Define the column indexes for elemI(:,:) array
    !% These are for the full array of all integers
    !%-------------------------------------------------------------------------
    enum, bind(c)
         enumerator :: ei_Lidx = 1                  !% local element index (static)
         enumerator :: ei_Gidx                      !% global element index  (static)
         !enumerator :: ei_main_idx_for_branch       !% idx of JM for a JB branch
         enumerator :: ei_elementType               !% KEY general element type  (static)
         enumerator :: ei_geometryType              !% KEY cross-sectional geometry type  (static)
         enumerator :: ei_barrels                   !% Integer number of barrels
         enumerator :: ei_HeqType                   !% KEY type of head equation (static)
         enumerator :: ei_link_Gidx_SWMM            !% link index from global SWMM network  (static)
         enumerator :: ei_link_Gidx_BIPquick        !% link index from global BIPquick network  (static)
         enumerator :: ei_link_pos                  !% position (elem from upstream = 1 to downstream = n) in link
         enumerator :: ei_Mface_uL                  !% map to upstream face local index  (static)
         enumerator :: ei_Mface_dL                  !% map to downstream face local index  (static)
         enumerator :: ei_node_Gidx_SWMM            !% node index from global SWMM network  (static)
         enumerator :: ei_node_Gidx_BIPquick        !% node index from global BIPquick network  (static)
         enumerator :: ei_QeqType                   !% KEY type of flow equation (static)
         ! enumerator :: ei_specificType              !% specific element type (static) NOT USED AS OF 20220626
         !% brh20211210s
         !enumerator :: ei_Subcatch_TableIdx         !% index in subcatchment table for linking to subcatchments to this element 
         !enumerator :: ei_Nsubcatch                 !% number of subcatchments feeding an element
         !% brh20211210e         
         enumerator :: ei_Temp01                    !% temporary array
         enumerator :: ei_Temp02
         enumerator :: ei_tmType                    !% KEY time march type (dynamic)
         enumerator :: ei_BoundaryArray_idx         !% if a boundary cell, then position in the elemB array
         enumerator :: ei_link_transect_idx         !% index of the link transect in link array
         enumerator :: ei_transect_idx              !% index of transect in transect array
         enumerator :: ei_lastplusone !% must be last enum item
    end enum
    integer, target :: Ncol_elemI = ei_lastplusone-1

    !%-------------------------------------------------------------------------
    !% Define the column indexes for elemR(:,:) array
    !% These are for the full arrays of all reals
    !%-------------------------------------------------------------------------
    enum, bind(c)
        enumerator :: er_Area = 1                   !% cross-sectional flow area (latest)
        enumerator :: er_Area_N0                    !% cross-sectional flow area (time N)
        enumerator :: er_Area_N1                    !% cross-sectional flow area (time N-1)
        enumerator :: er_AreaBelowBreadthMax        !% area below the max breadth in a conduit (static)
        enumerator :: er_AoverAfull                 !% ratio of are to full area -- used for tabular geometry
        !enumerator :: er_Beta                      !% bottom slope / roughness, so that Q/beta = section factor
        enumerator :: er_BottomSlope                !% bottom slope of the element
        enumerator :: er_BreadthMax                 !% maximum breadth of conduit (static)
        enumerator :: er_Depth                      !% actual maximum depth of open-channel flow
        enumerator :: er_DepthAtBreadthMax          !% depth below the point of maximum breadth
        enumerator :: er_dHdA                       !% geometric change in elevation with area (used in AC only)
        enumerator :: er_dSlotArea                  !% change in slot volume
        enumerator :: er_dSlotDepth                 !% change in slot depth
        enumerator :: er_dSlotVolume                !% change in slot volume
        enumerator :: er_EllDepth                        !% the ell (lower case L) modified hydraulic depth
        enumerator :: er_Flowrate                   !% flowrate (latest)
        enumerator :: er_FlowrateLimit               !% max flowrate from user.inp file (0 is no limit)
        enumerator :: er_Flowrate_N0                !% flowrate (time N)
        enumerator :: er_Flowrate_N1                !% flowrate (time N-1)
        enumerator :: er_FlowrateLateral            !% lateral inflow BC
        enumerator :: er_FlowrateStore              !% temporary storage used for adjacent AC and ETM elements
        enumerator :: er_FroudeNumber               !% froude number of flow
        enumerator :: er_FullArea                   !% cross-sectional area of a full conduit (static)
        enumerator :: er_FullDepth                  !% maximum possible flow depth in full conduit (static)
        !enumerator :: er_FullEllDepth                    !% ell of  full pipe
        !enumerator :: er_FullHydDepth               !% hydraulic (average) depth of full conduit (static)
        enumerator :: er_FullHydRadius              !% hydraulic radius of full conduit (static)
        enumerator :: er_FullPerimeter              !% wetted perimeter of full conduit or channel (static)
        enumerator :: er_FullTopwidth               !% Topwidth of full conduit or channel
        enumerator :: er_FullVolume                 !% Volume of a full conduit or channel (static)
        enumerator :: er_GammaC                     !% gamma continuity source term for AC solver
        enumerator :: er_GammaM                     !% gamma momentum source term for AC solver
        enumerator :: er_Head                       !% piezometric head (latest) -- water surface elevation in open channel
        enumerator :: er_Head_N0                    !% piezometric head (time N)
        !enumerator :: er_HeadAvg                    !% average of head on faces of an element.
        enumerator :: er_HeadLastAC                 !% piezometric head at start of last AC step
        enumerator :: er_HeadStore                  !% temporary storage used for adjacent AC and ETM elements
        !enumerator :: er_HydDepth                   !% hydraulic depth of flow
        enumerator :: er_HydRadius                  !% hydraulic radius of flow
        enumerator :: er_InterpWeight_uG            !% interpolation Weight, upstream, for geometry
        enumerator :: er_InterpWeight_dG            !% interpolation Weight, downstream, for geometry
        enumerator :: er_InterpWeight_uH            !% interpolation Weight, upstream for head
        enumerator :: er_InterpWeight_dH            !% interpolation Weight, downstream for head
        enumerator :: er_InterpWeight_uQ            !% interpolation Weight, upstream, for flowrate
        enumerator :: er_InterpWeight_dQ            !% interpolation Weight, downstream, for flowrate
        enumerator :: er_InterpWeight_uP            !% interpolation Weight, upstream, for Preissmann number
        enumerator :: er_InterpWeight_dP            !% interpolation Weight, downstream, for Preissmann number
        enumerator :: er_Ksource                    !% k source term for AC solver
        enumerator :: er_Kentry_MinorLoss           !% K factor for entry minor loss
        enumerator :: er_Kexit_MinorLoss            !% K factor for exit minor loss
        enumerator :: er_Kconduit_MinorLoss         !% K factor over the body of the conduit
        enumerator :: er_Length                     !% length of element (static)
        enumerator :: er_ones                       !% vector of ones (useful with sign function)
        enumerator :: er_Perimeter                  !% Wetted perimeter of flow
        enumerator :: er_Preissmann_Celerity        !% celerity due to Preissmann Slot
        enumerator :: er_Preissmann_Number          !% Preissmann number
        enumerator :: er_Preissmann_Number_N0       !% preissmann number at previous time step
        enumerator :: er_Preissmann_Number_initial  !% Initial Preissmann number before surcharge
        !enumerator :: er_Pressure_Head              !% Pressure head 
        enumerator :: er_ManningsN                  !% baseline Mannings N roughness value for friction model
        enumerator :: er_ManningsN_Dynamic          !% total ManningsN roughness, including dynamic adjustment (experimental)
        enumerator :: er_SedimentDepth
        enumerator :: er_SeepRate                   !% Local seepage rate in m/s
        enumerator :: er_Setting                    !% percent open setting for a link element, on (1) or off (0) for pump
        !enumerator :: er_SectionFactor              !% present value of Qn/S0 section factor
        !enumerator :: er_SectionFactor_Max          !% maximum value of section factor (for S0 = 0)
        enumerator :: er_SlotWidth                  !% slot width
        enumerator :: er_SlotDepth                  !% slot depth
        enumerator :: er_SlotDepth_N0
        enumerator :: er_SlotArea                   !% slot area
        enumerator :: er_SlotHydRadius              !% slot hydraulic radius 
        enumerator :: er_SlotVolume                 !% slot volume 
        enumerator :: er_SlotVolume_N0              !% old slot volume      
        enumerator :: er_SmallVolume                !% the value of a "small volume" for this element
        enumerator :: er_SmallVolume_CMvelocity     !% velocity by Chezy-Manning for a small volume
        enumerator :: er_SmallVolume_ManningsN      !% roughness used for computing Chezzy-Manning on small volume
        enumerator :: er_SmallVolumeRatio           !% blending ad hoc and solved velocity for small volume.
        enumerator :: er_SourceContinuity           !% source term for continuity equation
        enumerator :: er_SourceMomentum             !% source term for momentum equation
        enumerator :: er_Surcharge_Time             !% time a CC element is surcharged
        enumerator :: er_TargetSetting              !% target percent open setting for a link element in the next time step
        enumerator :: er_Temp01                     !% temporary array (use and set to null in a single procedure)
        enumerator :: er_Temp02                     !% temporary array (use and set to null in a single procedure)
        enumerator :: er_Temp03                     !% temporary array (use and set to null in a single procedure)
        enumerator :: er_Temp04                     !% temporary array (use and set to null in a single procedure)
        enumerator :: er_TimeLastSet                !% last time the er_Setting was changed
        enumerator :: er_Topwidth                   !% topwidth of flow at free surface
        enumerator :: er_Velocity                   !% velocity (latest)
        enumerator :: er_Velocity_N0                !% velocity time N
        enumerator :: er_Velocity_N1                !% velocity time N-1
        enumerator :: er_VelocityLastAC             !% velocity at start of last AC step
        enumerator :: er_Volume                     !% volume (latest)
        enumerator :: er_Volume_N0                  !% volume (time N)
        enumerator :: er_Volume_N1                  !% volume (time N-1)
        enumerator :: er_VolumeConservation         !% cumulative volume conservation
        enumerator :: er_VolumeLastAC               !% volume at start of last AC step
        enumerator :: er_VolumeOverFlow             !% volume lost for overflow in this time step.  20220124brh
        enumerator :: er_VolumeOverFlowTotal        !% total volume lost to overflow       20220124brh 
        enumerator :: er_VolumeArtificialInflow     !% artificial inflow in transition to zero depth in this step
        enumerator :: er_VolumeArtificialInflowTotal
        enumerator :: er_VolumePonded               !% total volume in ponding
        enumerator :: er_VolumeStore                !% temporary storage used for adjacent AC and ETM elements
        enumerator :: er_WaveSpeed                  !% wave speed in element
        enumerator :: er_YoverYfull                 !% ratio of depth to full depth for tabular geometry
        enumerator :: er_Zbottom                    !% bottom elevation of element (static)
        enumerator :: er_ZbreadthMax                !% elevation at maximum breadth
        enumerator :: er_Zcrown                     !% inside crown elevation of closed conduit (static)
        
        enumerator :: er_lastplusone !% must be last enum item
    end enum
    integer, target :: Ncol_elemR = er_lastplusone-1

    !%-------------------------------------------------------------------------
    !% Define the column indexes for elemYN(:,:) arrays
    !% These are the for the full arrays of logical
    !%-------------------------------------------------------------------------

    enum, bind(c)
        enumerator :: eYN_canSurcharge = 1              !% TRUE for element that can surcharge, FALSE where it cannot (static)
        enumerator :: eYN_hasSubcatchRunOff             !% TRUE if element connected to one or more subcatchments for Runoff
        enumerator :: eYN_hasFlapGate                   !% TRUE if 1-way flap gate is present
        enumerator :: eYN_isBoundary_up                 !% TRUE if the element is connected to a shared face upstream thus a boundary element of a partition
        enumerator :: eYN_isBoundary_dn                 !% TRUE if the element is connected to a shared face downstream thus a boundary element of a partition
        enumerator :: eYN_isCulvert                     !% TRUE if element is inlet, outlet or culvert barrel
        enumerator :: eYN_isDownstreamJB                !% TRUE if the element is downstream JB
        enumerator :: eYN_isDummy
        enumerator :: eYN_isElementDownstreamOfJB       !% TRUE if the element is immediate downstream of JB
        enumerator :: eYN_isForceMain                   !% TRUE if this is a force main element
        enumerator :: eYN_isOutput                      !% TRUE if the element is an output element
        enumerator :: eYN_isPSsurcharged                !% TRUE if Preissmann slot is present for this cell
        enumerator :: eYN_isSmallDepth                  !% TRUE is use small volume algorithm
        enumerator :: eYN_isSurcharged                 !% TRUE is a surcharged closed conduit, FALSE if non-surcharged 
        enumerator :: eYN_isZeroDepth                   !% TRUE if volume qualifies as "near zero"
        enumerator :: eYN_Temp01                        !% temporary logical space
        enumerator :: eYN_lastplusone !% must be last enum item
    end enum
    integer, target :: Ncol_elemYN = eYN_lastplusone-1

!%
!%==========================================================================
!% PACKED ELEMENTS
!%==========================================================================
!%  
    !%-------------------------------------------------------------------------
    !% Define the column indexes for elemP(:,:) array
    !% These are the for the packed arrays general elements
    !%-------------------------------------------------------------------------
    enum, bind(c)
        enumerator :: ep_AC = 1                     !% all AC elements
        enumerator :: ep_ALLtm                      !% all ETM, AC elements 
        enumerator :: ep_CC_AC                      !% all CC elements that are AC
        !enumerator :: ep_CC_ACsurcharged            !% all CC elements that are AC
        enumerator :: ep_CC_ALLtm                   !% all CC elements that are ETM or AC
        !enumerator :: ep_CC_ALLtm_ACsurcharged      !% all CC elements that are AC and surcharged
        enumerator :: ep_CC_ETM                     !% all CC elements that are ETM
        !enumerator :: ep_CC_ETM_PSsurcharged        !% CC elements that are ETM and surcharged
        enumerator :: ep_CC_H_ETM                   !% all CC elements that are ETM for H
        enumerator :: ep_CC_Q_AC                    !% all CC elements that are AC for Q
        enumerator :: ep_CC_Q_ETM                   !% all CC elements that are ETM for Q
        !enumerator :: ep_CCJB_ACsurcharged          !% all CC and JB elements that are AC
        !enumerator :: ep_CCJB_ALLtm                 !% all CC and JB elements that ar any TM
        !enumerator :: ep_CCJB_AC                    !% all CC and JB elements that are AC
        !enumerator :: ep_CCJB_ALLtm_ACsurcharged      !% all CC and JB elements that are AC and surcharged
        enumerator :: ep_CCJB_eETM_i_fAC            !% Any CC or JB element that is ETM and has an AC face.
        !enumerator :: ep_CCJB_ETM                   !% CC and JB elements that are ETM
        !enumerator :: ep_CCJB_ETM_PSsurcharged        !% CC and JB elements that are ETM and surcharged
        !enumerator :: ep_CCJM_H_AC_open             !% CC and JM elements that are AC for H and open channel
        enumerator :: ep_CCJM_H_ETM                 !% CC and JM elements that are ETM for H
        enumerator :: ep_CC_isClosedSetting          !% CC elements that have er_Setting = 0.0 indicating closed off
        enumerator :: ep_Culvert_Inlet              !% all CC elements that are also culvert inlets
        !enumerator :: ep_culvert_outlet             !% all CC elements that are also culvert outlets
        !enumerator :: ep_culvert_inout              !% all CC elements that are both culvert inlet and outlet
        enumerator :: ep_Diag                       !% diagnostic elements (static)
        enumerator :: ep_ETM                        !% all ETM elements
        enumerator :: ep_JM                         !% all JM elements
        enumerator :: ep_JM_AC                      !% junction mains using AC method
        enumerator :: ep_JM_ALLtm                   !% Junction mains with any time march (static)
        enumerator :: ep_JM_ETM                     !% junction mains using ETM method
        !enumerator :: ep_JB_AC                      !% junction branches using AC method
        !enumerator :: ep_JB_ALLtm                   !% Junction branches with any time march (static)
        !enumerator :: ep_JB_ETM                     !% junction branches using ETM method
        !enumerator :: ep_AC_ACnonSurcharged         !% all surcharged with AC
        !enumerator :: ep_ALLtm_NonSurcharged        !% all time march nonsurcharged (neither PS nor AC surcharge)
        !enumerator :: ep_ETM_PSnonSurcharged        !% all surcharged with ETM
        !enumerator :: ep_SmallDepth_CC_ALLtm_posSlope !% small depth conduit cells with any time march and positive bottom slope
        !enumerator :: ep_SmallDepth_CC_ALLtm_negSlope !% small depth conduit cells with any time march and negative (adverse) bottom slope
        !enumerator :: ep_SmallDepth_CC_ALLtm
        enumerator :: ep_SmallDepth_CC_ETM
        !enumerator :: ep_SmallDepth_CC_AC
        !enumerator :: ep_SmallDepth_JM_ALLtm  !% 20220122brh
        enumerator :: ep_SmallDepth_JM_ETM    !% 20220122brh
        !enumerator :: ep_SmallDepth_JM_AC     !% 20220122brh
        !enumerator :: ep_ZeroDepth_CC_ALLtm         !% zero depth with any time march
        enumerator :: ep_ZeroDepth_CC_ETM
        !enumerator :: ep_ZeroDepth_CC_AC
        !enumerator :: ep_ZeroDepth_JM_ALLtm         !% zero depth JM
        enumerator :: ep_ZeroDepth_JM_ETM
        !enumerator :: ep_ZeroDepth_JM_AC
        !enumerator :: ep_ACsurcharged              !% all surcharged with AC
        !enumerator :: ep_ALLtmSurcharged           !% all time march surcharged
        !enumerator :: ep_PSsurcharged             !% all surcharged with ETM
        !enumerator :: ep_CCJM_H_ACsurcharged       !% all CCJM surcharged for H and AC solution
        !enumerator :: ep_CCJM_H_AC                  !% all CCJM solved for head with AC
        !enumerator :: ep_CCJB_eAC_i_fETM            !% all AC next to ETM
        enumerator :: ep_BClat                      !% all elements with lateral BC
        enumerator :: ep_JB_Downstream            !% all the downstream JB elements 
        enumerator :: ep_CC_DownstreamJBadjacent    !% all CC element downstream of a JB 
        enumerator :: ep_CC_Open_Elements           !% all CC elements that are open channel
        enumerator :: ep_CC_Closed_Elements         !% all closed CC elements 
        !enumerator :: ep_JM_Closed_Elements         !% all closed JM elements
        enumerator :: ep_JB_Closed_Elements         !% all closed JB elements 
        enumerator :: ep_JB_Open_Elements           !% all open-channel JB elements  
        enumerator :: ep_Output_Elements            !% all output elements -- local index   
        enumerator :: ep_CC_NOTsmalldepth           !% all Conduits that have time-marching without small or zero depth
        enumerator :: ep_CC_NOTzerodepth            !% all Conduits that have time-marching and are above zero depth
        enumerator :: ep_JBJM_NOTsmalldepth         !% all JB JM elements used in CFL computation 
        enumerator :: ep_CCJBJM_NOTsmalldepth       !% all elements used in CFL computation
        enumerator :: ep_CCJM_NOTsmalldepth         !% alternate elements for CFL computation 
        enumerator :: ep_CC_Transect                !% all channel elements with irregular transect
        enumerator :: ep_FM_HW_all                  !% all Hazen-Williams Force Main elements
        enumerator :: ep_FM_HW_PSsurcharged      !% all Hazen-Williams Force Main elements Preissmann Slot method that are surcharged
        !enumerator :: ep_FM_HW_PS_NonSurcharged     !% all Hazen-Williams Force Main elements with Preissmann Slot that are not surcharged
        enumerator :: ep_FM_dw_PSsurcharged      !% all Darcy-Weisbach Force Main elements Preissmann Slot method that are surcharged
        enumerator :: ep_FM_dw_PSnonSurcharged     !% all Darcy-Weisbach Force Main elements with Preissmann Slot that are not surcharged
        enumerator :: ep_lastplusone !% must be last enum item
    end enum
    integer, target :: Ncol_elemP = ep_lastplusone-1

!%
!%==========================================================================
!% PACKED GEOMETRY ELEMENTS
!%==========================================================================
!%      
    !%-------------------------------------------------------------------------
    !% Define the column indexes for elemPGalltm(:,:), elemPGetm(:,:),
    !% and elemPGac(:,:) arrays
    !% These are the packed arrays of geometry
    !% Note the same columns are used in the three different arrays
    !%-------------------------------------------------------------------------

    enum, bind(c)
        !% --- open channels
        enumerator :: epg_CC_rectangular = 1          !% CC rectangular channels 
        enumerator :: epg_CC_trapezoidal              !% CC trapezoidal channels 
        enumerator :: epg_CC_triangular               !% CC triangular channels 
        enumerator :: epg_CC_parabolic                !% CC parabolic channels 
        enumerator :: epg_CC_power_function           !% CC power function channels
        enumerator :: epg_CC_irregular                !% CC irregular channels 
        !% --- closed conduits
        enumerator :: epg_CC_rectangular_closed       !% CC rectangular conduits
        enumerator :: epg_CC_rectangular_round        !% CC rectangular-round conduits
        enumerator :: epg_CC_rectangular_triangular   !% CC rectangular-triangular conduit
        enumerator :: epg_CC_circular                 !% CC circular conduits 
        enumerator :: epg_CC_semi_circular            !% CC semi-circular conduits
        enumerator :: epg_CC_filled_circular          !% CC filled circular conduits
        enumerator :: epg_CC_arch                     !% CC arch conduits
        enumerator :: epg_CC_basket_handle            !% CC basket-handle conduits
        enumerator :: epg_CC_catenary                 !% CC catenary conduits 
        enumerator :: epg_CC_egg_shaped               !% CC egg shaped conduits 
        enumerator :: epg_CC_gothic                   !% CC gothic conduits
        enumerator :: epg_CC_horse_shoe               !% CC horse shoe conduits 
        enumerator :: epg_CC_horiz_ellipse            !% CC horizontal ellipse conduits
        enumerator :: epg_CC_mod_basket               !% CC modified basked conduits
        enumerator :: epg_CC_semi_elliptical          !% CC semi-elliptical conduits
        enumerator :: epg_CC_vert_ellipse             !% CC vertical elliptical conduits
        !% --- junctions
        enumerator :: epg_JM_functionalStorage        !% JM functional geometry relationship 
        enumerator :: epg_JM_tabularStorage           !% JM tabular geometry relationship 
        enumerator :: epg_JM_impliedStorage           !% JM with artificial storage
        enumerator :: epg_lastplusone !% must be last enum item
    end enum
    integer, target :: Ncol_elemPGalltm =  epg_lastplusone-1
    integer, target :: Ncol_elemPGetm   =  epg_lastplusone-1
    integer, target :: Ncol_elemPGac    =  epg_lastplusone-1

!%
!%==========================================================================
!% SPECIAL FEATURE ELEMENTS
!%==========================================================================
!%  
    !%-------------------------------------------------------------------------
    !% Define the column indexes for elemSI(:,:) arrays
    !% These are the full arrays of special integer data
    !% NOTE: each esi type should ONLY write to rows for elements of the
    !% correct type, otherwise it could be ovewriting data for a
    !% different special element type. This idea is important when initializing
    !% esi values!
    !%-------------------------------------------------------------------------

    !% --- CULVERT
    enum, bind(c)
        !% define the column indexes for the elemSi(:,:)  culvert
        enumerator :: esi_Conduit_Culvert_Code  =1        !% culvert code number
        enumerator :: esi_Conduit_Culvert_EquationForm    !% first column of table H-2 in SWMM hydraulics manual
        enumerator :: esi_Conduit_Culvert_Part            !% type key for inlet, outlet in/out
        enumerator :: esi_Conduit_Culvert_OutletID        !% element # for outlet (only stored for inlet)
        enumerator :: esi_Conduit_ForceMain_Method        !% type key HazenWilliams or DarcyWeisbach
        enumerator :: esi_Conduit_lastplusone     !% must be last enum item
    end enum
    integer, parameter :: Ncol_elemSI_Conduit = esi_Conduit_lastplusone-1
    
!    ! % NOTE: storing force main data in esi or esr is
!     !% --- FORCE MAIN
!     enum, bind(c)
!         !% define the column indexes for the elemSi(:,:) force main elements
!         enumerator :: esi_Conduit_Forcemain_Method = 1       !% type key  HazenWilliams or DarcyWeisbach
!         !enumerator :: esi_ForceMain_isSubmerged    !% 0 = no, 1 = yes
!         enumerator :: esi_ForceMain_lastplusone    !% must be last enum item
!     end enum
!     integer, parameter :: Ncol_elemSI_ForceMain = esi_ForceMain_lastplusone-1


    !% --- JUNCTION
    enum, bind(c)
        !% define the column indexes for elemSI(:,:) junction branch elements
        !% Note that esi_JunctionMain, esi_JunctionBranch, and (if needed) esi_Storage will
        !% share the same column sets.
        enumerator ::  esi_JunctionMain_Type       = 1             !% KEY junction main type
        enumerator ::  esi_JunctionMain_Curve_ID                   !% id of the junction storage cure if exists
        enumerator ::  esi_JunctionBranch_Exists                   !% assigned 1 if branch exists
        enumerator ::  esi_JunctionBranch_Link_Connection          !% the link index connected to that junction branch
        enumerator ::  esi_JunctionBranch_Main_Index               !% elem idx of the junction main for this branch
        enumerator ::  esi_JunctionBranch_lastplusone !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is
    integer, parameter :: Ncol_elemSI_junction = esi_JunctionBranch_lastplusone-1

    enum, bind(c)
        !% define the column indexes for elemSi(:,:) weir elements
        enumerator :: esi_Weir_EndContractions = 1      !% number of endcontractions of the weir
        enumerator :: esi_Weir_FlowDirection            !% weir flow direction (-1, +1)
        enumerator :: esi_Weir_SpecificType             !% KEY specific weir type
        enumerator :: esi_Weir_GeometryType             !% KEY specific weir geometry type
        enumerator :: esi_Weir_RoadSurface              !% road surface type for roadway weir
        enumerator :: esi_Weir_lastplusone !% must be last enum item
    end enum

    integer, parameter :: Ncol_elemSI_weir = esi_Weir_lastplusone-1

    enum, bind(c)
        !% define the column indexes for elemSi(:,:) orifice elements
        enumerator :: esi_Orifice_FlowDirection = 1     !% orifice flow direction (-1, +1)
        enumerator :: esi_Orifice_SpecificType          !% KEY specific orifice type
        enumerator :: esi_Orifice_GeometryType          !% KEY specific orifice geometry type
        enumerator :: esi_Orifice_lastplusone !% must be last enum item
    end enum
    integer, parameter :: Ncol_elemSI_orifice = esi_Orifice_lastplusone-1


    enum, bind(c)
        !% define the column indexes for elemSi(:,:) outlet elements
        enumerator :: esi_Outlet_FlowDirection = 1     !% outlet flow direction (-1, +1)
        enumerator :: esi_Outlet_SpecificType          !% KEY specific outlet type
        enumerator :: esi_Outlet_CurveID               !% outlet curve id
        enumerator :: esi_Outlet_hasFlapGate           !% 1 if true, 0 if false
        enumerator :: esi_Outlet_lastplusone !% must be last enum item
    end enum
    integer, parameter :: Ncol_elemSI_outlet = esi_Orifice_lastplusone-1

    enum, bind(c)
        !% define the column indexes for elemSi(:,:) outlet elements
        enumerator :: esi_Pump_FlowDirection = 1     !% pump flow direction (-1, +1)
        enumerator :: esi_Pump_SpecificType          !% KEY specific pump type
        enumerator :: esi_Pump_CurveID               !% pump curve id
        enumerator :: esi_Pump_IsControlled          !% 1 for external control, 0 for upstream control
        !enumerator :: esi_Pump_Status                !% 1 = on, 0 =off !% 20220625brh removed -- use er_Setting
        enumerator :: esi_Pump_lastplusone !% must be last enum item
    end enum
    integer, parameter :: Ncol_elemSI_Pump = esi_Pump_lastplusone-1


    !% determine the largest number of columns for a special set
    integer, target :: Ncol_elemSI = max(&
                            Ncol_elemSI_Conduit,    &
                            Ncol_elemSI_Junction, &
                            Ncol_elemSI_Orifice, &
                            Ncol_elemSI_Outlet, &
                            Ncol_elemSI_Pump,  &
                            Ncol_elemSI_Weir &                            
                            )

    !%-------------------------------------------------------------------------
    !% Define the column indexes for elemSr(:,:) arrays
    !% These are the full arrays if special real data
    !% Note that different types of special elements (diagnostic, branches)
    !% share the same columns since a row can only have one type of element.
    !%-------------------------------------------------------------------------

    ! !% define the column indexes for elemSr(:,:) for geometry that has not yet been confirmed and assigned:
    ! enum, bind(c)
    !     enumerator ::  esr_JunctionBranch_lastplusone = 1 !% must be last enum item
    ! end enum

    ! integer, parameter :: Ncol_elemSR_JunctionBranch = esr_JunctionBranch_lastplusone-1

    !% define the column indexes for elemSR(:,:) for geometry that has not yet been confirmed and assigned:
    !% Note that esr_JunctionMain, esr_JunctionBranch and esr_Storage share the same column sets.

    !% --- Conduit data 
    enum, bind(c)
        enumerator :: esr_Conduit_Culvert_K =1      !% see Appendix H-1 in SWMM5 Hydraulics                
        enumerator :: esr_Conduit_Culvert_M 
        enumerator :: esr_Conduit_Culvert_C
        enumerator :: esr_Conduit_Culvert_Y
        enumerator :: esr_Conduit_Culvert_SCF
        enumerator :: esr_Conduit_ForceMain_Coef                  !% Hazen-Williams C or Darcy-Weisbach epsilon
        enumerator :: esr_Conduit_ForceMain_FrictionFactor         !% Darcy-Weisbach friction factor
        enumerator :: esr_Conduit_lastplusone             !% must be last enum item
    end enum
    integer, parameter :: Ncol_elemSR_Conduit = esr_Conduit_lastplusone-1

    !% --- FORCE MAIN 
    !%     This implies that Storage, Pump, Weir, Orifice, Outlet cannot also be Force Main
    !%     HACK -- if Storage needs to be defined as force main, then the coef will need to
    !%     be moved to the elemR array.
    ! enum, bind(c)
    !     enumerator :: esr_Conduit_ForceMain_Coef = 1               !% Hazen-Williams C or Darcy-Weisbach epsilon
    !     enumerator :: esr_ForceMain_FrictionFactor         !% Darcy-Weisbach friction factor
    !     enumerator :: esr_ForceMain_lastplusone            !% must be last enum item
    ! end enum
    ! integer, parameter :: Ncol_elemSR_ForceMain = esr_ForceMain_lastplusone-1
   
    !% --- JUNCTION MAIN and STORAGE
    enum, bind(c)
        enumerator ::  esr_JunctionMain_PondedArea = 1
        !enumerator ::  esr_JunctionMain_PondedVolume
        enumerator ::  esr_JunctionMain_SurchargeExtraDepth
        enumerator ::  esr_JunctionBranch_Kfactor
        enumerator ::  esr_Storage_Constant
        enumerator ::  esr_Storage_Coefficient
        enumerator ::  esr_Storage_Exponent
        enumerator ::  esr_Storage_Plane_Area
        enumerator ::  esr_Storage_FractionEvap   !% --- fraction of area over which evaporation acts
        enumerator ::  esr_Storage_lastplusone !% must be last enum item
    end enum
    integer, parameter :: Ncol_elemSR_Junction = esr_Storage_lastplusone-1

    !% --- ORIFICE
    enum, bind(c)
        enumerator ::  esr_Orifice_CriticalDepth = 1        !% critical depth bellow which the orifice acts like an weir
        enumerator ::  esr_Orifice_CriticalHead             !% critical head for weir flow through an orifice
        enumerator ::  esr_Orifice_FractionCriticalDepth    !% critical depth fracttion to distinct between weir and orifice flow
        enumerator ::  esr_Orifice_DischargeCoeff           !% discharge coefficient orifice
        enumerator ::  esr_Orifice_FullDepth                !% original orifice opening
        enumerator ::  esr_Orifice_FullArea                 !% original orifice opening area
        enumerator ::  esr_Orifice_EffectiveFullDepth       !% effective full depth after control intervention
        enumerator ::  esr_Orifice_EffectiveFullArea        !% effective full depth after control intervention
        enumerator ::  esr_Orifice_EffectiveHeadDelta       !% effective head delta across orifice
        enumerator ::  esr_Orifice_NominalDownstreamHead    !% nominal downstream head for orifice
        enumerator ::  esr_Orifice_Orate                    !% orifice time to operate (close the gate)
        enumerator ::  esr_Orifice_RectangularBreadth       !% rectangular orifice breadth
        enumerator ::  esr_Orifice_Zcrown                   !% orifice "crown" elevation - highest edge of orifice
        enumerator ::  esr_Orifice_Zcrest                   !% orifice "crest" elevation - lowest edge of orifice
        enumerator ::  esr_Orifice_lastplusone !% must be last enum item
    end enum
    integer, parameter :: Ncol_elemSR_Orifice = esr_Orifice_lastplusone-1

    !% --- OUTLET
    enum, bind(c)
        enumerator ::  esr_Outlet_DischargeCoeff = 1       !% discharge coefficient outlet
        enumerator ::  esr_Outlet_EffectiveHeadDelta       !% effective head delta across outlet
        enumerator ::  esr_Outlet_NominalDownstreamHead    !% nominal downstream head for outlet
        enumerator ::  esr_Outlet_Exponent                 !% exponent for outlet dishcharge relation
        enumerator ::  esr_Outlet_Coefficient              !% power for outlet dishcharge relation
        enumerator ::  esr_Outlet_Zcrest                   !% outlet "crest" elevation - lowest edge of outlet
        enumerator ::  esr_Outlet_lastplusone !% must be last enum item
    end enum
    integer, parameter :: Ncol_elemSR_Outlet = esr_Outlet_lastplusone-1

    !% --- PUMP
    enum, bind(c)
        enumerator ::  esr_Pump_EffectiveHeadDelta = 1     !% effective head delta across outlet
        enumerator ::  esr_Pump_NominalDownstreamHead      !% nominal downstream head for outlet
        enumerator ::  esr_Pump_yOn                        !% pump startup depth
        enumerator ::  esr_Pump_yOff                       !% pump shutoff depth
        enumerator ::  esr_Pump_xMin                       !% minimum pt. on pump curve 
        enumerator ::  esr_Pump_xMax                       !% maximum pt. on pump curve
        enumerator ::  esr_Pump_lastplusone                !% must be last enum item
    end enum
    integer, parameter :: Ncol_elemSR_Pump = esr_Pump_lastplusone-1

    !% --- WEIR
    enum, bind(c)
        enumerator ::  esr_Weir_Rectangular = 1         !% discharge coefficient for the rectangular portion
        enumerator ::  esr_Weir_Triangular              !% discharge coefficient for triangular weir part
        enumerator ::  esr_Weir_FullDepth               !% original weir opening
        enumerator ::  esr_Weir_FullArea                !% original weir opening area
        enumerator ::  esr_Weir_EffectiveFullDepth      !% effective full depth after control intervention
        enumerator ::  esr_Weir_EffectiveHeadDelta      !% effective head delta across weir
        enumerator ::  esr_Weir_NominalDownstreamHead   !% nominal downstream head
        enumerator ::  esr_Weir_RectangularBreadth      !% rectangular weir breadth
        enumerator ::  esr_Weir_TrapezoidalBreadth      !% trapezoidal weir breadth
        enumerator ::  esr_Weir_TrapezoidalLeftSlope    !% trapezoidal weir left slope
        enumerator ::  esr_Weir_TrapezoidalRightSlope   !% trapezoidal weir right slope
        enumerator ::  esr_Weir_TriangularSideSlope     !% triangular weir side slope
        enumerator ::  esr_Wier_RoadWidth               !% road width for roadway weir
        enumerator ::  esr_Weir_Zcrown                  !% weir crown elevation
        enumerator ::  esr_Weir_Zcrest                  !% weir crest elevation
        enumerator ::  esr_Weir_lastplusone !% must be last enum item
    end enum
    integer, parameter :: Ncol_elemSR_Weir = esr_Weir_lastplusone-1

    

    !% determine the largest number of columns for a special set
    integer, target :: Ncol_elemSR = max(&
                            Ncol_elemSR_Conduit,        &
                            Ncol_elemSR_Junction,        &
                            Ncol_elemSR_Orifice,        &
                            Ncol_elemSR_Outlet,         &
                            Ncol_elemSR_Pump,           &
                            Ncol_elemSR_Weir           &
                            ) !, &
                            ! Ncol_elemSR_Conduit)

    !% HACK: Ncol_elemSR must be updated when other special elements
    !% (i.e. orifice, pump, storage etc.) are added

!%
!%==========================================================================
!% SPECIAL GEOMETRY ELEMENTS
!%==========================================================================
!%                              
    !%-------------------------------------------------------------------------
    !% Define the column indexes for the elemSGR(:,:) arrays
    !% These are the full arrays of special, geometry, real data
    !%-------------------------------------------------------------------------

    !% ===== OPEN CHANNELS ======

    !% Define the column indexes for elemGSR(:,:) for irregular channel
    enum, bind(c)
        enumerator ::  esgr_Irregular_lastplusone = 1    !% dummy -- no data at this point
    end enum   
    integer, parameter :: Ncol_elemSGR_Irregular =  esgr_Irregular_lastplusone-1  

    !% Define the column indexes for elemGSR(:,:) for parabolic channel
    enum, bind(c)
        enumerator ::  esgr_Parabolic_Breadth = 1    !% breadth for parabolic geometry
        enumerator ::  esgr_Parabolic_Radius
        enumerator ::  esgr_Parabolic_lastplusone !% must be last enum item
    end enum    
    integer, parameter :: Ncol_elemSGR_Parabolic =  esgr_Parabolic_lastplusone-1    
    
    !% Define the column indexes for elemGSR(:,:) for powerfunction channel
    enum, bind(c)
        enumerator ::  esgr_PowerFunction_Breadth = 1    
        enumerator ::  esgr_PowerFunction_Radius
        enumerator ::  esgr_PowerFunction_lastplusone !% must be last enum item
    end enum    
    integer, parameter :: Ncol_elemSGR_PowerFunction =  esgr_PowerFunction_lastplusone-1     

    !% Define the column indexes for elemGSR(:,:) for rectangular pipe or channel
    enum, bind(c)
        enumerator ::  esgr_Rectangular_Breadth = 1    !% breadth for rectangular geometry
        enumerator ::  esgr_Rectangular_lastplusone !% must be last enum item
    end enum
    integer, parameter :: Ncol_elemSGR_Rectangular =  esgr_Rectangular_lastplusone-1

    !% Define the column indexes for elemGSR(:,:) for triangular channel
    enum, bind(c)
        enumerator ::  esgr_Triangular_Slope = 1    !% side-slope for triangular geometry
        enumerator ::  esgr_Triangular_TopBreadth  !% top breadth of triangular geometry
        enumerator ::  esgr_Triangular_lastplusone !% must be last enum item
    end enum
    integer, parameter :: Ncol_elemSGR_Triangular =  esgr_Triangular_lastplusone-1

    !% Define the column indexes for elemGSR(:,:) for trapezoidal pipe or channel
    enum, bind(c)
        enumerator ::  esgr_Trapezoidal_Breadth = 1    !% bottom breadth for trapezoidal geometry
        enumerator ::  esgr_Trapezoidal_LeftSlope      !% left slope for trapezoidal geometry
        enumerator ::  esgr_Trapezoidal_RightSlope     !% right slope for trapezoidal geometry
        enumerator ::  esgr_Trapezoidal_lastplusone !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSGR_Trapezoidal =  esgr_Trapezoidal_lastplusone-1



    !% ===== CLOSED CONDUITS =====

    !% Define the column indexes for elemGSR(:,:) for Arch shaped conduits
    enum, bind(c)
        enumerator ::  esgr_Arch_SoverSfull=1       !% S/Sfull for basket handle geometry
        enumerator ::  esgr_Arch_lastplusone      !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSGR_Arch = esgr_Arch_lastplusone-1

    !% Define the column indexes for elemGSR(:,:) for basket_handle_conduit
    enum, bind(c)
        enumerator ::  esgr_Basket_Handle_lastplusone = 1      !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSGR_Basket_Handle = esgr_Basket_Handle_lastplusone-1

    !% Define the column indexes for elemGSR(:,:) for Catenary shaped conduits
    enum, bind(c)
        enumerator ::  esgr_Catenary_SoverSfull = 1       !% S/Sfull for basket handle geometry
        enumerator ::  esgr_Catenary_lastplusone      !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSGR_Catenary = esgr_Catenary_lastplusone-1
   
    !% Define the column indexes for elemGSR(:,:) for circular conduit 
    enum, bind(c)
        enumerator ::  esgr_Circular_Diameter =1         !% diameter for circular geometry
        enumerator ::  esgr_Circular_Radius            !% radius for circular geometry
        enumerator ::  esgr_Circular_lastplusone !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSGR_Circular =  esgr_Circular_lastplusone-1

    !% Define the column indexes for elemGSR(:,:) for custom closed conduit
    enum, bind(c)
        enumerator ::  esgr_Custom_lastplusone = 1 !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSGR_Custom =  esgr_Custom_lastplusone-1

    !% Define the column indexes for elemGSR(:,:) for Egg_Shaped_conduit
    enum, bind(c)
         enumerator ::  esgr_Egg_Shaped_lastplusone = 1      !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSGR_Egg_Shaped = esgr_Egg_Shaped_lastplusone-1

    !% Define the column indexes for elemGSR(:,:) for filled circular pipe or channel
    enum, bind(c)
        enumerator ::  esgr_Filled_Circular_TotalPipeDiameter = 1   !% diameter for entire pipe (filled + flow)
        enumerator ::  esgr_Filled_Circular_TotalPipeArea       !% pipe cross-section area including filled area
        enumerator ::  esgr_Filled_Circular_TotalPipePerimeter  !% pipe perimeter including filled area
        enumerator ::  esgr_Filled_Circular_TotalPipeHydRadius  !% pipe hydradius including filled area
        enumerator ::  esgr_Filled_Circular_bottomArea          !% filled area of filled circular geometry
        enumerator ::  esgr_Filled_Circular_bottomPerimeter     !% filled wetted perimeter of filled circular geometry
        enumerator ::  esgr_Filled_Circular_bottomTopwidth      !% filled top-width of filled circular geometry
        enumerator ::  esgr_Filled_Circular_lastplusone     !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSGR_Filled_Circular =  esgr_Filled_Circular_lastplusone-1

    !% Define the column indexes for elemGSR(:,:) for Gothic shaped conduits
    enum, bind(c)
        enumerator ::  esgr_Gothic_SoverSfull = 1       !% S/Sfull for basket handle geometry
        enumerator ::  esgr_Gothic_lastplusone      !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSGR_Gothic = esgr_Gothic_lastplusone-1

    !% Define the column indexes for elemGSR(:,:) for HOrizontal ellipse shaped conduits
    enum, bind(c)
        enumerator ::  esgr_Horiz_Ellipse_SoverSfull = 1       !% S/Sfull for basket handle geometry
        enumerator ::  esgr_Horiz_Ellipse_lastplusone      !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSGR_Horiz_Ellipse = esgr_Horiz_Ellipse_lastplusone-1

    !% Define the column indexes for elemGSR(:,:) for Horse Shoe shaped conduits
    enum, bind(c)
        enumerator ::  esgr_Horse_Shoe_lastplusone = 1     !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSGR_Horse_Shoe = esgr_Horse_Shoe_lastplusone-1

    !% Define the column indexes for elemGSR(:,:) for mod_basket conduit
    enum, bind(c)
        enumerator ::  esgr_Mod_Basket_Ytop =1            !% height of top circular arc
        enumerator ::  esgr_Mod_Basket_Rtop             !% radius of top circular arc
        enumerator ::  esgr_Mod_Basket_Atop             !% area of top circular arc
        enumerator ::  esgr_Mod_Basket_ThetaTop         !% angle of top circular arc
        enumerator ::  esgr_Mod_Basket_lastplusone      !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSGR_Mod_Basket = esgr_Mod_Basket_lastplusone-1

     !% Define the column indexes for elemGSR(:,:) for Rectangular closedchannel
    enum, bind(c)
        enumerator ::  esgr_Rectangular_Closed_lastplusone = 1       !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSGR_Rectangular_Closed =  esgr_Rectangular_Closed_lastplusone-1


    !% Define the column indexes for elemGSR(:,:) for Rectangular round channel
    enum, bind(c)
        enumerator ::  esgr_Rectangular_Round_Ybot = 1              !% depth of bottom circular section
        enumerator ::  esgr_Rectangular_Round_Rbot              !% radius of the circular section
        enumerator ::  esgr_Rectangular_Round_Abot              !% area of the circular section
        enumerator ::  esgr_Rectangular_Round_ThetaBot          !% angle of the circular section
        enumerator ::  esgr_Rectangular_Round_lastplusone       !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSGR_Rectangular_Round =  esgr_Rectangular_Round_lastplusone-1

    !% Define the column indexes for elemGSR(:,:) for rectangular triangular channel
    enum, bind(c)
        enumerator ::  esgr_Rectangular_Triangular_BottomDepth = 1     !% depth of the triangular section
        enumerator ::  esgr_Rectangular_Triangular_BottomArea      !% area of the triangular section
        enumerator ::  esgr_Rectangular_Triangular_BottomSlope     !% side slope of the triangular section
        enumerator ::  esgr_Rectangular_Triangular_lastplusone     !% must be last enum item
    end enum
    integer, parameter :: Ncol_elemSGR_Rectangular_Triangular =  esgr_Rectangular_Triangular_lastplusone-1

    !% Define the column indexes for elemGSR(:,:) for Semi-Circular shaped conduits
    enum, bind(c)
        enumerator ::  esgr_Semi_Circular_SoverSfull  = 1     !% S/Sfull for basket handle geometry
        enumerator ::  esgr_Semi_Circular_lastplusone      !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSGR_Semi_Circular = esgr_Semi_Circular_lastplusone-1

    !% Define the column indexes for elemGSR(:,:) for Semi-Circular shaped conduits
    enum, bind(c)
        enumerator ::  esgr_Semi_Elliptical_SoverSfull = 1       !% S/Sfull for basket handle geometry
        enumerator ::  esgr_Semi_Elliptical_lastplusone      !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSGR_Semi_Elliptical = esgr_Semi_Elliptical_lastplusone-1


    !% Define the column indexes for elemGSR(:,:) for HOrizontal ellipse shaped conduits
    enum, bind(c)
        enumerator ::  esgr_Vert_Ellipse_SoverSfull = 1       !% S/Sfull for basket handle geometry
        enumerator ::  esgr_Vert_Ellipse_lastplusone      !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSGR_Vert_Ellipse = esgr_Vert_Ellipse_lastplusone-1

    !% Define the column indexes for elemSGR(:,:) for other geometry

    !% NEED OTHER GEOMETRY HERE

    !% determine the largest number of columns for a special set
    integer, target :: Ncol_elemSGR = max(&
                            Ncol_elemSGR_Irregular,              &
                            Ncol_elemSGR_Parabolic,              &
                            Ncol_elemSGR_PowerFunction,          &
                            Ncol_elemSGR_Rectangular,            &
                            Ncol_elemSGR_Trapezoidal,            &
                            Ncol_elemSGR_Arch,                   &
                            Ncol_elemSGR_Basket_Handle,          &
                            Ncol_elemSGR_Catenary,               &
                            Ncol_elemSGR_Circular,               &
                            Ncol_elemSGR_Custom,                 &
                            Ncol_elemSGR_Egg_Shaped,             &
                            Ncol_elemSGR_Filled_Circular,        &
                            Ncol_elemSGR_Gothic,                 &
                            Ncol_elemSGR_Horiz_Ellipse,          &
                            Ncol_elemSGR_Horse_Shoe,             &
                            Ncol_elemSGR_Mod_Basket,             &
                            Ncol_elemSGR_Rectangular_Closed,     &
                            Ncol_elemSGR_Rectangular_Round,      &
                            Ncol_elemSGR_Rectangular_Triangular, &
                            Ncol_elemSGR_Semi_Circular,          &
                            Ncol_elemSGR_Semi_Elliptical,        &
                            Ncol_elemSGR_Vert_Ellipse)


    !% HACK: Ncol_elemSR must be updated when other geometry types
    !% (i.e. triangular, circular etc.) are added for channel or
    !% conduit elements
!%
!%==========================================================================
!% SYSTEM CONTROL MONITORING AND ACTION ELEMENTS
!%==========================================================================
!%     
    !%-------------------------------------------------------------------------
    !% Define the column indexes for monitorI(:,:) array that stores
    !% monitoring point data that might be on other images
    !%-------------------------------------------------------------------------
    enum, bind(c)
        enumerator :: mi_idx = 1       !% unique ID for this monitoring element
        enumerator :: mi_image         !% image where the M element resides
        enumerator :: mi_elem_idx      !% element index on the image of the M point
        enumerator :: mi_linknodesimType  !% = 1 if link, 0 if node, -1 if simulation
        enumerator :: mi_linknode_idx  !% EPA SWMM link or node index for monitoring point 
        enumerator :: mi_lastplusone   !% must be last enum item
    end enum  
    integer, target :: Ncol_MonitoringPointI = mi_lastplusone - 1
    
    !%-------------------------------------------------------------------------
    !% Define the column indexes for monitorR(:,:) array that stores
    !% monitoring data that might be on other images
    !% NOTE if any other columns are added, the col must be
    !% updated with the corresponding elemR columns in util_allocate_monitor()
    !%-------------------------------------------------------------------------

    !{r_DEPTH, r_HEAD, r_VOLUME, r_INFLOW, r_FLOW, r_STATUS,
    !            r_SETTING, r_TIMEOPEN, r_TIMECLOSED, r_TIME, r_DATE,
    !            r_CLOCKTIME, r_DAYOFYEAR, r_DAY, r_MONTH};
    enum, bind(c)
        enumerator :: mr_Depth = 1    !% depth on control/monitoring element
        enumerator :: mr_Head     
        enumerator :: mr_Volume
        enumerator :: mr_Inflow
        enumerator :: mr_Flow
        enumerator :: mr_Setting      !% pump on/off status and link setting
        enumerator :: mr_TimeLastSet  !% last time (seconds) the setting was changed
        enumerator :: mr_lastplusone  !% must be last enum item
    end enum
    integer, target :: Ncol_MonitoringPointR = mr_lastplusone -1


    !%-------------------------------------------------------------------------
    !% Define the column indexes for actionI(:,:) array that stores
    !% control action point data that might be on other images
    !%-------------------------------------------------------------------------
    enum, bind(c)
        enumerator :: ai_idx = 1       !% unique ID for this action element
        enumerator :: ai_image         !% image where the action element resides
        enumerator :: ai_elem_idx      !% element index on the image of the action point
        enumerator :: ai_link_idx      !% EPA SWMM link index for action point
        enumerator :: ai_hasChanged    !% 1 if setting was changed, 0 if not
        enumerator :: ai_lastplusone   !% must be last enum item
    end enum  
    integer, target :: Ncol_ActionPointI = ai_lastplusone - 1

    !%-------------------------------------------------------------------------
    !% Define the column indexes for actionR(:,:) array that stores
    !% monitoring data that might be on other images
    !% NOTE if any other columns are added, the col must be
    !% updated with the corresponding elemR columns in util_allocate_action()
    !%-------------------------------------------------------------------------
    enum, bind(c)
        enumerator :: ar_dummy = 1         !% NOT SURE IF WE WILL NEED THIS STRUCTURE
        enumerator :: ar_lastplusone
    end enum
    integer, target :: Ncol_ActionPointR = ar_lastplusone -1


    ! !%-------------------------------------------------------------------------
    ! !% Define the column indexes for conmonYN(:,:) array that stores
    ! !% control and monitoring data that might be on other images
    ! !%-------------------------------------------------------------------------
    ! enum, bind(c)
    !     enumerator :: mpYN_lastplusone !% must be last enum item
    ! end enum
    ! integer, target :: Ncol_MonitoringPointYN = mpYN_lastplusone - 1

!%
!%==========================================================================
!% INTER-IMAGE BOUNDARY/GHOST ELEMENTS
!%==========================================================================
!%     
    !%-------------------------------------------------------------------------
    !% Define the column indexes for elemB%I/elemGI(:,:) arrays
    !% These arrays are used to store/transfer inter image data
    !%-------------------------------------------------------------------------
    enum, bind(c)
         enumerator :: ebgi_idx = 1                   !% row index of elemB%I/elemGI array 
         enumerator :: ebgi_elem_Lidx                 !% local element index (static)
         enumerator :: ebgi_elem_Gidx                 !% global element index  (static)
         enumerator :: ebgi_Mface_uL                  !% map to upstream face local index  (static)
         enumerator :: ebgi_Mface_dL                  !% map to downstream face local index  (static)
         enumerator :: ebgi_lastplusone !% must be last enum item
    end enum
    integer, target :: Ncol_elemBGI = ebgi_lastplusone-1

    !%-------------------------------------------------------------------------
    !% Define the column indexes for elemB%R/elemGR(:,:) arrays
    !% These arrays are used to store/transfer inter image data
    !%-------------------------------------------------------------------------
    enum, bind(c)
        enumerator ::  ebgr_Area = 1                  !% cross-sectional flow area (latest) boundary/ghost element 
        enumerator ::  ebgr_Topwidth                  !% topwidth of flow at free surfac boundary/ghost element
        enumerator ::  ebgr_Depth                     !% Depth of flow boundary/ghost element
        enumerator ::  ebgr_Head                      !% piezometric head (latest) -- water surface elevation in open channel boundary/ghost element
        enumerator ::  ebgr_Flowrate                  !% flowrate (latest) boundary/ghost element
        enumerator ::  ebgr_Preissmann_Number         !% preissmann number boundary/ghost element
        enumerator ::  ebgr_Volume                    !% volume at boundary/ghost element
        enumerator ::  ebgr_InterpWeight_dG           !% interpolation Weight, downstream, for geometry boundary/ghost element
        enumerator ::  ebgr_InterpWeight_uG           !% interpolation Weight, upstream, for geometry boundary/ghost element 
        enumerator ::  ebgr_InterpWeight_dH           !% interpolation Weight, downstream for head boundary/ghost element
        enumerator ::  ebgr_InterpWeight_uH           !% interpolation Weight, upstream for head boundary/ghost element 
        enumerator ::  ebgr_InterpWeight_dQ           !% interpolation Weight, downstream, for flowrate boundary/ghost element
        enumerator ::  ebgr_InterpWeight_uQ           !% interpolation Weight, upstream, for flowrate boundary/ghost element
        enumerator ::  ebgr_InterpWeight_dP           !% interpolation Weight, downstream, for preissman number boundary/ghost element
        enumerator ::  ebgr_InterpWeight_uP           !% interpolation Weight, upstream, for preissman number boundary/ghost element
        enumerator ::  ebgr_lastplusone               !% must be last enum item boundary/ghost element
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, target :: Ncol_elemBGR =  ebgr_lastplusone-1


!%
!%==========================================================================
!% FACES
!%==========================================================================
!%  
    !%-------------------------------------------------------------------------
    !% Define the column indexes for faceI(:,:) arrays
    !% These are the full arrays of face integer data
    !%-------------------------------------------------------------------------

    enum, bind(c)
        enumerator ::  fi_Lidx = 1                  !% local array index (row)
        enumerator ::  fi_Gidx                      !% global (unique) index
        enumerator ::  fi_BCtype                    !% KEY type of BC on face
        enumerator ::  fi_barrels                   !% number of barrels for the face
        enumerator ::  fi_jump_type                 !% KEY Type of hydraulic jump
        enumerator ::  fi_Melem_uL                  !% map to element upstream (local index)
        enumerator ::  fi_Melem_dL                  !% map to element downstream (local index)
        enumerator ::  fi_GhostElem_uL              !% map to upstream ghost element
        enumerator ::  fi_GhostElem_dL              !% map to downstream ghost element
        enumerator ::  fi_BoundaryElem_uL           !% map to upstream boundary/ghost element in the boundary/ghost array
        enumerator ::  fi_BoundaryElem_dL           !% map to dwonstream boundary/ghost element in the boundary/ghost array
        enumerator ::  fi_Connected_image           !% image number a shared face connected to
        enumerator ::  fi_node_idx_BIPquick         !% if the face is originated from a node, then the BQ idx
        enumerator ::  fi_link_idx_BIPquick         !% face connected to a BQ link element 
        enumerator ::  fi_node_idx_SWMM             !% if the face is originated from a node, then the SWMM idx
        enumerator ::  fi_link_idx_SWMM             !% face connected to a SWMM link element 
        enumerator ::  fi_zeroDepth                 !% 0 is no zerodepth, -1 is zerodepth upstream, +1 is zerodepth downstream, +2 is both
        !% HACK: THESE MIGHT NEED TO BE RESTORED
        ! enumerator ::  fi_Melem_uG                 !% map to element upstream (global index)
        ! enumerator ::  fi_Melem_dG                 !% map to element upstream (global index)
        ! enumerator ::  fi_eHeqType_u               !% type of H solution on element upstream
        ! enumerator ::  fi_eHeqType_d               !% type of H solution on element downstream
        ! enumerator ::  fi_eQeqType_u               !% type of Q solution on element upstream
        ! enumerator ::  fi_eQeqType_d               !% type of Q solution on element downstream
        enumerator :: fi_lastplusone !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, target :: Ncol_faceI =  fi_lastplusone-1

    !%-------------------------------------------------------------------------
    !% Define the column indexes for faceR(:,:) arrays
    !% These are the full arrays of face mapping data
    !%-------------------------------------------------------------------------
    enum, bind(c)
        enumerator :: fr_Area_d = 1             !% cross-sectional area on downstream side of face
        enumerator :: fr_Area_u                 !% cross-sectional area on upstream side of face
        enumerator :: fr_Depth_d 
        enumerator :: fr_Depth_u
        enumerator :: fr_Flowrate               !% flowrate through face (latest)
        enumerator :: fr_Flowrate_N0            !% flowrate through face (time N)    enumerator :: fr_Head_d  !% Piezometric head on downstream side of face
        enumerator :: fr_Flowrate_Conservative  !% the effective flow rate over the time step N to N+1
        enumerator :: fr_FlowrateMax            !% maximum flowrate based on upstream volume/timestep
        enumerator :: fr_FlowrateMin            !% minimum flowrate (negative maximum) based on downstream volume/timestep
        !enumerator :: fr_Flowrate_Max           !% maximum flowrate based on upstream element
        enumerator :: fr_Head_u                 !% piezometric head on upstream side of face
        enumerator :: fr_Head_d                 !% piezometric head on downstream side of face
        enumerator :: fr_Zbottom                !% zbottom of faces
        !enumerator :: fr_ZbreadthMax            !% elevation of maximum breadth
        !enumerator :: fr_HydDepth_d             !% hydraulic Depth on downstream side of face
        !3numerator :: fr_HydDepth_u             !% hydraulic Depth on upstream side of face
        !enumerator :: fr_EllDepth_d             !% modified hydraulic Depth on downstream side of face
        !enumerator :: fr_EllDepth_u             !% modified hydraulic Depth on upstream side of face
        !enumerator :: fr_Topwidth_d             !% topwidth on downstream side of face
        !enumerator :: fr_Topwidth_u             !% topwidth on upstream side of face
        enumerator :: fr_Velocity_d             !% velocity on downstream side of face
        enumerator :: fr_Velocity_u             !% velocity on upstream side of face
        enumerator :: fr_Preissmann_Number      !% preissmann number at face

        !% HACK: THE FOLLOWING MAY NEED TO BE RESTORED
        ! enumerator :: fr_Zbottom_u             !% Bottom elevation on upstream side of face
        ! enumerator :: fr_Zbottom_d             !% Bottom elevation on downstream side of face
        ! enumerator :: fr_X                     !% Linear X location in system
        enumerator :: fr_lastplusone !% must be last enum item
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, target :: Ncol_faceR =  fr_lastplusone-1

    !%-------------------------------------------------------------------------
    !% Define the column indexes for faceYN(:,:) arrays
    !% These are the full arrays of face logical data
    !%-------------------------------------------------------------------------
    enum, bind(c)
        enumerator :: fYN_isAC_adjacent = 1
        enumerator :: fYN_isInteriorFace
        enumerator :: fYN_isSharedFace
        enumerator :: fYN_isUpGhost
        enumerator :: fYN_isDnGhost
        enumerator :: fYN_isnull
        enumerator :: fYN_isPSsurcharged
        enumerator :: fYN_isDownstreamJbFace
        enumerator :: fYN_isFaceOut
        !% HACK: The following might not be needed
        enumerator :: fYN_isDiag_adjacent
        ! enumerator :: fYN_isETM_adjacent
        ! enumerator :: fYN_isBCface
        enumerator :: fYN_lastplusone !% must be last enum item
    end enum
    integer, target :: Ncol_faceYN =  fYN_lastplusone-1

!%
!%==========================================================================
!% PACKED FACES
!%==========================================================================
!%  
    !%-------------------------------------------------------------------------
    !% Define the column indexes for faceP(:,:) and facePS(:,:) arrays
    !% These are for the packed array of face data
    !%-------------------------------------------------------------------------
    enum, bind(c)
        enumerator :: fp_all = 1                !% all faces execpt boundary, null, and shared faces
        enumerator :: fp_AC                     !% face with adjacent AC element
        enumerator :: fp_Diag                   !% face with adjacent diagnostic element
        enumerator :: fp_elem_downstream_is_zero
        enumerator :: fp_elem_upstream_is_zero
        enumerator :: fp_elem_bothsides_are_zero
        enumerator :: fp_JumpDn                 !% face with hydraulic jump from nominal downstream to upstream
        enumerator :: fp_JumpUp                 !% face with hydraulic jump from nominal upstream to downstream
        enumerator :: fp_BCup
        enumerator :: fp_BCdn
        enumerator :: fp_J1                     !% faces that are dead-ends of link without inflow BC
        enumerator :: fp_J1_BCup                !% faces that are either J1 or BCup
        enumerator :: fp_Output_Faces           !% faces that are selected for output
        enumerator :: fp_lastplusone !% must be last enum item
    end enum
    integer, target :: Ncol_faceP =  fp_lastplusone-1

!%
!%==========================================================================
!% PROFILER
!%==========================================================================
!% 
    !% row indexes for profiler data
    enum, bind(c)
        enumerator :: pfr_thisstart = 1
        enumerator :: pfr_thisend
        enumerator :: pfr_cumulative
        enumerator :: pfr_lastplusone !% must be last enum item
    end enum
    integer, target :: Nrow_pf = pfr_lastplusone-1

    !% column indexes for profiler_data. Each is the name of a procedure preceded by pfc_
    enum, bind(c)
        enumerator :: pfc_initialize_all = 1
        enumerator :: pfc_init_partitioning
        enumerator :: pfc_init_BIPquick
        enumerator :: pfc_init_network_define_toplevel
        enumerator :: pfc_init_bc
        enumerator :: pfc_init_IC_setup
        enumerator :: pfc_init_IC_from_linkdata
        enumerator :: pfc_init_IC_get_depth_from_linkdata
        enumerator :: pfc_init_IC_get_flow_roughness_from_linkdata
        enumerator :: pfc_init_IC_get_elemtype_from_linkdata
        enumerator :: pfc_init_IC_get_geometry_from_linkdata
        enumerator :: pfc_init_IC_get_channel_geometry
        enumerator :: pfc_init_IC_get_conduit_geometry
        enumerator :: pfc_init_IC_get_weir_geometry
        enumerator :: pfc_init_IC_get_orifice_geometry
        enumerator :: pfc_init_IC_get_channel_conduit_velocity
        enumerator :: pfc_init_IC_from_nodedata
        enumerator :: pfc_init_IC_get_junction_data
        enumerator :: pfc_geo_assign_JB
        enumerator :: pfc_update_auxiliary_variables
        enumerator :: pfc_init_IC_set_SmallVolumes
        enumerator :: pfc_init_IC_diagnostic_interpolation_weights
        enumerator :: pfc_face_interpolation
        enumerator :: pfc_diagnostic_toplevel
        enumerator :: pfc_lastplusone  !% must be last enum item
    end enum
    integer, target :: Ncol_pf = pfc_lastplusone-1

!%
!%==========================================================================
!% OUTPUT ARRAYS
!%==========================================================================
!%  
    !% data types (columns) in the OutElemFixedI
    enum, bind(c)
        enumerator :: oefi_elem_Gidx = 1
        enumerator :: oefi_link_Gidx_SWMM
        enumerator :: oefi_node_Gidx_SWMM
        enumerator :: oefi_lastplusone !% must be the last enum item
    end enum
    integer, target :: Ncol_oefi = oefi_lastplusone-1

    !% data types (columns) in the OutFaceFixedI
    enum, bind(c)
        enumerator :: offi_face_Gidx = 1
        enumerator :: offi_node_Gidx_SWMM
        enumerator :: offi_lastplusone !% must be the last enum item
    end enum
    integer, target :: Ncol_offi = offi_lastplusone-1

!%
!%==========================================================================
!% SUBCATCHMENTS
!%==========================================================================
!%  
    !%-------------------------------------------------------------------------
    !% Define the column indexes for subcatchR(:,:) arrays
    !% These are arrays with setting%SWMMinput%N_subcatch rows.
    !%-------------------------------------------------------------------------
    enum, bind(c)
        enumerator :: sr_RunOffRate_baseline=1    !% hydrology step runoff rate from EPASWMM
        enumerator :: sr_lastplusone            !% must be the last enum item
    end enum
    integer, target :: Ncol_subcatchR = sr_lastplusone-1

    !% --- additional sr columns that are of length N_subcatch_runon
    integer, allocatable, target :: sr_RunOn_Volume(:)

    !%-------------------------------------------------------------------------
    !% Define the column indexes for subcatchI(:,:) arrays
    !% These are arrays with setting%SWMMinput%N_subcatch rows.
    !%-------------------------------------------------------------------------
    enum, bind(c)
        enumerator :: si_runoff_nodeIdx = 1   !% node index for runoff
        enumerator :: si_runoff_elemIdx       !% runoff element for this subcatchment  
        enumerator :: si_runoff_P_image       !% coarray image that holds element for runoff
        enumerator :: si_RunOn_count          !% Number of runons to this subcatchment
        enumerator :: si_lastplusone          !% must be the last enum item
    end enum
    integer, target :: Ncol_subcatchI = si_lastplusone-1


    !% --- additional si columns that are of length N_subcatch_runon
    integer, allocatable, target :: si_RunOn_nodeIdx(:)
    integer, allocatable, target :: si_RunOn_SWMMoutfallIdx(:)
    integer, allocatable, target :: si_RunOn_faceIdx(:)
    integer, allocatable, target :: si_RunOn_P_image(:)

    
    !%-------------------------------------------------------------------------
    !% Define the column indexes for subcatchYN(:,:) arrays
    !% These are arrays with setting%SWMMinput%N_subcatch rows.
    !%-------------------------------------------------------------------------
    enum, bind(c)
        enumerator :: sYN_hasRunOff            !% TRUE if subcatchment has runoff to an element
        enumerator :: sYN_hasRunOn             !% TRUE if at least one runon to subcatchment
        enumerator :: sYN_lastplusone          !% must be the last enum item
    end enum
    integer, target :: Ncol_subcatchYN = sYN_lastplusone-1

!%
!%==========================================================================
!% CURVE DATA
!%==========================================================================
!%  
    !% data types (columns) in the table
    !% define column indexes for storage curve types
    enum, bind(c)
        enumerator :: curve_storage_depth = 1
        enumerator :: curve_storage_area
        enumerator :: curve_storage_volume
        enumerator :: curve_storage_lastplusone !% must be the last enum item
    end enum
    integer, parameter :: Ncol_storage_curve = curve_storage_lastplusone-1

    !% define column indexes for pump curve types
    enum, bind(c)
        enumerator :: curve_pump_Xvar = 1    !% this can be volume, head or depth used for interpolation
        enumerator :: curve_pump_flowrate
        enumerator :: curve_pump_lastplusone !% must be the last enum item
    end enum

    integer, parameter :: Ncol_pump_curve = curve_pump_lastplusone-1
    
    !% define column indexes for outlet rating curve types
    enum, bind(c)
        enumerator :: curve_outlet_depth = 1
        enumerator :: curve_outlet_flowrate
        enumerator :: curve_outlet_lastplusone !% must be the last enum item
    end enum

    integer, parameter :: Ncol_outlet_curve = curve_outlet_lastplusone-1

    !% determine the largest number of columns for table data structure
    integer, target :: Ncol_curve = max(&
                            Ncol_storage_curve, &
                            Ncol_pump_curve, &
                            Ncol_outlet_curve)

!%
!%==========================================================================
!% TRANSECT ARRAYS
!%==========================================================================
!%  
    !% transect integer array indexes
    enum, bind(c)
        enumerator :: ti_idx = 1
        enumerator :: ti_lastplusone 
    end enum    

    integer, parameter :: Ncol_transectI = ti_lastplusone-1

    !% transect real array indexes
    enum, bind(c)
        enumerator :: tr_depthFull = 1     ! depth when full (yFull in EPA-SWMM)
        enumerator :: tr_areaFull          ! area when full (aFull in EPA_SWMM)
        enumerator :: tr_hydRadiusFull     ! hydradius when full (hFull in EPA-SWMM)
        enumerator :: tr_widthMax          ! max width (wMax in EPA-SWMM)
        enumerator :: tr_depthAtBreadthMax ! depth at max width (ywMax in EPA-SWMM)
        enumerator :: tr_sectionFactor     ! section factor at max flow (sMax in EPA-SWMM)
        enumerator :: tr_areaAtMaxFlow     ! area at max flow (aMax in EPA-SWMM)
        enumerator :: tr_lengthFactor      ! floodplain / channel length (lengthFactor in EPA-SWMM)
        enumerator :: tr_ManningsN         ! Mannings n (roughness in EPA-SWMM)
        enumerator :: tr_widthFull         ! not in EPA-SWMM - width at full depth
        enumerator :: tr_areaBelowBreadthMax ! not in EPA-SWMM
        enumerator :: tr_lastplusone
    end enum

    integer, parameter :: Ncol_transectR = tr_lastplusone-1

    !% transect table real array indexes
    enum, bind(c)
        enumerator :: tt_depth = 1    ! depth derived from uniform distribution in EPA-SWMM
        enumerator :: tt_area         ! stores EPA-SWMM Transect.areaTbl
        enumerator :: tt_hydradius    ! stores EPA-SWMM Transect.hradTbl
        enumerator :: tt_width        ! stores EPA-SWMM Transect.widthTbl
        enumerator :: tt_lastplusone
    end enum
    
    integer, parameter :: Ncol_transectTable = tt_lastplusone-1

!%
!%==========================================================================
!% UNIFORM TABLE ARRAYS (for individual elements)
!%==========================================================================
!%  
    !% uniformTableI column indexes
    enum, bind(c)
        enumerator :: uti_idx = 1           ! counter of uniformTableI
        enumerator :: uti_BChead_idx        ! index in the BChead array where uti is associated with a BC
        enumerator :: uti_elem_idx          ! element index
        enumerator :: uti_lastplusone
    end enum

    integer, parameter :: Ncol_uniformTableI = uti_lastplusone - 1

    !% uniformTableR column indexes
    enum, bind(c)
        enumerator :: utr_SFmax =1             !% maximum section factor, i.e max (A R_h^(2/3)
        enumerator :: utr_QcritMax             !% maximum flowrate for critical depth, i.e., max(A sqrt(gD))
        enumerator :: utr_DepthMax             !% max depth at cross-section
        enumerator :: utr_AreaMax              !% max area at cross-sectoin
        enumerator :: utr_lastplusone
    end enum

    integer, parameter :: Ncol_uniformTableR = utr_lastplusone - 1

    !% uniformTablDataR integer array indexes
    !% These are tables where uniformly-distributed value is NOT depth
    !% but is instead Section Factor or Q critical
    enum, bind(c)
        enumerator :: utd_SF_uniform = 1            ! uniform distribution of section factor
        enumerator :: utd_SF_depth_nonuniform       ! depth column corresponding to section factor value
        enumerator :: utd_SF_area_nonuniform        ! area column corresponding to section factor value  
        enumerator :: utd_Qcrit_uniform             ! uniform distribution of critical flow
        enumerator :: utd_Qcrit_depth_nonuniform    ! depth column corresponding to Qcritical value
        enumerator :: utd_Qcrit_area_nonuniform     ! area column corresponding to Qcritical value
        enumerator :: utd_lastplusone
    end enum
  
    integer, parameter :: Ncol_uniformTableDataR = utd_lastplusone - 1
!%
!%==========================================================================
!% GEOMETRY TABLE ARRAYS (Additions to standard xsections)
!%==========================================================================
!%      
!% HOLD FOR LATER USE 20220930
    ! !% --- data types for output indexed by 3rd column in geometryTableR(:,:,:)
    ! enum, bind(c) 
    !     enumerator :: gtr3_Area_from_Depth = 1    !% column of output based on uniformly-distributed input
    !     enumerator :: gtr3_Perimeter_from_Depth   !%
    !     enumerator :: gtr3_Topwidth_from_Depth
    !     enumerator :: gtr3_HydDepth_from_Depth
    !     enumerator :: gtr3_HydRadius_from_Depth
    !     enumerator :: gtr3_lastplusone
    ! end enum
    ! integer, parameter :: Ncol3_GeometryTableR = gtr3_lastplusone - 1
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
!%
end module define_indexes
