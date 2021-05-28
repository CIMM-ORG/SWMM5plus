! module assign_index
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

module assign_index

    use globals
    use iso_c_binding

    implicit none
    !
    !==========================================================================
    !==========================================================================
    !
    !%-------------------------------------------------------------------------
    !%  FIRST INDEXES (DO NOT CHANGE) -----------------------------------------
    !% In theory, we can use different first index values for the arrays.
    !% This was initially used in debugging, but it seemed the gfortran compiler
    !% had some behaviors with non-unity starting points that I couldn't figure out.
    !%-------------------------------------------------------------------------
    integer, parameter :: first_face_index  = 1
    integer, parameter :: first_elem_index  = 1

    !%-------------------------------------------------------------------------
    !% Number of maximum branches
    integer, parameter :: max_us_branch_per_node = 3
    integer, parameter :: max_ds_branch_per_node = 3
    integer, parameter :: max_branch_per_node = max_us_branch_per_node + max_ds_branch_per_node
    !%-------------------------------------------------------------------------
    enum, bind(c)
        enumerator :: Junction_main = 1       
        enumerator :: Junction_branch_1_in
        enumerator :: Junction_branch_2_out
        enumerator :: Junction_branch_3_in
        enumerator :: Junction_branch_4_out
        enumerator :: Junction_branch_5_in
        enumerator :: Junction_branch_6_out
    end enum
    integer, target :: Nelem_in_Junction = Junction_branch_6_out
    !%-------------------------------------------------------------------------
    !% Define the column indexes for nodeI(:,:) arrays
    !% These are the for the full arrays of integer data
    !%-------------------------------------------------------------------------
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
        enumerator :: li_length_adjusted  ! 0 = length was not adjusted, 1 = one side was adjusted, 2 = both side was adjusted
        enumerator :: li_BQ_image ! image number assigned from BIPquick
    end enum
    !% note, this must be changed to whatever the last enum element is
    integer, target :: Ncol_linkI = li_BQ_image

    !%------------------------------------------------------------------------- 
    !% Define the column indexes for nodeI(:,:) arrays
    !% These are the for the full arrays of integer data
    !%-------------------------------------------------------------------------
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
    integer, parameter :: ni_idx_base1 = ni_BQ_edge

    !% column indexes for multi-branch nodes
    integer, parameter :: ni_Mlink_u1   = ni_idx_base1+1 ! map to link of upstream branch 1
    integer, parameter :: ni_Mlink_u2   = ni_idx_base1+2 ! map to link up dowstream branch 1
    integer, parameter :: ni_Mlink_u3   = ni_idx_base1+3
    integer, parameter :: ni_idx_base2  = ni_idx_base1 + max_us_branch_per_node

    integer, parameter :: ni_Mlink_d1   = ni_idx_base2+1
    integer, parameter :: ni_Mlink_d2   = ni_idx_base2+2
    integer, parameter :: ni_Mlink_d3   = ni_idx_base2+3

    !% storage for link index for upstream and downstream links
    integer, dimension(max_us_branch_per_node) :: ni_MlinkUp = nullvalueI
    integer, dimension(max_ds_branch_per_node) :: ni_MlinkDn = nullvalueI

    integer, target :: Ncol_nodeI = ni_idx_base2 + max_ds_branch_per_node

    !%------------------------------------------------------------------------- 
    !% Define the column indexes for linkR(:,:) arrays
    !% These are the for the full arrays of real data
    !%------------------------------------------------------------------------- 
    enum, bind(c)
        enumerator :: lr_Length = 1
        enumerator :: lr_AdjustedLength ! lenght adjustment if multi-link junction is present
        enumerator :: lr_InletOffset    ! Every links should have a inlet and oulet offset
        enumerator :: lr_OutletOffset   ! to make it consistent with SWMM.
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
    !% note, this must be changed to whatever the last enum element is
    integer, target :: Ncol_linkR = lr_Capacity

    !%------------------------------------------------------------------------- 
    !% Define the column indexes for nodeYN(:,:) arrays
    !% These are the for the full arrays of real data
    !%------------------------------------------------------------------------- 
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

    !% column index for real data on multiple branches of a node
    integer, parameter :: nr_ElementLength_u1 = nr_idx_base1 + 1 ! used for subdividing junctions
    integer, parameter :: nr_ElementLength_u2 = nr_idx_base1 + 2 ! used for subdividing junctions
    integer, parameter :: nr_ElementLength_u3 = nr_idx_base1 + 3 ! used for subdividing junctions

    integer, parameter :: nr_idx_base2 = nr_idx_base1 + max_us_branch_per_node
    integer, parameter :: nr_ElementLength_d1 = nr_idx_base2 + 1 ! used for subdividing junctions
    integer, parameter :: nr_ElementLength_d2 = nr_idx_base2 + 2 ! used for subdividing junctions
    integer, parameter :: nr_ElementLength_d3 = nr_idx_base2 + 3 ! used for subdividing junctions

    !% storage of node indexes for multi-branch data
    integer, dimension(max_us_branch_per_node) :: nr_ElementLengthUp = nullvalueI 
    integer, dimension(max_ds_branch_per_node) :: nr_ElementLengthDn = nullvalueI

    integer, target :: Ncol_nodeR = nr_idx_base2 + max_ds_branch_per_node
    
    !%------------------------------------------------------------------------- 
    !% Define the column indexes for nodeYN(:,:) arrays
    !% These are the for the full arrays of logical
    !%------------------------------------------------------------------------- 
    enum, bind(c)
        enumerator :: nYN_temp1 = 1
    end enum 
    !% note, this must be changed to whatever the last enum element is
    integer, target :: Ncol_nodeYN  = nYN_temp1
    
    !%------------------------------------------------------------------------- 
    !% Define the column indexes for linkYN(:,:) arrays
    !% These are the for the full arrays of logical
    !%------------------------------------------------------------------------- 
    enum, bind(c)
        enumerator :: lYN_temp1 = 1
    end enum 
    !% note, this must be changed to whatever the last enum element is
    integer, target :: Ncol_linkYN  = lYN_temp1

    !%------------------------------------------------------------------------- 
    !% COLUMN INDEXES FOR INTEGER DATA IN P_nodeI partitioning arrays
    !%------------------------------------------------------------------------- 
    enum, bind(c)
        enumerator :: P_ni_idx_Partition = 1 ! the node index number
        enumerator :: P_ni_Partition_No ! the Partition number to which that node index belongs
        enumerator :: P_ni_is_boundary ! a binary marker that is 1 when the node is shared between partitions in the link-node paradigm
    end enum

    !%------------------------------------------------------------------------- 
    !% COLUMN INDEXES FOR INTEGER DATA IN P_linkI partitioning arrays
    !%------------------------------------------------------------------------- 
    enum, bind(c)
        enumerator :: P_li_idx_Partition = 1 ! the link index number
        enumerator :: P_li_Partition_No ! the Partition number to which that link index belongs
    end enum

    !%------------------------------------------------------------------------- 
    !% Define the column indexes for elemI(:,:) array   
    !% These are for the full array of all integers
    !%-------------------------------------------------------------------------   
    enum, bind(c)
         enumerator :: ei_Lidx = 1                  !% local element index (static)
         enumerator :: ei_Gidx                      !% global element index  (static)
         enumerator :: ei_elementType               !% general element type  (static)
         enumerator :: ei_geometryType              !% cross-sectional geometry type  (static)
         enumerator :: ei_HeqType                   !% type of head equation (static)
         enumerator :: ei_link_Gidx_SWMM            !% link index from global SWMM network  (static)
         enumerator :: ei_link_Gidx_BIPquick        !% link index from global BIPquick network  (static)
         enumerator :: ei_Mface_uL                  !% map to upstream face local index  (static)
         enumerator :: ei_Mface_dL                  !% map to downstream face local index  (static)
         enumerator :: ei_node_Gidx_SWMM            !% node index from global SWMM network  (static)
         enumerator :: ei_node_Gidx_BIPquick        !% node index from global BIPquick network  (static)
         enumerator :: ei_QeqType                   !% type of flow equation (static) 
         enumerator :: ei_specificType              !% specific element type (static) 
         enumerator :: ei_tmType                    !% time march type (dynamic)
    end enum
    !% note, this must be changed to whatever the last enum element is
    integer, target :: Ncol_elemI = ei_tmType 

    !%------------------------------------------------------------------------- 
    !% Define the column indexes for elemR(:,:) array  
    !% These are for the full arrays of all reals
    !%-------------------------------------------------------------------------  
    enum, bind(c)
        enumerator :: er_Area = 1                   !% cross-sectional flow area (latest)
        enumerator :: er_Area_N0                    !% cross-sectional flow area (time N)
        enumerator :: er_Area_N1                    !% cross-sectional flow area (time N-1)
        enumerator :: er_AreaBelowBreadthMax        !% area below the max breadth in a conduit (static)
        enumerator :: er_BreadthMax                 !% maximum breadth of conduit (static)
        enumerator :: er_Depth                      !% actual maximum depth of open-channel flow
        enumerator :: er_dHdA                       !% geometric change in elevation with area
        enumerator :: er_ell                        !% the ell (lower case L)  length scale in AC solver
        enumerator :: er_Flowrate                   !% flowrate (latest)
        enumerator :: er_Flowrate_N0                !% flowrate (time N)
        enumerator :: er_Flowrate_N1                !% flowrate (time N-1)
        enumerator :: er_FlowrateLateral            !% lateral inflow BC
        enumerator :: er_FlowrateStore              !% temporary storage used for adjacent AC and ETM elements
        enumerator :: er_FroudeNumber               !% froude number of flow
        enumerator :: er_FullArea                   !% cross-sectional area of a full conduit (static)
        enumerator :: er_FullDepth                  !% maximum possible flow depth in full conduit (static)
        enumerator :: er_FullPerimeter              !% wetted perimeter of full conduit (static)
        enumerator :: er_FullVolume                 !% Volume of a full conduit (static)
        enumerator :: er_GammaC                     !% gamma continuity source term for AC solver
        enumerator :: er_GammaM                     !% gamma momentum source term for AC solver   
        enumerator :: er_Head                       !% piezometric head (latest) -- water surface elevation in open channel
        enumerator :: er_Head_N0                    !% piezometric head (time N)     
        enumerator :: er_HeadLastAC                 !% piezometric head at start of last AC step
        enumerator :: er_HeadStore                  !% temporary storage used for adjacent AC and ETM elements
        enumerator :: er_HydDepth                   !% hydraulic depth of flow    
        enumerator :: er_HydRadius                  !% hydraulic radius of flow    
        enumerator :: er_InterpWeight_uG            !% interpolation Weight, upstream, for geometry
        enumerator :: er_InterpWeight_dG            !% interpolation Weight, downstream, for geometry
        enumerator :: er_InterpWeight_uH            !% interpolation Weight, upstream for head
        enumerator :: er_InterpWeight_dH            !% interpolation Weight, downstream for head
        enumerator :: er_InterpWeight_uQ            !% interpolation Weight, upstream, for flowrate
        enumerator :: er_InterpWeight_dQ            !% interpolation Weight, downstream, for flowrate
        enumerator :: er_Ksource                    !% k source term for AC solver
        enumerator :: er_Length                     !% length of element (static)
        enumerator :: er_ones                       !% vector of ones (useful with sign function)
        enumerator :: er_Perimeter                  !% Wetted perimeter of flow
        enumerator :: er_Roughness                  !% baseline roughness value for friction model 
        enumerator :: er_SmallVolume                !% the value of a "small volume" for this element 
        enumerator :: er_SmallVolume_CMvelocity     !% velocity by Chezy-Manning for a small volume
        enumerator :: er_SmallVolume_HeadSlope      !% head slope between faces for computing Chezy-Manning on small volume
        enumerator :: er_SmallVolume_ManningsN      !% roughness used for computing Chezzy-Manning on small volume
        enumerator :: er_SmallVolumeRatio           !% blending ad hoc and solved velocity for small volume.
        enumerator :: er_SourceContinuity           !% source term for continuity equation
        enumerator :: er_SourceMomentum             !% source term for momentum equation
        enumerator :: er_Topwidth                   !% topwidth of flow at free surface 
        enumerator :: er_Velocity                   !% velocity (latest)
        enumerator :: er_Velocity_N0                !% velocity time N
        enumerator :: er_Velocity_N1                !% velocity time N-1
        enumerator :: er_VelocityLastAC             !% velocity at start of last AC step
        enumerator :: er_Volume                     !% volume (latest)
        enumerator :: er_Volume_N0                  !% volume (time N)
        enumerator :: er_Volume_N1                  !% volume (time N-1)
        enumerator :: er_VolumeLastAC               !% volume at start of last AC step
        enumerator :: er_VolumeStore                !% temporary storage used for adjacent AC and ETM elements
        enumerator :: er_WaveSpeed                  !% wave speed in element
        enumerator :: er_Zbottom                    !% bottom elevation of element (static)
        enumerator :: er_ZbreadthMax                !% elevation at maximum breadth
        enumerator :: er_Zcrown                     !% inside crown elevation of closed conduit (static)
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, target :: Ncol_elemR = er_Zcrown 

    !%-------------------------------------------------------------------------
    !% Define the column indexes for elemP(:,:) array    
    !% These are the for the packed arrays general elements
    !%-------------------------------------------------------------------------  
    enum, bind(c)
        enumerator :: ep_AC = 1                     !% all AC elements
        enumerator :: ep_ALLtm                      !% all ETM, AC elements
        enumerator :: ep_CC_AC                      !% all CC elements that are AC
        enumerator :: ep_CC_ALLtm                   !% all CC elements that are ETM or AC
        enumerator :: ep_CC_ETM                     !% all CC elements that are ETM
        enumerator :: ep_CC_H_ETM                   !% all CC elements that are ETM for H
        enumerator :: ep_CC_Q_AC                    !% all CC elements that are AC for Q
        enumerator :: ep_CC_Q_ETM                   !% all CC elements that are ETM for Q
        enumerator :: ep_CCJB_AC_surcharged         !% all CC and JB elements that are AC       
        enumerator :: ep_CCJB_ALLtm                 !% all CC and JB elements that ar any TM
        enumerator :: ep_CCJB_AC                    !% all CC and JB elements that are AC
        enumerator :: ep_CCJB_ALLtm_surcharged      !% all CC and JB elements that are AC and surcharged
        enumerator :: ep_CCJB_eETM_i_fAC            !% Any CC or JB element that is ETM and has an AC face.
        enumerator :: ep_CCJB_ETM                   !% CC and JB elements that are ETM
        enumerator :: ep_CCJB_ETM_surcharged        !% CC and JB elements that are ETM and surcharged
        enumerator :: ep_CCJM_H_AC_open             !% CC and JM elements that are AC for H and open channel
        enumerator :: ep_CCJM_H_ETM                 !% CC and JM elements that are ETM for H
        enumerator :: ep_Diag                       !% diagnostic elements (static)
        enumerator :: ep_ETM                        !% all ETM elements
        enumerator :: ep_JM_AC                      !% junction mains using AC method
        enumerator :: ep_JM_ALLtm                   !% Junction mains with any time march (static)
        enumerator :: ep_JM_ETM                     !% junction mains using ETM method    
        enumerator :: ep_NonSurcharged_AC           !% all surcharged with AC 
        enumerator :: ep_NonSurcharged_ALLtm        !% all time march nonsurcharged     
        enumerator :: ep_NonSurcharged_ETM          !% all surcharged with ETM
        enumerator :: ep_smallvolume_AC             !% small volume cells with AC
        enumerator :: ep_smallvolume_ALLtm          !% small volume with any time march
        enumerator :: ep_smallvolume_ETM            !% small volume cells with ETM    
        enumerator :: ep_Surcharged_AC              !% all surcharged with AC
        enumerator :: ep_Surcharged_ALLtm           !% all time march surcharged
        enumerator :: ep_Surcharged_ETM             !% all surcharged with ETM
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, target :: Ncol_elemP = ep_Surcharged_ETM

    !%-------------------------------------------------------------------------  
    !% Define the column indexes for elemPGalltm(:,:), elemPGetm(:,:),
    !% and elemPGac(:,:) arrays    
    !% These are the packed arrays of geometry
    !%-------------------------------------------------------------------------  

    enum, bind(c)
        enumerator :: epg_CCJM_rectangular_nonsurcharged = 1 !% CC and JM rectangular channels without surcharge
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, target :: Ncol_elemPG =  epg_CCJM_rectangular_nonsurcharged

    !%-------------------------------------------------------------------------  
    !% Define the column indexes for elemYN(:,:) arrays    
    !% These are the for the full arrays of logical
    !%-------------------------------------------------------------------------  

    enum, bind(c)
        enumerator :: eYN_canSurcharge                  !% TRUE for element that can surcharge, FALSE where it cannot (static)
        enumerator :: eYN_isAdhocFlowrate               !% TRUE is use of ad hoc flowrate algorithm    
        enumerator :: eYN_isSmallVolume                 !% TRUE is use small volume algorithm    
        enumerator :: eYN_isSurcharged                  !% TRUE is a surcharged conduit, FALSE is open channel flow
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, target :: Ncol_elemYN = eYN_isSurcharged

    !%-------------------------------------------------------------------------  
    !% Define the column indexes for elemSI(:,:) arrays    
    !% These are the full arrays if special integer data
    !%------------------------------------------------------------------------- 

    enum, bind(c)
        !% define the column indexes for elemSI(:,:) junction elements
        enumerator ::  eSI_JunctionBranch_Exists = 1    !% assigned 1 if branch exists
    end enum
    !% note, this must be changed to whatever the last enum element is
    integer, parameter :: Ncol_elemSI_junction = eSI_JunctionBranch_Exists

    enum, bind(c)
        !% define the column indexes for elemSI(:,:) weir elements
        enumerator :: eSi_Weir_EndContractions = 1      !% number of endcontractions of the weir
        enumerator :: eSi_Weir_FlowDirection            !% weir flow direction
        enumerator :: eSi_Weir_SpecificType             !% specifice weir type
    end enum
    !% note, this must be changed to whatever the last enum element is
    integer, parameter :: Ncol_elemSI_weir = eSi_Weir_SpecificType

    !% determine the largest number of columns for a special set
    integer, target :: Ncol_elemSI = max(&
                            Ncol_elemSI_junction, &
                            Ncol_elemSI_weir)

    !%-------------------------------------------------------------------------  
    !% Define the column indexes for elemSR(:,:) arrays    
    !% These are the full arrays if special real data
    !%------------------------------------------------------------------------- 

    !% define the column indexes for elemSR(:,:) for geometry that has not yet been confirmed and assigned:
    enum, bind(c)
        enumerator ::  eSr_Weir_DischargeCoeff1 = 1     !% discharge coefficient for triangular weir 
        enumerator ::  eSr_Weir_DischargeCoeff2         !% discharge coefficient for rectangular weir part
        enumerator ::  eSr_Weir_EffectiveFullDepth      !% effective full depth after control intervention
        enumerator ::  eSr_Weir_EffectiveHead           !% effective head on weir
        enumerator ::  eSr_Weir_NominalDownstreamHead   !% nominal downstream head
        enumerator ::  eSr_Weir_RectangularBreadth      !% rectangular weir breadth
        enumerator ::  eSr_Weir_TrapezoidalBreadth      !% trapezoidal weir breadth
        enumerator ::  eSr_Weir_TrapezoidalLeftSlope    !% trapezoidal weir left slope
        enumerator ::  eSr_Weir_TrapezoidalRightSlope   !% trapezoidal weir right slope
        enumerator ::  eSr_Weir_TriangularSideSlope     !% triangular weir side slope
        enumerator ::  eSr_Weir_Zcrest                  !% weir crest elevation 
    end enum
    !% note, this must be changed to whatever the last enum element is
    integer, parameter :: Ncol_elemSR_Weir = eSr_Weir_Zcrest

    enum, bind(c)
        !% Define the column indexes for elemSR(:,:) for junction branches
        enumerator :: eSr_JunctionBranch_Kfactor = 1 !% friction head loss factor in branch
    end enum    
    !% note, this must be changed to whatever the last enum element is
    integer, parameter :: Ncol_elemSR_JunctionBranch = eSr_JunctionBranch_Kfactor

    !% NEED OTHER SPECIAL ELEMENTS HERE

    !% determine the largest number of columns for a special set
    integer, target :: Ncol_elemSR = max(&
                            Ncol_elemSR_JunctionBranch, &
                            Ncol_elemSR_Weir)

    !% HACK: Ncol_elemSR must be updated when other special elements
    !% (i.e. orifice, pump, storage etc.) are added
        
    !%------------------------------------------------------------------------- 
    !% Define the column indexes for the elemSGR(:,:) arrays
    !% These are the full arrays of special, geometry, real data
    !%------------------------------------------------------------------------- 

    !% Define the column indexes for elemGSR(:,:) for rectangular pipe or channel
    enum, bind(c)
         enumerator ::  eSGR_Rectangular_Breadth = 1  !% breadth for rectangular geometry
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, parameter :: Ncol_elemSGR_Rectangular =  eSGr_Rectangular_Breadth


    !% Define the column indexes for elemSGR(:,:) for other geometry

    !% NEED OTHER GEOMETRY HERE

    !% determine the largest number of columns for a special set
    integer, target :: Ncol_elemSGR = Ncol_elemSGR_Rectangular

    !% HACK: Ncol_elemSR must be updated when other geometry types
    !% (i.e. triangular, circular etc.) are added for channel or 
    !% conduit elements

    !%------------------------------------------------------------------------- 
    !% define the column indexes for elemWDR(:,:) 
    !% for width-depth pairs
    !%-------------------------------------------------------------------------

    !% HACK We are trying to reduce the amount of data stored as width-depth pairs.  
    !% This is still experimental and under development. 

    !% The elemWDI has one row for each element that has a width-depth pair, 
    !% and we provide an index to the elemI/elemR/elemYN arrays that contain 
    !% other data about this element (e.g., Mannings n). Note that we are 
    !% planning elemWDR will have more rows than elemWDI because we 
    !% need a row for each width-depth pair. We will probably need to modify 
    !% this to create a fast approach.

    !% define the column indexes for elemWDI(:,:) for width-depth pairs
    enum, bind(c)
        enumerator ::  eWDi_Melem_Lidx = 1      !% Map to local idx of element
        enumerator ::  eWDi_elemWDRidx_F        !% Location of first row in elemWDR array for this element
        enumerator ::  eWDi_elemWDRidx_L        !% Location of last row in elemWDR array for this element
        enumerator ::  eWDi_N_pair              !% Number of width-depth pairs (rows in elemWDR) for this element
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, target :: Ncol_elemWDI =  eWDi_Melem_Lidx

    !%------------------------------------------------------------------------- 
    !% define the column indexes for elemWDR(:,:) 
    !% for width-depth pairs
    !%-------------------------------------------------------------------------

    !% HACK: This is experimental for width-depth pairs. 
    !% We expect to have a row for each pair, so parsing 
    !% the data will require use of the elemWDI array.

    enum, bind(c)
        enumerator ::  eWDr_Width = 1               !% Width at a given depth
        enumerator ::  eWDr_Depth                   !% Depth at a given width 
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, target :: Ncol_elemWDR =  eWDr_Depth

    !%------------------------------------------------------------------------- 
    !% Define the column indexes for faceI(:,:) arrays 
    !% These are the full arrays of face integer data
    !%-------------------------------------------------------------------------  

    enum, bind(c)
        enumerator ::  fi_Lidx = 1                  !% local array index (row)
        enumerator ::  fi_Gidx                      !% global (unique) index
        enumerator ::  fi_BCtype                    !% type of BC on face
        enumerator ::  fi_jump_type                 !% Type of hydraulic jump    
        enumerator ::  fi_Melem_uL                  !% map to element upstream (local index)
        enumerator ::  fi_Melem_dL                  !% map to element downstream (local index)

        !% HACK: THESE MIGHT NEED TO BE RESTORED
        ! enumerator ::  fi_Melem_uG                 !% map to element upstream (global index)
        ! enumerator ::  fi_Melem_dG                 !% map to element upstream (global index)
        ! enumerator ::  fi_eHeqType_u               !% type of H solution on element upstream
        ! enumerator ::  fi_eHeqType_d               !% type of H solution on element downstream
        ! enumerator ::  fi_eQeqType_u               !% type of Q solution on element upstream
        ! enumerator ::  fi_eQeqType_d               !% type of Q solution on element downstream

    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, target :: Ncol_faceI =  fi_Melem_dL

    !%------------------------------------------------------------------------- 
    !% Define the column indexes for faceM(:,:) arrays 
    !% These are for the full arrays of face mapping data
    !%-------------------------------------------------------------------------  
    enum, bind(c)
        enumerator :: fm_all = 1
    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, target :: Ncol_faceM =  fm_all

    !%------------------------------------------------------------------------- 
    !% Define the column indexes for faceR(:,:) arrays 
    !% These are the full arrays of face mapping data
    !%------------------------------------------------------------------------- 
    enum, bind(c)
        enumerator :: fr_Area_d = 1             !% cross-sectional area on downstream side of face
        enumerator :: fr_Area_u                 !% cross-sectional area on upstream side of face
        enumerator :: fr_Flowrate               !% flowrate through face (latest)
        enumerator :: fr_Flowrate_N0            !% flowrate through face (time N)    enumerator :: fr_Head_d  !% Piezometric head on downstream side of face
        enumerator :: fr_Head_u                 !% piezometric head on upstream side of face
        enumerator :: fr_HydDepth_d             !% hydraulic Depth on downstream side of face
        enumerator :: fr_HydDepth_u             !% hydraulic Depth on upstream side of face
        enumerator :: fr_Topwidth_d             !% topwidth on downstream side of face
        enumerator :: fr_Topwidth_u             !% topwidth on upstream side of face
        enumerator :: fr_Velocity_d             !% velocity on downstream side of face
        enumerator :: fr_Velocity_u             !% velocity on upstream side of face

        !% HACK: THE FOLLOWING MAY NEED TO BE RESTORED
        ! enumerator :: fr_Zbottom_u             !% Bottom elevation on upstream side of face
        ! enumerator :: fr_Zbottom_d             !% Bottom elevation on downstream side of face
        ! enumerator :: fr_X                     !% Linear X location in system

    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, target :: Ncol_faceR =  fr_Velocity_u

    !%-------------------------------------------------------------------------
    !% Define the column indexes for faceP(:,:) arrays 
    !% These are for the packed array of face data
    !%-------------------------------------------------------------------------  
    enum, bind(c)
        enumerator :: fp_AC                     !% face with adjacent AC element
        enumerator :: fp_Diag                   !% face with adjacent diagnostic element
        enumerator :: fp_JumpDn                 !% face with hydraulic jump from nominal downstream to upstream
        enumerator :: fp_JumpUp                 !% face with hydraulic jump from nominal upstream to downstream

    end enum
    !% note, this must be changed to whatever the last enum element is!
    integer, target :: Ncol_faceP =  fp_JumpUp

    !%-------------------------------------------------------------------------  
    !% Define the column indexes for faceYN(:,:) arrays
    !% These are the full arrays of face logical data
    !%-------------------------------------------------------------------------  
    enum, bind(c)
        enumerator :: fYN_isAC_adjacent = 1
        enumerator :: fYN_isnull    

        !% HACK: The following might not be needed
        ! enumerator :: fYN_isDiag_adjacent 
        ! enumerator :: fYN_isETM_adjacent
        ! enumerator :: fYN_isBCface 
    end enum 
    !% note, this must be changed to whatever the last enum element is!
    integer, target :: Ncol_faceYN =  fYN_isnull
    !
    !==========================================================================
    ! END OF MODULE
    !==========================================================================
    !
end module assign_index