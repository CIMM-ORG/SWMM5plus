! module array_index
!
! This stores the parameters used as aliases for the column numbers in
! the elem2I, elemMI, elem2R, elemMR, elem2YN,, elemMYN, faceI, and faceR arrays.
!
! for example: elem2R(:,er_Volume) is the column with volume data for an element
! that only has 2 faces
!
! Note that this fixes the number of faces per element in elemM at 6.
! If more are needed then the code will have to be rewritten for more maps.
!
! Ben R. Hodges
! 20181105
!
!==========================================================================
module array_index

    use globals

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
    integer, parameter :: links_per_node = face_per_elemM

    !%  SIZE FOR elem2 (DO NOT CHANGE) -------------------------------------------
    !CODE NOTE: require elem_per_face = 2 unless code is restructured!
    integer, parameter :: elem_per_face           = 2  !maximum number of elements connected to face

    !%  elem2I ARRAY COLUMN INDEXES FOR INTEGER DATA -----------------------------
    integer, parameter :: e2i_idx             = 1    ! unique index of each element
    integer, parameter :: e2i_meta_elem_type  = 2    ! unique indes of each meta element type
    integer, parameter :: e2i_elem_type       = 3    ! type of element
    integer, parameter :: e2i_weir_elem_type  = 4    ! type of weir element
    integer, parameter :: e2i_orif_elem_type  = 5    ! type of orifice element
    integer, parameter :: e2i_pump_elem_type  = 6    ! type of pump element
    integer, parameter :: e2i_geometry        = 7    ! type of geometry
    integer, parameter :: e2i_roughness_type  = 8    ! roughness type
    integer, parameter :: e2i_link_ID         = 9    ! ID of link in the link/node space
    integer, parameter :: e2i_link_Pos        = 10   ! position (elem from downstream = 1 to upstream = n) in link
    integer, parameter :: e2i_Mface_u         = 11   ! map to upstream face
    integer, parameter :: e2i_Mface_d         = 12   ! map to downstream face
    integer, parameter :: e2i_solver          = 13   ! solver used SVE/AC
    integer, parameter :: e2i_idx_base1       = 14   ! number of base data stored

    integer, parameter :: e2i_temp1           = e2i_idx_base1 + 1
    integer, parameter :: e2i_temp2           = e2i_idx_base1 + 2

    integer, parameter :: e2i_n_temp          = 2    ! matching the number of e2i_tempX arrays
    integer, parameter :: e2i_idx_max  = e2i_idx_base1  +   e2i_n_temp

    ! storage for temp array index positions
    integer, dimension(e2i_n_temp) :: e2i_Temp = nullvalueI

    !%  elemMI ARRAY COLUMN INDEXES FOR INTEGER DATA -----------------------------
    ! column index for integer central data on a multi-branch junction
    integer, parameter :: eMi_idx             = 1    ! unique index of each element
    integer, parameter :: eMi_meta_elem_type  = 2    ! type of meta element
    integer, parameter :: eMi_elem_type       = 3    ! type of element
    integer, parameter :: eMi_geometry        = 4    ! type of geometry
    integer, parameter :: eMi_nfaces          = 5    ! number of faces for each element
    integer, parameter :: eMi_nfaces_u        = 6    ! number of upstream faces for each element
    integer, parameter :: eMi_nfaces_d        = 7    ! number of downstream faces for each element
    integer, parameter :: eMi_roughness_type  = 8    ! roughness type
    integer, parameter :: eMi_node_ID         = 9    ! ID of node in the link/node space
    integer, parameter :: eMi_curve_type      = 10   ! ID for storage surface area curve type. 1 for functional and 2 for tabular
    integer, parameter :: eMi_solver          = 11   ! ID for storage surface area curve type. 1 for functional and 2 for tabular
    integer, parameter :: eMi_temp1           = 12
    integer, parameter :: eMi_temp2           = 13
    integer, parameter :: eMi_idx_base1       = 13    ! number of base data stored

    integer, parameter :: eMi_n_temp          = 2     ! matching the the number of eMi_tempX arrays
    ! storage for temp array index positions
    integer, dimension(eMi_n_temp) :: eMi_Temp = nullvalueI

    ! column index for indiviual branches of multi-branch junction
    integer, parameter :: eMi_Mface_u1 = eMi_idx_base1+1 ! map to upstream face
    integer, parameter :: eMi_Mface_u2 = eMi_idx_base1+2
    integer, parameter :: eMi_Mface_u3 = eMi_idx_base1+3
    integer, parameter :: eMi_idx_base2= eMi_idx_base1 + upstream_face_per_elemM

    integer, parameter :: eMi_Mface_d1 = eMi_idx_base2+1 ! map to downstream face
    integer, parameter :: eMi_Mface_d2 = eMi_idx_base2+2
    integer, parameter :: eMi_Mface_d3 = eMi_idx_base2+3
    integer, parameter :: eMi_idx_max  = eMi_idx_base2 + dnstream_face_per_elemM

    integer, dimension(upstream_face_per_elemM) :: eMi_MfaceUp = nullvalueI
    integer, dimension(dnstream_face_per_elemM) :: eMi_MfaceDn = nullvalueI
    integer, dimension(upstream_face_per_elemM + dnstream_face_per_elemM) :: eMi_MfaceAll = nullvalueI

    !%  faceI ARRAY COLUMN INDEXES FOR INTEGER DATA -----------------------------
    ! column index for integer data in the faceI array
    integer, parameter :: fi_idx                = 1    ! unique index of each face
    integer, parameter :: fi_type               = 2    ! unique type for this face
    integer, parameter :: fi_Melem_u            = 3    ! unique map index of upstream element
    integer, parameter :: fi_Melem_d            = 4    ! unique map index of downstream element
    integer, parameter :: fi_meta_etype_u       = 5    ! type of element in nominal upstream direction
    integer, parameter :: fi_meta_etype_d       = 6    ! type of element in nominal downstream direction
    integer, parameter :: fi_etype_u            = 7    ! type of element in nominal upstream direction
    integer, parameter :: fi_etype_d            = 8    ! type of element in nominal downstream direction
    integer, parameter :: fi_branch_u           = 9    ! branch # if upstream is a junction
    integer, parameter :: fi_branch_d           = 10   ! branch # if downstream is a junciton
    integer, parameter :: fi_jump_type          = 11   ! type of hydraulic jump
    integer, parameter :: fi_node_ID            = 12   ! ID of node if face corresponds to a node in link/node space
    integer, parameter :: fi_link_ID            = 13   ! ID of link if face is interior of link
    integer, parameter :: fi_link_Pos           = 14   ! position of interior face by elem count (up to down) in link
    integer, parameter :: fi_BC_ID              = 15   ! ID of the BC if a BC face
    integer, parameter :: fi_temp1              = 16
    integer, parameter :: fi_temp2              = 17
    integer, parameter :: fi_idx_max            = 17

    integer, parameter :: fi_n_temp             = 2

    ! storage for temp array index positions
    integer, dimension(fi_n_temp) :: fi_Temp = nullvalueI

    !%  elem2R ARRAY COLUMN INDEXES FOR REAL DATA ----------------------------------
    ! column index for real data in the elem2R array
    integer, parameter :: e2r_Volume             = 1
    integer, parameter :: e2r_SmallVolume        = 2
    integer, parameter :: e2r_SmallVolumeRatio   = 3
    integer, parameter :: e2r_Flowrate           = 4    ! Flowrate being solved
    integer, parameter :: e2r_Flowrate_N0        = 5    ! Flowrate at time N
    integer, parameter :: e2r_Flowrate_N1        = 6    ! Flowrate at time N-1
    integer, parameter :: e2r_Velocity           = 7
    integer, parameter :: e2r_Timescale_Q_u      = 8
    integer, parameter :: e2r_Timescale_H_u      = 9
    integer, parameter :: e2r_Timescale_G_u      = 10
    integer, parameter :: e2r_Timescale_Q_d      = 11
    integer, parameter :: e2r_Timescale_H_d      = 12
    integer, parameter :: e2r_Timescale_G_d      = 13
    integer, parameter :: e2r_Friction           = 14
    integer, parameter :: e2r_Eta                = 15
    integer, parameter :: e2r_Eta_N0             = 16   ! Eta at time N
    integer, parameter :: e2r_Head               = 17
    integer, parameter :: e2r_Area               = 18
    integer, parameter :: e2r_Area_N0            = 19   ! Area at time N
    integer, parameter :: e2r_Area_N1            = 20   ! Area at time N-1
    integer, parameter :: e2r_FullArea           = 21   ! Full area for closed conduits
    integer, parameter :: e2r_Topwidth           = 22
    integer, parameter :: e2r_Perimeter          = 23
    integer, parameter :: e2r_Depth              = 24
    integer, parameter :: e2r_FullDepth          = 25   ! vertical opening of pipe, weir, orifice
    integer, parameter :: e2r_HydDepth           = 26
    integer, parameter :: e2r_HydRadius          = 27
    integer, parameter :: e2r_Radius             = 28   ! radius for circular geometry
    integer, parameter :: e2r_X                  = 29
    integer, parameter :: e2r_Length             = 30
    integer, parameter :: e2r_Zbottom            = 31
    integer, parameter :: e2r_Zcrown             = 32   
    integer, parameter :: e2r_BreadthScale       = 33   ! bttom breadth for trapezoidal and rectangular geometry
    integer, parameter :: e2r_Roughness          = 34
    integer, parameter :: e2r_VolumeConservation = 35
    integer, parameter :: e2r_FroudeNumber       = 36
    integer, parameter :: e2r_LeftSlope          = 37   ! for specialized geometry
    integer, parameter :: e2r_RightSlope         = 38   ! for specialized geometry
    integer, parameter :: e2r_ParabolaValue      = 39   ! for specialized geometry
    integer, parameter :: e2r_SideSlope          = 40   ! for specialized geometry
    integer, parameter :: e2r_InletOffset        = 41
    integer, parameter :: e2r_OutletOffset       = 42  
    integer, parameter :: e2r_DischargeCoeff1    = 43   ! discharge coefficient for triangular weir part or orifice element
    integer, parameter :: e2r_DischargeCoeff2    = 44   ! discharge coefficient for rectangular weir part
    integer, parameter :: e2r_EndContractions    = 45   ! End contractions for rectengular and trapezoidal weir
    integer, parameter :: e2r_elN                = 46   ! the ell (lower case L)  length scale in AC
    integer, parameter :: e2r_dHdA               = 47   ! geometric change in elevation with area
    integer, parameter :: e2r_CtestH0            = 48   ! Convergence test storage for AC: d/dtau of eta * A
    integer, parameter :: e2r_CtestQ0            = 49   ! Convergence test storage for AC: d/dtau of Q
    integer, parameter :: e2r_CtestH1            = 50   ! Convergence test storage for AC: d/dtau of eta * A
    integer, parameter :: e2r_CtestQ1            = 51   ! Convergence test storage for AC: d/dtau of Q
    integer, parameter :: e2r_idx_base1          = 51
    integer, parameter :: e2r_temp1              = e2r_idx_base1 + 1
    integer, parameter :: e2r_temp2              = e2r_idx_base1 + 2
    integer, parameter :: e2r_temp3              = e2r_idx_base1 + 3
    integer, parameter :: e2r_temp4              = e2r_idx_base1 + 4
    integer, parameter :: e2r_temp5              = e2r_idx_base1 + 5
    integer, parameter :: e2r_temp6              = e2r_idx_base1 + 6
    integer, parameter :: e2r_temp7              = e2r_idx_base1 + 7
    integer, parameter :: e2r_temp8              = e2r_idx_base1 + 8
    integer, parameter :: e2r_temp9              = e2r_idx_base1 + 9
    integer, parameter :: e2r_temp10             = e2r_idx_base1 + 10
    integer, parameter :: e2r_temp11             = e2r_idx_base1 + 11
    integer, parameter :: e2r_temp12             = e2r_idx_base1 + 12
    integer, parameter :: e2r_temp13             = e2r_idx_base1 + 13
    integer, parameter :: e2r_temp14             = e2r_idx_base1 + 14
    integer, parameter :: e2r_temp15             = e2r_idx_base1 + 15
    integer, parameter :: e2r_temp16             = e2r_idx_base1 + 16
    integer, parameter :: e2r_temp17             = e2r_idx_base1 + 17
    integer, parameter :: e2r_temp18             = e2r_idx_base1 + 18
    integer, parameter :: e2r_temp19             = e2r_idx_base1 + 19
    integer, parameter :: e2r_temp20             = e2r_idx_base1 + 20
    integer, parameter :: e2r_n_temp       = 20
    integer, parameter :: e2r_idx_max = e2r_idx_base1 + e2r_n_temp

    ! storage for temp array index positions
    integer, dimension(e2r_n_temp) :: e2r_Temp = nullvalueI

    !%  elemMR ARRAY COLUMN INDEXES FOR REAL DATA --------------------------------
    ! column index for real central data on a multi-branch junction
    integer, parameter :: eMr_Volume             = 1
    integer, parameter :: eMr_SmallVolume        = 2
    integer, parameter :: eMr_SmallVolumeRatio   = 3
    integer, parameter :: eMr_Flowrate           = 4
    integer, parameter :: eMr_Flowrate_N0        = 5    ! Flowrate at time N
    integer, parameter :: eMr_Flowrate_N1        = 6    ! Flowrate at time N-1
    integer, parameter :: eMr_Velocity           = 7
    integer, parameter :: eMr_Friction           = 8
    integer, parameter :: eMr_Eta                = 9
    integer, parameter :: eMr_Eta_N0             = 10   ! Eta at time N
    integer, parameter :: eMr_Head               = 11
    integer, parameter :: eMr_Area               = 12
    integer, parameter :: eMr_Area_N0            = 13   ! Area at time N
    integer, parameter :: eMr_Area_N1            = 14   ! Area at time N-1
    integer, parameter :: eMr_FullArea           = 15   
    integer, parameter :: eMr_SurfArea           = 16
    integer, parameter :: eMr_Topwidth           = 17
    integer, parameter :: eMr_Perimeter          = 18
    integer, parameter :: eMr_Depth              = 19
    integer, parameter :: eMr_HydDepth           = 20
    integer, parameter :: eMr_FullDepth          = 21
    integer, parameter :: eMr_HydRadius          = 22
    integer, parameter :: eMr_Radius             = 23
    integer, parameter :: eMr_X                  = 24
    integer, parameter :: eMr_Length             = 25
    integer, parameter :: eMr_Zbottom            = 26
    integer, parameter :: eMr_Zcrown             = 27
    integer, parameter :: eMr_BreadthScale       = 28
    integer, parameter :: eMr_Roughness          = 29
    integer, parameter :: eMr_VolumeConservation = 30
    integer, parameter :: eMr_FroudeNumber       = 31
    integer, parameter :: eMr_LeftSlope          = 32
    integer, parameter :: eMr_RightSlope         = 33
    integer, parameter :: eMr_ParabolaValue      = 34
    integer, parameter :: eMr_FullVolume         = 35     ! Storage full volume
    integer, parameter :: eMr_StorageConstant    = 36     ! Storage Constant for surface area
    integer, parameter :: eMr_StorageCoeff       = 37     ! Storage Coefficient for surface area
    integer, parameter :: eMr_StorageExponent    = 38     ! Storafe Exponent for surface area
    integer, parameter :: eMr_SurfaceArea        = 39     ! Storage Surface Area
    integer, parameter :: eMr_SurchargeDepth     = 40     ! Surcharge Depth
    integer, parameter :: eMr_elN                = 41     ! the elN (lower case L)  length scale in AC
    integer, parameter :: eMr_dHdA               = 42     ! geometric change in elevation with area
    integer, parameter :: eMr_CtestH0            = 43     ! Convergence test storage for AC: d/dtau of eta * A
    integer, parameter :: eMr_CtestQ0            = 44     ! Convergence test storage for AC: d/dtau of Q
    integer, parameter :: eMr_CtestH1            = 45     ! Convergence test storage for AC: d/dtau of eta * A
    integer, parameter :: eMr_CtestQ1            = 46     ! Convergence test storage for AC: d/dtau of Q
    integer, parameter :: eMr_idx_base1          = 47

    ! column indexes for real branch data on a multi-branch junction
    ! note that these indexes must be consecutive by type
    ! always the upstream u1,u2,u3 then the downstream d1, d2, d3

    integer, parameter :: eMr_Eta_u1        = eMr_idx_base1 + 1    ! Q in the u1 branch
    integer, parameter :: eMr_Eta_u2        = eMr_idx_base1 + 2
    integer, parameter :: eMr_Eta_u3        = eMr_idx_base1 + 3
    integer, parameter :: eMr_idx_base2     = eMr_idx_base1 + upstream_face_per_elemM
    integer, parameter :: eMr_Eta_d1        = eMr_idx_base2 + 1
    integer, parameter :: eMr_Eta_d2        = eMr_idx_base2 + 2
    integer, parameter :: eMr_Eta_d3        = eMr_idx_base2 + 3
    integer, parameter :: eMr_idx_base3     = eMr_idx_base2 + dnstream_face_per_elemM
    integer, parameter :: eMr_Flowrate_u1   = eMr_idx_base3 + 1    ! Q in the u1 branch
    integer, parameter :: eMr_Flowrate_u2   = eMr_idx_base3 + 2
    integer, parameter :: eMr_Flowrate_u3   = eMr_idx_base3 + 3
    integer, parameter :: eMr_idx_base4     = eMr_idx_base3 + upstream_face_per_elemM
    integer, parameter :: eMr_Flowrate_d1   = eMr_idx_base4 + 1
    integer, parameter :: eMr_Flowrate_d2   = eMr_idx_base4 + 2
    integer, parameter :: eMr_Flowrate_d3   = eMr_idx_base4 + 3
    integer, parameter :: eMr_idx_base5     = eMr_idx_base4 + dnstream_face_per_elemM
    integer, parameter :: eMr_Velocity_u1   = eMr_idx_base5 + 1    ! Q in the u1 branch
    integer, parameter :: eMr_Velocity_u2   = eMr_idx_base5 + 2
    integer, parameter :: eMr_Velocity_u3   = eMr_idx_base5 + 3
    integer, parameter :: eMr_idx_base6     = eMr_idx_base5 + upstream_face_per_elemM
    integer, parameter :: eMr_Velocity_d1   = eMr_idx_base6 + 1
    integer, parameter :: eMr_Velocity_d2   = eMr_idx_base6 + 2
    integer, parameter :: eMr_Velocity_d3   = eMr_idx_base6 + 3
    integer, parameter :: eMr_idx_base7     = eMr_idx_base6 + dnstream_face_per_elemM
    integer, parameter :: eMr_Timescale_u1  = eMr_idx_base7 + 1
    integer, parameter :: eMr_Timescale_u2  = eMr_idx_base7 + 2
    integer, parameter :: eMr_Timescale_u3  = eMr_idx_base7 + 3
    integer, parameter :: eMr_idx_base8     = eMr_idx_base7  + upstream_face_per_elemM
    integer, parameter :: eMr_Timescale_d1  = eMr_idx_base8 + 1
    integer, parameter :: eMr_Timescale_d2  = eMr_idx_base8 + 2
    integer, parameter :: eMr_Timescale_d3  = eMr_idx_base8 + 3
    integer, parameter :: eMr_idx_base9     = eMr_idx_base8 + dnstream_face_per_elemM
    integer, parameter :: eMr_Area_u1       = eMr_idx_base9 + 1
    integer, parameter :: eMr_Area_u2       = eMr_idx_base9 + 2
    integer, parameter :: eMr_Area_u3       = eMr_idx_base9 + 3
    integer, parameter :: eMr_idx_base10    = eMr_idx_base9 + upstream_face_per_elemM
    integer, parameter :: eMr_Area_d1       = eMr_idx_base10 + 1
    integer, parameter :: eMr_Area_d2       = eMr_idx_base10 + 2
    integer, parameter :: eMr_Area_d3       = eMr_idx_base10 + 3
    integer, parameter :: eMr_idx_base11    = eMr_idx_base10 + dnstream_face_per_elemM
    integer, parameter :: eMr_Topwidth_u1   = eMr_idx_base11 + 1
    integer, parameter :: eMr_Topwidth_u2   = eMr_idx_base11 + 2
    integer, parameter :: eMr_Topwidth_u3   = eMr_idx_base11 + 3
    integer, parameter :: eMr_idx_base12    = eMr_idx_base11 + upstream_face_per_elemM
    integer, parameter :: eMr_Topwidth_d1   = eMr_idx_base12 + 1
    integer, parameter :: eMr_Topwidth_d2   = eMr_idx_base12 + 2
    integer, parameter :: eMr_Topwidth_d3   = eMr_idx_base12 + 3
    integer, parameter :: eMr_idx_base13    = eMr_idx_base12 + dnstream_face_per_elemM
    integer, parameter :: eMr_HydDepth_u1   = eMr_idx_base13 + 1
    integer, parameter :: eMr_HydDepth_u2   = eMr_idx_base13 + 2
    integer, parameter :: eMr_HydDepth_u3   = eMr_idx_base13 + 3
    integer, parameter :: eMr_idx_base14    = eMr_idx_base13 + upstream_face_per_elemM
    integer, parameter :: eMr_HydDepth_d1   = eMr_idx_base14 + 1
    integer, parameter :: eMr_HydDepth_d2   = eMr_idx_base14 + 2
    integer, parameter :: eMr_HydDepth_d3   = eMr_idx_base14 + 3
    integer, parameter :: eMr_idx_base15    = eMr_idx_base14 + dnstream_face_per_elemM
    integer, parameter :: eMr_Length_u1     = eMr_idx_base15 + 1
    integer, parameter :: eMr_Length_u2     = eMr_idx_base15 + 2
    integer, parameter :: eMr_Length_u3     = eMr_idx_base15 + 3
    integer, parameter :: eMr_idx_base16    = eMr_idx_base15 + upstream_face_per_elemM
    integer, parameter :: eMr_Length_d1     = eMr_idx_base16 + 1
    integer, parameter :: eMr_Length_d2     = eMr_idx_base16 + 2
    integer, parameter :: eMr_Length_d3     = eMr_idx_base16 + 3
    integer, parameter :: eMr_idx_base17    = eMr_idx_base16 + dnstream_face_per_elemM
    integer, parameter :: eMr_Zbottom_u1    = eMr_idx_base17 + 1
    integer, parameter :: eMr_Zbottom_u2    = eMr_idx_base17 + 2
    integer, parameter :: eMr_Zbottom_u3    = eMr_idx_base17 + 3
    integer, parameter :: eMr_idx_base18    = eMr_idx_base17 + upstream_face_per_elemM
    integer, parameter :: eMr_Zbottom_d1    = eMr_idx_base18 + 1
    integer, parameter :: eMr_Zbottom_d2    = eMr_idx_base18 + 2
    integer, parameter :: eMr_Zbottom_d3    = eMr_idx_base18 + 3
    integer, parameter :: eMr_idx_base19    = eMr_idx_base18 + dnstream_face_per_elemM

    integer, parameter :: eMr_BreadthScale_u1    = eMr_idx_base19 + 1
    integer, parameter :: eMr_BreadthScale_u2    = eMr_idx_base19 + 2
    integer, parameter :: eMr_BreadthScale_u3    = eMr_idx_base19 + 3
    integer, parameter :: eMr_idx_base20         = eMr_idx_base19 + upstream_face_per_elemM

    integer, parameter :: eMr_BreadthScale_d1    = eMr_idx_base20 + 1
    integer, parameter :: eMr_BreadthScale_d2    = eMr_idx_base20 + 2
    integer, parameter :: eMr_BreadthScale_d3    = eMr_idx_base20 + 3
    integer, parameter :: eMr_idx_base21         = eMr_idx_base20 + dnstream_face_per_elemM

    integer, parameter :: eMr_temp1              = eMr_idx_base21 + 1
    integer, parameter :: eMr_temp2              = eMr_idx_base21 + 2
    integer, parameter :: eMr_temp3              = eMr_idx_base21 + 3
    integer, parameter :: eMr_temp4              = eMr_idx_base21 + 4
    integer, parameter :: eMr_temp5              = eMr_idx_base21 + 5
    integer, parameter :: eMr_temp6              = eMr_idx_base21 + 6
    integer, parameter :: eMr_temp7              = eMr_idx_base21 + 7
    integer, parameter :: eMr_temp8              = eMr_idx_base21 + 8
    integer, parameter :: eMr_temp9              = eMr_idx_base21 + 9
    integer, parameter :: eMr_temp10             = eMr_idx_base21 + 10
    integer, parameter :: eMr_n_temp       = 10
    integer, parameter :: eMr_idx_base22   = eMr_idx_base21 + eMr_n_temp
    integer, parameter :: eMr_idx_max      = eMr_idx_base22

    ! storage for temp array index positions
    integer, dimension(eMr_n_temp) :: eMr_Temp = nullvalueI

    ! storage arrays for all the column indexes on upstream branches
    integer, dimension(upstream_face_per_elemM) :: eMr_EtaUp            = nullvalueI
    integer, dimension(upstream_face_per_elemM) :: eMr_FlowrateUp       = nullvalueI
    integer, dimension(upstream_face_per_elemM) :: eMr_VelocityUp       = nullvalueI
    integer, dimension(upstream_face_per_elemM) :: eMr_TimescaleUp      = nullvalueI
    integer, dimension(upstream_face_per_elemM) :: eMr_AreaUp           = nullvalueI
    integer, dimension(upstream_face_per_elemM) :: eMr_TopwidthUp       = nullvalueI
    integer, dimension(upstream_face_per_elemM) :: eMr_HydDepthUp       = nullvalueI
    integer, dimension(upstream_face_per_elemM) :: eMr_LengthUp         = nullvalueI
    integer, dimension(upstream_face_per_elemM) :: eMr_ZbottomUp        = nullvalueI
    integer, dimension(upstream_face_per_elemM) :: eMr_BreadthScaleUp   = nullvalueI
    ! storage arrays for all the column indexes on downstream branches
    integer, dimension(dnstream_face_per_elemM) :: eMr_EtaDn            = nullvalueI
    integer, dimension(dnstream_face_per_elemM) :: eMr_FlowrateDn       = nullvalueI
    integer, dimension(dnstream_face_per_elemM) :: eMr_VelocityDn       = nullvalueI
    integer, dimension(dnstream_face_per_elemM) :: eMr_TimescaleDn      = nullvalueI
    integer, dimension(dnstream_face_per_elemM) :: eMr_AreaDn           = nullvalueI
    integer, dimension(dnstream_face_per_elemM) :: eMr_TopwidthDn       = nullvalueI
    integer, dimension(dnstream_face_per_elemM) :: eMr_HydDepthDn       = nullvalueI
    integer, dimension(dnstream_face_per_elemM) :: eMr_LengthDn         = nullvalueI
    integer, dimension(dnstream_face_per_elemM) :: eMr_ZbottomDn        = nullvalueI
    integer, dimension(upstream_face_per_elemM) :: eMr_BreadthScaleDn   = nullvalueI
    ! storage arrays for all the column indexes of all branches
    integer, dimension(face_per_elemM) :: eMr_EtaAll        = nullvalueI
    integer, dimension(face_per_elemM) :: eMr_FlowrateAll   = nullvalueI
    integer, dimension(face_per_elemM) :: eMr_VelocityAll   = nullvalueI
    integer, dimension(face_per_elemM) :: eMr_TimescaleAll  = nullvalueI
    integer, dimension(face_per_elemM) :: eMr_AreaAll       = nullvalueI
    integer, dimension(face_per_elemM) :: eMr_TopwidthAll   = nullvalueI
    integer, dimension(face_per_elemM) :: eMr_HydDepthAll   = nullvalueI
    integer, dimension(face_per_elemM) :: eMr_LengthAll     = nullvalueI
    integer, dimension(face_per_elemM) :: eMr_ZbottomAll    = nullvalueI
    integer, dimension(face_per_elemM) :: eMr_BreadthScaleAll    = nullvalueI

    !%  faceR COLUMN INDEXES FOR REAL DATA ON FACE -------------------------------
    ! column index for real faceR array
    integer, parameter :: fr_Area_d       = 1 ! area on downstream side of face (AreaP in SvePy)
    integer, parameter :: fr_Area_u       = 2 ! area on upstream side of face (AreaM in SvePy)
    integer, parameter :: fr_Eta_d        = 3 ! free surface on upstream side of face
    integer, parameter :: fr_Eta_u        = 4 ! free surface on downstream side of face
    integer, parameter :: fr_Flowrate     = 5
    integer, parameter :: fr_Flowrate_N0  = 6
    integer, parameter :: fr_Flowrate_net = 7
    integer, parameter :: fr_HydDepth_d   = 8
    integer, parameter :: fr_HydDepth_u   = 9
    integer, parameter :: fr_Topwidth     = 10
    integer, parameter :: fr_Velocity_d   = 11 ! velocity on downstream side of face (velocityP in SvePy)
    integer, parameter :: fr_Velocity_u   = 12 ! velocity on upstream side of face (velocityM in SvePy)
    integer, parameter :: fr_Zbottom      = 13
    integer, parameter :: fr_X            = 14
    integer, parameter :: fr_temp1        = 15
    integer, parameter :: fr_temp2        = 16
    integer, parameter :: fr_temp3        = 17
    integer, parameter :: fr_temp4        = 18
    integer, parameter :: fr_temp5        = 19
    integer, parameter :: fr_temp6        = 20
    integer, parameter :: fr_temp7        = 21
    integer, parameter :: fr_temp8        = 22
    integer, parameter :: fr_temp9        = 23
    integer, parameter :: fr_idx_max      = 23

    integer, parameter :: fr_n_temp       = 9
    ! storage for temp array index positions
    integer, dimension(fr_n_temp) :: fr_Temp = nullvalueI

    !%  elem2YN COLUMN INDEXES FOR LOGICAL DATA ON 2 FACE ELEMENT ----------------
    ! column index for logical data in elem2YN array
    integer, parameter :: e2YN_IsSmallVolume    = 1
    integer, parameter :: e2YN_IsAdhocFlowrate  = 2
    integer, parameter :: e2YN_CanSurcharge     = 3    ! for weir surcharge
    integer, parameter :: e2YN_IsSurcharged     = 4    ! Pipe is full
    integer, parameter :: e2YN_temp1            = 5
    integer, parameter :: e2YN_temp2            = 6
    integer, parameter :: e2YN_temp3            = 7
    integer, parameter :: e2YN_temp4            = 8
    integer, parameter :: e2YN_idx_max          = 8

    integer, parameter :: e2YN_n_temp = 4
    ! storage for temp array index positions
    integer, dimension(e2YN_n_temp) :: e2YN_Temp = nullvalueI

    !%  elemMYN COLUMN INDEXES FOR LOGICAL DATA ON MULTI-BRANCH JUNCTION ---------
    ! column index for logical data in elemMYN array
    integer, parameter :: eMYN_IsSmallVolume   = 1
    integer, parameter :: eMYN_IsAdhocFlowrate = 2
    integer, parameter :: eMYN_IsSurcharged    = 3 
    integer, parameter :: eMYN_temp1           = 4
    integer, parameter :: eMYN_temp2           = 5
    integer, parameter :: eMYN_temp3           = 6
    integer, parameter :: eMYN_temp4           = 7
    integer, parameter :: eMYN_idx_max         = 7

    integer, parameter :: eMYN_n_temp = 4
    ! storage for temp array index positions
    integer, dimension(eMYN_n_temp) :: eMYN_Temp = nullvalueI

    !%  faceYN COLUMN INDEXES FOR LOGICAL DATA ON FACE
    ! column index for logical data in elemYN array
    integer, parameter :: fYN_IsHQ2up            = 1
    integer, parameter :: fYN_IsHQ2dn            = 2
    integer, parameter :: fYN_IsHup              = 3
    integer, parameter :: fYN_IsHdn              = 4
    integer, parameter :: fYN_IsQup              = 5
    integer, parameter :: fYN_IsQdn              = 6
    integer, parameter :: fYN_temp1              = 7
    integer, parameter :: fYN_temp2              = 8
    integer, parameter :: fYN_idx_max            = 8

    integer, parameter :: fYN_n_temp = 2
    ! storage for temp array index positions
    integer, dimension(fYN_n_temp) :: fYN_Temp = nullvalueI

    !%  linkI COLUMN INDEXES FOR INTEGER DATA OF LINKS IN LINK/NODE SYSTEM -------
    ! column index for integer data in linkI array
    integer, parameter :: li_idx              = 1
    integer, parameter :: li_link_type        = 2
    integer, parameter :: li_weir_type        = 3  ! type of weir link
    integer, parameter :: li_orif_type        = 4  ! type of orifice link
    integer, parameter :: li_pump_type        = 5  ! type of pump link
    integer, parameter :: li_geometry         = 6
    integer, parameter :: li_roughness_type   = 7
    integer, parameter :: li_N_element        = 8  ! Number of elements in this link
    integer, parameter :: li_Mnode_u          = 9  ! map to upstream node connecting to link
    integer, parameter :: li_Mnode_d          = 10 ! map to downstram node connecting to link
    integer, parameter :: li_Melem_u          = 11 ! element ID of upstream element of link
    integer, parameter :: li_Melem_d          = 12 ! element ID of downstream element of link
    integer, parameter :: li_Mface_u          = 13 ! face ID of upstream face of link
    integer, parameter :: li_Mface_d          = 14 ! face ID of downstream face of link
    integer, parameter :: li_assigned         = 15 ! given 1 when link is assigned
    integer, parameter :: li_InitialDepthType = 16 ! 1=uniform, 2= lineary change, 3=exponential decay
    integer, parameter :: li_temp1            = 17
    integer, parameter :: li_idx_max          = 17

    !%  nodeI COLUMN INDEXES FOR INTEGER DATA OF NODES IN LINK/NODE SYSTEM -------
    ! column index for integer data in nodeI array
    integer, parameter :: ni_idx              = 1
    integer, parameter :: ni_node_type        = 2
    integer, parameter :: ni_N_link_u         = 3 ! number of upstream links at this node
    integer, parameter :: ni_N_link_d         = 4 ! number of downstram links at this node
    integer, parameter :: ni_curve_type       = 5 ! ID for nodal storage surface area curve type. 1 for functional and 2 for tabular
    integer, parameter :: ni_assigned         = 6 ! given 1 when node has been assigned to face/elem,
    integer, parameter :: ni_temp1            = 7
    integer, parameter :: ni_idx_base1        = 7

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
    integer, parameter :: lr_Length                 = 1
    integer, parameter :: lr_InletOffset            = 2   ! Every links should have a inlet and oulet offset
    integer, parameter :: lr_OutletOffset           = 3   ! to make it consistent with SWMM.
    integer, parameter :: lr_BreadthScale           = 4
    integer, parameter :: lr_TopWidth               = 5
    integer, parameter :: lr_ElementLength          = 6
    integer, parameter :: lr_Slope                  = 7
    integer, parameter :: lr_LeftSlope              = 8
    integer, parameter :: lr_RightSlope             = 9
    integer, parameter :: lr_Roughness              = 10
    integer, parameter :: lr_InitialFlowrate        = 11
    integer, parameter :: lr_InitialDepth           = 12
    integer, parameter :: lr_InitialUpstreamDepth   = 13
    integer, parameter :: lr_InitialDnstreamDepth   = 14
    integer, parameter :: lr_ParabolaValue          = 15
    integer, parameter :: lr_SideSlope              = 16   ! for weirs only
    integer, parameter :: lr_DischargeCoeff1        = 17   ! discharge coefficient for triangular weir part or orifice element
    integer, parameter :: lr_DischargeCoeff2        = 18   ! discharge coefficient for rectangular weir part
    integer, parameter :: lr_FullDepth              = 19   ! vertical opening of pipe, weir, orifice
    integer, parameter :: lr_EndContractions        = 20
    integer, parameter :: lr_temp1                  = 21
    integer, parameter :: lr_idx_max                = 21
    
    !%  nodeR COLUMN INDEXES FOR REAL DATA OF NODES IN LINK/NODE SYSTEM ----------
    ! column index for real data in the nodeR array
    integer, parameter :: nr_Zbottom                = 1
    integer, parameter :: nr_InitialDepth           = 2
    integer, parameter :: nr_FullDepth              = 3
    integer, parameter :: nr_StorageConstant        = 4
    integer, parameter :: nr_StorageCoeff           = 5
    integer, parameter :: nr_StorageExponent        = 6
    integer, parameter :: nr_PondedArea             = 7
    integer, parameter :: nr_SurchargeDepth         = 8
    integer, parameter :: nr_temp1                  = 9
    integer, parameter :: nr_idx_base1              = 9

    ! column index for real data on multiple branches of a node
    integer, parameter :: nr_ElementLength_u1 = nr_idx_base1 + 1 ! used for subdividing junctions
    integer, parameter :: nr_ElementLength_u2 = nr_idx_base1 + 2 ! used for subdividing junctions
    integer, parameter :: nr_ElementLength_u3 = nr_idx_base1 + 3 ! used for subdividing junctions
    integer, parameter :: nr_idx_base2        = nr_idx_base1 + upstream_face_per_elemM
    integer, parameter :: nr_ElementLength_d1 = nr_idx_base2 + 1 ! used for subdividing junctions
    integer, parameter :: nr_ElementLength_d2 = nr_idx_base2 + 2 ! used for subdividing junctions
    integer, parameter :: nr_ElementLength_d3 = nr_idx_base2 + 3 ! used for subdividing junctions
    integer, parameter :: nr_idx_max = nr_idx_base2 + dnstream_face_per_elemM

    ! storage of node indexes for multi-branch data
    integer, dimension(upstream_face_per_elemM) :: nr_ElementLengthUp = nullvalueI
    integer, dimension(dnstream_face_per_elemM) :: nr_ElementLengthDn = nullvalueI

    !%  nodeYN COLUMN INDEXES FOR LOGICAL DATA ON NODES --------------------------
    ! column index for logical data in nodeYN array
    integer, parameter :: nYN_temp1     = 1
    integer, parameter :: nYN_idx_max   = 1

    !%  linkYN COLUMN INDEXES FOR LOGICAL DATA ON LINKS ---------------------------
    ! column index for logical data in linkYN array
    integer, parameter :: lYN_temp1     = 1
    integer, parameter :: lYN_idx_max   = 1

    !% widthDepth COLUMN INDEXES
    integer, parameter :: wd_widthAtLayerTop         = 1
    integer, parameter :: wd_depthAtLayerTop         = 2
    integer, parameter :: wd_areaThisLayer           = 3
    integer, parameter :: wd_areaTotalBelowThisLayer = 4
    integer, parameter :: wd_depthTotalBelowThisLayer= 5
    integer, parameter :: wd_Dwidth                  = 6
    integer, parameter :: wd_Ddepth                  = 7
    integer, parameter :: wd_angle                   = 8
    integer, parameter :: wd_perimeterBelowThisLayer = 9
    integer, parameter :: wd_rHBelowThisLayer        = 10
    integer, parameter :: wd_gammaBelowThisLayer     = 11
    integer, parameter :: wd_area_difference         = 12
    integer, parameter :: wd_local_difference        = 13
    integer, parameter :: wd_idx_max                 = 13

    !==========================================================================
    ! END OF MODULE array_index
    !==========================================================================
end module array_index
