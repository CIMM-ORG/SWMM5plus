!
! module initial_condition
!
! Provides setup of initial conditions for geometry and dynamics
!
!==========================================================================
!
 module initial_condition
! 
! initial conditions on the elements and faces
!
    use array_index
    use bc
    use data_keys
    use element_geometry
    use element_dynamics
    use face_values
    use globals
    use junction
    use setting_definition
    use utility
    
    implicit none
    
    private
    
    public  :: initial_condition_setup

    integer :: debuglevel = 0
    
 contains
!
!========================================================================== 
!==========================================================================
!
 subroutine initial_condition_setup &
    (elem2R, elem2I, elem2YN, elemMR, elemMI, elemMYN, faceR, faceI, faceYN, &
     linkR, linkI, nodeR, nodeI, bcdataDn, bcdataUp, thisTime)
 
 character(64) :: subroutine_name = 'initial_condition_setup'
 
 real,      intent(in out)  :: elem2R(:,:),  elemMR(:,:), faceR(:,:)
 integer,   intent(in out)  :: elem2I(:,:),  elemMI(:,:), faceI(:,:)
 logical,   intent(in out)  :: elem2YN(:,:), elemMYN(:,:), faceYN(:,:)

 real,                intent(in)      :: linkR(:,:), nodeR(:,:)
 integer,   target,   intent(in)      :: linkI(:,:), nodeI(:,:)
 real,                intent(in)      :: thisTime

 type(bcType),        intent(in out)      :: bcdataDn(:), bcdataUp(:)  
 
 integer :: idx
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 !% get data that can be extracted from links
 call initial_conditions_from_linkdata &
    (elem2R, elem2I, elemMR, elemMI, linkR, linkI)
   
 call initial_junction_conditions &
    (faceR, faceI, elem2R, elem2I, elemMR, elemMI, nodeR, nodeI)   

 !% set the bc elements (outside of face) to null values
 call bc_nullify_ghost_elem (elem2R, bcdataDn)
 call bc_nullify_ghost_elem (elem2R, bcdataUp)

 !% update the geometry
 call element_geometry_update &
    (elem2R, elem2I, elem2YN, e2r_Volume, &
     elemMR, elemMI, elemMYN, eMr_Volume, &
     faceR, faceI, bcdataDn, bcdataUp, thisTime, 0)

 call meta_element_assign &
    (elem2I, e2i_elem_type, e2i_meta_elem_type) 

 call meta_element_assign &
    (elemMI, eMi_elem_type, eMi_meta_elem_type) 
 
 call element_dynamics_update &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     bcdataDn, bcdataUp, e2r_Velocity, eMr_Velocity, &
     e2r_Volume, eMr_Volume, thisTime)   

 call face_meta_element_assign &
    (faceI, elem2I, N_face, fi_Melem_u, fi_Melem_d, fi_meta_etype_u, &
     fi_meta_etype_d, e2i_Meta_elem_type)      

 call face_update &
    (elem2R, elem2I, elemMR, faceR, faceI, faceYN, &
     bcdataDn, bcdataUp, e2r_Velocity, eMr_Velocity,  &
     e2r_Volume, eMr_Volume, thisTime, 0)

 

 !% set the element-specific smallvolume value
 !% HACK - THIS IS ONLY FOR RECTANGULAR ELEMENTS
 if (setting%SmallVolume%UseSmallVolumes) then
    elem2R(:,e2r_SmallVolume) = setting%SmallVolume%DepthCutoff * elem2R(:,e2r_BreadthScale) * elem2R(:,e2r_Length) 
    elemMR(:,eMr_SmallVolume) = setting%SmallVolume%DepthCutoff * elemMR(:,eMr_BreadthScale) * elemMR(:,eMr_Length) 
 else

    elem2R(:,e2r_SmallVolume) = zeroR
    elemMR(:,eMr_SmallVolume) = zeroR
 endif

 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine initial_condition_setup
!
!========================================================================== 
! PRIVATE BELOW HERE
!==========================================================================
!
 subroutine initial_conditions_from_linkdata &
    (elem2R, elem2I, elemMR, elemMI, linkR, linkI)
!
! The link data structure can store a variety of geometric data.
! This will be expanded in the future
!
 character(64) :: subroutine_name = 'initial_conditions_from_linkdata'
 
 real,      intent(in out)  :: elem2R(:,:),  elemMR(:,:) 
 integer,   intent(in out)  :: elem2I(:,:),  elemMI(:,:)
    
 real,      target,   intent(in)      :: linkR(:,:)
 integer,   target,   intent(in)      :: linkI(:,:)
 
 real               :: kappa
 real,      pointer :: dup, ddn
 integer,   pointer :: Lindx, LdepthType
 integer :: ii, ei_max, mm
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 do ii=1,N_link
    Lindx      => linkI(ii,li_idx)
    LdepthType => linkI(ii,li_InitialDepthType)
    
    !% up and downstream depths on this link
    dup => linkR(ii,lr_InitialUpstreamDepth)
    ddn => linkR(ii,lr_InitialDnstreamDepth)

            
    select case (LdepthType)
    
        case (1)
    
            !%  Initial depth --------------------------------------------------
            if (linkR(ii,lr_InitialDepth) /= nullvalueR) then        
                !%  if the link has a uniform depth as an initial condition
                where (elem2I(:,e2i_link_ID) == Lindx)
                    elem2R(:,e2r_Depth) = linkR(ii,lr_InitialDepth)
                endwhere
            else
                where (elem2I(:,e2i_link_ID) == Lindx)
                    elem2R(:,e2r_Depth) = 0.5*(dup + ddn)
                endwhere
            endif
        case (2)        
            !% if the link has linearly-varying depth 
            !% depth at the downstream element (link position =1)
            where ( (elem2I(:,e2i_link_Pos) == 1) .and. (elem2I(:,e2i_link_ID) == Lindx) )
                elem2R(:,e2r_Depth) = ddn   
            endwhere
            
            !%  using a linear distribution over the links 
            ei_max = maxval(elem2I(:,e2i_link_Pos),1,elem2I(:,e2i_link_ID) == Lindx)
            do mm=2,ei_max
                !% find the element that is at the mm position in the link
                where ( (elem2I(:,e2i_link_Pos) == mm) .and. (elem2I(:,e2i_link_ID) == Lindx) )
                    ! use a linear interp
                    elem2R(:,e2r_Depth) = ddn + (dup - ddn) * real(mm-1) / real(ei_max-1)
                endwhere
            end do
            
        case (3)
            ! HACK - this needs work to make the exponent driven by dup and ddn
            ! to ensure decay to dup over length
            
            !% exponential decay


            !% if the link has linearly-varying depth 
            !% depth at the downstream element (link position =1)
            where ( (elem2I(:,e2i_link_Pos) == 1) .and. (elem2I(:,e2i_link_ID) == Lindx) )
                elem2R(:,e2r_Depth) = ddn   
            endwhere

            !%  using a linear distribution over the links 
            ei_max = maxval(elem2I(:,e2i_link_Pos),1,elem2I(:,e2i_link_ID) == Lindx)            

            do mm=2,ei_max
                        kappa = real(ei_max-1)
                    !% find the element that is at the mm position in the link
                if (ddn - dup > zeroR) then
                    !%  depth decreases exponentially going upstream
                    where ( (elem2I(:,e2i_link_Pos) == mm) .and. &
                        (elem2I(:,e2i_link_ID) == Lindx)        )
                         elem2R(:,e2r_HydDepth) = (ddn-dup) * exp(-real(mm-1)) + dup
                    endwhere
                elseif (ddn - dup < zeroR) then
                    !%  depth increases exponentially going upstream
                    where ( (elem2I(:,e2i_link_Pos) == mm) .and. &
                        (elem2I(:,e2i_link_ID) == Lindx)        )
                        elem2R(:,e2r_Depth) = dup - (dup-ddn) * exp(-real(mm-1))
                    endwhere
                else
                    !%  uniform depth
                    where ( (elem2I(:,e2i_link_Pos) == mm) .and. &
                        (elem2I(:,e2i_link_ID) == Lindx)        )
                        ! use a linear interp
                        elem2R(:,e2r_Depth) = ddn
                    endwhere
                     
                endif
                            
            end do
        end select
        !%  handle all the initial conditions that don't depend on geometry type
        !%
        where (elem2I(:,e2i_link_ID) == Lindx)
            elem2I(:,e2i_roughness_type) = linkI(ii,li_roughness_type)
            elem2R(:,e2r_Roughness)      = linkR(ii,lr_Roughness)
            elem2R(:,e2r_Flowrate)       = linkR(ii,lr_InitialFlowrate)
            elem2R(:,e2r_BreadthScale)   = linkR(ii,lr_BreadthScale)               
        endwhere
        
        ! print*, elem2R(2,e2r_Eta), 'e2r_Eta', elem2R(2, e2r_HydDepth), 'e2r_HydDepth'
        if (linkI(ii,li_geometry) == lRectangularChannel ) then
            !% handle rectangular elements
            
            where (elem2I(:,e2i_link_ID) == Lindx)
                elem2I(:,e2i_geometry)  = eRectangularChannel 
                elem2R(:,e2r_HydDepth)  = elem2R(:,e2r_Depth)
                elem2R(:,e2r_BreadthScale)  = linkR(ii,lr_BreadthScale)         
                elem2R(:,e2r_Topwidth)  = linkR(ii,lr_BreadthScale)
                elem2R(:,e2r_Eta)       = elem2R(:,e2r_Zbottom) + elem2R(:,e2r_HydDepth)
                elem2R(:,e2r_Area)      = elem2R(:,e2r_HydDepth) * elem2R(:,e2r_BreadthScale)
                elem2R(:,e2r_Volume)    = elem2R(:,e2r_Area)     * elem2R(:,e2r_Length)
                elem2R(:,e2r_Perimeter) = elem2R(:,e2r_BreadthScale) + twoR * elem2R(:,e2r_HydDepth)
            endwhere
            print*, elem2R(:,e2r_Volume)
        elseif (linkI(ii,li_geometry) == lVnotchWeir ) then
            !% handle triangular elements
            !% Talk to Ehsan about this
            where (elem2I(:,e2i_link_ID) == Lindx)
                ! All the geometry calculation here are similar to the geometry calculation in weir module
                elem2I(:,e2i_geometry)      = eVnotchWeir
                ! For weir element the hydraulic depth is the depth of water in the element 
                elem2R(:,e2r_HydDepth)      = onehalfR * elem2R(:, e2r_Depth)
                elem2R(:,e2r_BreadthScale)  = zeroR 
                elem2R(:,e2r_Area)          = setting%Weir%WeirSideSlope * elem2R(:,e2r_Depth) ** twoR 
                elem2R(:,e2r_Topwidth)      = twoR * setting%Weir%WeirSideSlope * elem2R(:,e2r_Depth)
                elem2R(:,e2r_Eta)           = elem2R(:,e2r_Zbottom) + elem2R(:,e2r_HydDepth)
                elem2R(:,e2r_Volume)        = elem2R(:,e2r_Area) * elem2R(:,e2r_Length)
                elem2R(:,e2r_Perimeter)     = twoR * elem2R(:,e2r_HydDepth) * sqrt(1 + setting%Weir%WeirSideSlope ** 2)
            endwhere
        else
            !% handle elements of other geometry types
            print *, 'error: initialization for non-rectangular elements needed in ',subroutine_name

        end if

        !%  Update velocity
        where (  (elem2I(:,e2i_link_ID) == Lindx) .and. (elem2R(:,e2r_Area) > zeroR) )
            elem2R(:,e2r_Velocity)  = elem2R(:,e2r_Flowrate) / elem2R(:,e2r_Area)
        endwhere
    
 enddo

 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine initial_conditions_from_linkdata
!
!========================================================================== 
!==========================================================================
!
 subroutine initial_junction_conditions &
    (faceR, faceI, elem2R, elem2I, elemMR, elemMI, nodeR, nodeI)
 
 character(64) :: subroutine_name = 'initial_junction_conditions'
 
 real,              intent(in out)  :: elemMR(:,:)
 real,      target, intent(in)      :: elem2R(:,:), nodeR(:,:), faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:), elemMI(:,:), nodeI(:,:), faceI(:,:)
 
 integer,   pointer :: tface, telem
 
 real   :: upvalue(upstream_face_per_elemM), dnvalue(dnstream_face_per_elemM)
  
 integer :: ii, mm 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 if (N_elemM > 0) then
    !% initialize the free surface from the average of the adjacent elements
    call junction_adjacent_element_average &
        (elem2R, elemMR, elemMI, faceI, e2r_Eta, eMr_Eta)
        
    !% initialize the branch areas to the values of the adjacent elements
    call junction_adjacent_element_values_to_branches &
        (elem2R, elemMR, elemMI, faceI, e2r_Area, eMr_AreaUp, eMr_AreaDn) 
        
    !% initialize the branch flowrates to the values of the adjacent elements
    call junction_adjacent_element_values_to_branches &
        (elem2R, elemMR, elemMI, faceI, e2r_Flowrate, eMr_FlowrateUp, eMr_FlowrateDn) 
        
    !% initialize element momentum to the average of the net upstream and downstrea
    !% fluxes
    call junction_branch_average_of_inflows_and_outflows (elemMR, elemMI)  
    
    !% here we assume the branch and junction topwidths are already initialized
    !% in a prior call to junction_geometry_setup    
    !print *, elemMR(:,eMr_Topwidth)
    !print *, elemMR(:,eMr_TopwidthAll)
    
    where (elemMI(:,eMi_elem_type) == eJunctionChannel)
        elemMR(:,eMr_HydDepth) = elemMR(:,eMr_Eta) - elemMR(:,eMr_Zbottom)
    endwhere
    
    ! HACK -- need other geometry types
    
    where ((elemMI(:,eMi_geometry) == eRectangularChannel) .and. &
           (elemMI(:,eMi_elem_type) == eJunctionChannel))
        elemMR(:,eMr_Area)      = elemMR(:,eMr_HydDepth) * elemMR(:,eMr_Topwidth)
        elemMR(:,eMr_Volume)    = elemMR(:,eMr_Area)     * elemMR(:,eMr_Length)
        elemMR(:,eMr_Perimeter) = elemMR(:,eMr_Breadthscale) + twoR * elemMR(:,eMr_HydDepth)
        elemMR(:,eMr_HydRadius) = elemMR(:,eMr_Area) / elemMR(:,eMr_Perimeter)
    endwhere
    
    !% velocities
    call junction_branch_velocities (elemMR, elemMI)
    
    where ((elemMI(:,eMi_elem_type) == eJunctionChannel) .and. &
           (elemMr(:,eMr_Area) > zeroR))
        elemMR(:,eMr_Velocity) = elemMR(:,eMr_Flowrate) / elemMR(:,eMr_Area)
    endwhere    
   
 end if

 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine initial_junction_conditions
!
!========================================================================== 
!==========================================================================
!
 subroutine meta_element_assign (elemI, ei_elem_type, ei_meta_elem_type)
!
! Assign meta element type to elements
!

 character(64) :: subroutine_name = 'meta_element_assign'
 
 integer,   target,     intent(inout)    :: elemI(:,:)
 integer,               intent(in)       :: ei_elem_type         

 integer,               intent(in)       :: ei_meta_elem_type

 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
    
 where ( (elemI(:,ei_elem_type) == eChannel)             .or. &
         (elemI(:,ei_elem_type) == ePipe)                .or. &
         (elemI(:,ei_elem_type) == eJunctionChannel)     .or. &
         (elemI(:,ei_elem_type) == eJunctionPipe)  )
        elemI(:,ei_meta_elem_type) = eHQ

 elsewhere ( (elemI(:,ei_elem_type) == eWeir)            .or. &
             (elemI(:,ei_elem_type) == eorifice)         .or. &
             (elemI(:,ei_elem_type) == ePump)  )   
        elemI(:,ei_meta_elem_type) = eQonly

 elsewhere ( (elemI(:,ei_elem_type) == eStorage) )
        elemI(:,ei_meta_elem_type) = eHonly
        
 elsewhere ( (elemI(:,ei_elem_type) == eBCup)            .or. &
             (elemI(:,ei_elem_type) == eBCdn)  )
        ! Assigning nonHQ meta elem type to boundary conditions. Confirm this!
        elemI(:,ei_meta_elem_type) = eNonHQ
 end where
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine meta_element_assign 
!
!========================================================================== 
!==========================================================================
 subroutine face_meta_element_assign &
    (faceI, elemI, N_face, fi_Melem_u, fi_Melem_d, fi_meta_etype_u, &
     fi_meta_etype_d, ei_Meta_elem_type)

 character(64) :: subroutine_name = 'interp_with_junction_upstream'
 
 integer,      target,     intent(in out)  :: faceI(:,:), elemI(:,:)
 integer,                  intent(in)      :: N_face, fi_Melem_u, fi_Melem_d
 integer,                  intent(in)      :: fi_meta_etype_u, fi_meta_etype_d
 integer,                  intent(in)      :: ei_Meta_elem_type 

 integer :: ii
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 do ii=1, N_face
    faceI(ii,fi_meta_etype_u) = elemI(faceI(ii,fi_Melem_u), ei_Meta_elem_type)
    faceI(ii,fi_meta_etype_d) = elemI(faceI(ii,fi_Melem_d), ei_Meta_elem_type)
 end do

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine face_meta_element_assign
!
!========================================================================== 
!==========================================================================
!
!
! subroutine initial_condition_setupOLD &
!    (elem2R, elemMR, elem2I, elemMI, elem2YN, elemMYN, &
!     faceR, faceI, faceYN, bcdataDn, bcdataUp)
! 
! character(64) :: subroutine_name = 'initial_condition_setup'
! 
! real,      intent(in out)  :: elem2R(:,:), elemMR(:,:), faceR(:,:)
! 
! integer,   intent(in out)  :: elem2I(:,:), elemMI(:,:)
! 
! logical,   intent(in out)  :: elem2YN(:,:), elemMYN(:,:), faceYN(:,:)
! 
! integer,   intent(in out)      :: faceI(:,:)
! 
! type(bcType),  intent(in)  :: bcdataDn(:), bcdataUp(:)
! 
! real  ::  uniform_water_depth, uniform_bottom_roughness, uniform_flowrate
!  
! integer :: mm 
!!-------------------------------------------------------------------------- 
! if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
! 
! select case (casename)
!    case ('1link_network')
!        uniform_water_depth = 0.5
!        uniform_bottom_roughness = 0.03
!        uniform_flowrate = 0.822225
!        
!        print *, 'in ',subroutine_name,'---------------------------'  
!        print *, 'setting initial uniform water depth of ',uniform_water_depth
!          
!!        call initial_condition_for_uniform_rectangular_channel &
!!            (uniform_water_depth, uniform_bottom_roughness, &
!!             elem2R, elemMR, elem2I, elemMI, elem2YN, elemMYN, &
!!             faceI, bcdataDn, bcdataUp)
!!
!        call bc_nullify_ghost_elem (elem2R, bcdataDn)
!        call bc_nullify_ghost_elem (elem2R, bcdataUp)
!        
!        call element_geometry_update &
!            (elem2R, elem2I, elem2YN, e2r_Volume, &
!             elemMR, elemMI, elemMYN, eMr_Volume, &
!             faceI, bcdataDn, bcdataUp)
!            
!        ! HACK hard code setup of flowrate
!    
!        elem2R(:,e2r_Flowrate) = uniform_flowrate
!        elem2R(:,e2r_Velocity) = elem2R(:,e2r_Flowrate) / elem2R(:,e2r_Area)
!        print *, 'setting initial flowrate of ',uniform_flowrate
!    
!        call element_dynamics_update &
!            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
!             bcdataDn, bcdataUp, e2r_Velocity, eMr_Velocity)             
!        
!        call face_update &
!            (elem2R, elemMR, faceR, faceI, faceYN, &
!             bcdataDn, bcdataUp, e2r_Velocity, eMr_Velocity)
!        
!    case ('6link_1_line_network')
!        print *, 'error initial condition not defined for ',casename,' in ',subroutine_name
!        stop
!        
!    case ('3link_Y_network')
!        print *, 'error initial condition not defined for ',casename,' in ',subroutine_name
!        stop
!        
!    case ('6link_Y_network')
!        print *, 'error initial condition not defined for ',casename,' in ',subroutine_name
!        stop
!    
!    case default
!        print *
!        print *, 'casename = ',casename
!        print *, 'error: valid casename not selected in ',subroutine_name
!        stop
! end select    
! 
! ! remainder of geometry based on branches
! call junction_geometry_from_branches (elemMR, elemMI)
!
! 
! if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
! end subroutine initial_condition_setupOLD
!!
!========================================================================== 
!==========================================================================
!
! subroutine initial_condition_for_uniform_rectangular_channelOLD &
!    (uniform_water_depth, uniform_bottom_roughness, &
!     elem2R, elemMR, elem2I, elemMI, elem2YN, elemMYN, &
!     faceI, bcdataDn, bcdataUp)
!!     
!! simple rectangular channels with initial uniform depth
!!
! character(64) :: subroutine_name = 'initial_condition_for_uniform_rectangular_channel'
! 
! real,      intent(in out)  :: elem2R(:,:), elemMR(:,:)
! 
! integer,   intent(in out)  :: elem2I(:,:), elemMI(:,:)
! 
! logical,   intent(in out)  :: elem2YN(:,:), elemMYN(:,:)
! 
! integer,   intent(in)      :: faceI(:,:)
! 
! type(bcType),  target,  intent(in)  :: bcdataDn(:), bcdataUp(:)
! 
! real,  intent(in)  :: uniform_water_depth, uniform_bottom_roughness
! 
! integer,   pointer :: eID
! 
! integer :: ii
!  
!!-------------------------------------------------------------------------- 
! if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
!
! ! 2-face elements - rectangular 
! ! note that topwidth, zbottom, length were extracted from link data in network processing
! elem2R(:,e2r_Area)         = uniform_water_depth * elem2R(:,e2r_Topwidth)
! 
! elem2R(:,e2r_Volume)       = elem2R(:,e2r_Area) * elem2R(:,e2r_Length)
! 
! elem2R(:,e2r_Perimeter)    = elem2R(:,e2r_Topwidth) + twoR * uniform_water_depth
! 
! elem2R(:,e2r_Eta)          = uniform_water_depth + elem2R(:,e2r_Zbottom)
! 
! elem2R(:,e2r_Roughness)    = uniform_bottom_roughness
! 
! elem2R(:,e2r_Flowrate)     = zeroR
! elem2R(:,e2r_Velocity)     = zeroR
!    
! ! multi-face elements - rectangular
! ! note that topwidth, zbottom, length were extracted from link data in network processing
! elemMR(:,eMr_AreaAll)      = uniform_water_depth * elemMR(:,eMr_TopwidthAll) 
!    
! elemMR(:,eMr_Eta)          = uniform_water_depth + elemMR(:,eMr_Zbottom)
! 
! elemMR(:,eMr_Perimeter)    = elemMR(:,eMr_Topwidth) + twoR * uniform_water_depth
!
! elem2R(:,eMr_Roughness)    = uniform_bottom_roughness
! 
! elemMR(:,eMr_Flowrate)     = zeroR
! elemMR(:,eMr_Velocity)     = zeroR
!    
!
! if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
! end subroutine initial_condition_for_uniform_rectangular_channelOLD
!!
!========================================================================== 
! END OF MODULE initial_condition
!==========================================================================
 end module initial_condition
