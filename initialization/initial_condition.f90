!
!% module initial_condition
!% Provides setup of initial conditions for geometry and dynamics
!
!
!==========================================================================
!
!==========================================================================
!
module initial_condition

    use define_indexes
    use define_keys
    use define_globals
    use define_settings
    use pack_mask_arrays
    use boundary_conditions
    use update
    use face
    use diagnostic_elements
    use geometry !BRHbugfix 20210813
    use circular_conduit
    use rectangular_channel, only: rectangular_area_from_depth
    use trapezoidal_channel, only: trapezoidal_area_from_depth
    use triangular_channel, only: triangular_area_from_depth
    use storage_geometry
    use adjust
    use xsect_tables
    use interface, only: interface_get_nodef_attribute
    use utility_profiler
    use utility_allocate
    use utility_deallocate
    use utility_crash, only: util_crashpoint
    use utility, only: util_CLprint

    implicit none

    public :: init_IC_toplevel

    private

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine init_IC_toplevel ()
        !%------------------------------------------------------------------
        !% Description:
        !% set up the initial conditions for all the elements
        !%------------------------------------------------------------------
        !% Declarations:
            integer          :: ii, iblank, whichTM
            integer, pointer :: whichSolver
            integer, allocatable :: tempP(:)  !% for debugging
            integer, pointer :: thisCol, npack, thisP(:)
            character(64)    :: subroutine_name = 'init_IC_toplevel'
        !%-------------------------------------------------------------------
        !% Preliminaries:
            if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% Aliases
            whichSolver => setting%Solver%SolverSelect
            select case (whichSolver)
            case (ETM)
                whichTM = ETM
            case (AC)
                whichTM = AC
            case (ETM_AC)
                whichTM = ALLtm
            case default
                print *, 'CODE ERROR: solver type unknown for # ', whichSolver
                print *, 'which has key ',trim(reverseKey(whichSolver))
                call util_crashpoint(478383)
                !return
            end select      
        !%-------------------------------------------------------------------    
        !% --- get data that can be extracted from links
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin init_IC_from_linkdata'
        !call util_CLprint('before init_IC_from_linkdata')
        !call sleep(1)
        call init_IC_from_linkdata ()

        !% --- set up background geometry for weir, orifice, etc.
        !%     from adjacent elements
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin init_IC_diagnostic_geometry_from_adjacent'
        !call util_CLprint('before init_IC_diagnostic_geometry_from_adjacent')
        !call sleep(1)
        call init_IC_diagnostic_geometry_from_adjacent ()

        !% --- sync after all the link data has been extracted
        !%     the junction branch data is read in from the elemR arrays which
        !%     may need inter-image communication, thus the sync call is needed
        sync all
        
        !% --- get data that can be extracted from nodes
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin init_IC_from_nodedata'
        !call util_CLprint('before init_IC_from_nodedata')
        !call sleep(1)
        call init_IC_from_nodedata ()

        !% --- set up the transect arrays
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin init_IC_elem_transect...'
        !call util_CLprint('before init_IC_elem_transect')
        !call sleep(1)
        call init_IC_elem_transect_arrays ()
        call init_IC_elem_transect_geometry ()

        !% --- identify the small and zero depths (must be done before pack)
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin adjust small/zero depth'
        !call util_CLprint('before adjust_smalldepth_identify_all')
        !call sleep(1)
        call adjust_smalldepth_identify_all ()
        call adjust_zerodepth_identify_all ()

        !% ---zero out the lateral inflow column
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin init_set_zero_lateral_inflow'
        !call util_CLprint('before init_IC_set_zero_lateral_inflow')
        !call sleep(1)
        call init_IC_set_zero_lateral_inflow ()

        !% --- update time marching type
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin init_IC_solver_select '
        !call util_CLprint(' before init_IC_solver_select')
        !call sleep(1)
        call init_IC_solver_select (whichSolver)

        !% --- set up all the static packs and masks
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin pack_mask arrays_all'
        !call util_CLprint('before pack_mask_arrays_all')
        !call sleep(1)
        call pack_mask_arrays_all ()

        !% --- initialize zerovalues for other than depth (must be done after pack)
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin IC_Zerovalues'
        !call util_CLprint('before init_IC_ZeroValues_nondepth')
        !call sleep(1)
        call init_IC_ZeroValues_nondepth ()

        !% --- set all the zero and small volumes
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin adjust small/zero depth 2'
        !call util_CLprint('before adjust_zero_and_small_depth_elem')
        !call sleep(1)
        call adjust_zero_and_small_depth_elem (ETM, .true.)
        !call util_CLprint ('here again')
        call adjust_zero_and_small_depth_face (ETM, .false.)

        !% --- get the bottom slope
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin IC bottom slope'
        !call util_CLprint('before init_IC_bottom_slope')
        !call sleep(1)
        call init_IC_bottom_slope ()

        !% --- set small volume values in elements
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin init_IC_set_SmallVolumes'
        !call util_CLprint('before init_IC_set_SmallVolumes')
        !call sleep(1)
        call init_IC_set_SmallVolumes ()

        !% --- initialize slots
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin init_IC_slot'
        !call util_CLprint('before init_IC_slot')
        !call sleep(1)
        call init_IC_slot ()

        !% --- get the velocity and any other derived data
        !%     These are data needed before bc and aux variables are updated
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin init_IC_derived_data'
        !call util_CLprint('before init_IC_derived_data')
        !call sleep(1)
        call init_IC_derived_data()

        !% --- set the reference head (based on Zbottom values)
        !%     this must be called before bc_update() so that
        !%     the timeseries for head begins correctly
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin init_reference_head'
        !call util_CLprint('before init_reference_head')
        !call sleep(1)
        call init_reference_head()

        !% --- remove the reference head values from arrays
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin init_subtract_reference_head'
        !call util_CLprint('before init_subtract_reference_head')
        !call sleep(1)
        call init_subtract_reference_head()

        !% --- initialize boundary conditions
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin init_bc'
        !call util_CLprint('before init_bc')
        !call sleep(1)
        call init_bc()
        if (crashI==1) return

        !% --- update the BC so that face interpolation works in update_aux...
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin bc_update'
        !call util_CLprint('before bc_update')
        !call sleep(1)
        call bc_update()
        if (crashI==1) return

        !% --- storing dummy values for branches that are invalid
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin branch dummy values'
        !call util_CLprint('before init_IC_branch_dummy_values')
        !call sleep(1)
        call init_IC_branch_dummy_values ()
  
        !% --- set all the auxiliary (dependent) variables
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin update_aux_variables'
        !call util_CLprint('before update_auxiliary_variables')
        !call sleep(1)
        call update_auxiliary_variables (whichTM)

        !% --- initialize old head 
        !%     HACK - make into a subroutine if more variables need initializing
        !%     after update_aux_var
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin setting old head'
        !call util_CLprint('before setting old head')
        !call sleep(1)
        elemR(:,er_Head_N0) = elemR(:,er_Head)

        !% --- update diagnostic interpolation weights
        !%     (the interpolation weights of diagnostic elements
        !%     stays the same throughout the simulation. Thus, they
        !%     are only needed to be set at the top of the simulation)
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,  'begin init_IC_diagnostic_interpolation_weights'
        !call util_CLprint('before init_IC_diagnostic_interpolation_weights')
        !call sleep(1)
        call init_IC_diagnostic_interpolation_weights()

        !% --- set small values to diagnostic element interpolation sets
        !%     Needed so that junk values does not mess up the first interpolation
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin  init_IC_small_values_diagnostic_elements'
        !call util_CLprint('before init_IC_small_values_diagnostic_elements')
        !call sleep(1)
        call init_IC_small_values_diagnostic_elements

        !% --- update faces
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin face_interpolation '
        !call util_CLprint('before face_interpolation')
        !call sleep(1)
        call face_interpolation (fp_all,ALLtm)

        !% --- update the initial condition in all diagnostic elements
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin diagnostic_toplevel'
        !call util_CLprint('before diagnostic_toplevel')
        !call sleep(1)
        call diagnostic_toplevel ()

        !% --- ensure that small and zero depth faces are correct
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin adjust small/zero depth 3'
        !call util_CLprint('before adjust_zero_and_small_depth_face')
        !call sleep(1)
        call adjust_zero_and_small_depth_face (ETM, .false.)

        

        !% ---populate er_ones columns with ones
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin init_IC_oneVectors'
        !call util_CLprint('before init_IC_oneVectors')
        !call sleep(1)
        call init_IC_oneVectors ()

        ! !% TEMPORARY TEST
        ! print *, '**************************************************************************'
        ! print *, 'TEMPORARY TEST INCREASING MANNINGS N'
        ! print *, '**************************************************************************'
        ! elemR(:,er_Roughness) = tenR * elemR(:,er_Roughness) 
    

        ! !% --- initialized the max face flowrate
        ! if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin face_flowrate_max...'
        ! !call util_CLprint()
        ! call sleep(1)
        ! call face_flowrate_max_interior (fp_all)
        ! call face_flowrate_max_shared   (fp_all)


        !% Notes on initial conditions brh20211215
        !% dHdA is not initialized in channels except where timemarch is AC

        !%-------------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%initial_condition) then
            print*, '----------------------------------------------------'
            print*, 'image = ', this_image()
            print*, '.....................elements.......................'
            print*, elemI(:,ei_elementType), 'element type'
            print*, elemI(:,ei_geometryType),'element geometry'
            print*, '-------------------Geometry Data--------------------'
            print*, elemR(:,er_Depth), 'depth'
            print*, elemR(:,er_Area), 'area'
            print*, elemR(:,er_Head), 'head'
            print*, elemR(:,er_Topwidth), 'topwidth'
            print*, elemR(:,er_HydDepth), 'hydraulic depth'
            print*, elemR(:,er_HydRadius), 'hydraulic radius'
            print*, elemR(:,er_Perimeter), 'wetted perimeter'
            print*, elemR(:,er_Volume),'volume'
            print*, '-------------------Dynamics Data--------------------'
            print*, elemR(:,er_Flowrate), 'flowrate'
            print*, elemR(:,er_Velocity), 'velocity'
            print*, elemR(:,er_FroudeNumber), 'froude Number'
            print*, elemR(:,er_InterpWeight_uQ), 'timescale Q up'
            print*, elemR(:,er_InterpWeight_dQ), 'timescale Q dn'
            print*, '..................faces..........................'
            print*, faceR(:,fr_Area_u), 'face area up'
            print*, faceR(:,fr_Area_d), 'face area dn'
            print*, faceR(:,fr_Head_u), 'face head up'
            print*, faceR(:,fr_Head_d), 'face head dn'
            print*, faceR(:,fr_Flowrate), 'face flowrate'
            print*, faceR(:,fr_Topwidth_u), 'face topwidth up'
            print*, faceR(:,fr_Topwidth_d), 'face topwidth dn'
            ! call execute_command_line('')
        end if

        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_IC_toplevel
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine init_IC_from_linkdata ()
        !%------------------------------------------------------------------
        !% Description:
        !% get the initial depth, flowrate, and geometry data from links
        !%------------------------------------------------------------------
        !% Declarations:
            integer                                     :: ii, pLink, npack
            integer, pointer                            :: thisLink, eIdx(:)
            integer, dimension(:), allocatable, target  :: packed_link_idx
            integer, dimension(:), allocatable, target  :: ePack
            integer           :: allocation_status, deallocation_status
            character(len=99) :: emsg
            character(64) :: subroutine_name = 'init_IC_from_linkdata'
        !%------------------------------------------------------------------
        !% Preliminaries
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% pack all the link indexes in an image
        packed_link_idx = pack(link%I(:,li_idx), (link%I(:,li_P_image) == this_image()))

        !% find the number of links in an image
        pLink = size(packed_link_idx)

        ! !% allocate an array for packed comparison
        ! allocate(ePack(N_elem(this_image())), stat=allocation_status, errmsg=emsg)
        ! call util_allocate_check(allocation_status, emsg, 'ePack')

        !% cycle through the links in an image
        do ii = 1,pLink
            !% necessary pointers
            thisLink    => packed_link_idx(ii)
    
            ! if (thisLink == 191) then
            !     print *, ii ,' in ',trim(subroutine_name)
            !     call util_CLprint ()
            ! end if
            
            ! if (thisLink == 191) then
            !     print *, 'calling IC_get_depth'
            !     call util_CLprint ()
            ! end if
            call init_IC_get_depth (thisLink)

            ! if (thisLink == 191) then
            !     print *, 'calling IC_get_flow_roughness_From_linkdata'
            !     call util_CLprint ()
            ! end if
            call init_IC_get_flow_roughness_from_linkdata (thisLink)

            ! if (thisLink == 191) then
            !     print *, 'calling IC_get_elemtype_from_linkdata'
            !     call util_CLprint ()
            ! end if
            call init_IC_get_elemtype_from_linkdata (thisLink)

            ! if (thisLink == 191) then
            !     print *, 'calling IC_get_geoemtry_from_linkdata'
                 !call util_CLprint ('before init_IC_get_geometry_from_linkdata')
            ! end if
            call init_IC_get_geometry_from_linkdata (thisLink)
        
            ! if (thisLink == 191) then
            !     print *, 'at end of init_IC_from_linkdata'
            !     call util_CLprint ()
            ! end if

            !%brh20211215 this stuff moved to init_IC_derived_data as it
            !% does not need to be done on a link-by-link basis.
            !call init_IC_get_channel_conduit_velocity (thisLink) 

            if ((setting%Output%Verbose) .and. (this_image() == 1)) then
                if (mod(ii,1000) == 0) then
                    print *, '... handling link ',ii
                end if
            end if

        end do
        
        !%------------------------------------------------------------------
        !% Closing
            !% deallocate the temporary array
            deallocate(packed_link_idx)
            ! deallocate(ePack, stat=deallocation_status, errmsg=emsg)
            ! call util_deallocate_check(deallocation_status, emsg, 'ePack')

            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_from_linkdata
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_get_depth (thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% get the initial depth data from links and nodes
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in)  :: thisLink
            integer              :: mm, ei_max, firstidx(1)
            integer, allocatable :: pElem(:)
            integer, pointer     :: LdepthType,  nUp, nDn
            integer, pointer     :: thisP(:)
            real(8), pointer     :: DepthUp, DepthDn
            real(8), pointer     :: zLinkUp, zLinkDn
            real(8), pointer     :: eDepth(:), eLength(:), eZbottom(:)
            real(8)              :: kappa,  headUp, headDn, linkLength, length2Here
            real(8)              :: dDelta
            
            character(64) :: subroutine_name = 'init_IC_get_depth_from_linkdata'
        !%-----------------------------------------------------------------
        !% Preliminaries
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------
        !% Aliases
            !% --- type of initial depth type
            LdepthType  => link%I(thisLink,li_InitialDepthType)
            !% --- upstream and downstream nodes
            nUp      => link%I(thisLink,li_Mnode_u)
            nDn      => link%I(thisLink,li_Mnode_d)
            !% --- link upstream and downstream bottom elevation
            zLinkUp  => link%R(thisLink,lr_ZbottomUp)
            zLinkDn  => link%R(thisLink,lr_ZbottomDn)
            !% --- depths upstream and downstream on link (not yet initialized)
            DepthUp     => link%R(thisLink,lr_InitialUpstreamDepth)
            DepthDn     => link%R(thisLink,lr_InitialDnstreamDepth)
            !% 
            eLength   => elemR(:,er_Length)
            eDepth    => elemR(:,er_Depth)
            eZbottom  => elemR(:,er_Zbottom)
        !%-----------------------------------------------------------------

        !% --- up and downstream heads on connected nodes
        headUp = node%R(nUp,nr_Zbottom) + node%R(nUp,nr_InitialDepth)
        headDn = node%R(nDn,nr_Zbottom) + node%R(nDn,nr_InitialDepth)

        !print *, 'headUP, Dn',headUp, headDn

        !% --- set upstream link depths including effects of offsets
        !%     where head upstream is less than zbottom, depth is zero
        DepthUp = max(headUp - zLinkUp, zeroR)

        !% --- set downstream link depths including effects of offsets
        !%     where downstream head is less than zbottom, depth is zero
        DepthDn = max(headDn - zLinkDn, zeroR)

        ! print *, ' '
        ! print *, 'in ',trim(subroutine_name)
        ! print *, 'thislink ',thisLink
        ! print *, 'DepthUp, Dn',DepthUp,DepthDn
        ! print *, 'ZlinkUp, Dn',zLinkUp,zLinkDn

        !% --- check for a downstream gate on the node
        !%     adjust depths and head as needed
        if (node%YN(nDn,nYN_hasFlapGate)) then
            !% --- if upstreamhead is lower than downstream head
            !%     then flap gate is closed
            if (headUp < headDn) then
                !% --- closed flap gate
                if (DepthUp == zeroR) then
                    !% --- if zero depth upstream, then downstream is also zero
                    DepthDn = zeroR
                    headDn  = zLinkDn
                else
                    !% --- if some positive depth is upstream, then
                    !%     set downstream at the same head (ponding at gate)
                    headDn  = headUp
                    DepthDn = headDn - zLinkDn
                    LdepthType = FixedHead
                end if
            end if
        end if

        ! print *, 'after: '
        ! print *, 'headUP, Dn',headUp, headDn
        ! print *, 'DepthUp, Dn',DepthUp,DepthDn


        !% --- pack the elements for this link
        pElem = pack(elemI(:,ei_Lidx), (elemI(:,ei_link_Gidx_BIPquick) == thisLink))

        ! print *, 'pElem ', pElem

        !% --- error checking, the upstream should be the first element in the pack
        firstidx = findloc(elemI(pElem,ei_link_Pos),1)
        if (firstidx(1) .ne. 1) then
            print *, 'CODE ERROR'
            print *, 'Possible problem in element ordeing in a link'
            print *, 'error with link ',trim(link%Names(thisLink)%str)
            print *, elemI(pElem,ei_link_Pos)
            call util_crashpoint(55872)
        end if

        !% --- total depth delta
        dDelta = DepthUp - DepthDn

        ! print *, 'dDelta ',dDelta

        !% --- total length of all elements in link
        linkLength = sum(eLength(pElem))

        ! print *, 'linkLength ',linkLength

        !% -- initialize length measure from the upper end of the link
        !%    to an interative element center
        length2Here = zeroR

        !% --- set single element of linearly varying or exponential to uniform depth
        if (size(pElem) == 1) then
            if (     (LdepthType == LinearlyVaryingDepth) &
                .or. (LdepthType == ExponentialDepth)      ) then
                LdepthType = UniformDepth
            end if
        end if

        !% ---set the depths in link elements from links
        select case (LdepthType)

            case (UniformDepth)
                !% --- uniform depth uses the average of upstream and downstream depths
                eDepth(pElem) = onehalfR * (DepthUp + DepthDn)

                ! print *, 'uniform depth eDepth(pElem) ',eDepth(pElem)

                ! !%  if the link has a uniform depth as an initial condition
                ! if (link%R(thisLink,lr_InitialDepth) .ne. nullvalueR) then

                !     where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                !         elemR(:,er_Depth) = link%R(thisLink,lr_InitialDepth)
                !     endwhere
                ! else
                !     where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                !         elemR(:,er_Depth) = onehalfR * (DepthUp + DepthDn)
                !     endwhere
                ! end if
        

            case (LinearlyVaryingDepth)
                !% --- linearly-varying depth distribution
                do mm=1,size(pElem)
                    !% --- use the length from upstream face to center of this element
                    length2Here       = length2Here + onehalfR * eLength(pElem(mm))
                    !% --- depth by linear interpolation
                    eDepth(pElem(mm)) = DepthUp - dDelta * length2Here /linkLength
                    !% --- add the remainder of this element to the length
                    length2Here       = length2Here + onehalfR * eLength(pElem(mm))
                end do

                ! print *, 'Linear varying eDepth(pElem)',eDepth(pElem)

                !where ( (elemI(:,ei_link_Pos) == 1) .and. (elemI(:,ei_link_Gidx_BIPquick) == thisLink) )
                !    elemR(:,er_Depth) = DepthUp
                !endwhere

                !%  using a linear distribution over the links
                !ei_max = maxval(elemI(:,ei_link_Pos), 1, elemI(:,ei_link_Gidx_BIPquick) == thisLink)

                ! do mm=2,ei_max
                !     !% find the element that is at the mm position in the link
                !     where ( (elemI(:,ei_link_Pos) == mm) .and. (elemI(:,ei_link_Gidx_BIPquick) == thisLink) )
                !         !% use a linear interpolation
                !         elemR(:,er_Depth) = DepthUp - (DepthUp - DepthDn) * real(mm - oneI) / real(ei_max - oneI)
                !     endwhere
                ! end do

            case (ExponentialDepth)
                !% --- if the link has exponentially increasing or decreasing depth

                do mm=1,size(pElem)
                    !% --- use the length from upstream face to center of this element
                    length2Here       = length2Here + onehalfR * eLength(pElem(mm))
                    !% --- normalized exponential decay
                    kappa = - exp(oneR) * length2Here / linkLength
                    !% --- depth by linear interpolation
                    eDepth(pElem(mm)) = DepthDn + dDelta * exp(-kappa)
                    !% --- add the remainder of this element to the length
                    length2Here       = length2Here + onehalfR * eLength(pElem(mm))
                end do

                ! print *, 'exponential eDepth(pElem)',eDepth(pElem)
                ! where ( (elemI(:,ei_link_Pos) == 1) .and. (elemI(:,ei_link_Gidx_BIPquick) == thisLink) )
                !     elemR(:,er_Depth) = DepthUp
                ! endwhere

                ! !% find the remaining elements in the link
                ! ei_max = maxval(elemI(:,ei_link_Pos), 1, elemI(:,ei_link_Gidx_BIPquick) == thisLink)

                ! do mm=2,ei_max
                !     kappa = real(mm - oneI)

                !     !%  depth decreases exponentially going downstream
                !     if (DepthUp - DepthDn > zeroR) then
                !         where ( (elemI(:,ei_link_Pos)           == mm      ) .and. &
                !                 (elemI(:,ei_link_Gidx_BIPquick) == thisLink) )
                !             elemR(:,er_Depth) = DepthUp - (DepthUp - DepthDn) * exp(-kappa)
                !         endwhere

                !     !%  depth increases exponentially going downstream
                !     elseif (DepthUp - DepthDn < zeroR) then
                !         where ( (elemI(:,ei_link_Pos)           == mm      ) .and. &
                !                 (elemI(:,ei_link_Gidx_BIPquick) == thisLink) )
                !             elemR(:,er_Depth) = DepthUp + (DepthDn - DepthUp) * exp(-kappa)
                !         endwhere

                !     !%  uniform depth
                !     else
                !         where ( (elemI(:,ei_link_Pos)           == mm      ) .and. &
                !                 (elemI(:,ei_link_Gidx_BIPquick) == thisLink) )
                !             elemR(:,er_Depth) = DepthUp
                !         endwhere
                !     end if
                ! end do

            case (FixedHead)    
                !% --- set the downstream depth as a fixed head (ponding)
                !%     over all the elements in the link.
                eDepth(pElem) = max(headDn - eZbottom(pElem), zeroR)

                ! print *, 'fixedhead eDepth(pElem)',eDepth(pElem)
            

            case default
                print *, 'In ', subroutine_name
                print *, 'CODE ERROR: unexpected initial depth type #', LdepthType,'  in link, ', thisLink
                print *, 'which has key ',trim(reverseKey(LdepthType)) 
                !stop 
                call util_crashpoint(83753)
                !return
        end select

        ! print *, 'at end ',eDepth(pElem)

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_depth
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_get_flow_roughness_from_linkdata (thisLink)
        !--------------------------------------------------------------------------
        !
        !% get the initial flowrate and roughness data from links
        !
        !--------------------------------------------------------------------------

            integer, intent(in) :: thisLink

            character(64) :: subroutine_name = 'init_IC_get_flow_roughness_from_linkdata'
        !--------------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !%  handle all the initial conditions that don't depend on geometry type
        where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
            elemR(:,er_Flowrate)       = link%R(thisLink,lr_InitialFlowrate)
            elemR(:,er_Flowrate_N0)    = link%R(thisLink,lr_InitialFlowrate)
            elemR(:,er_Flowrate_N1)    = link%R(thisLink,lr_InitialFlowrate)
            elemR(:,er_Roughness)      = link%R(thisLink,lr_Roughness)
        endwhere

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_flow_roughness_from_linkdata
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_get_elemtype_from_linkdata (thisLink)
        !--------------------------------------------------------------------------
        !% get the geometry data from links
        !--------------------------------------------------------------------------

            integer, intent(in) :: thisLink
            integer, pointer    :: linkType

            character(64) :: subroutine_name = 'init_IC_get_elemtype_from_linkdata'
        !--------------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% necessary pointers
        linkType      => link%I(thisLink,li_link_type)

        select case (linkType)

            case (lChannel)
                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    elemI(:,ei_elementType)     = CC
                    elemI(:,ei_HeqType)         = time_march
                    elemI(:,ei_QeqType)         = time_march
                endwhere


            case (lpipe)
                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    elemI(:,ei_elementType)     = CC
                    elemI(:,ei_HeqType)         = time_march
                    elemI(:,ei_QeqType)         = time_march
                    elemYN(:,eYN_canSurcharge)  =  .true.
                endwhere

            case (lweir)
                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    elemI(:,ei_elementType)         = weir
                    elemI(:,ei_QeqType)             = diagnostic
                    elemYN(:,eYN_canSurcharge)      = link%YN(thisLink,lYN_CanSurcharge)
                endwhere

            case (lOrifice)
                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    elemI(:,ei_elementType)            = orifice
                    elemI(:,ei_QeqType)                = diagnostic
                    elemYN(:,eYN_canSurcharge)         = .true.
                endwhere

            case (lPump)
                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    elemI(:,ei_elementType)            = pump
                    elemI(:,ei_QeqType)                = diagnostic
                    elemYN(:,eYN_canSurcharge)         = .true.
                endwhere

            case (lOutlet)
                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    elemI(:,ei_elementType)            = outlet
                    elemI(:,ei_QeqType)                = diagnostic
                    elemYN(:,eYN_canSurcharge)         = .true.
                endwhere

            case default

                print *, 'In ', subroutine_name
                print *, 'CODE ERROR: unexpected link type, ', linkType,'  in the network'
                print *, 'which has key ',trim(reverseKey(linkType))
                !stop 
                call util_crashpoint(65343)
                !return
        end select


        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_elemtype_from_linkdata
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_get_geometry_from_linkdata (thisLink)
        !--------------------------------------------------------------------------
        !% get the geometry data from links
        !--------------------------------------------------------------------------

            integer, intent(in) :: thisLink
            integer, pointer    :: linkType
            character(64) :: subroutine_name = 'init_IC_get_geometry_from_linkdata'
        !--------------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% necessary pointers
        linkType      => link%I(thisLink,li_link_type)

        ! if (thisLink == 191) then
        !     print *, 'link type',linkType,trim(reverseKey(linkType))
        ! end if


        select case (linkType)

            case (lChannel)
                !% get geomety data for channels
                call init_IC_get_channel_geometry (thisLink)

            case (lpipe)
                !% get geomety data for conduits
                call init_IC_get_conduit_geometry (thisLink)

            case (lweir)
                !% get geomety data for weirs
                call init_IC_get_weir_geometry (thisLink)

            case (lOrifice)
                !% get geomety data for orifices
                call init_IC_get_orifice_geometry (thisLink)

            case (lPump)
                !% get geomety data for pump
                call init_IC_get_pump_geometry (thisLink)

            case (lOutlet)
                !% get geomety data for outlets
                call init_IC_get_outlet_geometry (thisLink)

            case default

                print *, 'In ', subroutine_name
                print *, 'CODE ERROR: unexpected link type, ', linkType,'  in the network'
                print *, 'which has key ',trim(reverseKey(linkType))
                !stop 
                call util_crashpoint(99834)
                !return

        end select


        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_geometry_from_linkdata
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_get_channel_geometry (thisLink)
        !--------------------------------------------------------------------------
        !
        !% get the geometry data for open channel links
        !% and calculate element volumes
        !%
        !% Note that the "FullDepth" must be defined for open channels.    
        !--------------------------------------------------------------------------

            integer, intent(in) :: thisLink
            integer, pointer    :: geometryType, link_tidx
            integer :: ii, kk, thisTransectIdx, startT, endT
            integer, dimension(:), allocatable :: packIdx

            real(8), pointer :: depthnorm(:)

            character(64) :: subroutine_name = 'init_IC_get_channel_geometry'
        !--------------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% pointer to geometry type
        geometryType => link%I(thisLink,li_geometry)
        depthnorm    => elemR(:,er_Temp01)

        !print *, 'geometrytype ',geometryType,trim(reverseKey(geometryType))

        select case (geometryType)

            case (lRectangular)

                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)

                    elemI(:,ei_geometryType) = rectangular

                    !% --- independent data
                    elemSGR(:,esgr_Rectangular_Breadth) = link%R(thisLink,lr_BreadthScale)
                    elemR(:,er_BreadthMax)              = elemSGR(:,esgr_Rectangular_Breadth)
                    elemR(:,er_FullDepth)               = link%R(thisLink,lr_FullDepth)
                    elemR(:,er_FullHydDepth)            = link%R(thisLink,lr_FullDepth) 

                    !% --- dependent on full depth
                    elemR(:,er_FullPerimeter)           = twoR * link%R(thisLink,lr_FullDepth) + elemSGR(:,esgr_Rectangular_Breadth) !% 20220406brh
                    elemR(:,er_ZbreadthMax)             = elemR(:,er_FullDepth) + elemR(:,er_Zbottom)
                    elemR(:,er_Zcrown)                  = elemR(:,er_Zbottom) + elemR(:,er_FullDepth)
                    elemR(:,er_FullArea)                = elemSGR(:,esgr_Rectangular_Breadth) * elemR(:,er_FullDepth)

                    !% --- dependent on full area
                    elemR(:,er_FullVolume)              = elemR(:,er_FullArea) * elemR(:,er_Length)
                    elemR(:,er_AreaBelowBreadthMax)     = elemR(:,er_FullArea)
                    elemR(:,er_ell_max)                 = (elemR(:,er_Zcrown) - elemR(:,er_ZbreadthMax)) &
                                                         * elemR(:,er_BreadthMax)                      &
                                                         + elemR(:,er_AreaBelowBreadthMax) / elemR(:,er_BreadthMax) 
               
                    !% --- store IC data
                    elemR(:,er_Area)          = elemSGR(:,esgr_Rectangular_Breadth) * elemR(:,er_Depth)
                    elemR(:,er_Area_N0)       = elemR(:,er_Area)
                    elemR(:,er_Area_N1)       = elemR(:,er_Area)
                    elemR(:,er_Volume)        = elemR(:,er_Area) * elemR(:,er_Length)
                    elemR(:,er_Volume_N0)     = elemR(:,er_Volume)
                    elemR(:,er_Volume_N1)     = elemR(:,er_Volume)
                    
                endwhere

            case (lTrapezoidal)

                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)

                    elemI(:,ei_geometryType) = trapezoidal

                    !% --- independent data
                    elemSGR(:,esgr_Trapezoidal_Breadth)    = link%R(thisLink,lr_BreadthScale)
                    elemSGR(:,esgr_Trapezoidal_LeftSlope)  = link%R(thisLink,lr_LeftSlope)
                    elemSGR(:,esgr_Trapezoidal_RightSlope) = link%R(thisLink,lr_RightSlope)
                    elemR(:,er_FullDepth)                  = link%R(thisLink,lr_FullDepth)

                    !% --- dependent on FullDepth
                    elemR(:,er_BreadthMax)   = elemSGR(:,esgr_Trapezoidal_Breadth)        &
                                             + (  elemSGR(:,esgr_Trapezoidal_LeftSlope)    &
                                                + elemSGR(:,esgr_Trapezoidal_RightSlope)) * elemR(:,er_FullDepth)
                     
                    elemR(:,er_ZbreadthMax)  = elemR(:,er_FullDepth) + elemR(:,er_Zbottom)
                    elemR(:,er_Zcrown)       = elemR(:,er_Zbottom) + elemR(:,er_FullDepth)

                    elemR(:,er_FullArea)     = ( elemSGR(:,esgr_Trapezoidal_Breadth) &
                                                 + onehalfR  &
                                                  * (  elemSGR(:,esgr_Trapezoidal_LeftSlope) &
                                                     + elemSGR(:,esgr_Trapezoidal_RightSlope) &
                                                     ) * elemR(:,er_FullDepth) &
                                                ) * elemR(:,er_FullDepth)

                    elemR(:,er_FullPerimeter) = elemSGR(:,esgr_Trapezoidal_Breadth)                          &
                                               + elemR(:,er_FullDepth)                                        &
                                                * (   sqrt(oneR + elemSGR(:,esgr_Trapezoidal_LeftSlope )**2)  &
                                                    + sqrt(oneR + elemSGR(:,esgr_Trapezoidal_RightSlope)**2) )                             

                    !% --- depends on FullArea
                    elemR(:,er_AreaBelowBreadthMax) = elemR(:,er_FullArea)
                    elemR(:,er_FullVolume)          = elemR(:,er_FullArea) * elemR(:,er_Length)

                    !% --- depends on BreadthMax et al
                    elemR(:,er_ell_max)      = (elemR(:,er_Zcrown) - elemR(:,er_ZbreadthMax))  &
                                              * elemR(:,er_BreadthMax)                         &
                                              + elemR(:,er_AreaBelowBreadthMax) / elemR(:,er_BreadthMax) 

                    elemR(:,er_FullHydDepth) = elemR(:,er_FullArea) / elemR(:,er_BreadthMax) 

                
                    !% --- store IC data
                    elemR(:,er_Area)         = ( elemSGR(:,esgr_Trapezoidal_Breadth) &
                                                 + onehalfR  &
                                                 * (  elemSGR(:,esgr_Trapezoidal_LeftSlope)  &
                                                    + elemSGR(:,esgr_Trapezoidal_RightSlope) &
                                                    ) * elemR(:,er_Depth) &
                                               ) * elemR(:,er_Depth)

                    elemR(:,er_Area_N0)      = elemR(:,er_Area)
                    elemR(:,er_Area_N1)      = elemR(:,er_Area)
                    elemR(:,er_Volume)       = elemR(:,er_Area) * elemR(:,er_Length)
                    elemR(:,er_Volume_N0)    = elemR(:,er_Volume)
                    elemR(:,er_Volume_N1)    = elemR(:,er_Volume)
                               
                endwhere

            
            case (lTriangular)

                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)

                    elemI(:,ei_geometryType) = triangular
                    !% --- independent data
                    elemSGR(:,esgr_Triangular_TopBreadth)  = link%R(thisLink,lr_BreadthScale)
                    elemR(:,er_BreadthMax)                 = link%R(thisLink,lr_BreadthScale)
                    elemR(:,er_FullDepth)                  = link%R(thisLink,lr_FullDepth)

                    !% --- dependent on Full Depth
                    elemSGR(:,esgr_Triangular_Slope)       = elemSGR(:,esgr_Triangular_TopBreadth) / (twoR * elemR(:,er_FullDepth))
                    elemR(:,er_ZbreadthMax)                = elemR(:,er_FullDepth) + elemR(:,er_Zbottom)
                    elemR(:,er_Zcrown)                     = elemR(:,er_Zbottom) + elemR(:,er_FullDepth)
                    elemR(:,er_FullArea)                   = elemR(:,er_FullDepth) * elemR(:, er_FullDepth) * elemSGR(:,esgr_Triangular_Slope) 

                    !% --- dependent on Full Area
                    elemR(:,er_FullVolume)                 = elemR(:,er_FullArea) * elemR(:,er_Length)
                    elemR(:,er_AreaBelowBreadthMax)        = elemR(:,er_FullArea)
                    elemR(:,er_FullPerimeter)              = twoR * sqrt( ((elemSGR(:,esgr_Triangular_TopBreadth) ** twoR) / fourR) &
                                                                        + (elemR(:,er_FullDepth) ** twoR) )

                    !% -- depedendent on Breadth Max et al
                    elemR(:,er_ell_max)           = (elemR(:,er_Zcrown) - elemR(:,er_ZbreadthMax)) &
                                                   * elemR(:,er_BreadthMax)                        &
                                                   + elemR(:,er_AreaBelowBreadthMax) / elemR(:,er_BreadthMax) 

                    elemR(:,er_FullHydDepth)      = elemR(:,er_FullArea) / elemR(:,er_BreadthMax) 
                    
                    !% store IC data
                    elemR(:,er_Area)         = elemR(:,er_Depth) * elemR(:, er_Depth) * elemSGR(:,esgr_Triangular_Slope) 
                    elemR(:,er_Area_N0)      = elemR(:,er_Area)
                    elemR(:,er_Area_N1)      = elemR(:,er_Area)
                    elemR(:,er_Volume)       = elemR(:,er_Area) * elemR(:,er_Length)
                    elemR(:,er_Volume_N0)    = elemR(:,er_Volume)
                    elemR(:,er_Volume_N1)    = elemR(:,er_Volume)
                    
                endwhere

            case (lIrregular)
              
                link_tidx => link%I(thisLink,li_transect_idx)

                !% TESTING STUFF
                ! print *, 'link_tidx ',link_tidx  
                ! do ii=1,size(elemI,1)
                !     if  (elemI(ii,ei_link_Gidx_BIPquick) == thisLink) then
                !         elemR(ii,er_Depth) = (real(ii,8) / 30.d0) * link%transectR(link_tidx,tr_depthFull)
                !         print *,ii, 'depth ',elemR(ii,er_Depth)
                !     end if
                ! end do

                ! !% --- initialize transect table
                ! !%     get a packed index of elements in irregular link
                ! packIdx = pack(elemI(:,ei_Lidx),(elemI(:,ei_link_Gidx_BIPquick) == thisLink) )
                ! !% --- assign the EPA SWMM transect table index to these links
                ! elemI(packIdx,ei_link_transect_idx) = li_transect_idx

                ! !% --- store separate link transect data table for each element
                ! !%     although this is substantial extra storage, it allows us to
                ! !%     write vector processing for updating element depth from volume
                ! startT = lastTransectIdx+1
                ! thisTransectIdx = lastTransectIdx
                ! do kk=1,size(packIdx)
                !     thisTransectIdx = thisTransectIdx+1
                !     !% --- assign the element transect table index (unique for each element)
                !     elemI(packIdx(kk),ei_transect_idx) = thisTransectIdx
                !     transectTableAreaR (thisTransectIdx,:,:) = link%transectTableAreaR (link_tidx,:,:)
                !     transectTableDepthR(thisTransectIdx,:,:) = link%transectTableDepthR(link_tidx,:,:)
                ! end do
                ! lastTransectIdx = thisTransectIdx
                ! endT = lastTransectIdx
                
                !% assign non-table transect data
                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    elemI(:,ei_geometryType)  = irregular
                    elemI(:,ei_link_transect_idx)  = link_tidx
                    
                    !% store fixed data
                    elemR(:,er_BreadthMax)    = link%transectR(link_tidx,tr_widthMax)
                    elemR(:,er_FullDepth)     = link%transectR(link_tidx,tr_depthFull)
                    elemR(:,er_FullArea)      = link%transectR(link_tidx,tr_areaFull)
                    elemR(:,er_FullHydDepth)  = elemR(:,er_FullArea) / link%transectR(link_tidx,tr_widthFull)
                    elemR(:,er_ZbreadthMax)   = link%transectR(link_tidx,tr_depthAtBreadthMax) + elemR(:,er_Zbottom)
                    elemR(:,er_FullPerimeter) = elemR(:,er_FullArea) / link%transectR(link_tidx,tr_hydRadiusFull)
                    elemR(:,er_Zcrown)        = elemR(:,er_Zbottom)  + elemR(:,er_FullDepth)
                    elemR(:,er_FullVolume)    = elemR(:,er_FullArea) * elemR(:,er_Length)
                    elemR(:,er_AreaBelowBreadthMax)   =  link%transectR(link_tidx,tr_areaBelowBreadthMax)
                    elemR(:,er_ell_max) = (  (elemR(:,er_Zcrown) - elemR(:,er_ZbreadthMax)) * elemR(:,er_BreadthMax) &
                                            + elemR(:,er_AreaBelowBreadthMax )                                       &
                                          ) / elemR(:,er_BreadthMax)
                endwhere

                !% --- IC data for area and volume cannot be initialized until the transect tables are setup, which is
                !%     delayed until after the JB are initialized.

                ! !% --- store IC data
                ! !%     compute the normalized depth
                ! depthnorm(packIdx) =  elemR(packIdx,er_Depth)/elemR(packIdx,er_FullDepth)
                ! !% --- lookup the associated normalized area
                ! call xsect_table_lookup_array ( &
                !     elemR(:,er_Area), depthnorm, transectTableDepthR(:,:,tt_area), packIdx)
                ! !% --- convert normalized area to physical area   
                ! elemR(packIdx,er_Area) = elemR(packIdx,er_Area) * elemR(packIdx,er_FullArea)    
                ! !% --- store other areas and volumes
                ! elemR(packIdx,er_Area_N0)       = elemR(packIdx,er_Area)
                ! elemR(packIdx,er_Area_N1)       = elemR(packIdx,er_Area)
                ! elemR(packIdx,er_Volume)        = elemR(packIdx,er_Area) * elemR(packIdx,er_Length)
                ! elemR(packIdx,er_Volume_N0)     = elemR(packIdx,er_Volume)
                ! elemR(packIdx,er_Volume_N1)     = elemR(packIdx,er_Volume)

                ! ! do ii=1,size(packIdx)
                ! !     print *, 'here  in ',trim(subroutine_name)
                ! !     print *, ii, packIdx(ii)
                ! !     print *, elemR(packIdx(ii),er_Depth), elemR(packIdx(ii),er_Area)
                ! ! end do    
                
                ! depthnorm(packIdx) = nullvalueR
                ! deallocate(packIdx)

                ! !print *, 'working on irregular cross-sections'
                ! !call util_crashpoint(448792)
                ! !return
                
            case default

                print *, 'In, ', subroutine_name
                print *, 'CODE ERROR -- geometry type unknown for # ',geometryType
                print *, 'which has key ',trim(reverseKey(geometryType))
                !stop 
                call util_crashpoint(98734)
                !return

        end select

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_channel_geometry
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_get_conduit_geometry (thisLink)
        !--------------------------------------------------------------------------
        !
        !% get the geometry data for conduit links
        !% and calculate element volumes
        !
        !--------------------------------------------------------------------------
            integer :: ii
            integer, intent(in) :: thisLink
            integer, pointer    :: geometryType

            character(64) :: subroutine_name = 'init_IC_get_conduit_geometry'
        !--------------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% pointer to geometry type
        geometryType => link%I(thisLink,li_geometry)

        ! if (thisLink == 191) then
        !     print *, 'geometry Type',geometryType,' ', trim(reverseKey(geometryType))
        !     print *, 'depth ', elemR(2227,er_Depth)
        ! end if

        select case (geometryType)

        case (lRectangular_closed)

            where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)

                elemI(:,ei_geometryType)    = rectangular_closed
                !% store geometry specific data
                elemSGR(:,esgr_Rectangular_Breadth) = link%R(thisLink,lr_BreadthScale)
                elemR(:,er_BreadthMax)            = elemSGR(:,esgr_Rectangular_Breadth)
                elemR(:,er_Area)                  = elemSGR(:,esgr_Rectangular_Breadth) * elemR(:,er_Depth)
                elemR(:,er_Area_N0)               = elemR(:,er_Area)
                elemR(:,er_Area_N1)               = elemR(:,er_Area)
                elemR(:,er_Volume)                = elemR(:,er_Area) * elemR(:,er_Length)
                elemR(:,er_Volume_N0)             = elemR(:,er_Volume)
                elemR(:,er_Volume_N1)             = elemR(:,er_Volume)
                elemR(:,er_FullDepth)             = link%R(thisLink,lr_FullDepth)
                elemR(:,er_ZbreadthMax)           = elemR(:,er_FullDepth) + elemR(:,er_Zbottom)
                elemR(:,er_Zcrown)                = elemR(:,er_Zbottom) + elemR(:,er_FullDepth)
                elemR(:,er_FullArea)              = elemSGR(:,esgr_Rectangular_Breadth) * elemR(:,er_FullDepth)
                elemR(:,er_FullVolume)            = elemR(:,er_FullArea) * elemR(:,er_Length)
                elemR(:,er_AreaBelowBreadthMax)   = elemR(:,er_FullArea)
                elemR(:,er_ell_max)               = (elemR(:,er_Zcrown) - elemR(:,er_ZbreadthMax)) * elemR(:,er_BreadthMax) + &
                                                    elemR(:,er_AreaBelowBreadthMax) / elemR(:,er_BreadthMax)  
                elemR(:,er_FullHydDepth)          = elemR(:,er_FullDepth) 
                elemR(:,er_FullPerimeter)         = twoR * elemR(:,er_FullDepth) + elemSGR(:,esgr_Rectangular_Breadth)
            endwhere
            
        case (lCircular)

            where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)

                elemI(:,ei_geometryType)    = circular
                !% store geometry specific data
                elemSGR(:,esgr_Circular_Diameter) = link%R(thisLink,lr_BreadthScale)
                elemSGR(:,esgr_Circular_Radius)   = link%R(thisLink,lr_BreadthScale) / twoR
                elemR(:,er_FullDepth)             = link%R(thisLink,lr_FullDepth)
                elemR(:,er_Zcrown)                = elemR(:,er_Zbottom) + elemR(:,er_FullDepth)
                elemR(:,er_ZbreadthMax)           = elemR(:,er_FullDepth)/twoR + elemR(:,er_Zbottom)
                elemR(:,er_FullArea)              = pi * elemSGR(:,esgr_Circular_Radius) ** twoR
                elemR(:,er_FullVolume)            = elemR(:,er_FullArea) * elemR(:,er_Length)
                elemR(:,er_FullHydDepth)          = elemR(:,er_FullDepth)
                elemR(:,er_FullPerimeter)         = elemR(:,er_FullArea) / (onefourthR * elemR(:,er_FullDepth))
                elemR(:,er_BreadthMax)            = elemSGR(:,esgr_Circular_Diameter)
                elemR(:,er_AreaBelowBreadthMax)   = elemR(:,er_FullArea)  / twoR
                elemR(:,er_Area)                  = (elemSGR(:,esgr_Circular_Radius) **2) * &
                                (acos(1.0 - (elemR(:,er_Depth)/elemSGR(:,esgr_Circular_Radius))) - &
                                sin(2.0*acos(1.0 - (elemR(:,er_Depth)/elemSGR(:,esgr_Circular_Radius))))/2.0 )
                elemR(:,er_Area_N0)               = elemR(:,er_Area)
                elemR(:,er_Area_N1)               = elemR(:,er_Area)
                elemR(:,er_Volume)                = elemR(:,er_Area) * elemR(:,er_Length)
                elemR(:,er_Volume_N0)             = elemR(:,er_Volume)
                elemR(:,er_Volume_N1)             = elemR(:,er_Volume)
                elemR(:,er_ell_max)               = (elemR(:,er_Zcrown) - elemR(:,er_ZbreadthMax)) * elemR(:,er_BreadthMax) + &
                                                    elemR(:,er_AreaBelowBreadthMax) / elemR(:,er_BreadthMax) 
            end where

        case (lIrregular)
            print *, 'In ', trim(subroutine_name)
            print *, 'USER ERROR: Irregular cross-section geometry not allowed for closed conduits (open-channel only) in SWMM5+'
            call util_crashpoint(4409874)
            !return
        case default

            print *, 'In, ', trim(subroutine_name)
            print *, 'CODE ERROR geometry type unknown for # ', geometryType
            print *, 'which has key ',trim(reverseKey(geometryType))
            call util_crashpoint(887344)
            !return
        end select

        ! if (thisLink == 191) then
        !     print *, 'at end of init_IC_get_conduit_geometry'
        !     print *, 'geometry Type',geometryType,' ', trim(reverseKey(geometryType))
        !     print *, 'depth ', elemR(2227,er_Depth)
        ! end if

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_conduit_geometry
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_get_weir_geometry (thisLink)
        !--------------------------------------------------------------------------
        !
        !% get the geometry and other data data for weir links
        !
        !--------------------------------------------------------------------------

            integer, intent(in) :: thisLink
            integer, pointer    :: specificWeirType
            integer, allocatable :: thisPack(:)
            integer :: ii

            character(64) :: subroutine_name = 'init_IC_get_weir_geometry'
        !--------------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% pointer to specific weir type
        specificWeirType => link%I(thisLink,li_link_sub_type)

        select case (specificWeirType)
            !% --- set up weir specific data
            case (lTrapezoidalWeir)

                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    !% integer data
                    elemSI(:,esi_Weir_SpecificType)          = trapezoidal_weir
                    !% real data
                    elemSR(:,esr_Weir_EffectiveFullDepth)    = link%R(thisLink,lr_FullDepth)
                    elemSR(:,esr_Weir_Rectangular)           = link%R(thisLink,lr_DischargeCoeff1)
                    elemSR(:,esr_Weir_Triangular)            = link%R(thisLink,lr_DischargeCoeff2)
                    elemSR(:,esr_Weir_TrapezoidalBreadth)    = link%R(thisLink,lr_BreadthScale)
                    elemSR(:,esr_Weir_TrapezoidalLeftSlope)  = link%R(thisLink,lr_SideSlope)
                    elemSR(:,esr_Weir_TrapezoidalRightSlope) = link%R(thisLink,lr_SideSlope)
                    elemSR(:,esr_Weir_Zcrest)                = elemR(:,er_Zbottom) + link%R(thisLink,lr_InletOffset)
                    elemSR(:,esr_Weir_Zcrown)                = elemSR(:,esr_Weir_Zcrest) + link%R(thisLink,lr_FullDepth)

                    !% --- default channel geometry (overwritten later by adjacent CC shape)
                    !%     assumes channel is rectangular and  twice the breadth of weir and
                    !%     used weir crown and add the weir effective depth as the maximum overflow
                    elemI(:,ei_geometryType)            = rectangular
                    elemSGR(:,esgr_Rectangular_Breadth) = twoR * (  elemSR(:,esr_Weir_TrapezoidalBreadth) &
                                                               +    elemSR(:,esr_Weir_EffectiveFullDepth) &
                                                                *(  elemSR(:,esr_Weir_TrapezoidalLeftSlope) &
                                                                  + elemSR(:,esr_Weir_TrapezoidalRightSlope)) )
                    elemR(:,er_BreadthMax)              = elemSGR(:,esgr_Rectangular_Breadth)                                       
                    elemR(:,er_FullDepth)               = elemSR(:,esr_Weir_Zcrown) - elemR(:,er_Zbottom) &
                                                            + elemSR(:,esr_Weir_EffectiveFullDepth)  
                endwhere

            case (lSideFlowWeir)

                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    !% integer data
                    elemSI(:,esi_Weir_SpecificType)          = side_flow
                    elemSI(:,esi_Weir_EndContractions)       = link%I(thisLink,li_weir_EndContrations)
                    !% real data
                    elemSR(:,esr_Weir_EffectiveFullDepth)    = link%R(thisLink,lr_FullDepth)
                    elemSR(:,esr_Weir_Rectangular)           = link%R(thisLink,lr_DischargeCoeff1)
                    elemSR(:,esr_Weir_RectangularBreadth)    = link%R(thisLink,lr_BreadthScale)
                    elemSR(:,esr_Weir_Zcrest)                = elemR(:,er_Zbottom) + link%R(thisLink,lr_InletOffset)
                    elemSR(:,esr_Weir_Zcrown)                = elemSR(:,esr_Weir_Zcrest) + link%R(thisLink,lr_FullDepth)

                    !% --- default channel geometry (overwritten later by adjacent CC shape)
                    !%     assumes channel is rectangular and twice the breadth of weir and
                    !%     used weir crown and add the weir effective depth as the maximum overflow
                    elemI(:,ei_geometryType)            = rectangular
                    elemSGR(:,esgr_Rectangular_Breadth) = twoR * elemSR(:,esr_Weir_RectangularBreadth) 
                    elemR(:,er_BreadthMax)              = elemSGR(:,esgr_Rectangular_Breadth)                                       
                    elemR(:,er_FullDepth)               = elemSR(:,esr_Weir_Zcrown) - elemR(:,er_Zbottom) &
                                                            + elemSR(:,esr_Weir_EffectiveFullDepth)   
                endwhere

            case (lRoadWayWeir)

                print *, 'In ', subroutine_name
                print *, 'roadway weir is not handeled yet for #',specificWeirType
                print *, 'which has key ',trim(reverseKey(specificWeirType))
                !stop 
                call util_crashpoint(557834)
                !return

            case (lVnotchWeir)

                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    !% integer data
                    elemSI(:,esi_Weir_SpecificType)          = vnotch_weir
                    !% real data
                    elemSR(:,esr_Weir_EffectiveFullDepth)    = link%R(thisLink,lr_FullDepth)
                    elemSR(:,esr_Weir_Triangular)            = link%R(thisLink,lr_DischargeCoeff1)
                    elemSR(:,esr_Weir_TriangularSideSlope)   = link%R(thisLink,lr_SideSlope)
                    elemSR(:,esr_Weir_Zcrest)                = elemR(:,er_Zbottom) + link%R(thisLink,lr_InletOffset)
                    elemSR(:,esr_Weir_Zcrown)                = elemSR(:,esr_Weir_Zcrest) + link%R(thisLink,lr_FullDepth)

                    !% --- default channel geometry (overwritten later by adjacent CC shape)
                    !%     assumes channel is rectangular and twice the breadth of weir and
                    !%     used weir crown and add the weir effective depth as the maximum overflow
                    elemI(:,ei_geometryType)            = rectangular
                    elemSGR(:,esgr_Rectangular_Breadth) = fourR * elemSR(:,esr_Weir_EffectiveFullDepth) &
                                                                * elemSR(:,esr_Weir_TriangularSideSlope)
                    elemR(:,er_BreadthMax)              = elemSGR(:,esgr_Rectangular_Breadth)                                       
                    elemR(:,er_FullDepth)               = elemSR(:,esr_Weir_Zcrown) - elemR(:,er_Zbottom) &
                                                            + elemSR(:,esr_Weir_EffectiveFullDepth)   
                endwhere

            case (lTransverseWeir)

                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    !% integer data
                    elemSI(:,esi_Weir_SpecificType)          = transverse_weir
                    elemSI(:,esi_Weir_EndContractions)       = link%I(thisLink,li_weir_EndContrations)
                    !% real data
                    elemSR(:,esr_Weir_EffectiveFullDepth)    = link%R(thisLink,lr_FullDepth)
                    elemSR(:,esr_Weir_Rectangular)           = link%R(thisLink,lr_DischargeCoeff1)
                    elemSR(:,esr_Weir_RectangularBreadth)    = link%R(thisLink,lr_BreadthScale)
                    elemSR(:,esr_Weir_Zcrest)                = elemR(:,er_Zbottom)  + link%R(thisLink,lr_InletOffset)
                    elemSR(:,esr_Weir_Zcrown)                = elemSR(:,esr_Weir_Zcrest) + link%R(thisLink,lr_FullDepth)

                    !% --- default channel geometry (overwritten later by adjacent CC shape)
                    !%     assumes channel is rectangular and twice the breadth of weir and
                    !%     used weir crown and add the weir effective depth as the maximum overflow
                    elemI(:,ei_geometryType)            = rectangular
                    elemSGR(:,esgr_Rectangular_Breadth) = twoR * elemSR(:,esr_Weir_RectangularBreadth) 
                    elemR(:,er_BreadthMax)              = elemSGR(:,esgr_Rectangular_Breadth)                                       
                    elemR(:,er_FullDepth)               = elemSR(:,esr_Weir_Zcrown) - elemR(:,er_Zbottom) &
                                                            + elemSR(:,esr_Weir_EffectiveFullDepth)   
                endwhere

            case default

                print *, 'In ', trim(subroutine_name)
                print *, 'CODE ERROR: unknown weir type, ', specificWeirType,'  in network'
                print *, 'which has key ',trim(reverseKey(specificWeirType))
                !stop 
                call util_crashpoint(99834)
                !return

        end select

        !% --- set minimum crest height as 101% of the zero depth value for all weirs
        !%     this ensures that zero-height weir elements cannot cause flow for zerovalue depths
        thisPack = pack(elemI(:,ei_Lidx),(elemI(:,ei_link_Gidx_BIPquick) == thisLink) ) 
        do ii=1,size(thisPack)
            elemSR(thisPack(ii),esr_Weir_Zcrest) = &
                max( elemSR(thisPack(ii),esr_Weir_Zcrest), elemR(thisPack(ii),er_Zbottom) + setting%ZeroValue%Depth*1.d1  )
        end do
        deallocate(thisPack)

        !% --- initialize a default rectangular channel as the background of the weir
        call init_IC_diagnostic_default_geometry (thisLink)

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_IC_get_weir_geometry
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_get_orifice_geometry (thisLink)
        !%--------------------------------------------------------------------------
        !% Description:
        !% get the geometry and other data data for orifice links
        !%--------------------------------------------------------------------------
        !% Declarations
            integer, intent(in)  :: thisLink
            integer, pointer     :: specificOrificeType, OrificeGeometryType
            integer, allocatable :: thisPack(:)
            integer :: ii

            character(64) :: subroutine_name = 'init_IC_get_orifice_geometry'
        !--------------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% pointer to specific orifice type
        specificOrificeType => link%I(thisLink,li_link_sub_type)

        print *, 'CODE ERROR need orifice default geometry, similar to weir'
        call util_crashpoint(44982)

        select case (specificOrificeType)
        !% copy orifice specific data
        case (lBottomOrifice)
            where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                !% integer data
                elemSI(:,esi_Orifice_SpecificType)      = bottom_orifice
            endwhere
        case (lSideOrifice)
            where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                !% integer data
                elemSI(:,esi_Orifice_SpecificType)       = side_orifice
            endwhere
        case default
            print *, 'In ', subroutine_name
            print *, 'CODE ERROR: unknown orifice type, ', specificOrificeType,'  in network'
            print *, 'which has key ',trim(reverseKey(specificOrificeType))
            !stop 
            call util_crashpoint(8863411)
            !return
        end select

        !% pointer to specific orifice geometry
        OrificeGeometryType => link%I(thisLink,li_geometry)

        select case (OrificeGeometryType)
            !% copy orifice specific geometry data
        case (lRectangular_closed)  !% brh20211219 added Rect_closed
            where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                !% integer data
                elemI(:,ei_geometryType)          = rectangular_closed
                !% real data
                elemR(:,er_FullDepth)                    = link%R(thisLink,lr_FullDepth)
                elemSR(:,esr_Orifice_EffectiveFullDepth) = link%R(thisLink,lr_FullDepth)
                elemSR(:,esr_Orifice_DischargeCoeff)     = link%R(thisLink,lr_DischargeCoeff1)
                elemSR(:,esr_Orifice_Zcrest)             = elemR(:,er_Zbottom) + link%R(thisLink,lr_InletOffset)
                elemSR(:,esr_Orifice_Zcrown)             = elemSR(:,eSr_Orifice_Zcrest) + link%R(thisLink,lr_FullDepth)
                elemSR(:,esr_Orifice_RectangularBreadth) = link%R(thisLink,lr_BreadthScale)

                !% --- HACK the default background channel
                elemR(:,esgr_Rectangular_Breadth)   = twoR * elemSR(:,esr_Orifice_RectangularBreadth)
                elemR(:,er_BreadthMax)              = elemSGR(:,esgr_Rectangular_Breadth) 
                !% note -- full depth already used
            end where

        case (lCircular)
            where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                !% integer data
                elemI(:,ei_geometryType)    = circular
                !% real data
                elemSR(:,esr_Orifice_EffectiveFullDepth) = link%R(thisLink,lr_FullDepth)
                elemSR(:,esr_Orifice_DischargeCoeff)     = link%R(thisLink,lr_DischargeCoeff1)
                elemSR(:,esr_Orifice_Zcrest)             = elemR(:,er_Zbottom) + link%R(thisLink,lr_InletOffset)
                elemSR(:,esr_Orifice_Zcrown)             = elemSR(:,esr_Orifice_Zcrest) + link%R(thisLink,lr_FullDepth)

                !% --- HACK the default background channel
                elemR(:,esgr_Rectangular_Breadth)   = twoR * elemSR(:,esr_Orifice_EffectiveFullDepth)
                elemR(:,er_BreadthMax)              = elemSGR(:,esgr_Rectangular_Breadth) 
                elemR(:,er_FullDepth)               = elemSR(:,esr_Orifice_Zcrown) - elemR(:,er_Zbottom) &
                                                         + elemSR(:,esr_Orifice_EffectiveFullDepth)
                
            end where

        case default
            print *, 'In ', subroutine_name
            print *, 'CODE ERROR: unknown orifice geometry type, ', OrificeGeometryType,'  in network'
            print *, 'which has key ',trim(reverseKey(OrificeGeometryType))
            !stop 
            call util_crashpoint(8345553)
            !return
        end select

        !% --- set minimum crest height as 101% of the zero depth value for all orifices
        !%     this ensures that zero-height orifice elements cannot cause flow for zerovalue depths
        thisPack = pack(elemI(:,ei_Lidx),(elemI(:,ei_link_Gidx_BIPquick) == thisLink) ) 
        do ii=1,size(thisPack)
            elemSR(thisPack(ii),esr_Orifice_Zcrest) = &
                max( elemSR(thisPack(ii),esr_Orifice_Zcrest), elemR(thisPack(ii),er_Zbottom) + setting%ZeroValue%Depth*1.d1  )
        end do
        deallocate(thisPack)

        !% --- initialize a default rectangular channel as the background of the weir
        call init_IC_diagnostic_default_geometry (thisLink)


        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_orifice_geometry
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_get_pump_geometry (thisLink)
        !%-------------------------------------------------------------------
        !% Description:
        !% get the geometry for pumps
        !%-------------------------------------------------------------------
            integer             :: ii
            integer, intent(in) :: thisLink
            integer, pointer    :: specificPumpType, curveID, eIDx

            character(64) :: subroutine_name = 'init_IC_get_pump_geometry'
        !--------------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% pointer to specific pump type
        specificPumpType => link%I(thisLink,li_link_sub_type)
        curveID          => link%I(thisLink,li_curve_id)

        do ii = 1,N_elem(this_image())
            if (elemI(ii,ei_link_Gidx_BIPquick) == thisLink) then  
                !% real data
                elemSR(ii,esr_Pump_yOn)     = link%R(thisLink,lr_yOn)
                elemSR(ii,esr_Pump_yOff)    = link%R(thisLink,lr_yOff)

                if (curveID < zeroI) then
                    !% integer data
                    elemSI(ii,esi_Pump_SpecificType) = type_IdealPump 
                else
                    elemSI(ii,esi_Pump_CurveID) = curveID
                    Curve(curveID)%ElemIdx      = ii
                    !% copy pump specific data
                    if (specificPumpType == lType1Pump) then
                        !% integer data
                        elemSI(ii,esi_Pump_SpecificType) = type1_Pump

                    else if (specificPumpType == lType2Pump) then
                        !% integer data
                        elemSI(ii,esi_Pump_SpecificType) = type2_Pump

                    else if (specificPumpType == lType3Pump) then
                        !% integer data
                        elemSI(ii,esi_Pump_SpecificType) = type3_Pump

                    else if (specificPumpType == lType4Pump) then
                        !% integer data
                        elemSI(ii,esi_Pump_SpecificType) = type4_Pump

                    else if (specificPumpType == lTypeIdealPump) then
                        !% integer data
                        elemSI(ii,esi_Pump_SpecificType) = type_IdealPump 
                    else
                        print *, 'In ', subroutine_name
                        print *, 'CODE ERROR: unknown orifice type, ', specificPumpType,'  in network'
                        print *, 'which has key ',trim(reverseKey(specificPumpType))
                        !stop 
                        call util_crashpoint(8863411)
                        !return
                    end if
                end if 
            end if
        end do
     
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_pump_geometry
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_get_outlet_geometry (thisLink)
        !--------------------------------------------------------------------------
        !
        !% get the geometry and other data data for orifice links
        !
        !--------------------------------------------------------------------------
            integer             :: ii
            integer, intent(in) :: thisLink
            integer, pointer    :: specificOutletType, curveID, eIDx

            character(64) :: subroutine_name = 'init_IC_get_outlet_geometry'
        !--------------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% pointer to specific outlet type
        specificOutletType => link%I(thisLink,li_link_sub_type)
        curveID            => link%I(thisLink,li_curve_id)

        do ii = 1,N_elem(this_image())
            if (elemI(ii,ei_link_Gidx_BIPquick) == thisLink) then

                !% real data
                elemSR(ii,esr_Outlet_Coefficient) = link%R(thisLink,lr_DischargeCoeff1)
                elemSR(ii,esr_Outlet_Exponent)    = link%R(thisLink,lr_DischargeCoeff2)
                elemSR(ii,esr_Outlet_Zcrest)      = elemR(ii,er_Zbottom) + link%R(thisLink,lr_InletOffset)

                if ((specificOutletType == lNodeDepth) .and. (curveID == zeroI)) then
                    !% integer data
                    elemSI(ii,esi_Outlet_SpecificType)  = func_depth_outlet
                elseif ((specificOutletType == lNodeDepth) .and. (curveID /= zeroI)) then
                    !% integer data
                    elemSI(ii,esi_Outlet_SpecificType)  = tabl_depth_outlet
                    elemSI(ii,esi_Outlet_CurveID)       = curveID
                    Curve(curveID)%ElemIdx              = ii
                elseif ((specificOutletType == lNodeHead) .and. (curveID == zeroI)) then
                    !% integer data
                    elemSI(ii,esi_Outlet_SpecificType)  = func_head_outlet
                elseif ((specificOutletType == lNodeHead) .and. (curveID /= zeroI)) then
                    !% integer data
                    elemSI(ii,esi_Outlet_SpecificType)  = tabl_head_outlet
                    elemSI(ii,esi_Outlet_CurveID)       = curveID
                    Curve(curveID)%ElemIdx              = ii
                else
                    print*, 'In ', subroutine_name
                    print*, 'CODE ERROR: unknown outlet type, ', specificOutletType,'  in network'
                    print *, 'which has key ',trim(reverseKey(specificOutletType))
                    !stop 
                    call util_crashpoint(82564)
                    !return
                end if
            end if 
        end do

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_outlet_geometry
!%
!%=========================================================================
!%=========================================================================
!%
    subroutine init_IC_diagnostic_default_geometry (thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% Provides default background channel geometry for diagnostic
        !% elements. This geometry will be overwritten by that of an
        !% adjacent cell that has valid channel/conduit geometry
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisLink
        !%-----------------------------------------------------------------

        where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
            elemR(:,er_FullHydDepth)            = elemR(:,er_FullDepth)
            elemR(:,er_FullPerimeter)           = elemR(:,er_BreadthMax) + twoR * elemR(:,er_FullDepth)
            elemR(:,er_ZbreadthMax)             = elemR(:,er_FullDepth)
            elemR(:,er_Zcrown)                  = elemR(:,er_Zbottom) + elemR(:,er_FullDepth)
            elemR(:,er_FullArea)                = elemR(:,er_FullDepth) * elemR(:,er_BreadthMax)
            elemR(:,er_FullVolume)              = elemR(:,er_FullArea)  * elemR(:,er_Length)
            elemR(:,er_AreaBelowBreadthMax)     = elemR(:,er_FullArea)
            elemR(:,er_ell_max)                  = elemR(:,er_FullDepth)

            !% store IC data
            elemR(:,er_Area)          = elemSGR(:,esgr_Rectangular_Breadth) * elemR(:,er_Depth)
            elemR(:,er_Area_N0)       = elemR(:,er_Area)
            elemR(:,er_Area_N1)       = elemR(:,er_Area)
            elemR(:,er_Volume)        = elemR(:,er_Area) * elemR(:,er_Length)
            elemR(:,er_Volume_N0)     = elemR(:,er_Volume)
            elemR(:,er_Volume_N1)     = elemR(:,er_Volume)
        endwhere
     
    end subroutine init_IC_diagnostic_default_geometry
!%
!%=========================================================================
!%=========================================================================
!%
    subroutine init_IC_diagnostic_geometry_from_adjacent ()
        !%-----------------------------------------------------------------
        !% Description:   20220511brh
        !% Provides the additional "background" geometry of
        !% diagnostic elements based on their surroundings. This is the
        !% geometry of the channel/conduit in which the weir, orifice,
        !% pump, or outlet exists.  This ensures that a weir (etc) next to
        !% a JB branch has a valid geometry that can be used for the JB branch.
        !% Weir -- if upstream element is CC, the weir takes on the
        !%   geometry of the upstream CC element. If the upstream element is
        !%   other than CC, then the weir takes on the geometry of the
        !%   downstream CC element. If the downstream element is also other than 
        !%   CC then the default geometry from (e.g. init_IC_get_weir_geometry)
        !%   is used
        !% Orifice -- same as weir
        !% Pump --- tsame as weir
        !% Outlet --- to be determined
        !%-----------------------------------------------------------------
        !% Declarations
            integer, dimension(:), allocatable, target :: packIdx
            integer, pointer :: Fidx, Aidx, thisP
            integer, pointer :: linkIdx
            integer :: ii, jj, Ci
            integer :: thisCol(16)

            character(64) :: subroutine_name = 'init_IC_diagnostic_geometry'
        !%-----------------------------------------------------------------
        !% Preliminaries:
            !% --- get the set of weirs, orifices, pumps, and outlets
            packIdx = pack(elemI(:,ei_Lidx), &
                    ( (elemI(:,ei_elementType) .eq. weir) &
                    .or. &
                    (elemI(:,ei_elementType) .eq. orifice) &
                    .or. &
                    (elemI(:,ei_elementType) .eq. pump) &
                    .or. &
                    (elemI(:,ei_elementType) .eq. outlet) ) )

            !% --- set the column indexes of the standard geometry data that
            !%     are needed in a diagnostic element (i.e., copied to the JB branch)
            thisCol(1) = er_BreadthMax
            thisCol(2) = er_FullDepth
            thisCol(3) = er_ZbreadthMax
            thisCol(4) = er_Zcrown
            thisCol(5) = er_FullArea
            thisCol(6) = er_FullVolume
            thisCol(7) = er_AreaBelowBreadthMax
            thisCol(8) = er_ell_max
            thisCol(9) = er_FullHydDepth
            thisCol(10)= er_FullPerimeter
            thisCol(11)= er_Area
            thisCol(12)= er_Area_N0        
            thisCol(13)= er_Area_N1
            thisCol(14)= er_Volume
            thisCol(15)= er_Volume_N0
            thisCol(16)= er_Volume_N1
        !%-----------------------------------------------------------------
        !% --- cycle through to set geometry of diagnostic element
        !%     use the upstream geometry if it is CC
        do ii=1,size(packIdx)
            !% --- the present point
            thisP  => packIdx(ii)
            !% --- the link
            linkIdx => elemI(thisP,ei_link_Gidx_SWMM)
            !% --- the upstream face
            Fidx => elemI(thisP,ei_Mface_uL)
            !% --- the upstream element
            !%     which may be on a different image
            if (elemYN(thisP,eYN_isBoundary_up)) then
                Ci   =  faceI(Fidx,fi_Connected_image)
                Aidx => faceI(Fidx,fi_GhostElem_uL)
            else
                Ci   =  this_image()
                Aidx => faceI(Fidx,fi_Melem_uL)
            end if

            !% --- the element type upstream
            select case (elemI(Aidx,ei_elementType))
                case (CC)
                    !% --- if an upstream element is a channel/conduit, use this for the background channel
                    !%     geometry of the diagnostic element in which the weir/orifice/pump/outlet is embeded
                    elemI(thisP,ei_geometryType) = elemI(Aidx,ei_geometryType)[Ci]
                    elemI(thisP,thisCol)         = elemI(Aidx,thisCol)[Ci]
                    !% --- copy special geometry
                    call init_IC_diagnostic_special_geometry (thisP, Aidx, Ci)
                
                case default
                    !% --- if the upstream element is not CC, then use the downstream element geometry
                    !%     (assuming it is CC)
                    !% --- the downstream face
                    Fidx => elemI(thisP,ei_Mface_dL)
                    !% --- the downstream element
                    !%     which may be on a different image
                    if (elemYN(thisP,eYN_isBoundary_dn)) then
                        Ci   =  faceI(Fidx,fi_Connected_image)
                        Aidx => faceI(Fidx,fi_GhostElem_dL)
                    else
                        Ci   =  this_image()
                        Aidx => faceI(Fidx,fi_Melem_dL)
                    end if

                    !% --- the element type downstream
                    select case (elemI(Aidx,ei_elementType))
                        case (CC)
                            !% --- if a downstream element is a channel/conduit, use this for the
                            !%     geometry of the element in which the weir/orifice/pump is embeded
                            elemI(thisP,ei_geometryType) = elemI(Aidx,ei_geometryType)[Ci]
                            elemI(thisP,thisCol)         = elemI(Aidx,thisCol)[Ci]
                            !% --- copy over special geometry
                            call init_IC_diagnostic_special_geometry (thisP, Aidx, Ci)
                        
                        case default
                            !% --- if the downstream element is NOT channel/conduit, then accept
                            !%     the default geometry or fail
                            select case (elemI(thisP,ei_elementType))
                                case (weir)
                                    !% --- weirs are OK with default channel geometry, so no action
                                case (orifice)
                                    !% --- orifices are OK with default channel geometry, so no action
                                case (pump)
                                    !% --- pumps do not have default channel geometry, so they must
                                    !%     have a CC element upstream or downstream.
                                    print *, 'USER SYSTEM CONFIGURATION ERROR'
                                    print *, 'A pump requires at least one upstream or downstream link that is a'
                                    print *, 'conduit or channel. This condition violated for'
                                    print *, 'pump with name ',trim(link%Names(linkIdx)%str)
                                    call util_crashpoint(23987)
                                case (outlet)
                                    !% --- outlets do not have default channel geometry, so they must
                                    !%     have a CC element upstream or downstream.
                                    print *, 'USER SYSTEM CONFIGURATION ERROR'
                                    print *, 'An outelt requires at least one upstream or downstream link that is a'
                                    print *, 'conduit or channel. This condition violated for'
                                    print *, 'outlet with name ',trim(link%Names(linkIdx)%str)
                                    call util_crashpoint(92873)
                            end select
                    end select
            end select
        end do

        !%-----------------------------------------------------------------
        !% Closing:
            deallocate(packIdx)

    end subroutine init_IC_diagnostic_geometry_from_adjacent
!%
!%==========================================================================
!%==========================================================================
!% 
    subroutine init_IC_diagnostic_special_geometry (thisP, Aidx, Ci)
        !%-----------------------------------------------------------------
        !% Description:
        !% Copies the special fixed geoemtry (depends on element geometry type)
        !% from the adjacent cell (Aidx) to this cell (thisP) where
        !% Aidx is on the connected image (Ci). This is used to get the
        !% geometry for a JB junction branch
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisP, Aidx, Ci
            character(64) :: subroutine_name = 'init_IC_diagnostic_special_geometry'
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------
        !% --- copy over special geometry data depending on geometry type
        select case (elemI(thisP,ei_geometryType))
            case (rectangular)
                elemSGR(thisP,esgr_Rectangular_Breadth)    = elemSGR(Aidx,esgr_Rectangular_Breadth)[Ci]
            case (trapezoidal)
                elemSGR(thisP,esgr_Trapezoidal_Breadth)    = elemSGR(Aidx,esgr_Trapezoidal_Breadth)[Ci]
                elemSGR(thisP,esgr_Trapezoidal_LeftSlope)  = elemSGR(Aidx,esgr_Trapezoidal_LeftSlope)[Ci]
                elemSGR(thisP,esgr_Trapezoidal_RightSlope) = elemSGR(Aidx,esgr_Trapezoidal_RightSlope)[Ci]
            case (triangular)
                elemSGR(thisP,esgr_Triangular_TopBreadth)  = elemSGR(Aidx,esgr_Triangular_TopBreadth)[Ci]
                elemSGR(thisP,esgr_Triangular_Slope)       = elemSGR(Aidx,esgr_Triangular_Slope)[Ci] 
            case (irregular)
                elemI(thisP,ei_link_transect_idx)          = elemI(Aidx,ei_link_transect_idx)[Ci]
            case (rectangular_closed)
                elemSGR(thisP,esgr_Rectangular_Breadth)    = elemSGR(Aidx,esgr_Rectangular_Breadth)[Ci]
            case (circular)
                elemSGR(thisP,esgr_Circular_Diameter)      = elemSGR(Aidx,esgr_Circular_Diameter)[Ci]
                elemSGR(thisP,esgr_Circular_Radius)        = elemSGR(Aidx,esgr_Circular_Radius)[Ci]
            case default
                print *, 'CODE ERROR unexpected geometry'
                print *, 'ei_geometryType index # ',elemI(thisP,ei_geometryType)
                print *, 'which represents ',reverseKey(elemI(thisP,ei_geometryType))
                print *, 'is not handled in subroutine ',trim(subroutine_name)
                call util_crashpoint(99376)
        end select
                
    end subroutine init_IC_diagnostic_special_geometry
!%
!%==========================================================================
!%==========================================================================
!%    
    ! subroutine init_IC_get_channel_conduit_velocity (thisLink, ePack, npack)
    !     !% brh 20211216 obsolete -- replaced with init_IC_derived_values()
    !     !%-----------------------------------------------------------------
    !     !% Description:
    !     !% get the velocity of channel and conduits
    !     !% and sell all other velocity to zero
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: thisLink
    !         integer, pointer    :: specificWeirType
    !         character(64) :: subroutine_name = 'init_IC_get_channel_conduit_velocity'
    !     !%------------------------------------------------------------------
    !     !% Preliminaries:
    !         !if (crashYN) return
    !         if (setting%Debug%File%initial_condition) &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%------------------------------------------------------------------

    !     !% HACK: this might not be right
    !     where ( (elemI(:,ei_link_Gidx_BIPquick) == thisLink) .and. &
    !             (elemR(:,er_area)               .gt. zeroR ) .and. &
    !             (elemI(:,ei_elementType)        == CC      ) )

    !         elemR(:,er_Velocity)    = elemR(:,er_Flowrate) / elemR(:,er_Area)
    !         elemR(:,er_Velocity_N0) = elemR(:,er_Velocity)
    !         elemR(:,er_Velocity_N1) = elemR(:,er_Velocity)

    !     elsewhere ( (elemI(:,ei_link_Gidx_BIPquick) == thisLink) .and. &
    !                 (elemR(:,er_area)               .le. zeroR ) .and. &
    !                 (elemI(:,ei_elementType)        == CC    ) )

    !         elemR(:,er_Velocity)    = zeroR
    !         elemR(:,er_Velocity_N0) = zeroR
    !         elemR(:,er_Velocity_N1) = zeroR

    !     endwhere

    !     if (setting%Debug%File%initial_condition) &
    !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    ! end subroutine init_IC_get_channel_conduit_velocity
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_from_nodedata ()
        !--------------------------------------------------------------------------
        !% get the initial depth, and geometry data from nJm nodes
        !--------------------------------------------------------------------------

            integer                       :: ii, image, pJunction
            integer, pointer              :: thisJunctionNode
            integer, allocatable, target  :: packed_nJm_idx(:)

            character(64) :: subroutine_name = 'init_IC_from_nodedata'
        !--------------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% Setting the local image value
        image = this_image()

        !% pack all the link indexes in an image
        packed_nJm_idx = pack(node%I(:,ni_idx), &
                             ((node%I(:,ni_P_image) == image) .and. &
                              (node%I(:,ni_node_type) == nJm) ) )

        !% find the number of links in an image
        pJunction = size(packed_nJm_idx)

        !% cycle through the links in an image
        do ii = 1,pJunction
            !% necessary pointers
            thisJunctionNode => packed_nJm_idx(ii)
            call init_IC_get_junction_data (thisJunctionNode)
        end do

        !% deallocate the temporary array
        deallocate(packed_nJm_idx)

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_from_nodedata
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_get_junction_data (thisJunctionNode)        
        !--------------------------------------------------------------------------
        !% get data for the multi branch junction elements
        !--------------------------------------------------------------------------
        integer, intent(in) :: thisJunctionNode

        integer              :: ii, jj, JMidx, JBidx, Aidx, Ci
        integer, pointer     :: BranchIdx, JBgeometryType, JmType, curveID, NumRows
        integer, pointer     :: Fidx, F2idx
        integer              :: nbranches
        real(8), allocatable :: integrated_volume(:)

        character(64) :: subroutine_name = 'init_IC_get_junction_data'
        !--------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !print *, 'inside ',trim(subroutine_name)

        !%................................................................
        !% Junction main
        !%................................................................
        !% find the first element ID associated with that nJm
        !% masked on the global node number for this node.
        JMidx = minval(elemI(:,ei_Lidx), elemI(:,ei_node_Gidx_SWMM) == thisJunctionNode)

        !print *, 'JMidx ',JMidx

        !% the first element index is a junction main
        elemI(JMidx,ei_elementType)  = JM
        elemI(JMidx,ei_HeqType)      = time_march
        elemI(JMidx,ei_QeqType)      = notused

        !% set the type of junction main
        if (node%YN(thisJunctionNode,nYN_has_storage)) then
            if (node%I(thisJunctionNode,ni_curve_ID) .eq. 0) then
                elemSI(JMidx,esi_JunctionMain_Type)   = FunctionalStorage
                elemSR(JMidx,esr_Storage_Constant)    = node%R(thisJunctionNode,nr_StorageConstant)
                elemSR(JMidx,esr_Storage_Coefficient) = node%R(thisJunctionNode,nr_StorageCoeff)
                elemSR(JMidx,esr_Storage_Exponent)    = node%R(thisJunctionNode,nr_StorageExponent)

            else
                elemSI(JMidx,esi_JunctionMain_Type) = TabularStorage
                elemSI(JMidx,esi_JunctionMain_Curve_ID) = node%I(thisJunctionNode,ni_curve_ID)
            end if
        else
            !%-----------------------------------------------------------------------
            !% HACK: Junction main with implied storage are rectangular
            !%-----------------------------------------------------------------------
            elemSI(JMidx,esi_JunctionMain_Type) = ImpliedStorage
            elemI(JMidx,ei_geometryType)        = rectangular
        end if

        !% junction main depth and head from initial conditions
        elemR(JMidx,er_Depth)     = node%R(thisJunctionNode,nr_InitialDepth)
        elemR(JMidx,er_Head)      = elemR(JMidx,er_Depth) + elemR(JMidx,er_Zbottom)
        elemR(JMidx,er_FullDepth) = node%R(thisJunctionNode,nr_FullDepth)

        !% JM elements are not solved for momentum.
        elemR(JMidx,er_Flowrate)     = zeroR
        elemR(JMidx,er_Velocity)     = zeroR

        !% wave speed is the gravity wave speed for the depth
        elemR(JMidx,er_WaveSpeed)    = sqrt(setting%constant%gravity * elemR(JMidx,er_Depth))
        elemR(JMidx,er_FroudeNumber) = zeroR

        !% find if the node can surcharge
        if (node%R(thisJunctionNode,nr_SurchargeDepth) .ne. nullValueR) then
            elemYN(JMidx,eYN_canSurcharge)  = .true.
            elemR(JMidx,er_FullDepth)       = node%R(thisJunctionNode,nr_SurchargeDepth)
            ! elemI(JMidx,ei_geometryType)    = rectangular_closed
        else
            elemYN(JMidx,eYN_canSurcharge)  = .false.
        end if

        !%................................................................
        !% Junction Branches
        !%................................................................
        !% loopthrough all the branches
        !% HACK -- replace much of this with a call to the standard geometry after all other IC have
        !% been done. The branch depth should be based on the upstream or downstream depth of the
        !% adjacent element.
        do ii = 1,max_branch_per_node

            !% 20220406brh Rewritten to use adjacent element geometry initialization where possible.

            !% find the element id of junction branches
            JBidx = JMidx + ii

            !print *, 'JBidx ',JBidx

            !% set the adjacent element id where elemI and elemR data can be extracted
            !% note that elemSGR etc are not yet assigned
            if (mod(ii,2) == zeroI) then
                Fidx => elemI(JBidx,ei_MFace_dL)
                !% even are downstream branches
                if (elemYN(JBidx,eYN_isBoundary_dn)) then
                    Ci   = faceI(Fidx,fi_Connected_image)
                    Aidx = faceI(Fidx,fi_GhostElem_dL)
                else
                    Ci   = this_image()
                    Aidx = faceI(Fidx,fi_Melem_dL)
                end if
            !% odd are upstream branches
            else
                Fidx => elemI(JBidx,ei_MFace_uL)
                if (elemYN(JBidx,eYN_isBoundary_up)) then
                    Ci   = faceI(Fidx,fi_Connected_image)
                    Aidx = faceI(Fidx,fi_GhostElem_uL)
                else
                    Ci   = this_image()
                    Aidx = faceI(Fidx,fi_Melem_uL)
                end if
            end if

            !% set the junction branch element type
            elemI(JBidx,ei_elementType) = JB

            !% set the geometry for existing branches
            !% Note that elemSI(,...Exists) is set in init_network_handle_nJm
            if (.not. elemSI(JBidx,esi_JunctionBranch_Exists) == oneI) cycle

            BranchIdx      => elemSI(JBidx,esi_JunctionBranch_Link_Connection)
            JBgeometryType => link%I(BranchIdx,li_geometry)

            !print *, 'linkgeo', link%I(BranchIdx,li_geometry), reverseKey(link%I(BranchIdx,li_geometry))

            !% set the JB to time_march for use with splitting between AC
            !% and ETM in rk2_extrapolate_to_fullstep_ETM, rk2_restore_to_midstep_ETM
            !% rk2_interpolate_to_halfstep_AC, k2_restore_to_fullstep_AC
            elemI(JBidx,ei_HeqType) = time_march
            elemI(JBidx,ei_QeqType) = time_march

            !% Junction branch k-factor
            !% If the user does not input the K-factor for junction branches entrance/exit loses then
            !% use default from setting
            if (node%R(thisJunctionNode,nr_JunctionBranch_Kfactor) .ne. nullvalueR) then
                elemSR(JBidx,esr_JunctionBranch_Kfactor) = node%R(thisJunctionNode,nr_JunctionBranch_Kfactor)
            else
                elemSR(JBidx,esr_JunctionBranch_Kfactor) = setting%Junction%kFactor
            end if

            !% set the initial head to the same as the junction main
            elemR(JBidx,er_Head)    = elemR(JMidx,er_Head)
            elemR(JBidx,er_Depth)   = elemR(JBidx,er_Head) - elemR(JBidx,er_Zbottom)
            if (elemR(JBidx,er_Head) < elemR(JBidx,er_Zbottom)) then
                elemR(JBidx,er_Head) = elemR(JBidx,er_Zbottom)
                elemR(JBidx,er_Depth) = zeroR
            end if

            elemR(JBidx,er_VolumeOverFlow) = zeroR
            elemR(JBidx,er_VolumeOverFlowTotal) = zeroR

            !% JB elements initialized for momentum
            elemR(JBidx,er_Flowrate)     = elemR(Aidx,er_Flowrate)[Ci] !% flowrate of adjacent element
            elemR(JBidx,er_WaveSpeed)    = sqrt(setting%constant%gravity * elemR(JBidx,er_Depth))
            elemR(JBidx,er_FroudeNumber) = zeroR

            !% Set the face flowrates such that it does not blowup the initial interpolation
            if (elemI(JBidx, ei_Mface_uL) /= nullvalueI) then
                faceR(elemI(JBidx, ei_Mface_uL),fr_flowrate) = elemR(JBidx,er_Flowrate) 
            else if (elemI(JBidx, ei_Mface_dL) /= nullvalueI) then
                faceR(elemI(JBidx, ei_Mface_dL),fr_flowrate) = elemR(JBidx,er_Flowrate)
            end if

            !% Set the geometry from the adjacent elements
            elemI(JBidx,ei_geometryType)        = elemI(Aidx,ei_geometryType)[Ci]
            elemYN(JBidx,eYN_canSurcharge)      = elemYN(Aidx,eYN_canSurcharge)[Ci]

            !print *, 'JBidx',JBidx, elemI(JBidx,ei_geometryType)
            !print *, 'Aidx ',Aidx, elemI(Aidx,ei_geometryType)

            select case  (elemI(JBidx,ei_geometryType))

            case (rectangular, trapezoidal, triangular, rectangular_closed, circular, irregular)
                !% --- Copy all the geometry specific data from the adjacent element cell
                !%     Note that because irregular transect tables are not yet initialized, the
                !%     Area and Volume here will be junk for an irregular cross-section and will need to be
                !%     reset after transect tables are initialized. This occurs because we have
                !%     to cycle through all the CC, JM/JB before we can set the element transect
                !%     tables.
                elemR(JBidx,er_Area)                = elemR(Aidx,er_Area)[Ci]
                elemR(JBidx,er_AreaBelowBreadthMax) = elemR(Aidx,er_AreaBelowBreadthMax)[Ci]
                elemR(JBidx,er_BreadthMax)          = elemR(Aidx,er_BreadthMax)[Ci]
                elemR(JBidx,er_FullArea)            = elemR(Aidx,er_FullArea)[Ci]
                elemR(JBidx,er_FullDepth)           = elemR(Aidx,er_FullDepth)[Ci]
                elemR(JBidx,er_FullHydDepth)        = elemR(Aidx,er_FullHydDepth)[Ci]
                elemR(JBidx,er_FullPerimeter)       = elemR(Aidx,er_FullPerimeter)[Ci]
                elemR(JBidx,er_ZbreadthMax)         = elemR(Aidx,er_ZbreadthMax)[Ci]
                elemR(JBidx,er_Zcrown)              = elemR(Aidx,er_Zcrown)[Ci]         
                elemR(JBidx,er_Roughness)           = elemR(Aidx,er_Roughness)[Ci]
                elemI(JBidx,ei_link_transect_idx)   = elemI(Aidx,ei_link_transect_idx)[Ci]
                !% copy the entire row of the elemSGR array
                elemSGR(JBidx,:)                    = elemSGR(Aidx,:)[Ci]

            case (undefinedKey)
                print *, 'in ',trim(subroutine_name)
                print *, 'CODE ERROR: undefinedKey for ei_geometryType for junction'
                print *, 'at JBidx ',JBidx
                print * , ' '
                call util_crashpoint (23374)
                !return

            case default
                print *, 'in ',trim(subroutine_name)
                print *, 'CODE ERROR: unknown geometry type ',elemI(JBidx,ei_geometryType)
                print *, 'which has key ',trim(reverseKey(elemI(JBidx,ei_geometryType)))
                call util_crashpoint (4473)
                !return
            end select

            !% set the velocity
            if (elemR(JBidx,er_Area) .gt. setting%ZeroValue%Area) then ! BRHbugfix 20210813
                elemR(JBidx,er_Velocity) = elemR(JBidx,er_Flowrate) / elemR(JBidx,er_Area)
            else
                elemR(JBidx,er_Velocity) = zeroR
            end if

            !% Common geometry that do not depend on cross-section
            elemR(JBidx,er_Area_N0)      = elemR(JBidx,er_Area)
            elemR(JBidx,er_Area_N1)      = elemR(JBidx,er_Area)
            elemR(JBidx,er_FullVolume)   = elemR(JBidx,er_FullArea)  * elemR(JBidx,er_Length)
            elemR(JBidx,er_Volume)       = elemR(JBidx,er_Area) * elemR(JBidx,er_Length)
            elemR(JBidx,er_Volume_N0)    = elemR(JBidx,er_Volume)
            elemR(JBidx,er_Volume_N1)    = elemR(JBidx,er_Volume)

        end do
        
        !% get junction main geometry based on type
        JmType => elemSI(JMidx,esi_JunctionMain_Type)

        select case (JmType)

        case (ImpliedStorage)
            !% the JM characteristic length is the sum of the two longest branches
            elemR(JMidx,er_Length) = max(elemR(JMidx+1,er_Length), elemR(JMidx+3,er_Length), &
                                            elemR(JMidx+5,er_Length)) + &
                                        max(elemR(JMidx+2,er_Length), elemR(JMidx+4,er_Length), &
                                            elemR(JMidx+6,er_Length))

            !% --- Plane area is the sum of the branch plane area 
            !%     This uses simplified geometry approximations as the junction main is only
            !%     mass conservation only, which means its volume change can be approximated
            !%     as if it is a rectangular box of Storage_Plane_Area x Depth
            elemSR(JMidx,esr_Storage_Plane_Area) = zeroR

            do ii=1,max_branch_per_node
                JBidx = JMidx + ii
                if (.not. elemSI(JBidx,esi_JunctionBranch_Exists) == oneI) cycle

                BranchIdx      => elemSI(JBidx,esi_JunctionBranch_Link_Connection)
                JBgeometryType => link%I(BranchIdx,li_geometry)

                select case (JBgeometryType)
                case (lRectangular,lRectangular_closed)
                    elemSR(JMidx,esr_Storage_Plane_Area) = elemSR(JMidx,esr_Storage_Plane_Area)  &
                     +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)                       &
                          * elemR(  JBidx,er_Length)                                          &
                          * elemSGR(JBidx,esgr_Rectangular_Breadth) )

                case (lTrapezoidal)
                    elemSR(JMidx,esr_Storage_Plane_Area) = elemSR(JMidx,esr_Storage_Plane_Area)  &
                             +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)               &
                                  * elemR(  JBidx,er_Length)                                  &
                                  * elemSGR(JBidx,esgr_Trapezoidal_Breadth) )

                case (lTriangular)
                    elemSR(JMidx,esr_Storage_Plane_Area) = elemSR(JMidx,esr_Storage_Plane_Area)  &
                                +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)               &
                                    * elemR(  JBidx,er_Length)                                   &
                                    * (elemSGR(JBidx,esgr_Triangular_TopBreadth)/twoR) )

                case (lCircular)
                    elemSR(JMidx,esr_Storage_Plane_Area) = elemSR(JMidx,esr_Storage_Plane_Area)  &
                     +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)                       &
                          * elemR(  JBidx,er_Length)                                          &
                          * elemSGR(JBidx,esgr_Circular_Diameter) )

                case (lIrregular)    
                    !% --- for irregular geometry, we use the average breadth for the area below
                    !%     the maximum breadth as the characteristic width of the branch
                    elemSR(JMidx,esr_Storage_Plane_Area) = elemSR(JMidx,esr_Storage_Plane_Area)  &
                     +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)                       &
                          * elemR(  JBidx,er_Length)                                          &
                          * elemR(  JBidx,er_AreaBelowBreadthMax) / elemR(JBidx,er_ZbreadthMax))

                case default
                    print *, 'In, ', subroutine_name
                    print *, 'CODE ERROR: geometry type unknown for # ', JBgeometryType
                    print *, 'which has key ',trim(reverseKey(JBgeometryType))
                    !stop 
                    call util_crashpoint(308974)
                    !return
                end select
            end do

            !% Breadth is consistent with length and plane area
            elemSGR(JMidx,esgr_Rectangular_Breadth) =  elemSR(JMidx,esr_Storage_Plane_Area) &
                                                    /   elemR(JMidx,er_Length)

            !% Volume depends on plane area and depth
            elemR(JMidx,er_Volume)     = elemSR(JMidx,esr_Storage_Plane_Area) * elemR(JMidx,er_Depth)
            elemR(JMidx,er_Volume_N0)  = elemR(JMidx,er_Volume)
            elemR(JMidx,er_Volume_N1)  = elemR(JMidx,er_Volume)
            elemR(JMidx,er_FullVolume) = elemSR(JMidx,esr_Storage_Plane_Area) * elemR(JMidx,er_FullDepth)

        case (FunctionalStorage)
            elemR(JMidx,er_Volume)     = elemSR(JMidx,esr_Storage_Constant) * elemR(JMidx,er_Depth)      &
                + (elemSR(JMidx,esr_Storage_Coefficient) / (elemSR(JMidx,esr_Storage_Exponent) + oneR))  &
                    * elemR(JMidx,er_Depth) ** (elemSR(JMidx,esr_Storage_Exponent) + oneR)
            elemR(JMidx,er_Volume_N0)  = elemR(JMidx,er_Volume)
            elemR(JMidx,er_Volume_N1)  = elemR(JMidx,er_Volume)
            elemR(JMidx,er_FullVolume) = elemSR(JMidx,esr_Storage_Constant) * elemR(JMidx,er_FullDepth)  &
                + (elemSR(JMidx,esr_Storage_Coefficient) / (elemSR(JMidx,esr_Storage_Exponent) + oneR))  &
                    * elemR(JMidx,er_FullDepth) ** (elemSR(JMidx,esr_Storage_Exponent) + oneR)
            !% create a storage curve
            call storage_create_curve (JMidx)

        case (TabularStorage)
            CurveID => elemSI(JMidx,esi_JunctionMain_Curve_ID)
            NumRows => curve(CurveID)%NumRows 
            Curve(CurveID)%ElemIdx = JMidx
            !% SWMM5+ needs a volume vs depth relationship thus Trapezoidal rule is used
            !% to get to integrate the area vs depth curve
            call storage_integrate_area_vs_depth_curve (CurveID)

            !% now interpolate from the cure to get the volume
            call storage_interpolate_volume_from_depth_singular (JMidx)

            elemR(JMidx,er_Volume_N0)  = elemR(JMidx,er_Volume)
            elemR(JMidx,er_Volume_N1)  = elemR(JMidx,er_Volume)
            elemR(JMidx,er_FullVolume) = Curve(CurveID)%ValueArray(NumRows,curve_storage_volume)
        case default
            !% IMPORTANT -- if any other new type is defined, make sure that
            !% subroutine geo_depth_from_volume is updated
            print *, 'In, ', subroutine_name
            print *, 'CODE ERROR junction main type unknown for # ', JmType
            print *, 'which has key ',trim(reverseKey(JmType))
            !stop 
            call util_crashpoint(54895)
            !return

        end select

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_junction_data   
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_elem_transect_arrays ()
        !%------------------------------------------------------------------
        !% Description:
        !% initializes the transect tables for elements (as opposed to the
        !% link transect tables, which are already stored). Note that this
        !% must be done after both CC and JM/JB element geometry are 
        !% initialized so that we have a count of all the elements with
        !% irregular geometry.
        !% This presumes that the ei_link_transect_idx has already been 
        !% assigned
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: geometryType(:)
            integer, pointer :: elemTransectIdx(:), linkTransectIdx(:)
            integer          :: ii, thisTransectIdx

            character(64) :: subroutine_name = 'init_IC_elem_transect_arrays'
        !%------------------------------------------------------------------
        !% Aliases:
            geometryType    => elemI(:,ei_geometryType)
            elemTransectIdx => elemI(:,ei_transect_idx)
            linkTransectIdx => elemI(:,ei_link_transect_idx)
        !%------------------------------------------------------------------
        !% Preliminaries:
            !% --- count the total number of elements with irregular cross-sections
            N_transect = count(elemI(:,ei_geometryType) .eq. irregular)
            if (N_transect .le. 0) return
        !%------------------------------------------------------------------

        !% --- allocate the element transect arrays
        call util_allocate_element_transect ()

        thisTransectIdx = 0
        !% --- cycle through the elements
        do ii=1,N_elem(this_image())

            if (geometryType(ii) .ne. irregular) cycle
            !% --- increment the element transect index
            thisTransectIdx = thisTransectIdx +1
            !% --- store the element transect index
            elemTransectIdx(ii) = thisTransectIdx

            !% --- store the link transect data for each element
            transectTableDepthR(thisTransectIdx,:,:) = link%transectTableDepthR(linkTransectIdx(ii),:,:)
            transectTableAreaR (thisTransectIdx,:,:) = link%transectTableAreaR (linkTransectIdx(ii),:,:)
            transectI(thisTransectIdx,:)             = link%transectI          (linkTransectIdx(ii),:)
            transectR(thisTransectIdx,:)             = link%transectR          (linkTransectIdx(ii),:)
            transectID(thisTransectIdx)              = link%transectID         (linkTransectIdx(ii))%str
        end do
        
        !% --- HACK: we should develop an approach to allow smoothing of irregular geometry
        !%     tables where 2 links connect and there is an abrupt change of geometry that isn't
        !%     really supported by the data. This should be an optional algorithm that 
        !%     looks at the size of the adjacent links and smooths some fraction of the
        !%     elements on either side of the transition. Note that this should only bee
        !%     applied across a nj2 node (face) that joins exactly 2 links without storage.

        !%------------------------------------------------------------------
    end subroutine init_IC_elem_transect_arrays
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine init_IC_elem_transect_geometry ()
        !%------------------------------------------------------------------
        !% Description:
        !% initializes the time=0 area, volume etc. that depend on the 
        !% initial depth for irregular cross-sections. 
        !% This is delayed from all the other geometry initialization because
        !% the transect tables must be setup prior to computing initial
        !% geometry
        !%------------------------------------------------------------------
        !% Declarations
            integer, allocatable :: tpack(:)
            integer :: npack, ii
            real(8), pointer :: area(:), area0(:), area1(:)
            real(8), pointer :: depth(:), fulldepth(:), ell(:), length(:)
            real(8), pointer :: hydDepth(:), hydRadius(:), perimeter(:)
            real(8), pointer :: topwidth(:)
            real(8), pointer :: volume(:), volume0(:), volume1(:)
            real(8), pointer :: thisTable(:,:)
        !%------------------------------------------------------------------
        !% Aliases:
            area      => elemR(:,er_Area)
            area0     => elemR(:,er_Area_N0)
            area1     => elemR(:,er_Area_N1)
            depth     => elemR(:,er_Depth)
            fulldepth => elemR(:,er_FullDepth)
            ell       => elemR(:,er_ell)
            hydDepth  => elemR(:,er_HydDepth)
            hydRadius => elemR(:,er_HydRadius)
            length    => elemR(:,er_Length)
            perimeter => elemR(:,er_Perimeter)
            topwidth  => elemR(:,er_Topwidth)
            volume    => elemR(:,er_Volume)
            volume0   => elemR(:,er_Volume_N0)
            volume1   => elemR(:,er_Volume_N1)
        !%------------------------------------------------------------------
        !% Preliminaries
            npack = count(elemI(:,ei_geometryType) .eq. irregular)
            if (npack .eq. 0) return
        !%------------------------------------------------------------------
        !% --- get the packed irregular element list
        tpack = pack(elemI(:,ei_Lidx), elemI(:,ei_geometryType) .eq. irregular)

        !print *, depth(tpack)
        !print *, trim(reverseKey(elemI(tpack(1),ei_geometryType)))
        !stop 29837

        !% --- temporary store of the normalized depth
        depth(tpack) = depth(tpack) / fulldepth(tpack)

        !% --- get area from normalized depth
        thisTable => transectTableDepthR(:,:,tt_area)
        call xsect_table_lookup_array (area, depth, thisTable, tpack) 

        !% --- get hydraulic radius from normalized depth
        thisTable => transectTableDepthR(:,:,tt_hydradius)
        call xsect_table_lookup_array (hydRadius, depth, thisTable, tpack) 

        !% --- get top width from normalized depth
        thisTable => transectTableDepthR(:,:,tt_width)
        call xsect_table_lookup_array (topwidth, depth, thisTable, tpack) 

        !% --- restore the physical depth
        depth(tpack) = depth(tpack) * fulldepth(tpack)

        !% --- derived data
        area0(tpack)     = area(tpack)
        area1(tpack)     = area(tpack)
        hydDepth(tpack)  = area(tpack) / topwidth(tpack)
        ell(tpack)       = hydDepth(tpack)  !% HACK -- assumes x-sect is continually increasing in width with depth
        perimeter(tpack) = area(tpack) / hydRadius(tpack)
        volume(tpack)    = area(tpack) * length(tpack)
        volume0(tpack)   = volume(tpack)
        volume1(tpack)   = volume(tpack)
        
        !%------------------------------------------------------------------
        !% Closing:
            deallocate(tpack)

    end subroutine init_IC_elem_transect_geometry
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_derived_data ()
        !%------------------------------------------------------------------
        !% Description:
        !% Initial conditions for data derived from data already read from
        !% the input file
        !%------------------------------------------------------------------
        !% Declarations
            logical, pointer :: isSmallVol(:)
            integer, pointer :: npack, thisP(:)
            real(8), pointer :: area(:), flowrate(:), velocity(:)
            integer          :: thisCol
        !%------------------------------------------------------------------
        !% Preliminaries
        !%------------------------------------------------------------------
        !% Aliases
            thisCol    = ep_CC_NOTsmalldepth
            npack      => npack_elemP(thisCol)
            thisP      => elemP(1:npack,thisCol)
            area       => elemR(:,er_Area)
            flowrate   => elemR(:,er_Flowrate)
            velocity   => elemR(:,er_Velocity)
        !%------------------------------------------------------------------
        
        elemR(:,er_Velocity) = zeroR

        if (npack < 1) return
        velocity(thisP) = flowrate(thisP) / area(thisP)

        where (velocity(thisP) > setting%Limiter%Velocity%Maximum)
           ! velocity(thisP) = zeroR
            velocity(thisP) = 0.99d0
        end where
        
        elemR(:,er_Velocity_N0) = velocity
        elemR(:,er_Velocity_N1) = velocity

        !%------------------------------------------------------------------
        !% Closing
    end subroutine init_IC_derived_data
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_solver_select (whichSolver)
        !%------------------------------------------------------------------
        !% Desscription
        !% select the solver based on depth for all the elements
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: whichSolver
            character(64)       :: subroutine_name = 'init_IC_solver_select'
        !%------------------------------------------------------------------
        !% Preliminaries:
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        select case (whichSolver)
        case (ETM)
            where ( (elemI(:,ei_HeqType) == time_march) .or. &
                    (elemI(:,ei_QeqType) == time_march) )
                elemI(:,ei_tmType) = ETM
            endwhere
        case (AC)
            print*, 'In, ', subroutine_name
            print*, 'AC solver is not handeled at this moment'
            !stop 
            call util_crashpoint(83974)
            !return
        case (ETM_AC)
            print*, 'In, ', subroutine_name
            print*, 'ETM-AC solver is not handeled at this moment'
            !stop 
            call util_crashpoint(2975)
            !return
        case default
            print *, 'In, ', subroutine_name
            print *, 'CODE ERROR: unknown solver, ', whichSolver
            print *, 'which has key ',trim(reverseKey(whichSolver))
            !stop 
            call util_crashpoint(81878)
            !return
        end select

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_solver_select
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_small_values_diagnostic_elements ()
        !--------------------------------------------------------------------------
        !
        !% set the volume, area, head, other geometry, and flow to zero values
        !% for the diagnostic elements so no error is induced in the primary
        !% face update
        !
        !--------------------------------------------------------------------------

            character(64)       :: subroutine_name = 'init_IC_small_values_diagnostic_elements'

        !--------------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        where ( (elemI(:,ei_QeqType) == diagnostic) .or. (elemI(:,ei_HeqType) == diagnostic))
            !% HACK: settings%ZeroValues should be used here
            !% when the code is finalized
            elemR(:,er_Area)     = setting%ZeroValue%Area  ! 1.0e-6
            elemR(:,er_Topwidth) = setting%ZeroValue%Topwidth !1.0e-6
            elemR(:,er_HydDepth) = setting%ZeroValue%Depth  ! 1.0e-6
            elemR(:,er_Flowrate) = zeroR
            elemR(:,er_Head)     = setting%ZeroValue%Depth + elemR(:,er_Zbottom) !1.0e-6
        endwhere

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_small_values_diagnostic_elements
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_diagnostic_interpolation_weights ()
        !--------------------------------------------------------------------------
        !
        !% set the interpolation weights for diagnostic elements
        !
        !--------------------------------------------------------------------------

            character(64)       :: subroutine_name = 'init_IC_diagnostic_interpolation_weights'

            integer, pointer ::  Npack, thisP(:), tM
            integer :: ii, kk, tB
        !--------------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"


        !% Q-diagnostic elements will have minimum interp weights for Q
        !% and maximum interp values for G and H
        !% Theses serve to force the Q value of the diagnostic element to the faces
        !% the G and H values are obtained from adjacent elements.
        where (elemI(:,ei_QeqType) == diagnostic)
            elemR(:,er_InterpWeight_uQ) = setting%Limiter%InterpWeight%Minimum
            elemR(:,er_InterpWeight_dQ) = setting%Limiter%InterpWeight%Minimum
            elemR(:,er_InterpWeight_uG) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_dG) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_uH) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_dH) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_uP) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_dP) = setting%Limiter%InterpWeight%Maximum
        endwhere

        !% H-diagnostic elements will have minimum interp weights for H and G
        !% and maximum interp weights for Q
        !% These serve to force the G, and H values of the diagnostic element to the faces
        !% and the Q value is obtained from adjacent elements
        where (elemI(:,ei_HeqType) == diagnostic)
            elemR(:,er_InterpWeight_uQ) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_dQ) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_uG) = setting%Limiter%InterpWeight%Minimum
            elemR(:,er_InterpWeight_dG) = setting%Limiter%InterpWeight%Minimum
            elemR(:,er_InterpWeight_uH) = setting%Limiter%InterpWeight%Minimum
            elemR(:,er_InterpWeight_dH) = setting%Limiter%InterpWeight%Minimum
            elemR(:,er_InterpWeight_uP) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_dP) = setting%Limiter%InterpWeight%Maximum
        endwhere

        ! !% brh 20220204 -- ccommenting this approach
        ! !% Branch elements have invariant interpolation weights so are computed here
        ! !% These are designed so that the face of a JB gets the flowrate from the
        ! !% adjacent CC conduit or channel, but the geometry and head are from the JB.
        ! Npack => npack_elemP(ep_JM_ALLtm)
        ! if (Npack > 0) then
        !     thisP  => elemP(1:Npack,ep_JM_ALLtm)
        !     do ii=1,Npack
        !         tM => thisP(ii) !% junction main ID
        !         do kk=1,max_branch_per_node
        !             tB = tM + kk !% junction branch ID
        !             elemR(tB,er_InterpWeight_uQ) = setting%Limiter%InterpWeight%Maximum
        !             elemR(tB,er_InterpWeight_dQ) = setting%Limiter%InterpWeight%Maximum
        !             elemR(tB,er_InterpWeight_uG) = setting%Limiter%InterpWeight%Minimum
        !             elemR(tB,er_InterpWeight_dG) = setting%Limiter%InterpWeight%Minimum
        !             elemR(tB,er_InterpWeight_uH) = setting%Limiter%InterpWeight%Minimum
        !             elemR(tB,er_InterpWeight_dH) = setting%Limiter%InterpWeight%Minimum
        !         end do
        !     end do
        ! end if

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_diagnostic_interpolation_weights

!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_branch_dummy_values ()    
        !%------------------------------------------------------------------
        !% Description:
        !% assigns dummy values that non-zero and not excessivley large 
        !% to the velocity, depth, area, etc. of non-valid junction branches
        !% This allows these branches to be used in array computations 
        !% without causing either divide by zero or overflow/underflow.
        !%------------------------------------------------------------------
            integer, pointer :: npack, thisP(:), BranchExists(:)
            integer :: thisCol, ii
        !%------------------------------------------------------------------
            thisCol = ep_JM
            npack   => npack_elemP(thisCol)
            if (npack < 1) return
            thisP         => elemP(1:npack,thisCol)
            BranchExists  => elemSI(:,esi_JunctionBranch_Exists)
        !%------------------------------------------------------------------

        do ii=1,max_branch_per_node
            where (BranchExists(thisP+ii) == 0)
                elemR(thisP+ii,er_Area) = 0.33333
                elemR(thisP+ii,er_Depth) = 0.33333
                elemR(thisP+ii,er_Flowrate) = 0.33333
                elemR(thisP+ii,er_Head) = 0.33333
                elemR(thisP+ii,er_Velocity) = 0.33333
                elemR(thisP+ii,er_Volume) = 0.33333 
            end where
        end do

    end subroutine init_IC_branch_dummy_values
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_set_SmallVolumes ()
        !%------------------------------------------------------------------
        !% set the small volume values in elements that are used with
        !% the small depth cutoff to change the time-marching solution
        !% The ratio of the (fixed)small volume to the actual volume is 
        !% used to blend the CM small-volume solution with the SVE solution.
        !%------------------------------------------------------------------
        !% Declarations:
            character(64)       :: subroutine_name = 'init_IC_set_SmallVolumes'
            !logical, pointer    :: useSmallVolumes
            real(8), pointer    :: depthCutoff, smallVolume(:), length(:)
            real(8), pointer    :: theta(:), radius(:), rectB(:), tempDepth(:)
            real(8), pointer    :: trapB(:), trapL(:), trapR(:), depth(:)
            integer, pointer    :: geoType(:), tPack(:), eIdx(:)
            integer             :: npack, ii, indx
        !%------------------------------------------------------------------
        !% Preliminaries
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases
            !useSmallVolumes  => setting%SmallDepth%UseSmallVolumesYN
            depthCutoff      => setting%SmallDepth%DepthCutoff
            geoType          => elemI(1:N_elem(this_image()),ei_geometryType)
            tPack            => elemI(1:N_elem(this_image()),ei_Temp01)
            eIdx             => elemI(1:N_elem(this_image()),ei_Lidx)
            smallVolume      => elemR(1:N_elem(this_image()),er_SmallVolume)
            length           => elemR(1:N_elem(this_image()),er_Length)
            depth            => elemR(1:N_elem(this_image()),er_Depth)
            tempDepth        => elemR(1:N_elem(this_image()),er_Temp01)
            rectB            => elemSGR(1:N_elem(this_image()),esgr_Rectangular_Breadth)
            radius           => elemSGR(1:N_elem(this_image()),esgr_Circular_Radius)
            trapL            => elemSGR(1:N_elem(this_image()),esgr_Trapezoidal_LeftSlope)
            trapR            => elemSGR(1:N_elem(this_image()),esgr_Trapezoidal_RightSlope)
            trapB            => elemSGR(1:N_elem(this_image()),esgr_Trapezoidal_Breadth)
            theta            => elemR(1:N_elem(this_image()),er_Temp01)
        !%------------------------------------------------------------------
        !% More preliminaries
            elemR(:,er_SmallVolume) = zeroR

            !% error checking for circular pipes
            theta = zeroR ! temporary use of theta space for comparison, this isn't theta!
            where (geoType == circular)
                theta = radius - depthCutoff
            end where
            if (any(theta < zeroR)) then
                print *, 'FATAL ERROR'
                print *, 'Small Volume depth cutoff ',depthCutoff
                print *, 'is larger or equal to radius of the smallest pipe '
                print *, 'Small Volume depth cutoff must be smaller than the radius.'
                !stop 
                call util_crashpoint(398705)
                !return
            end if
        !%------------------------------------------------------------------

        !% --- temporarily store depth, and replace depth with the cutoff
        !%     so that we can use standard depth-to-area-functions.
        !%     must be reversed at end of subroutine
        tempDepth = Depth
        depth = depthCutoff
                
        !% --- rectangular conduit
        tPack = zeroI
        npack = count(geoType == rectangular)
        if (npack > 0) then
            tPack(1:npack) = pack(eIdx,geoType == rectangular)
            smallvolume(tPack(1:npack)) = rectangular_area_from_depth(tPack(1:npack)) * length(tPack(1:npack))
        end if
        !where (geoType == rectangular)
        !    !smallVolume = depthCutoff * length * rectB
        !end where
      
        !% --- trapezoidal conduit  
        tPack = zeroI
        npack = count(geoType == trapezoidal)
        if (npack > 0) then
            tPack(1:npack) = pack(eIdx,geoType == trapezoidal)
            smallvolume(tPack(1:npack)) = trapezoidal_area_from_depth(tPack(1:npack)) * length(tPack(1:npack))
        end if  
        ! where (geoType == trapezoidal)
        !     !rm smallVolume = trapB * depthCutoff  + onehalfR*( trapL + trapR ) * depthCutoff * depthCutoff
        !     smallVolume = (trapB * depthCutoff  + onehalfR*( trapL + trapR ) * depthCutoff * depthCutoff) * length !% 20220122brh
        ! end where

        !% --- triangular conduit  
        tPack = zeroI
        npack = count(geoType == triangular)
        if (npack > 0) then
            tPack(1:npack) = pack(eIdx,geoType == triangular)
            smallvolume(tPack(1:npack)) = triangular_area_from_depth(tPack(1:npack)) * length(tPack(1:npack))
        end if 

        !% ---  circular conduit


        tPack = zeroI
        npack = count(geoType == circular)
        if (npack > 0) then
            tPack(1:npack) = pack(eIdx,geoType == circular)
            !% HACK -- temporary until problem for small volumes in circular conduits is fixed
            smallvolume(tPack(1:npack)) = 7.d0 * elemSGR(tPack(1:npack),esgr_Circular_Diameter) * depthCutoff * length
        end if
           
            ! do ii=1,npack
            !     indx = tPack(ii)
            !     if (elemI(indx,ei_elementType) .eq. CC) then
            !         print *, indx, trim(reverseKey(elemI(indx,ei_elementType)))
            !         print *, elemSGR(indx,esgr_Circular_Diameter), elemSGR(indx,esgr_Circular_Radius)
            !         print *, elemSGR(indx, esgr_Circular_YoverYfull), elemSGR(indx,esgr_Circular_AoverAfull)
            !         call circular_depth_from_volume (elemPGetm, npack_elemPGetm(epg_CC_circular_nonsurcharged), epg_CC_circular_nonsurcharged)
            !         elemR(indx,er_Volume) = 1.0d6 * elemR(indx,er_Volume)
            !         print *, 'after depth',elemR(indx,er_Depth)
            !         print *, 'after volume',elemR(indx,er_Length) * circular_area_from_depth_singular (indx) 
            !     end if
            ! end do

            ! tPack(1:npack) = pack(eIdx,geoType == circular)
            ! !% HACK -- we need an array lookup that can use the temporary pack.
            ! do ii=1,npack
            !     smallvolume(tPack(ii)) = circular_area_from_depth_singular(tPack(ii)) * length(tPack(ii))
            !     !elemR(tPack(ii),er_Volume) = smallvolume(tPack(ii)) !% temporary for debugging
            !     !print *, ii, smallvolume(tPack(ii)), circular_depth_from_volume_singular(tPack(ii))
            ! end do
            ! elemR(tPack(1:npack),er_Volume) = smallvolume(tPack(ii))
            ! print *, 'before', elemR(tPack(1:npack),er_Depth)
            ! call circular_depth_from_volume (elemPGetm, npack_elemPGetm(epg_CC_circular_nonsurcharged), epg_CC_circular_nonsurcharged)
            ! print *, 'after ', elemR(tPack(1:npack),er_Depth)

        !stop 
        !call util_crashpoint(397803)

        !% restore the initial condition depth to the depth vector
        depth = tempDepth
 
        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_set_SmallVolumes
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_set_zero_lateral_inflow ()
        !%-----------------------------------------------------------------
        !% set all the lateral inflows to zero before start of a simulation
        !%-----------------------------------------------------------------

        elemR(:,er_FlowrateLateral) = zeroR

    end subroutine init_IC_set_zero_lateral_inflow
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_oneVectors ()
        !%-----------------------------------------------------------------
        !% set up a vector of real ones (useful in sign functions)
        !%-----------------------------------------------------------------

        elemR(:,er_ones) = oneR

    end subroutine init_IC_oneVectors
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_slot ()
        !%-----------------------------------------------------------------
        !% set all the slot values to zero before start of a simulation
        !%-----------------------------------------------------------------
            character(64)       :: subroutine_name = 'init_IC_slot'
        !------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !------------------------------------------------------------------
        elemR(1:size(elemR,1)-1,er_SlotVolume)       = zeroR
        ! elemR(1:20,er_SlotVolume)                    = 0.01 

        elemR(1:size(elemR,1)-1,er_SlotArea)              = zeroR
        ! elemR(1:20,er_SlotArea)                           = elemR(1:20,er_SlotVolume) / elemR(1:20,er_length)

        elemR(1:size(elemR,1)-1,er_SlotWidth)             = zeroR
        ! elemR(1:20,er_SlotWidth)                          = (setting%Constant%gravity * elemR(1:20,er_FullArea) )/(100.0**2.0)

        elemR(1:size(elemR,1)-1,er_SlotDepth)             = zeroR
        ! elemR(1:20,er_SlotDepth)                          = elemR(1:20,er_SlotArea)/elemR(1:20,er_SlotWidth) 
        ! elemR(1:20,er_SlotDepth)                          = 1.0 / (elemR(1:20,er_area)/(10.0 * elemR(1:20,er_SlotArea))-1.0)
        
        ! elemR(1:20,er_SlotWidth)                          = elemR(1:20,er_SlotArea) / elemR(1:20,er_SlotDepth) 
        
        elemR(1:size(elemR,1)-1,er_SlotHydRadius)         = zeroR
        elemR(1:size(elemR,1)-1,er_Preissmann_Celerity)   = zeroR
        elemR(1:size(elemR,1)-1,er_Preissmann_Number)     = setting%PreissmannSlot%TargetPreissmannCelerity / &
                                                            (setting%PreissmannSlot%PreissmannAlpha * sqrt(setting%Constant%gravity * &
                                                            elemR(1:size(elemR,1)-1,er_ell_max)))
        ! elemR(1:20,er_Preissmann_Celerity)                = sqrt(9.81 * elemR(1:20,er_area)/elemR(1:20,er_SlotWidth) )
        ! faceR(1:20,fr_Area_d) = elemR(1:20,er_SlotWidth) * elemR(1:20,er_Length) * elemR(1:20,er_SlotDepth)  
        !------------------------------------------------------------------
            if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_slot
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_reference_head ()
        !%------------------------------------------------------------------
        !% Description:
        !% computes the reference head 
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: Npack, thisP(:)
            integer :: er_set(5), esr_set(8), fr_set(3)
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------   
        !% Aliases
            !% --- use only the time-marching elements to set reference head
            Npack => npack_elemP(ep_ALLtm)
            thisP => elemP(1:Npack,ep_ALLtm)
        !%------------------------------------------------------------------  

        !% --- Get the reference head  
        if (Npack > zeroI) then      
            if (setting%Solver%SubtractReferenceHead) then
                setting%Solver%ReferenceHead  = minval(elemR(thisP,er_Zbottom))
                setting%Solver%ReferenceHead  &
                    = real(floor(setting%Solver%ReferenceHead),8) 
            else
                setting%Solver%ReferenceHead = zeroR
            end if
            setting%Solver%MaxZbottom     = maxval(elemR(thisP,er_Zbottom))
            setting%Solver%MinZbottom     = minval(elemR(thisP,er_Zbottom))
            setting%Solver%AverageZbottom =    sum(elemR(thisP,er_Zbottom)) / Npack
        end if
        !% set the min and max over all the processors.
        call co_min(setting%Solver%ReferenceHead)
        call co_max(setting%Solver%MaxZbottom)
        call co_min(setting%Solver%MinZbottom)
        call co_max(setting%Solver%AverageZbottom)

        !% --- bug checking
        ! if (this_image() == 1) then
        !     write(*,*) 'reference head   ',setting%Solver%ReferenceHead
        !     write(*,*) 'max z bottom     ',setting%Solver%MaxZbottom
        !     write(*,*) 'average z bottom ',setting%Solver%AverageZbottom
        !     write(*,*) 'min z bottom     ',setting%Solver%MinZbottom
        ! end if
  
        !%------------------------------------------------------------------    
        !% Closing 
    end subroutine init_reference_head
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_subtract_reference_head ()
        !%------------------------------------------------------------------
        !% Description:
        !% removes reference head from Z values in all arrays
        !% except for BC, which is done in bc_fetch()
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: Npack, thisP(:)
            integer :: er_set(5), esr_set(8), fr_set(3)
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------   
        !% Aliases
            !% --- use only the time-marching elements to set reference head
            Npack => npack_elemP(ep_ALLtm)
            thisP => elemP(1:Npack,ep_ALLtm)
        !%------------------------------------------------------------------   
        !% list of indexes using Z reference
                ! er_Head
                ! er_Head_N0
                ! er_Zbottom
                ! er_ZbreadthMax
                ! er_Zcrown
                ! esr_Weir_NominalDownstreamHead
                ! esr_Weir_Zcrown
                ! esr_Weir_Zcrest
                ! esr_Orifice_NominalDownstreamHead
                ! esr_Orifice_Zcrown
                ! esr_Orifce_Zcrest
                ! esr_Outlet_NominalDownstreamHead
                ! esr_Outlet_Zcrest
                ! fr_Head_u
                ! fr_Head_d
                ! fr_Zbottom

        !% --- subtract the reference head from elemR
        er_set = (/er_Head,        &
                  er_Head_N0,     &
                  er_Zbottom,      &
                  er_ZbreadthMax, &
                  er_Zcrown/)
        where (elemR(:,er_set) .ne. nullValueR)        
            elemR(:,er_set) = elemR(:,er_set) - setting%Solver%ReferenceHead
        endwhere

        ! !% --- subtract the reference head from elemSR
        ! esr_set = (/esr_Weir_NominalDownstreamHead,    &
        !            esr_Weir_Zcrown,                   &
        !            esr_Weir_Zcrest,                   &
        !            esr_Orifice_NominalDownstreamHead, &
        !            esr_Orifice_Zcrown,                &
        !            esr_Orifice_Zcrest,                &
        !            esr_Outlet_NominalDownstreamHead,  &
        !            esr_Outlet_Zcrest/)

        ! where (elemSR(:,esr_set) .ne. nullValueR)        
        !        elemSR(:,esr_set) = elemSR(:,esr_set) - setting%Solver%ReferenceHead
        ! endwhere
        
        !% substract the reference head from weir elemSR
        where (elemI(:,ei_elementType) == weir)      
               elemSR(:,esr_Weir_NominalDownstreamHead) = elemSR(:,esr_Weir_NominalDownstreamHead) &
                                                - setting%Solver%ReferenceHead
               elemSR(:,esr_Weir_Zcrown) = elemSR(:,esr_Weir_Zcrown) - setting%Solver%ReferenceHead
               elemSR(:,esr_Weir_Zcrest) = elemSR(:,esr_Weir_Zcrest) - setting%Solver%ReferenceHead
        endwhere

        !% substract the reference head from orifice elemSR
        where (elemI(:,ei_elementType) == orifice)      
               elemSR(:,esr_Orifice_NominalDownstreamHead) = elemSR(:,esr_Orifice_NominalDownstreamHead) &
                                                - setting%Solver%ReferenceHead
               elemSR(:,esr_Orifice_Zcrown) = elemSR(:,esr_Orifice_Zcrown) - setting%Solver%ReferenceHead
               elemSR(:,esr_Orifice_Zcrest) = elemSR(:,esr_Orifice_Zcrest) - setting%Solver%ReferenceHead
        endwhere

        !% substract the reference head from outlet elemSR
        where (elemI(:,ei_elementType) == outlet)      
               elemSR(:,esr_Outlet_NominalDownstreamHead) = elemSR(:,esr_Outlet_NominalDownstreamHead) &
                                                - setting%Solver%ReferenceHead
               elemSR(:,esr_Outlet_Zcrest) = elemSR(:,esr_Outlet_Zcrest) - setting%Solver%ReferenceHead
        endwhere

        !% --- subtract the refence head from faceR
        fr_set = (/fr_Head_u, &
                  fr_Head_d, &
                  fr_Zbottom/)
        where (faceR(:,fr_set) .ne. nullValueR)        
               faceR(:,fr_set) = faceR(:,fr_set) - setting%Solver%ReferenceHead
        endwhere          

        !%------------------------------------------------------------------    
        !% Closing 
    end subroutine init_subtract_reference_head
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_bc()
        !%-----------------------------------------------------------------------------
        !%
        !% Description:
        !%    Initializes boundary connditions
        !%
        !% Notes:
        !%    The structures are general enough to support 3 types of BCs:
        !%
        !%    BCup: updstream boundary condition which can be inflow or head BC
        !%    BCdn: downstream boundary condition which can be inflow or head BC
        !%    BClat: lateral inflow coming into and nJ2 or nJm node.
        !%
        !%    However, the code only supports inflow BCs for BCup and BClat,
        !%    and head BCs for BCdn, mimimcking EPA-SWMM 5.13 functionalities.
        !%    Further developments allowing other types of inflow and head BCs,
        !%    should store the respective BC in either the BC%inflowX or the
        !%    BC%headX arrays defining the corresponding type of BC (i.e., BCup,
        !%    BCdn, and BClat) in the BC%xI(:,bi_category) column.
        !%
        !%-----------------------------------------------------------------------------
        integer :: ii, nidx, ntype, counter_bc_er, outfallType
        integer :: SWMMtseriesIdx, SWMMbasepatType
        character(64) :: subroutine_name = "init_bc"
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        
        if (setting%Profile%useYN) call util_profiler_start (pfc_init_bc)

        call pack_nodes()
        call util_allocate_bc()

        ! do ii=1,size(BC%flowI,DIM=1)
        !     write(*,"(10i8)"), BC%flowI(ii,bi_idx), BC%flowI(ii,bi_node_idx), BC%flowI(ii,bi_face_idx), &
        !     BC%flowI(ii,bi_elem_idx), BC%flowI(ii,bi_category), BC%flowI(ii,bi_subcategory), BC%flowI(ii,bi_fetch)
        ! end do 
        ! stop 
        !call util_crashpoint(39766)

        !% --- Convention to denote that xR_timeseries arrays haven't been fetched
        if (N_flowBC > 0) then
            BC%flowI(:,bi_fetch) = 1
            BC%flowI(:,bi_TS_upper_idx) = 0  !% latest position of upper bound in flow table
            !% --- Convention to denote association between nodes and face/elements
            !%     BCup and BCdn BCs are associated with faces, thus bi_elem_idx is null
            !%    BClat BCs are associated with elements, thus bi_face_idx is null
            BC%flowI(:, bi_face_idx) = nullvalueI
            BC%flowI(:, bi_elem_idx) = nullvalueI
            BC%flowTimeseries = nullValueR
        end if
        !print *, 'here ddd'
        if (N_headBC > 0) then
            BC%headI = nullvalueI
            BC%headI(:,bi_fetch) = 1
            BC%headI(:,bi_TS_upper_idx) = 0
            BC%headTimeseries = nullValueR
        end if

        !% --- Initialize Inflow BCs
        if (N_flowBC > 0) then
            do ii = 1, N_flowBC
                nidx  = node%P%have_flowBC(ii)
                ntype = node%I(nidx, ni_node_type)

                !% Handle Inflow BCs (BCup and BClat only)
                if (node%YN(nidx, nYN_has_extInflow) .or. node%YN(nidx, nYN_has_dwfInflow)) then
                    if ((ntype == nJm) .or. (ntype == nJ2)) then
                        BC%flowI(ii, bi_category) = BClat
                        BC%flowI(ii, bi_elem_idx) = node%I(nidx, ni_elemface_idx) !% elem idx
                    else if (ntype == nBCup) then
                        BC%flowI(ii, bi_category) = BCup
                        BC%flowI(ii, bi_face_idx) = node%I(nidx, ni_elemface_idx) !% face idx
                    else
                        print *, "Error, BC type can't be an inflow BC for node " // node%Names(nidx)%str
                        !stop 
                        call util_crashpoint(739845)
                        !return
                    end if

                    BC%flowI(ii, bi_node_idx) = nidx
                    BC%flowI(ii, bi_idx) = ii
                    BC%flowYN(ii, bYN_read_input_file) = .true.

                    !print *, 'in ',trim(subroutine_name)
                    !print *, ii, trim(node%Names(nidx)%str), node%I(nidx, ni_elemface_idx)

                    !% HACK Pattern needs checking --- the following may be wrong! brh20211221
                    !% check whether there is a pattern (-1 is no pattern) for this inflow
                    SWMMbasepatType = &
                        interface_get_nodef_attribute(nidx, api_nodef_extInflow_basePat_type)
                    !% brh20211216 should be _type?    
                    !rm nbasepat = &
                    !rm    interface_get_nodef_attribute(nidx, api_nodef_extInflow_basePat)
                    
                    !% check whether there is a time series 
                    !% (-1 is none, >0 is index, API_NULL_VALUE_I is error, which crashes API)
                    SWMMtseriesIdx = &
                        interface_get_nodef_attribute(nidx, api_nodef_extInflow_tSeries)

                    !% BC does not have fixed value if its associated with dwfInflow
                    !% or if extInflow has tseries or pattern
                    BC%flowI(ii, bi_subcategory) = BCQ_tseries
                    if (.not. node%YN(nidx, nYN_has_dwfInflow)) then !% extInflow only
                        !% brh 20211216 modified the following for basepatType == -1 rather than basepat /= -1
                        if ((SWMMtseriesIdx == -1) .and. (SWMMbasepatType == -1)) then
                            BC%flowI(ii, bi_subcategory) = BCQ_fixed
                        end if
                    end if
                else
                    print *, "There is an error, only nodes with extInflow or dwfInflow can have inflow BC"
                    !stop 
                    call util_crashpoint(826549)
                    !return
                end if
            end do
        end if

        !% --- Initialize Head BCs
        if (N_headBC > 0) then
            do ii = 1, N_headBC
                nidx = node%P%have_headBC(ii)
                ntype = node%I(nidx, ni_node_type)

                if (ntype == nBCdn) then
                    BC%headI(ii, bi_category) = BCdn
                    BC%headI(ii, bi_face_idx) = node%I(nidx, ni_elemface_idx) !% face idx
                else
                    print *, "Error, BC type can't be a head BC for node " // node%Names(nidx)%str
                    !stop 
                    call util_crashpoint(57635)
                    !return
                end if

                BC%headI(ii, bi_idx) = ii
                BC%headI(ii, bi_node_idx) = nidx

                !% --- get the outfall type
                outfallType = int(interface_get_nodef_attribute(nidx, api_nodef_outfall_type))
                select case (outfallType)
                case (API_FREE_OUTFALL)
                    BC%headI(ii, bi_subcategory) = BCH_free
                    BC%headYN(ii, bYN_read_input_file) = .false.

                case (API_NORMAL_OUTFALL)
                    BC%headI(ii, bi_subcategory) = BCH_normal
                    BC%headYN(ii, bYN_read_input_file) = .false.

                case (API_FIXED_OUTFALL) 
                    BC%headI(ii, bi_subcategory) = BCH_fixed
                    BC%headYN(ii, bYN_read_input_file) = .false.

                case (API_TIDAL_OUTFALL)
                    BC%headI(ii, bi_subcategory) = BCH_tidal
                    BC%headYN(ii, bYN_read_input_file) = .true.

                case (API_TIMESERIES_OUTFALL)
                    BC%headI(ii, bi_subcategory) = BCH_tseries
                    BC%headYN(ii, bYN_read_input_file) = .true.
                case default
                    print *, 'CODE ERROR: unexpected case default'
                    call util_crashpoint(33875)
                end select

                !% --- check for a flap gate
                if (interface_get_nodef_attribute(nidx, api_nodef_hasFlapGate) == oneR) then
                    BC%headYN(ii,bYN_hasFlapGate) = .true.
                else
                    BC%headYN(ii,bYN_hasFlapGate) = .false.
                endif

            end do
        end if
    
        call bc_step()
        if (crashI==1) return

        call pack_bc()

        if (setting%Profile%useYN) call util_profiler_stop (pfc_init_bc)

        if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_bc
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_bottom_slope ()
        !%------------------------------------------------------------------ 
        !% Description:
        !% computes the bottom slope of all channel and conduit elements
        !%------------------------------------------------------------------
        integer, pointer :: npack, thisP(:), fup(:), fdn(:)
        integer          :: thisCol
        real(8), pointer :: slope(:), length(:), fZbottom(:)
        !%------------------------------------------------------------------
        !% Aliases
            thisCol = ep_CC_ALLtm
            npack   => npack_elemP(thisCol)
            if (npack < 1) return
            thisP   => elemP(1:npack,thisCol)
            fup     => elemI(:,ei_Mface_uL)
            fdn     => elemI(:,ei_Mface_dL)
            slope   => elemR(:,er_BottomSlope)
            length  => elemR(:,er_Length)
            fZbottom => faceR(:,fr_Zbottom)
        !%------------------------------------------------------------------
        
        slope(thisP) =  (fZbottom(fup(thisP)) - fZbottom(fdn(thisP))) / length(thisP)

        !%------------------------------------------------------------------
    end subroutine init_IC_bottom_slope    
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_ZeroValues_nondepth ()
        !%------------------------------------------------------------------
        !% Description:
        !% ensures consistent initialization of zero values. Uses the 
        !% detha the primary setting, then sets the other
        !% values for consistency
        !% Assumes that Topwidth value is set.
        !%------------------------------------------------------------------
        !% Declarations
            real(8), pointer :: area, topwidth, volume, depth, length
            integer, pointer :: Npack, thisP(:)
            integer :: ii
        !%------------------------------------------------------------------
        !% Aliases
            area     => setting%ZeroValue%Area
            topwidth => setting%ZeroValue%Topwidth
            volume   => setting%ZeroValue%Volume
            depth    => setting%ZeroValue%Depth
            length   => setting%Discretization%NominalElemLength
        !%------------------------------------------------------------------
        if (.not. setting%ZeroValue%UseZeroValues) return

        !% --- use set of all time-marching
        Npack => npack_elemP(ep_ALLtm)
        if (Npack > 0) then
            thisP => elemP(1:Npack,ep_ALLtm)

            !% --- temproary store of initial depth and replace with zero depth
            elemR(:,er_Temp02) = elemR(:,er_Depth)
            elemR(:,er_Depth) = depth

            !% --- compute the topwidths for zero depth
            !%     temporary store initial condition topwidth
            elemR(:,er_Temp01) = elemR(:,er_Topwidth)
            !% --- get the topwidth at zero depth using packed geometry arrays
            call geo_topwidth_from_depth (elemPGalltm, npack_elemPGalltm, col_elemPGalltm)
            !% --- use the minimum topwidth at zero depth as the smallest topwidth
            topwidth = minval(elemR(thisP,er_Topwidth))             
            !% --- return initial condition values to topwidth 
            elemR(:,er_Topwidth) = elemR(:,er_Temp01)
            !% --- return initial condition values to depth
            elemR(:,er_Depth)    = elemR(:,er_Temp02)
          
            !% OLD the zero topwidth is 5% of the max breadth        
            !OLD topwidth = minval(elemR(thisP,er_BreadthMax)) / twentyR

            !% the zerovalue area is 50% of the product of zerovalue depth and topwidth
            area = onehalfR * topwidth * depth

            !% the zero value volume uses 5% of the volume at minimum depth
            volume = area * minval(elemR(thisP,er_Length)) / twentyR

            !print *, topwidth, area, depth, volume, minval(elemR(thisP,er_Length))
        else
            print *, 'unexpected error -- no time-marching elements found '
            !stop 
            call util_crashpoint(398733)
            !return
        end if

        if (depth < 1e-16) then
            print *, 'error, setting%ZeroValue%Depth is too small'
            !stop 
            call util_crashpoint(3987095)
            !return
        end if

        if (topwidth < 1e-16) then
            print *, 'error, setting%ZeroValue%TopWidth is too small'
            !stop 
            call util_crashpoint(3987095)
            !return
        end if

        if (area < 1e-16) then
            print *, 'error, setting%ZeroValue%Area is too small'
            !stop 
            call util_crashpoint(93764)
            !return
        end if

        if (volume < 1e-16) then
            print *, 'error, setting%ZeroValue%Volume is too small'
            !stop 
            call util_crashpoint(77395)
            !return
        end if


        !%------------------------------------------------------------------
        !% Closing
    end subroutine init_IC_ZeroValues_nondepth
!%
!%==========================================================================    
!% END MODULE
!%==========================================================================
!%
end module initial_condition
