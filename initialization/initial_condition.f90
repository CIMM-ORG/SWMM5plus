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
    use update
    use face
    use diagnostic_elements
    use geometry !BRHbugfix 20210813
    use circular_conduit
    use storage_geometry

    implicit none

    public :: init_IC_toplevel

    private

contains
    !
    !==========================================================================
    ! PUBLIC
    !==========================================================================
    !
    subroutine init_IC_toplevel ()
    !--------------------------------------------------------------------------
    !
    !% set up the initial conditions for all the elements
    !
    !--------------------------------------------------------------------------

        integer          :: ii
        integer, pointer :: solver
        character(64)    :: subroutine_name = 'init_IC_toplevel'

    !--------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        solver => setting%Solver%SolverSelect

        !if (setting%Output%Verbose) print *,'begin init_IC_from_linkdata'

        !% get data that can be extracted from links
        call init_IC_from_linkdata ()

        !if (setting%Output%Verbose) print *,'begin init_IC_from_nodedata'

        !% get data that can be extracted from nodes
        call init_IC_from_nodedata ()

        !if (setting%Output%Verbose) print *,'begin init_set_zero_lateral_inflow'

        !% zero out the lateral inflow column
        call init_IC_set_zero_lateral_inflow ()

        !if (setting%Output%Verbose) print *, 'begin init_IC_solver_select '

        !% update time marching type
        call init_IC_solver_select (solver)

        !if (setting%Output%Verbose) print *, 'begin pack_mask arrays_all'

        !% set up all the static packs and masks
        call pack_mask_arrays_all ()

        !if (setting%Output%Verbose) print *, 'begin init_IC_set_SmallVolumes'

        !% set small volume values in elements
        call init_IC_set_SmallVolumes ()

        !if (setting%Output%Verbose) print *, 'begin update_auxiliary_variables'

        !% initialize slots
        call init_IC_slot ()

        !% update all the auxiliary variables
        call update_auxiliary_variables (solver)

        !if (setting%Output%Verbose) print *,  'begin init_IC_diagnostic_interpolation_weights'

        !% update diagnostic interpolation weights
        !% (the interpolation weights of diagnostic elements
        !% stays the same throughout the simulation. Thus, they
        !% are only needed to be set at the top of the simulation)
        call init_IC_diagnostic_interpolation_weights()

        !if (setting%Output%Verbose) print *, 'begin  init_IC_small_values_diagnostic_elements'

        !% set small values to diagnostic element interpolation sets
        !% so that junk values does not mess up the first interpolation
        call init_IC_small_values_diagnostic_elements

        !if (setting%Output%Verbose) print *, 'begin face_interpolation '

        !% update faces
        call face_interpolation (fp_all)

        !if (setting%Output%Verbose) print *, 'begin diagnostic_toplevel'

        !% update the initial condition in all diagnostic elements
        call diagnostic_toplevel ()

        !if (setting%Output%Verbose) print *, 'begin init_IC_oneVectors'

        !% populate er_ones columns with ones
        call init_IC_oneVectors ()

        ! if (setting%Profile%YN) call util_profiler_stop (pfc_init_IC_setup)

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
    !
    !==========================================================================
    ! PRIVATE
    !==========================================================================
    !
    subroutine init_IC_from_linkdata ()
    !--------------------------------------------------------------------------
    !
    !% get the initial depth, flowrate, and geometry data from links
    !
    !--------------------------------------------------------------------------

        integer                                     :: ii, image, pLink
        integer, pointer                            :: thisLink
        integer, dimension(:), allocatable, target  :: packed_link_idx

        character(64) :: subroutine_name = 'init_IC_from_linkdata'
    !--------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% Setting the local image value
        image = this_image()

        !% pack all the link indexes in an image
        packed_link_idx = pack(link%I(:,li_idx), (link%I(:,li_P_image) == image))

        !% find the number of links in an image
        pLink = size(packed_link_idx)

        !% cycle through the links in an image
        do ii = 1,pLink
            !% necessary pointers
            thisLink    => packed_link_idx(ii)

            call init_IC_get_depth_from_linkdata (thisLink)

            call init_IC_get_flow_roughness_from_linkdata (thisLink)

            call init_IC_get_elemtype_from_linkdata(thisLink)

            call init_IC_get_geometry_from_linkdata (thisLink)

            !% HACK we need to add a small/zero volume adjustment here

            call init_IC_get_channel_conduit_velocity (thisLink)

        end do

        !% deallocate the temporary array
        deallocate(packed_link_idx)

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_from_linkdata
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_IC_get_depth_from_linkdata (thisLink)
    !--------------------------------------------------------------------------
    !
    !% get the initial depth data from links
    !
    !--------------------------------------------------------------------------

        integer, intent(in) :: thisLink

        integer             :: mm, ei_max
        real(8)             :: kappa
        integer, pointer    :: LdepthType
        real(8), pointer    :: DepthUp, DepthDn

        character(64) :: subroutine_name = 'init_IC_get_depth_from_linkdata'
    !--------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% type of initial depth type
        LdepthType  => link%I(thisLink,li_InitialDepthType)

        !% up and downstream depths on this link
        DepthUp => link%R(thisLink,lr_InitialUpstreamDepth)
        DepthDn => link%R(thisLink,lr_InitialDnstreamDepth)

        !% set the depths in link elements from links
        select case (LdepthType)

            case (Uniform)

                !%  if the link has a uniform depth as an initial condition
                if (link%R(thisLink,lr_InitialDepth) .ne. nullvalueR) then

                    where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                        elemR(:,er_Depth) = link%R(thisLink,lr_InitialDepth)
                    endwhere
                else
                    where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                        elemR(:,er_Depth) = onehalfR * (DepthUp + DepthDn)
                    endwhere
                end if

            case (LinearlyVarying)

                !% if the link has linearly-varying depth
                !% depth at the upstream element (link position = 1)
                where ( (elemI(:,ei_link_Pos) == 1) .and. (elemI(:,ei_link_Gidx_BIPquick) == thisLink) )
                    elemR(:,er_Depth) = DepthUp
                endwhere

                !%  using a linear distribution over the links
                ei_max = maxval(elemI(:,ei_link_Pos), 1, elemI(:,ei_link_Gidx_BIPquick) == thisLink)

                do mm=2,ei_max
                    !% find the element that is at the mm position in the link
                    where ( (elemI(:,ei_link_Pos) == mm) .and. (elemI(:,ei_link_Gidx_BIPquick) == thisLink) )
                        !% use a linear interpolation
                        elemR(:,er_Depth) = DepthUp - (DepthUp - DepthDn) * real(mm - oneI) / real(ei_max - oneI)
                    endwhere
                end do

            case (ExponentialDecay)

                !% if the link has exponentially decayed depth
                !% depth at the upstream element (link position = 1)
                where ( (elemI(:,ei_link_Pos) == 1) .and. (elemI(:,ei_link_Gidx_BIPquick) == thisLink) )
                    elemR(:,er_Depth) = DepthUp
                endwhere

                !% find the remaining elements in the link
                ei_max = maxval(elemI(:,ei_link_Pos), 1, elemI(:,ei_link_Gidx_BIPquick) == thisLink)

                do mm=2,ei_max
                    kappa = real(mm - oneI)

                    !%  depth decreases exponentially going downstream
                    if (DepthUp - DepthDn > zeroR) then
                        where ( (elemI(:,ei_link_Pos)           == mm      ) .and. &
                                (elemI(:,ei_link_Gidx_BIPquick) == thisLink) )
                            elemR(:,er_Depth) = DepthUp - (DepthUp - DepthDn) * exp(-kappa)
                        endwhere

                    !%  depth increases exponentially going downstream
                    elseif (DepthUp - DepthDn < zeroR) then
                        where ( (elemI(:,ei_link_Pos)           == mm      ) .and. &
                                (elemI(:,ei_link_Gidx_BIPquick) == thisLink) )
                            elemR(:,er_Depth) = DepthUp + (DepthDn - DepthUp) * exp(-kappa)
                        endwhere

                    !%  uniform depth
                    else
                        where ( (elemI(:,ei_link_Pos)           == mm      ) .and. &
                                (elemI(:,ei_link_Gidx_BIPquick) == thisLink) )
                            elemR(:,er_Depth) = DepthUp
                        endwhere
                    end if
                end do

            case default
                print*, 'In ', subroutine_name
                print*, 'error: unexpected initial depth type, ', LdepthType,'  in link, ', thisLink
                stop

        end select

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_depth_from_linkdata
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
        if (icrash) return
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

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
    !
    !% get the geometry data from links
    !
    !--------------------------------------------------------------------------

        integer, intent(in) :: thisLink
        integer, pointer    :: linkType

        character(64) :: subroutine_name = 'init_IC_get_elemtype_from_linkdata'
    !--------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

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

                print*, 'In ', subroutine_name
                print*, 'pumps are not handeled yet'
                stop

            case default

                print*, 'In ', subroutine_name
                print*, 'error: unexpected link, ', linkType,'  in the network'
                stop

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
    !
    !% get the geometry data from links
    !
    !--------------------------------------------------------------------------

        integer, intent(in) :: thisLink
        integer, pointer    :: linkType

        character(64) :: subroutine_name = 'init_IC_get_flow_roughness_from_linkdata'
    !--------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% necessary pointers
        linkType      => link%I(thisLink,li_link_type)

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

                print*, 'In ', subroutine_name
                print*, 'pumps are not handeled yet'
                stop

            case default

                print*, 'In ', subroutine_name
                print*, 'error: unexpected link, ', linkType,'  in the network'
                stop

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
    !% get the geometry data for channel links
    !% and calculate element volumes
    !
    !--------------------------------------------------------------------------

        integer, intent(in) :: thisLink
        integer, pointer    :: geometryType

        character(64) :: subroutine_name = 'init_IC_get_channel_geometry'
    !--------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% pointer to geometry type
        geometryType => link%I(thisLink,li_geometry)

        select case (geometryType)

            case (lRectangular)

                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)

                    elemI(:,ei_geometryType) = rectangular

                    !% store geometry specific data
                    elemSGR(:,esgr_Rectangular_Breadth) = link%R(thisLink,lr_BreadthScale)

                    elemR(:,er_BreadthMax)   = elemSGR(:,esgr_Rectangular_Breadth)
                    elemR(:,er_Area)         = elemSGR(:,esgr_Rectangular_Breadth) * elemR(:,er_Depth)
                    elemR(:,er_Area_N0)      = elemR(:,er_Area)
                    elemR(:,er_Area_N1)      = elemR(:,er_Area)
                    elemR(:,er_Volume)       = elemR(:,er_Area) * elemR(:,er_Length)
                    elemR(:,er_Volume_N0)    = elemR(:,er_Volume)
                    elemR(:,er_Volume_N1)    = elemR(:,er_Volume)
                    !% the full depth of channel is set to a large depth so it
                    !% never surcharges. the large depth is set as, factor x width,
                    !% where the factor is an user controlled paratmeter.
                    elemR(:,er_FullDepth)    = setting%Limiter%Channel%LargeDepthFactor * &
                                                link%R(thisLink,lr_BreadthScale)
                    elemR(:,er_ZbreadthMax)  = elemR(:,er_FullDepth) + elemR(:,er_Zbottom)
                    elemR(:,er_Zcrown)       = elemR(:,er_Zbottom) + elemR(:,er_FullDepth)
                    elemR(:,er_FullArea)     = elemSGR(:,esgr_Rectangular_Breadth) * elemR(:,er_FullDepth)
                    elemR(:,er_FullVolume)   = elemR(:,er_FullArea) * elemR(:,er_Length)
                endwhere

            case (lTrapezoidal)

                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)

                    elemI(:,ei_geometryType) = trapezoidal

                    !% store geometry specific data
                    elemSGR(:,esgr_Trapezoidal_Breadth)    = link%R(thisLink,lr_BreadthScale)
                    elemSGR(:,esgr_Trapezoidal_LeftSlope)  = link%R(thisLink,lr_LeftSlope)
                    elemSGR(:,esgr_Trapezoidal_RightSlope) = link%R(thisLink,lr_RightSlope)

                    ! (Bottom width + averageSlope * Depth)*Depth
                    elemR(:,er_Area)         = (elemSGR(:,esgr_Trapezoidal_Breadth) + onehalfR * &
                                (elemSGR(:,esgr_Trapezoidal_LeftSlope) + elemSGR(:,esgr_Trapezoidal_RightSlope)) * &
                                elemR(:,er_Depth)) * elemR(:,er_Depth)

                    elemR(:,er_Area_N0)      = elemR(:,er_Area)
                    elemR(:,er_Area_N1)      = elemR(:,er_Area)
                    elemR(:,er_Volume)       = elemR(:,er_Area) * elemR(:,er_Length)
                    elemR(:,er_Volume_N0)    = elemR(:,er_Volume)
                    elemR(:,er_Volume_N1)    = elemR(:,er_Volume)

                    ! Bottom width + (lslope + rslope) * BankFullDepth
                    elemR(:,er_BreadthMax)   = elemSGR(:,esgr_Trapezoidal_Breadth) + (elemSGR(:,esgr_Trapezoidal_LeftSlope) + &
                                elemSGR(:,esgr_Trapezoidal_RightSlope)) * elemR(:,er_FullDepth)
                    !% the full depth of channel is set to a large depth so it
                    !% never surcharges. the large depth is set as, factor x width,
                    !% where the factor is an user controlled paratmeter.
                    elemR(:,er_FullDepth)    = setting%Limiter%Channel%LargeDepthFactor * &
                                                link%R(thisLink,lr_BreadthScale)
                    elemR(:,er_ZbreadthMax)  = elemR(:,er_FullDepth) + elemR(:,er_Zbottom)
                    elemR(:,er_Zcrown)       = elemR(:,er_Zbottom) + elemR(:,er_FullDepth)
                    elemR(:,er_FullArea)     = (elemSGR(:,esgr_Trapezoidal_Breadth) + onehalfR * &
                                (elemSGR(:,esgr_Trapezoidal_LeftSlope) + elemSGR(:,esgr_Trapezoidal_RightSlope)) * &
                                elemR(:,er_FullDepth)) * elemR(:,er_FullDepth)
                    elemR(:,er_FullVolume)   = elemR(:,er_FullArea) * elemR(:,er_Length)
                endwhere

            case default

                print*, 'In, ', subroutine_name
                print*, 'Only rectangular channel geometry is handeled at this moment'
                stop

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
        if (icrash) return
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% pointer to geometry type
        geometryType => link%I(thisLink,li_geometry)

        select case (geometryType)

        ! case (lRectangular)

        !     where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)

        !         elemI(:,ei_geometryType)    = rectangular_closed

        !         !% store geometry specific data
        !         elemSGR(:,esgr_Rectangular_Breadth) = link%R(thisLink,lr_BreadthScale)

        !         elemR(:,er_BreadthMax)      = elemSGR(:,esgr_Rectangular_Breadth)
        !         elemR(:,er_Area)            = elemSGR(:,esgr_Rectangular_Breadth) * elemR(:,er_Depth)
        !         elemR(:,er_Area_N0)         = elemR(:,er_Area)
        !         elemR(:,er_Area_N1)         = elemR(:,er_Area)
        !         elemR(:,er_Volume)          = elemR(:,er_Area) * elemR(:,er_Length)
        !         elemR(:,er_Volume_N0)       = elemR(:,er_Volume)
        !         elemR(:,er_Volume_N1)       = elemR(:,er_Volume)
        !         elemR(:,er_FullDepth)       = link%R(thisLink,lr_FullDepth)
        !         elemR(:,er_ZbreadthMax)     = elemR(:,er_FullDepth) + elemR(:,er_Zbottom)
        !         elemR(:,er_Zcrown)          = elemR(:,er_Zbottom) + elemR(:,er_FullDepth)
        !         elemR(:,er_FullArea)        = elemSGR(:,esgr_Rectangular_Breadth) * elemR(:,er_FullDepth)
        !         elemR(:,er_FullVolume)      = elemR(:,er_FullArea) * elemR(:,er_Length)
        !     endwhere

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

                ! elemR(ii,er_Area)                  = circular_area_from_depth_singular(ii)

                !% HACK: for initial condition steup, use analytical functions which works with where
                !% statement

                elemR(:,er_Area)                  = (elemSGR(:,esgr_Circular_Radius) **2) * &
                                (acos(1.0 - (elemR(:,er_Depth)/elemSGR(:,esgr_Circular_Radius))) - &
                                sin(2.0*acos(1.0 - (elemR(:,er_Depth)/elemSGR(:,esgr_Circular_Radius))))/2.0 )

                elemR(:,er_Area_N0)               = elemR(:,er_Area)
                elemR(:,er_Area_N1)               = elemR(:,er_Area)
                elemR(:,er_Volume)                = elemR(:,er_Area) * elemR(:,er_Length)
                elemR(:,er_Volume_N0)             = elemR(:,er_Volume)
                elemR(:,er_Volume_N1)             = elemR(:,er_Volume)
            end where

        case default

            print*, 'In, ', subroutine_name
            print*, 'Only rectangular, and circular conduit geometry is handeled at this moment'
            stop

        end select

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

        character(64) :: subroutine_name = 'init_IC_get_weir_geometry'
    !--------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% pointer to specific weir type
        specificWeirType => link%I(thisLink,li_weir_type)

        select case (specificWeirType)
            !% copy weir specific data
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
                endwhere

            case (lRoadWayWeir)

                print*, 'In ', subroutine_name
                print*, 'roadway weir is not handeled yet'
                stop

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
                endwhere

            case default

                print*, 'In ', subroutine_name
                print*, 'error: unknown weir type, ', specificWeirType,'  in network'
                stop

        end select

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_IC_get_weir_geometry
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_IC_get_orifice_geometry (thisLink)
    !--------------------------------------------------------------------------
    !
    !% get the geometry and other data data for orifice links
    !
    !--------------------------------------------------------------------------

        integer, intent(in) :: thisLink
        integer, pointer    :: specificOrificeType, OrificeGeometryType

        character(64) :: subroutine_name = 'init_IC_get_orifice_geometry'
    !--------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% pointer to specific orifice type
        specificOrificeType => link%I(thisLink,li_orif_type)

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

                print*, 'In ', subroutine_name
                print*, 'error: unknown orifice type, ', specificOrificeType,'  in network'
                stop

        end select

        !% pointer to specific orifice geometry
        OrificeGeometryType => link%I(thisLink,li_geometry)

        select case (OrificeGeometryType)
            !% copy orifice specific geometry data
            case (lRectangular)
                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    !% integer data
                    elemSI(:,ei_geometryType)          = rectangular

                    !% real data
                    elemR(:,er_FullDepth)                    = link%R(thisLink,lr_FullDepth)
                    elemSR(:,esr_Orifice_EffectiveFullDepth) = link%R(thisLink,lr_FullDepth)
                    elemSR(:,esr_Orifice_DischargeCoeff)     = link%R(thisLink,lr_DischargeCoeff1)
                    elemSR(:,esr_Orifice_Zcrest)             = elemR(:,er_Zbottom) + link%R(thisLink,lr_InletOffset)
                    elemSR(:,esr_Orifice_Zcrown)             = elemSR(:,eSr_Orifice_Zcrest) + link%R(thisLink,lr_FullDepth)
                    elemSR(:,esr_Orifice_RectangularBreadth) = link%R(thisLink,lr_BreadthScale)
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
                end where

            case default
                print*, 'In ', subroutine_name
                print*, 'error: unknown orifice geometry type, ', OrificeGeometryType,'  in network'
                stop
            end select


        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_orifice_geometry
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_IC_get_channel_conduit_velocity (thisLink)
    !--------------------------------------------------------------------------
    !
    !% get the velocity of channel and conduits
    !% and sell all other velocity to zero
    !
    !--------------------------------------------------------------------------

        integer, intent(in) :: thisLink
        integer, pointer    :: specificWeirType

        character(64) :: subroutine_name = 'init_IC_get_channel_conduit_velocity'
    !--------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% HACK: this might not be right
        where ( (elemI(:,ei_link_Gidx_BIPquick) == thisLink) .and. &
                (elemR(:,er_area)               .gt. zeroR ) .and. &
                (elemI(:,ei_elementType)        == CC      ) )

            elemR(:,er_Velocity)    = elemR(:,er_Flowrate) / elemR(:,er_Area)
            elemR(:,er_Velocity_N0) = elemR(:,er_Velocity)
            elemR(:,er_Velocity_N1) = elemR(:,er_Velocity)

        elsewhere ( (elemI(:,ei_link_Gidx_BIPquick) == thisLink) .and. &
                    (elemR(:,er_area)               .le. zeroR ) .and. &
                    (elemI(:,ei_elementType)        == CC    ) )

            elemR(:,er_Velocity)    = zeroR
            elemR(:,er_Velocity_N0) = zeroR
            elemR(:,er_Velocity_N1) = zeroR

        endwhere

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_IC_get_channel_conduit_velocity
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_IC_from_nodedata ()
    !--------------------------------------------------------------------------
    !
    !% get the initial depth, and geometry data from nJm nodes
    !
    !--------------------------------------------------------------------------

        integer                       :: ii, image, pJunction
        integer, pointer              :: thisJunctionNode
        integer, allocatable, target  :: packed_nJm_idx(:)

        character(64) :: subroutine_name = 'init_IC_from_nodedata'
    !--------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

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
    !
    !--------------------------------------------------------------------------
    !
    !% get data for the multi branch junction elements
    !
    !--------------------------------------------------------------------------
    !
        integer, intent(in) :: thisJunctionNode

        integer              :: ii, jj, JMidx, JBidx
        integer, pointer     :: BranchIdx, geometryType, JmType, curveID
        real(8), allocatable :: integrated_volume(:)

        character(64) :: subroutine_name = 'init_IC_get_junction_data'
    !--------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !%................................................................
        !% Junction main
        !%................................................................
        !% find the first element ID associated with that nJm
        !% masked on the global node number for this node.
        JMidx = minval(elemI(:,ei_Lidx), elemI(:,ei_node_Gidx_SWMM) == thisJunctionNode)

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
            !% HACK: Junction main with artifical storage are rectangular
            !%-----------------------------------------------------------------------
            elemSI(JMidx,esi_JunctionMain_Type) = ArtificalStorage
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
            elemR(JMidx,er_FullDepth) = node%R(thisJunctionNode,nr_SurchargeDepth)
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

            !% find the element id of junction branches
            JBidx = JMidx + ii

            !% set the junction branch element type
            elemI(JBidx,ei_elementType) = JB

            !% set the geometry for existing branches
            !% Note that elemSI(,...Exists) is set in init_network_handle_nJm
            if (elemSI(JBidx,esi_JunctionBranch_Exists) == oneI) then

                BranchIdx    => elemSI(JBidx,esi_JunctionBranch_Link_Connection)
                geometryType => link%I(BranchIdx,li_geometry)

                !% set the JB to time_march for use with splitting between AC
                !% and ETM in rk2_extrapolate_to_fullstep_ETM, rk2_restore_to_midstep_ETM
                !% rk2_interpolate_to_halfstep_AC, k2_restore_to_fullstep_AC
                elemI(JBidx,ei_HeqType) = time_march
                elemI(JBidx,ei_QeqType) = time_march

                !% set the initial head to the same as the junction main
                elemR(JBidx,er_Head)    = elemR(JMidx,er_Head)
                elemR(JBidx,er_Depth)   = elemR(JBidx,er_Head) - elemR(JBidx,er_Zbottom)
                if (elemR(JBidx,er_Depth) < setting%ZeroValue%Depth) then
                    elemR(JBidx,er_Depth) = setting%ZeroValue%depth
                    elemR(JBidx,er_Head)  = setting%ZeroValue%depth + elemR(JBidx,er_Zbottom)
                end if

                !% JB elements initialized for momentum
                elemR(JBidx,er_WaveSpeed)    = sqrt(setting%constant%gravity * elemR(JBidx,er_Depth))
                elemR(JBidx,er_FroudeNumber) = zeroR

                ! set the common geometry for conduits and non-conduits that are independent of cross-section shape
                if (link%I(BranchIdx,li_link_type) == lPipe) then
                    elemYN(JBidx,eYN_canSurcharge)  = .true.
                    elemR(JBidx,er_FullDepth)   = link%R(BranchIdx,lr_FullDepth)
                else
                    elemYN(JBidx,eYN_canSurcharge)  = .false.
                    elemR(JBidx,er_FullDepth)   = setting%Limiter%Channel%LargeDepthFactor * &
                                                    link%R(BranchIdx,lr_BreadthScale)
                end if
                elemR(JBidx,er_Zcrown)      = elemR(JBidx,er_Zbottom) + elemR(JBidx,er_FullDepth)

                !% Junction branch k-factor
                !% If the user does not input the K-factor for junction branches entrance/exit loses then
                !% use default from setting
                if (node%R(thisJunctionNode,nr_JunctionBranch_Kfactor) .ne. nullvalueR) then
                    elemSR(JBidx,esr_JunctionBranch_Kfactor) = node%R(thisJunctionNode,nr_JunctionBranch_Kfactor)
                else
                    elemSR(JBidx,esr_JunctionBranch_Kfactor) = setting%Junction%kFactor
                end if

                !% get the geometry data
                select case (geometryType)

                    case (lRectangular)
                        elemI(JBidx,ei_geometryType) = rectangular
                        elemR(JBidx,er_BreadthMax)   = link%R(BranchIdx,lr_BreadthScale)
                        elemR(JBidx,er_Area)         = elemR(JBidx,er_BreadthMax) * elemR(JBidx,er_Depth)

                        !% store geometry specific data
                        elemSGR(JBidx,esgr_Rectangular_Breadth) = link%R(BranchIdx,lr_BreadthScale)
                        elemR(JBidx,er_FullArea)    = elemR(JBidx,er_BreadthMax) * link%R(BranchIdx,lr_FullDepth)
                        elemR(JBidx,er_ZbreadthMax) = link%R(BranchIdx,lr_FullDepth) + elemR(JBidx,er_Zbottom)

                    case (lTrapezoidal)

                        print *, 'Need to rewrite trapezoidal inital conditions'
                        stop 3983
                        elemI(JBidx,ei_geometryType) = trapezoidal

                        !% store geometry specific data
                        elemSGR(JBidx,esgr_Trapezoidal_Breadth)    = link%R(BranchIdx,lr_BreadthScale)
                        elemSGR(JBidx,esgr_Trapezoidal_LeftSlope)  = link%R(BranchIdx,lr_LeftSlope)
                        elemSGR(JBidx,esgr_Trapezoidal_RightSlope) = link%R(BranchIdx,lr_RightSlope)


                        elemR(JBidx,er_ZBreadthMax) = link%R(BranchIdx,lr_FullDepth) + elemR(JBidx,er_Zbottom)

                        elemR(JBidx,er_BreadthMax)  = elemSGR(JBidx,esgr_Trapezoidal_Breadth) + &
                                (elemSGR(JBidx,esgr_Trapezoidal_LeftSlope) + &
                                elemSGR(JBidx,esgr_Trapezoidal_RightSlope)) * elemR(JBidx,er_ZbreadthMax)

                        elemR(JBidx,er_FullArea)    = (elemSGR(JBidx,esgr_Trapezoidal_Breadth) + onehalfR * &
                                (elemSGR(JBidx,esgr_Trapezoidal_LeftSlope) + elemSGR(JBidx,esgr_Trapezoidal_RightSlope)) * &
                                elemR(JBidx,er_FullDepth)) * elemR(JBidx,er_FullDepth) !  BRHbugfix 20210813

                    case (lCircular)
                        elemI(JBidx,ei_geometryType) = circular

                        !% store geometry specific data
                        elemSGR(JBidx,esgr_Circular_Diameter) = link%R(BranchIdx,lr_FullDepth)
                        elemSGR(JBidx,esgr_Circular_Radius)   = elemSGR(JBidx,esgr_Circular_Diameter) / twoR

                        elemR(JBidx,er_ZBreadthMax) = link%R(BranchIdx,lr_FullDepth) / twoR + elemR(JBidx,er_Zbottom)

                        elemR(JBidx,er_BreadthMax)  = link%R(BranchIdx,lr_FullDepth)

                        elemR(JBidx,er_FullArea)    = pi * elemSGR(JBidx,esgr_Circular_Radius) ** twoR

                    case default

                        print*, 'In, ', subroutine_name
                        print*, 'Only rectangular geometry is handeled at this moment'
                        stop

                end select

                !% get the flow data from links for junction branches
                !% this flowrate will always be lagged in junction branches
                elemR(JBidx,er_Flowrate) = link%R(BranchIdx,lr_InitialFlowrate)

                if (elemR(JBidx,er_Area) .gt. setting%ZeroValue%Area) then ! BRHbugfix 20210813
                    elemR(JBidx,er_Velocity) = elemR(JBidx,er_Flowrate) / elemR(JBidx,er_Area)
                else
                    elemR(JBidx,er_Velocity) = zeroR
                end if

                !% Common geometry that do not depend on cross-section
                elemR(JBidx,er_Area_N0)      = elemR(JBidx,er_Area)
                elemR(JBidx,er_Area_N1)      = elemR(JBidx,er_Area)
                elemR(JBidx,er_Volume)       = elemR(JBidx,er_Area) * elemR(JBidx,er_Length)
                elemR(JBidx,er_Volume_N0)    = elemR(JBidx,er_Volume)
                elemR(JBidx,er_Volume_N1)    = elemR(JBidx,er_Volume)

            end if
        end do

        !% HACK:
        !% set initial conditions for junction main from the junction branch data
        !% For the momentum we are simply using rectangular geometry as a damping pot for the junctions.
        !% Goal is to ensure consistency with the links and mass conservation.
        !% Need to replace how JM geometry is handled in timeloop before we change this.
        !% length -- here uses the largest 2 input and output to get a maximum length

        !% get junction main geometry based on type
        JmType => elemSI(JMidx,esi_JunctionMain_Type)

        select case (JmType)

            case (ArtificalStorage)

                elemR(JMidx,er_Length) = max(elemR(JMidx+1,er_Length), elemR(JMidx+3,er_Length), &
                                             elemR(JMidx+5,er_Length)) + &
                                         max(elemR(JMidx+2,er_Length), elemR(JMidx+4,er_Length), &
                                             elemR(JMidx+6,er_Length))

                !% HACK: finding the average breadth. This will not work for channels with other than rectangular geometry.
                !% we need to generalize this
                elemSGR(JMidx,esgr_Rectangular_Breadth) = (elemR(JMidx+1,er_Length)*elemSGR(JMidx+1,esgr_Rectangular_Breadth) + &
                                                           elemR(JMidx+2,er_Length)*elemSGR(JMidx+2,esgr_Rectangular_Breadth) + &
                                                           elemR(JMidx+3,er_Length)*elemSGR(JMidx+3,esgr_Rectangular_Breadth) + &
                                                           elemR(JMidx+4,er_Length)*elemSGR(JMidx+4,esgr_Rectangular_Breadth) + &
                                                           elemR(JMidx+5,er_Length)*elemSGR(JMidx+5,esgr_Rectangular_Breadth) + &
                                                           elemR(JMidx+6,er_Length)*elemSGR(JMidx+6,esgr_Rectangular_Breadth))/ &
                                                           elemR(JMidx,er_Length)

                !% Volume
                !% rectangular volume depends on characteristic length and breadth.
                elemR(JMidx,er_Volume) =   elemSGR(JMidx,esgr_Rectangular_Breadth) * elemR(JMidx,er_Length) * elemR(JMidx,er_Depth)

                elemR(JMidx,er_Volume_N0) = elemR(JMidx,er_Volume)
                elemR(JMidx,er_Volume_N1) = elemR(JMidx,er_Volume)

            case (FunctionalStorage)
                elemR(JMidx,er_Volume) = elemSR(JMidx,esr_Storage_Constant) * elemR(JMidx,er_Depth) +          &
                    (elemSR(JMidx,esr_Storage_Coefficient) / (elemSR(JMidx,esr_Storage_Exponent) + oneR)) * &
                    elemR(JMidx,er_Depth) ** (elemSR(JMidx,esr_Storage_Exponent) + oneR)

                elemR(JMidx,er_Volume_N0) = elemR(JMidx,er_Volume)
                elemR(JMidx,er_Volume_N1) = elemR(JMidx,er_Volume)

                !% create a storage curve
                call storage_create_curve (JMidx)

            case (TabularStorage)
                curveID => elemSI(JMidx,esi_JunctionMain_Curve_ID)
                Curve(curveID)%ElemIdx = JMidx
                !% SWMM5+ needs a volume vs depth relationship thus Trapezoidal rule is used
                !% to get to integrate the area vs depth curve
                call storage_integrate_area_vs_depth_curve (curveID)

                !% now interpolate from the cure to get the volume
                call storage_interpolate_volume_from_depth_singular (JMidx)

            case default
                print*, 'In, ', subroutine_name
                print*, 'error: unknown junction main type, ', JmType
                stop 54895

        end select

        ! call the standard geometry update for junction branches
        call geo_assign_JB (ALLtm, ep_JM_ALLtm)

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_junction_data
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_IC_solver_select (solver)
    !--------------------------------------------------------------------------
    !
    !% select the solver based on depth for all the elements
    !
    !--------------------------------------------------------------------------

        integer, intent(in) :: solver
        character(64)       :: subroutine_name = 'init_IC_solver_select'

    !--------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"


        select case (solver)

            case (ETM)

                where ( (elemI(:,ei_HeqType) == time_march) .or. &
                        (elemI(:,ei_QeqType) == time_march) )

                    elemI(:,ei_tmType) = ETM

                endwhere

            case (AC)

                print*, 'In, ', subroutine_name
                print*, 'AC solver is not handeled at this moment'
                stop 83974

            case (ETM_AC)

                print*, 'In, ', subroutine_name
                print*, 'ETM-AC solver is not handeled at this moment'
                stop 2975
            case default

                print*, 'In, ', subroutine_name
                print*, 'error: unknown solver, ', solver
                stop 81878

        end select



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
        if (icrash) return
        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        where ( (elemI(:,ei_QeqType) == diagnostic) .or. (elemI(:,ei_HeqType) == diagnostic))
            !% HACK: settings%ZeroValues should be used here
            !% when the code is finalized
            elemR(:,er_Area)     = 1.0e-6
            elemR(:,er_Topwidth) = 1.0e-6
            elemR(:,er_HydDepth) = 1.0e-6
            elemR(:,er_Flowrate) = 1.0e-6
            elemR(:,er_Head)     = 1.0e-6
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
        if (icrash) return
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
        endwhere

        !% Branch elements have invariant interpolation weights so are computed here
        !% These are designed so that the face of a JB gets the flowrate from the
        !% adjacent CC conduit or channel, but the geometry and head are from the JB.
        Npack => npack_elemP(ep_JM_ALLtm)
        if (Npack > 0) then
            thisP  => elemP(1:Npack,ep_JM_ALLtm)
            do ii=1,Npack
                tM => thisP(ii) !% junction main ID
                do kk=1,max_branch_per_node
                    tB = tM + kk !% junction branch ID
                    elemR(tB,er_InterpWeight_uQ) = setting%Limiter%InterpWeight%Maximum
                    elemR(tB,er_InterpWeight_dQ) = setting%Limiter%InterpWeight%Maximum
                    elemR(tB,er_InterpWeight_uG) = setting%Limiter%InterpWeight%Minimum
                    elemR(tB,er_InterpWeight_dG) = setting%Limiter%InterpWeight%Minimum
                    elemR(tB,er_InterpWeight_uH) = setting%Limiter%InterpWeight%Minimum
                    elemR(tB,er_InterpWeight_dH) = setting%Limiter%InterpWeight%Minimum
                end do
            end do
        end if

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_diagnostic_interpolation_weights
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_IC_set_SmallVolumes ()
    !--------------------------------------------------------------------------
    !
    !% set the small volume values in elements
    !
    !--------------------------------------------------------------------------

        character(64)       :: subroutine_name = 'init_IC_set_SmallVolumes'

    !--------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        if (setting%SmallVolume%UseSmallVolumes) then
            where (elemI(:,ei_geometryType) == rectangular)
                elemR(:,er_SmallVolume) = setting%SmallVolume%DepthCutoff * elemSGR(:,esgr_Rectangular_Breadth) * &
                    elemR(:,er_Length)

            elsewhere (elemI(:,ei_geometryType) == trapezoidal)
                elemR(:,er_SmallVolume) = (elemSGR(:,esgr_Trapezoidal_Breadth) + onehalfR * &
                    (elemSGR(:,esgr_Trapezoidal_LeftSlope) + elemSGR(:,esgr_Trapezoidal_RightSlope)) * &
                    setting%SmallVolume%DepthCutoff) * setting%SmallVolume%DepthCutoff

            elsewhere (elemI(:,ei_geometryType) == circular)
                ! HACK: not the correct small volume according to geometry. But it will work for now
                elemR(:,er_SmallVolume) = pi * (setting%SmallVolume%DepthCutoff / twoR) ** twoR
            endwhere
        else
            elemR(:,er_SmallVolume) = zeroR
        end if

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_set_SmallVolumes
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_IC_set_zero_lateral_inflow ()
    !--------------------------------------------------------------------------
    !
    !% set all the lateral inflows to zero before start of a simulation
    !
    !--------------------------------------------------------------------------

        character(64)       :: subroutine_name = 'init_IC_set_zero_lateral_inflow'

    !--------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        elemR(1:size(elemR,1)-1,er_FlowrateLateral) = zeroR

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_set_zero_lateral_inflow
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_IC_oneVectors ()
    !--------------------------------------------------------------------------
    !
    !% set all the lateral inflows to zero before start of a simulation
    !
    !--------------------------------------------------------------------------

        character(64)       :: subroutine_name = 'init_IC_oneVectors'

    !--------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        elemR(1:size(elemR,1)-1,er_ones) = oneR

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_oneVectors
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_IC_slot ()
    !--------------------------------------------------------------------------
    !
    !% set all the slot values to zero before start of a simulation
    !
    !--------------------------------------------------------------------------

        character(64)       :: subroutine_name = 'init_IC_slot'

    !--------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        elemR(1:size(elemR,1)-1,er_SlotWidth)             = zeroR
        elemR(1:size(elemR,1)-1,er_SlotVolume)            = zeroR
        elemR(1:size(elemR,1)-1,er_SlotDepth)             = zeroR
        elemR(1:size(elemR,1)-1,er_SlotArea)              = zeroR
        elemR(1:size(elemR,1)-1,er_SlotHydRadius)         = zeroR
        elemR(1:size(elemR,1)-1,er_Preissmann_Celerity)   = zeroR

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_slot
    !
    !==========================================================================
    !==========================================================================
    !
end module initial_condition
