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

    implicit none

    public :: init_IC_setup

    private

contains
    !
    !==========================================================================
    ! PUBLIC
    !==========================================================================
    !
    subroutine init_IC_setup ()
    !--------------------------------------------------------------------------
    !
    !% set up the initial conditions for all the elements
    !
    !--------------------------------------------------------------------------

        integer          :: ii
        integer, pointer :: solver
        character(64)    :: subroutine_name = 'init_IC_setup'

    !--------------------------------------------------------------------------
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

        solver => setting%Solver%SolverSelect

        !if (setting%Verbose) print *,'begin init_IC_from_linkdata'

        !% get data that can be extracted from links
        call init_IC_from_linkdata ()

        !if (setting%Verbose) print *,'begin init_IC_from_nodedata'

        !% get data that can be extracted from nodes
        call init_IC_from_nodedata ()

        !if (setting%Verbose) print *,'begin init_set_zero_lateral_inflow'

        !% zero out the lateral inflow column
        call init_IC_set_zero_lateral_inflow ()

        !if (setting%Verbose) print *, 'begin init_IC_solver_select '

        !% update time marching type
        call init_IC_solver_select (solver)

        !if (setting%Verbose) print *, 'begin pack_mask arrays_all'

        !% set up all the static packs and masks
        call pack_mask_arrays_all ()

        !if (setting%Verbose) print *, 'begin init_IC_set_SmallVolumes'

        !% set small volume values in elements 
        call init_IC_set_SmallVolumes ()

        !if (setting%Verbose) print *, 'begin update_auxiliary_variables'

        !% update all the auxiliary variables
        call update_auxiliary_variables (solver)
     
        !if (setting%Verbose) print *,  'begin init_IC_diagnostic_interpolation_weights'

        !% update diagnostic interpolation weights
        !% (the interpolation weights of diagnostic elements
        !% stays the same throughout the simulation. Thus, they
        !% are only needed to be set at the top of the simulation)
        call init_IC_diagnostic_interpolation_weights()
 
        !if (setting%Verbose) print *, 'begin  init_IC_small_values_diagnostic_elements'

        !% set small values to diagnostic element interpolation sets
        !% so that junk values does not mess up the first interpolation
        call init_IC_small_values_diagnostic_elements

        !if (setting%Verbose) print *, 'begin face_interpolation '

        !% update faces
        call face_interpolation (fp_all)

        !if (setting%Verbose) print *, 'begin diagnostic_toplevel'

        !% update the initial condition in all diagnostic elements
        call diagnostic_toplevel ()

        !if (setting%Verbose) print *, 'begin init_IC_oneVectors'

        !% populate er_ones columns with ones
        call init_IC_oneVectors ()

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
           call execute_command_line('')
        end if

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine init_IC_setup
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
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

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
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

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

                    where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
                        elemR(:,er_Depth) = link%R(thisLink,lr_InitialDepth)
                    endwhere
                else
                    where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
                        elemR(:,er_Depth) = onehalfR * (DepthUp + DepthDn)
                    endwhere
                end if

            case (LinearlyVarying)

                !% if the link has linearly-varying depth
                !% depth at the upstream element (link position = 1)
                where ( (elemI(:,ei_link_Pos) == 1) .and. (elemI(:,ei_link_Gidx_SWMM) == thisLink) )
                    elemR(:,er_Depth) = DepthUp
                endwhere

                !%  using a linear distribution over the links
                ei_max = maxval(elemI(:,ei_link_Pos), 1, elemI(:,ei_link_Gidx_SWMM) == thisLink)

                do mm=2,ei_max
                    !% find the element that is at the mm position in the link
                    where ( (elemI(:,ei_link_Pos) == mm) .and. (elemI(:,ei_link_Gidx_SWMM) == thisLink) )
                        !% use a linear interpolation
                        elemR(:,er_Depth) = DepthUp - (DepthUp - DepthDn) * real(mm - oneI) / real(ei_max - oneI)
                    endwhere
                end do

            case (ExponentialDecay)

                !% if the link has exponentially decayed depth
                !% depth at the upstream element (link position = 1)
                where ( (elemI(:,ei_link_Pos) == 1) .and. (elemI(:,ei_link_Gidx_SWMM) == thisLink) )
                    elemR(:,er_Depth) = DepthUp
                endwhere

                !% find the remaining elements in the link
                ei_max = maxval(elemI(:,ei_link_Pos), 1, elemI(:,ei_link_Gidx_SWMM) == thisLink)

                do mm=2,ei_max
                    kappa = real(mm - oneI)

                    !%  depth decreases exponentially going downstream
                    if (DepthUp - DepthDn > zeroR) then
                        where ( (elemI(:,ei_link_Pos)       == mm      ) .and. &
                                (elemI(:,ei_link_Gidx_SWMM) == thisLink) )
                            elemR(:,er_Depth) = DepthUp - (DepthUp - DepthDn) * exp(-kappa)
                        endwhere

                    !%  depth increases exponentially going downstream
                    elseif (DepthUp - DepthDn < zeroR) then
                        where ( (elemI(:,ei_link_Pos)       == mm      ) .and. &
                                (elemI(:,ei_link_Gidx_SWMM) == thisLink) )
                            elemR(:,er_Depth) = DepthUp + (DepthDn - DepthUp) * exp(-kappa)
                        endwhere

                    !%  uniform depth
                    else
                        where ( (elemI(:,ei_link_Pos)       == mm      ) .and. &
                                (elemI(:,ei_link_Gidx_SWMM) == thisLink) )
                            elemR(:,er_Depth) = DepthUp
                        endwhere
                    end if
                end do

            case default
                print*, 'In ', subroutine_name
                print*, 'error: unexpected initial depth type, ', LdepthType,'  in link, ', thisLink
                stop "in " // subroutine_name

        end select

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

        !%  handle all the initial conditions that don't depend on geometry type
        where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
            elemR(:,er_Flowrate)       = link%R(thisLink,lr_InitialFlowrate)
            elemR(:,er_Flowrate_N0)    = link%R(thisLink,lr_InitialFlowrate)
            elemR(:,er_Flowrate_N1)    = link%R(thisLink,lr_InitialFlowrate)
            elemR(:,er_Roughness)      = link%R(thisLink,lr_Roughness)
        endwhere

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

        !% necessary pointers
        linkType      => link%I(thisLink,li_link_type)

        select case (linkType)

            case (lChannel)

                where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
                    elemI(:,ei_elementType)     = CC
                    elemI(:,ei_HeqType)         = time_march
                    elemI(:,ei_QeqType)         = time_march
                endwhere


            case (lpipe)

                where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
                    elemI(:,ei_elementType)     = CC
                    elemI(:,ei_HeqType)         = time_march
                    elemI(:,ei_QeqType)         = time_march
                    elemYN(:,eYN_canSurcharge)  =  .true.
                endwhere

            case (lweir)

                where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
                    elemI(:,ei_elementType)     = weir
                    elemI(:,ei_QeqType)         = diagnostic
                    elemYN(:,eYN_canSurcharge)  = link%YN(thisLink,lYN_CanSurcharge)
                endwhere

            case (lOrifice)

                where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
                    elemI(:,ei_elementType)     = orifice
                    elemI(:,ei_QeqType)         = diagnostic
                    elemYN(:,eYN_canSurcharge)  = .true.
                endwhere

            case (lPump)

                print*, 'In ', subroutine_name
                print*, 'pumps are not handeled yet'
                stop "in " // subroutine_name

            case default

                print*, 'In ', subroutine_name
                print*, 'error: unexpected link, ', linkType,'  in the network'
                stop "in " // subroutine_name

        end select


        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

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
                stop "in " // subroutine_name

            case default

                print*, 'In ', subroutine_name
                print*, 'error: unexpected link, ', linkType,'  in the network'
                stop "in " // subroutine_name

        end select


        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

        !% pointer to geometry type
        geometryType => link%I(thisLink,li_geometry)

        select case (geometryType)

            case (lRectangular)

                where (elemI(:,ei_link_Gidx_SWMM) == thisLink)

                    elemI(:,ei_geometryType) = rectangular

                    !% store geometry specific data
                    elemSGR(:,eSGR_Rectangular_Breadth) = link%R(thisLink,lr_BreadthScale)

                    elemR(:,er_BreadthMax)   = elemSGR(:,eSGR_Rectangular_Breadth)
                    elemR(:,er_Area)         = elemSGR(:,eSGR_Rectangular_Breadth) * elemR(:,er_Depth)
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
                    elemR(:,er_FullArea)     = elemSGR(:,eSGR_Rectangular_Breadth) * elemR(:,er_FullDepth)
                    elemR(:,er_FullVolume)   = elemR(:,er_FullArea) * elemR(:,er_Length)
                endwhere

            case (lTrapezoidal)

                where (elemI(:,ei_link_Gidx_SWMM) == thisLink)

                    elemI(:,ei_geometryType) = trapezoidal

                    !% store geometry specific data
                    elemSGR(:,eSGR_Trapezoidal_Breadth)    = link%R(thisLink,lr_BreadthScale)
                    elemSGR(:,eSGR_Trapezoidal_LeftSlope)  = link%R(thisLink,lr_LeftSlope)
                    elemSGR(:,eSGR_Trapezoidal_RightSlope) = link%R(thisLink,lr_RightSlope)

                    ! (Bottom width + averageSlope * Depth)*Depth
                    elemR(:,er_Area)         = (elemSGR(:,eSGR_Trapezoidal_Breadth) + onehalfR * &
                                (elemSGR(:,eSGR_Trapezoidal_LeftSlope) + elemSGR(:,eSGR_Trapezoidal_RightSlope)) * &
                                elemR(:,er_Depth)) * elemR(:,er_Depth)

                    elemR(:,er_Area_N0)      = elemR(:,er_Area)
                    elemR(:,er_Area_N1)      = elemR(:,er_Area)
                    elemR(:,er_Volume)       = elemR(:,er_Area) * elemR(:,er_Length)
                    elemR(:,er_Volume_N0)    = elemR(:,er_Volume)
                    elemR(:,er_Volume_N1)    = elemR(:,er_Volume)

                    ! Bottom width + (lslope + rslope) * BankFullDepth
                    elemR(:,er_BreadthMax)   = elemSGR(:,eSGR_Trapezoidal_Breadth) + (elemSGR(:,eSGR_Trapezoidal_LeftSlope) + &
                                elemSGR(:,eSGR_Trapezoidal_RightSlope)) * elemR(:,er_FullDepth)
                    !% the full depth of channel is set to a large depth so it
                    !% never surcharges. the large depth is set as, factor x width,
                    !% where the factor is an user controlled paratmeter.
                    elemR(:,er_FullDepth)    = setting%Limiter%Channel%LargeDepthFactor * &
                                                link%R(thisLink,lr_BreadthScale)
                    elemR(:,er_ZbreadthMax)  = elemR(:,er_FullDepth) + elemR(:,er_Zbottom)
                    elemR(:,er_Zcrown)       = elemR(:,er_Zbottom) + elemR(:,er_FullDepth)
                    elemR(:,er_FullArea)     = (elemSGR(:,eSGR_Trapezoidal_Breadth) + onehalfR * &
                                (elemSGR(:,eSGR_Trapezoidal_LeftSlope) + elemSGR(:,eSGR_Trapezoidal_RightSlope)) * &
                                elemR(:,er_FullDepth)) * elemR(:,er_FullDepth)
                    elemR(:,er_FullVolume)   = elemR(:,er_FullArea) * elemR(:,er_Length)
                endwhere

            case default

                print*, 'In, ', subroutine_name
                print*, 'Only rectangular channel geometry is handeled at this moment'
                stop "in " // subroutine_name

        end select

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

        !% pointer to geometry type
        geometryType => link%I(thisLink,li_geometry)

        select case (geometryType)

        ! case (lRectangular)

        !     where (elemI(:,ei_link_Gidx_SWMM) == thisLink)

        !         elemI(:,ei_geometryType)    = rectangular_closed

        !         !% store geometry specific data
        !         elemSGR(:,eSGR_Rectangular_Breadth) = link%R(thisLink,lr_BreadthScale)

        !         elemR(:,er_BreadthMax)      = elemSGR(:,eSGR_Rectangular_Breadth)
        !         elemR(:,er_Area)            = elemSGR(:,eSGR_Rectangular_Breadth) * elemR(:,er_Depth)
        !         elemR(:,er_Area_N0)         = elemR(:,er_Area)
        !         elemR(:,er_Area_N1)         = elemR(:,er_Area)
        !         elemR(:,er_Volume)          = elemR(:,er_Area) * elemR(:,er_Length)
        !         elemR(:,er_Volume_N0)       = elemR(:,er_Volume)
        !         elemR(:,er_Volume_N1)       = elemR(:,er_Volume)
        !         elemR(:,er_FullDepth)       = link%R(thisLink,lr_FullDepth)
        !         elemR(:,er_ZbreadthMax)     = elemR(:,er_FullDepth) + elemR(:,er_Zbottom)
        !         elemR(:,er_Zcrown)          = elemR(:,er_Zbottom) + elemR(:,er_FullDepth)
        !         elemR(:,er_FullArea)        = elemSGR(:,eSGR_Rectangular_Breadth) * elemR(:,er_FullDepth)
        !         elemR(:,er_FullVolume)      = elemR(:,er_FullArea) * elemR(:,er_Length)
        !     endwhere

        case (lCircular)

            do ii = 1,N_elem(this_image())
                if (elemI(ii,ei_link_Gidx_SWMM) == thisLink) then

                    elemI(ii,ei_geometryType)    = circular

                    !% store geometry specific data
                    elemSGR(ii,eSGR_Circular_Diameter) = link%R(thisLink,lr_BreadthScale)
                    elemSGR(ii,eSGR_Circular_Radius)   = link%R(thisLink,lr_BreadthScale) / twoR

                    elemR(ii,er_FullDepth)             = link%R(thisLink,lr_FullDepth)
                    elemR(ii,er_Zcrown)                = elemR(ii,er_Zbottom) + elemR(ii,er_FullDepth)
                    elemR(ii,er_ZbreadthMax)           = elemR(ii,er_FullDepth)/twoR + elemR(ii,er_Zbottom)
                    elemR(ii,er_FullArea)              = pi * elemSGR(ii,eSGR_Circular_Radius) ** twoR
                    elemR(ii,er_FullVolume)            = elemR(ii,er_FullArea) * elemR(ii,er_Length)
                    elemR(ii,er_FullHydDepth)          = elemR(ii,er_FullDepth) 
                    elemR(ii,er_FullPerimeter)         = elemR(ii,er_FullArea) / (onefourthR * elemR(ii,er_FullDepth))
                    elemR(ii,er_BreadthMax)            = elemSGR(ii,eSGR_Circular_Diameter)
                    elemR(ii,er_AreaBelowBreadthMax)   = elemR(ii,er_FullArea)  / twoR

                    elemR(ii,er_Area)                  = circular_area_from_depth_singular(ii)
                    elemR(ii,er_Area_N0)               = elemR(ii,er_Area)
                    elemR(ii,er_Area_N1)               = elemR(ii,er_Area)
                    elemR(ii,er_Volume)                = elemR(ii,er_Area) * elemR(ii,er_Length)
                    elemR(ii,er_Volume_N0)             = elemR(ii,er_Volume)
                    elemR(ii,er_Volume_N1)             = elemR(ii,er_Volume)

                    elemSR(ii,eSr_conduit_SlotWidth)   = zeroR
                    elemSR(ii,eSr_conduit_SlotVolume)  = zeroR
                    elemSR(ii,eSr_conduit_SlotDepth)   = zeroR
                    elemSR(ii,eSr_conduit_SlotArea)    = zeroR
                end if
            end do
            
        case default

            print*, 'In, ', subroutine_name
            print*, 'Only rectangular, and circular conduit geometry is handeled at this moment'
            stop "in " // subroutine_name

        end select

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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

        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        !% pointer to specific weir type
        specificWeirType => link%I(thisLink,li_weir_type)

        select case (specificWeirType)
            !% copy weir specific data
            case (lTrapezoidalWeir)

                where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
                    !% integer data
                    elemSI(:,eSi_Weir_SpecificType)          = trapezoidal_weir

                    !% real data
                    elemSR(:,eSr_Weir_EffectiveFullDepth)    = link%R(thisLink,lr_FullDepth)
                    elemSR(:,eSr_Weir_DischargeCoeff1)       = link%R(thisLink,lr_DischargeCoeff1)
                    elemSR(:,eSr_Weir_DischargeCoeff2)       = link%R(thisLink,lr_DischargeCoeff2)
                    elemSR(:,eSr_Weir_TrapezoidalBreadth)    = link%R(thisLink,lr_BreadthScale)
                    elemSR(:,eSr_Weir_TrapezoidalLeftSlope)  = link%R(thisLink,lr_LeftSlope)
                    elemSR(:,eSr_Weir_TrapezoidalRightSlope) = link%R(thisLink,lr_RightSlope)
                    elemSR(:,eSr_Weir_Zcrest)                = elemR(:,er_Zbottom) + link%R(thisLink,lr_InletOffset)

                    !% HACK: I am not sure if we need to update the initial area or volume of an weir element
                    !% since they will all be set to zero values at the start of the simulation
                endwhere

            case (lSideFlowWeir)

                where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
                    !% integer data
                    elemSI(:,eSi_Weir_SpecificType)          = side_flow
                    elemSI(:,eSi_Weir_EndContractions)       = link%I(thisLink,li_weir_EndContrations)

                    !% real data
                    elemSR(:,eSr_Weir_EffectiveFullDepth)    = link%R(thisLink,lr_FullDepth)
                    elemSR(:,eSr_Weir_DischargeCoeff2)       = link%R(thisLink,lr_DischargeCoeff2)
                    elemSR(:,eSr_Weir_RectangularBreadth)    = link%R(thisLink,lr_BreadthScale)
                    elemSR(:,eSr_Weir_Zcrest)                = elemR(:,er_Zbottom) + link%R(thisLink,lr_InletOffset)

                    !% HACK: I am not sure if we need to update the initial area or volume of an weir element
                    !% since they will all be set to zero values at the start of the simulation
                endwhere

            case (lRoadWayWeir)

                print*, 'In ', subroutine_name
                print*, 'roadway weir is not handeled yet'
                stop "in " // subroutine_name

            case (lVnotchWeir)

                where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
                    !% integer data
                    elemSI(:,eSi_Weir_SpecificType)          = vnotch_weir

                    !% real data
                    elemSR(:,eSr_Weir_EffectiveFullDepth)    = link%R(thisLink,lr_FullDepth)
                    elemSR(:,eSr_Weir_DischargeCoeff1)       = link%R(thisLink,lr_DischargeCoeff1)
                    elemSR(:,eSr_Weir_TriangularSideSlope)   = link%R(thisLink,lr_SideSlope)
                    elemSR(:,eSr_Weir_Zcrest)                = elemR(:,er_Zbottom) + link%R(thisLink,lr_InletOffset)

                    !% HACK: I am not sure if we need to update the initial area or volume of an weir element
                    !% since they will all be set to zero values at the start of the simulation
                endwhere

            case (lTransverseWeir)

                where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
                    !% integer data
                    elemSI(:,eSi_Weir_SpecificType)          = transverse_weir
                    elemSI(:,eSi_Weir_EndContractions)       = link%I(thisLink,li_weir_EndContrations)

                    !% real data
                    elemSR(:,eSr_Weir_EffectiveFullDepth)    = link%R(thisLink,lr_FullDepth)
                    elemSR(:,eSr_Weir_DischargeCoeff2)       = link%R(thisLink,lr_DischargeCoeff2)
                    elemSR(:,eSr_Weir_RectangularBreadth)    = link%R(thisLink,lr_BreadthScale)
                    elemSR(:,eSr_Weir_Zcrest)                = elemR(:,er_Zbottom)  + link%R(thisLink,lr_InletOffset)

                    !% HACK: I am not sure if we need to update the initial area or volume of an weir element
                    !% since they will all be set to zero values at the start of the simulation
                endwhere

            case default

                print*, 'In ', subroutine_name
                print*, 'error: unknown weir type, ', specificWeirType,'  in network'
                stop "in " // subroutine_name

        end select

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

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

        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        !% pointer to specific orifice type
        specificOrificeType => link%I(thisLink,li_orif_type)

        select case (specificOrificeType)
            !% copy orifice specific data
            case (lBottomOrifice)

                where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
                    !% integer data
                    elemSI(:,eSi_specific_orifice_type)      = bottom_orifice

                endwhere

            case (lSideOrifice)

                where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
                    !% integer data
                    elemSI(:,eSi_specific_orifice_type)       = side_orifice

                endwhere

            case default

                print*, 'In ', subroutine_name
                print*, 'error: unknown orifice type, ', specificOrificeType,'  in network'
                stop "in " // subroutine_name

        end select

        !% pointer to specific orifice geometry
        OrificeGeometryType => link%I(thisLink,li_geometry)

        select case (OrificeGeometryType)
            !% copy orifice specific geometry data
            case (lRectangular)
                where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
                    !% integer data
                    elemSI(:,ei_geometryType)          = rectangular

                    !% real data
                    elemSR(:,eSr_Weir_EffectiveFullDepth)    = link%R(thisLink,lr_FullDepth)
                    elemSR(:,eSr_Orifice_DischargeCoeff)     = link%R(thisLink,lr_DischargeCoeff1)
                    elemSR(:,eSr_Orifice_Zcrest)             = elemR(:,er_Zbottom) + link%R(thisLink,lr_InletOffset)
                    elemSR(:,eSr_Orifice_RectangularBreadth) = link%R(thisLink,lr_BreadthScale)

                endwhere

            case (lCircular)
                print*, 'In ', subroutine_name
                print *, 'error, the circular orifice is not yet implemented'
                stop "in " // subroutine_name

            case default
                print*, 'In ', subroutine_name
                print*, 'error: unknown orifice geometry type, ', OrificeGeometryType,'  in network'
                stop "in " // subroutine_name
            end select


        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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

        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        !% HACK: this might not be right
        where ( (elemI(:,ei_link_Gidx_SWMM) == thisLink) .and. &
                (elemR(:,er_area)           .gt. zeroR   ) .and. &
                (elemI(:,ei_elementType)    == CC      ) )

            elemR(:,er_Velocity)    = elemR(:,er_Flowrate) / elemR(:,er_Area)
            elemR(:,er_Velocity_N0) = elemR(:,er_Velocity)
            elemR(:,er_Velocity_N1) = elemR(:,er_Velocity)

        elsewhere ( (elemI(:,ei_link_Gidx_SWMM) == thisLink) .and. &
                    (elemR(:,er_area)           .le. zeroR   ) .and. &
                    (elemI(:,ei_elementType)    == CC      ) )

            elemR(:,er_Velocity)    = zeroR
            elemR(:,er_Velocity_N0) = zeroR
            elemR(:,er_Velocity_N1) = zeroR

        endwhere

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

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
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

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
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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

        integer             :: ii, JMidx, JBidx
        integer, pointer    :: BranchIdx, geometryType

        character(64) :: subroutine_name = 'init_IC_get_junction_data'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        !%................................................................
        !% Junction main
        !%................................................................
        !% find the first element ID associated with that nJm
        !% masked on the global node number for this node.
        JMidx = minval(elemI(:,ei_Lidx), elemI(:,ei_node_Gidx_SWMM) == thisJunctionNode)

        !% the first element index is a junction main
        elemI(JMidx,ei_elementType)  = JM
        elemI(JMidx,ei_HeqType)      = time_march
        elemI(JMidx,ei_QeqType)      = none
        
        !%-----------------------------------------------------------------------
        !% HACK: Junction main are always rectangular
        elemI(JMidx,ei_geometryType) = rectangular
        !%-----------------------------------------------------------------------

        elemR(JMidx,er_Depth)        = node%R(thisJunctionNode,nr_InitialDepth)
        !% JM elements are not solved for momentum. 
        elemR(JMidx,er_Flowrate)     = zeroR
        elemR(JMidx,er_Velocity)     = zeroR
        !% wave speed is the gravity wave speed for the depth
        elemR(JMidx,er_WaveSpeed)    = sqrt(setting%constant%gravity * elemR(JMidx,er_Depth))
        elemR(JMidx,er_FroudeNumber) = zeroR

        !% find if the node can surcharge
        if (node%R(thisJunctionNode,nr_SurchargeDepth) .ne. nullValueR) then
            elemYN(JMidx,eYN_canSurcharge)  = .true.
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
            if (elemSI(JBidx,eSI_JunctionBranch_Exists) == oneI) then

                BranchIdx    => elemSI(JBidx,eSI_JunctionBranch_Link_Connection)
                geometryType => link%I(BranchIdx,li_geometry)

                !% set the head equation as as JB for now for existing branches
                !% only JM is time marched 
                !% HACK -- might need to get a realistic depth -- this perhaps should be done
                !% as an adjustment after all the other IC is done (e.g., we may want to compare
                !% with the upstream element depth.)
                elemI(JBidx,ei_HeqType) = JB
                elemR(JBidx,er_Depth)   = elemR(JMidx,er_Depth)
                elemR(JBidx,er_Head)    = elemR(JBidx,er_Depth) + elemR(JBidx,er_Zbottom) !% BRHbugfix 20210813
                

                !% HACK: JB elements are not solved for momentum. 
                !% set zero here and revise later if needed
                ! BRHbugfix 20210813 elemR(JBidx,er_WaveSpeed)    = sqrt(setting%constant%gravity * elemR(JBidx,er_Depth))
                ! BRHbugfix 20210813 elemR(JBidx,er_FroudeNumber) = zeroR
                ! BRHbugfix end

                !% BRHbugfix 20210813 start: moved these from inside select
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
                !% If the user does not input the k-factor for junction branches usee default from setting
                if (node%R(thisJunctionNode,nr_JunctionBranch_Kfactor) .ne. nullvalueR) then
                    elemSR(JBidx,eSr_JunctionBranch_Kfactor) = node%R(thisJunctionNode,nr_JunctionBranch_Kfactor)
                else
                    elemSR(JBidx,eSr_JunctionBranch_Kfactor) = setting%Junction%kFactor
                end if                
                !% BRHbugfix 20210813  end 

                !% get the geometry data
                select case (geometryType)

                    case (lRectangular)
                        elemI(JBidx,ei_geometryType) = rectangular
                        elemR(JBidx,er_BreadthMax)   = link%R(BranchIdx,lr_BreadthScale)
                        elemR(JBidx,er_Area)         = elemR(JBidx,er_BreadthMax) * elemR(JBidx,er_Depth)
                        ! BRHbugfix 20210813 elemR(JBidx,er_Area_N0)      = elemR(JBidx,er_Area)
                        ! BRHbugfix 20210813 elemR(JBidx,er_Area_N1)      = elemR(JBidx,er_Area)
                        ! BRHbugfix 20210813 elemR(JBidx,er_Volume)       = elemR(JBidx,er_Area) * elemR(JBidx,er_Length)
                        ! BRHbugfix 20210813 elemR(JBidx,er_Volume_N0)    = elemR(JBidx,er_Volume)
                        ! BRHbugfix 20210813 elemR(JBidx,er_Volume_N1)    = elemR(JBidx,er_Volume)

                        !% store geometry specific data
                        elemSGR(JBidx,eSGR_Rectangular_Breadth) = link%R(BranchIdx,lr_BreadthScale)
                        elemR(JBidx,er_FullArea)    = elemR(JBidx,er_BreadthMax) * elemR(JBidx,er_FullDepth)
                        elemR(JBidx,er_ZbreadthMax) = zeroR

                        !% HACK: not sure if we need surcharge condition for junction branches
                        !% ANSWER -- yes, we do need surcharge on JB for conduits, but should be more general
                        !% That is, does not need to be in rectangular sectoin
                        ! BRHbugfix 20210813 if (link%I(BranchIdx,li_link_type) == lPipe) then
                            ! BRHbugfix 20210813 elemYN(JBidx,eYN_canSurcharge)  = .true.
                            ! BRHbugfix 20210813 elemR(JBidx,er_FullDepth)   = link%R(BranchIdx,lr_FullDepth)
                            !elemR(JBidx,er_Zcrown)      = elemR(JBidx,er_Zbottom) + elemR(JBidx,er_FullDepth)
                            !elemR(JBidx,er_FullArea)    = elemR(JBidx,er_BreadthMax) * elemR(JBidx,er_FullDepth)
                            ! BRHbugfix 20210813 elemR(JBidx,er_FullVolume)  = elemR(JBidx,er_FullArea) * elemR(JBidx,er_Length)
                            ! BRHbugfix 20210813 elemR(JBidx,er_BreadthMax)  = zeroR
                            ! BRHbugfix 20210813 elemR(JBidx,er_ZbreadthMax) = zeroR
                        ! BRHbugfix 20210813 else
                            !% the full depth of channel is set to a large depth so it
                            !% never surcharges. the large depth is set as, factor x width,
                            !% where the factor is an user controlled paratmeter.
                            ! BRHbugfix 20210813 elemR(JBidx,er_FullDepth)   = setting%Limiter%Channel%LargeDepthFactor * &
                            ! BRHbugfix 20210813                             link%R(BranchIdx,lr_BreadthScale)
                            !elemR(JBidx,er_Zcrown)      = elemR(JBidx,er_Zbottom) + elemR(JBidx,er_FullDepth)
                            !elemR(JBidx,er_FullArea)    = elemR(JBidx,er_BreadthMax) * elemR(JBidx,er_FullDepth)
                            ! BRHbugfix 20210813 elemR(JBidx,er_FullVolume)  = elemR(JBidx,er_FullArea) * elemR(JBidx,er_Length)
                            ! BRHbugfix 20210813 elemR(JBidx,er_ZbreadthMax)  = link%R(BranchIdx,lr_FullDepth)
                        ! BRHbugfix 20210813 end if


                    case (lTrapezoidal)

                            print *, 'Need to rewrite trapezoidal inital conditions'
                            stop 3983
                            elemI(JBidx,ei_geometryType) = trapezoidal

                            !% store geometry specific data
                            elemSGR(JBidx,eSGR_Trapezoidal_Breadth)    = link%R(BranchIdx,lr_BreadthScale)
                            elemSGR(JBidx,eSGR_Trapezoidal_LeftSlope)  = link%R(BranchIdx,lr_LeftSlope)
                            elemSGR(JBidx,eSGR_Trapezoidal_RightSlope) = link%R(BranchIdx,lr_RightSlope)

                            ! (Bottom width + averageSlope * Depth)*Depth
                            ! BRHbugfix 20210813 elemR(JBidx,er_Area)         = (elemSGR(JBidx,eSGR_Trapezoidal_Breadth) + onehalfR * &
                            ! BRHbugfix 20210813             (elemSGR(JBidx,eSGR_Trapezoidal_LeftSlope) + elemSGR(JBidx,eSGR_Trapezoidal_RightSlope)) * &
                            ! BRHbugfix 20210813             elemR(JBidx,er_Depth)) * elemR(JBidx,er_Depth)

                            ! BRHbugfix 20210813 elemR(JBidx,er_Area_N0)      = elemR(JBidx,er_Area)
                            ! BRHbugfix 20210813 elemR(JBidx,er_Area_N1)      = elemR(JBidx,er_Area)
                            ! BRHbugfix 20210813 elemR(JBidx,er_Volume)       = elemR(JBidx,er_Area) * elemR(JBidx,er_Length)
                            ! BRHbugfix 20210813 elemR(JBidx,er_Volume_N0)    = elemR(JBidx,er_Volume)
                            ! BRHbugfix 20210813 elemR(JBidx,er_Volume_N1)    = elemR(JBidx,er_Volume)

                            !% Junction branch k-factor
                            !% HACK: if the user does not input the k-factor for junction brnaches,
                            !% get a default value from the setting
                            ! BRHbugfix 20210813 if (node%R(thisJunctionNode,nr_JunctionBranch_Kfactor) .ne. nullvalueR) then
                            ! BRHbugfix 20210813     elemSR(JBidx,eSr_JunctionBranch_Kfactor) = node%R(thisJunctionNode,nr_JunctionBranch_Kfactor)
                            ! BRHbugfix 20210813 else
                            ! BRHbugfix 20210813     elemSR(JBidx,eSr_JunctionBranch_Kfactor) = setting%Junction%kFactor
                            ! BRHbugfix 20210813 end if
                            
                            elemR(JBidx,er_ZBreadthMax)  =elemR(JBidx,er_FullDepth)  ! BRHbugfix 20210813 !

                            if (link%I(BranchIdx,li_link_type) == lPipe) then

                                ! BRHbugfix 20210813 elemYN(JBidx,eYN_canSurcharge)  = .true.

                                ! BRHbugfix 20210813 elemR(JBidx,er_FullDepth)   = link%R(BranchIdx,lr_FullDepth)
                                ! BRHbugfix 20210813 elemR(JBidx,er_Zcrown)      = elemR(JBidx,er_Zbottom) + elemR(JBidx,er_FullDepth)
                                ! BRHbugfix 20210813 elemR(JBidx,er_FullArea)    = elemR(JBidx,er_BreadthMax) * elemR(JBidx,er_FullDepth)
                                !  BRHbugfix 20210813 elemR(JBidx,er_FullArea)    = (elemSGR(JBidx,eSGR_Trapezoidal_Breadth) + onehalfR * &
                                !  BRHbugfix 20210813         (elemSGR(JBidx,eSGR_Trapezoidal_LeftSlope) + elemSGR(JBidx,eSGR_Trapezoidal_RightSlope)) * &
                                 !  BRHbugfix 20210813        elemR(JBidx,er_FullDepth)) * elemR(JBidx,er_FullDepth)
                                ! ! BRHbugfix 20210813 elemR(JBidx,er_BreadthMax)  = zeroR

                            else
                                elemR(JBidx,er_BreadthMax)  = elemSGR(JBidx,eSGR_Trapezoidal_Breadth) + &
                                            (elemSGR(JBidx,eSGR_Trapezoidal_LeftSlope) + &
                                            elemSGR(JBidx,eSGR_Trapezoidal_RightSlope)) * elemR(JBidx,er_ZbreadthMax)
                                !  BRHbugfix 20210813 elemR(JBidx,er_FullDepth)   = setting%Limiter%Channel%LargeDepthFactor * &
                                !    BRHbugfix 20210813                     link%R(BranchIdx,lr_BreadthScale)
                                !  BRHbugfix 20210813 elemR(JBidx,er_Zcrown)      = elemR(JBidx,er_Zbottom) + elemR(JBidx,er_FullDepth)
                           
                                !  BRHbugfix 20210813 elemR(JBidx,er_FullArea)    = (elemSGR(JBidx,eSGR_Trapezoidal_Breadth) + onehalfR * &
                                !  BRHbugfix 20210813         (elemSGR(JBidx,eSGR_Trapezoidal_LeftSlope) + elemSGR(JBidx,eSGR_Trapezoidal_RightSlope)) * &
                                 !  BRHbugfix 20210813        elemR(JBidx,er_FullDepth)) * elemR(JBidx,er_FullDepth)
                                !  BRHbugfix 20210813 elemR(JBidx,er_FullVolume)  = elemR(JBidx,er_FullArea) * elemR(JBidx,er_Length)
                            end if
                            elemR(JBidx,er_FullArea)    = (elemSGR(JBidx,eSGR_Trapezoidal_Breadth) + onehalfR * &
                                    (elemSGR(JBidx,eSGR_Trapezoidal_LeftSlope) + elemSGR(JBidx,eSGR_Trapezoidal_RightSlope)) * &
                                    elemR(JBidx,er_FullDepth)) * elemR(JBidx,er_FullDepth) !  BRHbugfix 20210813 

                    case default

                        print*, 'In, ', subroutine_name
                        print*, 'Only rectangular geometry is handeled at this moment'
                        stop "in " // subroutine_name

                end select

                !% get the flow data from links for junction branches
                !% this flowrate will always be lagged in junction branches
                elemR(JBidx,er_Flowrate) = link%R(BranchIdx,lr_InitialFlowrate)

                ! BRHbugfix 20210813  if (elemR(JBidx,er_Area) .gt. zeroR) then
                if (elemR(JBidx,er_Area) .gt. setting%ZeroValue%Area) then ! BRHbugfix 20210813
                    elemR(JBidx,er_Velocity) = elemR(JBidx,er_Flowrate) / elemR(JBidx,er_Area)
                else
                    elemR(JBidx,er_Velocity) = zeroR
                end if

                ! BRHbugfix 20210813 start
                !% Common geometry that do not depend on cross-section
                elemR(JBidx,er_Area_N0)      = elemR(JBidx,er_Area)
                elemR(JBidx,er_Area_N1)      = elemR(JBidx,er_Area)
                elemR(JBidx,er_Volume)       = elemR(JBidx,er_Area) * elemR(JBidx,er_Length)
                elemR(JBidx,er_Volume_N0)    = elemR(JBidx,er_Volume)
                elemR(JBidx,er_Volume_N1)    = elemR(JBidx,er_Volume)
                ! BRHbugfix 20210813 end

            end if
        end do

        !% HACK: 
        !% set initial conditions for junction main from the junction branch data
        !% For the momentum we are simply using rectangular geometry as a damping pot for the junctions.
        !% Goal is to ensure consistency with the links and mass conservation.
        !% Need to replace how JM geometry is handled in timeloop before we change this.
        !% length -- here uses the largest 2 input and output to get a maximum length
        elemR(JMidx,er_Length) = max(elemR(JMidx+1,er_Length), elemR(JMidx+3,er_Length), &
                                     elemR(JMidx+5,er_Length)) + &
                                 max(elemR(JMidx+2,er_Length), elemR(JMidx+4,er_Length), &
                                     elemR(JMidx+6,er_Length))

        !% HACK: finding the average breadth. This will not work for channels with other than rectangular geometry.
        !% we need to generalize this
        elemSGR(JMidx,eSGR_Rectangular_Breadth) = (elemR(JMidx+1,er_Length)*elemSGR(JMidx+1,eSGR_Rectangular_Breadth) + &
                                                   elemR(JMidx+2,er_Length)*elemSGR(JMidx+2,eSGR_Rectangular_Breadth) + &
                                                   elemR(JMidx+3,er_Length)*elemSGR(JMidx+3,eSGR_Rectangular_Breadth) + &
                                                   elemR(JMidx+4,er_Length)*elemSGR(JMidx+4,eSGR_Rectangular_Breadth) + &
                                                   elemR(JMidx+5,er_Length)*elemSGR(JMidx+5,eSGR_Rectangular_Breadth) + &
                                                   elemR(JMidx+6,er_Length)*elemSGR(JMidx+6,eSGR_Rectangular_Breadth))/ &   
                                                   elemR(JMidx,er_Length)

        !% Volume
        !% rectangular volume depends on characteristic length and breadth.
        elemR(JMidx,er_Volume) =   elemSGR(JMidx,eSGR_Rectangular_Breadth) * elemR(JMidx,er_Length) * elemR(JMidx,er_Depth)
                                                                                    
        !BRHbugfix 20210813 elemR(JBidx,er_Volume_N0) = elemR(JMidx,er_Volume)
        !BRHbugfix 20210813 elemR(JBidx,er_Volume_N1) = elemR(JMidx,er_Volume)
        elemR(JMidx,er_Volume_N0) = elemR(JMidx,er_Volume) !BRHbugfix 20210813 
        elemR(JMidx,er_Volume_N1) = elemR(JMidx,er_Volume) !BRHbugfix 20210813 

    
        ! call the standard geometry update for junction branches
        call geo_assign_JB (ALLtm, ep_JM_ALLtm) !BRHbugfix 20210813    

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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
        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"


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
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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
        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

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
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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
        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"


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
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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
        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        if (setting%SmallVolume%UseSmallVolumes) then
            where (elemI(:,ei_geometryType) == rectangular)
                elemR(:,er_SmallVolume) = setting%SmallVolume%DepthCutoff * elemSGR(:,eSGR_Rectangular_Breadth) * &
                    elemR(:,er_Length)

            elsewhere (elemI(:,ei_geometryType) == trapezoidal)
                elemR(:,er_SmallVolume) = (elemSGR(:,eSGR_Trapezoidal_Breadth) + onehalfR * &
                    (elemSGR(:,eSGR_Trapezoidal_LeftSlope) + elemSGR(:,eSGR_Trapezoidal_RightSlope)) * &
                    setting%SmallVolume%DepthCutoff) * setting%SmallVolume%DepthCutoff
            endwhere
        else
            elemR(:,er_SmallVolume) = zeroR
        end if

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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
        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        elemR(1:size(elemR,1)-1,er_FlowrateLateral) = zeroR

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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
        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        elemR(1:size(elemR,1)-1,er_ones) = oneR

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine init_IC_oneVectors
    !
    !==========================================================================
    !==========================================================================
    !
end module initial_condition
