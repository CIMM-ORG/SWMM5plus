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
    !use pack_mask_arrays

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

        integer         :: ii

        character(64)   :: subroutine_name = 'init_IC_setup'

    !--------------------------------------------------------------------------
        if (setting%Debug%File%initial_condition) print *, '*** leave ',subroutine_name

        !% get data that can be extracted from links
        call init_IC_from_linkdata ()

        !% get data that can be extracted from nodes
        ! call initial_condition_from_nodedata ()

        !% set up all the static packs and masks
        ! call pack_mask_arrays_all ()

        !% set up all the dynamic packs and masks
        ! call pack_dynamic_arrays ()

        !% update all the auxiliary variables
        ! call update_auxiliary_variables

        !% update faces
        ! call face_update()

        sync all

        if (setting%Debug%File%initial_condition) then
            !% only using the first processor to print results
            if (this_image() == 1) then

                do ii = 1,num_images()
                   print*, '----------------------------------------------------'
                   print*, 'image = ', ii
                   print*, '..................elements..........................'
                   print*, elemR(:,er_Depth)[ii], 'Depth'
                   call execute_command_line('')
                enddo

            endif
        endif

        if (setting%Debug%File%initial_condition) print *, '*** leave ',subroutine_name
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
        if (setting%Debug%File%initial_condition) print *, '*** leave ',subroutine_name

        !% Setting the local image value
        image = this_image()

        !% pack all the link indexes in an image
        packed_link_idx = pack(linkI(:,li_idx), (linkI(:,li_P_image) .eq. image))

        !% find the number of links in an image
        pLink = size(packed_link_idx)

        !% cycle through the links in an image
        do ii = 1,pLink
            !% necessary pointers
            thisLink    => packed_link_idx(ii)
            
            call init_IC_get_depth_from_linkdata (thisLink)

            call init_IC_get_flow_roughness_from_linkdata (thisLink)

            call init_IC_get_geometry_from_linkdata (thisLink)

            !% we need another call here to set the type of 
            !% time march to be used.

        end do

        if (setting%Debug%File%initial_condition) print *, '*** leave ',subroutine_name
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
        if (setting%Debug%File%initial_condition) print *, '*** leave ',subroutine_name

        !% type of initial depth type
        LdepthType  => linkI(thisLink,li_InitialDepthType)

        !% up and downstream depths on this link
        DepthUp => linkR(thisLink,lr_InitialUpstreamDepth)
        DepthDn => linkR(thisLink,lr_InitialDnstreamDepth)

        !% set the depths in link elements from links
        select case (LdepthType)

            case (Uniform)

                !%  if the link has a uniform depth as an initial condition
                if (linkR(thisLink,lr_InitialDepth) .ne. nullvalueR) then
                    
                    where (elemI(:,ei_link_Gidx_SWMM) .eq. thisLink)
                        elemR(:,er_Depth) = linkR(thisLink,lr_InitialDepth)
                    endwhere
                else
                    where (elemI(:,ei_link_Gidx_SWMM) .eq. thisLink)
                        elemR(:,er_Depth) = onehalfR * (DepthUp + DepthDn)
                    endwhere
                endif

            case (LinearlyVarying)

                !% if the link has linearly-varying depth
                !% depth at the upstream element (link position = 1)
                where ( (elemI(:,ei_link_Pos) .eq. 1) .and. (elemI(:,ei_link_Gidx_SWMM) .eq. thisLink) )
                    elemR(:,er_Depth) = DepthUp
                endwhere

                !%  using a linear distribution over the links
                ei_max = maxval(elemI(:,ei_link_Pos), 1, elemI(:,ei_link_Gidx_SWMM) .eq. thisLink)

                do mm=2,ei_max
                    !% find the element that is at the mm position in the link
                    where ( (elemI(:,ei_link_Pos) .eq. mm) .and. (elemI(:,ei_link_Gidx_SWMM) .eq. thisLink) )
                        !% use a linear interpolation
                        elemR(:,er_Depth) = DepthUp - (DepthUp - DepthDn) * real(mm - oneI) / real(ei_max - oneI)
                    endwhere
                end do

            case (ExponentialDecay)

                !% if the link has exponentially decayed depth
                !% depth at the upstream element (link position = 1)
                where ( (elemI(:,ei_link_Pos) .eq. 1) .and. (elemI(:,ei_link_Gidx_SWMM) .eq. thisLink) )
                    elemR(:,er_Depth) = DepthUp
                endwhere

                !% find the remaining elements in the link
                ei_max = maxval(elemI(:,ei_link_Pos), 1, elemI(:,ei_link_Gidx_SWMM) .eq. thisLink)

                do mm=2,ei_max
                    kappa = real(mm - oneI)

                    !%  depth decreases exponentially going downstream
                    if (DepthUp - DepthDn > zeroR) then
                        where ( (elemI(:,ei_link_Pos)       .eq. mm      ) .and. &
                                (elemI(:,ei_link_Gidx_SWMM) .eq. thisLink) )
                            elemR(:,er_Depth) = DepthUp - (DepthUp - DepthDn) * exp(-kappa)
                        endwhere

                    !%  depth increases exponentially going downstream
                    elseif (DepthUp - DepthDn < zeroR) then
                        where ( (elemI(:,ei_link_Pos)       .eq. mm      ) .and. &
                                (elemI(:,ei_link_Gidx_SWMM) .eq. thisLink) )
                            elemR(:,er_Depth) = DepthUp + (DepthDn - DepthUp) * exp(-kappa)
                        endwhere

                    !%  uniform depth
                    else
                        where ( (elemI(:,ei_link_Pos)       .eq. mm      ) .and. &
                                (elemI(:,ei_link_Gidx_SWMM) .eq. thisLink) )
                            elemR(:,er_Depth) = DepthUp
                        endwhere
                    endif
                end do

            case default
                print*, 'In ', subroutine_name
                print*, 'error: unexpected initial depth type, ', LdepthType,'  in link, ', thisLink
                stop

        end select

        if (setting%Debug%File%initial_condition) print *, '*** leave ',subroutine_name
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
        if (setting%Debug%File%initial_condition) print *, '*** leave ',subroutine_name

        !%  handle all the initial conditions that don't depend on geometry type
        where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
            elemR(:,er_Flowrate)       = linkR(thisLink,lr_InitialFlowrate)
            elemR(:,er_Flowrate_N0)    = linkR(thisLink,lr_InitialFlowrate)
            elemR(:,er_Flowrate_N1)    = linkR(thisLink,lr_InitialFlowrate)
            elemR(:,er_Roughness)      = linkR(thisLink,lr_Roughness)
        endwhere

        if (setting%Debug%File%initial_condition) print *, '*** leave ',subroutine_name
    end subroutine init_IC_get_flow_roughness_from_linkdata
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
        if (setting%Debug%File%initial_condition) print *, '*** leave ',subroutine_name

        !% necessary pointers
        linkType      => linkI(thisLink,li_link_type)

        select case (linkType)

            case (lChannel)
                !% get geomety data for channels
                call init_IC_get_channel_geometry (thisLink)

            case (lpipe)
                !% get geomety data for pipes
                call init_IC_get_pipe_geometry (thisLink)

            case (lweir)
                !% get geomety data for weirs
                call init_IC_get_weir_geometry(thisLink)

            case (lOrifice)

                print*, 'In ', subroutine_name
                print*, 'orifices are not handeled yet'
                stop

            case (lPump)

                print*, 'In ', subroutine_name
                print*, 'pumps are not handeled yet'
                stop

            case default

                print*, 'In ', subroutine_name
                print*, 'error: unexpected link, ', linkType,'  in the network'
                stop

        end select
        

        if (setting%Debug%File%initial_condition) print *, '*** leave ',subroutine_name
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
        if (setting%Debug%File%initial_condition) print *, '*** leave ',subroutine_name

        !% pointer to geometry type
        geometryType => linkI(thisLink,li_geometry)

        select case (geometryType)

        case (lRectangular)

            where (elemI(:,ei_link_Gidx_SWMM) == thisLink)

                elemI(:,ei_geometryType) = eRectangular
                elemR(:,er_BreadthMax)   = linkR(thisLink,lr_BreadthScale)
                elemR(:,er_Area)         = elemR(:,er_BreadthMax) * elemR(:,er_Depth)
                elemR(:,er_Area_N0)      = elemR(:,er_Area)
                elemR(:,er_Area_N1)      = elemR(:,er_Area)
                elemR(:,er_Volume)       = elemR(:,er_Area) * elemR(:,er_Length)
                elemR(:,er_Volume_N0)    = elemR(:,er_Volume)
                elemR(:,er_Volume_N1)    = elemR(:,er_Volume)
                elemR(:,er_ZbreadthMax)  = linkR(thisLink,lr_FullDepth)

                !% store geometry specific data
                elemSGR(:,eSGR_Rectangular_Breadth) = linkR(thisLink,lr_BreadthScale)
            endwhere

        case default

            print*, 'In, ', subroutine_name
            print*, 'Only rectangular channel geometry is handeled at this moment'
            stop

        end select

        if (setting%Debug%File%initial_condition) print *, '*** leave ',subroutine_name
    end subroutine init_IC_get_channel_geometry
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_IC_get_pipe_geometry (thisLink)
    !--------------------------------------------------------------------------
    !
    !% get the geometry data for pipe links 
    !% and calculate element volumes
    !
    !--------------------------------------------------------------------------

        integer, intent(in) :: thisLink
        integer, pointer    :: geometryType 

        character(64) :: subroutine_name = 'init_IC_get_pipe_geometry'
    !--------------------------------------------------------------------------
        if (setting%Debug%File%initial_condition) print *, '*** leave ',subroutine_name

        !% pointer to geometry type
        geometryType => linkI(thisLink,li_geometry)

        select case (geometryType)

        case (lRectangular)

            where (elemI(:,ei_link_Gidx_SWMM) == thisLink)

                elemI(:,ei_geometryType)    = eRectangular
                
                elemR(:,er_BreadthMax)      = linkR(thisLink,lr_BreadthScale)
                elemR(:,er_Area)            = elemR(:,er_BreadthMax) * elemR(:,er_Depth)
                elemR(:,er_Area_N0)         = elemR(:,er_Area)
                elemR(:,er_Area_N1)         = elemR(:,er_Area)
                elemR(:,er_Volume)          = elemR(:,er_Area) * elemR(:,er_Length)
                elemR(:,er_Volume_N0)       = elemR(:,er_Volume)
                elemR(:,er_Volume_N1)       = elemR(:,er_Volume)
                elemR(:,er_FullDepth)       = linkR(thisLink,lr_FullDepth)
                elemR(:,er_Zcrown)          = elemR(:,er_Zbottom) + elemR(:,er_FullDepth)
                elemR(:,er_FullArea)        = elemR(:,er_BreadthMax) * elemR(:,er_FullDepth)
                elemR(:,er_FullVolume)      = elemR(:,er_FullArea) * elemR(:,er_Length)

                elemYN(:,eYN_canSurcharge)  =  .true.

                !% store geometry specific data
                elemSGR(:,eSGR_Rectangular_Breadth) = linkR(thisLink,lr_BreadthScale)
            endwhere

        case default

            print*, 'In, ', subroutine_name
            print*, 'Only rectangular pipe geometry is handeled at this moment'
            stop

        end select

        if (setting%Debug%File%initial_condition) print *, '*** leave ',subroutine_name
    end subroutine init_IC_get_pipe_geometry
    !
    !==========================================================================
    !==========================================================================
    !
     subroutine init_IC_get_weir_geometry (thisLink)
    !--------------------------------------------------------------------------
    !
    !% get the geometry and other data data for weir links 
    !% and calculate element volumes
    !
    !--------------------------------------------------------------------------

        integer, intent(in) :: thisLink
        integer, pointer    :: specificWeirType 

        character(64) :: subroutine_name = 'init_IC_get_weir_geometry'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !% pointer to specific weir type
        specificWeirType => linkI(thisLink,li_weir_type)

        select case (specificWeirType)
            !% copy weir specific data
            case (lTrapezoidalWeir) 

                where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
                    !% integer data
                    elemI(:,ei_geometryType)                 = eTrapezoidal
                    elemSI(:,eSi_Weir_SpecificType)          = eTrapezoidalWeir

                    !% real data
                    elemSR(:,eSr_Weir_EffectiveFullDepth)    = linkR(thisLink,lr_FullDepth)
                    elemSR(:,eSr_Weir_DischargeCoeff1)       = linkR(thisLink,lr_DischargeCoeff1)
                    elemSR(:,eSr_Weir_DischargeCoeff2)       = linkR(thisLink,lr_DischargeCoeff2)
                    elemSR(:,eSr_Weir_TrapezoidalBreadth)    = linkR(thisLink,lr_BreadthScale)
                    elemSR(:,eSr_Weir_TrapezoidalLeftSlope)  = linkR(thisLink,lr_LeftSlope)
                    elemSR(:,eSr_Weir_TrapezoidalRightSlope) = linkR(thisLink,lr_RightSlope)
                    elemSR(:,eSr_Weir_Zcrest)                = elemR(:,er_Zbottom) + linkR(thisLink,lr_InletOffset)

                    !% HACK: I am not sure if we need to update the initial area or volume of an weir element
                    !% since they will all be set to zero values at the start of the simulation
                endwhere

            case (lSideFlowWeir) 

                where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
                    !% integer data
                    elemI(:,ei_geometryType)                 = eRectangular
                    elemSI(:,eSi_Weir_SpecificType)          = eSideFlowWeir
                    elemSI(:,eSi_Weir_EndContractions)       = linkI(thisLink,li_weir_EndContrations)

                    !% real data
                    elemSR(:,eSr_Weir_EffectiveFullDepth)    = linkR(thisLink,lr_FullDepth)
                    elemSR(:,eSr_Weir_DischargeCoeff2)       = linkR(thisLink,lr_DischargeCoeff2)
                    elemSR(:,eSr_Weir_RectangularBreadth)    = linkR(thisLink,lr_BreadthScale)
                    elemSR(:,eSr_Weir_Zcrest)                = elemR(:,er_Zbottom) + linkR(thisLink,lr_InletOffset)

                    !% HACK: I am not sure if we need to update the initial area or volume of an weir element
                    !% since they will all be set to zero values at the start of the simulation
                endwhere

            case (lRoadWayWeir)

                print*, 'In ', subroutine_name
                print*, 'roadway weir is not handeled yet'
                stop

            case (lVnotchWeir)

                where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
                    !% integer data
                    elemI(:,ei_geometryType)                 = eTriangular
                    elemSI(:,eSi_Weir_SpecificType)          = eVnotchWeir

                    !% real data
                    elemSR(:,eSr_Weir_EffectiveFullDepth)    = linkR(thisLink,lr_FullDepth)
                    elemSR(:,eSr_Weir_DischargeCoeff1)       = linkR(thisLink,lr_DischargeCoeff1)
                    elemSR(:,eSr_Weir_TriangularSideSlope)   = linkR(thisLink,lr_SideSlope)
                    elemSR(:,eSr_Weir_Zcrest)                = elemR(:,er_Zbottom) + linkR(thisLink,lr_InletOffset)

                    !% HACK: I am not sure if we need to update the initial area or volume of an weir element
                    !% since they will all be set to zero values at the start of the simulation
                endwhere

            case (lTransverseWeir)

                where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
                    !% integer data
                    elemI(:,ei_geometryType)                 = eRectangular
                    elemSI(:,eSi_Weir_SpecificType)          = eTransverseWeir
                    elemSI(:,eSi_Weir_EndContractions)       = linkI(thisLink,li_weir_EndContrations)

                    !% real data
                    elemSR(:,eSr_Weir_EffectiveFullDepth)    = linkR(thisLink,lr_FullDepth)
                    elemSR(:,eSr_Weir_DischargeCoeff2)       = linkR(thisLink,lr_DischargeCoeff2)
                    elemSR(:,eSr_Weir_RectangularBreadth)    = linkR(thisLink,lr_BreadthScale)
                    elemSR(:,eSr_Weir_Zcrest)                = elemR(:,er_Zbottom)  + linkR(thisLink,lr_InletOffset)

                    !% HACK: I am not sure if we need to update the initial area or volume of an weir element
                    !% since they will all be set to zero values at the start of the simulation
                endwhere

            case default

                print*, 'In ', subroutine_name
                print*, 'error: unknown weir type, ', specificWeirType,'  in network'
                stop

        end select
        
        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine init_IC_get_weir_geometry
    !
    !==========================================================================
    !==========================================================================
    !

    !
    !==========================================================================
    !==========================================================================
    !

    !
    !==========================================================================
    !==========================================================================
    !

    !
    !==========================================================================
    !==========================================================================
    !

    !
    !==========================================================================
    !==========================================================================
    !
end module initial_condition
