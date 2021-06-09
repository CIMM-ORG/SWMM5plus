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

    public :: initial_condition_setup

    private

contains
    !
    !==========================================================================
    ! PUBLIC
    !==========================================================================
    !
    subroutine initial_condition_setup ()
    !--------------------------------------------------------------------------
    !
    !% set up the initial conditions for all the elements
    !
    !--------------------------------------------------------------------------

        integer         :: ii

        character(64)   :: subroutine_name = 'initial_condition_setup'

    !--------------------------------------------------------------------------
        if (setting%Debug%File%initial_condition) print *, '*** leave ',subroutine_name

        !% get data that can be extracted from links
        call initial_condition_from_linkdata ()

        !% get data that can be extracted from nodes
        ! call initial_condition_from_nodedata ()

        !% set up all the static packs and masks
        ! call pack_mask_arrays_all ()

        !% set up all the dynamic packs and masks
        ! call pack_dynamic_arrays ()

        !% update all the auxiliary variables
        ! call update_auxiliary_variables

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
    end subroutine initial_condition_setup
    !
    !==========================================================================
    ! PRIVATE
    !==========================================================================
    !
    subroutine initial_condition_from_linkdata ()
    !--------------------------------------------------------------------------
    !
    !% get the initial depth, flowrate, and geometry data from links
    !
    !--------------------------------------------------------------------------

        integer                                     :: ii, mm, image, pLink
        integer                                     :: ei_max
        real(8)                                     :: kappa
        integer, pointer                            :: thisLink, LdepthType 
        real(8), pointer                            :: DepthUp, DepthDn
        integer, dimension(:), allocatable, target  :: packed_link_idx

        character(64) :: subroutine_name = 'initial_condition_from_linkdata'
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
            LdepthType  => linkI(thisLink,li_InitialDepthType)

            !% up and downstream depths on this link
            DepthUp => linkR(thisLink,lr_InitialUpstreamDepth)
            DepthDn => linkR(thisLink,lr_InitialDnstreamDepth)

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

        end do

        if (setting%Debug%File%initial_condition) print *, '*** leave ',subroutine_name
    end subroutine initial_condition_from_linkdata
    !
    !==========================================================================
    !==========================================================================
    ! 
end module initial_condition
