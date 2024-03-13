module discretization
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Performs discretization into multiple elements per link   
    !%
    !%==========================================================================

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use utility_crash, only: util_crashpoint

    implicit none

    public init_discretization_nominal
    private

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine init_discretization_nominal(link_idx)
    !%----------------------------------------------------------------------
    !% Description:
    !%   This subroutine sets the number of elements per link.  The element length
    !%   is adjusted so that an integer number of elements is assigned to each link.
    !%----------------------------------------------------------------------
    !% Declarations
        integer, intent(in) :: link_idx
        real(8) :: remainder
        real(8), pointer :: elem_nominal_length
        integer, pointer :: min_elem_per_link
        logical, pointer :: use_nominal_length
        character(64) :: subroutine_name = 'init_discretization_nominal'
    !%----------------------------------------------------------------------
    !% Preliminaries
        if (setting%Debug%File%discretization) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !%----------------------------------------------------------------------
    !% Aliases
        !use_nominal_length  => setting%Discretization%UseNominalElemLength
        elem_nominal_length => setting%Discretization%NominalElemLength
        min_elem_per_link   => setting%Discretization%MinElementPerLink
    !%----------------------------------------------------------------------

        !% --- treatment of for special links
        !%     note that equivalent orifices have type lOrifice
        if ((link%I(link_idx,li_link_type) == lWeir)    .or. &
            (link%I(link_idx,li_link_type) == lOrifice) .or. &
            (link%I(link_idx,li_link_type) == lOutlet)  .or. &
            (link%I(link_idx,li_link_type) == lPump)           ) then
            link%I(link_idx, li_N_element) = oneI
            link%R(link_idx, lr_ElementLength) = link%R(link_idx, lr_Length)
            return
        end if 

        !% --- for channel and pipe only
        if  ((link%I(link_idx,li_link_type) == lChannel)  .or. &
             (link%I(link_idx,li_link_type) == lPipe)           ) then
        
            select case (setting%Discretization%Method)

                case (EqualElements)
                    !% --- Adjusts the number of elements in a link based on the length so
                    !%     that element lengths are close to the nominal length

                    !% --- find remainder after division
                    remainder = mod(link%R(link_idx,lr_Length), elem_nominal_length)
                    
                    if ( remainder == zeroR ) then
                        !% --- the elements fit evenly into the length of link
                        link%I(link_idx, li_N_element)     = int(link%R(link_idx, lr_Length) / elem_nominal_length)
                        link%R(link_idx, lr_ElementLength) = link%R(link_idx, lr_Length) / link%I(link_idx, li_N_element)

                    elseif ( remainder .ge. onehalfR * elem_nominal_length ) then
                        !% --- the remainder is greater than half of an element length so the ceiling value is used
                        !%     and the elements will be slightly shorter than the nominal length
                        link%I(link_idx, li_N_element)     = ceiling(link%R(link_idx,lr_Length) / elem_nominal_length)
                        link%R(link_idx, lr_ElementLength) = link%R(link_idx, lr_Length) / link%I(link_idx, li_N_element)

                    else
                        !% --- the remainder is less than half of an element length so floor value is used and
                        !%     the elements will be slightly longer than the nominal length
                        link%I(link_idx, li_N_element)     = max(floor(link%R(link_idx,lr_Length) / elem_nominal_length), oneI)
                        link%R(link_idx, lr_ELementLength) = link%R(link_idx, lr_Length) / link%I(link_idx, li_N_element)

                    end if

                    !% --- Additional check to ensure that every link has at least one element
                    if ( link%R(link_idx, lr_Length) .le. elem_nominal_length ) then
                        link%I(link_idx, li_N_element) = oneI
                        link%R(link_idx, lr_ElementLength) = link%R(link_idx, lr_Length)
                    end if

                    if ((link%I(link_idx, li_N_element)     < min_elem_per_link)                .or. &
                        (link%R(link_idx, lr_ElementLength) < (onehalfR * elem_nominal_length))        ) then

                        select case (setting%Discretization%SmallElementHandling)

                            case (EquivalentOrifice,LengthenLink,FailLimiter) 
                                !% --- this should not be rearched as initialization should fix the problem
                                write(*,*) 'CODE ERROR: in Small Element Handling'
                                write(*,*) 'Link with insufficient length found in discretization'
                                write(*,*) 'but this should have been fixed in link initialization'
                                write(*,*) 'Link # is ',link_idx 
                                write(*,*) 'Link Name is ',trim(link%Names(link_idx)%str)
                                write(*,*) 'Link length    ',link%R(link_idx, lr_Length)
                                write(*,*) 'Element length ',link%R(link_idx, lr_ElementLength)
                                write(*,*) 'nominal length ',elem_nominal_length
                                write(*,*) 
                                call util_crashpoint(72120987)

                            case (AllowSmallLinks)
                                !% --- subdivide link into smaller elements to meet minimum
                                link%I(link_idx, li_N_element) = min_elem_per_link
                                link%R(link_idx, lr_ElementLength) = link%R(link_idx, lr_Length) / real(min_elem_per_link,real(8))

                            case default 
                                write(*,*) 'CODE ERROR: unexpected case default'
                                call util_crashpoint(2209874)

                        end select
                    else 
                        !% --- continue, sufficient elements per link
                    end if

                    ! !% --- only force a minimum number of elements per link if the minimum
                    ! !%     number of elements per link settings is set to greater than one
                    ! if (min_elem_per_link > oneI) then
                    !     !% check for links that has less elements than 
                    !     !% the specified minimum number of elements
                    !     if (link%I(link_idx, li_N_element) < min_elem_per_link) then
                    !         link%I(link_idx, li_N_element) = min_elem_per_link
                    !         link%R(link_idx, lr_ElementLength) = link%R(link_idx, lr_Length) / link%I(link_idx, li_N_element)
                    !     end if 
                    ! end if

                case (UnequalElements)
                    !% --- use the minimum number of elements per link in each link
                    !%     This results in different element sizes throughout system.
                    link%I(link_idx, li_N_element)     = min_elem_per_link
                    link%R(link_idx, lr_ElementLength) = link%R(link_idx, lr_Length) / real(min_elem_per_link,real(8))

                    if (link%R(link_idx, lr_ElementLength) < (onehalfR * elem_nominal_length)) then

                        select case (setting%Discretization%SmallElementHandling)
                            case (EquivalentOrifice,LengthenLink,FailLimiter) 
                                !% --- this should not be reached as initialization should fix the problem
                                write(*,*) 'CODE ERROR: in Small Element Handling'
                                write(*,*) 'Link with insufficient length found in discretization'
                                write(*,*) 'but this should have been fixed in link initialization'
                                write(*,*) 'Link # is ',link_idx 
                                write(*,*) 'Link Name is ',trim(link%Names(link_idx)%str)
                                write(*,*) 'Link length    ',link%R(link_idx, lr_Length)
                                write(*,*) 'Element length ',link%R(link_idx, lr_ElementLength)
                                write(*,*) 'nominal length ',elem_nominal_length
                                write(*,*) 
                                call util_crashpoint(829887)

                            case (AllowSmallLinks)
                                !% --- subdivide link into smaller elements to meet minimum
                                link%I(link_idx, li_N_element) = min_elem_per_link
                                link%R(link_idx, lr_ElementLength) = link%R(link_idx, lr_Length) / real(min_elem_per_link,real(8))

                            case default 
                                write(*,*) 'CODE ERROR: unexpected case default'
                                call util_crashpoint(2209874)
                        end select
                    else
                        !% --- continue, small element length not found
                    end if

                case default
                    write(*,*) 'CODE ERROR: unexpected case default'
                    write(*,*) 'missing handling of a link type with key #',link%I(link_idx,li_link_type)
                    write(*,*) trim(reverseKey(link%I(link_idx,li_link_type)))
                    call util_crashpoint(81109872)
            end select

        else 
            write(*,*) 'CODE ERROR -- unexpected else reached'
        end if

    end subroutine init_discretization_nominal
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_discretization_adjustlinklength()

        !% ARCHIVE -- this may be recalled in future code

        print *, 'OBSOLETE 20230507 brh'
        stop 539874
        !%-------------------------------------------------------------------------
        !% Description:
        !%   This subroutine computes the "Adusted_Length" that is a section of the
        !%   a channel/conduit that is used for a junction branch, which is completed
        !%   in init_network_nJm_branch_length()
        !% 
        !%   If "isAdjustLinkLength = .true., then this value is subtracted
        !%   from the channel/conduit length herein 
        !%-------------------------------------------------------------------------
        !% Declarations
        !     integer          :: ii, Adjustment_flag
        !     real(8)          :: temp_length
        !     logical, pointer :: isAdjustLinkLength
        !     real(8)          :: elem_nominal_length, elem_shorten_cof
        !     character(64)    :: subroutine_name = 'init_discretization_adjustlinklength'
        ! !%-------------------------------------------------------------------------
        ! if (setting%Debug%File%discretization) &
        !     write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        ! isAdjustLinkLength => setting%Discretization%AdjustLinkLengthForJunctionBranchYN

        ! ! if (isAdjustLinkLength) then 
        ! !     print *, 'CONFIGURATION ERROR'
        ! !     print *, 'setting.Discretization.AdjustLinkLengthForJunctionBranchYN = .true.'
        ! !     print *, 'SWMM5+ presently requires this to be .false.'
        ! !     print *, 'Please change the setting in your *.json file'
        ! !     stop 77098723
        ! ! end if
        
        ! do ii =1, N_link
        !     !% --- default shorting coefficient (reset for each link)
        !     elem_nominal_length = link%R(ii,lr_Length) / link%I(ii,li_N_element)
        !     elem_shorten_cof    = setting%Discretization%JunctionBranchLengthFactor

        !     temp_length = link%R(ii,lr_Length) ! length of link ii
        !     Adjustment_flag = oneI
           
        !     !% --- adjust the shortening for small link lengths 20220520brh 
        !     ! if (temp_length < (oneR + twoR*elem_shorten_cof) * elem_nominal_length) then
        !     !     !% --- limit the shortening to 1/4 of the total element length (1/8 on either side)
        !     !     elem_shorten_cof = temp_length / (elem_nominal_length * eightR)
        !     ! end if

        !     if ( node%I(link%I(ii,li_Mnode_u), ni_node_type) == nJm ) then
        !         temp_length = temp_length - elem_shorten_cof * elem_nominal_length ! make a cut for upstream M junction
        !         Adjustment_flag = Adjustment_flag + oneI
        !     end if

        !     if ( node%I(link%I(ii,li_Mnode_d), ni_node_type) == nJm ) then
        !         temp_length = temp_length - elem_shorten_cof * elem_nominal_length ! make a cut for downstream M junction
        !         Adjustment_flag = Adjustment_flag + oneI
        !     end if

        !     if ((link%I(ii,li_link_type) == lChannel) .or. (link%I(ii,li_link_type) == lPipe)) then
        !         link%I(ii,li_length_adjusted) = Adjustment_flag
        !         link%R(ii,lr_AdjustedLength) = temp_length
        !         !% set the new element length based on the adjusted link length if the user permits
        !         if (isAdjustLinkLength) link%R(ii,lr_ElementLength) = link%R(ii,lr_AdjustedLength)/link%I(ii,li_N_element)     
        !     else
        !         link%R(ii,lr_AdjustedLength) = link%R(ii,lr_ElementLength)
        !         link%I(ii,li_length_adjusted) = DiagAdjust
        !     end if
        ! end do

        ! if (setting%Debug%File%discretization)  &
        !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_discretization_adjustlinklength
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
!%
end module discretization
