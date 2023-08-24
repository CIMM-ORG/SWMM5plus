module air_tracking
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Procedures for dynamic Preissmann Slot following paper by
    !% Sharior, Hodges, and Vasconcelos (2023)
    !%==========================================================================
    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use utility_crash

    implicit none

    private
    ! public :: 

    !%    
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    !%
    subroutine air_constricted_links (pCol_Closed) 
        !%------------------------------------------------------------------
        !% Description:
        !% Find if a any element in a link is surcharge constricting 
        !% the airflow. 
        !%------------------------------------------------------------------
        integer, intent(in) :: 
        integer          :: ii
        integer, pointer :: thisP(:), npack, fup, fdn, elemLinkIdx(:), linkIdx(:)
        logical, pointer :: isSurcharged(:), linkConstrained(:)
        character(64) :: subroutine_name = 'air_constricted_links'
        !%------------------------------------------------------------------
        npackP  = npack_elemP(pCol_Closed)
        if (npackP < 1) return 
        thisP => elemP(1:npackP,pCol)

        !% Aliases
        isSurcharged => elemYN(:,eYN_isSurcharged)
        elemLinkIdx  => elemI(:,ei_link_Gidx_BIPquick)
        linkIdx      => link%I(:,li_idx)
        linkConstrained => link%YN(:,lYN_isConstricted)

        !% go through the links and check if any of the 
        !% element in that link is closed and surcharged
        !% and thus constraining the airflow
        do ii = 1,N_Link
            if (any((elemLinkIdx(thisP) == ii) .and. isSurcharged(thisP))) then
                linkConstrained(ii) = .true.
            end if
        end do

    end subroutine air_constricted_links
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    !%    
    !%==========================================================================
    !%==========================================================================
    !%

    contains
end module air_tracking