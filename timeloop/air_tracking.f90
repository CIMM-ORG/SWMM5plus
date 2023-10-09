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

contains
    !%    
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    !%
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine at_constricted_links (pCol_Closed) 
        !%------------------------------------------------------------------
        !% Description:
        !% Find if a any element in a link is surcharge constricting 
        !% the airflow. 
        !%------------------------------------------------------------------
        integer, intent(in) :: pCol_Closed
        integer          :: ii
        integer, pointer :: thisP(:), npack, fup, fdn, elemLinkIdx(:), linkIdx(:)
        logical, pointer :: isSurcharged(:), linkSurcharged(:)
        character(64) :: subroutine_name = 'at_constricted_links'
        !%------------------------------------------------------------------
        npackP  = npack_elemP(pCol_Closed)
        if (npackP < 1) return 
        thisP => elemP(1:npackP,pCol)

        !% Aliases
        isSurcharged => elemYN(:,eYN_isSurcharged)
        elemLinkIdx  => elemI(:,ei_link_Gidx_BIPquick)
        linkIdx      => link%I(:,li_idx)
        linkSurcharged => link%YN(:,lYN_isSurcharged)

        !% go through the links and check if any of the 
        !% element in that link is closed and surcharged
        !% and thus constraining the airflow
        do ii = 1,N_Link
            !% reset air constriction
            linkSurcharged(ii) = .false.
            !% if there is a surcharge element of a link, set the link as constricted
            if (any((elemLinkIdx(thisP) == ii) .and. isSurcharged(thisP))) then
                linkSurcharged(ii) = .true.
            end if
        end do

    end subroutine at_constricted_links
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine at_elem_air_volume (pCol_Closed) 
        !%------------------------------------------------------------------
        !% Description:
        !% Find the airvolume from the empty spaces in closed elements
        !%------------------------------------------------------------------
        integer, intent(in) :: pCol_Closed
        integer          :: ii
        integer, pointer :: thisP(:), npack
        real(8), pointer :: Volume(:), AirVolume(:), fullVolume(:)

        character(64) :: subroutine_name = 'at_elem_air_volume'
        !%------------------------------------------------------------------
        npackP  = npack_elemP(pCol_Closed)
        if (npackP < 1) return 
        thisP => elemP(1:npackP,pCol)

        !% Aliases
        Volume     => elemR(:,er_Volume)
        AirVolume  => elemR(:,er_Air_volume)
        fullVomume => elemR(:,er_FullVolume)

        !% assume the empty space will be occupied by air
        AirVolume(thisP) = fullVolume(thisP) - Volume(thisP)

        where (AirVolume(thisP) <= zeroR)
            AirVolume(thisP) = zeroR
        end where

    end subroutine at_elem_air_volume   
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine at_constricted_link_air_volume (pCol_Closed) 
        !%------------------------------------------------------------------
        !% Description:
        !% Find if a any element in a link is surcharge constricting 
        !% the airflow. 
        !%------------------------------------------------------------------
        integer, intent(in) :: pCol_Closed
        integer          :: ii
        integer, pointer :: thisP(:), npack, fup, fdn, elemLinkIdx(:), linkIdx(:)
        integer, pointer :: linkUpElem(:), linkDnElem(:)
        real(8), pointer :: linkAirVol(:), elemAirVol
        logical, pointer :: linkUpSurcharged(:), linkDnSurcharged(:), linkSurcharged(:)
        character(64) :: subroutine_name = 'at_constricted_link_air_volume'
        !%------------------------------------------------------------------
        npackP  = npack_elemP(pCol_Closed)
        if (npackP < 1) return 
        thisP => elemP(1:npackP,pCol)

        !% Aliases
        isSurcharged => elemYN(:,eYN_isSurcharged)
        elemLinkIdx  => elemI(:,ei_link_Gidx_BIPquick)
        elemAirVol   => elemR(:,er_Air_volume)
        linkIdx      => link%I(:,li_idx)
        linkAirVol   => link%R(:,lr_Air_Volume)
        linkUpElem   => link%YN(:,li_first_elem_idx)
        linkDnElem   => link%YN(:,li_last_elem_idx)
        linkSurcharged   => link%YN(:,lYN_isSurcharged)
        linkUpSurcharged => link%YN(:,lYN_isUpSurcharge)
        linkDnSurcharged => link%YN(:,lYN_isDnSurcharge)

        do ii - 1,N_link
            !% reset if the link is surcharged 
            linkSurcharged(ii) = .false.





            if (linkSurcharged(ii)) then
                linkSurcharged(ii) = sum(elemAirVol(thisP), elemLinkIdx(thisP) == ii)
            end if
        end do


    end subroutine at_constricted_link_air_volume
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    !%    
    !%==========================================================================
    !%==========================================================================
    !%

    
end module air_tracking