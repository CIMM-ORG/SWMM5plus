module air_entrapment
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

    
    public :: air_entrapment_toplevel 

contains
    !%    
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    !%
    subroutine air_entrapment_toplevel (istep)
        !%-----------------------------------------------------------------
        !% Description
        !%-----------------------------------------------------------------
        !% toplevel subroutine for air entrapment modeling
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: istep
            integer             ::  ii
        !%------------------------------------------------------------------

        call ae_elem_air_preliminaries ()

    end subroutine air_entrapment_toplevel
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    ! subroutine at_constricted_links (pCol_Closed) 
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Find if a any element in a link is surcharge constricting 
    !     !% the airflow. 
    !     !%------------------------------------------------------------------
    !     integer, intent(in) :: pCol_Closed
    !     integer          :: ii
    !     integer, pointer :: thisP(:), npack, fup, fdn, elemLinkIdx(:), linkIdx(:)
    !     logical, pointer :: isSurcharged(:), linkSurcharged(:)
    !     character(64) :: subroutine_name = 'at_constricted_links'
    !     !%------------------------------------------------------------------
    !     npackP  = npack_elemP(pCol_Closed)
    !     if (npackP < 1) return 
    !     thisP => elemP(1:npackP,pCol)

    !     !% Aliases
    !     isSurcharged => elemYN(:,eYN_isSurcharged)
    !     elemLinkIdx  => elemI(:,ei_link_Gidx_BIPquick)
    !     linkIdx      => link%I(:,li_idx)
    !     linkSurcharged => link%YN(:,lYN_isSurcharged)

    !     !% go through the links and check if any of the 
    !     !% element in that link is closed and surcharged
    !     !% and thus constraining the airflow
    !     do ii = 1,N_Link
    !         !% reset air constriction
    !         linkSurcharged(ii) = .false.
    !         !% if there is a surcharge element of a link, set the link as constricted
    !         if (any((elemLinkIdx(thisP) == ii) .and. isSurcharged(thisP))) then
    !             linkSurcharged(ii) = .true.
    !         end if
    !     end do

    ! end subroutine at_constricted_links
    ! !%    
    ! !%==========================================================================
    ! !%==========================================================================
    ! !%
    subroutine ae_elem_air_preliminaries () 
        !%------------------------------------------------------------------
        !% Description:
        !% Find the airvolume from the empty spaces in closed elements
        !% Copy the flows from upstream and downstream of the elements
        !%------------------------------------------------------------------
        !% Declarations:
            integer          :: ii 
            integer, pointer :: nElem, eIdx(:), fUp(:), fDn(:)
            real(8), pointer :: airVolume(:), volume(:), fullVolume(:)
            real(8), pointer :: flowUp(:), flowDn(:), faceFlow(:)
            logical, pointer :: upCons(:), dnCons(:), faceSurcharged(:) 
            logical, pointer :: elemCons(:), surcharged(:)
        !%------------------------------------------------------------------
        !% static pointers 
        fullVolume => elemR(:,er_FullVolume)
        volume     => elemR(:,er_Volume)
        surcharged => elemYN(:,eYN_isSurcharged)
        faceFlow   => faceR(:,fr_Flowrate)
        faceSurcharged => faceYN(:,fYN_isPSsurcharged)

        !% cycle through the links to find element air volumes
        do ii = 1,N_link
            !% additional pointers
            nElem     => link%I(ii,li_N_element)
            eIdx      => LinkElemMapsI(ii,1:nElem,lmi_elem_idx)
            fUp       => LinkElemMapsI(ii,1:nElem,lmi_elem_up_face)
            fDn       => LinkElemMapsI(ii,1:nElem,lmi_elem_dn_face)
            airVolume => elemAirR(ii,1:nElem,ear_air_volume)
            flowUp    => elemAirR(ii,1:nElem,ear_flowrate_up)
            flowDn    => elemAirR(ii,1:nElem,ear_flowrate_dn)
            elemCons  => elemAirYN(ii,1:nElem,eaYN_elem_pressurized)
            upCons    => elemAirYN(ii,1:nElem,eaYN_elem_up_constricted)
            dnCons    => elemAirYN(ii,1:nElem,eaYN_elem_dn_constricted)

            !% find the airvolumes in conduit elements
            airVolume = max(fullVolume(eIdx) - volume(eIdx), zeroR)
            !% find if the conduit element is surcharged
            elemCons  = surcharged(eIdx)
            !% copy the upstream and downstream flowrates from the faces
            flowUp    = faceFlow(fUp)
            flowDn    = faceFlow(fDn)
            !% find if upsteam or downstream of the element is constricted 
            !% due to water pressurization
            upCons    = faceSurcharged(fUp)
            dnCons    = faceSurcharged(fDn)

        end do

    end subroutine ae_elem_air_preliminaries   
    ! !%    
    ! !%==========================================================================
    ! !%==========================================================================
    ! !%
    ! subroutine at_constricted_link_air_volume (pCol_Closed) 
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Find if a any element in a link is surcharge constricting 
    !     !% the airflow. 
    !     !%------------------------------------------------------------------
    !     integer, intent(in) :: pCol_Closed
    !     integer          :: ii
    !     integer, pointer :: thisP(:), npack, fup, fdn, elemLinkIdx(:), linkIdx(:)
    !     integer, pointer :: linkUpElem(:), linkDnElem(:)
    !     real(8), pointer :: linkAirVol(:), elemAirVol
    !     logical, pointer :: linkUpSurcharged(:), linkDnSurcharged(:), linkSurcharged(:)
    !     character(64) :: subroutine_name = 'at_constricted_link_air_volume'
    !     !%------------------------------------------------------------------
    !     npackP  = npack_elemP(pCol_Closed)
    !     if (npackP < 1) return 
    !     thisP => elemP(1:npackP,pCol)

    !     !% Aliases
    !     isSurcharged => elemYN(:,eYN_isSurcharged)
    !     elemLinkIdx  => elemI(:,ei_link_Gidx_BIPquick)
    !     elemAirVol   => elemR(:,er_Air_volume)
    !     linkIdx      => link%I(:,li_idx)
    !     linkAirVol   => link%R(:,lr_Air_Volume)
    !     linkUpElem   => link%YN(:,li_first_elem_idx)
    !     linkDnElem   => link%YN(:,li_last_elem_idx)
    !     linkSurcharged   => link%YN(:,lYN_isSurcharged)
    !     linkUpSurcharged => link%YN(:,lYN_isUpSurcharge)
    !     linkDnSurcharged => link%YN(:,lYN_isDnSurcharge)

    !     do ii - 1,N_link
    !         !% reset if the link is surcharged 
    !         linkSurcharged(ii) = .false.





    !         if (linkSurcharged(ii)) then
    !             linkSurcharged(ii) = sum(elemAirVol(thisP), elemLinkIdx(thisP) == ii)
    !         end if
    !     end do


    ! end subroutine at_constricted_link_air_volume
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    !%    
    !%==========================================================================
    !%==========================================================================
    !%

    
end module air_entrapment