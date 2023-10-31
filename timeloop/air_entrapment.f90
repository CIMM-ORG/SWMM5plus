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

        call ae_elem_air_volume ()

        call ae_detect_air_bubbles ()

    end subroutine air_entrapment_toplevel
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine ae_elem_air_volume () 
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
            logical, pointer :: upSur(:), dnSur(:), faceSurcharged(:) 
            logical, pointer :: elemSur(:), surcharged(:)
        !%------------------------------------------------------------------
        !% static pointers 
        fullVolume => elemR(:,er_FullVolume)
        volume     => elemR(:,er_Volume)
        surcharged => elemYN(:,eYN_isSurcharged)
        faceFlow   => faceR(:,fr_Flowrate)
        faceSurcharged => faceYN(:,fYN_isPSsurcharged)

        !% cycle through the conduits to find element air volumes
        do ii = 1,N_conduit
            !% additional pointers
            nElem     => conduitAirI(ii,cai_N_elements)
            eIdx      => elemAirI(ii,1:nElem,eai_elem_idx)
            fUp       => elemAirI(ii,1:nElem,eai_elem_up_face)
            fDn       => elemAirI(ii,1:nElem,eai_elem_dn_face)
            airVolume => elemAirR(ii,1:nElem,ear_air_volume)
            flowUp    => elemAirR(ii,1:nElem,ear_flowrate_up)
            flowDn    => elemAirR(ii,1:nElem,ear_flowrate_dn)
            elemSur   => elemAirYN(ii,1:nElem,eaYN_elem_pressurized)
            upSur     => elemAirYN(ii,1:nElem,eaYN_elem_up_pressurized)
            dnSur     => elemAirYN(ii,1:nElem,eaYN_elem_dn_pressurized)

            !% find the airvolumes in conduit elements
            airVolume = max(fullVolume(eIdx) - volume(eIdx), zeroR)
            !% find if the conduit element is surcharged
            elemSur   = surcharged(eIdx)
            !% copy the upstream and downstream flowrates from the faces
            flowUp    = faceFlow(fUp)
            flowDn    = faceFlow(fDn)
            !% find if upsteam or downstream of the element is constricted 
            !% due to water pressurization
            upSur    = faceSurcharged(fUp)
            dnSur    = faceSurcharged(fDn)

        end do

    end subroutine ae_elem_air_volume   
    ! !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine ae_detect_air_bubbles () 
        !%------------------------------------------------------------------
        !% Description:
        !% Find air bubbles in an link
        !%------------------------------------------------------------------
            integer          :: ii, jj, startIdx, endIdx
            integer, pointer :: nElem, eIdx(:), airPocketType(:)
            logical, pointer :: upSur(:), dnSur(:), faceSurcharged(:) 
            logical, pointer :: elemSur(:), surcharged(:), entrappedAir(:)
            logical          :: possibleAirBubble, stopCheck
        !%------------------------------------------------------------------
        stopCheck = .false.
        !% cycle through the conduits to find air bubbles,
        do ii = 1,N_conduit
            !% pointers
            nElem     => conduitAirI(ii,cai_N_elements)
            eIdx      => elemAirI(ii,:,eai_elem_idx)
            elemSur   => elemAirYN(ii,:,eaYN_elem_pressurized)
            upSur     => elemAirYN(ii,:,eaYN_elem_up_pressurized)
            dnSur     => elemAirYN(ii,:,eaYN_elem_dn_pressurized)
            entrappedAir  => elemAirYN(ii,:,eaYN_elem_has_entrapped_air)
            airPocketType => elemAirI(ii,:,eai_airpocket_type)

            !% reset the possible air pocket detection logical
            possibleAirBubble = .false.
            !% reset the air pocket type
            airPocketType(1:nElem) = noAirPocket
            !% reset the airpocket logical
            entrappedAir(1:nElem)  = .false.

            !% initial detection of air bubbles
            !% if a conduit has a mixed of surcharged and non-surcharged cell,
            !% it will have to potential to entrap air pockets
            if ((any(elemSur(1:nElem))) .and. any(.not. elemSur(1:nElem))) then
            
                possibleAirBubble = .true.
                ! print*, 'possible air bubble detected at conduit  ', link%names(ii)%str
                ! print*
                
            end if

            !% if initially any airbubble is detected, map the location of that air bubble
            if (possibleAirBubble) then
                !% reset the indexes
                startIdx    = nullvalueI
                endIdx      = nullvalueI
                !% --- cycle through the conduit elements to 
                !%     find entrapped air pockets
                jj = 1
                do while (jj <= nElem)

                    !% the element is not surcharged and start counting 
                    !% the non-surcharged elements
                    if (.not. elemSur(jj)) then
                        !% set the starting element 
                        startIdx = jj

                        !% cycle through the next elements until a surcharge
                        !% element is encountered
                        do while (jj <= nElem .and. (.not. elemSur(jj)))
                            jj = jj + 1
                        end do

                        !% set the ending index
                        endIdx = jj - 1

                        !% set the entrapped airpocket logical to true
                        entrappedAir(startIdx:endIdx) = .true.
                        !% set the type of the airpocket 
                        if (startIdx == 1) then
                            airPocketType(startIdx:endIdx) = upReleaseAirpocket
                        else if (endIdx == nElem) then
                            airPocketType(startIdx:endIdx) = dnReleaseAirpocket
                        else
                            airPocketType(startIdx:endIdx) = entrappedAirpocket
                        end if

                        !% HACK: 
                        !% while surcharging, v-shape patterns can develop
                        !% thus eliminate the entrapped airpockets that
                        !% consnsts of only one element and both the upstream
                        !% and downstream elements are surcharged 
                        if ((startIdx - endIdx == zeroI) .and. (airPocketType(startIdx) == entrappedAirpocket)) then
                            airPocketType(startIdx:endIdx) = noAirPocket
                            entrappedAir(startIdx:endIdx)  = .false.
                        end if

                    !% progress the counter for surcharge elements
                    else
                        jj = jj + 1
                    end if

                end do
                ! if (any(airPocketType(1:nElem) == entrappedAirpocket)) then
                !     print*, 'Entrapped air at conduit  ',link%names(ii)%str
                !     print*, entrappedAir(1:nElem), 'entrappedAir(1:nElem)'
                !     print*, elemSur(1:nElem), 'elemSur(1:nElem)'
                !     print*, airPocketType(1:nElem), 'airPocketType(1:nElem)'
                ! end if
            end if

        end do

    end subroutine ae_detect_air_bubbles
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    !%    
    !%==========================================================================
    !%==========================================================================
    !%

    
end module air_entrapment