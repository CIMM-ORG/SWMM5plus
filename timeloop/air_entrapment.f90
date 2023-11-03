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

        call ae_detect_air_pockets ()

        call ae_initialize_air_pockets ()

    end subroutine air_entrapment_toplevel 
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine ae_detect_air_pockets () 
        !%------------------------------------------------------------------
        !% Description:
        !% Find air Pockets in an link
        !%------------------------------------------------------------------
            integer          :: ii, jj, startIdx, endIdx, airPocketIdx
            integer, pointer :: cIdx, nElem, eIdx(:), fUp(:), fDn(:)
            logical, pointer :: conAir, elemSur(:), faceSur(:) 
            logical          :: possibleAirPocket
        !%------------------------------------------------------------------
        !% cycle through the conduits to find air Pockets,
        do ii = 1,N_conduit
            !% pointers
            cIdx      => pConduitIdx(ii)
            nElem     => link%I(cIdx,li_N_element)
            conAir    => link%YN(cIdx,lYN_airPocketDetected)
            eIdx      => conduitElemMapsI(ii,1:nElem,cmi_elem_idx)
            fUp       => conduitElemMapsI(ii,1:nElem,cmi_elem_up_face)
            fDn       => conduitElemMapsI(ii,1:nElem,cmi_elem_dn_face)
            elemSur   => elemYN(:,eYN_isPSsurcharged)
            faceSur   => faceYN(:,fYN_isPSsurcharged)

            !% reset the possible air pocket detection logical
            possibleAirPocket = .false.

            !% ------------------------------------------------------------------
            !% initial air pockets detection
            !% search for mixed flow condition in a conduit
            !% if a conduit has a mixed of surcharged and non-surcharged cell,
            !% it will have to potential to entrap air pockets
            !% this reduces the need to cycle through all the conduits for
            !% possible air pocket
            if ((any(elemSur(eIdx))) .and. any(elemSur(eIdx))) then
                possibleAirPocket = .true.
            end if

            !% if initially any airPocket is detected, map the location of that air Pocket
            if (possibleAirPocket) then
                !% reset the counter for number of Pockets
                airPocketIdx = zeroI
                !% reset the indexes
                startIdx    = zeroI
                endIdx      = zeroI
                !% reset the air entrapment type
                airI(ii,:,airI_type) = noAirPocket
                !% reset the logicals
                airYN(ii,:,airYN_air_pocket_detected) = .false.
                !% initialize the index for the second do loop
                jj = oneI

                !% --- cycle through the conduit elements to 
                !%     find entrapped air pockets
                do while (jj <= nElem)
                    !% the element is not surcharged and start counting 
                    !% the non-surcharged elements
                    if (.not. elemSur(eIdx(jj))) then

                        !% set the starting index 
                        startIdx = jj

                        !% cycle through the next elements until a surcharge
                        !% element is encountered
                        do while (jj <= nElem .and. (.not. elemSur(eIdx(jj))))
                            jj = jj + oneI
                        end do

                        !% set the ending index
                        endIdx = jj - oneI
                        !% count the number of airpockets
                        airPocketIdx = airPocketIdx + oneI

                        !% --- HACK: 
                        !%     while surcharging, v-shape patterns can develop
                        !%     thus remove the saved entrapped airpockets that
                        !%     consnsts of only one element and both the upstream
                        !%     and downstream elements are surcharged (thus creating a v)
                        if ((startIdx - endIdx == zeroI)  .and. ((startIdx > oneI) .and. (endIdx < nElem)))    then
                            airPocketIdx = airPocketIdx - oneI
                            cycle
                        end if

                        !% --- air pocket limit
                        !%     if the airpocket limit doesnot exceed the permissible
                        !%     number of airpockets per conduit, store the data
                        if ((airPocketIdx > zeroI) .and. (airPocketIdx <= max_airpockets_per_conduit)) then
                            !% arrayI_pos selects the right column positions of airpockets
                            !% save the integer air pocket data
                            airI(ii,airPocketIdx,airI_idx)        = airPocketIdx
                            airI(ii,airPocketIdx,airI_elem_start) = eIdx(startIdx)
                            airI(ii,airPocketIdx,airI_face_up)    = fUp(startIdx)
                            airI(ii,airPocketIdx,airI_elem_end)   = eIdx(endIdx)
                            airI(ii,airPocketIdx,airI_face_dn)    = fDn(endIdx)

                            !% save the logical air pocket data
                            airYN(ii,airPocketIdx,airYN_air_pocket_detected) = .true.
                            !% save the airpocket detection at the conduit
                            conAir = .true.

                            !% --- set the type of the airpocket
                            !%     if the starting element is the first element in the conduit
                            !%     there will be an upstream release
                            if (startIdx == oneI) then
                                airI(ii,airPocketIdx,airI_type) = upReleaseAirpocket
                            
                            !%     else if the ending element is the last element in the conduit
                            !%     there will be a downstream release
                            else if (endIdx == nElem) then
                                airI(ii,airPocketIdx,airI_type) = dnReleaseAirpocket
                            
                            !%     else the pocket is entrapped
                            else
                                airI(ii,airPocketIdx,airI_type) = entrappedAirpocket
                            end if

                        !% if there is too many air pockets than permissible  
                        else if (airPocketIdx > max_airpockets_per_conduit) then
                            print*, "The conduit has, ", airPocketIdx, " airpockets "
                            print*, "which is more than maximum permissible of, ", max_airpockets_per_conduit
                            print*
                        end if

                    !% progress the counter for surcharge elements
                    else
                        jj = jj + oneI
                    end if
                end do

                !% debug printing
                ! if (any(airYN(ii,:,airYN_air_pocket_detected))) then
                !     print*, 'Entrapped air at conduit  ',link%names(ii)%str
                !     print*, 'airPocketIdx', airPocketIdx
                !     print*, 'idx 1     = ', airI(ii,1,airI_idx),         'idx 2     = ', airI(ii,2,airI_idx),        'idx 3     = ',airI(ii,3,airI_idx) 
                !     print*, 'type 1    = ', airI(ii,1,airI_type),        'type 2    = ', airI(ii,2,airI_type),       'type 3    = ', airI(ii,3,airI_type)
                !     print*, 'elem up 1 = ', airI(ii,1,airI_elem_start),  'elem up 2 = ', airI(ii,2,airI_elem_start), 'elem up 3 = ',airI(ii,3,airI_elem_start) 
                !     print*, 'elem dn 1 = ', airI(ii,1,airI_elem_end),    'elem dn 2 = ', airI(ii,2,airI_elem_end),   'elem dn 3 = ',airI(ii,3,airI_elem_end) 
                !     print* 
                ! end if
            end if
        end do

    end subroutine ae_detect_air_pockets
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine ae_initialize_air_pockets ()
        !%------------------------------------------------------------------
        !% Description:
        !% Find the initial volume when an air pocket is detected
        !%------------------------------------------------------------------
            integer          :: ii, jj, startIdx, endIdx, airPocketIdx
            integer, pointer :: cIdx, eStart, eEnd, fUp, fDn
            real(8), pointer :: airVol, flowUp, flowDn
            real(8), pointer :: elemVol(:), fullVol(:), faceFlow(:)
            logical, pointer :: conAir, airPocket
            logical          :: possibleAirPocket
        !%------------------------------------------------------------------
        !% static pointers
        elemVol   => elemR(:,er_Volume)
        fullVol   => elemR(:,er_FullVolume)
        faceFlow  => faceR(:,fr_Flowrate)
        !% cycle through the conduits to find air Pockets,
        do ii = 1,N_conduit
            !% pointers
            cIdx      => pConduitIdx(ii)
            conAir    => link%YN(cIdx,lYN_airPocketDetected)

            if (conAir) then
                do jj = 1,max_airpockets_per_conduit 
                    eStart    => airI(ii,jj,airI_elem_start)
                    eEnd      => airI(ii,jj,airI_elem_end)
                    fUp       => airI(ii,jj,airI_face_up)
                    fDn       => airI(ii,jj,airI_face_dn)
                    airVol    => airR(ii,jj,airR_volume)
                    flowUp    => airR(ii,jj,airR_flowUp)
                    flowDn    => airR(ii,jj,airR_flowDn)
                    airPocket => airYN(ii,jj,airYN_air_pocket_detected)
                    
                    !% calculate the air volume and 
                    !% save the water flowrate at the airpocket interface
                    if (airPocket) then
                        airVol = max(sum(fullVol(eStart:eEnd) - elemVol(eStart:eEnd)), zeroR)
                        flowUp = faceFlow(fUp)
                        flowDn = faceFlow(fDn)
                    else
                        airVol = zeroR
                        flowUp = zeroR
                        flowDn = zeroR
                    end if
                end do
            end if


            !% debug printing
            if (any(airYN(ii,:,airYN_air_pocket_detected))) then
                print*, 'Entrapped air at conduit  ',link%names(ii)%str
                print*, 'Vol 1     = ', airR(ii,1,airR_volume),         'Vol 2     = ', airR(ii,2,airR_volume),        'Vol 3     = ',airR(ii,3,airR_volume)
                
            end if
                
        end do
        

    end subroutine ae_initialize_air_pockets
    !%    
    !%==========================================================================
    !%==========================================================================
    !%

    
end module air_entrapment