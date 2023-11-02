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

        ! call ae_elem_air_volume ()

        call ae_detect_air_bubbles ()

    end subroutine air_entrapment_toplevel
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    ! subroutine ae_elem_air_volume () 
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Find the airvolume from the empty spaces in closed elements
    !     !% Copy the flows from upstream and downstream of the elements
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer          :: ii 
    !         integer, pointer :: nElem, eIdx(:), fUp(:), fDn(:)
    !         real(8), pointer :: airVolume(:), volume(:), fullVolume(:)
    !         real(8), pointer :: flowUp(:), flowDn(:), faceFlow(:)
    !         logical, pointer :: upSur(:), dnSur(:), faceSurcharged(:) 
    !         logical, pointer :: elemSur(:), surcharged(:)
    !     !%------------------------------------------------------------------
    !     !% static pointers 
    !     fullVolume => elemR(:,er_FullVolume)
    !     volume     => elemR(:,er_Volume)
    !     surcharged => elemYN(:,eYN_isSurcharged)
    !     faceFlow   => faceR(:,fr_Flowrate)
    !     faceSurcharged => faceYN(:,fYN_isPSsurcharged)

    !     !% cycle through the conduits to find element air volumes
    !     do ii = 1,N_conduit
    !         !% additional pointers
    !         nElem     => airI(ii,cairI_N_elements)
    !         eIdx      => conduitElemMapsI(ii,1:nElem,cmi_elem_idx)
    !         fUp       => conduitElemMapsI(ii,1:nElem,cmi_elem_up_face)
    !         fDn       => conduitElemMapsI(ii,1:nElem,cmi_elem_dn_face)
    !         airVolume => elemAirR(ii,1:nElem,ear_air_volume)
    !         flowUp    => elemAirR(ii,1:nElem,ear_flowrate_up)
    !         flowDn    => elemAirR(ii,1:nElem,ear_flowrate_dn)
    !         elemSur   => elemAirYN(ii,1:nElem,eairYN_elem_pressurized)
    !         upSur     => elemAirYN(ii,1:nElem,eairYN_elem_up_pressurized)
    !         dnSur     => elemAirYN(ii,1:nElem,eairYN_elem_dn_pressurized)

    !         !% find the airvolumes in conduit elements
    !         airVolume = max(fullVolume(eIdx) - volume(eIdx), zeroR)
    !         !% find if the conduit element is surcharged
    !         elemSur   = surcharged(eIdx)
    !         !% copy the upstream and downstream flowrates from the faces
    !         flowUp    = faceFlow(fUp)
    !         flowDn    = faceFlow(fDn)
    !         !% find if upsteam or downstream of the element is constricted 
    !         !% due to water pressurization
    !         upSur    = faceSurcharged(fUp)
    !         dnSur    = faceSurcharged(fDn)

    !     end do

    ! end subroutine ae_elem_air_volume   
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine ae_detect_air_bubbles () 
        !%------------------------------------------------------------------
        !% Description:
        !% Find air bubbles in an link
        !%------------------------------------------------------------------
            integer          :: ii, jj, startIdx, endIdx, airPocketIdx
            integer, pointer :: cIdx, nElem, eIdx(:), fUp(:), fDn(:)
            logical, pointer :: elemSur(:), faceSur(:)
            logical          :: possibleAirBubble
        !%------------------------------------------------------------------
        !% cycle through the conduits to find air bubbles,
        do ii = 1,N_conduit
            !% pointers
            cIdx      => pConduitIdx(ii)
            nElem     => link%I(cIdx,li_N_element)
            eIdx      => conduitElemMapsI(ii,1:nElem,cmi_elem_idx)
            fUp       => conduitElemMapsI(ii,1:nElem,cmi_elem_up_face)
            fDn       => conduitElemMapsI(ii,1:nElem,cmi_elem_dn_face)
            elemSur   => elemYN(:,eYN_isPSsurcharged)
            faceSur   => faceYN(:,fYN_isPSsurcharged)


            !% reset the possible air pocket detection logical
            possibleAirBubble = .false.

            !% ------------------------------------------------------------------
            !% initial air pockets detection
            !% search for mixed flow condition in a conduit
            !% if a conduit has a mixed of surcharged and non-surcharged cell,
            !% it will have to potential to entrap air pockets
            !% this reduces the need to cycle through all the conduits for
            !% possible air pocket
            if ((any(elemSur(eIdx))) .and. any(elemSur(eIdx))) then
                possibleAirBubble = .true.
            end if

            !% if initially any airbubble is detected, map the location of that air bubble
            if (possibleAirBubble) then
                !% reset the counter for number of bubbles
                airPocketIdx = zeroI
                !% reset the indexes
                startIdx    = zeroI
                endIdx      = zeroI
                !% reset the air entrapment type
                airI(ii,:,airI_type) = noAirPocket
                !% reset the logicals
                airYN(ii,:,airYN_is_air_pocket) = .false.
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
                            airYN(ii,airPocketIdx,airYN_is_air_pocket) = .true.

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
                ! if (any(airYN(ii,:,airYN_is_air_pocket))) then
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