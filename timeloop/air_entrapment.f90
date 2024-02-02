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
    use face, only: face_interpolation
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
            integer             :: ii
        !%------------------------------------------------------------------

        if (setting%AirTracking%StaticAirPocket) then
            !% find the airpocket in the network 
            !% (only one airpocket per superlink is allowed)
            call airpockets_detection_single ()
        else
            !% find the airpockets in the network
            call airpockets_detection ()
        end if

        !% find the pressure head due to airpocket
        call airpockets_calculation (istep)
        !% interpolate the faces again after air calculation
        !% to update the new heads to the faces (only head interp)
        call face_interpolation(fp_noBC_IorS, .false., .true., .false., .true., .true.)

    end subroutine air_entrapment_toplevel 
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine airpockets_calculation (istep)
        !%------------------------------------------------------------------
        !% Description:
        !% air pocket calculation
        !%------------------------------------------------------------------
            integer, intent(in) :: istep
            integer          :: ii, jj
            integer, pointer :: cIdx(:), max_airpockets
            logical, pointer :: conAir(:)
        !%------------------------------------------------------------------
        !% initialize if an element has pressurized air
        elemR(1:size(elemR,1)-1,er_Pressurized_Air)   =  zeroR
        elemR(1:size(elemR,1)-1,er_Air_Pressure_Head) =  zeroR
        elemYN(1:size(elemYN,1)-1,eYN_hasAirPocket)   =  .false.

        !% cycle through the continious conduits to find air Pockets,
        do ii = 1, N_super_conduit
            !% pointers
            cIdx      => sc_link_Idx(ii,1:links_per_sc(ii))
            conAir    => link%YN(:,lYN_airPocketDetected)

            max_airpockets => setting%AirTracking%NumberOfAirpocketsAllowed

            !% only go through the conduit airpocket calculation if the
            !% conduit contains one
            if (any(conAir(cIdx))) then
 
                ! print*, '------------------------------------------------------'

                do jj = 1,max_airpockets 
    
                    !% initialize the airpocket real values
                    call airpocket_initialization (ii,jj, istep)
      
                    !% calculate the net flowrate though the airpockets
                    call airpocket_netflowrate (ii,jj)
    
                    !% calculate air outflow rate form an airpocket if any vent is present
                    call airpocket_air_mass_outflow (ii, jj, istep)

                    !% update the air density at the air pocket
                    call airpocket_air_density_update (ii, jj, istep)  

                    !% calculate the air pressure head
                    call airpocket_pressure_head (ii, jj, istep)

                    !% add the heads back to the elements
                    call add_airpocket_heads_to_elem (ii, jj)

                end do
            
            end if
                
        end do
        

    end subroutine airpockets_calculation
        !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine airpockets_detection () 
        !%------------------------------------------------------------------
        !% Description:
        !% Find air Pockets in an link
        !%------------------------------------------------------------------
            integer          :: ii, jj, kk, startIdx, endIdx, airPocketIdx, nElem
            integer, pointer :: cIdx(:), eIdx(:), fUp(:), fDn(:), max_airpockets
            logical, pointer :: conAir(:), elemSur(:), fBlocked(:)
            logical          :: possibleAirPocket
        !%------------------------------------------------------------------
        !% cycle through the conduits to find air Pockets,
        do ii = 1,N_super_conduit
            !% pointers
            cIdx      => sc_link_Idx(ii,1:links_per_sc(ii))
            nElem     =  sum(link%I(cIdx,li_N_element))
            conAir    => link%YN(:,lYN_airPocketDetected)
            eIdx      => conduitElemMapsI(ii,1:nElem,cmi_elem_idx)
            fUp       => conduitElemMapsI(ii,1:nElem,cmi_elem_up_face)
            fDn       => conduitElemMapsI(ii,1:nElem,cmi_elem_dn_face)
            elemSur   => elemYN(:,eYN_isPSsurcharged)
            fBlocked  => faceYN(:,fYN_isAirflowBlocked)

            max_airpockets  => setting%AirTracking%NumberOfAirpocketsAllowed

            !% reset the possible air pocket detection logical
            possibleAirPocket = .false.
            !% reset the conduit air pocket detection logical
            conAir(cIdx) = .false.

            !% initialize the conduitElemMapsI at each detection for moving air pocket
            conduitElemMapsI(ii,:,cmi_airpocket_idx)  = nullvalueI
            conduitElemMapsI(ii,:,cmi_airpocket_type) = noAirPocket

            !% ------------------------------------------------------------------
            !% initial air pockets screening
            !% search for if any faces conduit elements that has blocked airflow
            if (any(elemSur(eIdx))) then
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
                !% initialize the index for the second do loop
                jj = oneI

                !% --- cycle through the conduit elements to 
                !%     find entrapped air pockets
                do while (jj <= nElem)
                    !% --- Find the airpocket based on mixed condition first
                    !%     the element is not surcharged and start counting 
                    !%     the non-surcharged elements
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
                        !%     consnsts of only one element 
                        if ((startIdx - endIdx == zeroI))    then
                            airPocketIdx = airPocketIdx - oneI
                            cycle
                        end if

                        !% --- air pocket limit
                        !%     if the airpocket limit doesnot exceed the permissible
                        !%     number of airpockets per conduit, store the data
                        if ((airPocketIdx > zeroI) .and. (airPocketIdx <= max_airpockets)) then
                
                            if (.not. airYN(ii,airPocketIdx,airYN_air_pocket_detected)) then
                                airYN(ii,airPocketIdx,airYN_new_air_pocket) =  .true.
                            end if
                            !% arrayI_pos selects the right column positions of airpockets
                            !% save the integer air pocket data
                            airI(ii,airPocketIdx,airI_idx)        = airPocketIdx
                            airI(ii,airPocketIdx,airI_elem_start) = eIdx(startIdx)
                            airI(ii,airPocketIdx,airI_face_up)    = fUp(startIdx)
                            airI(ii,airPocketIdx,airI_elem_end)   = eIdx(endIdx)
                            airI(ii,airPocketIdx,airI_face_dn)    = fDn(endIdx)

                            !% --- set the type of the airpocket
                            !%     if the starting element is the first element in the conduit
                            !%     there will be an upstream release
                            if (startIdx == oneI) then
                                airI(ii,airPocketIdx,airI_type) = upReleaseAirpocket
                            
                            !% else if the ending element is the last element in the conduit
                            !% there will be a downstream release
                            else if (endIdx == nElem) then
                                airI(ii,airPocketIdx,airI_type) = dnReleaseAirpocket
                            
                            !% else the pocket is entrapped
                            else
                                airI(ii,airPocketIdx,airI_type) = entrappedAirpocket
                            end if

                            !% save airpocket detection data in the conduitElemMapsI array
                            conduitElemMapsI(ii,startIdx:endIdx,cmi_airpocket_idx)  = airPocketIdx
                            conduitElemMapsI(ii,startIdx:endIdx,cmi_airpocket_type) = airI(ii,airPocketIdx,airI_type)

                            !% save the airpocket detection logical
                            airYN(ii,airPocketIdx,airYN_air_pocket_detected) = .true.

                            !% save the airpocket detection at the conduit
                            conAir(cIdx) = .true.
                            !% --- special conditions for upstream and downstream 
                            !%     airflow blockage
                            if (airI(ii,airPocketIdx,airI_type) == upReleaseAirpocket)  then

                                !% if the airflow of upstream release is blocked 
                                if (fBlocked(fUp(oneI))) then
                                    airI(ii,airPocketIdx,airI_type)  = entrappedAirpocket
                                end if
                            end if
                                
                            if (airI(ii,airPocketIdx,airI_type) == dnReleaseAirpocket) then

                                !% if the airflow of downstream release is blocked
                                if (fBlocked(fDn(nElem))) then
                                    airI(ii,airPocketIdx,airI_type)  = entrappedAirpocket
                                end if
                            end if

                        !% else if there is too many air pockets than permissible  
                        else if (airPocketIdx > max_airpockets) then
                            print*, "The conduit has, ", airPocketIdx, " airpockets "
                            print*, "which is more than maximum permissible of, ", max_airpockets
                            print*
                        end if

                    !% progress the counter for surcharge elements
                    else
                        jj = jj + oneI
                    end if
                end do

                !% reset the rest of the empty airpocket array locations
                if (airPocketIdx <= max_airpockets) then
                    airYN(ii,airPocketIdx+1:max_airpockets,:)                 = .false.
                    airI(ii,airPocketIdx+1:max_airpockets,:)                  = nullvalueI
                    !% reset the air entrapment type
                    airI(ii,airPocketIdx+1:max_airpockets,airI_type)          = noAirPocket
                    airR(ii,airPocketIdx+1:max_airpockets,:)                  = zeroR
                    airR(ii,airPocketIdx+1:max_airpockets,airR_absolute_head) = setting%AirTracking%AtmosphericPressureHead
                    airR(ii,airPocketIdx+1:max_airpockets,airR_density)       = setting%AirTracking%AirDensity
                end if

            !% if there is not any possible air pockets, reset the values of the air arrays of that correspondig conduit
            else
                conAir(cIdx)  = .false.
                airYN(ii,:,:) = .false.
                airI(ii,:,:)  = nullvalueI
                !% reset the air entrapment type
                airI(ii,:,airI_type) = noAirPocket
                airR(ii,:,:)  = zeroR
                airR(ii,:,airR_absolute_head) = setting%AirTracking%AtmosphericPressureHead
                airR(ii,:,airR_density)       = setting%AirTracking%AirDensity
                conduitElemMapsI(ii,:,cmi_airpocket_idx)  = nullvalueI
                conduitElemMapsI(ii,:,cmi_airpocket_type) = noAirPocket
            end if

            !% debug printing
            ! if (any(airYN(ii,:,airYN_air_pocket_detected))) then
            !     print*, 'air at conduit indexes ',cIdx, ' ', link%names(cIdx(1))%str
            !     print*, 'airPocketIdx', airPocketIdx
            !     print*, 'idx 1     = ', airI(ii,1,airI_idx),         'idx 2     = ', airI(ii,2,airI_idx),        'idx 3     = ',airI(ii,3,airI_idx) 
            !     print*, 'type 1    = ', airI(ii,1,airI_type),        'type 2    = ', airI(ii,2,airI_type),       'type 3    = ', airI(ii,3,airI_type)
            !     print*, 'elem up 1 = ', airI(ii,1,airI_elem_start),  'elem up 2 = ', airI(ii,5,airI_elem_start), 'elem up 3 = ',airI(ii,3,airI_elem_start) 
            !     print*, 'elem dn 1 = ', airI(ii,1,airI_elem_end),    'elem dn 2 = ', airI(ii,2,airI_elem_end),   'elem dn 3 = ',airI(ii,3,airI_elem_end) 
            !     print*, 'Head 1= ', airR(ii,1,airR_absolute_head), 'Head 2= ', airR(ii,2,airR_absolute_head), 'Head 3= ', airR(ii,3,airR_absolute_head)
            !     print* 
            ! end if

        end do

    end subroutine airpockets_detection
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine airpocket_initialization (sc_Idx,aIdx,istep)
        !%------------------------------------------------------------------
        !% Description:
        !%   newly detected airpocket: initialize the inflow, outflow,  
        !%   volume and heads
        !%   already detected airpockets: initialize the inflow and outflow 
        !%   no air pockets: zero the values 
        !%------------------------------------------------------------------
            integer, intent (in) :: sc_Idx, aIdx, istep
            integer, pointer     :: elemStartIdx, elemEndIdx, faceUp, faceDn
            real(8), pointer     :: airVolume, inflow, outflow, airDensity
            real(8), pointer     :: airMass, absHead, gaugeHead, atmHead, rho_a
            real(8), pointer     :: elemVol(:), fullVol(:), faceFlow(:)
            logical, pointer     :: isAirPocket, newAirPocket

            integer, allocatable :: pElem(:)
        !%------------------------------------------------------------------
        !% Aliases
            elemStartIdx  => airI(sc_Idx,aIdx,airI_elem_start)
            elemEndIdx    => airI(sc_Idx,aIdx,airI_elem_end)
            faceUp        => airI(sc_Idx,aIdx,airI_face_up)
            faceDn        => airI(sc_Idx,aIdx,airI_face_dn)
            airVolume     => airR(sc_Idx,aIdx,airR_volume)
            airDensity    => airR(sc_Idx,aIdx,airR_density)
            airMass       => airR(sc_Idx,aIdx,airR_mass)
            inflow        => airR(sc_Idx,aIdx,airR_inflow)
            outflow       => airR(sc_Idx,aIdx,airR_outflow)
            absHead       => airR(sc_Idx,aIdx,airR_absolute_head)
            gaugeHead     => airR(sc_Idx,aIdx,airR_gauge_head)
            isAirPocket   => airYN(sc_Idx,aIdx,airYN_air_pocket_detected)
            newAirPocket  => airYN(sc_Idx,aIdx,airYN_new_air_pocket)

        !% other aliases
            elemVol   => elemR(:,er_Volume)
            fullVol   => elemR(:,er_FullVolume)
            faceFlow  => faceR(:,fr_Flowrate)
            atmHead   => setting%AirTracking%AtmosphericPressureHead
            rho_a     => setting%AirTracking%AirDensity

        !% if a airpocket is detected, initialize the air pocket
        if (isAirPocket) then
            
            !% pack the elemnts in this air pocket
            pElem = pack(conduitElemMapsI(sc_Idx,:,cmi_elem_idx), conduitElemMapsI(sc_Idx,:,cmi_airpocket_idx)  == aIdx)

            !% in this algorithm, aiVolume = 0 means the air pocket is newly detected
            if (newAirPocket) then
                newAirPocket = .false.
                gaugeHead     = zeroR
                airDensity    = rho_a
                absHead       = atmHead
                !% calculate the air pocket volume from empty space
                airVolume  = max(sum(fullVol(pElem) - elemVol(pElem)), zeroR)
                !% for a new pocket formulation save the older values
                airMass      = airDensity * airVolume
                !% copy over the flow data
                inflow     = faceFlow(faceUp)
                outflow    = faceFlow(faceDn) 
            else
                !% calculate the air pocket volume from empty space
                airVolume  = max(sum(fullVol(pElem) - elemVol(pElem)), zeroR)
                !% copy over the flow data
                inflow     = faceFlow(faceUp)
                outflow    = faceFlow(faceDn)
            end if

            !% deallocate the temporary array
            deallocate(pElem)
            
        else
            !% if there is not an airpocket present, zero out the values
            airVolume   = zeroR
            airMass     = zeroR
            inflow      = zeroR
            outflow     = zeroR
            gaugeHead   = zeroR
            airDensity  = rho_a
            absHead     = atmHead
        end if

    end subroutine airpocket_initialization
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine airpocket_netflowrate (sc_Idx, aIdx)
        !%------------------------------------------------------------------
        !% Description:
        !%   calculate the net flowrate through an airpocket 
        !%------------------------------------------------------------------
            integer, intent (in) :: sc_Idx, aIdx
            real(8), pointer     :: dvdt, inflow, outflow, faceFlow(:)
            logical, pointer     :: isAirPocket
        !%------------------------------------------------------------------
        !% Aliases
            dvdt        => airR(sc_Idx,aIdx,airR_dvdt)
            inflow      => airR(sc_Idx,aIdx,airR_inflow)
            outflow     => airR(sc_Idx,aIdx,airR_outflow)
            isAirPocket => airYN(sc_Idx,aIdx,airYN_air_pocket_detected)
            faceFlow    => faceR(:,fr_Flowrate)
        
        !% if a airpocket is detected, calculate the net flowrate
        if (isAirPocket) then
            dvdt = outflow - inflow
        else
            dvdt = zeroR
        end if

    end subroutine airpocket_netflowrate
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine airpocket_air_mass_outflow (sc_Idx, aIdx, istep)
        !%------------------------------------------------------------------
        !% Description:
        !%   calculates the air outflow rate for airpocket with release
        !%------------------------------------------------------------------
            integer, intent (in) :: sc_Idx, aIdx, istep
            real(8)              :: ExpansionFactor, ratio, areaOpening
            integer, pointer     :: pocketType, elemStartIdx, elemEndIdx
            real(8), pointer     :: massOutflow, absHead, dt, kappa, atmHead, airDensity
            real(8), pointer     :: dishCoeff, rho_w, rho_a, grav, minVentArea, airMass
            logical, pointer     :: isAirPocket
        !%------------------------------------------------------------------
        !% Aliases
            pocketType   => airI(sc_Idx,aIdx,airI_type)
            elemStartIdx => airI(sc_Idx,aIdx,airI_elem_start)
            elemEndIdx   => airI(sc_Idx,aIdx,airI_elem_end)
            absHead      => airR(sc_Idx,aIdx,airR_absolute_head)
            massOutflow  => airR(sc_Idx,aIdx,airR_mass_flowrate)
            airDensity   => airR(sc_Idx,aIdx,airR_density)
            airMass      => airR(sc_Idx,aIdx,airR_mass)
            isAirPocket  => airYN(sc_Idx,aIdx,airYN_air_pocket_detected)

            dt           => setting%Time%Hydraulics%Dt
            grav         => setting%constant%gravity
            kappa        => setting%AirTracking%PolytropicExponent
            atmHead      => setting%AirTracking%AtmosphericPressureHead
            dishCoeff    => setting%AirTracking%AirDischargeCoefficient
            rho_w        => setting%AirTracking%WaterDensity
            rho_a        => setting%AirTracking%AirDensity
            minVentArea  => setting%AirTracking%MinimumVentArea

        !% calculate the new air reseale rate from the opening
        !% calculate the new airpocket pressure head
        if (isAirPocket) then
            select case (pocketType)

                case (entrappedAirpocket)

                    areaOpening = zeroR

                case (upReleaseAirpocket)
                    
                    areaOpening = min(minVentArea, max(elemR(elemStartIdx,er_FullArea) - elemR(elemStartIdx,er_Area), zeroR))

                    ! areaOpening = minVentArea
                    
                case (dnReleaseAirpocket)

                    areaOpening = min(minVentArea, max(elemR(elemEndIdx,er_FullArea) - elemR(elemEndIdx,er_Area), zeroR))

                    ! areaOpening = minVentArea
                    

                case (noAirPocket)

                    areaOpening = zeroR

                case default
                    print*, 'ERROR: this must not reach'
                    stop 8413354

            end select

            !% calculate the ratio between atmospheric vs. absolute pressue head at conduit
            ratio = atmHead / absHead

            if (absHead >= 1.893 * atmHead) then
                !% chocked air expulsion form the valve
                if (airDensity > zeroR) then
                    massOutflow = - dishCoeff * areaOpening * airDensity * sqrt(grav * (rho_w / airDensity) &
                                * absHead) * sqrt(kappa * (twoR / (kappa + oneR)) ** ((kappa + oneR) / (kappa - oneR)))
                else
                    massOutflow = zeroR
                end if

            else if (absHead > atmHead .and. absHead < 1.893 * atmHead) then
                !% subsonic air expulsion from the valve
                if (airDensity > zeroR) then
                    !% find the expansion factor
                    ExpansionFactor = sqrt((kappa / (kappa - oneR)) * (ratio ** (twoR / kappa)) &
                                        * ((oneR - ratio ** ((kappa - oneR) / kappa)) / (oneR - ratio)))
                    !% find the airflow
                    massOutflow = - dishCoeff * areaOpening * ExpansionFactor * airDensity &
                                * sqrt(twoR * grav * (rho_w / airDensity) * abs(absHead - atmHead))
                else
                    massOutflow = zeroR
                end if 
        
            else
                massOutflow = zeroR
            end if
        end if


    end subroutine airpocket_air_mass_outflow
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine airpocket_pressure_head (sc_Idx, aIdx, istep)
        !%------------------------------------------------------------------
        !% Description:
        !%   Solves the airpocket pressure head
        !%------------------------------------------------------------------
            integer, intent (in) :: sc_Idx, aIdx, istep
            real                 :: alpha, beta
            integer, pointer     :: pocketType
            real(8), pointer     :: airVolume, airMass, dvdt, absHead, gaugeHead
            real(8), pointer     :: absHead_N0, dmdt
            real(8), pointer     :: dt, kappa, atmHead, theta, crk(:)
            logical, pointer     :: isAirPocket
        !%------------------------------------------------------------------
        !% Aliases
            pocketType  => airI(sc_Idx,aIdx,airI_type)
            airVolume   => airR(sc_Idx,aIdx,airR_volume)
            airMass     => airR(sc_Idx,aIdx,airR_mass)
            absHead     => airR(sc_Idx,aIdx,airR_absolute_head)
            absHead_N0  => airR(sc_Idx,aIdx,airR_absolute_head_N0)
            gaugeHead   => airR(sc_Idx,aIdx,airR_gauge_head)
            dvdt        => airR(sc_Idx,aIdx,airR_dvdt)
            dmdt        => airR(sc_Idx,aIdx,airR_mass_flowrate)
            isAirPocket => airYN(sc_Idx,aIdx,airYN_air_pocket_detected)
            crk         => setting%Solver%crk2
            dt          => setting%Time%Hydraulics%Dt
            kappa       => setting%AirTracking%PolytropicExponent
            atmHead     => setting%AirTracking%AtmosphericPressureHead
            theta       => setting%AirTracking%theta
        
        !% calculate the new airpocket pressure head
        if (isAirPocket) then

            if (airVolume > zeroR .and. airMass > zeroR) then
                !% find the alpha and beta
                alpha = (kappa / airVolume) * dvdt
                beta  = (kappa / airMass)   * dmdt

                !% find the absolute head 
                absHead = absHead * (oneR + dt * onehalfR * (oneR - theta) * (- alpha + beta)) / (oneR - dt * onehalfR * theta * (- alpha + beta))
                
                !% find the gauge pressure head
                gaugeHead = absHead - atmHead
            else
                isAirPocket = .false.
                absHead   = atmHead
                gaugeHead = absHead - atmHead  
            end if

        end if
    
    end subroutine airpocket_pressure_head
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine airpocket_air_density_update (sc_Idx, aIdx, istep)
        !%------------------------------------------------------------------
        !% Description:
        !%   Update the air density at the airpocket
        !%------------------------------------------------------------------
            integer, intent (in) :: sc_Idx, aIdx, istep
            integer, pointer     :: pocketType
            real(8), pointer     :: airVolume, massOutflow, airMass, airMass_N0
            real(8), pointer     :: airDensity_N0, dvdt, airDensity, dt, rho_a, crk(:)
            logical, pointer     :: isAirPocket, colAirpocket
            real                 :: dRho
        !%------------------------------------------------------------------
        !% Aliases
            pocketType    => airI(sc_Idx,aIdx,airI_type)
            airVolume     => airR(sc_Idx,aIdx,airR_volume)
            dvdt          => airR(sc_Idx,aIdx,airR_dvdt)
            airMass       => airR(sc_Idx,aIdx,airR_mass)
            airMass_N0    => airR(sc_Idx,aIdx,airR_mass_N0)
            airDensity    => airR(sc_Idx,aIdx,airR_density)
            airDensity_N0 => airR(sc_Idx,aIdx,airR_density_N0)
            massOutflow   => airR(sc_Idx,aIdx,airR_mass_flowrate)
            isAirPocket   => airYN(sc_Idx,aIdx,airYN_air_pocket_detected)
            colAirpocket  => airYN(sc_Idx,aIdx,airYN_air_pocket_collapsed)
            dt            => setting%Time%Hydraulics%Dt
            crk           => setting%Solver%crk2
            rho_a         => setting%AirTracking%AirDensity

        !% update the new density based on mass outflow
        if (isAirPocket) then

            !% timemarch air mass through mass flowrate
            airMass = airMass + dt * onehalfR * massOutflow

            !% limit air mass
            airMass = max(airMass,zeroR)
            
            if (airMass > zeroR) then
                colAirpocket = .false.
                airDensity = airMass / airVolume
            else 
                airDensity  = zeroR
                !% no mass left thus the airpocket has been collasped
                isAirPocket = .false.
                colAirpocket = .true.
            end if

        end if

    end subroutine airpocket_air_density_update
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine airpockets_volume_update (sc_Idx, aIdx, istep)
        !%------------------------------------------------------------------
        !% Description:
        !%   Timemarch airpocket volume
        !%------------------------------------------------------------------
            integer, intent (in) :: sc_Idx, aIdx, istep
            real(8)              :: tempVol
            real(8), pointer     :: airVolume, dvdt, airVolume_N0, dt, crk(:)
            logical, pointer     :: isAirPocket
        !%------------------------------------------------------------------
        !% Aliases
            airVolume    => airR(sc_Idx,aIdx,airR_volume)
            dvdt         => airR(sc_Idx,aIdx,airR_dvdt)
            airVolume_N0 => airR(sc_Idx,aIdx,airR_volume_N0)
            isAirPocket  => airYN(sc_Idx,aIdx,airYN_air_pocket_detected)
            dt           => setting%Time%Hydraulics%Dt
            crk          => setting%Solver%crk2
        
        !% update the new volume based on inflow and outflow outflow
        if (isAirPocket) then

            airVolume = airVolume_N0 + crk(istep) * dt * dvdt
            !% limit air volume
            airVolume = max(airVolume, zeroR)

        end if

    end subroutine airpockets_volume_update
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine add_airpocket_heads_to_elem (sc_Idx, aIdx)
        !%------------------------------------------------------------------
        !% Description:
        !%  add the airpocket heads to corresponding elements
        !%------------------------------------------------------------------
        integer, intent (in) :: sc_Idx, aIdx
        integer, pointer     :: elemStartIdx, elemEndIdx
        real(8), pointer     :: gaugeHead, elemHead(:), hasAirPocket(:), AirHead(:)
        real(8), pointer     :: fPnumber(:), Pnumber(:), dt
        logical, pointer     :: isAirPocket
        integer, allocatable :: pElem(:)
        !%------------------------------------------------------------------
        !% Aliases
            fPnumber     => faceR(:,fr_Preissmann_Number)
            Pnumber      => elemR(:,er_Preissmann_Number)
            elemStartIdx => airI(sc_Idx,aIdx,airI_elem_start)
            elemEndIdx   => airI(sc_Idx,aIdx,airI_elem_end)
            gaugeHead    => airR(sc_Idx,aIdx,airR_gauge_head)
            isAirPocket  => airYN(sc_Idx,aIdx,airYN_air_pocket_detected)
            elemHead     => elemR(:,er_Head)
            AirHead      => elemR(:,er_Air_Pressure_Head)
            hasAirPocket => elemR(:,er_Pressurized_Air)
            dt           => setting%Time%Hydraulics%Dt

        !% add the gauge heads to elem heads
        if (isAirPocket) then
            !% pack the elements that contain air pocket
            pElem = pack(conduitElemMapsI(sc_Idx,:,cmi_elem_idx), conduitElemMapsI(sc_Idx,:,cmi_airpocket_idx)  == aIdx)
            !% add the air head to piezometric head
            elemHead(pElem)     = max(elemHead(pElem) + gaugeHead, elemHead(pElem))
            AirHead(pElem)      = gaugeHead
            hasAirPocket(pElem) = oneR
            elemYN(pElem,eYN_hasAirPocket) = .true.

            !% deallocate the temporary array
            deallocate(pElem)
        end if

    end subroutine add_airpocket_heads_to_elem
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine airpockets_detection_single () 
        !%------------------------------------------------------------------
        !% Description:
        !% Only one air pocket permitted per link
        !%------------------------------------------------------------------
            integer          :: ii, jj, kk, startIdx, endIdx, airPocketIdx, nElem
            integer, pointer :: cIdx(:), eIdx(:), fUp(:), fDn(:), max_airpockets
            logical, pointer :: conAir(:), elemSur(:), fBlocked(:), StaticAirPocket
            logical          :: possibleAirPocket
        !%------------------------------------------------------------------
        !% cycle through the conduits to find air Pockets,
        do ii = 1,N_super_conduit
            !% pointers
            cIdx      => sc_link_Idx(ii,1:links_per_sc(ii))
            nElem     =  sum(link%I(cIdx,li_N_element))
            conAir    => link%YN(:,lYN_airPocketDetected)
            eIdx      => conduitElemMapsI(ii,1:nElem,cmi_elem_idx)
            fUp       => conduitElemMapsI(ii,1:nElem,cmi_elem_up_face)
            fDn       => conduitElemMapsI(ii,1:nElem,cmi_elem_dn_face)
            elemSur   => elemYN(:,eYN_isPSsurcharged)
            fBlocked  => faceYN(:,fYN_isAirflowBlocked)

            !% reset the possible air pocket detection logical
            possibleAirPocket = .false.
            !% reset the conduit air pocket detection logical
            conAir(cIdx) = .false.

            !% ------------------------------------------------------------------
            !% initial air pockets screening
            !% search for if any faces conduit elements that has blocked airflow
            if (any(elemSur(eIdx))) then
                possibleAirPocket = .true.
            end if

            !% if initially any airPocket is detected, map the location of that air Pocket
            if (possibleAirPocket) then

                if (airYN(ii,oneI,airYN_air_pocket_detected)) then
                    !% air pocket already detected
                    conAir(cIdx) = .true.
                    !% --- special conditions for upstream and downstream 
                    !%     airflow blockage
                    if (airI(ii,oneI,airI_type) == upReleaseAirpocket)  then
                        !% if the airflow of upstream release is blocked 
                        if (fBlocked(fUp(oneI))) then
                            airI(ii,oneI,airI_type)  = entrappedAirpocket
                        end if
                        
                    else if (airI(ii,oneI,airI_type) == dnReleaseAirpocket) then
                        !% if the airflow of downstream release is blocked
                        if (fBlocked(fDn(nElem))) then
                            airI(ii,oneI,airI_type)  = entrappedAirpocket
                        end if

                    else if (airI(ii,oneI,airI_type) == entrappedAirpocket) then
                        !% if the airflow of downstream release is blocked
                        if (.not. fBlocked(fDn(nElem))) then
                            airI(ii,oneI,airI_type)  = dnReleaseAirpocket
                        else if (fBlocked(fUp(oneI))) then
                            airI(ii,oneI,airI_type)  = upReleaseAirpocket
                        end if
                    end if

                else
                    !% reset the counter for number of Pockets
                    airPocketIdx = zeroI
                    !% reset the indexes
                    startIdx    = zeroI
                    endIdx      = zeroI
                    !% initialize the conduitElemMapsI at each detection for moving air pocket
                    conduitElemMapsI(ii,:,cmi_airpocket_idx)  = nullvalueI
                    conduitElemMapsI(ii,:,cmi_airpocket_type) = noAirPocket
                    !% initialize the index for the second do loop
                    jj = oneI
                    !% --- cycle through the conduit elements to 
                    !%     find entrapped air pockets
                    do while (jj <= nElem)
                        !% --- Find the airpocket based on mixed condition first
                        !%     the element is not surcharged and start counting 
                        !%     the non-surcharged elements
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
                            !%     consnsts of only one element 
                            if ((endIdx - startIdx <= zeroI))    then
                                airPocketIdx = airPocketIdx - oneI
                                cycle
                            end if

                            !% --- air pocket limit
                            !%     if the airpocket limit doesnot exceed the permissible
                            !%     number of airpockets per conduit, store the data
                            if ((airPocketIdx > zeroI) .and. (airPocketIdx <= oneI)) then
                                !% new airpocket
                                airYN(ii,oneI,airYN_new_air_pocket) =  .true.
                                !% save the airpocket detection logical
                                airYN(ii,oneI,airYN_air_pocket_detected) = .true.
                                !% arrayI_pos selects the right column positions of airpockets
                                !% save the integer air pocket data
                                airI(ii,oneI,airI_idx)        = oneI
                                airI(ii,oneI,airI_elem_start) = eIdx(startIdx)
                                airI(ii,oneI,airI_face_up)    = fUp(startIdx)
                                airI(ii,oneI,airI_elem_end)   = eIdx(endIdx)
                                airI(ii,oneI,airI_face_dn)    = fDn(endIdx)

                                !% --- set the type of the airpocket
                                !%     if the starting element is the first element in the conduit
                                !%     there will be an upstream release
                                if (startIdx == oneI) then
                                    airI(ii,oneI,airI_type) = upReleaseAirpocket
                                
                                !% else if the ending element is the last element in the conduit
                                !% there will be a downstream release
                                else if (endIdx == nElem) then
                                    airI(ii,oneI,airI_type) = dnReleaseAirpocket
                                
                                !% else the pocket is entrapped
                                else
                                    airI(ii,oneI,airI_type) = entrappedAirpocket
                                end if

                                !% save airpocket detection data in the conduitElemMapsI array
                                conduitElemMapsI(ii,startIdx:endIdx,cmi_airpocket_idx)  = oneI
                                conduitElemMapsI(ii,startIdx:endIdx,cmi_airpocket_type) = airI(ii,oneI,airI_type)

                            !% else if there is too many air pockets than permissible  
                            else if (airPocketIdx > oneI) then
                                print*, "The conduit has, ", airPocketIdx, " airpockets "
                                print*, "which is more than maximum permissible of 1 air pockets for the static method"
                                print*
                            end if
                        !% progress the counter for surcharge elements
                        else
                            jj = jj + oneI
                        end if
                    end do
                end if
            !% if there is not any possible air pockets, reset the values of the air arrays of that correspondig conduit
            else
                conAir(cIdx)  = .false.
                airYN(ii,:,:) = .false.
                airI(ii,:,:)  = nullvalueI
                !% reset the air entrapment type
                airI(ii,:,airI_type) = noAirPocket
                airR(ii,:,:)  = zeroR
                airR(ii,:,airR_absolute_head) = setting%AirTracking%AtmosphericPressureHead
                airR(ii,:,airR_density)       = setting%AirTracking%AirDensity
                conduitElemMapsI(ii,:,cmi_airpocket_idx)  = nullvalueI
                conduitElemMapsI(ii,:,cmi_airpocket_type) = noAirPocket
            end if
        end do

    end subroutine airpockets_detection_single
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
end module air_entrapment