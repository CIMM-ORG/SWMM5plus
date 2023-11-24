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
            integer             :: ii
        !%------------------------------------------------------------------

        call airpockets_detection ()

        call airpockets_calculation (istep)

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
            integer, pointer :: cIdx, max_airpockets
            logical, pointer :: conAir
        !%------------------------------------------------------------------
        !% cycle through the conduits to find air Pockets,
        do ii = 1,N_conduit
            !% pointers
            cIdx      => pConduitIdx(ii)
            conAir    => link%YN(cIdx,lYN_airPocketDetected)

            max_airpockets => setting%AirTracking%NumberOfAirpocketsAllowed

            !% only go through the conduit airpocket calculation if the
            !% conduit contains one
            if (conAir) then

                ! print*, '------------------------------------------------------'

                do jj = 1,max_airpockets 
    
                    !% initialize the airpocket real values
                    call airpocket_initialization (ii,jj, istep)
      
                    !% calculate the net flowrate though the airpockets
                    call airpocket_netflowrate (ii,jj)
    
                    !% calculate air outflow rate form an airpocket if any vent is present
                    call airpocket_air_mass_outflow (ii, jj)

                    !% calculate the air pressure head
                    call airpocket_pressure_head (ii, jj, istep)

                    !% time-march volume
                    call airpockets_volume_update (ii, jj, istep)

                    !% update the air density at the air pocket
                    call airpocket_air_density_update (ii, jj, istep)  

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
            integer          :: ii, jj, startIdx, endIdx, airPocketIdx
            integer, pointer :: cIdx, nElem, eIdx(:), fUp(:), fDn(:), max_airpockets
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

            max_airpockets => setting%AirTracking%NumberOfAirpocketsAllowed

            !% reset the possible air pocket detection logical
            possibleAirPocket = .false.
            !% reset the conduit air pocket detection logical
            conAir = .false.

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
                        if ((airPocketIdx > zeroI) .and. (airPocketIdx <= max_airpockets)) then
                            !% arrayI_pos selects the right column positions of airpockets
                            !% save the integer air pocket data
                            airI(ii,airPocketIdx,airI_idx)        = airPocketIdx
                            airI(ii,airPocketIdx,airI_elem_start) = eIdx(startIdx)
                            airI(ii,airPocketIdx,airI_face_up)    = fUp(startIdx)
                            airI(ii,airPocketIdx,airI_elem_end)   = eIdx(endIdx)
                            airI(ii,airPocketIdx,airI_face_dn)    = fDn(endIdx)

                            !% save the airpocket detection logical
                            airYN(ii,airPocketIdx,airYN_air_pocket_detected) = .true.

                            !% save the airpocket detection at the conduit
                            conAir = .true.

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

            !% if there is not any possible air pockets, reset the values of the air arrays of that correspondig conduit
            else
                !% reset the logicals
                conAir = .false.
                airYN(ii,:,:) = .false.
                airI(ii,:,:)  = nullvalueI
                !% reset the air entrapment type
                airI(ii,:,airI_type) = noAirPocket
                airR(ii,:,:)  = zeroR
                airR(ii,:,airR_absolute_head) = setting%AirTracking%AtmosphericPressureHead
            end if

            !% debug printing
            ! if (any(airYN(ii,:,airYN_air_pocket_detected))) then
                ! print*, 'Entrapped air at conduit  ',link%names(ii)%str
                ! print*, 'airPocketIdx', airPocketIdx
                ! print*, 'idx 1     = ', airI(ii,1,airI_idx),         'idx 2     = ', airI(ii,2,airI_idx),        'idx 3     = ',airI(ii,3,airI_idx) 
                ! print*, 'type 1    = ', airI(ii,1,airI_type),        'type 2    = ', airI(ii,2,airI_type),       'type 3    = ', airI(ii,3,airI_type)
                ! print*, 'elem up 1 = ', airI(ii,1,airI_elem_start),  'elem up 2 = ', airI(ii,2,airI_elem_start), 'elem up 3 = ',airI(ii,3,airI_elem_start) 
                ! print*, 'elem dn 1 = ', airI(ii,1,airI_elem_end),    'elem dn 2 = ', airI(ii,2,airI_elem_end),   'elem dn 3 = ',airI(ii,3,airI_elem_end) 
                ! print* 
            ! end if

        end do

    end subroutine airpockets_detection
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine airpocket_initialization (cIdx,aIdx,istep)
        !%------------------------------------------------------------------
        !% Description:
        !%   newly detected airpocket: initialize the inflow, outflow,  
        !%   volume and heads
        !%   already detected airpockets: initialize the inflow and outflow 
        !%   no air pockets: zero the values 
        !%------------------------------------------------------------------
            integer, intent (in) :: cIdx, aIdx, istep
            integer, pointer     :: elemStartIdx, elemEndIdx, faceUp, faceDn
            real(8), pointer     :: airVolume, airVolume_N0, inflow, outflow 
            real(8), pointer     :: airMass, airMass_N0, airDensity, airDensity_N0
            real(8), pointer     :: absHead, absHead_N0, gaugeHead, atmHead, rho_a
            real(8), pointer     :: elemVol(:), fullVol(:), faceFlow(:)
            logical, pointer     :: isAirPocket, newAirPocket
        !%------------------------------------------------------------------
        !% Aliases
            elemStartIdx  => airI(cIdx,aIdx,airI_elem_start)
            elemEndIdx    => airI(cIdx,aIdx,airI_elem_end)
            faceUp        => airI(cIdx,aIdx,airI_face_up)
            faceDn        => airI(cIdx,aIdx,airI_face_dn)
            airVolume     => airR(cIdx,aIdx,airR_volume)
            airVolume_N0  => airR(cIdx,aIdx,airR_volume_N0)
            airDensity    => airR(cIdx,aIdx,airR_density)
            airDensity_N0 => airR(cIdx,aIdx,airR_density_N0)
            airMass       => airR(cIdx,aIdx,airR_mass)
            airMass_N0    => airR(cIdx,aIdx,airR_mass_N0)
            inflow        => airR(cIdx,aIdx,airR_inflow)
            outflow       => airR(cIdx,aIdx,airR_outflow)
            absHead       => airR(cIdx,aIdx,airR_absolute_head)
            absHead_N0    => airR(cIdx,aIdx,airR_absolute_head_N0)
            gaugeHead     => airR(cIdx,aIdx,airR_gauge_head)
            isAirPocket   => airYN(cIdx,aIdx,airYN_air_pocket_detected)
            newAirPocket  => airYN(cIdx,aIdx,airYN_new_air_pocket)

        !% other aliases
            elemVol   => elemR(:,er_Volume)
            fullVol   => elemR(:,er_FullVolume)
            faceFlow  => faceR(:,fr_Flowrate)
            atmHead   => setting%AirTracking%AtmosphericPressureHead
            rho_a     => setting%AirTracking%AirDensity

        !% if a airpocket is detected, initialize the air pocket
        if (isAirPocket) then
            !% in this algorithm, aiVolume = 0 means the air pocket is newly detected
            if (airVolume == zeroR) then
                newAirPocket = .true.
                gaugeHead     = zeroR
                airDensity    = rho_a
                airDensity_N0 = rho_a
                absHead       = atmHead
                absHead_N0    = atmHead 
                !% calculate the air pocket volume from empty space
                airVolume  = max(sum(fullVol(elemStartIdx:elemEndIdx) - elemVol(elemStartIdx:elemEndIdx)), zeroR)
                !% for a new pocket formulation save the older values
                airVolume_N0 = airVolume
                airMass      = airDensity * airVolume
                airMass_N0   = airMass
            end if
            
            !% copy over the flow data
            inflow     = faceFlow(faceUp)
            outflow    = faceFlow(faceDn)
            
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
    subroutine airpocket_netflowrate (cIdx, aIdx)
        !%------------------------------------------------------------------
        !% Description:
        !%   calculate the net flowrate through an airpocket 
        !%------------------------------------------------------------------
            integer, intent (in) :: cIdx, aIdx
            integer, pointer     :: fUp(:), fDn(:), elemStartIdx, elemEndIdx
            real(8), pointer     :: dvdt, inflow, outflow, faceFlow(:)
            logical, pointer     :: isAirPocket
        !%------------------------------------------------------------------
        !% Aliases
            dvdt        => airR(cIdx,aIdx,airR_dvdt)
            inflow      => airR(cIdx,aIdx,airR_inflow)
            outflow     => airR(cIdx,aIdx,airR_outflow)
            isAirPocket => airYN(cIdx,aIdx,airYN_air_pocket_detected)


            elemStartIdx => airI(cIdx,aIdx,airI_elem_start)
            elemEndIdx   => airI(cIdx,aIdx,airI_elem_end)
            faceFlow  => faceR(:,fr_Flowrate)
            fUp       => elemI(:,ei_Mface_uL)
            fDn       => elemI(:,ei_Mface_dL)
        
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
    subroutine airpocket_air_mass_outflow (cIdx, aIdx)
        !%------------------------------------------------------------------
        !% Description:
        !%   calculates the air outflow rate for airpocket with release
        !%------------------------------------------------------------------
            integer, intent (in) :: cIdx, aIdx
            real(8)              :: ExpansionFactor, ratio, areaOpening
            integer, pointer     :: pocketType, elemStartIdx, elemEndIdx, faceUp, faceDn
            real(8), pointer     :: massOutflow, absHead, dt, kappa, atmHead, airDensity
            real(8), pointer     :: dishCoeff, rho_w, rho_a, grav, minVentArea
            logical, pointer     :: isAirPocket, isVacuumed
        !%------------------------------------------------------------------
        !% Aliases
            pocketType   => airI(cIdx,aIdx,airI_type)
            elemStartIdx => airI(cIdx,aIdx,airI_elem_start)
            elemEndIdx   => airI(cIdx,aIdx,airI_elem_end)
            faceUp       => airI(cIdx,aIdx,airI_face_up)
            faceDn       => airI(cIdx,aIdx,airI_face_dn)
            absHead      => airR(cIdx,aIdx,airR_absolute_head)
            massOutflow  => airR(cIdx,aIdx,air_mass_flowrate)
            airDensity   => airR(cIdx,aIdx,airR_density)
            isAirPocket  => airYN(cIdx,aIdx,airYN_air_pocket_detected)
            isVacuumed   => airYN(cIdx,aIdx,airYN_air_pocket_vacuumed)
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

                    !% copy the area opening from face
                    ! areaOpening = max(elemR(elemStartIdx,er_FullArea) - faceR(faceUp,fr_Area_d), minVentArea)
                    areaOpening = minVentArea

                case (dnReleaseAirpocket)

                    !% copy the area opening from face
                    ! areaOpening = max(elemR(elemEndIdx,er_FullArea) - faceR(faceDn,fr_Area_u), minVentArea)
                    areaOpening = minVentArea

                case (noAirPocket)

                    areaOpening = zeroR

                case default

            end select

            !% calculate the ratio between atmospheric vs. absolute pressue head at conduit
            ratio = atmHead / absHead

            !% choked orifice airflow
            if (ratio > zeroR .and. ratio < 0.53 ) then
                if (airDensity > zeroR) then
                    massOutflow = dishCoeff * areaOpening * airDensity * sqrt(grav * (rho_w / airDensity) &
                                * absHead) * sqrt(kappa * (twoR / (kappa + oneR)) ** ((kappa + oneR) / (kappa - oneR)))
                else
                    massOutflow = zeroR
                end if
            !% normal orifice airflow
            else if (ratio >= 0.53 .and. ratio < 1.00) then

                !% find the expansion factor
                ExpansionFactor = sqrt((kappa / (kappa - oneR)) * (ratio ** (twoR / kappa)) &
                                    * ((oneR - ratio ** ((kappa - oneR) / kappa)) / (oneR - ratio)))

                if (airDensity > zeroR) then
                    

                    !% find the airflow
                    massOutflow = - dishCoeff * areaOpening * ExpansionFactor * airDensity &
                                * sqrt(twoR * grav * (rho_w / airDensity) * abs(absHead - atmHead))
                else
                    ! massOutflow = dishCoeff * areaOpening * ExpansionFactor * rho_a &
                    !             * sqrt(twoR * grav * (rho_w / rho_a) * abs(absHead - atmHead))

                    massOutflow = zeroR
                end if

            !% airflow into the airpocket
            else if (ratio > 1.00) then

                !% find the expansion factor
                ExpansionFactor = sqrt((kappa / (kappa - oneR)) * (ratio ** (twoR / kappa)) &
                                    * ((oneR - ratio ** ((kappa - oneR) / kappa)) / (oneR - ratio)))

                !% find the airflow
                massOutflow = dishCoeff * areaOpening * ExpansionFactor * airDensity &
                        * sqrt(twoR * grav * (rho_w / rho_a) * abs(atmHead - absHead))   
            else

                massOutflow = zeroR
            end if

            ! print*, '... ... ... ... ... ... ... ... ... ...'

            ! print*, 'airpocket_air_mass_outflow'
            ! print*, reverseKey(pocketType), ' at conduit  ',link%names(pConduitIdx(cIdx))%str
            ! print*, 'time   = ', setting%Time%Now
            ! print*, 'a idx  = ', aIdx
            ! print*, 'rato   = ', ratio
            ! print*, 'rho_a  = ', airDensity
            ! print*, 'EFctor = ', ExpansionFactor
            ! print*, 'masflo = ', massOutflow
            ! print*

        end if


    end subroutine airpocket_air_mass_outflow
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine airpocket_pressure_head (cIdx, aIdx, istep)
        !%------------------------------------------------------------------
        !% Description:
        !%   Solves the airpocket pressure head
        !%------------------------------------------------------------------
            integer, intent (in) :: cIdx, aIdx, istep
            real                 :: alpha
            integer, pointer     :: pocketType
            real(8), pointer     :: airVolume, airMass, dvdt, absHead, gaugeHead
            real(8), pointer     :: absHead_N0, dHdt, massOutflow
            real(8), pointer     :: dt, kappa, atmHead, crk(:)
            logical, pointer     :: isAirPocket
        !%------------------------------------------------------------------
        !% Aliases
            pocketType  => airI(cIdx,aIdx,airI_type)
            airVolume   => airR(cIdx,aIdx,airR_volume)
            airMass     => airR(cIdx,aIdx,airR_mass)
            absHead     => airR(cIdx,aIdx,airR_absolute_head)
            absHead_N0  => airR(cIdx,aIdx,airR_absolute_head_N0)
            gaugeHead   => airR(cIdx,aIdx,airR_gauge_head)
            dvdt        => airR(cIdx,aIdx,airR_dvdt)
            massOutflow => airR(cIdx,aIdx,air_mass_flowrate)
            dHdt        => airR(cIdx,aIdx,airR_temp01)
            isAirPocket => airYN(cIdx,aIdx,airYN_air_pocket_detected)
            crk         => setting%Solver%crk2
            dt          => setting%Time%Hydraulics%Dt
            kappa       => setting%AirTracking%PolytropicExponent
            atmHead     => setting%AirTracking%AtmosphericPressureHead
        
        !% calculate the new airpocket pressure head
        if (isAirPocket) then
            !% find the dvdt from the inflow and outflows
            ! dHdt = - kappa * (absHead/airVolume) * dvdt + kappa * (absHead/airMass) * massOutflow
            ! dHdt = - kappa * (absHead/airVolume) * dvdt

            if (airMass > zeroR) then
                ! dHdt  = dHdt + kappa * (absHead/airMass) * massOutflow

                dHdt = - kappa * (absHead/airVolume) * dvdt + kappa * (absHead/airMass) * massOutflow
            else
                dHdt = zeroR
            end if
            
            absHead   = absHead_N0 + crk(istep) * dHdt * dt 
            absHead   = max(absHead,atmHead)
            gaugeHead = absHead - atmHead


            ! print*, 'airpocket_pressure_head'
            ! print*, 'airVol = ', airVolume
            ! print*, 'airMas = ', airMass
            ! print*, 'masflo = ', massOutflow
            ! print*, 'dHdt   = ', dHdt
            ! print*, 'dvdt   = ', dvdt
            ! print*, 'aHead  = ', absHead
            ! print*, 'gHead  = ', gaugeHead
            ! print*

            !% save the new air absolute head and push down
            !% after the end of an RK step
            if (istep == twoI) then
                absHead_N0 = absHead
            end if

        end if
    
    end subroutine airpocket_pressure_head
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine airpocket_air_density_update (cIdx, aIdx, istep)
        !%------------------------------------------------------------------
        !% Description:
        !%   Update the air density at the airpocket
        !%------------------------------------------------------------------
            integer, intent (in) :: cIdx, aIdx, istep
            integer, pointer     :: pocketType
            real(8), pointer     :: airVolume, massOutflow, airMass, airMass_N0
            real(8), pointer     :: airDensity_N0, dvdt, airDensity, dt, crk(:)
            logical, pointer     :: isAirPocket, isVacuumed
            real                 :: dRhodT
        !%------------------------------------------------------------------
        !% Aliases
            pocketType    => airI(cIdx,aIdx,airI_type)
            airVolume     => airR(cIdx,aIdx,airR_volume)
            dvdt          => airR(cIdx,aIdx,airR_dvdt)
            airMass       => airR(cIdx,aIdx,airR_mass)
            airMass_N0    => airR(cIdx,aIdx,airR_mass_N0)
            airDensity    => airR(cIdx,aIdx,airR_density)
            airDensity_N0 => airR(cIdx,aIdx,airR_density_N0)
            massOutflow   => airR(cIdx,aIdx,air_mass_flowrate)
            isAirPocket   => airYN(cIdx,aIdx,airYN_air_pocket_detected)
            isVacuumed    => airYN(cIdx,aIdx,airYN_air_pocket_vacuumed)
            dt            => setting%Time%Hydraulics%Dt
            crk           => setting%Solver%crk2

        !% update the new density based on mass outflow
        if (isAirPocket) then

            !% timemarch air mass through mass flowrate
            airMass = airMass_N0 + dt * crk(istep) * massOutflow

            airMass = max(airMass,zeroR)


            if (airMass > zeroR) then
                isVacuumed = .false.
                airDensity = airMass / airVolume
            else 
                airMass = zeroR
                isVacuumed = .true.
                airDensity =  zeroR
            end if

            ! print*, 'airpocket_air_density_update'
            ! print*, 'airVol = ', airVolume
            ! print*, 'airMas = ', airMass
            ! print*, 'rho_a  = ', airDensity
            ! print*

            !% save the new mass and push down to old mass
            !% after the end of an RK step
            if (istep == twoI) then
                airDensity_N0 = airDensity 
                airMass_N0    = airMass
            end if

        end if


    end subroutine airpocket_air_density_update
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine airpockets_volume_update (cIdx, aIdx, istep)
        !%------------------------------------------------------------------
        !% Description:
        !%   Timemarch airpocket volume
        !%------------------------------------------------------------------
            integer, intent (in) :: cIdx, aIdx, istep
            integer, pointer     :: pocketType
            real(8)              :: tempVol
            real(8), pointer     :: airVolume, dvdt, airVolume_N0, dt, crk(:)
            logical, pointer     :: isAirPocket
        !%------------------------------------------------------------------
        !% Aliases
            pocketType   => airI(cIdx,aIdx,airI_type)
            airVolume    => airR(cIdx,aIdx,airR_volume)
            dvdt         => airR(cIdx,aIdx,airR_dvdt)
            airVolume_N0 => airR(cIdx,aIdx,airR_volume_N0)
            isAirPocket  => airYN(cIdx,aIdx,airYN_air_pocket_detected)
            dt           => setting%Time%Hydraulics%Dt
            crk          => setting%Solver%crk2
        
        !% update the new volume based on inflow and outflow outflow
        if (isAirPocket) then

            airVolume = airVolume_N0 + crk(istep) * dt * dvdt

            !% save the new volume and push down to old mass
            !% after the end of an RK step
            if (istep == twoI) then
                airVolume_N0 = airVolume
            end if

        end if

    end subroutine airpockets_volume_update
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine add_airpocket_heads_to_elem (cIdx, aIdx)
        !%------------------------------------------------------------------
        !% Description:
        !%  add the airpocket heads to corresponding elements
        !%------------------------------------------------------------------
        integer, intent (in) :: cIdx, aIdx
        integer, pointer     :: elemStartIdx, elemEndIdx
        real(8), pointer     :: gaugeHead, elemHead(:)
        logical, pointer     :: isAirPocket
        !%------------------------------------------------------------------
        !% Aliases
            elemStartIdx => airI(cIdx,aIdx,airI_elem_start)
            elemEndIdx   => airI(cIdx,aIdx,airI_elem_end)
            gaugeHead    => airR(cIdx,aIdx,airR_gauge_head)
            isAirPocket  => airYN(cIdx,aIdx,airYN_air_pocket_detected)
            elemHead     => elemR(:,er_Head)
        
        !% add the gauge heads to elem heads
        if (isAirPocket) then
            elemHead(elemStartIdx:elemEndIdx) = elemHead(elemStartIdx:elemEndIdx) + gaugeHead
        end if


    end subroutine add_airpocket_heads_to_elem
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    
end module air_entrapment