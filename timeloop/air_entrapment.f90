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
            ! call airpockets_detection_single ()
            
            !% only three airpockets per superconduit is allowed
            call airpockets_detection_static ()
        else
            !% find the airpockets in the network
            call airpockets_detection ()
        end if

        !% find the pressure head due to airpocket
        call airpockets_calculation (istep)

        if (setting%AirTracking%AirVentThroughJM) then
            !% --- compute the pressure head in a JM junction
            call airpockets_junction (istep)
            !% --- add the pressure head to CC elements adjacent to junctions
            !%     and update airR arrays
            call airpockets_update_air_pressure_from_JM (istep)
        end if

        !% pack the faces those are needed to be interpolated
        call airpocket_face_pack ()
        
        !% interpolate the faces again after air calculation
        !% to update the new heads to the faces (only head interp)
        call face_interpolation(fp_Airpockets, .false., .true., .false., .true., .true.)

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
        faceYN(:,fYN_isAirPressurized)                =  .false.

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

                    !% copy the mass outflow to JM
                    call airpockets_junction_air_maps (ii, jj)

                    !% update the air density at the air pocket
                    call airpocket_air_density_update (ii, jj, istep)  

                    !% calculate the air pressure head
                    call airpocket_pressure_head (ii, jj, istep)

                    !% add the heads back to the CC elements
                    call add_airpocket_heads_to_CC_elem (ii, jj)

                end do
            
            end if
                
        end do
        

    end subroutine airpockets_calculation
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine airpockets_junction (istep)
        !%------------------------------------------------------------------
        !% Description:
        !% air pocket calculation for a junction
        !% Requires the mass inflow rate from a conduit to be already 
        !% provided
        !%------------------------------------------------------------------
            integer, intent(in) :: istep
            integer, pointer    :: Npack, thisJM(:)
            integer             :: mm, JMidx
            real(8), pointer    :: kappa, dt, crk, theta, atmHead, rho_a
            real(8)             :: alpha, beta, dvdt, coef
            logical             :: JunctionAirPocket = .false.
        !%------------------------------------------------------------------
        !% Preliminaries
            Npack => npack_elemP(ep_JM)   
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
            thisJM  => elemP(1:Npack,ep_JM)
            kappa   => setting%AirTracking%PolytropicExponent
            dt      => setting%Time%Hydraulics%Dt
            crk     => setting%Solver%crk2(istep)
            theta   => setting%AirTracking%theta
            atmHead => setting%AirTracking%AtmosphericPressureHead
            rho_a   => setting%AirTracking%AirDensity
        !%------------------------------------------------------------------
        
        !% --- cycle through the JM junctions
        do mm = 1,Npack
            JMidx = thisJM(mm)

            !% --- find if this junction has airpocket
            JunctionAirPocket = airpocket_detection_JM (JMidx)
            
            !% if the junction doesnot contain air pocket, cycle
            if (.not. JunctionAirPocket) then
                elemR(JMidx,er_Pressurized_Air)          = zeroR
                elemR(JMidx,er_Air_Pressure_Head)        = zeroR
                elemSR(JMidx,esr_JM_Air_Volume)          = zeroR
                elemSR(JMidx,esr_JM_Air_Volume_N0)       = zeroR
                elemSR(JMidx,esr_JM_Air_HeadAbsolute)    = atmHead
                elemSR(JMidx,esr_JM_Air_HeadAbsolute_N0) = atmHead
                elemSR(JMidx,esr_JM_Air_Density)         = rho_a
                elemYN(JMidx,eYN_hasAirPocket)           = .false.
                cycle
            end if

            !% --- determine if the airpocket is new to initialize the airpocket
            if (elemSR(JMidx,esr_JM_Air_Volume) == zeroR) then
                !% --- compute the air volume 
                elemSR(JMidx,esr_JM_Air_Volume)    = airpockets_JM_airvolume(JMidx)
                elemSR(JMidx,esr_JM_Air_Volume_N0) = elemSR(JMidx,esr_JM_Air_Volume)
                elemSR(JMidx,esr_JM_Air_Mass)      = rho_a * elemSR(JMidx,esr_JM_Air_Volume)
                elemSR(JMidx,esr_JM_Air_Mass_N0)   = elemSR(JMidx,esr_JM_Air_Mass)
                elemSR(JMidx,esr_JM_Air_Density)   = rho_a
                elemSR(JMidx,esr_JM_Air_HeadAbsolute)    = atmHead
                elemSR(JMidx,esr_JM_Air_HeadAbsolute_N0) = atmHead
            else
                !% --- else only calculate the air volume
                elemSR(JMidx,esr_JM_Air_Volume)    = airpockets_JM_airvolume(JMidx)
            end if

            
            if (elemSR(JMidx,esr_JM_Air_Volume) .le. zeroR) then 
                !% --- air pocket has collapsed and JMidx is full of water (incipient surcharging)
                !%     zero out pocket and cycle to next JMidx
                elemSR(JMidx,esr_JM_Air_Volume)    = zeroR
            end if

            !% --- compute the mass outflow to atmosphere
            elemSR(JMidx,esr_JM_Air_MassOutflowRate) = airpockets_JM_outflow (JMidx) 

            !% --- compute the new air mass in the JMidx
            !%     NOTE THAT POSITIVE VALUES OF INFLOWS ARE INFLOWS AND POSITIVE OUTFLOWS ARE OUTFLOWS
            elemSR(JMidx,esr_JM_Air_Mass) = max(elemSR(JMidx,esr_JM_Air_Mass_N0) &
                + crk * dt * ( - elemSR(JMidx,esr_JM_Air_MassOutflowRate)), zeroR)

            !% update the air density
            elemSR(JMidx,esr_JM_Air_Density) = elemSR(JMidx,esr_JM_Air_Mass) / elemSR(JMidx,esr_JM_Air_Volume)

            if (elemSR(JMidx,esr_JM_Air_Mass) .le. zeroR) then 
                !% --- no air pressurization
                elemSR(JMidx,esr_JM_Air_Density)  = zeroR
                elemSR(JMidx,esr_JM_Air_Mass)     = zeroR
            end if

            !% --- compute the air volume rate of change (opposite of water volume change)
            !%     Negative is decreasing volume
            dvdt  = airpockets_JM_dvdt (JMidx, istep)

            !% --- beta calculation
            if (elemSR(JMidx,esr_JM_Air_Mass) > zeroR)  then                             
                !% --- compute the beta term, which should be negative for a net outflow
                beta = ( kappa / elemSR(JMidx,esr_JM_Air_Mass) ) * (- elemSR(JMidx,esr_JM_Air_MassOutflowRate))
            else 
                beta = zeroR
            end if

            !% --- alpha calculation
            if (elemSR(JMidx,esr_JM_Air_Volume) > zeroR) then
                !% --- compute the alpha term
                alpha = ( kappa / elemSR(JMidx,esr_JM_Air_Volume)) * dvdt 
            else
                alpha = zeroR
            end if

            if (elemSR(JMidx,esr_JM_Air_Volume) > zeroR .and. elemSR(JMidx,esr_JM_Air_Mass) > zeroR) then
                !% --- compute coef of theta method
                coef = (beta - alpha) * dt * crk

                !% --- Theta method to update the JM Head
                elemSR(JMidx,esr_JM_Air_HeadAbsolute) = elemSR(JMidx,esr_JM_Air_HeadAbsolute)   &
                                       * ((oneR + coef * (oneR - theta))  / (oneR - coef * theta))
                !% --- find the gauge head
                elemSR(JMidx,esr_JM_Air_HeadGauge)     =  elemSR(JMidx,esr_JM_Air_HeadAbsolute) - atmHead
            else
                elemSR(JMidx,esr_JM_Air_HeadAbsolute)  = atmHead
                elemSR(JMidx,esr_JM_Air_HeadGauge)     = zeroR
            end if

            !% --- check for head below atmospheric, which indicates the air pocket has emptied
            if (elemSR(JMidx,esr_JM_Air_HeadGauge) .le. zeroR) then   
                !% --- zero out the air head and mass
                elemSR(JMidx,esr_JM_Air_HeadGauge) = zeroR
                elemR(JMidx,er_Pressurized_Air)    = zeroR
                elemR(JMidx,er_Air_Pressure_Head)  = zeroR
                elemYN(JMidx,eYN_hasAirPocket)     = .false.
            else 
                !% --- update the air head (gauge)
                elemR(JMidx,er_Head) = elemR(JMidx,er_Head) + elemSR(JMidx,esr_JM_Air_HeadGauge)
                elemR(JMidx,er_Pressurized_Air)   = oneR
                elemR(JMidx,er_Air_Pressure_Head) = elemSR(JMidx,esr_JM_Air_HeadGauge) 
                elemYN(JMidx,eYN_hasAirPocket)    = .true.
            end if
        end do

    end subroutine airpockets_junction    
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
                            print*, "The conduit", ii," has, ", airPocketIdx, " airpockets "
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
                airYN(ii,:,airYN_air_pocket_detected)  = .false.
                airYN(ii,:,airYN_air_pocket_collapsed) = .true.
                airYN(ii,:,airYN_new_air_pocket)       = .false.
                airI(ii,:,airI_idx)                    = nullvalueI
                airI(ii,:,airI_elem_start)             = nullvalueI
                airI(ii,:,airI_elem_end)               = nullvalueI
                airI(ii,:,airI_face_up)                = nullvalueI
                airI(ii,:,airI_face_dn)                = nullvalueI
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
            real(8), pointer     :: airVolume, inflow, outflow, airDensity, airMass_N0
            real(8), pointer     :: airMass, absHead, gaugeHead, atmHead, rho_a, absHead_N0
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
            airMass_N0    => airR(sc_Idx,aIdx,airR_mass_N0)
            inflow        => airR(sc_Idx,aIdx,airR_inflow)
            outflow       => airR(sc_Idx,aIdx,airR_outflow)
            absHead       => airR(sc_Idx,aIdx,airR_absolute_head)
            absHead_N0    => airR(sc_Idx,aIdx,airR_absolute_head_N0)
            gaugeHead     => airR(sc_Idx,aIdx,airR_gauge_head)
            isAirPocket   => airYN(sc_Idx,aIdx,airYN_air_pocket_detected)
            newAirPocket  => airYN(sc_Idx,aIdx,airYN_new_air_pocket)

        !% other aliases
            elemVol   => elemR(:,er_Volume)
            fullVol   => elemR(:,er_FullVolume)
            faceFlow  => faceR(:,fr_Flowrate)
            atmHead   => setting%AirTracking%AtmosphericPressureHead
            rho_a     => setting%AirTracking%AirDensity

            !% flowrate pointer based on the RK step
            if (istep == oneI) then
                faceFlow  => faceR(:,fr_Flowrate)
            else
                faceFlow  => faceR(:,fr_Flowrate_Conservative)
            end if

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
                absHead_N0    = atmHead
                !% calculate the air pocket volume from empty space
                airVolume  = max(sum(fullVol(pElem) - elemVol(pElem)), zeroR)
                !% for a new pocket formulation save the older values
                airMass      = airDensity * airVolume
                airMass_N0   = airMass
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
            airMass_N0  = zeroR
            inflow      = zeroR
            outflow     = zeroR
            gaugeHead   = zeroR
            airDensity  = rho_a
            absHead     = atmHead
            absHead_N0  = atmHead
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
            logical, pointer     :: isAirPocket, JunctionVent, UseMinVentArea, JunctionAirPocket
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
            JunctionAirPocket => airYN(sc_Idx,aIdx,airYN_junction_air_pocket)

            dt           => setting%Time%Hydraulics%Dt
            grav         => setting%constant%gravity
            kappa        => setting%AirTracking%PolytropicExponent
            atmHead      => setting%AirTracking%AtmosphericPressureHead
            dishCoeff    => setting%AirTracking%AirDischargeCoefficient
            rho_w        => setting%AirTracking%WaterDensity
            rho_a        => setting%AirTracking%AirDensity
            minVentArea  => setting%AirTracking%MinimumVentArea
            JunctionVent => setting%AirTracking%AirVentThroughJM
            UseMinVentArea => setting%AirTracking%UseMinVentArea

        !% calculate the new air reseale rate from the opening
        !% calculate the new airpocket pressure head
        if (isAirPocket) then
            select case (pocketType)

                case (entrappedAirpocket)

                    areaOpening = zeroR

                case (upReleaseAirpocket)

                    if (JunctionVent) then
                        areaOpening =  zeroR
                        !% this airpocket vents into a junction. So the air discharge will be calculated at the junction
                        JunctionAirPocket = .true.
                    else
                        if (UseMinVentArea) then
                            areaOpening = minVentArea
                        else
                            areaOpening = max(elemR(elemStartIdx,er_FullArea) - elemR(elemStartIdx,er_Area), zeroR)
                        end if
                    end if
                    
                case (dnReleaseAirpocket)      

                    if (JunctionVent) then
                        areaOpening =  zeroR
                        !% this airpocket vents into a junction. So the air discharge will be calculated at the junction
                        JunctionAirPocket = .true.
                    else
                        if (UseMinVentArea) then
                            areaOpening = minVentArea
                        else
                            areaOpening = max(elemR(elemEndIdx,er_FullArea) - elemR(elemEndIdx,er_Area), zeroR)
                        end if
                    end if

                case (noAirPocket)

                    areaOpening = zeroR

                case default
                    print*, 'ERROR: this must not reach'
                    stop 8413354

            end select

            !% only calculate the mass outflows for the non-junction airpockets
            if (.not. JunctionAirPocket) then
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
        end if

    end subroutine airpocket_air_mass_outflow
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine airpockets_junction_air_maps (sc_Idx, aIdx)
        !%------------------------------------------------------------------
        !% Description:
        !%   calculates the air exchange rate between vlink and JM
        !%------------------------------------------------------------------
            integer, intent (in) :: sc_Idx, aIdx
            integer, pointer     :: pocketType, upJMidx, dnJMidx, upJBidx, dnJBidx
            logical, pointer     :: ventedUpJM, ventedDnJM, isAirPocket, JunctionAirPocket
            real(8), pointer     :: massOutflow
        !%------------------------------------------------------------------
        !% Aliases
            pocketType   => airI(sc_Idx,aIdx,airI_type)
            upJMidx      => airI(sc_Idx,aIdx,airI_Up_JM_idx)
            upJBidx      => airI(sc_Idx,aIdx,airI_Up_JB_idx)
            dnJMidx      => airI(sc_Idx,aIdx,airI_Dn_JM_idx)
            dnJBidx      => airI(sc_Idx,aIdx,airI_Dn_JB_idx)
            ventedUpJM   => airYN(sc_Idx,aIdx,airYN_air_vented_through_UpJM)
            ventedDnJM   => airYN(sc_Idx,aIdx,airYN_air_vented_through_DnJM)
            massOutflow  => airR(sc_Idx,aIdx,airR_mass_flowrate)
            isAirPocket  => airYN(sc_Idx,aIdx,airYN_air_pocket_detected)
            JunctionAirPocket => airYN(sc_Idx,aIdx,airYN_junction_air_pocket)

            if (isAirPocket) then
                select case (pocketType)

                    case (noAirPocket)
                        !% do nothing
                    case (entrappedAirpocket)
                        !% do nothing
                    case (upReleaseAirpocket)
                        if (ventedUpJM .and. JunctionAirPocket) then
                            ! elemSR(upJMidx,esr_JM_Air_MassInflowRate) = - massOutflow + elemSR(upJMidx,esr_JM_Air_MassInflowRate)
                            elemSI(upJBidx,esi_JB_vLink_Connection) = sc_Idx
                            elemSI(upJBidx,esi_JB_air_pocket_index) = aIdx
                            elemYN(upJBidx,eYN_hasAirPocket) = .true.
                        end if
                    case (dnReleaseAirpocket)
                        if (ventedDnJM .and. JunctionAirPocket) then
                            ! elemSR(dnJMidx,esr_JM_Air_MassInflowRate) = - massOutflow + elemSR(upJMidx,esr_JM_Air_MassInflowRate)
                            elemSI(dnJBidx,esi_JB_vLink_Connection) = sc_Idx
                            elemSI(dnJBidx,esi_JB_air_pocket_index) = aIdx
                            elemYN(dnJBidx,eYN_hasAirPocket) = .true.
                        end if

                    case default
                        print*, 'ERROR: this must not reach'
                        stop 8413354
                end select

            end if

    end subroutine airpockets_junction_air_maps
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
            real(8), pointer     :: dt, kappa, atmHead, theta, crk
            logical, pointer     :: isAirPocket, JunctionAirPocket
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
            JunctionAirPocket => airYN(sc_Idx,aIdx,airYN_junction_air_pocket)
            crk         => setting%Solver%crk2(istep)
            dt          => setting%Time%Hydraulics%Dt
            kappa       => setting%AirTracking%PolytropicExponent
            atmHead     => setting%AirTracking%AtmosphericPressureHead
            theta       => setting%AirTracking%theta
        
        !% calculate the new airpocket pressure head
        if (isAirPocket .and. .not. JunctionAirPocket) then

            if (airVolume > zeroR .and. airMass > zeroR) then
                !% find the alpha and beta
                alpha = (kappa / airVolume) * dvdt
                beta  = (kappa / airMass)   * dmdt

                !% find the absolute head  
                ! absHead = absHead_N0 * (oneR + dt * crk * (oneR - theta) * (- alpha + beta)) / (oneR - dt * crk * theta * (- alpha + beta))

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
            real(8), pointer     :: airDensity_N0, dvdt, airDensity, dt, rho_a, crk
            logical, pointer     :: isAirPocket, colAirpocket, JunctionAirPocket
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
            JunctionAirPocket => airYN(sc_Idx,aIdx,airYN_junction_air_pocket)
            colAirpocket  => airYN(sc_Idx,aIdx,airYN_air_pocket_collapsed)
            dt            => setting%Time%Hydraulics%Dt
            crk           => setting%Solver%crk2(istep)
            rho_a         => setting%AirTracking%AirDensity

        !% update the new density based on mass outflow
        if (isAirPocket .and. .not. JunctionAirPocket) then

            !% timemarch air mass through mass flowrate
            ! airMass = airMass_N0 + dt * crk * massOutflow

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
    ! subroutine airpockets_volume_update (sc_Idx, aIdx, istep)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !%   Timemarch airpocket volume
    !     !%------------------------------------------------------------------
    !         integer, intent (in) :: sc_Idx, aIdx, istep
    !         real(8)              :: tempVol
    !         real(8), pointer     :: airVolume, dvdt, airVolume_N0, dt, crk(:)
    !         logical, pointer     :: isAirPocket
    !     !%------------------------------------------------------------------
    !     !% Aliases
    !         airVolume    => airR(sc_Idx,aIdx,airR_volume)
    !         dvdt         => airR(sc_Idx,aIdx,airR_dvdt)
    !         airVolume_N0 => airR(sc_Idx,aIdx,airR_volume_N0)
    !         isAirPocket  => airYN(sc_Idx,aIdx,airYN_air_pocket_detected)
    !         dt           => setting%Time%Hydraulics%Dt
    !         crk          => setting%Solver%crk2
        
    !     !% update the new volume based on inflow and outflow outflow
    !     if (isAirPocket) then

    !         airVolume = airVolume_N0 + crk(istep) * dt * dvdt
    !         !% limit air volume
    !         airVolume = max(airVolume, zeroR)

    !     end if

    ! end subroutine airpockets_volume_update
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine add_airpocket_heads_to_CC_elem (sc_Idx, aIdx)
        !%------------------------------------------------------------------
        !% Description:
        !%  add the airpocket heads to corresponding elements
        !%------------------------------------------------------------------
        integer, intent (in) :: sc_Idx, aIdx
        real(8), pointer     :: gaugeHead, elemHead(:), hasAirPocket(:), AirHead(:)
        integer, pointer     :: fUp(:), fDn(:)
        logical, pointer     :: isAirPocket, JunctionAirPocket
        integer, allocatable :: pElem(:)
        !%------------------------------------------------------------------
        !% Aliases
            gaugeHead    => airR(sc_Idx,aIdx,airR_gauge_head)
            isAirPocket  => airYN(sc_Idx,aIdx,airYN_air_pocket_detected)
            elemHead     => elemR(:,er_Head)
            AirHead      => elemR(:,er_Air_Pressure_Head)
            hasAirPocket => elemR(:,er_Pressurized_Air)
            JunctionAirPocket => airYN(sc_Idx,aIdx,airYN_junction_air_pocket)
            fUp          => elemI(:,ei_Mface_uL)
            fDn          => elemI(:,ei_Mface_dL)

        !% add the gauge heads to elem heads
        if (isAirPocket .and. .not. JunctionAirPocket) then
            !% pack the elements that contain air pocket
            pElem = pack(conduitElemMapsI(sc_Idx,:,cmi_elem_idx), conduitElemMapsI(sc_Idx,:,cmi_airpocket_idx)  == aIdx)
            !% add the air head to piezometric head
            elemHead(pElem)     = max(elemHead(pElem) + gaugeHead, elemHead(pElem))
            AirHead(pElem)      = gaugeHead
            hasAirPocket(pElem) = oneR
            elemYN(pElem,eYN_hasAirPocket) = .true.
            !% faceYN for air pressurization
            faceYN(fUp(pElem),fYN_isAirPressurized) = .true.
            faceYN(fDn(pElem),fYN_isAirPressurized) = .true.
            !% this is solely for plotting 
            where (elemR(pElem,er_Volume) .ge. elemR(pElem,er_FullVolume)) 
                hasAirPocket(pElem) = zeroR
            end where

            !% deallocate the temporary array
            deallocate(pElem)
        end if

    end subroutine add_airpocket_heads_to_CC_elem
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine airpocket_face_pack ()
        !%------------------------------------------------------------------
        !% Description:
        !%  mask and pack the faces adjacent to element containing airpockets
        !%------------------------------------------------------------------
        integer, pointer     :: ptype, npack, fIdx(:)
        !%------------------------------------------------------------------
        !% Aliases
        fIdx   => faceI(:,fi_Lidx)
        ptype  => col_faceP(fp_Airpockets)
        npack  => npack_faceP(ptype)

        npack = count(faceYN(:,fYN_isAirPressurized))

        if (npack > 0) then 
            faceP(1:npack,ptype) = pack(fIdx,                        &
                      faceYN(:,fYN_isAirPressurized))
        end if

    end subroutine airpocket_face_pack
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

                                !% save the airpocket detection at the conduit
                                conAir(cIdx) = .true.
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
                airYN(ii,:,airYN_air_pocket_detected)  = .false.
                airYN(ii,:,airYN_air_pocket_collapsed) = .true.
                airYN(ii,:,airYN_new_air_pocket)       = .false.
                airI(ii,:,airI_idx)                    = nullvalueI
                airI(ii,:,airI_elem_start)             = nullvalueI
                airI(ii,:,airI_face_up)                = nullvalueI
                airI(ii,:,airI_face_dn)                = nullvalueI
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
    subroutine airpockets_detection_static () 
        !%------------------------------------------------------------------
        !% Description:
        !%      # airpockets per link
        !%------------------------------------------------------------------
            integer          :: ii, jj, kk, startIdx, endIdx, nElem, airPocketIdx, s1
            integer, pointer :: cIdx(:), eIdx(:), fUp(:), fDn(:), max_airpockets
            logical, pointer :: conAir(:), elemSur(:), fBlocked(:), StaticAirPocket
            logical          :: possibleAirPocket, contigious_pocket
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

            !% ------------------------------------------------------------------
            !% initial air pockets screening
            !% search for if any faces conduit elements that has blocked airflow
            if (any(elemSur(eIdx))) then
                possibleAirPocket = .true.
            end if

            !% if initially any airPocket is detected, map the location of that air Pocket
            if (possibleAirPocket) then

                !% check if airpockets in this vlink already been detected
                if (any(airYN(ii,:,airYN_air_pocket_detected))) then
                    !% air pocket already detected
                    conAir(cIdx) = .true.

                    do kk = 1,max_airpockets

                        !% --- special conditions for upstream and downstream 
                        !%     airflow blockage
                        if (airI(ii,kk,airI_type) == upReleaseAirpocket)  then
                            !% if the airflow of upstream release is blocked 
                            if (fBlocked(fUp(oneI))) then
                                airI(ii,kk,airI_type)  = entrappedAirpocket
                            end if
                            
                        else if (airI(ii,kk,airI_type) == dnReleaseAirpocket) then
                            !% if the airflow of downstream release is blocked
                            if (fBlocked(fDn(nElem))) then
                                airI(ii,kk,airI_type)  = entrappedAirpocket
                            end if

                        !% HACK: NOT SURE ABOUT THIS
                        else if (airI(ii,kk,airI_type) == entrappedAirpocket) then
                            !% if the airflow of downstream release is not blocked
                            if (.not. fBlocked(airI(ii,kk,airI_face_dn))) then
                                airI(ii,kk,airI_type)  = dnReleaseAirpocket
                            !% if the airflow of upstream release is not blocked
                            else if (.not. fBlocked(airI(ii,kk,airI_face_up))) then
                                airI(ii,kk,airI_type)  = upReleaseAirpocket

                            ! else if (.not. fBlocked(airI(ii,kk,airI_face_dn)) .and. &
                            !          .not. fBlocked(airI(ii,kk,airI_face_up))) then
                            !     airI(ii,kk,airI_type)  = noAirPocket  
                            !     airYN(ii,kk,airYN_air_pocket_detected) = .false. 
                            end if
                        end if

                    end do


                !% all the airpockets in this vlink has vented or collapsed
                !% detect/re-detect air pockets
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
                    !% initialize continious airpocket continuity
                    contigious_pocket = .false.
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

                            !% store type specific maps

                            !% --- set the type of the airpocket
                            !%     if the starting element is the first element in the conduit
                            !%     there will be an upstream release
                            if (startIdx == oneI) then 
                        
                                !% reset the continious entrapped airpocket identifier
                                if (contigious_pocket) then
                                    contigious_pocket = .false.
                                    !% reset the airpocket index
                                    airPocketIdx = airPocketIdx - oneI
                                end if

                                !% finally save data in airI, airYN arrays
                                !% new air pocket identifier
                                if (.not. airYN(ii,airPocketIdx,airYN_air_pocket_detected)) then
                                    airYN(ii,airPocketIdx,airYN_new_air_pocket) =  .true.
                                end if

                                !% save the airpocket detection logical
                                airYN(ii,airPocketIdx,airYN_air_pocket_detected) = .true.
                                airI(ii,airPocketIdx,airI_type)       = upReleaseAirpocket
                                !% arrayI_pos selects the right column positions of airpockets
                                !% save the integer air pocket data
                                airI(ii,airPocketIdx,airI_idx)        = airPocketIdx
                                airI(ii,airPocketIdx,airI_elem_start) = eIdx(oneI)
                                airI(ii,airPocketIdx,airI_face_up)    = fUp(oneI) 
                                airI(ii,airPocketIdx,airI_elem_end)   = eIdx(endIdx) 
                                airI(ii,airPocketIdx,airI_face_dn)    = fDn(endIdx)

                                !% save airpocket detection data in the conduitElemMapsI array
                                conduitElemMapsI(ii,startIdx:endIdx,cmi_airpocket_idx)  = airPocketIdx
                                conduitElemMapsI(ii,startIdx:endIdx,cmi_airpocket_type) = airI(ii,airPocketIdx,airI_type)

                            !% else if the ending element is the last element in the conduit
                            !% there will be a downstream release
                            else if (endIdx == nElem) then
                        
                                !% reset the continious entrapped airpocket identifier
                                if (contigious_pocket) then
                                    contigious_pocket = .false.
                                    !% reset the airpocket index
                                    airPocketIdx = airPocketIdx - oneI
                                end if

                                !% finally save data in airI, airYN arrays
                                !% new air pocket identifier
                                if (.not. airYN(ii,airPocketIdx,airYN_air_pocket_detected)) then
                                    airYN(ii,airPocketIdx,airYN_new_air_pocket) =  .true.
                                end if

                                !% save the airpocket detection logical
                                airYN(ii,airPocketIdx,airYN_air_pocket_detected) = .true.
                                airI(ii,airPocketIdx,airI_type)       = dnReleaseAirpocket
                                !% arrayI_pos selects the right column positions of airpockets
                                !% save the integer air pocket data
                                airI(ii,airPocketIdx,airI_idx)        = airPocketIdx
                                airI(ii,airPocketIdx,airI_elem_start) = eIdx(startIdx)
                                airI(ii,airPocketIdx,airI_face_up)    = fUp(startIdx) 
                                airI(ii,airPocketIdx,airI_elem_end)   = eIdx(nElem) 
                                airI(ii,airPocketIdx,airI_face_dn)    = fDn(nElem)

                                !% save airpocket detection data in the conduitElemMapsI array
                                conduitElemMapsI(ii,startIdx:endIdx,cmi_airpocket_idx)  = airPocketIdx
                                conduitElemMapsI(ii,startIdx:endIdx,cmi_airpocket_type) = airI(ii,airPocketIdx,airI_type)
                            
                            !% else the pocket is entrapped
                            else
                                if (.not.  contigious_pocket) then

                                    !% finally save data in airI, airYN arrays
                                    !% new air pocket identifier
                                    if (.not. airYN(ii,airPocketIdx,airYN_air_pocket_detected)) then
                                        airYN(ii,airPocketIdx,airYN_new_air_pocket) =  .true.
                                    end if

                                    !% save the airpocket detection logical
                                    airYN(ii,airPocketIdx,airYN_air_pocket_detected) = .true.
                                    airI(ii,airPocketIdx,airI_type)       = entrappedAirpocket
                                    !% arrayI_pos selects the right column positions of airpockets
                                    !% save the integer air pocket data
                                    airI(ii,airPocketIdx,airI_idx)        = airPocketIdx
                                    airI(ii,airPocketIdx,airI_elem_start) = eIdx(startIdx)
                                    airI(ii,airPocketIdx,airI_face_up)    = fUp(startIdx) 
                                    airI(ii,airPocketIdx,airI_elem_end)   = eIdx(endIdx) 
                                    airI(ii,airPocketIdx,airI_face_dn)    = fDn(endIdx)

                                    !% save airpocket detection data in the conduitElemMapsI array
                                    conduitElemMapsI(ii,startIdx:endIdx,cmi_airpocket_idx)  = airPocketIdx
                                    conduitElemMapsI(ii,startIdx:endIdx,cmi_airpocket_type) = airI(ii,airPocketIdx,airI_type)

                                    !% set the pocket continuity to true
                                    contigious_pocket = .true.
                                    s1 = startIdx
                                    
                                else
                                    !% reset the airpocket index
                                    airPocketIdx = airPocketIdx - oneI
                                    !% update the airpocket end maps
                                    airI(ii,airPocketIdx,airI_elem_end)   = eIdx(endIdx) 
                                    airI(ii,airPocketIdx,airI_face_dn)    = fDn(endIdx)

                                    !% save airpocket detection data in the conduitElemMapsI array
                                    conduitElemMapsI(ii,s1:endIdx,cmi_airpocket_idx)  = airPocketIdx
                                    conduitElemMapsI(ii,s1:endIdx,cmi_airpocket_type) = airI(ii,airPocketIdx,airI_type)

                                end if

                            end if

                            !% --- air pocket limit
                            !%     if the airpocket limit doesnot exceed the permissible
                            !%     number of airpockets per conduit, store the data
                            if ((airPocketIdx > zeroI) .and. (airPocketIdx <= threeI)) then

                                
                                !% save the airpocket detection at the conduit
                                conAir(cIdx) = .true.

                                !% --- special conditions for upstream and downstream 
                                !%     airflow blockage (gate or valves)
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
                                
                            else
                                print*, 'THE ALGORITHM HAS FAILED!!!'
                                stop 7897561
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
                !% airYN
                airYN(ii,:,airYN_air_pocket_detected)  = .false.
                airYN(ii,:,airYN_air_pocket_collapsed) = .true.
                airYN(ii,:,airYN_new_air_pocket)       = .false.
                airYN(ii,:,airYN_junction_air_pocket)  = .false.
                !% airI
                airI(ii,:,airI_idx)                    = nullvalueI
                airI(ii,:,airI_elem_start)             = nullvalueI
                airI(ii,:,airI_elem_end)               = nullvalueI
                airI(ii,:,airI_face_up)                = nullvalueI
                airI(ii,:,airI_face_dn)                = nullvalueI
                airI(ii,:,airI_type)                   = noAirPocket
                !% airR
                airR(ii,:,:)                  = zeroR
                airR(ii,:,airR_absolute_head) = setting%AirTracking%AtmosphericPressureHead
                airR(ii,:,airR_density)       = setting%AirTracking%AirDensity
                !% conduitElemMapsI
                conduitElemMapsI(ii,:,cmi_airpocket_idx)  = nullvalueI
                conduitElemMapsI(ii,:,cmi_airpocket_type) = noAirPocket
            end if

            !% debug printing
            ! if (any(airYN(ii,:,airYN_air_pocket_detected))) then
            !     print*, 'air at conduit indexes ',cIdx, ' '!, link%names(cIdx(ii))%str
            !     print*, 'airPocketIdx', airPocketIdx
            !     print*, 'idx 1     = ', airI(ii,1,airI_idx),         'idx 2     = ', airI(ii,2,airI_idx),        'idx 3     = ',airI(ii,3,airI_idx) 
            !     print*, 'type 1    = ', airI(ii,1,airI_type),        'type 2    = ', airI(ii,2,airI_type),       'type 3    = ', airI(ii,3,airI_type)
            !     print*, 'elem up 1 = ', airI(ii,1,airI_elem_start),  'elem up 2 = ', airI(ii,2,airI_elem_start), 'elem up 3 = ',airI(ii,3,airI_elem_start) 
            !     print*, 'elem dn 1 = ', airI(ii,1,airI_elem_end),    'elem dn 2 = ', airI(ii,2,airI_elem_end),   'elem dn 3 = ',airI(ii,3,airI_elem_end) 
            !     ! print*, 'Head 1= ', airR(2,1,airR_absolute_head), 'Head 2= ', airR(2,2,airR_absolute_head), 'Head 3= ', airR(2,3,airR_absolute_head)
            !     print* 
            ! end if

        end do

    end subroutine airpockets_detection_static
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function airpockets_JM_outflow (JMidx)
        !%------------------------------------------------------------------
        !% Description:
        !% computes the mass outflow rate from a JM junction
        !% POSITIVE VALUES ARE OUTFLOWS
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx
            real(8), pointer    :: atmHead, ventArea, dishCoeff, rho_a, rho_w
            real(8), pointer    :: kappa, grav
            real(8)             :: absHead, ratio, ExpansionFactor
        !%------------------------------------------------------------------
        !% --- pointers
        ventArea  => setting%AirTracking%MinimumVentArea
        atmHead   => setting%AirTracking%AtmosphericPressureHead
        dishCoeff => setting%AirTracking%AirDischargeCoefficient
        rho_a     => setting%AirTracking%AirDensity
        rho_w     => setting%AirTracking%WaterDensity
        kappa     => setting%AirTracking%PolytropicExponent
        grav      => setting%constant%gravity

        !% --- calculate the absolute head
        absHead = elemSR(JMidx,esr_JM_Air_HeadAbsolute)
        !% --- calculate the ratio to the atm head to pressure head at the junction
        ratio = atmHead / absHead
        ! !% ---- temporarily set to zero
        ! airpockets_JM_outflow = zeroR

        !% --- calculate he air outflow
        if (absHead >= 1.893 * atmHead) then
            !% chocked air expulsion form the valve
            if (elemSR(JMidx,esr_JM_Air_Density) > zeroR) then
                airpockets_JM_outflow =   dishCoeff * ventArea * elemSR(JMidx,esr_JM_Air_Density) &
                                        * sqrt(grav * (rho_w / elemSR(JMidx,esr_JM_Air_Density) ) &
                                        * absHead) * sqrt(kappa * (twoR / (kappa + oneR)) ** ((kappa + oneR) / (kappa - oneR)))
            else
                airpockets_JM_outflow = zeroR
            end if
        else if (absHead > atmHead .and. absHead < 1.893 * atmHead) then
            !% subsonic air expulsion from the valve
            if (elemSR(JMidx,esr_JM_Air_Density) > zeroR) then
                !% find the expansion factor
                ExpansionFactor = sqrt((kappa / (kappa - oneR)) * (ratio ** (twoR / kappa)) &
                                    * ((oneR - ratio ** ((kappa - oneR) / kappa)) / (oneR - ratio)))
                !% find the airflow
                airpockets_JM_outflow =   dishCoeff * ventArea * ExpansionFactor * elemSR(JMidx,esr_JM_Air_Density) &
                                      * sqrt(twoR * grav * (rho_w / elemSR(JMidx,esr_JM_Air_Density) ) * abs(absHead - atmHead))
            else
                airpockets_JM_outflow = zeroR
            end if 
    
        else
            airpockets_JM_outflow = zeroR
        end if

    end function airpockets_JM_outflow
!%    
!%==========================================================================
!%==========================================================================
!%
    logical function airpocket_detection_JM (JMidx)
        !%------------------------------------------------------------------
        !% Description:
        !% Detects if an air pocket is present from the junction branches
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx
            integer             :: ii
            real(8), pointer    :: depth(:), zbtm(:), zCrown(:)    
        !%------------------------------------------------------------------
        !% pointers
        depth  => elemR(:,er_Depth)
        zbtm   => elemR(:,er_Zbottom)
        zCrown => elemR(:,er_Zcrown)
        !%------------------------------------------------------------------
        airpocket_detection_JM = .true.

        do ii=1,max_branch_per_node

            if (elemSI(JMidx+ii,esi_JB_Exists) .ne. oneI) cycle

            !% if any of the junction is not surcharged, 
            !% or doesnot contain an airpocket, return false
            if ((depth(JMidx+ii) + zbtm(JMidx+ii) < zCrown(JMidx+ii)) .and. &
                (.not. elemYN(JMidx+ii,eYN_hasAirPocket))) then
                airpocket_detection_JM = .false.
                return
            end if 
            
        end do

    end function airpocket_detection_JM 
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function airpockets_JM_airvolume (JMidx)
        !%------------------------------------------------------------------
        !% Description:
        !% calculate volume in junction
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx
            integer             :: ii, sc_idx, aIdx
        !%------------------------------------------------------------------
        !% calculate the initial air volume in the junction
        airpockets_JM_airvolume = elemR(JMidx,er_FullVolume) - elemR(JMidx,er_Volume)
    
        !% add the air volumes from JM adjacent air pockets
        do ii=1,max_branch_per_node
            if (elemSI(JMidx+ii,esi_JB_Exists) .ne. oneI) cycle

            !% add the volume form adjacent air pockets
            if (elemYN(JMidx+ii,eYN_hasAirPocket)) then
                sc_idx = elemSI(JMidx+ii,esi_JB_vLink_Connection)
                aIdx   = elemSI(JMidx+ii,esi_JB_air_pocket_index)
                airpockets_JM_airvolume = airpockets_JM_airvolume + airR(sc_Idx,aIdx,airR_volume)  

            end if 
        end do

    end function airpockets_JM_airvolume
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function airpockets_JM_dvdt (JMidx, istep)
        !%------------------------------------------------------------------
        !% Description:
        !% calculate dvdt from air arrays
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx, istep
            integer             :: ii, sc_idx, aIdx
            real(8), pointer    :: dt, crk
        !%------------------------------------------------------------------
        dt      => setting%Time%Hydraulics%Dt
        crk     => setting%Solver%crk2(istep)
        !% --- compute the air volume rate of change in JM main first 
        !%     (opposite of water volume change) Negative is decreasing volume
        ! airpockets_JM_dvdt  = (elemR(JMidx,er_Volume_N0) - elemR(JMidx,er_Volume)) / (crk * dt)

        airpockets_JM_dvdt  = - elemR(JMidx,er_FlowrateLateral)

        !% --- handle the upstream branches
        do ii=1,max_branch_per_node,2
            if (elemSI(JMidx+ii,esi_JB_Exists) .ne. oneI) cycle
            
            !% add the upstream branch inflow 
            if (elemYN(JMidx+ii,eYN_hasAirPocket)) then
                sc_idx = elemSI(JMidx+ii,esi_JB_vLink_Connection)
                aIdx   = elemSI(JMidx+ii,esi_JB_air_pocket_index)
                airpockets_JM_dvdt = airpockets_JM_dvdt - airR(sc_Idx,aIdx,airR_inflow)  
            end if 
        end do

        !% --- handle the downstream branches
        do ii=2,max_branch_per_node,2
            if (elemSI(JMidx+ii,esi_JB_Exists) .ne. oneI) cycle
            
            !% add the upstream branch inflow 
            if (elemYN(JMidx+ii,eYN_hasAirPocket)) then
                sc_idx = elemSI(JMidx+ii,esi_JB_vLink_Connection)
                aIdx   = elemSI(JMidx+ii,esi_JB_air_pocket_index)
                airpockets_JM_dvdt = airpockets_JM_dvdt + airR(sc_Idx,aIdx,airR_outflow)  
            end if 
        end do
        
    end function airpockets_JM_dvdt
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine airpockets_update_air_pressure_from_JM (istep)
        !%------------------------------------------------------------------
        !% Description:
        !% add the air pressures from JM to airR
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in)  :: istep
            integer              :: mm, JMidx,ii, sc_idx, aIdx
            integer, pointer     :: Npack, thisJM(:)
            real(8), pointer     :: elemHead(:), hasAirPocket(:), AirHead(:)
            integer, pointer     :: fUp(:), fDn(:)
            integer, allocatable :: pElem(:)
        !%------------------------------------------------------------------
        !% Preliminaries
            Npack => npack_elemP(ep_JM)   
            if (Npack < 1) return
        !% pointers
            thisJM       => elemP(1:Npack,ep_JM)
            elemHead     => elemR(:,er_Head)
            AirHead      => elemR(:,er_Air_Pressure_Head)
            hasAirPocket => elemR(:,er_Pressurized_Air)
            fUp          => elemI(:,ei_Mface_uL)
            fDn          => elemI(:,ei_Mface_dL)

        !% --- cycle through the JM junctions
        do mm = 1,Npack
            JMidx = thisJM(mm)
            !% add the air volumes from JM adjacent air pockets
            do ii=1,max_branch_per_node
                if (elemSI(JMidx+ii,esi_JB_Exists) .ne. oneI) cycle

                !% add the volume form adjacent air pockets
                if (elemYN(JMidx+ii,eYN_hasAirPocket)) then
                    !% airpocket maps
                    sc_idx = elemSI(JMidx+ii,esi_JB_vLink_Connection)
                    aIdx   = elemSI(JMidx+ii,esi_JB_air_pocket_index)

                    airR(sc_Idx,aIdx,airR_gauge_head)    = elemSR(JMidx,esr_JM_Air_HeadGauge)
                    airR(sc_Idx,aIdx,airR_absolute_head) = elemSR(JMidx,esr_JM_Air_HeadAbsolute) 
                    airR(sc_Idx,aIdx,airR_mass)          = elemSR(JMidx,esr_JM_Air_Mass)
                    airR(sc_Idx,aIdx,airR_density)       = elemSR(JMidx,esr_JM_Air_Density) 


                    !% pack the elements that contain air pocket
                    pElem = pack(conduitElemMapsI(sc_Idx,:,cmi_elem_idx), conduitElemMapsI(sc_Idx,:,cmi_airpocket_idx)  == aIdx)
                    !% add the air head to piezometric head
                    elemHead(pElem)     = max(elemHead(pElem) + airR(sc_Idx,aIdx,airR_gauge_head), elemHead(pElem))
                    AirHead(pElem)      = airR(sc_Idx,aIdx,airR_gauge_head)
                    hasAirPocket(pElem) = oneR
                    elemYN(pElem,eYN_hasAirPocket) = .true.
                    !% faceYN for air pressurization
                    faceYN(fUp(pElem),fYN_isAirPressurized) = .true.
                    faceYN(fDn(pElem),fYN_isAirPressurized) = .true.
                    !% this is solely for plotting 
                    where (elemR(pElem,er_Volume) .ge. elemR(pElem,er_FullVolume)) 
                        hasAirPocket(pElem) = zeroR
                    end where

                    !% deallocate the temporary array
                    deallocate(pElem)
                end if 
            end do
        end do

    end subroutine airpockets_update_air_pressure_from_JM
!%    
!%==========================================================================
!%==========================================================================
!%
end module air_entrapment