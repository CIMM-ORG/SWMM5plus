module junction_lowlevel
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Lower-level procedures used for junction computation
    !%==========================================================================
    use define_globals
    use define_keys
    use define_indexes
    use define_xsect_tables
    use define_settings, only: setting
    use face, only: face_push_elemdata_to_face
    use geometry, only: geo_depth_from_volume_by_element_CC
    use geometry_lowlevel, only: llgeo_head_from_depth_pure
    use update, only: update_Froude_number_element, update_wavespeed_element, update_auxiliary_variables_CC
    use utility_crash, only: util_crashpoint

    implicit none

    private

    public :: lljunction_branch_velocity
    !public :: lljunction_branch_energy_outflow_OLD 
    public :: lljunction_branch_dQdH
    !public :: lljunction_branch_getface
    public :: lljunction_branch_Qnet
    public :: lljunction_branch_update_DeltaQ
    public :: lljunction_branch_update_flowrate
    
    public :: lljunction_CC_for_JBadjacent

    public :: lljunction_conservation_residual
    public :: lljunction_conservation_fix

    public :: lljunction_main_dHcompute
    public :: lljunction_main_dryingfix

    public :: lljunction_main_update_intermediate
    public :: lljunction_main_update_final

    public :: lljunction_main_dQdHoverflow
    ! public :: lljunction_main_dQdHstorage
    public :: lljunction_main_energyhead
    !public :: lljunction_main_head_bounds
    ! public :: lljunction_main_iscrossing_overflow_or_ponding
    ! public :: lljunction_main_iscrossing_surcharge
    public :: lljunction_main_netFlowrate
    public :: lljunction_main_overflow_conditions
    public :: lljunction_main_plan_area
    public :: lljunction_main_Qoverflow
    public :: lljunction_main_slotwidth
    public :: lljunction_main_sumBranches

    public :: lljunction_main_update_Qdependent_values
    !public :: lljunction_main_update_storage_rate 
    public :: lljunction_main_velocity
    public :: lljunction_main_volume_from_storageRate

    public :: lljunction_push_inflows_from_CC_to_JB_face
    public :: lljunction_push_adjacent_CC_elemdata_to_face


    integer :: printJM = 135
    integer :: printJB = 136
    
    integer :: stepCut = 57000
    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%    
subroutine lljunction_branch_velocity ()
    !%-----------------------------------------------------------------
        !% Description:
        !% Computes an energy-equation outflow for each outflow junction 
        !% branch. This is needed to ensure that zero outflow do not 
        !% become stuck.
        !%-----------------------------------------------------------------
        !% Declarations:
        integer, pointer :: JBidx, JMidx, Npack, thisP(:)
        integer, pointer :: fidx
        real(8), pointer :: HeadJM, HeadAdj
        real(8), pointer :: DepthAdj, ZbottomJB, Ke
        real(8), pointer :: VelocityJM, VelocityAdj, grav
        real(8), pointer :: BlendingFactor, ReverseDhFactor
        real(8), pointer :: EnergyHeadJM, EnergyHeadAdj
        real(8) :: deltaH, deltaE, bsign, eFlowrate, eVelocity, VelHead
        real(8) :: deltaEjmZ, deltaEAdjZ, VelHeadJM
        integer :: ii
        logical :: isOutflow, isUpstream

        !%------------------------------------------------------------------
        !% Aliases
            Npack => npack_elemP(ep_JB)
            grav  => setting%Constant%gravity
            BlendingFactor => setting%Junction%BlendingFactor
            ReverseDhFactor=> setting%Junction%ReverseDhFactor
        !%------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return   
        !%------------------------------------------------------------------

        thisP => elemP(1:Npack,ep_JB)

        !% --- cycle through branches
        do ii=1,Npack 
            JBidx => thisP(ii)
            if (elemSI(JBidx,esi_JB_Exists) .ne. oneI) cycle 

            JMidx        => elemSI(JBidx,esi_JB_Main_Index)
            HeadJM       => elemR (JMidx,er_Head)
            EnergyHeadJM => elemR (JMidx,er_EnergyHead)
            VelocityJM   => elemR (JMidx,er_Velocity)           
            
            ZbottomJB    => elemR (JBIdx,er_Zbottom)
            Ke           => elemSR(JBidx,esr_JB_Kfactor)

            VelHeadJM  = (VelocityJM**2) / (twoR * grav) !% always > 0
            
            !% --- get the up or down face for this JB
            !%     and set whether adjacent face indicates this as outflow or inflow
            if (elemSI(JBidx,esi_JB_IsUpstream) == oneI) then 
                isUpstream = .true.
                bsign = +oneR
                !% --- upstream branch
                fidx => elemI(JBidx,ei_Mface_uL)
                !% --- use face adjacent velocity to determine in/outflow
                if (faceR(fidx,fr_Velocity_Adjacent) .le. zeroR) then 
                    isOutflow = .true. 
                else
                    isOutflow = .false.
                end if
            else
                isUpstream = .false.
                bsign = -oneR
                !% --- downstream branch
                fidx => elemI(JBidx,ei_Mface_dL)
                !% --- use face adjacent velocity to determine in/outflow
                if (faceR(fidx,fr_Velocity_Adjacent) .ge. zeroR) then 
                    isOutflow = .true. 
                else
                    isOutflow = .false.
                end if
            endif
            HeadAdj       => faceR(fidx,fr_Head_Adjacent)
            EnergyHeadAdj => faceR(fidx,fr_EnergyHead_Adjacent)
            VelocityAdj   => faceR(fidx,fr_Velocity_Adjacent)

            deltaE     = EnergyHeadJM  - EnergyHeadAdj 
            deltaEjmZ  = EnergyHeadJM  - ZbottomJB
            deltaEAdjZ = EnergyHeadAdj - ZbottomJB

            if (isOutflow) then 
                if (deltaEjmZ .le. zeroR) then 
                    !% --- JM energy head below branch bottom, no outflow possible 
                    elemR(JBidx,er_Flowrate) = zeroR
                    elemR(JBidx,er_Velocity) = zeroR
                    cycle
                else
                    !% --- JM energy head above branch bottom allows an outflow
                    if (HeadAdj > ZbottomJB) then 
                        !% --- Adj Head is also above branch bottom so we have fluid connection
                        if (deltaE > zeroR) then 
                            !% --- simple outflow, which should be Vhead > 0
                            VelHead = EnergyHeadJM - HeadAdj - Ke * VelHeadJM
                            if (VelHead .le. zeroR) then 
                                !% --- VelocityJM must be an overestimate, so use a correction 
                                !%     Note that >0 is velocity head out
                                VelHead = ReverseDhFactor * deltaE
                            else
                                !% --- retain VelHead
                            end if
                        else 
                            !% --- inconsistent deltaE <= 0 implies inflow, but adjacent velocity implies outflow
                            !%     Note that <0 implies velocity head out occurs if JM < adj
                            VelHead = ReverseDhFactor * (EnergyHeadJM - HeadAdj)
                        end if
                    else
                        !% --- Adj Head is below branch bottom
                        !% --- Outflow from JM is a waterfall into Adj governed by deltaEjmZ
                        !%     Here >0 implies velocity head out
                        VelHead = deltaEjmZ - Ke * VelHeadJM
                        if (VelHead .le. zeroR) then 
                            !% --- VelocityJM is an overestimate, so use a correction 
                            !%     Here >0 implies velocity head out
                            VelHead = ReverseDhFactor * deltaEjmZ
                        else
                            !% --- retain VelHead
                        end if
                    end if
                end if
            else
                !% --- inflow
                if (deltaEAdjZ .le. zeroR) then 
                    !% --- adjacent energy head below branch bottom, so no inflow possible
                    elemR(JBidx,er_Flowrate) = zeroR
                    elemR(JBidx,er_Velocity) = zeroR
                    cycle
                else
                    !% --- adjacent energy head above branch bottom
                    if (HeadJM > ZbottomJB) then
                        !% --- JM head is above branch bottom so we have a fluid connection
                        !%     Key difference from Outflow algorithm is the the frictional loss depends
                        !%     on the velocity head adjacent, which we are solving for, so we have
                        !%     the 1/(1-Ke) form. Arguably, this should use EnergyHeadJM, but since
                        !%     the energy head in JM might be driven by another inflow, we will neglect
                        !%     it here. Might explore using a fraction of it depending on relavitve
                        !%     flowrate contributions to JM.
                        !%     Here we expect <0  as inflow
                        VelHead = -(HeadAdj - HeadJM) / (oneR - Ke)
                        if (VelHead .ge. zeroR) then 
                            !% --- inconsistent values
                            VelHead = ReverseDhFactor * deltaE
                        else
                            !% --- retain VelHead 
                        end if
                    else
                        !% --- HeadJM < Zbottom so a waterfall connection from adjacent inward
                        !%     >0 implies inflow
                        VelHead = -deltaEAdjZ / (oneR - Ke)
                        !% --- no need to check if deltaEAdjZ < 0 because that has already been removed 
                    end if
                end if
            end if

            !% --- velocity and flowrate implied by energy arguments
            eVelocity = - bsign * sign(oneR,VelHead) * sqrt(abs(VelHead) * twoR * grav)  
            eFlowrate = eVelocity * elemR(JBidx,er_AreaVelocity)  

            !% --- blend energy flowrate with flowrate stored in JB branch
            elemR(JBidx,er_Flowrate) =  (oneR - BlendingFactor) * eFlowrate                &
                                              + BlendingFactor  * elemR(JBidx,er_Flowrate)   

            !% --- compute the velocity
            if (elemR(JBidx,er_AreaVelocity) > setting%ZeroValue%Area) then
                elemR(JBidx,er_Velocity) = elemR(JBidx,er_Flowrate) / elemR(JBidx,er_AreaVelocity)
            else
                elemR(JBidx,er_Velocity) = zeroR
            end if

            !% --- apply velocity limiter (does not affect flowrate)
            if (abs(elemR(JBidx,er_Velocity)) > setting%Limiter%Velocity%Maximum) then 
                elemR(JBidx,er_Velocity) = sign(setting%Limiter%Velocity%Maximum * 0.99d0, &
                                                elemR(JBidx,er_Velocity))
            end if

        end do

        
    end subroutine lljunction_branch_velocity
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine lljunction_branch_energy_outflow_OLD ()
    !     !%-----------------------------------------------------------------
    !     !% Description:
    !     !% Computes an energy-equation outflow for each outflow junction 
    !     !% branch. This is needed to ensure that zero outflow do not 
    !     !% become stuck.
    !     !%-----------------------------------------------------------------
    !     !% Declarations:
    !     integer, pointer :: JBidx, JMidx, Npack, thisP(:)
    !     integer, pointer :: fidx
    !     real(8), pointer :: HeadJM, EnergyHeadJM, HeadAdj, EnergyHeadAdj
    !     real(8), pointer :: DepthAdj
    !     real(8), pointer :: VelocityJM, VelocityAdj, grav, DampingFactor
    !     real(8) :: deltaH, deltaE, bsign, eFlowrate, VelHead
    !     integer :: ii
    !     logical :: isOutflow, isUpstream

    !     !%------------------------------------------------------------------
    !     !% Aliases
    !         Npack => npack_elemP(ep_JB)
    !         grav  => setting%Constant%gravity
    !         DampingFactor => setting%Junction%BlendingFactor
    !     !%------------------------------------------------------------------
    !     !% Preliminaries
    !         if (Npack < 1) return   
    !     !%------------------------------------------------------------------

    !     thisP => elemP(1:Npack,ep_JB)

    !     !% --- cycle through branches
    !     do ii=1,Npack 
    !         JBidx => thisP(ii)
    !         if (elemSI(JBidx,esi_JB_Exists) .ne. oneI) cycle 

    !         JMidx        => elemSI(JBidx,esi_JB_Main_Index)
    !         HeadJM       => elemR (JMidx,er_Head)
    !         EnergyHeadJM => elemR (JMidx,er_EnergyHead)
    !         VelocityJM   => elemR (JMidx,er_Velocity)            
            
    !         !% --- get the up or down face for this JB
    !         if (elemSI(JBidx,esi_JB_IsUpstream) == oneI) then 
    !             isUpstream = .true.
    !             bsign = +oneR
    !             !% --- upstream branch
    !             fidx => elemI(JBidx,ei_Mface_uL)
    !             !% --- use face adjacent velocity to determine in/outflow
    !             if (faceR(fidx,fr_Velocity_Adjacent) .le. zeroR) then 
    !                 isOutflow = .true. 
    !             else
    !                 isOutflow = .false.
    !             end if
    !         else
    !             isUpstream = .false.
    !             bsign = -oneR
    !             !% --- downstream branch
    !             fidx => elemI(JBidx,ei_Mface_dL)
    !             !% --- use face adjacent velocity to determine in/outflow
    !             if (faceR(fidx,fr_Velocity_Adjacent) .ge. zeroR) then 
    !                 isOutflow = .true. 
    !             else
    !                 isOutflow = .false.
    !             end if
    !         endif
    !         HeadAdj       => faceR(fidx,fr_Head_Adjacent)
    !         EnergyHeadAdj => faceR(fidx,fr_EnergyHead_Adjacent)
    !         VelocityAdj   => faceR(fidx,fr_Velocity_Adjacent)
    !         DepthAdj      => faceR(fidx,fr_Depth_Adjacent)

    !         deltaH = HeadJM - HeadAdj
    !         deltaE = EnergyHeadJM - EnergyHeadAdj

    !         !% --- get the energy equation velocity head at the inlet/outlet
    !         if (isOutflow) then 
    !             if (deltaE .ge. zeroR) then 
    !                 !% --- outflow (positive) with consistent energy gradient
    !                 !%     outflow velocity head losing head based on K factor for approach velocity
    !                 !VelHead =  deltaE - elemSR(JBidx,esr_JB_Kfactor) &
    !                 !                    * (VelocityJM**2) / (twoR * grav)
    !                 VelHead = EnergyHeadJM - HeadAdj - elemSR(JBidx,esr_JB_Kfactor) &
    !                                                     * (VelocityJM**2) / (twoR * grav)                                 
    !                 if (VelHead .le. zeroR) then 
    !                     !% --- K factor head loss is too great, so ad hoc reduction of deltaE
    !                     VelHead = onehalfR * deltaE !% positive is outflow
    !                 end if     
    !             else
    !                 !% --- (deltaE < zeroR)
    !                 !% --- nominal outflow with inconsistent energy gradient (deltaE < 0)
    !                 !%     compute velocity head depending on deltaH
    !                 if (deltaH < zeroR) then 
    !                     !% --- velocity head driven solely by static head, reversing direction
    !                     !%     do NOT use energy because adjacent value has
    !                     !%     inconsistent direction
    !                     VelHead = onehalfR * deltaH  !% negative is inflow
    !                 else
    !                     !% --- inconsistent case that should not occur, deltaE < 0 and deltaH > 0
    !                     !VelHead = deltaH - elemSR(JBidx,esr_JB_Kfactor) &
    !                     VelHead = EnergyHeadJM - HeadAdj - elemSR(JBidx,esr_JB_Kfactor) &
    !                                 * (VelocityJM**2) / (twoR * grav)      
    !                     if (VelHead .le. zeroR) then 
    !                         !% --- K factor head loss is too great, so ad hoc reduction of deltaH
    !                         VelHead = onehalfR * deltaH  !% positive is outflow
    !                     end if
    !                 end if
    !             end if
    !         else
    !             !% --- nominal inflow
    !             if (deltaE .le. zeroR) then 
    !                 !% --- inflow with consistent energy gradient, velocity head is negative
    !                 VelHead = (EnergyHeadJM - HeadAdj) / (oneR - elemSR(JBidx,esr_JB_Kfactor) )
    !                 ! VelHead = -(EnergyHeadAdj - HeadJM - elemSR(JBidx,esr_JB_Kfactor) &
    !                 !                                    * (VelocityAdj**2) / (twoR * grav))
    !                 ! if (VelHead .ge. zeroR) then 
    !                 !     !% --- K factor head loss too great, so ad hoc reduction of delta E
    !                 !     VelHead = onehalfR * deltaE !% negative is inflow
    !                 ! end if
    !             else
    !                 !% --- (deltaE > zeroR) 
    !                 !% --- nominal inflow (dE should be negative) with inconsistent energy gradient
    !                 !%     compute reversed velocity head 
    !                 if (deltaH > zeroR) then 
    !                     !% --- velocity outflow head driven solely by static head
    !                     !%     do NOT use energy because adjacent value has
    !                     !%     inconsistent direction
    !                     VelHead = onehalfR * deltaH  !% positive is outflow 
    !                 else
    !                     !% --- case that should not occur, deltaE > 0 and deltaH < 0
    !                     !%     set inflow (negative) velocity based on deltaH and adjacent velocity losses
    !                     VelHead = deltaH + elemSR(JBidx,esr_JB_Kfactor) &
    !                                         * (VelocityAdj**2) / (twoR * grav)
    !                     if (VelHead .ge. zeroR) then 
    !                         !% --- K factor head loss is too great, so ad hoc reduction of deltaH
    !                         VelHead = onehalfR * deltaH !% negative is inflow
    !                     end if
    !                 end if
    !             end if

    !         endif

    !         !% --- get a flowrate based on the energy
    !         !% --- positive VelHead is outflow, which is negative velocity for upstream
    !         !%     but is positive velocity for downstream
    !         eFlowrate = - bsign * sign(oneR,VelHead) * sqrt(abs(VelHead) * twoR * grav) &
    !                         * elemR(JBidx,er_Area)

    !         ! if ((JBidx == 111) .and. (setting%Time%Step > 60484) ) then
    !         !     print *, 'EFLOW ',eFlowrate, elemR(JBidx,er_Flowrate)
    !         ! end if

    !         !% --- apply only on outflows 
    !         !% NOTE THE ABOVE COMPUTES ALL THE ENERGY FLOWRATES, BUT WE ONLY USE IT
    !         !% FOR OUTFLOWS. FUTURE -- simplify what is done above.
    !         if ((isUpstream .and. (eFlowrate .le. zeroR) ) & 
    !             .or. &
    !             (.not. isUpstream) .and. (eFlowrate .ge. zeroR))  then         

    !                 !% --- combine the energy flow rate with the prior flowrate
    !                 elemR(JBidx,er_Flowrate) =  (oneR - DampingFactor) * eFlowrate                &
    !                                                 + DampingFactor  * elemR(JBidx,er_Flowrate)                

    !                 ! if (eFlowrate .le. zeroR) then 
    !                 !     elemR(JBidx,er_Flowrate) = min(elemR(JBidx,er_Flowrate),eFlowrate)
    !                 ! else
    !                 !     elemR(JBidx,er_Flowrate) = max(elemR(JBidx,er_Flowrate),eFlowrate)
    !                 ! end if
                    
    !                 !% --- compute the velocity
    !                 if (elemR(JBidx,er_Area) > setting%ZeroValue%Area) then
    !                     elemR(JBidx,er_Velocity) = elemR(JBidx,er_Flowrate) / elemR(JBidx,er_Area)
    !                 else
    !                     elemR(JBidx,er_Velocity) = zeroR
    !                 end if

    !                 !% --- apply velocity limiter (does not affect flowrate)
    !                 if (abs(elemR(JBidx,er_Velocity)) > setting%Limiter%Velocity%Maximum) then 
    !                     elemR(JBidx,er_Velocity) = sign(setting%Limiter%Velocity%Maximum * 0.99d0, &
    !                                                     elemR(JBidx,er_Velocity))
    !                 end if

    !         end if
        
    !     end do

    ! end subroutine lljunction_branch_energy_outflow_OLD
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine lljunction_branch_energy_outflow_OLD2 ()
    !     !%-----------------------------------------------------------------
    !     !% Description:
    !     !% Computes an energy-equation outflow for each outflow junction 
    !     !% branch
    !     !%-----------------------------------------------------------------
    !     !% Declarations:
    !         integer, pointer :: Npack, JMar(:), thisJB(:), fidx(:)
    !         real(8), pointer :: deltaH(:), grav, energyQ(:), Vsq2g(:)

    !         logical :: isUpstreamBranch 
    !         real(8) :: bsign
    !         integer :: frHead, frArea, frHeadAdj, ii
    !     !%-----------------------------------------------------------------
    !     !% Aliases
    !         !% --- array for the JM index
    !         JMar   => elemSI(:,esi_JB_Main_Index)    
    !         grav   => setting%Constant%gravity
    !         energyQ=> elemR(:,er_Temp03)
    !         Vsq2g  => elemR(:,er_Temp04)
    !     !%-----------------------------------------------------------------

    !     !% --- cycle through nominal upstream and downstream JB
    !     !%     This should affect only outflow branches with consistent
    !     !%     pressure difference with adjacent CC
    !     do ii=1,2
    !         !% --- get the upstream or downstream JB elements
    !         if (ii==1) then 
    !             !% --- upstream JB
    !             isUpstreamBranch = .true.
    !             bsign = oneR
    !             Npack => npack_elemP(ep_JB_Upstream_CC_Adjacent)
    !             if (Npack > 0) then 
    !                 thisJB => elemP(1:Npack,ep_JB_Upstream_CC_Adjacent)
    !             else
    !                 cycle
    !             end if
    !             !% --- the JB-adjacent face is upstream
    !             fidx => elemI(:,ei_Mface_uL)
    !             frHead = fr_Head_u !% -- use u for jump purposes? QUESTION
    !             frArea = fr_Area_u !% QUESTION
    !             frHeadAdj = fr_Head_Adjacent
    !         else
    !             !% --- downstream JB
    !             bsign = -oneR
    !             isUpstreamBranch = .false.
    !             Npack => npack_elemP(ep_JB_Downstream_CC_Adjacent)
    !             if (Npack > 0) then 
    !                 thisJB => elemP(1:Npack,ep_JB_Downstream_CC_Adjacent)
    !             else 
    !                 cycle
    !             end if
    !             !% --- the JB-adjacent face is downstream
    !             fidx => elemI(:,ei_Mface_dL)
    !             frHead = fr_Head_d !% --- use d for jump purposes? QUESTION
    !             frArea = fr_Area_d !% QUESTION
    !             frHeadAdj = fr_Head_Adjacent
    !         end if

    !         !% --- head difference from junction main to element
    !         !%     note this is positive for any outflow
    !         deltaH  => elemR(:,er_Temp01)
    !         deltaH  = zeroR
    !         energyQ = zeroR
    !         Vsq2g   = zeroR

    !         where (elemR(thisJB,er_Depth) > setting%ZeroValue%Depth)
    !             where (faceR(fidx(thisJB),frHeadAdj) > faceR(fidx(thisJB),fr_Zbottom))
    !                 deltaH(thisJB) =  elemR(JMar(thisJB),er_Head) - faceR(fidx(thisJB),frHeadAdj)
    !             elsewhere
    !                 !% --- where adjacent head is lower than face zbottom
    !                 deltaH(thisJB) = elemR(JMar(thisJB),er_Head) - (faceR(fidx(thisJB),fr_Zbottom)+ setting%ZeroValue%Depth)
    !             endwhere
    !         endwhere

    !         where (deltaH(thisJB) > zeroR)
    !             Vsq2g(thisJB) =  deltaH(thisJB) &
    !                         + (elemR(JMar(thisJB),er_Velocity)**2) * (oneR - elemSR(thisJB,esr_JB_Kfactor))  &
    !                         /(twoR * grav)
    !         elsewhere 
    !             Vsq2g(thisJB) = deltaH(thisJB)
    !         endwhere

             
    !         !print *, 'vsq2g ',170, Vsq2g(170), deltaH(170)
        

    !         ! if (ii==2) then
    !         !     print *, ' '
    !         !     print *, 'in branch energy ', JMar(181)
    !         !     print *, 'heads  ', elemR(JMar(181),er_Head), faceR(fidx(181),frHeadAdj)
    !         !     !print *, fidx(181),frHeadAdj
    !         !     !print *, elemI(181,ei_Mface_dL)
    !         !     print *, 'deltaH ',deltaH(181), setting%Junction%ZeroHeadDiffValue
    !         !     print *, 'flow   ',elemR(181,er_Flowrate)
    !         ! end if


    !         !% --- For dH that increases outflow
    !         !where ((deltaH(thisJB) > setting%Junction%ZeroHeadDiffValue)                       &
    !         ! .and.                                                                       &
    !         !        (elemR(thisJB,er_Flowrate) * bsign < -setting%Junction%ZeroOutflowValue)    &
    !         !        )

    !         !% --- where energy equation gives an outflow
    !         where (Vsq2g(thisJB) .ge. zeroR)
    !                 !% --- outflow from JB  in upstream direction
    !                 !% --- or outflow from JB in downstream direction 
    !                 !%     applies junction main approach velocity brh20230829                    
    !                 energyQ(thisJB) = - bsign * elemR(thisJB,er_Area) * sqrt(twoR * grav * Vsq2g(thisJB))                
    !                 !% --- DAMPING: Average with existing flowrate
    !                 elemR(thisJB,er_Flowrate)  = (oneR - setting%Junction%BlendingFactor) * energyQ(thisJB)  &
    !                                                    + setting%Junction%BlendingFactor  * elemR(thisJB,er_Flowrate) 

    !                 !% --- handle small depths
    !                 where (elemR(thisJB,er_AreaVelocity) > setting%ZeroValue%Area)
    !                     elemR(thisJB,er_Velocity) = elemR(thisJB,er_Flowrate) / elemR(thisJB,er_AreaVelocity)
    !                 elsewhere 
    !                     elemR(thisJB,er_Velocity) = zeroR
    !                 endwhere
    
    !                 !% --- apply strict velocity limiter (does not affect flowrate)
    !                 where (abs(elemR(thisJB,er_Velocity)) > setting%Limiter%Velocity%Maximum)
    !                     elemR(thisJB,er_Velocity) = sign(setting%Limiter%Velocity%Maximum * 0.99d0, elemR(thisJB,er_Velocity))
    !                 endwhere       
    !         endwhere 

    !         !% --- where energy would reverse flow
    !         where ((Vsq2g(thisJB) < zeroR)      &
    !                     .and.                   &
    !                     (elemR(thisJB,er_Flowrate) * bsign < -setting%Junction%ZeroOutflowValue) )
    !             elemR(thisJB,er_Velocity) = zeroR
    !             elemR(thisJB,er_Flowrate) = zeroR
    !         endwhere



    !         ! !% --- for dH that decreases, but does not reverse outflow
    !         ! where ((deltaH(thisJB) < zeroR)                        &
    !         !        .and.                                           &
    !         !        (Vsq2g .ge. zeroR)                              &
    !         !        .and.                                           &
    !         !        (elemR(thisJB,er_Flowrate) * bsign < zeroR)     &
    !         !        )
    !         !        energyQ(thisJB) = - bsign * elemR(thisJB,er_Area) * sqrt(twoR * grav * Vsq2g)  
    !         ! endwhere

    !         ! !% --- for a trivial outflow with a trivial dH, set flow to zero
    !         ! !%     Note -- does not affect inflows or deltaH < zeroR
    !         ! where ((deltaH(thisJB) .ge. zeroR)                                                 &
    !         !         .and.                                                                      &
    !         !         (deltaH(thisJB)                   .le. setting%Junction%ZeroHeadDiffValue) &
    !         !         .and.                                                                      &
    !         !         (elemR(thisJB,er_Flowrate) * bsign  < -setting%Junction%ZeroOutflowValue)  &
    !         !        )
    !         !     elemR(thisJB,er_Flowrate) = zeroR
    !         !     elemR(thisJB,er_Velocity) = zeroR
    !         ! end where

    !         ! if (ii==2) then 
    !         !     print *, ' '
    !         !     print *, 'flowrate ',elemR(181,er_Flowrate)
    !         !     print *, ' '
    !         ! end if
    !         ! !% --- for inconsistent inflow pressure gradient with outflow on adjacent element
    !         ! where  (       (deltaH(thisJB)                     <   -setting%Junction%ZeroHeadDiffValue)  &
    !         !          .and. (elemR(thisJB,er_Flowrate) * bsign .le. -setting%Junction%ZeroOutflowValue) )
    !         !     !!% --- inconsistent flow direction and pressure gradient   
    !         !     !elemR(thisJB,er_Flowrate) = zeroR
                
    !         !     !% --- use approach velocity from manhole with K factor
    !         !     elemR(thisB,er_Velocity) = 
    !         ! endwhere

    !         !% --- push flowrate and velocity to faces
    !         faceR(fidx(thisJB),fr_Flowrate)   = elemR(thisJB,er_Flowrate)
    !         faceR(fidx(thisJB),fr_Velocity_u) = elemR(thisJB,er_Velocity)
    !         faceR(fidx(thisJB),fr_Velocity_d) = elemR(thisJB,er_Velocity)
    !         !% --- check if the elem data has either been pushed 
    !         !%     to a shared face. if so then mark that face
    !         where (faceYN(fidx(thisJB),fYN_isSharedFace))
    !             faceYN(fidx(thisJB),fYN_isSharedFaceDiverged) = .true.
    !         end where
    !     end do

    !     !% --- update JB froude number and wave speed
    !     Npack => npack_elemP(ep_JB)
    !     if (Npack > 0) then 
    !         thisJB => elemP(1:Npack, ep_JB)
    !         call update_Froude_number_element (thisJB) 
    !         call update_wavespeed_element (thisJB)
    !     end if
        
    ! end subroutine lljunction_branch_energy_outflow_OLD2
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine lljunction_branch_dQdH ()
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the branch dQdH from the
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: JBidx, JMidx, Npack, thisP(:)
            integer          :: fidx
            real(8), pointer :: fA, fH, fQ, fZ
            real(8), pointer :: Ladj, Tadj, Hadj, headJM
            real(8), pointer ::  crk, dt, grav
            real(8), pointer :: VelAdj , FrAdj, Dadj, ZBadj
            real(8)          :: bsign, denominator
            logical          :: isInflow, isDownstream
            integer          :: ii

            real(8) :: tempfactor !% used for experiments only 20230425
            integer :: istep=oneI  !% Required dQ/dH used prior to RK
        !%------------------------------------------------------------------
        !% Aliases
            crk  => setting%Solver%crk2(istep)
            dt   => setting%Time%Hydraulics%Dt
            grav => setting%Constant%gravity
            Npack => npack_elemP(ep_JB)
        !%------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return   
        !%------------------------------------------------------------------

        thisP => elemP(1:Npack,ep_JB)

        !% --- cycle through branches to set dQdH
        do ii=1,Npack 
            JBidx => thisP(ii)
            JMidx => elemSI(JBidx,esi_JB_Main_Index)

            headJM => elemR(JMidx,er_Head)

            !% --- get the face for this JB
            if (elemSI(JBidx,esi_JB_IsUpstream) == oneI) then 
                !% --- upstream JB
                isDownstream = .false.
                fidx  =  elemI(JBidx,ei_Mface_uL)
                fA    => faceR(fidx,fr_Area_u)
                fH    => faceR(fidx,fr_Head_u)
                fQ    => faceR(fidx,fr_Flowrate)
                fZ    => faceR(fidx,fr_Zbottom)
                FrAdj => faceR(fidx,fr_Froude_Adjacent)
                if (fQ .ge. 0) then 
                    isInflow = .true.
                else 
                    isInflow = .false.
                end if

            else
                !% --- downstream JB
                isDownstream = .true.
                fidx  =  elemI(JBidx,ei_Mface_dL)
                fA    => faceR(fidx,fr_Area_d)
                fH    => faceR(fidx,fr_Head_d)
                fQ    => faceR(fidx,fr_Flowrate)
                fZ    => faceR(fidx,fr_Zbottom)
                FrAdj => faceR(fidx,fr_Froude_Adjacent)
                if (fQ < 0) then 
                    isInflow = .true.
                else 
                    isInflow = .false.
                end if
            end if
            !% --- adjacent element data
            Ladj  => faceR(fidx,fr_Length_Adjacent)
            Tadj  => faceR(fidx,fr_Topwidth_Adjacent)
            Hadj  => faceR(fidx,fr_Head_Adjacent)
            VelAdj=> faceR(fidx,fr_Velocity_Adjacent)
            Dadj  => faceR(fidx,fr_Depth_Adjacent)
            ZBadj => faceR(fidx,fr_Zcrest_Adjacent)

            !% --- CC elements adjacent to JB
            if (elemSI(JBidx,esi_JB_CC_adjacent) == oneI) then 

                if (isDownstream) then 

                    !% --- downstream branch
                    if ((FrAdj .le. -oneR) .or. (elemR(JBidx,er_FroudeNumber) .le. -oneR)) then 
                        !% --- supercritical inflow
                        elemSR(JBidx,esr_JB_dQdH) = zeroR
                    else 
                        !% --- outflow or subcritical inflow
                        !elemSR(JBidx,esr_JB_dQdH) = + crk * grav * dt * fA / Ladj
                        elemSR(JBidx,esr_JB_dQdH) = + grav * dt * fA / (onehalfR * Ladj)
                    end if

                    !% --- handle waterfall inflow elements
                    if ((isInflow) .and. (headJM < fZ)) then 
                        elemSR(JBidx,esr_JB_dQdH) = zeroR
                    end if

                    !% --- handle uphill outflow
                    if ((.not. isInflow) .and. (headJM < Hadj)) then 
                        elemSR(JBidx,esr_JB_dQdH) = zeroR 
                    end if

                    ! !% --- Limit uphill outflow 
                    ! if ( (.not. isInflow) .and. (headJM < Hadj)) then 
                    !     !% --- outflow into an adverse pressure gradient,
                    !     !%     largest negative allowable dQdH is Q/deltaH
                    !     elemSR(JBidx,esr_JB_dQdH) = min(elemSR(JBidx,esr_JB_dQdH), elemR(JBidx,er_Flowrate) / (Hadj - headJM))
                    ! end if

                    !% --- limit outflow by 1/4 volume of JM
                    !%     Qmax = (1/4) V / dt;  V = Aplan * Depth
                    !%     Q/H = 1/4 (Aplan * Depth) / (dt * Depth) = (1/4) Aplan / dt
                    ! if (.not. isInflow) then   
                    !     elemSR(JBidx,esr_JB_dQdH) &
                    !         = min (elemSR(JBidx, esr_JB_dQdH),  &
                    !                 elemSR(JMidx, esr_Storage_Plan_Area) * onefourthR / dt)
                    ! end if

                else
                    !% --- upstream branch
                    if ((FrAdj .ge. +oneR) .or. (elemR(JBidx,er_FroudeNumber) .ge. +oneR)) then 
                        !% --- supercritical inflow
                        elemSR(JBidx,esr_JB_dQdH) = zeroR
                    else
                        !% --- outflow or subcritical inflow
                        !elemSR(JBidx,esr_JB_dQdH) = - crk* grav * dt * fA / Ladj
                        elemSR(JBidx,esr_JB_dQdH) = - grav * dt * fA / (onehalfR * Ladj)
                    end if

                    !% --- handle waterfall inflow elements
                    if ((isInflow) .and. (headJM < fZ)) then 
                        elemSR(JBidx,esr_JB_dQdH) = zeroR
                    end if

                    !% --- handle uphill outflow
                    if ((.not. isInflow) .and. (headJM < Hadj)) then 
                        elemSR(JBidx,esr_JB_dQdH) = zeroR 
                    end if

                    ! !% --- Limit uphill outflow 
                    ! if ( (.not. isInflow) .and. (headJM < Hadj)) then 
                    !     !% --- outflow into adverse pressure gradient
                    !     !%     largest (positive) dQ/dH is Q/deltaH
                    !     elemSR(JBidx,esr_JB_dQdH) = max(elemSR(JBidx,esr_JB_dQdH), elemR(JBidx,er_Flowrate) / (Hadj-headJM))
                    ! end if

                    !% --- Limit outflow dQdH by 1/4 volume of JM
                    !%     Qmax = (1/4) V / dt;  V = Aplan * Depth
                    !%     Q/H = 1/4 (Aplan * Depth) / (dt * Depth) = (1/4) Aplan / dt
                    !%     note upstream branch dQdH < 0
                    ! if (.not. isInflow) then  
                    !     elemSR(JBidx,esr_JB_dQdH) &
                    !         = max ( elemSR(JBidx, esr_JB_dQdH),  &
                    !                 -elemSR(JMidx, esr_Storage_Plan_Area) * onefourthR / dt)
                    ! end if

                end if

            !% --- Diagnostic element adjacent to JB
            elseif (elemSI(JBidx,esi_JB_Diag_adjacent)) then 
                !% --- dQdH for adjacent diagnostic element
                if (isInflow) then 
                    bsign = -oneR !% increasing head decreases flowrate for inflow
                else
                    bsign = +oneR !% increasing head increases flowrate for outflow
                endif
                if (elemR(JMidx,er_Head) < faceR(fidx,fr_Zcrest_Adjacent)) then 
                    !% -- insufficient head means the junction cannot affect the flow in 
                    !%    diagnostic branch (either in or out flow)
                    elemSR(JBidx,esr_JB_dQdH) = zeroR
                    !cycle
                else
                    elemSR(JBidx,esr_JB_dQdH) = bsign * faceR(fidx,fr_dQdH_Adjacent)
                    !cycle
                end if

            else 
                print *, 'CODE ERROR Unexpected else'
                call util_crashpoint(140987112)
            end if

        end do
        
    end subroutine lljunction_branch_dQdH
!%
!%==========================================================================
! !%==========================================================================
! !% 
!     pure subroutine lljunction_branch_getface (outdata, frCol, JMidx, fiIdx, kstart)
!         !%-----------------------------------------------------------------
!         !% Description
!         !% Stores face data of frCol on outdata element space
!         !% Operates either on upstream or downstream branches, but 
!         !% requires separate calls to for each.
!         !%-----------------------------------------------------------------
!             real(8), intent(inout) :: outdata(:)  !% element data for output
!             integer, intent(in)    :: JMidx       !% index of JM junction
!             integer, intent(in)    :: fiIdx       !%  index of map up or down to face
!             integer, intent(in)    :: frCol    !%  column in faceR array
!             integer, intent(in)    :: kstart       !% = 1 for upstream branches, 2 for down 
!             real(8) :: k1,k2   
!         !%-----------------------------------------------------------------
    
!         k1 = JMidx + kstart
!         k2 = JMidx + max_branch_per_node

!         where (elemSI(k1:k2:2,esi_JB_Exists) .eq. oneI)
!               outdata(k1:k2:2) = faceR(elemI(k1:k2:2,fiIdx),frCol) 
!         endwhere

!     end subroutine lljunction_branch_getface
! !%
! !%==========================================================================
!%==========================================================================
!% 
    real(8) pure function lljunction_branch_Qnet (JMidx,idir)    
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the net inflow (idir = 1) or outflow (idir=-1) from
        !% junction JMidx for branch flowrates (only)
        !%------------------------------------------------------------------
            integer, intent(in) :: JMidx,idir
            real(8), dimension(max_branch_per_node) :: Qdir
            integer :: kk
        !%------------------------------------------------------------------

        Qdir = zeroR
        do concurrent (kk=1:max_branch_per_node)
            if ((abs(elemR(JMidx+kk,er_Flowrate))               > zeroR) .and. &
                (elemSI(JMidx+kk,esi_JB_Exists)     == oneI) .and. &
                (elemSI(JMidx+kk,esi_JB_CanModifyQ) == oneI)         ) then

                Qdir(kk) =  branchsign(kk) *  elemR(JMidx+kk,er_Flowrate)      &
                    *onehalfR * (                                              &
                                    oneR + real(idir,8) * branchsign(kk)       &
                                         * elemR(JMidx+kk,er_Flowrate)         &
                                         / abs(elemR(JMidx+kk,er_Flowrate))    &
                                )
            end if
        end do

        lljunction_branch_Qnet = sum(Qdir)
        
    end function lljunction_branch_Qnet
!%    
!%==========================================================================
!%==========================================================================
!%     
    subroutine lljunction_branch_update_DeltaQ (JMidx, dH)   
        !%------------------------------------------------------------------
        !% Description:
        !% Updates the delta Q for a JB based on the dQ/dH
        !%------------------------------------------------------------------
            integer, intent(in) :: JMidx
            real(8), intent(in) :: dH
            integer :: ii
        !%------------------------------------------------------------------

        do ii=1,max_branch_per_node
            if (elemSI(JMidx+ii,esi_JB_Exists) .ne. oneI) cycle   
            
            elemR(JMidx+ii,er_DeltaQ) = &
                elemSR(JMidx+ii,esr_JB_dQdH) * dH 

        end do

    end subroutine lljunction_branch_update_DeltaQ
!%    
!%==========================================================================
!%==========================================================================
!% 
    subroutine lljunction_branch_update_flowrate (JMidx)   
        !%------------------------------------------------------------------
        !% Description:
        !% Updates the Qfor a JB based on the dQ/dH
        !%------------------------------------------------------------------
            integer, intent(in) :: JMidx
            integer             :: ii
        !%------------------------------------------------------------------

        do ii=1,max_branch_per_node
            if (elemSI(JMidx+ii,esi_JB_Exists) .ne. oneI) cycle
            elemR(JMidx+ii,er_Flowrate) &
                = elemR (JMidx+ii,er_Flowrate) + elemR(JMidx+ii,er_DeltaQ )
        end do

   end subroutine lljunction_branch_update_flowrate
!%    
!%==========================================================================
!%==========================================================================
!%     
    subroutine lljunction_CC_for_JBadjacent (thisColP, istep, isUpstreamYN)
        !%-----------------------------------------------------------------
        !% Description
        !% Adusts the values on the CC elements that are JB adjacent for the
        !% solution of the JM and JB elements. Note that this is called for
        !% the CC elements (thisColP). For parallelization, the DeltaQ from
        !% JB elements on another image must already be transferred to the
        !% intervening face and synced.
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisColP, istep
            logical, intent(in) :: isUpstreamYN

            integer, pointer    :: thisCC(:), Npack, fIdx(:)
            real(8), pointer    :: dt, crk, oldVolume(:)
            real(8) :: bsign

            integer :: mm
        !%-----------------------------------------------------------------   
        !% Aliases
            Npack       => npack_elemP(thisColP)
            if (Npack < 1) return
            thisCC      => elemP(1:Npack,thisColP)
            dt          => setting%Time%Hydraulics%Dt
            crk         => setting%Solver%crk2(istep)
            oldVolume   => elemR(:,er_Temp01)
        !%-----------------------------------------------------------------  

            print *, 'obsolete?  20230927'
            stop 2098734
            
        if (isUpstreamYN) then 
            !% --- pointer to the downstream JB face for CC upstream of JB
            fIdx => elemI(:,ei_Mface_dL)
            bsign = +oneR
        else
            !% --- pointer to the upstrea JB face for CC downstream of JB
            fIdx => elemI(:,ei_Mface_uL)
            bsign = -oneR
        end if    

        !% -- adjust the conservative flux on JBadjacent CC for deltaQ
        !if (istep == 2) then 
        ! faceR(fIdx(thisCC),fr_Flowrate_Conservative) &
        !         = faceR(fIdx(thisCC),fr_Flowrate_Conservative) &
        !         + faceR(fIdx(thisCC),fr_DeltaQ)
            

        ! !% PROBLEM brh 20230927 IS THIS NEEDED? 
        ! !% IF THIS IS SHARED, THE JB SHOULD ALREADY BE CHANGED, AND
        ! !% ABVOE WE ARE SETTING THE FACE ON THE ADJACENT PROCESSOR TO MATCH        
        ! !% --- check if the elem data has either been pushed 
        ! !%     to a shared face. if so then mark that face
        ! where (faceYN(fIdx(thisCC),fYN_isSharedFace))
        !     faceYN(fIdx(thisCC),fYN_isSharedFaceDiverged) = .true.
        ! end where
        !end if


        !% --- store present volume
        oldVolume(thisCC) = elemR(thisCC,er_Volume)

        !% --- changing the face flowrate changes the volume.
        !%     On an upstream element, a negative DeltaQ causes
        !%     an increase in the upstream element volume
        ! elemR(thisCC,er_Volume) = elemR(thisCC,er_Volume) &
        !     - bsign * crk * dt *  faceR(fIdx(thisCC),fr_DeltaQ)    

        elemR(thisCC,er_Volume) = elemR(thisCC,er_Volume) &
            - bsign * dt *  faceR(fIdx(thisCC),fr_DeltaQ)  

        !% --- changing the face flowrate and element volume changes the element velocity
        where (elemR(thisCC,er_Volume) > setting%ZeroValue%Volume)   
            elemR(thisCC,er_Velocity) &
                = (elemR(thisCC,er_Flowrate) * elemR(thisCC,er_Length) &
                   + faceR(fIdx(thisCC),fr_DeltaQ) * elemR(thisCC,er_Length) ) &
                / elemR(thisCC,er_Volume)
        elsewhere
            elemR(thisCC,er_Velocity) = zeroR 
        endwhere
  
        !% --- update the auxiliary variables
        !%     as we do not have separate packed open and closed for JB adjacent
        !%     we do this through the call with isSingularYN = .true. and the JM
        !%     index provided. Although this is inefficient, there should not be
        !%     that many elements in the upstream/downstream of junction sets.
        do mm=1,Npack
            call update_auxiliary_variables_CC (dummyIdx, dummyIdx, dummyIdx, &
                .false., .true., thisCC(mm))
        end do

    end subroutine lljunction_CC_for_JBadjacent
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine lljunction_conservation_fix (JMidx, resid, QnetIn, QnetOut)
        !%------------------------------------------------------------------
        !% Description:
        !% Ad hoc fix to junction flow rates to ensure mass conservation
        !% Typically required because dH adjustment was limited
        !% Modifies elemR(:,er_Flowrate) and Qoverflow
        !%------------------------------------------------------------------
            integer, intent(in)    :: JMidx
            real(8), intent(inout) :: resid,  QnetOut
            real(8), intent(inout) :: QnetIn

            real(8), pointer :: Qoverflow, Qstorage, dQdH(:)
            real(8), pointer :: OverflowDepth, PondedDepth
            real(8), pointer :: pi

            integer, pointer :: fup(:), fdn(:)

            integer :: kk
            real(8) :: dQoverflow, MinHeadForOverflow
            !%, OverflowDepth, PondedHead

            real(8), dimension(max_branch_per_node) :: dQ, dH, areaQ
            real(8), parameter :: localEpsilon = 1.0d-6
            real(8) :: Aout, Ain, AoutOverflow, AinPonded, Astorage
            integer :: bcount
            logical, dimension(max_branch_per_node) :: bFixYN

            logical :: repeatYN
        !%------------------------------------------------------------------
        !% Alias:
            Qoverflow     => elemSR(JMidx,esr_JM_OverflowPondingRate)
            Qstorage      => elemSR(JMidx,esr_JM_StorageRate)
            OverflowDepth => elemSR(JMidx,esr_JM_OverflowDepth)
            PondedDepth   => elemSR(JMidx,esr_JM_ExternalPondedDepth)
            dQdH          => elemSR(:,esr_JB_dQdH)

            fup       => elemI(:,ei_Mface_ul)
            fdn       => elemI(:,ei_Mface_dl)
            pi        => setting%Constant%pi
        !%------------------------------------------------------------------

        repeatYN = .true.

        dQ = zeroR
        dH = zeroR
        dQoverflow = zeroR

        Aout   = zeroR
        Ain    = zeroR

        !% --- note that by definition: QnetIn > 0 and QnetOut < 0 
        !%     resid > 0 implies too much inflow
        !%     resid < 0 implies too much outflow

        ! !% --- attempt to adjust residual with overflow/ponding alone
        ! if ((resid < zeroR) .and. (Qoverflow < zeroR)) then 
        !     !% --- reduce negative magnitude or eliminate negative Qoverflow rate
        !     !% --- note Qoverflow < 0 is an outflow
        !     if (Qoverflow < resid) then 
        !         !% --- resid can be fully accounted for by overflow reduction
        !         Qoverflow = Qoverflow - resid
        !         resid = zeroR
        !         return  !% --- no further residual processing needed
        !     else
        !         !% --- resid removes the overflow, but some discrepancy remains
        !         QnetOut   = QnetOut - Qoverflow
        !         resid     = resid - Qoverflow  
        !         Qoverflow = zeroR
        !         !% -- continue with residual processing
        !     end if
        ! elseif ((resid > zeroR) .and. (Qoverflow > zeroR)) then 
        !     !% --- reduce positive magnitude or elimenate positive Qoverflow rate
        !     !%     this is an inflow that only occurs due to ponding
        !     if (Qoverflow > zeroR) then 
        !         !% --- resid can be fully accounted for by overflowr eduction
        !         Qoverflow = Qoverflow - resid
        !         resid = zeroR
        !         return !% --- no further residual processing needed
        !     else
        !         !% --- resid reduces the overflow, but some resid remains
        !         QnetIn    = QnetIn - Qoverflow
        !         resid     = resid - Qoverflow 
        !         Qoverflow = zeroR
        !     end if
        ! else
        !     !% --- if Qoverflow == 0 no action required
        ! endif
            
        
        !% --- get flow areas of overflows or ponding
        if (Qoverflow > zeroR) then 
            !% --- inflow from ponding only (use external ponded depth)
            select case (elemSI(JMidx,esi_JM_OverflowType))
                case (PondedWeir)
                    AinPonded =  PondedDepth * twoR * sqrt( pi * elemSR(JMidx,esr_Storage_Plan_Area))
                case (PondedOrifice)
                    AinPonded =  PondedDepth * elemSR(JMidx,esr_JM_OverflowOrifice_Length)
                case default 
                    print *, 'CODE ERROR unexpected case default: Overflow with No overflow Type'
                    print *, 'JMidx ',JMIdx, ' ',trim(node%Names(elemI(JMidx,ei_node_Gidx_Bipquick))%str)
                    print *, 'Qoverflow ', Qoverflow 
                    print *, 'Head      ', elemR(JMidx,er_Head)
                    call util_crashpoint(77987233)
            end select
            AoutOverflow = zeroR

        elseif (Qoverflow < zeroR) then
            !% --- outflow either as overflow or to ponding (OverflowDepth > 0)
            select case (elemSI(JMidx,esi_JM_OverflowType))
                case (OverflowWeir, PondedWeir)
                    AoutOverflow = OverflowDepth * twoR * sqrt( pi * elemSR(JMidx,esr_Storage_Plan_Area))
                case (OverflowOrifice, PondedOrifice)
                    AoutOverflow = OverflowDepth * elemSR(JMidx,esr_JM_OverflowOrifice_Length)
                case default 
                    print *, 'CODE ERROR unexpected case default'
                    call util_crashpoint(7722366)
            end select
            AinPonded = zeroR

        else 
            AinPonded    = zeroR 
            AoutOverflow = zeroR
        end if

        !% --- Compute contributing areaQ(kk) in branches
        bFixYN = .false. !% --- whether or not a branch can be fixed for conservation
        areaQ  = zeroR   !% --- area for each branch
        bcount = zeroI   !% --- number of modifiable branches
        do kk=1,max_branch_per_node
            !% --- cycle if this branch cannot contribute to flow
            if ((elemSI(JMidx+kk,esi_JB_Exists)        .ne. oneI) .or.  &
                (elemSI(JMidx+kk,esi_JB_CanModifyQ)    .ne. oneI)        ) cycle

            if (elemSI(JMidx+kk,esi_JB_IsUpstream) == oneI) then 
                !% --- upstream branch
                !%     cycle if low junction head or high Fr in branch cannot be adjusted for Q
                if ((elemR(JMidx   ,er_Head) .le. faceR(fup(JMidx+kk),fr_Zbottom)) .or. &
                    (elemR(JMidx+kk,er_FroudeNumber) .ge. oneR)) cycle
                !% --- otherwise set bFix to modify this branch
                bFixYN(kk) = .true.  
                    bFixYN(kk) = .true.  
                bFixYN(kk) = .true.  
                bcount = bcount + oneI
                areaQ(kk) = max(faceR(fup(JMidx+kk),fr_Area_d),elemR(JMidx+kk,er_Area))

            else
                !% --- downstream branch
                !%     low junction head or high -Fr in branch cannot be adjusted for Q
                if ((elemR(JMidx  ,er_Head) .le. faceR(fdn(JMidx+kk),fr_Zbottom)) .or. &
                    (elemR(JMidx+kk,er_FroudeNumber) .le. -oneR)) cycle
                !% --- otherwise set bFix to modify this branch
                bFixYN(kk) = .true.
                bcount = bcount + oneI
                areaQ(kk) = max(faceR(fdn(JMidx+kk),fr_Area_u),elemR(JMidx+kk,er_Area))

            end if
            !% --- note, should not reach here unles bFix == .true.

            if (bcount < 1) then 
                print *, 'CODE ERROR bcount = 0; should not have reached this point '
                call util_crashpoint(6209873)
            end if

        end do

        !% --- if no branches with flow and no ponded inflow, then use all branches
        if ((bcount < 1) .and. (AinPonded == zeroR)) then 
            !% --- occurs during wetting drying, use all real branches for adjustment
            do kk=1,max_branch_per_node
                if ((elemSI(JMidx+kk,esi_JB_Exists)        .ne. oneI) .or.  &
                    (elemSI(JMidx+kk,esi_JB_CanModifyQ)    .ne. oneI)        ) cycle

                if (elemSI(JMidx+kk,esi_JB_IsUpstream) == oneI) then
                    !% --- upstream branch
                    bFixYN(kk) = .true.  
                        bFixYN(kk) = .true.  
                    bFixYN(kk) = .true.  
                    bcount = bcount + oneI
                    areaQ(kk) = max(faceR(fup(JMidx+kk),fr_Area_d),elemR(JMidx+kk,er_Area))

                else
                    !% --- downstream branch
                    bFixYN(kk) = .true.
                    bcount = bcount + oneI
                    areaQ(kk) = max(faceR(fdn(JMidx+kk),fr_Area_u),elemR(JMidx+kk,er_Area))

                end if
            end do

            if (bcount < 1) then 
                print *, 'CODE ERROR Mass residual in an unexpected condition'
                call util_crashpoint(2298522)
            endif
        else
            !% --- continue
        end if

        do while (repeatYN)

            !% --- Accumulate area weighting for inflows and outflows
            Ain  = AinPonded
            Aout = AoutOverflow
            do kk=1,max_branch_per_node
                if ((elemSI(JMidx+kk,esi_JB_Exists)        .ne. oneI) .or.  &
                    (elemSI(JMidx+kk,esi_JB_CanModifyQ)    .ne. oneI)        ) cycle
                !% --- accumulate the in and outflow areas for modifiable flows
                if ((real(branchsign(kk),8) * elemR(JMidx+kk,er_Flowrate)) > zeroR) then 
                    Ain  = Ain  + areaQ(kk)
                else
                    Aout = Aout + areaQ(kk)
                end if
            end do

            !% --- area associated with increased/decreased storage
            if (elemR(JMidx,er_Head) > elemR(JMidx,er_Zcrown)) then 
                Astorage = elemSR(JMidx,esr_JM_Present_PlanArea)
            else
                Astorage = elemSR(JMidx,esr_Storage_Plan_Area)
            end if

            if ((Ain < setting%ZeroValue%Area) .and. (Aout < setting%ZeroValue%Area)) then 
                !% degenerate condition 
                print *, 'CODE ERROR unexpected junction condition '
                print *, 'JMidx ',JMidx
                print *, Ain, Aout, setting%ZeroValue%Area
                call util_crashpoint(629784)
            else 
                !% --- distribute residual in proportion to area
                do kk=1,max_branch_per_node
                    if (.not. bFixYN(kk)) cycle 
                    !% --- note resid<0 is too much outflow and branchsign for
                    !%     downstream is -1, while the value of a downstream outflow
                    !%     is positive. We require dQ < 0 to decrease
                    !%     the positive outflow on a downstream branch. Similarly
                    !%     for an upstream branch the Q is inherently negative and the
                    !%     branch sign is positive so the dQ > 0 is needed to reduce
                    !%     the (negative) outflow on an upstream branch
                    dQ(kk) = - resid * real(branchsign(kk),8) * areaQ(kk) / (Ain + Aout + Astorage) 

                end do
                !% --- overflow rate adjustment
                if (Qoverflow < zeroR) then 
                    !% --- overflow/ponding is an outflow (negative Q)
                    !%     negative residual indicates too much outflow,
                    !%     which requires a positive dQ to reduce the
                    !%     outflow magnitude
                    dQoverflow =  - resid * AoutOverflow / (Ain + Aout + Astorage)

                    !%     For outflow, we limit the change in the outflow to zero
                    !%     i.e., we do not allow an outflow to become an inflow  

                    !% --- check to see if too much overflow was removed
                    !%     i.e., a positive residual (too much inflow)
                    !%     in which case reduce the residual set the overflow amd dQoverflow to zero
                    !%     and repeat without overflow
                    if (Qoverflow + dQoverflow > zeroR) then 
                        resid = resid - Qoverflow
                        Qoverflow    = zeroR
                        dQoverflow   = zeroR
                        AinPonded    = zeroR 
                        AoutOverflow = zeroR
                        repeatYN = .true.
                    else 
                        repeatYN = .false.
                    endif

                elseif (Qoverflow > zeroR) then 
                    !% --- Ponding with higher head outside is an inflow (positive Q)
                    !%     and a negative residual indicates too little inflow
                    !%     so dQ>0 is required to increase inflow rate.
                    !%     Conversely, a positive residual requires reduction of inflow
                    !%     so dQ < 0 is needed
                    dQoverflow = -resid * AinPonded / (Ain + Aout + Astorage)   
                    
                    !%     For inflow, we limit the change in the inflow value to zero
                    !%     i.e., we do not allow the inflow to become an outflow

                    !% --- check to see if too much ponding was removed
                    !%     i.e., a positive residual (too much inflow)
                    !%     in which case reduce the residual set the overflow amd dQoverflow to zero
                    !%     and repeat without overflow
                    if (Qoverflow + dQoverflow < zeroR) then 
                        resid = resid - Qoverflow
                        Qoverflow    = zeroR
                        dQoverflow   = zeroR
                        AinPonded    = zeroR 
                        AoutOverflow = zeroR
                        repeatYN = .true.
                    else 
                        repeatYN = .false.
                    endif

                else
                    !% Qoverflow == 0
                    dQoverflow = zeroR
                    repeatYN = .false.
                end if
            end if

        end do

        !% --- positive residual (too much inflow) causes increased storage rate
        Qstorage = Qstorage + resid *  Astorage / (Ain + Aout + Astorage)

        !% --- update flowrates
        do kk=1,max_branch_per_node
            if (elemSI(JMidx+kk,esi_JB_Exists) .ne. oneI) cycle
            if (bFixYN(kk) .ne. zeroI) then
                elemR(JMidx+kk,er_Flowrate) = elemR(JMidx+kk,er_Flowrate) + dQ(kk)
            end if
        end do

        !% --- update overflow
        Qoverflow = Qoverflow + dQoverflow

        !% --- recompute the residual
        resid = lljunction_conservation_residual (JMidx)

        if (abs(resid) > localEpsilon) then 
            print *, ' '
            print *, 'resid ',resid, ' at junction element ',JMidx

            print *, 'Node number ',elemI(JMidx,ei_node_Gidx_SWMM)
            if (elemI(Jmidx,ei_node_Gidx_SWMM) .ne. nullvalueI) then 
                print *, 'node name ',trim(node%Names(elemI(Jmidx,ei_node_Gidx_SWMM))%str)
            end if
            print *, 'on solution step ',setting%Time%Step
            print *, 'at time ',setting%Time%Now / (3600.d0), ' hours'
            print *, 'CODE ERROR unexpected flowrate residual'
            call util_crashpoint(6109873)
            return 
        end if

    end subroutine lljunction_conservation_fix
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8)  function lljunction_conservation_residual (JMidx) 
        !%------------------------------------------------------------------
        !% Description:
        !% computes flowrate residual, which is >0 for too much inflow 
        !% and < 0 for too much outflow
        !%------------------------------------------------------------------
            integer, intent(in) :: JMidx
            real(8) :: QnetBranches
        !%------------------------------------------------------------------

        QnetBranches = lljunction_main_sumBranches (JMidx,er_Flowrate, elemR)
        
        lljunction_conservation_residual = QnetBranches       &
             + elemSR(JMidx,esr_JM_OverflowPondingRate)    &
             - elemSR(JMidx,esr_JM_StorageRate)     &
             + elemR (JMidx,er_FlowrateLateral)

            !  if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then
            !     print *, ' '
            !     print *, 'resid components'
            !     print *, QnetBranches
            !     print *, elemSR(JMidx,esr_JM_OverflowPondingRate)
            !     print *, elemSR(JMidx,esr_JM_StorageRate)
            !     print *, elemR (JMidx,er_FlowrateLateral)
            !     print *, ' '

            !  end if     

    end function lljunction_conservation_residual
!%    
!%==========================================================================
!%==========================================================================
!% 
    subroutine lljunction_main_dHcompute  &
        (JMidx, dH, dQdHoverflow, dQdHstorage, Qnet, Hbound, istep, &
         isOverflow, isPonding, isCrossingIntoSurcharge, isCrossingOutofSurcharge)
        !%------------------------------------------------------------------
        !% Description
        !% Top level computation of dH in first step of junction solution
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in)    :: JMidx, istep
            logical, intent(in)    :: isOverflow, isPonding
            logical, intent(in)    :: isCrossingIntoSurcharge, isCrossingOutofSurcharge
            real(8), intent(inout) :: dH, dQdHoverflow, dQdHstorage
            real(8), intent(in)    :: Qnet
            real(8), dimension(2), intent(in) :: Hbound

            real(8), pointer :: Qstorage
            real(8)          :: dQdHbranches
            real(8)          :: divisor

            real(8), parameter :: localEpsilon = 1.0d-6
        !%------------------------------------------------------------------
        !% Aliases
            Qstorage    => elemSR(JMidx,esr_JM_StorageRate) !% positive is increasing storage
        !%------------------------------------------------------------------   
            
        !% --- if zero net flow, then storage and flowrates do not change.
        !%     Set storage rate to zero, volume to old volume, and
        !%     DeltaQ to zero. No need to do anything to flowrates
        !%     Then we're done with this junction
        if (Qnet == zeroR) then 
            Qstorage = zeroR
            dH = zeroR
            !% --- delta Q is NOT accumulative 
            
            return !% nothing more for this junction
        end if

        !% --- compute storage rate of change at the present head
        dQdHstorage = elemSR(JMidx,esr_JM_Present_PlanArea) / setting%Time%Hydraulics%Dt

        ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
        !     print *, 'dQdH storage, planarea ',dQdHstorage, elemSR(JMidx,esr_JM_Present_PlanArea)
        ! end if

        if (isOverflow .or. isPonding) then
            !% --- compute overflow rate of change with change in head
            dQdHoverflow = lljunction_main_dQdHoverflow (JMidx)
        else
            dQdHoverflow = zeroR
        end if

        ! if ((setting%Time%Step > 54165) .and. (JMidx == 109)) then
        !     print *, 'dQdHoverflow:    ',dQdHoverflow 
        ! end if

        !% --- compute net dQdH of branches
        dQdHbranches = lljunction_main_sumBranches(JMidx,esr_JB_dQdH, elemSR)

        !% --- divisor
        divisor =  dQdHstorage -  dQdHbranches - dQdHoverflow

        if (abs(divisor) > localEpsilon ) then 
            dH = Qnet / divisor
        else
            dH = zeroR
        end if

        ! if (printJM == JMidx) then 
        !     print *, ' '
        !     print *, 'dH here ',dH 
        !     print *, ' '
        ! end if

        ! !% --- limit dH
        ! if (dH < Hbound(1)) then 
        !     !% --- lower limit
        !     dH = Hbound(1)
        !     elemSI(JMidx,esi_JM_HeadLimit) = -oneI
        ! elseif (dH > Hbound(2)) then 
        !     dH = Hbound(2)
        !     elemSI(JMidx,esi_JM_HeadLimit) = +oneI
        ! else
        !     elemSI(JMidx,esi_JM_HeadLimit) = zeroI
        ! end if

    end subroutine lljunction_main_dHcompute
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine lljunction_main_dryingfix (JMidx, Qnet)
        !%------------------------------------------------------------------
        !% Description:
        !% Adjusts branch outflows and/or negative inflow rates when the
        !% net flowrate with provide a negative volume
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: JMidx
            real(8), intent(inout) :: Qnet

            real(8)          :: Qout, Qproportion, Vdelta
            real(8), pointer :: dt, Vzero
            logical          :: isLateraOutflow = .false.
            integer          :: kk
        !%------------------------------------------------------------------
        !% Aliases:
            dt    => setting%Time%Hydraulics%Dt
            Vzero => setting%ZeroValue%Volume
        !%------------------------------------------------------------------
        !% Preliminaries
            !% --- Compute change in volume minus a small amount
            !% --- note that this should only be called when Vdelta < Vzero
            !%     Qnet < 0 for net outflow and we will correct to 10*Vzero
            !%     Here -Vdelta is the volume increase required to maintain a
            !%     positive small volume when Volume is going negative
            Vdelta = elemR(JMidx,er_Volume_N0) + Qnet*dt - tenR*Vzero
            if (Vdelta .ge. Vzero) return  !% No adjustment required as Volume is increasing
        !%------------------------------------------------------------------   

        !% --- get the flowrates of all branch outflows (returns Qout < 0)   
        Qout = lljunction_main_sumBranches_InOrOutFlow &
                (JMidx,er_Flowrate,elemR,.false.,.true.)

        !% --- add Qlateral outflow to the net outflow
        if (elemR(Jmidx,er_FlowrateLateral) < zeroR) then  
            isLateraOutflow = .true.
            Qout = Qout + elemR(Jmidx,er_FlowrateLateral)
        end if     

        if (-Qout * dt < -Vdelta) then 
            !% --- If the magnitude of the outflow volume is smaller
            !%     than the volume change (-Vdelta) needed to bring
            !%     the volume back to a small positive value, then
            !%     reducing the outflow magnitude cannot bring the
            !%     volume back to positive, so we simply
            !%     eliminate all the outflows. The remaining negative
            !%     volume will need to be handled elsewhere      
            Qproportion = zeroR
            if (elemR(Jmidx,er_FlowrateLateral) < zeroR) then
                elemR(Jmidx,er_FlowrateLateral) = zeroR  
            end if
        else 
            !% --- the proportional fix required
            !%     Note: is > 0 as Vdelta < 0 and Qout < 0
            !%     expect abs(Vdelta) < abs(Qout * dt)
            !%     that is, the allowable volume change is less than out flow
            Qproportion = Vdelta / (Qout * dt)  
        end if
        
        do kk=1,max_branch_per_node
            !% --- ignore dummy branches
            if (elemSI(JMidx+kk,esi_JB_Exists) .ne. oneI) cycle
            !% --- ignore inflows
            if (branchsign(kk) * elemR(JMidx+kk,er_Flowrate) .ge. zeroR) cycle
            !% --- remove old outflow from Qnet
            Qnet = Qnet - elemR(JMidx+kk,er_Flowrate)
            !% --- reduce outflows
            elemR(JMidx+kk,er_Flowrate) =  elemR(JMidx+kk,er_Flowrate) * Qproportion
            !% --- add new flow to Qnet
            Qnet = Qnet + elemR(JMidx+kk,er_Flowrate)
        end do

        !% --- adjust lateral outflow
        if (isLateraOutflow) then 
            !% --- remove old lateral outflow
            Qnet = Qnet - elemR(JMidx,er_FlowrateLateral)
            !% --- adjust lateral outflow
            elemR(JMidx,er_FlowrateLateral) = elemR(JMidx,er_FlowrateLateral) * Qproportion
            !% --- add new lateral outflow
            Qnet = Qnet + elemR(JMidx,er_FlowrateLateral)
        end if
        
    


    end subroutine lljunction_main_dryingfix    
!%==========================================================================
!%==========================================================================
!% 
!% 
    !     subroutine lljunction_main_iscrossing_overflow_or_ponding ( JMidx,  dH, &
    !         isCrossingIntoOverflowOrPonding, isCrossingOutofOverflowOrPonding)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% updates depth, head, and deltaQ for a given dH on junction JMidx
    !         !%------------------------------------------------------------------
    !         !% Declarations:
    !             integer, intent(in)    :: JMidx
    !             real(8), intent(in)    :: dH
    !             logical, intent(inout) :: isCrossingIntoOverflowOrPonding
    !             logical, intent(inout) :: isCrossingOutofOverflowOrPonding
    !             real(8), pointer       :: MinHeadForOverflowPonding

    !         !%------------------------------------------------------------------
    !         !% Aliases
    !             MinHeadForOverflowPonding => elemSR(JMidx,esr_JM_MinHeadForOverflowPonding)
    !         !%------------------------------------------------------------------
            
    !         !% --- crossing from free surface to overflow/ponding
    !         if ((elemR(JMidx,er_Head)      .le. MinHeadForOverflowPonding)        & 
    !             .and.                                                    &
    !             ((elemR(JMidx,er_Head) + dH)  > MinHeadForOverflowPonding)      &
    !         ) then
    !             isCrossingIntoOverflowOrPonding  = .true.
    !             isCrossingOutofOverflowOrPonding = .false.
    !         else
    !             isCrossingIntoOverflowOrPonding = .false.
    !         end if

    !         !% --- crossing from overflow/ponding to free surface
    !         if ((elemR(JMidx,er_Head)        > MinHeadForOverflowPonding)        & 
    !             .and.                                                   &
    !             ((elemR(JMidx,er_Head) + dH) < MinHeadForOverflowPonding)      &
    !         ) then
    !             isCrossingOutofOverflowOrPonding  = .true.
    !             isCrossingIntoOverflowOrPonding   = .false.
    !         else
    !             isCrossingOutofOverflowOrPonding = .false.
    !         end if
    

    !     end subroutine lljunction_main_iscrossing_overflow_or_ponding    
    ! !%    
    ! !%==========================================================================
    ! !%==========================================================================
    ! !% 
    !     subroutine lljunction_main_iscrossing_surcharge ( JMidx,  dH, &
    !         isCrossingIntoSurcharge, isCrossingOutofSurcharge)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% updates depth, head, and deltaQ for a given dH on junction JMidx
    !         !%------------------------------------------------------------------
    !         !% Declarations:
    !             integer, intent(in)    :: JMidx
    !             real(8), intent(in)    :: dH
    !             logical, intent(inout) :: isCrossingIntoSurcharge
    !             logical, intent(inout) :: isCrossingOutofSurcharge
    !             real(8), pointer       :: MinHeadForOverflowPonding

    !         !%------------------------------------------------------------------
    !         !% Aliases
    !             MinHeadForOverflowPonding => elemSR(JMidx,esr_JM_MinHeadForOverflowPonding)
    !         !%------------------------------------------------------------------

    !         !% --- if head starts below or at crown and rises above the crown
    !         if ((elemR(JMidx,er_Head)      .le. elemR(JMidx,er_Zcrown))    &
    !             .and.                                                      &
    !             ((elemR(Jmidx,er_Head) + dH) >  elemR(Jmidx,er_Zcrown))    &
    !         ) then
    !             isCrossingIntoSurcharge  = .true.
    !             isCrossingOutofSurcharge = .false.
    !         else  
    !             isCrossingIntoSurcharge = .false.
    !         endif
    !         !% --- if head is above the crown, and drops below the crown
    !         if ((elemR(JMidx,er_Head)        > elemR(JMidx,er_Zcrown))  &
    !             .and.                                                   &
    !             ((elemR(Jmidx,er_Head) + dH) < elemR(Jmidx,er_Zcrown))  &
    !         ) then
    !                 isCrossingOutofSurcharge = .true.
    !                 isCrossingIntoSurcharge  = .false.
    !         else
    !             isCrossingOutofSurcharge = .false.
    !         end if

    !     end subroutine lljunction_main_iscrossing_surcharge
! !%    
!%==========================================================================
!%==========================================================================
!% 
    subroutine lljunction_main_netFlowrate &
        (JMidx, Qnet, canOverflowOrPond, isOverflow, isPonding)
        !%------------------------------------------------------------------
        !% Description
        !% Computes net flowrate in (positive) or out (negative) to
        !% a junction from all sources
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: JMidx
            real(8), intent(inout) :: Qnet
            logical, intent(in)    :: canOverflowOrPond
            logical, intent(inout) :: isOverflow, isPonding

            real(8), pointer       :: Qoverflow, Qlateral

            real(8)                :: QnetBranches
        !%------------------------------------------------------------------
        !% Aliases
            !% --- note that Qoverflow includes ponding rate
            Qoverflow   => elemSR(JMidx,esr_JM_OverflowPondingRate) !% negative is outflow
            Qlateral    => elemR (JMidx,er_FlowrateLateral) !% negative is outflow)
        !%------------------------------------------------------------------

        !% --- compute net flowrate from branches (both CC and Diag)
        QnetBranches = lljunction_main_sumBranches (JMidx,er_Flowrate, elemR)

        !% --- compute overflow/ponding rate (negative is outflow)
        !%     returns values for isOverflow and isPonding if present head
        !%     is overflow or ponding
        if (canOverflowOrPond) then
            Qoverflow = lljunction_main_Qoverflow &
                (JMidx,oneI,isOverflow,isPonding)
        else
            Qoverflow = zeroR
        end if

        !% --- net flowrate (Qnet > 0 is net inflow)
        Qnet = QnetBranches + Qoverflow + Qlateral 

        ! if (printJM == JMidx) then 
        !     print *, ' '
        !     print *, 'Qnet here ', Qnet 
        !     print *, ' '
        ! end if
        
    end subroutine lljunction_main_netFlowrate
!%    
!%==========================================================================
!%==========================================================================
!% 
    subroutine lljunction_main_update_intermediate &
        (JMidx, istep, dH, dQdHoverflow, dQdHstorage, MinHeadForOverflow, &
         isOverflow, isPonding, &
         isCrossingIntoOverflowOrPonding, isCrossingOutofOverflowOrPonding, &
         isCrossingIntoSurcharge, isCrossingOutofSurcharge)
        !%------------------------------------------------------------------
        !% Description
        !% updates depth, head, and deltaQ for a given dH on junction JMidx
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: JMidx, istep
            real(8), intent(in) :: dH, dQdHstorage, dQdHOverflow, MinHeadForOverflow
            logical, intent(in) :: isPonding, isOverflow
            logical, intent(in) :: isCrossingIntoOverflowOrPonding, isCrossingOutofOverflowOrPonding
            logical, intent(in) :: isCrossingIntoSurcharge, isCrossingOutofSurcharge

            real(8), pointer    :: Qstorage(:), Qoverflow(:)
            real(8)             :: QnetBranches
        !%------------------------------------------------------------------
        !% Aliases
            Qstorage    => elemSR(:,esr_JM_StorageRate)
            !% --- note that Qoverflow includes ponding rate
            Qoverflow   => elemSR(:,esr_JM_OverflowPondingRate) !% negative is outflow
        !%------------------------------------------------------------------

        !% --- update JM head and depth
        elemR(JMidx,er_Head)  = elemR(JMidx,er_Head)  + dH
        !% 20240209brh revised for compatibility with air trapping
        !elemR(JMidx,er_Depth) = min(elemR(JMidx,er_Head) - elemR(JMidx,er_Zbottom), elemR(JMidx,er_FullDepth))
        elemR(JMidx,er_Depth) = min(elemR(JMidx,er_Depth) + dH , elemR(JMidx,er_FullDepth))

        elemR(JMidx,er_Depth) = max(elemR(JMidx,er_Depth),0.99d0*setting%ZeroValue%Depth)
        elemR(JMidx,er_EllDepth) =  elemR(JMidx,er_Depth)

        if (isPonding .or. isOverflow) then
            !% --- update junction main overflow rate
            !%     We want Q = 0.5 (Q^n + Q^{n+1}) where
            !%     Q^{n+1} = Q^n + dH * dQdH 
            Qoverflow(JMidx) = Qoverflow(JMidx) + onehalfR * dH * dQdHoverflow
        end if

        ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then
        !     print *, 'dQdH ',elemSR(printJB,esr_JB_dQdH), dH
        ! end if

        !% --- compute JB element DeltaQ using dQdH * dH
        elemR((JMidx+1):(JMidx+max_branch_per_node), er_DeltaQ) = zeroR
        call lljunction_branch_update_DeltaQ (JMidx,dH)  

        ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then
        !     print *, 'DeltaQ ',elemR(printJB,er_DeltaQ)
        ! end if

        !% --- update the flowrates using DeltaQ
        call lljunction_branch_update_flowrate (JMidx) 

        !% --- update net Q branches (included CC and Diag)
        QnetBranches = lljunction_main_sumBranches (JMidx,er_Flowrate,elemR)

        !% --- update junction main storage flow rate
        Qstorage(JMidx) = dQdHstorage * dH   
         
        elemR(JMidx,er_Volume) =  lljunction_main_volume_from_storageRate (JMidx,istep)

        ! !% --- overwrite for threshold crossing for exact values
        ! if (isCrossingIntoSurcharge .or. isCrossingOutofSurcharge )then 
        !     elemR(JMidx,er_Volume)   = elemR(JMidx,er_FullVolume)
        !     elemR(JMidx,er_Depth)    = elemR(JMidx,er_FullDepth)
        !     elemR(JMidx,er_EllDepth) = elemR(Jmidx,er_FullDepth)
        !     elemR(JMidx,er_Head)     = elemR(JMidx,er_FullDepth) + elemR(JMidx,er_Zbottom)

        ! elseif (isCrossingIntoOverflowOrPonding .or. isCrossingOutofOverflowOrPonding) then
        !     !% PROBLEM SETTING VOLUME NEGLECTS SURCHARGE VOLUME IN SLOT THAT MAY EXIST!
        !     elemR(JMidx,er_Volume)   = elemR(JMidx,er_FullVolume) 
        !     elemR(JMidx,er_Depth)    = MinHeadForOverFlow - elemR(JMidx,er_Zbottom)
        !     elemR(JMidx,er_EllDepth) = elemR(JMidx,er_Depth)
        !     elemR(JMidx,er_Head)     = MinHeadForOverflow
        ! else 
        !     !% no action 
        ! end if

    end subroutine lljunction_main_update_intermediate
!%    
!%==========================================================================
!%==========================================================================
!% 
    subroutine lljunction_main_update_final &
        (JMidx, istep)
        !%------------------------------------------------------------------
        !% Description
        !% Top level computation of dH in first step of junction solution
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx, istep
            real(8), pointer    :: Qoverflow(:)
            real(8), pointer    :: dt, crk(:)
        !%------------------------------------------------------------------
        !% Aliases
            !% --- note that Qoverflow includes ponding rate
            Qoverflow   => elemSR(:,esr_JM_OverflowPondingRate) !% negative is outflow
            dt          => setting%Time%Hydraulics%Dt
            crk         => setting%Solver%crk2
        !%------------------------------------------------------------------   
 
        !% --- update the overflow volume based on rate and time step
        select case (elemSI(JMidx,esi_JM_OverflowType))
            case (OverflowWeir,OverflowOrifice)
                !% --- volume overflow (out is negative flowrate, gives positive overflow volume)
                !elemR(JMidx,er_VolumeOverflow) = -Qoverflow(JMidx)  * dt * crk(istep)
                elemR(JMidx,er_VolumeOverflow) = -Qoverflow(JMidx)  * dt 
            case (PondedWeir,PondedOrifice)
                !elemR(JMidx,er_VolumePonded)   = -Qoverflow(JMidx)  * dt * crk(istep)
                elemR(JMidx,er_VolumePonded)   = -Qoverflow(JMidx)  * dt 
            case (NoOverflow)
                !% no action
            case default
                print *, elemSI(JMidx,esi_JM_OverflowType), trim(reverseKey(elemSI(JMidx,esi_JM_OverflowType)))
                print *, 'CODE ERROR unexpected case default'
                call util_crashpoint(397894)
        end select

    end subroutine lljunction_main_update_final   
!%    
!%==========================================================================
!%==========================================================================
!% 
    real(8) function lljunction_main_dQdHoverflow (JMidx)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the overflow rate for  junction
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx
            real(8), pointer    :: Lorifice, OverflowDepth 
            real(8), pointer    :: coef2, coef4
            real(8), pointer    :: MinHeadForOverflowPonding,  ExternalPondedHeadDiff 
            real(8), pointer    :: Storage_Plan_Area
        !%------------------------------------------------------------------  
        !% Aliases
            Lorifice                  => elemSR(JMidx,esr_JM_OverflowOrifice_Length)
            OverflowDepth             => elemSR(JMidx,esr_JM_OverflowDepth)
            !ExternalPondedDepth       => elemSR(JMidx,esr_JM_ExternalPondedDepth)
            MinHeadForOverflowPonding => elemSR(JMidx,esr_JM_MinHeadForOverflowPonding)
            !ExternalPondedHead        => elemSR(JMidx,esr_JM_ExternalPondedDepth)
            ExternalPondedHeadDiff    => elemSR(JMidx,esr_JM_ExternalPondedHeadDiff)
            Storage_Plan_Area         => elemSR(JMidx,esr_JM_Present_PlanArea)       
            coef2    => setting%Junction%Overflow%coef2
            coef4    => setting%Junction%Overflow%coef4
        !%------------------------------------------------------------------  

        !% --- Head below surcharge, only possibility is a ponding inflow
        !%     But changing Q does not change the inflow
        if (elemR(JMidx,er_Head) .le. MinHeadForOverflowPonding) then 
            lljunction_main_dQdHoverflow = zeroR
            return 
        end if

        !% --- if possible overflow condition exists
        select case (elemSI(JMidx,esi_JM_OverflowType))  
            case (NoOverflow)
                lljunction_main_dQdHoverflow = zeroR 
                return 

            case (OverflowWeir)
                !% --- Using Brater and King Eq 5.10 dQ/dH = (3/2) C L H^{1/2}
                !%     Apply L = 2 (pi A)^{1/2} such that
                !%     coef2 = (3/2) (C sqrt(pi)) and
                !%     dQ/dH = coef2 * sqrt( A H )
                !%     minus sign as Q is outflow (negative) as H increases
                lljunction_main_dQdHoverflow = -coef2                            &
                    * sqrt(Storage_Plan_Area * OverflowDepth)
                return

            case (PondedWeir)
                !% --- see explanation for OverflowWeir, above
                !%     abs(OverflowDepth) required so that negative OverflowDepth
                !%     is allowed when ponding is an inflow
                lljunction_main_dQdHoverflow = -coef2                            &
                    * sqrt(Storage_Plan_Area * abs(ExternalPondedHeadDiff))
                return
                
            case (OverflowOrifice)
                !% --- use the supplied orifice length
                !%     Using Brater and King Eq. 4.16 or 4.17
                !%     dQ/dH = sqrt(2g) L sqrt(H)
                !%     coef4 = sqrt(2g)
                !%     minus sign as Q is outflow (negative) as H increases
                lljunction_main_dQdHoverflow = -coef4 * Lorifice * sqrt(OverflowDepth)   
                return
                        
            case (PondedOrifice)
                !% --- see explanation for OverflowOrifice, above
                !%     abs(OverflowDepth) required so that negative OverflowDepth
                !%     is allowed when ponding is an inflow
                lljunction_main_dQdHoverflow = -coef4 * Lorifice * sqrt(abs(OverflowDepth))   
                return

            case default
                !% --- should not reach here.
                print *, 'CODE ERROR unexpected case default'
                call util_crashpoint(1209874)

        end select 

    end function lljunction_main_dQdHoverflow
!%    
!%========================================================================== 
!%==========================================================================
!% 
    ! real(8) function lljunction_main_dQdHstorage &
    !         (JMidx,istep, isOverflow, isPonding, &
    !          isCrossingIntoSurcharge, isCrossingOutofSurcharge)  
    !     !%------------------------------------------------------------------
    !     !% Description
    !     !% Computes the storage rate of junction with tabular or functional
    !     !% storage
    !     !%------------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in) :: JMidx, istep
    !         logical, intent(in) :: isOverflow, isPonding
    !         logical, intent(in) :: isCrossingIntoSurcharge, isCrossingOutofSurcharge
    !         real(8)             :: planArea
    !     !%------------------------------------------------------------------  

    !     if (elemSI(JMidx,esi_JM_Type) .eq. NoStorage) then
    !         !% --- no storage rate for implied storage junctions is not finished
    !         print *, 'CODE ERROR NoStorage junction type is not supported'
    !         call util_crashpoint(11001093)
    !         lljunction_main_dQdHstorage = zeroR
    !         return
    !     endif
            
    !     lljunction_main_dQdHstorage = planArea / setting%Time%Hydraulics%Dt


    ! end function lljunction_main_dQdHstorage   
    !%    
!%========================================================================== 
!%==========================================================================
!%   
    subroutine lljunction_main_energyhead (thisColP)
        !%------------------------------------------------------------------
        !% Description
        !% computes total energy head at a junction main
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisColP
            integer, pointer    :: Npack, thisP(:)
        !%------------------------------------------------------------------
        !% Preliminaries
            Npack => npack_elemP(thisColP)
            if (Npack < 1) return
        !%-----------------------------------------------------------------
        !% Aliases
            thisP => elemP(1:Npack,thisColP)
        !%------------------------------------------------------------------

        elemR(thisP,er_EnergyHead) = elemR(thisP,er_Head) &
            + (elemR(thisP,er_Velocity)**2) / (twoR * setting%Constant%gravity)

    end subroutine lljunction_main_energyhead
!%    
!%==========================================================================
!%==========================================================================
!% 
    ! subroutine lljunction_main_head_bounds (JMidx, Hbound)
    !     !%-----------------------------------------------------------------
    !     !% Description
    !     !% Computes the minimum head, Hbound(1), and Hbound(2) 
    !     !%, maximum head for a junction
    !     !%     Hbound(2) is the maximum head in the surrounding elements
    !     !%     Hbound(1) is Zbottom of JM, or the lowest Z bottom of any 
    !     !%          branch if they are all higher than JM
    !     !% 20230912brh switched to using full energy head as bounds
    !     !%-----------------------------------------------------------------
    !     !% Declarations
    !         integer,               intent(in)    :: JMidx
    !         real(8), dimension(2), intent(inout) :: Hbound
    !         integer :: ii, JBidx
    !         integer, pointer :: fidx
    !         real(8), pointer :: grav, headJM, headAdj, EnergyHeadJM, EnergyHeadAdj
    !         real(8), dimension(max_branch_per_node,2) :: HbranchLimit
    !         real(8) :: lateralAdd
    !         logical :: jhead_lowlimit_TF
    !     !%-----------------------------------------------------------------
    !     !% Aliases
    !         grav => setting%Constant%gravity
    !     !%-----------------------------------------------------------------

    !     Hbound(1) = +huge(oneR)
    !     Hbound(2) = -huge(oneR)

    !     ! HbranchLimit(:,1) = +huge(oneR)
    !     ! HbranchLimit(:,2) = -huge(oneR)

    !     ! jhead_lowlimit_TF = .false.

    !     ! headJM       => elemR(JMidx,er_Head)
    !     ! EnergyHeadJM => elemR(JMidx,er_EnergyHead)

    !     ! !% FIND ALL BRANCHES WITH INFLOW HEAD THAT COULD LIMIT THE MAX, MIN
    !     ! !% JUNCTION HEAD
    !     ! !% --- cycle through branches (cannot be concurrent)
    !     ! !%     fadj* are faces for adjustment, zeroI is null value
    !     ! do ii=1,max_branch_per_node

    !     !     if (elemSI(JMidx+ii,esi_JB_Exists) .ne. oneI) cycle 
    !     !     !% --- diagnostic elements cannot be head limiters
    !     !     if (elemSI(JMidx+ii,esi_JB_Diag_adjacent) .eq. oneI) cycle
            
    !     !     if (mod(ii,2)== 0) then 
    !     !         !% --- downstream branch
    !     !         !% --- downstream face
    !     !         fidx => elemI(JMidx+ii,ei_Mface_dL)
    !     !         JBidx = JMidx+ ii
    !     !         headAdj       => faceR(fidx,fr_Head_Adjacent)
    !     !         EnergyHeadAdj => faceR(fidx,fr_EnergyHead_Adjacent)

    !     !         if (elemR(JBidx,er_Flowrate) > zeroR) then
    !     !             !% --- outflow on downstream branch 
    !     !             !%     provides only a lower limit
    !     !             if (headJM > headAdj) then
    !     !             ! if (EnergyHeadJM > EnergyHeadAdj) then
    !     !                 !% --- outflow with positive head gradient uses outside
    !     !                 !%     head as low limiter
    !     !                 HbranchLimit(ii,1) = headAdj

    !     !             ! elseif (headJM < headAdj) then
    !     !             ! ! elseif (EnergyHeadJM < EnergyHeadAdj) then
    !     !             !     !% --- outflow with inverse gradient cannot reduce below inside headJM
    !     !             !     HbranchLimit(ii,1) = headJM
    !     !             else 
    !     !                 !% --- zero gradient has no low limiter
    !     !             end if

    !     !         elseif (elemR(JBidx,er_Flowrate) < zeroR) then
    !     !             !% -- inflow on downstream branch
    !     !             !%    provides only an upper limit
    !     !             if (headAdj > headJM) then 
    !     !             ! if (EnergyHeadAdj > EnergyHeadJM) then 
    !     !                 !% --- inflow with positive head gradient uses outside
    !     !                 !%     head as high limiter
    !     !                 HbranchLimit(ii,2) = headAdj
    !     !                 !HbranchLimit(ii,2) = EnergyHeadAdj

    !     !             ! elseif (headAdj < headJM) then 
    !     !             ! ! elseif (EnergyHeadAdj < EnergyHeadJM) then     
    !     !             !     !% --- inflow with adverse head gradient cannot increase
    !     !             !     !%     JM head
    !     !             !     HbranchLimit(ii,2) = headJM 
    !     !             else 
    !     !                 !% --- zero gradient has no high limiter
    !     !             end if
    !     !         else
    !     !             !% no flowrate
    !     !         end if

    !     !     else
    !     !         !% --- upstream branch
    !     !         !% --- upstream face
    !     !         fidx => elemI(JMidx+ii,ei_Mface_uL)
    !     !         JBidx = JMidx+ ii
    !     !         headAdj => faceR(fidx,fr_Head_Adjacent)
    !     !         EnergyHeadAdj => faceR(fidx,fr_EnergyHead_Adjacent)

    !     !         if (elemR(JBidx,er_Flowrate) < zeroR) then
    !     !             !% --- outflow on upstream branch 
    !     !             !%     provides only a lower limit
    !     !             if (headJM > headAdj) then
    !     !             ! if (EnergyHeadJM > EnergyHeadAdj) then
    !     !                 !% --- outflow with positive head gradient uses outside
    !     !                 !%     head as low limiter
    !     !                 HbranchLimit(ii,1) = headAdj

    !     !             ! elseif (headJM < headAdj) then 
    !     !             ! ! elseif (EnergyHeadJM < EnergyHeadAdj) then 
    !     !             !     !% --- outflow with inverse gradient cannot reduce headJM
    !     !             !     HbranchLimit(ii,1) = headJM
    !     !             else 
    !     !                 !% --- zero gradient has no low limiter
    !     !             end if

    !     !         elseif (elemR(JBidx,er_Flowrate) > zeroR) then
    !     !             !% -- inflow on upstream branch
    !     !             !%    provides only an upper limit
    !     !             if (headAdj > headJM) then 
    !     !             ! if (EnergyHeadAdj > EnergyHeadJM) then 
    !     !                 !% --- inflow with positive head gradient uses outside
    !     !                 !%     head as high limiter
    !     !                 HbranchLimit(ii,2) = headAdj
    !     !                 !HbranchLimit(ii,2) = EnergyHeadAdj

    !     !             ! elseif (headAdj < headJM) then 
    !     !             ! ! elseif (EnergyHeadAdj < EnergyHeadJM) then 
    !     !             !     !% --- inflow with adverse head gradient cannot increase
    !     !             !     !%     JM head
    !     !             !     HbranchLimit(ii,2) = headJM
    !     !             else 
    !     !                 !% --- zero gradient has no high limiter
    !     !             end if

    !     !         else
    !     !             !% no flowrate
    !     !         end if
    !     !     end if
    !     ! end do

    !     ! !% --- select the lowest low limiter from all branches
    !     ! Hbound(1) = minval(HbranchLimit(:,1))

    !     ! !% --- select the highest limiter from all branches
    !     ! Hbound(2) = maxval(HbranchLimit(:,2))
        
    !     ! !% --- account for lateral inflow
    !     ! if (elemR(JMidx,er_FlowrateLateral) .ne. zeroR) then 
    !     !     if (elemSR(JMidx,esr_Storage_Plan_Area) > zeroR) then
    !     !         lateralAdd =  elemR(JMidx,er_FlowrateLateral) / elemSR(JMidx,esr_Storage_Plan_Area)
    !     !     else
    !     !         print *, 'Unexpected storage plan area of zero '
    !     !         print *, 'for junction element index ',JMidx
    !     !         print *, 'Node ',trim(node%Names(elemI(JMidx,ei_node_Gidx_SWMM))%str)
    !     !         call util_crashpoint(7109744)
    !     !     end if
    !     ! else
    !     !     lateralAdd = zeroR
    !     ! end if

    !     ! if ((elemR(JMidx,er_FlowrateLateral) < zeroR) .and. (Hbound(1) .ne. huge(oneR))) then 
    !     !     Hbound(1) = Hbound(1) + lateralAdd
    !     ! elseif ((elemR(JMidx,er_FlowrateLateral) > zeroR) .and. (Hbound(2) .ne. -huge(oneR))) then 
    !     !     Hbound(2) = Hbound(2) + lateralAdd
    !     ! end if

    !     if (Hbound(1) > 1000.d0) then
    !         Hbound(1) = -huge(oneR)
    !     end if
    !     if (Hbound(2) < -1000.d0) then 
    !         Hbound(2) = + huge(oneR)
    !     end if

    ! end subroutine lljunction_main_head_bounds  
!%    
!%========================================================================== 
!%==========================================================================
!%
    subroutine lljunction_main_overflow_conditions  (JMidx)
        !%-----------------------------------------------------------------
        !% Description
        !% Computes data needed for evaluating conditions for overflow
        !% and ponding
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx
            real(8), pointer    :: MinHeadForOverFlowPonding, OverflowDepth
            real(8), pointer    :: ExternalPondedArea, ExternalPondedDepth
            real(8), pointer    :: ExternalPondedHead, PondedHeadDiff
            real(8) :: WaterHead
        !%-----------------------------------------------------------------
        !% Aliases
            MinHeadForOverflowPonding => elemSR(JMidx,esr_JM_MinHeadForOverflowPonding)
            OverflowDepth             => elemSR(JMidx,esr_JM_OverflowDepth)
            ExternalPondedArea        => elemSR(JMidx,esr_JM_ExternalPondedArea)
            ExternalPondedDepth       => elemSR(JMidx,esr_JM_ExternalPondedDepth)
            ExternalPondedHead        => elemSR(JMidx,esr_JM_ExternalPondedHead)
            PondedHeadDiff            => elemSR(JMidx,esr_JM_ExternalPondedHeadDiff)
        !%-----------------------------------------------------------------

        if (setting%AirTracking%UseAirTrackingYN) then    
            WaterHead = elemR(JMidx,er_Head) - elemSR(JMidx,esr_JM_Air_HeadGauge)
        else 
            WaterHead = elemR(JMidx,er_Head)
        end if

        !% --- ponding junctions are treated differently than non-ponding
        if ( ExternalPondedArea > zeroR) then 

            !% --- used for conservation fix
            OverflowDepth = max(WaterHead - MinHeadForOverflowPonding, zeroR)    

            !% --- ponded head is the head available in the ponded area
            ExternalPondedDepth = elemR (JMidx,er_VolumePondedTotal)   &
                                / elemSR(JMidx,esr_JM_ExternalPondedArea)  

            ExternalPondedHead  = MinHeadForOverflowPonding  + ExternalPondedDepth 

            !% --- for ponding, the overflow depth is the difference                        
            !%     between the junction head and ponded head. A positive value
            !%     causes an outflow to ponding whereas a negative value
            !%     causes an inflow from ponding to the junciton
            PondedHeadDiff     = WaterHead - ExternalPondedHead

            !% --- limit negative ponded head difference to a waterfall condition
            !%     from the ponded depth
            if (PondedHeadDiff < zeroR) then
                PondedHeadDiff = max(PondedHeadDiff, -ExternalPondedDepth)
            else
                !% no action -- keep + PondedHeadDiff
            end if
        else 
            !% --- For non-ponding, the overflow depth is either positive or zero
            !%     This uses the present value of head
            OverflowDepth = max(WaterHead - MinHeadForOverflowPonding, zeroR)                   
        end if

        ! if ((setting%Time%Step > 54165) .and. (JMidx == 109)) then 
        !     print *, 'overflow depth: ',OverFlowDepth
        ! end if
        ! if (printJM == JMidx) then 
        !     print *, ' '
        !     print *, 'OverflowDepth, ExternalPondedDepth', OverflowDepth, ExternalPondedHead
        !     print *, ' '
        ! end if

    end subroutine lljunction_main_overflow_conditions
!%
!%========================================================================== 
!%==========================================================================
!% 
    real(8) function lljunction_main_plan_area (JMidx)
        !%------------------------------------------------------------------
        !% Description:
        !% sets the present plan area to either the storage area at the
        !% present head, or the latest surcharge length * width, or the
        !% special case of an overflow/ponding orifice that uses the 
        !% orifice dimensions
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) ::  JMidx
            real(8), pointer    :: SlotWidth, Length
            real(8)             :: WaterHead
        !%------------------------------------------------------------------
        !% Aliases
            SlotWidth => elemR(JMidx,er_SlotWidth)
            Length    => elemR(JMidx,er_Length)

        !%------------------------------------------------------------------

        if (setting%AirTracking%UseAirTrackingYN) then    
            WaterHead = elemR(JMidx,er_Head) - elemSR(JMidx,esr_JM_Air_HeadGauge)
        else 
            WaterHead = elemR(JMidx,er_Head)
        end if

        !% --- set the plan area used for storage
        if (elemYN(JMidx,eYN_canSurcharge)) then
            select case (elemSI(JMidx,esi_JM_OverflowType))
                case (NoOverflow)
                    if (elemR(JMidx,er_Head) < elemR(JMidx,er_Zcrown)) then
                        lljunction_main_plan_area = elemSR(JMidx,esr_Storage_Plan_Area) 
                        return
                    else    
                        !% --- preissmann slot surcharge plan area
                        !%     should be stored in present
                        lljunction_main_plan_area = elemSR(JMidx,esr_JM_Present_PlanArea)
                        return
                    end if

                case (PondedWeir,OverflowWeir)
                    if  (elemR(JMidx,er_Head) < elemR(JMidx,er_Zcrown)) then
                        lljunction_main_plan_area = elemSR(JMidx,esr_Storage_Plan_Area) 
                        return 
                    elseif ((elemR(JMidx,er_Head) .ge. elemR(JMidx,er_Zcrown)) .and. &
                            (           WaterHead  <   elemSR(JMidx,esr_JM_MinHeadForOverflowPonding)) ) then
                        !% --- preissmann slot surcharge plan area
                        lljunction_main_plan_area = elemSR(JMidx,esr_JM_Present_PlanArea)
                        return
                    else
                        !% --- overflow uses regular plan area
                        lljunction_main_plan_area = elemSR(JMidx,esr_Storage_Plan_Area) 
                        return
                    end if
                    

                case (PondedOrifice,OverflowOrifice)
                    if  (elemR(JMidx,er_Head) < elemR(JMidx,er_Zcrown)) then
                        lljunction_main_plan_area = elemSR(JMidx,esr_Storage_Plan_Area) 
                        return 
                    elseif ((elemR(JMidx,er_Head) .ge. elemR(JMidx,er_Zcrown)) .and. &
                            (           WaterHead  <   elemSR(JMidx,esr_JM_MinHeadForOverflowPonding)) ) then
                        !% --- preissmann slot surcharge plan area
                        lljunction_main_plan_area = elemSR(JMidx,esr_JM_Present_PlanArea)
                        return
                    else
                        !% --- overflow uses plan area similar to orifice
                        lljunction_main_plan_area = elemSR(JMidx,esr_JM_OverflowOrifice_Height) &
                                                  * elemSR(JMidx,esr_JM_OverflowDepth) 
                        return
                    end if

                case default
                    print *, 'CODE ERROR unexpected case default'
                    call util_crashpoint(22087445)
            end select
        else 
            lljunction_main_plan_area = elemSR(JMidx,esr_Storage_Plan_Area) 
        end if
    
    end function lljunction_main_plan_area
!%
!%==========================================================================
!%==========================================================================
!% 
    real(8) function lljunction_main_Qoverflow &
        (JMidx, istep, isOverflow, isPonding)  
        !%------------------------------------------------------------------
        !% Description
        !% Computes the overflow/ponding rate of a junction that is stored
        !% Negative is outflow, positive is inflow
        !% Weir overflow based on non-dimensional approach of Brater and King
        !% for broad-crested weir
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in)    :: JMidx, istep
            logical, intent(inout) :: isOverflow, isPonding
            real(8), pointer       :: Horifice, Lorifice, coef1, coef3, WeirFactor
            real(8), pointer       :: OverflowDepth, ExternalPondedDepth, Storage_Plan_Area
            real(8), pointer       :: VolumePondedTotal, MinHeadForOverflowPonding
            real(8), pointer       :: ExternalPondedArea, ExternalPondedHeadDiff
            real(8)                :: tempOverflow, WaterHead 

            real(8), pointer    :: dt, crk
        !%------------------------------------------------------------------  
        !% Aliases
            Horifice                  => elemSR(JMidx,esr_JM_OverflowOrifice_Height)
            Lorifice                  => elemSR(JMidx,esr_JM_OverflowOrifice_Length)
            MinHeadForOverFlowPonding => elemSR(JMidx,esr_JM_MinHeadForOverflowPonding)
            OverflowDepth             => elemSR(JMidx,esr_JM_OverflowDepth)
            ExternalPondedHeadDiff    => elemSR(JMidx,esr_JM_ExternalPondedHeadDiff)
            ExternalPondedArea        => elemSR(JMidx,esr_JM_ExternalPondedArea)
            ExternalPondedDepth       => elemSR(JMidx,esr_JM_ExternalPondedDepth)
            Storage_Plan_Area         => elemSR(JMidx,esr_Storage_Plan_Area)
            VolumePondedTotal         => elemR(Jmidx,er_VolumePondedTotal)
            WeirFactor => setting%Junction%Overflow%WeirLengthFactor
            coef1      => setting%Junction%Overflow%coef1
            coef3      => setting%Junction%Overflow%coef3
            dt         => setting%Time%Hydraulics%Dt
            crk        => setting%Solver%crk2(istep)
            
        !%------------------------------------------------------------------  

        !% --- Head below the surcharge, only possibility is a ponding inflow
        if (elemR(JMidx,er_Head) .le. MinHeadForOverflowPonding) then
            !% --- no outflowing overflow 
            isOverflow = .false.
            if (ExternalPondedArea == zeroR) then 
                isPonding = .false.
                !% --- no ponding inflow or outflow
                lljunction_main_Qoverflow = zeroR
                return

            else
                !% --- possible ponding inflow as waterfall into junction
                if (ExternalPondedDepth > zeroR) then 
                    isPonding  = .true.
                    isOverflow = .false.
                    !% --- ponding inflow as waterfall
                    select case (elemSI(JMidx,esi_JM_OverflowType))

                        case (PondedWeir)
                            !% --- see OverflowWeir below for explanation
                            !%     given + value as inflow
                            tempOverflow = coef1 *  WeirFactor * sqrt(Storage_Plan_Area) &
                                * ((abs(ExternalPondedHeadDiff))**threehalfR)
            
                        case (PondedOrifice)
                            !% ---- see OverflowOrifice below for explanation
                            !%      given + value as inflow
                            tempOverflow = coef3 * Lorifice &
                                    * ((abs(ExternalPondedHeadDiff))**threehalfR)

                        case default 
                            print *, 'CODE ERROR unexpected case default'
                            call util_crashpoint(7298723)
                    end select

                    !% --- limit inflow by ponded water available
                    !if ((tempOverflow * dt *crk)  > VolumePondedTotal) then
                    if ((tempOverflow * dt)  > VolumePondedTotal) then
                        !lljunction_main_Qoverflow = VolumePondedTotal / (dt * crk)
                        lljunction_main_Qoverflow = VolumePondedTotal / dt 
                    else 
                        lljunction_main_Qoverflow = tempOverflow
                    end if
                    return

                else
                    isPonding  = .false.
                    isOverflow = .false.
                    !% --- no ponding inflow or outflow
                    lljunction_main_Qoverflow = zeroR
                    return

                end if
            end if
        else 
            !% --- continue below for head above surcharge
            !%     Note that all of the conditional in the if..then above will result
            !%     in the subroutine returning. We only continue to the next step
            !%     if the junction head is above the minimum head for overflow,
            !%     which will require a overflow depth > 0 for a strict overflow
            !%     but may have a +/- overflow depth for ponding. Note that a negative
            !%     overflow depth indicates junction head is lower than ponded head
            !%     so the ponding causes and inflow
        end if

        !% --- possible inflow or outflow due to overflow or ponding
        !%     CONVENTIONS: 
        !%     (1) a flowrate is negative if it is an outflow
        !%      and positive if it is an inflow.
        !%     (2) a positive OverflowDepth causes an outflow (negative Q)
        !%      whereas a negative Overflow depth cannot occur
        !%     (3) a positive PondedHeadDiff causes and outflow (negative Q)
        !%      whereas a negative PondedHeadDiff causes an inflow (positive Q)
        !%     (4) A positive inflow is only allowed with ponding.
        !%     
        select case (elemSI(JMidx,esi_JM_OverflowType))

            case (NoOverflow)
                isPonding  = .false.
                isOverflow = .false.
                !% --- NoOverflow implies no ponding or overflow
                lljunction_main_Qoverflow = zeroR
                return

            case (OverflowWeir)
                isPonding  = .false.
                isOverflow = .true.
                !% --- weir overflow based on weir length as circumference of storage plan area
                !%     Q = C L H^{3/2}
                !%     L = 2 (pi A)^{1/2}
                !%     1.38 < C < 1.83 from Brater and King Table 5.1 where L in m
                !%     and C in m^{1/2} / S
                !%     Q = coef1 (sqrt(A)) H^(3/2)
                !%     minus sign as Q is outflow (negative)
                lljunction_main_Qoverflow                                &
                    = -coef1 * WeirFactor * sqrt(Storage_Plan_Area)      &
                       * (OverflowDepth**threehalfR)
                return

            case (PondedWeir)     
                isPonding  = .true.
                isOverflow = .false.      
                !% --- similar to OverflowWeir above, but allows in or outflow 
                if (ExternalPondedHeadDiff < zeroR) then 
                    !% --- inflow from ponding (+ value)
                    tempOverflow                                             &
                        = +coef1 * WeirFactor * sqrt(Storage_Plan_Area)      &
                            * (-ExternalPondedHeadDiff)**threehalfR

                    !% --- limit inflow by ponded water available
                    !if ((tempOverflow * dt *crk)  > VolumePondedTotal) then
                    if ((tempOverflow * dt)  > VolumePondedTotal) then
                        !lljunction_main_Qoverflow = VolumePondedTotal / (dt * crk)
                        lljunction_main_Qoverflow = VolumePondedTotal / dt
                    else 
                        lljunction_main_Qoverflow = tempOverflow
                    end if

                    return

                elseif (ExternalPondedHeadDiff > zeroR) then 
                    !% --- outflow to ponding (- value)
                    lljunction_main_Qoverflow                                &
                        = -coef1 * WeirFactor * sqrt(Storage_Plan_Area)      &
                            * (ExternalPondedHeadDiff**threehalfR)
                    return

                else 
                    !% --- no exchange
                    lljunction_main_Qoverflow = zeroR
                    return
                end if

            case (OverflowOrifice)
                isPonding  = .false.
                isOverflow = .true.
                !% --- orifice overflow assuming a single orifice
                !%     Using Brater and King Eq. 4.17 
                !%     Q = (2/3) L sqrt(2g) H^(3/2)
                !%     define coef3 = (2/3) sqrt(2g) so that
                !%     Q = coef3 * L * H^{3/2}
                !%     minus sign as Q is outflow (negative)
                if (elemR(JMidx,er_Head) .le. (MinHeadForOverflowPonding + Horifice )) then
                    !% --- head (water surface) below upper edge of orifice (H < Z + Horifice)
                    lljunction_main_Qoverflow = -coef3 * Lorifice * (OverflowDepth**threehalfR)
                    return 

                else
                    !% --- orifice is pressurized (H > Zcrown + Horifice)
                    !%     Using Brater and King eq. 4.16
                    !%     Q = (2/3) L sqrt(2g) ( (H-Zcrown)^(3/2) - (H-(zcrown+Horifice))^(3/2) )
                    !%     define coef3 = (2/3) sqrt(2g) so that
                    !%     Q = coef3 * L * ( (H-Zcrown)^(3/2) - (H-(zcrown+Horifice))^(3/2) )
                    !%     minus sign as Q is outflow (negative) 
                    lljunction_main_Qoverflow = -coef3 * Lorifice  &
                        * (                                                                                 &
                            +  ((OverflowDepth           )**threehalfR)   &
                            -  ((OverflowDepth - Horifice)**threehalfR)    &
                          )
                    return
                end if
             
            case (PondedOrifice)
                isPonding  = .true.
                isOverflow = .false.
                !% --- similar to OverflowOrifice but allows in or outflows with ponding
                if (ExternalPondedHeadDiff < zeroR) then
                    !% --- inflows from ponding require + value
                    if (ExternalPondedDepth .le. Horifice ) then
                        !% --- water surface below upper edge of orifice (H < Z + Horifice)
                        tempOverflow = +coef3 * Lorifice * ((-ExternalPondedHeadDiff)**threehalfR)

                    else
                        !% --- inflows from ponding with submerged orifice
                        tempOverflow = +coef3 * Lorifice                        &
                            * (                                                 &
                                +  ((-ExternalPondedHeadDiff           )**threehalfR)   &
                                -  ((-ExternalPondedHeadDiff - Horifice)**threehalfR)   &
                            )
                    end if

                    !% --- limit inflow by ponded water available
                    !if ((tempOverflow * dt *crk)  > VolumePondedTotal) then
                    if ((tempOverflow * dt *crk)  > VolumePondedTotal) then
                       !lljunction_main_Qoverflow = VolumePondedTotal / (dt * crk)
                        lljunction_main_Qoverflow = VolumePondedTotal / dt
                    else 
                        lljunction_main_Qoverflow = tempOverflow
                    end if
                    return 

                elseif (ExternalPondedHeadDiff > zeroR) then
                    !% --- outflows to ponding with - value
                    if (ExternalPondedDepth .le. Horifice) then
                        !% --- water surface below upper edge of orifice (H < Z + Horifice)
                        lljunction_main_Qoverflow = -coef3 * Lorifice &
                            * ((ExternalPondedHeadDiff)**threehalfR)
                        return 

                    else
                        lljunction_main_Qoverflow = -coef3 * Lorifice  &
                            * (                                                                                 &
                                +  ((ExternalPondedHeadDiff           )**threehalfR)   &
                                -  ((ExternalPondedHeadDiff - Horifice)**threehalfR)    &
                            )
                        return
                    end if

                else 
                    !% --- no exchange
                    lljunction_main_Qoverflow = zeroR
                    return
                end if

            case default
                !% --- should not reach here.
                !%     for debugging, change to impure function and uncomment
                !%     the following
                print *, 'CODE ERROR unexpected case default'
                call util_crashpoint(6209874)
        end select

    end function lljunction_main_Qoverflow   
!%    
!%==========================================================================
!%==========================================================================
!% 
    subroutine lljunction_main_slotwidth (thisColP) 
        !%-----------------------------------------------------------------
        !% Description
        !% ensures the Preissman Slot width is set correctly for non-surcharged
        !% junctions. This is needed so that incipient surcharge starts with
        !% the correct values
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisColP  !% should be only JM elements
            integer, pointer :: Npack, thisP(:)
        !%-----------------------------------------------------------------
        !% Preliminaries
            Npack => npack_elemP(thisColP)
            if (Npack < 1) return
            thisP => elemP(1:Npack,thisColP)
        !%-----------------------------------------------------------------

        where (.not. elemYN(thisP,eYN_isSurcharged))
            elemR(thisP,er_SlotWidth) = sqrt(elemSR(thisP,esr_Storage_Plan_Area))
        endwhere

    end subroutine lljunction_main_slotwidth
!%    
!%==========================================================================
!%==========================================================================
!% 
    real(8) function lljunction_main_sumBranches (JMidx,thisCol,thisArray)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the net values for a junction from all the branches
        !% Note -- assumes that all branches that do not exist have zero 
        !% values
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx, thisCol
            real(8), intent(in) :: thisArray(:,:)
            real(8) :: thissum
            integer :: kk
        !%------------------------------------------------------------------
        thissum = zeroR

        do kk = 1,max_branch_per_node 
            if (elemSI(JMidx+kk,esi_JB_Exists) .ne. oneI) cycle
            thissum = thissum + branchsign(kk) * thisArray(JMidx+kk,thisCol)
        end do

        lljunction_main_sumBranches = thissum

    end function lljunction_main_sumBranches   
!%    
!%========================================================================== 
!%==========================================================================
!% 
    real(8) function lljunction_main_sumBranches_InOrOutFlow &
         (JMidx,thisCol,thisArray,isInflow,isApplyBranchSign)
        !%------------------------------------------------------------------
        !% Description:
        !% sums real data in thisCol of thisArray for branches that are inflows
        !% if isInflow is true or outflows if isInflow is false.
        !% if isApplyBranchSign = true, the "thisCol" is a flowrate and the
        !% branch sign is applied in the summation to account for the different
        !% sign of inflows/outflows on upstream and downstream branches
        !% Returns > 0 for inflows and < 0 for outflows
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: JMidx, thisCol
            real(8), intent(in) :: thisArray(:,:)
            logical, intent(in) :: isInflow, isApplyBranchSign
            
            real(8) :: dsign, tbranchsign(max_branch_per_node), thissum 
            integer :: kk
        !%------------------------------------------------------------------
        !% Preliminaries
            if (isApplyBranchSign) then
                !% --- apply + to upstream branches and - to downstream
                !%     which is needed to make flowrates consistent 
                !%     inflows or outflows.
                tbranchsign = branchsign
            else 
                tbranchsign(:) = +oneR
            end if
        !%------------------------------------------------------------------
        thissum = zeroR

        !% --- set the sign for flowrate  
        if (isInflow) then 
            dsign = +oneR
        else 
            dsign = -oneR
        end if

        !% --- cycle through upstream branches
        do kk = 1,max_branch_per_node 
            !% --- check if branch exists
            if (elemSI(JMidx+kk,esi_JB_Exists) .ne. oneI) cycle
            !% --- check the screening based on flowrate
            !%     Note that branchsign*Q is positive for inflows and negative
            !%     for outflows, dsign = -1 switches outflows to positive and inflows to negative
            if ((dsign * branchsign(kk) * elemR(JMidx+kk,er_Flowrate)) > zeroR) then
                !% --- this is an inflow and we're looking for inflows
                !%     or this is an outflow and we are looking for outflows
                !%     Note tbranchsum is used so that the "thisCol" data is only reversed
                !%     in sign if it is a flowrate (e.g., does not apply to summing areas)
                thissum = thissum + tbranchsign(kk) * thisArray(JMidx+kk,thisCol)
            else
                !% --- is an outflow and we're looking for inflows or vice versa 
            end if
        end do

        lljunction_main_sumBranches_InOrOutFlow = thissum
        


    end function lljunction_main_sumBranches_InOrOutFlow
!%    
!%========================================================================== 
!%==========================================================================
!% 
    subroutine lljunction_main_update_Qdependent_values (JMidx,istep) 
        !%------------------------------------------------------------------
        !% Description:
        !% Updates Qrate values after Q branches etc. have changed
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx, istep
            real(8), pointer    :: Qstorage, Qoverflow, dt, crk
        !%------------------------------------------------------------------
        !% Alias
            Qstorage  => elemSR(JMidx,esr_JM_StorageRate)
            Qoverflow => elemSR(JMidx,esr_JM_OverflowPondingRate)
            dt        => setting%Time%Hydraulics%Dt
            crk       => setting%Solver%crk2(istep)
        !%------------------------------------------------------------------

        !% --- update storage volume
        if (elemSI(JMidx,esi_JM_Type) .ne. NoStorage) then
            elemR(JMidx,er_Volume) = elemR(JMidx,er_Volume_N0)                         &
                + dt *crk * Qstorage
        else 
            !% --- no storage does not change volume
            elemR(JMidx,er_Volume) = elemR(JMidx,er_Volume_N0)
        end if

        !% --- update this time step volume overflow and ponded
        select case (elemSI(JMidx,esi_JM_OverflowType))
            case (NoOverflow)
                !% -- no action

            case (OverflowWeir,OverflowOrifice)
                !% --- define the VolumeOverflow as positive
                !elemR(JMidx,er_VolumeOverFlow) = -Qoverflow * dt * crk
                elemR(JMidx,er_VolumeOverFlow) = -Qoverflow * dt

                ! if ((setting%Time%Step > 54165) .and. (JMidx == 109)) then 
                !     print *, 'Volume Overflow AAA', elemR(109,er_VolumeOverFlow)
                ! end if

            case (PondedWeir, PondedOrifice)
                !% --- define the VolumePonded as positive is increasing 
                !%     the ponded volume and negative decreasing it.
                !%     This changes the sign of Qoverflow, which is 
                !%     negative for an outflow rate
                !elemR(JMidx,er_VolumePonded)   = -Qoverflow * dt * crk
                elemR(JMidx,er_VolumePonded)   = -Qoverflow * dt

            case default
                print *, 'CODE ERROR unexpected case default'
                call util_crashpoint(711345)
        end select

    end subroutine lljunction_main_update_Qdependent_values
!%    
!%==========================================================================
!%==========================================================================
!% 
    ! real(8) function lljunction_main_update_storage_rate (JMidx, dH, QnetBranches, istep)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Updates the Qstorage for change in head dH
    !     !%------------------------------------------------------------------
    !         integer, intent(in) :: JMidx, istep
    !         real(8), intent(in) :: dH, QnetBranches
    !         real(8)             :: planArea
    !     !%------------------------------------------------------------------
    !     select case (elemSI(JMidx,esi_JM_Type))

    !         case (NoStorage)
    !             lljunction_main_update_storage_rate = zeroR
    !             return
    !         case (TabularStorage,FunctionalStorage,ImpliedStorage)
    !             !% --- compute the storage flowrate
    !             if (istep == 1) then
    !                 !% --- set the plan area
    !                 if (elemR(JMidx,er_Head) < elemR(JMidx,er_Zcrown)) then
    !                     !% --- for non-surcharged
    !                     planArea = elemSR(JMidx,esr_Storage_Plan_Area)
    !                 elseif (elemR(JMidx,er_Head) > elemR(JMidx,er_Zcrown)) then
    !                     !% --- for surcharged
    !                     planArea = elemSR(JMidx,esr_JM_Present_PlanArea)
    !                 else
    !                     !% --- head is exactly at crown
    !                     if (dH > zeroR) then 
    !                         !% --- rising head use the surcharge area
    !                         planArea = elemSR(JMidx,esr_JM_Present_PlanArea)
    !                     else
    !                         !% --- dropping head use the standard area
    !                         planArea = elemSR(JMidx,esr_Storage_Plan_Area)
    !                     end if
    !                 end if
    !                 !% --- on a rising surcharge, the surcharge plan area is not
    !                 !%     yet set, so use half of the storage plan area
    !                 if (planArea .eq. zeroR) then 
    !                     planArea = onehalfR * elemSR(JMidx,esr_Storage_Plan_Area)
    !                 end if

    !                 ! lljunction_main_update_storage_rate = planArea * dH &
    !                 !         /(setting%Solver%crk2(istep) * setting%Time%Hydraulics%Dt)

    !                 lljunction_main_update_storage_rate = planArea * dH &
    !                         / setting%Time%Hydraulics%Dt       

    !                         ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
    !                         !     print *, 'plan area, dH ', planArea, dH
    !                         ! end if
                
    !                 return

    !             elseif (istep == 2) then 
    !                 !% --- compute rate from mass conservation
    !                 lljunction_main_update_storage_rate = QnetBranches       &
    !                     + elemSR(JMidx,esr_JM_OverflowPondingRate) &
    !                     + elemR(JMIdx,er_FlowrateLateral)

    !             else 
    !                 !% should not be possible
    !                 print *, 'CODE ERROR unexpected else'
    !                 call util_crashpoint(6209874)
    !                 return
    !             end if

    !         case default
    !             lljunction_main_update_storage_rate = zeroR
    !             print *, 'CODE ERROR unexpected case default'
    !             call util_crashpoint(8852783)
    !             return

    !     end select

    ! end function lljunction_main_update_storage_rate   
!%    
!%========================================================================== 
!%==========================================================================
!%   
    subroutine lljunction_main_velocity (thisColP)
        !%-----------------------------------------------------------------
        !% Description
        !% computes an approximate velocity in a junction main for use
        !% in downstream fluxes.
        !% This should depend only on JB and not on face values
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisColP  !% should be only JM elements

            integer, pointer :: Npack
            integer :: mm, ii, kk, rr, JMidx, nInflow

            real(8), dimension(max_branch_per_node) :: inVelocity, inFlow
            real(8) :: weightedVelocity, massVelocity, totalQin
            real(8) :: junctionVelocity
        !%-----------------------------------------------------------------
        !% Preliminaries
            Npack => npack_elemP(thisColP)
            if (Npack < 1) return
        !%-----------------------------------------------------------------

        do mm=1,Npack 
            JMidx = elemP(mm,thisColP)

            if (elemYN(JMidx,eYN_isZeroDepth)) cycle

            elemR(JMidx,er_Flowrate) = zeroR 
            elemR(JMidx,er_Velocity) = zeroR

            inFlow     = zeroR
            inVelocity = zeroR

            !% --- cycle thru branches to get inflow velocities and flowrates
            !%     note that the net inflow velocity is directional (i.e. an inflow
            !%     from downstream against the flow cancels velocity of inflows
            !%     from upstream with the flow)
            kk=1
            do ii=1,max_branch_per_node
                if (elemSI(JMidx+ii,esi_JB_Exists) .ne. oneI)  cycle
                if (elemR(JMidx+ii,er_Depth) .le. setting%ZeroValue%Depth) cycle
                if (mod(ii,2) == 0) then 
                    !% --- downstream branch
                    if (elemR(JMidx+ii,er_Flowrate) < zeroR) then 
                        !% --- downstream branch inflow (do not reverse inflow velocity sign, see note above)
                        inVelocity(kk) =  elemR(JMidx+ii,er_Velocity)
                        inFlow(kk)     = -elemR(JMidx+ii,er_Flowrate)
                        kk=kk+1
                    else 
                        !% downstream outflow
                    end if
                else 
                    !% --- upstream branch
                    if (elemR(JMidx+ii,er_Flowrate) > zeroR) then 
                        !% --- upstream branch inflow
                        inVelocity(kk) = elemR(JMidx+ii,er_Velocity)
                        inFlow(kk)     = abs(elemR(JMidx+ii,er_Flowrate))
                        kk=kk+1
                    else 
                        !% upstream outflow
                    end if
                end if
            end do
            nInflow = kk-1

            !% --- sum the inflows
            if (nInflow > 0) then 
                totalQin = sum(inFlow(1:nInflow))
            else
                totalQin = zeroR
            end if

            ! if (JMidx == printJM) print *, 'total Q ',totalQin

            !% --- add lateral inflow
            if (elemR(JMidx,er_FlowrateLateral) > zeroR) then
                totalQin = totalQin + elemR(JMidx,er_FlowrateLateral)
            else 
                !% --- nothing to add 
            end if

            ! if (JMidx == printJM) print *, 'total Q ',totalQin

            if (totalQin > zeroR) then 
                !% --- weightedvelocity is the flow-weighted average velocity of the inflows
                !% --- HACK we need some branch reduction factors for more than 1 inflow branch
                weightedVelocity = sum(inFlow(1:nInflow) * inVelocity(1:nInflow) ) / totalQin

                ! if (JMidx == printJM) then 
                !     do rr = 1,nInflow
                !         print *, 'inFlow ',inFlow(rr), inVelocity(rr)
                !     end do
                ! end if
                ! if (JMidx == printJM) print *, 'weighted V ',weightedVelocity

                !% --- massVelocity is the velocity for the flow given the junction geometry
                !%     note that this is guaranteed positive
                massVelocity = totalQin / ( sqrt(elemSR(JMidx,esr_Storage_Plan_Area)) * elemR(JMidx,er_Depth))

                ! if (JMidx == printJM) print *, 'massVelocity ',massVelocity

                if (weightedVelocity == zeroR) then
                    junctionVelocity = massVelocity
                elseif (weightedVelocity < zeroR) then 
                    !% -- net upstream velocity
                    junctionVelocity = onehalfR * (weightedVelocity - massVelocity)
                else 
                    junctionVelocity = onehalfR * (weightedVelocity + massVelocity)
                end if

                !% --- excessive velocities are likely in error, so set to zero
                if (abs(junctionVelocity) > setting%Limiter%Velocity%Maximum) then
                    junctionVelocity = zeroR
                end if

                elemR(JMidx,er_Velocity) = junctionVelocity

                !% --- use the smaller magnitude of the total inflow rate or the flowrate implied by the junction velocity
                if (junctionVelocity < zeroR) then
                    elemR(JMidx,er_Flowrate) = max(-totalQin, &
                                                    junctionVelocity * elemR(JMidx,er_Depth) * sqrt(elemSR(JMidx,esr_Storage_Plan_Area)))
                elseif (junctionVelocity > zeroR) then
                    elemR(JMidx,er_Flowrate) = min(totalQin, &
                                                    junctionVelocity * elemR(JMidx,er_Depth) * sqrt(elemSR(JMidx,esr_Storage_Plan_Area)))
                else 
                    elemR(JMidx,er_Flowrate) = zeroR
                end if
                elemR(JMidx,er_FroudeNumber) = junctionVelocity / sqrt(setting%Constant%gravity * elemR(JMidx,er_Depth))


            else
                !% --- no inflows to provide momentum, 
                !%     no changes to flowrate, velocity or froude number
            end if

            ! if (JMidx == printJM) print *,  'flowrate ',elemR(JMidx,er_Flowrate)
        end do


    end subroutine lljunction_main_velocity
!%
!%==========================================================================
!%==========================================================================
!% 
    real(8) function lljunction_main_volume_from_storageRate (JMidx,istep)
        !%------------------------------------------------------------------
        !% Description: updates volume using the storage rate
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx, istep
            real(8), pointer :: dt, crk(:), Qstorage
        !%------------------------------------------------------------------
        !% Aliases
            Qstorage  => elemSR(JMidx,esr_JM_StorageRate)
            dt        => setting%Time%Hydraulics%Dt
            crk       => setting%Solver%crk2
        !%------------------------------------------------------------------

        if (elemSI(JMidx,esi_JM_Type) .ne. NoStorage) then
            lljunction_main_volume_from_storageRate = elemR(JMidx,er_Volume_N0)                         &
                + dt * crk(istep) * Qstorage
        else 
            !% --- no storage does not change volume
            elemR(JMidx,er_Volume) = elemR(JMidx,er_Volume_N0)
        end if

    end function lljunction_main_volume_from_storageRate    
!%    
!%========================================================================== 
!%==========================================================================
!%
    subroutine lljunction_push_inflows_from_CC_to_JB_face ()
        !%-----------------------------------------------------------------
        !% Description:
        !% Pushes the inflow flowrate from a CC element adjacent to a JB
        !% element to the CC/JB face. This is done without interpolation to
        !% ensure that inflows are driven by the upstream conditions.
        !% Note that Diagnostic elements adjacent to JB are handled 
        !% elsewhere
        !%-----------------------------------------------------------------
        !% Declarations
            integer, pointer :: Npack, fidx(:), thisP(:)
            integer :: ii, thisPCol, eiFace
            real(8) :: bsign
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------

        !% --- cycle through upstream and downstream
        do ii=1,2 
            if (ii==1) then 
                !% --- upstream CC elements
                thisPCol = ep_CC_UpstreamOfJunction
                !% --- CC/JB face relative to CC element
                eiFace   = ei_Mface_dL
                !% --- branch sign (sign of an inflow)
                bsign    = +oneR
            else 
                !% --- downstream CC elements
                thisPCol = ep_CC_DownstreamOfJunction
                !% --- CC/JB face relative to CC element
                eiFace   = ei_Mface_uL
                !% --- branch sign (sign of an inflow)
                bsign    = -oneR
            end if

            !% --- size of set of CC elements adjacent to JB
            Npack => npack_elemP(thisPCol)

            if (Npack > 0) then 
                !% --- set of CC elements adjacent to JB
                thisP => elemP(1:Npack,thisPCol)
                !% --- face indx for CC/JB face
                fidx  => elemI(:,eiFace)

                !% --- push element flowrate to face for inflows
                where (bsign * elemR(thisP,er_Flowrate) > zeroR)
                    faceR(fidx(thisP),fr_Flowrate) = elemR(thisP,er_Flowrate)
                    !% --- Note where a shared face has been diverged 
                    !%    from value on other image
                    where (faceYN(fidx(thisP),fYN_isSharedFace))
                        faceYN(fidx(thisP),fYN_isSharedFaceDiverged) = .true.
                    end where

                endwhere
            end if
        end do

    end subroutine lljunction_push_inflows_from_CC_to_JB_face
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine lljunction_push_adjacent_CC_elemdata_to_face ()
        !%-----------------------------------------------------------------
        !% Description:
        !% Pushes the elem data from the element upstream or downstream
        !% of a JB to the JB face for storage in the fr_..._Adjacent
        !% data columns
        !%-----------------------------------------------------------------
        !% Declarations
            integer :: ii, epCCcol
            logical :: isUpstreamFace
        !%-----------------------------------------------------------------

        do ii=1,2
            !% --- cycle over upstream and downstream faces of junction
            if (ii==1) then 
                !% --- an upstream face is CC downstream of junction
                isUpstreamFace = .true.
                epCCcol = ep_CC_DownstreamOfJunction
            else
                !% -- a downstream face is CC upstream of junction
                isUpstreamFace = .false.
                epCCcol = ep_CC_UpstreamOfJunction
            end if
            call face_push_elemdata_to_face (epCCcol, fr_Head_Adjacent,      er_Head,         elemR, isUpstreamface)
            call face_push_elemdata_to_face (epCCcol, fr_EnergyHead_Adjacent,er_EnergyHead,   elemR, isUpstreamface)
            call face_push_elemdata_to_face (epCCcol, fr_Topwidth_Adjacent,  er_Topwidth,     elemR, isUpstreamface)
            call face_push_elemdata_to_face (epCCcol, fr_Length_Adjacent,    er_Length,       elemR, isUpstreamface)
            call face_push_elemdata_to_face (epCCcol, fr_Zcrest_Adjacent,    er_Zbottom,      elemR, isUpstreamface)
            call face_push_elemdata_to_face (epCCcol, fr_Velocity_Adjacent,  er_Velocity,     elemR, isUpstreamface)
            call face_push_elemdata_to_face (epCCcol, fr_Froude_Adjacent,    er_FroudeNumber, elemR, isUpstreamface)
            call face_push_elemdata_to_face (epCCcol, fr_Depth_Adjacent,     er_Depth,        elemR, isUpstreamface)
        end do

    end subroutine lljunction_push_adjacent_CC_elemdata_to_face
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module junction_lowlevel