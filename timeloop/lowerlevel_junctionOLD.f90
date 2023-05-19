module lowerlevel_junction

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

!%----------------------------------------------------------------------------- 
!% Description:
!% Procedures called from withing junction_elements module
!%----------------------------------------------------------------------------- 

    implicit none

    private

    public :: lljunction_branch_energy_outflow 
    public :: lljunction_branch_dQdH
    public :: lljunction_branch_getface
    public :: lljunction_branch_Qnet
    public :: lljunction_branch_update_DeltaQ
    public :: lljunction_branch_update_flowrate
    
    public :: lljunction_CC_for_JBadjacent

    public :: lljunction_conservation_residual
    public :: lljunction_conservation_fix

    public :: lljunction_main_dQdHoverflow
    public :: lljunction_main_dQdHstorage
    public :: lljunction_main_head_bounds
    public :: lljunction_main_Qoverflow
    public :: lljunction_main_sumBranches
    public :: lljunction_main_update_Qdependent_values
    public :: lljunction_main_update_storage_rate 
    public :: lljunction_main_velocity
    public :: lljunction_main_volume_from_storageRate

    public :: lljunction_push_inflowCC_flowrates_to_face
    public :: lljunction_push_adjacent_elemdata_to_face

    integer :: printJM = 420
    
    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine lljunction_branch_energy_outflow ()
        !%-----------------------------------------------------------------
        !% Description:
        !% Computes an energy-equation outflow for each outflow junction 
        !% branch
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, pointer :: Npack, JMar(:), thisJB(:), fidx(:)
            real(8), pointer :: deltaH(:), grav

            logical :: isUpstreamBranch 

            real(8) :: bsign

            integer :: frHead, frArea, frHeadAdj, ii, mm
        !%-----------------------------------------------------------------
        !% Aliases
            !% --- array for the JM index
            JMar   => elemSI(:,esi_JunctionBranch_Main_Index)    
            grav   => setting%Constant%gravity
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------

        !% --- cycle through nominal upstream and downstream JB
        !%     This should affect only outflow branches with consistent
        !%     pressure difference with adjacent CC
        do ii=1,2
            !% --- get the upstream or downstream JB elements
            if (ii==1) then 
                !% --- upstream JB
                isUpstreamBranch = .true.
                bsign = oneR
                Npack => npack_elemP(ep_JB_Upstream_CC_Adjacent)
                if (Npack > 0) then 
                    thisJB => elemP(1:Npack,ep_JB_Upstream_CC_Adjacent)
                else
                    cycle
                end if
                !% --- the JB-adjacent face is upstream
                fidx => elemI(:,ei_Mface_uL)
                frHead = fr_Head_u !% -- use u for jump purposes? QUESTION
                frArea = fr_Area_u !% QUESTION
                frHeadAdj = fr_Head_Adjacent
            else
                !% --- downstream JB
                bsign = -oneR
                isUpstreamBranch = .false.
                Npack => npack_elemP(ep_JB_Downstream_CC_Adjacent)
                if (Npack > 0) then 
                    thisJB => elemP(1:Npack,ep_JB_Downstream_CC_Adjacent)
                else 
                    cycle
                end if
                !% --- the JB-adjacent face is downstream
                fidx => elemI(:,ei_Mface_dL)
                frHead = fr_Head_d !% --- use d for jump purposes? QUESTION
                frArea = fr_Area_d !% QUESTION
                frHeadAdj = fr_Head_Adjacent
            end if

            !% --- head difference from junction main to element
            deltaH => elemR(:,er_Temp01)
            deltaH = zeroR
            where (elemR(thisJB,er_Depth) > setting%ZeroValue%Depth)
                where (faceR(fidx(thisJB),frHeadAdj) > faceR(fidx(thisJB),fr_Zbottom))
                    deltaH(thisJB) =  elemR(JMar(thisJB),er_Head) - faceR(fidx(thisJB),frHeadAdj)
                elsewhere
                    !% --- where adjacent head is lower than face zbottom
                    deltaH(thisJB) = elemR(JMar(thisJB),er_Head) - (faceR(fidx(thisJB),fr_Zbottom)+ setting%ZeroValue%Depth)
                endwhere
                !deltaH(thisJB) =  elemR(JMar(thisJB),er_Head) - faceR(fidx(thisJB),frHead)
                
            endwhere

            ! do mm=1,size(thisJB)
            !     if (printJM == Jmar(thisJB(mm))) then
            !         if (.not. isUpstreamBranch) then 
            !             print *, 'AAA ', isUpstreamBranch
            !             print *, 'deltaH()  ', deltaH(thisJB(mm)), elemR(thisJB(mm),er_Depth)
            !             print *, 'JM, JB',  JMar(thisJB(mm)), thisJB(mm),  elemR(JMar(thisJB(mm)),er_Head)
            !             print *, ' '
            !            ! print *, 'fidx      ',  fidx(6), faceR(fidx(6),frHead)
            !         end if
            !     end if
            ! end do

            !% --- limit dH if it is more than 1/2 of depth in branch
            where (deltaH(thisJB) > onehalfR   * elemR(thisJB,er_Depth))
                deltaH(thisJB)    = onefourthR * elemR(thisJB,er_Depth)
            end where

            ! do mm=1,size(thisJB)
            !     if (printJM == Jmar(thisJB(mm))) then
            !         if (.not. isUpstreamBranch) then 
            !             print *, 'BBB ', isUpstreamBranch
            !             print *, 'deltaH()  ', deltaH(thisJB(mm)), elemR(thisJB(mm),er_Depth)
            !             !print *, 'JM, JB',  JMar(thisJB(mm)), thisJB(mm),  elemR(JMar(thisJB(mm)),er_Head)
            !             print *, ' '
            !            ! print *, 'fidx      ',  fidx(6), faceR(fidx(6),frHead)
            !         end if
            !     end if
            ! end do

            !% --- Subcritical outflow junction (super overwrites below)
            !%     Requires more than a trivial deltaH and there must be more than
            !%     a trivial outflowrate
            where ((deltaH(thisJB) > 1.0d-8)                        &
                   .and.                                            &
                   (elemR(thisJB,er_Flowrate) * bsign < -1.0d-8)    &
                   )
                    !% --- outflow from JB  in upstream direction
                    !% --- or outflow from JB in downstream direction
                    !%     NOTE: without junction main approach velocity
                    elemR(thisJB,er_Flowrate) = - bsign * faceR(fidx(thisJB),frArea)                  &
                            * sqrt(                                                                     &
                                    (twoR * grav * deltaH(thisJB))                                      &
                                    /(oneR + elemSR(thisJB,esr_JunctionBranch_Kfactor))                 &
                                    )   

                    where (elemR(thisJB,er_Area) > setting%ZeroValue%Area)
                        elemR(thisJB,er_Velocity) = elemR(thisJB,er_Flowrate) / elemR(thisJB,er_Area)
                    elsewhere 
                        elemR(thisJB,er_Velocity) = zeroR
                    endwhere
    
                    !% --- apply strict velocity limiter
                    where (abs(elemR(thisJB,er_Velocity)) > setting%Limiter%Velocity%Maximum)
                        elemR(thisJB,er_Velocity) = sign(setting%Limiter%Velocity%Maximum * 0.99d0, elemR(thisJB,er_Velocity))
                    endwhere

                    !% --- push flowrate and velocity to faces
                    faceR(fidx(thisJB),fr_Flowrate)   = elemR(thisJB,er_Flowrate)
                    faceR(fidx(thisJB),fr_Velocity_u) = elemR(thisJB,er_Velocity)
                    faceR(fidx(thisJB),fr_Velocity_d) = elemR(thisJB,er_Velocity)
            endwhere 

            ! do mm=1,size(thisJB)
            !     if (printJM == Jmar(thisJB(mm))) then
            !         if (.not. isUpstreamBranch) then 
            !             print *, 'CCC ', isUpstreamBranch
            !             print *, 'deltaH()  ', deltaH(thisJB(mm)), elemR(thisJB(mm),er_Depth)
            !             print *, 'Flowrate  ', elemR(thisJB(mm),er_Flowrate)
            !             print *, 'area      ', faceR(fidx(thisJB(mm)),frArea)
            !             !print *, 'JM, JB',  JMar(thisJB(mm)), thisJB(mm),  elemR(JMar(thisJB(mm)),er_Head)
            !             print *, ' '
            !            ! print *, 'fidx      ',  fidx(6), faceR(fidx(6),frHead)
            !         end if
            !     end if
            ! end do



            !% --- for a trivial outflow and a trivial dH, set flow to zero
            !%     Note -- does not affect inflows or deltaH < zeroR
            where ((deltaH(thisJB) .ge. zeroR) .and. (deltaH(thisJB) .le. 1.0d-8) .and. (elemR(thisJB,er_Flowrate) * bsign < -1.0d-8))
                elemR(thisJB,er_Flowrate) = zeroR
                elemR(thisJB,er_Velocity) = zeroR
            end where


            !% HACK --- NEED TO CHECK WHY THIS WORKS?  20230506

            !% --- for inflow pressure gradient and outflow
            where  ((deltaH(thisJB) < -1.0d-8) .and. (elemR(thisJB,er_Flowrate) * bsign .le. -1.0d-8))
                !% --- inconsistent flow direction and pressure gradient   
                elemR(thisJB,er_Flowrate) = zeroR
            endwhere


            !%STUFF THAT WAS TRIED BUT DIDN'T HELP:
            ! where ((deltaH(thisJB) > 1.0d-8) .and. (elemR(thisJB,er_Flowrate) * bsign > -1.0d-8))
            !     !% --- inconsistent flow direction and pressure gradient   
            !     elemR(thisJB,er_Flowrate) = zeroR
            ! endwhere 

            ! where  ((deltaH(thisJB) .le. 1.0d-8) .and. (deltaH(thisJB) .ge. -1.0d-8))
            !     !% --- effectively zero pressure gradient (applies to JB as zerodepth)
            !     elemR(thisJB,er_Flowrate) = zeroR
            ! endwhere 

            ! where ((deltaH(thisJB) < -1.0d-8) .and. (elemR(thisJB,er_Flowrate) * bsign > -1.0d-8))
            !     !% --- inflow using face flowrate
            !     elemR(thisJB,er_Flowrate) =  faceR(fidx(thisJB),fr_Flowrate) 
            ! endwhere

            !% THE SUPERCRITICAL CAUSED PROBLEMS FOR THIN LAYERS
            ! !% --- Supercritical junction (requires JM velocity and froude number)
            ! if (isUpstreamBranch) then
            !     where (elemR(JMar(thisJB),er_FroudeNumber) < -oneR)
            !         elemR(thisJB,er_Flowrate) = elemR(JMar(thisJB),er_Velocity) * faceR(fidx(thisJB),frArea)
            !         !elemR(thisJB,er_Flowrate) = elemR(thisJB,er_Velocity) * faceR(fidx(thisJB),frArea)
            !     endwhere     
            ! else
            !     where (elemR(JMar(thisJB),er_FroudeNumber) > oneR)
            !         elemR(thisJB,er_Flowrate) = elemR(JMar(thisJB),er_Velocity) * faceR(fidx(thisJB),frArea)
            !         !elemR(thisJB,er_Flowrate) = elemR(thisJB,er_Velocity) * faceR(fidx(thisJB),frArea)
            !     endwhere     
            ! end if  

            ! print *, 'after ',elemR(7,er_Flowrate)
            ! print *, ' '

            !% THIS CAUSED HUGE PROBLEMS
            ! ! !% --- volume limiter on outflow
            ! if (isUpstreamBranch) then 
            ! !     ! !% --- volume limit on outflow from upstream branch
            ! !     ! where ( elemR(JMar(thisJB),er_Volume) < -dt * elemR(thisJB,er_Flowrate) )
            ! !     !     elemR(thisJB,er_Flowrate) = - elemR(JMar(thisJB),er_Volume) / (dt * twoR)
            ! !     ! endwhere
            ! else 
            ! ! 20230430brh
            !     !% --- volume limit on outflow on downstream branch
            !     !where ( elemR(JMar(thisJB),er_Volume) < dt * elemR(thisJB,er_Flowrate) )
            !     !    elemR(thisJB,er_Flowrate) = elemR(JMar(thisJB),er_Volume) / (dt * twoR)
            !     do mm=1,size(thisJB)
            !         print *, thisJB(mm), elemR(thisJB(mm),er_Flowrate), elemR(JMar(thisJB(mm)),er_Volume) / (dt * twoR)
            !     end do
            !     !endwhere
            ! end if

            ! !% --- update velocity
            ! where (elemR(thisJB,er_Area) > setting%ZeroValue%Area)
            !     elemR(thisJB,er_Velocity) = elemR(thisJB,er_Flowrate) / elemR(thisJB,er_Area)
            ! elsewhere 
            !     elemR(thisJB,er_Velocity) = zeroR
            ! endwhere

            ! !% --- apply velocity limiter
            ! where (abs(elemR(thisJB,er_Velocity)) > setting%Limiter%Velocity%Maximum)
            !     elemR(thisJB,er_Velocity) = sign(setting%Limiter%Velocity%Maximum * 0.99d0, elemR(thisJB,er_Velocity))
            ! endwhere

        end do

        Npack => npack_elemP(ep_JB)
        if (Npack > 0) then 
            thisJB => elemP(1:Npack, ep_JB)
            call update_Froude_number_element (thisJB) 
            call update_wavespeed_element (thisJB)
            !call update_interpweights_JB (thisP, Npack, .true.) !% not needed here 20230504brh
        end if
        
    end subroutine lljunction_branch_energy_outflow
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
            real(8)          :: bsign, denominator, thisArea, FrFactor
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
            JMidx => elemSI(JBidx,esi_JunctionBranch_Main_Index)

            headJM => elemR(JMidx,er_Head)

            !% --- get the face for this JB
            if (elemSI(JBidx,esi_JunctionBranch_IsUpstream) == oneI) then 
                !% --- upstream JB
                isDownstream = .false.
                fidx  =  elemI(JBidx,ei_Mface_uL)
                fA    => faceR(fidx,fr_Area_u)
                fH    => faceR(fidx,fr_Head_u)
                fQ    => faceR(fidx,fr_Flowrate)
                fZ    => faceR(fidx,fr_Zbottom)
                FrAdj => faceR(fidx,fr_Froude_Adjacent)
                !bsign = +oneR
                !bsign = -oneR
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
                !bsign = -oneR
                !bsign = +oneR
                if (fQ < 0) then 
                    isInflow = .true.
                else 
                    isInflow = .false.
                end if

                ! print *, ' '
                ! print *, 'in dqdh ',isInflow, JBidx, fidx
                ! print *, ' '

            end if
            !% --- adjacent element data
            Ladj  => faceR(fidx,fr_Length_Adjacent)
            Tadj  => faceR(fidx,fr_Topwidth_Adjacent)
            Hadj  => faceR(fidx,fr_Head_Adjacent)
            VelAdj=> faceR(fidx,fr_Velocity_Adjacent)
            Dadj  => faceR(fidx,fr_Depth_Adjacent)
            ZBadj => faceR(fidx,fr_Zcrest_Adjacent)


            !% --- CC elements adjacent to JB
            if (elemSI(JBidx,esi_JunctionBranch_CC_adjacent) == oneI) then 

                if (isDownstream) then 

                    ! print *, 'FrAdj ',FrAdj

                    !% --- downstream branch
                    if ((FrAdj .le. -oneR) .or. (elemR(JBidx,er_FroudeNumber) .le. -oneR)) then 
                        !% --- supercritical inflow
                        elemSR(JBidx,esr_JunctionBranch_dQdH) = zeroR
                    else 
                        !% --- outflow or subcritical inflow
                        elemSR(JBidx,esr_JunctionBranch_dQdH) = + crk * grav * dt * fA / Ladj
                    end if

                    ! if (FrAdj .ge. oneR) then  
                    !     !% --- supercritical outflow  (HACK NEEDS BETTER TESTING BEFORE USING)
                    !     !elemSR(JBidx,esr_JunctionBranch_dQdH) =  VelAdj * Tadj * FrAdj  TEST MOD 20230430brh
                    !     elemSR(JBidx,esr_JunctionBranch_dQdH) = crk * grav * dt * fA / Ladj
                    ! elseif (FrAdj .le. -oneR) then 
                    !     !% --- supercritical inflow
                    !     elemSR(JBidx,esr_JunctionBranch_dQdH) = zeroR
                    ! else 
                    !     !% --- subcritical in/out flow
                    !     elemSR(JBidx,esr_JunctionBranch_dQdH) = crk * grav * dt * fA / Ladj
                    ! end if

                    ! print *, 'dqdh AA ',elemSR(JBidx,esr_JunctionBranch_dQdH)

                    !% --- handle waterfall inflow elements
                    if ((isInflow) .and. (headJM < fZ)) then 
                        elemSR(JBidx,esr_JunctionBranch_dQdH) = zeroR
                    end if

                    ! print *, 'dqdh BB ',elemSR(JBidx,esr_JunctionBranch_dQdH)

                    !% --- handle uphill outflow
                    if ((.not. isInflow) .and. (headJM < Hadj)) then 
                        elemSR(JBidx,esr_JunctionBranch_dQdH) = zeroR 
                    end if

                    !% --- TEST limit uphill outflow   20230506brh
                    if ( (.not. isInflow) .and. (headJM < Hadj)) then 
                        !% --- outflow into an adverse pressure gradient,
                        !%     largest negative allowable dQdH is Q/deltaH
                        elemSR(JBidx,esr_JunctionBranch_dQdH) = min(elemSR(JBidx,esr_JunctionBranch_dQdH), elemR(JBidx,er_Flowrate) / (Hadj - headJM))
                    end if

                    ! print *, 'dqdh CC ',elemSR(JBidx,esr_JunctionBranch_dQdH)
                    ! print *, JBidx, isInflow, headJM, Hadj
                    ! print *,' '

                    !% --- limit outflow by 1/4 volume of JM
                    !%     Qmax = (1/4) V / dt;  V = Aplan * Depth
                    !%     Q/H = 1/4 (Aplan * Depth) / (dt * Depth) = (1/4) Aplan / dt
                    if (.not. isInflow) then   !% CHANGE 20230506brh
                        elemSR(JBidx,esr_JunctionBranch_dQdH) &
                            = min (elemSR(JBidx, esr_JunctionBranch_dQdH),  &
                                    elemSR(JMidx, esr_Storage_Plan_Area) * onefourthR / dt)
                    end if

                    ! print *, 'dqdh DD ',elemSR(JBidx,esr_JunctionBranch_dQdH)
                    ! print *, JBidx, isInflow, headJM, Hadj
                    ! print *,' '

                else
                    !% --- downstream branch
                    if ((FrAdj .ge. +oneR) .or. (elemR(JBidx,er_FroudeNumber) .ge. +oneR)) then 
                        !% --- supercritical inflow
                        elemSR(JBidx,esr_JunctionBranch_dQdH) = zeroR
                    else
                        !% --- outflow or subcritical inflow
                        elemSR(JBidx,esr_JunctionBranch_dQdH) = - crk* grav * dt * fA / Ladj
                    end if


                    ! if (FrAdj .le. -oneR) then 
                    !     !% --- supercritical outflow (HACK NEEDS BETTER TESTING BEFORE USING)
                    !     !elemSR(JBidx,esr_JunctionBranch_dQdH) = VelAdj * Tadj * abs(FrAdj)  TEST MOD 20230430brh
                    !     elemSR(JBidx,esr_JunctionBranch_dQdH) = - crk* grav * dt * fA / Ladj
                    ! elseif (FrAdj .ge. +oneR) then 
                    !     !% --- supercritical inflow
                    !     elemSR(JBidx,esr_JunctionBranch_dQdH) = zeroR
                    ! else
                    !     !% --- subcritical in/out flow
                    !     elemSR(JBidx,esr_JunctionBranch_dQdH) = - crk* grav * dt * fA / Ladj
                    ! end if

                    !% --- handle waterfall inflow elements
                    if ((isInflow) .and. (headJM < fZ)) then 
                        elemSR(JBidx,esr_JunctionBranch_dQdH) = zeroR
                    end if

                    !% --- handle uphill outflow
                    ! if ((.not. isInflow) .and. (headJM < Hadj)) then 
                    !     elemSR(JBidx,esr_JunctionBranch_dQdH) = zeroR 
                    ! end if

                    !% --- TEST limit uphill outflow  20230506brh
                    if ( (.not. isInflow) .and. (headJM < Hadj)) then 
                        !% --- outflow into adverse pressure gradient
                        !%     largest (positive) dQ/dH is Q/deltaH
                        elemSR(JBidx,esr_JunctionBranch_dQdH) = max(elemSR(JBidx,esr_JunctionBranch_dQdH), elemR(JBidx,er_Flowrate) / (Hadj-headJM))
                    end if

                    !% --- limit outflow dQdH by 1/4 volume of JM
                    !%     Qmax = (1/4) V / dt;  V = Aplan * Depth
                    !%     Q/H = 1/4 (Aplan * Depth) / (dt * Depth) = (1/4) Aplan / dt
                    !%     note upstream branch dQdH < 0
                    if (.not. isInflow) then   !% CHANGE 20230506brh
                        elemSR(JBidx,esr_JunctionBranch_dQdH) &
                            = max ( elemSR(JBidx, esr_JunctionBranch_dQdH),  &
                                    -elemSR(JMidx, esr_Storage_Plan_Area) * onefourthR / dt)
                    end if

                end if

               
                !thisArea = fA + Tadj*(fH - Hadj)
                ! if (thisArea > zeroR) then 
                    
                    ! if ((elemSI(JBidx,esi_JunctionBranch_IsUpstream) .ne. oneI) .and. (FrAdj .ge. oneR)) then 
                    !     elemSR(JBidx,esr_JunctionBranch_dQdH) = FrAdj * (Tadj + onehalfR / (grav * Dadj))
                    ! else 
                    !     elemSR(JBidx,esr_JunctionBranch_dQdH) = tempfactor * (bsign * crk * dt / Ladj) * (grav * thisArea)
                    ! end if
                    ! if (JMidx == printJM) print *, JBidx, 'dQdH thisarea ',elemSR(JBidx,esr_JunctionBranch_dQdH)
                    ! if (JMidx == printJM) print *, JBidx, 'Froude number ',FrAdj
                    !elemSR(JBidx,esr_JunctionBranch_dQdH) = tempfactor * (bsign * crk * dt / Ladj) * (grav * thisArea + FrFactor * (VelAdj**2) * Tadj)
                    !if (JMidx == printJM) print *, JBidx, 'dqdH A ',(bsign * crk * dt / Ladj) * (grav * thisArea)
                    !if (JMidx == printJM) print *, JBidx, 'dQdH B ',(bsign * crk * dt / Ladj) * (grav * thisArea + FrFac * (VelAdj**2) * Tadj)
                    !if (JMidx == printJM) print *, JBidx, 'Vel1   ',VelAdj
                ! else 
                    ! elemSR(JBidx,esr_JunctionBranch_dQdH) = tempfactor * (bsign * crk * dt / Ladj) * (grav * fA)
                    ! if (JMidx == printJM) print *, JBidx, 'dQdH small A  ',elemSR(JBidx,esr_JunctionBranch_dQdH)
                    ! if (JMidx == printJM) print *, JBidx, 'Froude number ',FrAdj
                    !elemSR(JBidx,esr_JunctionBranch_dQdH) = tempfactor * (bsign * crk * dt / Ladj) * (grav * fA + FrFactor * (VelAdj**2) * Tadj)
                    !if (JMidx == printJM) print *, JBidx, 'dqdH A ',(bsign * crk * dt / Ladj) * (grav * thisArea)
                    !if (JMidx == printJM) print *, JBidx, 'dqdH B ',(bsign * crk * dt / Ladj) * (grav * thisArea + (VelAdj**2) * Tadj)
                    !if (JMidx == printJM) print *, JBidx, 'Vel2   ',VelAdj
                ! end if

            !% --- Diagnostic element adjacent to JB
            elseif (elemSI(JBidx,esi_JunctionBranch_Diag_adjacent)) then 
                !% --- dQdH for adjacent diagnostic element
                if (isInflow) then 
                    bsign = -oneR !% increasing head decreases flowrate for inflow
                else
                    bsign = +oneR !% increasing head increases flowrate for outflow
                endif
                if (elemR(JMidx,er_Head) < faceR(fidx,fr_Zcrest_Adjacent)) then 
                    !% -- insufficient head means the junction cannot affect the flow in 
                    !%    diagnostic branch (either in or out flow)
                    elemSR(JBidx,esr_JunctionBranch_dQdH) = zeroR
                    !cycle
                else
                    elemSR(JBidx,esr_JunctionBranch_dQdH) = bsign * faceR(fidx,fr_dQdH_Adjacent)
                    !cycle
                end if

            else 
                print *, 'CODE ERROR: Unexepected else'
                call util_crashpoint(140987112)
            end if

            ! print *, ' '
            ! print *, 'in junction_branch_dQdH'
            ! print *, JBidx, elemSR(JBidx,esr_JunctionBranch_dQdH)
            ! print *, ' '

        end do
        
    end subroutine lljunction_branch_dQdH
!%
!%==========================================================================
!%==========================================================================
!% 
    pure subroutine lljunction_branch_getface (outdata, frCol, JMidx, fiIdx, kstart)
        !%-----------------------------------------------------------------
        !% Description
        !% Stores face data of frCol on outdata element sapce
        !%-----------------------------------------------------------------
            real(8), intent(inout) :: outdata(:)  !% element data for output
            integer, intent(in)    :: JMidx       !% index of JM junction
            integer, intent(in)    :: fiIdx       !%  index of map up or down to face
            integer, intent(in)    :: frCol    !%  column in faceR array
            integer, intent(in)    :: kstart       !% = 1 for upstream branches, 2 for down 
            real(8) :: k1,k2   
        !%-----------------------------------------------------------------
    
        k1 = JMidx + kstart
        k2 = JMidx + max_branch_per_node

        !do concurrent (kk=kstart:max_branch_per_node:2)
        where (elemSI(k1:k2:2,esi_JunctionBranch_Exists) .eq. oneI)
              outdata(k1:k2:2) = faceR(elemI(k1:k2:2,fiIdx),frCol) 
        endwhere
        !end do

        !print *, 'in get face', outdata(k1:k2:2)

    end subroutine lljunction_branch_getface
!%
!%==========================================================================
!%==========================================================================
!% 
    real(8) pure function lljunction_branch_Qnet (JMidx,idir)    
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the net inflow (idir = 1) or outflow (idir=-1) from
        !% junction JMidx
        !%------------------------------------------------------------------
            integer, intent(in) :: JMidx,idir
            real(8), dimension(max_branch_per_node) :: Qdir
            integer :: kk
        !%------------------------------------------------------------------

        Qdir = zeroR
        do concurrent (kk=1:max_branch_per_node)
            if ((abs(elemR(JMidx+kk,er_Flowrate))               > zeroR) .and. &
                (elemSI(JMidx+kk,esi_JunctionBranch_Exists)     == oneI) .and. &
                (elemSI(JMidx+kk,esi_JunctionBranch_CanModifyQ) == oneI)         ) then

                Qdir(kk) =  branchsign(kk) *  elemR(JMidx+kk,er_Flowrate)      &
                    *onehalfR * (                                              &
                                    oneR + real(idir,8) * branchsign(kk)       &
                                         * elemR(JMidx+kk,er_Flowrate)         &
                                         / abs(elemR(JMidx+kk,er_Flowrate))    &
                                )
                ! print *, 'kk ',kk
                ! print *,  'A ',branchsign(kk) *  elemR(JMidx+kk,er_Flowrate)
                ! print *,  'B ',oneR + real(idir,8) * branchsign(kk) * elemR(JMidx+kk,er_Flowrate)   &
                ! / abs(elemR(JMidx+kk,er_Flowrate)) 
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
        !% Updates the Qfor a JB based on the dQ/dH
        !%------------------------------------------------------------------
            integer, intent(in) :: JMidx
            real(8), intent(in) :: dH
            integer :: ii
        !%------------------------------------------------------------------

        do ii=1,max_branch_per_node
            if (elemSI(JMidx+ii,esi_JunctionBranch_Exists) .ne. oneI) cycle   
            
            elemR(JMidx+ii,er_DeltaQ) = elemSR(JMidx+ii,esr_JunctionBranch_dQdH) * dH 

            ! print *, ' '
            ! print *, 'JB ',JMidx+ii 
            ! print *, 'DH dQdH ', dH, elemSR(JMidx+ii,esr_JunctionBranch_dQdH)
            ! print *, 'Delta Q ',elemR(JMidx+ii,er_DeltaQ)
            ! print *, ' '
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
           integer :: ii
       !%------------------------------------------------------------------

       do ii=1,max_branch_per_node
           if (elemSI(JMidx+ii,esi_JunctionBranch_Exists) .ne. oneI) cycle
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
        !% Adusts the values on the elements that are JB adjacent for the
        !% solution of the JM and JB elements
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisColP, istep
            logical, intent(in) :: isUpstreamYN

            integer, pointer    :: thisCC(:), Npack, fIdx(:)
            real(8), pointer    :: dt, crk, oldVolume(:)
            real(8) :: bsign

            integer :: mm, JMidx
        !%-----------------------------------------------------------------   
        !% Aliases
            Npack       => npack_elemP(thisColP)
            if (Npack < 1) return
            thisCC      => elemP(1:Npack,thisColP)
            dt          => setting%Time%Hydraulics%Dt
            crk         => setting%Solver%crk2(istep)
            oldVolume   => elemR(:,er_Temp01)
        !%-----------------------------------------------------------------  
        !% Preliminaries
        !%----------------------------------------------------------------- 
            
        if (isUpstreamYN) then 
            !% --- pointer to the downstream JB face for CC upstream of JB
                ! print *, ' '
                ! print *, 'UPSTREAM ================================='
            fIdx => elemI(:,ei_Mface_dL)
            bsign = +oneR
        else
            !% --- pointer to the upstrea JB face for CC downstream of JB
                ! print *, ' '
                ! print *, 'DOWNSTREAM ============================'
            fIdx => elemI(:,ei_Mface_uL)
            bsign = -oneR
        end if    

        !% -- adjust the conservative flux on JBadjacent CC for deltaQ
        if (istep == 2) then 
            faceR(fIdx(thisCC),fr_Flowrate_Conservative) &
                 = faceR(fIdx(thisCC),fr_Flowrate_Conservative) &
                 + faceR(fIdx(thisCC),fr_DeltaQ)

                !  if ((istep == 2) .and. (isUpstreamYN)) then     
                !     print *, ' '
                !     print *, 'faceQ after  ',faceR(6,fr_Flowrate_Conservative)
                !  end if
        end if

        !% --- store present volume
        oldVolume(thisCC) = elemR(thisCC,er_Volume)

        !% --- changing the face flowrate changes the volume
        !%     on an upstream element, a negative DeltaQ causes
        !%     an increase in the upstream element volume
        elemR(thisCC,er_Volume) = elemR(thisCC,er_Volume) &
            - bsign * crk * dt *  faceR(fIdx(thisCC),fr_DeltaQ)    

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
            integer, intent(in) :: JMidx
            real(8), intent(inout) :: resid,  QnetOut
            real(8), intent(inout) :: QnetIn

            real(8), pointer :: Qoverflow, Qstorage, dQdH(:)
            real(8), pointer :: pi

            integer, pointer :: fup(:), fdn(:)

            integer, dimension(1) :: JMar

            integer :: kk
            real(8) :: QratioIn, QratioOut, dQoverflow, Hinc
            real(8) :: MinHeadForOverflow, OverflowDepth, PondedHead

            real(8), dimension(max_branch_per_node) :: dQ, dH, areaQ
            real(8), dimension(1) :: Rarray
            real(8), parameter :: localEpsilon = 1.0d-6
            real(8) :: zbottom, Aout, Ain, AoutOverflow, AinPonded, Astorage
            real(8) :: dQstorage
            integer :: zcount, Qcount, bcount
            logical, dimension(max_branch_per_node) :: bFixYN

            logical :: repeatYN
        !%------------------------------------------------------------------
        !% Alias:
            Qoverflow => elemSR(JMidx,esr_JunctionMain_OverflowRate)
            Qstorage  => elemSR(JMidx,esr_JunctionMain_StorageRate)
            dQdH      => elemSR(:,    esr_JunctionBranch_dQdH)

            fup       => elemI(:,ei_Mface_ul)
            fdn       => elemI(:,ei_Mface_dl)
            pi        => setting%Constant%pi
        !%------------------------------------------------------------------

        repeatYN = .true.

            ! print *, ' '
            ! print *, '============FIX CONSERVATION'
            ! print *, 'JUNCTION ',JMidx
            ! print *, 'starting resid ',resid 

        dQ = zeroR
        dH = zeroR
        dQoverflow = zeroR

        Aout   = zeroR
        Ain    = zeroR

        if ((printJM == JMidx) .and. (Qoverflow .ne. zeroR)) then 
            print *, ' Step ',setting%Time%Step
            print *, 'Qoverflow AAA ',Qoverflow
            print *, 'HEAD          ',elemR(JMidx,er_Head)
        end if

        !% --- note that by definition: QnetIn > 0 and QnetOut < 0 
        !%     resid > 0 implies too much inflow
        !%     resid < 0 implies too much outflow

        !% --- attempt to adjust residual with overflow/ponding alone
        if ((resid < zeroR) .and. (Qoverflow < zeroR)) then 
            !% --- reduce negative magnitude or eliminate negative Qoverflow rate
            !% --- note Qoverflow < 0 is an outflow
            if (Qoverflow < resid) then 
                !% --- resid can be fully accounted for by overflow reduction
                Qoverflow = Qoverflow - resid
                resid = zeroR
                return  !% --- no further residual processing needed
            else
                !% --- resid removes the overflow, but some discrepancy remains
                QnetOut   = QnetOut - Qoverflow
                resid     = resid - Qoverflow  
                Qoverflow = zeroR
                !% -- continue with residual processing
            end if
        elseif ((resid > zeroR) .and. (Qoverflow > zeroR)) then 
            !% --- reduce positive magnitude or elimenate positive Qoverflow rate
            !%     this is an inflow that only occurs due to ponding
            if (Qoverflow > resid) then 
                !% --- resid can be fully accounted for by overflowr eduction
                Qoverflow = Qoverflow - resid
                resid = zeroR
                return !% --- no further residual processing needed
            else
                !% --- resid reduces the overflow, but some resid remains
                QnetIn    = QnetIn - Qoverflow
                resid     = resid - Qoverflow 
                Qoverflow = zeroR
            end if
        else
            !% --- if Qoverflow == 0 no action required
        endif

            ! print *, '0001  Resid    ',JMidx, resid
            ! print *, 'node name ',trim(node%Names(elemI(JMidx,ei_node_Gidx_Bipquick))%str)

            ! print *, 'head ',elemR(JMidx,er_Head)

        !% --- set the conditions for overflow and ponding
        call lljunction_main_overflow_conditions &
            (JMidx, MinHeadForOverflow, OverflowDepth, PondedHead)

            ! print *, MinHeadForOverflow, OverFlowDepth, PondedHead
            if ((printJM == JMidx) .and. (Qoverflow .ne. zeroR)) print *, 'Qoverflow here ',Qoverflow
            
        do while (repeatYN)


            !% --- get flow areas of overflows or ponding
            if (Qoverflow > zeroR) then 
                !% --- inflow from ponding only (OverflowDepth < 0)
                select case (elemSI(JMidx,esi_JunctionMain_OverflowType))
                    case (PondedWeir)
                        AinPonded =  -OverflowDepth * twoR * sqrt( pi * elemSR(JMidx,esr_JunctionMain_Surcharge_Plan_Area))
                    case (PondedOrifice)
                        AinPonded =  -OverFlowDepth * elemSR(JMidx,esr_JunctionMain_OverflowOrifice_Length)
                    case default 
                        print *, 'CODE ERROR: unexpected case default: Overflow with No overflow Type'
                        print *, 'JMidx ',JMIdx, ' ',trim(node%Names(elemI(JMidx,ei_node_Gidx_Bipquick))%str)
                        print *, 'Qoverflow ', Qoverflow 
                        print *, 'Head      ', elemR(JMidx,er_Head)
                        call util_crashpoint(77987233)
                end select
                AoutOverflow = zeroR

            elseif (Qoverflow < zeroR) then
                !% --- outflow either as overflow or to ponding (OverflowDepth > 0)
                select case (elemSI(JMidx,esi_JunctionMain_OverflowType))
                    case (OverflowWeir, PondedWeir)
                        AoutOverflow = OverflowDepth * twoR * sqrt( pi * elemSR(JMidx,esr_JunctionMain_Surcharge_Plan_Area))
                    case (OverflowOrifice, PondedOrifice)
                        AoutOverflow = OverFlowDepth * elemSR(JMidx,esr_JunctionMain_OverflowOrifice_Length)
                    case default 
                        print *, 'CODE ERROR: unexpected case default'
                        call util_crashpoint(7722366)
                end select
                AinPonded = zeroR
            else 
                AinPonded = zeroR 
                AoutOverflow = zeroR
            end if

            !% --- Compute contributing area in branches
            bFixYN = .false. !% --- whether or not a branch can be fixed for conservation
            areaQ  = zeroR   !% --- area for each branch
            bcount = zeroI   !% --- number of modifiable branches
            do kk=1,max_branch_per_node
                !% --- cycle if this branch cannot contribute to flow
                if ((elemSI(JMidx+kk,esi_JunctionBranch_Exists)        .ne. oneI) .or.  &
                    (elemSI(JMidx+kk,esi_JunctionBranch_CanModifyQ)    .ne. oneI)        ) cycle

                if (elemSI(JMidx+kk,esi_JunctionBranch_IsUpstream) == oneI) then 
                    !% --- upstream branch
                    !%     cycle if low junction head or high Fr in branch cannot be adjusted for Q
                    if ((elemR(JMidx   ,er_Head) .le. faceR(fup(JMidx+kk),fr_Zbottom)) .or. &
                        (elemR(JMidx+kk,er_FroudeNumber) .ge. oneR)) cycle
                    !% --- otherwise set bFix to modify this branch
                    bFixYN(kk) = .true.  
                    bcount = bcount + oneI
                    areaQ(kk) = max(faceR(fup(JMidx+kk),fr_Area_d),elemR(JMidx+kk,er_Area))

                        ! print *, 'ff 01 AREA for Flow ',kk, areaQ(kk)

                else
                    !% --- downstream branch
                    !%     low junction head or high -Fr in branch cannot be adjusted for Q
                    if ((elemR(JMidx  ,er_Head) .le. faceR(fdn(JMidx+kk),fr_Zbottom)) .or. &
                        (elemR(JMidx+kk,er_FroudeNumber) .le. -oneR)) cycle
                    !% --- otherwise set bFix to modify this branch
                    bFixYN(kk) = .true.
                    bcount = bcount + oneI
                    areaQ(kk) = max(faceR(fdn(JMidx+kk),fr_Area_u),elemR(JMidx+kk,er_Area))

                        ! print *, 'ff 02 Area for Flow',kk, areaQ(kk)
                end if
                !% --- note, should not reach here unles bFix == .true.

                if (bcount < 1) then 
                    print *, 'CODE ERROR: bcount = 0; should not have reached this point '
                    call util_crashpoint(6209873)
                end if

                    ! print *, 'QVAL ',kk, (real(branchsign(kk),8) * elemR(JMidx+kk,er_Flowrate))
            end do

                ! print *,' '
                ! print *, 'bcount ',bcount
                ! print *,' '

            !% --- if no branches with flow and no ponded inflow, then use all branches
            if ((bcount < 1) .and. (AinPonded == zeroR)) then 
                !% --- occurs during wetting drying, use all real branches for adjustment
                do kk=1,max_branch_per_node
                    if ((elemSI(JMidx+kk,esi_JunctionBranch_Exists)        .ne. oneI) .or.  &
                        (elemSI(JMidx+kk,esi_JunctionBranch_CanModifyQ)    .ne. oneI)        ) cycle

                    if (elemSI(JMidx+kk,esi_JunctionBranch_IsUpstream) == oneI) then
                        !% --- upstream branch
                        bFixYN(kk) = .true.  
                        bcount = bcount + oneI
                        areaQ(kk) = max(faceR(fup(JMidx+kk),fr_Area_d),elemR(JMidx+kk,er_Area))

                        !print *, 'ff 03 ',kk, areaQ(kk)

                    else
                        !% --- downstream branch
                        bFixYN(kk) = .true.
                        bcount = bcount + oneI
                        areaQ(kk) = max(faceR(fdn(JMidx+kk),fr_Area_u),elemR(JMidx+kk,er_Area))

                        !print *, 'ff 04 ',kk, areaQ(kk)

                    end if
                end do

                if (bcount < 1) then 
                    print *, 'CODE ERROR: Mass residual in an unexpected condition'
                    call util_crashpoint(2298522)
                endif
            else
                !% --- continue
            end if

            !% --- Accumulate area weighting for inflows and outflows
            Ain  = Ain  + AinPonded
            Aout = Aout + AoutOverflow
            do kk=1,max_branch_per_node
                if ((elemSI(JMidx+kk,esi_JunctionBranch_Exists)        .ne. oneI) .or.  &
                    (elemSI(JMidx+kk,esi_JunctionBranch_CanModifyQ)    .ne. oneI)        ) cycle
                !% --- accumulate the in and outflow areas for modifiable flows
                if ((real(branchsign(kk),8) * elemR(JMidx+kk,er_Flowrate)) > zeroR) then 
                    Ain  = Ain  + areaQ(kk)
                            ! print *, 'Ain ',kk, Ain
                else
                    Aout = Aout + areaQ(kk)
                            ! print *, 'Aout ',kk, Aout
                end if
            end do

                ! print *, ' '
                ! print *, 'Ain, Aout ',Ain, Aout
                ! print *, ' '

            !% --- area associated with increased/decreased storage
            Astorage = elemSR(JMidx,esr_Storage_Plan_Area)

            if ((Ain < setting%ZeroValue%Area) .and. (Aout < setting%ZeroValue%Area)) then 
                !% degenerate condition 
                print *, 'CODE ERROR: unexpected junction condition '
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

                     ! print *, 'dQ(kk) ',kk, dQ(kk)

                end do
                !% --- overflow rate adjustment
                 ! print *, 'Qoverflow ',Qoverflow
                if (Qoverflow < zeroR) then 
                    !% --- overflow/ponding is an outflow (negative Q)
                    !%     negative residual indicates too much outflow,
                    !%     which requires a positive dQ to reduce the
                    !%     outflow magnitude
                    dQoverflow =  - resid * AoutOverflow / (Ain + Aout + Astorage)

                    !% --- check to see if too much overflow was removed
                    !%     i.e., a positive residual (too much inflow)
                    !%     caused a
                    if (Qoverflow + dQoverflow > zeroR) then 
                        resid = resid - Qoverflow
                        Qoverflow = zeroR
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

                    !% --- check to see if too much ponding was removed
                    !%     i.e., a positive residual (too much inflow)
                    !%     caused a
                    if (Qoverflow + dQoverflow < zeroR) then 
                        resid = resid - Qoverflow
                        Qoverflow = zeroR
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

            ! print *, 'dQ : ', dQ

        end do

        !% --- positive residual (too much inflow) causes increased storage rate
        Qstorage = Qstorage + resid *  Astorage / (Ain + Aout + Astorage)

            ! print *, ' '
            ! print *, 'dQ '
            ! do  kk=1,max_branch_per_node
            !     if (elemSI(JMidx+kk,esi_JunctionBranch_Exists) .ne. oneI) cycle
            !     print *, dQ(kk), elemR(JMidx+kk,er_Flowrate)
            ! end do

        !% --- update flowrates
        do kk=1,max_branch_per_node
            if (elemSI(JMidx+kk,esi_JunctionBranch_Exists) .ne. oneI) cycle
            if (bFixYN(kk) .ne. zeroI) then
                elemR(JMidx+kk,er_Flowrate) = elemR(JMidx+kk,er_Flowrate) + dQ(kk)
            end if
        end do

        !% --- update overflow
        Qoverflow = Qoverflow + dQoverflow
        

        !% --- recompute the residual
        resid = lljunction_conservation_residual (JMidx)

            ! print *, 'NEW RESIDUAL, overflow ',resid, Qoverflow

            ! print *, abs(resid), localEpsilon

        if (abs(resid) > localEpsilon) then 
            print *, 'resid ',resid, ' at junction element ',JMidx
            print *, 'on solution step ',setting%Time%Step
            print *, 'at time ',setting%Time%Now / (3600.d0), ' hours'
            print *, 'CODE ERROR: unexpected flowrate residual'
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

        ! if (JMidx ==printJM) then
        !     print *, ' '
        !     print *, 'IN JUNCTION CONSERVATION RESID FOR ', JMidx
        !     print *, 'Qnetbranches: ',QnetBranches 
        !     print *, 'OverflowRate: ',elemSR(JMidx,esr_JunctionMain_OverflowRate)
        !     print *, 'StorageRate:  ',elemSR(JMidx,esr_JunctionMain_StorageRate)
        !     print *, 'QLat          ',elemR(JMidx,er_FlowrateLateral)
        !     print *, ' '
        ! end if
        
        lljunction_conservation_residual = QnetBranches         &
             + elemSR(JMidx,esr_JunctionMain_OverflowRate)    &
             - elemSR(JMidx,esr_JunctionMain_StorageRate)     &
             + elemR (JMidx,er_FlowrateLateral)

    end function lljunction_conservation_residual
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
            real(8), pointer    :: Lorifice(:), coef2, coef4
            real(8)             :: MinHeadForOverflow, OverFlowDepth, PondedHead
        !%------------------------------------------------------------------  
        !% Aliases
            Lorifice => elemSR(:,esr_JunctionMain_OverflowOrifice_Length)
            coef2    => setting%Junction%Overflow%coef2
            coef4    => setting%Junction%Overflow%coef4
        !%------------------------------------------------------------------  

        !% --- set the conditions for overflow or ponding
        call lljunction_main_overflow_conditions &
            (JMidx, MinHeadForOverflow, OverflowDepth, PondedHead)


        !% --- Head below surcharge, only possibility is a ponding inflow
        !%     But changing Q does not change the inflow
        if (elemR(JMidx,er_Head) .le. MinHeadForOverflow) then 
            lljunction_main_dQdHoverflow = zeroR
            return 
        end if

        !% --- if possible overflow condition exists
        select case (elemSI(JMidx,esi_JunctionMain_OverflowType))  
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
                    * sqrt(elemSR(JMidx,esr_Storage_Plan_Area) * OverflowDepth)
                return

            case (PondedWeir)
                !% --- see explanation for OverflowWeir, above
                !%     abs(OverflowDepth) required so that negative OverflowDepth
                !%     is allowed when ponding is an inflow
                lljunction_main_dQdHoverflow = -coef2                            &
                    * sqrt(elemSR(JMidx,esr_Storage_Plan_Area) * abs(OverflowDepth))
                return
                
            case (OverflowOrifice)
                !% --- use the supplied orifice length
                !%     Using Brater and King Eq. 4.16 or 4.17
                !%     dQ/dH = sqrt(2g) L sqrt(H)
                !%     coef4 = sqrt(2g)
                !%     minus sign as Q is outflow (negative) as H increases
                lljunction_main_dQdHoverflow = -coef4 * Lorifice(JMidx) * sqrt(OverflowDepth)   
                return
                        
            case (PondedOrifice)
                !% --- see explanation for OverflowOrifice, above
                !%     abs(OverflowDepth) required so that negative OverflowDepth
                !%     is allowed when ponding is an inflow
                lljunction_main_dQdHoverflow = -coef4 * Lorifice(JMidx) * sqrt(abs(OverflowDepth))   
                return

            case default
                !% --- should not reach here.
                !%     for debugging, change to impure function and uncomment
                !%     the following
                ! print *, 'CODE ERROR'
                ! print *, 'unexpected case default'
                ! call util_crashpoint(1209874)

        end select 

    end function lljunction_main_dQdHoverflow
    !%    
!%========================================================================== 
!%==========================================================================
!% 
    real(8) function lljunction_main_dQdHstorage (JMidx,iStep)  
        !%------------------------------------------------------------------
        !% Description
        !% Computes the storage rate of junction with tabular or functional
        !% storage
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx, iStep
        !%------------------------------------------------------------------  

        if (elemSI(JMidx,esi_JunctionMain_Type) .eq. NoStorage) then
            !% --- no storage rate for implied storage junctions
            lljunction_main_dQdHstorage = zeroR
            return
        else
            ! if (JMidx == printJM) then
            !     print *, ' '
            !     print *, 'in dQdHstorage ',JMidx, elemYN(JMidx,eYN_isSurcharged)
            !     print *, elemR(JMidx,er_Head), elemR(JMidx,er_Zcrown)
            !     print *, ' '
            ! end if

            if (.not.elemYN(JMidx,eYN_isSurcharged)) then
                !% --- Storage rate term based on time step
                lljunction_main_dQdHstorage = elemSR(JMidx,esr_Storage_Plan_Area) &
                    / (setting%Solver%crk2(iStep) * setting%Time%Hydraulics%Dt)
            else
                !% saz 20230504
                !% --- Storage rate term based on time step
                lljunction_main_dQdHstorage = elemSR(JMidx,esr_JunctionMain_Surcharge_Plan_Area) &
                                            / (setting%Solver%crk2(iStep) * setting%Time%Hydraulics%Dt)
            end if
        end if

    end function lljunction_main_dQdHstorage   
    !%    
!%========================================================================== 
!%==========================================================================
!% 
    subroutine lljunction_main_head_bounds (JMidx, Hbound)
        !%-----------------------------------------------------------------
        !% Description
        !% Computes the minimum head, Hbound(1), and Hbound(2) 
        !%, maximum head for a junction
        !%     Hbound(2) is the maximum head in the surrounding elements
        !%     Hbound(1) is Zbottom of JM, or the lowest Z bottom of any 
        !%          branch if they are all higher than JM
        !%-----------------------------------------------------------------
            integer,               intent(in)    :: JMidx
            real(8), dimension(2), intent(inout) :: Hbound
            integer :: ii, JBidx
            integer, pointer :: fidx
            real(8), pointer :: grav, headJM, headAdj
            real(8), dimension(max_branch_per_node,2) :: HbranchLimit
            logical :: jhead_lowlimit_TF
        !%-----------------------------------------------------------------
            grav => setting%Constant%gravity
            
        !%-----------------------------------------------------------------

            Hbound(1) = +huge(oneR)
            Hbound(2) = -huge(oneR)

            HbranchLimit(:,1) = +huge(oneR)
            HbranchLimit(:,2) = -huge(oneR)

            jhead_lowlimit_TF = .false.

            

            !% FIND ALL BRANCHES WITH INFLOW HEAD THAT COULD LIMIT THe
            !% JUNCTION HEAD
            !% --- cycle through branches (cannot be concurrent)
            !%     fadj* are faces for adjustment, zeroI is null value
            do ii=1,max_branch_per_node
                ! print *, ' '
                ! print *, 'diag adjacent ',JMidx+ii, elemSI(JMidx+ii,esi_JunctionBranch_Diag_adjacent)
                ! print *, ' '
                if (elemSI(JMidx+ii,esi_JunctionBranch_Exists) .ne. oneI) cycle 
                !% --- diagnostic elements cannot be head limiters
                if (elemSI(JMidx+ii,esi_JunctionBranch_Diag_adjacent) .eq. oneI) cycle
                headJM => elemR(JMidx,er_Head)
                if (mod(ii,2)== 0) then 
                    !% --- downstream branch
                    !% --- downstream face
                    fidx => elemI(JMidx+ii,ei_Mface_dL)
                    JBidx = JMidx+ ii
                    headAdj => faceR(fidx,fr_Head_Adjacent)

                    ! print *, ' '
                    ! print *, 'Head Adj ',JBidx, headAdj
                    ! print *, ' '

                    if (elemR(JBidx,er_Flowrate) > zeroR) then
                        !% --- outflow on downstream branch 
                        !%     provides only a lower limit
                        if (headJM > headAdj) then
                            !% --- outflow with positive head gradient uses outside
                            !%     head as low limiter
                            HbranchLimit(ii,1) = headAdj
                        elseif (HeadJM < headAdj) then
                            !% --- outflow with inverse gradient cannot reduce headJM
                            HbranchLimit(ii,1) = headJM
                        else 
                            !% --- zero gradient has no low limiter
                        end if
                    elseif (elemR(JBidx,er_Flowrate) < zeroR) then
                        !% -- inflow on downstream branch
                        !%    provides only an upper limit
                        if (headAdj > headJM) then 
                            !% --- inflow with positive head gradient uses outside
                            !%     head as high limiter
                            HbranchLimit(ii,2) = headAdj
                        elseif (headAdj < headJM) then 
                            !% --- inflow with adverse head gradient cannot increase
                            !%     JM head
                            HbranchLimit(ii,2) = headJM !% THIS SEEMS QUESTIONABLE IN MULTI-BRANCH?
                        else 
                            !% --- zero gradient has no high limiter
                        end if
                    else
                        !% no flowrate
                    end if
                else
                    !% --- upstream branch
                    !% --- upstream face
                    fidx => elemI(JMidx+ii,ei_Mface_uL)
                    JBidx = JMidx+ ii
                    headAdj => faceR(fidx,fr_Head_Adjacent)
                    if (elemR(JBidx,er_Flowrate) < zeroR) then
                        !% --- outflow on upstream branch 
                        !%     provides only a lower limit
                        if (headJM > headAdj) then
                            !% --- outflow with positive head gradient uses outside
                            !%     head as low limiter
                            HbranchLimit(ii,1) = headAdj
                        elseif (HeadJM < headAdj) then 
                            !% --- outflow with inverse gradient cannot reduce headJM
                            HbranchLimit(ii,1) = headJM
                        else 
                            !% --- zero gradient has no low limiter
                        end if
                    elseif (elemR(JBidx,er_Flowrate) > zeroR) then
                        !% -- inflow on upstream branch
                        !%    provides only an upper limit
                        if (headAdj > headJM) then 
                            !% --- inflow with positive head gradient uses outside
                            !%     head as high limiter
                            HbranchLimit(ii,2) = headAdj
                        elseif (headAdj < headJM) then 
                            !% --- inflow with adverse head gradient cannot increase
                            !%     JM head
                            HbranchLimit(ii,2) = headJM
                        else 
                            !% --- zero gradient has no high limiter
                        end if
                    else
                        !% no flowrate
                    end if
                end if
    
            end do

            !% --- select the lowest low limiter from all branches
            Hbound(1) = minval(HbranchLimit(:,1))
            !% --- select the highest limiter from all branches
            Hbound(2) = maxval(HbranchLimit(:,2))

            if (Hbound(1) > 1000.d0) then
                Hbound(1) = -huge(oneR)
            end if
            if (Hbound(2) < -1000.d0) then 
                Hbound(2) = + huge(oneR)
            end if


    end subroutine lljunction_main_head_bounds    
!%    
!%==========================================================================
!%==========================================================================
!% 
    real(8) function lljunction_main_Qoverflow (JMidx, istep)  
        !%------------------------------------------------------------------
        !% Description
        !% Computes the overflow/ponding rate of a junction that is stored
        !% Negative is outflow, positive is inflow
        !% Weir overflow based on non-dimensional approach of Brater and King
        !% for broad-crested weir
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx, istep
            real(8), pointer    :: Horifice(:), Lorifice(:), coef1, coef3, WeirFactor
            real(8)             :: MinHeadForOverflow, OverflowDepth, PondedHead
            real(8)             :: tempOverflow

            real(8), pointer    :: dt, crk
        !%------------------------------------------------------------------  
        !% Aliases
            Horifice   => elemSR(:,esr_JunctionMain_OverflowOrifice_Height)
            Lorifice   => elemSR(:,esr_JunctionMain_OverflowOrifice_Length)
            WeirFactor => setting%Junction%Overflow%WeirLengthFactor
            coef1      => setting%Junction%Overflow%coef1
            coef3      => setting%Junction%Overflow%coef3
            dt         => setting%Time%Hydraulics%Dt
            crk        => setting%Solver%crk2(istep)
        !%------------------------------------------------------------------  

        !% --- set the overflow conditions
        call lljunction_main_overflow_conditions &
            (JMidx,MinHeadForOverflow, OverflowDepth, PondedHead)

            ! if (JMidx==printJM) print *, ' in LLjunction head ',elemR(JMidx,er_Head), elemR(JMidx,er_Zcrown)

        !% --- Head below the surcharge, only possibility is a ponding inflow
        if (elemR(JMidx,er_Head) .le. MinHeadForOverflow) then
            !% --- no outflowing overflow 
            if ( elemSR(JMidx,esr_JunctionMain_PondedArea) == zeroR) then 
                !% --- no ponding inflow or outflow
                lljunction_main_Qoverflow = zeroR
                return

            else
                !% --- possible ponding inflow as waterfall into junction
                if (PondedHead > MinHeadForOverflow) then 
                    !% --- ponding inflow as waterfall
                    select case (elemSI(JMidx,esi_JunctionMain_OverflowType))

                        case (PondedWeir)
                            !% --- see OverflowWeir below for explanation
                            !%     given + value as inflow
                            tempOverflow = coef1 *  WeirFactor * sqrt(elemSR(JMidx,esr_Storage_Plan_Area)) &
                                * ((PondedHead - MinHeadForOverflow)**threehalfR)

            
                        case (PondedOrifice)
                            !% ---- see OverflowOrifice below for explanation
                            !%      given + value as inflow
                            tempOverflow = coef3 * Lorifice(JMidx) &
                                    * ((PondedHead - MinHeadForOverflow)**threehalfR)

                        case default 
                            print *, 'CODE ERROR: unexpected case default'
                            call util_crashpoint(7298723)
                    end select

                    !% --- limit inflow by ponded water available
                    if ((tempOverflow * dt *crk)  > elemR(Jmidx,er_VolumePondedTotal)) then
                        lljunction_main_Qoverflow = elemR(Jmidx,er_VolumePondedTotal) / (dt * crk)
                    else 
                        lljunction_main_Qoverflow = tempOverflow
                    end if
                    return

                else
                    !% --- no ponding inflow or outflow
                    lljunction_main_Qoverflow = zeroR
                    return

                end if
            end if
        else 
            !% --- continue
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
        !%     (1) a flowrate is positive if it is an outflow
        !%      and negative if it is an inflow.
        !%     (2) a positive Overflow Depth causes an outflow (negative Q)
        !%      whereas a negative Overflow depth causes a inflow (positive Q)
        !%     (3) A positive inflow is only allowed with ponding.
        !%     
        select case (elemSI(JMidx,esi_JunctionMain_OverflowType))

            case (NoOverflow)
                !% --- NoOverflow implies no ponding or overflow
                lljunction_main_Qoverflow = zeroR
                return

            case (OverflowWeir)
                !% --- weir overflow based on weir length as circumference of storage plan area
                !%     Q = C L H^{3/2}
                !%     L = 2 (pi A)^{1/2}
                !%     1.38 < C < 1.83 from Brater and King Table 5.1 where L in m
                !%     and C in m^{1/2} / S
                !%     Q = coef1 (sqrt(A)) H^(3/2)
                !%     minus sign as Q is outflow (negative)
                lljunction_main_Qoverflow                                &
                    = -coef1 * WeirFactor * sqrt(elemSR(JMidx,esr_Storage_Plan_Area)) &
                       * (OverFlowDepth**threehalfR)
                return

            case (PondedWeir)           
                !% --- similar to OverflowWeir above, but allows in or outflow 
                if (OverflowDepth < zeroR) then 
                    !% --- inflow from ponding (+ value)
                    tempOverflow                                             &
                        = +coef1 * WeirFactor * sqrt(elemSR(JMidx,esr_Storage_Plan_Area)) &
                            * (-OverflowDepth)**threehalfR

                    !% --- limit inflow by ponded water available
                    if ((tempOverflow * dt *crk)  > elemR(Jmidx,er_VolumePondedTotal)) then
                        lljunction_main_Qoverflow = elemR(Jmidx,er_VolumePondedTotal) / (dt * crk)
                    else 
                        lljunction_main_Qoverflow = tempOverflow
                    end if

                    return

                elseif (OverflowDepth > zeroR) then 
                    !% --- outflow to ponding (- value)
                    lljunction_main_Qoverflow                                &
                        = -coef1 * WeirFactor * sqrt(elemSR(JMidx,esr_Storage_Plan_Area)) &
                            * (OverflowDepth**threehalfR)
                    return

                else 
                    !% --- no exchange
                    lljunction_main_Qoverflow = zeroR
                    return
                end if

            case (OverflowOrifice)
                !% --- orifice overflow assuming a single orifice
                !%     Using Brater and King Eq. 4.17 
                !%     Q = (2/3) L sqrt(2g) H^(3/2)
                !%     define coef3 = (2/3) sqrt(2g) so that
                !%     Q = coef3 * L * H^{3/2}
                !%     minus sign as Q is outflow (negative)
                if (elemR(JMidx,er_Head) .le. (MinHeadForOverflow + Horifice(JMidx) )) then
                    !% --- head (water surface) below upper edge of orifice (H < Z + Horifice)
                    lljunction_main_Qoverflow = -coef3 * Lorifice(JMidx) * (OverFlowDepth**threehalfR)
                    return 

                else
                    !% --- orifice is pressurized (H > Zcrown + Horifice)
                    !%     Using Brater and King eq. 4.16
                    !%     Q = (2/3) L sqrt(2g) ( (H-Zcrown)^(3/2) - (H-(zcrown+Horifice))^(3/2) )
                    !%     define coef3 = (2/3) sqrt(2g) so that
                    !%     Q = coef3 * L * ( (H-Zcrown)^(3/2) - (H-(zcrown+Horifice))^(3/2) )
                    !%     minus sign as Q is outflow (negative) 
                    lljunction_main_Qoverflow = -coef3 * Lorifice(JMidx)  &
                        * (                                                                                 &
                            +  ((OverFlowDepth                  )**threehalfR)   &
                            -  ((OverFlowDepth - Horifice(JMidx))**threehalfR)    &
                          )
                    return
                end if
             
            case (PondedOrifice)
                !% --- similar to OverflowOrifice but allows in or outflows with ponding
                if (OverflowDepth < zeroR) then
                    !% --- inflows from ponding require + value
                    if (PondedHead .le. (MinHeadForOverflow + Horifice(JMidx) )) then
                        !% --- water surface below upper edge of orifice (H < Z + Horifice)
                        tempOverflow = +coef3 * Lorifice(JMidx) * ((-OverflowDepth)**threehalfR)

                    else
                        !% --- inflows from ponding
                        tempOverflow = +coef3 * Lorifice(JMidx)  &
                            * (                                                                                 &
                                +  ((-OverFlowDepth                  )**threehalfR)   &
                                -  ((-OverflowDepth - Horifice(JMidx))**threehalfR)    &
                            )
                    end if

                    !% --- limit inflow by ponded water available
                    if ((tempOverflow * dt *crk)  > elemR(Jmidx,er_VolumePondedTotal)) then
                        lljunction_main_Qoverflow = elemR(Jmidx,er_VolumePondedTotal) / (dt * crk)
                    else 
                        lljunction_main_Qoverflow = tempOverflow
                    end if
                    return 

                elseif (OverflowDepth > zeroR) then
                    !% --- outflows to ponding with - value
                    if (PondedHead .le. (MinHeadForOverflow + Horifice(JMidx) )) then
                        !% --- water surface below upper edge of orifice (H < Z + Horifice)
                        lljunction_main_Qoverflow = -coef3 * Lorifice(JMidx) &
                            * ((OverflowDepth)**threehalfR)
                        return 

                    else
                        lljunction_main_Qoverflow = -coef3 * Lorifice(JMidx)  &
                            * (                                                                                 &
                                +  ((OverFlowDepth                  )**threehalfR)   &
                                -  ((OverflowDepth - Horifice(JMidx))**threehalfR)    &
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
                print *, 'CODE ERROR'
                print *, 'unexpected case default'
                call util_crashpoint(6209874)
        end select

    end function lljunction_main_Qoverflow   
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
            integer :: k1, k2,kk
        !%------------------------------------------------------------------
        !k1 = JMidx + 1
        !k2 = JMidx + max_branch_per_node

        thissum = zeroR

        ! print *, ' '
        ! print *, 'junction_main_sumBranches ', JMidx
        ! print *, branchsign
        ! print *, ' '
        ! print *, thisArray(k1:k2,thisCol)
        ! print *, ' '

        !junction_main_sumBranches = sum(branchsign * thisArray(k1:k2,thisCol)) 
        !junction_main_sumBranches = sum( branchsign                                       &
        !                                * real(elemSI(k1:k2,esi_JunctionBranch_Exists),8) &
        !                                * thisArray(k1:k2,thisCol) )

        do kk = 1,max_branch_per_node 
            if (elemSI(JMidx+kk,esi_JunctionBranch_Exists) .ne. oneI) cycle
            thissum = thissum + branchsign(kk) * thisArray(JMidx+kk,thisCol)
            !print *, 'thissum ',kk, thissum
        end do

        lljunction_main_sumBranches = thissum

    end function lljunction_main_sumBranches   
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

            real(8) :: QnetBranches
        !%------------------------------------------------------------------
        !% Alias
            Qstorage  => elemSR(JMidx,esr_JunctionMain_StorageRate)
            Qoverflow => elemSR(JMidx,esr_JunctionMain_OverflowRate)
            dt        => setting%Time%Hydraulics%Dt
            crk       => setting%Solver%crk2(istep)
        !%------------------------------------------------------------------

        !% --- update storage volume
        if (elemSI(JMidx,esi_JunctionMain_Type) .ne. NoStorage) then
            elemR(JMidx,er_Volume) = elemR(JMidx,er_Volume_N0)                         &
                + dt *crk * Qstorage
        else 
            !% --- no storage does not change volume
            elemR(JMidx,er_Volume) = elemR(JMidx,er_Volume_N0)
        end if

        !% --- update this time step volume overflow and ponded
        select case (elemSI(JMidx,esi_JunctionMain_OverflowType))
            case (NoOverflow)
                !% -- no action

            case (OverflowWeir,OverflowOrifice)
                !% --- define the VolumeOverflow as positive
                elemR(JMidx,er_VolumeOverFlow) = -Qoverflow * dt * crk

            case (PondedWeir, PondedOrifice)
                !% --- define the VolumePonded as positive is increasing 
                !%     the ponded volume and negative decreasing it.
                !%     This changes the sign of Qoverflow, which is 
                !%     negative for an outflow rate
                elemR(JMidx,er_VolumePonded)   = -Qoverflow * dt * crk

            case default
                print *, 'CODE ERROR: unexpected case default'
                call util_crashpoint(711345)
        end select

    end subroutine lljunction_main_update_Qdependent_values
!%    
!%==========================================================================
!%==========================================================================
!% 
    real(8) function lljunction_main_update_storage_rate (JMidx, dH, QnetBranches, istep)
        !%------------------------------------------------------------------
        !% Description:
        !% Updates the Qstorage for change in head dH
        !%------------------------------------------------------------------
            integer, intent(in) :: JMidx, istep
            real(8), intent(in) :: dH, QnetBranches
            integer :: ii
        !%------------------------------------------------------------------
        select case (elemSI(JMidx,esi_JunctionMain_Type))

            case (NoStorage)
                lljunction_main_update_storage_rate = zeroR
                return
            case (TabularStorage,FunctionalStorage,ImpliedStorage)
                !% --- compute the storage flowrate
                if (istep == 1) then
                    !% --- compute rate from dH
                    if (.not. elemYN(JMidx,eYN_isSurcharged)) then
                        lljunction_main_update_storage_rate                                        &
                            = (elemSR(JMidx,esr_Storage_Plan_Area) * dH)                    &
                                / (setting%Solver%crk2(istep) * setting%Time%Hydraulics%Dt)
                    else
                        lljunction_main_update_storage_rate                                        &
                            = (elemSR(JMidx,esr_JunctionMain_Surcharge_Plan_Area) * dH)                    &
                                / (setting%Solver%crk2(istep) * setting%Time%Hydraulics%Dt)
                    end if
                    return
                elseif (istep == 2) then 
                    !% --- compute rate from mass conservation
                    lljunction_main_update_storage_rate = QnetBranches       &
                        + elemSR(JMidx,esr_JunctionMain_OverflowRate) &
                        + elemR(JMIdx,er_FlowrateLateral)
                else 
                    !% should not be possible
                    print *, 'CODE ERROR: unexpected else'
                    call util_crashpoint(6209874)
                    return
                end if

            case default
                lljunction_main_update_storage_rate = zeroR
                print *, 'CODE ERROR: unexpected case default'
                call util_crashpoint(8852783)
                return

        end select

    end function lljunction_main_update_storage_rate   
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
            integer :: mm, ii, kk, JMidx, nInflow

            real(8), dimension(max_branch_per_node) :: inVelocity, inFlow
            real(8) :: weightedVelocity, junctionVelocity, totalQ
            real(8) :: inletVelocity
        !%-----------------------------------------------------------------
        !% Preliminaries
            Npack => npack_elemP(thisColP)
            if (Npack < 1) return
        !%-----------------------------------------------------------------

        do mm=1,Npack 
            JMidx = elemP(mm,thisColP)

            elemR(JMidx,er_Flowrate) = zeroR 
            elemR(JMidx,er_Velocity) = zeroR

            ! ! print *, ' '
            ! ! print *, '============================'
            ! ! print *, 'in junction main velocity '
            ! ! print *, 'JMidx = ',JMidx

            inFlow     = zeroR
            inVelocity = zeroR

            !% --- cycle thru branches to get inflow velocities and flowrates
            !%     note that the net inflow velocity is directional (i.e. an inflow
            !%     from downstream against the flow cancels velocity of inflows
            !%     from upstream with the flow)
            kk=1
            do ii=1,max_branch_per_node
                if (elemSI(JMidx+ii,esi_JunctionBranch_Exists) .ne. oneI)  cycle
                if (elemR(JMidx+ii,er_Depth) .le. setting%ZeroValue%Depth) cycle
                if (mod(ii,2) == 0) then 
                    !% --- downstream branch
                    if (elemR(JMidx+ii,er_Flowrate) < zeroR) then 
                        !% --- downstream branch inflow
                        inVelocity(kk) = elemR(JMidx+ii,er_Velocity)
                        inFlow(kk)     = abs(elemR(JMidx+ii,er_Flowrate))
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
                totalQ = sum(inFlow(1:nInflow))
            else
                totalQ = zeroR
            end if

            !% --- add lateral inflow
            if (elemR(JMidx,er_FlowrateLateral) > zeroR) then
                totalQ = totalQ + elemR(JMidx,er_FlowrateLateral)
            else 
                !% --- nothing to add 
            end if


            if (totalQ > zeroR) then 
                !% --- weightedvelocity is the flow-weighted average velocity of the inflows
                !% --- HACK we need some branch reduction factors for more than 1 inflow branch
                weightedVelocity = sum(inFlow(1:nInflow) * inVelocity(1:nInflow) ) / totalQ

                !% --- junctionVelocity is the velocity for the flow given the junction geometry
                junctionVelocity = totalQ / ( sqrt(elemSR(JMidx,esr_Storage_Plan_Area)) * elemR(JMidx,er_Depth))

                if (weightedVelocity == zeroR) then
                    inletVelocity = junctionVelocity
                elseif (weightedVelocity < zeroR) then 
                    !% -- net upstream velocity
                    inletVelocity = max(weightedVelocity, -junctionVelocity)
                else 
                    inletVelocity = min(weightedVelocity,junctionVelocity)
                end if

                elemR(JMidx,er_Velocity) = inletVelocity

                if (abs(inletVelocity) > setting%Limiter%Velocity%Maximum) then
                    inletVelocity = sign(0.99d0 * setting%Limiter%Velocity%Maximum, inletVelocity) 
                end if

                elemR(JMidx,er_Flowrate)     = inletVelocity * elemR(JMidx,er_Depth) * sqrt(elemSR(JMidx,esr_Storage_Plan_Area))
                elemR(JMidx,er_FroudeNumber) = inletVelocity / sqrt(setting%Constant%gravity * elemR(JMidx,er_Depth))

            else
                !% --- no inflows to provide momentum, 
            end if
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
            Qstorage  => elemSR(JMidx,esr_JunctionMain_StorageRate)
            dt        => setting%Time%Hydraulics%Dt
            crk       => setting%Solver%crk2
        !%------------------------------------------------------------------

        if (elemSI(JMidx,esi_JunctionMain_Type) .ne. NoStorage) then
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
    subroutine lljunction_push_inflowCC_flowrates_to_face ()
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

                    !% --- adjust velocity where face area is sufficient
                    !% HACK --- THIS MAY NEED TO BE COMMENTED OUT 20230506brh
                    !% IT IS POSSIBLE THAT THESE VELOCITIES COULD BE TOO LARGE.
                    !% NEED TO TEST WHETHER IT IS BETTER TO PUSH THE ELEMENT VELOCITY 
                    !% TO THE FACE.
                    ! where (fr_Area_u > setting%ZeroValue%Area)
                    !     faceR(fidx(thisP),fr_Velocity_u) = elemR(thisP,er_Flowrate) / faceR(fidx(thisP),fr_Area_u)
                    ! endwhere
                    ! where (fr_Area_d > setting%ZeroValue%Area)
                    !     faceR(fidx(thisP),fr_Velocity_d) = elemR(thisP,er_Flowrate) / faceR(fidx(thisP),fr_Area_d)
                    ! endwhere

                    !% --- Note where a shared face has been diverged 
                    !%    from value on other image
                    where (faceYN(fidx(thisP),fYN_isSharedFace))
                        faceYN(fidx(thisP),fYN_isSharedFaceDiverged) = .true.
                    end where

                endwhere
            end if
        end do

        ! !% --- upstream CC elements
        ! Npack => npack_elemP(ep_CC_UpstreamOfJunction)
        ! if (Npack > 0) then 
        !     !% --- set of CC elements upstream of a JB
        !     thisP => elemP(1:Npack,ep_CC_UpstreamOfJunction)
        !     !% --- face indx for CC/JB face downstream of CC
        !     fidx  => elemI(:,ei_Mface_dL)
        !     where (elemR(thisP,er_Flowrate) > zeroR)
        !         !% --- inflows into junction
        !         faceR(fidx(thisP),fr_Flowrate)   = elemR(thisP,er_Flowrate)
        !         where (fr_Area_u > setting%ZeroValue%Area)
        !             faceR(fidx(thisP),fr_Velocity_u) = elemR(thisP,er_Flowrate) / faceR(fidx(thisP),fr_Area_u)
        !         endwhere
        !         where (fr_Area_d > setting%ZeroValue%Area)
        !             faceR(fidx(thisP),fr_Velocity_d) = elemR(thisP,er_Flowrate) / faceR(fidx(thisP),fr_Area_d)
        !         endwhere
        !         !% --- Note where a shared face has been diverged 
        !         !%    from value on other image
        !         where (faceYN(fidx(thisP),fYN_isSharedFace))
        !             faceYN(fidx(thisP),fYN_isSharedFaceDiverged) = .true.
        !         end where
        !     endwhere
        ! end if
        
        ! !% --- downstream CC elements
        ! Npack => npack_elemP(ep_CC_DownstreamOfJunction)
        ! if (Npack > 0) then 
        !     !% --- the set of CC elemenst downstream of a JB
        !     thisP => elemP(1:Npack,ep_CC_DownstreamOfJunction)
        !     fidx  => elemI(:,ei_Mface_uL)
        !     where (elemR(thisP,er_Flowrate) < zeroR)
        !         !% --- inflows into the junction
        !         faceR(fidx(thisP),fr_Flowrate)   = elemR(thisP,er_Flowrate)
        !         where (fr_Area_u > setting%ZeroValue%Area)
        !             faceR(fidx(thisP),fr_Velocity_u) = elemR(thisP,er_Flowrate) / faceR(fidx(thisP),fr_Area_u)
        !         endwhere 
        !         where (fr_Area_d > setting%ZeroValue%Area)
        !             faceR(fidx(thisP),fr_Velocity_d) = elemR(thisP,er_Flowrate) / faceR(fidx(thisP),fr_Area_d)
        !         endwhere
        !         !% check if any shared face has been diverged
        !         where (faceYN(fidx(thisP),fYN_isSharedFace))
        !             faceYN(fidx(thisP),fYN_isSharedFaceDiverged) = .true.
        !         end where
        !     endwhere
        ! end if

        ! ! ! ! 


    end subroutine lljunction_push_inflowCC_flowrates_to_face
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine lljunction_push_adjacent_elemdata_to_face ()
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
        !%-----------------------------------------------------------------
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
            call face_push_elemdata_to_face (epCCcol, fr_Head_Adjacent,     er_Head,         elemR, isUpstreamface)
            call face_push_elemdata_to_face (epCCcol, fr_Topwidth_Adjacent, er_Topwidth,     elemR, isUpstreamface)
            call face_push_elemdata_to_face (epCCcol, fr_Length_Adjacent,   er_Length,       elemR, isUpstreamface)
            call face_push_elemdata_to_face (epCCcol, fr_Zcrest_Adjacent,   er_Zbottom,      elemR, isUpstreamface)
            call face_push_elemdata_to_face (epCCcol, fr_Velocity_Adjacent, er_Velocity,     elemR, isUpstreamface)
            call face_push_elemdata_to_face (epCCcol, fr_Froude_Adjacent,   er_FroudeNumber, elemR, isUpstreamface)
            call face_push_elemdata_to_face (epCCcol, fr_Depth_Adjacent,    er_Depth,        elemR, isUpstreamface)
        end do

    end subroutine lljunction_push_adjacent_elemdata_to_face
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine lljunction_main_overflow_conditions &
            (JMidx, MinHeadForOverflow, OverflowDepth, PondedHead)
        !%-----------------------------------------------------------------
        !% Description
        !% Computes data needed for evaluating conditions for overflow
        !% and ponding
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in)    :: JMidx
            real(8), intent(inout) :: MinHeadForOverFlow, OverFlowDepth, PondedHead
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------

        !% --- minimum head at which overflow/poinding occurs
        MinHeadForOverflow = elemR(JMidx,er_Zcrown) &
             + elemSR(JMidx,esr_JunctionMain_OverflowHeightAboveCrown)

        !% --- ponding junctions are treated differently than non-ponding
        if ( elemSR(JMidx,esr_JunctionMain_PondedArea) > zeroR) then 
            !% --- ponded head is the head available in the ponded area
            PondedHead       = MinHeadForOverflow                           &
                                + (  elemR (JMidx,er_VolumePondedTotal)   &
                                    / elemSR(JMidx,esr_JunctionMain_PondedArea)  )
            !% --- for ponding, the overflow depth is the difference                        
            !%     between the junction head and ponded head. A positive value
            !%     causes an outflow to ponding whereas a negative value
            !%     causes an inflow from ponding to the junciton
            OverFlowDepth = elemR(JMidx,er_Head) - PondedHead
        else 
            !% --- by definition, the ponded head at an overflow is the minimum
            !%     head for overflow -- but this should never be needed
            PondedHead = MinHeadForOverflow 
            !% --- For non-ponding, the overflow depth is either positive or zero
            OverFlowDepth = max(elemR(JMidx,er_Head) - MinHeadForOverflow, zeroR)                   
        end if

    end subroutine lljunction_main_overflow_conditions
!%
!%==========================================================================

!%==========================================================================
!%
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------

!%
!%==========================================================================
!%==========================================================================
!% 
!%========================================================================== 
!% END OF MODULE
!%==========================================================================
end module lowerlevel_junction