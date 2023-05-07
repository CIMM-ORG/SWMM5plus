module lowerlevel_junction

    use define_globals
    use define_keys
    use define_indexes
    use define_xsect_tables
    use define_settings, only: setting
    use face, only: face_push_elemdata_to_face
    use update, only: update_Froude_number_element, update_wavespeed_element
    use utility_crash, only: util_crashpoint

!%----------------------------------------------------------------------------- 
!% Description:
!% Procedures called from withing junction_elements module
!%----------------------------------------------------------------------------- 

    implicit none

    private

    public :: lljunction_main_velocity
    public :: lljunction_push_inflowCC_flowrates_to_face
    public :: lljunction_push_adjacent_elemdata_to_face
    public :: lljunction_branch_energy_outflow 
    public :: lljunction_branch_dQdH


    contains
!%==========================================================================
!% PUBLIC
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

            integer :: frHead, frArea, ii
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
                end if
                !% --- the JB-adjacent face is upstream
                fidx => elemI(:,ei_Mface_uL)
                frHead = fr_Head_u !% -- use u for jump purposes? QUESTION
                frArea = fr_Area_u !% QUESTION
            else
                !% --- downstream JB
                bsign = -oneR
                isUpstreamBranch = .false.
                Npack => npack_elemP(ep_JB_Downstream_CC_Adjacent)
                if (Npack > 0) then 
                    thisJB => elemP(1:Npack,ep_JB_Downstream_CC_Adjacent)
                end if
                !% --- the JB-adjacent face is downstream
                fidx => elemI(:,ei_Mface_dL)
                frHead = fr_Head_d !% --- use d for jump purposes? QUESTION
                frArea = fr_Area_d !% QUESTION
            end if

            !% --- head difference from junction main to face
            deltaH => elemR(:,er_Temp01)
            deltaH = zeroR
            where (elemR(thisJB,er_Depth) > setting%ZeroValue%Depth)
                deltaH(thisJB) =  elemR(JMar(thisJB),er_Head) - faceR(fidx(thisJB),frHead)
            endwhere

            ! if (.not. isUpstreamBranch) then 
            !     print *, 'deltaH(6)',deltaH(6), elemR(6,er_Depth)
            !     print *, 'JM        ', JMar(6), elemR(JMar(6),er_Head)
            !     print *, 'fidx      ',fidx(6), faceR(fidx(6),frHead)
            ! end if

            !% --- limit dH if it is more than 1/2 of depth in branch
            where (deltaH(thisJB) > onehalfR   * elemR(thisJB,er_Depth))
                deltaH(thisJB)    = onefourthR * elemR(thisJB,er_Depth)
            end where

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