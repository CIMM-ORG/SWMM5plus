module junction_elements

    use define_globals
    use define_keys
    use define_indexes
    use define_xsect_tables
    use define_settings, only: setting
    use adjust
    use face, only: face_interpolation, face_force_JBadjacent_values
    use pack_mask_arrays, only: pack_small_and_zero_depth_elements, pack_zero_depth_interior_faces
    use update
    ! use utility_unit_testing, only: util_utest_CLprint
    use utility_crash, only: util_crashpoint

!%----------------------------------------------------------------------------- 
!% Description:
!% Computes junction elements
!%----------------------------------------------------------------------------- 

    implicit none

    private
    
    !public :: junction_toplevel
    public :: junction_toplevel_3
    public :: junction_force_Qinterpweights
    public :: junction_face_terms
    public :: junction_mass_conservation

    real(8), parameter :: Cbc = 0.46295d0  !% HACK dimensionless 1.45/sqrt(g) from Brater and King broad-crested weir coefficient
    real(8), parameter :: Horifice = 0.15d0  !% HACK fixed orifice height
    real(8), parameter :: Lorifice = 1.5d0   !% HACK fixed orifice length

    real(8) :: coef1, coef2, coef3, coef4

    integer :: printJM = 11 !% testing

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine junction_toplevel_3 (istep)
        !%-----------------------------------------------------------------
        !% Description:
        !% Controls computation of implicit junction element
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: istep
            integer, pointer    :: Npack, thisP(:)
            integer :: ii
        !%-----------------------------------------------------------------
        !% Preliminaries
            if (npack_elemP(ep_JM) == 0) return

            !% --- coefficients in orifice and weir eqquations for overflow
            coef1 = twoR * Cbc * sqrt(setting%Constant%gravity * setting%Constant%pi)
            coef2 = threehalfR * coef1
            coef3 = twothirdR  * sqrt(twoR * setting%Constant%gravity)
            coef4 = threehalfR * coef3        
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------
        !% --- Ensures JB interpweights are based on flow on both sides of face
        !%     This is needed so that when Q(JB) = 0 the adjacent element Q
        !%     will be interpolated to the intervening face
            
        Npack => npack_elemP(ep_JB)  
        thisP => elemP(1:Npack,ep_JB)
        call update_interpweights_JB (thisP, Npack, .false.)

            ! ! call util_utest_CLprint ('------- DDD08 after update_interpweights_JB')        

        !% --- interpolate latest flowrate (only) to face around JB
        !%     Qyn only, skipJump and skipZeroAdjust
        sync all
        call face_interpolation(fp_JB, .false., .false., .true., .true., .true.) 

            ! ! call util_utest_CLprint ('------- DDD09 after face_interpolation for JB') 

        !% --- computes JM values for head, VolumeOverflow, Volume, 
        !      and JB values  for head, dQdH, DeltaQ, flowrate, StorageRate, OverflowRate
        Npack => npack_elemP(ep_JM)
        thisP => elemP(1:Npack,ep_JM)
        call junction_calculation_3 (thisP, Npack, istep)

            ! ! call util_utest_CLprint ('------- DDD10 after junction_calculation') 

        !% --- Update the JM and JB auxiliary variables 
        !%     For JB this is geometry, Fr, C, interpwieghts
        !%     For JM this is depth, storage_plan_area, elldepth
        !%     .true. = force interp weights for Q to favor JB
        call update_auxiliary_variables_JMJB ( .true.)

            ! ! call util_utest_CLprint ('------- DDD11 after update_auxiliary_variables_JMJB') 

        call adjust_element_toplevel (JB)
        call adjust_element_toplevel (JM)
       

        call face_force_JBadjacent_values (ep_JM, .true.)
        call face_force_JBadjacent_values (ep_JM, .false.)

        ! ! call util_utest_CLprint ('------- FFF07  after face_force_JB... in junction')

        if (num_images() > oneI) then 
            print *, 'CODE ERROR: NEED UPDATE FOR FACE CHANGES DURING JB'
            stop 598743
        end if

        !% --- Adjust JB adjacent elements using new face values
        !%     Fixes flowrate, volume, and geometry
        call junction_CC_for_JBadjacent (ep_CC_UpstreamOfJunction,   istep, .true.)
        call junction_CC_for_JBadjacent (ep_CC_DownstreamOfJunction, istep, .false.)

            ! ! call util_utest_CLprint ('------- FFF08  after junction_CC_forJBadjacent')

        call adjust_face_toplevel(fp_JB)

    end subroutine junction_toplevel_3
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine junction_face_terms (thisColP, istep)
        !%------------------------------------------------------------------
        !% Description:
        !% computes the fa and fb terms from time n data for junction solution
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisColP   !% packed column of junction faces
            integer, intent(in) :: istep      !% rk2 substep
            integer, pointer    :: Npack, fidx, eUp, eDn

            real(8), pointer    :: crk, dt, grav
            real(8), pointer    :: fLengthAdj(:), fTopAdj(:), fHeadAdj(:)
            real(8), pointer    :: fHeadU(:), fHeadD(:), fAreaU(:), fAreaD(:)
            real(8), pointer    :: fVelU(:), fVelD(:)

            integer             :: mm
        !%------------------------------------------------------------------
        !% Aliases
            crk => setting%Solver%crk2(istep)
            dt  => setting%Time%Hydraulics%Dt
            grav=> setting%Constant%gravity

            fLengthAdj => faceR(:,fr_Length_Adjacent)
            fTopAdj    => faceR(:,fr_Topwidth_Adjacent)
            fHeadAdj   => faceR(:,fr_Head_Adjacent)
            fHeadU     => faceR(:,fr_Head_u)
            fHeadD     => faceR(:,fr_Head_d)
            fAreaU     => faceR(:,fr_Area_u)
            fAreaD     => faceR(:,fr_Area_d)
            fVelU      => faceR(:,fr_Velocity_u)
            fVelD      => faceR(:,fr_Velocity_u)
        !%------------------------------------------------------------------
        !% Preliminaries
            Npack => npack_faceP(thisColP)
            if (Npack == 0) return  !% no junctions
        !%------------------------------------------------------------------

        if (setting%Solver%MomentumSourceMethod .ne. T00) then 
            print *, 'CONFIGURATION ERROR'
            print *, 'SWMM5+ presently requires setting.Solver.MomentumSourceMethod = T00'
            print *, 'However, method ',trim(reverseKey(setting%Solver%MomentumSourceMethod)),' is used'
            call util_crashpoint(6698734)
        end if
        
            ! print *, 'in junction face terms'
        !% --- fa and fb terms for the JB elements adjacent to a JB face
        do mm=1,Npack

            !     print *, ' '
            !     print *, 'in junction_face_terms'
            !     print *, 'mm = ',mm

            fIdx => faceP(mm,thisColP)

            ! print *, 'fidx ',fidx

            eDn  => faceI(fIdx,fi_Melem_dL)
            eUp  => faceI(fIdx,fi_Melem_uL)

            !print *, 'fidx, eDn, eUp', fIdx, eDn, eUp

            !% --- error checking
            if ((elemI(eDn,ei_elementType) == JB) .and. (elemI(eUp,ei_elementType) == JB)) then 
                print *, 'CONFIGURATION ERROR: system has a junction/junction connection'
                print *, 'without a conduit, channel, weir, orifice, or pump between.'
                print *, 'This configuration is not allowed in SWMM5+'
                call util_crashpoint(639874)
            end if

            !% --- handle upstream branch (downstream JB element )
            if (elemI(eDn,ei_elementType) == JB) then
                !% --- for different upstream element types
                select case (elemI(eUp,ei_elementType))
                    case (CC)
                            !print *, 'is CC'
                        !% --- upstream conduit/channel, downstream JB
                        !% --- fa term
                        elemSR(eDn,esr_JunctionBranch_fa)                                    &
                            = - (crk * dt / fLengthAdj(fidx))                                &
                            * (                                                              &
                                + (fVelU(fidx)**twoI) * fTopAdj(fidx)                        &
                                + grav * (fAreaU(fidx)                                       &
                                            + fTopAdj(fidx) * (fHeadU(fidx) - fHeadAdj(fidx))) &
                                ) 

                            !print *, ' fa dn ',  elemSR(eDn,esr_JunctionBranch_fb)

                        !% --- fb term    
                        elemSR(eDn,esr_JunctionBranch_fb)                                     &
                            =  -crk * dt * grav * fTopAdj(fidx) / fLengthAdj(fidx)  

                            !print *, ' fb dn ',  elemSR(eDn,esr_JunctionBranch_fb) 

                        if (elemYN(eDn,eYN_isCulvert)) then 
                            print *, 'CODE ERROR: culvert adjacent to junction not handled'
                            call util_crashpoint(2209873)
                        end if

                         !% --- ensure high Fr inflows have dQ/dH for junction
                        if ((elemR(eDn,er_FroudeNumber) .ge. +oneR) .or. &
                            (elemR(eUp,er_FroudeNumber) .ge. +oneR)       ) then  
                            !print *, 'High Froude into JB from upstream ',elemR(eUp,er_FroudeNumber)   
                            elemSR(eDn,esr_JunctionBranch_fa) = zeroR
                            elemSR(eDn,esr_JunctionBranch_fb) = zeroR 
                        end if

                    case (weir)
                            ! print *, 'is Weir'
                        !% --- upstream weir, downstream JB
                        elemSR(eDn,esr_JunctionBranch_fb) = zeroR
                        if (elemSR(eUp,esi_Weir_FlowDirection) == oneI) then 
                                !print *, 'is + flow'
                            !% --- flow into JB
                            if (elemYN(eUp,eYN_isSurcharged)) then
                                !% --- weir surcharged flow, 
                                !%     dQdHe > 0 because Q > 0
                                !%     increasing H at JB decreases He of weir and decreases magnitude of the negative Q (towards zero)
                                elemSR(eDn,esr_JunctionBranch_fa) = -elemSR(eUp,esr_Weir_dQdHe) 
                            else 
                                !% --- weir not surcharged
                                !%     increase H at JB has no effect on Q
                                elemSR(eDn,esr_JunctionBranch_fa) = zeroR
                            end if
                        else 
                                !print *, 'is - flow', eUp
                                !print *, elemSR(eUp,esr_Weir_dQdHe)
                            !% --- flow out of JB (reverse flow)
                            !%     dQdHe < 0 because Q < 0
                            !%     increasing H at JB increases He and increases magnitude of the negative Q (towards negative infinity) 
                            !%     applies to boths surcharged and non-surcharged
                            elemSR(eDn,esr_JunctionBranch_fa) = elemSR(eUp,esr_Weir_dQdHe)
                        end if

                    case (orifice)
                        !% --- upstream orifice, downstream JB
                        print *, 'CODE ERROR: orifice not handled'
                        call util_crashpoint(53873403)
                    case (outlet)
                        !% --- upstream outlet, downstream JB
                        print *, 'CODE ERROR: outlet not handled'
                        call util_crashpoint(53873403)
                    case (pump)
                        !% --- upstream pump, downstream JB
                        print *, 'CODE ERROR: pump not handled'
                        call util_crashpoint(53873404)
                    case default
                        print *, 'CODE ERROR: unexpected case default'
                        call util_crashpoint(538734)
                end select

                !print *, 'JB downstream ',elemSR(eDn,esr_JunctionBranch_fa), elemSR(eDn,esr_JunctionBranch_fb)
            
            end if

            !% --- handle downstream branch (upstream JB element)
            if (elemI(eUp,ei_elementType) == JB) then
                !% --- for different downstream element types
                select case (elemI(eDn,ei_elementType))
                    case (CC)
                        !% --- downstream conduit/channel, upstream JB
                        !% --- fa term
                        elemSR(eUp,esr_JunctionBranch_fa)                                     &
                            = + (crk * dt / fLengthAdj(fidx))                                &
                            * (                                                              &
                                + (fVelD(fidx)**twoI) * fTopAdj(fidx)                        &
                                + grav * (fAreaD(fidx)                                       &
                                        + fTopAdj(fidx) * (fHeadD(fidx) - fHeadAdj(fidx)))   &
                            ) 

                            ! print *, ' fa up ',  elemSR(eUp,esr_JunctionBranch_fa) 

                        !% --- fb term    
                        elemSR(eUp,esr_JunctionBranch_fb)                                     &
                            =  +crk * dt * grav * fTopAdj(fidx) / fLengthAdj(fidx) 

                            ! print *, ' fb up ',  elemSR(eUp,esr_JunctionBranch_fb) 

                        !% --- ensure high Fr inflows have dQ/dH for junction
                        if ((elemR(eDn,er_FroudeNumber) .le. -oneR) .or. &
                            (elemR(eUp,er_FroudeNumber) .le. -oneR)         ) then  
                            !print *, 'High Froude into JB from downstream ',elemR(eDn,er_FroudeNumber)   
                            elemSR(eUp,esr_JunctionBranch_fa) = zeroR
                            elemSR(eUp,esr_JunctionBranch_fb) = zeroR 
                        end if

                        if (elemYN(eUp,eYN_isCulvert)) then 
                            print *, 'CODE ERROR: culvert adjacent to junction not handled'
                            call util_crashpoint(3209873)
                        end if    

                    case (weir)
                        !print *, 'is Weir'
                        elemSR(eUp,esr_JunctionBranch_fa) = elemSR(eDn,esr_Weir_dQdHe)
                        elemSR(eUp,esr_JunctionBranch_fb) = zeroR

                        !% --- downstream weir, upstream JB
                        elemSR(eDn,esr_JunctionBranch_fb) = zeroR
                        if (elemSR(eDn,esi_Weir_FlowDirection) == oneI) then 
                            !% --- flow out of JB
                            !%     dQdHe > 0 because Q > 0
                            !%     increasing H at JB increases He of weir and increases magnitude of positive Q (towards positive inifinity)
                            !%     applies to surcharged and non-surcharged
                            elemSR(eUp,esr_JunctionBranch_fa) = elemSR(eDn,esr_Weir_dQdHe)
                        else 
                            !% --- flow into JB (reverse flow)
                            if (elemYN(eDn,eYN_isSurcharged)) then
                                !% --- surcharged weir inflow
                                !%     dQdHe < 0 because Q < 0
                                !%     increasing H at JB decreases He of weir and decreases magnitude of negative Q (towards zero)
                                elemSR(eUp,esr_JunctionBranch_fa) = -elemSR(eDn,esr_Weir_dQdHe)
                            else 
                                !% --- non-surcharged weir inflow
                                !%     increase H at JB has no effect on Q
                                elemSR(eUp,esr_JunctionBranch_fa) = zeroR
                            end if
                        end if

                    case (orifice)
                        !% --- downstream orifice, upstream JB
                        print *, 'CODE ERROR: outlet not handled'
                        call util_crashpoint(598734102)
                    case (outlet)
                        !% --- downstream outlet, upstream JB
                        print *, 'CODE ERROR: outlet not handled'
                        call util_crashpoint(598734103)
                    case (pump)
                        !% --- downstream pump, upstream JB
                        print *, 'CODE ERROR: pump not handled'
                        call util_crashpoint(598734104)
                    case default
                        print *, 'CODE ERROR: unexpected case default'
                        call util_crashpoint(5987341)
                end select

                !print *, 'JB upstream    ',elemSR(eUp,esr_JunctionBranch_fa), elemSR(eUp,esr_JunctionBranch_fb)

            end if
        end do

        !stop 2309874

    end subroutine junction_face_terms    
!%
!%==========================================================================
!%==========================================================================
!% 
    subroutine junction_force_Qinterpweights (thisCol)
        !%------------------------------------------------------------------
        !% Description:
        !% sets the JB interpweights to maximum so that adjacent
        !% values are interpolated to face
        !%------------------------------------------------------------------
            integer, intent(in) :: thisCol
            integer, pointer    :: npack, thisP(:)
            integer             :: ii
        !%------------------------------------------------------------------
        !% Aliases
            npack => npack_elemP(thisCol)
            if (npack < 1) return
            thisP => elemP(1:npack,thisCol)
        !%------------------------------------------------------------------

        do ii=1,max_branch_per_node  
            elemR(thisP+ii,er_InterpWeight_uQ) = setting%Limiter%InterpWeight%Maximum
            elemR(thisP+ii,er_InterpWeight_dQ) = setting%Limiter%InterpWeight%Maximum
        end do
        
    end subroutine junction_force_Qinterpweights
!%
!%==========================================================================
!%==========================================================================
!% 
    subroutine junction_mass_conservation (thisColP, istep)
        !%-----------------------------------------------------------------
        !% Description:
        !% Finds residual of junctions and adjusts JB and face values
        !% to obtain mass conservation
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisColP, istep
            integer, pointer    :: thisJM(:), Npack, JMidx
            real(8), pointer    :: Qoverflow(:), Qstorage(:), Qlateral(:)

            integer :: mm,ii
            real(8) :: resid !% residual in/out flowrate
            real(8) :: QnetIn, QnetOut, QnetBranches

            real(8), parameter :: local_epsilon = 1.0d-15
        !%-----------------------------------------------------------------
        !% Aliases
            Npack   => npack_elemP(thisColP)
            if (Npack < 1) return
            thisJM  => elemP(1:Npack,thisColP)
            Qoverflow => elemSR(:,esr_JunctionMain_OverflowRate)
            Qstorage  => elemSR(:,esr_JunctionMain_StorageRate)
            Qlateral  => elemR (:,er_FlowrateLateral)
        !%-----------------------------------------------------------------

            if (istep == 2) then
                print *, 'skipping junction mass conservation '
                return
            end if

        do mm=1,Npack
            JMidx => thisJM(mm)
            ! !% --- compute junction residual
            resid = junction_conservation_residual (thisJM(mm)) 
                    ! if (JMidx == printJM) 

            ! print *, '    resid at 4     ',resid, JMidx 

            if (abs(resid) > local_epsilon) then 
                !% --- note that QnetIn > 0 and QnetOut < 0
                QnetIn  = junction_branch_Qnet (thisJM(mm),+oneI)
                QnetOut = junction_branch_Qnet (thisJM(mm),-oneI) + Qoverflow(thisJM(mm))
                    ! if (JMidx == printJM) print *, 'starting conservation resid flows '
                    ! if (JMidx == printJM) print *, QnetIn, QnetOut, elemR(thisJM(mm),er_FlowrateLateral)

                !% --- ad hoc adjustment of elemR(:,er_Flowrate) of JB and Qoverflow of JM
                call junction_fix_conservation(thisJM(mm),resid, QnetIn, QnetOut)
                        ! if (JMidx == printJM) then 
                        !     print *, 'Q update 2'
                        !     do ii=1,max_branch_per_node
                        !         print *, ii, elemR(thisJM(mm)+ii,er_Flowrate)
                        !     end do
                        ! end if
               
                !% --- update net Q branches for residual changes
                QnetBranches = junction_main_sumBranches (thisJM(mm),er_Flowrate,elemR)
                        ! if (JMidx == printJM) print *, 'QnetBranches (end)', QnetBranches

                !% --- volume conservative storage rate (not based on dH)
                Qstorage(thisJM(mm)) = QnetBranches + Qoverflow(thisJM(mm)) + Qlateral(thisJM(mm))

                        ! if (JMidx == printJM) print *, 'QStorage (end)', Qstorage

                !% --- update Volume, VolumeOverflow and JB face values
                call junction_update_Qdependent_values (thisJM(mm), istep)

            end if 
        end do

    end subroutine junction_mass_conservation
!%    
!%==========================================================================
!% PRIVATE -- PRIMARY CALLS
!%==========================================================================
!%
    subroutine junction_calculation_3 (thisJM, Npack, istep)
        !%-----------------------------------------------------------------
        !% Description:
        !% Solves for the junction head and flows based on quadratic
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisJM(:), Npack, istep
            integer             :: JMidx

            real(8), pointer :: Qstorage(:), Qoverflow(:), Qlateral(:)

            real(8) :: Qnet, QnetBranches, QnetIn, QnetOut
            real(8) :: aQuad, bQuad, cQuad, rQuad
            real(8) :: dQdHoverflow, dQdHstorage, dH, resid
            real(8) :: dHhi, dHlo
            real(8), dimension(3) :: Hbound ! low, average, high
            real(8), dimension(2) :: deltaH
            integer, dimension(1) :: minDHloc
            integer :: mm, ii, fadj, bcount
            !real(8), parameter :: localEpsilon = 1.0d-6
        !%-----------------------------------------------------------------   
        !% Aliases
            Qstorage    => elemSR(:,esr_JunctionMain_StorageRate)
            Qoverflow   => elemSR(:,esr_JunctionMain_OverflowRate)
            Qlateral    => elemR (:,er_FlowrateLateral)
        !%----------------------------------------------------------------- 

        do mm=1,Npack
            JMidx = thisJM(mm)
            ! if (JMidx==printJM) then
            !     print *, ' '
            !     print *, '-----------------------------------------------'
            !     print *, 'mm, JM ',mm, JMidx
            ! end if

            !% --- get the max and min heads allowable for the JM
            call junction_main_head_bounds (JMidx, Hbound)

                ! if (JMidx==printJM) print *, ' '
                ! if (JMidx==printJM) print *,  'HBound'
                ! if (JMidx==printJM) print *, Hbound(1), Hbound(2), Hbound(3)

            !% --- set the allowable change in junction main head
            call junction_main_dHlimits (JMidx, dHlo, dHhi, Hbound)

                ! if (JMidx==printJM) print *, 'dHlo, dHhi ',dHlo, dHhi

            !% --- store face flowrate in JB for upstream (1) and downstream (2)
            call junction_branch_getface (elemR(:,er_Flowrate),fr_Flowrate,JMidx,ei_Mface_uL,1)
            call junction_branch_getface (elemR(:,er_Flowrate),fr_Flowrate,JMidx,ei_Mface_dL,2)

                ! if (JMidx==printJM) then 
                !     do ii=1,max_branch_per_node
                !     print *, 'Qbranch ',ii, elemR(JMidx+ii,er_Flowrate)
                !     end do
                ! end if

            !% --- compute net flowrate from branches (both CC and Diag)
            QnetBranches = junction_main_sumBranches (JMidx,er_Flowrate, elemR)

                ! if (JMidx==printJM) print *, '   QnetBranches ',QnetBranches
               
            !% --- compute overflow rate
            Qoverflow(JMidx) = junction_main_Qoverflow (JMidx)

                ! if (JMidx==printJM) print *, '   Qoverflow   ',Qoverflow(JMidx)

            !% --- net flowrate (Qnet > 0 is inflow)
            Qnet = QnetBranches + Qoverflow(JMidx) + Qlateral(JMidx)   

                ! if (JMidx==printJM) print *, '   Qnet        ',Qnet

            !% --- compute storage rate of change with head
            dQdHstorage = junction_main_dQdHstorage (JMidx,iStep)

                ! if (JMidx==printJM) print *, '   dQdHstorage ',dQdHstorage

            !% --- compute overflow rate with change in head
            dQdHoverflow = junction_main_dQdHoverflow (JMidx)

                !if (JMidx==printJM) print *, '   dQdHoverflow',dQdHoverflow

            !% --- compute the net fa from branches (both CC and Diag) which is part of the 'B' quadratic coefficient
            bQuad = junction_main_sumBranches(JMidx,esr_JunctionBranch_fa, elemSR)

                ! if (JMidx==printJM) print *, '   bQuad       ',bQuad

            !% --- adding the remainder of the B quadratic term
            bQuad = bQuad + dQdHoverflow - dQdHstorage

                ! if (JMidx==printJM) print *, '   bQuad (rev) ',bQuad

            !% --- compute the net fb from branches (both CC and Diag) which is the 'A' quadratic coefficient
            aQuad = junction_main_sumBranches(JMidx,esr_JunctionBranch_fb, elemSR)

                ! if (JMidx==printJM) print *, '   aQuad       ',aQuad

            !% --- compute the junction continuity source term whichis the 'C' of the quadratic equation
            cQuad = Qlateral(JMidx) + Qoverflow(JMidx) + QnetBranches

                ! if (JMidx==printJM) print *, '   cQuad       ',cQuad

            !% --- quadratic radical term
            rQuad = (bQuad**twoI) - fourR * aQuad * cQuad

                ! if (JMidx==printJM) print *, '   rQuad       ',rQuad

            !% --- compute junction head change
            dH = junction_head_change (JMidx,Qnet,aQuad,bQuad,cQuad,rQuad, Hbound, dHlo, dHhi)
   
                !  if (JMidx==printJM) print *, '   dH  0      ',dH

                ! if (JMidx==printJM) print *, ' '
                ! if (JMidx==printJM) print *, 'DH values'
                ! if (JMidx==printJM) print *, dHlo, dH, dHhi
                ! if (JMidx==printJM) print *, ' '


            !% --- limit junction head change by geometry and adjacent head
            dH = junction_head_limiter (JMidx,dH, dHhi, dHlo)

                ! if (JMidx==printJM) print *, '   dH  5      ',dH

            !% --- update JM head
            elemR(JMidx,er_Head) = elemR(JMidx,er_Head) + dH


            !% --- compute JB element dQdH from fa, fb and dH
            call junction_update_branch_dQdH (JMidx,dH)

            !% --- compute JB element DeltaQ using dQdH
            call junction_update_branch_DeltaQ (JMidx,dH)  

            !% --- update JB elements Q using Delta Q
            call junction_update_branch_flowrate (JMidx)

            !% --- update JB element H where JM head > branch Zbottom
            call junction_update_branch_head (JMidx,dH)


                ! if (JMidx==printJM) then 
                !     print *, 'dQdH pieces'
                !     do ii=1,2 !max_branch_per_node
                !         print *, elemSR(JMidx+ii,esr_JunctionBranch_fa), elemSR(JMidx+ii,esr_JunctionBranch_fb)
                !     end do
                ! end if

                ! if (JMidx==printJM) then 
                !     print *, 'dQdH update'
                !     do ii=1,2 !max_branch_per_node
                !         print *, elemSI(JMidx+ii,esi_JunctionBranch_Exists), ii, elemSR(JMidx+ii,esr_JunctionBranch_dQdH)
                !     end do
                ! end if


                ! if (JMidx==printJM) then
                !     print *, 'Q update 1'
                !     do ii=1,max_branch_per_node
                !         print *, ii, elemR(JMidx+ii,er_Flowrate)
                !     end do
                ! end if

            !% --- update junction main overflow rate
            Qoverflow(JMidx) = Qoverflow(JMidx) + dH * dQdHoverflow

            !% --- update net Q branches (included CC and Diag)
            QnetBranches = junction_main_sumBranches (JMidx,er_Flowrate,elemR)

                    ! if (JMidx==printJM) print *, 'QnetBranches (end)', QnetBranches

            !% --- update junction main storage flow rate
            Qstorage(JMidx) = junction_update_storage_rate  &
                                    (JMidx, dH, QnetBranches,istep) 
            
            !% ---          
                                
                ! if (JMidx==printJM) print *,  'depth compare ', Qstorage(JMidx)  &
                !     * setting%Time%Hydraulics%Dt * setting%Solver%crk2(istep) / elemSR(JMidx,esr_Storage_Plan_Area) &
                !     + elemR(JMidx,er_Depth) &
                !     , elemR(JMidx,er_Head) - elemR(JMidx,er_Zbottom)

                ! if (JMidx==printJM) print *, 'Qstorage new ',Qstorage(JMidx)

                ! if (JMidx==printJM) print *, 'Qstorage() - QnetBranches ', Qstorage(JMidx) - QnetBranches      
  
            !% --- update Volume, VolumeOverflow and JB face values
            call junction_update_Qdependent_values (JMidx, istep)

            !% TEST
            ! if (JMidx==printJM) then 
            !     if (elemR(11,er_Head) > elemR(10,er_Head)) then 
            !         print *, ' '
            !         print *, 'WARNING -- UPSTREAM HEAD INCREMENT'
            !         print *, 'heads ',elemR(10,er_Head), elemR(11,er_Head)
            !         print *, ' '
            !         !stop 76098723
            !     end if
            ! end if

        end do

        !stop 5509873

    end subroutine junction_calculation_3
!%
!%==========================================================================
!%==========================================================================
!%     
    subroutine junction_CC_for_JBadjacent (thisColP, istep, isUpstreamYN)
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
            
        !% --- pointer to the JB face
        if (isUpstreamYN) then 
            ! print *, ' '
            ! print *, 'UPSTREAM ================================='
            fIdx => elemI(:,ei_Mface_dL)
            bsign = +oneR
        else
            ! print *, ' '
            ! print *, 'DOWNSTREAM ============================'
            fIdx => elemI(:,ei_Mface_uL)
            bsign = -oneR
        end if    

        !% --- store present volume
        oldVolume(thisCC) = elemR(thisCC,er_Volume)

        !% --- changing the face flowrate changes the volume
        !%     on an upstream element, a negative DeltaQ causes
        !%     an increase in the upstream element volume
        elemR(thisCC,er_Volume) = elemR(thisCC,er_Volume) &
            - bsign * crk * dt *  faceR(fIdx(thisCC),fr_DeltaQ)

        ! if (isUpstreamYN) then 
        !     print *, fIdx(10)
        !     print *, bsign, crk, dt
        !     print *, 'face    DN for Up',faceR(fIdx(10),fr_DeltaQ)
        !     print *, 'Volume  UP ', oldVolume(10),elemR(10,er_Volume)
        !     print *, ' '
        ! else 
        !     print *, fIdx(22)
        !     print *, bsign, crk, dt
        !     print *, 'face    UP for Dn',faceR(fIdx(22),fr_DeltaQ)
        !     print *, 'Volume  Dn ', oldVolume(22),elemR(22,er_Volume)
        !     print *, 'Volume  Dn ', oldVolume(22),elemR(22,er_Volume) 
        !     print *, ' '
        ! end if
        ! print *, ' '


        !% --- changing the face flowrate changes the element velocity
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
        !%     index provided
        do mm=1,Npack
            call update_auxiliary_variables_CC (dummyIdx, dummyIdx, dummyIdx, &
                .false., .true., thisCC(mm))
        end do


    end subroutine junction_CC_for_JBadjacent
!%
!%==========================================================================
!% PRIVATE -- SECONDARY CALLS
!%==========================================================================
!% 
    subroutine junction_main_head_bounds (JMidx, Hbound)
        !%-----------------------------------------------------------------
        !% Description
        !% Computes the minimum head, Hbound(1), average head, Hbound(2) 
        !% and maximum heead, Hbound(3) for a junction
        !%     Hbound(3) is the maximum head in the surrounding elements
        !%     Hbound(2) is the average head of all elements
        !%     Hbound(1) is Zbottom of JM, or the lowest Z bottom of any 
        !%          branch if they are all higher than JM
        !%-----------------------------------------------------------------
            integer,               intent(in)    :: JMidx
            real(8), dimension(3), intent(inout) :: Hbound
            integer :: ii, bcount, fadj
        !%-----------------------------------------------------------------

        !%-----------------------------------------------------------------

            Hbound(1) = +huge(oneR)
            Hbound(2) = zeroR
            Hbound(3) = -huge(oneR)
            bcount    = zeroI

            !% --- cycle through branches (cannot be concurrent)
            do ii=1,max_branch_per_node
                if (elemSI(JMidx+ii,esi_JunctionBranch_Exists)==oneI) then 
                    bcount = bcount+oneI
                    if (mod(ii,2)== 0) then 
                        !% --- downstream branch
                        fadj  = elemI(JMidx+ii,ei_Mface_dL)
                    else
                        !% --- upstream branch
                        fadj  = elemI(JMidx+ii,ei_Mface_uL)
                    end if
                    Hbound(1) = min(Hbound(1),faceR(fadj,fr_Head_Adjacent))
                    Hbound(2) = Hbound(2)   + faceR(fadj,fr_Head_Adjacent)
                    Hbound(3) = max(Hbound(3),faceR(fadj,fr_Head_Adjacent))
                end if
            end do

            Hbound(1) = max(Hbound(1),elemR(JMidx,er_Zbottom))

            Hbound(2) = Hbound(2) / real(bcount,8)


    end subroutine junction_main_head_bounds    
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine junction_main_dHlimits (JMidx, dHlo, dHhi, Hbound)
        !%-----------------------------------------------------------------
        !% Descriptions
        !% sets the allowable change in head for the junction main
        !%-----------------------------------------------------------------
            integer,               intent(in)     :: JMidx
            real(8), dimension(3), intent(in)     :: Hbound
            real(8),               intent(inout)  :: dHlo, dHhi
        !%-----------------------------------------------------------------

        !% --- positive allowed change in dH
        dHhi = Hbound(3) - elemR(JMidx,er_Head)
        if (dHhi < zeroR) dHhi = zeroR

        !% --- negative allowed change in dH
        dHlo = Hbound(1) - elemR(JMidx,er_Head)
        if (dHlo > zeroR ) dHlo = zeroR

    end subroutine junction_main_dHlimits    
!%    
!%==========================================================================
!%==========================================================================
!% 
    pure subroutine junction_branch_getface (outdata, frCol, JMidx, fiIdx, kstart)
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
            outdata(k1:k2:2) = faceR(elemI(k1:k2:2,fiIdx),frCol)            &
                             * real(elemSI(k1:k2:2,esi_JunctionBranch_Exists),8)
        !end do

        !print *, 'in get face', outdata(k1:k2:2)

    end subroutine junction_branch_getface
!%
!%==========================================================================
!%==========================================================================
!% 
    real(8) function junction_main_sumBranches (JMidx,thisCol,thisArray)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the net values for a junction from all the branches
        !% Note -- assumes that all branches that do not exist have zero 
        !% values
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx, thisCol
            real(8), intent(in) :: thisArray(:,:)
            integer :: k1, k2
        !%------------------------------------------------------------------
        k1 = JMidx + 1
        k2 = JMidx + max_branch_per_node

        ! print *, ' '
        ! print *, 'junction_main_sumBranches ', JMidx
        ! print *, branchsign
        ! print *, ' '
        ! print *, thisArray(k1:k2,thisCol)
        ! print *, ' '

        !junction_main_sumBranches = sum(branchsign * thisArray(k1:k2,thisCol)) 
        junction_main_sumBranches = sum( branchsign                                       &
                                        * real(elemSI(k1:k2,esi_JunctionBranch_Exists),8) &
                                        * thisArray(k1:k2,thisCol) )

    end function junction_main_sumBranches   
!%    
!%========================================================================== 
!%==========================================================================
!% 
    real(8) function junction_main_Qoverflow (JMidx)  
        !%------------------------------------------------------------------
        !% Description
        !% Computes the overflow rate of a junction
        !% HACK: for an overflow orifice we have hard-coded the orifice
        !%   height and length. Later these should be user inputs
        !% HACK: for an overflow weir we have hard-coded the weir coefficient.
        !%   Later this should be moved to the settings structure
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx
        !%------------------------------------------------------------------  

        !% --- return zero if the head is below the crown at the start.
        if (elemR(JMidx,er_Head) .le. elemR(JMidx,er_Zcrown)) then 
            junction_main_Qoverflow = zeroR
            return
        end if

        select case (elemSI(JMidx,esi_JunctionMain_OverflowType))
            case (NoOverflow,Ponded)
                !% --- ponding an NoOverflow do not have separate volume accounting
                !%     for overflow
                junction_main_Qoverflow = zeroR
            case (OverflowWeir)
                !% --- weir overflow based on estimated circumference of storage plan area
                junction_main_Qoverflow = -coef1 * sqrt(elemSR(JMidx,esr_Storage_Plan_Area)) &
                    * ((elemR(JMidx,er_Head) - elemR(JMidx,er_Zcrown))**threehalfR)
            case (OverflowOrifice)
                !% --- orifice overflow assuming a single orifice of standard dimensions
                if (elemR(JMidx,er_Head) .le. (elemR(JMidx,er_Zcrown) + Horifice )) then
                    !% --- water surface below upper edge of orifice
                    junction_main_Qoverflow = -coef3 * Lorifice &
                        * (elemR(JMidx,er_Head) - elemR(JMidx,er_Zcrown))
                else
                    !% --- orifice is pressurized (head above the Zcrown + Horifice)
                    junction_main_Qoverflow = -coef3 * Lorifice  &
                        * (                                                                                 &
                            +  ((elemR(JMidx,er_Head) -  elemR(JMidx,er_Zcrown)            )**threehalfR)   &
                            -  ((elemR(JMidx,er_Head) - (elemR(JMidx,er_Zcrown) + Horifice))**threehalfR)   &
                          )
                end if
            case default
                !% --- should not reach here.
                !%     for debugging, change to impure function and uncomment
                !%     the following
                print *, 'CODE ERROR'
                print *, 'unexpected case default'
                call util_crashpoint(6209874)
        end select

    end function junction_main_Qoverflow   
!%    
!%==========================================================================
!%==========================================================================
!% 
    real(8) pure function junction_main_dQdHstorage (JMidx,iStep)  
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
            junction_main_dQdHstorage = zeroR
            return
        else
            !% --- Storage rate term based on time step
            junction_main_dQdHstorage = elemSR(JMidx,esr_Storage_Plan_Area) &
                / (setting%Solver%crk2(iStep) * setting%Time%Hydraulics%Dt)
        end if

    end function junction_main_dQdHstorage   
    !%    
!%========================================================================== 
!%==========================================================================
!% 
    real(8) pure function junction_main_dQdHoverflow (JMidx)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the overflow rate for  junction
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx
        !%------------------------------------------------------------------  

        ! print *, ' '
        ! print *, 'OverFlowType ',elemSI(JMidx,esi_JunctionMain_OverflowType)
        ! print *, NoOverflow, Ponded, OverflowWeir, OverflowOrifice
        ! print *, elemR(JMidx,er_Head) - elemR(JMidx,er_Zcrown)

        !% --- if not an overflow at this time
        if (elemR(JMidx,er_Head) .le. elemR(JMidx,er_Zcrown)) then 
            junction_main_dQdHoverflow = zeroR
            return 
        end if

        !% --- if possible overflow condition exists
        select case (elemSI(JMidx,esi_JunctionMain_OverflowType))  
            case (NoOverflow,Ponded)
                junction_main_dQdHoverflow = zeroR
            case (OverflowWeir)
                junction_main_dQdHoverflow = -coef2                            &
                    * sqrt(   ( elemSR(JMidx,esr_Storage_Plan_Area))           &
                            * ( elemR(JMidx,er_Head - elemR(JMidx,er_Zcrown))) &
                          )
            case (OverflowOrifice)
                junction_main_dQdHoverflow = -coef4 * Lorifice                  &
                        * sqrt(elemR(JMidx,er_Head) - elemR(JMidx,er_Zcrown))     
            case default
                !% --- should not reach here.
                !%     for debugging, change to impure function and uncomment
                !%     the following
                ! print *, 'CODE ERROR'
                ! print *, 'unexpected case default'
                ! call util_crashpoint(1209874)
        end select 

    end function junction_main_dQdHoverflow
    !%    
!%==========================================================================  
!%==========================================================================
!% 
    real(8) function junction_head_change &
        (JMidx, Qnet, aQuad, bQuad, cQuad, rQuad, Hbound, dHlo, dHhi)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes solution for head change using dQdH to get flux 
        !% conservative in/out flows
        !% the aQuad, bQuad and cQuad are the a,b,c coefficients of quadratic
        !% rQuad is the radical term
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx
            real(8), intent(in) :: Qnet, aQuad, bQuad, cQuad, Hbound(:), dHlo, dHhi
            real(8), intent(inout) :: rQuad
            
            integer :: kk
            real(8), dimension(2) :: deltaH, testVal
            real(8), dimension(max_branch_per_node) :: dQtrial1, dQtrial2
            real(8) :: loDelta, hiDelta

            real(8), parameter  :: localEpsilon = 1.0d-6
        !%------------------------------------------------------------------

        dQtrial1 = zeroR
        dQtrial2 = zeroR

        !% --- check for small 'a' coefficent of quadratic
        if (abs(aQuad) < localEpsilon) then 
            !% --- small 'a' coefficient
            if (abs(bQuad) < +localEpsilon) then 
                !% --- small 'b' coefficient
                !% --- with small 'a' this implies zero value
                print *, 'SMALL bquad'
                stop 239874
                junction_head_change = zeroR
                return
            else
                !% --- small 'a' implies linear solution
                print *, 'Small aquad', JMidx
                print *,  'alternate -c/b ',-cQuad / bQuad
                !stop 669874
                junction_head_change = -cQuad / bQuad
                return
            end if
        else 
            !% --- quadratic solution
            if (rQuad < zeroR) then 
                !% --- quadratic solution fails
                !%     use linear
                if (aQuad .ne. zeroR) then 
                    print *, 'small Rquad '
                    stop 6698734
                    junction_head_change = -cQuad / aQuad 
                    !% -- set the fb terms to zero
                    do kk=1,max_branch_per_node
                        elemSR(JMidx+kk,esr_JunctionBranch_fb) = zeroR 
                    end do
                    return
                else 
                    print *, 'unexpected else condition'
                    call util_crashpoint(77098723)
                end if
            else
                !% --- two roots of quadratic
                rQuad = sqrt(rQuad) 
                deltaH(1) = (-bQuad  + rQuad) / (twoR * aQuad)
                deltaH(2) = (-bQuad  - rQuad) / (twoR * aQuad)

                    ! if (JMidx==printJM) print *, '  JMhead      ',elemR(JMidx,er_Head)
                    ! if (JMidx==printJM) print *, '    dlo, dhi  ', dHlo, dHhi
                    ! if (JMidx==printJM) print *, '    deltaH(:) ', deltaH(1), deltaH(2)
                    ! if (JMidx==printJM) print *, '         Qnet ',Qnet

                loDelta = minval(deltaH)
                hiDelta = maxval(deltaH)    
         
                !% --- compute JB element dQdH from fa, fb and dH
                call junction_update_branch_dQdH (JMidx,deltaH(1))
                !% --- possible branch flowrate
                do kk=1,max_branch_per_node
                    if (elemSI(JMidx+kk,esi_JunctionBranch_Exists)==oneI) then 
                        dQtrial1(kk) = elemSR(JMidx+kk,esr_JunctionBranch_dQdH) * deltaH(1)
                    end if
                end do

                do kk=1,max_branch_per_node
                    if (elemSI(JMidx+kk,esi_JunctionBranch_Exists)==oneI) then 
                        dQtrial2(kk) = elemSR(JMidx+kk,esr_JunctionBranch_dQdH) * deltaH(2)
                    end if
                end do

                    ! if (JMidx==printJM) print *, 'dQ trial ',sum(dQtrial1), sum(dQtrial2)

                if (abs(sum(dQtrial1)) < abs(sum(dQtrial2))) then
                    junction_head_change = deltaH(1)
                    return
                else
                    junction_head_change = deltaH(2)
                    return 
                end if

                ! if ((loDelta < dHlo) .and. (hiDelta > dHhi)) then 

                !     print *, 'degenerate condition JMIDX = ',JMidx
                !     print *, '-cQuad / aQuad ', -cQuad / aQuad
                !     print *, '-cQuad / bQuad ', -cQuad / bQuad
                !     stop 509874
                !     junction_head_change = zeroR
                !     !return
                ! else
                !     if (loDelta .ge. dHlo) then 
                !         junction_head_change = loDelta 
                !         return
                !     else
                !         junction_head_change = hiDelta
                !         return
                !     end if 
                ! end if


                ! if (deltaH(1)*deltaH(2) < zeroR) then
                !     !% --- one negative root and one positive
                !     !%     so there is only one solution that is consistent
                !     if (Qnet > zeroR) then
                !         !% --- for positive Qnet, dH should be positive
                !         !%     choose the positive root
                !         junction_head_change = maxval(deltaH)
                !         return
                !     elseif (Qnet < zeroR) then
                !         !% --- for negative Qnet, dH should be negative
                !         !%     choose the negative root
                !         junction_head_change = minval(deltaH)
                !         return
                !     else 
                !         !% --- for no Qnet, dH should be zero
                !         junction_head_change = zeroR
                !         return
                !     end if
                ! elseif (deltaH(1)* deltaH(2) > zeroR) then 
                !     !% --- two positive or two negative roots
                !     if (deltaH(1) > zeroR) then 
                !         !% --- two positive roots, implies Qnet > 0

                !         !% using smaller value
                !         junction_head_change = minval(deltaH)

                !         ! if (elemR(JMidx,er_Head) > Hbound(2)) then 
                !         !     !% --- depth is above average, so use the smaller value
                !         !     junction_head_change = minval(deltaH)
                !         ! else 
                !         !     !% --- depth is below the average, so use the larger value 
                !         !     junction_head_change = maxval(deltaH)
                !         ! end if
                !     else
                !         !% --- two negative roots, imlies Qnet < 0

                !         !% --- using larger (negative) magnitude 
                !         junction_head_change = minval(deltaH)
                        
                !         ! if (elemR(JMidx,er_Head) > Hbound(2)) then 
                !         !     !% --- depth is above average, so use the larger (negative value)
                !         !     junction_head_change = minval(deltaH)
                !         !     return
                !         ! else
                !         !     !% --- depth is below average, so use the smaller (negative value)
                !         !     junction_head_change = maxval(deltaH)
                !         !     return
                !         ! end if
                !     end if
                ! else 
                !     junction_head_change = zeroR
                ! end if


                ! if (deltaH(1)*deltaH(2) < zeroR) then
                !     !% --- one negative root and one positive
                !     !%     so there is only one solution that is consistent
                !     if (Qnet > zeroR) then
                !         !% --- for positive Qnet, dH should be positive
                !         !%     choose the positive root
                !         junction_head_change = maxval(deltaH)
                !         return
                !     elseif (Qnet < zeroR) then
                !         !% --- for negative Qnet, dH should be negative
                !         !%     choose the negative root
                !         junction_head_change = minval(deltaH)
                !         return
                !     else 
                !         !% --- for no Qnet, dH should be zero
                !         junction_head_change = zeroR
                !         return
                !     end if
                ! elseif (deltaH(1)*deltaH(2) > zeroR) then
                    ! !% --- either two negative roots or two positive roots
                    ! !%     choose dH that causes smallest net change in flowrates
                    ! do kk=1,max_branch_per_node 
                    !     dQtrial1(kk) = junction_branchDQ_from_dH (JMidx+kk,deltaH(1))
                    !     dQtrial2(kk) = junction_branchDQ_from_dH (JMidx+kk,deltaH(2))
                    ! end do
                    ! !% --- get the magnitude of flowrate changes required for two dH values
                    ! testVal(1) = sum(abs(dQtrial1))
                    ! testVal(2) = sum(abs(dQtrial2))

                    ! if (JMidx==printJM) print *, '   testVal(:) ', testVal(1), testVal(2)
                    ! !% -- select the dH causing smaller flowrate changes
                    ! if (testVal(1) < testVal(2)) then
                    !     junction_head_change = deltaH(1)
                    ! elseif (testVal(1) > testVal(2)) then
                    !     junction_head_change = deltaH(2)
                    ! else
                    !     !% --- if same flowrate change then choose the ??dH
                    !     if (abs(deltaH(1)) .le. abs(deltaH(2))) then 
                    !     !if (abs(deltaH(2)) .le. abs(deltaH(1))) then 
                    !         junction_head_change = deltaH(1)
                    !     else
                    !         junction_head_change = deltaH(2)
                    !     end if
                    ! end if

                !     if (deltaH(1) > zeroR) then 
                !         !% --- two positive roots choose the largest consistent with Hmax
                !         if (maxval(deltaH) .le. dHhi) then 
                !             junction_head_change = maxval(deltaH)
                !             return
                !         else 
                !             junction_head_change = minval(deltaH)
                !             return
                !         end if
                !     else
                !         !% --- two negative roots, choose largest negative magnitude 
                !         !%     (smallest negative number) consistent with with Hmin 
                !         if (minval(deltaH) .ge. dHlo) then 
                !             junction_head_change = minval(deltaH)
                !             return
                !         else
                !             junction_head_change = maxval(deltaH)
                !             return
                !         end if
                !     end if
                ! else 
                !     junction_head_change = zeroR
                !     return
                ! end if

            end if
        end if

    end function junction_head_change
!%    
!%==========================================================================  
!%==========================================================================
!% 
    real(8) function junction_head_limiter (JMidx, dHin, dHhi, dHlo) 
        !%------------------------------------------------------------------ 
        !% Description:
        !% limits the allowable head change based on geometry and inflow
        !% head available  
        !%------------------------------------------------------------------
            integer, intent(in) :: JMidx
            real(8), intent(in) :: dHin, dHhi, dHlo

            real(8) :: dH
        !%------------------------------------------------------------------
        dH = dHin    

        if (dH > dHhi) dH = dHhi
        if (dH < dHlo) dH = dHlo

        ! !% --- limit dH to prevent flow direction change
        ! if (dH > zeroR) then
        !     dH = min(dH, junction_dH_maxgain (JMidx,dH))
        ! elseif (dH < zeroR) then
        !     dH = max(dH, junction_dH_maxloss (JMidx,dH))
        ! else 
        !     !% -- if dH = zero, no change
        ! end if
        !     ! if (JMidx==31) print *, '   dH  1      ',dH

        ! !% --- limit dH to prevent negative head
        ! if ((elemR(JMidx,er_Head) - elemR(JMidx,er_Zbottom) + dH) < setting%ZeroValue%Depth) then 
        !     ! ! print *, 'limiters '
        !     ! ! print *, elemR(JMidx,,er_Head),  (elemR(JMidx,,er_Zbottom) + 0.99d0*setting%ZeroValue%Depth)
        !     dH = min(dH,(elemR(JMidx,er_Zbottom) + 0.99d0*setting%ZeroValue%Depth) - elemR(JMidx,er_Head))
        ! end if
        !     ! if (JMidx==31) print *, '   dH  2      ',dH

        ! !% --- limit dH dropping based on where overflow shuts off
        ! if (elemSR(JMidx,esr_JunctionMain_OverflowRate) .ne. zeroR) then 
        !     dH = max(dH, junction_dH_overflow_min(JMidx))
        ! end if
            ! if (JMidx==31) print *, '   dH  3      ',dH

        junction_head_limiter = dH

    end function junction_head_limiter
!%    
!%==========================================================================  
!%==========================================================================
!% 
    real(8) pure function junction_dH_overflow_min(JMidx)
        !%------------------------------------------------------------------
        !% Description:
        !% Limit the change of dH during overflow so that at most it goes
        !% an epsilon distance below the crown
        !%------------------------------------------------------------------
            integer, intent(in) :: JMidx

            real(8), parameter :: localEpsilon = 1.0d-4
        !%------------------------------------------------------------------

        junction_dH_overflow_min = &
             -(elemR(JMidx,er_Head) - elemR(JMidx,er_Zcrown)) - localEpsilon

    end function junction_dH_overflow_min
!%    
!%==========================================================================
!%==========================================================================
!% 
    real(8) pure function junction_branchDQ_from_dH (JBidx,dH)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes Q change for junction branch JBidx = JMidx + kk
        !% Computes dQdH = fa + fb*dH, but does not store.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JBidx
            real(8), intent(in) :: dH 
        !%------------------------------------------------------------------ 
        !% Preliminaries 
            if (elemSI(JBidx,esi_JunctionBranch_Exists .ne. oneI)) return
        !%------------------------------------------------------------------ 

        junction_branchDQ_from_dH =                             &
              (   elemSR(JBidx,esr_JunctionBranch_fa)           &
                + elemSR(JBidx,esr_JunctionBranch_fb) * dH ) * dH

    end function junction_branchDQ_from_dH
!%    
!%==========================================================================
!%==========================================================================
!% 
    subroutine junction_update_branch_dQdH (JMidx, dH)
        !%------------------------------------------------------------------
        !% Description:
        !% Updates the dQdH for a JB based on the fa, fb and dH values
        !%------------------------------------------------------------------
            integer, intent(in) :: JMidx
            real(8), intent(in) :: dH
            integer :: ii
        !%------------------------------------------------------------------

        do ii=1,max_branch_per_node
            elemSR(JMidx+ii,esr_JunctionBranch_dQdH)                  &
                 = real(elemSI(JMidx+ii,esi_JunctionBranch_Exists),8) &
                   * (  elemSR(JMidx+ii,esr_JunctionBranch_fa)        &
                      + elemSR(JMidx+ii,esr_JunctionBranch_fb) * dH)      
        end do
    
    end subroutine junction_update_branch_dQdH
!%    
!%==========================================================================
!%==========================================================================
!%     
    subroutine junction_update_branch_DeltaQ (JMidx, dH)   
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
        end do

    end subroutine junction_update_branch_DeltaQ
!%    
!%==========================================================================
!%==========================================================================
!% 
    subroutine junction_update_branch_flowrate (JMidx)   
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

    end subroutine junction_update_branch_flowrate
!%    
!%==========================================================================
 !%==========================================================================
!% 
    subroutine junction_update_branch_head (JMidx, dH)
        !%------------------------------------------------------------------
        !% Description
        !% Applies the updated Junction Main head to branches that are 
        !% influenced by Junction Main
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx
            real(8), intent(in) :: dH

            real(8), pointer :: head(:), velocity(:), Zbottom(:)
            real(8), pointer :: SedDepth(:), jbKfactor(:), grav

            integer, pointer :: jbExists(:)

            integer :: ii, jb
        !%------------------------------------------------------------------
        !% Aliases
            jbExists => elemSI(:,esi_JunctionBranch_Exists)

            head      => elemR(:,er_Head)
            velocity  => elemR(:,er_Velocity)
            Zbottom   => elemR(:,er_Zbottom)
            SedDepth  => elemR(:,er_SedimentDepth)
            jbKfactor => elemSR(:,esr_JunctionBranch_Kfactor)
            grav      => setting%Constant%gravity
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        
        !% -- cycle through branches
        do ii=1, max_branch_per_node
            jb = JMidx+ii
            !% --- skip branches that don't exist
            if (jbExists(jb) .ne. oneI) cycle 
            !% --- check if head of main is above bottom of branch
            if (head(JMidx) > (Zbottom(jb) + SedDepth(jb)) ) then 
                !% --- set junction head to main head
                head(jb) = head(JMidx)
                !% --- adjust for Kfactor
                if (jbKfactor(jb) > zeroR) then 
                    head(jb) = head(jb)                                        &
                        + branchsign(ii) * sign(oneR,Velocity(jb))             &
                        * jbKfactor(jb) * (Velocity(jb)**twoR)  / (twoR * grav) 
                else 
                    !% ---K=0 has no effect
                end if      
            else 
                !% --- JB head is unchanged by junction solution 
                !%     because mJM head is too low
            end if         
        end do

    end subroutine junction_update_branch_head    
!%    
!%==========================================================================  
!%==========================================================================
!% 
    real(8) function junction_update_storage_rate (JMidx, dH, QnetBranches, istep)
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
                junction_update_storage_rate = zeroR
                return
            case (TabularStorage,FunctionalStorage,ImpliedStorage)
                !% --- compute the storage flowrate
                if (istep == 1) then
                    !% --- compute rate from dH
                    junction_update_storage_rate                                        &
                        = (elemSR(JMidx,esr_Storage_Plan_Area) * dH)                    &
                            / (setting%Solver%crk2(istep) * setting%Time%Hydraulics%Dt)
                    return
                elseif (istep == 2) then 
                    !% --- compute rate from mass conservation
                    junction_update_storage_rate = QnetBranches       &
                        + elemSR(JMidx,esr_JunctionMain_OverflowRate) &
                        +  elemR(JMIdx,er_FlowrateLateral)
                else 
                    !% should not be possible
                    print *, 'CODE ERROR: unexpected else'
                    call util_crashpoint(6209874)
                    return
                end if

            case default
                junction_update_storage_rate = zeroR
                print *, 'CODE ERROR: unexpected case default'
                call util_crashpoint(8852783)
                return

        end select

    end function junction_update_storage_rate   
!%    
!%========================================================================== 
!%==========================================================================
!% 
    subroutine junction_update_Qdependent_values (JMidx,istep) 
        !%------------------------------------------------------------------
        !% Description:
        !% Updates Qrate values after Q branches etc. have changed
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx, istep
            real(8), pointer    :: Qstorage, Qoverflow

            real(8) :: QnetBranches
        !%------------------------------------------------------------------
        !% Alias
            Qstorage  => elemSR(JMidx,esr_JunctionMain_StorageRate)
            Qoverflow => elemSR(JMidx,esr_JunctionMain_OverflowRate)
        !%------------------------------------------------------------------

        !% --- update volume overflow 
        elemR(JMidx,er_VolumeOverflow) = Qoverflow &
                * setting%Time%Hydraulics%Dt * setting%Solver%crk2(istep) 

        !% --- update storage volume
        if (elemSI(JMidx,esi_JunctionMain_Type) .ne. NoStorage) then
            elemR(JMidx,er_Volume) = elemR(JMidx,er_Volume_N0)                         &
                + setting%Time%Hydraulics%Dt * setting%Solver%crk2(istep) * Qstorage
        else 
            !% --- no storage does not change volume
            elemR(JMidx,er_Volume) = elemR(JMidx,er_Volume_N0)
        end if

        

    end subroutine junction_update_Qdependent_values
!%    
!%==========================================================================
!%==========================================================================
! !%
!     subroutine junction_update_face_values (JMidx)
!         !%------------------------------------------------------------------
!         !% Descriptions
!         !% Copies JB element values to faces
!         !%  upstream called with "oneI" argument, downstream called with "twoI"
!         !%------------------------------------------------------------------
!             integer, intent(in) :: JMidx
!         !%------------------------------------------------------------------

!     end subroutine junction_update_face_values 
!%    
!%==========================================================================
!%==========================================================================
!% 
    subroutine junction_branchface_forceJBvalue (frCol, erCol, fiIdx, JMidx, kstart)
        !%------------------------------------------------------------------
        !% Description:
        !% Forces the JB element value on the adjacent face
        !%------------------------------------------------------------------
            integer, intent(in) :: frCol  !% column in faceR array for output
            integer, intent(in) :: erCol  !% column in elemR array for input
            integer, intent(in) :: fiIdx  !% face index column for up/dn map
            integer, intent(in) :: JMidx  !% junction main index
            integer, intent(in) :: kstart !% =1 for upstream, 2 for downstream
            integer :: k1, k2
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        k1 = JMidx + kstart
        k2 = JMidx + max_branch_per_node

        faceR(elemI(k1:k2:2,fiIdx),frCol) = elemR(k1:k2:2,erCol)
    
    end subroutine junction_branchface_forceJBvalue
!%    
!%==========================================================================
!% CONSERVATION
!%==========================================================================
!%
    real(8)  function junction_conservation_residual (JMidx) 
        !%------------------------------------------------------------------
        !% Description:
        !% computes flowrate residual, which is >0 for too much inflow 
        !% and < 0 for too much outflow
        !%------------------------------------------------------------------
            integer, intent(in) :: JMidx
            real(8) :: QnetBranches
        !%------------------------------------------------------------------

        QnetBranches = junction_main_sumBranches (JMidx,er_Flowrate, elemR)

        ! if (JMidx ==printJM) then
            ! print *, ' '
            ! print *, 'JUNCTION CONSERVATION RESID ', JMidx
            ! print *, 'QnetBranches: ',QnetBranches
            ! print *, 'lateral, storage, overflow'
            ! print *, elemR(JMidx,er_FlowrateLateral), elemSR(JMidx,esr_JunctionMain_StorageRate), elemSR(JMidx,esr_JunctionMain_OverflowRate)
            ! print *, ' '
        !end if
        
        junction_conservation_residual = QnetBranches         &
             + elemSR(JMidx,esr_JunctionMain_OverflowRate)    &
             - elemSR(JMidx,esr_JunctionMain_StorageRate)     &
             + elemR (JMidx,er_FlowrateLateral)

    end function junction_conservation_residual
!%    
!%==========================================================================
!%==========================================================================
!% 
    real(8) pure function junction_branch_Qnet (JMidx,idir)    
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

        junction_branch_Qnet = sum(Qdir)
        
    end function junction_branch_Qnet
!%    
!%==========================================================================
!%==========================================================================
!%    
    subroutine junction_fix_conservation (JMidx, resid, QnetIn, QnetOut)
        !%------------------------------------------------------------------
        !% Description:
        !% Ad hoc fix to junction flow rates to ensure mass conservation
        !% Typically required because dH adjustment was limited
        !% Modifies elemR(:,er_Flowrate) and Qoverflow
        !%------------------------------------------------------------------
            integer, intent(in) :: JMidx
            real(8), intent(inout) :: resid,  QnetOut
            real(8), intent(in) :: QnetIn

            real(8), pointer :: Qoverflow, Qstorage, dQdH(:)

            integer, pointer :: fup(:), fdn(:)

            integer :: kk
            real(8) :: QratioIn, QratioOut, dQoverflow, Hinc

            real(8), dimension(max_branch_per_node) :: dQ, dH, areaQ
            real(8), parameter :: localEpsilon = 1.0d-6
            real(8) :: zbottom, Aout, Ain
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
        !%------------------------------------------------------------------

        repeatYN = .true.

        ! print *, ' '
        ! print *, 'JUNCTION ',JMidx
        ! print *, 'starting resid ',resid 

        dQ = zeroR
        dH = zeroR
        dQoverflow = zeroR

        !% --- note that by definition: QnetIn > 0 and QnetOut < 0 
        !%     resid > 0 implies too much inflow
        !%     resid < 0 implies too much outflow

        !% --- reduce or eliminate Qoverflow rate for negative residual
        if ((resid < zeroR) .and. (Qoverflow < zeroR)) then 
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
        endif
            
        do while (repeatYN)
        ! select case (elemSI(JMidx,esi_JunctionMain_Type))
        !     case (TabularStorage,FunctionalStorage,ImpliedStorage)
        !         !% Conservation should NOT be an issue
        !         print *, 'RESIDUAL ',resid
        !         repeatYN = .false.
        !     case (NoStorage)
                ! print *, 'IN NO STORAGE'
                ! print *, 'QnetIn, QnetOut ',QnetIn, QnetOut
                ! !% --- check for a degenerate condition with no branch flows
                ! if ((QnetIn .le. zeroR) .and. (QnetOut .ge. zeroR)) then
                !     !% --- positive residual implies flow that must exit
                !     !%     Zero (or slightly negative) QnetIn implies no inflows
                !     !%     Zero (or slightly positive) QnetOut implies no outflows
                !     !%     typically occurs when lateral inflow exists but
                !     !%     there are no other inflows or outflows    
                !     if (resid > zeroR) then
                !         !% --- find the lowest channel(s) connected to junction
                !         !%     neglect all diagnostic adjacent
                !         zbottom = huge(oneR)
                !         zcount  = zeroI
                !         do kk=1,max_branch_per_node
                !             if ((elemSI(JMidx+kk,esi_JunctionBranch_Exists)        == oneI)  .and. &
                !                 (elemSI(JMidx+kk,esi_JunctionBranch_Diag_adjacent) == zeroI) .and. &
                !                 (elemSI(JMidx+kk,esi_JunctionBranch_CanModifyQ)    == oneI)          ) then
                !                 zcount = zcount + oneI
                !                 zbottom = min(zbottom,elemR(JMidx+kk,er_Zbottom))
                !             end if
                !         end do
                !         if (zcount .ge. oneI) then 
                !             !% --- spread the residual as outflows over lowest channels
                !             do kk=1,max_branch_per_node
                !                 if ((elemSI(JMidx+kk,esi_JunctionBranch_Exists)        == oneI)  .and. &
                !                     (elemSI(JMidx+kk,esi_JunctionBranch_Diag_adjacent) == zeroI) .and. &
                !                     (elemSI(JMidx+kk,esi_JunctionBranch_CanModifyQ)    == oneI)  .and. &
                !                     (zbottom == elemR(JMidx+kk,er_Zbottom))                              ) then
                !                     !% --- new outflow
                !                     dQ(kk) = - branchsign(kk) * resid / real(zcount,8)
                !                 end if
                !             end do
                !         else 
                !             !% --- cannot find a non-diagnostic element to fix residual, so it is lost
                !             print *, 'WARNING 928732: junction residual inflow that is ignored'
                !             print *, 'as there are no non-diagnostic branches'
                !             print *, 'Junction index in elemR array is',JMidx
                !         end if
                !     else !% (resid <= 0)
                !         !% --- residual < 0 implies too much outflow so a reduced outflow is needed (which 
                !         !%     is not possible) or an increased inflow is needed (also not possible)
                !         !%     since there are no viable branches.
                !         !%     Typically occurs when there is a negative lateral inflow and there is
                !         !%     inflow or outflow to the junction.
                !             print *, 'WARNING 928734: junction residual outflow that is ignored'
                !             print *, 'as there are no non-diagnostic branches'
                !             print *, 'Junction index in elemR array is',JMidx
                !     end if

                ! else !% --- for cases with some inflows and/or outflows

                    !% --- distribute change in in/outflows based on contributing area
                    bFixYN = .false.
                    Aout = zeroR
                    Ain  = zeroR
                    areaQ = zeroR
                    bcount = zeroI
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
                            print *, 'CODE ERROR: bcount = 0; should not have reached this point '
                            call util_crashpoint(6209873)
                        end if

                            ! print *, 'QVAL ',(real(branchsign(kk),8) * elemR(JMidx+kk,er_Flowrate))


                    end do

                    if (bcount < 1) then 
                        !% --- occurs during wetting drying, use all real branches for adjustment
                        do kk=1,max_branch_per_node
                            if ((elemSI(JMidx+kk,esi_JunctionBranch_Exists)        .ne. oneI) .or.  &
                                (elemSI(JMidx+kk,esi_JunctionBranch_CanModifyQ)    .ne. oneI)        ) cycle

                            if (elemSI(JMidx+kk,esi_JunctionBranch_IsUpstream) == oneI) then
                                !% --- upstream branch
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
                            print *, 'CODE ERROR: Mass residual in an unexpected condition'
                            call util_crashpoint(2298522)
                        endif
                    end if

                    do kk=1,max_branch_per_node
                        !% --- accumulate the in and outflow areas
                        if ((real(branchsign(kk),8) * elemR(JMidx+kk,er_Flowrate)) > zeroR) then 
                            Ain  = Ain  + areaQ(kk)
                                ! print *, 'Ain ',Ain
                        else
                            Aout = Aout + areaQ(kk)
                                ! print *, 'Aout ',Aout
                        end if
                       
                        if (Qoverflow < zeroR) then 
                            Aout = Aout + elemSR(JMidx,esr_JunctionMain_OverflowOrifice_Length) &
                                        * elemSR(JMidx,esr_JunctionMain_OverflowOrifice_Height)
                        end if
                    end do

                !     print *, ' '
                !    ! print *, 'bfix ',bFixYN
                !     print *, 'Ain, Aout ',Ain, Aout
                !     print *, ' '
                    

                    if ((Ain < setting%ZeroValue%Area) .and. (Aout < setting%ZeroValue%Area)) then 
                        !% degenerate condition 
                        print *, 'CODE ERROR: unexpected junction are condition '
                        print *, 'JMidx ',JMidx
                        print *, Ain, Aout, setting%ZeroValue%Area
                        call util_crashpoint(629784)
                    else 
                        !% --- distribute residual in proportion to area
                        do kk=1,max_branch_per_node
                            if (.not. bFixYN(kk)) cycle 
                            dQ(kk) = - resid * real(branchsign(kk),8) * areaQ(kk) / (Ain + Aout) 
                        end do
                        !% --- overflow rate adjustment
                        if (Qoverflow < zeroR) then 
                            dQoverflow = -resid * elemSR(JMidx,esr_JunctionMain_OverflowOrifice_Length) &
                                                * elemSR(JMidx,esr_JunctionMain_OverflowOrifice_Height) &
                                                / (Ain + Aout)

                            !% --- check to see if too much overflow was removed
                            if (Qoverflow + dQoverflow > zeroR) then 
                                resid = resid - Qoverflow
                                Qoverflow = zeroR
                                repeatYN = .true.
                            else 
                                repeatYN = .false.
                            endif
                        else
                            dQoverflow = zeroR
                            repeatYN = .false.
                        end if
                    end if

                    ! print *, 'dQ : ', dQ

                    ! if ((QnetOut .ge. zeroR) .and. (QnetIn > zeroR)) then 
                    !     !% --- only inflows to modify (no outflow)
                    !     !%     includes diagnostic except for those frozen
                    !     do kk=1,max_branch_per_node
                    !         !% --- select inflows to adjust
                    !         !%     signs based on QnetIn > 0
                    !         if (((branchsign(kk) * elemR(JMidx+kk,er_Flowrate)) > zeroR) .and. &
                    !             (elemSI(JMidx+kk,esi_JunctionBranch_Exists)     == oneI) .and. &
                    !             (elemSI(JMidx+kk,esi_JunctionBranch_CanModifyQ) == oneI)         ) then
                    !             !% --- adjust by proportion of total inflow
                    !             dQ(kk) = - resid * elemR(JMidx+kk,er_Flowrate) / QnetIn
                    !         end if
                    !     end do
                    !     !% --- no adjustment to overflow
                    !     dQoverflow = zeroR
            
                    ! elseif ((QnetOut < zeroR) .and. (QnetIn .le. zeroR)) then 
                    !     !% --- only outflows to modify (no inflow)
                    !     !%     note that signs based on QnetOut < 0
                    !     !%     includes diagnostic except for those frozen
                    !     do kk=1,max_branch_per_node
                    !         !% -- select outflows to adjust
                    !         if (((branchsign(kk) * elemR(JMidx+kk,er_Flowrate)) < zeroR) .and. & 
                    !             (elemSI(JMidx+kk,esi_JunctionBranch_Exists)     == oneI) .and. &            
                    !             (elemSI(JMidx+kk,esi_JunctionBranch_CanModifyQ) == oneI)         ) then
                    !             !% --- adjust by proportion of tototal outflow
                    !             dQ(kk) = - resid * elemR(JMidx+kk,er_Flowrate) / QnetOut
                    !         end if
                    !     end do
                    !     !% --- adjust overflow for proportional residual
                    !     dQoverflow = - resid * Qoverflow / QnetOut
            
                    ! else 
                    !     !% --- adjusting both inflows and outflows
                    !     QratioIn =  QnetIn  / (QnetIn - QnetOut)
                    !     QratioOut= -QnetOut / (QnetIn - QnetOut)
            
                    !     if (JMidx == printJM) print *, ' '
                    !     if (JMidx == printJM) print *, 'fixing conservation both inflow/outflow'
                    !     if (JMidx == printJM) print *, 'Qratio ',QratioIn, QratioOut
                    !     if (JMidx == printJM) print *, 'resid  ',resid
            
                    !     Qcount = zeroI
                    !     do kk=1,max_branch_per_node

                    !            if (JMidx == printJM) print *, 'kk, dQdH ',kk, dQdH(JMidx+kk)

                    !         if ((elemSI(JMidx+kk,esi_JunctionBranch_Exists)     == oneI) .and. &
                    !             (elemSI(JMidx+kk,esi_JunctionBranch_CanModifyQ) == oneI)         ) then 

                    !             if     ((branchsign(kk) * elemR(JMidx+kk,er_Flowrate)) > zeroR) then
                    !                 !% --- inflows affected in proportion
                    !                 dQ(kk) = - resid * QratioIn * elemR(JMidx+kk,er_Flowrate) / QnetIn
                    !             elseif ((branchsign(kk) * elemR(JMidx+kk,er_Flowrate)) < zeroR)  then
                    !                 !% --- outflows affected in proportion
                    !                 dQ(kk) = - resid * QratioOut * elemR(JMidx+kk,er_Flowrate) / QnetOut
                    !             elseif (branchsign(kk) * elemR(JMidx+kk,er_Flowrate) == zeroR) then
                    !                 dQ(kk) = zeroR
                    !             else 
                    !                 print *, 'CODE ERROR: unexpected else'
                    !                 call util_crashpoint(739874)
                    !                 return 
                    !             end if
                    !             !% --- associated change in head
                    !             if ((abs(dQdH(JMidx+kk)) > localEpsilon) .and. &
                    !                 (.not. elemYN(JMidx+kk,eYN_isZeroDepth))) then
                    !                 dH(kk) = dQ(kk) / dQdH(JMidx+kk)
                    !                 Qcount = Qcount + 1
                    !             else
                    !                 dH(kk) = zeroR
                    !             end if

                    !             ! if (JMidx == 31) print *, 'dQdH', elemSR(JMidx+kk,esr_JunctionBranch_dQdH)
                    !             ! if (JMidx == 31) print *, 'kk, dQ, dH ', dQ(kk), dH(kk)

                    !         else
                    !             !% no branch, no action
                    !         end if
                    !     end do
            
                        !% --- overflow rate adjustment
                        !dQoverflow = - resid * QratioOut * Qoverflow / QnetOut

                        ! if (JMidx == 31) then 
                        !     if (Qcount > 0) then
                        !     print *, 'Qcount ',Qcount
                        !     print *, 'sum(dH), sum(dH)/Qcount ' , sum(dH), sum(dH)/real(Qcount,8)
                        !     end if
                        ! end if

                        !% --- Junction head change implied by flowrate change
                        !%     head change is average of dH(kk) for branches with valid dQ/dH
                        ! if (Qcount > 0) then
                        !     Hinc = sum(dH) / real(Qcount,8)

                        !     ! if (JMidx == 31) then
                        !     !     print *, ' '
                        !     !     print *, 'JMidx, conservation DH ',JMidx, Hinc
                        !     ! end if
                            
                        !     ! do kk=1,max_branch_per_node
                        !     !    if (elemR(JMidx,er_Head)    > elemR(JMidx+kk,er_Zbottom)) then 
                        !     !         !% --- adjust branch head where junction head is above branch bottom
                        !     !        elemR(JMidx+kk,er_Head) = elemR(JMidx+kk,er_Head) + Hinc
                        !     !    end if
                        !     ! end do
                        !     ! elemR(JMidx,er_Head) = elemR(JMidx,er_Head) + Hinc
                        ! end if
                    

                ! end if

        !     case default
        !         print *, 'CODE ERROR: need to setup storage residual adjustment'
        !         call util_crashpoint(669873)
        !         return
        ! end select
        end do

            ! print *, ' '
            ! print *, 'dQ '
            ! do  kk=1,max_branch_per_node
            !     print *, dQ(kk), elemR(JMidx+kk,er_Flowrate)
            ! end do

        !% --- update flowrates
        do kk=1,max_branch_per_node
            ! if ((elemSI(JMidx+kk,esi_JunctionBranch_Exists)     == oneI) .and. &
            !     (elemSI(JMidx+kk,esi_JunctionBranch_CanModifyQ) == oneI)         ) then
            !     elemR(JMidx+kk,er_Flowrate) = elemR(JMidx+kk,er_Flowrate) + dQ(kk)
            ! end if
            if (bFixYN(kk) .ne. zeroI) then
                elemR(JMidx+kk,er_Flowrate) = elemR(JMidx+kk,er_Flowrate) + dQ(kk)
            end if
        end do

        !% --- update overflow
        Qoverflow = Qoverflow + dQoverflow

        !% --- recompute the residual
        resid = junction_conservation_residual (JMidx)

            !  print *, 'NEW RESIDUAL ',resid 

            ! print *, abs(resid), localEpsilon

        if (abs(resid) > localEpsilon) then 
            print *, 'resid ',resid
            print *, 'CODE ERROR: unexpected flowrate residual'
            call util_crashpoint(6109873)
            return 
        end if


    end subroutine junction_fix_conservation
!%    
!%==========================================================================
!%==========================================================================
!% 
    !     real(8) pure function junction_source_term (JMidx, Qlat, Qover, QnetBranches) 
    !         !%------------------------------------------------------------------
    !         !% Description:
    !         !% computes the junction source term for quadratic solution
    !         !%------------------------------------------------------------------
    !             real(8), intent(in) :: Qlat, Qover, QnetBranches
    !         !%------------------------------------------------------------------

    !         junction_source_term = Qlat + Qover + Qnetbranches

    !     end function junction_source_term
    ! !%    
!%==========================================================================
!%==========================================================================
!% 
    ! subroutine junction_branch_dQdH (dQdH, JMidx, fm, kstart)    
    !     !%-----------------------------------------------------------------
    !     !% Description
    !     !% Computes dQdH for junction branches
    !     !%-----------------------------------------------------------------
    !         real(8), intent(inout) :: dQdH(:)
    !         integer, intent(in)    :: JMidx !% index of JM junction
    !         integer, intent(in)    :: fm(:) !% map up or down to face
    !         integer, intent(in)    :: kstart !% = 1 for upstream branches, 2 for down

    !         real(8), pointer :: fpsiL(:), fEnergyHead(:), eHead(:), fQ(:)
    !         integer ::  kk, tB
    !         real(8) ::  oneL = 1.d0

    !         real(8), parameter :: localEpsilon = 1.0d-6
    !     !%-----------------------------------------------------------------
    !     !% Aliases
    !         eHead       => elemR(:,er_Head)
    !         fEnergyHead => faceR(:,fr_EnergyHead)
    !         fpsiL       => faceR(:,fr_2B_psiL)
    !         fQ          => faceR(:,fr_Flowrate)
    !     !%-----------------------------------------------------------------

    !     !% HACK -- this needs to be cleaned up and made concurrent
    !     !%  and the subroutine should be pure.
        
    !     do kk=kstart,max_branch_per_node,2
    !         !% --- this branch index
    !         tB = JMidx + kk
    !         dQdH(tB) = zeroR !% ensure that even dummy branches have a value

    !         ! print *, ' '
    !         ! print *, '================================================='
    !         ! print *, 'kk here ',kk

    !         ! print *, 'branch exists ', elemSI(tB,esi_JunctionBranch_Exists)

    !         !% --- cycle if not a valid branch
    !         if (elemSI(tB,esi_JunctionBranch_Exists) .ne. oneI) cycle

    !         !% --- for adjacent channel/conduit branches
    !         if (elemSI(tB,esi_JunctionBranch_CC_adjacent) .eq. oneI) then

    !             ! print *, ' '
    !             ! print *, 'dif 1',fEnergyHead(fm(tB)) - eHead(JMidx)
    !             ! print *, 'fpsil',fpsiL(fm(tB))

    !             !% --- limiter for small values
    !             ! if ( (abs(fEnergyHead(fm(tB)) - eHead(JMidx)) < setting%ZeroValue%Depth) .or. &
    !             !      (abs(fpsiL(fm(tB))) < localEpsilon) ) then 
    !             !     dQdH(tB) = zeroR
    !             ! end if


    !             ! print *, ' '
    !             ! print *, 'tB', tB
    !             ! print *, 'fpsil   ',fpsiL(fm(tB))
    !             ! print *, 'head dif',(fEnergyHead(fm(tB)) - eHead(JMidx))

    !             !% --- first computation for dQdH using face 2 * beta * psi * L
    !             !%     along with face E and junction H
    !             ! dQdH(tB) = oneR / ( fpsiL(fm(tB)) * (fEnergyHead(fm(tB)) - eHead(JMidx)) ) 

    !             ! print *, ' '
    !             ! print *, 'tB ',tB
    !             ! print *, 'fpsiL    ',fpsiL(fm(tB))
    !             ! print *, 'head dif ',(fEnergyHead(fm(tB)) - eHead(JMidx))

    !             !% --- error checking -- first computation should be positive
    !             !%     if not, use the Chezy-Manning approach
    !             ! if (dQdH(tB) .le. zeroR) then 

    !             ! if (elemR(tB,er_FroudeNumber) .ge. 0.99d0) then 
    !             !     dQdH(tB) = zeroR

    !             ! else

    !                 ! print *, ' '
    !                 ! print *, 'tB ,kk ',tB, kk
    !                 ! print *, 'branchsign',branchsign(kk)
    !                 ! print *, 'area, hyrad          ',elemR(tB,er_Area), elemR(tB,er_HydRadius)
    !                 ! print *, 'MN, Lenght           ', elemR(tB,er_ManningsN) , elemR(tB,er_Length)
    !                 ! print *, 'delta E              ',fEnergyHead(fm(tB)) - eHead(JMidx)


    !                 if (eHead(JMidx) < elemR(tB,er_Zbottom)) then 
    !                     !% junction head is too low such that change does not affect branch flow
    !                     dQdH(tB) = zeroR

    !                     ! if (JMidx == 31) then 
    !                     !     print *, 'dQdH zero for tB  ',tB
    !                     !  end if
    !                 else
    !                     !% --- Use CM approach
    !                     if (abs(fEnergyHead(fm(tB)) - eHead(JMidx)) > localEpsilon) then
    !                         dQdH(tB) = - branchsign(kk)                                                 &
    !                             * elemR(tB,er_Area) * (elemR(tB,er_HydRadius)**twothirdR)               &
    !                             / (                                                                     &
    !                                 twoR * elemR(tB,er_ManningsN) * (sqrt(elemR(tB,er_Length)))         &
    !                                 * sqrt(abs(fEnergyHead(fm(tB)) - eHead(JMidx)))  &
    !                                 )

    !                         !  if (JMidx == 31) then 
    !                         !     print *, 'tb, dQdH   ',tB, dQdH(tB)
    !                         !  end if
    !                     else 
    !                         dQdH(tB) = zeroR
    !                     end if
    !                 end if
    !             ! end if
    !                 ! print *, 'CODE ERROR'
    !                 ! print *, 'unexpected negative value for dQdH first part'
    !                 ! call util_crashpoint(629784)
    !             ! else 
    !             !     !% --- second computation for dQdH
    !             !     dQdH(tB) = -branchsign(kk) * sign(sqrt(dQdH(tB)),fQ(fm(tB)))
    !             ! end if

                    

                    
    !                 ! !% --- zero all inflows
    !                 ! if ((elemR(tb,er_Flowrate) * branchsign(kk)) > zeroR) then 
    !                 !     dQdH(tB) = zeroR
    !                 ! end if

    !                 ! print *, ' tb, dQdH ',tB, dQdH(tB)

        

    !         else
    !             !% --- for adjacent diagnostic branches
    !             print *, 'code error - unfinished - diagnostic branch next to junction'
    !             call util_crashpoint(229874)
    !         end if

    !     end do

    ! end subroutine junction_branch_dQdH
!%    
!%==========================================================================
!%==========================================================================
!%
    ! real(8) pure function junction_main_branchdQdHsum (JMidx)
    !     !%------------------------------------------------------------------
    !     !% Description
    !     !% Computes the net dQdH terms junction from all the branches
    !     !% Note -- assumes that all branches that do not exist have zero 
    !     !% dQdH
    !     !%------------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in) :: JMidx
    !         integer :: k1, k2
    !     !%------------------------------------------------------------------
    !     k1 = JMidx + 1
    !     k2 = JMidx + max_branch_per_node

    !     junction_main_branchdQdHsum = sum(branchsign * elemSR(k1:k2,esr_JunctionBranch_dQdH)) 

    ! end function junction_main_branchdQdHsum
!%    
!%==========================================================================
!%==========================================================================
!% 
    ! real(8) function junction_main_dH &
    !     (JMidx, QnetBranches, Qoverflow, dQdHstorage, dQdHoverflow, dQdHsum)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Computes the change in head for the junction
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: JMidx
    !         real(8), intent(in) :: Qoverflow, QnetBranches
    !         real(8), intent(in) :: dQdHoverflow, dQdHsum, dQdHstorage
    !         real(8), pointer :: Qlat
    !         real(8) :: denominator

    !         real(8), parameter :: localEpsilon = 1.0d-16
    !     !%------------------------------------------------------------------

    !     denominator = dQdHstorage - dQdHoverflow - dQdHsum

    !     ! if (JMidx == 634) then
    !     !     print *, ' '
    !     !     print *, ' flowrates ',elemR(JMidx,er_FlowrateLateral)
    !     !     print *, Qoverflow, QnetBranches
    !     !     print *, ' '
    !     !     print *, dQdHstorage, dQdHoverflow, dQdHsum
    !     !     print *, 'denominator ', denominator
    !     !     print *, ' '
    !     ! end if

    !     if (abs(denominator) > localEpsilon) then
    !         junction_main_dH = (elemR(JMidx,er_FlowrateLateral) + Qoverflow + QnetBranches) & 
    !                      / denominator
    !     else
    !         junction_main_dH = zeroR 
    !     end if

    ! end function junction_main_dH
!%    
!%==========================================================================






!%==========================================================================
!% 
    subroutine junction_branchface_forceElemValue (JMidx)
        !%------------------------------------------------------------------
        !% Description:
        !% Forces the JB face value on the adjacent element
        !% AS WRITTEN 20230316 THIS CAUSES BLOWUP: DO NOT USE!
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: JMidx
            integer             :: ii, fidx, eUp, eDn
        !%------------------------------------------------------------------

        print *, 'CODE ERROR junction_branchface_forceElemValue should not be used'
        call util_crashpoint(5390874)

        do ii=1,max_branch_per_node,2
            if (elemSI(Jmidx+ii,esi_JunctionBranch_Exists) == oneI) then 
                fidx = elemI(JMidx+ii,ei_Mface_uL)
                eUp  = faceI(fidx,fi_Melem_uL)
                elemR(eUp,er_Flowrate) = faceR(fidx,fr_Flowrate)
            end if
        end do

        do ii=2,max_branch_per_node,2
            if (elemSI(Jmidx+ii,esi_JunctionBranch_Exists) == oneI) then 
                fidx = elemI(JMidx+ii,ei_Mface_dL)
                eDn  = faceI(fidx,fi_Melem_dL)
                elemR(eDn,er_Flowrate) = faceR(fidx,fr_Flowrate)
            end if
        end do

    end subroutine junction_branchface_forceElemValue
!%    
!%==========================================================================
!%==========================================================================
!% 





















!%==========================================================================
!%    
!     subroutine junction_get_main_data (jMainR, jMainI thisJM, Npack)
!         !%------------------------------------------------------------------
!         !% Description:
!         !% stores the junction main data (not associated with branches)
!         !%------------------------------------------------------------------
!         !% Declarations
!             real(8), intent(inout) :: jMainR(:,:)
!             integer, intent(inout) :: jMainI(:,:)
!             integer, intent(in)    :: thisJM(:), Npack
!         !%------------------------------------------------------------------
!         !%------------------------------------------------------------------

!         jMainI(1:Npack, jmi_Jtype) = elemSI(thisJM,esi_JunctionMain_Type)    

!         jMainR(1:Npack,jmr_AreaPlan) = elemSI(thisJM,esi_Storage_Plan_Area)
!         jMainR(1:Npack,jmr_cJ)       = zeroR  !% initialized later
!         jMainR(1:Npack,jmr_Hstart)   = elemR(thisJM,er_Head) !% starting head
!         jMainR(1:Npack,jmr_Hdelta)   = zeroR !% initial dH 
!         jMainR(1:Npack,jmr_Hresid)   = zeroR !% initialized later
!         jMainR(1:Npack,jmr_Qlat)     = elemR(thisJM,er_FlowrateLateral)
!         jMainR(1:Npack,jmr_Sc)       = zeroR  !% initialized later

!     end subroutine junction_get_main_data
! !%
! !%==========================================================================
! !%==========================================================================
! !%    
!     subroutine junction_get_branch_data (jBranchR, Hstart)
    !         !%-----------------------------------------------------------------
    !         !% Description:
    !         !% Takes the data from faceR and elemR arrays and stores in the
    !         !% jBranchR 3D array
    !         !%-----------------------------------------------------------------
    !         !% Declarations:
    !             real(8), intent(inout) :: jBranchR(:,:,:)
    !             real(8), intent(in)    :: Hstart(:)
    !             integer, intent(in)    :: thisJM(:)

    !             integer :: kk
    !         !%-----------------------------------------------------------------

    !         !% --- upstream faces
    !         do kk=1,max_branch_per_node,2
    !             where (elemSI(thisJM+kk,esi_JunctionBranch_Exists) ==1)
    !                 jBranchR(:,kk,jbr_beta)     = +oneR !% +1 for upstream beta
    !                 jBranchR(:,kk,jbr_Edelta)   = (faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Head_u)               &
    !                                             + (faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Velocity_u)**twoI) &
    !                                             / (twoR * setting%constant%gravity)) - Hstart
    !                 jBranchR(:,kk,jbr_flowsign) = zeroR !% initialized later                            
    !                 jBranchR(:,kk,jbr_psiL2)    = faceR(elemI(thisJM+kk,ei_Mface_uL),fr_psiL2)
    !                 jBranchR(:,kk,jbr_Q)        = faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Flowrate)
    !                 jBranchR(:,kk,jbr_Qdelta)   = zeroR !% initial zero difference
    !                 jBranchR(:,kk,jbr_Qresid)   = zeroR !% initialized later
                    
    !             elsewhere
    !                 jBranchR(:,kk,jbr_beta)    = zeroR
    !                 jBranchR(:,kk,jbr_Edelta)  = zeroR
    !                 jBranchR(:,kk,jbr_psiL2)   = zeroR
    !                 jBranchR(:,kk,jbr_Q)       = zeroR
    !                 jBranchR(:,kk,jbr_Qdelta)  = zeroR
    !                 jBranchR(:,kk,jbr_Qresid)  = zeroR
                    
    !             endwhere 
    !         end do

    !         !% --- downstream faces
    !         do kk=2,max_branch_per_node,2
    !             where (elemSI(thisJM+kk,esi_JunctionBranch_Exists) ==1)
    !                 jBranchR(:,kk,jbr_beta)    = -oneR !% -1 for downstream beta
    !                 !% --- note, Edelta needs Hstart 
    !                 jBranchR(:,kk,jbr_Edelta)  = (faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Head_d)             &
    !                                             + (faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Velocity_d)**twoI)  &
    !                                             / (twoR * setting%constant%gravity)) - Hstart
    !                 jBranchR(:,kk,jbr_flowsign) = zeroR !% initialized later  
    !                 jBranchR(:,kk,jbr_psiL2)   = faceR(elemI(thisJM+kk,ei_Mface_dL),fr_psiL2)
    !                 jBranchR(:,kk,jbr_Q)       = faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Flowrate)
    !                 jBranchR(:,kk,jbr_Qdelta)  = zeroR !% initial zero difference
    !                 jBranchR(:,kk,jbr_Qresid)  = zeroR !% initialized later
                    
    !             elsewhere
    !                 jBranchR(:,kk,jbr_beta)    = zeroR
    !                 jBranchR(:,kk,jbr_Edelta)  = zeroR
    !                 jBranchR(:,kk,jbr_psiL2)   = zeroR
    !                 jBranchR(:,kk,jbr_Q)       = zeroR
    !                 jBranchR(:,kk,jbr_Qdelta)  = zeroR
    !                 jBranchR(:,kk,jbr_Qresid)  = zeroR
                    
    !             endwhere 
    !         end do

!     end subroutine junction_get_branch_data
! !%
! !%==========================================================================
! !%==========================================================================
! !%  
!     function junction_flowsign (Flowrate, dQ, nk)
    !         !%-----------------------------------------------------------------
    !         !% Description
    !         !% returns +- oneR depending on flow sign
    !         !%-----------------------------------------------------------------
    !         !% Declarations:
    !             integer, intent(in)    :: nk
    !             real(8), intent(in)    :: Flowrate(:), dQ(:)
    !             real(8), dimension(nk) :: oneOut = oneR

    !             real(8), dimension(nk) :: junction_flowsign
    !         !%-----------------------------------------------------------------

    !         junction_flowsign = sign(oneOut,Flowrate + dQ)
            
!     end function junction_flowsign
! !%
! !%==========================================================================
! !%==========================================================================
! !%    
!     subroutine junction_cJ (AreaPlan, jMainI, thisJM)
    !         !%-----------------------------------------------------------------
    !         !% Description:
    !         !% Sets the jmr_cJ value that depends on storage type for all junctions
    !         !%-----------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(inout) :: AreaPlan(:)
    !             integer, intent(in)    :: jMainI(:,:), thisJM(:)
    !         !%-----------------------------------------------------------------

    !         where (jMainI(:,jmi_Jtype) .eq. ImpliedStorage)
    !             AreaPlan = zeroR
    !         elsewhere
    !             AreaPlan = elemSR(thisJM,esr_Storage_Plan_Area)
    !         endwhere

!     end subroutine junction_cJ
! !%
! !%==========================================================================
! !%==========================================================================
! !%  
!     subroutine junction_Sc (Sc, beta, Qinitial)
    !         !%-----------------------------------------------------------------
    !         !% Description
    !         !% Computes the source term for continuity as the net sum
    !         !% of the initila flowrate
    !         !%-----------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(inout) :: Sc(:)
    !             real(8), intent(in)    :: beta(:,:)
    !             real(8), intent(in)    :: Qinitial(:,:)
    !         !%-----------------------------------------------------------------
    !         !%-----------------------------------------------------------------

    !         Sc = sum(beta * Qinitial,2)    

!     end subroutine junction_Sc
! !%
! !%==========================================================================
! !%==========================================================================
! !%  
!     subroutine junction_Sm (jBranchR)
    !         !%-----------------------------------------------------------------
    !         !% Description:
    !         !% sets the juction source term
    !         !%-----------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(inout) :: jBranchR(:,:,:)
    !         !%-----------------------------------------------------------------

    !         where (jBranchR(:,:,jbr_beta)  .ne. zeroI)
    !             jBranchR(:,:,jbr_Sm)  = jBranchR(:,:,jbr_a) * (jBranchR(:,:,jbr_Qinit)**twoI) &
    !                                   - jBranchR(:,:,jbr_Edelta)
    !         endwhere 

!     end subroutine junction_Sm    
! !%
! !%==========================================================================
! !%==========================================================================
! !%  
!     subroutine junction_Qresid_initial (Qresid, Sm)
    !         !%-----------------------------------------------------------------
    !         !% Description
    !         !% Sets the initial residual for the flowrate on the branches
    !         !%-----------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(inout) :: Qresid(:,:)
    !             real(8), intent(in)    :: Sm(:,:)
    !         !%-----------------------------------------------------------------

    !         Qresid = Sm

!     end subroutine junction_Qresid_initial
! !%
! !%==========================================================================
! !%==========================================================================
! !% 
!     subroutine junction_Hresid_initial (Hresid, Sc)
    !         !%-----------------------------------------------------------------
    !         !% Description
    !         !% Sets the initial residual for the flowrate on the branches
    !         !%-----------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(inout) :: Hresid(:,:)
    !             real(8), intent(in)    :: Sc(:,:)
    !         !%-----------------------------------------------------------------

    !         Hresid = Sc

!     end subroutine junction_Hresid_initial
! !%
! !%==========================================================================

! !%                              OLD BELOW

! !%==========================================================================
! !%
!     subroutine  junction_calculationOLD (Npack, thisColP, istep)
    !         !%-----------------------------------------------------------------
    !         !% Description:
    !         !% get the jacobean matrix for junction element following the 
    !         !% derivation by ...
    !         !%-----------------------------------------------------------------
    !         !% Declarations
    !         integer, intent(in) :: Npack, thisColP, istep
    !         integer, pointer    :: thisJM(:)

    !         real(8), pointer :: grav

    !         !% --- the junction data storage
    !         real(8) :: jDataR(Npack, max_branch_per_node, NCol_jDataR) 
    !         integer :: jDataI(Npack, max_branch_per_node, NCol_jDataI)

    !         !% --- solver data storage
    !         real(8), dimension(Npack) :: startEJ, EJ, Qlat, VolumeEJ, fullVolumeEJ
    !         real(8), dimension(Npack) :: Zcrown, Aplan, Aponded
    !         real(8), dimension(Npack) :: OverflowLength, OverflowHeight, Qoverflow

    !         real(8), dimension(Npack) :: alpha, normL2, Qin, Qout,  fcoef
    !         real(8), dimension(Npack) :: deltaEJ, residEJ, Volume

    !         integer, dimension(Npack) :: np, Niter, OverFlowType, JunctionType

    !         logical, dimension(Npack) :: isconverged, isRepack, isfailed, isFull, isSubmergedOrifice
    !         integer :: kk, mm

    !         !% --- epsilon used in junction to prevent zero values
    !         real(8) :: jEpsilon = 1.d-6

    !         !% --- maximum number of iterations (HACK -- make into global)
    !         integer :: maxIter = 20
    !         !% --- flowrate conservation convergence (HACK -- make into global)
    !         real(8) :: jConvergence = 1.d-12

    !         !%-----------------------------------------------------------------------------
    !         !%  Aliases:         
    !             grav         => setting%constant%gravity
    !             !% --- packed set of junctions
    !             thisJM => elemP(1:Npack,thisColP)
    !         !%-----------------------------------------------------------------------------
    !         !% Preliminaries
    !             !% --- iteration counter
    !             Niter = zeroI

    !             !% --- initialization
    !             alpha   = zeroR
    !             normL2  = zeroR 
    !             Qin     = zeroR 
    !             Qout    = zeroR 
    !             fcoef   = zeroR
    !             deltaEJ = zeroR
    !             residEJ = zeroR
    !             Qoverflow = zeroR

    !             jDataR = zeroR
    !             jDataI = zeroI


    !             !% --- stores the initial elemR k location
    !             !%     this is used to bring packed branch data back to element and face array
    !             do concurrent (kk=1:max_branch_per_node)
    !                 jDataI(:,kk,ji_kidx)  = kk  
    !             end do

    !             !% store the starting head in Temp array for conservation check at end
    !             elemR(thisJM,er_Temp01) = elemR(thisJM,er_Head)
    !         !%-----------------------------------------------------------------------------

    !         !% --- get the junction data 
    !         call junction_get_JM_dataOLD (startEJ, EJ, Qlat, volumeEJ, fullVolumeEJ, &
    !                 Zcrown, Aplan, Aponded, OverflowLength, OverFlowHeight, &
    !                 isFull, isSubmergedOrifice, JunctionType, OverflowType, thisJM, Npack)

    !         ! print *, ' '
    !         ! print *, 'AT START OF JUNCTION SOLUTION FOR step ',istep
    !         ! print *, 'StartEJ ',startEJ(1), elemR(thisJM(1),er_Head)
    !         ! print *, 'EJ      ',EJ(1), elemR(thisJM(1),er_Head)
    !         ! print *, 'Zcrown  ',Zcrown(1), elemR(thisJM(1),er_Zcrown)
    !         ! print *, 'Zbottom ','                    ',elemR(thisJM(1),er_Zbottom)
    !         ! print *, 'Depth   ','                    ',elemR(thisJM(1),er_Depth)
    !         ! print *, 'QLat    ',Qlat(1), elemR(thisJM(1),er_FlowrateLateral)
    !         ! print *, 'Vej     ',volumeEJ(1), elemR(thisJM(1),er_Volume)
    !         ! print *, 'fullV   ',fullVolumeEJ(1), elemR(thisJM(1),er_FullVolume)
    !         ! print *, 'Apond   ', Aponded(1), elemSR(thisJM(1),esr_JunctionMain_PondedArea)
    !         ! print *, 'Aplan   ', Aplan(1),elemSR(thisJM(1),esr_Storage_Plan_Area)
    !         ! print *, 'Jtype   ', JunctionType(1), trim(reverseKey(JunctionType(1)))
    !         ! print *, 'overflowtype ',OverflowType(1), trim(reverseKey(OverflowType(1)))
    !         ! print *, ' '        

    !         !% --- store the face and elem data in 3D array
    !         call junction_get_face_data(jDataR, thisJM)

    !         ! print *, 'Face  DATA'
    !         ! !print *, 'size jDataR ',size(JdataR,1), size(JdataR,2), size(JdataR,3)
    !         ! do kk=1,max_branch_per_node
    !         !     !print *, 'kk ',kk
    !         !     if (JdataR(1,kk,jr_beta) .ne. zeroI) then 
    !         !         print *, ' '
    !         !         print *, kk,'beta  ',JdataR(1,kk,jr_beta)
                    
    !         !         if (mod(kk,2) == 0) then 
    !         !             print *, 'Area  ',JdataR(1,kk,jr_Area),    faceR(elemI(thisJM(1)+kk,ei_Mface_dL),fr_Area_d)
    !         !             print *, 'E     ',jDataR(1,kk,jr_Ebranch), faceR(elemI(thisJM(1)+kk,ei_Mface_dL),fr_Head_d)
    !         !             print *, 'length',jDataR(1,kk,jr_Length),  faceR(elemI(thisJM(1)+kk,ei_Mface_dL),fr_Length_d)
    !         !             print *, 'gamma ',jDataR(1,kk,jr_Gamma),   faceR(elemI(thisJM(1)+kk,ei_Mface_dL),fr_GammaM)
    !         !             print *, 'Q     ',jDataR(1,kk,jr_Q),       faceR(elemI(thisJM(1)+kk,ei_Mface_dL),fr_Flowrate)
    !         !             print *, 'K     ',jDataR(1,kk,jr_K),       faceR(elemI(thisJM(1)+kk,ei_Mface_dL),fr_KJunction_MinorLoss)
                        
    !         !         else 
    !         !             print *, 'Area  ',JdataR(1,kk,jr_Area),    faceR(elemI(thisJM(1)+kk,ei_Mface_uL),fr_Area_u)
    !         !             print *, 'E     ',jDataR(1,kk,jr_Ebranch), faceR(elemI(thisJM(1)+kk,ei_Mface_uL),fr_Head_u)
    !         !             print *, 'length',jDataR(1,kk,jr_Length),  faceR(elemI(thisJM(1)+kk,ei_Mface_uL),fr_Length_u)
    !         !             print *, 'gamma ',jDataR(1,kk,jr_Gamma),   faceR(elemI(thisJM(1)+kk,ei_Mface_uL),fr_GammaM)
    !         !             print *, 'Q     ',jDataR(1,kk,jr_Q),       faceR(elemI(thisJM(1)+kk,ei_Mface_uL),fr_Flowrate)
    !         !             print *, 'K     ',jDataR(1,kk,jr_K),       faceR(elemI(thisJM(1)+kk,ei_Mface_uL),fr_KJunction_MinorLoss)
    !         !         end if
    !         !     end if
    !         ! end do

    !         ! print *, ' '
    !         ! print *, 'heads '
    !         ! print *, elemR(faceI(elemI(thisJM(1)+1,ei_Mface_uL),fi_Melem_uL),er_Head)
    !         ! print *, faceR(elemI(thisJM(1)+1,ei_Mface_uL),fr_Head_u)
    !         ! print *, faceR(elemI(thisJM(1)+1,ei_Mface_uL),fr_Head_d)
    !         ! print *, elemR(thisJM(1),er_Head)
    !         ! print *, faceR(elemI(thisJM(1)+2,ei_Mface_dL),fr_Head_u)
    !         ! print *, faceR(elemI(thisJM(1)+2,ei_Mface_dL),fr_Head_d)
    !         ! print *, elemR(faceI(elemI(thisJM(1)+2,ei_Mface_dL),fi_Melem_dL),er_Head)

    !         ! print *, ' '
    !         ! print *, 'flowrates'
    !         ! print *, elemR(faceI(elemI(thisJM(1)+1,ei_Mface_uL),fi_Melem_uL),er_Flowrate)
    !         ! print *, faceR(elemI(thisJM(1)+1,ei_Mface_uL),fr_Flowrate)
    !         ! !print *, faceR(elemI(thisJM(1)+1,ei_Mface_uL),fr_Head_d)
    !         ! print *, elemR(thisJM(1),er_Flowrate)
    !         ! print *, faceR(elemI(thisJM(1)+2,ei_Mface_dL),fr_Flowrate)
    !         ! !print *, faceR(elemI(thisJM(1)+2,ei_Mface_dL),fr_Head_d)
    !         ! print *, elemR(faceI(elemI(thisJM(1)+2,ei_Mface_dL),fi_Melem_dL),er_Flowrate)

    !         !stop 498734

    !         !% --- set the alpha at the junction
    !         call junction_alpha (alpha, Aplan, Aponded, isFull, JunctionType, &
    !              OverflowType, thisJM, Npack, istep)

    !         ! print *, ' '
    !         ! print *, 'Alpha ',alpha(1)
    !         ! print *, ' '

    !         !% --- set the fcoef for overflows
    !         call junction_fcoef(fcoef, Aplan, OverflowLength, OverflowType)

    !         ! print *, 'Fcoef ',fcoef(1)
    !         ! print *, ' '

    !         !% AFTER THIS POINT WE SHOULD NOT NEED TO USE thisJM(:) UNTIL
    !         !% MOVING DATA BACK TO faceR AND elemR

    !         !% --- identify branches with gamma = 0 and K=0 or gamma = 0 and Q = 0
    !         !%     requires setting beta to zero (thus, we cannot compute flow in 
    !         !%     this branch during this time step)
    !         call junction_initialize_beta (jDataR, jEpsilon)

    !         ! print *, ' '
    !         ! print *, '=========================================='
    !         ! print *, 'Beta '
    !         ! print *, jDataR(1,:,jr_beta)
    !         ! print *, ' '
    !         ! print *, 'Gamma '
    !         ! print *, jDataR(1,:,jr_Gamma)
    !         ! print *, ' '
    !         ! print *, 'K '
    !         ! print *, jDataR(1,:,jr_K)
    !         ! print *, ' '
    !         ! print *, 'Q '
    !         ! print *, jDataR(1,:,jr_Q)
    !         ! print *, ' '

    !         !% --- the expected total branches
    !         !%     This will be modified by packing if there are branches with beta=0.
    !         np = elemSI(thisJM,esi_JunctionMain_Total_Branches)     

    !         ! print *, 'np ',np
    !         ! print *, ' '
            
    !         !% --- code check for debugging
    !         call junction_check_beta (jDataR, np, nPack)

    !         !% --- pack the jDataR and jDataI to remove zero branches
    !         isRepack(:) = .true.
    !         call junction_pack_all (jDataR, jDataI, np, isRepack, Npack)

    !         ! print *, 'repack Ebranch'
    !         ! print *, jDataR(1,:,jr_Ebranch)
    !         ! print *, ' '
    !         ! print *, 'repack Area '
    !         ! print *, jDataR(1,:,jr_Area)

    !         !% --- compute the invariant terms that will not change in this time step
    !         call junction_invariant_terms (jDataR, np, Npack)

    !         ! print *, 'jA '
    !         ! print *, jDataR(1,:,jr_a)
    !         ! print *, ' '
    !         ! print *, 'jC' 
    !         ! print *, jDataR(1,:,jr_c)
    !         ! print *, ' '
    !         ! print *, 'lambda A'
    !         ! print *, jDataR(1,:,jr_LambdaA)
    !         ! print *, ' '
    !         ! print *, 'lambda B'
    !         ! print *, jDataR(1,:,jr_LambdaB)
    !         ! print *, ' '


    !         !% --- compute the residual in each branch
    !         call junction_all_residuals (jDataR, residEJ, EJ, &
    !             startEJ, alpha, fcoef, Qlat, Zcrown, OverflowHeight, &
    !             OverflowType, isFull, isSubmergedOrifice, np, Npack)

    !         ! print *, 'residEJ BB', residEJ(1)
    !         ! print *, ' '
            
    !         !stop 298734

    !         isconverged = .false. 
    !         isfailed    = .false.

    !         !do concurrent (mm=1:Npack)
    !         do mm=1,1 !%Npack
    !                 ! print *, 'in junction loop ',mm
    !                 ! print *, 'resid at start   ',residEJ(mm)
    !             do while (.not. isconverged(mm))
    !                     ! print *, '============================================================='
    !                     ! print *, 'Niter ',Niter(mm), '===',mm,'=================================='
    !                     ! print *, ' '
    !                     !% --- increment iteration counter for exit after maxIter loops
    !                 Niter(mm) = Niter(mm) + 1

    !                 !% --- compute varying lambdaC 
    !                 jDataR(mm,:,jr_LambdaC) = junction_lambdaC (jDataR(mm,:,:), np(mm))

    !                 !% --- look for small LambdaC and repack if needed
    !                 call junction_check_repack_oneJ (jDataR(mm,:,:), np(mm), jEpsilon, isRepack(mm))
            
    !                 ! print *, 'lambda C 1',jDataR(mm,1:2,jr_LambdaC)

    !                 if (isRepack(mm)) then 
    !                     call junction_pack_oneJ &
    !                         (jDataR(mm,:,:), jDataI(mm,:,:), np(mm), isRepack(mm))
    !                 end if

    !                 ! print *, 'lambda C 2',jDataR(mm,1:2,jr_LambdaC)

    !                 !% --- compute the junction iteration
    !                 call junction_iteration &
    !                     (jDataR(mm,:,:), deltaEJ(mm), EJ(mm),  residEJ(mm), alpha(mm), &
    !                      fcoef(mm), Zcrown(mm), OverflowHeight(mm), OverflowType(mm), isFull(mm), &
    !                      isSubmergedOrifice(mm), np(mm)) 

    !                 ! print *,'new Q     ', jDataR(mm,1,jr_Q),    jDataR(mm,2,jr_Q)
    !                 ! print *, 'oldresid ', jDataR(mm,1,jr_residQ), jDataR(mm,2,jr_residQ)

    !                 !% --- compute the residuals
    !                 jDataR(mm,1:np(mm),jr_residQ) = junction_Qresidual              &
    !                                              (jDataR(mm,:,:), EJ(mm), startEJ(mm), np(mm),mm)
        
    !                 ! print *, 'newQ resid',jDataR(mm,1,jr_residQ) ,jDataR(mm,2,jr_residQ)                                     

    !                 residEJ(mm) = junction_Eresidual                                           &
    !                              (jDataR(mm,:,:), EJ(mm), startEJ(mm), alpha(mm), fcoef(mm),   &
    !                               Qlat(mm), Zcrown(mm), OverflowHeight(mm), isFull(mm),        &
    !                               isSubmergedOrifice(mm), OverflowType(mm), np(mm))

    !                 !% --- compute residual L2 norm
    !                 normL2(mm) = junction_get_L2norm &
    !                             (jDataR(mm,:,jr_residQ), residEJ(mm), np(mm))

    !                 ! print *, 'residEJ ',residEJ(mm)
    !                 ! print *, 'normL2  ',normL2(mm), jConvergence

    !                 !% --- check norms for exit
    !                 if (normL2(mm) .le. jConvergence) then 
    !                     isconverged(mm) = .true.
    !                     ! print *, ' ' 
    !                     ! print *, 'CONVERGED'
    !                     ! print *, ' '
    !                 else
    !                     if (Niter(mm)+1 .ge. maxIter) then 
    !                         !% --- mark as failed and exit
    !                         isconverged(mm) = .true.
    !                         isfailed(mm)    = .false.
    !                         !print *, ' '
    !                         print *, 'FAILED '
    !                         !print *, ' '
    !                     end if
    !                 end if 


                    
    !             end do

    !             !% ---- CONVERGED OR EXITING FOR THIS JUNCTION------------------------------
    !             !% --- compute the net Q inflow
    !             Qin(mm) = junction_branch_Qnet (jDataR(mm,:,:), np(mm), +oneI)
    !             !% --- compute the net Q outflow
    !             Qout(mm)= junction_branch_Qnet (jDataR(mm,:,:), np(mm), -oneI)
    !             !% --- compute the Q overflow
    !             Qoverflow(mm) = junction_Qoverflow (EJ(mm), Zcrown(mm), fcoef(mm),            &
    !                 OverflowHeight(mm), OverflowType(mm), isFull(mm), isSubmergedOrifice(mm) )

    !             !% --- adjust Qs to ensure conservation
    !             jDataR(mm,1:np(mm),jr_Q) = junction_conservation &
    !                 (jDataR(mm,:,:), residEJ(mm), Qin(mm), Qout(mm), np(mm))

    !             !% --- restore data to the faceR ... arrays
    !             call junction_push_face_flowrate &
    !                 (faceR(:,fr_Flowrate), faceR(:,fr_Velocity_u), faceR(:,fr_Velocity_d), &
    !                  faceR(:,fr_Area_u), faceR(:,fr_Area_d), jDataR(mm,:,:), jDataI(mm,:,:), &
    !                  np(mm), thisJM(mm))

    !             !% --- head
    !             elemR(thisJM(mm),er_Head) = EJ(mm)

    !             !% --- volume
    !             elemR(thisJM(mm),er_Volume) = junction_volume &
    !                 (elemR(thisJM(mm),er_Volume_N0), Qin(mm), Qout(mm), Qlat(mm), Qoverflow(mm), istep) 

    !             !% --- set the junction overflow and ponded volume on the last step of RK
    !             if (istep == 2) then
    !                 !% --- overflow
    !                 elemR(thisJM(mm),er_VolumeOverFlow) = Qoverflow(mm) * setting%Time%Hydraulics%Dt

    !                 !% --- ponded volume
    !                 if ((OverflowType(mm) == Ponded) .and. (isFull(mm))) then 
    !                     elemR(thisJM(mm),er_VolumePonded) = elemR(thisJM(mm),er_Volume) - elemR(thisJM(mm),er_FullVolume)
    !                 end if
    !             end if

    !         end do

    !         !% --- check conservation -- can be commented after debugging
    !         call junction_check_conservation (thisJM, Qoverflow, jDataR, Npack ,istep)

    !         !stop 6098734

    !         !% HACK --- the new water surface elevation associated with the
    !         !% changing storage may NOT exactly balance conservation because of
    !         !% the use of a fixed planar area. We need to rederiv the 
    !         !% junction_conservation equations so that the exact change in volume
    !         !% implied by the change in head gives us the correct volume

    !         !% --- HACK: NEED TO DO SOMETHING FOR REPORTING CONVERGENCE FAILURE 

!     end subroutine junction_calculationOLD
! !%
! !%========================================================================== 
! !%========================================================================== 
! !% 
!     pure subroutine junction_get_JM_dataOLD &
    !         (startEJ, EJ, Qlat, volumeEJ, fullVolumeEJ, Zcrown, Aplan, Aponded, &
    !             OverflowLength, OverflowHeight, &
    !             isFull, isSubmergedOrifice, JunctionType, OverflowType, thisJM, Npack)
    !         !%------------------------------------------------------------------
    !         !% Description:
    !         !% Gets the JM data from the elemR array and stores in local
    !         !% vector
    !         !%------------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(inout) :: startEJ(:), EJ(:), Qlat(:)
    !             real(8), intent(inout) :: volumeEJ(:), fullVolumeEJ(:), Zcrown(:)
    !             real(8), intent(inout) :: Aplan(:), Aponded(:)
    !             real(8), intent(inout) :: OverflowLength(:), OverflowHeight(:)
    !             logical, intent(inout) :: isFull(:), isSubmergedOrifice(:)
    !             integer, intent(inout) :: OverflowType(:), JunctionType(:)
    !             integer, intent(in)    :: Npack, thisJM(:)
    !         !%------------------------------------------------------------------

    !         startEJ     (1:Npack) = elemR(thisJM,er_Head)
    !         EJ          (1:Npack) = elemR(thisJM,er_Head)
    !         Qlat        (1:Npack) = elemR(thisJM,er_FlowrateLateral)
    !         volumeEJ    (1:Npack) = elemR(thisJM,er_Volume)
    !         fullVolumeEJ(1:Npack) = elemR(thisJM,er_FullVolume)
    !         Zcrown      (1:Npack) = elemR(thisJM,er_Zcrown)
    !         OverflowType(1:Npack) = elemSI(thisJM,esi_JunctionMain_OverflowType)
    !         Aponded     (1:Npack) = elemSR(thisJM,esr_JunctionMain_PondedArea)
    !         Aplan       (1:Npack) = elemSR(thisJM,esr_Storage_Plan_Area)
    !         JunctionType(1:Npack) = elemSI(thisJM,esi_JunctionMain_Type)

    !         OverflowLength(1:Npack) = elemSR(thisJM,esr_JunctionMain_OverflowOrifice_Length)
    !         OverflowHeight(1:Npack) = elemSR(thisJM,esr_JunctionMain_OverflowOrifice_Height)

    !         where (JunctionType .eq. ImpliedStorage)
    !             !% --- The only use of Aplan with implied storage is for the
    !             !%     overflow weir effect.
    !             Aplan = 0.179d0 !% HACK Area gives Lw = 1.5 m as overflow weir length
    !             where (elemR(thisJM,er_Head) .ge. Zcrown)
    !                 isFull = .true.
    !             elsewhere
    !                 isFull = .false.
    !             endwhere
    !         elsewhere
    !             where (volumeEJ .ge. fullVolumeEJ)
    !                 isFull = .true.
    !             elsewhere 
    !                 isFull = .false.
    !             endwhere
    !         endwhere 

    !         where ((isFull) .and. (OverflowType .eq. OverflowOrifice) &
    !             .and. (EJ > Zcrown + OverflowHeight))
    !             isSubmergedOrifice = .true.
    !         elsewhere 
    !             isSubmergedOrifice = .false.
    !         endwhere

   
!     end subroutine junction_get_JM_dataOLD
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%    
!     pure subroutine junction_get_face_data (jDataR, thisJM )    
    !         !%------------------------------------------------------------------
    !         !% Takes the data from faceR and elemR arrays and stores in the
    !         !% jDataR 3D array
    !         !%------------------------------------------------------------------
    !             real(8), intent(inout) :: jDataR(:,:,:)
    !             integer, intent(in)    :: thisJM(:)

    !             integer :: kk
    !         !%------------------------------------------------------------------
    !         !%------------------------------------------------------------------

    !         !% --- store data for upstream faces
    !         do concurrent (kk=1:max_branch_per_node:2)
    !             where (elemSI(thisJM+kk,esi_JunctionBranch_Exists) == oneI)
    !                 jDataR(:,kk,jr_beta)    = +oneR !% +1 for upstream beta
    !                 jDataR(:,kk,jr_Area)    =  faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Area_u)
    !                 jDataR(:,kk,jr_Ebranch) =  faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Head_u)            &
    !                                         + (faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Velocity_u)**twoI) &
    !                                         / (twoR * setting%constant%gravity)
    !                 jDataR(:,kk,jr_Gamma)   =  faceR(elemI(thisJM+kk,ei_Mface_uL),fr_GammaM)
    !                 jDataR(:,kk,jr_Length)  =  faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Length_u)
    !                 jDataR(:,kk,jr_Q)       =  faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Flowrate)
    !                 jDataR(:,kk,jr_K)       =  faceR(elemI(thisJM+kk,ei_Mface_uL),fr_KJunction_MinorLoss)
    !             elsewhere
    !                 jDataR(:,kk,jr_beta)    = zeroR !%  0 if branch does not exist
    !             endwhere
    !         end do

    !         !% --- store data for downstream faces
    !         do concurrent (kk=2:max_branch_per_node:2)
    !             where (elemSI(thisJM+kk,esi_JunctionBranch_Exists) == oneI)
    !                 jDataR(:,kk,jr_beta)    = -oneR !% -1 for downstream beta
    !                 jDataR(:,kk,jr_Area)    =  faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Area_d)
    !                 jDataR(:,kk,jr_Ebranch) =  faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Head_d)             &
    !                                         + (faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Velocity_d)**twoI)  &
    !                                         / (twoR * setting%constant%gravity)
    !                 jDataR(:,kk,jr_Gamma)   =  faceR(elemI(thisJM+kk,ei_Mface_dL),fr_GammaM)
    !                 jDataR(:,kk,jr_Length)  =  faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Length_d)
    !                 jDataR(:,kk,jr_Q)       =  faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Flowrate)
    !                 jDataR(:,kk,jr_K)       =  faceR(elemI(thisJM+kk,ei_Mface_dL),fr_KJunction_MinorLoss)
    !             elsewhere
    !                 jDataR(:,kk,jr_beta)    = zeroR !%  0 if branch does not exist
    !             endwhere
    !         end do


!     end subroutine junction_get_face_data
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%      
!     subroutine junction_alpha (alpha, Aplan, Aponded, isFull, JunctionType,  &
    !                 OverflowType, thisJM, Npack, istep)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% Computes the alpha term that represents change in junction storage
    !         !%------------------------------------------------------------------
    !         !% Declarations
                
    !             real(8), intent(inout) :: alpha(:)
    !             real(8), intent(in)    :: Aplan(:), Aponded(:)
    !             integer, intent(in)    :: JunctionType(:), OverflowType(:), thisJM(:), Npack, istep
    !             logical, intent(in)    :: isFull(:)

    !             real(8), pointer :: crk(:), dt
    !             integer          :: mm
    !         !%------------------------------------------------------------------
    !         !% Aliases
    !             crk          => setting%Solver%crk2
    !             dt           => setting%Time%Hydraulics%Dt
    !         !%------------------------------------------------------------------

    !         do mm=1,Npack
    !             select case (JunctionType(mm))

    !                 case (ImpliedStorage)
    !                     !% --- no volume increase allowed
    !                     alpha(mm) = zeroR

    !                 case (FunctionalStorage, TabularStorage)
    !                     if (.not. isFull(mm)) then
    !                         !% --- below full volume
    !                         alpha(mm) = Aplan(mm) / (crk(istep)*dt)
    !                     else
    !                         !% --- at or above full volume
    !                         select case (OverflowType(mm))
    !                             case (NoOverflow, OverflowOrifice, OverflowWeir)
    !                                 !% --- alpha remains based on plan area
    !                                 alpha(mm) = Aplan(mm) / (crk(istep)*dt)
    !                             case (Ponded)
    !                                 !% --- alpha used ponded area instead of plan area
    !                                 alpha(mm) = Aponded(mm)  / (crk(istep)*dt)
    !                             case default 
    !                                 print *, 'CODE ERROR: Unexpected case default'
    !                                 call util_crashpoint(2287422)
    !                         end select
    !                     end if
    !                 case default
    !                     print *, 'unexpected case default'
    !                     call util_crashpoint(628733)
    !             end select  
    !         end do

!     end subroutine junction_alpha
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     pure subroutine junction_fcoef (fcoef, Aplan, OverflowLength, OverflowType)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% Computes the fcoef term for overflow weirs
    !         !%------------------------------------------------------------------
    !         !% Declarations   
    !             real(8), intent(inout) :: fcoef(:)
    !             real(8), intent(in)    :: Aplan(:), OverflowLength(:)
    !             integer, intent(in)    :: OverflowType(:)
    !             real(8), parameter :: Cbc = 1.45d0  !% weir coef. HACK MOVE TO SETTINGS
    !         !%------------------------------------------------------------------

    !         where (OverflowType == OverflowWeir)
    !             fcoef = twoR * Cbc * sqrt( setting%Constant%pi * Aplan)
    !         elsewhere (OverflowType == OverflowOrifice)
    !             fcoef = OverflowLength * sqrt(twoR * setting%Constant%pi)
    !         elsewhere
    !             fcoef = zeroR
    !         endwhere

!     end subroutine junction_fcoef    
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     pure subroutine junction_initialize_beta (jDataR, jEpsilon)
    !         !%------------------------------------------------------------------
    !         !% Description:
    !         !% sets the initial jDataR(:,jr_beta) = 0 for branches that
    !         !% cannot have a flow in this time step.
    !         !%------------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(inout) :: jDataR(:,:,:)
    !             real(8), intent(in)    :: jEpsilon
    !         !%------------------------------------------------------------------
    !         !%------------------------------------------------------------------

    !         where (((jDataR(:,:,jr_Gamma)     < jEpsilon) .and.   &
    !                 (jDataR(:,:,jr_K)         < jEpsilon)       ) &
    !             .or.                                              &
    !                ((jDataR(:,:,jr_Gamma)     < jEpsilon) .and.   & 
    !                 (abs(jDataR(:,:,jr_Q))    < jEpsilon)       ) ) 

    !             jDataR(:,:,jr_beta) = zeroR

    !         endwhere

!     end subroutine junction_initialize_beta
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     subroutine junction_check_beta (jR, np, nPack)
            ! !%------------------------------------------------------------------
            ! !% Description:
            ! !% Checks that np(:) is identical to the number of beta <> 0
            ! !% This can be commented out once the code is thoroughly debugged
            ! !%------------------------------------------------------------------
            !     real(8), intent(in) :: jR(:,:,:)
            !     integer, intent(in) :: np(:), nPack 
            !     integer :: mm
            !     real(8) :: nbeta

            ! !%------------------------------------------------------------------
            ! !%------------------------------------------------------------------

            ! do mm=1,Npack 
            !     nbeta = count(jR(mm,:,jr_beta) .ne. zeroR)
            !     if (np(mm) .ne. nbeta) then 
            !         print *, 'CODE ERROR:'
            !         print *, 'problem with np and nbeta in junction iteration'
            !         print *, 'mm     = ',mm 
            !         print *, 'np(mm) = ',np(mm)
            !         print *, 'n beta = ',nbeta
            !         call util_crashpoint(398233)
            !     end if
            ! end do

!     end subroutine junction_check_beta
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%
!     pure subroutine junction_invariant_terms (jR, np, Npack)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% computes the terms in the junction solution that do not change
    !         !% during a time step
    !         !%------------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(inout) :: jR(:,:,:)
    !             integer, intent(in)    :: np(:), nPack

    !             integer :: mm
    !         !%------------------------------------------------------------------

    !         do concurrent (mm=1:Npack)
    !             !% --- a
    !             jR(mm,1:np(mm),jr_a) = jR(mm,1:np(mm),jr_beta)             &
    !                 * jR(mm,1:np(mm),jr_Length) * jR(mm,1:np(mm),jr_Gamma) &
    !                 / (twoR * setting%constant%gravity * jR(mm,1:np(mm),jr_Area))

    !             !% --- c
    !             jR(mm,1:np(mm),jr_c) = jR(mm,1:np(mm),jr_beta)     &
    !                 * jR(mm,1:np(mm),jr_K)                             &
    !                 / (twoR * setting%constant%gravity * (jR(mm,1:np(mm),jr_Area)**twoI))   

    !             !% --- lambdaA
    !             jR(mm,1:np(mm),jr_LambdaA)                            &
    !                 = twoR * setting%constant%gravity  *  (jR(mm,1:np(mm),jr_Area)**twoI) 

    !             !% --- lambdaB
    !             jR(             mm,1:np(mm),jr_LambdaB)               &
    !                 =        jR(mm,1:np(mm),jr_Area)                  &
    !                        * jR(mm,1:np(mm),jr_Length)                &
    !                        * jR(mm,1:np(mm),jr_Gamma)
    !         end do

!     end subroutine junction_invariant_terms
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     pure function junction_lambdaC (jR, np)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% computes the lambdaC term that changes with each iteration
    !         !%------------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(in)    :: jR(:,:)
    !             integer, intent(in)    :: np

    !             real(8), dimension(np) :: junction_lambdaC

    !             integer :: mm
    !         !%------------------------------------------------------------------

    !         !% --- note that the LambdaC storage is later used for Lambda = LambdaA/LambdaC
    !         junction_lambdaC = jR(1:np,jr_LambdaB) + jR(1:np,jr_K) * abs(jR(1:np,jr_Q))

!     end function junction_lambdaC
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     subroutine junction_iteration (jR, deltaEJ, EJ, residEJ, alpha, &
    !             fcoef, Zcrown, OverflowHeight, OverflowType, isFull, isSubmergedOrifice, np) 
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% Iterative results for one step advancing Q and EJ
    !         !%------------------------------------------------------------------
    !             real(8), intent(inout) :: jR(:,:), deltaEJ, EJ
    !             real(8), intent(in)    :: residEJ, alpha, fcoef, Zcrown, OverflowHeight
    !             logical, intent(in)    :: isFull, isSubmergedOrifice
    !             integer, intent(in)    :: np, OverFlowType

    !             real(8) :: Bsum, Csum, Gterm

    !             integer :: mm
    !         !%------------------------------------------------------------------
    !         !%------------------------------------------------------------------

    !         !% --- compute lambda and store in lambdaC (recomputed each iteration)
    !         jR(1:np,jr_LambdaC) =  jR(1:np,jr_LambdaA)  / jR(1:np,jr_LambdaC)

    !         ! print *, 'Lambda ratio ',jR(1:np,jr_LambdaC)
            
    !         !% B summation
    !         Bsum = sum(jR(1:np,jr_LambdaC)) 

    !         ! print *, 'Bsum ',Bsum
            
    !         !% C summation
    !         Csum = sum(jR(1:np,jr_LambdaC) * jR(1:np,jr_residQ))

    !         ! print *, 'Csum ',Csum

    !         !% Gterm
    !         Gterm = junction_Gterm (alpha, EJ, Zcrown, fcoef, OverflowHeight, &
    !                                 OverflowType, isFull, isSubmergedOrifice)

    !         ! print *, 'Gterm ',Gterm

    !         !% change in Junction Energy
    !         deltaEJ = (-residEJ + Csum ) / (-Gterm - Bsum)

    !         ! print *, 'deltaEJ ', deltaEJ

    !         ! print *, ' '
    !         ! print *, '========='
    !         ! print *, 'lambdaC ',jR(1:np,jr_LambdaC)
    !         ! print *, 'residQ  ',jR(1:np,jr_residQ)
    !         ! print *, 'deltaEJ ',deltaEJ
    !         ! print *, 'beta    ',jR(1:np,jr_beta)
        
    !         !% change in flowrate (recall LambdaC was overwritten with LambdaI)
    !         jR(1:np,jr_DeltaQ) = jR(1:np,jr_LambdaC) * (-jR(1:np,jr_residQ) - deltaEJ) &
    !                            / jR(1:np,jr_beta)


    !         ! print *, 'dQ ', jR(1:np,jr_DeltaQ) 
    !         ! print *, '==============='
    !         ! print *, ' '

    !         ! print *, 'old Q ',jR(1:np,jr_Q)

    !         !%  New Q flowrate
    !         jR(1:np,jr_Q) = jR(1:np,jr_Q) + jR(1:np,jr_DeltaQ)

    !         ! print *, 'new Q ',jR(1:np,jr_Q)
    !         ! print *, ' '
    !         ! print *, 'old EJ ', EJ

    !         !% --- update junction head
    !         EJ = EJ + deltaEJ

            
    !         ! print *, 'new EJ ',EJ
    !         ! print *, ' '

!     end subroutine junction_iteration   
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%    
!     subroutine junction_all_residuals &
    !             (jR, residEJ, EJ, startEJ, alpha, fcoef, Qlat, &
    !              Zcrown, OverflowHeight, OverflowType, isFull, isSubmergedOrifice, np, Npack)
    !         !%------------------------------------------------------------------
    !         !% Description:
    !         !% Computes the residuals of the nonlinear junction equations for
    !         !% all junctions
    !         !%------------------------------------------------------------------
    !             real(8), intent(inout) :: jR(:,:,:), residEJ(:)
    !             real(8), intent(in)    :: EJ(:), startEJ(:), alpha(:), Qlat(:)
    !             real(8), intent(in)    :: Zcrown(:), fcoef(:), OverflowHeight(:)
    !             logical, intent(in)    :: isSubmergedOrifice(:), isFull(:)
    !             integer, intent(in)    :: OverflowType(:), np(:), Npack 

    !             integer :: mm
    !         !%------------------------------------------------------------------

    !         !do concurrent (mm=1:Npack)
    !         do mm=1,Npack
    !             !print *, 'mm here ',mm, '---------------------;'
    !             !% --- get the flowrate residuals
    !             jR(mm,1:np(mm),jr_residQ) = junction_Qresidual                  &
    !                 (jR(mm,:,:), EJ(mm), startEJ(mm),  np(mm), mm)

    !             ! if (mm==1) then
    !             !     print *, 'eJ       ',eJ(mm)
    !             !     print *, 'startEJ  ',startEJ(mm)
    !             !     print *, 'alpha    ',alpha(mm)
    !             !     print *, 'fcoef    ',fcoef(mm)
    !             !     print *, 'Qlat     ',Qlat(mm)
    !             !     print *, 'Zcrown   ',Zcrown(mm)
    !             !     print *, 'OverflowH',OverflowHeight(mm)
    !             !     print *, 'isFull   ', isFull(mm)
    !             !     print *, 'isSubmerO', isSubmergedOrifice(mm)
    !             !     print *, 'overFType', OverflowType(mm)
    !             !     print *, 'np       ', np(mm)
    !             !     print *, ' '
    !             ! end if

    !             residEJ(mm) = junction_Eresidual                                &
    !                 (jR(mm,:,:), EJ(mm), startEJ(mm), alpha(mm), fcoef(mm),     &
    !                 Qlat(mm), Zcrown(mm), OverflowHeight(mm), isFull(mm),       &
    !                 isSubmergedOrifice(mm), OverflowType(mm), np(mm))

    !             ! if (mm==1) then 
    !             !print *, 'residEJ AA',residEJ(mm), mm
    !             ! end if
    !         end do
            
!     end subroutine junction_all_residuals
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     real(8) function junction_Eresidual &
    !             (jR, EJ, startEJ, alpha, fcoef, Qlat, Zcrown, OverflowHeight, &
    !              isFull, isSubmergedOrifice, OverflowType, np)  
    !         !%------------------------------------------------------------------
    !         !% Description:
    !         !% Computes the E residual of the nonlinear junction equation for
    !         !% a single junction
    !         !%------------------------------------------------------------------
    !             real(8), intent(in)    :: jR(:,:)
    !             real(8), intent(in)    :: EJ, startEJ, alpha, fcoef, Qlat, Zcrown
    !             real(8), intent(in)    :: OverflowHeight
    !             logical, intent(in)    :: isSubmergedOrifice, isFull
    !             integer, intent(in)    :: np, OverflowType
    !         !%------------------------------------------------------------------

    !         !% --- get the junction energy residual
    !         if (.not. isFull) then
    !             !% --- junction that is not overflowing
    !             junction_Eresidual = -alpha * (EJ - startEJ) + Qlat &
    !                                 + sum( (jR(1:np,jr_beta) * jR(1:np,jr_Q)) )  
    !         else 
    !             select case (OverflowType)

    !                 case (NoOverflow)
    !                     !% --- Zcrown becomes volume limit
    !                     junction_Eresidual = -alpha * (Zcrown - startEJ) + Qlat     &
    !                                     + sum( (jR(1:np,jr_beta) * jR(1:np,jr_Q)) )

    !                 case (OverFlowWeir)
    !                     !% --- requires additional outflow term
    !                     junction_Eresidual = -alpha * (EJ - startEJ) + Qlat            & 
    !                                     - twoR * fcoef * (EJ - Zcrown)**threehalfR     &
    !                                     + sum( (jR(1:np,jr_beta) * jR(1:np,jr_Q)) )

    !                 case (OverflowOrifice)
    !                     if (isSubmergedOrifice) then 
    !                         junction_Eresidual = -alpha * (Zcrown - startEJ)                                     &
    !                                         - twothirdR * fcoef * ((EJ -  Zcrown                  )**threehalfR) &
    !                                         + twothirdR * fcoef * ((EJ - (Zcrown + OverflowHeight))**threehalfR) &
    !                                         + sum( (jR(1:np,jr_beta) * jR(1:np,jr_Q)) )
    !                     else
    !                         junction_Eresidual = -alpha * (Zcrown - startEJ)                  &
    !                                         - twothirdR * fcoef * ((EJ - Zcrown)**threehalfR) &
    !                                         + sum( (jR(1:np,jr_beta) * jR(1:np,jr_Q)) )
    !                     end if
    !                 case (Ponded)
    !                     !% --- same residual as not full (alpha has been changed)
    !                     junction_Eresidual = -alpha * (EJ - startEJ) + Qlat          &
    !                                     + sum( (jR(1:np,jr_beta) * jR(1:np,jr_Q)) )  
                
    !             case default
    !                 !% --- change to non-pure to have error checking here
    !             end select
    !         end if

!     end function junction_Eresidual
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     function junction_Qresidual (jR, EJ, startEJ, np, mm)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% Computes the flowrate residual for all the non-zero branches
    !         !% of a single junction
    !         !%------------------------------------------------------------------
    !             real(8), intent(in) :: jR(:,:), EJ, startEJ
    !             integer, intent(in) :: np, mm
    !             real(8), dimension(np) :: junction_Qresidual
    !         !%------------------------------------------------------------------

    !         ! if (mm==1) then 
    !         !             print *, ' '
    !         ! print *, 'in junction Qresid '
    !         ! print *, 'j_a ', jR(1:np,jr_a)
    !         ! print *, 'Q   ', jR(1:np,jr_Q)
    !         ! print *, 'j_c ', jR(1:np,jr_c)
    !         ! print *, 'signedQ^2 ', sign((jR(1:np,jr_Q)**twoI),jR(1:np,jr_Q))
    !         ! print *, 'EJ      ',EJ
    !         ! print *, 'ejBranch ',jR(1:np,jr_Ebranch)
    !         ! end if
            
        
    !         junction_Qresidual(1:np) = jR(1:np,jr_a) * jR(1:np,jr_Q)        &
    !             + jR(1:np,jr_c) * sign((jR(1:np,jr_Q)**twoI),jR(1:np,jr_Q)) &
    !             + EJ - jR(1:np,jr_Ebranch)

            
    !             ! if (mm==1 )then 
    !             !      print * ,'Qresid ',junction_Qresidual(1:np)
    !             !     print *, ' '
    !             ! end if

!     end function junction_Qresidual
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     real(8) pure function junction_Gterm &
    !             (alpha, EJ, Zcrown, fcoef, OverflowHeight, OverflowType, isFull, isSubmergedOrifice)    
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% Computes G term for various overflow (or not overflow) conditions
    !         !%------------------------------------------------------------------
    !             real(8), intent(in) :: alpha, EJ, Zcrown, fcoef, OverflowHeight
    !             integer, intent(in) :: OverflowType
    !             logical, intent(in) :: isFull, isSubmergedOrifice
    !         !%------------------------------------------------------------------
    !         if (.not. isFull) then 
    !             !% --- EJ < Zcrown
    !             junction_Gterm = alpha
    !         else 
    !             !% --- overflowing conditions
    !             select case (OverflowType)

    !                 case (NoOverflow)
    !                     !% --- surcharge without overflow
    !                     junction_Gterm = zeroR

    !                 case (OverflowOrifice)
    !                     if (isSubmergedOrifice) then 
    !                         junction_Gterm = fcoef                           &
    !                             * (                                         &
    !                                 + sqrt(EJ - Zcrown)                     &
    !                                 - sqrt(EJ - (Zcrown + OverflowHeight))  &
    !                                 )
    !                     else
    !                         junction_Gterm = fcoef * sqrt(EJ - Zcrown)
    !                     end if

    !                 case (OverflowWeir)
    !                     junction_Gterm = alpha + threeR * fcoef * sqrt(EJ - Zcrown)

    !                 case (Ponded)
    !                     !% --- note alpha is from the ponded area
    !                     junction_Gterm = alpha
    !                 case default 
    !                     !% --- code should not be here
    !                     !%     to have a failure for testing 
    !                     !%     the function must be unpure    
    !             end select      
    !         end if

!     end function junction_Gterm
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     ! real(8) pure function junction_Qoverflow_residual &
    !     !         (fcoef, EJ, Zcrown, Qoverflow)    
    !     !     !%------------------------------------------------------------------
    !     !     !% Description
    !     !     !% Computes the overflow flowrate residual\
    !     !     !%------------------------------------------------------------------
    !     !         real(8), intent(in) :: EJ, Zcrown, fcoef, Qoverflow
    !     !     !%------------------------------------------------------------------

    !     !     junction_Qoverflow_residual = fcoef * ((EJ - Zcrown)**(1.5d0)) + Qoverflow

!     ! end function junction_Qoverflow_residual
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     real(8) function junction_get_L2norm (residQ, residEJ, np)
    !         !%------------------------------------------------------------------
    !         !% Description:
    !         !% Computes the L2 norm of residQ for multiple branches and the
    !         !% single EJ for one junction
    !         !%------------------------------------------------------------------
    !             real(8), intent(in) :: residQ(:), residEJ
    !             integer, intent(in) :: np
    !         !%------------------------------------------------------------------

    !         ! print *, ' '
    !         ! print *, 'residQ, residEJ ',residQ(1:np), residEJ
    !         ! print *, 'sum    ' ,sum(residQ(1:np)**twoI)
    !         ! print *, 'square ',residEJ**twoI

    !         junction_get_L2norm = sqrt(sum(residQ(1:np)**twoI) + residEJ**twoI)

    !         !stop 5098734

!     end function junction_get_L2norm
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%      
!     real(8) pure function junction_branch_Qnet (jR, np, idir)
    !         !%------------------------------------------------------------------
    !         !% Description:
    !         !% Computes the net Q in (idir = 1) or net Q out (idir = -1)
    !         !% for a junction
    !         !%------------------------------------------------------------------
    !             real(8), intent(in) :: jR(:,:)
    !             integer, intent(in) :: np, idir
    !         !%------------------------------------------------------------------

    !         junction_branch_Qnet = onehalfR * sum                                      &
    !           (                                                                 &
    !             jR(1:np,jr_beta) * jR(1:np,jr_Q)                                &
    !             * (                                                             &
    !                 oneR + real(idir,8) * jR(1:np,jr_beta) * jR(1:np,jr_Q)      & 
    !                        / abs(jR(1:np,jr_Q))                                 &   
    !               )                                                             &
    !           )

!     end function junction_branch_Qnet
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%  
!     real(8) pure function junction_Qoverflow &
    !             (EJ, Zcrown, fcoef, OverflowHeight, OverflowType, &
    !              isFull, isSubmergedOrifice)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% Computes the overflow based on theory
    !         !%------------------------------------------------------------------
    !             real(8), intent(in) :: EJ, Zcrown, fcoef, OverflowHeight
    !             integer, intent(in) :: OverflowType
    !             logical, intent(in) :: isFull, isSubmergedOrifice
    !         !%------------------------------------------------------------------

    !         !% --- baseline is zero overflow
    !         junction_Qoverflow = zeroR
    !         if (.not. isFull) return 

    !         select case (OverflowType)
    !             case (NoOverflow)
    !                 !% --- retain zero
    !             case (OverflowOrifice)
    !                 if (isSubmergedOrifice) then 
    !                     junction_Qoverflow = - twothirdR * fcoef                        &
    !                         * (                                                &
    !                             ((EJ - Zcrown)**threehalfR)                    &
    !                            -((EJ - (Zcrown + OverflowHeight))**threehalfR) &
    !                         )

    !                 else
    !                     junction_Qoverflow = - twothirdR * fcoef * ((EJ - Zcrown)**(threehalfR))
    !                 end if
    !             case (OverflowWeir)
    !                 junction_Qoverflow = -twoR * fcoef * ((EJ - Zcrown)**(threehalfR))
    !             case (Ponded)
    !                 !% --- overflow volume is held in JM volume
    !                 !%     retain zero
    !             case Default
    !                 !% --- code should not be here
    !                 !%     to have a failure for testing 
    !                 !%     the function must be unpure    
    !         end select

!     end function junction_Qoverflow
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     pure function junction_conservation (jR, residEJ, QinTotal, QoutTotal, np)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% Proportionally adjust inflow/outflow Q to ensure flowrate
    !         !% and storage residual is zero to machine accuracy
    !         !% Applies to a single junction
    !         !%------------------------------------------------------------------
    !             real(8), intent(in)     :: jR(:,:), QInTotal, QoutTotal, residEJ
    !             integer, intent(in)     :: np

    !             real(8), dimension(np)  :: junction_conservation
    !         !%------------------------------------------------------------------

    !         junction_conservation = jR(1:np,jr_Q) &
    !             * (                                                                                  &
    !                oneR - onefourthR * residEJ                                                       &
    !                * (                                                                               &
    !                     (oneR + (jR(1:np,jr_beta) * jR(1:np,jr_Q)) / abs(jR(1:np,jr_Q)) ) / QinTotal  &
    !                    +                                                                             &
    !                     (oneR - (jR(1:np,jr_beta) * jR(1:np,jr_Q)) / abs(jR(1:np,jr_Q)) ) / QoutTotal &
    !                  )                                                                               &
    !                )
                
!     end function junction_conservation
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     pure subroutine junction_push_face_flowrate &
    !             (faceQ, faceVup, faceVdn, faceAup, faceAdn, jR, jI,  np, JMidx)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% Cycles through the solved branches (1:np) and restores their
    !         !% Q data to the faceR array 
    !         !%------------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(inout) :: faceQ(:),   faceVup(:), faceVdn(:)
    !             real(8), intent(in)    :: faceAup(:), faceAdn(:)
    !             real(8), intent(in)    :: jR(:,:)
    !             integer, intent(in)    :: np, JMidx
    !             integer, target, intent(in) ::  jI(:,:)

    !             integer :: mm, ii
    !             !integer, pointer :: tM, kidx, fIdx
    !             integer :: kidx(np), fIdx(np)
    !         !%------------------------------------------------------------------
    !         !%------------------------------------------------------------------

    !         !% --- cycle through the branches that were solved
    !         do mm=1,np
                
    !             !% --- branch kk in the 1:max_branch_per_node
    !             kidx(mm) = jI(mm,ji_kidx)
                
    !             !% --- handle upstream or downstream separately
    !             if (mod(kidx(mm)+1,twoI) .eq. zeroI) then 
    !                 !% --- upstream branch face
    !                 fIdx(mm) = elemI(JMidx+kidx(mm),ei_Mface_uL)
    !             else
    !                 !% --- downstream branch face
    !                 fIdx(mm) = elemI(JMidx+kidx(mm),ei_Mface_dL)
    !             end if
    !             !% --- set the flowrate
    !             faceQ(fIdx(mm))   = jR(mm,jr_Q)
    !             !% --- get consistent velocities (assumes area > epsilon)
    !             faceVup(fIdx(mm)) = faceQ(fIdx(mm)) / faceAup(fIdx(mm))
    !             faceVdn(fIdx(mm)) = faceQ(fIdx(mm)) / faceAdn(fIdx(mm))
    !         end do

!     end subroutine junction_push_face_flowrate   
! !%
! !%========================================================================== 
! !%========================================================================== 
! !% 
!     real(8) pure function junction_volume &
    !             (oldVolume, Qin, Qout, Qlat, Qoverflow, istep)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% Updates junction volume from Volume_N0 for an RK step
    !         !%------------------------------------------------------------------
    !         !% Declarations:
    !             real(8), intent(in) :: oldVolume, Qin, Qout, Qlat, Qoverflow
    !             integer, intent(in) :: istep
    !         !%------------------------------------------------------------------

    !         junction_volume = oldVolume                    &
    !             + setting%Solver%crk2(istep) * setting%Time%Hydraulics%Dt  &
    !             * (Qin - Qout + Qlat + Qoverflow)

!     end function junction_volume 
! !%
! !%========================================================================== 
! !%========================================================================== 
! !% 
!     subroutine junction_check_conservation (thisJM, Qoverflow, jR, Npack, istep)
    !         !%------------------------------------------------------------------
    !         !% Description:
    !         !% Checks flowrate and junction storage conservation
    !         !%------------------------------------------------------------------
    !         !% Declarations
    !             integer, target, intent(in) :: thisJM(:)
    !             integer, intent(in) :: Npack, istep
    !             real(8), intent(in) :: Qoverflow(:), jR(:,:,:)
    !             integer, pointer :: tM, upF, dnF
    !             integer          :: mm, kk
    !             real(8)          :: thisCons, thisVol, thisFlow
    !             character(64) :: subroutine_name = "junction_check_conservation"
    !         !%------------------------------------------------------------------
    !             ! print *, ' '
    !             ! print *, 'in ',trim(subroutine_name)

    !         do mm=1,Npack
    !             tM => thisJM(mm)

    !             ! print *, 'junction index ',tM
    !             ! print *, 'Qlateral ',elemR(tM,er_FlowrateLateral)
    !             ! print *, 'Qoverflow',Qoverflow(mm)
                
    !             !% --- lateral inflow
    !             thisCons = elemR(tM,er_FlowrateLateral) + Qoverflow(mm)
    !             thisFlow = abs(elemR(tM,er_FlowrateLateral))
                
    !             ! print *, 'thisCons 1',thisCons

    !             !% --- change in storage
    !             !%     assumes old head stored in er_Temp01
    !             select case (elemSI(tM,esi_JunctionMain_Type))
    !                 case (ImpliedStorage)
    !                     ! print *, 'in implied storage'
    !                     thisVol = sum(jR(mm,:,jr_Length) * jr(mm,:,jr_Area) * abs(jR(mm,:,jr_beta))) / twoR
    !                     !% --- no change to conservation by storage
    !                 case (TabularStorage,FunctionalStorage)
    !                     thisCons = thisCons &
    !                          - (elemSR(tM,esr_Storage_Plan_Area)    &
    !                               / (  setting%Solver%crk2(istep)    &
    !                                   * setting%Time%Hydraulics%Dt ) &
    !                            )                                     &
    !                             * (elemR(tM,er_Head) - elemR(tM,er_Temp01))
    !                     thisVol = elemR(tM,er_Volume)
    !                 case default
    !                     print *, 'CODE ERROR: unexpected case default'
    !             end select

    !             ! print *, 'thisCons 2',thisCons
    !             ! print *, 'thisVol  2',thisVol

    !             !% --- upstream branches
    !             do kk=1,max_branch_per_node,2
    !                 if (elemSI(tM+kk,esi_JunctionBranch_Exists) == oneI) then
    !                     upF => elemI(tM+kk,ei_Mface_uL)
    !                     thisCons = thisCons + faceR(upF,fr_Flowrate)
    !                     !print *, 'thisCons 3',thisCons
    !                     thisFlow = thisFlow + abs(faceR(upF,fr_Flowrate))
    !                 end if
    !             end do

    !             !% --- downstream branches
    !             do kk=2,max_branch_per_node,2
    !                 if (elemSI(tM+kk,esi_JunctionBranch_Exists) == oneI) then
    !                     dnF => elemI(tM+kk,ei_Mface_dL)
    !                     thisCons = thisCons - faceR(dnF,fr_Flowrate)
    !                     !print *, 'thisCons 4',thisCons
    !                     thisFlow = thisFlow + abs(faceR(dnF,fr_Flowrate))
    !                 end if
    !             end do

    !             ! print *, 'thisCons ',thisCons
    !             ! print *, 'thisFlow ',thisFlow
    !             ! print *, 'thisVol  ',thisVol
    !             ! print *, 'dt       ',setting%Time%Hydraulics%Dt

    !             !% volume based normalization
    !             if (thisVol > setting%ZeroValue%Volume) then
    !                 thisCons = thisCons * setting%Time%Hydraulics%Dt / thisVol
    !             else
    !                 !% flowrate based normalization
    !                 if (thisFlow > setting%ZeroValue%Volume) then
    !                     thisCons = thisCons / thisFlow
    !                 else 
    !                     !% no normalization
    !                 end if
    !             end if

    !             if (abs(thisCons) > 1.d-2) then 
    !                 print *, 'conservation error at ',thisJM(mm)
    !                 print *, 'magnitude ',thisCons
    !                 stop 5509873
    !             end if
    !         end do
            
!     end subroutine junction_check_conservation
! !%
! !%========================================================================== 
! !% JUNCTION PACKING BELOW
! !%========================================================================== 
! !%
!     pure subroutine junction_check_repack_all (isRepack, jR, np, Npack, jEpsilon)
!         !%------------------------------------------------------------------
!         !% Description
!         !% Checks if any branches are now zero and need repacking and 
!         !% adjusts beta to zero where branches are zero
!         !%------------------------------------------------------------------
!         !% Declarations
!             logical, intent(inout) :: isRepack(:)
!             real(8), intent(inout) :: jR(:,:,:)
!             integer, intent(in)    :: np(:), Npack
!             real(8), intent(in)    :: jEpsilon

!             integer :: mm
!         !%------------------------------------------------------------------
!         !%------------------------------------------------------------------
!         do concurrent (mm=1:Npack)

!             call junction_check_repack_oneJ (jR(mm,:,:), np(mm), jEpsilon, isRepack(mm))

!         end do

!     end subroutine junction_check_repack_all   
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%    
!     pure subroutine junction_check_repack_oneJ (jR, np, jEpsilon, isRepack)
!         !%------------------------------------------------------------------
!         !% Checks for zeros and resets beta if repack needed
!         !%------------------------------------------------------------------
!         !% Declarations
!             real(8), intent(inout) :: jR(:,:)
!             real(8), intent(in)    :: jEpsilon
!             integer, intent(in)    :: np 
!             logical, intent(inout) :: isRepack
!         !%------------------------------------------------------------------
!         !%------------------------------------------------------------------

!         if (any(jR(1:np,jr_LambdaC) < jEpsilon)) then 
!             isRepack = .true.
!             where (jR(1:np,jr_LambdaC) < jEpsilon)
!                 jR(1:np,jr_beta) = zeroR
!             endwhere
!         else
!             isRepack = .false.
!         end if

!     end subroutine junction_check_repack_oneJ
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%
!     pure subroutine junction_pack_all &
!         (jR, jI, np, isRepack, Npack)
!         !%------------------------------------------------------------------
!         !% Description
!         !% Packing on the 2nd dimension of 3D arrays, jR, jI
!         !% Cycles through the 1st dimension (packed JM index) in this procedure 
!         !% and then cycles thorugh through the 3rd dimension (data type)
!         !% in call to subsidiary procedures
!         !%------------------------------------------------------------------
!         !% Declarations
!             integer, intent(inout) :: np(:), jI(:,:,:)
!             real(8), intent(inout) :: jR(:,:,:)
!             logical, intent(inout) :: isRepack(:)
!             integer, intent(in)    :: Npack
        
!             integer :: mm
!         !%------------------------------------------------------------------
!         !% --- pack the jDataR and jDataI to remove zero branches
!         do concurrent (mm=1:Npack)
!             call junction_pack_oneJ(jR(mm,:,:), jI(mm,:,:), np(mm), isRepack(mm))
!         end do

!     end subroutine junction_pack_all
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     pure subroutine junction_pack_oneJ (jR, jI, np, isRepack)
!         !%------------------------------------------------------------------
!         !% Description
!         !% Packs data for one junction during implicit junction solution
!         !%------------------------------------------------------------------
!         !% Declarations
!             real(8), intent(inout) :: jR(:,:)
!             integer, intent(inout) :: jI(:,:), np
!             logical, intent(inout) :: isRepack
!         !%------------------------------------------------------------------
!         !%------------------------------------------------------------------
!         if (.not. isRepack) return 
    
!         !% --- count the nonzero branches
!         np = count(jR(:,jr_beta) .ne. zeroR)

!         !% --- perform the packs based on jr_beta (must delay packing jr_beta)
!         call junction_pack_all_jRtype (jR(:,:), np)
!         call junction_pack_all_jItype (jI(:,:), jR(:,:), np)

!         !% --- pack beta (always must be last)
!         jR(1:np,jr_beta)  = junction_pack_one_arrayR &
!                 (jR(:,jr_beta), jR(:,jr_beta), np)

!         !% --- repack finished, reset the logical
!         isRepack = .false.

!     end subroutine junction_pack_oneJ
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%        
!     pure subroutine junction_pack_all_jRtype (jR, np)  
!         !%------------------------------------------------------------------
!         !% Description:
!         !% Packs the first dimension of a 2D real array (jR) by cycling
!         !% through the data types in the 2nd dimension.
!         !% Note that np must equal the number of jR(:,jr_beta) that are 
!         !% nonzero
!         !%------------------------------------------------------------------
!         !% Declarations:
!             real(8), intent(inout) :: jR(:,:)
!             integer, intent(in)    :: np

!             !% --- jData elements that are packed
!             integer, parameter :: NjRset = 13
!             integer :: jRset(NjRset)
!             integer :: jj
!         !%------------------------------------------------------------------

!         !% --- define data sets used in packing 
!         !%     DO NOT INCLUDE jr_beta -- must be repacked last in a separate 
!         !%     call since it is themask used for packing.
!         jRset = [jr_Area, jr_a, jr_c, jr_Ebranch, jr_Gamma, jr_K, jr_Length, jr_Q, &
!                  jr_DeltaQ, jr_residQ, jr_LambdaA, jr_LambdaB, jr_LambdaC]

!         do concurrent (jj=1:NjRset)
!             !% --- pack the jRset(jj) data type for the thisJM(mm) branch
!             !%     to remove non-existent branches
!             jR(1:np,jRset(jj))  = junction_pack_one_arrayR                               &
!                         (jR(:,jRset(jj)), jR(:,jr_beta), np)
!         end do

!     end subroutine junction_pack_all_jRtype
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%        
!     pure subroutine junction_pack_all_jItype (jI, jR, np)  
!         !%------------------------------------------------------------------
!         !% Description:
!         !% Packs the first dimension of a 2D integer array (jI) by cycling
!         !% through the data types in the 2nd dimension.
!         !% Note that np must equal the number of jR(:,jr_beta) that are 
!         !% nonzero
!         !%------------------------------------------------------------------
!             integer, intent(inout) :: jI(:,:)
!             real(8), intent(in)    :: jR(:,:)
!             integer, intent(in)    :: np

!             integer, parameter :: NjIset = Ncol_jDataI
!             integer            :: jIset(NjIset)
!             integer            :: jj
!         !%------------------------------------------------------------------
!         !%------------------------------------------------------------------

!         !% --- define data sets used in packing
!         jIset = [ji_kidx]

!         do concurrent (jj=1:NjIset)
!             !% --- pack the jIset(jj) data type for the thisJM(mm) branch
!             !%     to remove non-existent branches
!             jI(1:np,jIset(jj))  = junction_pack_one_arrayI                               &
!                         (jI(:,jIset(jj)), jR(:,jr_beta), np)
!         end do

!     end subroutine junction_pack_all_jItype
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%      
!     pure function junction_pack_one_arrayR (inarray, beta, np)
!         !%------------------------------------------------------------------
!         !% Description
!         !% standard pack to remove any locations where beta = 0.0
!         !%------------------------------------------------------------------
!             real(8), intent(in)    :: inarray(:), beta(:)
!             integer, intent(in)    :: np
!             real(8), dimension(np) :: junction_pack_one_arrayR
!         !%------------------------------------------------------------------

!         junction_pack_one_arrayR(1:np) = pack (inarray, beta .ne. zeroR)

!     end function junction_pack_one_arrayR
! !%
! !%==========================================================================  
! !%========================================================================== 
! !%      
!     pure function junction_pack_one_arrayI (inarray, beta, np)
!         !%------------------------------------------------------------------
!         !% Description
!         !% standard pack to remove any locations where beta = 0.0
!         !%------------------------------------------------------------------
!             real(8), intent(in)    ::  beta(:)
!             integer, intent(in)    :: np, inarray(:)
!             integer, dimension(np) :: junction_pack_one_arrayI
!         !%------------------------------------------------------------------

!         junction_pack_one_arrayI(1:np) = pack (inarray, beta .ne. zeroR)

!     end function junction_pack_one_arrayI
! !%
! !%==========================================================================  
!     !%========================================================================== 
! !%
!     function matinv(A) result (invA)
!         !%-----------------------------------------------------------------
!         !% Description:
!         !% Find the inverse of an input matrix A
!         !%-----------------------------------------------------------------
!         !% Declarations
!         real(8),allocatable, intent(in) :: A(:,:)
!         real(8),allocatable             :: invA(:,:)
!         integer                         :: n   
!         !%-----------------------------------------------------------------
!         !% find the size of the matrix
!         n = size(A,1)

!         !% for square matrix size or 2-4, direct inversion is the fastest 
!         if (n == 2) then
!             invA = matinv2(A)
!         else if (n == 3) then
!             invA = matinv3(A)
!         else if (n == 4) then
!             invA = matinv4(A)
!         else
!         !% use the LAPACK for higher order matrices
!         !% talk with cesar about LAPACK installation recipie
!             print*, 'Matrix inversion for shape > 4 has not yet been developed'
!             stop 789451
!         end if     

!     end function 
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%
!     pure function matinv2(A) result(B)
!         !! Performs a direct calculation of the inverse of a 2×2 matrix.
!         real(8), intent(in) :: A(2,2)   !! Matrix
!         real(8)             :: B(2,2)   !! Inverse matrix
!         real(8)             :: detinv

!         ! Calculate the inverse determinant of the matrix
!         detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

!         ! Calculate the inverse of the matrix
!         B(1,1) = +detinv * A(2,2)
!         B(2,1) = -detinv * A(2,1)
!         B(1,2) = -detinv * A(1,2)
!         B(2,2) = +detinv * A(1,1)
!     end function
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%
!     pure function matinv3(A) result(B)
!         !! Performs a direct calculation of the inverse of a 3×3 matrix.
!         real(8), intent(in) :: A(3,3)   !! Matrix
!         real(8)             :: B(3,3)   !! Inverse matrix
!         real(8)             :: detinv

!         ! Calculate the inverse determinant of the matrix
!         detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
!                 - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
!                 + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

!         ! Calculate the inverse of the matrix
!         B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
!         B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
!         B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
!         B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
!         B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
!         B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
!         B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
!         B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
!         B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
!     end function
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%
!     pure function matinv4(A) result(B)
!         !! Performs a direct calculation of the inverse of a 4×4 matrix.
!         real(8), intent(in) :: A(4,4)   !! Matrix
!         real(8)             :: B(4,4)   !! Inverse matrix
!         real(8)             :: detinv

!         ! Calculate the inverse determinant of the matrix
!         detinv = &
!         1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
!         - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
!         + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
!         - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

!         ! Calculate the inverse of the matrix
!         B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
!         B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
!         B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
!         B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
!         B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
!         B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
!         B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
!         B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
!         B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
!         B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
!         B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
!         B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
!         B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
!         B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
!         B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
!         B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
!     end function 
!%
!%==========================================================================
    !%==========================================================================
!% 
    ! subroutine  junction_calculation (Npack, thisColP, istep)
    !     !%-----------------------------------------------------------------
    !     !% Description:
    !     !% get the jacobean matrix for junction element following the 
    !     !% derivation by ...
    !     !%-----------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in) :: Npack, thisColP, istep
    !         integer, pointer    :: thisJM(:)

    !         real(8) :: QnetBranches, dQdHsum, dQdHstorage, Qoverflow, dQdHoverflow
    !         real(8) :: dH, resid, Qstorage, QnetIn, QnetOut

    !         integer :: mm, kk

    !         real(8), parameter :: localEpsilon = 1.0d-6

    !     !%-----------------------------------------------------------------
    !     !% Preliminaries
    !         if (Npack < 1) return
    !     !%-----------------------------------------------------------------   
    !     !% Aliases
    !         thisJM      => elemP(1:Npack,thisColP)
    !     !%-----------------------------------------------------------------  
        
    !     do mm=1,Npack
    !         !% --- store face flowrate in JB for upstream (1) and downstream (2)
    !         call junction_branch_getface (elemR(:,er_Flowrate),fr_Flowrate,thisJM(mm),ei_Mface_uL,1)
    !         call junction_branch_getface (elemR(:,er_Flowrate),fr_Flowrate,thisJM(mm),ei_Mface_dL,2)

    !             ! if (thisJM(mm) == 634) then
    !             !     print *, ' '
    !             !     print *, '======================================================='
    !             !     print *, 'flowrates JM =', thisJM(mm)
    !             !     print *, elemR(thisJM(mm)+1,er_Flowrate), elemR(thisJM(mm)+2,er_Flowrate), elemR(thisJM(mm),er_FlowrateLateral)
    !             ! end if

    !         !% --- store face energy head in JB for upstream (1) and downstream (2)
    !         call junction_branch_getface (elemSR(:,er_EnergyHead),fr_EnergyHead,thisJM(mm),ei_Mface_uL,1)
    !         call junction_branch_getface (elemSR(:,er_EnergyHead),fr_EnergyHead,thisJM(mm),ei_Mface_dL,2)

    !         ! print *, ' '
    !         ! print *, 'energyHead'
    !         ! print *, elemR(thisJM+1,er_EnergyHead), elemR(thisJM+2,er_EnergyHead)

    !         !% --- compute net flowrate from branches
    !         QnetBranches = junction_main_sumBranches (thisJM(mm),er_Flowrate,elemR)

    !             ! if (thisJM(mm) == 634) then
    !             !     print *, ' '
    !             !     print *, 'Qnetbranches ',QnetBranches
    !             ! end if

    !             ! if (thisJM(mm) == 634) then 
    !             !     print *, ' '
    !             !     print *, 'Qlat         ',elemR(thisJM(mm),er_FlowrateLateral)
    !             ! end if

    !         !% --- compute overflow rate
    !         Qoverflow = junction_main_Qoverflow (thisJM(mm))

    !             ! if (thisJM(mm) == 634) then 
    !             !     print *, ' '
    !             !     print *, 'Qoverflow   ',Qoverflow
    !             ! end if

    !         !% --- compute storage rate of change with head
    !         dQdHstorage = junction_main_dQdHstorage (thisJM(mm),iStep)

    !             ! if (thisJM(mm) == 634) then 
    !             !     print *, ' '
    !             !     print *, 'dQdHstore   ',dQdHstorage
    !             ! end if

    !         !% --- compute overflow rate with change in head
    !         dQdHoverflow = junction_main_dQdHoverflow (thisJM(mm))

    !             ! if (thisJM(mm) == 634) then
    !             !     print *, ' '
    !             !     print *, 'dQdHoverflow ',dQdHoverflow
    !             ! end if

    !         !% --- compute dQdH for each branch -- NEEDS UPDATE to include diagnostic
    !         call junction_branch_dQdH (elemSR(:,esr_JunctionBranch_dQdH),thisJM(mm), elemI(:,ei_Mface_uL), 1)
    !         call junction_branch_dQdH (elemSR(:,esr_JunctionBranch_dQdH),thisJM(mm), elemI(:,ei_Mface_dL), 2)

    !             ! if (thisJM(mm) == 634) then
    !             !     print *, ' '
    !             !     print *, 'dQdH branches'
    !             !     print *, elemSR(thisJM(mm)+1,er_JunctionBranch_dQdH), elemSR(thisJM(mm)+2,esr_JunctionBranch_dQdH)
    !             ! end if

    !         !% --- compute net branch dQdH (product with beta)
    !         dQdHsum =  junction_main_branchdQdHsum (thisJM(mm))

    !             ! if (thisJM(mm) == 634) then
    !             !     print *, ' '
    !             !     print *, 'dQdHsum     ',dQdHsum
    !             ! end if

    !         !% --- compute the junction head change
    !         dH = junction_main_dH &
    !             (thisJM(mm), Qnetbranches, Qoverflow, dQdHstorage, dQdHoverflow, dQdHsum)

    !                 ! if (thisJM(mm) == 634) then
    !                 !     print *, ' '
    !                 !     print *, 'dH before limit = ',dH   
    !                 !     print *, ' '
    !                 ! end if

    !         ! print *, 'max gain ',junction_dH_maxgain (elemR(:,er_EnergyHead), elemR(thisJM(mm),er_Head), thisJM(mm))
    !         ! print *, 'max loss ',junction_dH_maxloss (elemR(:,er_EnergyHead), elemR(thisJM(mm),er_Head), thisJM(mm))
    !         !% --- compute overflow rate
    !             Qoverflow = junction_main_Qoverflow (thisJM(mm))
    !         !% limit dH to prevent flow direction change
    !         if (dH > zeroR) then
    !             dH = min(dH, junction_dH_maxgain (thisJM(mm),dH))
    !         elseif (dH < zeroR) then
    !             dH = max(dH, junction_dH_maxloss (thisJM(mm),dH))
    !         else 
    !             !% -- if dH = zero, no change
    !         end if

    !             ! if (thisJM(mm) == 634) then
    !             !     print *, ' '
    !             !     print *, 'dH after first limit = ',dH   
    !             !     print *, ' '
    !             ! end if

    !         !% limit dH to prevent negative head
    !         if ((elemR(thisJM(mm),er_Head) - elemR(thisJM(mm),er_Zbottom) + dH) < setting%ZeroValue%Depth) then 
    !             ! print *, 'limiters '
    !             ! print *, elemR(thisJM(mm),er_Head),  (elemR(thisJM(mm),er_Zbottom) + 0.99d0*setting%ZeroValue%Depth)
    !             dH = (elemR(thisJM(mm),er_Zbottom) + 0.99d0*setting%ZeroValue%Depth) - elemR(thisJM(mm),er_Head)
    !         end if

    !             ! if (thisJM(mm) == 634) then
    !             !     print *, ' '
    !             !     print *, 'dH after second limit = ',dH   
    !             !     print *, ' '
    !             ! end if

    !         !% --- limit dH dropping based on where overflow shuts off
    !         if (Qoverflow .ne. zeroR) then 
    !             dH = max(dH, junction_dH_overflow_min(thisJM(mm)))
    !         end if

    !             ! if (thisJM(mm) == 634) then
    !             !     print *, ' '
    !             !     print *, 'dH after third limit = ',dH   
    !             ! end if
    !             ! ! ! print *, 'available ', elemR(thisJM(mm),er_Head) - elemR(thisJM(mm),er_Zbottom)
    !             ! ! ! print *, ' '

    !         !% --- update junction head
    !         elemR(thisJM(mm),er_Head) = elemR(thisJM(mm),er_Head) + dH

    !         !% --- update branch Q
    !         elemR(thisJM(mm)+1:thisJM(mm)+max_branch_per_node,er_Flowrate) &
    !             = elemR(thisJM(mm)+1:thisJM(mm)+max_branch_per_node,er_Flowrate)  &
    !             + dH * elemSR(thisJM(mm)+1:thisJM(mm)+max_branch_per_node,esr_JunctionBranch_dQdH)


    !             ! if (thisJM(mm) == 634) then
    !             !     print *, ' '
    !             !     print *, 'flowrates after update', thisJM(mm)
    !             !     !do kk=1,max_branch_per_node
    !             !     do kk=1,2
    !             !         print *, 'Q : ',kk, elemR(thisJM(mm)+kk,er_Flowrate)
    !             !     end do
    !             ! end if

    !         !% --- update overflow
    !         Qoverflow = Qoverflow + dH * dQdHoverflow

    !         !% --- update storage flow rate
    !         select case (elemSI(thisJM(mm),esi_JunctionMain_Type))
    !             case (ImpliedStorage)
    !                 Qstorage = zeroR
    !             case (TabularStorage,FunctionalStorage)
    !                 !% --- compute the new storage volume
    !                 Qstorage = elemSR(thisJM(mm),esr_Storage_Plan_Area) * dH 
    !                 !% --- compute the storage flowrate
    !                 Qstorage = Qstorage / (setting%Solver%crk2(istep) * setting%Time%Hydraulics%Dt )

    !                 ! if (thisJM(mm) == 634) then
    !                 !     print *, ' '
    !                 !     print *, 'Qstorage rate ',Qstorage
    !                 !     print *, ' '
    !                 ! end if 
    !             case default
    !                 print *, 'CODE ERROR: unexpected case default'
    !                 call util_crashpoint(882873)
    !         end select

    !         !% --- compute final residual
    !         resid = junction_conservation_residual (thisJM(mm), Qoverflow, Qstorage) 

    !             ! if (thisJM(mm) == 634) then
    !             !     print *, 'conservation resid ',resid, ';   thisJM ',thisJM(mm)
    !             !     print *, ' '
    !             ! end if

    !         if (abs(resid) > 1.0d-15) then 
    !             !% --- net inflows and net outflows
    !             QnetIn  = junction_branch_Qnet (thisJM(mm),+oneI)
    !             QnetOut = junction_branch_Qnet (thisJM(mm),-oneI) + Qoverflow

    !             ! print *, 'starting conservation resid flows '
    !             ! print *, QnetIn, QnetOut, elemR(thisJM(mm),er_FlowrateLateral)

    !             call junction_fix_conservation(thisJM(mm),resid, Qoverflow, Qstorage, QnetIn, QnetOut)
    !         end if    

    !         !% --- update net Q branches
    !         QnetBranches = junction_main_sumBranches (thisJM(mm),er_Flowrate,elemR)

    !         !% --- update volume based on mass conservation for flux rate at end of
    !         !%     first step considered over the entire step. This is the mass 
    !         !%     flux used for transport in the second step of RK2
    !         if (istep == 1) then
    !              !% --- update volume overflow 
    !             elemR(thisJM(mm),er_VolumeOverflow) = Qoverflow * setting%Time%Hydraulics%Dt 

    !             elemR(thisJM(mm),er_Volume) = elemR(thisJM(mm),er_Volume_N0) &
    !                 + setting%Time%Hydraulics%Dt &
    !                 *(QnetBranches + Qoverflow  + elemR(thisJM(mm),er_FlowrateLateral))
    !         end if

    !             ! if (thisJM(mm) == 634) then
    !             !     print *, ' '
    !             !     print *, 'flowrates after conservation force', thisJM(mm)
    !             !     !do kk=1,max_branch_per_node
    !             !     do kk=1,2
    !             !         print *, 'Q : ',kk, elemR(thisJM(mm)+kk,er_Flowrate)
    !             !     end do

    !             !     print *, ' '
    !             !     print *, 'final conservation resid ',resid
    !             !     print *, ' '

    !             !     print *, elemR(thisJM(mm),er_Volume) - elemR(thisJM(mm),er_Volume_N0) &
    !             !     - setting%Time%Hydraulics%Dt * setting%Solver%crk2(istep) &
    !             !     * ( &
    !             !        QnetBranches &
    !             !     )
    !             ! end if


    !         !% --- push junction JB flowrate values back to face  
    !         call junction_branchface_forceJBvalue (fr_Flowrate, er_Flowrate, ei_Mface_uL, thisJM(mm), 1) 
    !         call junction_branchface_forceJBvalue (fr_Flowrate, er_Flowrate, ei_Mface_dL, thisJM(mm), 2) 

    !         !% --- reset conservative face fluxes
    !         if (istep == 1) then
    !             call junction_branchface_forceJBvalue (fr_Flowrate_Conservative, er_Flowrate, ei_Mface_uL, thisJM(mm), 1) 
    !             call junction_branchface_forceJBvalue (fr_Flowrate_Conservative, er_Flowrate, ei_Mface_dL, thisJM(mm), 2) 
    !         end if

    !             ! print *, ' '
    !             ! print *, 'faces flowrate should be JB values'
    !             ! print *, elemR(5,er_Flowrate), elemR(6,er_Flowrate)
    !             ! print *, faceR(4,fr_Flowrate), faceR(5,fr_Flowrate)

    !             ! if (abs(elemR(5,er_Flowrate) - elemR(6,er_Flowrate)) > localEpsilon) then 
    !             !     print *, 'flowrate mismatch '
    !             !     stop 669874
    !             ! end if 

    !             !% --- note the head values are handled in update_auxiliary_variables_JMJB


    !             ! print *, ' '
    !             ! print *, 'dH after limit = ',dH   
    !             ! print *, ' '

    !             ! if (abs(elemR(thisJM(mm),er_Head)) > 2000.d0) then 
    !             !     print *, 'Strange Head at ',thisJM(mm)
    !             !     print *, elemR(thisJM(mm),er_Head)
    !             !     call util_crashpoint(709874)
    !             ! end if
    !     end do

    ! end subroutine junction_calculation
    !%==========================================================================
!% 
    ! real(8) function junction_dH_maxgain (JMidx,dH)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Computes the maximum head change allowed to prevent reversing
    !     !% flow in all inflow branches
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: JMidx
    !         real(8), intent(in) :: dH

    !         integer :: fidx, ii
    !         real(8) :: dEmax
    !     !%------------------------------------------------------------------
    !     dEmax = zeroR   

    !     do concurrent (ii=1:max_branch_per_node:2)
    !         !% --- select upstream branches with inflows to junction
    !         if ((elemSI(JMidx+ii,esi_JunctionBranch_Exists) == oneI) &
    !             .and.                                                &
    !             ( elemR(JMidx+ii,er_Flowrate) > zeroR)) then 

    !             fidx  = elemI(JMidx+ii,ei_Mface_uL)
    !             dEmax = max(dEmax,faceR(fidx,fr_Head_Adjacent) - elemR(JMidx,er_Head))

    !         end if
    !     end do   

    !     do concurrent (ii=2:max_branch_per_node:2)
    !         !% --- select downstream branches with inflows to junction
    !         if ((elemSI(JMidx+ii,esi_JunctionBranch_Exists) == oneI) &
    !             .and.                                                &
    !             ( elemR(JMidx+ii,er_Flowrate) < zeroR) ) then 

    !             fidx  = elemI(JMidx+ii,ei_Mface_dL)
    !             dEmax = max(dEmax,faceR(fidx,fr_Head_Adjacent) - elemR(JMidx,er_Head))

    !         end if
    !     end do  

    !     junction_dH_maxgain = min(dEmax,dH)
                    
    ! end function junction_dH_maxgain  
!%    
!%==========================================================================
!%==========================================================================
!% 
    ! real(8) pure function junction_dH_maxloss (JMidx,dH)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Computes the largest negative head change allowed to prevent reversing
    !     !% outflow in all outflow branches
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: JMidx
    !         real(8), intent(in) :: dH
    !         integer :: fidx, ii
    !         real(8) :: dEmin
    !     !%------------------------------------------------------------------
    !     !% 
    !     !%------------------------------------------------------------------

    !     dEmin = zeroR

    !     do concurrent (ii=1:max_branch_per_node:2)
    !         !% --- select upstream branches with outflows from junciton
    !         if ((elemSI(JMidx+ii,esi_JunctionBranch_Exists) == oneI) &
    !             .and.                                                &
    !             ( elemR(JMidx+ii,er_Flowrate) < zeroR) ) then 

    !             fidx  = elemI(JMidx+ii,ei_Mface_uL)
    !             dEmin = min(dEmin,faceR(fidx,fr_Head_Adjacent) - elemR(JMidx,er_Head))

    !         end if
    !     end do   

    !     do concurrent (ii=2:max_branch_per_node:2)
    !         if ((elemSI(JMidx+ii,esi_JunctionBranch_Exists) == oneI) &
    !             .and.                                                &
    !             ( elemR(JMidx+ii,er_Flowrate) > zeroR)) then 

    !             fidx  = elemI(JMidx+ii,ei_Mface_dL)
    !             dEmin = min(dEmin,faceR(fidx,fr_Head_Adjacent) - elemR(JMidx,er_Head))

    !         end if
    !     end do  

    !     junction_dH_maxloss = max(dEmin,dH)
                    
    ! end function junction_dH_maxloss
!%    
!%==========================================================================
    !%==========================================================================
!%
    ! subroutine junction_toplevel (whichTM, istep)
    !     print *, 'OBSOLETE'
    !     stop 66908734
        ! !%-----------------------------------------------------------------
        ! !% Description:
        ! !% Controls computation of implicit junction element
        ! !%-----------------------------------------------------------------
        ! !% Declarations
        !     integer, intent(in) :: whichTM, istep
        !     integer             :: thisColP, Npack

        ! !%-----------------------------------------------------------------
        ! !% Preliminaries
        !     select case (whichTM)
        !         case (ALLtm)
        !             thisColP = col_elemP(ep_JM_ALLtm)
        !         case (ETM)
        !             thisColP = col_elemP(ep_JM_ETM)
        !         case (AC)
        !             thisColP = col_elemP(ep_JM_AC)
        !         case default
        !             print *, 'CODE ERROR: time march type unknown for # ', whichTM
        !             print *, 'which has key ',trim(reverseKey(whichTM))
        !             stop 7659
        !     end select

        !     Npack = npack_elemP(thisColP)

        !     if (Npack == 0) return

        !     !% --- coefficients in orifice and weir eqquations
        !     coef1 = twoR * Cbc * sqrt(setting%Constant%gravity * setting%Constant%pi)
        !     coef2 = threehalfR * coef1
        !     coef3 = twothirdR  * sqrt(twoR * setting%Constant%gravity)
        !     coef4 = threehalfR * coef3
        ! !%-----------------------------------------------------------------

        ! ! call junction_calculation (Npack, thisColP, istep)

        !     print *, 'CODE ERROR obsolete routine'
        !     call util_crashpoint(5209873)

        ! ! print *, ' '
        ! ! print *, 'after junction calc ', elemR(4,er_Volume), elemR(4,er_Head)
        ! ! print *, ' '

    ! end subroutine junction_toplevel
    !%
!%==========================================================================
!%    
!%========================================================================== 
!% END OF MODULE
!%==========================================================================
end module junction_elements