module lowlevel_rk2

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use geometry, only : geo_area_from_depth_singular
    use utility, only: util_sign_with_ones
    use utility_output
    use utility_crash, only: util_crashpoint
    ! use utility_unit_testing, only: util_utest_CLprint

    implicit none

    !%-----------------------------------------------------------------------------
    !% Description:
    !% Runge-Kutta method for time-marching (ETM and AC)
    !%
    !% METHOD:
    !%
    !%

    private

    public :: ll_alternate_JB
    public :: ll_continuity_netflowrate_CC
    public :: ll_continuity_netflowrate_JM
    public :: ll_continuity_volume_ETM
    public :: ll_continuity_volume_CCJM_AC_open
    public :: ll_continuity_head_CCJM_AC_surcharged
    public :: ll_continuity_add_gamma_CCJM_AC_open
    public :: ll_continuity_add_source_CCJM_AC_open
    public :: ll_continuity_add_source_CCJM_AC_surcharged
    public :: ll_momentum_Ksource_CC
    public :: ll_momentum_source_CC
    !public :: ll_momentum_lateral_source_CC
    public :: ll_momentum_gammaCM_CC
    public :: ll_momentum_gammaFM_CC
    public :: ll_momentum_solve_CC
    public :: ll_momentum_velocity_CC
    public :: ll_momentum_add_gamma_CC_AC
    public :: ll_momentum_add_source_CC_AC
    public :: ll_minorloss_friction_gamma_CC
    public :: ll_enforce_flapgate_CC
    public :: ll_store_in_temporary
    public :: ll_restore_from_temporary
    public :: ll_extrapolate_values
    public :: ll_interpolate_values
    public :: ll_flowrate_and_velocity_JB
    public :: ll_flowrate_and_velocity_JB_2
    !public :: ll_momentum_solve_JB
    ! public :: ll_CC_slot_computation_ETM
    ! public :: ll_JM_slot_computation_ETM
    public :: ll_get_dynamic_ManningsN

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine ll_continuity_netflowrate_CC (outCol, thisCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Compute net flowrates for channels, conduits and special elements
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: outCol, thisCol, Npack
        real(8), pointer :: fQ(:), eQlat(:)
        integer, pointer :: iup(:), idn(:), thisP(:)
        !%-----------------------------------------------------------------------------
        thisP => elemP(1:Npack,thisCol)
        fQ    => faceR(:,fr_Flowrate)
        eQlat => elemR(:,er_FlowrateLateral)
        iup   => elemI(:,ei_Mface_uL)
        idn   => elemI(:,ei_Mface_dL)
        !%-----------------------------------------------------------------------------

        elemR(thisP,outCol) = fQ(iup(thisP)) - fQ(idn(thisP)) + eQlat(thisP)

        ! print *, 'in ll_continuity_netflowrate_CC'
        ! print *, fQ(iup(1)), fQ(idn(1)), eQlat(1)
        !print *, fQ(iup(ietmp(3))), fQ(idn(ietmp(3))), eQlat(ietmp(3))
        !print *, elemR(ietmp(3),outCol)

    end subroutine ll_continuity_netflowrate_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_continuity_netflowrate_JM (outCol, thisCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% compute net flowrates for junction mains
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: outCol, thisCol, Npack
            real(8), pointer :: branchQ(:), eQlat(:), fQ(:)
            integer, pointer :: thisP(:), isbranch(:), fup(:), fdn(:)
            integer, pointer :: nBarrel(:)
            integer :: ii, jj
        !%------------------------------------------------------------------
        !% Aliases
            thisP    => elemP(1:Npack,thisCol)
            branchQ  => elemR(:,er_Flowrate)
            eQlat    => elemR(:,er_FlowrateLateral)
            isbranch => elemSI(:,esi_JunctionBranch_Exists)
            fQ       => faceR(:,fr_Flowrate)
            fup      => elemI(:,ei_Mface_uL)
            fdn      => elemI(:,ei_Mface_dL)
            nBarrel  => elemI(:,ei_barrels)
        !%------------------------------------------------------------------
        !% note that 1, 3 and 5 are nominal upstream branches and 2, 4, 6 are nominal
        !% downstream branches
        elemR(thisP,outCol) = eQlat(thisP)

        !  print *, 'in ll_continuity_netflowrate_JM'
        ! print *, elemR(iet(1),outCol)
        ! do ii=1,max_branch_per_node,2
        !     print *, fQ(fup(iet(1)+ii)), real(isbranch(iet(1)+ii  ),8)
        !     print *, fQ(fdn(iet(1)+ii+1)), real(isbranch(iet(1)+ii+1),8)
        ! end do

        !% --- testing approach using branch Q
        ! do ii = 1,max_branch_per_node,2
        !     elemR(thisP,outCol) = elemR(thisP,outCol)                 &
        !         + real(isbranch(thisP+ii  ),8) * branchQ(thisP+ii  )  &
        !         - real(isbranch(thisP+ii+1),8) * branchQ(thisP+ii+1)
        ! end do

        ! do ii=1,size(thisP)
        !     print *, thisP(ii),  fQ(fup(thisP(ii)+1  )) ,fQ(fdn(thisP(ii)+2))
        ! end do

        ! do ii=1,size(thisP)
        !     print *, thisP(ii), elemR(thisP(ii),outCol)
        ! end do

        !% --- using face Q up/dn of branch (mass conservative)
        !%     multiply Q by number of barrels of branch
        do ii = 1,max_branch_per_node,2
            elemR(thisP,outCol) = elemR(thisP,outCol) &
                + real(isbranch(thisP+ii  ),8) * fQ(fup(thisP+ii  )) * real(nBarrel(thisP+ii  ),8)  &
                - real(isbranch(thisP+ii+1),8) * fQ(fdn(thisP+ii+1)) * real(nBarrel(thisP+ii+1),8) 
        end do

        ! do ii=1,size(thisP)
        !     print *, thisP(ii), elemR(thisP(ii),outCol)
        ! end do

        ! if (this_image() == 2) then
        !     do ii=1,max_branch_per_node,2
        !         print *, fQ(fup(5428+ii)),   real(isbranch(5428+ii  ),8)
        !         print *, fQ(fdn(5428+ii+1)), real(isbranch(5428+ii+1),8)
        !     end do
        ! end if

        !%-----------------------------------------------------------------
    end subroutine ll_continuity_netflowrate_JM
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_continuity_volume_ETM (outCol, thisCol, Npack, istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Solve for volume from continuity in ETM step
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: outCol, thisCol, Npack, istep
        integer, pointer :: thisP(:)
        real(8), pointer :: Csource(:), VolumeN0(:)
        real(8), pointer :: crk(:), dt
        !%-----------------------------------------------------------------------------
        thisP    => elemP(1:Npack,thisCol)
        VolumeN0 => elemR(:,er_Volume_N0)
        Csource  => elemR(:,er_SourceContinuity)
        crk      => setting%Solver%crk2
        dt       => setting%Time%Hydraulics%Dt
        !%-----------------------------------------------------------------------------

        elemR(thisP,outCol) = VolumeN0(thisP) + crk(istep) * dt * Csource(thisP)

        !% reset the source
        Csource(thisP) = zeroR

        ! if (this_image() == 2) then
        !     write(*,"(A,4f12.4)") ' in ll source ', crk(istep) * dt * Csource(5428)
        ! end if

       !print *, 'in ll_continuity_volume'
       !print *, VolumeN0(1), crk(istep)* dt * Csource(1), elemR(1,outCol)

    end subroutine ll_continuity_volume_ETM
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_continuity_volume_CCJM_AC_open (outCol,  thisCol, Npack, istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Solves continuity for volume in AC method with open channel flow
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: outCol,  thisCol, Npack, istep
        integer, pointer :: thisP(:)
        real(8), pointer :: Csource(:), VolumeM(:), Cgamma(:)
        real(8), pointer :: crk, dtau
        !%-----------------------------------------------------------------------------
        thisP   => elemP(1:Npack,thisCol)
        VolumeM => elemR(:,er_VolumeLastAC) !% last complete AC solve
        Csource => elemR(:,er_SourceContinuity)
        Cgamma  => elemR(:,er_GammaC)
        crk     => setting%Solver%crk2(istep)
        dtau    => setting%ACmethod%dtau
        !%-----------------------------------------------------------------------------

        elemR(thisP,outCol) = &
              ( VolumeM(thisP) + crk * dtau * Csource(thisP) ) &
              / (oneR + crk * dtau * Cgamma(thisP) )

        !% reset the source
        Csource(thisP) = zeroR      

    end subroutine ll_continuity_volume_CCJM_AC_open
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_continuity_head_CCJM_AC_surcharged (outCol, thisCol, Npack, istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Solves continuity for head in AC method with surcharged flow
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: outCol,  thisCol, Npack, istep
        integer, pointer :: thisP(:)
        real(8), pointer :: Csource(:),  Cgamma(:), eHeadM(:)
        real(8), pointer :: crk, dtau
        !%-----------------------------------------------------------------------------
        thisP   => elemP(1:Npack,thisCol)
        eHeadM  => elemR(:,er_HeadLastAC)
        Csource => elemR(:,er_SourceContinuity)
        Cgamma  => elemR(:, er_GammaC)
        crk     => setting%Solver%crk2(istep)
        dtau    => setting%ACmethod%dtau
        !%-----------------------------------------------------------------------------

        elemR(thisP ,outCol) = &
              ( eHeadM(thisP) + crk * dtau * Csource(thisP ) ) &
              / (oneR + crk * dtau * Cgamma(thisP ) )

        !% reset the source
        Csource(thisP) = zeroR      

    end subroutine ll_continuity_head_CCJM_AC_surcharged
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_continuity_add_source_CCJM_AC_open (outCol, thisCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !%
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: outCol,  thisCol, Npack
        !%-----------------------------------------------------------------------------
        !%

        !% HACK -- needs revision for packing
        ! subroutine addto_source_continuity_CCJM_AC_open &
        !     (inoutCol, thisMaskCol)

        ! integer, intent(in) :: inoutCol, thisMaskCol

        ! real(8), pointer :: ell(:), area(:), dHdA(:), head(:)
        ! real(8), pointer :: Qnet(:), volumeN0(:), volumeN1(:)
        ! real(8), pointer :: Fr, a2, a3,
        ! !%-------------------------------------------------
        ! a2 => setting%ACmethod%ImplicitCoef%a2
        ! a3 => setting%ACmethod%ImplicitCoef%a3
        ! Fr => setting%ACmethod%Froude

        ! Qnet => elemR(:,inoutCol) !% used and updated

        ! ell => elemR(:,er_EllDepth)
        ! area => elemR(:,er_Area)
        ! dHdA => elemR(:,er_dHdA)
        ! head => elemR(:,er_Head)

        ! volumeN0 => elemR(:,Volume_N0)
        ! volumeN1 => elemR(:,Volume_N1)

        ! where elemM(:,thisMaskCol)
        !     elemR(:,inoutCol) =  &
        !         ell(:) *( Fr**twoR) / (area(:) * dHdA(:) + head(:))
        !         (Qnet(:) - (a2/dt) * volumeN0(:) - (a3/dt) * volumeN1(:))
        ! endwhere

        print *, "inside ll_continuity_add_source_CCJM_AC_open stub"
        stop 9366

    end subroutine ll_continuity_add_source_CCJM_AC_open
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_continuity_add_source_CCJM_AC_surcharged (outCol, thisCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: outCol,  thisCol, Npack
        !%-----------------------------------------------------------------------------
        !%
                !% HACK needs to be rewritten for packing

        ! subroutine addto_source_continuity_CCJM_AC_surcharged &
        !     (inoutCol, thisMaskCol)

        ! integer, intent(in) :: inoutCol, thisMaskCol

        ! real(8), pointer :: ell(:), area(:), Qnet(:)
        ! real(8), pointer :: Fr
        ! !%-------------------------------------------------
        ! Fr => setting%ACmethod%Froude

        ! Qnet => elemR(:,inoutCol) !% used and updated

        ! ell => elemR(:,er_EllDepth)
        ! area => elemR(:,er_Area)

        ! where elemM(:,thisMaskCol)
        !     elemR(:,inoutCol) =  &
        !         (ell(:) * (Fr**twoR) / area(:)) * Qnet(:)
        ! endwhere

        print *, "inside ll_continuity_add_source_CCMJM_AC_surcharged stub"
        stop 84792

    end subroutine ll_continuity_add_source_CCJM_AC_surcharged
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_continuity_add_gamma_CCJM_AC_open (outCol, thisCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !%
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: outCol,  thisCol, Npack
        !%-----------------------------------------------------------------------------
        !%
        ! subroutine gamma_continuity_CCJM_AC_open &
        !     (er_outCol, thisMaskCol)

        ! !% provides the gamma factor for continuity for the open-channel AC solution

        ! integer, intent(in) :: outCol, thisMaskCol

        ! real(8), pointer :: ell(:), area(:), dHdA(:), head(:)
        ! real(8), pointer :: Fr, a1
        ! !%-------------------------------------------------
        ! a1 => setting%ACmethod%ImplicitCoef%a1
        ! Fr => setting%ACmethod%Froude

        ! ell => elemR(:,er_EllDepth)
        ! area => elemR(:,er_Area)
        ! dHdA => elemR(:,er_dHdA)
        ! head => elemR(:,er_Head)

        ! where elemM(:,thisMaskCol)
        !     elemR(:,outCol) = &
        !         ell(:) * (Fr**twoR) * a1  &
        !         / ( dt * ( area(:) * dHdA(:) + head(:) ) )
        ! endwhere

        print *, "inside ll_continuity_add_gamma_CCJM_AC_open stub"
        stop 29870

    end subroutine ll_continuity_add_gamma_CCJM_AC_open
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_momentum_Ksource_CC (outCol, thisCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% momentum K source terms for different methods for ETM
        !% This is the K term common to AC and ETM momentum advance for
        !% different T00, T10, T20 methods
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: outCol, thisCol, Npack
        real(8), pointer :: fAdn(:), fAup(:), fHdn(:), fHup(:), eHead(:)
        real(8), pointer :: eArea(:), grav
        integer, pointer :: iup(:), idn(:), thisP(:)
        !%-----------------------------------------------------------------------------
        thisP  => elemP(1:Npack,thisCol)
        fAdn   => faceR(:,fr_Area_d)
        fAup   => faceR(:,fr_Area_u)
        fHdn   => faceR(:,fr_Head_d)
        fHup   => faceR(:,fr_Head_u)
        eHead  => elemR(:,er_Head)
        eArea  => elemR(:,er_Area)
        iup    => elemI(:,ei_Mface_uL)
        idn    => elemI(:,ei_Mface_dL)
        grav => setting%constant%gravity
        !%-----------------------------------------------------------------------------

        ! if (setting%Time%Step > 37466) then 
        !     print *, 'in ll_momentum_Ksource_CC '
        !     !print *, setting%Solver%MomentumSourceMethod, trim(reverseKey(setting%Solver%MomentumSourceMethod))
            
        !     print *, idn(1623), iup(1623), faceI(idn(1623),fi_Melem_dL)
        !     print *, fAup(idn(1623)), fAdn(iup(1623))
        !     print *, fHdn(iup(1623)), fHup(idn(1623))
        ! end if

        select case (setting%Solver%MomentumSourceMethod)
        case (T00)
            elemR(thisP,outCol) = grav * ( &
                ( fAup(idn(thisP)) - fAdn(iup(thisP)) ) * eHead(thisP) )

        case (T10)
            elemR(thisP,outCol) = grav * onehalfR *  ( &
                +fAup(idn(thisP)) * fHdn(iup(thisP))   &
                -fAdn(iup(thisP)) * fHup(idn(thisP)) )

            ! print *, ' '
            ! print *, 'in momentum source '
            ! print *, 'first term  ',   fAup(idn(61)) * fHdn(iup(61))
            ! print *, 'second term ',  -fAdn(iup(61)) * fHup(idn(61)) 

        case (T20)
            elemR(thisP,outCol) = grav * onesixthR *  (                       &
                +fAup(idn(thisP)) * ( fHdn(iup(thisP)) + fourR * eHead(thisP) )   &
                -fAdn(iup(thisP)) * ( fHup(idn(thisP)) + fourR * eHead(thisP) ) )

        case (T10s2)   
          
            elemR(thisP,outCol) = onehalfR * grav &
                                *(+( fAup(idn(thisP)) - fAdn(iup(thisP)) ) * eHead(thisP) &
                                  -( fHup(idn(thisP)) - fHdn(iup(thisP)) ) * eArea(thisP) &
                                 )
        case (TA1)
            print *, 'experimental code with Momentum SourceMethod = TA1 should not be used'
            !% note this generally gives a negative source if fAup < fAdn, which is a problem
            stop 59087342
            elemR(thisP,outCol) = grav * ( fAup(idn(thisP)) - fAdn(iup(thisP)) ) &
                                * (                                              &
                                -   onefourthR  * fHdn(iup(thisP))                 &
                                +   threehalfR  * eHead(thisP)                     &
                                -   onefourthR  * fHup(idn(thisP))                 &
                                )

        
            ! print *, 'fA ', fAup(idn(iet(4))) - fAdn(iup(iet(4)))
            ! print *, 'H  ', -onefourthR * fHdn(iup(iet(4))), threehalfR*eHead(iet(4)), -onefourthR*fHup(idn(iet(4)))
            ! print *, thisP
            ! print *, iet(3)
            ! print *, onefourthR, threehalfR
            ! print *, elemR(iet(4),outcol)
            !stop 98734
        case (TA2)
            print *, 'experimental code with Momentum SourceMethod = TA2 should not be used'
            stop 5908734
            !% EXPERIMENTAL, DO NOT USE
            ! elemR(thisP,outCol) = grav   &
            !                     * (                                                 &
            !                     + eArea(thisP) * ( fHdn(iup(thisP)) - eHead(thisP) ) &
            !                     + fAup(idn(thisP)) * eHead(thisP)                      &
            !                     - fAdn(iup(thisP)) * fHdn(iup(thisP))                     &
            !                     )                        
            ! elemR(thisP,outCol) = grav * ( fAup(idn(thisP)) - fAdn(iup(thisP)) ) &
            !                     * (                                              &
            !                     +   onefourthR  * fHdn(iup(thisP))                 &
            !                     +   onehalfR    * eHead(thisP)                     &
            !                     +   onefourthR  * fHup(idn(thisP))                 &
            !                     )
            ! elemR(thisP,outCol) = grav * ( fAup(idn(thisP)) - fAdn(iup(thisP)) ) &
            !                     * (                                              &
            !                     - onehalfR  *  fHdn(iup(thisP))                 &
            !                     + twoR      * eHead(thisP)                     &
            !                     - onehalfR  * fHup(idn(thisP))                 &
            !                     )   
            !  elemR(thisP,outCol) = grav * ( fAup(idn(thisP)) - fAdn(iup(thisP)) ) &
            !                     * (                                              &
            !                     +   onefourthR  * fHdn(iup(thisP))                 &
            !                     -   threehalfR    * eHead(thisP)                     &
            !                     +   onefourthR  * fHup(idn(thisP))                 &
            !                     )                                         
            ! elemR(thisP,outCol) = onehalfR * grav &
            !                     *(+( fAup(idn(thisP)) - fAdn(iup(thisP)) ) * eHead(thisP) &
            !                       -( fHup(idn(thisP)) - fHdn(iup(thisP)) ) * eArea(thisP) &
            !                       +  fAup(idn(thisP)) * fHup(idn(thisP))                  &
            !                       -  fAdn(iup(thisP)) * fHdn(iup(thisP))                  &
            !                     )
        case default
            print *, 'CODE ERROR setting.Solver.MomentumSourceMethod type unknown for # ', setting%Solver%MomentumSourceMethod
            print *, 'which has key ',trim(reverseKey(setting%Solver%MomentumSourceMethod))
            stop 2382
        end select

        ! !print *, ' Ksource ',elemR(780,outCol)
        !  print *, 'in ll_momentum_Ksource_CC'
        !  print *, elemR(iet(4),outCol)

    end subroutine ll_momentum_Ksource_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_momentum_source_CC (outCol, thisCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Common source for momentum on channels and conduits for ETM
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: outCol, thisCol, Npack
        real(8) :: delta
        real(8), pointer :: fQ(:), fUdn(:), fUup(:), fAdn(:), fAup(:)
        real(8), pointer :: fHdn(:), fHup(:), eKsource(:), grav
        integer, pointer :: iup(:), idn(:), thisP(:)
        character(64)    :: subroutine_name = "ll_momentum_source_CC"
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%lowlevel_rk2) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        thisP    => elemP(1:Npack,thisCol)
        fQ       => faceR(:,fr_Flowrate)
        fUdn     => faceR(:,fr_Velocity_d)
        fUup     => faceR(:,fr_Velocity_u)
        fAdn     => faceR(:,fr_Area_d)
        fAup     => faceR(:,fr_Area_u)
        fHdn     => faceR(:,fr_Head_d)
        fHup     => faceR(:,fr_Head_u)
        eKsource => elemR(:,er_Ksource)
        iup      => elemI(:,ei_Mface_uL)
        idn      => elemI(:,ei_Mface_dL)
        grav => setting%constant%gravity
        !%-----------------------------------------------------------------------------

        select case (setting%Solver%MomentumSourceMethod)
        case (T00)
            delta = zeroR
        case (T10)
            delta = onehalfR
        case (T20)
            delta = onesixthR
        case (T10s2)
            delta = onehalfR
        case (TA1)
            delta = zeroR
        case (TA2)
            delta = zeroR
        case default
            print *, 'CODE ERROR setting.Solver.MomentumSourceMethod type unknown for # ', setting%Solver%MomentumSourceMethod
            print *, 'which has key ',trim(reverseKey(setting%Solver%MomentumSourceMethod))
            stop 589794
        end select

        elemR(thisP,outCol) = &
            fQ(iup(thisP)) * fUdn(iup(thisP)) - fQ(idn(thisP)) * fUup(idn(thisP)) &
            + grav * (oneR - delta) &
                *(  &
                    + fAdn(iup(thisP)) * fHdn(iup(thisP))  &
                    - fAup(idn(thisP)) * fHup(idn(thisP))  &
                    ) &
                + eKsource(thisP)

        ! print *, ' '
        ! print *, 'in ll_momentum_source_cc'
        ! print *, 'faces: ',iup(58), idn(58)
        ! print *,  fQ(iup(58)), fUdn(iup(58))
        ! print *, '1st term ',fQ(iup(58)) * fUdn(iup(58))
        ! print *,  -fQ(idn(58)), fUup(idn(58))
        ! print *, '2nd term ',-fQ(idn(58)) * fUup(idn(58))
        ! print *, 'balance of 1-2 ',fQ(iup(58)) * fUdn(iup(58)) - fQ(idn(58)) * fUup(idn(58))
        ! print *, 'pieces ',fQ(iup(58)), fUdn(iup(58))
        ! print *, 'pieces ',fQ(idn(58)), fUup(idn(58))
        ! print *, 'pieces ',fAdn(iup(58)) , fHdn(iup(58))
        ! print *, 'pieces ',fAup(idn(58)) , fHup(idn(58))
        ! print *, '3rd term ',fAdn(iup(58)) * fHdn(iup(58))
        ! print *, '4th term ',-fAup(idn(58)) * fHup(idn(58))
        ! print *, 'balance of 1-2 with coef ',grav * (oneR - delta) * (fAdn(iup(58)) * fHdn(iup(58)) - fAup(idn(58)) * fHup(idn(58)))
        ! print *, 'source ',eKsource(58)
        ! print *, 'output ',elemR(58,outCol)
        ! print *, ' '

        if (setting%Debug%File%lowlevel_rk2) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine ll_momentum_source_CC
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine ll_momentum_lateral_source_CC (inoutCol, thisCol, Npack)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Adding lateral inflow source term to momentum
    !     !% EXPERIMENTAL 20220524 -- DO NOT USE
    !     !%------------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in) :: inoutCol, thisCol, Npack
    !         real(8), pointer :: Qlat(:), Area(:)
    !         integer, pointer :: thisP(:)
    !         character(64)    :: subroutine_name = "ll_momentum_lateral_source_CC"
    !     !%------------------------------------------------------------------
    !     !% Preliminaries:
    !     if (setting%Debug%File%lowlevel_rk2) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%------------------------------------------------------------------
    !     !% Aliases:
    !         thisP    => elemP(1:Npack,thisCol)
    !         Qlat     => elemR(:,er_FlowrateLateral)
    !         Area     => elemR(:,er_Area)
    !     !%------------------------------------------------------------------

    !         print *, 'CODE ERROR: momentum lateral source sould not be used'
    !         stop 559873

    !    ! print *, ' before qlat source ',elemR(780,inoutCol)    

    !     !% HACK the onehalfR should be replaced with a coefficient
    !     elemR(thisP,inoutCol) = elemR(thisP,inoutCol) &
    !         + (Qlat(thisP))**2 / Area(thisP)

    !     !print *, ' after qlat source ',elemR(780,inoutCol)     
    !     !%------------------------------------------------------------------
    !     !% Closing:
    !     if (setting%Debug%File%lowlevel_rk2) &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine ll_momentum_lateral_source_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_momentum_gammaCM_CC (outCol, thisCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% Common Gamma for momentum on channels and conduits 
        !% using the Chezy-Manning roughness approach.
        !% Computes the common part of the Gamma term, which
        !% is the implict friction used in both AC and ETM
        !%-------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: outCol, thisCol, Npack
            real(8), pointer :: velocity(:), mn(:), rh(:),  grav
            integer, pointer :: thisP(:)
            character(64) :: subroutine_name = 'll_momentum_gammaCM_CC'
        !%------------------------------------------------------------------
        !% Aliases
            thisP    => elemP(1:Npack,thisCol)
            velocity => elemR(:,er_velocity)
            rh       => elemR(:,er_HydRadius)
            grav     => setting%constant%gravity
            if (.not. setting%Solver%ManningsN%useDynamicManningsN) then
                mn   => elemR(:,er_ManningsN)
            else
                mn   => elemR(:,er_ManningsN_Dynamic)
            end if
            
        !%---------------------------------------------------------------------

        !% ---- standard Manning's n approach
        elemR(thisP,outCol) =                                       &
                grav * (mn(thisP)**twoR) * abs(velocity(thisP))     &
                /                                                   &
                ( rh(thisP)**fourthirdsR )                         
    
        !    print *, 'in ',trim(subroutine_name)
        !    print *, 'Gamma element values'
        !    print *, elemR(thisP,outCol)
        !    print *, ' '
        !    print *, '============================'
        !    print *, elemR(139,outCol)      
        !    print *, rh(139), mn(139),velocity(139)
        !    print *, elemR(139,er_ManningsN), elemR(139,er_ManningsN_Dynamic)
        !    print *, setting%Solver%ManningsN%useDynamicManningsN

                ! print *, 'in ', trim(subroutine_name)
                ! print *, mn(thisP)

    end subroutine ll_momentum_gammaCM_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_momentum_gammaFM_CC (outCol, thisCol, Npack, FMmethod)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the gamma term in momentum for Force Main roughness
        !% for surcharged pipes.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: outCol, thisCol, Npack, FMmethod
            real(8), pointer   :: grav, velocity(:), AFull(:), Pfull(:)
            real(8), pointer   :: FMcoef(:), DWf(:)
            integer, pointer   :: thisP(:)
            real(8), parameter :: HZfactor = 1.354d0
            real(8), parameter :: HZexpU   = 0.852d0
            real(8), parameter :: HZexpD1  = 1.852d0
            real(8), parameter :: HZexpD2  = 1.1667d0
            character(64) :: subroutine_name = 'll_momentum_gammaFM_CC'
        !%------------------------------------------------------------------
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
            thisP    => elemP(1:Npack,thisCol)
            velocity => elemR(:,er_velocity)
            Afull    => elemR(:,er_FullArea)
            Pfull    => elemR(:,er_FullPerimeter)
            FMcoef   => elemSR(:,esr_Conduit_ForceMain_Coef)
            DWf      => elemSR(:,esr_Conduit_ForceMain_FrictionFactor)
            grav => setting%constant%gravity
        !%------------------------------------------------------------------

        select case (FMmethod)
        case (HazenWilliams)
            elemR(thisP,outCol) =                                            &
                (HZfactor * grav * abs(velocity(thisP))**(HZexpU))           &
                / ( (FMcoef(thisP)**HZexpD1) * ((Afull(thisP) / Pfull(thisP))**HZexpD2) )

                ! print *, 'HazenWilliams in ',trim(subroutine_name)
                ! print *, elemR(thisP,outCol)

        case (DarcyWeisbach)
            elemR(thisP,outCol) = &
                (DWf(thisP) * abs(velocity(thisP)) ) &
                / ( eightR * Afull(thisP)/Pfull(thisP)  )

                ! print *, 'DarcyWeisbach in ',trim(subroutine_name) 
                ! print *, elemR(thisP,outCol)

        case default
            print *, 'CODE ERROR: unexpected case default'
            call util_crashpoint(5592283)
        end select
               
        !%------------------------------------------------------------------
    end subroutine ll_momentum_gammaFM_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_minorloss_friction_gamma_CC (inoutCol, thisCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% Adds the minor loss term for channels and conduits to elemR(:,inoutCol)
        !% note this term is g h_L / L = KU/2L
        !% REVISED 20230202 so that only conduit minor loss is considered
        !% the exit/entrance are now applied in the junction_elements only
        !% Note that for nJ2 and nBC connections the entry/exit losses have
        !% already been added to the Kconduit_MinorLoss for the first and
        !% last elements of the link.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: inoutCol, thisCol, Npack
            real(8), pointer :: velocity(:), oneVec(:)
            real(8), pointer :: Kentry(:), Kexit(:), Kconduit(:), length(:)
            integer, pointer :: thisP(:)
            character(64) :: subroutine_name = 'll_minorloss_friction_CC'
        !%--------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%--------------------------------------------------------------------
        !% Aliases
            thisP    => elemP(1:Npack,thisCol)
            velocity => elemR(:,er_velocity)
            !Kentry   => elemR(:,er_Kentry_MinorLoss)
            !Kexit    => elemR(:,er_Kexit_MinorLoss)
            Kconduit => elemR(:,er_Kconduit_MinorLoss)
            length   => elemR(:,er_Length)
            oneVec   => elemR(:,er_ones)
        !%------------------------------------------------------------------------------

        ! print *, ' '
        ! print *,  elemR(thisP,inoutCol)
        ! print *, ' '

        !% ---- minor loss term (without gravity, which cancels out in derivation)
        elemR(thisP,inoutCol) = elemR(thisP,inoutCol)               &
                + abs(velocity(thisP))                              & 
                * (Kconduit(thisP))                                 &
               /                                                    &
               (twoR * length(thisP)) 

        ! print *, 'in ', trim(subroutine_name)
        ! print *, 'Gamma element updated for minor loss'
        ! print *, elemR(thisP,inoutCol)
        ! !print *, elemR(thisP,inoutCol)
        ! print *, elemR(thisP,inoutCol)
        ! print *, ' '
        ! !print *, velocity(thisP)
        ! !print *, ' '
        ! print *, Kentry(thisP)
        ! print *, ' '
        ! print *, Kexit(thisP)
        ! print *, ' '
        ! print *, Kconduit(thisP)
        ! print *, ' '
        ! print *, length(thisP)
        ! stop 298734


    end subroutine ll_minorloss_friction_gamma_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_momentum_solve_CC (outCol, thisCol, Npack, thisMethod, istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Advance flowrate to n+1/2 for conduits and channels (either ETM or AC)
        !% momentum step that can handle either AC or ETM
        !% but cannot do both at the same time!
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: outCol, thisCol, Npack, thisMethod, istep
        integer, pointer :: thisP(:)
        real(8), pointer :: delt, crk(:)
        real(8), pointer :: volumeLast(:), velocityLast(:), Msource(:), GammaM(:)
        !%-----------------------------------------------------------------------------
        thisP => elemP(1:Npack,thisCOl)
        crk => setting%Solver%crk2
        !%-----------------------------------------------------------------------------

        if (thisMethod == AC) then
            delt         => setting%ACmethod%dtau
            volumeLast   => elemR(:,er_VolumeLastAC)
            velocityLast => elemR(:,er_VelocityLastAC)
        elseif (thisMethod == ETM) then !% real time march
            delt         => setting%Time%Hydraulics%Dt
            volumeLast   => elemR(:,er_Volume_N0)
            velocityLast => elemR(:,er_Velocity_N0)
        else
            print *, 'error, if-else that should not be reached'
            stop 38293
        end if

        Msource => elemR(:,er_SourceMomentum)
        GammaM  => elemR(:,er_GammaM)

        ! print *, ' '
        ! print *, 'in ll_momentum_solve_CC'
        ! print *, ' Msource ', Msource(2189)
        ! print *, ' Gamma   ', GammaM(2189)
        ! print *, 'velocity last ',velocityLast(2189)
        ! print *, ' Vprod   ',volumeLast(61) * velocityLast(2189)
        ! print *, 'crk,delt ', crk(istep),delt
        ! print *, ' '

        elemR(thisP,outCol) =  &
                ( volumeLast(thisP) * velocityLast(thisP) + crk(istep) * delt * Msource(thisP) ) &
                / ( oneR + crk(istep) * delt * GammaM(thisP) )

        ! print *, 'in ll_momentum_solve_CC'
        ! print *, elemR(2189,outCol) 
        ! print *, volumeLast(2189), velocityLast(2189), Msource(2198)
        ! print *, crk(istep), delt, GammaM(2189)  
        ! print *, volumeLast(2189) * velocityLast(2189), crk(istep) * delt * Msource(2189) 
        ! print *, ' '

    end subroutine ll_momentum_solve_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_momentum_velocity_CC (inoutCol, thisCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% velocity for ETM time march
        !%------------------------------------------------------------------
        integer, intent (in) :: inoutCol,  thisCol, Npack
        integer, pointer :: thisP(:)
        real(8), pointer :: momentum(:), volume(:)
        !%------------------------------------------------------------------
        thisP    => elemP(1:Npack,thisCol)
        !% the input integrated momentum is overwritten by the output velocity.
        momentum => elemR(:,inoutCol) !% overwritten
        volume   => elemR(:,er_Volume)
        !%------------------------------------------------------------------
        !% compute velocity

        !print *, ' flowrate ', momentum(780) / elemR(780,er_Length)

        elemR(thisP,inoutCol) = momentum(thisP) / volume(thisP)


        !print *, ' velocity ',elemR(780,inoutCol)
        ! print*
        ! print*, 'in ll_momentum_velocity_CC'
        ! print*, elemR(thisP,inoutCol), 'new velocity'
    end subroutine ll_momentum_velocity_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_enforce_flapgate_CC (inoutCol, thisCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% enforce zero back flow on flapgate elements
        !%------------------------------------------------------------------
        integer, intent (in) :: inoutCol,  thisCol, Npack
        integer, pointer :: thisP(:)
        real(8), pointer :: flow(:)
        !%------------------------------------------------------------------
        thisP    => elemP(1:Npack,thisCol)
        !% the inout column is either velocity, momentum, or flowrate
        flow => elemR(:,inoutCol) !% overwritten
        !%------------------------------------------------------------------

        where ((elemYN(thisP,eYN_hasFlapGate)) .and. (elemR(thisP,inoutCol) < zeroR))
            elemR(thisP,inoutCol) = zeroR
        endwhere

    end subroutine ll_enforce_flapgate_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_momentum_add_gamma_CC_AC (inoutCol, thisCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes the unique part of the Gamma term for AC method
        !% This adds to the existing (common) GammaM
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: inoutCol, thisCol, Npack
        integer, pointer :: thisP(:)
        real(8), pointer :: a1, dt, GammaM(:)
        !%-----------------------------------------------------------------------------
        !%
        thisP => elemP(1:Npack,thisCol)
        a1 => setting%ACmethod%ImplicitCoef%a1
        dt => setting%Time%Hydraulics%Dt
        GammaM => elemR(:,er_GammaM) ! used and updated
        !%-----------------------------------------------------------------------------

        elemR(thisP,inoutCol) = GammaM(thisP) + a1 / dt

    end subroutine ll_momentum_add_gamma_CC_AC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_momentum_add_source_CC_AC (inoutCol, thisCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% this adds to the common source term the unique terms from the AC method
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: inoutCol, thisCol, Npack
        integer, pointer :: thisP(:)
        real(8), pointer :: volumeN0(:), volumeN1(:), Msource(:)
        real(8), pointer :: velocityN0(:), velocityN1(:)
        real(8), pointer :: a2, a3, dt
        !%-----------------------------------------------------------------------------
        thisP => elemP(1:Npack,thisCol)
        dt => setting%Time%Hydraulics%Dt
        a2 => setting%ACmethod%ImplicitCoef%a2
        a3 => setting%ACmethod%ImplicitCoef%a3
        !%-----------------------------------------------------------------------------
        Msource    => elemR(:,inoutCol)
        volumeN0   => elemR(:,er_Volume_N0)
        volumeN1   => elemR(:,er_Volume_N1)
        velocityN0 => elemR(:,er_velocity_N0)
        velocityN1 => elemR(:,er_velocity_N1)

        elemR(thisP,inoutCol) = Msource(thisP) &
                - (+a2 * volumeN0(thisP) * velocityN0(thisP) &
                   +a3 * volumeN1(thisP) * velocityN1(thisP) ) / dt

    end subroutine ll_momentum_add_source_CC_AC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_store_in_temporary (thisCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Copies data to a temporary storage space.
        !% Used in RK2 with ETM/AC to handle time n+1/2 and n+1* communication issues
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisCol, Npack
        integer, pointer :: thisP(:)
        integer :: er_store(3), er_data(3)
        !%-----------------------------------------------------------------------------
        !%
        !% elements in pack
        thisP => elemP(1:Npack,thisCol)

        !% locations for storage
        !% HACK these should be temporary storage that can be re-used outside of
        !% the particular portion of the timeloop. However, this could be tricky
        !% because they need to be not written over!
        er_store= [er_FlowrateStore, er_VolumeStore, er_HeadStore]

        !% locations of data to store
        er_data= [er_Flowrate, er_Volume, er_Head]

        !% copy data to storage
        elemR(thisP,er_store) = elemR(thisP,er_data)

    end subroutine ll_store_in_temporary
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_restore_from_temporary (thisCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% restores data that was temporarily moved in ll_store_in_temporary
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisCol, Npack
        integer, pointer :: thisP(:)
        integer :: er_store(3), er_data(3)
        !%-----------------------------------------------------------------------------
        !% elements in pack
        thisP => elemP(1:Npack,thisCol)

        !% locations of storage
        er_store= [er_FlowrateStore, er_VolumeStore, er_HeadStore]

        !% locations of data to be overwritten
        er_data= [er_Flowrate, er_Volume, er_Head]

        !% overwrite data with storage
        elemR(thisP,er_data) = elemR(thisP,er_store)

    end subroutine ll_restore_from_temporary
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_extrapolate_values (thisCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Extrapolates from time 0 to time 1 using difference  (time 1/2 - time 0)
        !% Used for matching RK2 time levels between AC and ETM methods.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisCol, Npack
        integer, pointer :: thisP(:)
        integer :: eN0(3), eNow(3)
        !%-----------------------------------------------------------------------------
        !% elements in pack
        thisP => elemP(1:Npack,thisCol)

        !% NOTE eN0 and eNow should probably be in settings.
        !% values at time n
        eN0 = [er_Flowrate_N0, er_Volume_N0, er_Head_N0]

        !% values at time n+1/2
        eNow = [er_Flowrate, er_Volume, er_Head]

        !% linear extrapolation to n+1
        elemR(thisP,eNow) = elemR(thisP,eN0) &
                + twoR * ( elemR(thisP,eNow) - elemR(thisP,eN0) )

    end subroutine ll_extrapolate_values
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_interpolate_values (thisCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Interpolates to time n+1/2 from time n=0 and time n+1 data
        !% Used for matching RK2 time levels between AC and ETM methods
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisCol, Npack
        integer, pointer :: thisP(:)
        integer :: eN0(3), eNow(3)
        !%-----------------------------------------------------------------------------
        !% elements in pack
        thisP => elemP(1:Npack,thisCol)

        !% values at time n
        eN0 = [er_Flowrate_N0, er_Volume_N0, er_Head_N0]

        !% values at time n+1
        eNow = [er_Flowrate, er_Volume, er_Head]

        !% linear interpolation to n+1/2
        elemR(thisP,eNow) = onehalfR * ( elemR(thisP,eN0) + elemR(thisP,eNow) )

    end subroutine ll_interpolate_values
!%
!%==========================================================================
!%==========================================================================
!%   
    subroutine ll_flowrate_and_velocity_JB_2 (whichTM, istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Updates the flowrate and velocity on junction branches from face values
        !% obtained in the face interpolation
        !%-----------------------------------------------------------------------------
            integer, intent(in) :: whichTM, istep
            integer, pointer :: thisColP_JM, thisP(:), BranchExists(:), tM, iup(:), idn(:)
            integer, pointer :: Npack
            real(8), pointer :: eHead(:), fHead_u(:), fHead_d(:) !%, fFlowMax(:)
            real(8), pointer :: eFlow(:), fFlow(:), eArea(:), eVelocity(:), eRH(:), vMax
            real(8), pointer :: eVolume(:), eLength(:), dt, grav, epsH, crk(:)
            real(8), pointer :: eRough(:)
            real(8), pointer :: fZbottom(:), fDepth_d(:), fDepth_u(:)
            real(8), pointer :: fVelocity_d(:), fVelocity_u(:)
            integer :: ii, kk, tB
            real(8) :: dHead, gamma, FrFace, headloss
            integer, pointer :: iFaceUp(:), iFaceDn(:)
            integer, pointer :: tFup, tFdn
            logical, pointer :: isZeroDepth(:)
        !%-----------------------------------------------------------------------------
        !%
            BranchExists => elemSI(:,esi_JunctionBranch_Exists)
            eArea        => elemR(:,er_Area)
            eVelocity    => elemR(:,er_Velocity)
            eFlow        => elemR(:,er_Flowrate)
            eVolume      => elemR(:,er_Volume)
            eLength      => elemR(:,er_Length)
            eRH          => elemR(:,er_HydRadius)
            eRough       => elemR(:,er_ManningsN)
    
            fDepth_u     => faceR(:,fr_Depth_u)
            fDepth_d     => faceR(:,fr_Depth_d)
            fVelocity_u  => faceR(:,fr_Velocity_u)
            fVelocity_d  => faceR(:,fr_Velocity_d)
            fDepth_u     => faceR(:,fr_Depth_u)
            fDepth_d     => faceR(:,fr_Depth_d)
            fZbottom     => faceR(:,fr_Zbottom)
            fFlow        => faceR(:,fr_Flowrate)
            iFaceUp      => elemI(:,ei_Mface_uL)
            iFaceDn      => elemI(:,ei_Mface_dL)
    
            eHead        => elemR(:,er_Head)
            fHead_u      => faceR(:,fr_Head_u)
            fHead_d      => faceR(:,fr_Head_d)
    
            isZeroDepth  => elemYN(:,eYN_isZeroDepth)
    
            crk          => setting%Solver%crk2
            vMax         => setting%Limiter%Velocity%Maximum
            dt           => setting%Time%Hydraulics%Dt
            grav         => setting%constant%gravity
            epsH         => setting%Eps%Head
        !%-----------------------------------------------------------------------------
        !%
            select case (whichTM)
            case (ALLtm)
                thisColP_JM            => col_elemP(ep_JM_ALLtm)
            case (ETM)
                thisColP_JM            => col_elemP(ep_JM_ETM)
            case (AC)
                thisColP_JM            => col_elemP(ep_JM_AC)
            case default
                print *, 'CODE ERROR: time march type unknown for # ', whichTM
                print *, 'which has key ',trim(reverseKey(whichTM))
                stop 7659
            end select
    
            !% --- alternate form of junctions skips this
            if (setting%Junction%Method .ne. Explicit2) then 
                return
            end if
        !%-----------------------------------------------------------------------------
            !% --- original form of junctions:
            Npack => npack_elemP(thisColP_JM)
            if (Npack > 0) then
                thisP => elemP(1:Npack,thisColP_JM)
                do ii=1,Npack
                    tM => thisP(ii)  !% JM junction main ID
                    ! handle the upstream branches
                    do kk=1,max_branch_per_node,2
                        tB = tM + kk  !% JB branch ID
                        if (BranchExists(tB)==1) then
                            ! head difference across the branch
                            tFup => iFaceUp(tB)
                            !% use the downstream face, so that near-zero is handled
                            dHead = fHead_d(tFup) - eHead(tB) !% using elem to face

                            if (dHead < zeroR) then
                                !% --- head in JM is larger than upstream branch which implies reversed flow
                                !%     handle different depth conditions
                                if (eHead(tM) .le. (fZbottom(tFup) + setting%ZeroValue%Depth)) then
                                    !% --- JM head below the invert of the upstream conduit, so no flow 
                                    !%     Head on JB is the conduit face head
                                    !eHead(tB)     = fHead_d(tFup) 
                                    eFlow(tB)     = zeroR
                                else 
                                    !% --- JM head above invert of the upstream conduit
                                    !%     expected flow from JM reversed into the upstream branch
                                    !%     handle different velocity conditions
                                    if (fFlow(tFup) > zeroR) then 
                                        !% --- nominal flow against the head gradient is not allowed,
                                        !%     set flux to zero and head to JM head
                                        !eHead(tB)     = eHead(tM)
                                        eFlow(tB)     = zeroR
                                    else
                                        !% --- negative flow
                                        !%     nominal flow and head gradient are into upstream (reversed)
                                        !%     so use the head at JM for the JB head
                                        !eHead(tB)     = eHead(tM)
                                        !%     estimate flowrate driven by head gradient from last value
                                        !%     in JB
                                        gamma = oneR                                 &
                                                +   crk(istep) * dt * grav           &
                                                * abs(eFlow(tB)) * (eRough(tB)**2) &
                                                / (eArea(tB) * eRH(tB)**(fourthirdsR))
    
                                        !% --- note that the dHead is upstream - downstream
                                        eFlow(tB) = (   eFlow(tB)                                               &
                                                    + crk(istep) * dt * grav * eArea(tB) * dHead / eLength(tB) &
                                                    ) / gamma     
                                        
                                        !% --- use the larger magnitude (negative) of the downstream face flowrate of this flowrate
                                        !eFlow(tB)     = min(eFlow(tB),fFlow(tFup))
                                    end if
                                end if

                            elseif (dHead > zeroR) then     
                                !% --- head outside is larger than JM
                                !%     so expected flow into JM from upstream branch
                                !%     handle different depth conditions 
                                if (fHead_d(tFup).le. (fZbottom(tFup) + setting%ZeroValue%Depth)) then
                                    !% --- upstream head lower than invert of conduit at JM connection
                                    !%     so no flow into JM
                                    !eHead(tB)     = fHead_d(tFup)
                                    eFlow(tB)     = zeroR
                                else 
                                    !% --- expected flow into JM from upstream branch
                                    !%     handle waterfall or nonwaterfall conditions
                                    if (eHead(tM) < fZbottom(tFup)) then 
                                         !% --- low depth in JM causing waterfall into JM
                                        if (fDepth_d(tFup) > setting%ZeroValue%Depth) then 
                                            !% --- nonzero depth at face
                                            !% --- compute Froude number on face
                                            FrFace =  fVelocity_d(tFup) / (fDepth_d(tFup) * grav)
                                            if (FrFace .ge. oneR) then 
                                                !% --- supercritical flow
                                                !%     use velocity and flowrate of face
                                                !%     Head is smaller up face head or depth reduced by (HACK) 1/2
                                                !eHead(tB)     = min(fZbottom(tFup) + onehalfR * fDepth_d(tFup), fHead - onehalfR * fDepth_d(tFup))
                                                !% --- adjustment if above HACK makes head too small
                                                !if (eHead(tB) < fZbottom(tFup)) then 
                                                !    eHead(tB) = eHead(tB) + twoR * setting%ZeroValue%Depth
                                                !end if
                                                eFlow(tB)     = fFlow(tFup)
                                            elseif (FrFace < zeroR) then 
                                                !% --- nominally reversed flow
                                                !%     no flow allowed out, head approximated as just above Zbottom
                                                !%     of conduit at connection to JM
                                                !eHead(tB)     = fZbottom(tFup) + twoR *setting%ZeroValue%Depth
                                                eFlow(tB)     = zeroR
                                            else
                                                !% --- subctritical flow ending in waterfall
                                                !% --- HACK -- instead of computing the
                                                !%     critical depth we're just increasing
                                                !%     the velocity to get the headloss
                                                eVelocity(tB) = twoR * fVelocity_d(tFup)
                                                headloss      = (eVelocity(tB)**2 )/ (twoR * grav)
                                                eFlow(tB)     = fFlow(tFup)     
                                                !% --- JB head is diminished by head loss, but at least 
                                                !%     half of the head above the bottom of the conduit entrance
                                                !eHead(tB) = max(fHead - headloss, fZbottom(tFup) + onehalfR*fDepth_d(tFup))                                        
                                            end if
                                        else 
                                            !% --- zero depth at face 
                                            !eHead(tB)     = fHead 
                                            eFlow(tB)     = zeroR
                                        end if

                                    else 
                                        !% --- non waterfall case with hFace > hJM
                                        !%     handle different depth conditions
                                        if (fDepth_d(tFup) > setting%ZeroValue%Depth) then 
                                            !% --- handle different velocity directions

                                            ! print *, 'conditions 4'
                                            ! print *, fVelocity_d(tFup)
                                            ! print *, ' '

                                            if (fVelocity_d(tFup) > zeroR) then 
                                                !% --- velocity is nominal downstream direction
                                                !% --- inflow from face
                                                eFlow(tB)     = fFlow(tFup)
                                                headloss      = (fVelocity_d(tFup)**2 )/ (twoR * grav)
                                                !% --- head with limited head loss

                                                ! print *, 'tfHead ',fHead, headloss, fZbottom(tFup)
                                                ! print *,  fHead - headloss, fZbottom(tFup) + onehalfR*fDepth_d(tFup)
                                                ! print *, ' '

                                                !eHead(tB) = max(fHead_d(tFup)- headloss, fZbottom(tFup) + onehalfR*fDepth_d(tFup))
                                            else 
                                                !% --- outflow against the head gradient is not allowed
                                                !%     use JM head for JB
                                                !eHead(tB)     = eHead(tM) 
                                                eFlow(tB)     = zeroR
                                            end if
                                        else
                                            !% --- zero depth at face 
                                            !eHead(tB)     = eHead(tM)
                                            eFlow(tB)     = zeroR
                                        end if
                                    end if
                                end if
                            else
                                !% --- head difference is identically zero
                                !eHead(tB)     = eHead(tM)
                                eFlow(tB)     = zeroR 
                            end if

                            !% commenting 20230123
                            !% --- no JB driven inflow if JB is zerodepth
                            !%     note that flow across face can still be driven by upstream
                            ! if (isZeroDepth(tB)) then 
                            !     eFlow(tB) = zeroR
                            ! end if

                            !% commenting 20230123
                            !% --- prevent outflow is zero depth JM
                            ! if (isZeroDepth(tM) .and. (eFlow(tB) < zeroR )) then
                            !     eFlow(tB) = zeroR
                            ! end if

                            !% commenting 20230123
                            !% --- prevent outflow driven by JB if head JM <= head JB upstream face
                            ! if ((eHead(tM) .le. fHead_d(tFup)) .and. (eFlow(tB) < zeroR )) then 
                            !     eFlow(tB) = zeroR
                            ! end if

                            !% Fix for velocity blowup due to small areas
                            if (eArea(tB) <= setting%ZeroValue%Area) then
                                eVelocity(tB) = zeroR
                            else
                                eVelocity(tB) = eFlow(tB) / eArea(tB)
                            end if
    
                            if (abs(eVelocity(tB)) > vMax) then
                                eVelocity(tB) = sign( 0.99d0 * vMax, eVelocity(tB) )
                            end if
                        end if
                    end do

                    !% handle the downstream branches
                    do kk=2,max_branch_per_node,2
                        !print *, kk ,'in junction branch'
                        tB = tM + kk
                        if (BranchExists(tB)==1) then
                            !print *, kk, 'in junction branch'
                            tFdn => iFaceDn(tB)
                            !% use the upstream face so that near-zero is handled
                            dHead = eHead(tB) - fHead_u(tFdn) !% using elem to face

                            !% --- handle different head gradient directions
                            if (dHead > zeroR) then
                                !% --- head in JM is larger than upstream branch which implies standard flow
                                !%     handle different depth conditions

                                if (eHead(tM) .le. (fZbottom(tFdn) + setting%ZeroValue%Depth)) then
                                    !% --- JM head below the invert of the downstream conduit, so no flow 
                                    !%     Head on JB is the conduit face head
                                    !eHead(tB)     = fHead 
                                    eFlow(tB)     = zeroR
                                else
                                    !% --- JM head above invert of the upstream conduit
                                    !%     expected flow from JM reversed into the upstream branch
                                    !%     handle different velocity conditions
                                    if (fFlow(tFdn) < zeroR) then 
                                        !% --- nominal flow against the head gradient is not allowed,
                                        !%     set flux to zero and head to JM head
                                        !eHead(tB)     = eHead(tM)
                                        eFlow(tB)     = zeroR
                                    else 
                                        !% --- nominal flow and head gradient are into downstream 
                                        eHead(tB)     = eHead(tM)
                                        !%     estimate flowrate driven by head gradient from last value
                                        !%     in JB
                                        gamma = oneR &
                                                +   crk(istep) * dt * grav          &
                                                * abs(eFlow(tB)) * (eRough(tB)**2) &
                                                / (eArea(tB) * eRH(tB)**(fourthirdsR))
                
                                        eFlow(tB) = (   eFlow(tB)                                                &
                                                    +  crk(istep) * dt * grav * eArea(tB) * dHead / eLength(tB) &
                                                    ) / gamma      
                                      
                                        !% --- use the larger of the downstream face flowrate of this flowrate
                                        !eFlow(tB)     = max(eFlow(tB),fFlow(tFdn))
                                    end if
                                end if

                            elseif (dHead < zeroR) then 
                                !% --- head outside is larger than JM
                                !%     so expected flow into JM from downstream branch (reversed flow)
                                !%     handle different depth conditions 
                                if (fHead_u(tFdn) .le. (fZbottom(tFdn) + setting%ZeroValue%Depth)) then
                                    !% --- downstream head lower than invert of conduit at JM connection
                                    !%     so no flow into JM
                                    !eHead(tB)     = fHead
                                    eFlow(tB)     = zeroR
                                else 
                                    !% --- expected flow into JM from downstream branch
                                    !%     handle waterfall or nonwaterfall conditions
                                    if (eHead(tM) < fZbottom(tFdn)) then 
                                        !% --- low depth in JM causing waterfall into JM
                                        if (fDepth_u(tFdn) > setting%ZeroValue%Depth) then 
                                            !% --- nonzero depth at face
                                            !% --- compute Froude number on face
                                            FrFace =  fVelocity_u(tFdn) / (fDepth_u(tFdn) * grav)
                                            if (FrFace .le. oneR) then 
                                                !% --- upstream supercritical flow
                                                !%     use velocity and flowrate of face
                                                !%     Head is smaller up face head or depth reduced by (HACK) 1/2
                                                !eHead(tB)     = min(fZbottom(tFdn) + onehalfR * fDepth_u(tFdn), fHead - onehalfR * fDepth_u(tFdn))
                                                !% --- adjustment if above HACK makes head too small
                                                !if (eHead(tB) < fZbottom(tFdn)) then 
                                                !    eHead(tB) = eHead(tB) + twoR * setting%ZeroValue%Depth
                                                !end if
                                                eFlow(tB)     = fFlow(tFup)
                                            elseif (FrFace > zeroR) then 
                                                !% --- downstream flow not allowed with reversed head gradient
                                                !%     no flow allowed out, head approximated as just above Zbottom
                                                !%     of conduit at connection to JM
                                                !eHead(tB)     = fZbottom(tFdn) + twoR *setting%ZeroValue%Depth
                                                eFlow(tB)     = zeroR
                                            else
                                                !% --- negative Froude number > -1
                                                !%     subcritical flow upstream (reverse) ending in waterfall
                                                !% --- HACK -- instead of computing the
                                                !%     critical depth we're just increasing
                                                !%     the velocity to get the headloss
                                                eVelocity(tB) = twoR * fVelocity_u(tFdn)
                                                headloss      = (eVelocity(tB)**2 ) / (twoR * grav)
                                                eFlow(tB)     = fFlow(tFdn)     
                                                !% --- JB head is diminished by head loss, but at least 
                                                !%     half of the head above the bottom of the conduit entrance
                                                !eHead(tB) = max(fHead - headloss, fZbottom(tFdn) + onehalfR*fDepth_u(tFdn))                                      
                                            end if
                                        else 
                                            !% --- zero depth at face 
                                            !eHead(tB)     = fHead 
                                            eFlow(tB)     = zeroR
                                        end if
                                    else 
                                        !% --- non waterfall case with hface > hJM implying possible reverse flow
                                        !%     handle different depth conditions
                                        if (fDepth_u(tFdn) > setting%ZeroValue%Depth) then 
                                            !% --- handle different velocity directions
                                            if (fVelocity_u(tFdn) < zeroR) then 
                                                !% --- velocity is nominal upstream direction
                                                !% --- inflow from face
                                                eFlow(tB)     = fFlow(tFdn)
                                                headloss      = (fVelocity_u(tFdn)**2 ) / (twoR * grav)
                                                !% --- head with limited head loss
                                                !eHead(tB) = max(fHead- headloss, fZbottom(tFdn) + onehalfR*fDepth_u(tFdn))
                                            else 
                                                !% --- outflow against the head gradient is not allowed
                                                !%     use JM head for JB
                                                !eHead(tB)     = eHead(tM) 
                                                eFlow(tB)     = zeroR
                                            end if
                                        else
                                            !% --- zero depth at face 
                                            !eHead(tB)     = eHead(tM)
                                            eFlow(tB)     = zeroR
                                        end if
                                    end if
                                end if  
                            else                                      
                                !% --- head difference is identically zero
                                !eHead(tB)     = eHead(tM)
                                eFlow(tB)     = zeroR 
                            end if

                            !% Commenting 20230123
                            !% --- no JB driven inflow if JB is zerodepth
                            !%     note that flow across face can still be driven by upstream
                            ! if (isZeroDepth(tB)) then 
                            !     eFlow(tB) = zeroR
                            ! end if
    
                            ! ! print *, 'tB and flow BBB',tB, eFlow(tB)
    
                            ! !% --- prevent outflow from zero depth JM
                            ! if (isZeroDepth(tM) .and. (eFlow(tB) > zeroR )) then
                            !     eFlow(tB) = zeroR
                            ! end if

    
                            !% Commenting 20230123
                            !% --- prevent JB-driven outflow if head JM <= head JB downstream face
                            ! if ((eHead(tM) .le. fHead_u(tFdn)) .and. (eFlow(tB) > zeroR )) then 
                            !     eFlow(tB) = zeroR
                            ! end if
    
                            ! print *, 'tB and flow DDD',tB, eFlow(tB)
    
                            !% Fix for velocity blowup due to small areas
                            if (eArea(tB) <= setting%ZeroValue%Area) then
                                eVelocity(tB) = zeroR
                            else
                                eVelocity(tB) = eFlow(tB) / eArea(tB)
                            end if
                            ! print *, 'tB and flow DDD',tB, eFlow(tB)
    
                            !print *, 'CCCCCC ',tM, tB, eVelocity(tB)
    
                            if (abs(eVelocity(tB)) > vMax) then
                                eVelocity(tB) = sign( 0.99d0 * vMax, eVelocity(tB) )
                            end if
    
                            !print *, 'DDDDDD ',tM, tB, eVelocity(tB)
                            ! print *, 'tB and flow EEE',tB, eFlow(tB)    
    
                        end if
                    end do
    
                end do
            end if
    
            ! print *, ' '
            ! print *, 'in ll_flowrate_and_velocity at end'
            ! print *, elemR(51,er_Flowrate)
            ! print *, ' '
    
        end subroutine ll_flowrate_and_velocity_JB_2
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_flowrate_and_velocity_JB (whichTM, istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Updates the flowrate and velocity on junction branches from face values
        !% obtained in the face interpolation
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: whichTM, istep
        integer, pointer :: thisColP_JM, thisP(:), BranchExists(:), tM, iup(:), idn(:)
        integer, pointer :: Npack
        real(8), pointer :: eHead(:), fHead_u(:), fHead_d(:) !%, fFlowMax(:)
        real(8), pointer :: eFlow(:), fFlow(:), eArea(:), eVelocity(:), eRH(:), vMax
        real(8), pointer :: eVolume(:), eLength(:), dt, grav, epsH, crk(:)
        real(8), pointer :: eRough(:)
        integer :: ii, kk, tB
        real(8) :: dHead, gamma
        integer, pointer :: iFaceUp(:), iFaceDn(:)
        integer, pointer :: tFup, tFdn
        logical, pointer :: isZeroDepth(:)
        !%-----------------------------------------------------------------------------
        !%
        BranchExists => elemSI(:,esi_JunctionBranch_Exists)
        eArea        => elemR(:,er_Area)
        eVelocity    => elemR(:,er_Velocity)
        eFlow        => elemR(:,er_Flowrate)
        eVolume      => elemR(:,er_Volume)
        eLength      => elemR(:,er_Length)
        eRH          => elemR(:,er_HydRadius)
        eRough       => elemR(:,er_ManningsN)

        fFlow        => faceR(:,fr_Flowrate)
        !fFlowMax     => faceR(:,fr_Flowrate_Max)
        iFaceUp      => elemI(:,ei_Mface_uL)
        iFaceDn      => elemI(:,ei_Mface_dL)

        eHead        => elemR(:,er_Head)
        fHead_u      => faceR(:,fr_Head_u)
        fHead_d      => faceR(:,fr_Head_d)

        isZeroDepth  => elemYN(:,eYN_isZeroDepth)

        crk          => setting%Solver%crk2
        vMax         => setting%Limiter%Velocity%Maximum
        dt           => setting%Time%Hydraulics%Dt
        !rm 20220207brh headC        => setting%Junction%HeadCoef
        grav         => setting%constant%gravity
        epsH         => setting%Eps%Head
        !%-----------------------------------------------------------------------------
        !%
        select case (whichTM)
        case (ALLtm)
            thisColP_JM            => col_elemP(ep_JM_ALLtm)
        case (ETM)
            thisColP_JM            => col_elemP(ep_JM_ETM)
        case (AC)
            thisColP_JM            => col_elemP(ep_JM_AC)
        case default
            print *, 'CODE ERROR: time march type unknown for # ', whichTM
            print *, 'which has key ',trim(reverseKey(whichTM))
            stop 7659
        end select

        !% --- alternate form of junctions skip this
        if (setting%Junction%Method .ne. Explicit1) then
            !% do nothing
            return
        end if

        !% --- original form of junctions:
        Npack => npack_elemP(thisColP_JM)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisColP_JM)
            do ii=1,Npack
                tM => thisP(ii)  !% JM junction main ID
                ! handle the upstream branches
                do kk=1,max_branch_per_node,2
                    tB = tM + kk  !% JB branch ID
                    if (BranchExists(tB)==1) then

                        ! head difference across the branch
                        tFup => iFaceUp(tB)
                        !% use the downstream face, so that near-zero is handled
                        dHead = fHead_d(tFup) - eHead(tB) !% using elem to face

                        gamma = oneR                                 &
                                +   crk(istep) * dt * grav           &
                                  * abs(eFlow(tB)) * (eRough(tB)**2) &
                                  / (eArea(tB) * eRH(tB)**(fourthirdsR))

                        !% --- CHANGED 20230113 to represent the flowrate associated
                        !%     with the head gradient. This is later averaged
                        !%     with the face flowrate to get the final JB and face flowrate

                        !% --- note that the dHead is upstream - downstream
                        eFlow(tB) = (   eFlow(tB)                                                &
                                      + crk(istep) * dt * grav * eArea(tB) * dHead / eLength(tB) &
                                    ) / gamma     
                                    
                        ! eFlow(tB) = (                                                   &
                        !             + crk(istep) * dt * grav * eArea(tB) * dHead / eLength(tB) &
                        !           ) / gamma              
                                    
                                    ! if (tB == 51) then
                                    !     print *, ' '
                                    !     print *, '  here in JB lowlevel'
                                    !     print *, 'flowrate ',eFlow(51)
                                    !     print *, 'area     ' ,eArea(tB)
                                    !     print *, 'dhead    ',dHead
                                    !     print *, 'gamma    ',gamma
                                    !     print *, 'depth    ',elemR(51,er_Depth)
                                    !     print *, 'iszero   ',elemYN(51,eYN_isZeroDepth)
                                    !     print *, ' '
                                    ! end if           
                                    
                        !% --- no JB driven inflow if JB is zerodepth
                        !%     note that flow across face can still be driven by upstream
                        if (isZeroDepth(tB)) then 
                            eFlow(tB) = zeroR
                        end if

                        !% --- prevent outflow is zero depth JM
                        if (isZeroDepth(tM) .and. (eFlow(tB) < zeroR )) then
                            eFlow(tB) = zeroR
                        end if

                        !% --- prevent outflow driven by JB if head JM <= head JB upstream face
                        if ((eHead(tM) .le. fHead_d(tFup)) .and. (eFlow(tB) < zeroR )) then 
                            eFlow(tB) = zeroR
                        end if

                        !% Fix for velocity blowup due to small areas
                        if (eArea(tB) <= setting%ZeroValue%Area) then
                            eVelocity(tB) = zeroR
                        else
                            eVelocity(tB) = eFlow(tB) / eArea(tB)
                        end if

                        !print *, 'AAAAAA ',tM, tB, eVelocity(tB)

                        if (abs(eVelocity(tB)) > vMax) then
                            eVelocity(tB) = sign( 0.99d0 * vMax, eVelocity(tB) )
                        end if

                        !print *, 'BBBBBB ',tM, tB, eVelocity(tB)

                    end if
                end do
                !% handle the downstream branches
                do kk=2,max_branch_per_node,2
                    !print *, kk ,'in junction branch'
                    tB = tM + kk
                    if (BranchExists(tB)==1) then

                        !print *, kk, 'in junction branch'
                        tFdn => iFaceDn(tB)
                        !% use the upstream face so that near-zero is handled
                        dHead = eHead(tB) - fHead_u(tFdn) !% using elem to face
                       
                        gamma = oneR &
                                +   crk(istep) * dt * grav          &
                                  * abs(eFlow(tB)) * (eRough(tB)**2) &
                                  / (eArea(tB) * eRH(tB)**(fourthirdsR))
   
                        eFlow(tB) = (   eFlow(tB)                                                &
                                     +  crk(istep) * dt * grav * eArea(tB) * dHead / eLength(tB) &
                                    ) / gamma      
                                    
                        ! eFlow(tB) = (                                                  &
                        !             +  crk(istep) * dt * grav * eArea(tB) * dHead / eLength(tB) &
                        !            ) / gamma    

                        ! print *, 'tB and flow AAA',tB, eFlow(tB)            
                                   
                        !% --- no JB driven inflow if JB is zerodepth
                        !%     note that flow across face can still be driven by upstream
                        if (isZeroDepth(tB)) then 
                            eFlow(tB) = zeroR
                        end if

                        ! print *, 'tB and flow BBB',tB, eFlow(tB)

                        !% --- prevent outflow from zero depth JM
                        if (isZeroDepth(tM) .and. (eFlow(tB) > zeroR )) then
                            eFlow(tB) = zeroR
                        end if

                        ! print *, 'tB and flow CCC',tB, eFlow(tB)

                        !% --- prevent JB-driven outflow if head JM <= head JB downstream face
                        if ((eHead(tM) .le. fHead_u(tFdn)) .and. (eFlow(tB) > zeroR )) then 
                            eFlow(tB) = zeroR
                        end if

                        ! print *, 'tB and flow DDD',tB, eFlow(tB)

                        !% Fix for velocity blowup due to small areas
                        if (eArea(tB) <= setting%ZeroValue%Area) then
                            eVelocity(tB) = zeroR
                        else
                            eVelocity(tB) = eFlow(tB) / eArea(tB)
                        end if
                        ! print *, 'tB and flow DDD',tB, eFlow(tB)

                        !print *, 'CCCCCC ',tM, tB, eVelocity(tB)

                        if (abs(eVelocity(tB)) > vMax) then
                            eVelocity(tB) = sign( 0.99d0 * vMax, eVelocity(tB) )
                        end if

                        !print *, 'DDDDDD ',tM, tB, eVelocity(tB)
                        ! print *, 'tB and flow EEE',tB, eFlow(tB)    

                    end if
                end do

            end do
        end if

        ! print *, ' '
        ! print *, 'in ll_flowrate_and_velocity at end'
        ! print *, elemR(51,er_Flowrate)
        ! print *, ' '

    end subroutine ll_flowrate_and_velocity_JB
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_get_dynamic_ManningsN (thisP, dpnorm_col) 
        !%------------------------------------------------------------------
        !% Description:
        !% called to get the dynamic ManningsN for a set of points thisP(:)
        !% the dpnorm_col is the location where the normalized pressured 
        !% delta is stored.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisP(:), dpnorm_col
            real(8), pointer    :: dynamic_mn(:), mn(:), dp_norm(:)
            real(8), pointer    :: length(:)
            real(8), pointer    :: alpha, beta, dt, pi
            integer :: ii
            character(64) :: subroutine_name ='ll_get_dynamic_ManningsN'
        !%------------------------------------------------------------------  
        !% Aliases
            pi           => setting%Constant%pi
            dt           => setting%Time%Hydraulics%Dt
            alpha        => setting%Solver%ManningsN%alpha
            beta         => setting%Solver%ManningsN%beta
            mn           => elemR(:,er_ManningsN)
            dynamic_mn   => elemR(:,er_ManningsN_Dynamic)
            dp_norm      => elemR(:,dpnorm_col)
            length       => elemR(:,er_Length)
            
        !%------------------------------------------------------------------
            
       ! dynamic_mn(thisP) =  mn(thisP)
        ! dynamic_mn(thisP) =  mn(thisP) &
        !    +  alpha *  (dt / ((length(thisP))**(onethirdR))) * (exp(dp_norm(thisP)) - oneR )   


        dynamic_mn(thisP) = mn(thisP)                                     &
            * (oneR  + alpha                                              &
                * sin(                                                    &
                      max(zeroR, onehalfR * pi                                &
                             * min( (dp_norm(thisP))/beta, oneR )  &
                          )                                               &
                    )                                                     &
                ) 
            

        print *, 'DYNAMIC MANNINGS NCANNOT BE USED. PRODUCES PROBLEMS AT SMALL DEPTHS.'
        stop 1093874   

        ! do ii=1,size(thisP)
        !     if (dynamic_mn(thisP(ii)) > mn(thisP(ii))) then 
        !         print *, thisP(ii), dynamic_mn(thisP(ii)), dp_norm(thisP(ii))
        !         !stop 3098745
        !     end if
        ! end do

        ! print *, 'in ',trim(subroutine_name)
        ! print *, mn(139), alpha, dt
        ! print *, dp_norm(139)
        ! print *, exp(dp_norm(139))

        ! if (dynamic_mn(50) > 0.05d0) then
        !     print *, dynamic_mn(50)
        ! end if
            
        !% OTHER VERSIONS EXPERIMENTED WITH 20220802
             ! dynamic_mn(thisP) =  mn(thisP) &
           !     +  onehundredR *  (dt / volume**(oneninthR)) * (exp(dp_norm(thisP)) - oneR ) 

           ! dynamic_mn(thisP) =  mn(thisP) &
           !     +  onehundredR *  (dt / ((abs(eHead(thisP) - zBottom(thisP)))**(onethirdR))) * (exp(dp_norm(thisP)) - oneR ) 

    end subroutine ll_get_dynamic_ManningsN
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_alternate_JB (whichTM, istep)
        !%------------------------------------------------------------------
        !% Description:
        !% At the start of this routine, the face values of the JB should
        !% be identical to the elementvalues of the adjacent element. The
        !% values existing on the JB are holdovers from the last term 
        !% the relationship betweent hese face values and the JM are used 
        !% to set the JB values
        !%------------------------------------------------------------------
        !% Declarations:
        !%------------------------------------------------------------------
            integer, intent(in) :: whichTM, istep
            integer, pointer    :: thisColP_JM, Npack
            integer, pointer    :: thisP(:), BranchExists(:)
            integer, pointer    :: iFaceUp(:), iFaceDn(:)
            real(8), pointer    :: eHead(:), fHead_u(:), fHead_d(:)
            real(8), pointer    :: eVelocity(:), fVelocity_u(:), fVelocity_d(:)
            real(8), pointer    :: eFlow(:), fFlow(:)
            real(8), pointer    :: eZbottom(:), fZbottom(:)
            real(8), pointer    :: eDepth(:), fDepth_u(:), fDepth_d(:)
            real(8), pointer    :: eArea(:), fArea_u(:), fArea_d(:)
            real(8), pointer    :: eRH(:), eRough(:), eLength(:), eBottomSlope(:)
            real(8), pointer    :: dt, grav, crk(:)
            integer, pointer    :: tM, tFup, tFdn 

            real(8) :: dHead,  eHeadU,  headloss, FrFace, gamma, fHead, fDepth
            integer :: ii, kk, tB
        !%------------------------------------------------------------------    
        !% Aliases:
            BranchExists => elemSI(:,esi_JunctionBranch_Exists)
            iFaceUp      => elemI(:,ei_Mface_uL)
            iFaceDn      => elemI(:,ei_Mface_dL)
            eHead        => elemR(:,er_Head)
            fHead_u      => faceR(:,fr_Head_u)
            fHead_d      => faceR(:,fr_Head_d)
            eVelocity    => elemR(:,er_Velocity)
            fVelocity_u  => faceR(:,fr_Velocity_u)
            fVelocity_d  => faceR(:,fr_Velocity_d)
            eFlow        => elemR(:,er_Flowrate)
            fFlow        => faceR(:,fr_Flowrate)
            eZbottom     => elemR(:,er_Zbottom)
            fZbottom     => faceR(:,fr_Zbottom)
            eDepth       => elemR(:,er_Depth)
            fDepth_u     => faceR(:,fr_Depth_u)
            fDepth_d     => faceR(:,fr_Depth_d)
            eRH          => elemR(:,er_HydRadius)
            eRough       => elemR(:,er_ManningsN)
            eArea        => elemR(:,er_Area)
            fArea_u      => faceR(:,fr_Area_u)
            fArea_d      => faceR(:,fr_Area_d)
            eLength      => elemR(:,er_Length)
            eBottomSlope => elemR(:,er_BottomSlope)

            crk          => setting%Solver%crk2
            dt           => setting%Time%Hydraulics%Dt
            grav         => setting%Constant%gravity

        !%------------------------------------------------------------------
        !% Preliminaries:
            select case (whichTM)
            case (ALLtm)
                thisColP_JM            => col_elemP(ep_JM_ALLtm)
            case (ETM)
                thisColP_JM            => col_elemP(ep_JM_ETM)
            case (AC)
                thisColP_JM            => col_elemP(ep_JM_AC)
            case default
                print *, 'CODE ERROR: time march type unknown for # ', whichTM
                print *, 'which has key ',trim(reverseKey(whichTM))
                stop 7659
            end select 
        !%------------------------------------------------------------------
    
        Npack => npack_elemP(thisColP_JM)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisColP_JM)
            do ii=1,Npack
                tM => thisP(ii)
                
                ! print *, ' '
                ! print *, 'JB ii ',ii

                !% --- upstream branches
                do kk=1,max_branch_per_node,2
                    tB = tM + kk !% JB branchID 

                    if (BranchExists(tB) == 1) then 
                        tFup => iFaceUp(tB)  !% index of upstream face
                        
                        !% --- provisional head difference
                        dHead  = fHead_d(tFup) - eHead(tM)
                        fHead  = fHead_d(tFup)

                        ! print *, ' '
                        ! print *, fHead_d(tFup), fZbottom(tFup), eZbottom(tB)
                        ! print *, fDepth_d(tFup), fHead_d(tFup)- fZbottom(tFup)
                        ! print *, ' '

                        ! print *, 'provisional dHead ',dHead, fDepth_d(tFup)
                        

                        !% --- Check to see if a head based on uniform flow
                        !%     in the upstream should be used for computations
                        ! if (      (fFlow(tFup)      > zeroR) &
                        !     .and. (dHead            > zeroR) &
                        !     .and. (eBottomSlope(tB) > zeroR) &
                        !     .and. (fHead_d(tFup)    > (fZbottom(tFup) + tenR * setting%ZeroValue%Depth))) then

                        !     !% --- head at outlet assuming uniformly sloping water surface level
                        !     eHeadU  = fHead_d(tFup) - eBottomSlope(tB) * eLength(tB)

                        !     !% --- make sure eHeadU provides a significant depth
                        !     if (eHeadU > (fZbottom(tFup) + onethousandR * setting%ZeroValue%Depth)) then 
                        !         !% --- set values for head and depth based on uniform flow
                        !         fHead  = eHeadU
                        !     else
                        !         !% --- revert to upstream face head and depth
                        !         fHead  = fHead_d(tFup)
                        !     end if
                        ! else 
                        !     !% --- use upstream face head and depth
                        !     fHead  = fHead_d(tFup)

                        ! end if
                        !% --- reset the head difference
                        ! dHead = fHead - eHead(tM)

                        ! print *, 'final dHead ',dHead

                        !% --- handle different head gradient directions
                        if (dHead < zeroR) then
                            !% --- head in JM is larger than upstream branch which implies reversed flow
                            !%     handle different depth conditions
                            if (eHead(tM) .le. (fZbottom(tFup) + setting%ZeroValue%Depth)) then
                                !% --- JM head below the invert of the upstream conduit, so no flow 
                                !%     Head on JB is the conduit face head
                                eHead(tB)     = fHead_d(tFup) 
                                eFlow(tB)     = zeroR
                            else
                                !% --- JM head above invert of the upstream conduit
                                !%     expected flow from JM reversed into the upstream branch
                                !%     handle different velocity conditions
                                if (fFlow(tFup) > zeroR) then 
                                    !% --- nominal flow against the head gradient is not allowed,
                                    !%     set flux to zero and head to JM head
                                    eHead(tB)     = eHead(tM)
                                    eFlow(tB)     = zeroR
                                else 
                                    !% --- negative flow
                                    !%     nominal flow and head gradient are into upstream (reversed)
                                    !%     so use the head at JM for the JB head
                                    eHead(tB)     = eHead(tM)
                                    !%     estimate flowrate driven by head gradient from last value
                                    !%     in JB
                                    gamma = oneR &
                                        +   crk(istep) * dt * grav          &
                                        * abs(eFlow(tB)) * (eRough(tB)**2) &
                                        / (eArea(tB) * eRH(tB)**(fourthirdsR))

                                    !% --- note that the dHead is upstream - downstream, which is
                                    !%     negative so flow will be negative
                                    eFlow(tB) = (   eFlow(tB)                                                &
                                                    + crk(istep) * dt * grav * eArea(tB) * dHead / eLength(tB) &
                                                ) / gamma      
                                    
                                    !% --- use the larger magnitude (negative) of the downstream face flowrate of this flowrate
                                    eFlow(tB)     = min(eFlow(tB),fFlow(tFup))
                                end if
                            end if

                        elseif (dHead > zeroR) then 
                            !% --- head outside is larger than JM
                            !%     so expected flow into JM from upstream branch
                            !%     handle different depth conditions 

                            ! print *, 'conditions 1'
                            ! print *, fHead_d(tFup), fZbottom(tFup)
                            ! print *, ' '

                            if (fHead .le. (fZbottom(tFup) + setting%ZeroValue%Depth)) then
                                !% --- upstream head lower than invert of conduit at JM connection
                                !%     so no flow into JM
                                eHead(tB)     = fHead
                                eFlow(tB)     = zeroR
                            else 
                                ! print *, 'conditions 2'
                                ! print *, eHead(tM), fZbottom(tFup)
                                !print *, 'slope ', elemR(tB,er_BottomSlope)
                                !print *, ' '

                                !% --- expected flow into JM from upstream branch
                                !%     handle waterfall or nonwaterfall conditions
                                if (eHead(tM) < fZbottom(tFup)) then 
                                    !% --- low depth in JM causing waterfall into JM
                                    if (fDepth_d(tFup) > setting%ZeroValue%Depth) then 
                                        !% --- nonzero depth at face
                                        !% --- compute Froude number on face
                                        FrFace =  fVelocity_d(tFup) / (fDepth_d(tFup) * grav)
                                        if (FrFace .ge. oneR) then 
                                            !% --- supercritical flow
                                            !%     use velocity and flowrate of face
                                            !%     Head is smaller up face head or depth reduced by (HACK) 1/2
                                            eHead(tB)     = min(fZbottom(tFup) + onehalfR * fDepth_d(tFup), fHead - onehalfR * fDepth_d(tFup))
                                            !% --- adjustment if above HACK makes head too small
                                            if (eHead(tB) < fZbottom(tFup)) then 
                                                eHead(tB) = eHead(tB) + twoR * setting%ZeroValue%Depth
                                            end if
                                            eFlow(tB)     = fFlow(tFup)
                                        elseif (FrFace < zeroR) then 
                                            !% --- nominally reversed flow
                                            !%     no flow allowed out, head approximated as just above Zbottom
                                            !%     of conduit at connection to JM
                                            eHead(tB)     = fZbottom(tFup) + twoR *setting%ZeroValue%Depth
                                            eFlow(tB)     = zeroR
                                        else
                                            !% --- subctritical flow ending in waterfall
                                            !% --- HACK -- instead of computing the
                                            !%     critical depth we're just increasing
                                            !%     the velocity to get the headloss
                                            eVelocity(tB) = twoR * fVelocity_d(tFup)
                                            headloss      = (eVelocity(tB)**2 )/ (twoR * grav)
                                            eFlow(tB)     = fFlow(tFup)     
                                            !% --- JB head is diminished by head loss, but at least 
                                            !%     half of the head above the bottom of the conduit entrance
                                            eHead(tB) = max(fHead - headloss, fZbottom(tFup) + onehalfR*fDepth_d(tFup))                                        
                                        end if
                                    else 
                                        !% --- zero depth at face 
                                        eHead(tB)     = fHead 
                                        eFlow(tB)     = zeroR
                                    end if
                                else 
                                    ! print *, 'conditions 3'
                                    ! print *, fDepth_d(tFup), setting%ZeroValue%Depth 
                                    ! print *, ' '

                                    !% --- non waterfall case with hFace > hJM
                                    !%     handle different depth conditions
                                    if (fDepth_d(tFup) > setting%ZeroValue%Depth) then 
                                        !% --- handle different velocity directions

                                        ! print *, 'conditions 4'
                                        ! print *, fVelocity_d(tFup)
                                        ! print *, ' '

                                        if (fVelocity_d(tFup) > zeroR) then 
                                            !% --- velocity is nominal downstream direction
                                            !% --- inflow from face
                                            eFlow(tB)     = fFlow(tFup)
                                            headloss      = (fVelocity_d(tFup)**2 )/ (twoR * grav)
                                            !% --- head with limited head loss

                                            ! print *, 'tfHead ',fHead, headloss, fZbottom(tFup)
                                            ! print *,  fHead - headloss, fZbottom(tFup) + onehalfR*fDepth_d(tFup)
                                            ! print *, ' '

                                            eHead(tB) = max(fHead- headloss, fZbottom(tFup) + onehalfR*fDepth_d(tFup))
                                        else 
                                            !% --- outflow against the head gradient is not allowed
                                            !%     use JM head for JB
                                            eHead(tB)     = eHead(tM) 
                                            eFlow(tB)     = zeroR
                                        end if
                                    else
                                        !% --- zero depth at face 
                                        eHead(tB)     = eHead(tM)
                                        eFlow(tB)     = zeroR
                                    end if
                                end if
                            end if    
                        else
                            !% --- head difference is identically zero
                            eHead(tB)     = eHead(tM)
                            eFlow(tB)     = zeroR 
                        end if

                        !% --- set depth 
                        eDepth(tB) = max(eHead(tB) - eZbottom(tB), 0.99d0*setting%ZeroValue%Depth)

                        ! print *, ' '
                        ! print *, 'set depth ',fDepth_d(tFup), eDepth(tB)
                        ! print *, ' '

                        if ((fFlow(tB) > zeroR) .and. (fDepth_d(tFup) > onethousandR * setting%ZeroValue%Depth)) then 
                            eDepth(tB) = min(eDepth(tB), fDepth_d(tFup))
                        end if

                        !% --- set other variables -- note the JB only needs the variables that
                        !%     are used on the face, area, depth, flowrate, head, velocity.
                        eArea(tB)  = geo_area_from_depth_singular(tB,eDepth(tB),setting%ZeroValue%Depth)
                        if ((eFlow(tB) == zeroR) .or. (eArea(tB) .le. setting%ZeroValue%Depth)) then 
                            eVelocity(tB) = zeroR
                        else
                            eVelocity(tB) = sign(max( abs(eFlow(tB)) / eArea(tB), setting%Limiter%Velocity%Maximum),eFlow(tB))
                        end if
                      
                        !% --- push JB values back to faces
                        fArea_d(tFup)     = eArea(tB)
                        fArea_u(tFup)     = eArea(tB)
                        fDepth_d(tFup)    = eDepth(tB)
                        fDepth_u(tFup)    = eDepth(tB)
                        fFlow(tFup)       = eFlow(tB)
                        fHead_d(tFup)     = eHead(tB)
                        fHead_u(tFup)     = eHead(tB)
                        fVelocity_d(tFup) = eVelocity(tB)
                        fVelocity_u(tFup) = eVelocity(tB)
                    else 
                        cycle
                    end if
                end do


                !% --- downstream branches
                do kk=2,max_branch_per_node,2
                    tB = tM + kk !% JB branchID 

                    if (BranchExists(tB) == 1) then 
                        tFdn => iFaceDn(tB)  !% index of upstream face

                        !% --- provisional head difference
                        dHead = eHead(tM) - fHead_u(tFdn) 
                        fHead = fHead_u(tFdn)

                        !% --- Check to see if a head based on uniform flow
                        !%     from downstream to upstream in a reverse slope channel 
                        !%     (an unlikely condition)
                        ! if (      (fFlow(tFdn)     < zeroR) &
                        !     .and. (dHead           < zeroR) &
                        !     .and. (eBottomSlope(tB) < zeroR) &
                        !     .and. (fHead_u(tFdn) > (fZbottom(tFdn))+ tenR*setting%ZeroValue%Depth)) then 
                                
                        !     !% --- expected head at outlet assuming uniformly sloping water surface level
                        !     !%     note this is a reversed flow, since eBottomSlope < 0 this results
                        !     !%     in a smaller head
                        !     eHeadU  = fHead_u(tFdn) + eBottomSlope(tB) * eLength(tB)

                        !     !% --- make sure eHeadU provides a significant depth
                        !     if (eHeadU > (fZbottom(tFdn) + onethousandR * setting%ZeroValue%Depth)) then 
                        !         !% --- set values for head and depth based on uniform flow
                        !         fHead  = eHeadU
                        !     else
                        !         !% --- revert to upstream face head and depth
                        !         fHead  = fHead_u(tFdn)
                        !     end if
                        ! else 
                        !     !% --- use downstream face head and depth
                        !     fHead  = fHead_u(tFdn)
                        ! end if
                        !% --- reset the head difference
                        ! dHead = eHead(tM) - fHead

                        ! print *, 'conditions A'
                        ! print *, dHead, eHead(tM), fHead_u(tFdn)
                        ! print *, ' '

                        !% --- handle different head gradient directions
                        if (dHead > zeroR) then
                            !% --- head in JM is larger than upstream branch which implies standard flow
                            !%     handle different depth conditions

                            ! print *, 'conditions B'
                            ! print *, eHead(tM), fZbottom(tFdn), eHead(tM)-(fZbottom(tFdn) + setting%ZeroValue%Depth)
                            ! print *, ' '

                            if (eHead(tM) .le. (fZbottom(tFdn) + setting%ZeroValue%Depth)) then
                                !% --- JM head below the invert of the downstream conduit, so no flow 
                                !%     Head on JB is the conduit face head
                                eHead(tB)     = fHead 
                                eFlow(tB)     = zeroR
                            else

                                ! print *, 'conditions C'
                                ! print *, fFlow(tFdn)
                                ! print *, ' '

                                !% --- JM head above invert of the upstream conduit
                                !%     expected flow from JM reversed into the upstream branch
                                !%     handle different velocity conditions
                                if (fFlow(tFdn) < zeroR) then 
                                    !% --- nominal flow against the head gradient is not allowed,
                                    !%     set flux to zero and head to JM head
                                    eHead(tB)     = eHead(tM)
                                    eFlow(tB)     = zeroR
                                else 
                                    !% --- nominal flow and head gradient are into downstream 
                                    eHead(tB)     = eHead(tM)
                                    !%     estimate flowrate driven by head gradient from last value
                                    !%     in JB
                                    gamma = oneR &
                                        +   crk(istep) * dt * grav          &
                                        * abs(eFlow(tB)) * (eRough(tB)**2) &
                                        / (eArea(tB) * eRH(tB)**(fourthirdsR))

                                    !% --- note that the dHead is upstream - downstream
                                    eFlow(tB) = (   eFlow(tB)                                                &
                                                    + crk(istep) * dt * grav * eArea(tB) * dHead / eLength(tB) &
                                                ) / gamma      
                                    
                                    !% --- use the larger of the downstream face flowrate of this flowrate
                                    eFlow(tB)     = max(eFlow(tB),fFlow(tFdn))
                                end if
                            end if

                        elseif (dHead < zeroR) then 
                            !% --- head outside is larger than JM
                            !%     so expected flow into JM from downstream branch (reversed flow)
                            !%     handle different depth conditions 
                            if (fHead .le. (fZbottom(tFdn) + setting%ZeroValue%Depth)) then
                                !% --- downstream head lower than invert of conduit at JM connection
                                !%     so no flow into JM
                                eHead(tB)     = fHead
                                eFlow(tB)     = zeroR
                            else 
                                !% --- expected flow into JM from downstream branch
                                !%     handle waterfall or nonwaterfall conditions
                                if (eHead(tM) < fZbottom(tFdn)) then 
                                    !% --- low depth in JM causing waterfall into JM
                                    if (fDepth_u(tFdn) > setting%ZeroValue%Depth) then 
                                        !% --- nonzero depth at face
                                        !% --- compute Froude number on face
                                        FrFace =  fVelocity_u(tFdn) / (fDepth_u(tFdn) * grav)
                                        if (FrFace .le. oneR) then 
                                            !% --- upstream supercritical flow
                                            !%     use velocity and flowrate of face
                                            !%     Head is smaller up face head or depth reduced by (HACK) 1/2
                                            eHead(tB)     = min(fZbottom(tFdn) + onehalfR * fDepth_u(tFdn), fHead - onehalfR * fDepth_u(tFdn))
                                            !% --- adjustment if above HACK makes head too small
                                            if (eHead(tB) < fZbottom(tFdn)) then 
                                                eHead(tB) = eHead(tB) + twoR * setting%ZeroValue%Depth
                                            end if
                                            eFlow(tB)     = fFlow(tFup)
                                        elseif (FrFace > zeroR) then 
                                            !% --- downstream flow not allowed with reversed head gradient
                                            !%     no flow allowed out, head approximated as just above Zbottom
                                            !%     of conduit at connection to JM
                                            eHead(tB)     = fZbottom(tFdn) + twoR *setting%ZeroValue%Depth
                                            eFlow(tB)     = zeroR
                                        else
                                            !% --- negative Froude number > -1
                                            !%     subcritical flow upstream (reverse) ending in waterfall
                                            !% --- HACK -- instead of computing the
                                            !%     critical depth we're just increasing
                                            !%     the velocity to get the headloss
                                            eVelocity(tB) = twoR * fVelocity_u(tFdn)
                                            headloss      = (eVelocity(tB)**2 ) / (twoR * grav)
                                            eFlow(tB)     = fFlow(tFdn)     
                                            !% --- JB head is diminished by head loss, but at least 
                                            !%     half of the head above the bottom of the conduit entrance
                                            eHead(tB) = max(fHead - headloss, fZbottom(tFdn) + onehalfR*fDepth_u(tFdn))                                      
                                        end if
                                    else 
                                        !% --- zero depth at face 
                                        eHead(tB)     = fHead 
                                        eFlow(tB)     = zeroR
                                    end if
                                else 
                                    !% --- non waterfall case with hface > hJM implying possible reverse flow
                                    !%     handle different depth conditions
                                    if (fDepth_u(tFdn) > setting%ZeroValue%Depth) then 
                                        !% --- handle different velocity directions
                                        if (fVelocity_u(tFdn) < zeroR) then 
                                            !% --- velocity is nominal upstream direction
                                            !% --- inflow from face
                                            eFlow(tB)     = fFlow(tFdn)
                                            headloss      = (fVelocity_u(tFdn)**2 ) / (twoR * grav)
                                            !% --- head with limited head loss
                                            eHead(tB) = max(fHead- headloss, fZbottom(tFdn) + onehalfR*fDepth_u(tFdn))
                                        else 
                                            !% --- outflow against the head gradient is not allowed
                                            !%     use JM head for JB
                                            eHead(tB)     = eHead(tM) 
                                            eFlow(tB)     = zeroR
                                        end if
                                    else
                                        !% --- zero depth at face 
                                        eHead(tB)     = eHead(tM)
                                        eFlow(tB)     = zeroR
                                    end if
                                end if
                            end if    
                        else
                            !% --- head difference is identically zero
                            eHead(tB)     = eHead(tM)
                            eFlow(tB)     = zeroR 
                        end if

                        eDepth(tB) = max(eHead(tB) - eZbottom(tB), 0.99d0*setting%ZeroValue%Depth)
                        if ((fFlow(tB) < zeroR) .and. (fDepth_u(tFdn) > onethousandR * setting%ZeroValue%Depth)) then 
                            eDepth(tB) = min(eDepth(tB), fDepth_u(tFdn))
                        end if

                        !% --- set other variables -- note the JB only needs the variables that
                        !%     are used on the face, area, depth, flowrate, head, velocity.
                        eArea(tB)  = geo_area_from_depth_singular(tB,eDepth(tB),setting%ZeroValue%Depth)
                        if ((eFlow(tB) == zeroR) .or. (eArea(tB) .le. setting%ZeroValue%Depth)) then 
                            eVelocity(tB) = zeroR
                        else
                            eVelocity(tB) = sign(max( abs(eFlow(tB)) / eArea(tB), setting%Limiter%Velocity%Maximum),eFlow(tB))
                        end if

                        !% --- push JB values back to faces
                        fArea_d(tFdn)     = eArea(tB)
                        fArea_u(tFdn)     = eArea(tB)
                        fDepth_d(tFdn)    = eDepth(tB)
                        fDepth_u(tFdn)    = eDepth(tB)
                        fFlow(tFdn)       = eFlow(tB)
                        fHead_d(tFdn)     = eHead(tB)
                        fHead_u(tFdn)     = eHead(tB)
                        fVelocity_d(tFdn) = eVelocity(tB)
                        fVelocity_u(tFdn) = eVelocity(tB)
                    else 
                        cycle 
                    end if
                end do
            end do
        end if
    
        !%------------------------------------------------------------------
        !% Closing:

    end subroutine ll_alternate_JB

!%==========================================================================
!%==========================================================================
!%
        !%------------------------------------------------------------------
        !% Description:
        !%
        !%------------------------------------------------------------------
        !% Declarations:
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------
    
    
        !%------------------------------------------------------------------
        !% Closing:

!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module lowlevel_rk2