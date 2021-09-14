module lowlevel_rk2

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys

    implicit none

    real(8), pointer :: grav => setting%constant%gravity

    !%-----------------------------------------------------------------------------
    !% Description:
    !% Runge-Kutta method for time-marching (ETM and AC)
    !%
    !% METHOD:
    !%
    !%

    private

    public :: ll_continuity_netflowrate_CC
    public :: ll_continuity_netflowrate_JM
    public :: ll_continuity_volume_CCJM_ETM
    public :: ll_continuity_volume_CCJM_AC_open
    public :: ll_continuity_head_CCJM_AC_surcharged
    public :: ll_continuity_add_gamma_CCJM_AC_open
    public :: ll_continuity_add_source_CCJM_AC_open
    public :: ll_continuity_add_source_CCJM_AC_surcharged
    public :: ll_momentum_Ksource_CC
    public :: ll_momentum_source_CC
    public :: ll_momentum_gamma_CC
    public :: ll_momentum_solve_CC
    public :: ll_momentum_velocity_CC
    public :: ll_momentum_add_gamma_CC_AC
    public :: ll_momentum_add_source_CC_AC
    public :: ll_store_in_temporary
    public :: ll_restore_from_temporary
    public :: ll_extrapolate_values
    public :: ll_interpolate_values
    public :: ll_junction_branch_flowrate_and_velocity
    public :: ll_slot_computation_ETM

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

    end subroutine ll_continuity_netflowrate_CC
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine ll_continuity_netflowrate_JM (outCol, thisCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% compute net flowrates for junction mains
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: outCol, thisCol, Npack
        real(8), pointer :: fQ(:), eQlat(:)
        integer, pointer :: thisP(:), iup(:), idn(:)
        !%-----------------------------------------------------------------------------
        thisP => elemP(1:Npack,thisCol)
        fQ    => faceR(:,fr_Flowrate)
        eQlat => elemR(:,er_FlowrateLateral)
        iup   => elemI(:,ei_Mface_uL)
        idn   => elemI(:,ei_Mface_dL)
        !%-----------------------------------------------------------------------------

        !% note that 1, 3 and 5 are nominal upstream branches and 2, 4, 6 are nominal
        !% downstream branches

        elemR(thisP,outCol) =  &
            +fQ(iup(thisP+1)) - fQ(idn(thisP+2)) &
            +fQ(iup(thisP+3)) - fQ(idn(thisP+4)) &
            +fQ(iup(thisP+5)) - fQ(idn(thisP+6)) &
            +eQlat(thisP)
    end subroutine ll_continuity_netflowrate_JM
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine ll_continuity_volume_CCJM_ETM (outCol, thisCol, Npack, istep)
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

    end subroutine ll_continuity_volume_CCJM_ETM
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

        ! ell => elemR(:,er_ell)
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

        ! ell => elemR(:,er_ell)
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

        ! ell => elemR(:,er_ell)
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
        integer, pointer :: iup(:), idn(:), thisP(:)
        !%-----------------------------------------------------------------------------
        thisP  => elemP(1:Npack,thisCol)
        fAdn   => faceR(:,fr_Area_d)
        fAup   => faceR(:,fr_Area_u)
        fHdn   => faceR(:,fr_Head_d)
        fHup   => faceR(:,fr_Head_u)
        eHead  => elemR(:,er_Head)
        iup    => elemI(:,ei_Mface_uL)
        idn    => elemI(:,ei_Mface_dL)
        !%-----------------------------------------------------------------------------

        select case (setting%Solver%MomentumSourceMethod)
            case (T00)
                elemR(thisP,outCol) = grav * ( &
                    ( fAup(idn(thisP)) - fAdn(iup(thisP)) ) * eHead(thisP) )
            case (T10)
                elemR(thisP,outCol) = grav * onehalfR *  ( &
                    +fAup(idn(thisP)) * fHdn(iup(thisP))   &
                    -fAdn(iup(thisP)) * fHup(idn(thisP)) )
            case (T20)
                elemR(thisP,outCol) = grav * onesixthR *  (                       &
                    +fAup(idn(thisP)) * ( fHdn(iup(thisP)) + fourR * eHead(:) )   &
                    -fAdn(iup(thisP)) * ( fHup(idn(thisP)) + fourR * eHead(thisP) ) )
            case default
                print *, 'error, case default that should not be reached'
                stop 2382
            !% Error
        end select
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
        real(8), pointer :: fHdn(:), fHup(:), eKsource(:)
        integer, pointer :: iup(:), idn(:), thisP(:)
        character(64)    :: subroutine_name = "ll_momentum_source_CC"
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%lowlevel_rk2) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

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
        !%-----------------------------------------------------------------------------

        select case (setting%Solver%MomentumSourceMethod)
            case (T00)
                delta = zeroR
            case (T10)
                delta = onehalfR
            case (T20)
                delta = onesixthR
            case default
                stop "in " // subroutine_name
            !% Error
        end select

        elemR(thisP,outCol) = &
            fQ(iup(thisP)) * fUdn(iup(thisP)) - fQ(idn(thisP)) * fUup(idn(thisP)) &
            + grav * (oneR - delta) &
                *(  &
                    + fAdn(iup(thisP)) * fHdn(iup(thisP))  &
                    - fAup(idn(thisP)) * fHup(idn(thisP))  &
                    ) &
                + eKsource(thisP)

        if (setting%Debug%File%lowlevel_rk2) &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine ll_momentum_source_CC
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine ll_momentum_gamma_CC (outCol, thisCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Common Gamma for momentum on channels and conduits for  ETM
        !% Computes the common part of the Gamma term, which
        !% is the implict friction used in both AC and ETM
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: outCol, thisCol, Npack
        real(8), pointer :: velocity(:), mn(:), rh(:), oneVec(:)
        integer, pointer :: thisP(:)
        !%------------------------------------------------------------------------------
        thisP    => elemP(1:Npack,thisCol)
        velocity => elemR(:,er_velocity)
        mn       => elemR(:,er_Roughness)
        rh       => elemR(:,er_HydRadius)
        oneVec   => elemR(:,er_ones)
        !%------------------------------------------------------------------------------

        elemR(thisP,outCol) = &
                sign(oneVec(thisP), velocity(thisP)) &
                * grav * (mn(thisP)**twoR) * velocity(thisP)  &
                / &
                ( rh(thisP)**fourthirdsR )

    end subroutine ll_momentum_gamma_CC
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

        elemR(thisP,outCol) =  &
                ( volumeLast(thisP) * velocityLast(thisP) + crk(istep) * delt * Msource(thisP) ) &
                / ( oneR + crk(istep) * delt * GammaM(thisP) )

    end subroutine ll_momentum_solve_CC
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine ll_momentum_velocity_CC (inoutCol, thisCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% velocity for ETM time march
        !%-----------------------------------------------------------------------------
        integer, intent (in) :: inoutCol,  thisCol, Npack
        integer, pointer :: thisP(:)
        real(8), pointer :: momentum(:), volume(:)
        !%-----------------------------------------------------------------------------
        thisP    => elemP(1:Npack,thisCol)
        !% the input integrated momentum is overwritten by the output velocity.
        momentum => elemR(:,inoutCol) !% overwritten
        volume   => elemR(:,er_Volume)
        !%-----------------------------------------------------------------------------
        !% compute velocity

        elemR(thisP,inoutCol) = momentum(thisP) / volume(thisP)

        !% zero out velocities in those cells with near zero volumes
        where (elemYN(thisP,eYN_isNearZeroVolume))
            elemR(thisP,inoutCol) = zeroR
        endwhere

    end subroutine ll_momentum_velocity_CC
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
    subroutine ll_junction_branch_flowrate_and_velocity (whichTM)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Updates the flowrate and velocity on junction branches from face values
        !% obtained in the face interpolation
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: whichTM
        integer, pointer :: thisColP_JM, thisP(:), BranchExists(:), tM, iup(:), idn(:)
        integer, pointer :: Npack
        real(8), pointer :: eHead(:), fHead_u(:), fHead_d(:) ! BRHbugfix 20210811
        real(8), pointer :: eFlow(:), fFlow(:), eArea(:), eVelocity(:), vMax ! BRHbugfix 20210829
        real(8), pointer :: eVolume(:), dt  !BRHbugfix 20210829
        logical, pointer :: isAdhocFlowrate(:)
        integer :: ii, kk, tB
        !% BRHbugfix 20210812 start
        real(8) :: dHead
        integer, pointer :: iFaceUp(:), iFaceDn(:) !, iElemUp(:), iElemDn(:)
        integer, pointer :: tFup, tFdn !, tEup, tEdn
        !% BRHbugfix 20210812 end
        !%-----------------------------------------------------------------------------
        !%
        BranchExists => elemSI(:,eSI_JunctionBranch_Exists)
        eArea        => elemR(:,er_Area)
        eVelocity    => elemR(:,er_Velocity)
        eFlow        => elemR(:,er_Flowrate)
        eVolume      => elemR(:,er_Volume)       ! BRHbugfix 20210829

        fFlow        => faceR(:,fr_Flowrate)
        iFaceUp      => elemI(:,ei_Mface_uL) !% BRHbugfix 20210811
        iFaceDn      => elemI(:,ei_Mface_dL)!% BRHbugfix 20210811
        !iElemUp      => faceI(:,fi_Melem_uL) !% BRHbugfix 20210811
        !iElemDn       => faceI(:,fi_Melem_dL)!% BRHbugfix 20210811
        !% BRHbugfix 20210811 start
        eHead        => elemR(:,er_Head)
        fHead_u      => faceR(:,fr_Head_u)
        fHead_d      => faceR(:,fr_Head_d)
        vMax         => setting%Limiter%Velocity%Maximum
        isAdhocFlowrate => elemYN(:,eYN_IsAdhocFlowrate)
        !% BRHbugfix 20210811 end
        dt           => setting%Time%Hydraulics%Dt  ! BRHbugfix 20210829
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
            print *, 'error, case default should never be reached.'
            stop 7659
        end select

        Npack => npack_elemP(thisColP_JM)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisColP_JM)
            ! BRHbugfix 20210811 start
            do ii=1,Npack
                tM => thisP(ii)
                ! handle the upstream branches
                do kk=1,max_branch_per_node,2
                    tB = tM + kk
                    if (BranchExists(tB)==1) then
                        ! head difference across the branch
                        tFup => iFaceUp(tB)
                        !tEup => iElemUp(tFup)
                        !dHead = eHead(tEup) - eHead(tB) !% using elem to elem
                        dHead = fHead_u(tFup) - eHead(tB) !% using elem to face
                        if (dHead >= zeroR) then
                            ! downstream flow in an upstream branch use upstream values
                            !BRH bugfix 20210829eFlow(tB) = fFlow(tFup) !% using face
                            !eFlow(tB) = eFlow(tEup)!%  using elem
                            eFlow(tB) = eArea(tB) * sqrt(twoR * setting%Constant%gravity * dHead) !BRH bugfix 20210829
                        else
                            ! upstream flow in an upstream branch
                            eFlow(tB) = - eArea(tB) * sqrt(twoR * setting%Constant%gravity * (-dHead))
                            ! if outflow, limit negative flowrate by 1/3 main volume
                            eFlow(tB) = max(eFlow(tB), -eVolume(tM)/(threeR * dt) ) !BRHbugfix 20210829
                        end if

                        !% HACK: Fix for velocity blowup
                        if (eArea(tB) <= setting%ZeroValue%Area) then
                            eVelocity(tB) = zeroR
                        else
                            eVelocity(tB) = eFlow(tB) / eArea(tB)
                        end if

                        if (abs(eVelocity(tB)) > vMax) then
                            eVelocity(tB) = sign( 0.99 * vMax, eVelocity(tB) )
                            isAdhocFlowrate(tB) = .true.
                        end if

                    end if
                end do
                !% handle the downstream branches
                do kk=2,max_branch_per_node,2
                    tB = tM + kk
                    if (BranchExists(tB)==1) then
                        tFdn => iFaceDn(tB)
                        !tEdn => iElemDn(tFdn)
                        !dHead =  eHead(tB) - eHead(tEdn) !% using elem to elem
                        dHead = eHead(tB) - fHead_d(tFdn) !% using elem to face
                        if (dHead < zeroR) then
                            ! upstream flow in a downstream branch use downstream values
                            !BRH bugfix 20210829 eFlow(tB) = fFlow(tFdn) !% using face
                            !eFlow(tB) = eFlow(tEdn) !%  using elem
                            eFlow(tB) =  - eArea(tB) * sqrt(twoR * setting%Constant%gravity * (-dHead) ) !BRH bugfix 20210829
                        else
                            ! downstream flow in an downstream branch
                            eFlow(tB) =  + eArea(tB) * sqrt(twoR * setting%Constant%gravity * dHead )
                            ! if outflow, limit flowrate by 1/3 main volume
                            eFlow(tB) = min(eFlow(tB), eVolume(tM)/(threeR * dt) )  !BRHbugfix 20210829
                        end if

                        !% HACK: Fix for velocity blowup
                        if (eArea(tB) <= setting%ZeroValue%Area) then
                            eVelocity(tB) = zeroR
                        else
                            eVelocity(tB) = eFlow(tB) / eArea(tB)
                        end if

                        if (abs(eVelocity(tB)) > vMax) then
                            eVelocity(tB) = sign( 0.99 * vMax, eVelocity(tB) )
                            isAdhocFlowrate(tB) = .true.
                        end if

                    end if
                end do

            end do

            ! remove old stuff
            !     ! handle the upstream branches
            !     do kk=1,max_branch_per_node,2
            !         tB = tM + kk
            !         if (BranchExists(tB)==1) then
            !             eFlow(tB) = fFlow(iup(tB))
            !             eVelocity(tB) = eFlow(tB) / eArea(tB)
            !         end if
            !     end do
            !     !% handle the downstream branches
            !     do kk=2,max_branch_per_node,2
            !         tB = tM + kk
            !         if (BranchExists(tB)==1) then

            !             eFlow(tB) = fFlow(idn(tB))
            !             eVelocity(tB) = eFlow(tB) / eArea(tB)
            !         end if
            !     end do
            ! end do
            ! BRHbugfix 20210811 end
        end if

    end subroutine ll_junction_branch_flowrate_and_velocity
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine ll_slot_computation_ETM (thisCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Compute preissmann slot for conduits in ETM methods
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisCol, Npack
        integer, pointer    :: thisP(:), SlotMethod
        real(8), pointer    :: SlotWidth(:), SlotVolume(:), SlotDepth(:), SlotArea(:)
        real(8), pointer    :: volume(:), fullvolume(:), fullarea(:), ell(:), length(:)
        real(8), pointer    :: CelerityFactor

        character(64) :: subroutine_name = 'll_slot_computation_ETM'
        !%-----------------------------------------------------------------------------
        thisP => elemP(1:Npack,thisCol)
        volume     => elemR(:,er_Volume)
        fullvolume => elemR(:,er_FullVolume)
        fullarea   => elemR(:,er_FullArea)
        ell        => elemR(:,er_ell)
        length     => elemR(:,er_Length)
        SlotWidth  => elemSR(:,eSr_conduit_SlotWidth)
        SlotVolume => elemSR(:,eSr_conduit_SlotVolume)
        SlotDepth  => elemSR(:,eSr_conduit_SlotDepth)
        SlotArea   => elemSR(:,eSr_conduit_SlotArea)

        SlotMethod     => setting%PreissmannSlot%PreissmannSlotMethod
        CelerityFactor => setting%PreissmannSlot%CelerityFactor

        select case (SlotMethod)

            case (VariableSlot)

                SlotVolume(thisP) = max(volume(thisP) - fullvolume(thisP), zeroR)
                SlotWidth(thisP)  = fullarea(thisP) / (CelerityFactor * ell(thisP))
                SlotArea(thisP)   = SlotVolume(thisP) / length(thisP)
                SlotDepth(thisP)  = SlotArea(thisP) / SlotWidth(thisP)
                print*, thisP, 'thisP'
                print*, volume(thisP), 'volume(thisP)'
                print*, fullvolume(thisP), 'fullvolume(thisP)'
                print*, SlotVolume(thisP), 'SlotVolume(thisP)'
                print*, SlotWidth(thisP), 'SlotWidth(thisP)'
                print*, SlotArea(thisP), 'SlotArea(thisP)'
                print*, SlotDepth(thisP), 'SlotDepth(thisP)'

            case (StaticSlot)

                print*, 'In ', subroutine_name
                print*, 'Static Preissmann Slot is not handeled yet'
                stop "in " // subroutine_name

            case default
                !% should not reach this stage
                print*, 'In ', subroutine_name
                print*, 'error: unexpected Preissmann Slot Method, ', SlotMethod
                stop "in " // subroutine_name

        end select

    end subroutine ll_slot_computation_ETM
    !%
    !%==========================================================================
    !%==========================================================================
    !%
        !%-----------------------------------------------------------------------------
        !% Description:
        !%
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%

    !%==========================================================================
    !% END OF MODULE
    !%+=========================================================================
    !%==========================================================================
    !% PRIVATE
    !%==========================================================================

end module lowlevel_rk2