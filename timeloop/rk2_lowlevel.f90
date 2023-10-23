module rk2_lowlevel
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Provides the lower-level procedures for the Runge-Kutta 2 time march
    !%
    !%==========================================================================

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use geometry, only : geo_area_from_depth_singular
    use utility, only: util_sign_with_ones
    use utility_output
    use utility_crash, only: util_crashpoint

    implicit none
    private

    public :: ll_continuity_netflowrate_CC
    public :: ll_continuity_volume_CC
    public :: ll_momentum_Ksource_CC
    public :: ll_momentum_source_CC
    public :: ll_momentum_gammaCM_CC
    public :: ll_momentum_gammaFM_CC
    public :: ll_minorloss_friction_gamma_CC
    public :: ll_momentum_solve_CC
    public :: ll_momentum_velocity_CC
    public :: ll_enforce_flapgate_CC
    public :: ll_enforce_zerodepth_velocity

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine ll_continuity_netflowrate_CC (outCol, thisCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% Compute net flowrates for channels, conduits and special elements
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: outCol, thisCol, Npack
            real(8), pointer    :: fQ(:), eQlat(:)
            integer, pointer    :: iup(:), idn(:), thisP(:)
        !%------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
            thisP => elemP(1:Npack,thisCol)
            fQ    => faceR(:,fr_Flowrate)
            eQlat => elemR(:,er_FlowrateLateral)
            iup   => elemI(:,ei_Mface_uL)
            idn   => elemI(:,ei_Mface_dL)
        !%------------------------------------------------------------------

        elemR(thisP,outCol) = fQ(iup(thisP)) - fQ(idn(thisP)) + eQlat(thisP)

    end subroutine ll_continuity_netflowrate_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_continuity_volume_CC (outCol, thisCol, Npack, istep)
        !%------------------------------------------------------------------
        !% Description:
        !% Solve for volume from continuity for CC elements
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: outCol, thisCol, Npack, istep
            integer, pointer    :: thisP(:)
            real(8), pointer    :: Csource(:), VolumeN0(:)
            real(8), pointer    :: crk(:), dt
        !%------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
            thisP    => elemP(1:Npack,thisCol)
            VolumeN0 => elemR(:,er_Volume_N0)
            Csource  => elemR(:,er_SourceContinuity)
            crk      => setting%Solver%crk2
            dt       => setting%Time%Hydraulics%Dt
        !%------------------------------------------------------------------

        elemR(thisP,outCol) = VolumeN0(thisP) + crk(istep) * dt * Csource(thisP)

        !% --- reset the source
        Csource(thisP) = zeroR

    end subroutine ll_continuity_volume_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_momentum_Ksource_CC (outCol, thisCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% momentum K source terms for different methods for ETM
        !% This is the K term common to AC and ETM momentum advance for
        !% different T00, T10, T20 methods
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: outCol, thisCol, Npack
            real(8), pointer    :: fAdn(:), fAup(:), fHdn(:), fHup(:), eHead(:)
            real(8), pointer    :: eArea(:), grav
            integer, pointer    :: iup(:), idn(:), thisP(:)
        !%------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
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
        !%------------------------------------------------------------------

        select case (setting%Solver%MomentumSourceMethod)
            case (T00)
                elemR(thisP,outCol) = grav * (                               &
                    ( fAup(idn(thisP)) - fAdn(iup(thisP)) ) * eHead(thisP) )

            case (T10)
                elemR(thisP,outCol) = grav * onehalfR *  ( &
                    +fAup(idn(thisP)) * fHdn(iup(thisP))   &
                    -fAdn(iup(thisP)) * fHup(idn(thisP)) )

            case (T20)
                elemR(thisP,outCol) = grav * onesixthR *  (                           &
                    +fAup(idn(thisP)) * ( fHdn(iup(thisP)) + fourR * eHead(thisP) )   &
                    -fAdn(iup(thisP)) * ( fHup(idn(thisP)) + fourR * eHead(thisP) ) )

            case (T2L0)
                elemP(thisP,outCol) = grav * onefourthR * (   &
                    +fAup(idn(thisP)) * ( fHdn(iup(thisP)) + twoR * eHead(thisP) )   &
                    -fAdn(iup(thisP)) * ( fHup(idn(thisP)) + twoR * eHead(thisP) ) )       

            ! case (T10s2)   
            !     elemR(thisP,outCol) = onehalfR * grav &
            !                         *(+( fAup(idn(thisP)) - fAdn(iup(thisP)) ) * eHead(thisP) &
            !                         -( fHup(idn(thisP)) - fHdn(iup(thisP)) ) * eArea(thisP)   &
            !                         )
            ! case (TA1)
            !     print *, 'experimental code with Momentum SourceMethod = TA1 should not be used'
            !     !% --- note this generally gives a negative source if fAup < fAdn, which is a problem
            !     call util_crashpoint(1108744)
            !     elemR(thisP,outCol) = grav * ( fAup(idn(thisP)) - fAdn(iup(thisP)) )   &
            !                         * (                                                &
            !                         -   onefourthR  * fHdn(iup(thisP))                 &
            !                         +   threehalfR  * eHead(thisP)                     &
            !                         -   onefourthR  * fHup(idn(thisP))                 &
            !                         )
            ! case (TA2)
            !     print *, 'experimental code with Momentum SourceMethod = TA2 should not be used'
            !     call util_crashpoint(7119873)
            !     !% EXPERIMENTAL, DO NOT USE
            !     ! elemR(thisP,outCol) = grav   &
            !     !                     * (                                                 &
            !     !                     + eArea(thisP) * ( fHdn(iup(thisP)) - eHead(thisP) ) &
            !     !                     + fAup(idn(thisP)) * eHead(thisP)                      &
            !     !                     - fAdn(iup(thisP)) * fHdn(iup(thisP))                     &
            !     !                     )                        
            !     ! elemR(thisP,outCol) = grav * ( fAup(idn(thisP)) - fAdn(iup(thisP)) ) &
            !     !                     * (                                              &
            !     !                     +   onefourthR  * fHdn(iup(thisP))                 &
            !     !                     +   onehalfR    * eHead(thisP)                     &
            !     !                     +   onefourthR  * fHup(idn(thisP))                 &
            !     !                     )
            !     ! elemR(thisP,outCol) = grav * ( fAup(idn(thisP)) - fAdn(iup(thisP)) ) &
            !     !                     * (                                              &
            !     !                     - onehalfR  *  fHdn(iup(thisP))                 &
            !     !                     + twoR      * eHead(thisP)                     &
            !     !                     - onehalfR  * fHup(idn(thisP))                 &
            !     !                     )   
            !     !  elemR(thisP,outCol) = grav * ( fAup(idn(thisP)) - fAdn(iup(thisP)) ) &
            !     !                     * (                                              &
            !     !                     +   onefourthR  * fHdn(iup(thisP))                 &
            !     !                     -   threehalfR    * eHead(thisP)                     &
            !     !                     +   onefourthR  * fHup(idn(thisP))                 &
            !     !                     )                                         
            !     ! elemR(thisP,outCol) = onehalfR * grav &
            !     !                     *(+( fAup(idn(thisP)) - fAdn(iup(thisP)) ) * eHead(thisP) &
            !     !                       -( fHup(idn(thisP)) - fHdn(iup(thisP)) ) * eArea(thisP) &
            !     !                       +  fAup(idn(thisP)) * fHup(idn(thisP))                  &
            !     !                       -  fAdn(iup(thisP)) * fHdn(iup(thisP))                  &
            !     !                     )
            case default
                print *, 'CODE ERROR setting.Solver.MomentumSourceMethod type unknown for # ',&
                     setting%Solver%MomentumSourceMethod
                print *, 'which has key ',trim(reverseKey(setting%Solver%MomentumSourceMethod))
                call util_crashpoint(7298734)
        end select

    end subroutine ll_momentum_Ksource_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_momentum_source_CC (outCol, thisCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% Common source for momentum on channels and conduits
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: outCol, thisCol, Npack
            real(8)             :: delta
            real(8), pointer    :: fQ(:), fUdn(:), fUup(:), fAdn(:), fAup(:)
            real(8), pointer    :: fHdn(:), fHup(:), eKsource(:), grav
            integer, pointer    :: iup(:), idn(:), thisP(:)
            character(64)       :: subroutine_name = "ll_momentum_source_CC"
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (Npack < 1) return
            if (setting%Debug%File%lowlevel_rk2) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases
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
        !%------------------------------------------------------------------

        select case (setting%Solver%MomentumSourceMethod)
            case (T00)
                delta = zeroR
            case (T10)
                delta = onehalfR
            case (T20)
                delta = onesixthR
            case (T2L0)
                delta = onefourthR
            case (T10s2)
                delta = onehalfR
            case (TA1)
                delta = zeroR
            case (TA2)
                delta = zeroR
            case default
                print *, 'CODE ERROR setting.Solver.MomentumSourceMethod type unknown for # ',&
                            setting%Solver%MomentumSourceMethod
                print *, 'which has key ',trim(reverseKey(setting%Solver%MomentumSourceMethod))
                call util_crashpoint(7119223)
        end select

        elemR(thisP,outCol) = &
            fQ(iup(thisP)) * fUdn(iup(thisP)) - fQ(idn(thisP)) * fUup(idn(thisP)) &
            + grav * (oneR - delta) &
                *(  &
                    + fAdn(iup(thisP)) * fHdn(iup(thisP))  &
                    - fAup(idn(thisP)) * fHup(idn(thisP))  &
                    ) &
                + eKsource(thisP)

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%lowlevel_rk2) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine ll_momentum_source_CC
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
            real(8), pointer :: velocity(:), velocityold(:), mn(:), rh(:)
            real(8), pointer :: grav, ReversalFactor, SmallVelocity 
            integer, pointer :: thisP(:)
            character(64) :: subroutine_name = 'll_momentum_gammaCM_CC'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return    
        !%------------------------------------------------------------------
        !% Aliases
            thisP       => elemP(1:Npack,thisCol)
            velocity    => elemR(:,er_velocity)
            velocityOld => elemR(:,er_velocity_N0)
            rh          => elemR(:,er_HydRadius)
            mn          => elemR(:,er_ManningsN)
            grav        => setting%constant%gravity
            ReversalFactor => setting%Solver%ManningsN%FlowReversalFactor 
            SmallVelocity  => setting%Solver%ManningsN%SmallVelocity    

        !%---------------------------------------------------------------------
            
        where ((velocity(thisP) > SmallVelocity) .and. (velocityOld(thisP) > SmallVelocity))
            !% ---- standard Manning's n approach 
            !%      for consistent downstream flow
            elemR(thisP,outCol) =                                       &
                    grav * (mn(thisP)**twoR) * abs(velocity(thisP))     &
                    /                                                   &
                    ( rh(thisP)**fourthirdsR )    
        elsewhere
            !% --- HACK - increasing Manning's n on backflow or reversal
            !%     using average of last and this velocity along with an
            !%     increased Manning's n factor so that zero
            !%     velocity at one time step will still have a drag component.
            !%     This is necessary to damp backwards propagating waves with
            !%     small velocities that should have dissipation far in 
            !%     excess of that implied by Manning's n value
            elemR(thisP,outCol) =                                                      &
                    grav * ((ReversalFactor*mn(thisP))**twoR)                          &
                     * onehalfR * (abs(velocity(thisP)) + abs(velocityOld(thisP)))     &
                    /                                                                  &
                    ( rh(thisP)**fourthirdsR )
        endwhere                     

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
            !% ---- HACK the force main factors need to be moved to setting
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

            case (DarcyWeisbach)
                elemR(thisP,outCol) =                          &
                    (DWf(thisP) * abs(velocity(thisP)) )       &
                    / ( eightR * Afull(thisP)/Pfull(thisP)  )

            case default
                print *, 'CODE ERROR unexpected case default'
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
        !% Note this term is g h_L / L = KU/2L
        !% REVISED 20230202 so that only conduit minor loss is considered.
        !% The exit/entrance are now applied in the junction_elements only
        !% Note that for nJ2 and nBC connections the entry/exit losses have
        !% already been added to the Kconduit_MinorLoss for the first and
        !% last elements of the link.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: inoutCol, thisCol, Npack
            real(8), pointer    :: velocity(:), oneVec(:)
            real(8), pointer    :: Kconduit(:), length(:)
            integer, pointer    :: thisP(:)
            character(64)       :: subroutine_name = 'll_minorloss_friction_CC'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
            thisP    => elemP(1:Npack,thisCol)
            velocity => elemR(:,er_velocity)
            Kconduit => elemR(:,er_Kconduit_MinorLoss)
            length   => elemR(:,er_Length)
            oneVec   => elemR(:,er_ones)
        !%------------------------------------------------------------------
        !% ---- minor loss term (without gravity, which cancels out in derivation)
        elemR(thisP,inoutCol) = elemR(thisP,inoutCol)               &
                + abs(velocity(thisP))                              & 
                * (Kconduit(thisP))                                 &
               /                                                    &
               (twoR * length(thisP)) 

    end subroutine ll_minorloss_friction_gamma_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_momentum_solve_CC (outCol, thisCol, Npack, istep)
        !%------------------------------------------------------------------
        !% Description:
        !% Advance flowrate to n+1/2 for conduits and channel momentum 
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: outCol, thisCol, Npack, istep
            integer, pointer :: thisP(:)
            real(8), pointer :: delt, crk(:)
            real(8), pointer :: volumeLast(:), velocityLast(:), Msource(:), GammaM(:)
        !%------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
            thisP        => elemP(1:Npack,thisCOl)
            crk          => setting%Solver%crk2
            delt         => setting%Time%Hydraulics%Dt
            volumeLast   => elemR(:,er_Volume_N0)
            velocityLast => elemR(:,er_Velocity_N0)
            Msource      => elemR(:,er_SourceMomentum)
            GammaM       => elemR(:,er_GammaM)
        !%------------------------------------------------------------------

        elemR(thisP,outCol) =  &
                ( volumeLast(thisP) * velocityLast(thisP) + crk(istep) * delt * Msource(thisP) ) &
                / ( oneR + crk(istep) * delt * GammaM(thisP) )

    end subroutine ll_momentum_solve_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine ll_momentum_velocity_CC (inoutCol, thisCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% velocity for time march
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent (in) :: inoutCol,  thisCol, Npack
            integer, pointer     :: thisP(:)
            real(8), pointer     :: momentum(:), volume(:)
        !%------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
            thisP    => elemP(1:Npack,thisCol)
            !% --- the input integrated momentum is overwritten by the output velocity.
            momentum => elemR(:,inoutCol) !% overwritten
            volume   => elemR(:,er_Volume)
        !%------------------------------------------------------------------

        !% --- volume limiter preventing small volumes from generating flow spikes
        where ((abs(momentum(thisP) / elemR(thisP,er_Length))) &
                 .ge. (setting%Limiter%VolumeFractionInTimeStep  * volume(thisP) / setting%Time%Hydraulics%Dt))
            elemR(thisP,inoutCol) = sign((setting%Limiter%VolumeFractionInTimeStep  * volume(thisP) / setting%Time%Hydraulics%Dt), momentum(thisP) )
        elsewhere
            elemR(thisP,inoutCol) = momentum(thisP) / volume(thisP)
        endwhere

        !% --- minimum flowrate 
        if (setting%Limiter%Velocity%ZeroMinimumVelocitiesYN) then 
            where  ((abs(momentum(thisP) / elemR(thisP,er_Length))) < setting%Limiter%Velocity%Minimum)
                elemR(thisP,inoutCol) = zeroR 
            endwhere
        end if

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
            integer, pointer     :: thisP(:)
            real(8), pointer     :: flow(:)
        !%------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
            thisP    => elemP(1:Npack,thisCol)
            !% --- the inout column is either velocity, momentum, or flowrate
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
    subroutine ll_enforce_zerodepth_velocity (VelCol, thisCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% Enforces zero velocity on any element that is zerodepth, which
        !% prevents artificially-large thin-layer velocities from affecting
        !% face interpolation during filling of zerodepth
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: VelCol, thisCol, Npack
            integer, pointer    :: thisP(:)
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases:
            thisP => elemP(1:Npack,thisCol)
        !%------------------------------------------------------------------
        
        where (elemYN(thisP,eYN_isZeroDepth))
            elemR(thisP,VelCol) = zeroR
        endwhere
    
        !%------------------------------------------------------------------
        !% Closing:

    end subroutine ll_enforce_zerodepth_velocity
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module rk2_lowlevel