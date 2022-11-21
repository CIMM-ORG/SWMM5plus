module forcemain
    
    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use utility, only: util_sign_with_ones
    use utility_output
    use utility_crash, only: util_crashpoint

    implicit none

    !%-----------------------------------------------------------------------------
    !% Description:
    !% Functions for force main computations
    !%
    !% METHOD:
    !%
    !%

    private

    public forcemain_ManningsN

    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine forcemain_ManningsN ()
        !%------------------------------------------------------------------
        !% Description:
        !% Calls lower-level subroutines to set the equivalent Mannings n
        !% for face mains that have a free surface.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: thisPackCol, Npack, thisP(:)
            character(64)    :: subroutine_name = 'forcemain_ManningsN'
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (.not. setting%Solver%ForceMain%AllowForceMainTF) return
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------
        !%
        !% --- for Hazen-Williams
        if (.not. setting%Solver%ForceMain%HazenWilliams_equivalent_isSetTF) then
            !% --- the HW equivalent only needs to be called once and the equivalent
            !%     Manning's n is stored for re-use
            thisPackCol => col_elemP(ep_FM_HW_all)
            Npack       => npack_elemP(thisPackCol)
            if (Npack > 0) then
                call fm_equivalent_manningsN (thisPackCol,Npack,HazenWilliams)
            endif
            !% --- only check this once. The ep_FM_HW_all is a static pack, so
            !%     if Npack = 0 on the first time through it will always be zero
            setting%Solver%ForceMain%HazenWilliams_equivalent_isSetTF = .true.
        else 
            !% --- skip if already set
        end if

        !% --- for Darcy-Weisbach
        !%     Equivalent Mannings n is function of friction factor, so needs 
        !%     to be re-computed at each time step for the non-surcharged FM elements
        thisPackCol => col_elemP(ep_FM_dw_PSnonSurcharged)
        Npack       => npack_elemP(thisPackCol)
        if (Npack > 0) then
            call fm_dw_friction (thisPackCol, Npack)
            call fm_equivalent_manningsN (thisPackCol,Npack,DarcyWeisbach)
        endif

    end subroutine forcemain_ManningsN
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine fm_equivalent_manningsN (thisCol, Npack, fmMethod)  
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the equivalent Mannings N used in force mains that are
        !% not surcharged
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisCol, Npack, fmMethod
            integer, pointer    :: thisP(:)
            real(8), pointer    :: manningsN(:), slope(:), Afull(:), Pfull(:)
            real(8), pointer    :: HWcoef(:), Ffactor(:)
            real(8), pointer    :: slopeMin, grav
            real(8), parameter  :: HWfactorD = 0.8492d0
            real(8), parameter  :: HWslopeExp = 0.04d0
            real(8), parameter  :: HWhydradExp = 0.03667d0
            real(8), parameter  :: DWhydradExp = 0.1667d0 !% (1/6)
            character(64) :: subroutine_name = 'fm_equivalent_manningsN'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
            thisP      => elemP(1:Npack,thisCol)
            HWcoef     => elemSR(:,esr_Conduit_ForceMain_Coef)
            Ffactor    => elemSR(:,esr_Conduit_ForceMain_FrictionFactor)
            manningsN  => elemR(:,er_ManningsN)
            Afull      => elemR(:,er_FullArea)
            Pfull      => elemR(:,er_FullPerimeter)
            slope      => elemR(:,er_BottomSlope)
            slopeMin   => setting%Solver%ForceMain%minimum_slope
            grav       => setting%Constant%gravity
        !%------------------------------------------------------------------
        select case (fmMethod)
        case (HazenWilliams)
            !% N = Rh^0.037 / (0.85 * C * S^0.04)
            manningsN(thisP) = ( ((Afull(thisP) / Pfull(thisP)))**HWhydradExp ) &
                / (HWfactorD * HWcoef(thisP) * (max( slope(thisP), slopeMin )**HWslopeExp) )
        case (DarcyWeisbach)
            !% N = Rh^1/6  * sqrt( f / 8g )
            manningsN(thisP) = ( ((Afull(thisP) / Pfull(thisP)))**DWhydradExp ) &
                * sqrt( Ffactor(thisP) / (eightR * grav) )
        case default
            print *, 'CODE ERROR: unexpected case default'
            call util_crashpoint(779834)
        end select

        ! print *, 'in ',trim(subroutine_name)
        ! print *, elemR(1:4,er_ManningsN)

    end subroutine fm_equivalent_manningsN  
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine fm_dw_friction (thisCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% computes the friction factor for Darcy-Weisbach force mains
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisCol, Npack
            integer, pointer    :: thisP(:)
            real(8), pointer    :: Re(:), Re2(:), hydradius(:), Afull(:), Pfull(:)
            real(8), pointer    :: velocity(:), rough(:), Ffac(:), viscosity
            real(8), parameter  :: eCoef = 1.08108d0 !% roughness multiplier in SI (4/3.7)
            real(8), parameter  :: rCoef = 5.74d0  !% Reynolds number coef
            real(8), parameter  :: rExpon = 0.9d0  !% Reynolds number exponent
            character(64) :: subroutine_name = 'fm_dw_friction'
        !%--------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
            thisP     => elemP(1:Npack,thisCol)
            Re        => elemR(:,er_Temp01)
            Re2       => elemR(:,er_Temp02)
            hydradius => elemR(:,er_HydRadius)
            Afull     => elemR(:,er_FullArea)
            Pfull     => elemR(:,er_FullPerimeter)
            velocity  => elemR(:,er_Velocity)
            rough     => elemSR(:,esr_Conduit_ForceMain_Coef)
            Ffac      => elemSR(:,esr_Conduit_ForceMain_FrictionFactor)
            viscosity => setting%Constant%water_kinematic_viscosity
        !%------------------------------------------------------------------
        !% --- for Reynolds number, use Hydraulic Diameter = 4 * hydraulic radius
        !%     This uses the hydraulic radius associated with the full flow area
        !%     for consistency in the friction factor computation for equivalent 
        !%     manning's n computed at full pipe conditions.
        Re(thisP) = velocity(thisP) * fourR * (Afull(thisP) / Pfull(thisP)) / viscosity

        !% --- lower bound for Reynolds number following forcemain_getReynolds in EPA SWMM
        Re(thisP) = min(Re(thisP),tenR)

        !% --- setup for Re < 4000
        Re2(thisP) = max(Re(thisP),4000.d0)
     
        !% --- compute the friction factor for fully turbulent flow
        !%     This uses the full hydraulic radius for consistency in the derivation
        !%     of the equivalent Manning's n
        Ffac(thisP) = onefourthR                                             &
            / ( log10(                                                       &
                      (eCoef * rough(thisP) / (Afull(thisP) / Pfull(thisP))) &
                       + (rCoef / (Re2(thisP)**rExpon))                       &
                     )**2 ) 

        !% --- handle low Re following approach in EPA SWMM forcemain_getFricFactor
        where (Re(thisP) .le. 2000.d0)
            Ffac(thisP) = 64.d0 / Re(thisP)
        elsewhere ((Re(thisP) > 2000.d0) .and. (Re(thisP) < 4000.d0))
            Ffac(thisP) = 0.032d0 + (Ffac(thisP) - 0.032d0) * (Re(thisP) - 2000.d0) / 2000.d0
        elsewhere
            !% -- accept the f from full turbulence
        end where       
        
        !% --- reset temporary values
        Re = nullvalueR
        Re2 = nullvalueR

    end subroutine fm_dw_friction
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module forcemain
