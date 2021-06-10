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
        thisP => elemP(:1:Npack,thisCol)
        fQ    => faceR(:,fr_Flowrate)
        eQlat => elemR(:,er_FlowrateLateral)
        iup   => elemI(:,ei_Mface_uL)
        idn   => elemI(:,ei_Mface_dL)
        !%-----------------------------------------------------------------------------

        elemR(thisP,outCol) = fQ(iup(thisP)) - fQ(idn(thisP)) + eQlat(:)        

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
        !%-----------------------------------------------------------------------------
        
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
        !%-----------------------------------------------------------------------------  
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
                stop
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
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------
    
        !%-----------------------------------------------------------------------------
        !%             
    !%==========================================================================
    !% PRIVATE
    !%==========================================================================    
    
    
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------
    
        !%-----------------------------------------------------------------------------
        !%   

    
    
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------
    
        !%-----------------------------------------------------------------------------
        !%   

    !%==========================================================================
    !% END OF MODULE
    !%+=========================================================================
end module lowlevel_rk2