module adjust
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Ad hoc adjustments and limiters
    !%==========================================================================
    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use pack_mask_arrays, only: pack_small_or_zero_depth_elements, pack_CC_zeroDepth_interior_faces
    use preissmann_slot, only: slot_CC_Vshaped_adjust
    use utility
    use utility_crash

    implicit none

    private

    public :: adjust_element_toplevel
    public :: adjust_face_for_zero_setting_singular
    public :: adjust_limit_by_zerovalues  !% used in geometry
    public :: adjust_limit_by_zerovalues_singular  !% used in geometry   
    public :: adjust_Vfilter_CC
    public :: adjust_zero_and_small_depth_face    
    public :: adjust_zero_or_small_depth_identify_NEW

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine adjust_element_toplevel (elementType)
        !%------------------------------------------------------------------
        !% Description
        !% adjustments for zero depth or small depth
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: elementType !%, CC, JM, JB
            integer :: npack
        !%------------------------------------------------------------------

        !% --- CC ELEMENT AD HOC ADJUSTMENTS    
        !% --- identify zero depths (.true. is zero depth)
        call adjust_zero_or_small_depth_identify_NEW(elementType,.true.)

        if (setting%SmallDepth%useMomentumCutoffYN) then    
            !% --- identify small depths (.false. is small depth)
            call adjust_zero_or_small_depth_identify_NEW(elementType,.false.)                
        end if

        !% --- create packed arrays of zero and small depths
        call pack_small_or_zero_depth_elements (elementType,.true.)

        if (setting%SmallDepth%useMomentumCutoffYN) then   
            call pack_small_or_zero_depth_elements (elementType,.false.)
        end if

        !% --- adjust head, flowrate, and auxiliary values at zero depth
        call adjust_zerodepth_element_values (elementType) 

        select case (elementType)
            case (CC) 
                if (setting%SmallDepth%useMomentumCutoffYN) then
                    !% --- apply limiters to fluxes and velocity
                    !%     (.false. so that smalldepth fluxes are not set to zero)
                    call adjust_smalldepth_element_fluxes_CC (.false.)
                end if
                call adjust_limit_velocity_max (CC) 
            case (JB)
                call adjust_limit_velocity_max (JB) 
            case (JM)
                call adjust_limit_velocity_max (JM) 
            case default
                print *, 'CODE ERROR unexpected case default'
                call util_crashpoint(2250983)
        end select

    end subroutine adjust_element_toplevel
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine adjust_face_for_zero_setting_singular (iFidx)
        !%------------------------------------------------------------------
        !% Description:
        !% Input is a single face that is "closed" by a control action
        !% where the upstream element elemR(:,er_Setting) = 0.0
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: iFidx  !% index of the face
        !%------------------------------------------------------------------
        faceR(iFidx, fr_Flowrate)              = zeroR
        faceR(iFidx, fr_Flowrate_Conservative) = zeroR
        faceR(iFidx, fr_Velocity_d)            = zeroR
        faceR(iFidx, fr_Velocity_u)            = zeroR

    end subroutine adjust_face_for_zero_setting_singular
!%
!%==========================================================================
!%==========================================================================  
!%        
    subroutine adjust_limit_by_zerovalues (geocol, geozero, thisP, isVolume)
        !%------------------------------------------------------------------
        !% Description:
        !% This applies a zero value limiter for a packed array that is not
        !% a priori limited to zero value elements
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: geocol, thisP(:)
            real(8), intent(in) :: geozero
            logical, intent(in) :: isVolume
            real(8), pointer    :: geovalue(:), overflow(:)    
            character(64) :: subroutine_name = 'adjust_limit_by_zerovalues'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%adjust) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases
            geovalue     => elemR(:,geocol)
            overflow     => elemR(:,er_VolumeOverFlow)
        !%------------------------------------------------------------------

        if (isVolume) then
            !% --- we are gaining volume by resetting to the geozero (minimum),so
            !%    count this as a negative overflow
            where (geovalue(thisP) < geozero)
                overflow(thisP) = overflow(thisP) - (geozero - geovalue(thisP)) 
                geovalue(thisP) = geozero
            end where
        else
            where (geovalue(thisP) .le. geozero)
                geovalue(thisP) = geozero * 0.99d0
            endwhere
        end if
       
        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%adjust) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine adjust_limit_by_zerovalues
!%
!%==========================================================================  
!%==========================================================================  
!%    
    subroutine adjust_limit_by_zerovalues_singular &
        (eIdx, geocol, geozero, isVolume)
        !%------------------------------------------------------------------
        !% Description:
        !% Applies either the ZeroValue limiter (geozero) or zeroR as a lower limit to the
        !% geometry variable in elemR(:,geocol) for the single elemetn eIdx
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: geocol, eIdx
            real(8), intent(in) :: geozero
            logical, intent(in) :: isVolume
            real(8), pointer :: geovalue(:), overflow(:)        
            character(64) :: subroutine_name = 'adjust_limit_by_zerovalues_singular'
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (setting%Debug%File%adjust) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases:
            geovalue => elemR(:,geocol)
            overflow => elemR(:,er_VolumeOverFlow)
        !%------------------------------------------------------------------
        if (geovalue(eIdx) .le. geozero) then
            if (isVolume) then
                !% --- we are gaining volume by resetting to the geozero (minimum),so
                !%    count this as a negative overflow
                overflow(eIdx) = overflow(eIdx) - (geozero - geovalue(eIdx))
                geovalue(eIdx) = geozero
            else
                geovalue(eIdx) = geozero * 0.99d0
            end if
            
        end if

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%adjust) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
                
    end subroutine adjust_limit_by_zerovalues_singular
!%
!%========================================================================== 
!%==========================================================================
!%
    subroutine adjust_Vfilter_CC ()
        !%------------------------------------------------------------------
        !% Description:
        !% Performs ad-hoc adjustments that may be needed for stability
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: thisP(:), Npack
            real(8), pointer :: vMax

            character(64)    :: subroutine_name = 'adjust_Vfilter_CC'
        !%------------------------------------------------------------------
        !% Preleliminaries
            if (setting%Debug%File%adjust) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------   
        !% Aliases:
            Npack => npack_elemP(ep_CC)      
            if (Npack < 1) return 
            thisP => elemP(1:Npack,ep_CC)
            vMax       => setting%Limiter%Velocity%Maximum
        !%------------------------------------------------------------------

        !% --- ad hoc adjustments to flowrate 
        if (setting%Adjust%Flowrate%ApplyYN) then   
            select case (setting%Adjust%Flowrate%Approach)
                case (vshape)
                    !% --- suppress v-shape over face/element/face
                    call adjust_Vshaped_flowrate (thisP)
                case default
                    print *, 'CODE ERROR unknown setting.Adjust.Flowrate.Approach #',setting%Adjust%Flowrate%Approach
                    print *, 'which has key ',trim(reverseKey(setting%Adjust%Flowrate%Approach))
                    !stop 
                    call util_crashpoint( 4973)
                    !return
            end select
        else 
            !% --- no flow adjustment
        end if

        !% --- ad hoc adjustments to head
        !%     only affects head
        if (setting%Adjust%Head%ApplyYN) then   
            select case (setting%Adjust%Head%Approach)
                case (vshape_all_CC)
                    call adjust_Vshaped_head_CC(thisP,0,.false.)
                case (vshape_freesurface_CC)
                    call adjust_Vshaped_head_CC(thisP,eYN_isSurcharged,.false.)
                case (vshape_surcharge_CC)
                    call adjust_Vshaped_head_CC(thisP,eYN_isSurcharged,.true.)
                case default
                    print *,  'CODE ERROR unknown setting.Adjust.Head.Approach #',setting%Adjust%Head%Approach
                    print *, 'which has key ',trim(reverseKey(setting%Adjust%Head%Approach))
                    !stop 
                    call util_crashpoint( 9073)
                    !return
            end select
        else 
            !% --- nohead adjustment
        end if
 
        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%adjust) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine adjust_Vfilter_CC
!%
!%========================================================================== 
!%==========================================================================  
!%    
    subroutine adjust_zero_and_small_depth_face (ifixQCons)
        !%------------------------------------------------------------------
        !% Description:
        !% Top level control for face adjustment
        !% ifixQCons = .true. to use results to change the conservative face flux
        !%------------------------------------------------------------------
        !% Declarations:
            logical, intent(in) :: ifixQcons
        !%------------------------------------------------------------------

        if (setting%SmallDepth%useMomentumCutoffYN) then 
            call adjust_smalldepth_face_fluxes_CC      (ifixQCons)
            call adjust_smalldepth_face_fluxes_JMJB    (ifixQCons)
        end if

        call adjust_zerodepth_face_fluxes_CC   (ifixQCons)

        call adjust_zerodepth_face_fluxes_JMJB (ifixQCons)

        call adjust_JB_elem_flux_to_equal_face ()  
 
    end subroutine adjust_zero_and_small_depth_face
!%
!%========================================================================== 
!%==========================================================================
!%
    subroutine adjust_zero_or_small_depth_identify_NEW (elementType, isZero)    
        !%------------------------------------------------------------------
        !% Description
        !% identifies all the zero or small depths of element inType
        !% if isZero is true then zero depths are identified, else small depths
        !%------------------------------------------------------------------
        !% Declarations:
        !%------------------------------------------------------------------
            integer, intent(in) :: elementType 
            logical, intent(in) :: isZero
            logical, pointer :: thisDepth(:), otherDepth(:)
            integer, pointer :: elemType(:)
            real(8), pointer :: depth0
            real(8), pointer :: eDepth(:)
            real(8) :: lowcutoff
        !%------------------------------------------------------------------
        !% Aliases
            elemType     => elemI(:,ei_elementType)       
            eDepth       => elemR(:,er_Depth)

            if (isZero) then
                !% -- zero depth logical
                depth0     => setting%ZeroValue%Depth
                lowcutoff  = -huge(0.0d0)
                thisDepth  => elemYN(:,eYN_isZeroDepth)                
            else
                !% --- small depth logical
                depth0     => setting%SmallDepth%MomentumDepthCutoff
                lowcutoff  = setting%ZeroValue%Depth
                thisDepth  => elemYN(:,eYN_isSmallDepth)
            endif
        !%------------------------------------------------------------------

        where (elemType .eq. elementType)
            where ((eDepth .le. depth0) .and. (eDepth > lowcutoff))
                thisDepth  = .true.
            elsewhere
                thisDepth  = .false.
            endwhere
        endwhere

        ! print *, ' '
        ! print *, 'in adjust_zero...'
        ! print *, elemR(611,er_Depth)
        ! print *, ' '

    end subroutine adjust_zero_or_small_depth_identify_NEW
!%  
!%==========================================================================  
!% PRIVATE
!%==========================================================================  
!%
    subroutine adjust_limit_velocity_max (whichType) 
        !%------------------------------------------------------------------
        !% Description:
        !% employs velocity limiters and small volume treatments to limit 
        !% destabilizing large velocities.
        !%------------------------------------------------------------------  
        !% Declarations:
            integer, intent(in) :: whichType
            integer, pointer    :: thisCol, Npack, thisP(:)
            real(8), pointer    :: vMax, velocity(:)
            character(64)       :: subroutine_name = 'adjust_velocity'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (.not. setting%Limiter%Velocity%UseLimitMaxYN) return
            if (setting%Debug%File%adjust) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"          
        !%-------------------------------------------------------------------
        !% Aliases
            select case (whichType)
            case (CC)
                thisCol => col_elemP(ep_CC)
            case (JM)
                thisCol => col_elemP(ep_JM)
            case (JB)
                thisCol => col_elemP(ep_JB)
            case (Diag)
                thisCol => col_elemP(ep_Diag)
            case default
                print *, 'CODE ERROR -- unexpected case default'
            end select
            Npack     => npack_elemP(thisCol)
            if (Npack < 1) return
            thisP     => elemP(1:Npack,thisCol)
            velocity  => elemR(:,er_Velocity)
            vMax      => setting%Limiter%Velocity%Maximum
        !%------------------------------------------------------------------ 

        !% apply ad-hoc velocity limiter
        where (abs(velocity(thisP)) > vMax)
            velocity(thisP) = sign( 0.99d0 * vMax, velocity(thisP) )
        endwhere 

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%adjust) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine adjust_limit_velocity_max
!%
!%==========================================================================   
!%==========================================================================
!%
    subroutine adjust_zerodepth_element_values (whichType)  !% PRIVATE
        !% -----------------------------------------------------------------
        !% Description:
        !% thisCol must be one of the ZeroDepth packed arrays that identifies
        !% all the (near) zero depth locations.
        !% -----------------------------------------------------------------
            integer, intent(in)  :: whichType
            integer, pointer :: thisCol, npack, thisP(:)
        !% -----------------------------------------------------------------
        !% Preliminaries
            select case (whichType)
            case (CC)
                thisCol => col_elemP(ep_ZeroDepth_CC)
            case (JM)
                thisCol => col_elemP(ep_ZeroDepth_JM)
            case (JB)
                thisCol => col_elemP(ep_ZeroDepth_JB)
            case default
                print *, 'CODE ERROR -- unexpected case default'
                call util_crashpoint( 93287)
            end select
        !% -----------------------------------------------------------------
        !% Aliases
            npack   => npack_elemP(thisCol)
            if (npack < 1) return
            thisP   => elemP(1:npack,thisCol)
        !% -----------------------------------------------------------------
        !% --- only reset the depth if it is too small
        where (elemR(thisP,er_Depth) .le. setting%ZeroValue%Depth)
            elemR(thisP,er_Depth)    = setting%ZeroValue%Depth * 0.99d0
        end where

        elemR(thisP,er_Area)         = setting%ZeroValue%Area
        elemR(thisP,er_dHdA)         = oneR / setting%ZeroValue%TopWidth
        elemR(thisP,er_EllDepth)     = setting%ZeroValue%Depth * 0.99d0
        elemR(thisP,er_Flowrate)     = zeroR 
        elemR(thisP,er_FroudeNumber) = zeroR
        elemR(thisP,er_Perimeter)    = setting%ZeroValue%TopWidth + setting%ZeroValue%Depth
        elemR(thisP,er_HydRadius)    = setting%ZeroValue%Area / (setting%ZeroValue%TopWidth + setting%ZeroValue%Depth)
        elemR(thisP,er_TopWidth)     = setting%ZeroValue%TopWidth
        elemR(thisP,er_Velocity)     = zeroR
        elemR(thisP,er_WaveSpeed)    = zeroR
        elemR(thisP,er_Head)         = setting%ZeroValue%Depth*0.99d0 + elemR(thisP,er_Zbottom)
    
        !% --- only reset volume when it gets negative
        where (elemR(thisP,er_Volume) < zeroR)
            elemR(thisP,er_Volume) = setting%ZeroValue%Volume 
        end where
    
    end subroutine adjust_zerodepth_element_values 
!%  
!%==========================================================================   
!%==========================================================================
!%
    subroutine adjust_smalldepth_element_fluxes_CC (isZeroFlux)  
        !% -----------------------------------------------------------------
        !% Description
        !% uses the ep_SmallDepth_CC_ALLtm pack to set the velocity and
        !% flowrate on all small volumes
        !%
        !% HACK -- this needs to be reviewed and debugged before using
        !% -----------------------------------------------------------------
        !% Declarations:
            logical, intent(in) :: isZeroFlux !% sets small depth fluxes to zero
            integer, pointer :: thisCol
            integer, pointer :: npack, thisP(:), fdn(:), fup(:)
            real(8), pointer :: Area(:), CMvelocity(:), CMvelocity2(:)
            real(8), pointer :: Flowrate(:), HydRadius(:), ManningsN(:)
            real(8), pointer :: VelocityN0(:), Velocity(:), VelocityBlend(:), svRatio(:)
            real(8), pointer :: Head(:)
            real(8), pointer :: SmallVolume(:), Volume(:), faceFlow(:), faceFlowCons(:)
            real(8), pointer :: fHead_u(:), fHead_d(:), Length(:), oneArray(:)
            integer, target  :: pset(2)
            real(8)          :: psign(2)
            integer :: ii
            character(64) :: subroutine_name = 'adjust_smalldepth_element_fluxes_CC'
        !% -----------------------------------------------------------------
        !% Preliminaries:   
            if (.not. setting%SmallDepth%useMomentumCutoffYN) return
            thisCol => col_elemP(ep_SmallDepth_CC)
        !% -----------------------------------------------------------------    
        !% Aliases    
            Area          => elemR(:,er_Area)
            CMvelocity    => elemR(:,er_SmallVolume_CMvelocity) 
            Flowrate      => elemR(:,er_Flowrate)
            Head          => elemR(:,er_Head)
            HydRadius     => elemR(:,er_HydRadius)
            Length        => elemR(:,er_Length)
            ManningsN     => elemR(:,er_SmallVolume_ManningsN)  
            oneArray      => elemR(:,er_ones) 
            SmallVolume   => elemR(:,er_SmallVolume)  
            svRatio       => elemR(:,er_SmallVolumeRatio)
            !% --- use old velocity to prevent unreasonably large RK2 mid-step from affecting adjustment
            VelocityN0    => elemR(:,er_Velocity_N0) 
            Velocity      => elemR(:,er_Velocity)
            CMvelocity2   => elemR(:,er_Temp02)
            Volume        => elemR(:,er_Volume)
            fHead_d       => faceR(:,fr_Head_d)
            fHead_u       => faceR(:,fr_Head_u)
            fdn           => elemI(:,ei_Mface_dL)
            fup           => elemI(:,ei_Mface_uL)
        !% -----------------------------------------------------------------          
        npack   => npack_elemP(thisCol)
        if (npack < 1) return 
        thisP   => elemP(1:npack,thisCol)

        !% --- forced values to zero (typically for initial conditions)
        if (isZeroFlux) then
            elemR(thisP,er_Velocity) = zeroR
            elemR(thisP,er_Flowrate) = zeroR
            return 
        end if
        
        !% --- define the small volume ratio, 
        !%     limit to 1.0 needed for intermediate step where SV is being exceeded.
        svRatio(thisP) = min(Volume(thisP) / SmallVolume(thisP), oneR)  !% 20220122brh   
    
        !% use the larger of available ManningsN values
        ManningsN(thisP) = setting%SmallDepth%ManningsN
        ManningsN(thisP) = max(ManningsN(thisP), elemR(thisP,er_ManningsN))

        !% --- chezy-manning velocity based on head slope
        ! CMvelocity(thisP) = sign( oneArray(thisP), fHead_d(fup(thisP)) - fHead_u(fdn(thisP))) &
        !     * ( HydRadius(thisP)**(twothirdR) )                                              &
        !     * sqrt(abs(fHead_d(fup(thisP)) - fHead_u(fdn(thisP))) / Length(thisP) )           &
        !     / ManningsN(thisP)

        !% --- using two CM velocities
        !%     chezy-manning velocity in upper part of element using head slope
        CMvelocity(thisP) = sign( oneArray(thisP), fHead_d(fup(thisP)) - Head(thisP)) &
            * ( HydRadius(thisP)**(twothirdR) )                                              &
            * sqrt(abs(fHead_d(fup(thisP)) - Head(thisP)) / (onehalfR * Length(thisP)) )           &
            / ManningsN(thisP)

        CMvelocity2(thisP) = sign( oneArray(thisP), Head(thisP) - fHead_u(fdn(thisP))) &
                * ( HydRadius(thisP)**(twothirdR) )                                              &
                * sqrt(abs(Head(thisP) - fHead_u(fdn(thisP))) / (onehalfR * Length(thisP)) )           &
                / ManningsN(thisP)        

        !% --- if opposite signs use the sum, otherwise take the smaller magnitude
        where (CMvelocity(thisP) * CMvelocity(thisP) < zeroR)
            CMvelocity(thisP) = CMvelocity(thisP) + CMvelocity2(thisP)
        elsewhere
            CMvelocity(thisP) = sign(min(abs(CMVelocity(thisP)), abs(CMVelocity2(thisP))),CMvelocity(thisP))
        endwhere
                
        !% --- if svRatio > 1, then blend with existing flowrate
        Flowrate(thisP) = Flowrate(thisP) * (svRatio(thisP) - oneR)         &
                         + svRatio(thisP) * CMvelocity(thisP) * Area(thisP)                

        !% --- new velocity  
        elemR(thisP,er_Velocity) = Flowrate(thisP) / Area(thisP)

        !% -----------------------------------------------------------------
        !% ARCHIVE OF TRIAL METHODS
        !% --- blend the computed velocity with CM velocity
        ! VelocityBlend(thisP) = svRatio(thisP) * VelocityN0(thisP) &
        !                     + (oneR - svRatio(thisP)) * CMvelocity(thisP)

        !% --- use original RK2 velocity when its magnitude is smaller than blended CM                          
        !% --- use the smaller magnitude of RK2 velocity or blend if both are positive)
        ! where ((VelocityBlend(thisP) .ge. zeroR) .and. (Velocity(thisP) .ge. zeroR ))
        !     VelocityBlend(thisP) = min(VelocityBlend(thisP),Velocity(thisP))      
        ! endwhere       
        
        !% --- use the smaller magnitude whene both are negative
        ! where ((VelocityBlend(thisP) < zeroR) .and. (Velocity(thisP) < zeroR ))
        !     VelocityBlend(thisP) = max(VelocityBlend(thisP),Velocity(thisP))      
        ! endwhere   

        !% --- new velocity 
        ! elemR(thisP,er_Velocity) = VelocityBlend(thisP)

   
    end subroutine adjust_smalldepth_element_fluxes_CC  
!%  
!%========================================================================== 
!%==========================================================================
!%
    subroutine adjust_zerodepth_face_fluxes_CC (ifixQCons)  !% PRIVATE
        !% -----------------------------------------------------------------
        !% Description:
        !% thisCol must be one of the ZeroDepth packed arrays that identifies
        !% all the (near) zero depth element locations.
        !% only applicable to CC,
        !% if ifixQCons = .true. then the conservative fluxes are adjusted
        !% -----------------------------------------------------------------
            logical, intent(in)  :: ifixQCons
            integer, pointer :: npack, thisCol, thisP(:), fdn(:), fup(:)
            real(8), pointer :: fQ(:), fQCons(:), fVel_u(:), fVel_d(:)
            real(8), pointer :: fArea_u(:), fArea_d(:)
        !% -----------------------------------------------------------------
        !% Aliases
            thisCol => col_elemP(ep_ZeroDepth_CC)
            npack   => npack_elemP(thisCol)
            if (npack < 1) return
            thisP   => elemP(1:npack,thisCol)
            fdn     => elemI(:,ei_Mface_dL)
            fup     => elemI(:,ei_Mface_uL)
            fQ      => faceR(:,fr_Flowrate)
            fQCons  => faceR(:,fr_Flowrate_Conservative)
            fVel_u  => faceR(:,fr_Velocity_u)
            fVel_d  => faceR(:,fr_Velocity_d)
            fArea_u => faceR(:,fr_Area_u)
            fArea_d => faceR(:,fr_Area_d)
        !% -----------------------------------------------------------------
        !% --- choose either zero or an inflow
        fQ(fup(thisP)) = max(fQ(fup(thisP)), zeroR)
        fQ(fdn(thisP)) = min(fQ(fdn(thisP)), zeroR)

        !% --- Set inflow from adjacent cell based on head gradient
        call adjust_faceflux_for_headgradient (thisP, setting%SmallDepth%MomentumDepthCutoff)

        !% --- reset the conservative fluxes
        if (ifixQCons) then
            fQCons(fup(thisP)) = fQ(fup(thisP))
            fQCons(fdn(thisP)) = fQ(fdn(thisP))
        end if

        !% --- reset velocities 
        !%     HACK: it would be better if we could enusre that all the
        !%     areas are greater than zero so we wouldn't need the where
        !%     statements. But as of 20230114 there were areas that lead
        !%     to NAN values
        where (fArea_d(fup(thisP)) > setting%ZeroValue%Area)
            fVel_d(fup(thisP)) = fQ(fup(thisP)) /  fArea_d(fup(thisP))
        elsewhere
            fVel_d(fup(thisP)) = zeroR
        endwhere

        where (fArea_u(fup(thisP)) > setting%ZeroValue%Area)
            fVel_u(fup(thisP)) = fQ(fup(thisP)) /  fArea_u(fup(thisP))
        elsewhere
            fVel_u(fup(thisP)) = zeroR
        endwhere

        where (fArea_d(fdn(thisP)) > setting%ZeroValue%Area)
        fVel_d(fdn(thisP)) = fQ(fdn(thisP)) /  fArea_d(fdn(thisP))
        elsewhere
            fVel_d(fdn(thisP)) = zeroR
        endwhere
        where ( fArea_u(fdn(thisP)) > setting%ZeroValue%Area)
            fVel_u(fdn(thisP)) = fQ(fdn(thisP)) /  fArea_u(fdn(thisP))
        elsewhere 
            fVel_u(fdn(thisP)) = zeroR
        endwhere

        !% -----------------------------------------------------------------
    end subroutine adjust_zerodepth_face_fluxes_CC        
!%  
!%==========================================================================   
!%==========================================================================
!%    
    subroutine adjust_zerodepth_face_fluxes_JMJB (ifixQCons)
        !%------------------------------------------------------------------
        !% Description:
        !% Sets zero depth values on branches and JM. Input column must
        !% be a JM set
        !% if ifixQCons = .true. then the conservative fluxes are adjusted
        !%------------------------------------------------------------------
        !% Declarations:
            logical, intent(in)  :: ifixQCons
            integer, pointer :: npack, thisCol, thisP(:), fup(:), fdn(:), isBranch(:)
            real(8), pointer :: fQ(:), fQCons(:)
            integer :: ii
        !%------------------------------------------------------------------
        !% Aliases
            thisCol => col_elemP(ep_ZeroDepth_JM)
            npack => npack_elemp(thisCol)
            if (npack < 1) return
            thisP => elemP(1:npack,thisCol)
            fup   => elemI(:,ei_Mface_uL)
            fdn   => elemI(:,ei_Mface_dL)
            isbranch => elemSI(:,esi_JunctionBranch_Exists)
            fQ    => faceR(:,fr_Flowrate)
            fQCons=> faceR(:,fr_Flowrate_Conservative)
        !%------------------------------------------------------------------

        do ii=1,max_branch_per_node,2
            where (isbranch(thisP+ii) .eq. oneI) 
                fQ(fup(thisP+ii  )) = max(fQ(fup(thisP+ii  )),zeroR)
            endwhere
            where (isbranch(thisP+ii+1) .eq. oneI) 
                fQ(fdn(thisP+ii+1)) = min(fQ(fdn(thisP+ii+1)),zeroR)
            endwhere 
        end do

        if (ifixQCons) then
            do ii=1,max_branch_per_node,2
                where (isbranch(thisP+ii) .eq. oneI) 
                    fQCons(fup(thisP+ii  )) = fQ(fup(thisP+ii  ))
                endwhere
                where (isbranch(thisP+ii+1) .eq. oneI)
                    fQCons(fdn(thisP+ii+1)) = fQ(fdn(thisP+ii+1))
                endwhere
            end do
        end if

    end subroutine adjust_zerodepth_face_fluxes_JMJB
!%  
!%==========================================================================   
!%==========================================================================
!%
    subroutine adjust_JB_elem_flux_to_equal_face ()
        !%------------------------------------------------------------------
        !% Description:
        !% makes the JB flowrate equal to the face flowrate
        !%------------------------------------------------------------------
        !% Declarations
            integer :: thisCol, ii
            integer, pointer :: npack, thisP(:), fdn(:), fup(:), isbranch(:)
            real(8), pointer :: eQ(:), fQ(:)
        !%------------------------------------------------------------------
        !% Aliases
            thisCol = ep_JM
            npack => npack_elemP(thisCol)
            if (npack < 1) return 
            thisP => elemP(1:npack,thisCol)
            eQ    => elemR(:,er_Flowrate)
            fQ    => faceR(:,fr_Flowrate)
            fdn   => elemI(:,ei_Mface_dL)
            fup   => elemI(:,ei_Mface_uL)
            isbranch => elemSI(:,esi_JunctionBranch_Exists)
        !%------------------------------------------------------------------
        
        !% --- assign the upstream face flux to the JB 
        !%     note that JB and face must have consistent fluxes or we
        !%     get conservation errors
        do ii=1,max_branch_per_node,2
            !% --- inflows to JB
            where ((fQ(fup(thisP+ii)) > zeroR) .and. (isbranch(thisP+ii) .eq. oneI))
               eQ(thisP + ii) = fQ(fup(thisP+ii))
            endwhere
            !% --- outflows from JB (or zero on face)
            where ((fQ(fup(thisP+ii)) .le. zeroR) .and. (isbranch(thisP+ii) .eq. oneI))
               fQ(fup(thisP+ii)) = eQ(thisP + ii)
            endwhere
        end do

        !% --- assign the downstream face flux to the JB 
        do ii=2,max_branch_per_node,2
            !% --- inflows to JB
            where ((fQ(fdn(thisP+ii))  <   zeroR ) .and. (isbranch(thisP+ii) .eq. oneI))
                eQ(thisP + ii) = fQ(fdn(thisP+ii))
            endwhere
            !% --- outflows from JB
            where ((fQ(fdn(thisP+ii)) .ge. zeroR ) .and.  (isbranch(thisP+ii) .eq. oneI))
                fQ(fdn(thisP+ii)) =  eQ(thisP + ii)
            endwhere
        end do
        
    end subroutine adjust_JB_elem_flux_to_equal_face
!%  
!%==========================================================================   
!%==========================================================================
!%   
    subroutine adjust_smalldepth_face_fluxes_CC (ifixQCons) 
        !%------------------------------------------------------------------
        !% Description:
        !% Sets the face values around an element where the ad-hoc
        !% small depth algorithm is used. Must be done after the small depth
        !% element values have been set.
        !% if ifixQCons = .true. then the conservative fluxes are adjusted
        !%------------------------------------------------------------------
        !% Declarations:
            logical, intent(in) ::  ifixQCons
            integer, pointer :: fdn(:), fup(:), thisP(:), thisCol, npack
            integer, pointer :: thisColJM, thisJM(:), npackJM, isbranch(:) !% 20220122brh
            real(8), pointer :: faceQ(:), elemQ(:), fQCons(:)
            real(8), pointer :: fVel_u(:), fVel_d(:)
            real(8), pointer :: faceHu(:), faceHd(:), faceAu(:), faceAd(:)
            real(8), pointer :: faceDu(:), faceDd(:)
            real(8), pointer :: elemH(:), elemL(:), elemVol(:)
            real(8), pointer :: dt, grav
            integer :: ii
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (.not. setting%SmallDepth%useMomentumCutoffYN) return
        !%------------------------------------------------------------------
        !% Aliases:
            thisCol   => col_elemP(ep_SmallDepth_CC)
            faceQ     => faceR(:,fr_Flowrate)
            faceHu    => faceR(:,fr_Head_u)
            faceHd    => faceR(:,fr_Head_d)
            faceAu    => faceR(:,fr_Area_u)
            faceAd    => faceR(:,fr_Area_d)
            fQCons    => faceR(:,fr_Flowrate_Conservative)
            elemQ     => elemR(:,er_Flowrate)
            elemH     => elemR(:,er_Head)
            elemL     => elemR(:,er_Length)
            fVel_u    => faceR(:,fr_Velocity_u)
            fVel_d    => faceR(:,fr_Velocity_d)
            
            elemVol   => elemR(:,er_Volume_N0)
            isbranch  => elemSI(:,esi_JunctionBranch_Exists)
            dt        => setting%Time%Hydraulics%Dt
            grav      => setting%constant%gravity
            fdn       => elemI(:,ei_Mface_dL)
            fup       => elemI(:,ei_Mface_uL)
        !%------------------------------------------------------------------
        !% Handle the CC elements
        npack => npack_elemP(thisCol)
        if (npack > 0) then
            thisP  => elemP(1:npack,thisCol)
        
            where (elemQ(thisP) .ge. zeroR)
                !% --- flow in downstream direction
                !%     downstream face value is minimum of the face value or element value
                faceQ(fdn(thisP)) = min(elemQ(thisP)     , faceQ(fdn(thisP)) )      
                !%     upstream face value is either the inflow from face or zero
                faceQ(fup(thisP)) = max(faceQ(fup(thisP)), zeroR)
                !% --- downstream outflow is limited to 1/3 the element volume 
                faceQ(fdn(thisP)) = min(faceQ(fdn(thisP)), elemVol(thisP) / (threeR * dt) )   !% 20220122brh
            elsewhere
                !% --- flow in upstream direction
                !%     downstream face value is inflow (negative face flow) or zero
                faceQ(fdn(thisP)) = min(faceQ(fdn(thisP)), zeroR )
                !%     upstream face value is the inflow (faceQ > 0) or the larger (smaller magnitude, closer
                !%     to zero) of the negative flowrate at face or element
                faceQ(fup(thisP)) = max(elemQ(thisP)     , faceQ(fup(thisP)) ) 
                !% --- upstream out flow is limited to 1/3 the element volume !% 20220122brh
                faceQ(fup(thisP)) = max(faceQ(fup(thisP)), -elemVol(thisP) / (threeR * dt))
            endwhere


            !% ARCHIVE
            !% THIS CAUSES PROBLEMS WITH DRAINING JUNCTION
            !% --- provide inflow rate from large head differences with small volume cells
            !%     Derived from the SVE momentum neglecting all terms except dQ/dt and gA dH/dx
            ! call adjust_faceflux_for_headgradient (thisP, setting%SmallDepth%MomentumDepthCutoff)

            if (ifixQCons) then
                !% --- update the conservative face Q
                fQCons(fdn(thisP)) = faceQ(fdn(thisP))
                fQCons(fup(thisP)) = faceQ(fup(thisP))
            end if

            !% --- reset velocities
            where ( faceAd(fup(thisP)) > setting%ZeroValue%Area)
                fVel_d(fup(thisP)) = faceQ(fup(thisP)) /  faceAd(fup(thisP))
            elsewhere
                fVel_d(fup(thisP)) = zeroR
            endwhere

            where (faceAu(fup(thisP)) > setting%ZeroValue%Area)
                fVel_u(fup(thisP)) = faceQ(fup(thisP)) /  faceAu(fup(thisP))
            elsewhere
                fVel_u(fup(thisP)) = zeroR
            endwhere

            where (faceAd(fdn(thisP)) > setting%ZeroValue%Area)
                fVel_d(fdn(thisP)) = faceQ(fdn(thisP)) /  faceAd(fdn(thisP))
            elsewhere 
                fVel_d(fdn(thisP)) = zeroR
            end where

            where (faceAu(fdn(thisP)) > setting%ZeroValue%Area)
                fVel_u(fdn(thisP)) = faceQ(fdn(thisP)) /  faceAu(fdn(thisP))
            elsewhere 
                fVel_u(fdn(thisP)) = zeroR
            endwhere

        else
            !% --- no CC elements, so no action
        end if

    end subroutine adjust_smalldepth_face_fluxes_CC
!%  
!%==========================================================================   
!%==========================================================================
!%
    subroutine adjust_smalldepth_face_fluxes_JMJB (ifixQCons)
        !%------------------------------------------------------------------
        !% Description:
        !% Sets the face values around an element where the ad-hoc
        !% small depth algorithm is used. Must be done after the small depth
        !% element values have been set.
        !% if ifixQCons = .true. then the conservative fluxes are adjusted
        !%------------------------------------------------------------------
        !% Declarations:
            logical, intent(in) ::  ifixQCons
            integer, pointer :: fdn(:), fup(:), thisP(:), npack
            integer, pointer :: thisColJM, thisJM(:), npackJM, isbranch(:) !% 20220122brh
            real(8), pointer :: faceQ(:), elemQ(:), fQCons(:)
            real(8), pointer :: fVel_u(:), fVel_d(:)
            real(8), pointer :: faceHu(:), faceHd(:), faceAu(:), faceAd(:)
            real(8), pointer :: faceDu(:), faceDd(:)
            real(8), pointer :: elemH(:), elemL(:), elemVol(:)
            real(8), pointer :: dt, grav
            integer :: ii
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (.not. setting%SmallDepth%useMomentumCutoffYN) return
        !%------------------------------------------------------------------
        !% Aliases:
            thisColJM => col_elemP(ep_SmallDepth_JM) 
            faceQ     => faceR(:,fr_Flowrate)
            faceHu    => faceR(:,fr_Head_u)
            faceHd    => faceR(:,fr_Head_d)
            faceAu    => faceR(:,fr_Area_u)
            faceAd    => faceR(:,fr_Area_d)
            fQCons    => faceR(:,fr_Flowrate_Conservative)
            elemQ     => elemR(:,er_Flowrate)
            elemH     => elemR(:,er_Head)
            elemL     => elemR(:,er_Length)
            fVel_u    => faceR(:,fr_Velocity_u)
            fVel_d    => faceR(:,fr_Velocity_d)
  
            elemVol   => elemR(:,er_Volume_N0)
            isbranch  => elemSI(:,esi_JunctionBranch_Exists)
            dt        => setting%Time%Hydraulics%Dt
            grav      => setting%constant%gravity
            fdn       => elemI(:,ei_Mface_dL)
            fup       => elemI(:,ei_Mface_uL)
        !%------------------------------------------------------------------
        !% --- Handle the JM elements
        npackJM => npack_elemP(thisColJM)
        if (npackJM > 0) then
            thisJM => elemP(1:npackJM,thisColJM)

            do ii = 1,max_branch_per_node,2
                where ((elemQ(thisJM+ii) .ge. zeroR) .and. (isbranch(thisJM+ii) .eq. oneI))
                    !% --- flow in downstream direction in upstream branch
                    !%     upstream face value is either the inflow from face or zero if face is outflow
                    faceQ(fup(thisJM+ii)) = max(faceQ(fup(thisJM+ii)), zeroR)
                endwhere
                where ((elemQ(thisJM+ii)   <  zeroR) .and. (isbranch(thisJM+ii) .eq. oneI))
                    !% --- flow in upstream direction in upstream branch is the inflow (faceQ >0) or
                    !%     the larger (closer to zero of the negative flowrate at face or element)
                    faceQ(fup(thisJM+ii)) = max(faceQ(fup(thisJM+ii)), elemQ(thisJM+ii)) 
                    !% --- outflow (-faceQ) is limited to 1/3 volume in junction
                    faceQ(fup(thisJM+ii)) = max(faceQ(fup(thisJM+ii)), -elemVol(thisJM+ii) / (threeR * dt) )
                endwhere

                where ((elemQ(thisJM+1+ii) .ge. zeroR) .and. (isbranch(thisJM+1+ii) .eq. oneI))
                    !% --- flow downstream direction in downstream branch 
                    !%     downstream face value is the smaller of the face or branch values
                    faceQ(fdn(thisJM+1+ii)) = min(faceQ(fdn(thisJM+1+ii)), elemQ(thisJM+1+ii)) 
                    !% --- outflow (+faceQ) is limited to 1/3 volume in junction
                    faceQ(fdn(thisJM+1+ii)) = min(faceQ(fdn(thisJM+1+ii)), elemVol(thisJM+1+ii) / (threeR * dt) )
                endwhere
                where ((elemQ(thisJM+1+ii)  <   zeroR) .and. (isbranch(thisJM+1+ii) .eq. oneI))
                    !% --- flow upstream direction in downstream branch
                    !%     face flow is the inflow (negative direction) from downstream or zero.
                    faceQ(fdn(thisJM+1+ii)) = min(faceQ(fdn(thisJM+1+ii)),zeroR) * (real(isbranch(thisJM+1+ii),8))
                endwhere

                if (ifixQCons) then
                    !% --- update the conservative face Q
                    where (isbranch(thisJM+ii) .eq. oneI) 
                        fQCons(fup(thisJM+ii  )) = faceQ(fup(thisJM+ii))
                    endwhere
                    where (isbranch(thisJM+1+ii) .eq. oneI)
                        fQCons(fdn(thisJM+1+ii)) = faceQ(fdn(thisJM+1+ii))
                    endwhere
                end if

                !% --- reset velocities
                where ((faceAd(fup(thisJM+ii)) > setting%ZeroValue%Area) .and. (isbranch(thisJM+ii) .eq. oneI))
                    fVel_d(fup(thisJM+ii  )) = faceQ(fup(thisJM+ii  )) /  faceAd(fup(thisJM+ii  )) 
                endwhere
                where ((faceAd(fup(thisJM+ii)) .le. setting%ZeroValue%Area) .and. (isbranch(thisJM+ii) .eq. oneI))
                    fVel_d(fup(thisJM+ii  )) = zeroR
                endwhere

                where ((faceAu(fup(thisJM+ii)) > setting%ZeroValue%Area) .and. (isbranch(thisJM+ii) .eq. oneI))
                    fVel_u(fup(thisJM+ii  )) = faceQ(fup(thisJM+ii  )) /  faceAu(fup(thisJM+ii  ))
                endwhere
                where ((faceAu(fup(thisJM+ii)) .le. setting%ZeroValue%Area) .and. (isbranch(thisJM+ii) .eq. oneI))
                    fVel_u(fup(thisJM+ii  )) = zeroR
                endwhere

                where ((faceAd(fdn(thisJM+1+ii))  >   setting%ZeroValue%Area)  .and. (isbranch(thisJM+1+ii) .eq. oneI))
                    fVel_d(fdn(thisJM+1+ii)) = faceQ(fdn(thisJM+1+ii)) /  faceAd(fdn(thisJM+1+ii))
                endwhere
                where ((faceAd(fdn(thisJM+1+ii)) .le. setting%ZeroValue%Area)  .and. (isbranch(thisJM+1+ii) .eq. oneI))
                    fVel_d(fdn(thisJM+1+ii)) = zeroR
                endwhere

                where ((faceAu(fdn(thisJM+1+ii))  >   setting%ZeroValue%Area) .and. (isbranch(thisJM+1+ii) .eq. oneI))
                    fVel_u(fdn(thisJM+1+ii)) = faceQ(fdn(thisJM+1+ii)) /  faceAu(fdn(thisJM+1+ii))
                endwhere
                where ((faceAu(fdn(thisJM+1+ii)) .le. setting%ZeroValue%Area) .and. (isbranch(thisJM+1+ii) .eq. oneI))
                    fVel_u(fdn(thisJM+1+ii)) = zeroR
                endwhere
            end do
        end if

    end subroutine adjust_smalldepth_face_fluxes_JMJB
!%  
!%==========================================================================
!%==========================================================================   
!%  
    subroutine adjust_Vshaped_flowrate (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Removes V-shape between faces and element center by averaging
        !% the face fluxes
        !%------------------------------------------------------------------ 
            integer, intent(in) :: thisP(:)  
            integer, pointer :: mapUp(:), mapDn(:)
            real(8), pointer :: coef, vMax, Qlateral(:), Vcoef(:)
            real(8), pointer :: faceFlow(:), elemFlow(:), elemVel(:)
            real(8), pointer :: faceAu(:), faceAd(:)
            real(8), pointer ::  w_uQ(:), w_dQ(:), elemArea(:), Vvalue(:)
            real(8), pointer :: elemDepth(:), multiplier, smallDepth
            character(64) :: subroutine_name = 'adjust_Vshaped_flowrate'
        !%------------------------------------------------------------------
        !% Preliminaries    
            if (setting%Adjust%Flowrate%Coef .le. zeroR) return
            if (setting%Debug%File%adjust) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------
        !% Aliases        
            coef => setting%Adjust%Flowrate%Coef
            if (coef .le. zeroR) return
        
            mapUp    => elemI(:,ei_Mface_uL)
            mapDn    => elemI(:,ei_Mface_dL)   
            faceFlow => faceR(:,fr_Flowrate)
            faceAu   => faceR(:,fr_Area_u)
            faceAd   => faceR(:,fr_Area_d)  
            elemFlow => elemR(:,er_Flowrate)    
            elemVel  => elemR(:,er_Velocity)
            elemArea => elemR(:,er_Area)
            elemDepth=> elemR(:,er_Depth)
            w_uQ     => elemR(:,er_InterpWeight_uQ)
            w_dQ     => elemR(:,er_InterpWeight_dQ)
            Vvalue   => elemR(:,er_Temp01)
            Vcoef    => elemR(:,er_Temp02)
            Qlateral => elemR(:,er_FlowrateLateral)
            
            multiplier => setting%Adjust%Flowrate%SmallDepthMultiplier
            smallDepth => setting%SmallDepth%MomentumDepthCutoff
            vMax       => setting%Limiter%Velocity%Maximum
        !%-----------------------------------------------------------------
        !% --- setting%coef is the blending adjustment (between 0.0 and 1.0)
        !%     if coef == 1 then the V-shape element flowrate is replaced by 
        !%     average of its faces. If coef < 1 then average value is blended
        !%     with the present value of Q.

        !% --- Vcoef is the coefficient adjusted for local conditions
        Vcoef(thisP) = coef    

        !% --- find the cells that are deep enough to use the V filter
        Vvalue(thisP) = elemDepth(thisP) / (multiplier * smallDepth)

        where (Vvalue(thisP) > oneR)
            Vvalue(thisP) = oneR
        elsewhere
            Vvalue(thisP) = zeroR
            Vcoef(thisP)  = zeroR
        endwhere

        !% --- eliminate cells that have no upstream inflow  !% TEST 20221021 brh
        where (faceR(mapUp(thisP),fr_Flowrate) .eq. zeroR)
            Vcoef(thisP) = zeroR
        end where

        !% ARCHIVE METHOD
        !% --- Reducing V-filter when Qlateral is large  
        !%     HACK the fraction below should be replaced with a coefficient
        ! where (Qlateral(thisP) > onefourthR * abs(elemFlow(thisP)))
        !     Vcoef(thisP)  = Vcoef(thisP) * (onefourthR * abs(elemFlow(thisP)) / Qlateral(thisP))**2
        !     Vvalue(thisP) = zeroR     
        ! endwhere

        !% --- the Vvalue returns...
        !%    -1.0 if the element Q is between the face Q (not v-shaped)
        !%    +1.0 if the element Q is outside of the two face Q (v-shaped)
        !%     0.0 if the element Q is equal to one of the face Q (not v-shaped)
        !%     0.0 if the depth is too shallow
        Vvalue(thisP) =  (util_sign_with_ones_or_zero(faceFlow(mapUp(thisP)) - elemFlow(thisP)))      &
                        *(util_sign_with_ones_or_zero(faceFlow(mapDn(thisP)) - elemFlow(thisP)))      &
                        * Vvalue(thisP)     

        where (Vvalue(thisP) .le. zeroR)
            Vcoef(thisP) = zeroR
        endwhere

        !% --- blending velocity so that so that area does not matter
        ! elemVel(thisP)  =  (oneR - Vcoef(thisP)) * elemVel(thisP) &
        !         + Vcoef(thisP) * onehalfR * (faceR(mapDn(thisP),fr_Velocity_u) + faceR(mapUp(thisP),fr_Velocity_d))

        !% ARCHIVE METHOD
        !% --- blend the element and face-average flow rates
        elemFlow(thisP)  =  (oneR - Vcoef(thisP)) * elemFlow(thisP) &
                + Vcoef(thisP) * onehalfR * (faceflow(mapDn(thisP)) + faceflow(mapUp(thisP)))

        !% HACK: belnd the face areas as well
        ! elemArea(thisP)  =  (oneR - Vcoef(thisP)) * elemArea(thisP) &
        !         + Vcoef(thisP) * onehalfR * (faceAu(mapDn(thisP)) + faceAd(mapUp(thisP)))
        !% --- reset the velocity      
        elemVel(thisP) = elemFlow(thisP) / elemArea(thisP)   


        !% ARCHIVE METHOD
        ! where (Vvalue(thisP) > zeroR)
        !     !% simple linear interpolation
        !     elemFlow(thisP)  =  (oneR - coef) * elemFlow(thisP) &
        !         + coef * onehalfR * (faceflow(mapDn(thisP)) + faceflow(mapUp(thisP)))
        !     !% reset the velocity      
        !     elemVel(thisP) = elemFlow(thisP) / elemArea(thisP)   
        ! endwhere       
        !% --- reset for high velocity (typically due to small area)
        where ((abs(elemVel(thisP)) > vMax) .and. (Vcoef(thisP) > zeroR))
            elemVel(thisP)  = sign( 0.99d0 * vMax, elemVel(thisP) )
            elemFlow(thisP) = elemVel(thisP) * elemArea(thisP)
        endwhere 
        

        !%------------------------------------------------------------------
        !% Closing
            !% --- clear the temporary Vvalue array
            Vvalue(thisP) = nullvalueR
            if (setting%Debug%File%adjust) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine adjust_Vshaped_flowrate
!%
!%==========================================================================  
!%==========================================================================  
!%
    subroutine adjust_Vshaped_head_CC (thisP, eYNcol, eYNvalue)
        !%-------------------------------------------------------------------
        !% Description:
        !% Adjusts head for V-shaped conditions and fix depth/area. 
        !% Note that this breaks the 2-way relationship between depth and volume
        !% That is, after this point the depth and head are the "effective"
        !% values and not the average values associated with the volume.
        !% The eYNcolumn is a column in the elemYN that is used for screening
        !% which elements will be adjusted (typically, this is eYN_isSurcharged)
        !% If eYNcol = 0 then no screening is used. If a valid eYNcol is
        !% provided, then the eYNvalue is used to determine whether the .false.
        !% or the .true. elements are adjusted.
        !%-------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisP(:), eYNcol
            logical, intent(in) :: eYNvalue
            integer, pointer ::  fUp(:), fDn(:)
            real(8), pointer :: coef, fHup(:), fHdn(:), eHead(:), Vvalue(:)
            real(8), pointer :: fAup(:), fAdn(:), eArea(:)
        !%-------------------------------------------------------------------
        !% Preliminaries
        !%-------------------------------------------------------------------
        !% Aliases
            coef     => setting%Adjust%Head%Coef
            fUp      => elemI(:,ei_Mface_uL)
            fDn      => elemI(:,ei_Mface_dL)    
            fHup     => faceR(:,fr_Head_u)  
            fHdn     => faceR(:,fr_Head_d)  
            fAup     => faceR(:,fr_Area_u)  
            fAdn     => faceR(:,fr_Area_d)            
            eHead    => elemR(:,er_Head)    
            eArea    => elemR(:,er_Area)    
            Vvalue   => elemR(:,er_Temp01)
        !%-------------------------------------------------------------------

        !% --- discriminator that determines which cells are adjusted.
        Vvalue = zeroR 

        if (eYNcol == 0) then 
            !% --- check all elements in thisP
            Vvalue(thisP) = oneR
        else
            !% --- use the logical set to limit the checked elements
            where (elemYN(thisP,eYNcol) == eYNvalue)
                Vvalue(thisP) = oneR
            endwhere
        end if

        if (sum(Vvalue) == zeroR) return 

        !% --- identify the V-shape locations (Vvalue = 1)
        Vvalue(thisP) =  (util_sign_with_ones_or_zero(fHdn(fUp(thisP)) - eHead(thisP)))      &
                        *(util_sign_with_ones_or_zero(fHup(fDn(thisP)) - eHead(thisP)))      &
                        * Vvalue(thisP)   

        !% --- adjust where needed
        where (Vvalue(thisP) > zeroR)    
            !% --- simple linear combination depending on coef
            !%     note if coef==1 then this becomes the average of the face values
            eHead(thisP)  =  (oneR - coef) * eHead(thisP) &
               + coef * onehalfR * (fHup(fDn(thisP)) + fHdn(fUp(thisP)))

            !% ARCHIVE METHOD
            !% --- adjust depth
            ! elemR(thisP,er_Depth) = eHead(thisP) - elemR(thisP,er_Zbottom) 
            
            !% --- adjust area as a average rather than trying to use geometry
            ! eArea(thisP)  = (oneR - coef) * eArea(thisP)                    &
            !     + coef * onehalfR * (fAup(fDn(thisP)) + fAdn(fUp(thisP)))

            !% --- NOTE: volume is NOT adjusted.
        end where

        if (setting%Solver%PreissmannSlot%useSlotTF) then
            !% adjust v-shaped slots
            call slot_CC_Vshaped_adjust (thisP,er_Temp01)
        end if

      
    end subroutine adjust_Vshaped_head_CC
!%    
!%==========================================================================  
!%==========================================================================  
!%
    subroutine adjust_faceflux_for_headgradient (thisP, thisMomentumDepthCutoff)
        !%------------------------------------------------------------------
        !% Description
        !% Adjusts the face flux for a head gradient when the depth of the
        !% face is more than twice the input thisMomentumDepthCutoff
        !% Should only be applied to faces of zero depth or small depth cells  
        !% thisP should be a packed set of element indexes that are either
        !% small depth or zero depth elements.
        !% Uses the faceR(:,er_FlowrateMax) and faceR(:,er_FlowrateMin)
        !% to limit the allowable flowrate driven by head
        !%------------------------------------------------------------------
        !% Declarations
            integer,  intent(in) :: thisP(:)
            real(8),  intent(in) :: thisMomentumDepthCutoff
            integer, pointer :: fdn(:), fup(:)
            real(8), pointer :: faceQ(:), elemQ(:) 
            real(8), pointer :: faceHu(:), faceHd(:), faceAu(:), faceAd(:)
            real(8), pointer :: faceDu(:), faceDd(:), faceQmax(:), faceQmin(:)
            real(8), pointer :: elemH(:), elemL(:) 
            real(8), pointer :: dt, grav
            integer :: ii
        !%------------------------------------------------------------------
        !% Aliases:
            faceQ     => faceR(:,fr_Flowrate)
            faceHu    => faceR(:,fr_Head_u)
            faceHd    => faceR(:,fr_Head_d)
            faceAu    => faceR(:,fr_Area_u)
            faceAd    => faceR(:,fr_Area_d)
            faceDu    => faceR(:,fr_Depth_u)
            faceDd    => faceR(:,fr_Depth_d)
            faceQmax  => faceR(:,fr_FlowrateMax)
            faceQmin  => faceR(:,fr_FlowrateMin)
            elemH     => elemR(:,er_Head)
            elemL     => elemR(:,er_Length)
            
            dt        => setting%Time%Hydraulics%Dt
            grav      => setting%constant%gravity
            fdn       => elemI(:,ei_Mface_dL)
            fup       => elemI(:,ei_Mface_uL)
        !%------------------------------------------------------------------
        !% --- For the downstream face, dH/dx < 0 leads to a negative Q
        !%     Only applies where head gradient implies flow into the small volume and the
        !%     depth at the face is twice the small depth cutoff

        where ( (elemH(thisP) < faceHu(fdn(thisP)) ) &
                .and. &
                (faceDu(fdn(thisP)) > twoR * thisMomentumDepthCutoff) )

            !% --- value based on head gradient
            faceQ(fdn(thisP)) = min(                                                 &
                faceQ(fdn(thisP)),                                                   &
                dt * grav * faceAu(fdn(thisP)) * (elemH(thisP) - faceHu(fdn(thisP))) &
                / (onehalfR * (elemL(thisP)))                                        &
                )

            !% --- limit by available volume flowrate that empties downstream element
            faceQ(fdn(thisP)) = max( faceQ(fdn(thisP)), faceQmin(fdn(thisP)))   
        end where

        !% --- for the upstream face dH/dx > 0 leads to a positive Q
        !%     Only applies where head gradient implies flow into the small volume and the
        !%     depth on the face is twice the small depth cutoff
        where ( (elemH(thisP) < faceHd(fup(thisP)) ) &
                .and. &
                (faceDd(fup(thisP)) > twoR * thisMomentumDepthCutoff) )

            !% --- value based on head gradient
            faceQ(fup(thisP)) = max(                                                 &
                faceQ(fup(thisP)),                                                   &
                dt * grav * faceAd(fup(thisP)) * (faceHd(fup(thisP)) - elemH(thisP)) &
                / (onehalfR * (elemL(thisP)))                                        &
                )

           !% --- limit by available volume flowrate that empties upstream element
            faceQ(fup(thisP)) = min( faceQ(fup(thisP)), faceQmax(fup(thisP)))
        end where

    end subroutine adjust_faceflux_for_headgradient
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module adjust