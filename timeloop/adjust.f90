module adjust

    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use pack_mask_arrays, only: pack_small_and_zero_depth_elements
    use utility
    use utility_crash

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Makes ad hoc adjustments to ensure stability in hydraulic solution
    !%
    !% METHOD:
    !% 
    !%

    private

    public :: adjust_face_for_zero_setting
    public :: adjust_face_for_zero_setting_singular

    public :: adjust_zero_and_small_depth_elem
    public :: adjust_zero_and_small_depth_face
    public :: adjust_Vfilter

    public :: adjust_limit_by_zerovalues  !% used in geometry
    public :: adjust_limit_by_zerovalues_singular  !% used in geometry
    public :: adjust_limit_velocity_max
    !public :: adjust_JB_elem_flux_to_equal_face
    

    public :: adjust_zerodepth_identify_all
    ! public :: adjust_zerodepth_element_values 
    ! public :: adjust_zerodepth_face_fluxes_CC
    ! public :: adjust_zerodepth_face_fluxes_JMJB

    public :: adjust_smalldepth_identify_all
    ! public :: adjust_smalldepth_element_fluxes
    ! public :: adjust_smalldepth_face_fluxes
   

    contains
!%==========================================================================
!% PUBLIC
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
    subroutine adjust_face_for_zero_setting ()
        !%------------------------------------------------------------------
        !% Description:
        !% Sets all zero fluxes on downstream faces of "closed" CC elements
        !% i.e. where elemR(:,er_Setting) = 0.0
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: ptype, npack, thisP(:), dFace(:)
        !%------------------------------------------------------------------
        !% Aliases:    
            ptype => col_elemP(ep_CC_isClosedSetting)
            npack => npack_elemP(ptype)
            dFace => elemI(:,ei_Mface_dL)
        !%------------------------------------------------------------------    
        !% Preliminaries    
            if (npack < 1) return
        !%------------------------------------------------------------------ 
        !% --- elements that are closed (er_Setting = 0.0)        
        thisP => elemP(1:npack,ep_CC_isClosedSetting)

        !% --- force flows and velocities to zero
        faceR(dface(thisP), fr_Flowrate)              = zeroR
        faceR(dface(thisP), fr_Flowrate_Conservative) = zeroR
        faceR(dface(thisP), fr_Velocity_d)            = zeroR
        faceR(dface(thisP), fr_Velocity_u)            = zeroR

    end subroutine adjust_face_for_zero_setting
!%
!%==========================================================================
!%==========================================================================  
!%
    subroutine adjust_zero_and_small_depth_elem (whichTM, isreset)
        !%------------------------------------------------------------------
        !% Description:
        !% Top level adjustment routine for zero and small depth conditions
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: whichTM
            logical, intent(in) :: isreset !% true means the zero/small packs are reset
            integer, pointer   :: thisCol_CC, thisCol_JM
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:   
        !%------------------------------------------------------------------
     
        if (isreset) then
            call adjust_zerodepth_identify_all ()
           
            call adjust_smalldepth_identify_all ()
            
            call pack_small_and_zero_depth_elements (whichTM)
            
        end if

        call adjust_zerodepth_element_values (whichTM, CC) 

            ! call util_CLprint('-------------AAAA')
        
        call adjust_zerodepth_element_values (whichTM, JM) 

            ! call util_CLprint('-------------BBBB')
        
        call adjust_smalldepth_element_fluxes (whichTM)

            ! call util_CLprint('-------------CCCC')
        
        call adjust_limit_velocity_max (whichTM) 

            ! call util_CLprint('-------------DDDD')

        !%------------------------------------------------------------------
        !% Closing:

    end subroutine adjust_zero_and_small_depth_elem 
!%
!%==========================================================================
!%==========================================================================  
!%    
    subroutine adjust_zero_and_small_depth_face (whichTM, ifixQCons)
        !%------------------------------------------------------------------
        !% Description:
        !% Top level control for face adjustment
        !% ifixQCons = .true. to use results to change the conservative face flux
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: whichTM
            logical, intent(in) :: ifixQcons
            integer, pointer :: thisCol_CC, thisCol_JM
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------
        call adjust_smalldepth_face_fluxes     (whichTM,ifixQCons)
        call adjust_zerodepth_face_fluxes_CC   (whichTM,ifixQCons)
        call adjust_zerodepth_face_fluxes_JMJB (whichTM,ifixQCons)
        call adjust_JB_elem_flux_to_equal_face (whichTM) !% 20220123brh
    
        !%------------------------------------------------------------------
        !% Closing:
 
    end subroutine adjust_zero_and_small_depth_face
!%
!%========================================================================== 
!%==========================================================================
!%
    subroutine adjust_Vfilter (whichTM)
        !%------------------------------------------------------------------
        !% Description:
        !% Performs ad-hoc adjustments that may be needed for stability
        !%------------------------------------------------------------------
            integer, intent(in) :: whichTM  !% indicates which Time marching method
            character(64) :: subroutine_name = 'adjust_Vfilter'
        !%------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%adjust) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------   
                
        !% ad hoc adjustments to flowrate 
        if (setting%Adjust%Flowrate%ApplyYN) then   
            select case (setting%Adjust%Flowrate%Approach)
            case (vshape)
                !% suppress v-shape over face/element/face
                call adjust_Vshaped_flowrate (whichTM)
            case default
                print *, 'CODE ERROR: unknown setting.Adjust.Flowrate.Approach #',setting%Adjust%Flowrate%Approach
                print *, 'which has key ',trim(reverseKey(setting%Adjust%Flowrate%Approach))
                !stop 
                call util_crashpoint( 4973)
                !return
            end select
        end if
        
        !% ad hoc adjustments to head
        if (setting%Adjust%Head%ApplyYN) then          
            select case (setting%Adjust%Head%Approach)
            case (vshape_surcharge_only)
                call adjust_Vshaped_head_surcharged (whichTM)
            case default
                print *,  'CODE ERROR: unknown setting.Adjust.Head.Approach #',setting%Adjust%Head%Approach
                print *, 'which has key ',trim(reverseKey(setting%Adjust%Head%Approach))
                !stop 
                call util_crashpoint( 9073)
                !return
            end select
        end if
        
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine adjust_Vfilter
!%
!%========================================================================== 
!%==========================================================================  
!%        
    subroutine adjust_limit_by_zerovalues (geocol, geozero, thisCol, isVolume)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% This applies a zero value limiter for a packed array that is not
        !% a priori limited to zero value elements
        !%-----------------------------------------------------------------------------
        !logical, intent(in) :: isreset
        integer, intent(in) :: geocol, thisCol
        real(8), intent(in) :: geozero
        logical, intent(in) :: isVolume
        integer, pointer :: Npack, thisP(:)
        real(8), pointer :: geovalue(:), overflow(:)    
        !logical, pointer :: NearZeroDepth(:),   isSmallDepth(:)   
        character(64) :: subroutine_name = 'adjust_limit_by_zerovalues'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        Npack        => npack_elemP(thisCol)  
        geovalue     => elemR(:,geocol)
        overflow     => elemR(:,er_VolumeOverFlow)
        !isZeroDepth  => elemYN(:,eYN_isZeroDepth)
        !isSmallDepth => elemYN(:,eYN_isSmallDepth)
        !%-----------------------------------------------------------------------------

        if (Npack > 0) then
            thisP    => elemP(1:Npack,thisCol)
            if (isVolume) then
                !% --- we are gaining volume by resetting to the geozero (minimum),so
                !%    count this as a negative overflow
                where (geovalue(thisP) < geozero)
                    overflow(thisP) = overflow(thisP) - (geozero - geovalue(thisP)) 
                    geovalue(thisP) = geozero
                end where
            else
                where (geovalue(thisP) < geozero)
                    geovalue(thisP) = geozero
                endwhere
            end if
            ! if (isreset) then
            !     NearZeroDepth(thisP) = .false.
            !     where (geovalue(thisP) <= geozero)
            !         NearZeroDepth(thisP) = .true.
            !         isSmallDepth(thisP) = .false.
            !     end where
            ! end if
        end if    

        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine adjust_limit_by_zerovalues
!%
!%==========================================================================  
!%==========================================================================  
!%    
    subroutine adjust_limit_by_zerovalues_singular (eIdx, geocol, geozero, isVolume)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Applies either the ZeroValue limiter (geozero) or zeroR as a lower limit to the
        !% geometry variable in elemR(:,geocol) for the single elemetn eIdx
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: geocol, eIdx
        real(8), intent(in) :: geozero
        logical, intent(in) :: isVolume
        real(8), pointer :: geovalue(:), overflow(:)        
        character(64) :: subroutine_name = 'adjust_limit_by_zerovalues_singular'
        !if (crashYN) return
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        geovalue => elemR(:,geocol)
        overflow => elemR(:,er_VolumeOverFlow)
        !%-----------------------------------------------------------------------------
        if (geovalue(eIdx) < geozero) then
            if (isVolume) then
                !% --- we are gaining volume by resetting to the geozero (minimum),so
                !%    count this as a negative overflow
                overflow(eIdx) = overflow(eIdx) - (geozero - geovalue(eIdx))
            end if
            geovalue(eIdx) = geozero
        end if

        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine adjust_limit_by_zerovalues_singular
!%
!%==========================================================================  
!%==========================================================================  
!%
    subroutine adjust_limit_velocity_max (whichTM)
        !%------------------------------------------------------------------
        !% Description:
        !% employs velocity limiters and small volume treatments to limit 
        !% destabilizing large velocities.
        !%------------------------------------------------------------------  
            integer, intent(in) :: whichTM
            integer, pointer :: thisCol_all, Npack, thisP(:)
            real(8), pointer :: vMax, velocity(:)
            character(64) :: subroutine_name = 'adjust_velocity'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (.not. setting%Limiter%Velocity%UseLimitMaxYN) return
            !if (crashYN) return
            if (setting%Debug%File%adjust) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            !% small velocity adjustment should only be done for CC elements
            !% since the velocity is solved there,
            select case (whichTM)
                case (ALLtm)
                    thisCol_all => col_elemP(ep_CC_ALLtm)
                case (ETM)
                    thisCol_all => col_elemP(ep_CC_ETM)       
                case (AC)
                    thisCol_all => col_elemP(ep_CC_AC)        
                case default
                    print *, 'CODE ERROR: time march type unknown for # ', whichTM
                    print *, 'which has key ',trim(reverseKey(whichTM))
                    !stop 
                    call util_crashpoint( 8368)
                    !return
            end select
        !%-------------------------------------------------------------------
        !% Aliases
            Npack     => npack_elemP(thisCol_all)
            if (Npack < 1) return
            thisP     => elemP(1:Npack,thisCol_all)
            velocity  => elemR(:,er_Velocity)
            vMax      => setting%Limiter%Velocity%Maximum
        !%------------------------------------------------------------------ 
        !% apply ad-hoc velocity limiter
        where (abs(velocity(thisP)) > vMax)
            velocity(thisP) = sign( 0.99d0 * vMax, velocity(thisP) )
        endwhere 
        !%------------------------------------------------------------------
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine adjust_limit_velocity_max
!%
!%==========================================================================  
!%==========================================================================
!%
    ! subroutine adjust_smallvolume_facevalues (whichTM)
    !     !%----------------------------------------------------------------------
    !     !% Description
    !     !% sets the face fluxes to the CM for outflows from a smallvolume element
    !     !% call with thisCol=0 to do all elements
    !     !%----------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: whichTM
    !         logical, pointer :: isSmallVol(:)
    !         integer, pointer :: fdn(:), fup(:), thisP(:)
    !         integer, pointer :: Npack, thisCol
    !         integer :: ii
    !         real(8), pointer :: fQ(:), eQ(:)
    !     !%----------------------------------------------------------------------
    !     !% Preliminaries
    !     !%----------------------------------------------------------------------
    !     !% Aliases
    !         select case (whichTM)
    !         case (ALLtm)
    !             thisCol => col_elemP(ep_CC_ALLtm)
    !         case (ETM)
    !             thisCol => col_elemP(ep_CC_ETM)    
    !         case (AC)
    !             thisCol => col_elemP(ep_CC_AC)       
    !         case default
    !             print *, 'error, default case should not be reached.'
    !             stop 
    !               call util_crashpoint( 83655)
    !         end select
    !         isSmallVol => elemYN(:,eYN_isSmallDepth)
    !         fdn        => elemI(:,ei_Mface_dL)
    !         fup        => elemI(:,ei_Mface_uL)
    !         fQ         => faceR(:,fr_Flowrate)
    !         eQ         => elemR(:,er_Flowrate)
    !         Npack      => npack_elemP(thisCol)
    !         if (Npack .le. 0) return
    !         thisP => elemP(1:Npack,thisCol)
    !     !%----------------------------------------------------------------------
        
    !     !% --- a small volume can have an inflow from faces interpolated faces
    !     !%     but the outflow is limited to the CM outflow
    !     !%     Note the following will seg fault if used on JB or JM cells

    !     where ( isSmallVol(thisP) )
    !         !% set the downstream face flowrate to zero or keep the value if it is negative (inflow)
    !         fQ(fdn(thisP)) = onehalfR * ( fQ(fdn(thisP)) - abs( fQ(fdn(thisP)) ) )
    !     end where

    
    !     !% if the value is now zero, use the element value for the face flux
    !     where ( (isSmallVol(thisP)) .and. ( abs(faceR(fdn(thisP),fr_flowrate)) < setting%Eps%Velocity) )
    !         fQ(fdn(thisP)) = eQ(thisP)
    !     end where


    !     where (isSmallVol(thisP)) 
    !         !% set the upstream face flowrate to zero or keep the value if it is positive (inflow)    
    !         fQ(fup(thisP)) = onehalfR * ( fQ(fup(thisP)) + abs( fQ(fup(thisP))) )
    !     end where

    !     !% if the value is now zero, use the element value for the face flux
    !     where ( (isSmallVol(thisP)) .and. ( abs(faceR(fup(thisP),fr_flowrate)) < setting%Eps%Velocity) )
    !         fQ(fup(thisP)) = eQ(thisP)
    !     end where

    ! end subroutine adjust_smallvolume_facevalues
!%   
!%==========================================================================    
!%==========================================================================
!%
    subroutine adjust_zerodepth_identify_all()
        !%------------------------------------------------------------------
        !% Description
        !% identifies all the zero depths (without regard to TM or Diagnostic)
        !%------------------------------------------------------------------
        !% Declarations:
        !%------------------------------------------------------------------
            logical, pointer :: isZeroDepth(:), isAdHocQ(:), isSmallDepth(:)
            real(8), pointer :: depth0
            real(8), pointer :: eDepth(:)
        !%------------------------------------------------------------------
        !% Aliases
            depth0       => setting%ZeroValue%Depth
            eDepth       => elemR(:,er_Depth)
            isZeroDepth  => elemYN(:,eYN_isZeroDepth)
            isSmallDepth => elemYN(:,eYN_isSmallDepth)
        !%------------------------------------------------------------------

        isZeroDepth = .false.

        where (eDepth .le. depth0)
            isZeroDepth   = .true.
            isSmallDepth = .false.
        endwhere

        !%----------------------------------------------------------------------
        !% Closing
    end subroutine adjust_zerodepth_identify_all
!%  
!%==========================================================================   
!%==========================================================================
!%
    subroutine adjust_smalldepth_identify_all()
        !%------------------------------------------------------------------
        !% Description
        !% identifies all the small depths (without regard to TM or Diagnostic)
        !%------------------------------------------------------------------
        !% Declarations:
        !%------------------------------------------------------------------
            logical, pointer :: isZeroDepth(:), isAdHocQ(:), isSmallDepth(:)
            real(8), pointer :: depth0, depthS
            real(8), pointer :: eDepth(:)
        !%------------------------------------------------------------------
        !% Preliminaries
        !%------------------------------------------------------------------
        !% Aliases
            depth0       => setting%ZeroValue%Depth
            depthS       => setting%SmallDepth%DepthCutoff
            eDepth       => elemR(:,er_Depth)
            isSmallDepth => elemYN(:,eYN_isSmallDepth)
        !%------------------------------------------------------------------

        isSmallDepth = .false.

        where ((eDepth .le. depthS) .and. (eDepth > depth0) )
            isSmallDepth   = .true.
        endwhere

        !%----------------------------------------------------------------------
        !% Closing
    end subroutine adjust_smalldepth_identify_all
!%  
!%==========================================================================   
!%==========================================================================
!%
    subroutine adjust_zerodepth_element_values (whichTM, whichType)
        !% -----------------------------------------------------------------
        !% Description:
        !% thisCol must be one of the ZeroDepth packed arrays that identifies
        !% all the (near) zero depth locations.
        !% -----------------------------------------------------------------
            integer, intent(in)  :: whichTM, whichType
            integer, pointer :: thisCol, npack, thisP(:)
        !% -----------------------------------------------------------------
        !% Preliminaries
            select case (whichTM)
            case (ALLtm)
                select case (whichType)
                case (CC)
                    !thisCol => col_elemP(ep_ZeroDepth_CC_ALLtm)
                    print *, 'CODE ERROR: AC Not implemented'
                    call util_crashpoint(329874)
                case (JM)
                    !thisCol => col_elemP(ep_ZeroDepth_JM_ALLtm)
                    print *, 'CODE ERROR: AC Not implemented'
                call util_crashpoint(6698723)
                case default
                    print *, 'CODE ERROR -- unexpected case default'
                    !stop 
                    call util_crashpoint( 94733)
                    !return
                end select
            case (ETM)
                select case (whichType)
                case (CC)
                    thisCol => col_elemP(ep_ZeroDepth_CC_ETM)
                case (JM)
                    thisCol => col_elemP(ep_ZeroDepth_JM_ETM)
                case default
                    print *, 'CODE ERROR -- unexpected case default'
                    !stop 
                    call util_crashpoint( 93287)
                    !return
                end select
            case (AC)
                select case (whichType)
                case (CC)
                    !thisCol => col_elemP(ep_ZeroDepth_CC_AC)
                    print *, 'CODE ERROR: AC Not implemented'
                    call util_crashpoint(229843)
                case (JM)
                    !thisCol => col_elemP(ep_ZeroDepth_JM_AC)
                    print *, 'CODE ERROR: AC Not implemented'
                    call util_crashpoint(6698723)
                case default
                    print *, 'CODE ERROR -- unexpected case default'
                    !stop 
                    call util_crashpoint( 22733)
                    !return
                end select
            case default
                print *, 'CODE ERROR: time march type unknown for # ', whichTM
                print *, 'which has key ',trim(reverseKey(whichTM))
                !stop 
                call util_crashpoint( 55873)
                !return
            end select
        !% -----------------------------------------------------------------
        !% Aliases
            npack   => npack_elemP(thisCol)
            if (npack < 1) return
            thisP   => elemP(1:npack,thisCol)
        !% -----------------------------------------------------------------

        !% only reset the depth if it is too small
        where (elemR(thisP,er_Depth) .le. setting%ZeroValue%Depth)
            elemR(thisP,er_Depth)    = setting%ZeroValue%Depth
        end where

        elemR(thisP,er_Area)         = setting%ZeroValue%Area
        elemR(thisP,er_dHdA)         = oneR / setting%ZeroValue%TopWidth
        elemR(thisP,er_ell)          = setting%ZeroValue%Depth
        elemR(thisP,er_Flowrate)     = zeroR
        elemR(thisP,er_FroudeNumber) = zeroR
        elemR(thisP,er_HydDepth)     = setting%ZeroValue%Depth
        elemR(thisP,er_Perimeter)    = setting%ZeroValue%TopWidth + setting%ZeroValue%Depth
        elemR(thisP,er_HydRadius)    = setting%ZeroValue%Area / (setting%ZeroValue%TopWidth + setting%ZeroValue%Depth)
        elemR(thisP,er_TopWidth)     = setting%ZeroValue%TopWidth
        elemR(thisP,er_Velocity)     = zeroR
        elemR(thisP,er_WaveSpeed)    = zeroR
        elemR(thisP,er_Head)    = setting%ZeroValue%Depth + elemR(thisP,er_Zbottom)
        !elemR(thisP,er_HeadAvg) = setting%ZeroValue%Depth + elemR(thisP,er_Zbottom)

        !% only reset volume when it gets too small
        where (elemR(thisP,er_Volume) < setting%ZeroValue%Volume)
            elemR(thisP,er_Volume) = (oneR + onetenthR) * setting%ZeroValue%Volume 
        end where
    
    end subroutine adjust_zerodepth_element_values 
!%  
!%==========================================================================   
!%==========================================================================
!%
    subroutine adjust_smalldepth_element_fluxes(whichTM)
        !% -----------------------------------------------------------------
        !% Description
        !% uses the ep_SmallDepth_CC_ALLtm pack to set the velocity and
        !% flowrate on all small volumes
        !% -----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: whichTM
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
            character(64) :: subroutine_name = 'adjust_smalldepth_element_fluxes'
        !% -----------------------------------------------------------------
        !% Preliminaries:   
            select case (whichTM)
            case (ALLtm)
                !thisCol => col_elemP(ep_SmallDepth_CC_ALLtm)
                print *, 'CODE ERROR: AC Not implemented'
                call util_crashpoint(66987253)
            case (ETM)
                thisCol => col_elemP(ep_SmallDepth_CC_ETM)
            case (AC)
                !thisCol => col_elemP(ep_SmallDepth_CC_AC)
                print *, 'CODE ERROR: AC Not implemented'
                call util_crashpoint(6698723)
            case default
                print *, 'CODE ERROR: time march type unknown for # ', whichTM
                   print *, 'which has key ',trim(reverseKey(whichTM))
                !stop 
                call util_crashpoint( 557345)
                !return
            end select
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
            !% use old velocity to prevent unreasonably large RK2 mid-step from affecting adjustment
            VelocityN0    => elemR(:,er_Velocity_N0) 
            Velocity      => elemR(:,er_Velocity)
            VelocityBlend => elemR(:,er_Temp01)
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
        
        !% --- define the small volume ratio, 
        !%     limit to 1.0 needed for intermediate step where SV is being exceeded.
        svRatio(thisP) = min(Volume(thisP) / SmallVolume(thisP), oneR)  !% 20220122brh   

        !print *, 'svRatio ', svRatio(iet(1:2))
    
        !% use the larger of available ManningsN values
        ManningsN(thisP) = setting%SmallDepth%ManningsN
        ManningsN(thisP) = max(ManningsN(thisP), elemR(thisP,er_ManningsN))   
        if (setting%Solver%ManningsN%useDynamicManningsN) then
            ManningsN(thisP) = max(ManningsN(thisP), elemR(thisP,er_ManningsN_Dynamic))
        end if

        ! print *, 'in ',trim(subroutine_name)
        ! print *, ManningsN(139), elemR(139,er_ManningsN_Dynamic), elemR(139,er_ManningsN)

        !print *, 'mannings n', ManningsN(iet(1:2))

        !print *, 'head diff ',fHead_d(fup(iet(1:2))) - fHead_u(fdn(iet(1:2)))
        !print *, 'hydradius ',HydRadius(iet(1:2))
        ! print *, 'face Up',fup(iet(1:2))
        ! print *, 'face Dn',fdn(iet(1:2))
        ! print *, 'head Upstream ',fhead_d(fup(iet(1:2)))
        ! print *, 'head Dnstream ',fhead_u(fdn(iet(1:2)))

        !% chezy-manning velocity based on head slope
        ! CMvelocity(thisP) = sign( oneArray(thisP), fHead_d(fup(thisP)) - fHead_u(fdn(thisP))) &
        !     * ( HydRadius(thisP)**(twothirdR) )                                              &
        !     * sqrt(abs(fHead_d(fup(thisP)) - fHead_u(fdn(thisP))) / Length(thisP) )           &
        !     / ManningsN(thisP)

        

        !% Trial 20220716brh -- using two CM velocities
        !% chezy-manning velocity in upper part of element using head slope
        CMvelocity(thisP) = sign( oneArray(thisP), fHead_d(fup(thisP)) - Head(thisP)) &
            * ( HydRadius(thisP)**(twothirdR) )                                              &
            * sqrt(abs(fHead_d(fup(thisP)) - head(thisP)) / (onehalfR * Length(thisP)) )           &
            / ManningsN(thisP)

        CMvelocity2(thisP) = sign( oneArray(thisP), Head(thisP) - fHead_u(fdn(thisP))) &
                * ( HydRadius(thisP)**(twothirdR) )                                              &
                * sqrt(abs(Head(thisP) - fHead_u(fdn(thisP))) / (onehalfR * Length(thisP)) )           &
                / ManningsN(thisP)    

        !% if opposite signs use the sum, otherwise take the smaller magnitude
        where (CMvelocity(thisP) * CMvelocity(thisP) < zeroR)
            CMvelocity(thisP) = CMvelocity(thisP) + CMvelocity2(thisP)
        elsewhere
            CMvelocity(thisP) = sign(min(abs(CMVelocity(thisP)), abs(CMVelocity2(thisP))),CMvelocity(thisP))
        endwhere
                
                
        ! print *, 'head dif  ',fHead_d(fup(iet(8))) - fHead_u(fdn(iet(8)))    
        ! print *, 'hydRadius ',HydRadius(iet(8))
        ! print *, 'CMvelocity', CMvelocity(iet(8))
        ! print *, 'velocity  ', Velocity(iet(8)), VelocityN0(iet(8))
        ! print *, 'SVratio   ',svRatio(iet(8))

        !% blend the computed velocity with CM velocity
        VelocityBlend(thisP) = svRatio(thisP) * VelocityN0(thisP) &
                            + (oneR - svRatio(thisP)) * CMvelocity(thisP)

        !% 20220716brh
        !% --- use original RK2 velocity when its magnitude is smaller than blended CM                          
        !% --- use the smaller magnitude of RK2 velocity or blend if both are positive)
        where ((VelocityBlend(thisP) .ge. zeroR) .and. (Velocity(thisP) .ge. zeroR ))
            VelocityBlend(thisP) = min(VelocityBlend(thisP),Velocity(thisP))      
        endwhere       
        
        !% use the smaller magnitude whene both are negative
        where ((VelocityBlend(thisP) < zeroR) .and. (Velocity(thisP) < zeroR ))
            VelocityBlend(thisP) = max(VelocityBlend(thisP),Velocity(thisP))      
        endwhere   

        ! print *, 'Velblend ', VelocityBlend(iet(8))

        !% new flowrate
        Flowrate(thisP) = Area(thisP) * VelocityBlend(thisP)

        !% new velocity  % 20220712brh
        elemR(thisP,er_Velocity) = VelocityBlend(thisP)

        ! print *, ' '
        ! print *, 'in adjust small depth element fluxes'
        ! print *, 'Flowrate ',Flowrate(iet(7))
        ! print *, ' '

        !% reset the temporary storage
        VelocityBlend(thisP) = nullvalueR

        !% -----------------------------------------------------------------
    end subroutine adjust_smalldepth_element_fluxes  
!%  
!%========================================================================== 
!%==========================================================================
!%
    subroutine adjust_zerodepth_face_fluxes_CC (whichTM,ifixQCons)
        !% -----------------------------------------------------------------
        !% Description:
        !% thisCol must be one of the ZeroDepth packed arrays that identifies
        !% all the (near) zero depth element locations.
        !% only applicable to CC, not JM or JB
        !% if ifixQCons = .true. then the conservative fluxes are adjusted
        !% -----------------------------------------------------------------
            integer, intent(in)  :: whichTM
            logical, intent(in)  :: ifixQCons
            integer, pointer :: npack, thisCol, thisP(:), fdn(:), fup(:)
            real(8), pointer :: fQ(:), fQCons(:), fVel_u(:), fVel_d(:)
            real(8), pointer :: fArea_u(:), fArea_d(:)
        !% -----------------------------------------------------------------
        !% Preliminaries:
            select case (whichTM)
            case (ALLtm)
                !thisCol => col_elemP(ep_ZeroDepth_CC_ALLtm)
                print *, 'CODE ERROR: AC Not implemented'
                call util_crashpoint(1048366)
            case (ETM)
                thisCol => col_elemP(ep_ZeroDepth_CC_ETM)
            case (AC)
                !thisCol => col_elemP(ep_ZeroDepth_CC_AC)
                print *, 'CODE ERROR: AC Not implemented'
                call util_crashpoint(823453)
            case default
                print *, 'CODE ERROR: time march type unknown for # ', whichTM
                print *, 'which has key ',trim(reverseKey(whichTM))
                !stop 
                call util_crashpoint( 22487)
                !return
            end select
            npack   => npack_elemP(thisCol)
            if (npack < 1) return
        !% -----------------------------------------------------------------
        !% Aliases
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
        call adjust_faceflux_for_headgradient (thisP, setting%SmallDepth%DepthCutoff)

            !call util_CLprint ('-------- after adjust-faceflux-for-headgradient in adjust zerodepth')

        !% --- reset the conservative fluxes
        if (ifixQCons) then
            fQCons(fup(thisP)) = fQ(fup(thisP))
            fQCons(fdn(thisP)) = fQ(fdn(thisP))
        end if

        !% --- reset velocities
        fVel_d(fup(thisP)) = fQ(fup(thisP)) /  fArea_d(fup(thisP))
        fVel_u(fup(thisP)) = fQ(fup(thisP)) /  fArea_u(fup(thisP))

        fVel_d(fdn(thisP)) = fQ(fdn(thisP)) /  fArea_d(fdn(thisP))
        fVel_u(fdn(thisP)) = fQ(fdn(thisP)) /  fArea_u(fdn(thisP))



        !% -----------------------------------------------------------------
    end subroutine adjust_zerodepth_face_fluxes_CC        
!%  
!%==========================================================================   
!%==========================================================================
!%    
    subroutine adjust_zerodepth_face_fluxes_JMJB (whichTM, ifixQCons)
        !%------------------------------------------------------------------
        !% Description:
        !% Sets zero depth values on branches and JM. Input column must
        !% be a JM set
        !% if ifixQCons = .true. then the conservative fluxes are adjusted
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)  :: whichTM
            logical, intent(in)  :: ifixQCons
            integer, pointer :: npack, thisCol, thisP(:), fup(:), fdn(:), isBranch(:)
            real(8), pointer :: fQ(:), fQCons(:)
            integer :: ii
        !%------------------------------------------------------------------
        !% Preliminaries:
            select case (whichTM)
            case (ALLtm)
                !thisCol => col_elemP(ep_ZeroDepth_JM_ALLtm)
                print *, 'CODE ERROR: AC Not implemented'
                call util_crashpoint(6698723)
            case (ETM)
                thisCol => col_elemP(ep_ZeroDepth_JM_ETM)
            case (AC)
                !thisCol => col_elemP(ep_ZeroDepth_JM_AC)
                print *, 'CODE ERROR: AC Not implemented'
                call util_crashpoint(64438723)
            case default
                print *, 'CODE ERROR: time march type unknown for # ', whichTM
                print *, 'which has key ',trim(reverseKey(whichTM))
                !stop 
                call util_crashpoint( 224238)
                !return
            end select
            npack => npack_elemp(thisCol)
            if (npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases:
            thisP => elemP(1:npack,thisCol)
            fup   => elemI(:,ei_Mface_uL)
            fdn   => elemI(:,ei_Mface_dL)
            isbranch => elemSI(:,esi_JunctionBranch_Exists)
            fQ    => faceR(:,fr_Flowrate)
            fQCons=> faceR(:,fr_Flowrate_Conservative)
        !%------------------------------------------------------------------
        do ii=1,max_branch_per_node,2
            fQ(fup(thisP+ii  )) = max(fQ(fup(thisP+ii  )),zeroR) * real(isbranch(thisP+ii  ),8)
            fQ(fdn(thisP+ii+1)) = min(fQ(fdn(thisP+ii+1)),zeroR) * real(isbranch(thisP+ii+1),8)    
        end do
        
        if (ifixQCons) then
            do ii=1,max_branch_per_node,2
                fQCons(fup(thisP+ii  )) = fQ(fup(thisP+ii  ))
                fQCons(fdn(thisP+ii+1)) = fQ(fdn(thisP+ii+1))
            end do
        end if
    
        !%------------------------------------------------------------------
        !% Closing:

    end subroutine adjust_zerodepth_face_fluxes_JMJB
!%  
!%==========================================================================   
!%==========================================================================
!%
    subroutine adjust_JB_elem_flux_to_equal_face (whichTM)
        !%------------------------------------------------------------------
        !% Description:
        !% makes the JB flowrate equal to the face flowrate
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: whichTM
            integer :: thisCol, ii
            integer, pointer :: npack, thisP(:), fdn(:), fup(:), isbranch(:)
            real(8), pointer :: eQ(:), fQ(:)
        !%------------------------------------------------------------------
        !% Preliminaries
            select case (whichTM)
            case (ETM)
                thisCol = ep_JM_ETM
            case default
                print *, 'CODE ERROR: time march type unknown for # ', whichTM
                print *, 'which has key ',trim(reverseKey(whichTM))
                !stop 
                call util_crashpoint( 398703)
                !return
            end select
            npack => npack_elemP(thisCol)
            if (npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
            thisP => elemP(1:npack,thisCol)
            eQ    => elemR(:,er_Flowrate)
            fQ    => faceR(:,fr_Flowrate)
            fdn   => elemI(:,ei_Mface_dL)
            fup   => elemI(:,ei_Mface_uL)
            isbranch => elemSI(:,esi_JunctionBranch_Exists)
        !%------------------------------------------------------------------
        
        !% assign the upstream face flux to the JB
        do ii=1,max_branch_per_node,2
            eQ(thisP + ii) = fQ(fup(thisP+ii)) * real(isbranch(thisP+ii),8)
        end do

        !% assign the downstream face flux to the JB
        do ii=2,max_branch_per_node,2
            eQ(thisP + ii) = fQ(fdn(thisP+ii)) * real(isbranch(thisP+ii),8)
        end do
        

    end subroutine adjust_JB_elem_flux_to_equal_face
!%  
!%==========================================================================   
!%==========================================================================
!%
    subroutine adjust_smalldepth_face_fluxes (whichTM, ifixQCons)
        !%------------------------------------------------------------------
        !% Description:
        !% Sets the face values around an element where the ad-hoc
        !% small depth algorithm is used. Must be done after the small depth
        !% element values have been set.
        !% if ifixQCons = .true. then the conservative fluxes are adjusted
        !%------------------------------------------------------------------
        !% Declarations:
            logical, intent(in) ::  ifixQCons
            integer, intent(in) :: whichTM
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
            select case (whichTM)
            case (ALLtm)
                !thisCol   => col_elemP(ep_SmallDepth_CC_ALLtm)
                !thisColJM => col_elemP(ep_SmallDepth_JM_ALLtm)
                print *, 'CODE ERROR: AC Not implemented'
                call util_crashpoint(6634112)
            case (ETM)
                thisCol   => col_elemP(ep_SmallDepth_CC_ETM)
                thisColJM => col_elemP(ep_SmallDepth_JM_ETM) 
            case (AC)
                !thisCol   => col_elemP(ep_SmallDepth_CC_AC)
                !thisColJM => col_elemP(ep_SmallDepth_JM_AC)
                print *, 'CODE ERROR: AC Not implemented'
                call util_crashpoint(6634723)
            case default
                print *, 'CODE ERROR: time march type unknown for # ', whichTM
                print *, 'which has key ',trim(reverseKey(whichTM))
                !stop 
                call util_crashpoint( 447833)
                !return
            end select
        !%------------------------------------------------------------------
        !% Aliases:
            faceQ     => faceR(:,fr_Flowrate)
            faceHu    => faceR(:,fr_Head_u)
            faceHd    => faceR(:,fr_Head_d)
            faceAu    => faceR(:,fr_Area_u)
            faceAd    => faceR(:,fr_Area_d)
            faceDu    => faceR(:,fr_HydDepth_u)
            faceDd    => faceR(:,fr_HydDepth_d)
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

            !% 20220531brh
            !% --- provide inflow rate from large head differences with small volume cells
            !%     Derived from the SVE momentum neglecting all terms except dQ/dt and gA dH/dx
            call adjust_faceflux_for_headgradient (thisP, setting%SmallDepth%DepthCutoff)

            if (ifixQCons) then
                !% update the conservative face Q
                fQCons(fdn(thisP)) = faceQ(fdn(thisP))
                fQCons(fup(thisP)) = faceQ(fup(thisP))
            end if

            !% --- reset velocities
            fVel_d(fup(thisP)) = faceQ(fup(thisP)) /  faceAd(fup(thisP))
            fVel_u(fup(thisP)) = faceQ(fup(thisP)) /  faceAu(fup(thisP))

            fVel_d(fdn(thisP)) = faceQ(fdn(thisP)) /  faceAd(fdn(thisP))
            fVel_u(fdn(thisP)) = faceQ(fdn(thisP)) /  faceAu(fdn(thisP))

        
        else
            !% no CC elements
        end if

        !% Handle the JM elements
        npackJM => npack_elemP(thisColJM)
        if (npackJM > 0) then
            thisJM => elemP(1:npackJM,thisColJM)
            do ii = 1,max_branch_per_node,2
                where (elemQ(thisJM+ii) .ge. zeroR)
                    !% --- flow in downstream direction in upstream branch
                    !%     upstream face value is either the inflow from face or zero if face is outflow
                    faceQ(fup(thisJM+ii)) = max(faceQ(fup(thisJM+ii)), zeroR) * (real(isbranch(thisJM+ii),8))
                elsewhere
                    !% --- flow in upstream direction in upstream branch is the inflow (faceQ >0) or
                    !%     the larger (closer to zero of the negative flowrate at face or element)
                    faceQ(fup(thisJM+ii)) = max(faceQ(fup(thisJM+ii)), elemQ(thisJM+ii)) * (real(isbranch(thisJM+ii),8))
                    !% --- outflow (-faceQ) is limited to 1/3 volume in junction
                    faceQ(fup(thisJM+ii)) = max(faceQ(fup(thisJM+ii)), -elemVol(thisJM+ii) / (threeR * dt) )
                endwhere
                where (elemQ(thisJM+1+ii) .ge. zeroR)
                    !% --- flow downstream direction in downstream branch 
                    !%     downstream face value is the smaller of the face or branch values
                    faceQ(fdn(thisJM+1+ii)) = min(faceQ(fdn(thisJM+1+ii)), elemQ(thisJM+1+ii))  * (real(isbranch(thisJM+1+ii),8))
                    !% --- outflow (+faceQ) is limited to 1/3 volume in junction
                    faceQ(fdn(thisJM+1+ii)) = min(faceQ(fdn(thisJM+1+ii)), elemVol(thisJM+1+ii) / (threeR * dt) )
                elsewhere
                    !% --- flow upstream direction in downstream branch
                    !%     face flow is the inflow (negative direction) from downstream or zero.
                    faceQ(fdn(thisJM+1+ii)) = min(faceQ(fdn(thisJM+1+ii)),zeroR) * (real(isbranch(thisJM+1+ii),8))
                endwhere

                if (ifixQCons) then
                    !% update the conservative face Q
                    fQCons(fup(thisJM+ii  )) = faceQ(fup(thisJM+ii))
                    fQCons(fdn(thisJM+1+ii)) = faceQ(fdn(thisJM+1+ii))
                end if

                !% --- reset velocities
                fVel_d(fup(thisJM+ii  )) = faceQ(fup(thisJM+ii  )) /  faceAd(fup(thisJM+ii  ))  
                fVel_u(fup(thisJM+ii  )) = faceQ(fup(thisJM+ii  )) /  faceAu(fup(thisJM+ii  ))

                fVel_d(fdn(thisJM+1+ii)) = faceQ(fdn(thisJM+1+ii)) /  faceAd(fdn(thisJM+1+ii))
                fVel_u(fdn(thisJM+1+ii)) = faceQ(fdn(thisJM+1+ii)) /  faceAu(fdn(thisJM+1+ii))

            end do
        end if
        
        !%------------------------------------------------------------------
        !% Closing:
    end subroutine adjust_smalldepth_face_fluxes

!%  
!%==========================================================================
!% PRIVATE
!%==========================================================================   
!%  
    subroutine adjust_Vshaped_flowrate (whichTM)
        !%------------------------------------------------------------------
        !% Description:
        !% Removes V-shape between faces and element center by averaging
        !% the face fluxes
        !%------------------------------------------------------------------   
            integer, intent(in) :: whichTM
            integer, pointer :: thisCol, Npack
            integer, pointer :: thisP(:), mapUp(:), mapDn(:)
            real(8), pointer :: coef, vMax, Qlateral(:), Vcoef(:)
            real(8), pointer :: faceFlow(:), elemFlow(:), elemVel(:)
            real(8), pointer ::  w_uQ(:), w_dQ(:), elemArea(:), Vvalue(:)
            real(8), pointer :: elemDepth(:), multiplier, smallDepth
            !logical, pointer :: isSmallDepth(:), isNearZeroDepth(:), 
            character(64) :: subroutine_name = 'adjust_Vshaped_flowrate'
        !%------------------------------------------------------------------
        !% Preliminaries    
            if (setting%Adjust%Flowrate%Coef .le. zeroR) return
            !if (crashYN) return       
            if (setting%Debug%File%adjust) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------
        !% Aliases        
            select case (whichTM)
            case (ALLtm)
                thisCol => col_elemP(ep_CC_ALLtm)
            case (ETM)
                thisCol => col_elemP(ep_CC_ETM)
            case (AC)
                thisCol => col_elemP(ep_CC_AC)
            case default
                print *, 'CODE ERROR: time march type unknown for # ', whichTM
                print *, 'which has key ',trim(reverseKey(whichTM))
                !stop 
                call util_crashpoint( 9239)
                !return
            end select

            coef => setting%Adjust%Flowrate%Coef
            if (coef .le. zeroR) return

            Npack => npack_elemP(thisCol)
            if (Npack .le. 0) return
            
            thisP    => elemP(1:Npack,thisCol)
            mapUp    => elemI(:,ei_Mface_uL)
            mapDn    => elemI(:,ei_Mface_dL)   
            faceFlow => faceR(:,fr_Flowrate)  
            elemFlow => elemR(:,er_Flowrate)    
            elemVel  => elemR(:,er_Velocity)
            elemArea => elemR(:,er_Area)
            elemDepth=> elemR(:,er_Depth)
            w_uQ     => elemR(:,er_InterpWeight_uQ)
            w_dQ     => elemR(:,er_InterpWeight_dQ)
            Vvalue   => elemR(:,er_Temp01)
            Vcoef    => elemR(:,er_Temp02)
            Qlateral => elemR(:,er_FlowrateLateral)
            !isSmallDepth   => elemYN(:,eYN_isSmallDepth)
            !isNearZeroDepth => elemYN(:,eYN_isZeroDepth)
            
            multiplier => setting%Adjust%Flowrate%SmallDepthMultiplier
            smallDepth => setting%SmallDepth%DepthCutoff
            vMax       => setting%Limiter%Velocity%Maximum
        !%-----------------------------------------------------------------
        !% setting%coef is the blending adjustment (between 0.0 and 1.0)
        !% if coef == 1 then the V-shape element flowrate is replaced by 
        !% average of its faces. If coef < 1 then average value is blended
        !% with the present value of Q.

        !% Vcoef is the coefficient adjusted for local conditions
        Vcoef(thisP) = coef    

        !% find the cells that are deep enough to use the V filter
        Vvalue(thisP) = elemDepth(thisP) / (multiplier * smallDepth)
        where (Vvalue(thisP) > oneR)
            Vvalue(thisP) = oneR
        elsewhere
            Vvalue(thisP) = zeroR
            Vcoef(thisP)  = zeroR
        endwhere

        !% --- Reducing V-filter when Qlateral is large  20220524brh
        !%     HACK the fraction below should be replaced with a coefficient
        ! where (Qlateral(thisP) > onefourthR * abs(elemFlow(thisP)))
        !     Vcoef(thisP)  = Vcoef(thisP) * (onefourthR * abs(elemFlow(thisP)) / Qlateral(thisP))**2
        !     Vvalue(thisP) = zeroR     
        ! endwhere

        !% the Vvalue returns...
        !%  -1.0 if the element Q is between the face Q (not v-shaped)
        !%  +1.0 if the element Q is outside of the two face Q (v-shaped)
        !%  0.0 if the element Q is equal to one of the face Q (not v-shaped)
        !%  0.0 if the depth is too shallow
        Vvalue(thisP) =  (util_sign_with_ones_or_zero(faceFlow(mapUp(thisP)) - elemFlow(thisP)))      &
                        *(util_sign_with_ones_or_zero(faceFlow(mapDn(thisP)) - elemFlow(thisP)))      &
                        * Vvalue(thisP)     

        where (Vvalue(thisP) .le. zeroR)
            Vcoef(thisP) = zeroR
        endwhere

        !% blend the element and face-average flow rates
        elemFlow(thisP)  =  (oneR - Vcoef(thisP)) * elemFlow(thisP) &
                + Vcoef(thisP) * onehalfR * (faceflow(mapDn(thisP)) + faceflow(mapUp(thisP)))

        !% reset the velocity      
        elemVel(thisP) = elemFlow(thisP) / elemArea(thisP)   

        ! where (Vvalue(thisP) > zeroR)
        !     !% simple linear interpolation
        !     elemFlow(thisP)  =  (oneR - coef) * elemFlow(thisP) &
        !         + coef * onehalfR * (faceflow(mapDn(thisP)) + faceflow(mapUp(thisP)))
        !     !% reset the velocity      
        !     elemVel(thisP) = elemFlow(thisP) / elemArea(thisP)   
        ! endwhere
                
        !% reset for high velocity (typically due to small area)
        where (abs(elemVel(thisP)) > vMax)
            elemVel(thisP)  = sign( 0.99d0 * vMax, elemVel(thisP) )
            elemFlow(thisP) = elemVel(thisP) * elemArea(thisP)
        endwhere 
        
        !%------------------------------------------------------------------
        !% Closing
            !% clear the temporary Vvalue array
            Vvalue(thisP) = nullvalueR
            if (setting%Debug%File%adjust) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine adjust_Vshaped_flowrate
!%
!%==========================================================================  
!%==========================================================================  
!%
    subroutine adjust_Vshaped_head_surcharged (whichTM)
        !%------------------------------------------------------------------
        !% Description:
        !% Applies V-filter to surcharged head 
        !% HACK: ONLY APPLIES FOR PREISSMANN SLOT
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: whichTM
            integer, pointer :: thisCol, Npack
            integer, pointer :: thisP(:), mapUp(:), mapDn(:)
            real(8), pointer :: coef, multiplier, smallDepth
            real(8), pointer :: elemCrown(:), Vvalue(:), elemEllMax(:), Zbottom(:)
            real(8), pointer :: faceHeadUp(:), faceHeadDn(:), elemHead(:), elemVel(:)
            real(8), pointer :: w_uH(:), w_dH(:)
            logical, pointer :: isSlot(:)  !% Preissman Slot logical
            character(64) :: subroutine_name = 'adjust_Vshaped_head_surcharged'
        !%-------------------------------------------------------------------
        !% Preliminaries      
            if (setting%Debug%File%adjust) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% Aliases:
            select case (whichTM)
            case (ALLtm)
                !thisCol => col_elemP(ep_CC_ALLtm_ACsurcharged)
                print *, 'ALGORITHM DEVELOPMENT NEEDED FOR ALLtm with AC'
                call util_crashpoint(9587934)
            case (ETM)
                thisCol => col_elemP(ep_CC_Closed_Elements)
            case (AC)
                !thisCol => col_elemP(ep_CC_ACsurcharged)
                print *, 'ALGORITHM DEVLEOPMENT NEEDED FOR AC'
                call util_crashpoint(558723)
            case default
                print *, 'CODE ERROR: time march type unknown for # ', whichTM
                print *, 'which has key ',trim(reverseKey(whichTM))
                !stop 
                call util_crashpoint( 23943)
                return
            end select 

            !% coefficient for the blending adjustment (between 0.0 and 1.0)
            !% if coef == 1 then the V-shape element flowrate is replaced by 
            !% average of its faces.
            coef => setting%Adjust%Head%Coef
            
            if (coef .le. zeroR) return     
            Npack => npack_elemP(thisCol)
            if (Npack .le. 0) return

            thisP      => elemP(1:Npack,thisCol)
            mapUp      => elemI(:,ei_Mface_uL)
            mapDn      => elemI(:,ei_Mface_dL)    
            faceHeadUp => faceR(:,fr_Head_u)  
            faceHeadDn => faceR(:,fr_Head_d)          
            elemHead   => elemR(:,er_Head)    
            elemCrown  => elemR(:,er_Zcrown)
            elemEllMax => elemR(:,er_ell_max)
            w_uH       => elemR(:,er_InterpWeight_uH)
            w_dH       => elemR(:,er_InterpWeight_dH)
            Vvalue     => elemR(:,er_Temp01)
            Zbottom    => elemR(:,er_Zbottom)
            isSlot     => elemYN(:,eYN_isPSsurcharged)  !% Preissman slot

            multiplier => setting%Adjust%Head%FullDepthMultiplier

        !%-------------------------------------------------------------------
        !% --- For Preissman Slot, find the cells that are surcharged
        where (isSlot(thisP))
            Vvalue(thisP) = oneR
        elsewhere
            Vvalue(thisP) = zeroR 
        endwhere

        !% identify the V-shape locations
        Vvalue(thisP) =  (util_sign_with_ones_or_zero(faceHeadDn(mapUp(thisP)) - elemHead(thisP)))      &
                        *(util_sign_with_ones_or_zero(faceHeadUp(mapDn(thisP)) - elemHead(thisP)))      &
                        * Vvalue(thisP)   
                
        !% adjust where needed
        where (Vvalue(thisP) > zeroR)    
            !% simple linear interpolation
            elemHead(thisP)  =  (oneR - coef) * elemHead(thisP) &
               + coef * onehalfR * (faceHeadUp(mapDn(thisP)) + faceHeadDn(mapUp(thisP)))

        endwhere 

        !%-============================================================
        !% test 20220731
        !% --- for fully engaged preissmann slot -- use regular V-filter
        ! where ( isfSlot(mapUp(thisP)) .and. isfSlot(mapDn(thisP)) .and. (elemR(thisP,er_SlotDepth) > zeroR) )
        !     Vvalue(thisP) = oneR
        ! elsewhere
        !     Vvalue(thisP) = zeroR 
        ! endwhere

        ! !% identify the V-shape locations
        ! Vvalue(thisP) =  (util_sign_with_ones_or_zero(faceHeadDn(mapUp(thisP)) - elemHead(thisP)))      &
        !                 *(util_sign_with_ones_or_zero(faceHeadUp(mapDn(thisP)) - elemHead(thisP)))      &
        !                 * Vvalue(thisP)   
                
        ! !% adjust where needed
        ! where (Vvalue(thisP) > zeroR)    
        !     !% simple linear interpolation
        !     elemHead(thisP)  =  (oneR - coef) * elemHead(thisP) &
        !        + coef * onehalfR * (faceHeadUp(mapDn(thisP)) + faceHeadDn(mapUp(thisP)))
        ! endwhere   

        ! !% --- for cells where depth is below Zcrown, set the head to the minimum of neighbors
        ! where ( isfSlot(mapUp(thisP)) .and. isfSlot(mapDn(thisP)) .and. (elemR(thisP,er_SlotDepth) .le. zeroR) )
        !     Vvalue(thisP) = oneR
        ! elsewhere
        !     Vvalue(thisP) = zeroR 
        ! endwhere
        ! !% --- for these cells use the minimum of adjacent heads
        ! where (Vvalue(thisP) > zeroR)
        !     elemHead(thisP) = min(faceHeadUp(mapDn(thisP)),faceHeadDn(mapUp(thisP))) 
        ! end where

        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]" 
    end subroutine adjust_Vshaped_head_surcharged
!%    
!%==========================================================================  
!%==========================================================================  
!%
    subroutine adjust_faceflux_for_headgradient (thisP, thisDepthCutoff)
        !%------------------------------------------------------------------
        !% Description
        !% Adjusts the face flux for a head gradient when the depth of the
        !% face is more than twice the input thisDepthCutoff
        !% Should only be applied to faces of zero depth or small depth cells  
        !% thisP should be a packed set of element indexes that are either
        !% small depth or zero depth elements 
        !%------------------------------------------------------------------
        !% Declarations
            integer,  intent(in) :: thisP(:)
            real(8),  intent(in) :: thisDepthCutoff
            !logical, intent(in) ::  ifixQCons
            !integer, intent(in) :: whichTM
            integer, pointer :: fdn(:), fup(:)
            !integer, pointer :: thisColJM, thisJM(:), npackJM, isbranch(:) !% 20220122brh
            real(8), pointer :: faceQ(:), elemQ(:) !, fQCons(:)
            real(8), pointer :: faceHu(:), faceHd(:), faceAu(:), faceAd(:)
            real(8), pointer :: faceDu(:), faceDd(:)
            real(8), pointer :: elemH(:), elemL(:) !, elemVol(:)
            real(8), pointer :: dt, grav
            integer :: ii
        !%------------------------------------------------------------------
        !% Aliases:
            faceQ     => faceR(:,fr_Flowrate)
            faceHu    => faceR(:,fr_Head_u)
            faceHd    => faceR(:,fr_Head_d)
            faceAu    => faceR(:,fr_Area_u)
            faceAd    => faceR(:,fr_Area_d)
            faceDu    => faceR(:,fr_HydDepth_u)
            faceDd    => faceR(:,fr_HydDepth_d)
            !fQCons    => faceR(:,fr_Flowrate_Conservative)
            !elemQ     => elemR(:,er_Flowrate)
            elemH     => elemR(:,er_Head)
            elemL     => elemR(:,er_Length)
            
            !elemVol   => elemR(:,er_Volume_N0)
            !isbranch  => elemSI(:,esi_JunctionBranch_Exists)
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
                (faceDu(fdn(thisP)) > twoR * thisDepthCutoff) )
            faceQ(fdn(thisP)) = min(                                                 &
                faceQ(fdn(thisP)),                                                   &
                dt * grav * faceAu(fdn(thisP)) * (elemH(thisP) - faceHu(fdn(thisP))) &
                / (onehalfR * (elemL(thisP)))                                        &
                )
        end where

        !% --- for the upstream face dH/dx > 0 leads to a positive Q
        !%     Only applies where head gradient implies flow into the small volume and the
        !%     depth on the face is twice the small depth cutoff
        where ( (elemH(thisP) < faceHd(fup(thisP)) ) &
                .and. &
                (faceDd(fup(thisP)) > twoR * thisDepthCutoff) )
            faceQ(fup(thisP)) = max(                                                 &
                faceQ(fdn(thisP)),                                                   &
                dt * grav * faceAd(fup(thisP)) * (faceHd(fup(thisP)) - elemH(thisP)) &
                / (onehalfR * (elemL(thisP)))                                        &
                )
        end where

    end subroutine adjust_faceflux_for_headgradient
!%
!%==========================================================================
!%==========================================================================
!%
! subroutine adjust_smallvolumes_reset_old (Npack, thisCol)
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% nulls any prior storage of small volumes
    !     !%----------------------------------------------------------------------------- 
    !         integer, intent(in) :: Npack, thisCol  
    !         integer, pointer :: thisP(:)
    !     !%-----------------------------------------------------------------------------
    !         character(64) :: subroutine_name = 'adjust_smallvolumes_reset_old'
    !     !%-----------------------------------------------------------------------------
    !         if (Npack .le. 0) return
    !         !if (crashYN) return              
    !         if (setting%Debug%File%adjust) &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%-----------------------------------------------------------------------------
    !     !% Aliases       
    !         thisP => elemP(1:Npack,thisCol)
    !     !%----------------------------------------------------------------------------- 

    !     elemYN(thisP,eYN_isSmallDepth) = .false.
    !     elemR(thisP,er_SmallVolumeRatio) = nullvalueR

    !     !%-----------------------------------------------------------------------------
    !     if (setting%Debug%File%adjust) &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine  adjust_smallvolumes_reset_old
!%
!%==========================================================================
!%==========================================================================
!%
! subroutine adjust_smallvolumes_identify (Npack, thisCol, thisVolumeCol) 
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Identifies the small volumes out of a packed set
    !     !%------------------------------------------------------------------ 
    !     !% Declarations     
    !         integer, intent(in) :: Npack, thisCol, thisVolumeCol   
    !         integer, pointer :: thisP(:)
    !         real(8), pointer :: volume(:), smallvolume(:), svRatio(:)
    !         logical, pointer :: isSmallVol(:), isNearZeroDepth(:)
    !         character(64) :: subroutine_name = 'adjust_smallvolumes_identify'
    !     !%-----------------------------------------------------------------
    !     !% Preliminaries
    !         if (Npack .le. 0) return
    !         !if (crashYN) return              
    !         if (setting%Debug%File%adjust) &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%----------------------------------------------------------------- 
    !     !% Aliases
    !         thisP       => elemP(1:Npack,thisCol)
    !         volume      => elemR(:,thisVolumeCol)
    !         smallvolume => elemR(:,er_SmallVolume)
    !         svRatio     => elemR(:,er_SmallVolumeRatio)
    !         isSmallVol  => elemYN(:,eYN_isSmallDepth)
    !         isNearZeroDepth => elemYN(:,eYN_isZeroDepth)

    !     !%-------------------------------------------------------------------
 
    !     !% Find the small volume elements and set the SV ratio
    !     where (volume(thisP) < smallvolume(thisP))
    !         isSmallVol(thisP) = .true.
    !         svRatio(thisP) = volume(thisP) / smallvolume(thisP)
    !     endwhere
    
    !     where (isNearZeroDepth(thisP))
    !         !% for the elements that are near-zero let the zerovolume algorithm handle it
    !         isSmallVol(thisP) = .false.
    !     endwhere
        
    !     !%------------------------------------------------------------------
    !     if (setting%Debug%File%adjust) &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine adjust_smallvolumes_identify
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine adjust_smallvolumes_pack (Npack, thisColP, thisNewPackCol)
    !     !%  
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% packs an array for smallvolumes
    !     !%------------------------------------------------------------------
    !         integer, intent(in) :: Npack, thisColP, thisNewPackCol
    !         integer, pointer    :: thisP(:), eIdx(:)
    !         integer             :: newpack
    !         logical, pointer    :: isSmallVol(:)
    !         character(64) :: subroutine_name = 'adjust_smallvolumes_pack'
    !     !%------------------------------------------------------------------
    !         if (Npack .le. 0) return
    !         !if (crashYN) return              
    !         if (setting%Debug%File%adjust) &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%------------------------------------------------------------------
    !         thisP      => elemP(1:Npack,thisColP)
    !         eIdx       => elemI(:,ei_Lidx)
    !         isSmallVol => elemYN(:,eYN_isSmallDepth)
    !     !%------------------------------------------------------------------

    !     newpack = count(elemYN(thisP,eYN_isSmallDepth))
    !     npack_elemP(thisNewPackCol) = newpack  
    
    !     if (newpack > 0) then
    !         !% extract the set of small volumes.
    !         elemP(1:newpack,thisNewPackCol) = pack(eIdx(thisP), isSmallVol(thisP) )
    !     end if

    !     !% error checking
    !     if (any(elemR(elemP(1:newpack,thisNewPackCol) ,er_SmallVolume) == nullValueR)) then
    !         print *, 'FATAL ERROR -- CODE PROBLEM'
    !         print *, 'SmallVolume value that is set to nullvalueR.'
    !         print *, 'This is likely a failure to have one of the cross-sections shapes'
    !         print *, 'initialized in init_IC_set_SmallVolumes'
    !         print *, 'CODE FIX IS NEEDED/'
    !         stop 
    !           call util_crashpoint( 448795)
    !     end if

    !     !%------------------------------------------------------------------
    !     if (setting%Debug%File%adjust) &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine adjust_smallvolumes_pack
!%
!%==========================================================================
!%==========================================================================
!%
! subroutine adjust_velocity_smallvolume_blended &
    !     (Npack, thisCol)
    !     !%  
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% blends computed velocity with Chezy-Manning solution for small volumes 
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: Npack, thisCol
    !         integer, pointer :: thisP(:), fup(:), fdn(:)
    !         real(8), pointer :: fheadUp(:), fheadDn(:), length(:), area(:), HydRadius(:)
    !         real(8), pointer :: velocity(:), ManningsN(:), headslope(:), CMvelocity(:)
    !         real(8), pointer :: velocityBlend(:), svRatio(:), flowrate(:)
    !         character(64) :: subroutine_name = 'adjust_velocity_smallvolume_blended'
    !     !%------------------------------------------------------------------
    !     !% Preliminaries
    !         if (Npack .le. 0) return
    !         !if (crashYN) return              
    !         if (setting%Debug%File%adjust) &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%------------------------------------------------------------------
    !     !% Aliases
    !         thisP     => elemP(1:Npack,thisCol) !% only elements with small volumes
    !         fup       => elemI(:,ei_Mface_uL)
    !         fdn       => elemI(:,ei_Mface_dL)
    !         fheadUp   => faceR(:,fr_Head_u)
    !         fheadDn   => faceR(:,fr_Head_d)
    !         length    => elemR(:,er_Length)
    !         area      => elemR(:,er_Area)
    !         HydRadius => elemR(:,er_HydRadius)
    !         velocity  => elemR(:,er_Velocity)
    !         flowrate  => elemR(:,er_Flowrate)
    !         ManningsN => elemR(:,er_SmallVolume_ManningsN)
    !         headslope => elemR(:,er_SmallVolume_HeadSlope)
    !         CMvelocity => elemR(:,er_SmallVolume_CMvelocity)    
    !         svRatio    => elemR(:,er_SmallVolumeRatio)
    !         velocityBlend => elemR(:,er_Temp01)
    !     !%----------------------------------------------------------------------
    !     !% Adjust ManningsN for small volume CM velocity.
    !     !% Use the larger of the actual ManningsN or the setting% value

    !     ManningsN(thisP) = setting%SmallDepth%ManningsN
    !     where (ManningsN(thisP) < elemR(thisP,er_ManningsN))
    !         ManningsN(thisP) = elemR(thisP,er_ManningsN)
    !     endwhere

    !     !% slope of the piezometric head
    !     headslope(thisP) = (fheadDn(fup(thisP)) - fheadUp(fdn(thisP))) / length(thisP)

    !     !% absolute chezy-manning velocity based on slope
    !     CMvelocity(thisP) = ( HydRadius(thisP)**(twothirdR) ) * sqrt(abs(headslope(thisP))) / ManningsN(thisP)
                    
    !     !% assign direction to CM velocity.
    !     CMvelocity(thisP) = sign(CMvelocity(thisP), headslope(thisP))

    !     !% blend the computed velocity with CM velocity
    !     velocityBlend(thisP) = svRatio(thisP) * velocity(thisP) &
    !                         + (oneR - svRatio(thisP)) * CMvelocity(thisP)

    !     !% use the smaller velocity value, with the sign of the CM velocity
    !     velocity(thisP) = min(abs(CMvelocity(thisP)), velocity(thisP))
    !     velocity(thisP) = sign(velocity(thisP),CMvelocity(thisP))

    !     flowrate(thisP) = area(thisP) * velocity(thisP)

    !     !%----------------------------------------------------------------------
    !     !% Closing
    !     !% reset the temporary storage
    !         velocityBlend(thisP) = nullvalueR
    !         if (setting%Debug%File%adjust) &
    !             write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine adjust_velocity_smallvolume_blended
!%
!%==========================================================================
!%==========================================================================
!%
! subroutine adjust_zero_velocity_at_zero_depth &
    !     (Npack, thisCol)
    !     !%  
    !     !%----------------------------------------------------------------------
    !     !% Description:
    !     !% sets zeros to
    !     !%  velocity, flowrate, and head gradient across a zero-volume element
    !     !%----------------------------------------------------------------------  
    !         integer, intent(in) :: Npack, thisCol
    !         integer, pointer :: thisP(:), fup(:), fdn(:)
    !         real(8), pointer :: velocity(:), depth(:), flowrate(:), head(:)
    !         real(8), pointer :: fHeadUp(:), fHeadDn(:)
    !         character(64) :: subroutine_name = 'adjust_zero_velocity_at_zero_volume'
    !     !%----------------------------------------------------------------------
    !         if (Npack .le. 0) return
    !         !if (crashYN) return              
    !         if (setting%Debug%File%adjust) &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%--------------------------------------------------------------------
    !         thisP    => elemP(1:Npack,thisCol)
    !         depth    => elemR(:,er_Depth)
    !         velocity => elemR(:,er_Velocity)
    !         flowrate => elemR(:,er_Flowrate)
    !         head     => elemR(:,er_Head)
    !         fup      => elemI(:,ei_Mface_uL)
    !         fdn      => elemI(:,ei_Mface_dL)
    !         fHeadUp  => faceR(:,fr_Head_u)
    !         fHeadDn  => faceR(:,fr_Head_d) 
    !     !%------------------------------------------------------------------- 

    !     !% set the velocity and flowrate to zero and the head
    !     !% differenc across the element to zero.
    !     where (depth(thisP) .le. setting%ZeroValue%Depth)
    !         velocity(thisP) = zeroR
    !         flowrate(thisP) = zeroR  
    !         !% the head on the upstream side of the downstream face
    !         fHeadUp(fdn(thisP)) = head(thisP)
    !         !% the head on the downstream side of the usptream face
    !         fHeadDn(fup(thisP)) = head(thisP)
    !     endwhere    

    !     !%------------------------------------------------------------------- 
    !     if (setting%Debug%File%adjust) &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine adjust_zero_velocity_at_zero_depth
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module adjust