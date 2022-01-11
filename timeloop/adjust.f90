module adjust

    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use utility

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Makes ad hoc adjustments to ensure stability in hydraulic solution
    !%
    !% METHOD:
    !% 
    !%

    private

    public :: adjust_values

    public :: adjust_limit_by_zerovalues  !% used in geometry
    public :: adjust_limit_by_zerovalues_singular  !% used in geometry
    public :: adjust_limit_velocity_max



    !public :: adjust_smallvolume_face_values

    

    
    !public :: adjust_zerodepth_setvalues
    

    public :: adjust_zerodepth_identify_all
    public :: adjust_zerodepth_element_values 
    public :: adjust_zerodepth_face_fluxes_CC
    public :: adjust_zerodepth_face_fluxes_JMJB

    public :: adjust_smalldepth_identify_all
    public :: adjust_smalldepth_element_fluxes
    public :: adjust_smalldepth_face_fluxes
   

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine adjust_values (whichTM)
        !%------------------------------------------------------------------
        !% Description:
        !% Performs ad-hoc adjustments that may be needed for stability
        !%------------------------------------------------------------------
            integer, intent(in) :: whichTM  !% indicates which Time marching method
            character(64) :: subroutine_name = 'adjust_values'
        !%------------------------------------------------------------------
            if (icrash) return
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
                    print *, 'error, case default should not be reached'
                    print *, 'ad hoc flowrate adjust is .true., but approach is not supported'
                    stop 4973
                end select
        end if
        
        !% ad hoc adjustments to head
        if (setting%Adjust%Head%ApplyYN) then          
            select case (setting%Adjust%Head%Approach)
                case (vshape_surcharge_only)
                    call adjust_Vshaped_head_surcharged (whichTM)
                case default
                    print *, 'error, case default should not be reached'
                    print *, 'ad hoc head adjust is .true. but approach is not supported'
                    stop 9073
            end select
        end if
        
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine adjust_values
!%
!%========================================================================== 
!%==========================================================================  
!%        
    subroutine adjust_limit_by_zerovalues (geocol, geozero, thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% This applies a zero value limiter for a packed array that is not
        !% a priori limited to zero value elements
        !%-----------------------------------------------------------------------------
        !logical, intent(in) :: isreset
        integer, intent(in) :: geocol, thisCol
        real(8), intent(in) :: geozero
        integer, pointer :: Npack, thisP(:)
        real(8), pointer :: geovalue(:)    
        !logical, pointer :: NearZeroDepth(:),   isSmallDepth(:)   
        character(64) :: subroutine_name = 'adjust_limit_by_zerovalues'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        Npack        => npack_elemP(thisCol)  
        geovalue     => elemR(:,geocol)
        !isZeroDepth  => elemYN(:,eYN_isZeroDepth)
        !isSmallDepth => elemYN(:,eYN_isSmallDepth)
        !%-----------------------------------------------------------------------------

        if (Npack > 0) then
            thisP    => elemP(1:Npack,thisCol)
            where (geovalue(thisP) < geozero)
                    geovalue(thisP) = geozero
            endwhere
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
    subroutine adjust_limit_by_zerovalues_singular (eIdx, geocol, geozero)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Applies either the ZeroValue limiter (geozero) or zeroR as a lower limit to the
        !% geometry variable in elemR(:,geocol) for the single elemetn eIdx
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: geocol, eIdx
        real(8), intent(in) :: geozero
        real(8), pointer :: geovalue(:)        
        character(64) :: subroutine_name = 'adjust_limit_by_zerovalues_singular'
        if (icrash) return
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        geovalue => elemR(:,geocol)
        !%-----------------------------------------------------------------------------
        if (geovalue(eIdx) < geozero) then
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
            if (icrash) return
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
                    print *, 'error, default case should not be reached.'
                    stop 8368
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
            velocity(thisP) = sign( 0.99 * vMax, velocity(thisP) )
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
    !             stop 83655
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
    subroutine adjust_zerodepth_element_values (thisCol)
        !% -----------------------------------------------------------------
        !% Description:
        !% thisCol must be one of the ZeroDepth packed arrays that identifies
        !% all the (near) zero depth locations.
        !% -----------------------------------------------------------------
            integer, intent(in)  :: thisCol
            integer, pointer :: npack, thisP(:)
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
        elemR(thisP,er_HeadAvg) = setting%ZeroValue%Depth + elemR(thisP,er_Zbottom)
        

        !% only reset volume when it gets too small
        where (elemR(thisP,er_Volume) < setting%ZeroValue%Volume)
            elemR(thisP,er_Volume) = (oneR + onetenthR) * setting%ZeroValue%Volume 
        end where
        
        
    end subroutine adjust_zerodepth_element_values 
!%  
!%==========================================================================   
!%==========================================================================
!%
    subroutine adjust_smalldepth_element_fluxes()
        !% -----------------------------------------------------------------
        !% Description
        !% uses the ep_SmallDepth_CC_ALLtm pack to set the velocity and
        !% flowrate on all small volumes
        !% -----------------------------------------------------------------
        !% Declarations:
            integer, pointer :: thisCol
            integer, pointer :: npack, thisP(:), fdn(:), fup(:)
            real(8), pointer :: Area(:), BottomSlope(:), CMvelocity(:)
            real(8), pointer :: Flowrate(:), HydRadius(:), ManningsN(:)
            real(8), pointer :: Velocity(:), VelocityBlend(:), svRatio(:)
            real(8), pointer :: SmallVolume(:), Volume(:), faceFlow(:), faceFlowCons(:)
            integer, target  :: pset(2)
            real(8)          :: psign(2)
            integer :: ii
        !% -----------------------------------------------------------------
        !% Preliminaries:   
            !% local declarations to cycle through positive and negative slope sets
            pset(1) = ep_SmallDepth_CC_ALLtm_posSlope
            pset(2) = ep_SmallDepth_CC_ALLtm_negSlope
            psign(1) = oneR
            psign(2) = -oneR
        !% -----------------------------------------------------------------    
        !% Aliases    
            Area          => elemR(:,er_Area)
            BottomSlope   => elemR(:,er_BottomSlope)
            CMvelocity    => elemR(:,er_SmallVolume_CMvelocity) 
            Flowrate      => elemR(:,er_Flowrate)
            HydRadius     => elemR(:,er_HydRadius)
            ManningsN     => elemR(:,er_SmallVolume_ManningsN)   
            SmallVolume   => elemR(:,er_SmallVolume)  
            svRatio       => elemR(:,er_SmallVolumeRatio)
            Velocity      => elemR(:,er_Velocity)
            VelocityBlend => elemR(:,er_Temp01)
            Volume        => elemR(:,er_Volume)
        !% -----------------------------------------------------------------          
        do ii = 1,size(pset)
            thisCol => col_elemP(pset(ii))
            npack   => npack_elemP(thisCol)
            if (npack < 1) cycle
            thisP   => elemP(1:npack,thisCol)
            
            !% define the small volume ratio
            svRatio(thisP) = Volume(thisP) / SmallVolume(thisP)       
        
            !% use the larger of available roughness values
            ManningsN(thisP) = setting%SmallDepth%ManningsN
            ManningsN(thisP) = min(ManningsN(thisP), elemR(thisP,er_Roughness))
            ! where (ManningsN(thisP) < elemR(thisP,er_Roughness))
            !     ManningsN(thisP) = elemR(thisP,er_Roughness)
            ! endwhere    

            !% chezy-manning velocity based on bottom slope
            CMvelocity(thisP) = psign(ii) * ( HydRadius(thisP)**(twothirdR) ) &
                        * sqrt(abs(BottomSlope(thisP))) / ManningsN(thisP)
                    
            !% blend the computed velocity with CM velocity
            VelocityBlend(thisP) = svRatio(thisP) * velocity(thisP) &
                                + (oneR - svRatio(thisP)) * CMvelocity(thisP)

            !% use the smaller velocity value, with the sign of the CM velocity
            Velocity(thisP) = min(abs(CMvelocity(thisP)), velocity(thisP))
            Velocity(thisP) = sign(Velocity(thisP),CMvelocity(thisP))

            !% new flowrate
            Flowrate(thisP) = Area(thisP) * Velocity(thisP)

            !% reset the temporary storage
            VelocityBlend(thisP) = nullvalueR

        end do

        !% -----------------------------------------------------------------
    end subroutine adjust_smalldepth_element_fluxes  
!%  
!%========================================================================== 
!%==========================================================================
!%
    subroutine adjust_zerodepth_face_fluxes_CC (thisCol,ifixQCons)
        !% -----------------------------------------------------------------
        !% Description:
        !% thisCol must be one of the ZeroDepth packed arrays that identifies
        !% all the (near) zero depth element locations.
        !% only applicable to CC, not JM or JB
        !% if ifixQCons = .true. then the conservative fluxes are adjusted
        !% -----------------------------------------------------------------
            integer, intent(in)  :: thisCol
            logical, intent(in)  :: ifixQCons
            integer, pointer :: npack, thisP(:), fdn(:), fup(:)
            real(8), pointer :: fQ(:), fQCons(:)
        !% -----------------------------------------------------------------
        !% Aliases
            npack   => npack_elemP(thisCol)
            if (npack < 1) return
            thisP   => elemP(1:npack,thisCol)
            fdn     => elemI(:,ei_Mface_dL)
            fup     => elemI(:,ei_Mface_uL)
            fQ      => faceR(:,fr_Flowrate)
            fQCons  => faceR(:,fr_Flowrate_Conservative)
        !% -----------------------------------------------------------------
        !% choose either zero or an inflow
        fQ(fup(thisP)) = max(fQ(fup(thisP)), zeroR)
        fQ(fdn(thisP)) = min(fQ(fdn(thisP)), zeroR)

        !% reset the conservative fluxes
        if (ifixQCons) then
            fQCons(fup(thisP)) = fQ(fup(thisP))
            fQCons(fdn(thisP)) = fQ(fdn(thisP))
        end if

        !% -----------------------------------------------------------------
    end subroutine adjust_zerodepth_face_fluxes_CC        
!%  
!%==========================================================================   
!%==========================================================================
!%    
    subroutine adjust_zerodepth_face_fluxes_JMJB (thisCol, ifixQCons)
        !%------------------------------------------------------------------
        !% Description:
        !% Sets zero depth values on branches and JM. Input column must
        !% be a JM set
        !% if ifixQCons = .true. then the conservative fluxes are adjusted
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)  :: thisCol
            logical, intent(in)  :: ifixQCons
            integer, pointer :: npack, thisP(:), fup(:), fdn(:), isBranch(:)
            real(8), pointer :: fQ(:), fQCons(:)
            integer :: ii
        !%------------------------------------------------------------------
        !% Preliminaries:
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
    subroutine adjust_smalldepth_face_fluxes (ifixQCons)
        !%------------------------------------------------------------------
        !% Description:
        !% Sets the face values around an element where the ad-hoc
        !% small depth algorithm is used. Must be done after the small depth
        !% element values have been set.
        !% if ifixQCons = .true. then the conservative fluxes are adjusted
        !%------------------------------------------------------------------
        !% Declarations:
            logical, intent(in) :: ifixQCons
            integer, pointer :: fdn(:), fup(:), thisP(:), thisCol, npack
            real(8), pointer :: faceQ(:), elemQ(:), fQCons(:), slope(:)
            integer :: pset(2)
            integer :: ii
            logical :: istop = .false.
        !%------------------------------------------------------------------
        !% Preliminaries:
          ! % local declarations to cycle through positive and negative slope sets
            pset(1) = ep_SmallDepth_CC_ALLtm_posSlope
            pset(2) = ep_SmallDepth_CC_ALLtm_negSlope
        !%------------------------------------------------------------------
        !% Aliases:
            faceQ  => faceR(:,fr_Flowrate)
            fQCons => faceR(:,fr_Flowrate_Conservative)
            elemQ  => elemR(:,er_Flowrate)
            slope  => elemR(:,er_BottomSlope)
            fdn    => elemI(:,ei_Mface_dL)
            fup    => elemI(:,ei_Mface_uL)
        !%------------------------------------------------------------------
        do ii=1,size(pset)
            thisCol   => col_elemP(pset(ii))
            npack     => npack_elemP(thisCol)
            if (npack < 1) cycle 
            thisP     => elemP(1:npack,thisCol)

            !% set the downstream face as minimum of small volume flow or the
            !% (possibly negative) flow from the downstream element
            !% the upstream is zero or the inflow. Need where to handle
            !% elements with adverse slopes. Since this is a rare condition
            !% (and static, we can probably pack for this.)
            select case (ii)
            case (1)
                !% positive bottom slope
                faceQ(fdn(thisP)) = min(elemQ(thisP)     , faceQ(fdn(thisP)) )      
                faceQ(fup(thisP)) = max(faceQ(fup(thisP)), zeroR)
            case (2)
                !% negative bottom slope
                faceQ(fdn(thisP)) = min(faceQ(fdn(thisP)), zeroR )
                faceQ(fup(thisP)) = max(elemQ(thisP)     , faceQ(fup(thisP)) ) 
            case default
                print *, 'CODE ERROR -- default should not be reached'
                stop 3784848
            end select

            if (ifixQCons) then
                !% update the conservative face Q
                fQCons(fdn(thisP)) = faceQ(fdn(thisP))
                fQCons(fup(thisP)) = faceQ(fup(thisP))
            end if
       
            
        end do  
    
        !%------------------------------------------------------------------
        !% Closing:
    end subroutine adjust_smalldepth_face_fluxes

!%  
!%==========================================================================   
!%==========================================================================
!%
    ! subroutine adjust_zerodepth_element_values(thisCol)
    !     !% -----------------------------------------------------------------
    !     !% Description:
    !     !% thisCol must be one of the ZeroDepth packed arrays that identifies
    !     !% all the (near) zero depth locations.
    !     !% -----------------------------------------------------------------
    !         integer, intent(in)  :: thisCol
    !         integer, pointer :: npack, thisP(:)
    !     !% -----------------------------------------------------------------
    !     !% Aliases
    !         npack   => npack_elemP(thisCol)
    !         if (npack < 1) return
    !         thisP   => elemP(1:npack,thisCol)
    !     !% -----------------------------------------------------------------

    !     elemR(thisP,er_Velocity) = zeroR
    !     elemR(thisP,er_Flowrate) = zeroR


    ! end subroutine adjust_zerodepth_element_values
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine adjust_zerodepth_setvalues (thisCol)
    !     !%----------------------------------------------------------------------
    !     !% Description
    !     !% sets all the values in near zero-depth cells
    !     !%----------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: thisCol
    !         logical, pointer :: isZeroDepth(:)
    !         integer, pointer :: fdn(:), fup(:), thisP(:)
    !         integer          :: Npack
    !         integer :: ii
    !     !%----------------------------------------------------------------------
    !     !% Preliminaries
    !     !%----------------------------------------------------------------------
    !     !% Aliases
    !         isZeroDepth => elemYN(:,eYN_isZeroDepth)
    !         fdn         => elemI(:,ei_Mface_dL)
    !         fup         => elemI(:,ei_Mface_uL)
    !     !%----------------------------------------------------------------------

    !     !% --- set up for either the entire array or a packed section   
    !     if (thisCol == zeroI) then
    !         Npack = N_elem(this_image())
    !         thisP => elemI(1:Npack,ei_Lidx)
    !     else
    !         Npack = npack_elemP(thisCol)
    !         if (Npack > 0) then
    !             thisP => elemP(1:Npack,thisCol)
    !         else
    !             return
    !         end if
    !     end if
      
    !     !% HACK the following where statements were broken up as having too many
    !     !% statements between where's caused a segmentation fault for large systems.  
    !     where (isZeroDepth(thisP))
    !         elemR(thisP,er_Area)            = setting%ZeroValue%Area
    !     end where
    !     where (isZeroDepth(thisP))
    !         elemR(thisP,er_Depth)           = setting%ZeroValue%Depth
    !     end where
    !     where (isZeroDepth(thisP))
    !         elemR(thisP,er_dHdA)            = oneR / setting%ZeroValue%TopWidth
    !     end where
    !     where (isZeroDepth(thisP))
    !         elemR(thisP,er_ell)             = setting%ZeroValue%Depth
    !     end where
    !     where (isZeroDepth(thisP))
    !         elemR(thisP,er_Flowrate)        = zeroR
    !     end where
    !     where (isZeroDepth(thisP))
    !         elemR(thisP,er_FroudeNumber)    = zeroR
    !     end where
    !     where (isZeroDepth(thisP))
    !         elemR(thisP,er_HydDepth)        = setting%ZeroValue%Depth
    !     end where
    !     where (isZeroDepth(thisP))
    !         elemR(thisP,er_Perimeter)       = setting%ZeroValue%TopWidth + setting%ZeroValue%Depth
    !     end where
    !     where (isZeroDepth(thisP))
    !         elemR(thisP,er_HydRadius)       = setting%ZeroValue%Area  / (setting%ZeroValue%TopWidth + setting%ZeroValue%Depth)
    !     end where
    !     where (isZeroDepth(thisP))            
    !         elemR(thisP,er_Topwidth)        = setting%ZeroValue%TopWidth
    !     end where
    !     where (isZeroDepth(thisP))            
    !         elemR(thisP,er_Velocity)        = zeroR
    !     end where
    !     where (isZeroDepth(thisP))            
    !         elemR(thisP,er_WaveSpeed)       = zeroR
    !     end where
    !     where (isZeroDepth(thisP))            
    !         !% zerodepth forces the interpolation to use the neighbor, but later
    !         !% sets the fluxes so that only inward fluxes apply (see adjust_zerodepth_facevalues)
    !         !% Note that weights_uH and dH are unchanged on purpose!
    !         elemR(thisP,er_InterpWeight_uQ) = setting%Limiter%InterpWeight%Maximum
    !     end where
    !     where (isZeroDepth(thisP))           
    !         elemR(thisP,er_InterpWeight_dQ) = setting%Limiter%InterpWeight%Maximum
    !     end where
    !     where (isZeroDepth(thisP))            
    !         elemR(thisP,er_InterpWeight_uG) = setting%Limiter%InterpWeight%Maximum
    !     end where
    !     where (isZeroDepth(thisP))            
    !         elemR(thisP,er_InterpWeight_dG) = setting%Limiter%InterpWeight%Maximum
    !     end where
    !     where (isZeroDepth(thisP))
    !         elemR(thisP,er_Head)            = setting%ZeroValue%Depth + elemR(thisP,er_Zbottom)
    !     end where
    !     where (isZeroDepth(thisP))
    !         elemR(thisP,er_HeadAvg)         = setting%ZeroValue%Depth + elemR(thisP,er_Zbottom)
    !     end where

    !     !% don't reset the volume unless it gets smaller than 10% of zero value
    !     where (elemR(thisP,er_Volume) < setting%ZeroValue%Volume / tenR)
    !         elemR(thisP,er_Volume) = setting%ZeroValue%Volume
    !     end where

    !     call adjust_zerodepth_facevalues (thisCol)

    !     !%----------------------------------------------------------------------
    !     !% Closing    


    ! end subroutine adjust_zerodepth_setvalues    
!%     
!%==========================================================================
!%==========================================================================
!%   
    ! subroutine adjust_zerodepth_face_values (thisCol)
    !     !%----------------------------------------------------------------------
    !     !% Description
    !     !% sets the face fluxes to zero for outflows from a zerodepth element
    !     !% call with thisCol=0 to do all elements
    !     !%----------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: thisCol
    !         logical, pointer :: isZeroDepth(:)
    !         integer, pointer :: fdn(:), fup(:), thisP(:)
    !         integer          :: Npack
    !         integer :: ii
    !     !%----------------------------------------------------------------------
    !     !% Preliminaries
    !     !%----------------------------------------------------------------------
    !     !% Aliases
    !         isZeroDepth => elemYN(:,eYN_isZeroDepth)
    !         fdn => elemI(:,ei_Mface_dL)
    !         fup => elemI(:,ei_Mface_uL)
    !     !%----------------------------------------------------------------------

    !     !% --- set up for either the entire array or a packed section   
    !     if (thisCol == zeroI) then
    !         Npack = N_elem(this_image())
    !         thisP => elemI(1:Npack,ei_Lidx)
    !     else
    !         Npack = npack_elemP(thisCol)
    !         if (Npack > 0) then
    !             thisP => elemP(1:Npack,thisCol)
    !         else
    !             return
    !         end if
    !     end if

    !     !% On the upstream face set the downstream element values for area and head
    !     where ((isZeroDepth(thisP)) .and. (fup(thisP) .ne. nullValueI))
    !         faceR(fup(thisP),fr_Area_d)        = setting%ZeroValue%Area
    !         faceR(fup(thisP),fr_Head_d)        = setting%ZeroValue%Depth + elemR(thisP,er_Zbottom)
    !     end where

    !     !% On the downstream face set the upstream element values for area and head
    !     where ((isZeroDepth(thisP)) .and. (fdn(thisP) .ne. nullValueI))
    !         faceR(fdn(thisP),fr_Area_u)        = setting%ZeroValue%Area
    !         faceR(fdn(thisp),fr_Head_u)        = setting%ZeroValue%Depth + elemR(thisP,er_Zbottom)
    !     end where
            
    !     !% a zero volume can only have an inflow from faces (no outflow)...
    !     where ( (isZeroDepth(thisP)) .and. (fdn(thisP) .ne. nullvalueI) )
    !         !% set the downstream face flowrate to zero or keep the value if it is negative
    !         faceR(fdn(thisP),fr_Flowrate) &
    !             = onehalfR * ( faceR(fdn(thisP),fr_Flowrate) - abs(faceR(fdn(thisP),fr_Flowrate)) )
    !     end where

    !     where ( (isZeroDepth(thisP)) .and. (fup(thisP) .ne. nullvalueI) )
    !         !% set the upstream face flowrate to zero or keep the value if it is positive    
    !         faceR(fup(thisP),fr_Flowrate) &
    !             = onehalfR * ( faceR(fup(thisP),fr_Flowrate) + abs(faceR(fup(thisP),fr_Flowrate)) )
    !     end where

    ! end subroutine adjust_zerodepth_face_values
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
            real(8), pointer :: coef, vMax
            real(8), pointer :: faceFlow(:), elemFlow(:), elemVel(:)
            real(8), pointer ::  w_uQ(:), w_dQ(:), elemArea(:), Vvalue(:)
            logical, pointer :: isSmallDepth(:), isNearZeroDepth(:)
            character(64) :: subroutine_name = 'adjust_Vshaped_flowrate'
        !%------------------------------------------------------------------
        !% Preliminaries    
            if (setting%Adjust%Flowrate%Coef .le. zeroR) return
            if (icrash) return       
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
                print *, 'error, this default case should not be reached'
                stop 9239
            end select

            Npack => npack_elemP(thisCol)
            if (Npack .le. 0) return
            
            thisP    => elemP(1:Npack,thisCol)
            mapUp    => elemI(:,ei_Mface_uL)
            mapDn    => elemI(:,ei_Mface_dL)   
            faceFlow => faceR(:,fr_Flowrate)  
            elemFlow => elemR(:,er_Flowrate)    
            elemVel  => elemR(:,er_Velocity)
            elemArea => elemR(:,er_Area)
            w_uQ     => elemR(:,er_InterpWeight_uQ)
            w_dQ     => elemR(:,er_InterpWeight_dQ)
            Vvalue   => elemR(:,er_Temp01)
            isSmallDepth   => elemYN(:,eYN_isSmallDepth)
            isNearZeroDepth => elemYN(:,eYN_isZeroDepth)
            coef => setting%Adjust%Flowrate%Coef
            vMax => setting%Limiter%Velocity%Maximum
        !%-----------------------------------------------------------------
        !% coef is the blending adjustment (between 0.0 and 1.0)
        !% if coef == 1 then the V-shape element flowrate is replaced by 
        !% average of its faces. If coef < 1 then average value is blended
        !% with the present value of Q.

        !print *, 'made it here in adjust velocity'
    
        !% the Vvalue returns...
        !%  -1.0 if the element Q is between the face Q (not v-shaped)
        !%  +1.0 if the element Q is outside of the two face Q (v-shaped)
        !%  0.0 if the element Q is equal to one of the face Q (not v-shaped)
        Vvalue(thisP) =  (util_sign_with_ones_or_zero(faceFlow(mapUp(thisP)) - elemFlow(thisP)))      &
                        *(util_sign_with_ones_or_zero(faceFlow(mapDn(thisP)) - elemFlow(thisP)))    

        ! if (this_image() == debug_image) then
        !     print *, ' in adjust ',Vvalue(ietmp(2)), setting%Time%Now
        !     print *, faceFlow(mapUp(ietmp(2))), elemFlow(ietmp(2)), faceFlow(mapDn(ietmp(2)))
        ! end if

        where     ( (Vvalue(thisP) > zeroR)        &
            .and.   (.not. isSmallDepth  (thisP))   &
            .and.   (.not. isNearZeroDepth(thisP))  )

            !% averaging based on interpolation weights
            !% this had problems with lateral inflow conditions that
            ! !% change the interpolation weights -- brh 20220207
            ! elemFlow(thisP) =  (oneR - coef) * elemFlow(thisP) &
            !     + coef *                                       &
            !         (  w_uQ(thisP) * faceflow(mapDn(thisP))    &
            !          + w_dQ(thisP) * faceflow(mapUp(thisP)) )  &
            !     / ( w_uQ(thisP) + w_dQ(thisP) )

            !% simple linear interpolation
            elemFlow(thisP)  =  (oneR - coef) * elemFlow(thisP) &
                + coef * onehalfR * (faceflow(mapDn(thisP)) + faceflow(mapUp(thisP)))

            !% reset the velocity      
            elemVel(thisP) = elemFlow(thisP) / elemArea(thisP)   
        endwhere

        ! if (this_image() == debug_image) then
        !     print *, faceFlow(mapUp(ietmp(2))), elemFlow(ietmp(2)), faceFlow(mapDn(ietmp(2)))
        ! end if
                
        where (abs(elemVel(thisP)) > vMax)
            elemVel(thisP) = sign( 0.99 * vMax, elemVel(thisP) )
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
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: whichTM
        integer, pointer :: thisCol, Npack
        integer, pointer :: thisP(:), mapUp(:), mapDn(:)
        real(8), pointer :: coef
        real(8), pointer :: faceHeadUp(:), faceHeadDn(:), elemHead(:), elemVel(:)
        real(8), pointer :: w_uH(:), w_dH(:)
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'adjust_Vshaped_head_surcharged'
        !%-----------------------------------------------------------------------------
        if (icrash) return              
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        select case (whichTM)
            case (ALLtm)
                thisCol => col_elemP(ep_CC_ALLtm_surcharged)
            case (ETM)
                thisCol => col_elemP(ep_CC_ETM_surcharged)
            case (AC)
                thisCol => col_elemP(ep_CC_AC_surcharged)
            case default
                print *, 'error, this default case should not be reached'
                stop 2394
        end select 

        !% HACK: below is the old code. For similar reason mentioned above
        !% the code has been modified
        ! select case (whichTM)
        !     case (ALLtm)
        !         thisCol => col_elemP(ep_CCJB_ALLtm_surcharged)
        !     case (ETM)
        !         thisCol => col_elemP(ep_CCJB_ETM_surcharged)
        !     case (AC)
        !         thisCol => col_elemP(ep_CCJB_AC_surcharged)
        !     case default
        !         print *, 'error, this default case should not be reached'
        !         stop 2394
        ! end select   
    
        !% coefficient for the blending adjustment (between 0.0 and 1.0)
        !% if coef == 1 then the V-shape element flowrate is replaced by 
        !% average of its faces.
        coef => setting%Adjust%Head%Coef
        
        if (coef > zeroR) then       
            Npack => npack_elemP(thisCol)
            if (Npack > 0) then
                thisP      => elemP(1:Npack,thisCol)
                mapUp      => elemI(:,ei_Mface_uL)
                mapDn      => elemI(:,ei_Mface_dL)    
                faceHeadUp => faceR(:,fr_Head_u)  
                faceHeadDn => faceR(:,fr_Head_d)          
                elemHead   => elemR(:,er_Head)    
                w_uH       => elemR(:,er_InterpWeight_uH)
                w_dH       => elemR(:,er_InterpWeight_dH)
                
                !% identify the V-shape condition
                where  ( (util_sign_with_ones(faceHeadDn(mapUp(thisP)) - elemHead(thisP)))      &
                        *(util_sign_with_ones(faceHeadUp(mapDn(thisP)) - elemHead(thisP))) > 0)
                    
                    !% averaging based on interpolation weights
                    elemHead(thisP) =  (oneR - coef) * elemHead(thisP)  &
                        + coef *                                        &
                            (  w_uH(thisP) * faceHeadUp(mapDn(thisP))   &
                             + w_dH(thisP) * faceHeadDn(mapUp(thisP)) ) &
                        / ( w_uH(thisP) + w_dH(thisP) )
                        
                endwhere                       
            end if
        end if
        
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]" 
    end subroutine
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
    !         if (icrash) return              
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
    !         if (icrash) return              
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
    !         if (icrash) return              
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
    !         stop 448795
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
    !         if (icrash) return              
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
    !     !% Use the larger of the actual roughness or the setting% value

    !     ManningsN(thisP) = setting%SmallDepth%ManningsN
    !     where (ManningsN(thisP) < elemR(thisP,er_Roughness))
    !         ManningsN(thisP) = elemR(thisP,er_Roughness)
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
    !         if (icrash) return              
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