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
    public :: adjust_limit_by_zerovalues
    public :: adjust_limit_by_zerovalues_singular
    public :: adjust_velocity
    public :: adjust_face_dynamic_limit
   

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine adjust_values (whichTM)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs ad-hoc adjustments that may be needed for stability
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: whichTM  !% indicates which Time marching method
        character(64) :: subroutine_name = 'adjust_values'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
              
        !% ad hoc adjustments to flowrate 
        if (setting%Adjust%Flowrate%Apply) then   
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
        if (setting%Adjust%Head%Apply) then          
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
        !% Applies either the ZeroValue limiter (geozero) or zeroR as a lower limit to the
        !% geometry variable in elemR(:,geocol)
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: geocol, thisCol
        real(8), intent(in) :: geozero
        integer, pointer :: Npack, thisP(:)
        real(8), pointer :: geovalue(:)    
        logical, pointer :: NearZeroVolume(:)   
        character(64) :: subroutine_name = 'adjust_limit_by_zerovalues'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        Npack    => npack_elemP(thisCol)  
        geovalue => elemR(:,geocol)
        NearZeroVolume => elemYN(:,eYN_isNearZeroVolume)
        !%-----------------------------------------------------------------------------

        if (Npack > 0) then
            thisP    => elemP(1:Npack,thisCol)
            !% reset the NearZeroVolumes
            NearZeroVolume(thisP) = .false.
            
            if (setting%ZeroValue%UseZeroValues) then
                where (geovalue(thisP) < geozero)
                    geovalue(thisP) = geozero
                    NearZeroVolume(thisP) = .true.
                endwhere
            else
                where (geovalue(thisP) < zeroR)
                    geovalue(thisP) = zeroR
                    NearZeroVolume(thisP) = .true.
                endwhere  
            end if
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
        if (setting%ZeroValue%UseZeroValues) then
            if (geovalue(eIdx) < geozero) then
                geovalue(eIdx) = geozero
            end if
        else
            if (geovalue(eIdx) < zeroR) then
                geovalue(eIdx) = zeroR
            end if 
        end if 

        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine adjust_limit_by_zerovalues_singular
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine adjust_velocity (whichTM, velocityCol, volumeCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% employs velocity limiters and small volume treatments to limit 
        !% destabilizing velocities.
        !% The "velocityCol" is the velocity column in elemR that is adjusted
        !% The "volumeCol" is the volume column in elemR that is adjusted
        !%-----------------------------------------------------------------------------   
        integer, intent(in) ::whichTM, velocityCol, volumeCol
        integer, pointer :: thisCol_all, thisSmallVolumeCol, thisVelocityCol
        integer, pointer :: Npack
        character(64) :: subroutine_name = 'adjust_velocity'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        select case (whichTM)
            case (ALLtm)
                !% HACK: small velocity adjustment should only be done for CC elements
                !% since the velocity is solved there
                thisCol_all        => col_elemP(ep_CC_ALLtm)
                ! thisCol_all        => col_elemP(ep_ALLtm)
                thisSmallVolumeCol => col_elemP(ep_smallvolume_ALLtm)
                thisVelocityCol    => col_elemR(er_Velocity)
            case (ETM)
                !% HACK: small velocity adjustment should only be done for CC elements
                !% since the velocity is solved there
                thisCol_all        => col_elemP(ep_CC_ETM)
                ! thisCol_all        => col_elemP(ep_ETM)
                thisSmallVolumeCol => col_elemP(ep_smallvolume_ETM)
                thisVelocityCol    => col_elemR(er_Velocity)        
            case (AC)
                !% HACK: small velocity adjustment should only be done for CC elements
                !% since the velocity is solved there
                thisCol_all        => col_elemP(ep_CC_AC)
                ! thisCol_all        => col_elemP(ep_AC)
                thisSmallVolumeCol => col_elemP(ep_smallvolume_AC)    
                thisVelocityCol    => col_elemR(er_Velocity)        
            case default
                print *, 'error, default case should not be reached.'
                stop 8368
        end select
        
        !% reset the small volumes for flow/velocity limit computations
        if (setting%SmallVolume%UseSmallVolumes) then
            Npack => npack_elemP(thisCol_all) 
            if (Npack > 0) then   
                call adjust_smallvolumes_reset_old (Npack,thisCol_all)  
                call adjust_smallvolumes_identify (Npack, thisCol_all, volumeCol)        
                call adjust_smallvolumes_pack (Npack, thisCol_all, thisSmallVolumeCol)        
            end if
        end if
        
        !% apply ad-hoc velocity limiter
        if (setting%Limiter%Velocity%UseLimitMax) then
            Npack => npack_elemP(thiscol_all)
            if (Npack > 0) then 
                call adjust_velocity_limiter_reset_old (Npack, thiscol_all) 
                call adjust_velocity_limiter (Npack, thiscol_all, velocityCol) 
            end if
        end if
        
        !% For small volumes, compute a velocity that is blended from
        !% the update value and a Chezy-Manning computed using the 
        !% free surface slope of the element
        ! if (setting%SmallVolume%UseSmallVolumes) then
        !     Npack => npack_elemP(thisSmallVolumeCol)
        !     if (Npack > 0) then
        !         call adjust_velocity_smallvolume_blended    &
        !             (Npack, thisSmallVolumeCol, velocityCol)
        !     end if
        ! end if
        
        !% for extremely small volumes set velocity to zero
        if (setting%ZeroValue%UseZeroValues) then
             call adjust_zero_velocity_at_zero_volume    &
                (Npack, thiscol_all, velocityCol, volumeCol)
        end if
        
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine adjust_velocity
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine adjust_face_dynamic_limit (facePackCol, isInterior)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% This subroutine calculates the face valocity and adjusts for limiter
        !%-----------------------------------------------------------------------------   
        integer, intent(in) :: facePackCol
        logical, intent(in) :: isInterior
        integer, pointer :: Npack, thisP(:)
        real(8), pointer :: f_area_u(:), f_area_d(:), f_velocity_u(:), f_velocity_d(:)
        real(8), pointer :: f_flowrate(:), zeroValue, vMax
        character(64) :: subroutine_name = 'adjust_face_dynamic_limit'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        f_area_u     => faceR(:,fr_Area_u)
        f_area_d     => faceR(:,fr_Area_d)
        f_velocity_u => faceR(:,fr_Velocity_u)
        f_velocity_d => faceR(:,fr_Velocity_d)
        f_flowrate   => faceR(:,fr_Flowrate)
        zeroValue    => setting%ZeroValue%Area
        vMax         => setting%Limiter%Velocity%Maximum
        !%-----------------------------------------------------------------------------
        if (isInterior) then
            !% face velocity calculation at the interior faces
            Npack => npack_faceP(facePackCol)
            thisP => faceP(1:Npack,facePackCol)
        else
            !% face velocity calculation at the shared faces
            Npack => npack_facePS(facePackCol)
            thisP => facePS(1:Npack,facePackCol)
        end if

        if (Npack > 0) then
            if (setting%ZeroValue%UseZeroValues) then
                !% ensure face area_u is not smaller than zerovalue
                where (f_area_u(thisP) < zeroValue)
                    f_area_u(thisP) = zeroValue
                endwhere
                !% ensure face area_d is not smaller than zerovalue
                where (f_area_d(thisP) < zeroValue)
                    f_area_d(thisP) = zeroValue
                endwhere

                where (f_area_u(thisP) >= zeroValue)
                    f_velocity_u(thisP) = f_flowrate(thisP)/f_area_u(thisP)
                endwhere

                where (f_area_d(thisP) >= zeroValue)
                    f_velocity_d(thisP) = f_flowrate(thisP)/f_area_d(thisP)
                endwhere
            else
                where (f_area_u(thisP) < zeroR)
                    f_area_u(thisP) = zeroR
                endwhere 

                where (f_area_d(thisP) < zeroValue)
                    f_area_d(thisP) = zeroR
                endwhere 

                where (f_area_u(thisP) >= zeroR)
                    f_velocity_u(thisP) = f_flowrate(thisP)/f_area_u(thisP)
                endwhere

                where (f_area_d(thisP) >= zeroR)
                    f_velocity_d(thisP) = f_flowrate(thisP)/f_area_d(thisP)
                endwhere
            end if

            !%  limit high velocities
            if (setting%Limiter%Velocity%UseLimitMax) then
                where(abs(f_velocity_u(thisP))  > vMax)
                    f_velocity_u(thisP) = sign(0.99 * vMax, f_velocity_u(thisP))
                endwhere

                where(abs(f_velocity_d(thisP))  > vMax)
                    f_velocity_d(thisP) = sign(0.99 * vMax, f_velocity_d(thisP))
                endwhere
            end if
        end if
        
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine adjust_face_dynamic_limit
    !%
    !%==========================================================================
    !% PRIVATE
    !%==========================================================================   
    !%  
    subroutine adjust_Vshaped_flowrate (whichTM)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs 
        !%-----------------------------------------------------------------------------    
        integer, intent(in) :: whichTM
        integer, pointer :: thisCol, Npack
        integer, pointer :: thisP(:), mapUp(:), mapDn(:)
        real(8), pointer :: coef, vMax
        real(8), pointer :: faceFlow(:), elemFlow(:), elemVel(:), w_uQ(:), w_dQ(:), elemArea(:)
        logical, pointer :: isAdhocFlowrate(:)
        !%-----------------------------------------------------------------------------  
        character(64) :: subroutine_name = 'adjust_Vshaped_flowrate'
        !%-----------------------------------------------------------------------------
        if (icrash) return       
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------  
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
        
        !% HACK: below is the original code.
        !% I think the v-shape adjustment should only be called for CC elements
        !% because flowrate is only solved for CC elements
        !% calling the adjustment for JB elements will result in segmentation fault
        !% because JB elements will a ha nullvale for upstream/downstream face map 
        !% since if a JB will always connected to a JM which does not has any faces.

        ! select case (whichTM)
        !     case (ALLtm)
        !         thisCol => col_elemP(ep_CCJB_ALLtm)
        !     case (ETM)
        !         thisCol => col_elemP(ep_CCJB_ETM)
        !     case (AC)
        !         thisCol => col_elemP(ep_CCJB_AC)
        !     case default
        !         print *, 'error, this default case should not be reached'
        !         stop 9239
        ! end select
    
        !% coefficient for the blending adjustment (between 0.0 and 1.0)
        !% if coef == 1 then the V-shape element flowrate is replaced by the weighted
        !% average of its faces.
        coef => setting%Adjust%Flowrate%Coef
        vMax => setting%Limiter%Velocity%Maximum
        if (coef > zeroR) then      
            Npack => npack_elemP(thisCol)

            if (Npack > 0) then
                thisP    => elemP(1:Npack,thisCol)
                mapUp    => elemI(:,ei_Mface_uL)
                mapDn    => elemI(:,ei_Mface_dL)   
                faceFlow => faceR(:,fr_Flowrate)  
                elemFlow => elemR(:,er_Flowrate)    
                elemVel  => elemR(:,er_Velocity)
                elemArea => elemR(:,er_Area)
                w_uQ     => elemR(:,er_InterpWeight_uQ)
                w_dQ     => elemR(:,er_InterpWeight_dQ)
                isAdhocFlowrate => elemYN(:,eYN_IsAdhocFlowrate)


                !% identify the V-shape condition
                where  ( (util_sign_with_ones(faceFlow(mapUp(thisP)) - elemFlow(thisP)))      &
                        *(util_sign_with_ones(faceFlow(mapDn(thisP)) - elemFlow(thisP))) > 0)

                    !% averaging based on interpolation weights
                    elemFlow(thisP) =  (oneR - coef) * elemFlow(thisP) &
                        + coef *                                       &
                            (  w_uQ(thisP) * faceflow(mapDn(thisP))    &
                             + w_dQ(thisP) * faceflow(mapUp(thisP)) )  &
                        / ( w_uQ(thisP) + w_dQ(thisP) )

                    !% reset the velocity      
                    elemVel(thisP) = elemFlow(thisP) / elemArea(thisP)   
                endwhere
                
                where (abs(elemVel(thisP)) > vMax)
                    elemVel(thisP) = sign( 0.99 * vMax, elemVel(thisP) )
                    isAdhocFlowrate(thisP) = .true.
                endwhere 
            end if
        end if
        
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
        !% if coef == 1 then the V-shape element flowrate is replaced by the weighted
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
    subroutine adjust_smallvolumes_reset_old (Npack, thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% nulls any prior storage of small volumes
        !%----------------------------------------------------------------------------- 
        integer, intent(in) :: Npack, thisCol  
        integer, pointer :: thisP(:)
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'adjust_smallvolumes_reset_old'
        !%-----------------------------------------------------------------------------
        if (icrash) return              
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        thisP => elemP(1:Npack,thisCol)
        !%----------------------------------------------------------------------------- 
        elemYN(thisP,eYN_IsSmallVolume) = .false.
        elemR(thisP,er_SmallVolumeRatio) = nullvalueR
        
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine  adjust_smallvolumes_reset_old
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine adjust_smallvolumes_identify (Npack, thisCol, thisVolumeCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Identifies the small volumes
        !%-----------------------------------------------------------------------------        
        integer, intent(in) :: Npack, thisCol, thisVolumeCol   
        integer, pointer :: thisP(:)
        real(8), pointer :: volume(:), smallvolume(:), svRatio(:)
        logical, pointer :: isSmallVol(:)
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'adjust_smallvolumes_identify'
        !%-----------------------------------------------------------------------------
        if (icrash) return              
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------  
        thisP       => elemP(1:Npack,thisCol)
        volume      => elemR(:,thisVolumeCol)
        smallvolume => elemR(:,er_SmallVolume)
        svRatio     => elemR(:,er_SmallVolumeRatio)
        isSmallVol  => elemYN(:,eYN_isSmallVolume)
        !%----------------------------------------------------------------------------- 
 
        !% Find the small volume elements and set the SV ratio
        where (volume(thisP) < smallvolume(thisP))
            isSmallVol(thisP) = .true.
            svRatio(thisP) = volume(thisP) / smallvolume(thisP)
        endwhere
        !% for the elements that are near-zero, set the SV ratio to zero, This ensures only Chezy-Manning is used for solution
        where (volume(thisP) .le. setting%Zerovalue%Volume )
            sVratio(thisP) = zeroR
        endwhere
        
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine adjust_smallvolumes_identify
    !%  
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine adjust_smallvolumes_pack (Npack, thisColP, thisNewPackCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% packs an array for smallvolumes
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: Npack, thisColP, thisNewPackCol
        integer, pointer :: thisP(:), eIdx(:)
        integer :: newpack
        logical, pointer :: isSmallVol(:)
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'adjust_smallvolumes_pack'
        !%-----------------------------------------------------------------------------
        if (icrash) return              
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        thisP      => elemP(1:Npack,thisColP)
        eIdx       => elemI(:,ei_Lidx)
        isSmallVol => elemYN(:,eYN_isSmallVolume)
        !%-----------------------------------------------------------------------------
        newpack = count(elemYN(thisP,eYN_issmallvolume))
        npack_elemP(thisNewPackCol) = newpack  
    
        if (newpack > 0) then
            !% extract the set of small volumes.
            elemP(1:newpack,thisNewPackCol) = pack(eIdx(thisP), isSmallVol(thisP) )
        end if

        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine adjust_smallvolumes_pack
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine adjust_velocity_limiter_reset_old (Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
         !% removes ad-hoc flowrate designation
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: Npack, thisCol
        integer, pointer :: thisP(:)
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'adjust_velocity_limiter_reset_old'
        !%-----------------------------------------------------------------------------
        if (icrash) return              
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        thisP => elemP(1:Npack,thisCol)
        !%-----------------------------------------------------------------------------

        elemYN(thisP,eYN_IsAdhocFlowrate) = .false.
    
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine adjust_velocity_limiter_reset_old
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine adjust_velocity_limiter (Npack, thisPackCol, thisVelocityCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Ad hoc limit of velocity to 99% of maximum allowed
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: Npack, thisPackCol, thisVelocityCol
        integer, pointer :: thisP(:)
        real(8), pointer :: velocity(:), vMax
        logical, pointer :: isAdhocFlowrate(:)
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'adjust_velocity_limiter'
        !%-----------------------------------------------------------------------------
        if (icrash) return              
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        thisP           => elemP(1:Npack,thisPackCol)
        velocity        => elemR(:,thisVelocityCol)
        isAdhocFlowrate => elemYN(:,eYN_IsAdhocFlowrate)
        vMax            => setting%Limiter%Velocity%Maximum
        !%-----------------------------------------------------------------------------

        where (abs(velocity(thisP)) > vMax)
            velocity(thisP) = sign( 0.99 * vMax, velocity(thisP) )
            isAdhocFlowrate(thisP) = .true.
        endwhere 

        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine adjust_velocity_limiter
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine adjust_velocity_smallvolume_blended (Npack, thisCol, thisVelocityCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% blends computed velocity with Chezy-Manning solution for small volumes 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: Npack, thisCol, thisVelocityCol
        integer, pointer :: thisP(:)
        real(8), pointer :: fheadUp(:), fheadDn(:), length(:), area(:), HydRadius(:)
        real(8), pointer :: velocity(:), ManningsN(:), headslope(:), CMvelocity(:)
        real(8), pointer :: velocityBlend(:), svRatio(:)
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'adjust_velocity_smallvolume_blended'
        !%-----------------------------------------------------------------------------
        if (icrash) return              
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        thisP     => elemP(1:Npack,thisCol) !% only elements with small volumes
        fheadUp   => faceR(:,fr_Head_d)
        fheadDn   => faceR(:,fr_Head_u)
        length    => elemR(:,er_Length)
        area      => elemR(:,er_Area)
        HydRadius => elemR(:,er_HydRadius)
        velocity  => elemR(:,thisVelocityCol)
        ManningsN => elemR(:,er_SmallVolume_ManningsN)
        headslope => elemR(:,er_SmallVolume_HeadSlope)
        CMvelocity => elemR(:,er_SmallVolume_CMvelocity)    
        svRatio    => elemR(:,er_SmallVolumeRatio)
        velocityBlend => elemR(:,er_Temp01)
        !%-----------------------------------------------------------------------------
        !% Adjust ManningsN for small volume CM velocity.
        !% Use the larger of the actual roughness or the setting% value
        ManningsN(thisP) = setting%SmallVolume%ManningsN
        where (ManningsN(thisP) < elemR(thisP,er_Roughness))
            ManningsN(thisP) = elemR(thisP,er_Roughness)
        endwhere

        !% slope of the piezometric head
        headslope(thisP) = (fheadUp(thisP) - fheadDn(thisP)) / length(thisP)

        !% absolute chezy-manning velocity based on slope
        CMvelocity(thisP) = ( HydRadius(thisP)**(twothirdR) ) * sqrt(abs(headslope(thisP))) / ManningsN(thisP)
                    
        !% assign direction to CM velocity.
        CMvelocity(thisP) = sign(CMvelocity(thisP), headslope(thisP))

        !% blend the computed velocity with CM velocity
        velocityBlend(thisP) = svRatio(thisP) * velocity(thisP) &
                            + (oneR - svRatio(thisP)) * CMvelocity(thisP)

        !% use the smaller velocity value, with the sign of the CM velocity
        velocity(thisP) = min(abs(CMvelocity(thisP)), velocity(thisP))
        velocity(thisP) = sign(velocity(thisP),CMvelocity(thisP))

        elemYN(thisP,eYN_IsAdhocFlowrate) = .true.

        velocityBlend(thisP) = nullvalueR

        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine adjust_velocity_smallvolume_blended
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine adjust_zero_velocity_at_zero_volume &
        (Npack, thisCol, thisVelocityCol, thisVolumeCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% sets velocity to zero in near-zero volume s
        !%-----------------------------------------------------------------------------    
        integer, intent(in) :: Npack, thisCol, thisVelocityCol, thisVolumeCol
        integer, pointer :: thisP(:)
        real(8), pointer :: velocity(:), volume(:)
        logical, pointer :: isAdhocFlowrate(:)
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'adjust_zero_velocity_at_zero_volume'
        !%-----------------------------------------------------------------------------
        if (icrash) return              
        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%----------------------------------------------------------------------------- 
        thisP    => elemP(1:Npack,thisCol)
        volume   => elemR(:,thisVolumeCol)
        velocity => elemR(:,thisVelocityCol)
        isAdhocFlowrate => elemYN(:,eYN_IsAdhocFlowrate)
        !%----------------------------------------------------------------------------- 

        where (volume(thisP) <= setting%ZeroValue%Volume)
            velocity(thisP) = zeroR
            isAdhocFlowrate(thisP) = .true.
        endwhere    

        if (setting%Debug%File%adjust) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine adjust_zero_velocity_at_zero_volume
    !%
    !%==========================================================================
    !% END OF MODULE
    !%+=========================================================================
end module adjust