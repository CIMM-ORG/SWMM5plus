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
    public :: adjust_nearzero_volume

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
        !%-----------------------------------------------------------------------------
        
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
        endif
        
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
        endif
        
    end subroutine adjust_values    
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine adjust_nearzero_volume (er_Volume, thisCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: er_Volume, thisCol, Npack
        !%-----------------------------------------------------------------------------
        !%  
    end subroutine adjust_nearzero_volume
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
        real(8), pointer :: coef
        real(8), pointer :: faceFlow(:), elemFlow(:), elemVel(:), w_uQ(:), w_dQ(:), elemArea(:)
        !%-----------------------------------------------------------------------------    
        select case (whichTM)
            case (ALLtm)
                thisCol => col_elemP(ep_CCJB_ALLtm)
            case (ETM)
                thisCol => col_elemP(ep_CCJB_ETM)
            case (AC)
                thisCol => col_elemP(ep_CCJB_AC)
            case default
                print *, 'error, this default case should not be reached'
                stop 9239
        end select
    
        !% coefficient for the blending adjustment (between 0.0 and 1.0)
        !% if coef == 1 then the V-shape element flowrate is replaced by the weighted
        !% average of its faces.
        coef => setting%Adjust%Flowrate%Coef
        
        if (coef > zeroR) then      
            Npack => col_elemP(thisCol)
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
            endif
        endif
        
    end subroutine adjust_Vshaped_flowrate
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine adjust_Vshaped_head_surcharged (whichTM)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: whichTM
        integer, pointer :: thisCol, Npack
        integer, pointer :: thisP(:), mapUp(:), mapDn(:)
        real(8), pointer :: coef
        real(8), pointer :: faceHeadUp(:), faceHeadDn(:), elemHead(:), elemVel(:)
        real(8), pointer :: w_uH(:), w_dH(:)
        !%-----------------------------------------------------------------------------
        select case (whichTM)
            case (ALLtm)
                thisCol => col_elemP(ep_CCJB_ALLtm_surcharged)
            case (ETM)
                thisCol => col_elemP(ep_CCJB_ETM_surcharged)
            case (AC)
                thisCol => col_elemP(ep_CCJB_AC_surcharged)
            case default
                print *, 'error, this default case should not be reached'
                stop 2394
        end select    
    
        !% coefficient for the blending adjustment (between 0.0 and 1.0)
        !% if coef == 1 then the V-shape element flowrate is replaced by the weighted
        !% average of its faces.
        coef => setting%Adjust%Head%Coef
        
        if (coef > zeroR) then       
            Npack = col_elemP(thisCol)
            if (Npack > 0) then
                thisP      => elemP(1:Npack,thisCol)
                mapUp      => elemI(:,ei_Mface_uL)
                mapDn      => elemI(:,ei_Mface_dL)    
                faceHeadUp => faceR(:,fr_Head_u)  
                faceHeadDn => faceR(:,fr_Head_d)          
                elemHead   => elemR(:,er_Head)    
                w_uH => elemR(:,er_InterpWeight_uH)
                w_dH => elemR(:,er_InterpWeight_dH)
                
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
            endif
        endif
        
    end subroutine
    !%
    !%==========================================================================
    !% END OF MODULE
    !%+=========================================================================
end module adjust