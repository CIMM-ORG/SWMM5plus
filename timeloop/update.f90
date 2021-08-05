module update

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use geometry
    use adjust

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Updates values during timeloop of hydraulics.
    !%

    private
    
    public :: update_auxiliary_variables

    real(8), pointer :: grav => setting%constant%gravity

    contains
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    !%
    subroutine update_auxiliary_variables (whichTM)
        !%-----------------------------------------------------------------------------
        !% Description: 
        !% 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: whichTM  !% indicates which Time marching method
        integer, pointer :: thisCol_all
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'update_auxiliary_variables'
        if (setting%Debug%File%update) print *, '*** enter ',this_image(), subroutine_name     
        !%-----------------------------------------------------------------------------
        !%  
        !% update the head (non-surcharged) and geometry

        call geometry_toplevel (whichTM)

        !% adjust velocity with limiters and small volume treatment
        call adjust_velocity (whichTM, er_Velocity, er_Volume)

        !% set packed column for updated elements
        select case (whichTM)
            case (ALLtm)
                thisCol_all => col_elemP(ep_CC_ALLtm)
            case (ETM)
                thisCol_all => col_elemP(ep_CC_ETM)
            case (AC)
                thisCol_all => col_elemP(ep_CC_AC)
            case default
                print *, 'error, default case should not be reached'
                stop 7489
        end select

        !% Compute the flowrate on CC.
        !% Note that JM should have 0 flowrate and JB has lagged flowrate at this point.
        !% The JB flowrate is not updated until after face interpolation
        call update_CC_element_flowrate (thisCol_all)

        !% compute element Froude numbers for CC, JM
        call update_Froude_number_element (thisCol_all)

        !% compute element face interpolation weights on CC, JM
        call update_interpolation_weights_element (thisCol_all, whichTM)
   
        if (setting%Debug%File%update)  print *, '*** leave ', subroutine_name
    end subroutine update_auxiliary_variables
    !%
    !%==========================================================================
    !% PRIVATE
    !%==========================================================================   
    !%  
    subroutine update_CC_element_flowrate (thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisCol
        !%-----------------------------------------------------------------------------
        integer, pointer ::  Npack, thisP(:)
        real(8), pointer :: flowrate(:), velocity(:), area(:)
        !%-----------------------------------------------------------------------------   
        flowrate => elemR(:,er_Flowrate)
        velocity => elemR(:,er_Velocity)
        area     => elemR(:,er_Area)   
        !%-----------------------------------------------------------------------------    
        Npack => npack_elemP(thisCol)
        if (Npack > 0) then
            thisP    => elemP(1:Npack,thisCol)
            flowrate(thisP) = area(thisP) * velocity(thisP)
        endif

    end subroutine update_CC_element_flowrate
    !%
    !%==========================================================================   
    !%==========================================================================   
    !%  
    subroutine update_Froude_number_element (thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes Froude number on each element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisCol
        integer, pointer :: Npack, thisP(:)
        real(8), pointer :: Froude(:), velocity(:), depth(:)
        !%-----------------------------------------------------------------------------
        Froude   => elemR(:,er_FroudeNumber)
        velocity => elemR(:,er_Velocity)
        depth    => elemR(:,er_ell)  !% Use the ell value (modified hydraulic depth)
        !%-----------------------------------------------------------------------------
    
        Npack => npack_elemP(thisCol)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisCol)
            Froude(thisP) = velocity(thisP) / sqrt(grav * depth(thisP))
        endif
    
    end subroutine update_Froude_number_element
    !%   
    !%==========================================================================   
    !%==========================================================================   
    !%  
    subroutine update_interpolation_weights_element (thisCol, whichTM)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes the interpolation weights on each element form CC, JM
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisCol, whichTM
        integer, pointer :: Npack, Npack2, thisCol_AC,  thisP(:), thisP2(:)
        real(8), pointer :: velocity(:), wavespeed(:), depth(:), length(:)
        real(8), pointer :: w_uQ(:), w_dQ(:),  w_uG(:), w_dG(:),  w_uH(:), w_dH(:)
        !%-----------------------------------------------------------------------------
        velocity  => elemR(:,er_Velocity)
        wavespeed => elemR(:,er_WaveSpeed)
        depth     => elemR(:,er_ell)  !% modified hydraulic depth!
        length    => elemR(:,er_Length)
        w_uQ      => elemR(:,er_InterpWeight_uQ)
        w_dQ      => elemR(:,er_InterpWeight_dQ)
        w_uG      => elemR(:,er_InterpWeight_uG)
        w_dG      => elemR(:,er_InterpWeight_dG)
        w_uH      => elemR(:,er_InterpWeight_uH)
        w_dH      => elemR(:,er_InterpWeight_dH)       
        !%-----------------------------------------------------------------------------
        !% 2nd cases needed for handling surcharged AC elements and using the celerity
        !% multiplier of the AC method for the wavespeed
        select case (whichTM)
            case (ALLtm)
                thisCol_AC =>  col_elemP(ep_Surcharged_AC)
            case (ETM)
                !% no effect
            case (AC)
                thisCol_AC =>  col_elemP(ep_Surcharged_AC)
            case default
                print *, 'error, case default should not be reached.'
                stop 3987
        end select    
    
        Npack => npack_elemP(thisCol)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisCol)

            !% wavespeed at modified hydraulic depth (ell)    
            wavespeed(thisP) = sqrt(grav * depth(thisP))
        
            !% modify wavespeed for surcharged AC cells
            if (whichTM .ne. ETM) then
                Npack2 => npack_elemP(thisCol_AC)
                if (Npack2 > 0) then
                    thisP2 => elemP(1:Npack2,thisCol_AC)
                    wavespeed(thisP2) = wavespeed(thisP2) * setting%ACmethod%Celerity%RC
                endif    
            endif
            
            !% timescale interpolation weights for flowrate
            w_uQ(thisP) = - onehalfR * length(thisP)  / (velocity(thisP) - wavespeed(thisP))
            w_dQ(thisP) = + onehalfR * length(thisP)  / (velocity(thisP) + wavespeed(thisP))
            
            !% apply limiters to timescales
            where (w_uQ(thisP) < zeroR)
                w_uQ(thisP) = setting%Limiter%InterpWeight%Maximum
            endwhere
            where (w_uQ(thisP) < setting%Limiter%InterpWeight%Minimum)    
                w_uQ(thisP) = setting%Limiter%InterpWeight%Minimum
            endwhere
            where (w_uQ(thisP) > setting%Limiter%InterpWeight%Maximum)    
                w_uQ(thisP) = setting%Limiter%InterpWeight%Maximum
            endwhere    
            
            !% timescale interpolation for geometry are identical to flowrate
            !% but may be modified elsewhere
            w_uG(thisP) = w_uQ(thisP)
            w_dG(thisP) = w_dQ(thisP)
            
            !% head uses length scale interpolation
            !% This shouldn't need limiters.
            w_uH(thisP) = onehalfR * length(thisP)
            w_dH(thisP) = onehalfR * length(thisP)
            
        endif

    end subroutine update_interpolation_weights_element
    !%   !%
    !%==========================================================================
    !% END OF MODULE
    !%+=========================================================================
end module update