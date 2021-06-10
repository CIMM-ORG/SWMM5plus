module update

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use geometry

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Updates values during timeloop of hydraulics.
    !%

    private

    public :: update_auxiliary_variables
    public :: update_diagnostic_all

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
        !%-----------------------------------------------------------------------------
        !%  
        !% update the head (non-surcharged) and geometry
        call geometry_toplevel (whichTM)

        !% adjust velocity with limiters and small volume treatment
        !call adjust_velocity (whichTM)

        !% Compute the flowrate on CC.
        !% Note that JM should have 0 flowrate and JB has lagged flowrate
        !% at this point. 
        !% The JB flowrate is not updated until after face !% interpolation
        !call update_CC_element_flowrate (whichTM)

        !% compute element Froude numbers
        !call element_Froude_number

        !% compute element face interpolation weights
        !call element_interpolation_weights
   
    end subroutine update_auxiliary_variables  
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine update_diagnostic_all ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%  
    end subroutine update_diagnostic_all   
    !%
    !%==========================================================================
    !% PRIVATE
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
    !% END OF MODULE
    !%+=========================================================================
end module update