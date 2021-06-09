module runge_kutta2

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use adjust
    use update
    use face

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Runge-Kutta method for time-marching (ETM and AC)
    !%
    !% METHOD:
    !% 
    !%

    private

    public :: rk2_ETM_toplevel
    public :: rk2_AC_toplevel
    public :: rk2_ETMAC_toplevel

    contains
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    subroutine rk2_ETM_toplevel ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% single RK2 step for explicit time advance of SVE
        !%-----------------------------------------------------------------------------
        integer :: istep, whichTM  
        ! !%-----------------------------------------------------------------------------
        ! !% RK2 solution step 1   
        call rk2_ETM_step (istep=1)
        !% RK2 solution step 3 -- all aux variables for non-diagnostic
        call update_auxiliary_variables(whichTM=ETM)   

        !% RK2 solution step 4 -- all face interpolation
        if (num_images() == 1) then
            call face_interpolation_byMask(faceMaskCol=fm_all)
        else
            call face_interp_across_images()
            call face_interp_interior()
        end if
        
        !% RK2 solution step 5 -- update diagnostic elements and faces
        if (N_diag > 0) then
            call update_diagnostic_all()
        endif
        
        !% RK2 solution step X -- make ad hoc adjustments
        call adjust_values (whichTM=ETM)
        
        !% RK2 solution step 8 -- RK2 second step for ETM
        !% RK2 solution step 8(a)
        call rk2_ETM_step (istep=2)    
        !% RK2 solution step 8(c)
        call update_auxiliary_variables(whichTM=ETM)
        
        !% RK2 solution step 8(d,e) -- update all faces
        if (num_images() == 1) then      
            call face_interpolation_byMask (faceMaskCol=fm_all)
        else
            call face_interp_across_images()
            call face_interp_interior()
        endif
        
        !% RK2 solution step 9 -- update diagnostic elements and faces
        if (N_diag > 0) then 
            call update_diagnostic_all()
        endif

        !% RK2 solution step X -- make ad hoc adjustments
        call adjust_values (whichTM=ETM)
    
    end subroutine rk2_ETM_toplevel
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_AC_toplevel ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------
    
        !%-----------------------------------------------------------------------------
    end subroutine rk2_AC_toplevel    
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_ETMAC_toplevel ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------
        integer :: istep, whichTM, faceMaskCol, thisCol
        integer, pointer :: Npack

        !% N_ac is the number of AC elements
        !% N_etm is the number of ETM elements
        !% N_diag is the number of diagnoistic elements
        !%     
        !%-------------------------------------------------        
        !% step 1 -- RK2 step 1 for ETM
        istep=oneI
        if (N_etm > 0) then  
            call rk2_ETM_step (istep)
        endif
        
        !% step 2 -- RK2 step 1 for AC
        if (N_ac > 0) then
            !% step 2(a)
            call rk2_AC_step (istep)
            if (N_etm > 0) then        
                !% step 2(b,c) create time n+1(1*) consistency for AC 
                call rk2_ETM_extrapolate_to_fullstep()
            endif
        endif

        !% step 3 -- all aux variables for non-diagnostic
        call update_auxiliary_variables(whichTM=ALLtm)
        
        !% step 4 -- all face interpolation
        if (num_images() == 1) then
            call face_interpolation_byMask (faceMaskCol=fm_all)
        else
            call face_interp_across_images()
            call face_interp_interior()
        end if
        
        !% step 5 -- update diagnostic elements and faces
        if (N_diag > 0) then
            call update_diagnostic_all()
        endif
        
        !% step X -- make ad hoc adjustments
        call adjust_values (whichTM=ALLtm)
        
        !% step 6 -- RK2 step 2 for AC
        istep = 2
        if (N_ac > 0) then
            !% step 6(a)
            call rk2_AC_step (istep)
            if (N_etm > 0) then  
                !% step 6(b)
                call rk2_ETM_restore_to_midstep()
                !% step 6(c,d)
                call rk2_AC_interpolate_to_halfstep()
            endif
            !% step 6(e)
            call update_auxiliary_variables(whichTM=AC)      
            !% step 6(f,g) -- update faces for AC elements
            if (num_images() == 1) then 
                thisCol = fp_AC
                Npack => npack_faceP(thisCol)
                if (Npack > 0) then
                    call face_interpolation_byPack (thisCol, Npack)
                endif    
            else
                call face_interp_across_images()
                call face_interp_interior()
            endif
        endif  
        
        !% step 7 -- update diagnostic elements and faces
        if (N_diag > 0) then
            call update_diagnostic_all() 
        endif
        
        !% step 8 -- RK2 step 2 for ETM
        if (N_etm > 0) then
            !% step 8(a)
            call rk2_ETM_step (istep)
        
            if (N_ac > 0) then   
                !% step 8(b)
                call rk2_AC_restore_to_fullstep ()
            endif
            
            !% step 8(c)
            call update_auxiliary_variables(whichTM=ETM)
            
            !% step 8(d,e) -- update all faces
            if (num_images() == 1) then      
                call face_interpolation_byMask (faceMaskCol=fm_all)
            else
                call face_interp_across_images()
                call face_interp_interior()
            endif
        endif
        
        !% step 9 -- update diagnostic elements and faces      
        if (N_diag > 0) then   
            call update_diagnostic_all
        endif

        !% step X -- make ad hoc adjustments
        call adjust_values (whichTM=ALLtm)
                
        !%-----------------------------------------------------------------------------
    end subroutine rk2_ETMAC_toplevel    
    !%
    !%==========================================================================
    !% PRIVATE
    !%==========================================================================   
    !%  
    subroutine rk2_ETM_step (istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single rk2 step for ETM
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: istep
        integer :: tmType
        !%-----------------------------------------------------------------------------
        !%  
        !% perform the continuity step of the rk2 for ETM
        call rk2_ETM_continuity_step(istep)
        !% perform the momentum step of the rk2 for ETM
        call rk2_ETM_momentum_step(istep)      

    end subroutine rk2_ETM_step   
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_ETM_continuity_step (istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% perform the continuity step of the rk2 for ETM
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: istep   
        !%-----------------------------------------------------------------------------
        !%  
    end subroutine rk2_ETM_continuity_step       
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_ETM_momentum_step (istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% perform the momentum step of the rk2 for ETM
        integer, intent(in) :: istep

        !%-----------------------------------------------------------------------------
    
        !%-----------------------------------------------------------------------------
        !%      
    end subroutine rk2_ETM_momentum_step
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_AC_step (istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: istep
        !%-----------------------------------------------------------------------------
        !%          
        !%
    end subroutine rk2_AC_step
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_ETM_extrapolate_to_fullstep ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------
    
        !%-----------------------------------------------------------------------------
        !% 
    end subroutine rk2_ETM_extrapolate_to_fullstep      
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_ETM_restore_to_midstep ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------
    
        !%-----------------------------------------------------------------------------
        !%   
        !%
    end subroutine rk2_ETM_restore_to_midstep 

    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_AC_interpolate_to_halfstep()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------
    
        !%-----------------------------------------------------------------------------
        !%   
        !%
    end subroutine rk2_AC_interpolate_to_halfstep
    
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_AC_restore_to_fullstep()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------
    
        !%-----------------------------------------------------------------------------
        !%   
    end subroutine rk2_AC_restore_to_fullstep
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------
    
        !%-----------------------------------------------------------------------------
        !%      !%==========================================================================
    !% END OF MODULE
    !%+=========================================================================
end module runge_kutta2