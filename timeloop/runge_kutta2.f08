module runge_kutta2

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use adjust
    use update
    use face
    use lowlevel_rk2

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Runge-Kutta time-marching method for both ETM and AC approaches with RK2
    !%----------------------------------------------------------------------------- 

    private

    public :: rk2_toplevel_ETM
    public :: rk2_toplevel_AC
    public :: rk2_toplevel_ETMAC

    contains
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    subroutine rk2_toplevel_ETM ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% single RK2 step for explicit time advance of SVE
        !%-----------------------------------------------------------------------------
        integer :: istep, whichTM  
        ! !%-----------------------------------------------------------------------------
        ! !% RK2 solution step 1   
        istep=1
        call rk2_step_ETM (istep)
        !% RK2 solution step 3 -- all aux variables for non-diagnostic
        call update_auxiliary_variables (ETM)   

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
        call adjust_values (ETM)
        
        !% RK2 solution step 8 -- RK2 second step for ETM
        !% RK2 solution step 8(a)
        istep=2
        call rk2_step_ETM (istep)    
        !% RK2 solution step 8(c)
        call update_auxiliary_variables(ETM)
        
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
    
    end subroutine rk2_toplevel_ETM
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_toplevel_AC ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------
    
        !%-----------------------------------------------------------------------------
    end subroutine rk2_toplevel_AC    
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_toplevel_ETMAC ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------
        integer :: istep, whichTM, faceMaskCol, thisCol
        integer, pointer :: Npack
        !%-----------------------------------------------------------------------------    
        !% step 1 -- RK2 step 1 for ETM
        istep=1
        if (N_etm > 0) then  
            call rk2_step_ETM (istep)
        endif
        
        !% step 2 -- RK2 step 1 for AC
        if (N_ac > 0) then
            !% step 2(a)
            call rk2_step_AC (istep)
            if (N_etm > 0) then        
                !% step 2(b,c) create time n+1(1*) consistency for AC 
                call rk2_extrapolate_to_fullstep_ETM()
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
        istep=2
        if (N_ac > 0) then
            !% step 6(a)
            call rk2_step_AC (istep)
            if (N_etm > 0) then  
                !% step 6(b)
                call rk2_restore_to_midstep_ETM()
                !% step 6(c,d)
                call rk2_interpolate_to_halfstep_AC()
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
            call rk2_step_ETM (istep)
        
            if (N_ac > 0) then   
                !% step 8(b)
                call rk2_restore_to_fullstep_AC ()
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
    end subroutine rk2_toplevel_ETMAC    
    !%
    !%==========================================================================
    !% PRIVATE
    !%==========================================================================   
    !%  
    subroutine rk2_step_ETM (istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single rk2 step for ETM
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: istep
        integer :: tmType
        !%-----------------------------------------------------------------------------
        !%  
        !% perform the continuity step of the rk2 for ETM
        call rk2_continuity_step_ETM(istep)
        !% perform the momentum step of the rk2 for ETM
        call rk2_momentum_step_ETM(istep)      

    end subroutine rk2_step_ETM   
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_step_AC (istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Single step of RK2 for AC method
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: istep
        !%-----------------------------------------------------------------------------
        !% AC continuity
        call rk2_continuity_step_AC(istep)
        !% AC momentum
        call rk2_momentum_step_AC(istep)

    end subroutine rk2_step_AC
    !%       
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_continuity_step_ETM (istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% perform the continuity step of the rk2 for ETM
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: istep   
        integer, pointer :: thisPackCol, Npack
        !%-----------------------------------------------------------------------------
        !%          
        !% Baseline continuity source:
        !% Compute net flowrates for channels, conduits and special elements
        thisPackCol => col_elemP(ep_CC_H_ETM)
        Npack => npack_elemP(thisPackCol)
        if (Npack > 0) then
            call ll_continuity_netflowrate_CC (er_SourceContinuity, thisPackCol, Npack)
        endif
            
        !% compute net flowrates for junction mains   
        thisPackCol => col_elemP(ep_JM_ETM)
        Npack => npack_elemP(thisPackCol)
        if (Npack > 0) then
            call ll_continuity_netflowrate_JM (er_SourceContinuity, thisPackCol, Npack)
        endif
        
        !% Solve for volume in ETM step
        thisPackCol => col_elemP(ep_CCJM_H_ETM)
        Npack => npack_elemP(thisPackCol)
        if (Npack > 0) then
            call ll_continuity_volume_CCJM_ETM (er_Volume, thisPackCol, Npack, istep)
        endif
        
    end subroutine rk2_continuity_step_ETM       
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_continuity_step_AC (istep) 
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single RK2 continuity step for AC method
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: istep
        integer, pointer ::  thisMaskCol, thisPackCol, Npack
        !%-----------------------------------------------------------------------------
        !%          
        !% Compute net flowrates for channels, conduits and special elements
        thisPackCol => col_elemP(ep_CC_AC)
        Npack => npack_elemP(thisPackCol)
        if (Npack > 0) then
            call ll_continuity_netflowrate_CC (er_SourceContinuity, thisPackCol, Npack)
        endif
            
        !% compute net flowrates for junction mains    
        thisPackCol => col_elemP(ep_JM_AC)
        Npack => npack_elemP(thisPackCol)
        if (Npack > 0) then
            call ll_continuity_netflowrate_JM (er_SourceContinuity, thisPackCol, Npack)
        endif
        
        thisPackCol => col_elemP(ep_CCJM_H_AC_open)
        Npack => npack_elemP(thisPackCol)
        if (Npack > 0) then
            !% unique continuity source terms for AC open channel
            call ll_continuity_add_source_CCJM_AC_open (er_SourceContinuity, thisPackCol, Npack)
            !% unique continuity gamma terms for AC open channel
            call ll_continuity_add_gamma_CCJM_AC_open (er_GammaC, thisPackCol, Npack)        
            !% solve for volume in AC open step
            call ll_continuity_volume_CCJM_AC_open (er_Volume, thisPackCol, Npack, istep)
        endif
        
        thisPackCol => col_elemP(ep_CCJM_H_AC_surcharged)
        Npack => npack_elemP(thisPackCol)
        if (Npack > 0) then
            !% unique continuity source terms for AC surcharged channel head
            call ll_continuity_add_source_CCJM_AC_surcharged (er_SourceContinuity, thisPackCol, Npack)   
            !% solve for head in AC surcharged step 
            call ll_continuity_head_CCJM_AC_surcharged (er_Head, thisPackCol, Npack, istep)
        endif
    end subroutine rk2_continuity_step_AC
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_momentum_step_ETM (istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% perform the momentum step of the rk2 for ETM (or ETM part of ETMAC)
        integer, intent(in) :: istep
        integer, pointer :: thisPackCol, Npack
        integer, parameter :: thisMethod = ETM  
        !%-----------------------------------------------------------------------------
        thisPackCol => col_elemP(ep_CC_Q_ETM)
        Npack       => npack_elemP(thisPackCol)  

        !%-----------------------------------------------------------------------------
        !%    
        if (Npack > 0) then
            !% momentum K source terms for different methods for ETM
            call ll_momentum_Ksource_CC (er_Ksource, thisPackCol, Npack)           
            !% Common source for momentum on channels and conduits for ETM
            call ll_momentum_source_CC (er_SourceMomentum, thisPackCol, Npack)            
            !% Common Gamma for momentum on channels and conduits for  ETM 
            call ll_momentum_gamma_CC (er_GammaM, thisPackCol, Npack)
            !% Advance flowrate to n+1/2 for conduits and channels with ETM
            call ll_momentum_solve_CC (er_Velocity, thisPackCol, Npack, thisMethod, istep)
            !% velocity for ETM time march
            call ll_momentum_velocity_CC (er_Velocity, thisPackCol, Npack)  
        endif        

    end subroutine rk2_momentum_step_ETM
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_momentum_step_AC (istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: istep
        integer, pointer    :: thisCol, Npack
        integer, parameter  :: thisMethod = AC
        !%-----------------------------------------------------------------------------
        !%      
        thisCol => col_elemP(ep_CC_Q_AC)
        Npack => npack_elemP(thisCol)
        if (Npack > 0) then
            !% momentum K source terms for different methods
            call ll_momentum_Ksource_CC (er_Ksource, thisCol, Npack)         
            !% Source for momentum on channels and conduits
            call ll_momentum_source_CC (er_SourceMomentum, thisCol, Npack)           
            !% Gamma for momentum on channels and conduits
            call ll_momentum_gamma_CC (er_GammaM, thisCol, Npack)         
            !% additional momentum Gamma for AC time-march on channels and conduits
            call ll_momentum_add_gamma_CC_AC (er_GammaM, thisCol, Npack)            
            !% additional momentum source terms for AC time-march on channels and conduits
            call ll_momentum_add_source_CC_AC (er_SourceMomentum, thisCol, Npack)          
            !% AC elements advance flowrate to n+1(*) for conduits and channels
            call ll_momentum_solve_CC (er_Velocity, thisCol, Npack, thisMethod, istep)
        endif
           
        !% adjust zero elements
        thisCol => col_elemP(ep_CCJM_H_AC)
        Npack => npack_elemP(thisCol)
        if (Npack > 0) then
            call adjust_nearzero_volume (er_Volume, thisCol, Npack)
        endif
            
        thisCol => col_elemP(ep_CC_Q_AC)
        Npack => npack_elemP(thisCol)
        if (Npack > 0) then
            !% velocity for AC time march
            call ll_momentum_velocity_CC (er_Velocity, thisCol,Npack) 
        endif
          
    end subroutine rk2_momentum_step_AC
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_extrapolate_to_fullstep_ETM ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Finds ETM elements at time n+1/2 that are adjacent to fAC 
        !% Makes temporary store of data for Q, H, V at n+1/2
        !% overwrites the Q, H, V location with an extrapolation to n+1.
        !%-----------------------------------------------------------------------------
        integer, pointer :: thisCol, Npack
        !%-----------------------------------------------------------------------------
        thisCol => col_elemP(ep_CCJB_eETM_i_fAC)
        Npack => npack_elemP(thisPackCol)
            
        if (Npack > 0) then
            !% temporary storage of n+1/2 data
     !       call store_in_temporary (thisCol, Npack)
            
            !% extrapolation
     !       call extrapolate_values (thisCol, Npack)
            
     !       call update_auxiliary_variables_byPack (thisCol, Npack)
        endif

    end subroutine rk2_extrapolate_to_fullstep_ETM      
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_restore_to_midstep_ETM ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------
    
        !%-----------------------------------------------------------------------------
        !%   
        !%
    end subroutine rk2_restore_to_midstep_ETM

    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_interpolate_to_halfstep_AC()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------
    
        !%-----------------------------------------------------------------------------
        !%   
        !%
    end subroutine rk2_interpolate_to_halfstep_AC
    
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine rk2_restore_to_fullstep_AC()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------
    
        !%-----------------------------------------------------------------------------
        !%   
    end subroutine rk2_restore_to_fullstep_AC
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
end module runge_kutta2