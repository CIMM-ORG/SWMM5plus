module hydrology

    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use interface_
    !use utility

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Handles hydrology functions that work with SWMM-C interface
    !%
    !% METHOD: 
    !% 
    !%

    private

    public :: hydrology_runoff
   

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine hydrology_runoff()
    !     !%-------------------------------------------------------------------
    !     !% Description:
    !     !% Takes the hydrology runoff stored in the sbcatchR array and 
    !     !% applies it to the lateral inflows in the elemR array for use
    !     !% in hydraulics
    !     !%-------------------------------------------------------------------
    !     !% Declarations
    !         integer :: ii 
    !         integer, pointer :: thisE
    !         real(8), pointer :: latFlowRate(:), subcatchFlowrate(:)
    !         character(64) :: subroutine_name = 'hydrology_runoff'
    !     !%-------------------------------------------------------------------
    !     !% Preliminaries
    !         
    !         if (setting%Debug%File%adjust) &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%-------------------------------------------------------------------
    !     !% Aliases
    !     !%-------------------------------------------------------------------              
    !         latFlowRate      => elemR(:,er_FlowrateLateral)
    !         subcatchFlowRate => subcatchR(:,sr_RunoffRate_baseline) 
    !     !%-------------------------------------------------------------------  
    !     !% Multiple subcatchments can drain to a single element, but each
    !     !% subcatchment connects to only one element, so we cannot array
    !     !% process between subcatchR and elemR
    !     do ii =1,setting%SWMMinput%N_subcatch
    !         thisE => subcatchI(ii,si_runoff_elemIdx) 
    !         latFlowRate(thisE) = latFlowRate(thisE) + subcatchFlowRate(ii)
    !     end do
        
    !     !%-------------------------------------------------------------------  
    !     !% Closing                
    !     if (setting%Debug%File%adjust) &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine hydrology_runoff
!%
!%==========================================================================
!%==========================================================================  
!%    
    ! subroutine
    !     !%-------------------------------------------------------------------
    !     !% Description:
    !     !% 
    !     !%-------------------------------------------------------------------
    !     !% Declarations
    !         character(64) :: subroutine_name = 'hydro_'
    !     !%-------------------------------------------------------------------
    !     !% Preliminaries
    !        
    !         if (setting%Debug%File%adjust) &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%-------------------------------------------------------------------
    !     !% Aliases
    !     !%-------------------------------------------------------------------              
     
        
    !     !%-------------------------------------------------------------------  
    !     !% Closing                
    !     if (setting%Debug%File%adjust) &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine 
!%
!%==========================================================================
!% PRIVATE below here
!%==========================================================================  
!%    
    ! subroutine
    !     !%-------------------------------------------------------------------
    !     !% Description:
    !     !% 
    !     !%-------------------------------------------------------------------
    !     !% Declarations
    !         character(64) :: subroutine_name = 'hydro_'
    !     !%-------------------------------------------------------------------
    !     !% Preliminaries
    !         
    !         if (setting%Debug%File%adjust) &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%-------------------------------------------------------------------
    !     !% Aliases
    !     !%-------------------------------------------------------------------              
     
        
    !     !%-------------------------------------------------------------------  
    !     !% Closing                
    !     if (setting%Debug%File%adjust) &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine 
!%
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module hydrology