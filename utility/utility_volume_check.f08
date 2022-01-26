module utility_volume_check

    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting


!% OBSOLETE?

    implicit none

!-----------------------------------------------------------------------------
!
! Description:
!   Utility for check volume conservation
!
!-----------------------------------------------------------------------------

    private

    public :: util_volume_conservation

    
    contains
!%
!%==========================================================================  
!% PUBLIC
!%==========================================================================  
!%
    ! subroutine util_volume_conservation(diagnostic_type)
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% base on the given diagnostic_type, calculate the volume conservation
    !     !% diagnostic_type == 0 -> initialize all the variables 
    !     !% diagnostic_type == 1 -> initial condition volume
    !     !% disgnostic_type == 2 -> timeloop volume calculation (every timestep)
    !     !%-----------------------------------------------------------------------------
    !     integer, intent(in) :: diagnostic_type
    !     real(8) :: total_volume 
    !     real(8) :: channel_pipe_volume, weir_volume, orfice_volume
    !     real(8) :: JM_volume, JB_volume

    !     character(64) :: subroutine_name = 'util_volume_conservation'

    !     if (icrash) return
    !     if (setting%Debug%File%utility_volume_check) print *, '*** enter ',subroutine_name

    !     select case(diagnostic_type)
    !     case (0) ! this is for initialization
    !     case (1) ! Initial condition volume 
    !     case (2) ! Volume recorder in timeloop 
    !     case (3) ! Simply record max/min volume during the simulation
    !     case default
    !     end select       
             
    !     if (count() .and. count())

  
    !     channel_pipe_volume = sum()




    !     if (setting%Debug%File%utility_volume_check) print *, '*** leave ',subroutine_name

    ! end subroutine util_volume_conservation
!%
!%==========================================================================  
!% PRIVATE
!%==========================================================================  
!%  
    subroutine util_volume_bc(diagnostic_type)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Compute the volume from the boundary condition
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: diagnostic_type
        real(8) :: bc_volume

        character(64) :: subroutine_name = 'util_volume_bc'

        if (icrash) return
        if (setting%Debug%File%utility_volume_check) print *, '*** enter ',subroutine_name



        if (setting%Debug%File%utility_volume_check) print *, '*** leave ',subroutine_name
    end subroutine util_volume_bc


