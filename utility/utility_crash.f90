module utility_crash

    use define_globals

    implicit none

!%-----------------------------------------------------------------------------
!% Description:
!% Utility routines for code crash control
!% note that these cannot depend on anything except for define_globals so that it
!% can be called by any module except for defin_globals.f90
!%
!% USE: inside a conditional when a stop is needed use...
!%   call util_crashpoint(1093874) !% 1093874 is a unique index number for crash point
!%   This sets a "crash" condition where a single-processor code will immediately stop 
!%   with the index number displayed at the command line (e.g., stop 1093874).
!%   However, for multiprocessor the crashpoint sets crashI=1 on the local processor,
!%   which ensures that at the next call to util_crashstop() that all processors
!%   will be stopped.
!%
!% USE: outside of conditionals a stop is called with...
!%   call util_crashtop(93872) !% 93872 is a unique index number for the crash point.
!%   This will check to see if crashI==1 on ANY processor, and if so it will stop
!%   all the processors.
!%-----------------------------------------------------------------------------

    private

    public :: util_crashstop
    public :: util_crashpoint

    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine util_crashstop (indexnumber)
        !%------------------------------------------------------------------
        !% Description:
        !% checks for a crash condition across processors and stops all
        !% processors if any one has encountered a crash condition
        !% Note: requires a sync, so should NOT be called inside a conditional
        !% such as select-case or if-then.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: indexnumber !% unique ID of crash location
        !%------------------------------------------------------------------
        !% broadcast without sync
        call co_max(crashI)
        if (crashI==1) then
            print *, 'crash stop at...'
            stop indexnumber
        end if
        
    end subroutine util_crashstop
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_crashpoint (indexnumber)
        !%------------------------------------------------------------------
        !% Description:
        !% Creates immediate stop if this is a single processor run,
        !% otherwise prints the index number and sets the crashI=1 condition
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: indexnumber !% unique ID of crash location
        !%------------------------------------------------------------------
        crashI=1
        if (num_images()==1) then
            stop indexnumber
        else
            print *, 'crash point at...'
            print *, indexnumber
        end if
        
    end subroutine util_crashpoint
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module utility_crash