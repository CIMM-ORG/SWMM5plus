module utility_crash
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Utility routines for code crash control. These are preferred instead of
    !% a Fortran STOP command because the stop only operates on a single image
    !% which means a parallel simulation will hang at the nex SYNC command.
    !%
    !% Note that these cannot depend on anything except for define_globals 
    !% define_indexes, define_keys, and define_settings so that it 
    !% can be called by any module except these.
    !%
    !% Methods:
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
    !%==========================================================================
    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting

    implicit none

    private

    public :: util_crash_initialize
    public :: util_crashcheck
    public :: util_crashstop
    public :: util_crashpoint

    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine util_crash_initialize
        !%------------------------------------------------------------------
        !% Description:
        !% sets up the conditions for identifying blow-up crash of the
        !% simulation
        !% uses the size of integer*2 as a reasonable limit
        !%------------------------------------------------------------------
        !% Declarations:
            integer*2 :: ilimit = 0
        !%------------------------------------------------------------------

        setting%Crash%DepthMax               = real(huge(ilimit),8)
        setting%Crash%HeadMax                = setting%Crash%DepthMax + maxval(elemR(:,er_Zbottom))
        setting%Crash%FlowrateMax            = real(huge(ilimit),8)
        setting%Crash%PercentVelocityAtLimit = tenR

    end subroutine util_crash_initialize
!%    
!%==========================================================================
!%==========================================================================
!%    
    subroutine util_crashcheck (indexnumber)
        !%------------------------------------------------------------------
        !% Description:
        !% Checks for blow-up conditions and halts execution
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: indexnumber
            integer :: nlarge
            real(8) :: percentU
            logical :: iscrash
        !%------------------------------------------------------------------

        iscrash = .false.

        call util_crash_nantest(iscrash)

        !% --- blow up tests for individual cells
        call util_crash_blowuptest (er_Depth,    setting%Crash%DepthMax,    'depth',    iscrash)

        call util_crash_blowuptest (er_Head,     setting%Crash%HeadMax,     'head',     iscrash)

        call util_crash_blowuptest (er_Flowrate, setting%Crash%FlowrateMax, 'flowrate', iscrash)

        !% --- blow up test for too many cells close to velocity limit
        nlarge = count(      (abs(elemR(1:N_elem(this_image()),er_Velocity)) .ge. 0.98d0* setting%Limiter%Velocity%Maximum) &
                       .and. (abs(elemR(1:N_elem(this_image()),er_Velocity)) .ne. nullvalueR)  )
        percentU = real(nlarge,8) * onehundredR / real(N_elem(this_image()),8)

        if (percentU .ge. setting%Crash%PercentVelocityAtLimit) then
            iscrash = .true.
            print *, ' '
            print *, '========================================================== ', this_image()
            print *, '            SIMULATION BLOWING UP util_crashcheck'
            print *, 'Too many cells at artificial velocity limit'
            print *, nlarge, ' of ',N_elem(this_image()), ' cells on image ',this_image()
            print *, 'are at velocity near ',setting%Limiter%Velocity
        end if

        !% --- assign crashpoint
        if (iscrash) call util_crashpoint (indexnumber)

        !% --- set the crashI
        if (iscrash) crashI = 1
        
        !%------------------------------------------------------------------
    end subroutine util_crashcheck
 
!%==========================================================================
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

        !% --- identify crash on the crashing image
        if (crashI == 1) then
            print *, 'crash stop originating at', indexnumber,' on image ',this_image()
        end if

        sync all
        !% broadcast to all images and stop all
        call co_max(crashI)
        if (crashI==1) then
            print *, 'stopping image ',this_image()
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
        print *, 'crash point on image ',this_image(), ' at indexnumber...'
        print *, indexnumber
        if (num_images()==1) stop 
        
    end subroutine util_crashpoint
!%    
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine util_crash_blowuptest (thisCol, thislimit, blowupType, iscrash)
        !%------------------------------------------------------------------
        !% Description:
        !% Provides blowup test for elem(:,thisCol) with limit thislimit
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in)           :: thisCol
            real(8), intent(in)           :: thislimit
            logical, intent(inout)        :: iscrash
            character (len=*), intent(in) :: blowupType

            real(8), pointer              :: tempvalue(:)
            real(8) :: valuemax
            integer :: locOfMax(1), linknum, nodenum
        !%------------------------------------------------------------------
        !% Aliases 
            tempvalue => elemR(:,er_Temp01)
        !%------------------------------------------------------------------    
        !% --- store tempvalue with zeros for nullvalues
        tempvalue = abs(elemR(:,thisCol))
        where (tempvalue == nullvalueR)
            tempvalue = zeroR
        endwhere

        valuemax = maxval(tempvalue)
        locOfMax = maxloc(tempvalue)

        if (valuemax .ge. thislimit) then
            iscrash = .true.
            print *, ' '
            print *, '========================================================== ', this_image()
            print *, '      SIMULATION BLOWING UP util_crash_blowuptest'
            print *, trim(blowupType),' >    ',thislimit
            print *, 'at element index number', locOfMax
            print *, 'on image               ', this_image()
            linknum = elemI(locOfMax(1),ei_link_Gidx_SWMM)
            nodenum = elemI(locOfMax(1),ei_node_Gidx_SWMM)
            if (linknum .ne. nullvalueI) then
                print *, 'link index number is   ',linknum
                print *, 'SWMM input link ID is        ',trim(link%Names(linknum)%str)
            end if
            if (nodenum .ne. nullvalueI) then
                print *, 'node index number is  ',nodenum
                print *, 'SWMM input node ID is       ',trim(node%Names(nodenum)%str)
            end if
        end if
    
    end subroutine util_crash_blowuptest
!%   
!%==========================================================================
!%==========================================================================
!%
    subroutine util_crash_nantest (iscrash)
        !%------------------------------------------------------------------
        !% Description
        !% Tests for NaN values in the elemR array
        !% HACK -- should be expanded to faceR array as well
        !% HACK -- for comprehensive testing, we need to only do the static
        !% columns at the start, and then dynamic during the simulation.
        !%------------------------------------------------------------------
        !% Declarations:
            logical, intent(inout)        :: iscrash

            integer, pointer :: NN
            integer :: ii, kk

            character(len=14) :: ename(8)
            integer           :: eset(8)
        !%------------------------------------------------------------------
        !% Aliases
            NN => N_elem(this_image())
        !%------------------------------------------------------------------
        !% Preliminaries
            !% --- columns of elemR that will be checked for NaN
            eset = [er_Area, &
                    er_Depth, &
                    er_Flowrate, &
                    er_Head, &
                    er_Slotdepth, &
                    er_Velocity, &
                    er_Volume, &
                    er_WaveSpeed]

            ename = ['er_Area     ', &
                    'er_Depth    ', &
                    'er_Flowrate ', &
                    'er_Head     ', &
                    'er_Slotdepth', &
                    'er_Velocity ', &
                    'er_Volume   ', &
                    'er_WaveSpeed']
        !%------------------------------------------------------------------

        !% --- globally check for NaN where they shouldn't occur
        if (any(isnan(elemR(1:NN,eset)))) then
            iscrash = .true.
            !% --- cycle through the data types
            do ii=1,size(eset)
                if (any(isnan(elemR(1:NN,eset(ii))))) then
                    !% --- cycle through the individual elements
                    do kk=1,NN
                        if (isnan(elemR(kk,eset(ii)))) then
                            print *, 'NaN found in ',trim(ename(ii)),' on image ',this_image()
                            print *, 'element index ',kk
                        end if
                    end do
                end if
            end do
        end if
        
    end subroutine util_crash_nantest
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module utility_crash