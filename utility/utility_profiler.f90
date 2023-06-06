
module utility_profiler
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Computes cumulative elapsed CPU time on each processor for procedures with
    !% indexes in the pf_#subroutine name# enumerated in define_indexes
    !%==========================================================================
    use define_indexes
    use define_globals

    implicit none

    private
    public :: util_profiler_start
    public :: util_profiler_stop
    public :: util_profiler_print_summary 

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine util_profiler_start (this_procedure)
        !%------------------------------------------------------------------
        !% Description
        !% starts timing for the current processor and procedure  
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: this_procedure
        !%------------------------------------------------------------------    

        call cpu_time (profiler_data(pfr_thisstart,this_procedure))
        
    end subroutine util_profiler_start
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_profiler_stop (this_procedure)
        !%------------------------------------------------------------------
        !% Description
        !% stops timing for the current processor and procedure and records accumulation  
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: this_procedure
        !%------------------------------------------------------------------  
             
        call cpu_time (profiler_data(pfr_thisend,this_procedure))
        profiler_data(pfr_cumulative,this_procedure) =        &
              profiler_data(pfr_cumulative,this_procedure)    &
            + profiler_data(pfr_thisend   ,this_procedure)    &
            - profiler_data(pfr_thisstart ,this_procedure)

    end subroutine util_profiler_stop
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_profiler_print_summary ()
        !%------------------------------------------------------------------
        !% Description:
        !% provides command-line listing of profiling
        !%------------------------------------------------------------------ 
        !% Declarations:
            integer :: ii, mp,kk  
        !%------------------------------------------------------------------ 

        do mp = 1, num_images()
            if (this_image() == mp) then
                print *, '-------------- image ',mp ,' profile data------------------'
                do kk =1,3
                    print *, 'Level ',kk,' procedures'
                    do ii =1,Ncol_pf
                        if (profiler_procedure_level(ii) == kk) then
                            print *, profiler_data(pfr_cumulative,ii),  trim(profiler_procedure_name(ii))
                        endif
                    end do
                end do
                print *, 'unclassified levels'
                do ii =1,Ncol_pf
                    if ((profiler_procedure_level(ii) < 1) .or. (profiler_procedure_level(ii) > 3)) then
                        print *, profiler_data(pfr_cumulative,ii),  &
                                 trim(profiler_procedure_name(ii)), &
                                 profiler_procedure_level(ii)
                    end if    
                end do
            end if
        end do

    end subroutine util_profiler_print_summary
!%
!%==========================================================================
end module utility_profiler

