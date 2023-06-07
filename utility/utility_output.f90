module utility_output
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% common output functions
    !%==========================================================================
    use define_indexes
    use define_keys
    use define_globals
    use define_settings
    use define_types

    implicit none

    private

    public :: util_output_must_report

contains
!%
!%==========================================================================
!%
    function util_output_must_report() result(report)
        !%------------------------------------------------------------------
        !% Description:
        !% determines whether report is needed for Report.ThisStep report at 
        !% interval Report.TimeInterval
        !%------------------------------------------------------------------
            logical :: report
            integer, pointer :: reportStep
            real(8) :: timeNow, reportDt, startReport
        !%------------------------------------------------------------------
        reportStep  => setting%Output%Report%ThisStep
        timeNow     = setting%Time%Now
        reportDt    = setting%Output%Report%TimeInterval
        startReport = setting%Output%Report%StartTime

        if ((timeNow >= reportDt * (reportStep + 1)) .and. (timeNow > startReport))then
            report = .true.
        else if (timeNow == startReport) then
            report = .true.
            reportStep = -1
        else
            report = .false.
        end if

    end function util_output_must_report
!%
!%==========================================================================
!% END MODULE
!%==========================================================================
!%
end module utility_output
