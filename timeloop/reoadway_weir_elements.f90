module roadway_weir_elements

    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use utility, only: util_CLprint
    use utility_crash, only: util_crashpoint


    implicit none

    !% The discharge coefficients and submergence factors listed below were
    !% derived from Figure 10 in "Bridge Waterways Analysis Model: Research
    !% Report", U.S. Dept. of Transportation Federal Highway Administration
    !% Report No. FHWA/RD-86/108, McLean, VA, July 1986. 

    !% (these tables are copid over from SWMM5C and then moified.
    !% In SWMM5C, all the roadway weir calculations are done in CFS units.
    !% Since default units in SWMM5+ is SI, the 2nd colums (discharge coefficients) 
    !% of the tables were converted. The original unit of the discharge coefficients 
    !% were is ft^(1/2)/sec. Thus, the discharge coefficients were multipiled by
    !% 0.552 to convert them to m^(1/2)/sec.)

    !% Discharge Coefficients for (head / road width) <= 0.15
    integer, parameter :: N_Cd_LowPaved = 4
    real(8), dimension(4,2) :: Cd_LowPaved = (/ 0.0, 0.2, 0.7, 4.0, &
                                                1.57,1.63,1.67,1.68 &
                                             /)

    integer, parameter :: N_Cd_LowGravel = 8
    real(8), dimension(8,2) :: Cd_LowGravel = (/ 0.0, 0.5, 1.0, 1.5, &
                                                 2.0, 2.5, 3.0, 4.0, &
                                                 1.38,1.49,1.55,1.6, &
                                                 1.64,1.667,1.673,1.683 &
                                               /)
    !% Discharge Coefficients for (head / road width) > 0.15
    integer, parameter :: N_Cd_HighPaved = 2
    real(8), dimension(2,2) :: Cd_HighPaved = (/ 0.15,0.25, &
                                                 1.68,1.71  &
                                              /)

    integer, parameter :: N_Cd_HighGravel = 2
    real(8), dimension(2,2) :: Cd_HighGravel = (/ 0.15,0.30, &
                                                  1.63,1.71  &
                                               /)
    !% Submergence Factors
    !% (these are factor and unitless. Thus do not need correction)
    integer, parameter :: N_Sf_Paved = 9
    real(8), dimension(9,2) :: Sf_Paved = (/ 0.8, 0.85,0.9, 0.93, &
                                             0.95,0.97,0.98,0.99, &
                                             1.00,                &
                                             1.0, 0.98,0.92,0.85, &
                                             0.80,0.70,0.60,0.50, &
                                             0.40                 &
                                           /)
    integer, parameter :: N_Sf_Gravel = 12
    real(8), dimension(12,2) :: Sf_Gravel = (/ 0.75,0.80,0.83,0.86, &
                                               0.89,0.90,0.92,0.94, &
                                               0.96,0.98,0.99,1.00, &
                                               1., 0.985,0.97,0.93, &
                                               0.90,0.87,0.80,0.70, &
                                               0.60,0.50,0.40,0.24  &
                                             /)

    private

    public :: roadway_weir_flow

    contains

!%
!%==========================================================================
!% PUBLIC 
!%==========================================================================    
!% 
    subroutine roadway_weir_flow (eIdx)
    !%-----------------------------------------------------------------------------
    !% Description:
    !% 
    !% Find flow in roadway weirs
    !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx
        real(8), pointer    :: RoadHeight, RoadWidth, Head, NominalDSmHead 
        real(8), pointer    :: RectangularBreadth, disCoeff, Flowrate, dQdH
        integer, pointer    :: FlowDirection, RoadSurf
        real(8) :: HeadUp, HeadDn, cD
        logical :: useVariableCoeff
    !%-----------------------------------------------------------------------------
        dQdH                  => elemR(eIdx,er_dQdH)
        Flowrate              => elemR(eIdx,er_Flowrate)
        Head                  => elemR(eIdx,er_Head)
        RoadHeight            => elemSR(eIdx,esr_Weir_Zcrest) 
        RoadWidth             => elemSR(eIdx,esr_Wier_RoadWidth)
        disCoeff              => elemSR(eIdx,esr_Weir_Rectangular) 
        NominalDSmHead        => elemSR(eIdx,esr_Weir_NominalDownstreamHead)
        RectangularBreadth    => elemSR(eIdx,esr_Weir_RectangularBreadth)
        FlowDirection         => elemSI(eIdx,esi_Weir_FlowDirection)
        RoadSurf              => elemSI(eIdx, esi_Weir_RoadSurface)
    !%-----------------------------------------------------------------------------
        !% check if variable discharge coefficient is needed
        useVariableCoeff = .false.
        if ((RoadWidth > zeroR) .and. (RoadSurf .ne. NoRoadSurface)) useVariableCoeff = .true. 

        !% find water depths upstream and downstream of the roadway weir
        HeadUp = Head - RoadHeight
        HeadDn = NominalDSmHead - RoadHeight 
        cD     = disCoeff
        
        if (HeadUp > zeroR) then  
            if (useVariableCoeff) cD = get_roadway_discharge_coeff &
                (HeadUp, HeadDn, RoadWidth, RoadSurf)

            Flowrate = real(FlowDirection,8) * cD * RectangularBreadth * (HeadUp ** 1.5)
            !% find dQ/dH
            dQdH = 1.5 * Flowrate / HeadUp
        end if

    end subroutine roadway_weir_flow
!%
!%========================================================================== 
!% PRIVATE
!%==========================================================================    
!% 
    real(8) function get_roadway_discharge_coeff &
        (HeadUp, HeadDn, RoadWidth, RoadSurface) result (disCoeff)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: RoadSurface
        real(8), intent(in) :: HeadUp, HeadDn, RoadWidth
        real(8) :: hOverL, subFactor, dnHoverL
        !%-----------------------------------------------------------------------------
        !% set the initial submergence factor to 1
        subFactor = oneR
        !% if upstream head is zero, discharge coeff is zero
        if (HeadUp .le. zeroR) then
            disCoeff = zeroR
            return
        end if
        !% find the head/roadwidth ratio
        hOverL = HeadUp / RoadWidth

        !% find the discharge coeff from the tables based on head/roadwidth ratio
        if  (hOverL .le. 0.15) then
            if (RoadSurface == Paved) then
                disCoeff = get_table_val(hOverL, Cd_LowPaved, N_Cd_LowPaved)
            else
                disCoeff = get_table_val(hOverL, Cd_LowGravel, N_Cd_LowGravel)
            end if
        else
            if (RoadSurface == Paved) then
                disCoeff = get_table_val(hOverL, Cd_HighPaved, N_Cd_HighPaved)
            else
                disCoeff = get_table_val(hOverL, Cd_HighGravel, N_Cd_HighGravel)
            end if
        end if

        !% adjust for submergence
        if (HeadDn > zeroR) then
            dnHoverL = HeadDn / RoadWidth
            if (RoadSurface == Paved) then
                subFactor = get_table_val(dnHoverL, Sf_Paved, N_Sf_Paved)
            else
                subFactor = get_table_val(dnHoverL, Sf_Gravel, N_Sf_Gravel)
            end if
        end if

        !% multiply dis coeff with sub factor
        disCoeff = disCoeff * subFactor

    end function get_roadway_discharge_coeff
!%
!%========================================================================== 
!%==========================================================================    
!% 
    real(8) function get_table_val(xVal, table, nRow) result (yVal)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !%      This cide has been adapted from SWMM5C
        !%-----------------------------------------------------------------------------
        integer, intent(in) ::  nRow
        real(8), intent(in) :: xVal, table(:,:)
        integer :: ii
        real(8) :: x1, dx, y1, dy
        !%-----------------------------------------------------------------------------
        if (xVal .le. table(oneI,oneI)) then
            yVal = table(oneI,twoI)
        else if (xVal .ge. table(nRow,oneI)) then
             yVal = table(nRow,twoI)
        else
            do ii = 2,nRow
                if (xVal .le. table(ii,oneI)) then
                    x1 = table(ii-oneI,oneI)
                    dx = table(ii,oneI) - x1
                    y1 = table(ii-oneI,twoI)
                    dy = table(ii,twoI) - y1
                    yVal = y1 + (xVal - x1) * dy/dx
                    return
                end if
            end do
        end if

    
    end function get_table_val
!%
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module roadway_weir_elements