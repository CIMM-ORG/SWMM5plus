module utility_interpolate

    use define_indexes
    use define_indexes
    use define_keys
    use define_globals
    use define_settings
    use utility_allocate

    implicit none

    public :: util_curve_lookup_singular
    public :: util_interpolate_linear

    private

contains
!%  
!%==========================================================================
!%==========================================================================
!%

    real(8) function util_interpolate_linear(X, X1, X2, Y1, Y2) result (Y)
    !% This is a linear interpolation function
        real(8), intent(in) :: X, X1, X2, Y1, Y2
        if (Y1 == Y2) then
            Y = Y1
            return
        end if
        Y = Y1 + (X - X1) * (Y2 - Y1) / (X2 - X1)
    end function util_interpolate_linear
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine util_curve_lookup_singular(curveID, er_inCol, er_outCol, &
        xVal_col, yVal_col)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% look up from a curve and returns a interpolated value
        !% Inputs/Outputs:
        !%      curveID
        !%      er_inCol
        !%      er_outCol
        !%      xVal_col
        !%      yVal_col     
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: curveID, er_inCol, er_outCol, xVal_col, yVal_col
        real(8), pointer    :: x, y
        integer, pointer    :: nRows, ElemIdx
        real(8) :: x1, x2, y1, y2, slope
        integer :: ii

        character(64) :: subroutine_name = 'util_curve_lookup_singular'
        !%-----------------------------------------------------------------------------
        !% pointers:
        ElemIdx => curve(curveID)%ElemIdx
        nRows   => curve(curveID)%NumRows
        x       => elemR(ElemIdx,er_inCol)
        y       => elemR(ElemIdx,er_outCol)

        x1 = curve(curveID)%ValueArray(1,xVal_col)
        y1 = curve(curveID)%ValueArray(1,yVal_col)

        !% check if inVal is smaller than the first table entry
        if (x <= x1) then
            if (x1 > zeroR) then
                y = x/x1*y1
                return
            else
                y = y1
                return
            end if
        end if

        !% else loop through the table to find the position and interpolate
        do ii = 2,nRows
            x2 = curve(curveID)%ValueArray(ii,xVal_col)
            y2 = curve(curveID)%ValueArray(ii,yVal_col)
            if (x2 /= x1) slope = (y2 - y1) / (x2 - x1)

            if (x <= x2) then
                y = util_interpolate_linear(x, x1, x2, y1, y2)
                return
            end if
            x1 = x2
            y1 = y2
            if (slope < zeroR) slope = zeroR 
            y = y1 + slope * (x - x1)
        end do

    end subroutine util_curve_lookup_singular
!%  
!%==========================================================================
!%==========================================================================
!%
end module utility_interpolate