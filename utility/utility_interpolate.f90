module utility_interpolate

    use define_indexes
    use define_indexes
    use define_keys
    use define_globals
    use define_settings
    use utility_allocate
    use utility_crash, only: util_crashpoint

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
        ! print *, 'X-X1  ',X-X1
        ! print *, 'Y2-Y1 ', Y2-Y1
        ! print *, 'X2-X1 ',X2-X1
        ! print *, 'Y1    ',Y1
        ! print *, 'Y     ',Y
    end function util_interpolate_linear
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine util_curve_lookup_singular(&
            curveID, er_inCol, er_outCol,  xVal_col, yVal_col, interpFlag)
        !%------------------------------------------------------------------
        !% Description:
        !% look up from a curve and returns a interpolated value
        !% Inputs/Outputs:
        !%      curveID   -- CurveID for this element
        !%      er_inCol  -- column in elemR(:,er_inCol) with curve X data
        !%      er_outCol -- column in elemR(:,er_outCol) for curve Y data
        !%      xVal_col  -- x data column in the curve ValueArray
        !%      yVal_col  -- y data column in the curve Valuearray
        !%      interpFlag = 0 for step wise, 1 for linear interp
        !%--------------------------------------------------------------------
            integer, intent(in) :: curveID, er_inCol, er_outCol, xVal_col, yVal_col
            integer, intent(in) :: interpFlag
            real(8), pointer    :: x, y, x1, x2, y1, y2
            integer, pointer    :: nRows, ElemIdx
            real(8) :: slope
            integer :: ii

            character(64) :: subroutine_name = 'util_curve_lookup_singular'
        !%---------------------------------------------------------------------
        !% Aliases
            ElemIdx => curve(curveID)%ElemIdx
            nRows   => curve(curveID)%NumRows
            x       => elemR(ElemIdx,er_inCol)  !% this drives the curve
            y       => elemR(ElemIdx,er_outCol) !% this is the curve output
        !%---------------------------------------------------------------------

        !% --- starting position of the curve (these pointers change)
        x1 => curve(curveID)%ValueArray(1,xVal_col)
        y1 => curve(curveID)%ValueArray(1,yVal_col)

        !% --- maximum x, y values (these pointers change)
        x2 => curve(curveID)%ValueArray(nRows,xVal_col)
        y2 => curve(curveID)%ValueArray(nRows,yVal_col)

        !% --- parse for where x falls in range (x1,x2)
        if (x < x1) then
            !% --- if x input is smaller than lower bound
            if (x1 > zeroR) then
                !% --- We consider (0,0) to be the last point in the curve
                select case (interpFlag)
                case (0)
                    !% --- stepwise:use lowest value
                    y = y1
                case (1)
                    !% --- interpolation: set output by interpolating towards zero
                    y = (x/x1)*y1
                case default
                    print *, 'CODE ERROR: unexpected case default for interpFlag for ',trim(subroutine_name)
                    call util_crashpoint(698743)
                end select
                return !% output found and stored in y
            else
                !% --- if x1 < 0, there's no valid interpolation below x1
                !%     both stepwise and interpolation return the y1 value
                y = y1
                return !% output found and stored in y
            end if 
        elseif (x == x1) then
            !% --- x at lower bound
            y = y1
            return !% output found and stored in y
        elseif (x == x2) then
            !% --- x at upper bound
            y = y2
            return !% output found and stored in y
        elseif (x > x2) then
            !% --- if x is larger than upper bound
            select case (interpFlag)
            case(0)
                !% --- stepwise: return the largest y2 value
                y = y2
                return !% output found and stored in y
            case(1)
                !% --- linear: extrapolate based on slope at end of curve
                x1 => curve(curveID)%ValueArray(nRows-1,xVal_col)
                y1 => curve(curveID)%ValueArray(nRows-1,yVal_col)
                !% --- compute the slope for extrapolation
                if (x2 /= x1) then
                    slope = (y2 - y1) / (x2 - x1)
                    if (slope < zeroR) slope = zeroR !% -- don't allow negative slope extrapolation
                    !% --- interpolate: extrapolate using slope
                    y = y2 + slope * (x - x2)
                    return !% output found and stored in y
                else
                    !% --- last two points are identical
                    y = y2
                    return !% output found and stored in y
                end if
            case default
                print *, 'CODE ERROR: Unexpected case default in interpFlag for ',trim(subroutine_name)
                    call util_crashpoint(327768)
            end select
        else
            !% --- x is between lower and upper bound of curve
            !%     simple (but inefficient) step through curve
            do ii = 2,nRows
                !% --- get the next valuues from the curve
                x2 => curve(curveID)%ValueArray(ii,xVal_col)
                y2 => curve(curveID)%ValueArray(ii,yVal_col)
                if (x <= x2) then
                    select case (interpFlag)
                    case (0)
                        !% --- stepwise value
                        y = y1
                        return !% answer found and stored in y, so exit subroutine
                    case (1)
                        !% --- linear interpolation
                        y = util_interpolate_linear(x, x1, x2, y1, y2)
                        return !% answer found and stored in y, so exit subroutine
                    case default
                        print *, 'CODE ERROR: Unexpected case default in interpFlag for ',trim(subroutine_name)
                        call util_crashpoint(6682093)
                    end select
                end if
                !% --- unsuccessful, increment the x1,y1 position for next loop
                x1 => curve(curveID)%ValueArray(ii,xVal_col)
                y1 => curve(curveID)%ValueArray(ii,yVal_col)
            end do
            !% --- code should never reach here because we've parsed for x >= upper bound.
            print *, 'CODE ERROR: unexpected result -- reached end of loop without finding answer'
            call util_crashpoint(698743)
        end if


    end subroutine util_curve_lookup_singular
!%  
!%==========================================================================
!%==========================================================================
!%
end module utility_interpolate