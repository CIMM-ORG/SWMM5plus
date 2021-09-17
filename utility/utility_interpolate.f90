module utility_interpolate

    use define_indexes
    use define_indexes
    use define_keys
    use define_globals
    use define_settings
    use utility_allocate

    implicit none

    public

contains

    real(8) function util_interpolate_linear(X, X1, X2, Y1, Y2) result (Y)
    !% This is a linear interpolation function
        real(8), intent(in) :: X, X1, X2, Y1, Y2
        if (Y1 == Y2) then
            Y = Y1
            return
        end if
        Y = Y1 + (X - X1) * (Y2 - Y1) / (X2 - X1)
    end function util_interpolate_linear

end module utility_interpolate