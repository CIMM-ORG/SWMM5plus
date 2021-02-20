! module xsect_tables
!
! This module consists of tables of relative geometric properties for
! rounded cross-sections.
!
!==========================================================================
module xsect_tables

    ! use setting_definition
    ! use type_definitions
    use globals
    use setting_definition

    implicit none

    public


    !==========================================================================
    ! Circular Shape
    !==========================================================================
    !% Y/Yfull v. A/Afull
    integer, parameter      :: NYCirc = 51

    real(4), dimension(NYCirc) :: YCirc = (/0.0, 0.05236, 0.08369, 0.11025, 0.13423, 0.15643, 0.17755, 0.19772, 0.21704, 0.23581,&
        0.25412, 0.27194, 0.28948, 0.30653, 0.32349, 0.34017, 0.35666, 0.37298, 0.38915,      &
        0.40521, 0.42117, 0.43704, 0.45284, 0.46858, 0.48430, 0.50000, 0.51572, 0.53146,      &
        0.54723, 0.56305, 0.57892, 0.59487, 0.61093, 0.62710, 0.64342, 0.65991, 0.67659,      &
        0.69350, 0.71068, 0.72816, 0.74602, 0.76424, 0.78297, 0.80235, 0.82240, 0.84353,      &
        0.86563, 0.88970, 0.91444, 0.94749, 1.000 /)
    !% A/Afull v. Y/Yfull
    integer, parameter      :: NACirc = 51
    real(4), dimension(NACirc) :: ACirc = (/0.0, 0.00471, 0.0134, 0.024446, 0.0374, 0.05208, 0.0680, 0.08505, 0.1033, 0.12236,   &
        0.1423, 0.16310, 0.1845, 0.20665, 0.2292, 0.25236, 0.2759, 0.29985, 0.3242, 0.34874,  &
        0.3736, 0.39878, 0.4237, 0.44907, 0.4745, 0.5, 0.5255, 0.55093, 0.5763, 0.60135,      &
        0.6264, 0.65126, 0.6758, 0.70015, 0.7241, 0.74764, 0.7708, 0.79335, 0.8154, 0.8369,   &
        0.8576, 0.87764, 0.8967, 0.91495, 0.9320, 0.94792, 0.9626, 0.97555, 0.9866, 0.99516,  &
        1.000 /)
    !% R/Rfull v. Y/Yfull
    integer, parameter      :: NRCirc = 51
    real(4), dimension(NRCirc) :: RCirc = (/0.01, 0.0528, 0.1048, 0.1556, 0.2052, 0.2540, 0.3016, 0.3484, 0.3944, 0.4388, 0.4824,&
        0.5248, 0.5664, 0.6064, 0.6456, 0.6836, 0.7204, 0.7564, 0.7912, 0.8244, 0.8568,       &
        0.8880, 0.9176, 0.9464, 0.9736, 1.0000, 1.0240, 1.0480, 1.0700, 1.0912, 1.1100,       &
        1.1272, 1.1440, 1.1596, 1.1740, 1.1848, 1.1940, 1.2024, 1.2100, 1.2148, 1.2170,       &
        1.2172, 1.2150, 1.2104, 1.2030, 1.1920, 1.1780, 1.1584, 1.1320, 1.0940, 1.000 /)
    !%W/Wmax v. Y/Yfull
    integer, parameter      :: NWCirc = 51
    real(4), dimension(NWCirc) :: WCirc = (/0.0, 0.2800, 0.3919, 0.4750, 0.5426, 0.6000, 0.6499, 0.6940, 0.7332, 0.7684, 0.8000, &
        0.8285, 0.8542, 0.8773, 0.8980, 0.9165, 0.9330, 0.9474, 0.9600, 0.9708, 0.9798,       &
        0.9871, 0.9928, 0.9968, 0.9992, 1.0000, 0.9992, 0.9968, 0.9928, 0.9871, 0.9798,       &
        0.9708, 0.9600, 0.9474, 0.9330, 0.9165, 0.8980, 0.8773, 0.8542, 0.8285, 0.8000,       &
        0.7684, 0.7332, 0.6940, 0.6499, 0.6000, 0.5426, 0.4750, 0.3919, 0.2800, 0.0 /)
contains

    !
    !==========================================================================
    ! functions for table interpolation
    !==========================================================================
    !
    pure function table_lookup &
        (normalizedInput, table, nItems) result(normalizedOutput)
        !
        ! took up talbe values of circular types of geometry
        !
        real(4),      intent(in)      :: table(:)
        real(4),      intent(in)      :: normalizedInput
        integer,   intent(in)      :: nItems

        real(4)    :: normalizedOutput, normalizedOutput2
        real(4)    :: delta, startPos, endPos
        integer  :: ii

        !--------------------------------------------------------------------------
        !% find which segment of table contains x
        delta = oneR / (nItems - oneR)

        ii = int(normalizedInput / delta)

        if     ( ii .GE. (nItems - oneI) ) then

            normalizedOutput = table(nItems - 1)

        elseif ( ii .LE. zeroI) then

            normalizedOutput = zeroR

        else

            startPos = ii * delta
            endPos   = (ii + oneI) * delta

            normalizedOutput = table(ii) + (normalizedInput - startPos) * &
                (table(ii + oneI) - table(ii)) / delta

            if (ii == oneI) then
                ! use quadratic interpolation for low x value
                normalizedOutput2 = normalizedOutput + (normalizedInput - startPos) &
                    * (normalizedInput - endPos) / (delta*delta) * (table(ii)/2.0 - table(ii+1) &
                    + table(ii+2)/2.0)

                if ( normalizedOutput2 > 0.0 ) then
                    normalizedOutput = normalizedOutput2
                endif

            endif

        endif

    end function table_lookup
    !
    !==========================================================================
    !==========================================================================
    !
    pure function get_theta_of_alpha &
        (alpha) result(theta)
        !
        ! get the angle theta for small value of A/Afull (alpha) for circular geometry
        !
        real(4),      intent(in)      :: alpha

        real(4)    :: theta
        real(4)    :: theta1, d, ap
        integer  :: ii

        !--------------------------------------------------------------------------
        !% this code is adapted from SWMM 5.1 source code
        if     (alpha .GE. 1.0) then
            theta = 1.0
        elseif (alpha .LE. 0.0) then
            theta = 0.0
        elseif (alpha .LE. 1.0e-5) then
            theta = 37.6911 / 16.0 * alpha ** (onethirdR)
        else
            theta = 0.031715 - 12.79384 * alpha + 8.28479 * sqrt(alpha)
            theta1 = theta
            ap = twoR * pi *alpha
            do ii = 1,40
                d = - (ap - theta + sin(theta)) / (1.0 - cos(theta))
                if (d > 1.0) then
                    d = sign(real(1.0),d)
                endif
                theta = theta - d
                if ( abs(d) .LE. 0.0001 ) then
                    return
                endif
            enddo
            theta = theta1
            return
        endif

    end function get_theta_of_alpha
    !==========================================================================
    ! END OF MODULE xsect_tables
    !==========================================================================
end module xsect_tables
