! module xsect_tables
!
! This module consists of tables of relative geometric properties for
! rounded cross-sections.
!
!==========================================================================
module xsect_tables

    use define_indexes
    use define_globals
    use define_settings
    use define_xsect_tables

    implicit none

    public

contains

!
!==========================================================================
! functions for table interpolation
!==========================================================================
!
    ! pure function table_lookup &
    !     (normalizedInput, table, nItems) result(normalizedOutput)
        !
        ! table lookup function. This function is single operation
        !
        ! real(8),      intent(in)      :: table(:)
        ! real(8),      intent(in)      :: normalizedInput
        ! integer,   intent(in)      :: nItems

        ! real(8)     :: normalizedOutput, normalizedOutput2
        ! real(8)     :: delta, startPos, endPos
        ! integer  :: ii

        ! !--------------------------------------------------------------------------
        ! !% find which segment of table contains x
        ! delta = oneR / (nItems - oneR)

        ! ii = int(normalizedInput / delta)

        ! if     ( ii .GE. (nItems - oneI) ) then

        !     normalizedOutput = table(nItems)

        ! elseif ( ii .LE. zeroI) then

        !     normalizedOutput = zeroR

        ! else

        !     startPos = ii * delta
        !     endPos   = (ii + oneI) * delta

        !     normalizedOutput = table(ii) + (normalizedInput - startPos) * &
        !         (table(ii + oneI) - table(ii)) / delta

        !     if (ii == oneI) then
        !         ! use quadratic interpolation for low x value
        !         normalizedOutput2 = normalizedOutput + (normalizedInput - startPos) &
        !             * (normalizedInput - endPos) / (delta*delta) * (table(ii)/2.0 - table(ii+1) &
        !             + table(ii+2)/2.0)

        !         if ( normalizedOutput2 > 0.0 ) then
        !             normalizedOutput = normalizedOutput2
        !         endif

        !     endif

        ! endif

    ! end function table_lookup
!
!==========================================================================
!==========================================================================
!
    ! pure function get_theta_of_alpha &
    !     (alpha) result(theta)
        !
        ! get the angle theta for small value of A/Afull (alpha) for circular geometry
        !
        ! real(8),      intent(in)      :: alpha


        ! real(8)     :: theta
        ! real(8)     :: theta1, d, ap

        ! integer  :: ii

        ! !--------------------------------------------------------------------------
        ! !% this code is adapted from SWMM 5.1 source code
        ! if     (alpha .GE. 1.0) then
        !     theta = 1.0
        ! elseif (alpha .LE. 0.0) then
        !     theta = 0.0
        ! elseif (alpha .LE. 1.0e-5) then
        !     theta = 37.6911 / 16.0 * alpha ** (onethirdR)
        ! else
        !     theta = 0.031715 - 12.79384 * alpha + 8.28479 * sqrt(alpha)
        !     theta1 = theta
        !     ap = twoR * pi *alpha
        !     do ii = 1,40
        !         d = - (ap - theta + sin(theta)) / (1.0 - cos(theta))
        !         if (d > 1.0) then
        !             d = sign(oneR,d)
        !         endif
        !         theta = theta - d
        !         if ( abs(d) .LE. 0.0001 ) then
        !             return
        !         endif
        !     enddo
        !     theta = theta1
        !     return
        ! endif

    ! end function get_theta_of_alpha
!
!==========================================================================
!==========================================================================
!
    subroutine xsect_table_lookup &
        (inoutArray, normalizedInput, table, nItems, thisP)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% interpolates the normalized vaule from the lookup table.
        !% this subroutine is vectorized array operation.
        !%-----------------------------------------------------------------------------
        real(8), intent(inout)    :: inoutArray(:)
        real(8), intent(in)       :: normalizedInput(:), table(:)
        integer, intent(in)       :: nItems, thisP(:)
        integer, pointer          :: position(:)
        real(8)                   :: delta
        !%-----------------------------------------------------------------------------
        if (icrash) return

        !% pointer towards the position in the lookup table
        !% this is pointed towards temporary column
        position => elemI(:,ei_Temp01)

        delta = oneR / (nItems - oneR)

        !% this finds the position in the table for interpolation
        position(thisP) = int(normalizedInput(thisP) / delta)

        !% find the normalized output from the lookup table
        where (position(thisP) .LE. zeroI)
            inoutArray(thisP) = zeroR

        elsewhere ( (position(thisP) .GT. zeroI          ) .and. &
                    (position(thisP) .LT. (nItems - oneI)) )

            !%  Y = Y_a + (Y_b-Y_a)*(X_0-X_a)/(X_b-X_a)
            inoutArray(thisP) = table(position(thisP)+oneI) + &
                                (normalizedInput(thisP) - position(thisP) * delta) * &
                                (table(position(thisP) + twoI) - table(position(thisP)+oneI)) / delta

        elsewhere (position(thisP) .GE. (nItems - oneI))
            inoutArray(thisP) = table(nItems)
        endwhere

        !% quadratic interpolation for low value of normalizedInput
        where (position(thisP) .LT. twoI)
            inoutArray(thisP) = max(zeroR, &
                    (inoutArray(thisP) + (inoutArray(thisP) - delta) * &
                    (inoutArray(thisP) - twoI * delta) / (delta*delta) *         &
                    (table(oneI)/twoR - table(twoI)  +  table(threeI) / twoR)) )
        endwhere

        !% reset the temporary values to nullvalue
        position(thisP) = nullvalueI

    end subroutine xsect_table_lookup
!
!==========================================================================
!==========================================================================
!
    real(8) function xsect_table_lookup_singular &
        (normalizedInput, table, nItems) result (normalizedOutput)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% interpolates the normalized vaule from the lookup table.
        !% this function is singular operation.
        !%-----------------------------------------------------------------------------
        real(8), intent(in)       :: normalizedInput, table(:)
        integer, intent(in)       :: nItems
        integer                   :: position
        real(8)                   :: delta
        !%-----------------------------------------------------------------------------
        if (icrash) return

        delta = oneR / (nItems - oneR)

        !% this finds the position in the table for interpolation
        position = int(normalizedInput / delta)

        !% find the normalized output from the lookup table
        if (position .LE. zeroI) then
            normalizedOutput = zeroR

        else if ( (position .GT. zeroI          ) .and. &
                  (position .LT. (nItems - oneI)) ) then

            !%  Y = Y_a + (Y_b-Y_a)*(X_0-X_a)/(X_b-X_a)
            normalizedOutput = table(position+oneI) + (normalizedInput - position * delta) * &
                            (table(position+twoI) - table(position+oneI)) / delta

        else if (position .GE. (nItems - oneI)) then
            normalizedOutput = table(nItems)
        end if

        !% quadratic interpolation for low value of normalizedInput
        if (position .LT. twoI) then
            normalizedOutput = max(zeroR, &
                    (normalizedOutput + (normalizedOutput - delta) * &
                    (normalizedOutput - twoI * delta) / (delta*delta) *         &
                    (table(oneI)/twoR - table(twoI)  +  table(threeI) / twoR)) )
        end if

    end function xsect_table_lookup_singular
!%
!%==========================================================================
! END OF MODULE xsect_tables
!%==========================================================================
!%
end module xsect_tables
