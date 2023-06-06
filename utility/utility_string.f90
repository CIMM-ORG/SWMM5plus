module utility_string
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% string handling procedures
    !%==========================================================================

    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    elemental subroutine util_lower_case(word)
    !%------------------------------------------------------------------
    !% Description:
    !%   Convert a word to lower case
    !%------------------------------------------------------------------
    !% Declarations:
        character (len=*) , intent(in out) :: word
        integer                            :: ii, ic
    !%------------------------------------------------------------------

        do ii = 1, len(word)
            ic = ichar(word(ii:ii))
            if (ic >= 65 .and. ic < 90) word(ii:ii) = char(ic+32)
        end do

    end subroutine util_lower_case

end module utility_string
