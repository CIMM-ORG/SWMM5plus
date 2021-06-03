module utility_string

contains

    elemental subroutine utility_lower_case(word)
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   Convert a word to lower case
    !
    !-----------------------------------------------------------------------------

        character (len=*) , intent(in out) :: word
        integer                            :: ii, ic

    !-----------------------------------------------------------------------------

        do ii = 1, len(word)
            ic = ichar(word(ii:ii))
            if (ic >= 65 .and. ic < 90) word(ii:ii) = char(ic+32)
        end do

    end subroutine utility_lower_case

end module utility_string
