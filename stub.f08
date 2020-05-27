!==========================================================================
!
module stub
   !
   use array_index
   use data_keys
   use globals
   use setting_definition


   implicit none

   private

   integer :: debuglevel = 0

contains
   !
   !==========================================================================
   !==========================================================================
   !
   subroutine x &
      ()

      character(64) :: subroutine_name = 'x'


      !--------------------------------------------------------------------------
      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

      if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
   end subroutine x
   !
   !==========================================================================
   ! END OF MODULE stub
   !==========================================================================
end module stub
