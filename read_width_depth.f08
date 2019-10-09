!==========================================================================
! Reads a set of widthdepth pairs from a text file for each cell
! Text file is organized with in the form of "keyword value" with
! the width depth pairs as two items per line.  Typical layout and 
! allowable keyword value pairs are as follow (where xxxx is a number)
!     begin cell
!     ID  xxxx
!     Length xxxx
!     xDistance xxxx
!     zBottom xxxx
!     Breadth xxxx
!     TrapezoidAngle xxxx (in degrees)
!     cellType channel_WidthDepthPairs
!     numberPairs xxxx
!     WidthDepthPairs follow
!     xxxx xxxx
!     xxxx xxxx
!     end WidthDepthPairs
!     end cell
!     begin cell
!     ....
!     end cell
!==========================================================================
!
 module read_width_depth
! 
    use array_index
    use data_keys
    use globals
    use setting_definition

    
    implicit none
    
    private
    
    public  :: read_number_of_cells
    public  :: read_max_number_of_pairs
    public  :: read_widthdepth_pairs
    public  :: split_read_int

    integer :: debuglevel = 0
    
 contains
!
!========================================================================== 
!==========================================================================
!
 subroutine read_widthdepth_pairs &
    (NX, iunit, number_of_cells)
 
 character(64) :: subroutine_name = 'read_widthdepth_pairs'
 
 integer, intent(in)  :: iunit   ! the file unit number
 
 integer, intent(in)  :: NX
 integer, intent(out) :: number_of_cells
 
 logical :: next_data_is_pair  = .false.
 logical :: dont_quit          = .true.
 logical :: expecting_new_cell = .true.
 integer :: icell = -1
 
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 ! error checking
 number_of_cells = read_number_of_cells(iunit)
 
 do while (dont_quit .eqv. .true.)
 
 enddo
 
 
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine read_widthdepth_pairs
!
!========================================================================== 
!==========================================================================
!
 subroutine split_read_int (str, sep, pairValues)
 
 character(64) :: subroutine_name = 'split_read_int'
 
 character(len=*), intent(in)     :: str
 character(len=*), intent(in)     :: sep
 real, dimension(:), allocatable  :: pairValues
 integer :: i,n
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 n = 1
 do i=1, len(str)
     if (str(i:i) == sep) then
         n = n + 1
     endif
 end do
 allocate (pairValues(n))
 read (unit=str,fmt=*) pairValues
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine split_read_int
!
!========================================================================== 
!==========================================================================
!
 function read_number_of_cells (iunit) result(n_lines)
 
 character(64) :: subroutine_name = 'read_number_of_cells'
 
 integer,intent(in)  :: iunit   ! the file unit number
 integer             :: n_lines ! the number of lines in the file
 character(len=110)  :: tmp
 integer             :: istat
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 rewind(iunit)
 n_lines = 0
 do
    read(iunit,fmt='(A)',iostat=istat) tmp
    if (is_iostat_end(istat)) exit

    if (tmp == 'begin cell') then
        n_lines = n_lines + 1
    endif
 end do
 rewind(iunit)
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end function read_number_of_cells
!
!========================================================================== 
!==========================================================================
!
 function read_max_number_of_pairs (iunit) result(n_pairs)
 
 character(64) :: subroutine_name = 'read_max_number_of_pairs'
 
 integer,intent(in)  :: iunit   ! the file unit number
 integer             :: n_pairs ! the number of lines in the file
 character(len=110)  :: tmp
 integer             :: istat
 integer             :: tmpPair
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 rewind(iunit)
 n_pairs = 0
 do
    read(iunit,fmt='(A)',iostat=istat) tmp
    if (is_iostat_end(istat)) exit

    if (tmp(1:11) == 'numberPairs') then
        read(tmp(13:15),'(I3)') tmpPair
        n_pairs = max(n_pairs, tmpPair)
    endif
 end do
 rewind(iunit)
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end function read_max_number_of_pairs
!
!========================================================================== 
! END OF MODULE stub
!==========================================================================
 end module read_width_depth