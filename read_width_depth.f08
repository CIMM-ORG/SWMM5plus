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
    public  :: split_read

    integer :: debuglevel = 0
    
 contains
!
!========================================================================== 
!==========================================================================
!
 subroutine read_widthdepth_pairs &
    (iunit, ID, numberPairs, ManningsN, Length, zBottom, xDistance, &
     Breadth, widthDepthData, cellType)
 
 character(64) :: subroutine_name = 'read_widthdepth_pairs'
 
 integer, intent(in)  :: iunit   ! the file unit number
 
 integer :: number_of_cells
 integer :: max_number_of_pairs
 
 logical :: next_data_is_pair  = .false.
 logical :: dont_quit          = .true.
 logical :: expecting_new_cell = .true.
 integer :: icell = 0
 integer :: ipair = 0
 character(len=256)  :: tmp
 integer             :: istat
 
 
! These should be organized later
 integer, dimension(:), allocatable :: ID
 integer, dimension(:), allocatable :: numberPairs
 
 real, dimension(:), allocatable :: ManningsN
 real, dimension(:), allocatable :: Length
 real, dimension(:), allocatable :: zBottom
 real, dimension(:), allocatable :: xDistance
 real, dimension(:), allocatable :: Breadth
 
 real, dimension(:,:,:), allocatable :: widthDepthData
 
 character(len=:), allocatable :: cellType(:)
 
 character(len=256) :: value1
 character(len=256) :: value2
 
 
 integer :: allocation_status
 character(len=99) :: emsg
 real :: tmpID
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 ! error checking
 number_of_cells = read_number_of_cells(iunit)
 
 max_number_of_pairs = read_max_number_of_pairs (iunit)
 
 allocate(ID(number_of_cells), stat=allocation_status, errmsg=emsg)
 ID(:) = 0
 
 allocate(numberPairs(number_of_cells), stat=allocation_status, errmsg=emsg)
 numberPairs(:) = 0
 
 allocate(ManningsN(number_of_cells), stat=allocation_status, errmsg=emsg)
 ManningsN(:) = 0.0
 
 allocate(Length(number_of_cells), stat=allocation_status, errmsg=emsg)
 Length(:) = 0.0
 
 allocate(zBottom(number_of_cells), stat=allocation_status, errmsg=emsg)
 zBottom(:) = 0.0
 
 allocate(xDistance(number_of_cells), stat=allocation_status, errmsg=emsg)
 xDistance(:) = 0.0
 
 allocate(Breadth(number_of_cells), stat=allocation_status, errmsg=emsg)
 Breadth(:) = 0.0
 !+1 is a ghost cell for the chack that the width-depth pairs cover enough
 !depth and fixing of vertical walls
 allocate(widthDepthData(number_of_cells, max_number_of_pairs+1, wd_idx_max), stat=allocation_status, errmsg=emsg)
 widthDepthData(:,:,:) = 0.0
 
 allocate(character(100):: cellType(number_of_cells), stat=allocation_status, errmsg=emsg)
 
 rewind(iunit)
 do while (dont_quit .eqv. .true.)
    read(iunit,fmt='(A)',iostat=istat) tmp
    if (is_iostat_end(istat)) exit
    
    tmp = adjustl(tmp)
    if (tmp(1:1) == '#') then
        cycle
    endif
    
    tmp = trim(tmp)
    call split_read_two_strings (tmp, ' ', value1, value2)
    
    if (value1 == 'begin') then
        if(value2 == 'cell') then
            if (expecting_new_cell .neqv. .true.) then
                print*, tmp
                print*, 'error, likely misalignment in file'
                stop
            else
                !reset for new cell
                expecting_new_cell = .false.
                next_data_is_pair  = .false.
                icell = icell + 1
            endif
        else
            
            print*, tmp
            print*, 'error, not yet designed for features other than cells'
            stop
        endif
    endif
    
    call split_read_two_strings (tmp, ' ', value1, value2)
    
    if (value1 == 'ID') then
        read(value2 , *, iostat=istat) tmpID
        ID(icell) = int(tmpID)
    elseif (value1 == 'ManningsN') then
        read(value2 , *, iostat=istat) ManningsN(icell)
    elseif (value1 == 'Length') then
        read(value2 , *, iostat=istat) Length(icell)
    elseif (value1 == 'zBottom') then
        read(value2 , *, iostat=istat) zBottom(icell)
    elseif (value1 == 'xDistance') then
        read(value2 , *, iostat=istat) xDistance(icell)
    elseif (value1 == 'Breadth') then
        read(value2 , *, iostat=istat) Breadth(icell)
    elseif (value1 == 'cellType') then
        if (value2 == 'channel_WidthDepthPairs') then
            cellType(icell) = 'widthdepth_pair'
        else
            print*,'error: unknown value for cellType'
        endif
    elseif (value1 == 'numberPairs') then
        read(value2 , *, iostat=istat) numberPairs(icell)
    elseif (value1 == 'WidthDepthPairs') then
        if (value2 == 'follow') then
            next_data_is_pair = .true.
            ipair = 0
        else
            print*, 'error, not designed for 2nd argument other than follow'
        endif
    elseif (value1 == 'end') then
        if(value2 == 'WidthDepthPairs') then
            continue
        elseif(value2 == 'cell') then
            expecting_new_cell = .true.
        else
            print*, tmp
            print*, value2
            print*, 'error, unknown option'
            stop
        endif
    else
        if(ipair == 0) then 
            ipair = 1
        endif
        
        if(value1 == 'end') then
            next_data_is_pair = .false.
        else
            read(value1 , *, iostat=istat) widthDepthData(icell, ipair, wd_widthAtLayerTop)
            read(value2 , *, iostat=istat) widthDepthData(icell, ipair, wd_depthAtLayerTop)
            ipair = ipair+1
        endif
    endif
 enddo
 rewind(iunit)
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine read_widthdepth_pairs
!
!========================================================================== 
!==========================================================================
!
 subroutine split_read_two_strings (str, sep, value1, value2)
 
 character(64) :: subroutine_name = 'split_read_two_strings'
 
 character(len=*), intent(in)     :: str
 character(len=*), intent(in)     :: sep
 character(len=256), intent(inout) :: value1
 character(len=256), intent(inout) :: value2
 integer :: i,n
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 read (unit=str,fmt=*) value1, value2
 
 value1 = trim(value1)
 value2 = trim(value2)
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine split_read_two_strings
!
!========================================================================== 
!==========================================================================
!
 subroutine split_read (str, sep, pairValues)
 
 character(64) :: subroutine_name = 'split_read'
 
 character(len=*), intent(in)     :: str
 character(len=*), intent(in)     :: sep
 character(len=:), allocatable :: pairValues(:)
 integer :: i,n
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 n = 1
 do i=1, len(str)
     if (str(i:i) == sep) then
         n = n + 1
     endif
 end do
 allocate (character(256):: pairValues(n))
 read (unit=str,fmt=*) pairValues
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine split_read
!
!========================================================================== 
!==========================================================================
!
 function read_number_of_cells (iunit) result(n_lines)
 
 character(64) :: subroutine_name = 'read_number_of_cells'
 
 integer,intent(in)  :: iunit   ! the file unit number
 integer             :: n_lines ! the number of lines in the file
 character(len=256)  :: tmp
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
 character(len=256)  :: tmp
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
