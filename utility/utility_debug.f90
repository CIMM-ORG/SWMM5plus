module utility_debug

  use define_indexes
  use define_keys
  use define_globals
  use define_settings
  use, intrinsic :: iso_fortran_env, only: error_unit

  implicit none

  private

  public :: debug_2D_array_csv
  public :: debug_Nface_check

contains

  subroutine debug_2D_array_csv(file_name_input, type, header, arr_real, arr_int, arr_log)

    character(64) :: subroutine_name = 'debug_2D_array_csv'
    character(len = *), intent(in) :: file_name_input !% file name wanted
    character(len = 1), intent(in) :: type !% this should be I or i for int, R or r for real and L or l for log
    character(len = *), optional, intent(in) :: header !% header that is printed at the top of the file
    real(8), optional, intent(in) :: arr_real(:,:)
    integer, optional, intent(in) :: arr_int(:,:)
    logical, optional, intent(in) :: arr_log(:,:)
    character(len = :), allocatable :: file_name
    character(len = 5) :: str_image
    integer :: ii, jj, fu, rc, image

    if (setting%Debug%File%initialization) print *, '*** enter ', this_image(), subroutine_name


    !% fu stands for file unit which will be tied to the image and tells the system what file to open
    fu = this_image()

    !& here we create the file name that will be opened by the code
    write(str_image, '(i5.5)') fu
    file_name = 'debug/'//trim(file_name_input)//'_'//trim(str_image)//'.csv'

    !% opening the file, as well as error handing if the open fails
    open (action='write', file=file_name, status='replace', iostat=rc, newunit=fu)
    if (rc .ne. 0) then
       write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', rc
       stop "in " // subroutine_name
   end if

    write(fu,'(A)', advance = "no") header
    write(Fu, *)


    !% We write to the file depending on the array type selected.

    if (type == 'R' .or. type == 'r') then

       do ii = 1, size(arr_real(1,:))
          do jj = 1, size(arr_real(:,1))
             write (fu,'(F40.20)',advance = "no") arr_real(jj,ii)
             write (fu,'(A2)',advance = "no") ','
          end do
          write (fu, *)
       end do

    else if (type == 'I' .or. type == 'i') then

       do ii = 1, size(arr_int(1,:))
          do jj = 1, size(arr_int(:,1))
             write (fu,'(I10)',advance = "no") arr_int(jj,ii)
             write (fu,'(A2)',advance = "no") ','
          end do
          write (fu, *)
       end do


    else if (type == 'L' .or. type == 'l') then

       do ii = 1, size(arr_log(1,:))
          do jj = 1, size(arr_log(:,1))
             write (fu,'(L4)',advance = "no") arr_log(jj,ii)
             write (fu,'(A2)',advance = "no") ','
          end do
          write (fu, *)
       end do

    else

       print *, "INCORRECT TYPE INPUT, EMPTY CSV FILE CREATED"

    end if

    !% closing the file
    close(fu)


    if (setting%Debug%File%initialization)  print *, '*** leave ', this_image(), subroutine_name

  end subroutine debug_2D_array_csv



  subroutine debug_Nface_check
    !% debug to check that the number of faces on each processor is equal to N_Face

    integer :: ii, total_faces

    character(64) :: subroutine_name = 'debug_Nface_check'
    if (setting%Debug%File%initialization) print *, '*** enter ', this_image(), subroutine_name


    total_faces = 0

    !% Looping through elemI, and finding the last local index
    do ii=1, size(elemI(:,ei_Lidx))

       if (elemI(ii,ei_Mface_uL) /= nullvalueI .and.  total_faces < elemI(ii,ei_Mface_uL) ) then

          total_faces = elemI(ii, ei_Mface_uL)

       end if

       if (elemI(ii,ei_Mface_dL) /= nullvalueI .and. total_faces < elemI(ii,ei_Mface_dL)) then

          total_faces = elemI(ii, ei_Mface_dL)

       end if

    end do

    !% Then we compare it to N_Face and print if it is the same or not.

    if (total_faces == N_face(this_image())) then

       print *, "CORRECT NUMBER OF FACES ON IMAGE ::", this_image()

    else

       print *, "ERROR, INCORRECT NUMVER OF FACES ON IMAGE ::", this_image()
       print *, "total_faces = ", total_faces
       print *, "N_face(this_image()) =", N_face(this_image())
    end if

    if (setting%Debug%File%initialization)  print *, '*** leave ', this_image(), subroutine_name

  end subroutine debug_Nface_check

end module utility_debug
