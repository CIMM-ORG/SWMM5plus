module utility_deallocate

    use define_indexes
    use define_globals
    use define_settings, only: setting

    implicit none

!-----------------------------------------------------------------------------
!
! Description:
!   This has all the array and variables deallocation.
!   All the top-level storage should be deallocated in this module
!
!-----------------------------------------------------------------------------

    private

    integer           :: deallocation_status
    character(len=99) ::              emsg

    public :: util_deallocate_network_data
    public :: util_deallocate_partitioning_arrays
    public :: util_deallocate_check

contains

    subroutine util_deallocate_network_data()
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   deallocate all the data relate to defined network, including:
    !   linkX, nodeX, elemX, faceX
    !-----------------------------------------------------------------------------

        character(64) :: subroutine_name = 'util_deallocate_network_data'

    !-----------------------------------------------------------------------------
        if (setting%Debug%File%utility_deallocate) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        call util_deallocate_linknode()

        call util_deallocate_elemX_faceX()

        !stop 7504
        call util_deallocate_columns()

        call util_deallocate_bc()

        if (setting%Debug%File%utility_deallocate) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine util_deallocate_network_data
!
!==========================================================================
!==========================================================================
!
    subroutine util_deallocate_linknode()
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   deallocates the link and node storage used for the coarse representation
    !   of the network connectivity
    !
    ! Method:
    !   Release the memory used by the tables node%I, link%I, node%R, link%R, node%YN,
    !   link%YN.
    !   These are defined in globals.f08), and allocated in allocate_storage.f08
    !   Every time memory is deallocated, the utility_check_deallocation functionality
    !   (from utility.f08) is used to determine wheter or not there was an error during
    !   the deallocation.
    !
    !-----------------------------------------------------------------------------

        integer       :: ii
        character(64) :: subroutine_name = 'util_deallocate_linknode'

    !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_deallocate) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        deallocate(node%I, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'node%I')

        deallocate(link%I, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'link%I')

        deallocate(node%R, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'node%R')

        deallocate(node%YN, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'node%YN')

        if (allocated(node%P%have_flowBC)) then
            deallocate(node%P%have_flowBC, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'node%P%Have_flowBC')
        end if

        if (allocated(node%P%have_headBC)) then
            deallocate(node%P%have_headBC, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'node%P%have_headBC')
        end if

        !% Deallocate link/node names
        do ii = 1, size(link%Names)
            deallocate(link%Names(ii)%str, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'link%Names(ii)%str')
        end do
        do ii = 1, size(node%Names)
            deallocate(node%Names(ii)%str, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'node%Names(ii)%str')
        end do

        deallocate(node%Names, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'node%Names')

        deallocate(link%Names, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'link%Names')

        deallocate(node_output_idx,stat=deallocation_status,errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'node_output_idx')

        deallocate(link_output_idx,stat=deallocation_status,errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'link_output_idx')

        if (allocated(node%P%have_output)) then
            deallocate(node%P%have_output,stat=deallocation_status,errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'node%P%have_output')
        end if

        if (allocated(link%P%have_output)) then
            deallocate(link%P%have_output,stat=deallocation_status,errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'link%P%have_output')
        end if


        if (setting%Debug%File%utility_deallocate) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

    end subroutine util_deallocate_linknode

!
!==========================================================================
!==========================================================================
!

    subroutine util_deallocate_partitioning_arrays()
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   deallocates the partitioning arrays storage used for storing partitioned
    !   network information
    !
    !-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'util_deallocate_partitioning_arrays'
    !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_deallocate) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        deallocate(adjacent_links, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg,'adjacent_links')

        deallocate(elem_per_image, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elem_per_image')

        deallocate(image_full, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'image_full')

        if (setting%Debug%File%utility_deallocate) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine util_deallocate_partitioning_arrays

!
!==========================================================================
!==========================================================================
!


    subroutine util_deallocate_elemX_faceX ()
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   Simply deallocate the elemX and faceX thatwe assigned across all employed images
    !
    !-----------------------------------------------------------------------------

        character(64) :: subroutine_name = 'util_deallocate_elemX_faceX'

    !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_deallocate) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"


        !==== elem deallocation ====
        deallocate(elemR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elemR')

        deallocate(elemI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elemI')

        deallocate(elemYN, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elemYN')

         !% OK HERE

        deallocate(elemP, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elemP')

        deallocate(elemPGetm, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elemPGetm')

        deallocate(elemPGac, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elemPGac')


        deallocate(elemPGalltm, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elemPGalltm')

        deallocate(elemSI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elemSI')

        !% oK here

        !print *, elemSR(:,:)
        !deallocate(elemSR, stat=deallocation_status)
        !print *, deallocation_status
        !stop 76896

        deallocate(elemSR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elemSR')



        !==== face deallocation ====
        deallocate(faceR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'faceR')

        deallocate(faceI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'faceI')

        deallocate(faceYN, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'faceYN')

        deallocate(faceP, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'faceP')



        if (setting%Debug%File%utility_deallocate) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine util_deallocate_elemX_faceX

!
!==========================================================================
!==========================================================================
!

    subroutine util_deallocate_columns()
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   All variables are stored in col_elemX(:) arrays, release the memory here
    !   in all images
    !
    !-----------------------------------------------------------------------------

        character(64) :: subroutine_name = 'util_deallocate_columns'

    !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_deallocate) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        !==== col_elemI====
        deallocate(col_elemI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'col_elemI')
        !==== col_elemP====
        deallocate(npack_elemP, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'npack_elemP')

        deallocate(col_elemP, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'col_elemP')
        !==== col_elemPG====
        deallocate(npack_elemPGalltm, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'npack_elemPGalltm')

        deallocate(col_elemPGalltm, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'col_elemPGalltm')

        deallocate(npack_elemPGac, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'npack_elemPGac')

        deallocate(col_elemPGac, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'col_elemPGac')

        deallocate(npack_elemPGetm, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'npack_elemPGetm')

        deallocate(col_elemPGetm, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'col_elemPGetm')
        !==== col_elemR ====
        deallocate(col_elemR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'col_elemR')
        !==== col_elemSI ====
        deallocate(col_elemSI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'col_elemSI')
        !==== col_elemSR ====
        deallocate(col_elemSR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'col_elemSR')
        !==== col_elemSGR ====
        deallocate(col_elemSGR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'col_elemSGR')
        !==== col_elemWDI ====
        deallocate(col_elemWDI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'col_elemWDI')
        !==== col_elemWDR ====
        deallocate(col_elemWDR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'col_elemWDR')
        !==== col_elemYN ====
        deallocate(col_elemYN, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'col_elemYN')
        !==== col_faceI ====
        deallocate(col_faceI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'col_faceI')
        !==== col_faceM ====
        deallocate(col_faceM, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'col_faceM')
        !==== col_faceP ====
        deallocate(npack_faceP, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'npack_faceP')

        deallocate(col_faceP, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'col_faceP')
        !==== col_faceR ====
        deallocate(col_faceR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'col_faceR')
        !==== col_faceYN ====
        deallocate(col_faceYN, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'col_faceYN')


        if (setting%Debug%File%utility_deallocate) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine util_deallocate_columns


!
!==========================================================================
!==========================================================================
!

    subroutine util_deallocate_bc()

        character(64) :: subroutine_name = 'util_deallocate_bc'

        if (icrash) return
        if (setting%Debug%File%utility_deallocate) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        if (N_flowBC > 0) then
            deallocate(BC%flowI, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check (deallocation_status, emsg, 'BC%flowI')

            deallocate(BC%flowR_timeseries, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check (deallocation_status, emsg, 'BC%flowR_timeseries')

            deallocate(BC%flowIdx, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check (deallocation_status, emsg, 'BC%flowIdx')

            deallocate(BC%flowRI, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check (deallocation_status, emsg, 'BC%flowRI')

            if (allocated(BC%P%BCup)) then
                deallocate(BC%P%BCup, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check (deallocation_status, emsg, 'BC%P%BCup')
            end if

            if (allocated(BC%P%BClat)) then
                deallocate(BC%P%BClat, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check (deallocation_status, emsg, 'BC%P%BClat')
            end if
        end if

        if (N_headBC > 0) then
            deallocate(BC%headI, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check (deallocation_status, emsg, 'BC%headI')

            deallocate(BC%headR_timeseries, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check (deallocation_status, emsg, 'BC%headR_timeseries')

            deallocate(BC%headIdx, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check (deallocation_status, emsg, 'BC%headIdx')

            deallocate(BC%headRI, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check (deallocation_status, emsg, 'BC%headRI')

            if (allocated(BC%P%BCdn)) then
                deallocate(BC%P%BCdn, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check (deallocation_status, emsg, 'BC%P%BCdn')
            end if
        end if

        if (setting%Debug%File%utility_deallocate) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine util_deallocate_bc


!
!==========================================================================
!==========================================================================
!
    subroutine util_deallocate_check(deallocation_status, emsg, locationstring)
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   Checks deallocation status and stops if there is an error
        !
        !-----------------------------------------------------------------------------
            integer,            intent(in   ) :: deallocation_status
            character(len=*),   intent(in   ) :: emsg
            character(len=*),   intent(in   ) :: locationstring !% unique identifier of location

            character(64):: subroutine_name = 'util_deallocate_check'

        !-----------------------------------------------------------------------------
            if (icrash) return
            if (setting%Debug%File%utility) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

            if (deallocation_status > 0) then
                print *, trim(emsg)
                print *, 'variable trying to deallocate = ',trim(locationstring)
                stop
            end if

            if (setting%Debug%File%utility) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

    end subroutine util_deallocate_check

end module utility_deallocate
