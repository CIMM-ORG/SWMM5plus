module utility_deallocate

    use define_indexes
    use define_globals
    use define_settings, only: setting

    implicit none

!%-------------------------------------------------------------------------
!% Description:
!%   This has all the array and variables deallocation.
!%   All the top-level storage should be deallocated in this module
!%-------------------------------------------------------------------------

    private

    integer           :: deallocation_status
    character(len=99) ::              emsg

    public :: util_deallocate_network_data
    public :: util_deallocate_partitioning_arrays
    public :: util_deallocate_check

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine util_deallocate_network_data()
        !%------------------------------------------------------------------
        !% Description:
        !%   deallocate all the data relate to defined network, including:
        !%   linkX, nodeX, elemX, faceX
        !%-------------------------------------------------------------------
        !% Declarations:
            character(64) :: subroutine_name = 'util_deallocate_network_data'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%utility_deallocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------

        !print *, 'call util_deallocate_linknode() ',this_image()
        call util_deallocate_linknode()

        !print *, 'call util_deallocate_elemX_faceX() ', this_image()
        call util_deallocate_elemX_faceX()

        !print *, 'call util_deallocate_elem_boundary_ghost() ', this_image()
        call util_deallocate_elem_boundary_ghost()

        !print *, 'call util_deallocate_columns() ', this_image()
        call util_deallocate_columns()

        !print *, 'call util_deallocate_bc() ', this_image()
        call util_deallocate_bc()

        !%-------------------------------------------------------------------
        !% closing
            if (setting%Debug%File%utility_deallocate) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_deallocate_network_data
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_deallocate_partitioning_arrays()
        !%-------------------------------------------------------------------
        !% Description
        !%   deallocates the partitioning arrays storage used for storing partitioned
        !%   network information
        !%-------------------------------------------------------------------
        !% Declarations
            character(64) :: subroutine_name = 'util_deallocate_partitioning_arrays'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (crashYN) return
            if (setting%Debug%File%utility_deallocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------        

        deallocate(adjacent_links, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg,'adjacent_links')

        deallocate(elem_per_image, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elem_per_image')

        deallocate(image_full, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'image_full')

        !%-------------------------------------------------------------------
        !% closing
            if (setting%Debug%File%utility_deallocate) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_deallocate_partitioning_arrays
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_deallocate_check(deallocation_status, emsg, locationstring)
        !%-------------------------------------------------------------------
        !% Description:
        !%   Checks deallocation status and stops if there is an error
        !%-------------------------------------------------------------------
        !% Declaratoins
            integer,            intent(in   ) :: deallocation_status
            character(len=*),   intent(in   ) :: emsg
            character(len=*),   intent(in   ) :: locationstring !% unique identifier of location
            character(64):: subroutine_name = 'util_deallocate_check'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (crashYN) return
            if (setting%Debug%File%utility) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------

            if (deallocation_status > 0) then
                print *, trim(emsg)
                print *, 'variable trying to deallocate = ',trim(locationstring)
                stop
            end if

        !%-------------------------------------------------------------------    
            if (setting%Debug%File%utility) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_deallocate_check    
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine util_deallocate_linknode()
        !%-------------------------------------------------------------------
        !% Description:
        !%   deallocates the link and node storage used for the coarse representation
        !%   of the network connectivity
        !%
        !% Method:
        !%   Release the memory used by the tables node%I, link%I, node%R, link%R, node%YN,
        !%   link%YN.
        !%   These are defined in globals.f08), and allocated in allocate_storage.f08
        !%   Every time memory is deallocated, the utility_check_deallocation functionality
        !%   (from utility.f08) is used to determine wheter or not there was an error during
        !%   the deallocation.
        !%-------------------------------------------------------------------
        !% Declarations
            integer       :: ii
            character(64) :: subroutine_name = 'util_deallocate_linknode'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (crashYN) return
            if (setting%Debug%File%utility_deallocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
                
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
        ! do ii = 1, size(link%Names)
        !     deallocate(link%Names(ii)%str, stat=deallocation_status, errmsg=emsg)
        !     call util_deallocate_check(deallocation_status, emsg, 'link%Names(ii)%str')
        ! end do
        do ii = 1, size(node%Names)
            deallocate(node%Names(ii)%str, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'node%Names(ii)%str')
        end do

       
        deallocate(node%Names, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'node%Names')

        ! deallocate(link%Names, stat=deallocation_status, errmsg=emsg)
        ! call util_deallocate_check(deallocation_status, emsg, 'link%Names')

    
        !% brh20211217 -- 
        !% the followin commented asobsolete approach to storing output nodes and links
        !deallocate(node_output_idx,stat=deallocation_status,errmsg=emsg)
        !call util_deallocate_check(deallocation_status, emsg, 'node_output_idx')

        !deallocate(link_output_idx,stat=deallocation_status,errmsg=emsg)
        !call util_deallocate_check(deallocation_status, emsg, 'link_output_idx')

        ! if (allocated(node%P%have_output)) then
        !     deallocate(node%P%have_output,stat=deallocation_status,errmsg=emsg)
        !     call util_deallocate_check(deallocation_status, emsg, 'node%P%have_output')
        ! end if

        ! if (allocated(link%P%have_output)) then
        !     deallocate(link%P%have_output,stat=deallocation_status,errmsg=emsg)
        !     call util_deallocate_check(deallocation_status, emsg, 'link%P%have_output')
        ! end if

        !%-------------------------------------------------------------------
        !% closing
            if (setting%Debug%File%utility_deallocate) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_deallocate_linknode
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_deallocate_elemX_faceX ()
        !%-------------------------------------------------------------------
        ! Description:
        !   Simply deallocate the elemX and faceX thatwe assigned across all employed images
        !%-------------------------------------------------------------------
        !% Declarations:
            character(64) :: subroutine_name = 'util_deallocate_elemX_faceX'
        !%-------------------------------------------------------------------
        !% Preliminaries   
            if (crashYN) return
            if (setting%Debug%File%utility_deallocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------

        !==== elem deallocation ====
        !print *, 'elemR'
        deallocate(elemR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elemR')

        !print *, 'elemI'
        deallocate(elemI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elemI')

        !print *, 'elemYN'
        deallocate(elemYN, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elemYN')

        !print *, 'elemP'
        deallocate(elemP, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elemP')

        !print *, 'elemPGetm'
        deallocate(elemPGetm, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elemPGetm')

        !print *, 'elemPGac'
        deallocate(elemPGac, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elemPGac')

        !print *, 'elemPGalltm'
        deallocate(elemPGalltm, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elemPGalltm')

        !print *, 'elemSI'
        deallocate(elemSI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elemSI')

        !print *, 'elemSR'
        deallocate(elemSR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'elemSR')

        !==== face deallocation ====
        !print *, 'faceR'
        deallocate(faceR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'faceR')

        !print *, 'faceI'
        deallocate(faceI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'faceI')

        !print *, 'faceYN'
        deallocate(faceYN, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'faceYN')

        !print *, 'faceP'
        deallocate(faceP, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'faceP')

        !%-------------------------------------------------------------------
        !% closing
            if (setting%Debug%File%utility_deallocate) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_deallocate_elemX_faceX
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_deallocate_elem_boundary_ghost ()
        !%-------------------------------------------------------------------
        ! Description:
        !   Simply deallocate the elemB%R and elemGR that we assigned across all employed images
        !%-------------------------------------------------------------------
        !% Declarations:
            character(64) :: subroutine_name = 'util_deallocate_elem_boundary_ghost'
        !%-------------------------------------------------------------------
        !% Preliminaries   
            if (crashYN) return
            if (setting%Debug%File%utility_deallocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------

        !==== elemB deallocation ====
        if (num_images() > 1) then
            deallocate(elemB%R, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check (deallocation_status, emsg, 'elemB%R')

            deallocate(elemB%I, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check (deallocation_status, emsg, 'elemB%I')

            deallocate(elemGR, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'elemGR')
        end if
        !%-------------------------------------------------------------------
        !% closing
            if (setting%Debug%File%utility_deallocate) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_deallocate_elem_boundary_ghost
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_deallocate_columns()
        !%-------------------------------------------------------------------
        !% Description:
        !%   All variables are stored in col_elemX(:) arrays, release the memory here
        !%   in all images
        !%-------------------------------------------------------------------
            character(64) :: subroutine_name = 'util_deallocate_columns'
        !%-------------------------------------------------------------------
            if (crashYN) return
            if (setting%Debug%File%utility_deallocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------

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
       ! deallocate(col_faceM, stat=deallocation_status, errmsg=emsg)
       ! call util_deallocate_check(deallocation_status, emsg, 'col_faceM')
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

        !%-------------------------------------------------------------------
            if (setting%Debug%File%utility_deallocate) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_deallocate_columns
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_deallocate_bc()
        !%-------------------------------------------------------------------
            character(64) :: subroutine_name = 'util_deallocate_bc'
        !%-------------------------------------------------------------------
            if (crashYN) return
            if (setting%Debug%File%utility_deallocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------

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
        !%-------------------------------------------------------------------
        if (setting%Debug%File%utility_deallocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_deallocate_bc
!%
!%==========================================================================
!% END MODULE
!%==========================================================================
!%
end module utility_deallocate
