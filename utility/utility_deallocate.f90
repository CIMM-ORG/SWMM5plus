module utility_deallocate
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Provides deallocation of run-time arrays
    !%
    !%==========================================================================

    use define_indexes
    use define_globals
    use define_settings, only: setting
    use utility_crash, only: util_crashpoint

    implicit none

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

        ! print *, 'call util_deallocate_linknode() ',this_image()
        call util_deallocate_linknode()

        ! print *, 'call util_deallocate_elemX_faceX() ', this_image()
        call util_deallocate_elemX_faceX()

        ! print *, 'call util_deallocate_elem_boundary_ghost() ', this_image()
        call util_deallocate_elem_boundary_ghost()

        ! print *, 'call util_deallocate_columns() ', this_image()
        call util_deallocate_columns()

        ! print *, 'call util_deallocate_bc() ', this_image()
        call util_deallocate_bc()

        ! print *, 'call util_deallocate_entrapped_air_arrays() ', this_image()
        call util_deallocate_entrapped_air_arrays()

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
            !if (crashYN) return
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
        !% Closing
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
        !% Declarations
            integer,            intent(in   ) :: deallocation_status
            character(len=*),   intent(in   ) :: emsg
            character(len=*),   intent(in   ) :: locationstring !% unique identifier of location
            character(64):: subroutine_name = 'util_deallocate_check'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%utility) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
            if (deallocation_status > 0) then
                print *, trim(emsg)
                print *, 'variable trying to deallocate = ',trim(locationstring)
                call util_crashpoint(79873)
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

        if(this_image() .eq. 1 .and. allocated(output_profile_ids)) then 
            deallocate(output_profile_ids, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_profile_ids')
        end if

        ! deallocate(link%Names, stat=deallocation_status, errmsg=emsg)
        ! call util_deallocate_check(deallocation_status, emsg, 'link%Names')

        !%-------------------------------------------------------------------
        !% Closing
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
            !if (crashYN) return
            if (setting%Debug%File%utility_deallocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------

        !==== elem deallocation ====
        !print *, 'elemR'
        !% HACK: the deallocation of elemR is causing a segmentation fault
        !% for the arch geometry type. So, it has been commented out on 
        !% 11182022. Revisit and fix later
        ! deallocate(elemR, stat=deallocation_status, errmsg=emsg)
        ! call util_deallocate_check(deallocation_status, emsg, 'elemR')

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

        if(allocated(output_types_elemR)) then
            deallocate(output_types_elemR,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_types_elemR')

            deallocate(output_typeProcessing_elemR,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_typeProcessing_elemR')

            deallocate(output_typeMultiplyByBarrels_elemR,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_typeMultiplyByBarrels_elemR')

            deallocate(output_typeNames_elemR,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_typeNames_elemR')

            deallocate(output_typeUnits_elemR,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_typeUnits_elemR')

            deallocate(output_typeNames_withTime_elemR,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_typeNames_withTime_elemR')

            deallocate(output_typeUnits_withTime_elemR,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_typeUnits_withTime_elemR')
        end if

        if(allocated(output_static_types_elemR)) then

            deallocate(output_static_types_elemR,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_static_types_elemR')

            deallocate(output_static_typeProcessing_elemR,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_static_typeProcessing_elemR')

            deallocate(output_static_typeMultiplyByBarrels_elemR,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_static_typeMultiplyByBarrels_elemR')

            deallocate(output_static_typeNames_elemR,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_static_typeNames_elemR')

            deallocate(output_static_typeUnits_elemR,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_static_typeUnits_elemR')

            deallocate(output_static_elem,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_static_elem')

        end if 
        

        if(allocated(output_static_types_Link)) then

            deallocate(output_static_types_Link,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_static_types_Link')


            deallocate(output_static_typeProcessing_Link,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_static_typeProcessing_Link')


            deallocate(output_static_typeMultiplyByBarrels_Link,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_static_typeMultiplyByBarrels_Link')

            deallocate(output_static_typeNames_Link,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_static_typeNames_Link')

            deallocate(output_static_typeUnits_Link,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_static_typeUnits_Link')

            deallocate(output_static_Link,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_static_Link')

        end if 
        
        if(allocated(output_static_types_Node)) then
            
            deallocate(output_static_types_Node,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_static_types_Node')

            deallocate(output_static_typeProcessing_Node,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_static_typeProcessing_Node')

            deallocate(output_static_typeMultiplyByBarrels_Node,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_static_typeMultiplyByBarrels_Node')

            deallocate(output_static_typeNames_Node,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_static_typeNames_Node')

            deallocate(output_static_typeUnits_Node,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_static_typeUnits_Node')

            deallocate(output_static_Node,stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'output_static_Node')

        end if

        !    deallocate(output_types_elemR,stat=deallocation_status, errmsg=emsg)
        !    call util_deallocate_check(deallocation_status, emsg, 'output_types_elemR')

        !%-------------------------------------------------------------------
        !% closing
            if (setting%Debug%File%utility_deallocate) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_deallocate_elemX_faceX
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_deallocate_entrapped_air_arrays ()
        !%-------------------------------------------------------------------
        ! Description:
        !   deallocate the entrapped air arrays
        !%-------------------------------------------------------------------
        !% Declarations:
            character(64) :: subroutine_name = 'util_deallocate_entrapped_air_arrays'
        !%-------------------------------------------------------------------
        !% Preliminaries   
            if (setting%Debug%File%utility_deallocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------

        if (setting%AirTracking%UseAirTrackingYN) then

            deallocate(LinkElemMapsI, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'LinkElemMapsI')

            deallocate(linkAirR, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'linkAirR')

            deallocate(elemAirR, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'elemAirR')

            deallocate(linkAirYN, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'linkAirYN')

            deallocate(elemAirYN, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check(deallocation_status, emsg, 'elemAirYN')

        end if

        !%-------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%utility_deallocate) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_deallocate_entrapped_air_arrays
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_deallocate_elem_boundary_ghost ()
        !%-------------------------------------------------------------------
        ! Description:
        !   deallocate the elemB%R and elemGR that we assigned across all employed images
        !%-------------------------------------------------------------------
        !% Declarations:
            character(64) :: subroutine_name = 'util_deallocate_elem_boundary_ghost'
        !%-------------------------------------------------------------------
        !% Preliminaries   
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
        !% Closing
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
 
        !==== col_elemYN ====
        deallocate(col_elemYN, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'col_elemYN')
        !==== col_faceI ====
        deallocate(col_faceI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg, 'col_faceI')

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
        !% Closing
            if (setting%Debug%File%utility_deallocate) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_deallocate_columns
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_deallocate_bc()
        !%-------------------------------------------------------------------
        !% Description
        !% Deallcoate the boundary condition arrays in bc structure
        !%-------------------------------------------------------------------
        !% Declarations
            character(64) :: subroutine_name = 'util_deallocate_bc'
        !%-------------------------------------------------------------------
        !% Preliminaries:
            if (setting%Debug%File%utility_deallocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------

        if (N_flowBCnode > 0) then
            deallocate(BC%flowI, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check (deallocation_status, emsg, 'BC%flowI')

            deallocate(BC%flowTimeseries, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check (deallocation_status, emsg, 'BC%flowTimeseries')

            deallocate(BC%flowR, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check (deallocation_status, emsg, 'BC%flow_R')

            deallocate(BC%flowYN, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check (deallocation_status, emsg, 'BC%flowYN')

            if (allocated(BC%P%BCup)) then
                deallocate(BC%P%BCup, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check (deallocation_status, emsg, 'BC%P%BCup')
            end if

            if (allocated(BC%P%BClat)) then
                deallocate(BC%P%BClat, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check (deallocation_status, emsg, 'BC%P%BClat')
            end if
        end if

        if (N_headBCnode > 0) then
            deallocate(BC%headI, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check (deallocation_status, emsg, 'BC%headI')

            deallocate(BC%headTimeseries, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check (deallocation_status, emsg, 'BC%headTimeseries')

            deallocate(BC%headYN, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check (deallocation_status, emsg, 'BC%headYN')

            deallocate(BC%headR, stat=deallocation_status, errmsg=emsg)
            call util_deallocate_check (deallocation_status, emsg, 'BC%headR')

            if (allocated(BC%P%BCdn)) then
                deallocate(BC%P%BCdn, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check (deallocation_status, emsg, 'BC%P%BCdn')
            end if
        end if

        !%-------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%utility_deallocate) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_deallocate_bc
!%
!%==========================================================================
!% END MODULE
!%==========================================================================
!%
end module utility_deallocate
