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
        if (setting%Debug%File%utility_deallocate) print *, '*** enter ',subroutine_name

        call util_deallocate_linknode()
        call util_deallocate_elemX_faceX()
        call util_deallocate_columns()
        call util_deallocate_bc()

        if (setting%Debug%File%utility_deallocate) print *, '*** leave ',subroutine_name
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
    !   Release the memory used by the tables nodeI, linkI, nodeR, linkR, nodeYN,
    !   linkYN.
    !   These are defined in globals.f08), and allocated in allocate_storage.f08
    !   Every time memory is deallocated, the utility_check_deallocation functionality 
    !   (from utility.f08) is used to determine wheter or not there was an error during
    !   the deallocation.
    !
    !-----------------------------------------------------------------------------

        character(64) :: subroutine_name = 'util_deallocate_linknode_storage'

    !-----------------------------------------------------------------------------
        if (setting%Debug%File%utility_deallocate) print *, '*** enter ',subroutine_name

        deallocate(nodeI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(linkI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(nodeR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(linkR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(nodeYN, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(linkYN, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(nodeName, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(linkName, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        if (setting%Debug%File%utility_deallocate) print *, '*** leave ',subroutine_name

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
        if (setting%Debug%File%utility_deallocate) print *, '*** enter ',subroutine_name
    
        deallocate(adjacent_links, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
    
        deallocate(elem_per_image, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
    
        deallocate(image_full, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
    
        if (setting%Debug%File%utility_deallocate) print *, '*** leave ',subroutine_name
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
        if (setting%Debug%File%utility_deallocate) print *, '*** enter ',subroutine_name

        !==== elem deallocation ====
        deallocate(elemR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(elemI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(elemYN, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(elemP, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(elemPGetm, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(elemPGac, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
        
        deallocate(elemPGalltm, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(elemSI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(elemSR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        !==== face deallocation ====
        deallocate(faceR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
        
        deallocate(faceI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
        
        deallocate(faceYN, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
        
        deallocate(faceP, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
        
        
        if (setting%Debug%File%utility_deallocate) print *, '*** leave ',subroutine_name
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
        if (setting%Debug%File%utility_deallocate) print *, '*** enter ',subroutine_name
        
        !==== col_elemI====
        deallocate(col_elemI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
        !==== col_elemP====
        deallocate(npack_elemP, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(col_elemP, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
        !==== col_elemPG====
        deallocate(npack_elemPGalltm, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(col_elemPGalltm, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(npack_elemPGac, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(col_elemPGac, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(npack_elemPGetm, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(col_elemPGetm, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
        !==== col_elemR ====
        deallocate(col_elemR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
        !==== col_elemSI ====
        deallocate(col_elemSI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
        !==== col_elemSR ====
        deallocate(col_elemSR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
        !==== col_elemSGR ====
        deallocate(col_elemSGR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
        !==== col_elemWDI ====
        deallocate(col_elemWDI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
        !==== col_elemWDR ====
        deallocate(col_elemWDR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
        !==== col_elemYN ====
        deallocate(col_elemYN, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
        !==== col_faceI ====
        deallocate(col_faceI, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
        !==== col_faceM ====
        deallocate(col_faceM, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
        !==== col_faceP ====
        deallocate(npack_faceP, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)

        deallocate(col_faceP, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
        !==== col_faceR ====
        deallocate(col_faceR, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
        !==== col_faceYN ====
        deallocate(col_faceYN, stat=deallocation_status, errmsg=emsg)
        call util_deallocate_check(deallocation_status, emsg)
 

        if (setting%Debug%File%utility_deallocate) print *, '*** leave ',subroutine_name
    end subroutine util_deallocate_columns


!
!==========================================================================
!==========================================================================
!

    subroutine util_deallocate_bc()
        
        character(64) :: subroutine_name = 'util_deallocate_bc'
        
        if (setting%Debug%File%utility_deallocate) print *, '*** enter ',subroutine_name
        
        
        if (setting%Debug%File%utility_deallocate) print *, '*** leave ',subroutine_name
    end subroutine util_deallocate_bc


!
!==========================================================================
!==========================================================================
!
    
    
    subroutine util_deallocate_check(deallocation_status, emsg)
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   Checks deallocation status and stops if there is an error
        !
        !-----------------------------------------------------------------------------
            integer,            intent(in   ) :: deallocation_status
            character(len=*),   intent(in   ) :: emsg

            character(64):: subroutine_name = 'util_deallocate_check'
    
        !-----------------------------------------------------------------------------
    
            if (setting%Debug%File%utility) print *, '*** enter ',subroutine_name
    
            if (deallocation_status > 0) then
                print *, trim(emsg)
                stop
            end if
    
            if (setting%Debug%File%utility) print *, '*** leave ',subroutine_name
    
    end subroutine util_deallocate_check

    subroutine util_deallocate_interior_boundary_elements()
        deallocate(elemGR)
        deallocate(elemGI)
    end subroutine util_deallocate_interior_boundary_elements
end module utility_deallocate