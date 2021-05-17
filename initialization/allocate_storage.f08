module allocate_storage

    use array_index
    use globals
    use utility, only: utility_check_allocation
    use setting_definition, only: setting

    implicit none

!-----------------------------------------------------------------------------
!
! Description:
!   This has all the large array allocation. The only other allocation
!   should be in boundary conditions (module bc)
!   All the top-level storage should be allocated in this module
!
!-----------------------------------------------------------------------------

    private

    ! allocate_storage constants
    integer           :: allocation_status
    character(len=99) ::              emsg

    ! public members
    public :: allocate_linknode_storage

contains

    subroutine allocate_linknode_storage ()
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   Allocates the link and node storage used for the coarse representation
    !   of the network connectivity
    !
    ! Method:
    !   The tables nodeI, linkI, nodeR, linkR, nodeYN, linkYN, are allocated
    !   These are defined in globals.f08). Every time memory is allocated, the
    !   utility_check_allocation functionality (from utility.f08) is used to
    !   determine wheter or not there was an error during the allocation.
    !
    !-----------------------------------------------------------------------------

        character(64) :: subroutine_name = 'allocate_linknode_storage'

    !-----------------------------------------------------------------------------

        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        allocate(nodeI(N_node, ni_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        nodeI(:,:) = nullvalueI

        allocate(linkI(N_link, li_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        linkI(:,:) = nullvalueI

        allocate(B_nodeI(N_node, 3), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        B_nodeI(:,:) = nullvalueI

        allocate(B_linkI(N_link, 2), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        B_linkI(:,:) = nullvalueI

        allocate(nodeR(N_node, nr_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        nodeR(:,:) = nullvalueR

        allocate(linkR(N_link, lr_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        linkR(:,:) = nullvalueR

        allocate(nodeYN(N_node, nYN_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        nodeYN(:,:) = nullvalueL

        allocate(linkYN(N_link, lYN_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        linkYN(:,:) = nullvalueL

        allocate(nodeName(N_node), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)

        allocate(linkName(N_link), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name
    end subroutine allocate_linknode_storage


    subroutine coarray_storage_allocation()
        ! the max_caf_elem and max_caf_face are the maximum length of the coarray 
        ! across all employed images
        ! ==========================
        ! This will be excuted at parallel level
        ! ==========================

        character(64) :: subroutine_name = 'coarray_storage_allocation'
    
        !-----------------------------------------------------------------------------

        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !==== elem allocation ====
        !ncol => Ncol_elemR ! the maxmiumu number of columns
        allocate(elemR_caf(max_caf_elem_N, Ncol_elemR)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        elemR_caf(:,:) = nullvalueR

        !ncol => Ncol_elemI
        allocate(elemI_caf(max_caf_elem_N, Ncol_elemI)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        elemI_caf(:,:) = nullvalueI

        !ncol => Ncol_elemYN
        allocate(elemYN_caf(max_caf_elem_N, Ncol_elemYN)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        elemYN_caf(:,:) = nullvalueL

        !ncol => Ncol_elemP
        allocate(elemP_caf(max_caf_elem_N, Ncol_elemP)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        elemP_caf(:,:) = nullvalueI

        ! Ncol_elemPG not defined yet. Comment out
        !ncol => Ncol_elemPG
        !allocate(elemPG_caf(max_caf_elem_N, Ncol_elemPG)[*], stat=allocation_status, errmsg=emsg)
        !call utility_check_allocation(allocation_status, emsg)
        !elemPG_caf(:,:)[:] = nullvalueI

        !ncol => Ncol_elemSI
        allocate(elemSI_caf(max_caf_elem_N, Ncol_elemSI)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        elemSI_caf(:,:) = nullvalueI

        !ncol => Ncol_elemSR
        allocate(elemSR_caf(max_caf_elem_N, Ncol_elemSR)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        elemSR_caf(:,:) = nullvalueR

        !==== face allocation ====
        !ncol => Ncol_faceR 
        allocate(faceR_caf(max_caf_face_N, Ncol_faceR )[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        faceR_caf(:,:) = nullvalueR

        !ncol=> Ncol_faceI
        allocate(faceI_caf(max_caf_face_N, Ncol_faceI)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        faceI_caf(:,:) = nullvalueI

        !ncol=> Ncol_faceYN
        allocate(faceYN_caf(max_caf_face_N, Ncol_faceYN)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        faceYN_caf(:,:) = nullvalueL


        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name


    end subroutine coarray_storage_allocation


    subroutine allocate_columns()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   All the enumerated variables can not be used as pointers. Thus the 
        !   variables are stored in col_elemX(:) arrays that is a target
        !
        !-----------------------------------------------------------------------------

        character(64) :: subroutine_name = 'allocate_columns'

        !-----------------------------------------------------------------------------
        
        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !% allocation of the col_elemX and npack_elemX
        call allocate_col_elemI
        call allocate_col_elemP
        call allocate_col_elemPGalltm
        call allocate_col_elemPGac
        call allocate_col_elemPGetm
        call allocate_col_elemR
        call allocate_col_elemSI
        call allocate_col_elemSR
        call allocate_col_elemSGR
        call allocate_col_elemWDI
        call allocate_col_elemWDR
        call allocate_col_elemYN
        call allocate_col_faceI
        call allocate_col_faceM
        call allocate_col_faceP
        call allocate_col_faceR
        call allocate_col_faceYN

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name

    end subroutine allocate_columns

    subroutine allocate_col_elemI()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemI is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated ei_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'allocate_col_elemI'

        !-----------------------------------------------------------------------------
        
        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        ncol => Ncol_elemI

        !% allocate an array for storing the column 
        allocate( col_elemI(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemI(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name

    end subroutine allocate_col_elemI


    subroutine allocate_col_elemP()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemP is a vector of the columns in the elemP arrays
        !   that correspond to the enumerated ep_... array_index parameter
        !
        !   the npack_elemP(:) vector contains the number of packed elements
        !   for a given column.
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'allocate_col_elemP'

        !-----------------------------------------------------------------------------
        
        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemP

        !% allocate an array for storing the size of each packed type
        allocate( npack_elemP(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% allocate an array for storing the column of each packed type
        allocate( col_elemP(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemP(:) = [(ii,ii=1,ncol)]

        !% zero the number of packed items (to be defined in the packing)
        npack_elemP(:) = 0

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name

    end subroutine allocate_col_elemP


    subroutine allocate_col_elemPGalltm()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemPGalltm is a vector of the columns in the elemPGalltm array
        !   that correspond to the enumerated ePG_... array_index parameters
        !
        !   the npack_elemPGalltm(:) vector contains the number of packed elements
        !   for a given column.
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'allocate_col_elemPGalltm'

        !-----------------------------------------------------------------------------
        
        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemPG !% whatever the last item in the enumerator

        !% allocate an array for storing the size of each packed type
        allocate( npack_elemPGalltm(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% allocate an array for storing the enum type of each column
        allocate( col_elemPGalltm(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemPGalltm(:) = [(ii,ii=1,ncol)]

        !% zero the number of packed items (to be defined in the packing)
        npack_elemPGalltm(:) = 0

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name

    end subroutine allocate_col_elemPGalltm


    subroutine allocate_col_elemPGac()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemPGac is a vector of the columns in the elemPGac arrays
        !   that correspond to the enumerated epg_... array_index parameters
        !
        !   the npack_elemPGac(:) vector contains the number of packed elements
        !   for a given column.
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'allocate_col_elemPGac'

        !-----------------------------------------------------------------------------
        
        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemPG !% whatever the last item in the enumerator

        !% allocate an array for storing the size of each packed type
        allocate( npack_elemPGac(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% allocate an array for storing the enum type of each column
        allocate( col_elemPGac(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemPGac(:) = [(ii,ii=1,ncol)]

        !% zero the number of packed items (to be defined in the packing)
        npack_elemPGac(:) = 0

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name

    end subroutine allocate_col_elemPGac


    subroutine allocate_col_elemPGetm()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemPGetm is a vector of the columns in the elemPGetm arrays
        !   that correspond to the enumerated epg_... array_index parameters
        !    
        !   the npack_elemPGetm(:) vector contains the number of packed elements
        !   for a given column.
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'allocate_col_elemPGetm'

        !-----------------------------------------------------------------------------
        
        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemPG !% whatever the last item in the enumerator

        !% allocate an array for storing the size of each packed type
        allocate( npack_elemPGetm(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% allocate an array for storing the enum type of each column
        allocate( col_elemPGetm(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemPGetm(:) = [(ii,ii=1,ncol)]

        !% zero the number of packed items (to be defined in the packing)
        npack_elemPGetm(:) = 0

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name

    end subroutine allocate_col_elemPGetm


    subroutine allocate_col_elemR()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemR is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated er_... array_index parameter
        !
        !-----------------------------------------------------------------------------
        
        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'allocate_col_elemR'

        !-----------------------------------------------------------------------------
        
        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemR

        !% allocate an array for storing the column 
        allocate( col_elemR(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemR(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name

    end subroutine allocate_col_elemR


    subroutine allocate_col_elemSI()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemSI is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated eSI_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'allocate_col_elemSI'

        !-----------------------------------------------------------------------------
        
        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemSI

        !% allocate an array for storing the column 
        allocate( col_elemR(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemSI(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name

    end subroutine allocate_col_elemSI


    subroutine allocate_col_elemSR()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemSR is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated eSR_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'allocate_col_elemSR'

        !-----------------------------------------------------------------------------
        
        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemSR

        !% allocate an array for storing the column 
        allocate( col_elemSR(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemSR(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name

    end subroutine allocate_col_elemSR


    subroutine allocate_col_elemSGR()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemSGR is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated eSGR_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'allocate_col_elemSGR'

        !-----------------------------------------------------------------------------
        
        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemSGR

        !% allocate an array for storing the column 
        allocate( col_elemSGR(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemSGR(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name

    end subroutine allocate_col_elemSGR


    subroutine allocate_col_elemWDI()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemWDI is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated eWDI_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'allocate_col_elemWDI'

        !-----------------------------------------------------------------------------
        
        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemWDI

        !% allocate an array for storing the column 
        allocate( col_elemWDI(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemWDI(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name

    end subroutine allocate_col_elemWDI


    subroutine allocate_col_elemWDR()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemWDR is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated eWDR_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'allocate_col_elemWDI'

        !-----------------------------------------------------------------------------
        
        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemWDR

        !% allocate an array for storing the column 
        allocate( col_elemWDR(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemWDR(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name

    end subroutine allocate_col_elemWDR


    subroutine allocate_col_elemYN()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemYN is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated eYN_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'allocate_col_elemYN'

        !-----------------------------------------------------------------------------
        
        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemYN

        !% allocate an array for storing the column 
        allocate( col_elemYN(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemYN(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name

    end subroutine allocate_col_elemYN


    subroutine allocate_col_faceI()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_faceI is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated fi_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'allocate_col_faceI'

        !-----------------------------------------------------------------------------
        
        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_faceI

        !% allocate an array for storing the column 
        allocate( col_faceI(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_faceI(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name

    end subroutine allocate_col_faceI


    subroutine allocate_col_faceM()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_faceM is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated fM_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'allocate_col_faceM'

        !-----------------------------------------------------------------------------
        
        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_faceM

        !% allocate an array for storing the column 
        allocate( col_faceM(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_faceM(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name

    end subroutine allocate_col_faceM


    subroutine allocate_col_faceP()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   packed arrays for faces
        !   the col_faceP is a vector of the columns in the faceP arrays
        !   that correspond to the enumerated fp_... array_index parameters
        !
        !   the npack_faceP(:) vector contains the number of packed elements
        !   for a given column.
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'allocate_col_faceP'

        !-----------------------------------------------------------------------------
        
        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_faceP 

        !% allocate an array for storing the size of each packed type
        allocate( npack_faceP(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% allocate an array for storing the column of each packed type
        allocate( col_faceP(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_faceP(:) = [(ii,ii=1,ncol)]

        !% zero the number of packed items (to be defined in the packing)
        npack_faceP(:) = 0

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name

    end subroutine allocate_col_faceP


    subroutine allocate_col_faceR()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_faceR is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated fr_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'allocate_col_faceR'

        !-----------------------------------------------------------------------------
        
        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_faceR  !global

        !% allocate an array of column indexes that can be used as targets of pointers
        allocate( col_faceR(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_faceR(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name

    end subroutine allocate_col_faceR


    subroutine allocate_col_faceYN()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_faceYN is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated fYN_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'allocate_col_faceYN'

        !-----------------------------------------------------------------------------
        
        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_faceYN  !global

        !% allocate an array of column indexes that can be used as targets of pointers
        allocate( col_faceYN(ncol), stat=allocation_status, errmsg= emsg)
        call utility_check_allocation (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_faceYN(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name

    end subroutine allocate_col_faceYN

    ! subroutine allocate_bc()
    ! !
    ! ! allocate storage for boundary conditions.
    ! !
    ! !-----------------------------------------------------------------------------

    !     character(64)      :: subroutine_name = 'allocate_bc'
    !     integer            :: ii
    !     integer            :: allocation_status
    !     character(len=99)  :: emsg

    ! !-----------------------------------------------------------------------------

    !     if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

    !     !% the Upstream and Downstream bc structure
    !     allocate(bcdataUp(N_BCupstream), stat=allocation_status, errmsg=emsg)
    !     call utility_check_allocation (allocation_status, emsg)

    !     allocate(bcdataDn(N_BCdnstream), stat=allocation_status, errmsg=emsg)
    !     call utility_check_allocation (allocation_status, emsg)

    !     !% the downstream arrays - HACK default downstream is elevation
    !     do ii=1,N_BCdnstream
    !         bcdataDn(ii)%idx = ii
    !         bcdataDn(ii)%Updn          = bc_updn_downstream
    !         bcdataDn(ii)%Category      = bc_category_elevation
    !         bcdataDn(ii)%NodeID         = nullvalueI
    !         bcdataDn(ii)%FaceID         = nullvalueI
    !         bcdataDn(ii)%ElemGhostID    = nullvalueI
    !         bcdataDn(ii)%ElemInsideID   = nullvalueI
    !         bcdataDn(ii)%ThisValue      = nullvalueR
    !         bcdataDn(ii)%ThisTime       = nullvalueR
    !         bcdataDn(ii)%ThisFlowrate   = nullvalueR
    !         allocate(bcdataUp(ii)%TimeArray(setting%Constant%BCSlots))
    !         allocate(bcdataUp(ii)%ValueArray(setting%Constant%BCSlots))
    !     end do

    !     !% the upstream arrays = HACK default upstream is flowrate
    !     do ii=1,N_BCupstream
    !         bcdataUp(ii)%idx = ii
    !         bcdataUp(ii)%Updn       = bc_updn_upstream
    !         bcdataUp(ii)%category   = bc_category_inflowrate
    !         bcdataUp(ii)%NodeID         = nullvalueI
    !         bcdataUp(ii)%FaceID         = nullvalueI
    !         bcdataUp(ii)%ElemGhostID    = nullvalueI
    !         bcdataUp(ii)%ElemInsideID   = nullvalueI
    !         bcdataUp(ii)%ThisValue      = nullvalueR
    !         bcdataUp(ii)%ThisTime       = nullvalueR
    !         bcdataUp(ii)%ThisFlowrate   = nullvalueR
    !         allocate(bcdataUp(ii)%TimeArray(setting%Constant%BCSlots))
    !         allocate(bcdataUp(ii)%ValueArray(setting%Constant%BCSlots))
    !     end do

    !     if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name
    ! end subroutine allocate_bc

end module allocate_storage
