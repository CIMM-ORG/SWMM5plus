module utility

    use define_indexes
    use define_keys
    use define_globals
    use define_settings, only: setting
    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

!-----------------------------------------------------------------------------
!
! Description:
!   Utility routines that may be called in a number of places
!
!-----------------------------------------------------------------------------

    private

    public :: util_print_programheader
    public :: util_count_node_types
    public :: util_sign_with_ones
    public :: util_sign_with_ones_or_zero
    public :: util_print_warning
    public :: util_linspace
    
    public :: util_accumulate_volume_conservation
    public :: util_total_volume_conservation

    public :: util_find_elements_in_link
    public :: util_find_elements_in_junction_node
    public :: util_find_neighbors_of_CC_element
    public :: util_find_neighbors_of_JM_element

    public :: util_CLprint

    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine util_print_programheader ()
        if (this_image() .ne. 1) return
        write(*,*) " "
        write(*,*) "*********************************************************************"
        write(*,*) "*                            SWMM5+                                 *"
        write(*,*) "*                   beta 0.1 release 2022-01-10                     *"
        write(*,*) "*  A public-domain, finite-volume hydraulics engine for EPA SWMM.   *"
        write(*,*) "*                        developed by CIMM                          *"
        write(*,*) "*                                                                   *"
        write(*,*) "* CIMM is the Center for Infrastructure Modeling and Management     *"
        write(*,*) "* funded under US EPA Cooperative Agreement 83595001 awarded to the *" 
        write(*,*) "* University of Texas at Austin, 2017-23. PI: Prof. Ben R. Hodges   *"
        write(*,*) "*                                                                   *"
        write(*,*) "* Code authors:                                                     *"
        write(*,*) "*    2016-2022 Dr. Ben R. Hodges                                    *"
        write(*,*) "*    2016-2022 Edward Tiernan                                       *"
        write(*,*) "*    2019-2022 Sazzad Sharior                                       *"
        write(*,*) "*    2019-2022 Gerardo Riano-Briceno                                *"
        write(*,*) "*    2019-2022 Eric Jenkins                                         *"
        write(*,*) "*    2020-2022 Dr. Cheng-Wei (Justin) Yu                            *"
        write(*,*) "*    2021-2022 Christopher Brashear                                 *"
        write(*,*) "*    2021-2022 Abdulmuttalib Lokhandwala                            *"
        write(*,*) "*    2018-2020 Dr. Ehsan Madadi-Kandjani                            *"
        write(*,*) "*********************************************************************"
        write(*,*)
        write(*,"(A,i5,A)") "Simulation starts with ",num_images()," processors"
        write(*,*) ' '
        write(*,"(A)") 'Using the following files and folders:'
        write(*,"(A)") '  SWMM Input file   : '//trim(setting%File%inp_file)
        write(*,"(A)") '  SWMM Report file  : '//trim(setting%File%rpt_file)
        write(*,"(A)") '  SWMM Output file  : '//trim(setting%File%out_file)
        write(*,"(A)") '  Output folder     : '//trim(setting%File%output_timestamp_subfolder )
        write(*,"(A)") '  Settings file     : '//trim(setting%File%setting_file)
    end subroutine util_print_programheader  
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_count_node_types(N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2, N_nJ1)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% This subroutine uses the vectorized count() function to search the array for
        !% numberof instances of each node type
        !%-----------------------------------------------------------------------------
        integer, intent(in out) :: N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2, N_nJ1
        integer :: ii
        !%-----------------------------------------------------------------------------
        N_nBCup = count(node%I(:, ni_node_type) == nBCup)
        N_nBCdn = count(node%I(:, ni_node_type) == nBCdn)
        N_nJm = count(node%I(:, ni_node_type) == nJM)
        N_nStorage = count(node%I(:, ni_node_type) == nStorage)
        N_nJ2 = count(node%I(:, ni_node_type) == nJ2)
        N_nj1 = count(node%I(:, ni_node_type) == nJ1)

    end subroutine util_count_node_types
!%
!%==========================================================================
!%==========================================================================
!%
    pure elemental real(8) function util_sign_with_ones &
        (inarray) result (outarray)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% returns is an array of real ones with the sign of the inarray argument
        !%-----------------------------------------------------------------------------
        real(8),      intent(in)    :: inarray
        !%-----------------------------------------------------------------------------
        outarray = oneR
        outarray = sign(outarray,inarray)

    end function util_sign_with_ones

!%
!%==========================================================================
!%==========================================================================
!%
    function util_sign_with_ones_or_zero &
        (inarray) result (outarray)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% returns is an array of real ones with the sign of the inarray argument
        !% for non-zero inarray. For zero inarray returns zero.
        !%-----------------------------------------------------------------------------
        real(8),intent(in) :: inarray(:)
        real(8)            :: outarray(size(inarray,1))
        !%-----------------------------------------------------------------------------
        outarray = oneR
        outarray = sign(outarray,inarray)

        !% brh 20211227 
        where(inarray == zeroR)
            outarray = zeroI
        end where

    end function util_sign_with_ones_or_zero

!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_print_warning(msg,async)
        !%------------------------------------------------------------------
        !% Used for opening up the warning files and writing to the file
        !%------------------------------------------------------------------
            character(len = *), intent(in) :: msg
            logical, optional, intent(in) :: async
            logical :: async_actual
        !%------------------------------------------------------------------
        if (present(async)) then
            async_actual = async
        else
            async_actual = .true.
        end if
        if (this_image() == 1) then
            print *, "Warning: "//trim(msg)
        else if (async_actual) then
            print *, "Warning: "//trim(msg)
        end if

    end subroutine util_print_warning
!%
!%==========================================================================
!%==========================================================================
!%
    function util_linspace(startPoint,endPoint,N) result(outArray)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% similar to python/matlab linspace
        !%-----------------------------------------------------------------------------
        real(8), intent(in)  :: startPoint 
        real(8), intent(in)  :: endPoint
        integer, intent(in)  :: N
        real(8)              :: delta
        real(8), allocatable :: outArray(:)
        integer :: ii
        !%-----------------------------------------------------------------------------
        !% calculate step size
        delta = (endPoint - startPoint)/real(N-1,8)

        !% allocate the outArry based on number of samples
        allocate(outArray(N))

        do ii = 1, N
            outArray(ii) = startPoint + (ii-1)*delta
        end do

    end function util_linspace
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine util_accumulate_volume_conservation ()
        !%------------------------------------------------------------------
        !% Description:
        !% Computes/stores the cumulative mass conservation for each CC and JM
        !% element that is NOT small volume or zero depth
        !% HACK -- JM not completed
        !% HACK -- need a separate volume cons for the small/zero losses
        !%------------------------------------------------------------------
        !% Declarations:
            real(8), pointer :: eCons(:), fQ(:), eQLat(:), VolNew(:), VolOld(:), dt
            integer, pointer :: thisColCC, thisColJM, npack, thisP(:)
            integer, pointer :: fdn(:), fup(:), BranchExists(:)
            integer :: ii,kk
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:
            !thisColCC => col_elemP(ep_CC_Q_NOTsmalldepth)
            thisColCC => col_elemP(ep_CC_ALLtm)
            thisColJM => col_elemP(ep_JM_ALLtm)
            fQ      => faceR(:,fr_Flowrate_Conservative)
            eCons   => elemR(:,er_VolumeConservation)
            eQLat   => elemR(:,er_FlowrateLateral)
            VolNew  => elemR(:,er_Volume)  
            VolOld  => elemR(:,er_Volume_N0)  
            fup     => elemI(:,ei_Mface_uL)
            fdn     => elemI(:,ei_Mface_dL)
            dt      => setting%Time%Hydraulics%Dt
            BranchExists => elemSI(:,esi_JunctionBranch_Exists)
        !%------------------------------------------------------------------
        !% --- for the CC elements
        npack   => npack_elemP(thisColCC)
        if (npack > 0) then
            thisP => elemP(1:npack,thisColCC)

            !where (.not. elemYN(thisP,eYN_isZeroDepth))
            !% --- sum of the net inflow and lateral flow should be the change in volume from head
            !% --- for output, use an accumulator
            eCons(thisP) = eCons(thisP)                                            &
                         + dt * ( fQ(fup(thisP)) - fQ(fdn(thisP)) + eQlat(thisP) ) &
                         - (VolNew(thisP) - VolOld(thisP))

            !endwhere             
            !% --- for debugging, switch to using non-cumulative            
            !eCons(thisP) = dt * (fQ(fup(thisP)) - fQ(fdn(thisP)) + eQlat(thisP)) &
            !              - (VolNew(thisP) - VolOld(thisP))      
                          
            ! do ii = 1,size(thisP)
            !     if (abs(eCons(thisP(ii))) > 1.0e-4) then
            !         print *, 'CONSERVATION ISSUE CC', ii, thisP(ii), this_image()
            !         print *,  eCons(thisP(ii))
            !         print *, fQ(fup(thisP(ii))), fQ(fdn(thisP(ii))), eQlat(thisP(ii))
            !         print *, ' vol ',VolNew(thisP(ii)), VolOld(thisP(ii))
            !         do kk=1,num_images()
            !             stop 358783
            !         end do
            !     end if
            ! end do    

        end if

        !% HACK -- need an equivalent of the ep_CC_Q_NOTsmall depth for JM
        !% to make this work

        !% for the JM elements
        npack   => npack_elemP(thisColJM)
        if (npack > 0) then
            thisP => elemP(1:npack,thisColJM)

            eCons(thisP) = eCons(thisP) + dt * eQlat(thisP) - (VolNew(thisP) - VolOld(thisP))

            !% debug, compute just this time step conservation
            !eCons(thisP) = dt * eQlat(thisP) - (VolNew(thisP) - VolOld(thisP))

            do ii=1,max_branch_per_node,2
                eCons(thisP) = eCons(thisP)                                           &
                             + dt * ( real(BranchExists(thisP+ii  ),8) * fQ(fup(thisP+ii)) &
                                    - real(BranchExists(thisP+ii+1),8) * fQ(fdn(thisP+ii+1)) )
            end do
            
            ! do ii = 1,size(thisP)
            !     if (abs(eCons(thisP(ii))) > 1.0e-2) then
            !         print *, ' '
            !         print *, 'CONSERVATION ISSUE JM', ii, thisP(ii), this_image()
            !         print *,  'net cons ',eCons(thisP(ii))
            !         do kk = 1,max_branch_per_node,2
            !             print *, 'branch ',kk, fQ(fup(thisP(ii)+kk)  ) * real(BranchExists(thisP(ii)+kk  ),8) ,  &
            !                      fQ(fdn(thisP(ii)+kk+1)) * real(BranchExists(thisP(ii)+kk+1),8) 
            !         end do
            !         print *, ' vol ', VolNew(thisP(ii)), VolOld(thisP(ii))
            !         print *, 'd vol' , (VolNew(thisP(ii)) - VolOld(thisP(ii)))
            !         print *, 'net Q' , dt * (eQlat(thisP(ii)) + fQ(fup(thisP(ii)+1)  )+ fQ(fup(thisP(ii)+3)  )- fQ(fup(thisP(ii)+2)  ))
            !         print *, ' '
            !         do kk=1,num_images()
            !             stop 358783
            !         end do
            !     end if
            ! end do  

            !print *, eCons(thisP), (VolNew(thisP) - VolOld(thisP)), dt * eQlat(thisP)
            !print *, ' ', elemYN(thisP,eYN_isSmallDepth), elemYN(thisP,eYN_isZeroDepth)
            !write(*,"(8f12.4)") VolNew(thisP), VolOld(thisP), dt * fQ(fup(thisP+1)), dt * fQ(fup(thisP+3)), dt*fQ(fdn(thisP+2))
            !write(*,"(8f12.4)") VolNew(thisP) - VolOld(thisP), eCons(thisP)
        end if

        !%------------------------------------------------------------------

    end subroutine util_accumulate_volume_conservation
!%==========================================================================
!%==========================================================================
!%
    subroutine util_total_volume_conservation (volume_nonconservation)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes total volume non-conservation of all time-marching elements
        !%------------------------------------------------------------------
        !% Declarations:
            real(8), intent(inout) :: volume_nonconservation
            integer, pointer       :: npack, thisP(:), thisCol
            real(8), save          :: vstore[*]
            integer :: ii
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------
    
        !% for CC elements of any TM
        vstore = zeroR
        thisCol => col_elemP(ep_CC_ALLtm)
        npack   => npack_elemP(thisCol)
        if (npack > 0) then
            thisP => elemP(1:npack,thisCol)
            vstore = vstore + sum(elemR(thisP,er_VolumeConservation))
            !do ii = 1,npack
            !    if (abs(elemR(thisP(ii),er_VolumeConservation)) > 1.0) then
            !        print *, thisP(ii), elemR(thisP(ii),er_VolumeConservation)
            !    end if
            !end do
        end if
       ! print *, 'in util total_volume conservation ',vstore, this_image()

        !% for JM ETM elements
        thisCol =>col_elemP(ep_JM_ALLtm) 
        npack   => npack_elemP(thisCol)
        if (npack > 0) then
            thisP => elemP(1:npack,thisCol)
            vstore = vstore + sum(elemR(thisP,er_VolumeConservation))
            !do ii = 1,npack
            !    if (abs(elemR(thisP(ii),er_VolumeConservation)) > 1.0) then
            !        print *, thisP(ii), elemR(thisP(ii),er_VolumeConservation)
            !    end if
            !end do
        end if

        sync all
        call co_sum(vstore, result_image=1)

        volume_nonconservation = vstore
    
        !%------------------------------------------------------------------
        !% Closing:
    end subroutine util_total_volume_conservation
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_find_elements_in_link &
        (thislinkname, thislink_idx, thislink_image, elemInLink, nElemInLink)
        !%-------------------------------------------------------------------
        !% Description: 
        !% Finds the element indexes for the SWMM link with the name "thislinkname"
        !% stores the element indexes in the array elemInLink and
        !% returns the link index in link%I(:) and the image that host the link.
        !% Not efficiently written, but only purpose is for use in debugging.
        !%-------------------------------------------------------------------
            character(*), intent(in)   :: thislinkname
            integer, intent(inout)     :: elemInLink(:)
            integer, intent(inout)     :: thislink_idx, thislink_image
            integer, intent(inout)     :: nElemInLink
            integer ::  ii, jj
        !%-------------------------------------------------------------------
        elemInLink = 0
        thislink_idx = 0
        thislink_image = 0

        !% find the link idx and image for the input "thislinkname"
        if (this_image() == 1) then
            do ii = 1,size(link%I,dim=1)
                if (link%Names(ii)%str == thislinkname) then
                    write(*,"(A,A,A,i8,A,i6)")'link name ', trim(link%Names(ii)%str), ';  linkIdx= ', ii, ' On image = ',link%I(ii,li_P_image)
                    thislink_idx = ii
                    thislink_image = link%I(ii,li_P_image)
                end if
            end do
        end if
        !% broadcast result to all images
        call co_broadcast(thislink_idx,  source_image=1)
        call co_broadcast(thislink_image,source_image=1)
        sync all

    
        !% for the image that hosts the link, cycle through to find the elements in the link
        if (this_image() == thislink_image) then
            nElemInLink = 0
            do ii = 1,N_elem(this_image())
                if (elemI(ii,ei_link_Gidx_SWMM) == thislink_idx) then
                    nElemInLink = nElemInLink + 1
                    elemInLink(nElemInLink) = ii
                end if
            end do
            write(*,"(A,100i8)") 'elemIdx =',elemInLink(1:nElemInLink)
        end if
        call co_broadcast(elemInLink,source_image=thislink_image)
        sync all

    end subroutine util_find_elements_in_link
!%   
!%==========================================================================    
!%==========================================================================
!%
    subroutine util_find_elements_in_junction_node &
        (thisnodename, thisnode_idx, thisnode_image, elemJM_idx)
        !%------------------------------------------------------------------
        !% Description:
        !% Given the node name, this finds the node index, the image on
        !% which it is an element, and the JM element index
        !%
        !%------------------------------------------------------------------
        !% Declarations:
            character(*), intent(in)   :: thisnodename
            integer, intent(inout)     :: thisnode_idx, thisnode_image
            integer, intent(inout)     :: elemJM_idx
            integer ::  ii, jj
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------

        thisnode_idx = 0
        thisnode_image = 0

        !% find the node idx and image for the input "thisnodename"
        if (this_image() == 1) then
            do ii = 1,size(node%I,dim=1)
                !print *, ii, trim(node%Names(ii)%str)
                if (node%Names(ii)%str == thisnodename) then
                    write(*,"(A,A,A,i8,A,i6)")'node name ', trim(node%Names(ii)%str), ';  nodeIdx= ', ii, '; On image = ',node%I(ii,ni_P_image)
                    thisnode_idx = ii
                    thisnode_image = node%I(ii,ni_P_image)
                    exit !% the first node found should be the JM
                end if
            end do
        end if
        !% broadcast result to all images
        call co_broadcast(thisnode_idx,  source_image=1)
        call co_broadcast(thisnode_image,source_image=1)
        sync all

        elemJM_idx =0
        !% for the image that hosts the node, cycle through to find the element
        if (this_image() == thisnode_image) then
            do ii = 1,N_elem(this_image())
                !print *, ii, elemI(ii,ei_node_Gidx_SWMM)
                if (elemI(ii,ei_node_Gidx_SWMM) == thisnode_idx) then
                    write (*,"(A,i8,A,A)") 'elemIdx= ',ii,'; type = ',reverseKey(elemI(ii,ei_elementType))
                    elemJM_idx = ii
                    exit !% the first should be the correct nodes
                end if
            end do
        end if
        call co_broadcast(elemJM_idx,source_image=thisnode_image)
        sync all

        if (elemJM_idx == 0) then
            write(*,"(A,A)") 'No corresponding JM element found for the node ',trim(node%Names(thisnode_idx)%str)
            write(*,"(A)") 'This implies the node is a nJ1 or nJ2 and is either a boundary or a face'
        end if

        !%------------------------------------------------------------------
        !% Closing:
    end subroutine util_find_elements_in_junction_node
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_find_neighbors_of_CC_element (eIdx, iset)
        !%------------------------------------------------------------------
        !% Description:
        !% For debugging, given an element ID that is CC, find the upstream faces
        !% and elements. Output is in the form
        !% (upstream element, upstream face, element, downstream face, downstream element)\
        !% Note that this will not work across shared faces
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: eIdx
            integer, intent(inout) :: iset(5)
            integer :: ifaceUp, ifaceDn, ielemUp, ielemDn
            character(64) :: subroutine_name = 'util_find_neighbors_of_CC_element'
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (.not. (elemI(eIdx,ei_elementType) == CC)) then
                write(*,"(A,A,A,i8,A,A)") 'in ',trim(subroutine_name), ': element ',eIdx, 'is of type ',reverseKey(elemI(eIdx,ei_elementType))
                write(*,"(A)") 'but this procedure requires CC. Exiting with no result.'
                return
            end if
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------

        ifaceUp = elemI(eIdx,ei_Mface_uL)
        ifaceDn = elemI(eIdx,ei_Mface_dL)
        ielemUp = faceI(ifaceUp,fi_Melem_uL)
        ielemDn = faceI(ifaceDn,fi_Melem_dL)

        iset(1) = ielemUp
        iset(2) = iFaceUp
        iset(3) = eIdx
        iset(4) = ifaceDn
        iset(5) = ielemDn
    
        !%------------------------------------------------------------------
        !% Closing:
    end subroutine util_find_neighbors_of_CC_element
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_find_neighbors_of_JM_element &
        (eIdx, iUpSet, iDnSet, nUpBranch, nDnBranch)
        !%------------------------------------------------------------------
        !% Description:
        !% For debugging, given an element ID of type JM this finds the
        !% faces and neighbor elements
        !% Note that this will not work across shared faces
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: eIdx
            integer, intent(inout) :: iUpSet(max_up_branch_per_node,4)
            integer, intent(inout) :: iDnSet(max_dn_branch_per_node,4)
            integer, intent(inout) :: nUpBranch, nDnBranch
            integer, dimension(max_up_branch_per_node) :: ifaceUp, ielemUp, eIdx_BranchUp
            integer, dimension(max_dn_branch_per_node) :: ifaceDn, ielemDn, eIdx_BranchDn
            integer :: ii
            character(64) :: subroutine_name = 'util_find_neighbors_of_JM_element'
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (.not. (elemI(eIdx,ei_elementType) == JM)) then
                write(*,"(A,A,A,i8,A,A)") 'in ',trim(subroutine_name), ': element ',eIdx, 'is of type ',reverseKey(elemI(eIdx,ei_elementType))
                write(*,"(A)") 'but this procedure requires JM. Exiting with no result.'
                return
            end if
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------

        iUpSet = 0
        iDnSet = 0

        nUpBranch = 0
        do ii=1,max_branch_per_node,2
            if (elemSI(eIdx+ii,esi_JunctionBranch_Exists) > 0) then
                nUpBranch = nUpBranch + 1
                eIdx_BranchUp(nUpBranch) = eIdx + ii
                ifaceUp(nUpBranch) = elemI(eIdx_BranchUp(nUpBranch),ei_Mface_uL)
                ielemUp(nUpBranch) = faceI(ifaceUp(nUpBranch),fi_Melem_uL)
            end if
        end do    

        nDnBranch = 0
        do ii=2,max_branch_per_node,2
            if (elemSI(eIdx+ii,esi_JunctionBranch_Exists) > 0) then
                nDnBranch = nDnBranch + 1
                eIdx_BranchDn(nDnBranch) = eIdx + ii
                ifaceDn(nDnBranch) = elemI(eIdx_BranchDn(nDnBranch),ei_Mface_dL)
                ielemDn(nDnBranch) = faceI(ifaceUp(nDnBranch),fi_Melem_dL)
            end if
        end do    

        do ii=1,nUpBranch
            iUpSet(ii,1) = ielemUp(ii)
            iUpSet(ii,2) = ifaceUp(ii)
            iUpset(ii,3) = eIdx_BranchUp(ii)
            iUpSet(ii,4) = eIdx
        end do

        do ii=1,nDnBranch
            iDnSet(ii,1) = eIdx
            iDnset(ii,2) = eIdx_BranchDn(ii)
            iDnSet(ii,3) = ifaceDn(ii)
            iDnSet(ii,4) = ielemDn(ii)
        end do
            
    
        !%------------------------------------------------------------------
        !% Closing:
    end subroutine util_find_neighbors_of_JM_element
!%
!%==========================================================================    
!%==========================================================================
!%
    subroutine util_CLprint ()
        !%------------------------------------------------------------------
        !% Description:
        !% Used as a debugging write routine
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: ii
            integer, pointer :: fup(:), fdn(:)
            real(8), pointer :: dt
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:
            fup => faceI(:,fi_Melem_uL)
            fdn => faceI(:,fi_Melem_dL)
            dt  => setting%Time%Hydraulics%Dt
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        !% Closing:

        ! if (this_image() == 2) then
        !     print *, ' '
        !     ii = 5428
        !     print *, ii, reverseKey(elemI(ii,ei_elementType))
        !     print *, elemYN(ii,eYN_isSmallDepth), elemYN(ii,eYN_isZeroDepth)
        !     !print *, 'qlat ',elemR(ii,er_FlowrateLateral)
        !     !print *, elemR(ii+1,er_Flowrate), elemR(ii+3, er_Flowrate), elemR(ii+2, er_Flowrate)
        !     print *, faceR(fup(ii+1),fr_Flowrate), faceR(fup(ii+3),fr_Flowrate), faceR(fup(ii+2),fr_Flowrate)
        !     !print *, faceR(fup(ii+1),fr_Flowrate_Conservative), faceR(fup(ii+3), fr_Flowrate_Conservative), faceR(fdn(ii+2), fr_Flowrate_Conservative)
        !     print *, 'vol delta ',elemR(ii,er_Volume) - elemR(ii,er_Volume_N0)
        !     !print *,   dt * faceR(fup(ii+1), fr_Flowrate_Conservative) &
        !     !         + dt * faceR(fup(ii+3), fr_Flowrate_Conservative) &
        !     !        - dt * faceR(fdn(ii+2), fr_Flowrate_Conservative)
        !     print *,  dt * faceR(fup(ii+1), fr_Flowrate) &
        !             + dt * faceR(fup(ii+3), fr_Flowrate) &
        !             - dt * faceR(fdn(ii+2), fr_Flowrate)
        ! end if

        ! print *, elemYN(ieBup2,eYN_isSmallDepth)
        ! print *, elemYN(ieBup2,eYN_isZeroDepth)
        ! write(*,"(A,9(f10.4,A))") 'Q up2 ', &
        !                                 elemR(ieBup2(1),er_Flowrate), ' | ', &
        !                                 faceR(ifBup2(1),fr_Flowrate), ' | ', &
        !                                 elemR(ieBup2(2),er_flowrate), ' | ', &
        !                                 faceR(ifBup2(2),fr_Flowrate), ' | ', &
        !                                 elemR(ieBup2(3),er_flowrate), ' | ', &
        !                                 faceR(ifBup2(3),fr_Flowrate), ' | ', &
        !                                 elemR(ieBup2(4),er_Flowrate), ' | ', &
        !                                 faceR(ifBup2(4),fr_Flowrate), ' | ', &
        !                                 elemR(ieBup2(5),er_Flowrate)

        ! write(*,"(A,9(f10.4,A))") 'Q Con2', &
        !                                 elemR(ieBup2(1),er_FlowrateLateral), ' | ', &
        !                                 faceR(ifBup2(1),fr_Flowrate_Conservative), ' | ', &
        !                                 elemR(ieBup2(2),er_FlowrateLateral), ' | ', &
        !                                 faceR(ifBup2(2),fr_Flowrate_Conservative), ' | ', &
        !                                 elemR(ieBup2(3),er_FlowrateLateral), ' | ', &
        !                                 faceR(ifBup2(3),fr_Flowrate_Conservative), ' | ', &
        !                                 elemR(ieBup2(4),er_FlowrateLateral), ' | ', &
        !                                 faceR(ifBup2(4),fr_Flowrate_Conservative), ' | ', &
        !                                 elemR(ieBup2(5),er_FlowrateLateral)

        ! write(*,"(A,8(f10.4,A),2f10.4)") 'H     ', &
        !                             elemR(ieBup2(1),er_Head),   ' | ', &
        !                             faceR(ifBup2(1),fr_Head_u), ' | ', &
        !                             elemR(ieBup2(2),er_Head),   ' | ', &
        !                             faceR(ifBup2(2),fr_Head_u), ' | ', &
        !                             elemR(ieBup2(3),er_Head),   ' | ', &
        !                             faceR(ifBup2(3),fr_Head_u), ' | ', &
        !                             elemR(ieBup2(4),er_Head),   ' | ', &
        !                             faceR(ifBup2(4),fr_Head_u), ' | ', &
        !                             elemR(ieBup2(5),er_Head), elemR(612,er_Head)

        ! write(*,"(A,8(f10.4,A),2f10.4)") 'D     ', &
        !                             elemR(ieBup2(1),er_Depth),   ' | ', &
        !                             faceR(ifBup2(1),fr_HydDepth_u), ' | ', &
        !                             elemR(ieBup2(2),er_Depth),   ' | ', &
        !                             faceR(ifBup2(2),fr_HydDepth_u), ' | ', &
        !                             elemR(ieBup2(3),er_Depth),   ' | ', &
        !                             faceR(ifBup2(3),fr_HydDepth_u), ' | ', &
        !                             elemR(ieBup2(4),er_Depth),   ' | ', &
        !                             faceR(ifBup2(4),fr_HydDepth_u), ' | ', &
        !                             elemR(ieBup2(5),er_Depth), elemR(612,er_Depth)

    end subroutine util_CLprint    
!%
!%==========================================================================    
!%==========================================================================
!%
        !%------------------------------------------------------------------
        !% Description:
        !%
        !%------------------------------------------------------------------
        !% Declarations:
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------
    
    
        !%------------------------------------------------------------------
        !% Closing:
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module utility
