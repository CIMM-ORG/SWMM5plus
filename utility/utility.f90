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
    public :: util_find_elements_in_link
    public :: util_accumulate_volume_conservation
    public :: util_total_volume_conservation

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
    subroutine util_find_elements_in_link (thislinkname, thislink_idx, thislink_image, elemInLink)
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
            integer ::  ii, jj, ecount
        !%-------------------------------------------------------------------
        elemInLink = 0
        thislink_idx = 0
        thislink_image = 0

        !% find the link idx and image for the input "thislinkname"
        if (this_image() == 1) then
            do ii = 1,size(link%I,dim=1)
                if (link%Names(ii)%str == thislinkname) then
                    print *, ii, trim(link%Names(ii)%str), ' ' ,trim(thislinkname), ' ',link%I(ii,li_P_image)
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
            ecount = 1
            do ii = 1,N_elem(this_image())
                if (elemI(ii,ei_link_Gidx_SWMM) == thislink_idx) then
                    elemInLink(ecount) = ii
                    print *, ii, ' is elem #'
                    ecount = ecount + 1
                end if
            end do
        end if
        call co_broadcast(elemInLink,source_image=thislink_image)
        sync all

    end subroutine util_find_elements_in_link
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
            integer :: ii
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:
            thisColCC => col_elemP(ep_CC_Q_NOTsmalldepth)
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
            !% --- sum of the net inflow and lateral flow should be the change in volume from head
            !% --- for output, use an accumulator
            eCons(thisP) = eCons(thisP)                                            &
                         + dt * ( fQ(fup(thisP)) - fQ(fdn(thisP)) + eQlat(thisP) ) &
                         - (VolNew(thisP) - VolOld(thisP))

            !% --- for debugging, switch to using non-cumulative            
            !eCons(thisP) = dt * (fQ(fup(thisP)) - fQ(fdn(thisP)) + eQlat(thisP)) &
            !              - (VolNew(thisP) - VolOld(thisP))      
                          
            ! do ii = 1,size(thisP)
            !     if (abs(eCons(thisP(ii))) > 1.0) then
            !         print *, ii, thisP(ii), this_image()
            !         print *,  eCons(thisP(ii))
            !         print *, fQ(fup(thisP(ii))), fQ(fdn(thisP(ii))), eQlat(thisP(ii))
            !         print *, VolNew(thisP(ii)), VolOld(thisP(ii))
            !         stop 358783
            !     end if
            ! end do              
        end if

        !% HACK -- need an equivalent of the ep_CC_Q_NOTsmall depth for JM
        !% to make this work

        ! !% for the JM elements
        ! npack   => npack_elemP(thisColJM)
        ! if (npack > 0) then
        !     thisP => elemP(1:npack,thisColJM)
        !     eCons(thisP) = eCons(thisP) + eQlat(thisP) - (VolNew(thisP) - VolOld(thisP))
        !     do ii=1,max_branch_per_node,2
        !         eCons(thisP) = eCons(thisP)                                           &
        !                      + dt * ( real(BranchExists(thisP),8) * fQ(fup(thisP+ii)) &
        !                     - real(BranchExists(thisP),8) * fQ(fdn(thisP+ii+1)) )
        !     end do
        ! end if

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
            do ii = 1,npack
                if (abs(elemR(thisP(ii),er_VolumeConservation)) > 1.0) then
                    print *, thisP(ii), elemR(thisP(ii),er_VolumeConservation)
                end if
            end do
        end if
       ! print *, 'in util total_volume conservation ',vstore, this_image()

        !% HACK -- the JM elements aren't finished yet
        ! !% for JM ETM elements
        ! thisCol =>col_elemP(ep_JM_ETM) 
        ! npack   => npack_elemP(thisCol)
        ! if (npack > 0) then
        !     thisP => elemP(1:npack,thisCol)
        !     vstore = vstore + sum(elemR(thisP,er_VolumeConservation))
        ! end if

        ! !% for JJ AC elements
        ! thisCol =>col_elemP(ep_JM_AC) 
        ! npack   => npack_elemP(thisCol)
        ! if (npack > 0) then
        !     thisP => elemP(1:npack,thisCol)
        !     vstore = vstore + sum(elemR(thisP,er_VolumeConservation))
        ! end if

        ! sync all
        call co_sum(vstore, result_image=1)

        volume_nonconservation = vstore
    
        !%------------------------------------------------------------------
        !% Closing:
    end subroutine util_total_volume_conservation
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module utility
