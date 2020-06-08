! This module calculates the depth in storage nodes
!
!========================================================================== 
!
 module storage


    use adjustments
    use array_index
    use bc
    use data_keys
    use diagnostic
    use face_values
    use globals
    use setting_definition
    use utility

    implicit none

    private

    public :: storage_adjacent_element_average
    public :: storage_depth_volume

    integer :: debuglevel = 0

 contains
!
!========================================================================== 
!==========================================================================
!
subroutine storage_adjacent_element_average &
    (elem2R, elemMR, elemMI, faceI, e2r_data, eMr_out)
!
! this computes the average of values for all the elements upstream and
! downstream. Note that this should ONLY be use
! in setup routines (i.e. when initializing storage). This violates
! the "no-neighbor" rule and is time-consuming because it requires double
! mapping
! 
 character(64) :: subroutine_name = 'storage_adjacent_element_average'

 real,      target,     intent(in out)  :: elemMR(:,:)
 real,                  intent(in)      :: elem2R(:,:)
 integer,               intent(in)      :: elemMI(:,:), faceI(:,:)
 integer,               intent(in)      :: e2r_data, eMr_out
 
 real,      pointer :: Uvalue(:), Dvalue(:)
 integer :: eMr_tUp, eMr_tDn
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 eMr_tUp = eMr_Temp(next_eMr_temparray)
 Uvalue  => elemMR(:,eMr_tUp)
 next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)
 
 eMr_tDn = eMr_Temp(next_eMr_temparray)
 Dvalue  => elemMR(:,eMr_tDn)
 next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp) 
 
 Uvalue = zeroR
 call storage_summation_from_adjacent_elements_one_direction &
    (eMr_tUp, elem2R, elemMR, elemMI, faceI, &
     upstream_face_per_elemM, eMi_nfaces_u, eMi_MfaceUp, fi_Melem_u, e2r_data)

 Dvalue = zeroR
 call storage_summation_from_adjacent_elements_one_direction &
    (eMr_tDn, elem2R, elemMR, elemMI, faceI, &
     dnstream_face_per_elemM, eMi_nfaces_d, eMi_MfaceDn, fi_Melem_d, e2r_data)
 
 where (elemMI(:,eMi_elem_type) == eStorage)
    elemMR(:,eMr_out) = (Uvalue + Dvalue) / real( elemMI(:,eMi_nfaces_u) + elemMI(:,eMi_nfaces_d) )
 endwhere
 
 Dvalue = nullvalueR
 Uvalue = nullvalueR
 nullify(Dvalue,Uvalue)
 next_eMr_temparray = next_eMr_temparray-2
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine storage_adjacent_element_average
 !
!==========================================================================
!==========================================================================
!
 subroutine storage_summation_from_adjacent_elements_one_direction &
    (eMr_sumvalue, elem2R, elemMR, elemMI, faceI, &
     dir_face_per_elemM, eMi_nfaces_dir, eMi_MfaceDir, fi_Melem_dir, e2r_data)
      
 character(64) :: subroutine_name = 'storage_summation_from_adjacent_elements_one_direction'
!
! computes the sum of all the elements adjacent to a junctio in either
! the upstream or downstream direction
!
! THIS SHOULD ONLY BE USED IN SETUP AND INITIAL CONDITION ROUTINES
! 
 integer,           intent(in)      :: eMr_sumvalue
 real,      target, intent(in out)  :: elemMR(:,:)
 real,              intent(in)      :: elem2R(:,:)
 integer,   target, intent(in)      :: elemMI(:,:), faceI(:,:)
 
 integer,   intent(in)  :: dir_face_per_elemM, eMi_nfaces_dir, eMi_MfaceDir(:)
 integer,   intent(in)  :: fi_Melem_dir, e2r_data
 
 integer,   pointer :: tface, telem
 real   :: thisvalue(dir_face_per_elemM)
 integer :: mm, ii
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 
 do ii=1,N_elemM  
    thisvalue = zeroR
    do mm=1,dir_face_per_elemM
        if ((elemMI(ii,eMi_nfaces_dir) >= mm) .and. &
            (elemMI(ii,eMi_elem_type) == eStorage)) then
            !% the face on a branch in the direction specified by dir
            tface => elemMI(ii,eMi_MfaceDir(mm))
            !% the element upstream of the face
            telem => faceI(tface,fi_Melem_dir)
            !% the value at the element
            thisvalue(mm) = elem2R(telem,e2r_data)
        endif
    enddo
    elemMR(ii,eMr_sumvalue) = sum(thisvalue)
 enddo

 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine storage_summation_from_adjacent_elements_one_direction
!
!==========================================================================
!==========================================================================
!
 subroutine storage_depth_volume (elemMR, elemMI) 

 character(64) :: subroutine_name = 'storage_depth_volume'
!
! computes the volume in a storage unit by functional or tabular relationship
! 
 real,      target, intent(in out)  :: elemMR(:,:)
 integer,   target, intent(in)      :: elemMI(:,:)

 integer,   pointer :: cindx
 real,      pointer :: sVolume, sDepth, sArea, aConst, aCoeff, aExpon
 
 integer :: mm, ii
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 

 where (elemMI(:,eMi_elem_type) == eStorage)
        elemMR(:,eMr_HydDepth) = elemMR(:,eMr_Eta) - elemMR(:,eMr_Zbottom)
 endwhere

 do ii=1,N_elemM
    cindx    => elemMI(ii,eMi_curve_type)
    sDepth   => elemMR(ii,eMr_HydDepth)
    sVolume  => elemMR(ii,eMr_Volume)
    sArea    => elemMR(ii,eMr_SurfArea)
    aConst   => elemMR(ii,eMr_StorageConstant)
    aCoeff   => elemMR(ii,eMr_StorageCoeff)
    aExpon   => elemMR(ii,eMr_StorageExponent)

    if ((elemMI(ii,eMi_elem_type) == eStorage) .and. &
        (sDepth .LT. zeroR) ) then

        sDepth = zeroR

    elseif ((elemMI(ii,eMi_elem_type) == eStorage) .and. &
            (sDepth .GE. elemMR(ii,eMr_FullDepth))) then

        sDepth = elemMR(ii,eMr_FullDepth)

    endif
            
    select case (cindx)
        case (1)
        !% case where the relationship between depth and surface area is Functional
        if ((elemMI(ii,eMi_elem_type) == eStorage) .and. &
            (sDepth == zeroR) ) then

            sVolume = zeroR
            sArea   = zeroR
            
        elseif ((elemMI(ii,eMi_elem_type) == eStorage) .and. &
                (sDepth .GT. zeroR) ) then

            sVolume = aConst * sDepth + aCoeff / ((aExpon + oneR) &
                        * sDepth ** (aExpon + oneR))
            sArea   = aConst + aCoeff * sDepth ** aExpon

        endif

        case(2)
        !% case where the relationship between depth and surface area is Tabular
        if ((elemMI(ii,eMi_elem_type) == eStorage) .and. &
            (sDepth == zeroR) ) then

            print*, 'error: tabular storage curve type in not handeled at this moment'
            stop

        elseif ((elemMI(ii,eMi_elem_type) == eStorage) .and. &
                (sDepth .GT. zeroR) ) then

            print*, 'error: tabular storage curve type in not handeled at this moment'
            stop
        endif
    end select
 enddo

 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine storage_depth_volume
!
!==========================================================================
!==========================================================================
!
 end module storage