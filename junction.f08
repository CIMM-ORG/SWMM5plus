! module junction
!
! Procedures associated with junctions that have more than 2 faces
!
!==========================================================================
!
 module junction
! 
    use array_index
    use data_keys
    use globals
    use setting_definition

    
    implicit none
    
    private
    
    public  :: junction_geometry_from_branches
    public  :: junction_branches_assigned_to_faces

    integer :: debuglevel = 0
    
 contains
!
!========================================================================== 
!==========================================================================
!
 subroutine junction_geometry_from_branches  &
    (elemMR, elemMI)
!
! compute characteristic values for junction geometry
! using branch values
!
 character(64) :: subroutine_name = 'junction_geometry_from_branches'
 
 real,      intent(in out)  :: elemMR(:,:)
 
 integer,   intent(in)      :: elemMI(:,:)  
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 if (first_elemM_index /= 1) then
    print *, 'Check subroutine ',subroutine_name
    print *, 'functionality with first_elemM_index = ',first_elemM_index,&
             'is not guaranteed.'    
    !print *, 'Example using maxval:'         
    !print *, maxval( elemMR(:,eMr_BreadthScaleAll) , 2)
    !print *, maxval( elemMR(first_elemM_index:first_elemM_index+N_elemM-1,eMr_BreadthScaleAll) , 2)
    print *, 'Also need to examine topwidth vs. breadthscale'
    stop
 endif

 ! characteristic topwidth as the max of any branch
 elemMR(:,eMr_BreadthScale) = maxval( elemMR(:,eMr_BreadthScaleAll) , 2)

 ! volume as the sum of all the branch volumes 
 elemMR(:,eMr_Volume) = sum( (elemMR(:,eMr_AreaAll) * elemMR(:,eMr_LengthAll)) ,2)
 
 ! length as the sum of the average upstream and average downstream branch lengths
 elemMR(:,eMr_Length) = ( sum(elemMR(:,eMr_LengthUp),2) / elemMI(:,eMi_nfaces_u) ) &
                      + ( sum(elemMR(:,eMr_LengthDn),2) / elemMI(:,eMi_nfaces_d) )
 
 ! representative area based on volume/length
 elemMR(:,eMr_Area)  = elemMR(:,eMr_Volume) / elemMR(:,eMR_Length)
 
 ! perimeter based on rectangular cross section
 elemMR(:,eMr_Perimeter) = elemMR(:,eMr_BreadthScale) + twoR * (elemMR(:,eMr_Eta) - elemMR(:,eMr_Zbottom))
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine junction_geometry_from_branches
!
!==========================================================================
!==========================================================================
!
 subroutine junction_branches_assigned_to_faces &
    (faceI, elemMI)
 
 character(64) :: subroutine_name = 'junction_branches_assigned_to_faces'
 
 integer,   intent(in out)  :: faceI(:,:)
 
 integer,   intent(in)      :: elemMI(:,:)
 
 integer :: mm
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 do mm=1,upstream_face_per_elemM
    where (elemMI(:,eMi_nfaces_u) >= mm)
        faceI( elemMI(: , eMi_MfaceUp(mm)),fi_branch_d) = mm
    endwhere
 end do
 
 
 do mm=1,dnstream_face_per_elemM
    where (elemMI(:,eMi_nfaces_d) >= mm)
        faceI( elemMI(: , eMi_MfaceDn(mm)),fi_branch_u) = mm
    endwhere
 end do
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine junction_branches_assigned_to_faces
!
!========================================================================== 
! END OF MODULE junction
!==========================================================================
 end module junction