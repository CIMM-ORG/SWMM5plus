! module checking
!
! Perform various checks of the simulation setup
! In general, these are procedures that should only be called once.
!
!==========================================================================
!
 module checking
! 
! contains generalized checking algorithms
!
    use array_index
    use data_keys
    use globals
    
    implicit none
    private
    
    public :: checking_consistency
    public :: checking_smallvolume_consistency
    
    integer, private :: debuglevel = 0
    
 contains
!
!==========================================================================
!==========================================================================
!
 subroutine checking_consistency
!
! Checks consistency of parameter dimensions
!   The number of upstream and downstream faces per element in an elemM
!   array must be consistent with the total faces per element.
!   These are parameter dimensions that are separately set (rather than
!   variables) to prevent accidental change during the run, which could
!   be a catastrophic bug.
!
! TODO: consider adding checking for other parameters
!
    
 character(64) :: subroutine_name = 'checking_consistency'
 integer :: aa

!--------------------------------------------------------------------------  
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name  
 
 aa = upstream_face_per_elemM + dnstream_face_per_elemM
 if (aa /= face_per_elemM) then
    print *, 'upstream_face_per_elemM = ',upstream_face_per_elemM
    print *, 'dnstream_face_per_elemM = ',dnstream_face_per_elemM
    print *, 'face_per_elemM          = ',face_per_elemM
    print *, 'error: inconsistent parameters for faces per element in ',subroutine_name
    stop
 endif
 
 if (nullvalueI > 0) then
    print *, 'error: nullvalueI <=0 is required in ',subroutine_name
    STOP
 endif
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name  
 end subroutine checking_consistency
!
!==========================================================================
!==========================================================================
!
 subroutine checking_smallvolume_consistency &
    (elem2R, elemMR)
!
! checks the consistency of the smallvolume settings (if used) and the
! zero value settings
! 
 character(64) :: subroutine_name = 'checking_smallvolume_consistency'
 
 real,  intent(in)  :: elem2R(:,:), elemMR(:,:) 
 integer :: ii
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 if (.not. setting%SmallVolume%UseSmallVolumes) return
 
!  do ii= 1 , size(elem2R(:,e2r_SmallVolume),1)
!     print*, "ii=", ii, "elem2R(ii,e2r_SmallVolume)", elem2R(ii,e2r_SmallVolume)
!  enddo
 
 if (any(elem2R(:,e2r_SmallVolume) <= setting%ZeroValue%Volume) .or. &
     any(elemMR(:,eMr_SmallVolume) <= setting%ZeroValue%Volume)        ) then
    print *, 'setting%ZeroValue%Volume        = ',setting%ZeroValue%Volume
    print *, 'setting%SmallVolume%DepthCutoff = ',setting%SmallVolume%DepthCutoff
    print *, 'user error: the setting%ZeroValue%Volume is too large.'
    print *, 'It must be smaller than the smallest elemR(:,SmallVolume)'
    print *, 'which is computed using setting%SmallVolume%DepthCutoff.'
    print *, 'Recommend increasing SmallVolume%DepthCutoff or reducing the '
    print *, 'ZeroValue%Volume'
    stop
    
 endif
 
 if (setting%SmallVolume%DepthCutoff <= setting%ZeroValue%Depth ) then
    print *,'setting%SmallVolume%DepthCutoff = ',setting%SmallVolume%DepthCutoff
    print *,'setting%ZeroValue%Depth         = ',  setting%ZeroValue%Depth
    print *,'user error: setting%SmallVolume%DepthCutoff <= setting%ZeroValue%Depth'
    print *,'This is inconsistent. Increase the former or decrease the latter.'
    stop
 endif

 if (setting%SmallVolume%MinimumTopwidth < setting%ZeroValue%Topwidth ) then
    print *,'setting%SmallVolume%MinimumTopwidth = ',setting%SmallVolume%MinimumTopwidth
    print *,'setting%ZeroValue%Topwidth          = ',  setting%ZeroValue%Topwidth
    print *,'user error: setting%SmallVolume%MinimumTopwidth <= setting%ZeroValue%Topwidth'
    print *,'This is inconsistent. Increase the former or decrease the latter.'
    stop
 endif
 

 if (setting%SmallVolume%MinimumArea < setting%ZeroValue%Area ) then
    print *,'setting%SmallVolume%MinimumArea = ',setting%SmallVolume%MinimumArea
    print *,'setting%ZeroValue%Area          = ',  setting%ZeroValue%Area
    print *,'user error: setting%SmallVolume%MinimumArea <= setting%ZeroValue%Area'
    print *,'This is inconsistent. Increase the former or decrease the latter.'
    stop
 endif
 
 if (setting%SmallVolume%MinimumPerimeter < setting%ZeroValue%Topwidth ) then
    print *,'setting%SmallVolume%MinimumPerimeter = ',setting%SmallVolume%MinimumPerimeter
    print *,'setting%ZeroValue%Topwidth          = ',   setting%ZeroValue%Topwidth
    print *,'user error: setting%SmallVolume%MinimumPerimeter <= setting%ZeroValue%Topwidth'
    print *,'This is inconsistent. Increase the former or decrease the latter.'
    stop
 endif 
 
 if (setting%SmallVolume%MinimumHydRadius < setting%ZeroValue%Depth ) then
    print *,'setting%SmallVolume%MinimumHydRadius = ',setting%SmallVolume%MinimumHydRadius
    print *,'setting%ZeroValue%Depth              = ',  setting%ZeroValue%Depth
    print *,'user error: setting%SmallVolume%MinimumHydRadius <= setting%ZeroValue%Depth'
    print *,'This is inconsistent. Increase the former or decrease the latter.'
    stop
 endif  
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine checking_smallvolume_consistency
!
!========================================================================== 
! END OF MODULE checking
!==========================================================================
 end module checking
