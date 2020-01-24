! This module calculates the flow in a Weir Element
!
!========================================================================== 
!
 module weir


    use adjustments
    use array_index
    use bc
    use data_keys
    use diagnostic
    use globals
    use setting_definition
    use utility

    implicit none

    private

    public :: weir_step

    integer :: debuglevel = 0

 contains
!
!========================================================================== 
!==========================================================================
!
 subroutine weir_step &
    (e2r_Volume_old, e2r_Velocity_old, eMr_Volume_old, eMr_Velocity_old, &
     e2r_Volume_new, e2r_Velocity_new, eMr_Volume_new, eMr_Velocity_new, &
     elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     thiscoef)
!
 character(64) :: subroutine_name = 'weir_step'
 
! indexes for old/new volume and velocity storage
 integer,   intent(in) :: e2r_Volume_old, e2r_Velocity_old
 integer,   intent(in) :: eMr_Volume_old, eMr_Velocity_old
 integer,   intent(in) :: e2r_Volume_new, e2r_Velocity_new
 integer,   intent(in) :: eMr_Volume_new, eMr_Velocity_new

 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,           intent(in out)  :: elem2YN(:,:), elemMYN(:,:)
 real,              intent(in)      :: thiscoef

 real,  pointer ::  volume2old(:), volume2new(:), velocity2old(:), velocity2new(:)
 real,  pointer ::  volumeMold(:), volumeMnew(:), velocityMold(:), velocityMnew(:)
 real,  pointer ::  wCrest(:), wCrown(:), wZbottom(:), EffectiveHead(:)

 real,  pointer ::  wCoeff, wWidth, wHeight, wSideSlope, wInletoffset    !Weir discharge coefficient, Width, Height, sideslope

 integer, pointer                   :: iup(:), idn(:), dir(:)
 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


 volume2old   => elem2R(:,e2r_Volume_old)
 volume2new   => elem2R(:,e2r_Volume_new)

 velocity2old => elem2R(:,e2r_Velocity_old)
 velocity2new => elem2R(:,e2r_Velocity_new)

 volumeMold   => elemMR(:,eMr_Volume_old)
 volumeMnew   => elemMR(:,eMr_Volume_new)

 velocityMold => elemMR(:,eMr_Velocity_old)
 velocityMnew => elemMR(:,eMr_Velocity_new)

 wCrown       => elem2R(:,e2r_Depth)
 wZbottom     => elem2R(:,e2r_Zbottom)

!%  temporary space for elem2
 EffectiveHead => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 wCrest        => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 dir           => elem2I(:,e2i_Temp(next_e2i_temparray))
 next_e2i_temparray = utility_advance_temp_array (next_e2i_temparray,e2i_n_temp)

!%  zero temporary arrays
 EffectiveHead = zeroR
 dir = zeroI


!%  pointers for weir settings
 wWidth         => setting%Weir%WeirWidth
 wCoeff         => setting%Weir%WeirDischargeCoeff
 wHeight        => setting%Weir%WeirHeight
 wSideSlope     => setting%Weir%WeirSideSlope
 wInletoffset   => setting%Weir%WeirInletOffset 

 ! where ( ( elem2I(:,e2i_elem_type) == eWeir ) )
 !    wCrest = wInletoffset
 !    wCrown = wCrest + wHeight
 ! endwhere

 call weir_effective_head &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     wCrest, wCrown, dir, EffectiveHead)
 
 call weir_effective_length &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, & 
     wHeight, thiscoef)

 call weir_flow &
    (volume2old, velocity2old, volumeMold, velocityMold, &
     volume2new, velocity2new, volumeMnew, velocityMnew, &
     elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     wWidth, wHeight, wCoeff, wSideSlope, dir, EffectiveHead, &
     thiscoef)

 
 ! print *,'====================================================='
 ! print *, faceR(:,fr_HydDepth_d), 'H d'
 ! print *,'====================================================='
 ! print *, faceR(:,fr_HydDepth_u), 'H u'
 ! print *,'====================================================='
 ! print *, 'Weir Effective Head'
 ! print *, EffectiveHead



! stop
 ! release temporary arrays
 EffectiveHead  = nullvalueR
 wCrest         = nullvalueR

 dir            = nullvalueI
 nullify(EffectiveHead, wCrest)
 nullify(dir)
 next_e2r_temparray = next_e2r_temparray - 2
 next_e2i_temparray = next_e2i_temparray - 1

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_step
!
!==========================================================================
!==========================================================================
!
subroutine weir_effective_head &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     wCrest, wCrown, dir, EffectiveHead)
!
 character(64) :: subroutine_name = 'weir_effective_head'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,           intent(in out)  :: elem2YN(:,:), elemMYN(:,:)

 real,  pointer   ::  wCrest(:), wCrown(:), EffectiveHead(:)

 real,  pointer   ::  fEdn(:), fEup(:)
 integer, pointer ::  iup(:), idn(:), dir(:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!%  pointers for convenience in notation
 fEdn  => faceR(:,fr_Eta_d)
 fEup  => faceR(:,fr_Eta_u)

!% Calculate Effective Head (Assuming the flow is U/S to D/S needed to be fixed for flow revarsal
 iup   => elem2I(:,e2i_Mface_u)
 idn   => elem2I(:,e2i_Mface_d)

 where      ( (elem2I(:,e2i_elem_type) == eWeir) .and.  &
              (fEdn(iup) .GE. fEup(idn)) )

            dir = oneI
            EffectiveHead = fEdn(iup)  - fEup(idn)

 elsewhere ( (elem2I(:,e2i_elem_type) == eWeir) .and.  &
              (fEdn(iup) .LT. fEup(idn)) )

            dir = -oneI
            EffectiveHead = fEup(idn) - fEdn(iup)
 endwhere

 ! where     ( (elem2I(:,e2i_elem_type) == eWeir ) .and. &
 !             (fEup(idn) .GT. fEdn(iup)) )
 !    EffectiveHead  = fEup(idn) - wCrest

 ! elsewhere ( (elem2I(:,e2i_elem_type) == eWeir ) .and. &
 !             (fEup(idn) .LT. fEdn(iup)) )
 !    EffectiveHead  = wCrest - fEdn(iup)

 ! elsewhere ( (elem2I(:,e2i_elem_type) == eWeir ) .and. &
 !             (fEup(idn) .LT. wCrest) .and. &
 !             (fEdn(iup) .LT. wCrest) )
 !    EffectiveHead = zeroR

 ! elsewhere ( (elem2I(:,e2i_elem_type) == eWeir ) .and. &
 !             (fEup(idn) .GT. wCrown) .and. &
 !             (fEup(idn) .GT. fEdn(iup) ) )
 !    EffectiveHead = wCrown - wCrest

 ! elsewhere ( (elem2I(:,e2i_elem_type) == eWeir ) .and. &
 !             (fEdn(iup) .GT. wCrown) .and. &
 !             (fEdn(iup) .GT. fEup(idn) ) )
 !    EffectiveHead =  wCrest - wCrown

 ! elsewhere ( (elem2I(:,e2i_elem_type) == eWeir ) .and. &
 !             (fEup(idn) .EQ. fEdn(iup)) )
 !    EffectiveHead =  zeroR

!Need a fix for surcharge

 ! print*, EffectiveHead
 ! stop

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_effective_head
!
!========================================================================== 
!==========================================================================
!
subroutine weir_effective_length &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, wHeight, &
     thiscoef)
!
 character(64) :: subroutine_name = 'weir_effective_length'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,           intent(in out)  :: elem2YN(:,:), elemMYN(:,:)
 real,  pointer                     :: wLength(:)
 real,  pointer                     :: wHeight               !Weir Height
 real,              intent(in)      :: thiscoef

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

 wLength => elem2R(:,e2r_Length)

! Calculate Effective Length (This part is straight up from SWMM source code)
 where ( (elem2I(:,e2i_elem_type) == eWeir ) )
    wLength  = twoR*thiscoef*sqrt(grav*wHeight)
    wLength  = min(wLength, 200.0)
 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_effective_length
!
!==========================================================================
!==========================================================================
!
 subroutine weir_flow &
    (volume2old, velocity2old, volumeMold, velocityMold, &
     volume2new, velocity2new, volumeMnew, velocityMnew, &
     elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     wWidth, wHeight, wCoeff, wSideSlope, dir, EffectiveHead, &
     thiscoef)
!
 character(64) :: subroutine_name = 'weir_flow'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,           intent(in out)  :: elem2YN(:,:), elemMYN(:,:)
 real,              intent(in)      :: thiscoef

 
 real,  pointer ::  volume2old(:), volume2new(:), velocity2old(:), velocity2new(:)
 real,  pointer ::  volumeMold(:), volumeMnew(:), velocityMold(:), velocityMnew(:)
 real,  pointer ::  EffectiveHead(:), depth_element(:)
 real,  pointer ::  wCoeff, wWidth, wHeight, wSideSlope     !Weir discharge coefficient, Width, Height

 integer, pointer :: dir(:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

 depth_element => elem2R(:,e2r_Depth)

 where ( (elem2I(:,e2i_elem_type) == eWeir ).and. &
         (elem2I(:,e2i_geometry)  == eVnotchWeir) )

    ! Volume is weir flow equation * dt (this case dt = thiscoef)
    volume2new   = thiscoef * wCoeff * wSideSlope * EffectiveHead ** 2.5
  
    velocity2new = dir * wCoeff * sqrt(abs(EffectiveHead)) / wSideSlope

 endwhere
 ! print *,'====================================================='
 ! print *, faceR(:,fr_HydDepth_d), 'H d'
 ! print *,'====================================================='
 ! print *, faceR(:,fr_HydDepth_u), 'H u'
 ! print *,'====================================================='
 ! print *, 'Weir Effective Head'
 ! print *, EffectiveHead
 ! print*, '+++++++++++++++++++++++++'
 ! print *, 'Volume new'
 ! print *, volume2new
 ! print *, 'Velocity new'
 ! print *, velocity2new
 ! print *, 'RK2 coeff'
 ! print *, thiscoef
 ! print *, 'dir'
 ! print *, dir
 ! print *,'END OF STEP'

!%  pointers for volume and velocity storage (updating)
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_flow
!
!==========================================================================
!==========================================================================
!
 end module weir