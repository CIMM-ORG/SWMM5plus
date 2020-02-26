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
    use face_values
    use globals
    use setting_definition
    use utility

    implicit none

    private

    public :: weir_step
    public :: weir_freesurface_elevation
    public :: weir_provisional_geometry

    integer :: debuglevel = 0

 contains
!
!========================================================================== 
!==========================================================================
!
 subroutine weir_step &
    (e2r_Volume_old, e2r_Velocity_old, eMr_Volume_old, eMr_Velocity_old, &
     e2r_Volume_new, e2r_Velocity_new, eMr_Volume_new, eMr_Velocity_new, &
     elem2R, elemMR, faceI, faceR, faceYN, elem2I, elemMI, elem2YN, &
     elemMYN, thiscoef)
!
 character(64) :: subroutine_name = 'weir_step'
 
! indexes for old/new volume and velocity storage
 integer,   intent(in) :: e2r_Volume_old, e2r_Velocity_old
 integer,   intent(in) :: eMr_Volume_old, eMr_Velocity_old
 integer,   intent(in) :: e2r_Volume_new, e2r_Velocity_new
 integer,   intent(in) :: eMr_Volume_new, eMr_Velocity_new

 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 integer,           intent(in out)  :: faceI(:,:)
 real,      target, intent(in out)  :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,           intent(in out)  :: elem2YN(:,:), elemMYN(:,:)
 logical,           intent(in out)  :: faceYN(:,:)
 real,              intent(in)      :: thiscoef

 real,  pointer     ::  volume2old(:), volume2new(:), velocity2old(:), velocity2new(:)
 real,  pointer     ::  volumeMold(:), volumeMnew(:), velocityMold(:), velocityMnew(:)
 real,  pointer     ::  wFlow(:), wCrest(:), wCrown(:), wZbottom(:), wEta(:), EffectiveHead(:)

 real,  pointer     ::  wCoeff, wWidth, wHeight, wSideSlope, wInletoffset    !Weir discharge coefficient, Width, Height, sideslope
 real,  pointer     ::  fEdn(:), fEup(:)

 integer, pointer   ::  iup(:), idn(:), dir(:)
 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!%  pointers for convenience in notation
 fEdn         => faceR(:,fr_Eta_d)
 fEup         => faceR(:,fr_Eta_u)

 volume2old   => elem2R(:,e2r_Volume_old)
 volume2new   => elem2R(:,e2r_Volume_new)

 velocity2old => elem2R(:,e2r_Velocity_old)
 velocity2new => elem2R(:,e2r_Velocity_new)

 volumeMold   => elemMR(:,eMr_Volume_old)
 volumeMnew   => elemMR(:,eMr_Volume_new)

 velocityMold => elemMR(:,eMr_Velocity_old)
 velocityMnew => elemMR(:,eMr_Velocity_new)

 wflow        => elem2R(:,e2r_Flowrate)  
 wZbottom     => elem2R(:,e2r_Zbottom)
 wEta         => elem2R(:,e2r_eta)

 iup          => elem2I(:,e2i_Mface_u)
 idn          => elem2I(:,e2i_Mface_d)

!%  temporary space for elem2
 EffectiveHead => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 wCrest        => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 wCrown       => elem2R(:,e2r_Temp(next_e2r_temparray))
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
 
 where      ( (elem2I(:,e2i_elem_type) == eWeir) )
            wCrest =  wInletoffset + wZbottom
            wCrown =  wCrest + wHeight
 endwhere 

 call weir_provisional_geometry &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN) 
    
 call weir_freesurface_elevation &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     fEdn, fEup, iup, idn, wEta)

 call weir_effective_head &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     wCrest, wCrown, wEta, fEup, fEdn, iup, idn, dir, EffectiveHead)
 
 call weir_effective_length &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, & 
     wHeight, thiscoef)

 call weir_flow &
    (volume2old, velocity2old, volumeMold, velocityMold, &
     volume2new, velocity2new, volumeMnew, velocityMnew, &
     wFlow, elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, &
     elemMYN, wWidth, wHeight, wCoeff, wSideSlope, dir,     &
     EffectiveHead, thiscoef)

    ! print*, wFlow, 'wFlow before'
    ! print*,'**************************************'
    ! print*, faceR(:,fr_Flowrate), 'fr_Flowrate'
    ! print*,'**************************************'
    ! print*, faceR(:,fr_Velocity_u), 'fr_Velocity_u'
    ! print*,'**************************************'
    ! print*, faceR(:,fr_Velocity_d), 'fr_Velocity_d'
    ! print*,'**************************************'
 ! call flow_interp_for_upstream_weir_face (elem2R, faceR, faceI, faceYN)

 ! call flow_interp_for_downstream_weir_face (elem2R, faceR, faceI, faceYN)

    ! print*,'+++++++++++++++++++++++++++++++++++++'
    ! print*, faceR(:,fr_Flowrate), 'fr_Flowrate'
    ! print*,'+++++++++++++++++++++++++++++++++++++'
    ! print*, faceR(:,fr_Velocity_u), 'fr_Velocity_u'
    ! print*,'+++++++++++++++++++++++++++++++++++++'
    ! print*, faceR(:,fr_Velocity_d), 'fr_Velocity_d'
    ! print*,'+++++++++++++++++++++++++++++++++++++'

 ! release temporary arrays
 EffectiveHead  = nullvalueR
 wCrest         = nullvalueR
 dir            = nullvalueI

 nullify(EffectiveHead, wCrest, wCrown)
 nullify(dir)
 next_e2r_temparray = next_e2r_temparray - 3
 next_e2i_temparray = next_e2i_temparray - 1

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_step
!
!==========================================================================
!==========================================================================
!
subroutine weir_freesurface_elevation &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     fEdn, fEup, iup, idn, wEta)
!
 character(64) :: subroutine_name = 'weir_freesurface_elevation'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,           intent(in out)  :: elem2YN(:,:), elemMYN(:,:)

 real,  pointer   ::  fEdn(:), fEup(:)
 real,  pointer   ::  wEta(:)
 integer, pointer ::  iup(:), idn(:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


 where      ( (elem2I(:,e2i_elem_type) == eWeir) )
            wEta = max(fEdn(iup), fEup(idn))
 endwhere
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_freesurface_elevation
!
!========================================================================== 
!==========================================================================
!
subroutine weir_effective_head &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     wCrest, wCrown, wEta, fEup, fEdn, iup, idn, dir, EffectiveHead)
!
 character(64) :: subroutine_name = 'weir_effective_head'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,           intent(in out)  :: elem2YN(:,:), elemMYN(:,:)

 real,  pointer   ::  wCrest(:), wCrown(:), EffectiveHead(:), wEta(:)

 real,  pointer   ::  fEup(:), fEdn(:), nominalEup(:), nominalEdn(:)
 integer, pointer ::  iup(:), idn(:), dir(:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!% temporary allocation of pointers
 nominalEdn => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)
 
 nominalEup => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

!% nominal upstream and downstream head calculation
 where      ( (elem2I(:,e2i_elem_type) == eWeir) .and.  &
              (fEdn(iup) .GE. fEup(idn)) )
            dir = oneI
            nominalEup = fEdn(iup)
            nominalEdn = fEup(idn)
 elsewhere  ( (elem2I(:,e2i_elem_type) == eWeir) .and.   &
              (fEdn(iup) .LT. fEup(idn)) )
            dir = -oneI
            nominalEup = fEup(idn)
            nominalEdn = fEdn(iup)
 endwhere

!% effective head calculation
 where      ( (elem2I(:,e2i_elem_type) == eWeir) .and.  &
              (wEta .GT. wCrest) )

            EffectiveHead = min((wEta-wCrest), (nominalEup - nominalEdn))
 elsewhere  ( (elem2I(:,e2i_elem_type) == eWeir) .and.  &
              (wEta .LE. wCrest) )
            EffectiveHead = zeroR

 elsewhere  ( (elem2I(:,e2i_elem_type) == eWeir) .and.  &
              (wEta .GT. wCrown) )
            EffectiveHead = wCrown - wCrest
 endwhere

 nominalEdn = nullvalueR
 nominalEup = nullvalueR
 nullify(nominalEdn, nominalEup)
 next_e2r_temparray = next_e2r_temparray - 2
!Need a fix for surcharge

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
     wFlow, elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, &
     elemMYN, wWidth, wHeight, wCoeff, wSideSlope, dir,     &
     EffectiveHead, thiscoef)
!
 character(64) :: subroutine_name = 'weir_flow'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,           intent(in out)  :: elem2YN(:,:), elemMYN(:,:)
 real,              intent(in)      :: thiscoef

 
 real,  pointer ::  volume2old(:), volume2new(:), velocity2old(:), velocity2new(:)
 real,  pointer ::  volumeMold(:), volumeMnew(:), velocityMold(:), velocityMnew(:)
 real,  pointer ::  wFlow(:), EffectiveHead(:), depth_element(:)
 real,  pointer ::  wCoeff, wWidth, wHeight, wSideSlope     !Weir discharge coefficient, Width, Height

 integer, pointer :: dir(:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

 depth_element => elem2R(:,e2r_Depth)

 where ( (elem2I(:,e2i_elem_type) == eWeir ).and. &
         (elem2I(:,e2i_geometry)  == eVnotchWeir) )

    wFlow        = dir * wCoeff * wSideSlope * EffectiveHead ** 2.5
    velocity2new = dir * wCoeff * sqrt(abs(EffectiveHead))
    ! Volume is weir flow equation * dt (this case dt = thiscoef)
    volume2new   = thiscoef * wCoeff * wSideSlope * EffectiveHead ** 2.5  
 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_flow
!
!==========================================================================
!==========================================================================
!
subroutine weir_provisional_geometry &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN)
! this subroutine sets the weir geometry to zero.
 character(64) :: subroutine_name = 'weir_provisional_geometry'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,           intent(in out)  :: elem2YN(:,:), elemMYN(:,:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


 where      ( (elem2I(:,e2i_elem_type) == eWeir) )

            elem2R(:,e2r_Area)        = 1.0e-7 
            elem2R(:,e2r_Eta)         = 1.0e-7 
            elem2R(:,e2r_Perimeter)   = 1.0e-7 
            elem2R(:,e2r_HydDepth)    = 1.0e-7 
            elem2R(:,e2r_HydRadius)   = 1.0e-7 
            elem2R(:,e2r_Topwidth)    = 1.0e-7 
            elem2R(:,e2r_Depth)       = 1.0e-7 
 endwhere
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_provisional_geometry
!
!==========================================================================
!==========================================================================
!
 end module weir