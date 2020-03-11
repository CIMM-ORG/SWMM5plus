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
 real,  pointer     ::  wFlow(:), wZbottom(:), wEta(:), wEndContractions(:)
 real,  pointer     ::  wCrest(:), wCrown(:), EffectiveHead(:), EffectiveCrestLength(:) 
 real,  pointer     ::  wCoeff(:), wHeight(:), wSideSlope(:), wInletoffset(:)    
 real,  pointer     ::  fEdn(:), fEup(:), dir(:)

 integer, pointer   ::  iup(:), idn(:)
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

!%  pointers for weir settings
 wflow              => elem2R(:,e2r_Flowrate)  
 wZbottom           => elem2R(:,e2r_Zbottom)
 wHeight            => elem2R(:,e2r_FullDepth)
 wEta               => elem2R(:,e2r_eta)
 wCoeffTriangular   => elem2R(:,e2r_DischargeCoeff1)
 wCoeffRectangular  => elem2R(:,e2r_DischargeCoeff2)
 wInletoffset       => elem2R(:,e2r_InletOffset)
 wSideSlope         => elem2R(:,e2r_LeftSlope)
 wEndContractions   => elem2R(:,e2r_EndContractions)

 iup          => elem2I(:,e2i_Mface_u)
 idn          => elem2I(:,e2i_Mface_d)

!%  temporary space for elem2
 EffectiveHead => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 EffectiveCrestLength => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 wCrest        => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 wCrown       => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 dir           => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

!%  zero temporary arrays
 EffectiveHead = zeroR
 dir = zeroR
 
 where      ( (elem2I(:,e2i_elem_type) == eWeir) )
            wCrest =  wInletoffset + wZbottom
            wCrown =  wCrest + wHeight
 endwhere 


 call weir_provisional_geometry &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN) 
    
 call weir_freesurface_elevation &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     fEdn, fEup, iup, idn, dir,wEta)

 call weir_effective_head &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     wCrest, wCrown, wEta, fEup, fEdn, iup, idn, dir, EffectiveHead)
 
 call weir_effective_length &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, & 
     wHeight, wSideSlope, wEndContractions, EffectiveCrestLength, &
     EffectiveHead, thiscoef)

 call weir_flow &
    (volume2old, velocity2old, volumeMold, velocityMold, volume2new, &
     velocity2new, volumeMnew, velocityMnew, wFlow, elem2R, elemMR,  &
     faceR, elem2I, elemMI, elem2YN, elemMYN, wHeight, wCoeff,       &
     wSideSlope, dir, EffectiveHead, EffectiveCrestLength, thiscoef)

 ! release temporary arrays
 EffectiveHead          = nullvalueR
 EffectiveCrestLength   = nullvalueR
 wCrest                 = nullvalueR
 wCrown                 = nullvalueR
 dir                    = nullvalueI

 nullify(EffectiveHead, EffectiveCrestLength, wCrest, wCrown, dir)
 next_e2r_temparray = next_e2r_temparray - 5

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_step
!
!==========================================================================
!==========================================================================
!
subroutine weir_freesurface_elevation &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     fEdn, fEup, iup, idn, dir, wEta)
!
 character(64) :: subroutine_name = 'weir_freesurface_elevation'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,           intent(in out)  :: elem2YN(:,:), elemMYN(:,:)

 real,  pointer   ::  fEdn(:), fEup(:), wEta(:), dir(:)

 integer, pointer ::  iup(:), idn(:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


 where      ( (elem2I(:,e2i_elem_type) == eWeir) )
            wEta = max(fEdn(iup), fEup(idn))
            dir  = sign(oneR, ( fEdn(iup) - fEup(idn)))
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

 real,  pointer   ::  wCrest(:), wCrown(:), EffectiveHead(:)
 real,  pointer   ::  fEup(:), fEdn(:), wEta(:), dir(:)

 integer, pointer ::  iup(:), idn(:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


!% effective head calculation
 where      ( (elem2I(:,e2i_elem_type) == eWeir) .and. (wEta .GT. wCrest) )
    
            EffectiveHead = min((wEta-wCrest), dir*(fEdn(iup) - fEup(idn)))

 elsewhere  ( (elem2I(:,e2i_elem_type) == eWeir) .and. (wEta .LE. wCrest) )
            EffectiveHead = zeroR
            
 elsewhere  ( (elem2I(:,e2i_elem_type) == eWeir) .and. (wEta .GT. wCrown) )
            EffectiveHead = wCrown - wCrest
 endwhere

!Need a fix for surcharge

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_effective_head
!
!========================================================================== 
!==========================================================================
!
subroutine weir_effective_length &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, wHeight,  &
     wSideSlope, wEndContractions, EffectiveCrestLength, EffectiveHead, &
     thiscoef)
!
 character(64) :: subroutine_name = 'weir_effective_length'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,           intent(in out)  :: elem2YN(:,:), elemMYN(:,:)
 real,  pointer                     :: wLength(:), wHeight(:), wEndContractions(:), wSideSlope(:) 
 real,  pointer                     :: EffectiveHead(:), EffectiveCrestLength(:)
 real,              intent(in)      :: thiscoef

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

 wLength => elem2R(:,e2r_Length)

! Calculate Effective Length (This part is straight up from SWMM source code)
 where ( (elem2I(:,e2i_elem_type) == eWeir ) )
    wLength  = min(twoR*thiscoef*sqrt(grav*wHeight), 200.0)

! effective crest length is used in rectangular and trapezoidal weir flow calculation
 elsewhere ( (elem2I(:,e2i_elem_type) == eWeir ) .and. &
             (elem2I(:,e2i_geometry)  == eRectangular) )
    
    EffectiveCrestLength = max(twoR * wSideSlope * wHeight - 0.1 * &
            wEndContractions * EffectiveHead, 0.0)

 elsewhere ( (elem2I(:,e2i_elem_type) == eWeir ) .and. &
             (elem2I(:,e2i_geometry)  == eTrapezoidal) )

    EffectiveCrestLength = twoR*wSideSlope*wHeight

 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_effective_length
!
!==========================================================================
!==========================================================================
!
 subroutine weir_flow &
    (volume2old, velocity2old, volumeMold, velocityMold, volume2new, &
     velocity2new, volumeMnew, velocityMnew, wFlow, elem2R, elemMR,  &
     faceR, elem2I, elemMI, elem2YN, elemMYN, wHeight, wCoeff,       &
     wSideSlope, dir, EffectiveHead, EffectiveCrestLength, thiscoef)
!
 character(64) :: subroutine_name = 'weir_flow'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,           intent(in out)  :: elem2YN(:,:), elemMYN(:,:)
 real,              intent(in)      :: thiscoef

 
 real,  pointer ::  volume2old(:), volume2new(:), velocity2old(:), velocity2new(:)
 real,  pointer ::  volumeMold(:), volumeMnew(:), velocityMold(:), velocityMnew(:)
 real,  pointer ::  wFlow(:), wSideSlope(:), wCoeff(:), wHeight(:) , dir(:)
 real,  pointer ::  EffectiveHead(:), EffectiveCrestLength(:)   

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


 where ( (elem2I(:,e2i_elem_type) == eWeir ).and. &
         (elem2I(:,e2i_geometry)  == eTriangular) )        
    ! V-notch weir
    wFlow        = dir * wCoeff * wSideSlope * EffectiveHead ** 2.5
    velocity2new = dir * wCoeff * sqrt(abs(EffectiveHead))
    ! Volume is weir flow equation * dt (this case dt = thiscoef)
    volume2new   = thiscoef * wFlow 

 elsewhere ( (elem2I(:,e2i_elem_type) == eWeir ).and. &
           (elem2I(:,e2i_geometry)  == eTrapezoidal) )          
    ! Trapezoidal weir
    wFlow        = dir * wCoeff * wSideSlope * EffectiveHead ** 2.5
    velocity2new = dir * wCoeff * sqrt(abs(EffectiveHead))
    ! Volume is weir flow equation * dt (this case dt = thiscoef)
    volume2new   = thiscoef * wFlow 

 elsewhere ( (elem2I(:,e2i_elem_type) == eWeir ).and. &
          (elem2I(:,e2i_geometry)  == eRectangular) )         
    ! Transverse weir
    wFlow        = dir * wCoeff * EffectiveCrestLength * EffectiveHead ** 1.5
    velocity2new = dir * wCoeff * sqrt(abs(EffectiveHead))
    ! Volume is weir flow equation * dt (this case dt = thiscoef)
    volume2new   = thiscoef * wFlow 
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