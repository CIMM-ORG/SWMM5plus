! This module calculates the flow in orifice elements
!
!========================================================================== 
!
 module orifice


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

    ! public :: weir_step
    ! public :: weir_freesurface_elevation
    ! public :: weir_provisional_geometry

    integer :: debuglevel = 0

 contains
!
!========================================================================== 
!==========================================================================
!
 subroutine orifice_step &
    (e2r_Volume_old, e2r_Velocity_old, eMr_Volume_old, eMr_Velocity_old, &
     e2r_Volume_new, e2r_Velocity_new, eMr_Volume_new, eMr_Velocity_new, &
     elem2R, elemMR, faceI, faceR, faceYN, elem2I, elemMI, elem2YN, &
     elemMYN, thiscoef)
!
 character(64) :: subroutine_name = 'orifice_step'
 
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
 real,  pointer     ::  oFlow(:), oEta(:)
 real,  pointer     ::  oZbottom(:), oCrest(:), oCrown(:), oDischargeCoeff(:)
 real,  pointer     ::  oFullDepth(:), oInletoffset(:)
 real,  pointer     ::  oFullDepth(:), oInletoffset(:)      
 real,  pointer     ::  EffectiveHead(:)
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
 oFlow              => elem2R(:,e2r_Flowrate)
 oEta               => elem2R(:,e2r_eta)  
 oZbottom           => elem2R(:,e2r_Zbottom)
 oHeight            => elem2R(:,e2r_FullDepth)
 oDischargeCoeff    => elem2R(:,e2r_DischargeCoeff1)
 oInletoffset       => elem2R(:,e2r_InletOffset)

 iup          => elem2I(:,e2i_Mface_u)
 idn          => elem2I(:,e2i_Mface_d)

!%  temporary space for elem2
 EffectiveHead => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 oCrest        => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 oCrown       => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 cOrif        => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 cWeir        => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 dir           => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

!%  zero temporary arrays
 EffectiveHead  = zeroR
 cOrif          = zeroR
 cWeir          = zeroR
 dir            = zeroR
 
 !% the orifice flow calculation is only limited to non-sideflow orifice
 !% for now. the code will be updated for sideflow weir soon.
 where      ( (elem2I(:,e2i_elem_type) == eOrifice) )
            oCrest =  oInletoffset + oZbottom
            oCrown =  oCrest + oHeight
 endwhere 


 call orifice_provisional_geometry &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN) 
    
 call orifice_freesurface_elevation &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     fEdn, fEup, iup, idn, dir,oEta)

 ! call weir_effective_head &
 !    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
 !     wCrest, wCrown, wEta, fEup, fEdn, iup, idn, dir, EffectiveHead)
 
 ! call weir_effective_length &
 !    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, & 
 !     wHeight, wSideSlope, wEndContractions, EffectiveCrestLength, &
 !     EffectiveHead, thiscoef)

 ! call weir_flow &
 !    (volume2old, velocity2old, volumeMold, velocityMold, volume2new, &
 !     velocity2new, volumeMnew, velocityMnew, wFlow, elem2R, elemMR,  &
 !     faceR, elem2I, elemMI, elem2YN, elemMYN, wHeight, wSideSlope,   &
 !     wCoeffTriangular, wCoeffRectangular, wEndContractions, dir,     &
 !     EffectiveHead, EffectiveCrestLength, thiscoef)

 ! release temporary arrays
 EffectiveHead          = nullvalueR
 oCrest                 = nullvalueR
 oCrown                 = nullvalueR
 cOrif                  = nullvalueR
 cWeir                  = nullvalueR
 dir                    = nullvalueR

 nullify(EffectiveHead, oCrest, oCrown, cOrif, cWeir, dir)
 next_e2r_temparray = next_e2r_temparray - 6

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine orifice_step
!
!==========================================================================
!==========================================================================
!
subroutine orifice_provisional_geometry &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN)
! this subroutine sets the weir geometry to zero.
 character(64) :: subroutine_name = 'orifice_provisional_geometry'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,           intent(in out)  :: elem2YN(:,:), elemMYN(:,:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


 where      ( (elem2I(:,e2i_elem_type) == eOrifice) )

            elem2R(:,e2r_Area)        = 1.0e-7 
            elem2R(:,e2r_Eta)         = 1.0e-7 
            elem2R(:,e2r_Perimeter)   = 1.0e-7 
            elem2R(:,e2r_HydDepth)    = 1.0e-7 
            elem2R(:,e2r_HydRadius)   = 1.0e-7 
            elem2R(:,e2r_Topwidth)    = 1.0e-7 
            elem2R(:,e2r_Depth)       = 1.0e-7 
 endwhere
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine orifice_provisional_geometry
!
!==========================================================================
!==========================================================================
!
subroutine orifice_freesurface_elevation &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     fEdn, fEup, iup, idn, dir, wEta)
!
 character(64) :: subroutine_name = 'orifice_freesurface_elevation'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,           intent(in out)  :: elem2YN(:,:), elemMYN(:,:)

 real,  pointer   ::  fEdn(:), fEup(:), oEta(:), dir(:)

 integer, pointer ::  iup(:), idn(:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


 where      ( (elem2I(:,e2i_elem_type) == eOrifice) )
            oEta = max(fEdn(iup), fEup(idn))
            dir  = sign(oneR, ( fEdn(iup) - fEup(idn)))
 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine orifice_freesurface_elevation
!
!========================================================================== 
!==========================================================================
!
subroutine orifice_discharge_coefficients &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     fEdn, fEup, iup, idn, dir, wEta)
!
 character(64) :: subroutine_name = 'orifice_discharge_coefficients'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,           intent(in out)  :: elem2YN(:,:), elemMYN(:,:)

 real,  pointer   ::  fEdn(:), fEup(:), oEta(:), dir(:)

 integer, pointer ::  iup(:), idn(:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


 where      ( (elem2I(:,e2i_elem_type) == eOrifice) )
            oEta = max(fEdn(iup), fEup(idn))
            dir  = sign(oneR, ( fEdn(iup) - fEup(idn)))
 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine orifice_discharge_coefficients
!
!==========================================================================
!==========================================================================
!
 end module orifice