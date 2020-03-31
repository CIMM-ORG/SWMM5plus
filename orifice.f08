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
 logical,           intent(in)      :: elem2YN(:,:), elemMYN(:,:)
 logical,           intent(in out)  :: faceYN(:,:)
 real,              intent(in)      :: thiscoef

 real,  pointer     ::  volume2old(:), volume2new(:), velocity2old(:), velocity2new(:)
 real,  pointer     ::  volumeMold(:), volumeMnew(:), velocityMold(:), velocityMnew(:)
 real,  pointer     ::  oFlow(:), oEta(:), oZbottom(:), oWidth(:)
 real,  pointer     ::  oFullDepth(:), oInletoffset(:) , oDischargeCoeff(:)  
 real,  pointer     ::  hCrest(:), hCrown(:), hCrit(:)   
 real,  pointer     ::  hEffective(:), subFactor(:)
 real,  pointer     ::  fEdn(:), fEup(:)

 integer, pointer   ::  iup(:), idn(:), dir(:)

 logical, pointer   ::  maskarrayUpSubmerge(:), maskarrayDnSubmerge(:)
 logical, pointer   ::  maskarraySurcharge(:)

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
 oFullDepth         => elem2R(:,e2r_FullDepth)
 oWidth             => elem2R(:,e2r_BreadthScale)
 oDischargeCoeff    => elem2R(:,e2r_DischargeCoeff1)
 oInletoffset       => elem2R(:,e2r_InletOffset)

 iup          => elem2I(:,e2i_Mface_u)
 idn          => elem2I(:,e2i_Mface_d)

!%  temporary space for elem2
 hEffective    => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 hCrest        => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 hCrown        => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 hCrit         => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 cOrif         => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 cWeir         => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 subFactor     => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 dir           => elem2I(:,e2i_Temp(next_e2i_temparray))
 next_e2i_temparray = utility_advance_temp_array (next_e2i_temparray,e2i_n_temp)

 maskarrayDnSubmerge  => elem2YN(:,e2YN_Temp(next_e2YN_temparray) )
 next_e2YN_temparray  = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

 maskarrayUpSubmerge  => elem2YN(:,e2YN_Temp(next_e2YN_temparray) )
 next_e2YN_temparray  = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

 maskarraySurcharge   => elem2YN(:,e2YN_Temp(next_e2YN_temparray) )
 next_e2YN_temparray  = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

!%  zero temporary arrays
 hEffective     = nullvalueR
 hCrown         = nullvalueR
 hCrest         = nullvalueR
 hCrit          = nullvalueR
 cOrif          = nullvalueR
 cWeir          = nullvalueR
 subFactor      = nullvalueR

 dir            = nullvalueI

 maskarrayDnSubmerge  = nullvalueL
 maskarrayUpSubmerge  = nullvalueL
 maskarraySurcharge   = nullvalueL

! !% sets all orifice element geometry values ~ zero values 
!  call orifice_provisional_geometry &
!     (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN)

!% sets necessary orifice setting and find eta on orifice element
 call orifice_initialize &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,    &
     oInletoffset, oZbottom, oFullDepth, oLength, oEta, hCrown,  &
     hCrest, fEdn, fEup, iup, idn, dir, thiscoef) 

!% calculates the  equivalent orificeand weir discharge coefficients
 call orifice_equivalent_discharge_coefficient &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, oWidth, &
     oFullDepth, oDischargeCoeff, hCrit, cOrif, cWeir)

!% calculates effective head in orifice elements
 call orifice_effective_head &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     oEta, fEup, fEdn, iup, idn, dir, hCrest, hCrown, hcrit,  &
     hEffective, subFactor)

!% calculates flow in orifice elements
 call orifice_flow &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2new, &
     velocity2new, volumeMnew, velocityMnew, oFlow, cOrif, cWeir, oWidth, &
     oFullDepth, hEffective, dir, subFactor, thiscoef)

 ! release temporary arrays
 hEffective     = nullvalueR
 hCrown         = nullvalueR
 hCrest         = nullvalueR
 hCrit          = nullvalueR
 cOrif          = nullvalueR
 cWeir          = nullvalueR
 subFactor      = nullvalueR

 dir            = nullvalueI

 nullify(EffectiveHead, hCrest, hCrown, hCrit, cOrif, cWeir, subFactor, dir)
 next_e2r_temparray  = next_e2r_temparray  - 7
 next_e2i_temparray  = next_e2i_temparray  - 1
 next_e2YN_temparray = next_e2YN_temparray - 3

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
 logical,           intent(in)      :: elem2YN(:,:), elemMYN(:,:)

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
subroutine orifice_initialize &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,    &
     oInletoffset, oZbottom, oFullDepth, oLength, oEta, hCrown,  &
     hCrest, fEdn, fEup, iup, idn, dir, thiscoef)
!
 character(64) :: subroutine_name = 'orifice_initialize'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
 real,              intent(in)      :: thiscoef

 real,  pointer   :: oInletoffset(:), oZbottom(:)
 real,  pointer   :: oFullDepth(:), oLength(:), oEta(:)
 real,  pointer   :: hCrest(:), hCrown(:)
 real,  pointer   :: fEdn(:), fEup(:)

 integer, pointer :: iup(:), idn(:), dir(:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!% find orifice crest, crown, length , eta flow direction
 where ( (elem2I(:,e2i_elem_type) == eOrifice) )

        hCrest   = oInletoffset + oZbottom
        hCrown   = oCrest + oFullDepth
        ! find the effective weir length
        oLength  = min(twoR*thiscoef*sqrt(grav*oFullDepth), 200.0)
        ! set the free surface elevation at weir element
        oEta = max(fEdn(iup), fEup(idn)) 
        dir  = int(sign(oneR, ( fEdn(iup) - fEup(idn))))
 endwhere 

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine orifice_initialize
!
!==========================================================================
!==========================================================================
!
subroutine orifice_equivalent_discharge_coefficient &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, oWidth, &
     oFullDepth, oDischargeCoeff, hCrit, cOrif, cWeir)
!
 character(64) :: subroutine_name = 'orifice_equivalent_discharge_coefficient'
 
 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
 real,              intent(in)      :: thiscoef

 real,    pointer ::  oWidth(:), oFullDepth(:), oDischargeCoeff(:)
 real,    pointer ::  hCrit(:), cOrif(:), cWeir(:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name
 
 !% find effective orifice discharge coefficient.
 !% Co = CdAo(g)^0.5
 !% Cd = discharge coefficient
 !% Ao = Area of orifice opening
 where ( (elem2I(:, e2i_elem_type) == eOrifice) .and. 
         (elem2I(:e2i_geometry) == eCircular) )

        cOrif = oDischargeCoeff * (pi/4.0 *oFullDepth **2) * &
        sqrt(twoR * grav)

 elsewhere ( (elem2I(:, e2i_elem_type) == eOrifice) .and. 
         (elem2I(:e2i_geometry) == eRectangular) )

        cOrif = oDischargeCoeff * (oFullDepth * oWidth) * &
        sqrt(twoR * grav)
 endwhere   

 !% find critical height above opening where orifice flow
 !% turns into weir flow. It equals (Co/Cw)*(Area/Length)
 !% where Co is the orifice coeff., Cw is the weir coeff/sqrt(2g),
 !% Area is the area of the opening, and Length = circumference
 !% of the opening. For a basic sharp crested weir, Cw = 0.414.

 where ( (elem2I(:,e2i_orif_elem_type) == eBottomOrifice) .and. &
         (elem2I(:e2i_geometry) == eRectangular) )

        hCrit = oDischargeCoeff * (oFullDepth * oWidth) / &
        (0.414 * twoR * (oFullDepth + oWidth) )
        cWeir = oDischargeCoeff * (oFullDepth * oWidth) * &
        sqrt(twoR * grav * hCrit)

 elsewhere ( (elem2I(:,e2i_orif_elem_type) == eBottomOrifice) .and. &
             (elem2I(:e2i_geometry) == eCircular) )
        hcrit = oDischargeCoeff * oFullDepth / (0.414 * 4.0)
        cWeir = oDischargeCoeff * (pi/4.0 *oFullDepth **2) * &
        sqrt(twoR * grav * hCrit)

 elsewhere ( (elem2I(:,e2i_orif_elem_type) == eSideOrifice) .and. &
             (elem2I(:e2i_geometry) == eRectangular))

        hCrit = oFullDepth
        cWeir = oDischargeCoeff * (oFullDepth * oWidth) * &
        sqrt(grav * hCrit)

 elsewhere ( (elem2I(:,e2i_orif_elem_type) == eSideOrifice) .and. &
             (elem2I(:e2i_geometry) == eCircular))

        hCrit = oFullDepth
        cWeir = oDischargeCoeff * (pi/4.0 *oFullDepth **2) * &
        sqrt(grav * hCrit)

 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine orifice_equivalent_discharge_coefficient
!
!==========================================================================
!==========================================================================
!
subroutine orifice_effective_head &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     oEta, fEup, fEdn, iup, idn, dir, hCrest, hCrown, hcrit,  &
     hEffective, subFactor)
!
 character(64) :: subroutine_name = 'orifice_effective_head'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

 real,  pointer   :: hCrest(:), hCrown(:), hEffective(:), hCrit(:)
 real,  pointer   :: fEup(:), fEdn(:), oEta(:), subFactor(:)

 integer, pointer :: iup(:), idn(:), dir(:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!% effective head calculation for bottom orifice
 where      ( (elem2I(:,e2i_orif_elem_type) == eBottomOrifice) .and. &
              (oEta .LE. hCrest) )

            hEffective = zeroR

 elsewhere  ( (elem2I(:,e2i_orif_elem_type) == eBottomOrifice) .and. &
              (oEta .GT. hCrest) )

            hEffective = min( (oEta - hCrest), dir * ( fEdn(iup) - fEup(idn)) )
            ! find fraction of critical height for which weir flow occurs
            subFactor  = min(hEffective/hCrit, 1.0)
 endwhere
 
 !% find degree of submergence for side orifice
 where      ( (elem2I(:,e2i_orif_elem_type) == eSideOrifice) .and. &
              (oEta .LT. hCrown) .and. (hCrown .GT. hCrest))

            subFactor = (oEta - hCrown) / (hCrown - hCrest)

 elsewhere ( elem2I(:,e2i_orif_elem_type) == eSideOrifice ) 
            subFactor = OneR
 endwhere

!% effective head calculation for side orifice
 where      ( (elem2I(:,e2i_orif_elem_type) == eSideOrifice) .and. &
              (subFactor .LE. zeroR) )

            hEffective = zeroR

 elsewhere ( (elem2I(:,e2i_orif_elem_type) == eSideOrifice) .and. &
             (subFactor .GT. zeroR) .and. (subFactor .LT. oneR)) 

            hEffective = oEta - hCrest

 elsewhere (elem2I(:,e2i_orif_elem_type) == eSideOrifice)

            hEffective = min((oEta - (oCrest + oCrown)/twoR), &
                dir * ( fEdn(iup) - fEup(idn)))
 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine orifice_effective_head
!
!==========================================================================
!==========================================================================
!
 subroutine orifice_flow &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2new, &
     velocity2new, volumeMnew, velocityMnew, oFlow, cOrif, cWeir, oWidth, &
     oFullDepth, hEffective, dir, subFactor, thiscoef)
!
 character(64) :: subroutine_name = 'orifice_flow'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
 real,              intent(in)      :: thiscoef

 
 real,  pointer   ::  volume2new(:), velocity2new(:), volumeMnew(:), velocityMnew(:)
 real,  pointer   ::  oFlow(:), cOrif(:), cWeir(:), oWidth(:), oFullDepth(:) 
 real,  pointer   ::  subFactor(:), hEffective(:)

 integer, pointer :: dir(:)    

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!% flow and volume calculation 
 where     ( (elem2I(:,e2i_elem_type) == eOrifice) .and. &
             (subFactor .LE. zeroR) )     
    
        oFlow        = zeroR
        velocity2new = zeroR
        volume2new   = zeroR

 elsewhere ( (elem2I(:,e2i_elem_type) == eOrifice) .and. &
             (subFactor .GT. zeroR) .and. (subFactor .LT. oneR) ) 

        oFlow        = dir * cWeir * subFactor ** 1.5
        !% do the math to find the velocity (LATER) %!
        velocity2new = cWeir * subFactor ** 1.5
        volume2new   = oFlow * thiscoef

 elsewhere (elem2I(:,e2i_elem_type) == eOrifice) 

        oFlow        = dir * cOrif * sqrt(abs(hEffective)
        !% do the math to find the velocity (LATER) %!
        ! divide the flow by orifice opening (LATER)
        velocity2new = dir * cOrif * sqrt(abs(hEffective)
        volume2new   = oFlow * thiscoef
 endwhere
   
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine orifice_flow
!
!==========================================================================
!==========================================================================
!
 end module orifice