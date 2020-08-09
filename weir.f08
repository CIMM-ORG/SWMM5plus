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

    public :: weir_step

    private

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
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
 logical,           intent(in out)  :: faceYN(:,:)
 real,              intent(in)      :: thiscoef

 real,  pointer     ::  volume2old(:), volume2new(:), velocity2old(:), velocity2new(:)
 real,  pointer     ::  volumeMold(:), volumeMnew(:), velocityMold(:), velocityMnew(:)
 real,  pointer     ::  wBreadth(:), wInletoffset(:), cTriangular(:), cRectangular(:)
 real,  pointer     ::  wFullDepth(:), wZbottom(:), wSideSlope(:), wEndContractions(:)
 real,  pointer     ::  wFlow(:), wEta(:), wLength(:), wArea(:), hEffective(:)
 real,  pointer     ::  wPerimeter(:), wHyddepth(:), wHydradius(:), wTopwidth(:)
 real,  pointer     ::  lEffective(:), wCrest(:), wCrown(:), cOrif(:)       
 real,  pointer     ::  subFactor1(:), subFactor2(:)
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

!%  pointers for weir geometry and settings
!%  input
 wBreadth           => elem2R(:,e2r_BreadthScale)
 cTriangular        => elem2R(:,e2r_DischargeCoeff1)
 cRectangular       => elem2R(:,e2r_DischargeCoeff2)
 wInletoffset       => elem2R(:,e2r_InletOffset)
 wFullDepth         => elem2R(:,e2r_FullDepth)
 wSideSlope         => elem2R(:,e2r_SideSlope)
 wEndContractions   => elem2R(:,e2r_EndContractions)
 wZbottom           => elem2R(:,e2r_Zbottom)

!%  output
 wFlow              => elem2R(:,e2r_Flowrate)
 wEta               => elem2R(:,e2r_eta)  
 wLength            => elem2R(:,e2r_Length)
 wArea              => elem2R(:,e2r_Area)
 wPerimeter         => elem2R(:,e2r_Perimeter)
 wHyddepth          => elem2R(:,e2r_HydDepth)
 wHydradius         => elem2R(:,e2r_HydRadius)
 wTopwidth          => elem2R(:,e2r_Topwidth)
 hEffective         => elem2R(:,e2r_Depth)

!%  pointers for upstream and downstream faces
 iup          => elem2I(:,e2i_Mface_u)
 idn          => elem2I(:,e2i_Mface_d)

!%  temporary space for elem2r
 lEffective  => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 wCrest      => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 wCrown      => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)
 
 cOrif       => elem2R(:,e2r_Temp(next_e2r_temparray)) 
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 subFactor1  => elem2R(:,e2r_Temp(next_e2r_temparray)) 
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 subFactor2  => elem2R(:,e2r_Temp(next_e2r_temparray)) 
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

!%  temporary space for elem2I
 dir         => elem2I(:,e2i_Temp(next_e2i_temparray))
 next_e2i_temparray = utility_advance_temp_array (next_e2i_temparray,e2i_n_temp)

!%  temporary space for elem2YN
 maskarrayDnSubmerge  => elem2YN(:,e2YN_Temp(next_e2YN_temparray) )
 next_e2YN_temparray  = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

 maskarrayUpSubmerge  => elem2YN(:,e2YN_Temp(next_e2YN_temparray) )
 next_e2YN_temparray  = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

 maskarraySurcharge   => elem2YN(:,e2YN_Temp(next_e2YN_temparray) )
 next_e2YN_temparray  = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

!%  initializing temporary arrays
 subFactor1        = oneR
 subFactor2        = oneR
 dir               = zeroI

 maskarrayDnSubmerge  = nullvalueL
 maskarrayUpSubmerge  = nullvalueL
 maskarraySurcharge   = nullvalueL

!% set necessary weir setting and find eta on weir element    
call weir_initialize &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,     &
     wInletoffset, wZbottom, wCrown, wCrest, wFullDepth, wLength, &
     wEta, fEdn, fEup, iup, idn, dir, thiscoef)
 
!% calculate the  equivalent orifice discharge coefficient while surcharged
 call weir_surcharge_coefficient &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2new,    &
     velocity2new, volumeMnew, velocityMnew, wFlow, wSideSlope, cTriangular, &
     cRectangular, wBreadth, wArea, dir, wFullDepth, wEndContractions,       &
     lEffective, subFactor1, subFactor2, cOrif, thiscoef)

!% calculate effective head on weir element
 call weir_effective_head &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     wCrest, wCrown, wEta, fEup, fEdn, iup, idn, dir,         &
     hEffective, maskarrayDnSubmerge, maskarrayUpSubmerge,    &
     maskarraySurcharge)

!% calculates weir geometrices for effective head
call weir_geometry &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, wBreadth, &
     wSideSlope, wArea, wPerimeter, wHyddepth, wHydradius, wTopwidth,   &
     hEffective)

!% calculate weir length and effective crest length
 call weir_effective_crest_length &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,        &
     wEndContractions, wFullDepth, hEffective, wBreadth, wSideSlope, &
     lEffective)

!% Villemonte correction for downstream submergence
 call villemonte_weir_submergence_correction &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, wCrest, &
     subFactor1, subFactor2, fEdn, fEup, iup, idn, maskarrayDnSubmerge)

!% Villemonte correction for upstream submergence
 call villemonte_weir_submergence_correction &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, wCrest, &
     subFactor1, subFactor2, fEup, fEdn, idn, iup, maskarrayUpSubmerge)    

!% calculate weir flow
 call weir_flow &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2new, &
     velocity2new, volume2old, velocity2old, wFlow, wArea, wSideSlope,    &
     cTriangular, cRectangular, dir, hEffective, lEffective, subFactor1,  &
     subFactor2, thiscoef)

!%  flow calculataion when flow is surcharged
 call weir_surcharge_flow &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2new, &
     velocity2new, volumeMnew, velocityMnew, wCrest, wCrown, wEta, wFlow, &
     wArea, cOrif, hEffective, fEup, fEdn, iup, idn, dir, thiscoef,       &
     maskarraySurcharge)

 ! release temporary arrays
 lEffective     = nullvalueR
 subFactor1     = nullvalueR
 subFactor2     = nullvalueR
 cOrif          = nullvalueR
 wCrest         = nullvalueR
 wCrown         = nullvalueR
 dir            = nullvalueI

 nullify(lEffective, wCrest, wCrown, cOrif, subFactor1, subFactor2, dir, &
         maskarrayDnSubmerge, maskarrayUpSubmerge, maskarraySurcharge)

 next_e2r_temparray  = next_e2r_temparray  - 6
 next_e2i_temparray  = next_e2i_temparray  - 1
 next_e2YN_temparray = next_e2YN_temparray - 3

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_step
!
!==========================================================================
! PRIVATE BELOW
!==========================================================================
!
subroutine weir_initialize &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,   &
     inletoffset,zbottom, crown, crest, fulldepth, length, eta, &
     faceEtaDn, faceEtaUp, upFace, dnFace, dir, thiscoef)
!
 character(64) :: subroutine_name = 'weir_initialize'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
 real,              intent(in)      :: thiscoef

 real,  pointer   :: inletoffset(:), zbottom(:), crown(:)
 real,  pointer   :: crest(:), fullDepth(:), length(:), eta(:)
 real,  pointer   :: faceEtaDn(:), faceEtaUp(:)

 integer, pointer :: upFace(:), dnFace(:), dir(:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!% find weir crest, crown, length , eta flow direction
 where ( (elem2I(:,e2i_elem_type) == eWeir) )

        crest   =  inletoffset + zbottom
        crown   =  crest + fullDepth
        ! find the weir length
        length  = min(twoR*dt*sqrt(grav*fullDepth), 200.0)
        ! set the free surface elevation at weir element
        eta = max(faceEtaDn(upFace), faceEtaUp(dnFace)) 
        dir  = int(sign(oneR, ( faceEtaDn(upFace) - faceEtaUp(dnFace))))
 endwhere 

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_initialize
!
!==========================================================================
!==========================================================================
!
subroutine weir_surcharge_coefficient &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2new, &
     velocity2new, volumeMnew, velocityMnew, flow, sideslope, cTrig,      &
     cRect, breadth, area, dir, fulldepth, endcontractions, crestlength,  &
     submergenceFactor1, submergenceFactor2, corif, thiscoef)
!
!% when weir is surcharged, the flow becomes orifice flow. this subroutine
!% calculates the equivalent orifice discharge coefficient, cOrif
!
 character(64) :: subroutine_name = 'weir_surcharge_coefficient'
 
 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
 real,              intent(in)      :: thiscoef

 real,    pointer ::  volume2new(:), velocity2new(:), volumeMnew(:), velocityMnew(:)
 real,    pointer ::  flow(:), sideslope(:), cTrig(:), cRect(:), breadth(:), area(:)
 real,    pointer ::  fulldepth(:), endcontractions(:), crestlength(:)
 real,    pointer ::  submergenceFactor1(:), submergenceFactor2(:), corif(:)

 integer, pointer :: dir(:) 

!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!% get effective crest length for maximum weir opening
 call weir_effective_crest_length &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,   &
     endcontractions, fulldepth, fulldepth, breadth, sideslope, &
     crestlength)

!% get flow for maximum weir opening
 call weir_flow &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2new, &
     velocity2new, volumeMnew, velocityMnew, flow, area, sideslope,       &
     cTrig, cRect, dir, fulldepth, crestlength , submergenceFactor1,      &
     submergenceFactor2, thiscoef)
       
 where ( (elem2I(:,e2i_elem_type) == eWeir ) )  
        cOrif = flow / sqrt(fulldepth / twoR)
 endwhere 

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_surcharge_coefficient
!
!==========================================================================
!==========================================================================
!
subroutine weir_effective_head &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, crest,    &
     crown, eta, faceEtaUp, faceEtaDn, upFace, dnFace, dir,             &
     effectivehead, maskarray_dn_submergence, maskarray_up_submergence, &
     maskarray_surcharge)
!
 character(64) :: subroutine_name = 'weir_effective_head'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

 real,  pointer   :: crest(:), crown(:), effectivehead(:)
 real,  pointer   :: faceEtaUp(:), faceEtaDn(:), eta(:)

 integer, pointer :: upFace(:), dnFace(:), dir(:)

 logical, pointer :: maskarray_dn_submergence(:), maskarray_up_submergence(:)
 logical, pointer :: maskarray_surcharge(:)

 real             :: midpt

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!% effective head calculation
 where      ( (elem2I(:,e2i_elem_type) == eWeir) .and. (eta .LE. crest) )

            effectivehead = zeroR
 elsewhere  ( (elem2I(:,e2i_elem_type) == eWeir) .and. (eta .GT. crest) .and.&
              (eta .LT. crown) )

            effectivehead = eta - crest
            ! downstream submergence 
            maskarray_dn_submergence = ((dir .GT. zeroI) .and. &
                                        (faceEtaUp(dnFace) .GT. crest))
            ! upstream submergance
            maskarray_up_submergence = ((dir .LT. zeroI) .and. &
                                        (faceEtaDn(upFace) .GT. crest))
            
 elsewhere  ( (elem2I(:,e2i_elem_type) == eWeir) .and. (eta .GT. crown) )
            ! non surcharge weir flow
            effectivehead = crown - crest               
            ! surcharge conditon
            maskarray_surcharge = ((elem2I(:,e2i_elem_type) == eWeir) .and. &
                                   (elem2YN(:,e2YN_CanSurcharge))     .and. &
                                   (eta .GT. crown) )
 endwhere

!% surcharged weir flow is not added here. for surcharged weir, the effective head ...
!% calculation is commented out below. 

!  do mm = 1, N_elem2
!     if ( (elem2I(mm,e2i_elem_type) == eWeir) .and. (eta(mm) .GT. crest(mm)) ) then

!             effectivehead(mm) = (eta(mm)-crest(mm))

!     elseif ( (elem2I(mm,e2i_elem_type) == eWeir) .and. (eta(mm) .LE. crest(mm)) ) then
!             effectiveheadive(mm) = zeroR
            
!     elseif ( (elem2I(mm,e2i_elem_type) == eWeir) .and. (eta(mm) .GT. crown(mm)) ) then

!             if (elem2YN(mm,e2YN_CanSurcharge)) then
!                 ! weir surcharge condition

!                 midpt = (crown(mm) + crest(mm)) / twoR

!                 effectivehead(mm) = min(( eta(mm) - midpt ), &
!                                           dir(mm)*(faceEtaDn(upFace(mm))) - faceEtaUp(dnFace(mm))))
!             else

!                 effectivehead(mm) = crown(mm) - crest(mm)

!             endif
!     endif
! end do

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_effective_head
!
!==========================================================================
!==========================================================================
!
subroutine weir_geometry &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, breadth, &
     slope, area, perimeter, hyddepth, hydradius, topwidth, depth)
!
 character(64) :: subroutine_name = 'weir_geometry'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

 real,  pointer     ::  breadth(:), slope(:), area(:)
 real,  pointer     ::  perimeter(:), hyddepth(:), hydradius(:)
 real,  pointer     ::  topwidth(:), depth(:) 

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

 where      ( (elem2I(:,e2i_elem_type) == eWeir) .and. &
              (elem2I(:,e2i_geometry) == eRectangular ) )

        area        =   depth * breadth
        topwidth    =   breadth
        hyddepth    =   depth
        perimeter   =   breadth + twoR * hyddepth
        hydradius   =   area / perimeter

 elsewhere  ( (elem2I(:,e2i_elem_type) == eWeir) .and. &
              (elem2I(:,e2i_geometry) == eTrapezoidal ) )

        area        =   (breadth + slope * depth) * depth
        topwidth    =   breadth + twoR * slope * depth
        hyddepth    =   area / topwidth
        perimeter   =   breadth + twoR * depth * sqrt(oneR + slope ** twoR)
        hydradius   =   area / perimeter

 elsewhere  ( (elem2I(:,e2i_elem_type) == eWeir) .and. &
              (elem2I(:,e2i_geometry) == eTriangular ) ) 

        area        =   slope * depth ** twoR 
        topwidth    =   twoR * slope * depth 
        hyddepth    =   area / topwidth
        perimeter   =   twoR * depth * sqrt(oneR + slope ** twoR)
        hydradius   =   area / perimeter
 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_geometry
!
!==========================================================================
!==========================================================================
!
subroutine weir_effective_crest_length &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     endcontractions, fulldepth, depth, breadth, sideslope,   &
     crestlength)
!%  ths subroutine calculates effective creast length for 
!%  trapezoidal and rectangular weir
 character(64) :: subroutine_name = 'weir_effective_crest_length'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

 real,  pointer                     :: endcontractions(:), fulldepth(:)
 real,  pointer                     :: depth(:), breadth(:), sideslope(:)
 real,  pointer                     :: crestlength(:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

! effective crest length is used in rectangular and trapezoidal weir flow calculation
 where ( (elem2I(:,e2i_elem_type) == eWeir ) .and. &
         (elem2I(:,e2i_geometry)  == eRectangular) )
    
    crestlength = max(breadth - 0.1 * endcontractions * depth, zeroR)

 elsewhere ( (elem2I(:,e2i_elem_type) == eWeir ) .and. &
             (elem2I(:,e2i_geometry)  == eTrapezoidal) )

    crestlength = breadth + twoR * sideslope * fulldepth

 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_effective_crest_length
!
!==========================================================================
!==========================================================================
!
 subroutine villemonte_weir_submergence_correction &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, crest, &
     submergenceFactor1, submergenceFactor2, faceEtaDn, faceEtaUp,   &
     upFace, dnFace, maskarray_submergence)
!
 character(64) :: subroutine_name = 'villemonte_weir_submergence_correction'

 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

 real,  pointer ::  crest(:), submergenceFactor1(:), submergenceFactor2(:)       
 real,  pointer ::  faceEtaDn(:), faceEtaUp(:)

 integer, pointer :: upFace(:), dnFace(:)

 logical, pointer :: maskarray_submergence(:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!% calculate the submergance factor for different weirs according to Villemonte 1974
 where ( (elem2I(:,e2i_weir_elem_type) == eVnotchWeir) .and. &
         (maskarray_submergence) )     
    ! V-notch weir
    submergenceFactor1 = (oneR - ((faceEtaUp(dnFace) - crest) / (faceEtaDn(upFace) - crest)) &
        ** 2.5) ** 0.385

 elsewhere ( (elem2I(:,e2i_weir_elem_type) == eTrapezoidalWeir) .and. &
             (maskarray_submergence) )           
    ! Trapezoidal weir
    submergenceFactor1 = (oneR - ((faceEtaUp(dnFace) - crest) / (faceEtaDn(upFace) - crest)) &
        ** 2.5) ** 0.385
    submergenceFactor2 = (oneR - ((faceEtaUp(dnFace) - crest) / (faceEtaDn(upFace) - crest)) &
        ** 1.5) ** 0.385

 elsewhere ( (elem2I(:,e2i_weir_elem_type) == eTransverseWeir) .and. &
             (maskarray_submergence) )         
    ! Transverse weir
    submergenceFactor1 = (oneR - ((faceEtaUp(dnFace) - crest) / (faceEtaDn(upFace) - crest)) &
        ** 1.5) ** 0.385

 elsewhere ( (elem2I(:,e2i_weir_elem_type) == eSideFlowWeir) .and. &
             (maskarray_submergence) )     
   ! Side flow weir
   submergenceFactor1 = (oneR - ((faceEtaUp(dnFace) - crest) / (faceEtaDn(upFace) - crest)) &
        ** 1.67) ** 0.385

 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine villemonte_weir_submergence_correction
!
!========================================================================== 
!==========================================================================
!
 subroutine weir_flow &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2new, &
     velocity2new, volumeMnew, velocityMnew, flow, area, sideslope,       &
     cTrig, cRect, dir, effectivehead, crestlength , submergenceFactor1,  &
     submergenceFactor2, thiscoef)
!
 character(64) :: subroutine_name = 'weir_flow'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
 real,              intent(in)      :: thiscoef

 
 real,  pointer   ::  volume2new(:), velocity2new(:), volumeMnew(:), velocityMnew(:)
 real,  pointer   ::  flow(:), area(:), sideslope(:), cTrig(:), cRect(:) 
 real,  pointer   ::  submergenceFactor1(:), submergenceFactor2(:)
 real,  pointer   ::  effectivehead(:), crestlength(:)

 integer, pointer ::  dir(:)    

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

 where ( elem2I(:,e2i_weir_elem_type) == eVnotchWeir )     
    ! V-notch weir
    ! fs = submergance_factor
    ! d = direction
    ! Q = d*fs*(Cw1SH^2.5)  
    flow        = dir * submergenceFactor1 * cTrig * sideslope * &
        effectivehead ** 2.5

    velocity2new = flow / area

    ! Volume = Q * dt
    volume2new   = dt * flow 

 elsewhere ( elem2I(:,e2i_weir_elem_type) == eTrapezoidalWeir )          
    ! Trapezoidal weir
    ! Q = d*fs1*(Cw1SH^2.5) + d*fs2*(Cw2LH^1.5)
    flow        = dir * submergenceFactor1 * (cTrig * sideslope *    &
        effectivehead ** 2.5) + dir * submergenceFactor2 * (cRect *  &
        crestlength * effectivehead ** 1.5)

    velocity2new = flow / area

    ! Volume = Q * dt
    volume2new   = dt * flow 

 elsewhere ( elem2I(:,e2i_weir_elem_type) == eTransverseWeir )        
    ! Transverse weir
    ! Q = d*fs2*(Cw2LH^1.5)
    flow        = dir * submergenceFactor1 * cRect * &
        crestlength * effectivehead ** 1.5

    velocity2new = flow / area

    ! Volume = Q * dt
    volume2new   = dt * flow 

 elsewhere ( (elem2I(:,e2i_weir_elem_type) == eSideFlowWeir )  .and. &
             (dir .LE. zeroR) )    
    ! Side flow weir for reverse flow behaves like a Transverse weir
    ! Q = d*fs2*(Cw2LH^1.5)
    flow        = dir * submergenceFactor1 * cRect * &
        crestlength * effectivehead ** 1.5

    velocity2new = flow / area

    ! Volume = Q * dt
    volume2new   = dt * flow 

 elsewhere ( (elem2I(:,e2i_weir_elem_type) == eSideFlowWeir )  .and. &
             (dir .GT. zeroR) )    
    ! Corrected formula (see Metcalf & Eddy, Inc., Wastewater Engineering, McGraw-Hill, 1972 p. 164).
    flow        = dir * submergenceFactor1 * cRect * &
        (crestlength ** 0.83) * (effectivehead ** 1.67)

    velocity2new = flow / area

    ! Volume = Q * dt
    volume2new   = dt * flow 
 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_flow
!
!==========================================================================
!==========================================================================
!
 subroutine weir_surcharge_flow &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2new, &
     velocity2new, volume2old, velocity2old, crest, crown, eta, flow,     &
     area, cOrif, effectivehead, faceEtaup, faceEtadn, upFace, dnFace,    &
     dir, thiscoef, maskarray_surcharge)
!
 character(64) :: subroutine_name = 'weir_surcharge_flow'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
 real,              intent(in)      :: thiscoef

 
 real,    pointer ::  volume2new(:), velocity2new(:), volume2old(:), velocity2old(:)
 real,    pointer ::  crest(:), crown(:), eta(:), flow(:), area(:), cOrif(:)
 real,    pointer ::  effectivehead(:), faceEtaup(:), faceEtadn(:)

 integer, pointer ::  upFace(:), dnFace(:), dir(:)
 logical, pointer ::  maskarray_surcharge(:)

!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!% surchaged flow calculation
 where ( maskarray_surcharge ) 
    ! for surcharged flow, head is calculated from the midpoint of the weir opening    
    effectivehead = min(( eta - ((crown + crest) / twoR)), &
                          dir*(faceEtadn(upFace) - faceEtaup(dnFace)) )

    flow         = dir * cOrif * sqrt(abs(effectivehead))
    ! blend new velocity with old velocity -- needs further checking
    velocity2new  = flow / area
    ! Volume is weir flow equation * dt
    ! blend new volume with old volume -- needs further checking
    volume2new    =  dt * flow

 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_surcharge_flow
!
!==========================================================================
! END OF MODULE weir
!==========================================================================
!
 end module weir