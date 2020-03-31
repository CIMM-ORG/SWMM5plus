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
    public :: weir_provisional_geometry

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
 real,  pointer     ::  wFlow(:), wEta(:), wZbottom(:), wCrest(:), wCrown(:), wLength(:)
 real,  pointer     ::  wEndContractions(:), wCoeffTriangular(:), wCoeffRectangular(:)
 real,  pointer     ::  wHeight(:), wSideSlope(:), wInletoffset(:), wCoeffOrif(:)       
 real,  pointer     ::  EffectiveHead(:), EffectiveCrestLength(:), subFactor1(:)
 real,  pointer     ::  subFactor2(:), fEdn(:), fEup(:)

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
 wFlow              => elem2R(:,e2r_Flowrate)
 wEta               => elem2R(:,e2r_eta)  
 wZbottom           => elem2R(:,e2r_Zbottom)
 wLength            => elem2R(:,e2r_Length)
 wHeight            => elem2R(:,e2r_FullDepth)
 wCoeffTriangular   => elem2R(:,e2r_DischargeCoeff1)
 wCoeffRectangular  => elem2R(:,e2r_DischargeCoeff2)
 wInletoffset       => elem2R(:,e2r_InletOffset)
 wSideSlope         => elem2R(:,e2r_SideSlope)
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
 
 wCoeffOrif   => elem2R(:,e2r_Temp(next_e2r_temparray)) 
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 subFactor1   => elem2R(:,e2r_Temp(next_e2r_temparray)) 
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 subFactor2   => elem2R(:,e2r_Temp(next_e2r_temparray)) 
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
 EffectiveHead        = zeroR
 EffectiveCrestLength = zeroR

 subFactor1           = oneR
 subFactor2           = oneR

 maskarrayDnSubmerge  = nullvalueL
 maskarrayUpSubmerge  = nullvalueL
 maskarraySurcharge   = nullvalueL

 dir = zeroI

! !% set all weir element geometry values ~ zero values
!  call weir_provisional_geometry &
!     (elem2R, elemMR, faceR, elem2I, elemMI)

!% set necessary weir setting and find eta on weir element    
call weir_initialize &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,  &
     wInletoffset, wZbottom, wCrown, wCrest, wHeight, wLength, &
     wEta, fEdn, fEup, iup, idn, dir, thiscoef)
 
!% calculate the  equivalent orifice discharge coefficient while surcharged
 call weir_surcharge_coefficient &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2new, &
     velocity2new, volumeMnew, velocityMnew, wFlow, wSideSlope,           &
     wCoeffTriangular, wCoeffRectangular, dir, wHeight, wEndContractions, &
     EffectiveCrestLength, subFactor1, subFactor2, wCoeffOrif, thiscoef)

!% calculate effective head on weir element
 call weir_effective_head &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     wCrest, wCrown, wEta, fEup, fEdn, iup, idn, dir,         &
     EffectiveHead, maskarrayDnSubmerge, maskarrayUpSubmerge, &
     maskarraySurcharge)

!% calculate weir length and effective crest length
 call weir_effective_crest_length &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, & 
     wHeight, wSideSlope, wEndContractions, EffectiveCrestLength, &
     EffectiveHead, thiscoef)

!% Villemonte correction for downstream submergence
 call villemonte_submergence_correction &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, wCrest, &
     subFactor1, subFactor2, fEdn, fEup, iup, idn, maskarrayDnSubmerge)

!% Villemonte correction for upstream submergence
 call villemonte_submergence_correction &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, wCrest, &
     subFactor1, subFactor2, fEup, fEdn, idn, iup, maskarrayDnSubmerge)    

!% calculate weir flow
 call weir_flow &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2new, &
     velocity2new, volumeMnew, velocityMnew, wFlow, wSideSlope,           &
     wCoeffTriangular, wCoeffRectangular, dir, EffectiveHead,             &
     EffectiveCrestLength, subFactor1, subFactor2, thiscoef)

!%  flow calculataion when flow is surcharged
 call weir_surcharge_flow &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2new, &
     velocity2new, volumeMnew, velocityMnew, fEup, fEdn, iup, idn,        &
     wCrest, wCrown, wEta, wFlow, wCoeffOrif, dir,  EffectiveHead,        &
     thiscoef, maskarraySurcharge)

 ! release temporary arrays
 EffectiveHead          = nullvalueR
 EffectiveCrestLength   = nullvalueR
 subFactor1             = nullvalueR
 subFactor2             = nullvalueR
 wCoeffOrif             = nullvalueR
 wCrest                 = nullvalueR
 wCrown                 = nullvalueR

 dir                    = nullvalueI

 nullify(EffectiveHead, EffectiveCrestLength, wCoeffOrif, wCrest, &
    wCrown, subFactor1, subFactor2, dir, maskarrayDnSubmerge,     &
    maskarrayUpSubmerge, maskarraySurcharge)

 next_e2r_temparray  = next_e2r_temparray  - 7
 next_e2i_temparray  = next_e2i_temparray  - 1
 next_e2YN_temparray = next_e2YN_temparray - 3

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_step
!
!==========================================================================
!==========================================================================
!
subroutine weir_initialize &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,  &
     wInletoffset, wZbottom, wCrown, wCrest, wHeight, wLength, &
     wEta, fEdn, fEup, iup, idn, dir, thiscoef)
!
 character(64) :: subroutine_name = 'weir_initialize'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
 real,              intent(in)      :: thiscoef

 real,  pointer   :: wInletoffset(:), wZbottom(:), wCrown(:)
 real,  pointer   :: wCrest(:), wHeight(:), wLength(:),wEta(:)
 real,  pointer   :: fEdn(:), fEup(:)

 integer, pointer :: iup(:), idn(:), dir(:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!% find weir crest, crown, length , eta flow direction
 where ( (elem2I(:,e2i_elem_type) == eWeir) )

        wCrest   =  wInletoffset + wZbottom
        wCrown   =  wCrest + wHeight
        ! find the effective weir length
        wLength  = min(twoR*thiscoef*sqrt(grav*wHeight), 200.0)
        ! set the free surface elevation at weir element
        wEta = max(fEdn(iup), fEup(idn)) 
        dir  = int(sign(oneR, ( fEdn(iup) - fEup(idn))))
 endwhere 

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_initialize
!
!==========================================================================
!==========================================================================
!
subroutine weir_surcharge_coefficient &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2new, &
     velocity2new, volumeMnew, velocityMnew, wFlow, wSideSlope,           &
     wCoeffTriangular, wCoeffRectangular, dir, wHeight, wEndContractions, &
     EffectiveCrestLength, subFactor1, subFactor2, wCoeffOrif, thiscoef)
!
!% when weir is surcharged, the flow becomes orifice flow. this subroutine
!% calculates the equivalent orifice discharge coefficient, wCoeffOrif
!
 character(64) :: subroutine_name = 'weir_surcharge_coefficient'
 
 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
 real,              intent(in)      :: thiscoef

 real,    pointer ::  volume2new(:), velocity2new(:), volumeMnew(:), velocityMnew(:)
 real,    pointer ::  wFlow(:), wSideSlope(:), wCoeffTriangular(:), wCoeffRectangular(:)
 real,    pointer ::  wEndContractions(:), wCoeffOrif(:), wHeight(:), EffectiveCrestLength(:)
 real,    pointer ::  subFactor1(:), subFactor2(:)

 integer, pointer :: dir(:) 

!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!% get effective crest length for maximum weir opening
 call weir_effective_crest_length &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, & 
     wHeight, wSideSlope, wEndContractions, EffectiveCrestLength, &
     wHeight, thiscoef)

!% get flow for maximum weir opening
 call weir_flow &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2new, &
     velocity2new, volumeMnew, velocityMnew, wFlow, wSideSlope,           &
     wCoeffTriangular, wCoeffRectangular, dir, wHeight,                   &
     EffectiveCrestLength, subFactor1, subFactor2, thiscoef)

        
 where ( (elem2I(:,e2i_elem_type) == eWeir ) )  
        wCoeffOrif = wFlow / sqrt(wHeight / twoR)
 endwhere 

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_surcharge_coefficient
!
!==========================================================================
!==========================================================================
!
subroutine weir_effective_head &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     wCrest, wCrown, wEta, fEup, fEdn, iup, idn, dir,         &
     EffectiveHead, maskarray1, maskarray2, maskarray3)
!
 character(64) :: subroutine_name = 'weir_effective_head'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

 real,  pointer   :: wCrest(:), wCrown(:), EffectiveHead(:)
 real,  pointer   :: fEup(:), fEdn(:), wEta(:)

 integer, pointer :: iup(:), idn(:), dir(:)
 logical, pointer :: maskarray1(:), maskarray2(:), maskarray3(:)

 real             :: wMidPt

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!% effective head calculation
 where      ( (elem2I(:,e2i_elem_type) == eWeir) .and. (wEta .LE. wCrest) )

            EffectiveHead = zeroR
 elsewhere  ( (elem2I(:,e2i_elem_type) == eWeir) .and. (wEta .GT. wCrest) .and.&
              (wEta .LT. wCrown) )

            EffectiveHead = wEta - wCrest
            ! downstream submergence 
            maskarray1 = ((dir .GT. zeroI) .and. (fEup(idn) .GT. wCrest))
            ! upstream submergance
            maskarray2 = ((dir .LT. zeroI) .and. (fEdn(iup) .GT. wCrest))
            
 elsewhere  ( (elem2I(:,e2i_elem_type) == eWeir) .and. (wEta .GT. wCrown) )
            ! non surcharge weir flow
            EffectiveHead = wCrown - wCrest               
            ! surcharge conditon
            maskarray3 = ((elem2I(:,e2i_elem_type) == eWeir) .and. &
                         (elem2YN(:,e2YN_CanSurcharge))     .and. &
                         (wEta .GT. wCrown) )
 endwhere

!% surcharged weir flow is not added here. for surcharged weir, the effective head ...
!% calculation is commented out below. 

!  do mm = 1, N_elem2
!     if ( (elem2I(mm,e2i_elem_type) == eWeir) .and. (wEta(mm) .GT. wCrest(mm)) ) then

!             EffectiveHead(mm) = (wEta(mm)-wCrest(mm))

!     elseif ( (elem2I(mm,e2i_elem_type) == eWeir) .and. (wEta(mm) .LE. wCrest(mm)) ) then
!             EffectiveHead(mm) = zeroR
            
!     elseif ( (elem2I(mm,e2i_elem_type) == eWeir) .and. (wEta(mm) .GT. wCrown(mm)) ) then

!             if (elem2YN(mm,e2YN_CanSurcharge)) then
!                 ! weir surcharge condition

!                 wMidPt = (wCrown(mm) + wCrest(mm)) / twoR

!                 EffectiveHead(mm) = min(( wEta(mm) - wMidPt ), &
!                                           dir(mm)*(fEdn(iup(mm)) - fEup(idn(mm))))
!             else

!                 EffectiveHead(mm) = wCrown(mm) - wCrest(mm)

!             endif
!     endif
! end do

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_effective_head
!
!==========================================================================
!==========================================================================
!
subroutine weir_effective_crest_length &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, wHeight,  &
     wSideSlope, wEndContractions, EffectiveCrestLength, EffectiveHead, &
     thiscoef)
!
 character(64) :: subroutine_name = 'weir_effective_crest_length'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 real,              intent(in)      :: thiscoef
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

 real,  pointer                     :: wLength(:), wHeight(:), wEndContractions(:), wSideSlope(:) 
 real,  pointer                     :: wBreadth(:), EffectiveHead(:), EffectiveCrestLength(:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

 wLength    => elem2R(:,e2r_Length)
 wBreadth   => elem2R(:,e2r_BreadthScale)

! effective crest length is used in rectangular and trapezoidal weir flow calculation
 where ( (elem2I(:,e2i_elem_type) == eWeir ) .and. &
         (elem2I(:,e2i_geometry)  == eRectangular) )
    
    EffectiveCrestLength = max(wBreadth - 0.1 * wEndContractions * EffectiveHead, 0.0)

 elsewhere ( (elem2I(:,e2i_elem_type) == eWeir ) .and. &
             (elem2I(:,e2i_geometry)  == eTrapezoidal) )

    EffectiveCrestLength = wBreadth + twoR*wSideSlope*wHeight

 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_effective_crest_length
!
!==========================================================================
!==========================================================================
!
 subroutine villemonte_submergence_correction &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, wCrest, &
     subFactor1, subFactor2, Etadn, Etaup, upFace,  dnFace, maskarray)
!
 character(64) :: subroutine_name = 'villemonte_submergence_correction'

 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

 real,  pointer ::  wCrest(:), subFactor1(:), subFactor2(:)       
 real,  pointer ::  Etadn(:), Etaup(:)

 integer, pointer :: upFace(:), dnFace(:)
 logical, pointer :: maskarray(:)

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!% calculate the submergance factor for different weirs according to Villemonte 1974
 where ( (elem2I(:,e2i_weir_elem_type) == eVnotchWeir) .and. (maskarray) )     
    ! V-notch weir
    subFactor1 = (oneR - ((Etadn(upFace) - wCrest)/ (Etaup(dnFace) - wCrest)) &
        ** 2.5) ** 0.385

 elsewhere ( (elem2I(:,e2i_weir_elem_type) == eTrapezoidalWeir)  .and. (maskarray) )           
    ! Trapezoidal weir
    subFactor1 = (oneR - ((Etadn(upFace) - wCrest)/ (Etaup(dnFace) - wCrest)) &
        ** 2.5) ** 0.385
    subFactor2 = (oneR - ((Etadn(upFace) - wCrest)/ (Etaup(dnFace) - wCrest)) &
        ** 1.5) ** 0.385

 elsewhere ( (elem2I(:,e2i_weir_elem_type) == eTransverseWeir)  .and. (maskarray) )         
    ! Transverse weir
    subFactor1 = (oneR - ((Etadn(upFace) - wCrest)/ (Etaup(dnFace) - wCrest)) &
        ** 1.5) ** 0.385

 elsewhere ( (elem2I(:,e2i_weir_elem_type) == eSideFlowWeir)  .and. (maskarray) )     
   ! Side flow weir
   subFactor1 = (oneR - ((Etadn(upFace) - wCrest)/ (Etaup(dnFace) - wCrest)) &
        ** 1.67) ** 0.385

 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine villemonte_submergence_correction
!
!========================================================================== 
!==========================================================================
!
 subroutine weir_flow &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2new, &
     velocity2new, volumeMnew, velocityMnew, wFlow, wSideSlope,           &
     wCoeffTriangular, wCoeffRectangular, dir, EffectiveHead,             &
     EffectiveCrestLength, subFactor1, subFactor2, thiscoef)
!
 character(64) :: subroutine_name = 'weir_flow'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
 real,              intent(in)      :: thiscoef

 
 real,  pointer   ::  volume2new(:), velocity2new(:), volumeMnew(:), velocityMnew(:)
 real,  pointer   ::  wFlow(:), wSideSlope(:), wCoeffTriangular(:), wCoeffRectangular(:) 
 real,  pointer   ::  subFactor1(:), subFactor2(:)
 real,  pointer   ::  EffectiveHead(:), EffectiveCrestLength(:)

 integer, pointer :: dir(:)    

 integer :: mm
!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

 where ( elem2I(:,e2i_weir_elem_type) == eVnotchWeir )     
    ! V-notch weir
    ! fs = submergance_factor
    ! d = direction
    ! Q = d*fs*(Cw1SH^2.5)  
    wFlow        = dir * subFactor1 * wCoeffTriangular * wSideSlope * &
        EffectiveHead ** 2.5

    velocity2new = dir * subFactor1 * wCoeffTriangular * sqrt(abs(EffectiveHead))

    ! Volume = Q * dt
    volume2new   = thiscoef * wFlow 

 elsewhere ( elem2I(:,e2i_weir_elem_type) == eTrapezoidalWeir )          
    ! Trapezoidal weir
    ! Q = d*fs1*(Cw1SH^2.5) + d*fs2*(Cw2LH^1.5)
    wFlow        = dir * subFactor1 * (wCoeffTriangular * wSideSlope *   &
        EffectiveHead ** 2.5) + dir * subFactor2 * (wCoeffRectangular *  &
        EffectiveCrestLength * EffectiveHead ** 1.5)

    velocity2new = dir * sqrt(abs(EffectiveHead)) * (subFactor1 * &
        wCoeffTriangular + subFactor2 * wCoeffRectangular)

    ! Volume = Q * dt
    volume2new   = thiscoef * wFlow 

 elsewhere ( elem2I(:,e2i_weir_elem_type) == eTransverseWeir )        
    ! Transverse weir
    ! Q = d*fs2*(Cw2LH^1.5)
    wFlow        = dir * subFactor1 * wCoeffRectangular * &
        EffectiveCrestLength * EffectiveHead ** 1.5

    velocity2new = dir * subFactor1 * wCoeffRectangular * &
        sqrt(abs(EffectiveHead))

    ! Volume = Q * dt
    volume2new   = thiscoef * wFlow 

 elsewhere ( (elem2I(:,e2i_weir_elem_type) == eSideFlowWeir )  .and. &
             (dir .LE. zeroR) )    
    ! Side flow weir for reverse flow behaves like a Transverse weir
    ! Q = d*fs2*(Cw2LH^1.5)
    wFlow        = dir * subFactor1 * wCoeffRectangular * &
        EffectiveCrestLength * EffectiveHead ** 1.5

    velocity2new = dir * subFactor1 * wCoeffRectangular * &
        sqrt(abs(EffectiveHead))

    ! Volume = Q * dt
    volume2new   = thiscoef * wFlow 

 elsewhere ( (elem2I(:,e2i_weir_elem_type) == eSideFlowWeir )  .and. &
             (dir .GT. zeroR) )    
    ! Corrected formula (see Metcalf & Eddy, Inc., Wastewater Engineering, McGraw-Hill, 1972 p. 164).
    wFlow        = dir * subFactor1 * wCoeffRectangular * &
        (EffectiveCrestLength ** 0.83) * (EffectiveHead ** 1.67)

    velocity2new = dir * subFactor1 * wCoeffRectangular * &
        sqrt(abs(EffectiveHead))

    ! Volume = Q * dt
    volume2new   = thiscoef * wFlow 
 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_flow
!
!==========================================================================
!==========================================================================
!
 subroutine weir_surcharge_flow &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2new, &
     velocity2new, volumeMnew, velocityMnew, fEup, fEdn, iup, idn,        &
     wCrest, wCrown, wEta, wFlow, wCoeffOrif, dir,  EffectiveHead,        &
     thiscoef, maskarray)
!
 character(64) :: subroutine_name = 'weir_surcharge_flow'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
 real,              intent(in)      :: thiscoef

 
 real,    pointer ::  volume2new(:), velocity2new(:), volumeMnew(:), velocityMnew(:)
 real,    pointer ::  wFlow(:), wCoeffOrif(:)
 real,    pointer ::  wCrest(:), wCrown(:), wEta(:), EffectiveHead(:)
 real,    pointer ::  fEup(:), fEdn(:)

 integer, pointer ::  iup(:), idn(:), dir(:)
 logical, pointer ::  maskarray(:)

!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!% surchaged flow calculation
 where ( maskarray ) 
    ! for surcharged flow, head is calculated from the midpoint of the weir opening    
    EffectiveHead = min(( wEta - ((wCrown + wCrest) / twoR)), &
                          dir*(fEdn(iup) - fEup(idn)) )
    wFlow         = dir * wCoeffOrif * sqrt(abs(EffectiveHead))
    ! velocity equation needs correction
    velocity2new  = dir * wCoeffOrif * sqrt(abs(EffectiveHead))
    ! Volume is weir flow equation * dt (this case dt = thiscoef)
    volume2new    = thiscoef * wFlow 

 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

 end subroutine weir_surcharge_flow
!
!==========================================================================
!==========================================================================
!
subroutine weir_provisional_geometry &
    (elem2R, elemMR, faceR, elem2I, elemMI)
! this subroutine sets the weir geometry to zero.
 character(64) :: subroutine_name = 'weir_provisional_geometry'


 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)

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