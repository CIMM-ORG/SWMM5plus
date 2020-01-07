!==========================================================================
!
 module explicit_euler
!
    use adjustments
    use array_index
    use bc
    use data_keys
    use element_dynamics
    use element_geometry
    use face_values
    use globals
    use setting_definition
    use utility


    implicit none

    private

    public :: explicit_euler_advance
    public :: explicit_test_advance

    integer, parameter :: idummy = 0

    integer :: debuglevel = 0

 contains
!
!==========================================================================
!==========================================================================
!
 subroutine explicit_test_advance &
    (elem2R, elem2I, elem2YN, &
     elemMR, elemMI, elemMYN, &
     faceR,  faceI,  faceYN,  &
     bcdataDn, bcdataUp, thistime, dt)

 character(64) :: subroutine_name = 'explicit_euler_advance'

 real,      target, intent(in out) :: elem2R(:,:), elemMR(:,:), faceR(:,:)

 integer,   target, intent(in out) :: elem2I(:,:), elemMI(:,:), faceI(:,:)

 logical,   target, intent(in out) :: elem2YN(:,:), elemMYN(:,:), faceYN(:,:)

 type(bcType),  intent(in out)  :: bcdataDn(:), bcdataUp(:)

 real, intent(in) :: thistime, dt

 integer :: e2r_Volume_new, e2r_Velocity_new, eMr_Volume_new, eMr_Velocity_new

 real,  pointer :: newvolume(:), newvelocity(:)
 real,  pointer :: volume(:), velocity(:), eta(:)
 real,  pointer :: depth(:), zbottom(:), perimeter(:), mn(:), rh(:), area(:)
 real,  pointer :: timesUp(:), timesDn(:), breadth(:), length(:), flowrate(:)
 real,  pointer :: Qface(:), Vup(:), Vdn(:), Aup(:), Adn(:), Eup(:), Edn(:)

! TEST 20190102
! real,  pointer :: newarea(:), newflowrate(:)

 integer,   pointer :: fup(:), fdn(:), elemDn(:), elemUp(:), elemUp2(:)

 integer :: fitmp

!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

! set indexes for temporary space
 e2r_Volume_new = e2r_Temp(next_e2r_temparray)
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 e2r_Velocity_new = e2r_Temp(next_e2r_temparray)
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 fitmp = fi_Temp(next_fi_temparray)
 next_fi_temparray = utility_advance_temp_array (next_fi_temparray,fi_n_temp)

 eMr_Volume_new = 0
 eMr_Velocity_new = 0

!TEST 20190102
 newvolume   => elem2R(:,e2r_Volume_new)
 newvelocity => elem2R(:,e2r_Velocity_new)

!  newarea    => elem2R(:,e2r_Volume_new)
! newflowrate => elem2R(:,e2r_Velocity_new)


 elemUp2     => faceI(:,fitmp)

 volume   => elem2R(:,e2r_Volume)
 velocity => elem2R(:,e2r_Velocity)
 flowrate => elem2R(:,e2r_Flowrate)
 eta      => elem2R(:,e2r_Eta)
 mn       => elem2R(:,e2r_Roughness)
 rh       => elem2R(:,e2r_HydRadius)
 area     => elem2R(:,e2r_Area)
 depth    => elem2R(:,e2r_HydDepth)
 perimeter=> elem2R(:,e2r_Perimeter)
 zbottom  => elem2R(:,e2r_Zbottom)
 length   => elem2R(:,e2r_Length)
 breadth => elem2R(:,e2r_BreadthScale)
 timesUp => elem2R(:,e2r_Timescale_u)
 timesDn => elem2R(:,e2r_Timescale_d)

 Qface => faceR(:,fr_Flowrate)
 Vup   => faceR(:,fr_Velocity_u)
 Vdn   => faceR(:,fr_Velocity_d)
 Aup   => faceR(:,fr_Area_u)
 Adn   => faceR(:,fr_Area_d)
 Eup   => faceR(:,fr_Eta_u)
 Edn   => faceR(:,fr_Eta_d)

 fup => elem2I(:,e2i_Mface_u)
 fdn => elem2I(:,e2i_Mface_d)


!% using free surface as gA deta/dx
! where (elem2I(:,e2I_elem_type) == eChannel)
!    newvolume = volume + dt * (Qface(fup) - Qface(fdn))
!    newvelocity = (oneR / newvolume) &
!                * ( &
!                     velocity*volume &
!                    + dt * ( Qface(fup) * Vdn(fup) - Qface(fdn) * Vup(fdn)              &
!                            + grav * area * (Edn(fup) - Eup(fdn))                       &
!                            - grav * volume * (mn**2) * (velocity**2) / (rh**(4.0/3.0)) &
!                            ) &
!                  )
! endwhere
!
!% using free surface as T00
 where (elem2I(:,e2I_elem_type) == eChannel)
    newvolume = volume + dt * (Qface(fup) - Qface(fdn))
!    newvelocity = (oneR / newvolume) &
!                * ( &
!                     velocity*volume &
!                    + dt * ( Qface(fup) * Vdn(fup) - Qface(fdn) * Vup(fdn)              &
!                            + grav * ( Adn(fup)*Edn(fup) - Aup(fdn)*Eup(fdn))           &
!                            + grav * ( Aup(fdn) - Adn(fup) ) * eta                      &
!                            - grav * volume * (mn**2) * (velocity**2) / (rh**(4.0/3.0)) &
!                            ) &
!                  )
    newvelocity = (oneR / (newvolume*(oneR + dt * grav *  (mn**2) * velocity / (rh**(4.0/3.0)) ))) &
                * ( &
                     velocity*volume &
                    + dt * ( Qface(fup) * Vdn(fup) - Qface(fdn) * Vup(fdn)              &
                            + grav * ( Adn(fup)*Edn(fup) - Aup(fdn)*Eup(fdn))           &
                            + grav * ( Aup(fdn) - Adn(fup) ) * eta                      &
                            ) &
                   )
 endwhere
! AREA/FLOWRATE TEST 20180102
! print *, 'flow start :',flowrate(1:3)
! print *, 'qface start:',Qface(1:3)
! print *, 'Edn start:',Edn(1:3)
! print *, 'Eup start:',Eup(1:3)


! where (elem2I(:,e2I_elem_type) == eChannel)
!    newarea = area + dt * (Qface(fup) - Qface(fdn)) / length   ! now area
!    ! now flowrate
!    newflowrate =  flowrate + &
!                ( ( &
!                    + dt * ( Qface(fup) * Vdn(fup) - Qface(fdn) * Vup(fdn)              &
!                            + grav * ( Adn(fup)*Edn(fup) - Aup(fdn)*Eup(fdn))           &
!                            + grav * ( Aup(fdn) - Adn(fup) ) * eta                      &
!                            ) &
!                  ) / length) - dt * grav * (mn**2) * (flowrate**2) / ( area * rh**(4.0/3.0))
! endwhere
!
! where (elem2I(:,e2I_elem_type) == eChannel)
!    newarea = area + dt * (Qface(fup) - Qface(fdn)) / length   ! now area
!    ! now flowrate
!    newflowrate =  (flowrate + &
!                ( ( &
!                    + dt * ( Qface(fup) * Vdn(fup) - Qface(fdn) * Vup(fdn)              &
!                            + grav * ( Adn(fup)*Edn(fup) - Aup(fdn)*Eup(fdn))           &
!                            + grav * ( Aup(fdn) - Adn(fup) ) * eta                      &
!                            ) &
!                  ) / length)) / (1 + dt * grav * (mn**2) * (flowrate) / ( area * rh**(4.0/3.0))  )
! endwhere
!



! print *, 'new flow  ',newflowrate(1:3)
! !print *, 'qface:    ',Qface(size(Qface)-3:size(Qface))
! !print *, 'flowrate: ',flowrate(size(flowrate)-3:size(flowrate))
! where (elem2I(:,e2I_elem_type) == eChannel)
!    newarea = dt * ( Qface(fup) * Vdn(fup) - Qface(fdn) * Vup(fdn) ) / length
! endwhere
! print *, '1 line    ',newarea(1:3)
! where (elem2I(:,e2I_elem_type) == eChannel)
!   ! newarea = dt * grav * ( Adn(fup)*Edn(fup) - Aup(fdn)*Eup(fdn)) / length
!   newarea = Aup(fdn)*Eup(fdn)
! endwhere
!! print *, '2 line    ',newarea(1:3)
! print *, fup(1:3)
! print *, fdn(1:3)
! print *, 'Adn ',Adn(1:3)
! print *, 'Aup ',Aup(1:3)
! print *, 'Edn ',Edn(1:3)
! print *, 'Eup ',Eup(1:3)
!! print *, Adn(fup(1:3))
!! print *, Edn(fup(1:3))
!! print *, Aup(fdn(1:3))
!! print *, Eup(fdn(1:3))
!! print *, fdn(1:3)
!! print *, fup(1:3)
!stop
! where (elem2I(:,e2I_elem_type) == eChannel)
!    newarea = dt * grav * ( Aup(fdn) - Adn(fup) ) * eta  / length
! endwhere
! print *, '3 line    ',newarea(1:3)
! where (elem2I(:,e2I_elem_type) == eChannel)
!    newarea = dt * grav * area * (mn**2) * (velocity**2) / (rh**(4.0/3.0))
! endwhere
! print *, '4 line    ',newarea(1:3)
! print *, 'new flow  ',newflowrate(1:3)
!
! stop
!
!! AREA/FLOWRATE TEST 20180102

!
!%  using free surface as T10
! where (elem2I(:,e2I_elem_type) == eChannel)
!    newvolume = volume + dt * (Qface(fup) - Qface(fdn))
!    newvelocity = (oneR / newvolume) &
!                * ( &
!                     velocity*volume &
!                    + dt * ( Qface(fup) * Vdn(fup) - Qface(fdn) * Vup(fdn)              &
!                            + onehalfR * grav * ( Adn(fup)*Edn(fup) - Aup(fdn)*Eup(fdn) )           &
!                            + onehalfR * grav * ( Aup(fdn)*Edn(fup) - Adn(fup)*Eup(fdn) ) &
!                            - grav * volume * (mn**2) * (velocity**2) / (rh**(4.0/3.0)) &
!                            ) &
!                  )
! endwhere
!
!% update geometry and dynamics on element
 where (elem2I(:,e2I_elem_type) == eChannel)
    volume    = newvolume
    velocity  = newvelocity
    area      = volume / length
    flowrate  = velocity * area
    depth     = area / breadth
    perimeter = twoR * depth + breadth
    eta       = depth + zbottom
    rh        = area / perimeter
    !timesUp   = -onehalfR * length / (velocity - sqrt(grav * depth))
    !timesDn   = +onehalfR * length / (velocity + sqrt(grav * depth))
 endwhere

!! AREA/FLOWRATE TEST 20180102
! where (elem2I(:,e2I_elem_type) == eChannel)
!    volume    = newarea * length !test
!    velocity  = newflowrate / newarea
!    area      = newarea !renaming
!    flowrate  = newflowrate
!    depth     = area / breadth
!    perimeter = twoR * depth + breadth
!    eta       = depth + zbottom
!    rh        = area / perimeter
!    !timesUp   = -onehalfR * length / (velocity - sqrt(grav * depth))
!    !timesDn   = +onehalfR * length / (velocity + sqrt(grav * depth))
! endwhere
!!! AREA/FLOWRATE TEST 20180102

! BC updates
!call bc_applied_onelement (elem2R, bcdataDn, bcdataUp, thistime+dt, bc_category_elevation,idummy)
 call bc_applied_onelement (elem2R, bcdataDn, bcdataUp, thistime+dt, bc_category_inflowrate,e2r_Velocity)
 call bc_applied_onelement (elem2R, bcdataDn, bcdataUp, thistime+dt, bc_category_elevation, idummy)

! print *, 'after bc on element'
! print *, 'new flow  ',flowrate(size(flowrate)-3:size(flowrate))
! print *, 'Qface     ',Qface(size(Qface)-3:size(Qface))
!
 call test_face_update (elem2R, elem2I, faceR, faceI, elemUp2)

! print *, 'after face update'
! print *, 'new flow  ',flowrate(size(flowrate)-3:size(flowrate))
! print *, 'Qface     ',Qface(size(Qface)-3:size(Qface))

! call quadratic_face_update (elem2R, elem2I, faceR, faceI, elemUp2)

! BC applied on face
 call bc_applied_onface (faceR, faceI, elem2R, elem2I, bcdataDn, bcdataUp, e2r_Velocity_new, thistime+dt)

! print *, 'after bc on face'
! print *, 'new flow  ',flowrate(size(flowrate)-3:size(flowrate))
! print *, 'Qface     ',Qface(size(Qface)-3:size(Qface))

! BC extrapolated - moved into bc_applied_onface
! call bc_face_othervalues (faceR, faceI, elem2R, bcdataDn)
! call bc_face_othervalues (faceR, faceI, elem2R, bcdataUp)

 if (setting%Method%AdjustVshapedFlowrate%Apply) then
    call adjust_Vshaped_flowrate (elem2R, faceR, elem2I, elem2YN)
 endif

 ! set velocities and upstream values on faces (without hydraulic jump)
 call adjust_face_dynamic_limits &
    (faceR, faceI, newvolume, newvolume, &
     ( (faceI(:,fi_etype_u) == eChannel) .and. (faceI(:,fi_etype_d) == eChannel) ), .false. )

 !call bc_applied (faceR, bcdataDn, bcdataUp, thistime+dt)

 !print *, elem2R(elemDn ,e2r_Flowrate)
 !print *, faceR(:,fr_Flowrate)
 !print *, faceR(:,fr_Area_u)
 !print *, faceR(:,fr_Eta_u)
 !print *, faceR(:,fr_Velocity_u)
 !print *, elem2R(elemUp ,e2r_Flowrate)
 !print *, elem2R(elemUp2,e2r_Flowrate)


! reset temporary data space
 elem2R(:,e2r_Volume_new)    = nullvalueR
 elem2R(:,e2r_Velocity_new)  = nullvalueR

 next_e2r_temparray = next_e2r_temparray - 2
 next_fi_temparray  = next_fi_temparray - 1


 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine explicit_test_advance
!
!==========================================================================
!==========================================================================
!
 subroutine test_face_update &
    (elem2R, elem2I, faceR, faceI, elemUp2)

 character(64) :: subroutine_name = 'test_face_update'

 real,      target,     intent(in out)  :: elem2R(:,:), faceR(:,:)
 integer,   target,     intent(in out)  :: elem2I(:,:), faceI(:,:)

 integer,                  intent(in out)  :: elemUp2(:)

 integer,   pointer :: elemUp(:), elemDn(:)
 real,      pointer :: tscaleUp(:), tscaleDn(:), wavespeed(:)

 integer    :: e2r_wavespeed

!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

 e2r_wavespeed = e2r_Temp(next_e2r_temparray)
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

!% compute wavespeed
 wavespeed => elem2R(:,e2r_wavespeed)
 wavespeed = sqrt(grav * elem2R(:,e2r_HydDepth))

 tscaleUp => elem2R(:,e2r_Timescale_u)
 tscaleDn => elem2R(:,e2r_Timescale_d)

 tscaleUp = + (onehalfR * elem2R(:,e2r_Length) ) / (wavespeed - elem2R(:,e2r_Velocity))
 tscaleDn = + (onehalfR * elem2R(:,e2r_Length) ) / (wavespeed + elem2R(:,e2r_Velocity))

 where (tscaleUp < zeroR )
    tscaleUp = setting%Limiter%Timescale%Maximum
 endwhere
 where (tscaleDn < zeroR )
    tscaleDn = setting%Limiter%Timescale%Maximum
 endwhere

 where (tscaleUp < setting%Limiter%Timescale%Minimum)
    tscaleUp = setting%Limiter%Timescale%Minimum
 endwhere
 where (tscaleDn < setting%Limiter%Timescale%Minimum)
    tscaleDn = setting%Limiter%Timescale%Minimum
 endwhere

 where (tscaleUp > setting%Limiter%Timescale%Maximum)
    tscaleUp = setting%Limiter%Timescale%Maximum
 endwhere
 where (tscaleDn > setting%Limiter%Timescale%Maximum)
    tscaleDn = setting%Limiter%Timescale%Maximum
 endwhere

 !wavespeed = tscaleUp + tscaleDn

 !tscaleUp = tscaleUp/wavespeed
 !tscaleDn = tscaleDn/wavespeed

!% get map to upstream for a quadratic interp
 elemUp => faceI(:,fi_Melem_u)
 elemDn => faceI(:,fi_Melem_d)

 where (elemUp > size(elemUp))
    elemUp = size(elemUp)
 endwhere

 elemUp2 = elemUp

 where ((faceI(:,fi_etype_u) == fChannel) .and. (faceI(:,fi_etype_d) == fChannel))
    elemUp2 = elem2I(elemUp,e2i_Mface_u)
 endwhere

 where (elemUp2 > size(elemUp2))
    elemUp2 = size(elemUp2)
 endwhere

 where (elemUp2 < 1)
    elemUp2 = 1
 endwhere

 elemUp2 = faceI(elemUp2,fi_Melem_u)

 where (elemUp2 > size(elemUp2))
    elemUp2 = size(elemUp2)
 endwhere

! print *, elemDn
! print *, elemUp
! print *, elemUp2

! print *, size(faceI,1), size(elemDn), size(elemUp), size(elemUp2)

 where( (faceI(:,fi_etype_u) == fChannel) .and. (faceI(:,fi_etype_d) == fChannel))

    faceR(:,fr_Flowrate)   = 0.375 * elem2R(elemDn ,e2r_Flowrate) &
                           + 0.75  * elem2R(elemUp ,e2r_Flowrate) &
                           - 0.125 * elem2R(elemUp2,e2r_Flowrate)

!    faceR(:,fr_Flowrate)   = 0.5 * elem2R(elemDn  ,e2r_Flowrate) &
!                           + 0.5  * elem2R(elemUp ,e2r_Flowrate)

!    faceR(:,fr_Flowrate) = (    tscaleDn(elemUp) *  elem2R(elemDn ,e2r_Flowrate)  &
!                             +  tscaleUp(elemDn) *  elem2R(elemUp ,e2r_Flowrate)) &
!                             /( tscaleDn(elemUp) + tscaleUp(elemDn))


    faceR(:,fr_Velocity_d) = 0.375 * elem2R(elemDn ,e2r_Velocity) &
                           + 0.75  * elem2R(elemUp ,e2r_Velocity) &
                           - 0.125 * elem2R(elemUp2,e2r_Velocity)

!    faceR(:,fr_Velocity_d) = 0.5 * elem2R(elemDn ,e2r_Velocity) &
!                           + 0.5  * elem2R(elemUp ,e2r_Velocity)

!    faceR(:,fr_Velocity_d) = (  tscaleDn(elemDn) *  elem2R(elemDn ,e2r_Velocity)  &
!                              + tscaleUp(elemUp) *  elem2R(elemUp ,e2r_Velocity)) &
!                             /( tscaleDn(elemUp) + tscaleUp(elemDn))



    faceR(:,fr_Area_d)     = 0.375 * elem2R(elemDn ,e2r_Area) &
                           + 0.75  * elem2R(elemUp ,e2r_Area) &
                           - 0.125 * elem2R(elemUp2,e2r_Area)

!    faceR(:,fr_Area_d)     = 0.5 * elem2R(elemDn ,e2r_Area) &
!                           + 0.5  * elem2R(elemUp ,e2r_Area)

!    faceR(:,fr_Area_d)   =  (   tscaleDn(elemUp) *  elem2R(elemDn ,e2r_Area)  &
!                             +  tscaleUp(elemDn) *  elem2R(elemUp ,e2r_Area) ) &
!                             /( tscaleDn(elemUp) + tscaleUp(elemDn))


    !faceR(:,fr_Eta_d)      = 0.375 * elem2R(elemDn ,e2r_Eta) &
    !                       + 0.75  * elem2R(elemUp ,e2r_Eta) &
    !                       - 0.125 * elem2R(elemUp2,e2r_Eta)

    faceR(:,fr_Eta_d)      = 0.5  * elem2R(elemDn ,e2r_Eta) &
                           + 0.5  * elem2R(elemUp ,e2r_Eta)



    faceR(:,fr_Topwidth)   = 0.375 * elem2R(elemDn ,e2r_Topwidth) &
                           + 0.75  * elem2R(elemUp ,e2r_Topwidth) &
                           - 0.125 * elem2R(elemUp2,e2r_Topwidth)

!     faceR(:,fr_Topwidth)   = 0.5 * elem2R(elemDn ,e2r_Topwidth) &
!                           + 0.5  * elem2R(elemUp ,e2r_Topwidth)

!    faceR(:,fr_Topwidth)  = (   tscaleDn(elemUp) *  elem2R(elemDn ,e2r_Topwidth)   &
!                             +  tscaleUp(elemDn) *  elem2R(elemUp ,e2r_Topwidth) ) &
!                             /( tscaleDn(elemUp) + tscaleUp(elemDn))


    faceR(:,fr_Eta_u)      = faceR(:,fr_Eta_d)
    faceR(:,fr_Area_u)     = faceR(:,fr_Area_d)
    faceR(:,fr_Velocity_u) = faceR(:,fr_Velocity_d)
 endwhere

!% reset temporary space
 elem2R(:,e2r_wavespeed)  = nullvalueR
 next_e2r_temparray = next_e2r_temparray - 1

 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine test_face_update
!
!==========================================================================
!==========================================================================
!
 subroutine quadratic_face_update &
    (elem2R, elem2I, faceR, faceI, elemUp2)

 character(64) :: subroutine_name = 'quadratic_face_update'

 real,      target,     intent(in out)  :: elem2R(:,:), faceR(:,:)
 integer,   target,     intent(in out)  :: elem2I(:,:), faceI(:,:)

 integer,                  intent(in out)  :: elemUp2(:)

 integer,   pointer :: elemUp(:), elemDn(:)

!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

! get map to upstream for a quadratic interp
 elemUp => faceI(:,fi_Melem_u)
 elemDn => faceI(:,fi_Melem_d)

 where (elemUp > size(elemUp))
    elemUp = size(elemUp)
 endwhere

 elemUp2 = elemUp

 where ((faceI(:,fi_etype_u) == fChannel) .and. (faceI(:,fi_etype_d) == fChannel))
    elemUp2 = elem2I(elemUp,e2i_Mface_u)
 endwhere

 where (elemUp2 > size(elemUp2))
    elemUp2 = size(elemUp2)
 endwhere

 where (elemUp2 < 1)
    elemUp2 = 1
 endwhere

 elemUp2 = faceI(elemUp2,fi_Melem_u)

 where (elemUp2 > size(elemUp2))
    elemUp2 = size(elemUp2)
 endwhere

! print *, elemDn
! print *, elemUp
! print *, elemUp2

! print *, size(faceI,1), size(elemDn), size(elemUp), size(elemUp2)

 where( (faceI(:,fi_etype_u) == fChannel) .and. (faceI(:,fi_etype_d) == fChannel))

    faceR(:,fr_Flowrate)   = 0.375 * elem2R(elemDn ,e2r_Flowrate) &
                           + 0.75  * elem2R(elemUp ,e2r_Flowrate) &
                           - 0.125 * elem2R(elemUp2,e2r_Flowrate)

    faceR(:,fr_Velocity_d) = 0.375 * elem2R(elemDn ,e2r_Velocity) &
                           + 0.75  * elem2R(elemUp ,e2r_Velocity) &
                           - 0.125 * elem2R(elemUp2,e2r_Velocity)

    faceR(:,fr_Area_d)     = 0.375 * elem2R(elemDn ,e2r_Area) &
                           + 0.75  * elem2R(elemUp ,e2r_Area) &
                           - 0.125 * elem2R(elemUp2,e2r_Area)

    !faceR(:,fr_Eta_d)      = 0.375 * elem2R(elemDn ,e2r_Eta) &
    !                       + 0.75  * elem2R(elemUp ,e2r_Eta) &
    !                       - 0.125 * elem2R(elemUp2,e2r_Eta)
    faceR(:,fr_Eta_d)      = 0.5  * elem2R(elemDn ,e2r_Eta) &
                           + 0.5  * elem2R(elemUp ,e2r_Eta)

    faceR(:,fr_Topwidth)   = 0.375 * elem2R(elemDn ,e2r_Topwidth) &
                           + 0.75  * elem2R(elemUp ,e2r_Topwidth) &
                           - 0.125 * elem2R(elemUp2,e2r_Topwidth)

    faceR(:,fr_Eta_u)      = faceR(:,fr_Eta_d)
    faceR(:,fr_Area_u)     = faceR(:,fr_Area_d)
    faceR(:,fr_Velocity_u) = faceR(:,fr_Velocity_d)
 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine quadratic_face_update
!
!==========================================================================
!==========================================================================
!
 subroutine explicit_euler_advance &
    (elem2R, elem2I, elem2YN, &
     elemMR, elemMI, elemMYN, &
     faceR,  faceI,  faceYN,  &
     bcdataDn, bcdataUp, thistime, dt)

 character(64) :: subroutine_name = 'explicit_euler_advance'

 real,      target, intent(in out) :: elem2R(:,:), elemMR(:,:), faceR(:,:)

 integer,   target, intent(in out) :: elem2I(:,:), elemMI(:,:), faceI(:,:)

 logical,   target, intent(in out) :: elem2YN(:,:), elemMYN(:,:), faceYN(:,:)

 type(bcType),  intent(in out)  :: bcdataDn(:), bcdataUp(:)

 real, intent(in) :: thistime, dt

 integer :: e2r_Volume_new, e2r_Velocity_new, eMr_Volume_new, eMr_Velocity_new

 real,  pointer :: newvolume(:), newvelocity(:), volume(:), velocity(:), eta(:)
 real,  pointer :: depth(:), zbottom(:), perimeter(:), mn(:), rh(:), area(:)
 real,  pointer :: timesUp(:), timesDn(:), breadth(:), length(:), flowrate(:)
 real,  pointer :: Qface(:), Vup(:), Vdn(:), Aup(:), Adn(:), Eup(:), Edn(:)


 integer,   pointer :: fup(:), fdn(:), elemDn(:), elemUp(:), elemUp2(:)

 integer :: fitmp


!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

! set indexes for temporary space
 e2r_Volume_new = e2r_Temp(next_e2r_temparray)
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 e2r_Velocity_new = e2r_Temp(next_e2r_temparray)
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 fitmp = fi_Temp(next_fi_temparray)
 next_fi_temparray = utility_advance_temp_array (next_fi_temparray,fi_n_temp)


 eMr_Volume_new = 0
 eMr_Velocity_new = 0

 newvolume   => elem2R(:,e2r_Volume_new)
 newvelocity => elem2R(:,e2r_Velocity_new)

 elemUp2     => faceI(:,fitmp)

 volume   => elem2R(:,e2r_Volume)
 velocity => elem2R(:,e2r_Velocity)
 flowrate => elem2R(:,e2r_Flowrate)
 eta      => elem2R(:,e2r_Eta)
 mn       => elem2R(:,e2r_Roughness)
 rh       => elem2R(:,e2r_HydRadius)
 area     => elem2R(:,e2r_Area)
 depth    => elem2R(:,e2r_HydDepth)
 perimeter=> elem2R(:,e2r_Perimeter)
 zbottom  => elem2R(:,e2r_Zbottom)
 length   => elem2R(:,e2r_Length)
 breadth => elem2R(:,e2r_BreadthScale)
 timesUp => elem2R(:,e2r_Timescale_Q_u)
 timesDn => elem2R(:,e2r_Timescale_Q_d)

 Qface => faceR(:,fr_Flowrate)
 Vup   => faceR(:,fr_Velocity_u)
 Vdn   => faceR(:,fr_Velocity_d)
 Aup   => faceR(:,fr_Area_u)
 Adn   => faceR(:,fr_Area_d)
 Eup   => faceR(:,fr_Eta_u)
 Edn   => faceR(:,fr_Eta_d)

 fup => elem2I(:,e2i_Mface_u)
 fdn => elem2I(:,e2i_Mface_d)

 where (elem2I(:,e2I_elem_type) == eChannel)
    newvolume = volume + dt * (Qface(fup) - Qface(fdn))
    newvelocity = (oneR / newvolume) &
                * ( &
                     velocity*volume &
                    + dt * ( Qface(fup) * Vdn(fup) - Qface(fdn) * Vup(fdn)              &
                            + grav * area * (Edn(fup) - Eup(fdn))                       &
                            - grav * volume * (mn**2) * (velocity**2) / (rh**(4.0/3.0)) &
                            ) &
                  )
 endwhere

 where (elem2I(:,e2I_elem_type) == eChannel)
    volume    = newvolume
    velocity  = newvelocity
    area      = volume / length
    flowrate  = velocity * area
    depth     = area / breadth
    perimeter = twoR * depth + breadth
    eta       = depth + zbottom
    rh        = area / perimeter
    timesUp   = -onehalfR * length / (velocity - sqrt(grav * depth))
    timesDn   = +onehalfR * length / (velocity + sqrt(grav * depth))
 endwhere

 call bc_applied_onelement (elem2R, bcdataDn, bcdataUp, thistime+dt, bc_category_elevation,idummy)
 call bc_applied_onelement (elem2R, bcdataDn, bcdataUp, thistime+dt, bc_category_inflowrate,e2r_Velocity_new)

! get map to upstream for a quadratic interp
 elemUp => faceI(:,fi_Melem_u)
 elemDn => faceI(:,fi_Melem_d)

 where (elemUp > size(elemUp))
    elemUp = size(elemUp)
 endwhere

 elemUp2 = elemUp

 where ((faceI(:,fi_etype_u) == fChannel) .and. (faceI(:,fi_etype_d) == fChannel))
    elemUp2 = elem2I(elemUp,e2i_Mface_u)
 endwhere

 where (elemUp2 > size(elemUp2))
    elemUp2 = size(elemUp2)
 endwhere

 where (elemUp2 < 1)
    elemUp2 = 1
 endwhere

 elemUp2 = faceI(elemUp2,fi_Melem_u)

 where (elemUp2 > size(elemUp2))
    elemUp2 = size(elemUp2)
 endwhere

 where( (faceI(:,fi_etype_u) == fChannel) .and. (faceI(:,fi_etype_d) == fChannel))

    faceR(:,fr_Flowrate)   = 0.375 * elem2R(elemDn ,e2r_Flowrate) &
                           + 0.75  * elem2R(elemUp ,e2r_Flowrate) &
                           - 0.125 * elem2R(elemUp2,e2r_Flowrate)

    faceR(:,fr_Velocity_d) = 0.375 * elem2R(elemDn ,e2r_Velocity) &
                           + 0.75  * elem2R(elemUp ,e2r_Velocity) &
                           - 0.125 * elem2R(elemUp2,e2r_Velocity)

    faceR(:,fr_Area_d)     = 0.375 * elem2R(elemDn ,e2r_Area) &
                           + 0.75  * elem2R(elemUp ,e2r_Area) &
                           - 0.125 * elem2R(elemUp2,e2r_Area)

    !faceR(:,fr_Eta_d)      = 0.375 * elem2R(elemDn ,e2r_Eta) &
    !                       + 0.75  * elem2R(elemUp ,e2r_Eta) &
    !                       - 0.125 * elem2R(elemUp2,e2r_Eta)
    faceR(:,fr_Eta_d)      = 0.5  * elem2R(elemDn ,e2r_Eta) &
                           + 0.5  * elem2R(elemUp ,e2r_Eta)

    faceR(:,fr_Topwidth)   = 0.375 * elem2R(elemDn ,e2r_Topwidth) &
                           + 0.75  * elem2R(elemUp ,e2r_Topwidth) &
                           - 0.125 * elem2R(elemUp2,e2r_Topwidth)

    faceR(:,fr_Eta_u)      = faceR(:,fr_Eta_d)
    faceR(:,fr_Area_u)     = faceR(:,fr_Area_d)
    faceR(:,fr_Velocity_u) = faceR(:,fr_Velocity_d)
 endwhere

 call bc_applied_onface (faceR, faceI, elem2R, elem2I, bcdataDn, bcdataUp, e2r_Velocity_new, thistime+dt)

 if (setting%Method%AdjustVshapedFlowrate%Apply) then
    call adjust_Vshaped_flowrate (elem2R, faceR, elem2I, elem2YN)
 endif


 !call bc_applied (faceR, bcdataDn, bcdataUp, thistime+dt)

 !call bc_face_othervalues (faceR, faceI, elem2R, bcdataDn)
 !call bc_face_othervalues (faceR, faceI, elem2R, bcdataUp)


! reset temporary data space
 elem2R(:,e2r_Volume_new)    = nullvalueR
 elem2R(:,e2r_Velocity_new)  = nullvalueR

 next_e2r_temparray = next_e2r_temparray - 2
 next_fi_temparray  = next_fi_temparray - 1

 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine explicit_euler_advance
!
!==========================================================================
! END OF MODULE explicit_euler
!==========================================================================
 end module explicit_euler