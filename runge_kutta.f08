! module runge_kutta
!
! runge-kutta time advance for a single time step
!
!==========================================================================
!
 module runge_kutta

    use adjustments
    use array_index
    use bc
    use data_keys
    use diagnostic
    use element_geometry
    use element_dynamics
    use face_values
    use globals
    use setting_definition
    use utility
    use weir

    implicit none

    private

    public :: rk2

    integer :: debuglevel = 0

 contains
!==========================================================================
!==========================================================================
!
 subroutine rk2 &
    (elem2R, elemMR, elem2I, elemMI, faceR, faceI, elem2YN, elemMYN, faceYN, &
     bcdataDn, bcdataUp, thistime, dt)
!
! runge-kutta time advance for a single time step
!
 character(64) :: subroutine_name = 'rk2'

 real,      target, intent(in out) :: elem2R(:,:),  elemMR(:,:),  faceR(:,:)
 integer,   target, intent(in out) :: elem2I(:,:),  elemMI(:,:),  faceI(:,:)
 logical,   target, intent(in out) :: elem2YN(:,:), elemMYN(:,:), faceYN(:,:)
 type(bcType),      intent(in out) :: bcdataDn(:),  bcdataUp(:)
 real,              intent(in)     :: thistime, dt

 integer,   pointer :: fdn(:), fup(:)

 integer :: e2r_Volume_new, e2r_Velocity_new,  eMr_Volume_new, eMr_Velocity_new
 integer :: ii
 real    :: thiscoef(2), steptime

 integer :: ilink
    

!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!%  set indexes for temporary space
 e2r_Volume_new = e2r_Temp(next_e2r_temparray)
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 e2r_Velocity_new = e2r_Temp(next_e2r_temparray)
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 eMr_Volume_new = eMr_Temp(next_eMr_temparray)
 next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

 eMr_Velocity_new = eMr_Temp(next_eMr_temparray)
 next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

!%  zero out the temporary space
 elem2R(:,e2r_Volume_new)    = zeroR
 elem2R(:,e2r_Velocity_new)  = zeroR
 elemMR(:,eMr_Volume_new)    = zeroR
 elemMR(:,eMr_Velocity_new)  = zeroR
    
 if (  count(elem2I(:,e2i_elem_type) == eChannel) &
     + count(elemMI(:,eMi_elem_type) == eJunctionChannel) > zeroI) then
         ! + count(elem2I(:,e2i_elem_type) == eWeir) & 

    !%  coefficients for the rk2 steps
    thiscoef(1) = onehalfR
    thiscoef(2) = oneR    
    
    !%  step through the two steps of the RK2
    do ii=1,2
        steptime = thistime + thiscoef(ii) * dt



        call sve_rk2_step &
            (e2r_Volume, e2r_Velocity, eMr_Volume, eMr_Velocity, &
             e2r_Volume_new, e2r_Velocity_new, eMr_Volume_new, eMr_Velocity_new, &
             elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
             thiscoef(ii))

        call weir_step &
            (e2r_Volume, e2r_Velocity, eMr_Volume, eMr_Velocity, &
             e2r_Volume_new, e2r_Velocity_new, eMr_Volume_new, eMr_Velocity_new, &
             elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
             thiscoef(ii))
             
 
        call rk2_update_auxiliary_variables &
            (e2r_Velocity_new, eMr_Velocity_new, e2r_Volume_new, eMr_Volume_new, &
             elem2R, elem2I, elem2YN, &
             elemMR, elemMI, elemMYN, &
             faceR,  faceI,  faceYN,  &
             bcdataDn, bcdataUp, steptime, ii)
             
  
    
        if (ii==1) then
            !% store the net face fluxes that are used for volume advance.
            call diagnostic_element_volume_conservation_fluxes &
                (elem2R, elem2I, elemMR, elemMI, faceR)
        endif
    end do
    


    !% compute local element-based volume conservation
    call diagnostic_element_volume_conservation &
        (elem2R, elem2I, elemMR, elemMI, e2r_Volume_new, eMr_Volume_new)

    !%  update the velocity and volume
    call overwrite_old_values &
        (elem2R, elem2I, e2r_Velocity, e2r_Velocity_new, &
         e2r_Volume, e2r_Volume_new, e2i_elem_type, eChannel, .true.)

    call overwrite_old_values &
        (elemMR, elemMI, eMr_Velocity, eMr_Velocity_new, &
         eMr_Volume, eMr_Volume_new, eMi_elem_type, eJunctionChannel, .false.)

endif

!%  reset temporary data space
 elem2R(:,e2r_Volume_new)    = nullvalueR
 elem2R(:,e2r_Velocity_new)  = nullvalueR
 elemMR(:,eMr_Volume_new)    = nullvalueR
 elemMR(:,eMr_Velocity_new)  = nullvalueR

 next_e2r_temparray = next_e2r_temparray - 2
 next_eMr_temparray = next_eMr_temparray - 2

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine rk2
!
!==========================================================================
!==========================================================================
!
 subroutine sve_rk2_step &
    (e2r_Volume_old, e2r_Velocity_old, eMr_Volume_old, eMr_Velocity_old, &
     e2r_Volume_new, e2r_Velocity_new, eMr_Volume_new, eMr_Velocity_new, &
     elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     thiscoef)
!
! time advance - one step of RK2
!
 character(64) :: subroutine_name = 'sve_rk2_step'

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

 real,  pointer ::  fQ(:), fUdn(:), fUup(:), fAdn(:), fAup(:), fEdn(:), fEup(:)
 real,  pointer ::  kc2(:), ku2(:), ones2r(:)
 real,  pointer ::  kcM(:), kuM(:), onesMr(:)
 real,  pointer ::  volume2old(:), volume2new(:), velocity2old(:), velocity2new(:)
 real,  pointer ::  volumeMold(:), volumeMnew(:), velocityMold(:), velocityMnew(:)
 real,  pointer ::  eta2(:), etaM(:), rh2(:), mn2(:), rhM(:), mnM(:)

 integer, pointer       :: iup(:), idn(:)

 integer :: mm

!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!%  pointers for volume and velocity storage (updating)
 volume2old   => elem2R(:,e2r_Volume_old)
 volume2new   => elem2R(:,e2r_Volume_new)
 velocity2old => elem2R(:,e2r_Velocity_old)
 velocity2new => elem2R(:,e2r_Velocity_new)

 volumeMold   => elemMR(:,eMr_Volume_old)
 volumeMnew   => elemMR(:,eMr_Volume_new)
 velocityMold => elemMR(:,eMr_Velocity_old)
 velocityMnew => elemMR(:,eMr_Velocity_new)

 eta2      => elem2R(:,e2r_Eta)
 rh2       => elem2R(:,e2r_HydRadius)
 mn2       => elem2R(:,e2r_Roughness)

 etaM      => elemMR(:,eMr_Eta)
 rhM       => elemMR(:,eMr_HydRadius)
 mnM       => elemMR(:,eMr_Roughness)

!%  temporary space for channel elements
 kc2 => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 ku2 => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

 ones2r => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

!%  temporary space for juctions
 kcM => elemMR(:,eMr_Temp(next_eMr_temparray))
 next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

 kuM => elemMR(:,eMr_Temp(next_eMr_temparray))
 next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

 onesMr=> elemMR(:,eMr_Temp(next_eMr_temparray))
 next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

!%  zero temporary arrays
 kc2 = zeroR
 ku2 = zeroR

 kcM = zeroR
 kuM = zeroR

!%  initialize ones as 1.0
 ones2r = oneR
 onesMr = oneR

!%  pointers for convenience in notation
 fQ     => faceR(:,fr_Flowrate)
 fUdn   => faceR(:,fr_Velocity_d) ! velocity on the downstream side of a face
 fUup   => faceR(:,fr_Velocity_u) ! velocity on the upstream side of a face
 fAdn   => faceR(:,fr_Area_d)
 fAup   => faceR(:,fr_Area_u)
 fEdn   => faceR(:,fr_Eta_d)
 fEup   => faceR(:,fr_Eta_u)

!%  SOURCE TERMS ----------------------------------

!%  Channel elements (one upstream and one downstream face)
!%  Using the T00 method of Hodges & Liu, 2018 for momentum

 iup => elem2I(:,e2i_Mface_u)
 idn => elem2I(:,e2i_Mface_d)
 where ( elem2I(:,e2i_elem_type) == eChannel )
    kc2 = dt * ( fQ(iup) - fQ(idn) )
    ku2 = dt * ( fQ(iup) * fUdn(iup) - fQ(idn) * fUup(idn) &
                + grav * fAdn(iup) * (fEdn(iup) - eta2)   &
                - grav * fAup(idn) * (fEup(idn) - eta2) )
 endwhere

!%  Junctions (upstream faces)
!%  HACK -- needs factor for angle with main channel
 do mm=1,upstream_face_per_elemM
    iup   => elemMI(:,eMi_MfaceUp(mm))
    where ( (elemMI(:,eMi_elem_type) == eJunctionChannel).and. &
            (elemMI(:,eMi_nfaces_u) >= mm) )
        kcM = kcM + dt * fQ(iup)
        kuM = kuM + dt * ( fQ(iup) * fUdn(iup) &
                          + grav * fAdn(iup) * (fEdn(iup) - etaM) )
    endwhere
 enddo

!%  Junctions (downstream faces)
 do mm=1,dnstream_face_per_elemM
    idn   => elemMI(:,eMi_MfaceDn(mm))
    where ( (elemMI(:,eMi_elem_type) == eJunctionChannel) .and. &
            (elemMI(:,eMi_nfaces_d) >= mm) )
        kcM = kcM - dt * fQ(idn)
        kuM = kuM - dt * ( fQ(idn) * fUdn(idn) &
                          + grav * fAdn(idn) * (fEdn(idn) - etaM) )
    endwhere
 enddo

!%  UPDATES -----------------------------------------------

!%  Note the velocity2new is actually velocity*volume at this point
 where (elem2I(:,e2i_elem_type) == eChannel)
    volume2new   =  volume2old                + thiscoef * kc2
    velocity2new = (volume2old * velocity2old + thiscoef * ku2)  &
                   / (oneR + thiscoef * dt * grav *  (mn2**2) * velocity2old / (rh2**(4.0/3.0)) )
 endwhere


 where (elemMI(:,eMi_elem_type) == eJunctionChannel)
    volumeMnew   = volumeMold                + thiscoef * kcM
    velocityMnew = (volumeMold * velocityMold + thiscoef * kuM)  &
                   / (oneR + thiscoef * dt * grav *  (mnM**2) * velocityMold / (rhM**(4.0/3.0)) )
 endwhere

!%  CORRECTIONS ----------------------------------------------------------

!%  remove negative volumes to prevent problems in velocity computation
 call adjust_negative_volume_reset (volume2new)
 call adjust_negative_volume_reset (volumeMnew)

!%  VELOCITY - divide out the volume to get the ae2r_Tempctual velocity
 where (elem2I(:,e2i_elem_type) == eChannel)
    velocity2new = velocity2new / volume2new
 endwhere

 where (elemMI(:,eMi_elem_type) == eJunctionChannel)
    velocityMnew = velocityMnew / volumeMnew
 endwhere

! release temporary arrays
 kc2 = nullvalueR
 ku2 = nullvalueR
 ones2r = nullvalueR
 kcM = nullvalueR
 kuM = nullvalueR
 onesMr = nullvalueR
 nullify(kc2, ku2, ones2r, kcM, kuM, onesMr)
 next_e2r_temparray = next_e2r_temparray - 3
 next_eMr_temparray = next_eMr_temparray - 3

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine sve_rk2_step
!
!==========================================================================
!==========================================================================
!
 subroutine rk2_update_auxiliary_variables &
    (e2r_Velocity_new, eMr_Velocity_new, e2r_Volume_new, eMr_Volume_new, &
     elem2R, elem2I, elem2YN, &
     elemMR, elemMI, elemMYN, &
     faceR,  faceI,  faceYN,  &
     bcdataDn, bcdataUp, steptime, rkiteration)

 character(64) :: subroutine_name = 'rk2_update_auxiliary_variables'

 real,      target, intent(in out)  :: elemMR(:,:)
 real,              intent(in out)  :: elem2R(:,:), faceR(:,:)
 integer,           intent(in out)  :: elem2I(:,:), elemMI(:,:), faceI(:,:)
 logical,           intent(in out)  :: elem2YN(:,:),elemMYN(:,:),faceYN(:,:)
 type(bcType),      intent(in out)  :: bcdataDn(:), bcdataUp(:)
 integer,           intent(in)      :: e2r_Velocity_new, eMr_Velocity_new
 integer,           intent(in)      :: e2r_Volume_new,   eMr_Volume_new
 integer,           intent(in)      :: rkiteration
 real,              intent(in)      :: steptime

!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name  
    
!%  advance all geometry and dynamics
 call element_geometry_update &
    (elem2R, elem2I, elem2YN, e2r_Volume_new, &
     elemMR, elemMI, elemMYN, eMr_Volume_new, &
     faceR, faceI, bcdataDn, bcdataUp, steptime, 1)    
    
!%  at this point, the channels and the junction main sections have the correct
!%  geometry, but the junction branches have provisional geometry that is
!%  a functoin of the old face free surface elevation
 call element_dynamics_update &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     bcdataDn, bcdataUp, e2r_Velocity_new, eMr_Velocity_new,  &
     e2r_Volume_new, eMr_Volume_new, steptime)
        
!%  Updating the face values by interpolation from neighbor elements
!%  This uses the estimated values from the branches
 call face_update &
    (elem2R, elem2I, elemMR, faceR, faceI, faceYN, &
     bcdataDn, bcdataUp, e2r_Velocity_new, eMr_Velocity_new, &
     e2r_Volume_new, eMr_Volume_new, steptime, rkiteration)
    
!% fix the junction branches by interp with face values
 call element_geometry_branch_fix (elemMR, elemMI, faceR, faceI )
    
!%  Ad hoc adjustment for V-shaped flowrates across a channel element
 if (setting%Method%AdjustVshapedFlowrate%Apply) then
    call adjust_Vshaped_flowrate (elem2R, faceR, elem2I, elem2YN)
 endif

 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine rk2_update_auxiliary_variables
!
!==========================================================================
!==========================================================================
!
 subroutine overwrite_old_values &
    (elemR, elemI, er_Velocity, er_Velocity_new, er_Volume, er_Volume_new, &
     ei_elem_type, ThisElemType, overwriteGhost)
!
! overwrite a new velocity data into the old data space.
!
 character(64) :: subroutine_name = 'overwrite_old_values'

 real,      intent(in out)  :: elemR(:,:)

 integer,   intent(in)      :: elemI(:,:)

 integer,   intent(in)  :: er_Velocity, er_Velocity_new
 integer,   intent(in)  :: er_Volume, er_Volume_new
 integer,   intent(in)  :: ei_elem_type, ThisElemType
 logical,   intent(in)  :: overwriteGhost

!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

 where (elemI(:,ei_elem_type) == ThisElemType)
    elemR(:,er_Velocity) = elemR(:,er_Velocity_new)
    elemR(:,er_Volume)   = elemR(:,er_Volume_new)
 endwhere
 
 if (overwriteGhost) then
    where ( (elemI(:,ei_elem_type) == eBCup) .or. &
            (elemI(:,ei_elem_type) == eBCdn) )
        elemR(:,er_Velocity) = elemR(:,er_Velocity_new)
        elemR(:,er_Volume)   = elemR(:,er_Volume_new)
    endwhere        
 endif
 
! print *
! print *,trim(subroutine_name)
! print *, elemR(:,er_Velocity)
! print *, elemR(:,er_Velocity_new)
! print *
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine overwrite_old_values
!
!==========================================================================
! END OF MODULE runge_kutta
!==========================================================================
 end module runge_kutta