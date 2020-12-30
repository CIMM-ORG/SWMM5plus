! module runge_kutta
!
! runge-kutta time advance for a single time step
!
!==========================================================================
!
module runge_kutta

    use adjustments
    use array_index
    use artificial_compressibility
    use bc
    use data_keys
    use diagnostic
    use element_geometry
    use element_dynamics
    use face_values
    use globals
    use setting_definition
    use utility
    use storage
    use weir
    use orifice

    implicit none

    private

    public :: rk2
    public :: sve_rk2_step
    public :: overwrite_old_values

    integer :: debuglevel = 0

contains
    !==========================================================================
    !==========================================================================
    !
    subroutine rk2 &
        (elem2R, elemMR, elem2I, elemMI, faceR, faceI, elem2YN, elemMYN,      &
        faceYN, bcdataDn, bcdataUp, thistime, dt, ID, numberPairs, ManningsN, &
        Length, zBottom, xDistance, Breadth, widthDepthData, cellType,        &
        cycleSelect)
        !
        ! runge-kutta time advance for a single time step
        !
        character(64) :: subroutine_name = 'rk2'

        real,      target, intent(in out) :: elem2R(:,:),  elemMR(:,:),  faceR(:,:)
        integer,   target, intent(in out) :: elem2I(:,:),  elemMI(:,:),  faceI(:,:)
        logical,   target, intent(in out) :: elem2YN(:,:), elemMYN(:,:), faceYN(:,:)
        type(bcType),      intent(in out) :: bcdataDn(:),  bcdataUp(:)
        real,              intent(in)     :: thistime, dt
        logical,           intent(in)     :: cycleSelect(:)

        integer,   pointer :: fdn(:), fup(:)

        integer :: e2r_Volume_new, e2r_Velocity_new,  eMr_Volume_new, eMr_Velocity_new
        integer :: e2r_Flowrate_new, e2r_Flowrate_tmp, e2r_Area_new, e2r_Area_tmp
        integer :: e2r_Eta_new, e2r_Eta_tmp, fr_Flowrate_net_new
        integer :: eMr_Flowrate_new, eMr_Area_new, eMr_Eta_new
        integer :: ii
        real    :: thiscoef(2), steptime, af(3)

        integer :: ilink

        integer, intent(in out)    :: ID(:)
        integer, intent(in out)    :: numberPairs(:)
        real,    intent(in out)    :: ManningsN(:)
        real,    intent(in out)    :: Length(:)
        real,    intent(in out)    :: zBottom(:)
        real,    intent(in out)    :: xDistance(:)
        real,    intent(in out)    :: Breadth(:)
        real,    intent(in out)    :: widthDepthData(:,:,:)
        type(string), intent(in out)   :: cellType(:)

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !%  set indexes for temporary space
        e2r_Volume_new = e2r_Temp(next_e2r_temparray)
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        e2r_Velocity_new = e2r_Temp(next_e2r_temparray)
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        e2r_Flowrate_new = e2r_Temp(next_e2r_temparray)
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        e2r_Flowrate_tmp = e2r_Temp(next_e2r_temparray)
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        e2r_Area_new = e2r_Temp(next_e2r_temparray)
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        e2r_Area_tmp = e2r_Temp(next_e2r_temparray)
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        e2r_Eta_new = e2r_Temp(next_e2r_temparray)
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        e2r_Eta_tmp = e2r_Temp(next_e2r_temparray)
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        eMr_Volume_new = eMr_Temp(next_eMr_temparray)
        next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

        eMr_Velocity_new = eMr_Temp(next_eMr_temparray)
        next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

        eMr_Flowrate_new = eMr_Temp(next_eMr_temparray)
        next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

        eMr_Area_new = eMr_Temp(next_eMr_temparray)
        next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

        eMr_Eta_new = eMr_Temp(next_eMr_temparray)
        next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

        fr_Flowrate_net_new = fr_Temp(next_fr_temparray)
        next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

        !%  zero out the temporary space
        elem2R(:,e2r_Volume_new)    = zeroR
        elem2R(:,e2r_Velocity_new)  = zeroR
        elem2R(:,e2r_Flowrate_new)  = zeroR
        elem2R(:,e2r_Area_new)      = zeroR
        elem2R(:,e2r_Eta_new)       = zeroR
        elemMR(:,eMr_Volume_new)    = zeroR
        elemMR(:,eMr_Velocity_new)  = zeroR
        elemMR(:,eMr_Flowrate_new)  = zeroR
        elemMR(:,eMr_Area_new)      = zeroR
        elemMR(:,eMr_Eta_new)       = zeroR
        faceR(:,fr_Flowrate_net_new) = zeroR

        if (  count(elem2I(:,e2i_elem_type) == eChannel) &
            + count(elem2I(:,e2i_elem_type) == ePipe)    &
            + count(elem2I(:,e2i_elem_type) == eWeir)    &
            + count(elem2I(:,e2i_elem_type) == eOrifice) &
            + count(elemMI(:,eMi_elem_type) == eStorage) &
            + count(elemMI(:,eMi_elem_type) == eJunctionChannel) > zeroI) then

            !%  coefficients for the rk2 steps
            thiscoef(1) = onehalfR
            thiscoef(2) = oneR
            !%  Coefficients for the real time derivatives in AC
            if (setting%DefaultAC%TimeStencil == 'backwards3') then
                af = (/1.5, -2.0, 0.5/)
            elseif (setting%DefaultAC%TimeStencil == 'CN') then
                af = (/0.5, -0.5, 0.0/)
            else
                print*, 'error, unknown value for setting%DefaultAC%TimeStencil'
                print*, setting%DefaultAC%TimeStencil
                stop
            endif

            ! Convergence test storage for AC 
            ! checking convergence on d/dtau of eta * A and of Q.
            where (elem2I(:,e2i_solver) == AC )
                elem2R(:,e2r_CtestH1) = elem2R(:,e2r_Eta) * elem2R(:,e2r_Area)
                elem2R(:,e2r_CtestQ1) = elem2R(:,e2r_Flowrate)
                ! net flowrate on faces - this is the flux used for scalar
                ! transport (to be added later) and for volume conservation calculations
                ! fQNetNew = wrk * fQ + fQNet ! <= the mask condition is wrong, correct later
            endwhere
            ! print*, elem2I(48:53,e2i_solver), 'e2i_solver'
            !%  step through the two steps of the RK2
            do ii=1,2
                if (cycleSelect (1)) then
                    !%  cycleSelect(1) is only true when the solver is in normal time loop
                    steptime = thistime + thiscoef(ii) * dt
                else
                    !%  when cycleSelect(1) is false, the solver is in psuedo time
                    !%  so the steptime remains the same
                    steptime = thistime
                endif

                if ( (  count(elem2I(:,e2i_solver) == SVE) &
                      + count(elemMI(:,eMi_solver) == SVE)> zeroI)) then
                    !% cycleSelect is a logical which is [T T] for time loop and [F T] for 
                    !% psuedo time loop.

                    !% HACK: In the AC convergence loop, only the 2nd step of SVE RK2 is taken.
                    !% However, element center values are still needed for face reconstruction.
                    !% Since, the update algorithm is written in terms of new volume and velocity
                    !% column, the previously calculated SVE RK2 are transferred to the new value 
                    !% column for the first RK step in AC convergence loop. This ensures that the 
                    !% faces are updated properly. 

                    !% transferring old values to new volume and velocity column are very inefficient,
                    !% because the elegemt geomety and dynamics are calculated again. This indicates we
                    !% are almost taking the whole RK step. This can be bypassed by seperating the update 
                    !% routine based on SVE and AC solve. After an effective coupling between SVE and AC
                    !% solver, seperate update routine should be written.

                    !% HACK: This might not work when the solver switches. Needs further investigation.
                    
                    if (cycleSelect(ii)) then
                        call sve_rk2_step &
                            (e2r_Volume, e2r_Velocity, eMr_Volume, eMr_Velocity, &
                            e2r_Volume_new, e2r_Velocity_new, eMr_Volume_new, eMr_Velocity_new, &
                            elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
                            thiscoef(ii))
                    else
                        where (elem2I(:,e2i_solver) == SVE) 
                            elem2R(:,e2r_Volume_new)   = elem2R(:,e2r_Volume) 
                            elem2R(:,e2r_Velocity_new) = elem2R(:,e2r_Velocity)

                        elsewhere (elemMI(:,eMi_solver) == SVE) 
                            elemMR(:,eMr_Volume_new)    = elemMR(:,eMr_Volume)
                            elemMR(:,eMr_Velocity_new)  = elemMR(:,eMr_Velocity)
                        endwhere
                    endif                      
                endif

                if (  count(elem2I(:,e2i_solver) == AC) &
                    + count(elemMI(:,eMi_solver) == AC)> zeroI ) then

                    call ac_rk2_step &
                        (e2r_Flowrate, e2r_Area, e2r_Eta, eMr_Flowrate, eMr_Area, eMr_Eta, &
                        fr_Flowrate_net, e2r_Flowrate_new, e2r_Area_new, e2r_Eta_new,      &
                        e2r_Flowrate_tmp, e2r_Area_tmp, e2r_Eta_tmp, eMr_Flowrate_new,     &
                        eMr_Area_new, eMr_Eta_new, fr_Flowrate_net_new, elem2R, elemMR,    &
                        faceR, elem2I, elemMI, elem2YN, elemMYN, dt, af, thiscoef(ii)) 
                endif

                ! !%  Sets the Qonly element geometry to provisional values
                ! call QonlyElement_provisional_geometry &
                !     (elem2R, elemMR, faceR, elem2I, elemMI)

                ! if ( count(elemMI(:,eMi_elem_type) == eStorage) > zeroI) then
                !     ! call storage step if storage unit exists in the network
                !     call storage_step &
                !         (eMr_Volume, eMr_Velocity, eMr_Volume_new, eMr_Velocity_new,  &
                !         elemMR, faceR, elemMI, elemMYN, thiscoef(ii))
                ! endif

                call rk2_update_auxiliary_variables &
                    (e2r_Velocity_new, eMr_Velocity_new, e2r_Volume_new, eMr_Volume_new,  &
                    e2r_Flowrate_new, eMr_Flowrate_new, elem2R, elem2I, elem2YN, elemMR,  &
                    elemMI, elemMYN, faceR,  faceI, faceYN, bcdataDn, bcdataUp, steptime, &
                    ii, ID, numberPairs, ManningsN, Length, zBottom, xDistance, Breadth,  &
                    widthDepthData, cellType)

                ! !% advane Qonly elemnt
                ! call QonlyElement_step &
                !     (e2r_Volume, e2r_Velocity, eMr_Volume, eMr_Velocity, e2r_Volume_new, &
                !     e2r_Velocity_new, eMr_Volume_new, eMr_Velocity_new, elem2R, elemMR,  &
                !     faceI, faceR, faceYN, elem2I, elemMI, elem2YN, elemMYN, thiscoef(ii))

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

            call overwrite_old_values &
                (elem2R, elem2I, e2r_Velocity, e2r_Velocity_new, &
                e2r_Volume, e2r_Volume_new, e2i_elem_type, ePipe, .true.)

            call overwrite_old_values &
                (elem2R, elem2I, e2r_Velocity, e2r_Velocity_new, &
                e2r_Volume, e2r_Volume_new, e2i_elem_type, eWeir, .true.)

            call overwrite_old_values &
                (elem2R, elem2I, e2r_Velocity, e2r_Velocity_new, &
                e2r_Volume, e2r_Volume_new, e2i_elem_type, eOrifice, .true.)
                
            call overwrite_old_values &
                (elemMR, elemMI, eMr_Velocity, eMr_Velocity_new, &
                eMr_Volume, eMr_Volume_new, eMi_elem_type, eStorage, .false.)

            !%  assign the solver for the next RK step depending on area
            ! call assign_solver &
            !     (elem2I, elem2R, e2r_Area, e2r_FullArea, e2i_elem_type, ePipe,  &
            !     e2i_solver, e2r_Temp, e2r_n_temp, next_e2r_temparray)

            ! !%  HACK: AC solver for Junction Pipe has not derived yet
            ! call assign_solver &
            !     (elemMI, elemMR, eMr_Area, eMr_FullArea, eMi_elem_type, eJunctionPipe, &
            !     eMi_solver, eMr_Temp, eMr_n_temp, next_eMr_temparray)

            !%  Compute the dimensional change of the AC
            where (elem2I(:,e2i_solver) == AC)
                elem2R(:,e2r_CtestH1) = elem2R(:,e2r_Eta) * elem2R(:,e2r_Area) - &
                                                elem2R(:,e2r_CtestH1)
                elem2R(:,e2r_CtestQ1) = elem2R(:,e2r_Flowrate) - elem2R(:,e2r_CtestQ1)
            endwhere
        endif

        !%  reset temporary data space
        elem2R(:,e2r_Volume_new)    = nullvalueR
        elem2R(:,e2r_Velocity_new)  = nullvalueR
        elem2R(:,e2r_Flowrate_new)  = nullvalueR
        elem2R(:,e2r_Flowrate_tmp)  = nullvalueR
        elem2R(:,e2r_Area_new)      = nullvalueR
        elem2R(:,e2r_Area_tmp)      = nullvalueR
        elem2R(:,e2r_Eta_new)       = nullvalueR
        elem2R(:,e2r_Eta_tmp)       = nullvalueR
        elemMR(:,eMr_Volume_new)    = nullvalueR
        elemMR(:,eMr_Velocity_new)  = nullvalueR
        elemMR(:,eMr_Flowrate_new)  = nullvalueR
        elemMR(:,eMr_Area_new)      = nullvalueR
        elemMR(:,eMr_Eta_new)       = nullvalueR
        faceR(:,fr_Flowrate_net_new) = nullvalueR

        next_e2r_temparray = next_e2r_temparray - 8
        next_eMr_temparray = next_eMr_temparray - 5
        next_fr_temparray  = next_fr_temparray  - 1

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
        logical,   target, intent(in out)  :: elem2YN(:,:), elemMYN(:,:)
        real,              intent(in)      :: thiscoef

        real,  pointer ::  fQ(:), fUdn(:), fUup(:), fAdn(:), fAup(:), fEdn(:), fEup(:)
        real,  pointer ::  kc2(:), ku2(:), ones2r(:)
        real,  pointer ::  kcM(:), kuM(:), onesMr(:)
        real,  pointer ::  volume2old(:), volume2new(:), velocity2old(:), velocity2new(:)
        real,  pointer ::  volumeMold(:), volumeMnew(:), velocityMold(:), velocityMnew(:)
        real,  pointer ::  eta2(:), etaM(:), rh2(:), mn2(:), rhM(:), mnM(:)

        logical, pointer :: maskChannelPipeSVE(:), maskJunctionChannelPipeSVE(:)
        integer, pointer :: iup(:), idn(:)

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

        !%  temporary pointer for channel/pipe mask solved by SVE
        maskChannelPipeSVE => elem2YN(:,e2YN_Temp(next_e2YN_temparray))
        next_e2YN_temparray = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

        !%  temporary pointer for JunctionChannel/JunctionPipe mask solved by SVE
        maskJunctionChannelPipeSVE => elemMYN(:,eMYN_Temp(next_eMYN_temparray))
        next_eMYN_temparray = utility_advance_temp_array (next_eMYN_temparray,eMYN_n_temp)

        !%  zero temporary arrays
        kc2 = zeroR
        ku2 = zeroR

        kcM = zeroR
        kuM = zeroR

        !%  null temp array
        maskChannelPipeSVE = nullvalueL
        maskJunctionChannelPipeSVE = nullvalueL

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

        !%  finding the masks
        maskChannelPipeSVE = ( ( (elem2I(:,e2i_elem_type) == eChannel) .or.  &
                                 (elem2I(:,e2i_elem_type) == epipe) )  .and. &
                                 (elem2I(:,e2i_solver) == SVE) ) 


        maskJunctionChannelPipeSVE = ( ( (elemMI(:,eMi_elem_type) == eJunctionChannel) .or.  &
                                         (elemMI(:,eMi_elem_type) == eJunctionPipe) )  .and. &
                                         (elemMI(:,eMi_solver)    == SVE ) )
        !%  SOURCE TERMS ----------------------------------

        !%  Channel elements (one upstream and one downstream face)
        !%  Using the T00 method of Hodges & Liu, 2018 for momentum

        iup => elem2I(:,e2i_Mface_u)
        idn => elem2I(:,e2i_Mface_d)
        where (maskChannelPipeSVE)
            kc2 = dt * ( fQ(iup) - fQ(idn) )
            ! ku2 = dt * ( fQ(iup) * fUdn(iup) - fQ(idn) * fUup(idn) &
            !     + grav * fAdn(iup) * (fEdn(iup) - eta2)   &
            !     - grav * fAup(idn) * (fEup(idn) - eta2) )

            !% experimental: matching with SvePy code
            ku2 = dt * ( fQ(iup) * fUdn(iup) * sign(ones2r,fQ(iup)) &
                       - fQ(idn) * fUup(idn) * sign(ones2r,fQ(idn))   &
                       + grav * fAdn(iup) * (fEdn(iup) - eta2)   &
                       - grav * fAup(idn) * (fEup(idn) - eta2) )
        endwhere

        !%  Junctions (upstream faces)
        !%  HACK -- needs factor for angle with main channel
        do mm=1,upstream_face_per_elemM
            iup   => elemMI(:,eMi_MfaceUp(mm))
            where (maskJunctionChannelPipeSVE     .and. &
                  (elemMI(:,eMi_nfaces_u) >= mm)  )
                kcM = kcM + dt * fQ(iup)
                kuM = kuM + dt * ( fQ(iup) * fUdn(iup) &
                    + grav * fAdn(iup) * (fEdn(iup) - etaM) )
            endwhere
        enddo

        !%  Junctions (downstream faces)
        do mm=1,dnstream_face_per_elemM
            idn   => elemMI(:,eMi_MfaceDn(mm))
            where (maskJunctionChannelPipeSVE     .and. &
                  (elemMI(:,eMi_nfaces_d) >= mm)  )
                kcM = kcM - dt * fQ(idn)
                kuM = kuM - dt * ( fQ(idn) * fUdn(idn) &
                    + grav * fAdn(idn) * (fEdn(idn) - etaM) )
            endwhere
        enddo

        !%  UPDATES -----------------------------------------------

        !%  Note the velocity2new is actually velocity*volume at this point
        where (maskChannelPipeSVE)
            volume2new   =  volume2old                + thiscoef * kc2
            velocity2new = (volume2old * velocity2old + thiscoef * ku2)  &
                / (oneR + thiscoef * dt * grav *  (mn2**2) * abs(velocity2old) / (rh2**(4.0/3.0)) )
        endwhere


        where (maskJunctionChannelPipeSVE) 
            volumeMnew   = volumeMold                 + thiscoef * kcM
            velocityMnew = (volumeMold * velocityMold + thiscoef * kuM)  &
                / (oneR + thiscoef * dt * grav *  (mnM**2) * abs(velocityMold) / (rhM**(4.0/3.0)) )
        endwhere

        !%  CORRECTIONS ----------------------------------------------------------

        !%  remove negative volumes to prevent problems in velocity computation
        call adjust_negative_volume_reset (volume2new)
        call adjust_negative_volume_reset (volumeMnew)

        !%  VELOCITY - divide out the volume to get the e2r_Tempctual velocity
        where (maskChannelPipeSVE)
            velocity2new = velocity2new / volume2new
        endwhere

        where (maskJunctionChannelPipeSVE)
            velocityMnew = velocityMnew / volumeMnew
        endwhere

        ! release temporary arrays
        kc2 = nullvalueR
        ku2 = nullvalueR
        ones2r = nullvalueR
        kcM = nullvalueR
        kuM = nullvalueR
        onesMr = nullvalueR
        maskChannelPipeSVE = nullvalueL
        maskJunctionChannelPipeSVE = nullvalueL

        nullify(kc2, ku2, ones2r, kcM, kuM, onesMr, maskChannelPipeSVE, &
                maskJunctionChannelPipeSVE)
        
        next_e2r_temparray = next_e2r_temparray - 3
        next_eMr_temparray = next_eMr_temparray - 3
        next_e2YN_temparray = next_e2YN_temparray - 1
        next_eMYN_temparray = next_eMYN_temparray - 1

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine sve_rk2_step
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine rk2_update_auxiliary_variables &
        (e2r_Velocity_new, eMr_Velocity_new, e2r_Volume_new, eMr_Volume_new,  &
        e2r_Flowrate_new, eMr_Flowrate_new, elem2R, elem2I, elem2YN, elemMR,  &
        elemMI, elemMYN, faceR,  faceI, faceYN, bcdataDn, bcdataUp, steptime, &
        rkiteration, ID, numberPairs, ManningsN, Length, zBottom, xDistance,  &
        Breadth, widthDepthData, cellType)

        character(64) :: subroutine_name = 'rk2_update_auxiliary_variables'

        real,      target, intent(in out)  :: elemMR(:,:)
        real,              intent(in out)  :: elem2R(:,:), faceR(:,:)
        integer,           intent(in out)  :: elem2I(:,:), elemMI(:,:), faceI(:,:)
        logical,           intent(in out)  :: elem2YN(:,:),elemMYN(:,:),faceYN(:,:)
        type(bcType),      intent(in out)  :: bcdataDn(:), bcdataUp(:)
        integer,           intent(in)      :: e2r_Velocity_new, eMr_Velocity_new
        integer,           intent(in)      :: e2r_Volume_new,   eMr_Volume_new
        integer,           intent(in)      :: e2r_Flowrate_new, eMr_Flowrate_new
        integer,           intent(in)      :: rkiteration
        real,              intent(in)      :: steptime

        integer, intent(in out)    :: ID(:)
        integer, intent(in out)    :: numberPairs(:)
        real,    intent(in out)    :: ManningsN(:)
        real,    intent(in out)    :: Length(:)
        real,    intent(in out)    :: zBottom(:)
        real,    intent(in out)    :: xDistance(:)
        real,    intent(in out)    :: Breadth(:)
        real,    intent(in out)    :: widthDepthData(:,:,:)
        type(string), intent(in out)   :: cellType(:)

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !%  advance all geometry and dynamics
        call element_geometry_update &
            (elem2R, elem2I, elem2YN, e2r_Volume_new, elemMR, elemMI, elemMYN, &
            eMr_Volume_new, faceR, faceI, bcdataDn, bcdataUp, steptime, 1, ID, &
            numberPairs, ManningsN, Length, zBottom, xDistance, Breadth,       &
            widthDepthData, cellType)

        !%  at this point, the channels and the junction main sections have the correct
        !%  geometry, but the junction branches have provisional geometry that is
        !%  a functoin of the old face free surface elevation
        call element_dynamics_update &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, bcdataDn, &
            bcdataUp, e2r_Velocity_new, eMr_Velocity_new, e2r_Volume_new,       &
            eMr_Volume_new, e2r_Flowrate_new, eMr_Flowrate_new, steptime)

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
    !==========================================================================
    !
    subroutine QonlyElement_provisional_geometry &
        (elem2R, elemMR, faceR, elem2I, elemMI)
        ! this subroutine sets the Qonly element geometry to zero.
        character(64) :: subroutine_name = 'QonlyElement_provisional_geometry'


        real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real,      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)

        integer :: mm
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


        where      ( (elem2I(:,e2i_meta_elem_type) == eQonly) )

            elem2R(:,e2r_Area)        = 1.0e-7
            elem2R(:,e2r_Eta)         = 1.0e-7
            elem2R(:,e2r_Perimeter)   = 1.0e-7
            elem2R(:,e2r_HydDepth)    = 1.0e-7
            elem2R(:,e2r_HydRadius)   = 1.0e-7
            elem2R(:,e2r_Topwidth)    = 1.0e-7
            elem2R(:,e2r_Depth)       = 1.0e-7
            elem2R(:,e2r_Volume)      = 1.0e-7
            elem2R(:,e2r_Velocity)    = 1.0e-7
        endwhere

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine QonlyElement_provisional_geometry
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine QonlyElement_step &
        (e2r_Volume_old, e2r_Velocity_old, eMr_Volume_old, eMr_Velocity_old, &
        e2r_Volume_new, e2r_Velocity_new, eMr_Volume_new, eMr_Velocity_new, &
        elem2R, elemMR, faceI, faceR, faceYN, elem2I, elemMI, elem2YN, &
        elemMYN, thiscoef)
        !
        character(64) :: subroutine_name = 'QonlyElement_step'

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
        logical,   target, intent(in out)  :: faceYN(:,:)
        real,              intent(in)      :: thiscoef

        real,      pointer  :: valueUp(:), valueDn(:)
        real,      pointer  :: weightUpQ(:), weightDnQ(:)
        real,      pointer  :: faceQ(:)
        logical,   pointer  :: facemask(:)
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        valueUp => faceR(:,fr_Temp(next_fr_temparray))
        next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

        valueDn => faceR(:,fr_Temp(next_fr_temparray))
        next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

        weightUpQ => faceR(:,fr_Temp(next_fr_temparray))
        next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

        weightDnQ => faceR(:,fr_Temp(next_fr_temparray))
        next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

        facemask   => faceYN(:,fYN_Temp(next_fYN_temparray))
        next_fYN_temparray = utility_advance_temp_array (next_fYN_temparray,fYN_n_temp)

        faceQ => faceR(:,fr_Flowrate)

        !% advance flow, geometry in Qonly elements
        if ( count(elem2I(:,e2i_elem_type) == eWeir)  > zeroI) then
            !% call weir step if weirs exist in the network
            call weir_step &
                (e2r_Volume_old, e2r_Velocity_old, eMr_Volume_old, eMr_Velocity_old, &
                e2r_Volume_new, e2r_Velocity_new, eMr_Volume_new, eMr_Velocity_new, &
                elem2R, elemMR, faceI, faceR, faceYN, elem2I, elemMI, elem2YN, &
                elemMYN, thiscoef)
        endif

        if ( count(elem2I(:,e2i_elem_type) == eOrifice)  > zeroI) then
            !% call orifice step if orifices exist in the network
            call orifice_step &
                (e2r_Volume_old, e2r_Velocity_old, eMr_Volume_old, eMr_Velocity_old, &
                e2r_Volume_new, e2r_Velocity_new, eMr_Volume_new, eMr_Velocity_new, &
                elem2R, elemMR, faceI, faceR, faceYN, elem2I, elemMI, elem2YN, &
                elemMYN, thiscoef)
        endif

        !% face reconstruction -- only flow values
        !% update the flow to their faces
        facemask = ( (faceI(:,fi_meta_etype_u) == eQonly) .or. &
            (faceI(:,fi_meta_etype_d) == eQonly) )

        weightUpQ = setting%Limiter%Timescale%Maximum
        weightDnQ = setting%Limiter%Timescale%Maximum

        where (facemask)
            weightUpQ = elem2R(faceI(:,fi_Melem_u),e2r_Timescale_Q_d)
            weightDnQ = elem2R(faceI(:,fi_Melem_d),e2r_Timescale_Q_u)
            valueUp  = elem2R(faceI(:,fi_Melem_u),e2r_Flowrate)
            valueDn  = elem2R(faceI(:,fi_Melem_d),e2r_Flowrate)
            !% linear interpolation
            faceQ = (weightUpQ * valueDn + weightDnQ * valueUp) /(weightUpQ + weightDnQ)
        endwhere

        valueUp    = nullvalueR
        valueDn    = nullvalueR
        weightUpQ  = nullvalueR
        weightDnQ  = nullvalueR

        nullify(valueUp, valueDn, weightUpQ, weightDnQ, facemask)

        next_fr_temparray  = next_fr_temparray  - 4
        next_fYN_temparray = next_fYN_temparray - 1

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine QonlyElement_step
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine assign_solver &
        (elemI, elemR, er_Area, er_FullArea, ei_elem_type, ThisElemType,  &
        ei_solver, er_Temp, er_n_temp, next_er_temparray)
        !
        character(64) :: subroutine_name = 'assign_solver'

        real,      target, intent(inout)  :: elemR(:,:)
        integer,   target, intent(inout)  :: elemI(:,:)

        integer,   intent(in)      ::  er_Area, er_FullArea, er_n_temp
        integer,   intent(in)      ::  ei_elem_type, ei_solver, ThisElemType
        integer,   intent(in)      ::  er_Temp(:)
        integer,   intent(inout)   ::  next_er_temparray

        real,      pointer  :: AoverAfull(:), Area(:), FullArea(:)
        real                :: switchBufferP, switchBufferM
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


        !%  temporary space for pipe elements
        AoverAfull => elemR(:,er_Temp(next_er_temparray))
        next_er_temparray = utility_advance_temp_array (next_er_temparray,er_n_temp)

        Area     => elemR(:,er_Area)
        FullArea => elemR(:,er_FullArea)

        !%  values +/- buffer for the solver switch
        switchBufferP = setting%DefaultAC%Switch%Area + setting%DefaultAC%Switch%Buffer
        switchBufferM = setting%DefaultAC%Switch%Area - setting%DefaultAC%Switch%Buffer

        where ( elemI(:,ei_elem_type) == ThisElemType )
            AoverAfull = Area / FullArea
        endwhere

        ! selecting appropriate solver for pipe
        where ( (elemI(:,ei_elem_type) == ThisElemType) .and. &
                (elemI(:,ei_solver) == SVE)             .and. &
                (AoverAfull .GE. switchBufferP) )

            elemI(:,ei_solver) = AC

        elsewhere( (elemI(:,ei_elem_type) == ThisElemType) .and. &
                   (elemI(:,ei_solver) == AC)              .and. &
                   (AoverAfull .LE. switchBufferM) )

            elemI(:,ei_solver) = SVE
        endwhere

        ! print*, trim(subroutine_name)
        ! print*, 'AoverAfull', AoverAfull
        ! print*, 'Selected solver', elemI(:,ei_solver)
        
        AoverAfull = nullvalueR
        nullify(AoverAfull)
        next_er_temparray = next_er_temparray - 1

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine assign_solver
    !
    !==========================================================================
    ! END OF MODULE runge_kutta
    !==========================================================================
end module runge_kutta
