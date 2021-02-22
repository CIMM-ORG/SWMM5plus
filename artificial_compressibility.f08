! module artificial_compressibility
!
! Artificial Compressibility simulator for open/surcharged pipe flow.
! Adapted from PipeAC 2020
!
!==========================================================================
!
module artificial_compressibility

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
    use storage
    use weir
    use orifice
    use xsect_tables

    implicit none

    private

    public :: ac_rk2_step

    integer :: debuglevel = 0

contains

    !==========================================================================
    !==========================================================================
    !
    subroutine ac_rk2_step &
        (e2r_Volume_old, e2r_Velocity_old, e2r_Eta_old, eMr_Volume_old,  &
        eMr_Velocity_old, eMr_Eta_old, e2r_Volume_new, e2r_Velocity_new, &
        e2r_Eta_new, eMr_Volume_new, eMr_Velocity_new, eMr_Eta_new,      &
        elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, dt, af, &
        wrk)
        !
        ! A complete RK2 step of the AC method
        !
        character(64) :: subroutine_name = 'ac_rk2_step'

        ! indexes for old/new volume and velocity storage
        integer,   intent(in) :: e2r_Volume_old, e2r_Velocity_old, e2r_Eta_old 
        integer,   intent(in) :: eMr_Volume_old, eMr_Velocity_old, eMr_Eta_old
        integer,   intent(in) :: e2r_Volume_new, e2r_Velocity_new, e2r_Eta_new
        integer,   intent(in) :: eMr_Volume_new, eMr_Velocity_new, eMr_Eta_new

        real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real,      target, intent(in out)  :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in out)  :: elem2YN(:,:), elemMYN(:,:)
        real,              intent(in)      :: dt, wrk, af(:)
        
        real,  pointer ::  volume2old(:), volume2new(:), velocity2old(:), velocity2new(:)
        real,  pointer ::  volumeMold(:), volumeMnew(:), velocityMold(:), velocityMnew(:)
        real,  pointer ::  eta2old(:), eta2new(:), etaMold(:), etaMnew(:)
        real,  pointer ::  volume2n0(:), volume2n1(:), flowrate2n0(:), flowrate2n1(:)
        real,  pointer ::  VolumeMn0(:), VolumeMn1(:), flowrateMn0(:), flowrateMn1(:)
        real,  pointer ::  dHdA2(:), zcrown2(:), zbottom2(:), rh2(:), mn2(:)
        real,  pointer ::  fullVolume2(:), fullDepth2(:), fullVolumeM(:), dHdAM(:), rhM(:)
        real,  pointer ::  length2(:), elN2(:) , breadth2(:), lengthM(:), elNM(:), mnM(:)
        real,  pointer ::  fQ(:), fUdn(:), fUup(:), fAdn(:), fAup(:), fEdn(:), fEup(:)
        real,  pointer ::  kc2(:), ku2(:), kcM(:), kuM(:)

        integer,  pointer ::  iup(:), idn(:)
        logical,  pointer ::  isFull(:), fullPipeOpen(:)
        logical,  pointer ::  maskChannelPipeAC(:), maskJunctionAC(:)

        real              :: dtau, rc2
        integer :: mm

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 

        !% F^2 from derivation
        rc2  = setting%DefaultAC%Celerity%RC ** twoR        
        dtau = dt * setting%DefaultAC%dtauFactor%dtdtau

        !%  pointers for old/new volume, velocity and eta storage
        volume2old   => elem2R(:,e2r_Volume_old)
        volume2new   => elem2R(:,e2r_Volume_new)
        velocity2old => elem2R(:,e2r_Velocity_old)
        velocity2new => elem2R(:,e2r_Velocity_new)
        eta2old      => elem2R(:,e2r_Eta_old)
        eta2new      => elem2R(:,e2r_Eta_new)

        volumeMold   => elem2R(:,eMr_Volume_old)
        volumeMnew   => elem2R(:,eMr_Volume_new)
        velocityMold => elem2R(:,eMr_Velocity_old)
        velocityMnew => elem2R(:,eMr_Velocity_new)
        etaMold      => elem2R(:,eMr_Eta_old)
        etaMnew      => elem2R(:,eMr_Eta_new)

        !%  pointer required for volume, velocity and eta calculation of elem2
        volume2n0   => elem2R(:,e2r_Volume_N0)
        volume2n1   => elem2R(:,e2r_Volume_N1)
        flowrate2n0 => elem2R(:,e2r_Flowrate_N0)
        flowrate2n1 => elem2R(:,e2r_Flowrate_N1)
        rh2         => elem2R(:,e2r_HydRadius)
        mn2         => elem2R(:,e2r_Roughness)
        length2     => elem2R(:,e2r_Length)
        elN2        => elem2R(:,e2r_elN)
        zcrown2     => elem2R(:,e2r_Zcrown)
        zbottom2    => elem2R(:,e2r_Zbottom)
        dHdA2       => elem2R(:,e2r_dHdA)
        fullVolume2 => elem2R(:,e2r_FullVolume)
        fullDepth2  => elem2R(:,e2r_FullDepth)
        breadth2    => elem2R(:,e2r_BreadthScale)

        !%  pointer required for volume, velocity and eta calculation of elemM
        volumeMn0   => elemMR(:,eMr_Volume_N0)
        volumeMn1   => elemMR(:,eMr_Volume_N1)
        flowrateMn0 => elemMR(:,eMr_Flowrate_N0)
        flowrateMn1 => elemMR(:,eMr_Flowrate_N1)
        rhM         => elemMR(:,eMr_HydRadius)
        mnM         => elemMR(:,eMr_Roughness)
        lengthM     => elemMR(:,eMr_Length)
        elNM        => elemMR(:,eMr_elN)
        dHdAM       => elemMR(:,eMr_dHdA)
        fullVolumeM => elem2R(:,eMr_FullVolume)

        !%  pointers for convenience in notation
        fUdn   => faceR(:,fr_Velocity_d) 
        fUup   => faceR(:,fr_Velocity_u) 
        fAdn   => faceR(:,fr_Area_d)
        fAup   => faceR(:,fr_Area_u)
        fEdn   => faceR(:,fr_Eta_d)
        fEup   => faceR(:,fr_Eta_u)
        fQ     => faceR(:,fr_Flowrate)

        !%  pointers for U/S an D/S face indexes
        iup    => elem2I(:,e2i_Mface_u)
        idn    => elem2I(:,e2i_Mface_d)

        !%  pointer for identifying full pipe
        isFull => elem2YN(:,e2YN_IsSurcharged)

        !%  temporary space for pipe elements
        kc2 => elem2R(:,e2r_Temp(next_e2r_temparray))
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        ku2 => elem2R(:,e2r_Temp(next_e2r_temparray))
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        !%  temporary space for juctions-pipe elements
        kcM => elemMR(:,eMr_Temp(next_eMr_temparray))
        next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

        kuM => elemMR(:,eMr_Temp(next_eMr_temparray))
        next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

        !%  temporary pointer for channel/pipe mask solved by AC
        maskChannelPipeAC   => elem2YN(:,e2YN_Temp(next_e2YN_temparray))
        next_e2YN_temparray = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

        !%  temporary mask to find the full pipes that become open
        fullPipeOpen        => elem2YN(:,e2YN_Temp(next_e2YN_temparray))
        next_e2YN_temparray = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

        !%  temporary pointer for JunctionChannel/JunctionPipe mask solved by AC
        maskJunctionAC      => elemMYN(:,eMYN_Temp(next_eMYN_temparray))
        next_eMYN_temparray = utility_advance_temp_array (next_eMYN_temparray,eMYN_n_temp)


        !%  assign nullvalues to temporary pointers
        kc2 = nullvalueR
        ku2 = nullvalueR
        kcM = nullvalueR
        kuM = nullvalueR
        maskChannelPipeAC = nullvalueL
        maskJunctionAC    = nullvalueL
        fullPipeOpen      = nullvalueL

        !%  find mask for AC solvers
        maskChannelPipeAC = ( ( (elem2I(:,e2i_elem_type) == eChannel) .or.  &
                                (elem2I(:,e2i_elem_type) == epipe) )  .and. &
                                (elem2I(:,e2i_solver) == AC) ) 

        maskJunctionAC    = ( ( (elemMI(:,eMi_elem_type) == eJunctionChannel) .or.  &
                                (elemMI(:,eMi_elem_type) == eJunctionPipe) )  .and. &
                                (elemMI(:,eMi_solver)    == AC ) )

        if ( count(maskJunctionAC) > zeroI ) then
            print*, 'error: junctions are not handeled in AC solver yet'
            stop
        endif

        !% initialize the old volume/velocity/eta to new temporary column at the begining
        !% an rk step. this is necessary because C denominator in the kVolume uses 
        !% the updated volume/eta at the 2nd rk step. 
        if (wrk == onehalfR) then
            where (maskChannelPipeAC)
                volume2new   = volume2old
                velocity2new = velocity2old
                eta2new      = eta2old
            endwhere
        endif

        !%  AC RK2 volume term for elem2
        call Kvolume2AC &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,  kc2,   &
            volume2old, volume2new, eta2old, eta2new, volume2n0, volume2n1,   &
            length2, elN2, dHdA2, fQ, iup, idn, af, wrk, dt, dtau, rc2,       &
            maskChannelPipeAC, isFull)

        !%  AC RK2 velocity term for elem2
        call Kmomentum2AC &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, ku2,     &
            velocity2old, velocity2new, eta2old, eta2new, volume2old, length2, &
            flowrate2n0, flowrate2n1, mn2, rh2, fQ, fUdn, fUup, fAdn, fAup,    &
            fEdn, fEup, iup, idn, af, wrk, dt, dtau, rc2, maskChannelPipeAC,   &
            isFull)

        where ( maskChannelPipeAC .and. (isFull .eqv. .false.) )
            volume2new = volume2old + wrk * dtau * kc2
        elsewhere ( maskChannelPipeAC .and. (isFull .eqv. .true.) )
            eta2new = eta2old + wrk * dtau * kc2
        endwhere 

        ! print*, subroutine_name
        ! print*, '...................................'
        ! print*, volume2new, 'volume2new'
        
        !%  HACK: Need derivation for juction-pipe elements

        !% find the full pipes that become open to adjust negative eta2new
        fullPipeOpen = ( (elem2I(:,e2i_elem_type) == ePipe) .and. &
                         (isFull .eqv. .true. )             .and. &
                         (eta2new .lt. zcrown2)                   )
        
        !%  CORRECTIONS ----------------------------------------------------------
        !%  remove negative volumes to prevent problems in velocity computation
        call adjust_negative_volume_reset (volume2new)
        call adjust_negative_eta_reset (eta2new, zbottom2, fullPipeOpen)

        !%  All the full pipe basic geometry handling should be here
        call get_volume_from_eta &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, eta2new, &
            volume2new, length2, fullVolume2, fullDepth2, breadth2, zbottom2,  &
            zcrown2, eta2old, isFull)

        !%  VELOCITY - divide out the volume to get the e2r_Tempctual velocity
        where (maskChannelPipeAC)
            velocity2new = velocity2new / volume2new
        endwhere

        ! print*
        ! print*, velocity2new, 'velocity2new'
        ! print*, '...... ............... ...........'

        ! release temporary arrays
        kc2 = nullvalueR
        ku2 = nullvalueR
        kcM = nullvalueR
        kuM = nullvalueR
        maskChannelPipeAC = nullvalueL
        maskJunctionAC    = nullvalueL
        fullPipeOpen      = nullvalueL

        nullify(kc2, ku2, kcM, kuM, maskChannelPipeAC, maskJunctionAC, &
                fullPipeOpen)

        next_e2r_temparray = next_e2r_temparray - 2
        next_eMr_temparray = next_eMr_temparray - 2
        next_e2YN_temparray = next_e2YN_temparray - 2
        next_eMYN_temparray = next_eMYN_temparray - 1

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine ac_rk2_step
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine Kvolume2AC &
        (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,  kc2,   &
        volume2old, volume2new, eta2old, eta2new, volume2n0, volume2n1,   &
        length2, elN2, dHdA2, fQ, iup, idn, af, wrk, dt, dtau, rc2,       &
        maskChannelPipeAC, isFull)
        !%
        !%  The RHS of continuity for an AC RK2 step looped over the domain
        !%
        character(64) :: subroutine_name = 'Kvolume2AC'

        real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real,      target, intent(in out)  :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in out)  :: elem2YN(:,:), elemMYN(:,:)

        real,       intent(inout)  :: kc2(:)
        real,       intent(in)     :: volume2new(:), eta2new(:), volume2old(:)
        real,       intent(in)     :: eta2old(:), volume2n0(:), volume2n1(:)
        real,       intent(in)     :: length2(:), elN2(:),dHdA2(:), af(:), fQ(:)
        real,       intent(in)     :: dt, dtau, rc2, wrk
        integer,    intent(in)     :: iup(:), idn(:)
        logical,    intent(in)     :: maskChannelPipeAC(:), isFull(:)

        real                       :: invdt, gammaV, gammaH, lambdaV, lambdaH
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        invdt   = oneR / dt
        gammaV  = - af(1) * invdt
        gammaH  = zeroR 
        lambdaV = oneR / (oneR - wrk * dtau * gammaV)
        lambdaH = oneR

        !%  continuity source calculation
        where ( maskChannelPipeAC )
            !%  baseline continuity source 
            kc2 = fQ(iup) - fQ(idn)
        endwhere
        ! print*, subroutine_name
        ! print*, 'kc2 at continuity source calculation'
        ! print*, kc2
        ! print*

        !%  additional source terms
        where ( maskChannelPipeAC .and. (isFull .eqv. .false.) )
            !%  additional terms for open pipe continuity source
            kc2 = kc2 - invdt * af(2) * volume2n0 - invdt * af(3) * volume2n1
            !%  combine interior gamma and lambda
            kc2 = lambdaV * (kc2 + gammaV * volume2old)
        elsewhere ( maskChannelPipeAC .and. (isFull .eqv. .true.) )
            !%  combine interior gamma and lambda (irrelevant in H)
            kc2 = lambdaH * (kc2 + gammaH * eta2old)
        endwhere
        ! print*, 'kc2 at additional source terms'
        ! print*,kc2   
        ! print*

        !%  C term numerator calculation
        where ( maskChannelPipeAC )
            !%  baseline C term numerator
            kc2 = kc2 * elN2 * rc2
        endwhere

        ! print*, 'kc2 at baseline C term numerator'
        ! print*, kc2
        ! print*

        !%  additional C term numerator calculation
        where ( maskChannelPipeAC .and. (isFull .eqv. .false.) )
            !%  C denominator for open pipe (G)
            kc2 = kc2 / ((volume2new / length2) * dHdA2 + eta2new)
        elsewhere ( maskChannelPipeAC .and. (isFull .eqv. .true.) )
            !%  C denominator for closed pipe (V) 
            kc2 = kc2 / volume2new
        endwhere
        ! print*, 'kc2 at additional C term numerator calculation'
        ! print*, kc2
        ! print*

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine Kvolume2AC
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine Kmomentum2AC &
        (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, ku2,     &
        velocity2old, velocity2new, eta2old, eta2new, volume2old, length2, &
        flowrate2n0, flowrate2n1, mn2, rh2, fQ, fUdn, fUup, fAdn, fAup,    &
        fEdn, fEup, iup, idn, af, wrk, dt, dtau, rc2, maskChannelPipeAC,   &
        isFull)
        !
        !%  Momentum with baseline T(0,0) approach to pressure term and both the
        !%  d/dtau and friction handled with time m+1 stencil.
        !
        character(64) :: subroutine_name = 'Kmomentum2AC'

        real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real,      target, intent(in out)  :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in out)  :: elem2YN(:,:), elemMYN(:,:)

        real,       intent(inout) :: ku2(:), velocity2new(:)
        real,       intent(in)    :: velocity2old(:), eta2new(:),eta2old(:)
        real,       intent(in)    :: volume2old(:), length2(:), flowrate2n0(:)
        real,       intent(in)    :: flowrate2n1(:), mn2(:), rh2(:), fQ(:)
        real,       intent(in)    :: fUdn(:), fUup(:), fAdn(:), fAup(:), fEdn(:)
        real,       intent(in)    :: fEup(:), af(:)
        real,       intent(in)    :: wrk, dt, dtau, rc2
        integer,    intent(in)    :: iup(:), idn(:)
        logical,    intent(in)    :: maskChannelPipeAC(:), isFull(:)

        real,       pointer    :: gammaUV(:)     
        real                   :: invdt, tDelta
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !%  1/dt
        invdt = oneR / dt

        !%  temporary space for gammaQ elements
        gammaUV => elem2R(:,e2r_Temp(next_e2r_temparray))
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        !%  tDelta setection
        if (setting%DefaultAC%Tsource == 'T00') then
            !%  T(0,0) method
            tDelta = zeroR
        elseif (setting%DefaultAC%Tsource == 'T10') then
            !%  T(1,0) method
            tDelta = onehalfR
        elseif (setting%DefaultAC%Tsource == 'T20') then
            !%  T(2,0) method
            tDelta = oneR / sixR
        else
            print*, 'error, unknown value for setting%DefaultAC%Tsource of '
            print*,  setting%DefaultAC%Tsource
            stop
        endif

        !%  baseline source term calculation
        where (maskChannelPipeAC)
            ku2 = fQ(iup) * fUdn(iup) - fQ(idn) * fUup(idn) + grav * fAdn(iup) * fEdn(iup) * &
                    (oneR - tDelta) - grav * fAup(idn) * fEup(idn) * (oneR - tDelta)
        endwhere
        ! print*, 'ku2 at baseline source term calculation'
        ! print*, ku2
        ! print*

        !%  Tsource term calculation
        select case (setting%DefaultAC%Tsource)
            case ('T00')
                where (maskChannelPipeAC)
                    ku2 = ku2 + grav * (fAup(idn) - fAdn(iup)) * eta2new
                endwhere
            case ('T10')
                where (maskChannelPipeAC)
                    ku2 = ku2 + grav * tDelta * (fAup(idn) * fEdn(iup) - fAdn(iup) * fEup(idn)) 
                endwhere
            case ('T20')
                where (maskChannelPipeAC)
                    ku2 = ku2 + grav * tDelta * (fAup(idn) * (fEdn(iup) + fourR * eta2new) - fAdn(iup) &
                        * (fEup(idn) + fourR * eta2new))
                endwhere
            case default
                print*, 'error, unknown value for setting%DefaultAC%Tsource of '
                print*,  setting%DefaultAC%Tsource
                stop
        end select  
        ! print*, 'ku2 at Tsource term calculation'
        ! print*, ku2
        ! print*    

        !%  Other source term calculation
        where (maskChannelPipeAC)
            !%  adding real time levels to source
            ku2 = ku2 - (af(2) * flowrate2n0 + af(3) * flowrate2n1) * invdt * length2
            !%  gamma term
            gammaUV = - af(1) * invdt - grav * ((mn2**twoR) / (rh2**(fourR/threeR))) * abs(velocity2new)
            !%  adding gamma to source
            ku2 = ku2 + gammaUV * velocity2old * volume2old
            !%  multiplying by lambda (Note that C = 1 for Q)
            ku2 = ku2 / (oneR - wrk * dtau * gammaUV)
            !%  updated U (here U = velocity2new*volume2new)
            velocity2new = volume2old * velocity2old  + wrk * dtau * ku2
        endwhere
        ! print*, 'ku2 at Other source term calculation'
        ! print*, ku2
        ! print*
       
        !% release temporary arrays
        gammaUV = nullvalueR
        nullify(gammaUV)
        next_e2r_temparray = next_e2r_temparray - 1

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine Kmomentum2AC
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine get_volume_from_eta &
        (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, eta2new, &
        volume2new, length2, fullVolume2, fullDepth2, breadth2, zbottom2,  &
        zcrown2, eta2old, isFull)
        !
        !% Find the new volume and surcharge status if solved for eta
        !
        character(64) :: subroutine_name = 'get_volume_from_eta'

        real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real,      target, intent(in out)  :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in out)  :: elem2YN(:,:), elemMYN(:,:)

        real,       intent(inout) :: eta2new(:), volume2new(:), eta2old(:)
        real,       intent(in)    :: length2(:), fullVolume2(:), fullDepth2(:)
        real,       intent(in)    :: breadth2(:), zbottom2(:), zcrown2(:)
        logical,    intent(inout) :: isFull(:)

        real,       pointer       :: YoverYfull(:)
        logical,    pointer       :: OpenPipe(:), maskarray(:)
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name
        
        !% temporary array to calculate normalized depths for special geometry types
        YoverYfull => elem2R(:,e2r_Temp(next_e2r_temparray))
        next_e2r_temparray  = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)
        !% temporary mask array for detecting open pipes
        OpenPipe   => elem2YN(:,e2YN_Temp(next_e2YN_temparray))
        next_e2YN_temparray = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)
        !% temporary mask array for multiple purpose
        maskarray  => elem2YN(:,e2YN_Temp(next_e2YN_temparray))
        next_e2YN_temparray = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

        YoverYfull = nullvalueR
        OpenPipe   = nullvalueL
        maskarray  = nullvalueL

        !===========================================================================
        !% OPEN PIPES TRANSITION TO FULL
        !% Detect transition from open to full pipe
        where ( (elem2I(:,e2i_elem_type) == ePipe) .and. &
                (isfull .eqv. .false.)             .and. &
                (volume2new .GE. fullVolume2)            )
            eta2new    = zcrown2 + (volume2new - fullVolume2) / (breadth2 * length2)
            volume2new = fullVolume2
            isfull     = .true.
        endwhere
        !===========================================================================
        !% FULL PIPES
        !% Set the full pipe volume2new
        !% These cells already have eta2new directly updated from the time-stepping.
        where ( (elem2I(:,e2i_elem_type) == ePipe) .and. (isfull) )
            volume2new = fullVolume2
        endwhere
        !===========================================================================
        !% FULL PIPES TRANSITION TO OPEN 
        !% These have eta2new and need volume2new computed
        !% Note that these are not re-designated as open until after all
        !% the eta and volume computations are complete
        !% Detect full pipe that have become open
        OpenPipe = ( (elem2I(:,e2i_elem_type) == ePipe) .and. (isfull) .and. &
                     (eta2new .LT. zcrown2) )

        !% Open Rectangular Pipe.....................................................
        where ( OpenPipe .and. (elem2I(:,e2i_Geometry) == eRectangular) )
            volume2new = (eta2new - zbottom2) * breadth2 * length2
            isfull = .false.
        endwhere

        !% Open Circular Pipe........................................................
        maskarray = ( OpenPipe .and. (elem2I(:,e2i_Geometry) == eCircular) )
        where (maskarray)
            YoverYfull = (eta2new - zbottom2) / fullDepth2 
        endwhere

        !% find normalized area using Y/Yfull from lookup tables
        !% normalized area (A/Afull) is saved in the volume2new column
        call table_lookup_mask &
            (elem2I, elem2R, volume2new, YoverYfull, ACirc, NACirc, maskarray, &
            e2i_Temp, next_e2i_temparray, e2i_n_temp)

        !% finally find volume2new from normalized area and set isfull to false    
        where (maskarray)
            volume2new = volume2new * fullVolume2
            isfull     = .false.
        endwhere

        !% Overwrite old eta column
        !% In the AC RK2 derivation eta2old is only needed in combining gamma and lambda
        !% in source term of surcharged pipe. However for surcharged pipe gamma = 0, thus
        !% eta2old becomes irrelevant. As a result we can overwrite eta2old in surcharged 
        !% pipe cases here.
        where ( (elem2I(:,e2i_elem_type) == ePipe) .and. (isfull) )
            eta2old = eta2new
        endwhere

        if (count(isfull) >0) then
            print*, '------------------------------------'
            print*, 'full pipe deteted at ', subroutine_name
            print*, isfull, 'isFull'
            print*
            print*, fullVolume2, 'fullVolume2'
            print*
            print*, volume2new, 'volume2new'
            print*
            print*, eta2new, 'eta2new'
            print*, 'press return to continue'
            read(*,*)
        endif
            
        !% nullify temporary array
        YoverYfull = nullvalueR
        OpenPipe   = nullvalueL
        maskarray  = nullvalueL
        nullify(YoverYfull, OpenPipe, maskarray)
        next_e2r_temparray  = next_e2r_temparray  - 1
        next_e2YN_temparray = next_e2YN_temparray - 2

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine get_volume_from_eta
    !
    !==========================================================================
    ! END OF MODULE artificial_compressibility
    !==========================================================================
end module artificial_compressibility
