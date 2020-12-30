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

    implicit none

    private

    public :: ac_rk2_step
    ! public :: ac_rk2

    integer :: debuglevel = 0

contains

    !==========================================================================
    !==========================================================================
    !
    subroutine ac_rk2_step &
        (e2r_Flowrate_old, e2r_Area_old, e2r_Eta_old, eMr_Flowrate_old,   &
        eMr_Area_old, eMr_Eta_old, fr_Flowrate_net_old, e2r_Flowrate_new, &
        e2r_Area_new, e2r_Eta_new, e2r_Flowrate_tmp, e2r_Area_tmp,        &
        e2r_Eta_tmp, eMr_Flowrate_new, eMr_Area_new, eMr_Eta_new,         &
        fr_Flowrate_net_new, elem2R, elemMR, faceR, elem2I, elemMI,       &
        elem2YN, elemMYN, dt, af, wrk)
        !
        ! A complete RK2 step of the AC method
        !
        character(64) :: subroutine_name = 'ac_rk2_step'

        ! indexes for old/new volume and velocity storage
        integer,    intent(in) :: e2r_Flowrate_old, e2r_Area_old, e2r_Eta_old
        integer,    intent(in) :: eMr_Flowrate_old, eMr_Area_old, eMr_Eta_old
        integer,    intent(in) :: e2r_Flowrate_new, e2r_Area_new, e2r_Eta_new
        integer,    intent(in) :: e2r_Flowrate_tmp, e2r_Area_tmp, e2r_Eta_tmp
        integer,    intent(in) :: eMr_Flowrate_new, eMr_Area_new, eMr_Eta_new
        integer,    intent(in) :: fr_Flowrate_net_old, fr_Flowrate_net_new

        real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real,      target, intent(in out)  :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in out)  :: elem2YN(:,:), elemMYN(:,:)
        real,              intent(in)      :: dt, wrk, af(:)
        
        real,       pointer ::  flowrate2(:), flowrate2new(:), area2(:), area2new(:)
        real,       pointer ::  flowrateM(:), flowrateMnew(:), areaM(:), areaMnew(:)
        real,       pointer ::  eta2(:), eta2new(:), etaM(:), etaMnew(:), fQNet(:), fQNetNew(:)
        real,       pointer ::  area2n0(:), area2n1(:), flowrate2n0(:), flowrate2n1(:), velocity2(:) 
        real,       pointer ::  zcrown2(:), zbottom2(:), rh2(:), mn2(:), length2(:), elN2(:), dHdA2(:)
        real,       pointer ::  areaMn0(:), areaMn1(:), flowrateMn0(:), flowrateMn1(:), velocityM(:)
        real,       pointer ::  rhM(:), mnM(:), lengthM(:), elNM(:), dhDAM(:)
        real,       pointer ::  CtestH12(:), CtestQ12(:), CtestH1M(:), CtestQ1M(:)
        real,       pointer ::  kH2(:), kQ2(:), kHM(:), kQM(:)
        real,       pointer ::  fQ(:), fUdn(:), fUup(:), fAdn(:), fAup(:), fEdn(:), fEup(:)
        real,       pointer ::  area2tmp(:), eta2tmp(:), flowrate2tmp(:)

        integer,    pointer ::  iup(:), idn(:)
        logical,    pointer ::  maskPipeAc(:), isFull(:), maskarray(:)

        real                :: dtau, rc2

        integer :: mm

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 

        rc2  = setting%DefaultAC%Celerity%RC ** 2        ! F^2 in derivation
        dtau = dt * setting%DefaultAC%dtauFactor%dtdtau

        !%  pointers for Flowrate, area and eta storage (updating)
        flowrate2    => elem2R(:,e2r_Flowrate_old)
        flowrate2new => elem2R(:,e2r_Flowrate_new)
        flowrate2tmp => elem2R(:,e2r_Flowrate_tmp)
        area2        => elem2R(:,e2r_Area_old)
        area2new     => elem2R(:,e2r_Area_new)
        area2tmp     => elem2R(:,e2r_Area_tmp)
        eta2         => elem2R(:,e2r_Eta_old)
        eta2new      => elem2R(:,e2r_Eta_new)
        eta2tmp      => elem2R(:,e2r_Eta_tmp)

        flowrateM    => elem2R(:,eMr_Flowrate_old)
        flowrateMnew => elem2R(:,eMr_Flowrate_new)
        areaM        => elem2R(:,eMr_Area_old)
        areaMnew     => elem2R(:,eMr_Area_new)
        etaM         => elem2R(:,eMr_Eta_old)
        etaMnew      => elem2R(:,eMr_Eta_new)

        fQNet        => faceR(:,fr_Flowrate_net_old)
        fQNetNew     => faceR(:,fr_Flowrate_net_new)

        !%  pointer required for Flowrate, area and eta calculation of pipe elements
        area2n0     => elem2R(:,e2r_Area_N0)
        area2n1     => elem2R(:,e2r_Area_N1)
        flowrate2n0 => elem2R(:,e2r_Flowrate_N0)
        flowrate2n1 => elem2R(:,e2r_Flowrate_N1)
        velocity2   => elem2R(:,e2r_Velocity)
        rh2         => elem2R(:,e2r_HydRadius)
        mn2         => elem2R(:,e2r_Roughness)
        length2     => elem2R(:,e2r_Length)
        elN2        => elem2R(:,e2r_elN)
        zcrown2     => elem2R(:,e2r_Zcrown)
        zbottom2    => elem2R(:,e2r_Zbottom)
        dHdA2       => elem2R(:,e2r_dHdA)
        CtestH12    => elem2R(:,e2r_CtestH1)
        CtestQ12    => elem2R(:,e2r_CtestQ1)

        !%  pointer required for Flowrate, area and eta calculation of junction-pipe elements
        areaMn0     => elemMR(:,eMr_Area_N0)
        areaMn1     => elemMR(:,eMr_Area_N1)
        flowrateMn0 => elemMR(:,eMr_Flowrate_N0)
        flowrateMn1 => elemMR(:,eMr_Flowrate_N1)
        velocityM   => elemMR(:,eMr_Velocity)
        rhM         => elemMR(:,eMr_HydRadius)
        mnM         => elemMR(:,eMr_Roughness)
        lengthM     => elemMR(:,eMr_Length)
        elNM        => elemMR(:,eMr_elN)
        dHdAM       => elemMR(:,eMr_dHdA)
        CtestH1M    => elemMR(:,eMr_CtestH1)
        CtestQ1M    => elemMR(:,eMr_CtestH1)

        !%  temporary space for pipe elements
        kH2 => elem2R(:,e2r_Temp(next_e2r_temparray))
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        kQ2 => elem2R(:,e2r_Temp(next_e2r_temparray))
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        !%  temporary space for juctions-pipe elements
        kHM => elemMR(:,eMr_Temp(next_eMr_temparray))
        next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

        kQM => elemMR(:,eMr_Temp(next_eMr_temparray))
        next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

        !%  temporary pointer for pipe mask solving by AC
        maskPipeAc  => elem2YN(:,e2YN_Temp(next_e2YN_temparray))
        next_e2YN_temparray = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

        !%  temporary pointer to find the full pipes that become open
        maskarray   => elem2YN(:,e2YN_Temp(next_e2YN_temparray))
        next_e2YN_temparray = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

        !%  assign nullvalues to temporary pointers
        kH2 = nullvalueR
        kQ2 = nullvalueR
        kHM = nullvalueR
        kQM = nullvalueR
        maskPipeAc = nullvalueL
        maskarray  = nullvalueL

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

        maskPipeAc = ( (elem2I(:,e2i_elem_type) == ePipe) .and. (elem2I(:,e2i_solver) == AC) )
        ! store the original area, eta and flowrate going into RK2 steps
        ! these original values will be used to march the variables forward
        if (wrk == onehalfR) then
            where (maskPipeAc)
                area2tmp = area2
                eta2tmp = eta2
                flowrate2tmp = flowrate2
            endwhere
        endif

        !%  AC RK4 volume (area) term
        call KvolumePipe &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, kH2, &
            length2, area2new, area2tmp, area2, area2n0, area2n1, eta2new, &
            eta2tmp, eta2, elN2, dHdA2, fQ, iup, idn, af, wrk, dt, dtau,   &
            rc2, maskPipeAc, isFull)

        !%  AC RK4 flowrate term
        call Kmomentum3Pipe &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, kQ2, eta2, &
            length2, flowrate2new, flowrate2tmp, flowrate2, flowrate2n0,         &
            flowrate2n1, mn2, rh2, velocity2, fQ, fUdn, fUup, fAdn, fAup, fEdn,  &
            fEup, iup, idn, af, wrk, dt, dtau, rc2, maskPipeAc, isFull)

        where ( maskPipeAc .and. (isFull .eqv. .false.) )
            area2new = area2tmp + wrk * dtau * kH2
        elsewhere ( maskPipeAc .and. (isFull .eqv. .true.) )
            eta2new  = eta2tmp  + wrk * dtau * kH2
        endwhere 

        !% HACK: Need derivation for juction-pipe elements

        !%  CORRECTIONS ----------------------------------------------------------
        !%  remove negative area
        call adjust_negative_area_reset (area2new)
        !% find the full pipes that become open to adjust negative eta2new
        maskarray = ( maskPipeAc .and. (isFull .eqv. .true.) .and. &
                    ( eta2new .lt. zcrown2) )
        call adjust_negative_eta_reset (eta2new, zbottom2, maskarray)

        !% HACK: new area and new eta are needed to be transferred to the
        !% old value column. Though AC solver continuity is solved for 
        !% either area or eta, the SVE continuity solves for volume, and
        !% eta and area are derived variables which is used for face interplation 
        !% afterwards. To ensure a unified face update subroutine, new area and
        !% eta are now transferred to the old column. Talk with Dr. Hodges
        !% wheater this will cause problem.

        !% overwrite old area and eta values
        where (maskPipeAc)
            area2 = area2new
            eta2  = eta2new
        endwhere

        ! release temporary arrays
        kH2 = nullvalueR
        kQ2 = nullvalueR
        kHM = nullvalueR
        kQM = nullvalueR
        maskPipeAc = nullvalueL
        maskarray  = nullvalueL
        nullify(kH2, kQ2, kHM, kQM, maskPipeAc, maskarray)
        next_e2r_temparray = next_e2r_temparray - 2
        next_eMr_temparray = next_eMr_temparray - 2
        next_e2YN_temparray = next_e2YN_temparray - 2

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine ac_rk2_step
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine KvolumePipe &
        (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, kH2, &
        length2, area2new, area2old, area2, area2n0, area2n1, eta2new, &
        eta2old, eta2, elN2, dHdA2, fQ, iup, idn, af, wrk, dt, dtau,   &
        rc2, maskPipeAc, isFull)
        !%
        !%  The RHS of continuity for an AC RK2 step looped over the domain
        !%
        character(64) :: subroutine_name = 'KvolumePipe'

        real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real,      target, intent(in out)  :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in out)  :: elem2YN(:,:), elemMYN(:,:)

        real,       intent(inout)  :: kH2(:), area2new(:), eta2new(:), eta2old(:)
        real,       intent(in)     :: area2(:), area2old(:), area2n0(:), area2n1(:)
        real,       intent(in)     :: eta2(:), length2(:), elN2(:)
        real,       intent(in)     :: dHdA2(:), af(:), fQ(:)
        real,       intent(in)     :: dt, dtau, rc2, wrk
        integer,    intent(in)     :: iup(:), idn(:)
        logical,    intent(in)     :: maskPipeAc(:), isFull(:)

        real                    :: invdt, gammaA, gammaH, lambdaA, lambdaH
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        invdt   = oneR / dt
        gammaA  = - af(1) * invdt
        gammaH  = zeroR 
        lambdaA = oneR / (oneR - dtau * gammaA)
        lambdaH = oneR

        !%  continuity source calculation
        where ( maskPipeAc )
            !%  baseline continuity source 
            kH2 = ( fQ(iup) - fQ(idn) ) / length2
        endwhere

        !%  additional source terms
        where ( maskPipeAc .and. (isFull .eqv. .false.) )
            !%  additional terms for open pipe continuity source
            kH2 = kH2 - invdt * af(2) * area2n0 - invdt * af(3) * area2n1
            !%  combine interior gamma and lambda (irrelevant in H)
            kH2 = lambdaA * (kH2 + gammaA * area2)
        elsewhere ( maskPipeAc .and. (isFull .eqv. .true.) )
            !%  combine interior gamma and lambda (irrelevant in H)
            kH2 = lambdaH * (kH2 + gammaH * eta2)
        endwhere

        !% C term numerator calculation
        where ( maskPipeAc )
            !%  baseline C term numerator
            kH2 = kH2 * elN2 * rc2
        endwhere

        ! additional C term numerator calculation
        where ( maskPipeAc .and. (isFull .eqv. .false.) )
            !%  C denominator for open pipe (G)
            kH2 = kH2 / (area2 * dHdA2 + eta2)
        elsewhere ( maskPipeAc .and. (isFull .eqv. .true.) )
            !%  C denominator for closed pipe (A) 
            kH2 = kH2 / area2
        endwhere

        ! Update area and eta from volume solution (done after momentum!)
        
        ! output print for debug
        ! print*, '**************** AC SOLVER ****************'
        ! print*, trim(subroutine_name)
        ! print*, 'KH2         ',KH2(995:1000)
        
        ! print*, 'area    ',area2
        ! print*, 'area2new',area2new
        ! print*, 'eta2new ',eta2new

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine KvolumePipe
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine Kmomentum3Pipe &
        (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, kQ2, eta2, &
        length2, flowrate2new, flowrate2old, flowrate2, flowrate2n0,         &
        flowrate2n1, mn2, rh2, velocity2, fQ, fUdn, fUup, fAdn, fAup, fEdn,  &
        fEup, iup, idn, af, wrk, dt, dtau, rc2, maskPipeAc, isFull)
        !
        !%  Momentum with baseline T(0,0) approach to pressure term and both the
        !%  d/dtau and friction handled with time m+1 stencil.
        !
        character(64) :: subroutine_name = 'Kmomentum3Pipe'

        real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real,      target, intent(in out)  :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in out)  :: elem2YN(:,:), elemMYN(:,:)

        real,       intent(inout) :: kQ2(:),  flowrate2new(:), flowrate2old(:)
        real,       intent(in)    :: flowrate2n0(:),flowrate2n1(:),flowrate2(:)
        real,       intent(in)    :: mn2(:), rh2(:), velocity2(:), length2(:)
        real,       intent(in)    :: eta2(:), fQ(:), fUdn(:), fUup(:)
        real,       intent(in)    :: fAdn(:), fAup(:), fEdn(:), fEup(:), af(:)
        real,       intent(in)    :: wrk, dt, dtau, rc2
        integer,    intent(in)    :: iup(:), idn(:)
        logical,    intent(in)    :: maskPipeAc(:), isFull(:)

        real,       pointer    :: gammaQ(:)     
        real                   :: invdt, tDelta
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        invdt = 1.0 / dt

        !%  temporary space for gammaQ elements
        gammaQ => elem2R(:,e2r_Temp(next_e2r_temparray))
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        !%  tDelta setection
        if (setting%DefaultAC%Tsource == 'T00') then
            !%  T(0,0) method
            tDelta = 0.0
        elseif (setting%DefaultAC%Tsource == 'T10') then
            !%  T(1,0) method
            tDelta = 0.5
        elseif (setting%DefaultAC%Tsource == 'T20') then
            !%  T(2,0) method
            tDelta = 1.0 / 6.0
        else
            print*, 'error, unknown value for setting%DefaultAC%Tsource of '
            print*,  setting%DefaultAC%Tsource
            stop
        endif

        !%  baseline source term calculation
        where (maskPipeAc)
            kQ2 = fQ(iup) * fUdn(iup) - fQ(idn) * fUup(idn) + grav * fAdn(iup) * fEdn(iup) * &
                    (1.0 - tDelta) - grav * fAup(idn) * fEup(idn) * (1.0 - tDelta)
        endwhere
        ! print*, 'iup         ', iup(997:1001)
        ! print*, 'idn         ', idn(997:1001)
        ! print*, 'fQ          ', fQ(997:1001)
        ! print*, 'fUdn        ', fUdn(997:1001)
        ! print*, 'fUup        ', fUup(997:1001)
        ! print*, 'fAdn        ', fAdn(997:1001)
        ! print*, 'fAup        ', fAup(997:1001)
        ! print*, 'fEdn        ', fEdn(997:1001)
        ! print*, 'fEup        ', fEup(997:1001)
        ! print*, 'KQ2 base    ', KQ2(997:1001)
        ! print*, 'eta2        ', eta2(997:1001)
        !%  Tsource term calculation
        select case (setting%DefaultAC%Tsource)
            case ('T00')
                where (maskPipeAc)
                    kQ2 = kQ2 + grav * (fAup(idn) - fAdn(iup)) * eta2
                endwhere
            case ('T10')
                where (maskPipeAc)
                    kQ2 = kQ2 + grav * tDelta * (fAup(idn) * fEdn(iup) - fAdn(iup) * fEup(idn)) 
                endwhere
            case ('T20')
                where (maskPipeAc)
                    kQ2 = kQ2 + grav * tDelta * (fAup(idn) * (fEdn(iup) + 4.0 * eta2) - fAdn(iup) &
                        * (fEup(idn) + 4.0 * eta2))
                endwhere
            case default
                print*, 'error, unknown value for setting%DefaultAC%Tsource of '
                print*,  setting%DefaultAC%Tsource
                stop
        end select
        ! print*, 'KQ2 T00      ',KQ2(997:1001)
        !%  Other source term calculation
        where (maskPipeAc)
            kQ2 = (1.0 / length2) * kQ2
            !%  adding real time levels to source
            kQ2 = kQ2 - (af(2) * flowrate2n0 + af(3) * flowrate2n1) * invdt
            ! ! !%  gamma term
            gammaQ = - af(1) * invdt - grav * ((mn2**2) / (rh2**(4.0/3.0))) * abs(velocity2)
            ! !%  adding gamma to source
            kQ2 = kQ2 + gammaQ * flowrate2
            ! !%  multiplying by lambda (Note that C = 1 for Q)
            kQ2 = kQ2 / (1.0 - dtau * gammaQ)
            ! !%  updated Q
            flowrate2new = flowrate2old + wrk * dtau * kQ2
        endwhere
        ! output print for debug
        ! print*, trim(subroutine_name)
        ! print*, 'gammaQ      ',gammaQ
        ! print*, 'KQ2         ',KQ2(997:1001)
        ! print*, 'flowrate2   ',flowrate2
        ! print*, 'flowrate2new',flowrate2new(997:1001)
        ! print*, 'flowrate2old',flowrate2old
        ! if (wrk == oneR) then
        ! endif
        
        ! release temporary arrays
        gammaQ = nullvalueR
        nullify(gammaQ)
        next_e2r_temparray = next_e2r_temparray - 1

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine Kmomentum3Pipe
    !
    !==========================================================================
    !==========================================================================
    !

    !
    !==========================================================================
    !==========================================================================
    !
    ! subroutine ac_rk2_update_auxiliary_variables &
    !     (e2r_Velocity_new, eMr_Velocity_new, e2r_Volume_new, eMr_Volume_new,  &
    !     e2r_Flowrate_new, eMr_Flowrate_new, elem2R, elem2I, elem2YN, elemMR,  &
    !     elemMI, elemMYN, faceR,  faceI, faceYN, bcdataDn, bcdataUp, steptime, &
    !     rkiteration, ID, numberPairs, ManningsN, Length, zBottom, xDistance,  &
    !     Breadth, widthDepthData, cellType)

    !     character(64) :: subroutine_name = 'ac_rk2_update_auxiliary_variables'

    !     real,      target, intent(in out)  :: elemMR(:,:)
    !     real,              intent(in out)  :: elem2R(:,:), faceR(:,:)
    !     integer,           intent(in out)  :: elem2I(:,:), elemMI(:,:), faceI(:,:)
    !     logical,           intent(in out)  :: elem2YN(:,:),elemMYN(:,:),faceYN(:,:)
    !     type(bcType),      intent(in out)  :: bcdataDn(:), bcdataUp(:)
    !     integer,           intent(in)      :: e2r_Velocity_new, eMr_Velocity_new
    !     integer,           intent(in)      :: e2r_Volume_new,   eMr_Volume_new
    !     integer,           intent(in)      :: e2r_Flowrate_new, eMr_Flowrate_new
    !     integer,           intent(in)      :: rkiteration
    !     real,              intent(in)      :: steptime

    !     integer, intent(in out)    :: ID(:)
    !     integer, intent(in out)    :: numberPairs(:)
    !     real,    intent(in out)    :: ManningsN(:)
    !     real,    intent(in out)    :: Length(:)
    !     real,    intent(in out)    :: zBottom(:)
    !     real,    intent(in out)    :: xDistance(:)
    !     real,    intent(in out)    :: Breadth(:)
    !     real,    intent(in out)    :: widthDepthData(:,:,:)
    !     type(string), intent(in out)   :: cellType(:)
    !     !--------------------------------------------------------------------------
    !     if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

    !     !%  advance all geometry and dynamics
    !     call element_geometry_update &
    !         (elem2R, elem2I, elem2YN, e2r_Volume_new, elemMR, elemMI, elemMYN, &
    !         eMr_Volume_new, faceR, faceI, bcdataDn, bcdataUp, steptime, 1, ID, &
    !         numberPairs, ManningsN, Length, zBottom, xDistance, Breadth,       &
    !         widthDepthData, cellType)

    !     !%  at this point, the channels and the junction main sections have the correct
    !     !%  geometry, but the junction branches have provisional geometry that is
    !     !%  a functoin of the old face free surface elevation
    !     call element_dynamics_update &
    !         (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, bcdataDn, &
    !         bcdataUp, e2r_Velocity_new, eMr_Velocity_new, e2r_Volume_new,       &
    !         eMr_Volume_new, e2r_Flowrate_new, eMr_Flowrate_new, steptime)

    !     !%  Updating the face values by interpolation from neighbor elements
    !     !%  This uses the estimated values from the branches
    !     call face_update &
    !         (elem2R, elem2I, elemMR, faceR, faceI, faceYN, &
    !         bcdataDn, bcdataUp, e2r_Velocity_new, eMr_Velocity_new, &
    !         e2r_Volume_new, eMr_Volume_new, steptime, rkiteration)

    !     !% fix the junction branches by interp with face values
    !     call element_geometry_branch_fix (elemMR, elemMI, faceR, faceI )

    !     !%  Ad hoc adjustment for V-shaped flowrates across a channel element
    !     ! if (setting%Method%AdjustVshapedFlowrate%Apply) then
    !     !     call adjust_Vshaped_flowrate (elem2R, faceR, elem2I, elem2YN)
    !     ! endif

    !     if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    ! end subroutine ac_rk2_update_auxiliary_variables
    !==========================================================================
    ! END OF MODULE artificial_compressibility
    !==========================================================================
end module artificial_compressibility
