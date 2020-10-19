! module ac_convergence_loop
!
! pseudo time marching loop of artificial compressibility convergence
!
!==========================================================================
!
module ac_convergence_loop

    use array_index
    use artificial_compressibility
    use bc
    use data_keys
    use debug
    use diagnostic
    use explicit_euler
    use face_values
    use globals
    use output
    use runge_kutta
    use setting_definition
    use utility

    implicit none

    private

    public  :: pseudo_time_marching

    integer :: debuglevel = 0

    integer,   parameter :: idummy = 0

contains
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine pseudo_time_marching &
        (elem2R, elemMR, faceR, elem2I, elemMI, faceI, elem2YN, elemMYN,    &
        faceYN, bcdataDn, bcdataUp, linkI, thisStep, thisTime, AnormH,      &
        AnormQ, AnormHlast, AnormQlast, TnormH, TnormQ, RTnormH, RTnormQ,   &
        RLnormH,  RLnormQ, debugfile, diagnostic, threadedfile, ID,         &
        numberPairs, ManningsN, Length, zBottom, xDistance, Breadth,        &
        widthDepthData, cellType)
        !
        ! pseudo time marching loop of artificial compressibility convergence
        !
        character(64) :: subroutine_name = 'pseudo_time_marching'

        real,      target, intent(in out) :: elem2R(:,:),  elemMR(:,:),  faceR(:,:)
        integer,   target, intent(in out) :: elem2I(:,:),  elemMI(:,:),  faceI(:,:)
        logical,   target, intent(in out) :: elem2YN(:,:), elemMYN(:,:), faceYN(:,:)

        type(bcType),                  intent(in out) :: bcdataDn(:), bcdataUp(:)
        type(debugfileType),  target,  intent(in)     :: debugfile(:)
        type(diagnosticType), target,  intent(in out) :: diagnostic(:)
        type(threadedfileType),        intent(in)     :: threadedfile(:)

        integer,                       intent(in)     :: linkI(:,:)
        real,     intent(in)         :: thisTime
        integer,  intent(in)         :: thisStep   
        real,     intent(in out)     :: AnormH(:), AnormQ(:), AnormHlast(:), AnormQlast(:)
        real,     intent(in out)     :: TnormH(:), TnormQ(:)
        real,     intent(in out)     :: RTnormH(:), RTnormQ(:), RLnormH(:), RLnormQ(:)

        real    :: dtau, tauTime, NX
        integer :: thisIter, iterMax, iterMin
        logical :: isConverged

        logical, pointer :: maskAC(:)

        integer, intent(in out)    :: ID(:)
        integer, intent(in out)    :: numberPairs(:)
        real,    intent(in out)    :: ManningsN(:)
        real,    intent(in out)    :: Length(:)
        real,    intent(in out)    :: zBottom(:)
        real,    intent(in out)    :: xDistance(:)
        real,    intent(in out)    :: Breadth(:)
        real,    intent(in out)    :: widthDepthData(:,:,:)
        type(string), intent(in out)   :: cellType(:)

        integer :: ii
        logical :: acCycle(2)

        character(len=32) :: outdataName
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !%  temporary pointer to find the full pipes that become open
        maskAC   => elem2YN(:,e2YN_Temp(next_e2YN_temparray))
        next_e2YN_temparray = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

        maskAC = nullvalueL
        maskAC = (elem2I(:,e2i_solver) == AC)

        !%  set AC time step based on the dtau factor.
        setting%DefaultAC%dtau = dt * setting%DefaultAC%dtauFactor%dtdtau
        !%  initialize iteration counter for the AC loops
        thisIter = zeroI
        !%  initialization for AC loop
        isConverged = .false.
        tauTime = zeroR
        dtau = setting%DefaultAC%dtau
        !%  count the number of element in AC loop
        !%  HACK: Check with Dr. Hodges
        NX = count(elem2I(:,e2i_solver) == AC)

        if (thisStep == 1) then
            !%  Allows first step to take a large number of AC pseudo-time
            !%  iterations because the initial conditions typically are 
            !%  inconsistent
            iterMax = setting%DefaultAC%Iter%Firststep
        else
            iterMax = setting%DefaultAC%Iter%Max
        endif
        iterMin = setting%DefaultAC%Iter%Min

        !%  pseudo time loop
        do while ( (isConverged .eqv. .false.) .and. (thisIter < iterMax) )
            thisIter = thisIter + oneI
            !%  store norms from the last iteration
            AnormHlast = AnormH
            AnormQlast = AnormQ

            !%  acCycle mask ensures only the second RK step is take in AC convergence
            acCycle(1) = .false.
            acCycle(2) = .true.

            call rk2 & 
                (elem2R, elemMR, elem2I, elemMI, faceR, faceI, elem2YN, elemMYN, faceYN, &
                bcdataDn, bcdataUp, thistime, dt, ID, numberPairs, ManningsN, Length,   &
                zBottom, xDistance, Breadth, widthDepthData, cellType, acCycle)

            !%  baseline for tau convergence -- change from the first RK full step
            where (elem2I(:,e2i_solver) == AC)
                elem2R(:,e2r_CtestH0) = elem2R(:,e2r_eta) * elem2R(:,e2r_area) - &
                                        elem2R(:,e2r_Eta_N0) * elem2R(:,e2r_Area_N0)
                elem2R(:,e2r_CtestQ0) = elem2R(:,e2r_flowrate) - elem2R(:,e2r_Flowrate_N0)
            endwhere

            !%  norms for testing convergence
            AnormH(1) = sum(abs(elem2R(:,e2r_CtestH1)), maskAC) / Nx
            AnormH(2) = sqrt(sum(elem2R(:,e2r_CtestH1) ** twoR, maskAC)) / Nx
            AnormH(3) = maxval(abs(elem2R(:,e2r_CtestH1)), maskAC)

            AnormQ(1) = sum(abs(elem2R(:,e2r_CtestQ1)), maskAC) / Nx
            AnormQ(2) = sqrt(sum(elem2R(:,e2r_CtestQ1) ** twoR, maskAC)) / Nx
            AnormQ(3) = maxval(abs(elem2R(:,e2r_CtestQ1)), maskAC)

            TnormH(1) = sum(abs(elem2R(:,e2r_CtestH0)), maskAC) / Nx
            TnormH(2) = sqrt(sum(elem2R(:,e2r_CtestH0) ** twoR, maskAC)) / Nx
            TnormH(3) = maxval(abs(elem2R(:,e2r_CtestH0)), maskAC)

            TnormQ(1) = sum(abs(elem2R(:,e2r_CtestQ0)), maskAC) / Nx
            TnormQ(2) = sqrt(sum(elem2R(:,e2r_CtestQ0) ** twoR, maskAC)) / Nx
            TnormQ(3) = maxval(abs(elem2R(:,e2r_CtestQ0)), maskAC)
            

            RTnormH = AnormH / TnormH
            RTnormQ = AnormQ / TnormQ
            RLnormH = zeroR
            RLnormQ = zeroR

            where (AnormHlast > zeroR)
                RLnormH = abs(1.0 - AnormH / AnormHlast)
            endwhere

            where (AnormQlast > zeroR)
                RLnormQ = abs(1.0 - AnormQ / AnormQlast)
            endwhere
            !%  HACK: need same baseline for tau convergence for junctionpipe element


            !%  Norm testing for convergence
            !%  Either relative norm (slowing convergence) or the
            !%  absolute norm (dimensionally small change) may indicate
            !%  convergence
            if (thisIter .GE. iterMin) then
                if ( (      (RLnormH(2) .LE. setting%DefaultAC%Convergence%Hrelative) &
                      .and. (RLnormQ(2) .LE. setting%DefaultAC%Convergence%Qrelative) &
                      .and. (RLnormH(1) .LE. 5.0  * setting%DefaultAC%Convergence%Hrelative) &
                      .and. (RLnormQ(1) .LE. 5.0  * setting%DefaultAC%Convergence%Qrelative) &
                      .and. (RLnormH(3) .LE. 10.0 * setting%DefaultAC%Convergence%Hrelative) &
                      .and. (RLnormQ(3) .LE. 10.0 * setting%DefaultAC%Convergence%Qrelative) ) &
                .or. &
                     (      (AnormH(2) .LE. setting%DefaultAC%Convergence%Habsolute) &
                      .and. (AnormQ(2) .LE. setting%DefaultAC%Convergence%Qabsolute) &
                      .and. (AnormH(1) .LE. 5.0  * setting%DefaultAC%Convergence%Habsolute) &
                      .and. (AnormQ(1) .LE. 5.0  * setting%DefaultAC%Convergence%Qabsolute) &
                      .and. (AnormH(3) .LE. 10.0 * setting%DefaultAC%Convergence%Habsolute) &
                      .and. (AnormQ(3) .LE. 10.0 * setting%DefaultAC%Convergence%Qabsolute) ) ) then

                    isConverged = .true.
                    tauTime = zeroR
                else
                    tauTime = tauTime + dtau
                endif
            else
                tauTime = tauTime +dtau
            endif

            if (setting%DefaultAC%PrintConvergence) then
                if (isConverged) then
                    print*,'------------------------------------'
                    print*, 'AC converged at iteration ',thisiter
                    print*, 'H absolute norms:  ',AnormH(1), AnormH(2), AnormH(3)
                    print*, 'Q absolute norms:  ',AnormQ(1), AnormQ(2), AnormQ(3)
                    print*, 'H relative Lnorms: ',RLnormH(1), RLnormH(2), RLnormH(3)
                    print*, 'Q relative Lnorms: ',RLnormQ(1), RLnormQ(2), RLnormQ(3)
                    print*,'------------------------------------'

                else
                    print*,'------------------------------------'
                    print*, '**** AC residual larger than target'
                    print*, 'H absolute norms: ',AnormH(1), AnormH(2), AnormH(3)
                    print*, 'Q absolute norms: ',AnormQ(1), AnormQ(2), AnormQ(3)
                    print*, 'H relative Lnorms: ',RLnormH(1), RLnormH(2), RLnormH(3)
                    print*, 'Q relative Lnorms: ',RLnormQ(1), RLnormQ(2), RLnormQ(3)

                endif
            endif

            !%  HACK: add density anomay corrections

        end do

        ! release temporary arrays
        maskAC  = nullvalueL
        nullify(maskAC)
        next_e2YN_temparray = next_e2YN_temparray - 1

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
       
    end subroutine pseudo_time_marching
    !
    !==========================================================================
    ! END OF MODULE ac_convergence_loop
    !==========================================================================
end module ac_convergence_loop
