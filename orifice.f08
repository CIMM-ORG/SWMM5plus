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
    use xsect_tables

    implicit none

    public :: orifice_step

    private

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

        real(4),      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        integer,           intent(in out)  :: faceI(:,:)
        real(4),      target, intent(in out)  :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
        logical,           intent(in out)  :: faceYN(:,:)
        real(4),              intent(in)      :: thiscoef

        real(4),  pointer     ::  volume2old(:), volume2new(:), velocity2old(:), velocity2new(:)
        real(4),  pointer     ::  volumeMold(:), volumeMnew(:), velocityMold(:), velocityMnew(:)
        real(4),  pointer     ::  oBreadth(:), oDischargeCoeff(:), oInletoffset(:)
        real(4),  pointer     ::  oFullDepth(:), oZbottom(:)
        real(4),  pointer     ::  oFlow(:), oEta(:), oLength(:), oArea(:)
        real(4),  pointer     ::  oPerimeter(:), oHyddepth(:), oHydradius(:)
        real(4),  pointer     ::  oTopwidth(:), hEffective(:)
        real(4),  pointer     ::  hCrest(:), hCrown(:), hCrit(:)
        real(4),  pointer     ::  cOrif(:), cWeir(:), subFactor(:), subCorrection(:)
        real(4),  pointer     ::  fEdn(:), fEup(:)

        integer, pointer   ::  iup(:), idn(:), dir(:)

        logical, pointer   ::  maskarrayUpSubmerge(:), maskarrayDnSubmerge(:)

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

        !%  pointers for orifice geometry and settings
        !%  input
        oBreadth           => elem2R(:,e2r_BreadthScale)
        oDischargeCoeff    => elem2R(:,e2r_DischargeCoeff1)
        oInletoffset       => elem2R(:,e2r_InletOffset)
        oFullDepth         => elem2R(:,e2r_FullDepth)
        oZbottom           => elem2R(:,e2r_Zbottom)

        !%  output
        oFlow              => elem2R(:,e2r_Flowrate)
        oEta               => elem2R(:,e2r_eta)
        oLength            => elem2R(:,e2r_Length)
        oArea              => elem2R(:,e2r_Area)
        oPerimeter         => elem2R(:,e2r_Perimeter)
        oHyddepth          => elem2R(:,e2r_HydDepth)
        oHydradius         => elem2R(:,e2r_HydRadius)
        oTopwidth          => elem2R(:,e2r_Topwidth)
        hEffective         => elem2R(:,e2r_Depth)

        !%  pointers for upstream and downstream faces
        iup          => elem2I(:,e2i_Mface_u)
        idn          => elem2I(:,e2i_Mface_d)

        !%  temporary space for elem2
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

        subCorrection => elem2R(:,e2r_Temp(next_e2r_temparray))
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        dir           => elem2I(:,e2i_Temp(next_e2i_temparray))
        next_e2i_temparray = utility_advance_temp_array (next_e2i_temparray,e2i_n_temp)

        maskarrayDnSubmerge  => elem2YN(:,e2YN_Temp(next_e2YN_temparray) )
        next_e2YN_temparray  = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

        maskarrayUpSubmerge  => elem2YN(:,e2YN_Temp(next_e2YN_temparray) )
        next_e2YN_temparray  = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

        !%  zero temporary arrays
        hCrown         = nullvalueR
        hCrest         = nullvalueR
        hCrit          = nullvalueR
        cOrif          = nullvalueR
        cWeir          = nullvalueR
        subFactor      = nullvalueR

        dir            = nullvalueI
        subCorrection  = oneR

        maskarrayDnSubmerge  = nullvalueL
        maskarrayUpSubmerge  = nullvalueL

        !% sets necessary orifice setting and find eta on orifice element
        call orifice_initialize &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,    &
            oInletoffset, oZbottom, oFullDepth, oLength, oEta, hCrown,  &
            hCrest, fEdn, fEup, iup, idn, dir, thiscoef)

        !% calculates the  equivalent orificeand weir discharge coefficients
        call orifice_equivalent_discharge_coefficient &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, oBreadth, &
            oFullDepth, oDischargeCoeff, hCrit, cOrif, cWeir)

        !% calculates effective head in orifice elements
        call orifice_effective_head &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
            oEta, fEup, fEdn, iup, idn, dir, hCrest, hCrown, hcrit,  &
            hEffective, subFactor, maskarrayDnSubmerge,              &
            maskarrayUpSubmerge)

        !% updates geometry in orifice elements
        call orifice_geometry &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,         &
            oBreadth, oFullDepth, oArea, oPerimeter, oHyddepth, oHydradius,  &
            oTopwidth, hEffective)

        !% Villemonte correction for downstream submergence
        call villemonte_orifice_submergence_correction &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, hcrest, &
            subCorrection, fEdn, fEup, iup, idn, maskarrayDnSubmerge)

        !% Villemonte correction for upstream submergence
        call villemonte_orifice_submergence_correction &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, hcrest, &
            subCorrection, fEup, fEdn, idn, iup, maskarrayUpSubmerge)

        !% calculates flow in orifice elements
        call orifice_flow &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2new,   &
            velocity2new, volumeMnew, velocityMnew, oFlow, cOrif, cWeir, oBreadth, &
            oFullDepth, oArea, hEffective, dir, subFactor, subCorrection, thiscoef)

        ! release temporary arrays
        hCrown         = nullvalueR
        hCrest         = nullvalueR
        hCrit          = nullvalueR
        cOrif          = nullvalueR
        cWeir          = nullvalueR
        subFactor      = nullvalueR
        subCorrection  = nullvalueR

        dir            = nullvalueI

        maskarrayUpSubmerge = nullvalueL
        maskarrayDnSubmerge = nullvalueL

        nullify(hCrest, hCrown, hCrit, cOrif, cWeir, subFactor, subCorrection, &
            dir, maskarrayUpSubmerge, maskarrayDnSubmerge)

        next_e2r_temparray  = next_e2r_temparray  - 7
        next_e2i_temparray  = next_e2i_temparray  - 1
        next_e2YN_temparray = next_e2YN_temparray - 2

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

    end subroutine orifice_step
    !
    !==========================================================================
    ! PRIVATE BELOW
    !==========================================================================
    !
    subroutine orifice_initialize &
        (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,    &
        inletoffset, zbottom, fulldepth, length, eta, crown, crest, &
        faceEtaDn, faceEtaUp, upFace, dnFace, dir, thiscoef)
        !
        character(64) :: subroutine_name = 'orifice_initialize'

        real(4),      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real(4),      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
        real(4),              intent(in)      :: thiscoef

        real(4),  pointer   :: inletoffset(:), zbottom(:), fulldepth(:)
        real(4),  pointer   :: length(:), eta(:), crest(:), crown(:)
        real(4),  pointer   :: faceEtaDn(:), faceEtaUp(:)

        integer, pointer :: upFace(:), dnFace(:), dir(:)

        integer :: mm
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% find orifice crest, crown, length , eta flow direction
        where ( (elem2I(:,e2i_elem_type) == eOrifice) )

            crest   = inletoffset + zbottom
            crown   = crest + fulldepth
            ! find the effective orifice length
            length  = min(twoR * dt * sqrt(grav * fulldepth), 200.0)
            ! set the free surface elevation at orifice element
            eta     = max(faceEtaDn(upFace), faceEtaup(dnFace))
            dir     = int(sign(oneR, (faceEtaDn(upFace) - faceEtaup(dnFace))))
        endwhere

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

    end subroutine orifice_initialize
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine orifice_equivalent_discharge_coefficient &
        (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, breadth, &
        fulldepth, coeffDischarge, critDepth, coeffOrif, coeffWeir)
        !
        character(64) :: subroutine_name = 'orifice_equivalent_discharge_coefficient'

        real(4),      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real(4),      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

        real(4),    pointer ::  breadth(:), fulldepth(:), coeffDischarge(:)
        real(4),    pointer ::  critDepth(:), coeffOrif(:), coeffWeir(:)

        integer :: mm
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% find critical height above opening where orifice flow
        !% turns into weir flow. It equals (Co/Cw)*(Area/Length)
        !% where Co is the orifice coeff., Cw is the weir coeff/sqrt(2g),
        !% Area is the area of the opening, and Length = circumference
        !% of the opening. For a basic sharp crested weir, Cw = 0.414.

        where ( (elem2I(:,e2i_orif_elem_type) == eBottomOrifice) .and. &
            (elem2I(:,e2i_geometry) == eRectangular) )

            critDepth = coeffDischarge * (fulldepth * breadth) / &
                (0.414 * twoR  * (fullDepth + breadth) )

            coeffWeir = coeffDischarge * (fulldepth * breadth) * &
                sqrt(twoR * grav * critDepth)

        elsewhere ( (elem2I(:,e2i_orif_elem_type) == eBottomOrifice) .and. &
            (elem2I(:,e2i_geometry) == eCircular) )

            critDepth = coeffDischarge * fulldepth / (0.414 * fourR)

            coeffWeir = coeffDischarge * (pi / fourR * fulldepth ** twoR) * &
                sqrt(twoR * grav * critDepth)

        elsewhere ( (elem2I(:,e2i_orif_elem_type) == eSideOrifice) .and. &
            (elem2I(:,e2i_geometry) == eRectangular))

            critDepth = fulldepth

            coeffWeir = coeffDischarge * (fulldepth * breadth) * &
                sqrt(grav * critDepth)

        elsewhere ( (elem2I(:,e2i_orif_elem_type) == eSideOrifice) .and. &
            (elem2I(:,e2i_geometry) == eCircular))

            critDepth = fulldepth
            coeffWeir = coeffDischarge * ((pi / fourR) * fulldepth ** twoR) * &
                sqrt(grav * critDepth)

        endwhere
        !% find effective orifice discharge coefficient.
        !% Co = CdAo(g)^0.5
        !% Cd = discharge coefficient
        !% Ao = Area of orifice opening
        where     ( (elem2I(:,e2i_elem_type) == eOrifice) .and. &
            (elem2I(:,e2i_geometry) == eCircular) )

            coeffOrif = coeffDischarge * (pi / fourR * fulldepth ** twoR) * &
                sqrt(twoR * grav)

        elsewhere ( (elem2I(:,e2i_elem_type) == eOrifice) .and. &
            (elem2I(:,e2i_geometry) == eRectangular) )

            coeffOrif = coeffDischarge * (fulldepth * breadth) * &
                sqrt(twoR * grav)

        endwhere


        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

    end subroutine orifice_equivalent_discharge_coefficient
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine orifice_effective_head &
        (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, eta,        &
        faceEtaUp, faceEtaDn, upFace, dnFace, dir, crest, crown, critDepth,  &
        effectiveHead, submergenceFactor, maskarray_dn_submergence,          &
        maskarray_up_submergence)
        !
        character(64) :: subroutine_name = 'orifice_effective_head'

        real(4),      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real(4),      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

        real(4),  pointer   :: crest(:), crown(:), effectiveHead(:), critDepth(:)
        real(4),  pointer   :: faceEtaUp(:), faceEtaDn(:), eta(:), submergenceFactor(:)

        logical, pointer :: maskarray_dn_submergence(:), maskarray_up_submergence(:)

        integer, pointer :: upFace(:), dnFace(:), dir(:)


        integer :: mm
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% effective head calculation for bottom orifice
        where      ( (elem2I(:,e2i_orif_elem_type) == eBottomOrifice) .and. &
            (eta .LE. crest) )

            effectiveHead = zeroR

        elsewhere  ( (elem2I(:,e2i_orif_elem_type) == eBottomOrifice) .and. &
            (eta .GT. crest) )

            effectiveHead = min( (eta - crest), dir * &
                (faceEtaDn(upFace) - faceEtaUp(dnFace)) )
            ! find fraction of critical height for which weir flow occurs
            submergenceFactor  = min(effectiveHead / critDepth, oneR)

            ! downstream submergence
            maskarray_dn_submergence = ((dir .GT. zeroI)              .and. &
                (submergenceFactor .LT. oneR) .and. &
                (faceEtaUp(dnFace) .GT. crest)      )
            ! upstream submergance
            maskarray_up_submergence = ((dir .LT. zeroI)              .and. &
                (submergenceFactor .LT. oneR) .and. &
                (faceEtaDn(upFace) .GT. crest)      )
        endwhere

        !% find degree of submergence for side orifice
        where      ( (elem2I(:,e2i_orif_elem_type) == eSideOrifice) .and. &
            (eta .LT. crown) .and. (crown .GT. crest))

            submergenceFactor = (eta - crest) / (crown - crest)

            ! downstream submergence
            maskarray_dn_submergence = ((dir .GT. zeroI)              .and. &
                (submergenceFactor .LT. oneR) .and. &
                (faceEtaUp(dnFace) .GT. crest)      )
            ! upstream submergance
            maskarray_up_submergence = ((dir .LT. zeroI)              .and. &
                (submergenceFactor .LT. oneR) .and. &
                (faceEtaDn(upFace) .GT. crest)      )

        elsewhere ( elem2I(:,e2i_orif_elem_type) == eSideOrifice )
            submergenceFactor = OneR
        endwhere

        !% effective head calculation for side orifice
        where      ( (elem2I(:,e2i_orif_elem_type) == eSideOrifice) .and. &
            (submergenceFactor .LE. zeroR) )

            effectiveHead = zeroR

        elsewhere ( (elem2I(:,e2i_orif_elem_type) == eSideOrifice) .and. &
            (submergenceFactor .GT. zeroR) .and. (submergenceFactor .LT. oneR))

            effectiveHead = eta - crest

        elsewhere (elem2I(:,e2i_orif_elem_type) == eSideOrifice)

            effectiveHead = min((eta - (crest + crown) / twoR), &
                dir * (faceEtaDn(upFace) - faceEtaUp(dnFace)))
        endwhere

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

    end subroutine orifice_effective_head
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine orifice_geometry &
        (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,   &
        breadth, fullDepth, area, perimeter, hyddepth, hydradius,  &
        topwidth, depth)
        !
        character(64) :: subroutine_name = 'orifice_geometry'


        real(4),      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real(4),      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

        real(4),  pointer     ::  breadth(:), fullDepth(:), area(:)
        real(4),  pointer     ::  perimeter(:), hyddepth(:), hydradius(:)
        real(4),  pointer     ::  topwidth(:), depth(:)

        real(4)               ::  YoverYfull, Afull

        integer :: mm, ii
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        do ii=1, size(depth,1)
            if      ( (elem2I(ii,e2i_geometry)  == eCircular)          .and. &
                (elem2I(ii,e2i_elem_type) == eOrifice )          ) then

                YoverYfull  = depth(ii) / fulldepth(ii)
                Afull       = pi / fourR * fullDepth(ii) ** twoR

                if (YoverYfull .LE. zeroR) then
                    area(ii) = zeroR
                else
                    area(ii) = Afull * table_lookup(YoverYfull, ACirc, NACirc)
                endif
                topwidth(ii)   = fulldepth(ii) * table_lookup(YoverYfull, WCirc, NWCirc)
                hyddepth(ii)   = depth(ii)
                hydradius(ii)  = onefourthR * fulldepth (ii) * table_lookup(YoverYfull, RCirc, NRCirc)
                perimeter(ii)  = area(ii) / hydradius(ii)

            elseif ( (elem2I(ii,e2i_geometry)  == eRectangular) .and. &
                (elem2I(ii,e2i_elem_type) == eOrifice )       ) then

                area(ii)        = depth(ii) * breadth(ii)
                topwidth(ii)    = breadth(ii)
                hyddepth(ii)    = depth(ii)
                perimeter(ii)   = breadth(ii) + twoR * hyddepth(ii)
                hydradius(ii)   = area(ii) / perimeter(ii)
            endif
        enddo

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

    end subroutine orifice_geometry
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine villemonte_orifice_submergence_correction &
        (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, crest, &
        submerganceCorrection, faceEtaDn, faceEtaUp, upFace, dnFace,    &
        maskarray_submergence)
        !
        character(64) :: subroutine_name = 'villemonte_orifice_submergence_correction'

        real(4),      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real(4),      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

        real(4),  pointer ::  crest(:), submerganceCorrection(:)
        real(4),  pointer ::  faceEtaDn(:), faceEtaUp(:)

        integer, pointer :: upFace(:), dnFace(:)

        logical, pointer :: maskarray_submergence(:)

        integer :: mm
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% calculate the submergance factor for different orifice according to Villemonte 1974
        !% this only applies if the niminal d/s depth is higher than the crest
        where (maskarray_submergence)

            submerganceCorrection = (oneR - ((faceEtaUp(dnFace) - crest) / &
                (faceEtaDn(upFace) - crest)) ** 1.5) ** 0.385

        endwhere

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

    end subroutine villemonte_orifice_submergence_correction
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine orifice_flow &
        (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2new, &
        velocity2new, volumeMnew, velocityMnew, flow, coeffOrif, coeffWeir,  &
        breadth, fulldepth, area, effectiveHead, dir, submergenceFactor,     &
        submerganceCorrection, thiscoef)
        !
        character(64) :: subroutine_name = 'orifice_flow'


        real(4),      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real(4),      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
        real(4),              intent(in)      :: thiscoef


        real(4),  pointer   ::  volume2new(:), velocity2new(:), volumeMnew(:), velocityMnew(:)
        real(4),  pointer   ::  flow(:), coeffOrif(:), coeffWeir(:), breadth(:), fulldepth(:)
        real(4),  pointer   ::  area(:), submergenceFactor(:)
        real(4),  pointer   ::  effectiveHead(:),submerganceCorrection(:)

        integer, pointer :: dir(:)

        integer :: mm
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% flow and volume calculation
        where     ( (elem2I(:,e2i_elem_type) == eOrifice) .and. &
            (submergenceFactor .LE. zeroR) .or. (effectiveHead == zeroR))

            flow         = zeroR
            velocity2new = zeroR
            volume2new   = zeroR

        elsewhere ( (elem2I(:,e2i_elem_type) == eOrifice)      .and. &
            (submergenceFactor .LT. oneR) )

            flow         = dir * coeffWeir * submerganceCorrection * &
                submergenceFactor ** 1.5
            velocity2new = flow / area
            volume2new   = flow * dt

        elsewhere (elem2I(:,e2i_elem_type) == eOrifice )

            flow         = dir * coeffOrif * submerganceCorrection * &
                sqrt(abs(effectiveHead))
            velocity2new = flow / area
            volume2new   = flow * dt
        endwhere

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine orifice_flow
    !
    !==========================================================================
    ! END OF MODULE orifice
    !==========================================================================
    !
end module orifice
