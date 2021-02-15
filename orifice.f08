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
        (e2r_Volume_col, e2r_Velocity_col, elem2R, elemMR, faceI, faceR, &
        faceYN, elem2I, elemMI, elem2YN, elemMYN, thiscoef)
        !
        character(64) :: subroutine_name = 'orifice_step'

        ! indexes for old/new volume and velocity storage
        integer,   intent(in) :: e2r_Volume_col, e2r_Velocity_col

        real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        integer,           intent(in out)  :: faceI(:,:)
        real,      target, intent(in out)  :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
        logical,           intent(in out)  :: faceYN(:,:)
        real,              intent(in)      :: thiscoef

        real,  pointer     ::  volume2(:), velocity2(:), oBreadth(:) 
        real,  pointer     ::  oDischargeCoeff(:), oInletoffset(:)
        real,  pointer     ::  oFullDepth(:), oZbottom(:), oFlow(:), oDepth(:)
        real,  pointer     ::  oEta(:), oLength(:), oArea(:),oPerimeter(:)
        real,  pointer     ::  oHyddepth(:), oHydradius(:), oTopwidth(:)
        real,  pointer     ::  hEffective(:), oCrest(:), oCrown(:), oFullArea(:)
        real,  pointer     ::  hCrit(:), cOrif(:), cWeir(:)
        real,  pointer     ::  subFactor(:), subCorrection(:)
        real,  pointer     ::  fEdn(:), fEup(:)

        integer, pointer   ::  iup(:), idn(:), dir(:)

        logical, pointer   ::  maskarrayUpSubmerge(:), maskarrayDnSubmerge(:)

        integer :: mm
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !%  pointers for convenience in notation
        volume2   => elem2R(:,e2r_Volume_col)
        velocity2 => elem2R(:,e2r_Velocity_col)

        !%  pointers for orifice geometry and settings
        !%  input
        oBreadth           => elem2R(:,e2r_BreadthScale)
        oDischargeCoeff    => elem2R(:,e2r_DischargeCoeff1)
        oInletoffset       => elem2R(:,e2r_InletOffset)
        oFullDepth         => elem2R(:,e2r_FullDepth)
        oZbottom           => elem2R(:,e2r_Zbottom)
        oCrown             => elem2R(:,e2r_Zcrown)
        oFullArea          => elem2R(:,e2r_FullArea)
        fEdn               => faceR(:,fr_Eta_d)
        fEup               => faceR(:,fr_Eta_u)

        !%  output
        oFlow              => elem2R(:,e2r_Flowrate)
        oEta               => elem2R(:,e2r_eta)
        oLength            => elem2R(:,e2r_Length)
        oArea              => elem2R(:,e2r_Area)
        oPerimeter         => elem2R(:,e2r_Perimeter)
        oHyddepth          => elem2R(:,e2r_HydDepth)
        oHydradius         => elem2R(:,e2r_HydRadius)
        oTopwidth          => elem2R(:,e2r_Topwidth)
        oDepth             => elem2R(:,e2r_Depth)

        !%  pointers for upstream and downstream faces
        iup          => elem2I(:,e2i_Mface_u)
        idn          => elem2I(:,e2i_Mface_d)

        !%  temporary space for elem2
        hEffective    => elem2R(:,e2r_Temp(next_e2r_temparray))
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        oCrest        => elem2R(:,e2r_Temp(next_e2r_temparray))
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
        hEffective     = nullvalueR
        oCrest         = nullvalueR
        hCrit          = nullvalueR
        cOrif          = nullvalueR
        cWeir          = nullvalueR
        subFactor      = nullvalueR
        subCorrection  = oneR

        dir            = nullvalueI
    
        maskarrayDnSubmerge  = nullvalueL
        maskarrayUpSubmerge  = nullvalueL

        !% sets necessary orifice setting and find eta on orifice element
        call orifice_initialize &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,   &
            oInletoffset, oZbottom, oFullDepth, oLength, oEta, oCrown,  &
            oCrest, fEdn, fEup, iup, idn, dir)

        !% calculates the  equivalent orificeand weir discharge coefficients
        call orifice_equivalent_discharge_coefficient &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,   & 
            oBreadth, oFullDepth, oDischargeCoeff, hCrit, cOrif, cWeir)

        !% calculates effective head in orifice elements
        call orifice_effective_head &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
            oEta, fEup, fEdn, iup, idn, dir, oCrest, oCrown, hcrit,   &
            hEffective, subFactor, maskarrayDnSubmerge,               &
            maskarrayUpSubmerge)

        !% updates geometry in orifice elements
        call orifice_geometry &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,      &
            oBreadth, oFullDepth, oFullArea, oArea, oPerimeter, oHyddepth, &
            oHydradius, oTopwidth, oDepth, oEta, oCrest, oZbottom)

        !% Villemonte correction for downstream submergence
        call villemonte_orifice_submergence_correction &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, ocrest, &
            subCorrection, fEdn, fEup, iup, idn, maskarrayDnSubmerge)

        !% Villemonte correction for upstream submergence
        call villemonte_orifice_submergence_correction &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, ocrest, &
            subCorrection, fEup, fEdn, idn, iup, maskarrayUpSubmerge)

        !% calculates flow in orifice elements
        call orifice_flow &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2, &
            velocity2, oFlow, cOrif, cWeir, oBreadth, oFullDepth, oArea,       &
            hEffective, dir, subFactor, subCorrection, thiscoef)
        
        print*,'--------------------------------------------'
        print*,'Orifice values at ', subroutine_name
        print*
        print*, oEta(28), oEta(69), 'eta'
        print*
        print*, oZbottom(28), oZbottom(69), 'zbottom'
        print*
        print*, hEffective(28), hEffective(69), 'effective head'
        print*
        print*, oDepth(28), oDepth(69), 'depth'
        print*
        print*, oFullDepth(28), oFullDepth(69), 'fulldepth'
        print*
        print*, oFlow(28), oFlow(69), 'flow'
        print*
        print*, velocity2(28), velocity2(69), 'velocity'
        print*
        print*, oArea(28), oArea(69), 'area'
        print*
        print*, volume2(28), volume2(69), 'volume'
        print*
        ! print*, 'orifice debug: press return to continue'
        ! read(*,*)

        ! release temporary arrays
        hEffective     = nullvalueR
        oCrest         = nullvalueR
        hCrit          = nullvalueR
        cOrif          = nullvalueR
        cWeir          = nullvalueR
        subFactor      = nullvalueR
        subCorrection  = nullvalueR

        dir            = nullvalueI

        maskarrayUpSubmerge = nullvalueL
        maskarrayDnSubmerge = nullvalueL

        nullify(hEffective, oCrest, hCrit, cOrif, cWeir, subFactor, subCorrection, &
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
        (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,   &
        inletoffset, zbottom, fulldepth, length, eta, crown, crest, &
        faceEtaDn, faceEtaUp, upFace, dnFace, dir)
        !
        character(64) :: subroutine_name = 'orifice_initialize'

        real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real,      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

        real,    intent(inout)  :: crest(:), length(:), eta(:)
        real,    intent(in)     :: inletoffset(:), zbottom(:), fulldepth(:)
        real,    intent(in)     :: crown(:), faceEtaDn(:), faceEtaUp(:)
        integer, intent(inout)  :: dir(:)
        integer, intent(in)     :: upFace(:), dnFace(:)

        integer :: mm
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% find orifice crest, crown, length , eta flow direction
        where ( (elem2I(:,e2i_elem_type) == eOrifice) )

            crest   = inletoffset + zbottom
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
        !%  calculates equivalent discharge coefficients depending on
        !%  whether or not the flow is weir or orifice.
        !
        character(64) :: subroutine_name = 'orifice_equivalent_discharge_coefficient'

        real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real,      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

        real,    intent(inout) ::  critDepth(:), coeffOrif(:), coeffWeir(:)  
        real,    intent(in)    ::  breadth(:), fulldepth(:), coeffDischarge(:)

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% find critical height above opening where orifice flow
        !% turns into weir flow. It equals (Co/Cw)*(Area/Length)
        !% where Co is the orifice coeff., Cw is the weir coeff/sqrt(2g),
        !% Area is the area of the opening, and Length = circumference
        !% of the opening. For a basic sharp crested weir, Cw = 0.414.

        where ( (elem2I(:,e2i_orif_elem_type) == eBottomOrifice) .and. &
                (elem2I(:,e2i_geometry) == eRectangular        )       )

            critDepth = coeffDischarge * (fulldepth * breadth) / &
                (0.414 * twoR  * (fullDepth + breadth) )

            coeffWeir = coeffDischarge * (fulldepth * breadth) * &
                sqrt(twoR * grav * critDepth)

        elsewhere ( (elem2I(:,e2i_orif_elem_type) == eBottomOrifice) .and. &
                    (elem2I(:,e2i_geometry) == eCircular           )       )

            critDepth = coeffDischarge * fulldepth / (0.414 * fourR)

            coeffWeir = coeffDischarge * (pi / fourR * fulldepth ** twoR) * &
                sqrt(twoR * grav * critDepth)

        elsewhere ( (elem2I(:,e2i_orif_elem_type) == eSideOrifice) .and. &
                    (elem2I(:,e2i_geometry) == eRectangular      )       )

            critDepth = fulldepth

            coeffWeir = coeffDischarge * (fulldepth * breadth) * &
                sqrt(grav * critDepth)

        elsewhere ( (elem2I(:,e2i_orif_elem_type) == eSideOrifice) .and. &
                    (elem2I(:,e2i_geometry) == eCircular         )       )

            critDepth = fulldepth
            coeffWeir = coeffDischarge * ((pi / fourR) * fulldepth ** twoR) * &
                sqrt(grav * critDepth)

        endwhere

        !% find effective orifice discharge coefficient.
        !% Co = CdAo(g)^0.5
        !% Cd = discharge coefficient
        !% Ao = Area of orifice opening

        where ( (elem2I(:,e2i_elem_type) == eOrifice) .and. &
                (elem2I(:,e2i_geometry) == eCircular)       )

            coeffOrif = coeffDischarge * (pi / fourR * fulldepth ** twoR) * &
                sqrt(twoR * grav)

        elsewhere ( (elem2I(:,e2i_elem_type) == eOrifice   ) .and. &
                    (elem2I(:,e2i_geometry) == eRectangular)       )

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
        (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, eta,       &
        faceEtaUp, faceEtaDn, upFace, dnFace, dir, crest, crown, critDepth,  &
        effectiveHead, submergenceFactor, maskarray_dn_submergence,          &
        maskarray_up_submergence)
        !
        character(64) :: subroutine_name = 'orifice_effective_head'

        real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real,      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

        real,    intent(inout)   :: effectiveHead(:), submergenceFactor(:)
        real,    intent(in)      :: crest(:), crown(:), eta(:), critDepth(:)
        real,    intent(in)      :: faceEtaUp(:), faceEtaDn(:)

        logical, intent(inout)   :: maskarray_dn_submergence(:), maskarray_up_submergence(:)

        integer, intent(in)      :: upFace(:), dnFace(:), dir(:)

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% effective head calculation for bottom orifice
        where ( (elem2I(:,e2i_orif_elem_type) == eBottomOrifice) .and. &
                (eta .LE. crest                                )       )

            effectiveHead = zeroR

        elsewhere  ( (elem2I(:,e2i_orif_elem_type) == eBottomOrifice) .and. &
                     (eta .GT. crest                                )       )

            effectiveHead = min( (eta - crest), &
                    dir * (faceEtaDn(upFace) - faceEtaUp(dnFace)) )

            !% find fraction of critical height for which weir flow occurs
            submergenceFactor  = min(effectiveHead / critDepth, oneR)

            !% downstream submergence
            maskarray_dn_submergence = ((dir .GT. zeroI)               .and. &
                                        (submergenceFactor .LT. oneR ) .and. &
                                        (faceEtaUp(dnFace) .GT. crest)       )
            !% upstream submergance
            maskarray_up_submergence = ((dir .LT. zeroI)               .and. &
                                        (submergenceFactor .LT. oneR ) .and. &
                                        (faceEtaDn(upFace) .GT. crest)       )
        endwhere

        !% find degree of submergence for side orifice
        where ( (elem2I(:,e2i_orif_elem_type) == eSideOrifice) .and. &
                (eta .LT. Crown                              )       )

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
        where ( (elem2I(:,e2i_orif_elem_type ) == eSideOrifice) .and. &
                (submergenceFactor .LE. zeroR)                        )

            effectiveHead = zeroR

        elsewhere ( (elem2I(:,e2i_orif_elem_type ) == eSideOrifice) .and. &
                    (submergenceFactor .GT. zeroR)                  .and. &
                    (submergenceFactor .LT. oneR )                        )

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
        (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, breadth,   &
        fullDepth, fullArea, area, perimeter, hyddepth, hydradius, topwidth, &
        depth, eta, crest, zbottom)
        !
        !%  geometry handler for orifice elements
        !
        character(64) :: subroutine_name = 'orifice_geometry'


        real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real,      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)


        real,    intent(inout) ::  area(:), perimeter(:), hyddepth(:)
        real,    intent(inout) ::  hydradius(:), topwidth(:), depth(:)
        real,    intent(in)    ::  breadth(:), fullDepth(:), fullArea(:)
        real,    intent(in)    ::  eta(:), crest(:), zbottom(:)

        real,    pointer       ::  YoverYfull(:)
        logical, pointer       ::  maskCircularOrifice(:)

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% temporary pointer allocation fo circular for geometry update
        YoverYfull => elem2R(:,e2r_Temp(next_e2r_temparray))
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        !% temporary mask to find circular orifice elements
        maskCircularOrifice => elem2YN(:,e2YN_Temp(next_e2YN_temparray))
        next_e2YN_temparray = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

        YoverYfull          = nullvalueR
        maskCircularOrifice = nullvalueL

        !% ==========================================================================
        !%  rectangular orifice geometry handler
        !% ==========================================================================
        where ( (elem2I(:,e2i_elem_type) == eOrifice    ) .and. &
                (elem2I(:,e2i_geometry)  == eRectangular)       ) 

                depth       = min( max((eta - crest), zeroR), fullDepth)
                area        = depth * breadth
                topwidth    = breadth
                hyddepth    = depth
                perimeter   = breadth + twoR * hyddepth
                hydradius   = area / perimeter
        endwhere
        !% ==========================================================================
        !%  circular orifice geometry handler 
        !% ==========================================================================
        !%  mask circular orifice elements
        maskCircularOrifice = ( (elem2I(:,e2i_elem_type) == eOrifice ) .and. &
                                (elem2I(:,e2i_geometry)  == eCircular)       ) 

        !%  find Y/Yfull for table interpolation
        where (maskCircularOrifice)
            depth       = min( max((eta - crest), zeroR), fullDepth)
            YoverYfull  = depth / fullDepth
        endwhere

        !% find normalized area using Y/Yfull from lookup tables
        !% normalized area (A/Afull) is saved in the area column
        call table_lookup_mask &
            (elem2I, elem2R, area, YoverYfull, ACirc, NACirc, maskCircularOrifice, &
            e2i_Temp, next_e2i_temparray, e2i_n_temp)

        !% find normalized topwidth using Y/Yfull from lookup tables
        !% normalized topwidth (W/Wmax) is saved in the topwidth column
        call table_lookup_mask &
            (elem2I, elem2R, topwidth, YoverYfull, WCirc, NWCirc, maskCircularOrifice, &
            e2i_Temp, next_e2i_temparray, e2i_n_temp)

        !% find normalized hydraulic radius using Y/Yfull from lookup tables
        !% normalized hydraulic radius (R/Rmax) is saved in the hydradius column
        call table_lookup_mask &
            (elem2I, elem2R, hydradius, YoverYfull, RCirc, NRCirc, maskCircularOrifice, &
            e2i_Temp, next_e2i_temparray, e2i_n_temp)

        where (maskCircularOrifice)
            area      = fullArea  * area 
            topwidth  = fulldepth * topwidth
            hyddepth  = area / topwidth
            hydradius = onefourthR * fulldepth * hydradius
            perimeter = area / hydradius
        endwhere
        !% ==========================================================================
        !%  nullify and release temporary pointers
        YoverYfull          = nullvalueR
        maskCircularOrifice = nullvalueL
        nullify(YoverYfull, maskCircularOrifice)

        next_e2r_temparray  = next_e2r_temparray  - 1
        next_e2YN_temparray = next_e2YN_temparray - 1

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

        real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real,      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

        real,    intent(inout) ::  submerganceCorrection(:)
        real,    intent(in)    ::  crest(:), faceEtaDn(:), faceEtaUp(:)
        integer, intent(in)    ::  upFace(:), dnFace(:)
        logical, intent(in)    ::  maskarray_submergence(:)

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
        (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2, &
        velocity2, flow, coeffOrif, coeffWeir, breadth, fulldepth, area,   &
        effectiveHead, dir, submergenceFactor, submerganceCorrection,      &
        thiscoef)
        !
        !%  calculate orifice flow using effective head
        !
        character(64) :: subroutine_name = 'orifice_flow'


        real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real,      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
        real,              intent(in)      :: thiscoef


        real,    intent(inout)   ::  volume2(:), velocity2(:), flow(:)
        real,    intent(in)      ::  coeffOrif(:), coeffWeir(:), breadth(:)
        real,    intent(in)      ::  fulldepth(:), area(:), submergenceFactor(:)
        real,    intent(in)      ::  effectiveHead(:),submerganceCorrection(:)
        integer, intent(in)      ::  dir(:)

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% flow and volume calculation
        where ( (elem2I(:,e2i_elem_type) == eOrifice) .and. &
                (submergenceFactor .LE. zeroR       ) .or.  &
                (effectiveHead == zeroR             )       )

            flow      = zeroR
            velocity2 = zeroR
            volume2   = zeroR

        elsewhere ( (elem2I(:,e2i_elem_type) == eOrifice) .and. &
                    (submergenceFactor .LT. oneR        )       )

            flow      = dir * coeffWeir * submerganceCorrection * &
                submergenceFactor ** 1.5
            velocity2 = flow / area
            volume2   = flow * dt * thiscoef

        elsewhere (elem2I(:,e2i_elem_type) == eOrifice )

            flow      = dir * coeffOrif * submerganceCorrection * &
                sqrt(abs(effectiveHead))
            velocity2 = flow / area
            volume2   = flow * dt * thiscoef
        endwhere

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine orifice_flow
    !
    !==========================================================================
    ! END OF MODULE orifice
    !==========================================================================
    !
end module orifice
