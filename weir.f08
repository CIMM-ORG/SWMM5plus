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
        (e2r_Volume_col, e2r_Velocity_col, elem2R, elemMR, faceI, faceR, &
        faceYN, elem2I, elemMI, elem2YN, elemMYN, thiscoef)
        !
        character(64) :: subroutine_name = 'weir_step'

        ! indexes for old/new volume and velocity storage
        integer,   intent(in) :: e2r_Volume_col, e2r_Velocity_col

        real(8),      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        integer,           intent(in out)  :: faceI(:,:)
        real(8),      target, intent(in out)  :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
        logical,           intent(in out)  :: faceYN(:,:)
        real(8),              intent(in)      :: thiscoef

        real(8),  pointer     ::  volume2(:), velocity2(:), wFlow(:), wEta(:)
        real(8),  pointer     ::  wPerimeter(:), wHyddepth(:), wHydradius(:), wTopwidth(:)
        real(8),  pointer     ::  wArea(:), lEffective(:), hEffective(:), wCrest(:)
        real(8),  pointer     ::  cOrif(:), subFactor1(:), subFactor2(:)
        real(8),  pointer     ::  wBreadth(:), wInletoffset(:), wFullDepth(:), wZbottom(:)
        real(8),  pointer     ::  wSideSlope(:), wEndContractions(:), cTriangular(:), wCrown(:)
        real(8),  pointer     ::  cRectangular(:), wLength(:), fEdn(:), fEup(:)
        integer, pointer   ::  iup(:), idn(:), dir(:)

        logical, pointer   ::  maskarrayUpSubmerge(:), maskarrayDnSubmerge(:)
        logical, pointer   ::  IsSurcharged(:)


        integer :: mm
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !%  pointers for convenience in notation
        volume2   => elem2R(:,e2r_Volume_col)
        velocity2 => elem2R(:,e2r_Velocity_col)

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
        wCrown             => elem2R(:,e2r_Zcrown)
        fEdn               => faceR(:,fr_Eta_d)
        fEup               => faceR(:,fr_Eta_u)

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
        IsSurcharged       => elem2YN(:,e2YN_IsSurcharged)

        !%  pointers for upstream and downstream faces
        iup          => elem2I(:,e2i_Mface_u)
        idn          => elem2I(:,e2i_Mface_d)

        !%  temporary pointers
        !%  weir crest elevation
        wCrest       => elem2R(:,e2r_Temp(next_e2r_temparray))
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        !%  effective crest length for rectangular or trapezoidal weir
        lEffective   => elem2R(:,e2r_Temp(next_e2r_temparray))
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        !%  discharge coefficient for surcharged condition (orifice flow)
        cOrif        => elem2R(:,e2r_Temp(next_e2r_temparray))
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        !%  submergence factor for submergence
        subFactor1   => elem2R(:,e2r_Temp(next_e2r_temparray))
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        !  submergence factor for submergence
        subFactor2   => elem2R(:,e2r_Temp(next_e2r_temparray))
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)
        
        !%  flow direction
        dir          => elem2I(:,e2i_Temp(next_e2i_temparray))
        next_e2i_temparray = utility_advance_temp_array (next_e2i_temparray,e2i_n_temp)

        !%  mask for elem2YN d/s submergence
        maskarrayDnSubmerge  => elem2YN(:,e2YN_Temp(next_e2YN_temparray) )
        next_e2YN_temparray  = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

        !%  mask for elem2YN u/s submergence
        maskarrayUpSubmerge  => elem2YN(:,e2YN_Temp(next_e2YN_temparray) )
        next_e2YN_temparray  = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

        !%  initializing temporary arrays
        subFactor1        = oneR
        subFactor2        = oneR
        dir               = oneI

        maskarrayDnSubmerge  = nullvalueL
        maskarrayUpSubmerge  = nullvalueL

        !% set necessary weir setting and find eta on weir element
        call weir_initialize &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,     &
            wInletoffset, wZbottom, wCrown, wCrest, wFullDepth, wLength, &
            wEta, fEdn, fEup, iup, idn, dir, thiscoef)

        !% calculate the  equivalent orifice discharge coefficient while surcharged
        call weir_surcharge_coefficient &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2, &
            velocity2, wFlow, wSideSlope, cTriangular, cRectangular, wBreadth, &
            wArea, dir, wFullDepth, wEndContractions, lEffective, subFactor1,  &
            subFactor2, cOrif, thiscoef)

        !% calculate effective head on weir element
        call weir_effective_head &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
            wCrest, wCrown, wEta, fEup, fEdn, iup, idn, dir,         &
            hEffective, maskarrayDnSubmerge, maskarrayUpSubmerge,    &
            IsSurcharged)

        !% calculates weir geometrices for effective head
        call weir_geometry &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, wBreadth, &
            wSideSlope, wArea, wPerimeter, wHyddepth, wHydradius, wTopwidth,   &
            hEffective)

        !% calculate weir length and effective crest length
        call weir_effective_crest_length &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN,     &
            wEndContractions, hEffective, wBreadth, wSideSlope, lEffective)

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
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2,   &
            velocity2, wFlow, wArea, wSideSlope, cTriangular, cRectangular, dir, &
            hEffective, lEffective, subFactor1, subFactor2, thiscoef)

        !%  flow calculataion when flow is surcharged
        call weir_surcharge_flow &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2, &
            velocity2, wFlow, wCrest, wCrown, wEta, wArea, cOrif, hEffective,  &
            fEup, fEdn, iup, idn, dir, thiscoef, IsSurcharged)

        ! release temporary arrays
        lEffective     = nullvalueR
        subFactor1     = nullvalueR
        subFactor2     = nullvalueR
        cOrif          = nullvalueR
        wCrest         = nullvalueR
        dir            = nullvalueI

        nullify(lEffective, wCrest, wCrown, cOrif, subFactor1, subFactor2, &
                dir, maskarrayDnSubmerge, maskarrayUpSubmerge)

        next_e2r_temparray  = next_e2r_temparray  - 5
        next_e2i_temparray  = next_e2i_temparray  - 1
        next_e2YN_temparray = next_e2YN_temparray - 2

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


        real(8),      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real(8),      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
        real(8),              intent(in)      :: thiscoef

        real(8),    intent(inout) :: crest(:), length(:), eta(:)
        real(8),    intent(in)    :: inletoffset(:), zbottom(:), crown(:)
        real(8),    intent(in)    :: fullDepth(:), faceEtaDn(:), faceEtaUp(:)

        integer, intent(inout) :: dir(:)
        integer, intent(in)    :: upFace(:), dnFace(:)

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% find weir crest, length , eta flow direction
        where ( (elem2I(:,e2i_elem_type) == eWeir) )
            crest   =  inletoffset + zbottom
            ! find the effective weir length (changes depending on control setting)
            length  = min(twoR*dt*sqrt(grav*fullDepth), 200.0)
            ! set the free surface elevation at weir element
            eta  = max(faceEtaDn(upFace), faceEtaUp(dnFace))
            dir  = int(sign(oneR, ( faceEtaDn(upFace) - faceEtaUp(dnFace))))
        endwhere

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine weir_initialize
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine weir_surcharge_coefficient &
        (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2, &
        velocity2, flow, sideslope, cTrig, cRect, breadth, area, dir,      &
        fulldepth, endcontractions, crestlength, submergenceFactor1,       &
        submergenceFactor2, corif, thiscoef)
        !
        !% when weir is surcharged, the flow becomes orifice flow. this subroutine
        !% calculates the equivalent orifice discharge coefficient, cOrif
        !
        character(64) :: subroutine_name = 'weir_surcharge_coefficient'

        real(8),      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real(8),      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
        real(8),              intent(in)      :: thiscoef

        real(8),    intent(inout) ::  volume2(:), velocity2(:), flow(:), crestlength(:), corif(:)
        real(8),    intent(in)    ::  sideslope(:), cTrig(:), cRect(:), breadth(:), area(:)
        real(8),    intent(in)    ::  fulldepth(:), endcontractions(:)
        real(8),    intent(in)    ::  submergenceFactor1(:), submergenceFactor2(:)
        integer, intent(in)    ::  dir(:)

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% get effective crest length for maximum weir opening
        call weir_effective_crest_length &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
            endcontractions, fulldepth, breadth, sideslope, crestlength)

        !% get flow for maximum weir opening
        call weir_flow &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2, &
            velocity2, flow, area, sideslope, cTrig, cRect, dir, fulldepth,    &
            crestlength , submergenceFactor1, submergenceFactor2, thiscoef)

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
        (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, crest,     &
        crown, eta, faceEtaUp, faceEtaDn, upFace, dnFace, dir,effectivehead, &
        maskarray_dn_submergence, maskarray_up_submergence, is_surcharged)
        !
        !%  find the effective head on a weir
        !
        character(64) :: subroutine_name = 'weir_effective_head'


        real(8),      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real(8),      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

        real(8),    intent(inout)   :: effectivehead(:)
        real(8),    intent(in)      :: crest(:), crown(:), eta(:)
        real(8),    intent(in)      :: faceEtaUp(:), faceEtaDn(:)
        integer, intent(in)      :: upFace(:), dnFace(:), dir(:)
        logical, intent(inout)   :: maskarray_dn_submergence(:), maskarray_up_submergence(:)
        logical, intent(inout)   :: is_surcharged(:)

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% effective head calculation
        where ( (elem2I(:,e2i_elem_type) == eWeir) .and. &
                (eta .LE. crest)                         )

            effectivehead = zeroR

        elsewhere  ( (elem2I(:,e2i_elem_type) == eWeir) .and. &
                     (eta .GT. crest)                   .and. &
                     (eta .LT. crown)                         )

            effectivehead = eta - crest
            ! downstream submergence
            maskarray_dn_submergence = ( (dir .GT. zeroI)         .and. &
                                         (faceEtaUp(dnFace) .GT. crest) )
            ! upstream submergance
            maskarray_up_submergence = ( (dir .LT. zeroI)         .and. &
                                         (faceEtaDn(upFace) .GT. crest) )

        elsewhere  ( (elem2I(:,e2i_elem_type) == eWeir) .and. &
                     (eta .GT. crown)                         )

            ! non surcharge weir flow
            effectivehead = crown - crest
            ! mask for surcharged weirs
            is_surcharged = ( (elem2I(:,e2i_elem_type) == eWeir) .and. &
                              (elem2YN(:,e2YN_CanSurcharge))     .and. &
                              (eta .GT. crown)                         )
        endwhere

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


        real(8),      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real(8),      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

        real(8),  intent(inout)   ::  breadth(:), slope(:), area(:), topwidth(:)
        real(8),  intent(inout)   ::  perimeter(:), hyddepth(:), hydradius(:)
        real(8),  intent(in)      ::  depth(:)

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        where ( (elem2I(:,e2i_elem_type) == eWeir)              .and. &
                (elem2I(:,e2i_geometry) == eRectangular )             )

            area        =   depth * breadth
            topwidth    =   breadth
            hyddepth    =   depth
            perimeter   =   breadth + twoR * hyddepth
            hydradius   =   area / perimeter

        elsewhere  ( (elem2I(:,e2i_elem_type) == eWeir)         .and. &
                     (elem2I(:,e2i_geometry) == eTrapezoidal )        )

            area        =   (breadth + slope * depth) * depth
            topwidth    =   breadth + twoR * slope * depth
            hyddepth    =   area / topwidth
            perimeter   =   breadth + twoR * depth * sqrt(oneR + slope ** twoR)
            hydradius   =   area / perimeter

        elsewhere  ( (elem2I(:,e2i_elem_type) == eWeir)         .and. &
                     (elem2I(:,e2i_geometry) == eTriangular )         )

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
        endcontractions, head, breadth, sideslope, crestlength)
        !
        !%  ths subroutine calculates effective creast length for
        !%  trapezoidal and rectangular weir
        !
        character(64) :: subroutine_name = 'weir_effective_crest_length'


        real(8),      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real(8),      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

        real(8),  intent(inout)   :: crestlength(:) 
        real(8),  intent(in)      :: endcontractions(:), head(:)
        real(8),  intent(in)      :: breadth(:), sideslope(:)

        integer :: mm
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        ! effective crest length is used in rectangular and trapezoidal weir flow calculation
        where ( (elem2I(:,e2i_elem_type) == eWeir       ) .and. &
                (elem2I(:,e2i_geometry)  == eRectangular)       )

            crestlength = max(breadth - 0.1 * endcontractions * head, zeroR)

        elsewhere ( (elem2I(:,e2i_elem_type) == eWeir       ) .and. &
                    (elem2I(:,e2i_geometry)  == eTrapezoidal)       )

            crestlength = breadth + twoR * sideslope * head

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

        real(8),      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real(8),      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)

        real(8),    intent(inout) ::  submergenceFactor1(:), submergenceFactor2(:)
        real(8),    intent(in)    ::  crest(:), faceEtaDn(:), faceEtaUp(:)
        integer, intent(in)    ::  upFace(:), dnFace(:)
        logical, intent(in)    ::  maskarray_submergence(:)

        integer :: mm
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% calculate the submergance factor for different weirs according to Villemonte 1974
        where ( (elem2I(:,e2i_weir_elem_type) == eVnotchWeir) .and. &
                (maskarray_submergence)                             )
            ! V-notch weir
            submergenceFactor1 = (oneR - ((faceEtaUp(dnFace) - crest) / (faceEtaDn(upFace) - crest)) &
                ** 2.5) ** 0.385

        elsewhere ( (elem2I(:,e2i_weir_elem_type) == eTrapezoidalWeir) .and. &
                    (maskarray_submergence)                                  )
            ! Trapezoidal weir
            submergenceFactor1 = (oneR - ((faceEtaUp(dnFace) - crest) / (faceEtaDn(upFace) - crest)) &
                ** 2.5) ** 0.385
            submergenceFactor2 = (oneR - ((faceEtaUp(dnFace) - crest) / (faceEtaDn(upFace) - crest)) &
                ** 1.5) ** 0.385

        elsewhere ( (elem2I(:,e2i_weir_elem_type) == eTransverseWeir) .and. &
                    (maskarray_submergence)                                 )
            ! Transverse weir
            submergenceFactor1 = (oneR - ((faceEtaUp(dnFace) - crest) / (faceEtaDn(upFace) - crest)) &
                ** 1.5) ** 0.385

        elsewhere ( (elem2I(:,e2i_weir_elem_type) == eSideFlowWeir) .and. &
                    (maskarray_submergence)                               )
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
        (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2,  &
        velocity2, flow, area, sideslope, cTrig, cRect, dir, effectivehead, &
        crestlength , submergenceFactor1, submergenceFactor2, thiscoef)
        !
        !%  calculate flow through weir given an effective head
        !
        character(64) :: subroutine_name = 'weir_flow'


        real(8),      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real(8),      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
        real(8),              intent(in)      :: thiscoef

        real(8),    intent(inout)   ::  volume2(:), velocity2(:), flow(:)
        real(8),    intent(in)      ::  area(:), sideslope(:), cTrig(:), cRect(:)
        real(8),    intent(in)      ::  submergenceFactor1(:), submergenceFactor2(:)
        real(8),    intent(in)      ::  effectivehead(:), crestlength(:)
        integer, intent(in)      ::  dir(:)

        integer :: mm
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        where ( elem2I(:,e2i_weir_elem_type) == eVnotchWeir )
            ! V-notch weir
            ! fs = submergance_factor
            ! d = direction
            ! Q = d*fs*(Cw1SH^2.5)
            flow      = dir * submergenceFactor1 * cTrig * sideslope * &
                effectivehead ** 2.5
            velocity2 = flow / area
            ! Volume = Q * dt
            volume2   = thiscoef * dt * flow

        elsewhere ( elem2I(:,e2i_weir_elem_type) == eTrapezoidalWeir )
            ! Trapezoidal weir
            ! Q = d*fs1*(Cw1SH^2.5) + d*fs2*(Cw2LH^1.5)
            flow      = dir * submergenceFactor1 * (cTrig * sideslope *    &
                effectivehead ** 2.5) + dir * submergenceFactor2 * (cRect *  &
                crestlength * effectivehead ** 1.5)
            velocity2 = flow / area
            ! Volume = Q * dt
            volume2   = thiscoef * dt * flow

        elsewhere ( elem2I(:,e2i_weir_elem_type) == eTransverseWeir )
            ! Transverse weir
            ! Q = d*fs2*(Cw2LH^1.5)
            flow      = dir * submergenceFactor1 * cRect * &
                crestlength * effectivehead ** 1.5
            velocity2 = flow / area
            ! Volume = Q * dt
            volume2   = thiscoef * dt * flow

        elsewhere ( (elem2I(:,e2i_weir_elem_type) == eSideFlowWeir )  .and. &
                    (dir .LE. zeroR)                                        )
            ! Side flow weir for reverse flow behaves like a Transverse weir
            ! Q = d*fs2*(Cw2LH^1.5)
            flow      = dir * submergenceFactor1 * cRect * &
                crestlength * effectivehead ** 1.5
            velocity2 = flow / area
            ! Volume = Q * dt
            volume2   = thiscoef * dt * flow

        elsewhere ( (elem2I(:,e2i_weir_elem_type) == eSideFlowWeir )  .and. &
                    (dir .GT. zeroR)                                        )
            ! Corrected formula (see Metcalf & Eddy, Inc., Wastewater Engineering, McGraw-Hill, 1972 p. 164).
            flow      = dir * submergenceFactor1 * cRect * &
                (crestlength ** 0.83) * (effectivehead ** 1.67)
            velocity2 = flow / area
            ! Volume = Q * dt
            volume2   = thiscoef * dt * flow
        endwhere

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

    end subroutine weir_flow
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine weir_surcharge_flow &
        (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, volume2, &
        velocity2, flow, crest, crown, eta, area, cOrif, effectivehead,    &
        faceEtaup, faceEtadn, upFace, dnFace, dir, thiscoef, is_surcharged)
        !
        !%  find flow using orifice equation when weir is surcharged
        !
        character(64) :: subroutine_name = 'weir_surcharge_flow'


        real(8),      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        real(8),      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
        real(8),              intent(in)      :: thiscoef

        real(8),    intent(inout) ::  volume2(:), velocity2(:), flow(:), effectivehead(:)
        real(8),    intent(in)    ::  crest(:), crown(:), eta(:), area(:), cOrif(:)
        real(8),    Intent(in)    ::  faceEtaup(:), faceEtadn(:)
        integer, intent(in)    ::  upFace(:), dnFace(:), dir(:)
        logical, intent(in)    ::  is_surcharged(:)

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% surchaged flow calculation
        where ( (elem2I(:,e2i_elem_type) == eWeir) .and.  is_surcharged )
            !% for surcharged flow, head is calculated from the midpoint of the weir opening
            effectivehead = min(( eta - ((crown + crest) / twoR)), &
                dir*(faceEtadn(upFace) - faceEtaup(dnFace)) )

            flow       = dir * cOrif * sqrt(abs(effectivehead))
            !% blend new velocity with old velocity -- needs further checking
            velocity2  = flow / area
            !% Volume is weir flow equation * dt
            !% blend new volume with old volume -- needs further checking
            volume2    =  thiscoef * dt * flow
        endwhere

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine weir_surcharge_flow
    !
    !==========================================================================
    ! END OF MODULE weir
    !==========================================================================
    !
end module weir
