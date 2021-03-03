!
! module element_geometry
!
! Updates the geometry of elements
!
! Later this may be broken into separate modules for channels, pipes
! junctions, etc.
!
!==========================================================================
!
module element_geometry
    !
    use adjustments
    use array_index
    use bc
    use data_keys
    use junction
    use setting_definition
    use globals
    use utility
    use xsect_tables

    implicit none

    private

    public :: element_geometry_update
    public :: element_geometry_branch_fix

    integer :: debuglevel = 0

contains
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine element_geometry_update &
        (elem2R, elem2I, elem2YN, e2r_VolumeColumn, &
        elemMR, elemMI, elemMYN, eMr_VolumeColumn,  &
        faceR, faceI, bcdataDn, bcdataUp, thisTime, method_EtaM, ID, numberPairs, &
        ManningsN, Length, zBottom, xDistance, Breadth, widthDepthData, cellType)
        !
        ! Note that volume is handled as a separate temporary index location
        ! (rather than from the elemR(:,er_Volume) array) because we use
        ! this for the geometry update associated with an RK step where intermediate
        ! storage is used
        !
        character(64) :: subroutine_name = 'element_geometry_update'

        real,      target, intent(in out) :: elem2R(:,:),  elemMR(:,:)
        integer,   target, intent(in out) :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in out) :: elem2YN(:,:), elemMYN(:,:)
        integer,           intent(in)     :: faceI(:,:)
        integer,           intent(in)     :: e2r_VolumeColumn, eMr_VolumeColumn
        real,              intent(in)     :: faceR(:,:), thisTime
        type(bcType),      intent(in out) :: bcdataDn(:), bcdataUp(:)
        integer,           intent(in)     :: method_EtaM

        integer        :: eMr_EtaOld
        real,  pointer :: etaold(:)


        integer, parameter :: ilocaldummy = 0

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

        !%  assign temporary storage for old free surfac elevation on junctions
        eMr_EtaOld = eMr_Temp(next_eMr_temparray)
        etaold  => elemMR(:,eMr_EtaOld)
        next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)
        etaold = elemMR(:,eMr_Eta)

        !% reset values for small volume ahndling
        if (setting%SmallVolume%UseSmallVolumes) then
            elem2YN(:,e2YN_IsSmallVolume) = .false.
            elemMYN(:,eMYN_IsSmallVolume) = .false.
            elem2R(:,e2r_SmallVolumeRatio) = nullvalueR
            elemMR(:,eMr_SmallVolumeRatio) = nullvalueR
        endif

        !% apply any elevation bc and fix the element volume
        call bc_applied_onelement &
            (elem2R, bcdataDn, bcdataUp, thisTime, bc_category_elevation, ilocaldummy)

        !% rectangular geometry
        call geometry_update &
            (elem2R, elem2I, elem2YN, e2r_VolumeColumn, elemMR, elemMI, elemMYN,    &
            eMr_VolumeColumn, faceR, eMr_EtaOld, method_EtaM, ID, numberPairs,      &
            ManningsN, Length, zBottom, xDistance, Breadth, widthDepthData, cellType)
            
        !% HACK -- NEED OTHER GEOMETRY TYPES
        
        !% reset the computed geometry values where volumes are small
        if (setting%SmallVolume%UseSmallVolumes) then
            call adjust_smallvolumes &
                (elem2R, elem2I, elem2YN, e2r_VolumeColumn, &
                elemMR, elemMI, elemMYN, eMr_VolumeColumn    )
        endif

        !% reset the geometry (non-volume) where values are below minimums
        call adjust_for_zero_geometry (elem2R, elem2YN, elemMR, elemMI, elemMYN)

        !%  release the temp array
        etaold = nullvalueR
        nullify(etaold)
        next_eMr_temparray = next_eMr_temparray-1

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine element_geometry_update
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine element_geometry_branch_fix &
        (elemMR, elemMI, faceR, faceI )
        !
        ! called after face update to provide update to junction branches
        !
        character(64) :: subroutine_name = 'element_geometry_branch_fix'

        real,      target, intent(in out)  :: elemMR(:,:)
        real,              intent(in)      :: faceR(:,:)
        integer,           intent(in)      :: elemMI(:,:), faceI(:,:)

        integer                :: eMr_totalarea
        real,      pointer     :: totalarea(:)
        integer,   parameter   :: ilocaldummy = 0

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        eMr_totalarea = eMr_Temp(next_eMr_temparray)
        totalarea     => elemMR(:,eMr_totalarea)
        next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

        !% readjust branches geometry based on new face eta values -----------------
        !%  upstream branches. Note the use of fr_Eta_d for upstream values is logically
        !%  correct as you want the downstream value on the upstream face.
        !%  Note that ilocaldummy=0 because the old eta is not required for the
        !%  second update.
        call rectangular_junction_leg &
            (elemMR, elemMI, faceR, upstream_face_per_elemM, eMi_nfaces_u, eMi_MfaceUp,   &
            eMr_AreaUp, eMr_ZbottomUp, eMr_BreadthScaleUp, &
            eMr_TopwidthUp, eMr_EtaUp, eMr_HydDepthUp, ilocaldummy, fr_Eta_d,2)

        !%  downstream branches
        call rectangular_junction_leg &
            (elemMR, elemMI, faceR, dnstream_face_per_elemM, eMi_nfaces_d, eMi_MfaceDn,   &
            eMr_AreaDn, eMr_ZbottomDn, eMr_BreadthScaleDn, &
            eMr_TopwidthDn, eMr_EtaDn, eMr_HydDepthDn, ilocaldummy, fr_Eta_u,2)

        !% readjust branch flows for updated areas -----------------------------
        !%  get the total outflow branch areas
        totalarea = zeroR ! ensure temporary array is zero

        call junction_branch_sum_areas_by_direction &
            (eMR_totalarea, &
            dnstream_face_per_elemM, eMr_AreaDn, eMi_MfaceDn, eMi_nfaces_d, &
            upstream_face_per_elemM, eMr_AreaUp, eMi_MfaceUp, eMi_nfaces_u, &
            elemMR, elemMI, faceR)

        !%  distribute flow proportionally among outflows
        call junction_branch_velocity_and_flowrate_proportional_to_area &
            (eMR_totalarea, eMR_Flowrate, &
            dnstream_face_per_elemM, eMr_AreaDn, eMi_MfaceDn, eMi_nfaces_d, &
            eMr_FlowrateDn, eMr_VelocityDn, &
            upstream_face_per_elemM, eMr_AreaUp, eMi_MfaceUp, eMi_nfaces_u, &
            eMr_FlowrateUp, eMr_VelocityUp, &
            elemMR, elemMI, faceR)

        !%  enforce maximum velocities in junction branches
        call adjust_junction_branch_velocity_limit (elemMR, elemMI)

        totalarea = nullvalueR
        nullify(totalarea)
        next_eMr_temparray = next_eMr_temparray-1

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine element_geometry_branch_fix
    !
    !==========================================================================
    !
    ! PRIVATE BELOW HERE
    !
    !==========================================================================
    !
    subroutine geometry_update &
        (elem2R, elem2I, elem2YN, e2r_Volume_new, &
        elemMR, elemMI, elemMYN, eMr_Volume_new,  &
        faceR, eMr_EtaOld, method_EtaM, ID, numberPairs, ManningsN, Length,  &
        zBottom, xDistance, Breadth, widthDepthData, cellType)
        !
        ! Note that volume used is in a eTr storage location so that the update
        ! can be used on a temporary volume
        !
        character(64) :: subroutine_name = 'geometry_update'

        real,      intent(in out)  :: elem2R(:,:), elemMR(:,:)
        real,      intent(in)      :: faceR(:,:)
        integer,   intent(in)      :: elem2I(:,:), elemMI(:,:)
        logical,   intent(in out)  :: elem2YN(:,:), elemMYN(:,:)
        integer,   intent(in)      :: e2r_Volume_new, eMr_Volume_new, eMr_EtaOld
        integer,   intent(in)      :: method_EtaM

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

        !%  basic geometry update for rectangular channels and junctions
        !%   geometry for channels
        call channel_pipe_junction &
            (elem2R, elem2I, elem2YN, e2i_geometry, e2i_elem_type, eChannel, &
            e2r_Volume, e2r_Length, e2r_Zbottom, e2r_BreadthScale,           &
            e2r_LeftSlope, e2r_RightSlope, e2r_ParabolaValue, e2r_FullDepth, &
            e2r_FullArea, e2r_Area, e2r_Eta, e2r_Topwidth, e2r_Perimeter,    &
            e2r_Depth, e2r_HydDepth, e2r_HydRadius, e2r_dHdA, e2r_elN,       &
            e2YN_IsSurcharged, e2i_Temp, next_e2i_temparray, e2i_n_temp,     &
            e2r_Temp, next_e2r_temparray, e2r_n_temp, e2YN_Temp,             &
            next_e2YN_temparray, e2YN_n_temp, ID, numberPairs, ManningsN,    &
            Length, zBottom, xDistance, Breadth, widthDepthData, cellType)

        !%   geometry for pipes
        call channel_pipe_junction &
            (elem2R, elem2I, elem2YN, e2i_geometry, e2i_elem_type, ePipe,    &
            e2r_Volume, e2r_Length, e2r_Zbottom, e2r_BreadthScale,           &
            e2r_LeftSlope, e2r_RightSlope, e2r_ParabolaValue, e2r_FullDepth, &
            e2r_FullArea, e2r_Area, e2r_Eta, e2r_Topwidth, e2r_Perimeter,    &
            e2r_Depth, e2r_HydDepth, e2r_HydRadius, e2r_dHdA, e2r_elN,       &
            e2YN_IsSurcharged, e2i_Temp, next_e2i_temparray, e2i_n_temp,     &
            e2r_Temp, next_e2r_temparray, e2r_n_temp, e2YN_Temp,             &
            next_e2YN_temparray, e2YN_n_temp, ID, numberPairs, ManningsN,    &
            Length, zBottom, xDistance, Breadth, widthDepthData, cellType)

        !%   geomety for junctions
        call channel_pipe_junction &
            (elemMR, elemMI, elemMYN, eMi_geometry, eMi_elem_type, eJunctionChannel, &
            eMr_Volume, eMr_Length, eMr_Zbottom, eMr_BreadthScale,           &
            eMr_LeftSlope, eMr_RightSlope, eMr_ParabolaValue, eMr_FullDepth, &
            eMr_FullArea, eMr_Area, eMr_Eta, eMr_Topwidth, eMr_Perimeter,    &
            eMr_Depth, eMr_HydDepth, eMr_HydRadius, eMr_dHdA, eMr_elN,       &
            eMYN_IsSurcharged, eMi_Temp, next_eMi_temparray, eMi_n_temp,     &
            eMr_Temp, next_eMr_temparray, eMr_n_temp, eMYN_Temp,             &
            next_eMYN_temparray, eMYN_n_temp, ID, numberPairs, ManningsN,    &
            Length, zBottom, xDistance, Breadth, widthDepthData, cellType)


        ! HACK: We need to add junction-pipe element and geometry calculation of it


        !% upstream branches
        !% note the fr_Eta_d is used for the upstream face, whose downstream eta is
        !% seen by the upstream junction branch -- so the d is correct while other
        !% values are Up!
        call rectangular_junction_leg &
            (elemMR, elemMI, faceR, upstream_face_per_elemM, eMi_nfaces_u, eMi_MfaceUp,   &
            eMr_AreaUp, eMr_ZbottomUp, eMr_BreadthScaleUp, &
            eMr_TopwidthUp, eMr_EtaUp, eMr_HydDepthUp, eMr_EtaOld, fr_Eta_d, method_EtaM)


        !%  downstream branches
        call rectangular_junction_leg &
            (elemMR, elemMI, faceR, dnstream_face_per_elemM, eMi_nfaces_d, eMi_MfaceDn,   &
            eMr_AreaDn, eMr_ZbottomDn, eMr_BreadthScaleDn, &
            eMr_TopwidthDn, eMr_EtaDn, eMr_HydDepthDn, eMr_EtaOld, fr_Eta_u, method_EtaM)

        !print *, "check point 1222"

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine geometry_update
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine channel_pipe_junction &
        (elemR, elemI, elemYN, ei_geometry, ei_elem_type, elem_type_value,   &
        er_Volume, er_Length, er_Zbottom, er_BreadthScale, er_LeftSlope,     &
        er_RightSlope, er_ParabolaValue, er_FullDepth, er_FullArea, er_Area, &
        er_Eta, er_Topwidth, er_Perimeter, er_Depth, er_HydDepth,            &
        er_HydRadius, er_dHdA, er_elN, eYN_IsSurcharged, ei_Temp,            &
        next_ei_temparray, ei_n_temp, er_Temp, next_er_temparray, er_n_temp, &
        eYN_Temp, next_eYN_temparray, eYN_n_temp, wdID, wdnumberPairs,       &
        wdManningsN, wdLength, wdzBottom, wdxDistance, wdBreadth,            &
        widthDepthData, wdcellType)
        !
        ! computes element geometry for a rectangular channel or a channeljunction
        !
        character(64) :: subroutine_name = 'channel_pipe_junction'

        real,      target,     intent(in out)  :: elemR(:,:)
        logical,   target,     intent(in out)  :: elemYN(:,:)
        
        
        integer,   intent(in)      :: elemI(:,:)
        integer,   intent(in)      :: ei_geometry, ei_elem_type, elem_type_value
        integer,   intent(in)      :: er_Volume, er_Length, er_Zbottom
        integer,   intent(in)      :: er_BreadthScale, er_LeftSlope, er_RightSlope
        integer,   intent(in)      :: er_ParabolaValue, er_FullDepth, er_FullArea
        integer,   intent(in)      :: er_Area, er_Eta, er_Topwidth, er_Perimeter
        integer,   intent(in)      :: er_Depth, er_HydDepth, er_HydRadius, er_dHdA
        integer,   intent(in)      :: er_elN, eYN_IsSurcharged
        integer,   intent(in)      :: er_n_temp, ei_n_temp, eYN_n_temp 
        integer,   intent(in)      :: er_Temp(:), ei_Temp(:), eYN_Temp(:)

        integer,   intent(inout)   :: next_er_temparray, next_ei_temparray, next_eYN_temparray

        integer, target, intent(in out)    :: wdID(:)
        integer, target, intent(in out)    :: wdnumberPairs(:)
        real,    target, intent(in out)    :: wdManningsN(:)
        real,    target, intent(in out)    :: wdLength(:)
        real,    target, intent(in out)    :: wdzBottom(:)
        real,    target, intent(in out)    :: wdxDistance(:)
        real,    target, intent(in out)    :: wdBreadth(:)
        real,    target, intent(in out)    :: widthDepthData(:,:,:)
        type(string), target, intent(in out)   :: wdcellType(:)


        real,    pointer :: volume(:), length(:), zbottom(:), breadth(:)
        real,    pointer :: leftSlope(:), rightSlope(:), parabolaValue(:), fullDepth(:)
        real,    pointer :: fullArea(:), area(:), eta(:), topwidth(:), perimeter(:)
        real,    pointer :: depth(:), hyddepth(:), hydradius(:), dHdA(:), elN(:)
        
        logical, pointer :: isFull(:)

        real, pointer :: widthAtLayerTop(:,:), depthAtLayerTop(:,:), areaThisLayer(:,:)
        real, pointer :: areaTotalBelowThisLayer(:,:), dWidth(:,:)
        real, pointer :: dDepth(:,:), angle(:,:), perimeterBelowThisLayer(:,:)
        !real, dimension(:), allocatable :: area_difference, local_difference

        real, dimension(:), allocatable :: AA, BB, CC, DD, a_diff 
        ! w_d_variables are for solving qudratic function for width-depth geometry
        real, dimension(:), allocatable :: w_d_angle, w_d_widthAtLayerTop, w_d_depthAtLayerTop, w_d_perimeterBelowThisLayer 

        integer :: ii

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        ! inputs
        volume        => elemR(:,er_Volume)
        length        => elemR(:,er_Length)
        zbottom       => elemR(:,er_Zbottom)
        breadth       => elemR(:,er_BreadthScale) 
        leftSlope     => elemR(:,er_LeftSlope)
        rightSlope    => elemR(:,er_RightSlope)
        parabolaValue => elemR(:,er_ParabolaValue)
        fullDepth     => elemR(:,er_FullDepth)
        fullArea      => elemR(:,er_FullArea)
        isFull        => elemYN(:,eYN_IsSurcharged)

        ! outputs
        area       => elemR(:,er_Area)
        eta        => elemR(:,er_Eta)
        topwidth   => elemR(:,er_Topwidth)
        perimeter  => elemR(:,er_Perimeter)
        depth      => elemR(:,er_Depth)
        hyddepth   => elemR(:,er_HydDepth)
        hydradius  => elemR(:,er_HydRadius)
        dHdA       => elemR(:,er_dHdA)
        elN        => elemR(:,er_elN)

        widthAtLayerTop         => widthDepthData (:,:, wd_widthAtLayerTop)
        depthAtLayerTop         => widthDepthData (:,:, wd_depthAtLayerTop)
        areaThisLayer           => widthDepthData (:,:, wd_areaThisLayer)
        areaTotalBelowThisLayer => widthDepthData (:,:, wd_areaTotalBelowThisLayer)
        dWidth                  => widthDepthData (:,:, wd_Dwidth)
        dDepth                  => widthDepthData (:,:, wd_Ddepth)
        angle                   => widthDepthData (:,:, wd_angle)
        perimeterBelowThisLayer => widthDepthData (:,:, wd_perimeterBelowThisLayer)

        ! allocate these arrays for W-D geometry computation
        allocate (a_diff(size(volume,1)), AA(size(volume,1)), &
            BB(size(volume,1)), CC(size(volume,1)), DD(size(volume,1)), source = zeroR)
        !initialize the w_d_ array with all zeroR values
        allocate (w_d_angle(size(volume,1)), w_d_widthAtLayerTop(size(volume,1)), w_d_depthAtLayerTop(size(volume,1)), &
        w_d_perimeterBelowThisLayer(size(volume,1)), source = zeroR)


        !% compute the coefficients required for width-depth XS before we enter the where statement
        call width_depth_quadratic_function(elemI, elemR, ei_geometry, ei_elem_type, elem_type_value, &
            er_Length, er_Area, er_Volume, widthDepthData, w_d_angle, w_d_widthAtLayerTop, w_d_depthAtLayerTop, &
            w_d_perimeterBelowThisLayer, a_diff)

        !%  basic geometry types
        where ( (elemI(:,ei_elem_type) == elem_type_value) .and. &
                (elemI(:,ei_geometry)  == eRectangular)    .and. &
                (isFull .eqv. .false.)                           )

            area        = volume / length
            eta         = zbottom + (area / breadth)
            topwidth    = breadth
            perimeter   = breadth + 2.0 * ( eta - zbottom )
            hyddepth    = area / breadth
            depth       = hyddepth
            hydradius   = area / perimeter
            dHdA        = oneR / breadth
            elN         = hyddepth

        elsewhere ( (elemI(:,ei_elem_type) == elem_type_value) .and. &
                    (elemI(:,ei_geometry)  == eParabolic)      .and. &
                    (isFull .eqv. .false.)                           )
            area        = volume / length
            eta         = zbottom + parabolaValue ** oneThirdR &
                * (threefourthR * area) ** twothirdR
            hyddepth    = parabolaValue ** oneThirdR * (threefourthR * area) ** twothirdR
            depth       = (threeR / twoR) * hyddepth
            topwidth    = twoR * sqrt(depth/parabolaValue)
            perimeter   = onehalfR * topwidth &
                *( &
                sqrt( oneR + (fourR * depth/topwidth)**twoR )  &
                + (topwidth/fourR * depth) &
                *log &
                ( &
                fourR * depth/topwidth  &
                + sqrt( oneR + (fourR * depth/topwidth)**twoR ) &
                )  &
                )
            hydradius   = area / perimeter
            dHdA        = (parabolaValue ** oneThirdR) * (threefourthR ** twothirdR) &
                * twothirdR * (area) ** (-oneThirdR)
            elN         = hyddepth

        elsewhere ( (elemI(:,ei_elem_type) == elem_type_value) .and. &
                    (elemI(:,ei_geometry)  == eTrapezoidal)    .and. &
                    (isFull .eqv. .false.)                           )  

            area        = volume / length
            depth       = - onehalfR * (breadth/(onehalfR*(leftSlope + rightSlope)) &
                - sqrt((breadth/(onehalfR*(leftSlope + rightSlope))) ** twoR &
                + fourR * area/(onehalfR*(leftSlope + rightSlope))))

            eta         = zbottom + depth
            topwidth    = breadth + depth * (leftSlope + rightSlope)
            hyddepth    = area / topwidth
            perimeter   = breadth + depth &
                * (sqrt(oneR + leftSlope**twoR ) &
                + sqrt(oneR + rightSlope**twoR))
            hydradius   = area / perimeter
            dHdA        = oneR / topwidth
            elN         = hyddepth

        elsewhere ( (elemI(:,ei_elem_type) == elem_type_value) .and. &
                    (elemI(:,ei_geometry)  == eTriangular)     .and. &
                    (isFull .eqv. .false.)                           )  

            area        = volume / length
            depth       = sqrt(abs(area/(onehalfR*(leftSlope + rightSlope))))
            hyddepth    = onehalfR * depth
            eta         = zbottom + hyddepth
            topwidth    = (leftSlope + rightSlope) * depth
            perimeter   = depth * (sqrt(oneR + leftSlope**twoR) + sqrt(oneR + rightSlope**twoR))
            hydradius   = area / perimeter
            dHdA        = oneR / topwidth
            elN         = hyddepth

        elsewhere ( (elemI(:,ei_elem_type) == elem_type_value) .and. &
                    (elemI(:,ei_geometry)  == eWidthDepth)     .and. &
                    (isFull .eqv. .false.)                           )

            area        = volume / length
            AA          = oneR / tan( w_d_angle )
            BB          = w_d_widthAtLayerTop
            CC          = negoneR * a_diff
            DD          = ( negoneR * BB + sqrt( BB ** twoR - 4.0 * AA * CC) ) / ( twoR * AA )
            depth       = DD + w_d_depthAtLayerTop
            topwidth    = w_d_widthAtLayerTop + twoR * AA * DD
            hyddepth    = area / topwidth
            eta         = zbottom + depth
            perimeter   = w_d_perimeterBelowThisLayer + (twoR * DD / sin(w_d_angle))
            hydradius   = area / perimeter
            dHdA        = oneR / topwidth
            elN         = hyddepth 

        endwhere

        !%  specialized geometry types
        call circular_geometry &
            (elemI, elemR, elemYN, ei_geometry, ei_elem_type, elem_type_value,   &
            volume, length, zbottom, breadth, fulldepth, fullarea, area, eta,    &
            perimeter, depth, hyddepth, hydradius, topwidth, dHdA, elN, ei_Temp, &
            next_ei_temparray, ei_n_temp, er_Temp, next_er_temparray, er_n_temp, &
            eYN_Temp, next_eYN_temparray, eYN_n_temp, isFull)

        !% HACK: other specialized geometry types are needed

        !% geometry calculation for all surcharged elements
        call surcharged_element_geometry &
            (elemI, elemR, elemYN, ei_geometry, ei_elem_type, elem_type_value,   &
            volume, length, zbottom, breadth, fulldepth, fullarea, area, eta,    &
            perimeter, depth, hyddepth, hydradius, topwidth, dHdA, elN, isFull)
        
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine channel_pipe_junction
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine circular_geometry &
        (elemI, elemR, elemYN, ei_geometry, ei_elem_type, elem_type_value,   &
        volume, length, zbottom, breadth, fulldepth, fullarea, area, eta,    &
        perimeter, depth, hyddepth, hydradius, topwidth, dHdA, elN, ei_Temp, &
        next_ei_temparray, ei_n_temp, er_Temp, next_er_temparray, er_n_temp, &
        eYN_Temp, next_eYN_temparray, eYN_n_temp, isFull) 
        !
        ! Calculate circular geometry 
        !
        character(64) :: subroutine_name = 'circular_geometry'

        real,      target,     intent(inout)  :: elemR(:,:)
        logical,   target,     intent(inout)  :: elemYN(:,:)

        integer,   intent(in)       :: elemI(:,:)
        integer,   intent(in)       :: ei_geometry, ei_elem_type, elem_type_value
        real,      intent(in)       :: volume(:), length(:), zbottom(:), breadth(:)
        real,      intent(in)       :: fulldepth(:), fullarea(:)
        integer,   intent(in)       :: ei_n_temp, er_n_temp, eYN_n_temp
        integer,   intent(in)       :: ei_Temp(:), er_Temp(:), eYN_Temp(:)
        logical,   intent(in)       :: isFull(:)

        real,      intent(inout)    :: area(:), eta(:), perimeter(:), depth(:)
        real,      intent(inout)    :: hyddepth(:), hydradius(:), topwidth(:)
        real,      intent(inout)    :: dHdA(:), elN(:)
        integer,   intent(inout)    :: next_ei_temparray, next_er_temparray, next_eYN_temparray

        real,       pointer         :: AoverAfull(:), YoverYfull(:) 
        logical,    pointer         :: maskarray(:)
        real                        :: af, bf, cf
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% temporary pointer allocation fo circular for geometry update
        AoverAfull => elemR(:,er_Temp(next_er_temparray))
        next_er_temparray = utility_advance_temp_array (next_er_temparray,er_n_temp)

        YoverYfull => elemR(:,er_Temp(next_er_temparray))
        next_er_temparray = utility_advance_temp_array (next_er_temparray,er_n_temp)
        !% temporary mask to find open circular elements
        maskarray  => elemYN(:,eYN_Temp(next_eYN_temparray))
        next_eYN_temparray = utility_advance_temp_array (next_eYN_temparray,eYN_n_temp)

        maskarray  = nullvalueL
        AoverAfull = nullvalueR
        YoverYfull = nullvalueR

        maskarray = ( (elemI(:,ei_elem_type) == elem_type_value) .and. &
                      (elemI(:,ei_geometry)  == eCircular)       .and. &
                      (isFull .eqv. .false.)  )

        !% find Area and A/Afull for circular elements
        where (maskarray)
            area = volume / length
            AoverAfull = area / fullarea
        endwhere

        !% find normalized depth using A/Afull from lookup tables
        call table_lookup_mask &
            (elemI, elemR, YoverYfull, AoverAfull, YCirc, NYCirc, maskarray, &
            ei_Temp, next_ei_temparray, ei_n_temp)

        !% find normalized topwidth using Y/Yfull from lookup tables
        !% normalized topwidth (W/Wmax) is saved in the topwidth column
        call table_lookup_mask &
            (elemI, elemR, topwidth, YoverYfull, WCirc, NWCirc, maskarray, &
            ei_Temp, next_ei_temparray, ei_n_temp)

        !% find normalized hydraulic radius using Y/Yfull from lookup tables
        !% normalized hydraulic radius (R/Rmax) is saved in the hydradius column
        call table_lookup_mask &
            (elemI, elemR, hydradius, YoverYfull, RCirc, NRCirc, maskarray, &
            ei_Temp, next_ei_temparray, ei_n_temp)

        !% find all the circular geometric properties from normalized values
        where (maskarray)
           depth      = fulldepth * YoverYfull 
           eta        = zbottom + depth
           topwidth   = fulldepth * topwidth
           hydradius  = onefourthR * fulldepth * hydradius
           perimeter  = area / hydradius
        endwhere

        !% HACK: DONT LET TOPWIDTH FALL AT ZERO
        where (topwidth .LE. setting%SmallVolume%MinimumTopwidth)
            topwidth = setting%SmallVolume%MinimumTopwidth
        endwhere

        ! find modified hydraulic depth and for circular pipe from pipeAC Hodges 2020
        where ( (maskarray) .and. (AoverAfull .GT. onehalfR) )
            hyddepth   = eta - zbottom + (onehalfR*fulldepth) * (onefourthR*pi - oneR)
            elN        = hyddepth
        elsewhere ( (maskarray) .and. (AoverAfull .LE. onehalfR) .and. (AoverAfull .GT. zeroR) )
            hyddepth   = max(area / topwidth, zeroR)
            elN        = hyddepth
        endwhere
        ! if (count(maskarray)>zeroI) then
        !     print*,'..............................................'
        !     print*, subroutine_name
        !     print*, '..............................................'
        !     print*, eta, 'eta'
        !     print*
        !     print*, depth, 'depth'
        !     print*
        !     print*, topwidth, 'topwidth'
        !     print*
        !     print*, area, 'area'
        !     print*
        !     print*, elN, 'eln'
        !     print*,'...............................................'
        ! endif

        !% constants for circular dHdA (pipeAC2020)
        af = 1.29
        bf = 0.66
        cf = 0.34
        !% dHdA calculation
        where( (maskarray) .and. (AoverAfull .LT. onehalfR) ) 
            dHdA = (af * bf / ( (pi**bf) * ((fulldepth/twoR)**(twoR*bf - oneR)) ) ) & 
                    * area**(bf - oneR) + cf / (pi * (fulldepth/twoR))

        elsewhere( (maskarray) .and. (AoverAfull .GT. onehalfR) .and. &
                   (AoverAfull .LT. oneR) )
            dHdA = (af * bf / ( (pi**bf) * ((fulldepth/twoR)**(twoR*bf - oneR))) ) &
                    * (pi * ((fulldepth/twoR)**twoR) - area)**(bf - oneR)          &
                    + cf / (pi * (fulldepth/twoR))
        endwhere

        !% nullify temporary array
        AoverAfull = nullvalueR
        YoverYfull = nullvalueR

        nullify(AoverAfull, YoverYfull, maskarray)
        next_er_temparray  = next_er_temparray  - 2
        next_eYN_temparray = next_eYN_temparray - 1

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine circular_geometry
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine surcharged_element_geometry &
        (elemI, elemR, elemYN, ei_geometry, ei_elem_type, elem_type_value,   &
        volume, length, zbottom, breadth, fulldepth, fullarea, area, eta,    &
        perimeter, depth, hyddepth, hydradius, topwidth, dHdA, elN, isFull)
        !
        ! Calculates geometry for surcharged elements
        ! Surcharged elements has eta updated. Thus, these elements need special 
        ! Geometry handler.
        !
        character(64) :: subroutine_name = 'surcharged_element_geometry'

        real,      target,     intent(inout)  :: elemR(:,:)
        logical,   target,     intent(inout)  :: elemYN(:,:)

        integer,   intent(in)       :: elemI(:,:)
        integer,   intent(in)       :: ei_geometry, ei_elem_type, elem_type_value
        real,      intent(in)       :: volume(:), length(:), zbottom(:), breadth(:)
        real,      intent(in)       :: eta(:), fulldepth(:), fullarea(:)
        logical,   intent(in)       :: isFull(:)

        real,      intent(inout)    :: area(:), perimeter(:), depth(:)
        real,      intent(inout)    :: hyddepth(:), hydradius(:), topwidth(:)
        real,      intent(inout)    :: dHdA(:), elN(:)
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% surcharged elements already has eta updated. Only need to update
        !% other geometry values
        where ( (elemI(:,ei_elem_type) == elem_type_value) .and. &
                (elemI(:,ei_geometry)  == eRectangular)    .and. &
                (isFull) )
            area      = fullArea
            depth     = fullDepth 
            !%  minimum topwidth is equal to the breadth for rectangular full pipe
            topwidth  = breadth
            hyddepth  = eta - zbottom
            perimeter = twoR * (breadth + fulldepth)
            hydradius = area / perimeter
            dHdA      = zeroR
            elN       = hyddepth

        elsewhere ( (elemI(:,ei_elem_type) == elem_type_value) .and. &
                    (elemI(:,ei_geometry)  == eCircular)       .and. &
                    (isFull) )
            area      = fullArea
            depth     = fullDepth 
            !%  minimum topwidth at 5% of radius
            topwidth  = 0.05 * onehalfR * fullDepth
            hyddepth  = eta - zbottom + (onehalfR*fulldepth) * (onefourthR*pi - oneR)
            perimeter = pi * fulldepth
            hydradius = onefourthR * fulldepth
            dHdA      = zeroR 
            elN       = hyddepth

        endwhere

        ! if (count(( (elemI(:,ei_elem_type) == elem_type_value) .and. &
        !             (elemI(:,ei_geometry)  == eCircular)       .and. &
        !             (isFull) ))>0) then
        !     print*,'..............................................'
        !     print*, subroutine_name
        !     print*, '..............................................'
        !     print*, eta, 'eta'
        !     print*
        !     print*, depth, 'depth'
        !     print*
        !     print*, topwidth, 'topwidth'
        !     print*
        !     print*, area, 'area'
        !     print*
        !     print*, elN, 'eln'
        !     print*
        !     print*, 'pipe surcharged at element geom: press return to continue'
        !     print*,'...............................................'
        ! !     read(*,*)
        ! endif
        !% HACK: other special geometry types are needed

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine surcharged_element_geometry
    !
    !==========================================================================
    !==========================================================================
    !
    ! subroutine pipe_or_junction &
    !     (elemR, elemI, elemYN, ei_geometry, ei_elem_type, elem_type_value,   &
    !     ei_solver, er_Length, er_Zbottom, er_BreadthScale, er_Topwidth,      &
    !     er_Area, er_Eta, er_Perimeter, er_Depth, er_HydDepth, er_HydRadius,  &
    !     er_Radius, er_Volume, er_FullDepth, er_FullArea, er_Zcrown, er_dHdA, &
    !     er_elN, eYN_IsSurcharged, er_Temp, next_er_temparray, er_n_temp,     &
    !     ei_Temp, next_ei_temparray, ei_n_temp, eYN_Temp, next_eYN_temparray, &
    !     eYN_n_temp)
    !     !
    !     !% computes geometry for pipe and junction pipe
    !     !
    !     character(64) :: subroutine_name = 'pipe_or_junction'

    !     real,      target,     intent(inout)  :: elemR(:,:)
    !     logical,   target,     intent(inout)  :: elemYN(:,:)
    !     integer,   target,     intent(in)     :: elemI(:,:)

    !     integer,   intent(in)      :: ei_geometry, ei_elem_type, elem_type_value, ei_solver
    !     integer,   intent(in)      :: er_Length, er_Zbottom, er_BreadthScale, er_Perimeter 
    !     integer,   intent(in)      :: er_Topwidth, er_Area, er_Eta, er_dHdA, er_elN
    !     integer,   intent(in)      :: er_Depth, er_HydDepth, er_HydRadius, er_Radius
    !     integer,   intent(in)      :: er_Volume, er_FullDepth, er_FullArea, er_Zcrown
    !     integer,   intent(in)      :: er_n_temp, ei_n_temp, eYN_n_temp
    !     integer,   intent(in)      :: eYN_IsSurcharged 
    !     integer,   intent(in)      :: er_Temp(:), ei_Temp(:), eYN_Temp(:)

    !     integer,   intent(inout)   :: next_er_temparray, next_ei_temparray, next_eYN_temparray

    !     real,    pointer    :: volume(:), length(:), zbottom(:), breadth(:), hyddepth(:)
    !     real,    pointer    :: area(:),  eta(:), perimeter(:), depth(:), dHdA(:), elN(:)
    !     real,    pointer    :: hydradius(:), topwidth(:), fulldepth(:), zcrown(:)
    !     real,    pointer    :: radius(:), AoverAfull(:), YoverYfull(:), fullarea(:)
    !     integer, pointer    :: solver(:), geometry(:)
    !     logical, pointer    :: isfull(:), maskarray(:)

    !     integer :: ii

    !     !--------------------------------------------------------------------------
    !     if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

    !     !% inputs
        
    !     length     => elemR(:,er_Length)
    !     zbottom    => elemR(:,er_Zbottom)
    !     breadth    => elemR(:,er_BreadthScale)
    !     fulldepth  => elemR(:,er_FullDepth)
    !     zcrown     => elemR(:,er_Zcrown)
    !     fullarea   => elemR(:,er_FullArea)
    !     radius     => elemR(:,er_Radius)
    !     solver     => elemI(:,ei_solver)
    !     geometry   => elemI(:,ei_geometry)
    !     isfull     => elemYN(:,eYN_IsSurcharged)

    !     !% outputs
    !     area       => elemR(:,er_Area)
    !     eta        => elemR(:,er_Eta)
    !     perimeter  => elemR(:,er_Perimeter)
    !     depth      => elemR(:,er_Depth)
    !     hyddepth   => elemR(:,er_HydDepth)
    !     hydradius  => elemR(:,er_HydRadius)
    !     topwidth   => elemR(:,er_Topwidth)
    !     dHdA       => elemR(:,er_dHdA)
    !     elN        => elemR(:,er_elN)
    !     volume     => elemR(:,er_Volume)

    !     !% temporary array for geometry update
    !     AoverAfull => elemR(:,er_Temp(next_er_temparray))
    !     next_er_temparray = utility_advance_temp_array (next_er_temparray,er_n_temp)

    !     YoverYfull => elemR(:,er_Temp(next_er_temparray))
    !     next_er_temparray = utility_advance_temp_array (next_er_temparray,er_n_temp)

    !     maskarray  => elemYN(:,eYN_Temp(next_eYN_temparray))
    !     next_eYN_temparray = utility_advance_temp_array (next_eYN_temparray,eYN_n_temp)

    !     AoverAfull = nullvalueR
    !     YoverYfull = nullvalueR
    !     maskarray  = nullvalueL

    !     !% mask for circular pipes/junction pipe
    !     maskarray = ( (elemI(:,ei_elem_type) == elem_type_value) )

    !     !% For open pipe, set eta to zero for later error detection
    !     where (maskarray .and. (isfull .eqv. .false.)) 
    !           eta = zeroR
    !     endwhere
    !     ! print*,'----before-----'
    !     ! print*, area, 'area'
    !     ! print*, fullarea, 'fullarea'
    !     !% Open circular pipes that become full
    !     call open_pipe_transition_to_full &
    !         (elemI, elemR, elemYN, volume, length, zbottom, breadth, fulldepth,   &
    !         fullarea, zcrown, radius, area, eta, perimeter, depth, hyddepth,      &
    !         hydradius, topwidth, AoverAfull, YoverYfull, solver, geometry,        &
    !         ei_Temp, next_ei_temparray, ei_n_temp, eYN_Temp, next_eYN_temparray,  &
    !         eYN_n_temp, isfull, maskarray)

    !     !% Set the full pipe area
    !     !% These cells already have eta directly updated from the time-stepping.
    !     where (maskarray .and. (isfull .eqv. .true.)) 
    !         area   = fullarea
    !         volume = area * length
    !         depth  = fulldepth
    !         AoverAfull = oneR
    !         YoverYfull = oneR  
    !     endwhere

    !     !% NOTE: at this point, isfull will consist of those pipes that
    !     !% were designated as "full" coming into this routine and those that
    !     !% were detected as becoming full due to the increase in area. However
    !     !% the array ALSO is true (incorrectly) for those pipes that dropped
    !     !% from full to open during the prior time step. These must be handled
    !     !% separately.

    !     !% Full pipes that become open
    !     call full_pipe_transition_to_open &
    !         (elemI, elemR, elemYN, volume, length, zbottom, breadth, fulldepth,   &
    !         fullarea, zcrown, radius, area, eta, perimeter, depth, hyddepth,      &
    !         hydradius, topwidth, AoverAfull, YoverYfull, solver, geometry,        &
    !         ei_Temp, next_ei_temparray, ei_n_temp, eYN_Temp, next_eYN_temparray,  &
    !         eYN_n_temp, isfull, maskarray)

    !     !% Open pipes
    !     !% These has areas/volumes updated and needs eta to be computed
    !     call open_pipe &
    !         (elemI, elemR, elemYN, volume, length, zbottom, breadth, fulldepth,   &
    !         fullarea, zcrown, radius, area, eta, perimeter, depth, hyddepth,      &
    !         hydradius, topwidth, AoverAfull, YoverYfull, solver, geometry,        &
    !         ei_Temp, next_ei_temparray, ei_n_temp, eYN_Temp, next_eYN_temparray,  &
    !         eYN_n_temp, isfull, maskarray) 

    !     !% isfull reset
    !     !% set formerly full pipes that have become open to open
    !     where (maskarray .and. (isfull .eqv. .true.) .and. (eta < zcrown))
    !         isfull = .false.
    !     endwhere

    !     ! Now eta and area have been updated for all cases
    !     call pipe_additional_geometric_properties &
    !         (elemI, elemR, elemYN, volume, length, zbottom, breadth, fulldepth,   &
    !         fullarea, zcrown, radius, area, eta, perimeter, depth, hyddepth,      &
    !         hydradius, topwidth, dHdA, elN, AoverAfull, YoverYfull, solver,       &
    !         geometry, ei_Temp, next_ei_temparray, ei_n_temp, eYN_Temp,            &
    !         next_eYN_temparray, eYN_n_temp, isfull, maskarray)
    !     !% release temporary arrays
    !     AoverAfull = nullvalueR
    !     YoverYfull = nullvalueR
    !     maskarray  = nullvalueL
    !     ! print*, '-----after------'
    !     ! print*, volume, 'volume'
    !     ! print*, area, 'area'
    !     ! print*, depth, 'depth'
    !     ! print*, hyddepth, 'hyddepth'
    !     !% nullify temporary array
    !     nullify(AoverAfull, YoverYfull, maskarray)
    !     next_er_temparray  = next_er_temparray  - 2
    !     next_eYN_temparray = next_eYN_temparray - 1

    !     if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    ! end subroutine pipe_or_junction
    ! !
    ! !==========================================================================
    ! !==========================================================================
    ! !
    ! subroutine open_pipe_transition_to_full &
    !     (elemI, elemR, elemYN, volume, length, zbottom, breadth, fulldepth,   &
    !     fullarea, zcrown, radius, area, eta, perimeter, depth, hyddepth,      &
    !     hydradius, topwidth, AoverAfull, YoverYfull, solver, geometry,        &
    !     ei_Temp, next_ei_temparray, ei_n_temp, eYN_Temp, next_eYN_temparray,  &
    !     eYN_n_temp, isfull, maskarray)
    !     !
    !     ! OPEN PIPES THAT BECOME FULL ============
    !     ! Detect transition from open to full pipe
    !     ! area is advanced for open pipe. To eta is needed to be calculated here
    !     !
    !     character(64) :: subroutine_name = 'open_pipe_transition_to_full'

    !     real,      target,     intent(inout)  :: elemR(:,:)
    !     logical,   target,     intent(inout)  :: elemYN(:,:)

    !     integer,   intent(in)       :: elemI(:,:)
    !     real,      intent(in)       :: length(:), zbottom(:), breadth(:)
    !     real,      intent(in)       :: fulldepth(:), fullarea(:), zcrown(:), radius(:)
    !     integer,   intent(in)       :: ei_n_temp, eYN_n_temp
    !     integer,   intent(in)       :: solver(:), geometry(:), ei_Temp(:), eYN_Temp(:)
    !     logical,   intent(in)       :: maskarray(:)

    !     real,      intent(inout)    :: volume(:), area(:), eta(:), perimeter(:), depth(:)
    !     real,      intent(inout)    :: hyddepth(:), hydradius(:), topwidth(:)
    !     real,      intent(inout)    :: AoverAfull(:), YoverYfull(:)
    !     integer,   intent(inout)    :: next_ei_temparray, next_eYN_temparray 
    !     logical,   intent(inout)    :: isfull(:)

    !     logical,    pointer         :: maskarray_pipe_transition(:)

    !     !--------------------------------------------------------------------------
    !     if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

    !     !% temporary mask space for finding transiotional pipes
    !     maskarray_pipe_transition  => elemYN(:,eYN_Temp(next_eYN_temparray))
    !     next_eYN_temparray = utility_advance_temp_array (next_eYN_temparray,eYN_n_temp)

    !     maskarray_pipe_transition  = nullvalueL

    !     !% If a pipe become surcharged inbetween an SVE step, the solver will solve for volume
    !     !% this makes sure that area is updated beofre finding the transitional pipes. 
    !     where (maskarray .and. (isfull .eqv. .false.) .and. (solver == SVE) )
    !         area = volume / length
    !     endwhere

    !     maskarray_pipe_transition = (maskarray .and. (isfull .eqv. .false.) &
    !                                            .and. (area .GE. fullarea) )
    !     where (maskarray_pipe_transition) 
    !         ! Set head above the pipe crown based on the excess area in
    !         ! from the time advance divided by pipe width
    !         eta  = zcrown + (area - fullarea) / breadth
    !         area = fullarea
    !         AoverAfull = oneR
    !         YoverYfull = oneR
    !         volume = area * length
    !         isfull = .true.
    !     endwhere

    !     ! print*, 'geometry debug for ', trim(subroutine_name)
    !     ! print*, 'is full    ', isfull
    !     ! print*, 'Area       ', area
    !     ! print*, 'AoverAfull ', AoverAfull
    !     ! print*, 'depth      ', depth
    !     ! print*, 'eta        ', eta


    !     maskarray_pipe_transition  = nullvalueL
    !     !% nullify temporary array
    !     nullify(maskarray_pipe_transition)
    !     next_eYN_temparray = next_eYN_temparray - 1

    !     if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    ! end subroutine open_pipe_transition_to_full
    ! !
    ! !==========================================================================
    ! !==========================================================================
    ! !
    ! subroutine full_pipe_transition_to_open &
    !     (elemI, elemR, elemYN, volume, length, zbottom, breadth, fulldepth,   &
    !     fullarea, zcrown, radius, area, eta, perimeter, depth, hyddepth,      &
    !     hydradius, topwidth, AoverAfull, YoverYfull, solver, geometry,        &
    !     ei_Temp, next_ei_temparray, ei_n_temp, eYN_Temp, next_eYN_temparray,  &
    !     eYN_n_temp, isfull, maskarray)
    !     !
    !     ! FULL PIPES THAT TRANSITION TO OPEN ==================
    !     ! These have eta and need area computed
    !     ! Note that these are not re-designated as open until after all
    !     ! the eta and are computations are complete
    !     ! Detect full pipe that have become open
    !     !
    !     character(64) :: subroutine_name = 'full_pipe_transition_to_open'

    !     real,      target,     intent(inout)  :: elemR(:,:)
    !     logical,   target,     intent(inout)  :: elemYN(:,:)

    !     integer,   intent(in)       :: elemI(:,:)
    !     real,      intent(in)       :: length(:), zbottom(:), breadth(:)
    !     real,      intent(in)       :: fulldepth(:), fullarea(:), zcrown(:), radius(:)
    !     integer,   intent(in)       :: ei_n_temp, eYN_n_temp
    !     integer,   intent(in)       :: solver(:), geometry(:), ei_Temp(:), eYN_Temp(:)
    !     logical,   intent(in)       :: maskarray(:), isfull(:)

    !     real,      intent(inout)    :: volume(:), area(:), eta(:), perimeter(:), depth(:)
    !     real,      intent(inout)    :: hyddepth(:), hydradius(:), topwidth(:)
    !     real,      intent(inout)    :: AoverAfull(:), YoverYfull(:)
    !     integer,   intent(inout)    :: next_ei_temparray, next_eYN_temparray 

    !     logical,    pointer         :: maskarray_pipe_transition(:)
    !     logical,    pointer         :: maskarray_pipe_geometry(:)

    !     !--------------------------------------------------------------------------
    !     if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

    !     !% temporary mask space for finding transiotional pipes
    !     maskarray_pipe_transition  => elemYN(:,eYN_Temp(next_eYN_temparray))
    !     next_eYN_temparray = utility_advance_temp_array (next_eYN_temparray,eYN_n_temp)

    !     !% temporary mask space for geometry types that needs table interpolation
    !     maskarray_pipe_geometry  => elemYN(:,eYN_Temp(next_eYN_temparray))
    !     next_eYN_temparray = utility_advance_temp_array (next_eYN_temparray,eYN_n_temp)

    !     maskarray_pipe_transition = nullvalueL
    !     maskarray_pipe_geometry   = nullvalueL

    !     maskarray_pipe_transition = (maskarray .and. (isfull .eqv. .true.) &
    !                                            .and. (eta < zcrown) )
    !     where (maskarray_pipe_transition) 
    !         depth = eta - zbottom
    !         YoverYfull = depth / fulldepth
    !     endwhere

    !     !% set up the mask for circular geometry type as an input for table_lookup_mask subroutine
    !     maskarray_pipe_geometry = ( (maskarray_pipe_transition) .and. (geometry == eCircular) )

    !     !% get the normalized area from lookup table from full to open transitional pipe
    !     call table_lookup_mask &
    !         (elemI, elemR, AoverAfull, YoverYfull, ACirc, NACirc, maskarray_pipe_geometry, &
    !         ei_Temp, next_ei_temparray, ei_n_temp)

    !     !% release the temporary maskarray_pipe_geometry for other special geometry types
    !     maskarray_pipe_geometry = nullvalueL

    !     !% HACK: Other special pipe geometry types are needed
           
    !     where ( (maskarray_pipe_transition) .and. (geometry == eCircular) )
    !         !% get the pipe area by multiplying the normalized area with full area for circular pipe
    !         area   = AoverAfull * fullarea
    !         volume = area * length
    !     elsewhere ( (maskarray_pipe_transition) .and. (geometry == eRectangular) )
    !         area   = depth * breadth
    !         volume = area  * length
    !         AoverAfull = area / fullarea
    !     endwhere

    !     ! print*, 'geometry debug for ', trim(subroutine_name)
    !     ! print*, 'is full    ', isfull
    !     ! print*, 'Area       ', area
    !     ! print*, 'AoverAfull ', AoverAfull
    !     ! print*, 'depth      ', depth
    !     ! print*, 'eta        ', eta

    !     maskarray_pipe_transition = nullvalueL
    !     maskarray_pipe_geometry   = nullvalueL
    !     !% nullify temporary array
    !     nullify(maskarray_pipe_transition, maskarray_pipe_geometry)
    !     next_eYN_temparray = next_eYN_temparray - 2

    !     if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    ! end subroutine full_pipe_transition_to_open
    ! !
    ! !==========================================================================
    ! !==========================================================================
    ! !
    ! subroutine open_pipe &
    !     (elemI, elemR, elemYN, volume, length, zbottom, breadth, fulldepth,   &
    !     fullarea, zcrown, radius, area, eta, perimeter, depth, hyddepth,      &
    !     hydradius, topwidth, AoverAfull, YoverYfull, solver, geometry,        &
    !     ei_Temp, next_ei_temparray, ei_n_temp, eYN_Temp, next_eYN_temparray,  &
    !     eYN_n_temp, isfull, maskarray)   
    !     !
    !     ! this subroutine gets the value of circular pipe/junction-pipe geometry
    !     ! solved using the SVE/AC solver
    !     !
    !     character(64) :: subroutine_name = 'open_pipe'

    !     real,      target,     intent(inout)  :: elemR(:,:)
    !     logical,   target,     intent(inout)  :: elemYN(:,:)

    !     integer,   intent(in)       :: elemI(:,:)
    !     real,      intent(in)       :: length(:), zbottom(:), breadth(:)
    !     real,      intent(in)       :: fulldepth(:), fullarea(:), zcrown(:), radius(:)
    !     integer,   intent(in)       :: ei_n_temp, eYN_n_temp
    !     integer,   intent(in)       :: solver(:), geometry(:), ei_Temp(:), eYN_Temp(:)
    !     logical,   intent(in)       :: maskarray(:), isfull(:)

    !     real,      intent(inout)    :: volume(:), area(:), eta(:), perimeter(:), depth(:)
    !     real,      intent(inout)    :: hyddepth(:), hydradius(:), topwidth(:)
    !     real,      intent(inout)    :: AoverAfull(:), YoverYfull(:)
    !     integer,   intent(inout)    :: next_ei_temparray, next_eYN_temparray 

    !     logical,    pointer         :: maskarray_open_pipe(:)
    !     logical,    pointer         :: maskarray_pipe_geometry(:)

    !     !--------------------------------------------------------------------------
    !     if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

    !     !% temporary mask space for finding transiotional pipes
    !     maskarray_open_pipe  => elemYN(:,eYN_Temp(next_eYN_temparray))
    !     next_eYN_temparray = utility_advance_temp_array (next_eYN_temparray,eYN_n_temp)

    !     !% temporary mask space for geometry types that needs table interpolation
    !     maskarray_pipe_geometry  => elemYN(:,eYN_Temp(next_eYN_temparray))
    !     next_eYN_temparray = utility_advance_temp_array (next_eYN_temparray,eYN_n_temp)

    !     maskarray_open_pipe = nullvalueL
    !     maskarray_pipe_geometry = nullvalueL

    !     maskarray_open_pipe = (maskarray .and. (isfull .eqv. .false.))

    !     ! SVE solver solves for volume. Get the area from volume
    !     where (maskarray_open_pipe .and. (solver == SVE))
    !         area = volume / length
    !     endwhere

    !     ! get the normalized area to solve for other geometries
    !     where (maskarray_open_pipe)
    !         AoverAfull = area / fullarea
    !     endwhere
        
    !     !% set up the mask for circular geometry type as an input for table_lookup_mask subroutine
    !     maskarray_pipe_geometry = ( (maskarray_open_pipe) .and. (geometry == eCircular) )

    !     ! get normalized depth from the lookup table
    !     call table_lookup_mask &
    !         (elemI, elemR, YoverYfull, AoverAfull, YCirc, NYCirc, maskarray_pipe_geometry, &
    !         ei_Temp, next_ei_temparray, ei_n_temp)

    !     !% release the temporary maskarray_pipe_geometry for other special geometry types
    !     maskarray_pipe_geometry = nullvalueL

    !     !% HACK: Other special pipe geometry types are needed

    !     where (maskarray_open_pipe .and. (geometry == eCircular))
    !         ! get the depth by multiplying the normalized depth with fulldepth
    !         depth = fulldepth * YoverYfull
    !         eta   = zbottom + depth
    !     elsewhere (maskarray_open_pipe .and. (geometry == eRectangular))
    !         depth = area / breadth
    !         eta   = zBottom + depth
    !     endwhere 

    !     ! AC solver solves for area. Get the volume from area.
    !     where (maskarray_open_pipe .and. (solver == AC))
    !         volume = area * length
    !     endwhere

    !     ! print*, 'geometry debug for ', trim(subroutine_name)
    !     ! print*, 'Open pipe  ', maskarray_pipe_geometry
    !     ! print*, 'Area       ', area(997:1001)
    !     ! print*, 'AoverAfull ', AoverAfull
    !     ! print*, 'depth      ', depth
    !     ! print*, 'eta        ', eta

    !     maskarray_open_pipe  = nullvalueL
    !     maskarray_pipe_geometry = nullvalueL
    !     !% nullify temporary array
    !     nullify(maskarray_open_pipe, maskarray_pipe_geometry)
    !     next_eYN_temparray = next_eYN_temparray - 2

    !     if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    ! end subroutine open_pipe
    !
    !==========================================================================
    !==========================================================================
    !
    ! subroutine pipe_additional_geometric_properties &
    !     (elemI, elemR, elemYN, volume, length, zbottom, breadth, fulldepth,   &
    !     fullarea, zcrown, radius, area, eta, perimeter, depth, hyddepth,      &
    !     hydradius, topwidth, dHdA, elN, AoverAfull, YoverYfull, solver,       &
    !     geometry, ei_Temp, next_ei_temparray, ei_n_temp, eYN_Temp,            &
    !     next_eYN_temparray, eYN_n_temp, isfull, maskarray)
    !     !
    !     ! this subroutine gets the additional geometric properties of circular 
    !     ! pipe/junction-pipe geometry solved using the SVE/AC solver
    !     !
    !     character(64) :: subroutine_name = 'pipe_additional_geometric_properties'

    !     real,      target,     intent(inout)  :: elemR(:,:)
    !     logical,   target,     intent(inout)  :: elemYN(:,:)

    !     integer,   intent(in)       :: elemI(:,:)
    !     real,      intent(in)       :: volume(:), length(:), zbottom(:), breadth(:)
    !     real,      intent(in)       :: fulldepth(:), fullarea(:), zcrown(:), radius(:)
    !     integer,   intent(in)       :: ei_n_temp, eYN_n_temp
    !     integer,   intent(in)       :: solver(:), geometry(:), ei_Temp(:), eYN_Temp(:)
    !     logical,   intent(in)       :: maskarray(:), isfull(:)

    !     real,      intent(inout)    :: area(:), eta(:), perimeter(:), depth(:)
    !     real,      intent(inout)    :: hyddepth(:), hydradius(:), topwidth(:)
    !     real,      intent(inout)    :: dHdA(:), elN(:), AoverAfull(:), YoverYfull(:)
    !     integer,   intent(inout)    :: next_ei_temparray, next_eYN_temparray
    !     logical,   pointer          :: maskarray_pipe_geometry(:) 

    !     real :: af, bf, cf
    !     !--------------------------------------------------------------------------
    !     if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

    !     !% temporary mask space for geometry types that needs table interpolation
    !     maskarray_pipe_geometry  => elemYN(:,eYN_Temp(next_eYN_temparray))
    !     next_eYN_temparray = utility_advance_temp_array (next_eYN_temparray,eYN_n_temp)

    !     maskarray_pipe_geometry = nullvalueL

    !     !% set up the mask for circular geometry type as an input for table_lookup_mask subroutine
    !     maskarray_pipe_geometry = ( (maskarray) .and. (geometry == eCircular) )

    !     ! get normalized topwidth from lookup table for circular pipe
    !     call table_lookup_mask &
    !         (elemI, elemR, topwidth, YoverYfull, WCirc, NWCirc, maskarray_pipe_geometry, &
    !         ei_Temp, next_ei_temparray, ei_n_temp)
    !     ! get normalized hydradius from lookup table for circular pipe
    !     call table_lookup_mask &
    !         (elemI, elemR, hydradius, YoverYfull, RCirc, NRCirc, maskarray_pipe_geometry, &
    !         ei_Temp, next_ei_temparray, ei_n_temp)

    !     !% release the temporary maskarray_pipe_geometry for other special geometry types
    !     maskarray_pipe_geometry = nullvalueL

    !     !% HACK: Other special pipe geometry types are needed

    !     !%  additional geometric properties for OPEN PIPES
    !     where ( (maskarray) .and. (isfull .eqv. .false.) .and. (geometry == eCircular))
    !         ! get the depth by multiplying the normalized depth with fulldepth
    !         topwidth  = fulldepth * topwidth
    !         hydradius = onefourthR * fulldepth * hydradius
    !         perimeter = area / hydradius
    !     elsewhere ( (maskarray) .and. (isfull .eqv. .false.) .and. (geometry == eRectangular))
    !         topwidth  = breadth
    !         perimeter = breadth + twoR * depth 
    !         hydradius = area / perimeter
    !         hyddepth  = area / topwidth
    !     endwhere    

    !     !%  additional geometric properties for FULL PIPES
    !     where ((maskarray) .and. (isfull) .and. (geometry == eCircular))
    !         topwidth  = zeroR
    !         hydradius = onefourthR * fulldepth 
    !         perimeter = area / hydradius
    !     elsewhere ((maskarray) .and. (isfull) .and. (geometry == eRectangular))
    !         topwidth  = zeroR
    !         hydradius = (breadth * fulldepth) / (twoR * (breadth + fulldepth))
    !         perimeter = twoR * (breadth + fulldepth)
    !         hyddepth  = eta - zBottom
    !     endwhere 

    !     ! Get modified hydraulic depth for circular pipe from pipeAC Hodges 2020
    !     where ((maskarray) .and. (geometry == eCircular) .and. (AoverAfull .GT. onehalfR))
    !         hyddepth   = eta - zbottom + radius * (onefourthR * pi - oneR)
                             
    !     elsewhere ((maskarray) .and. (geometry == eCircular) .and. (AoverAfull .LT. onehalfR))
    !         hyddepth   = max(area / topwidth, zeroR)
    !     endwhere

    !     ! save the length scale values for pipes
    !     where(maskarray)
    !         elN = hyddepth
    !     endwhere

    !     !% Get dHdA (pipeAC2020)
    !     af = 1.29
    !     bf = 0.66
    !     cf = 0.34
    !     !% any kind of full pipes, dHdA = 0 
    !     where (maskarray .and. isfull)
    !         dHdA = zeroR
    !     endwhere
    !     !% RECTANGULAR PIPE
    !     where(maskarray .and. (geometry == eRectangular) .and. (isfull .eqv. .false.))
    !         dHdA = oneR / topwidth
    !     endwhere

    !     !% CIRCULAR PIPE
    !     where( (maskarray) .and. (geometry == eCircular) .and. &
    !            (AoverAfull .LT. onehalfR) ) 
    !         dHdA = (af * bf / ( (pi**bf) * (radius**(twoR*bf - oneR)) ) ) & 
    !                * area**(bf - oneR) + cf / (pi * radius)

    !     elsewhere( (maskarray) .and. (geometry == eCircular)   .and. &
    !                (AoverAfull .GT. onehalfR)                  .and. &
    !                (AoverAfull .LT. oneR) )
    !         dHdA = (af * bf / ( (pi**bf) * (radius**(twoR*bf - oneR))) ) &
    !                * (pi * (radius**twoR) - area)**(bf - oneR)  + cf / (pi * radius)
    !     endwhere

    !     ! print*, 'geometry debug for ', trim(subroutine_name)
    !     ! print*, 'is full    ', isfull
    !     ! print*, 'area        ', area(997:1001)
    !     ! print*, 'eta         ', eta(997:1001)
    !     ! print*,'----------------------------------------'
    !     ! print*, 'topwidth   ', topwidth
    !     ! print*, 'hydradius  ', hydradius
    !     ! print*, 'perimeter  ', perimeter
    !     ! print*, 'hyddepth   ', hyddepth(997:1001)
    !     ! print*, 'dHdA       ', dHdA

    !     maskarray_pipe_geometry = nullvalueL
    !     !% nullify temporary array
    !     nullify(maskarray_pipe_geometry)
    !     next_eYN_temparray = next_eYN_temparray - 1

    !     if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    ! end subroutine pipe_additional_geometric_properties
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine rectangular_junction_leg &
        (elemMR, elemMI, faceR, face_per_elemM, eMi_nfacesDir, eMi_MfaceDir, &
        eMr_AreaDir, eMr_ZbottomDir, eMr_BreadthScaleDir,  &
        eMr_TopwidthDir, eMr_EtaDir, eMr_HydDepthDir, eMr_EtaOld, fr_Eta_dir, method_EtaM)
        !
        ! computes geometry in the individual junction branches for a channel
        ! element with multiple connections
        !
        character(64) :: subroutine_name = 'rectangular_junction_leg'

        real,      target,     intent(in out)  :: elemMR(:,:)
        real,      target,     intent(in)      :: faceR(:,:)
        integer,   target,     intent(in)      :: elemMI(:,:)

        integer,               intent(in)      :: face_per_elemM, eMi_nfacesDir
        integer,               intent(in)      :: eMr_AreaDir(:)
        integer,               intent(in)      :: eMr_ZbottomDir(:)
        integer,               intent(in)      :: eMr_BreadthScaleDir(:)
        integer,               intent(in)      :: eMr_TopwidthDir(:)
        integer,               intent(in)      :: eMr_EtaDir(:)
        integer,               intent(in)      :: eMr_HydDepthDir(:)
        integer,               intent(in)      :: eMi_MfaceDir(:)
        integer,               intent(in)      :: eMr_EtaOld, fr_Eta_dir
        integer,               intent(in)      :: method_EtaM
        integer    :: mm

        integer,   pointer :: fdir(:)
        real,      pointer :: eta(:), etaM(:), etaMold(:), area(:), zbottom(:)
        real,      pointer :: depth(:), topwidth(:), breadth(:), etaFace(:)
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        etaM    => elemMR(:,eMr_Eta)
        etaMold => elemMR(:,eMr_EtaOld)

        !% cycle over branches. Dir is either {Up, Dn}
        do mm=1,face_per_elemM

            ! known
            zbottom     => elemMR(:,eMr_ZbottomDir(mm))
            breadth     => elemMR(:,eMr_BreadthScaleDir(mm))
            fdir        => elemMI(:,eMi_MfaceDir(mm))
            etaFace     => faceR (:,fr_Eta_dir)

            ! computed
            eta         => elemMR(:,eMr_EtaDir(mm))
            depth       => elemMR(:,eMr_HydDepthDir(mm))
            area        => elemMR(:,eMr_AreaDir(mm))
            topwidth    => elemMR(:,eMr_TopwidthDir(mm))

            select case (method_EtaM)
              case (0)
                !%  used for initiation when there are no face values
                !%  simple injection of junction eta
                where ( (elemMI(:,eMi_geometry)  == eRectangular)     .and. &
                    (elemMI(:,eMi_elem_type) == eJunctionChannel) .and. &
                    (elemMI(:,eMi_nfacesDir) >= mm)  )
                    eta = etaM
                endwhere
              case (1)
                !%  used before face values are updated (ie. etaF is old)
                !%  inject etaM with estimated correction for old gradient
                where ( (elemMI(:,eMi_geometry)  == eRectangular)     .and. &
                    (elemMI(:,eMi_elem_type) == eJunctionChannel) .and. &
                    (elemMI(:,eMi_nfacesDir) >= mm)  )
                    eta = etaM + onehalfR * (etaFace(fdir) - etaMold)
                endwhere
              case (2)
                !%  used after face values are updated (ie. etaF is new)
                !%  Average of face and junction value
                where ( (elemMI(:,eMi_geometry)  == eRectangular)     .and. &
                    (elemMI(:,eMi_elem_type) == eJunctionChannel) .and. &
                    (elemMI(:,eMi_nfacesDir) >= mm)  )
                    eta = onehalfR * (etaFace(fdir) + etaM)
                endwhere
              case default
                print *, method_EtaM
                print *, 'error: unexpected value of ',method_EtaM,' for method_EtaM in ',trim(subroutine_name)
            end select

            where ( (elemMI(:,eMi_geometry)  == eRectangular)     .and. &
                (elemMI(:,eMi_nfacesDir) >= mm)  )
                area = (eta - zbottom) * breadth
                topwidth = breadth
                depth = area / topwidth
            endwhere
        enddo

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine rectangular_junction_leg
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine width_depth_quadratic_function (elemI, elemR, ei_geometry, ei_elem_type, elem_type_value, &
        er_Length, er_Area, er_Volume, widthDepthData, w_d_angle, w_d_widthAtLayerTop, w_d_depthAtLayerTop, &
        w_d_perimeterBelowThisLayer, a_diff)
        ! This is the subroutine to update the irregular XS (width-depth system)
        ! Since the irregular XS is not linear, we cannot directly use interpolation to get the right geometry 
        ! between layers, solve the quadratic function here and get the depth DD
        
        character(64) :: subroutine_name = 'width_depth_quadratic_function'

        integer,                intent(in)      :: elemI(:,:)
        real,       target,     intent(in out)  :: elemR(:,:)
        real,       target,     intent(in)      :: widthDepthData(:,:,:)

        integer,                intent(in)      :: ei_geometry, ei_elem_type, elem_type_value
        integer,                intent(in)      :: er_Length
        integer,                intent(in)      :: er_Area
        integer,                intent(in)      :: er_Volume

        real,                   intent(in out)  :: w_d_angle(:)
        real,                   intent(in out)  :: w_d_widthAtLayerTop(:)
        real,                   intent(in out)  :: w_d_depthAtLayerTop(:)
        real,                   intent(in out)  :: w_d_perimeterBelowThisLayer(:)

        real,       pointer                     :: volume(:)
        real,       pointer                     :: area(:)
        real,       pointer                     :: length(:)

        real,       pointer                     :: areaTotalBelowThisLayer(:,:)
        real,       pointer                     :: widthAtLayerTop(:,:)
        real,       pointer                     :: depthAtLayerTop(:,:)
        real,       pointer                     :: angle(:,:)
        real,       pointer                     :: perimeterBelowThisLayer(:,:)
        
        
        real,       dimension(:), allocatable   :: area_difference

        real,                   intent(in out)  :: a_diff(:)
        integer                                 :: ind

        integer                                 :: ii, linkIDTemp
        real                                    :: local_diff

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        ! inputs
        volume        => elemR(:,er_Volume)
        length        => elemR(:,er_Length)
        area          => elemR(:,er_Area)

        areaTotalBelowThisLayer => widthDepthData (:,:, wd_areaTotalBelowThisLayer)
        widthAtLayerTop         => widthDepthData (:,:, wd_widthAtLayerTop)
        depthAtLayerTop         => widthDepthData (:,:, wd_depthAtLayerTop)
        angle                   => widthDepthData (:,:, wd_angle)
        perimeterBelowThisLayer => widthDepthData (:,:, wd_perimeterBelowThisLayer)



        ! area_difference is used for compute the area difference within each layer
        allocate (area_difference(size(widthDepthData,2)))


        !Calculate the area ahead
        where ( elemI(:, ei_geometry) == eWidthDepth .and. (elemI(:,ei_elem_type) == elem_type_value) )
            area        = volume / length
        endwhere
                
        !do loop to get the ind array and a_diff array, and use these two arrays to calculate AA, BB, CC, and DD for quadratic function 
        do ii=1, size(area,1)
            if ( (elemI(ii,ei_geometry)  == eWidthDepth) .and. (elemI(ii,ei_elem_type) == elem_type_value) ) then
                    linkIDTemp = elemI(ii,e2i_link_ID)

                    area_difference  = zeroR
                    area_difference(:) = area(ii) - areaTotalBelowThisLayer(linkIDTemp,:)
                    
                    local_diff = (areaTotalBelowThisLayer(linkIDTemp, minloc(abs(area_difference), DIM=1)) - area(ii) )

                     ! Find the closest area layer
                    if      ( local_diff > 0)   then !its between ind and ind-1
                        ind                             = minloc(abs(area_difference), DIM=1) -1 ! Since its between i & i-1, use the lower layer (i-1)
                        a_diff(ii)                      = (area(ii) - areaTotalBelowThisLayer(linkIDTemp,ind))
                        w_d_angle(ii)                   = angle(linkIDTemp, ind)
                        w_d_widthAtLayerTop(ii)         = widthAtLayerTop(linkIDTemp, ind)
                        w_d_depthAtLayerTop(ii)         = depthAtLayerTop(linkIDTemp, ind)
                        w_d_perimeterBelowThisLayer(ii) = perimeterBelowThisLayer(linkIDTemp, ind)
                    else if ( local_diff < 0)   then !its between ind and ind+1
                        ind                             = minloc(abs(area_difference), DIM=1) ! The floor is layer i
                        a_diff(ii)                      = (area(ii) - areaTotalBelowThisLayer(linkIDTemp, ind))
                        w_d_angle(ii)                   = angle(linkIDTemp, ind)
                        w_d_widthAtLayerTop(ii)         = widthAtLayerTop(linkIDTemp, ind)
                        w_d_depthAtLayerTop(ii)         = depthAtLayerTop(linkIDTemp, ind)
                        w_d_perimeterBelowThisLayer(ii) = perimeterBelowThisLayer(linkIDTemp, ind)
                    else if ( local_diff == 0)  then ! exactly at layer boundary
                        ind                             = minloc(abs(area_difference), DIM=1)
                        a_diff(ii)                      = 0
                        w_d_angle(ii)                   = angle(linkIDTemp, ind)
                        w_d_widthAtLayerTop(ii)         = widthAtLayerTop(linkIDTemp, ind)
                        w_d_depthAtLayerTop(ii)         = depthAtLayerTop(linkIDTemp, ind)
                        w_d_perimeterBelowThisLayer(ii) = perimeterBelowThisLayer(linkIDTemp, ind)
                    else
                        print *, "Error in Width-Depth Interpolation in element geometry"
                        print *, "local_diff = ", local_diff
                        print *, "areadifference = ", area_difference
                        print *, "=============="
                        print *, "area = ", area(ii)
                        print *, "=============="
                        print *, "areaTotalBelowThisLayer=", areaTotalBelowThisLayer(linkIDTemp,:)
                        stop
                
                    endif
            else
                a_diff(ii)  = nullvalueR
            endif
        enddo
            

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine width_depth_quadratic_function
    !
    !==========================================================================
    ! END OF MODULE element_geometry
    !==========================================================================
end module element_geometry