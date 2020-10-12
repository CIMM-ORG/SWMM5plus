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
        call channel_or_junction &
            (elem2R, elem2I, e2i_geometry, e2i_elem_type, eChannel, e2r_Length,   &
            e2r_Zbottom, e2r_BreadthScale, e2r_Topwidth, e2r_Area, e2r_Eta,      &
            e2r_Perimeter, e2r_Depth, e2r_HydDepth, e2r_HydRadius,               &
            e2r_FullDepth, e2r_Volume_new, e2r_LeftSlope, e2r_RightSlope,        &
            e2r_ParabolaValue, e2r_Temp, next_e2r_temparray, e2r_n_temp, ID,     &
            numberPairs, ManningsN, Length, zBottom, xDistance, Breadth,         &
            widthDepthData, cellType)

        !%   geomety for junctions
        call channel_or_junction &
            (elemMR, elemMI, eMi_geometry, eMi_elem_type, eJunctionChannel,       &
            eMr_Length, eMr_Zbottom, eMr_BreadthScale, eMr_Topwidth, eMr_Area,   &
            eMr_Eta, eMr_Perimeter, eMr_Depth, eMr_HydDepth, eMr_HydRadius,      &
            eMr_FullDepth, eMr_Volume_new, eMr_LeftSlope, eMr_RightSlope,        &
            eMr_ParabolaValue, eMr_Temp, next_eMr_temparray, eMr_n_temp, ID,     &
            numberPairs, ManningsN, Length, zBottom, xDistance, Breadth,         &
            widthDepthData, cellType)

        !%   geometry for pipe
        !%   pipe geometry is solved using area and eta. Thus, area and eta are the input
        !%   pipe volume is saved into the volume_new column for the consistancy of the code 
        !%   volume is updated finally after the rk2/ac-rk2 loop
        call pipe_or_junction &
            (elem2R, elem2I, elem2YN, e2i_geometry, e2i_elem_type, ePipe, e2i_solver,   &
            e2r_Length, e2r_Zbottom, e2r_BreadthScale, e2r_Topwidth, e2r_Area,          &
            e2r_Eta, e2r_Perimeter, e2r_Depth, e2r_HydDepth, e2r_HydRadius, e2r_Radius, &
            e2r_Volume_new, e2r_FullDepth, e2r_FullArea, e2r_Zcrown, e2r_dHdA, e2r_elN, &
            e2YN_IsSurcharged, e2r_Temp, next_e2r_temparray, e2r_n_temp, e2i_Temp,      &
            next_e2i_temparray, e2i_n_temp, e2YN_Temp, next_e2YN_temparray, e2YN_n_temp)


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


        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine geometry_update
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine channel_or_junction &
        (elemR, elemI, ei_geometry, ei_elem_type, elem_type_value, er_Length, &
        er_Zbottom, er_BreadthScale, er_Topwidth, er_Area, er_Eta,           &
        er_Perimeter, er_Depth, er_HydDepth, er_HydRadius, er_FullDepth,     &
        er_Volume, er_LeftSlope, er_RightSlope, er_ParabolaValue,            &
        er_Temp, next_er_temparray, er_n_temp, wdID, wdnumberPairs,          &
        wdManningsN, wdLength, wdzBottom, wdxDistance, wdBreadth,            &
        widthDepthData, wdcellType)
        !
        ! computes element geometry for a rectangular channel or a channeljunction
        !
        character(64) :: subroutine_name = 'channel_or_junction'

        real,      target,     intent(in out)  :: elemR(:,:)

        integer,   intent(in)      :: elemI(:,:)
        integer,   intent(in)      :: ei_geometry, ei_elem_type, elem_type_value
        integer,   intent(in)      :: er_Length, er_Zbottom, er_BreadthScale
        integer,   intent(in)      :: er_Area, er_Eta, er_Perimeter, er_Topwidth
        integer,   intent(in)      :: er_Depth, er_HydDepth, er_HydRadius
        integer,   intent(in)      :: er_Volume, er_FullDepth, er_n_temp
        integer,   intent(in)      :: er_LeftSlope, er_RightSlope, er_ParabolaValue
        integer,   intent(in)      :: er_Temp(:)
        integer,   intent(inout)   :: next_er_temparray

        integer, target, intent(in out)    :: wdID(:)
        integer, target, intent(in out)    :: wdnumberPairs(:)
        real,    target, intent(in out)    :: wdManningsN(:)
        real,    target, intent(in out)    :: wdLength(:)
        real,    target, intent(in out)    :: wdzBottom(:)
        real,    target, intent(in out)    :: wdxDistance(:)
        real,    target, intent(in out)    :: wdBreadth(:)
        real,    target, intent(in out)    :: widthDepthData(:,:,:)
        type(string), target, intent(in out)   :: wdcellType(:)


        real, pointer :: volume(:), length(:), zbottom(:), breadth(:)
        real, pointer :: area(:), eta(:), perimeter(:), depth(:), hyddepth(:)
        real, pointer :: hydradius(:), topwidth(:), fulldepth(:)
        real, pointer :: leftSlope(:), rightSlope(:), parabolaValue(:)
        real, pointer :: AoverAfull(:), YoverYfull(:)

        real, pointer :: widthAtLayerTop(:,:), depthAtLayerTop(:,:), areaThisLayer(:,:)
        real, pointer :: areaTotalBelowThisLayer(:,:), dWidth(:,:)
        real, pointer :: dDepth(:,:), angle(:,:), perimeterBelowThisLayer(:,:)
        real, dimension(:), allocatable :: area_difference, local_difference

        real :: AA, BB, CC, DD
        integer :: ii,ind, linkIDTemp

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        ! inputs
        volume        => elemR(:,er_Volume)
        length        => elemR(:,er_Length)
        zbottom       => elemR(:,er_Zbottom)
        breadth       => elemR(:,er_BreadthScale)
        fulldepth     => elemR(:,er_FullDepth)
        leftSlope     => elemR(:,er_LeftSlope)
        rightSlope    => elemR(:,er_RightSlope)
        parabolaValue => elemR(:,er_ParabolaValue)

        ! outputs
        area       => elemR(:,er_Area)
        eta        => elemR(:,er_Eta)
        perimeter  => elemR(:,er_Perimeter)
        depth      => elemR(:,er_Depth)
        hyddepth   => elemR(:,er_HydDepth)
        hydradius  => elemR(:,er_HydRadius)
        topwidth   => elemR(:,er_Topwidth)

        widthAtLayerTop         => widthDepthData (:,:, wd_widthAtLayerTop)
        depthAtLayerTop         => widthDepthData (:,:, wd_depthAtLayerTop)
        areaThisLayer           => widthDepthData (:,:, wd_areaThisLayer)
        areaTotalBelowThisLayer => widthDepthData (:,:, wd_areaTotalBelowThisLayer)
        dWidth                  => widthDepthData (:,:, wd_Dwidth)
        dDepth                  => widthDepthData (:,:, wd_Ddepth)
        angle                   => widthDepthData (:,:, wd_angle)
        perimeterBelowThisLayer => widthDepthData (:,:, wd_perimeterBelowThisLayer)

        allocate (area_difference(size(widthDepthData,2)))
        allocate (local_difference(size(widthDepthData,2)))

        !% temporary array for circular geometry update
        AoverAfull => elemR(:,er_Temp(next_er_temparray))
        next_er_temparray = utility_advance_temp_array (next_er_temparray,er_n_temp)

        YoverYfull => elemR(:,er_Temp(next_er_temparray))
        next_er_temparray = utility_advance_temp_array (next_er_temparray,er_n_temp)

        AoverAfull = nullvalueR
        YoverYfull = nullvalueR

        where ( (elemI(:,ei_geometry)  == eRectangular) .and. &
            (elemI(:,ei_elem_type) == elem_type_value    )         )
            area        = volume / length
            eta         = zbottom + (area / breadth)
            topwidth    = breadth
            perimeter   = breadth + 2.0 * ( eta - zbottom )
            hyddepth    = area / breadth
            depth       = hyddepth
            hydradius   = area / perimeter

        elsewhere ( (elemI(:,ei_geometry)  == eParabolic) .and. &
            (elemI(:,ei_elem_type) == elem_type_value    )         )
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

        elsewhere ( (elemI(:,ei_geometry)  == eTrapezoidal) .and. &
            (elemI(:,ei_elem_type) == elem_type_value    )         )
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

        elsewhere ( (elemI(:,ei_geometry)  == eTriangular) .and. &
            (elemI(:,ei_elem_type) == elem_type_value    )         )
            area        = volume / length
            depth       = sqrt(abs(area/(onehalfR*(leftSlope + rightSlope))))
            hyddepth    = onehalfR * depth
            eta         = zbottom + hyddepth
            topwidth    = (leftSlope + rightSlope) * depth
            perimeter   = depth * (sqrt(oneR + leftSlope**twoR) + sqrt(oneR + rightSlope**twoR))
            hydradius   = area / perimeter

        endwhere

        do ii=1, size(volume,1)
            if ( (elemI(ii,ei_geometry)  == eWidthDepth) .and. &
                (elemI(ii,ei_elem_type) == elem_type_value    )         ) then

                linkIDTemp = elemI(ii,e2i_link_ID)

                area_difference  = zeroR
                local_difference = zeroR
                area (ii) = volume(ii) / length(ii)
                area_difference(:) = area (ii) - areaTotalBelowThisLayer(linkIDTemp,:)
                local_difference(:) = area_difference(:) - areaThisLayer(linkIDTemp,:)
                ind = findloc(sign(oneR, area_difference(:)*local_difference(:)), -1.0, DIM=1)

                if (ind == 0) then
                    ind = size(widthAtLayerTop(linkIDTemp,:),1)
                endif

                AA = oneR/tan(angle(linkIDTemp,ind))
                BB = widthAtLayerTop(linkIDTemp,ind) - dWidth(linkIDTemp,ind)
                CC = - area_difference(ind)
                DD = (-BB + sqrt(BB**twoR - fourR*AA*CC))/(twoR*AA)

                hyddepth (ii)  = DD + depthAtLayerTop(linkIDTemp,ind) - dDepth(linkIDTemp,ind)
                eta (ii)       = zbottom (ii) + hyddepth (ii)
                depth (ii)     = hyddepth(ii)
                topwidth (ii)  = widthAtLayerTop(linkIDTemp,ind) - (dDepth(linkIDTemp,ind)-DD) &
                    *dWidth(linkIDTemp,ind)/dDepth(linkIDTemp,ind)
                perimeter (ii) = perimeterBelowThisLayer(linkIDTemp,ind) + twoR * DD/sin(angle(linkIDTemp,ind))
                hydradius (ii) = area(ii) / perimeter(ii)
            endif

            ! if ( (elemI(ii,ei_geometry)  == eCircular) .and. &
            !     (elemI(ii,ei_elem_type) == elem_type_value    )         ) then

            !     area(ii)        = volume(ii) / length(ii)
            !     AoverAfull(ii)  = area(ii) / (onefourthR * pi * fulldepth(ii) ** twoR)

            !     if (AoverAfull(ii) < 0.04) then
            !         depth(ii) = fulldepth(ii) * get_theta_of_alpha(AoverAfull(ii))
            !     else
            !         depth(ii) = fulldepth(ii) * table_lookup(AoverAfull(ii), YCirc, NYCirc)
            !     endif
            !     YoverYfull(ii) = depth(ii) / fulldepth(ii)
            !     topwidth(ii)   = fulldepth(ii) * table_lookup(YoverYfull(ii), WCirc, NWCirc)
            !     hyddepth(ii)   = area(ii) / topwidth(ii)
            !     eta(ii)        = zbottom(ii) + hyddepth(ii)
            !     hydradius(ii)  = onefourthR * fulldepth (ii) * table_lookup(YoverYfull(ii), RCirc, NRCirc)
            !     perimeter(ii)  = area(ii) / hydradius(ii)

            ! endif
        enddo

        AoverAfull = nullvalueR
        YoverYfull = nullvalueR

        !% nullify temporary array
        nullify(AoverAfull, YoverYfull)

        next_er_temparray  = next_er_temparray  - 2

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine channel_or_junction
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine pipe_or_junction &
        (elemR, elemI, elemYN, ei_geometry, ei_elem_type, elem_type_value,   &
        ei_solver, er_Length, er_Zbottom, er_BreadthScale, er_Topwidth,      &
        er_Area, er_Eta, er_Perimeter, er_Depth, er_HydDepth, er_HydRadius,  &
        er_Radius, er_Volume, er_FullDepth, er_FullArea, er_Zcrown, er_dHdA, &
        er_elN, eYN_IsSurcharged, er_Temp, next_er_temparray, er_n_temp,     &
        ei_Temp, next_ei_temparray, ei_n_temp, eYN_Temp, next_eYN_temparray, &
        eYN_n_temp)
        !
        !% computes geometry for pipe and junction pipe
        !
        character(64) :: subroutine_name = 'pipe_or_junction'

        real,      target,     intent(inout)  :: elemR(:,:)
        logical,   target,     intent(inout)  :: elemYN(:,:)
        integer,   target,     intent(in)     :: elemI(:,:)

        integer,   intent(in)      :: ei_geometry, ei_elem_type, elem_type_value, ei_solver
        integer,   intent(in)      :: er_Length, er_Zbottom, er_BreadthScale, er_Perimeter 
        integer,   intent(in)      :: er_Topwidth, er_Area, er_Eta, er_dHdA, er_elN
        integer,   intent(in)      :: er_Depth, er_HydDepth, er_HydRadius, er_Radius
        integer,   intent(in)      :: er_Volume, er_FullDepth, er_FullArea, er_Zcrown
        integer,   intent(in)      :: er_n_temp, ei_n_temp, eYN_n_temp
        integer,   intent(in)      :: eYN_IsSurcharged 
        integer,   intent(in)      :: er_Temp(:), ei_Temp(:), eYN_Temp(:)

        integer,   intent(inout)   :: next_er_temparray, next_ei_temparray, next_eYN_temparray

        real,    pointer    :: volume(:), length(:), zbottom(:), breadth(:), hyddepth(:)
        real,    pointer    :: area(:),  eta(:), perimeter(:), depth(:), dHdA(:), elN(:)
        real,    pointer    :: hydradius(:), topwidth(:), fulldepth(:), zcrown(:)
        real,    pointer    :: radius(:), AoverAfull(:), YoverYfull(:), fullarea(:)
        integer, pointer    :: solver(:)
        logical, pointer    :: isfull(:), maskarray(:)

        integer :: ii

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% inputs
        
        length     => elemR(:,er_Length)
        zbottom    => elemR(:,er_Zbottom)
        breadth    => elemR(:,er_BreadthScale)
        fulldepth  => elemR(:,er_FullDepth)
        zcrown     => elemR(:,er_Zcrown)
        fullarea   => elemR(:,er_FullArea)
        radius     => elemR(:,er_Radius)
        solver     => elemI(:,ei_solver)
        isfull     => elemYN(:,eYN_IsSurcharged)

        !% outputs
        area       => elemR(:,er_Area)
        eta        => elemR(:,er_Eta)
        perimeter  => elemR(:,er_Perimeter)
        depth      => elemR(:,er_Depth)
        hyddepth   => elemR(:,er_HydDepth)
        hydradius  => elemR(:,er_HydRadius)
        topwidth   => elemR(:,er_Topwidth)
        dHdA       => elemR(:,er_dHdA)
        elN        => elemR(:,er_elN)
        volume     => elemR(:,er_Volume)

        !% temporary array for geometry update
        AoverAfull => elemR(:,er_Temp(next_er_temparray))
        next_er_temparray = utility_advance_temp_array (next_er_temparray,er_n_temp)

        YoverYfull => elemR(:,er_Temp(next_er_temparray))
        next_er_temparray = utility_advance_temp_array (next_er_temparray,er_n_temp)

        maskarray  => elemYN(:,eYN_Temp(next_eYN_temparray))
        next_eYN_temparray = utility_advance_temp_array (next_eYN_temparray,eYN_n_temp)

        AoverAfull = nullvalueR
        YoverYfull = nullvalueR
        maskarray  = nullvalueL

        !% mask for circular pipes/junction pipe
        maskarray = ( (elemI(:,ei_geometry)  == eCircular)       .and. &
                      (elemI(:,ei_elem_type) == elem_type_value) )

        !% For open pipe, set eta to zero for later error detection
        where (maskarray .and. (isfull .eqv. .false.)) 
              eta = zeroR
        endwhere

        !% Open circular pipes that become full
        call open_circular_pipe_transition_to_full &
            (elemI, elemR, elemYN, volume, length, zbottom, breadth, fulldepth,   &
            fullarea, zcrown, radius, area, eta, perimeter, depth, hyddepth,      &
            hydradius, topwidth, AoverAfull, YoverYfull, solver, ei_Temp,         &
            next_ei_temparray, ei_n_temp, eYN_Temp, next_eYN_temparray,           &
            eYN_n_temp, isfull, maskarray)

        !% Set the full pipe area
        !% These cells already have eta directly updated from the time-stepping.
        where (maskarray .and. (isfull .eqv. .true.)) 
            area   = fullarea
            volume = area * length
            depth  = fulldepth
            AoverAfull = oneR  
        endwhere

        !% NOTE: at this point, isfull will consist of those pipes that
        !% were designated as "full" coming into this routine and those that
        !% were detected as becoming full due to the increase in area. However
        !% the array ALSO is true (incorrectly) for those pipes that dropped
        !% from full to open during the prior time step. These must be handled
        !% separately.

        !% Full pipes that become open
        call full_circular_pipe_transition_to_open &
            (elemI, elemR, elemYN, volume, length, zbottom, breadth, fulldepth,  &
            fullarea, zcrown, radius, area, eta, perimeter, depth, hyddepth,     &
            hydradius, topwidth, AoverAfull, YoverYfull, solver, ei_Temp,        &
            next_ei_temparray, ei_n_temp, eYN_Temp, next_eYN_temparray,          &
            eYN_n_temp, isfull, maskarray)

        !% Open pipes
        !% These have areas/volumes and need eta computed
        call open_circular_pipe &
            (elemI, elemR, elemYN, volume, length, zbottom, breadth, fulldepth,  &
            fullarea, zcrown, radius, area, eta, perimeter, depth, hyddepth,     &
            hydradius, topwidth, AoverAfull, YoverYfull, solver, ei_Temp,        &
            next_ei_temparray, ei_n_temp, eYN_Temp, next_eYN_temparray,          &
            eYN_n_temp, isfull, maskarray)

        !% isfull reset
        !% set formerly full pipes that have become open to open
        where (maskarray .and. (isfull .eqv. .true.) .and. (eta < zcrown))
            isfull = .false.
        endwhere

        ! Now eta and area have been updated for all cases
        call circular_pipe_additional_geometric_properties &
            (elemI, elemR, volume, length, zbottom, breadth, fulldepth, fullarea,   &
            zcrown, radius, area, eta, perimeter, depth, hyddepth, hydradius,       &
            topwidth, dHdA, elN, AoverAfull, YoverYfull, solver, ei_Temp,           &
            next_ei_temparray, ei_n_temp, eYN_Temp, next_eYN_temparray, eYN_n_temp, &
            isfull, maskarray) 
              
        !% release temporary arrays
        AoverAfull = nullvalueR
        YoverYfull = nullvalueR
        maskarray  = nullvalueL

        !% nullify temporary array
        nullify(AoverAfull, YoverYfull, maskarray)
        next_er_temparray  = next_er_temparray  - 2
        next_eYN_temparray = next_eYN_temparray - 1

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine pipe_or_junction
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine open_circular_pipe_transition_to_full &
        (elemI, elemR, elemYN, volume, length, zbottom, breadth, fulldepth,   &
        fullarea, zcrown, radius, area, eta, perimeter, depth, hyddepth,      &
        hydradius, topwidth, AoverAfull, YoverYfull, solver, ei_Temp,         &
        next_ei_temparray, ei_n_temp, eYN_Temp, next_eYN_temparray,           &
        eYN_n_temp, isfull, maskarray)
        !
        ! OPEN PIPES THAT BECOME FULL ============
        ! Detect transition from open to full pipe
        ! area is advanced for open pipe. To eta is needed to be calculated here
        !
        character(64) :: subroutine_name = 'open_circular_pipe_transition_to_full'

        real,      target,     intent(inout)  :: elemR(:,:)
        logical,   target,     intent(inout)  :: elemYN(:,:)

        integer,   intent(in)       :: elemI(:,:)
        real,      intent(in)       :: length(:), zbottom(:), breadth(:)
        real,      intent(in)       :: fulldepth(:), fullarea(:), zcrown(:), radius(:)
        integer,   intent(in)       :: ei_n_temp, eYN_n_temp
        integer,   intent(in)       :: solver(:), ei_Temp(:), eYN_Temp(:)
        logical,   intent(in)       :: maskarray(:)

        real,      intent(inout)    :: volume(:), area(:), eta(:), perimeter(:), depth(:)
        real,      intent(inout)    :: hyddepth(:), hydradius(:), topwidth(:)
        real,      intent(inout)    :: AoverAfull(:), YoverYfull(:)
        integer,   intent(inout)    :: next_ei_temparray, next_eYN_temparray 
        logical,   intent(inout)    :: isfull(:)

        logical,    pointer         :: maskarray_pipe_transition(:)

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% temporary mask space for finding transiotional pipes
        maskarray_pipe_transition  => elemYN(:,eYN_Temp(next_eYN_temparray))
        next_eYN_temparray = utility_advance_temp_array (next_eYN_temparray,eYN_n_temp)

        maskarray_pipe_transition  = nullvalueL

        maskarray_pipe_transition = (maskarray .and. (isfull .eqv. .false.) &
                                               .and. (area .GE. fullarea) )
        where (maskarray_pipe_transition) 
            ! Set head above the pipe crown based on the excess area in
            ! from the time advance divided by pipe width
            eta  = zcrown + (area - fullarea)/(twoR * radius)
            area = fullarea
            AoverAfull = oneR
            volume = area * length
            isfull = .true.
        endwhere

        maskarray_pipe_transition  = nullvalueL
        !% nullify temporary array
        nullify(maskarray_pipe_transition)
        next_eYN_temparray = next_eYN_temparray - 1

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine open_circular_pipe_transition_to_full
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine full_circular_pipe_transition_to_open &
        (elemI, elemR, elemYN, volume, length, zbottom, breadth, fulldepth,   &
        fullarea, zcrown, radius, area, eta, perimeter, depth, hyddepth,      &
        hydradius, topwidth, AoverAfull, YoverYfull, solver, ei_Temp,         &
        next_ei_temparray, ei_n_temp, eYN_Temp, next_eYN_temparray,           &
        eYN_n_temp, isfull, maskarray)
        !
        ! FULL PIPES THAT TRANSITION TO OPEN ==================
        ! These have eta and need area computed
        ! Note that these are not re-designated as open until after all
        ! the eta and are computations are complete
        ! Detect full pipe that have become open
        !
        character(64) :: subroutine_name = 'full_circular_pipe_transition_to_open'

        real,      target,     intent(inout)  :: elemR(:,:)
        logical,   target,     intent(inout)  :: elemYN(:,:)

        integer,   intent(in)       :: elemI(:,:)
        real,      intent(in)       :: length(:), zbottom(:), breadth(:)
        real,      intent(in)       :: fulldepth(:), fullarea(:), zcrown(:), radius(:)
        integer,   intent(in)       :: ei_n_temp, eYN_n_temp
        integer,   intent(in)       :: solver(:), ei_Temp(:), eYN_Temp(:)
        logical,   intent(in)       :: maskarray(:), isfull(:)

        real,      intent(inout)    :: volume(:), area(:), eta(:), perimeter(:), depth(:)
        real,      intent(inout)    :: hyddepth(:), hydradius(:), topwidth(:)
        real,      intent(inout)    :: AoverAfull(:), YoverYfull(:)
        integer,   intent(inout)    :: next_ei_temparray, next_eYN_temparray 

        logical,    pointer         :: maskarray_pipe_transition(:)

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% temporary mask space for finding transiotional pipes
        maskarray_pipe_transition  => elemYN(:,eYN_Temp(next_eYN_temparray))
        next_eYN_temparray = utility_advance_temp_array (next_eYN_temparray,eYN_n_temp)

        maskarray_pipe_transition  = nullvalueL

        maskarray_pipe_transition = (maskarray .and. (isfull .eqv. .true.) &
                                               .and. (eta < zcrown) )
        where (maskarray_pipe_transition) 
            depth = eta - zbottom
            YoverYfull = depth / fulldepth
        endwhere

        ! get the normalized area from lookup table from full to open transitional pipe
        call table_lookup_mask &
            (elemI, elemR, AoverAfull, YoverYfull, ACirc, NACirc, maskarray_pipe_transition, &
            ei_Temp, next_ei_temparray, ei_n_temp)

        ! get the pipe area by multiplying the normalized area with full area   
        where (maskarray_pipe_transition)
            area   = AoverAfull * fullarea
            volume = area * length
        endwhere

        maskarray_pipe_transition  = nullvalueL
        !% nullify temporary array
        nullify(maskarray_pipe_transition)
        next_eYN_temparray = next_eYN_temparray - 1

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine full_circular_pipe_transition_to_open
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine open_circular_pipe &
        (elemI, elemR, elemYN, volume, length, zbottom, breadth, fulldepth,   &
        fullarea, zcrown, radius, area, eta, perimeter, depth, hyddepth,      &
        hydradius, topwidth, AoverAfull, YoverYfull, solver, ei_Temp,         &
        next_ei_temparray, ei_n_temp, eYN_Temp, next_eYN_temparray,           &
        eYN_n_temp, isfull, maskarray)
        !
        ! this subroutine gets the value of circular pipe/junction-pipe geometry
        ! solved using the SVE/AC solver
        !
        character(64) :: subroutine_name = 'open_circular_pipe'

        real,      target,     intent(inout)  :: elemR(:,:)
        logical,   target,     intent(inout)  :: elemYN(:,:)

        integer,   intent(in)       :: elemI(:,:)
        real,      intent(in)       :: length(:), zbottom(:), breadth(:)
        real,      intent(in)       :: fulldepth(:), fullarea(:), zcrown(:), radius(:)
        integer,   intent(in)       :: ei_n_temp, eYN_n_temp
        integer,   intent(in)       :: solver(:), ei_Temp(:), eYN_Temp(:)
        logical,   intent(in)       :: maskarray(:), isfull(:)

        real,      intent(inout)    :: volume(:), area(:), eta(:), perimeter(:), depth(:)
        real,      intent(inout)    :: hyddepth(:), hydradius(:), topwidth(:)
        real,      intent(inout)    :: AoverAfull(:), YoverYfull(:)
        integer,   intent(inout)    :: next_ei_temparray, next_eYN_temparray 

        logical,    pointer         :: maskarray_open_pipe(:)

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% temporary mask space for finding transiotional pipes
        maskarray_open_pipe  => elemYN(:,eYN_Temp(next_eYN_temparray))
        next_eYN_temparray = utility_advance_temp_array (next_eYN_temparray,eYN_n_temp)

        maskarray_open_pipe  = nullvalueL

        maskarray_open_pipe = (maskarray .and. (isfull .eqv. .false.))

        ! SVE solver solves for volume. Get the area from volume
        where (maskarray_open_pipe .and. (solver == SVE))
            area = volume / length
        endwhere

        ! get the normalized area to solve for other geometries
        where (maskarray_open_pipe)
            AoverAfull = area / fullarea
        endwhere

        ! get normalized depth from the lookup table
        call table_lookup_mask &
            (elemI, elemR, YoverYfull, AoverAfull, YCirc, NYCirc, maskarray_open_pipe, &
            ei_Temp, next_ei_temparray, ei_n_temp)

        where (maskarray_open_pipe)
            ! get the depth by multiplying the normalized depth with fulldepth
            depth = fulldepth * YoverYfull
            eta   = zbottom + depth
        endwhere 

        ! AC solver solves for area. Get the volume from area.
        where (maskarray_open_pipe .and. (solver == AC))
            volume = area * length
        endwhere

        maskarray_open_pipe  = nullvalueL
        !% nullify temporary array
        nullify(maskarray_open_pipe)
        next_eYN_temparray = next_eYN_temparray - 1

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine open_circular_pipe
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine circular_pipe_additional_geometric_properties &
        (elemI, elemR, volume, length, zbottom, breadth, fulldepth, fullarea, &
        zcrown, radius, area, eta, perimeter, depth, hyddepth, hydradius,     &
        topwidth, dHdA, elN, AoverAfull, YoverYfull, solver, ei_Temp,         &
        next_ei_temparray, ei_n_temp, eYN_Temp, next_eYN_temparray,           &
        eYN_n_temp, isfull, maskarray)
        !
        ! this subroutine gets the additional geometric properties of circular 
        ! pipe/junction-pipe geometry solved using the SVE/AC solver
        !
        character(64) :: subroutine_name = 'circular_pipe_additional_geometric_properties'

        real,      target,     intent(inout)  :: elemR(:,:)

        integer,   intent(in)       :: elemI(:,:)
        real,      intent(in)       :: volume(:), length(:), zbottom(:), breadth(:)
        real,      intent(in)       :: fulldepth(:), fullarea(:), zcrown(:), radius(:)
        integer,   intent(in)       :: ei_n_temp, eYN_n_temp
        integer,   intent(in)       :: solver(:), ei_Temp(:), eYN_Temp(:)
        logical,   intent(in)       :: maskarray(:), isfull(:)

        real,      intent(inout)    :: area(:), eta(:), perimeter(:), depth(:)
        real,      intent(inout)    :: hyddepth(:), hydradius(:), topwidth(:)
        real,      intent(inout)    :: dHdA(:), elN(:), AoverAfull(:), YoverYfull(:)
        integer,   intent(inout)    :: next_ei_temparray, next_eYN_temparray 

        real :: af, bf, cf
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        ! get normalized topwidth from lookup table
        call table_lookup_mask &
            (elemI, elemR, topwidth, YoverYfull, WCirc, NWCirc, maskarray, &
            ei_Temp, next_ei_temparray, ei_n_temp)
        ! get normalized hydradius from lookup table
        call table_lookup_mask &
            (elemI, elemR, hydradius, YoverYfull, RCirc, NRCirc, maskarray, &
            ei_Temp, next_ei_temparray, ei_n_temp)

        where (maskarray)
            ! get the depth by multiplying the normalized depth with fulldepth
            topwidth  = fulldepth * topwidth
            hydradius = onefourthR * fulldepth * hydradius
            perimeter = area / hydradius
        endwhere 

        ! Get modified hydraulic depth from pipeAC Hodges 2020
        where ((maskarray) .and. (AoverAfull .GT. onehalfR))
            hyddepth   = eta - zbottom + radius * (onefourthR * pi - oneR)
                             
        elsewhere ((maskarray) .and. (AoverAfull .LT. onehalfR))
            hyddepth   = max(area / topwidth, zeroR)
        endwhere

        ! save the length scale values for circular pipe
        where(maskarray)
            elN = hyddepth
        endwhere

        ! Get dHdA (pipeAC2020)
        af = 1.29
        bf = 0.66
        cf = 0.34
        where( (maskarray) .and. (AoverAfull .LT. onehalfR) ) 
            dHdA = (af * bf / ( (pi**bf) * (radius**(twoR*bf - oneR)) ) ) & 
                   * area**(bf - oneR) + cf / (pi * radius)
        elsewhere( (maskarray) .and. (AoverAfull .GT. onehalfR) .and. &
                  (AoverAfull .LT. oneR) )
            dHdA = (af * bf / ( (pi**bf) * (radius**(twoR*bf - oneR))) ) &
                   * (pi * (radius**twoR) - area)**(bf - oneR)  + cf / (pi * radius)
        elsewhere( (maskarray) .and. (AoverAfull .GE. oneR) ) 
            dHdA = zeroR
        endwhere

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine circular_pipe_additional_geometric_properties
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
    ! END OF MODULE element_geometry
    !==========================================================================
end module element_geometry
