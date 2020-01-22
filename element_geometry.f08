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
     elemMR, elemMI, elemMYN, eMr_VolumeColumn, &
     faceR, faceI, bcdataDn, bcdataUp, thisTime, method_EtaM, &
     ID, numberPairs, ManningsN, Length, zBottom, xDistance, &
     Breadth, widthDepthData, cellType)
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
    (elem2R, elem2I, e2r_VolumeColumn, &
     elemMR, elemMI, eMr_VolumeColumn, faceR, eMr_EtaOld, method_EtaM, &
     ID, numberPairs, ManningsN, Length, zBottom, xDistance, &
     Breadth, widthDepthData, cellType)

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
    (elem2R, elem2I, e2r_Volume_new, &
     elemMR, elemMI, eMr_Volume_new, faceR, eMr_EtaOld, method_EtaM, &
     ID, numberPairs, ManningsN, Length, zBottom, xDistance, &
     Breadth, widthDepthData, cellType)
!
! Note that volume used is in a eTr storage location so that the update
! can be used on a temporary volume
!
 character(64) :: subroutine_name = 'geometry_update'

 real,      intent(in out)  :: elem2R(:,:), elemMR(:,:)
 real,      intent(in)      :: faceR(:,:)
 integer,   intent(in)      :: elem2I(:,:), elemMI(:,:)
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

!%  rectangular geometry for channels
 call channel_or_junction &
    (elem2R, elem2I, &
     e2i_geometry, e2i_elem_type, eChannel,  &
     e2r_Length, e2r_Zbottom, e2r_BreadthScale, e2r_Topwidth, e2r_Area, e2r_Eta, &
     e2r_Perimeter, e2r_Depth, e2r_HydDepth, e2r_HydRadius, e2r_Volume_new,    &
     e2r_LeftSlope, e2r_RightSlope, e2r_ParabolaValue, &
     ID, numberPairs, ManningsN, Length, zBottom, xDistance, &
     Breadth, widthDepthData, cellType)

!%  rectangular geomety for junctions
 call channel_or_junction &
    (elemMR, elemMI, &
     eMi_geometry, eMi_elem_type, eJunctionChannel,  &
     eMr_Length, eMr_Zbottom, eMr_BreadthScale, eMr_Topwidth, eMr_Area, eMr_Eta, &
     eMr_Perimeter, eMr_Depth, eMr_HydDepth, eMr_HydRadius, eMr_Volume_new,    &
     eMr_LeftSlope, eMr_RightSlope, eMr_ParabolaValue, &
     ID, numberPairs, ManningsN, Length, zBottom, xDistance, &
     Breadth, widthDepthData, cellType)

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
    (elemR, elemI, &
     ei_geometry, ei_elem_type, elem_type_value,  &
     er_Length, er_Zbottom, er_BreadthScale, er_Topwidth, er_Area, er_Eta, &
     er_Perimeter, er_Depth, er_HydDepth, er_HydRadius, er_Volume,    &
     er_LeftSlope, er_RightSlope, er_ParabolaValue, &
     wdID, wdnumberPairs, wdManningsN, wdLength, wdzBottom, wdxDistance, &
     wdBreadth, widthDepthData, wdcellType)
!
! computes element geometry for a rectangular channel or a channeljunction
!
 character(64) :: subroutine_name = 'channel_or_junction'

 real,      target,     intent(in out)  :: elemR(:,:)

 integer,   intent(in)      :: elemI(:,:)
 integer,   intent(in)      :: ei_geometry, ei_elem_type, elem_type_value
 integer,   intent(in)      :: er_Length, er_Zbottom, er_BreadthScale
 integer,   intent(in)      :: er_Area, er_Eta, er_Perimeter, er_Topwidth
 integer,   intent(in)      :: er_Depth, er_HydDepth, er_HydRadius, er_Volume
 integer,   intent(in)      :: er_LeftSlope, er_RightSlope, er_ParabolaValue
 
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
 real, pointer :: hydradius(:), topwidth(:)
 real, pointer :: leftSlope(:), rightSlope(:), parabolaValue(:)
 
 real, pointer :: widthAtLayerTop(:,:), depthAtLayerTop(:,:), areaThisLayer(:,:)
 real, pointer :: areaTotalBelowThisLayer(:,:), dWidth(:,:)
 real, pointer :: dDepth(:,:), angle(:,:), perimeterBelowThisLayer(:,:)
 real, pointer :: area_difference(:,:), local_difference(:,:)
 
 real :: AA, BB, CC, DD
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
 area_difference         => widthDepthData (:,:, wd_area_difference)
 local_difference        => widthDepthData (:,:, wd_local_difference)

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
    
 elsewhere ( (elemI(:,ei_geometry)  == eTriangle) .and. &
         (elemI(:,ei_elem_type) == elem_type_value    )         )
    area        = volume / length
    depth       = sqrt(abs(area/(onehalfR*(leftSlope + rightSlope))))
    hyddepth    = onehalfR * depth
    eta         = zbottom + hyddepth
    topwidth    = (leftSlope + rightSlope) * depth
    perimeter   = depth * (sqrt(oneR + leftSlope**twoR) + sqrt(oneR + rightSlope**twoR))
    hydradius   = area / perimeter
 
 endwhere
 
 do ii=1, size()
    if ( (elemI(ii,ei_geometry)  == eWidthDepth) .and. &
         (elemI(ii,ei_elem_type) == elem_type_value    )         ) then
         
            area            = volume / length
            area_difference = area
            ind = minloc(volume, 1.0)
            area            = volume / length
            hyddepth        = DD
            eta             = zbottom + hyddepth
            depth           = hyddepth
            topwidth        = 0.0
            perimeter       = 0.0
            hydradius       = area / perimeter
    endif
 enddo

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine channel_or_junction
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
