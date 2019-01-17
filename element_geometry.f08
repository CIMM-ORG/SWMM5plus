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
    use setting_definition
    use globals
    
    implicit none
    
    private
    
    public :: element_geometry_update
    
    integer :: debuglevel = 0
    
 contains
!    
!========================================================================== 
!==========================================================================
!
 subroutine element_geometry_update &
    (elem2R, elem2I, elem2YN, e2r_VolumeColumn, &
     elemMR, elemMI, elemMYN, eMr_VolumeColumn, &
     faceI, bcdataDn, bcdataUp, thisTime  )
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
 type(bcType),      intent(in out) :: bcdataDn(:), bcdataUp(:)
  real,             intent(in)     :: thisTime
    
 integer,   intent(in)  :: e2r_VolumeColumn, eMr_VolumeColumn 
 integer, parameter :: idummy = 0
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
!% reset values for small volume ahndling
 if (setting%SmallVolume%UseSmallVolumes) then
    elem2YN(:,e2YN_IsSmallVolume) = .false.
    elemMYN(:,eMYN_IsSmallVolume) = .false.
    elem2R(:,e2r_SmallVolumeRatio) = nullvalueR 
    elemMR(:,eMr_SmallVolumeRatio) = nullvalueR 
 endif
 
!% apply any elevation bc and fix the element volume 
 call bc_applied_onelement &
    (elem2R, bcdataDn, bcdataUp, thisTime, bc_category_elevation, idummy) 
    
!% rectangular geometry  
 call rectangular_geometry_update &
    (elem2R, elem2I, e2r_VolumeColumn, &
     elemMR, elemMI, eMr_VolumeColumn  )

!% HACK -- NEED OTHER GEOMETRY TYPES

!% reset the computed geometry values where volumes are small
 if (setting%SmallVolume%UseSmallVolumes) then
    call adjust_smallvolumes &
        (elem2R, elem2I, elem2YN, e2r_VolumeColumn, &
         elemMR, elemMI, elemMYN, eMr_VolumeColumn    )
 endif

!% reset the geometry (non-volume) where values are below minimums  
 call adjust_for_zero_geometry (elem2R, elemMR, elem2YN, elemMYN)

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine element_geometry_update
! 
!========================================================================== 
!
! PRIVATE BELOW HERE
!
!==========================================================================
!
 subroutine rectangular_geometry_update &
    (elem2R, elem2I, e2r_Volume_new, &
     elemMR, elemMI, eMr_Volume_new )
!
! Note that volume used is in a eTr storage location so that the update
! can be used on a temporary volume
!    
 character(64) :: subroutine_name = 'rectangular_geometry_update'
 
 real,      intent(in out)  :: elem2R(:,:), elemMR(:,:)
 integer,   intent(in)      :: elem2I(:,:), elemMI(:,:)

 integer,   intent(in)      :: e2r_Volume_new, eMr_Volume_new 

!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name
 
!%  basic geometry update for rectangular channels and junctions 
 
!%  rectangular geometry for channels
 call rectangular_channel_or_junction &
    (elem2R, elem2I, &
     e2i_geometry, e2i_elem_type, eChannel,  &
     e2r_Length, e2r_Zbottom, e2r_BreadthScale, e2r_Area, e2r_Eta, &
     e2r_Perimeter, e2r_HydDepth, e2r_HydRadius, e2r_Volume_new)

!%  rectangular geomety for junctions
 call rectangular_channel_or_junction &
    (elemMR, elemMI, &
     eMi_geometry, eMi_elem_type, eJunctionChannel,  &
     eMr_Length, eMr_Zbottom, eMr_BreadthScale, eMr_Area, eMr_Eta, &
     eMr_Perimeter, eMr_HydDepth, eMr_HydRadius, eMr_Volume_new)
     
!%  geometry for junction branches (area only - updated from eta)  
!%  upstream branches 
 call rectangular_junction_leg &
    (elemMR, elemMI, upstream_face_per_elemM, eMi_nfaces_u,   &
     eMr_AreaUp, eMr_ZbottomUp, eMr_BreadthScaleUp)
     
!%  downstream branches
 call rectangular_junction_leg &
    (elemMR, elemMI, dnstream_face_per_elemM, eMi_nfaces_d,   &
     eMr_AreaDn, eMr_ZbottomDn, eMr_BreadthScaleDn)
     
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine rectangular_geometry_update 
! 
!========================================================================== 
!========================================================================== 
!
 subroutine rectangular_channel_or_junction &
    (elemR, elemI, &
     ei_geometry, ei_elem_type, elem_typ_value,  &
     er_Length, er_Zbottom, er_BreadthScale, er_Area, er_Eta, &
     er_Perimeter, er_HydDepth, er_HydRadius, er_Volume)
!
! computes element geometry for a rectangular channel or a channeljunction
! 
 character(64) :: subroutine_name = 'rectangular_channel_or_junction'
 
 real,      target,     intent(in out)  :: elemR(:,:)
    
 integer,   intent(in)      :: elemI(:,:)  
 integer,   intent(in)      :: ei_geometry, ei_elem_type, elem_typ_value
 integer,   intent(in)      :: er_Length, er_Zbottom, er_BreadthScale
 integer,   intent(in)      :: er_Area, er_Eta, er_Perimeter
 integer,   intent(in)      :: er_HydDepth, er_HydRadius, er_Volume
 
 real,  pointer  :: volume(:), length(:), zbottom(:), breadth(:)
 real,  pointer  :: area(:), eta(:), perimeter(:), hyddepth(:), hydradius(:)
 
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 ! inputs
 volume     => elemR(:,er_Volume)
 length     => elemR(:,er_Length)
 zbottom    => elemR(:,er_Zbottom) 
 breadth    => elemR(:,er_BreadthScale)
 
! outputs
 area       => elemR(:,er_Area)
 eta        => elemR(:,er_Eta)
 perimeter  => elemR(:,er_Perimeter)
 hyddepth   => elemR(:,er_HydDepth)
 hydradius  => elemR(:,er_HydRadius)
 
 where ( (elemI(:,ei_geometry)  == eRectangular) .and. &
         (elemI(:,ei_elem_type) == elem_typ_value    )         )
    area        = volume / length
    eta         = zbottom + (area / breadth)
    perimeter   = breadth + 2.0 * ( eta - zbottom )
    hyddepth    = area / breadth
    hydradius   = area / perimeter
 endwhere
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine rectangular_channel_or_junction
!
!========================================================================== 
!==========================================================================
!
 subroutine rectangular_junction_leg &
    (elemMR, elemMI, eMnumberDir, eMi_nfacesDir, eMr_AreaDir,  &
     eMr_ZbottomDir, eMr_BreadthScaleDir)
!
! computes geometry in the individual junction branches for a channel
! element with multiple connections
! 
 character(64) :: subroutine_name = 'rectangular_junction_leg'
 
 real,      target,     intent(in out)  :: elemMR(:,:) 
 integer,               intent(in)      :: elemMI(:,:)
 
 integer,               intent(in)      :: eMnumberDir, eMi_nfacesDir
 integer,               intent(in)      :: eMr_AreaDir(:)  
 integer,               intent(in)      :: eMr_ZbottomDir(:) 
 integer,               intent(in)      :: eMr_BreadthScaleDir(:) 
 
 integer    :: mm
    
 real,      pointer   :: eta(:), area(:), zbottom(:), breadth(:)   
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 eta => elemMR(:,eMr_Eta)
 
 !% cycle over branches. Dir is either {Up, Dn}
 do mm=1,eMnumberDir
    area        => elemMR(:,eMr_AreaDir(mm))
    zbottom     => elemMR(:,eMr_ZbottomDir(mm))
    breadth     => elemMR(:,eMr_BreadthScaleDir(mm))
    
    !% HACK - rectangular geometry only
    where ( (elemMI(:,eMi_geometry)  == eRectangular)     .and. &
            (elemMI(:,eMi_elem_type) == eJunctionChannel) .and. &
            (elemMI(:,eMi_nfacesDir) >= mm)  )
        area = (eta - zbottom) * breadth
    endwhere
 enddo
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine rectangular_junction_leg
! 
!========================================================================== 
! END OF MODULE element_geometry
!========================================================================== 
 end module element_geometry