! module adjustments
!
! These are utility codes that make ad hoc adjustments to the simulation
! results. These are typically to prevent negative or zero values or
! to compensate for grid-scale oscillation.
!
!==========================================================================
!
 module adjustments
! 
    use array_index
    use data_keys
    use setting_definition
    use globals
    use utility
    
    implicit none
    
    private
    
    public :: adjust_channel_velocity_limiter
    public :: adjust_junction_branch_velocity_limit
    public :: adjust_face_dynamic_limits
    public :: adjust_for_zero_geometry
    public :: adjust_negative_volume_reset
    public :: adjust_smallvolumes
    public :: adjust_Vshaped_flowrate
    public :: adjust_zero_velocity_at_zero_volume

    integer :: debuglevel = 0
    
 contains
!
!========================================================================== 
!==========================================================================
!
 subroutine adjust_channel_velocity_limiter &
    (elemR, elemYN, elemI, &
     ei_elem_type, elemType, eYN_IsAdhocFlowrate, er_Velocity_new )
!
! Limits the velocity in channel elements to setting%Limiter%Velocity%Maximum
! In general, this should only be needed where the volumes start getting small,
! but also could accidentally mask (or delay) an instability.
!
 character(64) :: subroutine_name = 'adjust_channel_velocity_limiter'
 
 real,      target, intent(in out)  :: elemR(:,:)
 logical,           intent(in out)  :: elemYN(:,:)
 integer,           intent(in)      :: elemI(:,:)
 integer,           intent(in)      :: ei_elem_type, elemType, eYN_IsAdhocFlowrate
 integer,           intent(in)      :: er_Velocity_new
 
 real,  pointer ::  velocity(:)
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 if (setting%Limiter%Velocity%UseLimitMax) then
    velocity => elemR(:,er_Velocity_new)
    
    where ( (abs(velocity) > setting%Limiter%Velocity%Maximum) .and. &
            (elemI(:,ei_elem_type) == elemType) )
            
        velocity = sign( 0.99 * setting%Limiter%Velocity%Maximum, velocity )
        elemYN(:,eYN_IsAdhocFlowrate) = .true.
        
    endwhere 
    
 endif

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name 
 end subroutine adjust_channel_velocity_limiter
!
!========================================================================== 
!==========================================================================
!
 subroutine adjust_junction_branch_velocity_limit &
    (elemMR, elemMI)
! 
! This handles the velocity limiter only in junction branches,   
! Note that the channel elements and base junction velocity limits
! are enforced in a call to adjust_channel_velocity_limiter in the
! time marching.
!
! HACK - this could be cleaned up into 2 calls to a function, but be careful
! that we don't introduce pass-by-value.
!
 character(64) :: subroutine_name = 'adjust_junction_branch_velocity_limit'
 
 real,  target,     intent(in out)  :: elemMR(:,:)
 integer,           intent(in)      :: elemMI(:,:)
 
 real,  pointer :: velocity(:), flowrate(:), area(:)
 
 integer :: mm
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 if (setting%Limiter%Velocity%UseLimitMax) then
    do mm=1,upstream_face_per_elemM
        velocity => elemMR(:,eMr_VelocityUp(mm)) 
        flowrate => elemMR(:,eMr_FlowrateUp(mm))
        area     => elemMR(:,eMr_AreaUp(mm))
        where ( (abs(velocity) > setting%Limiter%Velocity%Maximum) .and. &
                (elemMI(:,eMi_elem_type) == eJunctionChannel) .and. &
                (elemMI(:,eMi_nfaces_u) >= mm) )
            velocity = sign( 0.99 * setting%Limiter%Velocity%Maximum, velocity )
            flowrate = velocity * area
        endwhere
    enddo
    do mm=1,dnstream_face_per_elemM
        velocity => elemMR(:,eMr_VelocityDn(mm)) 
        flowrate => elemMR(:,eMr_FlowrateDn(mm))
        area     => elemMR(:,eMr_AreaDn(mm))
        where ( (abs(velocity) > setting%limiter%velocity%Maximum)  .and. &
                (elemMI(:,eMi_elem_type) == eJunctionChannel) .and. &
                (elemMI(:,eMi_nfaces_d) >= mm) )
            velocity = sign( 0.99 * setting%Limiter%Velocity%Maximum, velocity )
            flowrate = velocity * area
        endwhere
    enddo
 endif
 
 nullify(velocity)
 nullify(flowrate)
 nullify(area)
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine adjust_junction_branch_velocity_limit
!
!========================================================================== 
!==========================================================================
!
 subroutine adjust_face_dynamic_limits &
    (faceR, faceI, volumeUp, volumeDn, facemask)
!
! ensures face velocity and areas are within limits
! 
 character(64) :: subroutine_name = 'adjust_face_dynamic_limits'
 
 real,                  intent(in out)  :: faceR(:,:)
 ! note, requires that elemR be provided separately for upstream and downstream
 ! these can be the samy arrays (e.g. in the case of elem2R for channel-channel)
 real,                  intent(in)      :: volumeUp(:), volumeDn(:)
 integer,   target,     intent(in)      :: faceI(:,:)
 logical,               intent(in)      :: facemask(:)
 
 real,  pointer :: volFrac
 integer,   pointer :: eUp(:), eDn(:)
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
! Ad hoc limit the volume that can be transported out of the upstream cell. 
 if (setting%Limiter%Flowrate%UseFaceVolumeTransport) then
    volFrac => setting%Limiter%flowrate%FaceVolumeTransport
    eUp => faceI(:,fi_Melem_u)
    eDn => faceI(:,fi_Melem_d)
    !%   for a downstream flow, limit flux from the upstream volume
    where ((facemask) .and. (faceR(:,fr_flowrate) * dt > volFrac * volumeUp(eUp)))
        faceR(:,fr_flowrate) =  volFrac * volumeUp(eUp) / dt
    endwhere
    
    where ((facemask) .and. (-faceR(:,fr_flowrate) * dt > volFrac * volumeDn(eDn)))
        faceR(:,fr_flowrate) =  volFrac * volumeDn(eDn) / dt
    endwhere

 endif
 
!%  ensure face area is not smaller than zerovalue     
 if (setting%ZeroValue%UseZeroValues) then 
    where ((facemask) .and. (faceR(:,fr_Area_d) < setting%Zerovalue%Area))
        faceR(:,fr_Area_d) = setting%Zerovalue%Area    
    endwhere
    where ((facemask) .and. (faceR(:,fr_Area_u) < setting%Zerovalue%Area))
        faceR(:,fr_Area_u) = setting%Zerovalue%Area    
    endwhere 

    where ((facemask) .and. (faceR(:,fr_Area_d) >= setting%Zerovalue%Area))
        faceR(:,fr_Velocity_d) = faceR(:,fr_Flowrate) / faceR(:,fr_Area_d)
    endwhere
    where ((facemask) .and. (faceR(:,fr_Area_u) >= setting%Zerovalue%Area))
        faceR(:,fr_Velocity_u) = faceR(:,fr_Flowrate) / faceR(:,fr_Area_u)
    endwhere
 endif
 
!%  limit high velocities
 if (setting%Limiter%Velocity%UseLimitMax) then
    where ( (facemask) .and.  &
             ( abs(faceR(:,fr_Velocity_d))  > setting%Limiter%Velocity%Maximum) )
        faceR(:,fr_Velocity_d) = sign( 0.99 * setting%Limiter%Velocity%Maximum, &
                                       faceR(:,fr_Velocity_d) )        
    endwhere
    where ( (facemask) .and.  &
             ( abs(faceR(:,fr_Velocity_u))  > setting%Limiter%Velocity%Maximum) )
        faceR(:,fr_Velocity_u) = sign( 0.99 * setting%Limiter%Velocity%Maximum, &
                                       faceR(:,fr_Velocity_u) )        
    endwhere
 endif

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine adjust_face_dynamic_limits
!
!========================================================================== 
!==========================================================================
!
 subroutine adjust_for_zero_geometry &
    (elem2R, elemMR, elem2YN, elemMYN)
!
! resets geometry to user setting%ZeroValues where the geometry is too small
! 
 character(64) :: subroutine_name = 'adjust_for_zero_geometry'
 
 real,          intent(in out)  ::  elem2R(:,:)
 real, target,  intent(in out)  ::  elemMR(:,:)
 logical,       intent(in)      ::  elem2YN(:,:), elemMYN(:,:)
 
 real, pointer  :: area(:), topwidth(:)
 
 integer :: mm
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
!%  Handle channels 
 call reset_element_for_zero_values &
    (elem2R, elem2YN, e2r_Area, e2r_Eta ,e2r_Zbottom, e2r_Topwidth,  &
     e2r_Perimeter, e2r_HydDepth, e2r_HydRadius, e2YN_IsSmallVolume)
     
!%  Handle central data for multi-branch junctions     
 call reset_element_for_zero_values &
    (elemMR, elemMYN, eMr_Area, eMr_Eta, eMr_Zbottom, eMr_Topwidth,  &
     eMr_Perimeter, eMr_HydDepth, eMr_HydRadius, eMYN_IsSmallVolume)    
     
!%  Handle branch data for multi-branch junctions
!%  Upstream branches    
 call reset_juctionbranches_for_zero_values &
    (elemMR, elemMYN, eMr_AreaUp, eMr_TopwidthUp, upstream_face_per_elemM)
    
!%  Downstream branches
 call reset_juctionbranches_for_zero_values &
    (elemMR, elemMYN, eMr_AreaDn, eMr_TopwidthDn, dnstream_face_per_elemM)
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine adjust_for_zero_geometry
!
!========================================================================== 
!==========================================================================
!
 subroutine adjust_negative_volume_reset &
    (volume)
!
! Ensures that any negative volumes are set to  setting%Zerovalue%Volume
! This is a limited routine called during time-stepping to ensure that
! we don't have a divide by zero (or small value) prior to the full 
! geometry reset.
!
 character(64) :: subroutine_name = 'adjust_negative_volume_reset'
 
 real,  intent(in out)  ::  volume(:)
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 if (setting%ZeroValue%UseZeroValues) then
    where (volume < setting%Zerovalue%Volume)
        volume = setting%Zerovalue%Volume
    endwhere
 else
    where (volume < zeroR)
        volume = zeroR
    endwhere
 endif

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine adjust_negative_volume_reset
!
!========================================================================== 
!==========================================================================
!
 subroutine adjust_smallvolumes &
    (elem2R, elem2I, elem2YN, e2r_VolumeColumn, &
     elemMR, elemMI, elemMYN, eMr_VolumeColumn)
!    
! Handle small volumes - ensures cells with depths that are at or below 
! setting%SmallVolume%DepthCutoff have geometry set to small finite values
! from setting%SmallVolume%...
!
! Note that input includes specific VolumeColumn indexes so that the adjustments
! can be made to a temporary data column during the time march.
!
 character(64) :: subroutine_name = 'adjust_smallvolumes'

 real,      intent(in out)      :: elem2R(:,:),  elemMR(:,:)
 integer,   intent(in)          :: elem2I(:,:),  elemMI(:,:)
 logical,   intent(in out)      :: elem2YN(:,:), elemMYN(:,:)
 integer,   intent(in)          :: e2r_VolumeColumn, eMr_VolumeColumn
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 if (.not. setting%SmallVolume%UseSmallVolumes) return 

!%  Identify the small volumes and update the small volume ratio
!%  which are stored in elemYN(:,isSmallVolume) and elemR(:,SmallVolumeRatio)
 call smallvolume_identification &
    (elem2R, elem2I, elem2YN, e2r_VolumeColumn, &
     elemMR, elemMI, elemMYN, eMr_VolumeColumn )
    
 call smallvolume_geometry (elem2R, elem2YN, elemMR, elemMYN)    
     
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine adjust_smallvolumes
!
!========================================================================== 
!==========================================================================
!
 subroutine adjust_Vshaped_flowrate &
    (elem2R, faceR, elem2I, elem2YN)
!
! Finds places where face flowrates and element flowrates are V-shaped
! The V is reduced by adjusting the element flowrate
!  
 character(64) :: subroutine_name = 'adjust_Vshaped_flowrate'
 
 real,      target,     intent(in out)  :: elem2R(:,:)
 real,      target,     intent(in)      :: faceR(:,:)
 integer,   target,     intent(in)      :: elem2I(:,:)
 logical,   target,     intent(in out)  :: elem2YN(:,:)
 
 real,      pointer     :: elemFlow(:), faceFlow(:), elemAdjust(:), elemVel(:), elemArea(:)
 real,      pointer     :: tscaleUp(:), tscaleDn(:)
 integer,   pointer     :: mapUp(:), mapDn(:)
 logical,   pointer     :: elemMask(:)
 
 integer :: e2r_adjustflow, e2YN_mask
 
 real,      pointer     :: coef
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 

 e2r_adjustflow = e2r_Temp(next_e2r_temparray)
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)
 
 e2YN_mask = e2YN_Temp(next_e2YN_temparray)
 next_e2YN_temparray = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

 elemFlow   => elem2R(:,e2r_Flowrate)
 elemVel    => elem2R(:,e2r_Velocity)
 elemArea   => elem2R(:,e2r_Area)
 tscaleUp   => elem2R(:,e2r_Timescale_u) 
 tscaleDn   => elem2R(:,e2r_Timescale_d) 
 elemAdjust => elem2R(:,e2r_adjustflow)
 
 elemMask   => elem2YN(:,e2YN_mask)
 elemMask = .false.
 
 faceFlow =>  faceR(:,fr_Flowrate)
 
 mapUp => elem2I(:,e2i_Mface_u)
 mapDn => elem2I(:,e2i_Mface_d)
 
 elemAdjust => elem2R(:, e2r_adjustflow)
 
 !% HACK - at this point only handling channel elements
 where (elem2I(:,e2i_elem_type) == eChannel)
    elemMask = .true.
 endwhere
 
!%  Mask is the cells where the flowrate in the element is higher than both
!%  of the faces.
 where (elemMask)
    elemMask = ( (utility_sign_with_ones(faceFlow(mapUp) - elemFlow))      &
                *(utility_sign_with_ones(faceFlow(mapDn) - elemFlow)) > 0)
 endwhere
   
!%  The adjusted element flowrate is a linear weighting of face flowrates 
!%  so that the face with the shorter timescale to reach has larger influence.
!%
!%  Arguably, this would likely be best using timescales, but we need to test
!%  this further  
! where (elemMask)  
!    elemAdjust =  (  tscaleUp * faceFlow(mapDn)   &
!                   + tscaleDn * faceFlow(mapUp) ) &
!                 / ( tscaleUp + tscaleDn )
! endwhere   
 where (elemMask)              
    elemAdjust =  (  0.5 * faceFlow(mapDn)   &
                   + 0.5 * faceFlow(mapUp) ) 
 endwhere   
    
!% apply a weighted combination based on the setting coefficient to damp the V 
 coef => setting%Method%AdjustVshapedFlowrate%Coef
 where (elemMask)   
    elemFlow = coef * elemAdjust + (oneR - coef) * elemFlow
    elemVel  = elemFlow / elemArea
 endwhere
 
!print *, elemFlow(size(elemFlow)-5:size(elemFlow))

 
!%  close up the temp arrays 
 elemAdjust = nullvalueR
 nullify(elemAdjust)
 next_e2r_temparray = next_e2r_temparray-1
 
 elemMask = nullvalueL
 nullify(elemMask)
 next_e2YN_temparray = next_e2YN_temparray-1 
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine adjust_Vshaped_flowrate
!
!========================================================================== 
!==========================================================================
!
 subroutine adjust_zero_velocity_at_zero_volume &
    (elem2R, elem2YN, e2r_VelocityColumn, e2r_VolumeColumn, &
     elemMR, elemMYN, eMr_VelocityColumn, eMr_VolumeColumn   )
!
! ensures that volumes smaller than the user limit have dynamics
! that are set to the user limit.
!
! This can create small non-zero momentum in the flow if the user minimums
! for setting%ZeroValue%Velocity, Flowrate are non-zero
!       
 character(64) :: subroutine_name = 'adjust_zero_velocity_at_zero_volume'
 
 real,      intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 logical,   intent(in out)  :: elem2YN(:,:), elemMYN(:,:)
 integer,   intent(in)      :: e2r_VelocityColumn, e2r_VolumeColumn
 integer,   intent(in)      :: eMr_VelocityColumn, eMr_VolumeColumn
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 if (setting%ZeroValue%Volume > zeroR) then
 
    call zero_velocity_at_zero_volume &
        (elem2R, elem2YN, e2r_VelocityColumn, e2r_VolumeColumn, e2YN_IsAdhocFlowrate )
    
    call zero_velocity_at_zero_volume &
        (elemMR, elemMYN, eMr_VelocityColumn, eMr_VolumeColumn, eMYN_IsAdhocFlowrate ) 
 endif
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name 
 end subroutine adjust_zero_velocity_at_zero_volume  
!
!========================================================================== 
!
! PRIVATE BELOW HERE
!
!==========================================================================
!
 subroutine smallvolume_identification &
    (elem2R, elem2I, elem2YN, eTr_Volume2, &
     elemMR, elemMI, elemMYN, eTr_VolumeM )
!
! updates small volume identification for different element types
!
 character(64) :: subroutine_name = 'smallvolume_identification'
 
 real,    target,   intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 integer, target,   intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical, target,   intent(in out)  :: elem2YN(:,:), elemMYN(:,:)
 integer,           intent(in)      :: eTr_Volume2,  eTr_VolumeM 
 
 real,      pointer :: smallvolumeratio(:), smallvolume(:), tvolume(:)
 integer,   pointer :: elemtype(:)
 
 logical, pointer :: issmallvolume(:)
 
 integer :: thiselemtype
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 

 if (.not. setting%SmallVolume%UseSmallVolumes) return
 
!%  for channel elements
 smallvolumeratio   => elem2R(:,e2r_SmallvolumeRatio)
 smallvolume        => elem2R(:,e2r_Smallvolume)
 issmallvolume      => elem2YN(:,e2YN_IsSmallVolume) 
 tvolume            => elem2R(:,eTr_Volume2)
 elemtype           => elem2I(:,e2i_elem_type)
 thiselemtype = eChannel
    
 call smallvolume_identification_for_element &
    (tvolume, smallvolumeratio, smallvolume, issmallvolume, &
     elemtype, thiselemtype)

!%  for junction elements  
 smallvolumeratio   => elemMR(:,eMr_SmallvolumeRatio)
 smallvolume        => elemMR(:,eMr_Smallvolume) 
 issmallvolume      => elemMYN(:,eMYN_IsSmallVolume)
 tvolume            => elemMR(:,eTr_VolumeM)
 elemtype           => elemMI(:,eMi_elem_type)
 thiselemtype = eJunctionChannel 
 
 call smallvolume_identification_for_element &
    (tvolume, smallvolumeratio, smallvolume, issmallvolume, &
     elemtype, thiselemtype)  
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine smallvolume_identification
!
!========================================================================== 
!========================================================================== 
!
 subroutine smallvolume_identification_for_element &
    (volume, smallvolumeratio, smallvolume, issmallvolume, &
     elemtype, thiselementtype )
!
! updates smallvolume identification for a single element type
!    
 character(64) :: subroutine_name = 'smallvolume_identification_for_element'
 
 real,      intent(in out)  :: volume(:), smallvolumeratio(:)
 real,      intent(in)      :: smallvolume(:)
 integer,   intent(in)      :: elemtype(:)
 logical,   intent(in out)  :: issmallvolume(:)
 integer,   intent(in)      :: thiselementtype
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 if (.not. setting%SmallVolume%UseSmallVolumes) return
  
!%  zero out small volumes (reset further below to setting%Zerovalue%Volume
!%  purpose is to ensure that volumes below the zero value will only use
!%  the Chezy-Manning approximation for flow.
!%  NOTE: this is a source of a small volume non-conservation of mass
 where ((volume < setting%Zerovalue%Volume) .and. (elemtype == thiselementtype))
    volume = zeroR
 endwhere
 
!%  Compute smallvolume_ratio where the local volume is small
!%  The ratio fraction used as a blending function in small volumes
!%  and is 1.0 where the element volume is exactly the small volume
 where ((volume < smallvolume) .and. (elemtype == thiselementtype))
    issmallvolume = .true.
    smallvolumeratio = volume / smallvolume
 endwhere
 
!%  reset zero volumes to small volume value
 where ((volume <= setting%Zerovalue%Volume) .and. (elemtype == thiselementtype))
    volume = setting%Zerovalue%Volume
 endwhere
  
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine smallvolume_identification_for_element   
!
!========================================================================== 
!========================================================================== 
!
 subroutine smallvolume_geometry &
    (elem2R, elem2YN, &
     elemMR, elemMYN)
!
! adjusts geometry that are designated as "smallvolume", which are derived
! from the setting%SmallDepth 
!
 character(64) :: subroutine_name = 'smallvolume_geometry'

 real,          intent(in out)  :: elem2R(:,:)
 real,  target, intent(in out)  :: elemMR(:,:)
 logical,       intent(in)      :: elem2YN(:,:), elemMYN(:,:)
 
 real,  pointer :: area(:), topwidth(:)
 integer :: mm
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 if (.not. setting%SmallVolume%UseSmallVolumes) return
 
!%  channel elements
 call smallvolume_element_geometry_reset &
    (elem2R, elem2YN, &
     e2r_Area, e2r_Eta, e2r_Perimeter, e2r_Zbottom, e2r_HydDepth, e2r_HydRadius, &
     e2r_Topwidth, e2YN_IsSmallVolume)

!%  junction elements
 call smallvolume_element_geometry_reset &
    (elemMR, elemMYN, &
     eMr_Area, eMr_Eta, eMr_Perimeter, eMr_Zbottom, eMr_HydDepth, eMr_HydRadius, &
     eMr_Topwidth, eMYN_IsSmallVolume) 

!%  junction branches  
 call smallvolume_junctionbranch_reset &
    (elemMR, elemMYN, eMr_AreaUp, eMr_TopwidthUp, upstream_face_per_elemM)  

 call smallvolume_junctionbranch_reset &
    (elemMR, elemMYN, eMr_AreaDn, eMr_TopwidthDn, dnstream_face_per_elemM)  

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine smallvolume_geometry
!
!========================================================================== 
!========================================================================== 
!
 subroutine smallvolume_element_geometry_reset &
    (elemR, elemYN, &
     er_Area, er_Eta, er_Perimeter, er_Zbottom, er_HydDepth, er_HydRadius, &
     er_Topwidth, eYN_IsSmallVolume)
!
! for small volumes, reset the geometry to the setting%smallvolume values
!   
 character(64) :: subroutine_name = 'smallvolume_element_geometry_reset'
 
 real,  target, intent(in out)  :: elemR(:,:)

 logical,   intent(in)  :: elemYN(:,:)
 
 integer,   intent(in)  :: er_Area, er_Eta, er_Perimeter, er_Zbottom
 integer,   intent(in)  :: er_HydDepth, er_HydRadius, er_Topwidth
 integer,   intent(in)  :: eYN_IsSmallVolume
 
 real, pointer  ::  area(:), eta(:), perimeter(:), zbottom(:)
 real, pointer  ::  hyddepth(:), hydradius(:)
 real, pointer  ::  topwidth(:) 

!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

 if (.not. setting%SmallVolume%UseSmallVolumes) return
 
 area               => elemR(:,er_Area)
 perimeter          => elemR(:,er_Perimeter)
 hyddepth           => elemR(:,er_Hyddepth)
 hydradius          => elemR(:,er_Hydradius)
 topwidth           => elemR(:,er_Topwidth)
 eta                => elemR(:,er_Eta)
 zbottom            => elemR(:,er_Zbottom)
 
 where (elemYN(:,eYN_IsSmallVolume))
    area        = setting%SmallVolume%MinimumArea
    perimeter   = setting%SmallVolume%MinimumTopwidth
    hyddepth    = setting%SmallVolume%DepthCutoff
    hydradius   = setting%SmallVolume%MinimumHydRadius
    topwidth    = setting%SmallVolume%MinimumTopwidth
    eta         = setting%SmallVolume%DepthCutoff + zbottom
 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine smallvolume_element_geometry_reset
!
!========================================================================== 
!==========================================================================
!
 subroutine smallvolume_junctionbranch_reset &
    (elemMR, elemMYN, eMr_AreaDir, eMr_TopwidthDir, &
     Dir_face_per_elem)
!
! reset the values for branches in a junction that has a smallvolume
! 
 character(64) :: subroutine_name = 'smallvolume_junctionbranch_reset'
 
 real,      target, intent(in out)  :: elemMR(:,:)
 logical,           intent(in)      :: elemMYN(:,:)
 
 ! dnstream_face_per_elemM or upstream_face_per_elemM
 integer,           intent(in)      :: Dir_face_per_elem 
 integer,           intent(in)      :: eMr_AreaDir(:), eMr_TopwidthDir(:)
 
 real,  pointer :: area(:), topwidth(:)
 integer        :: mm
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 if (.not. setting%SmallVolume%UseSmallVolumes) return
 
 do mm=1,Dir_face_per_elem
    area     => elemMR(:,eMr_AreaDir(mm))
    topwidth => elemMR(:,eMr_TopwidthDir(mm))   
    where (elemMYN(:,eMYN_IsSmallVolume))
        topwidth = setting%SmallVolume%MinimumTopwidth
        area     = setting%SmallVolume%MinimumArea
    endwhere    
 enddo
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine smallvolume_junctionbranch_reset
!
!========================================================================== 
!========================================================================== 
!
 subroutine reset_element_for_zero_values &
    (elemR, elemYN, er_Area, er_Eta, er_Zbottom, er_Topwidth, er_Perimeter, &
     er_HydDepth, er_HydRadius, eYN_IsSmallVolume)
!
! Adjusts miscellaneous geometry for small or negative values
! Note thismask is true only where elemYN(:,IsSmallVolume) is false). 
! That is, this only applies to elements NOT covered by the small volume 
! algorithm.
!
! In general, this should not be needed, but is provided in case the
! setting%SmallVolume algorithm is turned off or if some degenerate case
! occurs where the volume is large but the element geometry is small.
!    
 character(64) :: subroutine_name = 'reset_element_for_zero_values'
 
 real,  target,  intent(in out) :: elemR(:,:)
 logical,        intent(in)     :: elemYN(:,:)
 
 integer,        intent(in)     :: er_Area, er_Eta, er_Zbottom, er_Topwidth
 integer,        intent(in)     :: er_Perimeter, er_HydDepth, er_HydRadius
 integer,        intent(in)     :: eYN_IsSmallVolume
 
 real,  pointer :: area(:), eta(:), topwidth(:), perimeter(:) 
 real,  pointer :: zbottom(:), hyddepth(:), hydradius(:)
 
 integer :: ii

!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

 area       => elemR(:,er_Area)
 eta        => elemR(:,er_Eta)
 topwidth   => elemR(:,er_Topwidth)
 perimeter  => elemR(:,er_Perimeter)
 hyddepth   => elemR(:,er_HydDepth)
 hydradius  => elemR(:,er_HydRadius)
 zbottom    => elemR(:,er_Zbottom)
 
 if (setting%ZeroValue%UseZeroValues) then
    where (.not. elemYN(:,eYN_IsSmallVolume))
        where (area < setting%Zerovalue%Area) 
            area = setting%Zerovalue%Area
        endwhere

        where (eta - zbottom < setting%Zerovalue%Depth) 
            eta = zbottom + setting%Zerovalue%Depth
        endwhere
        
        where (topwidth < setting%Zerovalue%Topwidth) 
            topwidth = setting%Zerovalue%Topwidth
        endwhere
        
        where (perimeter < setting%Zerovalue%Topwidth) 
            perimeter = setting%Zerovalue%Topwidth
        endwhere
        
        where (hyddepth < setting%Zerovalue%Depth)
            hyddepth = setting%Zerovalue%Depth
        endwhere
        
        where (hydradius < setting%Zerovalue%Depth)
            hydradius = setting%Zerovalue%Depth
        endwhere
    endwhere
 else
    where (.not. elemYN(:,eYN_IsSmallVolume))
        where (area < zeroR) 
            area = zeroR
        endwhere

        where (eta - zbottom < zeroR) 
            eta = zbottom 
        endwhere
        
        where (topwidth < zeroR) 
            topwidth = zeroR
        endwhere
        
        where (perimeter < zeroR) 
            perimeter = zeroR
        endwhere
        
        where (hyddepth < zeroR)
            hyddepth = zeroR
        endwhere
        
        where (hydradius < zeroR)
            hydradius = zeroR
        endwhere
    endwhere
 endif
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine reset_element_for_zero_values   
!
!========================================================================== 
!==========================================================================
!
 subroutine reset_juctionbranches_for_zero_values &
    (elemMR, elemMYN, eMr_AreaDir, eMr_TopwidthDir, Dir_face_per_elemM)
!
! resets values below the zero value to the zero value in
! junction branches. Applied separately for upstream and downstream branches.
! 
 character(64) :: subroutine_name = 'reset_juctionbranches_for_zero_values'
 
 real,      target,     intent(in out)  :: elemMR(:,:)
 logical,               intent(in)      :: elemMYN(:,:)
 integer,               intent(in)      :: Dir_face_per_elemM
 integer,               intent(in)      :: eMr_AreaDir(:), eMr_TopwidthDir(:)
 real,  pointer :: area(:), topwidth(:)
 integer        :: mm
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 do mm=1,Dir_face_per_elemM
    area     => elemMR(:,eMr_AreaDir(mm))
    topwidth => elemMR(:,eMr_TopwidthDir(mm))
    if (setting%ZeroValue%UseZeroValues) then
        where (.not. elemMYN(:,eMYN_IsSmallVolume))
            where (area < setting%Zerovalue%Area)
                area = setting%Zerovalue%Area
            endwhere    
            where (topwidth < setting%Zerovalue%Topwidth)
                topwidth = setting%Zerovalue%Topwidth
            endwhere   
        endwhere
    else
        where (.not. elemMYN(:,eMYN_IsSmallVolume))
            where (area < zeroR)
                area = zeroR
            endwhere    
            where (topwidth < setting%Zerovalue%Topwidth)
                topwidth = zeroR
            endwhere   
        endwhere
    endif
 enddo
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine reset_juctionbranches_for_zero_values
!
!========================================================================== 
!==========================================================================
!
 subroutine zero_velocity_at_zero_volume &
    (elemR, elemYN, er_Velocity, er_Volume, eYN_IsAdhocFlowrate )
!
! sets velocity to minimal values if the volume is less
! than the user minimum. Note that this should only apply in volumes that
! are << elemR(:,Smallvolume)
!
 character(64) :: subroutine_name = 'zero_velocity_at_zero_volume'
 
 real,      target,     intent(in out)  :: elemR(:,:)
 logical,               intent(in out)  :: elemYN(:,:)
 integer,               intent(in)      :: er_Velocity, er_Volume, eYN_IsAdhocFlowrate
 
 real,  pointer :: velocity(:),  volume(:)
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 if (setting%ZeroValue%Volume > 0.0) then
    velocity => elemR(:,er_Velocity)
    volume   => elemR(:,er_Volume)
    where (volume <= setting%Zerovalue%Volume)
        velocity = sign(setting%Zerovalue%Velocity, velocity)
        elemYN(:,eYN_IsAdhocFlowrate) = .true.
    endwhere
 endif
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine zero_velocity_at_zero_volume
!
!========================================================================== 
! END OF MODULE adjustments
!========================================================================== 
 end module adjustments