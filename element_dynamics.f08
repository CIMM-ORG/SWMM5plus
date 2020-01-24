! module element_dynamics
!
! Updating values on elements during simulation time stepping
!
!==========================================================================
!
 module element_dynamics
! 
    use array_index
    use adjustments
    use bc
    use data_keys
    use globals
    use junction
    use setting_definition
    use utility
    
    implicit none
    
    private
    
    public  :: element_dynamics_update 

    integer :: debuglevel = 0
    
 contains
! 
!========================================================================== 
!==========================================================================
!
 subroutine element_dynamics_update &
    (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
     bcdataDn, bcdataUp, e2r_Velocity_new, eMr_Velocity_new, &
     e2r_Volume_new, eMr_Volume_new, thisTime)
!
! update the flow dynamics on an element given new velocity values stored
! in the array given by the column index e#r_Velocity_new 
!
 character(64) :: subroutine_name = 'element_dynamics_update'
 
 real,      target,     intent(in out)  ::  elem2R(:,:),  elemMR(:,:)
 real,      target,     intent(in)      ::  faceR(:,:)
 integer,               intent(in)      ::  elem2I(:,:),  elemMI(:,:)
 logical,   target,     intent(in out)  ::  elem2YN(:,:), elemMYN(:,:)
 type(bcType),          intent(in out)  ::  bcdataDn(:),  bcdataUp(:)
 real,                  intent(in)      ::  thisTime
 integer,               intent(in)      ::  e2r_Velocity_new, e2r_Volume_new
 integer,               intent(in)      ::  eMr_Velocity_new, eMr_Volume_new
  
 logical,   pointer :: isadhocflowrate(:) 
 integer            :: idx
    
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 

!% velocity limiter on channels
 call adjust_channel_velocity_limiter &
    (elem2R, elem2YN, elem2I, &
     e2i_elem_type, eChannel, e2YN_IsAdhocFlowrate, e2r_Velocity_new )

!% velocity limiter on junction (handled on main element only)     
 call adjust_channel_velocity_limiter &
    (elemMR, elemMYN, elemMI, &
     eMi_elem_type, eJunctionChannel, eMYN_IsAdhocFlowrate,  eMr_Velocity_new)    
    
!%  For small volumes, compute a velocity that is blended from the update value
!%  and a Chezy-Manning computed using the free surface slope of the element 
 if (setting%SmallVolume%UseSmallVolumes) then 
    call blended_smallvolume_velocity &
        (elem2R, elem2I, elem2YN, elemMR, elemMI, elemMYN, faceR,  &
        e2r_Velocity_new, eMr_Velocity_new)
 endif
    
!% for extremely small volumes - perform a separate reset
 if (setting%ZeroValue%Volume > zeroR) then
    call adjust_zero_velocity_at_zero_volume &
        (elem2R, elem2YN, e2r_Velocity_new, e2r_Volume_new, &
         elemMR, elemMYN, eMr_Velocity_new, eMr_Volume_new)   
 endif

!%  flowrate updated from velocity 
 call element_flowrate_update  &
     (elem2R, elemMR, faceR, elem2I, elemMI, e2r_Velocity_new, eMr_Velocity_new)
     
!%  apply the boundary conditions on velocity and flowrate     
 call bc_applied_onelement &
    (elem2R, bcdataDn, bcdataUp, thisTime, bc_category_inflowrate, e2r_Velocity_new)    
    
!% compute the timescales up and down
 call element_timescale &
    (elem2R, elem2I, elem2YN, elemMR, elemMI, elemMYN, bcdataDn, bcdataUp, &
     e2r_Velocity_new)
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine element_dynamics_update
!
!========================================================================== 
!
! PRIVATE BELOW HERE
!
!========================================================================== 
!
 subroutine element_flowrate_update &
    (elem2R, elemMR, faceR, elem2I, elemMI, &
     e2r_Velocity_new, eMr_Velocity_new)
!
! given a velocity and update areas in elemR, provide updates
! for flow rate of the element and junction branches
! 
 character(64) :: subroutine_name = 'element_flowrate_update'
 
 real,      target, intent(in out)  ::  elem2R(:,:), elemMR(:,:)
 real,      target, intent(in)      ::  faceR(:,:)
 integer,   target, intent(in)      ::  elem2I(:,:), elemMI(:,:)
 integer,           intent(in)      ::  e2r_Velocity_new, eMr_Velocity_new
 
 real,      pointer :: totalarea(:)
 
 integer :: mm, eMR_totalarea
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 

!%  assign total area from the temp space
 eMR_totalarea = eMr_Temp(next_eMr_temparray)
 totalarea  => elemMR(:,eMR_totalarea)
 next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

!%  update the channel flow rates with the velocity (may be from an RK step)    
 call flowrate_from_velocity &
    ( elem2R, elem2I, &
      e2r_Flowrate, e2r_Area, e2r_Velocity_new, e2i_elem_type, eChannel)
 
 call flowrate_from_velocity &
    ( elemMR, elemMI, &
      eMr_Flowrate, eMr_Area, eMr_Velocity_new, eMi_elem_type, eJunctionChannel)

 call flowrate_from_velocity &
    ( elem2R, elem2I, &
      e2r_Flowrate, e2r_Area, e2r_Velocity_new, e2i_elem_type, eWeir)


!%  FLOWS AND VELOCITIES IN JUNCTION BRANCHES -----------------
!%  The total flowrate is distributed among both outflow and inflow branches
!%  Note that this does NOT mean that the inflows and outflows are exactly 
!%  equal - the branch values are what are used in the face interpolation 
!%  to get the inflow/outflows at the faces.

!%  get the total outflow branch areas
 
!%  handle the downstream branches 
 totalarea = zeroR ! ensure temporary array is zero
 call junction_branch_sum_areas_by_direction &
    (eMR_totalarea, &
     dnstream_face_per_elemM, eMr_AreaDn, eMi_MfaceDn, eMi_nfaces_d, &
     upstream_face_per_elemM, eMr_AreaUp, eMi_MfaceUp, eMi_nfaces_u, &
     elemMR, elemMI, faceR)
       
!%  distribute flow proportionally among outflows downstream
 call junction_branch_velocity_and_flowrate_proportional_to_area &
    (eMR_totalarea, eMR_Flowrate, &
     dnstream_face_per_elemM, eMr_AreaDn, eMi_MfaceDn, eMi_nfaces_d, &
     eMr_FlowrateDn, eMr_VelocityDn, &
     upstream_face_per_elemM, eMr_AreaUp, eMi_MfaceUp, eMi_nfaces_u, &
     eMr_FlowrateUp, eMr_VelocityUp, &
     elemMR, elemMI, faceR)    

!%  handle upstream branches
 totalarea = zeroR ! ensure temporary array is zero          
 call junction_branch_sum_areas_by_direction &
    (eMR_totalarea, &
     upstream_face_per_elemM, eMr_AreaUp, eMi_MfaceUp, eMi_nfaces_u, &
     dnstream_face_per_elemM, eMr_AreaDn, eMi_MfaceDn, eMi_nfaces_d, &
     elemMR, elemMI, faceR)
     
!%  distribute flow proportionally among outflows upstream
 call junction_branch_velocity_and_flowrate_proportional_to_area &
    (eMR_totalarea, eMR_Flowrate, &
     upstream_face_per_elemM, eMr_AreaUp, eMi_MfaceUp, eMi_nfaces_u, &
     eMr_FlowrateUp, eMr_VelocityUp, &
     dnstream_face_per_elemM, eMr_AreaDn, eMi_MfaceDn, eMi_nfaces_d, &
     eMr_FlowrateDn, eMr_VelocityDn, &
     elemMR, elemMI, faceR)         
     
     
!print *, trim(subroutine_name)
!print *, elemMR(:,eMr_totalarea)
!print *, elemMR(:,eMr_Flowrate_u1), elemMR(:,eMr_Velocity_u1 )
!print *, elemMR(:,eMr_Flowrate_u2), elemMR(:,eMr_Velocity_u2 ) 
!print *, elemMR(:,eMr_Flowrate),    elemMR(:,eMr_Velocity_new)     
!print *, elemMR(:,eMr_Flowrate_d1), elemMR(:,eMr_Velocity_d1 )  

!%  enforce maximum velocities in junction branches
 call adjust_junction_branch_velocity_limit (elemMR, elemMI)

!%  release the temp array
 totalarea = nullvalueR
 nullify(totalarea)
 next_eMr_temparray = next_eMr_temparray-1

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine element_flowrate_update
!
!========================================================================== 
!==========================================================================
!
 subroutine flowrate_from_velocity &
    ( elemR, elemI, &
      er_Flowrate, er_Area, er_Velocity, ei_elem_type, eThisElemType)
!
! compute flowrate from velocity and area for eThisElemType
! 
 character(64) :: subroutine_name = 'flowrate_from_velocity'
 
 real,  target, intent(in out)  :: elemR(:,:)
 integer,       intent(in)      :: elemI(:,:)
 
 integer,       intent(in)  :: er_Flowrate, er_Area, er_Velocity, ei_elem_type
 integer,       intent(in)  :: eThisElemType
 
 real,  pointer :: flowrate(:), area(:), velocity(:)
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 flowrate => elemR(:,er_Flowrate)
 area     => elemR(:,er_Area)
 velocity => elemR(:,er_Velocity)
 
 where (elemI(:,ei_elem_type) == eThisElemType)
    flowrate = area * velocity
 endwhere
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine flowrate_from_velocity
!
!========================================================================== 
!==========================================================================
!
 subroutine blended_smallvolume_velocity &
    (elem2R, elem2I, elem2YN, elemMR, elemMI, elemMYN, faceR,  &
     e2r_Velocity_new, eMr_Velocity_new)
!
! any volume IsSmallVolume will be have its velocity blended with
! a Chezy velocity so that the small volumes reduce to a small finite value
! as the volume decreases.
!
! Note that the direction of the blended velocity will follow the 
! decrease in the free surface slope.
! 
 character(64) :: subroutine_name = 'blended_smallvolume_velocity'
 
 real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 real,      target, intent(in)      :: faceR(:,:)
 integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,           intent(in)      :: elem2YN(:,:), elemMYN(:,:)
 
 integer, intent(in) :: e2r_Velocity_new, eMr_Velocity_new
    
 integer :: e2r_tSlope, e2r_tManningsN, e2r_tSmallVelocity
 integer :: eMr_tSlope, eMr_tManningsN, eMr_tSmallVelocity
  
 integer :: mm
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 

!%  velocity blend for channel elements 
 call velocity_blend_with_mask &
    (elem2R, elem2I, elem2YN, faceR, &
     next_e2r_temparray, e2r_n_temp, e2r_Temp, &
     e2r_Velocity_new, e2r_Flowrate, e2r_Length,  &
     e2r_Area, e2r_HydRadius, e2r_SmallVolumeRatio, &
     e2YN_IsSmallVolume, e2i_roughness_type, e2r_Roughness, &
     e2i_elem_type, eChannel, e2i_Mface_u, e2i_Mface_d)

!%  velocity blend for junction elements 
!%  HACK - this arbitrarily uses the u1 and d1 for the slope
 call velocity_blend_with_mask &
    (elemMR, elemMI, elemMYN, faceR, &
     next_eMr_temparray, eMr_n_temp, eMr_Temp, &
     eMr_Velocity_new, eMr_Flowrate, eMr_Length,  &
     eMr_Area, eMr_HydRadius, eMr_SmallVolumeRatio, &
     eMYN_IsSmallVolume, eMi_roughness_type, eMr_Roughness, &
     eMi_elem_type, eJunctionChannel, eMi_Mface_u1, eMi_Mface_d1)
 
!%  perform small volume velocity bend for junction branches
 do mm=1,upstream_face_per_elemM
    !%  HACK this arbitrarily uses d1 for the slope
    !%  Not sure how to revise if d2 and d3 are non-zero
    call velocity_blend_with_mask &
        (elemMR, elemMI, elemMYN, faceR, &
         next_eMr_temparray, eMr_n_temp, eMr_Temp, &
         eMr_VelocityUp(mm), eMr_FlowrateUp(mm), eMr_LengthUp(mm),  &
         eMr_AreaUp(mm), eMr_HydRadius, eMr_SmallVolumeRatio, &
         eMYN_IsSmallVolume, eMi_roughness_type, eMr_Roughness, &
         eMi_elem_type, eJunctionChannel, eMi_MfaceUp(mm), eMi_Mface_d1)
 end do
 
 do mm=1,dnstream_face_per_elemM
    !%  HACK this arbitrarily uses u1 for the slope
    !%  PERHAPS REVISE TO USED Z at CENTER OF CHANNEL
    call velocity_blend_with_mask &
        (elemMR, elemMI, elemMYN, faceR, &
         next_eMr_temparray, eMr_n_temp, eMr_Temp, &
         eMr_VelocityDn(mm), eMr_FlowrateDn(mm), eMr_LengthDn(mm),  &
         eMr_AreaDn(mm), eMr_HydRadius, eMr_SmallVolumeRatio, &
         eMYN_IsSmallVolume, eMi_roughness_type, eMr_Roughness, &
         eMi_elem_type, eJunctionChannel, eMi_Mface_u1, eMi_MfaceDn(mm))
 end do

! HACK - the above applies a ManningsN approach to elements that might have
! other drag approaches.
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine blended_smallvolume_velocity
!
!========================================================================== 
!==========================================================================
!
 subroutine velocity_blend_with_mask &
    (elemR, elemI, elemYN, faceR, &
     next_er_temparray, er_n_temp, er_Temp, &
     er_Velocity_new, er_Flowrate, er_Length,  &
     er_Area, er_HydRadius, er_SmallVolumeRatio, &
     eYN_IsSmallVolume, ei_roughness_type, er_Roughness, &
     ei_elem_type, elemType, ei_Mface_u, ei_Mface_d)
!
! Performs the velocity blending between a Chezy-Manning solution and
! the simulated solution for small volumes
!
 character(64) :: subroutine_name = 'velocity_blend_with_mask'

 real,      target,     intent(in out)  :: elemR(:,:)
 real,                  intent(in)      :: faceR(:,:)
 integer,   target,     intent(in)      :: elemI(:,:)
 logical,               intent(in)      :: elemYN(:,:)
 
 integer,               intent(in out)  :: next_er_temparray
 integer,               intent(in)      :: er_n_temp, er_Temp(:)
    
 integer, intent(in) :: er_Velocity_new, er_Flowrate, er_Length
 integer, intent(in) :: er_HydRadius, er_SmallVolumeRatio, er_Area 
 integer, intent(in) :: eYN_IsSmallVolume, ei_roughness_type, er_Roughness
 integer, intent(in) :: ei_elem_type, elemType, ei_Mface_u, ei_Mface_d

 integer  :: er_tSlope, er_tSmallVelocity, er_tManningsN
 
 integer,   pointer :: fUp(:), fDn(:)
 real,      pointer :: Length(:), Velocity(:), Flowrate(:), HydRadius(:)
 real,      pointer :: SmallVolumeRatio(:), Area(:)
 real,      pointer :: Slope(:), smallVelocity(:), ManningsN(:)
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 er_tSlope = er_Temp(next_er_temparray)
 next_er_temparray = utility_advance_temp_array (next_er_temparray,er_n_temp)
 
 er_tManningsN = er_Temp(next_er_temparray)
 next_er_temparray = utility_advance_temp_array (next_er_temparray,er_n_temp)
    
 er_tSmallVelocity = er_Temp(next_er_temparray)
 next_er_temparray = utility_advance_temp_array (next_er_temparray,er_n_temp)

 fUp => elemI(:,ei_Mface_u)
 fDn => elemI(:,ei_Mface_d)
 
 Length           => elemR(:, er_Length)
 Velocity         => elemR(:, er_Velocity_new)
 Flowrate         => elemR(:, er_Flowrate)
 HydRadius        => elemR(:, er_HydRadius)
 SmallVolumeRatio => elemR(:, er_SmallVolumeRatio)
 Area             => elemR(:, er_Area)
 Slope            => elemR(:, er_tSlope)
 smallVelocity    => elemR(:, er_tSmallVelocity)
 ManningsN        => elemR(:, er_tManningsN)
 
 call smallvolume_ManningsN &
    (elemR, elemI, elemYN, er_tManningsN, er_tSmallVelocity, &
     eYN_IsSmallVolume, ei_roughness_type, er_Roughness)

 where ( (elemYN(:,eYN_IsSmallVolume)) .and. &
         (elemI(:,ei_elem_type) == elemType) .and. &
         (SmallVolumeRatio < oneR))
    Slope = ( faceR(fUp, fr_Eta_d) - faceR(fDn, fr_Eta_u) ) / Length
    
    !%  a small velocity based on Chezy-Manning and the free-surface slope
    smallVelocity = small_chezy_velocity (ManningsN,HydRadius,Slope)  
    
    !%  blend with the actual velocity, depending on the small volume ratio                  
    smallVelocity = velocity_blend (SmallVolumeRatio, Velocity, smallVelocity)
    
    !%  take the smaller value, with the sign of the C-M velocity            
    Velocity = sign(min(abs(smallVelocity), abs(Velocity)),smallVelocity)  
    Flowrate = Velocity * Area
 endwhere

 Slope          = nullvalueR
 ManningsN      = nullvalueR
 smallVelocity  = nullvalueR
 nullify(Slope, smallVelocity, ManningsN)
 next_er_temparray = next_er_temparray -3
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine velocity_blend_with_mask
!
!========================================================================== 
!==========================================================================
!
 subroutine smallvolume_ManningsN &
    (elemR, elemI, elemYN, er_tManningsN, er_tSmallVelocity, &
     eYN_IsSmallVolume, ei_roughness_type, er_Roughness)
 
 character(64) :: subroutine_name = 'smallvolume_ManningsN'
!
! Stores a temporary ManningsN and zero velocity for small volumes.
! We use a temporary storage so that we can set Manning's N to the larger value 
! of the actual Manning's N for the channel and a setting value (typically 
! larger)for small depths. 
!
 real,      target,     intent(in out)  :: elemR(:,:)  
 integer,               intent(in)      :: elemI(:,:)
 logical,               intent(in)      :: elemYN(:,:)
 
 integer,   intent(in) :: er_tManningsN, er_tSmallVelocity
 integer,   intent(in) :: eYN_IsSmallVolume, ei_roughness_type, er_Roughness
 
 real,  pointer :: ManningsN(:), smallVelocity(:)
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 ManningsN     => elemR(:,er_tManningsN)
 smallVelocity => elemR(:,er_tSmallVelocity)
 
!%  set the temporary ManningsN and zero velocity for small volumes
 where ( elemYN(:,eYN_IsSmallVolume) )
    ManningsN = setting%SmallVolume%ManningsN
    smallVelocity = zeroR
 endwhere
 
!%  reset the temporary ManningsN whereever the actual value is larger
 where ( (elemYN(:,eYN_IsSmallVolume))              .and. &
         (elemI(:,ei_roughness_type) == eManningsN) .and. &
         (elemR(:,er_Roughness) > ManningsN) )
    ManningsN = elemR(:,er_Roughness)
 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine smallvolume_ManningsN
!
!========================================================================== 
!========================================================================== 
!
 subroutine element_timescale &
    (elem2R, elem2I, elem2YN, elemMR, elemMI, elemMYN, bcdataDn, bcdataUp, &
     e2r_Velocity_new)
!
! timescale of channel and junction elements  
! Computed using new velocity for channel elements
! Assumes that new velocity has been stored in channel branches for junctions
!  
 character(64) :: subroutine_name = 'element_timescale'
 
 real,              intent(in out)  :: elem2R(:,:),  elemMR(:,:)
 integer,           intent(in)      :: elem2I(:,:),  elemMI(:,:)
 logical,           intent(in out)  :: elem2YN(:,:), elemMYN(:,:)
 type(bcType),      intent(in)      :: bcdataDn(:), bcdataUp(:)
 integer,           intent(in)      :: e2r_Velocity_new 

 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 call timescale_value_channel (elem2R, elem2I, elem2YN, e2r_Velocity_new)
    
! note that timescales are stored on branches only for a junction and the data
! in branch flowrate and velocities have already been updated in place, so
! no need to use a eMr..new index
 call timescale_value_junction (elemMR, elemMI, elemMYN)

 call timescale_value_HonlyElement (elem2R, elem2I, elem2YN)

 call timescale_value_QonlyElement (elem2R, elem2I, elem2YN)
 
 call bc_timescale_value (elem2R, bcdataDn)

 call bc_timescale_value (elem2R, bcdataUp)

 ! print*,'-------------------------------------------------'
 ! print*, elem2R(:, e2r_Timescale_Q_u), 'e2r_Timescale_Q_u'
 ! print*,'-------------------------------------------------'
 ! print*, elem2R(:, e2r_Timescale_Q_d), 'e2r_Timescale_Q_d'
 ! print*, '----------------------------------'
 ! print*, eMr_TimescaleUp, 'elemMR tscale up' 
 ! print*, '----------------------------------'
 ! print*, eMr_TimescaleDn, 'elemMR tscale dn'
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine element_timescale
!
!========================================================================== 
!==========================================================================
 subroutine timescale_value_HonlyElement &
    (elem2R, elem2I, elem2YN)

 character(64) :: subroutine_name = 'timescale_value_HonlyElement'
    
 real,      target, intent(in out)  :: elem2R(:,:)
 integer,           intent(in)      :: elem2I(:,:)
 logical,   target, intent(in out)  :: elem2YN(:,:)
 
 integer    ::  indx(2), maskindx1
 
 real,      pointer :: tscale_Q_up(:), tscale_Q_dn(:)
 real,      pointer :: tscale_H_up(:), tscale_H_dn(:), wavespeed(:)
 real,      pointer :: tscale_G_up(:), tscale_G_dn(:), velocity(:)
 real,      pointer :: length(:) 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

 tscale_Q_up => elem2R(:,e2r_Timescale_Q_u)
 tscale_Q_dn => elem2R(:,e2r_Timescale_Q_d)
 tscale_H_up => elem2R(:,e2r_Timescale_H_u)
 tscale_H_dn => elem2R(:,e2r_Timescale_H_d)
 tscale_G_up => elem2R(:,e2r_Timescale_G_u)
 tscale_G_dn => elem2R(:,e2r_Timescale_G_d)

!For Honly meta elements, the timescale for H is minimum. Timescale for Q and G is maximum
!TODO:Storage can be multi faced. So fix for multi faces. 

 where ( (elem2I(:,e2i_meta_elem_type) == eHonly) .and. &
         (elem2I(:,e2i_elem_type) == eStorage) )
    tscale_Q_up = setting%Limiter%Timescale%Maximum
    tscale_Q_dn = setting%Limiter%Timescale%Maximum
    tscale_H_up = setting%Limiter%Timescale%Minimum
    tscale_H_dn = setting%Limiter%Timescale%Minimum
    tscale_G_up = setting%Limiter%Timescale%Maximum
    tscale_G_dn = setting%Limiter%Timescale%Maximum
 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine timescale_value_HonlyElement
!
!==========================================================================
!==========================================================================
!
 subroutine timescale_value_QonlyElement &
    (elem2R, elem2I, elem2YN)

 character(64) :: subroutine_name = 'timescale_value_QonlyElement'
    
 real,      target, intent(in out)  :: elem2R(:,:)
 integer,           intent(in)      :: elem2I(:,:)
 logical,   target, intent(in out)  :: elem2YN(:,:)
 
 integer    ::  indx(2), maskindx1
 
 real,      pointer :: tscale_Q_up(:), tscale_Q_dn(:)
 real,      pointer :: tscale_H_up(:), tscale_H_dn(:), wavespeed(:)
 real,      pointer :: tscale_G_up(:), tscale_G_dn(:), velocity(:)
 real,      pointer :: length(:) 
 logical,   pointer :: maskarray1(:)
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

 tscale_Q_up => elem2R(:,e2r_Timescale_Q_u)
 tscale_Q_dn => elem2R(:,e2r_Timescale_Q_d)
 tscale_H_up => elem2R(:,e2r_Timescale_H_u)
 tscale_H_dn => elem2R(:,e2r_Timescale_H_d)
 tscale_G_up => elem2R(:,e2r_Timescale_G_u)
 tscale_G_dn => elem2R(:,e2r_Timescale_G_d)

!For Qonly meta elements, the timescale for Q is minimum. Timescale for H and G is maximum

 where ( (elem2I(:,e2i_meta_elem_type) == eQonly) )

    tscale_Q_up = setting%Limiter%Timescale%Minimum
    tscale_Q_dn = setting%Limiter%Timescale%Minimum

    tscale_H_up = setting%Limiter%Timescale%Maximum
    tscale_H_dn = setting%Limiter%Timescale%Maximum

    tscale_G_up = setting%Limiter%Timescale%Maximum
    tscale_G_dn = setting%Limiter%Timescale%Maximum


 endwhere

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine timescale_value_QonlyElement
!
!==========================================================================
!==========================================================================
!
 subroutine timescale_value_channel &
    (elem2R, elem2I, elem2YN, e2r_Velocity_new)
 
 character(64) :: subroutine_name = 'timescale_value_channel'
    
 integer,           intent(in)      :: e2r_Velocity_new   
 
 real,      target, intent(in out)  :: elem2R(:,:)
 integer,           intent(in)      :: elem2I(:,:)
 logical,   target, intent(in out)  :: elem2YN(:,:)
 
 integer    ::  indx(2), maskindx1, maskindx2
 
 real,      pointer :: wavespeed(:), velocity(:)
 real,      pointer :: tscale_Q_up(:), tscale_Q_dn(:)
 real,      pointer :: tscale_H_up(:), tscale_H_dn(:)
 real,      pointer :: tscale_G_up(:), tscale_G_dn(:)
 real,      pointer :: length(:) 
 logical,   pointer :: maskarray1(:), maskarray2(:)
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 

 maskindx1 = e2YN_Temp(next_e2YN_temparray) 
 maskarray1 => elem2YN(:,maskindx1)
 next_e2YN_temparray = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)

 maskindx2 = e2YN_Temp(next_e2YN_temparray) 
 maskarray2 => elem2YN(:,maskindx2)
 next_e2YN_temparray = utility_advance_temp_array (next_e2YN_temparray,e2YN_n_temp)
 
 wavespeed => elem2R(:,e2r_Temp(next_e2r_temparray))
 next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)
 
 wavespeed = zeroR
 
 tscale_Q_up => elem2R(:,e2r_Timescale_Q_u)
 tscale_Q_dn => elem2R(:,e2r_Timescale_Q_d)
 tscale_H_up => elem2R(:,e2r_Timescale_H_u)
 tscale_H_dn => elem2R(:,e2r_Timescale_H_d)
 tscale_G_up => elem2R(:,e2r_Timescale_G_u)
 tscale_G_dn => elem2R(:,e2r_Timescale_G_d)
 
 velocity  => elem2R(:,e2r_Velocity_new)
 length    => elem2R(:,e2r_Length)
 
!%  compute timescale 

 maskarray1 = ( (elem2I(:,e2i_elem_type) == eChannel) )

where (maskarray1)
    wavespeed = sqrt( grav * elem2R(:,e2r_HydDepth ))
    tscale_Q_up = - onehalfR * length / (velocity - wavespeed)
    tscale_Q_dn = + onehalfR * length / (velocity + wavespeed)
    tscale_G_up = tscale_Q_up
    tscale_G_dn = tscale_G_dn
    tscale_H_up = tscale_Q_up
    tscale_H_dn = tscale_Q_dn

 endwhere
 
! e2r_Timescale_G_u = e2r_Timescale_Q_u
! e2r_Timescale_G_d = e2r_Timescale_Q_d

 !TODO We have to add the limiter for each type of element base on their physic
 !TODO we have to add timescale calculation for e2r_Timescale_H_d and e2r_Timescale_H_u
 !from Hodges and Liu (2019)

!%  limiter for large, negative, and small values

 indx(1) = e2r_Timescale_Q_u
 indx(2) = e2r_Timescale_Q_d


 call timescale_limiter &
    (elem2R, elem2I, elem2YN, indx, maskindx1, maskindx2)
 
 wavespeed = nullvalueR
 maskarray1 = nullvalueL
 maskarray2 = nullvalueL
 nullify(wavespeed, maskarray1, maskarray2)
 next_e2r_temparray = next_e2r_temparray-1
 next_e2YN_temparray=next_e2YN_temparray-2
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine timescale_value_channel
!
!========================================================================== 
!==========================================================================
!
 subroutine timescale_value_junction &
    (elemMR, elemMI, elemMYN)
 
 character(64) :: subroutine_name = 'timescale_value_junction'
 
 real,  target,     intent(in out)  :: elemMR(:,:)
 integer,           intent(in)      :: elemMI(:,:)
 logical,           intent(in out)  :: elemMYN(:,:)
 
 real,  pointer ::  wavespeed(:), length(:), velocity(:), tscale(:)
    
 integer :: eMr_waveindx
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 eMr_waveindx = eMr_Temp(next_eMr_temparray)
 wavespeed => elemMR(:,eMr_waveindx)
 next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)
 
 wavespeed = zeroR

!%  compute wave speed
 where (elemMI(:,eMi_elem_type) == eJunctionChannel) 
    wavespeed = sqrt( grav * elemMR(:,eMr_HydDepth ))
 endwhere

!% compute timescale for junction in each direction.
 call timescale_junction_one_direction &
    (elemMR, elemMI, upstream_face_per_elemM, eMi_nfaces_u, &
     eMr_LengthUp, eMr_VelocityUp, eMr_TimescaleUp, eMr_waveindx, .true.)
 
 call timescale_junction_one_direction &
    (elemMR, elemMI, dnstream_face_per_elemM, eMi_nfaces_d, &
     eMr_LengthDn, eMr_VelocityDn, eMr_TimescaleDn, eMr_waveindx, .false.)

!%  apply limiters for negative, large, and small values 
 call timescale_limit_junction_one_direction &
    (elemMR, elemMI, elemMYN, upstream_face_per_elemM, eMI_nfaces_u, eMr_TimescaleUp)

 call timescale_limit_junction_one_direction &
    (elemMR, elemMI, elemMYN, dnstream_face_per_elemM, eMI_nfaces_d, eMr_TimescaleDn)
 
 wavespeed = nullvalueR
 nullify(wavespeed)
 next_eMr_temparray = next_eMr_temparray-1

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine timescale_value_junction
!
!========================================================================== 
!==========================================================================
!
 subroutine timescale_junction_one_direction &
    (elemMR, elemMI, dir_face_per_elemM, eMi_nfaces_dir, &
     eMr_LengthDir, eMr_VelocityDir, eMr_TimescaleDir, eMr_waveindx, isUp)
     
 character(64) :: subroutine_name = 'timescale_junction_one_direction'
 
 real,      target, intent(in out)  :: elemMR(:,:)
 integer,           intent(in)      :: elemMI(:,:)
 integer,           intent(in)      :: eMr_LengthDir(:), eMr_VelocityDir(:), eMr_TimescaleDir(:)
 integer,           intent(in)      :: dir_face_per_elemM, eMi_nfaces_dir, eMr_waveindx
 logical,           intent(in)      :: isUp
  
 real,  pointer :: length(:), velocity(:), tscale(:), wavespeed(:) 
 integer :: mm 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 wavespeed => elemMR(:,eMr_waveindx)
 
 do mm=1,dir_face_per_elemM
    length   => elemMR(:,eMr_LengthDir(mm))
    velocity => elemMR(:,eMr_VelocityDir(mm))
    tscale   => elemMR(:,eMr_TimescaleDir(mm))
    if (isUp) then
        where ( (elemMI(:,eMi_elem_type) == eJunctionChannel) .and. &
                (elemMI(:,eMi_nfaces_dir)  >= mm) )
            tscale = - length / (velocity - wavespeed)    
        endwhere
    else
        where ( (elemMI(:,eMi_elem_type) == eJunctionChannel) .and. &
                (elemMI(:,eMi_nfaces_dir)  >= mm) )
            tscale = + length / (velocity + wavespeed)    
        endwhere
    endif
 end do
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine timescale_junction_one_direction
!
!========================================================================== 
!==========================================================================
!
 subroutine timescale_limit_junction_one_direction &
    (elemMR, elemMI, elemMYN, dir_face_per_elemM, eMi_nfaces_dir, eMr_TimescaleDir)
 
 character(64) :: subroutine_name = 'timescale_limit_junction_one_direction'
 
 real,              intent(in out)  :: elemMR(:,:)
 integer,           intent(in)      :: elemMI(:,:)
 logical,   target, intent(in out)  :: elemMYN(:,:)
 integer,           intent(in)      :: dir_face_per_elemM, eMi_nfaces_dir
 integer,           intent(in)      :: eMr_TimescaleDir(:)
 
 
 logical,   pointer    :: maskarray1(:), maskarray2(:)
 integer               :: mm, maskindx1, maskindx2
 integer, dimension(1) :: indx
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 maskindx1  = eMYN_Temp(next_eMYN_temparray)
 maskarray1 => elemMYN(:,maskindx1)
 next_eMYN_temparray = utility_advance_temp_array(next_eMYN_temparray,eMYN_n_temp)

 maskindx2  = eMYN_Temp(next_eMYN_temparray)
 maskarray2 => elemMYN(:,maskindx2)
 next_eMYN_temparray = utility_advance_temp_array(next_eMYN_temparray,eMYN_n_temp)
 
 maskarray1 = nullvalueL
 maskarray2 = nullvalueL

 do mm=1,dir_face_per_elemM
    indx(1) = eMr_TimescaleDir(mm)
    maskarray1 = ((elemMI(:,eMi_nfaces_dir) >= mm) .and. &
                  (elemMI(:,eMi_elem_type) == eJunctionChannel) )
    call timescale_limiter &
        (elemMR, elemMI, elemMYN, indx, maskindx1, maskindx2)
 end do
 
 maskarray1 = nullvalueL
 maskarray2 = nullvalueL
 nullify(maskarray1, maskarray2)
 next_eMYN_temparray = next_eMYN_temparray-2
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine timescale_limit_junction_one_direction
!
!========================================================================== 
!==========================================================================
!
 subroutine timescale_limiter &
    (elemR, elemI, elemYN, indx, maskindx1, maskindx2)
!
! limits the timescales to prevent negatives, small values, or large values
! 
 character(64) :: subroutine_name = 'timescale_limiter'
 
 real,      target, intent(in out)  :: elemR(:,:)
 integer,           intent(in)      :: elemI(:,:)
 logical,   target, intent(in out)  :: elemYN(:,:)
 integer,           intent(in)      :: indx(:)
 integer,           intent(in)      :: maskindx1, maskindx2
 
 real,      pointer :: tscale(:)
 logical,   pointer :: maskarray1(:), maskarray2(:)
 integer            :: mm
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 maskarray1 => elemYN(:,maskindx1)
 maskarray2 => elemYN(:,maskindx2)
 
 do mm=1,size(indx)
 
    tscale => elemR(:,indx(mm))    
 
    !%  ensure negative time scales are stored as large values
    maskarray2 = (maskarray1 .and. ( tscale < 0.0 ))
    call apply_limiter_with_mask &
        (tscale, maskarray2, setting%Limiter%Timescale%Maximum)
        
    !%  ensure small time scales are not below the minimum
    maskarray2 = (maskarray1 .and. (tscale < setting%Limiter%Timescale%Minimum))
    call apply_limiter_with_mask &
        ( tscale, maskarray2, setting%Limiter%Timescale%Minimum)
         
    !%  ensure large time scales are not above the maximum
    maskarray2 = (maskarray1 .and. (tscale > setting%Limiter%Timescale%Maximum))
    call apply_limiter_with_mask &
        (tscale,  maskarray2, setting%Limiter%Timescale%Maximum)
 enddo

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine timescale_limiter
!
!========================================================================== 
!========================================================================== 
!
 pure subroutine apply_limiter_with_mask (inoutarray, maskarray, limitvalue) 
!
! limits the inoutarray to the limitvalue where maskarray is true
! 
 real,      intent(in out)  :: inoutarray(:)
 logical,   intent(in)      :: maskarray(:)
 real,      intent(in)      :: limitvalue
!-------------------------------------------------------------------------- 

 where (maskarray)
    inoutarray = limitvalue
 endwhere
 
 end subroutine apply_limiter_with_mask
!
!========================================================================== 
!========================================================================== 
!
 pure elemental function small_chezy_velocity (ManningsN,HydRadius,Slope)
!
! provides a chezy-manning velocity used for small volumes
!    
    real                :: small_chezy_velocity   
    real,   intent(in)  :: ManningsN, HydRadius, Slope
 
    small_chezy_velocity = sign(  (oneR / ManningsN)       &
                                * (HydRadius**(twothirdR))  &
                                * sqrt(abs(Slope)), Slope  )  
                              
 end function small_chezy_velocity
!
!========================================================================== 
!========================================================================== 
!
 pure elemental function velocity_blend (SmallVolumeRatio, velocity, smallVelocity)
!
! Blends two velocity solutions based on small volume ratio
!  
    real                :: velocity_blend
    real,   intent(in)  :: SmallVolumeRatio, velocity, smallVelocity
    
! HACK - need to ensure that 0 < SmallVolumeRatio < 1 
! Cannot do it here without removing the pure elemental nature.   

    velocity_blend = SmallVolumeRatio * velocity &
                    + (oneR - SmallVolumeRatio) * smallVelocity


 end function velocity_blend
!
!========================================================================== 
! END OF MODULE element_dynamics
!========================================================================== 
 end module element_dynamics
