!
! module face_values
!
! Interpolates values from adjacent elements to faces
!
!==========================================================================
!
 module face_values
! 
! compute values on faces from adjacent elements
!
    use adjustments
    use array_index
    use bc
    use data_keys
    use globals
    use setting_definition
    use utility
    
    implicit none
    
    private
    
    public :: face_update

    integer :: debuglevel = 0
    
 contains
!
!========================================================================== 
!==========================================================================
!
 subroutine face_update &
    (elem2R, elemMR, faceR, faceI, faceYN, &
     bcdataDn, bcdataUp, e2r_Velocity_new, eMr_Velocity_new, thisTime, thisIter)
!
! top-level computation of all face values from adjacent elements
! 
 character(64) :: subroutine_name = 'face_update'
 
 integer,       intent(in out)  :: faceI(:,:)
 real,          intent(in out)  :: faceR(:,:)
 real,          intent(in)      :: elem2R(:,:), elemMR(:,:)
 type(bcType),  intent(in out)  :: bcdataDn(:), bcdataUp(:)
 real,          intent(in)      :: thisTime
    
 integer,       intent(in)      :: e2r_Velocity_new, eMr_Velocity_new, thisIter
 
 logical,   intent(in out)  :: faceYN(:,:)
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 call face_interp_for_elem2 (elem2R, faceR, faceI, faceYN, bcdataDn, bcdataUp)
 
 if (N_elemM > 0) then
     print *, 'likely logical errors in the face interp in ',subroutine_name
     stop
 endif

 call bc_applied_onface (faceR, faceI, elem2R, bcdataDn, bcdataUp, thisTime)
 
 !print *, faceR(:,fr_flowrate)
 !print *, faceR(:,fr_Velocity_u)
 !print *, faceR(:,fr_Velocity_d)
 
 if (N_elemM > 0) then
    call face_interp_for_junctionchannel_down(elem2R, elemMR, faceR, faceI, faceYN)
 
    call face_interp_for_junctionchannel_up  (elem2R, elemMR, faceR, faceI, faceYN)
 endif
 
 call face_hydraulic_jump (elem2R, elemMR, faceR, faceI, e2r_Velocity_new, eMr_Velocity_new)
 
 call face_surface_elevation_interp (elem2R, elemMR, faceR, faceI, faceYN)
 
 if (thisIter == 1) then
    !% at end of first step of RK2, the face flow rate is the BC outflow for the step
    call face_bc_flowrate_update (bcdataDn, bcdataUp, faceR)
 endif
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine face_update
!
!========================================================================== 
! PRIVATE BELOW HERE
!==========================================================================
!
 subroutine face_bc_flowrate_update &
    (bcdataDn, bcdataUp, faceR)
!
! stores the last flowrate used on the boundary
! 
 character(64) :: subroutine_name = 'face_bc_flowrate_update'
 
 type(bcType),  target,     intent(in out)  :: bcdataDn(:), bcdataUp(:)
 real,                      intent(in)      :: faceR(:,:)
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 call face_bc_flowrate_update_onevalue (bcdataDn, faceR)
 call face_bc_flowrate_update_onevalue (bcdataUp, faceR)
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine face_bc_flowrate_update
!
!========================================================================== 
!==========================================================================
!
 subroutine face_bc_flowrate_update_onevalue &
    (bcdata, faceR)
 
 character(64) :: subroutine_name = 'face_bc_flowrate_update_onevalue'
 
 type(bcType),  target,     intent(in out)  :: bcdata(:)
 real,                      intent(in)      :: faceR(:,:)
 
 integer,   pointer :: fID
 integer            :: mm
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 do mm=1,size(bcdata)
    fID => bcdata(mm)%FaceID
    bcdata(mm)%ThisFlowrate = faceR(fID,fr_Flowrate)
    !print *,trim(subroutine_name),mm,bcdata(mm)%ThisFlowrate
 enddo
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine face_bc_flowrate_update_onevalue
!
!========================================================================== 
!==========================================================================
!
 subroutine face_interp_for_elem2 &
    (elem2R, faceR, faceI, faceYN, bcdataDn, bcdataUp)
!
! face interpolation between elements that have only 2 faces
! 
 character(64) :: subroutine_name = 'face_interp_for_elem2'
 
 integer,               intent(in)      :: faceI(:,:)
 real,      target,     intent(in out)  :: faceR(:,:)
 real,      target,     intent(in)      :: elem2R(:,:)
 logical,   target,     intent(in out)  :: faceYN(:,:)
 type(bcType),          intent(in)      :: bcdataDn(:), bcdataUp(:)
     
 logical,   pointer  :: facemask(:)
 real,      pointer  :: valueUp(:), valueDn(:), weightUp(:), weightDn(:)
 real,      pointer  :: inoutarray(:)
 
 integer :: mm
 
 integer,   dimension(3)    :: e2rset, frset

!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 valueUp => faceR(:,fr_Temp(next_fr_temparray))
 next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

 valueDn => faceR(:,fr_Temp(next_fr_temparray))
 next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)
 
 weightUp => faceR(:,fr_Temp(next_fr_temparray))
 next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

 weightDn => faceR(:,fr_Temp(next_fr_temparray))
 next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

 facemask   => faceYN(:,fYN_Temp(next_fYN_temparray))
 next_fYN_temparray = utility_advance_temp_array (next_fYN_temparray,fYN_n_temp) 

 facemask = ( (faceI(:,fi_etype_u) == eChannel) .and. &
              (faceI(:,fi_etype_d) == eChannel) )
                      
!%  use timescale for interpolation for Topwidth, Area, Flowrate
 weightUp = setting%Limiter%Timescale%Maximum
 weightDn = setting%Limiter%Timescale%Maximum
 where (facemask)
    weightUp = elem2R(faceI(:,fi_Melem_u),e2r_Timescale_d) !tscale acting downstream from upstream element
    weightDn = elem2R(faceI(:,fi_Melem_d),e2r_Timescale_u) !tscale acting upstream from downstream element
 endwhere
 
!%  set of face interpolations that use timescale weighting
 frset = (/fr_Topwidth,   fr_Area_d, fr_Flowrate /)
 e2rset= (/e2r_Topwidth, e2r_Area,  e2r_Flowrate /)
 
 do mm=1,size(frset)
    call interp_channel_onetype &
        (faceR, facemask, faceI, elem2R, &
         weightUp, weightDn, valueUp, valueDn, e2rset(mm), frset(mm))
 end do
 
!print *, trim(subroutine_name)
!print *, faceR(1:5,fr_flowrate) * setting%Time%DT
!print *, elem2R(2:4,e2r_Volume)
!
!do mm=2,(size(elem2R,1)-2)
!    print *, mm, faceR(mm,fr_flowrate)* setting%Time%DT,  elem2R(mm+1,e2r_Volume)
!    if (faceR(mm,fr_flowrate)* setting%Time%DT > elem2R(mm+1,e2r_Volume)) then
!        print *, mm
!        print *, 'here'
!        !stop
!    endif
!enddo
! 
!%  store the duplicate areas (later adjusted for jumps)
 where (facemask) 
    faceR(:,fr_Area_u) = faceR(:,fr_Area_d)  
 endwhere
 
!%  store simple extrapolation of the face - adjusted during jump procedure
 where (facemask)
    faceR(:,fr_Eta_u) = elem2R(faceI(:,fi_Melem_u),e2r_Eta)
    faceR(:,fr_Eta_d) = elem2R(faceI(:,fi_Melem_d),e2r_Eta)
 endwhere
 
!%  set velocities and upstream values on faces (without hydraulic jump)  

 call adjust_face_dynamic_limits (faceR, faceI, elem2R, facemask)
 
!print *,'after djust_face_dynamic_limit'
!do mm=2,(size(elem2R,1)-2)
!    print *, mm, faceR(mm,fr_flowrate)* setting%Time%DT,  elem2R(mm+1,e2r_Volume)
!    if (faceR(mm,fr_flowrate)* setting%Time%DT > elem2R(mm+1,e2r_Volume)) then
!        print *, mm
!        print *, 'here'
!        stop
!    endif
!enddo
!
  
!%  Store identical values for fr_XXX_u for the moment 
!%  These are later adjusted for hydraulic jumps
 where (facemask)
    faceR(:,fr_Area_u)     = faceR(:,fr_Area_d)
    faceR(:,fr_Velocity_u) = faceR(:,fr_Velocity_d)
 endwhere
  
 facemask = nullvalueL
 nullify(facemask)
 next_fYN_temparray = next_fYN_temparray - 1
 
 valueUp = nullvalueR
 valueDn = nullvalueR
 weightUp = nullvalueR
 weightDn = nullvalueR
 nullify(valueUp, valueDn, weightUp, weightDn)
 next_fr_temparray = next_fr_temparray - 4

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine face_interp_for_elem2
!
!========================================================================== 
!==========================================================================
!
 subroutine face_interp_for_junctionchannel_down &
    (elem2R, elemMR, faceR, faceI, faceYN)
!
! face interpolation with a junction downstream and channel upstream
! i.e. this is from an upstream junction branch to the channel
! 
 character(64) :: subroutine_name = 'face_interp_for_junctionchannel_down'
 
 integer,               intent(in)      :: faceI(:,:)
 real,      target,     intent(in out)  :: faceR(:,:)
 real,      target,     intent(in)      :: elem2R(:,:), elemMR(:,:)
 logical,   target,     intent(in out)  :: faceYN(:,:)
 
 logical,   pointer  :: facemask(:)
 real,      pointer  :: valueUp(:), valueDn(:), weightUp(:), weightDn(:)
 
 integer :: mm

!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 valueUp => faceR(:,fr_Temp(next_fr_temparray))
 next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

 valueDn => faceR(:,fr_Temp(next_fr_temparray))
 next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)
 
 weightUp => faceR(:,fr_Temp(next_fr_temparray))
 next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

 weightDn => faceR(:,fr_Temp(next_fr_temparray))
 next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

 facemask   => faceYN(:,fYN_Temp(next_fYN_temparray))
 next_fYN_temparray = utility_advance_temp_array (next_fYN_temparray,fYN_n_temp) 

 facemask = ( (faceI(:,fi_etype_u) == eChannel) .and. &
              (faceI(:,fi_etype_d) == eJunctionChannel) )
 
 where (facemask)
    weightUp = elem2R(faceI(:,fi_Melem_u),e2r_Timescale_d) !tscale acting downstream
 endwhere
 
 do mm=1,upstream_face_per_elemM
    where ( (facemask) .and. (faceI(:,fi_branch_d) == mm) )
        weightDn = elemMR(faceI(:,fi_Melem_d),eMr_TimescaleUp(mm)) !tscale acting upstream
    endwhere
 end do
 
!%  use timescale for interpolation for Topwidth, Area, Flowrate 
 call interp_junction_down &
    (faceR, facemask, faceI, elem2R, elemMR, &
     weightUp, weightDn, valueUp, valueDn, &
     e2r_Topwidth, eMr_TopwidthUp, fr_Topwidth)
 
 call interp_junction_down &
    (faceR, facemask, faceI, elem2R, elemMR, &
     weightUp, weightDn, valueUp, valueDn, &
     e2r_Area, eMr_AreaUp, fr_Area_d) 

 call interp_junction_down &
    (faceR, facemask, faceI, elem2R, elemMR, &
     weightUp, weightDn, valueUp, valueDn, &
     e2r_Flowrate, eMr_FlowrateUp, fr_Flowrate) 
 
!%  store identical areas (adjusted elsewhere for hyd jump)
 where (facemask)
    faceR(:,fr_Area_u) = faceR(:,fr_Area_d)
 endwhere
 
!%  Store the extraploated upstream and downstream free surface values
!%  Adjusted elsewhere with hydraulic jump   
 where (facemask) 
    faceR(:,fr_Eta_d) = elemMR(faceI(:,fi_Melem_d),eMr_Eta)
    faceR(:,fr_Eta_u) = elem2R(faceI(:,fi_Melem_u),e2r_Eta) 
 endwhere
  
! set velocities and upstream values on faces (without hydraulic jump)  
 print *, trim(subroutine_name),'error - this will not work for elemMR'
 stop
 call adjust_face_dynamic_limits (faceR, faceI, elemMR, facemask)    
    
 facemask = nullvalueL
 nullify(facemask)
 next_fYN_temparray = next_fYN_temparray - 1
 
 valueUp = nullvalueR
 valueDn = nullvalueR
 weightUp = nullvalueR
 weightDn = nullvalueR
 nullify(valueUp, valueDn, weightUp, weightDn)
 next_fr_temparray = next_fr_temparray - 4
     
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine face_interp_for_junctionchannel_down
!
!========================================================================== 
!==========================================================================
!
 subroutine face_interp_for_junctionchannel_up &
    (elem2R, elemMR, faceR, faceI, faceYN)
!
! face interpolation with a junction upstream and channel downstream
! i.e. this is from a downstream junction branch to the downstream channel
! 
 character(64) :: subroutine_name = 'face_interp_for_junctionchannel_up'
 
 integer,               intent(in)      :: faceI(:,:)
 real,      target,     intent(in out)  :: faceR(:,:)
 real,      target,     intent(in)      :: elem2R(:,:), elemMR(:,:)
 logical,   target,     intent(in out)  :: faceYN(:,:)
 
 logical,   pointer  :: facemask(:)
 real,      pointer  :: valueUp(:), valueDn(:), weightUp(:), weightDn(:)
 
 integer :: mm

!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 valueUp => faceR(:,fr_Temp(next_fr_temparray))
 next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

 valueDn => faceR(:,fr_Temp(next_fr_temparray))
 next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)
 
 weightUp => faceR(:,fr_Temp(next_fr_temparray))
 next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

 weightDn => faceR(:,fr_Temp(next_fr_temparray))
 next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

 facemask   => faceYN(:,fYN_Temp(next_fYN_temparray))
 next_fYN_temparray = utility_advance_temp_array (next_fYN_temparray,fYN_n_temp) 

 facemask = ( (faceI(:,fi_etype_u) == eChannel) .and. &
              (faceI(:,fi_etype_d) == eJunctionChannel) )
 
 where (facemask)
    weightDn = elem2R(faceI(:,fi_Melem_d),e2r_Timescale_u) !tscale acting upstream
 endwhere
 
 do mm=1,dnstream_face_per_elemM
    where ( (facemask) .and. (faceI(:,fi_branch_u) == mm) )
        weightUp = elemMR(faceI(:,fi_Melem_u),eMr_TimescaleDn(mm)) !tscale acting dnstream
    endwhere
 end do
 
!%  use timescale for interpolation for Topwidth, Area, Flowrate 
 call interp_junction_up &
    (faceR, facemask, faceI, elem2R, elemMR, &
     weightUp, weightDn, valueUp, valueDn, &
     e2r_Topwidth, eMr_TopwidthDn, fr_Topwidth)
 
 call interp_junction_up &
    (faceR, facemask, faceI, elem2R, elemMR, &
     weightUp, weightDn, valueUp, valueDn, &
     e2r_Area, eMr_AreaDn, fr_Area_d) 

 call interp_junction_up &
    (faceR, facemask, faceI, elem2R, elemMR, &
     weightUp, weightDn, valueUp, valueDn, &
     e2r_Flowrate, eMr_FlowrateDn, fr_Flowrate) 

!%  store identical areas (adjusted elsewhere for hyd jump)
 where (facemask)
    faceR(:,fr_Area_u) = faceR(:,fr_Area_d)
 endwhere
 
!%  Store the extraploated upstream and downstream free surface values
!%  Adjusted elsewhere with hydraulic jump   
 where (facemask) 
    faceR(:,fr_Eta_d) = elem2R(faceI(:,fi_Melem_d),e2r_Eta) 
    faceR(:,fr_Eta_u) = elemMR(faceI(:,fi_Melem_u),eMr_Eta) 
 endwhere
   
!%  set velocities and upstream values on faces (without hydraulic jump)  
 print *, trim(subroutine_name), 'error - this will not work for elemMR'
 stop 
 call adjust_face_dynamic_limits (faceR, faceI, elemMR, facemask)  
    
 facemask = nullvalueL
 nullify(facemask)
 next_fYN_temparray = next_fYN_temparray - 1
 
 valueUp = nullvalueR
 valueDn = nullvalueR
 weightUp = nullvalueR
 weightDn = nullvalueR
 nullify(valueUp, valueDn, weightUp, weightDn)
 next_fr_temparray = next_fr_temparray - 4
     
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine face_interp_for_junctionchannel_up
!
!========================================================================== 
!==========================================================================
!
 subroutine face_hydraulic_jump &
    (elem2R, elemMR, faceR, faceI, e2r_Velocity_new, eMr_Velocity_new)
!
! handles effects of hydraulic jump
! 
 character(64) :: subroutine_name = 'face_hydraulic_jump'
 
 real,      target,     intent(in out)  :: faceR(:,:)
 real,                  intent(in)      :: elem2R(:,:), elemMR(:,:)
 integer,   target,     intent(in out)  :: faceI(:,:)
 integer,               intent(in)      :: e2r_Velocity_new, eMr_Velocity_new
 
 integer,   pointer :: mapUp(:), mapDn(:), typUp(:), typDn(:), jumptype(:)
 real,      pointer :: froudeUp(:), froudeDn(:), etaUp(:), etaDn(:)
 real,      pointer :: areaUp(:), areaDn(:), velUp(:), velDn(:), flowrate(:)
 
 integer    :: fr_froudeUp, fr_froudeDn

 real :: feps
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 fr_froudeUp =  fr_Temp(next_fr_temparray)
 froudeUp    => faceR(:,fr_froudeUp)
 next_fr_temparray = utility_advance_temp_array (next_fr_temparray,fr_n_temp)

 fr_froudeDn =  fr_Temp(next_fr_temparray)
 froudeDn    => faceR(:,fr_froudeDn)
 next_fr_temparray = utility_advance_temp_array (next_fr_temparray,fr_n_temp)

 feps = setting%Eps%FroudeJump

 mapUp => faceI(:,fi_Melem_u)
 mapDn => faceI(:,fi_Melem_d)
 typUp => faceI(:,fi_etype_u)
 typDn => faceI(:,fi_etype_d)
 
 etaUp  => faceR(:,fr_Eta_u)
 etaDn  => faceR(:,fr_Eta_d)
 areaUp => faceR(:,fr_Area_u)
 areaDn => faceR(:,fr_Area_d)
 velUp  => faceR(:,fr_Velocity_u)
 velDn  => faceR(:,fr_Velocity_d)
 flowrate => faceR(:,fr_Flowrate)
 
 jumptype => faceI(:,fi_jump_type)

!%  compute the upstream and downstream froude number on each face 
!%  note the froude numbers are signed with + being downstream flow.
 where (typUp == fChannel)
    froudeUp = elem2R(mapUp,e2r_Velocity_new) / (sqrt(grav * elem2R(mapUp,e2r_HydDepth)))
 endwhere
 
 where (typDn == fChannel)
    froudeDn = elem2R(mapDn,e2r_Velocity_new) / (sqrt(grav * elem2R(mapDn,e2r_HydDepth)))
 endwhere 
 
!%  HACK - does not consider individual branches - only entire element values 
 where (typUp == fMultiple) 
    froudeUp = elemMR(mapUp,eMr_velocity_new) / (sqrt(grav * elemMR(mapUp,eMr_HydDepth)))
 endwhere
    
 where (typDn == fMultiple) 
    froudeDn = elemMR(mapDn,eMr_velocity_new) / (sqrt(grav * elemMR(mapDn,eMr_HydDepth)))
 endwhere
 

!%  Define whether the jump is upstream or downstream or none
!%  Note the mask for water surface level is required because with a coarse
!%  grid (relative to jump scales) on a steep slope it is possible to have 
!%  an eta transition that is not a step in the correct direction. In such 
!%  a case it is better to just use the time weighting approach (which give
!%  priority to the upstream side.
 where ( (froudeDn < oneR - feps) .and. (froudeUp > oneR + feps) .and. &
         (etaUp < etaDn) ) 
    jumptype = jump_downstream     
 endwhere
 
 where ( (froudeDn < -oneR - feps) .and. (froudeUp > -oneR + feps) .and. &
            (etaDn < etaUp) ) 
    jumptype = jump_upstream     
 endwhere
 
 where ((jumptype /= jump_upstream) .and. (jumptype /= jump_downstream))
    jumptype = jump_none
 endwhere
 
!%  Assign areas on either side of jump by extrapolation 
!%  Note that this depends on the fr_Eta_u and fr_Eta_d being initially 
!%  assigned by extrapolation
 where ( (jumptype /= jump_none) .and. (typUp == fChannel) )
    areaUp = elem2R(mapUp,e2r_Area)
 endwhere
 
 where ( (jumptype /= jump_none) .and. (typDn == fChannel) )
    areaDn = elem2R(mapDn,e2r_Area)
 endwhere
 
 where ( (jumptype /= jump_none) .and. (typUp == fMultiple) )
    areaUp = elemMR(mapUp,eMr_Area)
 endwhere
 
 where ( (jumptype /= jump_none) .and. (typDn == fMultiple) )
    areaDn = elemMR(mapDn,eMr_Area)
 endwhere  
 
!%  assign velocities on either side of jump
 where (jumptype /= jump_none)
    velUp  = flowrate / areaUp
    velDn  = flowrate / areaDn
 endwhere

!%  for no jumps, use a linear length interpolation for free surface (elsewhere)
 
 froudeUp = nullvalueR
 froudeDn = nullvalueR
 nullify(froudeUp, froudeDn)
 next_fr_temparray = next_fr_temparray-2
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine face_hydraulic_jump
!
!========================================================================== 
!==========================================================================
!
 subroutine face_surface_elevation_interp &
    (elem2R, elemMR, faceR, faceI, faceYN)
 
 character(64) :: subroutine_name = 'face_surface_elevation_interp'
 
 real,      target,     intent(in out)  :: faceR(:,:)
 real,                  intent(in)      :: elem2R(:,:), elemMR(:,:)
 integer,               intent(in)      :: faceI(:,:)
 logical,   target,     intent(in out)  :: faceYN(:,:)
 
 real,      pointer :: weightUp(:), weightDn(:), etaUp(:), etaDn(:)
 logical,   pointer :: facemask(:)
  
 integer :: mm 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
!%  temporary storage
 weightUp => faceR(:,fr_Temp(next_fr_temparray))
 next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

 weightDn => faceR(:,fr_Temp(next_fr_temparray))
 next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

 facemask   => faceYN(:,fYN_Temp(next_fYN_temparray))
 next_fYN_temparray = utility_advance_temp_array (next_fYN_temparray,fYN_n_temp) 

 etaUp => faceR(:,fr_Eta_u)
 etaDn => faceR(:,fr_Eta_d)
  
!%  use distance (length) for interpolation for free surface 
 where (faceI(:,fi_etype_d) == fChannel)
    weightDn = onehalfR * elem2R(faceI(:,fi_Melem_d),e2r_Length) 
 endwhere
 
 where (faceI(:,fi_etype_u) == fChannel)
     weightUp = onehalfR * elem2R(faceI(:,fi_Melem_u),e2r_Length) 
 endwhere

!%  HACK BELOW 
 if (N_elemM > 0) then
    print *, 'error: the following will not work in ',subroutine_name
    stop
 endif
 do mm=1,upstream_face_per_elemM
    where ( (faceI(:,fi_etype_d) == fMultiple) .and. (faceI(:,fi_branch_d) == mm) )
        weightDn = elemMR(faceI(:,fi_Melem_d),eMr_LengthUp(mm)) !legnth of upstream branch
    endwhere
 end do

 do mm=1,dnstream_face_per_elemM
    where ( (faceI(:,fi_etype_u) == fMultiple) .and. (faceI(:,fi_branch_u) == mm) )
        weightUp = elemMR(faceI(:,fi_Melem_u),eMr_LengthDn(mm)) !legnth of dnstream branch
    endwhere
 end do
!%  END HACK

!%  set the mask for channel and mulitple elements without a hyd jump
 facemask = ( ((faceI(:,fi_etype_d) == fChannel) .or. (faceI(:,fi_etype_d) == fMultiple)) &
             & .and. &
              ((faceI(:,fi_etype_u) == fChannel) .or. (faceI(:,fi_etype_u) == fMultiple)) &
             & .and. &
              (faceI(:,fi_jump_type) == jump_none) )
 
 where (facemask)
    etaDn = (weightUp * etaDn + weightDn * etaUp) / ( weightUp + weightDn )
 endwhere
 etaUp = etaDn

 weightUp = nullvalueR
 weightDn = nullvalueR
 facemask = nullvalueL
 nullify(weightUp, weightDn, facemask)
 next_fr_temparray = next_fr_temparray-2
 next_fYN_temparray = next_fYN_temparray-1
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine face_surface_elevation_interp
!
!========================================================================== 
!==========================================================================
!
 subroutine interp_channel_onetype &
    (faceR, facemask, faceI, elem2R, &
     weightUp, weightDn, valueUp, valueDn, e2r_ThisType, fr_ThisType)
!
! preforms weighted face interpolation for one type of channel
! 
 character(64) :: subroutine_name = 'interp_channel_onetype'
 
 real,      target,     intent(in out)  :: faceR(:,:)
 real,                  intent(in out)  :: valueUp(:),  valueDn(:)
 real,                  intent(in)      :: weightUp(:), weightDn(:)
 real,                  intent(in)      :: elem2R(:,:)
 integer,               intent(in)      :: faceI(:,:)
 
 integer,               intent(in)      :: e2r_ThisType, fr_ThisType
 logical,               intent(in)      :: facemask(:)
   
 real,  pointer :: inoutarray(:)  
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 where (facemask)
    valueUp  = elem2R(faceI(:,fi_Melem_u),e2r_ThisType)
    valueDn  = elem2R(faceI(:,fi_Melem_d),e2r_ThisType)
 endwhere 
    
 inoutarray => faceR(:,fr_ThisType)
 call linear_interpolation &
    (inoutarray, facemask, weightUp, weightDn, valueUp, valueDn)
    
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine interp_channel_onetype
!
!========================================================================== 
!==========================================================================
!
 subroutine interp_junction_down &
    (faceR, facemask, faceI, elem2R, elemMR, &
     weightUp, weightDn, valueUp, valueDn, &
     e2r_ThisType, eMr_ThisTypeUp, fr_ThisType)
!
! interp one value to a face when there is a junction downstream and channel upstream
! Uses branch values in elemMR array.
! 
 
 character(64) :: subroutine_name = 'interp_junction_down'
 
 real,      target,     intent(in out)  :: faceR(:,:)
 real,                  intent(in out)  :: valueUp(:),  valueDn(:)
 real,                  intent(in)      :: weightUp(:), weightDn(:)
 real,                  intent(in)      :: elem2R(:,:), elemMR(:,:)
 integer,               intent(in)      :: faceI(:,:)
 
 integer,               intent(in)      :: e2r_ThisType, fr_ThisType
 integer,               intent(in)      :: eMr_ThisTypeUp(:)
 logical,               intent(in)      :: facemask(:)
   
 real,  pointer :: inoutarray(:)  
    
 integer :: mm   
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 where (facemask)
    valueUp = elem2R(faceI(:,fi_Melem_u),e2r_ThisType) 
 endwhere
 
 do mm=1,upstream_face_per_elemM
    where ( (facemask) .and. (faceI(:,fi_branch_d) == mm) )
        valueDn = elemMR(faceI(:,fi_Melem_d),eMr_ThisTypeUp(mm)) 
    endwhere
 end do
 
 inoutarray => faceR(:,fr_ThisType)
 call linear_interpolation &
    (inoutarray, facemask, weightUp, weightDn, valueUp, valueDn)

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine interp_junction_down
!
!========================================================================== 
!==========================================================================
!
 subroutine interp_junction_up &
    (faceR, facemask, faceI, elem2R, elemMR, &
     weightUp, weightDn, valueUp, valueDn, &
     e2r_ThisType, eMr_ThisTypeDn, fr_ThisType)
!
! interp to a face when there is a junction upstream and channel downstream
! Uses branch values in elemMR array.
! 
 character(64) :: subroutine_name = 'interp_junction_up'
 
 real,      target,     intent(in out)  :: faceR(:,:)
 real,                  intent(in out)  :: valueUp(:),  valueDn(:)
 real,                  intent(in)      :: weightUp(:), weightDn(:)
 real,                  intent(in)      :: elem2R(:,:), elemMR(:,:)
 
 integer,               intent(in)      :: faceI(:,:)
 
 integer,               intent(in)      :: e2r_ThisType, fr_ThisType
 integer,               intent(in)      :: eMr_ThisTypeDn(:)
 
 logical,               intent(in)      :: facemask(:)
   
 real,  pointer :: inoutarray(:)  
    
 integer :: mm
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 where (facemask)
    valueDn = elem2R(faceI(:,fi_Melem_d),e2r_ThisType) 
 endwhere
 
 do mm=1,dnstream_face_per_elemM
    where ( (facemask) .and. (faceI(:,fi_branch_u) == mm) )
        valueUp = elemMR(faceI(:,fi_Melem_u),eMr_ThisTypeDn(mm)) 
    endwhere
 end do
 
 inoutarray => faceR(:,fr_ThisType)
 call linear_interpolation &
    (inoutarray, facemask, weightUp, weightDn, valueUp, valueDn)

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine interp_junction_up
!
!========================================================================== 
!==========================================================================
!
 pure subroutine linear_interpolation &
    (inoutarray, facemask, &
     upstreamWeight, downstreamWeight, upstreamValue, downstreamValue )
!
! linear weighting for interpolation where the upstreamWeight is applied to the
! downstream value - e.g., if length is the weighting, then the upstream value 
! being further away gives a larger emphasis to the downstream value 
! 
 real,      intent(in out)  :: inoutarray(:)
 logical,   intent(in)      :: facemask(:)
 real,      intent(in)      :: upstreamWeight(:), downstreamWeight(:)
 real,      intent(in)      :: upstreamValue(:),  downstreamValue(:)
  
!-------------------------------------------------------------------------- 
!    
 where (facemask)
    inoutarray = (  upstreamWeight   * downstreamValue   &
                  + downstreamWeight * upstreamValue   ) &
                 /( upstreamWeight + downstreamWeight )
 endwhere

 end subroutine linear_interpolation
!
!========================================================================== 
! END OF MODULE face_values
!==========================================================================
 end module face_values