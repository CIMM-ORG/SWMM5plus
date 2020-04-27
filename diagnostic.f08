! module diagnostic
!
! Utilities that are useful for diagnosing simulation behavior
!
!==========================================================================
!
 module diagnostic
! 
    use array_index
    use bc
    use data_keys
    use globals
    use setting_definition
    use type_definitions
    use utility

    implicit none
    
    public :: diagnostic_CFL
    public :: diagnostic_initialize
    public :: diagnostic_element_volume_conservation
    public :: diagnostic_element_volume_conservation_fluxes
    public :: diagnostic_froude_number
    public :: diagnostic_volume_conservation

    private

    integer :: debuglevel = 0
    
 contains
!
!==========================================================================
!==========================================================================
!
 pure function diagnostic_CFL &
    (elemR, er_Timescale_u, er_Timescale_d) result(cflmax)
!
! Computes thebarotropic + advection CFL from the timescales up and down
!
 real,      intent(in)      :: elemR(:,:)
 integer,   intent(in)      :: er_Timescale_u, er_Timescale_d
    
 real     :: cflmax
 real     :: cflu, cfld
  
!-------------------------------------------------------------------------- 
 !if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 cflu = maxval( setting%Time%Dt / elemR(:,er_Timescale_u), 1, &
               (elemR(:,er_Timescale_u) > setting%Limiter%Timescale%Minimum) )

 cfld = maxval( setting%Time%Dt / elemR(:,er_Timescale_d), 1, &
               (elemR(:,er_Timescale_d) > setting%Limiter%Timescale%Minimum) ) 
               
 cflmax = max(cflu, cfld)              
               
 !if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end function diagnostic_CFL
!
!========================================================================== 
!==========================================================================
!
 subroutine diagnostic_element_volume_conservation_fluxes &
    (elem2R, elem2I, elemMR, elemMI, faceR)
 
 character(64) :: subroutine_name = 'diagnostic_element_volume_conservation_fluxes'
 
 real,              intent(in out)  :: elem2R(:,:), elemMR(:,:)
 real,              intent(in)      :: faceR(:,:)
 integer, target,   intent(in)      :: elem2I(:,:), elemMI(:,:)
 
 integer,   pointer :: fup(:), fdn(:)
 
 integer :: mm
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
!% store the net volume flux through the face of continuity
 fup => elem2I(:,e2i_Mface_u)
 fdn => elem2I(:,e2i_Mface_d)
 where (elem2I(:,e2i_elem_type) == eChannel)
    elem2R(:,e2r_VolumeConservation) = - dt * (faceR(fup,fr_Flowrate) - faceR(fdn,fr_Flowrate) )
 endwhere 

! This needed to be fixed for weir element
 where (elem2I(:,e2i_elem_type) == eWeir)
    elem2R(:,e2r_VolumeConservation) = - dt * (faceR(fup,fr_Flowrate) - faceR(fdn,fr_Flowrate) )
 endwhere 

 ! This needed to be fixed for orifice element
 where (elem2I(:,e2i_elem_type) == eOrifice)
    elem2R(:,e2r_VolumeConservation) = - dt * (faceR(fup,fr_Flowrate) - faceR(fdn,fr_Flowrate) )
 endwhere 
 
 where ( elemMR(:,eMi_elem_type) == eJunctionChannel )
    elemMR(:,eMr_VolumeConservation) = zeroR
 endwhere
 
 do mm=1,upstream_face_per_elemM
    fup => elemMI(:,eMi_MfaceUp(mm))
    where ( (elemMR(:,eMi_elem_type) == eJunctionChannel) .and. &
            (elemMR(:,eMi_nfaces_u) >= mm) )
        elemMR(:,eMr_VolumeConservation) = elemMR(:,eMr_VolumeConservation) - dt * faceR(fup,fr_flowrate)
    endwhere
 enddo
 
 do mm=1,dnstream_face_per_elemM
    fdn => elemMI(:,eMi_MfaceDn(mm))
    where ( (elemMR(:,eMi_elem_type) == eJunctionChannel) .and. &
            (elemMR(:,eMi_nfaces_d) >= mm) )
        elemMR(:,eMr_VolumeConservation) = elemMR(:,eMr_VolumeConservation) + dt * faceR(fdn,fr_flowrate)
    endwhere
 enddo 
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine diagnostic_element_volume_conservation_fluxes
!
!========================================================================== 
!==========================================================================
!
 subroutine diagnostic_element_volume_conservation &
    (elem2R, elem2I, elemMR, elemMI, e2r_Volume_new, eMr_Volume_new)
!
! compute local volume conservation (subtracting net change from flux volume)
! this assumes that fluxes have beens tored in VolumeConservation column
! in a prior call to diagnostic_element_volume_conservation_fluxes  
!
 character(64) :: subroutine_name = 'diagnostic_element_volume_conservation'
 
 real,      intent(in out)  :: elem2R(:,:), elemMR(:,:)
 integer,   intent(in)      :: elem2I(:,:), elemMI(:,:)
 integer,   intent(in)      :: e2r_Volume_new, eMr_Volume_new
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 

 where (elem2I(:,e2i_elem_type) == eChannel)
    elem2R(:,e2r_VolumeConservation) = elem2R(:,e2r_VolumeConservation) &
         + (elem2R(:,e2r_Volume_new) - elem2R(:,e2r_Volume))
 endwhere   

! This needed to be fixed for weir element
 where (elem2I(:,e2i_elem_type) == eWeir)
    elem2R(:,e2r_VolumeConservation) = elem2R(:,e2r_VolumeConservation) &
         + (elem2R(:,e2r_Volume_new) - elem2R(:,e2r_Volume))
 endwhere

 ! This needed to be fixed for orifice element
 where (elem2I(:,e2i_elem_type) == eOrifice)
    elem2R(:,e2r_VolumeConservation) = elem2R(:,e2r_VolumeConservation) &
         + (elem2R(:,e2r_Volume_new) - elem2R(:,e2r_Volume))
 endwhere        
    

 where (elemMI(:,eMi_elem_type) == eJunctionChannel)
    elemMR(:,eMr_VolumeConservation) = elemMR(:,eMr_VolumeConservation) &
         + (elemMR(:,eMr_Volume_new) - elemMR(:,eMr_Volume))
 endwhere   

! print *, trim(subroutine_name)
! print *, elem2R(:,e2r_VolumeConservation) 
! print *, elemMR(:,eMr_VolumeConservation)
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine diagnostic_element_volume_conservation
!
!========================================================================== 
!==========================================================================
!
 subroutine diagnostic_froude_number &
    (elem2R, elem2I, elemMR, elemMI)
 
 character(64) :: subroutine_name = 'diagnostic_froude_number'
 
 real,      intent(in out)  :: elem2R(:,:), elemMR(:,:)
 integer,   intent(in)      :: elem2I(:,:), elemMI(:,:)
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 call diagnostic_froude_number_one &
    (elem2R, elem2I, e2r_FroudeNumber, e2r_Velocity, e2r_HydDepth, e2i_elem_type, eChannel)

!diagnostic_froude_number_one needed to be adapted for weir element
 call diagnostic_froude_number_one &
    (elem2R, elem2I, e2r_FroudeNumber, e2r_Velocity, e2r_HydDepth, e2i_elem_type, eWeir)

!diagnostic_froude_number_one needed to be adapted for Orifice element
 call diagnostic_froude_number_one &
    (elem2R, elem2I, e2r_FroudeNumber, e2r_Velocity, e2r_HydDepth, e2i_elem_type, eOrifice)
    
 call diagnostic_froude_number_one &
    (elemMR, elemMI, eMr_FroudeNumber, eMr_Velocity, eMr_HydDepth, eMi_elem_type, eJunctionChannel)    
    
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine diagnostic_froude_number
!
!========================================================================== 
!==========================================================================
!
 subroutine diagnostic_initialize &
    (diagnostic, elem2R, elem2I, elemMR, elemMI, faceR, &
     bcdataUp, bcdataDn)
 
 character(64) :: subroutine_name = 'diagnostic_initialize'
 
 type(diagnosticType), allocatable,    dimension(:), intent(out)    :: diagnostic
 
 real,                  intent(in out)  :: elem2R(:,:), elemMR(:,:), faceR(:,:)
 integer,   target,     intent(in)      :: elem2I(:,:), elemMI(:,:)
 type(bcType),          intent(in)      :: bcdataUp(:), bcdataDn(:)
 
 integer            :: allocation_status, total_steps
 character(len=99)  :: emsg

  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 

 total_steps = setting%Step%Final - setting%Step%First + 2
 
!%  create the storage of data
!%  note that position 1 is storage of initial conditions.
 allocate( diagnostic(total_steps), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 
 diagnostic%Volume%Step                 = 0
 
!% initialize all real values to zero 
 diagnostic%Volume%Time                 = zeroR
 diagnostic%Volume%Volume               = zeroR
 diagnostic%Volume%VolumeChange         = zeroR
 diagnostic%Volume%NetInflowVolume      = zeroR
 diagnostic%Volume%InflowRate           = zeroR
 diagnostic%Volume%OutflowRate          = zeroR
 diagnostic%Volume%ConservationThisStep = zeroR
 diagnostic%Volume%ConservationTotal    = zeroR
 
!% specific initializeation for volume conservation diagnostic (task = 0) 
 call diagnostic_volume_conservation &
    (diagnostic, elem2R, elem2I, elemMR, elemMI, faceR, &
     bcdataUp, bcdataDn, 1, 0) 

 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine diagnostic_initialize
!
!========================================================================== 
!==========================================================================
!
 subroutine diagnostic_volume_conservation &
    (diagnostic, elem2R, elem2I, elemMR, elemMI, faceR, &
     bcdataUp, bcdataDn, thisStep, diagnosticTask)
!
! handles volume conservation diagnostic
!   diagnosticTask = 0 -> initialize
!   diagnosticTask = 1 -> compute starting volume or transfer end volume to start
!   diagnosticTask = 2 -> compute volume change and BC flows 
!
 character(64) :: subroutine_name = 'diagnostic_volume_conservation'
 
 type(diagnosticType),  intent(in out)  :: diagnostic(:)
 real,                  intent(in out)  :: elem2R(:,:), elemMR(:,:), faceR(:,:)
 integer,   target,     intent(in)      :: elem2I(:,:), elemMI(:,:)
 type(bcType),          intent(in)      :: bcdataUp(:), bcdataDn(:)
 integer,               intent(in)      :: thisStep
 integer,               intent(in)      :: diagnosticTask
    
 integer,   pointer :: etype2(:), etypeM(:)
 real               :: channelVolume, junctionVolume, weirVolume, orificevolume, totalVolume
 real               :: inflowRate, outflowRate
 
 integer :: ii
    
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 etype2 => elem2I(:,e2i_elem_type)
 etypeM => elemMI(:,eMi_elem_type)
 
!  do ii= 1, size(elem2R(:,e2r_Volume),1)
!     print*, "elem2R(:,e2r_Volume) = ii =", ii, " == ", elem2R(ii,e2r_Volume)
!  enddo
 
 channelVolume  = sum(elem2R(:,e2r_Volume),1,etype2 == eChannel)
 weirVolume     = sum(elem2R(:,e2r_Volume),1,etype2 == eWeir)
 orificevolume  = sum(elem2R(:,e2r_Volume),1,etype2 == eOrifice)
 junctionVolume = sum(elemMR(:,eMr_Volume),1,etypeM == eJunctionChannel) 
 
 totalVolume = channelVolume + junctionVolume + weirVolume + orificevolume
  
 select case (diagnosticTask)
    case (0)
        !%  reset the accumulation storage
        diagnostic(thisStep)%Volume%ConservationTotal = zeroR
    case (1)
        !%  store only the volume (required prior to first step)
        diagnostic(thisStep)%Volume%Volume  = totalVolume
        diagnostic(thisStep)%Volume%Step = setting%Step%Current
        diagnostic(thisStep)%Volume%Time = setting%Time%ThisTime
    case (2)
        !% store switch the previous "this" volume to "last"
        !% update this volume
        diagnostic(thisStep)%Volume%Step = setting%Step%Current
        diagnostic(thisStep)%Volume%Time = setting%Time%ThisTime
        diagnostic(thisStep)%Volume%Volume = totalVolume
        diagnostic(thisStep)%Volume%VolumeChange = totalVolume &
            -  diagnostic(thisStep-1)%Volume%Volume
        
        call total_inout_flowrate &
            (inflowRate, outflowRate, faceR, bcdataUp, bcdataDn)
        
        diagnostic(thisStep)%Volume%InflowRate    = inflowRate
        diagnostic(thisStep)%Volume%OutflowRate   = outflowRate
        diagnostic(thisStep)%Volume%NetInflowVolume = &
            (diagnostic(thisStep)%Volume%InflowRate - diagnostic(thisStep)%Volume%OutflowRate) * dt
        
        diagnostic(thisStep)%Volume%ConservationThisStep =       &
                diagnostic(thisStep)%Volume%VolumeChange &
              - diagnostic(thisStep)%Volume%NetInflowVolume  
        
        diagnostic(thisStep)%Volume%ConservationTotal =  &
                diagnostic(thisStep-1)%Volume%ConservationTotal &
              + diagnostic(thisStep  )%Volume%ConservationThisStep
        
    case default
        print *,'error: unknown value for diagnosticTask in ',trim(subroutine_name)
        print *, diagnosticTask
        stop
 end select
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine diagnostic_volume_conservation
!
!========================================================================== 
!
! PRIVATE BELOW HERE
!
!==========================================================================
!
 subroutine diagnostic_froude_number_one &
    (elemR, elemI, er_FroudeNumber, er_Velocity, er_HydDepth, ei_elem_type, thisType)
 
 character(64) :: subroutine_name = 'diagnostic_froude_number_one'
 
 real,      intent(in out)  :: elemR(:,:)
 integer,   intent(in)      :: elemI(:,:)
 integer,   intent(in)      :: er_FroudeNumber, er_Velocity, er_HydDepth
 integer,   intent(in)      :: ei_elem_type, thisType
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 where (elemI(:,ei_elem_type) == thisType)
    elemR(:,er_FroudeNumber) = elemR(:,er_Velocity) / sqrt(grav * elemR(:,er_HydDepth))
 endwhere
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine diagnostic_froude_number_one
!
!========================================================================== !==========================================================================
!
 subroutine total_inout_flowrate &
    (inflowRate, outflowRate, faceR, bcdataUp, bcdataDn)
 
 character(64) :: subroutine_name = 'total_inout_flowrate'
 
 real,          target, intent(in)      :: faceR(:,:) 
 type(bcType),  target, intent(in)      :: bcdataUp(:), bcdataDn(:)
 real,                  intent(out)     :: inflowRate, outflowRate  
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 inflowRate  = zeroR
 outflowRate = zeroR
 
 call inout_flowrate_from_bcdata &
    (inflowRate, outflowRate, faceR, bcdataUp)
    
 call inout_flowrate_from_bcdata &
    (inflowRate, outflowRate, faceR, bcdataDn)   
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine total_inout_flowrate
!
!========================================================================== 
!==========================================================================
!
 subroutine inout_flowrate_from_bcdata &
    (inflowRate, outflowRate, faceR, bcdata)
!
! adds the inflow rate and outflow rate from bcdata to the
! accumulators inflowRate and outflowRate
 
 character(64) :: subroutine_name = 'inout_flowrate_from_bcdata'
 
 real,          target, intent(in)      :: faceR(:,:) 
 type(bcType),  target, intent(in)      :: bcdata(:)
 real,                  intent(in out)  :: inflowRate, outflowRate
 
 integer            :: mm
 integer,   pointer :: fID, bupdn
 real,      pointer :: flowrate
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 

 do mm=1,size(bcdata)
 
    flowrate => bcdata(mm)%ThisFlowrate
    bupdn    => bcdata(mm)%Updn

    if (bupdn == bc_updn_downstream) then
        !%  downstream boundaries are inflows if negative, outflows if positive
        if (flowrate < zeroR) then
            inflowRate  = inflowRate  - flowrate
        else
            outflowRate = outflowRate + flowrate
        endif
    elseif (bupdn == bc_updn_upstream) then
        !%  upstream boundaries are outflows if negative, inflows if positive
        if (flowrate < zeroR) then
            outflowRate = outflowRate - flowrate
        else
            inflowRate  = inflowRate  + flowrate
        endif
    else
        print *, 'error: unknown value for bcdata%Updn in ',trim(subroutine_name)
        print *, mm, bcdata(mm)%Updn
        stop
    endif
 enddo
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine inout_flowrate_from_bcdata
!
!========================================================================== 
! END OF MODULE diagnostic
!==========================================================================
 end module diagnostic
