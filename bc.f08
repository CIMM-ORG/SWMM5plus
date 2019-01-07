! module bc
!
! Setup and apply boundary conditions. All the BC functions and subroutines
! should be located in this module.
!
! Note that BC should only be applied on faces connected to an elem2 element.
! That is, a BC cannot be associated with a junction branch.
!
!==========================================================================
!
 module bc
!
! boundary condition definitions and enforcement
! 
    use array_index
    use data_keys
    use globals
    use setting_definition
    use utility
    
    implicit none
    
    private
    
    public :: bcType
    public :: bc_allocate
    public :: bc_checks
    public :: bc_applied_onface
    public :: bc_applied_onelement
    public :: bc_timescale_value
    public :: bc_nullify_ghost_elem
    
    
    type bcType
        integer :: Idx
        integer :: NodeID
        integer :: FaceID
        integer :: ElemGhostID
        integer :: ElemInsideID
        integer :: Updn      ! bc_updn_...  (0 = upstream,  1 = downstream)
        integer :: Category  ! bc_category_... (0 = elevation, 1 = inflowrate)
        real    :: FroudeInflowMaximum = 1.5 ! max value of Fr at inflow 
        real, dimension(:), allocatable :: TimeArray
        real, dimension(:), allocatable :: ValueArray
        real    :: ThisValue
        real    :: ThisTime
        real    :: ThisFlowrate
    end type bcType

    character(len=99), private  :: emsg
    
    integer :: debuglevel = 0
    
 contains
!
!========================================================================== 
!========================================================================== 
!
 subroutine bc_allocate &
    (bcdataDn, bcdataUp, ndnstreamBC, nupstreamBC, ntimepoint)
!
! allocate storage for boundary conditions.
! HACK - possibly move this to allocation_storage module
! 
 character(64) :: subroutine_name = 'bc_allocate'
 
 type(bcType), dimension(:),   allocatable,    intent(out) :: bcdataDn, bcdataUp
 
 integer,   intent(in)  :: ndnstreamBC, nupstreamBC, ntimepoint
 
 integer    :: ii

 integer            :: allocation_status
 character(len=99)  :: emsg
   
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
!% the Upstream and Downstream bc structure
 allocate( bcdataUp(nupstreamBC), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)

 allocate( bcdataDn(ndnstreamBC), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)

!% the downstream arrays - HACK default downstream is elevation
 bcdataDn%updn       = bc_updn_downstream
 bcdataDn%category   = bc_category_elevation
 bcdataDn%faceID     = nullvalueI
 bcdataDn%thisValue  = nullvalueR
 bcdataDn%thisTime   = nullvalueR
 
 do ii=1,ndnstreamBC
    allocate( bcdataDn(ii)%TimeArray(ntimepoint), stat=allocation_status, errmsg=emsg)
    call utility_check_allocation (allocation_status, emsg)
    
    allocate( bcdataDn(ii)%ValueArray(ntimepoint), stat=allocation_status, errmsg=emsg)
    call utility_check_allocation (allocation_status, emsg)
    
    bcdataDn(ii)%idx = ii
 end do

!% the upstream arrays = HACK default upstream is flowrate
 bcdataUp%updn       = bc_updn_upstream
 bcdataUp%category   = bc_category_inflowrate
 bcdataUp%faceID     = nullvalueI
 bcdataUp%thisValue  = nullvalueR
 bcdataUp%thisTime   = nullvalueR 
 
 do ii=1,nupstreamBC
    allocate( bcdataUp(ii)%TimeArray(ntimepoint), stat=allocation_status, errmsg=emsg)
    call utility_check_allocation (allocation_status, emsg)

    allocate( bcdataUp(ii)%ValueArray(ntimepoint), stat=allocation_status, errmsg=emsg)
    call utility_check_allocation (allocation_status, emsg)  
    
    bcdataUp(ii)%idx = ii   
 end do

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine bc_allocate
!
!========================================================================== 
!==========================================================================
!
 subroutine bc_applied_onface &
    (faceR, faceI, elem2R, bcdataDn, bcdataUp, thisTime)
!
! Apply boundary conditions to face arrays.
!
! This requires the elem2R data because an elevation BC is enforced on the
! ghost element so the face enforcement is an interpolation across the
! boundary.
!
 character(64) :: subroutine_name = 'bc_applied_onface'
 
 real,          intent(in out)  :: faceR(:,:)
 real,          intent(in)      :: elem2R(:,:)
 integer,       intent(in)      :: faceI(:,:)
 type(bcType),  intent(in out)  :: bcdataDn(:), bcdataUp(:)
 
 real,  intent(in)  :: thisTime 
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 call bc_onface(faceR, elem2R, bcdataDn, thisTime)
 call bc_onface(faceR, elem2R, bcdataUp, thisTime)
 
 call bc_face_othervalues (faceR, faceI, elem2R, bcdataDn) 
 call bc_face_othervalues (faceR, faceI, elem2R, bcdataUp) 
  
         
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 
 end subroutine bc_applied_onface
!
!========================================================================== 
!==========================================================================
!
 subroutine bc_applied_onelement &
    (elem2R, bcdataDn, bcdataUp, thisTime, thiscategory, e2r_VelocityColumn)
!
! Apply boundary condition on elements (not faces) for bc:Category = "thiscategory" 
! 
! The e2r_VelocityColumn is only used when thiscategory = bc_category_inflowrate
! in which case it may be either the e2r_Velocity or the e2r_Velocity_new column. 
! 
 character(64) :: subroutine_name = 'bc_applied_onelement'
 
 real,          intent(in out)  :: elem2R(:,:) 
 type(bcType),  intent(in out)  :: bcdataDn(:), bcdataUp(:)
 
 real,      intent(in)  :: thisTime 
 integer,   intent(in)  :: thiscategory
 integer,   intent(in)  :: e2r_VelocityColumn
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 call bc_onelement (elem2R, bcdataDn, thisTime, thiscategory, e2r_VelocityColumn)
 call bc_onelement (elem2R, bcdataUp, thisTime, thiscategory, e2r_VelocityColumn)
         
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine bc_applied_onelement
!
!========================================================================== 
!==========================================================================
!
 subroutine bc_checks &
    (bcdataUp, bcdataDn, elem2I, faceI, nodeI)
!
! checking of BC setup
!
 character(64) :: subroutine_name = 'bc_checks'
 
 integer,           intent(in)      :: elem2I(:,:), faceI(:,:)
 integer, target,   intent(in)      :: nodeI(:,:)
 
 type(bcType), intent(in out) :: bcdataUp(:), bcdataDn(:)
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
    
!% error checking for assignment of BC nodes to faces
 call bc_node_assignment_error_check (nodeI, faceI, nBCdn) 
 call bc_node_assignment_error_check (nodeI, faceI, nBCup)  
 
!% Assign bc face ID that corresponds to nodes
 call bc_assign_faceID (bcdataDn, faceI)
 call bc_assign_faceID (bcdataUp, faceI) 
 
!% assign initial values 
 call bc_updatevalue (bcdataDn, setting%Time%StartTime)
 call bc_updatevalue (bcdataUp, setting%Time%StartTime)
 
!% check to see that the estimated simulation time is within available BC data 
 call bc_adequate_coverage (bcdataUp, bcdataDn)
!  
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine bc_checks
!
!==========================================================================
!==========================================================================
!
 subroutine bc_timescale_value &
    (elem2R, bcdata)
!
! set timescales on ghost elements outside boundary to max value
! 
 character(64) :: subroutine_name = 'bc_timescale_value'
 
 real,                  intent(in out)  :: elem2R(:,:)
 type(bcType),  target, intent(in)      :: bcdata(:)
    
 integer :: ii
 
 integer,   pointer :: eID
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 do ii=1,size(bcdata)
    eID => bcdata(ii)%ElemGhostID
    elem2R(eID,e2r_Timescale_u) = setting%Limiter%Timescale%Maximum
    elem2R(eID,e2r_Timescale_d) = setting%Limiter%Timescale%Maximum
 end do
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine bc_timescale_value
!
!========================================================================== 
!==========================================================================
!
 subroutine bc_nullify_ghost_elem &
    (elem2R, bcdata)
!
! This sets all the ghost elements outside of boundary faces to nullvalues
! In general this should not be needed but is useful in debugging
! 
 character(64) :: subroutine_name = 'bc_nullify_ghost_elem'
 
 real,                      intent(in out) :: elem2R(:,:)
 type(bcType),  target,     intent(in)     :: bcdata(:)
 
 integer,   pointer :: eID
 integer            :: ii
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 do ii=1,size(bcdata)
    eID => bcdata(ii)%ElemGhostID
    elem2R(eID,:) = nullvalueR
 end do
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine bc_nullify_ghost_elem
!
!========================================================================== 
! PRIVATE BELOW HERE 
!==========================================================================
!
 subroutine bc_adequate_coverage &
    (bcdataUp, bcdataDn)
!
!check to see if the BC cover the entire simulation period
!
 character(64) :: subroutine_name = 'bc_adequate_coverage'

 type(bcType),  intent(in)  :: bcdataUp(:), bcdataDn(:) 

!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 call bc_adequate_coverage_onedir (bcdataUp, nBCup)    
 call bc_adequate_coverage_onedir (bcdataDn, nBCdn)  
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine bc_adequate_coverage
!
!========================================================================== 
!==========================================================================
!
 subroutine bc_onface &
    (faceR, elem2R, bcdata, thisTime)
!
! Apply boundary conditions on one face (not element)
! 
 character(64) :: subroutine_name = 'bc_onface'
 
 real,                      intent(in out)  :: faceR(:,:)
 real,                      intent(in)      :: elem2R(:,:)
 type(bcType),  target,     intent(in out)  :: bcdata(:)
 
 real,  intent(in)  :: thisTime 
 
 real,      pointer :: thisval
 integer,   pointer :: thisloc, thiscat, thisghost, thisinside
 integer :: ii
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 call bc_updatevalue (bcdata, thisTime)
 
 do ii=1,size(bcdata)
    thisloc => bcdata(ii)%faceID
    thisghost    => bcdata(ii)%ElemGhostID
    thisinside   => bcdata(ii)%ElemInsideID
    thiscat => bcdata(ii)%category
    thisval => bcdata(ii)%thisValue
    
    if (thiscat == bc_category_elevation) then
        !%  linear interpolation using ghost and interior cells
        faceR(thisloc,fr_Eta_u) = onehalfR *( elem2R(thisghost,e2r_Eta) + elem2R(thisinside,e2r_Eta))
        faceR(thisloc,fr_Eta_d) = faceR(thisloc,fr_Eta_u)
        
    elseif (thiscat == bc_category_inflowrate) then
        !%  direct enforcement of flowrate
        faceR(thisloc,fr_Flowrate) = thisval
    else
        print *, 'error: unexpected value for bcdata%category of ',thiscat,' in ',subroutine_name
        stop
    endif
 end do
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine bc_onface
!
!========================================================================== 
!==========================================================================
!
 subroutine bc_face_othervalues &
    (faceR, faceI, elem2R, bcdata)
!
! enforces the non-bc values on a bc face, e.g. flowrate on a face that
! has an elevation bc
! 
 character(64) :: subroutine_name = 'bc_face_othervalues'
 
 real,    target,   intent(in out)  :: faceR(:,:)
 real,    target,   intent(in)      :: elem2R(:,:)
 integer, target,   intent(in)      :: faceI(:,:)
 type(bcType),  target, intent(in)  :: bcdata(:)
 !integer,               intent(in)  :: fi_Melem_inside
 
 real                   :: thisFroudeNumber, thisDepth, sideslope
 real,      pointer     :: froudeMax, Qrate, Depth, Top, Dinc, Area
 integer, dimension(2)  :: e2rset, frset 
 integer,   pointer     :: fID, eID
 integer    :: ii
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 e2rset(1:2) = (/e2r_Topwidth, e2r_Area /)
 frset(1:2)  = (/fr_Topwidth,  fr_Area_d /)

!%  Extrapolation from the inside to the face at the boundary
!%  then overwrite with the actual bc
 do ii=1,size(bcdata)
    fID => bcdata(ii)%FaceID
    eID => bcdata(ii)%ElemInsideID
    froudeMax => bcdata(ii)%FroudeInflowMaximum
    !%  simple extrapolation of geometry from inside element
    !% for area and topwidth
    faceR(fID,frset)     = elem2R(eID,e2rset) 
    faceR(fID,fr_Area_u) =  faceR(fID,fr_Area_d)

    !%  elevation at non-elevation BC by extrapolating depth
    if (bcdata(ii)%Category /= bc_category_elevation) then
        !% INFLOW BC - need eta and velocity on face
        
        !faceR(fID,fr_Eta_d) = elem2R(eID,e2r_eta)
        !faceR(fID,fr_Eta_u) = faceR(fID,fr_Eta_d)  
        faceR(fID,fr_Eta_d) = faceR(fID,fr_Zbottom) + elem2R(eID,e2r_HydDepth)
        faceR(fID,fr_Eta_u) = faceR(fID,fr_Eta_d)
        if (faceR(fID,fr_Area_d) > zeroR) then
            faceR(fID,fr_Velocity_d) = faceR(fID,fr_flowrate) / faceR(fID,fr_Area_d)
            thisFroudeNumber = faceR(fID,fr_Velocity_d) / (sqrt(grav * elem2R(eID,e2r_HydDepth)) ) 
        else
            faceR(fID,fr_Velocity_d) = setting%Limiter%Velocity%Maximum
            thisFroudeNumber = twoR * setting%Eps%InflowDepthIncreaseFroudeLimit
        endif
        faceR(fID,fr_Velocity_u) = faceR(fID,fr_Velocity_d)
        
        !print *, trim(subroutine_name)
        !print *, fID ,faceR(fID,fr_flowrate), faceR(fID,fr_Velocity_d), faceR(fID,fr_Velocity_u)
        !print *, 'Fr=',thisFroudeNumber, 'Vel=',faceR(fID,fr_Velocity_d), 'D=',elem2R(eID,e2r_HydDepth)
        

        !%  Froude number limiter on inflow - changes area, velocity,topwidth, and eta
        !%  Note that depth is not stored on face, so we use interior depth
        if (thisFroudeNumber > froudeMax) then
            !% if max froude is exceeded, we need to reduce the inflow area and
            !% depth. Note we do not actually store geometry on the inflow, so
            !% we do this by approximation. We will assume the geometry devolves
            !% into a simple triangular area of Topwidth T and depth D, so that the
            !% area is given by 
            !% A = TD/2. 
            !% Then Fr^2 = 4Q^2 / (g T^2 D^3)
            !% The side slope is H = 2D / T so that
            !% T = 2D/H
            !% Thus, Fr^2 = 4Q^2 H^2 / (g 4D^2 D^3) = (QH)^2 / (g D^5)
            !% or D^5 = (QH)^2 / g Fr^2
            !% We assume that the high froude number conditions have the same
            !% effective side slop of H =2D(old) / T(old)
            Qrate       => faceR(fID,fr_Flowrate)
            Depth       => elem2R(eID,e2r_HydDepth) !depth not defined on face
            Top         => faceR(fID,fr_Topwidth)
            Area        => faceR(fID,fr_Area_d)
            sideslope   =  twoR * Depth / Top
            Dinc        => setting%Eps%InflowDepthIncreaseFroudeLimit
            
            thisDepth = ((((Qrate * sideslope) / froudeMax)**2) / grav )**(1.0/5.0)
            thisDepth = (oneR + Dinc)*thisDepth
            
            !print *, 'D=', thisDepth, elem2R(eID,e2r_HydDepth)
            !print *, 'T=',twoR * thisDepth / sideslope, faceR(fID,fr_Topwidth)
            
            faceR(fID,fr_Topwidth)  = twoR * thisDepth / sideslope
            
            !print *, 'A=',Top  * thisdepth / twoR, faceR(fID,fr_Area_d)
            
            faceR(fID,fr_Area_d)    = Top  * thisdepth / twoR
            
            faceR(fID,fr_Area_u)    = faceR(fID,fr_Area_d)
            
            !print *, 'E=',faceR(fID,fr_Zbottom) + thisDepth, faceR(fID,fr_Eta_d)
            
            faceR(fID,fr_Eta_d)     = faceR(fID,fr_Zbottom) + thisDepth
            faceR(fID,fr_Eta_u)     = faceR(fID,fr_Eta_d)
                                  
        endif
    end if
    
    !%  flowrate at non-flowrate BC from inside
    if (bcdata(ii)%Category /= bc_category_inflowrate) then
        !% ELEVATION BC
        
        faceR(fID,fr_Flowrate) = elem2R(eID,e2r_Flowrate)
        !%  reset the velocity
        if ( (setting%ZeroValue%UseZeroValues) .and.  &
             (faceR(fID,fr_Area_d) > setting%Zerovalue%Area)) then
            faceR(fID,fr_Velocity_d) = faceR(fID,fr_Flowrate) / faceR(fID,fr_Area_d)
            faceR(fID,fr_Velocity_u) = faceR(fID,fr_Velocity_d)
        else
            faceR(fID,fr_Velocity_d) = zeroR
            faceR(fID,fr_Velocity_u) = zeroR
        end if
    end if    
 end do
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine bc_face_othervalues
!
!========================================================================== 
!==========================================================================
!
 subroutine bc_onelement &
    (elem2R, bcdata, thisTime, thiscategory, e2r_VelocityColumn)
!
! Apply boundary conditions (either up or down) to an element.
!
! For elevation, this only affects a ghost cell (for use in interpolation to
! face in bc_onface). For flowrate this affects both ghost and interior cells
! for both flowrate and velocity.
!
! e2r_VelocityColumn only used when thiscategory = bc_category_inflowrate
! in which case it may be either e2r_Velocity or e2r_Velocity_new.

 character(64) :: subroutine_name = 'bc_onelement'
 
 real,                      intent(in out)  :: elem2R(:,:)
 type(bcType),  target,     intent(in out)  :: bcdata(:)
 real,                      intent(in)      :: thisTime
 integer,                   intent(in)      :: thiscategory, e2r_VelocityColumn
    
 real,      pointer :: thisval
 integer,   pointer :: thisloc, thiscat, thisghost
 integer :: ii  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 

 call bc_updatevalue (bcdata, thisTime)
 
 do ii=1,size(bcdata)
    thisghost => bcdata(ii)%ElemGhostID
    thisloc   => bcdata(ii)%ElemInsideID
    thiscat   => bcdata(ii)%category
    thisval   => bcdata(ii)%thisValue
    if (thiscat == thiscategory) then
    
        if (thiscat == bc_category_elevation) then
            !%  store elevation BC on ghost element
            elem2R(thisghost,e2r_Eta) = thisval
            
        elseif (thiscat == bc_category_inflowrate) then
            ! store the flowrate and compute the corresponding velocity on interior
            !elem2R(thisloc,e2r_Flowrate)       = thisval
            !elem2R(thisloc,e2r_VelocityColumn) = thisval / elem2R(thisloc,e2r_Area)
            
            ! store the same value on the ghost 
            elem2R(thisghost,e2r_Flowrate)       = thisval
            elem2R(thisghost,e2r_VelocityColumn) = thisval / elem2R(thisloc,e2r_Area)

        else
            print *, 'error: unexpected value for bcdata%category of ',thiscat,' in ',subroutine_name
            stop
        endif
    endif
 end do
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine bc_onelement
!
!========================================================================== 
!==========================================================================
!
 subroutine bc_updatevalue &
    (bcdata, thisTime)
!
! Interpolate the bcdata value for thisTime from the data array stored 
! in TimeArray:ValueArray
! 
 character(64) :: subroutine_name = 'bc_updatevalue'
 
 type(bcType), intent(in out) :: bcdata(:)
 
 real,  intent(in)  :: thisTime
 
 integer    :: ii
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 do ii=1,size(bcdata)
    bcdata%thisValue = utility_linear_interpolate_within_indexlist &
                        (thisTime, bcdata(ii)%TimeArray, bcdata(ii)%ValueArray )
    bcdata%thisTime  = thisTime
 end do
    
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 
 end subroutine bc_updatevalue
!
!========================================================================== 
!==========================================================================
!
 subroutine bc_node_assignment_error_check &
    (nodeI, faceI, nBCdir)
!
! check for node assignment errors in creating the link/node system
! 
 character(64) :: subroutine_name = 'bc_node_assignment_error_check'
 
 integer,   target,     intent(in)  :: nodeI(:,:)
 integer,               intent(in)  :: faceI(:,:)
 
 integer,               intent(in)  :: nBCdir
 
 character(len=8) :: cdir
 
 integer            :: ii
 integer,   pointer :: nidx
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 if (nBCdir == nBCdn) then
    cdir = 'dnstream'
 elseif (nBCdir == nBCup) then
    cdir = 'upstream'
 else
    print *, 'error: unknown value for nBCdir in ',subroutine_name
    stop
 endif
 
 do ii=1,size(nodeI,1)
    ! this node ID
    nidx => nodeI(ii,ni_idx)
    ! check to see if this node is a BC node
    if (nodeI(ii,ni_node_type) == nBCdir) then 
        ! count the number of faces that use this node as a BC
        ! there should be only 1.
        select case (count(faceI(:,fi_node_ID) == nidx))
            case (0)
                print *, 'error: BC ',cdir,' at node ',ii,' with ni_idx ',nidx, &
                 'does not have corresponding fi_node_ID in faceI array in ',subroutine_name
                stop
            case (1)
                ! no action - all OK
            case default
                print *, 'error: BC ',cdir,' at node ',ii,' with ni_idx ',nidx, &
                 'has more than one corresponding fi_node_ID in faceI array in ',subroutine_name
                stop
        end select
    endif    
 enddo
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine bc_node_assignment_error_check
!
!========================================================================== 
!==========================================================================
!
 subroutine bc_assign_faceID &
    (bcdata, faceI)
!
! Identify the faceID, ghostID and insideID for a bcdata set where the
! NodeID is already assigned
! 
 character(64) :: subroutine_name = 'bc_assign_faceID'
 
 type(bcType),  target, intent(in out) :: bcdata(:)
 
 integer,               intent(in)     :: faceI(:,:)
 
 integer,   pointer :: nID, fID, insideID, ghostID, updn
 integer :: ii 
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 

 do ii=1,size(bcdata)
    !%  known
    nID      => bcdata(ii)%NodeID  ! already assigned
    updn     => bcdata(ii)%Updn   
    
    !%  unknown
    fID      => bcdata(ii)%FaceID
    insideID => bcdata(ii)%ElemInsideID
    ghostID  => bcdata(ii)%ElemGhostID

    
    !%  location where the face nodeID is the same as the bc nodeID
    fID  = minloc(abs(faceI(:,fi_node_ID) - nID),1)
    
    !%  elem that is outside the face
    if     (updn == bc_updn_downstream) then
        ghostID  = faceI(fID,fi_Melem_d)
        insideID = faceI(fid,fi_Melem_u)
    elseif (updn == bc_updn_upstream)   then
        ghostID  = faceI(fID,fi_Melem_u)
        insideID = faceI(fID,fi_Melem_d)
    else
        print *, 'error: unexpected else in ',subroutine_name
        stop
    end if
 end do
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine bc_assign_faceID
!
!==========================================================================
!==========================================================================
!
 subroutine bc_adequate_coverage_onedir &
    (bcdata, nBCdir)
!
! check to see if the BC cover the entire simulation period
! 
 character(64) :: subroutine_name = 'bc_adequate_coverage_onedir'
 
 type(bcType),              intent(in)  :: bcdata(:)
 
 integer,                   intent(in)  :: nBCdir
 
 real :: timelow, timehigh
 
 integer :: ii
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 do ii=1,size(bcdata)
    timelow = bcdata(ii)%TimeArray(1)
    !print *, timelow
    timehigh= bcdata(ii)%TimeArray(size(bcdata(ii)%TimeArray,1))
    !print *, timehigh
    !print *, setting%Time%StartTime
    !print *, setting%Time%EndTime
    if (timelow > setting%Time%StartTime) then
        if (nBCdir == nBCup) then
            print *, 'error in bcdataUp(',ii,') in ',subroutine_name 
        elseif (nBCdir == nBCdn) then    
            print *, 'error in bcdataDn(',ii,') in ',subroutine_name 
        else
            print *, 'error: unexpected value of nBCdir ',nBCdir,' in ',subroutine_name 
        endif
        print *, 'start time of (',setting%Time%StartTime ,&
                 ' is below lowest available BC data time ',timelow
        stop
    endif
    if (timehigh < setting%Time%EndTime) then
        if (nBCdir == nBCup) then
            print *, 'error in bcdataUp(',ii,') in ',subroutine_name 
        elseif (nBCdir == nBCdn) then  
            print *, 'error in bcdataDn(',ii,') in ',subroutine_name 
        else
            print *, 'error: unexpected value of nBCdir ',nBCdir,' in ',subroutine_name 
        endif
        print *, 'end time of (',setting%Time%EndTime, &
                 ' is greater than highest available BC data time of ',timehigh
        stop
    endif
    
 enddo
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine bc_adequate_coverage_onedir
!
!==========================================================================
! END OF MODULE bc
!==========================================================================
 end module bc