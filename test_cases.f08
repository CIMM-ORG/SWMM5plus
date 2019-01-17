! module test_cases
!
! Calling routines for custom test cases, e.g. calls the case_simple_channel
! functions to setup a single channel reach.
!
!==========================================================================
!
 module test_cases
! 
    use array_index
    use bc
    use case_simple_channel
    use data_keys
    use globals
    use setting_definition
    use utility

    
    implicit none
    
    private
    
    public :: test_case_initiation

    integer :: debuglevel = 0
    
 contains
!
!========================================================================== 
!==========================================================================
!
 subroutine test_case_initiation &
    (linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName, &
     bcdataDn, bcdataUp)
 
 character(64) :: subroutine_name = 'test_case_initiation'
 
 integer,      dimension(:,:), allocatable, intent(out) :: linkI 
 integer,      dimension(:,:), allocatable, intent(out) :: nodeI 
 real,         dimension(:,:), allocatable, intent(out) :: linkR 
 real,         dimension(:,:), allocatable, intent(out) :: nodeR 
 logical,      dimension(:,:), allocatable, intent(out) :: linkYN
 logical,      dimension(:,:), allocatable, intent(out) :: nodeYN
 type(string), dimension(:),   allocatable, intent(out) :: linkName 
 type(string), dimension(:),   allocatable, intent(out) :: nodeName
 type(bcType), dimension(:),   allocatable, intent(out) :: bcdataUp, bcdataDn
 
 real :: Froude, CFL, breadth, ManningsN, total_length, lowerZ
 real :: subdivide_length, area, velocity, flowrate, upperZ
 real :: depth_dnstream, depth_upstream
 
 integer :: first_step, last_step, display_interval, idepth_type
 
 real :: climit, cvel
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 select case (setting%TestCase%TestName)
    
    !% Write a new case statement for each unique test case
    case ('simple_channel_001')
    
        ! step controls
        display_interval = 10000
        first_step = 1
        last_step  =  500000 ! note 1000 is good enough to show blow up or not, 10000 is smooth
    
        ! set up flow and time step for differen subcases
        ! tests that ran:  Fr = 0.25, 0.5
        Froude       = 1.5   ! determines flowrate and slope to get Froude
        CFL          = 0.025  ! determines dt from subdivide_length
    
        ! keep these physics fixed
        breadth         = 3.0
        depth_upstream  = 0.25
        depth_dnstream  = 25.0
        idepth_type     = 3  !1 = uniform, 2=linear, 3=exponential decay
        ManningsN       = 0.03
        total_length    = 10000.0   
        lowerZ          = 1.0
        subdivide_length = 250.0            

        call froude_driven_setup &
            (area, flowrate, velocity, upperZ, &
             Froude, depth_upstream, breadth, ManningsN, &
             lowerZ, total_length)

        call this_setting_for_time_and_steps &
            (CFL, velocity, depth_upstream, subdivide_length, first_step, last_step, & 
             display_interval)   
        
        call case_simple_channel_initialize &
            (total_length, breadth, subdivide_length, lowerZ, upperZ, &
             flowrate, depth_upstream, depth_dnstream, ManningsN,     &
             lManningsN, idepth_type,                                  &
             linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName,     &
             bcdataDn, bcdataUp)    
        
        if (.not. setting%Debugout%SuppressAllFiles) then
            call write_testcase_setup_file &
                (Froude, CFL, flowrate, velocity, depth_upstream,   &
                 depth_dnstream, breadth, area, total_length, subdivide_length, &
                 lowerZ, upperZ, ManningsN)
        endif         
        
    case default
        print *, setting%TestCase%TestName
        print *, 'error: no valid test case of ',&
            trim(setting%TestCase%TestName),' in ',subroutine_name
        stop
 end select
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine test_case_initiation
!
!========================================================================== 
!
! PRIVATE BELOW HERE
!
!==========================================================================
!
 subroutine this_setting_for_time_and_steps &
    (CFL, velocity, depth, subdivide_length, first_step, last_step, & 
     display_interval)
 
 character(64) :: subroutine_name = 'this_setting_for_time_and_steps'
 
 real,  intent(in) :: CFL, velocity, depth, subdivide_length
 
 integer, intent(in) :: first_step, last_step, display_interval
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 setting%Time%dt = get_dt_from_CFL (CFL, velocity, depth, subdivide_length)

 setting%Step%Current = 1

 setting%Step%First = first_step
 setting%Step%Final = last_step  

 setting%Debugout%DisplayInterval = display_interval
 
 setting%Time%StartTime = 0.0
 setting%Time%EndTime = setting%Time%StartTime  &
    + setting%Time%dt * (setting%Step%Final - setting%Step%First + 1)
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine this_setting_for_time_and_steps
!
!========================================================================== 
!==========================================================================
!
 subroutine froude_driven_setup &
    (area, flowrate, velocity, upperZ, &
     Froude, depth, breadth, ManningsN, &
     lowerZ, total_length)
 
 character(64) :: subroutine_name = 'froude_driven_setup'
  
 real,  intent(out)    :: area, flowrate, velocity, upperZ 
 real,  intent(in)     :: Froude,  breadth, ManningsN, lowerZ, total_length
 real,  intent(in)     :: depth 
 
 real :: perimeter, rh, slope
 
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 area = depth * breadth
 perimeter = 2.0 * depth + breadth
 rh = area / perimeter
 velocity = Froude * sqrt(grav * depth)
 flowrate = area * velocity
 slope = (velocity * ManningsN / (rh**(2.0/3.0)) )**2
 upperZ = lowerZ + slope * total_length
        

! print *, area
! print *, perimeter
! print *, rh
! print *, velocity
! print *, flowrate
! print *, slope
! print *, upperZ, lowerZ
! print *, total_length
! print *, slope*total_length
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine froude_driven_setup
!
!==========================================================================
!==========================================================================
!
 subroutine write_testcase_setup_file &
    (Froude, CFL, flowrate, velocity, depth_upstream, depth_dnstream, breadth,  &
     area, total_length, subdivide_length, lowerZ, upperZ, ManningsN)
 
 character(64) :: subroutine_name = ' write_testcase_setup_file'
 
 real,  intent(in)  :: Froude, CFL, flowrate, velocity,  breadth
 real,  intent(in)  :: area, total_length, subdivide_length, lowerZ, upperZ
 real,  intent(in)  :: ManningsN, depth_upstream, depth_dnstream
    
 integer        :: UnitNumber
 
 character(64)  :: thisFilePath, thisFileStatus, thisFileName
 character(256) :: thisFileWriteName
 
 logical        :: thisFileisOpen = .false.
 
 integer                :: open_status
 character(len=512)     :: emsg
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 open_status = 0
 
 UnitNumber = outputfile_next_unitnumber
 outputfile_next_unitnumber = outputfile_next_unitnumber+1
 
 thisFileName   = trim(setting%TestCase%TestName)
 thisFilePath   = trim(setting%DebugOut%FolderPath) &
                // trim(setting%Debugout%FolderName) // '/'
 thisFileStatus = 'new'
 thisFileIsOpen     = .true.
    
 thisFileWriteName  = trim(thisFilePath) // &
                      trim(thisFileName) // &
                      trim(setting%Time%DateTimeStamp) //&
                      '.txt'
 
! print *, trim(setting%TestCase%TestName)
! print *, trim(setting%DebugOut%FolderPath)
! print *, trim(setting%Debugout%FolderName)
! 
! print *, trim(thisFileName)
! print *, trim(thisFilePath)
! print *, trim(thisFileWriteName)
! stop
! 
 open(unit=UnitNumber, &
      file=trim(thisFileWriteName), &
      status = 'new', &
      access = 'sequential', &
      form   = 'formatted', &
      action = 'write', &
      iostat = open_status) 
     
 emsg = 'file exists: file open failed in '//trim(subroutine_name) &
       // '; filename = '//trim(thisFileWriteName) 
 call utility_check_fileopen (open_status, emsg)  
 
 write(UnitNumber,*) trim(setting%TestCase%TestName)
 write(UnitNumber,*) trim(setting%Time%DateTimeStamp)
 write(UnitNumber,*)
 write(UnitNumber,*) Froude  ,'=Froude'
 write(UnitNumber,*) CFL     ,'=CFL combined'
 write(UnitNumber,*) velocity * setting%Time%Dt / subdivide_length,'=CFL advective'
 write(UnitNumber,*) sqrt(grav * depth_upstream) * setting%Time%DT / subdivide_length,'=CFL barotropic'
 write(UnitNumber,*)
 write(UnitNumber,*) flowrate, '=flowrate'
 write(UnitNumber,*) velocity, '=velocity' 
 write(UnitNumber,*) setting%Time%Dt,' = dt'
 write(UnitNumber,*) 
 write(UnitNumber,*) depth_upstream   ,'=depth upstream'
 write(UnitNumber,*) depth_dnstream   ,'=depth downstream'
 write(UnitNumber,*) breadth ,'=breadth'
 write(UnitNumber,*) area    ,'=area'
 write(UnitNumber,*) total_length ,'=total_length'
 write(UnitNumber,*) subdivide_length ,'=subdivide_length'
 write(UnitNumber,*) area * subdivide_length,'=element_volume'
 write(UnitNumber,*) lowerZ,'=lowerZ'
 write(UnitNumber,*) upperZ,'=upperZ'
 write(UnitNumber,*) (upperZ - lowerZ )/ total_length,'=slope'
 write(UnitNumber,*) 
 write(UnitNumber,*) ManningsN,'=ManningsN'
 write(UnitNumber,*)
 write(UnitNumber,*) setting%Step%First,'=first step'
 write(UnitNumber,*) setting%Step%Final,'=last step'
 write(UnitNumber,*) setting%Time%StartTime,'=start time'
 write(UnitNumber,*) setting%Time%EndTime,'=end time'
 write(UnitNumber,*)

 close(UnitNumber)
 outputfile_next_unitnumber = outputfile_next_unitnumber-1       
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine  write_testcase_setup_file
!
!========================================================================== 
!==========================================================================
!
 function get_dt_from_CFL &
    (CFL, velocity, depth, element_length) &
    result (dt)
    
 character(64) :: subroutine_name = 'get_dt_from_CFL'   
 
 real,  intent(in) :: CFL, velocity, depth, element_length   
 real :: dt
 
!-------------------------------------------------------------------------- 
    
 dt = CFL * onehalfR * element_length / (velocity + sqrt(grav * depth))
 
! print *, subroutine_name
! print *, 'ad hoc cut of time step *******************************'
! dt = 0.25 * dt
 
 if (floor(dt) > 19) then
    dt = real(floor(dt))
 endif
 
 end function get_dt_from_CFL
!
!========================================================================== 
! END OF MODULE test_cases
!==========================================================================
 end module test_cases