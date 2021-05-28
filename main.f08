program main


   use globals
   use assign_index
   use initialization
   use setting_definition, only: setting
   use interface
   use coarray_partition
   use BIPquick
   use partitioning
   use network_define

   implicit none

   integer :: i
   logical :: arg_param = .false.
   character(len=8) :: param
   character(len=256) :: arg

   ! ---  Define paths

   call getcwd(setting%Paths%project)
   setting%Paths%setting = trim(setting%Paths%project) // '/initialization/settings.json'

   ! --- Read args

   do i = 1, iargc()
       call getarg(i, arg)
       if (.not. arg_param) then
           param = arg
           if (i == 1) then
               if (arg(:1) == '-') then
                   print *, "ERROR: it is necessary to define the path to the .inp file"
                   stop
               endif
               setting%Paths%inp = arg
           elseif ((trim(arg) == "-s") .or. & ! user provides settings file
               (trim(arg) == "-t")) then  ! hard coded test case
               arg_param = .true.
           else
               write(*, *) 'The argument ' // trim(arg) // ' is unsupported'
               stop
           end if
       else
           arg_param = .false.
           if (trim(param) == '-s') then
               setting%Paths%setting = arg
           elseif (trim(param) == '-t') then
               setting%TestCase%UseTestCase = .true.
               setting%TestCase%TestName = trim(arg)
               if (trim(arg) == 'simple_channel') then
               else if (trim(arg) == 'simple_orifice') then
               else if (trim(arg) == 'simple_pipe') then
               else if (trim(arg) == 'simple_weir') then
               else if (trim(arg) == 'swashes') then
               else if (trim(arg) == 'waller_creek') then
               else if (trim(arg) == 'y_channel') then
               else if (trim(arg) == 'y_storage_channel') then
               else
                   write(*, *) 'The test case ' // trim(arg) // ' is unsupported. Please use one of the following:'
                   print *, new_line('')
                   print *, "simple_channel, simple_orifice, simple_pipe"
                   print *, "simple_weir, swashes, waller_creek"
                   print *, "y_channel, y_storage_channel"
                   stop
               end if
           end if
       end if
   end do

   ! --- Load Settings

   call load_settings(setting%Paths%setting)
   call execute_command_line("if [ -d debug ]; then rm -r debug; fi && mkdir debug")

   ! --- Initialization

   call initialize_api()

   sync all
   
   call initialize_linknode_arrays()

   if(this_image()==1) then
     print*,'----------------------------------------------------------------'
     print*, 'upstream links to node 3'
     print*, nodeI(3,ni_Mlink_u1), 'Upstream link 1'
     print*, nodeI(3,ni_Mlink_u2), 'Upstream link 2'
     print*, nodeI(3,ni_Mlink_u3), 'Upstream link 3'
     print*, nodeI(3,ni_N_link_u), 'number of us link'
     print*
     print*, 'downstream links to node 3'
     print*
     print*, nodeI(3,ni_Mlink_d1), 'Downstream link 1'
     print*, nodeI(3,ni_Mlink_d2), 'Downstream link 2'
     print*, nodeI(3,ni_Mlink_d3), 'Downstream link 3'
     print*, nodeI(3,ni_N_link_d), 'number of ds link'
     print*,'----------------------------------------------------------------'
     print*, 'node types to understand how many us and ds boundary conditions'
     print*, nodeI(1,ni_node_type), '<== node 1 Typ'
     print*, nodeI(2,ni_node_type), '<== node 2 Typ'
     print*, nodeI(3,ni_node_type), '<== node 3 Typ'
     print*, nodeI(4,ni_node_type), '<== node 4 Typ'
     print*,'Up bc = 5, Dn bc = 4, nJm = 2'
     stop
    endif


  
   sync all

   call execute_partitioning()

   sync all
   
   call network_initiation()
   
   sync all
   
   ! --- Finalization
   call finalize_api() ! closes link with shared library

   print*, 'End of main'
   

   

end program main
