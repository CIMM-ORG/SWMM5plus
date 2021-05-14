program main

   use globals
   use array_index
   use initialization
   use setting_definition, only: setting
   use interface
   use BIPquick
   use partitioning

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
   call initialize_linknode_arrays()


   ! --- Testing Partitioning Module
   call partitioning_algorithm_check()

   if ( setting%Partitioning%UseDefault .eqv. .true. ) then
        call default_partitioning
   else if ( setting%Partitioning%UseBIPquick .eqv. .true. ) then
        call BIPquick_YJunction_Hardcode()
   end if

   do i = 1, size(P_nodeI, 1)
        print*, P_nodeI(i, :) 
   end do
   do i = 1, size(P_linkI, 1)
        print*, P_linkI(i, :)
   end do
   stop



   ! --- Finalization

   call finalize_api() ! closes link with shared library

end program main