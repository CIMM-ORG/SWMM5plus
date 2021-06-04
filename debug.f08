! module debug
!
! Provides simple text files output of arrays for debugging.
! This module needs significant restructuring if we want to maintain
! it within the code.
!
!==========================================================================
!
module debug
    !
    ! handles the debug output file opening, writing, and closing
    !
    use array_index
    use bc
    use data_keys
    use globals
    use setting_definition
    use type_definitions
    use utility


    implicit none

    public ::   debug_initialize
    public ::   debug_output
    public ::   debug_finalize

    private

    integer :: debuglevel = 0

contains
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine debug_initialize &
        (debugfile)
        !
        !   opens debug files for output (HACK - needs work to simplify)
        !
        character(64) :: subroutine_name = 'debug_initialize'

        type(debugfileType),  dimension(:),   allocatable, intent(out) :: debugfile

        integer             :: allocation_status
        character(len=512)   :: emsg

        integer :: ndebug
        integer :: ii

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !%  skip this if files are suppressed
        if (setting%Debugout%SuppressAllFiles) return

        !% HACK the following is rather brittle code. The numbers representing the
        !% number of debug files to be opened must be exactly equal to the calls to the
        !% debug_singlefile_open routine. This can be made much more compact and
        !% efficient - but it is not clear that this is a subroutine that will be
        !% needed for the long term code.
        ndebug = 0
        if (setting%Debugout%elem2R) then
            ndebug = ndebug + 23
        endif
        if (setting%Debugout%elem2I) then
            ndebug = ndebug + 10
        endif
        if (setting%Debugout%elem2YN) then
            ndebug = ndebug + 2
        endif

        if (setting%Debugout%elemMR) then
            ndebug = ndebug + 44
        endif
        if (setting%Debugout%elemMI) then
            ndebug = ndebug + 9
        endif
        if (setting%Debugout%elemMYN) then
            ndebug = ndebug + 2
        endif

        if (setting%Debugout%faceR) then
            ndebug = ndebug + 10
        endif
        if (setting%Debugout%faceI) then
            ndebug = ndebug + 1
            ! ndebug = ndebug + 12
        endif
        if (setting%Debugout%faceYN) then
            ndebug = ndebug + 0
        endif

        if (setting%Debugout%nodeR) then
            ndebug = ndebug + 6
        endif
        if (setting%Debugout%nodeI) then
            ndebug = ndebug + 5
        endif
        if (setting%Debugout%nodeYN) then
            ndebug = ndebug + 0
        endif

        if (setting%Debugout%linkR) then
            ndebug = ndebug + 5
        endif
        if (setting%Debugout%linkI) then
            ndebug = ndebug + 12
        endif
        if (setting%Debugout%linkYN) then
            ndebug = ndebug + 0
        endif

        allocate( debugfile(ndebug), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation (allocation_status, emsg)

        !%  HACK the following can be simplified by defining some string arrays
        !%  and an array of column indexes to use in an iterative call.
        ii=1
        if (setting%Debugout%elem2R) then
            call debug_singlefile_open (debugfile(ii), 'elem2R','Volume', e2r_Volume)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','SmallVolume', e2r_SmallVolume)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','SmallVolumeRatio', e2r_SmallVolumeRatio)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','Flowrate', e2r_Flowrate)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','Velocity',e2r_Velocity)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','Timescale_Q_u', e2r_Timescale_Q_u)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','Timescale_Q_d', e2r_Timescale_Q_d)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','Friction', e2r_Friction)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','Eta', e2r_Eta)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','Head', e2r_Head)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','Area', e2r_Area)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','Topwidth', e2r_Topwidth)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','Perimeter', e2r_Perimeter)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','HydDepth', e2r_HydDepth)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','HydRadius', e2r_HydRadius)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','X', e2r_X)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','Length', e2r_Length)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','Zbottom', e2r_Zbottom)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','BreadthScale', e2r_BreadthScale)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','Roughness', e2r_Roughness)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','VolumeConservation', e2r_VolumeConservation)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','FroudeNumber', e2r_FroudeNumber)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2R','FullDepth', e2r_FullDepth)
            ii=ii+1
        endif

        if (setting%Debugout%elem2I) then
            call debug_singlefile_open (debugfile(ii), 'elem2I','idx', e2i_idx)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2I','elem_type', e2i_meta_elem_type)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2I','elem_type', e2i_elem_type)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2I','geometry', e2i_geometry)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2I','roughness_type', e2i_roughness_type)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2I','link_ID', e2i_link_ID)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2I','link_Pos', e2i_link_Pos)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2I','Mface_u', e2i_Mface_u)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2I','Mface_d', e2i_Mface_d)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elem2I','solver', e2i_solver)
            ii=ii+1
        endif

        if (setting%Debugout%elem2YN) then
            print *, 'error: code for debug output of elem2YN not completed in ',subroutine_name
            stop
        endif

        if (setting%Debugout%elemMR) then
            call debug_singlefile_open (debugfile(ii), 'elemMR','Volume', eMr_Volume)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','SmallVolume', eMr_SmallVolume)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','SmallVolumeRatio', eMr_SmallVolumeRatio)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Flowrate', eMr_Flowrate)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Velocity', eMr_Velocity)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Friction', eMr_Friction)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Eta', eMr_Eta)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Head', eMr_Head)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Area', eMr_Area)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Topwidth', eMr_Topwidth)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Perimeter', eMr_Perimeter)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','HydDepth', eMr_HydDepth)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','X', eMr_X)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Length', eMr_Length)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Zbottom', eMr_Zbottom)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','BreadthScale', eMr_BreadthScale)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Roughness', eMr_Roughness)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','VolumeConservation', eMr_VolumeConservation)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','FroudeNumber', eMr_FroudeNumber)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Flowrate_u1', eMr_Flowrate_u1)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Flowrate_u2', eMr_Flowrate_u2)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Flowrate_d1', eMr_Flowrate_d1)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Velocity_u1', eMr_Velocity_u1)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Velocity_u2', eMr_Velocity_u2)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Velocity_d1', eMr_Velocity_d1)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Timescale_u1', eMr_Timescale_u1)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Timescale_u2', eMr_Timescale_u2)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Timescale_d1', eMr_Timescale_d1)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Timescale_d2', eMr_Timescale_d2)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Area_u1', eMr_Area_u1)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Area_u2', eMr_Area_u2)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Area_d1', eMr_Area_d1)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Topwidth_u1', eMr_Topwidth_u1)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Topwidth_u2', eMr_Topwidth_u2)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Topwidth_d1', eMr_Topwidth_d1)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Length_u1', eMr_Length_u1)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Length_u2', eMr_Length_u2)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Length_d1', eMr_Length_d1)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Zbottom_u1', eMr_Zbottom_u1)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Zbottom_u2', eMr_Zbottom_u2)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','Zbottom_d1', eMr_Zbottom_d1)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','BreadthScale_u1', eMr_BreadthScale_u1)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','BreadthScale_u2', eMr_BreadthScale_u2)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'elemMR','BreadthScale_d1', eMr_BreadthScale_d1)
            ii=ii+1
        endif

        if (setting%Debugout%elemMI) then
            print *, 'error: code for debug output of elemMI not completed in ',subroutine_name
            stop
        endif

        if (setting%Debugout%elemMYN) then
            print *, 'error: code for debug otuput of elemMYN not completed in ',subroutine_name
            stop
        endif

        if (setting%Debugout%faceR) then
            call debug_singlefile_open (debugfile(ii), 'faceR','Area_d', fr_Area_d)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'faceR','Area_u', fr_Area_u)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'faceR','Eta_d', fr_Eta_d)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'faceR','Eta_u', fr_Eta_u)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'faceR','Flowrate', fr_Flowrate)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'faceR','Topwidth', fr_Topwidth)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'faceR','Velocity_d', fr_Velocity_d)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'faceR','Velocity_u', fr_Velocity_u)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'faceR','Zbottom', fr_Zbottom)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'faceR','X', fr_X)
            ii=ii+1
        endif

        if (setting%Debugout%faceI) then
            ! call debug_singlefile_open (debugfile(ii), 'faceI','idx', fi_idx)
            ! ii=ii+1
            ! call debug_singlefile_open (debugfile(ii), 'faceI','fi_type', fi_type)
            ! ii=ii+1
            ! call debug_singlefile_open (debugfile(ii), 'faceI','Melem_u', fi_Melem_u)
            ! ii=ii+1
            ! call debug_singlefile_open (debugfile(ii), 'faceI','Melem_d', fi_Melem_d)
            ! ii=ii+1
            ! call debug_singlefile_open (debugfile(ii), 'faceI','etype_u', fi_etype_u)
            ! ii=ii+1
            ! call debug_singlefile_open (debugfile(ii), 'faceI','etype_d', fi_etype_d)
            ! ii=ii+1
            ! call debug_singlefile_open (debugfile(ii), 'faceI','branch_u', fi_branch_u)
            ! ii=ii+1
            ! call debug_singlefile_open (debugfile(ii), 'faceI','branch_d', fi_branch_d)
            ! ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'faceI','jump_type', fi_jump_type)
            ii=ii+1
            ! call debug_singlefile_open (debugfile(ii), 'faceI','node_ID', fi_node_ID)
            ! ii=ii+1
            ! call debug_singlefile_open (debugfile(ii), 'faceI','link_ID', fi_link_ID)
            ! ii=ii+1
            ! call debug_singlefile_open (debugfile(ii), 'faceI','link_Pos', fi_link_Pos)
            ! ii=ii+1
        endif

        if (setting%Debugout%faceYN) then
            print *, 'error: code for debug output of faceYN not completed in ',subroutine_name
            stop
        endif

        if (setting%Debugout%nodeR) then
            call debug_singlefile_open (debugfile(ii), 'nodeR','Eta', nr_Eta)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'nodeR','Depth', nr_Depth)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'nodeR','Volume', nr_Volume)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'nodeR','LateralInflow', nr_LateralInflow)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'nodeR','TotalInflow', nr_TotalInflow)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'nodeR','Flooding', nr_Flooding)
            ii=ii+1
        endif

        if (setting%Debugout%nodeI) then
            print *, 'error: code for debug output of nodeI not completed in ',subroutine_name
            stop
        endif

        if (setting%Debugout%nodeYN) then
            print *, 'error: code for debug output of nodeYN not completed in ',subroutine_name
            stop
        endif

        if (setting%Debugout%linkR) then
            call debug_singlefile_open (debugfile(ii), 'linkR','Flowrate', lr_Flowrate)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'linkR','Depth', lr_Depth)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'linkR','Volume', lr_Volume)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'linkR','Velocity', lr_Velocity)
            ii=ii+1
            call debug_singlefile_open (debugfile(ii), 'linkR','Capacity', lr_Capacity)
            ii=ii+1
        endif

        if (setting%Debugout%linkI) then
            print *, 'error: code for debug output of linkI not completed in ',subroutine_name
            stop
        endif

        if (setting%Debugout%linkYN) then
            print *, 'error: code for debug output of linkYN not completed in ',subroutine_name
            stop
        endif

        if (ii-1 /= ndebug) then
            print *, 'error: There were ',ii-1,' debug files open, but the debugfile array size was ',&
                ndebug, ' in ',subroutine_name
        endif

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine debug_initialize
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine debug_output &
        (debugfile, nodeR, linkR, elem2R, elem2I, elem2YN, elemMR, elemMI, &
        elemMYN, faceR, faceI, faceYN, bcdataUp, bcdataDn, thisstep)
        !
        !   writes the output to debug files
        !   call this subroutine anytime you want debug output
        !
        character(len=64) :: subroutine_name = 'debug_output'

        real(8),      target,     intent(in) :: nodeR(:,:),   linkR(:,:)
        real(8),      target,     intent(in) :: elem2R(:,:),  elemMR(:,:),  faceR(:,:)

        integer,   target,     intent(in) :: elem2I(:,:),  elemMI(:,:),  faceI(:,:)
        logical,   target,     intent(in) :: elem2YN(:,:), elemMYN(:,:), faceYN(:,:)


        type(bcType),  intent(in)  :: bcdataUp(:), bcdataDn(:)

        type(debugfileType), target,  intent(in)  :: debugfile(:)

        integer,       intent(in)  :: thisstep

        character(len=32), pointer ::  ArrayName

        integer,       pointer ::  CI, UnitNumber
        real(8),       pointer ::  thisdataR(:)
        integer,       pointer ::  thisdataI(:)

        integer :: arrayContains= 0  ! =1 for real(8), 2 for integer, 3 for logical

        integer, parameter :: dataR = 1
        integer, parameter :: dataI = 2
        integer, parameter :: dataYN = 3

        integer :: mm

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% skip this if file outpus are globally suppressed
        if (setting%Debugout%SuppressAllFiles) return

        if (setting%Debugout%DisplayInterval <= 0) return

        if (mod(thisstep,setting%Debugout%DisplayInterval) /= 0) return

        do mm=1,size(debugfile)

            UnitNumber => debugfile(mm)%FileInfo%UnitNumber
            ArrayName  => debugfile(mm)%ArrayName
            CI         => debugfile(mm)%ColumnIndex

            select case (ArrayName)
              case ('elem2R')
                thisdataR => elem2R(:,CI)
                arrayContains = dataR
              case ('elemMR')
                thisdataR => elemMR(:,CI)
                arrayContains = dataR
              case ('faceR')
                thisdataR => faceR(:,CI)
                arrayContains = dataR
              case ('elem2I')
                thisdataI => elem2I(:,CI)
                arrayContains = dataI
              case ('elemMI')
                print *, 'error: case not handled for ArrayName of ', ArrayName,' in ',subroutine_name
                stop
              case ('faceI')
                thisdataI => faceI(:,CI)
                arrayContains = dataI
              case ('elem2YN')
                print *, 'error: case not handled for ArrayName of ', ArrayName,' in ',subroutine_name
                stop
              case ('elemMYN')
                print *, 'error: case not handled for ArrayName of ', ArrayName,' in ',subroutine_name
                stop
              case ('faceYN')
                print *, 'error: case not handled for ArrayName of ', ArrayName,' in ',subroutine_name
                stop
              case ('nodeR')
                thisdataR => nodeR(:,CI)
                arrayContains = dataR
              case ('linkR')
                thisdataR => linkR(:,CI)
                arrayContains = dataR
              case default
                print *, 'error: unknown case for ArrayName of ',ArrayName,' in ',subroutine_name
                stop
            end select

            !%  write the time info for this slice
            if (.not. setting%Debugout%SuppressTimeStep)   write(UnitNumber,*) 'step=',setting%Step%Current
            if (.not. setting%Debugout%SuppressTimeValue)  write(UnitNumber,*) 'time=',setting%Time%ThisTime


            !%  write the number of data points and the data itself
            select case (arrayContains)
              case (dataR)
                if (.not. setting%Debugout%SuppressNdat) write(UnitNumber,*) 'ndat=',size(thisdataR)
                write(UnitNumber,*) thisdataR
              case (dataI)
                if (.not. setting%Debugout%SuppressNdat) write(UnitNumber,*) 'ndat=',size(thisdataI)
                write(UnitNumber,*) thisdataI
              case (dataYN)
                print *, 'error: case not handled for ArrayName of dataYN in ',subroutine_name
                stop
              case default
                print *, 'error: unknown case for arrayContains of ',arrayContains,' in ',subroutine_name
                stop
            end select

        end do

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine debug_output
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine debug_finalize &
        (debugfile)
        !
        !   closes debug files before finishing
        !
        character(64) :: subroutine_name = 'debug_finalize'

        type(debugfileType),  intent(in out) :: debugfile(:)

        integer :: mm

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% skip this if files are suppressed
        if (setting%Debugout%SuppressAllFiles) return

        do mm=1,size(debugfile)
            if (debugfile(mm)%FileInfo%IsOpen) then
                close(debugfile(mm)%FileInfo%UnitNumber)
                debugfile(mm)%FileInfo%IsOpen = .false.
            endif
        end do

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine debug_finalize
    !
    !==========================================================================
    !
    ! PRIVATE BELOW HERE
    !
    !==========================================================================
    !
    subroutine debug_singlefile_open &
        (debugfile, ArrayName, DataName, ColumnIndex)
        !
        ! Opens a single debug output file for writing.
        ! Will not overwrite any file - instead gives an error if file is found.
        ! Debug files are given a timestamp of yyyymmdd_hhmm, so overwrites failure
        ! typically only occurs when two runs are made quickly in the same directory.
        !
        ! On MacOS, we have seen random overwrite failures when the file actually
        ! does not exist. It seems likely there is a problem in the file directory of
        ! the Mac we were using.
        !
        character(64) :: subroutine_name = 'debug_singlefile_open'

        type(debugfileType), intent(in out) :: debugfile

        character(len=*),  intent(in)  :: ArrayName, DataName

        integer,           intent(in)  :: ColumnIndex

        integer                :: open_status
        character(len=512)     :: emsg

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        open_status = 0

        debugfile%FileInfo%UnitNumber = outputfile_next_unitnumber
        outputfile_next_unitnumber = outputfile_next_unitnumber+1

        debugfile%FileInfo%FileName   = trim(setting%Debugout%Filename)//'_'//trim(ArrayName)//'_'//trim(DataName)
        debugfile%FileInfo%FolderPath = trim(setting%Debugout%FolderPath)//trim(setting%Debugout%FolderName)//'/'
        debugfile%FileInfo%FileStatus = 'new'
        debugfile%FileInfo%IsOpen     = .true.
        debugfile%ArrayName           = trim(ArrayName)
        debugfile%ColumnIndex         = ColumnIndex

        debugfile%FileInfo%WriteName  = trim(debugfile%FileInfo%FolderPath) // &
            trim(debugfile%FileInfo%FileName) // &
            trim(setting%Time%DateTimeStamp) //&
            '.txt'

        open(unit=debugfile%FileInfo%Unitnumber, &
            file=trim(debugfile%FileInfo%WriteName), &
            status = 'new', &
            access = 'sequential', &
            form   = 'formatted', &
            action = 'write', &
            iostat = open_status)

        emsg = 'file exists or path/folder does not exist: file open failed in '//trim(subroutine_name) &
            // '; filename = '//trim(debugfile%FileInfo%WriteName)
        call utility_check_fileopen (open_status, emsg)


        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine debug_singlefile_open
    !
    !==========================================================================
    ! END OF MODULE debug
    !==========================================================================
end module debug
