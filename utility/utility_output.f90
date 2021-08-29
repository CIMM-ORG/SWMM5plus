Module utility_output

    use define_indexes
    use define_keys
    use define_globals
    use define_settings
    use define_types
    use interface
    use utility_datetime
    use output
    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    private

    public :: util_output_clean_folders
    public :: util_output_create_folder
    public :: util_output_create_elemR_files
    public :: util_output_create_faceR_files
    public :: util_output_create_summary_files
    public :: util_output_write_elemR_faceR
    public :: util_output_report
    public :: util_output_must_report

contains

    subroutine util_output_clean_folders
        character(64) :: subroutine_name = 'util_output_clean_folders'

        if (setting%Debug%File%utility_output) print *, "*** enter ", this_image(), subroutine_name

        if (this_image() == 1) then
            call system('mkdir -p debug_output')
            call system('rm -r debug_output')
            call system('mkdir -p debug')
            call system('rm -r debug')
        end if

        if (setting%Debug%File%utility_output) print *, "*** leave ", this_image(), subroutine_name

    end subroutine util_output_clean_folders

    subroutine util_output_create_folder
        character(64) :: subroutine_name = 'util_output_create_folder'

        if (setting%Debug%File%utility_output) print *, "*** enter ", this_image(), subroutine_name
        !creates and empties the folder before creating the debug files

        if( this_image() == 1) then
            call system('mkdir debug_output')
            call system('mkdir debug_output/partitioned')
            call system('mkdir debug_output/swmm5')
            call system('mkdir debug_output/partitioned/elemR')
            call system('mkdir debug_output/partitioned/faceR')
            call system('mkdir debug_output/partitioned/summary')
            call system('mkdir debug_output/partitioned/link')
            call system('mkdir debug_output/partitioned/node')
            call system('mkdir debug_output/swmm5/node')
            call system('mkdir debug_output/swmm5/link')
        end if

        sync all

        if (setting%Debug%File%utility_output) print *, "*** leave ", this_image(), subroutine_name
    end subroutine util_output_create_folder

    subroutine util_output_create_summary_files
        integer :: fu, open_status
        character(64) :: file_name
        character(64) :: subroutine_name = 'util_output_create_summary_files'

        if (setting%Debug%File%utility_output) print *, "*** enter ", this_image(), subroutine_name

        call system("mkdir -p debug_output")

        write(file_name, "(A,i1,A)") "debug_output/partitioned/summary/summary_", this_image(), ".csv"

        open(newunit=fu, file = file_name, status = 'replace',access = 'sequential', &
        form = 'formatted', action = 'write', iostat = open_status)

        write(fu, *) "In_Image,This_Time,CFL_max,dt,Velocity_Max,Wavespeed_Max"
        endfile(fu)
        close(fu)

        if (setting%Debug%File%utility_output) print *, "*** leave ", this_image(), subroutine_name
    end subroutine util_output_create_summary_files

    subroutine util_output_create_elemR_files

        integer :: fu, open_status, ii
        character(len = 250) :: file_name
        character(len = 40)  :: dir
        character(len = 4)   :: str_image
        character(len = 100) :: link_name
        character(len = 40)  :: str_elem_idx
        character(len = 10)  :: str_link_node_idx
        character(64) :: subroutine_name = 'util_output_create_elemR_files'

        if (setting%Debug%File%utility_output) print *, "*** enter ", this_image(), subroutine_name

        fu = this_image()

        write(str_image, '(i1)') fu

        !print *, "N_elem", N_elem(this_image())
        !print *, "size(elemI(:,ei_Gidx)",size(elemI(:,ei_Gidx))
        !print *, "size(link%names(:))", size(link%names(:))

        do ii = 1, N_elem(this_image())

            write(str_elem_idx,'(I10)') elemI(ii,ei_Gidx)

            if(elemI(ii,ei_elementType) == CC) then

                write(str_link_node_idx,'(I10)') elemI(ii,ei_link_Gidx_SWMM)

                file_name = "debug_output/partitioned/elemR/"//trim(str_image)//"_CC_" &
                    // trim(ADJUSTL(str_link_node_idx))// &
                    "_" // trim(ADJUSTL(str_elem_idx))//".csv"

                open(newunit=fu, file = file_name, status = 'replace',access = 'sequential', &
                form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                end if

            else if(elemI(ii,ei_elementType) == JM) then

                write(str_link_node_idx,'(I10)') elemI(ii,ei_node_Gidx_SWMM)

                file_name = "debug_output/partitioned/elemR/"//trim(str_image)//"_JM_" &
                    // trim(ADJUSTL(str_link_node_idx))// &
                    "_" // trim(ADJUSTL(str_elem_idx))//".csv"

                open(newunit=fu, file = file_name, status = 'replace',access = 'sequential', &
                form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                end if

            else if(elemI(ii,ei_elementType) == JB) then

                 write(str_link_node_idx,'(I10)') elemI(ii,ei_node_Gidx_SWMM)

                 file_name = "debug_output/partitioned/elemR/"//trim(str_image)//"_JB_" &
                     // trim(ADJUSTL(str_link_node_idx))// &
                     "_" // trim(ADJUSTL(str_elem_idx))//".csv"


                open(newunit=fu, file = file_name, status = 'replace',access = 'sequential', &
                    form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                end if

            end if

            write(fu, "(A)") "Timestamp,Time_In_Secs,Area,Area_N0,Area_N1,AreaBelowBreadthMax,BreadthMax,Depth,dHdA,ell," //&
                "Flowrate,Flowrate_N0,Flowrate_N1,FlowrateLateral,FlowrateStore,FroudeNumber," // &
                "FullArea,FullDepth,FullHydDepth,FullPerimeter,FullVolume,GammaC,GammaM,Head," // &
                "Head_N0,HeadLastAC,HeadStore,HydDepth,HydRadius,InterpWeight_uG,InterpWeight_dG," // &
                "InterpWeight_uH,InterpWeight_dH,InterpWeight_uQ,InterpWeight_dQ,Ksource,Length,ones," // &
                "Perimeter,Roughness,SmallVolume,SmallVolume_CMvelocity,SmallVolume_HeadSlope," // &
                "SmallVolume_ManningsN,SmallVolumeRatio,SourceContinuity,SourceMomentum,Temp01," // &
                "Topwidth,Velocity,Velocity_N0,Velocity_N1,VelocityLastAC,Volume,Volume_N0," // &
                "Volume_N1,VolumeLastAC,VolumeStore,WaveSpeed,Zbottom,ZbreadthMax,Zcrown"
            endfile(fu)
            close(fu)

        end do
        if (setting%Debug%File%utility_output) print *, "*** leave ", this_image(), subroutine_name

    end subroutine util_output_create_elemR_files

    subroutine util_output_create_faceR_files


        integer :: fu, open_status, ii
        character(len = 250) :: file_name
        character(len = 40)  :: dir
        character(len = 4)   :: str_image
        character(len = 100) :: link_name
        character(len = 40)  :: str_face_idx
        character(64) :: subroutine_name = 'util_output_create_faceR_files'


        if (setting%Debug%File%utility_output) print *, "*** enter ", this_image(), subroutine_name

        fu = this_image()

        write(str_image, '(i1)') fu

        do ii = 1, N_face(this_image())

            write(str_face_idx,'(I10)') faceI(ii,fi_Gidx)

            file_name = "debug_output/partitioned/faceR/"//trim(str_image)//"_face_" &
                // trim(ADJUSTL(str_face_idx))//".csv"

            open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
                form   = 'formatted', action = 'write', iostat = open_status)

            if (open_status /= 0) then
                write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
            end if


            write(fu, *) "Timestamp, Time_In_Secs, Area_d, Area_u, Flowrate, Flowrate_N0, Head_u, Head_d,"// &
                "Zbottom, HydDepth_d, HydDepth_u, Topwidth_d, Topwidth_u, Velocity_d, Velocity_u"
            endfile(fu)

            close(fu)

        end do

        if (setting%Debug%File%utility_output) print *, "*** leave ", this_image(), subroutine_name

    end subroutine util_output_create_faceR_files

    subroutine util_output_write_elemR_faceR

        integer :: fu, open_status, ii, yr, mnth, dy, hr, min, sec
        real(8) :: time_secs, time_epoch
        character(len = 250) :: file_name
        character(len = 40)  :: dir
        character(len = 4)   :: str_image
        character(len = 100) :: link_name
        character(len = 40)  :: str_elem_face_idx
        character(len = 10)  :: str_link_node_idx
        character(64) :: subroutine_name = 'util_output_write_elemR_faceR'

        if (setting%Debug%File%utility_output) print *, "*** enter ", this_image(), subroutine_name

        fu = this_image()
        time_secs = setting%Time%Now
        time_epoch = util_datetime_secs_to_epoch(time_secs)
        call util_datetime_decodedate(time_epoch, yr, mnth, dy)
        call util_datetime_decodetime(time_epoch, hr, min, sec)

        write(str_image, '(i1)') fu

        do ii = 1, N_elem(this_image())

            write(str_elem_face_idx,'(I10)') elemI(ii,ei_Gidx)

            if(elemI(ii,ei_elementType) == CC) then

                write(str_link_node_idx,'(I10)') elemI(ii,ei_link_Gidx_SWMM)

                file_name = "debug_output/partitioned/elemR/"//trim(str_image)//"_CC_" &
                    // trim(ADJUSTL(str_link_node_idx))// &
                    "_" // trim(ADJUSTL(str_elem_face_idx))//".csv"

                open(newunit=fu, file = file_name, status = 'old',access = 'append', &
                form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                end if


            else if(elemI(ii,ei_elementType) == JM) then

                write(str_link_node_idx,'(I10)') elemI(ii,ei_node_Gidx_SWMM)

                file_name = "debug_output/partitioned/elemR/"//trim(str_image)//"_JM_" &
                    // trim(ADJUSTL(str_link_node_idx))// &
                    "_" // trim(ADJUSTL(str_elem_face_idx))//".csv"


                open(newunit=fu, file = file_name, status = 'old',access = 'append', &
                form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                end if


            else if(elemI(ii,ei_elementType) == JB) then

                write(str_link_node_idx,'(I10)') elemI(ii,ei_node_Gidx_SWMM)

                file_name = "debug_output/partitioned/elemR/"//trim(str_image)//"_JB_" &
                    // trim(ADJUSTL(str_link_node_idx))// &
                    "_" // trim(ADJUSTL(str_elem_face_idx))//".csv"


                open(newunit=fu, file = file_name, status = 'old',access = 'append', &
                    form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                end if

            end if
            write(fu,fmt='(i4,2(a,i2.2))',advance = 'no') yr,"/",mnth,"/",dy
            write(fu,fmt = '(A)',advance = 'no') ' '
            write(fu,fmt='(2(i2.2,a),i2.2)',advance = 'no') hr,":",min,":",sec
            write(fu,'(A)', advance = 'no') ','
            write(fu, '(F0.16)', advance = 'no') time_secs
            write(fu,'(A)', advance = 'no') ','
            write(fu, '(*(G0.6,:,","))') elemR(ii,:)
            endfile(fu)
            close(fu)
        end do



        do ii = 1, N_face(this_image())


            write(str_elem_face_idx,'(I10)') faceI(ii,fi_Gidx)

            file_name = "debug_output/partitioned/faceR/"//trim(str_image)//"_face_" &
                // trim(ADJUSTL(str_elem_face_idx))//".csv"

            open(newunit=fu, file = file_name, status = 'old',access = 'Append', &
                form   = 'formatted', action = 'write', iostat = open_status)

            if (open_status /= 0) then
                write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
            end if

            !write the data to the file
            write(fu,fmt='(i4, 2(a,i2.2))',advance = 'no') yr,"/",mnth,"/",dy
            write(fu,fmt = '(A)',advance = 'no') ' '
            write(fu,fmt='(2(i2.2,a), i2.2)',advance = 'no') hr,":",min,":",sec
            write(fu,'(A)', advance = 'no') ','
            write(fu, '(F0.16)', advance = 'no') time_secs
            write(fu,'(A)', advance = 'no') ','
            write(fu, '(*(G0.6,:,","))') faceR(ii, :)

            endfile(fu)
            close(fu)

        end do
        if (setting%Debug%File%utility_output) print *, "*** leave ", this_image(), subroutine_name

    end subroutine util_output_write_elemR_faceR

    subroutine util_output_report
        character(64) :: subroutine_name = "util_output_report"

        if (setting%Debug%File%utility_output) print *, '*** enter ', this_image(), subroutine_name

        call util_output_report_summary()

        if (setting%Output%report .and. util_output_must_report()) then
            !call util_output_write_elemR_faceR()
            call output_write_link_files()
            call output_write_node_files()
        end if

        if (setting%Debug%File%utility_output) print *, '*** leave ', this_image(), subroutine_name
    end subroutine util_output_report

    subroutine util_output_report_summary()
        integer          :: fu, open_status, thisCol, Npack
        integer, pointer :: thisP(:)
        real(8)          :: thisCFL, max_velocity, max_wavespeed
        real(8), pointer :: dt, timeNow, velocity(:), wavespeed(:), length(:)
        character(512)    :: file_name
        character(64)    :: subroutine_name = "util_output_report_summary"

        if (setting%Debug%File%utility_output) print *, '*** enter ', this_image(), subroutine_name
        if (util_output_must_report() .and. setting%output%report) then

            write(file_name, "(A,i1,A)") "debug_output/partitioned/summary/summary_", this_image(), ".csv"

            thisCol   = col_elemP(ep_CC_ALLtm)
            Npack     = npack_elemP(thisCol)

            timeNow   => setting%Time%Now
            dt        => setting%Time%Dt
            velocity  => elemR(:,er_Velocity)
            wavespeed => elemR(:,er_WaveSpeed)
            length    => elemR(:,er_Length)
            thisP     => elemP(1:Npack,thisCol)

            thisCFL       = maxval((velocity(thisP) + wavespeed(thisP)) * dt / length(thisP))
            max_velocity  = maxval(abs(velocity(thisP)))
            max_wavespeed = maxval(abs(wavespeed(thisP)))

            open(newunit=fu, file = trim(file_name), status = 'old',access = 'Append', &
                form = 'formatted', action = 'write', iostat = open_status)
            write(fu, fmt='(*(G0.6 : ","))') &
                this_image(), timeNow, thisCFL, dt, max_velocity, max_wavespeed
            endfile(fu)
            close(fu)

            if (setting%verbose) then
                print*, '--------------------------------------'
                !% also print the summary in the terminal
                print('(*(G0.6))'), 'image = ', this_image(), ',  timeNow = ', timeNow, ',  dt = ', dt
                print('(*(G0.6))'), 'thisCFL = ',thisCFL, ',  max velocity = ', max_velocity, &
                ',  max wavespeed = ', max_wavespeed
            end if
        end if

        if (setting%Debug%File%utility_output) print *, '*** leave ', this_image(), subroutine_name
    end subroutine util_output_report_summary

    function util_output_must_report() result(report)
        logical :: report
        integer, pointer :: reportStep
        real(8) :: timeNow, reportDt, startReport

        reportStep  => setting%Output%reportStep
        timeNow     = setting%Time%Now
        reportDt    = setting%Output%reportDt
        startReport = setting%output%reportStartTime

        if ((timeNow >= reportDt * (reportStep + 1)) .and. (timeNow > startReport))then
            report = .true.
        else if (timeNow == startReport) then
            report = .true.
            reportStep = -1
        else
            report = .false.
        end if
    end function util_output_must_report


end Module utility_output
