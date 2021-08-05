Module utility_output

    use define_indexes
    use define_keys
    use define_globals
    use define_settings
    use define_types
    use interface
    use utility_datetime
    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    private

    public :: util_output_create_folder
    public :: util_output_create_elemR_files
    public :: util_output_create_faceR_files
    public :: util_output_create_summary_files
    public :: util_output_write_elemR_faceR
    public :: util_output_report
    public :: util_output_report_summary

contains

    subroutine util_output_create_folder

        !creates and empties the folder before creating the debug files

        if( this_image() == 1) then
            call system('rm -r debug_output')
            call system('mkdir debug_output')
            call system('mkdir debug_output/elemR')
            call system('mkdir debug_output/faceR')
            call system('mkdir debug_output/summary')
        end if

        sync all

    end subroutine util_output_create_folder

    subroutine util_output_create_summary_files
        integer :: fu, open_status
        character(64) :: file_name

        write(file_name, "(A,i1,A)") "debug_output/summary/summary_", this_image(), ".csv"

        open(newunit=fu, file = file_name, status = 'replace',access = 'sequential', &
        form   = 'formatted', action = 'write', iostat = open_status)

        write(fu, *) "In_Image,This_Time,CFL_max,dt,Velocity_Max,Wavespeed_Max"
            endfile(fu)
        close(fu)
    end subroutine util_output_create_summary_files

    subroutine util_output_create_elemR_files

        integer :: fu, open_status, ii
        character(len = 250) :: file_name
        character(len = 40)  :: dir
        character(len = 4)   :: str_image
        character(len = 100) :: link_name
        character(len = 40)  :: str_elem_idx

        fu = this_image()

        write(str_image, '(i1)') fu

        do ii = 1, N_elem(this_image())

            write(str_elem_idx,'(I10)') elemI(ii,ei_Gidx)

            if(elemI(ii,ei_elementType) == CC) then
                file_name = "debug_output/elemR/"//trim(str_image)//"_CC_" &
                    // trim(link%names(elemI(ii,ei_link_Gidx_SWMM))%str) // &
                    "_" // trim(ADJUSTL(str_elem_idx))//".csv"

                open(newunit=fu, file = file_name, status = 'replace',access = 'sequential', &
                form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                end if


            else if(elemI(ii,ei_elementType) == JM) then
                file_name = "debug_output/elemR/"//trim(str_image)//"_JM_" &
                    // trim(node%names(elemI(ii,ei_node_Gidx_SWMM))%str) // &
                    "_" // trim(ADJUSTL(str_elem_idx))//".csv"


                open(newunit=fu, file = file_name, status = 'replace',access = 'sequential', &
                form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                end if


            else if(elemI(ii,ei_elementType) == JB) then
                file_name = "debug_output/elemR/"//trim(str_image)//"_JB_" &
                    // trim(node%names(elemI(ii,ei_node_Gidx_SWMM))%str) // &
                    "_" // trim(ADJUSTL(str_elem_idx))//".csv"


                open(newunit=fu, file = file_name, status = 'replace',access = 'sequential', &
                    form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                end if

            end if

            write(fu, *) "Timestamp,Time_In_Secs,Area,Area_N0,Area_N1,AreaBelowBreadthMax,BreadthMax,Depth,dHdA,ell,&
                Flowrate,Flowrate_N0,Flowrate_N1,FlowrateLateral,FlowrateStore,FroudeNumber,&
                FullArea,FullDepth,FullHydDepth,FullPerimeter,FullVolume,GammaC,GammaM,Head,&
                Head_N0,HeadLastAC,HeadStore,HydDepth,HydRadius,InterpWeight_uG,InterpWeight_dG,&
                InterpWeight_uH,InterpWeight_dH,InterpWeight_uQ,InterpWeight_dQ,Ksource,Length,ones,&
                Perimeter,Roughness,SmallVolume,SmallVolume_CMvelocity,SmallVolume_HeadSlope,&
                SmallVolume_ManningsN,SmallVolumeRatio,SourceContinuity,SourceMomentum,Temp01,&
                Topwidth,Velocity,Velocity_N0,Velocity_N1,VelocityLastAC,Volume,Volume_N0,&
                Volume_N1,VolumeLastAC,VolumeStore,WaveSpeed,Zbottom,ZbreadthMax,Zcrown"
            endfile(fu)
            close(fu)

        end do

    end subroutine util_output_create_elemR_files

    subroutine util_output_create_faceR_files


        integer :: fu, open_status, ii
        character(len = 250) :: file_name
        character(len = 40)  :: dir
        character(len = 4)   :: str_image
        character(len = 100) :: link_name
        character(len = 40)  :: str_face_idx

        fu = this_image()

        write(str_image, '(i1)') fu

        do ii = 1, N_face(this_image())

            write(str_face_idx,'(I10)') faceI(ii,fi_Gidx)

            file_name = "debug_output/faceR/"//trim(str_image)//"_face_" &
                // trim(ADJUSTL(str_face_idx))//".csv"

            open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
                form   = 'formatted', action = 'write', iostat = open_status)

            if (open_status /= 0) then
                write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
            end if


            write(fu, *) "Timestamp, Time_In_Secs, Area_d, Area_u, Flowrate, Flowrate_N0, Head_u, Head_d,"// &
                "HydDepth_d, HydDepth_u, Topwidth_d, Topwidth_u, Velocity_d, Velocity_u"
            endfile(fu)

            close(fu)

        end do


    end subroutine util_output_create_faceR_files


    subroutine util_output_write_elemR_faceR

        integer :: fu, open_status, ii, yr, mnth, dy, hr, min, sec
        real(8) :: time_secs, time_epoch
        character(len = 250) :: file_name
        character(len = 40)  :: dir
        character(len = 4)   :: str_image
        character(len = 100) :: link_name
        character(len = 40)  :: str_elem_face_idx


        fu = this_image()
        time_secs = setting%Time%Hydraulics%timeNow
        time_epoch = util_datetime_secs_to_epoch(time_secs)
        call util_datetime_decodedate(time_epoch, yr, mnth, dy)
        call util_datetime_decodetime(time_epoch, hr, min, sec)

        write(str_image, '(i1)') fu


        do ii = 1, N_elem(this_image())

            write(str_elem_face_idx,'(I10)') elemI(ii,ei_Gidx)

            if(elemI(ii,ei_elementType) == CC) then
                file_name = "debug_output/elemR/"//trim(str_image)//"_CC_" &
                    // trim(link%names(elemI(ii,ei_link_Gidx_SWMM))%str) // &
                    "_" // trim(ADJUSTL(str_elem_face_idx))//".csv"

                open(newunit=fu, file = file_name, status = 'old',access = 'append', &
                form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                end if


            else if(elemI(ii,ei_elementType) == JM) then
                file_name = "debug_output/elemR/"//trim(str_image)//"_JM_" &
                    // trim(node%names(elemI(ii,ei_node_Gidx_SWMM))%str) // &
                    "_" // trim(ADJUSTL(str_elem_face_idx))//".csv"


                open(newunit=fu, file = file_name, status = 'old',access = 'append', &
                form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                end if


            else if(elemI(ii,ei_elementType) == JB) then
                file_name = "debug_output/elemR/"//trim(str_image)//"_JB_" &
                    // trim(node%names(elemI(ii,ei_node_Gidx_SWMM))%str) // &
                    "_" // trim(ADJUSTL(str_elem_face_idx))//".csv"


                open(newunit=fu, file = file_name, status = 'old',access = 'append', &
                    form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                end if

            end if
            write(fu,fmt='(i4, 2(a,i2.2))',advance = 'no') yr,"/",mnth,"/",dy
            write(fu,fmt = '(A)',advance = 'no') ' '
            write(fu,fmt='(2(i2.2,a), i2.2)',advance = 'no') hr,":",min,":",sec
            write(fu,'(A)', advance = 'no') ', '
            write(fu, '(F32.16)', advance = 'no') time_secs
            write(fu,'(A)', advance = 'no') ', '
            write(fu, '(*(G0.6,:,","))') elemR(ii,:)
            endfile(fu)
            close(fu)
        end do



        do ii = 1, N_face(this_image())


            write(str_elem_face_idx,'(I10)') faceI(ii,fi_Gidx)

            file_name = "debug_output/faceR/"//trim(str_image)//"_face_" &
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
            write(fu,'(A)', advance = 'no') ', '
            write(fu, '(F32.16)', advance = 'no') time_secs
            write(fu,'(A)', advance = 'no') ', '
            write(fu, '(*(G0.6,:,","))') faceR(ii, :)

            endfile(fu)
            close(fu)

        end do

    end subroutine util_output_write_elemR_faceR

    subroutine util_output_report
        character(64) :: subroutine_name = "util_output_report"

        if (setting%Debug%File%utility_output) print *, '*** enter ', subroutine_name

        if (util_output_must_report()) then
            call util_output_write_elemR_faceR()
        end if

        if (setting%Debug%File%utility_output) print *, '*** leave ', subroutine_name
    end subroutine util_output_report

    subroutine util_output_report_summary()
        integer :: fu, open_status
        real(8) :: thisCFL
        real(8), pointer :: dt, timeNow, velocity(:), wavespeed(:), length(:)
        integer, pointer :: thisCol, Npack, thisP(:)
        character(64) :: file_name
        character(64) :: subroutine_name = "util_output_report_summary"

        if (setting%Debug%File%utility_output) print *, '*** enter ', subroutine_name

        if (util_output_must_report()) then
            write(file_name, "(A,i1,A)") "debug_output/summary/summary_", this_image(), ".csv"
            timeNow   => setting%Time%Hydraulics%timeNow
            dt        => setting%Time%Hydraulics%Dt
            velocity  => elemR(:,er_Velocity)
            wavespeed => elemR(:,er_WaveSpeed)
            length    => elemR(:,er_Length)
            thisCol   => col_elemP(ep_CC_ALLtm)
            Npack     => npack_elemP(thisCol)
            thisP     => elemP(1:Npack,thisCol)
            thisCFL = maxval((velocity(thisP) + wavespeed(thisP)) * dt / length(thisP))

            if (setting%verbose) then
                print*, '--------------------------------------'
                print*, 'In image', this_image()
                print*, 'This Time = ', timeNow, 'dt = ', dt
                print*, 'CFL max = ', thisCFL, 'Velocity Max = ', maxval(abs(velocity(thisP))) , &
                'Wavespeed max = ', maxval(abs(wavespeed(thisP)))
            end if

            open(newunit=fu, file = trim(file_name), status = 'old',access = 'Append', &
                form = 'formatted', action = 'write', iostat = open_status)
            write(fu, fmt='(*(G0.6 : ","))') &
                this_image(), timeNow, thisCFL, dt, maxval(abs(velocity(thisP))), maxval(abs(wavespeed(thisP)))
            endfile(fu)
            close(fu)
        end if
        if (setting%Debug%File%utility_output) print *, '*** leave ', subroutine_name
    end subroutine util_output_report_summary

    function util_output_must_report() result(report)
        logical :: report
        real(8), pointer :: timeNow
        timeNow => setting%Time%Hydraulics%timeNow
        report = ((abs(mod(timeNow, setting%output%report_time) - &
             setting%output%report_time) <= setting%output%report_tol) .or. &
             (timeNow == 0))
    end function util_output_must_report

    subroutine util_output_debug_elemI

        integer :: fu, open_status
        character(len = 250) :: file_name
        character(len = 4) :: str_image

        fu = this_image()

        write(str_image, '(i1)') fu

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemI_ei_Lidx_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemI(:,ei_Lidx)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemI_ei_Gidx_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemI(:,ei_Gidx)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemI_ei_elementType_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemI(:,ei_elementType)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemI_ei_geometryType_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemI(:,ei_geometryType)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemI_ei_HeqType_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemI(:,ei_HeqType)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemI_ei_link_Gidx_SWMM_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemI(:,ei_link_Gidx_SWMM)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemI_ei_link_Gidx_BIPquick_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemI(:,ei_link_Gidx_BIPquick)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemI_ei_link_pos_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemI(:,ei_link_pos)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemI_ei_Mface_uL_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemI(:,ei_Mface_uL)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemI_ei_Mface_dL_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemI(:,ei_Mface_dL)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemI_ei_node_Gidx_SWMM_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemI(:, ei_node_Gidx_SWMM)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemI_ei_node_Gidx_BIPquick_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemI(:, ei_node_Gidx_BIPquick)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemI_ei_QeqType_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemI(:, ei_QeqType)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemI_ei_specificType_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemI(:, ei_specificType)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemI_ei_Temp01_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemI(:,ei_Temp01)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemI_ei_tmType_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemI(:,ei_tmType)

        close(fu)


    end subroutine util_output_debug_elemI

    subroutine util_output_debug_elemR

        integer fu, open_status
        character(len = 250) :: file_name
        character(len = 4) :: str_image

        fu = this_image()

        write(str_image, '(i1)') fu

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Area_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Area)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Area_N0_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Area_N0)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Area_N1_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Area_N1)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_AreaBelowBreadthMax_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_AreaBelowBreadthMax)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_BreadthMax_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_BreadthMax)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Depth_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Depth)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_dHdA_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_dHdA)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_ell_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_ell)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Flowrate_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Flowrate)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Flowrate_N0_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Flowrate_N0)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Flowrate_N1_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Flowrate_N1)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_FlowrateLateral_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_FlowrateLateral)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_FlowrateStore_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_FlowrateStore)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_FroudeNumber_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_FroudeNumber)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_FullArea_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_FullArea)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_FullDepth_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_FullDepth)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_FullHydDepth_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_FullHydDepth)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_FullPerimeter_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_FullPerimeter)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_FullVolume_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_FullVolume)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_GammaC_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_GammaC)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_GammaM_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_GammaM)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Head_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Head)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Head_N0_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Head_N0)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_HeadLastAC_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_HeadLastAC)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_HeadStore_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_HeadStore)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_HydDepth_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_HydDepth)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_HydRadius_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_HydRadius)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_InterpWeight_uG_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_InterpWeight_uG)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_InterpWeight_dG_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_InterpWeight_dG)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_InterpWeight_uH_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_InterpWeight_uH)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_InterpWeight_dH_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_InterpWeight_dH)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_InterpWeight_uQ_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_InterpWeight_uQ)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_InterpWeight_dQ_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_InterpWeight_dQ)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Ksource_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Ksource)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Length_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Length)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_ones_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_ones)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Perimeter_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Perimeter)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Roughness_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Roughness)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_SmallVolume_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_SmallVolume)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_SmallVolume_CMvelocity_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_SmallVolume_CMvelocity)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_SmallVolume_HeadSlope_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_SmallVolume_HeadSlope)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_SmallVolume_ManningsN_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_SmallVolume_ManningsN)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_SmallVolumeRatio_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_SmallVolumeRatio)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_SourceContinuity_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_SourceContinuity)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_SourceMomentum_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_SourceMomentum)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Temp01_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Temp01)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Topwidth_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Topwidth)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Velocity_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Velocity)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Velocity_N0_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Velocity_N0)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Velocity_N1_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Velocity_N1)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_VelocityLastAC_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_VelocityLastAC)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Volume_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Volume)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Volume_N0_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Volume_N0)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Volume_N1_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Volume_N1)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_VolumeLastAC_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_VolumeLastAC)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_VolumeStore_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_VolumeStore)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_WaveSpeed_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_WaveSpeed)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Zbottom_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Zbottom)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_ZbreadthMax_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_ZbreadthMax)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemR_er_Zcrown_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemR(:,er_Zcrown)

        close(fu)

    end subroutine util_output_debug_elemR


    subroutine util_output_debug_faceI

        integer fu, open_status
        character(len = 250) :: file_name
        character(len = 4) :: str_image

        fu = this_image()

        write(str_image, '(i1)') fu

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceI_fi_Lidx_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceI(:,fi_Lidx)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceI_fi_Gidx_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceI(:,fi_Gidx)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceI_fi_BCtype_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceI(:,fi_BCtype)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceI_fi_jump_type_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceI(:,fi_jump_type)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceI_fi_Melem_uL_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceI(:,fi_Melem_uL)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceI_fi_Melem_dL_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceI(:,fi_Melem_dL)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceI_fi_GhostElem_uL_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceI(:,fi_GhostElem_uL)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceI_fi_Melem_dL_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceI(:,fi_Melem_dL)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceI_fi_Connected_image_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceI(:,fi_Connected_image)

        close(fu)


    end subroutine util_output_debug_faceI

    subroutine util_output_debug_faceR

        integer fu, open_status
        character(len = 250) :: file_name
        character(len = 4) :: str_image

        fu = this_image()

        write(str_image, '(i1)') fu

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceR_fr_Area_d_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceR(:,fr_Area_d)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceR_fr_Area_u_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceR(:,fr_Area_u)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceR_fr_Flowrate_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceR(:,fr_Flowrate)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceR_fr_Flowrate_N0_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceR(:,fr_Flowrate_N0)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceR_fr_Head_u_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceR(:,fr_Head_u)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceR_fr_Head_d_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceR(:,fr_Head_d)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceR_fr_HydDepth_d_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceR(:,fr_HydDepth_d)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceR_fr_HydDepth_u_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceR(:,fr_HydDepth_u)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceR_fr_Topwidth_d_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceR(:,fr_Topwidth_d)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceR_fr_Topwidth_u_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceR(:,fr_Topwidth_u)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceR_fr_Velocity_d_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceR(:,fr_Velocity_d)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceR_fr_Velocity_u_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceR(:,fr_Velocity_u)

        close(fu)

    end subroutine util_output_debug_faceR


    subroutine util_output_debug_elemYN

        integer fu, open_status
        character(len = 250) :: file_name
        character(len = 4) :: str_image

        fu = this_image()

        write(str_image, '(i1)') fu

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemYN_eYN_canSurcharge_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemYN(:,eYN_canSurcharge)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemYN_eYN_isAdhocFlowrate_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemYN(:,eYN_isAdhocFlowrate)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemYN_eYN_isSmallVolume_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemYN(:,eYN_isSmallVolume)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemYN_eYN_isSurcharged_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemYN(:,eYN_isSurcharged)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemYN_eYN_isNearZeroVolume_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemYN(:,eYN_isNearZeroVolume)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_elemYN_eYN_isDummy_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) elemYN(:,eYN_isDummy)

        close(fu)


    end subroutine util_output_debug_elemYN

    subroutine util_output_debug_faceYN

        integer fu, open_status
        character(len = 250) :: file_name
        character(len = 4) :: str_image

        fu = this_image()

        write(str_image, '(i1)') fu

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceYN_fYN_isAC_adjacent_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceYN(:,fYN_isAC_adjacent)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceYN_fYN_isInteriorFace_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceYN(:,fYN_isInteriorFace)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceYN_fYN_isSharedFace_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceYN(:,fYN_isSharedFace)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceYN_fYN_isUpGhost_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceYN(:,fYN_isUpGhost)

        close(fu)

        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceYN_fYN_isDnGhost_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceYN(:,fYN_isDnGhost)

        close(fu)


        !-----------------------------------------------------------------------------------------------
        file_name = "debug_output/debug_faceYN_fYN_isnull_"//trim(str_image)//".csv"

        open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
            form   = 'formatted', action = 'write', iostat = open_status)

        if (open_status /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
        end if

        write(fu,*) faceYN(:,fYN_isnull)

        close(fu)

    end subroutine util_output_debug_faceYN





end Module utility_output
