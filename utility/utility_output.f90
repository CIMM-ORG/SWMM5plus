module utility_output

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
    public :: util_output_create_folders
    public :: util_output_export_linknode_input
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
            call system('mkdir -p debug_input')
            call system('mkdir -p debug_output')
            call system('mkdir -p swmm5_output')
            call system('rm -r debug_input')
            call system('rm -r debug_output')
            call system('rm -r swmm5_output')
        end if

        if (setting%Debug%File%utility_output) print *, "*** leave ", this_image(), subroutine_name

    end subroutine util_output_clean_folders

    subroutine util_output_create_folders
        character(64) :: subroutine_name = 'util_output_create_folders'

        if (setting%Debug%File%utility_output) print *, "*** enter ", this_image(), subroutine_name
        !creates and empties the folder before creating the debug files

        if ( this_image() == 1) then
            if (setting%Debug%Input) then
                call system('mkdir debug_input')
                call system('mkdir debug_input/node')
                call system('mkdir debug_input/link')
            end if
            if (setting%Debug%Output .or. setting%Output%report) then
                call system('mkdir debug_output')
                call system('mkdir debug_output/link')
                call system('mkdir debug_output/node')
                call system('mkdir swmm5_output')
                call system('mkdir swmm5_output/node')
                call system('mkdir swmm5_output/link')
            end if
            if (setting%Debug%Output) then
                call system('mkdir debug_output/elemR')
                call system('mkdir debug_output/faceR')
                call system('mkdir debug_output/summary')
                !% >>> BEGIN HACK
                call system('mkdir debug_output/swmm5')
                call system('mkdir debug_output/swmm5/link')
                call system('mkdir debug_output/swmm5/node')
                !% >>> END HACK
            end if
        end if

        sync all

        if (setting%Debug%File%utility_output) print *, "*** leave ", this_image(), subroutine_name
    end subroutine util_output_create_folders

    subroutine util_output_export_linknode_input()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Exports link and node tables as CSV files
        !%-----------------------------------------------------------------------------
        integer :: ii
        logical :: ex
        !%-----------------------------------------------------------------------------
        open(unit=1,file='debug_input/link/linkR.csv', status='unknown', action='write')
        open(unit=2,file='debug_input/link/linkI.csv', status='unknown', action='write')
        write(1, '(A)')                                                                    &
            "lr_Length,lr_AdjustedLength,lr_InletOffset,lr_OutletOffset,lr_BreadthScale,"                 // &
            "lr_TopWidth,lr_ElementLength,lr_Slope,lr_LeftSlope,lr_RightSlope,"         // &
            "lr_Roughness,lr_InitialFlowrate,lr_InitialDepth,lr_InitialUpstreamDepth,"  // &
            "lr_InitialDnstreamDepth,lr_ParabolaValue,lr_SideSlope,lr_DischargeCoeff1," // &
            "lr_DischargeCoeff2,lr_FullDepth,lr_Flowrate,lr_Depth,"  // &
            "lr_DepthUp,lr_DepthDn,lr_Volume,lr_Velocity,lr_Capacity"

        write(2, '(A)')                                                                    &
            "li_idx,li_link_type,li_weir_type,li_orif_type,li_pump_type,li_geometry,"    //&
            "li_roughness_type,li_N_element,li_Mnode_u,li_Mnode_d,li_assigned,"         // &
            "li_InitialDepthType,li_length_adjusted,li_P_image,li_parent_link,"     //&
            "li_weir_EndContrations,li_first_elem_idx,li_last_elem_idx"

        do ii = 1, size(link%R, 1)
            write(1,'(*(G0.6,:,","))') link%R(ii,:)
            write(2,'(*(G0.6,:,","))') link%I(ii,:)
        end do

        close(1)
        close(2)

        open(unit=3,file='debug_input/node/nodeR.csv', status='unknown', action='write')
        open(unit=4,file='debug_input/node/nodeI.csv', status='unknown', action='write')
        open(unit=5,file='debug_input/node/nodeYN.csv', status='unknown', action='write')

        write(3, '(A)')                                                                              &
            "nr_Zbottom,nr_InitialDepth,nr_FullDepth,nr_StorageConstant,nr_StorageCoeff,"         // &
            "nr_StorageExponent,nr_PondedArea,nr_SurchargeDepth,nr_MaxInflow,nr_Eta,"             // &
            "nr_Depth,nr_Volume,nr_LateralInflow,nr_TotalInflow,nr_Flooding,nr_ElementLength_u1," // &
            "nr_ElementLength_u2,nr_ElementLength_u3,nr_ElementLength_d1,nr_ElementLength_d2,"    // &
            "nr_ElementLength_d3"

        write(4, '(A)')                                                                                 &
            "ni_idx,ni_node_type,ni_N_link_u,ni_N_link_d,ni_curve_type,ni_assigned,"                 // &
            "ni_P_image,ni_P_is_boundary,ni_elemface_idx, ni_pattern_resolution,"                    // &
            "ni_Mlink_u1,ni_Mlink_u2,ni_Mlink_u3,ni_Mlink_d1,ni_Mlink_d2,ni_Mlink_d3"

        write(5, '(A)')                                                                                 &
            "nYN_has_inflow,nYN_has_extInflow,nYN_has_dwfInflow"

        do ii = 1, size(node%R, 1)
            write(3,'(*(G0.6,:,","))') node%R(ii,:)
            write(4,'(*(G0.6,:,","))') node%I(ii,:)
            write(5,'(*(G0.6,:,","))') node%YN(ii,:)
        end do

        close(3)
        close(4)
        close(5)

    end subroutine util_output_export_linknode_input

    subroutine util_output_create_elemR_files

        integer :: fu, open_status, ii
        character(len = 250) :: file_name
        character(len = 40)  :: dir
        character(len = 5)   :: str_image
        character(len = 100) :: link_name
        character(len = 40)  :: str_elem_idx
        character(len = 10)  :: str_link_node_idx
        character(64) :: subroutine_name = 'util_output_create_elemR_files'

        if (setting%Debug%File%utility_output) print *, "*** enter ", this_image(), subroutine_name

        fu = this_image()

        write(str_image, '(i5.5)') fu

        call system("mkdir -p debug_output")

        do ii = 1, N_elem(this_image())

            write(str_elem_idx,'(I10)') elemI(ii,ei_Gidx)

            if (elemI(ii,ei_elementType) == CC) then

                write(str_link_node_idx,'(I10)') elemI(ii,ei_link_Gidx_SWMM)

                file_name = "debug_output/elemR/"//trim(str_image)//"_CC_" &
                    // trim(ADJUSTL(str_link_node_idx))// &
                    "_" // trim(ADJUSTL(str_elem_idx))//".csv"

                open(newunit=fu, file = file_name, status = 'replace',access = 'sequential', &
                form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                    stop "in " // subroutine_name
                end if

            else if (elemI(ii,ei_elementType) == JM) then

                write(str_link_node_idx,'(I10)') elemI(ii,ei_node_Gidx_SWMM)

                file_name = "debug_output/elemR/"//trim(str_image)//"_JM_" &
                    // trim(ADJUSTL(str_link_node_idx))// &
                    "_" // trim(ADJUSTL(str_elem_idx))//".csv"

                open(newunit=fu, file = file_name, status = 'replace',access = 'sequential', &
                form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                    stop "in " // subroutine_name
                end if

            else if (elemI(ii,ei_elementType) == JB) then

                 write(str_link_node_idx,'(I10)') elemI(ii,ei_node_Gidx_SWMM)

                 file_name = "debug_output/elemR/"//trim(str_image)//"_JB_" &
                     // trim(ADJUSTL(str_link_node_idx))// &
                     "_" // trim(ADJUSTL(str_elem_idx))//".csv"


                open(newunit=fu, file = file_name, status = 'replace',access = 'sequential', &
                    form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                    stop "in " // subroutine_name
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
        character(len = 5)   :: str_image
        character(len = 100) :: link_name
        character(len = 40)  :: str_face_idx
        character(64) :: subroutine_name = 'util_output_create_faceR_files'


        if (setting%Debug%File%utility_output) print *, "*** enter ", this_image(), subroutine_name

        fu = this_image()

        write(str_image, '(i5.5)') fu

        do ii = 1, N_face(this_image())

            write(str_face_idx,'(I10)') faceI(ii,fi_Gidx)

            file_name = "debug_output/faceR/"//trim(str_image)//"_face_" &
                // trim(ADJUSTL(str_face_idx))//".csv"

            open(newunit=fu, file = file_name, status = 'new',access = 'sequential', &
                form   = 'formatted', action = 'write', iostat = open_status)

            if (open_status /= 0) then
                write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                stop "in " // subroutine_name
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
        character(len = 5)   :: str_image
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

        write(str_image, '(i5.5)') fu

        do ii = 1, N_elem(this_image())

            write(str_elem_face_idx,'(I10)') elemI(ii,ei_Gidx)

            if (elemI(ii,ei_elementType) == CC) then

                write(str_link_node_idx,'(I10)') elemI(ii,ei_link_Gidx_SWMM)

                file_name = "debug_output/elemR/"//trim(str_image)//"_CC_" &
                    // trim(ADJUSTL(str_link_node_idx))// &
                    "_" // trim(ADJUSTL(str_elem_face_idx))//".csv"

                open(newunit=fu, file = file_name, status = 'old',access = 'append', &
                form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                    stop "in " // subroutine_name
                end if


            else if (elemI(ii,ei_elementType) == JM) then

                write(str_link_node_idx,'(I10)') elemI(ii,ei_node_Gidx_SWMM)

                file_name = "debug_output/elemR/"//trim(str_image)//"_JM_" &
                    // trim(ADJUSTL(str_link_node_idx))// &
                    "_" // trim(ADJUSTL(str_elem_face_idx))//".csv"


                open(newunit=fu, file = file_name, status = 'old',access = 'append', &
                form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                    stop "in " // subroutine_name
                end if


            else if (elemI(ii,ei_elementType) == JB) then

                write(str_link_node_idx,'(I10)') elemI(ii,ei_node_Gidx_SWMM)

                file_name = "debug_output/elemR/"//trim(str_image)//"_JB_" &
                    // trim(ADJUSTL(str_link_node_idx))// &
                    "_" // trim(ADJUSTL(str_elem_face_idx))//".csv"


                open(newunit=fu, file = file_name, status = 'old',access = 'append', &
                    form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                    stop "in " // subroutine_name
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

            file_name = "debug_output/faceR/"//trim(str_image)//"_face_" &
                // trim(ADJUSTL(str_elem_face_idx))//".csv"

            open(newunit=fu, file = file_name, status = 'old',access = 'Append', &
                form   = 'formatted', action = 'write', iostat = open_status)

            if (open_status /= 0) then
                write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                stop "in " // subroutine_name
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

    subroutine util_output_create_summary_files
        integer :: fu, open_status
        character(64) :: file_name
        character(64) :: subroutine_name = 'util_output_create_summary_files'

        if (setting%Debug%File%utility_output) print *, "*** enter ", this_image(), subroutine_name

        write(file_name, "(A,i5.5,A)") "debug_output/summary/summary_", this_image(), ".csv"

        open(newunit=fu, file = file_name, status = 'replace',access = 'sequential', &
            form = 'formatted', action = 'write', iostat = open_status)

        write(fu, *) "In_Image,This_Time,CFL_max,dt,Velocity_Max,Wavespeed_Max"
        endfile(fu)
        close(fu)

        if (setting%Debug%File%utility_output) print *, "*** leave ", this_image(), subroutine_name
    end subroutine util_output_create_summary_files

    subroutine util_output_report
        character(64) :: subroutine_name = "util_output_report"

        if (setting%Debug%File%utility_output) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        if (setting%Debug%Output) call util_output_report_summary()

        if (setting%Output%report .and. util_output_must_report()) then
            if (setting%Debug%Output) call util_output_write_elemR_faceR()
            call output_write_link_files()
            call output_write_node_files()
        end if

        if (setting%Debug%File%utility_output) &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine util_output_report

    subroutine util_output_report_summary()
        integer          :: fu, open_status, thisCol, Npack
        integer, pointer :: thisP(:)
        real(8)          :: thisCFL, max_velocity, max_wavespeed
        real(8), pointer :: dt, timeNow, velocity(:), wavespeed(:), length(:)
        character(512)    :: file_name
        character(64)    :: subroutine_name = "util_output_report_summary"

        if (setting%Debug%File%utility_output) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"
        if (util_output_must_report() .and. setting%output%report) then

            write(file_name, "(A,i5.5,A)") "debug_output/summary/summary_", this_image(), ".csv"

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

        if (setting%Debug%File%utility_output) &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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

end module utility_output
