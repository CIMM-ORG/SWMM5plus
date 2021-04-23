module utility

    use array_index
    use data_keys
    use setting_definition
    use globals

    implicit none

!***********************************************************************************
!
! Description:
!   Utility routines that may be called in a number of places
!
!***********************************************************************************

    private

    integer :: debuglevel = 0

    public :: utility_check_allocation
    public :: utility_export_linknode_csv

contains

    subroutine utility_check_allocation(allocation_status, emsg)
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   Checks allocation status and STOPs if there is an error
    !
    !-----------------------------------------------------------------------------

        integer,           intent(in   ) :: allocation_status
        character(len=*),  intent(in   ) :: emsg

        character(64):: subroutine_name = 'utility_check_allocation'

    !-----------------------------------------------------------------------------

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        if (allocation_status > 0) then
            print *, trim(emsg)
            STOP
        end if

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

    end subroutine utility_check_allocation

    subroutine utility_export_linknode_csv()
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   Exports link and node tables as CSV files
    !
    !-----------------------------------------------------------------------------

        integer :: i
        logical :: ex

    !-----------------------------------------------------------------------------

        open(unit=1,file='debug/linkR.csv',status='new')
        open(unit=2,file='debug/linkI.csv',status='new')

        write(1, '(A)')                                                                    &
            "lr_Length,lr_InletOffset,lr_OutletOffset,lr_BreadthScale,"                 // &
            "lr_TopWidth,lr_ElementLength,lr_Slope,lr_LeftSlope,lr_RightSlope,"         // &
            "lr_Roughness,lr_InitialFlowrate,lr_InitialDepth,lr_InitialUpstreamDepth,"  // &
            "lr_InitialDnstreamDepth,lr_ParabolaValue,lr_SideSlope,lr_DischargeCoeff1," // &
            "lr_DischargeCoeff2,lr_FullDepth,lr_EndContractions,lr_Flowrate,lr_Depth,"  // &
            "lr_DepthUp,lr_DepthDn,lr_Volume,lr_Velocity,lr_Capacity"

        write(2, '(A)')                                                            &
            "li_idx,li_link_type,li_weir_type,li_orif_type,"                    // &
            "li_pump_type,li_geometry,li_roughness_type,"                       // &
            "li_N_element,li_Mnode_u,li_Mnode_d,li_Melem_u,"                    // &
            "li_Melem_d,li_Mface_u,li_Mface_d,li_assigned,li_InitialDepthType"

        do i = 1, N_link
            write(1,'(*(G0.6,:,","))') linkR(i,:)
            write(2,'(*(G0.6,:,","))') linkI(i,:)
        end do

        close(1)
        close(2)

        open(unit=3,file='debug/nodeR.csv',status='new')
        open(unit=4,file='debug/nodeI.csv',status='new')

        write(3, '(A)')                                                                              &
            "nr_Zbottom,nr_InitialDepth,nr_FullDepth,nr_StorageConstant,nr_StorageCoeff,"         // &
            "nr_StorageExponent,nr_PondedArea,nr_SurchargeDepth,nr_MaxInflow,nr_Eta,"             // &
            "nr_Depth,nr_Volume,nr_LateralInflow,nr_TotalInflow,nr_Flooding,nr_ElementLength_u1," // &
            "nr_ElementLength_u2,nr_ElementLength_u3,nr_ElementLength_d1,nr_ElementLength_d2,"    // &
            "nr_ElementLength_d3"

        write(4, '(A)')
            "ni_idx,ni_node_type,ni_N_link_u,ni_N_link_d,ni_curve_type,ni_assigned,ni_total_inflow," // &
            "ni_Mlink_u1,ni_Mlink_u2,ni_Mlink_u3,ni_Mlink_d1,ni_Mlink_d2,ni_Mlink_d3"

        do i = 1, N_node
            write(3,'(*(G0.6,:,","))') nodeR(i,:)
            write(4,'(*(G0.6,:,","))') nodeI(i,:)
        end do

        close(3)
        close(4)

    end subroutine utility_export_linknode_csv

end module utility