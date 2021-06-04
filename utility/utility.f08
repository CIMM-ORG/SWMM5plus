module utility

    use define_indexes
    use define_keys
    use define_globals
    use define_settings, only: setting

    implicit none

!-----------------------------------------------------------------------------
!
! Description:
!   Utility routines that may be called in a number of places
!
!-----------------------------------------------------------------------------

    private

    public :: util_export_linknode_csv
    public :: util_count_node_types

contains
    !
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    !
    subroutine util_export_linknode_csv()
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

        write(4, '(A)')                                                                                 &
            "ni_idx,ni_node_type,ni_N_link_u,ni_N_link_d,ni_curve_type,ni_assigned,ni_total_inflow," // &
            "ni_Mlink_u1,ni_Mlink_u2,ni_Mlink_u3,ni_Mlink_d1,ni_Mlink_d2,ni_Mlink_d3"

        do i = 1, N_node
            write(3,'(*(G0.6,:,","))') nodeR(i,:)
            write(4,'(*(G0.6,:,","))') nodeI(i,:)
        end do

        close(3)
        close(4)

    end subroutine util_export_linknode_csv
    !
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    !
    subroutine util_count_node_types(N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2)
        integer, intent(in out) :: N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2
        integer :: ii
    
        ! This subroutine uses the vectorized count() function to search the array for number of instances of each node type
        N_nBCup = count(nodeI(:, ni_node_type) == nBCup)
        N_nBCdn = count(nodeI(:, ni_node_type) == nBCdn)
        N_nJm = count(nodeI(:, ni_node_type) == nJM)
        N_nStorage = count(nodeI(:, ni_node_type) == nStorage)
        N_nJ2 = count(nodeI(:, ni_node_type) == nJ2)
    
    end subroutine util_count_node_types
    
end module utility