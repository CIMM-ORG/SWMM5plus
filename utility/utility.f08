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
    public :: util_sign_with_ones
    public :: util_head_and_flowdirection_for_diagnostic_element
    public :: util_velocity_from_flowrate_singular

    contains
    !%
    !%==========================================================================  
    !% PUBLIC
    !%==========================================================================  
    !%
    subroutine util_export_linknode_csv()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Exports link and node tables as CSV files
        !%-----------------------------------------------------------------------------
        integer :: ii
        logical :: ex
        !%-----------------------------------------------------------------------------
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
            "li_Melem_d,li_Mface_u,li_Mface_d,li_assigned,li_InitialDepthType," // &
            "li_length_adjusted,li_P_image,li_first_elem_idx,li_last_elem_idx"

        do ii = 1, N_link
            write(1,'(*(G0.6,:,","))') linkR(ii,:)
            write(2,'(*(G0.6,:,","))') linkI(ii,:)
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
            "ni_P_image,ni_P_is_boundary,ni_Mlink_u1,ni_Mlink_u2,ni_Mlink_u3,ni_Mlink_d1,"           // &
            "ni_Mlink_d2,ni_Mlink_d3"

        do ii = 1, N_node
            write(3,'(*(G0.6,:,","))') nodeR(ii,:)
            write(4,'(*(G0.6,:,","))') nodeI(ii,:)
        end do

        close(3)
        close(4)

    end subroutine util_export_linknode_csv
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine util_count_node_types(N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% This subroutine uses the vectorized count() function to search the array for 
        !% numberof instances of each node type
        !%-----------------------------------------------------------------------------
        integer, intent(in out) :: N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2
        integer :: ii
        !%-----------------------------------------------------------------------------
        N_nBCup = count(nodeI(:, ni_node_type) == nBCup)
        N_nBCdn = count(nodeI(:, ni_node_type) == nBCdn)
        N_nJm = count(nodeI(:, ni_node_type) == nJM)
        N_nStorage = count(nodeI(:, ni_node_type) == nStorage)
        N_nJ2 = count(nodeI(:, ni_node_type) == nJ2)
    
    end subroutine util_count_node_types
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    pure elemental real(8) function util_sign_with_ones &
        (inarray) result (outarray)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% returns is an array of real ones with the sign of the inarray argument
        !%-----------------------------------------------------------------------------
        real(8),      intent(in)    :: inarray
        !%-----------------------------------------------------------------------------
        outarray = oneR
        outarray = sign(outarray,inarray)

    end function util_sign_with_ones  
    !%========================================================================== 
    !%==========================================================================    
    !%  
    subroutine util_head_and_flowdirection_for_diagnostic_element (eIdx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes the maximum head, flow direction, and downstream head
        !% on any single diagnostic element. 
        !%
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx  !% must be a single element ID
        real(8), pointer :: Head, NominalDSHead
        real(8), pointer :: UpstreamFaceHead, DownstreamFaceHead, Zcrest
        integer, pointer :: FlowDirection, iupf, idnf
        !%-----------------------------------------------------------------------------
        !% outputs
        Head           => elemR(eIdx,er_Head)
        NominalDSHead  => elemSR(eIdx,eSr_Weir_NominalDownstreamHead)
        FlowDirection  => elemSI(eIdx,eSi_Weir_FlowDirection)
        !% element data used
        Zcrest => elemSR(eIdx,eSr_Weir_Zcrest)
        !% face locations
        iupf    => elemI(eIdx,ei_Mface_uL)
        idnf    => elemI(eIdx,ei_Mface_dL)
        ! face data used
        UpstreamFaceHead   => faceR(iupf,fr_Head_d)
        DownstreamFaceHead => faceR(idnf,fr_Head_u)
        !%-----------------------------------------------------------------------------        
        !% head on a diagnostic element as the maximum of upstream, downstream, or crest height.
        Head = max(UpstreamFaceHead , DownstreamFaceHead , Zcrest)
        
        !% flow direction on a diagnostic element assigned based up upstream and downstream heads
        FlowDirection = int(sign(oneR, (UpstreamFaceHead - DownstreamFaceHead)) )
        
        !% nominal downstream head on a diagnostic element
        NominalDSHead = min(UpstreamFaceHead, DownstreamFaceHead)
        
    end subroutine util_head_and_flowdirection_for_diagnostic_element
    !%
    !%========================================================================== 
    !%==========================================================================    
    !%  
    subroutine util_velocity_from_flowrate_singular (eIdx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx
        real(8), pointer :: Flowrate, Area, Velocity, Vmax
        logical, pointer :: isAdHocFlowrate
        !%-----------------------------------------------------------------------------
        Vmax     => setting%Limiter%Velocity%Maximum
        Velocity => elemR(eIdx,er_Velocity)
        Flowrate => elemR(eIdx,er_Flowrate)
        Area     => elemR(eIdx,er_Area)
        isAdHocFlowrate => elemYN(eIdx,eYN_IsAdHocFlowrate)
        !%-----------------------------------------------------------------------------

        Velocity = Flowrate / Area
        
        !% Velocity limiter
        elemYN(eIdx,eYN_IsAdHocFlowrate) = .true.
        if (setting%Limiter%Velocity%UseLimitMax) then
            if (abs(Velocity) > Vmax) then
                Velocity = sign(0.99*Vmax,Velocity)
                isAdHocFlowrate = .true.
            else
                isAdHocFlowrate = .false.
            endif
        endif
        
    end subroutine util_velocity_from_flowrate_singular
    !%      
    !%==========================================================================
    !% END OF MODULE
    !%+=========================================================================    
end module utility