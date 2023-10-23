module common_elements
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Common procedures for diagnostic elements
    !%==========================================================================
    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use utility_crash, only: util_crashpoint

    implicit none

    private

    public :: common_velocity_from_flowrate_singular
    public :: common_head_and_flowdirection_singular
    public :: common_outflow_energyhead_singular

    contains
!%
!%========================================================================== 
! PUBLIC
!%==========================================================================    
!%  
    subroutine common_velocity_from_flowrate_singular (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes velocity from flowrate for a single element
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx
            real(8), pointer :: Flowrate, Area, Velocity, Vmax
        !%------------------------------------------------------------------
        !% Aliases
            Vmax     => setting%Limiter%Velocity%Maximum
            Velocity => elemR(eIdx,er_Velocity)
            Flowrate => elemR(eIdx,er_Flowrate)
            Area     => elemR(eIdx,er_Area)
        !%------------------------------------------------------------------
        Velocity = Flowrate / Area
        
        !% Velocity limiter
        if (setting%Limiter%Velocity%UseLimitMaxYN) then
            if (abs(Velocity) > Vmax) then
                Velocity = sign(0.99d0*Vmax,Velocity)
            end if
        end if
        
    end subroutine common_velocity_from_flowrate_singular
!%      
!%==========================================================================
!%==========================================================================    
!%  
    subroutine common_head_and_flowdirection_singular &
        (eIdx, ZcrestCol, NominalDownstreamHeadCol, FlowDirectionCol)
        !%------------------------------------------------------------------
        !% Description:
        !% computes and stores the maximum head, flow direction, and downstream head
        !% on element. Designed for Weirs and Orifices for a single element
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx  !% must be a single element ID
            !% --- input is (e.g.) esr_Weir_Zcrest, esr_Weir_NominalDownstreamHead, esi_Weir_FlowDirection
            integer, intent(in) :: ZcrestCol, NominalDownstreamHeadCol, FlowDirectionCol
            real(8), pointer :: Head, NominalDSHead
            real(8), pointer :: UpstreamFaceHead, DownstreamFaceHead, Zcrest
            integer, pointer :: FlowDirection, iupf, idnf
        !%------------------------------------------------------------------
        !% Aliases
            !% --- outputs
            Head           => elemR(eIdx,er_Head)
            NominalDSHead  => elemSR(eIdx,NominalDownstreamHeadCol)
            FlowDirection  => elemSI(eIdx,FlowDirectionCol)
            !% --- element data used
            Zcrest => elemSR(eIdx,ZcrestCol)
            !% --- face locations
            iupf    => elemI(eIdx,ei_Mface_uL)
            idnf    => elemI(eIdx,ei_Mface_dL)
            !% --- face data used
            UpstreamFaceHead   => faceR(iupf,fr_Head_d)
            DownstreamFaceHead => faceR(idnf,fr_Head_u)
        !%------------------------------------------------------------------        
        !% head on a diagnostic element as the maximum of upstream, downstream, or crest height.
        Head = max(UpstreamFaceHead , DownstreamFaceHead)
        
        !% flow direction on a diagnostic element assigned based up upstream and downstream heads
        FlowDirection = int(sign(oneR, (UpstreamFaceHead - DownstreamFaceHead)))
        
        !% nominal downstream head on a diagnostic element
        NominalDSHead = min(UpstreamFaceHead, DownstreamFaceHead)
        
    end subroutine common_head_and_flowdirection_singular
!%      
!%==========================================================================
!%==========================================================================    
!%  
    subroutine common_outflow_energyhead_singular &
         (eIdx, NominalDownstreamHeadCol, FlowDirectionCol)
        !%------------------------------------------------------------------
        !% Description:
        !% computes and stores the energyhead available at the outflow
        !% of a diagnostic element (weir,pump,orifice)
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx  !% must be a single element ID
            integer, intent(in) :: NominalDownstreamHeadCol, FlowDirectionCol 
            integer, pointer :: fadj
        !%------------------------------------------------------------------
        !% Aliases
   
        !%------------------------------------------------------------------
        if (elemSI(eIdx,FlowDirectionCol) == oneI) then 
            !% --- downstream flow
            !%     pointer to downstream face
            fadj => elemI(eIdx,ei_Mface_dL)
            !% --- energy using upstream face values
            elemR(eIdx,er_EnergyHead) = elemSR(eIdx,NominalDownstreamHeadCol) &
                 + (faceR(fadj,fr_Velocity_u)**2) / (twoR * setting%Constant%gravity)

        elseif (elemSI(eIdx,FlowDirectionCol) == -oneI) then 
            !% --- upstream flow
            !%     pointer to upstream face
            fadj => elemI(eIdx,ei_Mface_uL)
            !% --- energy head using downstream face values (note that NominalDownstreamHead)
            !%     should contain the value at the upstream face.
            elemR(eIdx,er_EnergyHead) = elemSR(eIdx,NominalDownstreamHeadCol) &
                 + (faceR(fadj,fr_Velocity_d)**2) / (twoR * setting%Constant%gravity)
        else
            print *, 'Unexpected else'
            call util_crashpoint(710983)
        end if    
            
    end subroutine common_outflow_energyhead_singular    
!%
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module common_elements