module common_elements

    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Provides common computations useful over multiple element types
    !%
   
    private

    public :: common_velocity_from_flowrate_singular
    public :: common_head_and_flowdirection_singular


    contains
!%
!%========================================================================== 
! PUBLIC
!%==========================================================================    
!%  
    subroutine common_velocity_from_flowrate_singular (eIdx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx
        real(8), pointer :: Flowrate, Area, Velocity, Vmax
        !%-----------------------------------------------------------------------------
        if (crashYN) return
        Vmax     => setting%Limiter%Velocity%Maximum
        Velocity => elemR(eIdx,er_Velocity)
        Flowrate => elemR(eIdx,er_Flowrate)
        Area     => elemR(eIdx,er_Area)
        !%-----------------------------------------------------------------------------

        Velocity = Flowrate / Area
        
        !% Velocity limiter
        if (setting%Limiter%Velocity%UseLimitMaxYN) then
            if (abs(Velocity) > Vmax) then
                Velocity = sign(0.99*Vmax,Velocity)
            end if
        end if
        
    end subroutine common_velocity_from_flowrate_singular
!%      
!%==========================================================================
!%==========================================================================    
!%  
    subroutine common_head_and_flowdirection_singular &
        (eIdx, ZcrestCol, NominalDownstreamHeadCol, FlowDirectionCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes and stores the maximum head, flow direction, and downstream head
        !% on element. Designed for Weirs and Orifices
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx  !% must be a single element ID
        !% input is (e.g.) esr_Weir_Zcrest, esr_Weir_NominalDownstreamHead, esi_Weir_FlowDirection
        integer, intent(in) :: ZcrestCol, NominalDownstreamHeadCol, FlowDirectionCol
        real(8), pointer :: Head, NominalDSHead
        real(8), pointer :: UpstreamFaceHead, DownstreamFaceHead, Zcrest
        integer, pointer :: FlowDirection, iupf, idnf
        !%-----------------------------------------------------------------------------
        if (crashYN) return
        !% outputs
        Head           => elemR(eIdx,er_Head)
        NominalDSHead  => elemSR(eIdx,NominalDownstreamHeadCol)
        FlowDirection  => elemSI(eIdx,FlowDirectionCol)
        !% element data used
        Zcrest => elemSR(eIdx,ZcrestCol)
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
        
    end subroutine common_head_and_flowdirection_singular
!%
!%==========================================================================
!%==========================================================================    
!%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%  
!%
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module common_elements