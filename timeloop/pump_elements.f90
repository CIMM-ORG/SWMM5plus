module pump_elements

    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use common_elements
    use adjust
    use utility_interpolate
    use utility_crash, only: util_crashpoint


    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Computes diagnostic flow through orifice elements
    !%
    !% METHOD:
    !% 
    !%

    private

    public :: pump_toplevel

    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine pump_toplevel (eIdx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% We need a subroutine to get the new full depth (esr_EffectiveFullDepth)
        !% crown (esr_Zcrown) and crest (esr_Zcrest) elevation from control setting.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx  !% must be a single element ID

        character(64) :: subroutine_name = 'pump_toplevel'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%pump_elements) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

            call pump_flow (eIdx)

        if (setting%Debug%File%pump_elements)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine pump_toplevel
!%
!%========================================================================== 
!%==========================================================================    
!%  
    subroutine pump_flow (eIdx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx !% must be single element ID
        integer, pointer :: PumpType, iupf, idnf, CurveID
        real(8), pointer :: Head, FlowRate, UpFaceHead, DnFaceHead, PSetting
        real(8), pointer :: Depth, DnFaceZbottom
        !%-----------------------------------------------------------------------------
        ! if (crashYN) return
        PumpType => elemSI(eIdx,esi_Pump_SpecificType)
        CurveID  => elemSI(eIdx,esi_Pump_CurveID)
        Head     => elemR(eIdx,er_Head)
        Depth    => elemR(eIdx,er_Depth)
        FlowRate => elemR(eIdx,er_Flowrate)
        PSetting => elemR(eIdx,er_Setting)
        !% face locations
        iupf     => elemI(eIdx,ei_Mface_uL)
        idnf     => elemI(eIdx,ei_Mface_dL)
        !% face data used
        UpFaceHead => faceR(iupf,fr_Head_d)
        DnFaceHead => faceR(idnf,fr_Head_u)
        DnFaceZbottom => faceR(idnf,fr_Zbottom)

        !%-----------------------------------------------------------------------------
        select case (PumpType)
        case (type1_Pump)

                print *, 'CODE ERROR: ',trim(reverseKey(PumpType)), ' has not yet been developed'
                call util_crashpoint(561893)
                return

        case (type2_Pump)

                print *, 'CODE ERROR: ',trim(reverseKey(PumpType)), ' has not yet been developed'
                call util_crashpoint(23489)
                return

        case (type3_Pump)

                Head = max((DnFaceHead - UpFaceHead), zeroR)
                call util_curve_lookup_singular(CurveID, er_Head, er_Flowrate, &
                        curve_pump_Xvar, curve_pump_flowrate)

        case (type4_Pump)

                Depth = UpFaceHead - DnFaceZbottom
                call util_curve_lookup_singular(CurveID, er_Depth, er_Flowrate, &
                        curve_pump_Xvar, curve_pump_flowrate)

        case (type_IdealPump)  
                
                print *, 'CODE ERROR: ',trim(reverseKey(PumpType)), ' has not yet been developed'
                call util_crashpoint(784666)
                return

        end select

        !% reverse flow through pump are not allowed
        FlowRate = FlowRate * PSetting
        if (FlowRate < zeroR) FlowRate = zeroR 
  
    end subroutine pump_flow
! %
!%==========================================================================
!% PRIVATE
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
!%==========================================================================
end module pump_elements