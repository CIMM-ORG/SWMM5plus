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
    public :: pump_set_setting

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
        !% Preliminaries
            if (setting%Debug%File%pump_elements) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------    

        !% --- set pump setting

        call pump_flow (eIdx)

        !%-----------------------------------------------------------------------------
        !% Closing:
            if (setting%Debug%File%pump_elements)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine pump_toplevel
!%
!%========================================================================== 
!%==========================================================================    
!%     
    subroutine pump_set_setting (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% Updates the setting on a pump element, following EPA-SWMM
        !% link.c/link_setSetting where pump setting is immediately set to 
        !% targetsetting after a control change
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIdx
        !%------------------------------------------------------------------   

        elemR(eIdx,er_Setting) = elemR(eIdx,er_TargetSetting)

        !% --- error check
        !%     EPA-SWMM allows the pump setting to be only 0.0 or 1.0
        if (.not. ((elemR(eIdx,er_Setting) == 0.0) .or. (elemR(eIdx,er_Setting) == 1.0))) then
            print *, 'CODE ERROR: pump element has er_Setting that is not 0.0 or 1.0'
            call util_crashpoint(668723)
        end if

    end subroutine pump_set_setting
!%
!%========================================================================== 
!% PRIVATE
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



        case (type4_Pump)
            call pump_type4(eIdx)
               

        case (type_IdealPump)  
            call pump_ideal(eIdx)

        end select

        !% reverse flow through pump are not allowed
        FlowRate = FlowRate * PSetting
        if (FlowRate < zeroR) FlowRate = zeroR 
  
    end subroutine pump_flow
!%
!%========================================================================== 
!% PRIVATE
!%==========================================================================    
!%  
    subroutine pump_type1 (idx)
        !%------------------------------------------------------------------
        !% Description:
        !% EPA-SWMM type 1 pump that has constant flow rates for volume
        !% intervals
        !%------------------------------------------------------------------
        !% Declarations:
        integer, intent(in) :: idx  !% pump index
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------

        !%------------------------------------------------------------------
        !% Closing:
    end subroutine pump_type1
!%
!%==========================================================================   
!%==========================================================================    
!%  
    subroutine pump_type2 (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% 
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx  !% pump index
            integer, pointer :: PumpType, iupf, idnf, CurveID
            real(8), pointer :: FlowRate, UpFaceHead
            real(8), pointer :: Depth, DnFaceZbottom
        !%------------------------------------------------------------------
        !% Aliases:
            CurveID  => elemSI(eIdx,esi_Pump_CurveID)
            Depth    => elemR(eIdx,er_Depth)
            FlowRate => elemR(eIdx,er_Flowrate)
            !% face locations
            iupf     => elemI(eIdx,ei_Mface_uL)
            idnf     => elemI(eIdx,ei_Mface_dL)
            !% face data used
            UpFaceHead => faceR(iupf,fr_Head_d)
            DnFaceZbottom => faceR(idnf,fr_Zbottom)  
        !%------------------------------------------------------------------
        !% --- set the depth upstream that is the Type2 pump control point    
        Depth = UpFaceHead - DnFaceZbottom

        !% --- use the curve without interpolation
        call util_curve_lookup_singular( &
            CurveID, er_Depth, er_Flowrate, curve_pump_Xvar, curve_pump_flowrate,0)

    end subroutine pump_type2
!%
!%==========================================================================
!%==========================================================================    
!%
    subroutine pump_type3 (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% 
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: PumpType, iupf, idnf, CurveID
            real(8), pointer :: Head, FlowRate, UpFaceHead, DnFaceHead
            integer, intent(in) :: eIdx  !% pump index
        !%------------------------------------------------------------------
        !% Aliases:
            CurveID  => elemSI(eIdx,esi_Pump_CurveID)
            Head     => elemR(eIdx,er_Head)
            FlowRate => elemR(eIdx,er_Flowrate)
            !% face locations
            iupf     => elemI(eIdx,ei_Mface_uL)
            idnf     => elemI(eIdx,ei_Mface_dL)
            !% face data used
            UpFaceHead => faceR(iupf,fr_Head_d)
            DnFaceHead => faceR(idnf,fr_Head_u)
        !%------------------------------------------------------------------

        !% --- temporarily set the head to the head difference across the 
        !%     pump that controls a Type 3 (centrifugal) pump   
        Head = max((DnFaceHead - UpFaceHead), zeroR)

        !% --- use the curve for this pump and the head difference
        call util_curve_lookup_singular( &
            CurveID, er_Head, er_Flowrate, curve_pump_Xvar, curve_pump_flowrate,1)

        !% --- reset the pump element head to the upstream face value
        Head = UpFaceHead   

    end subroutine pump_type3
!%
!%==========================================================================
!%========================================================================== 
!%  
    subroutine pump_type4 (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% 
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx  !% pump index
            integer, pointer :: PumpType, iupf, idnf, CurveID
            real(8), pointer :: FlowRate, UpFaceHead
            real(8), pointer :: Depth, DnFaceZbottom
        !%------------------------------------------------------------------
        !% Aliases:
            CurveID  => elemSI(eIdx,esi_Pump_CurveID)
            Depth    => elemR(eIdx,er_Depth)
            FlowRate => elemR(eIdx,er_Flowrate)
            !% face locations
            iupf     => elemI(eIdx,ei_Mface_uL)
            idnf     => elemI(eIdx,ei_Mface_dL)
            !% face data used
            UpFaceHead => faceR(iupf,fr_Head_d)
            DnFaceZbottom => faceR(idnf,fr_Zbottom)  
        !%------------------------------------------------------------------

        !% --- set the depth upstream that is the Type4 pump control point    
        Depth = UpFaceHead - DnFaceZbottom

        !% --- use the curve for this pump to set the flowrate for this pump
        !%     based on the depth
        call util_curve_lookup_singular ( &
            CurveID, er_Depth, er_Flowrate, curve_pump_Xvar, curve_pump_flowrate,1)

        !%------------------------------------------------------------------
    end subroutine pump_type4
!%
!%==========================================================================
!%==========================================================================    
!%  
    subroutine pump_ideal (idx)
        !%------------------------------------------------------------------
        !% Description:
        !% EPA-SWMM Ideal Pump -- which simply repeats its inlet condition
        !% as a flow condition
        !%
        !% HACK---EPA-SWMM uses the upstream node flow plus any overflow rate
        !%        In the following we have not included any overflow as we
        !%        use the upstream face flow value
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: idx  !% pump index
            integer, pointer    :: fup
            real(8), pointer    :: flowrate, face_flowrate(:)
        !%------------------------------------------------------------------
        !% Aliases:
            fup           => elemI(idx,ei_Mface_uL)
            flowrate      => elemR(idx,er_Flowrate)
            face_flowrate => faceR(:  ,fr_Flowrate)
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
            flowrate = face_flowrate(fup)
        !%------------------------------------------------------------------
        !% Closing:
    end subroutine pump_ideal
!%
!%==========================================================================
!%==========================================================================   
!%  
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