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
        !%------------------------------------------------------------------
        !% Description:
        !% We need a subroutine to get the new full depth (esr_EffectiveFullDepth)
        !% crown (esr_Zcrown) and crest (esr_Zcrest) elevation from control setting.
        !%------------------------------------------------------------------
            integer, intent(in) :: eIdx  !% must be a single element ID
            integer, pointer :: PumpType
            real(8), pointer :: FlowRate,  PSetting
            character(64) :: subroutine_name = 'pump_toplevel'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%pump_elements) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------  
        !% Aliases
            PumpType => elemSI(eIdx,esi_Pump_SpecificType)
            FlowRate => elemR(eIdx,er_Flowrate)
            PSetting => elemR(eIdx,er_Setting)
        !%-------------------------------------------------------------------
        select case (PumpType)

            case (type1_Pump)
                call pump_type1(eIdx)

            case (type2_Pump)
                call pump_type2(eIdx)
            
            case (type3_Pump)
                call pump_type3(eIdx)

            case (type4_Pump)
                call pump_type4(eIdx)
                    
            case (type_IdealPump)  
                call pump_ideal(eIdx)

        end select
    
        !% --- linearly reduce flow by the setting value
        FlowRate = FlowRate * PSetting

        !% --- prohibit reverse flow through pump
        if (FlowRate < zeroR) FlowRate = zeroR 

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
    subroutine pump_type1 (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% EPA-SWMM type 1 pump whos flow is determined by volume up
        !% upstream element. This volume must be copied to the volume
        !% storage for the pump element
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx  !% pump index
            integer, pointer :: iupf, eup, CurveID, JMidx
            real(8), pointer :: FlowRate, Volume
            real(8), pointer :: Head, UpFaceHead
            logical, pointer :: isGhostUp
        !%------------------------------------------------------------------
        !% Aliases:
            CurveID  => elemSI(eIdx,esi_Pump_CurveID)
            FlowRate => elemR(eIdx,er_Flowrate)
            Head     => elemR(eIdx,er_Head)
            Volume   => elemR(eIdx,er_Volume)
            !% face locations
            iupf     => elemI(eIdx,ei_Mface_uL)
            !% face data used
            UpFaceHead    => faceR(iupf,fr_Head_d)
            eup           => faceI(iupf,fi_Melem_uL)
            isGhostUp     => faceYN(iupf,fYN_isUpGhost)
        !%------------------------------------------------------------------

        !% --- store volume of upstream element at pump
        !%     this volume isn't considered in volume conservation
        if (isGhostUp) then
            !% --- handle condition where upstream is on another image
            Volume = elemGR(eup,er_Volume)
        else
            !% --- check for a junction branch element upstream
            if (elemI(eUp,ei_elementType) == JB) then
                !% --- if JB, use the associated JM volume
                !JMidx => elemI(eup,ei_main_idx_for_branch)
                JMidx => elemSI(eup,esi_JunctionBranch_Main_Index)
                Volume = elemR(JMidx,er_Volume)
            else
                !% --- for any element other than JB
                Volume = elemR(eup,er_Volume)
            end if
        end if

        if (Volume .le. setting%ZeroValue%Volume) then
            !% --- no flow for effectively zero volume
            Flowrate = zeroR
        else
            !% --- use the curve without interpolation
            call util_curve_lookup_singular( &
                CurveID, er_Volume, er_Flowrate, curve_pump_Xvar, curve_pump_flowrate,0)
        end if

        !% --- set the head of the pump to the upstream face head
        Head = UpFaceHead

    end subroutine pump_type1
!%
!%==========================================================================   
!%==========================================================================    
!%  
    subroutine pump_type2 (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes flowrate for pump where flow is set by depth at upstream
        !% face
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx  !% pump index
            integer, pointer :: iupf, CurveID
            real(8), pointer :: FlowRate, UpFaceHead
            real(8), pointer :: Depth, Head, UpFaceZbottom
        !%------------------------------------------------------------------
        !% Aliases:
            CurveID  => elemSI(eIdx,esi_Pump_CurveID)
            Depth    => elemR(eIdx,er_Depth)
            FlowRate => elemR(eIdx,er_Flowrate)
            Head     => elemR(eIdx,er_Head)
            !% face locations
            iupf     => elemI(eIdx,ei_Mface_uL)
            !% face data used
            UpFaceHead    => faceR(iupf,fr_Head_d)
            UpFaceZbottom => faceR(iupf,fr_Zbottom)  
        !%------------------------------------------------------------------
        !% --- set the depth upstream that is the Type2 pump control point    
        Depth = UpFaceHead - UpFaceZbottom

        if (Depth .le. setting%ZeroValue%Depth) then
            !% --- no flow for effectively zero depth
            Flowrate = zeroR
        else
            !% --- use the curve without interpolation
            call util_curve_lookup_singular( &
                CurveID, er_Depth, er_Flowrate, curve_pump_Xvar, curve_pump_flowrate,0)
        end if

        !% --- set the head of the pump to the upstream face head
        Head = UpFaceHead


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
            integer, pointer :: iupf, idnf, CurveID
            real(8), pointer :: Head, Depth, FlowRate, UpFaceHead, DnFaceHead
            real(8), pointer :: UpFaceZbottom
            integer, intent(in) :: eIdx  !% pump index
        !%------------------------------------------------------------------
        !% Aliases:
            CurveID  => elemSI(eIdx,esi_Pump_CurveID)
            Head     => elemR(eIdx,er_Head)
            Depth    => elemR(eIdx,er_Depth)
            FlowRate => elemR(eIdx,er_Flowrate)
            !% face locations
            iupf     => elemI(eIdx,ei_Mface_uL)
            idnf     => elemI(eIdx,ei_Mface_dL)
            !% face data used
            UpFaceHead => faceR(iupf,fr_Head_d)
            DnFaceHead => faceR(idnf,fr_Head_u)
            UpFaceZbottom => faceR(iupf,fr_Zbottom)
        !%------------------------------------------------------------------

        !% --- temporarily set the head to the head difference across the 
        !%     pump that controls a Type 3 (centrifugal) pump   
        Head = max((DnFaceHead - UpFaceHead), zeroR)

        !% --- get the depth upstream
        Depth = UpFaceHead - UpFaceZbottom

        if ((Depth .le. setting%ZeroValue%Depth) .or. (Head .le. zeroR)) then
            !% --- no flow for effectively zero depth or no head difference
            Flowrate = zeroR
        else
            !% --- use the curve for this pump and the head difference
            call util_curve_lookup_singular( &
                CurveID, er_Head, er_Flowrate, curve_pump_Xvar, curve_pump_flowrate,1)
        end if

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
            integer, pointer :: iupf,  CurveID
            real(8), pointer :: FlowRate, UpFaceHead
            real(8), pointer :: Depth, Head, UpFaceZbottom
        !%------------------------------------------------------------------
        !% Aliases:
            CurveID  => elemSI(eIdx,esi_Pump_CurveID)
            Depth    => elemR(eIdx,er_Depth)
            Flowrate => elemR(eIdx,er_Flowrate)
            Head     => elemR(eIdx,er_Head)
            !% face locations
            iupf     => elemI(eIdx,ei_Mface_uL)
            !% face data used
            UpFaceHead    => faceR(iupf,fr_Head_d)
            UpFaceZbottom => faceR(iupf,fr_Zbottom)     
        !%------------------------------------------------------------------
        !% --- set the depth upstream that is the Type4 pump control point    
        Depth = UpFaceHead - UpFaceZbottom

        if (Depth .le. setting%ZeroValue%Depth) then
            !% --- no flow for effectively zero depth
            Flowrate = zeroR
        else
            !% --- use the curve for this pump to set the flowrate for this pump
            !%     based on the inlet depth
            !%     Note that CurveID stores the element index for the pump, and
            !%     the lookup table stores the result in elemR(idx,er_Flowrate)
            !%     for the idx of the pump.
            call util_curve_lookup_singular ( &
                CurveID, er_Depth, er_Flowrate, curve_pump_Xvar, curve_pump_flowrate,1)
        end if

         !% --- set the head of the pump to the upstream face head
        Head = UpFaceHead

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
            integer, pointer    :: iupf
            real(8), pointer    :: Flowrate, Head, UpFaceHead, UpFaceFlowrate
        !%------------------------------------------------------------------
        !% Aliases:
            Flowrate       => elemR(idx,er_Flowrate)
            Head           => elemR(idx,er_Head)
            iupf           => elemI(idx,ei_Mface_uL)
            !% face data used
            UpFaceHead     => faceR(iupf,fr_Head_d)
            UpFaceFlowrate => faceR(iupf,fr_Flowrate)
        !%------------------------------------------------------------------
        !%--- use the upstream face flowrate and head
        Flowrate = UpFaceFlowrate
        Head = UpFaceHead

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