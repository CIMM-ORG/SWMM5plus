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
                call pump_type2or4(eIdx,0)
            
            case (type3_Pump)
                call pump_type3(eIdx)

            case (type4_Pump)
                call pump_type2or4(eIdx,1)
                    
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
            integer, pointer    :: CurveID
            integer             :: Ci, Aidx
            real(8), pointer    :: FlowRate, Volume, Head, Depth, dQdH
            real(8)             :: upDepth, upHead, upVolume, upFlowrate, maxFlowrate
        !%------------------------------------------------------------------
        !% Aliases:
            CurveID  => elemSI(eIdx,esi_Pump_CurveID)
            FlowRate => elemR(eIdx,er_Flowrate)
            dQdH     => elemR(eIdx,er_dQdH)
            Head     => elemR(eIdx,er_Head)
            Depth    => elemR(eIdx,er_Depth)
            Volume   => elemR(eIdx,er_Volume)
        !%------------------------------------------------------------------
        !% Preliminaries
            upDepth=zeroR; upHead=zeroR; upVolume=zeroR; upFlowrate=zeroR
            maxFlowrate=zeroR
            Ci=1; Aidx=eIdx
        !%------------------------------------------------------------------

        !% -- get upstream data
        call pump_upstream_data &
            (type1_Pump, eIdx, Ci, Aidx, upDepth, upHead, upVolume, upFlowrate, maxFlowrate)

        !% --- store the upstream element depth as the pump depth
        Depth = upDepth

        !% --- set the head of the pump to the upstream element head
        Head  = upHead

        !% --- set the volume of the pump to the upstream element volume
        !%     which is necessary for curve lookup
        !%     NOTE: this volume should NOT be included in conservation calcs
        Volume = upVolume

        !% --- turn pump on or off
        call pump_turn_onoff(eIdx, Depth)

        !% --- exit if pump is off or depth is too small
        if ((elemR(eIdx,er_Setting) == zeroR) .or. &
            (upVolume .le. setting%ZeroValue%Volume)) then 
            Flowrate = zeroR
            Head     = zeroR
            return 
        end if

        ! print *, ' '
        ! print *, 'in pump '
        ! ! print *, 'CurveID ',CurveID
        ! ! print *, 'element ',eIdx, curve(CurveID)%ElemIdx
        ! print *, 'flowrate before ', Flowrate
        ! print *, 'volume before   ', Volume

        !% --- use the curve without interpolation
        call util_curve_lookup_singular( &
            CurveID, er_Volume, er_Flowrate, curve_pump_Xvar, curve_pump_flowrate,0)

         
        ! print *, 'flowrate after ', Flowrate

        !% --- flow limitation
        Flowrate = min(Flowrate,maxFlowrate)

        !% --- dQ/dH is zero for type 1 pump
        dQdH = zeroR

        ! print *, 'flowrate limit ', Flowrate

        !% --- reset pump depth and head to zero (no meaning)
        Depth = zeroR
        Head  = zeroR

    end subroutine pump_type1
!%
!%==========================================================================   
!%==========================================================================    
!%  
    subroutine pump_type2or4 (eIdx,interpType)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes flowrate for pump where flow is set by upstream depth
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx  !% pump index
            integer, intent(in) :: interpType !% 0 for type1, 1 for type4
            integer, pointer    :: CurveID
            integer             :: Ci, Aidx 
            real(8), pointer    :: FlowRate, Head, Depth, dQdH, tempFlowrate, tempDepth
            real(8)             :: upDepth, upHead, upVolume, upFlowrate, maxFlowrate, dH
        !%------------------------------------------------------------------
        !% Aliases:
            CurveID   => elemSI(eIdx,esi_Pump_CurveID)
            FlowRate  => elemR (eIdx,er_Flowrate)
            Head      => elemR (eIdx,er_Head)
            Depth     => elemR (eIdx,er_Depth)
            dQdH      => elemR (eIdx,er_dQdH)
            tempFlowrate => elemR(eIdx,er_Temp01)
            tempDepth    => elemR(eIdx,er_Temp02)
        !%-----------------------------------------------------------------
        !% Preliminaries
            upDepth=zeroR; upHead = zeroR; upVolume=zeroR; upFlowrate=zeroR
            maxFlowrate=zeroR
        !%------------------------------------------------------------------

        !% -- get upstream data
        call pump_upstream_data &
            (type2_Pump, eIdx, Ci, Aidx, upDepth, upHead, upVolume, upFlowrate, maxFlowrate)

        !% --- store the upstream element depth as the pump depth
        Depth = upDepth

        !% --- set the head of the pump to the upstream element head
        Head  = upHead
 
        !% --- turn pump on or off
        call pump_turn_onoff(eIdx, Depth)

        ! print *, 'depth, head ',Depth, Head
        ! print *, 'pump on/off ',elemR(eIdx,er_Setting)

        !% --- exit if pump is off or depth is too small
        if ((elemR(eIdx,er_Setting) == zeroR) .or. &
            (Depth .le. setting%ZeroValue%Depth)) then 
            Flowrate = zeroR
            Head     = zeroR
            Depth    = zeroR
            dQdH     = zeroR
            return 
        end if

        !% --- use the curve without interpolation
        call util_curve_lookup_singular( &
            CurveID, er_Depth, er_Flowrate, curve_pump_Xvar, curve_pump_flowrate,interpType)

        ! print *, 'Flowrate ',Flowrate, maxFlowrate

        !% --- flow limitation
        Flowrate = min(Flowrate,maxFlowrate)

        if (interpType == 0) then
            dQdH = zeroR
        else
            dH = 0.001
            tempDepth = Depth + dH
            
            !% interpolate for new temp depth and temp flowrate
            call util_curve_lookup_singular( &
                CurveID, er_Temp02, er_Temp01, curve_pump_Xvar, curve_pump_flowrate,interpType)
            
            dQdH = (tempFlowrate - Flowrate) / dH
        end if

        !% --- reset pump depth and head to zero (no meaning)
        Depth = zeroR
        Head  = zeroR


    end subroutine pump_type2or4
!%
!%==========================================================================
!%==========================================================================    
!%
    subroutine pump_type3 (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% Centrifugal pump whose flow depends on head difference across
        !% upstream and downstream elements
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx  !% pump index
            integer             :: Ci, Aidx
            integer, pointer    :: CurveID
            real(8), pointer    :: Head, Depth, Flowrate, dQdH, tempFlowrate, tempDepth
            real(8)             :: upDepth, upHead, upVolume, upFlowrate, maxFlowrate
            real(8)             :: dnHead, dH
        !%------------------------------------------------------------------
        !% Aliases:
            CurveID  => elemSI(eIdx,esi_Pump_CurveID)
            Head     => elemR(eIdx,er_Head)
            Depth    => elemR(eIdx,er_Depth)
            dQdH     => elemR(eIdx,er_dQdH)
            FlowRate => elemR(eIdx,er_Flowrate)
            tempFlowrate => elemR(eIdx,er_Temp01)
            tempDepth    => elemR(eIdx,er_Temp02)
        !%------------------------------------------------------------------
        !% Preliminaries
            upDepth=zeroR; upHead = zeroR; upVolume=zeroR; upFlowrate=zeroR
            maxFlowrate=zeroR; dnHead=zeroR
        !%------------------------------------------------------------------
        !% --- get upstream data
        call pump_upstream_data &
            (type3_Pump, eIdx, Ci, Aidx, upDepth, upHead, upVolume, upFlowrate, maxFlowrate)
        
        !% --- get the downstream data
        call pump_downstream_data (eIdx, Ci, Aidx, dnHead) 

        !% --- temporarily set the head at the pump to the head difference across the 
        !%     pump that controls a Type 3 (centrifugal) pump. This reset is needed
        !%     for the curve lookup  
        Head = max((dnHead - upHead), zeroR)

        !% --- store the depth upstream as the pump element depth
        Depth = upDepth

        ! print *, 'Delta head, depth ',Head, Depth

        !% --- turn pump on or off
        call pump_turn_onoff(eIdx, Depth)

        ! print *, 'setting ',elemR(eIdx,er_Setting)

        !% --- exit if pump is off or depth is too small
        if ((elemR(eIdx,er_Setting) == zeroR) .or. &
            (Depth .le. setting%ZeroValue%Depth)) then 
            Flowrate = zeroR
            Head     = zeroR
            Depth    = zeroR
            dQdH     = zeroR
            return 
        end if

        !% --- use the curve for this pump and the head difference
        call util_curve_lookup_singular( &
            CurveID, er_Head, er_Flowrate, curve_pump_Xvar, curve_pump_flowrate,1)

        ! print *, 'flowrate ',flowrate

        !% --- dQdh calculation
        dH = 0.001
        tempDepth = Depth + dH
        
        !% interpolate for new temp depth and temp flowrate
        call util_curve_lookup_singular( &
            CurveID, er_Temp02, er_Temp01, curve_pump_Xvar, curve_pump_flowrate,1)
        
        dQdH = (tempFlowrate - Flowrate) / dH

        !% --- reset the pump element head to the upstream value
        Head  = zeroR
        Depth = zeroR

        !% --- flow limitation
        Flowrate = min(Flowrate,maxFlowrate)

    end subroutine pump_type3
!%
!%==========================================================================
! !%========================================================================== 
! !%  
    !     subroutine pump_type4 (eIdx)
    !         !%------------------------------------------------------------------
    !         !% Description:
    !         !% 
    !         !%------------------------------------------------------------------
    !         !% Declarations:
    !             integer, intent(in) :: eIdx  !% pump index
    !             integer, pointer :: iupf,  CurveID
    !             real(8), pointer :: FlowRate, UpFaceHead
    !             real(8), pointer :: Depth, Head, UpFaceZbottom
    !         !%------------------------------------------------------------------
    !         !% Aliases:
    !             CurveID  => elemSI(eIdx,esi_Pump_CurveID)
    !             Depth    => elemR(eIdx,er_Depth)
    !             Flowrate => elemR(eIdx,er_Flowrate)
    !             Head     => elemR(eIdx,er_Head)
    !             !% face locations
    !             iupf     => elemI(eIdx,ei_Mface_uL)
    !             !% face data used
    !             UpFaceHead    => faceR(iupf,fr_Head_d)
    !             UpFaceZbottom => faceR(iupf,fr_Zbottom)     
    !         !%------------------------------------------------------------------
    !         !% -- get upstream data
    !             call pump_upstream_data &
    !             (type4_Pump, eIdx, Ci, Aidx, upDepth, upHead, upVolume, upFlowrate, maxFlowrate)

    !         !% --- store the upstream element depth as the pump depth
    !         Depth = upDepth

    !         !% --- set the head of the pump to the upstream element head
    !         Head  = upHead

    !         !% --- turn pump on or off
    !         call pump_turn_onoff(eIdx, Depth)

    !         !% --- exit if pump is off or depth is too small
    !         if ((elemR(eIdx,er_Setting) == zeroR) .or. &
    !             (Depth .le. setting%ZeroValue%Depth)) then 
    !             Flowrate = zeroR
    !             return 
    !         end if

    !         !% --- use the curve for this pump to set the flowrate for this pump
    !         !%     based on the inlet depth
    !         !%     Note that CurveID stores the element index for the pump, and
    !         !%     the lookup table stores the result in elemR(idx,er_Flowrate)
    !         !%     for the idx of the pump.
    !         call util_curve_lookup_singular ( &
    !             CurveID, er_Depth, er_Flowrate, curve_pump_Xvar, curve_pump_flowrate,1)


    !          !% --- set the head of the pump to the upstream face head
    !         Head = UpFaceHead

    !         !%------------------------------------------------------------------
!     end subroutine pump_type4
!%
!%==========================================================================
!%==========================================================================    
!%  
    subroutine pump_ideal (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% EPA-SWMM Ideal Pump -- which simply repeats its inlet condition
        !% as a flow condition
        !%
        !% HACK---EPA-SWMM uses the upstream node flow plus any overflow rate
        !%        In the following we have not included any overflow 
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx  !% pump index
            integer             :: Ci, Aidx 
            real(8), pointer    :: Flowrate, Head, Depth
            real(8)             :: upDepth, upHead, upVolume, upFlowrate, maxFlowrate
        !%------------------------------------------------------------------
        !% Aliases:
            Flowrate       => elemR(eIdx,er_Flowrate)
            Head           => elemR(eIdx,er_Head)
            Depth          => elemR(eIdx,er_Depth)
        !%------------------------------------------------------------------
       ! % Preliminaries
            upDepth=zeroR; upHead = zeroR; upVolume=zeroR; upFlowrate=zeroR
            maxFlowrate=zeroR
        !%------------------------------------------------------------------
        !% -- get upstream data
        call pump_upstream_data &
            (type_IdealPump, eIdx, Ci, Aidx, upDepth, upHead, upVolume, upFlowrate, maxFlowrate)

        !% --- store the upstream element depth as the pump depth
        Depth = upDepth

        !% --- set the head of the pump to the upstream element head
        Head  = upHead

        !% --- set the pump flowrate to the upstream element flowrate
        Flowrate = upFlowrate

        !% --- flow limitation
        Flowrate = min(Flowrate,maxFlowrate)

        !%------------------------------------------------------------------
        !% Closing:
    end subroutine pump_ideal
!%
!%==========================================================================
!%==========================================================================   
!%  
    subroutine pump_upstream_data &
        (pumpType, eIdx, Ci, Aidx, Depth, Head, Volume, Flowrate, maxFlowrate) 
        !%------------------------------------------------------------------
        !% Description:
        !% Upstream element data extracted across images as needed
        !% All pumps require Depth and Head.
        !% Volume is always computed for maxFlowrate calcuation
        !% Flowrate is only for Ideal pumps
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: pumpType, eIdx
            integer, intent(inout) :: Ci, Aidx
            real(8), intent(inout) :: Depth, Head, Volume, Flowrate, maxFlowrate
            integer, pointer       :: Fidx
            integer                :: AidxJB
        !%------------------------------------------------------------------
        !% Aliases
            Fidx     => elemI(eIdx,ei_Mface_uL) !% face upstream
        !%------------------------------------------------------------------
        !%  
             !print *, 'eIdx, Fidx ',eIdx,Fidx

        !% --- identify the upstream adjacent element and its image
        if (elemYN(eIdx,eYN_isBoundary_up)) then 
            Ci   = faceI(Fidx,fi_Connected_image)
            Aidx = faceI(Fidx,fi_GhostElem_uL)
            sync images(Ci)
        else
            Ci   =  this_image()
            Aidx =  faceI(Fidx,fi_Melem_uL)
        end if
        AidxJB = Aidx

        !% --- upstream switch Aidx to JM if junction
        if (elemI(Aidx,ei_elementType)[Ci] == JB) then
            !% --- for branches, switch to junction    
            Aidx   = elemSI(Aidx,esi_JunctionBranch_Main_Index)[Ci]
            Volume =  elemR(Aidx,er_Volume)[Ci] 
        end if

        ! print *, 'Aidx ',Aidx
        ! print *, 'node number ',elemI(Aidx,ei_node_Gidx_BIPquick)
        ! print *, 'node name   ',trim(node%Names(elemI(Aidx,ei_node_Gidx_BIPquick))%str)

        !% --- store the upstream element depth 
        Depth = elemR(Aidx,er_Depth)[Ci]

        !% --- get the upstream element head
        Head  = elemR(Aidx,er_Head)[Ci]

        !% --- get upstream element volume
        Volume = elemR(Aidx,er_Volume)[Ci]

        !% --- get upstream flowrate if needed
        select case (pumpType)
        case (type_IdealPump) !% ideal pump (use AidxJB when upstream is junction)
            Flowrate = elemR(AidxJB,er_Flowrate)[Ci]
        case default
            !% dummy return
            Flowrate = zeroR
        end select

        !% --- flowrate limitation for small depths
        if (elemYN(Aidx,eYN_isSmallDepth)[Ci]) then 
            !% --- limit flowrate by upstream volume for small depths
            !%     with a factor < 1.
            maxFlowrate = setting%SmallDepth%PumpVolumeFactor & 
                        * Volume / setting%Time%Hydraulics%Dt
        else 
            !% --- limit pump outflow by upstream volume (CFL limit)
            maxFlowrate = Volume / setting%Time%Hydraulics%Dt
        end if

    end subroutine pump_upstream_data
!%
!%========================================================================== 
!%==========================================================================    
!%  
    subroutine pump_downstream_data (eIdx, Ci, Aidx, Head) 
        !%------------------------------------------------------------------
        !% Description:
        !% Downstream element data
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: eIdx
            integer, intent(inout) :: Ci, Aidx
            real(8), intent(inout) :: Head
            integer, pointer       :: Fidx
        !%------------------------------------------------------------------
        !% Aliases
            Fidx     => elemI(eIdx,ei_Mface_dL) !% face dnstream
        !%------------------------------------------------------------------
        !%  
        !% --- identify the downstream adjacent element and its image
        if (elemYN(eIdx,eYN_isBoundary_dn)) then 
            Ci   = faceI(Fidx,fi_Connected_image)
            Aidx = faceI(Fidx,fi_GhostElem_dL)
            sync images(Ci)
        else
            Ci   =  this_image()
            Aidx =  faceI(Fidx,fi_Melem_dL)
        end if

        !% --- get the downstream element head
        Head  = elemR(Aidx,er_Head)[Ci]

    end subroutine pump_downstream_data
!%
!%========================================================================== 
!%==========================================================================    
!%  
    subroutine pump_turn_onoff (eIdx, Depth)
        !%------------------------------------------------------------------
        !% Description:
        !% Turns pump on or off depending on its yOn and yOff depth
        !% values. Sets flowrate to zero if pump is off
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: eIdx
            real(8), intent(in)    :: Depth
        !%------------------------------------------------------------------    

        ! print *, 'eidx ',eIdx
        ! print *, 'pump setting ', elemR(eIdx,er_Setting)
        ! print *, 'is controlled ',elemSI(eIdx,esi_Pump_IsControlled)
        ! print *, 'depth on      ',elemSR(eIdx,esr_Pump_yOn)

        !% --- check pump  off/on status and whether a change is required
        if (elemR(eIdx,er_Setting) == zeroR) then 
            !% --- pump is off
            if (elemSI(eIdx,esi_Pump_IsControlled) == zeroI) then
                !% --- if NOT externally controlled, 
                !%     check if pump should be turned on based on upstream condition
                if (Depth > elemSR(eIdx,esr_Pump_yOn))  then
                    !% -- upstream startup depth exceeded, turn pump on
                    elemR(eIdx,er_Setting) = oneR
                else 
                    !% --- no action because upstream depth is too small
                end if
            else 
                !% -- stays shutoff based on external control
                !%    regardless of upstream depth
            end if
        else 
            !% --- pump is on (due to prior upstream or external control)
            if (Depth .le. elemSR(eIdx,esr_Pump_yOff)) then 
                !% --- low inlet depth (safety shutoff)
                !%     turn pump off regardless of controls
                elemR(eIdx,er_Setting) = zeroR
            end if
        end if
        ! print *, 'pump setting ', elemR(eIdx,er_Setting)

    end subroutine pump_turn_onoff
!%
!%========================================================================== 
!%==========================================================================    
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