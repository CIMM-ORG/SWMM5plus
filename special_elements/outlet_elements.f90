module outlet_elements
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Computes outlet flow and head
    !%
    !%==========================================================================
    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use common_elements
    use adjust
    use utility_interpolate
    use utility_crash, only: util_crashpoint

    implicit none

    private

    public :: outlet_toplevel
    public :: outlet_set_setting

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine outlet_toplevel (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes diagnostic flow and head delta across a outlet.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx  !% must be a single element ID      
            character(64) :: subroutine_name = 'outlet_toplevel'
        !%------------------------------------------------------------------

        !% --- update the setting
        call outlet_set_setting (eIdx)

        !% --- get the flow direction and element head
        call  common_head_and_flowdirection_singular &
            (eIdx, esr_Outlet_Zcrest, esr_Outlet_NominalDownstreamHead, esi_Outlet_FlowDirection)
        
        !% --- find effective head difference accross outlet element
        call outlet_effective_head_delta (eIdx)
        
        !% --- find flow on outlet element
        call outlet_flow (eIdx)

        !% --- update outlet geometry
        call outlet_geometry_update (eIdx)
        
        !% --- update velocity from flowrate and area
        call common_velocity_from_flowrate_singular (eIdx)
        
    end subroutine outlet_toplevel
!%
!%==========================================================================
!%==========================================================================    
!%     
    subroutine outlet_set_setting (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% Updates the setting on an outlet element, following EPA-SWMM
        !% link.c/link_setSetting where outlet setting is immediately set to 
        !% targetsetting after a control change
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIdx
        !%------------------------------------------------------------------   

        elemR(eIdx,er_Setting) = elemR(eIdx,er_TargetSetting)

        !% --- error check
        !%     EPA-SWMM allows the outlet setting to be only 0.0 or 1.0
        if (.not. ((elemR(eIdx,er_Setting) == 0.0) .or. (elemR(eIdx,er_Setting) == 1.0))) then
            print *, 'CODE ERROR: outlet element has er_Setting that is not 0.0 or 1.0'
            call util_crashpoint(6387233)
        end if

    end subroutine outlet_set_setting
!%
!%========================================================================== 
!% PRIVATE
!%==========================================================================   
!%  
    subroutine outlet_effective_head_delta (eIdx)
        !%------------------------------------------------------------------   
        !% Description:
        !% Computes the effective head difference flowing over the top of a outlet
        !%------------------------------------------------------------------       
        !% Declarations
            integer, intent(in) :: eIdx !% single ID of element
            real(8), pointer :: EffectiveHeadDelta, Head, Zcrest
            real(8), pointer :: NominalDownstreamHead
            integer, pointer :: OutletType, FlowDirection
            logical, pointer :: hasFlapGate
        !%------------------------------------------------------------------ 
        !% Aliases
            OutletType            => elemSI(eIdx,esi_Outlet_SpecificType)
            FlowDirection         => elemSI(eIdx,esi_Outlet_FlowDirection)
            Head                  => elemR(eIdx,er_Head)
            Zcrest                => elemSR(eIdx,esr_Outlet_Zcrest)
            NominalDownstreamHead => elemSR(eIdx,esr_Outlet_NominalDownstreamHead)
            hasFlapGate           => elemYN(eIdx,eYN_hasFlapGate)
            !% --- output
            EffectiveHeadDelta    => elemSR(eIdx,esr_Outlet_EffectiveHeadDelta)
        !%------------------------------------------------------------------ 
        !% --- find the effective head delta
        if (hasFlapGate .and. (FlowDirection < zeroR)) then
            EffectiveHeadDelta = zeroR
        else
            select case (OutletType)
                case (func_head_outlet, tabl_head_outlet)
                    EffectiveHeadDelta = Head - max(NominalDownstreamHead,Zcrest)

                case (func_depth_outlet, tabl_depth_outlet) 
                    EffectiveHeadDelta = Head - Zcrest

                case default
                    print *, 'CODE ERROR outlet type unknown for # ', OutletType
                    print *, 'which has key ',trim(reverseKey(OutletType))
                    stop 5878
            end select
        end if
        
    end subroutine outlet_effective_head_delta
!%
!%========================================================================== 
!%==========================================================================    
!%  
    subroutine outlet_flow (eIdx)
        !%------------------------------------------------------------------ 
        !% Description:
        !% Computes flow at outlet
        !%------------------------------------------------------------------ 
            integer, intent(in) :: eIdx !% must be single element ID
            integer, pointer :: FlowDirection, CurveID, OutletType
            real(8), pointer :: Flowrate, Depth, CurrentSetting, EffectiveHeadDelta, qCoeff, qExpon, dQdH
            logical, pointer :: hasFlapGate
            real(8) :: CoeffOrifice
        !%------------------------------------------------------------------ 
        !% Aliases
            OutletType         => elemSI(eIdx,esi_Outlet_SpecificType)
            FlowDirection      => elemSI(eIdx,esi_Outlet_FlowDirection)
            CurveID            => elemSI(eIdx,esi_Outlet_CurveID)
            Depth              => elemR (eIdx,er_Depth)
            dQdH               => elemSR(eIdx,esr_Outlet_dQdHe)
            Flowrate           => elemR (eIdx,er_Flowrate) 
            CurrentSetting     => elemR (eIdx,er_Setting) 
            qCoeff             => elemSR(eIdx,esr_Outlet_Coefficient)
            qExpon             => elemSR(eIdx,esr_Outlet_Exponent)
            EffectiveHeadDelta => elemSR(eIdx,esr_Outlet_EffectiveHeadDelta)
            hasFlapGate        => elemYN(eIdx,eYN_hasFlapGate)
        !%------------------------------------------------------------------ 
 
        if (hasFlapGate .and. (FlowDirection < zeroR)) then
            Depth = zeroR
            Flowrate = zeroR
            dQdH = zeroR
        else
            select case (OutletType)
                case (func_head_outlet, func_depth_outlet) 
                    Depth = EffectiveHeadDelta
                    Flowrate = CurrentSetting * real(FlowDirection,8) * qCoeff * (effectiveHeadDelta ** qExpon)
                    !% --- dQdH is zero for outlet elements
                    dQdH = zeroR
                case (tabl_head_outlet, tabl_depth_outlet)

                    Depth = EffectiveHeadDelta
                    call util_curve_lookup_singular(CurveID, er_Depth, er_Flowrate, &
                        curve_outlet_depth, curve_outlet_flowrate,1)
                    Flowrate = Flowrate * CurrentSetting * real(FlowDirection,8)
                    !% --- dQdH is zero for outlet elements
                    dQdH = zeroR
                case default
                    print *, 'CODE ERROR outlet type unknown for # ', OutletType
                    print *, 'which has key ',trim(reverseKey(OutletType))
                    stop 7824
            end select
        end if

    end subroutine outlet_flow
!%
!%========================================================================== 
!%==========================================================================    
!%  
    subroutine outlet_geometry_update (eIdx)
        !%------------------------------------------------------------------ 
        !% Description:
        !% Outlet geometry is area, which uses the average of the upstream
        !% and downstream faces
        !%------------------------------------------------------------------ 
        !% Declarations
            integer, intent(in) :: eIdx
            integer, pointer :: fUp, fDn
            real(8), pointer :: fAUp, fADn, Area
        !%------------------------------------------------------------------ 
        !% Aliases
            fUp  => elemI(eIdx,ei_Mface_uL)
            fDn  => elemI(eIdx,ei_Mface_dL)
            fAUp => faceR(fUp,fr_area_d)
            fADn => faceR(fDn,fr_area_u)
            Area => elemR(eIdx,er_Area)
        !%------------------------------------------------------------------ 
        !% --- average area of upstream and downstream faces
        Area      =  (fAUp + fADn) / twoR

        !% --- apply geometry limiters
        call adjust_limit_by_zerovalues_singular (eIdx, er_Area, setting%ZeroValue%Area, .false.)

    end subroutine outlet_geometry_update
!%    
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module outlet_elements