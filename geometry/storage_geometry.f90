module storage_geometry
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Storage geometry computations
    !%==========================================================================
    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use utility_interpolate
    use utility_allocate
    use utility

    implicit none

    private

    public :: storage_functional_depth_from_volume
    public :: storage_tabular_depth_from_volume
    public :: storage_implied_depth_from_volume
    public :: storage_plan_area_from_volume
    public :: storage_create_integrated_volume_curve
    public :: storage_interpolate_volume_from_depth_singular
    public :: storage_create_curve_from_function
    public :: storage_volume_from_depth_singular
    !public :: storage_implied_volume_from_depth_singular
    !public :: storage_functional_volume_from_depth_singular
    !public :: storage_tabular_volume_from_depth_singular



    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================    
!%
    subroutine storage_functional_depth_from_volume (thisP, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on JM that has storage by functional depth
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: thisP(:), Npack
            integer, pointer :: curveID
            real(8), pointer :: depth, fulldepth, volume
            real(8), pointer :: aConst, aCoeff, aExpon, bb, vv
            real(8) :: aa
            integer :: ii
        !%------------------------------------------------------------------

        !% --- HACK: Find out a way to code this without do loop
        do ii = 1,Npack
            depth     => elemR(thisP(ii),er_Depth)
            fulldepth => elemR(thisP(ii),er_FullDepth)
            volume    => elemR(thisP(ii),er_Volume)
            aConst    => elemSR(thisP(ii),esr_Storage_Constant)
            aCoeff    => elemSR(thisP(ii),esr_Storage_Coefficient)
            aExpon    => elemSR(thisP(ii),esr_Storage_Exponent)
            curveID   => elemSI(thisP(ii),esi_JM_Curve_ID)

            if (aExpon == zeroR) then
                !% ---- if aExpon =  0, an explicit depth vs volume relation can be retrived
                depth = volume / (aConst + aCoeff)
            else
                if (aConst == zeroR) then
                    !% ----if aConst =  0, an explicit depth vs volume relation can be retrived
                    depth = (volume * ((aExpon + oneR)/ aCoeff) ) ** (oneR / (aExpon + oneR)) !% bugfix 20220720brh
                else
                    if (aExpon == oneR) then
                        !% --- for aExpon = 1 we have a quadratic relationship
                        aa =  aCoeff/(aExpon + oneR)
                        bb => aConst
                        vv => volume
                        !% --- quadratic formula, Note that the b^2 - 4av is written
                        !%     as b^2 + 4ac because v (the volume) is negative in the
                        !%     derivation but the pointer vv is to a positive value.
                        !%     This is guaranteed to be positive as long as aCoeff,
                        !%     aExpon, and volume are positive
                        depth = (- bb + sqrt(bb**2 + fourR * aa * vv) ) / (twoR * aa)
                    else            
                        !% --- interpolate from the curve created in initial_condition/init_IC_get_junction_data
                        !      using storage_create_curve_from_function()
                        call util_curve_lookup_singular(curveID, er_Volume, er_Depth, curve_storage_volume, &
                            curve_storage_depth, 1)
                    end if
                end if
            endif

            !% --- limit depth to the full depth of the element
            depth = min(depth,fulldepth)

        end do

    end subroutine storage_functional_depth_from_volume
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine storage_tabular_depth_from_volume (thisP, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on JM that has tabular storage (or non-surcharged trapezoidal conduits)
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: thisP(:), Npack
            integer, pointer :: curveID
            real(8), pointer :: depth, fulldepth, volume
            real(8), pointer :: aConst, aCoeff, aExpon
            integer :: ii
        !%------------------------------------------------------------------
        
        !% --- HACK: Find out a way to code this without do loop
        do ii = 1, Npack
            depth     => elemR(thisP(ii),er_Depth)
            fulldepth => elemR(thisP(ii),er_FullDepth)
            volume    => elemR(thisP(ii),er_Volume)
            curveID   => elemSI(thisP(ii),esi_JM_Curve_ID)

            !% --- interpolate from the curve
            call util_curve_lookup_singular(curveID, er_Volume, er_Depth, curve_storage_volume, &
                curve_storage_depth,1)

            !% --- limit depth to the full depth of the element
            depth = min(depth,fulldepth)

        end do

    end subroutine storage_tabular_depth_from_volume
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine storage_implied_depth_from_volume (thisP, Npack)
        !%-------------------------------------------------------------------
        !% Description:
        !% Computes depth from volume for a junction that does not have SWMM
        !% geometry specified using the storage plan area
        !%-------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: thisP(:), Npack
            real(8), pointer :: depth(:), fulldepth(:), pArea(:), volume(:)
            real(8), pointer :: slotVolume(:), slotDepth(:), slotArea(:)
            real(8), pointer :: fullVolume(:)
            logical, pointer :: isSlot(:)
        !%-------------------------------------------------------------------
        !% Aliases
            depth     => elemR(:,er_Depth)
            fulldepth => elemR(:,er_FullDepth)
            fullVolume=> elemR(:,er_FullVolume)
            pArea     => elemSR(:,esr_Storage_Plan_Area)
            volume    => elemR(:,er_Volume)
            slotVolume=> elemR(:,er_SlotVolume)
            slotDepth => elemR(:,er_SlotDepth)
            slotArea  => elemR(:,er_SlotArea)
            isSlot    => elemYN(:,eYN_isPSsurcharged)
        !%-------------------------------------------------------------------

        depth(thisP) = volume(thisP) / pArea(thisP)

        !% --- limit depth to the full depth of the element
        depth(thisP) = min(depth(thisP),fulldepth(thisP))
       
        !%-------------------------------------------------------------------
    end subroutine storage_implied_depth_from_volume
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine storage_plan_area_from_volume (thisP, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on JM that has storage by functional depth
        !%------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:), Npack
            integer, pointer :: curveID
            real(8), pointer :: planArea, tempPArea, volume
            integer :: ii
        !%------------------------------------------------------------------
        do ii = 1, Npack
            tempPArea => elemR(thisP(ii),er_Temp01)
            volume    => elemR(thisP(ii),er_Volume)
            planArea  => elemSR(thisP(ii),esr_Storage_Plan_Area)
            curveID   => elemSI(thisP(ii),esi_JM_Curve_ID)

            !% --- plan area is unchanged for implied or no storage
            if (elemSI(thisP(ii),esi_JM_type) == ImpliedStorage) return
            if (elemSI(thisP(ii),esi_JM_type) == NoStorage) return
           
            !% --- interpolate from the curve created in initial_condition/init_IC_get_junction_data
            !      using storage_create_curve_from_function()
            call util_curve_lookup_singular(curveID, er_Volume, er_Temp01, curve_storage_volume, &
                 curve_storage_area, 1)

            !% --- copy the data over from elemR temporary space to elemSR(thisP,esr_Storage_Plan_Area)
            planArea = tempPArea

        end do

    end subroutine storage_plan_area_from_volume
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine storage_create_integrated_volume_curve (CurveID)
        !%------------------------------------------------------------------
        !% Description:
        !% Integrate the SWMM5 area vs depth curve into volume vs depth curve
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: CurveID
            real(8), pointer    :: depth(:), area(:), integrated_volume(:)
            integer :: ii
        !%------------------------------------------------------------------
        !% Aliases
            depth => curve(curveID)%ValueArray(:,curve_storage_depth)
            area  => curve(curveID)%ValueArray(:,curve_storage_area)
            integrated_volume => curve(curveID)%ValueArray(:,curve_storage_volume)
            integrated_volume = nullvalueR
        !%------------------------------------------------------------------

        do ii = 1,size(curve(curveID)%ValueArray,1)

            if (ii == 1) then
                integrated_volume(ii) = onehalfR * depth(ii) * area(ii)
            else
                integrated_volume(ii) = integrated_volume(ii-oneI) &
                         + onehalfR * (depth(ii) - depth(ii-1)) * (area(ii) + area(ii-1))
            end if

        end do 

    end subroutine storage_create_integrated_volume_curve
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine storage_interpolate_volume_from_depth_singular (StorageIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% interpolate from the storage curve to get volume (stored in elemR)
        !% from the depth (stored in elemR)
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)  :: StorageIdx
            integer, pointer     :: curveID
            character(64) :: subroutine_name = 'storage_interpolate_from_curve'
        !%------------------------------------------------------------------
        !% Aliases
            curveID  => elemSI(StorageIdx,esi_JM_Curve_ID)
        !%------------------------------------------------------------------

        call util_curve_lookup_singular(curveID, er_Depth, er_Volume, curve_storage_depth, &
            curve_storage_volume, 1)

    end subroutine storage_interpolate_volume_from_depth_singular
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine storage_create_curve_from_function (StorageIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% create an artificial storage curve for the functional storage
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)  :: StorageIdx
            real(8), pointer     :: fullDepth, aConst, aCoeff, aExpon
            integer, pointer     :: CurveID, nRow
            character(64) :: subroutine_name = 'storage_create_curve_from_function'
        !%------------------------------------------------------------------
        !% Aliases
            fullDepth => elemR(StorageIdx,er_FullDepth)
            CurveID   => elemSI(StorageIdx,esi_JM_Curve_ID)
            aConst    => elemSR(StorageIdx,esr_Storage_Constant)
            aCoeff    => elemSR(StorageIdx,esr_Storage_Coefficient)  
            aExpon    => elemSR(StorageIdx,esr_Storage_Exponent)
        !%------------------------------------------------------------------

        !% --- num of rows (HACK: read from a setting file)
        nRow =>  setting%Junction%FunStorageN

        !% --- add a new curveID
        N_FunctionalStorage = N_FunctionalStorage + oneI
        CurveID  = setting%SWMMinput%N_curve + N_FunctionalStorage
      
        curve(CurveID)%ID       = CurveID
        curve(CurveID)%Type     = StorageCurve
        curve(CurveID)%NumRows  = nRow 
        curve(CurveID)%ElemIdx  = StorageIdx

        !% --- allocate the valueArray for the new curve
        call util_allocate_curve_entries (CurveID, nRow)

        Curve(CurveID)%ValueArray(:,curve_storage_depth)  = util_linspace(zeroR,fullDepth,nRow)
        Curve(CurveID)%ValueArray(:,curve_storage_area)   = aConst + aCoeff * &
                    Curve(CurveID)%ValueArray(:,curve_storage_depth) ** aExpon
        Curve(CurveID)%ValueArray(:,curve_storage_volume) =  aConst * Curve(CurveID)%ValueArray(:,curve_storage_depth) + &
                    (aCoeff/(aExpon+oneR)) * Curve(CurveID)%ValueArray(:,curve_storage_depth) ** (aExpon+oneR)

    end subroutine storage_create_curve_from_function
!%  
!%==========================================================================
!%==========================================================================
!%
    real(8) function storage_volume_from_depth_singular & 
        (idx, indepth) result(outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% computes the volume from the depth in a junction
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: idx
            real(8), intent(in) :: indepth
        !%------------------------------------------------------------------

        select case(elemSI(idx,esi_JM_Type))
            case (NoStorage)
                !% --- no storage is assigned
                outvalue = zeroR   

            case (ImpliedStorage)
                outvalue = storage_implied_volume_from_depth_singular(idx, indepth)

            case (TabularStorage)
                outvalue = storage_tabular_volume_from_depth_singular(idx, indepth)

            case (FunctionalStorage)
                outvalue = storage_functional_volume_from_depth_singular(idx, indepth)

            case default
                print *, 'CODE ERROR unexpected case default for storage '
                print *, 'problem in (elemSI(idx,esi_JM_Type) '
                print *, 'element index (idx)= ',idx
                print *, 'junction main type of ',(elemSI(idx,esi_JM_Type))
                print *, trim(reverseKey(elemSI(idx,esi_JM_Type)))
        end select

    end function storage_volume_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function storage_implied_volume_from_depth_singular &
        (idx, indepth) result(outvalue)
        !%-------------------------------------------------------------------
        !% Description:
        !% computes the volume from the depth in an implied storage junction
        !% ONLY APPLICABLE TO NON-SURCHARGED
        !%-------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: idx
            real(8), intent(in) :: indepth
        !%-------------------------------------------------------------------
        outvalue = elemSR(idx,esr_Storage_Plan_Area) * indepth
        !print *, 'in storage ',elemSR(idx,esr_Storage_Plan_Area), indepth
            
    end function storage_implied_volume_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function storage_functional_volume_from_depth_singular &
        (idx, indepth) result(outvalue)
        !%-------------------------------------------------------------------
        !% Description:
        !% computes the volume from the depth in a functional storage junction
        !%-------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: idx
            real(8), intent(in) :: indepth
            real(8), pointer :: aConst, aCoeff, aExpon
        !%-------------------------------------------------------------------
        !% Aliases
            aConst  => elemSR(idx,esr_Storage_Constant)
            aCoeff  => elemSR(idx,esr_Storage_Coefficient)
            aExpon  => elemSR(idx,esr_Storage_Exponent)
        !%-------------------------------------------------------------------
        !% --- integrate the functional form
        outvalue =  aConst * indepth      &
                    + (aCoeff / (aExpon + oneR)) * (indepth ** (aExpon+ oneR))
            
    end function storage_functional_volume_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%    
    real(8) function storage_tabular_volume_from_depth_singular &
        (idx, indepth)  result(outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% computes the volume from the depth in a tabular storage junction
        !% This uses the "storage_interpolate_volume_from_depth_singular that
        !% operates in place on elemR(:,er_Volume) based on elemR(:,er_Depth)
        !% However, this function is designed to leave both of the original
        !% values untouched.
        !% This should mainly be used in initial conditions, not elsewhere.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: idx
            real(8), intent(in) :: indepth
            real(8)             :: tdepth, tvolume
        !%------------------------------------------------------------------
        !% --- store current volume and depth
        tvolume = elemR(idx,er_Volume)   
        tdepth  = elemR(idx,er_Depth) 
      
        !% --- temporary switch of indepth into the current depth for interpolation
        elemR(idx,er_Depth) = indepth

        !% --- call the interpolation that works on elemR(:,er_Volume)
        !%     based on value in elemR(:,er_Depth)
        call storage_interpolate_volume_from_depth_singular (idx)

        !% --- get the output value
        outvalue = elemR(idx,er_Volume)

        !% --- switch the original incoming values back
        elemR(idx,er_Depth)  = tdepth
        elemR(idx,er_Volume) = tvolume

    end function storage_tabular_volume_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
end module storage_geometry