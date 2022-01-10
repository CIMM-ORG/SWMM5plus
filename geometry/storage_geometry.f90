module storage_geometry

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use utility_interpolate
    use utility_allocate
    use utility

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Trapezoidal channel geometry
    !%

    private

    public :: storage_functional_depth_from_volume
    public :: storage_tabular_depth_from_volume
    public :: storage_artificial_depth_from_volume
    public :: storage_integrate_area_vs_depth_curve
    public :: storage_interpolate_volume_from_depth_singular
    public :: storage_create_curve

    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine storage_functional_depth_from_volume (elemPGx, Npack, thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Only applies on JM that has storage (or non-surcharged trapezoidal conduits)
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
        integer, pointer :: thisP, curveID
        real(8), pointer :: depth, volume
        real(8), pointer :: aConst, aCoeff, aExpon
        integer :: ii
        !%-----------------------------------------------------------------------------

        !% HACK: Find out a way to code this without do loop
        do ii = 1, Npack
            thisP   => elemPGx(ii,thisCol) 
            depth   => elemR(thisP,er_Depth)
            volume  => elemR(thisP,er_Volume)
            aConst  => elemSR(thisP,esr_Storage_Constant)
            aCoeff  => elemSR(thisP,esr_Storage_Coefficient)
            aExpon  => elemSR(thisP,esr_Storage_Exponent)
            curveID => elemSI(thisP,esi_JunctionMain_Curve_ID)

            if (aExpon == zeroR) then
                !% if aExpon =  0, an explicit depth vs volume relation can be retrived
                depth = volume / (aConst + aCoeff)
            
            elseif (aConst == zeroR) then
                !% if aConst =  0, an explicit depth vs volume relation can be retrived
                depth = (volume / aCoeff * (oneR / (aExpon + oneR))) ** (oneR / (aExpon + oneR))

            else
                !% else interpolate from the curve
                call util_curve_lookup_singular(curveID, er_Volume, er_Depth, curve_storage_volume, &
                    curve_storage_depth)
            endif
        end do

    end subroutine storage_functional_depth_from_volume
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine storage_tabular_depth_from_volume (elemPGx, Npack, thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Only applies on JM that has storage (or non-surcharged trapezoidal conduits)
        !%-----------------------------------------------------------------------------
            integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
            integer, pointer :: thisP, curveID
            real(8), pointer :: depth, volume
            real(8), pointer :: aConst, aCoeff, aExpon
            integer :: ii
        !%-----------------------------------------------------------------------------

        !% HACK: Find out a way to code this without do loop
        do ii = 1, Npack
            thisP   => elemPGx(ii,thisCol) 
            depth   => elemR(thisP,er_Depth)
            volume  => elemR(thisP,er_Volume)
            curveID => elemSI(thisP,esi_JunctionMain_Curve_ID)

            !% interpolate from the curve
            call util_curve_lookup_singular(curveID, er_Volume, er_Depth, curve_storage_volume, &
                curve_storage_depth)
        end do

    end subroutine storage_tabular_depth_from_volume
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine storage_artificial_depth_from_volume (elemPGx, Npack, thisCol)
        !%  
        !%-------------------------------------------------------------------
        !% Description:
        !% Computes depth from volume for a junction that does not have SWMM
        !% geometry specified using the storage plane area
        !%-------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
            integer, pointer :: thisP(:)
            real(8), pointer :: depth(:), pArea(:), volume(:)
        !%-------------------------------------------------------------------
        !% Aliases
            thisP  => elemPGx(1:Npack,thisCol)
            depth  => elemR(:,er_Depth)
            pArea  => elemSR(:,esr_Storage_Plane_Area)
            volume => elemR(:,er_Volume)
        !%-------------------------------------------------------------------

        depth(thisP) = volume(thisP) / pArea(thisP)

        !print *, 'INSTORAGE ',volume(thisP)
       
        !%-------------------------------------------------------------------
    end subroutine storage_artificial_depth_from_volume
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine storage_integrate_area_vs_depth_curve (CurveID)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Intigrate the SWMM5 area vs depth curve into volume vs depth curve
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: CurveID
        real(8), pointer    :: depth(:), area(:), integrated_volume(:)
        integer :: ii
        !%-----------------------------------------------------------------------------
        
        !% pointers allocation
        depth => curve(curveID)%ValueArray(:,curve_storage_depth)
        area  => curve(curveID)%ValueArray(:,curve_storage_area)
        integrated_volume => curve(curveID)%ValueArray(:,curve_storage_volume)
        integrated_volume = nullvalueR

        do ii = 1,size(curve(curveID)%ValueArray,1)
            if (ii == 1) then
                integrated_volume(ii) = onehalfR * depth(ii) * area(ii)
                
            else
                integrated_volume(ii) = onehalfR * (depth(ii) - depth(ii-1)) * &
                    (area(ii) + area(ii-1))
            end if
        end do 

    end subroutine storage_integrate_area_vs_depth_curve
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine storage_interpolate_volume_from_depth_singular (StorageIdx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% interpolate from the storage curve
        !%-----------------------------------------------------------------------------
        integer, intent(in)  :: StorageIdx

        integer, pointer :: curveID

        character(64) :: subroutine_name = 'storage_interpolate_from_curve'
        !%-----------------------------------------------------------------------------
        curveID  => elemSI(StorageIdx,esi_JunctionMain_Curve_ID)

        call util_curve_lookup_singular(curveID, er_Depth, er_Volume, curve_storage_depth, &
            curve_storage_volume)

    end subroutine storage_interpolate_volume_from_depth_singular
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine storage_create_curve (StorageIdx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% create an artificial storage curve for the functional storage
        !%-----------------------------------------------------------------------------
        integer, intent(in)  :: StorageIdx

        real(8), pointer :: fullDepth, aConst, aCoeff, aExpon
        integer, pointer :: CurveID, nRow

        character(64) :: subroutine_name = 'storage_create_curve'
        !%-----------------------------------------------------------------------------
        !% pointer allocation
        fullDepth => elemR(StorageIdx,er_FullDepth)
        CurveID   => elemSI(StorageIdx,esi_JunctionMain_Curve_ID)
        aConst    => elemSR(StorageIdx,esr_Storage_Constant)
        aCoeff    => elemSR(StorageIdx,esr_Storage_Coefficient)  
        aExpon    => elemSR(StorageIdx,esr_Storage_Exponent)
        
        !% num of rows (HACK: read from a setting file)
        nRow =>  setting%Junction%FunStorageN

        !% add a new curveID
        CurveID  = SWMM_N_Curve + oneI

        curve(CurveID)%ID       = CurveID
        curve(CurveID)%Type     = StorageCurve
        curve(CurveID)%NumRows  = nRow 
        curve(CurveID)%ElemIdx  = StorageIdx

        !% allocate the valueArray for the new curve
        call util_allocate_curve_entries (CurveID, nRow)

        Curve(CurveID)%ValueArray(:,curve_storage_depth)  = util_linspace(zeroR,fullDepth,nRow)
        Curve(CurveID)%ValueArray(:,curve_storage_area)   = aConst + aCoeff * &
                    Curve(CurveID)%ValueArray(:,curve_storage_depth) ** aExpon
        Curve(CurveID)%ValueArray(:,curve_storage_volume) =  aConst * Curve(CurveID)%ValueArray(:,curve_storage_depth) + &
                    (aCoeff/(aExpon+oneR)) * Curve(CurveID)%ValueArray(:,curve_storage_depth) ** (aExpon+oneR)

    end subroutine storage_create_curve
!%  
!%==========================================================================

!%==========================================================================
!%
!%
!     real(8) function trapezoidal_area_from_depth_singular (indx) result (outvalue)
!         !%-----------------------------------------------------------------------------
!         !% Description:
!         !% Computes area from known depth for trapezoidal cross section of a single element
!         !% The input indx is the row index in full data 2D array.
!         !%-----------------------------------------------------------------------------
!         integer, intent(in) :: indx
!         real(8), pointer :: depth(:), breadth(:), lslope(:), rslope(:)
!         !%-----------------------------------------------------------------------------
!         depth   => elemR(:,er_Depth)
!         breadth => elemSGR(:,esgr_Trapezoidal_Breadth)
!         lslope  => elemSGR(:,esgr_Trapezoidal_LeftSlope)
!         rslope  => elemSGR(:,esgr_Trapezoidal_RightSlope)
!         !%-----------------------------------------------------------------------------
!         outvalue = (breadth(indx) + onehalfR * (lslope(indx) + rslope(indx)) * depth(indx)) * depth(indx)

!     end function trapezoidal_area_from_depth_singular

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
end module storage_geometry