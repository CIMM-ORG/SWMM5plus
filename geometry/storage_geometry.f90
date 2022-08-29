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
    public :: storage_implied_depth_from_volume
    public :: storage_integrate_area_vs_depth_curve
    public :: storage_interpolate_volume_from_depth_singular
    public :: storage_create_curve
    public :: storage_volume_from_depth_singular
    public :: storage_implied_volume_from_depth_singular
    public :: storage_functional_volume_from_depth_singular
    public :: storage_tabular_volume_from_depth_singular
    !public :: storage_implied_length

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
        real(8), pointer :: aConst, aCoeff, aExpon, bb, cc
        real(8) :: aa
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


        ! print *, 'in storage functional depth from volume'
        ! print *, volume


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
                        cc => volume
                        !% --- quadratic formula, Note that the b^2 - 4ac is written
                        !%     as b^2 + 4ac because c (the volume) is negative in the
                        !%     derivation but the pointer cc is to a positive value.
                        !%     This is guaranteed to be positive as long as aCoeff,
                        !%     aExpon, and volume are positive
                        depth = (- bb + sqrt(bb**2 + fourR * aa * cc) ) / (twoR * aa)
                    else            
                        !% --- interpolate from the curve created in initial_condition/init_IC_get_junction_data
                        !      using storage_create_curve()
                        call util_curve_lookup_singular(curveID, er_Volume, er_Depth, curve_storage_volume, &
                            curve_storage_depth, 1)
                    end if
                end if
            endif

            ! print *, volume

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

        ! print *, ' '
        ! print *, 'in storage tabular depth from volume'
        ! print *, ' '
        
        !% HACK: Find out a way to code this without do loop
        do ii = 1, Npack
            thisP   => elemPGx(ii,thisCol) 
            depth   => elemR(thisP,er_Depth)
            volume  => elemR(thisP,er_Volume)
            curveID => elemSI(thisP,esi_JunctionMain_Curve_ID)

            !% interpolate from the curve
            call util_curve_lookup_singular(curveID, er_Volume, er_Depth, curve_storage_volume, &
                curve_storage_depth,1)
        end do

    end subroutine storage_tabular_depth_from_volume
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine storage_implied_depth_from_volume (elemPGx, Npack, thisCol)
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
    end subroutine storage_implied_depth_from_volume
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
                integrated_volume(ii) = integrated_volume(ii-oneI) + onehalfR * (depth(ii) - depth(ii-1)) * &
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
        !% interpolate from the storage curve to get volume (stored in elemR)
        !% from the depth (stored in elemR)
        !%-----------------------------------------------------------------------------
        integer, intent(in)  :: StorageIdx

        integer, pointer :: curveID

        character(64) :: subroutine_name = 'storage_interpolate_from_curve'
        !%-----------------------------------------------------------------------------
        curveID  => elemSI(StorageIdx,esi_JunctionMain_Curve_ID)

        call util_curve_lookup_singular(curveID, er_Depth, er_Volume, curve_storage_depth, &
            curve_storage_volume, 1)

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
        N_FunctionalStorage = N_FunctionalStorage + oneI
        CurveID  = setting%SWMMinput%N_curve + N_FunctionalStorage
      
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

        select case(elemSI(idx,esi_JunctionMain_Type))
        case (ImpliedStorage)
            outvalue = storage_implied_volume_from_depth_singular(idx, indepth)
        case (TabularStorage)
            outvalue = storage_tabular_volume_from_depth_singular(idx, indepth)
        case (FunctionalStorage)
            outvalue = storage_functional_volume_from_depth_singular(idx, indepth)
        case default
            print *, 'CODE ERROR: unexpected case default for storage '
            print *, 'problem in (elemSI(idx,esi_JunctionMain_Type) '
            print *, 'element index (idx)= ',idx
            print *, 'junction main type of ',(elemSI(idx,esi_JunctionMain_Type))
            print *, trim(reverseKey(elemSI(idx,esi_JunctionMain_Type)))
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
        !%-------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: idx
            real(8), intent(in) :: indepth
        !%-------------------------------------------------------------------
        outvalue = elemSR(idx,esr_Storage_Plane_Area) * indepth
            
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
!%
    ! subroutine storage_implied_length (elemPGx, Npack, thisCol)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Computes an implied length scale of the storage element
    !     !% for tabular or functional storage elements
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
    !         integer, pointer :: thisP(:)
    !         real(8), pointer :: depth(:), volume(:), length(:)
    !     !%-----------------------------------------------------------------------------
    !     !% Aliases
    !         if (Npack < 1) return
    !         thisP   => elemPGx(1:Npack,thisCol) 
    !         depth   => elemR(:,er_Depth)
    !         volume  => elemR(:,er_Volume)
    !         length  => elemR(:,er_Length)
    !     !%------------------------------------------------------------------
        
    !     length(thisP) = sqrt(volume(thisP) / depth(thisP))    
        
    ! end subroutine storage_implied_length
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