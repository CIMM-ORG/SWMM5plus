module storage_geometry

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Trapezoidal channel geometry
    !%

    private

    public :: functional_storage_depth_from_volume

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine functional_storage_depth_from_volume (elemPGx, Npack, thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Only applies on JM that has storage (or non-surcharged trapezoidal conduits)
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:), volume(:)
        real(8), pointer :: aConst(:), aCoeff(:), aExpon(:)
        !%-----------------------------------------------------------------------------
        thisP   => elemPGx(1:Npack,thisCol) 
        depth   => elemR(:,er_Depth)
        volume  => elemR(:,er_Volume)
        aConst  => elemSR(:,esr_Storage_Constant)
        aCoeff  => elemSR(:,esr_Storage_Coefficient)
        aExpon  => elemSR(:,esr_Storage_Exponent)
            

        
    end subroutine functional_storage_depth_from_volume
! !%  
! !%==========================================================================
! !%==========================================================================
! !%
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