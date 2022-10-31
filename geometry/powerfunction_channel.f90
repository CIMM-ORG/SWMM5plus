module powerfunction_channel

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use geometry_lowlevel

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% powerfunction channel geometry
    !%

    private

    public :: powerfunction_depth_from_volume
    public :: powerfunction_topwidth_from_depth
    public :: powerfunction_perimeter_from_depth

    ! public :: powerfunction_area_from_depth_singular
    ! public :: powerfunction_topwidth_from_depth_singular 
    ! public :: powerfunction_perimeter_from_depth_singular

    contains

!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine powerfunction_depth_from_volume (elemPGx, Npack, thisCol)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !%--------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
            integer, pointer :: thisP(:)
            real(8), pointer :: depth(:), volume(:)
            real(8), pointer :: fulldepth(:), fullvolume(:)
        !%-------------------------------------------------------------------
            if (Npack < 1) return
        !%-------------------------------------------------------------------
            thisP      => elemPGx(1:Npack,thisCol) 
            depth      => elemR(:,er_Depth)
            volume     => elemR(:,er_Volume)   
            fulldepth  => elemR(:,er_FullDepth)
            fullvolume => elemR(:,er_FullVolume)
        !%--------------------------------------------------------------------

        if (setting%Discretization%AllowChannelOverflowTF) then     
            where (volume(thisP) >= fullvolume(thisP))
                !% --- truncate depth at full depth if there is overflow
                depth(thisP) = fulldepth(thisP)
            elsewhere
                !% --- standard powerfunction depth
                depth(thisP) = llgeo_powerfunction_depth_from_volume_pure &
                                (thisP, volume(thisP))
            endwhere
        else
            !% --- if overflow is NOT allowed
            where (volume(thisP) >= fullvolume(thisP))
                !% --- volume above max level is rectangular at max breadth
                depth(thisP) = llgeo_openchannel_depth_above_full_pure(thisP)
            elsewhere 
                !% --- standard powerfunction depth
                depth(thisP) = llgeo_powerfunction_depth_from_volume_pure &
                                (thisP, volume(thisP))
            endwhere
        end if
        
    end subroutine powerfunction_depth_from_volume
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine powerfunction_topwidth_from_depth (elemPGx, Npack, thisCol)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a powerfunction channel
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: elemPGx(:,:)
            integer, intent(in) ::  Npack, thisCol
            integer, pointer :: thisP(:)
            real(8), pointer :: topwidth(:), volume(:), fullvolume(:)
            real(8), pointer :: depth(:)
        !%-------------------------------------------------------------------
            thisP      => elemPGx(1:Npack,thisCol) 
            topwidth   => elemR(:,er_Topwidth)
            depth      => elemR(:,er_Depth)
            volume     => elemR(:,er_Volume)
            fullvolume => elemR(:,er_FullVolume)
        !%-------------------------------------------------------------------

        if (setting%Discretization%AllowChannelOverflowTF) then 
            !% --- depth is already limited to full depth
            topwidth(thisP) = llgeo_powerfunction_topwidth_from_depth_pure(thisP,depth(thisP))
        else
            where (volume(thisP) >= fullvolume(thisP))
                topwidth(thisP) = llgeo_openchannel_topwidth_above_full_pure(thisP)
            elsewhere
                !% --- use elemental form as depth <= fulldepth
                topwidth(thisP) = llgeo_powerfunction_topwidth_from_depth_pure(thisP,depth(thisP))
            endwhere
        end if
        
    end subroutine powerfunction_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine powerfunction_perimeter_from_depth (elemPGx, Npack, thisCol)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a powerfunction channel
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: elemPGx(:,:)
            integer, intent(in) ::  Npack, thisCol
            integer, pointer :: thisP(:)
            real(8), pointer :: perimeter(:), volume(:), fullvolume(:)
            real(8), pointer :: depth(:)
        !%-------------------------------------------------------------------
            thisP      => elemPGx(1:Npack,thisCol) 
            perimeter  => elemR(:,er_Perimeter)
            volume     => elemR(:,er_Volume)
            fullvolume => elemR(:,er_FullVolume)
            depth      => elemR(:,er_Depth)
        !%-------------------------------------------------------------------

        if (setting%Discretization%AllowChannelOverflowTF) then 
            !% --- depth is already limited to full depth
            perimeter(thisP) = llgeo_powerfunction_perimeter_from_depth_pure(thisP,depth(thisP))
        else
            where (volume(thisP) >= fullvolume(thisP))
                perimeter(thisP) = llgeo_openchannel_perimeter_above_full_pure(thisP)
            elsewhere
                !% --- use elemental form as depth <= fulldepth
                perimeter(thisP) = llgeo_powerfunction_perimeter_from_depth_pure(thisP,depth(thisP))
            endwhere
        end if

    end subroutine powerfunction_perimeter_from_depth
!%    
!%==========================================================================
!% SINGULAR
!%==========================================================================
!%    
!     pure real(8) function powerfunction_area_from_depth_singular &
!         (indx, depth) result (outvalue)
!         !%------------------------------------------------------------------
!         !% Description:
!         !% Computes area from known depth for powerfunction cross section
!         !% of a single element.
!         !% The input indx is the row index in full data 2D array.
!         !%------------------------------------------------------------------
!         !% Declarations
!             integer, intent(in) :: indx
!             real(8), intent(in) :: depth
!         !%------------------------------------------------------------------

!         outvalue = nullvalueR !% STUB

!     end function powerfunction_area_from_depth_singular 
! !%    
! !%==========================================================================
! !%==========================================================================
! !%
!     pure real(8) function powerfunction_topwidth_from_depth_singular &
!         (indx, depth) result (outvalue)
!         !%------------------------------------------------------------------
!         !% Description:
!         !% Computes topwidth from known depth for powerfunction cross section
!         !% of a single element.
!         !% The input indx is the row index in full data 2D array.
!         !%------------------------------------------------------------------
!         !% Declarations
!             integer, intent(in) :: indx
!             real(8), intent(in) :: depth
!      !%------------------------------------------------------------------

!         outvalue = nullvalueR !% STUB

!     end function powerfunction_topwidth_from_depth_singular 
! !%    
! !%==========================================================================
! !%==========================================================================
! !%        
!     pure real(8) function powerfunction_perimeter_from_depth_singular &
!         (indx, depth) result (outvalue)
!         !%------------------------------------------------------------------
!         !% Description:
!         !% Computes perimeter from known depth for powerfunction cross section
!         !% of a single element.
!         !% The input indx is the row index in full data 2D array.
!         !%------------------------------------------------------------------
!         !% Declarations
!             integer, intent(in) :: indx
!             real(8), intent(in) :: depth
!         !%------------------------------------------------------------------

!         outvalue = nullvalueR !% STUB

!     end function powerfunction_perimeter_from_depth_singular 
!%
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module powerfunction_channel