module trapezoidal_channel
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Geometry for trapezoidal open channel
    !%==========================================================================
    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use geometry_lowlevel

    implicit none

    private

    public :: trapezoidal_depth_from_volume
    public :: trapezoidal_topwidth_from_depth
    public :: trapezoidal_perimeter_from_depth

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine trapezoidal_depth_from_volume (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels 
        !% Input elemPGx is pointer (already assigned) for elemPGetm
        !% Assumes that volume > 0 is previuosly enforced in volume computations.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: depth(:), volume(:)
            real(8), pointer :: fullvolume(:), fulldepth(:)
        !%-------------------------------------------------------------------
        !% Aliases
            depth      => elemR(:,er_Depth)
            volume     => elemR(:,er_Volume)
            fullvolume => elemR(:,er_FullVolume)
            fulldepth  => elemR(:,er_FullDepth)
        !%------------------------------------------------------------------

        if (setting%Discretization%AllowChannelOverflowTF) then
            where (volume(thisP) >= fullvolume(thisP))
                !% --- truncate depth at full depth if there is overflow
                depth(thisP) = fulldepth(thisP)
            elsewhere
                !% --- standard trapezoidal depth
                depth(thisP) = llgeo_trapezoidal_depth_from_volume_pure &
                                    (thisP, volume(thisP))
            endwhere
        else
            !% --- if overflow is NOT allowed
            where (volume(thisP) >= fullvolume(thisP))
                !% --- volume above max level is rectangular at max breadth
                depth(thisP) = llgeo_openchannel_depth_above_full_pure(thisP)
            elsewhere
                !% --- standard trapezoidal depth
                depth(thisP) = llgeo_trapezoidal_depth_from_volume_pure &
                                 (thisP, volume(thisP))
            end where
        end if
 
    end subroutine trapezoidal_depth_from_volume
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine trapezoidal_topwidth_from_depth (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a trapezoidal channel
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: topwidth(:), volume(:), fullvolume(:)
            real(8), pointer :: depth(:)
        !%------------------------------------------------------------------
        !% Aliases    
            topwidth    => elemR(:,er_Topwidth)
            volume      => elemR(:,er_Volume)
            fullvolume  => elemR(:,er_FullVolume)
            depth       => elemR(:,er_Depth)
        !%----------------------------------------------------------------
        if (setting%Discretization%AllowChannelOverflowTF) then 
            !% --- depth is already limited to full depth
            topwidth(thisP) = llgeo_trapezoidal_topwidth_from_depth_pure(thisP,depth(thisP))
        else
            where (volume(thisP) >= fullvolume(thisP))
                topwidth(thisP) = llgeo_openchannel_topwidth_above_full_pure(thisP)
            elsewhere
                !% --- use elemental form as depth <= fulldepth
                topwidth(thisP) = llgeo_trapezoidal_topwidth_from_depth_pure(thisP,depth(thisP))
            endwhere

        end if

    end subroutine trapezoidal_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine trapezoidal_perimeter_from_depth (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a trapezoidal channel
        !%------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: perimeter(:), volume(:), fullvolume(:)
            real(8), pointer :: depth(:)
        !%------------------------------------------------------------------
        !% Aliases
            perimeter  => elemR(:,er_Perimeter)
            volume     => elemR(:,er_Volume)
            fullvolume => elemR(:,er_FullVolume)
            depth      => elemR(:,er_Depth)
        !%-----------------------------------------------------------------

        if (setting%Discretization%AllowChannelOverflowTF) then 
            !% --- depth is already limited to full depth
            perimeter(thisP) = llgeo_trapezoidal_perimeter_from_depth_pure(thisP,depth(thisP))
        else
            where (volume(thisP) >= fullvolume(thisP))
                perimeter(thisP) = llgeo_openchannel_perimeter_above_full_pure(thisP)
            elsewhere
                !% --- use elemental form as depth <= fulldepth
                perimeter(thisP) = llgeo_trapezoidal_perimeter_from_depth_pure(thisP,depth(thisP))
            endwhere
        end if
        

    end subroutine trapezoidal_perimeter_from_depth
!%    
!%==========================================================================
!% END MODULE
!%==========================================================================
end module trapezoidal_channel