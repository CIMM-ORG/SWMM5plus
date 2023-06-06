module rectangular_channel
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Geometry for rectangular open channel
    !%==========================================================================
    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use geometry_lowlevel

    implicit none

    private

    public :: rectangular_depth_from_volume
    public :: rectangular_topwidth_from_depth    
    public :: rectangular_perimeter_from_depth

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine rectangular_depth_from_volume (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels 
        !% Input elemPGx is pointer (already assigned) for elemPGetm
        !% Assumes that volume > 0 is previously enforced in volume computations.
        !%------------------------------------------------------------------
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: depth(:), volume(:)
            real(8), pointer :: fulldepth(:), fullvolume(:)
        !%------------------------------------------------------------------  
        !% Aliases:  
            depth      => elemR(:,er_Depth)
            volume     => elemR(:,er_Volume)
            fulldepth  => elemR(:,er_FullDepth)
            fullvolume => elemR(:,er_FullVolume)
        !%----------------------------------------------------------------- 

        if (setting%Discretization%AllowChannelOverflowTF) then
            where (volume(thisP) >= fullvolume(thisP))
                !% --- truncate depth at full depth if there is overflow
                depth(thisP) = fulldepth(thisP)
            elsewhere
                !% --- standard rectangular depth
                depth(thisP) = llgeo_rectangular_depth_from_volume_pure &
                                    (thisP, volume(thisP))
            endwhere
        else 
            !% --- if overflow is NOT allowed
            !% --- standard rectangular depth if overflow is suppressed
            depth(thisP) = llgeo_rectangular_depth_from_volume_pure &
                                    (thisP, volume(thisP))
        end if

    end subroutine rectangular_depth_from_volume
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_topwidth_from_depth (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a rectangular channel
        !%------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            real(8), pointer ::  topwidth(:), depth(:)
        !%------------------------------------------------------------------
        !% Aliases
            topwidth  => elemR(:,er_Topwidth)
            depth     => elemR(:,er_Depth)
        !%-----------------------------------------------------------------

        topwidth(thisP) =  llgeo_rectangular_topwidth_from_depth_pure(thisP,depth(thisP))

    end subroutine rectangular_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_perimeter_from_depth (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a rectangular channel
        !%------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: perimeter(:), depth(:)
        !%-------------------------------------------------------------------
        !% Aliases:
            perimeter => elemR(:,er_Perimeter)
            depth     => elemR(:,er_Depth)
        !%-------------------------------------------------------------------

        perimeter(thisP) = llgeo_rectangular_perimeter_from_depth_pure(thisP,depth(thisP))

    end subroutine rectangular_perimeter_from_depth
!%    
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module rectangular_channel