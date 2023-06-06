module parabolic_channel
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Geometry for parabolic open channel
    !%==========================================================================
    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use geometry_lowlevel

    implicit none

    private

    public :: parabolic_depth_from_volume
    public :: parabolic_topwidth_from_depth
    public :: parabolic_perimeter_from_depth

    contains

!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine parabolic_depth_from_volume (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels
        !% Input elemPGx is pointer (already assigned) for elemPGetm 
        !% Assumes that volume > 0 is previously enforced in volume computations.
        !%--------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: depth(:), volume(:)
            real(8), pointer :: fulldepth(:), fullvolume(:)
        !%-------------------------------------------------------------------
        !% Aliases
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
                !% --- standard parabolic depth
                depth(thisP) = llgeo_parabolic_depth_from_volume_pure &
                                (thisP, volume(thisP))
            endwhere
        else
            !% --- if overflow is NOT allowed
            where (volume(thisP) >= fullvolume(thisP))
                !% --- volume above max level is rectangular at max breadth
                depth(thisP) = llgeo_openchannel_depth_above_full_pure(thisP)
            elsewhere 
                !% --- standard parabolic depth
                depth(thisP) = llgeo_parabolic_depth_from_volume_pure &
                                (thisP, volume(thisP))
            endwhere
        end if
        
    end subroutine parabolic_depth_from_volume
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine parabolic_topwidth_from_depth (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a parabolic channel
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) ::  thisP(:)
            real(8), pointer :: topwidth(:), volume(:), fullvolume(:)
            real(8), pointer :: depth(:)
        !%-------------------------------------------------------------------
            topwidth   => elemR(:,er_Topwidth)
            volume     => elemR(:,er_Volume)
            fullvolume => elemR(:,er_FullVolume)
            depth      => elemR(:,er_Depth)
        !%-------------------------------------------------------------------

        if (setting%Discretization%AllowChannelOverflowTF) then 
            !% --- depth is already limited to full depth
            topwidth(thisP) = llgeo_parabolic_topwidth_from_depth_pure(thisP,depth(thisP))
        else
            where (volume(thisP) >= fullvolume(thisP))
                topwidth(thisP) = llgeo_openchannel_topwidth_above_full_pure(thisP)
            elsewhere
                !% --- use elemental form as depth <= fulldepth
                topwidth(thisP) = llgeo_parabolic_topwidth_from_depth_pure(thisP,depth(thisP))
            endwhere
        end if
        
    end subroutine parabolic_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine parabolic_perimeter_from_depth (thisP)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a parabolic channel
        !% NOTE THIS DOES NOT FOLLOW THE PATTERN OF OTHER SUBROUTINES AS IT
        !% DOES NOT USE ELEMENTAL FUNCTIONS
        !%-----------------------------------------------------------------------------
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: depth(:), perimeter(:), tempA(:)
            real(8), pointer :: volume(:), fullvolume(:)
        !%-----------------------------------------------------------------------------
            depth      => elemR(:,er_Depth)
            perimeter  => elemR(:,er_Perimeter)
            volume     => elemR(:,er_Volume)
            fullvolume => elemR(:,er_FullVolume)
            tempA      => elemR(:,er_Temp01)
        !%-----------------------------------------------------------------------------
        
        perimeter(thisP) = llgeo_parabolic_perimeter_from_depth_pure (thisP, depth(thisP))

        tempA(thisP) = zeroR

        !% --- correct for overfull channel
        !%     This could be made simpler if the function call for perimeter_from_depth
        !%     is elemental.
        if (.not. setting%Discretization%AllowChannelOverflowTF) then
            !% --- compute the corrected for all
            tempA(thisP) = llgeo_openchannel_perimeter_above_full_pure(thisP) 
            !%
            where (volume(thisP) > fullvolume(thisP))
                perimeter(thisP) = tempA(thisP)
            endwhere
        else 
            !% no changes needed if channel overflow is allowed
        endif

    end subroutine parabolic_perimeter_from_depth
!%    
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module parabolic_channel