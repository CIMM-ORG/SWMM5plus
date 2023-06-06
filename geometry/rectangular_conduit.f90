module rectangular_conduit
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Geometry for rectangular closed conduit
    !%==========================================================================
    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys

    implicit none

    private

    public :: rectangular_closed_depth_from_volume
    public :: rectangular_closed_topwidth_from_depth
    public :: rectangular_closed_perimeter_from_depth

    contains

!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine rectangular_closed_depth_from_volume (thisP)
        !%-------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels 
        !% Input elemPGx is pointer (already assigned) for elemPGetm
        !% Assumes that volume > 0 is enforced in volume computations.
        !%--------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: depth(:), volume(:), length(:), breadth(:)
            real(8), pointer :: fulldepth(:), fullvolume(:)
        !%--------------------------------------------------------------------
        !% Aliases
            depth       => elemR(:,er_Depth)
            volume      => elemR(:,er_Volume)
            fulldepth   => elemR(:,er_FullDepth)
            fullvolume  => elemR(:,er_FullVolume)
            length      => elemR(:,er_Length)
            breadth     => elemSGR(:,esgr_Rectangular_Breadth)
        !%---------------------------------------------------------------------  
        
        depth(thisP) = volume(thisP) / (length(thisP) * breadth(thisP))

        !% --- ensure the full depth is not exceeded
        depth(thisP) = min(depth(thisP),fulldepth(thisP))

    end subroutine rectangular_closed_depth_from_volume
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_closed_topwidth_from_depth (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a rectangular channel
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: thisP(:)
            integer, pointer :: GeomType(:)
            real(8), pointer :: breadth(:), topwidth(:), depth(:), fullDepth(:)
        !%-------------------------------------------------------------------
        !% Aliases:
            GeomType  => elemI(:,ei_geometryType)
            topwidth  => elemR(:,er_Topwidth)
            depth     => elemR(:,er_Depth)
            fullDepth => elemR(:,er_FullDepth)
            breadth   => elemSGR(:,esgr_Rectangular_Breadth)
        !%-------------------------------------------------------------------

        topwidth(thisP) = breadth(thisP)

        !% reset topwidth to zeroValue for depth higher than full depth
        where (depth(thisP) >= fullDepth(thisP))
               topwidth(thisP) = setting%ZeroValue%Topwidth
        end where

    end subroutine rectangular_closed_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_closed_perimeter_from_depth (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a rectangular channel
        !%------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: breadth(:), depth(:), perimeter(:), fulldepth(:), fullPerimeter(:)
        !%------------------------------------------------------------------
        !% Aliases:
            breadth   => elemSGR(:,esgr_Rectangular_Breadth)
            depth     => elemR(:,er_Depth)
            fulldepth => elemR(:,er_FullDepth)
            perimeter => elemR(:,er_Perimeter)
            fullPerimeter => elemR(:,er_FullPerimeter)
        !%-------------------------------------------------------------------
        
        where (depth(thisP) < fulldepth(thisP))
            perimeter(thisP) = twoR * depth(thisP) + breadth(thisP) 
        else where (depth(thisP) >= fulldepth(thisP))
            perimeter(thisP) = fullPerimeter(thisP)
        end where

    end subroutine rectangular_closed_perimeter_from_depth
!%    
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module rectangular_conduit