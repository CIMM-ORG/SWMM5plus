module rectangular_triangular_conduit
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Geometry for rectangular triangular closed conduit
    !%==========================================================================
    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys

    implicit none

    private

    public :: rectangular_triangular_depth_from_volume
    public :: rectangular_triangular_topwidth_from_depth
    public :: rectangular_triangular_perimeter_from_depth
    
    contains

!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine rectangular_triangular_depth_from_volume (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on rectangular_triangular closed conduit
        !% Input elemPGx is pointer (already assigned) for elemPGetm 
        !% Assumes that volume > 0 is enforced in volume computations.
        !% NOTE: this does NOT limit the depth by surcharge height at this point
        !% This will be done after the head is computed.
        !%-------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: depth(:), bottomDepth(:), bottomArea(:)
            real(8), pointer :: fulldepth(:)
            real(8), pointer :: volume(:), length(:), breadth(:), bottomSlope(:)

        !%-------------------------------------------------------------------
        !% Aliases
            depth       => elemR(:,er_Depth)
            fulldepth   => elemR(:,er_FullDepth)
            breadth     => elemR(:,er_FullTopWidth)
            bottomDepth => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth) 
            bottomArea  => elemSGR(:,esgr_Rectangular_Triangular_BottomArea)
            bottomSlope => elemSGR(:,esgr_Rectangular_Triangular_BottomSlope)
            volume      => elemR(:,er_Volume)
            length      => elemR(:,er_Length)
        !%----------------------------------------------------------------------

        where(volume(thisP) <= bottomArea(thisP) * length(thisP))
            depth(thisP) = sqrt(volume(thisP) / (length(thisP) * bottomSlope(thisP)))

        elsewhere(volume(thisP) > bottomarea(thisP)*length(thisP))
            depth(thisP) = bottomDepth(thisP) &
                + ((volume(thisP) / length(thisP)) - bottomarea(thisP)) / breadth(thisP)
        end where

        !% --- ensure the full depth is not exceeded
        depth(thisP) = min(depth(thisP),fulldepth(thisP))
                
    end subroutine rectangular_triangular_depth_from_volume
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_triangular_topwidth_from_depth (thisP)  
        !%-------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a rectangular_triangular channel
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: breadth(:), topwidth(:), bottomSlope(:), depth(:), bottomDepth(:)
        !%--------------------------------------------------------------------
        !% Aliases
            topwidth    => elemR(:,er_Topwidth)
            depth       => elemR(:,er_Depth)
            breadth     => elemR(:,er_FullTopWidth)
            bottomSlope => elemSGR(:,esgr_Rectangular_Triangular_BottomSlope)
            bottomDepth => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth) 
        !%---------------------------------------------------------------------

        where(depth(thisP) <= bottomDepth(thisP))
            topwidth(thisP) = twoR * bottomSlope(thisP) * depth(thisP)

        else where(depth(thisP) > bottomDepth(thisP))
            topwidth(thisP) = breadth(thisP)
        endwhere

    end subroutine rectangular_triangular_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_triangular_perimeter_from_depth (thisP)
        !%-------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a rectangular_triangular channel
        !%-------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: breadth(:), depth(:), bottomSlope(:), bottomDepth(:)
            real(8), pointer :: perimeter(:), fullPerimeter(:), fullDepth(:)
        !%-------------------------------------------------------------------
        !% Aliases:
            breadth       => elemR(:,er_FullTopWidth)
            depth         => elemR(:,er_Depth)
            bottomSlope   => elemSGR(:,esgr_Rectangular_Triangular_BottomSlope)
            bottomDepth   => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth)  
            perimeter     => elemR(:,er_Perimeter)
            fullperimeter => elemR(:,er_FullPerimeter)
            fulldepth     => elemR(:,er_FullDepth)
        !%--------------------------------------------------------------------

        where (depth(thisP) <= setting%ZeroValue%Depth)
            !% --- negligible water level
            perimeter(thisP) = setting%ZeroValue%Topwidth

        elsewhere ( (depth(thisP) > setting%ZeroValue%Depth) .and. (depth(thisP) <= bottomdepth(thisP)) )
            !% --- water level in lower triangular portion
            perimeter(thisP) = twoR * depth(thisP) * sqrt(oneR + bottomSlope(thisP) ** twoR)

        elsewhere ( (depth(thisP) > bottomdepth(thisP)) .and. (depth(thisP) < fulldepth(thisP)) )
            !% --- water level in upper rectangular portion
            perimeter(thisP) = twoR * bottomdepth(thisP) * sqrt(oneR + bottomSlope(thisP) ** twoR)  & !triangular section
                                 + twoR * (depth(thisP) - bottomDepth(thisP))                         !rectangular section
        elsewhere
            !% --- water level full
            perimeter(thisP) = fullperimeter(thisP)

        end where

    end subroutine rectangular_triangular_perimeter_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
end module rectangular_triangular_conduit