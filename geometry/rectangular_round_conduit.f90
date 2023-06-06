module rectangular_round_conduit
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Geometry for rectangular round closed conduit
    !%==========================================================================
    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use define_xsect_tables
    use xsect_tables
    use circular_conduit, only: circular_get_normalized_depth_from_area_analytical

    implicit none

    private

    public :: rect_round_depth_from_volume
    public :: rect_round_topwidth_from_depth
    public :: rect_round_perimeter_from_depth

    contains

!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine rect_round_depth_from_volume (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on closed conduit of rect_round cross section
        !% Input elemPGx is pointer (already assigned) for elemPGetm 
        !% Assumes that volume > 0 is enforced in volume computations.
        !%-------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: depth(:), fullArea(:), fullDepth(:)
            real(8), pointer :: length(:), breadth(:), AoverAfull(:), YoverYfull(:)
            real(8), pointer :: volume(:), aBot(:), rBot(:), yBot(:), pi
            logical, pointer :: botSection(:)
            integer, allocatable, target :: thisP_analytical(:), thisP_lookup(:)
            integer, target              :: Npack_analytical, Npack_lookup
        !%--------------------------------------------------------------------
        !% Aliases
            depth       => elemR(:,er_Depth)
            fullDepth   => elemR(:,er_FullDepth)
            fullArea    => elemR(:,er_FullArea) 
            breadth     => elemR(:,er_BreadthMax)
            aBot        => elemSGR(:,esgr_Rectangular_Round_Abot)
            rBot        => elemSGR(:,esgr_Rectangular_Round_Rbot)
            yBot        => elemSGR(:,esgr_Rectangular_Round_Ybot)
            volume      => elemR(:,er_Volume)
            length      => elemR(:,er_Length)
            AoverAfull  => elemR(:,er_Temp01)
            YoverYfull  => elemR(:,er_Temp02)
            botSection  => elemYN(:,eYN_Temp01)
            pi          => setting%Constant%pi
        !%----------------------------------------------------------------------

        !% --- initialize AoverAfull
        AoverAfull(thisP) = zeroR

        !% --- bottom circular section
        where (volume(thisP) > (aBot(thisP) * length(thisP)))
            depth(thisP) = yBot(thisP) + (volume(thisP) / length(thisP) - aBot(thisP)) / breadth(thisP)
            botSection(thisP) = .false.

        !% --- top rectangular part
        elsewhere
            !% --- find unfilled top-area/area of full circular top
            AoverAfull(thisP) = (volume(thisP) / length(thisP)) / (pi * rBot(thisP) ** twoR)
            botSection(thisP) = .true.
        end where

        !% --- pack where the circular top with AoverAfull <= 4% which will use analytical solution
        !%     from French, 1985 by using the central angle theta.
        !%     HACK -- this needs to be replaced with temporary storage rather than dynamic allocation
        Npack_analytical = count((AoverAfull(thisP) <= 0.04) .and. botSection(thisP))
        thisP_analytical = pack(thisP,(AoverAfull(thisP) <= 0.04 .and. botSection(thisP)))

        !% --- pack where the rest of the elements having AoverAfull > 0.04 which will use
        !%     lookup table for interpolation.
        !%     HACK -- this needs to be replaced with temporary storage rather than dynamic allocation
        Npack_lookup = count((AoverAfull(thisP) > 0.04) .and. botSection(thisP))
        thisP_lookup = pack(thisP,(AoverAfull(thisP) > 0.04 .and. botSection(thisP)))

        if (Npack_analytical > zeroI) then
            call circular_get_normalized_depth_from_area_analytical &
                (YoverYfull, AoverAfull, Npack_analytical, thisP_analytical)
        end if 

        if (Npack_lookup > zeroI) then        
            !% --- retrive the normalized Y/Yfull from the lookup table
            call xsect_table_lookup &
                (YoverYfull, AoverAfull, YCirc, thisP_lookup)  
        end if

        !% --- finally get the depth 
        where (botSection(thisP))
            depth(thisP) = twoR * rBot(thisP) * YoverYfull(thisP)
        end where

        !% --- ensure the full depth is not exceeded
        depth(thisP) = min(depth(thisP),fulldepth(thisP))
                
    end subroutine rect_round_depth_from_volume
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rect_round_topwidth_from_depth (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a rectangular_round channel
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: depth(:), topwidth(:), breadth(:), yBot(:), rBot(:)
        !%-------------------------------------------------------------------
        !% Aliases
            topwidth    => elemR(:,er_Topwidth)
            depth       => elemR(:,er_Depth)
            breadth     => elemR(:,er_BreadthMax)
            yBot        => elemSGR(:,esgr_Rectangular_Round_Ybot)
            rBot        => elemSGR(:,esgr_Rectangular_Round_Ybot) 
            
        !%--------------------------------------------------------------------

        where (depth(thisP) <= zeroR)
            topwidth(thisP) = setting%ZeroValue%Topwidth
        elsewhere (depth(thisP) > yBot(thisP))
            !% --- top rectangular section
            topwidth(thisP) = breadth(thisP)
        elsewhere (depth(thisP) <= yBot(thisP)) 
            !% --- bottom circular section
            topwidth(thisP) = twoR * sqrt(depth(thisP) * (twoR * rBot(thisP) - depth(thisP)))
        end where

    end subroutine rect_round_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine rect_round_perimeter_from_depth (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a rectangular_round channel
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: perimeter(:), breadth(:), depth(:)
            real(8), pointer :: fullDepth(:), theta(:), yBot(:), rBot(:)
            real(8), pointer :: fullPerimeter(:)
        !%------------------------------------------------------------------
        !% Aliases:
            depth        => elemR(:,er_Depth)
            perimeter    => elemR(:,er_Perimeter)
            fullDepth    => elemR(:,er_FullDepth)
            fullPerimeter=> elemR(:,er_FullPerimeter)
            breadth      => elemR(:,er_BreadthMax)
            yBot         => elemSGR(:,esgr_Rectangular_Round_Ybot)
            rBot         => elemSGR(:,esgr_Rectangular_Round_Rbot)  
            theta        => elemR(:,er_Temp01)
        !%-------------------------------------------------------------------

        where (depth(thisP) <= setting%ZeroValue%Depth)
            !% --- water level effectively zero depth
            perimeter(thisP) = setting%ZeroValue%Topwidth
        elsewhere ((depth(thisP) > setting%ZeroValue%Depth) .and.  (depth(thisP) <= yBot(thisP)))
            !% --- water level in bottom circular section
            theta(thisP)     = twoR * acos(oneR - depth(thisP) / rBot(thisP))
            perimeter(thisP) = rBot(thisP) * theta(thisP)
        elsewhere ((depth(thisP) > yBot(thisP)) .and. (depth(thisp) < fullDepth(thisP)))
            !% --- water level in top rectangular section
            theta(thisP)     = twoR * asin(breadth(thisP) / twoR / rBot(thisP))
            perimeter(thisP) = rBot(thisP) * theta(thisP)                                  &
                             + twoR * (min(depth(thisP), fullDepth(thisP)) - yBot(thisP))
        elsewhere
            !% --- if none of the above are met, then this must be full
            perimeter(thisP) = fullPerimeter(thisP)
        end where

    end subroutine rect_round_perimeter_from_depth
!%    
!%==========================================================================
!% END MODULE
!%==========================================================================
end module rectangular_round_conduit