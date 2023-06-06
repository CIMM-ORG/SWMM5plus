module mod_basket_conduit
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Geometry for modified basket conduit
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

    public :: mod_basket_depth_from_volume
    public :: mod_basket_topwidth_from_depth
    public :: mod_basket_perimeter_from_depth

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine mod_basket_depth_from_volume (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels 
        !% Input elemPGx is pointer (already assigned) for elemPGetm
        !% Assumes that volume > 0 is enforced in volume computations.
        !%-------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            integer, pointer :: thisP_analytical(:), thisP_lookup(:)
            integer, pointer :: thisPA(:), thisPL(:)
            real(8), pointer :: depth(:), fullArea(:), fullDepth(:), topArea(:), rTop(:)
            real(8), pointer :: length(:), breadth(:), AoverAfull(:), YoverYfull(:)
            real(8), pointer :: volume(:), pi
            logical, pointer :: topSection(:)
            integer, target              :: Npack_analytical, Npack_lookup
        !%---------------------------------------------------------------------
            depth       => elemR(:,er_Depth)
            fullDepth   => elemR(:,er_FullDepth)
            fullArea    => elemR(:,er_FullArea) 
            breadth     => elemR(:,er_BreadthMax)
            topArea     => elemSGR(:,esgr_Mod_Basket_Atop)
            rTop        => elemSGR(:,esgr_Mod_Basket_Rtop)
            volume      => elemR(:,er_Volume)
            length      => elemR(:,er_Length)
            AoverAfull  => elemR(:,er_AoverAfull)
            YoverYfull  => elemR(:,er_YoverYfull)
            topSection  => elemYN(:,eYN_Temp01)

            thisP_analytical => elemI(:,ei_Temp01)
            thisP_lookup     => elemI(:,ei_Temp02)
            pi               => setting%Constant%pi
        !%-----------------------------------------------------------------------------

        !% --- initialize AoverAfull
        AoverAfull(thisP) = zeroR 

        !% --- bottom rectangular section
        where(volume(thisP) <= (fullArea(thisP) - topArea(thisP)) * length(thisP))
            depth(thisP) = volume(thisP) / (length(thisP) * breadth(thisP))
            topSection(thisP) = .false.

        !% --- top circular part
        elsewhere
            !% --- find unfilled top-area/area of full circular top
            AoverAfull(thisP) = (fullArea(thisP) - volume(thisP) / length(thisP)) / (pi * rTop(thisP) ** twoR)
            topSection(thisP) = .true.
        end where

        !% --- pack where the circular top with AoverAfull <= 4% which will use analytical solution
        !%     from French, 1985 by using the central angle theta.
        !%     HACK -- this needs to be replaced with temporary storage rather than dynamic allocation
        Npack_analytical = count((AoverAfull(thisP) <= 0.04) .and. topSection(thisP))
        if (Npack_analytical > zeroI) then

            thisP_analytical = pack(thisP,(AoverAfull(thisP) <= 0.04 .and. topSection(thisP)))
            thisPA => thisP_analytical(1:Npack_analytical)

            call circular_get_normalized_depth_from_area_analytical &
                (YoverYfull, AoverAfull, Npack_analytical, thisPA)
        end if 

        !% --- pack where the rest of the elements having AoverAfull > 0.04 which will use
        !%     lookup table for interpolation.
        !%     HACK -- this needs to be replaced with temporary storage rather than dynamic allocation
        Npack_lookup = count((AoverAfull(thisP) > 0.04) .and. topSection(thisP))
        
        if (Npack_lookup > zeroI) then  

            thisP_lookup = pack(thisP,(AoverAfull(thisP) > 0.04 .and. topSection(thisP)))
            thisPL => thisP_lookup(1:Npack_lookup)       
            !% --- retrive the normalized Y/Yfull from the lookup table
            call xsect_table_lookup &
                (YoverYfull, AoverAfull, YCirc, thisPL)  
        end if

        !% --- get the depth by calculating the difference between full height & unfilled heigh
        where (topSection(thisP))
            depth(thisP) = fullDepth(thisP) - twoR * rTop(thisP) * YoverYfull(thisP)
        end where

        !% --- ensure the full depth is not exceeded
        depth(thisP) = min(depth(thisP),fulldepth(thisP))

        !% --- clear the temporary storage
        if (Npack_analytical > zeroI) thisPA = nullvalueI
        if (Npack_lookup     > zeroI) thisPL = nullvalueI
                
    end subroutine mod_basket_depth_from_volume
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine mod_basket_topwidth_from_depth (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a mod_basket channel
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: breadth(:), topwidth(:), fullDepth(:), depth(:)
            real(8), pointer :: yBreadthMax(:), rTop(:), emptyDepth(:)
        !%-------------------------------------------------------------------
            topwidth    => elemR(:,er_Topwidth)
            depth       => elemR(:,er_Depth)
            fullDepth   => elemR(:,er_FullDepth)
            yBreadthMax => elemR(:,er_DepthAtBreadthMax)
            breadth     => elemR(:,er_BreadthMax)
            rTop        => elemSGR(:,esgr_Mod_Basket_Rtop) 
            emptyDepth  => elemR(:,er_Temp01)
        !%------------------------------------------------------------------

        where(depth(thisP) <= zeroR)
            topwidth(thisP) = setting%ZeroValue%Topwidth

        elsewhere(depth(thisP) <= yBreadthMax(thisP))
            topwidth(thisP) = breadth(thisP)

        elsewhere (depth(thisP) > yBreadthMax(thisP))
            emptyDepth(thisP) = max(fullDepth(thisP) - depth(thisP), zeroR)
            topwidth(thisP)   = twoR * sqrt(emptyDepth(thisP) * (twoR * rTop(thisP) - emptyDepth(thisP)))
        end where

    end subroutine mod_basket_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine mod_basket_perimeter_from_depth (thisP) 
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a mod_basket channel
        !%------------------------------------------------------------------
        !% Declarastions
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: depth(:),fullDepth(:), perimeter(:)
            real(8), pointer :: breadth(:), yBreadthMax(:), rTop(:), thetaTop(:)
            real(8), pointer :: emptyDepth(:), emptyTheta(:)
        !%------------------------------------------------------------------
        !% Aliases:
            perimeter   => elemR(:,er_Perimeter)
            depth       => elemR(:,er_Depth)
            fullDepth   => elemR(:,er_FullDepth)
            yBreadthMax => elemR(:,er_DepthAtBreadthMax)
            breadth     => elemR(:,er_BreadthMax)
            rTop        => elemSGR(:,esgr_Mod_Basket_Rtop)
            thetaTop    => elemSGR(:,esgr_Mod_Basket_ThetaTop) 
            emptyDepth  => elemR(:,er_Temp01)
            emptyTheta  => elemR(:,er_Temp02)
        !%-----------------------------------------------------------------------------

        where(depth(thisP) <= yBreadthMax(thisP))
            perimeter(thisP) = twoR * depth(thisP) + breadth(thisP) 
        elsewhere
            !% --- find height of empty area
            emptyDepth(thisP) = max(fullDepth(thisP) - depth(thisP), zeroR)
            !% --- find angle of circular arc corresponding to this height
            emptyTheta(thisP) = twoR * acos(oneR - emptyDepth(thisP) / rTop(thisP))
            !% --- find perimeter of wetted portion of circular arc
            perimeter(thisP)  = (thetaTop(thisP) - emptyTheta(thisP)) * rTop(thisP)
            !% --- add on wetted perimeter of bottom rectangular area
            perimeter(thisP)  = perimeter(thisP) + twoR * yBreadthMax(thisP) + breadth(thisP) 
        end where

    end subroutine mod_basket_perimeter_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
end module mod_basket_conduit