module filled_circular_conduit
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Geometry for filled circular closed conduit
    !%
    !%==========================================================================
    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use define_xsect_tables
    use circular_conduit, only: circular_get_normalized_depth_from_area_analytical
    use xsect_tables

    implicit none

    private

    public :: filled_circular_depth_from_volume
    public :: filled_circular_topwidth_from_depth
    public :: filled_circular_hydradius_and_perimeter_from_depth

    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine filled_circular_depth_from_volume (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on conduits 
        !% Input elemPGx is pointer (already assigned) for elemPGetm 
        !% Assumes that volume > 0 is enforced in volume computations.
        !%-------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: depth(:), volume(:), length(:), AoverAfull(:)
            real(8), pointer :: YoverYfull(:), fullArea(:), fullDepth(:)
            real(8), pointer :: Ybottom(:), Abottom(:), tempYfull(:), tempAfull(:)

            integer, allocatable, target :: thisP_analytical(:), thisP_lookup(:)
            integer, target              :: Npack_analytical, Npack_lookup
            integer :: ii
        !%---------------------------------------------------------------------
        !% Aliases
            depth      => elemR(:,er_Depth)
            volume     => elemR(:,er_Volume)
            length     => elemR(:,er_Length)
            fullArea   => elemR(:,er_FullArea)
            fullDepth  => elemR(:,er_fullDepth)
            Ybottom    => elemR(:,er_SedimentDepth)
            Abottom    => elemSGR(:,esgr_Filled_Circular_bottomArea)
            AoverAfull => elemR(:,er_AoverAfull)
            YoverYfull => elemR(:,er_YoverYfull)
            tempYfull  => elemR(:,er_Temp01)
            tempAfull  => elemR(:,er_Temp02)
        !%----------------------------------------------------------------------

        !% --- calculate a temporary geometry by considering the whole cicrular cross-section
        tempYfull(thisP) = fullDepth(thisP) + Ybottom(thisP)
        tempAfull(thisP) = fullArea(thisP)  + Abottom(thisP)

        !% --- compute the relative volume
        AoverAfull(thisP) = (volume(thisP) / length(thisP) + Abottom(thisP)) / tempAfull(thisP)

        !% --- pack the filled circular elements with AoverAfull <= 4% which will use analytical solution
        !%     from French, 1985 by using the central angle theta.
        !% HACK -- this needs to be replaced with temporary storage rather than dynamic allocation
        Npack_analytical = count(AoverAfull(thisP) <= 0.04)
        thisP_analytical = pack(thisP,AoverAfull(thisP) <= 0.04)

        !% --- pack the rest of the filled circular elements having AoverAfull > 0.04 which will use
        !%     lookup table for interpolation.
        !% HACK -- this needs to be replaced with temporary storage rather than dynamic allocation
        Npack_lookup = count(AoverAfull(thisP) > 0.04)
        thisP_lookup = pack(thisP,AoverAfull(thisP) > 0.04)

        if (Npack_analytical > zeroI) then
            call circular_get_normalized_depth_from_area_analytical &
                (YoverYfull, AoverAfull, Npack_analytical, thisP_analytical)
        end if 

        if (Npack_lookup > zeroI) then        
            !% --- retrive the normalized Y/Yfull from the lookup table
            call xsect_table_lookup &
                (YoverYfull, AoverAfull, YCirc, thisP_lookup) 
        end if

        !% --- inally get the depth by multiplying the normalized depth with full depth
        !%     and substract the bottom depth
        depth(thisP) = YoverYfull(thisP) * tempYfull(thisP) - Ybottom(thisP)

        !% --- ensure the full depth is not exceeded
        depth(thisP) = min(depth(thisP),fulldepth(thisP))
        
    end subroutine filled_circular_depth_from_volume
!%
!%========================================================================== 
!%==========================================================================
!%
    subroutine filled_circular_topwidth_from_depth (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a filled_circular conduit
        !%-------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: depth(:), topwidth(:), YoverYfull(:)
            real(8), pointer :: fullDepth(:), Ybottom(:), tempYfull(:) 
            real(8), pointer :: breadthMax(:)
        !%-------------------------------------------------------------------
        !% Aliases
            depth      => elemR(:,er_Depth)
            topwidth   => elemR(:,er_Topwidth)
            fullDepth  => elemR(:,er_fullDepth)
            breadthMax => elemR(:,er_BreadthMax)
            Ybottom    => elemR(:,er_SedimentDepth)
            YoverYfull => elemR(:,er_YoverYfull)
            tempYfull  => elemR(:,er_Temp01)
        !%---------------------------------------------------------------------

        !% ---calculate a temporary depth by considering the whole circular
        !%    cross-section
        tempYfull(thisP) = fullDepth(thisP) + Ybottom(thisP)

        !% --- get the normalized depth
        YoverYfull(thisP) = (depth(thisP) + Ybottom(thisP)) / tempYfull(thisP)

        !% --- retrive the normalized T/Tmax from the lookup table for circular pipe
        !%     T/Tmax value is temporarily saved in the topwidth column
        call xsect_table_lookup &
            (topwidth, YoverYfull, TCirc, thisP)      

        !% --- get the topwidth by multiplying the T/Tmax with max breadth
        topwidth(thisP) = max(topwidth(thisP) * breadthMax(thisP), setting%ZeroValue%Topwidth)

    end subroutine filled_circular_topwidth_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine filled_circular_hydradius_and_perimeter_from_depth (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the hydraulic radius from a known depth in a filled_circular conduit
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: depth(:), area(:), hydRadius(:), perimeter(:)
            real(8), pointer :: tempArea(:)
            real(8), pointer :: totalFullDepth(:), totalFullArea(:), totalFullHydRadius(:)
            real(8), pointer :: sedimentDepth(:), sedimentPerimeter(:), sedimentTopwidth(:)
            real(8), pointer :: YoverYfull(:) 
        !%-------------------------------------------------------------------
        !% Aliases
            depth              => elemR(:,er_Depth)
            area               => elemR(:,er_Area)
            hydRadius          => elemR(:,er_HydRadius)
            perimeter          => elemR(:,er_Perimeter)
            sedimentDepth      => elemR(:,er_SedimentDepth)
            tempArea           => elemR(:,er_Temp01)
            totalFullDepth     => elemSGR(:,esgr_Filled_Circular_TotalPipeDiameter)
            totalFullArea      => elemSGR(:,esgr_Filled_Circular_TotalPipeArea)
            totalFullHydRadius => elemSGR(:,esgr_Filled_Circular_TotalPipeHydRadius)
            sedimentPerimeter  => elemSGR(:,esgr_Filled_Circular_bottomPerimeter)
            sedimentTopwidth     => elemSGR(:,esgr_Filled_Circular_bottomTopwidth)
            YoverYfull         => elemR(:,er_YoverYfull)  
        !%-----------------------------------------------------------------------------

        !% --- normalized depth based on entire circular cross-section
        YoverYfull(thisP) = (depth(thisP) + sedimentDepth(thisP)) &
                            / totalFullDepth(thisP)

        !% --- prevent overfull
        YoverYfull(thisP) = min(YoverYfull(thisP), oneR)

        !% --- prevent underfull (must be above sediment depth)
        YoverYfull(thisP) = max(YoverYfull(thisP), &
            (sedimentDepth(thisP) + setting%ZeroValue%Depth)/totalFullDepth )

        !% --- retrieve the normalized R/Rmax for the entire circular cross
        !%     section from the lookup table
        !%     R/Rmax value is temporarily saved in the hydRadius column
        !%     we can use this as temporary since we will overwrite it before
        !%     we are finished
        call xsect_table_lookup (hydRadius, YoverYfull, RCirc, thisP)  

        !% --- ensure no zeros in normalized hydraulic radius
        hydRadius(thisP) = max(hydRadius(thisP),setting%ZeroValue%Depth)     

        !% --- retrieve the normalized A/Afull for the entire circular cross
        !%     section from the lookup table. This must use a temp storage as
        !%     we don't want to overwrite the stored area
        call xsect_table_lookup (tempArea, YoverYfull, ACirc, thisP) 

        !% --- compute the perimeter of the flow section
        !%     get the total perimeter from hydradius and area of the flow and filled
        !%     section (must be unnormalized), then subtract the bottom perimeter 
        !%     and add the bottom topwidth
        perimeter(thisP) = (   ( totalFullArea(thisP)      * tempArea(thisP) )  &
                             / ( totalFullHydRadius(thisp) * hydRadius(thisP))  &
                           ) - sedimentPerimeter(thisP) + sedimentTopwidth(thisP)

        !% --- ensure no zeros in perimeter
        perimeter(thisP) = max(perimeter(thisP),setting%ZeroValue%Topwidth)                         

        !% --- compute hydraulic radius
        hydRadius(thisP) = area(thisP) / perimeter(thisP)

        !% --- clear temporary
        tempArea(thisP) = nullvalueR
     
    end subroutine filled_circular_hydradius_and_perimeter_from_depth
!%
!%=========================================================================
!% END OF MODULE
!%=========================================================================
end module filled_circular_conduit