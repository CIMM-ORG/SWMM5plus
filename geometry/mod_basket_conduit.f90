module mod_basket_conduit

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use define_xsect_tables
    use xsect_tables
    use circular_conduit, only: circular_get_normalized_depth_from_area_analytical

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% mod_basket channel geometry
    !%

    private

    public :: mod_basket_depth_from_volume
    public :: mod_basket_area_from_depth
    public :: mod_basket_area_from_depth_singular
    public :: mod_basket_topwidth_from_depth
    public :: mod_basket_topwidth_from_depth_singular 
    public :: mod_basket_perimeter_from_depth
    public :: mod_basket_perimeter_from_depth_singular
    public :: mod_basket_hyddepth_from_topwidth
    public :: mod_basket_hyddepth_from_topwidth_singular
    public :: mod_basket_hydradius_from_depth_singular

    contains

!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine mod_basket_depth_from_volume (elemPGx, Npack, thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels 
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !% NOTE: this does NOT limit the depth by surcharge height at this point
        !% This will be done after the head is computed.
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:), fullArea(:), fullDepth(:), topArea(:), rTop(:)
        real(8), pointer :: length(:), breadth(:), AoverAfull(:), YoverYfull(:)
        real(8), pointer :: volume(:), pi
        logical, pointer :: topSection(:)
        integer, allocatable, target :: thisP_analytical(:), thisP_lookup(:)
        integer, target              :: Npack_analytical, Npack_lookup
        !%-----------------------------------------------------------------------------
        thisP       => elemPGx(1:Npack,thisCol) 
        depth       => elemR(:,er_Depth)
        fullDepth   => elemR(:,er_FullDepth)
        fullArea    => elemR(:,er_FullArea) 
        topArea     => elemSGR(:,esgr_Mod_Basket_Atop)
        breadth     => elemSGR(:,esgr_Mod_Basket_BreadthMax)
        rTop        => elemSGR(:,esgr_Mod_Basket_Rtop)
        volume      => elemR(:,er_Volume)
        length      => elemR(:,er_Length)
        AoverAfull  => elemR(:,er_Temp01)
        YoverYfull  => elemR(:,er_Temp02)
        topSection  => elemYN(:,eYN_Temp01)
        pi          => setting%Constant%pi
        !%-----------------------------------------------------------------------------
        !% initialize AoverAfull
        AoverAfull(thisP) = zeroR 
        !% bottom rectangular section
        where(volume(thisP) <= (fullArea(thisP) - topArea(thisP)) * length(thisP))
            depth(thisP) = sqrt(volume(thisP) / (length(thisP) * breadth(thisP)))
            topSection(thisP) = .false.
        !% top circular part
        elsewhere
            !% find unfilled top-area/area of full circular top
            AoverAfull(thisP) = (fullArea(thisP) - volume(thisP) / length(thisP)) / (pi * rTop(thisP) ** twoR)
            topSection(thisP) = .true.
        end where

        !% --- pack where the circular top with AoverAfull <= 4% which will use analytical solution
        !%     from French, 1985 by using the central angle theta.
        !% HACK -- this needs to be replaced with temporary storage rather than dynamic allocation
        Npack_analytical = count((AoverAfull(thisP) <= 0.04) .and. topSection(thisP))
        thisP_analytical = pack(thisP,(AoverAfull(thisP) <= 0.04 .and. topSection(thisP)))

        !% --- pack where the rest of the elements having AoverAfull > 0.04 which will use
        !%     lookup table for interpolation.
        !% HACK -- this needs to be replaced with temporary storage rather than dynamic allocation
        Npack_lookup = count((AoverAfull(thisP) > 0.04) .and. topSection(thisP))
        thisP_lookup = pack(thisP,(AoverAfull(thisP) > 0.04 .and. topSection(thisP)))

        if (Npack_analytical > zeroI) then
            call circular_get_normalized_depth_from_area_analytical &
                (YoverYfull, AoverAfull, Npack_analytical, thisP_analytical)
        end if 

        if (Npack_lookup > zeroI) then        
            !% retrive the normalized Y/Yfull from the lookup table
            call xsect_table_lookup &
                (YoverYfull, AoverAfull, YCirc, thisP_lookup)  
        end if

        !% finally get the depth by calculating the difference between full height & unfilled heigh
        where (topSection(thisP))
            depth(thisP) = fullDepth(thisP) - twoR * rTop(thisP) * YoverYfull(thisP)
        end where
                
    end subroutine mod_basket_depth_from_volume
!%
!%==========================================================================
!%==========================================================================
!%
    elemental real(8) function mod_basket_area_from_depth (indx) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for rectangular cross section
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx  ! may be a packed array of indexes
        real(8) :: emptyDepth, emptyTheta, emptyArea
        !%-----------------------------------------------------------------------------
        
        if (elemR(indx,er_Depth) <= elemSGR(indx,esgr_Mod_Basket_YatMaxBreadth)) then
            outvalue = elemR(indx,er_Depth) * elemSGR(indx,esgr_Basket_Handle_BreadthMax)
        else
            emptyDepth = max(elemR(indx,er_FullDepth) - elemR(indx,er_Depth), zeroR)
            emptyTheta = twoR * acos(oneR - emptyDepth / elemSGR(indx,esgr_Mod_Basket_Rtop))
            emptyArea  = onehalfR * (elemSGR(indx,esgr_Mod_Basket_Rtop) ** twoR) * (emptyTheta - sin(emptyTheta))
            outvalue   = elemR(indx,er_FullArea) - emptyArea
        endif

    end function mod_basket_area_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function mod_basket_area_from_depth_singular (indx, depth) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for mod_basket cross section of a single element
        !% The input indx is the row index in full data 2D array.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer :: fulldepth(:), fullArea(:)
        real(8), pointer :: yBreadthMax(:), breadth(:), rTop(:)
        real(8) :: emptyDepth, emptyTheta, emptyArea
        !%-----------------------------------------------------------------------------
        fullArea    => elemR(:,er_FullArea)
        fulldepth   => elemR(:,er_FullDepth)
        yBreadthMax => elemSGR(:,esgr_Mod_Basket_YatMaxBreadth) 
        breadth     => elemSGR(:,esgr_Basket_Handle_BreadthMax)
        rTop        => elemSGR(:,esgr_Mod_Basket_Rtop)
        !%-----------------------------------------------------------------------------
        
        if(depth <= yBreadthMax(indx)) then
            outvalue = depth * breadth(indx)
        else
            emptyDepth = max(fulldepth(indx) - depth, zeroR) 
            emptyTheta = twoR * acos(oneR - emptyDepth / rTop(indx))
            emptyArea  = onehalfR * (rTop(indx) ** twoR) * (emptyTheta - sin(emptyTheta))
            outvalue   = fullArea(indx) - emptyArea  
        endif

    end function mod_basket_area_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine mod_basket_topwidth_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a mod_basket channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: breadth(:), topwidth(:), fullDepth(:), depth(:)
        real(8), pointer :: yBreadthMax(:), rTop(:), emptyDepth(:)
        !%-----------------------------------------------------------------------------
        thisP       => elemPGx(1:Npack,thisCol) 
        topwidth    => elemR(:,er_Topwidth)
        depth       => elemR(:,er_Depth)
        fullDepth   => elemR(:,er_FullDepth)
        yBreadthMax => elemSGR(:,esgr_Mod_Basket_YatMaxBreadth)
        breadth     => elemSGR(:,esgr_Mod_Basket_BreadthMax)
        rTop        => elemSGR(:,esgr_Mod_Basket_Rtop) 
        emptyDepth  => elemR(:,er_Temp01)
        !%-----------------------------------------------------------------------------

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
    real(8) function mod_basket_topwidth_from_depth_singular (indx, depth) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth for a mod_basket cross section of a single element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx 
        real(8), intent(in) :: depth
        real(8), pointer :: breadth(:), topwidth(:), fullDepth(:)
        real(8), pointer :: yBreadthMax(:), rTop(:)
        real(8) :: emptyDepth
        !%-----------------------------------------------------------------------------
        topwidth    => elemR(:,er_Topwidth)
        fullDepth   => elemR(:,er_FullDepth)
        yBreadthMax => elemSGR(:,esgr_Mod_Basket_YatMaxBreadth)
        breadth     => elemSGR(:,esgr_Mod_Basket_BreadthMax)
        rTop        => elemSGR(:,esgr_Mod_Basket_Rtop) 
        !%-----------------------------------------------------------------------------
         
        if (depth <= zeroR) then
            outvalue = setting%ZeroValue%Topwidth
        else if (depth <= yBreadthMax(indx)) then
            outvalue = breadth(indx)
        else
            emptyDepth = max(fullDepth(indx) - depth, zeroR)
            outvalue   = twoR * sqrt(emptyDepth * (twoR * rTop(indx) - emptyDepth))
        endif

    end function mod_basket_topwidth_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine mod_basket_perimeter_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a mod_basket channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:),fullDepth(:), perimeter(:)
        real(8), pointer :: breadth(:), yBreadthMax(:), rTop(:), thetaTop(:)
        real(8), pointer :: emptyDepth(:), emptyTheta(:)
        !%-----------------------------------------------------------------------------
        thisP       => elemPGx(1:Npack,thisCol)  
        perimeter   => elemR(:,er_Perimeter)
        fullDepth   => elemR(:,er_FullDepth)
        yBreadthMax => elemSGR(:,esgr_Mod_Basket_YatMaxBreadth)
        breadth     => elemSGR(:,esgr_Mod_Basket_BreadthMax)
        rTop        => elemSGR(:,esgr_Mod_Basket_Rtop)
        thetaTop    => elemSGR(:,esgr_Mod_Basket_ThetaTop) 
        emptyDepth  => elemR(:,er_Temp01)
        emptyTheta  => elemR(:,er_Temp02)
        !%-----------------------------------------------------------------------------

        where(depth(thisP) <= yBreadthMax(thisP))
            perimeter(thisP) = twoR * depth(thisP) + breadth(thisP) 
        elsewhere
            !% find height of empty area
            emptyDepth(thisP) = max(fullDepth(thisP) - depth(thisP), zeroR)
            !% find angle of circular arc corresponding to this height
            emptyTheta(thisP) = twoR * acos(oneR - emptyDepth(thisP) / rTop(thisP))
            !% find perimeter of wetted portion of circular arc
            perimeter(thisP)  = (thetaTop(thisP) - emptyTheta(thisP)) * rTop(thisP)
            !% add on wetted perimeter of bottom rectangular area
            perimeter(thisP)  = perimeter(thisP) + twoR * yBreadthMax(thisP) + breadth(thisP) 
        end where

    end subroutine mod_basket_perimeter_from_depth
!%    
!%==========================================================================    
!%==========================================================================
!%
    real(8) function mod_basket_perimeter_from_depth_singular (indx, depth) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes wetted perimeter from known depth for a mod_basket cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer :: fullDepth(:), breadth(:), yBreadthMax(:), rTop(:), thetaTop(:)
        real(8) :: emptyDepth, emptyTheta
        !%-----------------------------------------------------------------------------
        fullDepth   => elemR(:,er_FullDepth)
        yBreadthMax => elemSGR(:,esgr_Mod_Basket_YatMaxBreadth)
        breadth     => elemSGR(:,esgr_Mod_Basket_BreadthMax)
        rTop        => elemSGR(:,esgr_Mod_Basket_Rtop)
        thetaTop    => elemSGR(:,esgr_Mod_Basket_ThetaTop) 
        !%-----------------------------------------------------------------------------
        
        if(depth <= yBreadthMax(indx)) then
            outvalue = twoR * depth + breadth(indx) 
        else
            !% find height of empty area
            emptyDepth = max(fullDepth(indx) - depth, zeroR)
            !% find angle of circular arc corresponding to this height
            emptyTheta = twoR * acos(oneR - emptyDepth / rTop(indx))
            !% find perimeter of wetted portion of circular arc
            outvalue  = (thetaTop(indx) - emptyTheta) * rTop(indx)
            !% add on wetted perimeter of bottom rectangular area
            outvalue  = outvalue + twoR * yBreadthMax(indx) + breadth(indx)                     
        endif

    end function mod_basket_perimeter_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine mod_basket_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the hydraulic (average) depth from a known depth in a mod_basket channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: hyddepth(:), depth(:), area(:), topwidth(:), fullHydDepth(:)
        !%-----------------------------------------------------------------------------
        thisP       => elemPGx(1:Npack,thisCol) 
        depth       => elemR(:,er_Depth)
        area        => elemR(:,er_Area)
        topwidth    => elemR(:,er_Topwidth)
        hyddepth    => elemR(:,er_HydDepth)
        fullHydDepth => elemR(:,er_FullHydDepth)
        !%-----------------------------------------------------------------------------

        !% when conduit is empty
        where (depth(thisP) <= setting%ZeroValue%Depth)
            hyddepth(thisP) = setting%ZeroValue%Depth

        !% when conduit is not empty
        elsewhere (depth(thisP) > setting%ZeroValue%Depth)
            !% limiter for when the conduit is full
            hyddepth(thisP) = min(area(thisP) / topwidth(thisP), fullHydDepth(thisP))
        endwhere

    end subroutine mod_basket_hyddepth_from_topwidth
!%    
!%==========================================================================  
!%==========================================================================
!%
    real(8) function mod_basket_hyddepth_from_topwidth_singular (indx,topwidth,depth) result (outvalue)
    
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic depth from known depth for mod_basket cross section of 
        !% a single element
        !%-----------------------------------------------------------------------------   
        integer, intent(in) :: indx     
        real(8), intent(in) :: depth, topwidth
        real(8), pointer    :: area(:), fullHydDepth(:)
        !%-----------------------------------------------------------------------------
        area         => elemR(:,er_Area)
        fullHydDepth => elemR(:,er_FullHydDepth)
        !%--------------------------------------------------

        !% calculating hydraulic depth needs conditional since,
        !% topwidth can be zero in cross section for both
        !% full and empty condition.

        !% when conduit is empty
        if (depth <= setting%ZeroValue%Depth) then
            outvalue = setting%ZeroValue%Depth
        else
            !% limiter for when the conduit is full
            outvalue = min(area(indx) / topwidth, fullHydDepth(indx))
        endif

    end function mod_basket_hyddepth_from_topwidth_singular 
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function mod_basket_hydradius_from_depth_singular (indx, depth) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic radius from known depth for a mod_basket cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer :: fullDepth(:), fullArea(:) 
        real(8), pointer :: breadth(:), yBreadthMax(:), rTop(:), thetaTop(:)
        real(8) :: emptyDepth, emptyTheta, emptyArea, Perimeter, Area
        !%-----------------------------------------------------------------------------
        fullDepth   => elemR(:,er_FullDepth)
        yBreadthMax => elemSGR(:,esgr_Mod_Basket_YatMaxBreadth)
        breadth     => elemSGR(:,esgr_Mod_Basket_BreadthMax)
        rTop        => elemSGR(:,esgr_Mod_Basket_Rtop)
        thetaTop    => elemSGR(:,esgr_Mod_Basket_ThetaTop) 
        !%-----------------------------------------------------------------------------
        
        if(depth <= yBreadthMax(indx)) then
            outvalue = (depth * breadth(indx)) / (twoR * depth + breadth(indx) )
        else
            !% find height of empty area
            emptyDepth = max(fullDepth(indx) - depth, zeroR)
            !% find angle of circular arc corresponding to this height
            emptyTheta = twoR * acos(oneR - emptyDepth / rTop(indx))
            !% find the empty area
            emptyArea  = onehalfR * (rTop(indx) ** twoR) * (emptyTheta - sin(emptyTheta));
            !% find perimeter of wetted portion of circular arc
            Perimeter  = (thetaTop(indx) - emptyTheta) * rTop(indx)
            !% add on wetted perimeter of bottom rectangular area
            Perimeter  = Perimeter + twoR * yBreadthMax(indx) + breadth(indx) 
            !% find the area
            Area = fullArea(indx) - emptyArea 
            !% hydraulic radius
            outvalue = Area / Perimeter
        endif

    end function mod_basket_hydradius_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
end module mod_basket_conduit