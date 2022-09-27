module rectangular_round_conduit

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
    !% rectangular_round channel geometry
    !%

    private

    public :: rect_round_depth_from_volume
    public :: rect_round_area_from_depth
    public :: rect_round_area_from_depth_singular
    public :: rect_round_topwidth_from_depth
    public :: rect_round_topwidth_from_depth_singular 
    public :: rect_round_perimeter_from_depth
    public :: rect_round_perimeter_from_depth_singular
    public :: rect_round_hyddepth_from_topwidth
    public :: rect_round_hyddepth_from_topwidth_singular
    public :: rect_round_hydradius_from_depth_singular

    contains

!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine rect_round_depth_from_volume (elemPGx, Npack, thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels (or non-surcharged rectangular_round conduits)
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !% NOTE: this does NOT limit the depth by surcharge height at this point
        !% This will be done after the head is computed.
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:), fullArea(:), fullDepth(:)
        real(8), pointer :: length(:), breadth(:), AoverAfull(:), YoverYfull(:)
        real(8), pointer :: volume(:), aBot(:), rBot(:), yBot(:), pi
        logical, pointer :: botSection(:)
        integer, allocatable, target :: thisP_analytical(:), thisP_lookup(:)
        integer, target              :: Npack_analytical, Npack_lookup
        !%-----------------------------------------------------------------------------
        thisP       => elemPGx(1:Npack,thisCol) 
        depth       => elemR(:,er_Depth)
        fullDepth   => elemR(:,er_FullDepth)
        fullArea    => elemR(:,er_FullArea) 
        breadth     => elemSGR(:,esgr_Rectangular_Round_BreadthMax)
        aBot        => elemSGR(:,esgr_Rectangular_Round_Abot)
        rBot        => elemSGR(:,esgr_Rectangular_Round_Rbot)
        yBot        => elemSGR(:,esgr_Rectangular_Round_Ybot)
        volume      => elemR(:,er_Volume)
        length      => elemR(:,er_Length)
        AoverAfull  => elemR(:,er_Temp01)
        YoverYfull  => elemR(:,er_Temp02)
        botSection  => elemYN(:,eYN_Temp01)
        pi          => setting%Constant%pi
        !%-----------------------------------------------------------------------------
        !% initialize AoverAfull
        AoverAfull(thisP) = zeroR
        !% bottom circular section
        where (volume(thisP) > (aBot(thisP) * length(thisP)))
            depth(thisP) = yBot(thisP) + (volume(thisP) / length(thisP) - aBot(thisP)) / breadth(thisP)
            botSection(thisP) = .false.
        !% top rectangular part
        elsewhere
            !% find unfilled top-area/area of full circular top
            AoverAfull(thisP) = (volume(thisP) / length(thisP)) / (pi * rBot(thisP) ** twoR)
            botSection(thisP) = .true.
        end where

        !% --- pack where the circular top with AoverAfull <= 4% which will use analytical solution
        !%     from French, 1985 by using the central angle theta.
        !% HACK -- this needs to be replaced with temporary storage rather than dynamic allocation
        Npack_analytical = count((AoverAfull(thisP) <= 0.04) .and. botSection(thisP))
        thisP_analytical = pack(thisP,(AoverAfull(thisP) <= 0.04 .and. botSection(thisP)))

        !% --- pack where the rest of the elements having AoverAfull > 0.04 which will use
        !%     lookup table for interpolation.
        !% HACK -- this needs to be replaced with temporary storage rather than dynamic allocation
        Npack_lookup = count((AoverAfull(thisP) > 0.04) .and. botSection(thisP))
        thisP_lookup = pack(thisP,(AoverAfull(thisP) > 0.04 .and. botSection(thisP)))

        if (Npack_analytical > zeroI) then
            call circular_get_normalized_depth_from_area_analytical &
                (YoverYfull, AoverAfull, Npack_analytical, thisP_analytical)
        end if 

        if (Npack_lookup > zeroI) then        
            !% retrive the normalized Y/Yfull from the lookup table
            call xsect_table_lookup &
                (YoverYfull, AoverAfull, YCirc, thisP_lookup)  
        end if

        !% finally get the depth 
        where (botSection(thisP))
            depth(thisP) = twoR * rBot(thisP) * YoverYfull(thisP)
        end where
                
    end subroutine rect_round_depth_from_volume
!%
!%==========================================================================
!%==========================================================================
!%
    elemental real(8) function rect_round_area_from_depth (indx) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for rectangular cross section
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx  ! may be a packed array of indexes
        real :: theta
        !%-----------------------------------------------------------------------------
        
        if (elemR(indx,er_Depth) > elemSGR(indx,esgr_Rectangular_Round_Ybot)) then
            outvalue = elemSGR(indx,esgr_Rectangular_Round_Abot) + (elemR(indx,er_Depth) &
                     - elemSGR(indx,esgr_Rectangular_Round_Ybot)) * elemSGR(indx,esgr_Rectangular_Round_BreadthMax)
        else
            theta    = twoR * acos(oneR - elemR(indx,er_Depth) / elemSGR(indx,esgr_Rectangular_Round_Rbot))
            outvalue = onehalfR * (elemSGR(indx,esgr_Rectangular_Round_Rbot) ** twoR) * (theta - sin(theta))
        endif

    end function rect_round_area_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function rect_round_area_from_depth_singular (indx, depth) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for rectangular_round cross section of a single element
        !% The input indx is the row index in full data 2D array.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer :: yBot(:), rBot(:), aBot(:), breadth(:)
        real :: theta
        !%-----------------------------------------------------------------------------
        breadth => elemSGR(:,esgr_Rectangular_Round_BreadthMax)
        yBot    => elemSGR(:,esgr_Rectangular_Round_Ybot) 
        aBot    => elemSGR(:,esgr_Rectangular_Round_Abot)
        rBot    => elemSGR(:,esgr_Rectangular_Round_Rbot)
        !%-----------------------------------------------------------------------------
        if(depth > yBot(indx)) then
            outvalue = aBot(indx) + (depth - yBot(indx)) * breadth(indx)
        else
            theta    = twoR * acos(oneR - depth / rBot(indx))
            outvalue = onehalfR * (rBot(indx) ** twoR) * (theta - sin(theta))    
        endif

    end function rect_round_area_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rect_round_topwidth_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a rectangular_round channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:), topwidth(:), breadth(:), yBot(:), rBot(:)
        !%-----------------------------------------------------------------------------
        thisP       => elemPGx(1:Npack,thisCol) 
        topwidth    => elemR(:,er_Topwidth)
        depth       => elemR(:,er_Depth)
        yBot        => elemSGR(:,esgr_Rectangular_Round_Ybot)
        rBot        => elemSGR(:,esgr_Rectangular_Round_Ybot) 
        breadth     => elemSGR(:,esgr_Rectangular_Round_BreadthMax)
        !%-----------------------------------------------------------------------------

        where (depth(thisP) <= zeroR)
            topwidth(thisP) = setting%ZeroValue%Topwidth
        elsewhere (depth(thisP) > yBot(thisP))
            !% top rectangular section
            topwidth(thisP) = breadth(thisP)
        elsewhere (depth(thisP) <= yBot(thisP)) 
            !% bottom circular section
            topwidth(thisP) = twoR * sqrt(depth(thisP) * (twoR * rBot(thisP) - depth(thisP)))
        end where

    end subroutine rect_round_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function rect_round_topwidth_from_depth_singular (indx, depth) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth for a rectangular_round cross section of a single element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx 
        real(8), intent(in) :: depth
        real(8), pointer :: breadth(:), yBot(:), rBot(:)
        !%-----------------------------------------------------------------------------
        yBot        => elemSGR(:,esgr_Rectangular_Round_Ybot)
        rBot        => elemSGR(:,esgr_Rectangular_Round_Ybot) 
        breadth     => elemSGR(:,esgr_Rectangular_Round_BreadthMax)
        !%-----------------------------------------------------------------------------
         
        if (depth <= zeroR) then
            outvalue = setting%ZeroValue%Topwidth
        else if (depth > yBot(indx)) then
            !% top rectangular section
            outvalue = breadth(indx)
        else 
            !% bottom circular section
            outvalue = twoR * sqrt(depth * (twoR * rBot(indx) - depth))
        end if

    end function rect_round_topwidth_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rect_round_perimeter_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a rectangular_round channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: perimeter(:), breadth(:), depth(:)
        real(8), pointer :: fullDepth(:), theta(:), yBot(:), rBot(:)
        !%-----------------------------------------------------------------------------
        thisP       => elemPGx(1:Npack,thisCol) 
        depth       => elemR(:,er_Depth)
        perimeter   => elemR(:,er_Perimeter)
        fullDepth   => elemR(:,er_FullDepth)
        breadth     => elemSGR(:,esgr_Rectangular_Round_BreadthMax)
        yBot        => elemSGR(:,esgr_Rectangular_Round_Ybot)
        rBot        => elemSGR(:,esgr_Rectangular_Round_Rbot)  
        theta       => elemR(:,er_Temp01)
        !%-----------------------------------------------------------------------------

        where (depth(thisP) <= zeroR)
            perimeter(thisP) = zeroR

        elsewhere (depth(thisP) > yBot(thisP))
            !% top rectangular section
            theta(thisP)     = twoR * asin(breadth(thisP) / twoR / rBot(thisP))
            perimeter(thisP) = rBot(thisP) * theta(thisP) + twoR * (min(depth(thisP), &
                               fullDepth(thisP)) - yBot(thisP))

            where(depth(thisP) >= fullDepth(thisP))
                !% if the depth exceeds the full depth, the top part will contribute
                !% to overall wetted perimeter 
                perimeter(thisP) = perimeter(thisP) + breadth(thisP)
            end where

        elsewhere (depth(thisP) <= yBot(thisP))
            !% bottom circular section
            theta(thisP) = twoR * acos(oneR - depth(thisP) / rBot(thisP))
            perimeter(thisP) = rBot(thisP) * theta(thisP)
        end where

    end subroutine rect_round_perimeter_from_depth
!%    
!%==========================================================================    
!%==========================================================================
!%
    real(8) function rect_round_perimeter_from_depth_singular (indx, depth) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes wetted perimeter from known depth for a rectangular_round cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer :: breadth(:)
        real(8), pointer :: fullDepth(:), yBot(:), rBot(:)
        real(8) :: theta
        !%-----------------------------------------------------------------------------
        fullDepth   => elemR(:,er_FullDepth)
        breadth     => elemSGR(:,esgr_Rectangular_Round_BreadthMax)
        yBot        => elemSGR(:,esgr_Rectangular_Round_Ybot)
        rBot        => elemSGR(:,esgr_Rectangular_Round_Rbot)  
        !%-----------------------------------------------------------------------------
        
        if (depth <= zeroR) then
            outvalue = zeroR

        else if (depth > yBot(indx)) then
            !% top rectangular section
            theta    = twoR * asin(breadth(indx) / twoR / rBot(indx))
            outvalue = rBot(indx) * theta + twoR * (min(depth, fullDepth(indx)) - yBot(indx))

            if (depth >= fullDepth(indx)) then
                !% if the depth exceeds the full depth, the top part will contribute
                !% to overall wetted perimeter
                outvalue = outvalue + breadth(indx)
            end if
        else 
            !% bottom circular section
            theta    = twoR * acos(oneR - depth / rBot(indx))
            outvalue = rBot(indx) * theta
        end if

    end function rect_round_perimeter_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine rect_round_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the hydraulic (average) depth from a known depth in a rectangular_round channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: area(:), hyddepth(:), depth(:),  topwidth(:), fullHydDepth(:)
        !%-----------------------------------------------------------------------------
        thisP       => elemPGx(1:Npack,thisCol) 
        area        => elemR(:,er_Area)
        depth       => elemR(:,er_Depth)
        hyddepth    => elemR(:,er_HydDepth)
        topwidth    => elemR(:,er_Topwidth)
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

    end subroutine rect_round_hyddepth_from_topwidth
!%    
!%==========================================================================  
!%==========================================================================
!%
    real(8) function rect_round_hyddepth_from_topwidth_singular (indx, topwidth, depth) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic depth from known depth for rectangular_round cross section of 
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

    end function rect_round_hyddepth_from_topwidth_singular 
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function rect_round_hydradius_from_depth_singular (indx, depth) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic radius from known depth for a rectangular_round cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer :: breadth(:), fullDepth(:), yBot(:), rBot(:), aBot(:)
        real(8) :: theta, perimeter, area, topDepth
        !%-----------------------------------------------------------------------------
        fullDepth   => elemR(:,er_FullDepth)
        breadth     => elemSGR(:,esgr_Rectangular_Round_BreadthMax)
        aBot        => elemSGR(:,esgr_Rectangular_Round_Abot)
        yBot        => elemSGR(:,esgr_Rectangular_Round_Ybot)
        rBot        => elemSGR(:,esgr_Rectangular_Round_Rbot)  
        !%-----------------------------------------------------------------------------
        
        if (depth <= zeroR) then
            outvalue = zeroR

        else if (depth > yBot(indx)) then
            !% top rectangular section
            theta     = twoR * asin(breadth(indx) / twoR / rBot(indx))
            topDepth  = min(depth, fullDepth(indx)) - yBot(indx)
            perimeter = rBot(indx) * theta + twoR * topDepth

            if (depth >= fullDepth(indx)) then
                !% if the depth exceeds the full depth, the top part will contribute
                !% to overall wetted perimeter
                perimeter = perimeter + breadth(indx)
            end if

            !% find the area
            area     = aBot(indx) + topDepth * breadth(indx)
            outvalue = area / perimeter
        else 
            !% bottom circular section
            theta = twoR * acos(oneR - depth / rBot(indx))
            outvalue = onehalfR * rBot(indx) * (oneR - sin(theta)) / theta
        end if

    end function rect_round_hydradius_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
end module rectangular_round_conduit