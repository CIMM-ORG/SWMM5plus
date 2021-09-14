module circular_conduit

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use xsect_tables

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Circular conduit geometry
    !%

    private

    public :: circular_depth_from_volume
    public :: circular_area_from_depth_singular
    public :: circular_topwidth_from_depth
    public :: circular_topwidth_from_depth_singular 
    public :: circular_perimeter_from_depth
    public :: circular_perimeter_from_hydradius_singular
    public :: circular_hyddepth_from_topwidth
    public :: circular_hyddepth_from_topwidth_singular
    public :: circular_hydradius_from_depth_singular


    contains
    !%
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    !%
    subroutine circular_depth_from_volume (elemPGx, Npack, thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Only applies on conduits (or non-surcharged circular conduits)
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !% NOTE: this does NOT limit the depth by surcharge height at this point
        !% This will be done after the head is computed.
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:), volume(:), length(:), AoverAfull(:)
        real(8), pointer :: YoverYfull(:), fullArea(:), fulldepth(:)
        !%-----------------------------------------------------------------------------
        thisP      => elemPGx(1:Npack,thisCol) 
        depth      => elemR(:,er_Depth)
        volume     => elemR(:,er_Volume)
        length     => elemR(:,er_Length)
        fullArea   => elemR(:,er_FullArea)
        fulldepth  => elemR(:,er_FullDepth)
        AoverAfull => elemSGR(:,eSGR_Circular_AoverAfull)
        YoverYfull => elemSGR(:,eSGR_Circular_YoverYfull)
        !%-----------------------------------------------------------------------------  

        AoverAfull(thisP) = volume(thisP) / (length(thisP) * fullArea(thisP))

        !% HACK: when AoverAfull < 4%, SWMM5 uses a special function to get the 
        !% normalized depth using the central angle, theta 
        !% (Page 82 in the SWMM5 Hydraulic manual)
        !% Figure out a way to implement this later.

        !% retrive the normalized Y/Yfull from the lookup table
        call xsect_table_lookup &
            (YoverYfull, AoverAfull, YCirc, NYCirc, thisP)

        !% finally get the depth by multiplying the normalized depth with full depth
        depth(thisP) = YoverYfull(thisP) * fulldepth(thisP)

    end subroutine circular_depth_from_volume
    !%  
    !%==========================================================================
    !%==========================================================================
    !%
    real(8) function circular_area_from_depth_singular (indx) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for circular cross section of a single element
        !% The input indx is the row index in full data 2D array.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), pointer    :: depth(:), AoverAfull(:), YoverYfull(:)
        real(8), pointer    :: fullArea(:), fulldepth(:)
        !%-----------------------------------------------------------------------------
        depth      => elemR(:,er_Depth)
        fullArea   => elemR(:,er_FullArea)
        fulldepth  => elemR(:,er_FullDepth)
        AoverAfull => elemSGR(:,eSGR_Circular_AoverAfull)
        YoverYfull => elemSGR(:,eSGR_Circular_YoverYfull)
        !%----------------------------------------------------------------------------- 

        !% find Y/Yfull
        YoverYfull(indx) = depth(indx) / fulldepth(indx)

        !% get A/Afull from the lookup table using Y/Yfull
        AoverAfull(indx) = xsect_table_lookup_singular (YoverYfull(indx), ACirc, NACirc)

        !% finally get the area by multiplying the normalized area with full area 
        outvalue = AoverAfull(indx) * fullArea(indx)

    end function circular_area_from_depth_singular
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine circular_topwidth_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a circular conduit
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:), topwidth(:), YoverYfull(:), fulldepth(:)
        !%-----------------------------------------------------------------------------
        thisP      => elemPGx(1:Npack,thisCol) 
        depth      => elemR(:,er_Depth)
        topwidth   => elemR(:,er_Topwidth)
        fulldepth  => elemR(:,er_FullDepth)
        YoverYfull => elemSGR(:,eSGR_Circular_YoverYfull)
        !%-----------------------------------------------------------------------------  

        !% HACK: at this point, YoverYfull probably should be calculated already.
        YoverYfull(thisP) = depth(thisP) / fulldepth(thisP)

        !% retrive the normalized T/Tmax from the lookup table
        !% T/Tmax value is temporarily saved in the topwidth column
        call xsect_table_lookup &
            (topwidth, YoverYfull, TCirc, NTCirc, thisP)

        !% finally get the topwidth by multiplying the T/Tmax with full depth
        topwidth(thisP) = max (topwidth(thisP) * fulldepth(thisP), setting%ZeroValue%Topwidth)

    end subroutine circular_topwidth_from_depth
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    real(8) function circular_topwidth_from_depth_singular (indx) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth for a circular cross section of a single element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx 
        real(8), pointer    :: depth(:), YoverYfull(:), fulldepth(:)
        !%-----------------------------------------------------------------------------
        depth      => elemR(:,er_Depth)
        fulldepth  => elemR(:,er_FullDepth)
        YoverYfull => elemSGR(:,eSGR_Circular_YoverYfull)
        !%----------------------------------------------------------------------------- 

        !% find Y/Yfull
        YoverYfull(indx) = depth(indx) / fulldepth(indx)

        !% get topwidth by first retriving T/Tmax from the lookup table using Y/Yfull
        !% and then myltiplying it with Tmax (fullDepth for circular cross-section)
        outvalue = fulldepth(indx) * xsect_table_lookup_singular (YoverYfull(indx), TCirc, NTCirc)

        !% if topwidth <= zero, set it to zerovalue
        outvalue = max(outvalue, setting%ZeroValue%Topwidth)

    end function circular_topwidth_from_depth_singular
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine circular_perimeter_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a circular conduit
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:), hydRadius(:), YoverYfull(:) 
        real(8), pointer :: fulldepth(:), perimeter(:), area(:), fullperimeter(:)
        !%-----------------------------------------------------------------------------
        thisP      => elemPGx(1:Npack,thisCol) 
        depth      => elemR(:,er_Depth)
        area       => elemR(:,er_Area)
        hydRadius  => elemR(:,er_HydRadius)
        perimeter  => elemR(:,er_Perimeter)
        fulldepth  => elemR(:,er_FullDepth)
        YoverYfull => elemSGR(:,eSGR_Circular_YoverYfull)
        fullperimeter => elemR(:,er_FullPerimeter)
        !%-----------------------------------------------------------------------------  

        !% HACK: at this point, YoverYfull probably should be calculated already.
        YoverYfull(thisP) = depth(thisP) / fulldepth(thisP)

        !% retrive the normalized R/Rmax from the lookup table
        !% R/Rmax value is temporarily saved in the hydRadius column
        call xsect_table_lookup &
            (hydRadius, YoverYfull, RCirc, NRCirc, thisP)

        hydRadius(thisP) = onefourthR * fulldepth(thisP) * hydRadius(thisP)

        !% finally get the perimeter by dividing area by hydRadius
        perimeter(thisP) = min (area(thisP) / hydRadius(thisP), fullperimeter(thisP))

        !% HACK: perimeter correction is needed when the pipe is empty.
    end subroutine circular_perimeter_from_depth
    ! !%    
    ! !%==========================================================================    
    ! !%==========================================================================
    ! !%
    real(8) function circular_perimeter_from_hydradius_singular (indx) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes wetted perimeter from known depth for a circular cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), pointer :: hydRadius(:), area(:), fullperimeter(:)
        !%-----------------------------------------------------------------------------
        hydRadius     => elemR(:,er_HydRadius)
        area          => elemR(:,er_Area)
        fullperimeter => elemR(:,er_FullPerimeter)
        !%-----------------------------------------------------------------------------
        
        outvalue = min(area(indx) / hydRadius(indx), fullperimeter(indx))
        !% HACK: perimeter correction is needed when the pipe is empt
    end function circular_perimeter_from_hydradius_singular
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine circular_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the hydraulic (average) depth from a known depth in a circular conduit
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer    :: thisP(:)
        real(8), pointer    :: area(:), topwidth(:), fullHydDepth(:)
        real(8), pointer    :: depth(:), hyddepth(:)
        !%-----------------------------------------------------------------------------
        thisP        => elemPGx(1:Npack,thisCol) 
        area         => elemR(:,er_Area)
        topwidth     => elemR(:,er_Topwidth)
        depth        => elemR(:,er_Depth)
        hyddepth     => elemR(:,er_HydDepth)
        fullHydDepth => elemR(:,er_FullHydDepth)
        !%--------------------------------------------------

        !% calculating hydraulic depth needs conditional since,
        !% topwidth can be zero in circular cross section for both
        !% full and empty condition.

        !% when conduit is empty
        where (depth(thisP) <= setting%ZeroValue%Depth)
            hyddepth(thisP) = setting%ZeroValue%Depth

        !% when conduit is not empty
        elsewhere (depth(thisP) > setting%ZeroValue%Depth)
            !% limiter for when the conduit is full
            hyddepth(thisP) = min(area(thisP) / topwidth(thisP), fullHydDepth(thisP))
        endwhere

    end subroutine circular_hyddepth_from_topwidth
    !%    
    !%==========================================================================  
    !%==========================================================================
    !%
    real(8) function circular_hyddepth_from_topwidth_singular (indx) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic depth from known depth for circular cross section of 
        !% a single element
        !%-----------------------------------------------------------------------------   
        integer, intent(in) :: indx   
        real(8), pointer    :: area(:), topwidth(:), fullHydDepth(:), depth(:)
        !%-----------------------------------------------------------------------------
        depth        => elemR(:,er_Depth)
        area         => elemR(:,er_Area)
        topwidth     => elemR(:,er_Topwidth)
        fullHydDepth => elemR(:,er_FullHydDepth)
        !%--------------------------------------------------  

        !% calculating hydraulic depth needs conditional since,
        !% topwidth can be zero in circular cross section for both
        !% full and empty condition.

        !% when conduit is empty
        if (depth(indx) <= setting%ZeroValue%Depth) then
            outvalue = setting%ZeroValue%Depth
        else
            !% limiter for when the conduit is full
            outvalue = min(area(indx) / topwidth(indx), fullHydDepth(indx))
        endif

    end function circular_hyddepth_from_topwidth_singular
    ! !% 
    ! !%==========================================================================

    ! !%==========================================================================
    ! !%
    real(8) function circular_hydradius_from_depth_singular (indx) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic radius from known depth for a circular cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx 
        real(8), pointer    :: depth(:), YoverYfull(:), fulldepth(:)
        !%-----------------------------------------------------------------------------
        depth      => elemR(:,er_Depth)
        fulldepth  => elemR(:,er_FullDepth)
        YoverYfull => elemSGR(:,eSGR_Circular_YoverYfull)
        !%----------------------------------------------------------------------------- 

        !% find Y/Yfull
        YoverYfull(indx) = depth(indx) / fulldepth(indx)

        !% get hydRadius by first retriving R/Rmax from the lookup table using Y/Yfull
        !% and then myltiplying it with Rmax (fullDepth/4)
        outvalue = onefourthR * fulldepth(indx) * &
                xsect_table_lookup_singular (YoverYfull(indx), RCirc, NRCirc)

    end function circular_hydradius_from_depth_singular
    ! !%    
    ! !%==========================================================================

    ! !%
    ! !%    
    ! !%==========================================================================
    ! !%==========================================================================
    ! !%
    !     !%  
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% 
    !     !%-----------------------------------------------------------------------------

    !     !%-----------------------------------------------------------------------------
    !     !%  
    ! !%
    ! !%    
    ! !%==========================================================================
    ! !%==========================================================================
    ! !%
    !     !%  
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% 
    !     !%-----------------------------------------------------------------------------

    !     !%-----------------------------------------------------------------------------
    !     !%  
    ! !%
    ! !%    
    ! !%==========================================================================
    ! !%==========================================================================
    ! !%
    !     !%  
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% 
    !     !%-----------------------------------------------------------------------------

    !     !%-----------------------------------------------------------------------------
    !     !%  
    ! !%
    ! !%    
    ! !%==========================================================================
    ! !% PRIVATE
    ! !%==========================================================================   
    ! !%  
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% 
    !     !%-----------------------------------------------------------------------------

    !     !%-----------------------------------------------------------------------------
    !     !%  

    !    !%==========================================================================   
    ! ! !%
    ! ! subroutine circular_open_head_from_volume (elemPGx, Npack, thisCol)
    ! !     !%-----------------------------------------------------------------------------
    ! !     !% Description:
    ! !     !% Only applies on open conduits (or non-surcharged circular conduits)
    ! !     !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
    ! !     !% Assumes that volume > 0 is enforced in volume computations.
    ! !     !%-----------------------------------------------------------------------------
    ! !     integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
    ! !     integer, pointer :: thisP(:)
    ! !     real(8), pointer :: head(:), volume(:), length(:), breadth(:), zbottom(:)
    ! !     !%-----------------------------------------------------------------------------
    ! !     thisP   => elemPGx(1:Npack,thisCol) 
    ! !     head    => elemR(:,er_Head)
    ! !     volume  => elemR(:,er_Volume)
    ! !     length  => elemR(:,er_Length)
    ! !     breadth => elemSGR(:,eSGr_circular_Breadth)
    ! !     zbottom => elemR(:,er_Zbottom)
    ! !     !%-----------------------------------------------------------------------------   

    ! !     head(thisP) = zbottom(thisP) + volume(thisP) / (length(thisP) * breadth(thisP))
   
    ! ! end subroutine circular_open_head_from_volume
    ! !%  
    ! !%==========================================================================
    ! !%    !%==========================================================================
    ! !%
    ! ! subroutine circular_area_from_depth (elemPGx, Npack, thisCol)
    ! !     !%-----------------------------------------------------------------------------
    ! !     !% Description:
    ! !     !% Computes area of a circular open conduit given its depth
    ! !     !% Note, does NOT consider any closed top!
    ! !     !%-----------------------------------------------------------------------------
    ! !     integer, target, intent(in) :: elemPGx(:,:)
    ! !     integer, intent(in) ::  Npack, thisCol
    ! !     integer, pointer :: thisP(:)
    ! !     real(8), pointer :: area(:), depth(:), breadth(:)
    ! !     !%-----------------------------------------------------------------------------
    ! !     thisP   => elemPGx(1:Npack,thisCol) 
    ! !     area    => elemR(:,er_Area)
    ! !     depth   => elemR(:,er_Depth)
    ! !     breadth => elemSGR(:,eSGR_circular_Breadth)
    ! !     !%-----------------------------------------------------------------------------

    ! !     area(thisP) = depth(thisP) * breadth(thisP)

    ! ! end subroutine circular_area_from_depth
    ! ! !%
    ! ! !%==========================================================================
    ! !%==========================================================================
    ! !% END OF MODULE
    ! !%+=========================================================================
end module circular_conduit