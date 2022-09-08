module filled_circular_conduit

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use define_xsect_tables
    use circular_conduit, only: circular_get_normalized_depth_from_area_analytical
    use xsect_tables
    use utility, only: util_CLprint

    implicit none

    !%-----------------------------------------------------------------------------
    !% Description:
    !% filled_circular conduit geometry
    !%

    private

    public :: filled_circular_depth_from_volume
    public :: filled_circular_area_from_depth_singular
    public :: filled_circular_topwidth_from_depth
    public :: filled_circular_topwidth_from_depth_singular
    public :: filled_circular_perimeter_from_depth
    public :: filled_circular_perimeter_from_depth_singular
    public :: filled_circular_perimeter_from_hydradius_singular
    public :: filled_circular_hyddepth_from_topwidth
    public :: filled_circular_hyddepth_from_topwidth_singular
    public :: filled_circular_hydradius_from_depth_singular
    public :: filled_circular_normaldepth_from_sectionfactor_singular



    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine filled_circular_depth_from_volume (elemPGx, Npack, thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Only applies on conduits (or non-surcharged filled_circular conduits)
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !% NOTE: this does NOT limit the depth by surcharge height at this point
        !% This will be done after the head is computed.
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:), volume(:), length(:), AoverAfull(:)
        real(8), pointer :: YoverYfull(:), fullArea(:), fullDepth(:)
        real(8), pointer :: Ybottom(:), Abottom(:), tempYfull(:), tempAfull(:)

        integer, allocatable, target :: thisP_analytical(:), thisP_lookup(:)
        integer, target              :: Npack_analytical, Npack_lookup
        integer :: ii
        !%-----------------------------------------------------------------------------
        thisP      => elemPGx(1:Npack,thisCol)
        depth      => elemR(:,er_Depth)
        volume     => elemR(:,er_Volume)
        length     => elemR(:,er_Length)
        fullArea   => elemR(:,er_FullArea)
        fullDepth  => elemR(:,er_fullDepth)
        Abottom    => elemSGR(:,esgr_Filled_Circular_Abot)
        AoverAfull => elemSGR(:,esgr_Filled_Circular_AoverAfull)
        Ybottom    => elemSGR(:,esgr_Filled_Circular_Ybot)
        YoverYfull => elemSGR(:,esgr_Filled_Circular_YoverYfull)
        tempYfull  => elemR(:,er_Temp01)
        tempAfull  => elemR(:,er_Temp02)
        !%-----------------------------------------------------------------------------
        !% calculate a temporary geometry by considering the whole cicrular corss-section
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
            !% retrive the normalized Y/Yfull from the lookup table
            call xsect_table_lookup &
                (YoverYfull, AoverAfull, YCirc, thisP_lookup) 
        end if

        !% finally get the depth by multiplying the normalized depth with full depth
        !% and substract the bottom depth
        depth(thisP) = YoverYfull(thisP) * tempYfull(thisP) - Ybottom(thisP)
        

    end subroutine filled_circular_depth_from_volume
!%
!%==========================================================================      
!%==========================================================================
!%
    real(8) function filled_circular_area_from_depth_singular (indx, depth) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for filled_circular cross section of a single element
        !% The input indx is the row index in full data 2D array.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer    :: fullArea(:), fullDepth(:), Ybottom(:), Abottom(:)
        real(8), pointer    :: AoverAfull(:), YoverYfull(:)
        real(8) :: tempAfull, tempYfull
        !%-----------------------------------------------------------------------------
        !!if (crashYN) return
        fullArea   => elemR(:,er_FullArea)
        fullDepth  => elemR(:,er_fullDepth)
        Abottom    => elemSGR(:,esgr_Filled_Circular_Abot)
        AoverAfull => elemSGR(:,esgr_Filled_Circular_AoverAfull)
        Ybottom    => elemSGR(:,esgr_Filled_Circular_Ybot)
        YoverYfull => elemSGR(:,esgr_Filled_Circular_YoverYfull)
        !%-----------------------------------------------------------------------------
        !% calculate a temporary geometry by considering the whole cicrular corss-section
        tempAfull = fullArea(indx)  + Abottom(indx)
        tempYfull = fullDepth(indx) + Ybottom(indx)

        !% find Y/Yfull
        YoverYfull(indx) = (depth + Ybottom(indx)) / tempYfull

        !% get A/Afull from the lookup table using Y/Yfull
        AoverAfull(indx) = xsect_table_lookup_singular (YoverYfull(indx), ACirc)

        !% finally get the area by multiplying the normalized area with full area
        outvalue = AoverAfull(indx) * tempAfull - Abottom(indx)

    end function filled_circular_area_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine filled_circular_topwidth_from_depth (elemPGx, Npack, thisCol)
        !%
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a filled_circular conduit
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:), topwidth(:), YoverYfull(:)
        real(8), pointer :: fullDepth(:), Ybottom(:), tempYfull(:) 
        !%-----------------------------------------------------------------------------
        !!if (crashYN) return
        thisP      => elemPGx(1:Npack,thisCol)
        depth      => elemR(:,er_Depth)
        topwidth   => elemR(:,er_Topwidth)
        fullDepth  => elemR(:,er_fullDepth)
        Ybottom    => elemSGR(:,esgr_Filled_Circular_Ybot)
        YoverYfull => elemSGR(:,esgr_Filled_Circular_YoverYfull)
        tempYfull  => elemR(:,er_Temp01)
        !%-----------------------------------------------------------------------------
        !% calculate a temporary geometry by considering the whole cicrular corss-section
        tempYfull(thisP) = fullDepth(thisP) + Ybottom(thisP)

        !% HACK: at this point, YoverYfull probably should be calculated already.
        YoverYfull(thisP) = (depth(thisP) + Ybottom(thisP)) / tempYfull(thisP)

        !% retrive the normalized T/Tmax from the lookup table
        !% T/Tmax value is temporarily saved in the topwidth column
        call xsect_table_lookup &
            (topwidth, YoverYfull, TCirc, thisP) 

        !% finally get the topwidth by multiplying the T/Tmax with full depth
        topwidth(thisP) = max(topwidth(thisP) * tempYfull(thisP), setting%ZeroValue%Topwidth)

    end subroutine filled_circular_topwidth_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function filled_circular_topwidth_from_depth_singular (indx,depth) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth for a filled_circular cross section of a single element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer    ::  YoverYfull(:), fullDepth(:), Ybottom(:)
        real(8) :: tempYfull
        !%-----------------------------------------------------------------------------
        fullDepth  => elemR(:,er_fullDepth)
        Ybottom    => elemSGR(:,esgr_Filled_Circular_Ybot)
        YoverYfull => elemSGR(:,esgr_Filled_Circular_YoverYfull)
        !%-----------------------------------------------------------------------------
        !% calculate a temporary geometry by considering the whole cicrular corss-section
        tempYfull = fullDepth(indx) + Ybottom(indx)

        !% find Y/Yfull
        YoverYfull(indx) = (depth + Ybottom(indx)) / tempYfull

        !% get topwidth by first retriving T/Tmax from the lookup table using Y/Yfull
        !% and then myltiplying it with Tmax (fullDepth for filled_circular cross-section)
        outvalue = tempYfull * xsect_table_lookup_singular (YoverYfull(indx), TCirc) 

        !% if topwidth <= zero, set it to zerovalue
        outvalue = max(outvalue, setting%ZeroValue%Topwidth)

    end function filled_circular_topwidth_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine filled_circular_perimeter_from_depth (elemPGx, Npack, thisCol)
        !%
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a filled_circular conduit
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:), hydRadius(:), YoverYfull(:)
        real(8), pointer :: fullDepth(:), perimeter(:), area(:), fullperimeter(:)
        real(8), pointer :: Ybottom(:), Abottom(:), Pbottom(:), Tbottom(:), tempYfull(:)
        !%-----------------------------------------------------------------------------
        !!if (crashYN) return
        thisP      => elemPGx(1:Npack,thisCol)
        depth      => elemR(:,er_Depth)
        area       => elemR(:,er_Area)
        hydRadius  => elemR(:,er_HydRadius)
        perimeter  => elemR(:,er_Perimeter)
        fullDepth  => elemR(:,er_fullDepth)
        Abottom    => elemSGR(:,esgr_Filled_Circular_Abot)
        Ybottom    => elemSGR(:,esgr_Filled_Circular_Ybot)
        Pbottom    => elemSGR(:,esgr_Filled_Circular_Pbot)
        Tbottom    => elemSGR(:,esgr_Filled_Circular_Tbot)
        YoverYfull => elemSGR(:,esgr_Filled_Circular_YoverYfull)
        fullperimeter => elemR(:,er_FullPerimeter)
        tempYfull  => elemR(:,er_Temp01)
        !%-----------------------------------------------------------------------------
        !% calculate a temporary geometry by considering the whole cicrular corss-section
        tempYfull(thisP) = fullDepth(thisP) + Ybottom(thisP)

        YoverYfull(thisP) = (depth(thisP) + Ybottom(thisP)) / tempYfull(thisP)

        !% retrive the normalized R/Rmax from the lookup table
        !% R/Rmax value is temporarily saved in the hydRadius column
        call xsect_table_lookup &
            (hydRadius, YoverYfull, RCirc,  thisP)  

        hydRadius(thisP) = onefourthR * tempYfull(thisP) * hydRadius(thisP)

        !% get the perimeter by dividing area by hydRadius
        perimeter(thisP) = (area(thisP) + Abottom(thisP)) / hydRadius(thisP)

        !% finally get the perimeter by removing the perimeter of the filled portion
        perimeter(thisP) = min(perimeter(thisP) - Pbottom(thisP) + Tbottom(thisP), fullperimeter(thisP))

    end subroutine filled_circular_perimeter_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function filled_circular_perimeter_from_depth_singular &
        (idx, indepth) result(outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the filled_circular conduit perimeter for the given depth on
        !% the element idx
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: idx
            real(8), intent(in) :: indepth
            real(8), pointer :: fullDepth, fullarea, fullperimeter
            real(8), pointer :: Abottom, Ybottom, Pbottom, Tbottom
            real(8) :: hydRadius, YoverYfull, area
            real(8) :: tempYfull, tempAfull
        !%------------------------------------------------------------------
        !% Aliases
            fullDepth     => elemR(idx,er_fullDepth)
            fullarea      => elemR(idx,er_FullArea)
            fullperimeter => elemR(idx,er_FullPerimeter)
            Abottom       => elemSGR(idx,esgr_Filled_Circular_Abot)
            Pbottom       => elemSGR(idx,esgr_Filled_Circular_Pbot)
            Tbottom       => elemSGR(idx,esgr_Filled_Circular_Tbot)
            Ybottom       => elemSGR(idx,esgr_Filled_Circular_Ybot)
        !%------------------------------------------------------------------
        tempYfull  = fullDepth + Ybottom
        tempAfull  = fullarea  + Abottom

        YoverYfull = (indepth + Ybottom) / tempYfull

        !% 000 retrieve normalized A/Amax for this depth from lookup table
        area = xsect_table_lookup_singular (YoverYfull, ACirc)

        !% --- retrive the normalized R/Rmax for this depth from the lookup table
        hydradius =  xsect_table_lookup_singular (YoverYfull, RCirc)  !% 20220506 brh

        !% --- unnormalize
        hydRadius = onefourthR * tempYfull * hydRadius 
        area      = area * tempAfull

        !% --- get the perimeter by dividing area by hydRadius
        outvalue = min(area / hydRadius - Pbottom + Tbottom, fullperimeter)

    end function filled_circular_perimeter_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function filled_circular_perimeter_from_hydradius_singular (indx,hydradius) result (outvalue)
        !%
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes wetted perimeter from known depth for a filled_circular cross section of
        !% a single element
        !%-----------------------------------------------------------------------------
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: hydradius
        real(8), pointer ::  area(:), fullperimeter(:)
        !%-----------------------------------------------------------------------------
        area          => elemR(:,er_Area)
        fullperimeter => elemR(:,er_FullPerimeter)
        !%-----------------------------------------------------------------------------

        outvalue = min(area(indx) / hydRadius, fullperimeter(indx))

        !% HACK: perimeter correction is needed when the pipe is empty

    end function filled_circular_perimeter_from_hydradius_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine filled_circular_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        !%
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the hydraulic (average) depth from a known depth in a filled_circular conduit
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer    :: thisP(:)
        real(8), pointer    :: area(:), topwidth(:), fullHydDepth(:)
        real(8), pointer    :: depth(:), hyddepth(:)
        !%-----------------------------------------------------------------------------
        !!if (crashYN) return
        thisP        => elemPGx(1:Npack,thisCol)
        area         => elemR(:,er_Area)
        topwidth     => elemR(:,er_Topwidth)
        depth        => elemR(:,er_Depth)
        hyddepth     => elemR(:,er_HydDepth)
        fullHydDepth => elemR(:,er_FullHydDepth)
        !%--------------------------------------------------

        !% calculating hydraulic depth needs conditional since,
        !% topwidth can be zero in filled_circular cross section for both
        !% full and empty condition.

        !% when conduit is empty
        where (depth(thisP) <= setting%ZeroValue%Depth)
            hyddepth(thisP) = setting%ZeroValue%Depth

        !% when conduit is not empty
        elsewhere (depth(thisP) > setting%ZeroValue%Depth)
            !% limiter for when the conduit is full
            hyddepth(thisP) = min(area(thisP) / topwidth(thisP), fullHydDepth(thisP))
        endwhere

    end subroutine filled_circular_hyddepth_from_topwidth
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function filled_circular_hyddepth_from_topwidth_singular (indx,topwidth,depth) result (outvalue)
        !%
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic depth from known depth for filled_circular cross section of
        !% a single element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: topwidth, depth
        real(8), pointer    :: area(:), fullHydDepth(:)
        !%-----------------------------------------------------------------------------
        area         => elemR(:,er_Area)
        fullHydDepth => elemR(:,er_FullHydDepth)
        !%--------------------------------------------------

        !% calculating hydraulic depth needs conditional since,
        !% topwidth can be zero in filled_circular cross section for both
        !% full and empty condition.

        !% when conduit is empty
        if (depth <= setting%ZeroValue%Depth) then
            outvalue = setting%ZeroValue%Depth
        else
            !% limiter for when the conduit is full
            outvalue = min(area(indx) / topwidth, fullHydDepth(indx))
        endif

    end function filled_circular_hyddepth_from_topwidth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function filled_circular_hydradius_from_depth_singular (indx,depth) result (outvalue)
        !%
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic radius from known depth for a filled_circular cross section of
        !% a single element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer    :: YoverYfull(:), fullDepth(:), fullArea(:)
        real(8), pointer    :: Ybottom(:), Abottom(:), Pbottom(:), Tbottom(:)
        real(8) :: tempArea, tempYfull, tempAfull, tempPerimeter, tempRadius
        !%-----------------------------------------------------------------------------
        fullArea   => elemR(:,er_FullArea)
        fullDepth  => elemR(:,er_fullDepth)
        Abottom    => elemSGR(:,esgr_Filled_Circular_Abot)
        Pbottom    => elemSGR(:,esgr_Filled_Circular_Pbot)
        Tbottom    => elemSGR(:,esgr_Filled_Circular_Tbot)
        Ybottom    => elemSGR(:,esgr_Filled_Circular_Ybot)
        YoverYfull => elemSGR(:,esgr_Filled_Circular_YoverYfull)
        !%-----------------------------------------------------------------------------
        tempYfull = fullDepth(indx) + Ybottom(indx)
        tempAfull = fullArea(indx)  + Abottom(indx)
        !% find Y/Yfull
        YoverYfull(indx) = (depth + Ybottom(indx)) / tempYfull

        tempArea = xsect_table_lookup_singular (YoverYfull(indx), ACirc)
        tempArea = tempArea * tempAfull - Abottom(indx)
        !% get hydRadius by first retriving R/Rmax from the lookup table using Y/Yfull
        !% and then myltiplying it with Rmax (fullDepth/4)
        tempRadius = onefourthR * tempYfull * xsect_table_lookup_singular (YoverYfull(indx), RCirc)
        tempPerimeter = tempArea / tempRadius 
        
        outvalue = (tempArea - Abottom(indx)) / (tempPerimeter - Pbottom(indx) + Tbottom(indx))

    end function filled_circular_hydradius_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function filled_circular_normaldepth_from_sectionfactor_singular &
         (SFidx, inSF) result (outvalue)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the depth using the input section factor. This result is
        !% the normal depth of the flow. Note that the element MUST be in 
        !% the set of elements with tables stored in the uniformTableDataR
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: SFidx ! index in the section factor table
            real(8), intent(in) :: inSF
            integer, pointer    :: eIdx  !% element index
            real(8), pointer    :: thisTable(:)
            real(8)             :: normInput
        !%------------------------------------------------------------------
            thisTable => uniformTableDataR(SFidx,:,utd_SF_depth_nonuniform)
            eIdx      => uniformTableI(SFidx,uti_elem_idx)
        !%------------------------------------------------------------------
        !% --- normalize the input
        normInput = inSF / uniformTableR(SFidx,utr_SFmax)
        !% --- lookup the normalized depth
        outvalue = (xsect_table_lookup_singular(normInput,thisTable))
        !% --- unnormalize the depth for the output
        outvalue = outvalue * elemR(eIdx,er_fullDepth)

    end function filled_circular_normaldepth_from_sectionfactor_singular
!%
!%==========================================================================
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    !%----------------------------------------------------------------------
    !% Description:
    !%
    !%----------------------------------------------------------------------

    !%----------------------------------------------------------------------
    !%

    !    !%==========================================================================
    ! ! !%
    ! ! subroutine filled_circular_open_head_from_volume (elemPGx, Npack, thisCol)
    ! !     !%-----------------------------------------------------------------------------
    ! !     !% Description:
    ! !     !% Only applies on open conduits (or non-surcharged filled_circular conduits)
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
    ! !     breadth => elemSGR(:,esgr_Filled_Circular_Breadth)
    ! !     zbottom => elemR(:,er_Zbottom)
    ! !     !%-----------------------------------------------------------------------------

    ! !     head(thisP) = zbottom(thisP) + volume(thisP) / (length(thisP) * breadth(thisP))

    ! ! end subroutine filled_circular_open_head_from_volume
    ! !%
    ! !%==========================================================================
    ! !%    !%==========================================================================
    ! !%
    ! ! subroutine filled_circular_area_from_depth (elemPGx, Npack, thisCol)
    ! !     !%-----------------------------------------------------------------------------
    ! !     !% Description:
    ! !     !% Computes area of a filled_circular open conduit given its depth
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
    ! !     breadth => elemSGR(:,esgr_Filled_Circular_Breadth)
    ! !     !%-----------------------------------------------------------------------------

    ! !     area(thisP) = depth(thisP) * breadth(thisP)

    ! ! end subroutine filled_circular_area_from_depth
    ! ! !%
    ! ! !%==========================================================================
    ! !%==========================================================================
    ! !% END OF MODULE
    ! !%+=========================================================================
end module filled_circular_conduit