module arch_conduit

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use define_xsect_tables
    use xsect_tables
    use utility, only: util_CLprint

    implicit none

    !%-----------------------------------------------------------------------------
    !% Description:
    !% arch conduit geometry
    !%

    private

    public :: arch_depth_from_volume
    public :: arch_area_from_depth_singular
    public :: arch_topwidth_from_depth
    public :: arch_topwidth_from_depth_singular
    public :: arch_perimeter_from_depth
    public :: arch_perimeter_from_depth_singular
    public :: arch_perimeter_from_hydradius_singular
    public :: arch_hyddepth_from_topwidth
    public :: arch_hyddepth_from_topwidth_singular
    public :: arch_hydradius_from_depth_singular
    public :: arch_normaldepth_from_sectionfactor_singular



    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine arch_depth_from_volume (elemPGx, Npack, thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Only applies on conduits (or non-surcharged arch conduits)
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
        AoverAfull => elemSGR(:,esgr_Arch_AoverAfull)
        YoverYfull => elemSGR(:,esgr_Arch_YoverYfull)
        !%-----------------------------------------------------------------------------
        !% --- compute the relative volume
        AoverAfull(thisP) = volume(thisP) / (length(thisP) * fullArea(thisP))
      
        !% retrive the normalized Y/Yfull from the lookup table
        call xsect_table_lookup &
            (YoverYfull, AoverAfull, YArch, thisP)  

        !% finally get the depth by multiplying the normalized depth with full depth
        depth(thisP) = YoverYfull(thisP) * fulldepth(thisP)

    end subroutine arch_depth_from_volume
!%
!%==========================================================================      
!%==========================================================================
!%
    real(8) function arch_area_from_depth_singular (indx, depth) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for arch cross section of a single element
        !% The input indx is the row index in full data 2D array.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer    :: AoverAfull(:), YoverYfull(:)
        real(8), pointer    :: fullArea(:), fulldepth(:)
        !%-----------------------------------------------------------------------------
        !!if (crashYN) return
        fullArea   => elemR(:,er_FullArea)
        fulldepth  => elemR(:,er_FullDepth)
        AoverAfull => elemSGR(:,esgr_Arch_AoverAfull)
        YoverYfull => elemSGR(:,esgr_Arch_YoverYfull)
        !%-----------------------------------------------------------------------------

        !% find Y/Yfull
        YoverYfull(indx) = depth / fulldepth(indx)

        !% get A/Afull from the lookup table using Y/Yfull
        AoverAfull(indx) = xsect_table_lookup_singular (YoverYfull(indx), AArch)

        !% finally get the area by multiplying the normalized area with full area
        outvalue = AoverAfull(indx) * fullArea(indx)

    end function arch_area_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine arch_topwidth_from_depth (elemPGx, Npack, thisCol)
        !%
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a arch conduit
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:), topwidth(:), YoverYfull(:), fulldepth(:)
        !%-----------------------------------------------------------------------------
        !!if (crashYN) return
        thisP      => elemPGx(1:Npack,thisCol)
        depth      => elemR(:,er_Depth)
        topwidth   => elemR(:,er_Topwidth)
        fulldepth  => elemR(:,er_FullDepth)
        YoverYfull => elemSGR(:,esgr_Arch_YoverYfull)
        !%-----------------------------------------------------------------------------

        !% Calculate normalized depth
        YoverYfull(thisP) = depth(thisP) / fulldepth(thisP)

        !% retrive the normalized T/Tmax from the lookup table
        !% T/Tmax value is temporarily saved in the topwidth column
        call xsect_table_lookup &
            (topwidth, YoverYfull, TArch, thisP) 

        !% finally get the topwidth by multiplying the T/Tmax with full depth
        topwidth(thisP) = max (topwidth(thisP) * fulldepth(thisP), setting%ZeroValue%Topwidth)

    end subroutine arch_topwidth_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function arch_topwidth_from_depth_singular (indx,depth) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth for a arch cross section of a single element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer    ::  YoverYfull(:), fulldepth(:)
        !%-----------------------------------------------------------------------------
        fulldepth  => elemR(:,er_FullDepth)
        YoverYfull => elemSGR(:,esgr_Arch_YoverYfull)
        !%-----------------------------------------------------------------------------

        !% find Y/Yfull
        YoverYfull(indx) = depth / fulldepth(indx)

        !% get topwidth by first retriving T/Tmax from the lookup table using Y/Yfull
        !% and then myltiplying it with Tmax (fullDepth for arch cross-section)
        outvalue = fulldepth(indx) * xsect_table_lookup_singular (YoverYfull(indx), TArch) 

        !% if topwidth <= zero, set it to zerovalue
        outvalue = max(outvalue, setting%ZeroValue%Topwidth)

    end function arch_topwidth_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine arch_perimeter_from_depth (elemPGx, Npack, thisCol)
        !%
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a arch conduit
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:), hydRadius(:), YoverYfull(:)
        real(8), pointer :: fulldepth(:), perimeter(:), area(:), fullperimeter(:)
        !%-----------------------------------------------------------------------------
        !!if (crashYN) return
        thisP      => elemPGx(1:Npack,thisCol)
        depth      => elemR(:,er_Depth)
        area       => elemR(:,er_Area)
        hydRadius  => elemR(:,er_HydRadius)
        perimeter  => elemR(:,er_Perimeter)
        fulldepth  => elemR(:,er_FullDepth)
        YoverYfull => elemSGR(:,esgr_Arch_YoverYfull)
        fullperimeter => elemR(:,er_FullPerimeter)
        !%-----------------------------------------------------------------------------

        !% calculate normalized depth
        YoverYfull(thisP) = depth(thisP) / fulldepth(thisP)

        !% retrive the normalized R/Rmax from the lookup table
        !% R/Rmax value is temporarily saved in the hydRadius column
        call xsect_table_lookup &
            (hydRadius, YoverYfull, RArch, thisP)  

        hydRadius(thisP) = 0.2464 * fulldepth(thisP) * hydRadius(thisP)

        !% finally get the perimeter by dividing area by hydRadius
        perimeter(thisP) = min (area(thisP) / hydRadius(thisP), fullperimeter(thisP))

        !% HACK: perimeter correction is needed when the pipe is empty.
    end subroutine arch_perimeter_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function arch_perimeter_from_depth_singular &
        (idx, indepth) result(outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the arch conduit perimeter for the given depth on
        !% the element idx
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: idx
            real(8), intent(in) :: indepth
            real(8), pointer :: fulldepth, fullarea, fullperimeter
            real(8) :: hydRadius, YoverYfull, area

        !%------------------------------------------------------------------
        !% Aliases
            fulldepth     => elemR(idx,er_FullDepth)
            fullarea      => elemR(idx,er_FullArea)
            fullperimeter => elemR(idx,er_FullPerimeter)
        !%------------------------------------------------------------------
        YoverYfull = indepth / fulldepth

        !% --- retrieve normalized A/Amax for this depth from lookup table
        area = xsect_table_lookup_singular (YoverYfull, AArch)

        !% --- retrive the normalized R/Rmax for this depth from the lookup table
        hydradius =  xsect_table_lookup_singular (YoverYfull, RArch)  !% 20220506 brh

        !% --- unnormalize
        hydRadius = 0.2464 * fulldepth * hydradius
        area      = area * fullArea

        !% --- get the perimeter by dividing area by hydRadius
        outvalue = min(area / hydRadius, fullperimeter)

    end function arch_perimeter_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function arch_perimeter_from_hydradius_singular (indx,hydradius) result (outvalue)
        !%
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes wetted perimeter from known depth for a arch cross section of
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

    end function arch_perimeter_from_hydradius_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine arch_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        !%
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the hydraulic (average) depth from a known depth in a arch conduit
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
        !% topwidth can be zero in arch cross section for both
        !% full and empty condition.

        !% when conduit is empty
        where (depth(thisP) <= setting%ZeroValue%Depth)
            hyddepth(thisP) = setting%ZeroValue%Depth

        !% when conduit is not empty
        elsewhere (depth(thisP) > setting%ZeroValue%Depth)
            !% limiter for when the conduit is full
            hyddepth(thisP) = min(area(thisP) / topwidth(thisP), fullHydDepth(thisP))
        endwhere

    end subroutine arch_hyddepth_from_topwidth
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function arch_hyddepth_from_topwidth_singular (indx,topwidth,depth) result (outvalue)
        !%
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic depth from known depth for arch cross section of
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
        !% topwidth can be zero in arch cross section for both
        !% full and empty condition.

        !% when conduit is empty
        if (depth <= setting%ZeroValue%Depth) then
            outvalue = setting%ZeroValue%Depth
        else
            !% limiter for when the conduit is full
            outvalue = min(area(indx) / topwidth, fullHydDepth(indx))
        endif

    end function arch_hyddepth_from_topwidth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function arch_hydradius_from_depth_singular (indx,depth) result (outvalue)
        !%
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic radius from known depth for a arch cross section of
        !% a single element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer    ::  YoverYfull(:), fulldepth(:)
        !%-----------------------------------------------------------------------------
        fulldepth  => elemR(:,er_FullDepth)
        YoverYfull => elemSGR(:,esgr_Arch_YoverYfull)
        !%-----------------------------------------------------------------------------

        !% find Y/Yfull
        YoverYfull(indx) = depth / fulldepth(indx)

        !% get hydRadius by first retriving R/Rmax from the lookup table using Y/Yfull
        !% and then myltiplying it with Rmax (fullDepth/4)
        outvalue = 0.2464 * fulldepth(indx) * &
                xsect_table_lookup_singular (YoverYfull(indx), RArch)

    end function arch_hydradius_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function arch_normaldepth_from_sectionfactor_singular &
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
        outvalue = outvalue * elemR(eIdx,er_FullDepth)

    end function arch_normaldepth_from_sectionfactor_singular
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
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module arch_conduit