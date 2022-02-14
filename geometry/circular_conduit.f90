module circular_conduit

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use define_xsect_tables
    use xsect_tables

    implicit none

    !%-----------------------------------------------------------------------------
    !% Description:
    !% Circular conduit geometry
    !%

    private

    public :: circular_depth_from_volume
    !public :: circular_depth_from_volume_singular
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

        integer, allocatable, target :: thisP_analytical(:), thisP_lookup(:)
        integer, target              :: Npack_analytical, Npack_lookup
        !%-----------------------------------------------------------------------------
        thisP      => elemPGx(1:Npack,thisCol)
        depth      => elemR(:,er_Depth)
        volume     => elemR(:,er_Volume)
        length     => elemR(:,er_Length)
        fullArea   => elemR(:,er_FullArea)
        fulldepth  => elemR(:,er_FullDepth)
        AoverAfull => elemSGR(:,esgr_Circular_AoverAfull)
        YoverYfull => elemSGR(:,esgr_Circular_YoverYfull)
        !%-----------------------------------------------------------------------------
        if (icrash) return
        AoverAfull(thisP) = volume(thisP) / (length(thisP) * fullArea(thisP))

        !% when AoverAfull <= 4%, SWMM5 uses a special function to get the
        !% normalized depth using the central angle, theta

        !% pack the circular elements with AoverAfull <= 4% which will use analytical solution
        !% from French, 1985 by using the central angle theta.
        Npack_analytical = count(AoverAfull(thisP) <= 0.04)
        thisP_analytical = pack(thisP,AoverAfull(thisP) <= 0.04)

        !% pack the rest of the circular elements having AoverAfull > 0.04 which will use
        !% lookup table for interpolation.
        Npack_lookup = count(AoverAfull(thisP) > 0.04)
        thisP_lookup = pack(thisP,AoverAfull(thisP) > 0.04)

        if (Npack_analytical > zeroI) then
            call circular_get_normalized_depth_from_area_analytical &
                (YoverYfull, AoverAfull, Npack_analytical, thisP_analytical)
        end if 
    
        if (Npack_lookup > zeroI) then        
            !% retrive the normalized Y/Yfull from the lookup table
            call xsect_table_lookup &
                (YoverYfull, AoverAfull, YCirc, NYCirc, thisP_lookup)
        end if

        !% finally get the depth by multiplying the normalized depth with full depth
        depth(thisP) = YoverYfull(thisP) * fulldepth(thisP)

    end subroutine circular_depth_from_volume
!%
!%==========================================================================   
!%==========================================================================
!%
    ! real(8) function circular_depth_from_volume_singular (indx) result (outvalue)
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Only applies on conduits (or non-surcharged circular conduits)
    !     !%-----------------------------------------------------------------------------
    !     integer, target, intent(in) :: indx
    !     real(8), pointer :: depth(:), volume(:), length(:), AoverAfull(:)
    !     real(8), pointer :: YoverYfull(:), fullArea(:), fulldepth(:)
    !     !%-----------------------------------------------------------------------------
    !     depth      => elemR(:,er_Depth)
    !     volume     => elemR(:,er_Volume)
    !     length     => elemR(:,er_Length)
    !     fullArea   => elemR(:,er_FullArea)
    !     fulldepth  => elemR(:,er_FullDepth)
    !     AoverAfull => elemSGR(:,esgr_Circular_AoverAfull)
    !     YoverYfull => elemSGR(:,esgr_Circular_YoverYfull)
    !     !%-----------------------------------------------------------------------------
    !     if (icrash) return
    !     AoverAfull(indx) = volume(indx) / (length(indx) * fullArea(indx))

    !     !% HACK: when AoverAfull < 4%, SWMM5 uses a special function to get the
    !     !% normalized depth using the central angle, theta
    !     !% (Page 82 in the SWMM5 Hydraulic manual)
    !     !% Figure out a way to implement this later.

    !     !% retrive the normalized Y/Yfull from the lookup table
    !     !call xsect_table_lookup &
    !     !    (YoverYfull, AoverAfull, YCirc, NYCirc, thisP)

    !     YoverYfull(indx) = xsect_table_lookup_singular(AoverAfull(indx), YCirc, NYCirc)

    !     !% finally get the depth by multiplying the normalized depth with full depth
    !     outvalue = YoverYfull(indx) * fulldepth(indx)

    ! end function circular_depth_from_volume_singular
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
        if (icrash) return
        depth      => elemR(:,er_Depth)
        fullArea   => elemR(:,er_FullArea)
        fulldepth  => elemR(:,er_FullDepth)
        AoverAfull => elemSGR(:,esgr_Circular_AoverAfull)
        YoverYfull => elemSGR(:,esgr_Circular_YoverYfull)
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
        if (icrash) return
        thisP      => elemPGx(1:Npack,thisCol)
        depth      => elemR(:,er_Depth)
        topwidth   => elemR(:,er_Topwidth)
        fulldepth  => elemR(:,er_FullDepth)
        YoverYfull => elemSGR(:,esgr_Circular_YoverYfull)
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
        if (icrash) return
        depth      => elemR(:,er_Depth)
        fulldepth  => elemR(:,er_FullDepth)
        YoverYfull => elemSGR(:,esgr_Circular_YoverYfull)
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
        if (icrash) return
        thisP      => elemPGx(1:Npack,thisCol)
        depth      => elemR(:,er_Depth)
        area       => elemR(:,er_Area)
        hydRadius  => elemR(:,er_HydRadius)
        perimeter  => elemR(:,er_Perimeter)
        fulldepth  => elemR(:,er_FullDepth)
        YoverYfull => elemSGR(:,esgr_Circular_YoverYfull)
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
!%
!%==========================================================================
!%==========================================================================
!%
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
        if (icrash) return
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
        if (icrash) return
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
        if (icrash) return
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
!%
!%==========================================================================
!%==========================================================================
!%
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
        if (icrash) return
        depth      => elemR(:,er_Depth)
        fulldepth  => elemR(:,er_FullDepth)
        YoverYfull => elemSGR(:,esgr_Circular_YoverYfull)
        !%-----------------------------------------------------------------------------

        !% find Y/Yfull
        YoverYfull(indx) = depth(indx) / fulldepth(indx)

        !% get hydRadius by first retriving R/Rmax from the lookup table using Y/Yfull
        !% and then myltiplying it with Rmax (fullDepth/4)
        outvalue = onefourthR * fulldepth(indx) * &
                xsect_table_lookup_singular (YoverYfull(indx), RCirc, NRCirc)

    end function circular_hydradius_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine circular_get_normalized_depth_from_area_analytical &
        (normalizedDepth, normalizedArea, Npack, thisP)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% find the YoverYfull from AoverAfull only when AoverAfull <= 0.04
        !% This subroutine uses the analytical derivation from French, 1985 to
        !% calculate the central angle theta. THIS SUBROUTINE IS NOT VECTORIZED
        !%
        !% this piece of the code is adapted from SWMM5-c code
        !% for more reference check SWMM Reference Manual Volume II – Hydraulics, pp-81 
        !% and SWMM5 xsect.c code
        !%-----------------------------------------------------------------------------
        real(8), intent(inout)      :: normalizedDepth(:)
        real(8), intent(in)         :: normalizedArea(:)
        integer, target, intent(in) :: Npack, thisP(:)

        integer          :: ii, jj 
        integer, pointer :: eIdx
        real(8)          :: alpha, theta, theta1, dTheta
        !%-----------------------------------------------------------------------------
        if (icrash) return

        do ii = 1,Npack
            eIdx   => thisP(ii)
            alpha  = normalizedArea(eIdx)

            !% this piece of the code is adapted from SWMM5-c code
            if (alpha >= oneR) then
                normalizedDepth(eIdx) = oneR
            else if (alpha <= zeroR) then
                normalizedDepth(eIdx) = zeroR
            else if (alpha <= 1e-05) then 
                normalizedDepth(eIdx) = (((37.6911*alpha)**onethirdR)**twoR)/16.0
            else
                theta1  = 0.031715 - 12.79384 * alpha + 8.28479 * sqrt(alpha)
                theta = theta1
                do jj = 1,40
                    dTheta = - ((2.0 * pi) * alpha - theta1 + sin(theta1)) / (1.0 - cos(theta1))
                    if (dTheta > oneR) dTheta = sign(oneR,dtheta)
                    theta1 = theta1 - dTheta
                    if (abs(dTheta) <= 0.0001) then
                        theta = theta1
                        exit
                    end if    
                end do
                normalizedDepth(eIdx) = (oneR - cos(theta / twoR)) / twoR
            end if
        end do

    end subroutine circular_get_normalized_depth_from_area_analytical
!%
!%
!%==========================================================================
!%==========================================================================
!%
        !%
        !%-----------------------------------------------------------------------------
        !% Description:
        !%
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%
!%
!%
!%==========================================================================
!%==========================================================================
!%
        !%
        !%-----------------------------------------------------------------------------
        !% Description:
        !%
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%
!%
!%
!%==========================================================================
!%==========================================================================
!%
!%
    !%-----------------------------------------------------------------------------
    !% Description:
    !%
    !%-----------------------------------------------------------------------------

    !%-----------------------------------------------------------------------------
    !%
!%
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    !%-----------------------------------------------------------------------------
    !% Description:
    !%
    !%-----------------------------------------------------------------------------

    !%-----------------------------------------------------------------------------
    !%

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
    ! !     breadth => elemSGR(:,esgr_circular_Breadth)
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
    ! !     breadth => elemSGR(:,esgr_circular_Breadth)
    ! !     !%-----------------------------------------------------------------------------

    ! !     area(thisP) = depth(thisP) * breadth(thisP)

    ! ! end subroutine circular_area_from_depth
    ! ! !%
    ! ! !%==========================================================================
    ! !%==========================================================================
    ! !% END OF MODULE
    ! !%+=========================================================================
end module circular_conduit