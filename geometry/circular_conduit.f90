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
    
    !public :: circular_topwidth_from_depth
    !public :: circular_perimeter_from_depth
    
    ! public :: circular_area_from_depth_singular
    
    ! public :: circular_topwidth_from_depth_singular
    
    ! public :: circular_perimeter_from_depth_singular
    ! public :: circular_perimeter_from_hydradius_singular
    ! public :: circular_hyddepth_from_topwidth
    ! !public :: circular_hyddepth_from_topwidth_singular
    ! public :: circular_hydradius_from_depth_singular
    public :: circular_get_normalized_depth_from_area_analytical



    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine circular_depth_from_volume (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on conduits 
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !%-------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            integer, pointer :: thisP_lookup(:), thisP_analytical(:)
            integer, pointer :: thisPA(:), thisPL(:)
            real(8), pointer :: depth(:), volume(:), length(:), AoverAfull(:)
            real(8), pointer :: YoverYfull(:), fullArea(:), fulldepth(:)
            integer, target  :: Npack_analytical, Npack_lookup
        !%---------------------------------------------------------------------
        !%------------- --------------------------------------------------------
        !% Aliases
            depth       => elemR(:,er_Depth)
            volume      => elemR(:,er_Volume)
            length      => elemR(:,er_Length)
            fullArea    => elemR(:,er_FullArea)
            fulldepth   => elemR(:,er_FullDepth)
            AoverAfull  => elemR(:,er_AoverAfull)
            YoverYfull  => elemR(:,er_YoverYfull)
            
            thisP_analytical => elemI(:,ei_Temp01)
            thisP_lookup     => elemI(:,ei_Temp02)
        !%-----------------------------------------------------------------------------

        !% --- compute the relative volume, which is also the relative area
        AoverAfull(thisP) = volume(thisP) / (length(thisP) * fullArea(thisP))

        !% when AoverAfull <= 4%, SWMM5 uses a special function to get the
        !% normalized depth using the central angle, theta

        !% --- pack the circular elements with AoverAfull <= 4% which will use analytical solution
        !%     from French, 1985 by using the central angle theta.
        Npack_analytical = count(AoverAfull(thisP) <= 0.04)
        if (Npack_analytical > zeroI) then

            thisP_analytical(1:Npack_analytical) = pack(thisP,AoverAfull(thisP) <= 0.04)
            thisPA => thisP_analytical(1:Npack_analytical)

            call circular_get_normalized_depth_from_area_analytical &
                (YoverYfull, AoverAfull, Npack_analytical, thisPA)
        end if

        !% --- pack the rest of the circular elements having AoverAfull > 0.04 which will use
        !%     lookup table for interpolation.
        Npack_lookup = count(AoverAfull(thisP) > 0.04)
        if (Npack_lookup > zeroI) then
            
            thisP_lookup(1:Npack_lookup) = pack(thisP,AoverAfull(thisP) > 0.04)
            thisPL => thisP_lookup(1:Npack_lookup)

            call xsect_table_lookup &
               (YoverYfull, AoverAfull, YCirc, thisPL)

        endif

        !% --- unnormalize the depth 
        depth(thisP) = YoverYfull(thisP) * fulldepth(thisP)

        !% ensure the full depth is not exceeded
        depth(thisP) = min(depth(thisP),fulldepth(thisP))

        !% --- clear the temporary storage
        if (Npack_analytical > zeroI) thisPA = nullvalueI
        if (Npack_lookup     > zeroI) thisPL = nullvalueI

    end subroutine circular_depth_from_volume
!%
!%==========================================================================  
!%==========================================================================
!%
    subroutine circular_get_normalized_depth_from_area_analytical &
        (normalizedDepth, normalizedArea, Npack, thisP)
        !%------------------------------------------------------------------
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
        real(8), pointer :: pi
        !%-----------------------------------------------------------------------------
        pi => setting%Constant%pi

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
!%==========================================================================   
!%==========================================================================
!%
    ! subroutine circular_topwidth_from_depth (elemPGx, Npack, thisCol)
    !     !%
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Computes the topwidth from a known depth in a circular conduit
    !     !%-----------------------------------------------------------------------------
    !     integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
    !     integer, pointer :: thisP(:)
    !     real(8), pointer :: depth(:), topwidth(:), YoverYfull(:), fulldepth(:)
    !     !%-----------------------------------------------------------------------------
    !     !!if (crashYN) return
    !     thisP      => elemPGx(1:Npack,thisCol)
    !     depth      => elemR(:,er_Depth)
    !     topwidth   => elemR(:,er_Topwidth)
    !     fulldepth  => elemR(:,er_FullDepth)
    !     YoverYfull => elemSGR(:,esgr_Circular_YoverYfull)
    !     !%-----------------------------------------------------------------------------

    !     !% HACK: at this point, YoverYfull probably should be calculated already.
    !     YoverYfull(thisP) = depth(thisP) / fulldepth(thisP)

    !     !% retrive the normalized T/Tmax from the lookup table
    !     !% T/Tmax value is temporarily saved in the topwidth column
    !     call xsect_table_lookup &
    !         (topwidth, YoverYfull, TCirc, thisP)  !% 20220506brh
    !         !(topwidth, YoverYfull, TCirc, NTCirc, thisP)

    !     !% finally get the topwidth by multiplying the T/Tmax with full depth
    !     topwidth(thisP) = max (topwidth(thisP) * fulldepth(thisP), setting%ZeroValue%Topwidth)

  !  end subroutine circular_topwidth_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine circular_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     !%
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Computes the perimeter from a known depth in a circular conduit
    !     !%-----------------------------------------------------------------------------
    !     integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
    !     integer, pointer :: thisP(:)
    !     real(8), pointer :: depth(:), hydRadius(:), YoverYfull(:)
    !     real(8), pointer :: fulldepth(:), perimeter(:), area(:), fullperimeter(:)
    !     !%-----------------------------------------------------------------------------
    !     !!if (crashYN) return
    !     thisP      => elemPGx(1:Npack,thisCol)
    !     depth      => elemR(:,er_Depth)
    !     area       => elemR(:,er_Area)
    !     hydRadius  => elemR(:,er_HydRadius)
    !     perimeter  => elemR(:,er_Perimeter)
    !     fulldepth  => elemR(:,er_FullDepth)
    !     YoverYfull => elemSGR(:,esgr_Circular_YoverYfull)
    !     fullperimeter => elemR(:,er_FullPerimeter)
    !     !%-----------------------------------------------------------------------------

    !     !% HACK: at this point, YoverYfull probably should be calculated already.
    !     YoverYfull(thisP) = depth(thisP) / fulldepth(thisP)

    !     !% retrive the normalized R/Rmax from the lookup table
    !     !% R/Rmax value is temporarily saved in the hydRadius column
    !     call xsect_table_lookup &
    !         (hydRadius, YoverYfull, RCirc,  thisP)  !% 20220506 brh
    !         !hydRadius, YoverYfull, RCirc, NRCirc, thisP)

    !     hydRadius(thisP) = onefourthR * fulldepth(thisP) * hydRadius(thisP)

    !     !% finally get the perimeter by dividing area by hydRadius
    !     perimeter(thisP) = min (area(thisP) / hydRadius(thisP), fullperimeter(thisP))

    !     !% HACK: perimeter correction is needed when the pipe is empty.
    ! end subroutine circular_perimeter_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    !     subroutine circular_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
    !         !%
    !         !%-----------------------------------------------------------------------------
    !         !% Description:
    !         !% Computes the hydraulic (average) depth from a known depth in a circular conduit
    !         !%-----------------------------------------------------------------------------
    !         integer, target, intent(in) :: elemPGx(:,:)
    !         integer, intent(in) ::  Npack, thisCol
    !         integer, pointer    :: thisP(:)
    !         real(8), pointer    :: area(:), topwidth(:), fullHydDepth(:)
    !         real(8), pointer    :: depth(:), hyddepth(:)
    !         !%-----------------------------------------------------------------------------
    !         !!if (crashYN) return
    !         thisP        => elemPGx(1:Npack,thisCol)
    !         area         => elemR(:,er_Area)
    !         topwidth     => elemR(:,er_Topwidth)
    !         depth        => elemR(:,er_Depth)
    !         hyddepth     => elemR(:,er_HydDepth)
    !         fullHydDepth => elemR(:,er_FullHydDepth)
    !         !%--------------------------------------------------

    !         !% calculating hydraulic depth needs conditional since,
    !         !% topwidth can be zero in circular cross section for both
    !         !% full and empty condition.

    !         !% when conduit is empty
    !         where (depth(thisP) <= setting%ZeroValue%Depth)
    !             hyddepth(thisP) = setting%ZeroValue%Depth

    !         !% when conduit is not empty
    !         elsewhere (depth(thisP) > setting%ZeroValue%Depth)
    !             !% limiter for when the conduit is full
    !             hyddepth(thisP) = min(area(thisP) / topwidth(thisP), fullHydDepth(thisP))
    !         endwhere

!     end subroutine circular_hyddepth_from_topwidth
! !%
! !%==========================================================================
! !% SINGULAR
! !%==========================================================================
! !%
!     ! real(8) function circular_area_from_depth_singular &
!     !     (indx, depth) result (outvalue)
!     !     !%-----------------------------------------------------------------------------
!     !     !% Description:
!     !     !% Computes area from known depth for circular cross section of a single element
!     !     !% The input indx is the row index in full data 2D array.
!     !     !%-----------------------------------------------------------------------------
!     !     integer, intent(in) :: indx
!     !     real(8), intent(in) :: depth
!     !     real(8), pointer    :: fullArea(:), fulldepth(:)
!     !     real(8) :: AoverAfull, YoverYfull
!     !     !%-----------------------------------------------------------------------------
!     !     !!if (crashYN) return
!     !     fullArea   => elemR(:,er_FullArea)
!     !     fulldepth  => elemR(:,er_FullDepth)
!     !     !%-----------------------------------------------------------------------------

!     !     !% find Y/Yfull
!     !     YoverYfull = depth / fulldepth(indx)

!     !     !% get A/Afull from the lookup table using Y/Yfull
!     !     AoverAfull = xsect_table_lookup_singular (YoverYfull, ACirc) !% 20220506brh removed NACirc

!     !     !% finally get the area by multiplying the normalized area with full area
!     !     outvalue = AoverAfull * fullArea(indx)

!     ! end function circular_area_from_depth_singular
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     ! real(8) function circular_topwidth_from_depth_singular &
!     !     (indx,depth) result (outvalue)
!     !     !%-----------------------------------------------------------------------------
!     !     !% Description:
!     !     !% Computes the topwidth for a circular cross section of a single element
!     !     !%-----------------------------------------------------------------------------
!     !     integer, intent(in) :: indx
!     !     real(8), intent(in) :: depth
!     !     real(8), pointer    :: fulldepth(:)
!     !     real(8) :: YoverYfull
!     !     !%-----------------------------------------------------------------------------
!     !     fulldepth  => elemR(:,er_FullDepth)
!     !     !%-----------------------------------------------------------------------------

!     !     !% find Y/Yfull
!     !     YoverYfull = depth / fulldepth(indx)

!     !     !% get topwidth by first retriving T/Tmax from the lookup table using Y/Yfull
!     !     !% and then myltiplying it with Tmax (fullDepth for circular cross-section)
!     !     outvalue = fulldepth(indx) * xsect_table_lookup_singular (YoverYfull, TCirc) !% 20220506brh removed NTCirc

!     !     !% if topwidth <= zero, set it to zerovalue
!     !     outvalue = max(outvalue, setting%ZeroValue%Topwidth)

!     ! end function circular_topwidth_from_depth_singular
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     real(8) function circular_hydradius_from_depth_singular &
!          (indx,depth) result (outvalue)
!         !%
!         !%-----------------------------------------------------------------------------
!         !% Description:
!         !% Computes hydraulic radius from known depth for a circular cross section of
!         !% a single element
!         !%-----------------------------------------------------------------------------
!         integer, intent(in) :: indx
!         real(8), intent(in) :: depth
!         real(8), pointer    :: fulldepth(:)
!         real(8) :: YoverYfull
!         !%-----------------------------------------------------------------------------
!         fulldepth  => elemR(:,er_FullDepth)
!         !%-----------------------------------------------------------------------------

!         !% find Y/Yfull
!         YoverYfull = depth / fulldepth(indx)

!         !% get hydRadius by first retriving R/Rmax from the lookup table using Y/Yfull
!         !% and then myltiplying it with Rmax (fullDepth/4)
!         outvalue = onefourthR * fulldepth(indx) * &
!                 xsect_table_lookup_singular (YoverYfull, RCirc) !% 20220506brh removed NRCirc

!     end function circular_hydradius_from_depth_singular
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     real(8) function circular_perimeter_from_depth_singular  &
!         (idx, indepth) result(outvalue)
!         !%------------------------------------------------------------------
!         !% Description:
!         !% Computes the circular conduit perimeter for the given depth on
!         !% the element idx
!         !%------------------------------------------------------------------
!         !% Declarations:
!             integer, intent(in) :: idx
!             real(8), intent(in) :: indepth
!             real(8), pointer :: fulldepth, fullarea, fullperimeter
!             real(8) :: hydRadius, YoverYfull, area

!         !%------------------------------------------------------------------
!         !% Aliases
!             fulldepth     => elemR(idx,er_FullDepth)
!             fullarea      => elemR(idx,er_FullArea)
!             fullperimeter => elemR(idx,er_FullPerimeter)
!         !%------------------------------------------------------------------
!         YoverYfull = indepth / fulldepth

!         !% 000 retrieve normalized A/Amax for this depth from lookup table
!         area = xsect_table_lookup_singular (YoverYfull, ACirc)

!         !% --- retrive the normalized R/Rmax for this depth from the lookup table
!         hydradius =  xsect_table_lookup_singular (YoverYfull, RCirc)  !% 20220506 brh

!         !% --- unnormalize
!         hydRadius = hydradius * fullArea / fullPerimeter
!         area      = area * fullArea

!         !% --- get the perimeter by dividing area by hydRadius
!         outvalue = min(area / hydRadius, fullperimeter)

!     end function circular_perimeter_from_depth_singular
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     real(8) function circular_perimeter_from_hydradius_singular &
!          (indx,hydradius) result (outvalue)
!         !%
!         !%-----------------------------------------------------------------------------
!         !% Description:
!         !% Computes wetted perimeter from known depth for a circular cross section of
!         !% a single element
!         !%-----------------------------------------------------------------------------
!         !%-----------------------------------------------------------------------------
!         integer, intent(in) :: indx
!         real(8), intent(in) :: hydradius
!         real(8), pointer ::  area(:), fullperimeter(:)
!         !%-----------------------------------------------------------------------------
!         area          => elemR(:,er_Area)
!         fullperimeter => elemR(:,er_FullPerimeter)
!         !%-----------------------------------------------------------------------------

!         outvalue = min(area(indx) / hydRadius, fullperimeter(indx))

!         !% HACK: perimeter correction is needed when the pipe is empty

!     end function circular_perimeter_from_hydradius_singular
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     ! real(8) function circular_hyddepth_from_depth_singular &
!     !     (indx,depth) result (outvalue)
!     !     !%
!     !     !%-----------------------------------------------------------------------------
!     !     !% Description:
!     !     !% Computes hydraulic depth from known depth for circular cross section of
!     !     !% a single element
!     !     !%-----------------------------------------------------------------------------
!     !     integer, intent(in) :: indx
!     !     real(8), intent(in) :: depth
!     !     real(8), pointer    :: fullHydDepth, fullDepth
!     !     real(8)             :: topwidth, area
!     !     !%-----------------------------------------------------------------------------
!     !     !area         => elemR(:,er_Area)
!     !     fullHydDepth => elemR(indx,er_FullHydDepth)
!     !     fullDepth    => elemR(indx,er_FullDepth)
!     !     !%--------------------------------------------------

!     !     topwidth = circular_topwidth_from_depth_singular (indx,depth)
!     !     area     = circular_area_from_depth_singular (indx, depth)

!     !     if (depth <= setting%ZeroValue%Depth) then
!     !         !% --- empty
!     !         outvalue = setting%ZeroValue%Depth
!     !     elseif (depth >= fullHydDepth)
!     !         !% --- full
!     !         outvalue = fullHydDepth
!     !     else
!     !         !% --- otherwise
!     !         outvalue = area / topwidth
!     !     endif

!     ! end function circular_hyddepth_from_depth_singular
! !%
! !%==========================================================================

! !%==========================================================================
! !%
!     real(8) function circular_normaldepth_from_sectionfactor_singular &
!          (SFidx, inSF) result (outvalue)
!         !%------------------------------------------------------------------
!         !% Description
!         !% Computes the depth using the input section factor. This result is
!         !% the normal depth of the flow. Note that the element MUST be in 
!         !% the set of elements with tables stored in the uniformTableDataR
!         !%------------------------------------------------------------------
!         !% Declarations
!             integer, intent(in) :: SFidx ! index in the section factor table
!             real(8), intent(in) :: inSF
!             integer, pointer    :: eIdx  !% element index
!             real(8), pointer    :: thisTable(:)
!             real(8)             :: normInput
!         !%------------------------------------------------------------------
!             thisTable => uniformTableDataR(SFidx,:,utd_SF_depth_nonuniform)
!             eIdx      => uniformTableI(SFidx,uti_elem_idx)
!         !%------------------------------------------------------------------
!         !% --- normalize the input
!         normInput = inSF / uniformTableR(SFidx,utr_SFmax)
!         !% --- lookup the normalized depth
!         outvalue = (xsect_table_lookup_singular(normInput,thisTable))
!         !% --- unnormalize the depth for the output
!         outvalue = outvalue * elemR(eIdx,er_FullDepth)

!     end function circular_normaldepth_from_sectionfactor_singular
!%
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module circular_conduit