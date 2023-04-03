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
    public :: filled_circular_topwidth_from_depth
    public :: filled_circular_hydradius_and_perimeter_from_depth  


    
    
!     public :: filled_circular_area_from_depth_singular
    
!     public :: filled_circular_topwidth_from_depth_singular
!     public :: filled_circular_perimeter_from_depth
!     public :: filled_circular_perimeter_from_depth_singular
!     public :: filled_circular_perimeter_from_hydradius_singular
!     public :: filled_circular_hyddepth_from_topwidth
!    ! public :: filled_circular_hyddepth_from_topwidth_singular
!     public :: filled_circular_hydradius_from_depth_singular
!     public :: filled_circular_normaldepth_from_sectionfactor_singular



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
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
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

        !% calculate a temporary geometry by considering the whole cicrular cross-section
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

        !% ensure the full depth is not exceeded
        depth(thisP) = min(depth(thisP),fulldepth(thisP))
        

    end subroutine filled_circular_depth_from_volume
!%
!%========================================================================== 
!%==========================================================================
!%
    subroutine filled_circular_topwidth_from_depth (thisP)
        !%
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

        !% prevent overfull
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
!%==========================================================================
!%==========================================================================
! !%
!     subroutine filled_circular_perimeter_from_depth (elemPGx, Npack, thisCol)
!         !%
!         !%-----------------------------------------------------------------------------
!         !% Description:
!         !% Computes the perimeter from a known depth in a filled_circular conduit
!         !%-----------------------------------------------------------------------------
!         integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
!         integer, pointer :: thisP(:)
!         real(8), pointer :: depth(:), hydRadius(:), YoverYfull(:)
!         real(8), pointer :: fullDepth(:), perimeter(:), area(:), fullperimeter(:)
!         real(8), pointer :: Ybottom(:), Abottom(:), Pbottom(:), Tbottom(:), tempYfull(:)
!         !%-----------------------------------------------------------------------------
!         !!if (crashYN) return
!         thisP      => elemPGx(1:Npack,thisCol)
!         depth      => elemR(:,er_Depth)
!         area       => elemR(:,er_Area)
!         hydRadius  => elemR(:,er_HydRadius)
!         perimeter  => elemR(:,er_Perimeter)
!         fullDepth  => elemR(:,er_fullDepth)
!         Ybottom    => elemR(:,er_SedimentDepth)
!         Abottom    => elemSGR(:,esgr_Filled_Circular_bottomArea)
!         Pbottom    => elemSGR(:,esgr_Filled_Circular_bottomPerimeter)
!         Tbottom    => elemSGR(:,esgr_Filled_Circular_bottomTopwidth)
!         YoverYfull => elemSGR(:,esgr_Filled_Circular_YoverYfull)
!         fullperimeter => elemR(:,er_FullPerimeter)
!         tempYfull  => elemR(:,er_Temp01)
!         !%-----------------------------------------------------------------------------
!         !% calculate a temporary geometry by considering the whole cicrular corss-section
!         tempYfull(thisP) = fullDepth(thisP) + Ybottom(thisP)

!         YoverYfull(thisP) = (depth(thisP) + Ybottom(thisP)) / tempYfull(thisP)

!         !% retrive the normalized R/Rmax from the lookup table
!         !% R/Rmax value is temporarily saved in the hydRadius column
!         call xsect_table_lookup &
!             (hydRadius, YoverYfull, RCirc,  thisP)  

!         hydRadius(thisP) = onefourthR * tempYfull(thisP) * hydRadius(thisP)

!         !% get the perimeter by dividing area by hydRadius
!         perimeter(thisP) = (area(thisP) + Abottom(thisP)) / hydRadius(thisP)

!         !% finally get the perimeter by removing the perimeter of the filled portion
!         perimeter(thisP) = min(perimeter(thisP) - Pbottom(thisP) + Tbottom(thisP), fullperimeter(thisP))

!     end subroutine filled_circular_perimeter_from_depth
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     subroutine filled_circular_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
!         !%
!         !%-----------------------------------------------------------------------------
!         !% Description:
!         !% Computes the hydraulic (average) depth from a known depth in a filled_circular conduit
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
!         !% topwidth can be zero in filled_circular cross section for both
!         !% full and empty condition.

!         !% when conduit is empty
!         where (depth(thisP) <= setting%ZeroValue%Depth)
!             hyddepth(thisP) = setting%ZeroValue%Depth

!         !% when conduit is not empty
!         elsewhere (depth(thisP) > setting%ZeroValue%Depth)
!             !% limiter for when the conduit is full
!             hyddepth(thisP) = min(area(thisP) / topwidth(thisP), fullHydDepth(thisP))
!         endwhere

!     end subroutine filled_circular_hyddepth_from_topwidth
! !%
! !%==========================================================================
! !% SINGULAR

! !%==========================================================================
! !%

! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     real(8) function filled_circular_hydradius_from_depth_singular &
!         (indx,depth) result (outvalue)
!         !%
!         !%-----------------------------------------------------------------------------
!         !% Description:
!         !% Computes hydraulic radius from known depth for a filled_circular cross section of
!         !% a single element
!         !%-----------------------------------------------------------------------------
!         integer, intent(in) :: indx
!         real(8), intent(in) :: depth
!         real(8), pointer    :: YoverYfull(:), fullDepth(:), fullArea(:)
!         real(8), pointer    :: Ybottom(:), Abottom(:), Pbottom(:), Tbottom(:)
!         real(8) :: tempArea, tempYfull, tempAfull, tempPerimeter, tempRadius
!         !%-----------------------------------------------------------------------------
!         fullArea   => elemR(:,er_FullArea)
!         fullDepth  => elemR(:,er_fullDepth)
!         Ybottom    => elemR(:,er_SedimentDepth)
!         Abottom    => elemSGR(:,esgr_Filled_Circular_bottomArea)
!         Pbottom    => elemSGR(:,esgr_Filled_Circular_bottomPerimeter)
!         Tbottom    => elemSGR(:,esgr_Filled_Circular_bottomTopwidth)
!         YoverYfull => elemSGR(:,esgr_Filled_Circular_YoverYfull)
!         !%-----------------------------------------------------------------------------
!         tempYfull = fullDepth(indx) + Ybottom(indx)
!         tempAfull = fullArea(indx)  + Abottom(indx)
!         !% find Y/Yfull
!         YoverYfull(indx) = (depth + Ybottom(indx)) / tempYfull

!         tempArea = xsect_table_lookup_singular (YoverYfull(indx), ACirc)
!         tempArea = tempArea * tempAfull - Abottom(indx)
!         !% get hydRadius by first retriving R/Rmax from the lookup table using Y/Yfull
!         !% and then myltiplying it with Rmax (fullDepth/4)
!         tempRadius = onefourthR * tempYfull * xsect_table_lookup_singular (YoverYfull(indx), RCirc)
!         tempPerimeter = tempArea / tempRadius 
        
!         outvalue = (tempArea - Abottom(indx)) / (tempPerimeter - Pbottom(indx) + Tbottom(indx))

!     end function filled_circular_hydradius_from_depth_singular
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     real(8) function filled_circular_perimeter_from_depth_singular &
!         (idx, indepth) result(outvalue)
!         !%------------------------------------------------------------------
!         !% Description:
!         !% Computes the filled_circular conduit perimeter for the given depth on
!         !% the element idx
!         !%------------------------------------------------------------------
!         !% Declarations:
!             integer, intent(in) :: idx
!             real(8), intent(in) :: indepth
!             real(8), pointer :: fullDepth, fullarea, fullperimeter
!             real(8), pointer :: Abottom, Ybottom, Pbottom, Tbottom
!             real(8) :: hydRadius, YoverYfull, area
!             real(8) :: tempYfull, tempAfull
!         !%------------------------------------------------------------------
!         !% Aliases
!             fullDepth     => elemR(idx,er_fullDepth)
!             fullarea      => elemR(idx,er_FullArea)
!             fullperimeter => elemR(idx,er_FullPerimeter)
!             Ybottom       => elemR(idx,er_SedimentDepth)
!             Abottom       => elemSGR(idx,esgr_Filled_Circular_bottomArea)
!             Pbottom       => elemSGR(idx,esgr_Filled_Circular_bottomPerimeter)
!             Tbottom       => elemSGR(idx,esgr_Filled_Circular_bottomTopwidth)
!         !%------------------------------------------------------------------
!         tempYfull  = fullDepth + Ybottom
!         tempAfull  = fullarea  + Abottom

!         YoverYfull = (indepth + Ybottom) / tempYfull

!         !% 000 retrieve normalized A/Amax for this depth from lookup table
!         area = xsect_table_lookup_singular (YoverYfull, ACirc)

!         !% --- retrive the normalized R/Rmax for this depth from the lookup table
!         hydradius =  xsect_table_lookup_singular (YoverYfull, RCirc)  !% 20220506 brh

!         !% --- unnormalize
!         hydRadius = onefourthR * tempYfull * hydRadius 
!         area      = area * tempAfull

!         !% --- get the perimeter by dividing area by hydRadius
!         outvalue = min(area / hydRadius - Pbottom + Tbottom, fullperimeter)

!     end function filled_circular_perimeter_from_depth_singular
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     real(8) function filled_circular_perimeter_from_hydradius_singular &
!          (indx,hydradius) result (outvalue)
!         !%
!         !%-----------------------------------------------------------------------------
!         !% Description:
!         !% Computes wetted perimeter from known depth for a filled_circular cross section of
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

!     end function filled_circular_perimeter_from_hydradius_singular
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     ! real(8) function filled_circular_hyddepth_from_depth_singular &
!     !     (indx,depth) result (outvalue)
!     !     !%
!     !     !%-----------------------------------------------------------------------------
!     !     !% Description:
!     !     !% Computes hydraulic depth from known depth for filled_circular cross section of
!     !     !% a single element
!     !     !%-----------------------------------------------------------------------------
!     !         integer, intent(in) :: indx
!     !         real(8), intent(in) :: depth
!     !         real(8), pointer    :: fullDepth, fullHydDepth
!     !     !%-----------------------------------------------------------------------------
!     !         fullDepth    => elemR(indx,er_FullDepth)
!     !         fullHydDepth => elemR(indx,er_FullHydDepth)
!     !     !%--------------------------------------------------

!     !     topwidth = filled_circular_topwidth_from_depth_singular (indx,depth)
!     !     area     = filled_circular_area_from_depth_singular (indx, depth)

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

!     ! end function filled_circular_hyddepth_from_topwidth_singular
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     real(8) function filled_circular_normaldepth_from_sectionfactor_singular &
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
!         outvalue = outvalue * elemR(eIdx,er_fullDepth)

!     end function filled_circular_normaldepth_from_sectionfactor_singular
!%
!%==========================================================================
!% END OF MODULE
!%=========================================================================
end module filled_circular_conduit