module geometry_lowlevel

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use define_xsect_tables
    use xsect_tables

    use utility_profiler
    use utility_crash
    use utility, only: util_CLprint, util_syncwrite

    implicit none

!%-----------------------------------------------------------------------------
!% Description:
!% Geometry lowlevel computations
!% NOTES:
!% ..._pure functions (without singular) work on either a single element or 
!%    a packed array of elements in the "thisP" argument
!%    Requires that all values are stored in the elemR or elemSGR arrays
!%
!%..._singular functions work only on a single element and use the depth
!%    provided in the argument, which is NOT necessarily the depth in the
!%    corresponding elemR(indx,er_Depth) storage
!%

    private

    !% --- universal functions (apply to any cross-section)
    public :: llgeo_head_from_depth_pure
    public :: llgeo_area_from_volume_pure
    public :: llgeo_hydradius_from_area_and_perimeter_pure
    public :: llgeo_perimeter_from_hydradius_and_area_pure
    public :: llgeo_hyddepth_from_area_and_topwidth_pure
    public :: llgeo_FullEll_pure
    public :: llgeo_ell_pure
    !public :: llgeo_pressure_head_from_hyddepth_pure

    !% --- table look-ups (subroutines for arrays)
    public :: llgeo_tabular_depth_from_volume
    public :: llgeo_tabular_area_from_depth
    public :: llgeo_tabular_topwidth_from_depth
    public :: llgeo_tabular_hydradius_from_depth
    public :: llgeo_tabular_hydradius_from_sectionfactor_and_area

    !% --- table look-ups (singular)
    public :: llgeo_tabular_from_depth_singular
    public :: llgeo_tabular_hydradius_from_area_and_sectionfactor_singular

    !% --- above-full: extended openchannel geometry when overflow not allowed
    public :: llgeo_openchannel_depth_above_full_pure
    public :: llgeo_openchannel_topwidth_above_full_pure
    public :: llgeo_openchannel_perimeter_above_full_pure

    !% --- depth from volume -- OPEN CHANNEL
    public :: llgeo_parabolic_depth_from_volume_pure
    public :: llgeo_powerfunction_depth_from_volume_pure !% HACK -- NOT DONE
    public :: llgeo_rectangular_depth_from_volume_pure
    public :: llgeo_trapezoidal_depth_from_volume_pure
    public :: llgeo_triangular_depth_from_volume_pure

    !% --- area from depth -- OPEN CHANNEL
    public :: llgeo_parabolic_area_from_depth_pure
    public :: llgeo_powerfunction_area_from_depth_pure  !% HACK -- NOT DONE
    public :: llgeo_rectangular_area_from_depth_pure
    public :: llgeo_trapezoidal_area_from_depth_pure
    public :: llgeo_triangular_area_from_depth_pure

    !% --- area from depth -- CLOSED CONDUITS 
    public :: llgeo_filled_circular_area_from_depth_singular
    public :: llgeo_mod_basket_area_from_depth_singular
    public :: llgeo_rectangular_closed_area_from_depth_singular
    public :: llgeo_rect_round_area_from_depth_singular
    public :: llgeo_rectangular_triangular_area_from_depth_singular
    !public :: llgeo_semi_circular_area_from_depth_singular

    !% --- topwidth from depth -- OPEN CHANNEL
    public :: llgeo_parabolic_topwidth_from_depth_pure
    public :: llgeo_powerfunction_topwidth_from_depth_pure  !% HACK -- NOT DONE
    public :: llgeo_rectangular_topwidth_from_depth_pure
    public :: llgeo_trapezoidal_topwidth_from_depth_pure
    public :: llgeo_triangular_topwidth_from_depth_pure  

    !% --- topwidth from depth -- CLOSED CONDUITS
    public :: llgeo_filled_circular_topwidth_from_depth_singular
    public :: llgeo_mod_basket_topwidth_from_depth_singular
    public :: llgeo_rectangular_closed_topwidth_from_depth_singular
    public :: llgeo_rect_round_topwidth_from_depth_singular
    public :: llgeo_rectangular_triangular_topwidth_from_depth_singular
    
    !% --- perimeter from depth -- OPEN CHANNEL
    public :: llgeo_parabolic_perimeter_from_depth_pure
    public :: llgeo_powerfunction_perimeter_from_depth_pure
    public :: llgeo_rectangular_perimeter_from_depth_pure
    public :: llgeo_trapezoidal_perimeter_from_depth_pure
    public :: llgeo_triangular_perimeter_from_depth_pure

    !% --- perimeter from depth -- CLOSED CONDUITS
    !public :: llgeo_circular_perimeter_from_depth_singular
    public :: llgeo_filled_circular_perimeter_from_depth_singular
    public :: llgeo_mod_basket_perimeter_from_depth_singular
    public :: llgeo_rectangular_closed_perimeter_from_depth_singular
    public :: llgeo_rect_round_perimeter_from_depth_singular
    public :: llgeo_rectangular_triangular_perimeter_from_depth_singular
    !public :: llgeo_semi_circular_perimeter_from_depth_singular

    
    contains
!%    
!%==========================================================================
!% UNIVERSAL FUNCTIONS (do not depend on cross-section type)
!%==========================================================================
!%
    pure function llgeo_head_from_depth_pure &
            (thisP, depth) 
        !%------------------------------------------------------------------
        !% Description
        !% simple computation head from depth and z bottom
        !%------------------------------------------------------------------
            integer, intent(in)               :: thisP(:)
            real(8), intent(in)               :: depth(:)
            real(8), dimension(size(thisP))   :: llgeo_head_from_depth_pure
        !%------------------------------------------------------------------

            llgeo_head_from_depth_pure &
                 = depth                          &
                 + elemR(thisP,er_SedimentDepth)  &
                 + elemR(thisP,er_Zbottom)

    end function llgeo_head_from_depth_pure
!%
!%==========================================================================    
!%==========================================================================
!%
    pure function llgeo_area_from_volume_pure &
            (thisP, volume)
        !%------------------------------------------------------------------
        !% Description
        !% simple computation of cross-section area
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)
            real(8), intent(in) :: volume(:)
            real(8), dimension(1:size(thisP)) ::  llgeo_area_from_volume_pure
        !%------------------------------------------------------------------

        llgeo_area_from_volume_pure = volume / elemR(thisP,er_Length)

    end function llgeo_area_from_volume_pure
!%
!%==========================================================================    
!%==========================================================================
!%
    pure function llgeo_hydradius_from_area_and_perimeter_pure &
            (thisP, area, perimeter) 
        !%------------------------------------------------------------------
        !% Description
        !% simple computation of hydraulic radius
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)
            real(8), intent(in) :: area(:), perimeter(:)
            real(8), dimension(1:size(thisP)) :: llgeo_hydradius_from_area_and_perimeter_pure
        !%------------------------------------------------------------------

        llgeo_hydradius_from_area_and_perimeter_pure  = area / perimeter

    end function llgeo_hydradius_from_area_and_perimeter_pure
!%
!%==========================================================================
!%==========================================================================
!%
    pure function llgeo_perimeter_from_hydradius_and_area_pure &
            (thisP, hydradius, area)    
        !%------------------------------------------------------------------
        !% Description
        !% simple computation of hydraulic radius
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)
            real(8), intent(in) :: hydradius(:), area(:)
            real(8), dimension(1:size(thisP)) :: llgeo_perimeter_from_hydradius_and_area_pure
        !%------------------------------------------------------------------

            llgeo_perimeter_from_hydradius_and_area_pure &
                =  area / hydradius

    end function llgeo_perimeter_from_hydradius_and_area_pure  
!%
!%==========================================================================
!%==========================================================================
!%
    pure function llgeo_hyddepth_from_area_and_topwidth_pure &
            (thisP, area, topwidth)     
        !%------------------------------------------------------------------
        !% Description
        !% simple computation of hydraulic depth
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)
            real(8), intent(in) :: area(:), topwidth(:)
            real(8), dimension(1:size(thisP)) :: llgeo_hyddepth_from_area_and_topwidth_pure
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        llgeo_hyddepth_from_area_and_topwidth_pure &
             =  area / topwidth

    end function llgeo_hyddepth_from_area_and_topwidth_pure
!%
!%==========================================================================
!%==========================================================================
!%
    pure function llgeo_FullEll_pure &
            (thisP)   
        !%------------------------------------------------------------------
        !% Description
        !% computation of full value of ell as alternative hydraulic depth
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)
            real(8), dimension(1:size(thisP)) :: llgeo_FullEll_pure
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        llgeo_FullEll_pure &
                 = (                                                     &
                    elemR(thisP,er_Zcrown) - elemR(thisP,er_ZbreadthMax) &
                   )                                                     &
                   * elemR(thisP,er_BreadthMax)                          &
                   + elemR(thisP,er_AreaBelowBreadthMax) / elemR(thisP,er_BreadthMax) 

    end function llgeo_FullEll_pure
!%
!%==========================================================================
!%==========================================================================
!%
    pure function llgeo_ell_pure (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% computes the value of "ell" the modified hydraulic depth
        !%
        !% For simplicity, this only functions at the element location
        !% (i.e., you cannot call with different areas or topwidths)
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisP(:)
            real(8), dimension(size(thisP)) :: llgeo_ell_pure
            !real(8), pointer :: head, area, topwidth
            !real(8), pointer :: ZbreadthMax, breadthMax, areaBelowBreadthMax
        !%------------------------------------------------------------------    
        !% Aliases
            ! head                => elemR(indx,er_Head)
            ! area                => elemR(indx,er_Area)
            ! topwidth            => elemR(indx,er_Topwidth)
            ! ZbreadthMax         => elemR(indx,er_ZbreadthMax)
            ! breadthMax          => elemR(indx,er_BreadthMax)
            ! areaBelowBreadthMax => elemR(indx,er_AreaBelowBreadthMax)
        !%------------------------------------------------------------------

        where (elemR(thisP,er_Head) <  elemR(thisP,er_ZbreadthMax))
            llgeo_ell_pure = elemR(thisP,er_Area) / elemR(thisP,er_Topwidth)

        elsewhere
            llgeo_ell_pure = &
                       ( (elemR(thisP,er_Head) - elemR(thisP,er_ZbreadthMax)) &
                          * elemR(thisP,er_BreadthMax)                       &
                          + elemR(thisP,er_AreaBelowBreadthMax)             &
                       ) / elemR(thisP,er_BreadthMax)
        endwhere


        ! if (head .le. ZbreadthMax) then
        !     outvalue =  area / topwidth
        ! else
        !     outvalue = ( (head - ZbreadthMax) * breadthMax &
        !                     + areaBelowBreadthMax ) / breadthMax
        ! end if

    end function llgeo_ell_pure
!%
!%==========================================================================
!% TABLE LOOK-UP (PACKED)
!%==========================================================================
!%
    subroutine llgeo_tabular_depth_from_volume &
        (elemPGx, Npack, thisCol, Ytable)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on tabular look-ups (typically conduits) that have
        !% depth as a function of volume.
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !%-------------------------------------------------------------------
            integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
            real(8),         intent(in) :: Ytable(:)
            integer, pointer :: thisP(:)
            real(8), pointer :: depth(:), volume(:), length(:), AoverAfull(:)
            real(8), pointer :: YoverYfull(:), fullArea(:), fulldepth(:)
        !%---------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%--------------------------------------------------------------------
        !% Aliases
            thisP      => elemPGx(1:Npack,thisCol)
            depth      => elemR(:,er_Depth)
            volume     => elemR(:,er_Volume)
            length     => elemR(:,er_Length)
            fullArea   => elemR(:,er_FullArea)
            fulldepth  => elemR(:,er_FullDepth)
            AoverAfull => elemR(:,er_AoverAfull)
            YoverYfull => elemR(:,er_YoverYfull)
        !%-------------------------------------------------------------------

        !% --- compute the relative volume
        AoverAfull(thisP) = volume(thisP) / (length(thisP) * fullArea(thisP))
      
        !% retrive the normalized Y/Yfull from the lookup table
        call xsect_table_lookup &
            (YoverYfull, AoverAfull, Ytable, thisP)  

        !% finally get the depth by multiplying the normalized depth with full depth
        depth(thisP) = YoverYfull(thisP) * fulldepth(thisP)

        !% ensure the full depth is not exceeded
        depth(thisP) = min(depth(thisP),fulldepth(thisP))
        
    end subroutine llgeo_tabular_depth_from_volume   
!%
!%========================================================================== 
!%==========================================================================
!%
    subroutine llgeo_tabular_area_from_depth &
        (elemPGx, Npack, thisCol, Atable)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on tabular look-ups 
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !%-------------------------------------------------------------------
            integer, target, intent(in)   :: elemPGx(:,:)
            integer,  intent(in)   :: Npack, thisCol
            real(8),  intent(in)    :: Atable(:)
            integer, pointer :: thisP(:)
            real(8), pointer :: depth(:), volume(:),   area(:)
            real(8), pointer :: YoverYfull(:), AoverAfull(:)
            real(8), pointer :: fullArea(:), fulldepth(:)
        !%---------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%--------------------------------------------------------------------
        !% Aliases
            thisP      => elemPGx(1:Npack,thisCol)
            area       => elemR(:,er_Area)
            depth      => elemR(:,er_Depth)
            volume     => elemR(:,er_Volume)
            fullArea   => elemR(:,er_FullArea)
            fulldepth  => elemR(:,er_FullDepth)
            AoverAfull => elemR(:,er_AoverAfull)
            YoverYfull => elemR(:,er_YoverYfull)
        !%-------------------------------------------------------------------

        !% --- compute the relative area
        YoverYfull(thisP) =  depth(thisP)  / fulldepth(thisP)
      
        !% retrive the normalized Y/Yfull from the lookup table
        call xsect_table_lookup &
            (AoverAfull, YoverYfull, Atable, thisP)  

        !% finally get the area by multiplying the normalized area with full area
        area(thisP)  = AoverAfull(thisP) * fullarea(thisP)

        !% ensure the full area is not exceeded
        area(thisP) = min(area(thisP),fullarea(thisP))
        
    end subroutine llgeo_tabular_area_from_depth 
!%
!%==========================================================================     
!%========================================================================== 
!%
    subroutine llgeo_tabular_topwidth_from_depth &
        (elemPGx, Npack, thisCol, Ytable )
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on tabular look-ups (typically conduits) that have
        !% topwidth as a function of depth
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !%-------------------------------------------------------------------
            integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
            real(8),         intent(in) :: Ytable(:)
            integer, pointer :: thisP(:)
            real(8), pointer :: depth(:), topwidth(:), breadthMax(:)
            real(8), pointer :: YoverYfull(:), fulldepth(:)
        !%---------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%--------------------------------------------------------------------
        !% Aliases
            thisP      => elemPGx(1:Npack,thisCol)
            depth      => elemR(:,er_Depth)
            fulldepth  => elemR(:,er_FullDepth)
            topwidth   => elemR(:,er_TopWidth)
            breadthMax => elemR(:,er_BreadthMax)
            YoverYfull => elemR(:,er_YoverYfull)
        !%-------------------------------------------------------------------

        !% --- Calculate normalized depth
        YoverYfull(thisP) = depth(thisP) / fulldepth(thisP)

        !% --- ensure the topwidth fraction is max normalized depth
        YoverYfull(thisP) = min(YoverYfull(thisP),setting%Discretization%FullConduitTopwidthDepthFraction)

        !% --- retrive the normalized T/Tmax from the lookup table
        !%     T/Tmax value is temporarily saved in the topwidth column
        call xsect_table_lookup &
            (topwidth, YoverYfull, YTable, thisP) 

        !% --- get the topwidth by multiplying the T/Tmax with maximum breadth
        topwidth(thisP) = max(topwidth(thisP) * breadthMax(thisP), setting%ZeroValue%Topwidth)
    
    end subroutine llgeo_tabular_topwidth_from_depth
!%  
!%==========================================================================
!%==========================================================================      
!%
    subroutine llgeo_tabular_hydradius_from_depth &
        (elemPGx, Npack, thisCol, Ytable )
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on tabular look-ups (typically conduits) that have
        !% hydraulic radius as a function of depth
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !%-------------------------------------------------------------------
            integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
            real(8),         intent(in) :: Ytable(:)
            integer, pointer :: thisP(:)
            real(8), pointer :: depth(:), HydRadius(:), fullHydRadius(:)
            real(8), pointer :: YoverYfull(:), fulldepth(:)
        !%---------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%--------------------------------------------------------------------
        !% Aliases
            thisP         => elemPGx(1:Npack,thisCol)
            depth         => elemR(:,er_Depth)
            fulldepth     => elemR(:,er_FullDepth)
            HydRadius     => elemR(:,er_HydRadius)
            fullHydRadius => elemR(:,er_FullHydRadius)
            YoverYfull    => elemR(:,er_YoverYfull)
        !%-------------------------------------------------------------------

        !% --- Calculate normalized depth
        YoverYfull(thisP) = depth(thisP) / fulldepth(thisP)

        !% --- ensure 1.0 is max normalized depth
        YoverYfull(thisP) = min(YoverYfull(thisP),oneR)

        !% --- retrive the normalized R/Rmax from the lookup table
        !%     R/Rfull value is temporarily saved in the HydRadius column
        call xsect_table_lookup &
            (HydRadius, YoverYfull, YTable, thisP) 

        !% --- get the hydraulic radius by multiplying the R/Rfull with full hydraulic radius
        HydRadius(thisP) = max (HydRadius(thisP) * fullHydRadius(thisP), setting%ZeroValue%Depth)
    
    end subroutine llgeo_tabular_hydradius_from_depth
!%
!%==========================================================================
!%==========================================================================      
!%
    subroutine llgeo_tabular_hydradius_from_sectionfactor_and_area &
        (elemPGx, Npack, thisCol, Stable, SoverScol )
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on tabular look-ups (typically conduits) that have
        !% hydraulic radius as a function of depth
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !%-------------------------------------------------------------------
            integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
            integer,         intent(in) :: SoverScol
            real(8),         intent(in) :: Stable(:)
            integer, pointer :: thisP(:)
            real(8), pointer :: area(:), fullarea(:), fulldepth(:), fullhydradius(:)
            real(8), pointer :: AoverAfull(:), SoverSfull(:), HydRadius(:)
            real(8), pointer :: SectionFactor(:)
        !%---------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%--------------------------------------------------------------------
        !% Aliases
            thisP         => elemPGx(1:Npack,thisCol)
            area          => elemR(:,er_area)
            fullarea      => elemR(:,er_FullArea)
            fulldepth     => elemR(:,er_Fulldepth)
            fullhydradius => elemR(:,er_FullHydRadius)
            HydRadius     => elemR(:,er_HydRadius)
            AoverAfull    => elemR(:,er_AoverAfull)
            SoverSfull    => elemSGR(:,SoverScol)
            SectionFactor => elemR(:,er_Temp01)
        !%-------------------------------------------------------------------

        !% --- Calculate normalized area
        AoverAfull(thisP) = area(thisP) / fullarea(thisP)

        !% --- prevent overfull
        AoverAfull(thisP) = min(AoverAfull(thisP),oneR)

        !% --- prevent underfull
        AoverAfull(thisP) = max(AoverAfull(thisP), setting%ZeroValue%Area/fullarea(thisP))

        !% --- retrieve the normalized section factor
        call xsect_table_lookup (SoverSfull, AoverAfull, STable, thisP) 

        !% find the section factor 
        !% sF_full = Afull * (rFull) ^ (2/3)
        !% sF = sF_full * sF/sF_full
        SectionFactor(thisP) = &
                ( fullArea(thisP) & 
                  * ((fullhydradius(thisP))**twoThirdR) &
                ) * SoverSfull(thisP)

        !% --- compute hydraulic radius from section factor
        HydRadius(thisP) = (SectionFactor(thisP) / area(thisP)) ** threehalfR

        !% --- clear temporary storage
        SectionFactor(thisP) = nullvalueR
       
    end subroutine llgeo_tabular_hydradius_from_sectionfactor_and_area
!%
!%==========================================================================   
!% TABLE LOOK UP (SINGULAR)
!%========================================================================== 
!%
    real(8) function llgeo_tabular_from_depth_singular &
            (indx, depth, maxValue, zeroValue, Xtable) result (outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on tabular look-ups (typically conduits) as
        !% a function of input depth, which may not be stored in elemR(:,er_Depth)
        !% The zerovalue for depth is an argument so that it can be set
        !% for sediment-filled conduits to ensure the lookup depth is at or above
        !% the sediment layer.
        !% maxValue is the normalization scale for the Xtable type
        !% XTable is (e.g.) ACirc, TCirc, RHorseShoe, etc.
        !%-------------------------------------------------------------------
            integer, intent(in) :: indx
            real(8), intent(in) :: Xtable(:)
            real(8), intent(in) :: maxValue, zeroValue, depth
            real(8), pointer    :: fulldepth, breadthMax
            real(8)             :: YoverYfull
        !%--------------------------------------------------------------------
        !% Aliases
            fulldepth  => elemR(indx,er_FullDepth)
            breadthMax => elemR(indx,er_BreadthMax)
        !%-------------------------------------------------------------------

        !% --- Calculate normalized depth
        YoverYfull = depth / fulldepth

        !% --- prevent overfull
        YoverYfull = min(YoverYfull, oneR)

        !% --- prevent underfull
        YoverYfull = max(YoverYfull, zeroValue/fulldepth)

        !% --- look up
        outvalue = maxValue * xsect_table_lookup_singular(YoverYfull, XTable)
    
    end function llgeo_tabular_from_depth_singular
!%  
!%==========================================================================
!%==========================================================================  
!%   
    real(8) function llgeo_tabular_hydradius_from_area_and_sectionfactor_singular &
            (indx, area, fullHydRadius, zeroValue, Stable) result(outvalue)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the hydraulic radius using the section factor lookup
        !% table. This is the singular funciton version of the subroutine
        !% llgeo_tabular_hydradius_from_sectionfactor_and_area()    
        !% The argument zerovalue is for later use with sediment-filled elements
        !% to ensure the depth is at or above the sediment layer.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: indx
            real(8), intent(in) :: Stable(:),  area, fullHydRadius, zeroValue
            real(8), pointer    :: fullarea
            real(8)             :: AoverAfull, SoverSfull, SFfull
        !%------------------------------------------------------------------  
        !% Aliases
            fullarea => elemR(indx,er_FullArea)
        !%------------------------------------------------------------------ 
        !% --- Calculate normalized area
            AoverAfull = area / fullarea

            !% --- prevent overfull
            AoverAfull = min(AoverAfull,oneR)

            !% --- prevent underfull
            AoverAfull = max(AoverAfull, zeroValue/fullarea)
    
            !% --- retrieve the normalized section factor
            SoverSfull =  xsect_table_lookup_singular(AoverAfull, STable) 
    
            !%--- approach from lgeo_tabular_hydradius_from_sectionfactor_and_area()  
            !     that has been replaced by using the full hydraulic radius
            !% unnormalize section factor -- store in outvalue 
            !% sF_full = Afull * (rFull) ^ (2/3)
            !% sF = sF_full * sF/sF_full
            ! outvalue(thisP) = &
            !         ( fullArea & 
            !           * ((sfparam * fullDepth)**twoThirdR) &
            !         ) * SoverSfull

            !% --- full section factor
            SFfull = fullArea * (fullHydRadius**twoThirdR)
    
            !% --- R = (SF / Area)^(3/4) --- from unnormalized SF
            outvalue = (SFfull * SoverSfull / area) ** threehalfR 

    end function llgeo_tabular_hydradius_from_area_and_sectionfactor_singular
!%        
!%==========================================================================         
!% ABOVE FULL FUNCTIONS -- OPEN CHANNEL
!%==========================================================================
!%
    pure function llgeo_openchannel_depth_above_full_pure &
            (thisP)
        !%------------------------------------------------------------------
        !% Description
        !% computes the depth for any open channel when overflow is not 
        !% allowed. Channel banks are assumed to be vertical upwards
        !% from the maximum breadth, so the extra depth is based on
        !% the rectangular area
        !% NOTE: it is assumed that the set of points (thisP) are only
        !% those points where the volume > full volume and
        !% setting%Discretization%AllowChannelOverflowTF = .true.
        !% is set
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)
            real(8), dimension(size(thisP)) :: llgeo_openchannel_depth_above_full_pure
        !%------------------------------------------------------------------

        llgeo_openchannel_depth_above_full_pure                       &
            = elemR(thisP,er_FullDepth)                               &
            + ( elemR(thisP,er_Volume) - elemR(thisP,er_FullVolume) ) &
            / ( elemR(thisP,er_Length) * elemR(thisP,er_BreadthMax) )

    end function llgeo_openchannel_depth_above_full_pure
!%
!%==========================================================================
!%==========================================================================
!%
    pure function llgeo_openchannel_topwidth_above_full_pure &
            (thisP) 
        !%------------------------------------------------------------------
        !% Description
        !% computes the topwidth for any open channel when overflow is not 
        !% allowed. Channel banks are assumed to be vertical upwards
        !% from the maximum breadth, so the extra depth is based on
        !% the rectangular area
        !% NOTE: it is assumed that the set of points (thisP) are only
        !% those points where the volume > full volume and
        !% setting%Discretization%AllowChannelOverflowTF = .true.
        !% is set
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)
            real(8), dimension(size(thisP)) :: llgeo_openchannel_topwidth_above_full_pure
        !%------------------------------------------------------------------

        llgeo_openchannel_topwidth_above_full_pure = elemR(thisP,er_BreadthMax)


    end function llgeo_openchannel_topwidth_above_full_pure
!%
!%==========================================================================
!%==========================================================================
!%
    pure function llgeo_openchannel_perimeter_above_full_pure &
            (thisP) 
        !%------------------------------------------------------------------
        !% Description
        !% computes the perimeter for any open channel when overflow is not 
        !% allowed. Channel banks are assumed to be vertical upwards
        !% from the maximum breadth, so the extra depth is based on
        !% the rectangular area
        !% NOTE: it is assumed that the set of points (thisP) are only
        !% those points where the volume > full volume and
        !% setting%Discretization%AllowChannelOverflowTF = .true.
        !% is set
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)
            real(8), dimension(size(thisP)) :: llgeo_openchannel_perimeter_above_full_pure
        !%------------------------------------------------------------------

        !% --- new perimeter is full plus twice the depth of the extra volume
        llgeo_openchannel_perimeter_above_full_pure &
                = elemR(thisP,er_FullPerimeter) &
                + twoR * (elemR(thisP,er_Volume) - elemR(thisP,er_FullVolume)) &
                / (elemR(thisP,er_Length) * elemR(thisP,er_Breadthmax))

    end function llgeo_openchannel_perimeter_above_full_pure
!%
!%==========================================================================
!% DEPTH FROM VOLUME -- OPEN CHANNEL
!%==========================================================================
!%
    pure function llgeo_parabolic_depth_from_volume_pure &
            (thisP, volume)
        !%------------------------------------------------------------------
        !% Description
        !% computes parabolic open channel depth on input vector
        !% of points
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)  !% may be packed array of indexes
            real(8), intent(in) :: volume(:)
            real(8), dimension(size(thisP)) :: llgeo_parabolic_depth_from_volume_pure
        !%------------------------------------------------------------------

            llgeo_parabolic_depth_from_volume_pure                         &
                  = ( threefourthR                                         &
                    * ( volume / elemR(thisP,er_Length)    &
                        )                                                  &
                        / elemSGR(thisP,esgr_Parabolic_Radius)             &
                    ) ** twothirdR

    end function llgeo_parabolic_depth_from_volume_pure
!%
!%==========================================================================
!%==========================================================================
!%
    pure function llgeo_powerfunction_depth_from_volume_pure &
            (thisP, volume) 
        !%------------------------------------------------------------------
        !% Description
        !% computes power function open channel depth on input vector
        !% of points
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:) !% may be packed array of indexes
            real(8), intent(in) :: volume(:)
            real(8), dimension(size(thisP)) :: llgeo_powerfunction_depth_from_volume_pure
        !%------------------------------------------------------------------

        llgeo_powerfunction_depth_from_volume_pure = nullvalueR  !% STUB  HACK

    end function llgeo_powerfunction_depth_from_volume_pure
!%
!%==========================================================================
!%==========================================================================
!%
    pure function llgeo_rectangular_depth_from_volume_pure &
            (thisP, volume) 
        !%------------------------------------------------------------------
        !% Description
        !% computes rectangular open channel from depth on input vector
        !% of points
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)  !% may be packed array of indexes
            real(8), intent(in) :: volume(:)
            real(8), dimension(size(thisP)) :: llgeo_rectangular_depth_from_volume_pure
        !%------------------------------------------------------------------

        llgeo_rectangular_depth_from_volume_pure &
            = volume &
            / ( elemR(thisP,er_Length) * elemSGR(thisP,esgr_Rectangular_Breadth) )

    end function llgeo_rectangular_depth_from_volume_pure
!%
!%==========================================================================
!%==========================================================================
!%
    pure function llgeo_trapezoidal_depth_from_volume_pure &
            (thisP, volume)
        !%------------------------------------------------------------------
        !% Description
        !% computes depth from volume on packed elements thisP
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)
            real(8), intent(in) :: volume(:)
            real(8), dimension(size(thisP)) :: llgeo_trapezoidal_depth_from_volume_pure
        !%------------------------------------------------------------------

            llgeo_trapezoidal_depth_from_volume_pure                                         &
                = - onehalfR *                                                               &
                ( elemSGR(thisP,esgr_Trapezoidal_Breadth)                                    &
                      / (onehalfR * (  elemSGR(thisP,esgr_Trapezoidal_LeftSlope)             &
                                     + elemSGR(thisP,esgr_Trapezoidal_RightSlope)))          &
                       - sqrt( (elemSGR(thisP,esgr_Trapezoidal_Breadth)                      &
                               / (onehalfR * (  elemSGR(thisP,esgr_Trapezoidal_LeftSlope)    &
                                             + elemSGR(thisP,esgr_Trapezoidal_RightSlope)))  &
                               ) ** twoR                                                     &
                               + fourR * volume                              &
                                  / (onehalfR * elemR(thisP,er_Length)                       &
                                     * (  elemSGR(thisP,esgr_Trapezoidal_LeftSlope)          &
                                        + elemSGR(thisP,esgr_Trapezoidal_RightSlope))        &
                                    )                                                        &
                             )                                                               &
                )

    end function llgeo_trapezoidal_depth_from_volume_pure   
!%  
!%==========================================================================
!%==========================================================================
!%     
    pure function llgeo_triangular_depth_from_volume_pure  &
            (thisP, volume) 
        !%------------------------------------------------------------------
        !% Description
        !% computes depth from volume on packed elements thisP
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)
            real(8), intent(in) :: volume(:)
            real(8), dimension(size(thisP)) :: llgeo_triangular_depth_from_volume_pure
        !%------------------------------------------------------------------

        llgeo_triangular_depth_from_volume_pure                      &
                 = sqrt( volume                                      &
                         / (  elemR(thisP,er_Length)                 &
                            * elemSGR(thisP,esgr_Triangular_Slope)   &
                           )                                         &
                       )

    end function llgeo_triangular_depth_from_volume_pure 
!%  
!%==========================================================================
!% AREA FROM DEPTH -- OPEN CHANNEL   
!%==========================================================================
!%
    pure function llgeo_parabolic_area_from_depth_pure &
            (thisP, depth)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for parabolic cross section
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)  ! may be a packed array of indexes
            real(8), intent(in) :: depth(:)
            real(8), dimension(size(thisP)) :: llgeo_parabolic_area_from_depth_pure
        !%------------------------------------------------------------------
        
        llgeo_parabolic_area_from_depth_pure                               &
                = (fourR / threeR) * elemSGR(thisP, esgr_Parabolic_Radius) &
                   * depth *  sqrt(depth)
        
    end function llgeo_parabolic_area_from_depth_pure 
!%
!%==========================================================================
!%==========================================================================
!% 
    pure function llgeo_powerfunction_area_from_depth_pure &
            (thisP, depth) 
        !%------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for parabolic cross section
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)  ! may be a packed array of indexes
            real(8), intent(in) :: depth(:)
            real(8), dimension(size(thisP)) :: llgeo_powerfunction_area_from_depth_pure
        !%------------------------------------------------------------------

        llgeo_powerfunction_area_from_depth_pure = nullvalueR  !% STUB  HACK

    end function llgeo_powerfunction_area_from_depth_pure
!%
!%==========================================================================
!%==========================================================================
!%
    pure function llgeo_rectangular_area_from_depth_pure &
            (thisP, depth) 
        !%------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for rectangular cross section
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)  ! may be a packed array of indexes
            real(8), intent(in) :: depth(:)
            real(8), dimension(size(thisP)) :: llgeo_rectangular_area_from_depth_pure
        !%------------------------------------------------------------------

        llgeo_rectangular_area_from_depth_pure                        &
             = depth * elemR(thisP,er_BreadthMax)

    end function llgeo_rectangular_area_from_depth_pure
!%
!%==========================================================================
!%==========================================================================
!%
    pure function llgeo_trapezoidal_area_from_depth_pure  &
            (thisP, depth) 
        !%------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for trapezoidal cross section 
        !% input may be a scalar or a packed array of indexes
        !%-------------------------------------------------------------------
            integer, intent(in) :: thisP(:)
            real(8), intent(in) :: depth(:)
            real(8), dimension(size(thisP)) ::  llgeo_trapezoidal_area_from_depth_pure
        !%-------------------------------------------------------------------

        llgeo_trapezoidal_area_from_depth_pure                                     &
                 = (elemSGR(thisP,esgr_Trapezoidal_Breadth)                        &
                        + onehalfR * (  elemSGR(thisP,esgr_Trapezoidal_LeftSlope)  &
                                      + elemSGR(thisP,esgr_Trapezoidal_RightSlope) &
                                    )  * depth                    &
                    ) * depth

    end function llgeo_trapezoidal_area_from_depth_pure
!%
!%==========================================================================
!%==========================================================================
!%
    pure function llgeo_triangular_area_from_depth_pure &
            (thisP, depth) 
        !%----------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for triangular cross section
        !%----------------------------------------------------------------------
            integer, intent(in) :: thisP(:)  ! may be a packed array of indexes
            real(8), intent(in) :: depth(:)
            real(8), dimension(size(thisP)) :: llgeo_triangular_area_from_depth_pure
        !%----------------------------------------------------------------------

        llgeo_triangular_area_from_depth_pure           &
                 = elemSGR(thisP,esgr_Triangular_Slope) &
                    * (depth ** twoR) 

    end function llgeo_triangular_area_from_depth_pure
!%
!%==========================================================================
!%  AREA FROM DEPTH -- CLOSED CONDUITS
!%==========================================================================
!%
    real(8) function llgeo_filled_circular_area_from_depth_singular &
        (indx, depth) result (outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for filled_circular cross section 
        !% of a single element;
        !% The input indx is the row index in full data 2D array.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: indx
            real(8), intent(in) :: depth
            real(8), pointer    :: fullArea, fullDepth
            real(8), pointer    :: sedimentDepth, sedimentArea
            real(8)             :: AoverAfull, YoverYfull
            real(8)             :: totalFullArea, totalFullDepth
        !%-----------------------------------------------------------------
        !% Aliases:
            fullArea      => elemR(indx,er_FullArea)
            fullDepth     => elemR(indx,er_fullDepth)
            sedimentDepth => elemR(indx,er_SedimentDepth)
            sedimentArea  => elemSGR(indx,esgr_Filled_Circular_bottomArea)
        !%------------------------------------------------------------------

        !% calculate a temporary geometry by considering the whole cicrular cross-section
        totalFullArea  = fullArea  + sedimentArea
        totalFullDepth = fullDepth + sedimentDepth

        !% find Y/Yfull
        YoverYfull = (depth + sedimentDepth) / totalFullDepth

        !% get A/Afull from the lookup table using Y/Yfull
        AoverAfull = xsect_table_lookup_singular (YoverYfull, ACirc)

        !% --- get the area by multiplying the normalized area with full area
        !%     and subtracting the sediment area
        outvalue = AoverAfull * totalFullArea - sedimentArea

    end function llgeo_filled_circular_area_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function llgeo_mod_basket_area_from_depth_singular &
            (indx, depth) result (outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for mod_basket cross section of a single element
        !% The input indx is the row index in full data 2D array.
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: indx
            real(8), intent(in) :: depth
            real(8), pointer :: fulldepth, fullArea
            real(8), pointer :: depthAtBreadthMax, breadthMax, rTop
            real(8) :: emptyDepth, emptyTheta, emptyArea
        !%-----------------------------------------------------------------
        !% Aliases:
            fullArea          => elemR(indx,er_FullArea)
            fulldepth         => elemR(indx,er_FullDepth)
            depthAtBreadthMax => elemR(indx,er_DepthAtBreadthMax) 
            breadthMax        => elemR(indx,er_BreadthMax)
            rTop              => elemSGR(indx,esgr_Mod_Basket_Rtop)
        !%-----------------------------------------------------------------
        
        if (depth <= setting%ZeroValue%Depth) then
            outvalue = setting%ZeroValue%Area
        
        elseif ((depth > setting%ZeroValue%Depth) .and. &
                (depth <= depthAtBreadthMax)) then
            outvalue = depth * breadthMax

        elseif ((depth > depthAtBreadthMax) .and. &
                (depth < fullDepth)) then
            emptyDepth = max(fulldepth - depth, zeroR) 
            emptyTheta = twoR * acos(oneR - emptyDepth / rTop)
            emptyArea  = onehalfR * (rTop**2) * (emptyTheta - sin(emptyTheta))
            outvalue   = fullArea - emptyArea  

        else
            outvalue   = fullArea
        endif

    end function llgeo_mod_basket_area_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function llgeo_rectangular_closed_area_from_depth_singular &
            (indx, depth) result (outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for rectangular cross section of a single element
        !% The input indx is the row index in full data 2D array.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: indx
            real(8), intent(in) :: depth
            real(8), pointer ::  breadth, fulldepth
        !%------------------------------------------------------------------
        !% Aliases:
            breadth     =>   elemR(indx,er_BreadthMax)
            fulldepth   =>   elemR(indx,er_FullDepth)
        !%------------------------------------------------------------------

        if (depth <= setting%ZeroValue%Depth) then
            outvalue = setting%ZeroValue%Area

        elseif ((depth > setting%ZeroValue%Depth) .and. (depth < fulldepth)) then
            outvalue = depth * breadth
        else
            outvalue = fulldepth * breadth
        end if

    end function llgeo_rectangular_closed_area_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function llgeo_rect_round_area_from_depth_singular &
            (indx, depth) result (outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for rectangular_round cross section of a single element
        !% The input indx is the row index in full data 2D array.
        !%-----------------------------------------------------------------
            integer, intent(in) :: indx
            real(8), intent(in) :: depth
            real(8), pointer :: yBot, rBot, aBot, breadth, fulldepth, fullArea
            real :: theta
        !%------------------------------------------------------------------
            breadth   => elemR(indx,er_BreadthMax)
            yBot      => elemSGR(indx,esgr_Rectangular_Round_Ybot) 
            aBot      => elemSGR(indx,esgr_Rectangular_Round_Abot)
            rBot      => elemSGR(indx,esgr_Rectangular_Round_Rbot)
            fulldepth => elemR(indx,er_FullDepth)
            fullArea  => elemR(indx,er_FullArea)
        !%------------------------------------------------------------------

        if (depth <= setting%ZeroValue%Depth) then
            outvalue = setting%ZeroValue%Area

        elseif ((depth > setting%ZeroValue%Depth) .and. (depth <= yBot)) then
            theta    = twoR * acos(oneR - depth / rBot)
            outvalue = onehalfR * (rBot ** twoI) * (theta - sin(theta))  

        elseif ((depth > yBot) .and. (depth < fulldepth)) then
            outvalue = aBot + (depth - yBot) * breadth

        else
            !% --- truncate area at full depth
            !outvalue = aBot + (fulldepth - yBot) * breadth
            outvalue = fullArea
        endif

    end function llgeo_rect_round_area_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function llgeo_rectangular_triangular_area_from_depth_singular &
            (indx, depth) result (outvalue)
        !%-------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for rectangular_triangular cross section of a single element
        !% The input indx is the row index in full data 2D array.
        !%------------------------------------------------------------------
            integer, intent(in) :: indx
            real(8), intent(in) :: depth
            real(8), pointer :: bottomDepth, bottomSlope, bottomArea, breadth
            real(8), pointer :: fulldepth, fullArea
        !%------------------------------------------------------------------
            bottomDepth => elemSGR(indx,esgr_Rectangular_Triangular_BottomDepth) 
            bottomSlope => elemSGR(indx,esgr_Rectangular_Triangular_BottomSlope)
            bottomArea  => elemSGR(indx,esgr_Rectangular_Triangular_BottomArea)
            breadth     => elemR(indx,er_Breadthmax)
            fulldepth   => elemR(indx,er_FullDepth)
            fullArea    => elemR(indx,er_FullArea)
        !%-----------------------------------------------------------------
        
        if (depth <= setting%ZeroValue%Depth) then    
            outvalue = setting%ZeroValue%Area

        elseif ((depth > setting%ZeroValue%Depth) .and. (depth <= bottomDepth)) then
            outvalue = depth * depth * bottomSlope

        elseif ((depth > bottomDepth) .and. (depth < fulldepth)) then
            outvalue = bottomArea + (depth - bottomDepth) * breadth   

        else
            !% --- truncate at full area
            !outvalue = bottomArea + (fulldepth - bottomDepth) * breadth   
            outvalue = fullArea
        endif

    end function llgeo_rectangular_triangular_area_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    ! real(8) function llgeo_circular_area_from_depth_singular &
    !     (indx, depth) result (outvalue)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Computes area from known depth for circular cross section of a single element
    !     !% The input indx is the row index in full data 2D array.
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: indx
    !         real(8), intent(in) :: depth
    !         real(8), pointer    :: fullArea, fulldepth
    !         real(8) :: YoverYfull
    !     !%------------------------------------------------------------------
    !     !% Aliases
    !         fullArea   => elemR(indx,er_FullArea)
    !         fulldepth  => elemR(indx,er_FullDepth)
    !     !%------------------------------------------------------------------

    !     !% --- find Y/Yfull
    !     YoverYfull = depth / fulldepth

    !     !% --- prevent overfull
    !     YoverYfull = min(YoverYfull, oneR)

    !     !% --- prevent underfull
    !     YoverYfull = max(YoverYfull, setting%ZeroValue%Depth/fulldepth)

    !     !% --- get A/Afull from the lookup table using Y/Yfull and
    !     !%     unnormalize with full area
    !     outvalue = fullArea * xsect_table_lookup_singular (YoverYfull, ACirc) 


    ! end function llgeo_circular_area_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    ! real(8) function llgeo_semi_circular_area_from_depth_singular &
    !         (indx, depth) result (outvalue)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Computes area from known depth for semi_circular cross section of a single element
    !     !% The input indx is the row index in full data 2D array.
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: indx
    !         real(8), intent(in) :: depth
    !         real(8), pointer    :: fullArea, fulldepth
    !         real(8)             :: AoverAfull, YoverYfull
    !     !%-----------------------------------------------------------------
    !     !% Aliases
    !         fullArea   => elemR(indx,er_FullArea)
    !         fulldepth  => elemR(indx,er_FullDepth)
    !     !%-----------------------------------------------------------------

    !     !% find Y/Yfull
    !     YoverYfull(indx) = depth / fulldepth(indx)

    !     !% --- prevent overfull
    !     YoverYfull = min(YoverYfull, oneR)

    !     !% --- prevent underfull
    !     YoverYfull = max(YoverYfull, setting%ZeroValue%Depth/fulldepth)

    !     !% --- lookup area and unnormalize
    !     outvalue = fullArea * xsect_table_lookup_singular (YoverYfull, ASemiCircular)

    ! end function semi_circular_area_from_depth_singular
!%
!%==========================================================================
!% TOPWIDTH FROM DEPTH -- OPEN CHANNELS
!%==========================================================================
!%
    pure function llgeo_parabolic_topwidth_from_depth_pure &
            (thisP, depth) 
        !%------------------------------------------------------------------
        !% Description:
        !% Computes topwidth from known depth for parabolic cross section
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)  ! may be a packed array of indexes
            real(8), intent(in) :: depth(:)
            real(8), dimension(size(thisP)) :: llgeo_parabolic_topwidth_from_depth_pure
        !%------------------------------------------------------------------

        llgeo_parabolic_topwidth_from_depth_pure               &
                 = twoR * elemSGR(thisP,esgr_Parabolic_Radius) &
                         * sqrt(depth)

    end function llgeo_parabolic_topwidth_from_depth_pure
!%
!%==========================================================================
!%==========================================================================
!%
    pure function llgeo_powerfunction_topwidth_from_depth_pure &
            (thisP, depth) 
        !%------------------------------------------------------------------
        !% Description:
        !% Computes topwidth from known depth for powerfunction cross section
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)  ! may be a packed array of indexes
            real(8), intent(in) :: depth(:)
            real(8), dimension(size(thisP)) :: llgeo_powerfunction_topwidth_from_depth_pure
        !%------------------------------------------------------------------

        llgeo_powerfunction_topwidth_from_depth_pure = nullvalueR  !% STUB  HACK

    end function llgeo_powerfunction_topwidth_from_depth_pure
!%
!%==========================================================================    
!%==========================================================================
!%    
    pure function llgeo_rectangular_topwidth_from_depth_pure &
            (thisP, depth) 
        !%------------------------------------------------------------------
        !% Description:
        !% Computes topwidth from known depth for rectangular cross section
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)  ! may be a packed array of indexes
            real(8), intent(in) :: depth(:)
            real(8), dimension(size(thisP)) :: llgeo_rectangular_topwidth_from_depth_pure
        !%------------------------------------------------------------------

            llgeo_rectangular_topwidth_from_depth_pure = elemR(thisP,er_BreadthMax)

    end function llgeo_rectangular_topwidth_from_depth_pure
!%
!%==========================================================================
!%==========================================================================
!%
    pure function llgeo_trapezoidal_topwidth_from_depth_pure &
            (thisP, depth) 
        !%------------------------------------------------------------------
        !% Description:
        !% Computes topwidth from known depth for rectangular cross section
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)  ! may be a packed array of indexes
            real(8), intent(in) :: depth(:)
            real(8), dimension(size(thisP)) :: llgeo_trapezoidal_topwidth_from_depth_pure
        !%------------------------------------------------------------------

            llgeo_trapezoidal_topwidth_from_depth_pure            &
                 = elemSGR(thisP,esgr_Trapezoidal_Breadth)        &
                 + depth                          &
                 * (   elemSGR(thisP,esgr_Trapezoidal_LeftSlope)  &
                     + elemSGR(thisP,esgr_Trapezoidal_RightSlope) &
                    )

    end function llgeo_trapezoidal_topwidth_from_depth_pure
!%
!%==========================================================================
!%==========================================================================
!%
    pure function llgeo_triangular_topwidth_from_depth_pure &
            (thisP, depth) 
        !%------------------------------------------------------------------
        !% Description:
        !% Computes topwidth from known depth for triangular cross section
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)  ! may be a packed array of indexes
            real(8), intent(in) :: depth(:)
            real(8), dimension(size(thisP)) :: llgeo_triangular_topwidth_from_depth_pure
        !%------------------------------------------------------------------

        llgeo_triangular_topwidth_from_depth_pure &
             = twoR * elemSGR(thisP,esgr_Triangular_Slope) &
                        * depth

    end function llgeo_triangular_topwidth_from_depth_pure
!%
!%==========================================================================
!% TOPWIDTH FROM DEPTH --- CLOSED CONDUITS
!%==========================================================================
!%
    real(8) function llgeo_filled_circular_topwidth_from_depth_singular &
            (indx,depth) result (outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth for a filled_circular cross section of a single element
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: indx
            real(8), intent(in) :: depth
            real(8), pointer    :: fullDepth, sedimentDepth, breadthMax
            real(8), pointer    :: totalFullDepth
            real(8)             :: YoverYfull
        !%-------------------------------------------------------------------
        !% Aliases
            fullDepth      =>  elemR(indx,er_fullDepth)
            breadthMax     =>  elemR(indx,er_BreadthMax)
            sedimentDepth  =>  elemR(indx,er_SedimentDepth)
            totalFullDepth => elemSGR(indx,esgr_Filled_Circular_TotalPipeDiameter)
        !%-------------------------------------------------------------------

        !% --- normalized depth
        YoverYfull = (depth + sedimentDepth) / totalFullDepth

        !% --- prevent overfull
        YoverYfull = min(YoverYfull, oneR)

        !% --- prevent underfull (must be above sediment depth)
        YoverYfull = max(YoverYfull, (sedimentDepth + setting%ZeroValue%Depth)/totalFullDepth )

        !% --- get normalized topwidth and unnormalize
        outvalue =  breadthMax * xsect_table_lookup_singular (YoverYfull, TCirc) 

        !% --- if topwidth < zero value, set it to zerovalue
        outvalue = max(outvalue, setting%ZeroValue%Topwidth)

    end function llgeo_filled_circular_topwidth_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function llgeo_mod_basket_topwidth_from_depth_singular &
            (indx, depth) result (outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth for a mod_basket cross section of a single element
        !%-------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: indx 
            real(8), intent(in) :: depth
            real(8), pointer    :: breadthMax, fullDepth
            real(8), pointer    :: depthAtBreadthMax, rTop
            real(8)             :: emptyDepth, depthlimit
        !%-------------------------------------------------------------------
        !% Aliases:
            fullDepth         => elemR(indx,er_FullDepth)
            depthAtBreadthMax => elemR(indx,er_DepthAtBreadthMax)
            breadthMax        => elemR(indx,er_BreadthMax)
            rTop              => elemSGR(indx,esgr_Mod_Basket_Rtop) 
        !%--------------------------------------------------------------------
        
        depthlimit = fullDepth * setting%Discretization%FullConduitTopwidthDepthFraction

        if (depth <=  setting%ZeroValue%Depth) then
            outvalue = setting%ZeroValue%Topwidth

        elseif ((depth > setting%ZeroValue%Depth) .and. &
                (depth <= depthAtBreadthMax) ) then
            outvalue = breadthMax

        elseif ((depth > depthAtBreadthMax) .and. &
                (depth <= depthlimit) ) then
            emptyDepth = fullDepth - depth
            outvalue   = twoR * sqrt(emptyDepth * (twoR * rTop - emptyDepth))

        else 
            !% --- large depths are limited to prevent topwidth going to zero
            emptyDepth = fullDepth - depthlimit
            outvalue   = twoR * sqrt(emptyDepth * (twoR * rTop - emptyDepth))
        endif

        !% --- if topwidth <= zero value, set it to zerovalue
        outvalue = max(outvalue, setting%ZeroValue%Topwidth)

    end function llgeo_mod_basket_topwidth_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function llgeo_rectangular_closed_topwidth_from_depth_singular &
            (indx, depth) result (outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth for a rectangular cross section of a single element
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: indx 
            real(8), intent(in) :: depth
        !%------------------------------------------------------------------

        if (depth <= setting%ZeroValue%Depth) then
            outvalue = setting%ZeroValue%Topwidth
        else
            !% --- rectangular all the way to top (use breadth at full)
            outvalue = elemSGR(indx,esgr_Rectangular_Breadth)
        end if

    end function llgeo_rectangular_closed_topwidth_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function llgeo_rect_round_topwidth_from_depth_singular &
            (indx, depth) result (outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth for a rectangular_round cross section of a single element
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: indx 
            real(8), intent(in) :: depth
            real(8), pointer :: breadth, yBot, rBot
        !%------------------------------------------------------------------
        !% Aliases:
            yBot        => elemSGR(indx,esgr_Rectangular_Round_Ybot)
            rBot        => elemSGR(indx,esgr_Rectangular_Round_Ybot) 
            breadth     =>   elemR(indx,er_BreadthMax)
        !%------------------------------------------------------------------

        if (depth <= setting%ZeroValue%Depth) then
            outvalue = setting%ZeroValue%Topwidth

        elseif ((depth > setting%ZeroValue%Depth) .and. (depth < yBot)) then
            !% --- bottom circular section
            outvalue = twoR * sqrt(depth * (twoR * rBot - depth))

        else
            !% --- rectangular all the way to top
            outvalue = breadth
        end if
        
    end function llgeo_rect_round_topwidth_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function llgeo_rectangular_triangular_topwidth_from_depth_singular &
            (indx, depth) result (outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth for a rectangular_triangular cross section of a single element
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: indx 
            real(8), intent(in) :: depth
            real(8), pointer :: bottomDepth, bottomSlope, breadth
        !%------------------------------------------------------------------
        !% Aliases:
            bottomSlope => elemSGR(indx,esgr_Rectangular_Triangular_BottomSlope)
            bottomDepth => elemSGR(indx,esgr_Rectangular_Triangular_BottomDepth) 
            breadth     =>   elemR(indx,er_BreadthMax)
        !%-------------------------------------------------------------------
         
        if (depth <= setting%ZeroValue%Depth) then
            outvalue = setting%ZeroValue%Topwidth    

        elseif (( depth > setting%ZeroValue%Depth) .and. (depth < bottomDepth)) then
            outvalue = twoR * bottomSlope * depth

        else
            !% --- rectangular all the way to top
            outvalue = breadth
        endif

    end function llgeo_rectangular_triangular_topwidth_from_depth_singular
!%
!%========================================================================== 
!%==========================================================================
!%
    ! real(8) function llgeo_circular_topwidth_from_depth_singular &
    !     (indx,depth) result (outvalue)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Computes the topwidth for a circular cross section of a single element
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: indx
    !         real(8), intent(in) :: depth
    !         real(8), pointer    :: maxTopwidth, fulldepth
    !         real(8) :: YoverYfull
    !     !%------------------------------------------------------------------
    !     !% Aliases
    !         fulldepth   => elemR(indx,er_FullDepth)
    !         maxTopwidth => elemSGR(indx,esgr_Circular_Diameter)
    !     !%------------------------------------------------------------------

    !     !% --- find Y/Yfull
    !     YoverYfull = depth / fulldepth

    !     !% --- prevent overfull
    !     YoverYfull = min(YoverYfull, oneR)

    !     !% --- prevent underfull
    !     YoverYfull = max(YoverYfull, setting%ZeroValue%Depth/fulldepth)

    !     !% --- retriving T/Tmax from the lookup table using Y/Yfull and
    !     !%     unnormalizing with max topwidth
    !     outvalue = maxTopwidth * xsect_table_lookup_singular (YoverYfull, TCirc) !% 20220506brh removed NTCirc

    !     !% if topwidth <= zero, set it to zerovalue
    !     outvalue = max(outvalue, setting%ZeroValue%Topwidth)

    !     !% --- if topwidth < full value, set it to full value
    !     outvalue = max(outvalue, fullTopWidth)
    

    ! end function llgeo_circular_topwidth_from_depth_singular
!%
!%==========================================================================   
!%==========================================================================
!%
    ! real(8) function llgeo_semi_circular_topwidth_from_depth_singular &
    !     (indx,depth) result (outvalue)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Computes the topwidth for a semi_circular cross section of a single element
    !     !%------------------------------------------------------------------
    !         integer, intent(in) :: indx
    !         real(8), intent(in) :: depth
    !         real(8), pointer    :: fulldepth, maxTopwidth, fullTopwidth
    !         real(8)             :: YoverYfull
    !     !%-----------------------------------------------------------------
    !     fulldepth    => elemR(indx,er_FullDepth)
    !     maxTopwidth  => elemR(indx,er_BreadthMax)
    !     fullTopwidth => elemR(indx,er_FullTopwidth)
    !     !%-----------------------------------------------------------------

    !     !% find Y/Yfull
    !     YoverYfull(indx) = depth / fulldepth(indx)

    !     !% --- prevent overfull
    !     YoverYfull = min(YoverYfull, oneR)

    !     !% --- prevent underfull
    !     YoverYfull = max(YoverYfull, setting%ZeroValue%Depth/fulldepth)

    !     !% --- lookup table and unnoralize by maximum topwidth
    !     outvalue = maxTopwidth * xsect_table_lookup_singular (YoverYfull, TSemiCircular) 

    !     !% --- if topwidth <= zero, set it to zerovalue
    !     outvalue = max(outvalue, setting%ZeroValue%Topwidth)

    !     !% --- if topwidth < full value, set it to full value
    !     outvalue = max(outvalue, fullTopwidth)
    

    ! end function llgeo_semi_circular_topwidth_from_depth_singular
!%
!%==========================================================================
!% PERIMETER FROM DEPTH -- OPEN CHANNEL
!%==========================================================================
!% 
    pure function llgeo_parabolic_perimeter_from_depth_pure  &
            (thisP, depth)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes perimeter from known depth for parabolic cross section
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)  ! must be a packed array of indexes
            real(8), intent(in) :: depth(:)
            real(8), dimension(size(thisP)) :: temp1, temp2
            real(8), dimension(size(thisP)) :: llgeo_parabolic_perimeter_from_depth_pure
        !%------------------------------------------------------------------

        !% -- temporary arrays for clarity
        temp1 = twoR * sqrt(depth) &
                            / elemSGR(thisP,esgr_Parabolic_Radius)   

        !% -- temporary arrays for clarity
        temp2 = sqrt(oneR + (temp1**2))                     

        llgeo_parabolic_perimeter_from_depth_pure &
                = onehalfR * (elemSGR(thisP,esgr_Parabolic_Radius)**2)   &
                           * (                                           &
                               temp1 * temp2   + log(temp1 + temp2)      &
                             )                   

    end function llgeo_parabolic_perimeter_from_depth_pure                   
!%
!%==========================================================================
!%==========================================================================
!% 
    pure function llgeo_powerfunction_perimeter_from_depth_pure &
            (thisP, depth) 
        !%------------------------------------------------------------------
        !% Description:
        !% Computes perimeter from known depth for powerfunction cross section
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)  ! must be a packed array of indexes
            real(8), intent(in) :: depth(:)
            real(8), dimension(size(thisP)) :: llgeo_powerfunction_perimeter_from_depth_pure
        !%------------------------------------------------------------------

        llgeo_powerfunction_perimeter_from_depth_pure = nullvalueR  !% STUB HACK

    end function llgeo_powerfunction_perimeter_from_depth_pure                   
!%
!%==========================================================================
!%==========================================================================
!% 
    pure function llgeo_rectangular_perimeter_from_depth_pure &
            (thisP, depth) 
        !%------------------------------------------------------------------
        !% Description:
        !% Computes perimeter from known depth for rectangular cross section
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)  ! must be a packed array of indexes
            real(8), intent(in) :: depth(:)
            real(8), dimension(size(thisP)) :: llgeo_rectangular_perimeter_from_depth_pure
        !%------------------------------------------------------------------

        llgeo_rectangular_perimeter_from_depth_pure &
             = twoR * depth + elemSGR(thisP,esgr_Rectangular_Breadth)

    end function llgeo_rectangular_perimeter_from_depth_pure                   
!%
!%==========================================================================
!%==========================================================================
!% 
    pure function llgeo_trapezoidal_perimeter_from_depth_pure &
            (thisP, depth) 
        !%------------------------------------------------------------------
        !% Description:
        !% Computes perimeter from known depth for trapezoidal cross section
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)  ! must be a packed array of indexes
            real(8), intent(in) :: depth(:)
            real(8), dimension(size(thisP)) :: llgeo_trapezoidal_perimeter_from_depth_pure
        !%------------------------------------------------------------------

        llgeo_trapezoidal_perimeter_from_depth_pure                        &
            = elemSGR(thisP,esgr_Trapezoidal_Breadth)                      &
            + depth                                        &
            * (                                                            &
                sqrt(oneR + elemSGR(thisP,esgr_Trapezoidal_Leftslope )**2) &
               +sqrt(oneR + elemSGR(thisP,esgr_Trapezoidal_RightSlope)**2) &
              )

        !perimeter(thisP) = breadth(thisP) + depth(thisP) * (sqrt(oneR + lslope(thisP)**twoR) &
        !                + sqrt(oneR + rslope(thisP)**twoR))

    end function llgeo_trapezoidal_perimeter_from_depth_pure                   
!%
!%==========================================================================
!%==========================================================================
!% 
    pure function llgeo_triangular_perimeter_from_depth_pure &
            (thisP, depth) 
        !%------------------------------------------------------------------
        !% Description:
        !% Computes perimeter from known depth for triangular cross section
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:)  ! must be a packed array of indexes
            real(8), intent(in) :: depth(:)
            real(8), dimension(size(thisP)) :: llgeo_triangular_perimeter_from_depth_pure
        !%------------------------------------------------------------------

        llgeo_triangular_perimeter_from_depth_pure                   &
            = twoR * depth                          &
            * sqrt( oneR + elemSGR(thisP,esgr_Triangular_Slope)** twoI )

    end function llgeo_triangular_perimeter_from_depth_pure                   
!%
!%==========================================================================
!% PERIMETER FROM DEPTH -- CLOSED CONDUITS
!%==========================================================================
!%
    ! real(8) function llgeo_circular_perimeter_from_depth_singular  &
    !         (indx, depth) result(outvalue)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Computes the circular conduit perimeter for the given depth on
    !     !% the element idx
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: indx
    !         real(8), intent(in) :: depth
    !         real(8), pointer :: fulldepth, fullarea, fullperimeter
    !         real(8) :: hydRadius, YoverYfull, area

    !     !%------------------------------------------------------------------
    !     !% Aliases
    !         fulldepth     => elemR(indx,er_FullDepth)
    !         fullarea      => elemR(indx,er_FullArea)
    !         fullHydRadius => elemR(indx,er_FullHydRadius)
    !         fullperimeter => elemR(indx,er_FullPerimeter)
    !     !%------------------------------------------------------------------

    !     !% --- find Y/Yfull
    !     YoverYfull = depth / fulldepth

    !     !% --- prevent overfull
    !     YoverYfull = min(YoverYfull, oneR)

    !     !% --- prevent underfull
    !     YoverYfull = max(YoverYfull, setting%ZeroValue%Depth/fulldepth)

    !     !% --- retrieve normalized A/Amax for this depth from lookup table
    !     area = xsect_table_lookup_singular (YoverYfull, ACirc)

    !     !% --- retrive the normalized R/Rmax for this depth from the lookup table
    !     hydradius =  xsect_table_lookup_singular (YoverYfull, RCirc) 

    !     !% --- unnormalize
    !     hydRadius = hydradius * fullHydRadius
    !     area      = area * fullArea

    !     hydRadius = max(hydradius,setting%ZeroValue%Depth)

    !     !% --- get the perimeter by dividing area by hydRadius
    !     !%     limit by the full perimeter
    !     outvalue = min(area / hydRadius, fullperimeter)

    ! end function llgeo_circular_perimeter_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function llgeo_filled_circular_perimeter_from_depth_singular &
        (indx, depth) result (outvalue)   
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth for a filled_circular cross section of a single element
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: indx
            real(8), intent(in) :: depth
            real(8), pointer    :: fullDepth, sedimentDepth, breadthMax
            real(8), pointer    :: sedimentPerimeter, sedimentTopwidth
            real(8), pointer    :: totalFullDepth, totalFullArea, totalFullHydRadius
            real(8)             :: YoverYfull, tempHydRadius, tempArea
        !%-------------------------------------------------------------------
        !% Aliases
            fullDepth          =>   elemR(indx,er_fullDepth)
            breadthMax         =>   elemR(indx,er_BreadthMax)
            sedimentDepth      =>   elemR(indx,er_SedimentDepth)
            sedimentPerimeter  => elemSGR(indx,esgr_Filled_Circular_bottomPerimeter)     
            sedimentTopwidth   => elemSGR(indx,esgr_Filled_Circular_bottomTopWidth)
            totalFullDepth     => elemSGR(indx,esgr_Filled_Circular_TotalPipeDiameter)
            totalFullArea      => elemSGR(indx,esgr_Filled_Circular_TotalPipeArea)
            totalFullHydRadius => elemSGR(indx,esgr_Filled_Circular_TotalPipeHydRadius)
        !%-------------------------------------------------------------------
        !% --- normalized depth
        YoverYfull = (depth + sedimentDepth) / totalFullDepth

        !% --- prevent overfull
        YoverYfull = min(YoverYfull, oneR)

        !% --- prevent underfull (must be above sediment depth)
        YoverYfull = max(YoverYfull, (sedimentDepth + setting%ZeroValue%Depth)/totalFullDepth ) 
        
        !% --- get hydraulic radius of sediment + flow a
        tempHydRadius =  totalFullHydRadius * xsect_table_lookup_singular (YoverYfull, RCirc)

        !% --- ensure no zeros in normalized hydraulic radius
        tempHydRadius = max(temphydRadius,setting%ZeroValue%Depth) 

        !% --- compute area of sediment + flow
        tempArea  = totalFullArea * xsect_table_lookup_singular (YoverYfull, ACirc) 

        !% --- perimeter of the sediment + flow is the Area/Hyd Radius of 
        !%     sediment + flow, then with the sediment perimeter subtracted
        !%     and the sediment topwidth added

        outvalue = (tempArea / tempHydRadius) &
                    - sedimentPerimeter + sedimentTopwidth

        !% ensure perimeter is greater than zero
        outvalue = max(outvalue, setting%ZeroValue%Topwidth)            

    end function llgeo_filled_circular_perimeter_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function llgeo_mod_basket_perimeter_from_depth_singular &
            (indx, depth) result (outvalue)  
        !%------------------------------------------------------------------
        !% Description:
        !% Computes wetted perimeter from known depth for a mod_basket cross section of
        !% a single element 
        !%------------------------------------------------------------------
        !% Deckarations:
            integer, intent(in) :: indx
            real(8), intent(in) :: depth
            real(8), pointer :: fullDepth, breadth, yBreadthMax, rTop, thetaTop
            real(8), pointer :: fullPerimeter
            real(8) :: emptyDepth, emptyTheta
        !%------------------------------------------------------------------
        !% Aliases
            fullDepth     => elemR(indx,er_FullDepth)
            fullPerimeter => elemR(indx,er_FullPerimeter)
            yBreadthMax   => elemR(indx,er_DepthAtBreadthMax)
            breadth       => elemR(indx,er_BreadthMax)
            rTop        => elemSGR(indx,esgr_Mod_Basket_Rtop)
            thetaTop    => elemSGR(indx,esgr_Mod_Basket_ThetaTop) 
        !%------------------------------------------------------------------
        
        if (depth <= setting%ZeroValue%Depth) then
            outvalue = setting%ZeroValue%Topwidth

        elseif ((depth > setting%ZeroValue%Depth) .and. (depth <= yBreadthMax)) then
            outvalue = twoR * depth + breadth

        elseif ((depth > yBreadthMax) .and. (depth < fullDepth)) then
            !% find height of empty area
            emptyDepth = max(fullDepth - depth, zeroR)
            !% find angle of circular arc corresponding to this height
            emptyTheta = twoR * acos(oneR - emptyDepth / rTop)
            !% find perimeter of wetted portion of circular arc
            outvalue  = (thetaTop - emptyTheta) * rTop
            !% add on wetted perimeter of bottom rectangular area
            outvalue  = outvalue + twoR * yBreadthMax + breadth
        else 
            outvalue = fullPerimeter                 
        endif

    end function llgeo_mod_basket_perimeter_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function llgeo_rectangular_closed_perimeter_from_depth_singular &
            (indx, depth) result (outvalue)
        !%  
        !%------------------------------------------------------------------
        !% Description:
        !% Computes wetted perimeter from known depth for a rectangular cross section of
        !% a single element 
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: indx
            real(8), intent(in) :: depth
            real(8), pointer    ::  breadth, fulldepth
        !%-----------------------------------------------------------------
        !% Aliases:
            breadth       =>   elemR(indx,er_BreadthMax)
            fulldepth     =>   elemR(indx,er_FullDepth)
        !%-----------------------------------------------------------------

        if (depth < setting%ZeroValue%Depth) then
            outvalue = setting%ZeroValue%Topwidth

        elseif ((depth > setting%ZeroValue%Depth) .and. (depth < fulldepth)) then
            outvalue = twoR * depth + breadth

        else 
            !% --- full depth value
            outvalue = twoR * (fulldepth + breadth)

        end if

    end function llgeo_rectangular_closed_perimeter_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function llgeo_rect_round_perimeter_from_depth_singular &
            (indx, depth) result (outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes wetted perimeter from known depth for a rectangular_round cross section of
        !% a single element 
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: indx
            real(8), intent(in) :: depth
            real(8), pointer :: breadth
            real(8), pointer :: fullDepth, yBot, rBot
            real(8) :: theta
        !%------------------------------------------------------------------
        !% Aliases
            fullDepth   => elemR(indx,er_FullDepth)
            breadth     => elemR(indx,er_BreadthMax)
            yBot        => elemSGR(indx,esgr_Rectangular_Round_Ybot)
            rBot        => elemSGR(indx,esgr_Rectangular_Round_Rbot)  
        !%------------------------------------------------------------------
        
        if (depth <= setting%ZeroValue%Depth) then
            outvalue = setting%ZeroValue%Topwidth

        elseif ((depth > setting%ZeroValue%Depth) .and. (depth < yBot)) then
            !% --- bottom circular section
            theta    = twoR * acos(oneR - depth / rBot)
            outvalue = rBot * theta

        elseif ((depth > yBot) .and. (depth < fulldepth)) then
            !% --- top rectangular section
            theta    = twoR * asin(breadth / twoR / rBot)
            outvalue = rBot * theta + twoR * (depth - yBot)

        else 
            !% --- at or above full depth perimeter includes top
            theta    = twoR * asin(breadth / twoR / rBot)
            outvalue = rBot * theta + twoR * (fulldepth - yBot) + breadth
        end if

    end function llgeo_rect_round_perimeter_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function llgeo_rectangular_triangular_perimeter_from_depth_singular &
            (indx, depth) result (outvalue) 
        !%------------------------------------------------------------------
        !% Description:
        !% Computes wetted perimeter from known depth for a rectangular_triangular cross section of
        !% a single element 
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: indx
            real(8), intent(in) :: depth
            real(8), pointer :: bottomDepth, bottomSlope, breadth, fulldepth
        !%------------------------------------------------------------------
        !% Aliases:
            fullDepth   => elemR(indx,er_FullDepth)
            breadth     => elemR(indx,er_BreadthMax)
            bottomSlope => elemSGR(indx,esgr_Rectangular_Triangular_BottomSlope)
            bottomDepth => elemSGR(indx,esgr_Rectangular_Triangular_BottomDepth) 
        !%------------------------------------------------------------------
        
        if (depth <= setting%ZeroValue%Depth) then
            outvalue = setting%ZeroValue%Topwidth

        elseif ((depth > setting%ZeroValue%Depth) .and. (depth <= bottomDepth) ) then
            outvalue = twoR * depth * sqrt(oneR + bottomSlope**twoI)

        elseif ((depth >  bottomDepth) .and. (depth < fulldepth)) then
            outvalue = twoR * bottomDepth * sqrt(oneR + bottomSlope**twoI) &   !triangular section
                        + twoR * (depth - bottomDepth)                      !rectangular section

        else 
            !% --- at or above full depth perimeter includes top
            outvalue = twoR * bottomDepth * sqrt(oneR + bottomSlope**twoI) &   
                        + twoR * (fullDepth - bottomDepth) + breadth                   
        endif

    end function llgeo_rectangular_triangular_perimeter_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
    ! real(8) function llgeo_semi_circular_perimeter_from_depth_singular &
    !         (indx, depth) result(outvalue)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Computes the semi_circular conduit perimeter for the given depth on
    !     !% the element idx
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: indx
    !         real(8), intent(in) :: depth
    !         real(8), pointer :: fulldepth, fullarea, fullperimeter
    !         real(8) :: hydRadius, YoverYfull, AoverAfull, area, sf
    !         real(8), parameter :: sfparam = 0.2946d0
    !     !%------------------------------------------------------------------
    !     !% Aliases
    !         fulldepth     => elemR(indx,er_FullDepth)
    !         fullarea      => elemR(indx,er_FullArea)
    !         fullperimeter => elemR(indx,er_FullPerimeter)
    !     !%------------------------------------------------------------------

    !     !% --- set normalized depth   
    !     YoverYfull = depth / fulldepth

    !     !% --- prevent overfull
    !     YoverYfull = min(YoverYfull, oneR)

    !     !% --- prevent underfull
    !     YoverYfull = max(YoverYfull, setting%ZeroValue%Depth/fulldepth)

    !     !% --- retrieve normalized area for this depth from lookup table
    !     AoverAfull = xsect_table_lookup_singular (YoverYfull, ASemiCircular)

    !     !% --- retrieve normalized sectionfactor for this depth from lookup table
    !     sf = xsect_table_lookup_singular (AoverAfull, SSemiCircular)

    !     !% --- unnormalize
    !     sf   = (fullArea * (sfparam * fullDepth) ** twoThirdR) * sf
    !     area = AoverAfull * fullarea

    !     !% retrieve hyrdaulic radius from section factor
    !     hydRadius = (sF / area) ** threehalfR

    !     !% --- get the perimeter by dividing area by hydRadius
    !     outvalue = min(area / hydRadius, fullperimeter)

    ! end function llgeo_semi_circular_perimeter_from_depth_singular
!%    
!%==========================================================================

!% END MODULE
!%==========================================================================
end module geometry_lowlevel