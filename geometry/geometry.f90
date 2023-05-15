module geometry

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use define_xsect_tables
    use geometry_lowlevel
    use preissmann_slot
    use parabolic_channel
    use powerfunction_channel
    use irregular_channel
    use rectangular_channel
    use rectangular_conduit
    use trapezoidal_channel
    use triangular_channel
    use arch_conduit
    use basket_handle_conduit
    use catenary_conduit
    use circular_conduit
    use egg_shaped_conduit
    use filled_circular_conduit
    use gothic_conduit
    use horiz_ellipse_conduit
    use horse_shoe_conduit
    use mod_basket_conduit
    use rectangular_round_conduit
    use rectangular_triangular_conduit
    use semi_circular_conduit
    use semi_elliptical_conduit
    use vert_ellipse_conduit
    use storage_geometry
   
    use xsect_tables
    use adjust
    use utility_profiler
    use utility_crash

    ! use utility_unit_testing, only: util_utest_CLprint


    implicit none

!%-----------------------------------------------------------------------------
!% Description:
!% Geometry computations
!%

    private

    !public :: geometry_toplevel
    public :: geometry_toplevel_CC
    public :: geometry_toplevel_JMJB
    public :: geo_assign_JB_from_head
    public :: geo_JM_depth_area_from_volume
    public :: geo_common_initialize
    public :: geo_sectionfactor_from_depth_singular
    public :: geo_Qcritical_from_depth_singular
    public :: geo_critical_value_singular
    public :: geo_normaldepth_singular
    !public :: geo_topwidth_from_depth_by_type_CC
    ! public :: geo_hyddepth_from_area_and_topwidth_singular
    public :: geo_topwidth_from_depth_singular
    public :: geo_area_from_depth_singular
    public :: geo_perimeter_from_depth_singular
    !public :: geo_elldepth_pure
    !public :: geometry_table_initialize
    public :: geo_depth_from_volume_by_type_allCC
    public :: geo_depth_from_volume_by_element_CC
    public :: geo_depth_from_head
    !public :: geo_depth_from_volume_by_type_JM
    public :: geo_ZeroDepth_from_depth

    public :: geo_plan_area_from_volume_JM
    public :: geo_depth_from_volume_JM

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine geometry_toplevel_CC (                          &
            thisP, npackP, thisP_Open, npackP_Open,            &
            thisP_Closed, npackP_Closed, isSingularYN, isALLYN)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes geometry on channel/conduit (CC) elements for the  
        !% time-marching scheme of whichTM (ETM, AC, ALLtm)
        !%------------------------------------------------------------------
        !% Declarations
            !% packed data
            integer,  intent(in) :: thisP(:), thisP_Open(:), thisP_Closed(:)
            integer,  intent(in) :: npackP, npackP_Open, npackP_Closed
            !% singular element
            logical, intent(in) :: isSingularYN, isAllYN
        !%------------------------------------------------------------------
        !% Aliases
        !%------------------------------------------------------------------
        !% Preliminary
            if (npackP < 1) return
        !%------------------------------------------------------------------

        !% --- PREISSMAN SLOT    
        !% --- Handle Preissmann Slot for closed CC elements
        !%     with this time march type.
        if (npackP_Closed > 0) then
            call slot_CC_ETM (thisP_Closed)
        end if

        ! print *, 'in geometry AAAA',elemR(72,er_Depth), elemR(72,er_EllDepth)

        !% --- DEPTH
        !%     compute the depth on all elements of CC based on geometry.
        !%     If surcharged, this call returns the full depth of a closed conduit 
        !%     without adding Preissmann Slot depth.
        if (isAllYN) then
            call geo_depth_from_volume_by_type_allCC (elemPGetm, npack_elemPGetm, col_elemPGetm)
        else
            call geo_depth_from_volume_by_element_CC (thisP, npackP)
        end if

        ! print *, 'in geometry BBBB',elemR(72,er_Depth), elemR(72,er_EllDepth)

        !% --- ZERO DEPTH CC
        !%     reset all zero or near-zero depths in CC
        !%     Arguably this should not be needed as the individual depth computations
        !%     in geo_depth_from_volume_by_type_CC should use the zerovalues as minimums
        !%     but this needs to be confirmed.
        call adjust_limit_by_zerovalues &
            (er_Depth, setting%ZeroValue%Depth, thisP, .false.)

            ! print *, 'in geometry CCCC',elemR(5,er_Depth), elemR(5,er_EllDepth)

        !% --- PIEZOMETRIC HEAD
        !%     compute the head on all elements of CC
        !%     This sets head consistent with depth computed in geo_depth_from_volume
        !%     Head is strictly limited to the max depth + zbottom so it does not
        !%     include surcharge effects     
        elemR(thisP,er_Head) = llgeo_head_from_depth_pure &
                                    (thisP, elemR(thisP,er_Depth))

                                    ! print *, 'in geometry DDDD',elemR(5,er_Depth), elemR(5,er_EllDepth)

        !% --- OPEN CHANNEL OVERFLOW
        !%     Compute the overflow lost for CC open channels above
        !%     their maximum volume (no ponding allowed from open CC). 
        !%     Note that overflow or ponding for JM elements is handled 
        !%     in slot_JM_ETM.
        !%     Note, this is NOT standard in EPA-SWMM
        !% 20230508 This needs careful checking and integration
        if (npackP_Open > 0) then
            if (setting%Discretization%AllowChannelOverflowTF) then
                call geo_overflow_openchannels (thisP_Open)
            end if
        end if

        ! print *, 'in geometry DDDD',elemR(5,er_Depth), elemR(5,er_EllDepth)

        !% --- PREISSMAN SLOT VOLUME LIMIT CLOSED CONDUIT CC
        !%     limit the volume in closed element (CC) to the full volume
        !%     Note the excess volume has already been stored in the Preissman Slot
        if (npackP_Closed > 0) then
            call geo_volumelimit_closed (thisP_Closed)
        end if

        ! print *, 'in geometry EEEE',elemR(5,er_Depth), elemR(5,er_EllDepth)

        !% --- CROSS-SECTIONAL AREA
        !%     compute area from volume for CC
        !%     For closed conduits this is based on the volume limited by full volume.
        !%     For open channels the volume limit depends on if AllowChanneOverflowTF is false.
        elemR(thisP,er_Area) = llgeo_area_from_volume_pure(thisP,elemR(thisP,er_Volume))
        elemR(thisP,er_Area) = max(elemR(thisP,er_Area),setting%ZeroValue%Area)

        ! print *, 'in geometry FFFF',elemR(5,er_Depth), elemR(5,er_EllDepth)

        !% --- TOPWIDTH CC
        !%     compute topwidth from depth for all CC
        !%     Note: volume is limited to full depth UNLESS AllowChannelOverflowTF is false

        if (isAllYN) then
            call geo_topwidth_from_depth_by_type_allCC (elemPGetm, npack_elemPGetm, col_elemPGetm)
        else
            call geo_topwidth_from_depth_by_element_CC (thisP, npackP)
        end if

        ! print *, 'in geometry GGGG',elemR(5,er_Depth), elemR(5,er_EllDepth)

        !% --- PERIMETER AND HYDRAULIC RADIUS CC
        !%     compute hydraulic radius and perimeter
        !%     note these two are done together because for analytical cross-sections
        !%     we have equations for perimeter, whereas lookup cross-sections
        !%     have tables for hydraulic radius.
        if (isAllYN) then
            call geo_perimeter_and_hydradius_from_depth_by_type_CC (elemPGetm, npack_elemPGetm, col_elemPGetm) 
        else 
            call geo_perimeter_and_hydradius_from_depth_by_element_CC (thisP, npackP)
        end if

        ! print *, 'in geometry HHHH',elemR(5,er_Depth), elemR(5,er_EllDepth)

        !% --- ELLDEPTH MODIFIED HYDRAULIC DEPTH
        !%     the modified hydraulic depth "ell" is used for 
        !%     for Froude number computations on all CC elements
        call geo_elldepth_from_head_CC (thisP)

        ! print *, 'in geometry IIII',elemR(5,er_Depth), elemR(5,er_EllDepth)

        !% ---- ADJUST SLOT 
        !%      make adjustments for slots on closed elements only
        !%     These add slot values to volume, depth, head
        if (npackP_Closed > 0) then
            call slot_CC_adjustments (thisP_Closed)
        end if
        
        ! print *, 'in geometry JJJJ',elemR(5,er_Depth), elemR(5,er_EllDepth)

        !% --- check for crashpoint and stop here
        ! call util_crashstop(830984)

    end subroutine geometry_toplevel_CC
!% 
!%==========================================================================
!%==========================================================================
!%
    subroutine geometry_toplevel_JMJB ()
        !%------------------------------------------------------------------
        !% Description:
        !% Computes geometry on junction JM elements for the  
        !% time-marching scheme 
        !%------------------------------------------------------------------
        !% Declarations
        !%------------------------------------------------------------------

        print *, 'OBSOLETE geometry_toplevel_JMJB'
        stop 298372
        !% --- JB VALUES
        !%    assign the non-volume geometry on junction branches JB based on JM head
        !%    Values limited by full volume. Volume assigned is area * length
        call geo_assign_JB_from_head (ep_JM)

            ! call util_utest_CLprint ('------- in geometry after geo_assign_JB_from_head')
        
        !% --- JB slot adjustments
        !%     make the slot adjustments for JB after the geometry is assigned
        !%     this slot adjustment is based on the head on the JB
        call slot_JB_computation (ep_JM)


        !% --- JM values
       ! call geo_JM_depth_area_from_volume ()
        !% --- new junction plan area
        call geo_plan_area_from_volume_JM (elemPGetm, npack_elemPGetm, col_elemPGetm)
        !% --- ne junction depth 
        call geo_depth_from_volume_JM (elemPGetm, npack_elemPGetm, col_elemPGetm)


        end subroutine geometry_toplevel_JMJB 
!% 
!%==========================================================================
!%==========================================================================    
!%
    subroutine geo_common_initialize &
            (thisP, thisGeoType, ATableType, TTableType, RTableType, STableType)
        !%------------------------------------------------------------------
        !% Description
        !% Performs initializations that are similar for all geometryp
        !% ATableType = lookup table key for area
        !% TTableType = lookup table key for topwidth
        !% AoverAcol  = column of esgr_...AoverAfull
        !% YoverYcol  = column of esgr_...YoverYfull
        !% Values must already be defined for elemR(:,...) for
        !%    Depth, FullDepth, FullArea, FullHydRadius, FullTopwidth, Zbottom
        !%    BreadthMax, DepthAtMaxBreadth, Length
        !% NOTE: do not use setting%ZeroValue%... (except setting%ZeroValue%Depth)
        !%  as they have not been correctly set when this procedure is called
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisP(:), thisGeoType
            real(8), intent(in) :: ATableType(:), TTableType(:), RTableType(:), STableType(:)
            integer             :: ii, mm
            real(8), pointer    :: depth(:), tempDepth(:), fullArea(:)
            real(8), pointer    :: breadthMax(:), fullDepth(:), depthAtMaxBreadth(:)
            real(8), pointer    :: fullTopwidth(:), fullPerimeter(:), fullHydRadius(:)
            real(8) :: topwidthDepth
        !%------------------------------------------------------------------
        !% Aliases
            depth             => elemR(:,er_Depth)
            tempDepth         => elemR(:,er_Temp02)
            fullArea          => elemR(:,er_FullArea)
            fullDepth         => elemR(:,er_FullDepth)
            fullPerimeter     => elemR(:,er_FullPerimeter)
            fullTopwidth      => elemR(:,er_FullTopWidth)
            fullHydRadius     => elemR(:,er_FullHydRadius)
            breadthMax        => elemR(:,er_BreadthMax)
            depthAtMaxBreadth => elemR(:,er_DepthAtBreadthMax)
        !%------------------------------------------------------------------
        !% --- general analytical functions
        !%     Note that sediment depth is added to use this with filled circular
        !%     but sediment depth = 0 should apply for all other cross-sections
        elemR(thisP,er_FullPerimeter)  = elemR(thisP,er_FullArea)    &
                                       / elemR(thisP,er_FullHydRadius)
        
        elemR(thisP,er_Zcrown)         = elemR(thisP,er_Zbottom)      &
                                       + elemR(thisP,er_FullDepth)    &
                                       + elemR(thisP,er_SedimentDepth)

        elemR(thisP,er_FullVolume)     = elemR(thisP,er_FullArea)     &
                                       * elemR(thisP,er_Length)

        elemR(thisP,er_ZbreadthMax)    = elemR(thisP,er_DepthAtBreadthMax) &
                                       + elemR(thisP,er_SedimentDepth)     &
                                       + elemR(thisP,er_Zbottom)           
                                        
        !% --- get full topwidth depth and area below max breadth
        select case (thisGeoType)

            case (arch, basket_handle, catenary, circular, custom, eggshaped, gothic, &
                    horiz_ellipse, horseshoe, semi_circular, semi_elliptical, vert_ellipse)

                !% --- lookup functions using specific depths
                do ii=1,size(thisP)
                    mm = thisP(ii)

                    topwidthDepth = fullDepth(mm) &
                                    * setting%Discretization%FullConduitTopwidthDepthFraction

                    elemR(mm,er_FullTopwidth) = llgeo_tabular_from_depth_singular                  &
                        (mm, topwidthDepth, breadthMax(mm), setting%ZeroValue%Depth, zeroR, TTableType)   

                    elemR(mm,er_AreaBelowBreadthMax) =  llgeo_tabular_from_depth_singular           &
                        (mm, depthAtMaxBreadth(mm), fullArea(mm), setting%ZeroValue%Depth, zeroR, ATableType)   
                end do

            case (filled_circular)
                do ii=1,size(thisP)
                    mm = thisP(ii)

                    !% --- lookup functions using specific depths
                    topwidthDepth = fullDepth(mm) &
                                    * setting%Discretization%FullConduitTopwidthDepthFraction 
                
                    elemR(mm,er_FullTopwidth) = llgeo_filled_circular_topwidth_from_depth_singular     &
                                    (mm, topwidthDepth, zeroR)                                             

                    elemR(mm,er_AreaBelowBreadthMax) = llgeo_filled_circular_area_from_depth_singular  &   
                                    (mm, depthAtMaxBreadth(mm), zeroR)
                end do

            case (mod_basket)
                do ii=1,size(thisP)
                    mm = thisP(ii)

                    !% --- lookup functions using specific depths
                    topwidthDepth = fullDepth(mm) &
                                    * setting%Discretization%FullConduitTopwidthDepthFraction 

                    elemR(mm,er_FullTopwidth) = llgeo_mod_basket_topwidth_from_depth_singular     &
                                                (mm, topwidthDepth, zeroR)                                             

                    elemR(mm,er_AreaBelowBreadthMax) = llgeo_mod_basket_area_from_depth_singular  &   
                                                (mm, depthAtMaxBreadth(mm), zeroR)
                end do

            case (rectangular_closed)
                elemR(thisP,er_FullTopWidth)        = elemR(thisP,er_BreadthMax)
                elemR(thisP,er_AreaBelowBreadthMax) = onehalfR * elemR(thisP,er_FullArea) 

            case (rect_round)
                elemR(thisP,er_FullTopwidth) = elemR(thisP,er_BreadthMax)
                do ii=1,size(thisP)
                    mm = thisP(ii)
                    elemR(thisP(mm),er_AreaBelowBreadthMax) = llgeo_rect_round_area_from_depth_singular &
                                                               (mm, depthAtMaxBreadth(mm), zeroR)
                end do

            case (rect_triang)
                elemR(thisP,er_FullTopwidth) = elemR(thisP,er_BreadthMax)
                do ii=1,size(thisP)
                    mm = thisP(ii)
                    elemR(thisP(mm),er_AreaBelowBreadthMax) = llgeo_rectangular_triangular_area_from_depth_singular &
                                                               (mm, depthAtMaxBreadth(mm), zeroR)
                end do


            case default
                print *, 'CODE ERROR: Unexpected case default'
                call util_crashpoint(5298733)
        end select

        ! print *, ' '
        ! print *, 'in common geo'
        ! print *, 'thisP ',thisP
        ! print *, elemR(thisP,er_FullTopwidth)
        ! print *, ' '

        !% temporary store of depth
        tempDepth(thisP) = depth(thisP)

        !% --- temporary store or computing full values with general functions
        elemR(thisP,er_Depth)     = elemR(thisP,er_FullDepth)
        elemR(thisP,er_Area)      = elemR(thisP,er_FullArea)
        elemR(thisP,er_Perimeter) = elemR(thisP,er_FullPerimeter)
        elemR(thisP,er_Topwidth)  = elemR(thisP,er_FullTopwidth)

       

        !% --- standard functions using temporary store
        !elemR(thisP,er_FullEllDepth)       = llgeo_FullEll_pure(thisP) 
        !elemR(thisP,er_FullHydDepth)  = llgeo_hyddepth_from_area_and_topwidth_pure &
        !                                    (thisP, fullarea(thisP), fulltopwidth(thisP))
        fullperimeter(thisP) = max(fullperimeter(thisP),setting%ZeroValue%Depth)
        elemR(thisP,er_FullHydRadius) = llgeo_hydradius_from_area_and_perimeter_pure &
                                            (thisP, fullarea(thisP), fullperimeter(thisP))
        
        !% --- reset the initial depth
        depth(thisP) = tempDepth(thisP)
        
        !% --- get IC data for Area
        select case (thisGeoType)

            case (arch, basket_handle, circular, eggshaped, horiz_ellipse, horseshoe, vert_ellipse)

                do ii=1,size(thisP)
                    mm = thisP(ii)
                    elemR(mm,er_Area)        = llgeo_tabular_from_depth_singular &
                        (mm, depth(mm), fullArea(mm), setting%ZeroValue%Depth, zeroR, ATableType)

                    elemR(mm,er_Topwidth)    = llgeo_tabular_from_depth_singular &
                        (mm, depth(mm), breadthMax(mm), setting%ZeroValue%Depth, zeroR, TTableType)

                    elemR(mm,er_HydRadius)  = llgeo_tabular_from_depth_singular &
                        (mm, depth(mm), fullHydRadius(mm), setting%ZeroValue%Depth, zeroR, RtableType)

                    if (elemR(mm,er_HydRadius) > zeroR) then 
                        elemR(mm,er_Perimeter) = elemR(mm,er_Area) / elemR(mm,er_HydRadius)
                    else 
                        elemR(mm,er_Perimeter) = zeroR
                    end if
                end do        

            case (catenary, gothic, semi_circular, semi_elliptical)

                do ii=1,size(thisP)
                    mm = thisP(ii)
                    elemR(mm,er_Area)        = llgeo_tabular_from_depth_singular &
                        (mm, depth(mm), fullArea(mm), setting%ZeroValue%Depth, zeroR, ATableType)

                    elemR(mm,er_Topwidth)    = llgeo_tabular_from_depth_singular &
                        (mm, depth(mm), breadthMax(mm), setting%ZeroValue%Depth, zeroR, TTableType)

                    elemR(mm,er_HydRadius)  = llgeo_tabular_from_depth_singular &
                        (mm, depth(mm), fullHydRadius(mm), setting%ZeroValue%Depth, zeroR, RTableType)

                    elemR(mm,er_HydRadius) = llgeo_tabular_hydradius_from_area_and_sectionfactor_singular &
                        (mm, elemR(mm,er_Area) , fullhydradius(mm), zeroR, STableType)

                    if (elemR(mm,er_HydRadius) > zeroR) then 
                        elemR(mm,er_Perimeter) = elemR(mm,er_Area) / elemR(mm,er_HydRadius)
                    else 
                        elemR(mm,er_Perimeter) = zeroR
                    end if
                end do

            case (filled_circular)
                do ii=1,size(thisP)
                    mm = thisP(ii)
                    elemR(mm,er_Area)     = llgeo_filled_circular_area_from_depth_singular &
                                            (mm, depth(mm), zeroR)
                    elemR(mm,er_Topwidth) = llgeo_filled_circular_topwidth_from_depth_singular &
                                            (mm, depth(mm), zeroR)
                    elemR(mm,er_Perimeter) = llgeo_filled_circular_perimeter_from_depth_singular &
                                            (mm, depth(mm), zeroR)   
                    if (elemR(mm,er_Perimeter) > zeroR) then 
                        elemR(mm,er_HydRadius) = elemR(mm,er_Area) / elemR(mm,er_Perimeter)
                    else 
                        elemR(mm,er_HydRadius) = zeroR
                    end if         
                end do

            case (mod_basket)
                do ii = 1,size(thisP)
                    mm = thisP(ii)
                    elemR(mm,er_Area)      = llgeo_mod_basket_area_from_depth_singular &
                                            (mm, depth(mm), zeroR)
                    elemR(mm,er_Topwidth)  = llgeo_mod_basket_topwidth_from_depth_singular &
                                            (mm, depth(mm), zeroR)          
                    elemR(mm,er_Perimeter) = llgeo_mod_basket_perimeter_from_depth_singular &
                                            (mm, depth(mm), zeroR)  
                    if (elemR(mm,er_Perimeter) > zeroR) then 
                        elemR(mm,er_HydRadius) = elemR(mm,er_Area) / elemR(mm,er_Perimeter)
                    else 
                        elemR(mm,er_HydRadius) = zeroR
                    end if            
                end do

            case (rectangular_closed)
                do ii = 1,size(thisP)
                    mm = thisP(ii)
                    elemR(mm,er_Area)      = llgeo_rectangular_closed_area_from_depth_singular &
                                                (mm, depth(mm), zeroR)
                    elemR(mm,er_TopWidth)  = llgeo_rectangular_closed_topwidth_from_depth_singular &
                                                (mm, depth(mm), zeroR)
                    elemR(mm,er_Perimeter) = llgeo_rectangular_closed_perimeter_from_depth_singular &
                                                (mm, depth(mm), zeroR)
                    if (elemR(mm,er_Perimeter) > zeroR) then 
                        elemR(mm,er_HydRadius) = elemR(mm,er_Area) / elemR(mm,er_Perimeter)
                    else 
                        elemR(mm,er_HydRadius) = zeroR
                    end if
                end do
                
            case (rect_round)
                do ii=1,size(thisP)
                    mm = thisP(ii)
                    elemR(mm,er_Area)      = llgeo_rect_round_area_from_depth_singular &
                                            (mm, depth(mm), zeroR)
                    elemR(mm,er_Topwidth)  = llgeo_rect_round_topwidth_from_depth_singular &
                                            (mm, depth(mm), zeroR)
                    elemR(mm,er_Perimeter) = llgeo_rect_round_perimeter_from_depth_singular &
                                            (mm, depth(mm), zeroR)   
                    if (elemR(mm,er_Perimeter) > zeroR) then 
                        elemR(mm,er_HydRadius) = elemR(mm,er_Area) / elemR(mm,er_Perimeter)
                    else 
                        elemR(mm,er_HydRadius) = zeroR
                    end if                  
                end do
                !print*, elemR(thisP,er_Area), 'elemR(thisP,er_Area)'

            case (rect_triang)
                do ii=1,size(thisP)
                    mm = thisP(ii)
                    elemR(mm,er_Area)      = llgeo_rectangular_triangular_area_from_depth_singular &
                                            (mm, depth(mm), zeroR)
                    elemR(mm,er_Topwidth)  = llgeo_rectangular_triangular_topwidth_from_depth_singular &
                                            (mm, depth(mm), zeroR)   
                    elemR(mm,er_Perimeter) = llgeo_rectangular_triangular_perimeter_from_depth_singular &
                                            (mm, depth(mm), zeroR) 
                    if (elemR(mm,er_Perimeter) > zeroR) then 
                        elemR(mm,er_HydRadius) = elemR(mm,er_Area) / elemR(mm,er_Perimeter)
                    else 
                        elemR(mm,er_HydRadius) = zeroR
                    end if                                                                         
                end do
            case default
                print *, 'CODE ERROR: Unexpected case default'
                call util_crashpoint(5298722)
        end select

        elemR(thisP,er_Volume)        = elemR(thisP,er_Area) * elemR(thisP,er_Length)

        !% --- set areas, topwidth, perimeter for D = ZeroValue%Depth to zero
        !%     later these are corrected with setting%ZeroValue%...
        where (elemR(thisP,er_Depth) .le. setting%ZeroValue%Depth)
            elemR(thisP,er_Volume)   = zeroR !setting%ZeroValue%Volume
            elemR(thisP,er_Area)     = zeroR !setting%ZeroValue%Area
            elemR(thisP,er_Topwidth) = zeroR !setting%ZeroValue%Topwidth
            elemR(thisP,er_Perimeter)= zeroR !setting%ZeroValue%Topwidth
            elemR(thisP,er_HydRadius)= zeroR !setting%ZeroValue%Depth
        end where

        elemR(thisP,er_AoverAfull) = elemR(thisP,er_Area)  / elemR(thisP,er_FullArea)
        elemR(thisP,er_YoverYfull) = elemR(thisP,er_Depth) / elemR(thisP,er_FullDepth)

        call slot_initialize (thisP)

        !% store IC data
        elemR(thisP,er_Area_N0)       = elemR(thisP,er_Area)
        elemR(thisP,er_Area_N1)       = elemR(thisP,er_Area)
        elemR(thisP,er_Volume_N0)     = elemR(thisP,er_Volume)
        elemR(thisP,er_Volume_N1)     = elemR(thisP,er_Volume)

        !%--- clear the temp storage
        tempDepth(thisP) = nullvalueR

        !% note that er_Perimeter, er_Topwidth, er_HydRadius,, er_EllDepth are NOT initialized
    end subroutine geo_common_initialize    
!%
!%==========================================================================
!%==========================================================================    
!%
    real(8) function geo_sectionfactor_from_depth_singular &
         (eIdx, inDepth, ZeroValueArea, ZeroValuePerimeter) result (outvalue)  
        !%------------------------------------------------------------------
        !% Description
        !% computes the section factor for element with index eIdx for
        !% the depth "inDepth"
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in)  :: eIdx
            real(8), intent(in)  :: inDepth, ZeroValueArea, ZeroValuePerimeter
            real(8) :: thisPerimeter, thisArea
            character(64) :: subroutine_name = "geo_sectionfactor_from_depth_singular"
        !%------------------------------------------------------------------  
    
         !print *, 'input depth ',inDepth, setting%ZeroValue%Depth
        thisArea      = geo_area_from_depth_singular      (eIdx, inDepth, ZeroValueArea)
         !print *, '----- area     ',thisArea, setting%ZeroValue%Area
        thisPerimeter = geo_perimeter_from_depth_singular (eIdx, inDepth, ZeroValuePerimeter)
         !print *, '----- perimeter',thisPerimeter
        outvalue      = thisArea * ((thisArea / thisPerimeter)**twothirdR)
         !print *, '----- sf       ',outvalue

        !write(*,"(10f12.5)") inDepth, thisArea, thisPerimeter, outvalue 

    end function geo_sectionfactor_from_depth_singular
!%
!%==========================================================================    
!%==========================================================================   
!%
    real(8) function geo_Qcritical_from_depth_singular &
         (eIdx, inDepth, ZeroValue) result (outvalue)
        !%------------------------------------------------------------------
        !% computes the critical flow for element eIdx with depth "inDepth"
        !%------------------------------------------------------------------
         !% Declarations
         integer, intent(in)  :: eIdx
         real(8), intent(in)  :: inDepth, ZeroValue
         real(8), pointer     :: grav
         real(8)              :: thisArea
        !%------------------------------------------------------------------
            grav => setting%Constant%gravity
        !%------------------------------------------------------------------   

        thisArea      = geo_area_from_depth_singular (eIdx, inDepth, ZeroValue)
        outvalue      = thisArea * sqrt(inDepth * grav)

    end function geo_Qcritical_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function geo_critical_value_singular (UT_idx, utd_sheet) result (outvalue)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the critical depth or area for the uniformtable(UT_idx)
        !% using the flowrate in the associated element eIdx
        !% utd_sheet = utd_Qcrit_depth_nonuniform
        !%           = utd_Qcrit_area_nonuniform
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: UT_idx, utd_sheet
            integer, pointer    :: eIdx
            real(8), pointer    :: gravity, thistable(:)
            real(8)             :: QcritNormalized
            integer :: ii, utr_Max
        !%------------------------------------------------------------------
        !% Aliases
            eIdx      => uniformTableI(UT_idx,uti_elem_idx)
            thisTable => uniformTableDataR(UT_idx,:,utd_sheet)
        !%------------------------------------------------------------------
        !% Preliminaries
            select case (utd_sheet)
                case (utd_Qcrit_depth_nonuniform)
                    utr_Max = utr_DepthMax
                case (utd_Qcrit_area_nonuniform)
                    utr_Max = utr_AreaMax
                case default
                    print *, 'CODE ERROR: unexpected case default'
                    call util_crashpoint(6209873)
                end select
        !%------------------------------------------------------------------
        ! print *, ' '
        ! print *, 'in geo_critical_value_singular'
        ! print *, 'UT_idx',UT_idx
        ! print *, 'eIdx  ',eIdx
        ! print *, 'flowrate     ', elemR(eIdx,er_Flowrate)
        !print *, 'utr_Qcritmax ',utr_QcritMax
        !print *, 'table Qcritmax       ', uniformTableR(UT_idx,utr_QcritMax)
       !do ii=1,N_Elem(this_image())
       !   print *, ii, elemR(ii,er_Flowrate), elemR(ii,er_Head)
       ! end do
        ! print *, ' '
        
        !stop 2098374

        ! print *, 'Qcrit Max ',uniformTableR(UT_idx,utr_QcritMax)

        !% --- normalize the critical flowrate
        QcritNormalized = abs(elemR(eIdx,er_Flowrate) / uniformTableR(UT_idx,utr_QcritMax))
        
            !print *, 'QcritNormalized',QcritNormalized

        !% --- lookup the normalized critical depth for this critical flow
        outvalue = xsect_table_lookup_singular (QcritNormalized, thistable)

            !print *, 'normalized critical depth',outvalue

            !print *, 'DepthMax ',uniformTableR(UT_idx,utr_DepthMax)

        !% --- return depth to physical value
        outvalue = outvalue * uniformTableR(UT_idx,utr_Max)

            ! print *, 'critical depth ',outvalue

    end function geo_critical_value_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function geo_normaldepth_singular &
        (UT_idx) result (outvalue)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the normal depth for the location UT_idx in the
        !% uniform table array
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: UT_idx        ! index of element in the sectionfactonI/R arrays
            integer, pointer    :: eIdx
            real(8), pointer    :: thisTable(:)
            real(8)             :: sectionFactor, normSF
            character(64)       :: subroutine_name = 'geo_normaldepth_singular'
        !%------------------------------------------------------------------
        !% Aliases
            eIdx      => uniformTableI(UT_idx,uti_elem_idx)
            thisTable => uniformTableDataR(UT_idx,:,utd_SF_depth_nonuniform)  !% element index
        !%------------------------------------------------------------------
        !% --- section factor based on flowrate for the associated element
        if (elemR(eIdx,er_BottomSlope) > setting%ZeroValue%Slope) then
            sectionFactor = elemR(eIdx,er_Flowrate) * elemR(eIdx,er_ManningsN) / (sqrt(abs(elemR(eIdx,er_BottomSlope))))
        else
            sectionFactor = elemR(eIdx,er_Flowrate) * elemR(eIdx,er_ManningsN) / (sqrt(setting%ZeroValue%Slope))
        end if

        ! print *, 'sectionFactor ',sectionFactor
        ! print *, 'flowrate   ',elemR(eIdx,er_Flowrate)
        ! print *, 'mannings n ',elemR(eIdx,er_ManningsN)
        ! print *, 'slope      ',elemR(eIdx,er_BottomSlope)
        ! print *, 'beta       ',(elemR(eIdx,er_BottomSlope)**0.5) / elemR(eIdx,er_ManningsN)

        !% --- if flow is negative on a positive slope, or flow is positive on a negative slope,
        !%     then the section factor is negative, which implies an infinite normal depth
        if (sectionFactor .le. zeroR) then
            outvalue = setting%Limiter%NormalDepthInfinite
            return
        end if

        ! print *, 'uniform table ', uniformTableR(UT_idx,utr_SFmax)

        ! print *, size(uniformTableR,1), size(uniformTableR,2)

        !% --- normalize the section factor
        normSF   = sectionFactor / uniformTableR(UT_idx,utr_SFmax)

        ! print *, 'normSF ',normSF

        !% --- lookup the normalized normal depth
        outvalue = xsect_table_lookup_singular(normSF,thisTable)

        ! print *, 'normalized normal depth ',outvalue

        !% --- return normal depth to physical value
        outvalue = outvalue * uniformTableR(UT_idx,utr_DepthMax)

        ! print *, 'normal depth ',outvalue
    
    end function geo_normaldepth_singular
!%
!%==========================================================================
!% PRIVATE (except _singular functions)
!%==========================================================================
!%
    subroutine geo_ponding_inflow (thisColP_JM)
        !%------------------------------------------------------------------
        !% Description:
        !% Provides inflow volume from ponding outside of a JM element when
        !% the interior volume is below the maximum. Does NOT handle
        !% inflow into surcharge from ponding (see slot_JM_ETM)
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisColP_JM  !% packed column for JM elements
            integer, pointer    :: Npack, thisP(:)
            real(8), pointer    :: vPond(:), volume(:), vFull(:), vInflow(:)
        !%------------------------------------------------------------------
        !% Preliminaries
            if (.not. setting%SWMMinput%AllowPonding) return
            Npack => npack_elemP(thisColP_JM)
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases:
            thisp       => elemP(1:Npack,thisColP_JM)
            vPond       => elemR(:,er_VolumePonded)
            volume      => elemR(:,er_Volume)
            vFull       => elemR(:,er_FullVolume)
            vInflow     => elemR(:,er_Temp01)
        !%------------------------------------------------------------------
        !% --- volume available in the JM    
        vInflow(thisP) = vFull(thisP) - volume(thisP)
        !% --- inflow is the lesser of available volume and ponded volume
        vInflow(thisP)  = min(vInflow(thisp), vPond(thisP))

        where (vInflow(thisP) > zeroR)
            volume(thisP) = volume(thisP) + vInflow(thisP)
            vPond(thisP)  = vPond(thisp)  - vInflow(thisP)
        endwhere
  
    end subroutine geo_ponding_inflow
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_depth_from_volume_by_type_allCC (elemPGx, npack_elemPGx, col_elemPGx)
        !%------------------------------------------------------------------
        !% Description:
        !% This updates depths in nonsurcharged CC elements in PGx arrays
        !% The elemPGx determines whether this is ALLtm, ETM or AC elements
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: elemPGx(:,:), npack_elemPGx(:), col_elemPGx(:)
            integer, pointer :: Npack, thisCol
            character(64) :: subroutine_name = 'geo_depth_from_volume_by_type'
        !%-------------------------------------------------------------------
        !% cycle through different geometries  

                ! print *, ' '
                ! print *, 'at start of geo_depth_from_volume_by_type'
                ! print *, 54, elemR(54,er_Depth), elemR(54,er_Volume)
                ! print *, ' '
                
        !% --- open channels ------------------------------------------

        !% --- IRREGULAR
        thisCol => col_elemPGx(epg_CC_irregular)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call irregular_depth_from_volume (elemPGx(1:Npack,thisCol))
            call geo_ZeroDepth_from_volume   (elemPGx(1:Npack,thisCol))
        end if    
                
        !% -- POWER FUNCTION
        thisCol => col_elemPGx(epg_CC_power_function)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            print *, 'POWER FUNCTION CROSS-SECTION NOT COMPLETE'
            call util_crashpoint(5559872)
            !call powerfunction_depth_from_volume (elemPGx(1:Npack,thisCol))
            !call geo_ZeroDepth_from_volume  (elemPGx(1:Npack,thisCol))
        end if

        !% -- PARABOLIC
        thisCol => col_elemPGx(epg_CC_parabolic)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call parabolic_depth_from_volume (elemPGx(1:Npack,thisCol))
            call geo_ZeroDepth_from_volume   (elemPGx(1:Npack,thisCol))
        end if
                
        !% --- RECTANGULAR CHANNEL
        thisCol => col_elemPGx(epg_CC_rectangular)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
                ! print *, 'rect here Volume/length',elemR(15,er_Volume)/elemR(15,er_Length)
                ! print *,  'zero area           ',setting%ZeroValue%Area
            call rectangular_depth_from_volume (elemPGx(1:Npack,thisCol))
            call geo_ZeroDepth_from_volume     (elemPGx(1:Npack,thisCol))
        end if    

        !% --- TRAPEZOIDAL CHANNEL
        thisCol => col_elemPGx(epg_CC_trapezoidal)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call trapezoidal_depth_from_volume (elemPGx(1:Npack,thisCol))
            call geo_ZeroDepth_from_volume     (elemPGx(1:Npack,thisCol))
        end if

        !% --- TRIANGULAR CHANNEL
        thisCol => col_elemPGx(epg_CC_triangular)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call triangular_depth_from_volume (elemPGx(1:Npack,thisCol))
            call geo_ZeroDepth_from_volume    (elemPGx(1:Npack,thisCol))
        end if

        !% --- CLOSED CONDUITS  ---------------------------------------

        !% --  ARCH CONDUIT
        thisCol => col_elemPGx(epg_CC_arch)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume   &
                (elemPGx, Npack, thisCol, YArch)
            call geo_ZeroDepth_from_volume  (elemPGx(1:Npack,thisCol))
        end if

        !% --  BASKET_HANDLE
        thisCol => col_elemPGx(epg_CC_basket_handle)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume           &
                (elemPGx, Npack, thisCol, YBasketHandle)
            call geo_ZeroDepth_from_volume  (elemPGx(1:Npack,thisCol))   
        end if

        !% --  CATENARY
        thisCol => col_elemPGx(epg_CC_catenary)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume           &
                (elemPGx, Npack, thisCol, YCatenary)
            call geo_ZeroDepth_from_volume  (elemPGx(1:Npack,thisCol))    
        end if

        !% --- CIRCULAR CONDUIT
        thisCol => col_elemPGx(epg_CC_circular)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call circular_depth_from_volume (elemPGx(1:Npack,thisCol))
            call geo_ZeroDepth_from_volume  (elemPGx(1:Npack,thisCol))
        end if

        !% --  EGG_SHAPED
        thisCol => col_elemPGx(epg_CC_egg_shaped)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume           &
                (elemPGx, Npack, thisCol, YEgg)
            call geo_ZeroDepth_from_volume  (elemPGx(1:Npack,thisCol))
        end if

        !% --- FILLED CIRCULAR CONDUIT
        thisCol => col_elemPGx(epg_CC_filled_circular)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call filled_circular_depth_from_volume (elemPGx(1:Npack,thisCol))
            call geo_ZeroDepth_from_volume         (elemPGx(1:Npack,thisCol))
        end if

        !% --  GOTHIC
        thisCol => col_elemPGx(epg_CC_gothic)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume           &
                (elemPGx, Npack, thisCol, YGothic)
            call geo_ZeroDepth_from_volume  (elemPGx(1:Npack,thisCol))    
        end if

        !% --  HORIZONTAL ELLIPSE
        thisCol => col_elemPGx(epg_CC_horiz_ellipse)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume           &
                (elemPGx, Npack, thisCol, YHorizEllip)
            call geo_ZeroDepth_from_volume  (elemPGx(1:Npack,thisCol))    
        end if

        !% --  HORSESHOE
        thisCol => col_elemPGx(epg_CC_horse_shoe)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume           &
                (elemPGx, Npack, thisCol, YHorseShoe)
            call geo_ZeroDepth_from_volume  (elemPGx(1:Npack,thisCol))    
        end if

        !% --  MODIFIED BASKET HANDLE
        thisCol => col_elemPGx(epg_CC_mod_basket)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call mod_basket_depth_from_volume (elemPGx(1:Npack,thisCol))
            call geo_ZeroDepth_from_volume    (elemPGx(1:Npack,thisCol))
        end if

        !% --- RECTANGULAR CLOSED CONDUIT
        thisCol => col_elemPGx(epg_CC_rectangular_closed)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call rectangular_closed_depth_from_volume (elemPGx(1:Npack,thisCol))
            call geo_ZeroDepth_from_volume            (elemPGx(1:Npack,thisCol))
        end if        

        !% --- RECTANGULAR ROUND CONDUIT
        thisCol => col_elemPGx(epg_CC_rectangular_round)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call rect_round_depth_from_volume (elemPGx(1:Npack,thisCol))
            call geo_ZeroDepth_from_volume    (elemPGx(1:Npack,thisCol))
        end if

        !% -- RECTANGULAR TRIANGULAR
        thisCol => col_elemPGx(epg_CC_rectangular_triangular)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call rectangular_triangular_depth_from_volume (elemPGx(1:Npack,thisCol))
            call geo_ZeroDepth_from_volume                (elemPGx(1:Npack,thisCol))
        end if

        ! % --- SEMI-CIRCULAR CONDUIT
        thisCol => col_elemPGx(epg_CC_semi_circular)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume           &
                (elemPGx, Npack, thisCol, YSemiCircular)
            call geo_ZeroDepth_from_volume  (elemPGx(1:Npack,thisCol))    
        end if

        !% --  SEMI ELLIPTICAL
        thisCol => col_elemPGx(epg_CC_semi_elliptical)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume           &
                (elemPGx, Npack, thisCol, YSemiEllip)
            call geo_ZeroDepth_from_volume  (elemPGx(1:Npack,thisCol))    
        end if

        !% --  VERTICAL ELLIPSE
        thisCol => col_elemPGx(epg_CC_vert_ellipse)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume           &
                (elemPGx, Npack, thisCol, YVertEllip)
            call geo_ZeroDepth_from_volume  (elemPGx(1:Npack,thisCol))    
        end if

    end subroutine geo_depth_from_volume_by_type_allCC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_depth_from_volume_by_element_CC (thisP, Npack)
        !%------------------------------------------------------------------
        !% Description
        !% Companion to geo_depth_from_volume_by_type_allCC that
        !% computes for the entire set of cells of a given type.  This
        !% cycles through the set, so is less efficient, but is needed
        !% where only a subset of elements to be evaluated
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:), Npack
            integer :: mm
            integer, dimension(1) :: ap
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        do mm=1,Npack
            ap(1) = thisP(mm)
            select case (elemI(ap(1),ei_geometryType))
                case (irregular)
                    call irregular_depth_from_volume (ap)
                case (power_function)
                    print *, 'POWER FUNCTION CROSS-SECTION NOT COMPLETE'
                    call util_crashpoint(5559867)
                case (parabolic)
                    call parabolic_depth_from_volume (ap)
                case (rectangular)
                    call rectangular_depth_from_volume (ap)
                case (trapezoidal)
                    call trapezoidal_depth_from_volume (ap)
                case (triangular)
                    call triangular_depth_from_volume (ap)
                case (arch)
                    elemR(ap,er_Depth) &
                        = llgeo_tabular_depth_from_volume_singular (ap(1), YArch)
                case (basket_handle)
                    elemR(ap,er_Depth) &
                        = llgeo_tabular_depth_from_volume_singular (ap(1), YBasketHandle)
                case (catenary)
                    elemR(ap,er_Depth) &
                        = llgeo_tabular_depth_from_volume_singular (ap(1), YCatenary)
                case (circular)
                    call circular_depth_from_volume (ap)
                case (eggshaped)
                    elemR(ap,er_Depth) &
                        = llgeo_tabular_depth_from_volume_singular (ap(1), YEgg)
                case (filled_circular)
                    call filled_circular_depth_from_volume (ap)
                case (gothic)
                    elemR(ap,er_Depth) &
                        = llgeo_tabular_depth_from_volume_singular (ap(1), YGothic)
                case (horiz_ellipse)
                    elemR(ap,er_Depth) &
                        = llgeo_tabular_depth_from_volume_singular (ap(1), YHorizEllip)
                case (horseshoe)
                    elemR(ap,er_Depth) &
                        = llgeo_tabular_depth_from_volume_singular (ap(1), YHorseShoe)
                case (mod_basket)
                    call mod_basket_depth_from_volume (ap)
                case (rectangular_closed)
                    call rectangular_closed_depth_from_volume (ap)
                case (rect_round)
                    call rect_round_depth_from_volume (ap)
                case (rect_triang)
                    call rectangular_triangular_depth_from_volume (ap)
                case (semi_circular)
                    elemR(ap,er_Depth) &
                        = llgeo_tabular_depth_from_volume_singular (ap(1), YSemiCircular)
                case (semi_elliptical)
                    elemR(ap,er_Depth) &
                    = llgeo_tabular_depth_from_volume_singular (ap(1), YSemiEllip)
                case (vert_ellipse)
                    elemR(ap,er_Depth) &
                        = llgeo_tabular_depth_from_volume_singular (ap(1), YVertEllip)
                case (custom)
                    print *, 'CUSTOM CROSS-SECTION NOT COMPLETE'
                    call util_crashpoint(52498767)
                case default
                    print *, 'CODE ERROR: Unexpected case default'
                    call util_crashpoint(7209874)
            end select

            !% --- set zero depths
            call geo_ZeroDepth_from_volume (ap)

        end do

    end subroutine geo_depth_from_volume_by_element_CC    
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_depth_from_volume_JM (elemPGx, npack_elemPGx, col_elemPGx)
        !%------------------------------------------------------------------
        !% Description:
        !% This solves nonsurcharged CCJMJB elements because of PGx arrays
        !% The elemPGx determines whether this is ALLtm, ETM or AC elements
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: elemPGx(:,:), npack_elemPGx(:), col_elemPGx(:)
            integer, pointer :: Npack, thisP(:), thisCol
            character(64) :: subroutine_name = 'geo_depth_from_volume_by_type_JM'
        !%-------------------------------------------------------------------
        
        !% --- JUNCTIONS ---------------------------------------------------- 

        !% JM with functional geometry
        thisCol => col_elemPGx(epg_JM_functionalStorage)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            thisP => elemPGx(1:Npack,thisCol)
            call storage_functional_depth_from_volume (thisP,Npack)
            call geo_ZeroDepth_from_volume  (thisP)
        end if

        !% JM with tabular geometry
        thisCol => col_elemPGx(epg_JM_tabularStorage)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            thisP => elemPGx(1:Npack,thisCol)
            call storage_tabular_depth_from_volume (thisP, Npack)
            call geo_ZeroDepth_from_volume  (thisP)
        end if

        !% JM with implied storage 
        thisCol => col_elemPGx(epg_JM_impliedStorage)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            !print *, 'in implied storage'
            thisP => elemPGx(1:Npack,thisCol)
            call storage_implied_depth_from_volume (thisP, Npack)
            call geo_ZeroDepth_from_volume  (thisP)
        end if

        !% --- note that NoStorage junctions have no volume

    end subroutine geo_depth_from_volume_JM
 !%
!%==========================================================================

!%==========================================================================
!%
    subroutine geo_plan_area_from_volume_JM (elemPGx, npack_elemPGx, col_elemPGx)
        !%------------------------------------------------------------------
        !% Description:
        !% This calculates plan area for nonsurcharged JM elements because of PGx arrays
        !% The elemPGx determines whether this is ALLtm, ETM or AC elements
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: elemPGx(:,:), npack_elemPGx(:), col_elemPGx(:)
            integer, pointer :: Npack, thisP(:), thisCol
            character(64) :: subroutine_name = 'geo_depth_from_volume_by_type_JM'
        !%-------------------------------------------------------------------
        
        !% --- JUNCTIONS ---------------------------------------------------- 

        !% --- JM with functional geometry
        thisCol => col_elemPGx(epg_JM_functionalStorage)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            thisP => elemPGx(1:Npack,thisCol)
            call storage_plan_area_from_volume (thisP,Npack)
        end if

        !% --- JM with tabular geometry
        thisCol => col_elemPGx(epg_JM_tabularStorage)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            thisP => elemPGx(1:Npack,thisCol)
            call storage_plan_area_from_volume (thisP,Npack)
        end if

        !% --- JM with implied geometry
        !%     no action: plan area is fixed

    end subroutine geo_plan_area_from_volume_JM    
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_overflow_openchannels (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Adds the overflow of open channels that exceed their 
        !% Full Volume to the overflow accumulator for this step
        !% Reduces stored volume by the overflow amount.
        !% Note that open channels CANNOT pond in present version.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisP(:) !% must be open channels
            integer, pointer    :: nBarrels(:)
            real(8), pointer    :: volume(:), fullvolume(:), overflow(:)
            real(8), pointer    :: temp(:)
        !%------------------------------------------------------------------
        !% Preliminaries
            if (.not. setting%Discretization%AllowChannelOverflowTF) return
        !%------------------------------------------------------------------
        !% Aliases
            volume     => elemR(:,er_Volume)
            fullvolume => elemR(:,er_FullVolume)
            overflow   => elemR(:,er_VolumeOverFlow) 
            temp       => elemR(:,er_Temp01)
            nBarrels   => elemI(:,ei_barrels)
        !%------------------------------------------------------------------

        !% --- compute the potential overflow
        temp(thisP) = volume(thisP) - fullvolume(thisP)
        where (temp(thisP) > zeroR)
            !% --- overflow is cumulative within the time step
            !%     accounts for multiple barrels
            overflow(thisP) = overflow(thisP) * real(nBarrels(thisP),8) + temp(thisP) 
            !% --- set volume to full
            volume(thisP)   = fullvolume(thisP)
        endwhere
        !% --- reset the temp space
        temp(thisP) = zeroR
    
    end subroutine geo_overflow_openchannels
!%
!%==========================================================================    
!%==========================================================================
!%
    subroutine geo_volumelimit_closed (thisP) 
        !%------------------------------------------------------------------
        !% Description
        !% Sets closed element volumes to the full volume
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisP(:)
            real(8), pointer    :: volume(:), fullvolume(:)
        !%------------------------------------------------------------------
        !% Aliases:
            volume     => elemR(:,er_Volume)
            fullvolume => elemR(:,er_FullVolume)
        !%------------------------------------------------------------------

        !% limit the full volume by the full volume
        volume(thisP) = min( volume(thisP), fullvolume(thisP) )

    end subroutine geo_volumelimit_closed
!%   
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_assign_JB_from_head (thisColP_JM)
        !%------------------------------------------------------------------
        !% Description:
        !% Assigns geometry for head, depth, area and volume for JB (junction branches)
        !% When the main head is higher than the branch, this applies
        !% the main head to the branch with an adjustment for head loss.
        !% When the main head is below the branch, this sets the
        !% branch head to the bottom elevation plus a depth implied
        !%------------------------------------------------------------------
        !% by a Froude number of one.
        !%
        !% Note that the JB works in an inverse form from the other geometry computations.
        !% That is, for CC, JM we have volume a priori and then compute area, depth etc.
        !% However, for JB we get head then depth diagnostically and must compute area,
        !% etc. before we can get volume.
        !%
        !% 20210611 -- this is written in a simple loop form. See notes in draft SWMM5+
        !% NewCode Framework document on possible changes for a packed vector form.
        !% It is not clear that the number of junctions would make the change useful.
        !%
        !% HACK
        !% The following are NOT assigned on JB
        !% FullHydDepth, FullPerimeter, FullVolume, Roughness
        !%
        !%-------------------------------------------------------------------
            integer, intent(in) ::thisColP_JM

            integer, pointer :: Npack, thisP(:), BranchExists(:), thisSolve(:),  tM
            integer, pointer :: fup(:), fdn(:)
            real(8), pointer :: area(:), depth(:), head(:), hydradius(:)
            real(8), pointer :: length(:), perimeter(:), topwidth(:), velocity(:), flowrate(:)
            real(8), pointer :: volume(:), zBtm(:), Kfac(:), dHdA(:), ellDepth(:) !, ellMax(:)
            real(8), pointer :: zCrown(:), fullArea(:), fulldepth(:), fullperimeter(:)
            real(8), pointer :: sedimentDepth(:), thisTable(:,:)
            real(8), pointer :: fulltopwidth(:), breadthmax(:)
            real(8), pointer :: slotDepth(:), slotVolume(:), overflow(:), fullhydradius(:)
            real(8), pointer :: Atable(:), Ttable(:), Rtable(:), Stable(:)
            real(8), pointer :: grav, fVel_u(:), fVel_d(:) 
            real(8), pointer :: fZcrown_u(:), fZcrown_d(:), fHead_u(:), fHead_d(:) 
            logical, pointer :: isSlot(:)     

            real(8) :: depthnorm, zeroHydRadius
            integer :: tB, ii, kk, tBA(1)

            logical :: isUpBranch

            integer :: printJB = 5

        !% thisColP_JM is the column for the junction mains of a particular
        !% whichTM. For ALL ep_JM, for ETM, ep_JM_ETM, for AC ep_JM_AC
            integer, allocatable :: tempP(:)
            character(64) :: subroutine_name = 'geo_assign_JB_from_head'
            character(256) :: chartemp
        !%---------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%geometry) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            if (setting%Profile%useYN) call util_profiler_start (pfc_geo_assign_JB_from_head)
        !%----------------------------------------------------------------------
        !% Aliases
            Npack         => npack_elemP(thisColP_JM)
            area          => elemR(:,er_Area)
            breadthmax    => elemR(:,er_BreadthMax)
            depth         => elemR(:,er_Depth)
            dHdA          => elemR(:,er_dHdA)
            ellDepth      => elemR(:,er_EllDepth)
            flowrate      => elemR(:,er_Flowrate)
            head          => elemR(:,er_Head)
            hydradius     => elemR(:,er_HydRadius)
            length        => elemR(:,er_Length)
            perimeter     => elemR(:,er_Perimeter)
            sedimentDepth => elemR(:,er_SedimentDepth)
            topwidth      => elemR(:,er_Topwidth)
            velocity      => elemR(:,er_Velocity)
            volume        => elemR(:,er_Volume)
            zBtm          => elemR(:,er_Zbottom)
            zCrown        => elemR(:,er_Zcrown)
            fullArea      => elemR(:,er_FullArea)
            fulldepth     => elemR(:,er_FullDepth)
            fullTopWidth  => elemR(:,er_FullTopWidth)
            fullhydradius => elemR(:,er_FullHydRadius)
            fullperimeter => elemR(:,er_FullPerimeter)
            Kfac          => elemSR(:,esr_JunctionBranch_Kfactor)
            BranchExists  => elemSI(:,esi_JunctionBranch_Exists)
            thisSolve     => elemI(:,ei_tmType)

            fVel_d        => faceR(:,fr_velocity_d)
            fVel_u        => faceR(:,fr_velocity_u)
            fZcrown_d     => faceR(:,fr_Zcrown_d)
            fZcrown_u     => faceR(:,fr_Zcrown_u)
            fHead_d       => faceR(:,fr_Head_d)
            fHead_u       => faceR(:,fr_Head_u)

            fup           => elemI(:,ei_Mface_uL)
            fdn           => elemI(:,ei_Mface_dL)
            grav => setting%Constant%gravity
        !%------------------------------------------------------------------

        if (Npack > 0) then
            thisP  => elemP(1:Npack,thisColP_JM)

            !% cycle through the all the main junctions and each of its branches
            do ii=1,Npack             
                tM => thisP(ii) !% junction main ID

                ! print *, ' '
                ! print *, '===================================================='
                ! print *, 'tM ',tM
              
                !% cycle through the possible junction branches
                do kk=1,max_branch_per_node

                    if (mod(kk,2)==0) then 
                        isUpBranch = .false.
                    else
                        isUpBranch = .true.
                    endif
                    
                    tB = tM + kk !% junction branch ID
                    tBA(1) = tB  !% array for pure array functions

                    ! print *, 'tB ',tB 

                        ! if (tB == printJB) print *, 'HEAD AAA',head(tB)

                    if (BranchExists(tB) .ne. oneI) cycle

                    !% only when a branch exists.
                    !% --- head(tB) has already been updated in junction_calculation
                    !%     for any element that has head(tM) > bottom + sediment
                    !%     here we need to handle any JB elements that would have
                    !%     Froude =1 overflow
                    if ( head(tM) > (zBtm(tB) + sedimentDepth(tB)) ) then
                        !% 20230401 -- no action needed as this is done in junction_calculation

                        !% OBSOLETE IN THIS SECTION
                        ! % for main head above branch bottom entrance use a head
                        ! % loss approach. The branchsign and velocity sign ensure
                        ! % the headloss is added to an inflow and subtracted at
                        ! % an outflow
                        ! % Note this is a time-lagged velocity as the JB velocity
                        ! % is not updated until after face interpolation                                
                        ! head(tB) = head(tM)                                    &
                        !     + branchsign(kk) * sign(oneR,velocity(tB))         &
                        !     * (Kfac(tB) / (twoR * grav)) * (velocity(tB)**twoR)
                        
                        ! print *, 'HEAD A', head(tB)
                        ! print *,'Velocity ',velocity(tB)

                        head(tB) = head(tM)

                        ! if (tB == printJB) print *, 'HEAD BBB',head(tB)

                    else
                        !% for main head below the branch bottom entrance we assign a
                        !% Froude number of one on an inflow to the junction main. Note
                        !% an outflow from a junction main for this case gets head
                        !% of z_bottom of the branch (zero depth).
                        !% Note this is a time-lagged velocity as the JB velocity
                        !% is not updated until after face interpolation
                        ! head(tB) = zBtm(tB) + sedimentDepth(tB)                            &
                        !     + onehalfR * (oneR + branchsign(kk) * sign(oneR,velocity(tB))) &
                        !     *(velocity(tB)**twoR) / (grav) 

                        !% 20230322 brh using velocity on face  
                        if     ((      isUpBranch) .and. (fVel_d(fup(tB)) > zeroR)) then 
                            !% --- upstream branch with inflow

                            ! if (tB == printJB) then 
                            !     print *, ' '
                            !     print *, 'head up  ',fHead_d(fup(tB))
                            !     print *, 'zcrown_d ',fZcrown_d(fup(tB))
                            !     print *, 'Zbottom  ',zBtm(tB), sedimentDepth(tB)
                            !     print *, 'Kfac     ',Kfac(tB)
                            !     print *, 'Vel      ',fVel_d(fup(tB))
                            !     print *, ' '
                            ! end if

                            if (fHead_d(fup(tB)) > fZcrown_d(fup(tB))) then
                                !% --- surcharged inflow of upstream branch
                                !%     reduce by Kfactor
                                head(tB) = fHead_d(fup(tB)) &
                                    - Kfac(tB) * (fVel_d(fup(tB))**twoI) / (twoR * grav)
                            else
                                !% --- Fr=1 inflow from upstream branch
                                head(tB) = zBtm(tB) + sedimentDepth(tB)  &
                                    + (fVel_d(fup(tB))**twoI) / grav
                            endif

                            ! if (tB == printJB) print *, 'HEAD CCC',head(tB)

                        elseif ((.not. isUpbranch) .and. (fVel_u(fdn(tB)) < zeroR)) then
                            !% --- downstream branch with inflow
                            if (fHead_u(fdn(tB)) > fZcrown_u(fdn(tB))) then
                                !% --- surcharged inflow of downstream branch
                                !%     reduce by Kfactor
                                head(tB) = fHead_u(fdn(tB)) &
                                    - Kfac(tB) * (fVel_u(fdn(tB))**twoI) / (twoR * grav)
                            else
                                !% --- Fr=1 inflow from downstream branch
                                head(tB) = zBtm(tB) + sedimentDepth(tB)  &
                                    + (fVel_u(fdn(tB))**twoI) / grav
                            end if

                                ! if (tB == printJB) print *, 'HEAD DDD',head(tB)

                        else 
                            !% --- no flow, set below zerovalue
                            head(tB) = zBtm(tB) + sedimentDepth(tB) + 0.99d0 * setting%ZeroValue%Depth

                                ! if (tB == printJB) print *, 'HEAD EEE',head(tB)

                        end if

                            ! print *, 'HEAD B, velocity', head(tB), velocity(tB)
                            ! print *,'Velocity ',velocity(tB)
                    end if

                        ! if (tB == printJB) print *, 'HEAD ZZZ',head(tB)

                    !% compute provisional depth
                    depth(tB) = head(tB) - (zBtm(tB) + sedimentDepth(tB))

                        ! if (tB == printJB) print *, 'DEPTH ',depth(tB)

                    if (depth(tB) .ge. fulldepth(tB)) then
                        !% surcharged or incipient surcharged
                        depth(tB)        = fulldepth(tB)
                        area(tB)         = fullArea(tB)
                        perimeter(tB)    = fullperimeter(tB)
                        topwidth(tB)     = setting%ZeroValue%Topwidth
                        hydRadius(tB)    = fulldepth(tB) / fullperimeter(tB)
                        dHdA(tB)         = oneR / setting%ZeroValue%Topwidth
                        ellDepth(tBA)    = llgeo_elldepth_pure(tBA)

                    elseif (depth(tB) .le. setting%ZeroValue%Depth) then
                        !% negligible depth is treated with ZeroValues
                        depth(tB)        = setting%ZeroValue%Depth * 0.99d0
                        area(tB)         = setting%ZeroValue%Area
                        topwidth(tB)     = setting%ZeroValue%Topwidth
                        perimeter(tB)    = setting%ZeroValue%Topwidth + setting%ZeroValue%Depth
                        hydRadius(tB)    = setting%ZeroValue%Depth
                        dHdA(tB)         = oneR / setting%ZeroValue%Topwidth
                        ellDepth(tB)     = setting%ZeroValue%Depth * 0.99d0 

                    else
                        !% --- set lookup table names
                        select case (elemI(tB,ei_geometryType))
                            !% --- for tables using HydRadius
                            case (arch)
                                Atable => AArch
                                Rtable => RArch
                                Ttable => TArch
                            case (basket_handle)
                                Atable => ABasketHandle
                                Rtable => RBasketHandle
                                Ttable => TBasketHandle
                            case (circular)
                                Atable => ACirc
                                Rtable => RCirc
                                Ttable => TCirc
                            case (eggshaped)
                                Atable => AEgg
                                Rtable => REgg
                                Ttable => TEgg
                            case (horiz_ellipse)
                                Atable => AHorizEllip
                                Rtable => RHorizEllip
                                Ttable => THorizEllip
                            case (horseshoe)
                                Atable => AHorseShoe
                                Rtable => RHorseShoe
                                Ttable => THorseShoe
                            case (vert_ellipse)
                                Atable => AHorseShoe
                                Rtable => RHorseShoe
                                Ttable => THorseShoe

                            !% --- for tables using section factor
                            case (catenary)
                                Atable => ACatenary
                                Stable => SCatenary
                                Ttable => TCatenary
                            case (gothic)
                                Atable => AGothic
                                Stable => SGothic
                                Ttable => TGothic
                            case (semi_circular)
                                Atable => ASemiCircular
                                Stable => SSemiCircular
                                Ttable => TSemiCircular
                            case (semi_elliptical)
                                Atable => ASemiEllip
                                Stable => SSemiEllip
                                Ttable => TSemiEllip
                        end select

                        !% HACK --- dHdA, which is needed for AC method
                        !% is not computed in the following

                        !% --- not surcharged and non-negligible depth
                        select case (elemI(tB,ei_geometryType))

                            !% --- OPEN CHANNEL
                            !%     typical open channels have computations for area, topwidth, perimeter
                            !%     with standard geo_ function for hydraulic radius, hydraulic depth and ell

                            case (parabolic)   !% analytical
                                ! if (ii > 98) util_utest_CLprint('IIIIaa')
                                area(tBA)     = llgeo_parabolic_area_from_depth_pure         (tBA, depth(tBA))
                                area(tB )     = max(area(tB),setting%ZeroValue%Area)

                                topwidth(tBA) = llgeo_parabolic_topwidth_from_depth_pure     (tBA, depth(tBA))
                                topwidth(tB)  = max(topwidth(tB),setting%ZeroValue%Topwidth)

                                perimeter(tBA)= llgeo_parabolic_perimeter_from_depth_pure    (tBA, depth(tBA))
                                perimeter(tB) = max(perimeter(tB),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

                                hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                    (tBA, area(tBA), perimeter(tBA))
                                hydRadius(tB) = max(hydRadius(tB),setting%ZeroValue%Depth)

                                ellDepth(tBA) = llgeo_hyddepth_from_area_and_topwidth_pure   &
                                                    (tBA, area(tBA), topwidth(tBA))
                                ellDepth(tB)  = max(ellDepth(tB),setting%ZeroValue%Depth*0.99d0)

                            case (power_function) !% POSSIBLY LOOKUP
                                print *, 'CODE ERROR power function x-section not finished'
                                call util_crashpoint(6298349)

                            case (rectangular) !% analytical
                                ! if (ii > 98) util_utest_CLprint('IIIIb')
                                area(tBA)     = llgeo_rectangular_area_from_depth_pure       (tBA, depth(tBA))
                                area(tB)      = max(area(tB),setting%ZeroValue%Area)

                                topwidth(tBA) = llgeo_rectangular_topwidth_from_depth_pure   (tBA, depth(tBA))
                                topwidth(tB)  = max(topwidth(tB),setting%ZeroValue%Topwidth)

                                perimeter(tBA)= llgeo_rectangular_perimeter_from_depth_pure  (tBA, depth(tBA))
                                perimeter(tB) = max(perimeter(tB),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

                                hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                    (tBA, area(tBA), perimeter(tBA))
                                hydRadius(tB) = max(hydRadius(tB),setting%ZeroValue%Depth)

                                ellDepth(tBA) = llgeo_hyddepth_from_area_and_topwidth_pure  &
                                                    (tBA, area(tBA), topwidth(tBA))

                                        !print *, 'ell depth ', ellDepth(tB)

                                ellDepth(tB)  = max(ellDepth(tB),setting%ZeroValue%Depth*0.99d0)

                            case (trapezoidal) !% analytical
                                ! if (ii > 98) util_utest_CLprint('IIIIc')
                                area(tBA)     = llgeo_trapezoidal_area_from_depth_pure       (tBA, depth(tBA))
                                area(tB)      = max(area(tB),setting%ZeroValue%Area)

                                topwidth(tBA) = llgeo_trapezoidal_topwidth_from_depth_pure   (tBA, depth(tBA))
                                topwidth(tB)  = max(topwidth(tB),setting%ZeroValue%Topwidth)

                                perimeter(tBA)= llgeo_trapezoidal_perimeter_from_depth_pure  (tBA, depth(tBA))
                                perimeter(tB) = max(perimeter(tB),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

                                hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                    (tBA, area(tBA), perimeter(tBA))
                                hydRadius(tB) = max(hydRadius(tB),setting%ZeroValue%Depth)
                                
                                ellDepth(tBA)  = llgeo_hyddepth_from_area_and_topwidth_pure   &
                                                    (tBA, area(tBA), topwidth(tBA))
                                ellDepth(tB)  = max(ellDepth(tB),setting%ZeroValue%Depth*0.99d0)

                            case (triangular) !% analytical
                                ! if (ii > 98) util_utest_CLprint('IIIId')
                                area(tBA)     = llgeo_triangular_area_from_depth_pure        (tBA, depth(tBA))
                                area(tB)      = max(area(tB),setting%ZeroValue%Area)

                                topwidth(tBA) = llgeo_triangular_topwidth_from_depth_pure    (tBA, depth(tBA))
                                topwidth(tB)  = max(topwidth(tB),setting%ZeroValue%Topwidth)

                                perimeter(tBA)= llgeo_triangular_perimeter_from_depth_pure   (tBA, depth(tBA))
                                perimeter(tB) = max(perimeter(tB),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

                                hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                    (tBA, area(tBA), perimeter(tBA))
                                hydRadius(tB) = max(hydRadius(tB),setting%ZeroValue%Depth)

                                ellDepth(tBA) = llgeo_hyddepth_from_area_and_topwidth_pure   &
                                                    (tBA, area(tBA), topwidth(tBA))
                                ellDepth(tB)  = max(ellDepth(tB),setting%ZeroValue%Depth*0.99d0)

                            case (irregular)  !% lookup
                                area(tB)     = irregular_geometry_from_depth_singular ( &
                                                    tB,tt_area, depth(tB), elemR(tB,er_FullArea), setting%ZeroValue%Area)

                                ! if (ii > 98) util_utest_CLprint('IIIIe-1')
                                topwidth(tB) = irregular_geometry_from_depth_singular ( &
                                                    tB,tt_width, depth(tB), elemR(tB,er_FullTopWidth),setting%ZeroValue%TopWidth)

                                ! if (ii > 98) util_utest_CLprint('IIIIe-2')
                                !% --- note the irregular stores hyd radius rather than perimeter
                                hydRadius(tB) = irregular_geometry_from_depth_singular ( &
                                                    tB,tt_hydradius, depth(tB), elemR(tB,er_FullHydRadius), setting%ZeroValue%Depth)    

                                !% --- perimeter is derived geometry for irregular
                                ! if (ii > 98) util_utest_CLprint('IIIIe-4')
                                perimeter(tB) = area(tB) / hydRadius(tB)
                                perimeter(tB)= max(perimeter(tB),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

                                !% --- irregular must be continuously-increasing in width
                                ! if (ii > 98) util_utest_CLprint('IIIIe-5')
                                ellDepth(tB)  = geo_hyddepth_from_area_and_topwidth_singular (tB, area(tB), topwidth(tB), setting%ZeroValue%Depth*0.99d0) 

                            !% --- CLOSED CONDUITS
                            !%     closed conduits typically have look-up functions for area, topwidth and hydraulic
                            !%     radius, with standard geo_functions for perimeter, hydraulic depth, and ell
                            !%     However, where analytical functions are used, the perimeter is usually computed
                            !%     first and hydraulic radius is a geo_ function

                            !% --- lookups with Hydraulic Radius stored
                            case (arch, basket_handle, circular, eggshaped, horiz_ellipse, horseshoe, vert_ellipse) 
                                ! if (ii > 98) util_utest_CLprint('IIIIf')                                   
                                area(tB)     = llgeo_tabular_from_depth_singular &
                                    (tB, depth(tB), fullArea(tB), setting%ZeroValue%Depth, setting%ZeroValue%Area, Atable)

                                topwidth(tB) = llgeo_tabular_from_depth_singular &
                                    (tB, depth(tB), breadthmax(tB), setting%ZeroValue%Depth, setting%ZeroValue%Topwidth, Ttable)

                                topwidth(tB) = max(topwidth(tB), fulltopwidth(tB))

                                hydRadius(tB)= llgeo_tabular_from_depth_singular &
                                    (tB, depth(tB), fullHydRadius(tB), setting%ZeroValue%Depth, setting%ZeroValue%Depth, Rtable)

                                perimeter(tBA)= llgeo_perimeter_from_hydradius_and_area_pure &
                                                    (tBA, hydradius(tBA), area(tBA))
                                perimeter(tB) = max(perimeter(tB),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth) 
                                                
                                ellDepth(tBA) = llgeo_elldepth_pure (tBA) 
                                ellDepth(tB)  = max(ellDepth(tB),setting%ZeroValue%Depth*0.99d0)

                            !% --- lookups with SectionFactor stored
                            case (catenary, gothic, semi_circular, semi_elliptical)   
                                ! if (ii > 98) util_utest_CLprint('IIIIg') 
                                area(tB)     = llgeo_tabular_from_depth_singular &
                                    (tB, depth(tB), fullArea(tB), setting%ZeroValue%Depth, setting%ZeroValue%Area, Atable)

                                topwidth(tB) = llgeo_tabular_from_depth_singular &
                                    (tB, depth(tB), breadthmax(tB), setting%ZeroValue%Depth, setting%ZeroValue%Topwidth, Ttable)

                                topwidth(tB) = max(topwidth(tB), fulltopwidth(tB))

                                hydRadius(tB)= llgeo_tabular_hydradius_from_area_and_sectionfactor_singular &
                                    (tB, area(tB), fullhydradius(tB), setting%ZeroValue%Depth, Stable)

                                perimeter(tBA)= llgeo_perimeter_from_hydradius_and_area_pure &
                                                    (tBA, hydradius(tBA), area(tBA))
                                perimeter(tB) = max(perimeter(tB),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

                                ellDepth(tBA) = llgeo_elldepth_pure (tBA) 
                                ellDepth(tB)  = max(ellDepth(tB),setting%ZeroValue%Depth*0.99d0)

                            !% --- lookup with sediment
                            case (filled_circular)  
                                ! if (ii > 98) util_utest_CLprint('IIIIh')
                                area(tB)      = llgeo_filled_circular_area_from_depth_singular      (tB,depth(tB),setting%ZeroValue%Area)
                                topwidth(tB)  = llgeo_filled_circular_topwidth_from_depth_singular  (tB,depth(tB),setting%ZeroValue%Topwidth)
                                perimeter(tB) = llgeo_filled_circular_perimeter_from_depth_singular (tB,depth(tB),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                                
                                hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                    (tBA, area(tBA), perimeter(tBA))
                                hydRadius(tB) = max(hydRadius(tB),setting%ZeroValue%Depth)

                                ellDepth(tBA) = llgeo_elldepth_pure (tBA) 
                                ellDepth(tB)  = max(ellDepth(tB),setting%ZeroValue%Depth*0.99d0)

                            !% --- analytical closed-conduit cases
                            case (mod_basket)   !% analytical            
                                ! if (ii > 98) util_utest_CLprint('IIIIi')                     
                                area(tB)      = llgeo_mod_basket_area_from_depth_singular        (tB,depth(tB),setting%ZeroValue%Area)
                                topwidth(tB)  = llgeo_mod_basket_topwidth_from_depth_singular    (tB,depth(tB),setting%ZeroValue%Topwidth)
                                perimeter(tB) = llgeo_mod_basket_perimeter_from_depth_singular   (tB,depth(tB),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                                
                                hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                    (tBA, area(tBA), perimeter(tBA))
                                hydRadius(tB) = max(hydRadius(tB),setting%ZeroValue%Depth)

                                ellDepth(tBA) = llgeo_elldepth_pure (tBA)
                                ellDepth(tB)  = max(ellDepth(tB),setting%ZeroValue%Depth*0.99d0)

                            case (rectangular_closed) !% analytical
                                ! if (ii > 98) util_utest_CLprint('IIIIj')
                                area(tB)      = llgeo_rectangular_closed_area_from_depth_singular      (tB, depth(tB),setting%ZeroValue%Area)
                                topwidth(tB)  = llgeo_rectangular_closed_topwidth_from_depth_singular  (tB, depth(tB),setting%ZeroValue%Topwidth)
                                perimeter(tB) = llgeo_rectangular_closed_perimeter_from_depth_singular (tB, depth(tB),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                                
                                hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                    (tBA, area(tBA), perimeter(tBA))
                                hydRadius(tB) = max(hydRadius(tB),setting%ZeroValue%Depth)

                                ellDepth(tBA) = llgeo_elldepth_pure (tBA) 
                                ellDepth(tB)  = max(ellDepth(tB),setting%ZeroValue%Depth*0.99d0)

                            case (rect_round)  !% analytical         
                                ! if (ii > 98) util_utest_CLprint('IIIIk')                         
                                area(tB)      = llgeo_rect_round_area_from_depth_singular       (tB,depth(tB),setting%ZeroValue%Area)
                                topwidth(tB)  = llgeo_rect_round_topwidth_from_depth_singular   (tB,depth(tB),setting%ZeroValue%Topwidth)
                                perimeter(tB) = llgeo_rect_round_perimeter_from_depth_singular  (tB,depth(tB),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                                
                                hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                    (tBA, area(tBA), perimeter(tBA))
                                hydRadius(tB) = max(hydRadius(tB),setting%ZeroValue%Depth)

                                ellDepth(tBA) = llgeo_elldepth_pure (tBA) 
                                ellDepth(tB)  = max(ellDepth(tB),setting%ZeroValue%Depth*0.99d0)

                            case (rect_triang) !% analytical           
                                ! if (ii > 98) util_utest_CLprint('IIIIl')                        
                                area(tB)      = llgeo_rectangular_triangular_area_from_depth_singular      (tB,depth(tB),setting%ZeroValue%Area)
                                topwidth(tB)  = llgeo_rectangular_triangular_topwidth_from_depth_singular  (tB,depth(tB),setting%ZeroValue%Topwidth)
                                perimeter(tB) = llgeo_rectangular_triangular_perimeter_from_depth_singular (tB,depth(tB),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                                
                                hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                    (tBA, area(tBA), perimeter(tBA))
                                hydRadius(tB) = max(hydRadius(tB),setting%ZeroValue%Depth)

                                ellDepth(tBA) = llgeo_elldepth_pure (tBA) 
                                ellDepth(tB)  = max(ellDepth(tB),setting%ZeroValue%Depth*0.99d0)
                        
                            case default
                            print *, 'CODE ERROR: geometry type unknown for # ', elemI(tB,ei_geometryType)
                            print *, 'which has key ',trim(reverseKey(elemI(tB,ei_geometryType)))
                            print *, 'in ',trim(subroutine_name)
                            call util_crashpoint(399848)
                            !return
                            !stop 399848
                        end select

                        !% --- standard for all geometries
                        !hydDepth(tBA) = llgeo_hyddepth_from_area_and_topwidth_pure &
                        !                    (tBA, area(tBA), topwidth(tBA))
                        dHdA(tB)     = oneR / topwidth(tB)

                        ! if (ii > 98) util_utest_CLprint('KKKK')
                    end if
                    !% --- universal computation of volume
                    volume(tB) = area(tB) * length(tB)

                    ! print *, 'DEPTH B  ',depth(tB), setting%ZeroValue%Depth
                    ! print *, 'AREA  B  ',area(tB), setting%ZeroValue%Area
                    ! print *, 'FLOWRATE ',flowrate(tB)
                    ! print *, 'is small depth ',elemYN(tB,eYN_isSmallDepth)
                    ! print *, ' '

                    !% --- universal computation of velocity
                    if (area(tB) > setting%ZeroValue%Area) then 
                        velocity(tB) = flowrate(tB) / area(tB)

                        ! !% -- set slightly larger depths to velocity consistent with FR=1
                        ! if ((area(tB) < tenR * setting%ZeroValue%Area) .or. &
                        !     (elldepth(tB) < tenR * setting%ZeroValue%Depth) ) then 
                        !     velocity(tB) =sign(oneR,flowrate(tB)) * sqrt(grav * elldepth(tB))
                        ! end if

                    else
                        velocity(tB) = zeroR
                    end if


                    ! if (velocity(tB) > setting%Limiter%Velocity%Maximum) then
                    !     velocity(tB) = 0.99d0 * setting%Limiter%Velocity%Maximum 
                    !     ! print *, ' '
                    !     ! print *, 'SMALL VELOCITY JB -- AREA IS ',area(tB), setting%ZeroValue%Area
                    !     ! print *, 'Depth is ',depth(tB), ellDepth(tB)
                    !     ! print *, setting%ZeroValue%Depth
                    !     ! print *, ' '
                    ! end if
                            
                end do
            end do
        end if

        !% Note, the above can only be made a concurrent loop if we replace the tM
        !% with thisP(ii) and tB with thisP(ii)+kk, which makes the code
        !% difficult to read.

        if (setting%Profile%useYN) call util_profiler_stop (pfc_geo_assign_JB_from_head)

        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_assign_JB_from_head
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_JM_depth_area_from_volume ()
        !%------------------------------------------------------------------
        !% Description:
        !% Assigns depths and storage data for JM
        !% --- our approach is to use the volume (computed from conservative Q)
        !%     to get our functional and tabular plan areas.
        !%
        !%     The depth is simply the head - Zbottom, corrected for the maximum
        !%     allowed 
        !% 
        !%     Note that head, depth, and volume may be slightly inconsistent
        !%     
        !%------------------------------------------------------------------
        !% Declarations"
            integer, pointer    :: thisCol, Npack, thisP(:)
            integer, pointer    :: elemPGx(:,:), npack_elemPGx(:), col_elemPGx(:)
        !%------------------------------------------------------------------
        !% Aliases
        !% --- set packed column for updated elements
            elemPGx                => elemPGetm(:,:)
            npack_elemPGx          => npack_elemPGetm(:)
            col_elemPGx            => col_elemPGetm(:)
        !%------------------------------------------------------------------    

        print *, 'OBSOLETE'
        stop 29387043
        !% --- functional storage
        thisCol => col_elemPGx(epg_JM_functionalStorage)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            thisP => elemPGx(1:Npack, thisCol)
            call storage_plan_area_from_volume (thisP, Npack)
            call storage_functional_depth_from_volume (thisP, Npack)

            elemR(thisP,er_Depth) = min(elemR(thisP,er_Depth),   &
                                        elemR(thisP,er_FullDepth))

            elemR(thisP,er_Depth) = max(elemR(thisP,er_Depth),0.99d0*setting%ZeroValue%Depth) 

            elemR(thisP,er_EllDepth) = elemR(thisP,er_Depth)                           
        end if

        !% --- tabular storage
        thisCol => col_elemPGx(epg_JM_tabularStorage)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            thisP => elemPGx(1:Npack, thisCol)
            call storage_plan_area_from_volume (thisP, Npack)
            call storage_tabular_depth_from_volume (thisP, Npack)
            elemR(thisP,er_Depth) = min(elemR(thisP,er_Head) - elemR(thisP,er_Zbottom), &
                                        elemR(thisP,er_FullDepth))
            elemR(thisP,er_Depth) = max(elemR(thisP,er_Depth),0.99d0*setting%ZeroValue%Depth) 
            elemR(thisP,er_EllDepth) = elemR(thisP,er_Depth)  
        end if

        !% --- For implies storage with fixed plan area
        thisCol => col_elemPGx(epg_JM_impliedStorage)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            thisP => elemPGx(1:Npack, thisCol)
            where (elemSR(thisP,esr_Storage_Plan_Area) > zeroR)
                elemR(thisP,er_Depth) = elemR(thisP,er_Volume) / elemSR(thisP,esr_Storage_Plan_Area)
            endwhere
            elemR(thisP,er_Depth) = min(elemR(thisP,er_Head) - elemR(thisP,er_Zbottom), &
                                        elemR(thisP,er_FullDepth))
            elemR(thisP,er_Depth) = max(elemR(thisP,er_Depth),0.99d0*setting%ZeroValue%Depth)
            elemR(thisP,er_EllDepth) = elemR(thisP,er_Depth)    
        end if

        !% --- set the depths for no storage
        thisCol => col_elemPGx(epg_JM_noStorage)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then 
            thisP => elemPGx(1:Npack, thisCol)
            !% --- note that esr_Storage_Plan_Area is zero and never changes
            elemR(thisP,er_Depth) = min(elemR(thisP,er_Head) - elemR(thisP,er_Zbottom), &
                                        elemR(thisP,er_FullDepth))
            elemR(thisP,er_Depth) = max(elemR(thisP,er_Depth),0.99d0*setting%ZeroValue%Depth)
            elemR(thisP,er_EllDepth) = elemR(thisP,er_Depth)  
        end if

    end subroutine geo_JM_depth_area_from_volume
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_topwidth_from_depth_by_type_allCC &
        (elemPGx, npack_elemPGx, col_elemPGx)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth given depth of a non-surcharged element
        !%------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: elemPGx(:,:)
            integer, target, intent(in) :: npack_elemPGx(:), col_elemPGx(:)
            integer, pointer :: Npack, thisCol, thisP(:)

            real(8), pointer :: topwidth(:)

            character(64) :: subroutine_name = 'geo_topwidth_from_depth_by_type_CC'
        !%-------------------------------------------------------------------
        !% Aliases
            topwidth => elemR(:,er_TopWidth)
        !%-------------------------------------------------------------------
        !% Preliminaries
           if (setting%Debug%File%geometry) &
               write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% cycle through different geometries

        !% --- OPEN CHANNELS ----------------------------

        !% --- IRREGULAR
        Npack => npack_elemPGx(epg_CC_irregular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_irregular)
            thisP   => elemPGx(1:Npack,thisCol)
            call irregular_topwidth_from_depth (thisP)
            topwidth(thisP) = max(topwidth(thisP),setting%ZeroValue%Topwidth)
        end if

        !% -- PARABOLIC
        Npack => npack_elemPGx(epg_CC_parabolic)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_parabolic)
            thisP   => elemPGx(1:Npack,thisCol)
            call parabolic_topwidth_from_depth (thisP)
            topwidth(thisP) = max(topwidth(thisP),setting%ZeroValue%Topwidth)
        end if

        !% -- POWERFUNCTION
        Npack => npack_elemPGx(epg_CC_power_function)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_power_function)
            thisP   => elemPGx(1:Npack,thisCol)
            call powerfunction_topwidth_from_depth (thisP)
            topwidth(thisP) = max(topwidth(thisP),setting%ZeroValue%Topwidth)
        end if

        !% --- RECTANGULAR
        Npack => npack_elemPGx(epg_CC_rectangular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular)
            thisP   => elemPGx(1:Npack,thisCol)
            call rectangular_topwidth_from_depth (thisP)
            topwidth(thisP) = max(topwidth(thisP),setting%ZeroValue%Topwidth)
        end if

        !% --- TRAPEZOIDAL
        Npack => npack_elemPGx(epg_CC_trapezoidal)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_trapezoidal)
            thisP   => elemPGx(1:Npack,thisCol)
            call trapezoidal_topwidth_from_depth (thisP)
            topwidth(thisP) = max(topwidth(thisP),setting%ZeroValue%Topwidth)
        end if

        !% --- TRIANGULAR
        Npack => npack_elemPGx(epg_CC_triangular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_triangular)
            thisP   => elemPGx(1:Npack,thisCol)
            call triangular_topwidth_from_depth (thisP)
            topwidth(thisP) = max(topwidth(thisP),setting%ZeroValue%Topwidth)
        end if



        !% --- CLOSED CONDUITS -----------------------------------------------

        !% --  ARCH
        Npack   => npack_elemPGx(epg_CC_arch)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_arch)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, TArch, setting%ZeroValue%Topwidth)
        end if

        !% -- BASKET_HANDLE
        Npack => npack_elemPGx(epg_CC_basket_handle)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_basket_handle)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, TBasketHandle, setting%ZeroValue%Topwidth)
        end if

        !% -- CATENARY
        Npack => npack_elemPGx(epg_CC_catenary)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_catenary)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, TCatenary, setting%ZeroValue%Topwidth)
        end if

        !% --- CIRCULAR
        Npack => npack_elemPGx(epg_CC_circular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_circular)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, TCirc, setting%ZeroValue%Topwidth)
        end if

        !% -- EGG_SHAPED
        Npack => npack_elemPGx(epg_CC_egg_shaped)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_egg_shaped)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, TEgg, setting%ZeroValue%Topwidth)
        end if

        !% --- FILLED CIRCULAR
        Npack => npack_elemPGx(epg_CC_filled_circular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_filled_circular)
            thisP   => elemPGx(1:Npack,thisCol)
            call filled_circular_topwidth_from_depth (thisP)
            topwidth(thisP) = max(topwidth(thisP),setting%ZeroValue%Topwidth)
        end if

        !% -- GOTHIC
        Npack => npack_elemPGx(epg_CC_gothic)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_gothic)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, TGothic, setting%ZeroValue%Topwidth)
        end if

        !% --  HORIZONTAL ELLIPSE
        Npack   => npack_elemPGx(epg_CC_horiz_ellipse)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_horiz_ellipse)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, THorizEllip, setting%ZeroValue%Topwidth)
        end if

        !% -- HORSE_SHOE
        Npack => npack_elemPGx(epg_CC_horse_shoe)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_horse_shoe)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, THorseShoe, setting%ZeroValue%Topwidth)
        end if

        !% -- MODIFIED BASKET HANDLE
        Npack => npack_elemPGx(epg_CC_mod_basket)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_mod_basket)
            thisP   => elemPGx(1:Npack,thisCol)
            call mod_basket_topwidth_from_depth (thisP)
            topwidth(thisP) = max(topwidth(thisP),setting%ZeroValue%Topwidth)
        end if

        !% --- RECTANGULAR CLOSED
        Npack => npack_elemPGx(epg_CC_rectangular_closed)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_closed)
            thisP   => elemPGx(1:Npack,thisCol)
            call rectangular_closed_topwidth_from_depth (thisP)
            topwidth(thisP) = max(topwidth(thisP),setting%ZeroValue%Topwidth)
        end if

        !% --  RECTANGULAR ROUND
        Npack   => npack_elemPGx(epg_CC_rectangular_round)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_round)
            thisP   => elemPGx(1:Npack,thisCol)
            call rect_round_topwidth_from_depth (thisP)
            topwidth(thisP) = max(topwidth(thisP),setting%ZeroValue%Topwidth)
        end if
        
        !% -- RECTANGULAR TRIANGULAR
        Npack => npack_elemPGx(epg_CC_rectangular_triangular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_triangular)
            thisP   => elemPGx(1:Npack,thisCol)
            call rectangular_triangular_topwidth_from_depth (thisP)
            topwidth(thisP) = max(topwidth(thisP),setting%ZeroValue%Topwidth)
        end if

        !% -- SEMI-CIRCULAR
        Npack => npack_elemPGx(epg_CC_semi_circular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_semi_circular)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, TSemiCircular, setting%ZeroValue%Topwidth)
        end if

        !% -- SEMI-ELLIPTICAL
        Npack => npack_elemPGx(epg_CC_semi_elliptical)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_semi_elliptical)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, TSemiEllip, setting%ZeroValue%Topwidth)
        end if        

        !% --  VERTICAL ELLIPSE
        
        Npack   => npack_elemPGx(epg_CC_vert_ellipse)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_vert_ellipse)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, TVertEllip, setting%ZeroValue%Topwidth)

        end if

        !% TOPWIDTH NOT DEFINED FOR TABULAR, FUNCTIONAL, IMPLIED STORAGE

        !%-------------------------------------------------------------------
            if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_topwidth_from_depth_by_type_allCC
 !%
!%========================================================================== 
!%==========================================================================  
!%
    subroutine geo_topwidth_from_depth_by_element_CC (thisP, Npack)   
        !%------------------------------------------------------------------
        !% Description
        !% Companion to geo_topwidth_from_depth_by_type_allCC that
        !% computes for the entire set of cells of a given type.  This
        !% cycles through the set, so is less efficient, but is needed
        !% where only a subset of elements to be evaluated
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:), Npack
            integer :: mm
            integer, dimension(1) :: ap
            real(8), pointer :: breadthMax(:)
        !%------------------------------------------------------------------
            breadthMax => elemR(:,er_BreadthMax)
        !%------------------------------------------------------------------

            do mm=1,Npack
                ap(1) = thisP(mm)
                select case (elemI(ap(1),ei_geometryType))
                    case (irregular)
                        call irregular_topwidth_from_depth (ap)

                    case (power_function)
                        print *, 'POWER FUNCTION CROSS-SECTION NOT COMPLETE'
                        call util_crashpoint(5559867)

                    case (parabolic)
                        call parabolic_topwidth_from_depth (ap)

                    case (rectangular)
                        call rectangular_topwidth_from_depth (ap)

                    case (trapezoidal)
                        call trapezoidal_topwidth_from_depth (ap)

                    case (triangular)
                        call triangular_topwidth_from_depth (ap)

                    case (arch)
                        elemR(ap(1),er_Topwidth) &
                            = llgeo_tabular_from_depth_singular (                      &
                                ap(1), elemR(ap(1),er_Depth), breadthMax(ap(1)),       &
                                setting%ZeroValue%Depth, setting%ZeroValue%Topwidth,   &
                                TArch)

                    case (basket_handle)
                        elemR(ap(1),er_Topwidth) &
                            = llgeo_tabular_from_depth_singular (                      &
                                ap(1), elemR(ap(1),er_Depth), breadthMax(ap(1)),       &
                                setting%ZeroValue%Depth, setting%ZeroValue%Topwidth,   &
                                TBasketHandle)

                    case (catenary)
                        elemR(ap(1),er_Topwidth) &
                            = llgeo_tabular_from_depth_singular (                      &
                                ap(1), elemR(ap(1),er_Depth), breadthMax(ap(1)),       &
                                setting%ZeroValue%Depth, setting%ZeroValue%Topwidth,   &
                                TCatenary)

                    case (circular)
                        elemR(ap(1),er_Topwidth) &
                            = llgeo_tabular_from_depth_singular (                      &
                                ap(1), elemR(ap(1),er_Depth), breadthMax(ap(1)),       &
                                setting%ZeroValue%Depth, setting%ZeroValue%Topwidth,   &
                                TCirc)


                    case (eggshaped)
                        elemR(ap(1),er_Topwidth) &
                            = llgeo_tabular_from_depth_singular (                      &
                                ap(1), elemR(ap(1),er_Depth), breadthMax(ap(1)),       &
                                setting%ZeroValue%Depth, setting%ZeroValue%Topwidth,   &
                                TEgg)

                    case (filled_circular)
                        call filled_circular_topwidth_from_depth (ap)

                    case (gothic)
                        elemR(ap(1),er_Topwidth) &
                            = llgeo_tabular_from_depth_singular (                      &
                                ap(1), elemR(ap(1),er_Depth), breadthMax(ap(1)),       &
                                setting%ZeroValue%Depth, setting%ZeroValue%Topwidth,   &
                                TGothic)

                    case (horiz_ellipse)
                        elemR(ap(1),er_Topwidth) &
                            = llgeo_tabular_from_depth_singular (                      &
                                ap(1), elemR(ap(1),er_Depth), breadthMax(ap(1)),       &
                                setting%ZeroValue%Depth, setting%ZeroValue%Topwidth,   &
                                THorizEllip)

                    case (horseshoe)
                        elemR(ap(1),er_Topwidth) &
                            = llgeo_tabular_from_depth_singular (                      &
                                ap(1), elemR(ap(1),er_Depth), breadthMax(ap(1)),       &
                                setting%ZeroValue%Depth, setting%ZeroValue%Topwidth,   &
                                THorseShoe)

                    case (mod_basket)
                        call mod_basket_topwidth_from_depth (ap)

                    case (rectangular_closed)
                        call rectangular_closed_topwidth_from_depth (ap)

                    case (rect_round)
                        call rect_round_topwidth_from_depth (ap)

                    case (rect_triang)
                        call rectangular_triangular_topwidth_from_depth (ap)

                    case (semi_circular)
                        elemR(ap(1),er_Topwidth) &
                            = llgeo_tabular_from_depth_singular (                      &
                                ap(1), elemR(ap(1),er_Depth), breadthMax(ap(1)),       &
                                setting%ZeroValue%Depth, setting%ZeroValue%Topwidth,   &
                                TSemiCircular)

                    case (semi_elliptical)
                        elemR(ap(1),er_Topwidth) &
                            = llgeo_tabular_from_depth_singular (                      &
                                ap(1), elemR(ap(1),er_Depth), breadthMax(ap(1)),       &
                                setting%ZeroValue%Depth, setting%ZeroValue%Topwidth,   &
                                TSemiEllip)

                    case (vert_ellipse)
                        elemR(ap(1),er_Topwidth) &
                            = llgeo_tabular_from_depth_singular (                      &
                                ap(1), elemR(ap(1),er_Depth), breadthMax(ap(1)),       &
                                setting%ZeroValue%Depth, setting%ZeroValue%Topwidth,   &
                                TVertEllip)

                    case (custom)
                        print *, 'CUSTOM CROSS-SECTION NOT COMPLETE'
                        call util_crashpoint(52498767)
                    case default
                        print *, 'CODE ERROR: Unexpected case default'
                        call util_crashpoint(7209874)
                end select
    
            end do

    end subroutine geo_topwidth_from_depth_by_element_CC
!%
!%========================================================================== 
!%==========================================================================  
!%
    subroutine geo_perimeter_and_hydradius_from_depth_by_type_CC &
        (elemPGx, npack_elemPGx, col_elemPGx)

        !%-------------------------------------------------------------------
        !% Description:
        !% Computes the wetted perimeter and hydraulic radius for different
        !% types of cross-section. Note that for analytical cross sections
        !% the perimeter is computed first. For lookup cross-sections that
        !% hydraulic radius is computed first. The form of lookup depends
        !% on whether the hydraulic radius is tabulated directly using depth
        !% or if the section factor must be used to get hydraulic radius.
        !%-------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: elemPGx(:,:)
            integer, target, intent(in) :: npack_elemPGx(:), col_elemPGx(:)
            integer, pointer :: Npack, thisCol, thisP(:)
            real(8), pointer :: area(:), perimeter(:), hydradius(:)
            character(64) :: subroutine_name = 'geo_perimeter_and_hydradius_from_depth'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%geometry) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% Aliases
            area      => elemR(:,er_Area)
            hydradius => elemR(:,er_HydRadius)
            perimeter => elemR(:,er_Perimeter)
            
        !%-------------------------------------------------------------------
        !% cycle through different geometries

        !% --- OPEN CHANNELS -------------------------------

        !% --- IRREGULAR
        thisCol =>   col_elemPGx(epg_CC_irregular)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            thisP   => elemPGx(1:Npack,thisCol)
            !% --- solve both perimeter and hydraulic radius
            call irregular_perimeter_and_hydradius_from_depth (thisP,  &
                 setting%ZeroValue%TopWidth + setting%ZeroValue%Depth, &
                 setting%ZeroValue%Depth)
            !call irregular_perimeter_from_hydradius_area (elemPGx, Npack, thisCol)
        end if

        ! if (setting%Time%Step > 37466) &
        !   write(*,"(3(A,e12.5))") 'AA ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)


        !% --- PARABOLIC
        thisCol =>   col_elemPGx(epg_CC_parabolic)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            thisP   => elemPGx(1:Npack,thisCol)
            call parabolic_perimeter_from_depth (thisP)
            perimeter(thisP) = max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
            hydradius(thisP) = llgeo_hydradius_from_area_and_perimeter_pure &
                                        (thisP, area(thisP), perimeter(thisP))
            hydradius(thisP) = max(hydradius(thisP),setting%ZeroValue%Depth)
        end if
        ! if (setting%Time%Step > 37466) &
        ! write(*,"(3(A,e12.5))") 'BB ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)

        
        !% --- POWER FUNCTION
        thisCol => col_elemPGx(epg_CC_power_function)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            thisP   => elemPGx(1:Npack,thisCol)
            print *, 'POWER FUNCTION CROSS SECTION NOT COMPLETED'
            call util_crashpoint(54987)
            call powerfunction_perimeter_from_depth (thisP)
            perimeter(thisP) = max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
            hydRadius(thisP) = llgeo_hydradius_from_area_and_perimeter_pure &
                                         (thisP, area(thisP), perimeter(thisP))
            hydradius(thisP) = max(hydradius(thisP),setting%ZeroValue%Depth)
        end if

        !% --- RECTANGULAR
        thisCol => col_elemPGx(epg_CC_rectangular)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            thisP   => elemPGx(1:Npack,thisCol)
            call rectangular_perimeter_from_depth (thisP)
            perimeter(thisP) = max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
            hydRadius(thisP) = llgeo_hydradius_from_area_and_perimeter_pure &
                                        (thisP, area(thisP), perimeter(thisP))
            hydradius(thisP) = max(hydradius(thisP),setting%ZeroValue%Depth)
        end if

        ! if (setting%Time%Step > 37466) &
        ! write(*,"(3(A,e12.5))") 'CC ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)

        !% --- TRAPEZOIDAL
        Npack => npack_elemPGx(epg_CC_trapezoidal)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_trapezoidal)
            thisP   => elemPGx(1:Npack,thisCol)
            call trapezoidal_perimeter_from_depth (thisP)
            perimeter(thisP) = max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
            hydRadius(thisP) = llgeo_hydradius_from_area_and_perimeter_pure &
                                        (thisP, area(thisP), perimeter(thisP))
            hydradius(thisP) = max(hydradius(thisP),setting%ZeroValue%Depth)
        end if

        ! if (setting%Time%Step > 37466) &
        ! write(*,"(3(A,e12.5))") 'DD ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)

        !% --- TRIANGULAR
        Npack => npack_elemPGx(epg_CC_triangular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_triangular)
            thisP   => elemPGx(1:Npack,thisCol)
            call triangular_perimeter_from_depth (thisP)
            perimeter(thisP) = max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
            hydRadius(thisP) = llgeo_hydradius_from_area_and_perimeter_pure &
                                        (thisP, area(thisP), perimeter(thisP))
            hydradius(thisP) = max(hydradius(thisP),setting%ZeroValue%Depth)
        end if

        ! if (setting%Time%Step > 37466) &
        !   write(*,"(3(A,e12.5))") 'EE ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)

        
        !% --- CLOSED CONDUITS ---------------------------------

        !% --- ARCH
        Npack => npack_elemPGx(epg_CC_arch)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_arch)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_depth  &
                (thisP, RArch, setting%ZeroValue%Depth)
            perimeter(thisP) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            perimeter(thisP) = max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

        end if
        ! if (setting%Time%Step > 37466) &
        ! write(*,"(3(A,e12.5))") 'FF ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)


        !% -- BASKET_HANDLE
        Npack   => npack_elemPGx(epg_CC_basket_handle)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_basket_handle)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_depth  &
                (thisP, RBasketHandle, setting%ZeroValue%Depth)
            perimeter(thisP) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            perimeter(thisP)=  max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
 
        end if

        ! if (setting%Time%Step > 37466) &
        ! write(*,"(3(A,e12.5))") 'GG ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)

        !% -- CATENARY
        Npack => npack_elemPGx(epg_CC_catenary)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_catenary)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_sectionfactor_and_area &
                (thisP, SCatenary, esgr_Catenary_SoverSfull, &
                setting%ZeroValue%Depth, setting%ZeroValue%Area)
            perimeter(thisP) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            perimeter(thisP) = max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

        end if

        ! if (setting%Time%Step > 37466) &
        ! write(*,"(3(A,e12.5))") 'HH ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)

        !% --- CIRCULAR
        Npack => npack_elemPGx(epg_CC_circular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_circular)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_depth  &
                (thisP, RCirc, setting%ZeroValue%Depth)
            perimeter(thisP) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            perimeter(thisP) = max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
 
        end if

        ! if (setting%Time%Step > 37466) &
        ! write(*,"(3(A,e12.5))") 'II ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)


        !% -- EGG_SHAPED
        Npack => npack_elemPGx(epg_CC_egg_shaped)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_egg_shaped)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_depth  &
                (thisP, REgg, setting%ZeroValue%Depth)
            perimeter(thisP) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            perimeter(thisP)= max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
 
        end if

        ! if (setting%Time%Step > 37466) &
        ! write(*,"(3(A,e12.5))") 'JJ ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)

        !% --- FILLED CIRCULAR
        Npack => npack_elemPGx(epg_CC_filled_circular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_filled_circular)
            thisP   => elemPGx(1:Npack,thisCol)
            call filled_circular_hydradius_and_perimeter_from_depth (thisP)
            perimeter(thisP) = max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
            hydradius(thisP) = max(hydradius(thisP),setting%ZeroValue%Depth)

        end if
        ! if (setting%Time%Step > 37466) &
        ! write(*,"(3(A,e12.5))") 'KK ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)


        !% -- GOTHIC
        Npack => npack_elemPGx(epg_CC_gothic)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_gothic)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_sectionfactor_and_area &
                (thisP, SGothic, esgr_Gothic_SoverSfull, &
                setting%ZeroValue%Depth, setting%ZeroValue%Area)
            perimeter(thisP) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            perimeter(thisP) = max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

        end if

        ! if (setting%Time%Step > 37466) &
        ! write(*,"(3(A,e12.5))") 'LL ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)


        !% -- HORIZONTAL ELLIPSE
        Npack => npack_elemPGx(epg_CC_horiz_ellipse)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_horiz_ellipse)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_depth  &
                (thisP, RHorizEllip, setting%ZeroValue%Depth)
            perimeter(thisP) = llgeo_perimeter_from_hydradius_and_area_pure &
                                           (thisP, hydradius(thisP), area(thisP))
            perimeter(thisP) = max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

        end if

        ! if (setting%Time%Step > 37466) &
        ! write(*,"(3(A,e12.5))") 'MM ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)


        !% -- HORSE_SHOE
        Npack => npack_elemPGx(epg_CC_horse_shoe)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_horse_shoe)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_depth  &
                (thisP, RHorseShoe, setting%ZeroValue%Depth)
            perimeter(thisP) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            perimeter(thisP)= max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

        end if

        ! if (setting%Time%Step > 37466) &
        ! write(*,"(3(A,e12.5))") 'NN ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)


        !% -- MODIFIED BASKET HANDLE
        Npack => npack_elemPGx(epg_CC_mod_basket)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_mod_basket)
            thisP   => elemPGx(1:Npack,thisCol)
            call mod_basket_perimeter_from_depth (thisP)
            perimeter(thisP) = max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
            hydradius(thisP) = llgeo_hydradius_from_area_and_perimeter_pure &
                                         (thisP, area(thisP), perimeter(thisP))
            hydradius(thisP) = max(hydradius(thisP),setting%ZeroValue%Depth)
                                         
        end if

        ! if (setting%Time%Step > 37466) &
        ! write(*,"(3(A,e12.5))") 'OO ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)


        !% --- RECTANGULAR CLOSED
        Npack => npack_elemPGx(epg_CC_rectangular_closed)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_closed)
            thisP   => elemPGx(1:Npack,thisCol)
            call rectangular_closed_perimeter_from_depth (thisP)
            perimeter(thisP) = max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
            hydradius(thisP) = llgeo_hydradius_from_area_and_perimeter_pure &
                                         (thisP, area(thisP), perimeter(thisP))
            hydradius(thisP) = max(hydradius(thisP),setting%ZeroValue%Depth)
        end if

        ! if (setting%Time%Step > 37466) &
        ! write(*,"(3(A,e12.5))") 'PP ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)

        !% --  RECTANGULAR ROUND
        Npack   => npack_elemPGx(epg_CC_rectangular_round)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_round)
            thisP   => elemPGx(1:Npack,thisCol)
            call rect_round_perimeter_from_depth (thisP)
            perimeter(thisP) = max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
            hydradius(thisP) = llgeo_hydradius_from_area_and_perimeter_pure &
                                         (thisP, area(thisP), perimeter(thisP))
            hydradius(thisP) = max(hydradius(thisP),setting%ZeroValue%Depth)
        end if

        ! if (setting%Time%Step > 37466) &
        ! write(*,"(3(A,e12.5))") 'QQ ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)


        !% -- RECTANGULAR TRIANGULAR
        Npack => npack_elemPGx(epg_CC_rectangular_triangular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_triangular)
            thisP   => elemPGx(1:Npack,thisCol)
            call rectangular_triangular_perimeter_from_depth (thisP)
            perimeter(thisP) = max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
            hydradius(thisP) = llgeo_hydradius_from_area_and_perimeter_pure &
                                        (thisP, area(thisP), perimeter(thisP))
            hydradius(thisP) = max(hydradius(thisP),setting%ZeroValue%Depth)
        end if

        ! if (setting%Time%Step > 37466) &
        ! write(*,"(3(A,e12.5))") 'RR ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)

        !% -- SEMI-CIRCULAR
        Npack => npack_elemPGx(epg_CC_semi_circular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_semi_circular)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_sectionfactor_and_area &
                (thisP, SSemiCircular, esgr_Semi_Circular_SoverSfull, &
                setting%ZeroValue%Depth, setting%ZeroValue%Area)
            perimeter(thisP) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            perimeter(thisP) = max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

        end if

        ! if (setting%Time%Step > 37466) &
        ! write(*,"(3(A,e12.5))") 'SS ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)

        !% -- SEMI-ELLIPTICAL
        Npack => npack_elemPGx(epg_CC_semi_elliptical)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_semi_elliptical)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_sectionfactor_and_area &
                (thisP, SSemiEllip, esgr_Semi_Elliptical_SoverSfull, &
                setting%ZeroValue%Depth, setting%ZeroValue%Area)
            perimeter(thisP) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            perimeter(thisP)= max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

        end if

        ! if (setting%Time%Step > 37466) &
        ! write(*,"(3(A,e12.5))") 'TT ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)


        !% -- VERTICAL ELLIPSE
        Npack   => npack_elemPGx(epg_CC_vert_ellipse)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_vert_ellipse)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_depth  &
                (thisP, RVertEllip, setting%ZeroValue%Depth)
            perimeter(thisP) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            perimeter(thisP) = max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

        end if

        ! if (setting%Time%Step > 37466) &
        ! write(*,"(3(A,e12.5))") 'UU ',elemR(1624,er_HydRadius), ' ', elemR(1624,er_Depth), ' ', elemR(1624,er_Perimeter)

        !% TOPWIDTH is UNDEFINED FOR TABULAR, FUNCTIONAL, AND IMPLIED STORAGE

        !%-------------------------------------------------------------------
        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_perimeter_and_hydradius_from_depth_by_type_CC
!%
!%========================================================================== 
!%==========================================================================  
!%
    subroutine geo_perimeter_and_hydradius_from_depth_by_element_CC (thisP, Npack)   
        !%------------------------------------------------------------------
        !% Description
        !% Companion to geo_perimeter_and_hydradius_from_depth_by_type_allCC that
        !% computes for the entire set of cells of a given type.  This
        !% cycles through the set, so is less efficient, but is needed
        !% where only a subset of elements to be evaluated
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisP(:), Npack
            integer :: mm
            integer, dimension(1) :: ap

            real(8), pointer :: area(:), perimeter(:), hydradius(:)
        !%------------------------------------------------------------------
        !% Aliases
            area      => elemR(:,er_Area)
            hydradius => elemR(:,er_HydRadius)
            perimeter => elemR(:,er_Perimeter)
        !%------------------------------------------------------------------

        do mm=1,Npack
            ap(1) = thisP(mm)
            select case (elemI(ap(1),ei_geometryType))
                case (irregular)
                    call irregular_perimeter_and_hydradius_from_depth (ap, &
                        setting%ZeroValue%TopWidth + setting%ZeroValue%Depth,     &
                        setting%ZeroValue%Depth)

                case (power_function)
                    print *, 'POWER FUNCTION CROSS-SECTION NOT COMPLETE'
                    call util_crashpoint(5559867)

                case (parabolic)
                    call parabolic_perimeter_from_depth (ap)
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                    hydradius(ap) = llgeo_hydradius_from_area_and_perimeter_pure &
                                        (ap, area(ap), perimeter(ap))
                    hydradius(ap) = max(hydradius(ap),setting%ZeroValue%Depth)

                case (rectangular)
                    call rectangular_perimeter_from_depth (ap)
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                    hydRadius(ap) = llgeo_hydradius_from_area_and_perimeter_pure &
                                        (ap, area(ap), perimeter(ap))
                    hydradius(ap) = max(hydradius(ap),setting%ZeroValue%Depth)

                case (trapezoidal)
                    call trapezoidal_perimeter_from_depth (ap)
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                    hydRadius(ap) = llgeo_hydradius_from_area_and_perimeter_pure &
                                        (ap, area(ap), perimeter(ap))
                    hydradius(ap) = max(hydradius(ap),setting%ZeroValue%Depth)

                case (triangular)
                    call triangular_perimeter_from_depth (ap)
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                    hydRadius(ap) = llgeo_hydradius_from_area_and_perimeter_pure &
                                        (ap, area(ap), perimeter(ap))
                    hydradius(ap) = max(hydradius(ap),setting%ZeroValue%Depth)

                case (arch)
                    call llgeo_tabular_hydradius_from_depth (ap, RArch, setting%ZeroValue%Depth)
                    perimeter(ap) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (ap, hydradius(ap), area(ap))
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

                case (basket_handle)
                    call llgeo_tabular_hydradius_from_depth (ap, RBasketHandle, setting%ZeroValue%Depth)
                    perimeter(ap) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (ap, hydradius(ap), area(ap))
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

                case (catenary)
                    call llgeo_tabular_hydradius_from_sectionfactor_and_area &
                        (ap, SCatenary, esgr_Catenary_SoverSfull,       &
                        setting%ZeroValue%Depth, setting%ZeroValue%Area)
                    perimeter(ap) = llgeo_perimeter_from_hydradius_and_area_pure &
                                                    (ap, hydradius(thisP), area(ap))
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

                case (circular)
                    call llgeo_tabular_hydradius_from_depth (ap, RCirc, setting%ZeroValue%Depth)
                    perimeter(ap) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (ap, hydradius(ap), area(ap))
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

                case (eggshaped)
                    call llgeo_tabular_hydradius_from_depth (ap, REgg, setting%ZeroValue%Depth)
                    perimeter(ap) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (ap, hydradius(ap), area(ap))
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

                case (filled_circular)
                    call filled_circular_hydradius_and_perimeter_from_depth (ap)
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                    hydradius(ap) = max(hydradius(ap),setting%ZeroValue%Depth)

                case (gothic)
                    call llgeo_tabular_hydradius_from_sectionfactor_and_area &
                        (ap, SGothic, esgr_Gothic_SoverSfull,       &
                        setting%ZeroValue%Depth, setting%ZeroValue%Area)
                    perimeter(ap) = llgeo_perimeter_from_hydradius_and_area_pure &
                                             (ap, hydradius(thisP), area(ap))
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

                case (horiz_ellipse)
                    call llgeo_tabular_hydradius_from_depth (ap, RHorizEllip, setting%ZeroValue%Depth)
                    perimeter(ap) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (ap, hydradius(ap), area(ap))
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

                case (horseshoe)
                    call llgeo_tabular_hydradius_from_depth (ap, RHorseShoe, setting%ZeroValue%Depth)
                    perimeter(ap) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (ap, hydradius(ap), area(ap))
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

                case (mod_basket)
                    call mod_basket_perimeter_from_depth (ap)
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                    hydradius(ap) = llgeo_hydradius_from_area_and_perimeter_pure &
                                                (ap, area(ap), perimeter(ap))
                    hydradius(ap) = max(hydradius(ap),setting%ZeroValue%Depth)

                case (rectangular_closed)
                    call rectangular_closed_perimeter_from_depth (ap)
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                    hydradius(ap) = llgeo_hydradius_from_area_and_perimeter_pure &
                                         (ap, area(ap), perimeter(ap))
                    hydradius(ap) = max(hydradius(ap),setting%ZeroValue%Depth)

                case (rect_round)
                    call rect_round_perimeter_from_depth (ap)
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                    hydradius(ap) = llgeo_hydradius_from_area_and_perimeter_pure &
                                         (ap, area(ap), perimeter(ap))
                    hydradius(ap) = max(hydradius(ap),setting%ZeroValue%Depth)

                case (rect_triang)
                    call rectangular_triangular_perimeter_from_depth (ap)
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                    hydradius(ap) = llgeo_hydradius_from_area_and_perimeter_pure &
                                        (ap, area(ap), perimeter(ap))
                    hydradius(ap) = max(hydradius(ap),setting%ZeroValue%Depth)

                case (semi_circular)
                    call llgeo_tabular_hydradius_from_sectionfactor_and_area &
                        (ap, SSemiCircular, esgr_Semi_Circular_SoverSfull,       &
                         setting%ZeroValue%Depth, setting%ZeroValue%Area)
                    perimeter(ap) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (ap, hydradius(thisP), area(ap))
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)        

                case (semi_elliptical)
                    call llgeo_tabular_hydradius_from_sectionfactor_and_area &
                        (ap, SSemiEllip, esgr_Semi_Elliptical_SoverSfull,       &
                        setting%ZeroValue%Depth, setting%ZeroValue%Area)
                    perimeter(ap) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (ap, hydradius(thisP), area(ap))
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth) 

                case (vert_ellipse)
                    call llgeo_tabular_hydradius_from_depth (ap, RVertEllip, setting%ZeroValue%Depth)
                    perimeter(ap) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (ap, hydradius(ap), area(ap))
                    perimeter(ap) = max(perimeter(ap),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

                case (custom)
                    print *, 'CUSTOM CROSS-SECTION NOT COMPLETE'
                    call util_crashpoint(52498767)
                case default
                    print *, 'CODE ERROR: Unexpected case default'
                    call util_crashpoint(7209874)
            end select



        end do

    end subroutine geo_perimeter_and_hydradius_from_depth_by_element_CC
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine geo_perimeter_from_depth_by_type &
    !     (elemPGx, npack_elemPGx, col_elemPGx)
    !     !%-------------------------------------------------------------------
    !     !% Description:
    !     !% Computes the wetted perimeter given depth of a non-surcharged element
    !     !%-------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: elemPGx(:,:)
    !         integer, target, intent(in) :: npack_elemPGx(:), col_elemPGx(:)
    !         integer, pointer :: Npack, thisCol

    !         character(64) :: subroutine_name = 'geo_perimeter_from_depth'
    !     !%-------------------------------------------------------------------
    !     !% Preliminaries
    !         if (setting%Debug%File%geometry) &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%-------------------------------------------------------------------
    !     !% cycle through different geometries

    !     !% --- OPEN CHANNELS -------------------------------

    !     !% --- RECTANGULAR
    !     Npack => npack_elemPGx(epg_CC_rectangular)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_rectangular)
    !         call rectangular_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% --- TRAPEZOIDAL
    !     Npack => npack_elemPGx(epg_CC_trapezoidal)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_trapezoidal)
    !         call trapezoidal_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% --- TRIANGULAR
    !     Npack => npack_elemPGx(epg_CC_triangular)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_triangular)
    !         call triangular_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% --- PARABOLIC
    !     Npack => npack_elemPGx(epg_CC_parabolic)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_parabolic)
    !         call parabolic_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% --- POWER FUNCTION
    !     Npack => npack_elemPGx(epg_CC_power_function)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_power_function)
    !         print *, 'POWER FUNCTION CROSS SECTION NOT COMPLETED'
    !         call util_crashpoint(54987)
    !         !call power_function_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% --- IRREGULAR
    !     !%     note this requires first using the table lookup for hydraulic radius
    !     Npack => npack_elemPGx(epg_CC_irregular)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_irregular)
    !         call irregular_hydradius_from_depth (elemPGx, Npack, thisCol)
    !         call irregular_perimeter_from_hydradius_area (elemPGx, Npack, thisCol)
    !     end if

    !     !% --- CLOSED CONDUITS ---------------------------------

    !     !% --- RECTANGULAR CLOSED
    !     Npack => npack_elemPGx(epg_CC_rectangular_closed)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_rectangular_closed)
    !         call rectangular_closed_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if
        
    !     !% --  RECTANGULAR ROUND
    !     Npack   => npack_elemPGx(epg_CC_rectangular_round)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_rectangular_round)
    !         call rect_round_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% -- RECTANGULAR TRIANGULAR
    !     Npack => npack_elemPGx(epg_CC_rectangular_triangular)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_rectangular_triangular)
    !         call rectangular_triangular_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if
    
    !     !% --- CIRCULAR
    !     Npack => npack_elemPGx(epg_CC_circular)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_circular)
    !         call circular_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% -- SEMI-CIRCULAR
    !     Npack => npack_elemPGx(epg_CC_semi_circular)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_semi_circular)
    !         call semi_circular_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% --- FILLED CIRCULAR
    !     Npack => npack_elemPGx(epg_CC_filled_circular)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_filled_circular)
    !         call filled_circular_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !    !% --- ARCH
    !     Npack => npack_elemPGx(epg_CC_arch)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_arch)
    !         call arch_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% -- BASKET_HANDLE
    !     Npack   => npack_elemPGx(epg_CC_basket_handle)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_basket_handle)
    !         call basket_handle_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% -- CATENARY
    !     Npack => npack_elemPGx(epg_CC_catenary)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_catenary)
    !         call catenary_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% -- EGG_SHAPED
    !     Npack => npack_elemPGx(epg_CC_egg_shaped)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_egg_shaped)
    !         call egg_shaped_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% -- GOTHIC
    !     Npack => npack_elemPGx(epg_CC_gothic)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_gothic)
    !         call gothic_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% -- HORSE_SHOE
    !     Npack => npack_elemPGx(epg_CC_horse_shoe)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_horse_shoe)
    !         call horse_shoe_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% -- HORIZONTAL ELLIPSE
    !     Npack => npack_elemPGx(epg_CC_horiz_ellipse)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_horiz_ellipse)
    !         call horiz_ellipse_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% -- MODIFIED BASKET HANDLE
    !     Npack => npack_elemPGx(epg_CC_mod_basket)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_mod_basket)
    !         call mod_basket_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% -- SEMI-ELLIPTICAL
    !     Npack => npack_elemPGx(epg_CC_semi_elliptical)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_semi_elliptical)
    !         call semi_elliptical_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% -- VERTICAL ELLIPSE
    !     Npack   => npack_elemPGx(epg_CC_vert_ellipse)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_vert_ellipse)
    !         call vert_ellipse_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% HACK -- DO WE NEED TOPWIDTH FOR TABULAR, FUNCTIONAL, IMPLIED STORAGE?

    !     !%-------------------------------------------------------------------
    !         if (setting%Debug%File%geometry) &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine geo_perimeter_from_depth_by_type
!%
!%==========================================================================  
!%==========================================================================
!%
    ! subroutine geo_hyddepth_from_area_and_topwidth (thisColP)   
    !     !%-------------------------------------------------------------------
    !     !% Description:
    !     !% Computes Hydraulic Depth = Area/ Topwidth
    !     !% Assumes depth already computed for use in zeros or full
    !     !%-------------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in) :: thisColP
    !         integer, pointer    :: thisP(:), Npack
    !         real(8), pointer    :: hyddepth(:), area(:), topwidth(:)
    !         !% USE BELOW IF LIMITATIONS ARE NEEDED
    !         !real(8), pointer            :: depth(:), fulldepth(:), fullhyddepth(:)

    !         character(64) :: subroutine_name = 'geo_hyddepth_from_area_and_topwidth'
    !     !%-------------------------------------------------------------------
    !     !% Preliminaries
    !         Npack     => npack_elemP(thisColP)
    !         if (Npack < 1) return
    !         if (setting%Debug%File%geometry) &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%-------------------------------------------------------------------
    !     !% Aliases
    !             thisP    => elemP(1:Npack,thisColP)
    !         area         => elemR(:,er_Area)
    !         topwidth     => elemR(:,er_Topwidth)
    !         hyddepth     => elemR(:,er_HydDepth)
    !         !% USE BELOW IF LIMITATIONS ARE NEEDED
    !         ! depth        => elemR(:,er_Depth)
    !         ! fulldepth    => elemR(:,er_FullDepth)
    !         ! fullhyddepth => elemR(:,er_FullHydDepth)
            
    !     !%-------------------------------------------------------------------

    !     hyddepth(thisP) = area(thisP) / topwidth(thisP)

    !     !% USE BELOW IF LIMITATIONS ARE NEEDED
    !     ! where (depth(thisP) .le. setting%ZeroValue%Depth)
    !     !     hyddepth(thisP) = setting%ZeroValue%Depth
    !     ! elsewhere  (depth(thisP) .ge. fulldepth(thisP)) 
    !     !     hyddepth(thisP) = fullhyddepth(thisP)
    !     ! end where


    ! end subroutine geo_hyddepth_from_area_and_topwidth
!%
!%==========================================================================  
!%==========================================================================
!%
    ! subroutine geo_hyddepth_from_depth_or_topwidth (elemPGx, npack_elemPGx, col_elemPGx)
    !% OBSOLETE 20220930 brg
    !     !%-------------------------------------------------------------------
    !     !% Description:
    !     !% Note that hyddepth is the average depth, which is only area/topwidth
    !     !% for a simple open channel, and does not apply above midpoint in a
    !     !% conduit
    !     !%-------------------------------------------------------------------
    !         integer, intent(in) :: elemPGx(:,:)
    !         integer, target, intent(in) :: npack_elemPGx(:), col_elemPGx(:)
    !         integer, pointer :: Npack, thisCol

    !         character(64) :: subroutine_name = 'geo_hyddepth_from_depth'
    !     !%-------------------------------------------------------------------
    !     !% Preliminaries
    !         if (setting%Debug%File%geometry) &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%-------------------------------------------------------------------
    !     !% cycle through different geometries

    !     !% --- OPEN CHANNELS ----------------------------------

    !     !% --- RECTANGULAR
    !     Npack => npack_elemPGx(epg_CC_rectangular)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_rectangular)
    !         call rectangular_hyddepth_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% --- TRAPEZOIDAL
    !     Npack => npack_elemPGx(epg_CC_trapezoidal)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_trapezoidal)
    !         call trapezoidal_hyddepth_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% --- TRIANGULAR
    !     Npack => npack_elemPGx(epg_CC_triangular)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_triangular)
    !         call triangular_hyddepth_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% --- PARABOLIC
    !     Npack => npack_elemPGx(epg_CC_parabolic)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_parabolic)
    !         call parabolic_hyddepth_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% --- POWER FUNCTION
    !     Npack => npack_elemPGx(epg_CC_power_function)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_power_function)
    !         print *, 'POWER FUNCTION CROSS SECTION NOT COMPLETED'
    !         call util_crashpoint(5288987)
    !         !call power_function_hyddepth_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% --- IRREGULAR
    !     Npack => npack_elemPGx(epg_CC_irregular)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_irregular)
    !         call irregular_hyddepth_from_topwidth_area (elemPGx, Npack, thisCol)
    !     end if


    !     !% --- CLOSED CONDUITS -------------------------------------

    !     !% --- RECTANGULAR CLOSED
    !     Npack => npack_elemPGx(epg_CC_rectangular_closed)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_rectangular_closed)
    !         call rectangular_closed_hyddepth_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% --- RECTANGULAR ROUND
    !     Npack => npack_elemPGx(epg_CC_rectangular_round)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_rectangular_round)
    !         call rect_round_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
    !     end if

    !     !% -- RECTANGULAR TRIANGULAR
    !     Npack => npack_elemPGx(epg_CC_rectangular_triangular)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_rectangular_triangular)
    !         call rectangular_triangular_hyddepth_from_depth (elemPGx, Npack, thisCol)
    !     end if

    !     !% --- CIRCULAR
    !     Npack => npack_elemPGx(epg_CC_circular)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_circular)
    !         call circular_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
    !     end if

    !     !% -- SEMI-CIRCULAR
    !     Npack => npack_elemPGx(epg_CC_semi_circular)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_semi_circular)
    !         call semi_circular_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
    !     end if

    !     !% --- FILLED CIRCULAR
    !     Npack => npack_elemPGx(epg_CC_filled_circular)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_filled_circular)
    !         call filled_circular_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
    !     end if

    !     !% -- ARCH
    !     Npack => npack_elemPGx(epg_CC_arch)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_arch)
    !         call arch_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
    !     end if

    !     !% -- BASKET_HANDLE
    !     Npack => npack_elemPGx(epg_CC_basket_handle)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_basket_handle)
    !         call basket_handle_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
    !     end if

    !     !% --  CATENARY
    !     Npack   => npack_elemPGx(epg_CC_catenary)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_catenary)
    !         call catenary_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
    !     end if

    !     !% -- EGG_SHAPED
    !     Npack => npack_elemPGx(epg_CC_egg_shaped)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_egg_shaped)
    !         call egg_shaped_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
    !     end if

    !     !% --  GOTHIC
    !     Npack   => npack_elemPGx(epg_CC_gothic)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_gothic)
    !         call gothic_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
    !     end if

    !     !% -- HORSE_SHOE
    !     Npack => npack_elemPGx(epg_CC_horse_shoe)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_horse_shoe)
    !         call horse_shoe_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
    !     end if

    !     !% -- HORIZONTAL ELLIPSE
    !     Npack => npack_elemPGx(epg_CC_horiz_ellipse)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_horiz_ellipse)
    !         call horiz_ellipse_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
    !     end if

    !     !% --  MODIFIED BASKET HANDLE
    !     Npack   => npack_elemPGx(epg_CC_mod_basket)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_mod_basket)
    !         call mod_basket_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
    !     end if

    !     !% --  SEMI-ELLIPTICAL
    !     Npack   => npack_elemPGx(epg_CC_semi_elliptical)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_semi_elliptical)
    !         call semi_elliptical_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
    !     end if

    !     !% --  VERTICAL ELLIPSE
    !     Npack   => npack_elemPGx(epg_CC_vert_ellipse)
    !     if (Npack > 0) then
    !         thisCol => col_elemPGx(epg_CC_vert_ellipse)
    !         call vert_ellipse_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
    !     end if

    !     !% HACK -- DOES NOT HANDLE JM ELEMENTS (storage, junctions)

    !     !%-------------------------------------------------------------------
    !         if (setting%Debug%File%geometry) &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine geo_hyddepth_from_depth_or_topwidth
!%
!%==========================================================================
!%==========================================================================
!%
    !    subroutine geo_pressure_head_from_hyddepth (thisColP)
    !         !%------------------------------------------------------------------
    !         !% Description:
    !         !% calculates the pressure head
    !         !%------------------------------------------------------------------
    !         integer, intent(in) :: thisColP
    !         integer, pointer :: thisP(:), Npack
    !         real(8), pointer :: pressurehead(:), hyddepth(:), zbottom(:)

    !         character(64) :: subroutine_name = 'geo_pressure_head_from_hyddepth'
    !         !%------------------------------------------------------------------
    !         !% Preliminaries
    !             Npack     => npack_elemP(thisColP)
    !             if (Npack < 1) return
    !             if (setting%Debug%File%geometry) &
    !                 write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !         !%-------------------------------------------------------------------
    !         !% Aliases:
    !             thisP        => elemP(1:Npack,thisColP)
    !             pressurehead => elemR(:,er_Pressure_Head)
    !             hyddepth     => elemR(:,er_HydDepth)
    !             zbottom      => elemR(:,er_Zbottom)
    !         !%------------------------------------------------------------------
                
    !         pressurehead(thisP) = zbottom(thisP) + hyddepth(thisP)

    !         !%------------------------------------------------------------------
    !             if (setting%Debug%File%geometry) &
    !             write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     end subroutine geo_pressure_head_from_hyddepth 
!%
!%==========================================================================  
!%==========================================================================
!%
    ! subroutine geo_hydradius_from_area_perimeter (thisColP)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% sets area = volume/length which is common to all nonsurcharged elements
    !     !% Note this assumes volume has been limited by surcharge and zero values
    !     !%------------------------------------------------------------------
    !     integer, intent(in) :: thisColP
    !     integer, pointer :: thisP(:), Npack
    !     real(8), pointer :: area(:), hydradius(:), perimeter(:)

    !     character(64) :: subroutine_name = 'geo_hydradius_from_area_perimeter'
    !     !%------------------------------------------------------------------
    !     !% Preliminaries
    !         Npack     => npack_elemP(thisColP)
    !         if (Npack < 1) return
    !         if (setting%Debug%File%geometry) &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%-------------------------------------------------------------------
    !     !% Aliases:
    !         thisP     => elemP(1:Npack,thisColP)
    !         area      => elemR(:,er_Area)
    !         hydradius => elemR(:,er_HydRadius)
    !         perimeter => elemR(:,er_Perimeter)
    !     !%------------------------------------------------------------------
            
    !     hydradius(thisP) = area(thisP) / perimeter(thisP)

    !     !%------------------------------------------------------------------
    !         if (setting%Debug%File%geometry) &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine geo_hydradius_from_area_perimeter
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_elldepth_from_head_CC (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% computes the value of "ell" -- the modified hydraulic depth
        !% used as a length scale in AC method
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisP(:)
            real(8), pointer :: ellDepth(:), head(:), area(:), topwidth(:) , depth(:)
            real(8), pointer :: ZbreadthMax(:), breadthMax(:), areaBelowBreadthMax(:)
            real(8), pointer :: zcrown(:)
            integer :: ii

            character(64) :: subroutine_name = 'geo_elldepth_from_head_CC'
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (setting%Debug%File%geometry) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases:
            ellDepth            => elemR(:,er_EllDepth)
            head                => elemR(:,er_Head)
            area                => elemR(:,er_Area)
            depth               => elemR(:,er_Depth)
            topwidth            => elemR(:,er_Topwidth)
            ZbreadthMax         => elemR(:,er_ZbreadthMax)
            zcrown              => elemR(:,er_Zcrown)
            breadthMax          => elemR(:,er_BreadthMax)
            areaBelowBreadthMax => elemR(:,er_AreaBelowBreadthMax)
        !%-------------------------------------------------------------------

        ! print *, ' '
        ! print *, 'in elldepth ----------------------------------'
        ! print *, 'type ',trim(reverseKey(elemI(93,ei_geometryType)))
        ! print *, 'depth '
        ! print *, elemR(93,er_Depth)
        ! print *, ' '
        ! print *, 'head '
        ! print *, head(93)
        ! print *,' '
        ! print *, 'Zbreadthmax'
        ! print *, ZbreadthMax(93)
        ! print *, ' '
        ! print *, 'area '
        ! print *, area(93)
        ! print *, ' '
        ! print *, 'topwidth '
        ! print *, topwidth(93)
        ! print *, ' '
        ! print *, 'areaBelowBreadthMax '
        ! print *, areaBelowBreadthMax(93)
        ! print *, ' '
        ! print *, 'breadthMax '
        ! print *, breadthMax(93)    
        ! print *, ' '
            
        where ((head(thisP) .le. ZbreadthMax(thisP)) &
             .and. (topwidth(thisP) > setting%ZeroValue%Topwidth))
            !ell(thisP) =  hydDepth(thisP)
            ellDepth(thisP) =  area(thisP) / topwidth(thisP)
        elsewhere
            ellDepth(thisP) = ( (head(thisP) - ZbreadthMax(thisP)) * breadthMax(thisP) &
                  + areaBelowBreadthMax(thisP) ) / breadthMax(thisP)
        endwhere

        !% --- limiter to use for open channel flow
        !%     Needed because sometimes small depth topwidth too small,
        !%     so the area/topwidth provides ellDepth > depth, which is inconsistent.
        where (head(thisP) .le. zcrown(thisP))
            ellDepth(thisP) = min(depth(thisP), ellDepth(thisP))
        endwhere

        ellDepth(thisP) = max(ellDepth(thisP),setting%ZeroValue%Depth*0.99d0)

        !%-------------------------------------------------------------------
            if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_elldepth_from_head_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_JM_values ()
        !%------------------------------------------------------------------
        !% Description:
        !% The junction main (JM) values for HydDepth, ell,...
        !% Are not defined because geometry such as ZbreadthMax are not 
        !% defined. 
        !% Here we use the depth at the JM junctions so that we don't have
        !% nullvalueR stored here
        !%
        !% HACK
        !% the following variables are NOT defined on JM and perhaps need to
        !% be added:
        !% dHdA, FullArea, FroudeNumber, FullHydDepth, FullPerimeter,
        !% FullVolume, HydRadius, InterpWeight_xx, Length, Perimeter,
        !% Roughness, TopWidth, ZbreadthMax, Zcrown
        !%
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: thisCol, Npack, thisP(:)
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:
            thisCol => col_elemP(ep_JM)
            Npack   => npack_elemP(thisCol)
            thisP   => elemP(1:Npack,thisCol)
        !%------------------------------------------------------------------
        if (Npack > 0) then 
            !elemR(thisP,er_HydDepth)      = elemR(thisP,er_Depth)
            elemR(thisP,er_EllDepth)       = elemR(thisP,er_Depth)
            !elemR(thisP,er_Pressure_Head) = elemR(thisP,er_Depth)
        end if
        !%------------------------------------------------------------------
        !% Closing
    end subroutine geo_JM_values
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_dHdA (thisColP)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% This simply uses the inverse of the topwidth as dH/dA, which is an
        !% assumption of a small change. Arguably, for our known geometries we could be
        !% more precise, but it is not clear that it would be worth the effort.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisColP
        integer, pointer :: thisP(:), Npack
        real(8), pointer :: dHdA(:), topwidth(:)

        character(64) :: subroutine_name = 'geo_dHdA'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        Npack    => npack_elemP(thisColP)
        dHdA     => elemR(:,er_dHdA)
        topwidth => elemR(:,er_Topwidth)
        !%-----------------------------------------------------------------------------

        if (Npack > 0) then
            thisP    => elemP(1:Npack,thisColP)
            dHdA(thisP) = oneR / topwidth(thisP)
        end if

        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_dHdA
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_ACsurcharged (thisColP)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Sets volume, area, depth, perimeter, topwidth, hydraulic depth,
        !% and hydraulic radius for any surcharged conduit.
        !% Note that ell, and dHdA must be set elsewhere as they depend on specific geometry.
        !% Note the topwidth for surcharged is set to a small positive value to prevent
        !% division by zero in transition from surcharged to non-surcharged.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisColP
        integer, pointer :: Npack ,thisP(:)

        ! character(64) :: subroutine_name = 'geo_surcharged'
        ! !%-----------------------------------------------------------------------------
        ! !!if (crashYN) return
        ! Npack => npack_elemP(thisColP)
        ! !%-------------------------------------------------
        ! if (setting%Debug%File%geometry) &
        !     write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        ! if (Npack > 0) then
        !     thisP => elemP(1:Npack,thisColP)
        !     elemR(thisP,er_Volume)    = elemR(thisP,er_FullVolume)
        !     elemR(thisP,er_Area)      = elemR(thisP,er_FullArea)
        !     elemR(thisP,er_Depth)     = elemR(thisP,er_FullDepth)
        !     elemR(thisP,er_Perimeter) = elemR(thisP,er_FullPerimeter)
        !     elemR(thisP,er_HydDepth)  = elemR(thisP,er_FullHydDepth)
        !     elemR(thisP,er_HydRadius) = elemR(thisP,er_FullArea) / elemR(thisP,er_FullPerimeter)
        !     elemR(thisP,er_Topwidth)  = setting%ZeroValue%Topwidth
        ! end if

        ! if (setting%Debug%File%geometry) &
        !     write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        end subroutine geo_ACsurcharged
!%
!%==========================================================================
!%==========================================================================
!%    
    ! subroutine geo_limit_incipient_surcharge (geocol, fullcol, thisColP, isVolume) 
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Sets volume/depth limit to full volume for incipient surcharge.
    !     !% If input is volume, the excess is added to er_VolumeOverFlow
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: thisColP, geocol, fullcol
    !         logical, intent(in) :: isVolume !% 20220124brh
    !         integer, pointer :: Npack, thisP(:)
    !         real(8), pointer :: geovalue(:), fullvalue(:)
    !         real(8), pointer :: overflow(:)  !% 20220124brh

    !         character(64) :: subroutine_name = 'geo_limit_incipient_surcharge'
    !     !%------------------------------------------------------------------
    !     !% Preliminaries
    !         Npack      => npack_elemP(thisColP)
    !         if (Npack < 1) return
    !     !%-------------------------------------------------------------------
    !     !% Aliases
    !         thisP      => elemP(1:Npack,thisColP)
    !         eType      => elemR(:,ei_elementType)
    !         geovalue   => elemR(:,geocol)
    !         fullvalue  => elemR(:,fullcol)
    !         overflow   => elemR(:,er_VolumeOverFlow)  !% 20220124brh
    !         ponding    => elemR(:,er_VolumePonded)
    !     !%-------------------------------------------------------------------

    !         ! print *, 'in ',trim(subroutine_name),elemR(49,er_VolumeOverFlow)
    !         ! print *, geovalue(49), fullvalue(49), overflow(49)
        
    !     if (isVolume) then
    !         !% --- Note that Ponding and Overflow from JM has already been handled
    !         !%     in ll_JM_Slot_computation

    !         ! where (geovalue(thisP) > fullvalue(thisP))
    !         !     overflow(thisP) = geovalue(thisP) - fullvalue(thisP) + overflow(thisP)  
    !         !     geovalue(thisP) = fullvalue(thisP)
    !         ! endwhere
    !     else
    !         where (geovalue(thisP) > fullvalue(thisP))
    !             geovalue(thisP) = fullvalue(thisP)
    !         endwhere
    !     end if

    !     !print *, 'end of ',trim(subroutine_name),elemR(thisP,er_VolumeOverFlow)


    ! end subroutine geo_limit_incipient_surcharge
!%
!%==========================================================================
!% SINGULAR (PUBLIC) FUNCTIONS
!%==========================================================================
!%
    real(8) function geo_area_from_depth_singular &
        (idx, indepth, ZeroValueArea) result (outvalue)
        !%------------------------------------------------------------------
        !% Descriptions:
        !% computes the area for a given depth of a single element
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(in)   :: indepth, ZeroValueArea
            integer, intent(in)   :: idx
            integer, dimension(1) :: iA 
            real(8), dimension(1) :: depthA, outA
            real(8), pointer      :: depth(:), fullarea(:)
            real(8), pointer      :: ATable(:)
            character(64) :: subroutine_name = 'geo_area_from_depth_singular'
        !%------------------------------------------------------------------
        !% Aliases
            depth    => elemR(:,er_Depth)
            fullarea => elemR(:,er_FullArea)
        !%------------------------------------------------------------------
        !% size(1) arrays
            iA(1)     = idx
            depthA(1) = indepth
        !%------------------------------------------------------------------

        !% ensure that small depths have small perimeter
        if (indepth .le. setting%ZeroValue%Depth) then
            outvalue = setting%ZeroValue%Area
            return
        endif    

        !% --- set lookup table names
        select case (elemI(idx,ei_geometryType))
            !% --- for tables using HydRadius
            case (arch)
                Atable => AArch
            case (basket_handle)
                Atable => ABasketHandle
            case (circular)
                Atable => ACirc
            case (eggshaped)
                Atable => AEgg
            case (horiz_ellipse)
                Atable => AHorizEllip
            case (horseshoe)
                Atable => AHorseShoe
            case (vert_ellipse)
                Atable => AHorseShoe

            !% --- for tables using section factor
            case (catenary)
                Atable => ACatenary
            case (gothic)
                Atable => AGothic
            case (semi_circular)
                Atable => ASemiCircular
            case (semi_elliptical)
                Atable => ASemiEllip
        end select

        select case (elemI(idx,ei_geometryType))
        !% ----open channels  
        case (irregular)
            outvalue = irregular_geometry_from_depth_singular &
                (idx,tt_area, indepth, elemR(idx,er_FullArea), ZeroValueArea)

        case (parabolic)
            !outvalue = parabolic_area_from_depth_singular (idx, indepth)
            outA = llgeo_parabolic_area_from_depth_pure(iA, depthA)
            outvalue = outA(1)
            outvalue = max(outvalue,ZeroValueArea)

        case (power_function)
            print *, 'CODE ERROR powerfunction geometry not complete'
            call util_crashpoint(55098723)   

        case (rectangular)
            outA = llgeo_rectangular_area_from_depth_pure (iA, depthA)
            outvalue = outA(1)
            outvalue = max(outvalue,ZeroValueArea)

        case (trapezoidal)
            outA = llgeo_trapezoidal_area_from_depth_pure (iA, depthA)
            outvalue = outA(1)
            outvalue = max(outvalue,ZeroValueArea)

        case (triangular)
            outA= llgeo_triangular_area_from_depth_pure (iA, depthA)
            outvalue = max(outvalue,ZeroValueArea)
            outvalue = outA(1)

        !% --- closed conduits   
        case (arch, basket_handle, catenary, circular, eggshaped, gothic, &
            horiz_ellipse, horseshoe, semi_circular, semi_elliptical,     &
            vert_ellipse)

           ! outvalue = arch_area_from_depth_singular (idx, indepth)
            outvalue = llgeo_tabular_from_depth_singular &
                    (idx, indepth, fullArea(idx), setting%ZeroValue%Depth, ZeroValueArea, Atable)
            
        case (filled_circular)
            outvalue = llgeo_filled_circular_area_from_depth_singular (idx, indepth, ZeroValueArea)     
            
        case (mod_basket)
            outvalue = llgeo_mod_basket_area_from_depth_singular (idx, indepth, ZeroValueArea)

        case (rectangular_closed)
            outvalue = llgeo_rectangular_closed_area_from_depth_singular (idx, indepth, ZeroValueArea)

        case (rect_round )
            outvalue = llgeo_rect_round_area_from_depth_singular (idx, indepth, ZeroValueArea)

        case (rect_triang)
            outvalue = llgeo_rectangular_triangular_area_from_depth_singular (idx, indepth, ZeroValueArea)

        !case (circular )
        !    outvalue = circular_area_from_depth_singular (idx, indepth)
        !case (semi_circular)
        !    outvalue = semi_circular_area_from_depth_singular (idx, indepth)
        
        !case (basket_handle)
        !    outvalue = basket_handle_area_from_depth_singular (idx, indepth)
        !case (catenary)
        !    outvalue = catenary_area_from_depth_singular (idx, indepth)
        !case (eggshaped)
        !    outvalue = egg_shaped_area_from_depth_singular (idx, indepth)
        !case (gothic)
        !    outvalue = gothic_area_from_depth_singular (idx, indepth)
        !case (horseshoe)
        !    outvalue = horse_shoe_area_from_depth_singular (idx, indepth)
        !case (horiz_ellipse)
        !    outvalue = horiz_ellipse_area_from_depth_singular (idx, indepth)
        
        !case (semi_elliptical)
        !    outvalue = semi_elliptical_area_from_depth_singular (idx, indepth)
        !case (vert_ellipse)
        !    outvalue = vert_ellipse_area_from_depth_singular (idx, indepth)
        case (custom)
            print *, 'CODE ERROR: area for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(332341)

        case (force_main)
            print *, 'CODE ERROR: area for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'in ',trim(subroutine_name)   
            print *, 'This should never be reached as a force_main is not a valid geometryType'
            call util_crashpoint(332342)

        case default
            print *, 'CODE ERROR: area for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(332343)

        end select
           
    end function geo_area_from_depth_singular
!%
!%==========================================================================    
!%==========================================================================
!%
    real(8) function geo_topwidth_from_depth_singular &
        (idx, indepth, ZeroValueTopwidth) result (outvalue)
        !%------------------------------------------------------------------
        !% Descriptions:
        !% computes the topwidth for a given depth of a single element
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(in)  :: indepth, ZeroValueTopwidth
            integer, intent(in)  :: idx
            integer, dimension(1):: iA
            real(8), dimension(1):: depthA, outA
            real(8), pointer     :: depth(:), fullarea(:), fulltopwidth(:), breadthmax(:)
            real(8), pointer     :: TTable(:)
            character(64) :: subroutine_name = 'geo_topwidth_from_depth_singular'
        !%------------------------------------------------------------------
            depth        => elemR(:,er_Depth)
            breadthmax   => elemR(:,er_BreadthMax)
            fullarea     => elemR(:,er_FullArea)
            fulltopwidth => elemR(:,er_FullTopWidth)
        !%------------------------------------------------------------------
        !% size(1) arrays
            iA(1)     = idx
            depthA(1) = indepth
        !%------------------------------------------------------------------

        !% ensure that small depths have small topwidths
        if (indepth .le. setting%ZeroValue%Depth) then
            outvalue = setting%ZeroValue%Topwidth
            return
        endif

        select case (elemI(idx,ei_geometryType))
                !% --- for tables using HydRadius
            case (arch)
                Ttable => TArch
            case (basket_handle)
                Ttable => TBasketHandle
            case (circular)
                Ttable => TCirc
            case (eggshaped)
                Ttable => TEgg
            case (horiz_ellipse)
                Ttable => THorizEllip
            case (horseshoe)
                Ttable => THorseShoe
            case (vert_ellipse)
                Ttable => THorseShoe

            !% --- for tables using section factor
            case (catenary)
                Ttable => TCatenary
            case (gothic)
                Ttable => TGothic
            case (semi_circular)
                Ttable => TSemiCircular
            case (semi_elliptical)
                Ttable => TSemiEllip
            end select    
        select case (elemI(idx,ei_geometryType))
            
        !% --- OPEN CHANNELS    
        case (parabolic)
            outA = llgeo_parabolic_topwidth_from_depth_pure (iA, depthA)
            outvalue = outA(1)
            outvalue = max(outvalue,ZeroValueTopWidth)

        case (power_function)
            print *, 'CODE ERROR: topwidth for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)

        case (rectangular)
            outA = llgeo_rectangular_topwidth_from_depth_pure  (iA, depthA)
            outvalue = outA(1)
            outvalue = max(outvalue,ZeroValueTopWidth)

        case (trapezoidal)
            outA = llgeo_trapezoidal_topwidth_from_depth_pure (iA, depthA)
            outvalue = outA(1)
            outvalue = max(outvalue,ZeroValueTopWidth)

        case (triangular)
            outA = llgeo_triangular_topwidth_from_depth_pure (iA, depthA)
            outvalue = outA(1)
            outvalue = max(outvalue,ZeroValueTopWidth)
        
        case (irregular)
            
            outvalue = irregular_geometry_from_depth_singular &
                (idx,tt_width, indepth, elemR(idx,er_FullTopWidth), ZeroValueTopwidth)

            ! print *, 'calling irregular ',idx, tt_width, indepth
            ! print *, 'full topwidth ',elemR(idx,er_FullTopWidth), outvalue

        !% -----CLOSED CONDUITS ---------------------------------------------
        case (arch, basket_handle, catenary, circular, eggshaped, gothic,  &
                horiz_ellipse, horseshoe, semi_circular, semi_elliptical,  &
                vert_ellipse)

            outvalue = llgeo_tabular_from_depth_singular &
                (idx, depth(idx), breadthmax(idx), setting%ZeroValue%Depth,  ZeroValueTopwidth, Ttable)

            outvalue = min(outvalue, breadthMax(idx))

        case (filled_circular)
            outvalue = llgeo_filled_circular_topwidth_from_depth_singular  (idx, indepth, ZeroValueTopwidth)

        case (mod_basket)
            outvalue = llgeo_mod_basket_topwidth_from_depth_singular (idx, indepth, ZeroValueTopwidth)

        case (rectangular_closed)
            outvalue = llgeo_rectangular_closed_topwidth_from_depth_singular  (idx, indepth, ZeroValueTopwidth)

        case (rect_round)
            outvalue = llgeo_rect_round_topwidth_from_depth_singular (idx, indepth, ZeroValueTopwidth)

        case (rect_triang)
            outvalue = llgeo_rectangular_triangular_topwidth_from_depth_singular  (idx, indepth, ZeroValueTopwidth)
        
        !case (circular )
        !    outvalue = circular_topwidth_from_depth_singular  (idx, indepth)
        !case (semi_circular)
        !    outvalue = semi_circular_topwidth_from_depth_singular (idx, indepth)
        
        !case (arch)
        !    outvalue = arch_topwidth_from_depth_singular (idx, indepth)
        !case (basket_handle)
        !    outvalue = basket_handle_topwidth_from_depth_singular (idx, indepth)
        !case (catenary)
        !    outvalue = catenary_topwidth_from_depth_singular (idx, indepth)
        !case (eggshaped)
        !    outvalue = egg_shaped_topwidth_from_depth_singular (idx, indepth)
        !case (gothic)
         !   outvalue = gothic_topwidth_from_depth_singular (idx, indepth)
        !case (horseshoe)
        !    outvalue = horse_shoe_topwidth_from_depth_singular (idx, indepth)
        !case (horiz_ellipse)
        !    outvalue = horiz_ellipse_topwidth_from_depth_singular (idx, indepth)
        
        !case (semi_elliptical)
        !    outvalue = semi_elliptical_topwidth_from_depth_singular (idx, indepth)        
        !case (vert_ellipse)
         !   outvalue = vert_ellipse_topwidth_from_depth_singular (idx, indepth)
        
        case (custom)
            print *, 'CODE ERROR: topwidth for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)

        case (force_main)
            print *, 'CODE ERROR: topwidth for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'in ',trim(subroutine_name)   
            print *, 'This should never be reached as a force_main is not a valid geometryType'
            call util_crashpoint(4498734)

        case default
            print *, 'CODE ERROR: topwidth for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        end select

        
           
    end function geo_topwidth_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function geo_perimeter_from_depth_singular &
        (idx, indepth, ZeroValuePerimeter) result (outvalue)
        !%------------------------------------------------------------------
        !% Descriptions:
        !% computes the perimeter for a given depth of a single element
        !% Note that ALL values computed herein are based on the input depth
        !% which may NOT be the depth of the element idx. Thus, where we
        !% need more values than the depth to compute the perimeter we must 
        !% compute the other values as temporary values!
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(in)   :: indepth, ZeroValuePerimeter
            integer, intent(in)   :: idx
            integer, dimension(1) :: iA 
            real(8), dimension(1) :: indepthA, outA
            real(8), pointer      :: fullarea(:), fullHydRadius(:)
            real(8), pointer      :: Atable(:), Rtable(:), Stable(:), Ttable(:)
            real(8)               :: tempArea(1), tempHydRadius(1)
            character(64) :: subroutine_name = 'geo_perimeter_from_depth_singular'
        !%------------------------------------------------------------------
            fullarea      => elemR(:,er_FullArea)
            fullHydRadius => elemR(:,er_FullHydRadius)
        !%------------------------------------------------------------------

        !% ensure that small depths have small perimeter
        if (indepth .le. setting%ZeroValue%Depth) then
            outvalue = setting%ZeroValue%Topwidth
            return
        endif

        !% --- 1D array values of index and depth for pure functions
        iA(1)       = idx
        indepthA(1) = indepth

        select case (elemI(idx,ei_geometryType))
            !% --- for tables using HydRadius
            case (arch)
                Atable => AArch
                Rtable => RArch
                Ttable => TArch
            case (basket_handle)
                Atable => ABasketHandle
                Rtable => RBasketHandle
                Ttable => TBasketHandle
            case (circular)
                Atable => ACirc
                Rtable => RCirc
                Ttable => TCirc
            case (eggshaped)
                Atable => AEgg
                Rtable => REgg
                Ttable => TEgg
            case (horiz_ellipse)
                Atable => AHorizEllip
                Rtable => RHorizEllip
                Ttable => THorizEllip
            case (horseshoe)
                Atable => AHorseShoe
                Rtable => RHorseShoe
                Ttable => THorseShoe
            case (vert_ellipse)
                Atable => AHorseShoe
                Rtable => RHorseShoe
                Ttable => THorseShoe

            !% --- for tables using section factor
            case (catenary)
                Atable => ACatenary
                Stable => SCatenary
                Ttable => TCatenary
            case (gothic)
                Atable => AGothic
                Stable => SGothic
                Ttable => TGothic
            case (semi_circular)
                Atable => ASemiCircular
                Stable => SSemiCircular
                Ttable => TSemiCircular
            case (semi_elliptical)
                Atable => ASemiEllip
                Stable => SSemiEllip
                Ttable => TSemiEllip
        end select


        select case (elemI(idx,ei_geometryType))
            
        !% ---- OPEN CHANNELS 
        case (parabolic)
            outA = llgeo_parabolic_perimeter_from_depth_pure (iA, indepthA)
            outvalue = outA(1)
            outvalue = max(outvalue,ZeroValuePerimeter)
            !outvalue = parabolic_perimeter_from_depth_singular (idx, indepth)

        case (power_function)
            print *, 'CODE ERROR: perimeter for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(338234)    

        case (rectangular)
            !print *, 'in geo_perimeter_from_depth_singular ',iA, indepthA
            outA = llgeo_rectangular_perimeter_from_depth_pure (iA, indepthA)
            !print *, 'outA ',outA
            outvalue = outA(1)
            outvalue = max(outvalue,ZeroValuePerimeter)

        case (trapezoidal)
            outA = llgeo_trapezoidal_perimeter_from_depth_pure (iA, indepthA)
            outvalue = outA(1)
            outvalue = max(outvalue,ZeroValuePerimeter)

        case (triangular)
            outA = llgeo_triangular_perimeter_from_depth_pure (iA, indepthA)
            outvalue = outA(1)
            outvalue = max(outvalue,ZeroValuePerimeter)
        
        case (irregular)
            outvalue = irregular_geometry_from_depth_singular &
                (idx,tt_area, indepth, elemR(idx,er_FullArea), ZeroValuePerimeter)

        !% --- CLOSED CONDUITS   
        case (arch, basket_handle, circular, eggshaped, horiz_ellipse, horseshoe, vert_ellipse)
            tempArea(1)  = llgeo_tabular_from_depth_singular &
                    (idx, indepth, fullArea(idx), setting%ZeroValue%Depth, ZeroValuePerimeter, Atable)

            tempHydRadius(1) = llgeo_tabular_from_depth_singular &
                    (idx, indepth, fullHydRadius(idx), setting%ZeroValue%Depth, ZeroValuePerimeter, Rtable)

            outA = llgeo_perimeter_from_hydradius_and_area_pure &
                    (iA, tempHydRadius, tempArea)
            outvalue = outA(1)
            outvalue = max(outvalue,ZeroValuePerimeter)
        
        case (catenary, gothic, semi_circular, semi_elliptical)  

            tempArea(1)  = llgeo_tabular_from_depth_singular &
                    (idx, indepth, fullArea(idx), setting%ZeroValue%Depth, ZeroValuePerimeter, Atable)

            temphydRadius(1)= llgeo_tabular_hydradius_from_area_and_sectionfactor_singular &
                (idx, tempArea(1), fullhydradius(idx), setting%ZeroValue%Area, Stable)

            outA = llgeo_perimeter_from_hydradius_and_area_pure &
                        (iA, tempHydradius, tempArea)

            outvalue = outA(1)
            outvalue = max(outvalue,ZeroValuePerimeter)

        case (filled_circular)
            outvalue = llgeo_filled_circular_perimeter_from_depth_singular (idx, indepth, ZeroValuePerimeter)

        case (mod_basket)
            outvalue = llgeo_mod_basket_perimeter_from_depth_singular (idx, indepth, ZeroValuePerimeter)

        case (rectangular_closed)
            outvalue = llgeo_rectangular_closed_perimeter_from_depth_singular (idx, indepth, ZeroValuePerimeter) 

        case (rect_round )
            outvalue = llgeo_rect_round_perimeter_from_depth_singular (idx, indepth, ZeroValuePerimeter)

        case (rect_triang)
            outvalue =llgeo_rectangular_triangular_perimeter_from_depth_singular (idx, indepth, ZeroValuePerimeter)

        !case (circular )
        !    outvalue = circular_perimeter_from_depth_singular (idx, indepth)
        !case (semi_circular)
        !    outvalue = semi_circular_perimeter_from_depth_singular (idx, indepth)
        
        !case (arch)
        !    outvalue = arch_perimeter_from_depth_singular (idx, indepth)
        !case (basket_handle)
        !    outvalue = basket_handle_perimeter_from_depth_singular (idx, indepth)
        !case (catenary)
        !    outvalue = catenary_perimeter_from_depth_singular (idx, indepth)
        !case (eggshaped)
        !    outvalue = egg_shaped_perimeter_from_depth_singular (idx, indepth)
        !case (gothic)
        !    outvalue = gothic_perimeter_from_depth_singular (idx, indepth)
        !case (horseshoe)
        !    outvalue = horse_shoe_perimeter_from_depth_singular (idx, indepth)
        !case (horiz_ellipse)
            !print *, 'CODE ERROR: perimeter for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            !print *, 'has not been implemented in ',trim(subroutine_name)
            !call util_crashpoint(338234)
        !    outvalue = horiz_ellipse_perimeter_from_depth_singular (idx, indepth)
        
        !case (semi_elliptical)
         !   outvalue = semi_elliptical_perimeter_from_depth_singular (idx, indepth)
        !case (vert_ellipse)
        !    outvalue = vert_ellipse_perimeter_from_depth_singular (idx, indepth)
        
        case (custom)
            print *, 'CODE ERROR: perimeter for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(338234)
        case (force_main)
            print *, 'CODE ERROR: perimeter for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'in ',trim(subroutine_name)   
            print *, 'This should never be reached as a force_main is not a valid geometryType' 
            call util_crashpoint(338234)
        case default
            print *, 'CODE ERROR: perimeter for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(332344)
        end select



    end function geo_perimeter_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function geo_hyddepth_from_area_and_topwidth_singular &
        (idx, area, topwidth, ZeroValueHydDepth)  result (outvalue)
        !%------------------------------------------------------------------
        !% Descriptions:
        !% computes the hydraulic depth for area and topwidth of a single
        !% element
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(in)  :: area, topwidth, ZeroValueHydDepth
            integer, intent(in)  :: idx
            !real(8)              :: temp1, temp2
            character(64) :: subroutine_name = 'geo_hyddepth_from_area_and_topwidth_singular'
        !%------------------------------------------------------------------   
     
        if (topwidth > zeroR) then
            outvalue = area / topwidth
        else
            outvalue = zeroR
        end if

        outvalue = max(outvalue,ZeroValueHydDepth)
        
    end function geo_hyddepth_from_area_and_topwidth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    ! real(8) function geo_elldepth_singular &
    !     (head, area, topwidth, areaBelowMaxBreadth, maxBreadth, ZmaxBreadth) &
    !     result (outvalue)   
    !     !%------------------------------------------------------------------
    !     !% Descriptions:
    !     !% computes the modified hydraulic depth of a single
    !     !% element
    !     !%------------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in) :: idx
    !         real(8), intent(in) :: head, area, topwidth, areaBelowMaxBreadth,
    !         real(8), intent(in) :: maxBreadth, ZmaxBreadth
    !     !%------------------------------------------------------------------
        
    !     if (head .le. ZmaxBreadth) then 
    !         !% --- ell is simply the hydraulic depth below 
    !         outvalue = area / topwidth
    !     else
    !         !% -- ell is the modified hydraulic depth
    !         outvalue = (((head - ZmaxBreadth) * maxBreadth) + areaBelowMaxBreadth) &
    !                    / maxBreadth
    !     end if

    ! end function geo_elldepth_from_area_topwidth_and_maxbreadth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    ! real(8) function geo_hydradius_from_area_and_perimeter_singular &
    !     (idx, area, perimeter)  result (outvalue)
    !     !%------------------------------------------------------------------
    !     !% Descriptions:
    !     !% computes the hydraulic radius for area and perimeter of a single
    !     !% element
    !     !%------------------------------------------------------------------
    !     !% Declarations
    !         real(8), intent(in)  :: area, perimeter
    !         integer, intent(in)  :: idx
    !         !real(8)              :: temp1, temp2
    !         character(64) :: subroutine_name = 'geo_hydradius_from_area_and_perimeter_singular'
    !     !%------------------------------------------------------------------   
    
    !     outvalue = area / perimeter
    
    ! end function geo_hydradius_from_area_and_perimeter_singular
!%
!%==========================================================================
!%==========================================================================
!%
    ! real(8) function geo_perimeter_from_area_and_hydradius_singular &
    !     (idx, area, hydradius) result (outvalue)    
    !     !%------------------------------------------------------------------
    !     !% Descriptions:
    !     !% computes the perimeter for area and hydraulic radius of a single
    !     !% element
    !     !%------------------------------------------------------------------
    !     !% Declarations
    !         real(8), intent(in)  :: area, hydradius
    !         integer, intent(in)  :: idx
    !         character(64) :: subroutine_name = 'geo_perimeter_from_area_and_hydradius_singular'
    !     !%------------------------------------------------------------------   
    
    !     outvalue = area / hydradius
    
    ! end function geo_perimeter_from_area_and_hydradius_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_ZeroDepth_from_volume (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% ensures that if volume <= zeroArea * length the depth will be
        !% zeroDepth
        !%------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: depth(:), volume(:), length(:), area(:)
        !%------------- -----------------------------------------------------
        !% Aliases
            depth       => elemR(:,er_Depth)
            volume      => elemR(:,er_Volume)
            length      => elemR(:,er_Length)
        !%------------- ----------------------------------------------------

        where (volume(thisP)/length(thisP) .le. setting%ZeroValue%Area)
            depth(thisP) = setting%ZeroValue%Depth * 0.99d0
        end where



        ! print *, ' '
        ! print *, 'in geo_ZeroDepth_from_volume'
        ! print *, volume(54), length(54), setting%ZeroValue%Volume
        ! print *, volume(54)/length(54), setting%ZeroValue%Area
        ! print *, ' '

    end subroutine geo_ZeroDepth_from_volume
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_ZeroDepth_from_depth (thisP)    
        !%------------------------------------------------------------------
        !% Description:
        !% ensures that if depth <= zeroDepth the depth is reset to 99% 
        !% of zerodepth
        !%------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: depth(:)
        !%------------- -----------------------------------------------------
        !% Aliases
            depth       => elemR(:,er_Depth)
        !%------------- ----------------------------------------------------

        where (depth(thisP).le. setting%ZeroValue%Depth)
            depth(thisP) = setting%ZeroValue%Depth * 0.99d0
        end where

    end subroutine  geo_ZeroDepth_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_depth_from_head (epCol, Npack)
        !%------------------------------------------------------------------
        !% Description
        !% Computes depth from given head with 99% of Zerodepth as a minimum
        !% for a packed array of elements
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: epCol, Npack
            integer, pointer    :: thisP(:)
        !%------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return 
        !%------------------------------------------------------------------
        !% Aliases
            thisP => elemP(1:Npack,epCol)
        !%------------------------------------------------------------------

        elemR(thisP,er_Depth) = elemR(thisP,er_Head)                     &
             - (elemR(thisP,er_Zbottom) + elemR(thisP,er_SedimentDepth))
             
        where (elemR(thisP,er_Depth) .le. setting%ZeroValue%Depth)
            elemR(thisP,er_Depth) = 0.99d0* setting%ZeroValue%Depth
        endwhere

    end subroutine geo_depth_from_head
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function geo_depth_from_head_singular (Eidx)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes depth from a given head with 99% of Zerodepth as
        !% minimum for a single element
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: Eidx
        !%------------------------------------------------------------------

        geo_depth_from_head_singular                                                           &
            = max(                                                                             &
                elemR(Eidx,er_Head) - (elemR(Eidx,er_Zbottom) + elemR(Eidx,er_SedimentDepth)), &
                0.99d0 * setting%ZeroValue%Depth                                               &
            )

    end function geo_depth_from_head_singular
!%
!%==========================================================================
!%==========================================================================
!%
    ! real(8) function geo_hyddepth_from_depth_singular &
    !     (idx, indepth) result (outvalue)
    !     !% OBSOLETE 20220930 brh
    !     !%------------------------------------------------------------------
    !     !% Descriptions:
    !     !% computes the hydraulic depth for a given depth of a single element
    !     !%------------------------------------------------------------------
    !     !% Declarations
    !         real(8), intent(in)  :: indepth
    !         integer, intent(in)  :: idx
    !         real(8)              :: temp1, temp2
    !         character(64) :: subroutine_name = 'geo_hyddepth_from_depth_singular'
    !     !%------------------------------------------------------------------
    !     !%------------------------------------------------------------------
    !     select case (elemI(idx,ei_geometryType))
            
    !     !% --- OPEN CHANNELS -----------------------------   
    !     case (rectangular)
    !         outvalue = indepth
    !     case (trapezoidal)
    !         outvalue = trapezoidal_hyddepth_from_depth_singular (idx, indepth)
    !     case (triangular)
    !         outvalue = triangular_hyddepth_from_depth_singular (idx, indepth)
    !     case (parabolic)
    !         outvalue = parabolic_hyddepth_from_depth_singular (idx, indepth)
    !     case (power_function)
    !         print *, 'CODE ERROR: hyddepth for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
    !         print *, 'has not been implemented in ',trim(subroutine_name)
    !         call util_crashpoint(449734)
    !     case (irregular)
    !         !% --- get the area and topwidth, then compute the hydraulic depth
    !         temp1 = irregular_geometry_from_depth_singular (idx,tt_area,  indepth, setting%ZeroValue%Area)
    !         temp2 = irregular_geometry_from_depth_singular (idx,tt_width, indepth, setting%ZeroValue%TopWidth)
    !         outvalue = temp1 / temp2    

    !     !% --- CLOSED CONDUITS --------------------------------------    
    !     case (rectangular_closed)
    !         outvalue = rectangular_closed_hyddepth_from_depth_singular (idx, indepth)
    !     case (rect_round )
    !         outvalue = rect_round_hyddepth_from_depth_singular (idx,indepth)
    !     case (rect_triang)
    !         outvalue = rectangular_triangular_hyddepth_from_depth_singular (idx, indepth)
    !     case (circular )
    !         !% --- get the topwidth and use that to compute the hydraulic depth
    !         temp1    = circular_topwidth_from_depth_singular    (idx, indepth)
    !         outvalue = circular_hyddepth_from_topwidth_singular (idx,temp1,indepth)
    !     case (semi_circular)
    !         temp1    = semi_circular_topwidth_from_depth_singular    (idx, indepth)
    !         outvalue = semi_circular_hyddepth_from_topwidth_singular (idx,temp1,indepth)
    !     case (filled_circular)
    !         !% --- get the topwidth and use that to compute the hydraulic depth
    !         temp1    = filled_circular_topwidth_from_depth_singular    (idx, indepth)
    !         outvalue = filled_circular_hyddepth_from_topwidth_singular (idx,temp1,indepth)
    !     case (arch)
    !         !% --- get the topwidth and use that to compute the hydraulic depth
    !         temp1    = arch_topwidth_from_depth_singular    (idx, indepth)
    !         outvalue = arch_hyddepth_from_topwidth_singular (idx,temp1,indepth)
    !     case (basket_handle)
    !         !% --- get the topwidth and use that to compute the hydraulic depth
    !         temp1    = basket_handle_topwidth_from_depth_singular    (idx, indepth)
    !         outvalue = basket_handle_hyddepth_from_topwidth_singular (idx,temp1,indepth)
    !     case (catenary)
    !         temp1    = catenary_topwidth_from_depth_singular    (idx, indepth)
    !         outvalue = catenary_hyddepth_from_topwidth_singular (idx,temp1,indepth)
    !     case (eggshaped)
    !         temp1    = egg_shaped_topwidth_from_depth_singular    (idx, indepth)
    !         outvalue = egg_shaped_hyddepth_from_topwidth_singular (idx,temp1,indepth)
    !     case (gothic)
    !         temp1    = gothic_topwidth_from_depth_singular    (idx, indepth)
    !         outvalue = gothic_hyddepth_from_topwidth_singular (idx,temp1,indepth)
    !     case (horseshoe)
    !         !% --- get the topwidth and use that to compute the hydraulic depth
    !         temp1    = horse_shoe_topwidth_from_depth_singular    (idx, indepth)
    !         outvalue = horse_shoe_hyddepth_from_topwidth_singular (idx,temp1,indepth)
    !     case (horiz_ellipse)
    !         !% --- get the topwidth and use that to compute the hydraulic depth
    !         temp1    = horiz_ellipse_topwidth_from_depth_singular    (idx, indepth)
    !         outvalue = horiz_ellipse_hyddepth_from_topwidth_singular (idx,temp1,indepth)
    !     case (mod_basket)
    !         temp1    = mod_basket_topwidth_from_depth_singular    (idx, indepth)
    !         outvalue = mod_basket_hyddepth_from_topwidth_singular (idx,temp1,indepth)
    !     case (semi_elliptical)
    !         temp1    = semi_elliptical_topwidth_from_depth_singular    (idx, indepth)
    !         outvalue = semi_elliptical_hyddepth_from_topwidth_singular (idx,temp1,indepth)
    !     case (vert_ellipse)
    !         !% --- get the topwidth and use that to compute the hydraulic depth
    !         temp1    = vert_ellipse_topwidth_from_depth_singular    (idx, indepth)
    !         outvalue = vert_ellipse_hyddepth_from_topwidth_singular (idx,temp1,indepth)
        
    !     case (custom)
    !         print *, 'CODE ERROR: hyddepth for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
    !         print *, 'has not been implemented in ',trim(subroutine_name)
    !         call util_crashpoint(449734)
    !     case (force_main)
    !         print *, 'CODE ERROR: hyddepth for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
    !         print *, 'in ',trim(subroutine_name)   
    !         print *, 'This should never be reached as a force_main is not a valid geometryType'  
    !         call util_crashpoint(449734)
    !     case default
    !         print *, 'CODE ERROR: hyddepth for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
    !         print *, 'has not been implemented in ',trim(subroutine_name)
    !         call util_crashpoint(449734)
    !     end select
           
    ! end function geo_hyddepth_from_depth_singular
!%
!%==========================================================================

!%=========================================================================
!%
    ! subroutine geometry_table_initialize ()
    !     !%-----------------------------------------------------------------
    !     !% Description
    !     !% Creates element based geometry table
    !     !%-----------------------------------------------------------------
    !     !% Declarations:
    !         integer :: ii
    !         real(8), pointer :: table(:,:,:), fulldepth(:)
    !         !% local store of uniform distribution of depth
    !         real(8)          :: depth(Ncol2_GeometryTableR)
    !     !%-----------------------------------------------------------------
    !     !%-----------------------------------------------------------------
    !     !% Aliases
    !         table     =>  geometryTableR(:,:,:)
    !         eDepth    => elemR(:,er_Depth)
    !         fulldepth => elemR(:,er_FullDepth)
    !     !%-----------------------------------------------------------------
    !     !%-----------------------------------------------------------------

    !     !% --- cycle through elements (slow)
    !     do ii=1,N_elem(this_image())
    !         !% --- get a uniformly distributed depth array
    !         depth = real((/ (ii,ii=0,Ncol2_GeometryTableR-oneI) /),8)
    !         depth = fulldepth(ii) * depth / real(Ncol2_GeometryTableR-oneI,8)

    !         select case (elemI(ii,ei_geometryType))

    !         case (rectangular)
    !             !tdepth = eDepth(ii) !% temporary store
    !             !eDepth(ii) = depth(kk) !% set the local element depth
    !             !eDepth(ii) = tdepth !% restore the element depth
    !             do kk = 1,Ncol2_GeometryTableR
    !                 table(ii,kk,gtr_Area_from_Depth)      = rectangular_area_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Topwidth_from_Depth)  = rectangular_topwidth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Perimeter_from_Depth) = rectangular_perimeter_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydDepth_from_Depth)  = rectangular_hyddepth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydRadius_from_Depth) = rectangular_hydradius_from_depth_singular(ii,depth(kk))
    !             end do
                
    !         case (trapezoidal)
    !             do kk = 1,Ncol2_GeometryTableR
    !                 table(ii,kk,gtr_Area_from_Depth)      = trapezoidal_area_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Topwidth_from_Depth)  = trapezoidal_topwidth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Perimeter_from_Depth) = trapezoidal_perimeter_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydDepth_from_Depth)  = trapezoidal_hyddepth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydRadius_from_Depth) = trapezoidal_hydradius_from_depth_singular(ii,depth(kk))
    !             end do

    !         case (triangular)
    !             do kk = 1,Ncol2_GeometryTableR
    !                 table(ii,kk,gtr_Area_from_Depth)      = triangular_area_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Topwidth_from_Depth)  = triangular_topwidth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Perimeter_from_Depth) = triangular_perimeter_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydDepth_from_Depth)  = triangular_hyddepth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydRadius_from_Depth) = triangular_hydradius_from_depth_singular(ii,depth(kk))
    !             end do

    !         case (parabolic)
    !             do kk = 1,Ncol2_GeometryTableR
    !                 table(ii,kk,gtr_Area_from_Depth)      = parabolic_area_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Topwidth_from_Depth)  = parabolic_topwidth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Perimeter_from_Depth) = parabolic_perimeter_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydDepth_from_Depth)  = parabolic_hyddepth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydRadius_from_Depth) = parabolic_hydradius_from_depth_singular(ii,depth(kk))
    !             end do

    !         case (power_function)
    !             print *, 'CODE ERROR: POWER FUNCTION NOT COMPLETE'

    !         case (irregular)
    !             do kk = 1,Ncol2_GeometryTableR
    !                 table(ii,kk,gtr_Area_from_Depth)      = irregular_geometry_from_depth_singular(ii,tt_area,     depth(kk),setting%ZeroValue%Area)
    !                 table(ii,kk,gtr_Topwidth_from_Depth)  = irregular_geometry_from_depth_singular(ii,tt_width,    depth(kk),setting%ZeroValue%Topwidth)
    !                 table(ii,kk,gtr_HydRadius_from_Depth) = irregular_geometry_from_depth_singular(ii,tt_hydradius,depth(kk),setting%ZeroValue%Depth)
    !             end do
    !             table(ii,:,gtr_Perimeter_from_Depth) = table(ii,:,gtr_Area_from_Depth) / table(ii,:,gtr_HydRadius_from_Depth)
    !             table(ii,:,gtr_HydDepth_from_Depth)  = table(ii,:,gtr_Area_from_Depth) / table(ii,:,gtr_Topwidth_from_Depth)

    !         case (rectangular_closed)
    !             do kk = 1,Ncol2_GeometryTableR
    !                 table(ii,kk,gtr_Area_from_Depth)      = rectangular_closed_area_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Perimeter_from_Depth) = rectangular_closed_perimeter_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Topwidth_from_Depth)  = rectangular_closed_topwidth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydDepth_from_Depth)  = rectangular_closed_hyddepth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydRadius_from_Depth) = rectangular_closed_hydradius_from_depth_singular(ii,depth(kk))
    !             end do

    !         case (rect_round)
    !             do kk = 1,Ncol2_GeometryTableR
    !                 table(ii,kk,gtr_Area_from_Depth)      = rect_round_area_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Topwidth_from_Depth)  = rect_round_topwidth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Perimeter_from_Depth) = rect_round_perimeter_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydDepth_from_Depth)  = rect_round_hyddepth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydRadius_from_Depth) = rect_round_hydradius_from_depth_singular(ii,depth(kk))
    !             end do
                
    !         case (rect_triang)
    !             do kk = 1,Ncol2_GeometryTableR
    !                 table(ii,kk,gtr_Area_from_Depth)      = rectangular_triangular_area_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Topwidth_from_Depth)  = rectangular_triangular_topwidth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Perimeter_from_Depth) = rectangular_triangular_perimeter_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydDepth_from_Depth)  = rectangular_triangular_hyddepth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydRadius_from_Depth) = rectangular_triangular_hydradius_from_depth_singular(ii,depth(kk))
    !             end do

    !         case (circular)
    !             do kk = 1,Ncol2_GeometryTableR
    !                 table(ii,kk,gtr_Area_from_Depth)      = circular_area_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Topwidth_from_Depth)  = circular_topwidth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Perimeter_from_Depth) = circular_perimeter_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydDepth_from_Depth)  = circular_hyddepth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydRadius_from_Depth) = circular_hydradius_from_depth_singular(ii,depth(kk))
    !             end do

    !         case (semi_circular)
    !             do kk = 1,Ncol2_GeometryTableR
    !                 table(ii,kk,gtr_Area_from_Depth)      = semi_circular_area_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Topwidth_from_Depth)  = semi_circular_topwidth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Perimeter_from_Depth) = semi_circular_perimeter_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydDepth_from_Depth)  = semi_circular_hyddepth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydRadius_from_Depth) = semi_circular_hydradius_from_depth_singular(ii,depth(kk))
    !             end do

    !         case (filled_circular)
    !             do kk = 1,Ncol2_GeometryTableR
    !                 table(ii,kk,gtr_Area_from_Depth)      = filled_circular_area_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Topwidth_from_Depth)  = filled_circular_topwidth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Perimeter_from_Depth) = filled_circular_perimeter_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydDepth_from_Depth)  = filled_circular_hyddepth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydRadius_from_Depth) = filled_circular_hydradius_from_depth_singular(ii,depth(kk))
    !             end do

    !         case (arch)
    !             do kk = 1,Ncol2_GeometryTableR
    !                 table(ii,kk,gtr_Area_from_Depth)      = arch_area_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Topwidth_from_Depth)  = arch_topwidth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Perimeter_from_Depth) = arch_perimeter_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydDepth_from_Depth)  = arch_hyddepth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydRadius_from_Depth) = arch_hydradius_from_depth_singular(ii,depth(kk))
    !             end do

    !         case (basket_handle)
    !             do kk = 1,Ncol2_GeometryTableR
    !                 table(ii,kk,gtr_Area_from_Depth)      = basket_handle_area_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Topwidth_from_Depth)  = basket_handle_topwidth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Perimeter_from_Depth) = basket_handle_perimeter_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydDepth_from_Depth)  = basket_handle_hyddepth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydRadius_from_Depth) = basket_handle_hydradius_from_depth_singular(ii,depth(kk))
    !             end do

    !         case (catenary)
    !             do kk = 1,Ncol2_GeometryTableR
    !                 table(ii,kk,gtr_Area_from_Depth)      = catenary_area_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Topwidth_from_Depth)  = catenary_topwidth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Perimeter_from_Depth) = catenary_perimeter_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydDepth_from_Depth)  = catenary_hyddepth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydRadius_from_Depth) = catenary_hydradius_from_depth_singular(ii,depth(kk))
    !             end do

    !         case (eggshaped)
    !             do kk = 1,Ncol2_GeometryTableR
    !                 table(ii,kk,gtr_Area_from_Depth)      = egg_shaped_area_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Topwidth_from_Depth)  = egg_shaped_topwidth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Perimeter_from_Depth) = egg_shaped_perimeter_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydDepth_from_Depth)  = egg_shaped_hyddepth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydRadius_from_Depth) = egg_shaped_hydradius_from_depth_singular(ii,depth(kk))
    !             end do

    !         case (gothic)
    !             do kk = 1,Ncol2_GeometryTableR
    !                 table(ii,kk,gtr_Area_from_Depth)      = gothic_area_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Topwidth_from_Depth)  = gothic_topwidth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Perimeter_from_Depth) = gothic_perimeter_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydDepth_from_Depth)  = gothic_hyddepth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydRadius_from_Depth) = gothic_hydradius_from_depth_singular(ii,depth(kk))
    !             end do

    !         case (horseshoe)
    !             do kk = 1,Ncol2_GeometryTableR
    !                 table(ii,kk,gtr_Area_from_Depth)      = horse_shoe_area_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Topwidth_from_Depth)  = horse_shoe_topwidth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Perimeter_from_Depth) = horse_shoe_perimeter_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydDepth_from_Depth)  = horse_shoe_hyddepth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydRadius_from_Depth) = horse_shoe_hydradius_from_depth_singular(ii,depth(kk))
    !             end do

    !         case (horiz_ellipse)
    !             do kk = 1,Ncol2_GeometryTableR
    !                 table(ii,kk,gtr_Area_from_Depth)      = horiz_ellipse_area_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Topwidth_from_Depth)  = horiz_ellipse_topwidth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Perimeter_from_Depth) = horiz_ellipse_perimeter_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydDepth_from_Depth)  = horiz_ellipse_hyddepth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydRadius_from_Depth) = horiz_ellipse_hydradius_from_depth_singular(ii,depth(kk))
    !             end do

    !         case (mod_basket)
    !             do kk = 1,Ncol2_GeometryTableR
    !                 table(ii,kk,gtr_Area_from_Depth)      = mod_basket_area_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Topwidth_from_Depth)  = mod_basket_topwidth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_Perimeter_from_Depth) = mod_basket_perimeter_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydDepth_from_Depth)  = mod_basket_hyddepth_from_depth_singular(ii,depth(kk))
    !                 table(ii,kk,gtr_HydRadius_from_Depth) = mod_basket_hydradius_from_depth_singular(ii,depth(kk))
    !             end do

    !         case (semi_elliptical)


    !         case (vertical_ellipse)
    !         case (custom)
    !         case (force_main)
    !         end select
    !     end do


    ! end subroutine geometry_table_initialize
!%
!%=========================================================================
    !%        
!% 
!%==========================================================================
!%==========================================================================
!%
   ! subroutine geometry_toplevel()
        ! !%------------------------------------------------------------------
        ! !% Description:
        ! !% Input whichTM is one of ETM, AC, or ALLtm
        ! !% This should never be called for diagnostic arrays
        ! !% Note that the elemPGx arrays contain only time-marched elements so they
        ! !% will only handle CC and JM elements as the JB elements are not time-marched.
        ! !%------------------------------------------------------------------
        ! !% Declarations
        !     integer, intent(in) :: whichTM
        !     integer, pointer :: elemPGx(:,:), npack_elemPGx(:), col_elemPGx(:)
        !     !integer, pointer :: thisColP_surcharged, thisColP_NonSurcharged, 
        !     integer, pointer :: thisColP_all_TM, thisColP_Open_CC
        !     integer, pointer :: thisColP_JM, thisColP_JB, thisColP_CC
        !     integer, pointer :: thisColP_Closed_CC, thisColP_Closed_JB
        !     integer, pointer :: Npack, thisP(:)
        !     !integer, pointer ::  thisColP_Closed_JM
        !     logical :: isreset
        !     real(8), pointer :: depth(:), volume(:), area(:), topwidth(:)
        !     integer, allocatable :: tempP(:) !% debugging
        !     character(64) :: subroutine_name = 'geometry_toplevel'
        ! !%------------------------------------------------------------------
        ! !% Preliminaries
        !     !!if (crashYN) return
        !     if (setting%Debug%File%geometry) &
        !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        ! !%------------------------------------------------------------------
        ! !% Aliases
        !     area   => elemR(:,er_Area)
        !     depth  => elemR(:,er_Depth)
        !     topwidth => elemR(:,er_Topwidth)
        !     volume => elemR(:,er_Volume)
        ! !% set the packed geometry element array (elemPG) to use and columns of the
        ! !% packed elemP to use
        !     select case (whichTM)
        !         ! case (ALLtm)
        !         !     elemPGx                => elemPGalltm(:,:) 
        !         !     npack_elemPGx          => npack_elemPGalltm(:)
        !         !     col_elemPGx            => col_elemPGalltm(:)
        !         !     thisColP_CC            => col_elemP(ep_CC_ALLtm)
        !         !     thisColP_JM            => col_elemP(ep_JM_ALLtm)
        !         !     thisColP_JB            => col_elemP(ep_JB_ALLtm)
        !         !     !thisColP_surcharged    => col_elemP(ep_ALLtmSurcharged)
        !         !     !thisColP_NonSurcharged => col_elemP(ep_ALLtm_NonSurcharged)
        !         !     thisColP_all           => col_elemP(ep_ALLtm)
        !         !     thisColP_Open_CC       => col_elemP(ep_CC_Open_Elements)
        !         !     thisColP_Closed_CC     => col_elemP(ep_CC_Closed_Elements)
        !         !     thisColP_Closed_JB     => col_elemP(ep_JB_Closed_Elements)
        !         !     !thisColP_Closed_JM     => col_elemP(ep_JM_Closed_Elements)
        !         case (ETM)
        !             elemPGx                => elemPGetm(:,:)
        !             npack_elemPGx          => npack_elemPGetm(:)
        !             col_elemPGx            => col_elemPGetm(:)
        !             thisColP_CC            => col_elemP(ep_CC_ETM)
        !             thisColP_JM            => col_elemP(ep_JM_ETM)
        !             !thisColP_JB            => col_elemP(ep_JB_ETM)
        !             !thisColP_surcharged    => col_elemP(ep_PSsurcharged)
        !             !thisColP_NonSurcharged => col_elemP(ep_ETM_PSnonSurcharged)
        !             thisColP_all_TM        => col_elemP(ep_ETM)
        !             thisColP_Open_CC       => col_elemP(ep_CC_Open_Elements)
        !             thisColP_Closed_CC     => col_elemP(ep_CC_Closed_Elements)
        !             thisColP_Closed_JB     => col_elemP(ep_JB_Closed_Elements)
        !             !thisColP_Closed_JM     => col_elemP(ep_JM_Closed_Elements)
        !         ! case (AC)
        !         !     elemPGx                => elemPGac(:,:)
        !         !     npack_elemPGx          => npack_elemPGac(:)
        !         !     col_elemPGx            => col_elemPGac(:)
        !         !     thisColP_CC            => col_elemP(ep_CC_AC)
        !         !     thisColP_JM            => col_elemP(ep_JM_AC)
        !         !     thisColP_JB            => col_elemP(ep_JB_AC)
        !         !     !thisColP_surcharged    => col_elemP(ep_ACsurcharged)
        !         !     !thisColP_NonSurcharged => col_elemP(ep_AC_ACnonSurcharged)
        !         !     thisColP_all           => col_elemP(ep_AC)
        !         !     thisColP_Open_CC       => col_elemP(ep_CC_Open_Elements)
        !         !     thisColP_Closed_CC     => col_elemP(ep_CC_Closed_Elements)
        !         !     thisColP_Closed_JB     => col_elemP(ep_JB_Closed_Elements)
        !         !     !thisColP_Closed_JM     => col_elemP(ep_JM_Closed_Elements)
        !         case default
        !             print *, 'CODE ERROR: time march type unknown for # ', whichTM
        !             print *, 'which has key ',trim(reverseKey(whichTM))
        !             call util_crashpoint(7389)
        !             !return
        !             !stop 7389
        !     end select
        !     call util_crashstop(49872)
        ! !%--------------------------------------------------------------------
        !     ! ! ! call util_utestLprint ('in geometry at top==========================================')  

        ! !% STATUS: at this point we know volume and velocity on all elements
        ! !% from RK2 solution

        ! !% --- PONDING
        ! !%     adjust time-march volume for ponding inflow
        ! !%     This affects JM that have previous ponding but now have volumes
        ! !%     below the full volume
        ! if (setting%SWMMinput%AllowPonding) then
        !     call geo_ponding_inflow (thisColP_JM)   
        ! end if
        !     ! ! ! call util_utestLprint ('in geometry after ponding inflow') 
           
        ! !% --- PREISSMAN SLOT    
        ! !%     also adds to ponded volume or computes overflow volume for JM (only)
        ! !%     where slot depth + invert height exceeds maximum surcharge height
        ! !%     The Element Volume in JM is adjusted for any overflow or ponding
        ! !%     (but not for the slot)
        ! call slot_toplevel (thisColP_Closed_CC, thisColP_JM)

        !     ! ! ! call util_utestLprint ('in geometry after slot toplevel') 

        !     !stop 29387488

        ! !% STATUS: The Preissmann Slot values have been assigned for all CC and JM
        ! !% Overflow and Ponding have been assigned for JM (only). 
        ! !% JM element volumes have been adjusted for overflow or ponding, but
        ! !% not for slot volume. So both JM and CC have volume > fullvolume
        ! !% in surcharged elements.

        ! !% --- assign all geometry for surcharged elements CC, JM
        ! !%     Note: not used in Preissmann Slot (only AC)
        ! !%     DO NOT DELETE. HOLD THIS FOR LATER USE 20220909brh
        ! ! if ((whichTM .eq. ALLtm) .or. (whichTM .eq. AC)) then
        ! !     call geo_ACsurcharged (thisColP_surcharged)
        ! ! end if

        ! !% --- ZERO VOLUMES CC JM
        ! !%     reset all zero or near-zero volumes in all CC, JM
        ! !call adjust_limit_by_zerovalues (er_Volume, setting%ZeroValue%Volume, thisColP_NonSurcharged, .true.)
        ! call adjust_limit_by_zerovalues &
        !     (er_Volume, setting%ZeroValue%Volume, thisColP_all_TM, .true.)

        !     ! ! ! call util_utestLprint ('in geometry after limit_by_zerovalues (volume)') 

        ! !% --- DEPTH
        ! !%     compute the depth on all elements of CC JM based on geometry.
        ! !%     If surcharged, this call returns the full depth of a closed conduit 
        ! !%     without adding Preissmann Slot depth.
        ! call geo_depth_from_volume_by_type (elemPGx, npack_elemPGx, col_elemPGx)

        ! ! print *, 'after geo_depth_from_volume '
        ! ! print *, elemR(50,er_Depth), elemR(51,er_Depth)

        !     ! ! ! call util_utestLprint ('in geometry after depth_from_volume') 

        ! !% --- ZERO DEPTH CC JM -- now done in geo_depth_from_volume_by_type 20230113
        ! !%     reset all zero or near-zero depths in aa CC and JM
        ! !call adjust_limit_by_zerovalues (er_Depth, setting%ZeroValue%Depth, thisColP_NonSurcharged, .false.)
        ! call adjust_limit_by_zerovalues &
        !     (er_Depth, setting%ZeroValue%Depth, thisColP_all_TM, .false.)

        !     ! ! ! call util_utestLprint ('in geometry after limit_by_zerovalues (depth)') 

        ! !% --- PIEZOMETRIC HEAD
        ! !%     compute the head on all elements of CC and JM
        ! !%     This sets head consistent with depth computed in geo_depth_from_volume
        ! !%     Head is strictly limited to the max depth + zbottom so it does not
        ! !%     include surcharge effects     
        ! Npack     => npack_elemP(thisColP_all_TM)
        ! if (Npack > 0) then
        !     thisP => elemP(1:Npack,thisColP_all_TM)
        !     elemR(thisP,er_Head) = llgeo_head_from_depth_pure (thisP, depth(thisP))
        !     !elemR(thisP,er_PressureHead) = llgeo_head_from_depth_pure (thisP)
        !     !call geo_head_from_depth (thisColP_all_TM)
        ! end if
 
        ! ! print *, 'after geo_head_from_depth '
        ! ! print *, elemR(50,er_Head), elemR(51,er_Head)
        !     ! ! ! call util_utestLprint ('in geometry after head_from_depth')

        ! !% --- OPEN CHANNEL OVERFLOW
        ! !%     Compute the overflow lost for CC open channels above
        ! !%     their maximum volume (no ponding allowed from open CC). 
        ! !%     Note that overflow or ponding for JM elements is handled 
        ! !%     in slot_JM_ETM.
        ! !%     Note, this is NOT standard in EPA-SWMM
        ! if (setting%Discretization%AllowChannelOverflowTF) then
        !     call geo_overflow_openchannels (thisColP_Open_CC)
        ! end if

        !     ! ! ! call util_utestLprint ('in geometry after overflow_openchannels')

        ! !% --- PREISSMAN SLOT VOLUME LIMIT CLOSED CONDUIT CC JM
        ! !%     limit the volume in closed element (CC, JM) to the full volume
        ! !%     Note the excess volume has already been stored in the Preissman Slot
        ! call geo_volumelimit_closed (thisColP_Closed_CC)
        ! call geo_volumelimit_closed (thisColP_JM)

        !     ! ! ! call util_utestLprint ('in geometry after volumelimit_closed')

        ! !% REMOVED 20220909 brh
        ! !% --- limit volume for incipient surcharge. This is done after depth is computed
        ! !%     so that the "depth" algorithm can include depths greater than fulldepth
        ! !%     as a way to handle head for incipient surcharge.
        ! !call geo_limit_incipient_surcharge (er_Volume, er_FullVolume, thisColP_NonSurcharged,.true.) !% 20220124brh
        !     ! ! ! call util_utestLprint ('in geometry before geo_limit_incipient_surcharge (Depth)')  
        ! !% --- limit depth for surcharged on CC. This is done after head is computed
        ! !%     so that the depth algorithm can include depths greater than fulldepth where the 
        ! !%     geometry algorithm does not enforrce full depth
        ! !call geo_limit_incipient_surcharge (er_Depth, er_FullDepth, thisColP_NonSurcharged,.false.) 
        ! !@call geo_limit_incipient_surcharge (er_Depth, er_FullDepth, thisColP_all,.false.) 
        ! !% END REMOVE 20220909

        ! !% STATUS: At this point, the depths, heads and volumes of all CC, JM elements are
        ! !% at or below their full value.  For CC closed conduits and all JM the
        ! !% surcharged volume is stored in the er_SlotVolume and the surcharged extra
        ! !% depth is stored in the er_SlotDepth 

        ! !% --- PREISSMAN SLOT HEAD ADD IN JM 
        ! !%     adjust JM head to include Preissmann Slot Depth and ponding
        ! !%     This is needed before JB are computed
        ! call slot_JM_head_PSadd (thisColP_JM)

        !     ! ! ! call util_utestLprint ('in geometry after JM_head_PSadd') 
           
        ! !% --- JB VALUES
        ! !%    assign the non-volume geometry on junction branches JB based on JM head
        ! !%    Values limited by full volume. Volume assigned is area * length
        ! call geo_assign_JB_from_head (whichTM, thisColP_JM)

        !     ! ! ! call util_utestLprint ('in geometry after assign_JB') 

        ! !% --- JB CLOSED CONDUIT VOLUME LIMIT
        ! !%     further limiting the JB volume by full is probably not needed,
        ! !%     but might be useful if there's a numerical precision issues
        ! !%     with JB volume assigned by area * length.
        ! call geo_volumelimit_closed (thisColP_Closed_JB)

        !     ! ! ! call util_utestLprint ('in geometry after volumelimit_closed') 

        ! !% --- PREISSMANN SLOT HEAD REMOVE IN JM
        ! !%     we need to remove the PS and ponding from the JM cells so that we can easily
        ! !%     compute other geometry without full JM causing problems
        ! call slot_JM_head_PSremove (thisColP_JM)

        !     ! ! ! call util_utestLprint ('in geometry after JM_head_PSremove')  

        ! !% STATUS: at this point we have all geometry on CC, JM, JB that is
        ! !% limited by the full volume values. The CC and JM have slot values stored
        ! !% but no slot values have been computed for JB
        
        ! !% --- CROSS-SECTIONAL AREA
        ! !%    compute area from volume for CC, JM
        ! !%     For closed conduits this is based on the volume limited by full volume.
        ! !%     For open channels the volume limit depends on if AllowChanneOverflowTF is false.
        ! !%     Note that JB areas are already assigned in geo_assign_JB_from_head()
        ! Npack     => npack_elemP(thisColP_all_TM)
        ! if (Npack > 0) then
        !     thisP => elemP(1:Npack,thisColP_all_TM)
        !     elemR(thisP,er_Area) = llgeo_area_from_volume_pure (thisP, volume(thisP))
        !     elemR(thisP,er_Area) = max(elemR(thisP,er_Area),setting%ZeroValue%Area)
        !     !call geo_area_from_volume (thisColP_all_TM)
        ! end if

        !     ! ! ! call util_utestLprint ('in geometry after area_from_volume') 

        ! ! !% --- ZERO AREA CC JM
        ! ! !%     reset all zero or near-zero areas in CC and JM
        ! ! call adjust_limit_by_zerovalues &
        ! !      (er_Area, setting%ZeroValue%Area, thisColP_all_TM, .false.)

        !     ! ! ! call util_utestLprint ('in geometry after adjust_limit_by_zeroValues area')   

        ! !% --- TOPWIDTH CC
        ! !%     compute topwidth from depth for all CC
        ! !%     Note: Topwidth for JM is undefined in this subroutine
        ! !%     Note: volume is limited to full depth UNLESS AllowChannelOverflowTF is false
        ! call geo_topwidth_from_depth_by_type_CC (elemPGx, npack_elemPGx, col_elemPGx)

        !     ! ! ! call util_utestLprint ('in geometry after topwidth_from_depth') 

        ! ! !% --- ZERO TOPWIDTH CC
        ! ! !%     reset all zero or near-zero topwidth in CC 
        ! ! !%     but do not change the eYN(:,eYN_isZeroDepth) mask
        ! ! call adjust_limit_by_zerovalues &
        ! !      (er_Topwidth, setting%ZeroValue%Topwidth, thisColP_CC, .false.)

        !     ! ! ! call util_utestLprint ('in geometry after adjust_limit_by_zerovalues topwidth') 

        ! !% --- PERIMETER AND HYDRAULIC RADIUS CC
        ! !%     compute hydraulic radius and perimeter
        ! !%     note these two are done together because for analytical cross-sections
        ! !%     we have equations for perimeter, whereas lookup cross-sections
        ! !%     have tables for hydraulic radius.
        ! call geo_perimeter_and_hydradius_from_depth_by_type_CC (elemPGx, npack_elemPGx, col_elemPGx)  

        ! ! % --- compute perimeter from maximum depth for all CC
        ! ! %     Note: perimeter for JM is undefined in this subroutine
        ! !OBSOLETE  call geo_perimeter_from_depth (elemPGx, npack_elemPGx, col_elemPGx)

        !     ! ! ! call util_utestLprint ('in geometry after perimeter from depth') 

        ! !% --- compute hyddepth
        ! !call geo_hyddepth_from_depth_or_topwidth (elemPGx, npack_elemPGx, col_elemPGx)
        ! !% 20220930 replace with unified call
        ! ! Npack     => npack_elemP(thisColP_CC)
        ! ! if (Npack > 0) then
        ! !     thisP => elemP(1:Npack,thisColP_CC)
        ! !     elemR(thisP,er_HydDepth) = llgeo_hyddepth_from_area_and_topwidth_pure &
        ! !                                 (thisP, area(thisP), topwidth(thisP))
        ! ! end if
        ! ! OBSOLETE call geo_hyddepth_from_area_and_topwidth (thisColP_CC)

        ! !    util_utest_CLprint ('in geometry after hyddepth_from_depth')

        ! !% --- compute hydradius
        ! !%     Note: cannot be used for JM unless perimeter is defined prior.
        ! !OBSOLETE call geo_hydradius_from_area_perimeter (thisColP_CC)

        ! !    util_utest_CLprint ('in geometry after hydradius_from_area_perimeter') 

        ! !% --- the modified hydraulic depth "ell" is used for 
        ! !%     for Froude number computations on all CC elements
        ! !%     Note: ell for JM is undefined in this subroutine

        ! call geo_elldepth_from_head_CC (thisColP_CC)

        ! !% --- compute pressure head from the modified hydraulic depth
        ! ! Npack     => npack_elemP(thisColP_CC)
        ! ! if (Npack > 0) then
        ! !     thisP => elemP(1:Npack,thisColP_CC)
        ! !     !elemR(thisP,er_Pressure_Head) = llgeo_pressure_head_from_hyddepth_pure (thisP)
        ! !     !call geo_pressure_head_from_hyddepth (thisColP_CC)
        ! ! end if

        !     ! ! ! call util_utestLprint ('in geometry after ell_from_head') 

        ! !% --- make adjustments for slots on closed elements only
        ! !%     These add slot values to volume, depth, head
        ! call slot_CC_adjustments (thisColP_Closed_CC)

        !     ! ! ! call util_utestLprint ('in geometry after slot_CC_adustments') 

        ! call slot_JM_adjustments (thisColP_JM)

        !     ! ! ! call util_utestLprint ('in geometry after slot_JM_adjustments') 

        ! !% --- compute the SlotDepth, SlotArea, SlotVolume and update
        ! !%     elem volume and depth for JB. Note elem head on JB either with or
        ! !%     without surcharge is assigned in geo_assign_JB_from_head
        ! call slot_JB_computation (thisColP_JM)
        
        !     ! ! ! call util_utestLprint ('in geometry after slot_JB_computation') 

        ! !% Set JM values that are not otherwise defined
        ! !% HydDepth, ell. Note that topwidth, hydradius, perimeter are undefined.
        ! call geo_JM_values ()

        !     ! ! ! call util_utestLprint ('in geometry after JM_values') 

        ! !% HOLD UNTIL AC RE-VISITED
        ! ! !% compute the dHdA that are only for AC nonsurcharged
        ! ! if (whichTM .ne. ETM) then
        ! !     call geo_dHdA (ep_AC_ACnonSurcharged)
        ! ! end if

        !     ! ! ! call util_utestLprint ('in geometry at end') 


        ! !stop 2397843

        ! !% --- check for crashpoint and stop here
        ! call util_crashstop(322983)

        ! if (setting%Debug%File%geometry) &
        ! write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
 !   end subroutine geometry_toplevel
!%
!%==========================================================================
    !%==========================================================================
!%
    ! subroutine geo_area_from_volume (thisColP)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% sets area = volume/length which is common to all nonsurcharged elements
    !     !% Note this assumes volume has been limited by surcharge and zero values
    !     !%------------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in) :: thisColP
    !         integer, pointer :: thisP(:), Npack
    !         real(8), pointer :: area(:), volume(:), length(:)

    !         character(64) :: subroutine_name = 'geo_area_from_volume'
    !     !%--------------------------------------------------------------------
    !     !% Preliminaries
    !         Npack  => npack_elemP(thisColP)
    !         if (Npack < 1) return
    !         if (setting%Debug%File%geometry) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%--------------------------------------------------------------------
    !     !% Aliases
    !         thisP  => elemP(1:Npack,thisColP)
    !         area   => elemR(:,er_Area)
    !         volume => elemR(:,er_Volume)
    !         length => elemR(:,er_Length)
    !     !%--------------------------------------------------------------------
     
    !     !% --- note, this could cause issues if length = 0, which should not happen
    !     area(thisP) = volume(thisP) / length(thisP)

    !     !%--------------------------------------------------------------------
    !         if (setting%Debug%File%geometry) &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine geo_area_from_volume
!%
!%==========================================================================
    !%==========================================================================
!%
    ! subroutine geo_depth_from_volume_by_type (elemPGx, npack_elemPGx, col_elemPGx)

        !% OBSOLETE WITH JUNCTION IMPLICIT. REPLACED BY ..._CC and ..._JM

        !%------------------------------------------------------------------
        !% Description:
        !% This solves nonsurcharged CCJMJB elements because of PGx arrays
        !% The elemPGx determines whether this is ALLtm, ETM or AC elements
        !%------------------------------------------------------------------
!         !% Declarations:
!     integer, target, intent(in) :: elemPGx(:,:), npack_elemPGx(:), col_elemPGx(:)
!     integer, pointer :: Npack, thisCol
!     character(64) :: subroutine_name = 'geo_depth_from_volume_by_type'
! !%-------------------------------------------------------------------
! !% Preliminaries
!     !!if (crashYN) return
!     if (setting%Debug%File%geometry) &
!         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
!%-------------------------------------------------------------------    
!% cycle through different geometries  

        ! print *, ' '
        ! print *, 'at start of geo_depth_from_volume_by_type'
        ! print *, 54, elemR(54,er_Depth), elemR(54,er_Volume)
        ! print *, ' '

        
!% --- open channels ------------------------------------------

! !% --- IRREGULAR
! thisCol => col_elemPGx(epg_CC_irregular)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call irregular_depth_from_volume (elemPGx, Npack, thisCol)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)
! end if    
        
! !% -- POWER FUNCTION
! thisCol => col_elemPGx(epg_CC_power_function)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     print *, 'POWER FUNCTION CROSS-SECTION NOT COMPLETE'
!     call util_crashpoint(5559872)
!     !call powerfunction_depth_from_volume (elemPGx, Npack, thisCol)
!     !call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)
! end if

! !% -- PARABOLIC
! thisCol => col_elemPGx(epg_CC_parabolic)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call parabolic_depth_from_volume (elemPGx, Npack, thisCol)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)
! end if
        
! !% --- RECTANGULAR CHANNEL
! thisCol => col_elemPGx(epg_CC_rectangular)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!         ! print *, 'rect here Volume/length',elemR(15,er_Volume)/elemR(15,er_Length)
!         ! print *,  'zero area           ',setting%ZeroValue%Area
!     call rectangular_depth_from_volume (elemPGx, Npack, thisCol)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)
! end if    

! !% --- TRAPEZOIDAL CHANNEL
! thisCol => col_elemPGx(epg_CC_trapezoidal)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call trapezoidal_depth_from_volume (elemPGx, Npack, thisCol)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)
! end if

! !% --- TRIANGULAR CHANNEL
! thisCol => col_elemPGx(epg_CC_triangular)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call triangular_depth_from_volume (elemPGx, Npack, thisCol)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)
! end if

! !% --- CLOSED CONDUITS  ---------------------------------------

! !% --  ARCH CONDUIT
! thisCol => col_elemPGx(epg_CC_arch)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call llgeo_tabular_depth_from_volume   &
!         (elemPGx, Npack, thisCol, YArch)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)
! end if

! !% --  BASKET_HANDLE
! thisCol => col_elemPGx(epg_CC_basket_handle)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call llgeo_tabular_depth_from_volume           &
!         (elemPGx, Npack, thisCol, YBasketHandle)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)   
! end if

! !% --  CATENARY
! thisCol => col_elemPGx(epg_CC_catenary)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call llgeo_tabular_depth_from_volume           &
!         (elemPGx, Npack, thisCol, YCatenary)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)    
! end if

! !% --- CIRCULAR CONDUIT
! thisCol => col_elemPGx(epg_CC_circular)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call circular_depth_from_volume (elemPGx, Npack, thisCol)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)
! end if

! !% --  EGG_SHAPED
! thisCol => col_elemPGx(epg_CC_egg_shaped)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call llgeo_tabular_depth_from_volume           &
!         (elemPGx, Npack, thisCol, YEgg)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)
! end if

! !% --- FILLED CIRCULAR CONDUIT
! thisCol => col_elemPGx(epg_CC_filled_circular)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call filled_circular_depth_from_volume (elemPGx, Npack, thisCol)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)
! end if

! !% --  GOTHIC
! thisCol => col_elemPGx(epg_CC_gothic)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call llgeo_tabular_depth_from_volume           &
!         (elemPGx, Npack, thisCol, YGothic)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)    
! end if

! !% --  HORIZONTAL ELLIPSE
! thisCol => col_elemPGx(epg_CC_horiz_ellipse)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call llgeo_tabular_depth_from_volume           &
!         (elemPGx, Npack, thisCol, YHorizEllip)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)    
! end if

! !% --  HORSESHOE
! thisCol => col_elemPGx(epg_CC_horse_shoe)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call llgeo_tabular_depth_from_volume           &
!         (elemPGx, Npack, thisCol, YHorseShoe)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)    
! end if

! !% --  MODIFIED BASKET HANDLE
! thisCol => col_elemPGx(epg_CC_mod_basket)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call mod_basket_depth_from_volume (elemPGx, Npack, thisCol)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)
! end if

! !% --- RECTANGULAR CLOSED CONDUIT
! thisCol => col_elemPGx(epg_CC_rectangular_closed)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call rectangular_closed_depth_from_volume (elemPGx, Npack, thisCol)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)
! end if        

! !% --- RECTANGULAR ROUND CONDUIT
! thisCol => col_elemPGx(epg_CC_rectangular_round)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call rect_round_depth_from_volume (elemPGx, Npack, thisCol)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)
! end if

! !% -- RECTANGULAR TRIANGULAR
! thisCol => col_elemPGx(epg_CC_rectangular_triangular)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call rectangular_triangular_depth_from_volume (elemPGx, Npack, thisCol)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)
! end if

! ! % --- SEMI-CIRCULAR CONDUIT
! thisCol => col_elemPGx(epg_CC_semi_circular)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call llgeo_tabular_depth_from_volume           &
!         (elemPGx, Npack, thisCol, YSemiCircular)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)    
! end if

! !% --  SEMI ELLIPTICAL
! thisCol => col_elemPGx(epg_CC_semi_elliptical)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call llgeo_tabular_depth_from_volume           &
!         (elemPGx, Npack, thisCol, YSemiEllip)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)    
! end if

! !% --  VERTICAL ELLIPSE
! thisCol => col_elemPGx(epg_CC_vert_ellipse)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call llgeo_tabular_depth_from_volume           &
!         (elemPGx, Npack, thisCol, YVertEllip)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)    
! end if


! !% --- JUNCTIONS ---------------------------------------------------- 

! !% JM with functional geometry
! thisCol => col_elemPGx(epg_JM_functionalStorage)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call storage_functional_depth_from_volume (elemPGx, Npack, thisCol)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)
! end if

! !% JM with tabular geometry
! thisCol => col_elemPGx(epg_JM_tabularStorage)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call storage_tabular_depth_from_volume (elemPGx, Npack, thisCol)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)
! end if

! !% JM with implied storage 
! thisCol => col_elemPGx(epg_JM_impliedStorage)
! Npack   => npack_elemPGx(thisCol)
! if (Npack > 0) then
!     call storage_implied_depth_from_volume (elemPGx, Npack, thisCol)
!     call geo_ZeroDepth_from_volume  (elemPGx, Npack, thisCol)
! end if

! !%-------------------------------------------------------------------
! !% Closing
!     if (setting%Debug%File%geometry) &
        ! write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
! end subroutine geo_depth_from_volume_by_type
!%
!%==========================================================================
!%==========================================================================
!%
! subroutine geo_head_from_depth (thisColP)
!     !%------------------------------------------------------------------
!     !% Description:
!     !% Computes head from depth for elements of CC, JM
!     !%------------------------------------------------------------------
!         integer, intent(in) :: thisColP
!         integer, pointer :: Npack, thisP(:)
!         real(8), pointer :: depth(:), fulldepth(:), head(:), Zbtm(:)
!         !real(8), pointer :: pressurehead(:)

!         character(64) :: subroutine_name = 'geo_head_from_depth'
!     !%------------------------------------------------------------------
!     !% Preliminaries
!         Npack     => npack_elemP(thisColP)
!         if (Npack < 1) return
!     !%------------------------------------------------------------------
!     !% Aliases    
!         thisP     => elemP(1:Npack,thisColP)
!         depth     => elemR(:,er_Depth)
!         fulldepth => elemR(:,er_FullDepth)
!         head      => elemR(:,er_Head)
!         !pressurehead => elemR(:,er_Pressure_Head)
!         Zbtm      => elemR(:,er_Zbottom)
!     !%------------------------------------------------------------------
!     !%

!     head(thisP) = depth(thisP) + Zbtm(thisP)
!     !% set the pressure head to piezometric head
!     !% this will be fixed later for CC and JB elements
!     !pressurehead(thisP) = head(thisP)

!     if (setting%Debug%File%geometry) &
!     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
! end subroutine geo_head_from_depth
!%
!%==========================================================================
!% END OF MODULE
!%=========================================================================
end module geometry