module geometry
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Toplevel methods for geometry computation and
    !% procedures needed during initialization and time marching
    !%
    !%==========================================================================

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
    ! use arch_conduit
    ! use basket_handle_conduit
    ! use catenary_conduit
    use circular_conduit
    ! use egg_shaped_conduit
    use filled_circular_conduit
    ! use gothic_conduit
    ! use horiz_ellipse_conduit
    ! use horse_shoe_conduit
    use mod_basket_conduit
    use rectangular_round_conduit
    use rectangular_triangular_conduit
    ! use semi_circular_conduit
    ! use semi_elliptical_conduit
    ! use vert_ellipse_conduit
    use storage_geometry
   
    use xsect_tables
    use adjust
    use utility_profiler
    use utility_crash

    ! use utility_unit_testing, only: util_utest_CLprint

    implicit none

    private

    public :: geometry_toplevel_CC
    public :: geo_assign_JB_from_head
    public :: geo_common_initialize
    public :: geo_sectionfactor_from_depth_singular
    public :: geo_Qcritical_from_depth_singular
    public :: geo_critical_value_singular
    public :: geo_normaldepth_singular
    public :: geo_topwidth_from_depth_singular
    public :: geo_area_from_depth_singular
    public :: geo_depth_from_volume_by_type_allCC
    public :: geo_depth_from_volume_by_element_CC

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
        !% Preliminary
            if (npackP < 1) return
        !%------------------------------------------------------------------

            ! if (.not. isSingularYN) call util_utest_CLprint('       111 geo  - - - - - - - - - - ')

            ! if ((.not. isSingularYN) .and. (setting%Time%Step > 102057) ) then 
            !     print *,' '
            !     print *,'geo 111 step ',setting%Time%Step ,'=============================='
            !     print *, elemR(113,er_Volume), elemR(113,er_Depth), elemR(113,er_Head)
            !     print *, elemR(113,er_FullVolume), elemR(113,er_FullDepth), elemR(113,er_Head) - elemR(113,er_Zcrown)
            !     print *, elemR(113,er_SlotVolume), elemR(113,er_SlotDepth), elemR(113,er_SlotWidth)
            !     print *, elemYN(113,eYN_isPSsurcharged), elemYN(113,eYN_isSurcharged)
            ! end if
        !% --- PREISSMAN SLOT    
        !% --- Handle Preissmann Slot for closed CC elements
        !%     with this time march type.
        if (npackP_Closed > 0) then
            call slot_CC (thisP_Closed, isSingularYN)
        end if

        ! if (.not. isSingularYN) call util_utest_CLprint('       222 geo  - - - - - - - - - - ')
        ! if (.not. isSingularYN) print *, elemR(113,er_Volume), elemR(113,er_Depth), elemR(113,er_Head)

        ! if ((.not. isSingularYN) .and. (setting%Time%Step > 102057) ) then 
        !     print *,' '
        !     print *,'geo 222'
        !     print *, elemR(113,er_Volume), elemR(113,er_Depth), elemR(113,er_Head)
        !     print *, elemR(113,er_FullVolume), elemR(113,er_FullDepth), elemR(113,er_Head) - elemR(113,er_Zcrown)
        !     print *, elemR(113,er_SlotVolume), elemR(113,er_SlotDepth), elemR(113,er_SlotWidth)
        !     print *, elemYN(113,eYN_isPSsurcharged), elemYN(113,eYN_isSurcharged)
        ! end if
        !% --- DEPTH
        !%     compute the depth on all elements of CC based on geometry.
        !%     If surcharged, this call returns the full depth of a closed conduit 
        !%     without adding Preissmann Slot depth.
        if (isAllYN) then
            call geo_depth_from_volume_by_type_allCC (elemPGetm, npack_elemPGetm, col_elemPGetm)

            ! if (.not. isSingularYN)  call util_utest_CLprint('       333 geo  - - - - - - - - - - ')
        else
            call geo_depth_from_volume_by_element_CC (thisP, npackP)

            ! if (.not. isSingularYN) call util_utest_CLprint('       444 geo  - - - - - - - - - - ')
        end if

        ! if ((.not. isSingularYN) .and. (setting%Time%Step > 102057)) then 
        !     print *,' '
        !     print *,'geo 444'
        !     print *, elemR(113,er_Volume), elemR(113,er_Depth), elemR(113,er_Head)
        !     print *, elemR(113,er_FullVolume), elemR(113,er_FullDepth), elemR(113,er_Head) - elemR(113,er_Zcrown)
        !     print *, elemR(113,er_SlotVolume), elemR(113,er_SlotDepth), elemR(113,er_SlotWidth)
        !     print *, elemYN(113,eYN_isPSsurcharged), elemYN(113,eYN_isSurcharged)
        ! end if

        !% --- ZERO DEPTH CC
        !%     reset all zero or near-zero depths in CC
        !%     Arguably this should not be needed as the individual depth computations
        !%     in geo_depth_from_volume_by_type_CC should use the zerovalues as minimums
        !%     but this needs to be confirmed.
        call adjust_limit_by_zerovalues &
            (er_Depth, setting%ZeroValue%Depth, thisP, .false.)

            ! if (.not. isSingularYN) call util_utest_CLprint('       555 geo  - - - - - - - - - - ')
            ! if ((.not. isSingularYN) .and. (setting%Time%Step > 102057) ) then 
            !     print *,' '
            !     print *,'geo 555'
            !     print *, elemR(113,er_Volume), elemR(113,er_Depth), elemR(113,er_Head)
            !     print *, elemR(113,er_FullVolume), elemR(113,er_FullDepth), elemR(113,er_Head) - elemR(113,er_Zcrown)
            !     print *, elemR(113,er_SlotVolume), elemR(113,er_SlotDepth), elemR(113,er_SlotWidth)
            !     print *, elemYN(113,eYN_isPSsurcharged), elemYN(113,eYN_isSurcharged)
            ! end if

        !% --- PIEZOMETRIC HEAD
        !%     compute the head on all elements of CC
        !%     This sets head consistent with depth computed in geo_depth_from_volume
        !%     Head is strictly limited to the max depth + zbottom so it does not
        !%     include surcharge effects     
        elemR(thisP,er_Head) = llgeo_head_from_depth_pure &
                                    (thisP, elemR(thisP,er_Depth))

            ! if (.not. isSingularYN) call util_utest_CLprint('       666 geo  - - - - - - - - - - ')

            ! if ((.not. isSingularYN) .and. (setting%Time%Step > 102058) ) then 
            !     print *,' '
            !     print *,'geo 666'
            !     print *, elemR(113,er_Volume), elemR(113,er_Depth), elemR(113,er_Head)
            !     print *, elemR(113,er_FullVolume), elemR(113,er_FullDepth), elemR(113,er_Head) - elemR(113,er_Zcrown)
            !     print *, elemR(113,er_SlotVolume), elemR(113,er_SlotDepth), elemR(113,er_SlotWidth)
            !     print *, elemYN(113,eYN_isPSsurcharged), elemYN(113,eYN_isSurcharged)
            ! end if
        !% --- OPEN CHANNEL OVERFLOW
        !%     Compute the overflow lost for CC open channels above
        !%     their maximum volume (no ponding allowed from open CC). 
        !%     Note that overflow or ponding for JM elements is handled 
        !%     in slot_JM.
        !%     Note, this is NOT standard in EPA-SWMM
        !%     HACK 20230508 This needs further checking
        if (npackP_Open > 0) then
            if (setting%Discretization%AllowChannelOverflowTF) then
                call geo_overflow_openchannels (thisP_Open)
            end if
        end if

        ! if (.not. isSingularYN) call util_utest_CLprint('       777 geo  - - - - - - - - - - ')

        !% --- PREISSMAN SLOT VOLUME LIMIT CLOSED CONDUIT CC
        !%     limit the volume in closed element (CC) to the full volume
        !%     Note the excess volume has already been stored in the Preissman Slot
        if (npackP_Closed > 0) then
            call geo_volumelimit_closed (thisP_Closed)
        end if

        ! if (.not. isSingularYN)  call util_utest_CLprint('       888 geo  - - - - - - - - - - ')

        !% --- CROSS-SECTIONAL AREA
        !%     compute area from volume for CC
        !%     For closed conduits this is based on the volume limited by full volume.
        !%     For open channels the volume limit depends on if AllowChanneOverflowTF is false.
        elemR(thisP,er_Area) = llgeo_area_from_volume_pure(thisP,elemR(thisP,er_Volume))
        elemR(thisP,er_Area) = max(elemR(thisP,er_Area),setting%ZeroValue%Area)

        ! if (.not. isSingularYN) call util_utest_CLprint('       999 geo  - - - - - - - - - - ')

        !% --- TOPWIDTH CC
        !%     compute topwidth from depth for all CC
        !%     Note: volume is limited to full depth UNLESS AllowChannelOverflowTF is false
        if (isAllYN) then
            call geo_topwidth_from_depth_by_type_allCC (elemPGetm, npack_elemPGetm, col_elemPGetm)
        else
            call geo_topwidth_from_depth_by_element_CC (thisP, npackP)
        end if

        ! if (.not. isSingularYN) call util_utest_CLprint('       aaa geo  - - - - - - - - - - ')

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

        ! if (.not. isSingularYN) call util_utest_CLprint('       bbb geo  - - - - - - - - - - ')

        !% --- ELLDEPTH MODIFIED HYDRAULIC DEPTH
        !%     the modified hydraulic depth "ell" is used for 
        !%     for Froude number computations on all CC elements
        call geo_elldepth_from_head_CC (thisP)

        ! if (.not. isSingularYN) call util_utest_CLprint('       ccc geo  - - - - - - - - - - ')

        !% ---- ADJUST SLOT 
        !%      make adjustments for slots on closed elements only
        !%     These add slot values to volume, depth, head
        if (npackP_Closed > 0) then
            call slot_CC_adjustments (thisP_Closed)
        end if

        ! if (.not. isSingularYN)  call util_utest_CLprint('       ddd geo  - - - - - - - - - - ')

    end subroutine geometry_toplevel_CC
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
                    elemR(mm,er_AreaBelowBreadthMax) = llgeo_rect_round_area_from_depth_singular &
                                                        (mm, depthAtMaxBreadth(mm), zeroR)
                end do

            case (rect_triang)
                elemR(thisP,er_FullTopwidth) = elemR(thisP,er_BreadthMax)
                do ii=1,size(thisP)
                    mm = thisP(ii)
                    elemR(mm,er_AreaBelowBreadthMax) = llgeo_rectangular_triangular_area_from_depth_singular &
                                                        (mm, depthAtMaxBreadth(mm), zeroR)
                end do

            case default
                print *, 'CODE ERROR Unexpected case default'
                call util_crashpoint(5298733)
        end select

        !% --- temporary store of depth
        tempDepth(thisP) = depth(thisP)

        !% --- temporary store or computing full values with general functions
        elemR(thisP,er_Depth)     = elemR(thisP,er_FullDepth)
        elemR(thisP,er_Area)      = elemR(thisP,er_FullArea)
        elemR(thisP,er_Perimeter) = elemR(thisP,er_FullPerimeter)
        elemR(thisP,er_Topwidth)  = elemR(thisP,er_FullTopwidth)

        !% --- standard functions using temporary store
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
                print *, 'CODE ERROR Unexpected case default'
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
    
        thisArea      = geo_area_from_depth_singular      (eIdx, inDepth, ZeroValueArea)
        thisPerimeter = geo_perimeter_from_depth_singular (eIdx, inDepth, ZeroValuePerimeter)
        outvalue      = thisArea * ((thisArea / thisPerimeter)**twothirdR)

    end function geo_sectionfactor_from_depth_singular
!%
!%==========================================================================    
!%==========================================================================   
!%
    real(8) function geo_Qcritical_from_depth_singular &
         (eIdx, inDepth, ZeroValue) result (outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% computes the critical flow for element eIdx with depth "inDepth"
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in)  :: eIdx
            real(8), intent(in)  :: inDepth, ZeroValue
            real(8), pointer     :: grav
            real(8)              :: thisArea
        !%------------------------------------------------------------------
        !% Aliases
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
            real(8), pointer    :: thistable(:)
            real(8)             :: QcritNormalized
            integer :: utr_Max
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
                    print *, 'CODE ERROR unexpected case default'
                    call util_crashpoint(6209873)
            end select
        !%------------------------------------------------------------------
    
        !% --- normalize the critical flowrate
        QcritNormalized = abs(elemR(eIdx,er_Flowrate) / uniformTableR(UT_idx,utr_QcritMax))

        !% --- lookup the normalized critical depth for this critical flow
        outvalue = xsect_table_lookup_singular (QcritNormalized, thistable)

        !% --- return depth to physical value
        outvalue = outvalue * uniformTableR(UT_idx,utr_Max)

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

        !% --- if flow is negative on a positive slope, or flow is positive on a negative slope,
        !%     then the section factor is negative, which implies an infinite normal depth
        if (sectionFactor .le. zeroR) then
            outvalue = setting%Limiter%NormalDepthInfinite
            return
        end if

        !% --- normalize the section factor
        normSF   = sectionFactor / uniformTableR(UT_idx,utr_SFmax)

        !% --- lookup the normalized normal depth
        outvalue = xsect_table_lookup_singular(normSF,thisTable)

        !% --- return normal depth to physical value
        outvalue = outvalue * uniformTableR(UT_idx,utr_DepthMax)
    
    end function geo_normaldepth_singular
!%
!%==========================================================================
!% PRIVATE (except _singular functions)
!%==========================================================================
!%
    !% ARCHIVE METHOD  May need to reinstate
    ! subroutine geo_ponding_inflow (thisColP_JM)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Provides inflow volume from ponding outside of a JM element when
    !     !% the interior volume is below the maximum. Does NOT handle
    !     !% inflow into surcharge from ponding (see slot_JM)
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: thisColP_JM  !% packed column for JM elements
    !         integer, pointer    :: Npack, thisP(:)
    !         real(8), pointer    :: vPond(:), volume(:), vFull(:), vInflow(:)
    !     !%------------------------------------------------------------------
    !     !% Preliminaries
    !         if (.not. setting%SWMMinput%AllowPonding) return
    !         Npack => npack_elemP(thisColP_JM)
    !         if (Npack < 1) return
    !     !%------------------------------------------------------------------
    !     !% Aliases:
    !         thisp       => elemP(1:Npack,thisColP_JM)
    !         vPond       => elemR(:,er_VolumePonded)
    !         volume      => elemR(:,er_Volume)
    !         vFull       => elemR(:,er_FullVolume)
    !         vInflow     => elemR(:,er_Temp01)
    !     !%------------------------------------------------------------------
    !     !% --- volume available in the JM    
    !     vInflow(thisP) = vFull(thisP) - volume(thisP)
    !     !% --- inflow is the lesser of available volume and ponded volume
    !     vInflow(thisP)  = min(vInflow(thisp), vPond(thisP))

    !     where (vInflow(thisP) > zeroR)
    !         volume(thisP) = volume(thisP) + vInflow(thisP)
    !         vPond(thisP)  = vPond(thisp)  - vInflow(thisP)
    !     endwhere
  
    ! end subroutine geo_ponding_inflow
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
        !% Declarations
            integer, intent(in) :: thisP(:), Npack
            integer :: mm
            integer, dimension(1) :: ap
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
                    print *, 'CODE ERROR Unexpected case default'
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

        !% --- limit the full volume by the full volume
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
        !% by a Froude number of one.
        !%------------------------------------------------------------------
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
            !% --- thisColP_JM is the column for the junction mains, which should be ep_JM
            integer, intent(in) ::thisColP_JM

            integer, pointer :: Npack, thisP(:), BranchExists(:), thisSolve(:),  tM
            integer, pointer :: fup(:), fdn(:)
            real(8), pointer :: area(:), depth(:), head(:), hydradius(:), AreaVelocity(:)
            real(8), pointer :: length(:), perimeter(:), topwidth(:), velocity(:), flowrate(:)
            real(8), pointer :: volume(:), zBtm(:), Kfac(:), dHdA(:), ellDepth(:)
            real(8), pointer :: zCrown(:), fullArea(:), fulldepth(:), fullperimeter(:)
            real(8), pointer :: sedimentDepth(:), fulltopwidth(:), breadthmax(:)
            real(8), pointer :: fullhydradius(:), Atable(:), Ttable(:), Rtable(:), Stable(:)
            real(8), pointer :: grav, fVel_u(:), fVel_d(:) 
            real(8), pointer :: fZcrown_u(:), fZcrown_d(:), fHead_u(:), fHead_d(:)  

            integer :: tB, ii, kk, tBA(1)
            logical :: isUpBranch, isWaterfall

            character(64)  :: subroutine_name = 'geo_assign_JB_from_head'
        !%---------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%geometry) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            if (setting%Profile%useYN) call util_profiler_start (pfc_geo_assign_JB_from_head)
        !%----------------------------------------------------------------------
        !% Aliases
            Npack         => npack_elemP(thisColP_JM)
            area          => elemR(:,er_Area)
            AreaVelocity  => elemR(:,er_AreaVelocity)
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
            Kfac          => elemSR(:,esr_JB_Kfactor)
            BranchExists  => elemSI(:,esi_JB_Exists)
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
              
                !% cycle through the possible junction branches
                do kk=1,max_branch_per_node

                    if (mod(kk,2)==0) then 
                        isUpBranch = .false.
                    else
                        isUpBranch = .true.
                    endif
                    
                    tB = tM + kk !% junction branch ID
                    tBA(1) = tB  !% array for pure array functions

                    if (BranchExists(tB) .ne. oneI) cycle

                    if ( head(tM) > (zBtm(tB) + sedimentDepth(tB)) ) then

                        select case (setting%Junction%HeadMethodJB)
                        case (use_JM)
                            head(tB) = head(tM)
                        case (linear_interp)
                            !% --- head(tB) is linear interp from face value.
                            if (isUpBranch) then 
                                !% --- upstream branch
                                if (fHead_d(fup(tB)) > (zBtm(tB) + sedimentDepth(tB) + setting%ZeroValue%Depth)) then
                                    head(tB) = head(tM)  &
                                    + onehalfR * (fHead_d(fup(tB)) - head(tM))
                                else 
                                    head(tB) = head(tM) &
                                    + onehalfR * ((zBtm(tB) + sedimentDepth(tB)) - head(tM))
                                end if
                            else
                                !% --- downstream branch
                                if (fHead_u(fdn(tB)) > (zBtm(tB) + sedimentDepth(tB) + setting%ZeroValue%Depth)) then
                                    head(tB) = head(tM)  &
                                        + onehalfR * (fHead_u(fdn(tB)) - head(tM))
                                else 
                                    head(tB) = head(tM) !%&
                                    !+ onehalfR * ((zBtm(tB) + sedimentDepth(tB)) - head(tM))
                                end if
                            end if
                        case default 
                            print *, 'CODE ERROR: unknown case for setting.Junction.HeadMethodJB of ',setting%Junction%HeadMethodJB
                            call util_crashpoint(7209874)
                        end select

                        iswaterfall = .false.

                    else
                        !% --- for main head below the branch bottom entrance we assign a
                        !%     Froude number of one on an inflow to the junction main. Note
                        !%     an outflow from a junction main for this case gets head
                        !%     of z_bottom of the branch (zero depth).
                        !%     Note this is a time-lagged velocity as the JB velocity
                        !%      is not updated until after face interpolation

                        !% ARCHIVE METHOD
                        ! head(tB) = zBtm(tB) + sedimentDepth(tB)                            &
                        !     + onehalfR * (oneR + branchsign(kk) * sign(oneR,velocity(tB))) &
                        !     *(velocity(tB)**twoR) / (grav) 

                        !% --- using velocity on face  
                        if     (      isUpBranch) then 

                            if (fVel_d(fup(tB)) > zeroR) then
                                !% --- upstream branch with waterfall inflow to junction main
                                iswaterfall = .true.

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
                            else
                                iswaterfall = .false.
                                !% --- no flow, set below zerovalue
                                head(tB) = zBtm(tB) + sedimentDepth(tB) + 0.99d0 * setting%ZeroValue%Depth
                            end if

                        elseif (.not. isUpbranch) then

                            if (fVel_u(fdn(tB)) < zeroR) then
                                !% --- downstream branch with inflow waterfall into junction main

                                iswaterfall = .true.
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
                            else 
                                iswaterfall = .false.
                                !% --- no flow, set below zerovalue
                                head(tB) = zBtm(tB) + sedimentDepth(tB) + 0.99d0 * setting%ZeroValue%Depth

                            end if
                        end if
                    end if

                    !% compute provisional depth
                    depth(tB) = head(tB) - (zBtm(tB) + sedimentDepth(tB))

                    if (depth(tB) .ge. fulldepth(tB)) then
                        !%--- surcharged or incipient surcharged
                        depth(tB)        = fulldepth(tB)
                        area(tB)         = fullArea(tB)
                        AreaVelocity(tB) = fullArea(tB)
                        perimeter(tB)    = fullperimeter(tB)
                        topwidth(tB)     = setting%ZeroValue%Topwidth
                        hydRadius(tB)    = fulldepth(tB) / fullperimeter(tB)
                        dHdA(tB)         = oneR / setting%ZeroValue%Topwidth
                        ellDepth(tBA)    = llgeo_elldepth_pure(tBA)
                        elemYN(tB,eYN_isSurcharged) = .true.

                    elseif (depth(tB) .le. setting%ZeroValue%Depth) then
                        !% --- negligible depth is treated with ZeroValues
                        depth(tB)        = setting%ZeroValue%Depth * 0.99d0
                        area(tB)         = setting%ZeroValue%Area
                        AreaVelocity(tB) = setting%ZeroValue%Area
                        topwidth(tB)     = setting%ZeroValue%Topwidth
                        perimeter(tB)    = setting%ZeroValue%Topwidth + setting%ZeroValue%Depth
                        hydRadius(tB)    = setting%ZeroValue%Depth
                        dHdA(tB)         = oneR / setting%ZeroValue%Topwidth
                        ellDepth(tB)     = setting%ZeroValue%Depth * 0.99d0 
                        elemYN(tB,eYN_isSurcharged) = .false.

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

                                ellDepth(tB)  = max(ellDepth(tB),setting%ZeroValue%Depth*0.99d0)

                            case (trapezoidal) !% analytical
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

                                topwidth(tB) = irregular_geometry_from_depth_singular ( &
                                                    tB,tt_width, depth(tB), elemR(tB,er_FullTopWidth),setting%ZeroValue%TopWidth)

                                !% --- note the irregular stores hyd radius rather than perimeter
                                hydRadius(tB) = irregular_geometry_from_depth_singular ( &
                                                    tB,tt_hydradius, depth(tB), elemR(tB,er_FullHydRadius), setting%ZeroValue%Depth)    

                                !% --- perimeter is derived geometry for irregular
                                perimeter(tB) = area(tB) / hydRadius(tB)
                                perimeter(tB)= max(perimeter(tB),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)

                                !% --- irregular must be continuously-increasing in width
                                ellDepth(tB)  = geo_hyddepth_from_area_and_topwidth_singular (area(tB), topwidth(tB), setting%ZeroValue%Depth*0.99d0) 

                            !% --- CLOSED CONDUITS
                            !%     closed conduits typically have look-up functions for area, topwidth and hydraulic
                            !%     radius, with standard geo_functions for perimeter, hydraulic depth, and ell
                            !%     However, where analytical functions are used, the perimeter is usually computed
                            !%     first and hydraulic radius is a geo_ function

                            !% --- lookups with Hydraulic Radius stored
                            case (arch, basket_handle, circular, eggshaped, horiz_ellipse, horseshoe, vert_ellipse)                                 
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
                                area(tB)      = llgeo_mod_basket_area_from_depth_singular        (tB,depth(tB),setting%ZeroValue%Area)
                                topwidth(tB)  = llgeo_mod_basket_topwidth_from_depth_singular    (tB,depth(tB),setting%ZeroValue%Topwidth)
                                perimeter(tB) = llgeo_mod_basket_perimeter_from_depth_singular   (tB,depth(tB),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                                
                                hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                    (tBA, area(tBA), perimeter(tBA))
                                hydRadius(tB) = max(hydRadius(tB),setting%ZeroValue%Depth)

                                ellDepth(tBA) = llgeo_elldepth_pure (tBA)
                                ellDepth(tB)  = max(ellDepth(tB),setting%ZeroValue%Depth*0.99d0)

                            case (rectangular_closed) !% analytical
                                area(tB)      = llgeo_rectangular_closed_area_from_depth_singular      (tB, depth(tB),setting%ZeroValue%Area)
                                topwidth(tB)  = llgeo_rectangular_closed_topwidth_from_depth_singular  (tB, depth(tB),setting%ZeroValue%Topwidth)
                                perimeter(tB) = llgeo_rectangular_closed_perimeter_from_depth_singular (tB, depth(tB),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                                
                                hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                    (tBA, area(tBA), perimeter(tBA))
                                hydRadius(tB) = max(hydRadius(tB),setting%ZeroValue%Depth)

                                ellDepth(tBA) = llgeo_elldepth_pure (tBA) 
                                ellDepth(tB)  = max(ellDepth(tB),setting%ZeroValue%Depth*0.99d0)

                            case (rect_round)  !% analytical                            
                                area(tB)      = llgeo_rect_round_area_from_depth_singular       (tB,depth(tB),setting%ZeroValue%Area)
                                topwidth(tB)  = llgeo_rect_round_topwidth_from_depth_singular   (tB,depth(tB),setting%ZeroValue%Topwidth)
                                perimeter(tB) = llgeo_rect_round_perimeter_from_depth_singular  (tB,depth(tB),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                                
                                hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                    (tBA, area(tBA), perimeter(tBA))
                                hydRadius(tB) = max(hydRadius(tB),setting%ZeroValue%Depth)

                                ellDepth(tBA) = llgeo_elldepth_pure (tBA) 
                                ellDepth(tB)  = max(ellDepth(tB),setting%ZeroValue%Depth*0.99d0)

                            case (rect_triang) !% analytical                                
                                area(tB)      = llgeo_rectangular_triangular_area_from_depth_singular      (tB,depth(tB),setting%ZeroValue%Area)
                                topwidth(tB)  = llgeo_rectangular_triangular_topwidth_from_depth_singular  (tB,depth(tB),setting%ZeroValue%Topwidth)
                                perimeter(tB) = llgeo_rectangular_triangular_perimeter_from_depth_singular (tB,depth(tB),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                                
                                hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                    (tBA, area(tBA), perimeter(tBA))
                                hydRadius(tB) = max(hydRadius(tB),setting%ZeroValue%Depth)

                                ellDepth(tBA) = llgeo_elldepth_pure (tBA) 
                                ellDepth(tB)  = max(ellDepth(tB),setting%ZeroValue%Depth*0.99d0)
                        
                            case default
                                print *, 'CODE ERROR geometry type unknown for # ', elemI(tB,ei_geometryType)
                                print *, 'which has key ',trim(reverseKey(elemI(tB,ei_geometryType)))
                                print *, 'in ',trim(subroutine_name)
                                call util_crashpoint(399848)

                        end select

                        elemYN(tB,eYN_isSurcharged) = .false.
                        AreaVelocity(tB) = area(tB)

                        !% --- standard for all geometries
                        dHdA(tB)     = oneR / topwidth(tB)

                    end if
                    !% --- universal computation of volume
                    volume(tB) = area(tB) * length(tB)

                    if (iswaterfall) then 
                        !% --- limit velocity and flowrate by Fr=1
                        !%     Note this may be inconsistent with flowrate
                        !%     set by mass conservation of junction
                        if (isUpBranch) then
                            velocity(tB) = sqrt(grav * depth(tB))
                        else
                            velocity(tB) = -sqrt(grav * depth(tB))
                        end if
                    else
                        !% --- universal computation of velocity
                        if (area(tB) > setting%ZeroValue%Area) then 
                            velocity(tB) = flowrate(tB) / AreaVelocity(tB)

                            !% ARCHIVE METHOD
                            !% -- set slightly larger depths to velocity consistent with FR=1
                            ! if ((area(tB) < tenR * setting%ZeroValue%Area) .or. &
                            !     (elldepth(tB) < tenR * setting%ZeroValue%Depth) ) then 
                            !     velocity(tB) =sign(oneR,flowrate(tB)) * sqrt(grav * elldepth(tB))
                            ! end if

                        else
                            velocity(tB) = zeroR
                        end if
                    end if


                    if (abs(velocity(tB)) > setting%Limiter%Velocity%Maximum) then
                        velocity(tB) = sign(0.99d0 * setting%Limiter%Velocity%Maximum, velocity(tB)) 
                    end if
                            
                end do
            end do
        end if

        !% --- Note, the above can only be made a concurrent loop if we replace the tM
        !%     with thisP(ii) and tB with thisP(ii)+kk, which makes the code
        !%     difficult to read.

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Profile%useYN) call util_profiler_stop (pfc_geo_assign_JB_from_head)

            if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine geo_assign_JB_from_head
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
        !% Closing
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
                        print *, 'CODE ERROR Unexpected case default'
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

        !% --- FILLED CIRCULAR
        Npack => npack_elemPGx(epg_CC_filled_circular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_filled_circular)
            thisP   => elemPGx(1:Npack,thisCol)
            call filled_circular_hydradius_and_perimeter_from_depth (thisP)
            perimeter(thisP) = max(perimeter(thisP),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
            hydradius(thisP) = max(hydradius(thisP),setting%ZeroValue%Depth)

        end if

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

        !% TOPWIDTH is UNDEFINED FOR TABULAR, FUNCTIONAL, AND IMPLIED STORAGE

        !%-------------------------------------------------------------------
        !% Closing
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
                    print *, 'CODE ERROR Unexpected case default'
                    call util_crashpoint(7209874)
            end select

        end do

    end subroutine geo_perimeter_and_hydradius_from_depth_by_element_CC
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
       
        where ((head(thisP) .le. ZbreadthMax(thisP)) &
             .and. (topwidth(thisP) > setting%ZeroValue%Topwidth))
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
        !% Closing
            if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine geo_elldepth_from_head_CC
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

        case (custom)
            print *, 'CODE ERROR area for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(332341)

        case (force_main)
            print *, 'CODE ERROR area for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'in ',trim(subroutine_name)   
            print *, 'This should never be reached as a force_main is not a valid geometryType'
            call util_crashpoint(332342)

        case default
            print *, 'CODE ERROR area for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
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
            print *, 'CODE ERROR topwidth for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
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
        
        case (custom)
            print *, 'CODE ERROR topwidth for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)

        case (force_main)
            print *, 'CODE ERROR topwidth for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'in ',trim(subroutine_name)   
            print *, 'This should never be reached as a force_main is not a valid geometryType'
            call util_crashpoint(4498734)

        case default
            print *, 'CODE ERROR topwidth for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        end select

        
           
    end function geo_topwidth_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function geo_perimeter_from_depth_singular &
        (idx, indepth, ZeroValuePerimeter) result (outvalue)  !% PRIVATE
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
            print *, 'CODE ERROR perimeter for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
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
        
        case (custom)
            print *, 'CODE ERROR perimeter for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(338234)
        case (force_main)
            print *, 'CODE ERROR perimeter for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'in ',trim(subroutine_name)   
            print *, 'This should never be reached as a force_main is not a valid geometryType' 
            call util_crashpoint(338234)
        case default
            print *, 'CODE ERROR perimeter for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(332344)
        end select

    end function geo_perimeter_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function geo_hyddepth_from_area_and_topwidth_singular &
        (area, topwidth, ZeroValueHydDepth)  result (outvalue)
        !%------------------------------------------------------------------
        !% Descriptions:
        !% computes the hydraulic depth for area and topwidth of a single
        !% element
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(in)  :: area, topwidth, ZeroValueHydDepth
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
    subroutine geo_ZeroDepth_from_volume (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% ensures that if volume <= zeroArea * length the depth will be
        !% zeroDepth
        !%------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: depth(:), volume(:), length(:)
        !%------------- -----------------------------------------------------
        !% Aliases
            depth       => elemR(:,er_Depth)
            volume      => elemR(:,er_Volume)
            length      => elemR(:,er_Length)
        !%------------- ----------------------------------------------------

        where (volume(thisP)/length(thisP) .le. setting%ZeroValue%Area)
            depth(thisP) = setting%ZeroValue%Depth * 0.99d0
        end where

    end subroutine geo_ZeroDepth_from_volume
!%
!%==========================================================================
!% END OF MODULE
!%=========================================================================
end module geometry