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
    use utility, only: util_CLprint, util_syncwrite


    implicit none

!%-----------------------------------------------------------------------------
!% Description:
!% Geometry computations
!%

    private

    public :: geometry_toplevel
    public :: geo_common_initialize
    public :: geo_sectionfactor_from_depth_singular
    public :: geo_Qcritical_from_depth_singular
    public :: geo_criticaldepth_singular
    public :: geo_normaldepth_singular
    !public :: geo_topwidth_from_depth_by_type
    public :: geo_hyddepth_from_area_and_topwidth_singular
    public :: geo_topwidth_from_depth_singular
    public :: geo_area_from_depth_singular
    public :: geo_perimeter_from_depth_singular
    !public :: geo_ell_pure
    !public :: geometry_table_initialize

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine geometry_toplevel (whichTM)
        !%------------------------------------------------------------------
        !% Description:
        !% Input whichTM is one of ETM, AC, or ALLtm
        !% This should never be called for diagnostic arrays
        !% Note that the elemPGx arrays contain only time-marched elements so they
        !% will only handle CC and JM elements as the JB elements are not time-marched.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: whichTM
            integer, pointer :: elemPGx(:,:), npack_elemPGx(:), col_elemPGx(:)
            !integer, pointer :: thisColP_surcharged, thisColP_NonSurcharged, 
            integer, pointer :: thisColP_all_TM, thisColP_Open_CC
            integer, pointer :: thisColP_JM, thisColP_JB, thisColP_CC
            integer, pointer :: thisColP_Closed_CC, thisColP_Closed_JB
            integer, pointer :: Npack, thisP(:)
            !integer, pointer ::  thisColP_Closed_JM
            logical :: isreset
            real(8), pointer :: depth(:), volume(:), area(:), topwidth(:)
            integer, allocatable :: tempP(:) !% debugging
            character(64) :: subroutine_name = 'geometry_toplevel'
        !%------------------------------------------------------------------
        !% Preliminaries
            !!if (crashYN) return
            if (setting%Debug%File%geometry) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases
            area   => elemR(:,er_Area)
            depth  => elemR(:,er_Depth)
            topwidth => elemR(:,er_Topwidth)
            volume => elemR(:,er_Volume)
        !% set the packed geometry element array (elemPG) to use and columns of the
        !% packed elemP to use
            select case (whichTM)
                ! case (ALLtm)
                !     elemPGx                => elemPGalltm(:,:) 
                !     npack_elemPGx          => npack_elemPGalltm(:)
                !     col_elemPGx            => col_elemPGalltm(:)
                !     thisColP_CC            => col_elemP(ep_CC_ALLtm)
                !     thisColP_JM            => col_elemP(ep_JM_ALLtm)
                !     thisColP_JB            => col_elemP(ep_JB_ALLtm)
                !     !thisColP_surcharged    => col_elemP(ep_ALLtmSurcharged)
                !     !thisColP_NonSurcharged => col_elemP(ep_ALLtm_NonSurcharged)
                !     thisColP_all           => col_elemP(ep_ALLtm)
                !     thisColP_Open_CC       => col_elemP(ep_CC_Open_Elements)
                !     thisColP_Closed_CC     => col_elemP(ep_CC_Closed_Elements)
                !     thisColP_Closed_JB     => col_elemP(ep_JB_Closed_Elements)
                !     !thisColP_Closed_JM     => col_elemP(ep_JM_Closed_Elements)
                case (ETM)
                    elemPGx                => elemPGetm(:,:)
                    npack_elemPGx          => npack_elemPGetm(:)
                    col_elemPGx            => col_elemPGetm(:)
                    thisColP_CC            => col_elemP(ep_CC_ETM)
                    thisColP_JM            => col_elemP(ep_JM_ETM)
                    !thisColP_JB            => col_elemP(ep_JB_ETM)
                    !thisColP_surcharged    => col_elemP(ep_PSsurcharged)
                    !thisColP_NonSurcharged => col_elemP(ep_ETM_PSnonSurcharged)
                    thisColP_all_TM        => col_elemP(ep_ETM)
                    thisColP_Open_CC       => col_elemP(ep_CC_Open_Elements)
                    thisColP_Closed_CC     => col_elemP(ep_CC_Closed_Elements)
                    thisColP_Closed_JB     => col_elemP(ep_JB_Closed_Elements)
                    !thisColP_Closed_JM     => col_elemP(ep_JM_Closed_Elements)
                ! case (AC)
                !     elemPGx                => elemPGac(:,:)
                !     npack_elemPGx          => npack_elemPGac(:)
                !     col_elemPGx            => col_elemPGac(:)
                !     thisColP_CC            => col_elemP(ep_CC_AC)
                !     thisColP_JM            => col_elemP(ep_JM_AC)
                !     thisColP_JB            => col_elemP(ep_JB_AC)
                !     !thisColP_surcharged    => col_elemP(ep_ACsurcharged)
                !     !thisColP_NonSurcharged => col_elemP(ep_AC_ACnonSurcharged)
                !     thisColP_all           => col_elemP(ep_AC)
                !     thisColP_Open_CC       => col_elemP(ep_CC_Open_Elements)
                !     thisColP_Closed_CC     => col_elemP(ep_CC_Closed_Elements)
                !     thisColP_Closed_JB     => col_elemP(ep_JB_Closed_Elements)
                !     !thisColP_Closed_JM     => col_elemP(ep_JM_Closed_Elements)
                case default
                    print *, 'CODE ERROR: time march type unknown for # ', whichTM
                    print *, 'which has key ',trim(reverseKey(whichTM))
                    call util_crashpoint(7389)
                    !return
                    !stop 7389
            end select
            call util_crashstop(49872)
        !%--------------------------------------------------------------------
            !   call util_CLprint ('in geometry at top==========================================')  

        !% STATUS: at this point we know volume and velocity on all elements
        !% from RK2 solution

        !% --- PONDING
        !%     adjust time-march volume for ponding inflow
        !%     This affects JM that have previous ponding but now have volumes
        !%     below the full volume
        if (setting%SWMMinput%AllowPonding) then
            call geo_ponding_inflow (thisColP_JM)   
        end if

            ! call util_CLprint ('in geometry after ponding inflow') 
           
        !% --- PREISSMAN SLOT    
        !%     also adds to ponded volume or computes overflow volume for JM (only)
        !%     where slot depth + invert height exceeds maximum surcharge height
        !%     The Element Volume in JM is adjusted for any overflow or ponding
        !%     (but not for the slot)
        call slot_toplevel (whichTM, thisColP_Closed_CC, thisColP_JM)

            ! call util_CLprint ('in geometry after slot toplevel') 

        !% STATUS: The Preissmann Slot values have been assigned for all CC and JM
        !% Overflow and Ponding have been assigned for JM (only). 
        !% JM element volumes have been adjusted for overflow or ponding, but
        !% not for slot volume. So both JM and CC have volume > fullvolume
        !% in surcharged elements.

        !% --- assign all geometry for surcharged elements CC, JM
        !%     Note: not used in Preissmann Slot (only AC)
        !%     DO NOT DELETE. HOLD THIS FOR LATER USE 20220909brh
        ! if ((whichTM .eq. ALLtm) .or. (whichTM .eq. AC)) then
        !     call geo_ACsurcharged (thisColP_surcharged)
        ! end if

        !% --- ZERO VOLUMES CC JM
        !%     reset all zero or near-zero volumes in all CC, JM
        !call adjust_limit_by_zerovalues (er_Volume, setting%ZeroValue%Volume, thisColP_NonSurcharged, .true.)
        call adjust_limit_by_zerovalues &
            (er_Volume, setting%ZeroValue%Volume, thisColP_all_TM, .true.)

            ! call util_CLprint ('in geometry after limit_by_zerovalues (volume)') 

        !% --- DEPTH
        !%     compute the depth on all elements of CC JM based on geometry.
        !%     If surcharged, this call returns the full depth of a closed conduit 
        !%     without adding Preissmann Slot depth.
        call geo_depth_from_volume_by_type (elemPGx, npack_elemPGx, col_elemPGx)

         !    call util_CLprint ('in geometry after depth_from_volume') 

        !% --- ZERO DEPTH CC JM
        !%     reset all zero or near-zero depths in aa CC and JM
        !call adjust_limit_by_zerovalues (er_Depth, setting%ZeroValue%Depth, thisColP_NonSurcharged, .false.)
        call adjust_limit_by_zerovalues &
            (er_Depth, setting%ZeroValue%Depth, thisColP_all_TM, .false.)

            !call util_CLprint ('in geometry after limit_by_zerovalues (depth)') 

        !% --- PIEZOMETRIC HEAD
        !%     compute the head on all elements of CC and JM
        !%     This sets head consistent with depth computed in geo_depth_from_volume
        !%     Head is strictly limited to the max depth + zbottom so it does not
        !%     include surcharge effects     
        Npack     => npack_elemP(thisColP_all_TM)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisColP_all_TM)
            elemR(thisP,er_Head) = llgeo_head_from_depth_pure (thisP, depth(thisP))
            !elemR(thisP,er_PressureHead) = llgeo_head_from_depth_pure (thisP)
            !call geo_head_from_depth (thisColP_all_TM)
        end if
 
            ! call util_CLprint ('in geometry after head_from_depth')
        
        !% --- OPEN CHANNEL OVERFLOW
        !%     Compute the overflow lost for CC open channels above
        !%     their maximum volume (no ponding allowed from open CC). 
        !%     Note that overflow or ponding for JM elements is handled 
        !%     in slot_JM_ETM.
        !%     Note, this is NOT standard in EPA-SWMM
        if (setting%Discretization%AllowChannelOverflowTF) then
            call geo_overflow_openchannels (thisColP_Open_CC)
        end if

            ! call util_CLprint ('in geometry after overflow_openchannels')

        !% --- PREISSMAN SLOT VOLUME LIMIT CLOSED CONDUIT CC JM
        !%     limit the volume in closed element (CC, JM) to the full volume
        !%     Note the excess volume has already been stored in the Preissman Slot
        call geo_volumelimit_closed (thisColP_Closed_CC)
        call geo_volumelimit_closed (thisColP_JM)

            ! call util_CLprint ('in geometry after volumelimit_closed')

        !% REMOVED 20220909 brh
        !% --- limit volume for incipient surcharge. This is done after depth is computed
        !%     so that the "depth" algorithm can include depths greater than fulldepth
        !%     as a way to handle head for incipient surcharge.
        !call geo_limit_incipient_surcharge (er_Volume, er_FullVolume, thisColP_NonSurcharged,.true.) !% 20220124brh
            ! call util_CLprint ('in geometry before geo_limit_incipient_surcharge (Depth)')  
        !% --- limit depth for surcharged on CC. This is done after head is computed
        !%     so that the depth algorithm can include depths greater than fulldepth where the 
        !%     geometry algorithm does not enforrce full depth
        !call geo_limit_incipient_surcharge (er_Depth, er_FullDepth, thisColP_NonSurcharged,.false.) 
        !@call geo_limit_incipient_surcharge (er_Depth, er_FullDepth, thisColP_all,.false.) 
        !% END REMOVE 20220909

        !% STATUS: At this point, the depths, heads and volumes of all CC, JM elements are
        !% at or below their full value.  For CC closed conduits and all JM the
        !% surcharged volume is stored in the er_SlotVolume and the surcharged extra
        !% depth is stored in the er_SlotDepth 

        !% --- PREISSMAN SLOT HEAD ADD IN JM 
        !%     adjust JM head to include Preissmann Slot Depth and ponding
        !%     This is needed before JB are computed
        call slot_JM_head_PSadd (thisColP_JM)

            ! call util_CLprint ('in geometry after JM_head_PSadd') 
           
        !% --- JB VALUES
        !%    assign the non-volume geometry on junction branches JB based on JM head
        !%    Values limited by full volume. Volume assigned is area * length
        call geo_assign_JB (whichTM, thisColP_JM)

            ! call util_CLprint ('in geometry after assign_JB') 

        !% --- JB CLOSED CONDUIT VOLUME LIMIT
        !%     further limiting the JB volume by full is probably not needed,
        !%     but might be useful if there's a numerical precision issues
        !%     with JB volume assigned by area * length.
        call geo_volumelimit_closed (thisColP_Closed_JB)

            ! call util_CLprint ('in geometry after volumelimit_closed') 

        !% --- PREISSMANN SLOT HEAD REMOVE IN JM
        !%     we need to remove the PS and ponding from the JM cells so that we can easily
        !%     compute other geometry without full JM causing problems
        call slot_JM_head_PSremove (thisColP_JM)

            ! call util_CLprint ('in geometry after JM_head_PSremove')  

        !% STATUS: at this point we have all geometry on CC, JM, JB that is
        !% limited by the full volume values. The CC and JM have slot values stored
        !% but no slot values have been computed for JB
        
        !% --- CROSS-SECTIONAL AREA
        !%    compute area from volume for CC, JM
        !%     For closed conduits this is based on the volume limited by full volume.
        !%     For open channels the volume limit depends on if AllowChanneOverflowTF is false.
        !%     Note that JB areas are already assigned in geo_assign_JB()
        Npack     => npack_elemP(thisColP_all_TM)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisColP_all_TM)
            elemR(thisP,er_Area) = llgeo_area_from_volume_pure (thisP, volume(thisP))
            !call geo_area_from_volume (thisColP_all_TM)
        end if

            ! call util_CLprint ('in geometry after area_from_volume') 

        !% --- ZERO AREA CC JM
        !%     reset all zero or near-zero areas in CC and JM
        call adjust_limit_by_zerovalues &
             (er_Area, setting%ZeroValue%Area, thisColP_all_TM, .false.)

            ! call util_CLprint ('in geometry after adjust_limit_by_zeroValues area')   

        !% --- TOPWIDTH CC
        !%     compute topwidth from depth for all CC
        !%     Note: Topwidth for JM is undefined in this subroutine
        !%     Note: volume is limited to full depth UNLESS AllowChannelOverflowTF is false
        call geo_topwidth_from_depth_by_type (elemPGx, npack_elemPGx, col_elemPGx)

            ! call util_CLprint ('in geometry after topwidth_from_depth') 

        !% --- ZERO TOPWIDTH CC
        !%     reset all zero or near-zero topwidth in CC 
        !%     but do not change the eYN(:,eYN_isZeroDepth) mask
        call adjust_limit_by_zerovalues &
             (er_Topwidth, setting%ZeroValue%Topwidth, thisColP_CC, .false.)

            ! call util_CLprint ('in geometry after adjust_limit_by_zerovalues topwidth') 

        !% --- PERIMETER AND HYDRAULIC RADIUS CC
        !%     compute hydraulic radius and perimeter
        !%     note these two are done together because for analytical cross-sections
        !%     we have equations for perimeter, whereas lookup cross-sections
        !%     have tables for hydraulic radius.
        call geo_perimeter_and_hydradius_from_depth_by_type (elemPGx, npack_elemPGx, col_elemPGx)  

        ! % --- compute perimeter from maximum depth for all CC
        ! %     Note: perimeter for JM is undefined in this subroutine
        !OBSOLETE  call geo_perimeter_from_depth (elemPGx, npack_elemPGx, col_elemPGx)

            ! call util_CLprint ('in geometry after perimeter from depth') 

        !% --- compute hyddepth
        !call geo_hyddepth_from_depth_or_topwidth (elemPGx, npack_elemPGx, col_elemPGx)
        !% 20220930 replace with unified call
        Npack     => npack_elemP(thisColP_CC)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisColP_CC)
            elemR(thisP,er_HydDepth) = llgeo_hyddepth_from_area_and_topwidth_pure &
                                        (thisP, area(thisP), topwidth(thisP))
        end if
        ! OBSOLETE call geo_hyddepth_from_area_and_topwidth (thisColP_CC)

            ! call util_CLprint ('in geometry after hyddepth_from_depth')

        !% --- compute hydradius
        !%     Note: cannot be used for JM unless perimeter is defined prior.
        !OBSOLETE call geo_hydradius_from_area_perimeter (thisColP_CC)

            ! call util_CLprint ('in geometry after hydradius_from_area_perimeter') 

        !% --- the modified hydraulic depth "ell" is used for 
        !%     for Froude number computations on all CC elements
        !%     Note: ell for JM is undefined in this subroutine
        call geo_ell_from_head (thisColP_CC)

        !% --- compute pressure head from the modified hydraulic depth
        Npack     => npack_elemP(thisColP_CC)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisColP_CC)
            !elemR(thisP,er_Pressure_Head) = llgeo_pressure_head_from_hyddepth_pure (thisP)
            !call geo_pressure_head_from_hyddepth (thisColP_CC)
        end if

            ! call util_CLprint ('in geometry after ell_from_head') 

        !% --- make adjustments for slots on closed elements only
        !%     These add slot values to volume, depth, head
        call slot_CC_adjustments (thisColP_Closed_CC)

            ! call util_CLprint ('in geometry after slot_CC_adustments') 

        call slot_JM_adjustments (thisColP_JM)

            ! call util_CLprint ('in geometry after slot_JM_adjustments') 

        !% --- compute the SlotDepth, SlotArea, SlotVolume and update
        !%     elem volume and depth for JB. Note elem head on JB either with or
        !%     without surcharge is assigned in geo_assign_JB
        call slot_JB_computation (thisColP_JM)
        
            ! call util_CLprint ('in geometry after slot_JB_computation') 

        !% Set JM values that are not otherwise defined
        !% HydDepth, ell. Note that topwidth, hydradius, perimeter are undefined.
        call geo_JM_values ()

            ! call util_CLprint ('in geometry after JM_values') 

        !% HOLD UNTIL AC RE-VISITED
        ! !% compute the dHdA that are only for AC nonsurcharged
        ! if (whichTM .ne. ETM) then
        !     call geo_dHdA (ep_AC_ACnonSurcharged)
        ! end if

            ! call util_CLprint ('in geometry at end') 

        !% --- check for crashpoint and stop here
        call util_crashstop(322983)

        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geometry_toplevel
!%
!%==========================================================================
!%==========================================================================    
!%
    subroutine geo_common_initialize &
            (thisP, thisGeoType, ATableType, TTableType)
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
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisP(:), thisGeoType
            real(8), intent(in) :: ATableType(:), TTableType(:)
            integer             :: ii, mm
            real(8), pointer    :: depth(:), tempDepth(:), fullArea(:)
            real(8), pointer    :: breadthMax(:), fullDepth(:), depthAtMaxBreadth(:)
            real(8), pointer    :: fullTopwidth(:), fullPerimeter(:)
            real(8)             :: topwidthDepth
        !%------------------------------------------------------------------
        !% Aliases
            depth             => elemR(:,er_Depth)
            tempDepth         => elemR(:,er_Temp02)
            fullArea          => elemR(:,er_FullArea)
            fullDepth         => elemR(:,er_FullDepth)
            fullPerimeter     => elemR(:,er_FullPerimeter)
            fullTopwidth      => elemR(:,er_FullTopWidth)
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
                    (mm, topwidthDepth, breadthMax(mm), setting%ZeroValue%Topwidth, TTableType)   

                elemR(mm,er_AreaBelowBreadthMax) =  llgeo_tabular_from_depth_singular           &
                    (mm, depthAtMaxBreadth(mm), fullArea(mm), setting%ZeroValue%Area, ATableType)
            end do

        case (filled_circular)
            do ii=1,size(thisP)
                mm = thisP(ii)

                !% --- lookup functions using specific depths
                topwidthDepth = fullDepth(mm) &
                                * setting%Discretization%FullConduitTopwidthDepthFraction 
            
                elemR(mm,er_FullTopwidth) = llgeo_filled_circular_topwidth_from_depth_singular     &
                                (mm, topwidthDepth)                                             

                elemR(mm,er_AreaBelowBreadthMax) = llgeo_filled_circular_area_from_depth_singular  &   
                                (mm, depthAtMaxBreadth(mm))
            end do

        case (mod_basket)
            do ii=1,size(thisP)
                mm = thisP(ii)

                !% --- lookup functions using specific depths
                topwidthDepth = fullDepth(mm) &
                                * setting%Discretization%FullConduitTopwidthDepthFraction 

                elemR(mm,er_FullTopwidth) = llgeo_mod_basket_topwidth_from_depth_singular     &
                                            (mm, topwidthDepth)                                             

                elemR(mm,er_AreaBelowBreadthMax) = llgeo_mod_basket_area_from_depth_singular  &   
                                            (mm, depthAtMaxBreadth(mm))
            end do

        case (rectangular_closed)
            elemR(thisP,er_FullTopWidth)        = elemR(thisP,er_BreadthMax)
            elemR(thisP,er_AreaBelowBreadthMax) = elemR(thisP,er_FullArea) 

        case (rect_round)
            elemR(thisP,er_FullTopwidth)        = elemR(thisP,er_BreadthMax)
            elemR(thisP,er_AreaBelowBreadthMax) = elemR(thisP,er_FullArea)

        case (rect_triang)
            elemR(thisP,er_FullTopwidth)        = elemR(thisP,er_BreadthMax)
            elemR(thisP,er_AreaBelowBreadthMax) = elemR(thisP,er_FullArea)


        case default
            print *, 'CODE ERROR: Unexpected case default'
            call util_crashpoint(5298733)
        end select

        !% temporary store of depth
        tempDepth(thisP) = depth(thisP)

        !% --- temporary store or computing full values with general functions
        elemR(thisP,er_Depth)     = elemR(thisP,er_FullDepth)
        elemR(thisP,er_Area)      = elemR(thisP,er_FullArea)
        elemR(thisP,er_Perimeter) = elemR(thisP,er_FullPerimeter)
        elemR(thisP,er_Topwidth)  = elemR(thisP,er_FullTopwidth)

        !% --- standard functions using temporary store
        elemR(thisP,er_FullEll)       = llgeo_FullEll_pure(thisP) 
        elemR(thisP,er_FullHydDepth)  = llgeo_hyddepth_from_area_and_topwidth_pure &
                                            (thisP, fullarea(thisP), fulltopwidth(thisP))
        elemR(thisP,er_FullHydRadius) = llgeo_hydradius_from_area_and_perimeter_pure &
                                            (thisP, fullarea(thisP), fullperimeter(thisP))
        
        !% --- reset the initial depth
        depth(thisP) = tempDepth(thisP)
        
        !% --- get IC data for Area
        select case (thisGeoType)
        case (arch, basket_handle, catenary, circular, eggshaped, gothic, &
                horiz_ellipse, horseshoe, semi_circular, semi_elliptical, &
                vert_ellipse)
            do ii=1,size(thisP)
                mm = thisP(ii)
                elemR(mm,er_Area) = llgeo_tabular_from_depth_singular &
                    (mm, depth(mm), fullArea(mm), setting%ZeroValue%Area, ATableType)
            end do

        case (filled_circular)
            do ii=1,size(thisP)
                mm = thisP(ii)
                elemR(mm,er_Area) = llgeo_filled_circular_area_from_depth_singular &
                                    (mm, depth(mm))
            end do

        case (mod_basket)
            do ii = 1,size(thisP)
                mm = thisP(ii)
                elemR(mm,er_Area) = llgeo_mod_basket_area_from_depth_singular &
                                    (mm, depth(mm))
            end do

        case (rectangular_closed)
            elemR(thisP,er_Area) = elemR(thisP,er_Depth) * elemR(thisP,er_Breadthmax)

        case (rect_round)
            do ii=1,size(thisP)
                mm = thisP(ii)
                elemR(mm,er_Area) = llgeo_rect_round_area_from_depth_singular &
                                    (mm, depth(mm))
            end do
            print*, elemR(thisP,er_Area), 'elemR(thisP,er_Area)'

        case (rect_triang)
            do ii=1,size(thisP)
                mm = thisP(ii)
                elemR(mm,er_Area) = llgeo_rectangular_triangular_area_from_depth_singular &
                                        (mm, depth(mm))
            end do
        case default
            print *, 'CODE ERROR: Unexpected case default'
            call util_crashpoint(5298722)
        end select

        elemR(thisP,er_AoverAfull) = elemR(thisP,er_Area)  / elemR(thisP,er_FullArea)
        elemR(thisP,er_YoverYfull) = elemR(thisP,er_Depth) / elemR(thisP,er_FullDepth)

        call slot_initialize (thisP)

        !% store IC data
        elemR(thisP,er_Area_N0)       = elemR(thisP,er_Area)
        elemR(thisP,er_Area_N1)       = elemR(thisP,er_Area)
        elemR(thisP,er_Volume)        = elemR(thisP,er_Area) * elemR(thisP,er_Length)
        elemR(thisP,er_Volume_N0)     = elemR(thisP,er_Volume)
        elemR(thisP,er_Volume_N1)     = elemR(thisP,er_Volume)

        !%--- clear the temp storage
        tempDepth(thisP) = nullvalueR

        !% note that er_Perimeter, er_Topwidth, er_HydRadius,, er_Ell are NOT initialized
    end subroutine geo_common_initialize    
!%
!%==========================================================================
!%==========================================================================    
!%
    real(8) function geo_sectionfactor_from_depth_singular &
         (eIdx,inDepth) result (outvalue)  
        !%------------------------------------------------------------------
        !% Description
        !% computes the section factor for element with index eIdx for
        !% the depth "inDepth"
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in)  :: eIdx
            real(8), intent(in)  :: inDepth
            real(8) :: thisPerimeter, thisArea
            character(64) :: subroutine_name = "geo_sectionfactor_from_depth_singular"
        !%------------------------------------------------------------------  
        ! print *, 'in ',trim(subroutine_name)    
        ! print *, 'input depth ',inDepth
        thisArea      = geo_area_from_depth_singular      (eIdx,inDepth)
        ! print *, '----- area     ',thisArea
        thisPerimeter = geo_perimeter_from_depth_singular (eIdx,inDepth)
        ! print *, '----- perimeter',thisPerimeter
        outvalue      = thisArea * ((thisArea / thisPerimeter)**twothirdR)
        ! print *, '----- sf       ',outvalue

    end function geo_sectionfactor_from_depth_singular
!%
!%==========================================================================    
!%==========================================================================   
!%
    real(8) function geo_Qcritical_from_depth_singular &
         (eIdx,inDepth) result (outvalue)
        !%------------------------------------------------------------------
        !% computes the critical flow for element eIdx with depth "inDepth"
        !%------------------------------------------------------------------
         !% Declarations
         integer, intent(in)  :: eIdx
         real(8), intent(in)  :: inDepth
         real(8), pointer     :: grav
         real(8)              :: thisArea
        !%------------------------------------------------------------------
            grav => setting%Constant%gravity
        !%------------------------------------------------------------------     
        thisArea      = geo_area_from_depth_singular (eIdx, inDepth)
        outvalue      = thisArea * sqrt(inDepth * grav)

    end function geo_Qcritical_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function geo_criticaldepth_singular (UT_idx) result (outvalue)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the critical depth for the uniformtable(UT_idx)
        !% using the flowrate in the associated element eIdx
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: UT_idx
            integer, pointer    :: eIdx
            real(8), pointer    :: gravity, thistable(:)
            real(8)             :: normFlowrate
            integer :: ii
        !%------------------------------------------------------------------
        !% Aliases
            eIdx      => uniformTableI(UT_idx,uti_elem_idx)
            thisTable => uniformTableDataR(UT_idx,:,utd_Qcrit_depth_nonuniform)
        !%------------------------------------------------------------------
        ! print *, 'UT_idx',UT_idx
        ! print *, 'eIdx  ',eIdx
        ! print *, 'flowrate     ', elemR(eIdx,er_Flowrate)
        ! print *, 'utr_Qcritmax ',utr_QcritMax
        ! print *, 'table        ', uniformTableR(UT_idx,utr_QcritMax)
        ! do ii=1,N_Elem(this_image())
        !    print *, ii, elemR(ii,er_Flowrate)
        ! end do
        
        !stop 2098374

        !% --- normalize the critical flowrate
        normFlowrate = abs(elemR(eIdx,er_Flowrate) / uniformTableR(UT_idx,utr_QcritMax))
        
        ! print *, 'normflowrate ',normFlowrate

        !% --- lookup the normalized critical depth for this critical flow
        outvalue = xsect_table_lookup_singular (normFlowrate, thistable)

        ! print *, 'outvalue 1',outvalue

        !% --- return depth to physical value
        outvalue = outvalue * uniformTableR(UT_idx,utr_DepthMax)

        ! print *, 'outvalue 2 ',outvalue

    end function geo_criticaldepth_singular
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
            !print *, 'here AAAA '
            eIdx      => uniformTableI(UT_idx,uti_elem_idx)
            !print *, 'here BBBB '
            thisTable => uniformTableDataR(UT_idx,:,utd_SF_depth_nonuniform)  !% element index
        !%------------------------------------------------------------------
        !% --- section factor for the associated element
        sectionFactor = elemR(eIdx,er_Flowrate) * elemR(eIdx,er_ManningsN) / elemR(eIdx,er_BottomSlope)

        ! print *, 'sectionFactor ',sectionFactor
        ! print *, 'flowrate   ',elemR(eIdx,er_Flowrate)
        ! print *, 'mannings n ',elemR(eIdx,er_ManningsN)
        ! print *, 'slope      ',elemR(eIdx,er_BottomSlope)

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

        !print *, 'normSF ',normSF

        !% --- lookup the normalized normal depth
        outvalue = xsect_table_lookup_singular(normSF,thisTable)

        ! print *, 'outvalue 1 ',outvalue

        !% --- return normal depth to physical value
        outvalue = outvalue * uniformTableR(UT_idx,utr_DepthMax)

        ! print *, 'outvalue 2 ',outvalue
    
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
    subroutine geo_depth_from_volume_by_type (elemPGx, npack_elemPGx, col_elemPGx)
        !%------------------------------------------------------------------
        !% Description:
        !% This solves nonsurcharged CCJMJB elements because of PGx arrays
        !% The elemPGx determines whether this is ALLtm, ETM or AC elements
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: elemPGx(:,:), npack_elemPGx(:), col_elemPGx(:)
            integer, pointer :: Npack, thisCol
            character(64) :: subroutine_name = 'geo_depth_from_volume_by_type'
        !%-------------------------------------------------------------------
        !% Preliminaries
            !!if (crashYN) return
            if (setting%Debug%File%geometry) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------    
        !% cycle through different geometries  
                
        !% --- open channels ------------------------------------------

        !% --- IRREGULAR
        thisCol => col_elemPGx(epg_CC_irregular)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call irregular_depth_from_volume (elemPGx, Npack, thisCol)
        end if        
                
        !% -- POWER FUNCTION
        thisCol => col_elemPGx(epg_CC_power_function)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            print *, 'POWER FUNCTION CROSS-SECTION NOT COMPLETE'
            call util_crashpoint(5559872)
            !call powerfunction_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% -- PARABOLIC
        thisCol => col_elemPGx(epg_CC_parabolic)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call parabolic_depth_from_volume (elemPGx, Npack, thisCol)
        end if
                
        !% --- RECTANGULAR CHANNEL
        thisCol => col_elemPGx(epg_CC_rectangular)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call rectangular_depth_from_volume (elemPGx, Npack, thisCol)
        end if    

        !% --- TRAPEZOIDAL CHANNEL
        thisCol => col_elemPGx(epg_CC_trapezoidal)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call trapezoidal_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% --- TRIANGULAR CHANNEL
        thisCol => col_elemPGx(epg_CC_triangular)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call triangular_depth_from_volume (elemPGx, Npack, thisCol)
        end if


        !% --- CLOSED CONDUITS  ---------------------------------------

        !% --  ARCH CONDUIT
        thisCol => col_elemPGx(epg_CC_arch)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume   &
                (elemPGx, Npack, thisCol, YArch)
        end if

        !% --  BASKET_HANDLE
        thisCol => col_elemPGx(epg_CC_basket_handle)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume           &
                (elemPGx, Npack, thisCol, YBasketHandle)
        end if

        !% --  CATENARY
        thisCol => col_elemPGx(epg_CC_catenary)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume           &
                (elemPGx, Npack, thisCol, YCatenary)
        end if

        !% --- CIRCULAR CONDUIT
        thisCol => col_elemPGx(epg_CC_circular)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call circular_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% --  EGG_SHAPED
        thisCol => col_elemPGx(epg_CC_egg_shaped)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume           &
                (elemPGx, Npack, thisCol, YEgg)
        end if

        !% --- FILLED CIRCULAR CONDUIT
        thisCol => col_elemPGx(epg_CC_filled_circular)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call filled_circular_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% --  GOTHIC
        thisCol => col_elemPGx(epg_CC_gothic)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume           &
                (elemPGx, Npack, thisCol, YGothic)
        end if

        !% --  HORIZONTAL ELLIPSE
        thisCol => col_elemPGx(epg_CC_horiz_ellipse)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume           &
                (elemPGx, Npack, thisCol, YHorizEllip)
        end if

        !% --  HORSESHOE
        thisCol => col_elemPGx(epg_CC_horse_shoe)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume           &
                (elemPGx, Npack, thisCol, YHorseShoe)
        end if

        !% --  MODIFIED BASKET HANDLE
        thisCol => col_elemPGx(epg_CC_mod_basket)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call mod_basket_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% --- RECTANGULAR CLOSED CONDUIT
        thisCol => col_elemPGx(epg_CC_rectangular_closed)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call rectangular_closed_depth_from_volume (elemPGx, Npack, thisCol)
        end if        

        !% --- RECTANGULAR ROUND CONDUIT
        thisCol => col_elemPGx(epg_CC_rectangular_round)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call rect_round_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% -- RECTANGULAR TRIANGULAR
        thisCol => col_elemPGx(epg_CC_rectangular_triangular)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call rectangular_triangular_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        ! % --- SEMI-CIRCULAR CONDUIT
        thisCol => col_elemPGx(epg_CC_semi_circular)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume           &
                (elemPGx, Npack, thisCol, YSemiCircular)
        end if

        !% --  SEMI ELLIPTICAL
        thisCol => col_elemPGx(epg_CC_semi_elliptical)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume           &
                (elemPGx, Npack, thisCol, YSemiEllip)
        end if

        !% --  VERTICAL ELLIPSE
        thisCol => col_elemPGx(epg_CC_vert_ellipse)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_depth_from_volume           &
                (elemPGx, Npack, thisCol, YVertEllip)
        end if

        
        !% --- JUNCTIONS ---------------------------------------------------- 

        !% JM with functional geometry
        thisCol => col_elemPGx(epg_JM_functionalStorage)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call storage_functional_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% JM with tabular geometry
        thisCol => col_elemPGx(epg_JM_tabularStorage)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call storage_tabular_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% JM with implied storage 
        thisCol => col_elemPGx(epg_JM_impliedStorage)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call storage_implied_depth_from_volume (elemPGx, Npack, thisCol)
        end if
  
        !%-------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%geometry) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_depth_from_volume_by_type
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
!%==========================================================================
!%
    subroutine geo_overflow_openchannels (thisColP_Open_CC)
        !%------------------------------------------------------------------
        !% Description:
        !% Adds the overflow of open channels that exceed their 
        !% Full Volume to the overflow accumulator for this step
        !% Reduces stored volume by the overflow amount.
        !% Note that open channels CANNOT pond in present version.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisColP_Open_CC !% must be open channels
            integer, pointer    :: Npack, thisP(:), nBarrels(:)
            real(8), pointer    :: volume(:), fullvolume(:), overflow(:)
            real(8), pointer    :: temp(:)
        !%------------------------------------------------------------------
        !% Preliminaries
            Npack      => npack_elemP(thisColP_Open_CC)
            if (Npack < 1) return
            if (.not. setting%Discretization%AllowChannelOverflowTF) return
        !%------------------------------------------------------------------
        !% Aliases
            thisP      => elemP(1:Npack,thisColP_Open_CC)
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
    subroutine geo_volumelimit_closed (thisColP_Closed) 
        !%------------------------------------------------------------------
        !% Description
        !% Sets closed element volumes to the full volume
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisColP_Closed
            integer, pointer    :: Npack, thisP(:)
            real(8), pointer    :: volume(:), fullvolume(:)
        !%------------------------------------------------------------------
        !% Preliminaries
            Npack      => npack_elemP(thisColP_Closed)
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases:
            thisP      => elemP(1:Npack,thisColP_Closed)
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
    subroutine geo_assign_JB (whichTM, thisColP_JM)
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
            integer, intent(in) :: whichTM, thisColP_JM

            integer, pointer ::  Npack, thisP(:), BranchExists(:), thisSolve(:),  tM
            real(8), pointer :: area(:), depth(:), head(:), hyddepth(:), hydradius(:)
            real(8), pointer :: length(:), perimeter(:), topwidth(:), velocity(:)
            real(8), pointer :: volume(:), zBtm(:), Kfac(:), dHdA(:), ell(:), ellMax(:)
            real(8), pointer :: zCrown(:), fullArea(:), fulldepth(:), fullperimeter(:)
            real(8), pointer :: fullhyddepth(:), sedimentDepth(:), thisTable(:,:)
            real(8), pointer :: fulltopwidth(:)
            !real(8), pointer :: pressurehead(:)
            real(8), pointer :: slotDepth(:), slotVolume(:), overflow(:), fullhydradius(:)
            real(8), pointer :: Atable(:), Ttable(:), Rtable(:), Stable(:)
            real(8), pointer :: grav  
            logical, pointer :: isSlot(:)     

            real(8) :: depthnorm, zeroHydRadius
            integer :: tB, ii, kk, tBA(1)


        !% thisColP_JM is the column for the junction mains of a particular
        !% whichTM. For ALL ep_JM, for ETM, ep_JM_ETM, for AC ep_JM_AC
            integer, allocatable :: tempP(:)
            character(64) :: subroutine_name = 'geo_assign_JB'
        !%---------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%geometry) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            if (setting%Profile%useYN) call util_profiler_start (pfc_geo_assign_JB)
        !%----------------------------------------------------------------------
        !% Aliases
            Npack         => npack_elemP(thisColP_JM)
            area          => elemR(:,er_Area)
            depth         => elemR(:,er_Depth)
            dHdA          => elemR(:,er_dHdA)
            ell           => elemR(:,er_ell)
            ellMax        => elemR(:,er_FullEll)
            head          => elemR(:,er_Head)
            hyddepth      => elemR(:,er_HydDepth)
            hydradius     => elemR(:,er_HydRadius)
            length        => elemR(:,er_Length)
            perimeter     => elemR(:,er_Perimeter)
            !pressurehead  => elemR(:,er_Pressure_Head)
            sedimentDepth => elemR(:,er_SedimentDepth)
            topwidth      => elemR(:,er_Topwidth)
            velocity      => elemR(:,er_Velocity)
            volume        => elemR(:,er_Volume)
            zBtm          => elemR(:,er_Zbottom)
            zCrown        => elemR(:,er_Zcrown)
            fullArea      => elemR(:,er_FullArea)
            fulldepth     => elemR(:,er_FullDepth)
            fullTopWidth  => elemR(:,er_FullTopWidth)
            fullhyddepth  => elemR(:,er_FullHydDepth)
            fullhydradius => elemR(:,er_FullHydRadius)
            fullperimeter => elemR(:,er_FullPerimeter)
            !overflow      => elemR(:,er_VolumeOverFlow)
            !slotDepth     => elemR(:,er_SlotDepth)
            !slotVolume    => elemR(:,er_SlotVolume)
            Kfac          => elemSR(:,esr_JunctionBranch_Kfactor)
            BranchExists  => elemSI(:,esi_JunctionBranch_Exists)
            thisSolve     => elemI(:,ei_tmType)
            !isSlot        => elemYN(:,eYN_isPSsurcharged)
            grav => setting%Constant%gravity
        !%------------------------------------------------------------------

        if (Npack > 0) then
            thisP  => elemP(1:Npack,thisColP_JM)

            !% cycle through the all the main junctions and each of its branches
            do ii=1,Npack
                
                tM => thisP(ii) !% junction main ID

                !% moved 20220909brh
                ! !% if a slot present, add the slot depth and volume back to JM
                ! if (isSlot(tM)) then
                !     volume(tM)   = volume(tM)  + SlotVolume(tM) 
                !     depth(tM)    = depth(tM)   + SlotDepth(tM)
                !     head(tM)     = head(tM)    + SlotDepth(tM)
                !     ell(tM)      = ellMax(tM)
                !     !% Overflow(tM) = zeroR
                ! end if 

                !% only execute for whichTM of ALL or thisSolve (of JM) matching input whichTM
                if ((whichTM == ALLtm) .or. (thisSolve(tM) == whichTM)) then
                    !% cycle through the possible junction branches
                    do kk=1,max_branch_per_node
                        
                        tB = tM + kk !% junction branch ID
                        tBA(1) = tB  !% array for pure array functions

                        ! print *, kk, tB
                        ! print *, BranchExists(tB)

                        if (BranchExists(tB) == 1) then
                            !% only when a branch exists.
                            ! print *, head(tM), zBtm(tB)
                            ! print *, kk, branchsign(kk)
                            ! print *, velocity(tB)
                            ! print *, Kfac(tB)
                            if ( head(tM) > (zBtm(tB) + sedimentDepth(tB)) ) then
                                !% for main head above branch bottom entrance use a head
                                !% loss approach. The branchsign and velocity sign ensure
                                !% the headloss is added to an inflow and subtracted at
                                !% an outflow
                                !% Note this is a time-lagged velocity as the JB velocity
                                !% is not updated until after face interpolation                                
                                head(tB) = head(tM) + sedimentDepth(tB)                &
                                    + branchsign(kk) * sign(oneR,velocity(tB))         &
                                    * (Kfac(tB) / (twoR * grav)) * (velocity(tB)**twoR)
                               
                            else
                                !% for main head below the branch bottom entrance we assign a
                                !% Froude number of one on an inflow to the junction main. Note
                                !% an outflow from a junction main for this case gets head
                                !% of z_bottom of the branch (zero depth).
                                !% Note this is a time-lagged velocity as the JB velocity
                                !% is not updated until after face interpolation
                                head(tB) = zBtm(tB) + sedimentDepth(tB)                            &
                                    + onehalfR * (oneR + branchsign(kk) * sign(oneR,velocity(tB))) &
                                    *(velocity(tB)**twoR) / (grav) 
                            end if

                            !% HACK -- the above uses a Froude number argument for head(TM) < zBtm(tB)
                            !%      however, when the JB is surcharged we probably should be using the
                            !%      K factor approach and require K=1.
                           
                            !% compute provisional depth
                            depth(tB) = head(tB) - (zBtm(tB) + sedimentDepth(tB))

                            ! print *, 'in geo_assign_JB  ',trim(reverseKey(elemI(tB,ei_geometryType)))
                            ! print *, 'depth ',depth(tB), fulldepth(tB), setting%ZeroValue%Depth
                            
                            if (depth(tB) .ge. fulldepth(tB)) then
                                !% surcharged or incipient surcharged
                                depth(tB)        = fulldepth(tB)
                                area(tB)         = fullArea(tB)
                                hyddepth(tB)     = fullhyddepth(tB)
                                perimeter(tB)    = fullperimeter(tB)
                                topwidth(tB)     = setting%ZeroValue%Topwidth
                                hydRadius(tB)    = fulldepth(tB) / fullperimeter(tB)
                                dHdA(tB)         = oneR / setting%ZeroValue%Topwidth
                                ell(tBA)          = llgeo_ell_pure(tBA)
                                !pressurehead(tB) = zBtm(tB) + ell(tB)

                                ! write(*,"(A,i5,10f12.5)") 'AAA ell ',tB, ell(tB), depth(tB), hydDepth(tB), fulldepth(tB)

                            elseif ((depth(tB) < setting%ZeroValue%Depth) .and. (setting%ZeroValue%UseZeroValues)) then
                                !% negligible depth is treated with ZeroValues
                                depth(tB)        = setting%ZeroValue%Depth
                                area(tB)         = setting%ZeroValue%Area
                                topwidth(tB)     = setting%ZeroValue%Topwidth
                                hyddepth(tB)     = setting%ZeroValue%Depth !% setting%ZeroValue%Area / topwidth(tB) 20220712brh
                                perimeter(tB)    = topwidth(tB) + setting%ZeroValue%Depth
                                hydRadius(tB)    = setting%ZeroValue%Area / perimeter(tB)
                                dHdA(tB)         = oneR / topwidth(tB)
                                ell(tB)          = setting%ZeroValue%Depth !%hydDepth(tB)  20220712 brh
                                !pressurehead(tB) = zBtm(tB) + ell(tB)

                                ! write(*,"(A,i5,10f12.5)"), 'BBB ell ',tB, ell(tB), depth(tB), hydDepth(tB), fulldepth(tB)

                            elseif ((depth(tB) .le. zeroR) .and. (.not. setting%ZeroValue%UseZeroValues)) then
                                !% negative depth without zero value treatment (not recommended!) is treated as exactly zero
                                depth(tB) = zeroR
                                area(tB)  = zeroR
                                topwidth(tB) = zeroR
                                hydDepth(tB) = zeroR
                                perimeter(tB) = zeroR
                                hydRadius(tB) = zeroR
                                dHdA(tB)      = oneR / setting%ZeroValue%Topwidth
                                ell(tB)       = zeroR
                                !pressurehead(tB) = zBtm(tB) + ell(tB)

                                ! write(*,"(A,i5,10f12.5)") 'CCC ell ',tB, ell(tB), depth(tB), hydDepth(tB), fulldepth(tB)

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
                                    area(tBA)     = llgeo_parabolic_area_from_depth_pure         (tBA, depth(tBA))
                                    topwidth(tBA) = llgeo_parabolic_topwidth_from_depth_pure     (tBA, depth(tBA))
                                    perimeter(tBA)= llgeo_parabolic_perimeter_from_depth_pure    (tBA, depth(tBA))
                                    hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                        (tBA, area(tBA), perimeter(tBA))
                                    ell(tBA)      = llgeo_hyddepth_from_area_and_topwidth_pure   &
                                                        (tBA, area(tBA), topwidth(tBA))
 
                                case (power_function) !% POSSIBLY LOOKUP
                                    print *, 'CODE ERROR power function x-section not finished'
                                    call util_crashpoint(6298349)

                                case (rectangular) !% analytical
                                    area(tBA)     = llgeo_rectangular_area_from_depth_pure       (tBA, depth(tBA))
                                    topwidth(tBA) = llgeo_rectangular_topwidth_from_depth_pure   (tBA, depth(tBA))
                                    perimeter(tBA)= llgeo_rectangular_perimeter_from_depth_pure  (tBA, depth(tBA))
                                    hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                        (tBA, area(tBA), perimeter(tBA))
                                    ell(tBA)      = llgeo_hyddepth_from_area_and_topwidth_pure  &
                                                        (tBA, area(tBA), topwidth(tBA))

                                case (trapezoidal) !% analytical
                                    area(tBA)     = llgeo_trapezoidal_area_from_depth_pure       (tBA, depth(tBA))
                                    topwidth(tBA) = llgeo_trapezoidal_topwidth_from_depth_pure   (tBA, depth(tBA))
                                    perimeter(tBA)= llgeo_trapezoidal_perimeter_from_depth_pure  (tBA, depth(tBA))
                                    hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                        (tBA, area(tBA), perimeter(tBA))
                                    ell(tBA)      = llgeo_hyddepth_from_area_and_topwidth_pure   &
                                                        (tBA, area(tBA), topwidth(tBA))

                                case (triangular) !% analytical
                                    area(tBA)     = llgeo_triangular_area_from_depth_pure        (tBA, depth(tBA))
                                    topwidth(tBA) = llgeo_triangular_topwidth_from_depth_pure    (tBA, depth(tBA))
                                    perimeter(tBA)= llgeo_triangular_perimeter_from_depth_pure   (tBA, depth(tBA))
                                    hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                         (tBA, area(tBA), perimeter(tBA))
                                    ell(tBA)      = llgeo_hyddepth_from_area_and_topwidth_pure   &
                                                        (tBA, area(tBA), topwidth(tBA))

                                case (irregular)  !% lookup
                                    area(tB)     = irregular_geometry_from_depth_singular ( &
                                                        tB,tt_area, depth(tB), setting%ZeroValue%Depth)

                                    topwidth(tB) = irregular_geometry_from_depth_singular ( &
                                                        tB,tt_width, depth(tB), setting%ZeroValue%TopWidth)

                                    !% --- note the irregular stores hyd radius rather than perimeter
                                    zeroHydRadius = setting%ZeroValue%Area / (setting%ZeroValue%TopWidth + setting%ZeroValue%Depth)
                                    hydRadius(tB) = irregular_geometry_from_depth_singular ( &
                                                        tB,tt_hydradius, depth(tB), zeroHydRadius)                
                                    !% --- perimeter is derived geometry for irregular
                                    perimeter(tB) = area(tB) / hydRadius(tB)
                                    !% --- irregular must be continuously-increasing in width
                                    ell(tB)       = geo_hyddepth_from_area_and_topwidth_singular (tB, area(tB), topwidth(tB)) 

                                !% --- CLOSED CONDUITS
                                !%     closed conduits typically have look-up functions for area, topwidth and hydraulic
                                !%     radius, with standard geo_functions for perimeter, hydraulic depth, and ell
                                !%     However, where analytical functions are used, the perimeter is usually computed
                                !%     first and hydraulic radius is a geo_ function

                                !% --- lookups with Hydraulic Radius stored
                                case (arch, basket_handle, circular, eggshaped, horiz_ellipse, horseshoe, vert_ellipse)                                    
                                    area(tB)     = llgeo_tabular_from_depth_singular &
                                        (tB, depth(tB), fullArea(tB), setting%ZeroValue%Depth, Atable)

                                    topwidth(tB) = llgeo_tabular_from_depth_singular &
                                        (tB, depth(tB), fullTopWidth(tB), setting%ZeroValue%Depth, Ttable)

                                    hydRadius(tB)= llgeo_tabular_from_depth_singular &
                                        (tB, depth(tB), fullHydRadius(tB), setting%ZeroValue%Depth, Rtable)

                                    perimeter(tBA)= llgeo_perimeter_from_hydradius_and_area_pure &
                                                        (tBA, hydradius(tBA), area(tBA))
                                    ell(tBA)      = llgeo_ell_pure (tBA) 

                                !% --- lookups with SectionFactor stored
                                case (catenary, gothic, semi_circular, semi_elliptical)    
                                    area(tB)     = llgeo_tabular_from_depth_singular &
                                         (tB, depth(tB), fullArea(tB), setting%ZeroValue%Depth, Atable)

                                    topwidth(tB) = llgeo_tabular_from_depth_singular &
                                         (tB, depth(tB), fullTopWidth(tB), setting%ZeroValue%Depth, Ttable)

                                    hydRadius(tB)= llgeo_tabular_hydradius_from_area_and_sectionfactor_singular &
                                         (tB, area(tB), fullhydradius(tB), setting%ZeroValue%Area, Stable)

                                    perimeter(tBA)= llgeo_perimeter_from_hydradius_and_area_pure &
                                                        (tBA, hydradius(tBA), area(tBA))
                                    ell(tBA)      = llgeo_ell_pure (tBA) 

                                !% --- lookup with sediment
                                case (filled_circular)  
                                    area(tB)      = llgeo_filled_circular_area_from_depth_singular      (tB,depth(tB))
                                    topwidth(tB)  = llgeo_filled_circular_topwidth_from_depth_singular  (tB,depth(tB))
                                    perimeter(tB) = llgeo_filled_circular_perimeter_from_depth_singular (tB,depth(tB))
                                    hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                        (tBA, area(tBA), perimeter(tBA))
                                    ell(tBA)      = llgeo_ell_pure (tBA) 

                                !% --- analytical closed-conduit cases
                                case (mod_basket)   !% analytical                                 
                                    area(tB)      = llgeo_mod_basket_area_from_depth_singular        (tB,depth(tB))
                                    topwidth(tB)  = llgeo_mod_basket_topwidth_from_depth_singular    (tB,depth(tB))
                                    perimeter(tB) = llgeo_mod_basket_perimeter_from_depth_singular   (tB,depth(tB))
                                    hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                        (tBA, area(tBA), perimeter(tBA))
                                    ell(tBA)      = llgeo_ell_pure (tBA) 

                                case (rectangular_closed) !% analytical
                                    area(tB)      = llgeo_rectangular_closed_area_from_depth_singular      (tB, depth(tB))
                                    topwidth(tB)  = llgeo_rectangular_closed_topwidth_from_depth_singular  (tB, depth(tB))
                                    perimeter(tB) = llgeo_rectangular_closed_perimeter_from_depth_singular (tB, depth(tB))
                                    hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                        (tBA, area(tBA), perimeter(tBA))
                                    ell(tBA)      = llgeo_ell_pure (tBA) 

                                case (rect_round)  !% analytical                                  
                                    area(tB)      = llgeo_rect_round_area_from_depth_singular       (tB,depth(tB))
                                    topwidth(tB)  = llgeo_rect_round_topwidth_from_depth_singular   (tB,depth(tB))
                                    perimeter(tB) = llgeo_rect_round_perimeter_from_depth_singular  (tB,depth(tB))
                                    hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure&
                                                        (tBA, area(tBA), perimeter(tBA))
                                    ell(tBA)      = llgeo_ell_pure (tBA) 

                                case (rect_triang) !% analytical                                   
                                    area(tB)      = llgeo_rectangular_triangular_area_from_depth_singular      (tB,depth(tB))
                                    topwidth(tB)  = llgeo_rectangular_triangular_topwidth_from_depth_singular  (tB,depth(tB))
                                    perimeter(tB) = llgeo_rectangular_triangular_perimeter_from_depth_singular (tB,depth(tB))
                                    hydRadius(tBA)= llgeo_hydradius_from_area_and_perimeter_pure &
                                                        (tBA, area(tBA), perimeter(tBA))
                                    ell(tBA)      = llgeo_ell_pure (tBA) 
                            
                                case default
                                    print *, 'CODE ERROR: geometry type unknown for # ', elemI(tB,ei_geometryType)
                                    print *, 'which has key ',trim(reverseKey(elemI(tB,ei_geometryType)))
                                    print *, 'in ',trim(subroutine_name)
                                    call util_crashpoint(399848)
                                    !return
                                    !stop 399848
                                end select
                                !% --- standard for all geometries
                                hydDepth(tBA) = llgeo_hyddepth_from_area_and_topwidth_pure &
                                                    (tBA, area(tBA), topwidth(tBA))
                                dHdA(tB)     = oneR / topwidth(tB)
                            end if

                            !% --- universal computation of volume
                            volume(tB) = area(tB) * length(tB)
                        end if
                    end do
                end if
            end do
        end if

        !% Note, the above can only be made a concurrent loop if we replace the tM
        !% with thisP(ii) and tB with thisP(ii)+kk, which makes the code
        !% difficult to read.

        if (setting%Profile%useYN) call util_profiler_stop (pfc_geo_assign_JB)

        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_assign_JB
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
    subroutine geo_topwidth_from_depth_by_type &
        (elemPGx, npack_elemPGx, col_elemPGx)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth given depth of a non-surcharged element
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: elemPGx(:,:)
            integer, target, intent(in) :: npack_elemPGx(:), col_elemPGx(:)
            integer, pointer :: Npack, thisCol

            character(64) :: subroutine_name = 'geo_topwidth_from_depth_by_type'
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
            call irregular_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- PARABOLIC
        Npack => npack_elemPGx(epg_CC_parabolic)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_parabolic)
            call parabolic_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- POWERFUNCTION
        Npack => npack_elemPGx(epg_CC_power_function)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_power_function)
            call powerfunction_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- RECTANGULAR
        Npack => npack_elemPGx(epg_CC_rectangular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular)
            call rectangular_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- TRAPEZOIDAL
        Npack => npack_elemPGx(epg_CC_trapezoidal)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_trapezoidal)
            call trapezoidal_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- TRIANGULAR
        Npack => npack_elemPGx(epg_CC_triangular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_triangular)
            call triangular_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if



        !% --- CLOSED CONDUITS -----------------------------------------------

        !% --  ARCH
        thisCol => col_elemPGx(epg_CC_arch)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, TArch)
        end if

        !% -- BASKET_HANDLE
        Npack => npack_elemPGx(epg_CC_basket_handle)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_basket_handle)
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, TBasketHandle)
        end if

        !% -- CATENARY
        Npack => npack_elemPGx(epg_CC_catenary)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_catenary)
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, TCatenary)
        end if

        !% --- CIRCULAR
        Npack => npack_elemPGx(epg_CC_circular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_circular)
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, TCirc)
        end if

        !% -- EGG_SHAPED
        Npack => npack_elemPGx(epg_CC_egg_shaped)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_egg_shaped)
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, TEgg)
        end if

        !% --- FILLED CIRCULAR
        Npack => npack_elemPGx(epg_CC_filled_circular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_filled_circular)
            call filled_circular_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- GOTHIC
        Npack => npack_elemPGx(epg_CC_gothic)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_gothic)
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, TGothic)
        end if

        !% --  HORIZONTAL ELLIPSE
        thisCol => col_elemPGx(epg_CC_horiz_ellipse)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, THorizEllip)
        end if

        !% -- HORSE_SHOE
        Npack => npack_elemPGx(epg_CC_horse_shoe)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_horse_shoe)
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, THorseShoe)
        end if

        !% -- MODIFIED BASKET HANDLE
        Npack => npack_elemPGx(epg_CC_mod_basket)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_mod_basket)
            call mod_basket_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- RECTANGULAR CLOSED
        Npack => npack_elemPGx(epg_CC_rectangular_closed)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_closed)
            call rectangular_closed_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --  RECTANGULAR ROUND
        Npack   => npack_elemPGx(epg_CC_rectangular_round)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_round)
            call rect_round_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if
        
        !% -- RECTANGULAR TRIANGULAR
        Npack => npack_elemPGx(epg_CC_rectangular_triangular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_triangular)
            call rectangular_triangular_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- SEMI-CIRCULAR
        Npack => npack_elemPGx(epg_CC_semi_circular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_semi_circular)
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, TSemiCircular)
        end if

        !% -- SEMI-ELLIPTICAL
        Npack => npack_elemPGx(epg_CC_semi_elliptical)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_semi_elliptical)
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, TSemiEllip)
        end if        

        !% --  VERTICAL ELLIPSE
        thisCol => col_elemPGx(epg_CC_vert_ellipse)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call llgeo_tabular_topwidth_from_depth  &
                (elemPGx, Npack, thisCol, TVertEllip)
        end if

        !% TOPWIDTH NOT DEFINED FOR TABULAR, FUNCTIONAL, IMPLIED STORAGE

        !%-------------------------------------------------------------------
            if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_topwidth_from_depth_by_type
!%
!%========================================================================== 
!%==========================================================================  
!%
    subroutine geo_perimeter_and_hydradius_from_depth_by_type &
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
        Npack => npack_elemPGx(epg_CC_irregular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_irregular)
            !% --- solve both perimeter and hydraulic radius
            call irregular_perimeter_and_hydradius_from_depth (elemPGx, Npack, thisCol)
            !call irregular_perimeter_from_hydradius_area (elemPGx, Npack, thisCol)
        end if

        !% --- PARABOLIC
        Npack => npack_elemPGx(epg_CC_parabolic)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_parabolic)
            thisP   => elemPGx(1:Npack,thisCol)
            call parabolic_perimeter_from_depth (elemPGx, Npack, thisCol)
            elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure &
                                        (thisP, area(thisP), perimeter(thisP))
        end if
                
        !% --- POWER FUNCTION
        Npack => npack_elemPGx(epg_CC_power_function)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_power_function)
            thisP   => elemPGx(1:Npack,thisCol)
            print *, 'POWER FUNCTION CROSS SECTION NOT COMPLETED'
            call util_crashpoint(54987)
            call powerfunction_perimeter_from_depth (elemPGx, Npack, thisCol)
            elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure &
                                         (thisP, area(thisP), perimeter(thisP))
        end if

        !% --- RECTANGULAR
        Npack => npack_elemPGx(epg_CC_rectangular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular)
            thisP   => elemPGx(1:Npack,thisCol)
            call rectangular_perimeter_from_depth (elemPGx, Npack, thisCol)
            elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure &
                                        (thisP, area(thisP), perimeter(thisP))
        end if

        !% --- TRAPEZOIDAL
        Npack => npack_elemPGx(epg_CC_trapezoidal)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_trapezoidal)
            thisP   => elemPGx(1:Npack,thisCol)
            call trapezoidal_perimeter_from_depth (elemPGx, Npack, thisCol)
            elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure &
                                        (thisP, area(thisP), perimeter(thisP))
        end if

        !% --- TRIANGULAR
        Npack => npack_elemPGx(epg_CC_triangular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_triangular)
            thisP   => elemPGx(1:Npack,thisCol)
            call triangular_perimeter_from_depth (elemPGx, Npack, thisCol)
            elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure &
                                        (thisP, area(thisP), perimeter(thisP))
        end if


        
        !% --- CLOSED CONDUITS ---------------------------------

        !% --- ARCH
        Npack => npack_elemPGx(epg_CC_arch)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_arch)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_depth  &
                (elemPGx, Npack, thisCol, RArch)
            elemR(thisP,er_Perimeter) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            !call arch_perimeter_from_depth (elemPGx, Npack, thisCol)
            !elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure(thisP)
        end if

        !% -- BASKET_HANDLE
        Npack   => npack_elemPGx(epg_CC_basket_handle)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_basket_handle)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_depth  &
                (elemPGx, Npack, thisCol, RBasketHandle)
            elemR(thisP,er_Perimeter) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            !call basket_handle_perimeter_from_depth (elemPGx, Npack, thisCol)
            !elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure(thisP)
        end if

        !% -- CATENARY
        Npack => npack_elemPGx(epg_CC_catenary)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_catenary)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_sectionfactor_and_area &
                (elemPGx, Npack, thisCol, SCatenary, esgr_Catenary_SoverSfull)
            elemR(thisP,er_Perimeter) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            !call catenary_perimeter_from_depth (elemPGx, Npack, thisCol)
            !elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure(thisP)
        end if

        !% --- CIRCULAR
        Npack => npack_elemPGx(epg_CC_circular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_circular)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_depth  &
                (elemPGx, Npack, thisCol, RCirc)
            elemR(thisP,er_Perimeter) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            !call circular_perimeter_from_depth (elemPGx, Npack, thisCol)
            !elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure(thisP)
        end if

        !% -- EGG_SHAPED
        Npack => npack_elemPGx(epg_CC_egg_shaped)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_egg_shaped)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_depth  &
                (elemPGx, Npack, thisCol, REgg)
            elemR(thisP,er_Perimeter) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            !call egg_shaped_perimeter_from_depth (elemPGx, Npack, thisCol)
            !elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure(thisP)
        end if

        !% --- FILLED CIRCULAR
        Npack => npack_elemPGx(epg_CC_filled_circular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_filled_circular)
            thisP   => elemPGx(1:Npack,thisCol)
            call filled_circular_hydradius_and_perimeter_from_depth (elemPGx, Npack, thisCol)
            !call filled_circular_perimeter_from_depth (elemPGx, Npack, thisCol)
            !elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure(thisP)
        end if

        !% -- GOTHIC
        Npack => npack_elemPGx(epg_CC_gothic)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_gothic)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_sectionfactor_and_area &
                (elemPGx, Npack, thisCol, SGothic, esgr_Gothic_SoverSfull)
            elemR(thisP,er_Perimeter) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            !call gothic_perimeter_from_depth (elemPGx, Npack, thisCol)
            !elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure(thisP)
        end if

        !% -- HORIZONTAL ELLIPSE
        Npack => npack_elemPGx(epg_CC_horiz_ellipse)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_horiz_ellipse)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_depth  &
                (elemPGx, Npack, thisCol, RHorizEllip)
            elemR(thisP,er_Perimeter) = llgeo_perimeter_from_hydradius_and_area_pure &
                                           (thisP, hydradius(thisP), area(thisP))
            !call horiz_ellipse_perimeter_from_depth (elemPGx, Npack, thisCol)
            !elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure(thisP)
        end if

        !% -- HORSE_SHOE
        Npack => npack_elemPGx(epg_CC_horse_shoe)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_horse_shoe)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_depth  &
                (elemPGx, Npack, thisCol, RHorseShoe)
            elemR(thisP,er_Perimeter) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            !call horse_shoe_perimeter_from_depth (elemPGx, Npack, thisCol)
            !elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure(thisP)
        end if

        !% -- MODIFIED BASKET HANDLE
        Npack => npack_elemPGx(epg_CC_mod_basket)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_mod_basket)
            thisP   => elemPGx(1:Npack,thisCol)
            call mod_basket_perimeter_from_depth (elemPGx, Npack, thisCol)
            elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure &
                                         (thisP, area(thisP), perimeter(thisP))
        end if

        !% --- RECTANGULAR CLOSED
        Npack => npack_elemPGx(epg_CC_rectangular_closed)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_closed)
            thisP   => elemPGx(1:Npack,thisCol)
            call rectangular_closed_perimeter_from_depth (elemPGx, Npack, thisCol)
            elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure &
                                         (thisP, area(thisP), perimeter(thisP))
        end if

        !% --  RECTANGULAR ROUND
        Npack   => npack_elemPGx(epg_CC_rectangular_round)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_round)
            thisP   => elemPGx(1:Npack,thisCol)
            call rect_round_perimeter_from_depth (elemPGx, Npack, thisCol)
            elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure &
                                         (thisP, area(thisP), perimeter(thisP))
        end if

        !% -- RECTANGULAR TRIANGULAR
        Npack => npack_elemPGx(epg_CC_rectangular_triangular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_triangular)
            thisP   => elemPGx(1:Npack,thisCol)
            call rectangular_triangular_perimeter_from_depth (elemPGx, Npack, thisCol)
            elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure &
                                        (thisP, area(thisP), perimeter(thisP))
        end if

        !% -- SEMI-CIRCULAR
        Npack => npack_elemPGx(epg_CC_semi_circular)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_semi_circular)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_sectionfactor_and_area &
                (elemPGx, Npack, thisCol, SSemiCircular, esgr_Semi_Circular_SoverSfull)
            elemR(thisP,er_Perimeter) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            !call semi_circular_perimeter_from_depth (elemPGx, Npack, thisCol)
            !elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure(thisP)
        end if

        !% -- SEMI-ELLIPTICAL
        Npack => npack_elemPGx(epg_CC_semi_elliptical)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_semi_elliptical)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_sectionfactor_and_area &
                (elemPGx, Npack, thisCol, SSemiEllip, esgr_Semi_Elliptical_SoverSfull)
            elemR(thisP,er_Perimeter) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            !call semi_elliptical_perimeter_from_depth (elemPGx, Npack, thisCol)
            !elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure(thisP)
        end if

        !% -- VERTICAL ELLIPSE
        Npack   => npack_elemPGx(epg_CC_vert_ellipse)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_vert_ellipse)
            thisP   => elemPGx(1:Npack,thisCol)
            call llgeo_tabular_hydradius_from_depth  &
                (elemPGx, Npack, thisCol, RVertEllip)
            elemR(thisP,er_Perimeter) = llgeo_perimeter_from_hydradius_and_area_pure &
                                            (thisP, hydradius(thisP), area(thisP))
            !call vert_ellipse_perimeter_from_depth (elemPGx, Npack, thisCol)
            !elemR(thisP,er_HydRadius) = llgeo_hydradius_from_area_and_perimeter_pure(thisP)
        end if

        !% TOPWIDTH is UNDEFINED FOR TABULAR, FUNCTIONAL, AND IMPLIED STORAGE

        !%-------------------------------------------------------------------
        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_perimeter_and_hydradius_from_depth_by_type
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
    subroutine geo_ell_from_head (thisColP)
        !%------------------------------------------------------------------
        !% Description:
        !% computes the value of "ell" -- the modified hydraulic depth
        !% used as a length scale in AC method
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisColP
            integer, pointer :: thisP(:), Npack
            real(8), pointer :: ell(:), head(:), area(:), topwidth(:), hydDepth(:)
            real(8), pointer :: ZbreadthMax(:), breadthMax(:), areaBelowBreadthMax(:)
            integer :: ii

            character(64) :: subroutine_name = 'geo_ell'
        !%------------------------------------------------------------------
        !% Preliminaries:
            Npack => npack_elemP(thisColP)
            if (Npack < 1) return
            if (setting%Debug%File%geometry) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases:
            thisP               => elemP(1:Npack,thisColP)
            ell                 => elemR(:,er_ell)
            head                => elemR(:,er_Head)
            hydDepth            => elemR(:,er_HydDepth)
            area                => elemR(:,er_Area)
            topwidth            => elemR(:,er_Topwidth)
            ZbreadthMax         => elemR(:,er_ZbreadthMax)
            breadthMax          => elemR(:,er_BreadthMax)
            areaBelowBreadthMax => elemR(:,er_AreaBelowBreadthMax)
        !%-------------------------------------------------------------------

        where (head(thisP) .le. ZbreadthMax(thisP))
            ell(thisP) =  hydDepth(thisP)
        elsewhere
            ell(thisP) = ( (head(thisP) - ZbreadthMax(thisP)) * breadthMax(thisP) &
                            + areaBelowBreadthMax(thisP) ) / breadthMax(thisP)
        endwhere

        !%-------------------------------------------------------------------
            if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_ell_from_head
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
            elemR(thisP,er_HydDepth)      = elemR(thisP,er_Depth)
            elemR(thisP,er_ell)           = elemR(thisP,er_Depth)
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
        (idx, indepth) result (outvalue)
        !%------------------------------------------------------------------
        !% Descriptions:
        !% computes the area for a given depth of a single element
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(in)   :: indepth
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
                (idx,tt_area, indepth, setting%ZeroValue%Depth)

        case (parabolic)
            !outvalue = parabolic_area_from_depth_singular (idx, indepth)
            outA = llgeo_parabolic_area_from_depth_pure(iA, depthA)
            outvalue = outA(1)

        case (power_function)
            print *, 'CODE ERROR powerfunction geometry not complete'
            call util_crashpoint(55098723)   

        case (rectangular)
            outA = llgeo_rectangular_area_from_depth_pure (iA, depthA)
            outvalue = outA(1)

        case (trapezoidal)
            outA = llgeo_trapezoidal_area_from_depth_pure (iA, depthA)
            outvalue = outA(1)

        case (triangular)
            outA= llgeo_triangular_area_from_depth_pure (iA, depthA)
            outvalue = outA(1)

        !% --- closed conduits   
        case (arch, basket_handle, catenary, circular, eggshaped, gothic, &
            horiz_ellipse, horseshoe, semi_circular, semi_elliptical,     &
            vert_ellipse)

           ! outvalue = arch_area_from_depth_singular (idx, indepth)
            outvalue = llgeo_tabular_from_depth_singular &
                    (idx, indepth, fullArea(idx), setting%ZeroValue%Depth, Atable)
            
        case (filled_circular)
            outvalue = llgeo_filled_circular_area_from_depth_singular (idx, indepth)     
            
        case (mod_basket)
            outvalue = llgeo_mod_basket_area_from_depth_singular (idx, indepth)

        case (rectangular_closed)
            outvalue = llgeo_rectangular_closed_area_from_depth_singular (idx, indepth)

        case (rect_round )
            outvalue = llgeo_rect_round_area_from_depth_singular (idx, indepth)

        case (rect_triang)
            outvalue = llgeo_rectangular_triangular_area_from_depth_singular (idx, indepth)

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
            call util_crashpoint(33234)

        case (force_main)
            print *, 'CODE ERROR: area for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'in ',trim(subroutine_name)   
            print *, 'This should never be reached as a force_main is not a valid geometryType'
            call util_crashpoint(33234)

        case default
            print *, 'CODE ERROR: area for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(33234)

        end select
           
    end function geo_area_from_depth_singular
!%
!%==========================================================================    
!%==========================================================================
!%
    real(8) function geo_topwidth_from_depth_singular &
        (idx, indepth) result (outvalue)
        !%------------------------------------------------------------------
        !% Descriptions:
        !% computes the topwidth for a given depth of a single element
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(in)  :: indepth
            integer, intent(in)  :: idx
            integer, dimension(1):: iA
            real(8), dimension(1):: depthA, outA
            real(8), pointer     :: depth(:), fullarea(:), fulltopwidth(:)
            real(8), pointer     :: TTable(:)
            character(64) :: subroutine_name = 'geo_topwidth_from_depth_singular'
        !%------------------------------------------------------------------
            depth        => elemR(:,er_Depth)
            fullarea     => elemR(:,er_FullArea)
            fulltopwidth => elemR(:,er_FullTopWidth)
        !%------------------------------------------------------------------
        !% size(1) arrays
            iA(1)     = idx
            depthA(1) = indepth
        !%------------------------------------------------------------------
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

        case (power_function)
            print *, 'CODE ERROR: topwidth for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)

        case (rectangular)
            outA = llgeo_rectangular_topwidth_from_depth_pure  (iA, depthA)
            outvalue = outA(1)

        case (trapezoidal)
            outA = llgeo_trapezoidal_topwidth_from_depth_pure (iA, depthA)
            outvalue = outA(1)

        case (triangular)
            outA = llgeo_triangular_topwidth_from_depth_pure (iA, depthA)
            outvalue = outA(1)
        
        case (irregular)
            outvalue = irregular_geometry_from_depth_singular &
                (idx,tt_width, indepth, setting%ZeroValue%TopWidth)

        !% -----CLOSED CONDUITS ---------------------------------------------
        case (arch, basket_handle, catenary, circular, eggshaped, gothic,  &
                horiz_ellipse, horseshoe, semi_circular, semi_elliptical,  &
                vert_ellipse)

            outvalue = llgeo_tabular_from_depth_singular &
                (idx, depth(idx), fullTopWidth(idx), setting%ZeroValue%Depth, Ttable)

        case (filled_circular)
            outvalue = llgeo_filled_circular_topwidth_from_depth_singular  (idx, indepth)

        case (mod_basket)
            outvalue = llgeo_mod_basket_topwidth_from_depth_singular (idx, indepth)

        case (rectangular_closed)
            outvalue = llgeo_rectangular_closed_topwidth_from_depth_singular  (idx, indepth)

        case (rect_round)
            outvalue = llgeo_rect_round_topwidth_from_depth_singular (idx, indepth)

        case (rect_triang)
            outvalue = llgeo_rectangular_triangular_topwidth_from_depth_singular  (idx, indepth)
        
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
        (idx, indepth) result (outvalue)
        !%------------------------------------------------------------------
        !% Descriptions:
        !% computes the perimeter for a given depth of a single element
        !% Note that ALL values computed herein are based on the input depth
        !% which may NOT be the depth of the element idx. Thus, where we
        !% need more values than the depth to compute the perimeter we must 
        !% compute the other values as temporary values!
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(in)   :: indepth
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

        case (trapezoidal)
            outA = llgeo_trapezoidal_perimeter_from_depth_pure (iA, indepthA)
            outvalue = outA(1)

        case (triangular)
            outA = llgeo_triangular_perimeter_from_depth_pure (iA, indepthA)
            outvalue = outA(1)
        
        case (irregular)
            outvalue = irregular_geometry_from_depth_singular &
                (idx,tt_area, indepth, setting%ZeroValue%Depth)

        !% --- CLOSED CONDUITS   
        case (arch, basket_handle, circular, eggshaped, horiz_ellipse, horseshoe, vert_ellipse)
            tempArea(1)  = llgeo_tabular_from_depth_singular &
                    (idx, indepth, fullArea(idx), setting%ZeroValue%Depth, Atable)

            tempHydRadius(1) = llgeo_tabular_from_depth_singular &
                    (idx, indepth, fullHydRadius(idx), setting%ZeroValue%Depth, Rtable)

            outA = llgeo_perimeter_from_hydradius_and_area_pure &
                    (iA, tempHydRadius, tempArea)
            outvalue = outA(1)
        
        case (catenary, gothic, semi_circular, semi_elliptical)  

            tempArea(1)  = llgeo_tabular_from_depth_singular &
                    (idx, indepth, fullArea(idx), setting%ZeroValue%Depth, Atable)

            temphydRadius(1)= llgeo_tabular_hydradius_from_area_and_sectionfactor_singular &
                (idx, tempArea(1), fullhydradius(idx), setting%ZeroValue%Area, Stable)

            outA = llgeo_perimeter_from_hydradius_and_area_pure &
                        (iA, tempHydradius, tempArea)

            outvalue = outA(1)

        case (filled_circular)
            outvalue = llgeo_filled_circular_perimeter_from_depth_singular (idx, indepth)

        case (mod_basket)
            outvalue = llgeo_mod_basket_perimeter_from_depth_singular (idx, indepth)

        case (rectangular_closed)
            outvalue = llgeo_rectangular_closed_perimeter_from_depth_singular (idx, indepth) 

        case (rect_round )
            outvalue = llgeo_rect_round_perimeter_from_depth_singular (idx, indepth)

        case (rect_triang)
            outvalue =llgeo_rectangular_triangular_perimeter_from_depth_singular (idx, indepth)

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
            call util_crashpoint(33234)
        end select

    end function geo_perimeter_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function geo_hyddepth_from_area_and_topwidth_singular &
        (idx, area, topwidth)  result (outvalue)
        !%------------------------------------------------------------------
        !% Descriptions:
        !% computes the hydraulic depth for area and topwidth of a single
        !% element
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(in)  :: area, topwidth
            integer, intent(in)  :: idx
            !real(8)              :: temp1, temp2
            character(64) :: subroutine_name = 'geo_hyddepth_from_area_and_topwidth_singular'
        !%------------------------------------------------------------------   
     
        outvalue = area / topwidth
        
    end function geo_hyddepth_from_area_and_topwidth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function geo_hydradius_from_area_and_perimeter_singular &
        (idx, area, perimeter)  result (outvalue)
        !%------------------------------------------------------------------
        !% Descriptions:
        !% computes the hydraulic radius for area and perimeter of a single
        !% element
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(in)  :: area, perimeter
            integer, intent(in)  :: idx
            !real(8)              :: temp1, temp2
            character(64) :: subroutine_name = 'geo_hydradius_from_area_and_perimeter_singular'
        !%------------------------------------------------------------------   
    
        outvalue = area / perimeter
    
    end function geo_hydradius_from_area_and_perimeter_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function geo_perimeter_from_area_and_hydradius_singular &
        (idx, area, hydradius) result (outvalue)    
        !%------------------------------------------------------------------
        !% Descriptions:
        !% computes the perimeter for area and hydraulic radius of a single
        !% element
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(in)  :: area, hydradius
            integer, intent(in)  :: idx
            character(64) :: subroutine_name = 'geo_perimeter_from_area_and_hydradius_singular'
        !%------------------------------------------------------------------   
    
        outvalue = area / hydradius
    
    end function geo_perimeter_from_area_and_hydradius_singular
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
!% END OF MODULE
!%=========================================================================
end module geometry