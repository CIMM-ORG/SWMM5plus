module geometry

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use rectangular_channel
    use trapezoidal_channel
    use circular_conduit
    use storage_geometry
    use adjust
    use utility_profiler


    implicit none

!%-----------------------------------------------------------------------------
!% Description:
!% Geometry computations
!%

    private

    public :: geometry_toplevel
    public :: geo_assign_JB  !BRHbugfix 20210813

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine geometry_toplevel (whichTM)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Input whichTM is one of ETM, AC, or ALLtm
        !% This should never be called for diagnostic arrays
        !% Note that the elemPGx arrays contain only time-marched elements so they
        !% will only handle CC and JM elements as the JB elements are not time-marched.
        !%-----------------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: whichTM
            integer, pointer :: elemPGx(:,:), npack_elemPGx(:), col_elemPGx(:)
            integer, pointer :: thisColP_surcharged, thisColP_NonSurcharged, thisColP_all
            integer, pointer :: thisColP_JM, thisColP_JB, thisColP_ClosedElems
            logical :: isreset
            integer, allocatable :: tempP(:) !% debugging
            character(64) :: subroutine_name = 'geometry_toplevel'
        !%-----------------------------------------------------------------------------
        !% Preliminaries
            if (icrash) return
            if (setting%Debug%File%geometry) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        !% Aliases
        !% set the packed geometry element array (elemPG) to use and columns of the
        !% packed elemP to use
            select case (whichTM)
                case (ALLtm)
                    elemPGx                => elemPGalltm(:,:)
                    npack_elemPGx          => npack_elemPGalltm(:)
                    col_elemPGx            => col_elemPGalltm(:)
                    thisColP_JM            => col_elemP(ep_JM_ALLtm)
                    thisColP_JB            => col_elemP(ep_JB_ALLtm)
                    thisColP_surcharged    => col_elemP(ep_Surcharged_ALLtm)
                    thisColP_NonSurcharged => col_elemP(ep_NonSurcharged_ALLtm)
                    thisColP_all           => col_elemP(ep_ALLtm)
                    thisColP_ClosedElems   => col_elemP(ep_Closed_Elements)
                case (ETM)
                    elemPGx                => elemPGetm(:,:)
                    npack_elemPGx          => npack_elemPGetm(:)
                    col_elemPGx            => col_elemPGetm(:)
                    thisColP_JM            => col_elemP(ep_JM_ETM)
                    thisColP_JB            => col_elemP(ep_JB_ETM)
                    thisColP_surcharged    => col_elemP(ep_Surcharged_ETM)
                    thisColP_NonSurcharged => col_elemP(ep_NonSurcharged_ETM)
                    thisColP_all           => col_elemP(ep_ETM)
                    thisColP_ClosedElems   => col_elemP(ep_Closed_Elements)
                case (AC)
                    elemPGx                => elemPGac(:,:)
                    npack_elemPGx          => npack_elemPGac(:)
                    col_elemPGx            => col_elemPGac(:)
                    thisColP_JM            => col_elemP(ep_JM_AC)
                    thisColP_JB            => col_elemP(ep_JB_AC)
                    thisColP_surcharged    => col_elemP(ep_Surcharged_AC)
                    thisColP_NonSurcharged => col_elemP(ep_NonSurcharged_AC)
                    thisColP_all           => col_elemP(ep_AC)
                    thisColP_ClosedElems   => col_elemP(ep_Closed_Elements)
                case default
                    print *, 'CODE ERROR: time march type unknown for # ', whichTM
                    print *, 'which has key ',trim(reverseKey(whichTM))
                    stop 7389
            end select
        !%-----------------------------------------------------------------------------
        !% STATUS: at this point we know volume on Non-surcharged CC, JM,
        !% elements and head on all surcharged CC, JM elements

        !% assign all geometry for surcharged elements CC, JM (and JB?)
        call geo_surcharged (thisColP_surcharged)

        !% reset all zero or near-zero volumes in non-surcharged CC and JM
        call adjust_limit_by_zerovalues (er_Volume, setting%ZeroValue%Volume, thisColP_NonSurcharged)

        !% compute the depth on all non-surcharged elements of CC and JM
        call geo_depth_from_volume (elemPGx, npack_elemPGx, col_elemPGx)

        !% reset all zero or near-zero depths in non-surcharged CC and JM
        call adjust_limit_by_zerovalues (er_Depth, setting%ZeroValue%Depth, thisColP_NonSurcharged)

        !% compute the head on all non-surcharged elements of CC and JM
        call geo_head_from_depth (thisColP_NonSurcharged)

        !% limit volume for incipient surcharge. This is done after depth is computed
        !% so that the "depth" algorithm can include depths greater than fulldepth
        !% as a way to handle head for incipient surcharge.
        !call geo_limit_incipient_surcharge (er_Volume, er_FullVolume, thisColP_NonSurcharged)
        call geo_limit_incipient_surcharge (er_Volume, er_FullVolume, thisColP_NonSurcharged,.true.) !% 20220124brh

        !% limit depth for incipient surcharged. This is done after head is computed
        !% so that the depth algorithm can include depths greater than fulldepth to
        !% handle incipient surcharge
        !call geo_limit_incipient_surcharge (er_Depth, er_FullDepth, thisColP_NonSurcharged)
        call geo_limit_incipient_surcharge (er_Depth, er_FullDepth, thisColP_NonSurcharged,.false.) !% 20220124brh

        !% STATUS: at this point we know depths and heads in all CC, JM elements
        !% (surcharged and nonsurcharged) with limiters for conduit depth and zero depth
           
        !% assign the head, depth, geometry on junction branches JB based on JM head
        call geo_assign_JB (whichTM, thisColP_JM)

        !% STATUS at this point we know geometry on all JB and all surcharged, with
        !% depth, head, volume on all non-surcharged or incipient surcharge.

        !% compute area from volume for CC, JM nonsurcharged
        call geo_area_from_volume (thisColP_NonSurcharged)

        !% reset all zero or near-zero areas in non-surcharged CC and JM
        call adjust_limit_by_zerovalues (er_Area, setting%ZeroValue%Area, thisColP_NonSurcharged)

        !% compute topwidth from depth for all CC, JM nonsurcharged
        call geo_topwidth_from_depth (elemPGx, npack_elemPGx, col_elemPGx)

        !% reset all zero or near-zero topwidth in non-surcharged CC and JM
        !% but do not change the eYN(:,eYN_isZeroDepth) mask
        call adjust_limit_by_zerovalues (er_Topwidth, setting%ZeroValue%Topwidth, thisColP_NonSurcharged)

        !% compute perimeter from maximum depth for all CC, JM nonsurcharged
        call geo_perimeter_from_depth (elemPGx, npack_elemPGx, col_elemPGx)

        !% compute hyddepth
        call geo_hyddepth (elemPGx, npack_elemPGx, col_elemPGx)

        !% compute hydradius
        call geo_hydradius_from_area_perimeter (thisColP_NonSurcharged)

        !% make adjustments for slots on closed elements only for ETM
        if (whichTM .eq. ETM) then
            call geo_slot_adjustments (thisColP_ClosedElems)
        end if

        !% the modified hydraulic depth "ell" is used for AC computations and
        !% for Froude number computations on all elements, whether ETM or AC.
        call geo_ell (thisColP_all)

        !% Set JM values that are not otherwise defined
        call geo_JM_values ()

        !% compute the dHdA that are only for AC nonsurcharged
        if (whichTM .ne. ETM) then
            call geo_dHdA (ep_NonSurcharged_AC)
        end if

        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geometry_toplevel

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
        !% The following initializations are unknown as of 20211215
        !% Preissmann_Celerity, SlotVolume, SlotArea, SlotWidth, SlotDepth
        !% SmallVolume_xxx,
        !%-------------------------------------------------------------------
            integer, intent(in) :: whichTM, thisColP_JM
            integer, pointer ::  Npack, thisP(:), BranchExists(:), thisSolve(:),  tM
            real(8), pointer :: area(:), depth(:), head(:), hyddepth(:), hydradius(:)
            real(8), pointer :: length(:), perimeter(:), topwidth(:), velocity(:)
            real(8), pointer :: volume(:), zBtm(:), Kfac(:), dHdA(:), ell(:)
            real(8), pointer :: zCrown(:), fullarea(:), fulldepth(:), fullperimeter(:)
            real(8), pointer :: fullhyddepth(:)
            real(8), pointer :: grav
            integer :: tB, ii, kk
        !% branchsign assume branches are ordered as nominal inflow, outflow, inflow...
        !real(8) :: branchsign(6) = [+oneR,-oneR,+oneR,-oneR,+oneR,-oneR]

        !% thisColP_JM is the column for the junction mains of a particular
        !% whichTM. For ALL ep_JM, for ETM, ep_JM_ETM, for AC ep_JM_AC
            integer, allocatable :: tempP(:)
            character(64) :: subroutine_name = 'geo_assign_JB'
        !%---------------------------------------------------------------------
        !% Preliminaries
            if (icrash) return
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
            head          => elemR(:,er_Head)
            hyddepth      => elemR(:,er_HydDepth)
            hydradius     => elemR(:,er_HydRadius)
            length        => elemR(:,er_Length)
            perimeter     => elemR(:,er_Perimeter)
            topwidth      => elemR(:,er_Topwidth)
            velocity      => elemR(:,er_Velocity)
            volume        => elemR(:,er_Volume)
            zBtm          => elemR(:,er_Zbottom)
            zCrown        => elemR(:,er_Zcrown)
            fullarea      => elemR(:,er_FullArea)
            fulldepth     => elemR(:,er_FullDepth)
            fullhyddepth  => elemR(:,er_FullHydDepth)
            fullperimeter => elemR(:,er_FullPerimeter)
            Kfac          => elemSR(:,esr_JunctionBranch_Kfactor)
            BranchExists  => elemSI(:,esi_JunctionBranch_Exists)
            thisSolve     => elemI(:,ei_tmType)
            grav => setting%Constant%gravity
        !%------------------------------------------------------------------

        if (Npack > 0) then
            thisP  => elemP(1:Npack,thisColP_JM)
            !% cycle through the all the main junctions and each of its branches
            do ii=1,Npack
                tM => thisP(ii) !% junction main ID
                !% only execute for whichTM of ALL or thisSolve (of JM) matching input whichTM
                if ((whichTM == ALLtm) .or. (thisSolve(tM) == whichTM)) then
                    !% cycle through the possible junction branches
                    do kk=1,max_branch_per_node
                        tB = tM + kk !% junction branch ID
                        if (BranchExists(tB) == 1) then
                            !% only when a branch exists.
                            if ( head(tM) > zBtm(tB) ) then
                                !% for main head above branch bottom entrance use a head
                                !% loss approach. The branchsign and velocity sign ensure
                                !% the headloss is added to an inflow and subtracted at
                                !% an outflow
                                !% Note this is a time-lagged velocity as the JB velocity
                                !% is not updated until after face interpolation                                
                                head(tB) = head(tM)  + branchsign(kk) * sign(oneR,velocity(tB)) &
                                    * (Kfac(tB) / (twoR * grav)) * (velocity(tB)**twoR)
                            else
                                !% for main head below the branch bottom entrance we assign a
                                !% Froude number of one on an inflow to the junction main. Note
                                !% an outflow from a junction main for this case gets head
                                !% of z_bottom of the branch (zero depth).
                                !% Note this is a time-lagged velocity as the JB velocity
                                !% is not updated until after face interpolation
                                head(tB) = zBtm(tB)  &
                                    + onehalfR * (oneR + branchsign(kk) * sign(oneR,velocity(tB))) &
                                    *(velocity(tB)**twoR) / grav
                            end if
                           
                            
                            !% compute provisional depth
                            depth(tB) = head(tB) - zBtm(tB)
                            
                            if (depth(tB) .ge. fulldepth(tB)) then
                                !% surcharged or incipient surcharged
                                depth(tB)     = fulldepth(tB)
                                area(tB)      = fullarea(tB)
                                hyddepth(tB)  = fullhyddepth(tB)
                                perimeter(tB) = fullperimeter(tB)
                                topwidth(tB)  = setting%ZeroValue%Topwidth
                                hydRadius(tB) = fulldepth(tB) / fullperimeter(tB)
                                dHdA(tB)      = oneR / setting%ZeroValue%Topwidth
                            elseif ((depth(tB) < setting%ZeroValue%Depth) .and. (setting%ZeroValue%UseZeroValues)) then
                                !% negligible depth is treated with ZeroValues
                                depth(tB)     = setting%ZeroValue%Depth
                                area(tB)      = setting%ZeroValue%Area
                               
                                ! HACK fix
                                !f (elemI(tB,ei_geometryType) == rectangular)  then
                                !    topwidth(tB) = elemSGR(tB,esgr_Rectangular_Breadth)
                                !else
                                    topwidth(tB)  = setting%ZeroValue%Topwidth
                                !end if
                                !% HACK
                                hyddepth(tB)  = setting%ZeroValue%Area / topwidth(tB)
                                ! hyddepth(tB)  = setting%ZeroValue%Area / setting%ZeroValue%Topwidth
                                !% HACK
                                perimeter(tB) = topwidth(tB) + setting%ZeroValue%Depth
                                ! perimeter(tB) = setting%ZeroValue%Topwidth + setting%ZeroValue%Depth
                                hydRadius(tB) = setting%ZeroValue%Area / perimeter(tB)
                                !% HACK
                                dHdA(tB)      = oneR / topwidth(tB)
                                ! dHdA(tB)      = oneR / setting%ZeroValue%Topwidth
                            elseif ((depth(tB) .le. zeroR) .and. (.not. setting%ZeroValue%UseZeroValues)) then
                                !% negative depth is treated as exactly zero
                                depth(tB) = zeroR
                                area(tB)  = zeroR
                                topwidth(tB) = zeroR
                                hydDepth(tB) = zeroR
                                perimeter(tB) = zeroR
                                hydRadius(tB) = zeroR
                                dHdA(tB)      = oneR / setting%ZeroValue%Topwidth
                                !dHdA(tB)      = setting%ZeroValue%Topwidth
                            else
                                !% not surcharged and non-negligible depth
                                select case (elemI(tB,ei_geometryType))
                                case (rectangular)
                                    area(tB)     = rectangular_area_from_depth_singular (tB)
                                    topwidth(tB) = rectangular_topwidth_from_depth_singular (tB)
                                    hydDepth(tB) = rectangular_hyddepth_from_depth_singular (tB)
                                    perimeter(tB)= rectangular_perimeter_from_depth_singular (tB)
                                    hydRadius(tB)= rectangular_hydradius_from_depth_singular (tB)
                                    ell(tB)      = hydDepth(tB) !geo_ell_singular (tB) !BRHbugfix 20210812 simpler for rectangle
                                    dHdA(tB)     = oneR / topwidth(tB)
                                case (trapezoidal)
                                    area(tB)     = trapezoidal_area_from_depth_singular (tB)
                                    topwidth(tB) = trapezoidal_topwidth_from_depth_singular (tB)
                                    hydDepth(tB) = trapezoidal_hyddepth_from_depth_singular (tB)
                                    perimeter(tB)= trapezoidal_perimeter_from_depth_singular (tB)
                                    hydRadius(tB)= trapezoidal_hydradius_from_depth_singular (tB)
                                    ell(tB)      = hydDepth(tB) !geo_ell_singular (tB) !BRHbugfix 20210812 simpler for trapezoid
                                    dHdA(tB)     = oneR / topwidth(tB)
                                case (circular)
                                    area(tB)     = circular_area_from_depth_singular (tB)
                                    topwidth(tB) = circular_topwidth_from_depth_singular (tB)
                                    hydDepth(tB) = circular_hyddepth_from_topwidth_singular (tB)
                                    hydRadius(tB)= circular_hydradius_from_depth_singular (tB)
                                    perimeter(tB)= circular_perimeter_from_hydradius_singular (tB)
                                    ell(tB)      = hydDepth(tB) !geo_ell_singular (tB) !BRHbugfix 20210812 simpler for circular
                                    dHdA(tB)     = oneR / topwidth(tB)
                                case default
                                    print *, 'CODE ERROR: geometry type unknown for # ', elemI(tB,ei_geometryType)
                                    print *, 'which has key ',trim(reverseKey(elemI(tB,ei_geometryType)))
                                    print *, 'in ',trim(subroutine_name)
                                    stop 399848
                                end select
                            end if
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
!% PRIVATE
!%==========================================================================
!%
    subroutine geo_surcharged (thisColP)
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

        character(64) :: subroutine_name = 'geo_surcharged'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        Npack => npack_elemP(thisColP)
        !%-------------------------------------------------
        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        if (Npack > 0) then
            thisP => elemP(1:Npack,thisColP)
            elemR(thisP,er_Volume)    = elemR(thisP,er_FullVolume)
            elemR(thisP,er_Area)      = elemR(thisP,er_FullArea)
            elemR(thisP,er_Depth)     = elemR(thisP,er_FullDepth)
            elemR(thisP,er_Perimeter) = elemR(thisP,er_FullPerimeter)
            elemR(thisP,er_HydDepth)  = elemR(thisP,er_FullHydDepth)
            elemR(thisP,er_HydRadius) = elemR(thisP,er_FullArea) / elemR(thisP,er_FullPerimeter)
            elemR(thisP,er_Topwidth)  = setting%ZeroValue%Topwidth
        end if

        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        end subroutine geo_surcharged
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_depth_from_volume (elemPGx, npack_elemPGx, col_elemPGx)
        !%------------------------------------------------------------------
        !% Description:
        !% This solves nonsurcharged CCJM elements because of PGx arrays
        !% The elemPGx determines whether this is ALLtm, ETM or AC elements
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: elemPGx(:,:), npack_elemPGx(:), col_elemPGx(:)
            integer, pointer :: Npack, thisCol
            character(64) :: subroutine_name = 'geo_depth_from_volume'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (icrash) return
            if (setting%Debug%File%geometry) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------    
        !% cycle through different geometries
        !% RECTANGULAR
        thisCol => col_elemPGx(epg_CC_rectangular_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call rectangular_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% TRAPEZOIDAL
        thisCol => col_elemPGx(epg_CC_trapezoidal_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call trapezoidal_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% CIRCULAR
        thisCol => col_elemPGx(epg_CC_circular_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call circular_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% HACK Needs additional geometries

        !% JM with functional geomtery
        thisCol => col_elemPGx(epg_JM_functional_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call storage_functional_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% JM with tabular geomtery
        thisCol => col_elemPGx(epg_JM_tabular_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call storage_tabular_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% JM with artificial storage
        thisCol => col_elemPGx(epg_JM_artificial_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call storage_artificial_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !%-------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%geometry) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_depth_from_volume
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_limit_incipient_surcharge (geocol, fullcol, thisColP, isVolume)  !% 20220124brh
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Sets volume/depth limit to full volume for incipient surcharge.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisColP, geocol, fullcol
        logical, intent(in) :: isVolume !% 20220124brh
        integer, pointer :: Npack, thisP(:)
        real(8), pointer :: geovalue(:), fullvalue(:)
        real(8), pointer :: overflow(:)  !% 20220124brh

        character(64) :: subroutine_name = 'geo_limit_incipient_surcharge'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        Npack      => npack_elemP(thisColP)
        geovalue   => elemR(:,geocol)
        fullvalue  => elemR(:,fullcol)
        overflow   => elemR(:,er_VolumeOverFlow)  !% 20220124brh
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        if (Npack > 0) then
            thisP      => elemP(1:Npack,thisColP)
             !% 20220124brh REWRITE START
            if (isVolume) then
                where (geovalue(thisP) > fullvalue(thisP))
                    overflow(thisP) = geovalue(thisP) - fullvalue(thisP) + overflow(thisP)  !% 20220124brh
                    geovalue(thisP) = fullvalue(thisP)
                endwhere
            else
                where (geovalue(thisP) > fullvalue(thisP))
                    geovalue(thisP) = fullvalue(thisP)
                endwhere
            end if
            !% 20220124brh REWRITE END
            !where (geovalue(thisP) > fullvalue(thisP))
            !    geovalue(thisP) = fullvalue(thisP)
            !endwhere
        end if

        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_limit_incipient_surcharge
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_head_from_depth (thisColP)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes head from depth for non-surcharged elements of CC, JM
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisColP
        integer, pointer :: Npack, thisP(:)
        real(8), pointer :: depth(:), fulldepth(:), head(:), Zbtm(:)

        character(64) :: subroutine_name = 'geo_head_from_depth'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        Npack     => npack_elemP(thisColP)
        depth     => elemR(:,er_Depth)
        fulldepth => elemR(:,er_FullDepth)
        head      => elemR(:,er_Head)
        Zbtm      => elemR(:,er_Zbottom)
        !%-----------------------------------------------------------------------------
        !%
        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        if (Npack > 0) then
            thisP     => elemP(1:Npack,thisColP)
            head(thisP) = depth(thisP) + Zbtm(thisP)
        end if

        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_head_from_depth
!%
!%==========================================================================
!
!%==========================================================================
!%
    subroutine geo_area_from_volume (thisColP)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% sets area = volume/length which is common to all nonsurcharged elements
        !% Note this assumes volume has been limited by surcharge and zero values
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisColP
        integer, pointer :: thisP(:), Npack
        real(8), pointer :: area(:), volume(:), length(:)

        character(64) :: subroutine_name = 'geo_area_from_volume'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        Npack  => npack_elemP(thisColP)
        area   => elemR(:,er_Area)
        volume => elemR(:,er_Volume)
        length => elemR(:,er_Length)
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        if (Npack > 0) then
            thisP  => elemP(1:Npack,thisColP)
            area(thisP) = volume(thisP) / length(thisP)
        end if

        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_area_from_volume
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_topwidth_from_depth &
        (elemPGx, npack_elemPGx, col_elemPGx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth given depth of a non-surcharged element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: elemPGx(:,:)
        integer, target, intent(in) :: npack_elemPGx(:), col_elemPGx(:)
        integer, pointer :: Npack, thisCol

        character(64) :: subroutine_name = 'geo_topwidth_from_depth'
        !%-----------------------------------------------------------------------------
        !% cycle through different geometries
        if (icrash) return
        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        Npack => npack_elemPGx(epg_CC_rectangular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_nonsurcharged)
            call rectangular_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        Npack => npack_elemPGx(epg_CC_trapezoidal_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_trapezoidal_nonsurcharged)
            call trapezoidal_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        Npack => npack_elemPGx(epg_CC_circular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_circular_nonsurcharged)
            call circular_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% HACK NEED OTHER GEOMETRIES
        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_topwidth_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_perimeter_from_depth &
        (elemPGx, npack_elemPGx, col_elemPGx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the wetted perimeter given depth of a non-surcharged element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: elemPGx(:,:)
        integer, target, intent(in) :: npack_elemPGx(:), col_elemPGx(:)
        integer, pointer :: Npack, thisCol

        character(64) :: subroutine_name = 'geo_perimeter_from_depth'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% cycle through different geometries
        Npack => npack_elemPGx(epg_CC_rectangular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_nonsurcharged)
            call rectangular_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% cycle through different geometries
        Npack => npack_elemPGx(epg_CC_trapezoidal_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_trapezoidal_nonsurcharged)
            call trapezoidal_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        Npack => npack_elemPGx(epg_CC_circular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_circular_nonsurcharged)
            call circular_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% HACK NEED OTHER GEOMETRIES
        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_perimeter_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_hyddepth (elemPGx, npack_elemPGx, col_elemPGx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Note that hyddepth is the average depth, which is only area/topwidth
        !% for a simple open channel, and does not apply above midpoint in a
        !% conduit
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: elemPGx(:,:)
        integer, target, intent(in) :: npack_elemPGx(:), col_elemPGx(:)
        integer, pointer :: Npack, thisCol

        character(64) :: subroutine_name = 'geo_hyddepth_from_depth'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% cycle through different geometries
        Npack => npack_elemPGx(epg_CC_rectangular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_nonsurcharged)
            call rectangular_hyddepth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% cycle through different geometries
        Npack => npack_elemPGx(epg_CC_trapezoidal_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_trapezoidal_nonsurcharged)
            call trapezoidal_hyddepth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% cycle through different geometries
        Npack => npack_elemPGx(epg_CC_circular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_circular_nonsurcharged)
            call circular_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        end if

        !% HACK need other geometries
        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_hyddepth
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_hydradius_from_area_perimeter (thisColP)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% sets area = volume/length which is common to all nonsurcharged elements
        !% Note this assumes volume has been limited by surcharge and zero values
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisColP
        integer, pointer :: thisP(:), Npack
        real(8), pointer :: area(:), hydradius(:), perimeter(:)

        character(64) :: subroutine_name = 'geo_hydradius_from_area_perimeter'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        Npack     => npack_elemP(thisColP)
        area      => elemR(:,er_Area)
        hydradius => elemR(:,er_HydRadius)
        perimeter => elemR(:,er_Perimeter)
        !%-----------------------------------------------------------------------------

        if (Npack > 0) then
            thisP     => elemP(1:Npack,thisColP)
            hydradius(thisP) = area(thisP) / perimeter(thisP)
        end if

        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_hydradius_from_area_perimeter
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_ell (thisColP)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes the value of "ell" -- the modified hydraulic depth
        !% used as a length scale in AC method
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisColP
        integer, pointer :: thisP(:), Npack
        real(8), pointer :: ell(:), head(:), area(:), topwidth(:), hydDepth(:)
        real(8), pointer :: ZbreadthMax(:), breadthMax(:), areaBelowBreadthMax(:)

        character(64) :: subroutine_name = 'geo_ell'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        Npack               => npack_elemP(thisColP)
        ell                 => elemR(:,er_ell)
        head                => elemR(:,er_Head)
        hydDepth            => elemR(:,er_HydDepth)
        area                => elemR(:,er_Area)
        topwidth            => elemR(:,er_Topwidth)
        ZbreadthMax         => elemR(:,er_ZbreadthMax)
        breadthMax          => elemR(:,er_BreadthMax)
        areaBelowBreadthMax => elemR(:,er_AreaBelowBreadthMax)
        !%-----------------------------------------------------------------------------

        if (Npack > 0) then
            thisP               => elemP(1:Npack,thisColP)
            where (head(thisP) .le. ZbreadthMax(thisP))
                ell(thisP) =  hydDepth(thisP)
            elsewhere
                ell(thisP) = ( (head(thisP) - ZbreadthMax(thisP)) * breadthMax(thisP) &
                                + areaBelowBreadthMax(thisP) ) / breadthMax(thisP)
            endwhere
        end if

        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_ell
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function geo_ell_singular (indx) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes the value of "ell" -- the length scale for the AC method for
        !% a single index point
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), pointer :: head(:), area(:), topwidth(:)
        real(8), pointer :: ZbreadthMax(:), breadthMax(:), areaBelowBreadthMax(:)

        character(64) :: subroutine_name = 'geo_ell_singular'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        head                => elemR(:,er_Head)
        area                => elemR(:,er_Area)
        topwidth            => elemR(:,er_Topwidth)
        ZbreadthMax         => elemR(:,er_ZbreadthMax)
        breadthMax          => elemR(:,er_BreadthMax)
        areaBelowBreadthMax => elemR(:,er_AreaBelowBreadthMax)
        !%-----------------------------------------------------------------------------
        if (head(indx) .le. ZbreadthMax(indx)) then
            outvalue =  area(indx) / topwidth(indx)
        else
            outvalue = ( (head(indx) - ZbreadthMax(indx)) * breadthMax(indx) &
                            + areaBelowBreadthMax(indx) ) / breadthMax(indx)
        end if

        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end function geo_ell_singular
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
        if (icrash) return
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
    subroutine geo_slot_adjustments (thisColP)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% This subroutine adds back the slot geometry in all the closed elements
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisColP
        integer, pointer    :: thisP(:), Npack
        real(8), pointer    :: SlotWidth(:), SlotVolume(:), SlotDepth(:), SlotArea(:)
        real(8), pointer    :: volume(:), depth(:), area(:), head(:), SlotHydRadius(:)
        real(8), pointer    :: hydRadius(:), ell(:), zbottom(:), fullDepth(:)

        character(64) :: subroutine_name = 'geo_slot_adjustments'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        Npack      => npack_elemP(thisColP)
        volume     => elemR(:,er_Volume)
        depth      => elemR(:,er_Depth)
        fullDepth  => elemR(:,er_FullDepth)
        area       => elemR(:,er_Area)
        head       => elemR(:,er_Head)
        ell        => elemR(:,er_ell)
        hydRadius  => elemR(:,er_HydRadius)
        zbottom    => elemR(:,er_Zbottom)
        SlotWidth  => elemR(:,er_SlotWidth)
        SlotVolume => elemR(:,er_SlotVolume)
        SlotDepth  => elemR(:,er_SlotDepth)
        SlotArea   => elemR(:,er_SlotArea)
        SlotHydRadius => elemR(:,er_SlotHydRadius)
        !%-----------------------------------------------------------------------------

        if (Npack > 0) then
            thisP    => elemP(1:Npack,thisColP)
            where (SlotVolume(thisP) .gt. zeroR) 
                volume(thisP) = volume(thisP) + SlotVolume(thisP)
                area(thisP)   = area(thisP)   + SlotArea(thisP)
                depth(thisP)  = depth(thisP)  + SlotDepth(thisP)
                head(thisP)   = zbottom(thisP) + fullDepth(thisP) + SlotDepth(thisP)
            end where 
        end if

        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_slot_adjustments
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
        !% The following initializations are unknown as of 20211215
        !% Preissmann_Celerity, SlotVolume, SlotArea, SlotWidth, SlotDepth
        !% SmallVolume_xxx
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
            elemR(thisP,er_HydDepth) = elemR(thisP,er_Depth)
            elemR(thisP,er_ell)      = elemR(thisP,er_Depth)
        end if
        !%------------------------------------------------------------------
        !% Closing
    end subroutine geo_JM_values
!%
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module geometry