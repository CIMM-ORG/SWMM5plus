module geometry

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use rectangular_channel
    use rectangular_conduit
    use trapezoidal_channel
    use triangular_channel
    use rectangular_triangular_conduit
    use rect_round_conduit
    use arch_conduit
    use basket_handle_conduit
    use mod_basket_conduit
    use catenary_conduit
    use egg_shaped_conduit
    use circular_conduit
    use semi_circular_conduit
    use filled_circular_conduit
    use semi_elliptical_conduit
    use gothic_conduit
    use horiz_ellipse_conduit
    use vert_ellipse_conduit
    use horse_shoe_conduit
    use irregular_channel
    use parabolic_channel
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
    public :: geo_sectionfactor_from_depth_singular
    public :: geo_Qcritical_from_depth_singular
    public :: geo_criticaldepth_singular
    public :: geo_normaldepth_singular
    public :: geo_assign_JB  !BRHbugfix 20210813
    public :: geo_topwidth_from_depth
    public :: geo_hyddepth_from_depth_singular
    public :: geo_topwidth_from_depth_singular
    public :: geo_area_from_depth_singular
    public :: geo_perimeter_from_depth_singular
    public :: geo_ell_singular

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
            integer, pointer :: thisColP_JM, thisColP_JB
            integer, pointer :: thisColP_Closed_CC, thisColP_Closed_JB, thisColP_Closed_JM
            logical :: isreset
            integer, allocatable :: tempP(:) !% debugging
            character(64) :: subroutine_name = 'geometry_toplevel'
        !%-----------------------------------------------------------------------------
        !% Preliminaries
            !!if (crashYN) return
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
                    thisColP_Closed_CC     => col_elemP(ep_CC_Closed_Elements)
                    thisColP_Closed_JB     => col_elemP(ep_Closed_JB_Elements)
                    thisColP_Closed_JM     => col_elemP(ep_JM_Closed_Elements)
                case (ETM)
                    elemPGx                => elemPGetm(:,:)
                    npack_elemPGx          => npack_elemPGetm(:)
                    col_elemPGx            => col_elemPGetm(:)
                    thisColP_JM            => col_elemP(ep_JM_ETM)
                    thisColP_JB            => col_elemP(ep_JB_ETM)
                    thisColP_surcharged    => col_elemP(ep_Surcharged_ETM)
                    thisColP_NonSurcharged => col_elemP(ep_NonSurcharged_ETM)
                    thisColP_all           => col_elemP(ep_ETM)
                    thisColP_Closed_CC     => col_elemP(ep_CC_Closed_Elements)
                    thisColP_Closed_JB     => col_elemP(ep_Closed_JB_Elements)
                    thisColP_Closed_JM     => col_elemP(ep_JM_Closed_Elements)
                case (AC)
                    elemPGx                => elemPGac(:,:)
                    npack_elemPGx          => npack_elemPGac(:)
                    col_elemPGx            => col_elemPGac(:)
                    thisColP_JM            => col_elemP(ep_JM_AC)
                    thisColP_JB            => col_elemP(ep_JB_AC)
                    thisColP_surcharged    => col_elemP(ep_Surcharged_AC)
                    thisColP_NonSurcharged => col_elemP(ep_NonSurcharged_AC)
                    thisColP_all           => col_elemP(ep_AC)
                    thisColP_Closed_CC     => col_elemP(ep_CC_Closed_Elements)
                    thisColP_Closed_JB     => col_elemP(ep_Closed_JB_Elements)
                    thisColP_Closed_JM     => col_elemP(ep_JM_Closed_Elements)
                case default
                    print *, 'CODE ERROR: time march type unknown for # ', whichTM
                    print *, 'which has key ',trim(reverseKey(whichTM))
                    call util_crashpoint(7389)
                    !return
                    !stop 7389
            end select
            call util_crashstop(49872)
        !%-----------------------------------------------------------------------------
        !% STATUS: at this point we know volume on Non-surcharged CC, JM,
        !% elements and head on all surcharged CC, JM elements

            ! call util_CLprint ('in geometry at top')    

        !% --- assign all geometry for surcharged elements CC, JM
        !%     Note: not used in Preissmann Slot
        call geo_surcharged (thisColP_surcharged)

            ! call util_CLprint ('in geometry before adjust_limit_by_zerovalues') 

        !% --- reset all zero or near-zero volumes in non-surcharged CC, JM
        call adjust_limit_by_zerovalues (er_Volume, setting%ZeroValue%Volume, thisColP_NonSurcharged, .true.)

            ! call util_CLprint ('in geometry before geo_depth_from_volume') 

        !% --- compute the depth on all non-surcharged elements of CC, JM
        call geo_depth_from_volume (elemPGx, npack_elemPGx, col_elemPGx)

            ! call util_CLprint ('in geometry before adjust_limit_by_zerovalues (2)') 

        !% reset all zero or near-zero depths in non-surcharged CC and JM
        call adjust_limit_by_zerovalues (er_Depth, setting%ZeroValue%Depth, thisColP_NonSurcharged, .false.)

            ! call util_CLprint ('in geometry before geo_head_from_depth') 

        !% --- compute the head on all non-surcharged elements of CC and JM
        !%     This sets head consistent with depth
        call geo_head_from_depth (thisColP_NonSurcharged)
 
            ! call util_CLprint ('in geometry before geo_limit_incipient_surcharge (Volume)') 

        !% --- limit volume for incipient surcharge. This is done after depth is computed
        !%     so that the "depth" algorithm can include depths greater than fulldepth
        !%     as a way to handle head for incipient surcharge.
        call geo_limit_incipient_surcharge (er_Volume, er_FullVolume, thisColP_NonSurcharged,.true.) !% 20220124brh

            ! call util_CLprint ('in geometry before geo_limit_incipient_surcharge (Depth)')  

        !% limit depth for incipient surcharged. This is done after head is computed
        !% so that the depth algorithm can include depths greater than fulldepth to
        !% handle incipient surcharge
        !call geo_limit_incipient_surcharge (er_Depth, er_FullDepth, thisColP_NonSurcharged)
        call geo_limit_incipient_surcharge (er_Depth, er_FullDepth, thisColP_NonSurcharged,.false.) !% 20220124brh

            ! call util_CLprint ('in geometry before geo_assign_JB') 

        !% STATUS: at this point we know depths and heads in all CC, JM elements
        !% (surcharged and nonsurcharged) with limiters for conduit depth and zero depth
           
        !% assign the head, depth, geometry on junction branches JB based on JM head
        call geo_assign_JB (whichTM, thisColP_JM)

            ! call util_CLprint ('in geometry before geo_area_from_volume')  

        !% STATUS at this point we know geometry on all JB and all surcharged, with
        !% depth, head, volume on all non-surcharged or incipient surcharge.

        !% compute area from volume for CC, JM nonsurcharged
        call geo_area_from_volume (thisColP_NonSurcharged)

            ! call util_CLprint ('in geometry before adjust_limit_by_zerovalues') 

        !% reset all zero or near-zero areas in non-surcharged CC and JM
        call adjust_limit_by_zerovalues (er_Area, setting%ZeroValue%Area, thisColP_NonSurcharged, .false.)

            ! call util_CLprint ('in geometry before topwidth_from_depth')   

        !% compute topwidth from depth for all CC, JM nonsurcharged
        call geo_topwidth_from_depth (elemPGx, npack_elemPGx, col_elemPGx)

            ! call util_CLprint ('in geometry before adjust_limit_by_zerovalues') 

        !% reset all zero or near-zero topwidth in non-surcharged CC and JM
        !% but do not change the eYN(:,eYN_isZeroDepth) mask
        call adjust_limit_by_zerovalues (er_Topwidth, setting%ZeroValue%Topwidth, thisColP_NonSurcharged, .false.)

            ! call util_CLprint ('in geometry before perimeter_from_depth') 

        !% compute perimeter from maximum depth for all CC, JM nonsurcharged
        call geo_perimeter_from_depth (elemPGx, npack_elemPGx, col_elemPGx)

            ! call util_CLprint ('in geometry before hyddepth_from_depth') 

        !% compute hyddepth
        call geo_hyddepth_from_depth (elemPGx, npack_elemPGx, col_elemPGx)

            ! call util_CLprint ('in geometry before hydradius_from_area_perimeter')   

        !% compute hydradius  (applies to all nonsurcharged)
        call geo_hydradius_from_area_perimeter (thisColP_NonSurcharged)

            ! call util_CLprint ('in geometry before ell_from_head') 

        !% the modified hydraulic depth "ell" is used for AC computations and
        !% for Froude number computations on all elements, whether ETM or AC.
        call geo_ell_from_head (thisColP_all)

        !% correct the head by adding zbottom and modified hydraulic depth
        call geo_head_from_ell (thisColP_all)

            ! call util_CLprint ('in geometry before slot_adjustments') 

        !% make adjustments for slots on closed elements only for ETM
        if (whichTM .eq. ETM) then
            call geo_CC_slot_adjustments (thisColP_Closed_CC)

            call geo_JB_slot_computation_ETM(thisColP_JM)
        end if
            ! call util_CLprint ('in geometry before JM_values') 

        !% Set JM values that are not otherwise defined
        call geo_JM_values ()

            ! call util_CLprint ('in geometry before dHdA') 

        !% compute the dHdA that are only for AC nonsurcharged
        if (whichTM .ne. ETM) then
            call geo_dHdA (ep_NonSurcharged_AC)
        end if

            ! call util_CLprint ('in geometry at end') 

        call util_crashstop(322983)

        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geometry_toplevel
!%
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
        !print *, 'in ',trim(subroutine_name)    
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
        !%------------------------------------------------------------------
        !% Aliases
            eIdx      => uniformTableI(UT_idx,uti_elem_idx)
            thisTable => uniformTableDataR(UT_idx,:,utd_Qcrit_depth_nonuniform)
        !%------------------------------------------------------------------
        !% --- normalize the critical flowrate
        normFlowrate = abs(elemR(eIdx,er_Flowrate) / uniformTableR(UT_idx,utr_QcritMax))

        !% --- lookup the normalized critical depth for this critical flow
        outvalue = xsect_table_lookup_singular (normFlowrate, thistable)

        !% --- return depth to physical value
        outvalue = outvalue * uniformTableR(UT_idx,utr_DepthMax)

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
            eIdx      => uniformTableI(UT_idx,uti_elem_idx)
            thisTable => uniformTableDataR(UT_idx,:,utd_SF_depth_nonuniform)  !% element index
        !%------------------------------------------------------------------
        !% --- section factor for the associated element
        sectionFactor = elemR(eIdx,er_Flowrate) * elemR(eIdx,er_ManningsN) / elemR(eIdx,er_BottomSlope)

        !print *, 'sectionFactor ',sectionFactor

        !% --- if flow is negative on a positive slope, or flow is positive on a negative slope,
        !%     then the section factor is negative, which implies an infinite normal depth
        if (sectionFactor .le. zeroR) then
            outvalue = setting%Limiter%NormalDepthInfinite
            return
        end if

        !% --- normalize the section factor
        normSF   = sectionFactor / uniformTableR(UT_idx,utr_SFmax)

        !print *, 'normSF ',normSF

        !% --- lookup the normalized normal depth
        outvalue = xsect_table_lookup_singular(normSF,thisTable)


        !print *, 'outvalue 1 ',outvalue

        !% --- return normal depth to physical value
        outvalue = outvalue * uniformTableR(UT_idx,utr_DepthMax)

        !print *, 'outvalue 2 ',outvalue
    
    end function geo_normaldepth_singular
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
            real(8), pointer :: fullhyddepth(:), thisTable(:,:)
            real(8), pointer :: slotDepth(:), slotVolume(:), overflow(:)
            real(8), pointer :: grav  
            logical, pointer :: isSlot(:)     

            real(8) :: depthnorm, zeroHydRadius
            integer :: tB, ii, kk

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
            ellMax        => elemR(:,er_ell_max)
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
            fullArea      => elemR(:,er_FullArea)
            fulldepth     => elemR(:,er_FullDepth)
            fullhyddepth  => elemR(:,er_FullHydDepth)
            fullperimeter => elemR(:,er_FullPerimeter)
            overflow      => elemR(:,er_VolumeOverFlow)
            slotDepth     => elemR(:,er_SlotDepth)
            slotVolume    => elemR(:,er_SlotVolume)
            Kfac          => elemSR(:,esr_JunctionBranch_Kfactor)
            BranchExists  => elemSI(:,esi_JunctionBranch_Exists)
            thisSolve     => elemI(:,ei_tmType)
            isSlot        => elemYN(:,eYN_isSlot)
            grav => setting%Constant%gravity
        !%------------------------------------------------------------------

        if (Npack > 0) then
            thisP  => elemP(1:Npack,thisColP_JM)

            !% cycle through the all the main junctions and each of its branches
            do ii=1,Npack
                
                tM => thisP(ii) !% junction main ID

                !% if a slot present, add the slot depth and volume back to JM
                if (isSlot(tM)) then
                    volume(tM)   = volume(tM)  + SlotVolume(tM) 
                    depth(tM)    = depth(tM)   + SlotDepth(tM)
                    head(tM)     = head(tM)    + SlotDepth(tM)
                    ell(tM)      = ellMax(tM)
                    Overflow(tM) = zeroR
                end if 

                !% only execute for whichTM of ALL or thisSolve (of JM) matching input whichTM
                if ((whichTM == ALLtm) .or. (thisSolve(tM) == whichTM)) then
                    !% cycle through the possible junction branches
                    do kk=1,max_branch_per_node
                        
                        tB = tM + kk !% junction branch ID

                        ! print *, kk, tB
                        ! print *, BranchExists(tB)

                        if (BranchExists(tB) == 1) then
                            !% only when a branch exists.
                            ! print *, head(tM), zBtm(tB)
                            ! print *, kk, branchsign(kk)
                            ! print *, velocity(tB)
                            ! print *, Kfac(tB)
                            if ( head(tM) > zBtm(tB) ) then
                                !% for main head above branch bottom entrance use a head
                                !% loss approach. The branchsign and velocity sign ensure
                                !% the headloss is added to an inflow and subtracted at
                                !% an outflow
                                !% Note this is a time-lagged velocity as the JB velocity
                                !% is not updated until after face interpolation                                
                                head(tB) = head(tM) + branchsign(kk) * sign(oneR,velocity(tB)) &
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
                                    *(velocity(tB)**twoR) / (grav)   !% 20220307 brh ADDED 2 to factor -- removed 20220615
                            end if

                            !% HACK -- the above uses a Froude number argument for head(TM) < zBtm(tB)
                            !%      however, when the JB is surcharged we probably should be using the
                            !%      K factor approach and require K=1.
                           
                            !% compute provisional depth
                            depth(tB) = head(tB) - zBtm(tB)

                            ! print *, 'in geo_assign_JB  ',trim(reverseKey(elemI(tB,ei_geometryType)))
                            ! print *, 'depth ',depth(tB), fulldepth(tB), setting%ZeroValue%Depth
                            
                            if (depth(tB) .ge. fulldepth(tB)) then
                                !% surcharged or incipient surcharged
                                depth(tB)     = fulldepth(tB)
                                area(tB)      = fullArea(tB)
                                hyddepth(tB)  = fullhyddepth(tB)
                                perimeter(tB) = fullperimeter(tB)
                                topwidth(tB)  = setting%ZeroValue%Topwidth
                                hydRadius(tB) = fulldepth(tB) / fullperimeter(tB)
                                dHdA(tB)      = oneR / setting%ZeroValue%Topwidth
                                ell(tB)       = geo_ell_singular(tB)

                                ! write(*,"(A,i5,10f12.5)") 'AAA ell ',tB, ell(tB), depth(tB), hydDepth(tB), fulldepth(tB)

                            elseif ((depth(tB) < setting%ZeroValue%Depth) .and. (setting%ZeroValue%UseZeroValues)) then
                                !% negligible depth is treated with ZeroValues
                                depth(tB)     = setting%ZeroValue%Depth
                                area(tB)      = setting%ZeroValue%Area
                                topwidth(tB)  = setting%ZeroValue%Topwidth
                                hyddepth(tB)  = setting%ZeroValue%Depth !% setting%ZeroValue%Area / topwidth(tB) 20220712brh
                                perimeter(tB) = topwidth(tB) + setting%ZeroValue%Depth
                                hydRadius(tB) = setting%ZeroValue%Area / perimeter(tB)
                                dHdA(tB)      = oneR / topwidth(tB)
                                ell(tB)       = setting%ZeroValue%Depth !%hydDepth(tB)  20220712 brh

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

                                ! write(*,"(A,i5,10f12.5)") 'CCC ell ',tB, ell(tB), depth(tB), hydDepth(tB), fulldepth(tB)

                            else
                                !% not surcharged and non-negligible depth
                                select case (elemI(tB,ei_geometryType))
                                case (rectangular)
                                    area(tB)     = rectangular_area_from_depth_singular      (tB, depth(tB))
                                    topwidth(tB) = rectangular_topwidth_from_depth_singular  (tB, depth(tB))
                                    hydDepth(tB) = rectangular_hyddepth_from_depth_singular  (tB, depth(tB))
                                    perimeter(tB)= rectangular_perimeter_from_depth_singular (tB, depth(tB))
                                    hydRadius(tB)= rectangular_hydradius_from_depth_singular (tB, depth(tB))
                                    ell(tB)      = hydDepth(tB) !geo_ell_singular (tB) !BRHbugfix 20210812 simpler for rectangle
                                    dHdA(tB)     = oneR / topwidth(tB)
                                
                                case (rectangular_closed)
                                    area(tB)     = rectangular_closed_area_from_depth_singular      (tB, depth(tB))
                                    topwidth(tB) = rectangular_closed_topwidth_from_depth_singular  (tB, depth(tB))
                                    hydDepth(tB) = rectangular_closed_hyddepth_from_depth_singular  (tB, depth(tB))
                                    perimeter(tB)= rectangular_closed_perimeter_from_depth_singular (tB, depth(tB))
                                    hydRadius(tB)= rectangular_closed_hydradius_from_depth_singular (tB, depth(tB))
                                    ell(tB)      = hydDepth(tB) !geo_ell_singular (tB) !BRHbugfix 20210812 simpler for rectangle
                                    dHdA(tB)     = oneR / topwidth(tB) 

                                !    write(*,"(A,i5,10f12.5)") 'DDD ell ',tB, ell(tB), depth(tB), hydDepth(tB), fulldepth(tB)

                                case (triangular)
                                    area(tB)     = triangular_area_from_depth_singular      (tB,depth(tB))
                                    topwidth(tB) = triangular_topwidth_from_depth_singular  (tB,depth(tB))
                                    hydDepth(tB) = triangular_hyddepth_from_depth_singular  (tB,depth(tB))
                                    perimeter(tB)= triangular_perimeter_from_depth_singular (tB,depth(tB))
                                    hydRadius(tB)= triangular_hydradius_from_depth_singular (tB,depth(tB))
                                    ell(tB)      = geo_ell_singular (tB) 
                                    dHdA(tB)     = oneR / topwidth(tB)

                                !    write(*,"(A,i5,10f12.5)") 'EEE ell ',tB, ell(tB), depth(tB), hydDepth(tB), fulldepth(tB)
                                    
                                case (rect_triang)                                    
                                    area(tB)     = rectangular_triangular_area_from_depth_singular      (tB,depth(tB))
                                    topwidth(tB) = rectangular_triangular_topwidth_from_depth_singular  (tB,depth(tB))
                                    hydDepth(tB) = rectangular_triangular_hyddepth_from_depth_singular  (tB,depth(tB))
                                    perimeter(tB)= rectangular_triangular_perimeter_from_depth_singular (tB,depth(tB))
                                    hydRadius(tB)= rectangular_triangular_hydradius_from_depth_singular (tB,depth(tB))
                                    ell(tB)      = hydDepth(tB) !geo_ell_singular (tB) !BRHbugfix 20210812 simpler for rect_triang
                                    dHdA(tB)     = oneR / topwidth(tB)
                                
                                case (rect_round)                                    
                                    area(tB)     = rect_round_area_from_depth_singular        (tB,depth(tB))
                                    topwidth(tB) = rect_round_topwidth_from_depth_singular    (tB,depth(tB))
                                    hydDepth(tB) = rect_round_hyddepth_from_topwidth_singular (tB,topwidth(tB),depth(tB))
                                    perimeter(tB)= rect_round_perimeter_from_depth_singular   (tB,depth(tB))
                                    hydRadius(tB)= rect_round_hydradius_from_depth_singular   (tB,depth(tB))
                                    ell(tB)      = hydDepth(tB) !geo_ell_singular (tB) !BRHbugfix 20210812 simpler for rect_triang
                                    dHdA(tB)     = oneR / topwidth(tB)
                                
                                case (basket_handle)                                    
                                    area(tB)     = basket_handle_area_from_depth_singular        (tB,depth(tB))
                                    topwidth(tB) = basket_handle_topwidth_from_depth_singular    (tB,depth(tB))
                                    hydDepth(tB) = basket_handle_hyddepth_from_topwidth_singular (tB,topwidth(tB),depth(tB))
                                    perimeter(tB)= basket_handle_perimeter_from_depth_singular   (tB,depth(tB))
                                    hydRadius(tB)= basket_handle_hydradius_from_depth_singular   (tB,depth(tB))
                                    ell(tB)      = hydDepth(tB) !geo_ell_singular (tB) !BRHbugfix 20210812 simpler for rect_triang
                                    dHdA(tB)     = oneR / topwidth(tB)
                                
                                case (arch)                                    
                                    area(tB)     = arch_area_from_depth_singular        (tB,depth(tB))
                                    topwidth(tB) = arch_topwidth_from_depth_singular    (tB,depth(tB))
                                    hydDepth(tB) = arch_hyddepth_from_topwidth_singular (tB,topwidth(tB),depth(tB))
                                    perimeter(tB)= arch_perimeter_from_depth_singular   (tB,depth(tB))
                                    hydRadius(tB)= arch_hydradius_from_depth_singular   (tB,depth(tB))
                                    ell(tB)      = hydDepth(tB) !geo_ell_singular (tB) !BRHbugfix 20210812 simpler for rect_triang
                                    dHdA(tB)     = oneR / topwidth(tB)
                                
                                case (horiz_ellipse)                                    
                                    area(tB)     = horiz_ellipse_area_from_depth_singular        (tB,depth(tB))
                                    topwidth(tB) = horiz_ellipse_topwidth_from_depth_singular    (tB,depth(tB))
                                    hydDepth(tB) = horiz_ellipse_hyddepth_from_topwidth_singular (tB,topwidth(tB),depth(tB))
                                    perimeter(tB)= horiz_ellipse_perimeter_from_depth_singular   (tB,depth(tB))
                                    hydRadius(tB)= horiz_ellipse_hydradius_from_depth_singular   (tB,depth(tB))
                                    ell(tB)      = hydDepth(tB) !geo_ell_singular (tB) !BRHbugfix 20210812 simpler for rect_triang
                                    dHdA(tB)     = oneR / topwidth(tB)
                                
                                case (vert_ellipse)                                    
                                    area(tB)     = vert_ellipse_area_from_depth_singular        (tB,depth(tB))
                                    topwidth(tB) = vert_ellipse_topwidth_from_depth_singular    (tB,depth(tB))
                                    hydDepth(tB) = vert_ellipse_hyddepth_from_topwidth_singular (tB,topwidth(tB),depth(tB))
                                    perimeter(tB)= vert_ellipse_perimeter_from_depth_singular   (tB,depth(tB))
                                    hydRadius(tB)= vert_ellipse_hydradius_from_depth_singular   (tB,depth(tB))
                                    ell(tB)      = hydDepth(tB) !geo_ell_singular (tB) !BRHbugfix 20210812 simpler for rect_triang
                                    dHdA(tB)     = oneR / topwidth(tB)
                                
                                case (eggshaped)                                    
                                    area(tB)     = egg_shaped_area_from_depth_singular        (tB,depth(tB))
                                    topwidth(tB) = egg_shaped_topwidth_from_depth_singular    (tB,depth(tB))
                                    hydDepth(tB) = egg_shaped_hyddepth_from_topwidth_singular (tB,topwidth(tB),depth(tB))
                                    perimeter(tB)= egg_shaped_perimeter_from_depth_singular   (tB,depth(tB))
                                    hydRadius(tB)= egg_shaped_hydradius_from_depth_singular   (tB,depth(tB))
                                    ell(tB)      = hydDepth(tB) !geo_ell_singular (tB) !BRHbugfix 20210812 simpler for rect_triang
                                    dHdA(tB)     = oneR / topwidth(tB)
                                
                                case (horseshoe)                                    
                                    area(tB)     = horse_shoe_area_from_depth_singular        (tB,depth(tB))
                                    topwidth(tB) = horse_shoe_topwidth_from_depth_singular    (tB,depth(tB))
                                    hydDepth(tB) = horse_shoe_hyddepth_from_topwidth_singular (tB,topwidth(tB),depth(tB))
                                    perimeter(tB)= horse_shoe_perimeter_from_depth_singular   (tB,depth(tB))
                                    hydRadius(tB)= horse_shoe_hydradius_from_depth_singular   (tB,depth(tB))
                                    ell(tB)      = hydDepth(tB) !geo_ell_singular (tB) !BRHbugfix 20210812 simpler for rect_triang
                                    dHdA(tB)     = oneR / topwidth(tB)
                                
                                case (catenary)                                    
                                    area(tB)     = catenary_area_from_depth_singular        (tB,depth(tB))
                                    topwidth(tB) = catenary_topwidth_from_depth_singular    (tB,depth(tB))
                                    hydDepth(tB) = catenary_hyddepth_from_topwidth_singular (tB,topwidth(tB),depth(tB))
                                    perimeter(tB)= catenary_perimeter_from_depth_singular   (tB,depth(tB))
                                    hydRadius(tB)= catenary_hydradius_from_depth_singular   (tB,depth(tB))
                                    ell(tB)      = hydDepth(tB) !geo_ell_singular (tB) !BRHbugfix 20210812 simpler for rect_triang
                                    dHdA(tB)     = oneR / topwidth(tB)
                                
                                case (gothic)                                    
                                    area(tB)     = gothic_area_from_depth_singular        (tB,depth(tB))
                                    topwidth(tB) = gothic_topwidth_from_depth_singular    (tB,depth(tB))
                                    hydDepth(tB) = gothic_hyddepth_from_topwidth_singular (tB,topwidth(tB),depth(tB))
                                    perimeter(tB)= gothic_perimeter_from_depth_singular   (tB,depth(tB))
                                    hydRadius(tB)= gothic_hydradius_from_depth_singular   (tB,depth(tB))
                                    ell(tB)      = hydDepth(tB) !geo_ell_singular (tB) !BRHbugfix 20210812 simpler for rect_triang
                                    dHdA(tB)     = oneR / topwidth(tB)
                                
                                case (mod_basket)                                    
                                    area(tB)     = mod_basket_area_from_depth_singular        (tB,depth(tB))
                                    topwidth(tB) = mod_basket_topwidth_from_depth_singular    (tB,depth(tB))
                                    hydDepth(tB) = mod_basket_hyddepth_from_topwidth_singular (tB,topwidth(tB),depth(tB))
                                    perimeter(tB)= mod_basket_perimeter_from_depth_singular   (tB,depth(tB))
                                    hydRadius(tB)= mod_basket_hydradius_from_depth_singular   (tB,depth(tB))
                                    ell(tB)      = hydDepth(tB) !geo_ell_singular (tB) !BRHbugfix 20210812 simpler for rect_triang
                                    dHdA(tB)     = oneR / topwidth(tB)
                                
                                case (semi_elliptical)                                    
                                    area(tB)     = semi_elliptical_area_from_depth_singular        (tB,depth(tB))
                                    topwidth(tB) = semi_elliptical_topwidth_from_depth_singular    (tB,depth(tB))
                                    hydDepth(tB) = semi_elliptical_hyddepth_from_topwidth_singular (tB,topwidth(tB),depth(tB))
                                    perimeter(tB)= semi_elliptical_perimeter_from_depth_singular   (tB,depth(tB))
                                    hydRadius(tB)= semi_elliptical_hydradius_from_depth_singular   (tB,depth(tB))
                                    ell(tB)      = hydDepth(tB) !geo_ell_singular (tB) !BRHbugfix 20210812 simpler for rect_triang
                                    dHdA(tB)     = oneR / topwidth(tB)
                                
                                case (trapezoidal)
                                    area(tB)     = trapezoidal_area_from_depth_singular      (tB,depth(tB))
                                    topwidth(tB) = trapezoidal_topwidth_from_depth_singular  (tB,depth(tB))
                                    hydDepth(tB) = trapezoidal_hyddepth_from_depth_singular  (tB,depth(tB))
                                    perimeter(tB)= trapezoidal_perimeter_from_depth_singular (tB,depth(tB))
                                    hydRadius(tB)= trapezoidal_hydradius_from_depth_singular (tB,depth(tB))
                                    ell(tB)      = geo_ell_singular (tB) 
                                    dHdA(tB)     = oneR / topwidth(tB)

                                !    write(*,"(A,i5,10f12.5)") 'FFF ell ',tB, ell(tB), depth(tB), hydDepth(tB), fulldepth(tB)

                                case (circular)
                                    area(tB)     = circular_area_from_depth_singular          (tB,depth(tB))
                                    topwidth(tB) = circular_topwidth_from_depth_singular      (tB,depth(tB))
                                    hydDepth(tB) = circular_hyddepth_from_topwidth_singular   (tB,topwidth(tB),depth(tB))
                                    hydRadius(tB)= circular_hydradius_from_depth_singular     (tB,depth(tB))
                                    perimeter(tB)= circular_perimeter_from_hydradius_singular (tB,hydRadius(tB))
                                    ell(tB)      = geo_ell_singular (tB) 
                                    dHdA(tB)     = oneR / topwidth(tB)
                                
                                case (semi_circular)                                    
                                    area(tB)     = semi_circular_area_from_depth_singular        (tB,depth(tB))
                                    topwidth(tB) = semi_circular_topwidth_from_depth_singular    (tB,depth(tB))
                                    hydDepth(tB) = semi_circular_hyddepth_from_topwidth_singular (tB,topwidth(tB),depth(tB))
                                    perimeter(tB)= semi_circular_perimeter_from_depth_singular   (tB,depth(tB))
                                    hydRadius(tB)= semi_circular_hydradius_from_depth_singular   (tB,depth(tB))
                                    ell(tB)      = hydDepth(tB) !geo_ell_singular (tB) !BRHbugfix 20210812 simpler for rect_triang
                                    dHdA(tB)     = oneR / topwidth(tB)

                                    ! write(*,"(A,i5,10f12.5)"), 'GGG ell ',tB, ell(tB), depth(tB), hydDepth(tB), fulldepth(tB)
                                case (filled_circular)
                                    area(tB)     = filled_circular_area_from_depth_singular          (tB,depth(tB))
                                    topwidth(tB) = filled_circular_topwidth_from_depth_singular      (tB,depth(tB))
                                    hydDepth(tB) = filled_circular_hyddepth_from_topwidth_singular   (tB,topwidth(tB),depth(tB))
                                    hydRadius(tB)= filled_circular_hydradius_from_depth_singular     (tB,depth(tB))
                                    perimeter(tB)= filled_circular_perimeter_from_hydradius_singular (tB,hydRadius(tB))
                                    ell(tB)      = geo_ell_singular (tB) 
                                    dHdA(tB)     = oneR / topwidth(tB)

                                case (parabolic)
                                    area(tB)     = parabolic_area_from_depth_singular      (tB, depth(tB))
                                    topwidth(tB) = parabolic_topwidth_from_depth_singular  (tB, depth(tB))
                                    hydDepth(tB) = parabolic_hyddepth_from_depth_singular  (tB, depth(tB))
                                    perimeter(tB)= parabolic_perimeter_from_depth_singular (tB, depth(tB))
                                    hydRadius(tB)= parabolic_hydradius_from_depth_singular (tB, depth(tB))
                                    ell(tB)      = hydDepth(tB) !geo_ell_singular (tB) !BRHbugfix 20210812 simpler for rectangle
                                    dHdA(tB)     = oneR / topwidth(tB)

                                case (irregular)
                                    area(tB)    = irregular_geometry_from_depth_singular ( &
                                                        tB,tt_area, depth(tB), setting%ZeroValue%Depth)

                                    topwidth(tB) = irregular_geometry_from_depth_singular ( &
                                                        tB,tt_width, depth(tB), setting%ZeroValue%TopWidth)

                                    zeroHydRadius = setting%ZeroValue%Area / (setting%ZeroValue%TopWidth + setting%ZeroValue%Depth)
                                    hydRadius(tB) = irregular_geometry_from_depth_singular ( &
                                                        tB,tt_hydradius, depth(tB), zeroHydRadius)                

                                    hydDepth(tB)  = area(tB) / topwidth(tB)     
                                    
                                    perimeter(tB) = area(tB) / hydRadius(tB)
                                    ell(tB)       = hydDepth(tB)  !% HACK -- assumes irregular is continuously-increasing in width
                                    dHdA(tB)      = oneR / topwidth(tB)

                                !    write(*,"(A,i5,10f12.5)") 'HHH ell ',tB, ell(tB), depth(tB), hydDepth(tB), fulldepth(tB)

                                    ! !% get the transect by depth table 
                                    ! thisTable => transectTableDepthR(elemI(tB,ei_transect_idx),:,:)
                                    ! depthnorm     = depth(tB)/fulldepth(tB)

                                    ! area(tB)      = xsect_table_lookup_singular (depthnorm, thisTable(:,tt_area))
                                    ! !% xsect quadratic interp for small values can produce zero
                                    ! area(tB)      = max(area(tB), setting%ZeroValue%Area)

                                    ! topwidth(tB)  = xsect_table_lookup_singular (depthnorm, thisTable(:,tt_width))
                                    ! !% xsect quadratic interp for small values can produce zero
                                    ! topwidth(tB)  = max(topwidth(tB),setting%ZeroValue%TopWidth)

                                    ! hydRadius(tB) = xsect_table_lookup_singular (depthnorm, thisTable(:,tt_hydradius))
                                    ! !% xsect quadratic interp for small values can produce zero
                                    ! hydRadius(tB) = max(hydRadius(tB),setting%ZeroValue%Area / (setting%ZeroValue%TopWidth + setting%ZeroValue%Depth))
                                    
                                    

                                case default
                                    print *, 'CODE ERROR: geometry type unknown for # ', elemI(tB,ei_geometryType)
                                    print *, 'which has key ',trim(reverseKey(elemI(tB,ei_geometryType)))
                                    print *, 'in ',trim(subroutine_name)
                                    call util_crashpoint(399848)
                                    !return
                                    !stop 399848
                                end select
                            end if

                            ! print *, 'at bottom '

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
        !!if (crashYN) return
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
        !% This solves nonsurcharged CCJMJB elements because of PGx arrays
        !% The elemPGx determines whether this is ALLtm, ETM or AC elements
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: elemPGx(:,:), npack_elemPGx(:), col_elemPGx(:)
            integer, pointer :: Npack, thisCol
            character(64) :: subroutine_name = 'geo_depth_from_volume'
        !%-------------------------------------------------------------------
        !% Preliminaries
            !!if (crashYN) return
            if (setting%Debug%File%geometry) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------    
        !% cycle through different geometries

        ! call util_CLprint('start of geo depth from volume')        

        !% --- RECTANGULAR CC
        thisCol => col_elemPGx(epg_CC_rectangular_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call rectangular_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% --- RECTANGULAR CLOSED
        thisCol => col_elemPGx(epg_CC_rectangular_closed_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call rectangular_closed_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        ! call util_CLprint('after rectangular') 

        !% --- TRAPEZOIDAL CC
        thisCol => col_elemPGx(epg_CC_trapezoidal_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call trapezoidal_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        ! call util_CLprint('after trapezoidal') 

        !% --- TRIANGULAR CC
        thisCol => col_elemPGx(epg_CC_triangular_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call triangular_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        ! call util_CLprint('after triangular') 

        !% --- CIRCULAR CC
        thisCol => col_elemPGx(epg_CC_circular_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call circular_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% --  SEMI-CIRCULAR
        thisCol => col_elemPGx(epg_CC_semi_circular_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call semi_circular_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% --- FILLED CIRCULAR CC
        thisCol => col_elemPGx(epg_CC_filled_circular_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call filled_circular_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% -- PARABOLIC
        thisCol => col_elemPGx(epg_CC_parabolic_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call parabolic_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% -- RECT_TRIANG
        thisCol => col_elemPGx(epg_CC_rectangular_triangular_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call rectangular_triangular_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% --  BASKET_HANDLE
        thisCol => col_elemPGx(epg_CC_basket_handle_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call basket_handle_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% --  ARCH
        thisCol => col_elemPGx(epg_CC_arch_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call arch_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% --  Horizontal Ellipse
        thisCol => col_elemPGx(epg_CC_horiz_ellipse_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call horiz_ellipse_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% --  Vertical Ellipse
        thisCol => col_elemPGx(epg_CC_vert_ellipse_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call vert_ellipse_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% --  EGG_SHAPED
        thisCol => col_elemPGx(epg_CC_egg_shaped_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call egg_shaped_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% --  HORSE_SHOE
        thisCol => col_elemPGx(epg_CC_horse_shoe_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call horse_shoe_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% --  CATENARY
        thisCol => col_elemPGx(epg_CC_catenary_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call catenary_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% --  GOTHIC
        thisCol => col_elemPGx(epg_CC_gothic_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call gothic_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% --  MODIFIED BASKET HANDLE
        thisCol => col_elemPGx(epg_CC_mod_basket_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call mod_basket_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% --  RECTANGULAR ROUND
        thisCol => col_elemPGx(epg_CC_rectangular_round_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call rect_round_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !% --  SEMI-ELLIPTICAL
        thisCol => col_elemPGx(epg_CC_semi_elliptical_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call semi_elliptical_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !call util_CLprint('after circular') 

        !% --- IRREGULAR
        thisCol => col_elemPGx(epg_CC_irregular_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call irregular_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        ! call util_CLprint('after irregular') 
        !% HACK Needs additional geometries

        !% JM with functional geometry
        thisCol => col_elemPGx(epg_JM_functionalStorage_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call storage_functional_depth_from_volume (elemPGx, Npack, thisCol)
            !call storage_implied_length(elemPGx, Npack, thisCol)
        end if

        ! call util_CLprint('after functional storage') 

        !% JM with tabular geomtery
        thisCol => col_elemPGx(epg_JM_tabularStorage_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call storage_tabular_depth_from_volume (elemPGx, Npack, thisCol)
            !call storage_implied_length(elemPGx, Npack, thisCol)
        end if

        !call util_CLprint('after tabular storage') 

        !% JM with implied storage (note that length is already defined)
        thisCol => col_elemPGx(epg_JM_impliedStorage_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call storage_implied_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !call util_CLprint('after implied storage') 
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
        !!if (crashYN) return
        Npack      => npack_elemP(thisColP)
        geovalue   => elemR(:,geocol)
        fullvalue  => elemR(:,fullcol)
        overflow   => elemR(:,er_VolumeOverFlow)  !% 20220124brh
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

            ! print *, 'in ',trim(subroutine_name),elemR(49,er_VolumeOverFlow)
            ! print *, geovalue(49), fullvalue(49), overflow(49)

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

        ! print *, 'end of ',trim(subroutine_name),elemR(48,er_VolumeOverFlow)

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
        !!if (crashYN) return
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

        ! print *, 'thisP in geo_head_from_depth'
        ! print *, thisP

        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_head_from_depth
!%
!%==========================================================================
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
        !!if (crashYN) return
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
    real(8) function geo_area_from_depth_singular &
        (idx, indepth) result (outvalue)
        !%------------------------------------------------------------------
        !% Descriptions:
        !% computes the area for a given depth of a single element
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(in)  :: indepth
            integer, intent(in)  :: idx
            character(64) :: subroutine_name = 'geo_area_from_depth_singular'
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        select case (elemI(idx,ei_geometryType))
            
        case (rectangular)
            outvalue = rectangular_area_from_depth_singular (idx, indepth)
        case (trapezoidal)
            outvalue = trapezoidal_area_from_depth_singular (idx, indepth)
        case (triangular)
            outvalue = triangular_area_from_depth_singular (idx, indepth)
        case (parabolic)
            outvalue = parabolic_area_from_depth_singular (idx, indepth)
        case (power_function)
            print *, 'CODE ERROR: area for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(33234)
        case (rect_triang)
            outvalue = rectangular_triangular_area_from_depth_singular (idx, indepth)
        case (rect_round )
            outvalue = rect_round_area_from_depth_singular (idx, indepth)
        case (mod_basket)
            outvalue = mod_basket_area_from_depth_singular (idx, indepth)
        case (irregular)
            outvalue = irregular_geometry_from_depth_singular (idx,tt_area, indepth, setting%ZeroValue%Depth)
        case (circular )
            outvalue = circular_area_from_depth_singular (idx, indepth)
        case (filled_circular)
            outvalue = filled_circular_area_from_depth_singular (idx, indepth)
        case (rectangular_closed)
            outvalue = rectangular_closed_area_from_depth_singular (idx, indepth)
        case (horiz_ellipse)
            outvalue = horiz_ellipse_area_from_depth_singular (idx, indepth)
        case (vert_ellipse)
            outvalue = vert_ellipse_area_from_depth_singular (idx, indepth)
        case (arch)
            outvalue = arch_area_from_depth_singular (idx, indepth)
        case (eggshaped)
            outvalue = egg_shaped_area_from_depth_singular (idx, indepth)
        case (horseshoe)
            outvalue = horse_shoe_area_from_depth_singular (idx, indepth)
        case (gothic)
            outvalue = gothic_area_from_depth_singular (idx, indepth)
        case (catenary)
            outvalue = catenary_area_from_depth_singular (idx, indepth)
        case (semi_elliptical)
            outvalue = semi_elliptical_area_from_depth_singular (idx, indepth)
        case (basket_handle)
            outvalue = basket_handle_area_from_depth_singular (idx, indepth)
        case (semi_circular)
            outvalue = semi_circular_area_from_depth_singular (idx, indepth)
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
        !!if (crashYN) return
        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% --- RECTANGULAR
        Npack => npack_elemPGx(epg_CC_rectangular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_nonsurcharged)
            call rectangular_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- RECTANGULAR CLOSED
        Npack => npack_elemPGx(epg_CC_rectangular_closed_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_closed_nonsurcharged)
            call rectangular_closed_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- TRAPEZOIDAL
        Npack => npack_elemPGx(epg_CC_trapezoidal_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_trapezoidal_nonsurcharged)
            call trapezoidal_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- TRIANGULAR
        Npack => npack_elemPGx(epg_CC_triangular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_triangular_nonsurcharged)
            call triangular_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- CIRCULAR
        Npack => npack_elemPGx(epg_CC_circular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_circular_nonsurcharged)
            call circular_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- SEMI-CIRCULAR
        Npack => npack_elemPGx(epg_CC_semi_circular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_semi_circular_nonsurcharged)
            call semi_circular_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- FILLED CIRCULAR
        Npack => npack_elemPGx(epg_CC_filled_circular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_filled_circular_nonsurcharged)
            call filled_circular_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- PARABOLIC
        Npack => npack_elemPGx(epg_CC_parabolic_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_parabolic_nonsurcharged)
            call parabolic_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- RECT_TRIANG
        Npack => npack_elemPGx(epg_CC_rectangular_triangular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_triangular_nonsurcharged)
            call rectangular_triangular_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- BASKET_HANDLE
        Npack => npack_elemPGx(epg_CC_basket_handle_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_basket_handle_nonsurcharged)
            call basket_handle_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- ARCH
        Npack => npack_elemPGx(epg_CC_arch_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_arch_nonsurcharged)
            call arch_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- HORIZONTAL ELLIPSE
        Npack => npack_elemPGx(epg_CC_horiz_ellipse_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_horiz_ellipse_nonsurcharged)
            call horiz_ellipse_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- VERTICAL ELLIPSE
        Npack => npack_elemPGx(epg_CC_vert_ellipse_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_vert_ellipse_nonsurcharged)
            call vert_ellipse_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- EGG_SHAPED
        Npack => npack_elemPGx(epg_CC_egg_shaped_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_egg_shaped_nonsurcharged)
            call egg_shaped_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- HORSE_SHOE
        Npack => npack_elemPGx(epg_CC_horse_shoe_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_horse_shoe_nonsurcharged)
            call horse_shoe_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- CATENARY
        Npack => npack_elemPGx(epg_CC_catenary_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_catenary_nonsurcharged)
            call catenary_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- GOTHIC
        Npack => npack_elemPGx(epg_CC_gothic_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_gothic_nonsurcharged)
            call gothic_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- MODIFIED BASKET HANDLE
        Npack => npack_elemPGx(epg_CC_mod_basket_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_mod_basket_nonsurcharged)
            call mod_basket_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- RECTANGULAR ROUND
        Npack => npack_elemPGx(epg_CC_rectangular_round_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_round_nonsurcharged)
            call rect_round_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- SEMI-ELLIPTICAL
        Npack => npack_elemPGx(epg_CC_semi_elliptical_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_semi_elliptical_nonsurcharged)
            call semi_elliptical_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- IRREGULAR
        Npack => npack_elemPGx(epg_CC_irregular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_irregular_nonsurcharged)
            call irregular_topwidth_from_depth (elemPGx, Npack, thisCol)
        end if
        !% HACK NEED OTHER GEOMETRIES
        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_topwidth_from_depth
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
            character(64) :: subroutine_name = 'geo_topwidth_from_depth_singular'
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        select case (elemI(idx,ei_geometryType))
            
        case (rectangular)
            outvalue = rectangular_topwidth_from_depth_singular  (idx, indepth)
        case (trapezoidal)
            outvalue = trapezoidal_topwidth_from_depth_singular (idx, indepth)
        case (triangular)
            outvalue = triangular_topwidth_from_depth_singular  (idx, indepth)
        case (parabolic)
            outvalue = parabolic_topwidth_from_depth_singular (idx, indepth)
        case (power_function)
            print *, 'CODE ERROR: topwidth for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (rect_triang)
            outvalue = rectangular_triangular_topwidth_from_depth_singular  (idx, indepth)
        case (rect_round )
            outvalue = rect_round_topwidth_from_depth_singular (idx, indepth)
        case (mod_basket)
            outvalue = mod_basket_topwidth_from_depth_singular (idx, indepth)
        case (irregular)
            outvalue = irregular_geometry_from_depth_singular (idx,tt_width, indepth, setting%ZeroValue%TopWidth)
        case (circular )
            outvalue = circular_topwidth_from_depth_singular  (idx, indepth)
        case (filled_circular)
            outvalue = filled_circular_topwidth_from_depth_singular  (idx, indepth)
        case (rectangular_closed)
            outvalue = rectangular_closed_topwidth_from_depth_singular  (idx, indepth)
        case (horiz_ellipse)
            outvalue = horiz_ellipse_topwidth_from_depth_singular (idx, indepth)
        case (vert_ellipse)
            outvalue = vert_ellipse_topwidth_from_depth_singular (idx, indepth)
        case (arch)
            outvalue = arch_topwidth_from_depth_singular (idx, indepth)
        case (eggshaped)
            outvalue = egg_shaped_topwidth_from_depth_singular (idx, indepth)
        case (horseshoe)
            outvalue = horse_shoe_topwidth_from_depth_singular (idx, indepth)
        case (gothic)
            outvalue = gothic_topwidth_from_depth_singular (idx, indepth)
        case (catenary)
            outvalue = catenary_topwidth_from_depth_singular (idx, indepth)
        case (semi_elliptical)
            outvalue = semi_elliptical_topwidth_from_depth_singular (idx, indepth)
        case (basket_handle)
            outvalue = basket_handle_topwidth_from_depth_singular (idx, indepth)
        case (semi_circular)
            outvalue = semi_circular_topwidth_from_depth_singular (idx, indepth)
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
        !!if (crashYN) return
        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% cycle through different geometries

        !% --- RECTANGULAR
        Npack => npack_elemPGx(epg_CC_rectangular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_nonsurcharged)
            call rectangular_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- RECTANGULAR CLOSED
        Npack => npack_elemPGx(epg_CC_rectangular_closed_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_closed_nonsurcharged)
            call rectangular_closed_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- TRAPEZOIDAL
        Npack => npack_elemPGx(epg_CC_trapezoidal_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_trapezoidal_nonsurcharged)
            call trapezoidal_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- TRIANGULAR
        Npack => npack_elemPGx(epg_CC_triangular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_triangular_nonsurcharged)
            call triangular_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- CIRCULAR
        Npack => npack_elemPGx(epg_CC_circular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_circular_nonsurcharged)
            call circular_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- SEMI-CIRCULAR
        Npack => npack_elemPGx(epg_CC_semi_circular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_semi_circular_nonsurcharged)
            call semi_circular_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- FILLED CIRCULAR
        Npack => npack_elemPGx(epg_CC_filled_circular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_filled_circular_nonsurcharged)
            call filled_circular_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- PARABOLIC
        Npack => npack_elemPGx(epg_CC_parabolic_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_parabolic_nonsurcharged)
            call parabolic_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- RECT_TRIANG
        Npack => npack_elemPGx(epg_CC_rectangular_triangular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_triangular_nonsurcharged)
            call rectangular_triangular_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- BASKET_HANDLE
        Npack => npack_elemPGx(epg_CC_basket_handle_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_basket_handle_nonsurcharged)
            call basket_handle_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- ARCH
        Npack => npack_elemPGx(epg_CC_arch_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_arch_nonsurcharged)
            call arch_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- HORIZONTAL ELLIPSE
        Npack => npack_elemPGx(epg_CC_horiz_ellipse_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_horiz_ellipse_nonsurcharged)
            call horiz_ellipse_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- VERTICAL ELLIPSE
        Npack => npack_elemPGx(epg_CC_vert_ellipse_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_vert_ellipse_nonsurcharged)
            call vert_ellipse_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- EGG_SHAPED
        Npack => npack_elemPGx(epg_CC_egg_shaped_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_egg_shaped_nonsurcharged)
            call egg_shaped_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- HORSE_SHOE
        Npack => npack_elemPGx(epg_CC_horse_shoe_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_horse_shoe_nonsurcharged)
            call horse_shoe_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- CATENARY
        Npack => npack_elemPGx(epg_CC_catenary_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_catenary_nonsurcharged)
            call catenary_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- GOTHIC
        Npack => npack_elemPGx(epg_CC_gothic_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_gothic_nonsurcharged)
            call gothic_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- MODIFIED BASKET HANDLE
        Npack => npack_elemPGx(epg_CC_mod_basket_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_mod_basket_nonsurcharged)
            call mod_basket_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- RECTANGULAR ROUND
        Npack => npack_elemPGx(epg_CC_rectangular_round_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_round_nonsurcharged)
            call rect_round_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- SEMI-ELLIPTICAL
        Npack => npack_elemPGx(epg_CC_semi_elliptical_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_semi_elliptical_nonsurcharged)
            call semi_elliptical_perimeter_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- IRREGULAR
        !%     note this requires first using the table lookup for hydraulic radius
        Npack => npack_elemPGx(epg_CC_irregular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_irregular_nonsurcharged)
            call irregular_hydradius_from_depth (elemPGx, Npack, thisCol)
            call irregular_perimeter_from_hydradius_area (elemPGx, Npack, thisCol)
        end if

        !% HACK NEED OTHER GEOMETRIES
        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_perimeter_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function geo_perimeter_from_depth_singular &
        (idx, indepth) result (outvalue)
        !%------------------------------------------------------------------
        !% Descriptions:
        !% computes the perimeter for a given depth of a single element
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(in)  :: indepth
            integer, intent(in)  :: idx
            character(64) :: subroutine_name = 'geo_perimeter_from_depth_singular'
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        select case (elemI(idx,ei_geometryType))
            
        case (rectangular)
            outvalue = rectangular_perimeter_from_depth_singular (idx, indepth)
        case (trapezoidal)
            outvalue = trapezoidal_perimeter_from_depth_singular (idx, indepth)
        case (triangular)
            outvalue = triangular_perimeter_from_depth_singular (idx, indepth)
        case (parabolic)
            outvalue = parabolic_perimeter_from_depth_singular (idx, indepth)
        case (power_function)
            print *, 'CODE ERROR: perimeter for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(338234)
        case (rect_triang)
            outvalue = rectangular_triangular_perimeter_from_depth_singular (idx, indepth)
        case (rect_round )
            outvalue = rect_round_perimeter_from_depth_singular (idx, indepth)
        case (mod_basket)
            outvalue = mod_basket_perimeter_from_depth_singular (idx, indepth)
        case (irregular)
            outvalue = irregular_geometry_from_depth_singular (idx,tt_area, indepth, setting%ZeroValue%Depth)
        case (circular )
            outvalue = circular_perimeter_from_depth_singular (idx, indepth)
        case (filled_circular)
            outvalue = filled_circular_perimeter_from_depth_singular (idx, indepth)
        case (rectangular_closed)
            outvalue = rectangular_closed_perimeter_from_depth_singular (idx, indepth)
        case (horiz_ellipse)
            outvalue = horiz_ellipse_perimeter_from_depth_singular (idx, indepth)
        case (vert_ellipse)
            outvalue = vert_ellipse_perimeter_from_depth_singular (idx, indepth)
        case (arch)
            outvalue = arch_perimeter_from_depth_singular (idx, indepth)
        case (eggshaped)
            outvalue = egg_shaped_perimeter_from_depth_singular (idx, indepth)
        case (horseshoe)
            outvalue = horse_shoe_perimeter_from_depth_singular (idx, indepth)
        case (gothic)
            outvalue = gothic_perimeter_from_depth_singular (idx, indepth)
        case (catenary)
            outvalue = catenary_perimeter_from_depth_singular (idx, indepth)
        case (semi_elliptical)
            outvalue = semi_elliptical_perimeter_from_depth_singular (idx, indepth)
        case (basket_handle)
            outvalue = basket_handle_perimeter_from_depth_singular (idx, indepth)
        case (semi_circular)
            outvalue = semi_circular_perimeter_from_depth_singular (idx, indepth)
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
    subroutine geo_hyddepth_from_depth (elemPGx, npack_elemPGx, col_elemPGx)
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
        !!if (crashYN) return
        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% cycle through different geometries

        !% --- RECTANGULAR
        Npack => npack_elemPGx(epg_CC_rectangular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_nonsurcharged)
            call rectangular_hyddepth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- RECTANGULAR CLOSED
        Npack => npack_elemPGx(epg_CC_rectangular_closed_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_closed_nonsurcharged)
            call rectangular_closed_hyddepth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- TRAPEZOIDAL
        Npack => npack_elemPGx(epg_CC_trapezoidal_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_trapezoidal_nonsurcharged)
            call trapezoidal_hyddepth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- TRIANGULAR
        Npack => npack_elemPGx(epg_CC_triangular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_triangular_nonsurcharged)
            call triangular_hyddepth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% --- CIRCULAR
        Npack => npack_elemPGx(epg_CC_circular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_circular_nonsurcharged)
            call circular_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        end if

        !% -- SEMI-CIRCULAR
        Npack => npack_elemPGx(epg_CC_semi_circular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_semi_circular_nonsurcharged)
            call semi_circular_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        end if

        !% --- FILLED CIRCULAR
        Npack => npack_elemPGx(epg_CC_filled_circular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_filled_circular_nonsurcharged)
            call filled_circular_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        end if

        !% --- PARABOLIC
        Npack => npack_elemPGx(epg_CC_parabolic_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_parabolic_nonsurcharged)
            call parabolic_hyddepth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- RECT_TRIANG
        Npack => npack_elemPGx(epg_CC_rectangular_triangular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_triangular_nonsurcharged)
            call rectangular_triangular_hyddepth_from_depth (elemPGx, Npack, thisCol)
        end if

        !% -- BASKET_HANDLE
        Npack => npack_elemPGx(epg_CC_basket_handle_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_basket_handle_nonsurcharged)
            call basket_handle_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        end if

        !% -- ARCH
        Npack => npack_elemPGx(epg_CC_arch_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_arch_nonsurcharged)
            call arch_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        end if

        !% -- HORIZONTAL ELLIPSE
        Npack => npack_elemPGx(epg_CC_horiz_ellipse_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_horiz_ellipse_nonsurcharged)
            call horiz_ellipse_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        end if

        !% -- VERTICAL ELLIPSE
        Npack => npack_elemPGx(epg_CC_vert_ellipse_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_vert_ellipse_nonsurcharged)
            call vert_ellipse_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        end if

        !% -- EGG_SHAPED
        Npack => npack_elemPGx(epg_CC_egg_shaped_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_egg_shaped_nonsurcharged)
            call egg_shaped_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        end if

        !% -- HORSE_SHOE
        Npack => npack_elemPGx(epg_CC_horse_shoe_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_horse_shoe_nonsurcharged)
            call horse_shoe_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        end if

        !% -- CATENARY
        Npack => npack_elemPGx(epg_CC_catenary_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_catenary_nonsurcharged)
            call catenary_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        end if

        !% -- GOTHIC
        Npack => npack_elemPGx(epg_CC_gothic_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_gothic_nonsurcharged)
            call gothic_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        end if

        !% -- MODIFIED BASKET HANDLE
        Npack => npack_elemPGx(epg_CC_mod_basket_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_mod_basket_nonsurcharged)
            call mod_basket_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        end if

        !% -- RECTANGULAR ROUND
        Npack => npack_elemPGx(epg_CC_rectangular_round_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_rectangular_round_nonsurcharged)
            call rect_round_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        end if

        !% -- SEMI-ELLIPTICAL
        Npack => npack_elemPGx(epg_CC_semi_elliptical_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_semi_elliptical_nonsurcharged)
            call semi_elliptical_hyddepth_from_topwidth (elemPGx, Npack, thisCol)
        end if

        !% --- IRREGULAR
        Npack => npack_elemPGx(epg_CC_irregular_nonsurcharged)
        if (Npack > 0) then
            thisCol => col_elemPGx(epg_CC_irregular_nonsurcharged)
            call irregular_hyddepth_from_topwidth_area (elemPGx, Npack, thisCol)
        end if

        !% HACK need other geometries
        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_hyddepth_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function geo_hyddepth_from_depth_singular &
        (idx, indepth) result (outvalue)
        !%------------------------------------------------------------------
        !% Descriptions:
        !% computes the hydraulic depth for a given depth of a single element
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(in)  :: indepth
            integer, intent(in)  :: idx
            real(8)              :: temp1, temp2
            character(64) :: subroutine_name = 'geo_hyddepth_from_depth_singular'
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        select case (elemI(idx,ei_geometryType))
            
        case (rectangular)
            outvalue = indepth
        case (trapezoidal)
            outvalue = trapezoidal_hyddepth_from_depth_singular (idx, indepth)
        case (triangular)
            outvalue = triangular_hyddepth_from_depth_singular (idx, indepth)
        case (parabolic)
            outvalue = parabolic_hyddepth_from_depth_singular (idx, indepth)
        case (power_function)
            print *, 'CODE ERROR: hyddepth for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(449734)
        case (rect_triang)
            outvalue = rectangular_triangular_hyddepth_from_depth_singular (idx, indepth)
        case (rect_round )
            temp1    = rect_round_topwidth_from_depth_singular    (idx, indepth)
            outvalue = rect_round_hyddepth_from_topwidth_singular (idx,temp1,indepth)
        case (mod_basket)
            temp1    = mod_basket_topwidth_from_depth_singular    (idx, indepth)
            outvalue = mod_basket_hyddepth_from_topwidth_singular (idx,temp1,indepth)
        case (irregular)
            !% --- get the area and topwidth, then compute the hydraulic depth
            temp1 = irregular_geometry_from_depth_singular (idx,tt_area,  indepth, setting%ZeroValue%Area)
            temp2 = irregular_geometry_from_depth_singular (idx,tt_width, indepth, setting%ZeroValue%TopWidth)
            outvalue = temp1 / temp2
        case (circular )
            !% --- get the topwidth and use that to compute the hydraulic depth
            temp1    = circular_topwidth_from_depth_singular    (idx, indepth)
            outvalue = circular_hyddepth_from_topwidth_singular (idx,temp1,indepth)
        case (filled_circular)
            !% --- get the topwidth and use that to compute the hydraulic depth
            temp1    = filled_circular_topwidth_from_depth_singular    (idx, indepth)
            outvalue = filled_circular_hyddepth_from_topwidth_singular (idx,temp1,indepth)
        case (rectangular_closed)
            outvalue = rectangular_closed_hyddepth_from_depth_singular (idx, indepth)
        case (horiz_ellipse)
            !% --- get the topwidth and use that to compute the hydraulic depth
            temp1    = horiz_ellipse_topwidth_from_depth_singular    (idx, indepth)
            outvalue = horiz_ellipse_hyddepth_from_topwidth_singular (idx,temp1,indepth)
        case (vert_ellipse)
            !% --- get the topwidth and use that to compute the hydraulic depth
            temp1    = vert_ellipse_topwidth_from_depth_singular    (idx, indepth)
            outvalue = vert_ellipse_hyddepth_from_topwidth_singular (idx,temp1,indepth)
        case (arch)
            !% --- get the topwidth and use that to compute the hydraulic depth
            temp1    = arch_topwidth_from_depth_singular    (idx, indepth)
            outvalue = arch_hyddepth_from_topwidth_singular (idx,temp1,indepth)
        case (eggshaped)
            temp1    = egg_shaped_topwidth_from_depth_singular    (idx, indepth)
            outvalue = egg_shaped_hyddepth_from_topwidth_singular (idx,temp1,indepth)
        case (horseshoe)
            !% --- get the topwidth and use that to compute the hydraulic depth
            temp1    = horse_shoe_topwidth_from_depth_singular    (idx, indepth)
            outvalue = horse_shoe_hyddepth_from_topwidth_singular (idx,temp1,indepth)
        case (gothic)
            temp1    = gothic_topwidth_from_depth_singular    (idx, indepth)
            outvalue = gothic_hyddepth_from_topwidth_singular (idx,temp1,indepth)
        case (catenary)
            temp1    = catenary_topwidth_from_depth_singular    (idx, indepth)
            outvalue = catenary_hyddepth_from_topwidth_singular (idx,temp1,indepth)
        case (semi_elliptical)
            temp1    = semi_elliptical_topwidth_from_depth_singular    (idx, indepth)
            outvalue = semi_elliptical_hyddepth_from_topwidth_singular (idx,temp1,indepth)
        case (basket_handle)
            !% --- get the topwidth and use that to compute the hydraulic depth
            temp1    = basket_handle_topwidth_from_depth_singular    (idx, indepth)
            outvalue = basket_handle_hyddepth_from_topwidth_singular (idx,temp1,indepth)
        case (semi_circular)
            temp1    = semi_circular_topwidth_from_depth_singular    (idx, indepth)
            outvalue = semi_circular_hyddepth_from_topwidth_singular (idx,temp1,indepth)
        case (custom)
            print *, 'CODE ERROR: hyddepth for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(449734)
        case (force_main)
            print *, 'CODE ERROR: hyddepth for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'in ',trim(subroutine_name)   
            print *, 'This should never be reached as a force_main is not a valid geometryType'  
            call util_crashpoint(449734)
        case default
            print *, 'CODE ERROR: hyddepth for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(449734)
        end select
           
    end function geo_hyddepth_from_depth_singular
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
        !if (crashYN) return
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
    subroutine geo_ell_from_head (thisColP)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes the value of "ell" -- the modified hydraulic depth
        !% used as a length scale in AC method
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisColP
        integer, pointer :: thisP(:), Npack
        real(8), pointer :: ell(:), head(:), area(:), topwidth(:), hydDepth(:)
        real(8), pointer :: ZbreadthMax(:), breadthMax(:), areaBelowBreadthMax(:)
        integer :: ii

        character(64) :: subroutine_name = 'geo_ell'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
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
    end subroutine geo_ell_from_head
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function geo_ell_singular (indx) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes the value of "ell" the modified hydraulic depth
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), pointer :: head(:), area(:), topwidth(:)
        real(8), pointer :: ZbreadthMax(:), breadthMax(:), areaBelowBreadthMax(:)

        character(64) :: subroutine_name = 'geo_ell_singular'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
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
   subroutine geo_head_from_ell (thisColP)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes the value of "ell" -- the modified hydraulic depth
        !% used as a length scale in AC method
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisColP
        integer, pointer :: thisP(:), Npack
        real(8), pointer :: ell(:), head(:), zbottom(:)
        character(64) :: subroutine_name = 'geo_head_from_ell'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        Npack   => npack_elemP(thisColP)
        ell     => elemR(:,er_ell)
        head    => elemR(:,er_Head)
        zbottom => elemR(:,er_Zbottom)
        !%-----------------------------------------------------------------------------

        if (Npack > 0) then
            thisP       => elemP(1:Npack,thisColP)        
            head(thisP) = zbottom(thisP) + ell(thisP)
        end if

        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_head_from_ell 
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
    subroutine geo_CC_slot_adjustments (thisColP_closed_CC)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% This subroutine adds back the slot geometry in all the closed elements
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisColP_closed_CC
        integer, pointer    :: thisP(:), Npack, SlotMethod
        real(8), pointer    :: SlotWidth(:), SlotVolume(:), SlotDepth(:), dSlotDepth(:)
        real(8), pointer    :: volume(:), ell(:), depth(:), area(:), SlotArea(:)
        real(8), pointer    :: head(:), fullVolume(:), fullArea(:), fullDepth(:)
        real(8), pointer    :: Overflow(:), zcrown(:), ellMax(:), SlotHydRad(:)
        logical, pointer    :: isSlot(:)

        character(64) :: subroutine_name = 'geo_CC_slot_adjustments'
        !%-----------------------------------------------------------------------------

        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        Npack      => npack_elemP(thisColP_closed_CC)
        area       => elemR(:,er_Area)
        depth      => elemR(:,er_Depth)
        dSlotDepth => elemR(:,er_dSlotDepth)
        ell        => elemR(:,er_ell)
        ellMax     => elemR(:,er_ell_max)
        fullDepth  => elemR(:,er_FullDepth)
        fullvolume => elemR(:,er_FullVolume)
        fullArea   => elemR(:,er_FullArea)
        head       => elemR(:,er_Head)
        Overflow   => elemR(:,er_VolumeOverFlow)
        SlotWidth  => elemR(:,er_SlotWidth)
        SlotVolume => elemR(:,er_SlotVolume)
        SlotDepth  => elemR(:,er_SlotDepth)
        SlotArea   => elemR(:,er_SlotArea)
        SlotHydRad => elemR(:,er_SlotHydRadius)
        volume     => elemR(:,er_Volume)
        zcrown     => elemR(:,er_Zcrown)
        isSlot     => elemYN(:,eYN_isSlot)

        !% pointer to necessary settings struct
        SlotMethod => setting%PreissmannSlot%PreissmannSlotMethod
        !%-----------------------------------------------------------------------------

        !% CC slot adjustment
        if (Npack > 0) then
            thisP    => elemP(1:Npack,thisColP_closed_CC)

            select case (SlotMethod)
                case (StaticSlot)
                    where (isSlot(thisP)) 
                        volume(thisP) = volume(thisP)  + SlotVolume(thisP)
                        ! area(thisP)   = area(thisP)    + SlotArea(thisP) 
                        depth(thisP)  = depth(thisP)   + SlotDepth(thisP)
                        head(thisP)   = head(thisP)    + SlotDepth(thisP)
                        ell(thisP)    = ellMax(thisP)
                        Overflow(thisP) = zeroR
                    end where 

                case (DynamicSlot)
                    where (isSlot(thisP)) 
                        volume(thisP)    = volume(thisP)  + SlotVolume(thisP)
                        ! area(thisP)      = area(thisP)    + SlotArea(thisP)
                        SlotDepth(thisP) = max(SlotDepth(thisP) + dSlotDepth(thisP), zeroR) !% HACK: TESTING STUFF
                        head(thisP)      = head(thisP)  + SlotDepth(thisP)
                        ell(thisP)       = ellMax(thisP)
                        Overflow(thisP)  = zeroR
                    elsewhere
                        SlotDepth(thisP)  = zeroR
                        ! SlotVolume(thisP) = zeroR
                    end where 
            end select
        end if

        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine geo_CC_slot_adjustments
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
            elemR(thisP,er_HydDepth) = elemR(thisP,er_Depth)
            elemR(thisP,er_ell)      = elemR(thisP,er_Depth)
        end if
        !%------------------------------------------------------------------
        !% Closing
    end subroutine geo_JM_values
    !%
!%==========================================================================
!%==========================================================================
!%
    subroutine geo_JB_slot_computation_ETM (thisColP_JM)
        !%------------------------------------------------------------------
        !% Description:
        !%      Slot computation for Junction Branches
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisColP_JM
            integer, pointer :: Npack, thisP(:), tM, BranchExists(:)
            real(8), pointer :: area(:), depth(:), head(:), length(:), volume(:), zcrown(:), zbottom(:)
            real(8), pointer :: fullDepth(:), fullArea(:), fPNumber(:), PNumber(:), PCelerity(:), ell(:)
            real(8), pointer :: SlotWidth(:), SlotVolume(:), SlotDepth(:), SlotArea(:), ellMax(:)
            real(8), pointer :: overflow(:), grav, TargetPCelerity, PreissmannAlpha
            logical, pointer :: isSlot(:) , fSlot(:), isDnJB(:)
            integer, pointer :: SlotMethod, fUp(:), fDn(:)
            integer :: tB, ii, kk
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases
            Npack         => npack_elemP(thisColP_JM)
            area          => elemR(:,er_Area)
            depth         => elemR(:,er_Depth)
            head          => elemR(:,er_Head)
            length        => elemR(:,er_Length)
            fullArea      => elemR(:,er_FullArea)
            fullDepth     => elemR(:,er_FullDepth)
            overflow      => elemR(:,er_VolumeOverFlow)
            volume        => elemR(:,er_Volume)
            zcrown        => elemR(:,er_Zcrown)
            zbottom       => elemR(:,er_Zbottom)
            ell           => elemR(:,er_ell)
            ellMax        => elemR(:,er_ell_max)
            fUp           => elemI(:,ei_Mface_uL)
            fDn           => elemI(:,ei_Mface_dL)
            BranchExists  => elemSI(:,esi_JunctionBranch_Exists)
            grav          => setting%Constant%gravity
        !% Slot Aliases
            PNumber    => elemR(:,er_Preissmann_Number)
            PCelerity  => elemR(:,er_Preissmann_Celerity)
            SlotWidth  => elemR(:,er_SlotWidth)
            SlotVolume => elemR(:,er_SlotVolume)
            SlotDepth  => elemR(:,er_SlotDepth)
            SlotArea   => elemR(:,er_SlotArea)
            fPNumber   => faceR(:,fr_Preissmann_Number)
            isSlot     => elemYN(:,eYN_isSlot)
            isDnJB     => elemYN(:,eYN_isDownstreamJB)
            fSlot      => faceYN(:,fYN_isSlot)
            SlotMethod      => setting%PreissmannSlot%PreissmannSlotMethod
            TargetPCelerity => setting%PreissmannSlot%TargetPreissmannCelerity
            PreissmannAlpha => setting%PreissmannSlot%PreissmannAlpha
        !%------------------------------------------------------------------

        !% JB slot adjustment
        if (Npack > 0) then
            thisP  => elemP(1:Npack,thisColP_JM)
            !% cycle through the all the main junctions and each of its branches
            do ii=1,Npack
                tM => thisP(ii) !% junction main ID
                ! handle the upstream branches
                do kk=1,max_branch_per_node,2
                    tB = tM + kk  !% JB branch ID
                    if (BranchExists(tB)==1) then
                        !% initialize slot
                        isSlot(tB)     = .false.
                        SlotDepth(tB)  = zeroR
                        SlotArea(tB)   = zeroR
                        SlotWidth(tB)  = zeroR
                        SlotVolume(tB) = zeroR
                        PCelerity(tB)  = zeroR

                        !% assuming a slot if the head is above the crown
                        !% or the upstream CC is in a slot
                        if (head(tB) .gt. zcrown(tB)) then
                            isSlot(tB)     = .true.
                            fSlot(fUp(tB)) = .true.
                            PNumber(tB)    = fPNumber(fUp(tB))
                            PCelerity(tB)  = min(TargetPCelerity / PNumber(tB), TargetPCelerity)
                            SlotDepth(tB)  = max(depth(tB) - fulldepth(tB), zeroR)   
                            SlotArea(tB)   = (SlotDepth(tB) * (PNumber(tB)**twoR) * grav * &
                                                fullArea(tB)) / (TargetPCelerity ** twoR)
                            SlotVolume(tB) = SlotArea(tB) * length(tB)
                            
                            !% add the slot geometry back to previously solved geometry
                            volume(tB) = volume(tB)  + SlotVolume(tB)
                            area(tB)   = area(tB)    + SlotArea(tB)
                            depth(tB)  = depth(tB)   + SlotDepth(tB)
                            Overflow(tB) = zeroR
                        end if  
                    end if
                end do
                !% handle the downstream branches
                do kk=2,max_branch_per_node,2
                    tB = tM + kk
                    if (BranchExists(tB)==1) then
                        !% initialize slot
                        isSlot(tB)     = .false.
                        SlotDepth(tB)  = zeroR
                        SlotArea(tB)   = zeroR
                        SlotWidth(tB)  = zeroR
                        SlotVolume(tB) = zeroR
                        PCelerity(tB)  = zeroR

                        !% assuming a slot if the head is above the crown
                        !% or the downstream CC is in a slot
                        if (head(tB) .gt. zcrown(tB)) then
                            isSlot(tB)     = .true.
                            fSlot(fDn(tB)) = .true.
                            PNumber(tB)    = fPNumber(fDn(tB))
                            PCelerity(tB)  = min(TargetPCelerity / PNumber(tB), TargetPCelerity)
                            SlotDepth(tB)  = max(depth(tB) - fulldepth(tB), zeroR)    
                            SlotArea(tB)   = (SlotDepth(tB) * (PNumber(tB)**twoR) * grav * &
                                                fullArea(tB)) / (TargetPCelerity ** twoR)
                            SlotVolume(tB) = SlotArea(tB) * length(tB)

                            !% add the slot geometry back to previously solved geometry
                            volume(tB) = volume(tB)  + SlotVolume(tB)
                            area(tB)   = area(tB)    + SlotArea(tB)
                            depth(tB)  = depth(tB)   + SlotDepth(tB)
                            Overflow(tB) = zeroR
                        end if
                    end if
                end do
            end do
        end if
                  
        !%------------------------------------------------------------------
        !% Closing
    end subroutine geo_JB_slot_computation_ETM
!%
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module geometry