module geometry

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use rectangular_channel
    use rectangular_conduit
    use trapezoidal_channel
    use triangular_channel
    use circular_conduit
    use irregular_channel
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
    public :: geo_assign_JB  !BRHbugfix 20210813
    public :: geo_topwidth_from_depth
    public :: geo_hyddepth_from_depth_singular
    public :: geo_topwidth_from_depth_singular
    public :: geo_area_from_depth_singular
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
            integer, pointer :: thisColP_JM, thisColP_JB, thisColP_ClosedElems
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
                    call util_crashpoint(7389)
                    !return
                    !stop 7389
            end select
            call util_crashstop(49872)
        !%-----------------------------------------------------------------------------
        !% STATUS: at this point we know volume on Non-surcharged CC, JM,
        !% elements and head on all surcharged CC, JM elements

            ! call util_CLprint ('in geometry at top')    

        !% --- assign all geometry for surcharged elements CC, JM and JB
        !%     Note: not used in Preissmann Slot
        call geo_surcharged (thisColP_surcharged)

            ! call util_CLprint ('in geometry before adjust_limit_by_zerovalues') 

        !% --- reset all zero or near-zero volumes in non-surcharged CC, JM, and JB
        call adjust_limit_by_zerovalues (er_Volume, setting%ZeroValue%Volume, thisColP_NonSurcharged, .true.)

            ! print *, this_image(),  '    geomTL ccc',this_image(),setting%Time%Step
            ! call util_CLprint ('in geometry before geo_depth_from_volume') 

        !% --- compute the depth on all non-surcharged elements of CC, JM and JB
        call geo_depth_from_volume (elemPGx, npack_elemPGx, col_elemPGx)

            ! print *, this_image(),  '    geomTL  ddd',this_image(),setting%Time%Step
            ! call util_CLprint ('in geometry before adjust_limit_by_zerovalues (2)') 

        !% reset all zero or near-zero depths in non-surcharged CC and JM and JB
        call adjust_limit_by_zerovalues (er_Depth, setting%ZeroValue%Depth, thisColP_NonSurcharged, .false.)

            !print *,this_image(),  '     geomTL  eee',this_image(),setting%Time%Step
            ! call util_CLprint ('in geometry before geo_head_from_depth') 

        !% --- compute the head on all non-surcharged elements of CC and JM and JB
        !%     This sets head consistent with depth
        call geo_head_from_depth (thisColP_NonSurcharged)

            !print *, this_image(),  '    geomTL  fff',this_image(),setting%Time%Step
            ! call util_CLprint ('in geometry before geo_limit_incipient_surcharge (Volume)') 

        !% --- limit volume for incipient surcharge. This is done after depth is computed
        !%     so that the "depth" algorithm can include depths greater than fulldepth
        !%     as a way to handle head for incipient surcharge.
        call geo_limit_incipient_surcharge (er_Volume, er_FullVolume, thisColP_NonSurcharged,.true.) !% 20220124brh

            !print *, this_image(),  '    geomTL  ggg',this_image(),setting%Time%Step
            ! call util_CLprint ('in geometry before geo_limit_incipient_surcharge (Depth)')  

        ! print *, 'in ',trim(subroutine_name),elemR(48,er_VolumeOverFlow)

        !% limit depth for incipient surcharged. This is done after head is computed
        !% so that the depth algorithm can include depths greater than fulldepth to
        !% handle incipient surcharge
        !call geo_limit_incipient_surcharge (er_Depth, er_FullDepth, thisColP_NonSurcharged)
        call geo_limit_incipient_surcharge (er_Depth, er_FullDepth, thisColP_NonSurcharged,.false.) !% 20220124brh

            !print *, this_image(),  '    geomTL  hhh',setting%Time%Step
            !  call util_CLprint ('in geometry before geo_assign_JB') 

        !% STATUS: at this point we know depths and heads in all CC, JM elements
        !% (surcharged and nonsurcharged) with limiters for conduit depth and zero depth
           
        !% assign the head, depth, geometry on junction branches JB based on JM head
        call geo_assign_JB (whichTM, thisColP_JM)

            !print *, this_image(),  '    geomTL  iii',setting%Time%Step
            !  call util_CLprint ('in geometry before geo_area_from_volume')  

        !% STATUS at this point we know geometry on all JB and all surcharged, with
        !% depth, head, volume on all non-surcharged or incipient surcharge.

        !% compute area from volume for CC, JM nonsurcharged
        call geo_area_from_volume (thisColP_NonSurcharged)

            ! print *, this_image(),  '    geomTL  jjj',this_image()
            ! call util_CLprint ('in geometry before adjust_limit_by_zerovalues') 

        !% reset all zero or near-zero areas in non-surcharged CC and JM
        call adjust_limit_by_zerovalues (er_Area, setting%ZeroValue%Area, thisColP_NonSurcharged, .false.)

            ! print *, this_image(),  '    geomTL kkk',this_image()
            ! call util_CLprint ('in geometry before topwidth_from_depth')   

        !% compute topwidth from depth for all CC, JM nonsurcharged
        call geo_topwidth_from_depth (elemPGx, npack_elemPGx, col_elemPGx)

            ! print *, this_image(),  '    geomTL  lll', this_image()
            ! call util_CLprint ('in geometry before adjust_limit_by_zerovalues') 

        !% reset all zero or near-zero topwidth in non-surcharged CC and JM
        !% but do not change the eYN(:,eYN_isZeroDepth) mask
        call adjust_limit_by_zerovalues (er_Topwidth, setting%ZeroValue%Topwidth, thisColP_NonSurcharged, .false.)

            ! print *, this_image(),  '    geomTL  mmm',this_image()
            ! call util_CLprint ('in geometry before perimeter_from_depth') 

        !% compute perimeter from maximum depth for all CC, JM nonsurcharged
        call geo_perimeter_from_depth (elemPGx, npack_elemPGx, col_elemPGx)

            ! print *, this_image(),  '    geomTL  nnn',this_image()
            ! call util_CLprint ('in geometry before hyddepth_from_depth') 

        !% compute hyddepth
        call geo_hyddepth_from_depth (elemPGx, npack_elemPGx, col_elemPGx)

            ! print *, this_image(),  '    geomTL  ooo',this_image()
            ! call util_CLprint ('in geometry before hydradius_from_area_perimeter')   

        !% compute hydradius  (applies to all nonsurcharged)
        call geo_hydradius_from_area_perimeter (thisColP_NonSurcharged)


            ! print *, this_image(),  '    geomTL  qqq',this_image()
            ! call util_CLprint ('in geometry before ell_from_head') 

        !% the modified hydraulic depth "ell" is used for AC computations and
        !% for Froude number computations on all elements, whether ETM or AC.
        call geo_ell_from_head (thisColP_all)

            ! print *,  this_image(),  '    geomTL  rrr',this_image()
            ! call util_CLprint ('in geometry before slot_adjustments') 

        !% make adjustments for slots on closed elements only for ETM
        if (whichTM .eq. ETM) then
            call geo_slot_adjustments (thisColP_ClosedElems)
        end if

            ! print *,  this_image(),  '    geomTL  sss',this_image()
            ! call util_CLprint ('in geometry before JM_values') 

        !% Set JM values that are not otherwise defined
        call geo_JM_values ()

            ! print *, this_image(),  '    geomTL ttt',this_image()
            ! call util_CLprint ('in geometry before dHdA') 

        !% compute the dHdA that are only for AC nonsurcharged
        if (whichTM .ne. ETM) then
            call geo_dHdA (ep_NonSurcharged_AC)
        end if

            ! print *,  this_image(),  '    geomTL uuu',this_image()
            ! call util_CLprint ('in geometry at end') 

        call util_crashstop(322983)

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
            real(8), pointer :: zCrown(:), fullArea(:), fulldepth(:), fullperimeter(:)
            real(8), pointer :: fullhyddepth(:), thisTable(:,:)
            real(8), pointer :: grav       

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

                                !write(*,"(A,i5,10f12.5)") 'AAA ell ',tB, ell(tB), depth(tB), hydDepth(tB), fulldepth(tB)

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

                                !write(*,"(A,i5,10f12.5)"), 'BBB ell ',tB, ell(tB), depth(tB), hydDepth(tB), fulldepth(tB)

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

                                !write(*,"(A,i5,10f12.5)") 'CCC ell ',tB, ell(tB), depth(tB), hydDepth(tB), fulldepth(tB)

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

                                    ! print *, 'in geo_assign_JB  for rect element'
                                    ! print *, 'area ',area(tB), depth(tB)

                                   ! write(*,"(A,i5,10f12.5)") 'DDD ell ',tB, ell(tB), depth(tB), hydDepth(tB), fulldepth(tB)

                                case (triangular)
                                    area(tB)     = triangular_area_from_depth_singular      (tB,depth(tB))
                                    topwidth(tB) = triangular_topwidth_from_depth_singular  (tB,depth(tB))
                                    hydDepth(tB) = triangular_hyddepth_from_depth_singular  (tB,depth(tB))
                                    perimeter(tB)= triangular_perimeter_from_depth_singular (tB,depth(tB))
                                    hydRadius(tB)= triangular_hydradius_from_depth_singular (tB,depth(tB))
                                    ell(tB)      = geo_ell_singular (tB) 
                                    dHdA(tB)     = oneR / topwidth(tB)

                                   ! write(*,"(A,i5,10f12.5)") 'EEE ell ',tB, ell(tB), depth(tB), hydDepth(tB), fulldepth(tB)
                                    
                                case (trapezoidal)                                    
                                    area(tB)     = trapezoidal_area_from_depth_singular      (tB,depth(tB))
                                    topwidth(tB) = trapezoidal_topwidth_from_depth_singular  (tB,depth(tB))
                                    hydDepth(tB) = trapezoidal_hyddepth_from_depth_singular  (tB,depth(tB))
                                    perimeter(tB)= trapezoidal_perimeter_from_depth_singular (tB,depth(tB))
                                    hydRadius(tB)= trapezoidal_hydradius_from_depth_singular (tB,depth(tB))
                                    ell(tB)      = geo_ell_singular (tB) 
                                    dHdA(tB)     = oneR / topwidth(tB)

                                   ! write(*,"(A,i5,10f12.5)") 'FFF ell ',tB, ell(tB), depth(tB), hydDepth(tB), fulldepth(tB)

                                case (circular)
                                    area(tB)     = circular_area_from_depth_singular          (tB,depth(tB))
                                    topwidth(tB) = circular_topwidth_from_depth_singular      (tB,depth(tB))
                                    hydDepth(tB) = circular_hyddepth_from_topwidth_singular   (tB,topwidth(tB),depth(tB))
                                    hydRadius(tB)= circular_hydradius_from_depth_singular     (tB,depth(tB))
                                    perimeter(tB)= circular_perimeter_from_hydradius_singular (tB,hydRadius(tB))
                                    ell(tB)      = geo_ell_singular (tB) 
                                    dHdA(tB)     = oneR / topwidth(tB)

                                    !write(*,"(A,i5,10f12.5)"), 'GGG ell ',tB, ell(tB), depth(tB), hydDepth(tB), fulldepth(tB)

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

                                   ! write(*,"(A,i5,10f12.5)") 'HHH ell ',tB, ell(tB), depth(tB), hydDepth(tB), fulldepth(tB)

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

                            ! print *, 'in geo_assign_JB at bottom'
                            ! write(*,"(A,i5,10f12.5)") 'III ell ',tB, ell(tB), depth(tB), hydDepth(tB), fulldepth(tB)
                            !write(*,"(A,10f12.5)") 'hyd depth', hydDepth(tB)
                            !print *, area(tB), length(tB)

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

        !call util_CLprint('start of geo depth from volume')        

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

        !call util_CLprint('after rectangular') 

        !% --- TRAPEZOIDAL CC
        thisCol => col_elemPGx(epg_CC_trapezoidal_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call trapezoidal_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !call util_CLprint('after trapezoidal') 

        !% --- TRIANGULAR CC
        thisCol => col_elemPGx(epg_CC_triangular_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call triangular_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !call util_CLprint('after triangular') 

        !% --- CIRCULAR CC
        thisCol => col_elemPGx(epg_CC_circular_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call circular_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !call util_CLprint('after circular') 
 
        !% --- IRREGULAR CC
        thisCol => col_elemPGx(epg_CC_irregular_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call irregular_depth_from_volume (elemPGx, Npack, thisCol)
        end if

        !call util_CLprint('after irregular') 
        !% HACK Needs additional geometries

        !% JM with functional geometry
        thisCol => col_elemPGx(epg_JM_functionalStorage_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call storage_functional_depth_from_volume (elemPGx, Npack, thisCol)
            !call storage_implied_length(elemPGx, Npack, thisCol)
        end if

        !call util_CLprint('after functional storage') 

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
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(33234)
        case (power_function)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(33234)
        case (rect_triang)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(33234)
        case (rect_round )
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(33234)
        case (mod_basket)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(33234)
        case (irregular)
            outvalue = irregular_geometry_from_depth_singular (idx,tt_area, indepth, setting%ZeroValue%Depth)
        case (circular )
            outvalue = circular_area_from_depth_singular (idx, indepth)
        case (filled_circular)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(33234)
        case (rectangular_closed)
            outvalue = rectangular_closed_area_from_depth_singular (idx, indepth)
        case (horiz_ellipse)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(33234)
        case (vert_ellipse)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(33234)
        case (arch)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(33234)
        case (eggshaped)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(33234)
        case (horseshoe)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(33234)
        case (gothic)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(33234)
        case (catenary)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(33234)
        case (semi_elliptical)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(33234)
        case (basket_handle)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(33234)
        case (semi_circular)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(33234)
        case (custom)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(33234)
        case (force_main)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)   
            call util_crashpoint(33234)
        case default
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
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
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (power_function)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (rect_triang)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (rect_round )
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (mod_basket)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (irregular)
            outvalue = irregular_geometry_from_depth_singular (idx,tt_width, indepth, setting%ZeroValue%TopWidth)
        case (circular )
            outvalue = circular_topwidth_from_depth_singular  (idx, indepth)
        case (filled_circular)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (rectangular_closed)
            outvalue = rectangular_closed_topwidth_from_depth_singular  (idx, indepth)
        case (horiz_ellipse)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (vert_ellipse)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (arch)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (eggshaped)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (horseshoe)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (gothic)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (catenary)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (semi_elliptical)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (basket_handle)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (semi_circular)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (custom)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (force_main)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)   
            call util_crashpoint(4498734)
        case default
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
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
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (power_function)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (rect_triang)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (rect_round )
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (mod_basket)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
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
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (rectangular_closed)
            outvalue = rectangular_closed_hyddepth_from_depth_singular (idx, indepth)
        case (horiz_ellipse)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (vert_ellipse)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (arch)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (eggshaped)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (horseshoe)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (gothic)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (catenary)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (semi_elliptical)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (basket_handle)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (semi_circular)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (custom)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
        case (force_main)
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)   
            call util_crashpoint(4498734)
        case default
            print *, 'CODE ERROR: geometry code for cross-section ',trim(reverseKey(elemI(idx,ei_geometryType)))
            print *, 'has not been implemented in ',trim(subroutine_name)
            call util_crashpoint(4498734)
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


        ! print *, 'in geo ell_from_head'
        ! do ii=1,size(thisP)
        !     write(*,"(i5,10f12.5)")  thisP(ii), head(thisP(ii)), elemR(thisP(ii),er_Zbottom), elemR(thisP(ii),er_Depth), &
        !         ZbreadthMax(thisP(ii)), breadthMax(thisP(ii)), areaBelowBreadthMax(thisP(ii)),  ell(thisP(ii)), hydDepth(thisP(ii))
        ! end do
        ! print *, 'thisP ',thisP
        ! print *, 'head  ',head(thisP)
        ! print *, 'Zbmax ',ZbreadthMax(thisP)
        ! print *, 'ell   ',ell(thisP)

        ! print *, 'at end of geo_ell_from_head', elemR(15,er_ell)

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
    subroutine geo_slot_adjustments (thisColP)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% This subroutine adds back the slot geometry in all the closed elements
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisColP
        integer, pointer    :: thisP(:), Npack
        real(8), pointer    :: SlotWidth(:), SlotVolume(:), SlotDepth(:), SlotArea(:)
        real(8), pointer    :: volume(:), volumeN0(:), depth(:), area(:)
        real(8), pointer    :: head(:), headN0(:), fullVolume(:), fullArea(:), fullDepth(:)
        real(8), pointer    :: Overflow(:), zbottom(:), ellMax(:)

        character(64) :: subroutine_name = 'geo_slot_adjustments'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        Npack      => npack_elemP(thisColP)
        area       => elemR(:,er_Area)
        volume     => elemR(:,er_Volume)
        volumeN0   => elemR(:,er_Volume_N0)
        Overflow   => elemR(:,er_VolumeOverFlow)
        depth      => elemR(:,er_Depth)
        ellMax     => elemR(:,er_ell_max)
        fullDepth  => elemR(:,er_FullDepth)
        fullvolume => elemR(:,er_FullVolume)
        fullArea   => elemR(:,er_FullArea)
        head       => elemR(:,er_Head)
        headN0     => elemR(:,er_Head_N0)
        zbottom    => elemR(:,er_Zbottom)
        SlotWidth  => elemR(:,er_SlotWidth)
        SlotVolume => elemR(:,er_SlotVolume)
        SlotDepth  => elemR(:,er_SlotDepth)
        SlotArea   => elemR(:,er_SlotArea)

        !%-----------------------------------------------------------------------------

        
        if (Npack > 0) then
            thisP    => elemP(1:Npack,thisColP)

            !print *, 'in geo_slot',this_image(), thisP

            where (SlotVolume(thisP) .gt. zeroR) 
                volume(thisP) = volume(thisP)  + SlotVolume(thisP)
                area(thisP)   = area(thisP)    + SlotArea(thisP)
                depth(thisP)  = depth(thisP)   + SlotDepth(thisP)
                head(thisP)   = zbottom(thisP) + fullDepth(thisP) + SlotDepth(thisP)
                Overflow(thisP) = zeroR
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