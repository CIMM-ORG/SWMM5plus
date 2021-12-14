module boundary_conditions

    use interface
    use define_indexes
    use define_keys
    use define_globals
    use utility, only: util_print_warning
    use utility_interpolate
    use define_settings, only: setting
    use face, only: face_interpolate_bc
    use define_xsect_tables
    use xsect_tables

    implicit none

    private

    public :: bc_step
    public :: bc_update

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine bc_update()
        !%-----------------------------------------------------------------------------
        integer :: ii
        character(64) :: subroutine_name = "bc_update"
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%boundary_conditions)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        call bc_step()

        !% face interpolation for all BC faces
        if (N_flowBC > 0 .or. N_headBC > 0) then
            call bc_interpolate() ! computes interpolation
            call face_interpolate_bc() ! broadcast interpolation to face & elem arrays
        end if

        if (setting%Debug%File%boundary_conditions) then
            print *, "INFLOW BC"
            print *, "BC times"
            do ii = 1, setting%BC%slots
                print *, BC%flowR_timeseries(:, ii, br_time)
            end do
            print *, "BC values"
            do ii = 1, setting%BC%slots
                print *, BC%flowR_timeseries(:, ii, br_value)
            end do
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

            print *, "HEAD BC"
            print *, "BC times"
            do ii = 1, setting%BC%slots
                print *, BC%headR_timeseries(:, ii, br_time)
            end do
            print *, "BC values"
            do ii = 1, setting%BC%slots
                print *, BC%headR_timeseries(:, ii, br_value)
            end do
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        end if
    end subroutine bc_update
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine bc_step()
        !%-----------------------------------------------------------------------------
        integer :: ii, nidx, lidx
        integer :: tstep_larger_than_resolution
        real(8) :: ttime, tnow, tend
        character(64) :: subroutine_name = "bc_step"
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%boundary_conditions)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        tnow = setting%Time%Now
        tend = setting%Time%End

        if (N_flowBC > 0) then
            do ii = 1, N_flowBC
                nidx = BC%flowI(ii, bi_node_idx)
                if (BC%flowYN(ii,bYN_read_input_file)) then
                    if (BC%flowIdx(ii) == 0) then ! First fetch
                        call bc_fetch_flow(ii)
                    else
                        ttime = BC%flowR_timeseries(ii, BC%flowIdx(ii), br_time) ! Current time slot (upper bound of time interval)
                        if (tnow > ttime) then ! Needs update
                            if (BC%flowIdx(ii) == setting%BC%slots) then
                                call bc_fetch_flow(ii)
                            else
                                tstep_larger_than_resolution = -1
                                do while((tnow > ttime) .and. (BC%flowIdx(ii) < setting%BC%slots))
                                    BC%flowIdx(ii) = BC%flowIdx(ii) + 1
                                    ttime = BC%flowR_timeseries(ii, BC%flowIdx(ii), br_time)
                                    if ((BC%flowIdx(ii) == setting%BC%slots) .and. (tnow > ttime)) then
                                        call bc_fetch_flow(ii)
                                        ttime = BC%flowR_timeseries(ii, BC%flowIdx(ii), br_time)
                                    end if
                                    tstep_larger_than_resolution = tstep_larger_than_resolution + 1
                                end do
                                if ((tstep_larger_than_resolution > 0) .and. setting%Output%Warning) then
                                    call util_print_warning("Warning (bc.f08): The flow boundary condition for node " &
                                    // node%Names(nidx)%str // " has a higher resolution than the time step")
                                end if
                            end if
                        end if
                    end if
                else
                    call util_print_warning("Error (bc.f08): The flow boundary condition for node " &
                    // node%Names(nidx)%str // " should always read from an input file")
                    stop 87453
                end if
            end do
        end if

        if (N_headBC > 0) then
            do ii = 1, N_headBC
                nidx = BC%headI(ii, bi_node_idx)
                if (BC%headYN(ii,bYN_read_input_file)) then
                    if (BC%headIdx(ii) == 0) then ! First fetch
                        call bc_fetch_head(ii)
                        !% HACK - we are assuming that outfalls can only have one link upstream
                        lidx = node%I(nidx, ni_Mlink_u1)
                        link%R(lidx, lr_InitialDnstreamDepth) = BC%headR_timeseries(ii, 1, br_value) - node%R(nidx,nr_Zbottom)
                    else
                        ttime = BC%headR_timeseries(ii, BC%headIdx(ii), br_time) ! Current time slot (upper bound of time interval)
                        if (tnow > ttime) then ! Needs update
                            if (BC%headIdx(ii) == setting%BC%slots) then
                                call bc_fetch_head(ii)
                            else
                                tstep_larger_than_resolution = -1
                                do while((tnow > ttime) .and. (BC%headIdx(ii) < setting%BC%slots))
                                    BC%headIdx(ii) = BC%headIdx(ii) + 1
                                    ttime = BC%headR_timeseries(ii, BC%headIdx(ii), br_time)
                                    if ((BC%headIdx(ii) == setting%BC%slots) .and. (tnow > ttime)) then
                                        call bc_fetch_head(ii)
                                        ttime = BC%headR_timeseries(ii, BC%headIdx(ii), br_time)
                                    end if
                                    tstep_larger_than_resolution = tstep_larger_than_resolution + 1
                                end do
                                if ((tstep_larger_than_resolution > 0) .and. setting%Output%Warning) then
                                    call util_print_warning("Warning (bc.f08): The head boundary condition for node " &
                                    // node%Names(nidx)%str // " has a higher resolution than the time step")
                                end if
                            end if
                        end if
                    end if  
                end if
            end do
        end if

        if (setting%Debug%File%boundary_conditions)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine bc_step
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine bc_fetch_flow(bc_idx)
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: bc_idx
        integer             :: ii, NN
        real(8)             :: new_inflow_time
        real(8)             :: new_inflow_value
        character(64)       :: subroutine_name = "bc_fetch_flow"
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%boundary_conditions)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        NN = setting%BC%slots

        if (BC%flowIdx(bc_idx) == 0) then ! First fetch
            BC%flowR_timeseries(bc_idx, 1, br_time) = setting%Time%Start
            BC%flowR_timeseries(bc_idx, 1, br_value) = interface_get_flowBC(bc_idx, setting%Time%Start)
        else ! last value becomes first
            BC%flowR_timeseries(bc_idx, 1, br_time) = BC%flowR_timeseries(bc_idx, NN, br_time)
            BC%flowR_timeseries(bc_idx, 1, br_value) = BC%flowR_timeseries(bc_idx, NN, br_value)
        end if

        new_inflow_time = setting%Time%Start
        do ii = 2, NN
            new_inflow_time = min(setting%Time%End, interface_get_next_inflow_time(bc_idx, new_inflow_time))
            BC%flowR_timeseries(bc_idx, ii, br_time) = new_inflow_time
            BC%flowR_timeseries(bc_idx, ii, br_value) = interface_get_flowBC(bc_idx, new_inflow_time)
            if (new_inflow_time == setting%Time%End) exit
        end do
        BC%flowIdx(bc_idx) = 2

        if (setting%Debug%File%boundary_conditions) then
            do ii = 1, NN
                write(*, "(*(G0.4 : ','))") BC%flowR_timeseries(bc_idx, ii, :)
            end do
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        end if

    end subroutine bc_fetch_flow
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine bc_fetch_head(bc_idx)
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: bc_idx
        integer             :: ii, NN
        real(8)             :: new_head_time
        real(8)             :: new_head_value, normDepth, critDepth
        integer, pointer    :: nodeIdx, faceIdx, elemUpIdx
        character(64)       :: subroutine_name = "bc_fetch_head"
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%boundary_conditions)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        NN = setting%BC%slots

        if (BC%headIdx(bc_idx) == 0) then ! First fetch

            BC%headR_timeseries(bc_idx, 1, br_time) = setting%Time%Start
            BC%headR_timeseries(bc_idx, 1, br_value) = interface_get_headBC(bc_idx, setting%Time%Start)

        else ! last value becomes first
            BC%headR_timeseries(bc_idx, 1, br_time) = BC%headR_timeseries(bc_idx, NN, br_time)
            BC%headR_timeseries(bc_idx, 1, br_value) = BC%headR_timeseries(bc_idx, NN, br_value)
        end if

        do ii = 2, NN
            new_head_time = min(setting%Time%End, interface_get_next_head_time(bc_idx, setting%Time%Start))
            BC%headR_timeseries(bc_idx, ii, br_time) = new_head_time
            BC%headR_timeseries(bc_idx, ii, br_value) = interface_get_headBC(bc_idx, new_head_time)
            if (new_head_time == setting%Time%End) exit
        end do
        BC%headIdx(bc_idx) = 2

        if (setting%Debug%File%boundary_conditions) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine bc_fetch_head
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine bc_interpolate()
    !%-----------------------------------------------------------------------------
    !% Description:
    !% This subroutine is for boundary condition interpolation.
    !% Base on the time passed from the time loop, we interpolate (linear interpolation for now)
    !% the boundary condition to get the corresponding value.
    !%-----------------------------------------------------------------------------
        real(8) :: tnow, normDepth, critDepth
        integer :: ii, slot_idx, upper_idx, lower_idx 
        integer, pointer :: nodeIdx, faceIdx, elemUpIdx
        character(64) :: subroutine_name = 'bc_interpolate'
    !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%boundary_conditions)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        tnow = setting%Time%Now

        do ii=1, N_flowBC
            slot_idx = BC%flowIdx(ii)
            upper_idx = slot_idx
            lower_idx = slot_idx - 1
            !% Find the cloest index first, assign it to lower_idx for now
            if (BC%flowR_timeseries(ii, lower_idx, br_time) == tnow) then
                !% no need to do the interpolation, directly take the existing BC data
                BC%flowRI(ii) = BC%flowR_timeseries(ii, lower_idx, br_value)
            else if (lower_idx > 0) then
                if ( BC%flowR_timeseries(ii, lower_idx, br_value) == BC%flowR_timeseries(ii, upper_idx, br_value)) then
                    BC%flowRI(ii) = BC%flowR_timeseries(ii, lower_idx, br_value)
                    !% constant value, no need to do the interpolation
                else
                    !% interpolation step
                    if (.not. setting%BC%disableInterpolation) then
                        BC%flowRI(ii) = util_interpolate_linear( &
                            tnow, &
                            BC%flowR_timeseries(ii, lower_idx, br_time), &
                            BC%flowR_timeseries(ii, upper_idx, br_time), &
                            BC%flowR_timeseries(ii, lower_idx, br_value), &
                            BC%flowR_timeseries(ii, upper_idx, br_value))
                    else
                        BC%flowRI(ii) = BC%flowR_timeseries(ii, upper_idx, br_value) !% will be zero
                    end if
                end if
            end if
        end do

        do ii=1, N_headBC
            nodeIdx     => BC%headI(ii,bi_node_idx)
            faceIdx     => BC%headI(ii,bi_face_idx)
            elemUpIdx   => faceI(faceIdx,fi_Melem_uL)

            if ((BC%headI(ii,bi_subcategory) == BCH_fixed)   .or. &
                (BC%headI(ii,bi_subcategory) == BCH_tseries) .or. &
                (BC%headI(ii,bi_subcategory) == BCH_tidal)) then
                slot_idx = BC%headIdx(ii)
                upper_idx = slot_idx
                lower_idx = slot_idx - 1
                !% Find the cloest index first, assign it to lower_idx just for now
                if (BC%headR_timeseries(ii, lower_idx, br_time) == tnow) then
                    !% no need to do the interpolation, directly take the existing BC data
                    BC%headRI(ii) = BC%headR_timeseries(ii, lower_idx, br_value)
                else if (lower_idx .ne. 0) then
                    if ( BC%headR_timeseries(ii, lower_idx, br_value) == BC%headR_timeseries(ii, upper_idx, br_value)) then
                        BC%headRI(ii) = BC%headR_timeseries(ii, lower_idx, br_value)
                        !% constant value, no need to do the interpolation
                    else
                        BC%headRI(ii) = util_interpolate_linear( &
                            tnow, &
                            BC%headR_timeseries(ii, lower_idx, br_time), &
                            BC%headR_timeseries(ii, upper_idx, br_time), &
                            BC%headR_timeseries(ii, lower_idx, br_value), &
                            BC%headR_timeseries(ii, upper_idx, br_value))
                    end if
                end if

            else if (BC%headI(ii,bi_subcategory) == BCH_normal) then

                !% for normal dnBC, if the node has a invert,
                !% the depth in the node is zero
                if (node%R(nodeIdx,nr_Zbottom) > zeroR) then
                    BC%headRI(ii) = faceR(faceIdx,fr_Zbottom)
                else
                    if (elemI(elemUpIdx,ei_elementType) == CC) then
                        BC%headRI(ii) = elemR(elemUpIdx,er_Depth) + faceR(faceIdx,fr_Zbottom)
                    else
                        !% for normal dnBC, if the upstream link is not CC (i.e. weir, orifice etc)
                        !% the depth in the node is zero
                        BC%headRI(ii) =  faceR(faceIdx,fr_Zbottom)
                    end if
                end if

            else if (BC%headI(ii,bi_subcategory) == BCH_free) then

                !% for free dnBC, if the node has a invert,
                !% the depth in the node is zero
                if (node%R(nodeIdx,nr_Zbottom) > zeroR) then
                    BC%headRI(ii) = faceR(faceIdx,fr_Zbottom)
                else
                    if (elemI(elemUpIdx,ei_elementType) == CC) then
                        normDepth = elemR(elemUpIdx,er_Depth)
                        critDepth = bc_get_CC_critical_depth(elemUpIdx)

                    else
                        !% for free dnBC, if the upstream link is not CC (i.e. weir, orifice etc)
                        !% the depth in the node is zero
                        BC%headRI(ii) =  faceR(faceIdx,fr_Zbottom)
                    end if

                end if

                stop 55474
            else
                call util_print_warning("Error (bc.f08): Unknown downstream boundary condition type at " &
                    // node%Names(nodeIdx)%str // " node")
                    stop 86474
            end if
        end do

        if (setting%Debug%File%boundary_conditions) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine bc_interpolate
!%
!%==========================================================================
!%==========================================================================
!%
function bc_get_CC_critical_depth(elemIdx) result (criticalDepth)
    !%-----------------------------------------------------------------------------
    !% Description:
    !% gets the critical depth of a CC element
    !% this code has been adapted from SWMM5-C code
    !%-----------------------------------------------------------------------------
        integer, intent(in)  :: elemIdx
        real(8)              :: criticalDepth
        integer, pointer     :: elemGeometry 
        real(8), pointer     :: Flowrate, FullDepth, BreadthMax, grav, FullArea
        real(8)              :: Q2g, critDepthEstimate, ratio
        character(64) :: subroutine_name = 'bc_get_CC_critical_depth'
    !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%boundary_conditions)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        elemGeometry => elemI(elemIdx,ei_geometryType)
        Flowrate     => elemR(elemIdx,er_Flowrate)
        FullDepth    => elemR(elemIdx,er_FullDepth)
        FullArea     => elemR(elemIdx,er_FullArea)
        BreadthMax   => elemR(elemIdx,er_BreadthMax)
        grav         => setting%constant%gravity

        Q2g = Flowrate**twoR / grav

        if (Q2g == zeroR) then
            criticalDepth =  zeroR
            
        else
            select case (elemGeometry)

                case (rectangular)
                    criticalDepth = (Q2g/BreadthMax**twoR)**onethirdR
                    
                case default
                    !% first estimate Critical Depth for an equivalent circular conduit
                    critDepthEstimate = min(1.01*(Q2g/FullDepth)**onefourthR, FullDepth)

                    !% find ratio of conduit area to equiv. circular area
                    ratio = FullArea/((pi/fourR)*FullDepth**twoR)

                    if ((ratio >= onehalfR) .and. (ratio <= twoR)) then
                        

                        criticalDepth = bc_critDepth_enum (elemIdx, Flowrate, critDepthEstimate)
                        
                    else
                        print*, 'riddlers method has not been developed yet'
                        stop 123541
                        
                    end if

            end select
        end if
        print*, criticalDepth, 'criticalDepth'
        stop 121345

        if (setting%Debug%File%boundary_conditions) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

end function bc_get_CC_critical_depth
!%
!%==========================================================================
!%==========================================================================
!%
function bc_critDepth_enum (elemIdx, Q, yEstimate) result (criticalDepth)
    !%-----------------------------------------------------------------------------
    !% Description:
    !% gets the critical depth of a CC by enumeration method
    !% this code has been adapted from SWMM5-C code
    !%-----------------------------------------------------------------------------
        integer, intent(in)  :: elemIdx
        real(8), intent(in)  :: Q, yEstimate
        real(8)              :: criticalDepth
        real(8), pointer     :: FullDepth
        real(8)              :: dY, Y, Q0, QC
        integer              :: i1, ii
        character(64) :: subroutine_name = 'bc_critDepth_enum'
    !%-----------------------------------------------------------------------------

    FullDepth => elemR(elemIdx,er_FullDepth)

    !% find the critical depth from SWMM5 ennumeration method
    dY = FullDepth/25.0
    i1 = int(yEstimate/dY)
    Q0 = bc_get_critical_flow(elemIdx, i1*dY)

    if (Q0 < Q) then
        criticalDepth = FullDepth
        do  ii = i1+1,25
            QC = bc_get_critical_flow(elemIdx,ii*dY)
            if (QC > Q) then
                criticalDepth = ((Q-Q0)/(QC-Q0)+ real((ii-1),8))*dY
                return
            end if
            Q0 = QC
        end do
    else
        criticalDepth = zeroR
        do  ii = i1-1,0,-1
            QC = bc_get_critical_flow(elemIdx,ii*dY)
            if (QC < Q) then
                criticalDepth = ((Q-QC)/(Q0-QC)+ real((ii),8))*dY
                return
            end if
            Q0 = QC
        end do
    end if

end function bc_critDepth_enum
!%
!%==========================================================================
!%==========================================================================
!%
function bc_get_critical_flow (elemIdx, Depth) result (Flow)
    !%-----------------------------------------------------------------------------
    !% Description:
    !% gets the critical flow of a CC element for a given depth
    !%-----------------------------------------------------------------------------
        integer, intent(in)  :: elemIdx
        real(8), intent(in)  :: Depth 
        real(8)              :: Flow
        integer, pointer     :: elemGeometry
        real(8), pointer     :: fullDepth, fullArea, grav
        real(8)              :: YoverYfull, Area, Topwidth
        character(64) :: subroutine_name = 'bc_get_critical_flow'
    !%-----------------------------------------------------------------------------

    elemGeometry => elemI(elemIdx,ei_geometryType)
    fullDepth    => elemR(elemIdx,er_FullDepth)
    fullArea     => elemR(elemIdx,er_FullArea)
    grav         => setting%constant%gravity

    select case (elemGeometry)

        case (circular)
            YoverYfull = Depth/fullDepth
            Area     = fullArea * xsect_table_lookup_singular (YoverYfull, ACirc, NACirc)
            Topwidth = max(fullDepth * xsect_table_lookup_singular (YoverYfull, TCirc, NTCirc), &
                        setting%ZeroValue%Topwidth)

            Flow  = Area * sqrt(grav*Area/Topwidth)

        case default
            print*, 'This default case should not reach'
            stop 84198
     end select

end function bc_get_critical_flow
!%
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module boundary_conditions