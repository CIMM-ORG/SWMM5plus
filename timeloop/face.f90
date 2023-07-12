module face
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Gets face values based on adjacent elements
    !%
    !% Methods:
    !% Time-scale interpolation for geometry and flowrate, linear space
    !% interpolation for head. Junction JB and BC have special forcing.
    !%==========================================================================

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use adjust
    use geometry
    use jump
    use pack_mask_arrays, only: pack_CC_zeroDepth_interior_faces, pack_CC_zeroDepth_shared_faces
    use utility_profiler
    use utility, only: util_sign_with_ones
    use utility_crash, only: util_crashpoint

    implicit none

    private

    public :: face_pull_facedata_to_JBelem
    public :: face_push_elemdata_to_face
    public :: face_push_diag_adjacent_data_to_face
    public :: face_interpolation
    public :: face_interpolate_bc
    public :: face_force_JBelem_to_face
    public :: face_shared_face_sync

    public :: face_flowrate_for_openclosed_elem

    public :: face_zeroDepth


    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine face_zeroDepth (fp_downstream, fp_upstream, fp_bothsides)
        !%------------------------------------------------------------------
        !% Description:
        !% sets zero depth geometry values on faces depending on the
        !% adjacent elements
        !% --- For CC
        !% arguments are fp_CC_downstream_is_zero_IorS,
        !% fp_CC_upstream_is_zero_IorS, fp_CC_bothsides_are_zero_IorS
        !% --- For JB
        !% arguments are fp_JB_downstream_is_zero_IorS,
        !% fp_JB_upstream_is_zero_IorS, fp_JB_bothsides_are_zero_IorS
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: fp_downstream, fp_upstream, fp_bothsides
        !%------------------------------------------------------------------

        !% --- geometry (interior)
        call face_zeroDepth_geometry_interior(fp_downstream)
        call face_zeroDepth_geometry_interior(fp_upstream)
        call face_zeroDepth_geometry_interior(fp_bothsides)

        !% --- geometry (shared)
        !% sync all processors before the starting
        sync all

        call face_zeroDepth_geometry_shared(fp_downstream)
        call face_zeroDepth_geometry_shared(fp_upstream)
        call face_zeroDepth_geometry_shared(fp_bothsides)

        !% --- flowrate (interior)
        call face_zeroDepth_flowrates_interior(fp_downstream)
        call face_zeroDepth_flowrates_interior(fp_upstream)
        call face_zeroDepth_flowrates_interior(fp_bothsides)

        !% --- flowrate (shared)
        call face_zeroDepth_flowrates_shared(fp_downstream)
        call face_zeroDepth_flowrates_shared(fp_upstream)
        call face_zeroDepth_flowrates_shared(fp_bothsides)

        !% --- sync all processors before exiting
        sync all

    end subroutine face_zeroDepth
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine face_push_elemdata_to_face &
        (thisPCol, frCol, erCol, elemXR, UpstreamFaceTF)
        !%------------------------------------------------------------------
        !% Description
        !% pushes one column (erCol) of elemXR(:,:) array to one fr_.. column for the
        !% elemP packed column. If UpstreamFaceTF is true the push is to the
        !% upstream face, otherwise to downstream 
        !% NOTE this can give a segmentation fault if the face map for the
        !% any of the packed elements gives a nullvalueI. To prevent this from
        !% happening, make sure that calls to this including JB elements are
        !% done either using upstream only or downstream only.
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(in) :: elemXR(:,:)
            integer, intent(in) :: thisPCol, frCol, erCol
            logical, intent(in) :: UpstreamFaceTF
            integer, pointer :: Npack, thisP(:)
            integer :: eiMface
        !%------------------------------------------------------------------
        !% Preliminaries
            Npack => npack_elemP(thisPCol)
            if (Npack < 1) return
            thisP => elemP(1:Npack,thisPCol)
            if (UpstreamFaceTF) then 
                eiMface = ei_Mface_uL 
            else 
                eiMface = ei_Mface_dL 
            end if
        !%------------------------------------------------------------------
        !% --- push data
        faceR(elemI(thisP,eiMface),frCol) = elemXR(thisP,erCol)

        !% --- check if the elem data has either been pushed 
        !%     to a shared face. if so then mark that face
        where (faceYN(elemI(thisP,eiMface),fYN_isSharedFace))
            faceYN(elemI(thisP,eiMface),fYN_isSharedFaceDiverged) = .true.
        end where

    end subroutine face_push_elemdata_to_face
!%
!%==========================================================================    
!%==========================================================================
!%
    subroutine face_push_diag_adjacent_data_to_face (thisPCol)
        !%------------------------------------------------------------------
        !% Description
        !% Pushes element data into the fr_..._adjacent data fields
        !% thsiPCol must be a packed set of diagnostic elements that are
        !% adjacent to JB elements
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisPCol
            integer, pointer    :: Npack, thisP(:)
            integer :: ii, kk, ff
        !%------------------------------------------------------------------
            Npack => npack_elemP(thisPCol)
            thisP => elemP(1:Npack,thisPCol)
            if (Npack < 1) return 
        !%------------------------------------------------------------------   
        !% --- cycle through a set of diagnostic elements
        do ii=1,Npack
            !% -- cycle through upstream and downstream faces
            do kk=1,2
                if (kk==1) then !% -- upstream face
                    ff = elemI(thisP(ii),ei_Mface_uL)
                    if (.not. faceYN(ff,fYN_isUpstreamJBFace)) cycle !% if not JB adjacent
                else !% -- downstream face
                    ff = elemI(thisP(ii),ei_Mface_dL)
                    if (.not. faceYN(ff,fYN_isDownstreamJBFace)) cycle !% if not JB adjacent
                end if

                !% --- set the adjacent element value storage on the face
                select case (elemI(thisP(ii),ei_elementType))
                    case (weir)
                        faceR(ff,fr_Zcrest_Adjacent) = elemSR(thisP(ii),esr_Weir_Zcrest)
                        faceR(ff,fr_dQdH_Adjacent)   = elemSR(thisP(ii),esr_Weir_dQdHe)
                    case (orifice)
                        faceR(ff,fr_Zcrest_Adjacent) = elemSR(thisP(ii),esr_Orifice_Zcrest)
                        faceR(ff,fr_dQdH_Adjacent)   = elemSR(thisP(ii),esr_Orifice_dQdHe)
                    case (outlet)
                        faceR(ff,fr_Zcrest_Adjacent) = elemSR(thisP(ii),esr_Outlet_Zcrest)
                        faceR(ff,fr_dQdH_Adjacent)   = elemSR(thisP(ii),esr_Outlet_dQdHe)
                    case (pump)
                        faceR(ff,fr_Zcrest_Adjacent) = elemSR(thisP(ii),esr_Pump_Zcrest)
                        faceR(ff,fr_dQdH_Adjacent)   = elemSR(thisP(ii),esr_Pump_dQdHp)
                    case default 
                        print *, 'CODE ERROR: unexpected case default'
                        call util_crashpoint(3111987)
                end select
                !% --- check if the elem data has either been pushed 
                !%     to a shared face. if so then mark that face
                if (faceYN(ff,fYN_isSharedFace)) then
                    faceYN(ff,fYN_isSharedFaceDiverged) = .true.
                end if
            end do
        end do    

    end  subroutine face_push_diag_adjacent_data_to_face
!%
!%==========================================================================    
!%==========================================================================
!%    
    subroutine face_pull_facedata_to_JBelem (thisPcol, frCol, eDataOut)
        !%------------------------------------------------------------------
        !% Description:
        !% Pulls data from faceR(:,frCol) to eDataOut for the elements
        !% in elemP(:,thisPcol). eDataOut should be (e.g.) elemR(:,er_Flowrate)
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in)    :: thisPcol, frCol
            real(8), intent(inout) :: eDataOut(:)
            integer, pointer       :: NPack, thisP(:)
        !%------------------------------------------------------------------
        !% Preliminaries
            Npack => npack_elemP(thisPcol)
            if (Npack < 1) return 
            thisP => elemP(1:Npack,thisPcol)
        !%------------------------------------------------------------------

        where (elemSI(thisP,esi_JunctionBranch_IsUpstream) == oneI)
            !% --- upstream JB branch
            eDataOut(thisP) = faceR(elemI(thisP,ei_Mface_uL),frCol)
        elsewhere
            !% --- downstream JB branch
            eDataOut(thisP) = faceR(elemI(thisP,ei_Mface_dL),frCol)
        endwhere

    end subroutine face_pull_facedata_to_JBelem    
!%
!%==========================================================================    
!%==========================================================================
!%     
    subroutine face_interpolation &
        (facecol, Gyn, Hyn, Qyn, skipJump, skipZeroAdjust)
        !%------------------------------------------------------------------
        !% Description:
        !% Interpolates faces from elements
        !% NOTE -- calls to this subroutine CANNOT be within a conditional
        !% that would prevent it from being called by all images. That is,
        !% this subroutine MUST be called by all images (even with a null)
        !% to make sure that the images can be synced before sharing.
        !% Gyn, Hyn, Qyn = true will interpolate geometry, head, flowrate, respectively
        !% skipZeroAdjust = true skips the adjustment for zero cells, which
        !%    must be skipped for interpolation of JM to ensure mass conservation
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in)  :: faceCol !, whichTM
            logical, intent(in)  :: Qyn, Hyn, Gyn
            logical, intent(in)  :: skipZeroAdjust, skipJump
            logical              :: isBConly !, isTM
            
            character(64) :: subroutine_name = 'face_interpolation'
        !%-------------------------------------------------------------------
        !% Preliminaries    
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            
            if (setting%Profile%useYN) call util_profiler_start (pfc_face_interpolation)
        !%--------------------------------------------------------------------
        isBConly = .false.

        !% --- face reconstruction of all the interior faces
        call face_interpolation_interior (faceCol, Gyn, Hyn, Qyn, skipJump)

        !% --- face reconstruction of all the shared faces
        call face_interpolation_shared (faceCol, Gyn, Hyn, Qyn, .true.)

        call face_interpolate_bc (isBConly)

        !%-------------------------------------------------------------------
        !% Closing
            if (setting%Profile%useYN) call util_profiler_stop (pfc_face_interpolation)

            if (setting%Debug%File%face)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine face_interpolation
!%    
!%==========================================================================    
!%==========================================================================
!%
    subroutine face_interpolate_bc(isBConly)
        !%------------------------------------------------------------------
        !% Description:
        !% Interpolates all data to upstream and downstream faces
        !% when isBConly == true, this only handles BC enforcement
        !%-------------------------------------------------------------------
        !% Declarations:
            logical, intent(in) :: isBConly
            character (64) :: subroutine_name = 'face_interpolate_bc'
        !%-------------------------------------------------------------------
            if (setting%Debug%File%face)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------

        if ((N_nBCup > 0) .or. (N_nJ1 > 0)) call face_interpolation_upBC(isBConly)

        if (N_nBCdn > 0) call face_interpolation_dnBC(isBConly)

        !%-------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine face_interpolate_bc
!%    
!%==========================================================================    
!%==========================================================================
!%    
    subroutine face_force_JBelem_to_face (thisColP, isJBupstreamYN)
        !%------------------------------------------------------------------
        !% Description
        !% Forces faces adjacent to JB to JB values without interpolation
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisColP
            logical, intent(in) :: IsJBupstreamYN !% true if JB is upstream of JM
            integer, pointer    :: thisJM(:), Npack
            integer :: mm, ei_Mface, kstart
        !%------------------------------------------------------------------
        !% Preliminaries:
            Npack => npack_elemP(thisColP)
            if (Npack < 1) return
            thisJM => elemP(1:Npack,thisColP) 
            if (isJBupstreamYN) then 
                ei_Mface = ei_Mface_uL
                kstart = oneI !% JB upstream of JM
            else
                ei_Mface = ei_Mface_dL
                kstart = twoI !% JB downstream of JM
            end if
        !%------------------------------------------------------------------

        do mm=1,Npack
            !% --- push junction JB flowrate values back to faces  
            call face_force_JBvalues (fr_Flowrate, er_Flowrate, ei_Mface, thisJM(mm), kstart) 

            !% ---add junction JB Delta Q values to face value caused by interpolation
            call face_force_JBvalues (fr_DeltaQ, er_DeltaQ, ei_Mface, thisJM(mm), kstart) 
        end do

    end subroutine face_force_JBelem_to_face
!%    
!%==========================================================================    
!% PRIVATE
!%==========================================================================
!%
    subroutine face_interpolation_upBC(isBConly)
        !%------------------------------------------------------------------
        !% Description:
        !% Interpolates data to all upstream BC faces
        !% When called with isBConly == true, only does the BC update
        !%------------------------------------------------------------------
        !% Declarations:
            logical, intent(in) :: isBConly
            integer :: fGeoSetU(2), fGeoSetD(2), eGeoSet(2)
            integer :: ii
            integer, pointer :: edn(:), idx_P(:), fdn(:)
            integer, pointer :: idx_fBC(:), idx_fJ1(:), idx_fBoth(:)
            integer, pointer :: npackBC, npackJ1, npackBoth
            character(64) :: subroutine_name = 'face_interpolation_upBC'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%boundary_conditions)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% Aliases:
        !% For the head/geometry at the upstream faces, we directly take the dnwnstream element
        !% So there is no eup for upstream BCs
            edn       => faceI(:,fi_Melem_dL)
            fdn       => elemI(:,ei_Mface_dL)   
            npackBC   => npack_faceP(fp_BCup)
            idx_fBC   => faceP(1:npackBC,fp_BCup)
            npackJ1   => npack_faceP(fp_J1)
            idx_fJ1   => faceP(1:npackJ1,  fp_J1)
            npackBoth => npack_faceP(fp_J1_BCup)
            idx_fBoth => faceP(1:npackBoth,  fp_J1_BCup)
            idx_P     => BC%P%BCup(:)
        !%-------------------------------------------------------------------
        !% --- enforce stored inflow BC    
        if (npackBC > 0) then
            faceR(idx_fBC, fr_Flowrate) = BC%flowR(idx_P,br_value)
        end if

        !% --- enforce zero flow on J1 faces
        if (npackJ1 > 0) then
            faceR(idx_fJ1, fr_Flowrate)  = zeroR
            faceR(idx_fJ1, fr_Velocity_u) = zeroR
            faceR(idx_fJ1, fr_Velocity_d) = zeroR
        end if

        !% --- no upstream BC
        if (npackBoth < 1) return 

        !% --- update geometry data (don't do on a BC-only call)
        if (.not. isBConly) then
            !% --- Define sets of points for the interpolation, we are going from
            !%     the elements to the faces.       
            fGeoSetU = [fr_Area_u, fr_Depth_u]
            fGeoSetD = [fr_Area_d, fr_Depth_d]
            eGeoSet  = [er_Area,   er_Depth]

            !% --- Copying geometric data from elements to the BC/J1 faces
            do ii = 1,size(fGeoSetU)
                faceR(idx_fBoth,fGeoSetD(ii)) = elemR(edn(idx_fBoth),eGeoSet(ii)) 
                !% --- upstream side of face matches downstream (no jump)
                faceR(idx_fBoth,fGeoSetU(ii)) = faceR(idx_fBoth,fGeoSetD(ii))
            end do

            !% --- HACK: copy the preissmann number as well
            !%     Note that for static slot this will always be unity
            faceR(idx_fBoth,fr_Preissmann_Number) = elemR(edn(idx_fBoth),er_Preissmann_Number) 
            
            !% --- gradient extrapolation for head at infow
            if (npackBC > 0) then 
                faceR(idx_fBC, fr_Head_d) = elemR(edn(idx_fBC),er_Head)        &
                                          + elemR(edn(idx_fBC),er_Head)        &
                                          - faceR(fdn(edn(idx_fBC)),fr_Head_d) 
                !% --- use larger depth of the direct injection (from GeoSet above)
                !%     and depth implied by head gradient extrapolation.                          
                faceR(idx_fBC,fr_Depth_d) = max(faceR(idx_fBC,fr_Depth_d),  &
                                                  faceR(idx_fBC,fr_Head_d)    &
                                                - faceR(idx_fBC,fr_Zbottom))                          
            end if

            !% --- zero head gradient at J1 cells
            if (npackJ1 > 0) then
                faceR(idx_fJ1,fr_Head_d) = elemR(edn(idx_fJ1),er_Head) 
            end if

            !% --- head/ depth on downstream side of face is copied to upstream side (no jump)
            faceR(idx_fBoth,fr_Head_u)  = faceR(idx_fBoth,fr_Head_d)
            faceR(idx_fBoth,fr_Depth_u) = faceR(idx_fBoth,fr_Depth_d)
           
            !% --- ensure face area_u is not smaller than zerovalue
            if (npackBC > 0) then
                where (faceR(idx_fBC,fr_Area_d) <= setting%ZeroValue%Area)
                    faceR(idx_fBC,fr_Area_d)     = setting%ZeroValue%Area
                end where
                where (faceR(idx_fBC,fr_Area_u) <= setting%ZeroValue%Area)
                    faceR(idx_fBC,fr_Area_u)     = setting%ZeroValue%Area
                endwhere
            end if

        end if
        
        !% --- inflow velocity (provides momentum transported into downstream element)
        if (npackBC > 0) then 
            !% --- set velocity based on flowrate
            faceR(idx_fBC,fr_Velocity_d) = faceR(idx_fBC,fr_Flowrate)/faceR(idx_fBC,fr_Area_d)
            !% --- limit Velocity of inflow to Froude Number = 1
            faceR(idx_fBC,fr_Velocity_d) = min(faceR(idx_fBC,fr_Velocity_d), &
                                               sqrt( setting%Constant%gravity &
                                                     * faceR(idx_fBC,fr_Depth_d) ))

            !% --- set upstream velocity to downstream
            faceR(idx_fBC,fr_Velocity_u) = faceR(idx_fBC,fr_Velocity_d)  

            !%  --- reset velocity using high velocity limiter
            if (setting%Limiter%Velocity%UseLimitMaxYN) then            
                where(abs(faceR(idx_fBC,fr_Velocity_d))  > setting%Limiter%Velocity%Maximum)
                    faceR(idx_fBC,fr_Velocity_d) = sign(0.99d0 * setting%Limiter%Velocity%Maximum, &
                        faceR(idx_fBC,fr_Velocity_d))
                endwhere
                where(abs(faceR(idx_fBC,fr_Velocity_u))  > setting%Limiter%Velocity%Maximum)
                    faceR(idx_fBC,fr_Velocity_u) = sign(0.99d0 * setting%Limiter%Velocity%Maximum, &
                        faceR(idx_fBC,fr_Velocity_u))
                endwhere
            end if

        end if

        !%-------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%boundary_conditions) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine face_interpolation_upBC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine face_interpolation_dnBC(isBConly)
        !%------------------------------------------------------------------
        !% Description:
        !% Interpolates data to all downstream boundary faces
        !% When called with isBConly == true, only does the BC update
        !%------------------------------------------------------------------
        !% Declarations
            logical, intent(in) :: isBConly !% = true when only head BC is enforced
            
            integer, pointer :: nBC, eup, fidx(:), thisF, nodeIdx
            real(8), pointer :: eHead(:), eDepth(:), eArea(:)
        
            integer :: ii
            real(8) :: maxQin, maxQout, headdif, critDepth, critArea, normDepth
        !%------------------------------------------------------------------
        !% Aliases
            nBC           => npack_faceP(fp_BCdn)
            fidx          => faceP(1:nBC,fp_BCdn)    
            eHead         => elemR(:,er_Head)
            eArea         => elemR(:,er_Area)
            eDepth        => elemR(:,er_Depth)
        !%------------------------------------------------------------------

        !% --- cycle through the BC to set depth and head on face
        do ii=1,nBC
            thisF   => fidx(ii)
            nodeIdx => BC%HeadI(ii,bi_node_idx)
            eup     => faceI(thisF,fi_Melem_uL)

            select case (BC%headI(ii,bi_subcategory))
                case (BCH_tidal, BCH_tseries, BCH_fixed)
                    !% --- The head valeu is set in bc_interpolate_head

                    !% --- upstream face head is the BC, mustat least Zbottom
                    faceR(thisF, fr_Head_u)  = max(BC%headR(ii,br_value),faceR(thisF,fr_Zbottom)) 
                    !% --- depth at face
                    faceR(thisF, fr_Depth_u) = max(BC%headR(ii,br_value) - faceR(thisF,fr_Zbottom), zeroR)
                    !% --- adjust for reference head
                    faceR(thisF, fr_Head_u)  = faceR(thisF, fr_Head_u) - setting%Solver%ReferenceHead

                    !% --- the downstream side of face is the same as the upstream face (unless gate, see below)
                    faceR(thisF, fr_Head_d)  = faceR(thisF, fr_Head_u)
                    faceR(thisF, fr_Depth_d) = faceR(thisF, fr_Depth_u)

                    !% --- for a flap gate on a BC with higher head downstream
                    if ( (BC%headYN(ii,bYN_hasFlapGate)) .and. (eHead(eup)  < faceR(thisF,fr_Head_u))) then
                    !% --- reset the head on the upstream side of face for closed gate
                        faceR(thisF, fr_Head_u)  = max(eHead(eup)+setting%Solver%ReferenceHead, faceR(thisF,fr_Zbottom))
                        faceR(thisF, fr_Depth_u) = eHead(eup) - faceR(thisF,fr_Zbottom)
                        !% --- adjust for reference head
                        faceR(thisF, fr_Head_u)  = faceR(thisF, fr_Head_u) - setting%Solver%ReferenceHead
                        !% faceR(thisF, fr_Head_d)  is unchanged
                        !% faceR(thisF, fr_Depth_d) is unchanged
                    endif

                    !% --- consistent areas
                    faceR(thisF,fr_Area_u)  = geo_area_from_depth_singular   &
                                    (eUp, faceR(thisF, fr_Depth_u), setting%ZeroValue%Area)
                    faceR(thisF,fr_Area_d)  = faceR(thisF,fr_Area_u)

                case (BCH_normal)

                    if (elemI(eup,ei_elementType) == CC) then
                        !% --- Error check, normal depth is infinite for adverse slope
                        !%     Use setting%Eps%Machine so that slope must be greater than precision 
                        if (elemR(eup,er_BottomSlope) .le. onehundredR*setting%Eps%Machine) then
                            print *, 'USER CONFIGURATION ERROR: a NORMAL OUTFALL must be connected to an...'
                            print *, '...conduit/channel element with non-zero, positive bottom slope.'
                            print *, 'Problem for Outfall ',trim(node%Names(BC%headI(ii,bi_node_idx))%str)
                            print *, 'Connected to element ', eUp
                            print *, 'Part of link ',trim(  link%Names(elemI(eup,ei_link_Gidx_BIPquick))%str)
                            print *, 'Bottom slope is ',elemR(eup,er_BottomSlope)
                            call util_crashpoint(728474)
                        end if

                        !% --- get the normal depth
                        faceR(thisF,fr_Depth_u) = geo_normaldepth_singular (BC%HeadI(ii,bi_UTidx))

                        if (faceR(thisF,fr_Depth_u) > elemR(eup,er_Head)) then
                            !% --- use upstream head if the normal depth is too large (no backflow allowed)
                            faceR(thisF,fr_Depth_u) = elemR(eup,er_Head) - faceR(thisF,fr_Zbottom)
                        end if

                        !% --- head is the normal depth + Zbottom - referencehead
                        faceR(thisF,fr_Head_u)  = faceR(thisF,fr_Zbottom)                   &
                                                + faceR(thisF,fr_Depth_u)                   &
                                                - setting%Solver%ReferenceHead

                        !% --- the downstream side of face is the same as the upstream face
                        faceR(thisF, fr_Head_d)  = faceR(thisF, fr_Head_u)
                        faceR(thisF, fr_Depth_d) = faceR(thisF, fr_Depth_u)

                        !% --- consistent areas
                        faceR(thisF,fr_Area_u)  = geo_area_from_depth_singular   &
                                                    (eUp, faceR(thisF, fr_Depth_u), setting%ZeroValue%Area)
                        faceR(thisF,fr_Area_d)  = faceR(thisF,fr_Area_u)

                    else
                        print *, 'CODE ERROR: NEED ALGORITHM DESIGN FOR OUTFALL WITH UPSTREAM DIAGNOSTIC ELEMENT'
                        call util_crashpoint(792873)
                        !% for free dnBC, if the upstream link is not CC (i.e. weir, orifice etc)
                        !% the depth in the node is zero
                        !headValue(ii) =  faceR(fIdx,fr_Zbottom)
                    end if

                case (BCH_free)

                    normDepth = geo_normaldepth_singular   (BC%HeadI(ii,bi_UTidx))
                    critDepth = geo_critical_value_singular(BC%HeadI(ii,bi_UTidx),utd_Qcrit_depth_nonuniform)
                    critArea  = geo_critical_value_singular(BC%HeadI(ii,bi_UTidx),utd_Qcrit_area_nonuniform)

                    !% --- select between normal and critical depth depending on conditions
                    if (normDepth > critDepth) then 
                        !% --- GVF mild slope
                        if (eDepth(eup) .ge. normDepth) then 
                            !% --- use normal depth when deep drawdown
                            faceR(thisF,fr_Depth_u) = normDepth
                            faceR(thisF,fr_Area_u)  = geo_area_from_depth_singular   &
                                    (eUp, normDepth, setting%ZeroValue%Area)
                        else
                            if (eDepth(eup) .ge. critDepth) then
                                !% --- use critical depth for drawdown below normal
                                faceR(thisF,fr_Depth_u) = critDepth
                                faceR(thisF,fr_Area_u)  = critArea
                            else
                                !% --- supercritical flow retains same depth
                                faceR(thisF,fr_Depth_u) = eDepth(eup)
                                faceR(thisF,fr_Depth_u) = eArea(eup)
                            end if
                        end if
                    else
                        !% --- GVF steep slope (normDepth < critDepth)
                        if (eDepth(eup) .ge. critDepth) then 
                            !% --- use critical depth for drawdown
                            faceR(thisF,fr_Depth_u) = critDepth
                            faceR(thisF,fr_Area_u)  = critArea
                        else
                            if (eDepth(eup) .ge. normDepth) then 
                                faceR(thisF,fr_Depth_u) = normDepth
                                faceR(thisF,fr_Area_u)  = geo_area_from_depth_singular   &
                                    (eUp, normDepth, setting%ZeroValue%Area)
                            else
                                !% --- supercritical flow retains same depth
                                faceR(thisF,fr_Depth_u) = eDepth(eup)
                                faceR(thisF,fr_Depth_u) = eArea(eup)
                            end if
                        end if
                    end if

                    faceR(thisF,fr_Head_u) = faceR(thisF,fr_Depth_u) &
                                           + faceR(thisF,fr_Zbottom) &
                                           - setting%Solver%ReferenceHead

                    !% --- the downstream side of face is the same as the upstream face
                    faceR(thisF, fr_Head_d)  = faceR(thisF, fr_Head_u)
                    faceR(thisF, fr_Depth_d) = faceR(thisF, fr_Depth_u)
                    faceR(thisF, fr_Area_d)  = faceR(thisF, fr_Area_u)

                case default
                    print *, 'CODE ERROR: unexpected case default'
                    if ((nodeIdx < N_node) .and. (nodeIdx > 0)) then
                        print *, 'Unknown downstream boundary condition type for node'  &
                                // trim(node%Names(nodeIdx)%str)
                    end if
                    call util_crashpoint(864741)

            end select
            
        end do

        !% --- get geometry for face from upstream element
        if (.not. isBConly) then
            do ii=1,nBC 
                thisF   => fidx(ii)
                eup     => faceI(thisF,fi_Melem_uL)
                !% --- ensure face area isnot smaller than zero value
                if (faceR(thisF,fr_Area_u) < setting%ZeroValue%Area) then 
                    faceR(thisF,fr_Area_u) = setting%ZeroValue%Area
                end if
                if (faceR(thisF,fr_Area_d) < setting%ZeroValue%Area) then 
                    faceR(thisF,fr_Area_d) = setting%ZeroValue%Area
                end if
            end do
        end if

        !% --- set the flowrate at the boundary
        if (.not. isBConly) then
            do ii=1,nBC 
                thisF   => fidx(ii)
                eup     => faceI(thisF,fi_Melem_uL)

                !% --- default face flowrate is the upstream element flowrate
                faceR(thisF,fr_Flowrate) = elemR(eup,er_Flowrate)

                !% --- FLAPGATE check for closed flap gate that sets flow to zero
                if (BC%headYN(ii,bYN_hasFlapGate)) then
                        ! print *, 'has flapgate'
                    if (faceR(thisF,fr_Head_u) < faceR(thisF,fr_Head_d) ) then
                        !% --- set BC flow to zero for closed flap gate
                        !%     no outflow until head at face upstream exceeds the gate BC head
                        faceR(thisF,fr_Flowrate) = zeroR
                        ! print *, 'setting face flowrate to zero'
                    else
                        !% --- no change to Q
                    end if
                end if

                !% --- BACKFLOW handle possible backflow
                if (faceR(thisF,fr_Flowrate) < zeroR) then
                    select case (BC%headI(ii,bi_subcategory))
                        case (BCH_tidal, BCH_tseries, BCH_fixed)
                            !% --- backflow is allowed
                        case (BCH_normal, BCH_free)
                            !% -- backflow not allowed
                            faceR(thisF,fr_Flowrate) = zeroR
                        case default
                            print *, 'CODE ERROR: unexpected case default'
                            call util_crashpoint(61098755)
                    end select
                else
                    !% --- no requirement for outflow
                end if

                !% --- FLOW LIMITERS

                !% --- head difference near boundary
                headdif = eHead(eup) - faceR(thisF,fr_Head_u) !% -- should be negative for inflow

                !% --- flowrate limiter for inflow
                !%     to prevent oscillations, limit inflows to 1/4 of the volume 
                !%     that brings the upstream element up to the head BC value
                !%     NOTE this should never occur with normal or free boundaries
                if (faceR(thisF,fr_flowrate) < zeroR) then
                    !% --- head difference driving the boundary flow
                    if (headdif .ge. zeroR) then 
                        !% --- head gradient and inflow are inconsistent, so use zero flow
                        faceR(thisF,fr_Flowrate) = zeroR
                    else
                        !% --- for an inflow
                        if (elemYN(eup,eYN_isZeroDepth)) then 
                            !% --- for zerodepth the head dif is limited to 1/2 the face depth 
                            headdif = -min(-headdif, onehalfR*faceR(thisF,fr_Depth_d))
                        end if
                        !% --- flow that brings upstream element 1/4 of the way to the driving head
                        !%     this should be a negative (inflow) value
                        maxQin  = onefourthR * headdif * elemR(eup,er_Length) &
                                * elemR(eup,er_Topwidth) / setting%Time%Hydraulics%Dt
                        !% --- select the smaller inflow (negative) magnitude.
                        faceR(thisF,fr_Flowrate) = max(maxQin,faceR(thisF,fr_Flowrate)) 
                    end if

                elseif (faceR(thisF,fr_Flowrate) > zeroR) then
                    !% --- flowrate limiter for outflow (only small/zero depths)    
                    if (elemYN(eup,eYN_isZeroDepth)) then 
                        !% --- no outflow from zero depth
                        faceR(thisF,fr_Flowrate) = zeroR 
                    elseif (elemYN(eup,eYN_isSmallDepth)) then 
                        if (headdif .le. zeroR) then 
                            !% --- head gradient and outflow are inconsistent, so use zero flow
                            faceR(thisF,fr_Flowrate) = zeroR 
                        else
                            !% --- headdif limited by depth
                            headdif = min(headdif, eDepth(eup))
                            !% --- no outflow if less than zerodepth
                            if (headdif .le. setting%ZeroValue%Depth) then 
                                headdif = zeroR
                            end if
                            !% --- max outflow is 1/4 of available volume implied by headdif
                            maxQout  = onefourthR * headdif * elemR(eup,er_Length) &
                                * elemR(eup,er_Topwidth) / setting%Time%Hydraulics%Dt
                            !% --- select smaller outflow (positive) magnitude
                            faceR(thisF,fr_Flowrate) = min(maxQout,faceR(thisF,fr_Flowrate))
                        end if
                    else 
                        !% --- no outflow limiter for deeper then small depths
                    end if 
                else 
                    !% -- no change to zero flowrate
                end if
            end do
        end if

        !% --- set the Preissmann number
        if (.not. isBConly) then 
            do ii=1,nBC
                thisF   => fidx(ii)
                eup     => faceI(thisF,fi_Melem_uL)
                faceR(thisF,fr_Preissmann_Number) = elemR(eup,er_Preissmann_Number)
            end do
        end if

        !% --- set the velocity
        if (.not. isBConly) then 
            do ii=1,nBC
                thisF   => fidx(ii)
                eup     => faceI(thisF,fi_Melem_uL)
                faceR(thisF,fr_Velocity_u) = faceR(thisF,fr_Flowrate) / faceR(thisF,fr_Area_u)
                faceR(thisF,fr_Velocity_d) = faceR(thisF,fr_Flowrate) / faceR(thisF,fr_Area_d)

                !% --- apply limiter
                if (setting%Limiter%Velocity%UseLimitMaxYN) then 
                    if (abs(faceR(thisF,fr_Velocity_u)) > setting%Limiter%Velocity%Maximum) then
                        faceR(thisF,fr_Velocity_u) = sign(0.99d0 * setting%Limiter%Velocity%Maximum, &
                                                     faceR(thisF,fr_Velocity_u))
                    end if
                    if (abs(faceR(thisF,fr_Velocity_d)) > setting%Limiter%Velocity%Maximum) then
                        faceR(thisF,fr_Velocity_d) = sign(0.99d0 * setting%Limiter%Velocity%Maximum, &
                                                     faceR(thisF,fr_Velocity_d))
                    end if
                end if
            end do
        end if

    end subroutine face_interpolation_dnBC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine face_interpolation_interior (facePackCol, Gyn, Hyn, Qyn, skipJump)
        !%------------------------------------------------------------------
        !% Description:
        !% Interpolates all faces using a pack
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: facePackCol  !% Column in faceP array for needed pack
            logical, intent(in) :: Gyn, Hyn, Qyn
            logical, intent(in) :: skipJump !% this is true if jump_compute is skipped 
            integer, pointer    ::  Npack        !% expected number of packed rows in faceP.
            integer :: fGeoSetU(2), fGeoSetD(2), eGeoSet(2)
            integer :: fHeadSetU(1), fHeadSetD(1), eHeadSet(1)
            integer :: fFlowSet(2), eFlowSet(2)

            real(8), pointer :: fQold(:)

            character(64) :: subroutine_name = 'face_interpolation_interior'
        !%------------------------------------------------------------------
        !% Aliases       
            Npack => npack_faceP(facePackCol)
            if (Npack < 1) return
            fQold => faceR(:,fr_Temp01)
        !%------------------------------------------------------------------  
        !% Preliminaries    
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"  
        !%------------------------------------------------------------------     
        !% --- interpolate to ..._u
        !%     identify hydraulic jumps
        !%     set .._u and ..d based on jumps

        if (Gyn) then
            fGeoSetU = [fr_Area_u, fr_Depth_u]
            fGeoSetD = [fr_Area_d, fr_Depth_d]
            eGeoSet  = [er_Area,   er_Depth]

            call face_interp_interior_set &
                (fGeoSetU, eGeoSet,   er_InterpWeight_dG, er_InterpWeight_uG, facePackCol, Npack) 

            !% --- copy upstream to downstream storage at a face
            !%     (only for Head and Geometry types)
            !%     note that these might be reset by hydraulic jump
            call face_copy_upstream_to_downstream_interior &
                (fGeoSetD, fGeoSetU, facePackCol, Npack)
        end if

        if (Hyn) then
            fHeadSetU = [fr_Head_u]
            fHeadSetD = [fr_Head_d]
            eHeadSet  = [er_Head]

            call face_interp_interior_set &
                (fHeadSetU, eHeadSet, er_InterpWeight_dH, er_InterpWeight_uH, facePackCol, Npack)

            !% --- copy upstream to downstream storage at a face
            !%     (only for Head and Geometry types)
            !%     note that these might be reset by hydraulic jump
            call face_copy_upstream_to_downstream_interior &
                (fHeadSetD, fHeadSetU, facePackCol, Npack)
        end if

        if (Qyn) then
            fFlowSet = [fr_Flowrate, fr_Preissmann_Number]
            eFlowSet = [er_Flowrate, er_Preissmann_Number]

            call face_interp_interior_set &
                (fFlowSet, eFlowSet, er_InterpWeight_dQ, er_InterpWeight_uQ, facePackCol, Npack)

            !% --- calculate the velocity in faces and put limiter
            call face_velocities (facePackCol, .true.)
            
            !% --- compute volume-based limits on flowrate
            call face_flowrate_limits_interior (facePackCol)
        end if

        !% --- reset all the hydraulic jump interior faces
        if (.not. skipJump) then
            call jump_compute
        end if

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine face_interpolation_interior
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine face_interpolation_shared (facePackCol, Gyn, Hyn, Qyn, skipJump)
        !%------------------------------------------------------------------
        !% Description:
        !% Interpolates all the shared faces
        !% NOTE -- we do NOT put Npack conditionals on the subroutines that
        !% are called herein so that we can effectively use sync all across
        !% images
        !%-------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: facePackCol  !% Column in faceP array for needed pack
            logical, intent(in) :: Gyn, Hyn, Qyn !% = true for interp geometry, head, flowrate, respectively
            logical, intent(in ):: skipJump
            integer, pointer    :: Npack        !% expected number of packed rows in faceP.
            integer :: fGeoSetU(2), fGeoSetD(2), eGeoSet(2)
            integer :: fHeadSetU(1), fHeadSetD(1), eHeadSet(1)
            integer :: fFlowSet(2), eFlowSet(2)
            integer(kind=8) :: crate, cmax, cval
            character(64) :: subroutine_name = 'face_interpolation_shared'
        !%-------------------------------------------------------------------
        !% Aliases
            Npack => npack_facePS(facePackCol)
        !%-------------------------------------------------------------------
        !% Preliminaries   
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

            !% start the shared timer    
            sync all    
            if (this_image()==1) then
                call system_clock(count=cval,count_rate=crate,count_max=cmax)
                setting%Time%WallClock%SharedStart = cval
                setting%Time%WallClock%SharedStart_C = cval
            end if
        !%-------------------------------------------------------------------

        !% --- General approach
        !%     interpolate to ..._u
        !%     identify hydraulic jumps
        !%     set .._u and ..d based on jumps

        !% --- transfer all the local elemR data needed for face interpolation into elemB data structure
        call local_data_transfer_to_boundary_array (facePackCol, Npack)

        !% --- use elemB to transfer remote data to local elemG array for interpolation
        call inter_image_data_transfer (facePackCol, Npack)

        !% --- set the matching sets
        if (Gyn) then
            fGeoSetU = [fr_Area_u, fr_Depth_u]
            fGeoSetD = [fr_Area_d, fr_Depth_d]
            eGeoSet  = [er_Area, er_Depth]

            call face_interp_shared_set &
                (fGeoSetU, eGeoSet, er_InterpWeight_dG, er_InterpWeight_uG, facePackCol, Npack)

            !% --- copy upstream to downstream storage at a face
            !%     (only for Head and Geometry types)
            !%     note that these might be reset by hydraulic jump
            call face_copy_upstream_to_downstream_shared &
                (fGeoSetD, fGeoSetU, facePackCol, Npack)
        end if

        if (Hyn) then
            fHeadSetU = [fr_Head_u]
            fHeadSetD = [fr_Head_d]
            eHeadSet  = [er_Head]

            call face_interp_shared_set &
                (fHeadSetU, eHeadSet, er_InterpWeight_dH, er_InterpWeight_uH, facePackCol, Npack)

            !% --- copy upstream to downstream storage at a face
            !%     (only for Head and Geometry types)
            !%     note that these might be reset by hydraulic jump
            call face_copy_upstream_to_downstream_shared &
                (fHeadSetD, fHeadSetU, facePackCol, Npack)
        end if

        if (Qyn) then
            fFlowSet = [fr_Flowrate, fr_Preissmann_Number]
            eFlowSet = [er_Flowrate, er_Preissmann_Number]

            call face_interp_shared_set &
                (fFlowSet, eFlowSet, er_InterpWeight_dQ, er_InterpWeight_uQ, facePackCol, Npack)

            call face_velocities (facePackCol, .false.)

            !% --- compute volume-based limits on flowrate
            call face_flowrate_limits_shared (facePackCol)
        end if
            
        !% --- reset all the hydraulic jump interior faces
        if (.not. skipJump) then
            print *, 'ERROR: need to set up a test case for a hydraulic jump across a shared image'
            print *, 'need a version of the jump_compute for shared faces'
            call util_crashpoint(7209885)
        end if


        !%-------------------------------------------------------------------
        !% closing   
            !% stop the shared timer
            sync all
            if (this_image()==1) then
                call system_clock(count=cval,count_rate=crate,count_max=cmax)
                setting%Time%WallClock%SharedStop = cval
                setting%Time%WallClock%SharedCumulative &
                        = setting%Time%WallClock%SharedCumulative &
                        + setting%Time%WallClock%SharedStop &
                        - setting%Time%WallClock%SharedStart

                setting%Time%WallClock%SharedStop_C = cval
                setting%Time%WallClock%SharedCumulative_C &
                        = setting%Time%WallClock%SharedCumulative_C &
                        + setting%Time%WallClock%SharedStop_C &
                        - setting%Time%WallClock%SharedStart_C            
            end if 

            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine face_interpolation_shared
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine face_interp_interior_set &
        (fset, eset, eWdn, eWup, facePackCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% Interpolates to a face for a set of variables 
        !% NOTE cannot sync all in this subroutine
        !%-------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: fset(:), eset(:), eWdn, eWup, facePackCol, Npack
            integer, pointer :: thisP(:), eup(:), edn(:)
            integer :: ii
            character(64) :: subroutine_name = 'face_interp_interior_set'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% Aliases        
            thisP => faceP(1:Npack,facePackCol)
            eup   => faceI(:,fi_Melem_uL)
            edn   => faceI(:,fi_Melem_dL)
        !%--------------------------------------------------------------------
        !% --- cycle interpolation through each type in the set.

        do ii=1,size(fset)

            faceR(thisP,fset(ii)) = &
                (+elemR(eup(thisP),eset(ii)) * elemR(edn(thisP),eWup) &
                 +elemR(edn(thisP),eset(ii)) * elemR(eup(thisP),eWdn) &
                ) / &
                ( elemR(edn(thisP),eWup) + elemR(eup(thisP),eWdn))

        end do

        !% NOTES
        !% elemR(eup(thisP),eset(ii)) is the element value upstream of the face
        !% elemR(edn(thisP),eset(ii) is the element value downstream of the face.
        !% elemR(eup(thisp),eWdn) is the downstream weighting of the upstream element

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine face_interp_interior_set
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine face_interp_shared_set &
        (fset, eset, eWdn, eWup, facePackCol, Npack)
        !%-------------------------------------------------------------------
        !% Description:
        !% Interpolates faces shared between processor
        !% NOTE cannot "sync all" in this subroutine
        !%-------------------------------------------------------------------
        !% Declarations
            integer             :: ii, jj
            integer, intent(in) :: fset(:), eset(:), eWdn, eWup
            integer, intent(in) :: facePackCol, Npack
            integer, pointer    :: thisP, eup, edn, BUpIdx, BDnIdx
            logical, pointer    :: isGhostUp, isGhostDn
            integer(kind=8)     :: crate, cmax, cval
            character(64)       :: subroutine_name = 'face_interp_shared_set'
        !%--------------------------------------------------------------------
        !%  Preliminaries
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

            if (this_image()==1) then
                call system_clock(count=cval,count_rate=crate,count_max=cmax)
                setting%Time%WallClock%SharedStart_A = cval
            end if    
        !%--------------------------------------------------------------------

        !% --- cycle through all the shared faces
        do ii = 1,Npack
            !%-----------------------------------------------------------------
            !% Aliases
                thisP           => facePS(ii,facePackCol)
                eup             => faceI(thisP,fi_Melem_uL)
                edn             => faceI(thisP,fi_Melem_dL)
                BUpIdx          => faceI(thisP,fi_BoundaryElem_uL)
                BDnIdx          => faceI(thisP,fi_BoundaryElem_dL)
                isGhostUp       => faceYN(thisP,fYN_isUpGhost)
                isGhostDn       => faceYN(thisP,fYN_isDnGhost)
            !%-----------------------------------------------------------------
            !% --- cycle through each element in the set.
            !%     This is designed for fset and eset being vectors, but it
            !%     is not clear that this is needed.
            do jj=1,size(fset)

                !% --- condition for upstream element of the shared face is ghost and in a different image
                if (isGhostUp) then
                    faceR(thisP,fset(jj)) = &
                        (+elemGR(ii,eset(jj))  * elemB%R(ii,eWup)  &
                         +elemB%R(ii,eset(jj)) * elemGR(ii,eWdn) &
                        ) / &
                        ( elemB%R(ii,eWup) + elemGR(ii,eWdn) )

                !% --- condition for downstream element of the shared face is ghost and in a different image
                elseif (isGhostDn) then

                    faceR(thisP,fset(jj)) = &
                        (+elemB%R(ii,eset(jj)) * elemGR(ii,eWup) &
                         +elemGR(ii,eset(jj))  * elemB%R(ii,eWdn)  &
                        ) / &
                        ( elemGR(ii,eWup) + elemB%R(ii,eWdn) )
                else
                    write(*,*) 'CODE ERROR: unexpected else'
                    call util_crashpoint( 487874)
                end if      
            end do
        end do

        !% NOTES
        !% elemB%R(ii,eset(jj)) is the element value of the boundary element
        !% elemGR(ii,eset(jj)) is the element value of the ghost element
        !% elemB%R(ii,eWdn) is the downstream weighting of the boundary element
        !% elemGR(ii,eWup)) is the upstream weighting of the ghost element

        !%--------------------------------------------------------------------
        !% Closing
            sync all
            if (this_image()==1) then
                !% ---stop the shared timer
                call system_clock(count=cval,count_rate=crate,count_max=cmax)
                setting%Time%WallClock%SharedStop_A = cval
                setting%Time%WallClock%SharedCumulative_A &
                        = setting%Time%WallClock%SharedCumulative_A &
                        + setting%Time%WallClock%SharedStop_A &
                        - setting%Time%WallClock%SharedStart_A                    
            end if 
                if (setting%Debug%File%face) &
                    write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine face_interp_shared_set
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine local_data_transfer_to_boundary_array &
        (facePackCol, Npack)
        !%-------------------------------------------------------------------
        !% Description:
        !% transfers local data from elemR to elemB%R
        !%-------------------------------------------------------------------
        !% Declarations
            integer             :: ii !, eColumns(Ncol_elemBGR) 
            integer, intent(in) :: facePackCol, Npack
            integer, pointer    :: thisP, eUp, eDn, JMidx
            logical, pointer    :: isGhostUp, isGhostDn
            character(64)       :: subroutine_name = 'local_data_transfer_to_boundary_array'
        !%--------------------------------------------------------------------
        !%  Preliminaries
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"  
        !%--------------------------------------------------------------------
        !% cycle through all the shared faces
        sync all
        do ii = 1,Npack
            
            !%-----------------------------------------------------------------
            !% Aliases
                thisP      => facePS(ii,facePackCol)
                isGhostUp  => faceYN(thisP,fYN_isUpGhost)
                isGhostDn  => faceYN(thisP,fYN_isDnGhost)
                eUp        => faceI(thisP,fi_Melem_uL)
                eDn        => faceI(thisP,fi_Melem_dL)
            !%-----------------------------------------------------------------

            !% --- condition for upstream element is ghost
            if (isGhostUp) then
                elemB%R(ii,:) = elemR(eDn,:)
            !% --- condition for downstream element is ghost
            elseif (isGhostDn) then
                elemB%R(ii,:) = elemR(eUp,:)
            else
                write(*,*) 'CODE ERROR: unexpected else'
                call util_crashpoint( 487874)
            end if     
            !% --- handle special case for volume used by Pump Type 1 when
            !%     the upstream element is a JB
            if ((isGhostUp) .and. (elemI(eUp,ei_elementType) == JB)) then
                JMidx => elemSI(eup,esi_JunctionBranch_Main_Index)
                elemB%R(ii,er_Volume) = elemR(JMidx,er_Volume)
            end if
        end do

        !%--------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine local_data_transfer_to_boundary_array
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine inter_image_data_transfer &
        (facePackCol, Npack)
        !%-------------------------------------------------------------------
        !% Description:
        !% transfers data from connected images
        !%-------------------------------------------------------------------
        !% Declarations
            integer             :: ii  
            integer, intent(in) :: facePackCol, Npack
            integer, pointer    :: thisP, ci, BUpIdx, BDnIdx
            logical, pointer    :: isGhostUp, isGhostDn
            integer(kind=8)     :: crate, cmax, cval
            character(64)       :: subroutine_name = 'inter_image_data_transfer'
        !%--------------------------------------------------------------------
        !%  Preliminaries
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

            sync all
            if (this_image()==1) then
                call system_clock(count=cval,count_rate=crate,count_max=cmax)
                setting%Time%WallClock%SharedStart_A = cval
            end if 
          
        !%--------------------------------------------------------------------

        !% ---cycle through all the shared faces
        do ii = 1,Npack
            !%-----------------------------------------------------------------
            !% Aliases
                thisP      => facePS(ii,facePackCol)
                ci         => faceI(thisP,fi_Connected_image)
                BUpIdx     => faceI(thisP,fi_BoundaryElem_uL)
                BDnIdx     => faceI(thisP,fi_BoundaryElem_dL)
                isGhostUp  => faceYN(thisP,fYN_isUpGhost)
                isGhostDn  => faceYN(thisP,fYN_isDnGhost)
            !%-----------------------------------------------------------------
            !% --- condition for upstream element of the shared face is ghost and in a different image
            if (isGhostUp) then
                elemGR(ii,:) = elemB[ci]%R(BUpIdx,:) 

            !% ---condition for downstream element of the shared face is ghost and in a different image
            elseif (isGhostDn) then
                elemGR(ii,:) = elemB[ci]%R(BDnIdx,:)

            else
                write(*,*) 'CODE ERROR: unexpected else'
                call util_crashpoint(487874)

            end if        
        end do
        
        !%--------------------------------------------------------------------
        !% Closing
            sync all
            
            if (this_image()==1) then
                !% stop the shared timer
                call system_clock(count=cval,count_rate=crate,count_max=cmax)
                setting%Time%WallClock%SharedStop_A = cval
                setting%Time%WallClock%SharedCumulative_A &
                        = setting%Time%WallClock%SharedCumulative_A &
                        + setting%Time%WallClock%SharedStop_A &
                        - setting%Time%WallClock%SharedStart_A                    
            end if 
                if (setting%Debug%File%face) &
                    write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
                    
    end subroutine inter_image_data_transfer
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine face_copy_upstream_to_downstream_interior &
        (downstreamSet, upstreamSet, facePackCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% Copies the interpolated value on the upstream side to the downstream side
        !% These values are later adjusted for hydraulic jumps
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: facePackCol, Npack, downstreamSet(:), upstreamSet(:)
            integer, pointer :: thisP(:)
            character(64) :: subroutine_name = 'face_copy_upstream_to_downstream_interior'
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------

        if (Npack > 0) then
            thisP => faceP(1:Npack,facePackCol)
            faceR(thisP,downstreamSet) = faceR(thisP,upstreamSet)
        end if

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine face_copy_upstream_to_downstream_interior
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine face_copy_upstream_to_downstream_shared &
        (downstreamSet, upstreamSet, facePackCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% Copies the interpolated value on the upstream side to the downstream side
        !% These values are later adjusted for hydraulic jumps
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: facePackCol, Npack, downstreamSet(:), upstreamSet(:)
            integer, pointer :: thisP(:)
            integer(kind=8) :: crate, cmax, cval
            character(64) :: subroutine_name = 'face_copy_upstream_to_downstream'
        !%-------------------------------------------------------------------
        !% Preliminaries
            sync all
            if (this_image()==1) then
                call system_clock(count=cval,count_rate=crate,count_max=cmax)
                setting%Time%WallClock%SharedStart_B = cval
            end if 
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------

        if (Npack > 0) then
            thisP => facePS(1:Npack,facePackCol)
            faceR(thisP,downstreamSet) = faceR(thisP,upstreamSet)
        end if

        !%-------------------------------------------------------------------
        !% Closing
            sync all
            if (this_image()==1) then
                !% stop the shared timer
                call system_clock(count=cval,count_rate=crate,count_max=cmax)
                setting%Time%WallClock%SharedStop_B = cval
                setting%Time%WallClock%SharedCumulative_B &
                        = setting%Time%WallClock%SharedCumulative_B &
                        + setting%Time%WallClock%SharedStop_B &
                        - setting%Time%WallClock%SharedStart_B                    
            end if 
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine face_copy_upstream_to_downstream_shared
!%
!%==========================================================================
!%==========================================================================
!%   
    subroutine face_velocities (facePackCol, isInterior)
        !%------------------------------------------------------------------
        !% Description:
        !% This subroutine calculates the face valocity and adjusts for limiter
        !% NOTE: this subroutine CANNOT address adjacent element data. The
        !% facePackCol is allowed to include fp_..._all data, which contains
        !% both interior and shared faces.
        !%-------------------------------------------------------------------  
            integer, intent(in) :: facePackCol
            logical, intent(in) :: isInterior
            integer, pointer :: Npack, thisP(:)
            real(8), pointer :: f_area_u(:), f_area_d(:), f_velocity_u(:), f_velocity_d(:)
            real(8), pointer :: f_flowrate(:), zeroValue, vMax
            character(64) :: subroutine_name = 'adjust_face_dynamic_limit'
        !%-------------------------------------------------------------------
            if (setting%Debug%File%adjust) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% Aliases
            f_area_u     => faceR(:,fr_Area_u)
            f_area_d     => faceR(:,fr_Area_d)
            f_velocity_u => faceR(:,fr_Velocity_u)
            f_velocity_d => faceR(:,fr_Velocity_d)
            f_flowrate   => faceR(:,fr_Flowrate)
            zeroValue    => setting%ZeroValue%Area
            vMax         => setting%Limiter%Velocity%Maximum
        !%----------------------------------------------------------------------
        if (isInterior) then
            !% face velocity calculation at the interior faces
            Npack => npack_faceP(facePackCol)
            thisP => faceP(1:Npack,facePackCol)
        else
            !% face velocity calculation at the shared faces
            Npack => npack_facePS(facePackCol)
            thisP => facePS(1:Npack,facePackCol)
        end if

        if (Npack > 0) then
   
            !% ensure face area_u is not smaller than zerovalue
            where (f_area_u(thisP) < zeroValue)
                f_area_u(thisP) = zeroValue
            endwhere
            !% ensure face area_d is not smaller than zerovalue
            where (f_area_d(thisP) < zeroValue)
                f_area_d(thisP) = zeroValue
            endwhere

            f_velocity_u(thisP) = f_flowrate(thisP)/f_area_u(thisP)

            f_velocity_d(thisP) = f_flowrate(thisP)/f_area_d(thisP)

            !%  limit high velocities
            if (setting%Limiter%Velocity%UseLimitMaxYN) then
                where(abs(f_velocity_u(thisP))  > vMax)
                    f_velocity_u(thisP) = sign(0.99d0 * vMax, f_velocity_u(thisP))
                endwhere

                where(abs(f_velocity_d(thisP))  > vMax)
                    f_velocity_d(thisP) = sign(0.99d0 * vMax, f_velocity_d(thisP))
                endwhere
            end if

        end if    
        
        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%adjust) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine face_velocities
!%  
!%========================================================================== 
!%==========================================================================
!%
    subroutine face_flowrate_limits_interior (facePackCol)   
        !%------------------------------------------------------------------
        !% Description:
        !% stores maximum and minimum flowrates on the face based on
        !% emptying the volume of the upstream element (maximum) or the 
        !% downstream element (reversed flow maximum negative or minimum flow)
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: facePackCol
            integer, pointer :: Npack
            integer, pointer :: thisP(:), eup(:), edn(:)
            integer, pointer :: idx_fBCdn(:), idx_fBCup(:)
            real(8), pointer :: dt, eVolume(:)
        !%------------------------------------------------------------------
        !% Aliases
            Npack => npack_faceP(facePackCol)
            if (Npack < 1) return
            thisP   => faceP(1:Npack,facePackCol)
            eup     => faceI(:,fi_Melem_uL)
            edn     => faceI(:,fi_Melem_dL)
            eVolume => elemR(:,er_Volume)
            dt      => setting%Time%Hydraulics%Dt

            !% --- downstream BC on faces
            idx_fBCdn       => faceP(1:npack_faceP(fp_BCdn),fp_BCdn)
            !% --- upstream BC on faces
            idx_fBCup       => faceP(1:npack_faceP(fp_BCup),fp_BCup)
        !%------------------------------------------------------------------ 

        faceR(thisP,fr_FlowrateMax) =  eVolume(eup(thisP)) / dt
        faceR(thisP,fr_FlowrateMin) = -eVolume(edn(thisP)) / dt

        !% --- set downstream BC to allow any level of inflow
        faceR(idx_fBCdn,fr_FlowrateMin) = -nullvalueR

        !% --- set upstream BC face to allow any level of inflow
        faceR(idx_fBCup,fr_FlowrateMax) = nullvalueR

    end subroutine face_flowrate_limits_interior
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine face_flowrate_limits_shared (facePackCol)
        !%-------------------------------------------------------------------
        !% Description:
        !% stores maximum and minimum flowrates on the face based on
        !% emptying the volume of the upstream element (maximum) or the 
        !% downstream element (reversed flow maximum negative or minimum flow)
        !% Must be conducted AFTER  inter_image_data_transfer
        !%-------------------------------------------------------------------
        !% Declarations
            integer             :: ii
            integer, intent(in) :: facePackCol  !% Column in faceP array for needed pack
            integer, pointer    :: Npack        !% expected number of packed rows in faceP.
            integer, pointer    :: thisP, eup, edn, BUpIdx, BDnIdx
            real(8), pointer    :: dt
            logical, pointer    :: isGhostUp, isGhostDn
        !%-------------------------------------------------------------------
        !% Preliminaries   
            Npack => npack_facePS(facePackCol)
            if (Npack < 1) return
            dt => setting%Time%Hydraulics%Dt
    
        !%-------------------------------------------------------------------   
        do ii=1,Npack
            !%---------------------------------------------------------------
            !% Local Aliases
                thisP           => facePS(ii,facePackCol)
                eup             => faceI(thisP,fi_Melem_uL)
                edn             => faceI(thisP,fi_Melem_dL)
                BUpIdx          => faceI(thisP,fi_BoundaryElem_uL)
                BDnIdx          => faceI(thisP,fi_BoundaryElem_dL)
                isGhostUp       => faceYN(thisP,fYN_isUpGhost)
                isGhostDn       => faceYN(thisP,fYN_isDnGhost)
            !%---------------------------------------------------------------

            if (isGhostUp) then 
                faceR(thisP,fr_FlowrateMax) =  elemGR(ii,er_Volume) / dt
            elseif (isGhostDn) then 
                faceR(thisP,fr_FlowrateMin) = -elemGR(ii,er_Volume) / dt 
            end if

        end do


    end subroutine face_flowrate_limits_shared
!%
!%==========================================================================
!%==========================================================================
!% 
    subroutine face_force_JBvalues (frCol, erCol, fiIdx, JMidx, kstart)
        !%------------------------------------------------------------------
        !% Description:
        !% Forces the JB element value to the adjacent face
        !%------------------------------------------------------------------
            integer, intent(in) :: frCol  !% column in faceR array for output
            integer, intent(in) :: erCol  !% column in elemR array for input
            integer, intent(in) :: fiIdx  !% face index column for up/dn map
            integer, intent(in) :: JMidx  !% junction main index
            integer, intent(in) :: kstart !% =1 for upstream, 2 for downstream
            integer :: kk
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        do kk=kstart,max_branch_per_node,2
            if (elemSI(JMidx+kk,esi_JunctionBranch_Exists).ne. oneI) cycle

            faceR(elemI(JMidx+kk,fiIdx),frCol) = elemR(JMidx+kk,erCol)
            !% --- set the shared face diverge to true
            if (faceYN(elemI(JMidx+kk,fiIdx),fYN_isSharedFace)) then
                faceYN(elemI(JMidx+kk,fiIdx),fYN_isSharedFaceDiverged) = .true.
            end if
        end do
    
    end subroutine face_force_JBvalues
!%    
!%==========================================================================  
!%==========================================================================
!%   
    subroutine face_zeroDepth_geometry_interior (facePcol)
        !% -----------------------------------------------------------------
        !% Description:
        !% sets the face geoemtry values for interior faces with neighbor elements
        !% that are zero depth
        !% --- For CC
        !%  facePcol must be one of fp_CC_downstream_is_zero_IorS,
        !%  fp_CC_upstream_is_zero_IorS, fp_CC_bothsides_are_zero_IorS
        !% --- For JB
        !%  facePcol must be one of fp_JB_downstream_is_zero_IorS,
        !%  fp_JB_upstream_is_zero_IorS, fp_JB_bothsides_are_zero_IorS
        !% -----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: facePcol
            integer, pointer :: npack, thisP(:), edn(:), eup(:)
        !% -----------------------------------------------------------------
        !% Preliminaries
        !% -----------------------------------------------------------------
        npack => npack_faceP(facePCol)
        if (npack < 1) return
        
        thisP => faceP(1:npack,facePcol)
        edn   => faceI(:,fi_Melem_dL)
        eup   => faceI(:,fi_Melem_uL)

        select case (facePcol)
            case (fp_CC_downstream_is_zero_IorS, fp_JB_downstream_is_zero_IorS)

                !% ---- head is the smaller of the recent face value, the upstream element value,
                !%      or the depth upstream applied to the face zbottom (new 20230430brh)
                faceR(thisP,fr_Head_u) = min(faceR(thisP,fr_Head_u), elemR(eup(thisP),er_Head), elemR(eup(thisP),er_Depth) + faceR(thisP,fr_Zbottom))
                !% --- store on downstream face
                faceR(thisP,fr_Head_d) = faceR(thisP,fr_Head_u)

                !% --- depth is the computed depth from head or zeroDepth
                faceR(thisP,fr_Depth_u) = max(faceR(thisP,fr_Head_u) - faceR(thisP,fr_Zbottom), 0.99d0 * setting%ZeroValue%Depth)
                faceR(thisP,fr_Depth_d) = faceR(thisP,fr_Depth_u)

                !% --- area is larger of the existing or the upstream
                faceR(thisP,fr_Area_u) = max(faceR(thisP,fr_Area_u), elemR(eup(thisP),er_Area))
                faceR(thisP,fr_Head_d) = faceR(thisP,fr_Head_u)

                !% NOTE: the "consistent" approach does not work because at small depths the
                !% computation of area from depth is not consistent with the depth obtained from
                !% volume

            case (fp_CC_upstream_is_zero_IorS,fp_JB_upstream_is_zero_IorS)
                !% --- head is the smaller of the recent value or the downstream element value
                !%     or the depth upstream applied to the face zbottom (new 20230430brh)
                faceR(thisP,fr_Head_d) = min(faceR(thisP,fr_Head_d), elemR(edn(thisP),er_Head), elemR(edn(thisP),er_Depth)+ faceR(thisP,fr_Zbottom))
                faceR(thisP,fr_Head_u) = faceR(thisP,fr_Head_d)

                !% --- depth is the computed depth from head or zeroDepth
                faceR(thisP,fr_Depth_d) = max(faceR(thisP,fr_Head_d) - faceR(thisP,fr_Zbottom), 0.99d0 * setting%ZeroValue%Depth)
                faceR(thisP,fr_Depth_u) = faceR(thisP,fr_Depth_d)

                !% --- area is larger of the existing or the downstream
                faceR(thisP,fr_Area_d) = max(faceR(thisP,fr_Area_d), elemR(edn(thisP),er_Area))
                faceR(thisP,fr_Head_u) = faceR(thisP,fr_Head_d)

                !% NOTE: the "consistent" approach does not work because at small depths the
                !% computation of area from depth is not consistent with the depth obtained from
                !% volume

            case (fp_CC_bothsides_are_zero_IorS,fp_JB_bothsides_are_zero_IorS)
                faceR(thisP,fr_Depth_u) = 0.99d0 * setting%ZeroValue%Depth
                faceR(thisP,fr_Depth_d) = 0.99d0 * setting%ZeroValue%Depth
                faceR(thisP,fr_Head_u) = faceR(thisP,fr_Depth_u) + faceR(thisP,fr_Zbottom)
                faceR(thisP,fr_Head_d) = faceR(thisP,fr_Depth_d) + faceR(thisP,fr_Zbottom)
                faceR(thisP,fr_Area_u) = setting%ZeroValue%Area
                faceR(thisP,fr_Area_d) = setting%ZeroValue%Area
                
            case default
                print *, 'CODE ERROR: unexpected case default'
                call util_crashpoint(619873)
        end select

    end subroutine face_zeroDepth_geometry_interior
!%  
!%========================================================================== 
!%==========================================================================
!%   
    subroutine face_zeroDepth_geometry_shared (facePcol)
        !% -----------------------------------------------------------------
        !% Description:
        !% sets the face geoemtry values for shared faces with neighbor elements
        !% that are zero depth 
        !% --- For CC
        !%  facePcol must be one of fp_CC_downstream_is_zero_IorS,
        !%  fp_CC_upstream_is_zero_IorS, fp_CC_bothsides_are_zero_IorS
        !% --- For JB
        !%  facePcol must be one of fp_JB_downstream_is_zero_IorS,
        !%  fp_JB_upstream_is_zero_IorS, fp_JB_bothsides_are_zero_IorS
        !% -----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: facePcol
            integer, pointer :: npack, thisP, edn, eup, GUp, GDn, ci, Ifidx
            logical, pointer :: isGhostUp, isGhostDn
            integer :: ii
        !% -----------------------------------------------------------------
        !% Preliminaries
            npack => npack_facePS(facePCol)
            if (npack < 1) return
        !% -----------------------------------------------------------------

        do ii = 1,Npack
            !%-----------------------------------------------------------------
            !% Aliases
                thisP           => facePS(ii,facePCol)
                ci              => faceI(thisP,fi_Connected_image)
                eup             => faceI(thisP,fi_Melem_uL)
                edn             => faceI(thisP,fi_Melem_dL)
                GUp             => faceI(thisP,fi_GhostElem_uL)
                GDn             => faceI(thisP,fi_GhostElem_dL)
                Ifidx           => faceI(thisP,fi_Identical_Lidx)
                isGhostUp       => faceYN(thisP,fYN_isUpGhost)
                isGhostDn       => faceYN(thisP,fYN_isDnGhost)
            
            select case (facePcol)
                case (fp_CC_downstream_is_zero_IorS, fp_JB_downstream_is_zero_IorS)

                    if (.not. isGhostUp) then
                        !% --- head is the smaller of the recent value or the downstream element value
                        !%     or the depth upstream applied to the face zbottom
                        faceR(thisP,fr_Head_u) = min(faceR(thisP,fr_Head_u), elemR(eup,er_Head), elemR(eup,er_Depth) + faceR(thisP,fr_Zbottom)) 

                        faceR(thisP,fr_Head_d) = faceR(thisP,fr_Head_u)
                        faceR(thisP,fr_Depth_u) = max(faceR(thisP,fr_Head_u) - faceR(thisP,fr_Zbottom), 0.99d0 * setting%ZeroValue%Depth)
                        faceR(thisP,fr_Depth_d) = faceR(thisP,fr_Depth_u)

                        ! if (faceR(thisP,fr_Depth_u) > setting%ZeroValue%Depth) then 
                        !     faceR(thisP,fr_Area_u) = geo_area_from_depth_singular(eup,faceR(thisP,fr_Depth_u),setting%ZeroValue%Area)
                        ! else
                        !     faceR(thisP,fr_Area_u) = setting%ZeroValue%Area
                        ! end if
                        faceR(thisP,fr_Area_u) = max(faceR(thisP,fr_Area_u), elemR(eup,er_Area))
                        faceR(thisP,fr_Area_d) = faceR(thisP,fr_Area_u) 

                        !% --- the face values should be identical apart from the newly adjusted values
                        !%     transfer the whole data column to the indetical array 
                        faceR(Ifidx,:)[ci] = faceR(thisP,:)
                    
                    else 
                        !% --- no action
                    end if

                case (fp_CC_upstream_is_zero_IorS,fp_JB_upstream_is_zero_IorS)
                    if (.not. isGhostDn) then
                        !% --- head is the smaller of the recent value or the downstream element value
                        !%     or the depth upstream applied to the face zbottom (new 20230430brh)
                        faceR(thisP,fr_Head_d) = min(faceR(thisP,fr_Head_d), elemR(edn,er_Head), elemR(edn,er_Depth)+ faceR(thisP,fr_Zbottom))
                        faceR(thisP,fr_Head_u) = faceR(thisP,fr_Head_d)

                        !% --- depth is the computed depth from head or zeroDepth
                        faceR(thisP,fr_Depth_d) = max(faceR(thisP,fr_Head_d) - faceR(thisP,fr_Zbottom), 0.99d0 * setting%ZeroValue%Depth)
                        faceR(thisP,fr_Depth_u) = faceR(thisP,fr_Depth_d)
                    
                        ! if (faceR(thisP,fr_Depth_d) > setting%ZeroValue%Depth) then 
                        !     faceR(thisP,fr_Area_d) = geo_area_from_depth_singular(edn,faceR(thisP,fr_Depth_d),setting%ZeroValue%Area)
                        ! else
                        !     faceR(thisP,fr_Area_d) = setting%ZeroValue%Area
                        ! end if
                        faceR(thisP,fr_Area_d) = max(faceR(thisP,fr_Area_d), elemR(edn,er_Area))
                        faceR(thisP,fr_Area_u) = faceR(thisP,fr_Area_d)

                        !% --- the face values should be identical apart from the newly adjusted values
                        !%     transfer the whole data column to the indetical array 
                        faceR(Ifidx,:)[ci] = faceR(thisP,:)
                    else 
                        !% --- no action
                    end if
                
                case (fp_CC_bothsides_are_zero_IorS,fp_JB_bothsides_are_zero_IorS)
                    faceR(thisP,fr_Depth_u) = 0.99d0 * setting%ZeroValue%Depth
                    faceR(thisP,fr_Depth_d) = 0.99d0 * setting%ZeroValue%Depth
                    faceR(thisP,fr_Head_u) = faceR(thisP,fr_Depth_u) + faceR(thisP,fr_Zbottom)
                    faceR(thisP,fr_Head_d) = faceR(thisP,fr_Depth_d) + faceR(thisP,fr_Zbottom)
                    faceR(thisP,fr_Area_u) = setting%ZeroValue%Area
                    faceR(thisP,fr_Area_d) = setting%ZeroValue%Area

                case default
                    print *, 'CODE ERROR: unexpected case default'
                    call util_crashpoint(2548712)

            end select
        
        end do 

    end subroutine face_zeroDepth_geometry_shared
!%  
!%========================================================================== 
!%==========================================================================
!%
    subroutine face_zeroDepth_flowrates_interior (facePcol)
        !% -----------------------------------------------------------------
        !% Description:
        !% sets the face flux values for interior faces with neighbor elements
        !% that are zero depth 
        !% --- For CC
        !%  facePcol must be one of fp_CC_downstream_is_zero_IorS,
        !%  fp_CC_upstream_is_zero_IorS, fp_CC_bothsides_are_zero_IorS
        !% --- For JB
        !%  facePcol must be one of fp_JB_downstream_is_zero_IorS,
        !%  fp_JB_upstream_is_zero_IorS, fp_JB_bothsides_are_zero_IorS
        !% -----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: facePcol
            integer, pointer :: npack, thisP(:), edn(:), eup(:)
        !% -----------------------------------------------------------------
        !% Preliminaries
            npack => npack_faceP(facePCol)
            if (npack < 1) return
        !% -----------------------------------------------------------------
        !% Aliases
            thisP => faceP(1:npack,facePcol)
            edn   => faceI(:,fi_Melem_dL)
            eup   => faceI(:,fi_Melem_uL)
        !% -----------------------------------------------------------------

        select case (facePcol)
            case (fp_CC_downstream_is_zero_IorS,fp_JB_downstream_is_zero_IorS)
                !% --- use the downstream flow from the upstream element (or zero if upstream flow)
                !%     NOTE: cutting out limit on Head>Zbottom as it caused issues with Pump that has zero head
                !%     If this revision causes other problems, we may need to have a distinction for Diag elements upstream
                where (elemR(eup(thisP),er_Flowrate) .ge. zeroR)
                    faceR(thisP,fr_Flowrate) = elemR(eup(thisP),er_Flowrate)
                elsewhere
                    faceR(thisP,fr_Flowrate) = zeroR
                endwhere

            case (fp_CC_upstream_is_zero_IorS,fp_JB_upstream_is_zero_IorS)
                !% --- use the upstream flow from the downstream element (or zero if downstream flow)
                where (elemR(edn(thisP),er_Flowrate) .le. zeroR) 
                    faceR(thisP,fr_Flowrate) = elemR(edn(thisP),er_Flowrate)
                elsewhere
                    faceR(thisP,fr_Flowrate) = zeroR
                endwhere

            case (fp_CC_bothsides_are_zero_IorS,fp_JB_bothsides_are_zero_IorS)
                faceR(thisP,fr_Flowrate)   = zeroR
                faceR(thisP,fr_Velocity_d) = zeroR 
                faceR(thisP,fr_Velocity_u) = zeroR

            case default
                print *, 'CODE ERROR: unexpected case default'
                call util_crashpoint(619873)
        end select

        !% ---reset velocities
        !%    HACK: it would be better if we could enusre that all the
        !%    face areas are greater than zero so we wouldn't need the where
        !%    statements. But as of 20230114 there were areas that lead
        !%    to NAN values
        where (faceR(thisP,fr_Area_d) > setting%ZeroValue%Area)
            faceR(thisP,fr_Velocity_d) = faceR(thisP,fr_Flowrate) /  faceR(thisP,fr_Area_d)
        elsewhere
            faceR(thisP,fr_Velocity_d) = zeroR
        endwhere

        where (faceR(thisP,fr_Area_u) > setting%ZeroValue%Area)
            faceR(thisP,fr_Velocity_u) = faceR(thisP,fr_Flowrate) /  faceR(thisP,fr_Area_u)
        elsewhere
            faceR(thisP,fr_Velocity_u) = zeroR
        endwhere

    end subroutine face_zeroDepth_flowrates_interior
!%  
!%========================================================================== 
!%==========================================================================
!%    
    subroutine face_zeroDepth_flowrates_shared (facePcol)
        !% -----------------------------------------------------------------
        !% Description:
        !% sets the face flux values for shared faces with neighbor elements
        !% that are zero depth 
        !% --- For CC
        !%  facePcol must be one of fp_CC_downstream_is_zero_IorS,
        !%  fp_CC_upstream_is_zero_IorS, fp_CC_bothsides_are_zero_IorS
        !% --- For JB
        !%  facePcol must be one of fp_JB_downstream_is_zero_IorS,
        !%  fp_JB_upstream_is_zero_IorS, fp_JB_bothsides_are_zero_IorS
        !% -----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: facePcol
            integer, pointer :: npack, thisP, edn, eup, gUp, gDn, ci, Ifidx
            logical, pointer :: isGhostUp, isGhostDn
            integer :: ii
        !% -----------------------------------------------------------------
            npack => npack_facePS(facePCol)
            if (npack < 1) return
        
        do ii = 1,Npack
            !%-----------------------------------------------------------------
            !% Aliases
                thisP           => facePS(ii,facePCol)
                ci              => faceI(thisP,fi_Connected_image)
                eup             => faceI(thisP,fi_Melem_uL)
                edn             => faceI(thisP,fi_Melem_dL)
                gUp             => faceI(thisP,fi_GhostElem_uL)
                gDn             => faceI(thisP,fi_GhostElem_dL)
                Ifidx           => faceI(thisP,fi_Identical_Lidx)
                isGhostUp       => faceYN(thisP,fYN_isUpGhost)
                isGhostDn       => faceYN(thisP,fYN_isDnGhost)
            !% -----------------------------------------------------------------
            select case (facePcol)

                case (fp_CC_downstream_is_zero_IorS,fp_JB_downstream_is_zero_IorS)

                    if (.not. isGhostUp) then

                        if (elemR(eup,er_Flowrate) .ge. zeroR) then
                            faceR(thisP,fr_Flowrate) = elemR(eup,er_Flowrate)
                        else 
                            faceR(thisP,fr_Flowrate) = zeroR
                        end if
                    else 
                        !% --- no action
                    end if

                case (fp_CC_upstream_is_zero_IorS,fp_JB_upstream_is_zero_IorS)
                    
                    if (.not. isGhostDn) then

                        if (elemR(edn,er_Flowrate) .ge. zeroR) then
                            faceR(thisP,fr_Flowrate) = elemR(edn,er_Flowrate)
                        else 
                            faceR(thisP,fr_Flowrate) = zeroR
                        end if
                    else 
                        !% --- no action
                    end if

                case (fp_CC_bothsides_are_zero_IorS,fp_JB_bothsides_are_zero_IorS)
                    faceR(thisP,fr_Flowrate)   = zeroR
                    faceR(thisP,fr_Velocity_d) = zeroR 
                    faceR(thisP,fr_Velocity_u) = zeroR

                case default
                    print *, 'CODE ERROR: unexpected case default'
                    call util_crashpoint(6198732)
            end select        
        end do

    end subroutine face_zeroDepth_flowrates_shared
!%  
!%========================================================================== 
!%==========================================================================
!%
    subroutine face_shared_face_sync (facePcol)
        !% -----------------------------------------------------------------
        !% Description:
        !% sync data between shared faces
        !% 
        !% -----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: facePcol
            integer, pointer :: npack, thisP, ci, Ifidx
            logical, pointer :: isDiverged
            integer :: ii
        !% -----------------------------------------------------------------
        npack => npack_facePS(facePCol)
        if (npack < 1) return

        do ii = 1,Npack
            !%----------------------------------------------------------
            !% Aliases
                thisP           => facePS(ii,facePCol)
                ci              => faceI(thisP,fi_Connected_image)
                Ifidx           => faceI(thisP,fi_Identical_Lidx)
                isDiverged      => faceYN(thisP,fYN_isSharedFaceDiverged)
            !%----------------------------------------------------------
                
            !% --- check if the shared face has been diverged
            if (isDiverged) then
                !% --- if the face has been diverged, copy the whole shared face 
                !%     column in the identical face location at the connected image
                faceR(Ifidx,:)[ci] = faceR(thisP,:)

                !% --- after the transfer, set the divergence check to false
                isDiverged = .false.
            end if
        end do

    end subroutine face_shared_face_sync
!%  
!%========================================================================== 
!%==========================================================================  
!%  
    subroutine face_shared_face_sync_single (facePcol,col)
        !% -----------------------------------------------------------------
        !% Description:
        !% sync data between shared faces
        !% 
        !% -----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: facePcol,col
            integer, pointer :: npack, thisP, ci, Ifidx
            logical, pointer :: isDiverged
            integer :: ii
        !% -----------------------------------------------------------------
        npack => npack_facePS(facePCol)
        if (npack < 1) return

        do ii = 1,Npack
            !%----------------------------------------------------------
            !% Aliases
                thisP           => facePS(ii,facePCol)
                ci              => faceI(thisP,fi_Connected_image)
                Ifidx           => faceI(thisP,fi_Identical_Lidx)
                isDiverged      => faceYN(thisP,fYN_isSharedFaceDiverged)
            !%----------------------------------------------------------
                
            !% --- check if the shared face has been diverged
            if (isDiverged) then
                !% --- if the face has been diverged, copy the whole shared face 
                !%     column in the identical face location at the connected image
                faceR(Ifidx,col)[ci] = faceR(thisP,col)

                !% --- after the transfer, set the divergence check to false
                isDiverged = .false.
            end if
        end do

    end subroutine face_shared_face_sync_single
!%  
!%========================================================================== 
!%==========================================================================  
!% 
    subroutine face_flowrate_for_openclosed_elem (elemPcol)
        !%------------------------------------------------------------------
        !% Description:
        !% Sets all zero fluxes on downstream faces of "closed" CC elements
        !% i.e. where elemR(:,er_Setting) = 0.0
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: elemPcol
            integer, pointer :: npack, thisP(:), fdn(:)
        !%------------------------------------------------------------------  
        npack => npack_elemP(elemPcol)
        if (npack < 1) return 
        thisP => elemP(1:npack,elemPcol)
        fdn   => elemI(:,ei_Mface_dL)

        where (elemR(thisP,er_Setting) == zeroR)
            faceR(fdn(thisP),fr_Flowrate) = zeroR
            faceR(fdn(thisP),fr_Velocity_d) = zeroR 
            faceR(fdn(thisP),fr_Velocity_u) = zeroR

            where (faceYN(fdn(thisP),fYN_isSharedFace))
                faceYN(fdn(thisP),fYN_isSharedFaceDiverged) = .true.
            end where
        endwhere
   
       end subroutine face_flowrate_for_openclosed_elem    
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module face
