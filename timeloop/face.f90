module face

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use adjust
    use geometry
    use jump
    use pack_mask_arrays, only: pack_CC_zeroDepth_interior_faces, pack_CC_zeroDepth_shared_faces
    use utility_profiler
    use utility, only: util_sign_with_ones, util_syncwrite
    use utility_crash, only: util_crashpoint
   !use utility_unit_testing, only: util_utest_CLprint


    implicit none

    !%----------------------------------------------------------------------
    !% Description:
    !% Provides computation of face values for timeloop of hydraulics
    !%
    !% Methods:
    !% The faces values are determined by interpolation.
    !%
    !%----------------------------------------------------------------------
    private

    public :: face_pull_facedata_to_JBelem
    public :: face_push_elemdata_to_face
    public :: face_push_diag_adjacent_data_to_face
    public :: face_interpolation
    public :: face_interpolate_bc
    public :: face_force_JBadjacent_values
    public :: face_velocities


    public :: face_flowrate_for_openclosed_elem

    public :: face_zeroDepth

    !% --- these later will be private
    public :: face_zeroDepth_geometry_interior
    public :: face_zeroDepth_geometry_shared
    public :: face_zeroDepth_flowrates_interior
    public :: face_zeroDepth_flowrates_shared

    
    ! public :: face_FluxCorrection_interior
    ! public :: face_flowrate_max_interior
    ! public :: face_flowrate_max_shared

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%==========================================================================
!%
    subroutine face_zeroDepth (fp_downstream, fp_upstream, fp_bothsides)
        !%------------------------------------------------------------------
        !% Description:
        !% sets zero depth geometry values on faces depending on the
        !%  adjacent elements
        !% For CC
        !%  arguments are fp_CC_downstream_is_zero_IorS,
        !%  fp_CC_upstream_is_zero_IorS, fp_CC_bothsides_are_zero_IorS
        !% For JB
        !%  rguments are fp_JB_downstream_is_zero_IorS,
        !%  fp_JB_upstream_is_zero_IorS, fp_JB_bothsides_are_zero_IorS
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: fp_downstream, fp_upstream, fp_bothsides
        !%------------------------------------------------------------------

        !% --- geometry (interior)
        call face_zeroDepth_geometry_interior(fp_downstream)
        call face_zeroDepth_geometry_interior(fp_upstream)
        call face_zeroDepth_geometry_interior(fp_bothsides)

        !% --- geometry (shared)
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

    end subroutine face_zeroDepth
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine face_push_elemdata_to_face (thisPCol, frCol, erCol, elemXR, UpstreamFaceTF)
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

    end subroutine face_push_elemdata_to_face
!%
!%==========================================================================    
!%==========================================================================
!%
    subroutine face_push_diag_adjacent_data_to_face (thisPCol )
        !%------------------------------------------------------------------
        !% Description
        !% Pushes element data into the fr_..._adjacent data fields
        !% thsiPCol must be a packed set of diagnostic elements that are
        !% adjacent to JB elements
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisPCol
            integer, pointer    :: Npack, thisP(:), fup, fdn
            integer :: ii, kk, ff
        !%------------------------------------------------------------------
            Npack => npack_elemP(thisPCol)
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

        ! print *, ' '
        ! print *, 'in face_pull '
        ! print *, 'thisP ',thisP   
        ! print *, ' ' 
        ! print *, elemSI(thisP,esi_JunctionBranch_IsUpstream)
        ! print *, ' '
        ! print *, faceR(elemI(thisP(1),ei_Mface_uL),frCol)
        ! print *, faceR(elemI(thisP(2),ei_Mface_dL),frCol)
        ! print *, ' '
        ! print *, elemI(thisP(1),ei_Mface_uL), elemI(thisP(2),ei_Mface_dL)
        ! print *, ' '
        ! print *, trim(reverseKey(elemI(7,ei_elementType)))
        ! print *, 'upstream?   ', elemSI(7,esi_JunctionBranch_IsUpstream)
        ! print *, 'face map up ',elemI(7,ei_Mface_uL)
        ! print *, 'elem map up ',faceI(elemI(7,ei_Mface_uL),fi_Melem_uL)
        ! print *, ' '
        ! print *, trim(reverseKey(elemI(8,ei_elementType)))
        ! print *, 'upstream?   ', elemSI(8,esi_JunctionBranch_IsUpstream)
        ! print *, 'face map dn ',elemI(8,ei_Mface_dL)
        ! print *, 'elem map dn ',faceI(elemI(8,ei_Mface_dL),fi_Melem_dL)
        ! print *, ' '

    

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
    subroutine face_interpolation (facecol, Gyn, Hyn, Qyn, skipJump, skipZeroAdjust)
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
            integer, pointer :: Npack
            logical :: isBConly !, isTM
            integer :: iblank
            
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

            ! call util_utest_CLprint ('    XXX01 face after face_interpolation_interior')

        !% MOVED OUT OF THE INTERPOLATION ROUTINE 20230419
        ! !% --- force zero fluxes on downstream faces of closed CC element
        ! !%     This is the EPA-SWMM "Setting" value for open/closing elements
        ! !%     and should not be confused with
        ! !%     the setting%ZeroDepth algorithm in SWMM5+
        ! !%     note this does not require a "faceCol" argument as we
        ! !%     will execute this for both fp_noBC_IorS and fp_Diag_IorS calls
        ! if (Qyn) then
        !     call adjust_face_for_zero_setting ()
        ! end if

          !!!  ! call util_utest_CLprint ('    XXX02 face after adjust face for zero setting')

        !% REMOVE THIS FROM HERE 20230420
        ! if (.not. skipZeroAdjust) then

        !     print *, 'calling in not skipZeroAdjust'

        !     call face_zerodepth_interior(fp_elem_downstream_is_zero)

        !         ! call util_utest_CLprint ('    XXX03 face after face_zerodepth_interior 1')

        !     call face_zerodepth_interior(fp_elem_upstream_is_zero)
                
        !         ! call util_utest_CLprint ('    XXX04 face after face_zerodepth_interior 2')

        !         ! print *, ' '
        !         ! print *, 'in zero depth interior bothsides are zero index AAAA'
        !         ! print *, 'THISF FACE COL ',faceCol, fp_noBC_IorS
        !         ! print *, faceP(1:3,fp_elem_bothsides_are_zero)
        !         ! print *, ' '
        !         ! print *, faceP(1:npack_faceP(fp_elem_bothsides_are_zero),fp_elem_bothsides_are_zero)
        !         ! print *, ' '

        !     !call face_zerodepth_interior(fp_elem_bothsides_are_zero)
        ! end if

        !% --- face reconstruction of all the shared faces
        call face_interpolation_shared (faceCol, Gyn, Hyn, Qyn, .true.)

            ! call util_utest_CLprint ('    XXX06 face after face interpolation_shared')

        call face_interpolate_bc (isBConly)

            ! call util_utest_CLprint ('    XXX07 face after face_interpolate_BC')

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

        !% brh20211211 MOVED -- this is an element update
        !rm if (N_nBClat > 0) call face_interpolation_latBC_byPack()

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
    subroutine face_force_JBadjacent_values (thisColP, isJBupstreamYN)
        !%------------------------------------------------------------------
        !% Description
        !% Forces faces adjacent to JB to JB values without interpolation
        !%------------------------------------------------------------------
        !% Declaration
            integer, intent(in) :: thisColP
            logical, intent(in) :: IsJBupstreamYN !% true if JB is upstream of JM
            integer, pointer    :: thisJM(:), Npack
            integer :: mm, ei_Mface, kstart, kk
        !%------------------------------------------------------------------
        !% Preliminaries
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

            ! print *, 'mm',mm
            ! print *, 'thisJM(mm)', thisJM(mm)
            ! print *, 'IS upstream ? ',isJBUpstreamYN


            !% --- push junction JB flowrate values back to faces  
            call face_force_JBvalues (fr_Flowrate, er_Flowrate, ei_Mface, thisJM(mm), kstart) 

            !% ---add junction JB Delta Q values to face value caused by interpolation
            !call face_add_JBvalues (fr_DeltaQ, er_DeltaQ, ei_Mface, thisJM(mm), kstart) 
            call face_force_JBvalues (fr_DeltaQ, er_DeltaQ, ei_Mface, thisJM(mm), kstart) 


            !% 20230427 APPEARS TO WORK BETTER WITHOUT FORCING HEAD
            !% INSTEAD THE INTERPOLATION WILL DOMINATE
            !% ALSO CAUSED NONCONSERVATION ISSUES
            !% --- push JB head to adjacent face
            ! call face_force_JBvalues (fr_Head_d, er_Head, ei_Mface, thisJM(mm), kstart) 
            ! call face_force_JBvalues (fr_Head_u, er_Head, ei_Mface, thisJM(mm), kstart) 

            !% --- fix the face depth and area
            !% 20230427 CAUSES MASS CONSERVATION PROBLEMS IN CASE T007b_RO_Free-dx0010.inp
            !% eventually leads to segmenation fault
            ! do kk=kstart,max_branch_per_node,2
            !     if (elemSI(thisJM(mm)+kk,esi_JunctionBranch_Exists).ne. oneI) cycle
        
            !     !% --- get the depth
            !     faceR(elemI(thisJM(mm)+kk,ei_Mface),fr_Depth_u) &
            !          = faceR(elemI(thisJM(mm)+kk,ei_Mface),fr_Head_u)  &
            !          - faceR(elemI(thisJM(mm)+kk,ei_Mface),fr_Zbottom)

            !     !% --- use zerovalue%depth if depth is small
            !     faceR(elemI(thisJM(mm)+kk,ei_Mface),fr_Depth_u) &
            !         = max(faceR(elemI(thisJM(mm)+kk,ei_Mface),fr_Depth_u),0.99d0*setting%ZeroValue%Depth)  

            !     !% --- store depth on other side of face
            !     faceR(      elemI(thisJM(mm)+kk,ei_Mface),fr_Depth_d) &
            !         = faceR(elemI(thisJM(mm)+kk,ei_Mface),fr_Depth_u)

            !     !% --- compute area 
            !     if (faceR(elemI(thisJM(mm)+kk,ei_Mface),fr_Depth_u) > setting%ZeroValue%Depth) then
                    
            !         faceR(elemI(thisJM(mm)+kk,ei_Mface),fr_Area_u) &
            !             = geo_area_from_depth_singular (thisJM(mm)+kk,                                   &
            !                                             faceR(elemI(thisJM(mm)+kk,ei_Mface),fr_Depth_u), &
            !                                             setting%ZeroValue%Area)
            !     else 
            !         faceR(elemI(thisJM(mm)+kk,ei_Mface),fr_Area_u) = setting%ZeroValue%Area
            !     end if

            !     !% --- store depth on other side of face
            !     faceR(      elemI(thisJM(mm)+kk,ei_Mface),fr_Area_d) &
            !         = faceR(elemI(thisJM(mm)+kk,ei_Mface),fr_Area_u)
            ! end do

        end do

    end subroutine face_force_JBadjacent_values
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
            character(64) :: subroutine_name = 'face_interpolation_upBC'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%boundary_conditions)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% Aliases
        !% For the head/geometry at the upstream faces, we directly take the dnwnstream element
        !% So there is no eup for upstream BCs
            edn       => faceI(:,fi_Melem_dL)
            fdn       => elemI(:,ei_Mface_dL)   
            idx_fBC   => faceP(1:npack_faceP(fp_BCup),fp_BCup)
            idx_fJ1   => faceP(1:npack_faceP(fp_J1),  fp_J1)
            idx_fBoth => faceP(1:npack_faceP(fp_J1_BCup),  fp_J1_BCup)
            idx_P     => BC%P%BCup(:)
        !%-------------------------------------------------------------------
        !% enforce stored inflow BC    
        faceR(idx_fBC, fr_Flowrate) = BC%flowR(idx_P,br_value)

        ! print *, 'in face_interpolation_upBC'
        ! print *, faceR(idx_fBC,fr_Area_d), faceR(idx_fBC, fr_Flowrate)

        !% enforce zero flow on J1 faces
        faceR(idx_fJ1, fr_Flowrate) = zeroR
        faceR(idx_fJ1, fr_Velocity_u) = zeroR
        faceR(idx_fJ1, fr_Velocity_d) = zeroR

        !% update geometry data (don't do on a BC-only call)
        if (.not. isBConly) then
            !% Define sets of points for the interpolation, we are going from
            !% the elements to the faces.       
            !fGeoSetU = [fr_Area_u, fr_Topwidth_u, fr_HydDepth_u]
            !fGeoSetD = [fr_Area_d, fr_Topwidth_d, fr_HydDepth_d]
            !eGeoSet  = [er_Area,   er_Topwidth,   er_EllDepth]

            fGeoSetU = [fr_Area_u, fr_Depth_u]
            fGeoSetD = [fr_Area_d, fr_Depth_d]
            eGeoSet  = [er_Area,   er_Depth]

            !% Copying geometric data from elements to the BC/J1 faces
            do ii = 1,size(fGeoSetU)
                faceR(idx_fBoth,fGeoSetD(ii)) = elemR(edn(idx_fBoth),eGeoSet(ii)) 
                !% upstream side of face matches downstream (no jump)
                faceR(idx_fBoth,fGeoSetU(ii)) = faceR(idx_fBoth,fGeoSetD(ii))
            end do

            !% --- HACK: copy the preissmann number as well
            !%     Note that for static slot this will always be unity
            faceR(idx_fBoth,fr_Preissmann_Number) = elemR(edn(idx_fBoth),er_Preissmann_Number) 
            
            !% gradient extrapolation for head at infow
            faceR(idx_fBC, fr_Head_d) = elemR(edn(idx_fBC),er_Head)       &
                                      + elemR(edn(idx_fBC),er_Head)        &
                                      - faceR(fdn(edn(idx_fBC)),fr_Head_d) 

            !% zero head gradient at J1 cells
            faceR(idx_fJ1,fr_Head_d) = elemR(edn(idx_fJ1),er_Head) 

            !% head on downstream side of face is copied to upstream side (no jump)
            faceR(idx_fBoth,fr_Head_u) = faceR(idx_fBoth,fr_Head_d)
           
            !% ensure face area_u is not smaller than zerovalue
            where (faceR(idx_fBC,fr_Area_d) <= setting%ZeroValue%Area)
                faceR(idx_fBC,fr_Area_d)     = setting%ZeroValue%Area
            end where
            where (faceR(idx_fBC,fr_Area_u) <= setting%ZeroValue%Area)
                faceR(idx_fBC,fr_Area_u)     = setting%ZeroValue%Area
            endwhere

        end if

        !% set velocity based on flowrate
        faceR(idx_fBC,fr_Velocity_d) = faceR(idx_fBC,fr_Flowrate)/faceR(idx_fBC,fr_Area_d)
        faceR(idx_fBC,fr_Velocity_u) = faceR(idx_fBC,fr_Velocity_d)  

        !%  If an inflow has a high velocity, reset the value of the high velocity limiter
        !%  
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

        !stop 209873

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
        !% Preliminaries
            
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
                            print *, 'CONFIGURATION ERROR: a NORMAL OUTFALL must be connected to an...'
                            print *, '...conduit/channel element with non-zero, positive bottom slope.'
                            print *, 'Problem for Outfall ',trim(node%Names(BC%headI(ii,bi_node_idx))%str)
                            print *, 'Connected to element ', eUp
                            print *, 'Part of link ',trim(  link%Names(elemI(eup,ei_link_Gidx_BIPquick))%str)
                            print *, 'Bottom slope is ',elemR(eup,er_BottomSlope)
                            call util_crashpoint(728474)
                        end if

                        !% --- get the normal depth
                        faceR(thisF,fr_Depth_u) = geo_normaldepth_singular (BC%HeadI(ii,bi_UTidx))

                        ! print *, ' '
                        ! print *, ' normal depth ',faceR(thisF,fr_Depth_u) 
                        ! print *, ' '

                        if (faceR(thisF,fr_Depth_u) > elemR(eup,er_Head)) then
                            !% --- use upstream head if the normal depth is too large (no backflow allowed)
                            faceR(thisF,fr_Depth_u) = elemR(eup,er_Head) - faceR(thisF,fr_Zbottom)
                        end if

 

                        !% --- head is the normal depth + Zbottom - referencehead
                        faceR(thisF,fr_Head_u)  = faceR(thisF,fr_Zbottom)                   &
                                                + faceR(thisF,fr_Depth_u)                  &
                                                - setting%Solver%ReferenceHead

                        !% --- Outfall head should not be larger than the upstream for normal BC
                        !faceR(thisF,fr_Head_u) = min(faceR(thisF,fr_Head_u),eHead(eup))    

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

                    ! print *, ' '
                    ! print *, 'upstream   ', elemR(eup,er_Depth), eDepth(eup)
                    ! print *, 'norm depth ', normDepth 
                    ! print *, 'crit depth ', critDepth 
                    ! print *, 'crit area  ', critArea
                    ! print *, 'Flowrate   ', elemR(eup,er_Flowrate), critArea * sqrt(setting%Constant%gravity * critDepth)
                    ! print *, ' '
                    
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


            ! print *, ' '
            ! print *, 'depth, area ',faceR(thisF,fr_Depth_u), faceR(thisF,fr_Area_u)
            ! print *, 'flowrate ',elemR(eup,er_Flowrate), faceR(thisF,fr_Flowrate)
            ! print *, 'Asqrt(gh)',faceR(thisF,fr_Area_u) * sqrt(setting%Constant%gravity * faceR(thisF,fr_Depth_u))
            ! print *, 'FR ', (faceR(thisF,fr_Flowrate) / faceR(thisF,fr_Area_u)) / sqrt(setting%Constant%gravity * faceR(thisF,fr_Depth_u))
            ! print *, ' '

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
            integer, intent(in) :: facePackCol  !% Column in faceP array for needed pack
            logical, intent(in) :: Gyn, Hyn, Qyn
            logical, intent(in) :: skipJump !% this is true if jump_compute is skipped 
            integer, pointer    ::  Npack        !% expected number of packed rows in faceP.
            integer :: fGeoSetU(2), fGeoSetD(2), eGeoSet(2)
            integer :: fHeadSetU(1), fHeadSetD(1), eHeadSet(1)
            integer :: fFlowSet(2), eFlowSet(2)

            real(8), pointer :: fQold(:)
            integer, pointer :: thisP(:)

            !integer :: fOtherSet(1), eOtherSet(1)
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
        !% interpolate to ..._u
        !% identify hydraulic jumps
        !% set .._u and ..d based on jumps

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

            ! if (facePackCol == fp_JB_IorS) then 
            !     !% --- store the old flowrate so we can compute the deltaQ
            !     thisP => faceP(1:Npack,facePackCol)
            !     fQold(thisP) = faceR(thisP,fr_Flowrate)
            ! end if
            
            ! print *, ' '
            ! print *, 'in face before '
            ! print *, elemR(10,er_Flowrate), faceR(11,fr_Flowrate), elemR(12,er_Flowrate)
            ! print *, ' '
            ! print *, ' '
            ! print *, 'in face before face interp '
            ! print *, elemR(50,er_InterpWeight_dQ), elemR(52,er_InterpWeight_uQ)
            ! print *, elemR(50,er_Flowrate), faceR(51,fr_Flowrate), elemR(52,er_Flowrate)
            ! print *, ' '

            call face_interp_interior_set &
                (fFlowSet, eFlowSet, er_InterpWeight_dQ, er_InterpWeight_uQ, facePackCol, Npack)

                ! print *, ' '
                ! print *, 'in face after face interp '
                ! print *, elemR(50,er_InterpWeight_dQ), elemR(52,er_InterpWeight_uQ)
                ! print *, elemR(50,er_Flowrate), faceR(51,fr_Flowrate), elemR(52,er_Flowrate)
                ! print *, 'value ',(elemR(50,er_InterpWeight_dQ) * elemR(52,er_Flowrate) &
                !                + elemR(52,er_InterpWeight_uQ) * elemR(50,er_Flowrate)) &
                !                /(elemR(50,er_InterpWeight_dQ)+ elemR(52,er_InterpWeight_uQ))
                ! print *, ' '

            ! if (facePackCol == fp_JB_IorS) then 
            !         ! print *, ' '
            !         ! print *, '                                  CALLING face_DeltaQ =============='
            !     call face_deltaQ (facePackCol,.true.,fQold)
            ! end if

            !% --- calculate the velocity in faces and put limiter
            call face_velocities (facePackCol, .true.)
            
            !% --- compute volume-based limits on flowrate
            call face_flowrate_limits_interior (facePackCol)
        end if

        !% --- reset all the hydraulic jump interior faces
        if (.not. skipJump) then
            call jump_compute
        end if

        !% REMOVE FOR RK_ETM_5
        !% --- for JB faces (only) store the adjacent head, topwidth, length values
        !if (facePackCol .ne. fp_notJB_interior) then
     !       call face_junction_adjacent_values (fp_JB_IorS)
        !end if

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
            !integer :: fOtherSet(1), eGhostOtherSet(1)
            integer(kind=8) :: crate, cmax, cval
            character(64) :: subroutine_name = 'face_interpolation_shared'
        !%-------------------------------------------------------------------
        !% Aliases
            Npack => npack_facePS(facePackCol)
        !%-------------------------------------------------------------------
        !% Preliminaries   
            if (Npack < 1) return

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

        !% General approach
        !% interpolate to ..._u
        !% identify hydraulic jumps
        !% set .._u and ..d based on jumps

        !% t--- transfer all the local elemR data needed for face interpolation into elemB data structure
        call local_data_transfer_to_boundary_array (facePackCol, Npack)

        !% --- use elemB to transfer remote data to local elemG array for interpolation
        call inter_image_data_transfer (facePackCol, Npack)

        !% set the matching sets
        if (Gyn) then
            fGeoSetU = [fr_Area_u, fr_Depth_u]
            fGeoSetD = [fr_Area_d, fr_Depth_d]
            eGeoSet  = [er_Area, er_Depth]

            call face_interp_shared_set &
                (fGeoSetU, eGeoSet, er_InterpWeight_dG, er_InterpWeight_uG, facePackCol, Npack)

            !% --- copy upstream to downstream storage at a face
            !%    (only for Head and Geometry types)
            !%    note that these might be reset by hydraulic jump
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

            ! print *,'  '
            ! print *, ' NEED the shared routine for fr_DeltaQ '
            ! stop 2098734

            call face_interp_shared_set &
                (fFlowSet, eFlowSet, er_InterpWeight_dQ, er_InterpWeight_uQ, facePackCol, Npack)

            call face_velocities (facePackCol, .false.)

            !% --- compute volume-based limits on flowrate
            call face_flowrate_limits_shared (facePackCol)
        end if

        !%fOtherSet      = [  fr_Preissmann_Number,   fr_GammaM,   fr_KJunction_MinorLoss,   fr_2B_psiL,   fr_EnergyHead]
        !eGhostOtherSet = [ebgr_Preissmann_Number, ebgr_GammaM, ebgr_KJunction_MinorLoss, ebgr_2B_psiL, ebgr_EnergyHead]
        !fOtherSet      = [  fr_Preissmann_Number]
        !eGhostOtherSet = [ebgr_Preissmann_Number]

        !call face_interp_shared_set &
        !    (fOtherSet, eGhostOtherSet, ebgr_InterpWeight_dQ, ebgr_InterpWeight_uQ, facePackCol, Npack)
            
        !% --- reset all the hydraulic jump interior faces
        if (.not. skipJump) then
            print *, 'ERROR: need to set up a test case for a hydraulic jump across a shared image'
            print *, 'need a version of the jump_compute for shared faces'
            call util_crashpoint(7209885)
        end if

        ! !% --- for JB faces (only) store the adjacent head, topwidth, length values
        ! if (facePackCol .ne. fp_notJB_interior) then
        !     !call ??? face_junction_adjacent_values (fp_JB_IorS)
        !     print *, 'ERROR: shared faces not complete'
        !     print *, 'need to handle the fr_..._Adjacent values as in face_junction_adjacent_values'
        !     call util_crashpoint(7229873)
        ! end if

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
            integer :: ii, jj
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
        !% cycle interpolation through each type in the set.

        ! print *, 'G ',elemR(49,er_InterpWeight_uG),elemR(49,er_InterpWeight_dG)

            ! print *, ' ', size(fset)
            ! print *, 'thisP in fset ',thisP


        do ii=1,size(fset)

            ! if (ii == 1) then 
            !     do jj=1,Npack
            !         if (thisP(jj) == 11) then
            !             print *, ' '
            !             print *, 'in face interpolation set'
            !             print *, thisP(jj), eWup, eWdn
            !             print *, 'eup, edn ',eup(thisP(jj)),edn(thisP(jj))
            !             print *, 'elemUp, elemD',elemR(eup(thisP(jj)),eset(ii)),elemR(edn(thisP(jj)),eset(ii))
            !             print *, 'weights E up, E dn  ',elemR(eup(thisP(jj)),eWdn),    elemR(edn(thisP(jj)),eWup)
            !             !print *, 'G ',elemR(thisP(),er_InterpWeight_uG),elemR(49,er_InterpWeight_dG)
            !             print *, ' '
            !         end if
            !     end do
            ! end if

      
            ! if (ii==1) then 
            !     print *, 'in interp for Q?'
            !     print *, 'working on: ',eup(51), 51, edn(51)
            !     print *, elemR(eup(51),eset(ii)), elemR(edn(51),eWup)
            !     print *, elemR(edn(51),eset(ii)), elemR(eup(51),eWdn)
            !     print *, ' ' 
            ! end if
      

            faceR(thisP,fset(ii)) = &
                (+elemR(eup(thisP),eset(ii)) * elemR(edn(thisP),eWup) &
                 +elemR(edn(thisP),eset(ii)) * elemR(eup(thisP),eWdn) &
                ) / &
                ( elemR(edn(thisP),eWup) + elemR(eup(thisP),eWdn))

            ! print *, 'face ',faceR(thisP(1),fset(ii))    
            ! print *, ' '
        end do

        !% NOTES
        !% elemR(eup(thisP),eset(ii)) is the element value upstream of the face
        !% elemR(edn(thisP),eset(ii) is the element value downstream of the face.
        !% elemR(eup(thisp),eWdn) is the downstream weighting of the upstream element----------------------------------------------------------

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
        !% NOTE cannot sync all in this subroutine
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
        !% cycle through all the shared faces
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
            !% cycle through each element in the set.
            !% This is designed for fset and eset being vectors, but it
            !%   is not clear that this is needed.
            do jj=1,size(fset)

                !% condition for upstream element of the shared face is ghost and in a different image
                if (isGhostUp) then
                    faceR(thisP,fset(jj)) = &
                        (+elemGR(ii,eset(jj))  * elemB%R(ii,eWup)  &
                         +elemB%R(ii,eset(jj)) * elemGR(ii,eWdn) &
                        ) / &
                        ( elemB%R(ii,eWup) + elemGR(ii,eWdn) )

                !% condition for downstream element of the shared face is ghost and in a different image
                elseif (isGhostDn) then

                    faceR(thisP,fset(jj)) = &
                        (+elemB%R(ii,eset(jj)) * elemGR(ii,eWup) &
                         +elemGR(ii,eset(jj))  * elemB%R(ii,eWdn)  &
                        ) / &
                        ( elemGR(ii,eWup) + elemB%R(ii,eWdn) )
                else
                    write(*,*) 'CODE ERROR: unexpected else'
                    !stop 
                    call util_crashpoint( 487874)
                    !return
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

            !% saz 20230414 
            !% HACK: elemB%R/elemGR(:,:) will share the same number of columns as elemR
            !% thus we can readily use er_.... column indexs

            ! !% HACK: this eset has to be exactly mimic the indexes for ebgr_... 
            ! eColumns = [er_Area, er_Topwidth, er_Depth, er_DeltaQ, &
            !             er_Head, er_Flowrate, er_Preissmann_Number, er_Volume,     &
            !             er_Velocity, er_GammaM, er_Length, er_KJunction_MinorLoss, &
            !             er_InterpWeight_dG, er_InterpWeight_uG,                    &
            !             er_InterpWeight_dH, er_InterpWeight_uH,                    &
            !             er_InterpWeight_dQ, er_InterpWeight_uQ,                    &
            !             er_InterpWeight_dP, er_InterpWeight_uP] 

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

            !print *, 'xxAA ',this_image(), ii, thisP, isGhostUp, isGhostDn, eUp, eDn
            !% condition for upstream element is ghost
            if (isGhostUp) then
                elemB%R(ii,:) = elemR(eDn,:)
            !% condition for downstream element is ghost
            elseif (isGhostDn) then
                elemB%R(ii,:) = elemR(eUp,:)
            else
                write(*,*) 'CODE ERROR: unexpected else'
                !stop 
                call util_crashpoint( 487874)
                !return
            end if     
            !% --- handle special case for volume used by Pump Type 1 when
            !%     the upstream element is a JB
            if ((isGhostUp) .and. (elemI(eUp,ei_elementType) == JB)) then
                !JMidx => elemI(eup,ei_main_idx_for_branch)
                JMidx => elemSI(eup,esi_JunctionBranch_Main_Index)
                elemB%R(ii,er_Volume) = elemR(JMidx,er_Volume)
            end if
        end do

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
            integer, pointer    :: thisP, ci, BUpIdx, BDnIdx, eUp, eDn
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
        !% cycle through all the shared faces
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
            !% condition for upstream element of the shared face is ghost and in a different image
            if (isGhostUp) then
                elemGR(ii,:) = elemB[ci]%R(BUpIdx,:) 
                ! print*, elemGR(ii,:), 'elemGR(ii,:)'
                ! print*
                ! print*, elemB[ci]%R(BUpIdx,:) , 'elemB[ci]%R(BUpIdx,:) '
                ! print*
                ! print*, reverseKey(elemI(faceI(thisP,fi_GhostElem_uL), ei_elementType)[ci])
                ! print*, elemI(faceI(thisP,fi_GhostElem_uL), :)[ci], 'elem row'
            !% condition for downstream element of the shared face is ghost and in a different image
            elseif (isGhostDn) then
                elemGR(ii,:) = elemB[ci]%R(BDnIdx,:)
                ! print*, elemGR(ii,:), 'elemGR(ii,:)'
                ! print*
                ! print*, elemB[ci]%R(BUpIdx,:) , 'elemB[ci]%R(BUpIdx,:) '
                ! print*
                ! print*, reverseKey(elemI(faceI(thisP,fi_GhostElem_dL), ei_elementType)[ci])
                ! print*, elemR(faceI(thisP,fi_GhostElem_dL), :)[ci], 'elem row'
            else
                write(*,*) 'CODE ERROR: unexpected else'
                !stop 
                call util_crashpoint( 487874)
                !return
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
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Copies the interpolated value on the upstream side to the downstream side
        !% These values are later adjusted for hydraulic jumps
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: facePackCol, Npack, downstreamSet(:), upstreamSet(:)
        integer, pointer :: thisP(:)
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'face_copy_upstream_to_downstream_interior'
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%face) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------

        if (Npack > 0) then
            thisP => faceP(1:Npack,facePackCol)
            faceR(thisP,downstreamSet) = faceR(thisP,upstreamSet)
        end if

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
        !%-----------------------------------------------------------------------------

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
    subroutine face_deltaQ (facePackCol, isInterior, faceQold)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the change in face DQ caused by face interpolation for
        !% JB as part of junctions. This must be included in the total deltaQ
        !% change caused by the junction
        !%-------------------------------------------------------------------  
        !% Declarations:
            integer, intent(in) :: facePackCol
            logical, intent(in) :: isInterior
            real(8), intent(in) :: faceQold(:)

            integer, pointer :: Npack, thisP(:)

            integer :: ii
        !%-------------------------------------------------------------------  
        !% Preliminaries
            if (isInterior) then
                !% calculation at the interior faces
                Npack => npack_faceP(facePackCol)
                thisP => faceP(1:Npack,facePackCol)
            else
                !% calculation at the shared faces
                Npack => npack_facePS(facePackCol)
                thisP => facePS(1:Npack,facePackCol)
            end if
            if (Npack < 1) return
        !%-------------------------------------------------------------------      

        faceR(thisP,fr_DeltaQ) = faceR(thisP,fr_Flowrate) - faceQold(thisP) 

        ! print *, ' '
        ! print *, 'in face Delta Q'
        ! do ii=1,size(thisP)
        !     print *, ii, thisP(ii), faceR(thisP(ii),fr_Flowrate), faceQold(thisP(ii))
        ! end do

    end subroutine face_deltaQ
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

        ! print *, 'thisP ',thisP
        ! print *, 'eup(thisP)', eup(thisP)
        ! print *, 'edn(thisP)', edn(thisP)
        ! print *, faceR(thisP,fr_FlowrateMax)
        ! print *, faceR(thisP,fr_FlowrateMin)
        ! stop 43534
        
   

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
                faceR(thisP,fr_FlowrateMax) = elemGR(ii,er_Volume) / dt
            elseif (isGhostDn) then 
                faceR(thisP,fr_FlowrateMin) = -elemGR(ii,er_Volume) / dt 
            end if

        end do


    end subroutine face_flowrate_limits_shared
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine face_head_limited (facePackCol)
        !%-------------------------------------------------------------------
        !% Description:
        !% Finds where Head < Zbottom on face and sets Depth and Area to zero
        !%-------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: facePackCol
            integer, pointer :: Npack
            integer, pointer :: thisP(:)
        !%-------------------------------------------------------------------
        !% Aliases
            Npack => npack_faceP(facePackCol)
            if (Npack < 1) return
            thisP   => faceP(1:Npack,facePackCol)
        !%-------------------------------------------------------------------

        ! where ((onehalfR * (faceR(thisP,fr_Head_u) + faceR(thisP,fr_Head_d))) &
        !         .le. faceR(thisP,fr_Zbottom))
        !     faceR(thisP,fr_Depth_d)    = setting%ZeroValue%Depth
        !     faceR(thisP,fr_Depth_u)    = setting%ZeroValue%Depth
        !     faceR(thisP,fr_Area_d)     = setting%ZeroValue%Area
        !     faceR(thisP,fr_Area_u)     = setting%ZeroValue%Area
        !     faceR(thisP,fr_Velocity_d) = zeroR
        !     faceR(thisP,fr_Velocity_u) = zeroR
        !     faceR(thisP,fr_Flowrate)   = zeroR
        ! end where


    end subroutine face_head_limited
!%  
!%========================================================================== 
!%==========================================================================
!%
    subroutine face_CC_zerodepth_interior (facePackCol)
        !%------------------------------------------------------------------
        !% Description:
        !% where one side has a zero depth element, the face head is
        !% adjusted to the smaller of (1) the head computed by interpolation
        !% or (2) the head on the non-zero depth element.
        !% Input: column in the faceP(:,:) array containing the packed 
        !% indexes of faces with zero elements on upstream, downstream, o
        !% both sides.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: facePackCol
            real(8), pointer :: fHeadDn(:), fHeadUp(:), fDepthDn(:), fDepthUp(:)
            real(8), pointer :: fAreaDn(:), fAreaUp(:), fFlowrate(:)
            real(8), pointer :: fVelocityDn(:), fVelocityUp(:), fZbottom(:)
            real(8), pointer :: eHead(:), eFlowrate(:), eArea(:), eVelocity(:)
            integer, pointer :: eDn(:), eUp(:)
            integer, pointer :: mm, Npack, thisP(:)
            integer :: ii
        !%------------------------------------------------------------------
        !% Aliases:
            fHeadDn     => faceR(:,fr_Head_d)
            fHeadUp     => faceR(:,fr_Head_u)
            fDepthDn    => faceR(:,fr_Depth_d)
            fDepthUp    => faceR(:,fr_Depth_u)
            fAreaDn     => faceR(:,fr_Area_d)
            fAreaUp     => faceR(:,fr_Area_u)
            fVelocityDn => faceR(:,fr_Velocity_d)
            fVelocityUp => faceR(:,fr_Velocity_u)
            fFlowrate   => faceR(:,fr_Flowrate)
            fZbottom    => faceR(:,fr_Zbottom)
            eHead    => elemR(:,er_Head)
            eArea    => elemR(:,er_Area)
            eFlowrate=> elemR(:,er_Flowrate)
            eVelocity=> elemR(:,er_Velocity)
            eDn      => faceI(:,fi_Melem_dL)
            eUp      => faceI(:,fi_Melem_uL)
        !%------------------------------------------------------------------

        !% --- CASE: Face has zero element upstream
        Npack => npack_faceP(facePackCol)

        ! print *, 'in face zerodepth interior'
        ! print *, 'faces: ', faceP(1:Npack,facePackCol)
        ! print *, 'facePackCol ',facePackCol, fp_elem_downstream_is_zero, fp_elem_upstream_is_zero, fp_elem_bothsides_are_zero
        
        if (Npack > 0) then 
            thisP => faceP(1:Npack,facePackCol) 

            ! print *, ' '
            ! print *, 'in face_zerodepth_interior'
            ! print *, 'thisP ',thisP, eUp(thisP)

            !% --- set the head on the face for elements that have adjacent zero
            select case (facePackCol)

                case (fp_CC_downstream_is_zero_IorS)

                        ! print *, 'in fp_elem_downstream_is_zero'
                    !% ---set head to the smaller of the face head and the non-zero element upstream

                        ! print *, 'Head ',fHeadUp(553), eHead(eUp(553)), fHeadDn(553)

                    fHeadUp(thisP) = min(fHeadUp(thisP), eHead(eUp(thisP)))
                    fHeadDn(thisP) = fHeadUp(thisP)

                        ! print *, 'fHead(thisP)',fHeadUp(553), fHeadDn(553)
                        ! print *, 'zbottom     ',fZbottom(553)

                    !% --- get a face depth consistent with this head
                    !%     note this depth might be negative
                    fDepthUp(thisP) = max(fHeadUp(thisP) - fZbottom(thisP), zeroR)
                    fDepthDn(thisP) = fDepthUp(thisP)

                        ! print *, 'fDepthUp ',fDepthUp(553)
                        ! print *, 'fDepthDn ',fDepthDn(553)

                    !% --- get face area consistent with this depth
                    do ii=1,Npack
                        mm => thisP(ii) 
                            ! print *, 'mm ', mm
                            ! print *, 'eUp',eUp(mm)
                            ! print *, 'fDepthUp ',fDepthUp(mm)
                        fAreaUp(mm) = geo_area_from_depth_singular(eUp(mm), fDepthUp(mm), zeroR)
                    end do
                    fAreaDn(thisP) = fAreaUp(thisP)

                        ! print *, 'fDepthUp(thisP)   ', fDepthUp(thisP)
                        ! print *, 'fAreaUp(thisP)    ', fAreaUp(thisP) 
                        ! print *, 'eArea(eUp(thisP)) ',eArea(eUp(thisP))

                    !% --- set the flowrate through the face
                    !%     flowrate can only be from upstream to downstream (positive flow allowed) 
                    where ((fDepthUp(thisP) > setting%ZeroValue%Depth) .and. &
                           (fAreaUp(thisP) .ge. eArea(eUp(thisP))))
                        !%--- if the depth > 0 and the face area is greater than the element area
                        !%     then the face flow can support the entire flow rate (if an outflow)
                        fFlowrate(thisP)   = max(zeroR,eFlowrate(eUp(thisP)))

                    elsewhere ((fDepthUp(thisP) > setting%ZeroValue%Depth) .and. &
                               (fAreaUp(thisP)  < eArea(eUp(thisP))))
                        !% --- the outflow is reduced in proportion with the area
                        !%     i.e., the velocity of the element is used for the outflow  
                        !%     note velocity must be negative or outflow is zero 
                        fFlowrate(thisP) = max(zeroR, eVelocity(eUp(thisP)) * fAreaUp(thisP))
                        
                    elsewhere !% for face depth < zero there is no possibility of a face flow
                        fFlowrate(thisP) = zeroR
                    endwhere
                    
                    ! print *, ' '
                    ! print *, 'eflowrate ',eFlowrate(eUp(5))
                    ! print *, 'vel rate  ',eVelocity(eUp(5)) * fAreaUp(5)
                    ! print *, 'vel       ',eVelocity(eUp(5))
                    ! print *, 'area      ',fAreaUp(5)
                    ! print *, 'flowrate   ',fFlowrate(5)
                    ! print *, ' '

                    !% --- zero the face velocities to prevent small areas from being a problem
                    !%     this also means the advection of momentum into the element on the
                    !%     next time step is discounted until the element is no longer zero depth.
                    fVelocityUp(thisP) = zeroR
                    fVelocityDn(thisP) = zeroR

                case (fp_CC_upstream_is_zero_IorS)

                    ! print *, 'in fp_elem_upstream_is_zero'

                    !% ---set head to the smaller of the face head and the non-zero element downstream
                    fHeadDn(thisP) = min(fHeadDn(thisP), eHead(eDn(thisP)))
                    fHeadUp(thisP) = fHeadDn(thisP)

                    !% --- get a face depth consistent with this head
                    !%     note this depth might be negative
                    fDepthDn(thisP) = max(fHeadDn(thisP) - fZbottom(thisP), zeroR)
                    fDepthUp(thisP) = fDepthDn(thisP)

                    !% --- get face area consistent with this depth
                    do ii=1,Npack
                        mm => thisP(ii) 
                        fAreaDn(mm) = geo_area_from_depth_singular(eDn(mm), fDepthDn(mm), zeroR)
                    end do
                    fAreaUp(thisP) = fAreaDn(thisP)

                    !% --- set the flowrate through the face
                    !%     flowrate can only be from downstream to upstream (negative flow allowed)     
                    where ((fDepthDn(thisP) > setting%ZeroValue%Depth) .and. &
                           (fAreaDn(thisP) .ge. eArea(eDn(thisP))))
                        
                        fFlowrate(thisP)   = min(zeroR,eFlowrate(eDn(thisP)))

                    elsewhere ((fDepthDn(thisP) > setting%ZeroValue%Depth) .and. &
                             (fAreaDn(thisP) < eArea(eDn(thisP))))
                        !% --- the outflow is reduced in proportion with the area
                        !%     i.e., the velocity of the element is used for the outflow  
                        !%     note velocity must be negative or outflow is zero 
                        fFlowrate(thisP) = min(zeroR, eVelocity(eDn(thisP)) * fAreaDn(thisP))

                    elsewhere !% for face depth < zero there is no possibility of a face flow
                        fFlowrate(thisP) = zeroR
                    endwhere

                    !% --- zero the face velocities to prevent small areas from being a problem
                    !%     this also means the advection of momentum into the element on the
                    !%     next time step is discounted until the element is no longer zero depth.
                    fVelocityUp(thisP) = zeroR
                    fVelocityDn(thisP) = zeroR

                    ! print *, fFlowrate(43), fDepthDN(thisP)
                    ! print *, eFlowrate(eDn(43))

                case (fp_CC_bothsides_are_zero_IorS)

                    !% --- keep the previously interpolated head and depth, 
                    !%     but ensure area and fluxes are zero
                    fAreaDn(thisP)     = zeroR
                    fAreaUp(thisP)     = zeroR
                    fFlowrate(thisP)   = zeroR
                    fVelocityUp(thisP) = zeroR
                    fVelocityDn(thisP) = zeroR

                case default
            end select
        end if

    end subroutine face_CC_zerodepth_interior
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine face_junction_adjacent_values  (facePackCol)
        !%------------------------------------------------------------------
        !% Description:
        !% stores the upstream element head (for upstream JB) or downstream
        !% element head (for downstream JB) that is needed in junction 
        !% solutions
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: facePackCol
            integer, pointer    :: fidx, eidx, Npack
            integer :: ii, mm
        !%------------------------------------------------------------------
        !% Preliminaries
            if (facePackCol .ne. fp_JB_IorS) return !% --- only used for JB faces
            Npack => npack_faceP(facePackCol)
            if (Npack < oneI) return !% --- only used if such faces exist
        !%------------------------------------------------------------------ 

        !% HACK violates no-neighbor rule
            
        do mm=1,Npack
            fidx => faceP(mm,facePackCol)  !% --- face index
            !% -- does this have JB downstream?
            if (elemSI(faceI(fidx,fi_Melem_dL),esi_JunctionBranch_Exists) == oneI) then
                eidx => faceI(fidx,fi_Melem_uL) !% --- upstream element NO NEIGHBOR VIOLATION
            end if 
            !% --- does this have JB upstream?
            if (elemSI(faceI(fidx,fi_Melem_uL),esi_JunctionBranch_Exists) == oneI) then 
                eidx => faceI(fidx,fi_Melem_dL) !% -- downstream element NO NEIGHBOR VIOLATION
            end if

            faceR(fidx,fr_Head_Adjacent)     = elemR(eidx,er_Head)
            faceR(fidx,fr_Topwidth_Adjacent) = elemR(eidx,er_TopWidth)
            faceR(fidx,fr_Length_Adjacent)   = elemR(eidx,er_Length)

            select case (elemI(eidx,ei_elementType))
                case (CC,outlet,pump)
                    faceR(fidx,fr_Zcrest_Adjacent) = elemR(eidx,er_Zbottom)
                case (weir)
                    faceR(fidx,fr_Zcrest_Adjacent) = elemSR(eidx,esr_Weir_Zcrest) 
                case (orifice)
                    faceR(fidx,fr_Zcrest_Adjacent) = elemSR(eidx,esr_Orifice_Zcrest)
                case default
                    print *, 'CODE ERROR: unexpected case default'
                    call util_crashpoint(629873)
            end select


        end do
 
    end subroutine face_junction_adjacent_values    
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

        !do concurrent (kk=kstart:max_branch_per_node:2)
        do kk=kstart,max_branch_per_node,2
            if (elemSI(JMidx+kk,esi_JunctionBranch_Exists).ne. oneI) cycle

            faceR(elemI(JMidx+kk,fiIdx),frCol) = elemR(JMidx+kk,erCol)

        end do
    
    end subroutine face_force_JBvalues
!%    
!%==========================================================================
!%==========================================================================
!% 
    subroutine face_add_JBvalues (frCol, erCol, fiIdx, JMidx, kstart)
        !%------------------------------------------------------------------
        !% Description:
        !% Adds the JB element value to the adjacent face value
        !%------------------------------------------------------------------
            integer, intent(in) :: frCol  !% column in faceR array for output
            integer, intent(in) :: erCol  !% column in elemR array for input
            integer, intent(in) :: fiIdx  !% face index column for up/dn map
            integer, intent(in) :: JMidx  !% junction main index
            integer, intent(in) :: kstart !% =1 for upstream, 2 for downstream
            integer :: kk
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        !do concurrent (kk=kstart:max_branch_per_node:2)
        do kk=kstart,max_branch_per_node,2
            if (elemSI(JMidx+kk,esi_JunctionBranch_Exists).ne. oneI) cycle

            faceR(elemI(JMidx+kk,fiIdx),frCol) = faceR(elemI(JMidx+kk,fiIdx),frCol) + elemR(JMidx+kk,erCol)

        end do
    
    end subroutine face_add_JBvalues
!%    
!%==========================================================================    
!%==========================================================================
!%   
    subroutine face_zeroDepth_geometry_interior (facePcol)
        !% -----------------------------------------------------------------
        !% Description:
        !% sets the face geoemtry values for interior faces with neighbor elements
        !% that are zero depth
        !% For CC
        !%  facePcol must be one of fp_CC_downstream_is_zero_IorS,
        !%  fp_CC_upstream_is_zero_IorS, fp_CC_bothsides_are_zero_IorS
        !% For JB
        !%  facePcol must be one of fp_JB_downstream_is_zero_IorS,
        !%  fp_JB_upstream_is_zero_IorS, fp_JB_bothsides_are_zero_IorS
        !% -----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: facePcol
            integer, pointer :: npack, thisP(:), edn(:), eup(:), mm
            integer :: ii
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

            ! print *, ' '
            ! print *, 'here AA ',faceR(5,fr_Head_u), elemR(eup(5),er_Head)
            ! print *, 'here zbt',faceR(5,fr_Zbottom),  elemR(eup(5),er_Depth)
            ! print *, 'here    ', elemR(eup(5),er_Depth) + faceR(5,fr_Zbottom)
            ! print *, ' '

            !% ---- head is the smaller of the recent face value, the upstream element value,
            !%      or the depth upstream applied to the face zbottom (new 20230430brh)
            faceR(thisP,fr_Head_u) = min(faceR(thisP,fr_Head_u), elemR(eup(thisP),er_Head), elemR(eup(thisP),er_Depth) + faceR(thisP,fr_Zbottom))
            !% --- store on downstream face
            faceR(thisP,fr_Head_d) = faceR(thisP,fr_Head_u)

            ! print *, ' '
            ! print *, 'here BB ',faceR(5,fr_Head_d)
            ! print *, ' '
            ! print *, 'thisP',thisP, eup(thisP)
            ! print *, min(faceR(5,fr_Head_u),elemR(4,er_Head),elemR(4,er_Depth)+faceR(5,fr_Zbottom))

            !% --- depth is the computed depth from head or zeroDepth
            faceR(thisP,fr_Depth_u) = max(faceR(thisP,fr_Head_u) - faceR(thisP,fr_Zbottom), 0.99d0 * setting%ZeroValue%Depth)
            faceR(thisP,fr_Depth_d) = faceR(thisP,fr_Depth_u)

            !% --- compute face area consistent with depth
            do ii=1,npack 
                mm => thisP(ii)
                if (faceR(mm,fr_Depth_u) > setting%ZeroValue%Depth) then 
                    faceR(mm,fr_Area_u) = geo_area_from_depth_singular(eup(mm),faceR(mm,fr_Depth_u),setting%ZeroValue%Area)
                else
                    faceR(mm,fr_Area_u) = setting%ZeroValue%Area
                end if
                faceR(mm,fr_Area_d) = faceR(mm,fr_Area_u)
            end do

        case (fp_CC_upstream_is_zero_IorS,fp_JB_upstream_is_zero_IorS)
            !% --- head is the smaller of the recent value or the downstream element value
            !% or the depth upstream applied to the face zbottom (new 20230430brh)
            faceR(thisP,fr_Head_d) = min(faceR(thisP,fr_Head_d), elemR(edn(thisP),er_Head), elemR(edn(thisP),er_Depth)+ faceR(thisP,fr_Zbottom))
            faceR(thisP,fr_Head_u) = faceR(thisP,fr_Head_d)

            !% --- depth is the computed depth from head or zeroDepth
            faceR(thisP,fr_Depth_d) = max(faceR(thisP,fr_Head_d) - faceR(thisP,fr_Zbottom), 0.99d0 * setting%ZeroValue%Depth)
            faceR(thisP,fr_Depth_u) = faceR(thisP,fr_Depth_d)

            !% --- compute face area consistent with depth
            do ii=1,npack 
                mm => thisP(ii)
                if (faceR(mm,fr_Depth_d) > setting%ZeroValue%Depth) then 
                    faceR(mm,fr_Area_d) = geo_area_from_depth_singular(edn(mm),faceR(mm,fr_Depth_d),setting%ZeroValue%Area)
                else
                    faceR(mm,fr_Area_d) = setting%ZeroValue%Area
                end if
                faceR(mm,fr_Area_u) = faceR(mm,fr_Area_d)
            end do

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
        !% For CC
        !%  facePcol must be one of fp_CC_downstream_is_zero_IorS,
        !%  fp_CC_upstream_is_zero_IorS, fp_CC_bothsides_are_zero_IorS
        !% For JB
        !%  facePcol must be one of fp_JB_downstream_is_zero_IorS,
        !%  fp_JB_upstream_is_zero_IorS, fp_JB_bothsides_are_zero_IorS
        !% -----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: facePcol
            integer, pointer :: npack, thisP, edn, eup, GUp, GDn, ci
            logical, pointer :: isGhostUp, isGhostDn
            integer :: ii, Ifidx
        !% -----------------------------------------------------------------
        !% -----------------------------------------------------------------
        !% Preliminaries
        !% -----------------------------------------------------------------
        npack => npack_facePS(facePCol)
        if (npack < 1) return

        !% sync all processors before the start of the subroutine
        sync all
        do ii = 1,Npack
            !%-----------------------------------------------------------------
            !% Aliases
            thisP           => facePS(ii,facePCol)
            ci              => faceI(thisP,fi_Connected_image)
            eup             => faceI(thisP,fi_Melem_uL)
            edn             => faceI(thisP,fi_Melem_dL)
            GUp             => faceI(thisP,fi_GhostElem_uL)
            GDn             => faceI(thisP,fi_GhostElem_dL)
            isGhostUp       => faceYN(thisP,fYN_isUpGhost)
            isGhostDn       => faceYN(thisP,fYN_isDnGhost)
            
            select case (facePcol)
            case (fp_CC_downstream_is_zero_IorS, fp_JB_downstream_is_zero_IorS)

                if (.not. isGhostUp) then
                    !% --- head is the smaller of the recent value or the downstream element value
                    !% or the depth upstream applied to the face zbottom (new 20230430brh)
                    faceR(thisP,fr_Head_u) = min(faceR(thisP,fr_Head_u), elemR(eup,er_Head), elemR(eup,er_Depth) + faceR(thisP,fr_Zbottom)) 

                    faceR(thisP,fr_Head_d) = faceR(thisP,fr_Head_u)
                    faceR(thisP,fr_Depth_u) = max(faceR(thisP,fr_Head_u) - faceR(thisP,fr_Zbottom), 0.99d0 * setting%ZeroValue%Depth)
                    faceR(thisP,fr_Depth_d) = faceR(thisP,fr_Depth_u)

                    if (faceR(thisP,fr_Depth_u) > setting%ZeroValue%Depth) then 
                        faceR(thisP,fr_Area_u) = geo_area_from_depth_singular(eup,faceR(thisP,fr_Depth_u),setting%ZeroValue%Area)
                    else
                        faceR(thisP,fr_Area_u) = setting%ZeroValue%Area
                    end if
                    faceR(thisP,fr_Area_d) = faceR(thisP,fr_Area_u) 

                    !% find the index of the indentical face in the connected image
                    Ifidx = elemI(GDn,ei_Mface_uL)[ci]

                    !% the face values should be identical apart from the newly adjusted values
                    !% transfer the whole data column to the indetical array 
                    faceR(Ifidx,:)[ci] = faceR(thisP,:)
                
                end if

            case (fp_CC_upstream_is_zero_IorS,fp_JB_upstream_is_zero_IorS)
                if (.not. isGhostDn) then
                    !% --- head is the smaller of the recent value or the downstream element value
                !% or the depth upstream applied to the face zbottom (new 20230430brh)
                faceR(thisP,fr_Head_d) = min(faceR(thisP,fr_Head_d), elemR(edn,er_Head), elemR(edn,er_Depth)+ faceR(thisP,fr_Zbottom))
                faceR(thisP,fr_Head_u) = faceR(thisP,fr_Head_d)

                !% --- depth is the computed depth from head or zeroDepth
                faceR(thisP,fr_Depth_d) = max(faceR(thisP,fr_Head_d) - faceR(thisP,fr_Zbottom), 0.99d0 * setting%ZeroValue%Depth)
                faceR(thisP,fr_Depth_u) = faceR(thisP,fr_Depth_d)
                
                if (faceR(thisP,fr_Depth_d) > setting%ZeroValue%Depth) then 
                    faceR(thisP,fr_Area_d) = geo_area_from_depth_singular(edn,faceR(thisP,fr_Depth_d),setting%ZeroValue%Area)
                else
                    faceR(thisP,fr_Area_d) = setting%ZeroValue%Area
                end if
                faceR(thisP,fr_Area_u) = faceR(thisP,fr_Area_d)

                !% find the index of the indentical face in the connected image
                Ifidx = elemI(GUp,ei_Mface_dL)[ci]

                !% the face values should be identical apart from the newly adjusted values
                !% transfer the whole data column to the indetical array 
                faceR(Ifidx,:)[ci] = faceR(thisP,:)

                end if
            
            case (fp_CC_bothsides_are_zero_IorS,fp_JB_bothsides_are_zero_IorS)
                faceR(thisP,fr_Depth_u) = 0.99d0 * setting%ZeroValue%Depth
                faceR(thisP,fr_Depth_d) = 0.99d0 * setting%ZeroValue%Depth
                faceR(thisP,fr_Head_u) = faceR(thisP,fr_Depth_u) + faceR(thisP,fr_Zbottom)
                faceR(thisP,fr_Head_d) = faceR(thisP,fr_Depth_d) + faceR(thisP,fr_Zbottom)
                faceR(thisP,fr_Area_u) = setting%ZeroValue%Area
                faceR(thisP,fr_Area_d) = setting%ZeroValue%Area

            end select
        
        end do 

        !% HACK: not sure if we need this sync
        sync all

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
        !% For CC
        !%  facePcol must be one of fp_CC_downstream_is_zero_IorS,
        !%  fp_CC_upstream_is_zero_IorS, fp_CC_bothsides_are_zero_IorS
        !% For JB
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
        !%     HACK: it would be better if we could enusre that all the
        !%     face areas are greater than zero so we wouldn't need the where
        !%     statements. But as of 20230114 there were areas that lead
        !%     to NAN values
        where (faceR(thisP,fr_Area_d) > setting%ZeroValue%Area)
            faceR(thisP,fr_Velocity_d) = faceR(thisP,fr_Flowrate) /  faceR(thisP,fr_Area_d)
        elsewhere
            faceR(thisP,fr_Velocity_d) = zeroR
        endwhere

            ! call util_utest_CLprint ('------- LLL.02D  in face zerodepth ')

        where (faceR(thisP,fr_Area_u) > setting%ZeroValue%Area)
            faceR(thisP,fr_Velocity_u) = faceR(thisP,fr_Flowrate) /  faceR(thisP,fr_Area_u)
        elsewhere
            faceR(thisP,fr_Velocity_u) = zeroR
        endwhere

            ! call util_utest_CLprint ('------- LLL.02E  in face zerodepth ')

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
        !% For CC
        !%  facePcol must be one of fp_CC_downstream_is_zero_IorS,
        !%  fp_CC_upstream_is_zero_IorS, fp_CC_bothsides_are_zero_IorS
        !% For JB
        !%  facePcol must be one of fp_JB_downstream_is_zero_IorS,
        !%  fp_JB_upstream_is_zero_IorS, fp_JB_bothsides_are_zero_IorS
        !% -----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: facePcol
            integer, pointer :: npack, thisP, edn, eup, gUp, gDn, ci
            logical, pointer :: isGhostUp, isGhostDn
            integer :: ii, Ifidx
        !% -----------------------------------------------------------------
            npack => npack_facePS(facePCol)
            if (npack < 1) return
        
        !% sync all processors before the start of the subroutine
        sync all
        do ii = 1,Npack
            !%-----------------------------------------------------------------
            !% Aliases
            thisP           => facePS(ii,facePCol)
            ci              => faceI(thisP,fi_Connected_image)
            eup             => faceI(thisP,fi_Melem_uL)
            edn             => faceI(thisP,fi_Melem_dL)
            gUp             => faceI(thisP,fi_GhostElem_uL)
            gDn             => faceI(thisP,fi_GhostElem_dL)
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
                    !% find the index of the indentical face in the connected image
                    Ifidx = elemI(gDn,ei_Mface_uL)[ci]
                    !% the face values should be identical apart from the newly adjusted values
                    !% transfer only the flowrate data column to the indetical location 
                    faceR(Ifidx,fr_Flowrate)[ci] = faceR(thisP,fr_Flowrate)

                end if

            case (fp_CC_upstream_is_zero_IorS,fp_JB_upstream_is_zero_IorS)
                
                if (.not. isGhostDn) then

                    if (elemR(edn,er_Flowrate) .ge. zeroR) then
                        faceR(thisP,fr_Flowrate) = elemR(edn,er_Flowrate)
                    else 
                        faceR(thisP,fr_Flowrate) = zeroR
                    end if
                    !% find the index of the indentical face in the connected image
                    Ifidx = elemI(gUp,ei_Mface_dL)[ci]
                    !% the face values should be identical apart from the newly adjusted values
                    !% transfer only the flowrate data column to the indetical location 
                    faceR(Ifidx,fr_Flowrate)[ci] = faceR(thisP,fr_Flowrate)

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

        !% not sure if we need this sync
        sync all

    end subroutine face_zeroDepth_flowrates_shared
!%  
!%========================================================================== 
!%==========================================================================
!%
    ! subroutine face_zerodepth_shared (facePackCol)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% where one side has a zero depth element, the face head is
    !     !% adjusted to the smaller of (1) the head computed by interpolation
    !     !% or (2) the head on the non-zero depth element.
    !     !% Input: column in the facePS(:,:) array containing the packed 
    !     !% indexes of faces with zero elements on upstream, downstream, o
    !     !% both sides.
    !     !%------------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in) :: facePackCol
    !         integer             :: ii
    !         integer, pointer    :: Npack, thisP, eup, edn, BUpIdx, BDnIdx
    !         real(8), pointer    :: fHeadDn, fHeadUp, fDepthDn, fDepthUp
    !         real(8), pointer    :: fAreaDn, fAreaUp, fFlowrate
    !         real(8), pointer    :: fVelocityDn, fVelocityUp, fZbottom
    !         real(8), pointer    :: eHead, eFlowrate, eArea, eVelocity
    !         logical, pointer    :: isGhostUp, isGhostDn
    !     !%------------------------------------------------------------------
    !     !% Aliases:
    !         ! fHeadDn     => faceR(:,fr_Head_d)
    !         ! fHeadUp     => faceR(:,fr_Head_u)
    !         ! fDepthDn    => faceR(:,fr_Depth_d)
    !         ! fDepthUp    => faceR(:,fr_Depth_u)
    !         ! fAreaDn     => faceR(:,fr_Area_d)
    !         ! fAreaUp     => faceR(:,fr_Area_u)
    !         ! fVelocityDn => faceR(:,fr_Velocity_d)
    !         ! fVelocityUp => faceR(:,fr_Velocity_u)
    !         ! fFlowrate   => faceR(:,fr_Flowrate)
    !         ! fZbottom    => faceR(:,fr_Zbottom)
    !         ! eHead    => elemR(:,er_Head)
    !         ! eArea    => elemR(:,er_Area)
    !         ! eFlowrate=> elemR(:,er_Flowrate)
    !         ! eVelocity=> elemR(:,er_Velocity)
    !         ! eDn      => faceI(:,fi_Melem_dL)
    !         ! eUp      => faceI(:,fi_Melem_uL)
    !     !%------------------------------------------------------------------
    !     !%------------------------------------------------------------------
    !         !% --- CASE: Face has zero element upstream
    !         Npack => npack_facePS(facePackCol)

    !         if (Npack > 0) then
    !             do ii = 1,Npack
    !                 !% Aliases
    !                 thisP       => facePS(ii,facePackCol)
    !                 eup         => faceI(thisP,fi_Melem_uL)
    !                 edn         => faceI(thisP,fi_Melem_dL)
    !                 BUpIdx      => faceI(thisP,fi_BoundaryElem_uL)
    !                 BDnIdx      => faceI(thisP,fi_BoundaryElem_dL)
    !                 isGhostUp   => faceYN(thisP,fYN_isUpGhost)
    !                 isGhostDn   => faceYN(thisP,fYN_isDnGhost)
    !                 fHeadDn     => faceR(thisP,fr_Head_d)
    !                 fHeadUp     => faceR(thisP,fr_Head_u)
    !                 fDepthDn    => faceR(thisP,fr_Depth_d)
    !                 fDepthUp    => faceR(thisP,fr_Depth_u)
    !                 fAreaDn     => faceR(thisP,fr_Area_d)
    !                 fAreaUp     => faceR(thisP,fr_Area_u)
    !                 fVelocityDn => faceR(thisP,fr_Velocity_d)
    !                 fVelocityUp => faceR(thisP,fr_Velocity_u)
    !                 fFlowrate   => faceR(thisP,fr_Flowrate)
    !                 fZbottom    => faceR(thisP,fr_Zbottom)
    !                 eDn         => faceI(thisP,fi_Melem_dL)
    !                 eUp         => faceI(thisP,fi_Melem_uL)

    !                 select case (facePackCol)

    !                     case (fp_CC_downstream_is_zero_IorS)

    !                         !% set up aliases to upstream (non-zero depth) values
    !                         if (isGhostUp) then
    !                             eHead     => elemGR(ii,ebgr_Head)
    !                             eArea     => elemGR(ii,ebgr_Area)
    !                             eFlowrate => elemGR(ii,ebgr_Flowrate)
    !                             eVelocity => elemGR(ii,ebgr_Velocity)
    !                         else
    !                             eHead     => elemB%R(ii,ebgr_Head)
    !                             eArea     => elemB%R(ii,ebgr_Area)
    !                             eFlowrate => elemB%R(ii,ebgr_Flowrate)
    !                             eVelocity => elemB%R(ii,ebgr_Velocity)
    !                         end if

    !                         !% ---set head to the smaller of the face head and the non-zero element upstream
    !                         fHeadUp = min(fHeadUp, eHead)  
    !                         fHeadDn = fHeadUp 

    !                         !% --- get a face depth consistent with this head
    !                         !%     note this depth might be negative
    !                         fDepthUp = max(fHeadUp - fZbottom, zeroR)
    !                         fDepthDn = fDepthUp
    !                         !% --- get face area consistent with this depth

    !                     case (fp_CC_upstream_is_zero_IorS)
    !                         print *, 'SHARED NOT COMPLETED'
    !                         call util_crashpoint(60987433)
                            
    !                     case (fp_CC_bothsides_are_zero_IorS)
    !                         print *, 'SHARED NOT COMPLETED'
    !                         call util_crashpoint(60987431)

    !                     case default
    !                         print *, 'CODE ERROR: unexpected case default'
    !                         call util_crashpoint(6091723)
    !                 end select

    !             end do
    !         else 
    !             !% no shared faces
    !         end if

    ! end subroutine face_zerodepth_shared
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
           endwhere
   
       end subroutine face_flowrate_for_openclosed_elem    
   !%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine face_flowrate_max_interior (facePackCol)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% stores the upstream element flowrate on the face.
    !     !% this is used in the JB algorithms
    !     !% "facePackCol" must be a face pack (faceP) array
    !     !% Stores the flowrate of the actual upstream element (i.e., where
    !     !% the flow is coming from if both flows are the same direction. 
    !     !% Stores the difference between the flows if the flows are in 
    !     !% opposite directions.
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: facePackCol
    !         integer, pointer :: npack, thisF(:), eup(:), edn(:)
    !         real(8), pointer :: fQmax(:), eQ(:)
    !     !%------------------------------------------------------------------
    !     !% Aliases:
    !         npack => npack_faceP(facePackCol)
    !         if (npack < 1) return
    !         thisF => faceP(1:npack,facePackCol)
    !         fQmax => faceR(:,fr_Flowrate_Max)
    !         eQ    => elemR(:,er_Flowrate)
    !         eup   => faceI(:,fi_Melem_uL)
    !         edn   => faceI(:,fi_Melem_dL)
    !     !%------------------------------------------------------------------
    !     !% use (1 + sign(U)) and (1 - sign(U)) trick to discriminate between
    !     !% flows downstream from the upstream element and flows upstream
    !     !% from the downstream element.
    !     fQmax(thisF) = onehalfR * (                                                         &
    !                 (oneR + util_sign_with_ones( eQ( eup(thisF) ) ) ) *  eQ( eup(thisF) )   &
    !               + (oneR - util_sign_with_ones( eQ (edn(thisF) ) ) ) *  eQ( edn(thisF) ) )


    !     ! print *, ' '
    !     ! print *, '-------------------------------------------------'
    !     ! print *, 'in face flowrate'
    !     ! print *, thisF
    !     ! print *, ' '
    !     ! print *, eup(thisF)
    !     ! print *, ' '
    !     ! print *, edn(thisF)
    !     ! print *, ' '
    !     ! print *, eQ( eup(thisF) )
    !     ! print *, ' '
    !     ! print *,  eQ( edn(thisF) )
    !     ! print *, ' '
    !     ! print *, fQmax(thisF)
    !     ! print *, ' '
    !     ! print *, fQmax(:)
    !     ! print *, '-------------------------------------------------'
    !     ! print *, ' '


    !     !%------------------------------------------------------------------
    ! end subroutine face_flowrate_max_interior
!%    
!%==========================================================================    
!%==========================================================================
!%
    ! subroutine face_flowrate_max_shared (facePackCol)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% stores the upstream element flowrate on the face.
    !     !% this is used in the JB algorithms
    !     !% "facePackCol" must be a face pack (facePS) array
    !     !% Stores the flowrate of the actual upstream element (i.e., where
    !     !% the flow is coming from if both flows are the same direction. 
    !     !% Stores the difference between the flows if the flows are in 
    !     !% opposite directions.
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: facePackCol
    !         integer, pointer    :: npack, thisF, ghostUp, ghostDn, ci
    !         integer, pointer    :: eup(:), edn(:)
    !         logical, pointer    :: isGhostUp, isGhostDn
    !         real(8), pointer    :: eQ(:), fQmax(:)
    !         integer             :: ii
    !         integer(kind=8) :: crate, cmax, cval
    !     !%------------------------------------------------------------------
    !     !% Preliminaries:
    !         !% start the shared timer
    !         sync all
    !         if (this_image()==1) then
    !             call system_clock(count=cval,count_rate=crate,count_max=cmax)
    !             setting%Time%WallClock%SharedStart = cval
    !             !setting%Time%WallClock%SharedStart_A = cval
    !         end if 
    !     !%------------------------------------------------------------------    
    !     !% Aliases
    !     !% note that aliases cannot be used where coarrays are invoked
    !         npack       => npack_facePS(facePackCol)
    !         if (npack < 1) return
    !         eup         => faceI(:,fi_Melem_uL)
    !         edn         => faceI(:,fi_Melem_dL)
    !         fQmax       => faceR(:,fr_Flowrate_Max)
    !     !%------------------------------------------------------------------
    !     !% cycle through the shared faces (does not readily vectorize)  
    !     do ii = 1,npack
    !         thisF       => facePS(ii,facePackCol)
    !         ci          => faceI(thisF,fi_Connected_image)
    !         ghostUp     => faceI(thisF,fi_GhostElem_uL)
    !         ghostDn     => faceI(thisF,fi_GhostElem_dL)
    !         isGhostUp   => faceYN(thisF,fYN_isUpGhost)
    !         isGhostDn   => faceYN(thisF,fYN_isDnGhost)

    !         !% condition for upstream element of the shared face is ghost and in a different image
    !         if (isGhostUp) then
    !             fQmax(thisF) = onehalfR * (                                                     &
    !                 (oneR + util_sign_with_ones( elemR(ghostUp,   er_Flowrate)[ci] ) ) * elemR(ghostUp,   er_Flowrate)[ci]  &
    !               + (oneR - util_sign_with_ones( elemR(edn(thisF),er_Flowrate)     ) ) * elemR(edn(thisF),er_Flowrate)     )
    !         !% condition for downstream element of the shared face is ghost and in a different image
    !         elseif (isGhostDn) then
    !             fQmax(thisF) = onehalfR * (                                                      &
    !                 (oneR + util_sign_with_ones( elemR(eup(thisF),er_Flowrate)     ) ) * elemR(eup(thisF),er_Flowrate)      &
    !               + (oneR - util_sign_with_ones( elemR(ghostDn,   er_Flowrate)[ci] ) ) * elemR(ghostDn,   er_Flowrate)[ci] )
    !         else
    !             write(*,*) 'CODE ERROR: unexpected else'
    !             !stop 
    !             call util_crashpoint( 88355)
    !             !return
    !         end if 
    !     end do

    !     !%------------------------------------------------------------------
    !     !% Closing
    !     sync all
    !     if (this_image()==1) then
    !         !% stop the shared timer
    !         call system_clock(count=cval,count_rate=crate,count_max=cmax)
    !         setting%Time%WallClock%SharedStop = cval
    !         setting%Time%WallClock%SharedCumulative &
    !                 = setting%Time%WallClock%SharedCumulative &
    !                 + setting%Time%WallClock%SharedStop &
    !                 - setting%Time%WallClock%SharedStart

    !         ! setting%Time%WallClock%SharedStop_A = cval
    !         ! setting%Time%WallClock%SharedCumulative_A &
    !         !         = setting%Time%WallClock%SharedCumulative_A &
    !         !         + setting%Time%WallClock%SharedStop_A &
    !         !         - setting%Time%WallClock%SharedStart_A                    
    !     end if 
    ! end subroutine face_flowrate_max_shared
!%
!%==========================================================================
    !%==========================================================================
!%
    ! subroutine face_interp_set_byMask &
    !     (fset, eset, eWdn, eWup, faceMaskCol)
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Interpolates to a face for a set of variables using a mask
    !     !%-----------------------------------------------------------------------------
    !     integer, intent(in) :: fset(:), eset(:), eWdn, eWup, faceMaskCol
    !     integer, pointer :: eup(:), edn(:)
    !     integer :: ii
        
    !     character(64) :: subroutine_name = 'face_interp_set_byMask'
    !     !%-----------------------------------------------------------------------------
    !     if (setting%Debug%File%face) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%-----------------------------------------------------------------------------
    !     eup => faceI(:,fi_Melem_uL)
    !     edn => faceI(:,fi_Melem_dL)
    !     !%-----------------------------------------------------------------------------
    !     !% cycle through each element in the set.
    !     do ii=1,size(fset)
    !         where (faceM(:,faceMaskCol))
    !             faceR(:,fset(ii)) = &
    !                 (+elemR(eup(:),eset(ii)) * elemR(edn(:),eWup) &
    !                  +elemR(edn(:),eset(ii)) * elemR(eup(:),eWdn) &
    !                 ) / &
    !                 ( elemR(edn(:),eWup) + elemR(eup(:),eWdn))
    !         endwhere
    !     end do

    !     print *, 'in face_interp_set_byMask -- may be obsolete'
    !     stop 87098

    !     if (setting%Debug%File%face) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine face_interp_set_byMask
    !%==========================================================================
!%
!% brh 20211211 moved and renamed -- this doesn't belong in Face routines
! subroutine face_interpolation_latBC_byPack()
        
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Interpolates all boundary faces using a pack arrays -- base on bi_category
    !     !%-----------------------------------------------------------------------------
    !     !integer :: fGeoSetU(3), fGeoSetD(3), eGeoSet(3)
    !     !integer :: fFlowSet(1), 
    !     integer :: ii
    !     integer, pointer :: elem_P(:), idx_P(:)
    !     integer :: eFlowSet(1)
    !     !integer :: fHeadSetU(1), fHeadSetD(1), eHeadSet(1)
    !     character(64) :: subroutine_name = 'face_interpolation_latBC_byPack'

    !     !%-----------------------------------------------------------------------------
    !     if (setting%Debug%File%boundary_conditions)  &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"


    !     elem_P => elemP(1:npack_elemP(ep_BClat),ep_BClat)
    !     idx_P  => BC%P%BClat

    !     !fGeoSetU = [fr_Area_u, fr_Topwidth_u, fr_HydDepth_u]
    !     !fGeoSetD = [fr_Area_d, fr_Topwidth_d, fr_HydDepth_d]
    !     !eGeoSet  = [er_Area,   er_Topwidth,   er_HydDepth]

    !     !fHeadSetU = [fr_Head_u]
    !     !fHeadSetD = [fr_Head_d]
    !     !eHeadSet = [er_Head]

    !     !fFlowSet = [fr_Flowrate]
    !     eFlowSet = [er_FlowrateLateral]

    !     do ii=1,size(eFlowSet)
    !         elemR(elem_P,eFlowSet(ii)) = BC%flowR(idx_P,br_value)
    !     end do
    !     !% For lateral flow, just update the flow at the element >> elemR(flow) + BC_lateral_flow

    !     if (setting%Debug%File%boundary_conditions) &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine face_interpolation_latBC_byPack
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine face_interpolation_dnBC_OLD(isBConly)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Interpolates data to all downstream boundary faces
    !     !% When called with isBConly == true, only does the BC update
    !     !%-------------------------------------------------------------------
    !     !% Declarations
    !         logical, intent(in) :: isBConly
    !         !integer :: fGeoSetU(3), fGeoSetD(3), eGeoSet(3)
    !         integer :: ii
    !         integer, pointer :: idx_fBC(:), eup(:), idx_P(:)
    !         integer, pointer :: elemUpstream
    !         real(8), pointer :: depthBC(:), headBC(:), eHead(:),  eFlow(:)
    !         real(8), pointer :: eVelocity(:), eZbottom(:), eLength(:)
    !         real(8), pointer :: eVolume(:), grav, eDepth(:), eTopwidth(:)
    !         logical, pointer :: isZeroDepth(:), hasFlapGateBC(:)
    !         !real(8), pointer :: eArea(:), eHydDepth(:), eTopWidth(:)
    !         real(8) :: Vtemp, headdif, thisDepth, thisVolume, thisQ
    !         character(64) :: subroutine_name = 'face_interpolation_dnBC'
    !     !%--------------------------------------------------------------------
    !     !% Preliminaries
    !         if (setting%Debug%File%boundary_conditions)  &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%--------------------------------------------------------------------
        
    !         eup           => faceI(:,fi_Melem_uL)
    !         idx_fBC       => faceP(1:npack_faceP(fp_BCdn),fp_BCdn)
    !         idx_P         => BC%P%BCdn
    !         depthBC       => BC%headR(:,br_Temp01)
    !         headBC        => BC%headR(:,br_value)
    !         hasFlapGateBC => BC%headYN(:,bYN_hasFlapGate) 
    !         eDepth        => elemR(:,er_Depth)
    !         eHead         => elemR(:,er_Head)
    !         eFlow         => elemR(:,er_Flowrate)
    !         eTopwidth     => elemR(:,er_TopWidth)
    !         eVelocity     => elemR(:,er_Velocity)
    !         eZbottom      => elemR(:,er_Zbottom)
    !         eVolume       => elemR(:,er_Volume)
    !         eLength       => elemR(:,er_Length)
    !         isZeroDepth   => elemYN(:,eYN_isZeroDepth)
    !         grav        => setting%Constant%gravity
    !         ! eArea     => elemR(:,er_Area)
    !         ! eHydDepth => elemR(:,er_HydDepth)
    !         ! eTopWidth => elemR(:,er_TopWidth)
    !     !%--------------------------------------------------------------------   

    !     !% For fixed, tidal, and timeseries BC
    !     !% The BC is imagined as enforced on a ghost cell outside the boundary
    !     !% so the face value is given by linear interpolation using ghost and interior cells 

    !     !% --- upstream face head is the BC
    !     !faceR(idx_fBC, fr_Head_u) = 0.5 * (eHead(eup(idx_fBC)) + headBC) 
    !     faceR(idx_fBC, fr_Head_u)  = headBC(idx_P)
    !     faceR(idx_fBC, fr_Depth_u) = max(headBC(idx_P) - faceR(idx_fBC,fr_Zbottom), setting%ZeroValue%Depth*0.99d0)

    !     print*, ' '
    !     print *, 'in ',trim(subroutine_name)
    !     print *, 'BC Head idx ', BC%headI(:, bi_idx)
    !     print *, 'BC head cat ', BC%headI(:, bi_subcategory)
    !     print *, 'BC head key ', trim(reverseKey(BC%headI(1, bi_subcategory))), ' ',trim(reverseKey(BC%headI(1, bi_subcategory)))
    !     print *, 'face ID ', idx_fBC
    !     print *, 'faceR head value ', faceR(idx_fBC, fr_Head_u)
    !     print *, 'BC head          ',headBC(idx_P)
    !     print *, 'BC depth         ',headBC(idx_P) - faceR(idx_fBC,fr_Zbottom)

    

    !     !% --- the downstream side of face is the same as the upstream face (unless gate, see below)
    !     faceR(idx_fBC, fr_Head_d)  = faceR(idx_fBC, fr_Head_u)
    !     faceR(idx_fBC, fr_Depth_d) = faceR(idx_fBC, fr_Depth_u)

    !     print *, 'e head, head bc ', eHead(eup(idx_fBC)), headBC(idx_P)
        
    !     !% --- for a flap gate on a BC with higher head downstream
    !     where ( hasFlapGateBC(idx_P) .and. (eHead(eup(idx_fBC))  < headBC(idx_P) ) )
    !         !% --- reset the head on the upstream and downstream side of face for closed gate
    !         faceR(idx_fBC, fr_Head_u)  = eHead(eup(idx_fBC))
    !         faceR(idx_fBC, fr_Depth_u) = max(eHead(eup(idx_fBC)) - faceR(idx_fBC,fr_Zbottom), setting%ZeroValue%Depth*0.99d0)
    !         faceR(idx_fBC, fr_Head_d)  = headBC(idx_P)
    !         faceR(idx_fBC, fr_Depth_d) = max(headBC(idx_P)       - faceR(idx_fBC,fr_Zbottom), setting%ZeroValue%Depth*0.99d0)
    !     endwhere

    !     print *, 'faceR value head ',faceR(idx_fBC, fr_Head_u), faceR(idx_fBC, fr_Head_d)

    !     ! ! ! ! ! call util_utest_CLprint ('    face interpolate down BC')
        
    !     !% --- get geometry for face from upstream element shape
    !     if (.not. isBConly) then
    !         !% --- get the depth on the face (upstream) from BC (temporary store)
    !         !%     this ensures that if gate is closed the depth is the upstream depth
    !         depthBC(idx_P) = faceR(idx_fBC, fr_Depth_u) !faceR(idx_fBC, fr_Head_u) - faceR(idx_fBC,fr_Zbottom)

    !         !% --- compute elldepth, topwidth, area geometry from depth based on relationship for upstream element
    !         !%     but using the depthBC at the upstream side of the face (which may be closed gate)   
    !         do ii=1,size(idx_fBC)
    !             elemUpstream => eup(idx_fBC(ii))
    !             !faceR(idx_fBC(ii),fr_HydDepth_u) = geo_hyddepth_from_depth_singular(elemUpstream,depthBC(idx_P(ii)))
    !             !faceR(idx_fBC(ii),fr_Length_u)   = elemR(elemupstream,er_Length)
    !             !faceR(idx_fBC(ii),fr_Topwidth_u) = geo_topwidth_from_depth_singular(elemUpstream,depthBC(idx_P(ii)))
    !             faceR(idx_fBC(ii),fr_Area_u)     = geo_area_from_depth_singular   &
    !                 (elemUpstream, depthBC(idx_P(ii)), setting%ZeroValue%Area)
    !             !faceR(idx_fBC(ii),fr_HydDepth_u) = geo_hyddepth_from_area_and_topwidth_singular(elemUpstream, faceR(idx_fBC(ii),fr_Area_u),faceR(idx_fBC(ii),fr_Topwidth_u) )
    !             ! faceR(idx_fBC(ii),fr_EllDepth_u) = geo_elldepth_singular &
    !             !     (faceR(idx_fBC, fr_Head_u), faceR(idx_fBC(ii),fr_Area_u), faceR(idx_fBC(ii),fr_Topwidth_u), &
    !             !      elemR(elemUpstream,er_AreaBelowBreadthMax), elemR(elemUpstream,er_BreadthMax), elemR(elemUpstream,er_ZbreadthMax))
                
    !             ! !% TEST 20220712--- apply simple linear interpolation to prevent large downstream area from causing numerical problems
    !             ! !% 20220712brh
    !             ! faceR(idx_fBC(ii),fr_Area_u)     =  (faceR(idx_fBC(ii),fr_Area_u)     + eArea(elemUpstream)    ) * onehalfR
    !             ! faceR(idx_fBC(ii),fr_HydDepth_u) =  (faceR(idx_fBC(ii),fr_HydDepth_u) + eHydDepth(elemUpstream)) * onehalfR
    !             ! faceR(idx_fBC(ii),fr_Topwidth_u) =  (faceR(idx_fBC(ii),fr_Topwidth_u) + eTopWidth(elemUpstream)) * onehalfR
    !         end do
    !         !% --- store downstream side of face
    !         !faceR(idx_fBC,fr_Topwidth_d) = faceR(idx_fBC,fr_Topwidth_u) 
    !         faceR(idx_fBC,fr_Area_d)     = faceR(idx_fBC,fr_Area_u) 
    !         !faceR(idx_fBC,fr_HydDepth_d) = faceR(idx_fBC,fr_HydDepth_u)

    !         print *, ' '
    !         print *, 'in face_interpolation_dnBC'
    !         print *, 'depthBC ',depthBC(idx_P(1))
    !         print *, 'area    ',faceR(idx_fBC(1),fr_Area_u)
    !         print *, ' '
    !         ! !stop 2987355

    !         !% --- set the flowrate at the boundary
    !         do ii = 1,size(idx_fBC)
    !             !% --- upstream element ID
    !             elemUpstream => eup(idx_fBC(ii))

    !             !% --- default face flowrate is the upstream element flowrate
    !             faceR(idx_fBC(ii),fr_Flowrate) = eFlow(elemUpstream)

    !             !% --- check for closed flap gate
    !             if (hasFlapGateBC(idx_P(ii))) then
    !                 ! print *, 'has flapgate'
    !                 if (faceR(idx_fBC(ii), fr_Head_u) < headBC(idx_P(ii))) then
    !                     !% --- set BC flow to zero for closed flap gate
    !                     !%     no outflow until head at face upstream exceeds the gate BC head
    !                     faceR(idx_fBC(ii), fr_Flowrate) = zeroR
    !                     ! print *, 'setting face flowrate to zero'
    !                 else
    !                     !% --- set BC flow to the flow at the upstream element center
    !                     !%     note that if the upstream value is negative the BC is set to zero
    !                     !%     to prevent backflow.
    !                     faceR(idx_fBC(ii), fr_Flowrate) = max(eFlow(elemUpstream),zeroR)
    !                     ! print *, 'face flowrate ',faceR(idx_fBC(ii), fr_Flowrate)
    !                 end if
    !             end if

    !             !% --- limit the flowrate based on volume to bring upstream element
    !             !%     to the same elevation as the head BC
    !             !%     if headdif > 0 then thisQ >0 represents a limit on the outflow
    !             !%     if headdif < 0 then thisQ <0 represents a limit on the inflow
    !             headdif = eHead(elemUpstream) - faceR(idx_fBC(ii), fr_Head_u)
    !             !% --- positive thisQ is an outflow
    !             thisQ   = onefourthR * headdif * eLength(elemUpstream)  &
    !                         * eTopwidth(elemUpstream) / setting%Time%Hydraulics%Dt

    !             print *, 'HeadDif     ', headdif
    !             print *, 'fr_Q, thisQ ', faceR(idx_fBC(1),fr_Flowrate), thisQ

    !             !% --- flow limiters
    !             if (faceR(idx_fBC(ii), fr_Flowrate) .ge. zeroR) then 
    !                 !% --- nominal outflow
    !                 if (thisQ > zeroR) then
    !                     !% --- limit outflow to the smaller value
    !                     faceR(idx_fBC(ii), fr_Flowrate) = min(faceR(idx_fBC(ii), fr_Flowrate), thisQ) 
    !                 else 
    !                     !% --- outflow across a boundary of higher head should not occur
    !                     faceR(idx_fBC(ii), fr_Flowrate) = zeroR        
    !                 end if 
    !             else 
    !                 !% --- nominal inflow (cannot occur with flap gate)
    !                 if (thisQ > zeroR)  then 
    !                     !% --- inflow should not occur if thisQ > 0
    !                     faceR(idx_fBC(ii), fr_Flowrate) = zeroR  
    !                 else
    !                     !% headdif < 0 and thisQ < 0 implies backflow
    !                     if (isZeroDepth(elemUpstream))then
    !                         !% --- for zero depth, limit headdif by 50% of the face depth as the driving head
    !                         headdif = -min(-headdif, onehalfR * faceR(idx_fBC(ii),fr_Depth_d))
    !                         thisQ   = onefourthR * headdif * eLength(elemUpstream)  &
    !                                     * eTopwidth(elemUpstream) / setting%Time%Hydraulics%Dt
    !                     end if
    !                     if (eFlow(elemUpstream) < zeroR) then
    !                         !% --- nominal inflow over BC
    !                         !%     use the smaller magnitude (more positive value)
    !                         faceR(idx_fBC(ii), fr_Flowrate) = max(thisQ,eFlow(elemUpstream))
    !                     else
    !                         !% -- - use the filling value
    !                         faceR(idx_fBC(ii), fr_Flowrate) = thisQ
    !                     end if
    !                 end if
    !             end if

    !             print *, 'flowrate DD',faceR(idx_fBC(ii),fr_Flowrate)
    !         end do

    !             !% --- handle possible backflows---------------------
                
    !             ! !% ---compute adverse head difference (from downstreamface to upstream element center)
    !             ! !%    that drives could drive an inflow on a downstream head BC        
    !             ! headdif = faceR(idx_fBC(ii), fr_Head_u) - eHead(elemUpstream)
    !             ! if (isZeroDepth(elemUpstream))then
    !             !     !% --- for zero depth, limit headdif by 50% of the face depth as the driving head
    !             !     !headdif = min(headdif, onehalfR * faceR(idx_fBC(ii),fr_HydDepth_d))
    !             !     headdif = min(headdif, onehalfR * faceR(idx_fBC(ii),fr_Depth_d))
    !             ! end if

    !             !  print *, 'ii, headdif ',ii, headdif
    !             !  print *, faceR(idx_fBC(ii), fr_Head_u), eHead(elemUpstream), eZbottom(elemUpstream)

    !             !% --- check for adverse head gradient going upstream from boundary
    !             !%     when element flow is downstream or zero
    !             !      note that when element velocity is already negative we do not add the headdif term
    !         !     if (headdif > zeroR) then 
    !         !         !% --- adverse head gradient at outlet
    !         !         if (eFlow(elemUpstream) .ge. zeroR) then
    !         !             !% --- potential for inflow from downstream boundary despite downstream
    !         !             !%     flow on element.
    !         !             !%     Piezometric head difference minus upstream velocity
    !         !             !%     head provides the velocity head of inflow
    !         !             !%     If positive, then this is the upstream velocity at face
    !         !             !%     If negative, then the downstream flow in element is
    !         !             !%     able to overcome the adverse piezometric head gradient,
    !         !             !%     so the difference provides a reduced outflow at face
    !         !             !% --- The following is the estimated Bernoulli velocity squared. Note that
    !         !             !%     this assumes Velocity sign is consistent with flowrate sign.
    !         !             !Vtemp = twoR * grav * headdif - (eVelocity(elemUpstream)**twoI)
    !         !             !%     Use the velocity head if the upstream flow is distributed over the downstream area
    !         !             Vtemp = twoR * grav * headdif - ((eFlow(elemUpstream) / faceR(idx_fBC(ii),fr_Area_u))**twoI)
    !         !             !% --- The +Vtemp is flow in the upstream (negative) direction
    !         !             Vtemp = -sign(sqrt(abs(Vtemp)),Vtemp)
    !         !             !% --- take the smaller of the Q implied by Vtemp and the downstream Q of the element
    !         !             !%     Note Vtemp < 0 is upstream flow and will automatically be selected here since the upstream
    !         !             !%     elem flowrate is guaranteed to be positive. The following statement simply
    !         !             !%     ensures that the head balance approach for a downstream flow Vtemp > 0 does not exceed existing
    !         !             !%     downstream flow on the element imlied by the (flowrate >= 0) conditional above.
    !         !             !%     NOTE this flowrate may be either upstream or downstream!
    !         !             faceR(idx_fBC(ii),fr_Flowrate) = min(Vtemp * faceR(idx_fBC(ii),fr_Area_u),eFlow(elemUpstream))

    !         !             !  print *, 'naive face flowrate ', faceR(idx_fBC(ii),fr_Flowrate)
    !         !         end if

    !         !         !% --- Limit inflow magnitude to that the brings level up to the BC elevation
    !         !         if (faceR(idx_fBC(ii),fr_Flowrate) < zeroR) then
    !         !             !% --- volume rate adjustment: negative flowrate cannot provide more volume than 
    !         !             !%     fills the pipe to the depth equivalent to the downstream BC head in a single time step.
    !         !             !%     This reduces oscillatory behavior by preventing over filling on back flow
    !         !             thisDepth  = faceR(idx_fBC(ii), fr_Head_u) - eZbottom(elemUpstream)

    !         !             ! print *, 'thisDepth ',thisDepth, elemR(elemUpstream,er_Depth)
    !         !             ! print *, 'this head ',faceR(idx_fBC(ii), fr_Head_u), elemR(elemUpstream,er_Head)
    !         !             ! print *, ' '
    !         !             ! print *, 'elemUpstream ',elemUpstream
    !         !             ! print *, 'head,z ', elemR(elemUpstream,er_Head), elemR(elemUpstream,er_Zbottom)
    !         !             ! print *, 'depth  ', elemR(elemUpstream,er_Head) - elemR(elemUpstream,er_Zbottom), elemR(elemUpstream,er_Depth)

    !         !             if (thisDepth > zeroR) then
    !         !                 ! !% --- get the volume of the upstream element if filled to the BC head level
    !         !                 ! thisVolume = geo_area_from_depth_singular(elemUpstream, thisDepth, setting%ZeroValue%Area) * eLength(elemUpstream)
    !         !                 ! !% --- flowrate to fill to this volume in one time step
    !         !                 ! thisQ      = (thisVolume - eVolume(elemUpstream)) / setting%Time%Hydraulics%Dt

    !         !                 ! print *, 'thisQ ',thisQ
    !         !                 ! print *, thisVolume, eVolume(elemUpstream), eVolume(73)

    !         !                 !% --- using volume caused problem with tables not being precisely invertable
    !         !                 !%     from depth -> area  and from volume -> depth -> area
    !         !                 !% --- approximate the volume flowrate into the upstream element
    !         !                 !%     by the volume associated with the depth difference
    !         !                 thisQ = onehalfR * (thisDepth - eDepth(elemUpstream)) * eLength(elemUpstream)  &
    !         !                         * eTopwidth(elemUpstream) / setting%Time%Hydraulics%Dt

    !         !                 if (thisQ > zeroR) then
    !         !                     !% --- inflow is allowed, but use the smaller magnitude flow (larger negative value)
    !         !                     faceR(idx_fBC(ii),fr_Flowrate) = max(faceR(idx_fBC(ii),fr_Flowrate),-thisQ)
    !         !                 else
    !         !                     !% --- if thisVolume <= 0, then upstream is already to the maximum level
    !         !                     faceR(idx_fBC(ii),fr_Flowrate) = zeroR
    !         !                 end if
                            
    !         !                 ! print *, 'face flowrate with depth ',faceR(idx_fBC(ii),fr_Flowrate)
    !         !             else
    !         !                 !% --- negative upstream depth implies Zbottom of upstream element is higher
    !         !                 !%     than the downstream BC.  We shouldn't get here because it should only 
    !         !                 !%     occur if headdif < 0, 
    !         !                 print *, 'CODE ERROR -- unexpected else condition reached'
    !         !                 call util_crashpoint(5593341)
    !         !             end if
    !         !         else
    !         !             !% --- if flowrate is > 0, no volume rate adjustment needed and we accept 
    !         !             !%     previously computed value
    !         !         end if
    !         !     else
    !         !         !% --- limit inflow to that which brings upstream element up to the head of the BC
    !         !         thisDepth  = faceR(idx_fBC(ii), fr_Head_u) - eZbottom(elemUpstream)
    !         !         if (thisDepth > zeroR) then
    !         !         else 
    !         !             !% --- negative upstream depth implies Zbottom of upstream element is higher
    !         !                 !%     than the downstream BC.  We shouldn't get here because it should only 
    !         !                 !%     occur if headdif < 0, 
    !         !             print *, 'CODE ERROR -- unexpected else condition reached'
    !         !                 call util_crashpoint(5593341)
    !         !         end if

    !         !         !% --- if there is no adverse head gradient or if element velocity is already negative
    !         !         !%     in/outflow rate is upstream element flowrate                        
    !         !         faceR(idx_fBC(ii), fr_Flowrate) = eFlow(elemUpstream)
    !         !     end if

    !         !     ! print *, 'final face flowrate ', faceR(idx_fBC(ii),fr_Flowrate)

    !         ! end do

    !         ! ! ! ! ! call util_utest_CLprint ('    face YYYY')

    !         !% --- set the Preissmann number to the upstream element value
    !         faceR(idx_fBC, fr_Preissmann_Number) = elemR(eup(idx_fBC), er_Preissmann_Number) 

    !         ! !% --- set to zero flow for closed gate  !% MOVED UP INTO DO LOOP 20220716brh
    !         ! where ( BC%headYN(idx_P, bYN_hasFlapGate) .and. (faceR(idx_fBC, fr_Head_u) < BC%headR(idx_P,br_value)) )
    !         !     faceR(idx_fBC, fr_Flowrate) = zeroR
    !         ! end where

    !         ! print *, 'face area u ', faceR(idx_fBC,fr_Area_u)
    !         ! print *, 'face area d ', faceR(idx_fBC,fr_Area_d)

    !         !% --- ensure face area_u is not smaller than zerovalue
    !         where (faceR(idx_fBC,fr_Area_d) < setting%ZeroValue%Area)
    !             faceR(idx_fBC,fr_Area_d) = setting%ZeroValue%Area
    !         endwhere
    !         where (faceR(idx_fBC,fr_Area_u) < setting%ZeroValue%Area)
    !             faceR(idx_fBC,fr_Area_u) = setting%ZeroValue%Area
    !         endwhere

    !         faceR(idx_fBC,fr_Velocity_u) = faceR(idx_fBC,fr_Flowrate)/faceR(idx_fBC,fr_Area_u)
    !         faceR(idx_fBC,fr_Velocity_d) = faceR(idx_fBC,fr_Flowrate)/faceR(idx_fBC,fr_Area_d)  

    !         print *, 'face velocity u ', faceR(idx_fBC,fr_Velocity_u)
    !         print *, 'face velocity d ', faceR(idx_fBC,fr_Velocity_d)

    !         !%  limit high velocities
    !         if (setting%Limiter%Velocity%UseLimitMaxYN) then
    !             where(abs(faceR(idx_fBC,fr_Velocity_u))  > setting%Limiter%Velocity%Maximum)
    !                 faceR(idx_fBC,fr_Velocity_u) = sign(0.99d0 * setting%Limiter%Velocity%Maximum, &
    !                     faceR(idx_fBC,fr_Velocity_u))
    !             endwhere 
    !             where(abs(faceR(idx_fBC,fr_Velocity_d))  > setting%Limiter%Velocity%Maximum)
    !                 faceR(idx_fBC,fr_Velocity_d) = sign(0.99d0 * setting%Limiter%Velocity%Maximum, &
    !                     faceR(idx_fBC,fr_Velocity_d))
    !             end where                    
    !         end if

    !         ! ! ! ! ! call util_utest_CLprint ('    face ZZZZ')
    !     else
    !         !% continue
    !     end if

    !     !%--------------------------------------------------------------------
    !     !% Closing
    !         if (setting%Debug%File%boundary_conditions) &
    !             write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine face_interpolation_dnBC_OLD
!%
!%==========================================================================
    !%==========================================================================
!%
    ! subroutine face_interp_shared_set_old &
    !     (fset, eset, eWdn, eWup, facePackCol, Npack)
    !     !%-------------------------------------------------------------------
    !     !% Description:
    !     !% Interpolates faces shared between processor
    !     !%-------------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in) :: fset(:), eset(:), eWdn, eWup, facePackCol, Npack
    !         integer, pointer :: thisP, eup, edn, connected_image, ghostUp, ghostDn
    !         logical, pointer :: isGhostUp, isGhostDn
    !         integer :: ii, jj   
    !         integer(kind=8) :: crate, cmax, cval
    !         character(64) :: subroutine_name = 'face_interp_shared_set_old'
    !     !%--------------------------------------------------------------------
    !     !%  Preliminaries
    !         if (setting%Debug%File%face) &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !         sync all
    !         if (this_image()==1) then
    !             call system_clock(count=cval,count_rate=crate,count_max=cmax)
    !             setting%Time%WallClock%SharedStart_A = cval
    !         end if    
    !     !%--------------------------------------------------------------------
    !     !% cycle through all the shared faces
    !     do ii = 1,Npack
    !         !%-----------------------------------------------------------------
    !         !% Aliases
    !         thisP           => facePS(ii,facePackCol)
    !         connected_image => faceI(thisP,fi_Connected_image)
    !         eup             => faceI(thisP,fi_Melem_uL)
    !         edn             => faceI(thisP,fi_Melem_dL)
    !         ghostUp         => faceI(thisP,fi_GhostElem_uL)
    !         ghostDn         => faceI(thisP,fi_GhostElem_dL)
    !         isGhostUp       => faceYN(thisP,fYN_isUpGhost)
    !         isGhostDn       => faceYN(thisP,fYN_isDnGhost)
    !         !%-----------------------------------------------------------------
    !         !% cycle through each element in the set.
    !         !% This is designed for fset and eset being vectors, but it
    !         !%   is not clear that this is needed.
    !         do jj=1,size(fset)

    !             !% condition for upstream element of the shared face is ghost and in a different image
    !             if (isGhostUp) then

    !                 faceR(thisP,fset(jj)) = &
    !                     (+elemR(ghostUp,eset(jj))[connected_image] * elemR(edn,eWup) &
    !                      +elemR(edn,eset(jj)) * elemR(ghostUp,eWdn)[connected_image] &
    !                     ) / &
    !                     ( elemR(edn,eWup) + elemR(ghostUp,eWdn)[connected_image] )

    !             !% condition for downstream element of the shared face is ghost and in a different image
    !             elseif (isGhostDn) then

    !                 faceR(thisP,fset(jj)) = &
    !                     (+elemR(eup,eset(jj)) * elemR(ghostDn,eWup)[connected_image] &
    !                      +elemR(ghostDn,eset(jj))[connected_image] * elemR(eup,eWdn) &
    !                     ) / &
    !                     ( elemR(ghostDn,eWup)[connected_image] + elemR(eup,eWdn) )

    !             else
    !                 write(*,*) 'CODE ERROR: unexpected else'
    !                 !stop 
    !                 call util_crashpoint( 487874)
    !                 !return
    !             end if        
    !         end do
    !     end do

    !     !% NOTES
    !     !% elemR(eup,eset(jj)) is the element value upstream of the face
    !     !% elemR(edn,eset(jj) is the element value downstream of the face.
    !     !% elemR(eup,eWdn) is the downstream weighting of the upstream element
    !     !% elemR(edn,eWup)) is the upstream weighting of the downstream element

    !     !% elemR(ghostUp,eset(jj))[connected_image] is the elem value from the upstream image of the face
    !     !% elemR(ghostDn,eset(jj))[connected_image] is the elem value from the downstream image of the face
    !     !% elemR(ghostUp,eWdn)[connected_image] is the downstream weighting of the upstream image element
    !     !% elemR(ghostDn,eWup))[connected_image] is the upstream weighting of the downstream image element

    !     !%--------------------------------------------------------------------
    !     !% Closing
    !     sync all
    !     if (this_image()==1) then
    !         !% stop the shared timer
    !         call system_clock(count=cval,count_rate=crate,count_max=cmax)
    !         setting%Time%WallClock%SharedStop_A = cval
    !         setting%Time%WallClock%SharedCumulative_A &
    !                 = setting%Time%WallClock%SharedCumulative_A &
    !                 + setting%Time%WallClock%SharedStop_A &
    !                 - setting%Time%WallClock%SharedStart_A                    
    !     end if 
    !         if (setting%Debug%File%face) &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine face_interp_shared_set_old
!%
!%==========================================================================
!%==========================================================================
!%  
    ! subroutine face_head_average_on_element &
    !     (whichTM)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% computes the average head of the faces on an element
    !     !%-------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: whichTM 
    !         integer, pointer :: Npack, thisP(:), thisCol
    !         integer, pointer :: mapUp(:), mapDn(:)
    !         real(8), pointer :: fHeadU(:), fHeadD(:), eHeadAvg(:)
    !         character(64) :: subroutine_name = 'face_head_average_on_element'
    !     !%-------------------------------------------------------------------
    !     !% Preliminaries
    !         !if (.not. setting%Solver%QinterpWithLocalHeadGradient) return  
    !         if (setting%Debug%File%face) &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !         select case (whichTM)
    !             case (ALLtm)
    !                 thisCol => col_elemP(ep_CC_ALLtm)
    !             case (ETM)
    !                 thisCol => col_elemP(ep_CC_ETM)
    !             case (AC)
    !                 thisCol => col_elemP(ep_CC_AC)
    !             case (dummy)
    !                 ! print *, 'error - the FluxCorrection has not been coded for diagnostic elements'
    !                 ! stop 58704
    !                 return
    !             case default
    !                 print *, 'error, this default case should not be reached'
    !                 stop 2394
    !         end select         
    !     !%-------------------------------------------------------------------
    !     !% Aliases             
    !         Npack    => npack_elemP(thisCol)
    !         if (Npack .le. 0) return
    !         thisP    => elemP(1:Npack,thisCol)
    !         mapUp    => elemI(:,ei_Mface_uL)
    !         mapDn    => elemI(:,ei_Mface_dL)   
    !         fHeadU   => faceR(:,fr_Head_u)  
    !         fHeadD   => faceR(:,fr_Head_d)
    !         eHeadAvg => elemR(:,er_HeadAvg)    
    !     !%-------------------------------------------------------------------   
    !     !% The map up must use the downstream head on the face.
    !     !% The map dn must use the upstream head on the face.
    !     eHeadAvg(thisP) = onehalfR * (fHeadU(mapDn(thisP)) + fHeadD(mapUp(thisP)))

    !     ! print *, 'in ',trim(subroutine_name)
    !     ! print *, thisP
    !     ! print *, eHeadAvg(thisP)
    !     ! print *, ' '
    !     ! print *, eHeadAvg(:)
    !     ! print *, ' '
    !     ! print *, fHeadU(:)
    !     ! print *, ' '
    !     ! print *, fHeadD(:)
    !     ! print *, ' '
    !     !%-------------------------------------------------------------------
    !     !% Closing
    !         if (setting%Debug%File%face) &
    !             write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine face_head_average_on_element
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine face_FluxCorrection_interior &
        !     (faceCol, whichTM)
        !     !%------------------------------------------------------------------
        !     !% Description:
        !     !% Adds the head gradient term to the face flowrate for interior faces
        !     !% should be done after Q and H are interpolated to face
        !     !% and element HeadAvg is computed.
        !     !%-------------------------------------------------------------------
        !     !% Declarations
        !         integer, intent(in) :: faceCol, whichTM
        !         integer, pointer :: Npack, thisF(:), eup(:), edn(:), elist(:)
        !         real(8), pointer :: fQ(:), eArea(:), eHead(:), eHeadAvg(:)
        !         real(8), pointer :: eLength(:) !, qLateral(:), qChannel(:)
        !         real(8), pointer ::  dt, grav !, qfac, qratio 
        !         logical          :: isBConly
        !         character(64) :: subroutine_name = 'face_FluxCorrection_interior'
        !     !%-------------------------------------------------------------------
        !     !% Preliminaries
        !         !if (.not. setting%Solver%QinterpWithLocalHeadGradient) return  
        !         if (setting%Debug%File%face) &
        !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !     !%-------------------------------------------------------------------
        !     !% Aliases   
        !         eup      => faceI(:,fi_Melem_uL)
        !         edn      => faceI(:,fi_Melem_dL)
        !         fQ       => faceR(:,fr_Flowrate)
        !         eArea    => elemR(:,er_Area)
        !         eHead    => elemR(:,er_Head)
        !         !eHeadAvg => elemR(:,er_HeadAvg)
        !         eLength  => elemR(:,er_Length)
        !         eList    => elemI(:,ei_Temp01)
        !         dt       => setting%Time%Hydraulics%Dt
        !         grav     => setting%constant%gravity

        !         Npack => npack_faceP(faceCol)
        !         if (Npack .le. 0) return
        !         thisF    => faceP(1:Npack,faceCol)
        !     !%-------------------------------------------------------------------
        !     !% --- compute the average head for the elements 
        !     !call face_head_average_on_element (whichTM)

        !     !% -- we need a custom selector array because we don't have a packed array
        !     !%    that handles the combined face/element condition needed
        !     eList = zeroI
        !     select case (whichTM)
        !     case (ALLtm)
        !         eList(ep_CC_ALLtm) = oneI
        !     case (ETM)
        !         eList(ep_CC_ETM) = oneI
        !     case (AC)
        !         eList(ep_CC_AC) = oneI
        !     case (dummy)
        !         ! print *, 'error - the FluxCorrection has not been coded for diagnostic elements'
        !         return
        !     case default
        !         print *, 'error, this default case should not be reached'
        !         stop 239483
        !     end select

        !     !% --- adds term dt * grav A [ (dh/dx) - (dh_avg/dx) ] where not zerovolume
        !     where (      (.not. elemYN(eup(thisF),eYN_isZeroDepth   )        ) & 
        !            .and. (.not. elemYN(eup(thisF),eYN_isSmallDepth  )        ) &
        !            .and. (       elemR(eup(thisF),er_FlowrateLateral) > zeroR) &
        !            .and. (       eList(eup(thisF))                    == oneI) )
        !         fQ(thisF) = fQ(thisF) + dt * grav *                                         &
        !             (                                                                       &
        !                 +( eArea(eup(thisF)) * ( eHead(eup(thisF)) - eHeadAvg(eup(thisF)) ) &
        !                     / ( onehalfR * eLength(eup(thisF) ) ) )                         &
        !             )                        
        !     end where

        !     where (       (.not. elemYN(edn(thisF),eYN_isZeroDepth   )        ) & 
        !             .and. (.not. elemYN(edn(thisF),eYN_isSmallDepth  )        ) &
        !             .and. (       elemR(edn(thisF),er_FlowrateLateral) > zeroR) &
        !             .and. (       eList(edn(thisF))                    == oneI) )
        !         fQ(thisF) = fQ(thisF) + dt * grav *                                         &
        !             (                                                                       &
        !                 -( eArea(edn(thisF)) * ( eHead(edn(thisF)) - eHeadAvg(edn(thisF)) ) &
        !                     / ( onehalfR * eLength(edn(thisF) ) ) )                         &
        !             ) 
        !     end where

        !     !% --- need another call to face_interpolate so that the Q_HeadGradient
        !         !%     does not change the upper boundary inflow condition
        !     isBConly = .true.
        !     call face_interpolate_bc (isBConly)

        !     !% for lateral inflows upstream of a face with downstream flow 
        !     !% note: null set for negative inflow   
        !     ! where ( qLateral(eup(thisF)) > qratio * abs(qChannel(eup(thisF))) )
        !     !     fQ(thisF) = fQ(thisF) + qfac * dt * grav                                       &
        !     !         * ( util_sign_with_ones(fQ(thisF)) + oneR ) * onehalfR                     &
        !     !         *(                                                                         &
        !     !             +( eArea(eup(thisF)) * ( eHead(eup(thisF)) - eHeadAvg(eup(thisF)) ) )  &
        !     !             / ( onehalfR * eLength(eup(thisF)) )                                   &
        !     !         )
        !     ! end where

        !     ! !% for lateral inflows downstream of a face with an upstream flow
        !     ! !% note: null set for negative inflow
        !     ! where ( qLateral(edn(thisF)) > qratio * abs(qChannel(edn(thisF))) )
        !     !     fQ(thisF) = fQ(thisF) + qfac * dt * grav                                       &
        !     !         * ( util_sign_with_ones(fQ(thisF)) - oneR ) * onehalfR                     &
        !     !         *(                                                                         &   
        !     !           +( eArea(edn(thisF)) * ( eHead(edn(thisF)) - eHeadAvg(edn(thisF)) ) )    &
        !     !             / ( onehalfR * eLength(edn(thisF)) )                                   &
        !     !         )
        !     ! end where

        !     !% --- for downstream flow
        !     ! where (.not. elemYN(eup(thisP),eYN_isZeroDepth))
        !     !     fQ(thisP) = fQ(thisP) + qfac * dt * grav                                    &
        !     !        *( util_sign_with_ones(fQ(thisP)) + oneR ) * onehalfR                    &
        !     !        *(                                                                       &
        !     !             +( eArea(eup(thisP)) * ( eHead(eup(thisP)) - eHeadAvg(eup(thisP)) ) &
        !     !                 / ( onehalfR * eLength(eup(thisP) ) ) )                         &
        !     !         ) 
        !     ! end where

        !     ! !% --- for upstream flow
        !     ! where (.not. elemYN(edn(thisP),eYN_isZeroDepth))
        !     !     fQ(thisP) = fQ(thisP) + qfac * dt * grav                                    &
        !     !        *( util_sign_with_ones(fQ(thisP)) - oneR ) * onehalfR                    &
        !     !        *(                                                                       &
        !     !             -( eArea(edn(thisP)) * ( eHead(edn(thisP)) - eHeadAvg(edn(thisP)) ) &
        !     !                 / ( onehalfR * eLength(edn(thisP) ) ) )                         &
        !     !         ) 
        !     ! end where
        
        !     !%-------------------------------------------------------------------
        !     !% Closing
        !         if (setting%Debug%File%face) &
        !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        ! end subroutine face_FluxCorrection_interior
!%
!%==========================================================================
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module face
