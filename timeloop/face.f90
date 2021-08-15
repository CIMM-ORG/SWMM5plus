module face

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use adjust
    use jump

    implicit none

    !%-----------------------------------------------------------------------------
    !% Description:
    !% Provides computation of face values for timeloop of hydraulics
    !%
    !% METHOD:
    !%
    !%

    private

    public :: face_interpolation
    public :: face_interpolate_bc

    contains
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    !%
    subroutine face_interpolation (facecol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Interpolates faces
        !%-----------------------------------------------------------------------------
        integer, intent(in)  :: faceCol
        integer, pointer :: Npack
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'face_interpolation'
        if (setting%Debug%File%face) print *, '*** enter ', this_image(), subroutine_name
        !%-----------------------------------------------------------------------------

        !% face reconstruction of all the interior faces
        Npack => npack_faceP(faceCol)
        if (Npack > 0) then
            call face_interpolation_interior_byPack (faceCol, Npack)
        endif

        sync all

        !% face reconstruction of all the shared faces
        Npack => npack_facePS(faceCol)
        if (Npack > 0) then
            call face_interpolation_shared_byPack (faceCol, Npack)
        endif

        !% wait for all the processors to finish face interpolation across images
        sync all

        call face_interpolate_bc ()

        if (setting%Debug%File%face)  print *, '*** leave ', this_image(), subroutine_name
    end subroutine face_interpolation
    !%
    !%==========================================================================
    !% PRIVATE
    !%==========================================================================
    !%
    subroutine face_interpolate_bc()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Upper level BC interpolation face data population
        !%-----------------------------------------------------------------------------
        character (64) :: subroutine_name = 'face_interpolate_bc'
        !%-----------------------------------------------------------------------------

        if (setting%Debug%File%face)  print *, '*** enter ', this_image(), subroutine_name

        if (N_nBCup > 0) call face_interpolation_upBC_byPack()

        if (N_nBClat > 0) call face_interpolation_latBC_byPack()

        if (N_nBCdn > 0) call face_interpolation_dnBC_byPack()

        if (setting%Debug%File%face) print *, '*** leave ', this_image(), subroutine_name

    end subroutine face_interpolate_bc
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine face_interpolation_upBC_byPack()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Interpolates all boundary faces using a pack arrays -- base on bi_category
        !%-----------------------------------------------------------------------------

        integer :: fGeoSetU(3), fGeoSetD(3), eGeoSet(3)
        integer :: fFlowSet(1), eFlowSet(1)
        integer :: fHeadSetU(1), fHeadSetD(1), eHeadSet(1)
        character(64) :: subroutine_name = 'face_interpolation_upBC_byPack'
        integer :: ii
        integer, pointer :: face_P(:), edn(:), idx_P(:)

        !%-----------------------------------------------------------------------------

        if (setting%Debug%File%boundary_conditions)  print *, '*** enter ', this_image(), subroutine_name

        !% For the head/geometry at the upstream faces, we directly take the dnwnstream element
        !% So there is no eup for upstream BCs
        edn => faceI(:,fi_Melem_dL)

        face_P => faceP(1:npack_faceP(fp_BCup),fp_BCup)
        idx_P  => BC%P%BCup(:)

        fGeoSetU = [fr_Area_u, fr_Topwidth_u, fr_HydDepth_u]
        fGeoSetD = [fr_Area_d, fr_Topwidth_d, fr_HydDepth_d]
        eGeoSet  = [er_Area,   er_Topwidth,   er_HydDepth]

        fHeadSetU = [fr_Head_u]
        fHeadSetD = [fr_Head_d]
        eHeadSet = [er_Head]

        fFlowSet = [fr_Flowrate]
        eFlowSet = [er_Flowrate]

        do ii=1,size(fFlowSet)
            faceR(face_P, fFlowSet(ii)) = BC%flowRI(idx_P)
        end do

        !% Copying other data to the BC faces
        do ii = 1,size(fGeoSetU)
            faceR(face_P,fGeoSetD(ii)) = elemR(edn(face_P),eGeoSet(ii)) ! Copying the Geo
            faceR(face_P,fGeoSetU(ii)) = faceR(face_P,fGeoSetD(ii))
        end do

        faceR(face_P,fr_Head_d) = elemR(edn(face_P),er_HydDepth) + faceR(face_P,fr_Zbottom)!Copying the Head
        faceR(face_P,fr_Head_u) = faceR(face_P,fr_Head_d)
        
        !% HACK: This is needed to be revisited later
        if (setting%ZeroValue%UseZeroValues) then
            !% ensure face area_u is not smaller than zerovalue
            where (faceR(face_P,fr_Area_d) <= setting%ZeroValue%Area)
                faceR(face_P,fr_Area_d)     = setting%ZeroValue%Area
                faceR(face_P,fr_Area_u)     = setting%ZeroValue%Area
                faceR(face_P,fr_Velocity_d) = zeroR
                faceR(face_P,fr_Velocity_u) = zeroR
            endwhere

            where (faceR(face_P,fr_Area_d) > setting%ZeroValue%Area)
                faceR(face_P,fr_Velocity_d) = faceR(face_P,fr_Flowrate)/faceR(face_P,fr_Area_d)
                faceR(face_P,fr_Velocity_u) = faceR(face_P,fr_Velocity_d)  
            endwhere
        else
            !% ensure face area_u is not smaller than zerovalue
            where (faceR(face_P,fr_Area_d) <= zeroR)
                faceR(face_P,fr_Area_d)     = zeroR
                faceR(face_P,fr_Area_u)     = zeroR
                faceR(face_P,fr_Velocity_d) = zeroR
                faceR(face_P,fr_Velocity_u) = zeroR
            endwhere

            where (faceR(face_P,fr_Area_d) > zeroR)
                faceR(face_P,fr_Velocity_d) = faceR(face_P,fr_Flowrate)/faceR(face_P,fr_Area_d)
                faceR(face_P,fr_Velocity_u) = faceR(face_P,fr_Velocity_d)
            endwhere
        endif

        !%  limit high velocities
        if (setting%Limiter%Velocity%UseLimitMax) then
            where(abs(faceR(face_P,fr_Velocity_d))  > setting%Limiter%Velocity%Maximum)
                faceR(face_P,fr_Velocity_d) = sign(0.99 * setting%Limiter%Velocity%Maximum, &
                    faceR(face_P,fr_Velocity_d))

                faceR(face_P,fr_Velocity_u) = faceR(face_P,fr_Velocity_d)
            endwhere
        endif

        if (setting%Debug%File%boundary_conditions) print *, '*** leave ', this_image(), subroutine_name
    end subroutine face_interpolation_upBC_byPack
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine face_interpolation_latBC_byPack()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Interpolates all boundary faces using a pack arrays -- base on bi_category
        !%-----------------------------------------------------------------------------
        integer :: fGeoSetU(3), fGeoSetD(3), eGeoSet(3)
        integer :: fFlowSet(1), eFlowSet(1)
        integer :: fHeadSetU(1), fHeadSetD(1), eHeadSet(1)
        character(64) :: subroutine_name = 'face_interpolation_latBC_byPack'
        integer :: ii
        integer, pointer :: elem_P(:), idx_P(:)

        !%-----------------------------------------------------------------------------

        if (setting%Debug%File%boundary_conditions)  print *, '*** enter ', this_image(), subroutine_name


        elem_P => elemP(1:npack_elemP(ep_BClat),ep_BClat)
        idx_P  => BC%P%BClat

        fGeoSetU = [fr_Area_u, fr_Topwidth_u, fr_HydDepth_u]
        fGeoSetD = [fr_Area_d, fr_Topwidth_d, fr_HydDepth_d]
        eGeoSet  = [er_Area,   er_Topwidth,   er_HydDepth]

        fHeadSetU = [fr_Head_u]
        fHeadSetD = [fr_Head_d]
        eHeadSet = [er_Head]

        fFlowSet = [fr_Flowrate]
        eFlowSet = [er_FlowrateLateral]

        do ii=1,size(eFlowSet)
            elemR(elem_P,eFlowSet(ii)) = BC%flowRI(idx_P)
        end do
        !% For lateral flow, just update the flow at the element >> elemR(flow) + BC_lateral_flow

        if (setting%Debug%File%boundary_conditions) print *, '*** leave ', this_image(), subroutine_name
    end subroutine face_interpolation_latBC_byPack
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine face_interpolation_dnBC_byPack()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Interpolates all boundary faces using a pack arrays -- base on bi_category
        !%-----------------------------------------------------------------------------
        integer :: fGeoSetU(3), fGeoSetD(3), eGeoSet(3)
        integer :: fFlowSet(1), eFlowSet(1)
        integer :: fHeadSetU(1), fHeadSetD(1), eHeadSet(1)
        character(64) :: subroutine_name = 'face_interpolation_dnBC_byPack'
        integer :: ii
        integer, pointer :: face_P(:), eup(:), idx_P(:)

        !%-----------------------------------------------------------------------------

        if (setting%Debug%File%boundary_conditions)  print *, '*** enter ', this_image(), subroutine_name


        eup => faceI(:,fi_Melem_uL)

        face_P => faceP(1:npack_faceP(fp_BCdn),fp_BCdn)
        idx_P  => BC%P%BCdn

        fGeoSetU = [fr_Area_u, fr_Topwidth_u, fr_HydDepth_u]
        fGeoSetD = [fr_Area_d, fr_Topwidth_d, fr_HydDepth_d]
        eGeoSet  = [er_Area,   er_Topwidth,   er_HydDepth]

        fHeadSetU = [fr_Head_u]
        fHeadSetD = [fr_Head_d]
        eHeadSet = [er_Head]

        fFlowSet = [fr_Flowrate]
        eFlowSet = [er_Flowrate]


        do ii=1,size(fHeadSetD)
            !%  linear interpolation using ghost and interior cells
            faceR(face_P, fHeadSetU(ii)) = 0.5 * (elemR(eup(face_P), er_Head) + BC%headRI(idx_P)) !% downstream head update
            faceR(face_P, fHeadSetD(ii)) = faceR(face_P, fHeadSetU(ii))
        end do

        do ii=1,size(fFlowSet)
            faceR(face_P, fFlowSet(ii)) = elemR(eup(face_P), eFlowSet(ii)) !% Copying the flow from the upstream element
        end do

        do ii=1,size(fGeoSetD)
            faceR(face_P, fGeoSetD(ii)) = elemR(eup(face_P), eGeoSet(ii)) !% Copying other geo factors from the upstream element
            faceR(face_P, fGeoSetU(ii)) = faceR(face_P, fGeoSetD(ii))
        end do

        !% HACK: This is needed to be revisited later
        if (setting%ZeroValue%UseZeroValues) then
            !% ensure face area_u is not smaller than zerovalue
            where (faceR(face_P,fr_Area_d) < setting%ZeroValue%Area)
                faceR(face_P,fr_Area_d) = setting%ZeroValue%Area
                faceR(face_P,fr_Area_u) = setting%ZeroValue%Area
            endwhere

            where (faceR(face_P,fr_Area_d) >= setting%ZeroValue%Area)
                faceR(face_P,fr_Velocity_d) = faceR(face_P,fr_Flowrate)/faceR(face_P,fr_Area_d)
                faceR(face_P,fr_Velocity_u) = faceR(face_P,fr_Velocity_d)  
            endwhere
        else
            !% ensure face area_u is not smaller than zerovalue
            where (faceR(face_P,fr_Area_d) < zeroR)
                faceR(face_P,fr_Area_d) = zeroR
            endwhere

            where (faceR(face_P,fr_Area_d) >= zeroR)
                faceR(face_P,fr_Velocity_d) = faceR(face_P,fr_Flowrate)/faceR(face_P,fr_Area_d)
                faceR(face_P,fr_Velocity_u) = faceR(face_P,fr_Velocity_d)
            endwhere
        endif

        !%  limit high velocities
        if (setting%Limiter%Velocity%UseLimitMax) then
            where(abs(faceR(face_P,fr_Velocity_d))  > setting%Limiter%Velocity%Maximum)
                faceR(face_P,fr_Velocity_d) = sign(0.99 * setting%Limiter%Velocity%Maximum, &
                    faceR(face_P,fr_Velocity_d))

                faceR(face_P,fr_Velocity_u) = faceR(face_P,fr_Velocity_d)
            endwhere
        endif

        if (setting%Debug%File%boundary_conditions) print *, '*** leave ', this_image(), subroutine_name
    end subroutine face_interpolation_dnBC_byPack

    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine face_interpolation_interior_byPack (facePackCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Interpolates all faces using a pack
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: facePackCol  !% Column in faceP array for needed pack
        integer, intent(in) :: Npack        !% expected number of packed rows in faceP.
        integer :: fGeoSetU(3), fGeoSetD(3), eGeoSet(3)
        integer :: fHeadSetU(1), fHeadSetD(1), eHeadSet(1)
        integer :: fFlowSet(1), eFlowSet(1)
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'face_interpolation_byPack'
        if (setting%Debug%File%face) print *, '*** enter ', this_image(), subroutine_name
        !%-----------------------------------------------------------------------------
        !% Face values are needed for
        !% Area_u, Area_d, Head_u, Head_d, Flowrate,

        !% not sure if we need
        !% Topwidth_u, Topwidth_d, HydDepth_u, HydDepth_d
        !% Velocity_u, Velocity_d

        !% General approach
        !% interpolate to ..._u
        !% identify hydraulic jumps
        !% set .._u and ..d based on jumps

        !% set the matching sets
        !% THESE SHOULD BE DONE IN A GLOBAL -- MAYBE SETTINGS
        !% Note these can be expanded for other terms to be interpolated.
        fGeoSetU = [fr_Area_u, fr_Topwidth_u, fr_HydDepth_u]
        fGeoSetD = [fr_Area_d, fr_Topwidth_d, fr_HydDepth_d]
        eGeoSet  = [er_Area,   er_Topwidth,   er_HydDepth]

        fHeadSetU = [fr_Head_u]
        fHeadSetD = [fr_Head_d]
        eHeadSet = [er_Head]

        fFlowSet = [fr_Flowrate]
        eFlowSet = [er_Flowrate]

        !% two-sided interpolation to using the upstream face set
        call face_interp_interior_set_byPack &
            (fGeoSetU, eGeoSet, er_InterpWeight_dG, er_InterpWeight_uG, facePackCol, Npack)
        call face_interp_interior_set_byPack &
            (fHeadSetU, eHeadSet, er_InterpWeight_dH, er_InterpWeight_uH, facePackCol, Npack)
        call face_interp_interior_set_byPack &
            (fFlowSet, eFlowSet, er_InterpWeight_dQ, er_InterpWeight_uQ, facePackCol, Npack)

        !% copy upstream to downstream storage at a face
        !% (only for Head and Geometry types)
        !% note that these might be reset by hydraulic jump
        call face_copy_upstream_to_downstream_interior_byPack &
            (fGeoSetD, fGeoSetU, facePackCol, Npack)

        call face_copy_upstream_to_downstream_interior_byPack &
            (fHeadSetD, fHeadSetU, facePackCol, Npack)

        !% calculate the velocity in faces and put limiter
        call adjust_face_dynamic_limit (facePackCol, .true.)

        !% reset all the hydraulic jump interior faces
        call jump_compute

        if (setting%Debug%File%face) print *, '*** leave ', this_image(), subroutine_name
    end subroutine face_interpolation_interior_byPack
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine face_interp_across_images ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !%
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%
        print *, "HACK in face_interp_across_images stub -- may be obsolete"
        stop 7387

    end subroutine face_interp_across_images
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine face_interpolation_shared_byPack (facePackCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Interpolates all the shared faces
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: facePackCol  !% Column in faceP array for needed pack
        integer, intent(in) :: Npack        !% expected number of packed rows in faceP.
        integer :: fGeoSetU(3), fGeoSetD(3), eGeoSet(3)
        integer :: fHeadSetU(1), fHeadSetD(1), eHeadSet(1)
        integer :: fFlowSet(1), eFlowSet(1)
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'face_interpolation_shared_byPack'
        if (setting%Debug%File%face) print *, '*** enter ', this_image(), subroutine_name
        !%-----------------------------------------------------------------------------
        !% Face values are needed for
        !% Area_u, Area_d, Head_u, Head_d, Flowrate,

        !% not sure if we need
        !% Topwidth_u, Topwidth_d, HydDepth_u, HydDepth_d
        !% Velocity_u, Velocity_d

        !% General approach
        !% interpolate to ..._u
        !% identify hydraulic jumps
        !% set .._u and ..d based on jumps

        !% set the matching sets
        !% THESE SHOULD BE DONE IN A GLOBAL -- MAYBE SETTINGS
        !% Note these can be expanded for other terms to be interpolated.
        fGeoSetU = [fr_Area_u, fr_Topwidth_u, fr_HydDepth_u]
        fGeoSetD = [fr_Area_d, fr_Topwidth_d, fr_HydDepth_d]
        eGeoSet  = [er_Area,   er_Topwidth,   er_HydDepth]

        fHeadSetU = [fr_Head_u]
        fHeadSetD = [fr_Head_d]
        eHeadSet = [er_Head]

        fFlowSet = [fr_Flowrate]
        eFlowSet = [er_Flowrate]

        !% two-sided interpolation to using the upstream face set
        call face_interp_shared_set_byPack &
            (fGeoSetU, eGeoSet, er_InterpWeight_dG, er_InterpWeight_uG, facePackCol, Npack)
        call face_interp_shared_set_byPack &
            (fHeadSetU, eHeadSet, er_InterpWeight_dH, er_InterpWeight_uH, facePackCol, Npack)
        call face_interp_shared_set_byPack &
            (fFlowSet, eFlowSet, er_InterpWeight_dQ, er_InterpWeight_uQ, facePackCol, Npack)

        !% copy upstream to downstream storage at a face
        !% (only for Head and Geometry types)
        !% note that these might be reset by hydraulic jump
        call face_copy_upstream_to_downstream_shared_byPack &
            (fGeoSetD, fGeoSetU, facePackCol, Npack)

        call face_copy_upstream_to_downstream_shared_byPack &
            (fHeadSetD, fHeadSetU, facePackCol, Npack)

        call adjust_face_dynamic_limit (facePackCol, .false.)
        !% HACK needs jump computation for across shared faces
        ! print *, "HACK missing hydraulic jump that occurs on shared faces 36987"

        if (setting%Debug%File%face) print *, '*** leave ', this_image(), subroutine_name
    end subroutine face_interpolation_shared_byPack
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine face_interp_set_byMask &
        (fset, eset, eWdn, eWup, faceMaskCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Interpolates to a face for a set of variables using a mask
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: fset(:), eset(:), eWdn, eWup, faceMaskCol
        integer, pointer :: eup(:), edn(:)
        integer :: ii
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'face_interp_set_byMask'
        if (setting%Debug%File%face) print *, '*** enter ', this_image(), subroutine_name
        !%-----------------------------------------------------------------------------
        eup => faceI(:,fi_Melem_uL)
        edn => faceI(:,fi_Melem_dL)
        !%-----------------------------------------------------------------------------
        !% cycle through each element in the set.
        do ii=1,size(fset)
            where (faceM(:,faceMaskCol))
                faceR(:,fset(ii)) = &
                    (+elemR(eup(:),eset(ii)) * elemR(edn(:),eWup) &
                     +elemR(edn(:),eset(ii)) * elemR(eup(:),eWdn) &
                    ) / &
                    ( elemR(edn(:),eWup) + elemR(eup(:),eWdn))
            endwhere
        enddo

        print *, 'in face_interp_set_byMask -- may be obsolete'
        stop 87098

        if (setting%Debug%File%face) print *, '*** enter ', this_image(), subroutine_name
    end subroutine face_interp_set_byMask
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine face_interp_interior_set_byPack &
        (fset, eset, eWdn, eWup, facePackCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Interpolates to a face for a set of variables using a mask
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: fset(:), eset(:), eWdn, eWup, facePackCol, Npack
        integer, pointer :: thisP(:), eup(:), edn(:)
        integer :: ii
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'face_interp_interior_set_byPack'
        if (setting%Debug%File%face) print *, '*** enter ', this_image(), subroutine_name
        !%-----------------------------------------------------------------------------
        thisP => faceP(1:Npack,facePackCol)

        eup => faceI(:,fi_Melem_uL)
        edn => faceI(:,fi_Melem_dL)
        !%-----------------------------------------------------------------------------
        !% cycle through each element in the set.
        !% This is designed for fset and eset being vectors, but it
        !%   is not clear that this is needed.
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
        !% elemR(edn(thisp),eWup)) is the upstream weighting of the downstream element

        if (setting%Debug%File%face) print *, '*** enter ', this_image(), subroutine_name
    end subroutine face_interp_interior_set_byPack
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine face_interp_shared_set_byPack &
        (fset, eset, eWdn, eWup, facePackCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Interpolates to a face for a set of variables using a mask
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: fset(:), eset(:), eWdn, eWup, facePackCol, Npack
        integer, pointer :: thisP, eup, edn, connected_image, ghostUp, ghostDn
        logical, pointer :: isGhostUp, isGhostDn
        integer :: ii, jj
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'face_interp_shared_set_byPack'
        if (setting%Debug%File%face) print *, '*** enter ', this_image(), subroutine_name
        !%-----------------------------------------------------------------------------
        !% cycle through all the shared faces
        do ii = 1,Npack
            thisP           => facePS(ii,facePackCol)
            connected_image => faceI(thisP,fi_Connected_image)
            eup             => faceI(thisP,fi_Melem_uL)
            edn             => faceI(thisP,fi_Melem_dL)
            ghostUp         => faceI(thisP,fi_GhostElem_uL)
            ghostDn         => faceI(thisP,fi_GhostElem_dL)
            isGhostUp       => faceYN(thisP,fYN_isUpGhost)
            isGhostDn       => faceYN(thisP,fYN_isDnGhost)
            !%-----------------------------------------------------------------------------
            !% cycle through each element in the set.
            !% This is designed for fset and eset being vectors, but it
            !%   is not clear that this is needed.
            do jj=1,size(fset)

                !% condition for upstream element of the shared face is ghost and in a different image
                if (isGhostUp) then

                    faceR(thisP,fset(jj)) = &
                        (+elemR(ghostUp,eset(jj))[connected_image] * elemR(edn,eWup) &
                         +elemR(edn,eset(jj)) * elemR(ghostUp,eWdn)[connected_image] &
                        ) / &
                        ( elemR(edn,eWup) + elemR(ghostUp,eWdn)[connected_image] )

                !% condition for downstream element of the shared face is ghost and in a different image
                elseif (isGhostDn) then

                    faceR(thisP,fset(jj)) = &
                        (+elemR(eup,eset(jj)) * elemR(ghostDn,eWup)[connected_image] &
                         +elemR(ghostDn,eset(jj))[connected_image] * elemR(eup,eWdn) &
                        ) / &
                        ( elemR(ghostDn,eWup)[connected_image] + elemR(eup,eWdn) )
                endif
            end do
        end do

        !% NOTES
        !% elemR(eup,eset(jj)) is the element value upstream of the face
        !% elemR(edn,eset(jj) is the element value downstream of the face.
        !% elemR(eup,eWdn) is the downstream weighting of the upstream element
        !% elemR(edn,eWup)) is the upstream weighting of the downstream element

        !% elemR(ghostUp,eset(jj))[connected_image] is the elem value from the upstream image of the face
        !% elemR(ghostDn,eset(jj))[connected_image] is the elem value from the downstream image of the face
        !% elemR(ghostUp,eWdn)[connected_image] is the downstream weighting of the upstream image element
        !% elemR(ghostDn,eWup))[connected_image] is the upstream weighting of the downstream image element

        if (setting%Debug%File%face) print *, '*** enter ', this_image(), subroutine_name
    end subroutine face_interp_shared_set_byPack
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine face_copy_upstream_to_downstream_interior_byPack &
        (downstreamSet, upstreamSet, facePackCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Copies the interpolated value on the upstream side to the downstream side
        !% These values are later adjusted for hydraulic jumps
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: facePackCol, Npack, downstreamSet(:), upstreamSet(:)
        integer, pointer :: thisP(:)
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'face_copy_upstream_to_downstream_interior_byPack'
        if (setting%Debug%File%face) print *, '*** enter ', this_image(), subroutine_name
        !%-----------------------------------------------------------------------------

        thisP => faceP(1:Npack,facePackCol)

        faceR(thisP,downstreamSet) = faceR(thisP,upstreamSet)

        if (setting%Debug%File%face) print *, '*** leave ', this_image(), subroutine_name
    end subroutine face_copy_upstream_to_downstream_interior_byPack
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine face_copy_upstream_to_downstream_shared_byPack &
        (downstreamSet, upstreamSet, facePackCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Copies the interpolated value on the upstrea side to the downstream side
        !% These values are later adjusted for hydraulic jumps
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: facePackCol, Npack, downstreamSet(:), upstreamSet(:)
        integer, pointer :: thisP(:)
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'face_copy_upstream_to_downstream_byPack'
        if (setting%Debug%File%face) print *, '*** enter ', this_image(), subroutine_name
        !%-----------------------------------------------------------------------------

        thisP => facePS(1:Npack,facePackCol)

        faceR(thisP,downstreamSet) = faceR(thisP,upstreamSet)

        if (setting%Debug%File%face) print *, '*** leave ', this_image(), subroutine_name
    end subroutine face_copy_upstream_to_downstream_shared_byPack
    !%
    !%==========================================================================
    !% END OF MODULE
    !%==========================================================================
end module face
