module face

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
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
        if (setting%Debug%File%face) print *, '*** enter ', subroutine_name
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

        if (setting%Debug%File%face)  print *, '*** leave ', subroutine_name
    end subroutine face_interpolation
    !%
    !%==========================================================================
    !% PRIVATE
    !%==========================================================================
    !%
    subroutine face_interpolation_byMask (faceMaskCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Interpolates all faces using a mask -- assumes single processor
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: faceMaskCol !% Column in face array containing mask for all valid faces
        integer :: fGeoSetU(3), fGeoSetD(3), eGeoSet(3)
        integer :: fHeadSetU(1), fHeadSetD(1), eHeadSet(1)
        integer :: fFlowSet(1), eFlowSet(1)
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'face_interpolation_byMask'
        if (setting%Debug%File%face) print *, '*** enter ', subroutine_name
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
        call face_interp_set_byMask &
            (fGeoSetU, eGeoSet, er_InterpWeight_dG, er_InterpWeight_uG, faceMaskCol)
        call face_interp_set_byMask &
            (fHeadSetU, eHeadSet, er_InterpWeight_dH, er_InterpWeight_uH, faceMaskCol)
        call face_interp_set_byMask &
            (fFlowSet, eFlowSet, er_InterpWeight_dQ, er_InterpWeight_uQ, faceMaskCol)

        !% copy upstream to downstream storage at a face
        !% (only for Head and Geometry types as flow has a single value)
        !% note that these might be reset by hydraulic jump
        call face_copy_upstream_to_downstream_byMask (fGeoSetD,  fGeoSetU,  faceMaskCol)
        call face_copy_upstream_to_downstream_byMask (fHeadSetD, fHeadSetU, faceMaskCol)

        !% reset all the hydraulic jump faces
        call jump_compute

        if (setting%Debug%File%face) print *, '*** leave ', subroutine_name
    end subroutine face_interpolation_byMask
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
        if (setting%Debug%File%face) print *, '*** enter ', subroutine_name
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

        !% reset all the hydraulic jump faces
        call jump_compute

        if (setting%Debug%File%face) print *, '*** leave ', subroutine_name
    end subroutine face_interpolation_interior_byPack
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine face_interp_across_images ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%
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
        if (setting%Debug%File%face) print *, '*** enter ', subroutine_name
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

        !% HACK needs junp computation for across shared faces

        if (setting%Debug%File%face) print *, '*** leave ', subroutine_name
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
        if (setting%Debug%File%face) print *, '*** enter ', subroutine_name
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

        if (setting%Debug%File%face) print *, '*** enter ', subroutine_name
    end subroutine face_interp_set_byMask
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine face_copy_upstream_to_downstream_byMask &
        (downstreamSet, upstreamSet, faceMaskCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Copies the interpolated value on the upstrea side to the downstream side
        !% These values are later adjusted for hydraulic jumps
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: faceMaskCol, downstreamSet(:), upstreamSet(:)
        integer :: ii
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'face_copy_upstream_to_downstream_byMask'
        if (setting%Debug%File%face) print *, '*** enter ', subroutine_name
        !%-----------------------------------------------------------------------------

        do ii=1,size(downstreamset)
            where (faceM(:,faceMaskCol))
                faceR(:,downstreamSet(ii)) = faceR(:,upstreamSet(ii))
            endwhere
        enddo

        if (setting%Debug%File%face) print *, '*** leave ', subroutine_name
    end subroutine face_copy_upstream_to_downstream_byMask
    !%
    !%==========================================================================
    !%==========================================================================
    !%
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
        if (setting%Debug%File%face) print *, '*** enter ', subroutine_name
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

        if (setting%Debug%File%face) print *, '*** enter ', subroutine_name
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
        integer :: ii, jj
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'face_interp_shared_set_byPack'
        if (setting%Debug%File%face) print *, '*** enter ', subroutine_name
        !%-----------------------------------------------------------------------------
        !% cycle through all the shared faces
        do ii = 1,Npack
            thisP           => facePS(ii,facePackCol)
            connected_image => faceI(thisP,fi_Connected_image)
            eup             => faceI(thisP,fi_Melem_uL)
            edn             => faceI(thisP,fi_Melem_dL)
            ghostUp         => faceI(thisP,fi_GhostElem_uL)
            ghostDn         => faceI(thisP,fi_GhostElem_dL)
            !%-----------------------------------------------------------------------------
            !% cycle through each element in the set.
            !% This is designed for fset and eset being vectors, but it
            !%   is not clear that this is needed.
            do jj=1,size(fset)

                !% condition for upstream element of the shared face is in a different image
                if (eup .eq. nullValueI) then

                    faceR(thisP,fset(jj)) = &
                        (+elemR(ghostUp,eset(jj))[connected_image] * elemR(edn,eWup) &
                         +elemR(edn,eset(jj)) * elemR(ghostUp,eWdn)[connected_image] &
                        ) / &
                        ( elemR(edn,eWup) + elemR(ghostUp,eWdn)[connected_image] )
                print*, this_image(), 'image'
                print*, connected_image, 'connected_image'
                print*, thisP, 'face idx'  
                print*, faceR(thisP,fset(jj)), 'faceR(thisP,fset(jj))'  
                !% condition for downstream element of the shared face is in a different image
                elseif (edn .eq. nullValueI) then

                    faceR(thisP,fset(jj)) = &
                        (+elemR(eup,eset(jj)) * elemR(ghostDn,eWup)[connected_image] &
                         +elemR(ghostDn,eset(jj))[connected_image] * elemR(eup,eWdn) &
                        ) / &
                        ( elemR(ghostDn,eWup)[connected_image] + elemR(eup,eWdn) )
                print*, this_image(), 'image'
                print*, connected_image, 'connected_image'
                print*, thisP, 'face idx'  
                print*, faceR(thisP,fset(jj)), 'faceR(thisP,fset(jj))' 
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

        if (setting%Debug%File%face) print *, '*** enter ', subroutine_name
    end subroutine face_interp_shared_set_byPack
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine face_copy_upstream_to_downstream_interior_byPack &
        (downstreamSet, upstreamSet, facePackCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Copies the interpolated value on the upstrea side to the downstream side
        !% These values are later adjusted for hydraulic jumps
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: facePackCol, Npack, downstreamSet(:), upstreamSet(:)
        integer, pointer :: thisP(:)
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'face_copy_upstream_to_downstream_interior_byPack'
        if (setting%Debug%File%face) print *, '*** enter ', subroutine_name
        !%-----------------------------------------------------------------------------

        thisP => faceP(1:Npack,facePackCol)

        faceR(thisP,downstreamSet) = faceR(thisP,upstreamSet)

        if (setting%Debug%File%face) print *, '*** leave ', subroutine_name
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
        if (setting%Debug%File%face) print *, '*** enter ', subroutine_name
        !%-----------------------------------------------------------------------------

        thisP => facePS(1:Npack,facePackCol)

        faceR(thisP,downstreamSet) = faceR(thisP,upstreamSet)

        if (setting%Debug%File%face) print *, '*** leave ', subroutine_name
    end subroutine face_copy_upstream_to_downstream_shared_byPack
    !%
    !%==========================================================================
    !% END OF MODULE
    !%==========================================================================
end module face
