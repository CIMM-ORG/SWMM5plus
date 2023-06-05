module jump
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Computes jump conditions across a face
    !%
    !% Methods:
    !% Extrapolates heads from adjacent elements at supercritical
    !% to subcritical transition
    !%==========================================================================
    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use utility_crash

    implicit none

    private

    public :: jump_compute

    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================  
!%  
    subroutine jump_compute
        !%------------------------------------------------------------------
        !% Description:
        !% Computes a hydraulic jump at a face 
        !%------------------------------------------------------------------
        !% Declarations
            integer, pointer :: facePackCol, Npack
            character(64) :: subroutine_name = 'jump_compute'
        !%------------------------------------------------------------------
        !% Preliminareis
            if (setting%Debug%File%jump) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]" 
        !%------------------------------------------------------------------
        !%  
        !% --- identify hydraulic jump (create pack facemap in global)
        call jump_face_identify
        
        !% --- enforce hydraulic jump on downstream jumps
        facePackCol => col_faceP(fp_JumpDn_IorS)
        Npack       => npack_faceP(fp_JumpDn_IorS)
        if (Npack > 0) then
            call jump_enforce (facePackCol, Npack, jump_from_downstream)
        end if
        
        !% --- enforce hydraulic jump on upstream jumps
        facePackCol => col_faceP(fp_JumpUp_IorS)
        Npack       => npack_faceP(fp_JumpUp_IorS)
        if (Npack > 0) then
            call jump_enforce (facePackCol, Npack, jump_from_upstream)
        end if

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%jump)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine jump_compute   
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%    
    subroutine jump_face_identify
        !%------------------------------------------------------------------
        !% Description:
        !% This identifies and packs the hydraulic jump locations
        !% in the faceP array. Because the Froude number is computed
        !% for all elements (including surcharged), we need to limit
        !% this selection to the non-surcharged elements. This means
        !% that when a supercritical open channel flow meets a 
        !% surcharged pipe the face elevation will use a standard 
        !% interpolation (which will likely make the face surcharged).
        !%
        !% This subroutine should be modified to be two calls to a single function 
        !% that packs for either the upstream or downstream elements. However, do not
        !% rewrite this until we are sure that we actually need to discriminate between
        !% upstream and downstream jumps. We might be able to use a single function that
        !% simply identifies any jump condition.
        !%
        !% Note that limiting the jump conditions to non-surcharged elements on both 
        !% sides is a slight change from the prototype code where the jump was only 
        !% limited when the surcharged was on the downstream side.
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: image
            integer, pointer :: faceIdx(:), eup(:), edn(:), thisP(:), jumptype(:)
            integer, pointer :: Npack_jumpUp, Npack_JumpDn, Nfaces
            logical, pointer :: isSurcharged(:), isInterior(:), isZeroDepth(:), isSmallDepth(:)
            real(8), pointer :: feps, Fr(:), Head(:)
        !%------------------------------------------------------------------
        !% Aliases
            image  = this_image() !% not a pointer!
            Nfaces => N_face(image)
            Fr           => elemR(:,er_FroudeNumber)
            Head         => elemR(:,er_Head)
            isSurcharged => elemYN(:,eYN_isSurcharged)
            isZeroDepth  => elemYN(:,eYN_isZeroDepth)
            isSmallDepth => elemYN(:,eYN_isSmallDepth)
        !%------------------------------------------------------------------
            isInterior   => faceYN(1:Nfaces,fYN_isInteriorFace)
            eup          => faceI(1:Nfaces,fi_Melem_uL)
            faceIdx      => faceI(1:Nfaces,fi_Lidx) 
            edn          => faceI(1:Nfaces,fi_Melem_dL)
            jumptype     => faceI(1:Nfaces,fi_jump_type)
            feps         => setting%EPS%FroudeJump
        !%------------------------------------------------------------------

        !% --- zero out old jump
        faceI(1:Nfaces,fi_jump_type) = jump_none

        !% --- count the number of faces with flow in nominal downstream and a jump 
        !%     from upstream (supercritical) to downstream (subcritical)
        !%     and open-channel flow on either side
        npack_faceP(fp_JumpUp_IorS) = count( &
            isInterior &
            .and. &
            (Head(eup) < Head(edn)) &
            .and. &
            (Fr(eup) > oneR + feps)  &  ! supercritical downstream flow in upstream
            .and. &
            (Fr(edn) < oneR - feps)   &  ! subcritical in downstream
            .and. &
            (Fr(edn) >= zeroR)        &  ! not a reverse flow downstream
            .and. &
            (.not. isSurcharged(eup)) &
            .and. &
            (.not. isSurcharged(edn))  & 
            .and. &
            (.not. isSmallDepth(eup)) &
            .and. & 
            (.not. isSmallDepth(edn)) &
            .and. &
            (.not. isZeroDepth(eup)) &
            .and. &
            (.not. isZeroDepth(edn)) )

        Npack_JumpUp => npack_faceP(fp_JumpUp_IorS) 

        !%--- pack the indexes 
        if (Npack_JumpUp > 0) then
            faceP(1:Npack_JumpUp, fp_JumpUp_IorS) = pack(faceIdx, &
                isInterior &
                .and. &
                (Head(eup) < Head(edn)) &
                .and. &
                (Fr(eup) > oneR + feps)  &   !% supercritical downstream flow in upstream
                .and. &
                (Fr(edn) < oneR - feps)   &  !% subcritical in downstream
                .and. &
                (Fr(edn) >= zeroR)        &  !% not a reverse flow downstream
                .and. &
                (.not. isSurcharged(eup)) &
                .and. &
                (.not. isSurcharged(edn)) & 
                .and. &
                (.not. isSmallDepth(eup)) &
                .and. & 
                (.not. isSmallDepth(edn)) &
                .and. &
                (.not. isZeroDepth(eup)) &
                .and. &
                (.not. isZeroDepth(edn)) )

            !% --- pointer to the packed indexes
            thisP => faceP(1:Npack_JumpUp,fp_JumpUp_IorS)

            !% --- designate these as an upstream jump
            jumptype(thisP) = jump_from_upstream
        end if



        !% --- count the number of faces with reverse flow and a jump from 
        !%     downstream (supercritical) to upstream (subcritical)
        !%     and open channel flow on either side
        npack_faceP(fp_JumpDn_IorS)  = count( &
            isInterior &
            .and. &
            (Head(eup) > Head(edn)) &
            .and. &
            (Fr(eup) > -oneR + feps) &   !% subcritical reverse flow in upstream
            .and. &
            (Fr(eup) <= zeroR) &         !% not a downstream flow in upstream
            .and. &
            (Fr(edn) < -oneR - feps) &   !% supercritical reverse flow in downstream
            .and. &
            (.not. isSurcharged(eup)) &
            .and. &
            (.not. isSurcharged(edn))  & 
            .and. &
            (.not. isSmallDepth(eup)) &
            .and. & 
            (.not. isSmallDepth(edn)) &
            .and. &
            (.not. isZeroDepth(eup)) &
            .and. &
            (.not. isZeroDepth(edn)) )

        !% --- assign the above count to the npack storage for later use
        Npack_JumpDn => npack_faceP(fp_JumpDn_IorS)

        !% --- pack the indexes 
        if (Npack_JumpDn > 0) then
            faceP(1:Npack_JumpDn, fp_JumpDn_IorS) = pack(faceIdx, &
                isInterior &
                .and. &
                (Head(eup) > Head(edn)) &
                .and. & 
                (Fr(eup) > -oneR + feps) &  !% subcritical reverse flow in upstream
                .and. &
                (Fr(eup) <= zeroR) &        !% not a downstream flow in upstream
                .and. &
                (Fr(edn) < -oneR - feps) &  !% supercritical reverse flow in downstream
                .and. &
                (.not. isSurcharged(eup)) &
                .and. &
                (.not. isSurcharged(edn))  & 
                .and. &
                (.not. isSmallDepth(eup)) &
                .and. & 
                (.not. isSmallDepth(edn)) &
                .and. &
                (.not. isZeroDepth(eup)) &
                .and. &
                (.not. isZeroDepth(edn)) )
            
            !% --- pointer to thee packed indexes
            thisP => faceP(1:Npack_JumpDn,fp_JumpDn_IorS)

            !% --- designate these as an upstream jump
            jumptype(thisP) = jump_from_downstream

        end if

    end subroutine jump_face_identify
!%
!%========================================================================== 
!%==========================================================================    
!%  
    subroutine jump_enforce (facePackCol, Npack, jump_from)
        !%------------------------------------------------------------------
        !% Description:
        !% Enforces the hydraulic jump condition at a face
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: facePackCol, Npack, jump_from
            integer, pointer    :: thisP(:), eUp(:), eDn(:)
            integer             :: fsetUp(2), fsetDn(2), eset(2)
        !%------------------------------------------------------------------ 
        !% Aliases
            thisP => faceP(1:Npack,facePackCol)     
            eUp => faceI(:,fi_Melem_uL)
            eDn => faceI(:,fi_Melem_dL)
        !%------------------------------------------------------------------ 

        !% --- columns for face data that requires jump conditions
        fsetUp = [fr_Area_u, fr_Head_u]
        fsetDn = [fr_Area_d, fr_Head_d]
        !% --- columns for elem data stored on either side of jump
        eset = [er_Area,  er_Head]

        select case (jump_from)
            case (jump_from_upstream)
                !% --- enforce jump from upstream faces
                !%     so that downstream face is the downstream element value
                faceR(thisP,fsetDn) = elemR(eDn(thisP),eset)
                faceR(thisP,fr_Velocity_u) = faceR(thisP,fr_flowrate) / faceR(thisP,fr_Area_u)
                faceR(thisP,fr_Velocity_d) = faceR(thisP,fr_flowrate) / faceR(thisP,fr_Area_d)
                
            case (jump_from_downstream)
                !% --- enforce jump from downstream faces
                !%     so that upstream face is the upstream element value
                faceR(thisP,fsetUp) = elemR(eUp(thisP),eset)
                faceR(thisP,fr_Velocity_u) = faceR(thisP,fr_flowrate) / faceR(thisP,fr_Area_u)
                faceR(thisP,fr_Velocity_d) = faceR(thisP,fr_flowrate) / faceR(thisP,fr_Area_d)

            case default
                print *, 'CODE ERROR geometry type unknown for # ', jump_from
                print *, 'which has key ',trim(reverseKey(jump_from))
                call util_crashpoint(33782)
        end select
        
    end subroutine jump_enforce
!%
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module jump