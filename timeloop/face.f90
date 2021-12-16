module face

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use adjust
    use jump
    use utility_profiler


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

    public :: face_interpolation
    public :: face_interpolate_bc

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine face_interpolation (facecol,whichTM)
        !%------------------------------------------------------------------
        !% Description:
        !% Interpolates faces from elements
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in)  :: faceCol, whichTM
            integer, pointer :: Npack
            logical :: isBConly
            
            character(64) :: subroutine_name = 'face_interpolation'
        !%-------------------------------------------------------------------
        !% Preliminaries    
            if (icrash) return
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            
            if (setting%Profile%useYN) call util_profiler_start (pfc_face_interpolation)
        !%--------------------------------------------------------------------
        isBConly = .false.

        !print *, 'in ',trim(subroutine_name),' -------------------'    
        ! print *, faceR(:,fr_Flowrate)
        ! print *, ' '
        !print *, size(faceR,DIM=1)
        !print *, size(elemR,DIM=1)

        Npack => npack_faceP(faceCol)
        if (Npack > 0) then
            !% face reconstruction of all the interior faces
            call face_interpolation_interior (faceCol, Npack)
        end if

        Npack => npack_facePS(faceCol)
        if (Npack > 0) then
            sync all
            !% face reconstruction of all the shared faces
            call face_interpolation_shared (faceCol, Npack)
            sync all
        end if

        ! print *, ' '
        ! print *, 'before BC'
        ! print *, faceR(:,fr_Flowrate)
        ! print *, ' '

        call face_interpolate_bc (isBConly)

        !print *, 'before Qinterp'
        !print *, faceR(:,fr_Flowrate)
        !print *, ' '

        !% adjust face flowrate for QinterpWithLocalHeadGradient method
        if (setting%Solver%QinterpWithLocalHeadGradient) then 
        
            ! print *, 'in ',trim(subroutine_name)
            ! print *,' aaaa'
            ! !print *, faceR(8,fr_Head_u)
            ! !stop 39780
            ! ! print *, faceR(1:8,fr_Head_u)
            ! ! print *, ' '
            ! ! print *, faceR(1:8,fr_Head_d)
            ! ! print *, ' '
            ! print *, elemR(7,er_HeadAvg), elemR(7,er_Head)

            !% compute the average head for the elements updated before call to face_interpolation
            call face_head_average_on_element (whichTM)

            ! print *, ' bbbb'
            ! print *, elemR(7,er_HeadAvg), elemR(7,er_Head)
            ! stop 39705

            Npack => npack_faceP(faceCol)
            if (Npack > 0) then 
                !% Add the Q gradient term to the flowrate for interior faces
                call face_Q_HeadGradient_interior (faceCol, Npack)      
            end if

            !print *, ' cccc'
            !print *, faceR(:,fr_Flowrate)

            Npack => npack_facePS(faceCol)
            if (Npack > 0) then
                sync all
                !% add the Q gradient term to the flowrate for shared faces
                call face_Q_HeadGradient_shared (faceCol, Npack)
                sync all
            end if

            !% need another call to face_interpolate so that the Q_HeadGradient
            !% does not change the upper boundarin inflow conditon
            isBConly = .true.
            call face_interpolate_bc (isBConly)

            !print *, ' dddd'
            !print *, faceR(:,fr_Flowrate)
            !print *, ' '
            !stop 398705
        end if
        
        !print *, 'at end'
        !print *, faceR(:,fr_Flowrate)
        !stop 877386

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
            if (icrash) return
            if (setting%Debug%File%face)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------

        !print *, 'in 763764 ',trim(subroutine_name)
        !write(*,'(8f9.2)') faceR(1:N_elem(1)+1,fr_Flowrate)

        if (N_nBCup > 0) call face_interpolation_upBC(isBConly)

        !print *, 'in 9353 ',trim(subroutine_name)
        !write(*,'(8f9.2)') faceR(1:N_elem(1)+1,fr_Flowrate)

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
            integer :: fGeoSetU(3), fGeoSetD(3), eGeoSet(3)
            integer :: ii
            integer, pointer :: face_P(:), edn(:), idx_P(:), fdn(:)
            character(64) :: subroutine_name = 'face_interpolation_upBC'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (icrash) return
            if (setting%Debug%File%boundary_conditions)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% Aliases
        !% For the head/geometry at the upstream faces, we directly take the dnwnstream element
        !% So there is no eup for upstream BCs
            edn => faceI(:,fi_Melem_dL)
            fdn => elemI(:,ei_Mface_dL)
        
            face_P => faceP(1:npack_faceP(fp_BCup),fp_BCup)
            idx_P  => BC%P%BCup(:)
        !%-------------------------------------------------------------------
        !% enforce BC    
        faceR(face_P, fr_Flowrate) = BC%flowRI(idx_P)

        !print *, 'in 22975 ',trim(subroutine_name)
        !print *, idx_P
        !print *, face_P
        !print *, BC%flowRI(idx_P(1))
        !print *, faceR(face_P(1), fr_Flowrate)
        !print *, BC%flowRI(:)

        !% update geometry data (don't do on a BC-only call)
        if (.not. isBConly) then
            !% Define sets of points for the interpolation, we are going from
            !% the elements to the faces.       
            fGeoSetU = [fr_Area_u, fr_Topwidth_u, fr_HydDepth_u]
            fGeoSetD = [fr_Area_d, fr_Topwidth_d, fr_HydDepth_d]
            eGeoSet  = [er_Area,   er_Topwidth,   er_HydDepth]

            !% Copying other data to the BC faces
            do ii = 1,size(fGeoSetU)
                faceR(face_P,fGeoSetD(ii)) = elemR(edn(face_P),eGeoSet(ii)) ! Copying the Geo
                faceR(face_P,fGeoSetU(ii)) = faceR(face_P,fGeoSetD(ii))
            end do

            !% gradient extrapolation for head from PipeAC20
            faceR(face_P,fr_Head_d) = elemR(edn(face_P),er_Head) + elemR(edn(face_P),er_Head) - &
                        faceR(fdn(edn(face_P)),fr_Head_d) !% Bugfix SAZ 09212021

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
            else
                !% ensure face area_u is not smaller than zerovalue
                where (faceR(face_P,fr_Area_d) <= zeroR)
                    faceR(face_P,fr_Area_d)     = zeroR
                    faceR(face_P,fr_Area_u)     = zeroR
                    faceR(face_P,fr_Velocity_d) = zeroR
                    faceR(face_P,fr_Velocity_u) = zeroR
                endwhere
            end if  

        end if

        !% set velocity based on flowrate an limit
        if (setting%ZeroValue%UseZeroValues) then
            where (faceR(face_P,fr_Area_d) > setting%ZeroValue%Area)
                faceR(face_P,fr_Velocity_d) = faceR(face_P,fr_Flowrate)/faceR(face_P,fr_Area_d)
                faceR(face_P,fr_Velocity_u) = faceR(face_P,fr_Velocity_d)  
            endwhere
        else     
            where (faceR(face_P,fr_Area_d) > zeroR)
                faceR(face_P,fr_Velocity_d) = faceR(face_P,fr_Flowrate)/faceR(face_P,fr_Area_d)
                faceR(face_P,fr_Velocity_u) = faceR(face_P,fr_Velocity_d)
            endwhere
        end if

        !%  limit high velocities
        if (setting%Limiter%Velocity%UseLimitMaxYN) then
            where(abs(faceR(face_P,fr_Velocity_d))  > setting%Limiter%Velocity%Maximum)
                faceR(face_P,fr_Velocity_d) = sign(0.99 * setting%Limiter%Velocity%Maximum, &
                    faceR(face_P,fr_Velocity_d))

                faceR(face_P,fr_Velocity_u) = faceR(face_P,fr_Velocity_d)
            endwhere
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
    !     if (icrash) return
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
    !         elemR(elem_P,eFlowSet(ii)) = BC%flowRI(idx_P)
    !     end do
    !     !% For lateral flow, just update the flow at the element >> elemR(flow) + BC_lateral_flow

    !     if (setting%Debug%File%boundary_conditions) &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine face_interpolation_latBC_byPack
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine face_interpolation_dnBC(isBConly)
        !%------------------------------------------------------------------
        !% Description:
        !% Interpolates data to all downstream boundary faces
        !% When called with isBConly == true, only does the BC update
        !%-------------------------------------------------------------------
        !% Declarations
            logical, intent(in) :: isBConly
            integer :: fGeoSetU(3), fGeoSetD(3), eGeoSet(3)
            integer :: ii
            integer, pointer :: face_P(:), eup(:), idx_P(:)
            real :: DownStreamBcHead
            character(64) :: subroutine_name = 'face_interpolation_dnBC'
        !%--------------------------------------------------------------------
        !% Preliminaries
            if (icrash) return
            if (setting%Debug%File%boundary_conditions)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%--------------------------------------------------------------------
        !% Aliases
            eup    => faceI(:,fi_Melem_uL)
            face_P => faceP(1:npack_faceP(fp_BCdn),fp_BCdn)
            idx_P  => BC%P%BCdn
        !%--------------------------------------------------------------------   
        !% For fixed, tidal, and timeseries BC
        !% The BC is imagined as enforced on a ghost cell outside the boundary
        !% so the face value is given by linear interpolation using ghost and interior cells 
        faceR(face_P, fr_Head_u) = 0.5 * (elemR(eup(face_P), er_Head) + BC%headRI(idx_P)) !% downstream head update
        faceR(face_P, fr_Head_d) = faceR(face_P, fr_Head_u)

        if (.not. isBConly) then
            
            fGeoSetU = [fr_Area_u, fr_Topwidth_u, fr_HydDepth_u]
            fGeoSetD = [fr_Area_d, fr_Topwidth_d, fr_HydDepth_d]
            eGeoSet  = [er_Area,   er_Topwidth,   er_HydDepth]

            faceR(face_P, fr_Flowrate) = elemR(eup(face_P), er_Flowrate) !% Copying the flow from the upstream element

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
            end if

            !%  limit high velocities
            if (setting%Limiter%Velocity%UseLimitMaxYN) then
                where(abs(faceR(face_P,fr_Velocity_d))  > setting%Limiter%Velocity%Maximum)
                    faceR(face_P,fr_Velocity_d) = sign(0.99 * setting%Limiter%Velocity%Maximum, &
                        faceR(face_P,fr_Velocity_d))

                    faceR(face_P,fr_Velocity_u) = faceR(face_P,fr_Velocity_d)
                endwhere
            end if
        else
            !% continue
        end if
      
        !% endif
        !%--------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%boundary_conditions) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine face_interpolation_dnBC
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine face_interpolation_dnBC(isBConly)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Interpolates all boundary faces using a pack arrays -- base on bi_category
    !     !%-----------------------------------------------------------------------------
    !     integer :: fGeoSetU(3), fGeoSetD(3), eGeoSet(3)
    !     integer :: fFlowSet(1), eFlowSet(1)
    !     integer :: fHeadSetU(1), fHeadSetD(1), eHeadSet(1)
    !     character(64) :: subroutine_name = 'face_interpolation_dnBC_byPack'
    !     integer :: ii
    !     integer, pointer :: face_P(:), eup(:), idx_P(:), bcType
    !     real :: DownStreamBcHead

    !     !%-----------------------------------------------------------------------------
    !     if (icrash) return
    !     if (setting%Debug%File%boundary_conditions)  &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"


    !     eup => faceI(:,fi_Melem_uL)

    !     face_P => faceP(1:npack_faceP(fp_BCdn),fp_BCdn)
    !     idx_P  => BC%P%BCdn

    !     fGeoSetU = [fr_Area_u, fr_Topwidth_u, fr_HydDepth_u]
    !     fGeoSetD = [fr_Area_d, fr_Topwidth_d, fr_HydDepth_d]
    !     eGeoSet  = [er_Area,   er_Topwidth,   er_HydDepth]

    !     fHeadSetU = [fr_Head_u]
    !     fHeadSetD = [fr_Head_d]
    !     eHeadSet = [er_Head]

    !     fFlowSet = [fr_Flowrate]
    !     eFlowSet = [er_Flowrate]


    !     do ii=1,size(fHeadSetD)
    !         !%  linear interpolation using ghost and interior cells
    !         faceR(face_P, fHeadSetU(ii)) = 0.5 * (elemR(eup(face_P), er_Head) + BC%headRI(idx_P)) !% downstream head update
    !         faceR(face_P, fHeadSetD(ii)) = faceR(face_P, fHeadSetU(ii))
    !     end do

    !     do ii=1,size(fFlowSet)
    !         faceR(face_P, fFlowSet(ii)) = elemR(eup(face_P), eFlowSet(ii)) !% Copying the flow from the upstream element
    !     end do

    !     do ii=1,size(fGeoSetD)
    !         faceR(face_P, fGeoSetD(ii)) = elemR(eup(face_P), eGeoSet(ii)) !% Copying other geo factors from the upstream element
    !         faceR(face_P, fGeoSetU(ii)) = faceR(face_P, fGeoSetD(ii))
    !     end do

    !     !% HACK: This is needed to be revisited later
    !     if (setting%ZeroValue%UseZeroValues) then
    !         !% ensure face area_u is not smaller than zerovalue
    !         where (faceR(face_P,fr_Area_d) < setting%ZeroValue%Area)
    !             faceR(face_P,fr_Area_d) = setting%ZeroValue%Area
    !             faceR(face_P,fr_Area_u) = setting%ZeroValue%Area
    !         endwhere

    !         !% HACK: This is needed to be revisited later
    !         if (setting%ZeroValue%UseZeroValues) then
    !             !% ensure face area_u is not smaller than zerovalue
    !             where (faceR(face_P,fr_Area_d) < setting%ZeroValue%Area)
    !                 faceR(face_P,fr_Area_d) = setting%ZeroValue%Area
    !                 faceR(face_P,fr_Area_u) = setting%ZeroValue%Area
    !             endwhere

    !             where (faceR(face_P,fr_Area_d) >= setting%ZeroValue%Area)
    !                 faceR(face_P,fr_Velocity_d) = faceR(face_P,fr_Flowrate)/faceR(face_P,fr_Area_d)
    !                 faceR(face_P,fr_Velocity_u) = faceR(face_P,fr_Velocity_d)  
    !             endwhere
    !         else
    !             !% ensure face area_u is not smaller than zerovalue
    !             where (faceR(face_P,fr_Area_d) < zeroR)
    !                 faceR(face_P,fr_Area_d) = zeroR
    !             endwhere

    !             where (faceR(face_P,fr_Area_d) >= zeroR)
    !                 faceR(face_P,fr_Velocity_d) = faceR(face_P,fr_Flowrate)/faceR(face_P,fr_Area_d)
    !                 faceR(face_P,fr_Velocity_u) = faceR(face_P,fr_Velocity_d)
    !             endwhere
    !         end if

    !         !%  limit high velocities
    !         if (setting%Limiter%Velocity%UseLimitMaxYN) then
    !             where(abs(faceR(face_P,fr_Velocity_d))  > setting%Limiter%Velocity%Maximum)
    !                 faceR(face_P,fr_Velocity_d) = sign(0.99 * setting%Limiter%Velocity%Maximum, &
    !                     faceR(face_P,fr_Velocity_d))

    !                 faceR(face_P,fr_Velocity_u) = faceR(face_P,fr_Velocity_d)
    !             endwhere
    !         end if
    !     else
    !         !% continue
    !     end if
      
    !     !% endif
    !     !%--------------------------------------------------------------------
    !     !% Closing
    !         if (setting%Debug%File%boundary_conditions) &
    !             write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine face_interpolation_dnBC

!%
!%==========================================================================
!%==========================================================================
!%
    subroutine face_interpolation_interior (facePackCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% Interpolates all faces using a pack
        !%------------------------------------------------------------------
            integer, intent(in) :: facePackCol  !% Column in faceP array for needed pack
            integer, intent(in) :: Npack        !% expected number of packed rows in faceP.
            integer :: fGeoSetU(3), fGeoSetD(3), eGeoSet(3)
            integer :: fHeadSetU(1), fHeadSetD(1), eHeadSet(1)
            integer :: fFlowSet(1), eFlowSet(1)
            character(64) :: subroutine_name = 'face_interpolation_interior'
        !%------------------------------------------------------------------
        !% Preliminaries    
            if (icrash) return
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Face values are needed for
        !% Area_u, Area_d, Head_u, Head_d, Flowrate,

        !% HACK: not sure if we need
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
        call face_interp_interior_set &
            (fGeoSetU, eGeoSet, er_InterpWeight_dG, er_InterpWeight_uG, facePackCol, Npack)
        call face_interp_interior_set &
            (fHeadSetU, eHeadSet, er_InterpWeight_dH, er_InterpWeight_uH, facePackCol, Npack)
        call face_interp_interior_set &
            (fFlowSet, eFlowSet, er_InterpWeight_dQ, er_InterpWeight_uQ, facePackCol, Npack)

        !% copy upstream to downstream storage at a face
        !% (only for Head and Geometry types)
        !% note that these might be reset by hydraulic jump
        call face_copy_upstream_to_downstream_interior &
            (fGeoSetD, fGeoSetU, facePackCol, Npack)

        call face_copy_upstream_to_downstream_interior &
            (fHeadSetD, fHeadSetU, facePackCol, Npack)

        !% calculate the velocity in faces and put limiter
        call adjust_face_dynamic_limit (facePackCol, .true.)

        !% reset all the hydraulic jump interior faces
        call jump_compute

        !%------------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%face) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine face_interpolation_interior
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine face_interpolation_shared (facePackCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% Interpolates all the shared faces
        !%-------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: facePackCol  !% Column in faceP array for needed pack
            integer, intent(in) :: Npack        !% expected number of packed rows in faceP.
            integer :: fGeoSetU(3), fGeoSetD(3), eGeoSet(3)
            integer :: fHeadSetU(1), fHeadSetD(1), eHeadSet(1)
            integer :: fFlowSet(1), eFlowSet(1)
            character(64) :: subroutine_name = 'face_interpolation_shared'
        !%-------------------------------------------------------------------
        !% Preliminaries   
            if (icrash) return
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% Face values are needed for
        !% Area_u, Area_d, Head_u, Head_d, Flowrate,

        !% HACK not sure if we need
        !% Topwidth_u, Topwidth_d, HydDepth_u, HydDepth_d
        !% Velocity_u, Velocity_d

        !% General approach
        !% interpolate to ..._u
        !% identify hydraulic jumps
        !% set .._u and ..d based on jumps

        !% set the matching sets
        !% HACK THESE SHOULD BE DONE IN A GLOBAL -- MAYBE SETTINGS
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
        call face_interp_shared_set &
            (fGeoSetU, eGeoSet, er_InterpWeight_dG, er_InterpWeight_uG, facePackCol, Npack)
        call face_interp_shared_set &
            (fHeadSetU, eHeadSet, er_InterpWeight_dH, er_InterpWeight_uH, facePackCol, Npack)   
        call face_interp_shared_set &
            (fFlowSet, eFlowSet, er_InterpWeight_dQ, er_InterpWeight_uQ, facePackCol, Npack)

        !% copy upstream to downstream storage at a face
        !% (only for Head and Geometry types)
        !% note that these might be reset by hydraulic jump
        call face_copy_upstream_to_downstream_shared &
            (fGeoSetD, fGeoSetU, facePackCol, Npack)

        call face_copy_upstream_to_downstream_shared &
            (fHeadSetD, fHeadSetU, facePackCol, Npack)

        call adjust_face_dynamic_limit (facePackCol, .false.)
        !% HACK needs jump computation for across shared faces
        ! print *, "HACK missing hydraulic jump that occurs on shared faces 36987"

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
        !%-------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: fset(:), eset(:), eWdn, eWup, facePackCol, Npack
            integer, pointer :: thisP(:), eup(:), edn(:)
            integer :: ii
            character(64) :: subroutine_name = 'face_interp_interior_set'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (icrash) return
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% Aliases        
            thisP => faceP(1:Npack,facePackCol)
            eup   => faceI(:,fi_Melem_uL)
            edn   => faceI(:,fi_Melem_dL)
        !%--------------------------------------------------------------------
        !% cycle interpolation through each type in the set.
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
        !%-------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: fset(:), eset(:), eWdn, eWup, facePackCol, Npack
            integer, pointer :: thisP, eup, edn, connected_image, ghostUp, ghostDn
            logical, pointer :: isGhostUp, isGhostDn
            integer :: ii, jj   
            character(64) :: subroutine_name = 'face_interp_shared_set'
        !%--------------------------------------------------------------------
        !%  Preliminaries
            if (icrash) return
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%--------------------------------------------------------------------
        !% cycle through all the shared faces
        do ii = 1,Npack
            !%-----------------------------------------------------------------
            !% Aliases
            thisP           => facePS(ii,facePackCol)
            connected_image => faceI(thisP,fi_Connected_image)
            eup             => faceI(thisP,fi_Melem_uL)
            edn             => faceI(thisP,fi_Melem_dL)
            ghostUp         => faceI(thisP,fi_GhostElem_uL)
            ghostDn         => faceI(thisP,fi_GhostElem_dL)
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

                else
                    write(*,*) 'CODE ERROR: unexpected else'
                    stop 487874
                end if        
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

        !%--------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine face_interp_shared_set
!%
!%==========================================================================
!%==========================================================================
!%  
    subroutine face_head_average_on_element &
        (whichTM)
        !%------------------------------------------------------------------
        !% Description:
        !% computes the average head of the faces on an element
        !%-------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: whichTM 
            integer, pointer :: Npack, thisP(:), thisCol
            integer, pointer :: mapUp(:), mapDn(:)
            real(8), pointer :: fHeadU(:), fHeadD(:), eHeadAvg(:)
            character(64) :: subroutine_name = 'face_head_average_on_element'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (.not. setting%Solver%QinterpWithLocalHeadGradient) return  
            if (icrash) return
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            select case (whichTM)
                case (ALLtm)
                    thisCol => col_elemP(ep_CC_ALLtm)
                case (ETM)
                    thisCol => col_elemP(ep_CC_ETM)
                case (AC)
                    thisCol => col_elemP(ep_CC_AC)
                case (dummy)
                    print *, 'error - the QinterpWithLocalHeadGradient has not been coded for diagnostic elements'
                    stop 58704
                case default
                    print *, 'error, this default case should not be reached'
                    stop 2394
            end select         
        !%-------------------------------------------------------------------
        !% Aliases             
            Npack    => npack_elemP(thisCol)
            thisP    => elemP(1:Npack,thisCol)
            mapUp    => elemI(:,ei_Mface_uL)
            mapDn    => elemI(:,ei_Mface_dL)   
            fHeadU   => faceR(:,fr_Head_u)  
            fHeadD   => faceR(:,fr_Head_d)
            eHeadAvg => elemR(:,er_HeadAvg)    
        !%-------------------------------------------------------------------   
        !% The map up must use the downstream head on the face.
        !% The map dn must use the upstream head on the face.
        eHeadAvg(thisP) = onehalfR * (fHeadU(mapDn(thisP)) + fHeadD(mapUp(thisP)))

        ! print *, 'in ',trim(subroutine_name)
        ! print *, thisP
        ! print *, eHeadAvg(thisP)
        ! print *, ' '
        ! print *, eHeadAvg(:)
        ! print *, ' '
        ! print *, fHeadU(:)
        ! print *, ' '
        ! print *, fHeadD(:)
        ! print *, ' '
        !%-------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine face_head_average_on_element
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine face_Q_HeadGradient_interior &
        (faceCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% Adds the head gradient term to the face flowrate for interior faces
        !% should be done after Q and H are interpolated to face
        !% and element HeadAvg is computed.
        !%-------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: faceCol, Npack
            integer, pointer :: thisP(:), eup(:), edn(:)
            real(8), pointer :: fQ(:), eArea(:), eHead(:), eHeadAvg(:)
            real(8), pointer :: eLength(:), dt, grav 
            character(64) :: subroutine_name = 'face_interp_Q_HeadGradeint'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (.not. setting%Solver%QinterpWithLocalHeadGradient) return  
            if (icrash) return
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% Aliases   
            thisP    => faceP(1:Npack,faceCol)
            eup      => faceI(:,fi_Melem_uL)
            edn      => faceI(:,fi_Melem_dL)
            fQ       => faceR(:,fr_Flowrate)
            eArea    => elemR(:,er_Area)
            eHead    => elemR(:,er_Head)
            eHeadAvg => elemR(:,er_HeadAvg)
            eLength  => elemR(:,er_Length)
            dt       => setting%Time%Hydraulics%Dt
            grav => setting%constant%gravity

        ! print *, ' '
        ! print *, 'in dQ term'
        ! print *,  'head comparison'
        ! print *, ' 1-3'
        ! print *,  eHead(1:3)
        ! print *,  eHeadAvg(1:3)
        ! print *, ' 4-6' 
        ! print *,  eHead(4:6)
        ! print *,  eHeadAvg(4:6)
        ! print *, ' 7'
        ! print *,  eHead(7)
        ! print *,  eHeadAvg(7)
        ! print *, ' ' 
        !%-------------------------------------------------------------------
        !% adds term dt * grav A [ (dh/dx) - (dh_avg/dx) ]
        fQ(thisP) = fQ(thisP) + dt * grav *                                         &
            (                                                                       &
                +( eArea(eup(thisP)) * ( eHead(eup(thisP)) - eHeadAvg(eup(thisP)) ) &
                    / ( onehalfR * eLength(eup(thisP) ) ) )                         &
                -( eArea(edn(thisP)) * ( eHead(edn(thisP)) - eHeadAvg(edn(thisP)) ) &
                    / ( onehalfR * eLength(edn(thisP) ) ) )                         &
            ) 

        !%-------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine face_Q_HeadGradient_interior
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine face_Q_HeadGradient_shared &
        (faceCol, Npack) 
        !%------------------------------------------------------------------
        !% Description:
        !% Adds the head gradient term to the face flowrate for shared faces
        !% should be done after Q and H are interpolated to face
        !% and element HeadAvg is computed.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: faceCol, Npack
            integer, pointer :: thisP, eup, edn, ci, ghostUp, ghostDn
            logical, pointer :: isGhostUp, isGhostDn
            integer :: ii, jj   
            real(8), pointer :: fQ(:)
            real(8), pointer :: eLength(:), dt, grav
            character(64) :: subroutine_name = 'face_Q_HeadGradient_shared'
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (.not. setting%Solver%QinterpWithLocalHeadGradient) return  
            if (icrash) return
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases
                fQ   => faceR(:,fr_Flowrate)
                dt   => setting%Time%Hydraulics%Dt     
                grav => setting%constant%gravity   
        !%------------------------------------------------------------------
        !% Cycle over individual faces is required because of Coarrays        
        do ii = 1,Npack
            !%-----------------------------------------------------------------
            !% Local Aliases
                thisP       => facePS(ii,faceCol)
                ci          => faceI(thisP,fi_Connected_image) !% Connected_image
                eup         => faceI(thisP,fi_Melem_uL)
                edn         => faceI(thisP,fi_Melem_dL)
                ghostUp     => faceI(thisP,fi_GhostElem_uL)
                ghostDn     => faceI(thisP,fi_GhostElem_dL)
                isGhostUp   => faceYN(thisP,fYN_isUpGhost)
                isGhostDn   => faceYN(thisP,fYN_isDnGhost)
            !%-----------------------------------------------------------------
            !% upstream element of the shared face is ghost and in a different image
            if (isGhostUp) then
                fQ(thisP) = fQ(thisP) + dt * grav *                                           &
                    (                                                                         &
                        +( elemR(ghostUp,er_Area)[ci]                                         &
                            * ( elemR(ghostUp,er_Head)[ci] - elemR(ghostUp,er_HeadAvg)[ci] )  &
                            / ( onehalfR * elemR(ghostUp,er_Length)[ci] ) )                   &
                        -( elemR(edn,er_Area)                                                 &
                            * ( elemR(edn,er_Head)         - elemR(edn,er_HeadAvg)         )  &
                            / ( onehalfR * elemR(edn,er_Length)         ) )                   &
                    )

            !% downstream element of the shared face is ghost and in a different image
            elseif (isGhostDn) then
                fQ(thisP) = fQ(thisP) + dt * grav *                                           &
                    (                                                                         &
                        +( elemR(eup,er_Area)                                                 &
                            * ( elemR(eup,er_Head)         - elemR(eup,er_HeadAvg)         )  &
                            / ( onehalfR * elemR(eup,er_Length) ) )                           &
                        -( elemR(ghostDn,er_Area)[ci]                                         &
                            * ( elemR(ghostDn,er_Head)[ci] - elemR(ghostDn,er_HeadAvg)[ci] )  &
                            / ( onehalfR * elemR(ghostDn,er_Length)[ci] ) )                   &
                    )
            else
                write(*,*) 'CODE ERROR: unexpected else'
                stop 488758
            end if        
        end do        
        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%face) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine face_Q_HeadGradient_shared
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
        if (icrash) return
        if (setting%Debug%File%face) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------

        thisP => faceP(1:Npack,facePackCol)

        faceR(thisP,downstreamSet) = faceR(thisP,upstreamSet)

        if (setting%Debug%File%face) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine face_copy_upstream_to_downstream_interior
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine face_copy_upstream_to_downstream_shared &
        (downstreamSet, upstreamSet, facePackCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Copies the interpolated value on the upstrea side to the downstream side
        !% These values are later adjusted for hydraulic jumps
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: facePackCol, Npack, downstreamSet(:), upstreamSet(:)
        integer, pointer :: thisP(:)
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'face_copy_upstream_to_downstream'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%face) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------

        thisP => facePS(1:Npack,facePackCol)

        faceR(thisP,downstreamSet) = faceR(thisP,upstreamSet)

        if (setting%Debug%File%face) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine face_copy_upstream_to_downstream_shared
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
        
        character(64) :: subroutine_name = 'face_interp_set_byMask'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%face) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
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
        end do

        print *, 'in face_interp_set_byMask -- may be obsolete'
        stop 87098

        if (setting%Debug%File%face) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine face_interp_set_byMask
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module face
