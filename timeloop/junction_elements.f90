module junction_elements

    use define_globals
    use define_keys
    use define_indexes
    use define_xsect_tables
    use define_settings, only: setting
    use utility, only: util_CLprint
    use utility_crash, only: util_crashpoint

!%----------------------------------------------------------------------------- 
!% Description:
!% Computes junction elements
!%----------------------------------------------------------------------------- 

    implicit none

    private
    
    public :: junction_toplevel

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine junction_toplevel (whichTM, istep)
        !%-----------------------------------------------------------------
        !% Description:
        !% Controls computation of implicit junction element
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: whichTM, istep
            integer             :: thisColP, Npack

        !%-----------------------------------------------------------------
        !% Preliminaries
            select case (whichTM)
                case (ALLtm)
                    thisColP = col_elemP(ep_JM_ALLtm)
                case (ETM)
                    thisColP = col_elemP(ep_JM_ETM)
                case (AC)
                    thisColP = col_elemP(ep_JM_AC)
                case default
                    print *, 'CODE ERROR: time march type unknown for # ', whichTM
                    print *, 'which has key ',trim(reverseKey(whichTM))
                    stop 7659
            end select

            Npack = npack_elemP(thisColP)

            if (Npack == 0) return
        !%-----------------------------------------------------------------

        call junction_calculation (Npack, thisColP, istep)

    end subroutine junction_toplevel
!%    
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine  junction_calculation (Npack, thisColP, istep)
        !%-----------------------------------------------------------------
        !% Description:
        !% get the jacobean matrix for junction element following the 
        !% derivation by ...
        !%-----------------------------------------------------------------
        !% Declarations
        integer, intent(in) :: Npack, thisColP, istep
        integer, pointer    :: thisJM(:)

        real(8), pointer :: grav

        !% --- the junction data storage
        real(8) :: jDataR(Npack, max_branch_per_node, NCol_jDataR) 
        integer :: jDataI(Npack, max_branch_per_node, NCol_jDataI)

        !% --- solver data storage
        real(8), dimension(Npack) :: startEJ, EJ, Qlat, VolumeEJ, fullVolumeEJ
        real(8), dimension(Npack) :: Zcrown, Aplan, Aponded
        real(8), dimension(Npack) :: OverflowLength, OverflowHeight, Qoverflow

        real(8), dimension(Npack) :: alpha, normL2, Qin, Qout,  fcoef
        real(8), dimension(Npack) :: deltaEJ, residEJ, Volume

        integer, dimension(Npack) :: np, Niter, OverFlowType, JunctionType

        logical, dimension(Npack) :: isconverged, isRepack, isfailed, isFull, isSubmergedOrifice
        integer :: kk, mm

        !% --- epsilon used in junction to prevent zero values
        real(8) :: jEpsilon = 1.d-6

        !% --- maximum number of iterations (HACK -- make into global)
        integer :: maxIter = 10
        !% --- flowrate conservation convergence (HACK -- make into global)
        real(8) :: jConvergence = 1.d-4

        !%-----------------------------------------------------------------------------
        !%  Aliases:         
            grav         => setting%constant%gravity
            !% --- packed set of junctions
            thisJM => elemP(1:Npack,thisColP)
        !%-----------------------------------------------------------------------------
        !% Preliminaries
            !% --- iteration counter
            Niter = zeroI

            !% --- initialization
            alpha   = zeroR
            normL2  = zeroR 
            Qin     = zeroR 
            Qout    = zeroR 
            fcoef   = zeroR
            deltaEJ = zeroR
            residEJ = zeroR
            Qoverflow = zeroR


            !% --- stores the initial elemR k location
            !%     this is used to bring packed branch data back to element and face array
            do concurrent (kk=1:max_branch_per_node)
                jDataI(:,kk,ji_kidx)  = kk  
            end do

            !% store the starting head in Temp array for conservation check at end
            elemR(thisJM,er_Temp01) = elemR(thisJM,er_Head)
        !%-----------------------------------------------------------------------------

        !% --- get the junction data 
        call junction_get_JM_data (startEJ, EJ, Qlat, volumeEJ, fullVolumeEJ, &
                Zcrown, Aplan, Aponded, OverflowLength, OverFlowHeight, &
                isFull, isSubmergedOrifice, JunctionType, OverflowType, thisJM, Npack)

        !% --- store the face and elem data in 3D array
        call junction_get_face_data(jDataR, thisJM)

        !% --- set the alpha at the junction
        call junction_alpha (alpha, Aplan, Aponded, isFull, JunctionType, &
             OverflowType, thisJM, Npack, istep)

        !% --- set the fcoef for overflows
        call junction_fcoef(fcoef, Aplan, OverflowLength, OverflowType)

        !% AFTER THIS POINT WE SHOULD NOT NEED TO USE thisJM(:) UNTIL
        !% MOVING DATA BACK TO faceR AND elemR

        !% --- identify branches with gamma = 0 and K=0 or gamma = 0 and Q = 0
        !%     requires setting beta to zero (thus, we cannot compute flow in 
        !%     this branch during this time step)
        call junction_initialize_beta (jDataR, jEpsilon)

        !% --- the expected total branches
        !%     This will be modified by packing if there are branches with beta=0.
        np = elemSI(thisJM,esi_JunctionMain_Total_Branches)     
        
        !% --- code check for debugging
        call junction_check_beta (jDataR, np, nPack)

        !% --- pack the jDataR and jDataI to remove zero branches
        isRepack(:) = .true.
        call junction_pack_all (jDataR, jDataI, np, isRepack, Npack)

        !% --- compute the invariant terms that will not change in this time step
        call junction_invariant_terms (jDataR, np, Npack)

        !% --- compute the residual in each branch
        call junction_all_residuals (jDataR, residEJ, EJ, &
            startEJ, alpha, fcoef, Qlat, Zcrown, OverflowHeight, &
            OverflowType, isFull, isSubmergedOrifice, np, Npack)

        isconverged = .false. 
        isfailed    = .false.

        do concurrent (mm=1:Npack)
            do while (.not. isconverged(mm))
                !% --- increment iteration counter for exit after maxIter loops
                Niter(mm) = Niter(mm) + 1

                !% --- compute varying lambdaC 
                jDataR(mm,:,jr_LambdaC) = junction_lambdaC (jDataR(mm,:,:), np(mm))

                !% --- look for small LambdaC and repack if needed
                call junction_check_repack_oneJ (jDataR(mm,:,:), np(mm), jEpsilon, isRepack(mm))
          
                if (isRepack(mm)) then 
                    call junction_pack_oneJ &
                        (jDataR(mm,:,:), jDataI(mm,:,:), np(mm), isRepack(mm))
                end if

                !% --- compute the junction iteration
                call junction_iteration &
                    (jDataR(mm,:,:), deltaEJ(mm), EJ(mm),  residEJ(mm), alpha(mm), &
                     fcoef(mm), Zcrown(mm), OverflowHeight(mm), OverflowType(mm), isFull(mm), &
                     isSubmergedOrifice(mm), np(mm)) 

                !% --- compute the residuals
                jDataR(mm,1:np(mm),jr_residQ) = junction_Qresidual              &
                                             (jDataR(mm,:,:), EJ(mm), startEJ(mm), np(mm))
                         
                residEJ(mm) = junction_Eresidual                                           &
                             (jDataR(mm,:,:), EJ(mm), startEJ(mm), alpha(mm), fcoef(mm),   &
                              Qlat(mm), Zcrown(mm), OverflowHeight(mm), isFull(mm),        &
                              isSubmergedOrifice(mm), OverflowType(mm), np(mm))

                !% --- compute residual L2 norm
                normL2(mm) = junction_get_L2norm &
                            (jDataR(mm,:,jr_residQ), residEJ(mm), np(mm))

                !% --- check norms for exit
                if (normL2(mm) .le. jConvergence) then 
                    isconverged(mm) = .true.
                else
                    if (Niter(mm)+1 .ge. maxIter) then 
                        !% --- mark as failed and exit
                        isconverged(mm) = .true.
                        isfailed(mm)    = .false.
                    end if
                end if 
            end do

            !% ---- CONVERGED OR EXITING FOR THIS JUNCTION------------------------------
            !% --- compute the net Q inflow
            Qin(mm) = junction_branch_Qnet (jDataR(mm,:,:), np(mm), +oneI)
            !% --- compute the net Q outflow
            Qout(mm)= junction_branch_Qnet (jDataR(mm,:,:), np(mm), -oneI)
            !% --- compute the Q overflow
            Qoverflow(mm) = junction_Qoverflow (EJ(mm), Zcrown(mm), fcoef(mm),            &
                OverflowHeight(mm), OverflowType(mm), isFull(mm), isSubmergedOrifice(mm) )

            !% --- adjust Qs to ensure conservation
            jDataR(mm,1:np(mm),jr_Q) = junction_conservation &
                (jDataR(mm,:,:), residEJ(mm), Qin(mm), Qout(mm), np(mm))

            !% --- restore data to the faceR ... arrays
            call junction_push_face_flowrate &
                (faceR(:,fr_Flowrate), faceR(:,fr_Velocity_u), faceR(:,fr_Velocity_d), &
                 faceR(:,fr_Area_u), faceR(:,fr_Area_d), jDataR(mm,:,:), jDataI(mm,:,:), &
                 np(mm), thisJM(mm))

            !% --- head
            elemR(thisJM(mm),er_Head) = EJ(mm)

            !% --- volume
            elemR(thisJM(mm),er_Volume) = junction_volume &
                (elemR(thisJM(mm),er_Volume_N0), Qin(mm), Qout(mm), Qlat(mm), Qoverflow(mm), istep) 

            !% --- set the junction overflow and ponded volume on the last step of RK
            if (istep == 2) then
                !% --- overflow
                elemR(thisJM(mm),er_VolumeOverFlow) = Qoverflow(mm) * setting%Time%Hydraulics%Dt

                !% --- ponded volume
                if ((OverflowType(mm) == Ponded) .and. (isFull(mm))) then 
                    elemR(thisJM(mm),er_VolumePonded) = elemR(thisJM(mm),er_Volume) - elemR(thisJM(mm),er_FullVolume)
                end if
            end if

        end do

        !% --- check conservation -- can be commented after debugging
        call junction_check_conservation (thisJM, Qoverflow, Npack ,istep)

        !% HACK --- the new water surface elevation associated with the
        !% changing storage may NOT exactly balance conservation because of
        !% the use of a fixed planar area. We need to rederiv the 
        !% junction_conservation equations so that the exact change in volume
        !% implied by the change in head gives us the correct volume

        !% --- HACK: NEED TO DO SOMETHING FOR REPORTING CONVERGENCE FAILURE 

    end subroutine junction_calculation
!%
!%========================================================================== 
!% PRIVATE
!%========================================================================== 
!% 
    pure subroutine junction_get_JM_data &
        (startEJ, EJ, Qlat, volumeEJ, fullVolumeEJ, Zcrown, Aplan, Aponded, &
            OverflowLength, OverflowHeight, &
            isFull, isSubmergedOrifice, JunctionType, OverflowType, thisJM, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% Gets the JM data from the elemR array and stores in local
        !% vector
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(inout) :: startEJ(:), EJ(:), Qlat(:)
            real(8), intent(inout) :: volumeEJ(:), fullVolumeEJ(:), Zcrown(:)
            real(8), intent(inout) :: Aplan(:), Aponded(:)
            real(8), intent(inout) :: OverflowLength(:), OverflowHeight(:)
            logical, intent(inout) :: isFull(:), isSubmergedOrifice(:)
            integer, intent(inout) :: OverflowType(:), JunctionType(:)
            integer, intent(in)    :: Npack, thisJM(:)
        !%------------------------------------------------------------------

        startEJ     (1:Npack) = elemR(thisJM,er_Head)
        EJ          (1:Npack) = elemR(thisJM,er_Head)
        Qlat        (1:Npack) = elemR(thisJM,er_FlowrateLateral)
        volumeEJ    (1:Npack) = elemR(thisJM,er_Volume)
        fullVolumeEJ(1:Npack) = elemR(thisJM,er_FullVolume)
        Zcrown      (1:Npack) = elemR(thisJM,er_Zcrown)
        OverflowType(1:Npack) = elemSI(thisJM,esi_JunctionMain_OverflowType)
        Aponded     (1:Npack) = elemSR(thisJM,esr_JunctionMain_PondedArea)
        Aplan       (1:Npack) = elemSR(thisJM,esr_Storage_Plane_Area)
        JunctionType(1:Npack) = elemSI(thisJM,esi_JunctionMain_Type)

        OverflowLength(1:Npack) = elemSR(thisJM,esr_JunctionMain_OverflowOrifice_Length)
        OverflowHeight(1:Npack) = elemSR(thisJM,esr_JunctionMain_OverflowOrifice_Height)

        where (JunctionType .eq. ImpliedStorage)
            !% --- The only use of Aplan with implied storage is for the
            !%     overflow weir effect.
           Aplan = 0.179d0 !% Area gives Lw = 1.5 m as overflow weir length
        endwhere

        where (volumeEJ .ge. fullVolumeEJ)
            isFull = .true.
        elsewhere 
            isFull = .false.
        endwhere

        where ((isFull) .and. (OverflowType .eq. OverflowOrifice) &
            .and. (EJ > Zcrown + OverflowHeight))
            isSubmergedOrifice = .true.
        elsewhere 
            isSubmergedOrifice = .false.
        endwhere

   
    end subroutine junction_get_JM_data
!%
!%========================================================================== 
!%========================================================================== 
!%    
    pure subroutine junction_get_face_data (jDataR, thisJM )    
        !%------------------------------------------------------------------
        !% Takes the data from faceR and elemR arrays and stores in the
        !% jDataR 3D array
        !%------------------------------------------------------------------
            real(8), intent(inout) :: jDataR(:,:,:)
            integer, intent(in)    :: thisJM(:)

            integer :: kk
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        !% --- store data for upstream faces
        do concurrent (kk=1:max_branch_per_node:2)
            where (elemSI(thisJM+kk,esi_JunctionBranch_Exists) == oneI)
                jDataR(:,kk,jr_beta)    = +oneR !% +1 for upstream beta
                jDataR(:,kk,jr_Area)    =  faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Area_u)
                jDataR(:,kk,jr_Ebranch) =  faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Head_u)             &
                                        + (faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Velocity_u)**twoI) &
                                        / (twoR * setting%constant%gravity)
                jDataR(:,kk,jr_Gamma)   =  faceR(elemI(thisJM+kk,ei_Mface_uL),fr_GammaM)
                jDataR(:,kk,jr_Length)  =  faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Length_u)
                jDataR(:,kk,jr_Q)       =  faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Flowrate)
                jDataR(:,kk,jr_K)       =  faceR(elemI(thisJM+kk,ei_Mface_uL),fr_KJunction_MinorLoss)
            elsewhere
                jDataR(:,kk,jr_beta)    = zeroR !%  0 if branch does not exist
            endwhere
        end do

        !% --- store data for downstream faces
        do concurrent (kk=2:max_branch_per_node:2)
            where (elemSI(thisJM+kk,esi_JunctionBranch_Exists) == oneI)
                jDataR(:,kk,jr_beta)    = -oneR !% -1 for downstream beta
                jDataR(:,kk,jr_Area)    =  faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Area_d)
                jDataR(:,kk,jr_Ebranch) =  faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Head_d)             &
                                        + (faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Velocity_d)**twoI)  &
                                        / (twoR * setting%constant%gravity)
                jDataR(:,kk,jr_Gamma)   =  faceR(elemI(thisJM+kk,ei_Mface_dL),fr_GammaM)
                jDataR(:,kk,jr_Length)  =  faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Length_d)
                jDataR(:,kk,jr_Q)       =  faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Flowrate)
                jDataR(:,kk,jr_K)       =  faceR(elemI(thisJM+kk,ei_Mface_dL),fr_KJunction_MinorLoss)
            elsewhere
                jDataR(:,kk,jr_beta)    = zeroR !%  0 if branch does not exist
            endwhere
        end do


    end subroutine junction_get_face_data
!%
!%========================================================================== 
!%========================================================================== 
!%      
    subroutine junction_alpha (alpha, Aplan, Aponded, isFull, JunctionType,  &
                OverflowType, thisJM, Npack, istep)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the alpha term that represents change in junction storage
        !%------------------------------------------------------------------
        !% Declarations
            
            real(8), intent(inout) :: alpha(:)
            real(8), intent(in)    :: Aplan(:), Aponded(:)
            integer, intent(in)    :: JunctionType(:), OverflowType(:), thisJM(:), Npack, istep
            logical, intent(in)    :: isFull(:)

            real(8), pointer :: crk(:), dt
            integer          :: mm
        !%------------------------------------------------------------------
        !% Aliases
            crk          => setting%Solver%crk2
            dt           => setting%Time%Hydraulics%Dt
        !%------------------------------------------------------------------

        do mm=1,Npack
            select case (JunctionType(mm))

                case (ImpliedStorage)
                    !% --- no volume increase allowed
                    alpha(mm) = zeroR

                case (FunctionalStorage, TabularStorage)
                    if (.not. isFull(mm)) then
                        !% --- below full volume
                        alpha(mm) = Aplan(mm) / (crk(istep)*dt)
                    else
                        !% --- at or above full volume
                        select case (OverflowType(mm))
                            case (NoOverflow, OverflowOrifice, OverflowWeir)
                                !% --- alpha remains based on plane area
                                alpha(mm) = Aplan(mm) / (crk(istep)*dt)
                            case (Ponded)
                                !% --- alpha used ponded area instead of plane area
                                alpha(mm) = Aponded(mm)  / (crk(istep)*dt)
                            case default 
                                print *, 'CODE ERROR: Unexpected case default'
                                call util_crashpoint(2287422)
                        end select
                    end if
                case default
                    print *, 'unexpected case default'
                    call util_crashpoint(628733)
            end select  
        end do

    end subroutine junction_alpha
!%
!%========================================================================== 
!%========================================================================== 
!%   
    pure subroutine junction_fcoef (fcoef, Aplan, OverflowLength, OverflowType)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the fcoef term for overflow weirs
        !%------------------------------------------------------------------
        !% Declarations   
            real(8), intent(inout) :: fcoef(:)
            real(8), intent(in)    :: Aplan(:), OverflowLength(:)
            integer, intent(in)    :: OverflowType(:)
            real(8), parameter :: Cbc = 1.45d0  !% weir coef. HACK MOVE TO SETTINGS
        !%------------------------------------------------------------------

        where (OverflowType == OverflowWeir)
            fcoef = twoR * Cbc * sqrt( setting%Constant%pi * Aplan)
        elsewhere (OverflowType == OverflowOrifice)
            fcoef = OverflowLength * sqrt(twoR * setting%Constant%pi)
        elsewhere
            fcoef = zeroR
        endwhere

    end subroutine junction_fcoef    
!%
!%========================================================================== 
!%========================================================================== 
!%   
    pure subroutine junction_initialize_beta (jDataR, jEpsilon)
        !%------------------------------------------------------------------
        !% Description:
        !% sets the initial jDataR(:,jr_beta) = 0 for branches that
        !% cannot have a flow in this time step.
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(inout) :: jDataR(:,:,:)
            real(8), intent(in)    :: jEpsilon
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        where (((jDataR(:,:,jr_Gamma)     < jEpsilon) .and.   &
                (jDataR(:,:,jr_K)         < jEpsilon)       ) &
            .or.                                              &
               ((jDataR(:,:,jr_Gamma)     < jEpsilon) .and.   & 
                (abs(jDataR(:,:,jr_Q))    < jEpsilon)       ) ) 

            jDataR(:,:,jr_beta) = zeroR

        endwhere

    end subroutine junction_initialize_beta
!%
!%========================================================================== 
!%========================================================================== 
!%   
    subroutine junction_check_beta (jR, np, nPack)
        !%------------------------------------------------------------------
        !% Description:
        !% Checks that np(:) is identical to the number of beta <> 0
        !% This can be commented out once the code is thoroughly debugged
        !%------------------------------------------------------------------
            real(8), intent(in) :: jR(:,:,:)
            integer, intent(in) :: np(:), nPack 
            integer :: mm
            real(8) :: nbeta

        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        do mm=1,Npack 
            nbeta = count(jR(mm,:,jr_beta) .ne. zeroR)
            if (np(mm) .ne. nbeta) then 
                print *, 'CODE ERROR:'
                print *, 'problem with np and nbeta in junction iteration'
                print *, 'mm     = ',mm 
                print *, 'np(mm) = ',np(mm)
                print *, 'n beta = ',nbeta
                call util_crashpoint(398233)
            end if
        end do

    end subroutine junction_check_beta
!%
!%========================================================================== 
!%========================================================================== 
!%
    pure subroutine junction_invariant_terms (jR, np, Npack)
        !%------------------------------------------------------------------
        !% Description
        !% computes the terms in the junction solution that do not change
        !% during a time step
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(inout) :: jR(:,:,:)
            integer, intent(in)    :: np(:), nPack

            integer :: mm
        !%------------------------------------------------------------------

        do concurrent (mm=1:Npack)
            !% --- a
            jR(mm,1:np(mm),jr_a) = jR(mm,1:np(mm),jr_beta)             &
                * jR(mm,1:np(mm),jr_Length) * jR(mm,1:np(mm),jr_Gamma) &
                / (twoR * setting%constant%gravity * jR(mm,1:np(mm),jr_Area))

            !% --- c
            jR(mm,1:np(mm),jr_c) = jR(mm,1:np(mm),jr_beta)     &
                * jR(mm,1:np(mm),jr_K)                             &
                / (twoR * setting%constant%gravity * (jR(mm,1:np(mm),jr_Area)**twoI))   

            !% --- lambdaA
            jR(mm,1:np(mm),jr_LambdaA)                            &
                = twoR * setting%constant%gravity  *  (jR(mm,1:np(mm),jr_Area)**twoI) 

            !% --- lambdaB
            jR(             mm,1:np(mm),jr_LambdaB)               &
                =        jR(mm,1:np(mm),jr_Area)                  &
                       * jR(mm,1:np(mm),jr_Length)                &
                       * jR(mm,1:np(mm),jr_Gamma)
        end do

    end subroutine junction_invariant_terms
!%
!%========================================================================== 
!%========================================================================== 
!%   
    pure function junction_lambdaC (jR, np)
        !%------------------------------------------------------------------
        !% Description
        !% computes the lambdaC term that changes with each iteration
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(in)    :: jR(:,:)
            integer, intent(in)    :: np

            real(8), dimension(np) :: junction_lambdaC

            integer :: mm
        !%------------------------------------------------------------------

        !% --- note that the LambdaC storage is later used for Lambda = LambdaA/LambdaC
        junction_lambdaC = jR(1:np,jr_LambdaB) + jR(1:np,jr_K) * abs(jR(1:np,jr_Q))

    end function junction_lambdaC
!%
!%========================================================================== 
!%========================================================================== 
!%   
    pure subroutine junction_iteration (jR, deltaEJ, EJ, residEJ, alpha, &
            fcoef, Zcrown, OverflowHeight, OverflowType, isFull, isSubmergedOrifice, np) 
        !%------------------------------------------------------------------
        !% Description
        !% Iterative results for one step advancing Q and EJ
        !%------------------------------------------------------------------
            real(8), intent(inout) :: jR(:,:), deltaEJ, EJ
            real(8), intent(in)    :: residEJ, alpha, fcoef, Zcrown, OverflowHeight
            logical, intent(in)    :: isFull, isSubmergedOrifice
            integer, intent(in)    :: np, OverFlowType

            real(8) :: Bsum, Csum, Gterm

            integer :: mm
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        !% --- compute lambda and store in lambdaC (recomputed each iteration)
        jR(1:np,jr_LambdaC) =  jR(1:np,jr_LambdaA)  / jR(1:np,jr_LambdaC)
        
        !% B summation
        Bsum = sum(jR(1:np,jr_LambdaC)) 
        
        !% C summation
        Csum = sum(jR(1:np,jr_LambdaC) * jR(1:np,jr_residQ))

        !% Gterm
        Gterm = junction_Gterm (alpha, EJ, Zcrown, fcoef, OverflowHeight, &
                                OverflowType, isFull, isSubmergedOrifice)

        !% change in Junction Energy
        deltaEJ = (residEJ - Csum ) / (-Gterm - Bsum)
    
        !% change in flowrate (recall LambdaC was overwritten with LambdaI)
        jR(1:np,jr_DeltaQ) = jR(1:np,jr_LambdaC) * (jR(1:np,jr_residQ) - deltaEJ) &
                           / jR(1:np,jr_beta)

        !%  New Q flowrate
        jR(1:np,jr_Q) = jR(1:np,jr_Q) + jR(1:np,jr_DeltaQ)

        !% --- update junction head
        EJ = EJ + deltaEJ

    end subroutine junction_iteration   
!%
!%========================================================================== 
!%========================================================================== 
!%    
    pure subroutine junction_all_residuals &
            (jR, residEJ, EJ, startEJ, alpha, fcoef, Qlat, &
             Zcrown, OverflowHeight, OverflowType, isFull, isSubmergedOrifice, np, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the residuals of the nonlinear junction equations for
        !% all junctions
        !%------------------------------------------------------------------
            real(8), intent(inout) :: jR(:,:,:), residEJ(:)
            real(8), intent(in)    :: EJ(:), startEJ(:), alpha(:), Qlat(:)
            real(8), intent(in)    :: Zcrown(:), fcoef(:), OverflowHeight(:)
            logical, intent(in)    :: isSubmergedOrifice(:), isFull(:)
            integer, intent(in)    :: OverflowType(:), np(:), Npack 

            integer :: mm
        !%------------------------------------------------------------------

        do concurrent (mm=1:Npack)
            !% --- get the flowrate residuals
            jR(mm,1:np(mm),jr_residQ) = junction_Qresidual                  &
                (jR(mm,:,:), EJ(mm), startEJ(mm),  np(mm))

            residEJ(mm) = junction_Eresidual                                &
                (jR(mm,:,:), EJ(mm), startEJ(mm), alpha(mm), fcoef(mm),     &
                Qlat(mm), Zcrown(mm), OverflowHeight(mm), isFull(mm),       &
                isSubmergedOrifice(mm), OverflowType(mm), np(mm))
        end do
        
    end subroutine junction_all_residuals
!%
!%========================================================================== 
!%========================================================================== 
!%   
    real(8) pure function junction_Eresidual &
            (jR, EJ, startEJ, alpha, fcoef, Qlat, Zcrown, OverflowHeight, &
             isFull, isSubmergedOrifice, OverflowType, np)  
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the E residual of the nonlinear junction equation for
        !% a single junction
        !%------------------------------------------------------------------
            real(8), intent(in)    :: jR(:,:)
            real(8), intent(in)    :: EJ, startEJ, alpha, fcoef, Qlat, Zcrown
            real(8), intent(in)    :: OverflowHeight
            logical, intent(in)    :: isSubmergedOrifice, isFull
            integer, intent(in)    :: np, OverflowType
        !%------------------------------------------------------------------

        !% --- get the junction energy residual
        if (.not. isFull) then
            !% --- junction that is not overflowing
            junction_Eresidual = -alpha * (EJ - startEJ) + Qlat &
                                + sum( (jR(1:np,jr_beta) * jR(1:np,jr_Q)) )  
        else 
            select case (OverflowType)

                case (NoOverflow)
                    !% --- Zcrown becomes volume limit
                    junction_Eresidual = -alpha * (Zcrown - startEJ) + Qlat     &
                                    + sum( (jR(1:np,jr_beta) * jR(1:np,jr_Q)) )

                case (OverFlowWeir)
                    !% --- requires additional outflow term
                    junction_Eresidual = -alpha * (EJ - startEJ) + Qlat            & 
                                    - twoR * fcoef * (EJ - Zcrown)**threehalfR     &
                                    + sum( (jR(1:np,jr_beta) * jR(1:np,jr_Q)) )

                case (OverflowOrifice)
                    if (isSubmergedOrifice) then 
                        junction_Eresidual = -alpha * (Zcrown - startEJ)                                     &
                                        - twothirdR * fcoef * ((EJ -  Zcrown                  )**threehalfR) &
                                        + twothirdR * fcoef * ((EJ - (Zcrown + OverflowHeight))**threehalfR) &
                                        + sum( (jR(1:np,jr_beta) * jR(1:np,jr_Q)) )
                    else
                        junction_Eresidual = -alpha * (Zcrown - startEJ)                  &
                                        - twothirdR * fcoef * ((EJ - Zcrown)**threehalfR) &
                                        + sum( (jR(1:np,jr_beta) * jR(1:np,jr_Q)) )
                    end if
                case (Ponded)
                    !% --- same residual as not full (alpha has been changed)
                    junction_Eresidual = -alpha * (EJ - startEJ) + Qlat          &
                                    + sum( (jR(1:np,jr_beta) * jR(1:np,jr_Q)) )  
            
            case default
                !% --- change to non-pure to have error checking here
            end select
        end if

    end function junction_Eresidual
!%
!%========================================================================== 
!%========================================================================== 
!%   
    pure function junction_Qresidual (jR, EJ, startEJ, np)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the flowrate residual for all the non-zero branches
        !% of a single junction
        !%------------------------------------------------------------------
            real(8), intent(in) :: jR(:,:), EJ, startEJ
            integer, intent(in) :: np
            real(8), dimension(np) :: junction_Qresidual
        !%------------------------------------------------------------------
    
        junction_Qresidual(1:np) = jR(1:np,jr_a) * jr(1:np,jr_Q)        &
            + jr(1:np,jr_c) * sign((jR(1:np,jr_Q)**twoI),jR(1:np,jr_Q)) &
            + EJ - startEJ

    end function junction_Qresidual
!%
!%========================================================================== 
!%========================================================================== 
!%   
    real(8) pure function junction_Gterm &
            (alpha, EJ, Zcrown, fcoef, OverflowHeight, OverflowType, isFull, isSubmergedOrifice)    
        !%------------------------------------------------------------------
        !% Description
        !% Computes G term for various overflow (or not overflow) conditions
        !%------------------------------------------------------------------
            real(8), intent(in) :: alpha, EJ, Zcrown, fcoef, OverflowHeight
            integer, intent(in) :: OverflowType
            logical, intent(in) :: isFull, isSubmergedOrifice
        !%------------------------------------------------------------------
        if (.not. isFull) then 
            !% --- EJ < Zcrown
            junction_Gterm = alpha
        else 
            !% --- overflowing conditions
            select case (OverflowType)

                case (NoOverflow)
                    !% --- surcharge without overflow
                    junction_Gterm = zeroR

                case (OverflowOrifice)
                    if (isSubmergedOrifice) then 
                        junction_Gterm = fcoef                           &
                            * (                                         &
                                + sqrt(EJ - Zcrown)                     &
                                - sqrt(EJ - (Zcrown + OverflowHeight))  &
                                )
                    else
                        junction_Gterm = fcoef * sqrt(EJ - Zcrown)
                    end if

                case (OverflowWeir)
                    junction_Gterm = alpha + threeR * fcoef * sqrt(EJ - Zcrown)

                case (Ponded)
                    !% --- note alpha is from the ponded area
                    junction_Gterm = alpha
                case default 
                    !% --- code should not be here
                    !%     to have a failure for testing 
                    !%     the function must be unpure    
            end select      
        end if

    end function junction_Gterm
!%
!%========================================================================== 
!%========================================================================== 
!%   
    ! real(8) pure function junction_Qoverflow_residual &
    !         (fcoef, EJ, Zcrown, Qoverflow)    
    !     !%------------------------------------------------------------------
    !     !% Description
    !     !% Computes the overflow flowrate residual\
    !     !%------------------------------------------------------------------
    !         real(8), intent(in) :: EJ, Zcrown, fcoef, Qoverflow
    !     !%------------------------------------------------------------------

    !     junction_Qoverflow_residual = fcoef * ((EJ - Zcrown)**(1.5d0)) + Qoverflow

    ! end function junction_Qoverflow_residual
!%
!%========================================================================== 
!%========================================================================== 
!%   
    real(8) pure function junction_get_L2norm (residQ, residEJ, np)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the L2 norm of residQ for multiple branches and the
        !% single EJ for one junction
        !%------------------------------------------------------------------
            real(8), intent(in) :: residQ(:), residEJ
            integer, intent(in) :: np
        !%------------------------------------------------------------------

        junction_get_L2norm = sqrt(sum(residQ(1:np)**twoI) + residEJ**twoI)

    end function junction_get_L2norm
!%
!%========================================================================== 
!%========================================================================== 
!%      
    real(8) pure function junction_branch_Qnet (jR, np, idir)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the net Q in (idir = 1) or net Q out (idir = -1)
        !% for a junction
        !%------------------------------------------------------------------
            real(8), intent(in) :: jR(:,:)
            integer, intent(in) :: np, idir
        !%------------------------------------------------------------------

        junction_branch_Qnet = onehalfR * sum                                      &
          (                                                                 &
            jR(1:np,jr_beta) * jR(1:np,jr_Q)                                &
            * (                                                             &
                oneR + real(idir,8) * jR(1:np,jr_beta) * jR(1:np,jr_Q)      & 
                       / abs(jR(1:np,jr_Q))                                 &   
              )                                                             &
          )

    end function junction_branch_Qnet
!%
!%========================================================================== 
!%========================================================================== 
!%  
    real(8) pure function junction_Qoverflow &
            (EJ, Zcrown, fcoef, OverflowHeight, OverflowType, &
             isFull, isSubmergedOrifice)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the overflow based on theory
        !%------------------------------------------------------------------
            real(8), intent(in) :: EJ, Zcrown, fcoef, OverflowHeight
            integer, intent(in) :: OverflowType
            logical, intent(in) :: isFull, isSubmergedOrifice
        !%------------------------------------------------------------------

        !% --- baseline is zero overflow
        junction_Qoverflow = zeroR
        if (.not. isFull) return 

        select case (OverflowType)
            case (NoOverflow)
                !% --- retain zero
            case (OverflowOrifice)
                if (isSubmergedOrifice) then 
                    junction_Qoverflow = - twothirdR * fcoef                        &
                        * (                                                &
                            ((EJ - Zcrown)**threehalfR)                    &
                           -((EJ - (Zcrown + OverflowHeight))**threehalfR) &
                        )

                else
                    junction_Qoverflow = - twothirdR * fcoef * ((EJ - Zcrown)**(threehalfR))
                end if
            case (OverflowWeir)
                junction_Qoverflow = -twoR * fcoef * ((EJ - Zcrown)**(threehalfR))
            case (Ponded)
                !% --- overflow volume is held in JM volume
                !%     retain zero
            case Default
                !% --- code should not be here
                !%     to have a failure for testing 
                !%     the function must be unpure    
        end select

    end function junction_Qoverflow
!%
!%========================================================================== 
!%========================================================================== 
!%   
    pure function junction_conservation (jR, residEJ, QinTotal, QoutTotal, np)
        !%------------------------------------------------------------------
        !% Description
        !% Proportionally adjust inflow/outflow Q to ensure flowrate
        !% and storage residual is zero to machine accuracy
        !% Applies to a single junction
        !%------------------------------------------------------------------
            real(8), intent(in)     :: jR(:,:), QInTotal, QoutTotal, residEJ
            integer, intent(in)     :: np

            real(8), dimension(np)  :: junction_conservation
        !%------------------------------------------------------------------

        junction_conservation = jR(1:np,jr_Q) &
            * (                                                                                  &
               oneR - onefourthR * residEJ                                                       &
               * (                                                                               &
                    (oneR + (jR(1:np,jr_beta) * jR(1:np,jr_Q)) / abs(jR(1:np,jr_Q)) ) / QinTotal  &
                   +                                                                             &
                    (oneR - (jR(1:np,jr_beta) * jR(1:np,jr_Q)) / abs(jR(1:np,jr_Q)) ) / QoutTotal &
                 )                                                                               &
               )
            
    end function junction_conservation
!%
!%========================================================================== 
!%========================================================================== 
!%   
    pure subroutine junction_push_face_flowrate &
            (faceQ, faceVup, faceVdn, faceAup, faceAdn, jR, jI,  np, JMidx)
        !%------------------------------------------------------------------
        !% Description
        !% Cycles through the solved branches (1:np) and restores their
        !% Q data to the faceR array 
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(inout) :: faceQ(:),   faceVup(:), faceVdn(:)
            real(8), intent(in)    :: faceAup(:), faceAdn(:)
            real(8), intent(in)    :: jR(:,:)
            integer, intent(in)    :: np, JMidx
            integer, target, intent(in) ::  jI(:,:)

            integer :: mm, ii
            !integer, pointer :: tM, kidx, fIdx
            integer :: kidx(np), fIdx(np)
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        !% --- cycle through the branches that were solved
        do mm=1,np
            
            !% --- branch kk in the 1:max_branch_per_node
            kidx(mm) = jI(mm,ji_kidx)
            
            !% --- handle upstream or downstream separately
            if (mod(kidx(mm)+1,twoI) .eq. zeroI) then 
                !% --- upstream branch face
                fIdx(mm) = elemI(JMidx+kidx(mm),ei_Mface_uL)
            else
                !% --- downstream branch face
                fIdx(mm) = elemI(JMidx+kidx(mm),ei_Mface_dL)
            end if
            !% --- set the flowrate
            faceQ(fIdx(mm))   = jR(mm,jr_Q)
            !% --- get consistent velocities (assumes area > epsilon)
            faceVup(fIdx(mm)) = faceQ(fIdx(mm)) / faceAup(fIdx(mm))
            faceVdn(fIdx(mm)) = faceQ(fIdx(mm)) / faceAdn(fIdx(mm))
        end do

    end subroutine junction_push_face_flowrate   
!%
!%========================================================================== 
!%========================================================================== 
!% 
    real(8) pure function junction_volume &
            (oldVolume, Qin, Qout, Qlat, Qoverflow, istep)
        !%------------------------------------------------------------------
        !% Description
        !% Updates junction volume from Volume_N0 for an RK step
        !%------------------------------------------------------------------
        !% Declarations:
            real(8), intent(in) :: oldVolume, Qin, Qout, Qlat, Qoverflow
            integer, intent(in) :: istep
        !%------------------------------------------------------------------

        junction_volume = oldVolume                    &
            + setting%Solver%crk2(istep) * setting%Time%Hydraulics%Dt  &
            * (Qin - Qout + Qlat + Qoverflow)

    end function junction_volume 
!%
!%========================================================================== 
!%========================================================================== 
!% 
    subroutine junction_check_conservation (thisJM, Qoverflow, Npack, istep)
        !%------------------------------------------------------------------
        !% Description:
        !% Checks flowrate and junction storage conservation
        !%------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisJM(:)
            integer, intent(in) :: Npack, istep
            real(8), intent(in) :: Qoverflow(:)
            integer, pointer :: tM, upF, dnF
            integer          :: mm, kk
            real(8)          :: thisCons
        !%------------------------------------------------------------------

        do mm=1,Npack
            tM => thisJM(mm)
            
            !% --- lateral inflow
            thisCons = elemR(tM,er_FlowrateLateral) + Qoverflow(mm)

            !% --- change in storage
            !%     assumes old head stored in er_Temp01
            select case (elemSI(tM,esi_JunctionMain_Type))
                case (ImpliedStorage)
                    !% --- no change to conservation by storage
                case (TabularStorage,FunctionalStorage)
                    thisCons = thisCons &
                         - (elemSR(tM,esr_Storage_Plane_Area)    &
                              / (  setting%Solver%crk2(istep)    &
                                  * setting%Time%Hydraulics%Dt ) &
                           )                                     &
                            * (elemR(tM,er_Head) - elemR(tM,er_Temp01))

                case default
                    print *, 'CODE ERROR: unexpected case default'
            end select

            !% --- upstream branches
            do kk=1,max_branch_per_node,2
                if (elemSI(tM+kk,esi_JunctionBranch_Exists) == oneI) then
                    upF => elemI(tM+kk,ei_Mface_uL)
                    thisCons = thisCons + faceR(upF,fr_Flowrate)
                end if
            end do

            !% --- downstream branches
            do kk=2,max_branch_per_node,2
                if (elemSI(tM+kk,esi_JunctionBranch_Exists) == oneI) then
                    dnF => elemI(tM+kk,ei_Mface_dL)
                    thisCons = thisCons + faceR(dnF,fr_Flowrate)
                end if
            end do

            if (abs(thisCons) > 1.d-6) then 
                print *, 'conservation error at ',thisJM(mm)
                print *, 'magnitude ',thisCons
                stop 5509873
            end if
        end do
        
    end subroutine junction_check_conservation
!%
!%========================================================================== 
!% JUNCTION PACKING BELOW
!%========================================================================== 
!%
    pure subroutine junction_check_repack_all (isRepack, jR, np, Npack, jEpsilon)
        !%------------------------------------------------------------------
        !% Description
        !% Checks if any branches are now zero and need repacking and 
        !% adjusts beta to zero where branches are zero
        !%------------------------------------------------------------------
        !% Declarations
            logical, intent(inout) :: isRepack(:)
            real(8), intent(inout) :: jR(:,:,:)
            integer, intent(in)    :: np(:), Npack
            real(8), intent(in)    :: jEpsilon

            integer :: mm
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        do concurrent (mm=1:Npack)

            call junction_check_repack_oneJ (jR(mm,:,:), np(mm), jEpsilon, isRepack(mm))

        end do

    end subroutine junction_check_repack_all   
!%
!%========================================================================== 
!%========================================================================== 
!%    
    pure subroutine junction_check_repack_oneJ (jR, np, jEpsilon, isRepack)
        !%------------------------------------------------------------------
        !% Checks for zeros and resets beta if repack needed
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(inout) :: jR(:,:)
            real(8), intent(in)    :: jEpsilon
            integer, intent(in)    :: np 
            logical, intent(inout) :: isRepack
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        if (any(jR(1:np,jr_LambdaC) < jEpsilon)) then 
            isRepack = .true.
            where (jR(1:np,jr_LambdaC) < jEpsilon)
                jR(1:np,jr_beta) = zeroR
            endwhere
        else
            isRepack = .false.
        end if

    end subroutine junction_check_repack_oneJ
!%
!%========================================================================== 
!%========================================================================== 
!%
    pure subroutine junction_pack_all &
        (jR, jI, np, isRepack, Npack)
        !%------------------------------------------------------------------
        !% Description
        !% Packing on the 2nd dimension of 3D arrays, jR, jI
        !% Cycles through the 1st dimension (packed JM index) in this procedure 
        !% and then cycles thorugh through the 3rd dimension (data type)
        !% in call to subsidiary procedures
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(inout) :: np(:), jI(:,:,:)
            real(8), intent(inout) :: jR(:,:,:)
            logical, intent(inout) :: isRepack(:)
            integer, intent(in)    :: Npack
        
            integer :: mm
        !%------------------------------------------------------------------
        !% --- pack the jDataR and jDataI to remove zero branches
        do concurrent (mm=1:Npack)
            call junction_pack_oneJ(jR(mm,:,:), jI(mm,:,:), np(mm), isRepack(mm))
        end do

    end subroutine junction_pack_all
!%
!%========================================================================== 
!%========================================================================== 
!%   
    pure subroutine junction_pack_oneJ (jR, jI, np, isRepack)
        !%------------------------------------------------------------------
        !% Description
        !% Packs data for one junction during implicit junction solution
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(inout) :: jR(:,:)
            integer, intent(inout) :: jI(:,:), np
            logical, intent(inout) :: isRepack
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        if (.not. isRepack) return 
    
        !% --- count the nonzero branches
        np = count(jR(:,jr_beta) .ne. zeroR)

        !% --- perform the packs based on jr_beta (must delay packing jr_beta)
        call junction_pack_all_jRtype (jR(:,:), np)
        call junction_pack_all_jItype (jI(:,:), jR(:,:), np)

        !% --- pack beta (always must be last)
        jR(1:np,jr_beta)  = junction_pack_one_arrayR &
                (jR(:,jr_beta), jR(:,jr_beta), np)

        !% --- repack finished, reset the logical
        isRepack = .false.

    end subroutine junction_pack_oneJ
!%
!%========================================================================== 
!%========================================================================== 
!%        
    pure subroutine junction_pack_all_jRtype (jR, np)  
        !%------------------------------------------------------------------
        !% Description:
        !% Packs the first dimension of a 2D real array (jR) by cycling
        !% through the data types in the 2nd dimension.
        !% Note that np must equal the number of jR(:,jr_beta) that are 
        !% nonzero
        !%------------------------------------------------------------------
        !% Declarations:
            real(8), intent(inout) :: jR(:,:)
            integer, intent(in)    :: np

            !% --- jData elements that are packed
            integer, parameter :: NjRset = 13
            integer :: jRset(NjRset)
            integer :: jj
        !%------------------------------------------------------------------

        !% --- define data sets used in packing 
        !%     DO NOT INCLUDE jr_beta -- must be repacked last in a separate 
        !%     call since it is themask used for packing.
        jRset = [jr_Area, jr_a, jr_c, jr_Ebranch, jr_Gamma, jr_K, jr_Length, jr_Q, &
                 jr_DeltaQ, jr_residQ, jr_LambdaA, jr_LambdaB, jr_LambdaC]

        do concurrent (jj=1:NjRset)
            !% --- pack the jRset(jj) data type for the thisJM(mm) branch
            !%     to remove non-existent branches
            jR(1:np,jRset(jj))  = junction_pack_one_arrayR                               &
                        (jR(:,jRset(jj)), jR(:,jr_beta), np)
        end do

    end subroutine junction_pack_all_jRtype
!%
!%========================================================================== 
!%========================================================================== 
!%        
    pure subroutine junction_pack_all_jItype (jI, jR, np)  
        !%------------------------------------------------------------------
        !% Description:
        !% Packs the first dimension of a 2D integer array (jI) by cycling
        !% through the data types in the 2nd dimension.
        !% Note that np must equal the number of jR(:,jr_beta) that are 
        !% nonzero
        !%------------------------------------------------------------------
            integer, intent(inout) :: jI(:,:)
            real(8), intent(in)    :: jR(:,:)
            integer, intent(in)    :: np

            integer, parameter :: NjIset = Ncol_jDataI
            integer            :: jIset(NjIset)
            integer            :: jj
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        !% --- define data sets used in packing
        jIset = [ji_kidx]

        do concurrent (jj=1:NjIset)
            !% --- pack the jIset(jj) data type for the thisJM(mm) branch
            !%     to remove non-existent branches
            jI(1:np,jIset(jj))  = junction_pack_one_arrayI                               &
                        (jI(:,jIset(jj)), jR(:,jr_beta), np)
        end do

    end subroutine junction_pack_all_jItype
!%
!%========================================================================== 
!%========================================================================== 
!%      
    pure function junction_pack_one_arrayR (inarray, beta, np)
        !%------------------------------------------------------------------
        !% Description
        !% standard pack to remove any locations where beta = 0.0
        !%------------------------------------------------------------------
            real(8), intent(in)    :: inarray(:), beta(:)
            integer, intent(in)    :: np
            real(8), dimension(np) :: junction_pack_one_arrayR
        !%------------------------------------------------------------------

        junction_pack_one_arrayR(1:np) = pack (inarray, beta .ne. zeroR)

    end function junction_pack_one_arrayR
!%
!%==========================================================================  
!%========================================================================== 
!%      
    pure function junction_pack_one_arrayI (inarray, beta, np)
        !%------------------------------------------------------------------
        !% Description
        !% standard pack to remove any locations where beta = 0.0
        !%------------------------------------------------------------------
            real(8), intent(in)    ::  beta(:)
            integer, intent(in)    :: np, inarray(:)
            integer, dimension(np) :: junction_pack_one_arrayI
        !%------------------------------------------------------------------

        junction_pack_one_arrayI(1:np) = pack (inarray, beta .ne. zeroR)

    end function junction_pack_one_arrayI
!%
!%==========================================================================  
    !%========================================================================== 
!%
    function matinv(A) result (invA)
        !%-----------------------------------------------------------------
        !% Description:
        !% Find the inverse of an input matrix A
        !%-----------------------------------------------------------------
        !% Declarations
        real(8),allocatable, intent(in) :: A(:,:)
        real(8),allocatable             :: invA(:,:)
        integer                         :: n   
        !%-----------------------------------------------------------------
        !% find the size of the matrix
        n = size(A,1)

        !% for square matrix size or 2-4, direct inversion is the fastest 
        if (n == 2) then
            invA = matinv2(A)
        else if (n == 3) then
            invA = matinv3(A)
        else if (n == 4) then
            invA = matinv4(A)
        else
        !% use the LAPACK for higher order matrices
        !% talk with cesar about LAPACK installation recipie
            print*, 'Matrix inversion for shape > 4 has not yet been developed'
            stop 789451
        end if     

    end function 
!%
!%========================================================================== 
!%========================================================================== 
!%
    pure function matinv2(A) result(B)
        !! Performs a direct calculation of the inverse of a 2×2 matrix.
        real(8), intent(in) :: A(2,2)   !! Matrix
        real(8)             :: B(2,2)   !! Inverse matrix
        real(8)             :: detinv

        ! Calculate the inverse determinant of the matrix
        detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

        ! Calculate the inverse of the matrix
        B(1,1) = +detinv * A(2,2)
        B(2,1) = -detinv * A(2,1)
        B(1,2) = -detinv * A(1,2)
        B(2,2) = +detinv * A(1,1)
    end function
!%
!%========================================================================== 
!%========================================================================== 
!%
    pure function matinv3(A) result(B)
        !! Performs a direct calculation of the inverse of a 3×3 matrix.
        real(8), intent(in) :: A(3,3)   !! Matrix
        real(8)             :: B(3,3)   !! Inverse matrix
        real(8)             :: detinv

        ! Calculate the inverse determinant of the matrix
        detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
                - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
                + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

        ! Calculate the inverse of the matrix
        B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
        B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
        B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
        B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
        B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
        B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
        B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
        B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
        B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    end function
!%
!%========================================================================== 
!%========================================================================== 
!%
    pure function matinv4(A) result(B)
        !! Performs a direct calculation of the inverse of a 4×4 matrix.
        real(8), intent(in) :: A(4,4)   !! Matrix
        real(8)             :: B(4,4)   !! Inverse matrix
        real(8)             :: detinv

        ! Calculate the inverse determinant of the matrix
        detinv = &
        1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
        - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
        + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
        - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

        ! Calculate the inverse of the matrix
        B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
        B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
        B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
        B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
        B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
        B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
        B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
        B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
        B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
        B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
        B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
        B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
        B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
        B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
        B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
        B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
    end function 
!%
!%========================================================================== 
!% END OF MODULE
!%==========================================================================
end module junction_elements