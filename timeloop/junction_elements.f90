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
    

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine  junction_calculation (whichTM, istep)
        !%-----------------------------------------------------------------
        !% Description:
        !% get the jacobean matrix for junction element following the 
        !% derivation by ...
        !%-----------------------------------------------------------------
        !% Declarations
        integer, intent(in) :: whichTM, istep
        integer, pointer :: thisColP_JM, thisP(:), BranchExists(:), tM, iup(:), idn(:)
        integer, pointer :: nBranches(:), iFaceUp(:), iFaceDn(:)
        integer, pointer :: Npack, tFup, tFdn
        real(8), pointer :: eFlow(:), eHead(:), dfdQ(:), fFlow(:), fHead_u(:), fHead_d(:)
        real(8), pointer :: planArea(:), dt, grav, epsH, crk(:)
        real(8), pointer :: tempFlow(:), tempHead(:)
        integer :: ii, jj, kk, tB, jDim, info
        real(8) :: alpha, beta, gamma, kFactor
        real(8), allocatable :: jMatrix(:,:), inv_jMatrix(:,:), fVector(:), fSolution(:)
        !%-----------------------------------------------------------------------------
        !%
        nBranches    => elemSI(:,esi_JunctionMain_Total_Branches)
        BranchExists => elemSI(:,esi_JunctionBranch_Exists)
        eFlow        => elemR(:,er_Flowrate)
        eHead        => elemR(:,er_Head)
        dfdQ         => elemSR(:,esr_JunctionBranch_dfdQ)
        planArea     => elemSR(:,esr_Storage_Plane_Area)
        tempFlow     => elemR(:,er_Temp01)
        tempHead     => elemR(:,er_Temp01)

        iFaceUp      => elemI(:,ei_Mface_uL)
        iFaceDn      => elemI(:,ei_Mface_dL)

        fFlow        => faceR(:,fr_Flowrate)
        fHead_u      => faceR(:,fr_Head_u)
        fHead_d      => faceR(:,fr_Head_d)

        crk          => setting%Solver%crk2
        dt           => setting%Time%Hydraulics%Dt
        grav         => setting%constant%gravity
        epsH         => setting%Eps%Head
        !%-----------------------------------------------------------------------------
        !%
        select case (whichTM)
        case (ALLtm)
            thisColP_JM            => col_elemP(ep_JM_ALLtm)
        case (ETM)
            thisColP_JM            => col_elemP(ep_JM_ETM)
        case (AC)
            thisColP_JM            => col_elemP(ep_JM_AC)
        case default
            print *, 'CODE ERROR: time march type unknown for # ', whichTM
            print *, 'which has key ',trim(reverseKey(whichTM))
            stop 7659
        end select


        Npack => npack_elemP(thisColP_JM)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisColP_JM)

            !% ================================================================================
            !% first loop through the junctions to branch flowrate heads
            !% and initial guess to the solution
            do ii=1,Npack
                tM => thisP(ii)  !% JM junction main ID

                tempHead(tM) = eHead(tM)
                !% --- Junction branch setup
                !% set the junction branch head and flowrate from the faces
                !% also set the first approximation from the previously solved value
                do kk=1,max_branch_per_node,2
                    tB = tM + kk  !% JB branch ID
                    if (BranchExists(tB)==1) then
                        !% pointer to the face upstream of the branch
                        tFup => iFaceUp(tB)
                        !% copy the face flow over to junction branch as the first guess
                        eFlow(tB) = fFlow(tFup) 
                        !% copy the downstream face head over to junction branch as the first guess
                        eHead(tB) = fHead_d(tFup) 
                        !% copy the head and the flow as a first guess
                        tempHead(tB) = eHead(tB)
                        tempFlow(tB) = eHead(tB) 
                    end if
                end do
                !% handle the downstream branches
                do kk=2,max_branch_per_node,2
                    !print *, kk ,'in junction branch'
                    tB = tM + kk
                    if (BranchExists(tB)==1) then
                        !% pointer to the face downstream of the branch
                        tFdn => iFaceDn(tB)
                        !% copy the face flow over to junction element as a first guess
                        eFlow(tB) = fFlow(tFdn)
                        !% copy the upstream face head over to junction branch as the first guess
                        eHead(tB) = fHead_u(tFdn)
                        !% copy the head and the flow as a first guess
                        tempHead(tB) = eHead(tB)
                        tempFlow(tB) = eHead(tB)
                    end if
                end do
            end do
            !% ================================================================================


            !% ================================================================================
            !% secondly loop through all the junctions to set up the jacobean matrices
            do ii=1,Npack
                tM => thisP(ii)  !% JM junction main ID

                !% --- jacobean matrix allocation
                !% find the dimension of the square jacobean matrix
                !% this will also be the dimension of the function vector
                jDim = nBranches(tM) + oneI
                !% allocate the jacobean martix for this junction main
                allocate(jMatrix(jDim, jDim))
                !% populate the jacobean matrix with zero
                jMatrix(:,:) =  zeroR
                !% allocate the inverse jacobean martix for this junction main
                allocate(inv_jMatrix(jDim, jDim))
                !% populate the jacobean matrix with zero
                inv_jMatrix(:,:) =  zeroR
                !% counter for jMatrix row index advance
                jj = zeroI

                !% --- function vector allocation
                !% allocate the function vector for this junction main
                allocate(fVector(jDim))
                !% populate the function vector with zero
                fVector(:) =  zeroR
                !% allocate the function solution for this junction main
                allocate(fSolution(jDim))
                !% populate the function solution with zero
                fSolution(:) =  zeroR

                !% --- Junction branch setup
                !% find the jacobean derivative
                do kk=1,max_branch_per_node,2
                    tB = tM + kk  !% JB branch ID
                    if (BranchExists(tB)==1) then
                        !% for a valid branch advance the row index of the jMatrix
                        jj = jj + oneI

                        !% assuming at this point the junction branch already had the flow and
                        !% heads copied over from the subsequent faces

                        !% --- HACK: hard coded for testing
                        !% --- gamma is one for upstream branches
                        gamma = oneR
                        !% --- beta is gamm muliplied by the flow direction
                        beta = gamma * sign(oneR, eHead(tB) - eHead(tM))
                        !% --- k is hardcoded to a constant value of 1 for u/s branch
                        kFactor = oneR
                        !% --- find the functional derivative of flow for the junction branch
                        dfdQ(tB) = twoR * beta * kFactor * abs(eFlow(tB))  

                        !% --- populate the jMatrix
                        !% the functional derivative will be at the diagonal
                        jMatrix(jj,jj) = dfdQ(tB)  
                        !% the last element of this row will be 1
                        jMatrix(jj,jDim) = oneR 
                        !% the last row of the jacobean matrix will beta for the consequent columns
                        jMatrix(jDim,jj) = beta   

                        !% --- populate the function vector
                        fVector(jj) = beta * kFactor * (eFlow(tB) ** twoI) + tempHead(tM) - eHead(tB)    

                        !% --- populate the solution vector with the first gussess
                        fSolution(jj) = tempFlow(tB)

                        !% the last element of the function vector is the sum of all the branch flowrates with 
                        !% some additional term, which will be dealt with later 
                        fVector(jDim) = fVector(jDim) + gamma * eFlow(tB)       
                    end if
                end do

                !% handle the downstream branches
                do kk=2,max_branch_per_node,2
                    !print *, kk ,'in junction branch'
                    tB = tM + kk
                    if (BranchExists(tB)==1) then
                        !% for a valid branch advance the row index of the jMatrix
                        jj = jj + oneI

                        !% assuming at this point the junction branch already had the flow and
                        !% heads copied over from the subsequent faces

                        !% --- HACK: hard coded for testing
                        !% --- gamma is negatice one for upstream branches
                        gamma = - oneR
                        !% --- beta is gamma muliplied by the flow direction
                        beta = gamma * sign(oneR, eHead(tM) - eHead(tB))
                        !% --- k is hardcoded to a constant value of 0.5 for u/s branch
                        kFactor = onehalfR
                        !% --- find the functional derivative of flow for the junction branch
                        dfdQ(tB) = twoR * beta * kFactor * abs(eFlow(tB))  

                        !% --- populate the jMatrix
                        !% the functional derivative will be at the diagonal
                        jMatrix(jj,jj) = dfdQ(tB)  
                        !% the last element of this row will be 1
                        jMatrix(jj,jDim) = oneR  
                        !% the last row of the jacobean matrix will beta for the consequent columns
                        jMatrix(jDim,jj) = beta  

                        !% --- populate the function vector
                        fVector(jj) = beta * kFactor * (eFlow(tB) ** twoI) + tempHead(tM) - eHead(tB)

                        !% --- populate the solution vector with the first gussess
                        fSolution(jj) = tempFlow(tB)

                        !% the last element of the function vector is the sum of all the branch flowrates with 
                        !% some additional term, which will be dealt with later 
                        fVector(jDim) = fVector(jDim) + gamma * eFlow(tB)                            
                    end if
                end do

                alpha =  planArea(tM) / crk(istep) * dt
                !% ---  after finishing looping through branches, populate the last elemet of the jacobean matrix
                jMatrix(jDim,jDim) = alpha

                !% --- after finishing looping through branches, add the alpha terms to the function vector
                fVector(jDim) = fVector(jDim) + alpha * (tempHead(tM) -  eHead(tM))

                !% --- after finishing looping through branches, add the head of the JM at the last element space
                fSolution(jDim) = tempHead(tM)

                ! !% find the inverse of the jacobean matrix
                inv_jMatrix = matinv(jMatrix)
    

                !% find the new solution
                fSolution = fSolution - matmul(inv_jMatrix,fVector)
            end do
            !% ================================================================================

        end if

    end subroutine junction_calculation
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