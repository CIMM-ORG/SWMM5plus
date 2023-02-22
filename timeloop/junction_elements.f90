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
    public :: junction_force_Qinterpweights

    real(8), parameter :: Cbc = 0.46295d0  !% HACK dimensionless 1.45/sqrt(g) from Brater and King broad-crested weir coefficient
    real(8), parameter :: Horifice = 0.15d0  !% HACK fixed orifice height
    real(8), parameter :: Lorifice = 1.5d0   !% HACK fixed orifice length

    real(8) :: coef1, coef2, coef3, coef4

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

            !% --- coefficients in orifice and weir eqquations
            coef1 = twoR * Cbc * sqrt(setting%Constant%gravity * setting%Constant%pi)
            coef2 = threehalfR * coef1
            coef3 = twothirdR  * sqrt(twoR * setting%Constant%gravity)
            coef4 = threehalfR * coef3
        !%-----------------------------------------------------------------

        call junction_calculation (Npack, thisColP, istep)

        ! print *, ' '
        ! print *, 'after junction calc ', elemR(4,er_Volume), elemR(4,er_Head)
        ! print *, ' '

    end subroutine junction_toplevel
!%
!%==========================================================================
!%==========================================================================
!% 
    subroutine junction_force_Qinterpweights (thisCol)
        !%------------------------------------------------------------------
        !% Description:
        !% sets the JB interpweights to maximum so that adjacent
        !% values are interpolated to face
        !%------------------------------------------------------------------
            integer, intent(in) :: thisCol
            integer, pointer    :: npack, thisP(:)
            integer             :: ii
        !%------------------------------------------------------------------
        !% Aliases
            npack => npack_elemP(thisCol)
            if (npack < 1) return
            thisP => elemP(1:npack,thisCol)
        !%------------------------------------------------------------------

        do ii=1,max_branch_per_node  
            elemR(thisP+ii,er_InterpWeight_uQ) = setting%Limiter%InterpWeight%Maximum
            elemR(thisP+ii,er_InterpWeight_dQ) = setting%Limiter%InterpWeight%Maximum
        end do
        
    end subroutine junction_force_Qinterpweights
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

            real(8) :: QnetBranches, dQdHsum, dQdHstorage, Qoverflow, dQdHoverflow
            real(8) :: dH, resid, Qstorage, QnetIn, QnetOut

            integer :: mm, kk

            real(8), parameter :: localEpsilon = 1.0d-6

        !%-----------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%-----------------------------------------------------------------   
        !% Aliases
            thisJM      => elemP(1:Npack,thisColP)
        !%-----------------------------------------------------------------  
        
        do mm=1,Npack
            !% --- store face flowrate in JB for upstream (1) and downstream (2)
            call junction_branch_getface (elemR(:,er_Flowrate),fr_Flowrate,thisJM(mm),ei_Mface_uL,1)
            call junction_branch_getface (elemR(:,er_Flowrate),fr_Flowrate,thisJM(mm),ei_Mface_dL,2)

            ! print *, ' '
            ! print *, '======================================================='
            ! print *, 'flowrates JM =', thisJM(mm)
            ! print *, elemR(thisJM(mm)+1,er_Flowrate), elemR(thisJM(mm)+2,er_Flowrate), elemR(thisJM(mm),er_FlowrateLateral)

            !% --- store face energy head in JB for upstream (1) and downstream (2)
            call junction_branch_getface (elemR(:,er_EnergyHead),fr_EnergyHead,thisJM(mm),ei_Mface_uL,1)
            call junction_branch_getface (elemR(:,er_EnergyHead),fr_EnergyHead,thisJM(mm),ei_Mface_dL,2)

            ! print *, ' '
            ! print *, 'energyHead'
            ! print *, elemR(thisJM+1,er_EnergyHead), elemR(thisJM+2,er_EnergyHead)

            !% --- compute net flowrate from branches
            QnetBranches = junction_main_QnetBranches (thisJM(mm))

            ! print *, ' '
            ! print *, 'Qnetbranches ',QnetBranches

            !% --- compute overflow rate
            Qoverflow = junction_main_Qoverflow (thisJM(mm))

            !% --- compute storage rate of change with head
            dQdHstorage = junction_main_dQdHstorage (thisJM(mm),iStep)

            !% --- compute overflow rate with change in head
            dQdHoverflow = junction_main_dQdHoverflow (thisJM(mm))

            ! print *, ' '
            ! print *, 'dQdHoverflow ',dQdHoverflow

            !% --- compute dQdH for each branch -- must include diagnostic
            call junction_branch_dQdH (elemR(:,er_dQdH),thisJM(mm), elemI(:,ei_Mface_uL), 1)
            call junction_branch_dQdH (elemR(:,er_dQdH),thisJM(mm), elemI(:,ei_Mface_dL), 2)

            ! print *, ' '
            ! print *, 'dQdH branches'
            ! print *, elemR(thisJM+1,er_dQdH), elemR(thisJM+2,er_dQdH)

            !% --- compute net branch dQdH (product with beta)
            dQdHsum =  junction_main_branchdQdHsum (thisJM(mm))

            ! print *, ' '
            ! print *, 'dQdHsum ',dQdHsum

            !% --- compute the junction head change
            dH = junction_main_dH &
                (thisJM(mm), Qnetbranches, Qoverflow, dQdHstorage, dQdHoverflow, dQdHsum)

            ! print *, ' '
            ! print *, 'dH before limit = ',dH   
            ! print *, ' '

            ! print *, 'max gain ',junction_dH_maxgain (elemR(:,er_EnergyHead), elemR(thisJM(mm),er_Head), thisJM(mm))
            ! print *, 'max loss ',junction_dH_maxloss (elemR(:,er_EnergyHead), elemR(thisJM(mm),er_Head), thisJM(mm))

            !% limit dH to prevent flow direction change
            if (dH > zeroR) then
                dH = min(dH, junction_dH_maxgain (elemR(:,er_EnergyHead), elemR(thisJM(mm),er_Head), thisJM(mm)))
            elseif (dH < zeroR) then
                dH = max(dH, junction_dH_maxloss (elemR(:,er_EnergyHead), elemR(thisJM(mm),er_Head), thisJM(mm)))
            else 
                !% -- if dH = zero, no change
            end if

            ! print *, ' '
            ! print *, 'dH after first limit = ',dH   
            ! print *, ' '

            !% limit dH to prevent negative head
            if ((elemR(thisJM(mm),er_Head) - elemR(thisJM(mm),er_Zbottom) + dH) < setting%ZeroValue%Depth) then 
                ! print *, 'limiters '
                ! print *, elemR(thisJM(mm),er_Head),  (elemR(thisJM(mm),er_Zbottom) + 0.99d0*setting%ZeroValue%Depth)
                dH = (elemR(thisJM(mm),er_Zbottom) + 0.99d0*setting%ZeroValue%Depth) - elemR(thisJM(mm),er_Head)
            end if

            ! print *, ' '
            ! print *, 'dH after second limit = ',dH   
            ! print *, ' '

            !% --- limit dH dropping based on where overflow shuts off
            if (Qoverflow .ne. zeroR) then 
                dH = max(dH, junction_dH_overflow_min(thisJM(mm)))
            end if

            ! print *, ' '
            ! print *, 'dH after third limit = ',dH   
            ! ! print *, 'available ', elemR(thisJM(mm),er_Head) - elemR(thisJM(mm),er_Zbottom)
            ! ! print *, ' '

            !% --- update junction head
            elemR(thisJM(mm),er_Head) = elemR(thisJM(mm),er_Head) + dH

            !% --- update branch Q
            elemR(thisJM(mm)+1:thisJM(mm)+max_branch_per_node,er_Flowrate) &
                = elemR(thisJM(mm)+1:thisJM(mm)+max_branch_per_node,er_Flowrate)  &
                + dH * elemR(thisJM(mm)+1:thisJM(mm)+max_branch_per_node,er_dQdH)


                ! print *, ' '
                ! print *, 'flowrates after update'
                ! !do kk=1,max_branch_per_node
                ! do kk=1,3
                !     print *, 'Q : ',kk, elemR(thisJM(mm)+kk,er_Flowrate)
                ! end do

            !% --- update overflow
            Qoverflow = Qoverflow + dH * dQdHoverflow

            !% --- update storage flow rate
            select case (elemSI(thisJM(mm),esi_JunctionMain_Type))
                case (ImpliedStorage)
                    Qstorage = zeroR
                case (TabularStorage,FunctionalStorage)
                    !% --- compute the new storage volume
                    Qstorage = elemSR(thisJM(mm),esr_Storage_Plan_Area) * dH 
                    elemR(thisJM(mm),er_Volume) = elemR(thisJM(mm),er_Volume_N0)+ Qstorage

                    !% --- compute the storage flowrate
                    Qstorage = Qstorage / (setting%Solver%crk2(istep) * setting%Time%Hydraulics%Dt )
                case default
                    print *, 'CODE ERROR: unexpected case default'
                    call util_crashpoint(882873)
            end select

            !% --- adjust for non-conservation
            resid = junction_conservation_residual (thisJM(mm), Qoverflow, Qstorage) 

                ! print *, 'conservation resid ',resid
                ! print *, ' '

            if (abs(resid) > 1.0d-16) then 
                QnetIn  = junction_branch_Qnet (thisJM(mm),+oneI)
                QnetOut = junction_branch_Qnet (thisJM(mm),-oneI) + Qoverflow

                ! print *, 'starting conservation resid flows '
                ! print *, QnetIn, QnetOut, elemR(thisJM(mm),er_FlowrateLateral)

                if (elemR(thisJM(mm),er_FlowrateLateral) .ge. zeroR) then 
                    QnetIn = QnetIn + elemR(thisJM(mm),er_FlowrateLateral)
                else
                    QnetOut = QnetOut + elemR(thisJM(mm),er_FlowrateLateral)
                end if


                call junction_fix_conservation(thisJM(mm),resid, Qoverflow, Qstorage, QnetIn, QnetOut)
            end if    

            !% --- update volume overflow 
            elemR(thisJM(mm),er_VolumeOverflow) = Qoverflow * setting%Time%Hydraulics%Dt

            ! print *, ' '
            ! print *, 'flowrates after conservation force', thisJM(mm)
            ! !do kk=1,max_branch_per_node
            ! do kk=1,3
            !     print *, 'Q : ',kk, elemR(thisJM(mm)+kk,er_Flowrate)
            ! end do

            ! print *, 'final conservation resid ',resid
            ! print *, ' '



            !% --- push junction JB flowrate values back to face  
            call junction_branchface_forceJBvalue (fr_Flowrate, er_Flowrate, ei_Mface_uL, thisJM(mm), 1) 
            call junction_branchface_forceJBvalue (fr_Flowrate, er_Flowrate, ei_Mface_dL, thisJM(mm), 2) 

            ! print *, ' '
            ! print *, 'faces flowrate should be JB values'
            ! print *, elemR(5,er_Flowrate), elemR(6,er_Flowrate)
            ! print *, faceR(4,fr_Flowrate), faceR(5,fr_Flowrate)

            ! if (abs(elemR(5,er_Flowrate) - elemR(6,er_Flowrate)) > localEpsilon) then 
            !     print *, 'flowrate mismatch '
            !     stop 669874
            ! end if 

            !% --- note the head values are handled in update_auxiliary_variables_JM


            ! print *, ' '
            ! print *, 'dH after limit = ',dH   
            ! print *, ' '

            ! if (abs(elemR(thisJM(mm),er_Head)) > 2000.d0) then 
            !     print *, 'Strange Head at ',thisJM(mm)
            !     print *, elemR(thisJM(mm),er_Head)
            !     call util_crashpoint(709874)
            ! end if
        end do

    end subroutine junction_calculation
!%    
!%==========================================================================
!%==========================================================================
!% 
    pure subroutine junction_branch_getface (outdata, frCol, JMidx, fiIdx, kstart)
        !%-----------------------------------------------------------------
        !% Description
        !% Stores face data on JB element
        !%-----------------------------------------------------------------
            real(8), intent(inout) :: outdata(:)  !% element data for output
            integer, intent(in)    :: JMidx       !% index of JM junction
            integer, intent(in)    :: fiIdx       !%  index of map up or down to face
            integer, intent(in)    :: frCol    !%  column in faceR array
            integer, intent(in)    :: kstart       !% = 1 for upstream branches, 2 for down 
            real(8) :: k1,k2   
        !%-----------------------------------------------------------------
    
        k1 = JMidx + kstart
        k2 = JMidx + max_branch_per_node

        !do concurrent (kk=kstart:max_branch_per_node:2)
            outdata(k1:k2:2) = faceR(elemI(k1:k2:2,fiIdx),frCol)            &
                             * real(elemSI(k1:k2:2,esi_JunctionBranch_Exists),8)
        !end do

        !print *, 'in get face', outdata(k1:k2:2)

    end subroutine junction_branch_getface
!%
!%==========================================================================
!%==========================================================================
!% 
    real(8) pure function junction_main_QnetBranches (JMidx)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the net flow rate into a junction from all the branches
        !% Note -- assumes that all branches that do not exist have zero 
        !% flowrate
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx
            integer :: k1, k2
        !%------------------------------------------------------------------
        k1 = JMidx + 1
        k2 = JMidx + max_branch_per_node

        ! print *, ' '
        ! print *, 'junction_main_QnetBranches ', JMidx
        ! print *, branchsign
        ! print *, ' '
        ! print *, elemR(k1:k2,er_Flowrate)


        junction_main_QnetBranches = sum(branchsign * elemR(k1:k2,er_Flowrate)) 

    end function junction_main_QnetBranches   
!%    
!%==========================================================================
!%==========================================================================
!% 
    real(8) function junction_main_Qoverflow (JMidx)  
        !%------------------------------------------------------------------
        !% Description
        !% Computes the overflow rate of a junction
        !% HACK: for an overflow orifice we have hard-coded the orifice
        !%   height and length. Later these should be user inputs
        !% HACK: for an overflow weir we have hard-coded the weir coefficient.
        !%   Later this should be moved to the settings structure
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx
        !%------------------------------------------------------------------  

        !% --- return zero if the head is below the crown at the start.
        if (elemR(JMidx,er_Head) .le. elemR(JMidx,er_Zcrown)) then 
            junction_main_Qoverflow = zeroR
            return
        end if

        select case (elemSI(JMidx,esi_JunctionMain_OverflowType))
            case (NoOverflow,Ponded)
                !% --- ponding an NoOverflow do not have separate volume accounting
                !%     for overflow
                junction_main_Qoverflow = zeroR
            case (OverflowWeir)
                !% --- weir overflow based on estimated circumference of storage plan area
                junction_main_Qoverflow = -coef1 * sqrt(elemSR(JMidx,esr_Storage_Plan_Area)) &
                    * ((elemR(JMidx,er_Head) - elemR(JMidx,er_Zcrown))**threehalfR)
            case (OverflowOrifice)
                !% --- orifice overflow assuming a single orifice of standard dimensions
                if (elemR(JMidx,er_Head) .le. (elemR(JMidx,er_Zcrown) + Horifice )) then
                    !% --- water surface below upper edge of orifice
                    junction_main_Qoverflow = -coef3 * Lorifice &
                        * (elemR(JMidx,er_Head) - elemR(JMidx,er_Zcrown))
                else
                    !% --- orifice is pressurized (head above the Zcrown + Horifice)
                    junction_main_Qoverflow = -coef3 * Lorifice  &
                        * (                                                                                 &
                            +  ((elemR(JMidx,er_Head) -  elemR(JMidx,er_Zcrown)            )**threehalfR)   &
                            -  ((elemR(JMidx,er_Head) - (elemR(JMidx,er_Zcrown) + Horifice))**threehalfR)   &
                          )
                end if
            case default
                !% --- should not reach here.
                !%     for debugging, change to impure function and uncomment
                !%     the following
                print *, 'CODE ERROR'
                print *, 'unexpected case default'
                call util_crashpoint(6209874)
        end select

    end function junction_main_Qoverflow   
!%    
!%==========================================================================
!%==========================================================================
!% 
    real(8) pure function junction_main_dQdHstorage (JMidx,iStep)  
        !%------------------------------------------------------------------
        !% Description
        !% Computes the storage rate of junction with tabular or functional
        !% storage
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx, iStep
        !%------------------------------------------------------------------  

        if (elemSI(JMidx,esi_JunctionMain_Type) .eq. ImpliedStorage) then
            !% --- no storage rate for implied storage junctions
            junction_main_dQdHstorage = zeroR
            return
        else
            !% --- Storage rate term based on time step
            junction_main_dQdHstorage = elemSR(JMidx,esr_Storage_Plan_Area) &
                / (setting%Solver%crk2(iStep) * setting%Time%Hydraulics%Dt)
        end if

    end function junction_main_dQdHstorage   
    !%    
!%==========================================================================    
!%==========================================================================
!% 
    real(8) pure function junction_main_dQdHoverflow (JMidx)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the overflow rate for  junction
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx
        !%------------------------------------------------------------------  

        ! print *, ' '
        ! print *, 'OverFlowType ',elemSI(JMidx,esi_JunctionMain_OverflowType)
        ! print *, NoOverflow, Ponded, OverflowWeir, OverflowOrifice
        ! print *, elemR(JMidx,er_Head) - elemR(JMidx,er_Zcrown)

        !% --- if not an overflow at this time
        if (elemR(JMidx,er_Head) .le. elemR(JMidx,er_Zcrown)) then 
            junction_main_dQdHoverflow = zeroR
            return 
        end if

        !% --- if possible overflow condition exists
        select case (elemSI(JMidx,esi_JunctionMain_OverflowType))  
            case (NoOverflow,Ponded)
                junction_main_dQdHoverflow = zeroR
            case (OverflowWeir)
                junction_main_dQdHoverflow = -coef2                            &
                    * sqrt(   ( elemSR(JMidx,esr_Storage_Plan_Area))           &
                            * ( elemR(JMidx,er_Head - elemR(JMidx,er_Zcrown))) &
                          )
            case (OverflowOrifice)
                junction_main_dQdHoverflow = -coef4 * Lorifice                  &
                        * sqrt(elemR(JMidx,er_Head) - elemR(JMidx,er_Zcrown))     
            case default
                !% --- should not reach here.
                !%     for debugging, change to impure function and uncomment
                !%     the following
                ! print *, 'CODE ERROR'
                ! print *, 'unexpected case default'
                ! call util_crashpoint(1209874)
        end select 

    end function junction_main_dQdHoverflow
!%    
!%==========================================================================    
!%==========================================================================
!% 
    subroutine junction_branch_dQdH (dQdH, JMidx, fm, kstart)    
        !%-----------------------------------------------------------------
        !% Description
        !% Computes dQdH for junction branches
        !%-----------------------------------------------------------------
            real(8), intent(inout) :: dQdH(:)
            integer, intent(in)    :: JMidx !% index of JM junction
            integer, intent(in)    :: fm(:) !% map up or down to face
            integer, intent(in)    :: kstart !% = 1 for upstream branches, 2 for down

            real(8), pointer :: fpsiL(:), fEnergyHead(:), eHead(:), fQ(:)
            integer ::  kk, tB
            real(8) ::  oneL = 1.d0

            real(8), parameter :: localEpsilon = 1.0d-6
        !%-----------------------------------------------------------------
        !% Aliases
            eHead       => elemR(:,er_Head)
            fEnergyHead => faceR(:,fr_EnergyHead)
            fpsiL       => faceR(:,fr_2B_psiL)
            fQ          => faceR(:,fr_Flowrate)
        !%-----------------------------------------------------------------

        !% HACK -- this needs to be cleaned up and made concurrent
        !%  and the subroutine should be pure.
        
        do kk=kstart,max_branch_per_node,2
            !% --- this branch index
            tB = JMidx + kk
            dQdH(tB) = zeroR !% ensure that even dummy branches have a value

            ! print *, ' '
            ! print *, '================================================='
            ! print *, 'kk here ',kk

            ! print *, 'branch exists ', elemSI(tB,esi_JunctionBranch_Exists)

            !% --- cycle if not a valid branch
            if (elemSI(tB,esi_JunctionBranch_Exists) .ne. oneI) cycle

            !% --- for adjacent channel/conduit branches
            if (elemSI(tB,esi_JunctionBranch_CC_adjacent) .eq. oneI) then

                ! print *, ' '
                ! print *, 'dif 1',fEnergyHead(fm(tB)) - eHead(JMidx)
                ! print *, 'fpsil',fpsiL(fm(tB))

                !% --- limiter for small values
                ! if ( (abs(fEnergyHead(fm(tB)) - eHead(JMidx)) < setting%ZeroValue%Depth) .or. &
                !      (abs(fpsiL(fm(tB))) < localEpsilon) ) then 
                !     dQdH(tB) = zeroR
                ! end if


                ! print *, ' '
                ! print *, 'tB', tB
                ! print *, 'fpsil   ',fpsiL(fm(tB))
                ! print *, 'head dif',(fEnergyHead(fm(tB)) - eHead(JMidx))

                !% --- first computation for dQdH using face 2 * beta * psi * L
                !%     along with face E and junction H
                ! dQdH(tB) = oneR / ( fpsiL(fm(tB)) * (fEnergyHead(fm(tB)) - eHead(JMidx)) ) 

                ! print *, ' '
                ! print *, 'tB ',tB
                ! print *, 'fpsiL    ',fpsiL(fm(tB))
                ! print *, 'head dif ',(fEnergyHead(fm(tB)) - eHead(JMidx))

                !% --- error checking -- first computation should be positive
                !%     if not, use the Chezy-Manning approach
                ! if (dQdH(tB) .le. zeroR) then 

                ! if (elemR(tB,er_FroudeNumber) .ge. 0.99d0) then 
                !     dQdH(tB) = zeroR

                ! else

                    ! print *, ' '
                    ! print *, 'tB ,kk ',tB, kk
                    ! print *, 'branchsign',branchsign(kk)
                    ! print *, 'area, hyrad          ',elemR(tB,er_Area), elemR(tB,er_HydRadius)
                    ! print *, 'MN, Lenght           ', elemR(tB,er_ManningsN) , elemR(tB,er_Length)
                    ! print *, 'delta E              ',fEnergyHead(fm(tB)) - eHead(JMidx)

                    !% --- Use CM approach
                    if (abs(fEnergyHead(fm(tB)) - eHead(JMidx)) > localEpsilon) then
                        dQdH(tB) = - branchsign(kk)                                                 &
                            * elemR(tB,er_Area) * (elemR(tB,er_HydRadius)**twothirdR)               &
                            / (                                                                     &
                                twoR * elemR(tB,er_ManningsN) * (sqrt(elemR(tB,er_Length)))         &
                                * sqrt(abs(fEnergyHead(fm(tB)) - eHead(JMidx)))  &
                                )

                        ! print *, 'dQdH                  ', dQdH(tB)
                    else 
                        dQdH(tB) = zeroR
                    end if
                ! end if
                    ! print *, 'CODE ERROR'
                    ! print *, 'unexpected negative value for dQdH first part'
                    ! call util_crashpoint(629784)
                ! else 
                !     !% --- second computation for dQdH
                !     dQdH(tB) = -branchsign(kk) * sign(sqrt(dQdH(tB)),fQ(fm(tB)))
                ! end if

                    

                    
                    ! !% --- zero all inflows
                    ! if ((elemR(tb,er_Flowrate) * branchsign(kk)) > zeroR) then 
                    !     dQdH(tB) = zeroR
                    ! end if

                    ! print *, ' tb, dQdH ',tB, dQdH(tB)

        

            else
                !% --- for adjacent diagnostic branches
                print *, 'code error - unfinished - diagnostic branch next to junction'
                call util_crashpoint(229874)
            end if

        end do

    end subroutine junction_branch_dQdH
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) pure function junction_main_branchdQdHsum (JMidx)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the net dQdH terms junction from all the branches
        !% Note -- assumes that all branches that do not exist have zero 
        !% dQdH
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMidx
            integer :: k1, k2
        !%------------------------------------------------------------------
        k1 = JMidx + 1
        k2 = JMidx + max_branch_per_node

        junction_main_branchdQdHsum = sum(branchsign * elemR(k1:k2,er_dQdH)) 

    end function junction_main_branchdQdHsum
!%    
!%==========================================================================
!%==========================================================================
!% 
    real(8) pure function junction_main_dH &
        (JMidx, QnetBranches, Qoverflow, dQdHstorage, dQdHoverflow, dQdHsum)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the change in head for the junction
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: JMidx
            real(8), intent(in) :: Qoverflow, QnetBranches
            real(8), intent(in) :: dQdHoverflow, dQdHsum, dQdHstorage
            real(8), pointer :: Qlat
            real(8) :: denominator

            real(8), parameter :: localEpsilon = 1.0d-16
        !%------------------------------------------------------------------

        denominator = dQdHstorage - dQdHoverflow - dQdHsum

        ! print *, ' '
        ! print *, ' flowrates ',elemR(JMidx,er_FlowrateLateral)
        ! print *, Qoverflow, QnetBranches
        ! print *, ' '
        ! print *, dQdHstorage, dQdHoverflow, dQdHsum
        ! print *, 'denominator ', denominator
        ! print *, ' '

        if (abs(denominator) > localEpsilon) then
            junction_main_dH = (elemR(JMidx,er_FlowrateLateral) + Qoverflow + QnetBranches) & 
                         / denominator
        else
            junction_main_dH = zeroR 
        end if

    end function junction_main_dH
!%    
!%==========================================================================
!%==========================================================================
!% 
    real(8) pure function junction_dH_maxgain (energyB, headJ, JMidx)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the maximum head change allowed to prevent reversing
        !% flow in a branch based on the smallest (non-zero) head difference
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: JMidx
            real(8), intent(in) :: energyB(:), headJ
            integer :: k1, k2
            real(8) :: dE(max_branch_per_node)
        !%------------------------------------------------------------------
        !% 
        !%------------------------------------------------------------------

        k1 = JMidx + 1
        k2 = JMidx + max_branch_per_node

        !% --- note the + abs() provides only E >= H branches.
        dE = onehalfR * (energyB(k1:k2) - headJ + abs(energyB(k1:k2) - headJ))  &
                         * real(elemSI(k1:k2,esi_JunctionBranch_Exists),8)  
        
        where (dE .le. setting%ZeroValue%Depth)
            dE = abs(nullvalueR)
        endwhere

        junction_dH_maxgain = minval(dE)
                    
    end function junction_dH_maxgain  
!%    
!%==========================================================================
!%==========================================================================
!% 
    real(8) pure function junction_dH_maxloss (energyB, headJ, JMidx)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the largest negative head change allowed to prevent reversing
        !% flow in a branch based on the smallest (non-zero) head difference
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: JMidx
            real(8), intent(in) :: energyB(:), headJ
            integer :: k1, k2
            real(8) :: dE(max_branch_per_node)
        !%------------------------------------------------------------------
        !% 
        !%------------------------------------------------------------------

        k1 = JMidx + 1
        k2 = JMidx + max_branch_per_node

        !% --- note the - abs() provides only E <= H branches.
        dE = onehalfR * (energyB(k1:k2) - headJ - abs(energyB(k1:k2) - headJ))  &
                         * real(elemSI(k1:k2,esi_JunctionBranch_Exists),8)  

        ! print *, ' '
        ! print *, 'in junction dh maxloss'
        ! print *, energyB(k1), headJ
        ! print *, energyB(k1) - headJ, abs(energyB(k1) - headJ)
        ! print *, dE(1),dE(2)
        ! print *, ' '
        
        where (dE .ge. -setting%ZeroValue%Depth)
            dE = -abs(nullvalueR)
        endwhere


        junction_dH_maxloss = maxval(dE)
                    
    end function junction_dH_maxloss
    !%    
!%==========================================================================
!%==========================================================================
!% 
    real(8) pure function junction_dH_overflow_min(JMidx)
        !%------------------------------------------------------------------
        !% Description:
        !% Limit the droppin of dH during overflow so that at most it goes
        !% an epsilon distance below the crown
        !%------------------------------------------------------------------
            integer, intent(in) :: JMidx

            real(8), parameter :: localEpsilon = 1.0d-4
        !%------------------------------------------------------------------

        junction_dH_overflow_min = &
             -(elemR(JMidx,er_Head) - elemR(JMidx,er_Zcrown)) - localEpsilon

    end function junction_dH_overflow_min
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) pure function junction_conservation_residual (JMidx, Qoverflow, Qstorage) 
        !%------------------------------------------------------------------
        !% Description:
        !% computes flowrate residual, which is >0 for too much inflow 
        !% and < 0 for too much outflow
        !%------------------------------------------------------------------
            integer, intent(in) :: JMidx
            real(8), intent(in) :: Qoverflow, Qstorage
            real(8) :: QnetBranches
        !%------------------------------------------------------------------

        QnetBranches = junction_main_QnetBranches (JMidx)

        ! print *, ' '
        ! print *, 'junction cons ', JMidx
        ! print *, QnetBranches, Qoverflow, Qstorage
        
        junction_conservation_residual =  &
            QnetBranches + Qoverflow + Qstorage + elemR(JMidx,er_FlowrateLateral)

    end function junction_conservation_residual
!%    
!%==========================================================================
!%==========================================================================
!% 
    real(8) pure function junction_branch_Qnet (JMidx,idir)    
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the net inflow (idir = 1) or outflow (idir=-1) from
        !% junction JMidx
        !%------------------------------------------------------------------
            integer, intent(in) :: JMidx,idir
            real(8), dimension(max_branch_per_node) :: Qdir
            integer :: kk
        !%------------------------------------------------------------------

        Qdir = zeroR
        do concurrent (kk=1:max_branch_per_node)
            if ((abs(elemR(JMidx+kk,er_Flowrate)) > zeroR) .and. &
                (elemSI(JMidx+kk,esi_JunctionBranch_exists)== oneI)) then

                Qdir(kk) =  branchsign(kk) *  elemR(JMidx+kk,er_Flowrate)                    &
                    *onehalfR * (                                                                       &
                                    oneR + real(idir,8) * branchsign(kk) * elemR(JMidx+kk,er_Flowrate)   &
                                    / abs(elemR(JMidx+kk,er_Flowrate))                                   &
                                )
                ! print *, 'kk ',kk
                ! print *,  'A ',branchsign(kk) *  elemR(JMidx+kk,er_Flowrate)
                ! print *,  'B ',oneR + real(idir,8) * branchsign(kk) * elemR(JMidx+kk,er_Flowrate)   &
                ! / abs(elemR(JMidx+kk,er_Flowrate)) 
            end if
        end do

        junction_branch_Qnet = sum(Qdir)
        

    end function junction_branch_Qnet
!%    
!%==========================================================================
!%==========================================================================
!%    
    subroutine junction_fix_conservation (JMidx, resid, Qoverflow, Qstorage, QnetIn, QnetOut)
        !%------------------------------------------------------------------
        !% Description:
        !% Ad hoc fix to junction flow rates to ensure mass conservation
        !% Typically required because dH adjustment was limited
        !%------------------------------------------------------------------
            integer, intent(in) :: JMidx
            real(8), intent(inout) :: Qoverflow, resid
            real(8), intent(in) :: QnetIn, QnetOut, Qstorage
            integer :: kk
            real(8) :: QratioIn, QratioOut, dQoverflow
            real(8), dimension(max_branch_per_node) :: dQ
            real(8), parameter :: localEpsilon = 1.0d-6
        !%------------------------------------------------------------------

        ! print *, ' '
        ! print *, 'JUNCTION ',JMidx
        ! print *, 'starting resid ',resid 

        dQ = zeroR
        dQoverflow = zeroR

        !% --- check for a degenerate condition
        if ((QnetIn .le. zeroR) .and. (QnetOut .ge. zeroR)) then 
            print *, 'CODE ERROR: unexpected zero fluxes and non-zero residual'
            call util_crashpoint(229873)
            return
        end if

        !% --- get the changes to flow rates
        if ((QnetOut .ge. zeroR) .and. (QnetIn > zeroR)) then 
            !% --- only inflows (no outflow)
            do kk=1,max_branch_per_node
                !% --- select inflows to adjust
                !%     signs based on QnetIn > 0
                if ((branchsign(kk) * elemR(JMidx+kk,er_Flowrate)) > zeroR) then
                    dQ(kk) = - resid * elemR(JMidx+kk,er_Flowrate) / QnetIn
                end if
            end do

        elseif ((QnetOut < zeroR) .and. (QnetIn .le. zeroR)) then 
            !% --- only outflows (no inflow)
            !%     note that signs based on QnetOut < 0
            do kk=1,max_branch_per_node
                !% -- select outflows to adjust
                if ((branchsign(kk) * elemR(JMidx+kk,er_Flowrate)) < zeroR) then
                    dQ(kk) = - resid * elemR(JMidx+kk,er_Flowrate) / QnetOut
                end if
            end do
            !% --- adjust overflow
            dQoverflow = - resid * Qoverflow / QnetOut

        else 
            QratioIn =  QnetIn  / (QnetIn - QnetOut)
            QratioOut= -QnetOut / (QnetIn - QnetOut)

            ! print *, ' '
            ! print *, 'Qratio ',QratioIn, QratioOut

            do kk=1,max_branch_per_node
                !print *, kk, elemR(JMidx+kk,er_Flowrate)
                if     ((branchsign(kk) * elemR(JMidx+kk,er_Flowrate)) > zeroR) then
                    !% --- inflows
                    dQ(kk) = - resid * QratioIn * elemR(JMidx+kk,er_Flowrate) / QnetIn
                elseif ((branchsign(kk) * elemR(JMidx+kk,er_Flowrate)) < zeroR) then
                    !% --- outflows
                    dQ(kk) = - resid * QratioOut * elemR(JMidx+kk,er_Flowrate) / QnetOut
                elseif (branchsign(kk) * elemR(JMidx+kk,er_Flowrate) == zeroR) then
                    dQ(kk) = zeroR
                else 
                    print *, 'CODE ERROR: unexpected else'
                    call util_crashpoint(739874)
                    return 
                end if
            end do

            !% --- adjust overflow
            dQoverflow = - resid * QratioOut * Qoverflow / QnetOut
        end if

        ! print *, ' '
        ! print *, 'dQ '
        ! do  kk=1,max_branch_per_node
        !     print *, dQ(kk), elemR(JMidx+kk,er_Flowrate)
        ! end do

        !% --- update flowrates
        do kk=1,max_branch_per_node
            elemR(JMidx+kk,er_Flowrate) = elemR(JMidx+kk,er_Flowrate) + dQ(kk)
        end do
        !% --- update overflow
        Qoverflow = Qoverflow + dQoverflow

        !% --- recompute the residual
        resid = junction_conservation_residual (JMidx, Qoverflow, Qstorage)

        ! print *, 'NEW RESIDUAL ',resid 

        ! print *, abs(resid), localEpsilon

        if (abs(resid) > localEpsilon) then 
            print *, 'CODE ERROR: unexpected volume residual'
            call util_crashpoint(6109873)
            return 
        end if

    end subroutine junction_fix_conservation
!%    
!%==========================================================================
!%==========================================================================
!% 
    subroutine junction_branchface_forceJBvalue (frCol, erCol, fiIdx, JMidx, kstart)
        !%------------------------------------------------------------------
        !% Description:
        !% Forces the JB element value on the adjacent face
        !%------------------------------------------------------------------
            integer, intent(in) :: frCol  !% column in faceR array for output
            integer, intent(in) :: erCol  !% column in elemR array for input
            integer, intent(in) :: fiIdx  !% face index column for up/dn map
            integer, intent(in) :: JMidx  !% junction main index
            integer, intent(in) :: kstart !% =1 for upstream, 2 for downstream
            integer :: k1, k2
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        k1 = JMidx + kstart
        k2 = JMidx + max_branch_per_node

        faceR(elemI(k1:k2:2,fiIdx),frCol) = elemR(k1:k2:2,erCol)
    
    end subroutine junction_branchface_forceJBvalue
!%    
!%==========================================================================
!%==========================================================================
!% 





















!%==========================================================================
!%    
!     subroutine junction_get_main_data (jMainR, jMainI thisJM, Npack)
!         !%------------------------------------------------------------------
!         !% Description:
!         !% stores the junction main data (not associated with branches)
!         !%------------------------------------------------------------------
!         !% Declarations
!             real(8), intent(inout) :: jMainR(:,:)
!             integer, intent(inout) :: jMainI(:,:)
!             integer, intent(in)    :: thisJM(:), Npack
!         !%------------------------------------------------------------------
!         !%------------------------------------------------------------------

!         jMainI(1:Npack, jmi_Jtype) = elemSI(thisJM,esi_JunctionMain_Type)    

!         jMainR(1:Npack,jmr_AreaPlan) = elemSI(thisJM,esi_Storage_Plan_Area)
!         jMainR(1:Npack,jmr_cJ)       = zeroR  !% initialized later
!         jMainR(1:Npack,jmr_Hstart)   = elemR(thisJM,er_Head) !% starting head
!         jMainR(1:Npack,jmr_Hdelta)   = zeroR !% initial dH 
!         jMainR(1:Npack,jmr_Hresid)   = zeroR !% initialized later
!         jMainR(1:Npack,jmr_Qlat)     = elemR(thisJM,er_FlowrateLateral)
!         jMainR(1:Npack,jmr_Sc)       = zeroR  !% initialized later

!     end subroutine junction_get_main_data
! !%
! !%==========================================================================
! !%==========================================================================
! !%    
!     subroutine junction_get_branch_data (jBranchR, Hstart)
    !         !%-----------------------------------------------------------------
    !         !% Description:
    !         !% Takes the data from faceR and elemR arrays and stores in the
    !         !% jBranchR 3D array
    !         !%-----------------------------------------------------------------
    !         !% Declarations:
    !             real(8), intent(inout) :: jBranchR(:,:,:)
    !             real(8), intent(in)    :: Hstart(:)
    !             integer, intent(in)    :: thisJM(:)

    !             integer :: kk
    !         !%-----------------------------------------------------------------

    !         !% --- upstream faces
    !         do kk=1,max_branch_per_node,2
    !             where (elemSI(thisJM+kk,esi_JunctionBranch_Exists) ==1)
    !                 jBranchR(:,kk,jbr_beta)     = +oneR !% +1 for upstream beta
    !                 jBranchR(:,kk,jbr_Edelta)   = (faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Head_u)               &
    !                                             + (faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Velocity_u)**twoI) &
    !                                             / (twoR * setting%constant%gravity)) - Hstart
    !                 jBranchR(:,kk,jbr_flowsign) = zeroR !% initialized later                            
    !                 jBranchR(:,kk,jbr_psiL2)    = faceR(elemI(thisJM+kk,ei_Mface_uL),fr_psiL2)
    !                 jBranchR(:,kk,jbr_Q)        = faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Flowrate)
    !                 jBranchR(:,kk,jbr_Qdelta)   = zeroR !% initial zero difference
    !                 jBranchR(:,kk,jbr_Qresid)   = zeroR !% initialized later
                    
    !             elsewhere
    !                 jBranchR(:,kk,jbr_beta)    = zeroR
    !                 jBranchR(:,kk,jbr_Edelta)  = zeroR
    !                 jBranchR(:,kk,jbr_psiL2)   = zeroR
    !                 jBranchR(:,kk,jbr_Q)       = zeroR
    !                 jBranchR(:,kk,jbr_Qdelta)  = zeroR
    !                 jBranchR(:,kk,jbr_Qresid)  = zeroR
                    
    !             endwhere 
    !         end do

    !         !% --- downstream faces
    !         do kk=2,max_branch_per_node,2
    !             where (elemSI(thisJM+kk,esi_JunctionBranch_Exists) ==1)
    !                 jBranchR(:,kk,jbr_beta)    = -oneR !% -1 for downstream beta
    !                 !% --- note, Edelta needs Hstart 
    !                 jBranchR(:,kk,jbr_Edelta)  = (faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Head_d)             &
    !                                             + (faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Velocity_d)**twoI)  &
    !                                             / (twoR * setting%constant%gravity)) - Hstart
    !                 jBranchR(:,kk,jbr_flowsign) = zeroR !% initialized later  
    !                 jBranchR(:,kk,jbr_psiL2)   = faceR(elemI(thisJM+kk,ei_Mface_dL),fr_psiL2)
    !                 jBranchR(:,kk,jbr_Q)       = faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Flowrate)
    !                 jBranchR(:,kk,jbr_Qdelta)  = zeroR !% initial zero difference
    !                 jBranchR(:,kk,jbr_Qresid)  = zeroR !% initialized later
                    
    !             elsewhere
    !                 jBranchR(:,kk,jbr_beta)    = zeroR
    !                 jBranchR(:,kk,jbr_Edelta)  = zeroR
    !                 jBranchR(:,kk,jbr_psiL2)   = zeroR
    !                 jBranchR(:,kk,jbr_Q)       = zeroR
    !                 jBranchR(:,kk,jbr_Qdelta)  = zeroR
    !                 jBranchR(:,kk,jbr_Qresid)  = zeroR
                    
    !             endwhere 
    !         end do

!     end subroutine junction_get_branch_data
! !%
! !%==========================================================================
! !%==========================================================================
! !%  
!     function junction_flowsign (Flowrate, dQ, nk)
    !         !%-----------------------------------------------------------------
    !         !% Description
    !         !% returns +- oneR depending on flow sign
    !         !%-----------------------------------------------------------------
    !         !% Declarations:
    !             integer, intent(in)    :: nk
    !             real(8), intent(in)    :: Flowrate(:), dQ(:)
    !             real(8), dimension(nk) :: oneOut = oneR

    !             real(8), dimension(nk) :: junction_flowsign
    !         !%-----------------------------------------------------------------

    !         junction_flowsign = sign(oneOut,Flowrate + dQ)
            
!     end function junction_flowsign
! !%
! !%==========================================================================
! !%==========================================================================
! !%    
!     subroutine junction_cJ (AreaPlan, jMainI, thisJM)
    !         !%-----------------------------------------------------------------
    !         !% Description:
    !         !% Sets the jmr_cJ value that depends on storage type for all junctions
    !         !%-----------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(inout) :: AreaPlan(:)
    !             integer, intent(in)    :: jMainI(:,:), thisJM(:)
    !         !%-----------------------------------------------------------------

    !         where (jMainI(:,jmi_Jtype) .eq. ImpliedStorage)
    !             AreaPlan = zeroR
    !         elsewhere
    !             AreaPlan = elemSR(thisJM,esr_Storage_Plan_Area)
    !         endwhere

!     end subroutine junction_cJ
! !%
! !%==========================================================================
! !%==========================================================================
! !%  
!     subroutine junction_Sc (Sc, beta, Qinitial)
    !         !%-----------------------------------------------------------------
    !         !% Description
    !         !% Computes the source term for continuity as the net sum
    !         !% of the initila flowrate
    !         !%-----------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(inout) :: Sc(:)
    !             real(8), intent(in)    :: beta(:,:)
    !             real(8), intent(in)    :: Qinitial(:,:)
    !         !%-----------------------------------------------------------------
    !         !%-----------------------------------------------------------------

    !         Sc = sum(beta * Qinitial,2)    

!     end subroutine junction_Sc
! !%
! !%==========================================================================
! !%==========================================================================
! !%  
!     subroutine junction_Sm (jBranchR)
    !         !%-----------------------------------------------------------------
    !         !% Description:
    !         !% sets the juction source term
    !         !%-----------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(inout) :: jBranchR(:,:,:)
    !         !%-----------------------------------------------------------------

    !         where (jBranchR(:,:,jbr_beta)  .ne. zeroI)
    !             jBranchR(:,:,jbr_Sm)  = jBranchR(:,:,jbr_a) * (jBranchR(:,:,jbr_Qinit)**twoI) &
    !                                   - jBranchR(:,:,jbr_Edelta)
    !         endwhere 

!     end subroutine junction_Sm    
! !%
! !%==========================================================================
! !%==========================================================================
! !%  
!     subroutine junction_Qresid_initial (Qresid, Sm)
    !         !%-----------------------------------------------------------------
    !         !% Description
    !         !% Sets the initial residual for the flowrate on the branches
    !         !%-----------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(inout) :: Qresid(:,:)
    !             real(8), intent(in)    :: Sm(:,:)
    !         !%-----------------------------------------------------------------

    !         Qresid = Sm

!     end subroutine junction_Qresid_initial
! !%
! !%==========================================================================
! !%==========================================================================
! !% 
!     subroutine junction_Hresid_initial (Hresid, Sc)
    !         !%-----------------------------------------------------------------
    !         !% Description
    !         !% Sets the initial residual for the flowrate on the branches
    !         !%-----------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(inout) :: Hresid(:,:)
    !             real(8), intent(in)    :: Sc(:,:)
    !         !%-----------------------------------------------------------------

    !         Hresid = Sc

!     end subroutine junction_Hresid_initial
! !%
! !%==========================================================================

! !%                              OLD BELOW

! !%==========================================================================
! !%
!     subroutine  junction_calculationOLD (Npack, thisColP, istep)
    !         !%-----------------------------------------------------------------
    !         !% Description:
    !         !% get the jacobean matrix for junction element following the 
    !         !% derivation by ...
    !         !%-----------------------------------------------------------------
    !         !% Declarations
    !         integer, intent(in) :: Npack, thisColP, istep
    !         integer, pointer    :: thisJM(:)

    !         real(8), pointer :: grav

    !         !% --- the junction data storage
    !         real(8) :: jDataR(Npack, max_branch_per_node, NCol_jDataR) 
    !         integer :: jDataI(Npack, max_branch_per_node, NCol_jDataI)

    !         !% --- solver data storage
    !         real(8), dimension(Npack) :: startEJ, EJ, Qlat, VolumeEJ, fullVolumeEJ
    !         real(8), dimension(Npack) :: Zcrown, Aplan, Aponded
    !         real(8), dimension(Npack) :: OverflowLength, OverflowHeight, Qoverflow

    !         real(8), dimension(Npack) :: alpha, normL2, Qin, Qout,  fcoef
    !         real(8), dimension(Npack) :: deltaEJ, residEJ, Volume

    !         integer, dimension(Npack) :: np, Niter, OverFlowType, JunctionType

    !         logical, dimension(Npack) :: isconverged, isRepack, isfailed, isFull, isSubmergedOrifice
    !         integer :: kk, mm

    !         !% --- epsilon used in junction to prevent zero values
    !         real(8) :: jEpsilon = 1.d-6

    !         !% --- maximum number of iterations (HACK -- make into global)
    !         integer :: maxIter = 20
    !         !% --- flowrate conservation convergence (HACK -- make into global)
    !         real(8) :: jConvergence = 1.d-12

    !         !%-----------------------------------------------------------------------------
    !         !%  Aliases:         
    !             grav         => setting%constant%gravity
    !             !% --- packed set of junctions
    !             thisJM => elemP(1:Npack,thisColP)
    !         !%-----------------------------------------------------------------------------
    !         !% Preliminaries
    !             !% --- iteration counter
    !             Niter = zeroI

    !             !% --- initialization
    !             alpha   = zeroR
    !             normL2  = zeroR 
    !             Qin     = zeroR 
    !             Qout    = zeroR 
    !             fcoef   = zeroR
    !             deltaEJ = zeroR
    !             residEJ = zeroR
    !             Qoverflow = zeroR

    !             jDataR = zeroR
    !             jDataI = zeroI


    !             !% --- stores the initial elemR k location
    !             !%     this is used to bring packed branch data back to element and face array
    !             do concurrent (kk=1:max_branch_per_node)
    !                 jDataI(:,kk,ji_kidx)  = kk  
    !             end do

    !             !% store the starting head in Temp array for conservation check at end
    !             elemR(thisJM,er_Temp01) = elemR(thisJM,er_Head)
    !         !%-----------------------------------------------------------------------------

    !         !% --- get the junction data 
    !         call junction_get_JM_dataOLD (startEJ, EJ, Qlat, volumeEJ, fullVolumeEJ, &
    !                 Zcrown, Aplan, Aponded, OverflowLength, OverFlowHeight, &
    !                 isFull, isSubmergedOrifice, JunctionType, OverflowType, thisJM, Npack)

    !         ! print *, ' '
    !         ! print *, 'AT START OF JUNCTION SOLUTION FOR step ',istep
    !         ! print *, 'StartEJ ',startEJ(1), elemR(thisJM(1),er_Head)
    !         ! print *, 'EJ      ',EJ(1), elemR(thisJM(1),er_Head)
    !         ! print *, 'Zcrown  ',Zcrown(1), elemR(thisJM(1),er_Zcrown)
    !         ! print *, 'Zbottom ','                    ',elemR(thisJM(1),er_Zbottom)
    !         ! print *, 'Depth   ','                    ',elemR(thisJM(1),er_Depth)
    !         ! print *, 'QLat    ',Qlat(1), elemR(thisJM(1),er_FlowrateLateral)
    !         ! print *, 'Vej     ',volumeEJ(1), elemR(thisJM(1),er_Volume)
    !         ! print *, 'fullV   ',fullVolumeEJ(1), elemR(thisJM(1),er_FullVolume)
    !         ! print *, 'Apond   ', Aponded(1), elemSR(thisJM(1),esr_JunctionMain_PondedArea)
    !         ! print *, 'Aplan   ', Aplan(1),elemSR(thisJM(1),esr_Storage_Plan_Area)
    !         ! print *, 'Jtype   ', JunctionType(1), trim(reverseKey(JunctionType(1)))
    !         ! print *, 'overflowtype ',OverflowType(1), trim(reverseKey(OverflowType(1)))
    !         ! print *, ' '        

    !         !% --- store the face and elem data in 3D array
    !         call junction_get_face_data(jDataR, thisJM)

    !         ! print *, 'Face  DATA'
    !         ! !print *, 'size jDataR ',size(JdataR,1), size(JdataR,2), size(JdataR,3)
    !         ! do kk=1,max_branch_per_node
    !         !     !print *, 'kk ',kk
    !         !     if (JdataR(1,kk,jr_beta) .ne. zeroI) then 
    !         !         print *, ' '
    !         !         print *, kk,'beta  ',JdataR(1,kk,jr_beta)
                    
    !         !         if (mod(kk,2) == 0) then 
    !         !             print *, 'Area  ',JdataR(1,kk,jr_Area),    faceR(elemI(thisJM(1)+kk,ei_Mface_dL),fr_Area_d)
    !         !             print *, 'E     ',jDataR(1,kk,jr_Ebranch), faceR(elemI(thisJM(1)+kk,ei_Mface_dL),fr_Head_d)
    !         !             print *, 'length',jDataR(1,kk,jr_Length),  faceR(elemI(thisJM(1)+kk,ei_Mface_dL),fr_Length_d)
    !         !             print *, 'gamma ',jDataR(1,kk,jr_Gamma),   faceR(elemI(thisJM(1)+kk,ei_Mface_dL),fr_GammaM)
    !         !             print *, 'Q     ',jDataR(1,kk,jr_Q),       faceR(elemI(thisJM(1)+kk,ei_Mface_dL),fr_Flowrate)
    !         !             print *, 'K     ',jDataR(1,kk,jr_K),       faceR(elemI(thisJM(1)+kk,ei_Mface_dL),fr_KJunction_MinorLoss)
                        
    !         !         else 
    !         !             print *, 'Area  ',JdataR(1,kk,jr_Area),    faceR(elemI(thisJM(1)+kk,ei_Mface_uL),fr_Area_u)
    !         !             print *, 'E     ',jDataR(1,kk,jr_Ebranch), faceR(elemI(thisJM(1)+kk,ei_Mface_uL),fr_Head_u)
    !         !             print *, 'length',jDataR(1,kk,jr_Length),  faceR(elemI(thisJM(1)+kk,ei_Mface_uL),fr_Length_u)
    !         !             print *, 'gamma ',jDataR(1,kk,jr_Gamma),   faceR(elemI(thisJM(1)+kk,ei_Mface_uL),fr_GammaM)
    !         !             print *, 'Q     ',jDataR(1,kk,jr_Q),       faceR(elemI(thisJM(1)+kk,ei_Mface_uL),fr_Flowrate)
    !         !             print *, 'K     ',jDataR(1,kk,jr_K),       faceR(elemI(thisJM(1)+kk,ei_Mface_uL),fr_KJunction_MinorLoss)
    !         !         end if
    !         !     end if
    !         ! end do

    !         ! print *, ' '
    !         ! print *, 'heads '
    !         ! print *, elemR(faceI(elemI(thisJM(1)+1,ei_Mface_uL),fi_Melem_uL),er_Head)
    !         ! print *, faceR(elemI(thisJM(1)+1,ei_Mface_uL),fr_Head_u)
    !         ! print *, faceR(elemI(thisJM(1)+1,ei_Mface_uL),fr_Head_d)
    !         ! print *, elemR(thisJM(1),er_Head)
    !         ! print *, faceR(elemI(thisJM(1)+2,ei_Mface_dL),fr_Head_u)
    !         ! print *, faceR(elemI(thisJM(1)+2,ei_Mface_dL),fr_Head_d)
    !         ! print *, elemR(faceI(elemI(thisJM(1)+2,ei_Mface_dL),fi_Melem_dL),er_Head)

    !         ! print *, ' '
    !         ! print *, 'flowrates'
    !         ! print *, elemR(faceI(elemI(thisJM(1)+1,ei_Mface_uL),fi_Melem_uL),er_Flowrate)
    !         ! print *, faceR(elemI(thisJM(1)+1,ei_Mface_uL),fr_Flowrate)
    !         ! !print *, faceR(elemI(thisJM(1)+1,ei_Mface_uL),fr_Head_d)
    !         ! print *, elemR(thisJM(1),er_Flowrate)
    !         ! print *, faceR(elemI(thisJM(1)+2,ei_Mface_dL),fr_Flowrate)
    !         ! !print *, faceR(elemI(thisJM(1)+2,ei_Mface_dL),fr_Head_d)
    !         ! print *, elemR(faceI(elemI(thisJM(1)+2,ei_Mface_dL),fi_Melem_dL),er_Flowrate)

    !         !stop 498734

    !         !% --- set the alpha at the junction
    !         call junction_alpha (alpha, Aplan, Aponded, isFull, JunctionType, &
    !              OverflowType, thisJM, Npack, istep)

    !         ! print *, ' '
    !         ! print *, 'Alpha ',alpha(1)
    !         ! print *, ' '

    !         !% --- set the fcoef for overflows
    !         call junction_fcoef(fcoef, Aplan, OverflowLength, OverflowType)

    !         ! print *, 'Fcoef ',fcoef(1)
    !         ! print *, ' '

    !         !% AFTER THIS POINT WE SHOULD NOT NEED TO USE thisJM(:) UNTIL
    !         !% MOVING DATA BACK TO faceR AND elemR

    !         !% --- identify branches with gamma = 0 and K=0 or gamma = 0 and Q = 0
    !         !%     requires setting beta to zero (thus, we cannot compute flow in 
    !         !%     this branch during this time step)
    !         call junction_initialize_beta (jDataR, jEpsilon)

    !         ! print *, ' '
    !         ! print *, '=========================================='
    !         ! print *, 'Beta '
    !         ! print *, jDataR(1,:,jr_beta)
    !         ! print *, ' '
    !         ! print *, 'Gamma '
    !         ! print *, jDataR(1,:,jr_Gamma)
    !         ! print *, ' '
    !         ! print *, 'K '
    !         ! print *, jDataR(1,:,jr_K)
    !         ! print *, ' '
    !         ! print *, 'Q '
    !         ! print *, jDataR(1,:,jr_Q)
    !         ! print *, ' '

    !         !% --- the expected total branches
    !         !%     This will be modified by packing if there are branches with beta=0.
    !         np = elemSI(thisJM,esi_JunctionMain_Total_Branches)     

    !         ! print *, 'np ',np
    !         ! print *, ' '
            
    !         !% --- code check for debugging
    !         call junction_check_beta (jDataR, np, nPack)

    !         !% --- pack the jDataR and jDataI to remove zero branches
    !         isRepack(:) = .true.
    !         call junction_pack_all (jDataR, jDataI, np, isRepack, Npack)

    !         ! print *, 'repack Ebranch'
    !         ! print *, jDataR(1,:,jr_Ebranch)
    !         ! print *, ' '
    !         ! print *, 'repack Area '
    !         ! print *, jDataR(1,:,jr_Area)

    !         !% --- compute the invariant terms that will not change in this time step
    !         call junction_invariant_terms (jDataR, np, Npack)

    !         ! print *, 'jA '
    !         ! print *, jDataR(1,:,jr_a)
    !         ! print *, ' '
    !         ! print *, 'jC' 
    !         ! print *, jDataR(1,:,jr_c)
    !         ! print *, ' '
    !         ! print *, 'lambda A'
    !         ! print *, jDataR(1,:,jr_LambdaA)
    !         ! print *, ' '
    !         ! print *, 'lambda B'
    !         ! print *, jDataR(1,:,jr_LambdaB)
    !         ! print *, ' '


    !         !% --- compute the residual in each branch
    !         call junction_all_residuals (jDataR, residEJ, EJ, &
    !             startEJ, alpha, fcoef, Qlat, Zcrown, OverflowHeight, &
    !             OverflowType, isFull, isSubmergedOrifice, np, Npack)

    !         ! print *, 'residEJ BB', residEJ(1)
    !         ! print *, ' '
            
    !         !stop 298734

    !         isconverged = .false. 
    !         isfailed    = .false.

    !         !do concurrent (mm=1:Npack)
    !         do mm=1,1 !%Npack
    !                 ! print *, 'in junction loop ',mm
    !                 ! print *, 'resid at start   ',residEJ(mm)
    !             do while (.not. isconverged(mm))
    !                     ! print *, '============================================================='
    !                     ! print *, 'Niter ',Niter(mm), '===',mm,'=================================='
    !                     ! print *, ' '
    !                     !% --- increment iteration counter for exit after maxIter loops
    !                 Niter(mm) = Niter(mm) + 1

    !                 !% --- compute varying lambdaC 
    !                 jDataR(mm,:,jr_LambdaC) = junction_lambdaC (jDataR(mm,:,:), np(mm))

    !                 !% --- look for small LambdaC and repack if needed
    !                 call junction_check_repack_oneJ (jDataR(mm,:,:), np(mm), jEpsilon, isRepack(mm))
            
    !                 ! print *, 'lambda C 1',jDataR(mm,1:2,jr_LambdaC)

    !                 if (isRepack(mm)) then 
    !                     call junction_pack_oneJ &
    !                         (jDataR(mm,:,:), jDataI(mm,:,:), np(mm), isRepack(mm))
    !                 end if

    !                 ! print *, 'lambda C 2',jDataR(mm,1:2,jr_LambdaC)

    !                 !% --- compute the junction iteration
    !                 call junction_iteration &
    !                     (jDataR(mm,:,:), deltaEJ(mm), EJ(mm),  residEJ(mm), alpha(mm), &
    !                      fcoef(mm), Zcrown(mm), OverflowHeight(mm), OverflowType(mm), isFull(mm), &
    !                      isSubmergedOrifice(mm), np(mm)) 

    !                 ! print *,'new Q     ', jDataR(mm,1,jr_Q),    jDataR(mm,2,jr_Q)
    !                 ! print *, 'oldresid ', jDataR(mm,1,jr_residQ), jDataR(mm,2,jr_residQ)

    !                 !% --- compute the residuals
    !                 jDataR(mm,1:np(mm),jr_residQ) = junction_Qresidual              &
    !                                              (jDataR(mm,:,:), EJ(mm), startEJ(mm), np(mm),mm)
        
    !                 ! print *, 'newQ resid',jDataR(mm,1,jr_residQ) ,jDataR(mm,2,jr_residQ)                                     

    !                 residEJ(mm) = junction_Eresidual                                           &
    !                              (jDataR(mm,:,:), EJ(mm), startEJ(mm), alpha(mm), fcoef(mm),   &
    !                               Qlat(mm), Zcrown(mm), OverflowHeight(mm), isFull(mm),        &
    !                               isSubmergedOrifice(mm), OverflowType(mm), np(mm))

    !                 !% --- compute residual L2 norm
    !                 normL2(mm) = junction_get_L2norm &
    !                             (jDataR(mm,:,jr_residQ), residEJ(mm), np(mm))

    !                 ! print *, 'residEJ ',residEJ(mm)
    !                 ! print *, 'normL2  ',normL2(mm), jConvergence

    !                 !% --- check norms for exit
    !                 if (normL2(mm) .le. jConvergence) then 
    !                     isconverged(mm) = .true.
    !                     ! print *, ' ' 
    !                     ! print *, 'CONVERGED'
    !                     ! print *, ' '
    !                 else
    !                     if (Niter(mm)+1 .ge. maxIter) then 
    !                         !% --- mark as failed and exit
    !                         isconverged(mm) = .true.
    !                         isfailed(mm)    = .false.
    !                         !print *, ' '
    !                         print *, 'FAILED '
    !                         !print *, ' '
    !                     end if
    !                 end if 


                    
    !             end do

    !             !% ---- CONVERGED OR EXITING FOR THIS JUNCTION------------------------------
    !             !% --- compute the net Q inflow
    !             Qin(mm) = junction_branch_Qnet (jDataR(mm,:,:), np(mm), +oneI)
    !             !% --- compute the net Q outflow
    !             Qout(mm)= junction_branch_Qnet (jDataR(mm,:,:), np(mm), -oneI)
    !             !% --- compute the Q overflow
    !             Qoverflow(mm) = junction_Qoverflow (EJ(mm), Zcrown(mm), fcoef(mm),            &
    !                 OverflowHeight(mm), OverflowType(mm), isFull(mm), isSubmergedOrifice(mm) )

    !             !% --- adjust Qs to ensure conservation
    !             jDataR(mm,1:np(mm),jr_Q) = junction_conservation &
    !                 (jDataR(mm,:,:), residEJ(mm), Qin(mm), Qout(mm), np(mm))

    !             !% --- restore data to the faceR ... arrays
    !             call junction_push_face_flowrate &
    !                 (faceR(:,fr_Flowrate), faceR(:,fr_Velocity_u), faceR(:,fr_Velocity_d), &
    !                  faceR(:,fr_Area_u), faceR(:,fr_Area_d), jDataR(mm,:,:), jDataI(mm,:,:), &
    !                  np(mm), thisJM(mm))

    !             !% --- head
    !             elemR(thisJM(mm),er_Head) = EJ(mm)

    !             !% --- volume
    !             elemR(thisJM(mm),er_Volume) = junction_volume &
    !                 (elemR(thisJM(mm),er_Volume_N0), Qin(mm), Qout(mm), Qlat(mm), Qoverflow(mm), istep) 

    !             !% --- set the junction overflow and ponded volume on the last step of RK
    !             if (istep == 2) then
    !                 !% --- overflow
    !                 elemR(thisJM(mm),er_VolumeOverFlow) = Qoverflow(mm) * setting%Time%Hydraulics%Dt

    !                 !% --- ponded volume
    !                 if ((OverflowType(mm) == Ponded) .and. (isFull(mm))) then 
    !                     elemR(thisJM(mm),er_VolumePonded) = elemR(thisJM(mm),er_Volume) - elemR(thisJM(mm),er_FullVolume)
    !                 end if
    !             end if

    !         end do

    !         !% --- check conservation -- can be commented after debugging
    !         call junction_check_conservation (thisJM, Qoverflow, jDataR, Npack ,istep)

    !         !stop 6098734

    !         !% HACK --- the new water surface elevation associated with the
    !         !% changing storage may NOT exactly balance conservation because of
    !         !% the use of a fixed planar area. We need to rederiv the 
    !         !% junction_conservation equations so that the exact change in volume
    !         !% implied by the change in head gives us the correct volume

    !         !% --- HACK: NEED TO DO SOMETHING FOR REPORTING CONVERGENCE FAILURE 

!     end subroutine junction_calculationOLD
! !%
! !%========================================================================== 
! !%========================================================================== 
! !% 
!     pure subroutine junction_get_JM_dataOLD &
    !         (startEJ, EJ, Qlat, volumeEJ, fullVolumeEJ, Zcrown, Aplan, Aponded, &
    !             OverflowLength, OverflowHeight, &
    !             isFull, isSubmergedOrifice, JunctionType, OverflowType, thisJM, Npack)
    !         !%------------------------------------------------------------------
    !         !% Description:
    !         !% Gets the JM data from the elemR array and stores in local
    !         !% vector
    !         !%------------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(inout) :: startEJ(:), EJ(:), Qlat(:)
    !             real(8), intent(inout) :: volumeEJ(:), fullVolumeEJ(:), Zcrown(:)
    !             real(8), intent(inout) :: Aplan(:), Aponded(:)
    !             real(8), intent(inout) :: OverflowLength(:), OverflowHeight(:)
    !             logical, intent(inout) :: isFull(:), isSubmergedOrifice(:)
    !             integer, intent(inout) :: OverflowType(:), JunctionType(:)
    !             integer, intent(in)    :: Npack, thisJM(:)
    !         !%------------------------------------------------------------------

    !         startEJ     (1:Npack) = elemR(thisJM,er_Head)
    !         EJ          (1:Npack) = elemR(thisJM,er_Head)
    !         Qlat        (1:Npack) = elemR(thisJM,er_FlowrateLateral)
    !         volumeEJ    (1:Npack) = elemR(thisJM,er_Volume)
    !         fullVolumeEJ(1:Npack) = elemR(thisJM,er_FullVolume)
    !         Zcrown      (1:Npack) = elemR(thisJM,er_Zcrown)
    !         OverflowType(1:Npack) = elemSI(thisJM,esi_JunctionMain_OverflowType)
    !         Aponded     (1:Npack) = elemSR(thisJM,esr_JunctionMain_PondedArea)
    !         Aplan       (1:Npack) = elemSR(thisJM,esr_Storage_Plan_Area)
    !         JunctionType(1:Npack) = elemSI(thisJM,esi_JunctionMain_Type)

    !         OverflowLength(1:Npack) = elemSR(thisJM,esr_JunctionMain_OverflowOrifice_Length)
    !         OverflowHeight(1:Npack) = elemSR(thisJM,esr_JunctionMain_OverflowOrifice_Height)

    !         where (JunctionType .eq. ImpliedStorage)
    !             !% --- The only use of Aplan with implied storage is for the
    !             !%     overflow weir effect.
    !             Aplan = 0.179d0 !% HACK Area gives Lw = 1.5 m as overflow weir length
    !             where (elemR(thisJM,er_Head) .ge. Zcrown)
    !                 isFull = .true.
    !             elsewhere
    !                 isFull = .false.
    !             endwhere
    !         elsewhere
    !             where (volumeEJ .ge. fullVolumeEJ)
    !                 isFull = .true.
    !             elsewhere 
    !                 isFull = .false.
    !             endwhere
    !         endwhere 

    !         where ((isFull) .and. (OverflowType .eq. OverflowOrifice) &
    !             .and. (EJ > Zcrown + OverflowHeight))
    !             isSubmergedOrifice = .true.
    !         elsewhere 
    !             isSubmergedOrifice = .false.
    !         endwhere

   
!     end subroutine junction_get_JM_dataOLD
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%    
!     pure subroutine junction_get_face_data (jDataR, thisJM )    
    !         !%------------------------------------------------------------------
    !         !% Takes the data from faceR and elemR arrays and stores in the
    !         !% jDataR 3D array
    !         !%------------------------------------------------------------------
    !             real(8), intent(inout) :: jDataR(:,:,:)
    !             integer, intent(in)    :: thisJM(:)

    !             integer :: kk
    !         !%------------------------------------------------------------------
    !         !%------------------------------------------------------------------

    !         !% --- store data for upstream faces
    !         do concurrent (kk=1:max_branch_per_node:2)
    !             where (elemSI(thisJM+kk,esi_JunctionBranch_Exists) == oneI)
    !                 jDataR(:,kk,jr_beta)    = +oneR !% +1 for upstream beta
    !                 jDataR(:,kk,jr_Area)    =  faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Area_u)
    !                 jDataR(:,kk,jr_Ebranch) =  faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Head_u)            &
    !                                         + (faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Velocity_u)**twoI) &
    !                                         / (twoR * setting%constant%gravity)
    !                 jDataR(:,kk,jr_Gamma)   =  faceR(elemI(thisJM+kk,ei_Mface_uL),fr_GammaM)
    !                 jDataR(:,kk,jr_Length)  =  faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Length_u)
    !                 jDataR(:,kk,jr_Q)       =  faceR(elemI(thisJM+kk,ei_Mface_uL),fr_Flowrate)
    !                 jDataR(:,kk,jr_K)       =  faceR(elemI(thisJM+kk,ei_Mface_uL),fr_KJunction_MinorLoss)
    !             elsewhere
    !                 jDataR(:,kk,jr_beta)    = zeroR !%  0 if branch does not exist
    !             endwhere
    !         end do

    !         !% --- store data for downstream faces
    !         do concurrent (kk=2:max_branch_per_node:2)
    !             where (elemSI(thisJM+kk,esi_JunctionBranch_Exists) == oneI)
    !                 jDataR(:,kk,jr_beta)    = -oneR !% -1 for downstream beta
    !                 jDataR(:,kk,jr_Area)    =  faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Area_d)
    !                 jDataR(:,kk,jr_Ebranch) =  faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Head_d)             &
    !                                         + (faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Velocity_d)**twoI)  &
    !                                         / (twoR * setting%constant%gravity)
    !                 jDataR(:,kk,jr_Gamma)   =  faceR(elemI(thisJM+kk,ei_Mface_dL),fr_GammaM)
    !                 jDataR(:,kk,jr_Length)  =  faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Length_d)
    !                 jDataR(:,kk,jr_Q)       =  faceR(elemI(thisJM+kk,ei_Mface_dL),fr_Flowrate)
    !                 jDataR(:,kk,jr_K)       =  faceR(elemI(thisJM+kk,ei_Mface_dL),fr_KJunction_MinorLoss)
    !             elsewhere
    !                 jDataR(:,kk,jr_beta)    = zeroR !%  0 if branch does not exist
    !             endwhere
    !         end do


!     end subroutine junction_get_face_data
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%      
!     subroutine junction_alpha (alpha, Aplan, Aponded, isFull, JunctionType,  &
    !                 OverflowType, thisJM, Npack, istep)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% Computes the alpha term that represents change in junction storage
    !         !%------------------------------------------------------------------
    !         !% Declarations
                
    !             real(8), intent(inout) :: alpha(:)
    !             real(8), intent(in)    :: Aplan(:), Aponded(:)
    !             integer, intent(in)    :: JunctionType(:), OverflowType(:), thisJM(:), Npack, istep
    !             logical, intent(in)    :: isFull(:)

    !             real(8), pointer :: crk(:), dt
    !             integer          :: mm
    !         !%------------------------------------------------------------------
    !         !% Aliases
    !             crk          => setting%Solver%crk2
    !             dt           => setting%Time%Hydraulics%Dt
    !         !%------------------------------------------------------------------

    !         do mm=1,Npack
    !             select case (JunctionType(mm))

    !                 case (ImpliedStorage)
    !                     !% --- no volume increase allowed
    !                     alpha(mm) = zeroR

    !                 case (FunctionalStorage, TabularStorage)
    !                     if (.not. isFull(mm)) then
    !                         !% --- below full volume
    !                         alpha(mm) = Aplan(mm) / (crk(istep)*dt)
    !                     else
    !                         !% --- at or above full volume
    !                         select case (OverflowType(mm))
    !                             case (NoOverflow, OverflowOrifice, OverflowWeir)
    !                                 !% --- alpha remains based on plane area
    !                                 alpha(mm) = Aplan(mm) / (crk(istep)*dt)
    !                             case (Ponded)
    !                                 !% --- alpha used ponded area instead of plane area
    !                                 alpha(mm) = Aponded(mm)  / (crk(istep)*dt)
    !                             case default 
    !                                 print *, 'CODE ERROR: Unexpected case default'
    !                                 call util_crashpoint(2287422)
    !                         end select
    !                     end if
    !                 case default
    !                     print *, 'unexpected case default'
    !                     call util_crashpoint(628733)
    !             end select  
    !         end do

!     end subroutine junction_alpha
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     pure subroutine junction_fcoef (fcoef, Aplan, OverflowLength, OverflowType)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% Computes the fcoef term for overflow weirs
    !         !%------------------------------------------------------------------
    !         !% Declarations   
    !             real(8), intent(inout) :: fcoef(:)
    !             real(8), intent(in)    :: Aplan(:), OverflowLength(:)
    !             integer, intent(in)    :: OverflowType(:)
    !             real(8), parameter :: Cbc = 1.45d0  !% weir coef. HACK MOVE TO SETTINGS
    !         !%------------------------------------------------------------------

    !         where (OverflowType == OverflowWeir)
    !             fcoef = twoR * Cbc * sqrt( setting%Constant%pi * Aplan)
    !         elsewhere (OverflowType == OverflowOrifice)
    !             fcoef = OverflowLength * sqrt(twoR * setting%Constant%pi)
    !         elsewhere
    !             fcoef = zeroR
    !         endwhere

!     end subroutine junction_fcoef    
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     pure subroutine junction_initialize_beta (jDataR, jEpsilon)
    !         !%------------------------------------------------------------------
    !         !% Description:
    !         !% sets the initial jDataR(:,jr_beta) = 0 for branches that
    !         !% cannot have a flow in this time step.
    !         !%------------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(inout) :: jDataR(:,:,:)
    !             real(8), intent(in)    :: jEpsilon
    !         !%------------------------------------------------------------------
    !         !%------------------------------------------------------------------

    !         where (((jDataR(:,:,jr_Gamma)     < jEpsilon) .and.   &
    !                 (jDataR(:,:,jr_K)         < jEpsilon)       ) &
    !             .or.                                              &
    !                ((jDataR(:,:,jr_Gamma)     < jEpsilon) .and.   & 
    !                 (abs(jDataR(:,:,jr_Q))    < jEpsilon)       ) ) 

    !             jDataR(:,:,jr_beta) = zeroR

    !         endwhere

!     end subroutine junction_initialize_beta
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     subroutine junction_check_beta (jR, np, nPack)
            ! !%------------------------------------------------------------------
            ! !% Description:
            ! !% Checks that np(:) is identical to the number of beta <> 0
            ! !% This can be commented out once the code is thoroughly debugged
            ! !%------------------------------------------------------------------
            !     real(8), intent(in) :: jR(:,:,:)
            !     integer, intent(in) :: np(:), nPack 
            !     integer :: mm
            !     real(8) :: nbeta

            ! !%------------------------------------------------------------------
            ! !%------------------------------------------------------------------

            ! do mm=1,Npack 
            !     nbeta = count(jR(mm,:,jr_beta) .ne. zeroR)
            !     if (np(mm) .ne. nbeta) then 
            !         print *, 'CODE ERROR:'
            !         print *, 'problem with np and nbeta in junction iteration'
            !         print *, 'mm     = ',mm 
            !         print *, 'np(mm) = ',np(mm)
            !         print *, 'n beta = ',nbeta
            !         call util_crashpoint(398233)
            !     end if
            ! end do

!     end subroutine junction_check_beta
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%
!     pure subroutine junction_invariant_terms (jR, np, Npack)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% computes the terms in the junction solution that do not change
    !         !% during a time step
    !         !%------------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(inout) :: jR(:,:,:)
    !             integer, intent(in)    :: np(:), nPack

    !             integer :: mm
    !         !%------------------------------------------------------------------

    !         do concurrent (mm=1:Npack)
    !             !% --- a
    !             jR(mm,1:np(mm),jr_a) = jR(mm,1:np(mm),jr_beta)             &
    !                 * jR(mm,1:np(mm),jr_Length) * jR(mm,1:np(mm),jr_Gamma) &
    !                 / (twoR * setting%constant%gravity * jR(mm,1:np(mm),jr_Area))

    !             !% --- c
    !             jR(mm,1:np(mm),jr_c) = jR(mm,1:np(mm),jr_beta)     &
    !                 * jR(mm,1:np(mm),jr_K)                             &
    !                 / (twoR * setting%constant%gravity * (jR(mm,1:np(mm),jr_Area)**twoI))   

    !             !% --- lambdaA
    !             jR(mm,1:np(mm),jr_LambdaA)                            &
    !                 = twoR * setting%constant%gravity  *  (jR(mm,1:np(mm),jr_Area)**twoI) 

    !             !% --- lambdaB
    !             jR(             mm,1:np(mm),jr_LambdaB)               &
    !                 =        jR(mm,1:np(mm),jr_Area)                  &
    !                        * jR(mm,1:np(mm),jr_Length)                &
    !                        * jR(mm,1:np(mm),jr_Gamma)
    !         end do

!     end subroutine junction_invariant_terms
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     pure function junction_lambdaC (jR, np)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% computes the lambdaC term that changes with each iteration
    !         !%------------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(in)    :: jR(:,:)
    !             integer, intent(in)    :: np

    !             real(8), dimension(np) :: junction_lambdaC

    !             integer :: mm
    !         !%------------------------------------------------------------------

    !         !% --- note that the LambdaC storage is later used for Lambda = LambdaA/LambdaC
    !         junction_lambdaC = jR(1:np,jr_LambdaB) + jR(1:np,jr_K) * abs(jR(1:np,jr_Q))

!     end function junction_lambdaC
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     subroutine junction_iteration (jR, deltaEJ, EJ, residEJ, alpha, &
    !             fcoef, Zcrown, OverflowHeight, OverflowType, isFull, isSubmergedOrifice, np) 
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% Iterative results for one step advancing Q and EJ
    !         !%------------------------------------------------------------------
    !             real(8), intent(inout) :: jR(:,:), deltaEJ, EJ
    !             real(8), intent(in)    :: residEJ, alpha, fcoef, Zcrown, OverflowHeight
    !             logical, intent(in)    :: isFull, isSubmergedOrifice
    !             integer, intent(in)    :: np, OverFlowType

    !             real(8) :: Bsum, Csum, Gterm

    !             integer :: mm
    !         !%------------------------------------------------------------------
    !         !%------------------------------------------------------------------

    !         !% --- compute lambda and store in lambdaC (recomputed each iteration)
    !         jR(1:np,jr_LambdaC) =  jR(1:np,jr_LambdaA)  / jR(1:np,jr_LambdaC)

    !         ! print *, 'Lambda ratio ',jR(1:np,jr_LambdaC)
            
    !         !% B summation
    !         Bsum = sum(jR(1:np,jr_LambdaC)) 

    !         ! print *, 'Bsum ',Bsum
            
    !         !% C summation
    !         Csum = sum(jR(1:np,jr_LambdaC) * jR(1:np,jr_residQ))

    !         ! print *, 'Csum ',Csum

    !         !% Gterm
    !         Gterm = junction_Gterm (alpha, EJ, Zcrown, fcoef, OverflowHeight, &
    !                                 OverflowType, isFull, isSubmergedOrifice)

    !         ! print *, 'Gterm ',Gterm

    !         !% change in Junction Energy
    !         deltaEJ = (-residEJ + Csum ) / (-Gterm - Bsum)

    !         ! print *, 'deltaEJ ', deltaEJ

    !         ! print *, ' '
    !         ! print *, '========='
    !         ! print *, 'lambdaC ',jR(1:np,jr_LambdaC)
    !         ! print *, 'residQ  ',jR(1:np,jr_residQ)
    !         ! print *, 'deltaEJ ',deltaEJ
    !         ! print *, 'beta    ',jR(1:np,jr_beta)
        
    !         !% change in flowrate (recall LambdaC was overwritten with LambdaI)
    !         jR(1:np,jr_DeltaQ) = jR(1:np,jr_LambdaC) * (-jR(1:np,jr_residQ) - deltaEJ) &
    !                            / jR(1:np,jr_beta)


    !         ! print *, 'dQ ', jR(1:np,jr_DeltaQ) 
    !         ! print *, '==============='
    !         ! print *, ' '

    !         ! print *, 'old Q ',jR(1:np,jr_Q)

    !         !%  New Q flowrate
    !         jR(1:np,jr_Q) = jR(1:np,jr_Q) + jR(1:np,jr_DeltaQ)

    !         ! print *, 'new Q ',jR(1:np,jr_Q)
    !         ! print *, ' '
    !         ! print *, 'old EJ ', EJ

    !         !% --- update junction head
    !         EJ = EJ + deltaEJ

            
    !         ! print *, 'new EJ ',EJ
    !         ! print *, ' '

!     end subroutine junction_iteration   
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%    
!     subroutine junction_all_residuals &
    !             (jR, residEJ, EJ, startEJ, alpha, fcoef, Qlat, &
    !              Zcrown, OverflowHeight, OverflowType, isFull, isSubmergedOrifice, np, Npack)
    !         !%------------------------------------------------------------------
    !         !% Description:
    !         !% Computes the residuals of the nonlinear junction equations for
    !         !% all junctions
    !         !%------------------------------------------------------------------
    !             real(8), intent(inout) :: jR(:,:,:), residEJ(:)
    !             real(8), intent(in)    :: EJ(:), startEJ(:), alpha(:), Qlat(:)
    !             real(8), intent(in)    :: Zcrown(:), fcoef(:), OverflowHeight(:)
    !             logical, intent(in)    :: isSubmergedOrifice(:), isFull(:)
    !             integer, intent(in)    :: OverflowType(:), np(:), Npack 

    !             integer :: mm
    !         !%------------------------------------------------------------------

    !         !do concurrent (mm=1:Npack)
    !         do mm=1,Npack
    !             !print *, 'mm here ',mm, '---------------------;'
    !             !% --- get the flowrate residuals
    !             jR(mm,1:np(mm),jr_residQ) = junction_Qresidual                  &
    !                 (jR(mm,:,:), EJ(mm), startEJ(mm),  np(mm), mm)

    !             ! if (mm==1) then
    !             !     print *, 'eJ       ',eJ(mm)
    !             !     print *, 'startEJ  ',startEJ(mm)
    !             !     print *, 'alpha    ',alpha(mm)
    !             !     print *, 'fcoef    ',fcoef(mm)
    !             !     print *, 'Qlat     ',Qlat(mm)
    !             !     print *, 'Zcrown   ',Zcrown(mm)
    !             !     print *, 'OverflowH',OverflowHeight(mm)
    !             !     print *, 'isFull   ', isFull(mm)
    !             !     print *, 'isSubmerO', isSubmergedOrifice(mm)
    !             !     print *, 'overFType', OverflowType(mm)
    !             !     print *, 'np       ', np(mm)
    !             !     print *, ' '
    !             ! end if

    !             residEJ(mm) = junction_Eresidual                                &
    !                 (jR(mm,:,:), EJ(mm), startEJ(mm), alpha(mm), fcoef(mm),     &
    !                 Qlat(mm), Zcrown(mm), OverflowHeight(mm), isFull(mm),       &
    !                 isSubmergedOrifice(mm), OverflowType(mm), np(mm))

    !             ! if (mm==1) then 
    !             !print *, 'residEJ AA',residEJ(mm), mm
    !             ! end if
    !         end do
            
!     end subroutine junction_all_residuals
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     real(8) function junction_Eresidual &
    !             (jR, EJ, startEJ, alpha, fcoef, Qlat, Zcrown, OverflowHeight, &
    !              isFull, isSubmergedOrifice, OverflowType, np)  
    !         !%------------------------------------------------------------------
    !         !% Description:
    !         !% Computes the E residual of the nonlinear junction equation for
    !         !% a single junction
    !         !%------------------------------------------------------------------
    !             real(8), intent(in)    :: jR(:,:)
    !             real(8), intent(in)    :: EJ, startEJ, alpha, fcoef, Qlat, Zcrown
    !             real(8), intent(in)    :: OverflowHeight
    !             logical, intent(in)    :: isSubmergedOrifice, isFull
    !             integer, intent(in)    :: np, OverflowType
    !         !%------------------------------------------------------------------

    !         !% --- get the junction energy residual
    !         if (.not. isFull) then
    !             !% --- junction that is not overflowing
    !             junction_Eresidual = -alpha * (EJ - startEJ) + Qlat &
    !                                 + sum( (jR(1:np,jr_beta) * jR(1:np,jr_Q)) )  
    !         else 
    !             select case (OverflowType)

    !                 case (NoOverflow)
    !                     !% --- Zcrown becomes volume limit
    !                     junction_Eresidual = -alpha * (Zcrown - startEJ) + Qlat     &
    !                                     + sum( (jR(1:np,jr_beta) * jR(1:np,jr_Q)) )

    !                 case (OverFlowWeir)
    !                     !% --- requires additional outflow term
    !                     junction_Eresidual = -alpha * (EJ - startEJ) + Qlat            & 
    !                                     - twoR * fcoef * (EJ - Zcrown)**threehalfR     &
    !                                     + sum( (jR(1:np,jr_beta) * jR(1:np,jr_Q)) )

    !                 case (OverflowOrifice)
    !                     if (isSubmergedOrifice) then 
    !                         junction_Eresidual = -alpha * (Zcrown - startEJ)                                     &
    !                                         - twothirdR * fcoef * ((EJ -  Zcrown                  )**threehalfR) &
    !                                         + twothirdR * fcoef * ((EJ - (Zcrown + OverflowHeight))**threehalfR) &
    !                                         + sum( (jR(1:np,jr_beta) * jR(1:np,jr_Q)) )
    !                     else
    !                         junction_Eresidual = -alpha * (Zcrown - startEJ)                  &
    !                                         - twothirdR * fcoef * ((EJ - Zcrown)**threehalfR) &
    !                                         + sum( (jR(1:np,jr_beta) * jR(1:np,jr_Q)) )
    !                     end if
    !                 case (Ponded)
    !                     !% --- same residual as not full (alpha has been changed)
    !                     junction_Eresidual = -alpha * (EJ - startEJ) + Qlat          &
    !                                     + sum( (jR(1:np,jr_beta) * jR(1:np,jr_Q)) )  
                
    !             case default
    !                 !% --- change to non-pure to have error checking here
    !             end select
    !         end if

!     end function junction_Eresidual
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     function junction_Qresidual (jR, EJ, startEJ, np, mm)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% Computes the flowrate residual for all the non-zero branches
    !         !% of a single junction
    !         !%------------------------------------------------------------------
    !             real(8), intent(in) :: jR(:,:), EJ, startEJ
    !             integer, intent(in) :: np, mm
    !             real(8), dimension(np) :: junction_Qresidual
    !         !%------------------------------------------------------------------

    !         ! if (mm==1) then 
    !         !             print *, ' '
    !         ! print *, 'in junction Qresid '
    !         ! print *, 'j_a ', jR(1:np,jr_a)
    !         ! print *, 'Q   ', jR(1:np,jr_Q)
    !         ! print *, 'j_c ', jR(1:np,jr_c)
    !         ! print *, 'signedQ^2 ', sign((jR(1:np,jr_Q)**twoI),jR(1:np,jr_Q))
    !         ! print *, 'EJ      ',EJ
    !         ! print *, 'ejBranch ',jR(1:np,jr_Ebranch)
    !         ! end if
            
        
    !         junction_Qresidual(1:np) = jR(1:np,jr_a) * jR(1:np,jr_Q)        &
    !             + jR(1:np,jr_c) * sign((jR(1:np,jr_Q)**twoI),jR(1:np,jr_Q)) &
    !             + EJ - jR(1:np,jr_Ebranch)

            
    !             ! if (mm==1 )then 
    !             !      print * ,'Qresid ',junction_Qresidual(1:np)
    !             !     print *, ' '
    !             ! end if

!     end function junction_Qresidual
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     real(8) pure function junction_Gterm &
    !             (alpha, EJ, Zcrown, fcoef, OverflowHeight, OverflowType, isFull, isSubmergedOrifice)    
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% Computes G term for various overflow (or not overflow) conditions
    !         !%------------------------------------------------------------------
    !             real(8), intent(in) :: alpha, EJ, Zcrown, fcoef, OverflowHeight
    !             integer, intent(in) :: OverflowType
    !             logical, intent(in) :: isFull, isSubmergedOrifice
    !         !%------------------------------------------------------------------
    !         if (.not. isFull) then 
    !             !% --- EJ < Zcrown
    !             junction_Gterm = alpha
    !         else 
    !             !% --- overflowing conditions
    !             select case (OverflowType)

    !                 case (NoOverflow)
    !                     !% --- surcharge without overflow
    !                     junction_Gterm = zeroR

    !                 case (OverflowOrifice)
    !                     if (isSubmergedOrifice) then 
    !                         junction_Gterm = fcoef                           &
    !                             * (                                         &
    !                                 + sqrt(EJ - Zcrown)                     &
    !                                 - sqrt(EJ - (Zcrown + OverflowHeight))  &
    !                                 )
    !                     else
    !                         junction_Gterm = fcoef * sqrt(EJ - Zcrown)
    !                     end if

    !                 case (OverflowWeir)
    !                     junction_Gterm = alpha + threeR * fcoef * sqrt(EJ - Zcrown)

    !                 case (Ponded)
    !                     !% --- note alpha is from the ponded area
    !                     junction_Gterm = alpha
    !                 case default 
    !                     !% --- code should not be here
    !                     !%     to have a failure for testing 
    !                     !%     the function must be unpure    
    !             end select      
    !         end if

!     end function junction_Gterm
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     ! real(8) pure function junction_Qoverflow_residual &
    !     !         (fcoef, EJ, Zcrown, Qoverflow)    
    !     !     !%------------------------------------------------------------------
    !     !     !% Description
    !     !     !% Computes the overflow flowrate residual\
    !     !     !%------------------------------------------------------------------
    !     !         real(8), intent(in) :: EJ, Zcrown, fcoef, Qoverflow
    !     !     !%------------------------------------------------------------------

    !     !     junction_Qoverflow_residual = fcoef * ((EJ - Zcrown)**(1.5d0)) + Qoverflow

!     ! end function junction_Qoverflow_residual
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     real(8) function junction_get_L2norm (residQ, residEJ, np)
    !         !%------------------------------------------------------------------
    !         !% Description:
    !         !% Computes the L2 norm of residQ for multiple branches and the
    !         !% single EJ for one junction
    !         !%------------------------------------------------------------------
    !             real(8), intent(in) :: residQ(:), residEJ
    !             integer, intent(in) :: np
    !         !%------------------------------------------------------------------

    !         ! print *, ' '
    !         ! print *, 'residQ, residEJ ',residQ(1:np), residEJ
    !         ! print *, 'sum    ' ,sum(residQ(1:np)**twoI)
    !         ! print *, 'square ',residEJ**twoI

    !         junction_get_L2norm = sqrt(sum(residQ(1:np)**twoI) + residEJ**twoI)

    !         !stop 5098734

!     end function junction_get_L2norm
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%      
!     real(8) pure function junction_branch_Qnet (jR, np, idir)
    !         !%------------------------------------------------------------------
    !         !% Description:
    !         !% Computes the net Q in (idir = 1) or net Q out (idir = -1)
    !         !% for a junction
    !         !%------------------------------------------------------------------
    !             real(8), intent(in) :: jR(:,:)
    !             integer, intent(in) :: np, idir
    !         !%------------------------------------------------------------------

    !         junction_branch_Qnet = onehalfR * sum                                      &
    !           (                                                                 &
    !             jR(1:np,jr_beta) * jR(1:np,jr_Q)                                &
    !             * (                                                             &
    !                 oneR + real(idir,8) * jR(1:np,jr_beta) * jR(1:np,jr_Q)      & 
    !                        / abs(jR(1:np,jr_Q))                                 &   
    !               )                                                             &
    !           )

!     end function junction_branch_Qnet
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%  
!     real(8) pure function junction_Qoverflow &
    !             (EJ, Zcrown, fcoef, OverflowHeight, OverflowType, &
    !              isFull, isSubmergedOrifice)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% Computes the overflow based on theory
    !         !%------------------------------------------------------------------
    !             real(8), intent(in) :: EJ, Zcrown, fcoef, OverflowHeight
    !             integer, intent(in) :: OverflowType
    !             logical, intent(in) :: isFull, isSubmergedOrifice
    !         !%------------------------------------------------------------------

    !         !% --- baseline is zero overflow
    !         junction_Qoverflow = zeroR
    !         if (.not. isFull) return 

    !         select case (OverflowType)
    !             case (NoOverflow)
    !                 !% --- retain zero
    !             case (OverflowOrifice)
    !                 if (isSubmergedOrifice) then 
    !                     junction_Qoverflow = - twothirdR * fcoef                        &
    !                         * (                                                &
    !                             ((EJ - Zcrown)**threehalfR)                    &
    !                            -((EJ - (Zcrown + OverflowHeight))**threehalfR) &
    !                         )

    !                 else
    !                     junction_Qoverflow = - twothirdR * fcoef * ((EJ - Zcrown)**(threehalfR))
    !                 end if
    !             case (OverflowWeir)
    !                 junction_Qoverflow = -twoR * fcoef * ((EJ - Zcrown)**(threehalfR))
    !             case (Ponded)
    !                 !% --- overflow volume is held in JM volume
    !                 !%     retain zero
    !             case Default
    !                 !% --- code should not be here
    !                 !%     to have a failure for testing 
    !                 !%     the function must be unpure    
    !         end select

!     end function junction_Qoverflow
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     pure function junction_conservation (jR, residEJ, QinTotal, QoutTotal, np)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% Proportionally adjust inflow/outflow Q to ensure flowrate
    !         !% and storage residual is zero to machine accuracy
    !         !% Applies to a single junction
    !         !%------------------------------------------------------------------
    !             real(8), intent(in)     :: jR(:,:), QInTotal, QoutTotal, residEJ
    !             integer, intent(in)     :: np

    !             real(8), dimension(np)  :: junction_conservation
    !         !%------------------------------------------------------------------

    !         junction_conservation = jR(1:np,jr_Q) &
    !             * (                                                                                  &
    !                oneR - onefourthR * residEJ                                                       &
    !                * (                                                                               &
    !                     (oneR + (jR(1:np,jr_beta) * jR(1:np,jr_Q)) / abs(jR(1:np,jr_Q)) ) / QinTotal  &
    !                    +                                                                             &
    !                     (oneR - (jR(1:np,jr_beta) * jR(1:np,jr_Q)) / abs(jR(1:np,jr_Q)) ) / QoutTotal &
    !                  )                                                                               &
    !                )
                
!     end function junction_conservation
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     pure subroutine junction_push_face_flowrate &
    !             (faceQ, faceVup, faceVdn, faceAup, faceAdn, jR, jI,  np, JMidx)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% Cycles through the solved branches (1:np) and restores their
    !         !% Q data to the faceR array 
    !         !%------------------------------------------------------------------
    !         !% Declarations
    !             real(8), intent(inout) :: faceQ(:),   faceVup(:), faceVdn(:)
    !             real(8), intent(in)    :: faceAup(:), faceAdn(:)
    !             real(8), intent(in)    :: jR(:,:)
    !             integer, intent(in)    :: np, JMidx
    !             integer, target, intent(in) ::  jI(:,:)

    !             integer :: mm, ii
    !             !integer, pointer :: tM, kidx, fIdx
    !             integer :: kidx(np), fIdx(np)
    !         !%------------------------------------------------------------------
    !         !%------------------------------------------------------------------

    !         !% --- cycle through the branches that were solved
    !         do mm=1,np
                
    !             !% --- branch kk in the 1:max_branch_per_node
    !             kidx(mm) = jI(mm,ji_kidx)
                
    !             !% --- handle upstream or downstream separately
    !             if (mod(kidx(mm)+1,twoI) .eq. zeroI) then 
    !                 !% --- upstream branch face
    !                 fIdx(mm) = elemI(JMidx+kidx(mm),ei_Mface_uL)
    !             else
    !                 !% --- downstream branch face
    !                 fIdx(mm) = elemI(JMidx+kidx(mm),ei_Mface_dL)
    !             end if
    !             !% --- set the flowrate
    !             faceQ(fIdx(mm))   = jR(mm,jr_Q)
    !             !% --- get consistent velocities (assumes area > epsilon)
    !             faceVup(fIdx(mm)) = faceQ(fIdx(mm)) / faceAup(fIdx(mm))
    !             faceVdn(fIdx(mm)) = faceQ(fIdx(mm)) / faceAdn(fIdx(mm))
    !         end do

!     end subroutine junction_push_face_flowrate   
! !%
! !%========================================================================== 
! !%========================================================================== 
! !% 
!     real(8) pure function junction_volume &
    !             (oldVolume, Qin, Qout, Qlat, Qoverflow, istep)
    !         !%------------------------------------------------------------------
    !         !% Description
    !         !% Updates junction volume from Volume_N0 for an RK step
    !         !%------------------------------------------------------------------
    !         !% Declarations:
    !             real(8), intent(in) :: oldVolume, Qin, Qout, Qlat, Qoverflow
    !             integer, intent(in) :: istep
    !         !%------------------------------------------------------------------

    !         junction_volume = oldVolume                    &
    !             + setting%Solver%crk2(istep) * setting%Time%Hydraulics%Dt  &
    !             * (Qin - Qout + Qlat + Qoverflow)

!     end function junction_volume 
! !%
! !%========================================================================== 
! !%========================================================================== 
! !% 
!     subroutine junction_check_conservation (thisJM, Qoverflow, jR, Npack, istep)
    !         !%------------------------------------------------------------------
    !         !% Description:
    !         !% Checks flowrate and junction storage conservation
    !         !%------------------------------------------------------------------
    !         !% Declarations
    !             integer, target, intent(in) :: thisJM(:)
    !             integer, intent(in) :: Npack, istep
    !             real(8), intent(in) :: Qoverflow(:), jR(:,:,:)
    !             integer, pointer :: tM, upF, dnF
    !             integer          :: mm, kk
    !             real(8)          :: thisCons, thisVol, thisFlow
    !             character(64) :: subroutine_name = "junction_check_conservation"
    !         !%------------------------------------------------------------------
    !             ! print *, ' '
    !             ! print *, 'in ',trim(subroutine_name)

    !         do mm=1,Npack
    !             tM => thisJM(mm)

    !             ! print *, 'junction index ',tM
    !             ! print *, 'Qlateral ',elemR(tM,er_FlowrateLateral)
    !             ! print *, 'Qoverflow',Qoverflow(mm)
                
    !             !% --- lateral inflow
    !             thisCons = elemR(tM,er_FlowrateLateral) + Qoverflow(mm)
    !             thisFlow = abs(elemR(tM,er_FlowrateLateral))
                
    !             ! print *, 'thisCons 1',thisCons

    !             !% --- change in storage
    !             !%     assumes old head stored in er_Temp01
    !             select case (elemSI(tM,esi_JunctionMain_Type))
    !                 case (ImpliedStorage)
    !                     ! print *, 'in implied storage'
    !                     thisVol = sum(jR(mm,:,jr_Length) * jr(mm,:,jr_Area) * abs(jR(mm,:,jr_beta))) / twoR
    !                     !% --- no change to conservation by storage
    !                 case (TabularStorage,FunctionalStorage)
    !                     thisCons = thisCons &
    !                          - (elemSR(tM,esr_Storage_Plan_Area)    &
    !                               / (  setting%Solver%crk2(istep)    &
    !                                   * setting%Time%Hydraulics%Dt ) &
    !                            )                                     &
    !                             * (elemR(tM,er_Head) - elemR(tM,er_Temp01))
    !                     thisVol = elemR(tM,er_Volume)
    !                 case default
    !                     print *, 'CODE ERROR: unexpected case default'
    !             end select

    !             ! print *, 'thisCons 2',thisCons
    !             ! print *, 'thisVol  2',thisVol

    !             !% --- upstream branches
    !             do kk=1,max_branch_per_node,2
    !                 if (elemSI(tM+kk,esi_JunctionBranch_Exists) == oneI) then
    !                     upF => elemI(tM+kk,ei_Mface_uL)
    !                     thisCons = thisCons + faceR(upF,fr_Flowrate)
    !                     !print *, 'thisCons 3',thisCons
    !                     thisFlow = thisFlow + abs(faceR(upF,fr_Flowrate))
    !                 end if
    !             end do

    !             !% --- downstream branches
    !             do kk=2,max_branch_per_node,2
    !                 if (elemSI(tM+kk,esi_JunctionBranch_Exists) == oneI) then
    !                     dnF => elemI(tM+kk,ei_Mface_dL)
    !                     thisCons = thisCons - faceR(dnF,fr_Flowrate)
    !                     !print *, 'thisCons 4',thisCons
    !                     thisFlow = thisFlow + abs(faceR(dnF,fr_Flowrate))
    !                 end if
    !             end do

    !             ! print *, 'thisCons ',thisCons
    !             ! print *, 'thisFlow ',thisFlow
    !             ! print *, 'thisVol  ',thisVol
    !             ! print *, 'dt       ',setting%Time%Hydraulics%Dt

    !             !% volume based normalization
    !             if (thisVol > setting%ZeroValue%Volume) then
    !                 thisCons = thisCons * setting%Time%Hydraulics%Dt / thisVol
    !             else
    !                 !% flowrate based normalization
    !                 if (thisFlow > setting%ZeroValue%Volume) then
    !                     thisCons = thisCons / thisFlow
    !                 else 
    !                     !% no normalization
    !                 end if
    !             end if

    !             if (abs(thisCons) > 1.d-2) then 
    !                 print *, 'conservation error at ',thisJM(mm)
    !                 print *, 'magnitude ',thisCons
    !                 stop 5509873
    !             end if
    !         end do
            
!     end subroutine junction_check_conservation
! !%
! !%========================================================================== 
! !% JUNCTION PACKING BELOW
! !%========================================================================== 
! !%
!     pure subroutine junction_check_repack_all (isRepack, jR, np, Npack, jEpsilon)
!         !%------------------------------------------------------------------
!         !% Description
!         !% Checks if any branches are now zero and need repacking and 
!         !% adjusts beta to zero where branches are zero
!         !%------------------------------------------------------------------
!         !% Declarations
!             logical, intent(inout) :: isRepack(:)
!             real(8), intent(inout) :: jR(:,:,:)
!             integer, intent(in)    :: np(:), Npack
!             real(8), intent(in)    :: jEpsilon

!             integer :: mm
!         !%------------------------------------------------------------------
!         !%------------------------------------------------------------------
!         do concurrent (mm=1:Npack)

!             call junction_check_repack_oneJ (jR(mm,:,:), np(mm), jEpsilon, isRepack(mm))

!         end do

!     end subroutine junction_check_repack_all   
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%    
!     pure subroutine junction_check_repack_oneJ (jR, np, jEpsilon, isRepack)
!         !%------------------------------------------------------------------
!         !% Checks for zeros and resets beta if repack needed
!         !%------------------------------------------------------------------
!         !% Declarations
!             real(8), intent(inout) :: jR(:,:)
!             real(8), intent(in)    :: jEpsilon
!             integer, intent(in)    :: np 
!             logical, intent(inout) :: isRepack
!         !%------------------------------------------------------------------
!         !%------------------------------------------------------------------

!         if (any(jR(1:np,jr_LambdaC) < jEpsilon)) then 
!             isRepack = .true.
!             where (jR(1:np,jr_LambdaC) < jEpsilon)
!                 jR(1:np,jr_beta) = zeroR
!             endwhere
!         else
!             isRepack = .false.
!         end if

!     end subroutine junction_check_repack_oneJ
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%
!     pure subroutine junction_pack_all &
!         (jR, jI, np, isRepack, Npack)
!         !%------------------------------------------------------------------
!         !% Description
!         !% Packing on the 2nd dimension of 3D arrays, jR, jI
!         !% Cycles through the 1st dimension (packed JM index) in this procedure 
!         !% and then cycles thorugh through the 3rd dimension (data type)
!         !% in call to subsidiary procedures
!         !%------------------------------------------------------------------
!         !% Declarations
!             integer, intent(inout) :: np(:), jI(:,:,:)
!             real(8), intent(inout) :: jR(:,:,:)
!             logical, intent(inout) :: isRepack(:)
!             integer, intent(in)    :: Npack
        
!             integer :: mm
!         !%------------------------------------------------------------------
!         !% --- pack the jDataR and jDataI to remove zero branches
!         do concurrent (mm=1:Npack)
!             call junction_pack_oneJ(jR(mm,:,:), jI(mm,:,:), np(mm), isRepack(mm))
!         end do

!     end subroutine junction_pack_all
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%   
!     pure subroutine junction_pack_oneJ (jR, jI, np, isRepack)
!         !%------------------------------------------------------------------
!         !% Description
!         !% Packs data for one junction during implicit junction solution
!         !%------------------------------------------------------------------
!         !% Declarations
!             real(8), intent(inout) :: jR(:,:)
!             integer, intent(inout) :: jI(:,:), np
!             logical, intent(inout) :: isRepack
!         !%------------------------------------------------------------------
!         !%------------------------------------------------------------------
!         if (.not. isRepack) return 
    
!         !% --- count the nonzero branches
!         np = count(jR(:,jr_beta) .ne. zeroR)

!         !% --- perform the packs based on jr_beta (must delay packing jr_beta)
!         call junction_pack_all_jRtype (jR(:,:), np)
!         call junction_pack_all_jItype (jI(:,:), jR(:,:), np)

!         !% --- pack beta (always must be last)
!         jR(1:np,jr_beta)  = junction_pack_one_arrayR &
!                 (jR(:,jr_beta), jR(:,jr_beta), np)

!         !% --- repack finished, reset the logical
!         isRepack = .false.

!     end subroutine junction_pack_oneJ
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%        
!     pure subroutine junction_pack_all_jRtype (jR, np)  
!         !%------------------------------------------------------------------
!         !% Description:
!         !% Packs the first dimension of a 2D real array (jR) by cycling
!         !% through the data types in the 2nd dimension.
!         !% Note that np must equal the number of jR(:,jr_beta) that are 
!         !% nonzero
!         !%------------------------------------------------------------------
!         !% Declarations:
!             real(8), intent(inout) :: jR(:,:)
!             integer, intent(in)    :: np

!             !% --- jData elements that are packed
!             integer, parameter :: NjRset = 13
!             integer :: jRset(NjRset)
!             integer :: jj
!         !%------------------------------------------------------------------

!         !% --- define data sets used in packing 
!         !%     DO NOT INCLUDE jr_beta -- must be repacked last in a separate 
!         !%     call since it is themask used for packing.
!         jRset = [jr_Area, jr_a, jr_c, jr_Ebranch, jr_Gamma, jr_K, jr_Length, jr_Q, &
!                  jr_DeltaQ, jr_residQ, jr_LambdaA, jr_LambdaB, jr_LambdaC]

!         do concurrent (jj=1:NjRset)
!             !% --- pack the jRset(jj) data type for the thisJM(mm) branch
!             !%     to remove non-existent branches
!             jR(1:np,jRset(jj))  = junction_pack_one_arrayR                               &
!                         (jR(:,jRset(jj)), jR(:,jr_beta), np)
!         end do

!     end subroutine junction_pack_all_jRtype
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%        
!     pure subroutine junction_pack_all_jItype (jI, jR, np)  
!         !%------------------------------------------------------------------
!         !% Description:
!         !% Packs the first dimension of a 2D integer array (jI) by cycling
!         !% through the data types in the 2nd dimension.
!         !% Note that np must equal the number of jR(:,jr_beta) that are 
!         !% nonzero
!         !%------------------------------------------------------------------
!             integer, intent(inout) :: jI(:,:)
!             real(8), intent(in)    :: jR(:,:)
!             integer, intent(in)    :: np

!             integer, parameter :: NjIset = Ncol_jDataI
!             integer            :: jIset(NjIset)
!             integer            :: jj
!         !%------------------------------------------------------------------
!         !%------------------------------------------------------------------

!         !% --- define data sets used in packing
!         jIset = [ji_kidx]

!         do concurrent (jj=1:NjIset)
!             !% --- pack the jIset(jj) data type for the thisJM(mm) branch
!             !%     to remove non-existent branches
!             jI(1:np,jIset(jj))  = junction_pack_one_arrayI                               &
!                         (jI(:,jIset(jj)), jR(:,jr_beta), np)
!         end do

!     end subroutine junction_pack_all_jItype
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%      
!     pure function junction_pack_one_arrayR (inarray, beta, np)
!         !%------------------------------------------------------------------
!         !% Description
!         !% standard pack to remove any locations where beta = 0.0
!         !%------------------------------------------------------------------
!             real(8), intent(in)    :: inarray(:), beta(:)
!             integer, intent(in)    :: np
!             real(8), dimension(np) :: junction_pack_one_arrayR
!         !%------------------------------------------------------------------

!         junction_pack_one_arrayR(1:np) = pack (inarray, beta .ne. zeroR)

!     end function junction_pack_one_arrayR
! !%
! !%==========================================================================  
! !%========================================================================== 
! !%      
!     pure function junction_pack_one_arrayI (inarray, beta, np)
!         !%------------------------------------------------------------------
!         !% Description
!         !% standard pack to remove any locations where beta = 0.0
!         !%------------------------------------------------------------------
!             real(8), intent(in)    ::  beta(:)
!             integer, intent(in)    :: np, inarray(:)
!             integer, dimension(np) :: junction_pack_one_arrayI
!         !%------------------------------------------------------------------

!         junction_pack_one_arrayI(1:np) = pack (inarray, beta .ne. zeroR)

!     end function junction_pack_one_arrayI
! !%
! !%==========================================================================  
!     !%========================================================================== 
! !%
!     function matinv(A) result (invA)
!         !%-----------------------------------------------------------------
!         !% Description:
!         !% Find the inverse of an input matrix A
!         !%-----------------------------------------------------------------
!         !% Declarations
!         real(8),allocatable, intent(in) :: A(:,:)
!         real(8),allocatable             :: invA(:,:)
!         integer                         :: n   
!         !%-----------------------------------------------------------------
!         !% find the size of the matrix
!         n = size(A,1)

!         !% for square matrix size or 2-4, direct inversion is the fastest 
!         if (n == 2) then
!             invA = matinv2(A)
!         else if (n == 3) then
!             invA = matinv3(A)
!         else if (n == 4) then
!             invA = matinv4(A)
!         else
!         !% use the LAPACK for higher order matrices
!         !% talk with cesar about LAPACK installation recipie
!             print*, 'Matrix inversion for shape > 4 has not yet been developed'
!             stop 789451
!         end if     

!     end function 
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%
!     pure function matinv2(A) result(B)
!         !! Performs a direct calculation of the inverse of a 2×2 matrix.
!         real(8), intent(in) :: A(2,2)   !! Matrix
!         real(8)             :: B(2,2)   !! Inverse matrix
!         real(8)             :: detinv

!         ! Calculate the inverse determinant of the matrix
!         detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

!         ! Calculate the inverse of the matrix
!         B(1,1) = +detinv * A(2,2)
!         B(2,1) = -detinv * A(2,1)
!         B(1,2) = -detinv * A(1,2)
!         B(2,2) = +detinv * A(1,1)
!     end function
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%
!     pure function matinv3(A) result(B)
!         !! Performs a direct calculation of the inverse of a 3×3 matrix.
!         real(8), intent(in) :: A(3,3)   !! Matrix
!         real(8)             :: B(3,3)   !! Inverse matrix
!         real(8)             :: detinv

!         ! Calculate the inverse determinant of the matrix
!         detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
!                 - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
!                 + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

!         ! Calculate the inverse of the matrix
!         B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
!         B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
!         B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
!         B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
!         B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
!         B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
!         B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
!         B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
!         B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
!     end function
! !%
! !%========================================================================== 
! !%========================================================================== 
! !%
!     pure function matinv4(A) result(B)
!         !! Performs a direct calculation of the inverse of a 4×4 matrix.
!         real(8), intent(in) :: A(4,4)   !! Matrix
!         real(8)             :: B(4,4)   !! Inverse matrix
!         real(8)             :: detinv

!         ! Calculate the inverse determinant of the matrix
!         detinv = &
!         1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
!         - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
!         + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
!         - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

!         ! Calculate the inverse of the matrix
!         B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
!         B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
!         B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
!         B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
!         B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
!         B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
!         B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
!         B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
!         B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
!         B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
!         B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
!         B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
!         B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
!         B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
!         B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
!         B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
!     end function 
!%
!%========================================================================== 
!% END OF MODULE
!%==========================================================================
end module junction_elements