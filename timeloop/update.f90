module update

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use geometry
    use adjust

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Updates values during timeloop of hydraulics.
    !%

    private
    
    public :: update_auxiliary_variables
    public :: update_Froude_number_junction_branch

    real(8), pointer :: grav => setting%constant%gravity

    contains
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    !%
    subroutine update_auxiliary_variables (whichTM)
        !%-----------------------------------------------------------------------------
        !% Description: 
        !% 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: whichTM  !% indicates which Time marching method
        integer, pointer :: thisCol_all
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'update_auxiliary_variables'
        if (setting%Debug%File%update) print *, '*** enter ', this_image(), subroutine_name     
        !%-----------------------------------------------------------------------------
        !%  
        !% update the head (non-surcharged) and geometry

        !print *, '---- in ',subroutine_name,'   y01'
        !write(*,'(7F9.4,A15)') elemR(ietmp,er_Head),' Head elem '

        call geometry_toplevel (whichTM)

        !print *, '---- in ',subroutine_name,'   y02'
        !write(*,'(7F9.4,A15)') elemR(ietmp,er_Head),' Head elem '

        !% adjust velocity with limiters and small volume treatment
        call adjust_velocity (whichTM, er_Velocity, er_Volume)

        !print *, '---- in ',subroutine_name,'   y03'
        !write(*,'(7F9.4,A15)') elemR(ietmp,er_Head),' Head elem '

        !% set packed column for updated elements
        select case (whichTM)
            case (ALLtm)
                thisCol_all => col_elemP(ep_CC_ALLtm)
            case (ETM)
                thisCol_all => col_elemP(ep_CC_ETM)
            case (AC)
                thisCol_all => col_elemP(ep_CC_AC)
            case default
                print *, 'error, default case should not be reached'
                stop 7489
        end select

        !% Compute the flowrate on CC.
        !% Note that JM should have 0 flowrate and JB has lagged flowrate at this point.
        !% The JB flowrate is not updated until after face interpolation
        call update_CC_element_flowrate (thisCol_all)

        !print *, '---- in ',subroutine_name,'   y04'
        !write(*,'(7F9.4,A15)') elemR(ietmp,er_Head),' Head elem '

        !% compute element Froude numbers for CC, JM
        call update_Froude_number_element (thisCol_all)

        !print *, '---- in ',subroutine_name,'   y05'
        !write(*,'(7F9.4,A15)') elemR(ietmp,er_Head),' Head elem '

        !print *, '---- in ',subroutine_name,'   y07'
        !write(*,'(7F9.4,A15)') elemR(ietmp,er_Head),' Head elem '

        !% compute element face interpolation weights on CC, JM
        call update_interpolation_weights_element (thisCol_all, whichTM)

        !print *, '---- in ',subroutine_name,'   y06'
        !write(*,'(7F9.4,A15)') elemR(ietmp,er_Head),' Head elem '
   
        if (setting%Debug%File%update)  print *, '*** leave ', this_image(), subroutine_name
    end subroutine update_auxiliary_variables
    !%
    !%==========================================================================
    !% PRIVATE
    !%==========================================================================   
    !%  
    subroutine update_CC_element_flowrate (thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisCol
        !%-----------------------------------------------------------------------------
        integer, pointer ::  Npack, thisP(:)
        real(8), pointer :: flowrate(:), velocity(:), area(:)
        !%-----------------------------------------------------------------------------   
        flowrate => elemR(:,er_Flowrate)
        velocity => elemR(:,er_Velocity)
        area     => elemR(:,er_Area)   
        !%-----------------------------------------------------------------------------    
        Npack => npack_elemP(thisCol)
        if (Npack > 0) then
            thisP    => elemP(1:Npack,thisCol)
            flowrate(thisP) = area(thisP) * velocity(thisP)
        endif

    end subroutine update_CC_element_flowrate
    !%
    !%==========================================================================   
    !%==========================================================================   
    !%  
    subroutine update_Froude_number_element (thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes Froude number on each element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisCol
        integer, pointer :: Npack, thisP(:)
        real(8), pointer :: Froude(:), velocity(:), depth(:)
        !%-----------------------------------------------------------------------------
        Froude   => elemR(:,er_FroudeNumber)
        velocity => elemR(:,er_Velocity)
        depth    => elemR(:,er_ell)  !% Use the ell value (modified hydraulic depth)
        !%-----------------------------------------------------------------------------
    
        Npack => npack_elemP(thisCol)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisCol)
            Froude(thisP) = velocity(thisP) / sqrt(grav * depth(thisP))
        endif
    
    end subroutine update_Froude_number_element
    !%   
    !%==========================================================================   
    !%==========================================================================   
    !%  
    subroutine update_Froude_number_junction_branch (thisCol_JM)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes Froude number on each junction branch element
        !% BRHbugfix 20210812
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'update_Froude_number_junction_branch'
        integer, intent(in) :: thisCol_JM
        integer, pointer :: Npack, thisP(:), tM, BranchExists(:)
        real(8), pointer :: Froude(:), velocity(:), depth(:)
        integer :: ii, kk, tB
        !%-----------------------------------------------------------------------------
        Froude   => elemR(:,er_FroudeNumber)
        velocity => elemR(:,er_Velocity)
        depth    => elemR(:,er_ell)  !% Use the ell value (modified hydraulic depth)
        BranchExists => elemSI(:,eSI_JunctionBranch_Exists)
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%update) print *, '*** enter ', this_image(), subroutine_name
    
        Npack => npack_elemP(thisCol_JM)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisCol_JM)
            do ii=1,Npack
                tM => thisP(ii)
                do kk=1,max_branch_per_node
                    tB = tM + kk
                    if (BranchExists(tB)==1) then
                        Froude(tB) = velocity(tB) / sqrt(grav * depth(tB))
                        !print *, kk, tB, Froude(tB), velocity(tB),'  Froude JB'
                    end if
                end do
            end do
        end if

        if (setting%Debug%File%update)  print *, '*** leave ', this_image(), subroutine_name
    end subroutine update_Froude_number_junction_branch
    !%   
    !%==========================================================================   
    !%==========================================================================   
    !%  
    subroutine update_interpolation_weights_element (thisCol, whichTM)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes the interpolation weights on each element form CC, JM
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'update_interpolation_weights_element'
        integer, intent(in) :: thisCol, whichTM
        integer, pointer :: Npack, Npack2, thisCol_AC,  thisP(:), thisP2(:)
        real(8), pointer :: velocity(:), wavespeed(:), depth(:), length(:)
        real(8), pointer :: w_uQ(:), w_dQ(:),  w_uG(:), w_dG(:),  w_uH(:), w_dH(:)
        real(8), pointer :: Fr(:) !BRHbugfix20210811 test
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%update) print *, '*** enter ', this_image(), subroutine_name

        velocity  => elemR(:,er_Velocity)
        wavespeed => elemR(:,er_WaveSpeed)
        depth     => elemR(:,er_ell)  !% modified hydraulic depth!
        length    => elemR(:,er_Length)
        w_uQ      => elemR(:,er_InterpWeight_uQ)
        w_dQ      => elemR(:,er_InterpWeight_dQ)
        w_uG      => elemR(:,er_InterpWeight_uG)
        w_dG      => elemR(:,er_InterpWeight_dG)
        w_uH      => elemR(:,er_InterpWeight_uH)
        w_dH      => elemR(:,er_InterpWeight_dH)   
        
        Fr        => elemR(:,er_FroudeNumber)  !BRHbugfix20210811 test
        !%-----------------------------------------------------------------------------
        !% 2nd cases needed for handling surcharged AC elements and using the celerity
        !% multiplier of the AC method for the wavespeed
        select case (whichTM)
            case (ALLtm)
                thisCol_AC =>  col_elemP(ep_Surcharged_AC)
            case (ETM)
                !% no effect
            case (AC)
                thisCol_AC =>  col_elemP(ep_Surcharged_AC)
            case default
                print *, 'error, case default should not be reached.'
                stop 3987
        end select    
    
        Npack => npack_elemP(thisCol)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisCol)

            !% wavespeed at modified hydraulic depth (ell)    
            wavespeed(thisP) = sqrt(grav * depth(thisP))
        
            !% modify wavespeed for surcharged AC cells
            if (whichTM .ne. ETM) then
                Npack2 => npack_elemP(thisCol_AC)
                if (Npack2 > 0) then
                    thisP2 => elemP(1:Npack2,thisCol_AC)
                    wavespeed(thisP2) = wavespeed(thisP2) * setting%ACmethod%Celerity%RC
                endif    
            endif
    

            !% timescale interpolation weights for flowrate
            !% Modified from original approach by Froude number weighting
            !% Note that Fr is +/- depending on flow direction, so if the Fr is an odd power
            !% it needs to have an abs() e.g, abs(Fr(thisp)**3) *
            w_uQ(thisP) = - onehalfR * length(thisP)  / ( abs(Fr(thisp)**10) * velocity(thisP) - wavespeed(thisP)) !BRHbugfix 20210813 testing Fr
            w_dQ(thisP) = + onehalfR * length(thisP)  / ( abs(Fr(thisp)**10) * velocity(thisP) + wavespeed(thisP)) !BRHbugfix 20210813 testing Fr
            
            !% apply limiters to timescales
            where (w_uQ(thisP) < zeroR)
                w_uQ(thisP) = setting%Limiter%InterpWeight%Maximum
            endwhere
            where (w_uQ(thisP) < setting%Limiter%InterpWeight%Minimum)    
                w_uQ(thisP) = setting%Limiter%InterpWeight%Minimum
            endwhere
            where (w_uQ(thisP) > setting%Limiter%InterpWeight%Maximum)    
                w_uQ(thisP) = setting%Limiter%InterpWeight%Maximum
            endwhere    
            
            !% timescale interpolation for geometry are identical to flowrate
            !% but may be modified elsewhere
            w_uG(thisP) = w_uQ(thisP)
            w_dG(thisP) = w_dQ(thisP)
            
            !% head uses length scale interpolation
            !% This shouldn't need limiters.
            
            w_uH(thisP) = onehalfR * length(thisP)
            w_dH(thisP) = onehalfR * length(thisP)
         
        endif

        if (setting%FaceInterp%DownJBFaceInterp == dynamic) then
            if (num_images() > oneI) then
                print*, 'error: dynamic face interpolation for ds JB does not support multiple processors yet'
            else
                !% testin a new branch interp technique
                call update_interpolation_weights_ds_JB ()
            endif
        endif

        !print *
        !print *,'--- in ',trim(subroutine_name),' ----------------------------------------- end'
        !write(*,'(7e11.4,A15)') elemR(ietmp,er_InterpWeight_dQ),' InterpWeight_dQ'
        !write(*,'(7e11.4,A15)') elemR(ietmp,er_InterpWeight_uQ),' InterpWeight_uQ'
        ! print *, elemR(ietmp(1), er_InterpWeight_dQ)
        ! print *, elemR(ietmp(2), er_InterpWeight_dQ)
        ! print *, elemR(ietmp(3), er_InterpWeight_dQ)
        ! print *, elemR(ietmp(4), er_InterpWeight_dQ)
        ! print *, elemR(ietmp(5), er_InterpWeight_dQ)
        ! print *, elemR(ietmp(6), er_InterpWeight_dQ)
        ! print *, elemR(ietmp(7), er_InterpWeight_dQ)

        if (setting%Debug%File%update)  print *, '*** leave ', this_image(), subroutine_name
    end subroutine update_interpolation_weights_element
    !%   
    !%==========================================================================   
    !%==========================================================================   
    !%  
    subroutine update_interpolation_weights_ds_JB ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% This is a test subroutine that violates no neighbour algorithm
        !% This subroutine sets the interpolation wights in ds JB to its 
        !% conneceted link element
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'update_interpolation_weights_ds_JB'
        integer, pointer :: thisColP_JM, fUp(:), fDn(:), eUp(:), eDn(:), tM
        integer, pointer :: Npack, Npack2, thisCol_AC,  thisP(:), thisP2(:), BranchExists(:)
        real(8), pointer :: w_uQ(:), w_dQ(:),  w_uG(:), w_dG(:),  w_uH(:), w_dH(:)
        integer :: ii, kk, tB
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%update)  print *, '*** enter ', subroutine_name
        w_uQ      => elemR(:,er_InterpWeight_uQ)
        w_dQ      => elemR(:,er_InterpWeight_dQ)
        w_uG      => elemR(:,er_InterpWeight_uG)
        w_dG      => elemR(:,er_InterpWeight_dG)
        w_uH      => elemR(:,er_InterpWeight_uH)
        w_dH      => elemR(:,er_InterpWeight_dH)
        fUp       => elemI(:,ei_Mface_uL)    
        fDn       => elemI(:,ei_Mface_dL) 
        eUp       => faceI(:,fi_Melem_uL) 
        eDn       => faceI(:,fi_Melem_dL)
        BranchExists => elemSI(:,eSI_JunctionBranch_Exists)
        !%-----------------------------------------------------------------------------

        thisColP_JM  => col_elemP(ep_JM_ALLtm)
        Npack        => npack_elemP(thisColP_JM)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisColP_JM)
            do ii=1,Npack
                tM => thisP(ii) !% junction main ID
                !% only execute for whichTM of ALL or thisSolve (of JM) matching input whichTM
                !% setting the interp weight of ds JB same as its ds link element
                !% handle the downstram branches
                do kk=2,max_branch_per_node,2
                    tB = tM + kk
                    if (BranchExists(tB)==1) then
                        !% Baseline is all commented
                        !% case 1  Q of downstream JB equal with Q upstream of next element down
                        w_dQ(tB) = w_uQ(eDn(fDn(tB)))
                        ! w_uQ(tB) = w_dQ(eUp(fUp(tB)))

                        !% case 2  G of downstream JB equal with G upstream of next element down
                        w_dG(tB) = w_uG(eDn(fDn(tB)))
                        ! w_uG(tB) = w_dG(eUp(fUp(tB)))

                        !% case 3 H of downstream JB equal with H of upstream of next element down  
                        w_dH(tB) = w_uH(eDn(fDn(tB)))
                        ! w_uH(tB) = w_dH(eUp(fUp(tB)))


                        !% case 4 Q,G,H all changed 
                         !w_dQ(tB) = w_uQ(eDn(fDn(tB)))
                         !w_dG(tB) = w_uG(eDn(fDn(tB)))
                         !w_dH(tB) = w_uH(eDn(fDn(tB)))

                        !  print *,'---'
                        !  print *, tM,'tM'
                        !  print *, fUp(tB) ,'fUp(tB)'
                        !  print *, tb ,'tB'
                        !  print *, fDn(tB) ,'fDn(tB)'
                        !  print *, eDn(fDn(tB)), 'eDn(fDn(tB))'
                        !  print *, fDn(eDn(fDn(tB))), 'fDn(eDn(fDn(tB)))'
                        !  print *, eUp(fDn(tB)), 'eUp(fDn(tB))'
                        !  print *, w_uQ(eDn(fDn(tB))), 'w_uQ(eDn(fDn(tB)))'
                        !  print *, w_dQ(tB), 'w_dQ(tB)'
                        !  print *, w_uG(eDn(fDn(tB))), 'w_uG(eDn(fDn(tB)))'
                        !  print *, w_dG(tB), 'w_dG(tB)'
                        !  print *, w_uH(eDn(fDn(tB))), 'w_uH(eDn(fDn(tB)))'
                        !  print *, w_dH(tB), 'w_dH(tB)'                        
                        !  print *,'-----'

                        !w_uQ(tB) = w_uQ(eDn(fDn(tB)))
                        !w_dQ(tB) = w_dQ(eDn(fDn(tB)))
                        !w_uG(tB) = w_uG(eDn(fDn(tB)))
                        !w_dG(tB) = w_dG(eDn(fDn(tB)))
                        !w_uH(tB) = w_uH(eDn(fDn(tB)))
                        !w_dH(tB) = w_dH(eDn(fDn(tB)))
                    end if
                end do
            end do
        endif

        if (setting%Debug%File%update)  print *, '*** leave ', subroutine_name
    end subroutine update_interpolation_weights_ds_JB
    !% 
    !%==========================================================================
    !% END OF MODULE
    !%+=========================================================================
end module update