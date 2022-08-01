module control_hydraulics
    use define_keys
    use define_globals
    use define_settings
    use define_indexes
    use interface_
    use orifice_elements
    use weir_elements
    use pump_elements
    use outlet_elements
    use adjust,           only: adjust_face_for_zero_setting_singular
    use utility_allocate, only: util_allocate_monitor_points, util_allocate_action_points
    use utility,          only: util_unique_rank
    use utility_crash,    only: util_crashpoint

    implicit none
    !%-----------------------------------------------------------------------------
    !% Description:
    !%    Handles controls/monitoring algorithms for hydraulics based
    !%    on EPA-SWMM controls.c code
    !%-----------------------------------------------------------------------------

    private

    public :: control_update
    public :: control_init_monitoring_and_action_from_EPASWMM

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine control_update ()
        !%------------------------------------------------------------------
        !% Description:
        !% Implements the control functionality in EPA-SWMM.
        !% This extracts the monitoring data across the elements of SWMM5+
        !% and stores these in the link/node system of EPA-SWMM that has
        !% previously been initialized on each image. The EPA-SWMM code
        !% for parsing the controls and adjusting the Link[].settings is
        !% then called through the API interface. Once the controls have
        !% been evaluated, the action changes are retrieved from EPA-SWMM
        !% back into SWMM5+ and the elem(:,er_Setting) and 
        !% elem(:,er_TargetSetting) values are updated.
        !%------------------------------------------------------------------
        !% Declarations:
            character(64) :: subroutine_name = 'control_update'
        !%------------------------------------------------------------------

        !% --- store monitoring element data across all images
        call control_update_monitor_across_images()

        !% --- transfer monitoring data to EPA-SWMM to all images
        call control_update_EPASWMM_monitor_data()

        !% --- execute control algorithm in EPA-SWMM on all images
        call interface_controls_execute ()

        !% --- retrieve action changes elemR(:,er_TargetSettings) and
        !%     elemR(:,er_TimeLastSet)
        call control_update_actions ()

        !% --- adjust settings
        !%     here we need to change the elemR(:,er_settings) for orifices, 
        !%     weirs, pumps, outlets and conduits.
        call control_update_setting ()

        !% --- adjust the values in the elemR array for new elemR(:,er_Setting)
        call control_update_element_values ()

        !%------------------------------------------------------------------

    end subroutine control_update
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine control_init_monitoring_and_action_from_EPASWMM()
        !%------------------------------------------------------------------
        !% Description
        !% initializes the storage for controls and monitoring from the
        !% EPASWMM code
        !%------------------------------------------------------------------
            integer :: nRules, nPremise, nThenAction, nElseAction, ii, rr
            integer :: iLeft, iRight, thisPremiseLevel, success, numUnique
            integer :: npoint
            integer, allocatable :: location(:), linknodesimType(:), attribute(:)
            integer, allocatable :: irank(:)
            character(64) :: subroutine_name = 'control_init_monitoring_and_action_from_EPASWMM'
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        !% --- get the number of locations in link, node arrays that are
        !%     involved with control and monitoring
        call interface_controls_count(nRules, nPremise, nThenAction, nElseAction)

            ! print *, 'in ',trim(subroutine_name)
            ! print *, '       nRules, nPremise, nThenAction, nElseAction'
            ! print *, nRules, nPremise, nThenAction, nElseAction

        !% --- allocate and initialize the monitorI(:,:) array from control premise data
        call control_init_monitor_points(nPremise, nRules)

        !% --- allocate and initialize the actionI(:,:) array from control action data
        call control_init_action_points(nThenAction, nElseAction, nRules)

        !% --- set the elements and images for the monitoring points
        call control_init_monitor_elements()

        !% --- set the elements and images for the action points
        call control_init_action_elements()

        ! do ii=1,N_actionPoint
        !     print *, ii, actionI(ii,ai_idx)
        !     print *, 'elem idx   ',actionI(ii,ai_elem_idx)
        !     print *, 'image      ',actionI(ii,ai_image)
        !     print *, 'link idx   ',actionI(ii,ai_link_idx)
        !     print *, 'haschanged ',actionI(ii,ai_hasChanged)
        ! end do
        ! stop 298734


    end subroutine control_init_monitoring_and_action_from_EPASWMM
!%    
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine control_update_monitor_across_images()
        !%------------------------------------------------------------------
        !% Description:
        !% stores the latest element data across all images in the
        !% monitorR array.
        !%------------------------------------------------------------------
        !% Declarations
            integer :: ii, kk
            real(8), target :: rdummy
            integer, pointer :: Eidx
            character(64) :: subroutine_name = 'control_update_monitor_across_images'
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        rdummy = nullvalueR

        ! print *, 'in ',trim(subroutine_name)
        ! do ii=1,size(monitorI,1)
        !     print *, ' '
        !     print *, 'ii = ',ii
        !     print *, monitorI(ii,mi_idx)
        !     print *, monitorI(ii,mi_image)
        !     print *, monitorI(ii,mi_elem_idx)
        !     print *, monitorI(ii,mi_linknodesimType)
        !     print *, monitorI(ii,mi_linknode_idx)
        ! end do

        !% --- cycle through the monitor points to communicate data
        do ii=1,N_MonitorPoint
            !print *, ii

            if (monitorI(ii,mi_image) == this_image()) then
                Eidx => monitorI(ii,mi_elem_idx)
                !print *, 'Eidx ',Eidx
                !% -- store data in coarray to pass through to other images
                !%    This approach is taken because we cannot broadcast only
                !%    part of an array.
                !%    HACK we might be able to neaten this up by copying only
                !%    to the images that have action, rather than all images.
                monitorPassR(1) = elemR(Eidx,er_Depth)
                monitorPassR(2) = elemR(Eidx,er_Head)
                monitorPassR(3) = elemR(Eidx,er_Volume)  !% HACK -- this will not be the link volume!
                monitorPassR(4) = elemR(Eidx,er_FlowrateLateral) 
                monitorPassR(5) = elemR(Eidx,er_Flowrate)
                monitorPassR(6) = elemR(Eidx,er_Setting)
                monitorPassR(7) = elemR(Eidx,er_TimeLastSet)
                !% --- broadcast to all images
                call co_broadcast(monitorPassR, source_image = ii)     
            else
                !% continue
            end if

            ! print *, 'monitorPassR'
            ! print *, monitorPassR

            !% --- pause all the images and wait for the image corresponding
            !%     to this monitoring point to finish its broadcast
            sync all

            !% --- store the broadcast data to local non-coarrays
            monitorR(ii,mr_Depth)       = monitorPassR(1)
            monitorR(ii,mr_Head)        = monitorPassR(2)
            monitorR(ii,mr_Volume)      = monitorPassR(3)
            monitorR(ii,mr_Inflow)      = monitorPassR(4)
            monitorR(ii,mr_Flow)        = monitorPassR(5)
            monitorR(ii,mr_Setting)     = monitorPassR(6)
            monitorR(ii,mr_TimeLastSet) = monitorPassR(7)
        end do

     
        !% --- at this point, the monitorR array on every image has exactly the same data
        !%     so a control action can be conducted on any image.                
        

    end subroutine control_update_monitor_across_images
 !%    
!%==========================================================================    
!%==========================================================================    
!%     
    subroutine control_update_EPASWMM_monitor_data()
        !%------------------------------------------------------------------
        !% Description
        !% updates the EPA-SWMM monitor data (links/nodes) with the SWMM5+ 
        !% Note that this updates EPA-SWMM on every image, which is a
        !% shotgun approach because the EPA-SWMM controls_evaluate will
        !% evaluate ALL the controls, whether or not they are actually
        !% on the image. If we don't shotgun all the images then we would
        !% have to find a way to parse the different controls to different
        !% images, which would likely cause bugs.
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: ii
            real(8), pointer :: tDepth, tHead, tVolume, tInflow, tFlow
            real(8), pointer :: tSetting, tTimeLastSet
            integer, pointer :: LinkNodeNum, linknodesimType
            character(64) :: subroutine_name = 'control_update_EPASWMM_monitor_data'
        !%------------------------------------------------------------------

        !% DEBUGGING OUTPUT ================================================
            ! print *, ' '
            ! print *, 'in ',trim(subroutine_name), ': MONITOR DATA SENT TO EPA SWMM'    
        !% DEBUGGING OUTPUT ================================================

        do ii=1,N_MonitorPoint
            !% --- get the monitoring data     
            tDepth       => monitorR(ii,mr_Depth)
            tHead        => monitorR(ii,mr_Head)
            tVolume      => monitorR(ii,mr_Volume)
            tInflow      => monitorR(ii,mr_Inflow)
            tFlow        => monitorR(ii,mr_Flow)
            tSetting     => monitorR(ii,mr_Setting)
            tTimeLastSet => monitorR(ii,mr_TimeLastSet)

            LinkNodeNum     => monitorI(ii,mi_linknode_idx)
            linknodesimType => monitorI(ii,mi_linknodesimType)

            !% --- send monitor data into EPA-SWMM
            call interface_controls_transfer_monitor_data &
                (tDepth, tHead, tVolume, tInflow, tFlow, tSetting, &
                 tTimeLastSet, LinkNodeNum, linknodesimType)

            !% DEBUGGING OUTPUT ================================================
                ! print *, 'monitor ii  ',ii
                ! print *, 'Depth       ',tDepth
                ! print *, 'Head        ',tHead
                ! print *, 'Volume      ',tVolume
                ! print *, 'Inflow      ',tInflow
                ! print *, 'Flow        ',tFlow
                ! print *, 'Setting     ',tSetting
                ! print *, 'timeLastSet ',tTimeLastSet
                ! print *, 'linkNodeNum ',LinkNodeNum
                ! print *, 'linknodesimType      ',linknodesimType
            !% DEBUGGING OUTPUT ================================================

        end do

    end subroutine control_update_EPASWMM_monitor_data
!%    
!%==========================================================================     
!%==========================================================================    
!%  
    subroutine control_update_actions ()
        !%------------------------------------------------------------------
        !% Description
        !% gets the updated targetsettings and timelastset from EPA-SWMM after
        !%  control rules have been applied
        !%------------------------------------------------------------------
        !% Declarations
            integer :: ii
            integer, pointer :: Eidx, Lidx
            real(8), pointer :: targetsetting(:), timelastset(:)
            character(64) :: subroutine_name = 'control_update_actions'
        !%------------------------------------------------------------------
        !% Aliases:
            targetsetting => elemR(:,er_TargetSetting)
            timelastset   => elemR(:,er_TimeLastSet)
        !%------------------------------------------------------------------
        !% DEBUGGING OUTPUT ================================================
            ! print *, ' '
            ! print *, 'in ',trim(subroutine_name)
        !% DEBUGGING OUTPUT ================================================

        !% --- cycle through the action points
        do ii=1,N_ActionPoint

                ! print *, 'action point ',ii
            
            !% --- only modify elemR data for the image associated with the action point
            if (actionI(ii,ai_image) == this_image()) then
                !% --- set the element and link indexes
                Eidx => actionI(ii,ai_elem_idx)
                Lidx => actionI(ii,ai_link_idx)
                    ! print *, 'Eidx, Lidx, ',Eidx, Lidx
                    ! print *, 'link name = ', trim(link%Names(Lidx)%str)
                !% --- get the new target setting and time last set
                call interface_controls_get_action_results &
                    (targetsetting(Eidx),timelastset(Eidx),Lidx)

                !% DEBUGGING OUTPUT ================================================
                    ! print *, 'Action Point idx ',Eidx
                    ! print *, 'Link idx         ',Lidx, '; link name = ', trim(link%Names(Lidx)%str)
                    ! print *, 'targetSetting    ',targetsetting(Eidx)
                    ! print *, 'timelastset      ',timelastset(Eidx)
                    ! print *, 'elem setting before action taken ',elemR(Eidx,er_Setting)
                !% DEBUGGING OUTPUT ================================================

            end if  

        end do

    end subroutine control_update_actions
!%    
!%==========================================================================     
!%==========================================================================    
!%   
    subroutine control_update_setting ()
        !%------------------------------------------------------------------
        !% Description
        !% Update of elemR(:,er_Setting) based on elemR(:,er_TargetSetting)
        !% and element type for controls
        !%------------------------------------------------------------------
        !% Declarations
            integer :: ii
            integer, pointer :: Eidx, elemType, hasChanged
            real(8), pointer :: thisSetting, targetSetting, timeLastSet
            logical, pointer :: isClosedConduit
            character(64) :: subroutine_name = 'control_update_setting'
        !%------------------------------------------------------------------

        !% DEBUGGING OUTPUT ================================================
            ! print *, ' '
            ! print *, 'in ',trim(subroutine_name)
        !% DEBUGGING OUTPUT ================================================

        !% --- cycle through the action points (few)
        do ii = 1,N_ActionPoint

            !% --- only conduct actions updates on the appropriate images
            if (actionI(ii,ai_image) .ne. this_image()) cycle

            !% --- Aliases
            Eidx            => actionI(ii,ai_elem_idx)
            hasChanged      => actionI(ii,ai_hasChanged)
            thisSetting     => elemR( Eidx,er_Setting)
            targetSetting   => elemR( Eidx,er_TargetSetting)
            timeLastSet     => elemR( Eidx,er_TimeLastSet)
            elemType        => elemI( Eidx,ei_elementType)
            isClosedConduit => elemYN(Eidx,eYN_canSurcharge)

            !% DEBUGGING OUTPUT ================================================
                ! print *, 'ACTION POINTS before action'
                ! print *, 'Eidx           ',Eidx
                ! print *, 'hasChanged      ', hasChanged
                ! print *, 'thisSetting     ',thisSetting
                ! print *, 'targetSetting   ',targetSetting
                ! print *, 'timeLastSet     ',timeLastSet
                ! print *, 'elemType        ',trim(reverseKey(elemType))
                ! print *, 'isClosedConduit ',isClosedConduit
            !% DEBUGGING OUTPUT ================================================

            !% --- error checking
            if ((targetSetting < zeroR) .or. (targetSetting > oneR)) then
                print *, 'CODE ERROR: target setting from controls is less than 0 or greater than 1'
                call util_crashpoint(778734)
            end if

            !% --- set the timeLastSet to reflect only changes to or from closed
            !%     i.e., it does not count incremental changes (probably should be renamed timeLastClosedOrOpened)
            !%     NOTE that in EPA-SWMM, a closed conduit link should only be
            !%     at 1 or 0 (open or closed), but this is written to allow 
            !%     an incremental change, which is generally applied across all element types
            if (targetSetting .ne. thisSetting) then
                if (targetSetting * thisSetting == zeroR) then   
                    timeLastSet = setting%Time%Now
                end if
            end if

            !% --- update thisSetting for different element types
            !%     Note this does not update any other data.
            select case (elemType)
            case (CC)  !--- channels and conduits
                !% --- NOTE: closing a conduit changes the setting to 0, but
                !%     does NOT change the output at this time step.
                if (isClosedConduit) then
                    !% --- closed conduit that can use a setting of 1 (open) or 0 closed
                    if (targetSetting .ne. thisSetting) then
                        !% --- immediate open/close change for a link
                        thisSetting = targetSetting 
                        hasChanged = oneI
                    else
                        hasChanged = zeroI
                    end if
                    !% --- error check
                    if ((thisSetting .ne. zeroR) .or. (thisSetting .ne. oneR)) then
                        print *, 'CODE ERROR: a closed conduit element has an er_Setting of other than 0.0 or 1.0'
                        call util_crashpoint(598723)
                    end if
                else
                    !% --- not a closed conduit (open channel)
                    print *, 'CODE/USER ERROR: control action point on an open channel, which is not allowed'
                    call util_crashpoint(587223)
                end if

            case (pump)
                !% --- reset pump setting based on targetsetting
                if (targetSetting .ne. thisSetting) then
                    call pump_set_setting (Eidx)
                    hasChanged = oneI
                else
                    hasChanged = zeroI
                end if

            case (weir)
                !% --- reset weir setting based on targetsetting
                if (targetSetting .ne. thisSetting) then
                    call weir_set_setting (Eidx)
                    hasChanged = oneI
                else
                    hasChanged = zeroI
                end if

            case (orifice)
                !% --- reset orifice setting based on targetsetting
                if (targetSetting .ne. thisSetting) then
                    call orifice_set_setting (Eidx)
                    hasChanged = oneI
                else
                    hasChanged = zeroI
                end if

            case (outlet)
                !% --- reset orifice based on targetsetting
                if (targetSetting .ne. thisSetting) then
                    call outlet_set_setting (Eidx)
                    hasChanged = oneI
                else
                    hasChanged = zeroI
                end if

            case (JM)
                print *, 'CODE/USER ERROR: control action point on a junction main, which is not allowed'
                call util_crashpoint(558273)
                
            case (JB)
                print *, 'CODE/USER ERROR control action point on a junction branch, which is not allowed'
                call util_crashpoint(682734)

            case default
                print *, 'CODE ERROR: unexpected element type # ',elemType
                if (elemType .ne. nullvalueI) print *, 'with keyword of ',trim(reverseKey(elemType))
                call util_crashpoint(343723)

            end select


            !% DEBUGGING OUTPUT ================================================
                ! print *, 'ACTION POINTS after action'
                ! print *, 'Eidx           ',Eidx
                ! print *, 'hasChanged      ', hasChanged
                ! print *, 'thisSetting     ',thisSetting
                ! print *, 'targetSetting   ',targetSetting
                ! print *, 'timeLastSet     ',timeLastSet
                ! print *, 'elemType        ',trim(reverseKey(elemType))
                ! print *, 'isClosedConduit ',isClosedConduit
            !% DEBUGGING OUTPUT ================================================           
       
        end do

    end subroutine control_update_setting
!%    
!%==========================================================================     
!%==========================================================================    
!%   
    subroutine control_update_element_values ()
        !%------------------------------------------------------------------
        !% Description:
        !% Update of values in elemR array for changes in elemR(:,er_Setting)
        !% due to controls.
        !% This is handled separately from control_update_setting so that we
        !% can update settings without changing the actual flow conditions
        !% NOTE THIS IS LIMITED TO CONTROL ACTION POINTS ONLY
        !%------------------------------------------------------------------
        !% Declarastions:
            integer :: ii
            integer, pointer :: Eidx, elemType, hasChanged, dface
            real(8), pointer :: thisSetting
            logical, pointer :: isClosedConduit
            character(64) :: subroutine_name = 'control_update_element_values'
        !%------------------------------------------------------------------

            !% DEBUGGING OUTPUT ================================================
                ! print *, ' '
                ! print *, 'in ',trim(subroutine_name)
            !% DEBUGGING OUTPUT ================================================

        !% --- cycle through the action points    
        do ii = 1,N_ActionPoint
            !print *, ii ,' in ',trim(subroutine_name)
            !print *, 'actionI(ii,ai_image) ',actionI(ii,ai_image)

            !% --- only conduct actions updates on the appropriate images
            if (actionI(ii,ai_image) .ne. this_image()) cycle

            !% --- Aliases
            Eidx            => actionI(ii,ai_elem_idx)
            hasChanged      => actionI(ii,ai_hasChanged)
            elemType        => elemI(Eidx,ei_elementType)
            thisSetting     => elemR(Eidx,er_Setting)

            ! print *, 'Eidx ',Eidx
            ! print *, 'hasChanged ',hasChanged
            ! print *, 'elemType   ',trim(reverseKey(elemType))
            ! print *, 'thisSetting',thisSetting
            
            select case (elemType)

            case (CC)
                !% --- for regular CC elements, whether or not the setting has changed is irrelevant
                !%     we force zero flux on downstream faces that have Setting = 0
                !%     note that this is also done as part of face interpolation.
                if (thisSetting == zeroR) then
                    dface => elemI(Eidx,ei_Mface_dL)
                    call adjust_face_for_zero_setting_singular(dface)
                end if

            case (pump)
                if (hasChanged == oneI) then 
                    call pump_toplevel(Eidx)
                else
                    !% no change 
                end if

            case (weir)
                if (hasChanged == oneI) then 
                    call weir_toplevel(Eidx)
                else
                    !% no change 
                end if

            case (orifice)
                !print *, 'in here 98273'
                if (hasChanged == oneI) then 
                    call orifice_toplevel(Eidx)
                else
                    !% no change 
                end if
                !print *, 'done here 559873'

            case (outlet)
                if (hasChanged == oneI) then 
                    call outlet_toplevel(Eidx)
                else
                    !% no change 
                end if

            case (JM)
                print *, 'CODE/USER ERROR: control action point on a junction main, which is not allowed'
                call util_crashpoint(558273)
                
            case (JB)
                print *, 'CODE/USER ERROR control action point on a junction branch, which is not allowed'
                call util_crashpoint(682734)

            case default
                print *, 'CODE ERROR: unexpected element type # ',elemType
                if (elemType .ne. nullvalueI) print *, 'with keyword of ',trim(reverseKey(elemType))
                call util_crashpoint(343723)

            end select

            !% DEBUGGING OUTPUT ================================================
                ! print *, ' Eidx, haschanged ',Eidx, hasChanged, trim(reverseKey(elemType))
                ! print *, ' thissetting ',thissetting
                ! print *, 'elem setting ',elemR(Eidx,er_setting)
            !% DEBUGGING OUTPUT ================================================

        end do

        !print *, 'leaving ',trim(subroutine_name)
        
    end subroutine control_update_element_values
!%    
!%==========================================================================     
!%==========================================================================    
!%         
    subroutine control_init_monitor_points (nPremise, nRules)
        !%------------------------------------------------------------------
        !% Description:
        !% uses premises in the control logic of EPA-SWMM to identify
        !% the unique locations in the link and node arrays that are used
        !% as monitoring points for controls
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: nRules, nPremise
            integer :: ii, rr
            integer :: iLeft, iRight, thisPremiseLevel, success, numUnique
            integer :: npoint
            integer, allocatable :: location(:), linknodesimType(:), attribute(:)
            integer, allocatable :: irank(:), simulation(:)
            character(512) :: thisName
            character(64) :: subroutine_name = 'control_init_monitor_points'
        !%------------------------------------------------------------------
        !% --- each premise may add two monitoring points (LHS and RHS)
        npoint = twoI * nPremise

        !% --- if no premise points, exit
        if (npoint <= 0) return

        !% --- allocate temporary arrays
        allocate(location(npoint))
        allocate(linknodesimType(npoint))
        allocate(attribute(npoint))
        allocate(irank(npoint))
        location(:)        = nullvalueI
        irank(:)           = nullvalueI
        linknodesimType(:) = nullvalueI
        attribute(:)       = nullvalueI

        !% --- initial indexes for LHS and RHS of premise data
        iLeft   = 1
        iRight  = 2
        !% --- cycle through the rules (they index from 0 in C)
        do rr = 0, nRules-1
            !% --- each rule starts at the 0 premise level
            thisPremiseLevel = 0
                ! print *, ' '
                ! print *, 'rr ',rr, thisPremiseLevel
            success = 1
            do while (success == 1)
                    ! print *, ' '
                    ! print *, ' calling interface_controls_get_premise_data'
                !% --- get the premise level data, and (if successful) increment thisPremiseLevel
                call interface_controls_get_premise_data (                    &
                    location       (iLeft), location       (iRight),          &
                    linknodesimType(iLeft), linknodesimType(iRight),          &
                    attribute      (iLeft), attribute      (iRight),          &
                    thisPremiseLevel, rr, success)
                !% --- increment the storage location for the next premise data    
                !%     if not successful the last location was not used
                    !   print *, 'success ',success
                if (success) then
                        ! print *, 'thisPremiseLevel', thisPremiseLevel
                        ! print *, iLeft, iRight
                        ! print *, location(iLeft), location (iRight)
                        ! print *, linknodesimType(iLeft), linknodesimType(iRight)
                        ! print *, attribute(iLeft),attribute(iRight)
                    iLeft  = iLeft  + twoI
                    iRight = iRight + twoI

                    !% --- error checking
                    if (iRight > size(location)) then
                        print *, 'CODE ERROR: more premises found that space allocated for storage '
                        call util_crashpoint(3298734)
                    end if
                else
                    !% --- reset values
                    location(iLeft)  = nullvalueI
                    location(iRight) = nullvalueI
                    linknodesimType(iLeft)    = nullvalueI
                    linknodesimType(iRight)   = nullvalueI
                    attribute(iLeft) = nullvalueI
                    attribute(iRight)= nullvalueI
                end if
                    ! print *, 'end of while'
            end do
        end do

            ! print *, ' '
            ! print *, 'CONTROL monitoring points'
            ! do ii=1,size(location)
            !     print *, location(ii), linknodesimType(ii), attribute(ii)
            ! end do

        !% --- find unique locations for premises, 
        !      note the max number is nullvalueI
        call util_unique_rank(location,irank,numUnique)

            ! !% --- testing: printout the monitoring points
            ! print *, 'CONTROL Unique monitoring points'
            ! do ii=1,numUnique
            !     print *, location(irank(ii)), linknodesimType(irank(ii)), attribute(irank(ii))
            ! end do

            ! print *, ' LINKS'
            ! do ii = 1,N_link
            !     print *, ii, link%I(ii,li_idx), trim(reverseKey(link%I(ii,li_link_type))), ' , ', trim(link%Names(ii)%str)
            ! end do

            ! print *, ' NODES'
            ! do ii = 1,N_node
            !     print *, ii, node%I(ii,ni_idx), trim(reverseKey(node%I(ii,ni_node_type))), ' , ', trim(node%Names(ii)%str)
            ! end do

        !% --- set the number of monitoring points (removing the nullvalueI)
        !%     check for nullvalue in the unique location
        if (location(irank(numUnique)) == nullValueI) then
            N_MonitorPoint = numUnique-1
        else
            N_MonitorPoint = numUnique
        end if
        
        !% DEBUGGING OUTPUT ============================================
            ! print *, ' '
            ! print *, 'in control_init_monitor_points: MONITOR POINTS'
            ! do ii=1,N_monitorPoint
            !     if (linknodesimType(irank(ii))) then
            !         thisName = trim(link%Names(location(irank(ii)))%str)
            !     else
            !         thisname = trim(node%Names(location(irank(ii)))%str)
            !     end if
            !     write(*,"(A,i6,A,4i6)"),' LinkNode #, linkNode Name, linknodesimType, attr: ',location(irank(ii)),trim(thisName),  linknodesimType(irank(ii)), attribute(irank(ii))
            ! end do
        !% END DEBUGGING ================================================

        if (N_MonitorPoint < 1) then
            !% --- this occurs if only SIMULATION premises are found
            !%     i.e., all the linknodesimTyp = -1
            !% --- temporary set of # monitor points to 1 to allocate array
            N_MonitorPoint = 1
            call util_allocate_monitor_points()
            N_MonitorPoint = 0
        else
            !% --- allocate the monitor point storage
            call util_allocate_monitor_points()

            !% --- store the unique monitoring points
            !%     note, attributes are NOT stored because non-unique points have different
            !%     attributes. If we later need attributes we will need non-unique storage
            monitorI(:,mi_linknode_idx)    = location(irank(1:N_MonitorPoint))
            monitorI(:,mi_linknodesimType) = linknodesimType(irank(1:N_MonitorPoint))

            !% --- error checking
            !      print *, 'checking that link/node location is in index set'
            do ii=1,N_MonitorPoint
                !print *, 'monitorI linknodesimType: ',int(monitorI(ii,mi_linknodesimType))
                select case (monitorI(ii,mi_linknodesimType))
                case (0)  !node
                    rr = count(node%I(:,ni_idx) == monitorI(ii,mi_linknode_idx))
                case (1)  ! link
                    rr = count(link%I(:,li_idx) == monitorI(ii,mi_linknode_idx))
                case default
                    print *, 'CODE ERROR: unexpected default case'
                    call util_crashpoint(448229)
                end select
                if (rr .ne. 1) then
                    print *, 'CODE OR DATA ERROR: monitor point does not match node or link indexes'
                    call util_crashpoint(598272)
                end if
            end do
        end if

        !%------------------------------------------------------------------
        !% Closing:
            deallocate(irank)
            deallocate(location)
            deallocate(linknodesimType)
            deallocate(attribute)

    end subroutine control_init_monitor_points
!%    
!%==========================================================================
!%==========================================================================
!%  
    subroutine control_init_action_points (nThenAction, nElseAction, nRules)
        !%------------------------------------------------------------------
        !% Description:
        !% initializes the actionI array from then/else actions in control
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: nThenAction, nElseAction, nRules
            integer :: npoint, rr, ii, kIdx, thisActionLevel, success, numUnique
            integer :: isThen
            integer, allocatable :: location(:), linknodesimType(:), attribute(:)
            integer, allocatable :: irank(:)
            character(64) :: subroutine_name = 'control_init_action_points'
        !%------------------------------------------------------------------
        if (nRules < 1) return

        !% --- total number of action points that are possible
        npoint = nThenAction + nElseAction

        if (npoint < 1) return

        !% --- allocate temporary arrays
        allocate(location(npoint))
        allocate(linknodesimType(npoint))
        allocate(attribute(npoint))
        allocate(irank(npoint))
        location(:)  = nullvalueI
        irank(:)     = nullvalueI
        linknodesimType(:)    = nullvalueI
        attribute(:) = nullvalueI

        kIdx = 1 !% index for the action storage

        !% --- cycle through the rules (they index from 0 in C)
        do rr = 0, nRules-1
            !% --- get the "then" actions
            !% --- each rule starts at the 0 action level
            thisActionLevel = 0
                ! print *, ' '
                ! print *, 'rr ',rr, thisActionLevel
            success = 1
            do while (success == 1)
                isThen = 1
                call interface_controls_get_action_data (   &
                    location (kIdx),                             &
                    attribute(kIdx),                             &
                    thisActionLevel, rr, success,isThen)  
                if (success) kIdx = kIdx +1  
            end do    

            !% --- get the "else" actions
            thisActionLevel = 0
            success = 1
            do while (success == 1)
                isThen = 0
                call interface_controls_get_action_data (   &
                    location (kIdx),                             &
                    attribute(kIdx),                             &
                    thisActionLevel, rr, success,isThen)  
                if (success) kIdx = kIdx +1  
            end do   
        end do

            ! print *, 'all action points in ',trim(subroutine_name)
            ! do ii=1,npoint
            !     print *, ii, location(ii), attribute(ii)
            ! end do

        !% --- find unique locations for actions, 
        !      note the max number is nullvalueI
        call util_unique_rank(location,irank,numUnique)

            ! !% --- testing: printout the action points
            ! print *, 'Unique action points'
            ! do ii=1,numUnique
            !     print *, ii, location(irank(ii))
            ! end do
            ! print *, 'end unique'

        !% --- set the number of action points (removing the nullvalueI)
        !%     check for nullvalue in the unique location
        if (location(irank(numUnique)) == nullValueI) then
            N_ActionPoint = numUnique-1
        else
            N_ActionPoint = numUnique
        end if

            ! print *, 'N_actionPoint ',N_ActionPoint

            ! !% --- testing: printout the action points
            ! print *, 'Unique action points '
            ! do ii=1,N_ActionPoint
            !     print *, location(irank(ii)),  attribute(irank(ii))
            ! end do

        !% --- error checking
        if (N_ActionPoint < 1) then
            print *, 'CODE OR DATA ERROR: expected at least 1 action point from controls'
            call util_crashpoint(77632)
        end if
       
        call util_allocate_action_points()

        !% --- store the unique action points
        actionI(:,ai_link_idx) = location(irank(1:N_ActionPoint))
        actionI(1:N_ActionPoint,ai_idx) = (/ 1:N_ActionPoint /)

        !% DEBUGGING OUTPUT ================================================
            ! print *, ' '
            ! print *, 'in control_init_action_points: ACTION POINTS'
            ! print *, 'all action points'
            ! do ii = 1,N_ActionPoint
            !     write(*,"(A,3i6)") ' actionID, idx ',ii, actionI(ii,ai_idx)
            !     print *, 'Link #, link Name ', actionI(ii,ai_link_idx), '; ', trim(link%Names(actionI(ii,ai_link_idx))%str)
            ! end do
        !% END DEBUGGING ===================================================o

        !% --- error checking
            ! print *, 'checking that link location is in index set'
        do ii=1,N_ActionPoint
            rr = count(link%I(:,li_idx) == actionI(ii,ai_link_idx))
            if (rr .ne. 1) then
                print *, 'CODE OR DATA ERROR: action point does not match link indexes'
                call util_crashpoint(398743)
            end if
        end do

        !%------------------------------------------------------------------
        !% Closing
            deallocate(irank)
            deallocate(location)
            deallocate(linknodesimType)
            deallocate(attribute)

    end subroutine control_init_action_points
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine control_init_monitor_elements ()
        !%------------------------------------------------------------------
        !% Description:
        !% sets the element and image locations of monitor 
        !% elements of EPA-SWMM controls on link/nodes
        !% If a link is composed of more than one element, this algorithm will
        !% choose the central element as the monitor location.
        !% This algorithm uses the EPA-SWMM link indexes, which are the
        !% parent links in the partitioning algorithm. Because we do not
        !% enquire of parent/child status, the monitor point will always be
        !% on the parent link, which retains the original EPA-SWMM link
        !% index number. If more refined control is needed for a monitoring
        !% point, the user must provide further subdivision of the links in 
        !% the EPA-SWMM *.inp file or choose a node for the control.
        !%------------------------------------------------------------------
        !% Declarations
            integer :: ii
            integer, pointer :: linknodesimType(:), LNidx(:), eIdx(:)
            integer, pointer :: nodeType(:), numElement(:)
            integer, pointer :: monitorImage(:), linkImage(:), nodeImage(:)
            integer, pointer :: elemStart(:), elemEnd(:), Lidx, Nidx
            character(64) :: subroutine_name = 'control_init_monitor_elements'
        !%------------------------------------------------------------------
        !% Aliases
            linknodesimType => monitorI(:,mi_linknodesimType)
            LNidx  => monitorI(:,mi_linknode_idx)
            eIdx   => monitorI(:,mi_elem_idx)
            monitorImage => monitorI(:,mi_image)
            numElement   => link%I(:,li_N_element)
            elemStart    => link%I(:,li_first_elem_idx)
            elemEnd      => link%I(:,li_last_elem_idx)
            nodeType     => node%I(:,ni_node_type)
            linkImage    => link%I(:,li_P_image)
            nodeImage    => node%I(:,ni_P_image)
        !%------------------------------------------------------------------

        !% --- cycle through the monitoring points    
        do ii=1,N_MonitorPoint
            monitorI(ii,mi_idx) = ii

            !% --- each point is either associated with a link or node
            select case (linknodesimType(ii))
            case (1) !% is link
                Lidx => LNidx(ii) !% --- EPA-SWMM link index
                monitorImage(ii) = linkImage(Lidx)
                if (numElement(Lidx) == 1) then
                    !% --- only one element in link, so use that as monitor element
                    eIdx(ii)  = elemStart(Lidx)
                else
                    !% --- choose the central element 
                    !%     note that integer division gives bias to the
                    !%     smaller of two central values (upstream) if an even number
                    !%     of elements
                    eIdx(ii) = (elemStart(Lidx) + elemEnd(Lidx)) / twoI
                end if
            case (0) !% is node
                Nidx => LNidx(ii)
                !print *, 'nidx ',Nidx
                !print *, 'node type ',trim(reverseKey(nodeType(Nidx)))
                monitorImage(ii) = nodeImage(Nidx)
                select case (nodeType(Nidx))

                case (nJ1,nBCup)
                    !% --- connect monitor for node to the first element of downstream link
                    Lidx => node%I(Nidx,ni_N_link_d)
                    if (Lidx == nullvalueI) then
                        print *, 'CODE/SYSTEM ERROR: unexpected nullvalue for link index'
                        call util_crashpoint(598723)
                    else
                        eIdx(ii) = elemStart(Lidx)
                    end if

                case (nJ2, nBCdn)
                    !% --- connect monitor for node to the last element of upstream link
                    Lidx => node%I(Nidx,ni_N_link_u)
                    if (Lidx == nullvalueI) then
                        print *, 'CODE/SYSTEM ERROR: unexpected nullvalue for link index'
                        call util_crashpoint(98273)
                    else
                        eIdx(ii) = elemEnd(Lidx)
                    end if

                case (nJM,nStorage)
                    !if (node%I(Nidx,ni_elemface_idx) == nullvalueI) then
                    if (node%I(Nidx,ni_elem_idx) == nullvalueI) then
                        print *, 'CODE/SYSTEM ERROR: unexpected nullvalue for node index'
                        call util_crashpoint(429933)
                    else
                        !eIdx(ii) =  node%I(Nidx,ni_elemface_idx)
                        eIdx(ii) =  node%I(Nidx,ni_elem_idx)
                    end if

                case default
                    print *, 'CODE ERROR: Unexpected case default, ni_node_type of ',trim(reverseKey(nodeType(Nidx)))
                end select

            case default
                print *, 'CODE ERROR: Unexpected case default, monitor(:,mi_linknodesimType) unsupported value of ',linknodesimType(ii)
                call util_crashpoint(58723)
            end select

            !% DEBUGGING OUTPUT ================================================
                ! print *, ' '
                ! print *, 'in ',trim(subroutine_name), ': MONITOR ELEMENTS'
                ! print *, 'monitor ii=  ',ii
                ! print *, 'linknodesimType       ',linknodesimType(ii)
                ! print *, 'element #    ',eIdx(ii)
                ! print *, 'element Type ',trim(reverseKey(elemI(eIdx(ii),ei_elementType)))
                ! if (linknodesimType(ii)) then
                !     print *, 'link #, type ',LNidx(ii), trim(link%Names(LnIdx(ii))%str)
                ! else
                !     print *, 'node #, type ',LNidx(ii), trim(node%Names(LnIdx(ii))%str)
                ! end if
            !% END DEBUGGING ===================================================

        end do

    end subroutine control_init_monitor_elements
!%    
!%==========================================================================
!%==========================================================================
!%    
    subroutine control_init_action_elements ()
        !%------------------------------------------------------------------
        !% Description
        !% assigns the elements and images for the action links of control
        !% actions.
        !% If a link is composed of more than one element, this algorithm will
        !% choose the central element as the monitor location.
        !% This algorithm uses the EPA-SWMM link indexes, which are the
        !% parent links in the partitioning algorithm. Because we do not
        !% enquire of parent/child status, the monitor point will always be
        !% on the parent link, which retains the original EPA-SWMM link
        !% index number. If more refined control is needed for an action
        !% point, the user must provide further subdivision of the links in 
        !% the EPA-SWMM *.inp file
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: ii
            integer, pointer :: eIdx(:), Lidx, actionImage(:), linkImage(:)
            integer, pointer :: numElement(:), elemStart(:), elemEnd(:)
            character(64) :: subroutine_name = 'control_init_action_elements'
        !%------------------------------------------------------------------
        !% Aliases
            eIdx         => actionI(:,ai_elem_idx)
            actionImage  => actionI(:,ai_image)
            numElement   => link%I(:,li_N_element)
            elemStart    => link%I(:,li_first_elem_idx)
            elemEnd      => link%I(:,li_last_elem_idx)
            linkImage    => link%I(:,li_P_image)
        !%------------------------------------------------------------------

        !% DEBUGGING OUTPUT ================================================
            ! print *, 'in ',trim(subroutine_name), ': ACTION ELEMENTS'
        !% DEBUGGING OUTPUT ================================================

        !% --- cycle through the action points
        do ii=1,N_ActionPoint
            !print *, ii, actionI(ii,ai_link_idx)

            !% --- get the link for this action point
            Lidx => actionI(ii,ai_link_idx)
            !print *, 'link index ',Lidx

            !% --- store the image that this link is partitioned to
            actionImage(ii) = linkImage(Lidx)
            !print *, 'actionImage ',actionImage(ii)
            !print *, 'numElement  ',numElement(Lidx)

            !% --- find the element on this link to use
            if (numElement(Lidx) == 1) then
                !% --- only one element in link, so use that as action element
                eIdx(ii) = elemStart(Lidx)
                !print *, 'elemstart ',elemStart(Lidx)
            else
                !% --- choose central element
                eIdx(ii) = (elemStart(Lidx) + elemEnd(Lidx)) / twoI
            end if

            !% DEBUGGING OUTPUT ================================================
                ! print *, 'action ii      ',ii
                ! print *, 'action Link     ',Lidx
                ! print *, 'element #, type ',eIdx(ii),trim(reverseKey(elemI(eIdx(ii),ei_elementType)))
            !% END DEBUGGING ===================================================

        end do

    end subroutine control_init_action_elements
!%    
!%==========================================================================
!% END OF MODULE
!%==========================================================================
!%
end module control_hydraulics