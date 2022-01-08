module output

    use define_indexes
    use define_keys
    use define_globals
    use define_settings
    use define_types
    use interface
    use utility_datetime
    use utility_allocate
    use utility_deallocate


    !use, intrinsic :: iso_fortran_env, only: *

    !%-----------------------------------------------------------------------------
    !% Description
    !% Provides main output for two different types of output schemes
    !% (1) outputML_ is a multi-level output scheme of limited data
    !% (2) outputD_ is a csv -only scheme that dumps everything at every time step - obsolete
    !% (3) output_COMMON are common to both schemes -- obsolete
    implicit none

    private

    !% public subroutines common to both
    !rm brh20211207 public :: output_COMMON_nodes_selection
    !rm brh20211207 public :: output_COMMON_links_selection

    public :: outputML_selection
    public :: outputML_setup

    !public :: outputML_element_selection
    !public :: outputML_face_selection
    !public :: outputML_size_OutElem_by_image
    !public :: outputML_size_OutFace_by_image
    !public :: outputML_element_outtype_selection
    !public :: outputML_face_outtype_selection
    public :: outputML_store_data
    public :: outputML_write_control_file
    public :: outputML_combine_and_write_data
    public :: outputML_convert_elements_to_linknode_and_write
    

    ! public :: outputD_read_csv_link_names
    ! public :: outputD_read_csv_node_names
    ! public :: outputD_create_link_files
    ! public :: outputD_create_node_files
    ! public :: outputD_write_link_files
    ! public :: outputD_write_node_files
    ! public :: outputD_combine_links
    ! public :: outputD_move_node_files
    ! public :: outputD_update_swmm_out

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine outputML_selection ()
        !%-------------------------------------------------------------------
        !% Description:
        !% provides selection of elements/faces of the multi-level output
        !%-------------------------------------------------------------------
        !%-------------------------------------------------------------------

        !% --- designate the corresponding elements for output
        call outputML_element_selection ()
        !% --- deisgnate the corresponding face to output
        call outputML_face_selection ()
      
    end subroutine outputML_selection
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine outputML_setup ()
        !%-------------------------------------------------------------------
        !% Description:
        !% provides setup of the multi-level output
        !%-------------------------------------------------------------------
            integer :: nMaxElem, nMaxFace
            integer :: allocation_status
            character(len=99) :: emsg
        !%-------------------------------------------------------------------

        !% --- compute the N_OutElem for each image
        !print *, 'outputML_size_OutElem ', this_image()
        call outputML_size_OutElem_by_image ()

        !% --- compute the N_OutFace for each image
        !print *, 'outputML_size_OutFace ', this_image()
        call outputML_size_OutFace_by_image ()

        !% --- setup the output element data types
        !print *, 'outputML_element_outtype_selection ', this_image()
        call outputML_element_outtype_selection ()

        !% -- setup the output face data types
        !print *, 'outputML_face_outtype_selection ',this_image()
        call outputML_face_outtype_selection ()

        !% --- create storage space for multi-level output data
        !print *, 'util allocate_outputML_storage ', this_image()
        call util_allocate_outputML_storage ()

        !% --- create storage for output times
        !print *, 'util_allocate_outputML_times ', this_image()
        call util_allocate_outputML_times ()

        !% --- create storage for output binary filenames
        !print *, 'util allocate outputML_filenames ', this_image()
        call util_allocate_outputML_filenames ()

        ! !% allocations that were previously made in output.f90  brh20211223
        ! nMaxElem = maxval(N_OutElem)
        ! allocate(thisElementOut(nMaxElem), stat=allocation_status, errmsg=emsg)
        ! call util_allocate_check(allocation_status, emsg, 'thisElementOut')

        ! nMaxFace = maxval(N_OutFace)
        ! allocate(thisFaceOut(nMaxFace), stat=allocation_status, errmsg=emsg)
        ! call util_allocate_check(allocation_status, emsg, 'thisFaceOut')

    end subroutine outputML_setup
!%
!%==========================================================================    
!%==========================================================================
!%
!    subroutine output_COMMON_nodes_selection ()
        !% brh 20211207 -- obsolete
        !% Now using the api_node_rptFlag to initialize
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Gets the nodes that are output from the SWMM.inp file
        !% Stores T/F in the YN in the node%YN(:,nYN_isOutput) column
        !%-----------------------------------------------------------------------------
         !   integer :: ii
        !%-----------------------------------------------------------------------------
        !% HACK -- presently only handles ALL
        !node%YN(1:N_node,nYN_isOutput) = .true.

        !print *, 'printing the node output'
        !do ii = 1,N_node
        !    print *, ii, node%YN(ii,nYN_isOutput)
        !end do

        !print *, 'TESTING TO SEE WHAT HAPPENS IF ONLY ONE NODE IS OUTPUT'
        !node%YN(1:N_node,nYN_isOutput) = .false.

        !node%YN(1,nYN_isOutput)  = .true.
        !node%YN(2,nYN_isOutput)  = .true.
        !node%YN(3,nYN_isOutput)  = .true.
        !node%YN(4,nYN_isOutput)  = .true.

!    end subroutine output_COMMON_nodes_selection
!%
!%==========================================================================
!%==========================================================================
!%
!   subroutine output_COMMON_links_selection ()
        !% brh 20211207 -- obsolete
        !% Now using the api_link_rptFlag to initialize
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Gets the links that are output from the SWMM.inp file
        !% Stores T/F in the YN in the link%YN(:,lYN_isOutput) column
        !%-----------------------------------------------------------------------------
        ! integer :: ii
        !%-----------------------------------------------------------------------------
        !% HACK -- presently only handles ALL
        !link%YN(1:N_link,lYN_isOutput) = .true.

        !print *, 'printing the link output'
        !do ii = 1,N_link
        !    print *, ii, link%YN(ii,lYN_isOutput)
        !end do

        !print *, 'TESTING TO SEE WHAT HAPPENS IF ONLY ONE LINK IS OUTPUT'

        !link%YN(1:N_link,lYN_isOutput) = .false.

        !link%YN(1,lYN_isOutput)  = .true.
        ! link%YN(2,lYN_isOutput)  = .true.
        ! link%YN(3,lYN_isOutput)  = .true.

!    end subroutine output_COMMON_links_selection
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine outputML_element_selection ()
        !%------------------------------------------------------------------
        !% Description:
        !% Gets the elements corresponding to the nodes and links of SWMM output
        !% Stores T/F in the YN in the elemYN(:,eYN_isOutput) column
        !%------------------------------------------------------------------
            integer :: ii
            integer, pointer :: elementType(:), link_idx(:), node_idx(:), tlink, tnode
            logical, pointer :: isLinkOut(:), isNodeOut(:), isElemOut(:)
            character(64) :: subroutine_name = 'outputML_element_selection'
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (setting%Debug%File%output) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases:
            elementType => elemI(:,ei_elementType)
            link_idx    => elemI(:,ei_link_Gidx_SWMM)
            node_idx    => elemI(:,ei_node_Gidx_SWMM)
            isLinkOut   => link%YN(:,lYN_isOutput)
            isNodeOut   => node%YN(:,nYN_isOutput)
            isElemOut   => elemYN(:,eYN_isOutput)
        !%------------------------------------------------------------------
        !% ensure we start with all output elements to false
        isElemOut(:) = .false.

        !% --- Translate the link%YN to the elemYN
        !% --- HACK brute force with do loop
        !% --- note that each image has a different number of elements
        do ii =1,N_elem(this_image())
            tlink => link_idx(ii)
            select case (elementType(ii))
                case (CC)
                    !% --- conduits always correspond to links
                    if (isLinkOut(tlink)) then
                        isElemOut(ii) = .true.
                    else
                        isElemOut(ii) = .false.
                    end if
                case (orifice)
                    !% --- orifices always correspond to links
                    if (isLinkOut(tlink)) then
                        isElemOut(ii) = .true.
                    else
                        isElemOut(ii) = .false.
                    end if
                case (weir)
                    !% --- weirs always correspond to links
                    if (isLinkOut(tlink)) then
                        isElemOut(ii) = .true.
                    else
                        isElemOut(ii) = .false.
                    end if
                case (outlet)
                    !% --- outlets that correspond to links
                    if (isLinkOut(tlink)) then
                        isElemOut(ii) = .true.
                    else
                        isElemOut(ii) = .false.
                    end if
                case (JM)
                    tnode => node_idx(ii)
                    !% --- junction mains correspond to nodes
                    if (isNodeOut(tnode)) then
                        isElemOut(ii) = .true.
                    else
                        isElemOut(ii) = .false.
                    end if
                case (JB)
                    !% --- no output on junction branches
                    isElemOut(ii) = .false.
                case default
                    write (*,"(A)") 'ERROR (code): statement should be unreachable'
                    write (*,"(A)") ' invalid element type was ',elementType(ii)
                    stop
            end select
        end do
        !%------------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine outputML_element_selection
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine outputML_face_selection ()
        !%------------------------------------------------------------------
        !% Description:
        !% Gets the faces corresponding to the nJ2 nodes of SWMM output
        !% Stores T/F in the YN in the elemYN(:,eYN_isOutput) column
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: ii
            integer, pointer :: node_idx(:), face_idx(:)
            logical, pointer :: isNodeOut(:), isFaceOut(:)
            character(64) :: subroutine_name = 'outputML_face_selection'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%output) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases
            face_idx    => faceI(:,fi_Gidx)
            node_idx    => faceI(:,fi_node_idx_SWMM)
            isNodeOut   => node%YN(:,nYN_isOutput)
            isFaceOut   => faceYN(:,fYN_isFaceOut)
        !%------------------------------------------------------------------
        !% --- Translate the node%YN to the faceYN
        !% --- HACK brute force with do loop
        !% --- note that each image has a different number of faces
        do ii =1,N_face(this_image())
            !% --- check for valid nodes
            !% --- this should only be nJ1, nJ2, nBCup and nBCdn nodes
            !% --- (note that nJm nodes are already handled with elements)
            if (node_idx(ii) .ne. nullvalueI) then
                !% --- check that node is an output node
                if (isNodeOut(node_idx(ii))) then
                    !% --- check that this node has a SWMM name
                    if (.not. (node%Names(node_idx(ii))%str == "")) then
                        !% --- assign face output
                        isFaceOut(ii) = .true.
                    end if
                end if
            end if
        end do
        !%------------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine outputML_face_selection
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine outputML_size_OutElem_by_image ()
        !%-----------------------------------------------------------------------------
        !% Description
        !% computes and allocates N_OutElem for each image (number of output elements)
        !%-----------------------------------------------------------------------------
        !% Declarations
            integer :: ii
            character(64)    :: subroutine_name = 'outputML_size_OutElem_by_image'
        !%-----------------------------------------------------------------------------
        !% Preliminaries:
            if (setting%Output%Report%suppress_MultiLevel_Output) return  
            if (setting%Debug%File%output) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"       
        !%----------------------------------------------------------------------------
        !% initialized the logical controls
        setting%Output%ElementsExist_global = .false.

        !% --- the packed number of elements is the number with "true" in eYN_isOutput
        !%     for each image. Note that the array is initialized to 0, so the value
        !%     will only be non-zero for the operating image.
        N_OutElem(this_image()) = npack_elemP(ep_Output_Elements)
        sync all
        !%--- distribute the max number of output elements over all images
        !%    this ensures each image knows how many elements are output in the others
        call co_max(N_OutElem)
        sync all

        !% use the N_OutElem to set the logical array
        do ii=1,size(N_OutElem)
            if (N_OutElem(ii) > 0) then
                setting%Output%ElementsExist_byImage(ii) = .true.
            else
                setting%Output%ElementsExist_byImage(ii) = .false.
            end if
        end do

        !% set the global logical
        if (any(setting%Output%ElementsExist_byImage)) setting%Output%ElementsExist_global = .true.

        ! print *, setting%Output%ElementsExist_byImage
        ! print *, setting%Output%ElementsExist_global
        ! print *, N_OutElem
    

        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine outputML_size_OutElem_by_image
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine outputML_size_OutFace_by_image ()
        !%------------------------------------------------------------------
        !% Description
        !% computes and allocates N_OutFace for each image (number of output faces)
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: ii
            character(64)    :: subroutine_name = 'outputML_size_OutFace_by_image'
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (setting%Output%Report%suppress_MultiLevel_Output) return
            if (setting%Debug%File%output) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% initialized the logical controls
        setting%Output%FacesExist_global = .false.

        !% the packed number of faces is the number with "true" in eYN_isOutput
        !% for each image
        N_OutFace(this_image()) = npack_faceP(fp_Output_Faces)
        sync all
        !% distribute the number of output faces over all images
        call co_max(N_OutFace)
        sync all

        !% use the N_outFace to set the logical array
        do ii=1,size(N_OutFace)
            if (N_OutFace(ii) > 0) then
                setting%Output%FacesExist_byImage(ii) = .true.
            else 
                setting%Output%FacesExist_byImage(ii) = .false.
            end if
        end do

        !% set the global logical
        if (any(setting%Output%FacesExist_byImage)) setting%Output%FacesExist_global = .true.

        ! print *, setting%Output%FacesExist_byImage
        ! print *, setting%Output%FacesExist_global
        ! print *, N_OutFace
        ! stop 3648053

        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine outputML_size_OutFace_by_image
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine outputML_element_outtype_selection ()
        !%------------------------------------------------------------------
        !% Description
        !% initializes the types of output data provided for elements
        !% Note that this must be done on every image, including those that
        !% don't have any output elements.
        !%------------------------------------------------------------------
            integer :: ii
            character(64)        :: subroutine_name = 'outputML_element_outtype_selection'
        !%------------------------------------------------------------------
            if (setting%Output%Report%suppress_MultiLevel_Output) return
            if (.not. setting%Output%ElementsExist_global) return 
            if (setting%Debug%File%output) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        N_OutTypeElem = 0
        !% --- count the number of true output element types
        if (setting%Output%DataOut%isAreaOut)           N_OutTypeElem =  N_OutTypeElem + 1
        if (setting%Output%DataOut%isDepthOut)          N_OutTypeElem =  N_OutTypeElem + 1
        if (setting%Output%DataOut%isFlowrateOut)       N_OutTypeElem =  N_OutTypeElem + 1
        if (setting%Output%DataOut%isFluxConsOut)       N_OutTypeElem =  N_OutTypeElem + 1
        if (setting%Output%DataOut%isFroudeNumberOut)   N_OutTypeElem =  N_OutTypeElem + 1
        if (setting%Output%DataOut%isHeadOut)           N_OutTypeElem =  N_OutTypeElem + 1
        if (setting%Output%DataOut%isHydRadiusOut)      N_OutTypeElem =  N_OutTypeElem + 1
        if (setting%Output%DataOut%isPerimeterOut)      N_OutTypeElem =  N_OutTypeElem + 1
        if (setting%Output%DataOut%isSlotWidthOut)      N_OutTypeElem =  N_OutTypeElem + 1
        if (setting%Output%DataOut%isSlotDepthOut)      N_OutTypeElem =  N_OutTypeElem + 1
        if (setting%Output%DataOut%isTopWidthOut)       N_OutTypeElem =  N_OutTypeElem + 1
        if (setting%Output%DataOut%isVelocityOut)       N_OutTypeElem =  N_OutTypeElem + 1
        if (setting%Output%DataOut%isVolumeOut)         N_OutTypeElem =  N_OutTypeElem + 1
        if (setting%Output%DataOut%isWaveSpeedOut)      N_OutTypeElem =  N_OutTypeElem + 1

        if (N_OutTypeElem == 0) then
            !% --- if no outputtypes are specified, then suppress the element output
            setting%Output%ElementsExist_global = .false.
            write(*,*) '***********************************************************'
            write(*,*) '** WARNING: Output requested on links or nodes, but the  **'
            write(*,*) '** setting%Output%DataOut%... values did not include any **'
            write(*,*) '** valid data types for elements. Run is proceeding      **'
            write(*,*) '** without any  output data for finite-volume elements   **'
            write(*,*) '***********************************************************'
            return
        end if

        !% --- allocate space for the vector of indexes
        call util_allocate_outputML_elemtypes ()

        !% --- store the true element column indexes
        ii = 0
        !% --- Area
        if (setting%Output%DataOut%isAreaOut) then
            ii = ii+1
            output_types_elemR(ii) = er_Area
            output_typenames_elemR(ii) = 'Area'
            output_typeUnits_elemR(ii) = 'm^2'
            output_typeProcessing_elemR(ii) = AverageElements
        end if
        !% --- Depth
        if (setting%Output%DataOut%isDepthOut) then
            ii = ii+1
            output_types_elemR(ii) = er_Depth
            output_typenames_elemR(ii) = 'Depth'
            output_typeUnits_elemR(ii) = 'm'
            output_typeProcessing_elemR(ii) = AverageElements
        end if
        !% --- Flow rate
        if (setting%Output%DataOut%isFlowrateOut) then
            ii = ii+1
            output_types_elemR(ii) = er_Flowrate
            output_typenames_elemR(ii) = 'Flowrate'
            output_typeUnits_elemR(ii) = 'm^3/s'
            output_typeProcessing_elemR(ii) = AverageElements
        end if
        !% --- Conservative Flux rates, on elements this is the lateral flows
        if (setting%Output%DataOut%isFluxConsOut) then
            ii = ii+1
            output_types_elemR(ii) = er_FlowrateLateral
            output_typenames_elemR(ii) = 'FlowrateLateral'
            output_typeUnits_elemR(ii) = 'm^3/s'
            output_typeProcessing_elemR(ii) = SumElements
        end if
        !% --- Froude Number
        if (setting%Output%DataOut%isFroudeNumberOut) then
            ii = ii+1
            output_types_elemR(ii) = er_FroudeNumber
            output_typenames_elemR(ii) = 'FroudeNumber'
            output_typeUnits_elemR(ii) = 'unitless'
            output_typeProcessing_elemR(ii) = MaximumValue
        end if
        !% --- Head
        if (setting%Output%DataOut%isHeadOut) then
            ii = ii+1
            output_types_elemR(ii) = er_Head
            output_typenames_elemR(ii) = 'PiezometricHead'
            output_typeUnits_elemR(ii) = 'm'
            output_typeProcessing_elemR(ii) = AverageElements
            setting%Output%ElemHeadIndex = ii
        end if
        !% --- HydRadius
        if (setting%Output%DataOut%isHydRadiusOut) then
            ii = ii+1
            output_types_elemR(ii) = er_HydRadius
            output_typenames_elemR(ii) = 'HydraulicRadius'
            output_typeUnits_elemR(ii) = 'm'
            output_typeProcessing_elemR(ii) = AverageElements
        end if
        !% --- Perimeter
        if (setting%Output%DataOut%isPerimeterOut) then
            ii = ii+1
            output_types_elemR(ii) = er_Perimeter
            output_typenames_elemR(ii) = 'WettedPerimeter'
            output_typeUnits_elemR(ii) = 'm'
            output_typeProcessing_elemR(ii) = AverageElements
        end if
        !% --- SlotWidth
        if (setting%Output%DataOut%isSlotWidthOut) then
            ii = ii+1
            output_types_elemR(ii) = er_SlotWidth
            output_typenames_elemR(ii) = 'PreissmannSlotWidth'
            output_typeUnits_elemR(ii) = 'm'
            output_typeProcessing_elemR(ii) = AverageElements
        end if
        !% --- SlotDepth
        if (setting%Output%DataOut%isSlotDepthOut) then
            ii = ii+1
            output_types_elemR(ii) = er_SlotDepth
            output_typenames_elemR(ii) = 'PreissmannSlotDepth'
            output_typeUnits_elemR(ii) = 'm'
            output_typeProcessing_elemR(ii) = AverageElements
        end if
        !% --- TopWidth
        if (setting%Output%DataOut%isTopWidthOut) then
            ii = ii+1
            output_types_elemR(ii) = er_TopWidth
            output_typenames_elemR(ii) = 'FreeSurfaceTopWidth'
            output_typeUnits_elemR(ii) = 'm'
            output_typeProcessing_elemR(ii) = AverageElements
        end if
        !% --- Velocity
        if (setting%Output%DataOut%isVelocityOut) then
            ii = ii+1
            output_types_elemR(ii) = er_Velocity
            output_typenames_elemR(ii) = 'Velocity'
            output_typeUnits_elemR(ii) = 'm/s'
            output_typeProcessing_elemR(ii) = AverageElements
        end if
        !% --- Volume
        if (setting%Output%DataOut%isVolumeOut) then
            ii = ii+1
            output_types_elemR(ii) = er_Volume
            output_typenames_elemR(ii) = 'Volume'
            output_typeUnits_elemR(ii) = 'm^3'
            output_typeProcessing_elemR(ii) = SumElements
        end if
        !% --- WaveSpeed
        if (setting%Output%DataOut%isWaveSpeedOut) then
            ii = ii+1
            output_types_elemR(ii) = er_WaveSpeed
            output_typenames_elemR(ii) = 'EffectiveWaveSpeed'
            output_typeUnits_elemR(ii) = 'm/s'
            output_typeProcessing_elemR(ii) = AverageElements
        end if

        !% -- store 'time' for use in output
        output_typeNames_withTime_elemR(2:ii+1) = output_typeNames_elemR(:)
        output_typeNames_withTime_elemR(1) = 'Time'

        output_typeUnits_withTime_elemR(2:ii+1) = output_typeUnits_elemR(:)
        output_typeUnits_withTime_elemR(1) = 'seconds'

        !% --- set the type of time units for the output
        select case (setting%Output%Report%TimeUnits)
            case (InSeconds)
                output_typeUnits_withtime_elemR(1) = 'seconds'
            case (InMinutes)
                output_typeUnits_withtime_elemR(1) = 'minutes'
            case (InHours)
                output_typeUnits_withtime_elemR(1) = 'hours'
            case (InDays)
                output_typeUnits_withtime_elemR(1) = 'days'
            case default
                write(*,'(A)') 'ERROR (code, user) unknown value forsetting.Output.Report.TimeUnits of...'
                write(*,*) setting%Output%Report%TimeUnits
                stop
        end select

        !%------------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine outputML_element_outtype_selection
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine outputML_face_outtype_selection ()
        !%------------------------------------------------------------------
        !% Description
        !% initializes the types of output data provided for faces
        !% Note that this must be done on every image, including those that
        !% don't have any output faces.
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: ii
            character(64)        :: subroutine_name = 'output_face_outtype_selection'
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (setting%Output%Report%suppress_MultiLevel_Output) return
            if (.not. setting%Output%FacesExist_global) return
            if (setting%Debug%File%output) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        N_OutTypeFace = 0
        !% ---This is tuned to the type of data that exist at faces
        !% count the number of true elements
        !% --- There are 2 areas at a face
        if (setting%Output%DataOut%isAreaOut)           N_OutTypeFace =  N_OutTypeFace + 2

        !% --- Depth does not exist at a face
        !if (setting%Output%DataOut%isDepthOut)         N_OutTypeFace =  N_OutTypeFace + 1

        !% --- Flowrate
        if (setting%Output%DataOut%isFlowrateOut)       N_OutTypeFace =  N_OutTypeFace + 1

        !% --- Conservative fluxes
        if (setting%Output%DataOut%isFluxConsOut)       N_OutTypeFace =  N_OutTypeFace + 1

        !% --- there is no Froude number computed at a face
        !if (setting%Output%DataOut%isFroudeNumberOut)   N_OutTypeFace =  N_OutTypeFace + 1

        !% --- There are 2 Heads at a face
        if (setting%Output%DataOut%isHeadOut)           N_OutTypeFace =  N_OutTypeFace + 2

        !% HydRadius, Perimeter, SlotWidth, and SlotDepth do not exist at a face
        !if (setting%Output%DataOut%isHydRadiusOut)      N_OutTypeFace =  N_OutTypeFace + 1
        !if (setting%Output%DataOut%isPerimeterOut)      N_OutTypeFace =  N_OutTypeFace + 1
        !if (setting%Output%DataOut%isSlotWidthOut)      N_OutTypeFace =  N_OutTypeFace + 1
        !if (setting%Output%DataOut%isSlotDepthOut)      N_OutTypeFace =  N_OutTypeFace + 1

         !% --- There are 2 TopWidths at a face
        if (setting%Output%DataOut%isTopWidthOut)       N_OutTypeFace =  N_OutTypeFace + 2

        !% --- there are 2 velocities at a face
        if (setting%Output%DataOut%isVelocityOut)       N_OutTypeFace =  N_OutTypeFace + 2

        !% --- Volume does not exist at a face
        !if (setting%Output%DataOut%isVolumeOut)         N_OutTypeFace =  N_OutTypeFace + 1

        !% --- Wave speed does not exist at a face
       ! if (setting%Output%DataOut%isWaveSpeedOut)      N_OutTypeFace =  N_OutTypeFace + 1

        if (N_OutTypeFace == 0) then
            !% --- if no outputtypes are specified, then suppress the face output
            setting%Output%ElementsExist_global = .false.
            write(*,*) '***********************************************************'
            write(*,*) '** WARNING: Output requested on links or nodes, but the  **'
            write(*,*) '** setting%Output%DataOut%... values did not include any **'
            write(*,*) '** valid data types for faces. Run is proceeding without **'
            write(*,*) '** any utput data for finite-volume faces                **'
            write(*,*) '***********************************************************'
            return
        endif

        !% allocate space for the vector of indexes
        call util_allocate_outputML_facetypes ()

        !% store the true element column indexes
        ii = 0
        !% --- Area
        if (setting%Output%DataOut%isAreaOut) then
            ii = ii+1
            output_types_faceR(ii) = fr_Area_u
            output_typeNames_faceR(ii) = 'AreaUpstream'
            output_typeUnits_faceR(ii) = 'm^2'
            output_typeProcessing_faceR(ii) = SingleValue
            ii = ii+1
            output_types_faceR(ii) = fr_Area_d
            output_typeNames_faceR(ii) = 'AreaDownstream'
            output_typeUnits_faceR(ii) = 'm^2'
            output_typeProcessing_faceR(ii) = SingleValue
        end if
        !% --- Flowrate
        if (setting%Output%DataOut%isFlowrateOut) then
            ii = ii+1
            output_types_faceR(ii) = fr_Flowrate
            output_typeNames_faceR(ii) = 'Flowrate'
            output_typeUnits_faceR(ii) = 'm^3/s'
            output_typeProcessing_faceR(ii) = SingleValue
        end if
        !% --- Conservative Fluxes
        if (setting%Output%DataOut%isFluxConsOut) then
            ii = ii+1
            output_types_faceR(ii) = fr_Flowrate_Conservative
            output_typeNames_faceR(ii) = 'FlowrateConservative'
            output_typeUnits_faceR(ii) = 'm^3/s'
            output_typeProcessing_faceR(ii) = SingleValue
        end if
        !% --- Head
        if (setting%Output%DataOut%isHeadOut) then
            ii = ii+1
            output_types_faceR(ii) = fr_Head_u
            output_typeNames_faceR(ii) = 'PiezometricHeadUpstream'
            output_typeUnits_faceR(ii) = 'm'
            output_typeProcessing_faceR(ii) = SingleValue
            setting%Output%FaceUpHeadIndex = ii
            ii = ii+1
            output_types_faceR(ii) = fr_Head_d
            output_typeNames_faceR(ii) = 'PiezometricHeadDownstream'
            output_typeUnits_faceR(ii) = 'm'
            output_typeProcessing_faceR(ii) = SingleValue
            setting%Output%FaceDnHeadIndex = ii
        end if
        !% --- TopWidth
        if (setting%Output%DataOut%isTopWidthOut) then
            ii = ii+1
            output_types_faceR(ii) = fr_TopWidth_u
            output_typeNames_faceR(ii) = 'FreeSurfaceTopWidthUpstream'
            output_typeUnits_faceR(ii) = 'm'
            output_typeProcessing_faceR(ii) = SingleValue
            ii = ii+1
            output_types_faceR(ii) = fr_TopWidth_d
            output_typeNames_faceR(ii) = 'FreeSurfaceTopWidthDownstream'
            output_typeUnits_faceR(ii) = 'm'
            output_typeProcessing_faceR(ii) = SingleValue
        end if
        !% -- Velocity
        if (setting%Output%DataOut%isVelocityOut) then
            ii = ii+1
            output_types_faceR(ii) = fr_Velocity_u
            output_typeNames_faceR(ii) = 'VelocityUpstream'
            output_typeUnits_faceR(ii) = 'm/s'
            output_typeProcessing_faceR(ii) = SingleValue
            ii = ii+1
            output_types_faceR(ii) = fr_Velocity_d
            output_typeNames_faceR(ii) = 'VelocityDownstream'
            output_typeUnits_faceR(ii) = 'm/s'
            output_typeProcessing_faceR(ii) = SingleValue
        end if

        !% --- store 'time' for use in output
        output_typeNames_withtime_faceR(2:ii+1) = output_typenames_faceR(:)
        output_typeNames_withtime_faceR(1) = 'Time'

        output_typeUnits_withTime_faceR(2:ii+1) = output_typeUnits_faceR(:)
        output_typeUnits_withTime_faceR(1) = 'seconds'

        !% --- set the type of time units for the output
        select case (setting%Output%Report%TimeUnits)
            case (InSeconds)
                output_typeUnits_withtime_faceR(1) = 'seconds'
            case (InMinutes)
                output_typeUnits_withtime_faceR(1) = 'minutes'
            case (InHours)
                output_typeUnits_withtime_faceR(1) = 'hours'
            case (InDays)
                output_typeUnits_withtime_faceR(1) = 'days'
            case default
                write(*,'(A)') 'ERROR (code, user) unknown value setting.Output.Report.TimeUnits of...'
                write(*,*) setting%Output%Report%TimeUnits
                stop
        end select

        !%------------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine outputML_face_outtype_selection
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine outputML_store_data (isLastStep)
        !%------------------------------------------------------------------
        !% Description -- stores the multi-level output data in memory
        !% Hiearchy of data for outpuot
        !% --- elemR(:,:) is working data at the current time level on each image (coarray)
        !% --- elemOutR(:,:,:) is multi-time-level storage data of elemR for each image (coarray)
        !% --- OutElemDataR(:,:,:) combines the elements from different images into a global set (not coarray)
        !% --- OutElemFixedI(:,:) is static integer data (e.g. indexes) that are needed for output (not coarray)
        !% --- OutLink_ElemDataR(:,:,:,:) organizes OutElemDataR() by links, and add time as a data type (not coarray)
        !% --- OutLink_ProcessedDataR(:,:,:) averages the OutLink_ElemDataR over the links (not coarray)
        !%-------------------------------------------------------------------
            logical, intent(in) :: isLastStep
            integer, pointer :: thisLevel, Npack, thisP(:), thisType(:)
            character(64)    :: subroutine_name = 'outputML_store_data'
        !%--------------------------------------------------------------------
        !% Preliminaries
            if (setting%Output%Report%suppress_MultiLevel_Output) return
            if (setting%Debug%File%output) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%--------------------------------------------------------------------
        !% Aliases
            thisLevel => setting%Output%LastLevel
        !%--------------------------------------------------------------------
        !% --- increment the stored time level counter
        setting%Output%LastLevel              = setting%Output%LastLevel +1
        setting%Output%NumberOfTimeLevelSaved = setting%Output%NumberOfTimeLevelSaved +1
        
        !% --- store the time for this data
        output_times(thisLevel) = setting%Time%Now

        !% --- null the element storage (array exists in all images)
        elemOutR(:,:,thisLevel) = nullvalueR

        !% --- store the element data
        if (setting%Output%ElementsExist_byImage(this_image())) then
            !% --- get the pack size of output elements
            Npack => npack_elemP(ep_Output_Elements)
            !% --- set of output elements
            thisP => elemP(1:Npack,ep_Output_Elements)
            !% --- set of output types
            thisType => output_types_elemR(:)
            !% --- vector store
            elemOutR(1:Npack,:,thisLevel) = elemR(thisP,thisType)
            !% --- correct head to report in same reference base as the *.inp file
            if (setting%Solver%SubtractReferenceHead) then
                if (setting%Output%ElemHeadIndex > 0) then
                    elemOutR(1:Npack,setting%Output%ElemHeadIndex,thisLevel)  &
                  = elemOutR(1:Npack,setting%Output%ElemHeadIndex,thisLevel)  &
                  + setting%Solver%ReferenceHead
                end if
            end if
            !% no action
        end if

        !% null the storage  (array exists in all images)
        faceOutR(:,:,thisLevel) = nullvalueR 

        !% --- store the face data
        if (setting%Output%FacesExist_byImage(this_image())) then
            !% --- get the pack size of faces
            Npack => npack_faceP(fp_Output_Faces)
            !% --- set of output faces
            thisP => faceP(1:Npack,fp_Output_Faces)
            !% --- set of output types
            thisType => output_types_faceR(:)
            !% --- vector store
            faceOutR(1:Npack,:,thisLevel) = faceR(thisP,thisType)
            !% --- correct head to report in same reference base as the *.inp file
            if (setting%Solver%SubtractReferenceHead) then
                if (setting%Output%FaceUpHeadIndex > 0) then
                    faceOutR(1:Npack,setting%Output%FaceUpHeadIndex,thisLevel)  &
                  = faceOutR(1:Npack,setting%Output%FaceUpHeadIndex,thisLevel)  &
                  + setting%Solver%ReferenceHead
                end if
                if (setting%Output%FaceDnHeadIndex > 0) then
                    faceOutR(1:Npack,setting%Output%FaceDnHeadIndex,thisLevel)  &
                  = faceOutR(1:Npack,setting%Output%FaceDnHeadIndex,thisLevel)  &
                  + setting%Solver%ReferenceHead
                end if
            end if
        else
            !% no action
        end if

        !if (setting%Output%Verbose) write(*,"(A,i5)") &
        !    '**************************************************** Storing Time Level #',thisLevel

        !% --- if storage limit is reached, combine the output, write to file,
        !% --- and reset the storage time level
        if ((thisLevel == setting%Output%StoredLevels) .or. (isLastStep)) then
            setting%Output%NumberOfWriteSteps = setting%Output%NumberOfWriteSteps+1
            call outputML_combine_and_write_data (thisLevel, isLastStep)
            thisLevel = 0
        end if

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%output) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine outputML_store_data
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine outputML_combine_and_write_data (nLevel,isLastStep)
        !%------------------------------------------------------------------
        !% Description:
        !% combine all the stored output data
        !% This operates only on image 1, but cycles through each of the
        !% images to collect the data for output
        !%------------------------------------------------------------------
        !% Declarations:
        !%    OpenCoarrays has problems with array-based transfers, which 
        !%    is sometimes slow. Use .false. here unless there are problems 
        !%    causing seg faults in the output.
            logical :: OpenCoarraysMethod = .false.

            integer, intent(in) :: nLevel      !% number of time levels to be written in this file
            logical, intent(in) :: isLastStep  !% is this the last step of model run
            integer :: nTotalElem , nTotalFace !% total number of elements/faces to be written
            integer :: nTypeElem  , nTypeFace  !% total number of element/face types to be written
            integer, pointer ::  nWritten   !% number of files written (# of calls to this procedure)
            integer, pointer ::  nTotalTimeLevels  !% sum of all the time levels written in all files

            integer, pointer :: thiSE(:), thisF(:)

            integer :: nMaxElem, nMaxFace

            integer ::  ios
            integer :: Lasti, firstIdx, lastIdx, npack, ii,  kk, mm, pp
            integer ::  thisUnit  !% file unit numbers

            integer :: dimvector(3)    !% used to store (nTotalElem, nType, nLevel)
            character(len=256) :: file_name
            character(len=5) :: thisnum
            character(len=16) :: str_header
            character(64)    :: subroutine_name = 'outputML_combine_and_write_data'
        !%------------------------------------------------------------------
        !% Preliminaries
            !%=============================
            if (this_image() .ne. 1) return !% --- only run serial
            !%=============================
            if (setting%Output%Report%suppress_MultiLevel_Output) return
            if (setting%Debug%File%output) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases
            nMaxElem = maxval(N_OutElem) !% not an alias, but needed here
            thisE => thisElementOut(1:nMaxElem) !% temporary output element indexes

            nMaxFace = maxval(N_OutFace)  !% not an alias, but needed here
            thisF => thisFaceOut(1:nMaxFace) !% temporary output face indexes

            !% --- point to latest data on number of files written and time levels writte
            nWritten         => setting%File%outputML_Ncombined_file_written
            nTotalTimeLevels => setting%File%outputML_total_timelevels_written
        !%------------------------------------------------------------------
        !% --------------------------------------
        !% --- ELEMENT INDEX DATA
        !%
        nTotalElem = sum(N_OutElem(:))
        nTypeElem  = size(output_types_elemR)

        !% null the index array
        thisE(:) = nullvalueI

        !% --- combine the data from all images
        Lasti = 0 !% counter of the last element or face that has been stored
        !% ---
        do ii = 1,num_images()

            !% cycle if no output elements on this image
            if (.not. setting%Output%ElementsExist_byImage(ii)) cycle

            !% --- get the number of packed LOCAL element indexes
            npack = npack_elemP(ep_Output_Elements)[ii]

            !if (npack == 0) cycle !% if no elements to output, Lasti is unchanged

            !% --- store the coarray elements on this image
            thisE(1:npack) = elemP(1:npack,ep_Output_Elements)[ii]

            !% --- store the GLOBAL indexes for elements, SWMM links and SWMM nodes
            !%     Open Coarrays doesn't support these type of array statements
            ! if (OpenCoarraysMethod) then
                do kk=Lasti+1,Lasti+npack
                    OutElemFixedI(kk,oefi_elem_Gidx)      = elemI(thisE(kk-Lasti),ei_Gidx)[ii]
                    OutElemFixedI(kk,oefi_link_Gidx_SWMM) = elemI(thisE(kk-Lasti),ei_link_Gidx_SWMM)[ii]
                    OutElemFixedI(kk,oefi_node_Gidx_SWMM) = elemI(thisE(kk-Lasti),ei_node_Gidx_SWMM)[ii]
                end do !%kk
            ! else
                !% the following section gives seg fault on large systems, even with only a 
                !% single processor
            !     print *, 'C ',size(OutElemFixedI,DIM=1), size(OutElemFixedI,DIM=2)
            !     print *, '  ',size(elemI,DIM=1), size(elemI,DIM=2)
            !     print *, '  ',size(thisE), npack
            !     print *, '  ',Lasti, npack
            !     print *, '  ' ,size(OutElemFixedI,DIM=2),oefi_elem_Gidx,oefi_link_Gidx_SWMM,oefi_node_Gidx_SWMM
            !     print *, '  ' , size(elemI,DIM=2),ei_Gidx, ei_link_Gidx_SWMM,ei_node_Gidx_SWMM
            !     print *, thisE(1:npack)
            !     OutElemFixedI(Lasti+1:Lasti+npack,oefi_elem_Gidx)      = elemI(thisE(1:npack),ei_Gidx)[ii]
            !     print *, 'here A'
            !     OutElemFixedI(Lasti+1:Lasti+npack,oefi_link_Gidx_SWMM) = elemI(thisE(1:npack),ei_link_Gidx_SWMM)[ii]
            !     print *, 'here B'
            !     OutElemFixedI(Lasti+1:Lasti+npack,oefi_node_Gidx_SWMM) = elemI(thisE(1:npack),ei_node_Gidx_SWMM)[ii]
            !     print *, 'here C'
            ! end if

            !% --- store the real data
            !%     Open Coarrays doesn't support these type of array statements
            if (OpenCoarraysMethod) then
                do pp=1,nLevel
                    do mm=1,nTypeElem
                        do kk=Lasti+1,Lasti+npack
                            OutElemDataR(kk,mm,pp) =  elemOutR(kk-Lasti,mm,pp)[ii]
                        end do !% kk
                    end do !% mm
                end do !% pp
            else
                OutElemDataR(Lasti+1:Lasti+npack,1:nTypeElem,1:nLevel) &
                 =  elemOutR(      1:npack      ,1:nTypeElem,1:nLevel)[ii]
            end if
            !% increment Lasti for the next image
            Lasti = Lasti + npack
        end do !% images

        !% --------------------------------------
        !% --- FACE INDEX DATA
        !%
        nTotalFace = sum(N_OutFace(:))
        nTypeFace = size(output_types_faceR)

        !% null the index array
        thisF(:) = nullvalueI

        Lasti = 0 !% counter of the last element or face that has been stored
        do ii = 1,num_images()

            if (.not. setting%Output%FacesExist_byImage(ii)) cycle

            !% --- get the packed LOCAL face indexes
            npack = npack_faceP(fp_Output_Faces)[ii]

            !if (npack == 0) cycle !% if no elements to output, Lasti is unchanged

            !% --- store the coarray elements on image(1)
            thisF(1:npack) = faceP(1:npack,fp_Output_Faces)[ii]

            !% --- store the GLOBAL indexes for faces, SWMM links and SWMM nodes
            !%     Open Coarrays doesn't support these type of array statements
            ! if (OpenCoarraysMethod) then
                do kk=Lasti+1,Lasti+npack
                    OutFaceFixedI(kk,offi_face_Gidx)      = faceI(thisF(kk-Lasti),fi_Gidx)[ii]
                    OutFaceFixedI(kk,offi_node_Gidx_SWMM) = faceI(thisF(kk-Lasti),fi_node_idx_SWMM)[ii]
                end do !% kk
            ! else
            !     OutFaceFixedI(Lasti+1:Lasti+npack,offi_face_Gidx)      = faceI(thisF(1:npack),fi_Gidx)[ii]
            !     OutFaceFixedI(Lasti+1:Lasti+npack,offi_node_Gidx_SWMM) = faceI(thisF(1:npack),fi_node_idx_SWMM)[ii]
            ! end if


            !% --- store the real data
            !%     Open Coarrays doesn't support these type of array statements
            if (OpenCoarraysMethod) then
                do pp=1,nLevel
                    do mm=1,nTypeFace
                        do kk=Lasti+1,Lasti+npack
                            OutFaceDataR(kk,mm,pp) =  faceOutR(kk-Lasti,mm,pp)[ii]
                        end do !% kk
                    end do !% mm
                end do !% pp
            else
                OutFaceDataR(Lasti+1:Lasti+npack            ,1:nTypeFace,1:nLevel) &
                 =  faceOutR(Lasti+1-Lasti:Lasti+npack-Lasti,1:nTypeFace,1:nLevel)[ii]
            end if

            !% --- increment Lasti for the next image
            Lasti = Lasti + npack
        end do !% images()

        !% --- keeping track of number of files written
        !% --- begins at 0 so nWritten=1 is when writing first file
        nWritten = nWritten + 1
        !% --- HACK the following limits number of files -- fails if nWritten > 99999
        !% --- Need a method for (1) larger string (2)) warning if the number
        !% --- of files gets large and (3) soft stop when max is hit
        if (nWritten > 99999) then
            write(*,"(A)") 'ERROR (code, user): the intermediate output files have reached ...'
            write(*,"(A,I6,A)") '...the code the limit of',nWritten,' files...'
            write(*,"(A)") 'Suggest increasing setting.Output.StoredLevels.'
            stop
        end if

        !% -----------------------------------------------
        !% --- WRITING DATA
        !%
        write(thisnum,"(I5.5)") nWritten

        !% --- file name to write data to
        file_name = trim(setting%File%outputML_combinedfile_kernel) &
            //'_'//thisnum//'.bin'
        !% --- store the filename for later use
        call outputML_store_binary_output_filenames (nWritten, file_name)

        !% --- open unformatted file for data writing
        open(newunit=thisUnit, file=trim(file_name), form='unformatted', &
                action='write',status='new', iostat = ios)
        if (ios /= 0) then
            write(*,"(A)") 'ERROR (CODE) file could not be opened for writing...'
            write(*,"(A)") 'filename is ...'
            write(*,"(A)") trim(file_name)
            stop
        end if

        !% ----------------------------------
        !% --- TIME LEVELS ARE COMMON TO FACES AND ELEMENT OUTPUT
        !% --- the time levels
        write(thisUnit) nLevel
        write(thisUnit) output_times(1:nLevel)
        nTotalTimeLevels = nTotalTimeLevels+nLevel

        !% ----------------------------------
        !% --- WRITING ELEMENTS
        !% --- dimension vector
        dimvector(1) = nTotalElem
        dimvector(2) = nTypeElem
        dimvector(3) = nLevel

        if (nTotalElem > 0) then
            !% --- the output types
            write(thisUnit) nTypeElem
            write(thisUnit) output_types_elemR(1:nTypeElem)

            !% --- the fixed integer data
            write(thisUnit) (/ nTotalElem, Ncol_oefi/)
            write(thisUnit) OutElemFixedI(:,:)

            !% --- the combined elemR output array
            write(thisUnit) dimvector(:)
            write(thisUnit) OutElemDataR(1:nTotalElem,1:nTypeElem,1:nLevel)
        end if

        !% ----------------------------------
        !% --- WRITING FACES
        !% --- dimension vector
        dimvector(1) = nTotalFace
        dimvector(2) = nTypeFace
        dimvector(3) = nLevel

        if (nTotalFace > 0) then
            !% --- the output face types
            write(thisUnit) nTypeFace
            write(thisUnit) output_types_faceR(1:nTypeFace)

            !% --- the fixed integer data
            write(thisUnit) (/ nTotalFace, Ncol_offi/)
            write(thisUnit) OutFaceFixedI(:,:)

            !% --- the combined elemR output array
            write(thisUnit) dimvector(:)
            write(thisUnit) OutFaceDataR(1:nTotalFace,1:nTypeFace,1:nLevel)
        end if

        if (setting%Output%Verbose) write(*,"(A,i5)") &
            '************************************************* finished writing file #', nWritten

        !% --- close the unformatted unit
        close(thisUnit)

        !% --- note that fnunit for setting%File%outputML_filename_file remains open
        !% --- this is so that subsequent calls can write to it.

        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine outputML_combine_and_write_data
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine outputML_write_control_file ()
        !%------------------------------------------------------------------
        !% Description
        !% stores the global data to a file so that the outputML_convert_elements_to_linknode_and_write
        !% can be made independent of the run
        !%------------------------------------------------------------------
            integer, pointer :: nTotalTimeLevels  !% sum of all the time levels written in all files
            integer          :: nTotalElem            !% total number of elements to be written
            integer          :: nTypeElem             !% total number of element types to be written (not including time)

            integer, pointer :: thisunit
            integer          :: ios

            character(len=256) :: file_name
            character(64) :: subroutine_name = 'outputML_write_control_file'
        !%------------------------------------------------------------------
        !% Preliminaries:
            !%=============================
            if (this_image() .ne. 1) return !% --- only run serial
            !%=============================
            if (setting%Output%Report%suppress_MultiLevel_Output) return
            if (setting%Debug%File%output) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% --- open the control file
        thisunit => setting%File%UnitNumber%outputML_control_file
        file_name =       trim(setting%File%outputML_control_file)

        open(unit=thisUnit, file=trim(file_name), form='unformatted', &
                action='write',status='new', iostat = ios)
        if (ios /= 0) then
            write(*,"(A)") 'ERROR (CODE) file could not be opened for writing...'
            write(*,"(A)") 'filename is ...'
            write(*,"(A)") trim(file_name)
            stop
        end if

        !% --- BEGIN WRITING
        !% --- number of files written (# of calls to outputML_combine_and_write_data)
        write(thisUnit) setting%File%outputML_Ncombined_file_written
        !% --- number of time levels written
        write(thisUnit) setting%File%outputML_total_timelevels_written
        !% --- the model starting time
        write(thisUnit) setting%Time%StartEpoch
        !% --- integer key for report time units (e.g. IsHours)
        write(thisUnit) setting%Output%Report%TimeUnits
        !% --- total number number of output elements
        write(thisUnit) sum(N_OutElem(:))
        !% --- number of output element types (excluding time)
        write(thisUnit) size(output_types_elemR)
         !% --- total number number of output faces
        write(thisUnit) sum(N_OutFace(:))
        !% --- number of output face types (excluding time)
        write(thisUnit) size(output_types_faceR)
        !% --- maximum number of stored levels in a file
        write(thisUnit) setting%Output%StoredLevels

        close(thisunit)

        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine outputML_write_control_file
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine outputML_convert_elements_to_linknode_and_write ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% uses stored binary files or multi-level memory to create and write linknode
        !% data for output
        !%
        !% Understanding Arrays
        !% OutElem... leading size is the total number of output elements (NtotalOutputElements)
        !% OutLink...  leading size is total number of links with output
        !% OutNodeElem...  leading size is total number of nodes from elements
        !% OutNodeFace...   leading size is total number of nodes from faces
        !% SWMMlink... leading size is total number of SWMM links (SWMM_N_link)
        !% SWMMnode... leading size is total nuber of SWMM nodes (SWMM_N_node)
        !%
        !% The compliexity of this is due to the need to separate node reporting that
        !% is associated with a finite volume element (i.e., a junction) -- which is the
        !% set "NodeElem" -- from node reporting from finite volume faces (i.e., a
        !% junction between two links) -- which is the set "NodeFace'
        !%
        !% DEVELOPER NOTE: Our goal is for this subroutine to get all of its data from
        !% the files that it reads (i.e. no globals used) so that it is easy to make
        !% this a stand-along program at a later date.
        !% To make this stand alone, we will need to write some of the global arrays
        !% to the control file, then allocate new space here when we read this in.
        !% These include: output_types_elemR, output_times, elemI(:,ei_link_Gidx_SWMM).
        !% Also for independence, below the global OutElemDataR(:,:,:),
        !% which are used in the ML write, needs to be allocated as local null array.
        !% This doesn't need to be written in the control array, we can simply create
        !% locally rather than use from global.
        !%-----------------------------------------------------------------------------
        !% Declarations
            integer :: nWritten
            integer :: nTotalTimeLevels
            integer, pointer   :: thiselem, thislink, thisface, thisnode, thisType
            integer, pointer   :: SWMMlink, SWMMnode
            integer, pointer   :: swmmIdx(:), tlink(:), tnode(:)
            integer :: lasttimestart, lasttimeread !% last timelevel started, read and processed
            integer :: ii, kk, mm, pp,  mminc, ios, allocation_status, fu
            integer :: npackElem     !% number of element items in a pack
            integer :: npackFace     !% number of face items in a pack
            integer :: nPackVolume   !% number of output types of volume

            integer :: nTypeElem, nTypeFace     !% number of data types
            integer :: nTypeElemWtime, nTypeFaceWtime !% number of data types with time included
            integer :: nTotalElem, oldnTotalElem  !% total number of output elements (and prior value)
            integer :: nTotalFace, oldnTotalFace  !% total number of output faces (and prior value)
            integer :: nLevel    !% number of time levels
            integer :: nMax_elemInLink     !% maximum number of elements in any output link
            integer :: nMax_elemInNodeElem ! maximum number of elements in any output node
            integer :: nMax_faceInNodeFace ! maximum number of faces in any output node
            integer :: nOutLink  !% total number of output links
            integer :: nOutNodeElem  !% total number of output nodes from elements
            integer :: nOutNodeFace  !% total number of output nodes from faces
            integer :: nOutElemFixedColumns !% number of columns in the OutElem_FixedI(:,:) array
            integer :: nOutFaceFixedColumns !% number of columns in the OutFace_FixedI(:,:) array
            integer :: dimvector(3), olddimvector(3)  !% size of 3D array (and prior value)
            integer :: additional_rows  !% number of additional rows in link and node arrays due to Bquick

            !integer, allocatable :: SWMMlink_num_elements(:) !% number of elements in each output link
            !integer, allocatable :: SWMMnode_num_elements(:) !% number of elements in each output link
            
            integer, allocatable :: OutLink_N_elem_in_link(:)  !% number of elements in links
            integer, allocatable :: OutNodeElem_N_elem_in_node(:)  !% number of elements in nodes
            integer, allocatable :: OutNodeFace_N_face_in_node(:) !% number of faces in nodes

            integer, allocatable, target :: OutLink_pSWMMidx(:)       !% Global link index packed for output link size
            integer, allocatable, target :: OutNodeElem_pSWMMidx(:)   !% Global node index packed for output node/elem size
            integer, allocatable, target :: OutNodeFace_pSWMMidx(:)           !% Global node index packed for output node/face size

            !integer, allocatable         :: OutElem_LinearIdx(:)        !% 1..n indexes corresponding to 1:nTotalElem for packing
            !integer, allocatable         :: OutFace_LinearIdx(:)        !% 1..n indexes corresponding to 1:nTotalElem for packing
            integer, allocatable         :: pOutLinkElem(:)             !% local packed locations of elements that are links
            integer, allocatable         :: pOutNodeElem(:)             !% local packed locations of elements that are nodes
            integer, allocatable         :: pOutNodeFace(:)              !% local packed locations o faces that are nodes

            !% full list of packed Elem ID for each output link (link, list of elemID) note the valid length of columns is SWMMlink_num_elements(kk)
            integer, allocatable, target :: OutLink_pOutElemIdx(:,:) !
            integer, allocatable, target :: OutNodeElem_pOutElemIdx(:,:)        !% packed locations of 1:nTotalElem with elements that are nodes
            integer, allocatable, target :: OutNodeFace_pOutFaceIdx(:,:)   !% packed locations of OutFaceGidx with faces that are nodes

            real(8), allocatable         :: OutLink_ElemDataR(:,:,:,:) !% (link, elements in link, data types, time levels)
            real(8), allocatable         :: OutNodeElem_ElemDataR(:,:,:,:) !% (node, elements in node, data types, time levels)
            real(8), allocatable         :: OutNodeFace_FaceDataR(:,:,:,:) !% (node, faces in node, data types, time levels)

            real(8), allocatable         :: OutLink_ProcessedDataR(:,:,:)       !%  (link, data types, time levels
            real(8), allocatable         :: OutNodeElem_ProcessedDataR(:,:,:)   !%  (node, data types, time levels)
            real(8), allocatable         :: OutNodeFace_ProcessedDataR(:,:,:)   !%  (node, data types, time levels)

            !% aliases
            integer, pointer             :: pOutElem_Gidx(:)
            integer, pointer             :: pOutElem_Link_SWMM_idx(:)
            integer, pointer             :: pOutElem_Node_SWMM_idx(:)

            integer, pointer             :: pOutFace_Gidx(:)
            integer, pointer             :: pOutFace_Node_SWMM_idx(:)

            integer, pointer             :: pElem(:), pFace(:)

            ! logical to suppress face file writing of processed file (used for junction to avoid processing faces)
            logical, allocatable         :: isOutNodeFaceWriteFVonly(:)
            logical, allocatable         :: isOutNodeElemWriteFVOnly(:)
            logical, allocatable         :: isOutLinkWriteFVOnly(:)

            integer :: rlimits(2)  ! reshaping array

            logical :: isopen = .false.

            integer            :: deallocation_status
            integer            :: thisUnit
            character(len=256) :: thisFile
            character(len=32)  :: tlinkname, tnodename
            character(len=256) :: fn_link_unf, fn_link_csv, fn_linkFV_csv
            character(len=256) :: fn_nodeElem_unf, fn_nodeElem_csv, fn_nodeElemFV_csv
            character(len=256) :: fn_nodeFace_unf, fn_nodeFace_csv, fn_nodeFaceFV_csv
            integer            :: fU_link_unf, fU_link_csv, fU_linkFV_csv
            integer            :: fU_nodeElem_unf, fU_nodeElem_csv, fU_nodeElemFV_csv
            integer            :: fU_nodeFace_unf, fU_nodeFace_csv, fU_nodeFaceFV_csv
            character(len=99)  :: emsg
            character(len=8)   :: tstatus

            character(len=16)  :: time_units_str

            real(8) :: StartTimeEpoch
            integer :: reportTimeUnits
            integer :: NtotalOutputElements
            integer :: NtotalOutputFaces
            integer :: StoredLevels
            logical :: verbose

            integer :: dummyarrayI(1) = 1
            integer :: dummyI = 1

            real(8) :: time_secs, time_epoch, time_scale_for_output
            integer :: startdate(6) !% yr, month, day, hr, min, sec
            character(64)      :: subroutine_name = 'outputML_convert_elements_to_linknode_and_write'
        !%-----------------------------------------------------------------------------
        !% Preliminaries
            !%=============================
            if (this_image() .ne. 1) return !% --- only run serial
            !%=============================
            if (setting%Output%Report%suppress_MultiLevel_Output) return
            if (setting%Debug%File%output) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"               
            !% --- HACK -- need to change this to an input variable if making this independent
            verbose = setting%Output%Verbose
            if (verbose) write(*,*) '**** Beginning unformatted file conversion and multi-level output processing ****'
        !% -------------------------------------------------------------
        !% --- CONTROL FILE
        !% --- open and read the control file
        thisFile = trim(setting%File%outputML_control_file)

        open(newunit=thisUnit, file=trim(thisFile), form='unformatted', &
                action='read',status='old', iostat = ios)
        if (ios /= 0) then
            write(*,"(A)") 'ERROR (CODE) file could not be opened for writing...'
            write(*,"(A)") 'filename is ...'
            write(*,"(A)") trim(thisFile)
            stop 309870
        end if

        !% --- get the total number of combined files written
        read(thisUnit) nWritten
        !% --- get the total number of time levels written
        read(thisUnit) nTotalTimeLevels
        !% --- get the Epoch(days) start time
        read(thisUnit) StartTimeEpoch
        !% --- get the integer key for report time units (e.g., IsHours)
        read(thisUnit) reportTimeUnits
        !% --- get the total number of elements written in all linlks
        read(thisUnit) NtotalOutputElements
        !% --- get the total number of element data types written (excluding time)
        read(thisUnit) nTypeElem
        !% --- get the total number offaces written in all linlks
        read(thisUnit) NtotalOutputFaces
        !% --- get the total number of data types written (excluding time)
        read(thisUnit) nTypeFace
        !% --- get the maximum number of stored levels in a file
        read(thisUnit) StoredLevels

        ! print *, nWritten, 'nWritten'
        ! print *, nTotalTimeLevels,'nTotalTimeLevels'
        ! print *, StartTimeEpoch, 'StartTimeEpoch'
        ! print *, reportTimeUnits, 'reportTimeUnits'
        ! print *, NtotalOutputElements, 'NtotalOutputElements'
        ! print *, nTypeElem, 'nTypeElem'
        ! print *, NtotalOutputFaces, 'NtotalOutputFaces'
        ! print *, nTypeFace, 'nTypeFace'
        ! print *, StoredLevels, 'StoredLevels'

        !% --- close the control file
        close(thisUnit)

        !% --- get the start time as "thisd
        call util_datetime_decodedate(StartTimeEpoch, startdate(1), startdate(2), startdate(3))
        call util_datetime_decodetime(StartTimeEpoch, startdate(4), startdate(5), startdate(6))

        !% --- set the type of output time units
        select case (reportTimeUnits)
        case (InSeconds)
            time_units_str = 'seconds'
            time_scale_for_output = oneR
        case (InMinutes)
            time_units_str = 'minutes'
            time_scale_for_output = seconds_per_minute
        case (InHours)
            time_units_str = 'hours'
            time_scale_for_output = seconds_per_hour
        case (InDays)
            time_units_str = 'days'
            time_scale_for_output = seconds_per_day
        case default
            write(*,'(A)') 'ERROR (code) unknown value of setting%Output%Report%TimeUnits of ...'
            write(*,*) setting%Output%Report%TimeUnits
            stop
        end select

        !% --- HACK to make this independent of globals, this call will have to be changed and files always written/read.
        call outputML_get_all_output_binary_filenames (nWritten)

        lasttimeread = 0
        lasttimestart = 0

        !% -----------------------------------------
        !% --- OUTER LOOP CYCLES THROUGH ALL THE MULTI-LEVEL FILES THAT HAVE BEEN WRITTEN
        do ii=1,nWritten
            !% -----------------------------------
            !% --- PART 1 --- ALLOCATING ELEMENT AND FACE ARRAYS AND READING DATA
            !% -----------------------------------
                thisFile = trim(output_binary_filenames_all(ii))

                if (verbose) write(*,"(A)") 'reading file :',trim(thisFile)

                open(newunit=thisUnit, file=thisFile, form='unformatted', &
                    action='read', iostat=ios)
                if (ios .ne. 0) then
                    write(*,"(A)") 'ERROR (file): iostat /=0 for open() file named ...'
                    write(*,"(A)") trim(thisFile)
                    write(*,"(A)") '... file is an unformated file of output data...'
                    write(*,"(A,i5)") '... iostat value = ',ios
                    stop
                end if

                !% -------------------------------
                !% --- read and store the time levels
                read(thisUnit) nLevel
                lasttimestart = lasttimeread + 1
                lasttimeread  = lasttimeread+nLevel
                if (nLevel .le. size(output_times)) then
                    output_times(:) = nullvalueR
                    read(thisUnit) output_times(1:nLevel)
                else
                    write(*,"(A)") 'ERROR (code, file): unexpected array size problem...'
                    write(*,"(A,i5)") '...output_times has size ',shape(output_times)
                    write(*,"(A,i5)") '...but needs to read in ',nLevel
                    stop
                end if

                !% -------------------------------------------
                !% --- BELOW HERE FOR ELEMENTS
                if (NtotalOutputElements >0 ) then
                    !% -------------------------------
                    !% --- READ THE ELEMENTS
                    olddimvector(1) = NtotalOutputElements
                    olddimvector(2) = NtypeElem
                    olddimvector(3) = StoredLevels

                    !% --- read and store the vector of output types
                    read(thisUnit) nTypeElem
                    if (nTypeElem .le. size(output_types_elemR)) then
                        output_types_elemR(:) = nullvalueI
                        read(thisUnit) output_types_elemR(1:nTypeElem)
                    else
                        write(*,"(A)") 'ERROR (code, file), unexpected array size problem...'
                        write(*,"(A,i5)") '...output_types_elemR has size ',size(output_types_elemR)
                        write(*,"(A,i5)") '...but needs to read in ',nTypeElem
                        stop
                    end if

                    !% --- read and store the fixed integer data
                    read(thisUnit) nTotalElem, nOutElemFixedColumns
                    read(thisUnit) OutElemFixedI(:,:)

                    !% --- set the aliases
                    pOutElem_Gidx          => OutElemFixedI(:,oefi_elem_Gidx)
                    pOutElem_Link_SWMM_idx => OutElemFixedI(:,oefi_link_Gidx_SWMM)
                    pOutElem_Node_SWMM_idx => OutElemFixedI(:,oefi_node_Gidx_SWMM)

                    !% --- read and store the 3D array
                    read(thisUnit) dimvector
                    !% --- error checking
                    if (dimvector(1) .ne. nTotalElem) then
                        write(*,"(A)") 'ERROR (code, file): unexpected array size problem...'
                        write(*,"(A,i5)") '...dimvector(1) has value ',dimvector(1)
                        write(*,"(A,i5)") '...but nTotalElem is ',nTotalElem
                        stop
                    end if
                    if (dimvector(2) .ne. nTypeElem) then
                        write(*,"(A)") 'ERROR (code, file): unexpected array size problem...'
                        write(*,"(A,i5)") '...dimvector(2) has value ',dimvector(2)
                        write(*,"(A,i5)") '...but nTotalElem is ',nTypeElem
                        stop
                    end if
                    !% --- ensure that dimvector(3) is less than or equal to allocated
                    if (dimvector(3) > olddimvector(3)) then
                        write(*,"(A)") 'ERROR (code, file): unexpected array size problem...'
                        write(*,"(A,i5)") '...dimvector(3) has value ',dimvector(3)
                        write(*,"(A,i5)") '...but olddimvector(3) is ',olddimvector(3)
                        stop
                    end if
                    !% --- check the dimvector(3) is consistent with lasttimestart and lasttimeread
                    if (dimvector(3) .ne. lasttimeread + 1 - lasttimestart) then
                        write(*,"(A)") 'ERROR (code, file): unexpected array size problem...'
                        write(*,"(A,i5)") '...dimvector(3) has value ',dimvector(3)
                        write(*,"(A,i5)") '...but lasttimeread is  ',lasttimeread
                        write(*,"(A,i5)") '...and lasttimestart is ',lasttimestart
                        stop
                    end if
                    nTotalElem = dimvector(1)
                    nTypeElem  = dimvector(2)
                    nLevel = dimvector(3)
                    OutElemDataR(:,:,:) = nullvalueR

                    nTypeElemWtime = nTypeElem + 1
                    read(thisUnit) OutElemDataR(1:nTotalElem,1:nTypeElem,1:nLevel)
                end if !% NtotalOutputElements > 0

                !% -------------------------------------------
                !% --- BELOW HERE FOR FACES
                if (NtotalOutputFaces > 0) then
                    !% -------------------------------
                    !% --- READ THE FACES
                    olddimvector(1) = NtotalOutputFaces
                    olddimvector(2) = NtypeFace
                    olddimvector(3) = StoredLevels

                    !% --- read and store the vector of output types
                    read(thisUnit) nTypeFace
                    if (nTypeFace .le. size(output_types_faceR)) then
                        output_types_faceR(:) = nullvalueI
                        read(thisUnit) output_types_faceR(1:nTypeFace)
                    else
                        write(*,"(A)") 'ERROR (code, file), unexpected array size problem...'
                        write(*,"(A,i5)") '...output_types_faceR has size ',size(output_types_faceR)
                        write(*,"(A,i5)") '...but needs to read in ',nTypeFace
                        stop
                    end if

                    !% --- read and store the fixed integer data
                    read(thisUnit) nTotalFace, nOutFaceFixedColumns
                    read(thisUnit) OutFaceFixedI(:,:)

                    !% --- set the aliases
                    pOutFace_Gidx          => OutFaceFixedI(:,offi_face_Gidx)
                    pOutFace_Node_SWMM_idx => OutFaceFixedI(:,offi_node_Gidx_SWMM)

                    !% --- read and store the 3D array
                    read(thisUnit) dimvector
                    !% --- error checking
                    if (dimvector(1) .ne. nTotalFace) then
                        write(*,"(A)") 'ERROR (code, file): unexpected array size problem...'
                        write(*,"(A,i5)") '...dimvector(1) has value ',dimvector(1)
                        write(*,"(A,i5)") '...but nTotalFace is ',nTotalFace
                        stop
                    end if
                    if (dimvector(2) .ne. nTypeFace) then
                        write(*,"(A)") 'ERROR (code, file): unexpected array size problem...'
                        write(*,"(A,i5)") '...dimvector(2) has value ',dimvector(2)
                        write(*,"(A,i5)") '...but nTotalElem is ',nTypeFace
                        stop
                    end if
                    !% --- ensure that dimvector(3) is less than or equal to allocated
                    if (dimvector(3) > olddimvector(3)) then
                        write(*,"(A)") 'ERROR (code, file): unexpected array size problem...'
                        write(*,"(A,i5)") '...dimvector(3) has value ',dimvector(3)
                        write(*,"(A,i5)") '...but olddimvector(3) is ',olddimvector(3)
                        stop
                    end if
                    !% --- check the dimvector(3) is consistent with lasttimestart and lasttimeread
                    if (dimvector(3) .ne. lasttimeread + 1 - lasttimestart) then
                        write(*,"(A)") 'ERROR (code, file): unexpected array size problem...'
                        write(*,"(A,i5)") '...dimvector(3) has value ',dimvector(3)
                        write(*,"(A,i5)") '...but lasttimeread is  ',lasttimeread
                        write(*,"(A,i5)") '...and lasttimestart is ',lasttimestart
                        stop
                    end if
                    nTotalFace = dimvector(1)
                    nTypeFace  = dimvector(2)
                    nLevel = dimvector(3)
                    OutFaceDataR(:,:,:) = nullvalueR

                    nTypeFaceWtime = nTypeFace + 1
                    read(thisUnit) OutFaceDataR(1:nTotalFace,1:nTypeFace,1:nLevel)
                end if !% NtotalOutputFaces > 0
                !% -- done reading this file
                close(thisUnit)

            !% -----------------------------------
            !% --- PART 2a --- COUNT THE NUMBER OF ELEMENTS PER LINK AND ELEMENTS PER NODE
            !% -----------------------------------
                if (NtotalOutputElements > 0) then
                    if (ii==1) then
                        !% Moved to til_allocate_outputML_storage  brh 20220202
                        ! !% --- allocate space for storing the number of elements in each link (including non-output links)
                        ! !% --- note this MUST be the size of the SWMM_N_link so that we can later pack indexes
                        ! allocate(SWMMlink_num_elements(SWMM_N_link), stat=allocation_status, errmsg=emsg)
                        ! call util_allocate_check(allocation_status, emsg, 'SWMMlink_num_elements')
                        ! SWMMlink_num_elements(:) = 0

                        ! !% --- allocate space for storing the number of elements in each node (including non-output nodes)
                        ! !% --- note this MUST be the size of the SWMM_N_node so that we can later pack indexes
                        ! allocate(SWMMnode_num_elements(SWMM_N_node), stat=allocation_status, errmsg=emsg)
                        ! call util_allocate_check(allocation_status, emsg, 'SWMMnode_num_elements')
                        ! SWMMnode_num_elements(:) = 0

                        do kk=1,nTotalElem
                            !% --- get the number of elements in each output link
                            !% --- separate elements that are nodes from links
                            !% --- cycle through the output elements to get a count of links involved
                            !% --- note that there is always 1:1 correspondence between elem:nodes

                            !% --- the SWMM link idx for the kk output element
                            SWMMlink => pOutElem_Link_SWMM_idx(kk)
                            if (SWMMlink == nullvalueI) then !% then this is a node
                                SWMMnode => pOutElem_Node_SWMM_idx(kk)
                                if (SWMMnode == nullvalueI) then
                                    write(*,'(A)') 'ERROR (code) unexpected nullvalue for an output element...'
                                    write(*,'(A)') '... appears to be neither a link nor a node.'
                                    write(*,'(A,i8)') '... kk = ',kk
                                    write(*,'(A,i8)') '... Global Element Index = ',pOutElem_Gidx(kk)
                                    stop
                                end if
                                !% -- store the node index for each of the output elements
                                !OutElem_SWMMnodeIdx(kk) = SWMMnode
                                SWMMnode_num_elements(SWMMnode) = SWMMnode_num_elements(SWMMnode) + 1
                            else
                                !% -- store the link index for each of the output elements
                                !OutElem_SWMMlinkIdx(kk) = SWMMlink
                                SWMMlink_num_elements(SWMMlink) = SWMMlink_num_elements(SWMMlink) + 1
                            end if

                            !print *, 'kk', kk, OutElem_SWMMlinkIdx(kk), pOutElem_Link_SWMM_idx(kk)

                        end do
                    end if !% ii=1
                end if !% NtotalOutputElements > 0

            !print *, 'EEE FacesExist_byImage',setting%Output%FacesExist_byImage
            !% -----------------------------------
            !% --- PART 2b --- COUNT THE NUMBER OF FACES PER NODE
            !% -----------------------------------
                if (NtotalOutputFaces > 0) then
                    if (ii==1) then
                        !% Moved to til_allocate_outputML_storage  brh 20220202
                        ! !% --- allocate space for storing the number of faces in each node (including non-output node)
                        ! !% --- this should be one for all faces. It is allocate to be able to create common output routines
                        ! allocate(SWMMnode_num_faces(SWMM_N_node), stat=allocation_status, errmsg=emsg)
                        ! call util_allocate_check(allocation_status, emsg, 'SWMMnode_num_faces')
                        ! SWMMnode_num_faces(:) = 0

                        do kk=1,nTotalFace
                            !% --- get the number of faces in each output node
                            !% --- cycle through the output faces to get a count of nodes involved
                            !% --- the corresponding node idx for the face
                            SWMMnode => pOutFace_Node_SWMM_idx(kk)

                            if (SWMMnode == nullvalueI) then !% then this is a face not pointing at a node
                                write(*,'(A)') 'ERROR (code) unexpected nullvalue for an output face...'
                                write(*,'(A)') '... appears to be not part of the SWMM node set.'
                                stop
                            else
                                !% -- store the link index for each of the output elements
                                !OutFace_SWMMnodeIdx(kk) = SWMMnode
                                SWMMnode_num_faces(SWMMnode) = SWMMnode_num_faces(SWMMnode) + 1
                            end if
                        end do
                    end if !% ii=1
                end if !% NtotalOutputFaces > 0

            !print *, 'FFF ElementsExist_byImage', setting%Output%ElementsExist_byImage
            !% -----------------------------------
            !% --- PART 3a --- STORAGE ELEM->LINK CONVERSION
            !% -----------------------------------
                if (NtotalOutputElements > 0) then
                    !% finding the number of additional rows added to the link and node
                    !% arrays to accomodate phantom links and nodes
                    if (setting%Partitioning%PartitioningMethod == BQuick) then
                        additional_rows = num_images() - 1
                    else
                        additional_rows = 0
                    end if
                    !% --- only on first pass through with first file
                    if (ii==1) then
                        !% --- the maximum number of elements in any link
                        nMax_elemInLink = maxval(SWMMlink_num_elements)
                        !% -- count the number of links that have non-zero elements
                        nOutLink = count(SWMMlink_num_elements > 0)

                        !% code error check
                        if (SWMM_N_link .ne. (size(link%I(:,li_idx))-additional_rows)) then
                            write(*,"(A)") 'ERROR (code): we assumed size of SWMM_N_link and size of li_idx are identical...'
                            write(*,"(A)") '... they are not, which is a mismatch for the output. Need code rewrite.'
                            write(*,"(A)") '... SWMM_N_link is ',SWMM_N_link
                            write(*,"(A)") ',... size(link%I(:,li_idx)) is ',(size(link%I(:,li_idx))-additional_rows)
                            ! stop
                        end if

                        if (nOutLink > 0) then

                            !% HACK -- these needs to be moved into an allocation that occurs during initialization
                            !% However, this requires that some of the above computations of faces and elements per
                            !% node must be done during initialization as well.  brh 20220102

                            !% --- create the packed index list of SWMM output links
                            allocate(OutLink_pSWMMidx(nOutLink), stat=allocation_status, errmsg=emsg)
                            call util_allocate_check(allocation_status, emsg, 'OutLink_pSWMMidx')
                            OutLink_pSWMMidx(:) = pack( (/ (mm, mm=1,SWMM_N_link) /) ,(SWMMlink_num_elements > 0))

                            !% --- SPACE FOR DATA
                            !% --- Create a space to store all the elem data for each link and the set of time levels
                            !% --- note we use nTypeElem+1 so that we can also store the time as a datatype
                            allocate(OutLink_ElemDataR(nOutLink,nMax_elemInLink,nTypeElem+1,nLevel), &
                                stat=allocation_status, errmsg=emsg)
                            call util_allocate_check(allocation_status, emsg, 'OutLink_ElemDataR')
                            OutLink_ElemDataR(:,:,:,:) = nullValueR

                            !% -- create storage space for processed link data
                            !% --- note we use nTypeElem+1 so that we can also store the time as a data type
                            allocate(OutLink_ProcessedDataR(nOutLink,nTypeElem+1,nLevel), stat=allocation_status, errmsg=emsg)
                            call util_allocate_check(allocation_status, emsg, 'OutLink_ProcessedDataR')
                            OutLink_ProcessedDataR(:,:,:) = nullValueI

                            !% --- SPACE FOR MAPS
                            !% --- allocate storage of packed element indexes for each link
                            allocate(OutLink_pOutElemIdx(nOutLink,nMax_elemInLink), &
                                stat=allocation_status, errmsg=emsg)
                            call util_allocate_check(allocation_status, emsg, 'OutLink_pOutElemIdx')
                            OutLink_pOutElemIdx(:,:) = nullvalueI

                            !% --- space for logical that limits printing of multi-face nodes to FV output
                            allocate(isOutLinkWriteFVonly(nOutLink), stat=allocation_status, errmsg=emsg)
                            call util_allocate_check(allocation_status, emsg, 'isOutLinkWriteFVonly')
                            isOutLinkWriteFVonly(:) = .false.

                            !% --- create storage for the number of elements in each output link
                            allocate(OutLink_N_elem_in_link(nOutLink), stat=allocation_status, errmsg=emsg)
                            call util_allocate_check(allocation_status, emsg, 'OutLink_N_elem_in_link')
                            !% set the packed number of elements to the value previously computed in the global array
                            OutLink_N_elem_in_link = SWMMlink_num_elements(OutLink_pSWMMidx)

                        end if   !% nOutlink > 0
                    end if ! ii=1
                end if !% NtotalOutputElements > 0

            !print *, 'GGG ElementsExist_byImage', setting%Output%ElementsExist_byImage
            !% -----------------------------------
            !% --- PART 3b --- STORAGE FOR ELEM->NODE CONVERSION
            !% -----------------------------------
                if (NtotalOutputElements > 0) then
                    !% allocation only at the first file read
                    if (ii==1) then
                        !% --- the maximum number of elements in any node/elem
                        nMax_elemInNodeElem = maxval(SWMMnode_num_elements)
                        nOutNodeElem = count(SWMMnode_num_elements > 0)

                        !% code error check
                        if (SWMM_N_node .ne. (size(node%I(:,ni_idx))-additional_rows)) then
                            write(*,"(A)") 'ERROR (code): we assumed size of SWMM_N_node and size of ni_idx are identical...'
                            write(*,"(A)") '... they are not, which is a mismatch for the output. Need code rewrite.'
                            write(*,"(A)") '... SWMM_N_node is ',SWMM_N_node
                            write(*,"(A)") ',... size(node%I(:,ni_idx)) is ',(size(node%I(:,ni_idx))-additional_rows)
                            ! stop
                        end if

                        if (nOutNodeElem > 0) then

                            !% --- create the packed index list of SWMM output nodes
                            allocate(OutNodeElem_pSWMMidx(nOutNodeElem), stat=allocation_status, errmsg=emsg)
                            call util_allocate_check(allocation_status, emsg, 'OutNodeElem_pSWMMidx')
                            OutNodeElem_pSWMMidx(:) = pack( (/ (mm, mm=1,SWMM_N_node)/),(SWMMnode_num_elements > 0))

                            !% --- SPACE FOR DATA
                            !% --- allocate storage of packed node indexes for each link
                            allocate(OutNodeElem_ElemDataR(nOutNodeElem,nMax_elemInNodeElem,1:nTypeElem+1,nLevel), &
                                stat=allocation_status, errmsg=emsg)
                            call util_allocate_check(allocation_status, emsg, 'OutNodeElem_ElemDataR')
                            OutNodeElem_ElemDataR(:,:,:,:) = nullValueR

                            !% -- create storage space for processed node data
                            !% --- note we use nTypeElem+1 so that we can also store the time as a data type
                            allocate(OutNodeElem_ProcessedDataR(nOutNodeElem,nTypeElem+1,nLevel), &
                                stat=allocation_status, errmsg=emsg)
                            call util_allocate_check(allocation_status, emsg, 'OutNodeElem_ProcessedDataR')
                            OutNodeElem_ProcessedDataR(:,:,:) = nullValueR

                            !% --- SPACE FOR MAPS
                            !% --- allocate the element indexes for node indexes in 1:nMax of nTotalNode
                            allocate(OutNodeElem_pOutElemIdx(nOutNodeElem,nMax_elemInNodeElem), &
                                stat=allocation_status, errmsg=emsg)
                            call util_allocate_check(allocation_status, emsg, 'OutNodeElem_pOutElemIdx')
                            OutNodeElem_pOutElemIdx(:,:) = nullvalueI

                            !% --- space for logical that limits printing of multi-face nodes to FV output
                            allocate(isOutNodeElemWriteFVonly(nOutNodeElem), stat=allocation_status, errmsg=emsg)
                            call util_allocate_check(allocation_status, emsg, 'isOutNodeElemWriteFVonly')
                            isOutNodeElemWriteFVonly(:) = .false.

                            !% --- create storage for the number of elements in each output  ode
                            allocate(OutNodeElem_N_elem_in_node(nOutNodeElem), stat=allocation_status, errmsg=emsg)
                            call util_allocate_check(allocation_status, emsg, 'OutNodeElem_N_elem_in_node')
                            OutNodeElem_N_elem_in_node = SWMMnode_num_elements(OutNodeElem_pSWMMidx)

                        end if !% nOutNodeElem > 0
                    end if !% ii=1
                end if !% NtotalOutputElements > 0

                !print *, 'HHH %FacesExist_byImage', setting%Output%FacesExist_byImage
            !% -----------------------------------
            !% --- PART 3c --- STORAGE FOR NODE->FACE CONVERSION
            !% -----------------------------------
                if (NtotalOutputFaces > 0) then
                    !% --- only on first pass through with first file
                    if (ii==1) then
                        !% --- the maximum number of faces in any node
                        nMax_faceInNodeFace = maxval(SWMMnode_num_faces)
                        !% -- count the number of faces that have non-zero elements
                        nOutNodeFace = count(SWMMnode_num_faces > 0)

                        !% code error check
                        if (SWMM_N_node .ne. (size(node%I(:,ni_idx))- additional_rows)) then
                            write(*,"(A)") 'ERROR (code): we assumed size of SWMM_N_node and size of ni_idx are identical...'
                            write(*,"(A)") '... they are not, which is a mismatch for the output. Need code rewrite.'
                            write(*,"(A)") '... SWMM_N_node is ',SWMM_N_node
                            write(*,"(A)") ',... size(node%I(:,ni_idx)) is ',(size(node%I(:,ni_idx))-additional_rows)
                            ! stop
                        end if

                        if (nOutNodeFace > 0) then

                            !% --- create the packed index list of SWMM output nodes
                            allocate(OutNodeFace_pSWMMidx(nOutNodeFace), stat=allocation_status, errmsg=emsg)
                            call util_allocate_check(allocation_status, emsg, 'OutNodeFace_pSWMMidx')
                            OutNodeFace_pSWMMidx(:) = pack((/ (mm,mm=1,SWMM_N_node)/),(SWMMnode_num_faces > 0))

                            !% --- SPACE FOR DATA
                            !% --- Create a space to store all the face data for each node and the set of time levels
                            !% --- note we use nTypeFace+1 so that we can also store the time as a datatype
                            allocate(OutNodeFace_FaceDataR(nOutNodeFace,nMax_faceInNodeFace,nTypeFace+1,nLevel), &
                                stat=allocation_status, errmsg=emsg)
                            call util_allocate_check(allocation_status, emsg, 'OutNodeFace_FaceDataR')
                            OutNodeFace_FaceDataR(:,:,:,:) = nullValueR

                            !% -- create storage space for processing of node data
                            !% --- note we use nTypeFace+1 so that we can also store the time as a data type
                            allocate(OutNodeFace_ProcessedDataR(nOutNodeFace,nTypeFace+1,nLevel), &
                                stat=allocation_status, errmsg=emsg)
                            call util_allocate_check(allocation_status, emsg, 'OutNodeFace_ProcessedDataR')
                            OutNodeFace_ProcessedDataR(:,:,:) = nullValueR

                            !% --- SPACE FOR MAPS
                            !% --- allocate storage of packed node indexes for each node
                            allocate(OutNodeFace_pOutFaceIdx(nOutNodeFace,nMax_faceInNodeFace), &
                                stat=allocation_status, errmsg=emsg)
                            call util_allocate_check(allocation_status, emsg, 'OutNodeFace_pOutFaceIdx')
                            OutNodeFace_pOutFaceIdx(:,:) = nullvalueI

                            !% --- space for logical that limits printing of multi-face nodes to FV output
                            allocate(isOutNodeFaceWriteFVonly(nOutNodeFace), stat=allocation_status, errmsg=emsg)
                            call util_allocate_check(allocation_status, emsg, 'isOutNodeFaceWriteFVonly')
                            isOutNodeFaceWriteFVonly(:) = .false.

                            !% --- create storage for the number of faces in each output node
                            allocate(OutNodeFace_N_face_in_node(nOutNodeFace), stat=allocation_status, errmsg=emsg)
                            call util_allocate_check(allocation_status, emsg, 'OutNodeFace_N_face_in_node')
                            !% set the packed number of faces to the value previously computed in the global array
                            OutNodeFace_N_face_in_node(:) = SWMMnode_num_faces(OutNodeFace_pSWMMidx)

                        end if  !% nOutNodeFace > 0
                    end if ! ii=1
                end if !% NtotalOutputFaces > 0

                !print *, 'III nOutLink ElementsExist_byImage ',setting%Output%ElementsExist_byImage
            !% -----------------------------------
            !% --- PART 4a --- PERFORM ELEM->LINK CONVERSION
            !% -----------------------------------
                !% HACK -- consider making the first element in the arrays 1:nLevel and the last
                !%  as the link index -- this should be faster. But the change could create
                !%  lots of bugs.

                if ( NtotalOutputElements > 0) then
                    do kk=1,nOutLink
                        !% --- Each link must be handled separately because they each
                        !% --- have different numbers of elements.

                        !% --- get the global swmm link index for this kk
                        SWMMlink => OutLink_pSWMMidx(kk)

                        !% --- get the element indexes that match this link
                        npackElem = count(pOutElem_Link_SWMM_idx == SWMMlink)

                        !% pack the OutElem Indexes for the elements in a link
                        OutLink_pOutElemIdx(kk,1:npackElem) &
                            = pack((/ (mm, mm=1,nTotalElem) /), pOutElem_Link_SWMM_idx == SWMMlink)

                        !% --- select the current portion of the pack for use in storage
                        pElem => Outlink_pOutElemIdx(kk,1:npackElem)

                        !% --- store the element data by link (start with nTypeElem=2 to save space for time)
                        OutLink_ElemDataR  (kk,1:npackElem,2:nTypeElem+1,1:nLevel) &
                             = OutElemDataR(pElem         ,1:nTypeElem  ,1:nLevel)

                        !% --- if there is only one element in the link, then store that value for all types
                        if (npackElem == 1) then
                             OutLink_ProcessedDataR(kk,  2:nTypeElemWtime,1:nLevel) &
                                = OutLink_ElemDataR(kk,1,2:nTypeElemWtime,1:nLevel)
                        else
                            !% --- cycle through the data types for different processing (e.g. average, sum, max)
                            !% --- reshape() seems necessary to remove singleton dimensions and sum without seg fault
                            rlimits = (/ npackElem, nLevel/) !% limits for reshaping
                            do pp=2,nTypeElemWtime !% starts at 2 to skip the time column
                                select case (output_typeProcessing_elemR(pp-1)) !% -1 needed for correct index excluding time
                                    case (AverageElements)
                                        !% --- first sum the elements
                                        OutLink_ProcessedDataR(kk,pp,1:nLevel) = &
                                            sum(reshape(OutLink_ElemDataR(kk,1:npackElem,pp,1:nLevel),rlimits(1:2)),dim=1)
                                        !% --- then divide by number of elements in link
                                        OutLink_ProcessedDataR(kk,pp,1:nLevel) = &
                                            OutLink_ProcessedDataR(kk,pp,1:nLevel) / OutLink_N_elem_in_link(kk)
                                    case (SumElements)
                                        !% --- simple sum of the elements in a link
                                        OutLink_ProcessedDataR(kk,pp,1:nLevel) = &
                                            sum(reshape(OutLink_ElemDataR(kk,1:npackElem,pp,1:nLevel),rlimits(1:2)),dim=1)
                                    case (MaximumValue)
                                        !% --- get the maximum value in the link
                                        OutLink_ProcessedDataR(kk,pp,1:nLevel) = &
                                            maxval(reshape(OutLink_ElemDataR(kk,1:npackElem,pp,1:nLevel),rlimits(1:2)),dim=1)
                                    case (SingleValue)
                                        if (npackElem > 1) then
                                            !% if there is more than one value and a single value is desired, only
                                            !% print the FV file
                                            isOutLinkWriteFVonly(kk) = .true.
                                        else
                                            OutLink_ProcessedDataR(kk,pp,1:nLevel) &
                                                = OutLink_ElemDataR(kk,1,pp,1:nLevel)
                                        end if
                                    case default
                                        write(*,'(A)') 'ERROR (code) unknown key index for output_typeProcessing_elemR of ...'
                                        write(*,*) output_typeProcessing_elemR(pp)
                                        stop
                                end select
                            end do
                        end if

                        !% --- add in the time levels (in user-selected time unit)
                        do mm=1,npackElem
                            OutLink_ElemDataR(kk,mm,1,1:nLevel) = output_times(1:nLevel) /  time_scale_for_output
                        end do  !% mm
                        OutLink_ProcessedDataR(kk,1,1:nLevel)   = output_times(1:nLevel) / time_scale_for_output
                    end do !% kk
                end if !% NtotalOutputElements > 0

                !print *, 'JJJ nOutNodeElem ElementsExist_byImage',setting%Output%ElementsExist_byImage
            !% -----------------------------------
            !% --- PART 4b --- PERFORM ELEM->NODE CONVERSION
            !% -----------------------------------
                if (NtotalOutputElements > 0) then
                    do kk=1,nOutNodeElem
                        !% --- Each node must be handled separately because they each
                        !% --- have different numbers of elements.

                        !% --- get the global swmm node index for this kk
                        SWMMnode => OutNodeElem_pSWMMidx(kk)

                        !% --- get the element indexes in SWMM IDx that are for this node
                        npackElem = count(pOutElem_Node_SWMM_idx  == SWMMnode)

                        !% --- pack the global node indexes used for output header
                        OutNodeElem_pOutElemIdx(kk,1:npackElem) &
                            = pack((/ (mm,mm=1,nTotalElem)/), pOutElem_Node_SWMM_idx == SWMMnode)

                        !% --- select the current portion of the pack for use in storage
                        pElem => OutNodeElem_pOutElemIdx(kk,1:npackElem)

                        !% --- store the element data by node (start with nTypeElem=2 to save space for time)
                        OutNodeElem_ElemDataR(kk,1:npackElem,2:nTypeElem+1,1:nLevel) &
                            =    OutElemDataR(pElem         ,1:nTypeElem  ,1:nLevel)

                        !% --- if there is only one element in the node, then store that value for all types
                        if (npackElem == 1) then
                            OutNodeElem_ProcessedDataR (kk  ,2:nTypeElemWtime,1:nLevel)  &
                                = OutNodeElem_ElemDataR(kk,1,2:nTypeElemWtime,1:nLevel)
                        else
                            !% --- cycle through the data types for different processing (e.g. average, sum, max)
                            !% --- reshape() seems necessary to remove singleton dimensions and sum without seg fault
                            rlimits = (/ npackElem, nLevel/) !% limits for reshaping
                            do pp=2,nTypeElemWtime !% starts at 2 to skip the time column
                                select case (output_typeProcessing_elemR(pp-1)) !% -1 needed for correct index excluding time
                                    case (AverageElements)
                                        !% --- first sum the elements
                                        OutNodeElem_ProcessedDataR(kk,pp,1:nLevel) = &
                                            sum(reshape(OutNodeElem_ElemDataR(kk,1:npackElem,pp,1:nLevel),rlimits(1:2)),dim=1)
                                        !% --- then divide by number of elements in link
                                        OutNodeElem_ProcessedDataR(kk,pp,1:nLevel) = &
                                            OutNodeElem_ProcessedDataR(kk,pp,1:nLevel) / OutNodeElem_N_elem_in_node(kk)
                                    case (SumElements)
                                        !% --- simple sum of the elements in a link
                                        OutNodeElem_ProcessedDataR(kk,pp,1:nLevel) = &
                                            sum(reshape(OutNodeElem_ElemDataR(kk,1:npackElem,pp,1:nLevel),rlimits(1:2)),dim=1)
                                    case (MaximumValue)
                                        !% --- get the maximum value in the link
                                        OutNodeElem_ProcessedDataR(kk,pp,1:nLevel) = &
                                            maxval(reshape(OutNodeElem_ElemDataR(kk,1:npackElem,pp,1:nLevel),rlimits(1:2)),dim=1)
                                    case (SingleValue)
                                        if (npackElem > 1) then
                                            !% if there is more than one value and a single value is desired, only
                                            !% print the FV file
                                            isOutNodeElemWriteFVonly(kk) = .true.
                                        else
                                            OutNodeElem_ProcessedDataR(kk,pp,1:nLevel) &
                                                = OutNodeElem_ElemDataR(kk,1,pp,1:nLevel)
                                        end if
                                    case default
                                        write(*,'(A)') 'ERROR (code) unknown key index for output_typeProcessing_elemR of ...'
                                        write(*,*) output_typeProcessing_elemR(pp)
                                        stop
                                end select
                            end do
                        end if
                        !% --- add in the time levels (in user-selected time unit)
                        do mm=1,npackElem
                            OutNodeElem_ElemDataR(kk,mm,1,1:nLevel) = output_times(1:nLevel) /  time_scale_for_output
                        end do !% mm
                        OutNodeElem_ProcessedDataR(kk,1,1:nLevel)   = output_times(1:nLevel) / time_scale_for_output
                    end do !% kk
                end if !% NtotalOutputElements > 0

                !print *, 'KKK FacesExist_byImage (nOutNodeFace)',setting%Output%FacesExist_byImage
            !% -----------------------------------
            !% --- PART 4c --- PERFORM FACE->NODE CONVERSION
            !% -----------------------------------
                if (NtotalOutputFaces > 0) then
                    do kk=1,nOutNodeFace
                        !% --- Each node must be handled separately because they each
                        !% --- have different numbers of faces.

                        !% --- get the global swmm node index for this kk
                        SWMMnode => OutNodeFace_pSWMMidx(kk)

                        !% --- get the face indexes that are for this node
                        npackFace = count(pOutFace_Node_SWMM_idx == SWMMnode)

                        !% pack the global face indexes used for output header
                        OutNodeFace_pOutFaceIdx(kk,1:npackFace) &
                            = pack((/ (mm,mm=1,nTotalFace)/), pOutFace_Node_SWMM_idx == SWMMnode)

                        !% --- select the current portion of the pack for use in storage
                        pFace => OutNodeFace_pOutFaceIdx(kk,1:npackFace)

                        !% --- store the face data by node (start with nTypeElem=2 to save space for time)
                        OutNodeFace_FaceDataR(kk,1:npackFace,2:nTypeFace+1,1:nLevel) &
                            =    OutFaceDataR(pFace         ,1:nTypeFace  ,1:nLevel)

                        !% --- a node that has multiple faces is a junction, and only gets FV output
                        !% --- CUSTOM FOR NODEFACE ONLY
                        if (npackFace > 1) then
                            isOutNodeFaceWriteFVonly(kk) = .true.
                        end if

                        !% --- if there is only one face in the node, then store that value for all types
                        if (npackFace == 1) then
                            OutNodeFace_ProcessedDataR (kk  ,2:nTypeFaceWtime,1:nLevel) &
                                = OutNodeFace_FaceDataR(kk,1,2:nTypeFaceWtime,1:nLevel)
                        else
                            !% --- cycle through the data types for different processing (e.g. average, sum, max)
                            !% --- reshape() seems necessary to remove singleton dimensions and sum without seg fault
                            rlimits = (/ npackFace, nLevel/) !% limits for reshaping
                            do pp=2,nTypeFaceWtime !% starts at 2 to skip the time column
                                select case (output_typeProcessing_faceR(pp-1)) !% -1 needed for correct index excluding time
                                    case (AverageElements)
                                        !% --- first sum the elements
                                        OutNodeFace_ProcessedDataR(kk,pp,1:nLevel) = &
                                            sum(reshape(OutNodeFace_FaceDataR(kk,1:npackFace,pp,1:nLevel),rlimits(1:2)),dim=1)
                                        !% --- then divide by number of elements in link
                                        OutNodeFace_ProcessedDataR(kk,pp,1:nLevel) = &
                                            OutNodeFace_ProcessedDataR(kk,pp,1:nLevel) / OutNodeFace_N_face_in_node(kk)
                                    case (SumElements)
                                        !% --- simple sum of the elements in a link
                                        OutNodeFace_ProcessedDataR(kk,pp,1:nLevel) = &
                                            sum(reshape(OutNodeFace_FaceDataR(kk,1:npackFace,pp,1:nLevel),rlimits(1:2)),dim=1)
                                    case (MaximumValue)
                                        !% --- get the maximum value in the link
                                        OutNodeFace_ProcessedDataR(kk,pp,1:nLevel) = &
                                            maxval(reshape(OutNodeFace_FaceDataR(kk,1:npackFace,pp,1:nLevel),rlimits(1:2)),dim=1)
                                    case (SingleValue)
                                        if (npackFace > 1) then
                                            !% if there is more than one value and a single value is desired, only
                                            !% print the FV file
                                            isOutNodeFaceWriteFVonly(kk) = .true.
                                        else
                                            OutNodeFace_ProcessedDataR(kk,pp,1:nLevel) &
                                                = OutNodeFace_FaceDataR(kk,1,pp,1:nLevel)
                                        end if
                                    case default
                                        write(*,'(A)') 'ERROR (code) unknown key index for output_typeProcessing_faceR of ...'
                                        write(*,*) output_typeProcessing_faceR(pp)
                                        stop
                                end select
                            end do
                        end if

                        !% --- add in the time levels (in user-selected time unit)
                        do mm=1,npackFace
                            OutNodeFace_FaceDataR(kk,mm,1,1:nLevel) = output_times(1:nLevel) /  time_scale_for_output
                        end do  !% mm
                        OutNodeFace_ProcessedDataR(kk,1,1:nLevel)   = output_times(1:nLevel) / time_scale_for_output
                    end do !% kk
                end if !% NtotalOutputFaces > 0

                !print *, 'LLL ', nOutLink
            !% -----------------------------------
            !% --- PART 5 --- WRITE TO OUTPUT FILES (open and close each)
            !% -----------------------------------
            !% --- 6 sets of output files (at this point we are not doing finite-volume unformatted)
            !% node_####.unf -> all data and time steps for the node with string ID #### in fortran unformatted
            !% link_####.unf -> all data and time steps for the link with string ID #### in fortran unformatted
            !% node_####.csv -> same as node_####.unf but in csv file
            !% link_####.csv -> same as link_####.unf but in csv file
            !% linkFV_####_xxxx.csv -> data type xxxx finite-volume element-level data for link with string ID #### in csv
            !% nodeFV_####_xxx.csv -. data type xxxx finite-volume node data (elem or face) for node with string ID #### and data type xxxx in csv

            !% --- HACK --- THE UNFORMATTED FILES HAVE NOT BEEN VERIFIED

            !% -----------------------------------
            !% --- PART 5a --- WRITE OUTPUT FOR LINKS
            !% -----------------------------------
                if (NtotalOutputElements > 0) then
                    do kk=1,nOutLink
                        !% --- Cycle through the links to create the individual link output files

                        !% get the global SWMM index for this link
                        SWMMlink => OutLink_pSWMMidx(kk)

                        !% -- error checking
                        !% -- HACK -- need to write a check routine with VERIFY() to make sure link%Names(SWMMlink)%str is a valid string
                        if (.not. allocated(link%Names(SWMMlink)%str)) then
                            write (*,"(A,i8)") 'ERROR (code): link%Name(SWMMlink)%str not allocated for SWMMlink=',SWMMlink
                            stop
                        end if

                        if (len(link%Names(SWMMlink)%str) == 0) then
                            write(*,"(A,i8)") 'ERROR (code)): link%Name(kk)%str is empty for SWMMlink= ',SWMMlink
                            stop
                        end if

                        if (len(link%Names(SWMMlink)%str) > len(tlinkname)) then
                            write(*,"(A)") 'ERROR (user): User link name is too long...'
                            write(*,"(A,i8)") '... in link%Name(kk)%str is too long for SWMMlink= ',SWMMlink
                            write(*,"(A,i8)") '... max length is: ',len(link%Names(SWMMlink)%str)
                            write(*,"(A)") '... link name in SWMM is ...'
                            write(*,"(A)") trim(link%Names(SWMMlink)%str)
                            stop
                        end if

                        !% --- use a temporary name for convenience
                        tlinkname = trim(link%Names(SWMMlink)%str)

                        !% -----------------------------------------
                        !% --- LINK FILES, CSV AND UNF (all types in 1 file)
                        !%
                        if (.not. isOutLinkWriteFVOnly(kk)) then
                            !% --- set the filenames for output of SWMM links

                            fn_link_unf = trim(setting%File%outputML_Link_kernel) // '_' //trim(tlinkname) //'.unf'
                            fn_link_csv = trim(setting%File%outputML_Link_kernel) // '_' //trim(tlinkname) //'.csv'

                            if (ii==1) then  !% --- Create new link output files and write headers for first file read
                                ! !% --- open unformatted link file
                                ! open(newunit=fU_link_unf, file=trim(fn_link_unf), form='unformatted', &
                                !     action='write', access='append', status='new')
                                ! !% --- write header to unformated link file
                                ! call outputML_unf_header( &
                                !     fU_link_unf, nTypeElemWtime, nTotalTimeLevels, &
                                !     startdate, setting%Time%StartEpoch, &
                                !     output_typeNames_withTime_elemR, output_typeUnits_withTime_elemR, &
                                !     dummyarrayI,    &
                                !     tlinkname, setting%Time%DateTimeStamp, time_units_str, .false.)

                                !% --- open formatted csv link file
                                open(newunit=fU_link_csv, file=trim(fn_link_csv), form='formatted', &
                                    action='write', access='append')
                                !% --- write header to csv link file
                                call outputML_csv_header( &
                                    fU_link_csv, nTypeElem, nTotalTimeLevels, dummyI, &
                                    OutLink_pSWMMidx(kk), &
                                    startdate, setting%Time%StartEpoch, &
                                    output_typeNames_withTime_elemR, output_typeUnits_withTime_elemR, &
                                    dummyarrayI,    &
                                    tlinkname, setting%Time%DateTimeStamp, time_units_str, LinkOut, .false.)

                                !% --- finished with the headers
                            else !% --- for ii > 2, link csv and unf files exists so we just need to open
                                ! open(newunit=fU_link_unf, file=trim(fn_link_unf), form='unformatted', &
                                !     action='write', access='append', status='old')
                                open(newunit=fU_link_csv, file=trim(fn_link_csv), form='formatted',  &
                                    action='write', access='append')
                            end if

                            !% --- write link data to the unformatted file
                            ! call outputML_unf_writedata (&
                            !       fU_link_unf, nTypeElemWtime, nTotalTimeLevels, kk, OutLink_ProcessedDataR)
                            !% --- write link data to the csv formatted data for these nLevels
                            call outputML_csv_writedata ( &
                                fU_link_csv, kk, nTypeElemWtime, dummyI, nLevel,  &
                                OutLink_ProcessedDataR, OutLink_ElemDataR, .false. )

                            ! close(fU_link_unf) !% close the unformatted link file
                            close(fU_link_csv) !% close the csv link file
                        end if !% write link files

                        !% -----------------------------------------
                        !% --- LINK FINITE-VOLUME FILES CSV (1 type per file)
                        !%
                        do mm=1,nTypeElem
                            !% --- cycle through the types
                            mminc = mm+1 !% increment to skip time level
                            !% --- create the filename
                            fn_linkFV_csv = trim(setting%File%outputML_Link_kernel) &
                                // 'FV_' //trim(tlinkname) //'_'//trim(output_typeNames_elemR(mm)) //'.csv'

                            if (ii==1) then
                                !% --- open a new file for this type and set the header
                                open(newunit=fU_linkFV_csv, file=trim(fn_linkFV_csv), form='formatted', &
                                    action='write', access='append')
                                !% --- write the header
                                call outputML_csv_header( &
                                    fU_linkFV_csv, OutLink_N_elem_in_link(kk), nLevel, mminc, &
                                    OutLink_pSWMMidx(kk), &
                                    startdate, setting%Time%StartEpoch, &
                                    output_typeNames_withTime_elemR, output_typeUnits_withTime_elemR, &
                                    pOutElem_Gidx(OutLink_pOutElemIdx(kk,1:OutLink_N_elem_in_link(kk))), &
                                    tlinkname, setting%Time%DateTimeStamp, time_units_str, LinkOut, .true. )
                                !% --- finished writing headers
                            else !% --- for ii> 2, open the existing FV file for this type and link
                                open(newunit=fU_linkFV_csv, file=trim(fn_linkFV_csv), form='formatted', &
                                    action='write', position='append')
                            end if
                            !% --- write the csv FV output for elements of kk link with mm type and the latest 1:nLevel
                            call outputML_csv_writedata ( &
                            fU_linkFV_csv, kk, OutLink_N_elem_in_link(kk), mminc, nLevel,  &
                            OutLink_ProcessedDataR, OutLink_ElemDataR, .true.)
                            !% --- close this file so that the unit# can be used for a new file
                            close(fU_linkFV_csv)
                        end do !% mm
                        !% --- finished writing all Link files for link ii
                    end do !% kk
                end if !% NtotalOutputElements > 0

                !print *, 'MMM ', nOutNodeElem
            !% -----------------------------------
            !% --- PART 5b --- WRITE OUTPUT FOR NODES THAT ARE ELEMENTS
            !% -----------------------------------
                if (NtotalOutputElements > 0) then
                    do kk=1,nOutNodeElem
                        !% --- Cycle through the nodes to create the individual nodes output files
                        !% get the global SWMM index for this node
                        SWMMnode => OutNodeElem_pSWMMidx(kk)

                        !% -- error checking
                        !% -- HACK -- need to write a check routine with VERIFY() to make sure link%Names(SWMMlink)%str is a valid string
                        if (.not. allocated(node%Names(SWMMnode)%str)) then
                            write (*,"(A,i8)") 'ERROR (code): node%Name(SWMMnode)%str not allocated for SWMMnode=',SWMMnode
                            stop
                        end if

                        if (len(node%Names(SWMMnode)%str) == 0) then
                            write(*,"(A,i8)") 'ERROR (code)): node%Name(SWMMnode)%str is empty for SWMMnode= ',SWMMnode
                            stop
                        end if

                        if (len(node%Names(SWMMnode)%str) > len(tnodename)) then
                            write(*,"(A)") 'ERROR (user): User node name is too long...'
                            write(*,"(A,i8)") '... in node%Name(SWMMnode)%str is too long for SWMMnode= ',SWMMnode
                            write(*,"(A,i8)") '... max length is: ',len(node%Names(SWMMnode)%str)
                            write(*,"(A)") '... node name in SWMM is ...'
                            write(*,"(A)") trim(node%Names(SWMMnode)%str)
                            stop
                        end if

                        !% --- use a temporary name for convenience
                        tnodename = trim(node%Names(SWMMnode)%str)

                        !% -----------------------------------------
                        !% --- NODE-ELEM FILES, CSV AND UNF (all types in 1 file)
                        !%
                        if (.not. isOutNodeElemWriteFVOnly(kk)) then
                            !% --- set the filenames for output of SWMM links
                            fn_nodeElem_unf = trim(setting%File%outputML_Node_kernel) // '_' //trim(tnodename) //'.unf'
                            fn_nodeElem_csv = trim(setting%File%outputML_Node_kernel) // '_' //trim(tnodename) //'.csv'

                            if (ii==1) then  !% --- Create new node output files and write headers for first file read
                                ! !% --- open unformatted node file
                                ! open(newunit=fU_nodeElem_unf, file=trim(fn_nodeElem_unf), form='unformatted', &
                                !     action='write', access='append', status='new')
                                ! !% --- write header to unformated node file

                                ! call outputML_unf_header( &
                                !     fU_nodeElem_unf, nTypeElemWtime, nTotalTimeLevels, &
                                !     startdate, setting%Time%StartEpoch, &
                                !     output_typeNames_withTime_elemR, output_typeUnits_withTime_elemR, &
                                !     dummyarrayI,    &
                                !     tlinkname, setting%Time%DateTimeStamp, time_units_str, .false.)

                                !% --- open formatted csv node file
                                open(newunit=fU_nodeElem_csv, file=trim(fn_nodeElem_csv), form='formatted', &
                                    action='write', access='append')
                                !% --- write header to csv node file
                                call outputML_csv_header( &
                                    fU_nodeElem_csv, nTypeElem, nTotalTimeLevels, dummyI, &
                                    OutNodeElem_pSWMMidx(kk), &
                                    startdate, setting%Time%StartEpoch, &
                                    output_typeNames_withTime_elemR, output_typeUnits_withTime_elemR, &
                                    pOutElem_Gidx(OutNodeElem_pOutElemIdx(kk,1:OutNodeElem_N_elem_in_node(kk))), &
                                    tnodename, setting%Time%DateTimeStamp, time_units_str, NodeElemOut, .false.)
                                !% --- finished with the headers
                            else
                                ! open(newunit=fU_nodeElem_unf, file=trim(fn_nodeElem_unf), form='unformatted', &
                                !     action='write', access='append', status='old')
                                open(newunit=fU_nodeElem_csv, file=trim(fn_nodeElem_csv), form='formatted',  &
                                    action='write', access='append')
                            end if

                            !% --- write node data to the unformatted file
                            ! call outputML_unf_writedata (&
                            !     fU_nodeElem_unf, nTypeElemWtime, nTotalTimeLevels, kk, &
                            !     OutNodeElem_ProcessedDataR)
                            !% --- write node data to the csv formatted data for these nLevels
                            call outputML_csv_writedata ( &
                                fU_nodeElem_csv, kk, nTypeElemWtime, dummyI, nLevel,  &
                                OutNodeElem_ProcessedDataR, OutNodeElem_ElemDataR, .false. )

                            ! close(fU_nodeElem_unf) !% close the unformatted link file
                            close(fU_nodeElem_csv) !% close the csv link file
                        end if !% write node-elem files

                        !% -----------------------------------------
                        !% --- NODE-ELEM FINITE-VOLUME FILES CSV (1 type per file)
                        !%
                        do mm=1,nTypeElem
                            !% --- cycle through the types
                            mminc = mm+1 !% increment to skip time level
                            !% --- create the filename
                            fn_nodeElemFV_csv = trim(setting%File%outputML_Node_kernel) &
                                // 'FV_' //trim(tnodename) //'_'//trim(output_typeNames_elemR(mm)) //'.csv'

                            if (ii==1) then
                                !% --- open a new file for this type and set the header
                                open(newunit=fU_nodeElemFV_csv, file=trim(fn_nodeElemFV_csv), form='formatted', &
                                    action='write', access='append')
                                !% --- write the header
                                call outputML_csv_header( &
                                    fU_nodeElemFV_csv, OutNodeElem_N_elem_in_node(kk), nLevel, mminc, &
                                    OutNodeElem_pSWMMidx(kk), &
                                    startdate, setting%Time%StartEpoch, &
                                    output_typeNames_withTime_elemR, output_typeUnits_withTime_elemR, &
                                    pOutElem_Gidx(OutNodeElem_pOutElemIdx(kk,1:OutNodeElem_N_elem_in_node(kk))), &
                                    tnodename, setting%Time%DateTimeStamp, time_units_str, NodeElemOut, .true.)
                                !% --- finished writing headers
                            else !% --- for ii> 2, open the existing FV file for this type and node
                                open(newunit=fU_nodeElemFV_csv, file=trim(fn_nodeElemFV_csv), form='formatted', &
                                    action='write', position='append')
                            end if
                            call outputML_csv_writedata ( &
                                fU_nodeElemFV_csv, kk, OutNodeElem_N_elem_in_node(kk), mminc, nLevel,  &
                                OutNodeElem_ProcessedDataR, OutNodeElem_ElemDataR, .true.)
                            !% --- close this file so that the unit# can be used for a new file
                            close(fU_nodeElemFV_csv)
                        end do !% mm
                    end do !% kk
                end if !% NtotalOutputElements > 0
                !% --- finished writing all Node output files for NodeElem

                !print *, 'NNN ', nOutNodeFace
            !% -----------------------------------
            !% --- PART 5c --- WRITE OUTPUT FOR NODES THAT ARE FACES
            !% -----------------------------------
                if (NtotalOutputFaces > 0) then
                    do kk=1,nOutNodeFace
                        !% --- Cycle through the nodes to create the individual nodes output files
                        !% get the global SWMM index for this node
                        SWMMnode => OutNodeFace_pSWMMidx(kk)

                        !% -- error checking
                        !% -- HACK -- need to write a check routine with VERIFY() to make sure link%Names(SWMMnode)%str is a valid string
                        if (.not. allocated(node%Names(SWMMnode)%str)) then
                            write (*,"(A,i8)") 'ERROR (code): node%Name(SWMMnode)%str not allocated for SWMMnode=',SWMMnode
                            stop
                        end if

                        if (len(node%Names(SWMMnode)%str) == 0) then
                            write(*,"(A,i8)") 'ERROR (code)): node%Name(SWMMnode)%str is empty for SWMMnode= ',SWMMnode
                            stop
                        end if

                        if (len(node%Names(SWMMnode)%str) > len(tnodename)) then
                            write(*,"(A)") 'ERROR (user): User node name is too long...'
                            write(*,"(A,i8)") '... in node%Name(SWMMnode)%str is too long for SWMMnode= ',SWMMnode
                            write(*,"(A,i8)") '... max length is: ',len(node%Names(SWMMnode)%str)
                            write(*,"(A)") '... node name in SWMM is ...'
                            write(*,"(A)") trim(node%Names(SWMMnode)%str)
                            stop
                        end if

                        !% --- use a temporary name for convenience
                        tnodename = trim(node%Names(SWMMnode)%str)

                        !% -----------------------------------------
                        !% --- NODE-FACE FILES, CSV AND UNF (all types in 1 file)
                        !%
                        if (.not. isOutNodeFaceWriteFVonly(kk)) then
                            !% --- set the filenames for output of SWMM links
                            fn_nodeFace_unf = trim(setting%File%outputML_Node_kernel) &
                                // '_face_' //trim(tnodename) //'.unf'
                            fn_nodeFace_csv = trim(setting%File%outputML_Node_kernel) &
                                // '_face_' //trim(tnodename) //'.csv'

                            if (ii==1) then  !% --- Create new node output files and write headers for first file read
                                ! !% --- open unformatted node file
                                ! open(newunit=fU_nodeFace_unf, file=trim(fn_nodeFace_unf), form='unformatted', &
                                !     action='write', access='append', status='new')
                                ! !% --- write header to unformated node file
                                ! call outputML_unf_header( &
                                !     fU_nodeFace_unf, nTypeFaceWtime, nTotalTimeLevels, &
                                !     startdate, setting%Time%StartEpoch, &
                                !     output_typeNames_withTime_faceR, output_typeUnits_withTime_faceR, &
                                !     dummyarrayI,    &
                                !     tnodename, setting%Time%DateTimeStamp, time_units_str, .false.)

                                !% --- open formatted csv node file
                                open(newunit=fU_nodeFace_csv, file=trim(fn_nodeFace_csv), form='formatted', &
                                    action='write', access='append')
                                !% --- write header to csv node file
                                call outputML_csv_header( &
                                    fU_nodeFace_csv, nTypeFace, nTotalTimeLevels, dummyI, &
                                    OutNodeFace_pSWMMidx(kk), &
                                    startdate, setting%Time%StartEpoch, &
                                    output_typeNames_withTime_faceR, output_typeUnits_withTime_faceR, &
                                    pOutFace_Gidx(OutNodeFace_pOutFaceIdx(kk,1:OutNodeFace_N_face_in_node(kk))), &
                                    tnodename, setting%Time%DateTimeStamp, time_units_str, NodeFaceOut, .false.)
                                !% --- finished with the headers
                            else
                                !print *,' in here'
                                ! open(newunit=fU_nodeFace_unf, file=trim(fn_nodeFace_unf), form='unformatted', &
                                !     action='write', access='append', status='old')
                                open(newunit=fU_nodeFace_csv, file=trim(fn_nodeFace_csv), form='formatted',  &
                                    action='write', access='append')
                            end if

                            !% --- write node data to the unformatted file
                            ! call outputML_unf_writedata ( &
                            !     fU_nodeFace_unf, nTypeFaceWtime, nTotalTimeLevels, kk, OutNodeFace_ProcessedDataR)

                            !% --- write node data to the csv formatted data for these nLevels
                            call outputML_csv_writedata ( &
                                fU_nodeFace_csv, kk, nTypeFaceWtime, dummyI, nLevel,  &
                                OutNodeFace_ProcessedDataR, OutNodeFace_FaceDataR, .false. )

                            ! close(fU_nodeFace_unf) !% close the unformatted link file
                            close(fU_nodeFace_csv) !% close the csv link file
                        end if

                        !% -----------------------------------------
                        !% --- NODE-FACE FINITE-VOLUME FILES CSV (1 type per file)
                        !%
                        do mm=1,nTypeFace
                            !% --- cycle through the types
                            mminc = mm+1 !% increment to skip time level
                            !% --- create the filename
                            fn_nodeFaceFV_csv = trim(setting%File%outputML_Node_kernel) &
                                // 'FV_face_' //trim(tnodename) //'_'//trim(output_typeNames_faceR(mm)) //'.csv'

                            if (ii==1) then
                                !% --- open a new file for this type and set the header
                                open(newunit=fU_nodeFaceFV_csv, file=trim(fn_nodeFaceFV_csv), form='formatted', &
                                    action='write', access='append')
                                !% --- write the header
                                call outputML_csv_header( &
                                    fU_nodeFaceFV_csv, OutNodeFace_N_face_in_node(kk), nLevel, mminc, &
                                    OutNodeFace_pSWMMidx(kk), &
                                    startdate, setting%Time%StartEpoch, &
                                    output_typeNames_withTime_faceR, output_typeUnits_withTime_faceR, &
                                    pOutFace_Gidx(OutNodeFace_pOutFaceIdx(kk,1:OutNodeFace_N_face_in_node(kk))), &
                                    tnodename, setting%Time%DateTimeStamp, time_units_str, NodeFaceOut, .true. )
                                !% --- finished writing headers
                            else !% --- for ii> 2, open the existing FV file for this type and node
                                open(newunit=fU_nodeFaceFV_csv, file=trim(fn_nodeFaceFV_csv), form='formatted', &
                                    action='write', position='append')
                            end if
                            call outputML_csv_writedata ( &
                                fU_nodeFaceFV_csv, kk, OutNodeFace_N_face_in_node(kk), mminc, nLevel,  &
                                OutNodeFace_ProcessedDataR, OutNodeFace_FaceDataR, .true.)
                            !% --- close this file so that the unit# can be used for a new file
                            close(fU_nodeFaceFV_csv)
                        end do !% mm

                    end do !% kk
                end if !% NtotalOutputFaces > 0
            !% --- finished writing all Node output files for NodeFace

                !print *, 'OOO '
        end do !% ii
        if (verbose) write(*,*) '**** Finished writing output files ****'

        !% ----------------------------
        !% --- DEALLOCATE LOCAL STORAGE
        if (NtotalOutputElements > 0) then
            if (nOutLink > 0) then
                deallocate(OutLink_pSWMMidx, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check(deallocation_status, emsg, 'OutLink_pSWMMidx')

                deallocate(OutLink_ElemDataR, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check(deallocation_status, emsg, 'OutLink_ElemDataR')

                deallocate(OutLink_ProcessedDataR, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check(deallocation_status, emsg, 'OutLink_ProcessedDataR')

                deallocate(Outlink_pOutElemIdx, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check(deallocation_status, emsg, 'Outlink_pOutElemIdx')

                deallocate(isOutLinkWriteFVonly, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check(deallocation_status, emsg, 'isOutLinkWriteFVonly')

                deallocate(OutLink_N_elem_in_link, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check(deallocation_status, emsg, 'OutLink_N_elem_in_link')

            end if

            if (nOutNodeElem > 0) then

                deallocate(OutNodeElem_pSWMMidx, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check(deallocation_status, emsg, 'OutNodeElem_pSWMMidx')

                deallocate(OutNodeElem_ElemDataR, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check(deallocation_status, emsg, 'OutNodeElem_ElemDataR')

                deallocate(OutNodeElem_ProcessedDataR, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check(deallocation_status, emsg, 'OutNodeElem_ProcessedDataR')

                deallocate(OutNodeElem_pOutElemIdx, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check(deallocation_status, emsg, 'OutNodeElem_pOutElemIdx')

                deallocate(isOutNodeElemWriteFVonly, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check(deallocation_status, emsg, 'isOutNodeElemWriteFVonly')

                deallocate(OutNodeElem_N_elem_in_node, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check(deallocation_status, emsg, 'OutNodeElem_N_elem_in_node')

            end if
        end if

        if (NtotalOutputFaces > 0) then
            if (nOutNodeFace > 0) then

                deallocate(OutNodeFace_pSWMMidx, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check(deallocation_status, emsg, 'OutNodeFace_pSWMMidx')

                deallocate(OutNodeFace_FaceDataR, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check(deallocation_status, emsg, 'OutNodeFace_FaceDataR')

                deallocate(OutNodeFace_ProcessedDataR, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check(deallocation_status, emsg, 'OutNodeFace_ProcessedDataR')

                deallocate(OutNodeFace_pOutFaceIdx, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check(deallocation_status, emsg, 'OutNodeFace_pOutFaceIdx')

                deallocate(isOutNodeFaceWriteFVonly, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check(deallocation_status, emsg, 'isOutNodeFaceWriteFVonly')

                deallocate(OutNodeFace_N_face_in_node, stat=deallocation_status, errmsg=emsg)
                call util_deallocate_check(deallocation_status, emsg, 'OutNodeFace_N_face_in_node')

            end if
        end if

        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine outputML_convert_elements_to_linknode_and_write
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine outputML_csv_header &
        (funitIn, nType, nTotalTimeLevels, thistype, thisIndex,   &
        startdate, startEpoch,                &
        output_type_names, output_type_units, &
        elementsInLink,                       &
        tlinkname, ModelRunID,time_units_str, &
        FeatureType, isFV)
        !%-----------------------------------------------------------------------------
        !% standard header for csv output file
        !% Should be independent of global data
        !%-----------------------------------------------------------------------------
        integer, intent(in)             :: funitIn  !% formatted file unit, must be open
        integer, intent(in)             :: nType    !% number of data types (excluding time)
        integer, intent(in)             :: nTotalTimeLevels !% expected number of time rows in file
        integer, intent(in)             :: startdate(6)      !% yr, month, day, hr, min sec of model start date
        integer, intent(in)             :: thistype !% type index for FV output
        integer, intent(in)             :: thisIndex !% global index for link or node
        real(8), intent(in)             :: startEpoch    !% start date in Epoch days

        character(len=*), intent(in)    :: output_type_names(:) !% must be size nType
        character(len=*), intent(in)    :: output_type_units(:) !% must be size nType
        character(len=*), intent(in)    :: time_units_str

        integer, intent(in)             :: elementsInLink(:)

        character(len=*), intent(in)    :: tlinkname !% this linkID from SWMM
        character(len=*), intent(in)    :: ModelRunID   !% datetime stamp from run

        integer, intent(in)             :: FeatureType !% (e.g., LinkOut, NodeOut)

        logical, intent(in)             :: isFV      !% true for a FV file

        integer :: mm
        character(64) :: subroutine_name = 'outputML_csv_header'
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% -- ROW 1 --- LINK OR NODE ID (keyword, string)
        write(funitIn,fmt='(2a)') 'SWMM_ID: ,', trim(tlinkname)

        !% --- ROW 2 --- FEATURE TYPE (keyword, string)
        select case (FeatureType)
            case (LinkOut)
                write(funitIn,fmt='(2a)') 'FeatureType: ,', 'Link'
            case (NodeElemOut)
                write(funitIn,fmt='(2a)') 'FeatureType: ,', 'Node(FVelement)'
            case (NodeFaceOut)
                write(funitIn,fmt='(2a)') 'FeatureType: ,', 'Node(FVface)'
            case default
                write(*,'(A)') 'ERROR (code): Unknown FeatureType of ',FeatureType
                stop
        end select

        !% --- ROW 3 --- SWMM INDEX NUMBER IN CODE
        select case (FeatureType)
            case (LinkOut)
                write(funitIn,fmt='(a,i8)') 'CODE(...link_Gidx_SWMM): ,', thisIndex
            case (NodeElemOut)
                write(funitIn,fmt='(a,i8)') 'CODE(...node_Gidx_SWMM): ,', thisIndex
            case (NodeFaceOut)
                write(funitIn,fmt='(a,i8)') 'CODE(...node_Gidx_SWMM): ,', thisIndex
            case default
                write(*,'(A)') 'ERROR (code): Unknown FeatureType of ',FeatureType
                stop
        end select

        !% --- ROW 4 --- MODEL RUN ID (keyword, string)
        write(funitIn,fmt='(2a)') 'ModelRunID: ,',trim(ModelRunID)

        !% --- ROW 5 -- MODEL START DAY BY SWMM EPOCH DAYS (keyword, real)
        write(funitIn,fmt='(a,G0.16,a,i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')  &
            'Model_start_day(epoch):,',startEpoch

        !% --- ROW 6 -- MODEL START DAY BY yyyy_mm_dd (keyword, string)
        write(funitIn,fmt='(a,i4,a,i2.2,a,i2.2)') &
             'Model_start_day(yyyy_mm_ss):,', &
            startdate(1),'_',startdate(2),'_',startdate(3)

        !% --- ROW 7 -- MODEL START TIME BY hh:mm:ss (keyword, string)
        write(funitIn,fmt='(a,i2.2,a,i2.2,a,i2.2)')  &
            'Model_start_time(hh:mm:ss):,', &
            startdate(4),":",startdate(5),":",startdate(6)

        !% --- ROW 8 -- NUMBER OF DATA HEADER ROWS (keyword integer)
        write(funitIN,fmt='(2a,i8)') 'NumberOfDataHeaderRows:',',',3

        !% --- ROW 9 -- HEADER ROW CONTENTS (keyword, string, string ... to NumberOfHeaderRows)
        select case (FeatureType)
            case (LinkOut)
                write(funitIn,fmt='(*(a))') 'HeaderRowsContain: ',',', 'ElementID',',','DataType',',','Units'
            case (NodeElemOut)
                write(funitIn,fmt='(*(a))') 'HeaderRowsContain: ',',', 'ElementID',',','DataType',',','Units'
            case (NodeFaceOut)
                write(funitIn,fmt='(*(a))') 'HeaderRowsContain: ',',', 'FaceID',',','DataType',',','Units'
        end select

        !% --- ROW 10 --- EXPECTED NUMBER OF DATA ROWS (time levels)
        write(funitIn,fmt='(a,i8)') 'NumberDataRows: ,',nTotalTimeLevels

        !% --- ROW 11 --- EXPECTED NUMBER OF DATA COLUMNS
        write(funitIn,fmt='(a,i8)') 'NumberDataColumns: ,',nType +1 !%(including time column))

        !% --- ROW 12 --- BEGIN STATEMENT
        write(funitIn,fmt='(a)') 'BEGIN_HEADERS_AND_DATA'

        !% --- ROW 13 --- 1st HEADER ROW -- ELEMENT INDEX
        if (isFV) then
            write(funitIn,fmt='(*(i8,a))',advance='no') 0,','  !% time column
            do mm=1,nType-1
                write(funitIn,fmt='(*(i8,a))',advance='no')  elementsInLink(mm),','
            end do
            write(funitIn,fmt='(i8)') elementsInLink(nType)
        else
            select case (FeatureType)
                case (LinkOut)
                    write(funitIn,fmt='(*(i8,a))',advance='no') 0,','
                    do mm=1,nType-1
                        write(funitIn,fmt='(*(i8,a))',advance='no') 0,','
                    end do
                    write(funitIn,fmt='(i8)') 0
                case (NodeElemOut,NodeFaceOut)
                    write(funitIn,fmt='(*(i8,a))',advance='no') 0,','
                    do mm=1,nType-1
                        write(funitIn,fmt='(*(i8,a))',advance='no') elementsInLink(thisType),','
                    end do
                    write(funitIn,fmt='(i8)') elementsInLink(thisType)
                case default
                    write(*,*) 'ERROR (code): unexpected case default 389705'
                    stop
            end select
        end if

        !% --- ROW 14 -- 2nd HEADER ROW -- DATA TYPE
        if (isFV) then
            write(funitIn,fmt='(2a)',advance='no') trim(output_type_names(1)),','
            do mm=1,nType-1
                write(funitIn,fmt='(2a)',advance='no')  trim(output_type_names(thisType)),','
            end do
            write(funitIn,fmt='(a)') trim(output_type_names(thisType))
        else
            do mm=1,nType
                write(funitIn,fmt='(2a)',advance='no')  trim(output_type_names(mm)),','
            end do
            write(funitIn,fmt='(a)') trim(output_type_names(nType+1))
        end if

        !% --- ROW 15 --- 3rd HEADER ROW -- DATA UNITS
        if (isFV) then
            write(funitIn,fmt='(2a)',advance='no') output_type_units(1),','
            do mm=2,nType
                write(funitIn,fmt='(2a)',advance='no')  trim(output_type_units(thisType)),','
            end do
            write(funitIn,fmt='(a)') trim(output_type_units(thisType))
        else
            do mm=1,nType
                write(funitIn,fmt='(2a)',advance='no') trim(output_type_units(mm)),','
            end do
            write(funitIn,fmt='(a)') trim(output_type_units(nType+1))
        end if

        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine outputML_csv_header
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine outputML_unf_header &
        (funitIn, nType, nTotalTimeLevels,    &
         startdate, startEpoch,                &
         output_type_names, output_type_units, &
         elementsInLink,                       &
         tlinkname, ModelRunID, time_units_str, isFV)
        !%-----------------------------------------------------------------------------
        !% standard header for unf output file
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: funitIn  !% formatted file unit, must be open
        integer, intent(in) :: nType    !% number of data types (excluding time)
        integer, intent(in) :: nTotalTimeLevels !% expected number of time rows in file
        integer, intent(in) :: startdate(6)
        real(8), intent(in) :: startEpoch    !% start date in Epoch days

        character(len=*), intent(in)    :: output_type_names(:) !% must be size nType
        character(len=*), intent(in)    :: output_type_units(:) !% must be size nType
        character(len=*), intent(in)    :: time_units_str

        integer, intent(in)             :: elementsInLink(:)

        character(len=*), intent(in)    :: tlinkname  !% this linkID from SWMM
        character(len=*), intent(in)    :: ModelRunID   !% datetime stamp from run

        logical, intent(in)             :: isFV      !% true for a FV file

        character(64) :: subroutine_name = 'outputML_unf_header'
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% HACK --- THIS NEEDS TO BE REVISED

        !% --- write the link ID
        write(funitIn) 'LinkID:', trim(tlinkname)

        !% --- write timestamp for model run
        write(funitIn) 'ModelRunID:',ModelRunID

        !% --- write the initial epoch days and date/time at top of file
        write(funitIn)  &
            'Model_start_day(epoch):,',startEpoch

        write(funitIn) 'Model_start_date(y,m,d,h,m,s)', startdate(1:6)

        !% --- write the time units (column1)
        write(funitIn) 'TimeUnits(Column01)in',time_units_str

        !% --- write the expected number of time levels
        write(funitIn) 'ExpectedTimeLevels(rows):',nTotalTimeLevels

        !% ---- write the number of data types
        write(funitIn) 'NumberDataTypes(columns):',nType

        if (isFV) then
            write(funitIn) 'FirstRowAreElementIndexes(0=TimeColumn)'
            write(funitIn) 0, elementsInLink(:)
        else
            !% --- write the headers for the 2D array of data
            write(funitIn) 'FirstRowAreDataNames_SecondRowAreDataUnits'
            write(funitIn) output_type_names(:)
        end if

        if (.not. isFV) then
            !% --- write the units used
            write(funitIn) output_type_units(:)
        end if

        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine outputML_unf_header
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine outputML_csv_writedata &
        (funitIn, idx1, nIdx2, idx3, nLevel, Out_ProcessedDataR, Out_ElemDataR, isFV)
        !%-----------------------------------------------------------------------------
        !% writes data set of rows (time levels) and columns (data type)
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: funitIn !% file unit number to write to
        integer, intent(in) :: nLevel  !% number of time levels in this data set
        integer, intent(in) :: idx1    !% the link being output (kk)
        integer, intent(in) :: nIdx2   !% N items in columns (nType or SWMMlink_num_elements)
        integer, intent(in) :: idx3    !% the single data type being processed (FV only)
        real(8), intent(in) :: Out_ProcessedDataR(:,:,:) !% (link/node,type,timelevel)
        real(8), intent(in) :: Out_ElemDataR(:,:,:,:)    !% (link/node,element,type,timelevel)
        logical, intent(in) :: isFV   !% is finite volume output

        character(64) :: subroutine_name = 'outputML_csv_writedata'

        integer :: mm
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        if (isFV) then
            !% --- FV write is columns that are elements of the idx1 link
            do mm=1,nLevel
                write(funitIn,'(*(G0.6 : ","))') &
                 Out_ElemDataR(idx1,1,1,mm),   &
                 Out_ElemDataR(idx1,1:nIdx2,idx3,mm)
            end do
        else
            !% --- nonFV write is columns of data types (e.g., Area, Velocity)
            do mm=1,nLevel
                write(funitIn,'(*(G0.6 : ","))') Out_ProcessedDataR(idx1,1:nIdx2,mm)
            end do
        end if

        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine  outputML_csv_writedata
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine outputML_unf_writedata &
        (funitIn, nType, nLevel, klink, OutLink_ProcessedDataR)
        !%-----------------------------------------------------------------------------
        !% writes data set of rows (time levels) and columns (data type)
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: funitIn !% file unit number to write to
        integer, intent(in) :: nType   !% number of data types (excluding time)
        integer, intent(in) :: nLevel  !% number of time levels in this data set
        integer, intent(in) :: klink   !% the link being output
        real(8), intent(in) :: OutLink_ProcessedDataR(:,:,:)

        integer :: mm

        character(64) :: subroutine_name = 'outputML_unf_writedata'
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        do mm = 1,nLevel
            write(funitIn) OutLink_ProcessedDataR(klink,1:nType,mm)
        end do

        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine outputML_unf_writedata
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine outputD_update_swmm_out()
        ! !%-----------------------------------------------------------------------------
        ! !% Description:
        ! !%
        ! !%-----------------------------------------------------------------------------
        ! character(len = 250) :: fname
        ! character(len = 24) :: timestamp
        ! logical :: wrote_all_links = .false.
        ! logical :: wrote_all_nodes = .false.
        ! integer :: ii, rc, node_idx, link_idx
        ! integer, allocatable :: fus_nodes(:), fus_links(:)
        ! real(8) :: node_head, node_result
        ! real(8) :: link_flowrate, link_result
        ! real(8) :: timesecs
        ! character(64) :: subroutine_name = 'outputD_update_swmm_out'
        ! !%--------------------------------------------------------------------------
        ! if (setting%Debug%File%output) &
        !     write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        ! if (this_image() == 1) then
        !     allocate(fus_nodes(size(node%P%have_output)))
        !     allocate(fus_links(size(link%P%have_output)))
        !     fus_links = nullvalueI
        !     fus_nodes = nullvalueI

        !     do while(.not. (wrote_all_links .and. wrote_all_nodes))

        !         if (size(node%P%have_output) < 1) then
        !             wrote_all_nodes = .true.
        !         else
        !             do ii=1, size(node%P%have_output)
        !                 !print *, 'node ',ii,  size(node%P%have_output),  node%P%have_output(ii)
        !                 node_idx = node%P%have_output(ii)
        !                 if ((node_idx > 0) .and. (node_idx < nullvalueI)) then
        !                     if (fus_nodes(ii) == nullvalueI) then
        !                         !% open files to process .out
        !                         !fname = "swmm5_output/node/"//trim(node%names(node_idx)%str)//".csv"
        !                         fname = trim(setting%File%swmm5_output_node_folder) &
        !                             //trim(node%names(node_idx)%str)//".csv"
        !                         open(action='read', file=trim(fname), iostat=rc, newunit=fus_nodes(ii))
        !                         read(fus_nodes(ii), *, iostat = rc) timestamp
        !                     end if

        !                     read(fus_nodes(ii), "(A,2F10.8)", iostat = rc) timestamp, timesecs, node_head
        !                     if (rc /= 0) then
        !                         wrote_all_nodes = .true.
        !                         close(fus_nodes(ii))
        !                         !% Write line of .out
        !                         !exit
        !                     else
        !                         wrote_all_nodes = .false.
        !                         node_result = node_head - node%R(node_idx,nr_Zbottom)
        !                         !% stage values in .out
        !                         call interface_update_nodeResult(node_idx, api_output_node_depth, node_result)
        !                     end if
        !                 else
        !                     if (ii .eq. size(node%P%have_output)) then
        !                         wrote_all_nodes = .true.
        !                     end if
        !                 end if
        !                 if (ii .lt. size(node%P%have_output)) then
        !                     wrote_all_nodes = .false.
        !                 end if
        !             end do
        !         end if

        !         if (size(link%P%have_output) < 1) then
        !             wrote_all_links = .true.
        !         else
        !             do ii=1, size(link%P%have_output)
        !                 !print *, 'link ',ii,  size(link%P%have_output),  link%P%have_output(ii)
        !                 link_idx = link%P%have_output(ii)
        !                 if ((link_idx > 0) .and. (link_idx < nullvalueI)) then
        !                     if (fus_links(ii) == nullvalueI) then
        !                         !% open files to process .out
        !                         !fname = "swmm5_output/link/"//trim(link%names(link_idx)%str)//".csv"
        !                         fname =  trim(setting%File%swmm5_output_link_folder) &
        !                             //trim(link%names(link_idx)%str)//".csv"
        !                         open(action='read', file=trim(fname), iostat=rc, newunit=fus_links(ii))
        !                         read(fus_links(ii), *, iostat = rc) timestamp
        !                     end if

        !                     read(fus_links(ii), "(A,2F10.8)", iostat = rc) timestamp, timesecs, link_flowrate
        !                     if (rc /= 0) then
        !                         wrote_all_links = .true.
        !                         close(fus_links(ii))
        !                         !% Write line of .out
        !                         !exit
        !                     else
        !                         wrote_all_links = .false.
        !                         link_result = link_flowrate
        !                         !% stage values in .out
        !                         call interface_update_linkResult(link_idx, api_output_link_flow, link_result)
        !                     end if
        !                 else
        !                     if (ii .eq. size(link%P%have_output)) then
        !                         wrote_all_links = .true.
        !                     end if
        !                 end if
        !                 if (ii .lt. size(link%P%have_output)) then
        !                     wrote_all_links = .false.
        !                 end if
        !             end do
        !         end if
        !         call interface_write_output_line(timesecs)
        !     end do

        !     deallocate(fus_links)
        !     deallocate(fus_nodes)

        !     ! do ii = 1, size(link%P%have_output)
        !     !     link_idx = link%P%have_output(ii)
        !     !     call interface_export_link_results(link_idx)
        !     ! end do
        ! end if

        ! if (setting%Debug%File%output) &
        ! write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"


    !         !% --- check to see if file exists
    !         inquire (FILE=setting%File%links_input_file, EXIST=file_exists)

    !         if (file_exists) then
    !             !% open the csv file of link names
    !             open(unit= setting%File%UnitNumber%links_input_file, &
    !                  file=trim(setting%File%links_input_file), &
    !                  action='read', &
    !                  iostat=rc)
    !             if (rc /= 0) then
    !                 write (*, '(3a, i0)') 'ERROR (user): Opening file ', trim(setting%File%links_input_file), ' failed: ', rc
    !                 stop
    !             end if
    !             pp = 1 ! parent link
    !             endoffile = .false.
    !             !% ---loop through till the end of the file and save the valid links
    !             do while ((.not. endoffile) .and. (pp .le. setting%Output%max_links_csv))
    !                 inquire (UNIT=setting%File%UnitNumber%links_input_file, position=thispos)
    !                 if (thispos .eq. 'APPEND') then
    !                     endoffile = .true.
    !                     pp = pp-1 !% so that pp=0 indicates nothing read (empty file)
    !                     exit !end the do loop
    !                 end if
    !                 !% --- read in the link name from the csv
    !                 read(setting%File%UnitNumber%links_input_file, *,  iostat = rc) link_name
    !                 !% --- crash on error
    !                 if (rc /= 0) then
    !                     close(setting%File%UnitNumber%links_input_file)
    !                     !exit
    !                     write(*,"(A)") 'ERROR (user): reading file ', trim(setting%File%links_input_file)
    !                     write(*,"(A,i5)") 'failed before end of file with error ',rc
    !                     stop
    !                 end if

    !                 !% --- converting link name to link idx using the interface
    !                 link_idx = interface_find_object(object_type=API_LINK, object_name = link_name)
    !                 !% --- crash on error
    !                 if (link_idx == 0) then
    !                     write(*, "(A)") "ERROR (user): Link " // trim(link_name) // " in " // &
    !                         trim(setting%File%links_input_file) // " couldn't be found"
    !                     !exit
    !                     stop
    !                 end if

    !                 !% --- store index of link for output and increase index
    !                 link_output_idx(pp) = link_idx
    !                 pp = pp + 1

    !                 !% checking if the link is split across processors if
    !                 !% so then store the id of the phantom link for output
    !                 phantom_counter = 0
    !                 if (link%I(link_idx, li_parent_link) == link_idx) then
    !                     do jj = SWMM_N_link+1, N_link
    !                         if (link%I(jj, li_parent_link) == link_idx) then
    !                             link_output_idx(pp+phantom_counter+1) = jj
    !                             phantom_counter = phantom_counter + 1
    !                         end if
    !                     end do
    !                 end if
    !                 pp = pp + phantom_counter + 1
    !                 link%I(link_idx,li_num_phantom_links) = phantom_counter

    !                 !% --- check for stopping due to max links setting.
    !                 if (pp .ge. setting%Output%max_links_csv) then
    !                     if (setting%Output%Warning) then
    !                         write(*,"(A)") 'WARNING: stopped reading links_input_file file due to excessive number of links'
    !                         write(*,"(A)") 'Filename = ',trim(setting%File%links_input_file)
    !                         write(*,"(A,i5)") 'Maximum links set by setting.Output.max_links_csv as: ',setting%Output%max_links_csv
    !                     end if
    !                 end if
    !             end do
    !             !% --- if exited without reading (empty file) send warning
    !             if (pp == 0) then
    !                 if (setting%Output%Warning) then
    !                     write(*,"(A)") 'WARNING: did not find an links in the links_input_file'
    !                     write(*,"(A)") 'Filename = ',trim(setting%File%links_input_file)
    !                 end if
    !             end if
    !         end if

    !         if (.not. file_exists) then
    !             !% --- if links_input_file is not specified we output all the links up to the maximum allowed
    !             pp = 1 !% parent link
    !             !do link_idx = 1, SWMM_N_link
    !             do while ( (pp <= SWMM_N_link) .and. (pp .le. setting%Output%max_links_csv))
    !                 phantom_counter = 0
    !                 link_output_idx(pp) = link_idx
    !                 !% --- only parent links have associated phantoms
    !                 if (link%I(link_idx, li_parent_link) == link_idx) then
    !                     do jj = SWMM_N_link+1, N_link
    !                         if (link%I(jj, li_parent_link) == link_idx) then
    !                             link_output_idx(pp+phantom_counter+1) = jj
    !                             phantom_counter = phantom_counter + 1
    !                         end if
    !                     end do
    !                 end if
    !                 pp = pp + phantom_counter + 1
    !                 link%I(link_idx,li_num_phantom_links) = phantom_counter

    !                 !% --- check for stopping due to max links setting.
    !                 if (pp .ge. setting%Output%max_links_csv) then
    !                     if (setting%Output%Warning) then
    !                         write(*,"(A)") 'WARNING: stopped selecting links for csv output due to excessive number of links'
    !                         write(*,"(A,i5)") 'Maximum links set by setting.Output.max_links_csv as: ',setting%Output%max_links_csv
    !                     end if
    !                 end if
    !             end do
    !         end if

    !         if (setting%Debug%File%output) &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     end subroutine outputD_read_csv_link_names
    ! !%
    ! !%==========================================================================
    ! !%==========================================================================
    ! !%
    !     subroutine outputD_read_csv_node_names()
    !         !%-----------------------------------------------------------------------------
    !         !% Description:
    !         !% Reading the node input file and store nodes in note_output_idx.
    !         !% If file does not exist, then choose the first 1:setting%Output%max_nodes_csv
    !         !% nodes for output.
    !         !%-----------------------------------------------------------------------------
    !         character(len = 250) :: node_name
    !         integer :: rc, fu, ii, node_idx, maxnode
    !         character(16) :: thispos
    !         logical :: file_exists, endoffile
    !         character(64) :: subroutine_name = 'outputD_read_csv_node_names'
    !         !%--------------------------------------------------------------------------
    !         if (setting%Debug%File%output) &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !         !% --- abandon procedure if printout of nodes not needed
    !         node_output_idx = nullvalueI
    !         if (.not. setting%Output%print_nodes_csv) return

    !         !% --- check to see if file exists
    !         inquire (FILE=setting%File%nodes_input_file, EXIST=file_exists)

    !         !% --- read the file that exists
    !         if (file_exists) then
    !             !% --- open the csv file of link names
    !             open(unit= setting%File%UnitNumber%nodes_input_file, &
    !                  file=trim(setting%File%nodes_input_file), &
    !                  action='read', &
    !                  iostat=rc)
    !             if (rc /= 0) then
    !                 write (*, '(3a, i0)') 'ERROR (user): Opening file ', trim(setting%File%nodes_input_file), ' failed: ', rc
    !                 stop
    !             end if
    !             ii = 1
    !             endoffile = .false.
    !             do while ((.not. endoffile) .and. (ii .le. setting%Output%max_nodes_csv))
    !                 inquire (UNIT=setting%File%UnitNumber%nodes_input_file, position=thispos)
    !                 if (thispos .eq. 'APPEND') then
    !                     endoffile = .true.
    !                     ii = ii-1 !% so that ii=0 indicates nothing read (empty file)
    !                     exit !end the do loop
    !                 end if
    !                 !% --- read in nodes
    !                 read(setting%File%UnitNumber%nodes_input_file, *, iostat = rc) node_name
    !                 if (rc /= 0) then
    !                     node_output_idx(ii:) = nullvalueI
    !                     close(setting%File%UnitNumber%nodes_input_file)
    !                     !exit
    !                     write(*,"(A)") 'ERROR (user): reading file ', trim(setting%File%nodes_input_file)
    !                     write(*,"(A,i5)") 'failed before end of file with error ',rc
    !                     stop
    !                 end if
    !                 !% --- converting node name to node idx using the interface
    !                 node_idx = interface_find_object(object_type=API_NODE, object_name = node_name)
    !                 if (node_idx == 0) then
    !                     write(*, "(A)") "Node " // trim(node_name) // " in " // &
    !                     trim(setting%File%nodes_input_file) // " couldn't be found"
    !                     stop
    !                 end if
    !                 node_output_idx(ii) = node_idx
    !                 ii = ii + 1
    !             end do
    !             !% --- if exited without reading (empty file) send warning
    !             if (ii == 0) then
    !                 if (setting%Output%Warning) then
    !                     write(*,"(A)") 'WARNING: did not find an nodes in the nodes_input_file'
    !                     write(*,"(A)") 'Filename = ',trim(setting%File%nodes_input_file)
    !                 end if
    !             end if
    !         end if

    !         !% --- store 1:setting%Output%max_nodes_csv nodes for output when file doesn't exist
    !         if (.not. file_exists) then
    !             !% --- Output all nodes (up to max) if CSV file
    !             maxnode = max( SWMM_N_node, setting%Output%max_nodes_csv)
    !             node_output_idx = (/ (ii, ii =1, maxnode)/)
    !             !% --- send warning if output is truncated
    !             if (SWMM_N_node > maxnode) then
    !                 write(*,"(A)") 'WARNING: stopped selecting nodes for csv output due to excessive number of nodes'
    !                 write(*,"(A,i5)") 'Maximum nodes set by setting.Output.max_nodes_csv as: ',setting%Output%max_links_csv
    !             end if
    !         end if

    !         if (setting%Debug%File%output) &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
!     end subroutine outputD_read_csv_node_names
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine outputD_create_link_files
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !%   Creation of link files and header for the files
    !     !%-----------------------------------------------------------------------------
    !     integer :: ii,fu, open_status, temp_link_idx
    !     character(len = 250) :: file_name
    !     character(len = 100) :: link_name
    !     character(len = 5)   :: str_image
    !     character(len = 10)  :: str_idx
    !     character(64) :: subroutine_name = 'outputD_create_link_files'
    !     !%--------------------------------------------------------------------------
    !     if (setting%Debug%File%output) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !     write(str_image, '(i5.5)') this_image()

    !     do ii=1, size(link%P%have_output)

    !         !% HACK needed to ensure that a valid index is provided
    !         if (link%P%have_output(ii) > 0) then

    !             !% check if the link is a phantom link and if so find original link name and open correct file for the correct processor
    !             !% otherwise open file in the usual format of "link name_imageID.csv"
    !             if (link%P%have_output(ii) > size(link%names(:))) then
    !                 write(str_idx, '(i5.5)') link%P%have_output(ii)
    !                 temp_link_idx = link%P%have_output(ii)
    !                 !file_name = "debug_output/link/"//trim(link%names(link%I(temp_link_idx,li_parent_link))%str) &
    !                 !    //"_"//trim(str_image)//"_"//trim(str_idx)//".csv"
    !                 file_name = trim(setting%File%debug_output_link_folder) &
    !                     //trim(link%names(link%I(temp_link_idx,li_parent_link))%str) &
    !                     //"_"//trim(str_image)//"_"//trim(str_idx)//".csv"
    !             else
    !                 !file_name = "debug_output/link/"//trim(link%names(link%P%have_output(ii))%str) &
    !                 !    //"_"//trim(str_image)//".csv"
    !                 file_name = trim(setting%File%debug_output_link_folder) &
    !                     //trim(link%names(link%P%have_output(ii))%str) &
    !                     //"_"//trim(str_image)//".csv"
    !             end if

    !             open(newunit=fu, file = file_name, status = 'replace',access = 'sequential', &
    !                 form   = 'formatted', action = 'write', iostat = open_status)

    !             if (open_status /= 0) then
    !                 write (*, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
    !                 stop
    !             end if

    !             !% Write the header of the file, set end for next write and then close file
    !             write(fu, *) "Timestamp,Time_In_Secs,flowrate"
    !             endfile(fu)
    !             close(fu)
    !         end if
    !     end do

    !     if (setting%Debug%File%output) &
    !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine outputD_create_link_files
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine outputD_create_node_files
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !%   Creation of node files and header for the files
    !     !%-----------------------------------------------------------------------------
    !     integer :: ii,fu, open_status
    !     character(len = 250) :: file_name
    !     character(len = 100) :: node_name
    !     character(len = 5)   :: str_image
    !     character(64) :: subroutine_name = 'outputD_create_node_files'
    !     !%--------------------------------------------------------------------------
    !     if (setting%Debug%File%output) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !     !% Get current image as a string
    !     write(str_image, '(i5.5)') this_image()

    !     do ii=1, size(node%P%have_output)

    !         !% HACK -- needed because node%P%have_output(ii)= 0 is assigned somehwere
    !         if ((node%P%have_output(ii) > 0) .and. (node%P%have_output(ii) < nullvalueI )) then

    !             !% Open the node file
    !             !file_name = "debug_output/node/"//trim(node%names(node%P%have_output(ii))%str) &
    !             !    //"_"//trim(str_image)//".csv"
    !             file_name = trim(setting%File%debug_output_node_folder) &
    !                 //trim(node%names(node%P%have_output(ii))%str) &
    !                 //"_"//trim(str_image)//".csv"

    !             print *, 'a :',ii, node%P%have_output(ii)
    !             print *, 'b :',size(node%Names)
    !             print *, 'c : ',node%Names(node%P%have_output(ii))%str

    !             print *, 'd :', file_name

    !             open(newunit=fu, file = file_name, status = 'replace',access = 'sequential', &
    !                 form   = 'formatted', action = 'write', iostat = open_status)

    !             if (open_status /= 0) then
    !                 write (*, '(A)') 'Opening file failed'
    !                 write (*, '(A)') trim(FILE_NAME)
    !                 write (*, '(i0)') open_status
    !                 stop
    !             end if

    !             !% Write the header, this endfile and close the file
    !             write(fu, *) "Timestamp,Time_In_Secs,Head"
    !             endfile(fu)
    !             close(fu)
    !         end if
    !     end do

    !     if (setting%Debug%File%output) &
    !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine outputD_create_node_files
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine outputD_write_link_files
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Called at the report time step to calculate the flowrate in the link and
    !     !% write it to the file
    !     !%-----------------------------------------------------------------------------
    !     integer :: ii, fu, open_status, yr, mnth, dy, hr, min, sec
    !     integer :: start_elem, end_elem, temp_link_idx
    !     real(8) :: time_secs, time_epoch, avg_flowrate
    !     character(len = 250) :: file_name
    !     character(len = 100) :: link_name
    !     character(len = 5)   :: str_image
    !     character(len = 10)  :: str_idx
    !     character(64) :: subroutine_name = 'outputD_write_link_files'
    !     !%--------------------------------------------------------------------------
    !     if (setting%Debug%File%output) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !     write(str_image, '(i5.5)') this_image()
    !     time_secs = setting%Time%Now
    !     time_epoch = util_datetime_secs_to_epoch(time_secs)
    !     call util_datetime_decodedate(time_epoch, yr, mnth, dy)
    !     call util_datetime_decodetime(time_epoch, hr, min, sec)

    !     do ii=1, size(link%P%have_output)

    !         if (link%P%have_output(ii) > 0 ) then

    !             !% store the store the location of the start and end elem for easier reading
    !             start_elem = link%I(link%P%have_output(ii),li_first_elem_idx)
    !             end_elem = link%I(link%P%have_output(ii),li_last_elem_idx)

    !             !% calculate average flowrate by summing up the elems and diving about the number of elems
    !             avg_flowrate = sum(elemR(start_elem:end_elem,er_Flowrate))/(end_elem-start_elem)

    !             if (start_elem == end_elem) then
    !                 avg_flowrate = elemR(start_elem,er_flowrate)
    !             end if
    !             !% check if the link is a phantom link and if so find original link name and open correct file for the correct processor
    !             !% otherwise open file in the usual format of "link name_imageID.csv"
    !             if (link%P%have_output(ii) > size(link%names(:))) then
    !                 write(str_idx, '(i5.5)') link%P%have_output(ii)
    !                 temp_link_idx = link%P%have_output(ii)

    !                 !file_name = "debug_output/link/"//trim(link%names(link%I(temp_link_idx,li_parent_link))%str) &
    !                 !    //"_"//trim(str_image)//"_"//trim(str_idx)//".csv"
    !                 file_name = trim(setting%File%debug_output_link_folder) &
    !                     //trim(link%names(link%I(temp_link_idx,li_parent_link))%str) &
    !                     //"_"//trim(str_image)//"_"//trim(str_idx)//".csv"
    !             else
    !                 !file_name = "debug_output/link/"//trim(link%names(link%P%have_output(ii))%str) &
    !                 !    //"_"//trim(str_image)//".csv"
    !                 file_name = trim(setting%File%debug_output_link_folder) &
    !                     //trim(link%names(link%P%have_output(ii))%str) &
    !                     //"_"//trim(str_image)//".csv"
    !             end if

    !             open(newunit=fu, file = file_name, status = 'old',access = 'append', &
    !                 form   = 'formatted', action = 'write', iostat = open_status)

    !             !% writing timestamped output to file for average flowrate across the link
    !             write(fu,fmt='(i4, 2(a,i2.2))',advance = 'no') yr,"_",mnth,"_",dy
    !             write(fu,fmt = '(A)',advance = 'no') '_'
    !             write(fu,fmt='(2(i2.2,a), i2.2)',advance = 'no') hr,":",min,":",sec
    !             write(fu,'(A)', advance = 'no') ','
    !             write(fu, '(F0.16)', advance = 'no') time_secs
    !             write(fu,'(A)', advance = 'no') ','
    !             write(fu, '(*(G0.6 : ","))') avg_flowrate

    !             !% set the end of the file for next write and close file
    !             endfile(fu)
    !             close(fu)

    !         end if
    !     end do

    !     if (setting%Debug%File%output) &
    !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine outputD_write_link_files
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine outputD_write_node_files
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Called at the report time step to calculate the head in the node and
    !     !% write it to the file
    !     !%-----------------------------------------------------------------------------
    !     integer :: ii, fu, open_status, yr, mnth, dy, hr, min, sec
    !     integer :: temp_node_idx
    !     real(8) :: time_secs, time_epoch, avg_head
    !     character(len = 250) :: file_name
    !     character(len = 100) :: link_name
    !     character(len = 5)   :: str_image
    !     character(64) :: subroutine_name = 'outputD_write_node_files'

    !     !%--------------------------------------------------------------------------
    !     if (setting%Debug%File%output) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !     !% converter image ID to string, as well as get current time
    !     write(str_image, '(i5.5)') this_image()
    !     time_secs = setting%Time%Now
    !     time_epoch = util_datetime_secs_to_epoch(time_secs)
    !     call util_datetime_decodedate(time_epoch, yr, mnth, dy)
    !     call util_datetime_decodetime(time_epoch, hr, min, sec)

    !     !% loop through nodes we have to output
    !     do ii=1, size(node%P%have_output)

    !         if ((node%P%have_output(ii) > 0) .and. (node%P%have_output(ii) < nullvalueI )) then

    !             !% open node file
    !             !file_name = "debug_output/node/"//trim(node%names(node%P%have_output(ii))%str) &
    !             !    //"_"//trim(str_image)//".csv"
    !             file_name = trim(setting%File%debug_output_node_folder) &
    !                 //trim(node%names(node%P%have_output(ii))%str) &
    !                 //"_"//trim(str_image)//".csv"

    !             open(newunit=fu, file = file_name, status = 'old',access = 'append', &
    !                 form   = 'formatted', action = 'write', iostat = open_status)

    !             if (open_status /= 0) then
    !                 write (*, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
    !                 stop
    !             end if

    !             !% temp value for easier to read code
    !             temp_node_idx = node%P%have_output(ii)

    !             !% check if the node is a BC type node or a junction type node, then write appropriate data
    !             if (node%I(temp_node_idx,ni_node_type) == nBCup .or. node%I(temp_node_idx,ni_node_type) == nBCdn) then
    !                 write(fu,fmt='(i4, 2(a,i2.2))',advance = 'no') yr,"/",mnth,"/",dy
    !                 write(fu,fmt = '(A)',advance = 'no') ' '
    !                 write(fu,fmt='(2(i2.2,a), i2.2)',advance = 'no') hr,":",min,":",sec
    !                 write(fu,'(A)', advance = 'no') ','
    !                 write(fu, '(F0.16)', advance = 'no') time_secs
    !                 write(fu,'(A)', advance = 'no') ','
    !                 write(fu, '(*(G0.6 : ","))') faceR(node%I(temp_node_idx,ni_elemface_idx),fr_Head_d)

    !             else if (node%I(temp_node_idx,ni_node_type) == nJ2 .or. node%I(temp_node_idx,ni_node_type) == nJm) then
    !                 write(fu,fmt='(i4, 2(a,i2.2))',advance = 'no') yr,"/",mnth,"/",dy
    !                 write(fu,fmt = '(A)',advance = 'no') ' '
    !                 write(fu,fmt='(2(i2.2,a), i2.2)',advance = 'no') hr,":",min,":",sec
    !                 write(fu,'(A)', advance = 'no') ','
    !                 write(fu, '(F0.16)', advance = 'no') time_secs
    !                 write(fu,'(A)', advance = 'no') ','
    !                 write(fu, '(*(G0.6 : ","))') elemR(node%I(temp_node_idx,ni_elemface_idx),er_Head)
    !             !% if the node type is neither a BC or junction type then print the warning
    !             else
    !                 print *, "WARNING: node selected is neither BCup, BCdn, nJ2 or nJM node, no output will be written"
    !                 print *, "temp_node_idx",temp_node_idx
    !                 print *, "node%I(temp_node_idx,ni_node_type)", node%I(temp_node_idx,ni_node_type)
    !             end if

    !             !% call endfile and close the file
    !             endfile(fu)
    !             close(fu)
    !         end if
    !     end do

    !     if (setting%Debug%File%output) &
    !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    ! end subroutine outputD_write_node_files
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine outputD_combine_links
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Combines phantom links into single output links
    !     !%-----------------------------------------------------------------------------
    !     integer :: ii, jj, pp, rc, open_status, N_parents, N_phantoms
    !     integer :: temp_link_idx, temp_phantom_link ,link_output_idx_length
    !     integer :: start_elem, end_elem,num_elems
    !     integer, allocatable :: file_idx(:)
    !     real(8)              :: time_secs, time_epoch, flowrate
    !     real(8), pointer     :: avg_flowrate(:)
    !     real(8), allocatable :: full_length(:), tt(:)
    !     logical, allocatable :: first_iteration(:)
    !     logical, allocatable :: its_over(:)
    !     character(len = 250) :: parent_file_name, phantom_file_name
    !     character(len = 250) :: final_file_name
    !     character(len = 100) :: link_name
    !     character(len = 5)   :: str_image
    !     character(len = 10)  :: str_idx
    !     character(len = 19)  :: str_time
    !     character(64) :: subroutine_name = 'outputD_write_link_files'
    !     !%-----------------------------------------------------------------------------

    !     link_output_idx_length = count(link_output_idx(:) /= nullvalueI)
    !     N_phantoms = sum(link%I(:, li_num_phantom_links))
    !     N_parents = link_output_idx_length - N_phantoms

    !     !print *, link_output_idx_length
    !     !print *, N_phantoms
    !     !print *, N_parents


    !     allocate(first_iteration(N_parents))
    !     allocate(full_length(N_parents))
    !     allocate(its_over(N_parents))
    !     allocate(tt(N_parents))
    !     allocate(file_idx(link_output_idx_length+N_parents))

    !     avg_flowrate => link%R(:, lr_flowrate)

    !     first_iteration(:) = .true.
    !     full_length(:) = zeroR
    !     avg_flowrate(:) = zeroR
    !     its_over(:) = .false.
    !     tt(:) = 1
    !     file_idx = 0

    !     !% only execute if there are actually links to output
    !     if (sum(link_output_idx(:)) > 0) then

    !         do while (any(.not. its_over))
    !             ii = 1
    !             pp = 1
    !             do while (pp <= N_parents)
    !                 !% Set some starting values
    !                 !% we use temp_link_idx as the index of the parent link we are getting the output for
    !                 temp_link_idx = link_output_idx(ii)

    !                 if (first_iteration(pp)) then
    !                     !% define parent filename
    !                     write(str_image, '(i5.5)') link%I(temp_link_idx,li_P_image)

    !                     !parent_file_name = "debug_output/link/"// &
    !                     !    trim(link%names(link%I(temp_link_idx,li_parent_link))%str) &
    !                     !        //"_"//trim(str_image)//".csv"
    !                     parent_file_name = trim(setting%File%debug_output_link_folder)// &
    !                             trim(link%names(link%I(temp_link_idx,li_parent_link))%str) &
    !                                 //"_"//trim(str_image)//".csv"

    !                     !final_file_name = "swmm5_output/link/"//trim(link%names(temp_link_idx)%str)//".csv"
    !                     final_file_name = trim(setting%File%swmm5_output_link_folder) &
    !                         //trim(link%names(temp_link_idx)%str)//".csv"

    !                     !% Open parent file
    !                     open(newunit=file_idx(ii), action='read', file=parent_file_name, iostat=rc)
    !                     if (rc /= 0) then
    !                         write (*, '(3a, i0)') 'Opening file "', trim(parent_file_name), '" failed: ', rc
    !                         stop
    !                     end if
    !                     read (file_idx(ii), *, iostat=rc) str_time ! advance one line (skip header)

    !                     ! start calculating the full_length of the link for the phantom links
    !                     full_length(pp) = full_length(pp) + link%R(temp_link_idx,lr_AdjustedLength)

    !                     !Now we open the Final file for the link output which is in a different location and just the name of the link
    !                     !We also write the header
    !                     open(newunit=file_idx(link_output_idx_length+pp), &
    !                         file = final_file_name, status = 'replace',access = 'sequential', &
    !                         form = 'formatted', action = 'write', iostat = open_status)
    !                     if (open_status /= 0) then
    !                         write (*, '(3a, i0)') 'Opening file "', trim(Final_File_NAME), '" failed: ', open_status
    !                         stop
    !                     end if
    !                     write(file_idx(link_output_idx_length+pp), *) "Timestamp,Time_In_Secs,flowrate"

    !                     !Now we check if the parent link has phantoms related to it and loop through those values
    !                     !While looping we open each of the phantom link files and calculate the full_length
    !                     if (link%I(temp_link_idx, li_num_phantom_links) > 0) then
    !                         do jj = 1, link%I(temp_link_idx, li_num_phantom_links)
    !                             temp_phantom_link = link_output_idx(ii+jj)
    !                             write(str_image, '(i5.5)') link%I(temp_phantom_link,li_P_image)
    !                             write(str_idx, '(i5.5)')   temp_phantom_link

    !                             !phantom_file_name = "debug_output/link/"// &
    !                             !    trim(link%names(link%I(temp_phantom_link,li_parent_link))%str) &
    !                             !   //"_"//trim(str_image)//"_"//trim(str_idx)//".csv"
    !                             phantom_file_name =  trim(setting%File%debug_output_link_folder)// &
    !                                 trim(link%names(link%I(temp_phantom_link,li_parent_link))%str) &
    !                                 //"_"//trim(str_image)//"_"//trim(str_idx)//".csv"
    !                             full_length(pp) = full_length(pp) + link%R(temp_phantom_link,lr_AdjustedLength)

    !                             open(newunit=file_idx(ii+jj), action='read', file=phantom_file_name, iostat=rc)
    !                             if (rc /= 0) then
    !                                 write (*, '(3a, i0)') 'Opening file "', trim(phantom_file_name), '" failed: ', rc
    !                                 stop
    !                             end if
    !                             read (file_idx(ii+jj), *, iostat=rc) str_time ! advance one line (skip header)
    !                         end do
    !                         tt(pp) = tt(pp) + 1
    !                     else
    !                         !If the link doesn't have any phantom links related to it we close the parent and final file and simply rename and move the parent file to the final file location.
    !                         !Then increment one.
    !                         close(file_idx(ii))
    !                         close(file_idx(link_output_idx_length+pp))
    !                         call rename(trim(parent_file_name), trim(final_file_name))
    !                         its_over(pp) = .true.
    !                     end if
    !                     first_iteration(pp) = .false.
    !                 else
    !                     if (.not. its_over(pp)) then ! here we will only encounter parents with phantom links
    !                         ! Now that we have the parent file, the final file and the phantom files open
    !                         ! we can start recombing the flowrates and writing the final file
    !                         read (file_idx(ii), *, iostat=rc) str_time, time_secs, flowrate
    !                         if (rc /= 0) then
    !                             its_over(pp) = .true.
    !                             ii = ii + link%I(temp_link_idx, li_num_phantom_links) + 1
    !                             pp = pp + 1
    !                             cycle
    !                         end if

    !                         avg_flowrate(pp) = avg_flowrate(pp) + &
    !                             (flowrate * link%R(temp_link_idx,lr_AdjustedLength) / full_length(pp))

    !                         do jj = 1, link%I(temp_link_idx, li_num_phantom_links)
    !                             temp_phantom_link = link_output_idx(ii+jj)
    !                             read (file_idx(ii+jj), *, iostat=rc) str_time, time_secs, flowrate
    !                             avg_flowrate(pp) = avg_flowrate(pp) + &
    !                                 (flowrate * link%R(temp_phantom_link,lr_AdjustedLength) / full_length(pp))
    !                         end do

    !                         write(file_idx(link_output_idx_length+pp), '(A)',     advance = 'no') str_time
    !                         write(file_idx(link_output_idx_length+pp), '(A)',     advance = 'no') ','
    !                         write(file_idx(link_output_idx_length+pp), '(F0.16)', advance = 'no') time_secs
    !                         write(file_idx(link_output_idx_length+pp), '(A)',     advance = 'no') ','
    !                         write(file_idx(link_output_idx_length+pp), '(*(G0.6 : ","))') avg_flowrate(pp)

    !                         avg_flowrate(pp) = 0
    !                         tt(pp) = tt(pp) + 1
    !                     end if
    !                 end if
    !                 !% Now we increment based off of how many phantom links there where related to the parent link that we combined the output
    !                 ii = ii + link%I(temp_link_idx, li_num_phantom_links) + 1
    !                 pp = pp + 1
    !             end do
    !         end do
    !         do ii = 1, size(file_idx)
    !             close(file_idx(ii))
    !         end do
    !     end if

    !     deallocate(first_iteration)
    !     deallocate(full_length)
    !     deallocate(its_over)
    !     deallocate(file_idx)

    ! end subroutine outputD_combine_links
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine outputD_move_node_files
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Moves the node files to the swmm5 folder
    !     !%-----------------------------------------------------------------------------
    !     integer :: ii, fu, open_status
    !     character(len = 250) :: file_name, file_name_new
    !     character(len = 100) :: node_name
    !     character(len = 5)   :: str_image
    !     character(24) :: timestamp
    !     character(64) :: subroutine_name = 'outputD_move_node_files'
    !     !%--------------------------------------------------------------------------
    !     if (setting%Debug%File%output) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !     !% Get current image as a string
    !     write(str_image, '(i5.5)') this_image()

    !     do ii=1, size(node%P%have_output)

    !         if ((node%P%have_output(ii) > 0) .and. (node%P%have_output(ii) < nullvalueI )) then

    !             !% Open the node file
    !             !file_name = "debug_output/node/"//trim(node%names(node%P%have_output(ii))%str) &
    !             !    //"_"//trim(str_image)//".csv"
    !             file_name = trim(setting%File%debug_output_node_folder) &
    !                 //trim(node%names(node%P%have_output(ii))%str) &
    !                 //"_"//trim(str_image)//".csv"

    !             !file_name_new = "swmm5_output/node/"//trim(node%names(node%P%have_output(ii))%str)//".csv"
    !             file_name_new = trim(setting%File%swmm5_output_node_folder) &
    !                 //trim(node%names(node%P%have_output(ii))%str)//".csv"

    !             call rename(file_name, file_name_new)
    !         end if
    !     end do

    !     if (setting%Debug%File%output) &
    !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    ! end subroutine outputD_move_node_files
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine outputD_update_swmm_out()
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !%
    !     !%-----------------------------------------------------------------------------
    !     character(len = 250) :: fname
    !     character(len = 24) :: timestamp
    !     logical :: wrote_all_links = .false.
    !     logical :: wrote_all_nodes = .false.
    !     integer :: ii, rc, node_idx, link_idx
    !     integer, allocatable :: fus_nodes(:), fus_links(:)
    !     real(8) :: node_head, node_result
    !     real(8) :: link_flowrate, link_result
    !     real(8) :: timesecs
    !     character(64) :: subroutine_name = 'outputD_update_swmm_out'
    !     !%--------------------------------------------------------------------------
    !     if (setting%Debug%File%output) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !     if (this_image() == 1) then
    !         allocate(fus_nodes(size(node%P%have_output)))
    !         allocate(fus_links(size(link%P%have_output)))
    !         fus_links = nullvalueI
    !         fus_nodes = nullvalueI

    !         do while(.not. (wrote_all_links .and. wrote_all_nodes))

    !             if (size(node%P%have_output) < 1) then
    !                 wrote_all_nodes = .true.
    !             else
    !                 do ii=1, size(node%P%have_output)
    !                     !print *, 'node ',ii,  size(node%P%have_output),  node%P%have_output(ii)
    !                     node_idx = node%P%have_output(ii)
    !                     if ((node_idx > 0) .and. (node_idx < nullvalueI)) then
    !                         if (fus_nodes(ii) == nullvalueI) then
    !                             !% open files to process .out
    !                             !fname = "swmm5_output/node/"//trim(node%names(node_idx)%str)//".csv"
    !                             fname = trim(setting%File%swmm5_output_node_folder) &
    !                                 //trim(node%names(node_idx)%str)//".csv"
    !                             open(action='read', file=trim(fname), iostat=rc, newunit=fus_nodes(ii))
    !                             read(fus_nodes(ii), *, iostat = rc) timestamp
    !                         end if

    !                         read(fus_nodes(ii), "(A,2F10.8)", iostat = rc) timestamp, timesecs, node_head
    !                         if (rc /= 0) then
    !                             wrote_all_nodes = .true.
    !                             close(fus_nodes(ii))
    !                             !% Write line of .out
    !                             !exit
    !                         else
    !                             wrote_all_nodes = .false.
    !                             node_result = node_head - node%R(node_idx,nr_Zbottom)
    !                             !% stage values in .out
    !                             call interface_update_nodeResult(node_idx, api_output_node_depth, node_result)
    !                         end if
    !                     else
    !                         if (ii .eq. size(node%P%have_output)) then
    !                             wrote_all_nodes = .true.
    !                         end if
    !                     end if
    !                     if (ii .lt. size(node%P%have_output)) then
    !                         wrote_all_nodes = .false.
    !                     end if
    !                 end do
    !             end if

    !             if (size(link%P%have_output) < 1) then
    !                 wrote_all_links = .true.
    !             else
    !                 do ii=1, size(link%P%have_output)
    !                     !print *, 'link ',ii,  size(link%P%have_output),  link%P%have_output(ii)
    !                     link_idx = link%P%have_output(ii)
    !                     if ((link_idx > 0) .and. (link_idx < nullvalueI)) then
    !                         if (fus_links(ii) == nullvalueI) then
    !                             !% open files to process .out
    !                             !fname = "swmm5_output/link/"//trim(link%names(link_idx)%str)//".csv"
    !                             fname =  trim(setting%File%swmm5_output_link_folder) &
    !                                 //trim(link%names(link_idx)%str)//".csv"
    !                             open(action='read', file=trim(fname), iostat=rc, newunit=fus_links(ii))
    !                             read(fus_links(ii), *, iostat = rc) timestamp
    !                         end if

    !                         read(fus_links(ii), "(A,2F10.8)", iostat = rc) timestamp, timesecs, link_flowrate
    !                         if (rc /= 0) then
    !                             wrote_all_links = .true.
    !                             close(fus_links(ii))
    !                             !% Write line of .out
    !                             !exit
    !                         else
    !                             wrote_all_links = .false.
    !                             link_result = link_flowrate
    !                             !% stage values in .out
    !                             call interface_update_linkResult(link_idx, api_output_link_flow, link_result)
    !                         end if
    !                     else
    !                         if (ii .eq. size(link%P%have_output)) then
    !                             wrote_all_links = .true.
    !                         end if
    !                     end if
    !                     if (ii .lt. size(link%P%have_output)) then
    !                         wrote_all_links = .false.
    !                     end if
    !                 end do
    !             end if
    !             call interface_write_output_line(timesecs)
    !         end do

    !         deallocate(fus_links)
    !         deallocate(fus_nodes)

    !         ! do ii = 1, size(link%P%have_output)
    !         !     link_idx = link%P%have_output(ii)
    !         !     call interface_export_link_results(link_idx)
    !         ! end do
    !     end if

    !     if (setting%Debug%File%output) &
    !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    ! end subroutine outputD_update_swmm_out
!%
!%==========================================================================
!%==========================================================================
!%    
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine outputML_store_binary_output_filenames (nWritten, file_name)
        !%-----------------------------------------------------------------------------
        !% Description
        !% Stores the binary output filenames in memory or writes to a file
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: nWritten !% # of file writing about to begin
        character(len=256), intent(in) :: file_name !% name of file to be written
        integer :: kk
        integer :: ios
        integer, pointer ::  fnunit  !% pointers for file unit numbers
        character(64)    :: subroutine_name = 'outputML_store_binary_output_filenames'
        !%-----------------------------------------------------------------------------
        fnunit   => setting%File%UnitNumber%outputML_filename_file

        !% --- store the file names
        !% default is into memory until StoredFileNames limit is exceeded.
        !% then write all to a file (name stored in output_binary_filenames)
        !% and clear the memory to prevent confusion
        if (setting%Output%UseFileNameFile) then
            !% --- forced to use filename file
            if (nWritten == 1) then
                !% --- open the output filename file
                open(unit=fnunit, &
                file=trim(setting%File%outputML_filename_file), &
                form='formatted',action='write',status='new', iostat=ios)
                if (ios /= 0) then
                    write(*,"(A)") 'ERROR (CODE) file could not be opened for writing...'
                    write(*,"(A)") 'filename is ...'
                    write(*,"(A)") trim(file_name)
                    stop
                end if
                !% --- write the filename
                write(fnunit,"(A)") trim(file_name)
                !print *, '1 writing in filename file ',trim(file_name)
            else
                !% --- file stays open, so just write the next file name when ii > 1
                write(fnunit,"(A)") trim(file_name)
                !print *, '2 writing in filename file ',trim(file_name)
            end if
        else
            !% --- if not forced to use filename file, use local memory storage
            if (nWritten .le. setting%Output%StoredFileNames ) then
                output_binary_filenames(nWritten) = trim(file_name)
            elseif (nWritten == (setting%Output%StoredFileNames+1) ) then
                !% --- when allocated memory is exceeded...
                !% --- open the filenames file to store all the filenames
                open(unit=fnunit, &
                    file=trim(setting%File%outputML_filename_file), &
                    form='formatted',action='write',status='new', iostat=ios)
                    if (ios /= 0) then
                        write(*,"(A)") 'ERROR (CODE) file could not be opened for writing...'
                        write(*,"(A)") 'filename is ...'
                        write(*,"(A)") trim(file_name)
                        stop
                    end if
                !% --- write the prior filenames from memory to the file and the delete
                do kk=1,setting%Output%StoredFileNames
                    write(fnunit,"(A)") trim(output_binary_filenames(kk))
                    !print *, '3 writing in filename file ',trim(output_binary_filenames(kk))
                    output_binary_filenames(kk) = ""
                end do
                !% --- set flag to recover file names by reading
                setting%Output%UseFileNameFile = .true.
                write(fnunit,"(A)") trim(file_name)
                !print *, '4 writing in filename file ',trim(file_name)
            else
                write(fnunit,"(A)") trim(file_name)
                !print *, '5 writing in filename file ',trim(file_name)
            end if
        end if

        !% note, do not close(fnunit) as it may need to be written into on later steps

    end subroutine outputML_store_binary_output_filenames
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine outputML_get_all_output_binary_filenames (nWritten)
        !%------------------------------------------------------------------
        !% Description
        !% gets binary output filenames in memory from a file file
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: nWritten  !% number of files that were written
            integer             :: allocation_status, ios, ii
            integer, pointer    :: fnunit
            logical             :: isopen, doesexist
            character(len=99)   :: emsg
            character(64)       :: subroutine_name = 'outputML_get_output_binary_filenames'
        !%-------------------------------------------------------------------

        !% --- create filename storage that is large enough for all the files
        allocate(output_binary_filenames_all(nWritten), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'output_binary_filenames_all')

        if (setting%Output%UseFileNameFile) then
            !% --- reading in the file of filenames
            fnunit   => setting%File%UnitNumber%outputML_filename_file
            inquire(UNIT=fnunit,EXIST=doesexist)
            if (.not. doesexist) then 
                write(*,"(A)") 'ERROR (file), expected the file of output file names to already exist.'
                write(*,"(A)") 'but it cannot be found. The filename is'
                write(*,"(A)") trim(setting%File%outputML_filename_file)
                write(*,"(A,i6)") 'and the unit number is ',fnunit 
                stop 339182
            end if     
            inquire(UNIT=fnunit,OPENED=isopen)

            !% --- if previously open for writing, we want to close to switch to reading
            if (isopen) close(fnunit)
            open(unit=fnunit, &
                file=trim(setting%File%outputML_filename_file), &
                form='formatted',action='read', iostat=ios)
            if (ios /= 0) then
                write(*,"(A)") 'ERROR (file): iostat /=0 for open() file named ....'
                write(*,"(A)") trim(setting%File%outputML_filename_file)
                write(*,"(A)") '... file is the outputML_filename_file ...'
                write(*,"(A,i5)") '... iostat value = ',ios
                stop 89075
            end if
            rewind(unit=fnunit)
            do ii=1,nWritten
                read(fnunit,"(A)") output_binary_filenames_all(ii)
            end do
            close(fnunit)
        else
            output_binary_filenames_all(:) = output_binary_filenames(1:nWritten)
        end if

    end subroutine outputML_get_all_output_binary_filenames
!%
!%==========================================================================
!% END MODULE
!%==========================================================================
!%
end module output
