!
!% module pack_mask_arrays
!
!
!==========================================================================
!
!==========================================================================
!
module pack_mask_arrays

    use define_globals
    use define_indexes
    use define_keys
    use define_settings
    use utility_crash
    
    implicit none

    private

    public :: pack_mask_arrays_all
    public :: pack_dynamic_arrays
    public :: pack_nodes
    public :: pack_bc
    public :: pack_element_outputML
    public :: pack_face_outputML
    public :: pack_small_and_zero_depth_elements

contains
!
!==========================================================================
! PUBLIC
!==========================================================================
!
    subroutine pack_mask_arrays_all()
        !--------------------------------------------------------------------------
        !% set all the static packs and masks
        !--------------------------------------------------------------------------
        integer :: ii
        character(64) :: subroutine_name = 'pack_mask_arrays_all'
        !--------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%pack_mask_arrays) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        
        !call mask_faces_whole_array_static()
        
        call pack_geometry_alltm_elements()
        
        call pack_geometry_etm_elements()
        
        call pack_geometry_ac_elements()
        
        call pack_nongeometry_static_elements()
        
        call pack_nongeometry_dynamic_elements()
        
        call pack_static_interior_faces()
        
        call pack_static_shared_faces()
        
        call pack_dynamic_interior_faces()
        
        call pack_dynamic_shared_faces()

        call pack_small_and_zero_depth_elements(ALLtm)
        call pack_small_and_zero_depth_elements(ETM)

        if (setting%Debug%File%initial_condition) then
            !% only using the first processor to print results
            if (this_image() == 1) then

                do ii = 1,num_images()
                   print*, '----------------------------------------------------'
                   print*, 'image = ', ii
                   print*, '..........packed local element indexes of...........'
                   print*, elemP(:,ep_ALLtm)[ii], 'all ETM, AC elements'
                   print*, elemP(:,ep_ETM)[ii], 'all ETM elements'
                   print*, elemP(:,ep_CC_ETM)[ii], 'all CC elements that are ETM'
                   print*, elemP(:,ep_Diag)[ii], 'all diagnostic elements'
                   print*, '.................face logicals......................'
                   print*, faceP(:,fp_all)[ii], 'all the interior faces'
                   print*, facePS(:,fp_all)[ii], 'all the shared faces'
                   ! call execute_command_line('')
                end do

            end if
        end if

        if (setting%Debug%File%pack_mask_arrays) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine pack_mask_arrays_all
!
!==========================================================================
!==========================================================================
!
    subroutine pack_dynamic_arrays()
        !--------------------------------------------------------------------------
        !% set all the dynamic packs and masks
        !--------------------------------------------------------------------------
        character(64) :: subroutine_name = 'pack_dynamic_arrays'
        !--------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%pack_mask_arrays) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        call pack_geometry_etm_elements()
        call pack_geometry_ac_elements()
        call pack_nongeometry_dynamic_elements()
        call pack_dynamic_interior_faces()
        call pack_dynamic_shared_faces()

        if (setting%Debug%File%pack_mask_arrays) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine pack_dynamic_arrays
!    
!==========================================================================    
!==========================================================================
!
    subroutine pack_nodes()
        !--------------------------------------------------------------------------
        !% This allocates and packs the node data in the arrays of node%P.
        !% With this approach using the P type, each of the arrays on the images
        !% are allocated to the size needed.
        !--------------------------------------------------------------------------
        character(64)    :: subroutine_name = 'pack_nodes'
        !--------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%pack_mask_arrays) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        N_flowBC = count(node%YN(:,nYN_has_inflow) .and. &
                        (node%I(:,ni_P_image) == this_image()))

        if (N_flowBC > 0) then
            allocate(node%P%have_flowBC(N_flowBC))
            node%P%have_flowBC = pack(node%I(:,ni_idx), &
                node%YN(:,nYN_has_inflow) .and. (node%I(:,ni_P_image) == this_image()))
        end if

        !% HACK -- this assumes that a head BC is always a downstream BC.
        N_headBC = count((node%I(:, ni_node_type) == nBCdn) .and. &
                        (node%I(:,ni_P_image) == this_image()))

        if (N_headBC > 0) then
            allocate(node%P%have_headBC(N_headBC))
            node%P%have_headBC = pack(node%I(:,ni_idx), &
            (node%I(:, ni_node_type) == nBCdn) .and. &
            (node%I(:,ni_P_image) == this_image()))
        end if
        if (setting%Debug%File%pack_mask_arrays) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine pack_nodes
!
!==========================================================================
!==========================================================================
!
    subroutine pack_bc
        !--------------------------------------------------------------------------
        !% 
        !--------------------------------------------------------------------------
        integer :: psize
        character(64) :: subroutine_name = 'pack_bc'
        !--------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%pack_mask_arrays) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% BC packs
        !% zero out the number of upBC to get a new count of how many is in a given partition
        N_nBCup = 0
        if (N_flowBC > 0) then
            N_nBCup = count(BC%flowI(:, bi_category) == BCup)
            if (N_nBCup > 0) then
                allocate(BC%P%BCup(N_nBCup))
                BC%P%BCup = pack(BC%flowI(:, bi_idx), BC%flowI(:, bi_category) == BCup)
                !% Face packs
                npack_faceP(fp_BCup) = N_nBCup
                faceP(1:N_nBCup,fp_BCup) = BC%flowI(BC%P%BCup, bi_face_idx)
            end if

            N_nBClat = count(BC%flowI(:, bi_category) == BClat)
            if (N_nBClat > 0) then
                allocate(BC%P%BClat(N_nBClat))
                BC%P%BClat = pack(BC%flowI(:, bi_idx), BC%flowI(:, bi_category) == BClat)
                !% Elem Packs
                npack_elemP(ep_BClat) = N_nBClat
                elemP(1:N_nBClat,ep_BClat) = BC%flowI(BC%P%BClat, bi_elem_idx)
            end if
        end if

        !% BC packs
        !% zero out the number of dnBC to get a new count of how many is in a given partition
        N_nBCdn = 0
        if (N_headBC > 0) then
            N_nBCdn = count(BC%headI(:, bi_category) == BCdn)
            if (N_nBCdn > 0) then
                allocate(BC%P%BCdn(N_nBCdn))
                BC%P%BCdn = pack(BC%headI(:, bi_idx), BC%headI(:, bi_category) == BCdn)
                !% Face packs
                npack_faceP(fp_BCdn) = N_nBCdn
                faceP(1:N_nBCdn,fp_BCdn) = BC%headI(BC%P%BCdn, bi_face_idx)
            end if
        end if
        if (setting%Debug%File%pack_mask_arrays) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine pack_bc
!
!==========================================================================
!==========================================================================
!
    subroutine pack_element_outputML ()
        !%-----------------------------------------------------------------
        !% Description:
        !% creates a pack for elemR, elemI, elemYN for output elements
        !%-----------------------------------------------------------------
        !% Declarations:
            logical, pointer :: isElemOut(:)
            integer, pointer :: eIdx(:), ptype, npack
            character(64) :: subroutine_name = 'pack_element_outputML'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Output%Report%suppress_MultiLevel_Output) return
            !if (crashYN) return
        !%------------------------------------------------------------------
        !% Aliases
            eIdx => elemI(:,ei_Lidx)
            !% logical control on output for each element
            isElemOut => elemYN(:,eYN_isOutput)
            ptype => col_elemP(ep_Output_Elements)
            npack => npack_elemP(ptype)
        !%------------------------------------------------------------------
        !% count the output elements to be packed
        npack = count(isElemOut)

        !% output the true element indexes into a pack
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,isElemOut)
        else 
            !% continue     
        end if    

    end subroutine pack_element_outputML        
!
!==========================================================================
!==========================================================================
!    
    subroutine pack_face_outputML ()
        !%-----------------------------------------------------------------
        !% Description:
        !% packed arrays faceR, faceI, faceYN output elements
        !%-----------------------------------------------------------------
        !% Declarations:
            logical, pointer :: isFaceOut(:)
            integer, pointer :: fIdx(:), ptype, npack
            character(64) :: subroutine_name = 'pack_face_outputML'
        !%-----------------------------------------------------------------
        !% Preliminaries
            if (setting%Output%Report%suppress_MultiLevel_Output) return
            !if (crashYN) return
        !%-----------------------------------------------------------------
        !% Aliases:
            fIdx => faceI(:,fi_Lidx)
            !% logical control on output for each element
            isFaceOut => faceYN(:,fYN_isFaceOut)
            !% pointers for storage of type and npack storage
            ptype => col_faceP(fp_Output_Faces)
            npack => npack_faceP(ptype)
        !%-----------------------------------------------------------------
        !% count the output faces to be packed
        npack = count(isFaceOut)

        !% output the true element indexes into a pack
        if (npack > 0) then
            faceP(1:npack,ptype) = pack(fIdx,isFaceOut)
        else 
            !% continue    
        end if  

    end subroutine pack_face_outputML
!
!==========================================================================
! PRIVATE
!==========================================================================
!
    ! subroutine mask_faces_whole_array_static()
    !     !--------------------------------------------------------------------------
    !     !% find all the faces except boundary and null faces
    !     !--------------------------------------------------------------------------
    !     integer, pointer :: mcol
    !     character(64) :: subroutine_name = 'mask_faces_whole_array_static'
    !     !--------------------------------------------------------------------------
    !     !if (crashYN) return
    !     if (setting%Debug%File%pack_mask_arrays) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !     mcol => col_faceM(fm_all)

    !     faceM(:,mcol) = ( &
    !         (faceI(:,fi_BCtype) == doesnotexist) &
    !         .and. &
    !         (.not. faceYN(:,fYN_isnull))  &
    !         .and. &
    !         (.not. faceYN(:,fYN_isSharedFace)) &
    !         )

    !     if (setting%Debug%File%pack_mask_arrays) &
    !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine mask_faces_whole_array_static
!
!==========================================================================
!==========================================================================
!
    subroutine pack_geometry_alltm_elements()
        !--------------------------------------------------------------------------
        !% packed arrays for geometry types in elemPGalltm
        !--------------------------------------------------------------------------
        integer, pointer :: ptype, npack, eIDx(:)
        character(64) :: subroutine_name = 'pack_geometry_alltm_elements'
        !--------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%pack_mask_arrays) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        eIdx => elemI(:,ei_Lidx)

        !% --- rectangular channels 
        ptype => col_elemPGalltm(epg_CC_rectangular_nonsurcharged)
        npack => npack_elemPGalltm(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                ( &
                    (elemI(:,ei_geometryType) == rectangular) &
                ) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))

        if (npack > 0) then
            elemPGalltm(1:npack, ptype) = pack(eIdx, &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                ( &
                    (elemI(:,ei_geometryType) == rectangular) &
                ) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))
        end if

        !% --- rectangular conduits 
        ptype => col_elemPGalltm(epg_CC_rectangular_closed_nonsurcharged)
        npack => npack_elemPGalltm(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                ( &
                    (elemI(:,ei_geometryType) == rectangular_closed) &
                ) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))

        if (npack > 0) then
            elemPGalltm(1:npack, ptype) = pack(eIdx, &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                ( &
                    (elemI(:,ei_geometryType) == rectangular_closed) &
                ) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))
        end if

        !% --- trapezoidal channels 
        ptype => col_elemPGalltm(epg_CC_trapezoidal_nonsurcharged)
        npack => npack_elemPGalltm(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == trapezoidal) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))

        if (npack > 0) then
            elemPGalltm(1:npack, ptype) = pack(eIdx, &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == trapezoidal) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))
        end if

        !% --- triangular channels 
        ptype => col_elemPGalltm(epg_CC_triangular_nonsurcharged)
        npack => npack_elemPGalltm(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == triangular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))

        if (npack > 0) then
            elemPGalltm(1:npack, ptype) = pack(eIdx, &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == triangular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))
        end if

        !% --- rectangular triangular channels 
        ptype => col_elemPGalltm(epg_CC_rectangular_triangular_nonsurcharged)
        npack => npack_elemPGalltm(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == rect_triang) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))

        if (npack > 0) then
            elemPGalltm(1:npack, ptype) = pack(eIdx, &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == rect_triang) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))
        end if

        !% --- circular conduits 
        ptype => col_elemPGalltm(epg_CC_circular_nonsurcharged)
        npack => npack_elemPGalltm(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == circular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))

        if (npack > 0) then
            elemPGalltm(1:npack, ptype) = pack(eIdx, &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == circular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))
        end if

        !% --- irregular channels
        ptype => col_elemPGalltm(epg_CC_irregular_nonsurcharged)
        npack => npack_elemPGalltm(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == irregular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))

        if (npack > 0) then
            elemPGalltm(1:npack, ptype) = pack(eIdx, &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == irregular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))
        end if

        !% --- junction main with functional geometry relationship
        ptype => col_elemPGalltm(epg_JM_functionalStorage_nonsurcharged)
        npack => npack_elemPGalltm(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemSI(:,esi_JunctionMain_Type) == FunctionalStorage) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))

        if (npack > 0) then
            elemPGalltm(1:npack, ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemSI(:,esi_JunctionMain_Type) == FunctionalStorage) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))
        end if

        !% --- junction main with tabular geometry relationship
        ptype => col_elemPGalltm(epg_JM_tabularStorage_nonsurcharged)
        npack => npack_elemPGalltm(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemSI(:,esi_JunctionMain_Type) == TabularStorage) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))

        if (npack > 0) then
            elemPGalltm(1:npack, ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemSI(:,esi_JunctionMain_Type) == TabularStorage) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))
        end if

        !% --- junction main with artificial storage relationship -- for ALL tm
        ptype => col_elemPGalltm(epg_JM_impliedStorage_nonsurcharged)
        npack => npack_elemPGalltm(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemSI(:,esi_JunctionMain_Type) == ImpliedStorage) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))
      

        if (npack > 0) then
            elemPGalltm(1:npack, ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemSI(:,esi_JunctionMain_Type) == ImpliedStorage) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))
        end if

        if (setting%Debug%File%pack_mask_arrays) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine pack_geometry_alltm_elements
!
!==========================================================================
!==========================================================================
!
    subroutine pack_geometry_ac_elements()
        !--------------------------------------------------------------------------
        !% packed arrays for geometry types for AC elements
        !--------------------------------------------------------------------------
        integer, pointer :: ptype, npack, eIDx(:)
        character(64) :: subroutine_name = 'pack_geometry_alltm_elements'
        !--------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%pack_mask_arrays) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        eIdx => elemI(:,ei_Lidx)

        !% rectangular channels 
        ptype => col_elemPGac(epg_CC_rectangular_nonsurcharged)
        npack => npack_elemPGac(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                ( &
                    (elemI(:,ei_geometryType) == rectangular) &
                ) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )

        if (npack > 0) then
            elemPGac(1:npack, ptype) = pack(eIdx, &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                ( &
                    (elemI(:,ei_geometryType) == rectangular) &
                ) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged))&
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )
        end if

        !% rectangular conduits 
        ptype => col_elemPGac(epg_CC_rectangular_closed_nonsurcharged)
        npack => npack_elemPGac(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                ( &
                    (elemI(:,ei_geometryType) == rectangular_closed) &
                ) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )

        if (npack > 0) then
            elemPGac(1:npack, ptype) = pack(eIdx, &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                ( &
                    (elemI(:,ei_geometryType) == rectangular_closed) &
                ) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged))&
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )
        end if

        !% trapezoidal channels 
        ptype => col_elemPGac(epg_CC_trapezoidal_nonsurcharged)
        npack => npack_elemPGac(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == trapezoidal) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )

        if (npack > 0) then
            elemPGac(1:npack, ptype) = pack(eIdx, &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == trapezoidal) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged))&
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )
        end if

        !% triangular channels 
        ptype => col_elemPGac(epg_CC_triangular_nonsurcharged)
        npack => npack_elemPGac(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == triangular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )

        if (npack > 0) then
            elemPGac(1:npack, ptype) = pack(eIdx, &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == triangular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged))&
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )
        end if

       !% irregular channels 
        ptype => col_elemPGac(epg_CC_irregular_nonsurcharged)
        npack => npack_elemPGac(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == irregular) &
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )

        if (npack > 0) then
            elemPGac(1:npack, ptype) = pack(eIdx, &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == irregular) &
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )
        end if

        
        !% rectangular triangular channels 
        ptype => col_elemPGac(epg_CC_rectangular_triangular_nonsurcharged)
        npack => npack_elemPGac(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == rect_triang) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )

        if (npack > 0) then
            elemPGac(1:npack, ptype) = pack(eIdx, &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == rect_triang) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged))&
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )
        end if

        !% circular conduits 
        ptype => col_elemPGac(epg_CC_circular_nonsurcharged)
        npack => npack_elemPGac(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == circular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )

        if (npack > 0) then
            elemPGac(1:npack, ptype) = pack(eIdx, &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == circular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged))&
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )
        end if

        !% junction main with functional geometry relationship
        ptype => col_elemPGac(epg_JM_functionalStorage_nonsurcharged)
        npack => npack_elemPGac(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemSI(:,esi_JunctionMain_Type) == FunctionalStorage) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )

        if (npack > 0) then
            elemPGac(1:npack, ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemSI(:,esi_JunctionMain_Type) == FunctionalStorage) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )
        end if

        !% junction main with functional geometry relationship
        ptype => col_elemPGac(epg_JM_tabularStorage_nonsurcharged)
        npack => npack_elemPGac(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemSI(:,esi_JunctionMain_Type) == TabularStorage) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )

        if (npack > 0) then
            elemPGac(1:npack, ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemSI(:,esi_JunctionMain_Type) == TabularStorage) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )
        end if

        !% junction main with artificial storage relationship -- for AC
        ptype => col_elemPGac(epg_JM_impliedStorage_nonsurcharged)
        npack => npack_elemPGac(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemSI(:,esi_JunctionMain_Type) == ImpliedStorage) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( elemI(:,ei_tmType) == AC) &
                )
      

        if (npack > 0) then
            elemPGac(1:npack, ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemSI(:,esi_JunctionMain_Type) == ImpliedStorage) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( elemI(:,ei_tmType) == AC) &
                )
        end if       

        if (setting%Debug%File%pack_mask_arrays) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine pack_geometry_ac_elements
!
!==========================================================================
!==========================================================================
!
    subroutine pack_geometry_etm_elements()
        !%--------------------------------------------------------------------------
        !% packed arrays for geometry types that are explicit time marching (ETM)
        !%--------------------------------------------------------------------------
            integer, pointer :: ptype, npack, eIDx(:)
            character(64) :: subroutine_name = 'pack_geometry_etm_elements'
        !%--------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%pack_mask_arrays) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        eIdx => elemI(:,ei_Lidx)

        !% --- rectangular channels nonsurcharged
        ptype => col_elemPGetm(epg_CC_rectangular_nonsurcharged)
        npack => npack_elemPGetm(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                ( &
                    (elemI(:,ei_geometryType) == rectangular) &
                ) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )

        if (npack > 0) then
            elemPGetm(1:npack, ptype) = pack(eIdx, &
                 (elemI(:,ei_elementType) == CC)  &
                .and. &
                ( &
                    (elemI(:,ei_geometryType) == rectangular) &
                ) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )
        end if


        !% --- rectangular conduits nonsurcharged
        ptype => col_elemPGetm(epg_CC_rectangular_closed_nonsurcharged)
        npack => npack_elemPGetm(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                ( &
                    (elemI(:,ei_geometryType) == rectangular_closed) &
                ) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )

        if (npack > 0) then
            elemPGetm(1:npack, ptype) = pack(eIdx, &
                 (elemI(:,ei_elementType) == CC)  &
                .and. &
                ( &
                    (elemI(:,ei_geometryType) == rectangular_closed) &
                ) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )
        end if

        !% --- trapezoidal channels, conduits nonsurcharged
        ptype => col_elemPGetm(epg_CC_trapezoidal_nonsurcharged)
        npack => npack_elemPGetm(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == trapezoidal) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )

        if (npack > 0) then
            elemPGetm(1:npack, ptype) = pack(eIdx, &
                    (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == trapezoidal) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )
        end if

        !% --- triangular channels, conduits nonsurcharged
        ptype => col_elemPGetm(epg_CC_triangular_nonsurcharged)
        npack => npack_elemPGetm(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == triangular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )

        if (npack > 0) then
            elemPGetm(1:npack, ptype) = pack(eIdx, &
                    (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == triangular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )
        end if

        !% --- rectangular triangular channels, conduits 
        ptype => col_elemPGetm(epg_CC_rectangular_triangular_nonsurcharged)
        npack => npack_elemPGetm(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == rect_triang) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )

        if (npack > 0) then
            elemPGetm(1:npack, ptype) = pack(eIdx, &
                    (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == rect_triang) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )
        end if

        !% --- circular conduits, channels
        ptype => col_elemPGetm(epg_CC_circular_nonsurcharged)
        npack => npack_elemPGetm(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_geometryType) == circular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )

        if (npack > 0) then
            elemPGetm(1:npack, ptype) = pack(eIdx, &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == circular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged))&
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )
        end if

        !% --- irregular channels nonsurcharged
        ptype => col_elemPGetm(epg_CC_irregular_nonsurcharged)
        npack => npack_elemPGalltm(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == irregular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )

        if (npack > 0) then
            elemPGalltm(1:npack, ptype) = pack(eIdx, &
                (elemI(:,ei_elementType) == CC)  &
                .and. &
                (elemI(:,ei_geometryType) == irregular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )
        end if

        !% --- junction main with functional geometry relationship, nonsurcharged
        ptype => col_elemPGetm(epg_JM_functionalStorage_nonsurcharged)
        npack => npack_elemPGetm(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemSI(:,esi_JunctionMain_Type) == FunctionalStorage) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( elemI(:,ei_tmType) == ETM) &
                )

        if (npack > 0) then
            elemPGetm(1:npack, ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemSI(:,esi_JunctionMain_Type) == FunctionalStorage) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( elemI(:,ei_tmType) == ETM) &
                )
        end if

        !% --- junction main with functional geometry relationship, nonsurcharged
        ptype => col_elemPGetm(epg_JM_tabularStorage_nonsurcharged)
        npack => npack_elemPGetm(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemSI(:,esi_JunctionMain_Type) == TabularStorage) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( elemI(:,ei_tmType) == ETM) &
                )

        if (npack > 0) then
            elemPGetm(1:npack, ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemSI(:,esi_JunctionMain_Type) == TabularStorage) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( elemI(:,ei_tmType) == ETM) &
                )
        end if


        !% --- junction main with artificial storage relationship -- for ETM, nonsurcharged
        ptype => col_elemPGetm(epg_JM_impliedStorage_nonsurcharged)
        npack => npack_elemPGetm(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemSI(:,esi_JunctionMain_Type) == ImpliedStorage) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( elemI(:,ei_tmType) == ETM) &
                )
      

        if (npack > 0) then
            elemPGetm(1:npack, ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemSI(:,esi_JunctionMain_Type) == ImpliedStorage) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( elemI(:,ei_tmType) == ETM) &
                )
        end if 

        if (setting%Debug%File%pack_mask_arrays) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine pack_geometry_etm_elements
!
!==========================================================================
!==========================================================================
!
    subroutine pack_nongeometry_static_elements()
        !--------------------------------------------------------------------------
        !
        !% packed arrays for non geometry static elements
        !%
        !% HACK: Note that in this approach, an element that is diagnostic in Q and
        !% time march in H would appear in both the ep_Diag and ep_ALLtm packed arrays.
        !% Is this something we want to allow? It seems like we might need to require
        !% that an element that is diagnostic in Q must be either diagnostic in H or
        !% doesnotexist. Similarly, an element that is time march in H must be time-march
        !% in Q or doesnotexist
        !
        !--------------------------------------------------------------------------
        integer, pointer :: ptype, npack, eIDx(:)
        character(64) :: subroutine_name = 'pack_nongeometry_static_elements'
        !--------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%pack_mask_arrays) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        eIdx => elemI(:,ei_Lidx)

        !% ep_ALLtm
        !% - all elements that have a time march
        ptype => col_elemP(ep_ALLtm)
        npack => npack_elemP(ptype)
        npack = count( &
                (elemI(:,ei_QeqType) == time_march ) &
                .or. &
                (elemI(:,ei_HeqType) == time_march ) &
                )

        if (npack > 0) then
            elemP(1:npack,ptype) = pack( eIdx, &
                (elemI(:,ei_QeqType) == time_march ) &
                .or. &
                (elemI(:,ei_HeqType) == time_march ) &
                )
        end if

        !% ep_CC_ALLtm
        !% - all time march elements that are CC
        ptype => col_elemP(ep_CC_ALLtm)
        npack => npack_elemP(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_QeqType) == time_march ) &
                    .or. &
                    (elemI(:,ei_HeqType) == time_march ) &
                ) &
                .and. &
                (elemI(:,ei_elementType) == CC) &
                )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack( eIdx, &
                ( &
                    (elemI(:,ei_QeqType) == time_march ) &
                    .or. &
                    (elemI(:,ei_HeqType) == time_march ) &
                ) &
                .and. &
                (elemI(:,ei_elementType) == CC) &
                )
        end if

        !% ep_CCJB_ALLtm
        !% - all time march elements that are CC or JB
        ptype => col_elemP(ep_CCJB_ALLtm)
        npack => npack_elemP(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_QeqType) == time_march ) &
                    .or. &
                    (elemI(:,ei_HeqType) == time_march ) &
                ) &
                .and. &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB) &
                ))
        if (npack > 0) then
            elemP(1:npack,ptype) = pack( eIdx, &
                ( &
                    (elemI(:,ei_QeqType) == time_march ) &
                    .or. &
                    (elemI(:,ei_HeqType) == time_march ) &
                ) &
                .and. &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB) &
                ))
        end if

        !% ep_Diag
        !% - all elements that are diagnostic
        ptype => col_elemP(ep_Diag)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemI(:,ei_QeqType) == diagnostic ) &
                .or. &
                (elemI(:,ei_HeqType) == diagnostic ) &
                )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack( eIdx, &
                (elemI(:,ei_QeqType) == diagnostic ) &
                .or. &
                (elemI(:,ei_HeqType) == diagnostic ) &
                )
        end if

        !print *, 'JM'
        !% ep_JM
        !% - all elements that are JM
        ptype => col_elemP(ep_JM)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemI(:,ei_elementType) == JM))
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx, &
                (elemI(:,ei_elementType) == JM))
        end if

        !% ep_JM_ALLtm
        !% - all junction main elements that are time march
        ptype => col_elemP(ep_JM_ALLtm)
        npack => npack_elemP(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_QeqType) == time_march ) &
                    .or. &
                    (elemI(:,ei_HeqType) == time_march ) &
                ) &
                .and. &
                (elemI(:,ei_elementType) == JM) &
                )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack( eIdx, &
                ( &
                    (elemI(:,ei_QeqType) == time_march ) &
                    .or. &
                    (elemI(:,ei_HeqType) == time_march ) &
                ) &
                .and. &
                (elemI(:,ei_elementType) == JM) &
                )
        end if

        !% ep_JB_ALLtm
        !% - all junction main elements that are time march
        ptype => col_elemP(ep_JB_ALLtm)
        npack => npack_elemP(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_QeqType) == time_march ) &
                    .or. &
                    (elemI(:,ei_HeqType) == time_march ) &
                ) &
                .and. &
                (elemI(:,ei_elementType) == JB) &
                )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack( eIdx, &
                ( &
                    (elemI(:,ei_QeqType) == time_march ) &
                    .or. &
                    (elemI(:,ei_HeqType) == time_march ) &
                ) &
                .and. &
                (elemI(:,ei_elementType) == JB) &
                )
        end if

        !% ep_JB_DownStreamJB
        !% - all downstream JB elements
        ptype => col_elemP(ep_JB_DownStreamJB)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemI(:,ei_elementType) == JB) &
                .and. &
                (elemYN(:,eYN_isDownstreamJB)) &
                )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack( eIdx, &
                ( &
                (elemI(:,ei_elementType) == JB) &
                .and. &
                (elemYN(:,eYN_isDownstreamJB)) &
                ))
        endif

        !% ep_CC_DownstreamJbAdjacent
        !% - all CC element downstream of a JB
        ptype => col_elemP(ep_CC_DownstreamJbAdjacent)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemYN(:,eYN_isElementDownstreamOfJB)) &
                )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack( eIdx, &
                ( &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemYN(:,eYN_isElementDownstreamOfJB)) &
                ))
        endif

        !% ep_CC_Closed_Elements
        !% - all the closed time-marching elements
        ptype => col_elemP(ep_CC_Closed_Elements)
        npack => npack_elemP(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                ) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ) &
                .and. &
                ( &
                    (elemI(:,ei_geometryType) == circular) &
                    .or. &
                    (elemI(:,ei_geometryType) == filled_circular) &
                    .or. &
                    (elemI(:,ei_geometryType) == rectangular_closed) &
                    .or. &
                    (elemI(:,ei_geometryType) == horiz_ellipse) &
                    .or. &
                    (elemI(:,ei_geometryType) == arch) &
                    .or. &
                    (elemI(:,ei_geometryType) == eggshaped) &
                    .or. &
                    (elemI(:,ei_geometryType) == horseshoe) &
                    .or. &
                    (elemI(:,ei_geometryType) == gothic) &
                    .or. &
                    (elemI(:,ei_geometryType) == catenary) &
                    .or. &
                    (elemI(:,ei_geometryType) == semi_elliptical) &
                    .or. &
                    (elemI(:,ei_geometryType) == basket_handle) &
                    .or. &
                    (elemI(:,ei_geometryType) == semi_circular) &
                    .or. &
                    (elemI(:,ei_geometryType) == custom) &
                    .or. &
                    (elemI(:,ei_geometryType) == force_main) &
                ) )

        if (npack > 0) then
            elemP(1:npack, ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                ) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ) &
                .and. &
                ( &
                    (elemI(:,ei_geometryType) == circular) &
                    .or. &
                    (elemI(:,ei_geometryType) == filled_circular) &
                    .or. &
                    (elemI(:,ei_geometryType) == rectangular_closed) &
                    .or. &
                    (elemI(:,ei_geometryType) == horiz_ellipse) &
                    .or. &
                    (elemI(:,ei_geometryType) == arch) &
                    .or. &
                    (elemI(:,ei_geometryType) == eggshaped) &
                    .or. &
                    (elemI(:,ei_geometryType) == horseshoe) &
                    .or. &
                    (elemI(:,ei_geometryType) == gothic) &
                    .or. &
                    (elemI(:,ei_geometryType) == catenary) &
                    .or. &
                    (elemI(:,ei_geometryType) == semi_elliptical) &
                    .or. &
                    (elemI(:,ei_geometryType) == basket_handle) &
                    .or. &
                    (elemI(:,ei_geometryType) == semi_circular) &
                    .or. &
                    (elemI(:,ei_geometryType) == custom) &
                    .or. &
                    (elemI(:,ei_geometryType) == force_main) &
                ) )
        end if

        !% ep_Closed_Elements_JB
        !% - all the closed time-marching elements
        ptype => col_elemP(ep_Closed_Elements_JB)
        npack => npack_elemP(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == JB) &
                ) &
                .and. &
                ( &
                    (elemSI(:,esi_JunctionBranch_Exists) == 1) &
                ) &
                .and. &
                ( &
                    (elemI(:,ei_geometryType) == circular) &
                    .or. &
                    (elemI(:,ei_geometryType) == filled_circular) &
                    .or. &
                    (elemI(:,ei_geometryType) == rectangular_closed) &
                    .or. &
                    (elemI(:,ei_geometryType) == horiz_ellipse) &
                    .or. &
                    (elemI(:,ei_geometryType) == arch) &
                    .or. &
                    (elemI(:,ei_geometryType) == eggshaped) &
                    .or. &
                    (elemI(:,ei_geometryType) == horseshoe) &
                    .or. &
                    (elemI(:,ei_geometryType) == gothic) &
                    .or. &
                    (elemI(:,ei_geometryType) == catenary) &
                    .or. &
                    (elemI(:,ei_geometryType) == semi_elliptical) &
                    .or. &
                    (elemI(:,ei_geometryType) == basket_handle) &
                    .or. &
                    (elemI(:,ei_geometryType) == semi_circular) &
                    .or. &
                    (elemI(:,ei_geometryType) == custom) &
                    .or. &
                    (elemI(:,ei_geometryType) == force_main) &
                ) )

        if (npack > 0) then
            elemP(1:npack, ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == JB) &
                ) &
                .and. &
                ( &
                    (elemSI(:,esi_JunctionBranch_Exists) == 1) &
                ) &
                .and. &
                ( &
                    (elemI(:,ei_geometryType) == circular) &
                    .or. &
                    (elemI(:,ei_geometryType) == filled_circular) &
                    .or. &
                    (elemI(:,ei_geometryType) == rectangular_closed) &
                    .or. &
                    (elemI(:,ei_geometryType) == horiz_ellipse) &
                    .or. &
                    (elemI(:,ei_geometryType) == arch) &
                    .or. &
                    (elemI(:,ei_geometryType) == eggshaped) &
                    .or. &
                    (elemI(:,ei_geometryType) == horseshoe) &
                    .or. &
                    (elemI(:,ei_geometryType) == gothic) &
                    .or. &
                    (elemI(:,ei_geometryType) == catenary) &
                    .or. &
                    (elemI(:,ei_geometryType) == semi_elliptical) &
                    .or. &
                    (elemI(:,ei_geometryType) == basket_handle) &
                    .or. &
                    (elemI(:,ei_geometryType) == semi_circular) &
                    .or. &
                    (elemI(:,ei_geometryType) == custom) &
                    .or. &
                    (elemI(:,ei_geometryType) == force_main) &
                ) )
        end if
        if (setting%Debug%File%pack_mask_arrays) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine pack_nongeometry_static_elements
!
!==========================================================================
!==========================================================================
!
    subroutine pack_nongeometry_dynamic_elements()
        !--------------------------------------------------------------------------
        !% packed arrays for non geometry dynamic elements
        !--------------------------------------------------------------------------
        integer          :: ii
        integer, pointer :: ptype, npack, fup, fdn, eIDx(:)
        character(64) :: subroutine_name = 'pack_nongeometry_dynamic_elements'
        !--------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%pack_mask_arrays) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        eIdx => elemI(:,ei_Lidx)

        !% ep_AC ================================================
        !% - all elements that use AC
        ptype => col_elemP(ep_AC)
        npack => npack_elemP(ptype)
        npack = count( &
                (elemI(:,ei_tmType) == AC))
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemI(:,ei_tmType) == AC))
        end if

        !print *, 'CC_AC'
        !% ep_CC_AC ================================================
        !% - all channel conduit elements that use AC
        ptype => col_elemP(ep_CC_AC)
        npack => npack_elemP(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_tmType) == AC)     )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_tmType) == AC)     )
        end if

        !print *, 'CC_ETM'
        !% ep_CC_ETM ================================================
        !% - all channel conduit elements that use ETM
        ptype => col_elemP(ep_CC_ETM)
        npack => npack_elemP(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_tmType) == ETM)     )

        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_tmType) == ETM)     )
        end if

        !print *, 'CC_H_ETM'
        !% ep_CC_H_ETM ================================================
        !% - all channel conduit elements that have head time march using ETM
        ptype => col_elemP(ep_CC_H_ETM)
        npack => npack_elemP(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_HeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == ETM)     )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
               (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_HeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == ETM)     )
        end if

        !print *, 'CC_Q_AC'
        !% ep_CC_Q_AC ================================================
        !% - all channel conduit elements that have flow time march using AC
        ptype => col_elemP(ep_CC_Q_AC)
        npack => npack_elemP(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_QeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == AC)     )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
               (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_QeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == AC)     )
        end if

        !print *, 'CC_Q_ETM'
        !% ep_CC_Q_ETM ================================================
        !% - all channel conduit elements elements that have flow time march using ETM
        ptype => col_elemP(ep_CC_Q_ETM)
        npack => npack_elemP(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_QeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == ETM)     )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
               (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_QeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == ETM)     )
        end if

        !print *, 'CCJB_AC'
        !% ep_CCJB_AC ================================================
        !% - all channel conduit or junction branch elements elements that are AC
        ptype => col_elemP(ep_CCJB_AC)
        npack => npack_elemP(ptype)

        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or.&
                    (elemI(:,ei_elementType) == JB) &
                ) &
                .and. &
                (elemI(:,ei_tmType) == AC)     )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or.&
                    (elemI(:,ei_elementType) == JB) &
                ) &
                .and. &
                (elemI(:,ei_tmType) == AC)     )
        end if

        !print *, 'CC_AC_surcharged'
        !% ep_CC_AC_surcharged ================================================
        !% - all channel conduit elements elements that are AC and surcharged
        ptype => col_elemP(ep_CC_AC_surcharged)
        npack => npack_elemP(ptype)

        npack = count( &
                (  &
                    (elemI(:,ei_elementType) == CC) &
                ) &
                .and. &
                (elemI(:,ei_tmType) == AC)        &
                .and. &
                (elemYN(:,eYN_isSurcharged)))
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (  &
                    (elemI(:,ei_elementType) == CC) &
                ) &
                .and. &
                (elemI(:,ei_tmType) == AC)        &
                .and. &
                (elemYN(:,eYN_isSurcharged)))
        end if

        !print *, 'CCJB_AC_surcharged'
        !% ep_CCJB_AC_surcharged ================================================
        !% - all channel conduit or junction branch elements elements that are AC and surcharged
        ptype => col_elemP(ep_CCJB_AC_surcharged)
        npack => npack_elemP(ptype)

        npack = count( &
                (  &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB) &
                ) &
                .and. &
                (elemI(:,ei_tmType) == AC)        &
                .and. &
                (elemYN(:,eYN_isSurcharged)))
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (  &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB) &
                ) &
                .and. &
                (elemI(:,ei_tmType) == AC)        &
                .and. &
                (elemYN(:,eYN_isSurcharged)))
        end if

        !print *, 'CC_ALLtm_surcharged'
        !% ep_CC_ALLtm_surcharged ================================================
        !% - all channel conduit elements with any time march and surcharged
        ptype => col_elemP(ep_CC_ALLtm_surcharged)
        npack => npack_elemP(ptype)

        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                ) &
                .and. &
                ( &
                    (elemI(:,ei_tmType) == AC)        &
                    .or. &
                    (elemI(:,ei_tmType) == ETM)  &
                ) &
                .and. &
                (elemYN(:,eYN_isSurcharged)))
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                ) &
                .and. &
                ( &
                    (elemI(:,ei_tmType) == AC)        &
                    .or. &
                    (elemI(:,ei_tmType) == ETM)  &
                ) &
                .and. &
                (elemYN(:,eYN_isSurcharged)))
        end if

        ! print *, ' '
        ! print *, 'in ',trim(subroutine_name), ' ',ptype
        ! print *, 'thisP '
        ! print *, elemP(1:npack,ptype)
        ! print *, ' '
        ! print *, 'is_surcharged'
        ! print *, elemYN(elemP(1:npack,ptype),eYN_isSurcharged)
        ! print *, ' '

        !print *, 'CCJB_ALLtm_surcharged'
        !% ep_CCJB_ALLtm_surcharged ================================================
        !% - all channel conduit or junction branch elements with any time march and surcharged
        ptype => col_elemP(ep_CCJB_ALLtm_surcharged)
        npack => npack_elemP(ptype)

        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB) &
                ) &
                .and. &
                ( &
                    (elemI(:,ei_tmType) == AC)        &
                    .or. &
                    (elemI(:,ei_tmType) == ETM)  &
                ) &
                .and. &
                (elemYN(:,eYN_isSurcharged)))
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB) &
                ) &
                .and. &
                ( &
                    (elemI(:,ei_tmType) == AC)        &
                    .or. &
                    (elemI(:,ei_tmType) == ETM)  &
                ) &
                .and. &
                (elemYN(:,eYN_isSurcharged)))
        end if

        !% ep_CCJB_eETM_i_fAC ================================================
        !% conduits, channels, and junction branches that are ETM and have
        !% an adjacent face that is AC
        ptype => col_elemP(ep_CCJB_eETM_i_fAC)
        npack => npack_elemP(ptype)

        !% HACK: due to elements having u/s or d/s null faces, the code below
        !% is written in do loop untill other fix is figured out

        npack = zeroI
        do ii = 1,N_elem(this_image())
            fup => elemI(ii,ei_Mface_uL)
            fdn => elemI(ii,ei_Mface_dL)
            if ((fup /= nullvalueI) .and. (fdn /= nullValueI)) then
                if (((elemI(ii,ei_elementType) == CC) .or. (elemI(ii,ei_elementType) == JB)) &
                    .and. &
                    (elemI(ii,ei_tmType) == ETM) &
                    .and. &
                    ((faceYN(fup,fYN_isAC_adjacent)) .or. (faceYN(fdn,fYN_isAC_adjacent)))) then
                    
                    npack = npack + oneI
                    elemP(npack,ptype) = ii
                end if
            end if
        end do

        !% ep_CCJB_ETM ================================================
        !% - all channel conduit or junction branch that are ETM
        ptype => col_elemP(ep_CCJB_ETM)
        npack => npack_elemP(ptype)

        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB)  &
                 ) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB)  &
                 ) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )
        end if

        !print *, 'CC_ETM_surcharged'
        !% ep_CC_ETM_surcharged ================================================
        !% - all channel conduit or junction branch that are ETM and surcharged
        ptype => col_elemP(ep_CC_ETM_surcharged)
        npack => npack_elemP(ptype)

        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                 ) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                .and. &
                (elemYN(:,eYN_isSurcharged)) &
                )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                 ) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                .and. &
                (elemYN(:,eYN_isSurcharged)) &
                )
        end if

        !print *, 'CCJB_ETM_surcharged'
        !% ep_CCJB_ETM_surcharged ================================================
        !% - all channel conduit or junction branch that are ETM and surcharged
        ptype => col_elemP(ep_CCJB_ETM_surcharged)
        npack => npack_elemP(ptype)

        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB)  &
                 ) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                .and. &
                (elemYN(:,eYN_isSurcharged)) &
                )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB)  &
                 ) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                .and. &
                (elemYN(:,eYN_isSurcharged)) &
                )
        end if

        !print *, 'CCJM_H_AC_open'
        !% ep_CCJM_H_AC_open ================================================
        !% - all channel conduit or junction main elements solving head with AC and are non-surcharged
        ptype => col_elemP(ep_CCJM_H_AC_open)
        npack => npack_elemP(ptype)

        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM)  &
                 ) &
                .and. &
                (elemI(:,ei_HeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == AC) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM)  &
                 ) &
                .and. &
                (elemI(:,ei_HeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == AC) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                )
        end if

        !print *, 'CCJM_H_ETM'
        !% ep_CCJM_H_ETM ================================================
        !% - all channel conduit or junction main that use head solution with ETM
        ptype => col_elemP(ep_CCJM_H_ETM)
        npack => npack_elemP(ptype)

        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM)  &
                 ) &
                .and. &
                (elemI(:,ei_HeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM)  &
                 ) &
                .and. &
                (elemI(:,ei_HeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )
        end if

        !print *, 'ETM'
        !% ep_ETM ================================================
        !% - all elements that use ETM
        ptype => col_elemP(ep_ETM)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemI(:,ei_tmType) == ETM))
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx, &
                (elemI(:,ei_tmType) == ETM))
        end if

        !print *, 'JM_AC'
        !% ep_JM_AC ================================================
        !% - all elements that are junction mains and use AC
        ptype => col_elemP(ep_JM_AC)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemI(:,ei_elementType) == JM ) &
                .and. &
                (elemI(:,ei_tmType) == AC))
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemI(:,ei_elementType) == JM ) &
                .and. &
                (elemI(:,ei_tmType) == AC))
        end if

        !print *, 'JB_AC'
        !% ep_JB_AC ================================================
        !% - all elements that are junction mains and use AC
        ptype => col_elemP(ep_JB_AC)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemI(:,ei_elementType) == JB ) &
                .and. &
                (elemI(:,ei_tmType) == AC))
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemI(:,ei_elementType) == JB ) &
                .and. &
                (elemI(:,ei_tmType) == AC))
        end if

        !print *, 'JM_ETM'
        !% ep_JM_ETM ================================================
        !% - all elements that are junction mains and ETM
        ptype => col_elemP(ep_JM_ETM)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemI(:,ei_elementType) == JM ) &
                .and. &
                (elemI(:,ei_tmType) == ETM))

        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemI(:,ei_elementType) == JM ) &
                .and. &
                (elemI(:,ei_tmType) == ETM))
        end if

        !print *, 'JB_ETM'
        !% ep_JB_ETM ================================================
        !% - all elements that are junction mains and ETM
        ptype => col_elemP(ep_JB_ETM)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemI(:,ei_elementType) == JB ) &
                .and. &
                (elemI(:,ei_tmType) == ETM))

        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemI(:,ei_elementType) == JB ) &
                .and. &
                (elemI(:,ei_tmType) == ETM))
        end if

        !print *, 'NonSurcharged_AC'
        !% ep_NonSurcharged_AC ================================================
        !% - all AC elements that are not surcharged
        ptype => col_elemP(ep_NonSurcharged_AC)
        npack => npack_elemP(ptype)

        npack = count( &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == AC))
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == AC))
        end if

        !print *,'NonSurcharged_ALLtm'
        !% ep_NonSurcharged_ALLtm ================================================
        !% -- elements with any time march that are not surcharged
        ptype => col_elemP(ep_NonSurcharged_ALLtm)
        npack => npack_elemP(ptype)

        npack = count( &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_tmType) == AC) &
                    .or.&
                    (elemI(:,ei_tmType) == ETM) &
                ))

        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_tmType) == AC) &
                    .or.&
                    (elemI(:,ei_tmType) == ETM) &
                ))
        end if

        !print *, 'NonSurcharged_ETM'
        !% ep_NonSurcharged_ETM ================================================
        !% -- elements with ETM time march that are not surcharged 
        ptype => col_elemP(ep_NonSurcharged_ETM)
        npack => npack_elemP(ptype)

        npack = count( &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM))

        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM))
        end if

       

        !print *,'Surcharged_AC'
        !% ep_Surcharged_AC ================================================
        !% - all AC elements that are surcharged
        ptype => col_elemP(ep_Surcharged_AC)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == AC))

        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == AC))
        end if

        !print *, 'Surcharged_ALLtm'
        !% ep_Surcharged_ALLtm ================================================
        !% - all elements of any time march that are surcharged
        ptype => col_elemP(ep_Surcharged_ALLtm)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_tmType) == AC) &
                    .or. &
                    (elemI(:,ei_tmType) == ETM) &
                ))

        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_tmType) == AC) &
                    .or. &
                    (elemI(:,ei_tmType) == ETM) &
                ))
        end if

        !print *, 'Surcharged_ETM'
        !% ep_Surcharged_ETM ================================================
        !% - all ETM elements that are surcharged
        ptype => col_elemP(ep_Surcharged_ETM)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM))

        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM))
        end if

        !print *, 'CC closed elements'
        !% ep_CC_isclosed =====================================================
        !% --- all CC elements that are closed by elemR(:,er_Setting) = 0.0
        ptype => col_elemP(ep_CC_isclosed)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemI(:,ei_elementType) == CC ) &
                .and. &
                (elemR(:,er_Setting) == zeroR ) )

        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx, &
                (elemI(:,ei_elementType) == CC ) &
                .and. &
                (elemR(:,er_Setting) == zeroR ) )        
        end if

        if (setting%Debug%File%pack_mask_arrays) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine pack_nongeometry_dynamic_elements
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine pack_small_and_zero_depth_elements (whichTM)
        !%-----------------------------------------------------------------
        !% dynamic packing of elements with small and zero depths
        !%-----------------------------------------------------------------
        integer, intent(in)  :: whichTM
        integer, pointer :: ptype, npack, eIdx(:)
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------
        
        eIdx => elemI(:,ei_Lidx)

        select case (whichTM)
        case (ALLtm)
            !% ep_SmallDepth_CC_ALLtm ====================================
            !% - all Small depth that are CC and any time march
            ptype => col_elemP(ep_SmallDepth_CC_ALLtm)
            npack => npack_elemP(ptype)

            npack = count( &
                    (elemYN(:,eYN_isSmallDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == CC) &
                    .and. &
                    ( &
                        (elemI(:,ei_tmType) == ETM) &
                        .or. &
                        (elemI(:,ei_tmType) == AC) &
                    ) )

            if (npack > 0) then
                elemP(1:npack,ptype) = pack(eIdx,  &
                    (elemYN(:,eYN_isSmallDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == CC) &
                    .and. &
                    ( &
                        (elemI(:,ei_tmType) == ETM) &
                        .or. &
                        (elemI(:,ei_tmType) == AC) &
                    ) )
            end if
            !% ep_SmallDepth_JM_ALLtm ====================================
            !% - all Small depth that are JM and any time march
            ptype => col_elemP(ep_SmallDepth_JM_ALLtm)
            npack => npack_elemP(ptype)

            npack = count( &
                    (elemYN(:,eYN_isSmallDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == JM) &
                    .and. &
                    ( &
                        (elemI(:,ei_tmType) == ETM) &
                        .or. &
                        (elemI(:,ei_tmType) == AC) &
                    ) )

            if (npack > 0) then
                elemP(1:npack,ptype) = pack(eIdx,  &
                    (elemYN(:,eYN_isSmallDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == JM) &
                    .and. &
                    ( &
                        (elemI(:,ei_tmType) == ETM) &
                        .or. &
                        (elemI(:,ei_tmType) == AC) &
                    ) )
            end if
            
            !% ep_ZeroDepth_CC_ALLtm ====================================
            !% - all zero depth that are CC and any time march
            ptype => col_elemP(ep_ZeroDepth_CC_ALLtm)
            npack => npack_elemP(ptype)

            npack = count( &
                    (elemYN(:,eYN_isZeroDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == CC) &
                    .and. &
                    ( &
                        (elemI(:,ei_tmType) == ETM) &
                        .or. &
                        (elemI(:,ei_tmType) == AC) &
                    ) )

            if (npack > 0) then
                elemP(1:npack,ptype) = pack(eIdx,  &
                    (elemYN(:,eYN_isZeroDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == CC) &
                    .and. &
                    ( &
                        (elemI(:,ei_tmType) == ETM) &
                        .or. &
                        (elemI(:,ei_tmType) == AC) &
                    ) )
            end if

            !% ep_ZeroDepth_JM_ALLtm ====================================
            !% - all Zero depth that are JM for any TM
            ptype => col_elemP(ep_ZeroDepth_JM_ALLtm)
            npack => npack_elemP(ptype)

            npack = count( &
                    (elemYN(:,eYN_isZeroDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == JM) &
                    .and. &
                    ( &
                        (elemI(:,ei_tmType) == ETM) &
                        .or. &
                        (elemI(:,ei_tmType) == AC) &
                    ) )

            if (npack > 0) then
                elemP(1:npack,ptype) = pack(eIdx,  &
                    (elemYN(:,eYN_isZeroDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == JM) &
                    .and. &
                    ( &
                        (elemI(:,ei_tmType) == ETM) &
                        .or. &
                        (elemI(:,ei_tmType) == AC) &
                    ) )
            end if

        case (ETM)
            !% ep_SmallDepth_CC_ETM ====================================
            !% - all Small depth that are CC and ETM time march
            ptype => col_elemP(ep_SmallDepth_CC_ETM)
            npack => npack_elemP(ptype)

            npack = count( &
                    (elemYN(:,eYN_isSmallDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == CC) &
                    .and. &
                    (elemI(:,ei_tmType) == ETM) )

            if (npack > 0) then
                elemP(1:npack,ptype) = pack(eIdx,  &
                    (elemYN(:,eYN_isSmallDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == CC) &
                    .and. &
                    (elemI(:,ei_tmType) == ETM))
            end if

             !% BeginNew 20220122brh
            !% ep_SmallDepth_JM_ETM ====================================
            !% - all Small depth that are JM and ETM time march
            ptype => col_elemP(ep_SmallDepth_JM_ETM)
            npack => npack_elemP(ptype)

            npack = count( &
                    (elemYN(:,eYN_isSmallDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == JM) &
                    .and. &
                    (elemI(:,ei_tmType) == ETM) )

            if (npack > 0) then
                elemP(1:npack,ptype) = pack(eIdx,  &
                    (elemYN(:,eYN_isSmallDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == JM) &
                    .and. &
                    (elemI(:,ei_tmType) == ETM))
            end if
            !% EndNew 20220122brh

            !% ep_ZeroDepth_CC_ETM ====================================
            !% - all zero depth that are CC and ETM time march
            ptype => col_elemP(ep_ZeroDepth_CC_ETM)
            npack => npack_elemP(ptype)

            npack = count( &
                    (elemYN(:,eYN_isZeroDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == CC) &
                    .and. &
                    (elemI(:,ei_tmType) == ETM) )

            if (npack > 0) then
                elemP(1:npack,ptype) = pack(eIdx,  &
                    (elemYN(:,eYN_isZeroDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == CC) &
                    .and. &
                    (elemI(:,ei_tmType) == ETM) )
            end if

            !% ep_ZeroDepth_JM_ETM ====================================
            !% - all Zero depth that are JM for ETM
            ptype => col_elemP(ep_ZeroDepth_JM_ETM)
            npack => npack_elemP(ptype)

            npack = count( &
                    (elemYN(:,eYN_isZeroDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == JM) &
                    .and. &
                    (elemI(:,ei_tmType) == ETM))

            if (npack > 0) then
                elemP(1:npack,ptype) = pack(eIdx,  &
                    (elemYN(:,eYN_isZeroDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == JM) &
                    .and. &
                    (elemI(:,ei_tmType) == ETM))
            end if

        case (AC)
                !% BeginNew 20220122brh
            !% ep_SmallDepth_CC_AC ====================================
            !% - all Small depth that are CC and AC time march
            ptype => col_elemP(ep_SmallDepth_CC_AC)
            npack => npack_elemP(ptype)

            npack = count( &
                    (elemYN(:,eYN_isSmallDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == CC) &
                    .and. &
                    (elemI(:,ei_tmType) == AC) )

            if (npack > 0) then
                elemP(1:npack,ptype) = pack(eIdx,  &
                    (elemYN(:,eYN_isSmallDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == CC) &
                    .and. &
                    (elemI(:,ei_tmType) == AC))
            end if

            !% ep_SmallDepth_JM_AC ====================================
            !% - all Small depth that are JM and AC time march
            ptype => col_elemP(ep_SmallDepth_JM_AC)
            npack => npack_elemP(ptype)

            npack = count( &
                    (elemYN(:,eYN_isSmallDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == JM) &
                    .and. &
                    (elemI(:,ei_tmType) == AC) )

            if (npack > 0) then
                elemP(1:npack,ptype) = pack(eIdx,  &
                    (elemYN(:,eYN_isSmallDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == JM) &
                    .and. &
                    (elemI(:,ei_tmType) == AC))
            end if
            
            !% ep_ZeroDepth_CC_AC ====================================
            !% - all zero depth that are CC and AC time march
            ptype => col_elemP(ep_ZeroDepth_CC_AC)
            npack => npack_elemP(ptype)

            npack = count( &
                    (elemYN(:,eYN_isZeroDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == CC) &
                    .and. &
                    (elemI(:,ei_tmType) == AC) )

            if (npack > 0) then
                elemP(1:npack,ptype) = pack(eIdx,  &
                    (elemYN(:,eYN_isZeroDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == CC) &
                    .and. &
                    (elemI(:,ei_tmType) == AC) )
            end if

            !% ep_ZeroDepth_JM_AC ====================================
            !% - all Zero depth that are JM for AC
            ptype => col_elemP(ep_ZeroDepth_JM_AC)
            npack => npack_elemP(ptype)

            npack = count( &
                    (elemYN(:,eYN_isZeroDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == JM) &
                    .and. &
                    (elemI(:,ei_tmType) == AC))

            if (npack > 0) then
                elemP(1:npack,ptype) = pack(eIdx,  &
                    (elemYN(:,eYN_isZeroDepth)) &
                    .and. &
                    (elemI(:,ei_elementType) == JM) &
                    .and. &
                    (elemI(:,ei_tmType) == AC))
            end if

            !% EndNew 20220122brh
        case default
            print *, 'CODE ERROR: time march type unknown for # ', whichTM
            print *, 'which has key ',trim(reverseKey(whichTM))
            !stop 
            call util_crashpoint(3987053)
            !return
        end select
        

        !     % THE THE POS/NEG SLOPE SHOULD BE OBSOLETE
        ! !% ep_SmallDepth_CC_ALLtm_posSlope
        !     !% - all small depth that are CC and any time march with positive bottom slope
        !     ptype => col_elemP(ep_SmallDepth_CC_ALLtm_posSlope)
        !     npack => npack_elemP(ptype)

        !     npack = count( &
        !             (elemYN(:,eYN_isSmallDepth)) &
        !             .and. &
        !             (elemI(:,ei_elementType) == CC) &
        !             .and. &
        !             (elemR(:,er_BottomSlope) .ge. zeroR) &
        !             .and. &
        !             ( &
        !                 (elemI(:,ei_tmType) == ETM) &
        !                 .or. &
        !                 (elemI(:,ei_tmType) == AC) &
        !             ) )

        !     if (npack > 0) then
        !         elemP(1:npack,ptype) = pack(eIdx,  &
        !             (elemYN(:,eYN_isSmallDepth)) &
        !             .and. &
        !             (elemI(:,ei_elementType) == CC) &
        !             .and. &
        !             (elemR(:,er_BottomSlope) .ge. zeroR) &
        !             .and. &
        !             ( &
        !                 (elemI(:,ei_tmType) == ETM) &
        !                 .or. &
        !                 (elemI(:,ei_tmType) == AC) &
        !             ) )
        !     end if

        !     !% ep_SmallDepth_CC_ALLtm_negSlope
        !     !% - all small depth that are CC and any time march with negative (adverse) bottom slope
        !     ptype => col_elemP(ep_SmallDepth_CC_ALLtm_negSlope)
        !     npack => npack_elemP(ptype)

        !     npack = count( &
        !             (elemYN(:,eYN_isSmallDepth)) &
        !             .and. &
        !             (elemI(:,ei_elementType) == CC) &
        !             .and. &
        !             (elemR(:,er_BottomSlope) < zeroR) &
        !             .and. &
        !             ( &
        !                 (elemI(:,ei_tmType) == ETM) &
        !                 .or. &
        !                 (elemI(:,ei_tmType) == AC) &
        !             ) )

        !     if (npack > 0) then
        !         elemP(1:npack,ptype) = pack(eIdx,  &
        !             (elemYN(:,eYN_isSmallDepth)) &
        !             .and. &
        !             (elemI(:,ei_elementType) == CC) &
        !             .and. &
        !             (elemR(:,er_BottomSlope) < zeroR) &
        !             .and. &
        !             ( &
        !                 (elemI(:,ei_tmType) == ETM) &
        !                 .or. &
        !                 (elemI(:,ei_tmType) == AC) &
        !             ) )
        !     end if



         
        !print *, 'CC_NOTsmalldepth'
        !% ep_CC_NOTsmalldepth  ====================================
        !% Flow solution that are NOT small volume or zero depth
        !% -- needed to limit where CFL is computed and volume conservation
        ptype => col_elemP(ep_CC_NOTsmalldepth)
        npack => npack_elemP(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_QeqType) == time_march) &
                .and. &
                (.not. elemYN(:,eYN_isSmallDepth)) &
                .and. &
                (.not. elemYN(:,eYN_isZeroDepth))     )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_QeqType) == time_march) &
                .and. &
                (.not. elemYN(:,eYN_isSmallDepth))   &
                .and. &
                (.not. elemYN(:,eYN_isZeroDepth))     )
        end if

        !print *, 'JBJM_NOTsmalldepth'
        !% ep_JBJM_NOTsmalldepth  ====================================
        !% Flow solution that are NOT small volume or zero depth
        !% -- needed to limit where CFL is computed
        ptype => col_elemP(ep_JBJM_NOTsmalldepth)
        npack => npack_elemP(ptype)
        npack = count( &
                (      &
                      ( elemI(:,ei_elementType) == JM) &
                 .or. ((elemI(:,ei_elementType) == JB) .and. (elemSI(:,esi_JunctionBranch_Exists) == oneR) ) &   
                ) &
                .and. &
                (.not. elemYN(:,eYN_isSmallDepth)) &
                .and. &
                (.not. elemYN(:,eYN_isZeroDepth))     )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (      &
                      ( elemI(:,ei_elementType) == JM) &
                 .or. ((elemI(:,ei_elementType) == JB) .and. (elemSI(:,esi_JunctionBranch_Exists) == oneR) ) &   
                ) &
                .and. &
                (.not. elemYN(:,eYN_isSmallDepth)) &
                .and. &
                (.not. elemYN(:,eYN_isZeroDepth))     ) 
        end if

        !print *, 'CCJBJM_NOTsmalldepth'
        !% ep_CCJBJM_NOTsmalldepth  ====================================
        !% Flow solution that are NOT small volume or zero depth
        !% -- needed to limit where CFL is computed
        ptype => col_elemP(ep_CCJBJM_NOTsmalldepth)
        npack => npack_elemP(ptype)
        npack = count( &
                (      &
                      ( elemI(:,ei_elementType) == CC) &
                 .or. ( elemI(:,ei_elementType) == JM) &
                 .or. ((elemI(:,ei_elementType) == JB) .and. (elemSI(:,esi_JunctionBranch_Exists) == oneR) ) &   
                ) &
                .and. &
                (.not. elemYN(:,eYN_isSmallDepth)) &
                .and. &
                (.not. elemYN(:,eYN_isZeroDepth))     )
        ! print *, ' '
        ! print *, ' in pack_mask_arrays checking CCJBJM_NOTsmalldepth'
        ! print *, 'npack ',npack        
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (      &
                      ( elemI(:,ei_elementType) == CC) &
                 .or. ( elemI(:,ei_elementType) == JM) &
                 .or. ((elemI(:,ei_elementType) == JB) .and. (elemSI(:,esi_JunctionBranch_Exists) == oneR) ) &   
                ) &
                .and. &
                (.not. elemYN(:,eYN_isSmallDepth)) &
                .and. &
                (.not. elemYN(:,eYN_isZeroDepth))     ) 
        end if

        !print *, 'CCJM_NOTsmalldepth'
        !% ep_CCJM_NOTsmalldepth  ====================================
        !% Flow solution that are NOT small volume or zero depth
        !% -- alternate needed to limit where CFL is computed
        ptype => col_elemP(ep_CCJM_NOTsmalldepth)
        npack => npack_elemP(ptype)
        npack = count( &
                (      &
                      ( elemI(:,ei_elementType) == CC) &
                 .or. ( elemI(:,ei_elementType) == JM) &   
                ) &
                .and. &
                (.not. elemYN(:,eYN_isSmallDepth)) &
                .and. &
                (.not. elemYN(:,eYN_isZeroDepth))     )
        ! print *, ' '
        ! print *, ' in pack_mask_arrays checking CCJM_NOTsmalldepth'
        ! print *, 'npack ',npack        
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (      &
                      ( elemI(:,ei_elementType) == CC) &
                 .or. ( elemI(:,ei_elementType) == JM) &   
                ) &
                .and. &
                (.not. elemYN(:,eYN_isSmallDepth)) &
                .and. &
                (.not. elemYN(:,eYN_isZeroDepth))     ) 
        end if


    end subroutine pack_small_and_zero_depth_elements
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine pack_static_interior_faces()
        !--------------------------------------------------------------------------
        !% packed arrays for static faces
        !--------------------------------------------------------------------------
        integer :: ii, image
        integer, pointer :: Nfaces, ptype, npack, fIdx(:), eup(:), edn(:)

        character(64) :: subroutine_name = 'pack_static_interior_faces'

        !--------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%pack_mask_arrays) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% pointing to the number of faces in this image
        image  = this_image()
        Nfaces => N_face(image)

        fIdx => faceI(1:Nfaces,fi_Lidx)
        eup  => faceI(1:Nfaces,fi_Melem_uL)
        edn  => faceI(1:Nfaces,fi_Melem_dL)

        !% fp_all
        !% - all faces execpt boundary, null, and shared faces
        ptype => col_faceP(fp_all)
        npack => npack_faceP(ptype)

        npack = count(faceYN(1:Nfaces,fYN_isInteriorFace))

        if (npack > 0) then
            faceP(1:npack,ptype) = pack(fIdx, &
                faceYN(1:Nfaces,fYN_isInteriorFace) &
                )
        end if

        !% fp_J1
        !% - faces with only one link that are not inflow BC
        ptype => col_faceP(fp_J1)
        npack => npack_faceP(ptype)

        npack = count(faceI(1:Nfaces,fi_BCtype)==BCnone)

        if (npack > 0) then 
            faceP(1:npack, ptype) = pack( fIdx, &
                    faceI(1:Nfaces,fi_BCtype)==BCnone)
        end if

        !% fp_J1_BCup
        !% - faces with only one link that are not outfalls
        !%   these are either J1 or BCup
        ptype => col_faceP(fp_J1_BCup)
        npack => npack_faceP(ptype)

        npack = count(                              &
                (faceI(1:Nfaces,fi_BCtype)==BCnone) &
                .or.                                &
                (faceI(1:Nfaces,fi_BCtype)==BCup) )

        if (npack > 0) then 
            faceP(1:npack, ptype) = pack( fIdx,    &
                (faceI(1:Nfaces,fi_BCtype)==BCnone) &
                .or.                                &
                (faceI(1:Nfaces,fi_BCtype)==BCup) )
        end if

        !% fp_Diag
        !% - all faces adjacent to a diagnostic element
        ptype => col_faceP(fp_Diag)
        npack => npack_faceP(ptype)

        npack =  count( &
                faceYN(1:Nfaces,fYN_isInteriorFace)   &
                .and. &
                (elemI(edn,ei_HeqType) == diagnostic) &
                .or.  &
                (elemI(edn,ei_QeqType) == diagnostic) &
                .or.  &
                (elemI(eup,ei_HeqType) == diagnostic) &
                .or.  &
                (elemI(eup,ei_QeqType) == diagnostic))

        if (npack > 0) then
            faceP(1:npack, ptype) = pack( fIdx, &
                    faceYN(1:Nfaces,fYN_isInteriorFace)   &
                    .and. &
                    (elemI(edn,ei_HeqType) == diagnostic) &
                    .or.  &
                    (elemI(edn,ei_QeqType) == diagnostic) &
                    .or.  &
                    (elemI(eup,ei_HeqType) == diagnostic) &
                    .or.  &
                    (elemI(eup,ei_QeqType) == diagnostic))
        end if

        if (setting%Debug%File%pack_mask_arrays) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine pack_static_interior_faces
!
!==========================================================================
!==========================================================================
!
    subroutine pack_dynamic_interior_faces()
        !--------------------------------------------------------------------------
        !% packed arrays for dynamic faces
        !%
        !% HACK: Should the jump packing be called after all jump conditions are
        !% changed? or can it wait until the end of a time step? Note that this
        !% simply packs what is stored in faceI(:,fi_jump_type) as the actual
        !% computation of what is a jump is in the identify_hydraulic_jump subroutine.
        !--------------------------------------------------------------------------
        integer          :: ii, image
        integer, pointer :: Nfaces, ptype, npack, fIdx(:), eup(:), edn(:)
        character(64) :: subroutine_name = 'pack_dynamic_interior_faces'
        !--------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%pack_mask_arrays) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% pointing to the number of faces in this image
        image  = this_image()
        Nfaces => N_face(image)

        fIdx => faceI(1:Nfaces,fi_Lidx)
        eup  => faceI(1:Nfaces,fi_Melem_uL)
        edn  => faceI(1:Nfaces,fi_Melem_dL)

        !% fp_AC
        !% - faces with any AC adjacent
        ptype => col_faceP(fp_AC)
        npack => npack_faceP(ptype)

        npack = count( &
                faceYN(1:Nfaces,fYN_isInteriorFace) &
                .and.&
                (elemI(edn,ei_tmType) == AC) &
                .or. &
                (elemI(eup,ei_tmType) == AC))

        if (npack > 0) then
            faceP(1:npack, ptype) = pack( fIdx, &
                faceYN(1:Nfaces,fYN_isInteriorFace) &
                .and.&
                (elemI(edn,ei_tmType) == AC)    &
                .or. &
                (elemI(eup,ei_tmType) == AC)    )
        end if

        !% fp_JumpUp
        !% -Hydraulic jump from nominal upstream to downstream
        ptype => col_faceP(fp_JumpUp)
        npack => npack_faceP(ptype)

        npack = count( &
            faceYN(1:Nfaces,fYN_isInteriorFace) &
            .and. &
            faceI(1:Nfaces,fi_jump_type) == jump_from_upstream )

        if (npack > 0) then
            faceP(1:npack, ptype) = pack( fIdx, &
                faceYN(1:Nfaces,fYN_isInteriorFace) &
                .and.&
                faceI(1:Nfaces,fi_jump_type) == jump_from_upstream )
        end if

        !% fp_JumpDn
        !% --Hydraulic jump from nominal downstream to upstream
        ptype => col_faceP(fp_JumpDn)
        npack => npack_faceP(ptype)

        npack = count( &
            faceYN(1:Nfaces,fYN_isInteriorFace) &
            .and. &
            faceI(1:Nfaces,fi_jump_type) == jump_from_downstream )

        if (npack > 0) then
            faceP(1:npack, ptype) = pack( fIdx, &
                faceYN(1:Nfaces,fYN_isInteriorFace) &
                .and. &
                faceI(1:Nfaces,fi_jump_type) == jump_from_downstream )
        end if

        if (setting%Debug%File%pack_mask_arrays) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine pack_dynamic_interior_faces
!
!==========================================================================
!==========================================================================
!
    subroutine pack_static_shared_faces()
        !--------------------------------------------------------------------------
        !% packed arrays for static shared faces
        !--------------------------------------------------------------------------
        integer :: ii
        integer, pointer :: ptype, npack, fIdx(:), eup, edn, gup, gdn, Nfaces
        integer, pointer :: c_image, N_shared_faces, thisP
        logical, pointer :: isUpGhost, isDnGhost
        integer(kind=8) :: crate, cmax, cval
        character(64) :: subroutine_name = 'pack_static_shared_faces'
        !--------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%pack_mask_arrays) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        ! !% start the shared timer    
        ! sync all    
        ! if (this_image()==1) then
        !     call system_clock(count=cval,count_rate=crate,count_max=cmax)
        !     setting%Time%WallClock%SharedStart = cval
        ! end if

        !% pointing to the number of faces in this image
        Nfaces => N_face(this_image())

        fIdx   => faceI(1:Nfaces,fi_Lidx)

        ! % fp_all (shared faces)
        !% - all faces that are shared across images (Internal Boundary Faces)
        ptype => col_facePS(fp_all)
        npack => npack_facePS(ptype)

        npack = count( &
                faceYN(1:Nfaces,fYN_isSharedFace))

        if (npack > 0) then
            facePS(1:npack, ptype) = pack( fIdx, &
                faceYN(1:Nfaces,fYN_isSharedFace))
        end if

        sync all

        !% fp_Diag (shared faces)
        !% - all faces adjacent to a diagnostic element which is shared across images
        ptype => col_facePS(fp_Diag)
        npack => npack_facePS(ptype)

        !% pointer towards the total number of shared faces in an image
        N_shared_faces  => npack_facePS(fp_all)

        if (N_shared_faces > 0) then
            do ii = 1,N_shared_faces
                thisP       => facePS(ii,fp_all)
                eup         => faceI(thisP,fi_Melem_uL)
                edn         => faceI(thisP,fi_Melem_dL)
                isUpGhost   => faceYN(thisP,fYN_isUpGhost)
                gup         => faceI(thisP,fi_GhostElem_uL)
                isDnGhost   => faceYN(thisP,fYN_isDnGhost)
                gdn         => faceI(thisP,fi_GhostElem_dL)
                c_image     => faceI(thisP,fi_Connected_image)

                if (isUpGhost) then
                   if ((elemI(gup,ei_HeqType)[c_image] == diagnostic) .or.  &
                       (elemI(gup,ei_QeqType)[c_image] == diagnostic) .or.  &
                       (elemI(edn,ei_HeqType)          == diagnostic) .or.  &
                       (elemI(edn,ei_QeqType)          == diagnostic))   then

                        !% advance the number of pack value
                        npack = npack + oneI
                        !% save the face index
                        facePS(npack,ptype) = thisP
                    end if

                elseif (isDnGhost) then
                    if ((elemI(gdn,ei_HeqType)[c_image] == diagnostic) .or.  &
                        (elemI(gdn,ei_QeqType)[c_image] == diagnostic) .or.  &
                        (elemI(eup,ei_HeqType)          == diagnostic) .or.  &
                        (elemI(eup,ei_QeqType)          == diagnostic))  then

                        !% advance the number of pack value
                        npack = npack + oneI
                        !% save the face index
                        facePS(npack,ptype) = thisP
                    end if
                end if
            end do
        end if
        faceYN(faceP(1:npack, ptype),fYN_isDiag_adjacent) = .true.
        !% stop the shared timer
        sync all
        ! if (this_image()==1) then
        !     call system_clock(count=cval,count_rate=crate,count_max=cmax)
        !     setting%Time%WallClock%SharedStop = cval
        !     setting%Time%WallClock%SharedCumulative &
        !             = setting%Time%WallClock%SharedCumulative &
        !             + setting%Time%WallClock%SharedStop &
        !             - setting%Time%WallClock%SharedStart
        ! end if 
        if (setting%Debug%File%pack_mask_arrays) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine pack_static_shared_faces
!
!==========================================================================
!==========================================================================
!
    subroutine pack_dynamic_shared_faces()
        !--------------------------------------------------------------------------
        !% packed arrays for dynamic shared faces
        !%
        !% HACK: Should the jump packing be called after all jump conditions are
        !% changed? or can it wait until the end of a time step? Note that this
        !% simply packs what is stored in faceI(:,fi_jump_type) as the actual
        !% computation of what is a jump is in the identify_hydraulic_jump subroutine.
        !--------------------------------------------------------------------------
        integer          :: ii, image
        integer, pointer :: ptype, npack, fIdx(:), Nfaces
        integer, pointer :: N_shared_faces, thisP, eup, edn, gup, gdn, c_image
        logical, pointer :: isUpGhost, isDnGhost
        integer(kind=8) :: crate, cmax, cval
        character(64)    :: subroutine_name = 'pack_dynamic_shared_faces'
        !--------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%pack_mask_arrays) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            
        !% start the shared timer    
        sync all    
        if (this_image()==1) then
            call system_clock(count=cval,count_rate=crate,count_max=cmax)
            setting%Time%WallClock%SharedStart = cval
            setting%Time%WallClock%SharedStart_B = cval
        end if

        !% pointing to the number of faces in this image
        image  = this_image()
        Nfaces => N_face(image)

        fIdx => faceI(1:Nfaces,fi_Lidx)

        !% fp_AC (shared faces)
        !% - faces with any AC adjacent which is shared across images
        ptype => col_facePS(fp_AC)
        npack => npack_facePS(ptype)

        !% pointer towards the total number of shared faces in an image
        N_shared_faces  => npack_facePS(fp_all)

        if (N_shared_faces > 0) then
            do ii = 1,N_shared_faces
                thisP       => facePS(ii,fp_all)
                eup         => faceI(thisP,fi_Melem_uL)
                edn         => faceI(thisP,fi_Melem_dL)
                isUpGhost   => faceYN(thisP,fYN_isUpGhost)
                gup         => faceI(thisP,fi_GhostElem_uL)
                isDnGhost   => faceYN(thisP,fYN_isDnGhost)
                gdn         => faceI(thisP,fi_GhostElem_dL)
                c_image     => faceI(thisP,fi_Connected_image)

                if (isUpGhost) then
                    if ((elemI(gup,ei_tmType)[c_image] == AC) .or.  &
                        (elemI(edn,ei_tmType)          == AC))  then

                        !% advance the number of pack value
                        npack = npack + oneI
                        !% save the face index
                        facePS(npack,ptype) = thisP
                    end if

                elseif (isDnGhost) then
                    if ((elemI(gdn,ei_tmType)[c_image] == AC) .or.  &
                        (elemI(eup,ei_tmType)          == AC))  then

                        !% advance the number of pack value
                        npack = npack + oneI
                        !% save the face index
                        facePS(npack,ptype) = thisP
                    end if
                end if
            end do
        end if

        !% fp_JumpUp (shared faces)
        !% -Hydraulic jump from nominal upstream to downstream
        ptype => col_facePS(fp_JumpUp)
        npack => npack_facePS(ptype)

        npack = count( &
                faceYN(1:Nfaces,fYN_isSharedFace)              .and. &
                (faceI(1:Nfaces,fi_jump_type) == jump_from_upstream))

        if (npack > 0) then
            facePS(1:npack, ptype) = pack( fIdx, &
                faceYN(1:Nfaces,fYN_isSharedFace)              .and. &
                (faceI(1:Nfaces,fi_jump_type) == jump_from_upstream))
        end if

        !% fp_JumpDn
        !% --Hydraulic jump from nominal downstream to upstream
        ptype => col_facePS(fp_JumpDn)
        npack => npack_facePS(ptype)

        npack = count( &
                faceYN(1:Nfaces,fYN_isSharedFace)                .and. &
                (faceI(1:Nfaces,fi_jump_type) == jump_from_downstream))

        if (npack > 0) then
            facePS(1:npack, ptype) = pack( fIdx, &
                faceYN(1:Nfaces,fYN_isSharedFace)                .and. &
                (faceI(1:Nfaces,fi_jump_type) == jump_from_downstream))
        end if

        !% stop the shared timer
        sync all
        if (this_image()==1) then
            call system_clock(count=cval,count_rate=crate,count_max=cmax)
            setting%Time%WallClock%SharedStop = cval
            setting%Time%WallClock%SharedCumulative &
                    = setting%Time%WallClock%SharedCumulative &
                    + setting%Time%WallClock%SharedStop &
                    - setting%Time%WallClock%SharedStart

            setting%Time%WallClock%SharedStop_B = cval
            setting%Time%WallClock%SharedCumulative_B &
                    = setting%Time%WallClock%SharedCumulative_B &
                    + setting%Time%WallClock%SharedStop_B &
                    - setting%Time%WallClock%SharedStart_B            
        end if 
        if (setting%Debug%File%pack_mask_arrays) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine pack_dynamic_shared_faces
!%
!%==========================================================================
!%==========================================================================
!%
!% OBSOLETE?
    ! subroutine pack_link_node_output
    !     integer :: ii, jj, link_output_idx_length, node_output_idx_length
    !     character(64)    :: subroutine_name = 'pack_link_node_output'
    !     !% --------------------------------------------------------------------------
    !     !if (crashYN) return
    !     if (setting%Debug%File%pack_mask_arrays) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !     !% count the amount of valid output links
    !     link_output_idx_length = count(link_output_idx(:) /= nullvalueI)
    !     node_output_idx_length = count(node_output_idx(:) /= nullvalueI)

    !     !% allocate the pack
    !     allocate(link%P%have_output(link_output_idx_length))
    !     allocate(node%P%have_output(node_output_idx_length))

    !     !% fill the pack
    !     link%P%have_output = pack(link%I(link_output_idx(1:link_output_idx_length), li_idx), &
    !         link%I(link_output_idx(1:link_output_idx_length), li_P_image) == this_image())
    !     node%P%have_output = pack(node%I(node_output_idx(1:node_output_idx_length), ni_idx), &
    !         node%I(node_output_idx(1:node_output_idx_length), ni_P_image) == this_image())

    !     if (setting%Debug%File%pack_mask_arrays) &
    !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine pack_link_node_output
!
!==========================================================================
! END MODULE
!==========================================================================
!    
end module pack_mask_arrays
