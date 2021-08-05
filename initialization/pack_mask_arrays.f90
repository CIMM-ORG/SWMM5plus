!
!% module pack_mask_arrays
!
!
!==========================================================================
!
!==========================================================================
!
module pack_mask_arrays

    use define_indexes
    use define_keys
    use define_globals
    use define_settings

    implicit none

    private

    public :: pack_mask_arrays_all
    public :: pack_dynamic_arrays
    public :: pack_nodes
    public :: pack_bc

contains
    !
    !==========================================================================
    ! PUBLIC
    !==========================================================================
    !
    subroutine pack_mask_arrays_all()
        !--------------------------------------------------------------------------
        !
        !% set all the static packs and masks
        !
        !--------------------------------------------------------------------------

        integer :: ii

        character(64) :: subroutine_name = 'pack_mask_arrays_all'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

        call mask_faces_whole_array_static()
        call pack_geometry_alltm_elements()
        call pack_geometry_etm_elements()
        call pack_geometry_ac_elements()
        call pack_nongeometry_static_elements()
        call pack_nongeometry_dynamic_elements()
        call pack_static_interior_faces()
        call pack_static_shared_faces()
        call pack_dynamic_interior_faces()
        call pack_dynamic_shared_faces()

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
                   call execute_command_line('')
                enddo

            endif
        endif

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_mask_arrays_all
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine pack_dynamic_arrays()
        !--------------------------------------------------------------------------
        !
        !% set all the dynamic packs and masks
        !
        !--------------------------------------------------------------------------

        character(64) :: subroutine_name = 'pack_dynamic_arrays'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

        call pack_geometry_etm_elements()
        call pack_geometry_ac_elements()
        call pack_nongeometry_dynamic_elements()
        call pack_dynamic_interior_faces()
        call pack_dynamic_shared_faces()

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_dynamic_arrays
    !
    !==========================================================================
    ! PRIVATE
    !==========================================================================
    !
    subroutine mask_faces_whole_array_static()
        !--------------------------------------------------------------------------
        !
        !% find all the faces except boundary and null faces
        !
        !--------------------------------------------------------------------------

        integer, pointer :: mcol

        character(64) :: subroutine_name = 'mask_faces_whole_array_static'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

        mcol => col_faceM(fm_all)

        faceM(:,mcol) = ( &
            (faceI(:,fi_BCtype) == doesnotexist) &
            .and. &
            (.not. faceYN(:,fYN_isnull))  &
            .and. &
            (.not. faceYN(:,fYN_isSharedFace)) &
            )

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine mask_faces_whole_array_static
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine pack_geometry_alltm_elements()
        !--------------------------------------------------------------------------
        !
        !% packed arrays for geometry types in elemPGalltm
        !
        !--------------------------------------------------------------------------

        integer, pointer :: ptype, npack, eIDx(:)

        character(64) :: subroutine_name = 'pack_geometry_alltm_elements'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

        eIdx => elemI(:,ei_Lidx)

        !% rectangular channels, conduits and junction main
        ptype => col_elemPGalltm(epg_CCJM_rectangular_nonsurcharged)
        npack => npack_elemPGalltm(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemI(:,ei_geometryType) == rectangular) &
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
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemI(:,ei_geometryType) == rectangular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ))
        endif

        !% trapezoidal channels, conduits and junction main
        ptype => col_elemPGalltm(epg_CCJM_trapezoidal_nonsurcharged)
        npack => npack_elemPGalltm(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM) &
                ) &
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
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM) &
                ) &
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
        endif

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_geometry_alltm_elements
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine pack_geometry_ac_elements()
        !--------------------------------------------------------------------------
        !
        !% packed arrays for geometry types for AC elements
        !
        !--------------------------------------------------------------------------

        integer, pointer :: ptype, npack, eIDx(:)

        character(64) :: subroutine_name = 'pack_geometry_alltm_elements'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

        eIdx => elemI(:,ei_Lidx)

        !% rectangular channels, conduits and junction main
        ptype => col_elemPGac(epg_CCJM_rectangular_nonsurcharged)
        npack => npack_elemPGac(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemI(:,ei_geometryType) == rectangular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )

        if (npack > 0) then
            elemPGac(1:npack, ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemI(:,ei_geometryType) == rectangular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged))&
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )
        endif

        !% trapezoidal channels, conduits and junction main
        ptype => col_elemPGac(epg_CCJM_trapezoidal_nonsurcharged)
        npack => npack_elemPGac(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemI(:,ei_geometryType) == trapezoidal) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )

        if (npack > 0) then
            elemPGac(1:npack, ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemI(:,ei_geometryType) == trapezoidal) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged))&
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )
        endif

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_geometry_ac_elements
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine pack_geometry_etm_elements()
        !--------------------------------------------------------------------------
        !
        !% packed arrays for geometry types
        !
        !--------------------------------------------------------------------------

        integer, pointer :: ptype, npack, eIDx(:)

        character(64) :: subroutine_name = 'pack_geometry_etm_elements'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

        eIdx => elemI(:,ei_Lidx)

        !% rectangular channels, conduits and junction main
        ptype => col_elemPGetm(epg_CCJM_rectangular_nonsurcharged)
        npack => npack_elemPGetm(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemI(:,ei_geometryType) == rectangular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )

        if (npack > 0) then
            elemPGetm(1:npack, ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemI(:,ei_geometryType) == rectangular) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )
        endif

        !% trapezoidal channels, conduits and junction main
        ptype => col_elemPGetm(epg_CCJM_trapezoidal_nonsurcharged)
        npack => npack_elemPGetm(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemI(:,ei_geometryType) == trapezoidal) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )

        if (npack > 0) then
            elemPGetm(1:npack, ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemI(:,ei_geometryType) == trapezoidal) &
                .and. &
                (.not. elemYN(:,eYN_isSurcharged)) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )
        endif

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
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
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

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
        endif

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
        endif

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
        endif

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
        endif

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
        endif

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
        endif

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_nongeometry_static_elements
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine pack_nongeometry_dynamic_elements()
        !--------------------------------------------------------------------------
        !
        !% packed arrays for non geometry dynamic elements
        !
        !--------------------------------------------------------------------------

        integer          :: ii

        integer, pointer :: ptype, npack, eIDx(:)
        integer, allocatable :: fup(:), fdn(:)

        character(64) :: subroutine_name = 'pack_nongeometry_dynamic_elements'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

        eIdx => elemI(:,ei_Lidx)

        fup = pack(elemI(:,ei_Mface_uL), elemI(:,ei_Mface_uL) /= nullvalueI)
        fdn = pack(elemI(:,ei_Mface_dL), elemI(:,ei_Mface_dL) /= nullvalueI)

        !% ep_AC
        !% - all elements that use AC
        ptype => col_elemP(ep_AC)
        npack => npack_elemP(ptype)
        npack = count( &
                (elemI(:,ei_tmType) == AC))
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemI(:,ei_tmType) == AC))
        endif

        !% ep_CC_AC
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
        endif

        !% ep_CC_ETM
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
        endif

        !% ep_CC_H_ETM
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
        endif

        !% ep_CC_Q_AC
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
        endif

        !% ep_CC_Q_ETM
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
        endif

        !% ep_CCJB_AC
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
        endif


        !% ep_CC_AC_surcharged
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
        endif

        !% ep_CCJB_AC_surcharged
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
        endif

        !% ep_CC_ALLtm_surcharged
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
        endif

        !% ep_CCJB_ALLtm_surcharged
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
        endif

        !% ep_CCJB_eETM_i_fAC
        !% conduits, channels, and junction branches that are ETM and have
        !% an adjacent face that is AC
        ptype => col_elemP(ep_CCJB_eETM_i_fAC)
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
                ( &
                    (faceYN(fup,fYN_isAC_adjacent)) &
                    .or. &
                    (faceYN(fdn,fYN_isAC_adjacent)) &
                ))

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
                ( &
                    (faceYN(fup,fYN_isAC_adjacent)) &
                    .or. &
                    (faceYN(fdn,fYN_isAC_adjacent)) &
                ))
        endif

        !% ep_CCJB_ETM
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
        endif

        !% ep_CC_ETM_surcharged
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
        endif

        !% ep_CCJB_ETM_surcharged
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
        endif

        !% ep_CCJM_H_AC_open
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
        endif

        !% ep_CCJM_H_ETM
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
        endif

        !% ep_ETM
        !% - all elements that use ETM
        ptype => col_elemP(ep_ETM)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemI(:,ei_tmType) == ETM))
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx, &
                (elemI(:,ei_tmType) == ETM))
        endif

        !% ep_JM_AC
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
        endif

        !% ep_JB_AC
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
        endif

        !% ep_JM_ETM
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
        endif

        !% ep_JB_ETM
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
        endif

        !% ep_NonSurcharged_AC
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
        endif

        !% ep_NonSurcharged_ALLtm
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
        endif

        !% ep_NonSurcharged_ETM
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
        endif

        !NOT SURE IF THIS SHOULD BE DONE HERE OR WHERE SMALL VOLUMES ARE DECLARED
        !% ep_smallvolume_AC
        !% - all small volumes that are AC
        ptype => col_elemP(ep_smallvolume_AC)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemYN(:,eYN_isSmallVolume)) &
                .and. &
                (elemI(:,ei_tmType) == AC))

        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemYN(:,eYN_isSmallVolume)) &
                .and. &
                (elemI(:,ei_tmType) == AC))
        endif

        !NOT SURE IF THIS SHOULD BE DONE HERE OR WHERE SMALL VOLUMES ARE DECLARED
        !% ep_smallvolume_ALLtm
        !% - all small volumes that are any time march
        ptype => col_elemP(ep_smallvolume_ALLtm)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemYN(:,eYN_isSmallVolume)) &
                .and. &
                (   &
                    (elemI(:,ei_tmType) == AC) &
                    .or. &
                    (elemI(:,ei_tmType) == ETM) &
                ))

        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemYN(:,eYN_isSmallVolume)) &
                .and. &
                (   &
                    (elemI(:,ei_tmType) == AC) &
                    .or. &
                    (elemI(:,ei_tmType) == ETM) &
                ))
        endif

        !NOT SURE IF THIS SHOULD BE DONE HERE OR WHERE SMALL VOLUMES ARE DECLARED
        !% ep_smallvolume_ETM
        !% - all small volumes that are ETM
        ptype => col_elemP(ep_smallvolume_ETM)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemYN(:,eYN_isSmallVolume)) &
                .and. &
                (elemI(:,ei_tmType) == ETM))

        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemYN(:,eYN_isSmallVolume)) &
                .and. &
                (elemI(:,ei_tmType) == ETM))
        endif

        !% ep_Surcharged_AC
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
        endif

        !% ep_Surcharged_ALLtm
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
        endif

        !% ep_Surcharged_ETM
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
        endif

        if allocated(fup) deallocate(fup)
        if allocated(fdn) deallocate(fdn)

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_nongeometry_dynamic_elements
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine pack_static_interior_faces()
        !--------------------------------------------------------------------------
        !
        !% packed arrays for static faces
        !
        !--------------------------------------------------------------------------
        integer :: ii, image

        integer, pointer :: Nfaces, ptype, npack, fIdx(:), eup(:), edn(:)

        character(64) :: subroutine_name = 'pack_static_interior_faces'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

        !% pointing to the number of faces in this image
        image  = this_image()
        Nfaces => N_face(image)

        fIdx => faceI(1:Nfaces,fi_Lidx)
        eup  => faceI(1:Nfaces,fi_Melem_uL)
        edn  => faceI(1:Nfaces,fi_Melem_dL)

        ! % fp_all
        ! % - all faces execpt boundary, null, and shared faces
        ptype => col_faceP(fp_all)
        npack => npack_faceP(ptype)

        npack = count(faceYN(1:Nfaces,fYN_isInteriorFace))

        if (npack > 0) then
            faceP(1:npack,ptype) = pack(fIdx, &
                faceYN(1:Nfaces,fYN_isInteriorFace) &
                )
        endif


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
        endif

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_static_interior_faces
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine pack_dynamic_interior_faces()
        !--------------------------------------------------------------------------
        !
        !% packed arrays for dynamic faces
        !% HACK: Should the jump packing be called after all jump conditions are
        !% changed? or can it wait until the end of a time step? Note that this
        !% simply packs what is stored in faceI(:,fi_jump_type) as the actual
        !% computation of what is a jump is in the identify_hydraulic_jump subroutine.
        !
        !--------------------------------------------------------------------------

        integer          :: ii, image
        integer, pointer :: Nfaces, ptype, npack, fIdx(:), eup(:), edn(:)

        character(64) :: subroutine_name = 'pack_dynamic_interior_faces'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

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
        endif

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
        endif

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
        endif

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_dynamic_interior_faces
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine pack_static_shared_faces()
        !--------------------------------------------------------------------------
        !
        !% packed arrays for static shared faces
        !
        !--------------------------------------------------------------------------
        integer :: ii, image
        integer, pointer :: ptype, npack, fIdx(:), eup, edn, gup, gdn, Nfaces
        integer, pointer :: c_image, N_shared_faces, thisP
        logical, pointer :: isUpGhost, isDnGhost

        character(64) :: subroutine_name = 'pack_static_shared_faces'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

        !% pointing to the number of faces in this image
        image  = this_image()
        Nfaces => N_face(image)

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
        endif

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
                    endif

                elseif (isDnGhost) then
                    if ((elemI(gdn,ei_HeqType)[c_image] == diagnostic) .or.  &
                        (elemI(gdn,ei_QeqType)[c_image] == diagnostic) .or.  &
                        (elemI(eup,ei_HeqType)          == diagnostic) .or.  &
                        (elemI(eup,ei_QeqType)          == diagnostic))  then

                        !% advance the number of pack value
                        npack = npack + oneI
                        !% save the face index
                        facePS(npack,ptype) = thisP
                    endif
                endif
            end do
        endif

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_static_shared_faces
    !
    !==========================================================================
    !==========================================================================
    !
     subroutine pack_dynamic_shared_faces()
        !--------------------------------------------------------------------------
        !
        !% packed arrays for dynamic shared faces
        !% HACK: Should the jump packing be called after all jump conditions are
        !% changed? or can it wait until the end of a time step? Note that this
        !% simply packs what is stored in faceI(:,fi_jump_type) as the actual
        !% computation of what is a jump is in the identify_hydraulic_jump subroutine.
        !
        !--------------------------------------------------------------------------

        integer          :: ii, image
        integer, pointer :: ptype, npack, fIdx(:), Nfaces
        integer, pointer :: N_shared_faces, thisP, eup, edn, gup, gdn, c_image
        logical, pointer :: isUpGhost, isDnGhost
        character(64)    :: subroutine_name = 'pack_dynamic_shared_faces'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

        sync all

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
                    endif

                elseif (isDnGhost) then
                    if ((elemI(gdn,ei_tmType)[c_image] == AC) .or.  &
                        (elemI(eup,ei_tmType)          == AC))  then

                        !% advance the number of pack value
                        npack = npack + oneI
                        !% save the face index
                        facePS(npack,ptype) = thisP
                    endif
                endif
            end do
        endif

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
        endif

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
        endif

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_dynamic_shared_faces
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
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

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
        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_nodes

    subroutine pack_bc
        integer :: psize
        character(64) :: subroutine_name = 'pack_bc'
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

        !% BC packs
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
        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_bc
end module pack_mask_arrays
