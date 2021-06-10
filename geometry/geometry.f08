module geometry

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use rectangular


    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Geometry computations
    !%

    private

    public :: geometry_toplevel 

    contains
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    !%
    subroutine geometry_toplevel (whichTM)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Here whichTM is one of ETM, AC, or ALLtm
        !% This should never be called for diagnostic arrays
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: whichTM
        integer, pointer :: elemPGx(:,:), npack_elemPGx(:), col_elemPGx(:)
        integer, pointer :: thisCol_surcharged, thisCol_all, thisCol_JM, thisCol_NonSurcharged
        integer, pointer :: Npack, thisCol
        !%-----------------------------------------------------------------------------
        !% set the packed element array to use
        select case (whichTM)
            case (ALLtm)
                elemPGx               => elemPGalltm
                npack_elemPGx         => npack_elemPGalltm
                col_elemPGx           => col_elemPGalltm
                thisCol_surcharged    => col_elemPGalltm(ep_Surcharged_ALLtm)
                thisCol_all           => col_elemPGalltm(ep_ALLtm)
                thisCol_JM            => col_elemPGalltm(ep_JM_ALLtm)  
                thisCol_NonSurcharged => col_elemPGalltm(ep_NonSurcharged_ALLtm)
            case (ETM)
                elemPGx               => elemPGetm
                npack_elemPGx         => npack_elemPGetm
                col_elemPGx           => col_elemPGetm
                thisCol_surcharged    => col_elemPGetm(ep_Surcharged_ETM)        
                thisCol_all           => col_elemPGetm(ep_ETM)
                thisCol_JM            => col_elemPGetm(ep_JM_ETM)    
                thisCol_NonSurcharged => col_elemPGetm(ep_NonSurcharged_ETM)        
            case (AC)
                elemPGx               => elemPGac
                npack_elemPGx         => npack_elemPGac 
                col_elemPGx           => col_elemPGac         
                thisCol_surcharged    => col_elemPGac(ep_Surcharged_AC)
                thisCol_all           => col_elemPGac(ep_AC)
                thisCol_JM            => col_elemPGac(ep_JM_AC)
                thisCol_NonSurcharged => col_elemPGac(ep_NonSurcharged_AC)
            case default
                print *, 'error, case default should never be reached.'
                stop 7389
        end select

        !% compute the head on all non-surcharged elements of CC and JM
        call geo_head_from_volume_CCJM (elemPGx, npack_elemPGx, col_elemPGx)

        !% assign the head on junction branches JB
        Npack => npack_elemPGx(thisCol_JM)
        if (Npack > 0) then
    !        call geometry_head_assign_JB (Npack, thisCol_JM)
        endif    

        !% assign all geometry for surcharged elements CC, JM, JB
        Npack => npack_elemPGx(thisCol_surcharged)
        if (Npack > 0) then
    !        call geometry_surcharged_common (Npack, thisCol_surcharged)
        endif

        !% compute area from volume and max depth from head for all Nonsurcharged
        Npack => npack_elemPGx(thisCol_NonSurcharged)
        if (Npack > 0) then
    !        call geometry_area_from_volume_common (Npack, thisCol_NonSurcharged)
    !        call geometry_depth_from_head_common (Npack, thisCol_NonSurcharged)
        endif

        !% compute topwidth from maximum depth for all CC, JM, JB
    !    call geometry_topwidth_from_depth (elemPGx, npack_elemPGx, col_elemPGx)

        !% compute perimeter from maximum depth for all CC, JM, JB
    !    call geometry_perimeter_from_depth (elemPGx, npack_elemPGx, col_elemPGx)

        !% compute hyddepth and hydradius for all Nonsurcharged
        Npack => npack_elemPGx(thisCol_NonSurcharged)
        if (Npack > 0) then
    !        call geometry_hyddepth_common (Npack, thisCol_NonSurcharged)
    !        call geometry_hydradius_common (Npack, thisCol_NonSurcharged)   
        endif

        !% the modified hydraulic depth "ell" is used for AC computations and
        !% for Froude number computations on all elements, whether ETM or AC.
        Npack => npack_elemPGx(thisCol_all)
        if (Npack > 0) then 
    !        call geometry_ell (Npack, thisCol_all)
        endif

        !% compute the dHdA that are only required for the AC method.
        !% Note that dHdA only for non-surcharged elements
        if (whichTM .ne. ETM) then
            !% computing dHdA for AC method non-surcharged
            thisCol => col_elemPGx(ep_NonSurcharged_AC)
            Npack   => npack_elemPGx(thisCol)
            if (Npack > 0) then 
    !            call geometry_dHdA (Npack, thisCol)
            endif    
        endif
  
    end subroutine geometry_toplevel
    !%
    !%==========================================================================
    !% PRIVATE
    !%==========================================================================   
    !%
    subroutine geo_head_from_volume_CCJM (elemPGx, npack_elemPGx, col_elemPGx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% This is limited to solving nonsurcharged CCJM elements
        !% The elemPGx determines whether this is ALLtm, ETM or AC elements
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:), npack_elemPGx(:), col_elemPGx(:)
        integer, pointer :: Npack, thisCol
        !%-----------------------------------------------------------------------------
        !% cycle through different geometries
        thisCol => col_elemPGx(epg_CCJM_rectangular_nonsurcharged)
        Npack   => npack_elemPGx(thisCol)
        if (Npack > 0) then
            call rectangular_open_head_from_volume (elemPGx, Npack, thisCol)    
        endif

        !HACK Needs additional geometries, including surcharged rectangular conduits
        ! with and without Preissman slot.

    end subroutine geo_head_from_volume_CCJM
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%   
    !subroutine 
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%  
    !end subroutine  
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    !subroutine
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------

    !%     
        !%-----------------------------------------------------------------------------
        !% Description:
        !%
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%  
    !%
    !%==========================================================================
    !% END OF MODULE
    !%+=========================================================================
end module geometry