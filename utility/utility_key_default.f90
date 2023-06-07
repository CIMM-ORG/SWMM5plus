module utility_key_default
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% sets the value for any KEY column in an integer array to the undefinedKey
    !% value. 
    !% See define_indexes for KEY columns of ***I arrays.
    !% See define_keys for the list of keys
    !%==========================================================================

    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting

    implicit none
    
    private

    public :: util_key_default_linknode
    public :: util_key_default_bc
    public :: util_key_default_elemX
    public :: util_key_default_face

    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine util_key_default_linknode ()
        !%------------------------------------------------------------------
        !% Description
        !% sets default undefinedKey for link%I and node%I arrays
        !%------------------------------------------------------------------
            integer, allocatable :: keylist(:)
        !%------------------------------------------------------------------
        !% --- for links
        !%     keylist must match the number of KEY columns for li_ indexes
        allocate(keylist(4))  
        keylist = (/    li_link_type, &
                        li_link_sub_type, &
                        li_geometry, &
                        li_InitialDepthType /)

        link%I(:,keylist) = undefinedKey
        deallocate(keylist)

        !% --- for nodes
        !%     keylist must match the number of KEY columns for ni_ indexes
        allocate(keylist(1))
        keylist(1) = ni_node_type
        node%I(:,keylist) = undefinedKey
        deallocate(keylist)

    end subroutine util_key_default_linknode
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_key_default_bc()
        !%------------------------------------------------------------------
        !% Description
        !% sets default undefinedKey for BC%headI and BC%flowI arrays
        !%------------------------------------------------------------------
        !% Declarations
            integer, allocatable :: keylist(:)
        !%------------------------------------------------------------------
        !% --- keylist must match the number of KEY columns for bi_ indexes
        allocate(keylist(2))  
        keylist = (/    bi_category, &
                        bi_subcategory /)

        BC%headI(:,keylist) = undefinedKey
        BC%flowI(:,keylist) = undefinedKey
        deallocate(keylist)

    end subroutine util_key_default_bc
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_key_default_elemX ()
        !%------------------------------------------------------------------
        !% Description
        !% sets default undefinedKey for key items in elemI elemSI arrays
        !%------------------------------------------------------------------
        !% Declarations
            integer, allocatable :: keylist(:)
        !%------------------------------------------------------------------
        !% --- keylist must match the number of KEY columns for ei_ indexes
        allocate(keylist(5))  
        keylist = (/    ei_elementType, &
                        ei_geometryType, &
                        ei_HeqType, &
                        ei_QeqType, &
                        ei_tmType /)
        elemI(:,keylist) = undefinedKey
        deallocate(keylist)

        !% --- keylist must match the number of KEY columns for esi_JunctionMain indexes
        allocate(keylist(1))  
        keylist(1) = esi_JunctionMain_Type 
        elemSI(:,keylist) = undefinedKey
        deallocate(keylist)

        !% --- keylist must match the number of KEY columns for esi_Weir indexes
        allocate(keylist(2))  
        keylist = (/ esi_Weir_SpecificType, &
                     esi_Weir_GeometryType /)
        elemSI(:,keylist) = undefinedKey
        deallocate(keylist)

        !% --- keylist must match the number of KEY columns for esi_Orifice indexes
        allocate(keylist(2))  
        keylist = (/ esi_Orifice_SpecificType, &
                     esi_Orifice_GeometryType /)
        elemSI(:,keylist) = undefinedKey
        deallocate(keylist)

        !% --- keylist must match the number of KEY columns for esi_Outlet indexes
        allocate(keylist(1))  
        keylist(1) = esi_Outlet_SpecificType 
        elemSI(:,keylist) = undefinedKey
        deallocate(keylist)

        !% --- keylist must match the number of KEY columns for esi_Pump indexes
        allocate(keylist(1))  
        keylist(1) = esi_Pump_SpecificType 
        elemSI(:,keylist) = undefinedKey
        deallocate(keylist)

    end subroutine util_key_default_elemX
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_key_default_face ()
        !%------------------------------------------------------------------
        !% Description
        !% sets default undefinedKey for faceI arrays
        !%------------------------------------------------------------------
        integer, allocatable :: keylist(:)
        !%------------------------------------------------------------------
        !% --- keylist must match the number of KEY columns for fi_ indexes
        allocate(keylist(2))  
        keylist = (/    fi_BCtype, &
                        fi_jump_type /)

        faceI(:,keylist) = undefinedKey
        deallocate(keylist)

    end subroutine util_key_default_face
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
!%
end module utility_key_default