module irregular_channel

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use xsect_tables

    implicit none

    !%-----------------------------------------------------------------------
    !% Description:
    !% Irregular channel geometry
    !%

    private

    public :: irregular_depth_from_volume
    public :: irregular_topwidth_from_depth
    public :: irregular_hydradius_from_depth
    public :: irregular_perimeter_from_hydradius_area
    public :: irregular_hyddepth_from_topwidth_area

contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine irregular_depth_from_volume(elemPGx, Npack, thisCol)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the depth from volume for irregular channel
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !% NOTE: this does NOT limit the depth by surcharge height at this point
        !% This will be done after the head is computed.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
            integer, pointer :: thisP(:), tidx(:)
            real(8), pointer :: depth(:), fulldepth(:), volume(:), fullvolume(:)
            real(8), pointer :: thisTable(:,:), normInput(:)
        !%------------------------------------------------------------------
        !% Aliases
            thisP      => elemPGx(1:Npack,thisCol)
            depth      => elemR(:,er_Depth)
            fulldepth  => elemR(:,er_FullDepth)
            volume     => elemR(:,er_Volume)
            fullvolume => elemR(:,er_FullVolume)
            normInput  => elemR(:,er_Temp01)
            thisTable  => transectTableAreaR(:,:,tt_depth)
            tidx       => elemI(:,ei_transect_idx)
        !%------------------------------------------------------------------
        !% Preliminaries
            normInput(:) = zeroR
        !%------------------------------------------------------------------

        !% --- normalize the input
        normInput(thisP) = volume(thisP) / fullvolume(thisP)

        !% --- lookup the normalized depth
        call xsect_table_lookup_array (depth, normInput, thisTable, thisP)

        !% --- convert to physical depth
        depth(thisP) = depth(thisP) * transectR(tidx(thisP),tr_depthFull)

        !%------------------------------------------------------------------
        !% Closing
        !% --- reset the temp array
            normInput(:) = zeroR

    end subroutine irregular_depth_from_volume
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine irregular_topwidth_from_depth(elemPGx, Npack, thisCol)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from depth for irregular channel
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !% NOTE: this does NOT limit the depth by surcharge height at this point
        !% This will be done after the head is computed.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
            integer, pointer :: thisP(:), tidx(:)
            real(8), pointer :: depth(:), fulldepth(:), topwidth(:)
            real(8), pointer :: thisTable(:,:), normInput(:)
        !%------------------------------------------------------------------
        !% Aliases
            thisP      => elemPGx(1:Npack,thisCol)
            depth      => elemR(:,er_Depth)
            fulldepth  => elemR(:,er_FullDepth)
            topwidth   => elemR(:,er_TopWidth)
            normInput  => elemR(:,er_Temp01)
            thisTable  => transectTableDepthR(:,:,tt_width)
            tidx       => elemI(:,ei_transect_idx)
        !%------------------------------------------------------------------
        !% Preliminaries
            normInput(:) = zeroR
        !%------------------------------------------------------------------

        !% --- normalize the input
        normInput(thisP) = depth(thisP) / fulldepth(thisP)

        !% --- lookup the topwidth
        call xsect_table_lookup_array (topwidth, normInput, thisTable, thisP)

        !% --- convert to physical topwidth
        topwidth(thisP) = topwidth(thisP) * transectR(tidx(thisP),tr_widthMax)

        !%------------------------------------------------------------------
        !% Closing
        !% --- reset the temp array
            normInput(:) = zeroR

    end subroutine irregular_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine irregular_hydradius_from_depth(elemPGx, Npack, thisCol)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the hydraulic radius from depth for irregular channel
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !% NOTE: this does NOT limit the depth by surcharge height at this point
        !% This will be done after the head is computed.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
            integer, pointer :: thisP(:), tidx(:)
            real(8), pointer :: depth(:), fulldepth(:), hydradius(:)
            real(8), pointer :: thisTable(:,:), normInput(:)
        !%------------------------------------------------------------------
        !% Aliases
            thisP          => elemPGx(1:Npack,thisCol)
            depth          => elemR(:,er_Depth)
            fulldepth      => elemR(:,er_FullDepth)
            hydradius      => elemR(:,er_HydRadius)
            normInput      => elemR(:,er_Temp01)
            thisTable      => transectTableDepthR(:,:,tt_width)
            tidx           => elemI(:,ei_transect_idx)
        !%------------------------------------------------------------------
        !% Preliminaries
            normInput(:) = zeroR
        !%------------------------------------------------------------------

        !% --- normalize the input
        normInput(thisP) = depth(thisP) / fulldepth(thisP)

        !% --- lookup the topwidth
        call xsect_table_lookup_array (hydradius, normInput, thisTable, thisP)

        !% --- convert to physical topwidth
        hydradius(thisP) = hydradius(thisP) * transectR(tidx(thisP),tr_hydRadiusFull)

        !%------------------------------------------------------------------
        !% Closing
        !% --- reset the temp array
            normInput(:) = zeroR

    end subroutine irregular_hydradius_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine irregular_perimeter_from_hydradius_area (elemPGx, Npack, thisCol)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from hydraulic radius and area for
        !% irregular channel
        !% This does NOT require a lookup because we use lookups of the
        !% perimeter and area before
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
            integer, pointer :: thisP(:)
            real(8), pointer :: perimeter(:), area(:), hydradius(:)
        !%------------------------------------------------------------------
        !% Aliases
            thisP     => elemPGx(1:Npack,thisCol)
            perimeter => elemR(:,er_Perimeter)
            area      => elemR(:,er_area)
            hydradius => elemR(:,er_HydRadius)
        !%------------------------------------------------------------------
            perimeter(thisP) = hydradius(thisP) / area(thisP)

    end subroutine irregular_perimeter_from_hydradius_area
!%    
!%==========================================================================
!%==========================================================================
!%  
    subroutine irregular_hyddepth_from_topwidth_area (elemPGx, Npack, thisCol)  
        !%------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic depth as A/T
        !%------------------------------------------------------------------
        !% Declarations:
        integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: topwidth(:), area(:), hyddepth(:)
        !%------------------------------------------------------------------
        !% Aliases
            thisP     => elemPGx(1:Npack,thisCol)
            topwidth  => elemR(:,er_TopWidth)
            area      => elemR(:,er_area)
            hyddepth  => elemR(:,er_HydDepth)
        !%------------------------------------------------------------------
            hyddepth(thisP) = area(thisP) / topwidth(thisP)

    end subroutine irregular_hyddepth_from_topwidth_area
!%    
!%==========================================================================
!%==========================================================================
!%    
!%==========================================================================
!% END MODULE
!%==========================================================================
!%
!%  
    !%----------------------------------------------------------------------
    !% Description:
    !% 
    !%----------------------------------------------------------------------

    !%----------------------------------------------------------------------
    !%  
end module irregular_channel