module irregular_channel

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use geometry_lowlevel
    use xsect_tables

    implicit none

    !%-----------------------------------------------------------------------
    !% Description:
    !% Irregular channel geometry
    !%

    private

    public :: irregular_depth_from_volume
    public :: irregular_topwidth_from_depth
    public :: irregular_perimeter_and_hydradius_from_depth


    public :: irregular_perimeter_from_hydradius_area
    ! public :: irregular_hyddepth_from_topwidth_area
    public :: irregular_geometry_from_depth_singular

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
        !% Assumes that volume > 0 is previously enforced in volume computations.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
            integer, pointer :: thisP(:), tidx(:)
            real(8), pointer :: depth(:), fulldepth(:), volume(:), fullvolume(:)
            real(8), pointer :: thisTable(:,:), normInput(:)
        !%------------------------------------------------------------------
            if (Npack < 1) return
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

        !% --- correct the depth for overfull channel
        if (setting%Discretization%AllowChannelOverflowTF) then
            !% --- volume is already limited to full volume
            !%     so no correction needed
        else
            !% --- if overflow is not allowed, use rectangular volume abovev
            !%     for overfull elements
            where  (volume(thisP) >= fullvolume(thisP))
                !% --- volume above max level is rectangular at max breadth
                depth(thisP) = llgeo_openchannel_depth_above_full_pure(thisP)
            endwhere
        end if

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
        !% Assumes that volume > 0 is previousy enforced in volume computations.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
            integer, pointer :: thisP(:), tidx(:)
            real(8), pointer :: depth(:), fulldepth(:), topwidth(:)
            real(8), pointer :: volume(:), fullvolume(:)
            real(8), pointer :: thisTable(:,:), normInput(:)
        !%------------------------------------------------------------------
        !% Aliases
            thisP      => elemPGx(1:Npack,thisCol)
            depth      => elemR(:,er_Depth)
            fulldepth  => elemR(:,er_FullDepth)
            topwidth   => elemR(:,er_TopWidth)
            volume     => elemR(:,er_Volume)
            fullvolume => elemR(:,er_FullVolume)
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

        !% --- correct the depth for overfull channel
        if (setting%Discretization%AllowChannelOverflowTF) then
            !% --- depth is already limited to full depth
            !%     so no correction needed
        else
            !% --- if overflow is not allowed, then volume > full volume is possible
            !%     use rectangular shape to correct
            !%     topwidth for depth greater than full depth
            where  (volume(thisP) >= fullvolume(thisP))
                !% --- volume above max level is rectangular at max breadth
                topwidth(thisP) = llgeo_openchannel_depth_above_full_pure(thisP)
            endwhere
        end if

        !%------------------------------------------------------------------
        !% Closing
        !% --- reset the temp array
            normInput(:) = zeroR

    end subroutine irregular_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine irregular_perimeter_and_hydradius_from_depth(elemPGx, Npack, thisCol)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the hydraulic radius from depth for irregular channel
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is previuosly enforced in volume computations.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
            integer, pointer :: thisP(:), tidx(:)
            real(8), pointer :: depth(:), fulldepth(:), hydradius(:)
            real(8), pointer :: perimeter(:), area(:), volume(:), fullvolume(:)
            real(8), pointer :: thisTable(:,:), normInput(:)
        !%------------------------------------------------------------------
        !% Aliases
            thisP          => elemPGx(1:Npack,thisCol)
            depth          => elemR(:,er_Depth)
            fulldepth      => elemR(:,er_FullDepth)
            hydradius      => elemR(:,er_HydRadius)
            perimeter      => elemR(:,er_Perimeter)
            volume         => elemR(:,er_Volume)
            fullvolume     => elemR(:,er_FullVolume)
            area           => elemR(:,er_Area)
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

        !% --- correct for overfull channel
        if (setting%Discretization%AllowChannelOverflowTF) then
            !% --- depth is already limited to full depth
            !%     so no correction needed to hydradius
            perimeter(thisP) = llgeo_perimeter_from_hydradius_and_area_pure &
                            (thisP, hydradius(thisP), area(thisP))
        else
            !% --- if overflow is not allowed, then volume > full volume is possible
            !%     use rectangular shape to correct
            !%     thydraidius and perimeter for depth greater than full depth
            where  (volume(thisP) >= fullvolume(thisP))
                !% --- use the full perimeter and volume excess to get the new perimeter
                perimeter(thisP) = llgeo_openchannel_perimeter_above_full_pure(thisP)
                !% --- note that the area with an extra volume already includes the
                !%     extra volume effect.
                hydradius(thisP) = llgeo_hydradius_from_area_and_perimeter_pure &
                                        (thisP, area(thisP), perimeter(thisP))
            endwhere
        end if

        !%------------------------------------------------------------------
        !% Closing
        !% --- reset the temp array
            normInput(:) = zeroR

    end subroutine irregular_perimeter_and_hydradius_from_depth
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
    ! subroutine irregular_hyddepth_from_topwidth_area (elemPGx, Npack, thisCol)  
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Computes hydraulic depth as A/T
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !     integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
    !     integer, pointer :: thisP(:)
    !     real(8), pointer :: topwidth(:), area(:), hyddepth(:)
    !     !%------------------------------------------------------------------
    !     !% Aliases
    !         thisP     => elemPGx(1:Npack,thisCol)
    !         topwidth  => elemR(:,er_TopWidth)
    !         area      => elemR(:,er_area)
    !         hyddepth  => elemR(:,er_HydDepth)
    !     !%------------------------------------------------------------------
    !         hyddepth(thisP) = area(thisP) / topwidth(thisP)

    ! end subroutine irregular_hyddepth_from_topwidth_area
!%    
!%==========================================================================
!%==========================================================================
!%   
    real(8) function irregular_geometry_from_depth_singular &
                (indx, table_idx, depth, ZeroValueGeometry) result (outvalue)
        !%----------------------------------------------------------------------
        !% Description:
        !%  Computes cross-sectional geometry for a given depth for a single element (indx)
        !%  depth is the unormalized depth of interest
        !%  table_idx is in {tt_area, tt_width, tt_depth, tt_hydradius}
        !%  zerovalue is the result used in place of zeros
        !%----------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: indx, table_idx
            real(8), intent(in) :: depth, ZeroValueGeometry
            real(8), pointer    :: fullDepth(:), thisTable(:)
            real(8)             :: depthnorm
        !%----------------------------------------------------------------------
        !% Aliases
            !% --- get the transect by depth table 
            thisTable => transectTableDepthR(elemI(indx,ei_transect_idx),:,table_idx)
            fullDepth => elemR(:,er_FullDepth)
        !%----------------------------------------------------------------------
        !% --- normalized depth
        depthnorm     = depth/fulldepth(indx)

        !% --- max is used because xsect quadratic interp for small values can produce zero
        outvalue = max( xsect_table_lookup_singular (depthnorm, thisTable(:)),ZeroValueGeometry)

    end function irregular_geometry_from_depth_singular    
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