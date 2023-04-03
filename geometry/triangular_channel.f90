module triangular_channel

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use geometry_lowlevel

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% triangular channel geometry
    !%

    private

    public :: triangular_depth_from_volume
    public :: triangular_topwidth_from_depth
    public :: triangular_perimeter_from_depth

    ! public :: triangular_area_from_depth_singular
    ! public :: triangular_topwidth_from_depth_singular 
    ! public :: triangular_perimeter_from_depth_singular
    
    !public :: triangular_hyddepth_from_depth
    
    !public :: triangular_hyddepth_from_depth_singular
    !public :: triangular_hydradius_from_depth_singular

    contains

!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine triangular_depth_from_volume (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels 
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is previuosly enforced in volume computations.
        !%------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: depth(:), volume(:)
            real(8), pointer ::fullvolume(:), fulldepth(:)
        !%------------------------------------------------------------------
        !% Aliases:
            depth      => elemR(:,er_Depth)
            volume     => elemR(:,er_Volume)
            fulldepth  => elemR(:,er_FullDepth)
            fullvolume => elemR(:,er_FullVolume)
        !%------------------------------------------------------------------  

        if (setting%Discretization%AllowChannelOverflowTF) then
            where (volume(thisP) >= fullvolume(thisP))
                !% --- truncate depth at full depth if there is overflow
                depth(thisP) = fulldepth(thisP)
            elsewhere
                !% --- standard triangular depth
                depth(thisP) = llgeo_triangular_depth_from_volume_pure &
                                    (thisP, volume(thisP))
            endwhere
        else
            !% --- if overflow is NOT allowed
            where (volume(thisP) >= fullvolume(thisP))
                !% --- volume above max level is rectangular at max breadth
                depth(thisP) = llgeo_openchannel_depth_above_full_pure(thisP)
            elsewhere
                !% --- standard triangular depth
                depth(thisP) = llgeo_triangular_depth_from_volume_pure &
                                    (thisP, volume(thisP))
            endwhere
        end if

    end subroutine triangular_depth_from_volume
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine triangular_topwidth_from_depth (thisP)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a triangular channel
        !%-----------------------------------------------------------------------------
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: topwidth(:), depth(:)
            real(8), pointer :: volume(:), fullvolume(:)
        !%------------------------------------------------------------------
        !% Aliases    
            topwidth    => elemR(:,er_Topwidth)
            volume      => elemR(:,er_Volume)
            fullvolume  => elemR(:,er_FullVolume)
            depth       => elemR(:,er_Depth)
        !%-----------------------------------------------------------------------------

        if (setting%Discretization%AllowChannelOverflowTF) then 
            !% --- depth is already limited to full depth
            topwidth(thisP) = llgeo_triangular_topwidth_from_depth_pure(thisP,depth(thisP))
        else
            where (volume(thisP) >= fullvolume(thisP))
                topwidth(thisP) = llgeo_openchannel_topwidth_above_full_pure(thisP)
            elsewhere
                !% --- use elemental form as depth <= fulldepth
                topwidth(thisP) = llgeo_triangular_topwidth_from_depth_pure(thisP,depth(thisP))
            endwhere
        endif

        !topwidth(thisP) = twoR * sideslope(thisP) * depth(thisP) 

    end subroutine triangular_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine triangular_perimeter_from_depth (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a triangular channel
        !%------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: volume(:), fullvolume(:), perimeter(:)
            real(8), pointer :: depth(:)
        !%------------------------------------------------------------------
        !% Aliases
            fullvolume => elemR(:,er_FullVolume)
            volume     => elemR(:,er_Volume)
            perimeter  => elemR(:,er_Perimeter)
            depth      => elemR(:,er_Depth)
        !%------------------------------------------------------------------

        if (setting%Discretization%AllowChannelOverflowTF) then 
            !% --- depth is already limited to full depth
            perimeter(thisP) = llgeo_triangular_perimeter_from_depth_pure(thisP,depth(thisP))
        else
            where (volume(thisP) >= fullvolume(thisP))
                perimeter(thisP) = llgeo_openchannel_perimeter_above_full_pure(thisP)
            elsewhere
                !% --- use elemental form as depth <= fulldepth
                perimeter(thisP) = llgeo_triangular_perimeter_from_depth_pure(thisP,depth(thisP))
            endwhere
        end if

    end subroutine triangular_perimeter_from_depth
!%    
!%==========================================================================  
!% SINGULAR
!%==========================================================================
!%
!     pure real(8) function triangular_area_from_depth_singular &
!         (indx,depth) result (outvalue)
!         !%------------------------------------------------------------------
!         !% Description:
!         !% Computes area from known depth for triangular cross section of a single element
!         !% The input indx is the row index in full data 2D array.
!         !%------------------------------------------------------------------
!             integer, intent(in) :: indx
!             real(8), intent(in) :: depth
!         !%------------------------------------------------------------------

!         outvalue = elemSGR(indx,esgr_Triangular_Slope) * (depth ** 2)

!     end function triangular_area_from_depth_singular
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     pure real(8) function triangular_topwidth_from_depth_singular &
!         (indx,depth) result (outvalue)
!         !%------------------------------------------------------------------
!         !% Description:
!         !% Computes the topwidth for a triangular cross section of a single element
!         !%------------------------------------------------------------------
!             integer, intent(in) :: indx 
!             real(8), intent(in) :: depth
!         !%------------------------------------------------------------------
 
!         outvalue = twoR * elemSGR(indx,esgr_Triangular_Slope) * depth

!     end function triangular_topwidth_from_depth_singular
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     pure real(8) function triangular_perimeter_from_depth_singular &
!         (indx,depth) result (outvalue)
!         !%  
!         !%-----------------------------------------------------------------
!         !% Description:
!         !% Computes wetted perimeter from known depth for a triangular cross section of
!         !% a single element 
!         !%------------------------------------------------------------------
!             integer, intent(in) :: indx
!             real(8), intent(in) :: depth
!         !%------------------------------------------------------------------
        
!         outvalue = twoR * depth * sqrt(oneR + elemSGR(indx,esgr_Triangular_Slope) ** 2)

!     end function triangular_perimeter_from_depth_singular
! !%    
! !%==========================================================================
! !%==========================================================================
! !%
!     ! subroutine triangular_hyddepth_from_depth (elemPGx, Npack, thisCol)
!     !     !%  
!     !     !%-----------------------------------------------------------------------------
!     !     !% Description:
!     !     !% Computes the hydraulic (average) depth from a known depth in a triangular channel
!     !     !%-----------------------------------------------------------------------------
!     !     integer, target, intent(in) :: elemPGx(:,:)
!     !     integer, intent(in) ::  Npack, thisCol
!     !     integer, pointer :: thisP(:)
!     !     real(8), pointer :: hyddepth(:), depth(:)
!     !     !%-----------------------------------------------------------------------------
!     !     thisP     => elemPGx(1:Npack,thisCol) 
!     !     depth     => elemR(:,er_Depth)
!     !     hyddepth  => elemR(:,er_HydDepth)
!     !     !%-----------------------------------------------------------------------------

!     !     hyddepth(thisP) = depth(thisP) / twoR

!     ! end subroutine triangular_hyddepth_from_depth
! !%    
! !%==========================================================================  
! !%==========================================================================
! !%
!     ! real(8) function triangular_hyddepth_from_depth_singular (indx,depth) result (outvalue)
!     !     !%  
!     !     !%-----------------------------------------------------------------------------
!     !     !% Description:
!     !     !% Computes hydraulic depth from known depth for triangular cross section of 
!     !     !% a single element
!     !     !%-----------------------------------------------------------------------------   
!     !     integer, intent(in) :: indx     
!     !     real(8), intent(in) :: depth
!     !     !%-----------------------------------------------------------------------------  

!     !     outvalue = depth / twoR

!     ! end function triangular_hyddepth_from_depth_singular 
! !%    
! !%==========================================================================
! !%==========================================================================
! !%
!     ! real(8) function triangular_hydradius_from_depth_singular (indx,depth) result (outvalue)
!     !     !%  
!     !     !%-----------------------------------------------------------------------------
!     !     !% Description:
!     !     !% Computes hydraulic radius from known depth for a triangular cross section of
!     !     !% a single element 
!     !     !%-----------------------------------------------------------------------------
!     !     integer, intent(in) :: indx
!     !     real(8), intent(in) :: depth
!     !     real(8), pointer    :: breadth(:), sideslope(:)
!     !     !%-----------------------------------------------------------------------------
!     !     sideslope => elemSGR(:,esgr_Triangular_Slope)
!     !     breadth   => elemSGR(:,esgr_Triangular_TopBreadth)
!     !     !%-----------------------------------------------------------------------------
        
!     !     outvalue = (sideslope(indx) * depth) / (twoR * sqrt(oneR + sideslope(indx) ** twoR))

!     ! end function triangular_hydradius_from_depth_singular
    !%    
!%==========================================================================

end module triangular_channel