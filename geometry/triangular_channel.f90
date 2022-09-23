module triangular_channel

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% triangular channel geometry
    !%

    private

    public :: triangular_depth_from_volume
    public :: triangular_area_from_depth
    public :: triangular_area_from_depth_singular
    public :: triangular_topwidth_from_depth
    public :: triangular_topwidth_from_depth_singular 
    public :: triangular_perimeter_from_depth
    public :: triangular_perimeter_from_depth_singular
    public :: triangular_hyddepth_from_depth
    public :: triangular_hyddepth_from_depth_singular
    public :: triangular_hydradius_from_depth_singular

    contains

!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine triangular_depth_from_volume (elemPGx, Npack, thisCol)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels 
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !% NOTE: this does NOT limit the depth by surcharge height at this point
        !% This will be done after the head is computed.
        !%------------------------------------------------------------------
            integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
            integer, pointer :: thisP(:)
            real(8), pointer :: depth(:), volume(:), length(:), breadth(:)
            real(8), pointer :: sideslope(:), fullvolume(:), fulldepth(:)
        !%-----------------------------------------------------------------
        thisP      => elemPGx(1:Npack,thisCol) 
        depth      => elemR(:,er_Depth)
        volume     => elemR(:,er_Volume)
        length     => elemR(:,er_Length)
        fulldepth  => elemR(:,er_FullDepth)
        fullvolume => elemR(:,er_FullVolume)
        breadth    => elemSGR(:,esgr_Triangular_TopBreadth)
        sideslope  => elemSGR(:,esgr_Triangular_Slope)
        !%-----------------------------------------------------------------------------  
        where (volume(thisP) >= fullvolume(thisP))
            depth(thisP) = fulldepth(thisP)
        elsewhere
            depth(thisP) = sqrt((volume(thisP) / length(thisP)) / sideslope(thisP))
        endwhere

    end subroutine triangular_depth_from_volume
    !%  
!%==========================================================================
!%==========================================================================
!%
    elemental real(8) function triangular_area_from_depth (indx) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for triangular cross section
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx  ! may be a packed array of indexes
        !%-----------------------------------------------------------------------------
        outvalue = elemSGR(indx,esgr_Triangular_Slope) * (elemR(indx,er_Depth) ** twoR)

    end function triangular_area_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function triangular_area_from_depth_singular (indx,depth) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for triangular cross section of a single element
        !% The input indx is the row index in full data 2D array.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer ::  breadth(:), sideslope(:)
        !%-----------------------------------------------------------------------------
        sideslope => elemSGR(:,esgr_Triangular_Slope)
        breadth => elemSGR(:,esgr_Triangular_TopBreadth)
        !%-----------------------------------------------------------------------------
        outvalue = sideslope(indx) * (depth ** twoR)

    end function triangular_area_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine triangular_topwidth_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a triangular channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: breadth(:), topwidth(:), sideslope(:), depth(:)
        !%-----------------------------------------------------------------------------
        thisP    => elemPGx(1:Npack,thisCol) 
        topwidth => elemR(:,er_Topwidth)
        sideslope => elemSGR(:,esgr_Triangular_Slope)
        depth    => elemR(:,er_Depth)
        !%-----------------------------------------------------------------------------

        topwidth(thisP) = twoR * sideslope(thisP) * depth(thisP) 

    end subroutine triangular_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function triangular_topwidth_from_depth_singular (indx,depth) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth for a triangular cross section of a single element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx 
        real(8), intent(in) :: depth
        real(8), pointer    :: sideslope(:)
        !%-----------------------------------------------------------------------------
        sideslope => elemSGR(:,esgr_Triangular_Slope)
        !%-----------------------------------------------------------------------------
        !%  
        outvalue = twoR * sideslope(indx) * depth

    end function triangular_topwidth_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine triangular_perimeter_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a triangular channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: slope(:), depth(:), perimeter(:)
        !%-----------------------------------------------------------------------------
        thisP     => elemPGx(1:Npack,thisCol) 
        slope     => elemSGR(:,esgr_Triangular_Slope)
        depth     => elemR(:,er_Depth)
        perimeter => elemR(:,er_Perimeter)
        !%-----------------------------------------------------------------------------

        perimeter(thisP) = twoR * depth(thisP) * sqrt(oneR + slope(thisP) ** twoR)

    end subroutine triangular_perimeter_from_depth
!%    
!%==========================================================================    
!%==========================================================================
!%
    real(8) function triangular_perimeter_from_depth_singular (indx,depth) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes wetted perimeter from known depth for a triangular cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer    :: slope(:)
        !%-----------------------------------------------------------------------------
        slope => elemSGR(:,esgr_Triangular_Slope)
        !%-----------------------------------------------------------------------------
        
        outvalue = twoR * depth * sqrt(1 + slope(indx) ** twoR)

    end function triangular_perimeter_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine triangular_hyddepth_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the hydraulic (average) depth from a known depth in a triangular channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: hyddepth(:), depth(:)
        !%-----------------------------------------------------------------------------
        thisP     => elemPGx(1:Npack,thisCol) 
        depth     => elemR(:,er_Depth)
        hyddepth  => elemR(:,er_HydDepth)
        !%-----------------------------------------------------------------------------

        hyddepth(thisP) = depth(thisP) / twoR

    end subroutine triangular_hyddepth_from_depth
!%    
!%==========================================================================  
!%==========================================================================
!%
    real(8) function triangular_hyddepth_from_depth_singular (indx,depth) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic depth from known depth for triangular cross section of 
        !% a single element
        !%-----------------------------------------------------------------------------   
        integer, intent(in) :: indx     
        real(8), intent(in) :: depth
        !%-----------------------------------------------------------------------------  

        outvalue = depth / twoR

    end function triangular_hyddepth_from_depth_singular 
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function triangular_hydradius_from_depth_singular (indx,depth) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic radius from known depth for a triangular cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer    :: breadth(:), sideslope(:)
        !%-----------------------------------------------------------------------------
        sideslope => elemSGR(:,esgr_Triangular_Slope)
        breadth   => elemSGR(:,esgr_Triangular_TopBreadth)
        !%-----------------------------------------------------------------------------
        
        outvalue = (sideslope(indx) * depth) / (twoR * sqrt(oneR + sideslope(indx) ** twoR))

    end function triangular_hydradius_from_depth_singular
    !%    
!%==========================================================================

end module triangular_channel