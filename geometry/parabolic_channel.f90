module parabolic_channel

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% parabolic channel geometry
    !%

    private

    public :: parabolic_depth_from_volume
    public :: parabolic_area_from_depth
    public :: parabolic_area_from_depth_singular
    public :: parabolic_topwidth_from_depth
    public :: parabolic_topwidth_from_depth_singular 
    public :: parabolic_perimeter_from_depth
    public :: parabolic_perimeter_from_depth_singular
    public :: parabolic_hyddepth_from_depth
    !public :: parabolic_hyddepth_from_depth_singular
    public :: parabolic_hydradius_from_depth_singular

    contains

!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine parabolic_depth_from_volume (elemPGx, Npack, thisCol)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !% NOTE: this does NOT limit the depth by surcharge height at this point
        !% This will be done after the head is computed.
        !%--------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
            integer, pointer :: thisP(:)
            real(8), pointer :: depth(:), volume(:), length(:), breadth(:)
            real(8), pointer :: fulldepth(:), fullvolume(:), rbot(:)
        !%-------------------------------------------------------------------
            if (Npack < 1) return
        !%-------------------------------------------------------------------
            thisP      => elemPGx(1:Npack,thisCol) 
            depth      => elemR(:,er_Depth)
            volume     => elemR(:,er_Volume)
            length     => elemR(:,er_Length)
            fulldepth  => elemR(:,er_FullDepth)
            fullvolume => elemR(:,er_FullVolume)
            breadth    => elemSGR(:,esgr_Parabolic_Breadth)
            rbot       => elemSGR(:, esgr_Parabolic_Radius)
        !%-------------------------------------------------------------------- 
        where (volume(thisP) >= fullvolume(thisP))
            depth(thisP) = fulldepth(thisP)
        elsewhere
            depth(thisP) = ( threeR/fourR * (volume(thisP)/length(thisP)) &
                         / rbot(thisP) ) ** twothirdR
        endwhere
    end subroutine parabolic_depth_from_volume
    !%  
!%==========================================================================
!%==========================================================================
!%
    elemental real(8) function parabolic_area_from_depth (indx) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for parabolic cross section
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx  ! may be a packed array of indexes
        !%-----------------------------------------------------------------------------
        
        outvalue = (fourR / threeR) * elemSGR(indx, esgr_Parabolic_Radius) * elemR(indx, er_Depth) *  sqrt(elemR(indx, er_Depth))
        
    end function parabolic_area_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine parabolic_topwidth_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a parabolic channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: topwidth(:), depth(:), rbot(:)
        !%-----------------------------------------------------------------------------
        thisP     => elemPGx(1:Npack,thisCol) 
        topwidth  => elemR(:,er_Topwidth)
        depth     => elemR(:,er_Depth)
        rbot      => elemSGR(:, esgr_Parabolic_Radius)
        !%-----------------------------------------------------------------------------

        topwidth(thisP) = twoR * rbot(thisP) * sqrt(depth(thisP))
        
    end subroutine parabolic_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine parabolic_perimeter_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a parabolic channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:), perimeter(:), rbot(:), xx(:), tt(:)
        !%-----------------------------------------------------------------------------
        thisP     => elemPGx(1:Npack,thisCol) 
        depth     => elemR(:,er_Depth)
        perimeter => elemR(:,er_Perimeter)
        rbot      => elemSGR(:,esgr_Parabolic_Radius)
        xx         => elemR(:, er_Temp02)
        tt         => elemr(:, er_Temp03)
        !%-----------------------------------------------------------------------------
        xx(thisP) = twoR * sqrt(depth(thisP)) / rbot(thisP)
        tt(thisP) = sqrt(oneR + xx(thisP) * xx(thisP))
        !%-----------------------------------------------------------------------------

        perimeter(thisP) = onehalfR * rbot(thisP) * rbot(thisP) * (xx(thisP) * tt(thisP) + log(xx(thisP) + tt(thisP)))

    end subroutine parabolic_perimeter_from_depth
!%    
!%==========================================================================  
!%==========================================================================
!%
    subroutine parabolic_hyddepth_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the hydraulic (average) depth from a known depth in a parabolic channel
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

        hyddepth(thisP) = (twoR/threeR) * depth(thisP)

    end subroutine parabolic_hyddepth_from_depth
!%    
!%==========================================================================  
!% SINGULAR
!%==========================================================================
!%
    real(8) function parabolic_area_from_depth_singular (indx, depth) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for parabolic cross section of a single element
        !% The input indx is the row index in full data 2D array.
        !%-----------------------------------------------------------------------------
            integer, intent(in) :: indx
            real(8), intent(in) :: depth
            real(8), pointer    ::  rbot
        !%-----------------------------------------------------------------------------
            rbot => elemSGR(indx, esgr_Parabolic_Radius)
        !%-----------------------------------------------------------------------------
        outvalue = (fourR / threeR) * rbot * depth *  sqrt(depth)

    end function parabolic_area_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function parabolic_topwidth_from_depth_singular (indx, depth) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth for a parabolic cross section of a single element
        !%-----------------------------------------------------------------------------
            integer, intent(in) :: indx 
            real(8), intent(in) :: depth
            real(8), pointer    ::  rbot
        !%-----------------------------------------------------------------------------
            rbot => elemSGR(indx, esgr_Parabolic_Radius)
        !%-----------------------------------------------------------------------------

        outvalue = twoR * rbot * sqrt(depth)

    end function parabolic_topwidth_from_depth_singular
!%
!%==========================================================================  
!%==========================================================================
!%
    real(8) function parabolic_perimeter_from_depth_singular (indx, depth) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes wetted perimeter from known depth for a parabolic cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer    ::  rbot
        real(8) :: xx, tt
        !%-----------------------------------------------------------------------------
            rbot => elemSGR(indx, esgr_Parabolic_Radius)
        !%-----------------------------------------------------------------------------
        xx = twoR * sqrt(depth) / rbot
        tt = sqrt(oneR + xx * xx)

        outvalue = onehalfR * rbot * rbot * (xx * tt + log(xx + tt))

    end function parabolic_perimeter_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
    ! real(8) function parabolic_hyddepth_from_depth_singular (indx,depth) result (outvalue)
    !     !%  
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Computes hydraulic depth from known depth for parabolic cross section of 
    !     !% a single element
    !     !%-----------------------------------------------------------------------------   
    !     integer, intent(in) :: indx   
    !     real(8), intent(in) :: depth  
    !     !%-----------------------------------------------------------------------------  

    !     outvalue = (twoR/threeR) * depth

    ! end function parabolic_hyddepth_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function parabolic_hydradius_from_depth_singular &
        (indx, depth) result (outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic radius from known depth for a parabolic cross section of
        !% a single element 
        !%------------------------------------------------------------------
            integer, intent(in) :: indx
            real(8), intent(in) :: depth
            real(8) :: xx, tt, area, perimeter
            real(8), pointer :: rbot
            !real(8), pointer :: breadth, fulldepth, rbot
        !%------------------------------------------------------------------
        !breadth   => elemSGR(indx,esgr_parabolic_Breadth)
        !fulldepth => elemR(indx, er_FullDepth)
        rbot      => elemSGR(indx, esgr_Parabolic_Radius)
        !%------------------------------------------------------------------
        xx = (twoR * sqrt(depth) / rbot)
        tt = sqrt(oneR + (xx * xx))
        area = (fourR / threeR) * rbot * depth *  sqrt(depth)

        perimeter = onehalfR * rbot * rbot * (xx * tt + log(xx + tt))
        
        outvalue = area / perimeter

    end function parabolic_hydradius_from_depth_singular
!%      
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module parabolic_channel