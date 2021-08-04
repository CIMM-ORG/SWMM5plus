module trapezoidal_channel

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Trapezoidal channel geometry
    !%

    private

    public :: trapezoidal_depth_from_volume
    public :: trapezoidal_area_from_depth_singular
    public :: trapezoidal_topwidth_from_depth
    public :: trapezoidal_topwidth_from_depth_singular 
    public :: trapezoidal_perimeter_from_depth
    public :: trapezoidal_perimeter_from_depth_singular
    public :: trapezoidal_hyddepth_from_depth
    public :: trapezoidal_hyddepth_from_depth_singular
    public :: trapezoidal_hydradius_from_depth_singular


    contains
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    !%
    subroutine trapezoidal_depth_from_volume (elemPGx, Npack, thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels (or non-surcharged trapezoidal conduits)
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !% NOTE: this does NOT limit the depth by surcharge height at this point
        !% This will be done after the head is computed.
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:), volume(:), length(:), breadth(:)
        real(8), pointer :: lslope(:), rslope(:)
        !%-----------------------------------------------------------------------------
        thisP   => elemPGx(1:Npack,thisCol) 
        depth   => elemR(:,er_Depth)
        volume  => elemR(:,er_Volume)
        length  => elemR(:,er_Length)
        breadth => elemSGR(:,eSGR_Trapezoidal_Breadth)
        lslope  => elemSGR(:,eSGR_Trapezoidal_LeftSlope)
        rslope  => elemSGR(:,eSGR_Trapezoidal_RightSlope)
        !%-----------------------------------------------------------------------------  

        depth(thisP)       = - onehalfR * (breadth(thisP)/(onehalfR*(lslope(thisP) + rslope(thisP))) &
                - sqrt((breadth(thisP)/(onehalfR*(lslope(thisP) + rslope(thisP)))) ** twoR &
                + fourR * volume(thisP)/(onehalfR*length(thisP)*(lslope(thisP) + rslope(thisP)))))

    end subroutine trapezoidal_depth_from_volume
    !%  
    !%==========================================================================
    !%==========================================================================
    !%
    real(8) function trapezoidal_area_from_depth_singular (indx) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for trapezoidal cross section of a single element
        !% The input indx is the row index in full data 2D array.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), pointer :: depth(:), breadth(:), lslope(:), rslope(:)
        !%-----------------------------------------------------------------------------
        depth   => elemR(:,er_Depth)
        breadth => elemSGR(:,eSGR_Trapezoidal_Breadth)
        lslope  => elemSGR(:,eSGR_Trapezoidal_LeftSlope)
        rslope  => elemSGR(:,eSGR_Trapezoidal_RightSlope)
        !%-----------------------------------------------------------------------------
        outvalue = (breadth(indx) + onehalfR * (lslope(indx) + rslope(indx)) * depth(indx)) * depth(indx)

    end function trapezoidal_area_from_depth_singular
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine trapezoidal_topwidth_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a rectangular channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: topwidth(:), depth(:), breadth(:), lslope(:), rslope(:)
        !%-----------------------------------------------------------------------------
        thisP    => elemPGx(1:Npack,thisCol) 
        topwidth => elemR(:,er_Topwidth)
        depth    => elemR(:,er_Depth)
        breadth  => elemSGR(:,eSGR_Trapezoidal_Breadth)
        lslope   => elemSGR(:,eSGR_Trapezoidal_LeftSlope)
        rslope   => elemSGR(:,eSGR_Trapezoidal_RightSlope)
        !%-----------------------------------------------------------------------------

        topwidth(thisP) = breadth(thisP) + depth(thisP) * (lslope(thisP) + rslope(thisP))

    end subroutine trapezoidal_topwidth_from_depth
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    real(8) function trapezoidal_topwidth_from_depth_singular (indx) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth for a trapezoidal cross section of a single element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx 
        real(8), pointer :: depth(:), breadth(:), lslope(:), rslope(:)
        !%-----------------------------------------------------------------------------
        depth   => elemR(:,er_Depth)
        breadth => elemSGR(:,eSGR_Trapezoidal_Breadth)
        lslope  => elemSGR(:,eSGR_Trapezoidal_LeftSlope)
        rslope  => elemSGR(:,eSGR_Trapezoidal_RightSlope)
        !%-----------------------------------------------------------------------------
        outvalue = breadth(indx) + depth(indx) * (lslope(indx) + rslope(indx))

    end function trapezoidal_topwidth_from_depth_singular
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine trapezoidal_perimeter_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a trapezoidal channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: perimeter(:), depth(:), breadth(:), lslope(:), rslope(:)
        !%-----------------------------------------------------------------------------
        thisP     => elemPGx(1:Npack,thisCol) 
        perimeter => elemR(:,er_Perimeter)
        depth     => elemR(:,er_Depth)
        breadth   => elemSGR(:,eSGR_Trapezoidal_Breadth)
        lslope    => elemSGR(:,eSGR_Trapezoidal_LeftSlope)
        rslope    => elemSGR(:,eSGR_Trapezoidal_RightSlope)
        !%-----------------------------------------------------------------------------

        perimeter(thisP) = breadth(thisP) + depth(thisP) * (sqrt(oneR + lslope(thisP)**twoR) &
                        + sqrt(oneR + rslope(thisP)**twoR))

    end subroutine trapezoidal_perimeter_from_depth
    !%    
    !%==========================================================================    
    !%==========================================================================
    !%
    real(8) function trapezoidal_perimeter_from_depth_singular (indx) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes wetted perimeter from known depth for a trapezoidal cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), pointer :: depth(:), breadth(:), lslope(:), rslope(:)
        !%-----------------------------------------------------------------------------
        depth   => elemR(:,er_Depth)
        breadth => elemSGR(:,eSGR_Trapezoidal_Breadth)
        lslope  => elemSGR(:,eSGR_Trapezoidal_LeftSlope)
        rslope  => elemSGR(:,eSGR_Trapezoidal_RightSlope)
        !%----------------------------------------------------------------------------- 
        
        outvalue =  breadth(indx) + depth(indx) * (sqrt(oneR + lslope(indx)**twoR) + &
                    sqrt(oneR + rslope(indx)**twoR))

    end function trapezoidal_perimeter_from_depth_singular
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine trapezoidal_hyddepth_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the hydraulic (average) depth from a known depth in a rectangular channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: hyddepth(:), depth(:), breadth(:), lslope(:), rslope(:)
        !%-----------------------------------------------------------------------------
        thisP     => elemPGx(1:Npack,thisCol) 
        hyddepth  => elemR(:,er_HydDepth)
        depth   => elemR(:,er_Depth)
        breadth => elemSGR(:,eSGR_Trapezoidal_Breadth)
        lslope  => elemSGR(:,eSGR_Trapezoidal_LeftSlope)
        rslope  => elemSGR(:,eSGR_Trapezoidal_RightSlope)
        !%----------------------------------------------------------------------------- 

        hyddepth(thisP) = ((breadth(thisP) + onehalfR * (lslope(thisP) + rslope(thisP)) * &
                    depth(thisP)) * depth(thisP)) / (breadth(thisP) + depth(thisP) * &
                    (lslope(thisP) + rslope(thisP)))

    end subroutine trapezoidal_hyddepth_from_depth
    !%    
    !%==========================================================================  
    !%==========================================================================
    !%
    real(8) function trapezoidal_hyddepth_from_depth_singular (indx) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic depth from known depth for trapezoidal cross section of 
        !% a single element
        !%-----------------------------------------------------------------------------   
        integer, intent(in) :: indx 
        real(8), pointer :: depth(:), breadth(:), lslope(:), rslope(:)
        !%-----------------------------------------------------------------------------
        depth   => elemR(:,er_Depth)
        breadth => elemSGR(:,eSGR_Trapezoidal_Breadth)
        lslope  => elemSGR(:,eSGR_Trapezoidal_LeftSlope)
        rslope  => elemSGR(:,eSGR_Trapezoidal_RightSlope)
        !%-----------------------------------------------------------------------------     

        outvalue = ((breadth(indx) + onehalfR * (lslope(indx) + rslope(indx)) * &
                    depth(indx)) * depth(indx)) / (breadth(indx) + depth(indx) * &
                    (lslope(indx) + rslope(indx)))

    end function trapezoidal_hyddepth_from_depth_singular
    !%    
    !%==========================================================================

    !%==========================================================================
    !%
    real(8) function trapezoidal_hydradius_from_depth_singular (indx) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic radius from known depth for a trapezoidal cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), pointer :: depth(:), breadth(:), lslope(:), rslope(:)
        !%-----------------------------------------------------------------------------
        depth   => elemR(:,er_Depth)
        breadth => elemSGR(:,eSGR_Trapezoidal_Breadth)
        lslope  => elemSGR(:,eSGR_Trapezoidal_LeftSlope)
        rslope  => elemSGR(:,eSGR_Trapezoidal_RightSlope)
        !%----------------------------------------------------------------------------- 
        
        outvalue = ((breadth(indx) + onehalfR * (lslope(indx) + rslope(indx)) *  &
                    depth(indx)) * depth(indx)) / (breadth(indx) + depth(indx) * &
                    (sqrt(oneR + lslope(indx)**twoR) + sqrt(oneR + rslope(indx)**twoR)))

    end function trapezoidal_hydradius_from_depth_singular
    !%    
    !%==========================================================================

    !%
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%  
    !%
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%  
    !%
    !%    
    !%==========================================================================
    !%==========================================================================
    !%
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%  
end module trapezoidal_channel