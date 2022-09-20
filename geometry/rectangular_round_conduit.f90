module rectangular_round_conduit

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% rectangular_round channel geometry
    !%

    private

    public :: get_Y
    public :: get_theta_of_alpha
    public :: rectangular_round_depth_from_volume
    public :: rectangular_round_area_from_depth
    public :: rectangular_round_area_from_depth_singular
    public :: rectangular_round_topwidth_from_depth
    public :: rectangular_round_topwidth_from_depth_singular 
    public :: rectangular_round_perimeter_from_depth
    public :: rectangular_round_perimeter_from_depth_singular
    public :: rectangular_round_hyddepth_from_depth
    public :: rectangular_round_hyddepth_from_depth_singular
    public :: rectangular_round_hydradius_from_depth_singular

    contains

!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    real(8) subroutine get_Y (alpha) result (outvalue)
        real(8), intent(in) :: alpha
        real(8) :: theta, outvalue

        if (alpha >= oneR) then 
            outvalue = oneR
        else if (alpha <= zeroR) then
            outvalue = zeroR
        else if (alpha <= 1.0e-5) then
            theta = 37.6911 ** onethirdR
            outvalue = theta * theta / sixteenR
        else 
            theta = get_theta_of_alpha(alpha)
            outvalue = (oneR - cos(theta/twoR)) / twoR
    end subroutine get_Y

    subroutine get_theta_of_alpha(alpha) result (outvalue)
        real(8), intent(in) :: alpha
        integer :: k
        real(8) :: theta, theta1, ap, d

        if (alpha > 0.04d) then
            theta = 1.2d + 5.08d * (alpha - 0.04d) / 0.96d
        else 
            theta = 0.031715 - 12.79384 * alpha + 8.28479 * sqrt(alpha);
        theta1 = theta
        ap = (twoR * pi) * alpha
        

    end subroutine get_theta_of_alpha


    subroutine rectangular_round_depth_from_volume (elemPGx, Npack, thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels (or non-surcharged rectangular_round conduits)
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !% NOTE: this does NOT limit the depth by surcharge height at this point
        !% This will be done after the head is computed.
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:), bottomDepth(:), bottomArea(:), volume(:), length(:), breadth(:), bottomSlope(:)
        real(8) :: alpha
        !%-----------------------------------------------------------------------------
        thisP       => elemPGx(1:Npack,thisCol) 
        depth       => elemR(:,er_Depth)
        bottomDepth => elemSGR(:,esgr_Rectangular_Round_BottomDepth) 
        bottomArea  => elemSGR(:,esgr_Rectangular_Round_BottomArea)
        breadth     => elemSGR(:,esgr_Rectangular_Round_TopBreadth)
        bottomRadius => elemSGR(:,esgr_Rectangular_Round_BottomRadius)
        volume      => elemR(:,er_Volume)
        length      => elemR(:,er_Length)
        !%-----------------------------------------------------------------------------  

        where(volume(thisP) <= bottomArea(thisP) * length(thisP))
            alpha = volume(thisP) / length(thisP) / (pi * bottomRadius(thisP) * bottomRadius(thisP) )
            if (alpha < 0.04) then
                depth(thisP) = twoR * bottomRadius(thisP) * get_Y(alpha)
            else
                depth(thisP) = () * 




        elsewhere(volume(thisP) > bottomarea(thisP)*length(thisP))
            depth(thisP) = bottomDepth(thisP) + ((volume(thisP) / length(thisP)) - bottomarea(thisP)) / breadth(thisP)
        end where
        

        if (setting%Debug%File%geometry) &
                print *, 'depth =   ' , depth(thisP)
                
    end subroutine rectangular_round_depth_from_volume
!%
!%==========================================================================
!%==========================================================================
!%
    elemental real(8) function rectangular_round_area_from_depth (indx) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for rectangular cross section
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx  ! may be a packed array of indexes
        !%-----------------------------------------------------------------------------
        
        if (elemR(indx,er_Depth) <= elemSGR(indx,esgr_Rectangular_Round_BottomDepth)) then
            outvalue = elemR(indx,er_Depth) * elemR(indx,er_Depth) * elemSGR(indx,esgr_Rectangular_Round_BottomSlope)
        else
            outvalue =  elemSGR(indx,esgr_Rectangular_Round_BottomArea) &    !round section
                        + ((elemR(indx,er_Depth)  - elemSGR(indx,esgr_Rectangular_Round_BottomDepth)) * &
                        elemSGR(indx,esgr_rectangular_Round_TopBreadth))       !rectangular section
        endif

    end function rectangular_round_area_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function rectangular_round_area_from_depth_singular (indx, depth) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for rectangular_round cross section of a single element
        !% The input indx is the row index in full data 2D array.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer :: bottomDepth(:), bottomSlope(:), bottomArea(:), breadth(:)
        !%-----------------------------------------------------------------------------
        bottomDepth => elemSGR(:,esgr_Rectangular_Round_BottomDepth) 
        bottomSlope => elemSGR(:,esgr_Rectangular_Round_BottomSlope)
        bottomArea  => elemSGR(:,esgr_Rectangular_Round_BottomArea)
        breadth     => elemSGR(:,esgr_Rectangular_Round_TopBreadth)
        !%-----------------------------------------------------------------------------
        
        if(depth <= bottomDepth(indx)) then
            outvalue = depth * depth * bottomSlope(indx)
        else
            outvalue = bottomArea(indx) + (depth - bottomDepth(indx)) * breadth(indx)    
        endif

        if (setting%Debug%File%geometry) &
            print *, 'area = ' , outvalue

    end function rectangular_round_area_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_round_topwidth_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a rectangular_round channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: breadth(:), topwidth(:), bottomSlope(:), depth(:), bottomDepth(:)
        !%-----------------------------------------------------------------------------
        thisP       => elemPGx(1:Npack,thisCol) 
        topwidth    => elemR(:,er_Topwidth)
        depth       => elemR(:,er_Depth)
        bottomSlope => elemSGR(:,esgr_Rectangular_Round_BottomSlope)
        bottomDepth => elemSGR(:,esgr_Rectangular_Round_BottomDepth) 
        breadth     => elemSGR(:,esgr_Rectangular_Round_TopBreadth)
        !%-----------------------------------------------------------------------------

        where(depth(thisP) <= bottomDepth(thisP))
            topwidth(thisP) = twoR * bottomSlope(thisP) * depth(thisP)

        else where(depth(thisP) > bottomDepth(thisP))
            topwidth(thisP) = breadth(thisP)
        endwhere

        if (setting%Debug%File%geometry) &
            print *, 'topwidth = ' , topwidth(thisP)

    end subroutine rectangular_round_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function rectangular_round_topwidth_from_depth_singular (indx, depth) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth for a rectangular_round cross section of a single element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx 
        real(8), intent(in) :: depth
        real(8), pointer :: bottomDepth(:), bottomSlope(:), breadth(:)
        !%-----------------------------------------------------------------------------
        bottomSlope => elemSGR(:,esgr_Rectangular_Round_BottomSlope)
        bottomDepth => elemSGR(:,esgr_Rectangular_Round_BottomDepth) 
        breadth     => elemSGR(:,esgr_Rectangular_Round_TopBreadth)
        !%-----------------------------------------------------------------------------
         
        if(depth <= bottomDepth(indx)) then
            outvalue = twoR * bottomSlope(indx) * depth
        else
            outvalue = breadth(indx)
        endif

        if (setting%Debug%File%geometry) &
            print *, 'topwidth = ' , outvalue
    end function rectangular_round_topwidth_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_round_perimeter_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a rectangular_round channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: breadth(:), depth(:), bottomSlope(:), bottomDepth(:), perimeter(:)
        !%-----------------------------------------------------------------------------
        thisP       => elemPGx(1:Npack,thisCol) 
        breadth     => elemSGR(:,esgr_Rectangular_Round_TopBreadth)
        depth       => elemR(:,er_Depth)
        bottomSlope => elemSGR(:,esgr_Rectangular_Round_BottomSlope)
        bottomDepth => elemSGR(:,esgr_Rectangular_Round_BottomDepth)  
        perimeter   => elemR(:,er_Perimeter)
        !%-----------------------------------------------------------------------------

        where(depth(thisP) <= bottomdepth(thisP))
            perimeter(thisP) = twoR * depth(thisP) * sqrt(oneR + bottomSlope(thisP) ** twoR)

        elsewhere
            perimeter(thisP) = twoR * bottomdepth(thisP) * sqrt(oneR + bottomSlope(thisP) ** twoR)  & !round section
                                 + twoR * (depth(thisP) - bottomDepth(thisP))                         !rectangular section
        end where

        if (setting%Debug%File%geometry) &
                print *, 'perimeter = ' , perimeter(thisP)

    end subroutine rectangular_round_perimeter_from_depth
!%    
!%==========================================================================    
!%==========================================================================
!%
    real(8) function rectangular_round_perimeter_from_depth_singular (indx, depth) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes wetted perimeter from known depth for a rectangular_round cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer :: bottomDepth(:), bottomSlope(:), breadth(:)
        !%-----------------------------------------------------------------------------
        breadth     => elemSGR(:,esgr_Rectangular_Round_TopBreadth)
        bottomSlope => elemSGR(:,esgr_Rectangular_Round_BottomSlope)
        bottomDepth => elemSGR(:,esgr_Rectangular_Round_BottomDepth) 
        !%-----------------------------------------------------------------------------
        
        if(depth <= bottomDepth(indx)) then
            outvalue = twoR * depth * sqrt(oneR + bottomSlope(indx) ** twoR)
        else
            outvalue = twoR * bottomDepth(indx) * sqrt(oneR + bottomSlope(indx) ** twoR) & !round section
                        + twoR * (depth - bottomDepth(indx))                               !rectangular section

        endif

        if (setting%Debug%File%geometry) &
            print *, 'perimeter = ' , outvalue

    end function rectangular_round_perimeter_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_round_hyddepth_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the hydraulic (average) depth from a known depth in a rectangular_round channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: hyddepth(:), depth(:),  bottomdepth(:)
        !%-----------------------------------------------------------------------------
        thisP       => elemPGx(1:Npack,thisCol) 
        depth       => elemR(:,er_Depth)
        hyddepth    => elemR(:,er_HydDepth)
        bottomdepth => elemSGR(:,esgr_Rectangular_Round_BottomDepth) 
        !%-----------------------------------------------------------------------------

        where(depth(thisP) <= bottomdepth(thisP))
            hyddepth(thisP) = depth(thisP) / twoR

        else where(depth(thisP) > bottomdepth(thisP))
            hyddepth(thisP) = (bottomdepth(thisP) / twoR) &         !round section
                            + (depth(thisP) - bottomdepth(thisP))   !rectangular Section
        endwhere

        if (setting%Debug%File%geometry) &
            print *, 'hydepth = ' , hyddepth(thisP)

    end subroutine rectangular_round_hyddepth_from_depth
!%    
!%==========================================================================  
!%==========================================================================
!%
    real(8) function rectangular_round_hyddepth_from_depth_singular (indx, depth) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic depth from known depth for rectangular_round cross section of 
        !% a single element
        !%-----------------------------------------------------------------------------   
        integer, intent(in) :: indx     
        real(8), intent(in) :: depth
        real(8), pointer :: bottomdepth(:)
        !%-----------------------------------------------------------------------------
        bottomdepth => elemSGR(:,esgr_Rectangular_Round_BottomDepth) 
        !%-----------------------------------------------------------------------------  

        if(depth <= bottomdepth(indx)) then
            outvalue = depth / twoR
        else
            outvalue = (bottomdepth(indx) / twoR) &         !round section
                     + (depth - bottomdepth(indx))          !rectangular section
        endif

        if (setting%Debug%File%geometry) &
            print *, 'hyddepth = ' , outvalue

    end function rectangular_round_hyddepth_from_depth_singular 
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function rectangular_round_hydradius_from_depth_singular (indx, depth) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic radius from known depth for a rectangular_round cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer :: bottomdepth(:), breadth(:), bottomSlope(:)
        !%-----------------------------------------------------------------------------
        bottomSlope => elemSGR(:,esgr_Rectangular_Round_BottomSlope)
        breadth     => elemSGR(:,esgr_Rectangular_Round_TopBreadth)
        bottomdepth => elemSGR(:,esgr_Rectangular_Round_BottomDepth) 
        !%-----------------------------------------------------------------------------
        
        if(depth <= bottomdepth(indx)) then
            outvalue = (bottomSlope(indx) * depth) / (twoR * sqrt(oneR + (bottomSlope(indx) ** twoR)))
        else
            outvalue = ((bottomSlope(indx) * bottomdepth(indx)) / (twoR * sqrt(oneR + (bottomSlope(indx) ** twoR)))) &          !round section
                     + (((depth - bottomdepth(indx)) * breadth(indx)) / (twoR * (depth - bottomdepth(indx))))                   !rectangular section
        endif

        if (setting%Debug%File%geometry) &
            print *, 'hydradius = ' , outvalue

    end function rectangular_round_hydradius_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
end module rectangular_round_conduit