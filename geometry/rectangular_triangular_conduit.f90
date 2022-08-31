module rectangular_triangular_conduit

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% rectangular_triangular channel geometry
    !%

    private

    public :: rectangular_triangular_depth_from_volume
    public :: rectangular_triangular_area_from_depth
    public :: rectangular_triangular_area_from_depth_singular
    public :: rectangular_triangular_topwidth_from_depth
    public :: rectangular_triangular_topwidth_from_depth_singular 
    public :: rectangular_triangular_perimeter_from_depth
    public :: rectangular_triangular_perimeter_from_depth_singular
    public :: rectangular_triangular_hyddepth_from_depth
    public :: rectangular_triangular_hyddepth_from_depth_singular
    public :: rectangular_triangular_hydradius_from_depth_singular

    contains

!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine rectangular_triangular_depth_from_volume (elemPGx, Npack, thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels (or non-surcharged rectangular_triangular conduits)
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !% NOTE: this does NOT limit the depth by surcharge height at this point
        !% This will be done after the head is computed.
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:), bottomDepth(:), bottomArea(:), volume(:), length(:), breadth(:), bottomSlope(:)
        !%-----------------------------------------------------------------------------
        thisP       => elemPGx(1:Npack,thisCol) 
        depth       => elemR(:,er_Depth)
        bottomDepth => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth) 
        bottomArea  => elemSGR(:,esgr_Rectangular_Triangular_BottomArea)
        breadth     => elemSGR(:,esgr_Rectangular_Triangular_TopBreadth)
        bottomSlope => elemSGR(:,esgr_Rectangular_Triangular_BottomSlope)
        volume      => elemR(:,er_Volume)
        length      => elemR(:,er_Length)
        !%-----------------------------------------------------------------------------  

        where(volume(thisP) <= bottomArea(thisP) * length(thisP))
            depth(thisP) = sqrt(volume(thisP) / (length(thisP) * bottomSlope(thisP)))

        elsewhere(volume(thisP) > bottomarea(thisP)*length(thisP))
            depth(thisP) = bottomDepth(thisP) + ((volume(thisP) / length(thisP)) - bottomarea(thisP)) / breadth(thisP)
        end where

        if (setting%Debug%File%geometry) &
                print *, 'depth =   ' , depth(thisP)
                
    end subroutine rectangular_triangular_depth_from_volume
!%
!%==========================================================================
!%==========================================================================
!%
    elemental real(8) function rectangular_triangular_area_from_depth (indx) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for rectangular cross section
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx  ! may be a packed array of indexes
        !%-----------------------------------------------------------------------------
        
        if (elemR(indx,er_Depth) <= elemSGR(indx,esgr_Rectangular_Triangular_BottomDepth)) then
            outvalue = elemR(indx,er_Depth) * elemR(indx,er_Depth) * elemSGR(indx,esgr_Rectangular_Triangular_BottomSlope)
        else
            outvalue =  elemSGR(indx,esgr_Rectangular_Triangular_BottomArea) &    !triangular section
                        + ((elemR(indx,er_Depth)  - elemSGR(indx,esgr_Rectangular_Triangular_BottomDepth)) * &
                        elemSGR(indx,esgr_rectangular_Triangular_TopBreadth))       !rectangular section
        endif

    end function rectangular_triangular_area_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function rectangular_triangular_area_from_depth_singular (indx, depth) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for rectangular_triangular cross section of a single element
        !% The input indx is the row index in full data 2D array.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer :: bottomDepth(:), bottomSlope(:), bottomArea(:), breadth(:)
        !%-----------------------------------------------------------------------------
        bottomDepth => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth) 
        bottomSlope => elemSGR(:,esgr_Rectangular_Triangular_BottomSlope)
        bottomArea  => elemSGR(:,esgr_Rectangular_Triangular_BottomArea)
        breadth     => elemSGR(:,esgr_Rectangular_Triangular_TopBreadth)
        !%-----------------------------------------------------------------------------
        
        if(depth <= bottomDepth(indx)) then
            outvalue = depth * depth * bottomSlope(indx)
        else
            outvalue = bottomArea(indx) + (depth - bottomDepth(indx)) * breadth(indx)    
        endif

        if (setting%Debug%File%geometry) &
            print *, 'area = ' , outvalue

    end function rectangular_triangular_area_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_triangular_topwidth_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a rectangular_triangular channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: breadth(:), topwidth(:), bottomSlope(:), depth(:), bottomDepth(:)
        !%-----------------------------------------------------------------------------
        thisP       => elemPGx(1:Npack,thisCol) 
        topwidth    => elemR(:,er_Topwidth)
        depth       => elemR(:,er_Depth)
        bottomSlope => elemSGR(:,esgr_Rectangular_Triangular_BottomSlope)
        bottomDepth => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth) 
        breadth     => elemSGR(:,esgr_Rectangular_Triangular_TopBreadth)
        !%-----------------------------------------------------------------------------

        where(depth(thisP) <= bottomDepth(thisP))
            topwidth(thisP) = twoR * bottomSlope(thisP) * depth(thisP)

        else where(depth(thisP) > bottomDepth(thisP))
            topwidth(thisP) = breadth(thisP)
        endwhere

        if (setting%Debug%File%geometry) &
            print *, 'topwidth = ' , topwidth(thisP)

    end subroutine rectangular_triangular_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function rectangular_triangular_topwidth_from_depth_singular (indx, depth) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth for a rectangular_triangular cross section of a single element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx 
        real(8), intent(in) :: depth
        real(8), pointer :: bottomDepth(:), bottomSlope(:), breadth(:)
        !%-----------------------------------------------------------------------------
        bottomSlope => elemSGR(:,esgr_Rectangular_Triangular_BottomSlope)
        bottomDepth => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth) 
        breadth     => elemSGR(:,esgr_Rectangular_Triangular_TopBreadth)
        !%-----------------------------------------------------------------------------
         
        if(depth <= bottomDepth(indx)) then
            outvalue = twoR * bottomSlope(indx) * depth
        else
            outvalue = breadth(indx)
        endif

        if (setting%Debug%File%geometry) &
            print *, 'topwidth = ' , outvalue
    end function rectangular_triangular_topwidth_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_triangular_perimeter_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a rectangular_triangular channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: breadth(:), depth(:), bottomSlope(:), bottomDepth(:), perimeter(:)
        !%-----------------------------------------------------------------------------
        thisP       => elemPGx(1:Npack,thisCol) 
        breadth     => elemSGR(:,esgr_Rectangular_Triangular_TopBreadth)
        depth       => elemR(:,er_Depth)
        bottomSlope => elemSGR(:,esgr_Rectangular_Triangular_BottomSlope)
        bottomDepth => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth)  
        perimeter   => elemR(:,er_Perimeter)
        !%-----------------------------------------------------------------------------

        where(depth(thisP) <= bottomdepth(thisP))
            perimeter(thisP) = twoR * depth(thisP) * sqrt(oneR + bottomSlope(thisP) ** twoR)

        elsewhere
            perimeter(thisP) = twoR * bottomdepth(thisP) * sqrt(oneR + bottomSlope(thisP) ** twoR)  & !triangular section
                                 + twoR * (depth(thisP) - bottomDepth(thisP))                         !rectangular section
        end where

        if (setting%Debug%File%geometry) &
                print *, 'perimeter = ' , perimeter(thisP)

    end subroutine rectangular_triangular_perimeter_from_depth
!%    
!%==========================================================================    
!%==========================================================================
!%
    real(8) function rectangular_triangular_perimeter_from_depth_singular (indx, depth) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes wetted perimeter from known depth for a rectangular_triangular cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer :: bottomDepth(:), bottomSlope(:), breadth(:)
        !%-----------------------------------------------------------------------------
        breadth     => elemSGR(:,esgr_Rectangular_Triangular_TopBreadth)
        bottomSlope => elemSGR(:,esgr_Rectangular_Triangular_BottomSlope)
        bottomDepth => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth) 
        !%-----------------------------------------------------------------------------
        
        if(depth <= bottomDepth(indx)) then
            outvalue = twoR * depth * sqrt(oneR + bottomSlope(indx) ** twoR)
        else
            outvalue = twoR * bottomDepth(indx) * sqrt(oneR + bottomSlope(indx) ** twoR) & !triangular section
                        + twoR * (depth - bottomDepth(indx))                               !rectangular section

        endif

        if (setting%Debug%File%geometry) &
            print *, 'perimeter = ' , outvalue

    end function rectangular_triangular_perimeter_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_triangular_hyddepth_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the hydraulic (average) depth from a known depth in a rectangular_triangular channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: hyddepth(:), depth(:),  bottomdepth(:)
        !%-----------------------------------------------------------------------------
        thisP       => elemPGx(1:Npack,thisCol) 
        depth       => elemR(:,er_Depth)
        hyddepth    => elemR(:,er_HydDepth)
        bottomdepth => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth) 
        !%-----------------------------------------------------------------------------

        where(depth(thisP) <= bottomdepth(thisP))
            hyddepth(thisP) = depth(thisP) / twoR

        else where(depth(thisP) > bottomdepth(thisP))
            hyddepth(thisP) = (bottomdepth(thisP) / twoR) &         !triangular section
                            + (depth(thisP) - bottomdepth(thisP))   !rectangular Section
        endwhere

        if (setting%Debug%File%geometry) &
            print *, 'hydepth = ' , hyddepth(thisP)

    end subroutine rectangular_triangular_hyddepth_from_depth
!%    
!%==========================================================================  
!%==========================================================================
!%
    real(8) function rectangular_triangular_hyddepth_from_depth_singular (indx, depth) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic depth from known depth for rectangular_triangular cross section of 
        !% a single element
        !%-----------------------------------------------------------------------------   
        integer, intent(in) :: indx     
        real(8), intent(in) :: depth
        real(8), pointer :: bottomdepth(:)
        !%-----------------------------------------------------------------------------
        bottomdepth => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth) 
        !%-----------------------------------------------------------------------------  

        if(depth <= bottomdepth(indx)) then
            outvalue = depth / twoR
        else
            outvalue = (bottomdepth(indx) / twoR) &         !triangular section
                     + (depth - bottomdepth(indx))          !rectangular section
        endif

        if (setting%Debug%File%geometry) &
            print *, 'hyddepth = ' , outvalue

    end function rectangular_triangular_hyddepth_from_depth_singular 
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function rectangular_triangular_hydradius_from_depth_singular (indx, depth) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic radius from known depth for a rectangular_triangular cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer :: bottomdepth(:), breadth(:), bottomSlope(:)
        !%-----------------------------------------------------------------------------
        bottomSlope => elemSGR(:,esgr_Rectangular_Triangular_BottomSlope)
        breadth     => elemSGR(:,esgr_Rectangular_Triangular_TopBreadth)
        bottomdepth => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth) 
        !%-----------------------------------------------------------------------------
        
        if(depth <= bottomdepth(indx)) then
            outvalue = (bottomSlope(indx) * depth) / (twoR * sqrt(oneR + (bottomSlope(indx) ** twoR)))
        else
            outvalue = ((bottomSlope(indx) * bottomdepth(indx)) / (twoR * sqrt(oneR + (bottomSlope(indx) ** twoR)))) &          !triangular section
                     + (((depth - bottomdepth(indx)) * breadth(indx)) / (twoR * (depth - bottomdepth(indx))))                   !rectangular section
        endif

        if (setting%Debug%File%geometry) &
            print *, 'hydradius = ' , outvalue

    end function rectangular_triangular_hydradius_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
end module rectangular_triangular_conduit