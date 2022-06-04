module rectangular_triangular_channel

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
        real(8), pointer :: depth(:), bottomdepth(:), volume(:), length(:), breadth(:), sideslope(:)
        real(8) :: bottomvolume
        !%-----------------------------------------------------------------------------
        thisP       => elemPGx(1:Npack,thisCol) 
        depth       => elemR(:,er_Depth)
        bottomdepth => elemR(:,er_BottomDepth)
        volume      => elemR(:,er_Volume)
        length      => elemR(:,er_Length)
        breadth     => elemSGR(:,esgr_rectangular_Triangular_TopBreadth)
        sideslope   => elemSGR(:,esgr_rectangular_Triangular_Slope)
        !%-----------------------------------------------------------------------------  

        bottomvolume = (sideslope(thisP) * (bottomdepth(thisP) ** twoR)) * length(thisP)

        if(volume(thisP) <= bottomvolume) then
            depth(thisP) = sqrt((volume(thisP)/length(thisP)) / sideslope(thisP))
        else
            depth(thisP) = (volume(thisP) - bottomvolume) / (length(thisP) * breadth(thisP)) + bottomdepth(thisP)
        endif

    end subroutine rectangular_triangular_depth_from_volume
    !%  
!%==========================================================================
!%==========================================================================
!%
    elemental real(8) function rectangular_triangular_area_from_depth (indx) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for rectangular_triangular cross section
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx  ! may be a packed array of indexes
        real(8), pointer :: depth(:), bottomdepth(:), breadth(:), sideslope(:)
        !%-----------------------------------------------------------------------------

        if(elemR(indx,er_Depth) <= elemR(indx,er_BottomDepth)) then
            outvalue = elemR(indx,er_Depth) * elemR(indx,er_Depth) * elemSGR(indx,esgr_rectangular_Triangular_Slope)
        else
            outvalue = (elemR(indx,er_Depth) * elemR(indx,er_Depth) * elemSGR(indx,esgr_rectangular_Triangular_Slope)) & !triangular section
                        + ((elemR(indx,er_Depth) - elemR(indx,er_BottomDepth) * elemSGR(indx,esgr_rectangular_Triangular_TopBreadth))) !rectangular section
        end if

    end function rectangular_triangular_area_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function rectangular_triangular_area_from_depth_singular (indx) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for rectangular_triangular cross section of a single element
        !% The input indx is the row index in full data 2D array.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), pointer :: depth(:), bottomdepth(:), breadth(:), sideslope(:)
        !%-----------------------------------------------------------------------------
        depth       => elemR(:,er_Depth)
        bottomdepth => elemR(:,er_BottomDepth)
        sideslope   => elemSGR(:,esgr_rectangular_Triangular_Slope)
        breadth     => elemSGR(:,esgr_rectangular_Triangular_TopBreadth)
        !%-----------------------------------------------------------------------------
        
        if(depth(indx) <= bottomdepth(indx)) then
            outvalue = depth(indx) * depth(indx) * sideslope(indx)
        else
            outvalue = (bottomdepth(indx) * bottomdepth(indx) * sideslope(indx)) &  !triangular section
                        + ((depth(indx) - bottomdepth(indx)) * breadth(indx))       !rectangular section
        endif

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
        real(8), pointer :: breadth(:), topwidth(:), sideslope(:), depth(:), bottomdepth(:)
        !%-----------------------------------------------------------------------------
        thisP       => elemPGx(1:Npack,thisCol) 
        topwidth    => elemR(:,er_Topwidth)
        sideslope   => elemSGR(:,esgr_rectangular_Triangular_Slope)
        depth       => elemR(:,er_Depth)
        bottomdepth => elemR(:,er_BottomDepth)
        breadth     => elemSGR(:,esgr_rectangular_Triangular_TopBreadth)
        !%-----------------------------------------------------------------------------

        if(elemR(thisP,er_Depth) <= elemR(thisP,er_BottomDepth)) then
            topwidth(thisP) = twoR * sideslope(thisP) * depth(thisP)
        else
            topwidth(thisP) = breadth(thisP)
        endif

    end subroutine rectangular_triangular_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function rectangular_triangular_topwidth_from_depth_singular (indx) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth for a rectangular_triangular cross section of a single element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx 
        real(8), pointer :: depth(:), bottomdepth(:), sideslope(:), breadth(:)
        !%-----------------------------------------------------------------------------
        depth       => elemR(:,er_Depth)
        sideslope   => elemSGR(:,esgr_rectangular_Triangular_Slope)
        bottomdepth => elemR(:,er_BottomDepth)
        breadth     => elemSGR(:,esgr_rectangular_Triangular_TopBreadth)
        !%-----------------------------------------------------------------------------
         
        if(depth(indx) <= bottomdepth(indx)) then
            outvalue = twoR * sideslope(indx) * depth(indx)
        else
            outvalue = breadth(indx)
        endif

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
        real(8), pointer :: breadth(:), depth(:), bottomdepth(:), perimeter(:)
        !%-----------------------------------------------------------------------------
        thisP       => elemPGx(1:Npack,thisCol) 
        breadth     => elemSGR(:,esgr_rectangular_Triangular_TopBreadth)
        depth       => elemR(:,er_Depth)
        bottomdepth => elemR(:,er_BottomDepth)
        perimeter   => elemR(:,er_Perimeter)
        !%-----------------------------------------------------------------------------

        if(depth(thisP) <= bottomdepth(thisP)) then
            perimeter(thisP) = twoR * sqrt(((breadth(thisP) ** twoR) / fourR) + (depth(thisP) ** twoR))
        else
            perimeter(thisP) = (twoR * sqrt(((breadth(thisP) ** twoR) / fourR) + (bottomdepth(thisP) ** twoR))) & !triangular section
                                + (twoR * ((depth(thisP) - bottomdepth(thisP)) + breadth(thisP)))           !rectangular section
        endif

    end subroutine rectangular_triangular_perimeter_from_depth
!%    
!%==========================================================================    
!%==========================================================================
!%
    real(8) function rectangular_triangular_perimeter_from_depth_singular (indx) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes wetted perimeter from known depth for a rectangular_triangular cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), pointer :: depth(:), bottomdepth(:), breadth(:)
        !%-----------------------------------------------------------------------------
        breadth     => elemSGR(:,esgr_rectangular_Triangular_TopBreadth)
        depth       => elemR(:,er_Depth)
        bottomdepth => elemR(:,er_BottomDepth)
        !%-----------------------------------------------------------------------------
        
        if(depth(indx) <= bottomdepth(indx)) then
            outvalue = twoR * sqrt(((breadth(indx) ** twoR) / fourR) + (depth(indx) ** twoR))
        else
            outvalue = (twoR * sqrt(((breadth(indx) ** twoR) / fourR) + (bottomdepth(indx) ** twoR))) & !triangular section
                        + (twoR * ((depth(indx) - bottomdepth(indx)) + breadth(indx)))            !rectangular section

        endif

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
        bottomdepth => elemR(:,er_BottomDepth)
        !%-----------------------------------------------------------------------------

        if (depth(thisP) <= bottomdepth(thisP)) then
            hyddepth(thisP) = depth(thisP) / twoR
        else
            hyddepth(thisP) = (bottomdepth(thisP) / twoR) &         !triangular section
                            + (depth(thisP) - bottomdepth(thisP))   !rectangular Section
        endif

    end subroutine rectangular_triangular_hyddepth_from_depth
!%    
!%==========================================================================  
!%==========================================================================
!%
    real(8) function rectangular_triangular_hyddepth_from_depth_singular (indx) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic depth from known depth for rectangular_triangular cross section of 
        !% a single element
        !%-----------------------------------------------------------------------------   
        integer, intent(in) :: indx     
        real(8), pointer :: depth(:), bottomdepth(:)
        !%-----------------------------------------------------------------------------
        depth       => elemR(:,er_Depth)
        bottomdepth => elemR(:,er_BottomDepth)
        !%-----------------------------------------------------------------------------  

        if(depth(indx) <= bottomdepth(indx)) then
            outvalue = depth(indx) / twoR
        else
            outvalue = (bottomdepth(indx) / twoR) &         !triangular section
                     + (depth(indx) - bottomdepth(indx))    !rectangular section
        endif

    end function rectangular_triangular_hyddepth_from_depth_singular 
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function rectangular_triangular_hydradius_from_depth_singular (indx) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic radius from known depth for a rectangular_triangular cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), pointer :: depth(:), bottomdepth(:), breadth(:), sideslope(:)
        !%-----------------------------------------------------------------------------
        depth       => elemR(:,er_Depth)
        sideslope   => elemSGR(:,esgr_rectangular_Triangular_Slope)
        breadth     => elemSGR(:,esgr_rectangular_Triangular_TopBreadth)
        bottomdepth => elemR(:,er_BottomDepth)
        !%-----------------------------------------------------------------------------
        
        if(depth(indx) <= bottomdepth(indx)) then
            outvalue = (sideslope(indx) * depth(indx)) / (twoR * sqrt(oneR + (sideslope(indx) ** twoR)))
        else
            outvalue = ((sideslope(indx) * bottomdepth(indx)) / (twoR * sqrt(oneR + (sideslope(indx) ** twoR)))) &                          !triangular section
                     + (((depth(indx) - bottomdepth(indx)) * breadth(indx)) / (breadth(indx) + (twoR * (depth(indx) - bottomdepth(indx))))) !rectangular section
        endif

    end function rectangular_triangular_hydradius_from_depth_singular
    !%    
!%==========================================================================

end module rectangular_triangular_channel