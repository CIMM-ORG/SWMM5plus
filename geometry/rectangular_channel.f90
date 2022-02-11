module rectangular_channel

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Rectangular channel geometry
    !%

    private

    public :: rectangular_depth_from_volume
    public :: rectangular_area_from_depth
    public :: rectangular_area_from_depth_singular
    public :: rectangular_topwidth_from_depth
    public :: rectangular_topwidth_from_depth_singular 
    public :: rectangular_perimeter_from_depth
    public :: rectangular_perimeter_from_depth_singular
    public :: rectangular_hyddepth_from_depth
    public :: rectangular_hyddepth_from_depth_singular
    public :: rectangular_hydradius_from_depth_singular

    contains

!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine rectangular_depth_from_volume (elemPGx, Npack, thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels (or non-surcharged rectangular conduits)
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !% NOTE: this does NOT limit the depth by surcharge height at this point
        !% This will be done after the head is computed.
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:), volume(:), length(:), breadth(:)
        !%-----------------------------------------------------------------------------
        thisP   => elemPGx(1:Npack,thisCol) 
        depth   => elemR(:,er_Depth)
        volume  => elemR(:,er_Volume)
        length  => elemR(:,er_Length)
        breadth => elemSGR(:,esgr_Rectangular_Breadth)
        !%-----------------------------------------------------------------------------  

        depth(thisP) = volume(thisP) / (length(thisP) * breadth(thisP))

    end subroutine rectangular_depth_from_volume
    !%  
!%==========================================================================
!%==========================================================================
!%
    elemental real(8) function rectangular_area_from_depth (indx) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for rectangular cross section
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx  ! may be a packed array of indexes
        !%-----------------------------------------------------------------------------
        outvalue = elemR(indx,er_Depth) * elemSGR(indx,esgr_Rectangular_Breadth)

    end function rectangular_area_from_depth
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function rectangular_area_from_depth_singular (indx) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for rectangular cross section of a single element
        !% The input indx is the row index in full data 2D array.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), pointer :: depth(:), breadth(:)
        !%-----------------------------------------------------------------------------
        depth   => elemR(:,er_Depth)
        breadth => elemSGR(:,esgr_Rectangular_Breadth)
        !%-----------------------------------------------------------------------------
        outvalue = depth(indx) * breadth(indx)

    end function rectangular_area_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_topwidth_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a rectangular channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:), elemType(:)
        real(8), pointer :: breadth(:), topwidth(:), depth(:), fullDepth(:)
        !%-----------------------------------------------------------------------------
        thisP     => elemPGx(1:Npack,thisCol) 
        elemType  => elemI(:,ei_elementType)
        topwidth  => elemR(:,er_Topwidth)
        depth     => elemR(:,er_Depth)
        fullDepth => elemR(:,er_FullDepth)
        breadth   => elemSGR(:,esgr_Rectangular_Breadth)
        !%-----------------------------------------------------------------------------

        topwidth(thisP) = breadth(thisP)

        !% HACK code: testing if rectangular open and closed can be
        !% incorporated in a single piece of code
        where ((elemType(thisP) == rectangular_closed) .and. &
               (depth(thisP) >= fullDepth(thisP)))
               topwidth = setting%ZeroValue%Topwidth
        end where

    end subroutine rectangular_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function rectangular_topwidth_from_depth_singular (indx) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth for a rectangular cross section of a single element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx 
        !%-----------------------------------------------------------------------------
        !%  
        outvalue = elemSGR(indx,esgr_Rectangular_Breadth)

    end function rectangular_topwidth_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_perimeter_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a rectangular channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: breadth(:), depth(:), perimeter(:)
        !%-----------------------------------------------------------------------------
        thisP     => elemPGx(1:Npack,thisCol) 
        breadth   => elemSGR(:,esgr_Rectangular_Breadth)
        depth     => elemR(:,er_Depth)
        perimeter => elemR(:,er_Perimeter)
        !%-----------------------------------------------------------------------------

        perimeter(thisP) = twoR * depth(thisP) + breadth(thisP) 

    end subroutine rectangular_perimeter_from_depth
!%    
!%==========================================================================    
!%==========================================================================
!%
    real(8) function rectangular_perimeter_from_depth_singular (indx) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes wetted perimeter from known depth for a rectangular cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), pointer :: depth(:), breadth(:)
        !%-----------------------------------------------------------------------------
        depth   => elemR(:,er_Depth)
        breadth => elemSGR(:,esgr_Rectangular_Breadth)
        !%-----------------------------------------------------------------------------
        
        outvalue = twoR * depth(indx) + breadth(indx)

    end function rectangular_perimeter_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_hyddepth_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the hydraulic (average) depth from a known depth in a rectangular channel
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

        hyddepth(thisP) = depth(thisP)

    end subroutine rectangular_hyddepth_from_depth
!%    
!%==========================================================================  
!%==========================================================================
!%
    real(8) function rectangular_hyddepth_from_depth_singular (indx) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic depth from known depth for rectangular cross section of 
        !% a single element
        !%-----------------------------------------------------------------------------   
        integer, intent(in) :: indx     
        !%-----------------------------------------------------------------------------  

        outvalue = elemR(indx,er_Depth)

    end function rectangular_hyddepth_from_depth_singular 
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function rectangular_hydradius_from_depth_singular (indx) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic radius from known depth for a rectangular cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), pointer :: depth(:), breadth(:)
        !%-----------------------------------------------------------------------------
        depth   => elemR(:,er_Depth)
        breadth => elemSGR(:,esgr_Rectangular_Breadth)
        !%-----------------------------------------------------------------------------
        
        outvalue = (depth(indx) * breadth(indx)) / ( twoR * depth(indx) + breadth(indx) )

    end function rectangular_hydradius_from_depth_singular
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
!%
!%    
!%==========================================================================
!% PRIVATE
!%==========================================================================   
!%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%  

       !%==========================================================================   
    ! !%
    ! subroutine rectangular_open_head_from_volume (elemPGx, Npack, thisCol)
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Only applies on open channels (or non-surcharged rectangular conduits)
    !     !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
    !     !% Assumes that volume > 0 is enforced in volume computations.
    !     !%-----------------------------------------------------------------------------
    !     integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
    !     integer, pointer :: thisP(:)
    !     real(8), pointer :: head(:), volume(:), length(:), breadth(:), zbottom(:)
    !     !%-----------------------------------------------------------------------------
    !     thisP   => elemPGx(1:Npack,thisCol) 
    !     head    => elemR(:,er_Head)
    !     volume  => elemR(:,er_Volume)
    !     length  => elemR(:,er_Length)
    !     breadth => elemSGR(:,esgr_Rectangular_Breadth)
    !     zbottom => elemR(:,er_Zbottom)
    !     !%-----------------------------------------------------------------------------   

    !     head(thisP) = zbottom(thisP) + volume(thisP) / (length(thisP) * breadth(thisP))
   
    ! end subroutine rectangular_open_head_from_volume
    !%  
    !%==========================================================================
    !%    !%==========================================================================
    !%
    ! subroutine rectangular_area_from_depth (elemPGx, Npack, thisCol)
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Computes area of a rectangular open channel given its depth
    !     !% Note, does NOT consider any closed top!
    !     !%-----------------------------------------------------------------------------
    !     integer, target, intent(in) :: elemPGx(:,:)
    !     integer, intent(in) ::  Npack, thisCol
    !     integer, pointer :: thisP(:)
    !     real(8), pointer :: area(:), depth(:), breadth(:)
    !     !%-----------------------------------------------------------------------------
    !     thisP   => elemPGx(1:Npack,thisCol) 
    !     area    => elemR(:,er_Area)
    !     depth   => elemR(:,er_Depth)
    !     breadth => elemSGR(:,esgr_Rectangular_Breadth)
    !     !%-----------------------------------------------------------------------------

    !     area(thisP) = depth(thisP) * breadth(thisP)

    ! end subroutine rectangular_area_from_depth
    ! !%
    ! !%==========================================================================
    !%==========================================================================
    !% END OF MODULE
    !%+=========================================================================
end module rectangular_channel