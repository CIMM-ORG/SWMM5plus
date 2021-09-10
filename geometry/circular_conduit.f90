module circular_conduit

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use xsect_tables

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Circular conduit geometry
    !%

    private

    public :: circular_depth_from_volume
    ! public :: circular_area_from_depth_singular
    ! public :: circular_topwidth_from_depth
    ! public :: circular_topwidth_from_depth_singular 
    ! public :: circular_perimeter_from_depth
    ! public :: circular_perimeter_from_depth_singular
    ! public :: circular_hyddepth_from_depth
    ! public :: circular_hyddepth_from_depth_singular
    ! public :: circular_hydradius_from_depth_singular


    contains
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    !%
    subroutine circular_depth_from_volume (elemPGx, Npack, thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Only applies on open conduits (or non-surcharged circular conduits)
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !% NOTE: this does NOT limit the depth by surcharge height at this point
        !% This will be done after the head is computed.
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:), volume(:), length(:), AoverAfull(:), YoverYfull(:)
        real(8), pointer :: fullArea(:)
        !%-----------------------------------------------------------------------------
        thisP      => elemPGx(1:Npack,thisCol) 
        depth      => elemR(:,er_Depth)
        volume     => elemR(:,er_Volume)
        length     => elemR(:,er_Length)
        fullArea   => elemR(:,er_FullArea)
        AoverAfull => elemSGR(:,eSGR_Circular_AoverAfull)
        YoverYfull => elemSGR(:,eSGR_Circular_YoverYfull)
        !%-----------------------------------------------------------------------------  

        AoverAfull(thisP) = volume(thisP) / (length(thisP) * fullArea(thisP)) 

    end subroutine circular_depth_from_volume
    !%  
    !%==========================================================================
    !%==========================================================================
    !%
    ! real(8) function circular_area_from_depth_singular (indx) result (outvalue)
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Computes area from known depth for circular cross section of a single element
    !     !% The input indx is the row index in full data 2D array.
    !     !%-----------------------------------------------------------------------------
    !     integer, intent(in) :: indx
    !     real(8), pointer :: depth(:), breadth(:)
    !     !%-----------------------------------------------------------------------------
    !     depth   => elemR(:,er_Depth)
    !     breadth => elemSGR(:,eSGR_circular_Breadth)
    !     !%-----------------------------------------------------------------------------
    !     outvalue = depth(indx) * breadth(indx)

    ! end function circular_area_from_depth_singular
    ! !%
    ! !%==========================================================================
    ! !%==========================================================================
    ! !%
    ! subroutine circular_topwidth_from_depth (elemPGx, Npack, thisCol)
    !     !%  
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Computes the topwidth from a known depth in a circular conduit
    !     !%-----------------------------------------------------------------------------
    !     integer, target, intent(in) :: elemPGx(:,:)
    !     integer, intent(in) ::  Npack, thisCol
    !     integer, pointer :: thisP(:)
    !     real(8), pointer :: breadth(:), topwidth(:)
    !     !%-----------------------------------------------------------------------------
    !     thisP    => elemPGx(1:Npack,thisCol) 
    !     topwidth => elemR(:,er_Topwidth)
    !     breadth  => elemSGR(:,eSGR_circular_Breadth)
    !     !%-----------------------------------------------------------------------------

    !     topwidth(thisP) = breadth(thisP)

    ! end subroutine circular_topwidth_from_depth
    ! !%    
    ! !%==========================================================================
    ! !%==========================================================================
    ! !%
    ! real(8) function circular_topwidth_from_depth_singular (indx) result (outvalue)
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Computes the topwidth for a circular cross section of a single element
    !     !%-----------------------------------------------------------------------------
    !     integer, intent(in) :: indx 
    !     !%-----------------------------------------------------------------------------
    !     !%  
    !     outvalue = elemSGR(indx,eSGR_circular_Breadth)

    ! end function circular_topwidth_from_depth_singular
    ! !%
    ! !%==========================================================================
    ! !%==========================================================================
    ! !%
    ! subroutine circular_perimeter_from_depth (elemPGx, Npack, thisCol)
    !     !%  
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Computes the perimeter from a known depth in a circular conduit
    !     !%-----------------------------------------------------------------------------
    !     integer, target, intent(in) :: elemPGx(:,:)
    !     integer, intent(in) ::  Npack, thisCol
    !     integer, pointer :: thisP(:)
    !     real(8), pointer :: breadth(:), depth(:), perimeter(:)
    !     !%-----------------------------------------------------------------------------
    !     thisP     => elemPGx(1:Npack,thisCol) 
    !     breadth   => elemSGR(:,eSGR_circular_Breadth)
    !     depth     => elemR(:,er_Depth)
    !     perimeter => elemR(:,er_Perimeter)
    !     !%-----------------------------------------------------------------------------

    !     perimeter(thisP) = twoR * depth(thisP) + breadth(thisP) 

    ! end subroutine circular_perimeter_from_depth
    ! !%    
    ! !%==========================================================================    
    ! !%==========================================================================
    ! !%
    ! real(8) function circular_perimeter_from_depth_singular (indx) result (outvalue)
    !     !%  
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Computes wetted perimeter from known depth for a circular cross section of
    !     !% a single element 
    !     !%-----------------------------------------------------------------------------
    !     !%-----------------------------------------------------------------------------
    !     integer, intent(in) :: indx
    !     real(8), pointer :: depth(:), breadth(:)
    !     !%-----------------------------------------------------------------------------
    !     depth   => elemR(:,er_Depth)
    !     breadth => elemSGR(:,eSGR_circular_Breadth)
    !     !%-----------------------------------------------------------------------------
        
    !     outvalue = twoR * depth(indx) + breadth(indx)

    ! end function circular_perimeter_from_depth_singular
    ! !%    
    ! !%==========================================================================
    ! !%==========================================================================
    ! !%
    ! subroutine circular_hyddepth_from_depth (elemPGx, Npack, thisCol)
    !     !%  
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Computes the hydraulic (average) depth from a known depth in a circular conduit
    !     !%-----------------------------------------------------------------------------
    !     integer, target, intent(in) :: elemPGx(:,:)
    !     integer, intent(in) ::  Npack, thisCol
    !     integer, pointer :: thisP(:)
    !     real(8), pointer :: hyddepth(:), depth(:)
    !     !%-----------------------------------------------------------------------------
    !     thisP     => elemPGx(1:Npack,thisCol) 
    !     depth     => elemR(:,er_Depth)
    !     hyddepth  => elemR(:,er_HydDepth)
    !     !%-----------------------------------------------------------------------------

    !     hyddepth(thisP) = depth(thisP)

    ! end subroutine circular_hyddepth_from_depth
    ! !%    
    ! !%==========================================================================  
    ! !%==========================================================================
    ! !%
    ! real(8) function circular_hyddepth_from_depth_singular (indx) result (outvalue)
    !     !%  
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Computes hydraulic depth from known depth for circular cross section of 
    !     !% a single element
    !     !%-----------------------------------------------------------------------------   
    !     integer, intent(in) :: indx     
    !     !%-----------------------------------------------------------------------------  

    !     outvalue = elemR(indx,er_Depth)

    ! end function circular_hyddepth_from_depth_singular
    ! !% 
    ! !%==========================================================================

    ! !%==========================================================================
    ! !%
    ! real(8) function circular_hydradius_from_depth_singular (indx) result (outvalue)
    !     !%  
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Computes hydraulic radius from known depth for a circular cross section of
    !     !% a single element 
    !     !%-----------------------------------------------------------------------------
    !     integer, intent(in) :: indx
    !     real(8), pointer :: depth(:), breadth(:)
    !     !%-----------------------------------------------------------------------------
    !     depth   => elemR(:,er_Depth)
    !     breadth => elemSGR(:,eSGR_circular_Breadth)
    !     !%-----------------------------------------------------------------------------
        
    !     outvalue = (depth(indx) * breadth(indx)) / ( twoR * depth(indx) + breadth(indx) )

    ! end function circular_hydradius_from_depth_singular
    ! !%    
    ! !%==========================================================================

    ! !%
    ! !%    
    ! !%==========================================================================
    ! !%==========================================================================
    ! !%
    !     !%  
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% 
    !     !%-----------------------------------------------------------------------------

    !     !%-----------------------------------------------------------------------------
    !     !%  
    ! !%
    ! !%    
    ! !%==========================================================================
    ! !%==========================================================================
    ! !%
    !     !%  
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% 
    !     !%-----------------------------------------------------------------------------

    !     !%-----------------------------------------------------------------------------
    !     !%  
    ! !%
    ! !%    
    ! !%==========================================================================
    ! !%==========================================================================
    ! !%
    !     !%  
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% 
    !     !%-----------------------------------------------------------------------------

    !     !%-----------------------------------------------------------------------------
    !     !%  
    ! !%
    ! !%    
    ! !%==========================================================================
    ! !% PRIVATE
    ! !%==========================================================================   
    ! !%  
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% 
    !     !%-----------------------------------------------------------------------------

    !     !%-----------------------------------------------------------------------------
    !     !%  

    !    !%==========================================================================   
    ! ! !%
    ! ! subroutine circular_open_head_from_volume (elemPGx, Npack, thisCol)
    ! !     !%-----------------------------------------------------------------------------
    ! !     !% Description:
    ! !     !% Only applies on open conduits (or non-surcharged circular conduits)
    ! !     !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
    ! !     !% Assumes that volume > 0 is enforced in volume computations.
    ! !     !%-----------------------------------------------------------------------------
    ! !     integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
    ! !     integer, pointer :: thisP(:)
    ! !     real(8), pointer :: head(:), volume(:), length(:), breadth(:), zbottom(:)
    ! !     !%-----------------------------------------------------------------------------
    ! !     thisP   => elemPGx(1:Npack,thisCol) 
    ! !     head    => elemR(:,er_Head)
    ! !     volume  => elemR(:,er_Volume)
    ! !     length  => elemR(:,er_Length)
    ! !     breadth => elemSGR(:,eSGr_circular_Breadth)
    ! !     zbottom => elemR(:,er_Zbottom)
    ! !     !%-----------------------------------------------------------------------------   

    ! !     head(thisP) = zbottom(thisP) + volume(thisP) / (length(thisP) * breadth(thisP))
   
    ! ! end subroutine circular_open_head_from_volume
    ! !%  
    ! !%==========================================================================
    ! !%    !%==========================================================================
    ! !%
    ! ! subroutine circular_area_from_depth (elemPGx, Npack, thisCol)
    ! !     !%-----------------------------------------------------------------------------
    ! !     !% Description:
    ! !     !% Computes area of a circular open conduit given its depth
    ! !     !% Note, does NOT consider any closed top!
    ! !     !%-----------------------------------------------------------------------------
    ! !     integer, target, intent(in) :: elemPGx(:,:)
    ! !     integer, intent(in) ::  Npack, thisCol
    ! !     integer, pointer :: thisP(:)
    ! !     real(8), pointer :: area(:), depth(:), breadth(:)
    ! !     !%-----------------------------------------------------------------------------
    ! !     thisP   => elemPGx(1:Npack,thisCol) 
    ! !     area    => elemR(:,er_Area)
    ! !     depth   => elemR(:,er_Depth)
    ! !     breadth => elemSGR(:,eSGR_circular_Breadth)
    ! !     !%-----------------------------------------------------------------------------

    ! !     area(thisP) = depth(thisP) * breadth(thisP)

    ! ! end subroutine circular_area_from_depth
    ! ! !%
    ! ! !%==========================================================================
    ! !%==========================================================================
    ! !% END OF MODULE
    ! !%+=========================================================================
end module circular_conduit