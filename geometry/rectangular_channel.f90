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
    !public :: rectangular_normaldepth_from_depth_singular
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
    ! real(8) function rectangular_depth_from_sectionfactor_singular &
    !     (sectionFactor, sectionFactorMax) result (outvalue)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Computes the depth for a given section factor for a
    !     !% rectangular channel
    !     !%------------------------------------------------------------------
    !     real(8), intent(in) :: sectionFactor, sectionFactorMax
    !     !%-------------------------------------------------------------------
    !     !%-------------------------------------------------------------------

    !     sfNorm = sectionFactor / sectionFactorMax
        
        

    ! end function rectangular_normaldepth_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_depth_from_volume (elemPGx, Npack, thisCol)
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
            real(8), pointer :: fulldepth(:), fullvolume(:)
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (Npack < 1) return
        !%------------------------------------------------------------------  
        !% Aliases:  
            thisP      => elemPGx(1:Npack,thisCol) 
            depth      => elemR(:,er_Depth)
            volume     => elemR(:,er_Volume)
            length     => elemR(:,er_Length)
            fulldepth  => elemR(:,er_FullDepth)
            fullvolume => elemR(:,er_FullVolume)
            breadth    => elemSGR(:,esgr_Rectangular_Breadth)
        !%----------------------------------------------------------------- 

        where (volume(thisP) >= fullvolume(thisP))
            depth(thisP) = fulldepth(thisP)
        elsewhere
            depth(thisP) = volume(thisP) / (length(thisP) * breadth(thisP))
        endwhere

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
    subroutine rectangular_topwidth_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a rectangular channel
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:), GeomType(:)
        real(8), pointer :: breadth(:), topwidth(:), depth(:), fullDepth(:)
        !%-----------------------------------------------------------------------------
        thisP     => elemPGx(1:Npack,thisCol) 
        GeomType  => elemI(:,ei_geometryType)
        topwidth  => elemR(:,er_Topwidth)
        depth     => elemR(:,er_Depth)
        fullDepth => elemR(:,er_FullDepth)
        breadth   => elemSGR(:,esgr_Rectangular_Breadth)
        !%-----------------------------------------------------------------------------

        topwidth(thisP) = breadth(thisP)

    end subroutine rectangular_topwidth_from_depth
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
!% SINGULAR
!%==========================================================================
!%
    real(8) function rectangular_area_from_depth_singular (indx, depth) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes area from known depth for rectangular cross section of a single element
        !% The input indx is the row index in full data 2D array.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer ::  breadth(:)
        !%-----------------------------------------------------------------------------
        breadth => elemSGR(:,esgr_Rectangular_Breadth)
        !%-----------------------------------------------------------------------------
       ! print *, 'in rectangular_area_from_depth_singular'
       ! print *, indx, breadth(indx), depth
        outvalue = depth * breadth(indx)

    end function rectangular_area_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function rectangular_topwidth_from_depth_singular (indx, depth) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth for a rectangular cross section of a single element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx 
        real(8), intent(in) :: depth
        !%-----------------------------------------------------------------------------
        !%  
        outvalue = elemSGR(indx,esgr_Rectangular_Breadth)

    end function rectangular_topwidth_from_depth_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function rectangular_perimeter_from_depth_singular (indx, depth) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes wetted perimeter from known depth for a rectangular cross section of
        !% a single element 
        !%-----------------------------------------------------------------------------
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: indx
        real(8), intent(in) :: depth
        real(8), pointer ::  breadth(:)
        !%-----------------------------------------------------------------------------
        breadth => elemSGR(:,esgr_Rectangular_Breadth)
        !%-----------------------------------------------------------------------------
        
        outvalue = twoR * depth + breadth(indx)

    end function rectangular_perimeter_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function rectangular_hyddepth_from_depth_singular (indx,depth) result (outvalue)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic depth from known depth for rectangular cross section of 
        !% a single element
        !%-----------------------------------------------------------------------------   
        integer, intent(in) :: indx   
        real(8), intent(in) :: depth  
        !%-----------------------------------------------------------------------------  

        outvalue = depth

    end function rectangular_hyddepth_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
    real(8) function rectangular_hydradius_from_depth_singular &
        (indx, depth) result (outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes hydraulic radius from known depth for a rectangular cross section of
        !% a single element 
        !%------------------------------------------------------------------
            integer, intent(in) :: indx
            real(8), intent(in) :: depth
            real(8), pointer :: breadth(:)
        !%------------------------------------------------------------------
        breadth => elemSGR(:,esgr_Rectangular_Breadth)
        !%------------------------------------------------------------------
        
        outvalue = (depth * breadth(indx)) / ( twoR * depth + breadth(indx) )

    end function rectangular_hydradius_from_depth_singular
!%      
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module rectangular_channel