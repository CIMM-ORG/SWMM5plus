module rectangular_conduit

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

    public :: rectangular_closed_depth_from_volume
    public :: rectangular_closed_topwidth_from_depth
    public :: rectangular_closed_perimeter_from_depth

    ! public :: rectangular_closed_area_from_depth
    ! !public :: rectangular_closed_area_from_depth_singular
    
    ! public :: rectangular_closed_topwidth_from_depth_singular 
    
    ! public :: rectangular_closed_perimeter_from_depth_singular
    ! public :: rectangular_closed_hyddepth_from_depth
    ! !public :: rectangular_closed_hyddepth_from_depth_singular
    ! public :: rectangular_closed_hydradius_from_depth_singular

    contains

!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine rectangular_closed_depth_from_volume (elemPGx, Npack, thisCol)
        !%-------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels 
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !%--------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
            integer, pointer :: thisP(:)
            real(8), pointer :: depth(:), volume(:), length(:), breadth(:)
            real(8), pointer :: fulldepth(:), fullvolume(:)
        !%--------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%--------------------------------------------------------------------
        !% Aliases
            thisP       => elemPGx(1:Npack,thisCol) 
            depth       => elemR(:,er_Depth)
            volume      => elemR(:,er_Volume)
            fulldepth   => elemR(:,er_FullDepth)
            fullvolume  => elemR(:,er_FullVolume)
            length      => elemR(:,er_Length)
            breadth     => elemSGR(:,esgr_Rectangular_Breadth)
        !%---------------------------------------------------------------------  
        
        depth(thisP) = volume(thisP) / (length(thisP) * breadth(thisP))

        !% --- ensure the full depth is not exceeded
        depth(thisP) = min(depth(thisP),fulldepth(thisP))
        

        ! where (volume(thisP) < fullvolume(thisP))
         !   depth(thisP) = volume(thisP) / (length(thisP) * breadth(thisP))
        ! else where (volume(thisP) >= fullvolume(thisP))
        !     depth(thisP) = fulldepth(thisP)
        ! end where

    end subroutine rectangular_closed_depth_from_volume
    !%  
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_closed_topwidth_from_depth (elemPGx, Npack, thisCol)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a rectangular channel
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: elemPGx(:,:)
            integer, intent(in) ::  Npack, thisCol
            integer, pointer :: thisP(:), GeomType(:)
            real(8), pointer :: breadth(:), topwidth(:), depth(:), fullDepth(:)
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%-------------------------------------------------------------------
        !% Aliases:
            thisP     => elemPGx(1:Npack,thisCol) 
            GeomType  => elemI(:,ei_geometryType)
            topwidth  => elemR(:,er_Topwidth)
            depth     => elemR(:,er_Depth)
            fullDepth => elemR(:,er_FullDepth)
            breadth   => elemSGR(:,esgr_Rectangular_Breadth)
        !%-------------------------------------------------------------------

        topwidth(thisP) = breadth(thisP)

        !% reset topwidth to zeroValue for depth higher than full depth
        where (depth(thisP) >= fullDepth(thisP))
               topwidth(thisP) = setting%ZeroValue%Topwidth
        end where

    end subroutine rectangular_closed_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_closed_perimeter_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a rectangular channel
        !%------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: elemPGx(:,:)
            integer, intent(in) ::  Npack, thisCol
            integer, pointer :: thisP(:)
            real(8), pointer :: breadth(:), depth(:), perimeter(:), fulldepth(:), fullPerimeter(:)
        !%------------------------------------------------------------------       
        !% Preliminaries
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases:
            thisP     => elemPGx(1:Npack,thisCol) 
            breadth   => elemSGR(:,esgr_Rectangular_Breadth)
            depth     => elemR(:,er_Depth)
            fulldepth => elemR(:,er_FullDepth)
            perimeter => elemR(:,er_Perimeter)
            fullPerimeter => elemR(:,er_FullPerimeter)
        !%-------------------------------------------------------------------
        
        where (depth(thisP) < fulldepth(thisP))
            perimeter(thisP) = twoR * depth(thisP) + breadth(thisP) 
        else where (depth(thisP) >= fulldepth(thisP))
            perimeter(thisP) = fullPerimeter(thisP)
        end where

    end subroutine rectangular_closed_perimeter_from_depth
!%    
!%========================================================================== 
!%==========================================================================
!%
!     elemental real(8) function rectangular_closed_area_from_depth &
!         (indx) result (outvalue)
!         !%-----------------------------------------------------------------------------
!         !% Description:
!         !% Computes area from known depth for rectangular cross section
!         !%-----------------------------------------------------------------------------
!         integer, intent(in) :: indx  ! may be a packed array of indexes
!         !%-----------------------------------------------------------------------------

!         if (elemR(indx,er_Depth) < elemR(indx,er_FullDepth)) then
!             outvalue = elemR(indx,er_Depth) * elemSGR(indx,esgr_Rectangular_Breadth)
!         else 
!             outvalue = elemR(indx,er_FullArea)
!         end if

!     end function rectangular_closed_area_from_depth
! !%
! !%==========================================================================

 
! !%==========================================================================
! !%
!     subroutine rectangular_closed_hyddepth_from_depth (elemPGx, Npack, thisCol)
!         !%  
!         !%-----------------------------------------------------------------------------
!         !% Description:
!         !% Computes the hydraulic (average) depth from a known depth in a rectangular channel
!         !%-----------------------------------------------------------------------------
!         integer, target, intent(in) :: elemPGx(:,:)
!         integer, intent(in) ::  Npack, thisCol
!         integer, pointer :: thisP(:)
!         real(8), pointer :: hyddepth(:), depth(:), fulldepth(:)
!         !%-----------------------------------------------------------------------------
!         thisP     => elemPGx(1:Npack,thisCol) 
!         depth     => elemR(:,er_Depth)
!         hyddepth  => elemR(:,er_HydDepth)
!         fulldepth => elemR(:,er_FullDepth)
!         !%-----------------------------------------------------------------------------
!         where (depth(thisP) < fulldepth(thisP))
!             hyddepth(thisP) = depth(thisP)
!         else where (depth(thisP) >= fulldepth(thisP))
!             hyddepth(thisP) = fulldepth(thisP)
!         end where

!     end subroutine rectangular_closed_hyddepth_from_depth
! !%    
! !%==========================================================================  
! !% SINGULAR
! !%==========================================================================
! !%
!     ! real(8) function rectangular_closed_area_from_depth_singular &
!     !     (indx, depth) result (outvalue)
!     !     !%------------------------------------------------------------------
!     !     !% Description:
!     !     !% Computes area from known depth for rectangular cross section of a single element
!     !     !% The input indx is the row index in full data 2D array.
!     !     !%-------------------------------------------------------------------
!     !     integer, intent(in) :: indx
!     !     real(8), intent(in) :: depth
!     !     real(8), pointer ::  breadth(:), fulldepth(:)
!     !     !%--------------------------------------------------------------------
!     !     breadth     => elemSGR(:,esgr_Rectangular_Breadth)
!     !     fulldepth   => elemR(:,er_FullDepth)
!     !     !%-----------------------------------------------------------------------------
!     !     if (depth < fulldepth(indx)) then
!     !         outvalue = depth * breadth(indx)
!     !     else
!     !         outvalue = fulldepth(indx)
!     !     end if

!     ! end function rectangular_closed_area_from_depth_singular
! !%
! !%==========================================================================


! !%==========================================================================
! !%
!     ! real(8) function rectangular_closed_hyddepth_from_depth_singular &
!     !     (indx,depth) result (outvalue)
!     !     !%  
!     !     !%-----------------------------------------------------------------------------
!     !     !% Description:
!     !     !% Computes hydraulic depth from known depth for rectangular cross section of 
!     !     !% a single element
!     !     !%-----------------------------------------------------------------------------   
!     !     integer, intent(in) :: indx   
!     !     real(8), intent(in) :: depth  
!     !     !%-----------------------------------------------------------------------------  

!     !     if (depth <= setting%ZeroValue%Depth) then
!     !         !% --- empty
!     !         outvalue = setting%ZeroValue%Depth
!     !     elseif (depth >= elemR(indx,er_FullDepth))
!     !         !% --- full
!     !         outvalue = elemR(indx,er_FullHydDepth)
!     !     else
!     !         outvalue = depth
!     !     end if


!     ! end function rectangular_closed_hyddepth_from_depth_singular 
! !%    
! !%==========================================================================
! !%==========================================================================
! !%
!     real(8) function rectangular_closed_hydradius_from_depth_singular &
!         (indx, depth) result (outvalue)
!         !%  
!         !%-----------------------------------------------------------------------------
!         !% Description:
!         !% Computes hydraulic radius from known depth for a rectangular cross section of
!         !% a single element 
!         !%-----------------------------------------------------------------------------
!         integer, intent(in) :: indx
!         real(8), intent(in) :: depth
!         real(8), pointer :: breadth(:), fulldepth(:)
!         !%-----------------------------------------------------------------------------
!         breadth   => elemSGR(:,esgr_Rectangular_Breadth)
!         fulldepth => elemR(:,er_FullDepth)
!         !%-----------------------------------------------------------------------------
!         if (depth < fulldepth(indx)) then
!             outvalue = (depth * breadth(indx)) / ( twoR * depth + breadth(indx) )
!         else 
!             outvalue = (fulldepth(indx) * breadth(indx)) / ( twoR * fulldepth(indx) + breadth(indx) )
!         end if

!     end function rectangular_closed_hydradius_from_depth_singular
!%  
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module rectangular_conduit