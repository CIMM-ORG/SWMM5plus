module rectangular_channel

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use geometry_lowlevel

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Rectangular channel geometry
    !%

    private

    public :: rectangular_depth_from_volume
    public :: rectangular_topwidth_from_depth    
    public :: rectangular_perimeter_from_depth

    ! public :: rectangular_area_from_depth_singular
    ! public :: rectangular_topwidth_from_depth_singular 
    ! public :: rectangular_perimeter_from_depth_singular

    !public :: rectangular_hyddepth_from_depth

    !public :: rectangular_hyddepth_from_depth_singular
   !public :: rectangular_hydradius_from_depth_singular

    contains

!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine rectangular_depth_from_volume (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels 
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is previously enforced in volume computations.
        !%------------------------------------------------------------------
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: depth(:), volume(:)
            real(8), pointer :: fulldepth(:), fullvolume(:)
        !%------------------------------------------------------------------  
        !% Aliases:  
            depth      => elemR(:,er_Depth)
            volume     => elemR(:,er_Volume)
            fulldepth  => elemR(:,er_FullDepth)
            fullvolume => elemR(:,er_FullVolume)
        !%----------------------------------------------------------------- 

        if (setting%Discretization%AllowChannelOverflowTF) then
            where (volume(thisP) >= fullvolume(thisP))
                !% --- truncate depth at full depth if there is overflow
                depth(thisP) = fulldepth(thisP)
            elsewhere
                !% --- standard rectangular depth
                depth(thisP) = llgeo_rectangular_depth_from_volume_pure &
                                    (thisP, volume(thisP))
            endwhere
        else 
            !% --- if overflow is NOT allowed
            !% --- standard rectangular depth if overflow is suppressed
            depth(thisP) = llgeo_rectangular_depth_from_volume_pure &
                                    (thisP, volume(thisP))
        end if

    end subroutine rectangular_depth_from_volume
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_topwidth_from_depth (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a rectangular channel
        !%------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            real(8), pointer ::  topwidth(:), depth(:)
        !%------------------------------------------------------------------
        !% Aliases
            topwidth  => elemR(:,er_Topwidth)
            depth     => elemR(:,er_Depth)
        !%-----------------------------------------------------------------

        topwidth(thisP) =  llgeo_rectangular_topwidth_from_depth_pure(thisP,depth(thisP))

    end subroutine rectangular_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_perimeter_from_depth (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a rectangular channel
        !%------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: perimeter(:), depth(:)
        !%-------------------------------------------------------------------
        !% Aliases:
            perimeter => elemR(:,er_Perimeter)
            depth     => elemR(:,er_Depth)
        !%-------------------------------------------------------------------

        perimeter(thisP) = llgeo_rectangular_perimeter_from_depth_pure(thisP,depth(thisP))

    end subroutine rectangular_perimeter_from_depth
!%    
!%==========================================================================  
!%==========================================================================
!%
    ! subroutine rectangular_hyddepth_from_depth (elemPGx, Npack, thisCol)
    !     !%  
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Computes the hydraulic (average) depth from a known depth in a rectangular channel
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

    ! end subroutine rectangular_hyddepth_from_depth
!%    
!%==========================================================================   
!% SINGULAR
!%==========================================================================
!%
!     pure real(8) function rectangular_area_from_depth_singular (indx, depth) result (outvalue)
!         !%-------------------------------------------------------------------
!         !% Description:
!         !% Computes area from known depth for rectangular cross section of a single element
!         !% The input indx is the row index in full data 2D array.
!         !%--------------------------------------------------------------------
!         !% Declarations
!             integer, intent(in) :: indx
!             real(8), intent(in) :: depth
!         !%---------------------------------------------------------------------
!         !%---------------------------------------------------------------------

!         outvalue = depth * elemSGR(indx,esgr_Rectangular_Breadth)

!     end function rectangular_area_from_depth_singular
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     pure real(8) function rectangular_topwidth_from_depth_singular &
!         (indx, depth) result (outvalue)
!         !%------------------------------------------------------------------
!         !% Description:
!         !% Computes the topwidth for a rectangular cross section of a single element
!         !%------------------------------------------------------------------
!         !% Declarations
!             integer, intent(in) :: indx 
!             real(8), intent(in) :: depth
!         !%------------------------------------------------------------------

!         outvalue = elemSGR(indx,esgr_Rectangular_Breadth)

!     end function rectangular_topwidth_from_depth_singular
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     pure real(8) function rectangular_perimeter_from_depth_singular &
!          (indx, depth) result (outvalue)
!         !%------------------------------------------------------------------
!         !% Description:
!         !% Computes wetted perimeter from known depth for a rectangular cross section of
!         !% a single element 
!         !%-------------------------------------------------------------------
!         !% Declarations
!             integer, intent(in) :: indx
!             real(8), intent(in) :: depth
!         !%-------------------------------------------------------------------
        
!         outvalue = twoR * depth + elemSGR(indx,esgr_Rectangular_Breadth)

!     end function rectangular_perimeter_from_depth_singular
! !%    
! !%==========================================================================
! !%==========================================================================
! !%
!     ! real(8) function rectangular_hyddepth_from_depth_singular (indx,depth) result (outvalue)
!     !     !%  
!     !     !%-----------------------------------------------------------------------------
!     !     !% Description:
!     !     !% Computes hydraulic depth from known depth for rectangular cross section of 
!     !     !% a single element
!     !     !%-----------------------------------------------------------------------------   
!     !     integer, intent(in) :: indx   
!     !     real(8), intent(in) :: depth  
!     !     !%-----------------------------------------------------------------------------  

!     !     outvalue = depth

!     ! end function rectangular_hyddepth_from_depth_singular
! !%    
! !%==========================================================================
! !%==========================================================================
! !%
! !     real(8) function rectangular_hydradius_from_depth_singular &
! !         (indx, depth) result (outvalue)
! !         !%------------------------------------------------------------------
! !         !% Description:
! !         !% Computes hydraulic radius from known depth for a rectangular cross section of
! !         !% a single element 
! !         !%------------------------------------------------------------------
! !             integer, intent(in) :: indx
! !             real(8), intent(in) :: depth
! !             real(8), pointer :: breadth(:)
! !         !%------------------------------------------------------------------
! !         breadth => elemSGR(:,esgr_Rectangular_Breadth)
! !         !%------------------------------------------------------------------
        
! !         outvalue = (depth * breadth(indx)) / ( twoR * depth + breadth(indx) )

! !     end function rectangular_hydradius_from_depth_singular
! !%      
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module rectangular_channel