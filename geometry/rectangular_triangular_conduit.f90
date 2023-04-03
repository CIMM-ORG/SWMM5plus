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
    public :: rectangular_triangular_topwidth_from_depth
    public :: rectangular_triangular_perimeter_from_depth
    
    ! public :: rectangular_triangular_area_from_depth
    ! public :: rectangular_triangular_area_from_depth_singular
    
    ! public :: rectangular_triangular_topwidth_from_depth_singular 
    
    ! public :: rectangular_triangular_perimeter_from_depth_singular
    ! public :: rectangular_triangular_hyddepth_from_depth
    ! !public :: rectangular_triangular_hyddepth_from_depth_singular
    ! public :: rectangular_triangular_hydradius_from_depth_singular

    contains

!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine rectangular_triangular_depth_from_volume (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels 
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is enforced in volume computations.
        !% NOTE: this does NOT limit the depth by surcharge height at this point
        !% This will be done after the head is computed.
        !%-------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: depth(:), bottomDepth(:), bottomArea(:)
            real(8), pointer :: fulldepth(:)
            real(8), pointer :: volume(:), length(:), breadth(:), bottomSlope(:)

        !%-------------------------------------------------------------------
        !% Aliases
            depth       => elemR(:,er_Depth)
            fulldepth   => elemR(:,er_FullDepth)
            breadth     => elemR(:,er_FullTopWidth)
            bottomDepth => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth) 
            bottomArea  => elemSGR(:,esgr_Rectangular_Triangular_BottomArea)
            bottomSlope => elemSGR(:,esgr_Rectangular_Triangular_BottomSlope)
            volume      => elemR(:,er_Volume)
            length      => elemR(:,er_Length)
        !%----------------------------------------------------------------------

        where(volume(thisP) <= bottomArea(thisP) * length(thisP))
            depth(thisP) = sqrt(volume(thisP) / (length(thisP) * bottomSlope(thisP)))

        elsewhere(volume(thisP) > bottomarea(thisP)*length(thisP))
            depth(thisP) = bottomDepth(thisP) &
                + ((volume(thisP) / length(thisP)) - bottomarea(thisP)) / breadth(thisP)
        end where

        !% --- ensure the full depth is not exceeded
        depth(thisP) = min(depth(thisP),fulldepth(thisP))
                
    end subroutine rectangular_triangular_depth_from_volume
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_triangular_topwidth_from_depth (thisP)
        !%  
        !%-------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a rectangular_triangular channel
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: breadth(:), topwidth(:), bottomSlope(:), depth(:), bottomDepth(:)
        !%--------------------------------------------------------------------
        !% Aliases
            topwidth    => elemR(:,er_Topwidth)
            depth       => elemR(:,er_Depth)
            breadth     => elemR(:,er_FullTopWidth)
            bottomSlope => elemSGR(:,esgr_Rectangular_Triangular_BottomSlope)
            bottomDepth => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth) 
        !%---------------------------------------------------------------------

        where(depth(thisP) <= bottomDepth(thisP))
            topwidth(thisP) = twoR * bottomSlope(thisP) * depth(thisP)

        else where(depth(thisP) > bottomDepth(thisP))
            topwidth(thisP) = breadth(thisP)
        endwhere

    end subroutine rectangular_triangular_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine rectangular_triangular_perimeter_from_depth (thisP)
        !%-------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a rectangular_triangular channel
        !%-------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: breadth(:), depth(:), bottomSlope(:), bottomDepth(:)
            real(8), pointer :: perimeter(:), fullPerimeter(:), fullDepth(:)
        !%-------------------------------------------------------------------
        !% Aliases:
            breadth       => elemR(:,er_FullTopWidth)
            depth         => elemR(:,er_Depth)
            bottomSlope   => elemSGR(:,esgr_Rectangular_Triangular_BottomSlope)
            bottomDepth   => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth)  
            perimeter     => elemR(:,er_Perimeter)
            fullperimeter => elemR(:,er_FullPerimeter)
            fulldepth     => elemR(:,er_FullDepth)
        !%--------------------------------------------------------------------

        where (depth(thisP) <= setting%ZeroValue%Depth)
            !% --- negligible water level
            perimeter(thisP) = setting%ZeroValue%Topwidth

        elsewhere ( (depth(thisP) > setting%ZeroValue%Depth) .and. (depth(thisP) <= bottomdepth(thisP)) )
            !% --- water level in lower triangular portion
            perimeter(thisP) = twoR * depth(thisP) * sqrt(oneR + bottomSlope(thisP) ** twoR)

        elsewhere ( (depth(thisP) > bottomdepth(thisP)) .and. (depth(thisP) < fulldepth(thisP)) )
            !% --- water level in upper rectangular portion
            perimeter(thisP) = twoR * bottomdepth(thisP) * sqrt(oneR + bottomSlope(thisP) ** twoR)  & !triangular section
                                 + twoR * (depth(thisP) - bottomDepth(thisP))                         !rectangular section
        elsewhere
            !% --- water level full
            perimeter(thisP) = fullperimeter(thisP)

        end where

    end subroutine rectangular_triangular_perimeter_from_depth
!%    
!%==========================================================================  


!%==========================================================================
!%
!     elemental real(8) function rectangular_triangular_area_from_depth &
!         (indx) result (outvalue)
!         !%-----------------------------------------------------------------------------
!         !% Description:
!         !% Computes area from known depth for rectangular cross section
!         !%-----------------------------------------------------------------------------
!         integer, intent(in) :: indx  ! may be a packed array of indexes
!         !%-----------------------------------------------------------------------------
        
!         if (elemR(indx,er_Depth) <= elemSGR(indx,esgr_Rectangular_Triangular_BottomDepth)) then
!             outvalue = elemR(indx,er_Depth) * elemR(indx,er_Depth) * elemSGR(indx,esgr_Rectangular_Triangular_BottomSlope)
!         else
!             outvalue =  elemSGR(indx,esgr_Rectangular_Triangular_BottomArea) &    !triangular section
!                         + ((elemR(indx,er_Depth)  - elemSGR(indx,esgr_Rectangular_Triangular_BottomDepth)) * &
!                         elemSGR(indx,esgr_rectangular_Triangular_TopBreadth))       !rectangular section
!         endif

!     end function rectangular_triangular_area_from_depth
! !%
! !%==========================================================================


! !%==========================================================================
! !%
!     subroutine rectangular_triangular_hyddepth_from_depth (elemPGx, Npack, thisCol)
!         !%  
!         !%-----------------------------------------------------------------------------
!         !% Description:
!         !% Computes the hydraulic (average) depth from a known depth in a rectangular_triangular channel
!         !%-----------------------------------------------------------------------------
!         integer, target, intent(in) :: elemPGx(:,:)
!         integer, intent(in) ::  Npack, thisCol
!         integer, pointer :: thisP(:)
!         real(8), pointer :: hyddepth(:), depth(:),  bottomdepth(:)
!         !%-----------------------------------------------------------------------------
!         thisP       => elemPGx(1:Npack,thisCol) 
!         depth       => elemR(:,er_Depth)
!         hyddepth    => elemR(:,er_HydDepth)
!         bottomdepth => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth) 
!         !%-----------------------------------------------------------------------------

!         where(depth(thisP) <= bottomdepth(thisP))
!             hyddepth(thisP) = depth(thisP) / twoR

!         else where(depth(thisP) > bottomdepth(thisP))
!             hyddepth(thisP) = (bottomdepth(thisP) / twoR) &         !triangular section
!                             + (depth(thisP) - bottomdepth(thisP))   !rectangular Section
!         endwhere

!         if (setting%Debug%File%geometry) &
!             print *, 'hydepth = ' , hyddepth(thisP)

!     end subroutine rectangular_triangular_hyddepth_from_depth
! !%    
! !%==========================================================================  
! !% SINGULAR
! !%==========================================================================
! !%
!     ! real(8) function rectangular_triangular_area_from_depth_singular &
!     !     (indx, depth) result (outvalue)
!     !     !%-----------------------------------------------------------------------------
!     !     !% Description:
!     !     !% Computes area from known depth for rectangular_triangular cross section of a single element
!     !     !% The input indx is the row index in full data 2D array.
!     !     !%-----------------------------------------------------------------------------
!     !     integer, intent(in) :: indx
!     !     real(8), intent(in) :: depth
!     !     real(8), pointer :: bottomDepth(:), bottomSlope(:), bottomArea(:), breadth(:)
!     !     !%-----------------------------------------------------------------------------
!     !     bottomDepth => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth) 
!     !     bottomSlope => elemSGR(:,esgr_Rectangular_Triangular_BottomSlope)
!     !     bottomArea  => elemSGR(:,esgr_Rectangular_Triangular_BottomArea)
!     !     breadth     => elemSGR(:,esgr_Rectangular_Triangular_TopBreadth)
!     !     !%-----------------------------------------------------------------------------
        
!     !     if(depth <= bottomDepth(indx)) then
!     !         outvalue = depth * depth * bottomSlope(indx)
!     !     else
!     !         outvalue = bottomArea(indx) + (depth - bottomDepth(indx)) * breadth(indx)    
!     !     endif

!     !     if (setting%Debug%File%geometry) &
!     !        print *, 'area = ' , outvalue

!     ! end function rectangular_triangular_area_from_depth_singular
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     ! real(8) function rectangular_triangular_topwidth_from_depth_singular &
!     !      (indx, depth) result (outvalue)
!     !     !%-----------------------------------------------------------------------------
!     !     !% Description:
!     !     !% Computes the topwidth for a rectangular_triangular cross section of a single element
!     !     !%-----------------------------------------------------------------------------
!     !     integer, intent(in) :: indx 
!     !     real(8), intent(in) :: depth
!     !     real(8), pointer :: bottomDepth(:), bottomSlope(:), breadth(:)
!     !     !%-----------------------------------------------------------------------------
!     !     bottomSlope => elemSGR(:,esgr_Rectangular_Triangular_BottomSlope)
!     !     bottomDepth => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth) 
!     !     breadth     => elemSGR(:,esgr_Rectangular_Triangular_TopBreadth)
!     !     !%-----------------------------------------------------------------------------
         
!     !     if(depth <= bottomDepth(indx)) then
!     !         outvalue = twoR * bottomSlope(indx) * depth
!     !     else
!     !         outvalue = breadth(indx)
!     !     endif

!     !     if (setting%Debug%File%geometry) &
!     !         print *, 'topwidth = ' , outvalue
!     ! end function rectangular_triangular_topwidth_from_depth_singular
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     ! real(8) function rectangular_triangular_perimeter_from_depth_singular &
!     !     (indx, depth) result (outvalue)
!     !     !%  
!     !     !%-----------------------------------------------------------------------------
!     !     !% Description:
!     !     !% Computes wetted perimeter from known depth for a rectangular_triangular cross section of
!     !     !% a single element 
!     !     !%-----------------------------------------------------------------------------
!     !     !%-----------------------------------------------------------------------------
!     !     integer, intent(in) :: indx
!     !     real(8), intent(in) :: depth
!     !     real(8), pointer :: bottomDepth(:), bottomSlope(:), breadth(:)
!     !     !%-----------------------------------------------------------------------------
!     !     breadth     => elemSGR(:,esgr_Rectangular_Triangular_TopBreadth)
!     !     bottomSlope => elemSGR(:,esgr_Rectangular_Triangular_BottomSlope)
!     !     bottomDepth => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth) 
!     !     !%-----------------------------------------------------------------------------
        
!     !     if(depth <= bottomDepth(indx)) then
!     !         outvalue = twoR * depth * sqrt(oneR + bottomSlope(indx) ** twoR)
!     !     else
!     !         outvalue = twoR * bottomDepth(indx) * sqrt(oneR + bottomSlope(indx) ** twoR) & !triangular section
!     !                     + twoR * (depth - bottomDepth(indx))                               !rectangular section

!     !     endif

!     !     if (setting%Debug%File%geometry) &
!     !         print *, 'perimeter = ' , outvalue

!     ! end function rectangular_triangular_perimeter_from_depth_singular
! !%    
! !%==========================================================================
! !%==========================================================================
! !%
!     ! real(8) function rectangular_triangular_hyddepth_from_depth_singular &
!     !     (indx, depth) result (outvalue)
!     !     !%-----------------------------------------------------------------------------
!     !     !% Description:
!     !     !% Computes hydraulic depth from known depth for rectangular_triangular cross section of 
!     !     !% a single element
!     !     !%-----------------------------------------------------------------------------   
!     !     integer, intent(in) :: indx     
!     !     real(8), intent(in) :: depth
!     !     real(8) :: topwidth, area
!     !     !real(8), pointer :: bottomdepth, fulldepth
!     !     !%-----------------------------------------------------------------------------
!     !     !bottomdepth => elemSGR(indx,esgr_Rectangular_Triangular_BottomDepth) 
!     !     !fulldepth   => elemR(indx,er_FullDepth)
!     !     !%-----------------------------------------------------------------------------  

!     !     topwidth = rect_round_topwidth_from_depth_singular (indx, indepth)
!     !     area     = rect_round_area_from_depth_singular (indx, depth)

!     !     if (depth <= setting%ZeroValue%Depth) then
!     !         !% --- empty
!     !         outvalue = setting%ZeroValue%Depth
!     !     elseif (depth >= elemR(indx,er_FullDepth))
!     !         !% --- full
!     !         outvalue = elemR(indx,er_FullHydDepth)
!     !     else
!     !         outvalue = area / topwidth
!     !     endif
        
!     !     ! if (depth <= setting%ZeroValue%Depth) then
!     !     !     !% --- empty
!     !     !     outvalue = setting%ZeroValue%Depth
!     !     ! elseif ((depth > setting%ZeroValue%Depth) .and. &
!     !     !         (depth <= bottomdepth) ) then
!     !     !     !% --- in triangular section
!     !     !     outvalue = depth / twoR
!     !     ! elseif (depth < fullDepth)
!     !     !     !% --- in rectangular section
!     !     !     outvalue = (bottomdepth / twoR) &         !triangular section
!     !     !              + (depth - bottomdepth)          !rectangular section
!     !     ! elseif (depth >= elemR(indx,er_FullDepth))
!     !     !     !% --- full
!     !     !     outvalue = elemR(indx,er_FullHydDepth)
!     !     ! endif
            
!     !     if (setting%Debug%File%geometry) &
!     !         print *, 'hyddepth = ' , outvalue

!     ! end function rectangular_triangular_hyddepth_from_depth_singular 
! !%    
! !%==========================================================================
! !%==========================================================================
! !%
!     real(8) function rectangular_triangular_hydradius_from_depth_singular &
!         (indx, depth) result (outvalue)
!         !%  
!         !%-----------------------------------------------------------------------------
!         !% Description:
!         !% Computes hydraulic radius from known depth for a rectangular_triangular cross section of
!         !% a single element 
!         !%-----------------------------------------------------------------------------
!         integer, intent(in) :: indx
!         real(8), intent(in) :: depth
!         real(8), pointer :: bottomdepth(:), breadth(:), bottomSlope(:)
!         !%-----------------------------------------------------------------------------
!         bottomSlope => elemSGR(:,esgr_Rectangular_Triangular_BottomSlope)
!         breadth     => elemSGR(:,esgr_Rectangular_Triangular_TopBreadth)
!         bottomdepth => elemSGR(:,esgr_Rectangular_Triangular_BottomDepth) 
!         !%-----------------------------------------------------------------------------
        
!         if(depth <= bottomdepth(indx)) then
!             outvalue = (bottomSlope(indx) * depth) / (twoR * sqrt(oneR + (bottomSlope(indx) ** twoR)))
!         else
!             outvalue = ((bottomSlope(indx) * bottomdepth(indx)) / (twoR * sqrt(oneR + (bottomSlope(indx) ** twoR)))) &          !triangular section
!                      + (((depth - bottomdepth(indx)) * breadth(indx)) / (twoR * (depth - bottomdepth(indx))))                   !rectangular section
!         endif

!         if (setting%Debug%File%geometry) &
!             print *, 'hydradius = ' , outvalue

!     end function rectangular_triangular_hydradius_from_depth_singular
!%    
!%==========================================================================
!%==========================================================================
!%
end module rectangular_triangular_conduit