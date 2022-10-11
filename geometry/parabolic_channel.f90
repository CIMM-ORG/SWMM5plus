module parabolic_channel

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use geometry_lowlevel

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% parabolic channel geometry
    !%

    private

    public :: parabolic_depth_from_volume
    public :: parabolic_topwidth_from_depth
    public :: parabolic_perimeter_from_depth

    ! public :: parabolic_area_from_depth_singular
    ! public :: parabolic_topwidth_from_depth_singular 
    ! public :: parabolic_perimeter_from_depth_singular
    
    
    !public :: parabolic_hyddepth_from_depth
    
    !public :: parabolic_hyddepth_from_depth_singular
    !public :: parabolic_hydradius_from_depth_singular

    contains

!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine parabolic_depth_from_volume (elemPGx, Npack, thisCol)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is previously enforced in volume computations.
        !%--------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: elemPGx(:,:), Npack, thisCol
            integer, pointer :: thisP(:)
            real(8), pointer :: depth(:), volume(:)
            real(8), pointer :: fulldepth(:), fullvolume(:)
        !%-------------------------------------------------------------------
            if (Npack < 1) return
        !%-------------------------------------------------------------------
            thisP      => elemPGx(1:Npack,thisCol) 
            depth      => elemR(:,er_Depth)
            volume     => elemR(:,er_Volume)
            fulldepth  => elemR(:,er_FullDepth)
            fullvolume => elemR(:,er_FullVolume)
        !%--------------------------------------------------------------------

        if (setting%Discretization%AllowChannelOverflowTF) then     
            where (volume(thisP) >= fullvolume(thisP))
                !% --- truncate depth at full depth if there is overflow
                depth(thisP) = fulldepth(thisP)
            elsewhere
                !% --- standard parabolic depth
                depth(thisP) = llgeo_parabolic_depth_from_volume_pure &
                                (thisP, volume(thisP))
            endwhere
        else
            !% --- if overflow is NOT allowed
            where (volume(thisP) >= fullvolume(thisP))
                !% --- volume above max level is rectangular at max breadth
                depth(thisP) = llgeo_openchannel_depth_above_full_pure(thisP)
            elsewhere 
                !% --- standard parabolic depth
                depth(thisP) = llgeo_parabolic_depth_from_volume_pure &
                                (thisP, volume(thisP))
            endwhere
        end if
        
    end subroutine parabolic_depth_from_volume
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine parabolic_topwidth_from_depth (elemPGx, Npack, thisCol)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a parabolic channel
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: elemPGx(:,:)
            integer, intent(in) ::  Npack, thisCol
            integer, pointer :: thisP(:)
            real(8), pointer :: topwidth(:), volume(:), fullvolume(:)
            real(8), pointer :: depth(:)
        !%-------------------------------------------------------------------
            thisP      => elemPGx(1:Npack,thisCol) 
            topwidth   => elemR(:,er_Topwidth)
            volume     => elemR(:,er_Volume)
            fullvolume => elemR(:,er_FullVolume)
            depth      => elemR(:,er_Depth)
        !%-------------------------------------------------------------------

        if (setting%Discretization%AllowChannelOverflowTF) then 
            !% --- depth is already limited to full depth
            !topwidth(thisP) = twoR * rbot(thisP) * sqrt(depth(thisP))
            topwidth(thisP) = llgeo_parabolic_topwidth_from_depth_pure(thisP,depth(thisP))
        else
            where (volume(thisP) >= fullvolume(thisP))
                topwidth(thisP) = llgeo_openchannel_topwidth_above_full_pure(thisP)
            elsewhere
                !% --- use elemental form as depth <= fulldepth
                topwidth(thisP) = llgeo_parabolic_topwidth_from_depth_pure(thisP,depth(thisP))
            endwhere
        end if
        
    end subroutine parabolic_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine parabolic_perimeter_from_depth (elemPGx, Npack, thisCol)
        !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a parabolic channel
        !% NOTE THIS DOES NOT FOLLOW THE PATTERN OF OTHER SUBROUTINES AS IT
        !% DOES NOT USE ELEMENTAL FUNCTIONS
        !%-----------------------------------------------------------------------------
        integer, target, intent(in) :: elemPGx(:,:)
        integer, intent(in) ::  Npack, thisCol
        integer, pointer :: thisP(:)
        real(8), pointer :: depth(:), perimeter(:), tempA(:)
        real(8), pointer :: volume(:), fullvolume(:)
        !%-----------------------------------------------------------------------------
        thisP      => elemPGx(1:Npack,thisCol) 
        depth      => elemR(:,er_Depth)
        perimeter  => elemR(:,er_Perimeter)
        volume     => elemR(:,er_Volume)
        fullvolume => elemR(:,er_FullVolume)
        tempA      => elemR(:,er_Temp01)
        ! rbot      => elemSGR(:,esgr_Parabolic_Radius)
        ! xx         => elemR(:, er_Temp02)
        ! tt         => elemr(:, er_Temp03)
        ! !%-----------------------------------------------------------------------------
        ! xx(thisP) = twoR * sqrt(depth(thisP)) / rbot(thisP)
        ! tt(thisP) = sqrt(oneR + xx(thisP) * xx(thisP))
        !%-----------------------------------------------------------------------------
        
        perimeter(thisP) = llgeo_parabolic_perimeter_from_depth_pure (thisP, depth(thisP))

        tempA(thisP) = zeroR

        !% --- correct for overfull channel
        !%     This could be made simpler if the function call for perimeter_from_depth
        !%     is elemental.
        if (.not. setting%Discretization%AllowChannelOverflowTF) then
            !% --- compute the corrected for all
            tempA(thisP) = llgeo_openchannel_perimeter_above_full_pure(thisP) 
            !%
            where (volume(thisP) > fullvolume(thisP))
                perimeter(thisP) = tempA(thisP)
            endwhere
        else 
            !% no changes needed if channel overflow is allowed
        endif

        !perimeter(thisP) = onehalfR * rbot(thisP) * rbot(thisP) * (xx(thisP) * tt(thisP) + log(xx(thisP) + tt(thisP)))

    end subroutine parabolic_perimeter_from_depth
!%    
!%==========================================================================  
!%==========================================================================
!%
    ! subroutine parabolic_hyddepth_from_depth (elemPGx, Npack, thisCol)
    !     !%  
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Computes the hydraulic (average) depth from a known depth in a parabolic channel
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

    !     hyddepth(thisP) = (twoR/threeR) * depth(thisP)

    ! end subroutine parabolic_hyddepth_from_depth
!%    
!%==========================================================================  
!% SINGULAR
!%==========================================================================
!%
!     pure real(8) function parabolic_area_from_depth_singular &
!         (indx, depth) result (outvalue)
!         !%------------------------------------------------------------------
!         !% Description:
!         !% Computes area from known depth for parabolic cross section of a single element
!         !% The input indx is the row index in full data 2D array.
!         !%------------------------------------------------------------------
!             integer, intent(in) :: indx
!             real(8), intent(in) :: depth
!         !%------------------------------------------------------------------

!         outvalue = fourthirdsR * elemSGR(indx,esgr_Parabolic_Radius) * depth *  sqrt(depth)

!     end function parabolic_area_from_depth_singular
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     pure real(8) function parabolic_topwidth_from_depth_singular &
!         (indx, depth) result (outvalue)
!         !%------------------------------------------------------------------
!         !% Description:
!         !% Computes the topwidth for a parabolic cross section of a single element
!         !%------------------------------------------------------------------
!             integer, intent(in) :: indx 
!             real(8), intent(in) :: depth
!         !%------------------------------------------------------------------

!         outvalue = twoR * elemSGR(indx, esgr_Parabolic_Radius) * sqrt(depth)

!     end function parabolic_topwidth_from_depth_singular
! !%
! !%==========================================================================  
! !%==========================================================================
! !%
!     pure real(8) function parabolic_perimeter_from_depth_singular &
!         (indx, depth) result (outvalue)
!         !%------------------------------------------------------------------
!         !% Description:
!         !% Computes wetted perimeter from known depth for a parabolic cross section of
!         !% a single element 
!         !%------------------------------------------------------------------
!         !%------------------------------------------------------------------
!         integer, intent(in) :: indx
!         real(8), intent(in) :: depth
!         real(8), pointer    ::  rbot
!         real(8) :: xx, tt
!         !%-----------------------------------------------------------------
!         xx = twoR * sqrt(depth) / elemSGR(indx, esgr_Parabolic_Radius)
!         tt = sqrt(oneR + xx * xx)

!         outvalue = onehalfR * (elemSGR(indx, esgr_Parabolic_Radius)**2) &
!                  * (xx * tt + log(xx + tt))

!     end function parabolic_perimeter_from_depth_singular
! !%    
! !%==========================================================================
! !%==========================================================================
! !%
!     ! real(8) function parabolic_hyddepth_from_depth_singular (indx,depth) result (outvalue)
!     !     !%  
!     !     !%-----------------------------------------------------------------------------
!     !     !% Description:
!     !     !% Computes hydraulic depth from known depth for parabolic cross section of 
!     !     !% a single element
!     !     !%-----------------------------------------------------------------------------   
!     !     integer, intent(in) :: indx   
!     !     real(8), intent(in) :: depth  
!     !     !%-----------------------------------------------------------------------------  

!     !     outvalue = (twoR/threeR) * depth

!     ! end function parabolic_hyddepth_from_depth_singular
! !%    
! !%==========================================================================
! !%==========================================================================
! !%
!     ! real(8) function parabolic_hydradius_from_depth_singular &
!     !     (indx, depth) result (outvalue)
!     !     !%------------------------------------------------------------------
!     !     !% Description:
!     !     !% Computes hydraulic radius from known depth for a parabolic cross section of
!     !     !% a single element 
!     !     !%------------------------------------------------------------------
!     !         integer, intent(in) :: indx
!     !         real(8), intent(in) :: depth
!     !         real(8) :: xx, tt, area, perimeter
!     !         real(8), pointer :: rbot
!     !         !real(8), pointer :: breadth, fulldepth, rbot
!     !     !%------------------------------------------------------------------
!     !     !breadth   => elemSGR(indx,esgr_parabolic_Breadth)
!     !     !fulldepth => elemR(indx, er_FullDepth)
!     !     rbot      => elemSGR(indx, esgr_Parabolic_Radius)
!     !     !%------------------------------------------------------------------
!     !     xx = (twoR * sqrt(depth) / rbot)
!     !     tt = sqrt(oneR + (xx * xx))
!     !     area = (fourR / threeR) * rbot * depth *  sqrt(depth)

!     !     perimeter = onehalfR * rbot * rbot * (xx * tt + log(xx + tt))
        
!     !     outvalue = area / perimeter

!     ! end function parabolic_hydradius_from_depth_singular
!%      
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module parabolic_channel