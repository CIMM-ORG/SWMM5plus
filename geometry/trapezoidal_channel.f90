module trapezoidal_channel

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use geometry_lowlevel

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Trapezoidal channel geometry
    !%

    private

    public :: trapezoidal_depth_from_volume
    public :: trapezoidal_topwidth_from_depth
    public :: trapezoidal_perimeter_from_depth

    ! public :: trapezoidal_area_from_depth_singular
    ! public :: trapezoidal_topwidth_from_depth_singular 
    ! public :: trapezoidal_perimeter_from_depth_singular
    
    !public :: trapezoidal_hyddepth_from_depth
    !public :: trapezoidal_hyddepth_from_depth_singular
    
    !public :: trapezoidal_hydradius_from_depth_singular


    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine trapezoidal_depth_from_volume (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on open channels 
        !% Input elemPGx is pointer (already assigned) for elemPGalltm, elemPGetm or elemPGac
        !% Assumes that volume > 0 is previuosly enforced in volume computations.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: depth(:), volume(:)
            real(8), pointer :: fullvolume(:), fulldepth(:)
        !%-------------------------------------------------------------------
        !% Aliases
            depth      => elemR(:,er_Depth)
            volume     => elemR(:,er_Volume)
            fullvolume => elemR(:,er_FullVolume)
            fulldepth  => elemR(:,er_FullDepth)
        !%------------------------------------------------------------------

        if (setting%Discretization%AllowChannelOverflowTF) then
            where (volume(thisP) >= fullvolume(thisP))
                !% --- truncate depth at full depth if there is overflow
                depth(thisP) = fulldepth(thisP)
            elsewhere
                !% --- standard trapezoidal depth
                depth(thisP) = llgeo_trapezoidal_depth_from_volume_pure &
                                    (thisP, volume(thisP))
            endwhere
        else
            !% --- if overflow is NOT allowed
            where (volume(thisP) >= fullvolume(thisP))
                !% --- volume above max level is rectangular at max breadth
                depth(thisP) = llgeo_openchannel_depth_above_full_pure(thisP)
            elsewhere
                !% --- standard trapezoidal depth
                depth(thisP) = llgeo_trapezoidal_depth_from_volume_pure &
                                 (thisP, volume(thisP))
            end where
        end if
 
    end subroutine trapezoidal_depth_from_volume
!%  
!%==========================================================================
!%==========================================================================
!%
    subroutine trapezoidal_topwidth_from_depth (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the topwidth from a known depth in a trapezoidal channel
        !%------------------------------------------------------------------
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: topwidth(:), volume(:), fullvolume(:)
            real(8), pointer :: depth(:)
        !%------------------------------------------------------------------
        !% Aliases    
            topwidth    => elemR(:,er_Topwidth)
            volume      => elemR(:,er_Volume)
            fullvolume  => elemR(:,er_FullVolume)
            depth       => elemR(:,er_Depth)
        !%----------------------------------------------------------------
        if (setting%Discretization%AllowChannelOverflowTF) then 
            !% --- depth is already limited to full depth
            topwidth(thisP) = llgeo_trapezoidal_topwidth_from_depth_pure(thisP,depth(thisP))
        else
            where (volume(thisP) >= fullvolume(thisP))
                topwidth(thisP) = llgeo_openchannel_topwidth_above_full_pure(thisP)
            elsewhere
                !% --- use elemental form as depth <= fulldepth
                topwidth(thisP) = llgeo_trapezoidal_topwidth_from_depth_pure(thisP,depth(thisP))
            endwhere

        end if

        !topwidth(thisP) = breadth(thisP) + depth(thisP) * (lslope(thisP) + rslope(thisP))

    end subroutine trapezoidal_topwidth_from_depth
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine trapezoidal_perimeter_from_depth (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the perimeter from a known depth in a trapezoidal channel
        !%------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            real(8), pointer :: perimeter(:), volume(:), fullvolume(:)
            real(8), pointer :: depth(:)
        !%------------------------------------------------------------------
        !% Aliases
            perimeter  => elemR(:,er_Perimeter)
            volume     => elemR(:,er_Volume)
            fullvolume => elemR(:,er_FullVolume)
            depth      => elemR(:,er_Depth)
        !%-----------------------------------------------------------------

        if (setting%Discretization%AllowChannelOverflowTF) then 
            !% --- depth is already limited to full depth
            perimeter(thisP) = llgeo_trapezoidal_perimeter_from_depth_pure(thisP,depth(thisP))
        else
            where (volume(thisP) >= fullvolume(thisP))
                perimeter(thisP) = llgeo_openchannel_perimeter_above_full_pure(thisP)
            elsewhere
                !% --- use elemental form as depth <= fulldepth
                perimeter(thisP) = llgeo_trapezoidal_perimeter_from_depth_pure(thisP,depth(thisP))
            endwhere
        end if
        

    end subroutine trapezoidal_perimeter_from_depth
!%    
!%==========================================================================   
!%==========================================================================
!%
    ! subroutine trapezoidal_hyddepth_from_depth (elemPGx, Npack, thisCol)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Computes the hydraulic (average) depth from a known depth in a rectangular channel
    !     !%----------------------------------------------------------------
    !     integer, target, intent(in) :: elemPGx(:,:)
    !     integer, intent(in) ::  Npack, thisCol
    !     integer, pointer :: thisP(:)
    !     real(8), pointer :: hyddepth(:), depth(:), breadth(:), lslope(:), rslope(:)
    !     !%-----------------------------------------------------------------------------
    !     thisP     => elemPGx(1:Npack,thisCol) 
    !     hyddepth  => elemR(:,er_HydDepth)
    !     depth   => elemR(:,er_Depth)
    !     breadth => elemSGR(:,esgr_Trapezoidal_Breadth)
    !     lslope  => elemSGR(:,esgr_Trapezoidal_LeftSlope)
    !     rslope  => elemSGR(:,esgr_Trapezoidal_RightSlope)
    !     !%----------------------------------------------------------------------------- 

    !     hyddepth(thisP) = ((breadth(thisP) + onehalfR * (lslope(thisP) + rslope(thisP)) * &
    !                 depth(thisP)) * depth(thisP)) / (breadth(thisP) + depth(thisP) * &
    !                 (lslope(thisP) + rslope(thisP)))

    ! end subroutine trapezoidal_hyddepth_from_depth
!%    
!%========================================================================== 
!% SINGULAR 
!%==========================================================================
!%
!     pure real(8) function trapezoidal_area_from_depth_singular &
!         (indx, depth) result (outvalue)
!         !%------------------------------------------------------------------
!         !% Description:
!         !% Computes area from known depth for trapezoidal cross section of a single element
!         !% The input indx is the row index in full data 2D array.
!         !%------------------------------------------------------------------
!         integer, intent(in) :: indx
!         real(8), intent(in) :: depth
!         !%------------------------------------------------------------------
!         outvalue = ( elemSGR(indx,esgr_Trapezoidal_Breadth)                     &
!                     + onehalfR * (   elemSGR(indx,esgr_Trapezoidal_LeftSlope)   &
!                                    + elemSGR(indx,esgr_Trapezoidal_RightSlope)  &
!                                  ) * depth                                      &
!                    ) * depth

!     end function trapezoidal_area_from_depth_singular
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     pure real(8) function trapezoidal_topwidth_from_depth_singular &
!         (indx, depth) result (outvalue)
!         !%------------------------------------------------------------------
!         !% Description:
!         !% Computes the topwidth for a trapezoidal cross section of a single element
!         !%------------------------------------------------------------------
!             integer, intent(in) :: indx 
!             real(8), intent(in) :: depth
!         !%------------------------------------------------------------------
!         outvalue = elemSGR(indx,esgr_Trapezoidal_Breadth)             &
!              + depth * (  elemSGR(indx,esgr_Trapezoidal_LeftSlope)    &
!                         + elemSGR(indx,esgr_Trapezoidal_RightSlope)   &
!                        )

!     end function trapezoidal_topwidth_from_depth_singular
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     pure real(8) function trapezoidal_perimeter_from_depth_singular &
!         (indx, depth) result (outvalue)
!         !%  
!         !%------------------------------------------------------------------
!         !% Description:
!         !% Computes wetted perimeter from known depth for a trapezoidal cross section of
!         !% a single element 
!         !%-----------------------------------------------------------------
!             integer, intent(in) :: indx
!             real(8), intent(in) :: depth
!         !%-----------------------------------------------------------------------------
        
!         outvalue =  elemSGR(indx,esgr_Trapezoidal_Breadth) &
!              + depth * (  sqrt(oneR + elemSGR(indx,esgr_Trapezoidal_LeftSlope )**2) &
!                         + sqrt(oneR + elemSGR(indx,esgr_Trapezoidal_RightSlope)**2) &
!                         )

!     end function trapezoidal_perimeter_from_depth_singular
! !%    
! !%==========================================================================
! !%==========================================================================
! !%
!     ! real(8) function trapezoidal_hyddepth_from_depth_singular &
!     !     (indx, depth) result (outvalue)
!     !     !%  
!     !     !%-----------------------------------------------------------------------------
!     !     !% Description:
!     !     !% Computes hydraulic depth from known depth for trapezoidal cross section of 
!     !     !% a single element
!     !     !%-----------------------------------------------------------------------------   
!     !     integer, intent(in) :: indx 
!     !     real(8), intent(in) :: depth
!     !     real(8), pointer    ::  breadth(:), lslope(:), rslope(:)
!     !     !%-----------------------------------------------------------------------------
!     !     breadth => elemSGR(:,esgr_Trapezoidal_Breadth)
!     !     lslope  => elemSGR(:,esgr_Trapezoidal_LeftSlope)
!     !     rslope  => elemSGR(:,esgr_Trapezoidal_RightSlope)
!     !     !%-----------------------------------------------------------------------------     

!     !     outvalue = ( (breadth(indx)                                             &    
!     !                     + onehalfR * (lslope(indx) + rslope(indx)) * depth      &
!     !                  ) * depth                                                  &
!     !                  / (breadth(indx) + depth * (lslope(indx) + rslope(indx)) ) &
!     !                 )

!     ! end function trapezoidal_hyddepth_from_depth_singular
! !%    
! !%==========================================================================
! !%==========================================================================
! !%
!     real(8) function trapezoidal_hydradius_from_depth_singular (indx, depth) result (outvalue)
!         !%  
!         !%-----------------------------------------------------------------------------
!         !% Description:
!         !% Computes hydraulic radius from known depth for a trapezoidal cross section of
!         !% a single element 
!         !%-----------------------------------------------------------------------------
!         integer, intent(in) :: indx
!         real(8), intent(in) :: depth
!         real(8), pointer ::  breadth(:), lslope(:), rslope(:)
!         !%-----------------------------------------------------------------------------
!         breadth => elemSGR(:,esgr_Trapezoidal_Breadth)
!         lslope  => elemSGR(:,esgr_Trapezoidal_LeftSlope)
!         rslope  => elemSGR(:,esgr_Trapezoidal_RightSlope)
!         !%----------------------------------------------------------------------------- 
        
!         outvalue = ((breadth(indx) + onehalfR * (lslope(indx) + rslope(indx)) &
!                      *  depth) * depth                                        &
!                    )                                                          &
!                    / ( breadth(indx)                                          &
!                        + depth * (  sqrt(oneR + lslope(indx)**twoR)           &
!                                   + sqrt(oneR + rslope(indx)**twoR))          &
!                      )

!     end function trapezoidal_hydradius_from_depth_singular
!%    
!%==========================================================================
!% END MODULE
!%==========================================================================
end module trapezoidal_channel