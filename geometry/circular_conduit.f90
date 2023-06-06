module circular_conduit
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Geometry for circular closed conduit
    !%
    !%==========================================================================

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use define_xsect_tables
    use xsect_tables

    implicit none

    private

    public :: circular_depth_from_volume
    public :: circular_get_normalized_depth_from_area_analytical

    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine circular_depth_from_volume (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Only applies on conduits 
        !% Input elemPGx is pointer (already assigned) for elemPGetm
        !% Assumes that volume > 0 is enforced in volume computations.
        !%-------------------------------------------------------------------
        !% Declarations
            integer, target, intent(in) :: thisP(:)
            integer, pointer :: thisP_lookup(:), thisP_analytical(:)
            integer, pointer :: thisPA(:), thisPL(:)
            real(8), pointer :: depth(:), volume(:), length(:), AoverAfull(:)
            real(8), pointer :: YoverYfull(:), fullArea(:), fulldepth(:)
            integer, target  :: Npack_analytical, Npack_lookup
        !%---------------------------------------------------------------------
        !% Aliases
            depth       => elemR(:,er_Depth)
            volume      => elemR(:,er_Volume)
            length      => elemR(:,er_Length)
            fullArea    => elemR(:,er_FullArea)
            fulldepth   => elemR(:,er_FullDepth)
            AoverAfull  => elemR(:,er_AoverAfull)
            YoverYfull  => elemR(:,er_YoverYfull)
            
            thisP_analytical => elemI(:,ei_Temp01)
            thisP_lookup     => elemI(:,ei_Temp02)
        !%-----------------------------------------------------------------------------

        !% --- compute the relative volume, which is also the relative area
        AoverAfull(thisP) = volume(thisP) / (length(thisP) * fullArea(thisP))

        !% --- when AoverAfull <= 4%, SWMM5 uses a special function to get the
        !%     normalized depth using the central angle, theta

        !% --- pack the circular elements with AoverAfull <= 4% which will use analytical solution
        !%     from French, 1985 by using the central angle theta.
        Npack_analytical = count(AoverAfull(thisP) <= 0.04)
        if (Npack_analytical > zeroI) then

            thisP_analytical(1:Npack_analytical) = pack(thisP,AoverAfull(thisP) <= 0.04)
            thisPA => thisP_analytical(1:Npack_analytical)

            call circular_get_normalized_depth_from_area_analytical &
                (YoverYfull, AoverAfull, Npack_analytical, thisPA)
        end if

        !% --- pack the rest of the circular elements having AoverAfull > 0.04 which will use
        !%     lookup table for interpolation.
        Npack_lookup = count(AoverAfull(thisP) > 0.04)
        if (Npack_lookup > zeroI) then
            
            thisP_lookup(1:Npack_lookup) = pack(thisP,AoverAfull(thisP) > 0.04)
            thisPL => thisP_lookup(1:Npack_lookup)

            call xsect_table_lookup &
               (YoverYfull, AoverAfull, YCirc, thisPL)

        endif

        !% --- unnormalize the depth 
        depth(thisP) = YoverYfull(thisP) * fulldepth(thisP)

        !% ensure the full depth is not exceeded
        depth(thisP) = min(depth(thisP),fulldepth(thisP))

        !% --- clear the temporary storage
        if (Npack_analytical > zeroI) thisPA = nullvalueI
        if (Npack_lookup     > zeroI) thisPL = nullvalueI

    end subroutine circular_depth_from_volume
!%
!%==========================================================================  
!%==========================================================================
!%
    subroutine circular_get_normalized_depth_from_area_analytical &
        (normalizedDepth, normalizedArea, Npack, thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% find the YoverYfull from AoverAfull only when AoverAfull <= 0.04
        !% This subroutine uses the analytical derivation from French, 1985 to
        !% calculate the central angle theta. THIS SUBROUTINE IS NOT VECTORIZED
        !%
        !% this piece of the code is adapted from SWMM5-c code
        !% for more reference check SWMM Reference Manual Volume II – Hydraulics, pp-81 
        !% and SWMM5 xsect.c code
        !%-----------------------------------------------------------------------------
        !% Declarations
            real(8), intent(inout)      :: normalizedDepth(:)
            real(8), intent(in)         :: normalizedArea(:)
            integer, target, intent(in) :: Npack, thisP(:)

            integer          :: ii, jj 
            integer, pointer :: eIdx
            real(8)          :: alpha, theta, theta1, dTheta
            real(8), pointer :: pi
        !%-----------------------------------------------------------------------------
        pi => setting%Constant%pi

        do ii = 1,Npack
            eIdx   => thisP(ii)
            alpha  = normalizedArea(eIdx)

            !% --- this piece of the code is adapted from EPA-SWMM-C code
            if (alpha >= oneR) then
                normalizedDepth(eIdx) = oneR
            else if (alpha <= zeroR) then
                normalizedDepth(eIdx) = zeroR
            else if (alpha <= 1e-05) then 
                normalizedDepth(eIdx) = (((37.6911*alpha)**onethirdR)**twoR)/16.0
            else
                theta1  = 0.031715 - 12.79384 * alpha + 8.28479 * sqrt(alpha)
                theta = theta1
                do jj = 1,40
                    dTheta = - ((2.0 * pi) * alpha - theta1 + sin(theta1)) / (1.0 - cos(theta1))
                    if (dTheta > oneR) dTheta = sign(oneR,dtheta)
                    theta1 = theta1 - dTheta
                    if (abs(dTheta) <= 0.0001) then
                        theta = theta1
                        exit
                    end if    
                end do
                normalizedDepth(eIdx) = (oneR - cos(theta / twoR)) / twoR
            end if
        end do

    end subroutine circular_get_normalized_depth_from_area_analytical
!%
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module circular_conduit