module culvert_elements

    use define_globals
    use define_keys
    use define_indexes
    use define_xsect_tables
    use define_settings, only: setting
    use geometry_lowlevel
    use xsect_tables
    use utility_crash, only: util_crashpoint

!%----------------------------------------------------------------------------- 
!% Description:
!% Computes culverts
!%----------------------------------------------------------------------------- 

    implicit none

    private

    !% --- unit conversion (BG to SI) used in FHWA culvert equations
    real(8), parameter  :: kUnit = 1.811d0

    public :: culvert_parameter_values
    public :: culvert_toplevel
    

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine culvert_parameter_values()
        !%------------------------------------------------------------------
        !% Description:
        !% Stores the culvert parameter values from EPA-SWMM culvert.c
        !% These are (in order) FORM, K, M, C, Y, SCF
        !% NOTE SCF is slope correction factor, added to table for SWMM5+
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        !% Circular concrete
        culvertValue(1,:) = (/1.0d0, 0.0098d0, 2.00d0, 0.0398d0, 0.67d0, -0.5d0 /) !% Square edge w/headwall
        culvertValue(2,:) = (/1.0d0, 0.0018d0, 2.00d0, 0.0292d0, 0.74d0, -0.5d0/)  !% Groove end w/headwall
        culvertValue(3,:) = (/1.0d0, 0.0045d0, 2.00d0, 0.0317d0, 0.69d0, -0.5d0/)  !% Groove end projecting
    
        !% Circular Corrugated Metal Pipe
        culvertValue(4,:) = (/1.0d0, 0.0078d0, 2.00d0, 0.0379d0, 0.69d0, -0.5d0/)  !% Headwall
        culvertValue(5,:) = (/1.0d0, 0.0210d0, 1.33d0, 0.0463d0, 0.75d0, +0.7d0/) !% Mitered to slope
        culvertValue(6,:) = (/1.0d0, 0.0340d0, 1.50d0, 0.0553d0, 0.54d0, -0.5d0/)  !% Projecting
    
        !% Circular Pipe, Beveled Ring Entrance
        culvertValue(7,:) = (/1.0d0, 0.0018d0, 2.50d0, 0.0300d0, 0.74d0, -0.5d0/)  !% Beveled ringd0, 45 deg bevels
        culvertValue(8,:) = (/1.0d0, 0.0018d0, 2.50d0, 0.0243d0, 0.83d0, -0.5d0/)  !% Beveled ringd0, 33.7 deg bevels
    
        !% Rectangular Box with Flared Wingwalls
        culvertValue(9,:)  = (/1.0d0, 0.026d0, 1.0d0,   0.0347d0, 0.81d0, -0.5d0/)  !% 30-75 deg. wingwall flares
        culvertValue(10,:) = (/1.0d0, 0.061d0, 0.75d0,  0.0400d0, 0.80d0, -0.5d0/)  !% 90 or 15 deg. wingwall flares
        culvertValue(11,:) = (/1.0d0, 0.061d0, 0.75d0,  0.0423d0, 0.82d0, -0.5d0/)  !% 0 deg. wingwall flares (striaght sides)
    
        !% Rectanglar Box with Flared Wingwalls & Top Edge Bevel
        culvertValue(12,:) = (/2.0d0, 0.510d0, 0.667d0, 0.0309d0, 0.80d0, -0.5d0/)  !% 45 deg. flare; 0.43D top edge bevel
        culvertValue(13,:) = (/2.0d0, 0.486d0, 0.667d0, 0.0249d0, 0.83d0, -0.5d0/)  !% 18-33.7 deg flare; 0.083D top edge bevel
    
        !% Rectangular Box; 90-deg Headwall; Chamfered or Beveled Inlet Edges
        culvertValue(14,:) = (/2.0d0, 0.515d0, 0.667d0, 0.0375d0, 0.79d0, -0.5d0/)  !% chamfered 3/4-in
        culvertValue(15,:) = (/2.0d0, 0.495d0, 0.667d0, 0.0314d0, 0.82d0, -0.5d0/)  !% beveled 1/2-in/ft at 45 deg (1:1)
        culvertValue(16,:) = (/2.0d0, 0.486d0, 0.667d0, 0.0252d0, 0.865d0, -0.5d0/) !% beveled 1-in/ft at 33.7 deg (1:1.5)
    
        !% Rectangular Box; Skewed Headwall; Chamfered or Beveled Inlet Edges
        culvertValue(17,:) = (/2.0d0, 0.545d0, 0.667d0, 0.04505d0,0.73d0, -0.5d0/)  !% 3/4" chamfered edged0, 45 deg skewed headwall
        culvertValue(18,:) = (/2.0d0, 0.533d0, 0.667d0, 0.0425d0, 0.705d0, -0.5d0/) !% 3/4" chamfered edged0, 30 deg skewed headwall
        culvertValue(19,:) = (/2.0d0, 0.522d0, 0.667d0, 0.0402d0, 0.68d0, -0.5d0/)  !% 3/4" chamfered edged0, 15 deg skewed headwall
        culvertValue(20,:) = (/2.0d0, 0.498d0, 0.667d0, 0.0327d0, 0.75d0, -0.5d0/)  !% 45 deg beveled edged0, 10-45 deg skewed headwall
    
        !% Rectangular box, Non-offset Flared Wingwalls; 3/4" Chamfer at Top of Inlet
        culvertValue(21,:) = (/2.0d0, 0.497d0, 0.667d0, 0.0339d0, 0.803d0, -0.5d0/) !% 45 deg (1:1) wingwall flare
        culvertValue(22,:) = (/2.0d0, 0.493d0, 0.667d0, 0.0361d0, 0.806d0, -0.5d0/) !% 18.4 deg (3:1) wingwall flare
        culvertValue(23,:) = (/2.0d0, 0.495d0, 0.667d0, 0.0386d0, 0.71d0, -0.5d0/)  !% 18.4 deg (3:1) wingwall flared0, 30 deg inlet skew
    
        !% Rectangular box, Offset Flared Wingwallsd0, Beveled Edge at Inlet Top
        culvertValue(24,:) = (/2.0d0, 0.497d0, 0.667d0, 0.0302d0, 0.835d0, -0.5d0/)  !% 45 deg (1:1) flared0, 0.042D top edge bevel
        culvertValue(25,:) = (/2.0d0, 0.495d0, 0.667d0, 0.0252d0, 0.881d0, -0.5d0/)  !% 33.7 deg (1.5:1) flared0, 0.083D top edge bevel
        culvertValue(26,:) = (/2.0d0, 0.493d0, 0.667d0, 0.0227d0, 0.887d0, -0.5d0/)  !% 18.4 deg (3:1) flared0, 0.083D top edge bevel
    
        !%  Corrugated Metal Box
        culvertValue(27,:) = (/1.0d0, 0.0083d0, 2.00d0, 0.0379d0, 0.69d0, -0.5d0/)  !% 90 deg headwall
        culvertValue(28,:) = (/1.0d0, 0.0145d0, 1.75d0, 0.0419d0, 0.64d0, -0.5d0/)  !% Thick wall projecting
        culvertValue(29,:) = (/1.0d0, 0.0340d0, 1.50d0, 0.0496d0, 0.57d0, -0.5d0/)  !% Thin wall projecting
    
        !%  Horizontal Ellipse Concrete
        culvertValue(30,:) = (/1.0d0, 0.0100d0, 2.00d0, 0.0398d0, 0.67d0, -0.5d0/)  !% Square edge w/headwall
        culvertValue(31,:) = (/1.0d0, 0.0018d0, 2.50d0, 0.0292d0, 0.74d0, -0.5d0/)  !% Grooved end w/headwall
        culvertValue(32,:) = (/1.0d0, 0.0045d0, 2.00d0, 0.0317d0, 0.69d0, -0.5d0/)  !% Grooved end projecting
    
        !%  Vertical Ellipse Concrete
        culvertValue(33,:) = (/1.0d0, 0.0100d0, 2.00d0, 0.0398d0, 0.67d0, -0.5d0/)  !% Square edge w/headwall
        culvertValue(34,:) = (/1.0d0, 0.0018d0, 2.50d0, 0.0292d0, 0.74d0, -0.5d0/)  !% Grooved end w/headwall
        culvertValue(35,:) = (/1.0d0, 0.0095d0, 2.00d0, 0.0317d0, 0.69d0, -0.5d0/)  !% Grooved end projecting
    
        !%  Pipe Arch, 18" Corner Radiusd0, Corrugated Metal
        culvertValue(36,:) = (/1.0d0, 0.0083d0, 2.00d0, 0.0379d0, 0.69d0, -0.5d0/)  !% 90 deg headwall
        culvertValue(37,:) = (/1.0d0, 0.0300d0, 1.00d0, 0.0463d0, 0.75d0, +0.7d0/)  !% Mitered to slope
        culvertValue(38,:) = (/1.0d0, 0.0340d0, 1.50d0, 0.0496d0, 0.57d0, -0.5d0/)  !% Projecting
    
        !%  Pipe Arch, 18" Corner Radiusd0, Corrugated Metal
        culvertValue(39,:) = (/1.0d0, 0.0300d0, 1.50d0, 0.0496d0, 0.57d0, -0.5d0/)  !% Projecting
        culvertValue(40,:) = (/1.0d0, 0.0088d0, 2.00d0, 0.0368d0, 0.68d0, -0.5d0/)  !% No bevels
        culvertValue(41,:) = (/1.0d0, 0.0030d0, 2.00d0, 0.0269d0, 0.77d0, -0.5d0/)  !% 33.7 deg bevels
    
        !%  Pipe Arch, 31" Corner Radiusd0, Corrugated Metal
        culvertValue(42,:) = (/1.0d0, 0.0300d0, 1.50d0, 0.0496d0, 0.57d0, -0.5d0/)  !% Projecting
        culvertValue(43,:) = (/1.0d0, 0.0088d0, 2.00d0, 0.0368d0, 0.68d0, -0.5d0/)  !% No bevels
        culvertValue(44,:) = (/1.0d0, 0.0030d0, 2.00d0, 0.0269d0, 0.77d0, -0.5d0/)  !% 33.7 deg. bevels
    
        !%  Arch, Corrugated Metal
        culvertValue(45,:) = (/1.0d0, 0.0083d0, 2.00d0, 0.0379d0, 0.69d0, -0.5d0/)  !% 90 deg headwall
        culvertValue(46,:) = (/1.0d0, 0.0300d0, 1.00d0, 0.0473d0, 0.75d0, +0.7d0/)  !% Mitered to slope                     !% (5.1.013)
        culvertValue(47,:) = (/1.0d0, 0.0340d0, 1.50d0, 0.0496d0, 0.57d0, -0.5d0/)  !% Thin wall projecting
    
        !%  Circular Culvert
        culvertValue(48,:) = (/2.0d0, 0.534d0, 0.555d0, 0.0196d0, 0.90d0, -0.5d0/)  !% Smooth tapered inlet throat
        culvertValue(49,:) = (/2.0d0, 0.519d0, 0.640d0, 0.0210d0, 0.90d0, -0.5d0/)  !% Rough tapered inlet throat
    
        !%  Elliptical Inlet Face
        culvertValue(50,:) = (/2.0d0, 0.536d0, 0.622d0, 0.0368d0, 0.83d0, -0.5d0/)  !% Tapered inletd0, beveled edges
        culvertValue(51,:) = (/2.0d0, 0.5035d0,0.719d0, 0.0478d0, 0.80d0, -0.5d0/)  !% Tapered inletd0, square edges
        culvertValue(52,:) = (/2.0d0, 0.547d0, 0.800d0, 0.0598d0, 0.75d0, -0.5d0/)  !% Tapered inletd0, thin edge projecting
    
        !%  Rectangular
        culvertValue(53,:) = (/2.0d0, 0.475d0, 0.667d0, 0.0179d0, 0.97d0, -0.5d0/)  !% Tapered inlet throat
    
        !%  Rectangular Concrete
        culvertValue(54,:) = (/2.0d0, 0.560d0, 0.667d0, 0.0446d0, 0.85d0, -0.5d0/)  !% Side taperedd0, less favorable edges
        culvertValue(55,:) = (/2.0d0, 0.560d0, 0.667d0, 0.0378d0, 0.87d0, -0.5d0/)  !% Side taperedd0, more favorable edges
    
        !%  Rectangular Concrete
        culvertValue(56,:) = (/2.0d0, 0.500d0, 0.667d0, 0.0446d0, 0.65d0, -0.5d0/) !% Slope taperedd0, less favorable edges
        culvertValue(57,:) = (/2.0d0, 0.500d0, 0.667d0, 0.0378d0, 0.71d0, -0.5d0/)  !% Slope taperedd0, more favorable edges
    

    end subroutine culvert_parameter_values
!%
!%==========================================================================
!%==========================================================================    
!%
    subroutine culvert_toplevel ()
        !%------------------------------------------------------------------
        !% Description:
        !% checks for and enforces culvert behaviors
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: thisCol, npack, eIn, eOut, EqForm, geoType, fInlet
            integer, pointer :: fUp(:), fDn(:), thisE(:)
            real(8), pointer :: fHeadD(:), fHeadU(:),  Flowrate(:)
            real(8), pointer :: Zbtm(:), Atable(:), Ttable(:)
            integer          :: ii

            logical :: isConverged     = .false.
            logical :: isReversedFlow  = .false.
            logical :: isSubmerged     = .false.
            logical :: isTabular       = .false.
            logical :: isTransition    = .false.
            logical :: isInconsistentFlow = .false.
            real(8) :: QIC, Dinlet, PsiOut, H1s, H1u

        !%------------------------------------------------------------------
        !% Aliases
        
            thisCol  => col_elemP(ep_Culvert_Inlet)
            npack    => npack_elemP(ep_Culvert_Inlet)
            if (npack < 1) return
            
            thisE    => elemP(:,thisCol)  !% culvert inlet element
            fUp      => elemI(:,ei_Mface_uL) !% face upstream of inlet
            fDn      => elemI(:,ei_Mface_dL) !% face downstream of outlet
            fHeadD   => faceR(:,fr_Head_d)
            fHeadU   => faceR(:,fr_Head_u)
            Flowrate => elemR(:,er_Flowrate)
            Zbtm     => elemR(:,er_Zbottom)
            

            Atable => ACirc !% Dummy to prevent unallocated pointer
            Ttable => TCirc !% Dummy to prevent unallocated pointer

        !%------------------------------------------------------------------

        !print *, 'starting culvert_toplevel'

        !% --- cycle through the culverts
        do ii=1,npack
            !print *, 'ii ',ii, npack
            !% --- inlet and outlet
            eIn    => thisE(ii)
            fInlet => elemI(thisE(ii),ei_Mface_uL)
            eOut   => elemSI(eIn,esi_Conduit_Culvert_OutletID)
            EqForm => elemSI(eIn,esi_Conduit_Culvert_EquationForm)
            isConverged = .false.
            isInconsistentFlow = .false.

            !% --- geometry type
            geoType => elemI(eIn,ei_geometryType)

            !print *, 'calling culvert table pointers'
            !% --- get the correct table pointers for the geometry types
            call culvert_table_pointers(geoType, isTabular, Atable, Ttable)
            
            !print *, 'calling culvert_isReversedFlow'
            !% --- check for flow reversal (upstream flow) and reset eIn=eOut if 
            !%     reversed flow occurs.
            isReversedFlow = culvert_isReversedFlow(eIn, eOut, fInlet, isInconsistentFlow)

            !% --- where flow is into culvert from both ends or out of culvert from b
            !%     both ends, do not use culvert flow limitation
            if (isInconsistentFlow) cycle
                
            !print *, 'calling culvert_inlet_depth'
            !% --- effective inlet flow depth
            Dinlet = culvert_inlet_depth(eIn,isReversedFlow)

            !print *, 'Dinlet , dculvert ',Dinlet, elemR(eIn,er_Depth)

            !% --- if negative inlet depth then there is no QIC
            if (Dinlet .le. zeroR) cycle

            !print *, 'calling culvert_isSubmerged'
            !% --- check for submerged culvert and store H1s
            call culvert_isSubmerged (eIn, isReversedFlow, isSubmerged, H1s)
            
            !% --- if submerged, then get QIC directly
            if (isSubmerged) then
                !print *, 'isSubmerged ',isSubmerged
                QIC =  culvert_QIC_submerged_eq (eIn, Dinlet, isReversedFlow)
            else
                !print *, 'calling culvert_istransition'
                !% --- check for transition status (between submerged and unsubmerged)
                !%     and store H1u
                call culvert_isTransition (eIn,Dinlet,isTransition,H1u)

                !if (isTransition) print *, 'isTransition ',isTransition
                !print *, 'EqForm       ',EqForm

                !% --- compute Psi = D_c/D_full (only for unsumberged form 1)
                if (EqForm == 1) then      
                    !% --- get the normalized critical depth
                    PsiOut =  culvert_Psi_unsubmerged_form1 (isConverged, isReversedFlow, &
                                isTabular, geoType, eIn, Dinlet, Atable, Ttable)
                else 
                    !% --- PsiOut is a dummy that will not affect the solution
                    PsiOut = zeroR
                end if   

                !% --- compute QIC for either form 1 or 2 unsubmerged
                !%     Note this assigns a negative value for a reversed flow
                QIC = culvert_QIC_unsubmerged &
                        (isTabular,isReversedFlow, EqForm, eIn, geoType, PsiOut, Dinlet, &
                         Atable, TTable)
         
                !print *, 'out of culvert QIC_unsubmerged'
                !% ---  handle transition
                if (isTransition) then 
                    !print *, 'calling culvert_QIC_tranistion'
                    QIC = culvert_QIC_transition (eIn, Dinlet, QIC, H1s, H1u, &
                                                  isTransition, isReversedFlow)
                end if

            endif

            ! print *, 'QIC, flowrate: ',QIC, elemR(eIn,er_Flowrate)

            !% --- reset the flowrate on element and upstrem face
            !%     if inlet control is less than time-advance flowrate
            if (abs(QIC) < abs(faceR(fInlet,fr_Flowrate))) then 
                elemR(eIn,er_Flowrate) = QIC
                faceR(fInlet,fr_Flowrate) = QIC
            end if
            
        end do


    end subroutine culvert_toplevel
!%
!%==========================================================================  
!%==========================================================================  
!%
    subroutine culvert_table_pointers (geoType, isTabular, Atable, Ttable)    
        !%------------------------------------------------------------------
        !% Description
        !% provides pointers to the correct tables for the geometry type
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in)    :: geoType
            logical, intent(inout) :: isTabular
            real(8), pointer, intent(inout) :: Atable(:), Ttable(:)
        !%------------------------------------------------------------------

        select case (geoType)

        case (arch, basket_handle, catenary, circular, custom, eggshaped, gothic, &
            horiz_ellipse, horseshoe, semi_circular, semi_elliptical, vert_ellipse)
            !% --- set the table for tabular geometry
            isTabular = .true.
            select case (geoType)
            case (arch)
                Atable => AArch 
                Ttable => TArch

            case (basket_handle)
                Atable => ABasketHandle
                Ttable => TBasketHandle

            case (catenary)
                Atable => ACatenary
                Ttable => TCatenary

            case (circular)
                Atable => ACirc
                Ttable => TCirc

            case (custom)
                !Atable => A
                !Ttable => T
                print *, 'CODE ERROR: Custom conduit cross sections not completed'
                call util_crashpoint(6229873)

            case (eggshaped)
                Atable => AEgg
                Ttable => TEgg

            case (gothic)
                Atable => AGothic
                Ttable => TGothic

            case (horiz_ellipse)
                Atable => AHorizEllip
                Ttable => THorizEllip

            case (horseshoe)
                Atable => AHorseShoe
                Ttable => THorseShoe

            case (semi_circular)
                Atable => ASemiCircular
                Ttable => TSemiCircular

            case (semi_elliptical)
                Atable => ASemiEllip
                Ttable => TSemiEllip

            case (vert_ellipse)
                Atable => AVertEllip
                Ttable => TVertEllip

            case default
                print *, 'CODE ERROR: unexpected case default'
                call util_crashpoint(7298733)
            end select

        case (filled_circular)
            isTabular = .false.
            Atable => ACirc
            Ttable => TCirc

        case (mod_basket)
            isTabular = .false.
            Atable => ACirc  !% Dummy to prevent unallocated pointer
            Ttable => TCirc  !% Dummy to prevent unallcoated pointer

        case (rectangular_closed)
            isTabular = .false.
            Atable => ACirc  !% Dummy to prevent unallocated pointer
            Ttable => TCirc  !% Dummy to prevent unallcoated pointer

        case (rect_round)
            isTabular = .false.
            Atable => ACirc  !% Dummy to prevent unallocated pointer
            Ttable => TCirc  !% Dummy to prevent unallcoated pointer

        case (rect_triang)
            isTabular = .false.
            Atable => ACirc  !% Dummy to prevent unallocated pointer
            Ttable => TCirc  !% Dummy to prevent unallcoated pointer

        case default
            print *, 'CODE ERROR: unexpected case default'
            call util_crashpoint(6098723)

        end select

    end subroutine culvert_table_pointers    
!%
!%==========================================================================  
!%==========================================================================  
!%
    logical function culvert_isReversedFlow &
            (eIn, eOut, fInlet, isInconsistentFlow) result (isReversedFlow)
        !%------------------------------------------------------------------
        !% Description:
        !% Determines whether the flow in culvert is in the nominal downstream
        !% direction, reversed, or inconsistent
        !%------------------------------------------------------------------
        !% Declarations:
            logical, intent(inout)          :: isInconsistentFlow
            integer, pointer, intent(inout) :: eIn, fInlet
            integer, pointer, intent(in)    :: eOut
            real(8), pointer                :: Flowrate(:)
        !%------------------------------------------------------------------
        !% Aliases:
            Flowrate => elemR(:,er_Flowrate)
        !%------------------------------------------------------------------

        !% --- check for reversed flow direction
        if ((Flowrate(eIn) < zeroR) .and. (Flowrate(eOut) < zeroR)) then 
            !% --- upstream (reversed) flow
            isReversedFlow = .true.
            isInconsistentFlow = .false.
            !% --- make the inlet point at the outlet
            eIn => eOut
            !% --- make the face point at the outlet downstream face
            fInlet => elemI(eOut,ei_MFace_dL)
        elseif ((Flowrate(eIn) > zeroR) .and. (Flowrate(eOut) > zeroR)) then  
            !% --- downstream (normal direction) flow
            isReversedFlow = .false.
            isInconsistentFlow = .false.
        else
            !% --- inconsistent flow directions, ignore culvert equations
            isReversedFlow = .false.
            isInconsistentFlow = .true.
        end if

    end function culvert_isReversedFlow
!%
!%==========================================================================  
!%==========================================================================  
!%
    real(8) function culvert_inlet_depth (eIn,isReversedFlow) result (outvalue)
        !%------------------------------------------------------------------
        !% Description
        !% Sets the inlet depth as difference between the head on the
        !% upstream face and the Zbottom of the first element in the culvert
        !% Note this assumes that eIn is the nominal outlet end of culvert if
        !% the culvert flow is reversed.
        !% Note for a sediment-filled culvert this is the depth WITHOUT
        !% considering the sediment.
        !%------------------------------------------------------------------
        !% Declarations: 
            logical, intent(in) :: isReversedFlow
            integer, intent(in) :: eIn  !% inlet  elements of culvert
            integer, pointer    :: fup
        !%------------------------------------------------------------------   
        if (isReversedFlow) then 
            fup      => elemI(eIn,ei_Mface_dL)
            outvalue = faceR(fup,fr_Head_u) - elemR(eIn,er_Zbottom)
        else
            fup      => elemI(eIn,ei_Mface_uL)
            outvalue = faceR(fup,fr_Head_d) - elemR(eIn,er_Zbottom)   
        endif
        
    
    end function culvert_inlet_depth
!%
!%==========================================================================  
!%==========================================================================  
!%   
    subroutine culvert_isSubmerged &
            (eIn, isReversedFlow, isSubmerged, H1s )
        !%------------------------------------------------------------------
        !% Description
        !% determines whether the culvert meets the submerged flow condition
        !% of inlet head greater than crown of culvert or equation 7-48
        !% of the SWMM5 Reference Manual Vol II Hydraulics.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: eIn
            logical, intent(in)    :: isReversedFlow
            logical, intent(inout) :: isSubmerged
            real(8), intent(inout) :: H1s
            real(8), pointer :: Dfull,  Zbtm, Zcrown, slope
            real(8), pointer :: Ccoef, Ycoef, SCF, Hup
        !%------------------------------------------------------------------
        !% Aliases
            Dfull  =>  elemR(eIn,er_FullDepth)
            if (.not. isReversedFlow) then
                Hup    =>  faceR(elemI(eIn,ei_Mface_uL),fr_Head_d)
            else
                Hup    =>  faceR(elemI(eIn,ei_Mface_dL),fr_Head_u)
            end if
            Zbtm   =>  elemR(eIn,er_Zbottom)
            Zcrown =>  elemR(eIn,er_Zcrown)
            slope  =>  elemR(eIn,er_BottomSlope)
            Ccoef  => elemSR(eIn,esr_Conduit_Culvert_C)
            Ycoef  => elemSR(eIn,esr_Conduit_Culvert_Y)
            SCF    => elemSR(eIn,esr_Conduit_Culvert_SCF)
        !%------------------------------------------------------------------
        if (.not. isReversedFlow) then
            H1s = (ZBtm + Dfull * (16.d0 * Ccoef + Ycoef + SCF*slope))
            if ( (Hup > Zcrown) .or. (Hup > H1s) ) then 
                isSubmerged = .true.
            else 
                isSubmerged = .false.
            end if
        else 
            !% --- reverse flow neglects slope term
            H1s = (ZBtm + Dfull * (16.d0 * Ccoef + Ycoef))
            if ( (Hup > Zcrown) .or. (Hup > H1s)) then 
                isSubmerged = .true.
            else 
                isSubmerged = .false.
            end if   
        end if

    end subroutine culvert_isSubmerged
!%
!%==========================================================================     
!%==========================================================================  
!% 
    subroutine culvert_isTransition (eIn, Dinlet, isTransition, H1u)     
        !%------------------------------------------------------------------
        !% Description
        !% Checks whether an unsubmerged culvert is in the transition zone
        !% Should only be called on nominally unsubmerged elements (i.e.
        !% that have already passed the submergence test).
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in)    :: eIn
            real(8), intent(in)    :: Dinlet
            real(8), intent(inout) :: H1u
            logical, intent(inout) :: isTransition
            real(8), parameter  :: dfactor = 0.95d0 !% from EPA-SWMM Hydraulics manual eq 7-49
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        !% --- store the cutoff depth for the transition regime
        H1u =  dfactor * elemR(eIn,er_FullDepth)

        !% --- set logical if in transition
        if (Dinlet .ge. H1u) then
            isTransition = .true.
        else
            isTransition = .false.
        end if

        !% --- output is cutoff head for transition
        H1u = H1u + elemR(eIn,er_Zbottom)

    end subroutine culvert_isTransition
!%
!%==========================================================================  
!%==========================================================================  
!%
    real(8) function culvert_QIC_submerged_eq &
            (eIn, Dinlet, isReversedFlow) result (outvalue)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the inlet submerged form of the culvert equations
        !% Output is the inlet-controlled flowrate
        !% eIn is the index of the true inlet to the culvert (i.e., the
        !% nominal outlet in reversed flow)
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIn
            real(8), intent(in) :: Dinlet
            logical, intent(in) :: isReversedFlow
            real(8), pointer :: Afull, Dfull, slope
            real(8), pointer :: Ccoef, Ycoef, SCF
        !%------------------------------------------------------------------
            Afull  =>  elemR(eIn,er_FullArea)
            Dfull  =>  elemR(eIn,er_FullDepth)
            slope  =>  elemR(eIn,er_BottomSlope)
            Ccoef  => elemSR(eIn,esr_Conduit_Culvert_C)
            Ycoef  => elemSR(eIn,esr_Conduit_Culvert_Y)
            SCF    => elemSR(eIn,esr_Conduit_Culvert_SCF)
        !%------------------------------------------------------------------

        if (.not. isReversedFlow) then
            !% --- conventional downstream flow
            outvalue = ( Afull * sqrt(Dfull) / kUnit )                  &
                    * sqrt(                                             &
                            ( ( Dinlet / Dfull) - Ycoef - SCF * slope ) &
                            / Ccoef                                     &
                        ) 
        else
            !% --- reversed flow in culvert (neglects slope term)
            !%     This applies the negative for the flowrate direction
            outvalue = - ( Afull * sqrt(Dfull) / kUnit )                &
                    * sqrt(                                             &
                            ( ( Dinlet / Dfull) - Ycoef )               &
                            / Ccoef                                     &
                        ) 
        end if


    end function culvert_QIC_submerged_eq
!%
!%==========================================================================     
!%==========================================================================  
!%
    real(8) function culvert_Psi_unsubmerged_form1              &
            (isConverged, isReversedFlow, isTabular,        &
             geoType, eIn, Dinlet, Atable, Ttable ) result (outvalue)
        !%------------------------------------------------------------------
        !% Description
        !% Solves form1 unsubmerged equation using SWMM5+ approach
        !%------------------------------------------------------------------
        !% Declarations:
            logical, intent(inout) :: isConverged
            logical, intent(in) :: isReversedFlow, isTabular
            integer, intent(in) :: eIn, geoType
            real(8), intent(in) :: Dinlet, Atable(:), Ttable(:)
            
            real(8), pointer :: sedimentDepth, totalPipeDiameter

            real(8) :: Bhat, Khat, Dhat, Delta, halfM
            real(8) :: Psi(4), resid(4), Omega(3)

            integer :: ii

            real(8), parameter :: epsConverged = 0.0003d0
            real(8), parameter :: smalldepth = 0.001d0
        !%------------------------------------------------------------------
        !% Aliases:
            if (geoType == filled_circular) then
                sedimentDepth     =>   elemR(eIn,er_SedimentDepth)
                totalPipeDiameter => elemSGR(eIn,esgr_Filled_Circular_TotalPipeDiameter)
            end if
        !%------------------------------------------------------------------    
        !% --- STEP 1: get the culvert constants
        halfM = elemSR(eIn,esr_Conduit_Culvert_M) * onehalfR
        Bhat = culvert_Bhat (eIn)
        Khat = culvert_Khat (eIn, Bhat)
        Dhat = culvert_Dhat (eIn, Dinlet, isReversedFlow)

        ! print *, ' '
        ! print *, 'halfM ',halfM
        ! print *, 'Bhat  ',Bhat 
        ! print *, 'Khat  ',Khat 
        ! print *, 'Dhat  ',Dhat
        !print *, 'smalldepth ',smalldepth

        !% --- STEP 2: initial values
        !%     Psi is the nondimensionalized depth (solution variable)
        Psi(1) = smalldepth
        Psi(2) = min(oneR, Dinlet/elemR(eIn,er_FullDepth) )

        !% --- for very small Dinlet/Dconduit, set Q = 0
        if (Psi(2) < Psi(1)) then
            outvalue = zeroR
            return
        end if

        ! print *, 'Dinlet FullDepth ',Dinlet, elemR(eIn,er_FullDepth)

        ! print *, 'Psi(1:2) ', Psi(1), Psi(2)

        if (geoType == filled_circular) then 
            !% --- Psi must be based on total pipe diameter, not full depth
            !%     This is needed for later table lookups
            Psi(1) = (sedimentDepth + smalldepth) / totalPipeDiameter 
            Psi(2) = min(oneR, (Dinlet + sedimentDepth) / totalPipeDiameter)
        end if
        Psi(3) = onehalfR * (Psi(1) + Psi(2))

        ! print *, 'Psi(3) ',Psi(3)

        !% --- STEP 3: define Delta
        Delta = abs(Psi(1)-Psi(2))

        ! print *, 'Delta ',Delta

        !% --- STEPS 4, 5, 6: define gamma, phi, omega, and residual
        do ii=1,2
            Omega(ii) = culvert_omega (isTabular, geoType, eIn, Psi(ii), Atable, Ttable)
            resid(ii) = culvert_residual_form1 (Psi(ii), Omega(ii), halfM, Dhat, Khat, Bhat)
        end do

        ! print *, 'resid(1:2) ',resid(1), resid(2)
        ! print *, 'Omega(:)   ',Omega(1), Omega(2)


        !% --- STEP 7: check to see if we're done
        if ( abs(resid(1)) < epsConverged ) then
            outvalue = Psi(1)
            return
        end if
        if ( abs(resid(2)) < epsConverged ) then 
            outvalue = Psi(2)
            return
        end if

        !% --- STEP 8: check for root being constrained by residuals
        !%     they must be opposite sign
        if ( resid(1) * resid(2) > zeroR ) then 
            outvalue = nullvalueR  !% failure
            print *, resid(1), resid(2)
            print *, 'CODE ERROR: failure of residual constraint in culvert_Psi_unsubmerged_form1 '
            call util_crashpoint(62987662)
            !return
        end if

        isConverged = .false.

        do while (.not. isConverged)
            !% --- STEPS 4, 5, and 6: Update for Psi(3)
            Omega(3) = culvert_omega (isTabular, geoType, eIn, Psi(3), Atable, Ttable)
            resid(3) = culvert_residual_form1 (Psi(3), Omega(3), halfM, Dhat, Khat, Bhat)

            !print *, 'iteration Step 6 resid(3) ', resid(3)

            !% --- STEP 7: check for solution
            if ( abs(resid(3)) < epsConverged ) then 
                outvalue = Psi(3)
                isConverged = .true.
                return
            end if

            !% --- STEP 9: Compute Psi(4) and Omega(4)
            Psi(4) = Psi(3)                                                     &
                + (resid(3) * (Psi(3) - Psi(1)) * sign(oneR,resid(1)-resid(2))) &
                / sqrt( (resid(3)**2) - resid(1) * resid(2)  )

            !print *, 'Psi(4) ',Psi(4)   

            !% --- STEP 10: check for solution
            if (abs(Psi(4)-Psi(3)) .le. epsConverged) then 
                outvalue = Psi(3)
                isConverged = .true.
                return
            end if

            !% --- STEP 12,13: Compute Omega(4) and resid(4)
            Omega(4) = culvert_omega (isTabular, geoType, eIn, Psi(4), Atable, Ttable)  
            resid(4) = culvert_residual_form1 (Psi(4), Omega(4), halfM, Dhat, Khat, Bhat)

            !% --- STEP 13a: check for solution
            if ( abs(resid(4)) < epsConverged) then 
                outvalue = Psi(4)
                isConverged = .true.
                return
            end if

            ! print *, 'step 12, 13',Omega(4), resid(4)
            ! print *, 'resid(3),resid(4)', resid(3), resid(4)

            !% --- STEPS 14,15,16,17,18
            if (resid(3) * resid(4) < zeroR) then 
                Psi(1)   = Psi(3)
                Psi(2)   = Psi(4)
                resid(1) = resid(3)
                resid(2) = resid(4)
            elseif (resid(1) * resid(4) < zeroR) then 
                Psi(2)   = Psi(4)
                resid(2) = resid(4)
            elseif (resid(2) * resid(4) < zeroR) then 
                Psi(1)   = Psi(4)
                resid(1) = resid(4)
            else
                !% FAILURE -- should not reach this point
                outvalue = nullvalueR ! 
                print *, 'CODE ERROR: diverging solution in culvert'
                call util_crashpoint(598734)
                !return
            end if

            !% --- STEP 19: update Psi(3)
            Psi(3) = onehalfR * (Psi(1) + Psi(2))

            !% --- STEP 20,21,22: check if finished and update Delta
            if (abs(Psi(1)- Psi(2)) .le. epsConverged) then 
                outvalue = Psi(3)
                isConverged = .true.
                return
            elseif ((abs(Psi(1)- Psi(2)) > Delta)) then
                !% FAILURE -- diverging solution
                outvalue = nullvalueR
                print *, 'CODE ERROR: Diverging solution in culvert'
                call util_crashpoint(7798743)
                !return
            else 
                Delta = abs(Psi(1) - Psi(2))
                !% continue do loop
            end if

            !% TESTING TO PREVENT INFINITE LOOP
            !isConverged = .true.

        end do
        
    end function culvert_Psi_unsubmerged_form1 
!%
!%==========================================================================  
!%==========================================================================  
!%   
    real(8) function culvert_QIC_unsubmerged                       & 
            (isTabular, isReversedFlow, EqnForm, eIn, geoType, Psi, Dinlet, &
             Atable, Ttable ) result (outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes inlet-constrained flowrate for culvert
        !%------------------------------------------------------------------
        !% Declarations
            logical, intent(in) :: isTabular, isReversedFlow
            integer, intent(in) :: EqnForm, eIn, geoType
            real(8), intent(in) :: Psi, Dinlet, Atable(:), Ttable(:)

            real(8), pointer    :: grav, Afull, Dfull, Tmax, sedimentArea
            real(8)             :: Tvalue, Avalue
        !%------------------------------------------------------------------
        !% Aliases:
            Afull  =>  elemR(eIn,er_FullArea)
            Dfull  =>  elemR(eIn,er_FullDepth)
            Tmax   =>  elemR(eIn,er_BreadthMax)
            grav   =>  setting%Constant%gravity
            sedimentArea =>  elemSGR(eIn,esgr_Filled_Circular_bottomArea)
        !%------------------------------------------------------------------
        select case (EqnForm)
        case (1)
            if (isTabular) then
                !% --- QIC based on normalized value lookups for unsumberged EqForm=1
                !%     where Psi has been computed
                outvalue = sqrt(                                                       &
                                ( grav * (Afull**3)                                    &
                                * (xsect_table_lookup_singular(Psi,Atable)**3)         &
                                )                                                      &
                                / ( Tmax * xsect_table_lookup_singular(Psi,Ttable) )   &
                                )
            else 
                select case (geoType)
                !% --- QIC based on computations for physical values    
                case(filled_circular)
                    !% --- requires Psi to be based on total pipe diameter
                    Tvalue = Tmax  * xsect_table_lookup_singular(Psi,Ttable)
                    Avalue = Afull * xsect_table_lookup_singular(Psi,Atable) - sedimentArea
                    outvalue = sqrt(grav * (Avalue**3) / Tvalue)

                case (mod_basket)
                    Tvalue = llgeo_mod_basket_topwidth_from_depth_singular (eIn, Psi*Dfull, setting%ZeroValue%Topwidth)
                    Avalue = llgeo_mod_basket_area_from_depth_singular     (eIn, Psi*Dfull, setting%ZeroValue%Area)
                    outvalue = sqrt( grav * (Avalue**3) / Tvalue )

                case (rectangular_closed)
                    outvalue = sqrt( grav * ((Psi*Dfull*Tmax)**3) / Tmax )

                case (rect_round)
                    !% --- using physical values
                    Tvalue = llgeo_rect_round_topwidth_from_depth_singular (eIn, Psi*Dfull, setting%ZeroValue%Topwidth)
                    Avalue = llgeo_rect_round_area_from_depth_singular     (eIn, Psi*Dfull, setting%ZeroValue%Area)
                    outvalue = sqrt( grav * (Avalue**3) / Tvalue )

                case (rect_triang)
                    !% --- using physical values
                    Tvalue = llgeo_rectangular_triangular_topwidth_from_depth_singular (eIn, Psi*Dfull, setting%ZeroValue%Topwidth)
                    Avalue = llgeo_rectangular_triangular_area_from_depth_singular     (eIn, Psi*Dfull, setting%ZeroValue%Area)
                    outvalue = sqrt( grav * (Avalue**3) / Tvalue )

                case default
                    print *, 'CODE ERROR: unexpected case default'
                    call util_crashpoint(5098723)
                end select
            end if

        case (2)
            !% --- simpler form 2 equation:
            outvalue = culvert_QIC_unsubmerged_form2 (eIn, Dinlet)

        case default     
            print *, 'CODE ERROR: unexpected default case'  
            call util_crashpoint(598723) 
        end select

        !% --- provide negative flowrate for reversed flow
        if (isReversedFlow) outvalue = -outvalue

    end function culvert_QIC_unsubmerged
!%
!%==========================================================================     
!%==========================================================================  
!% 
    real(8) function culvert_QIC_transition &
         (eIn, Dinlet, QICu, H1s, H1u,&
          isTransition, isReversedFlow) result (outvalue)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the transition regime for culvert, in between the
        !% sumberged and unsubmerged equations
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIn
            real(8), intent(in) :: QICu, Dinlet, H1s, H1u
            logical, intent(in) :: isTransition, isReversedFlow
            real(8), pointer    :: Zbtm
            real(8)             :: QICs
        !%------------------------------------------------------------------
        !% Preliminaries
            if (.not. isTransition) return
        !%------------------------------------------------------------------
            Zbtm => elemR(eIn,er_Zbottom)
        !%------------------------------------------------------------------
        
        if (H1s > H1u) then
            !% --- get the submerged flowrate    
            QICs =  culvert_QIC_submerged_eq (eIn, Dinlet, isReversedFlow)
            outvalue = QICu + (QICs - QICu) * (Dinlet + Zbtm - H1u) &
                                             / (H1s - H1u)
        else
            !% --- inconsistent result in H1s and H1u, so retain QICu value
            outvalue = QICu
        end if

end function culvert_QIC_transition    
!%
!%==========================================================================     
!%==========================================================================  
!%  
    real(8) function culvert_omega &
            (isTabular, geoType, eIn, Psi, Atable, Ttable) result (Omega)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes local Omega = Gamma / Phi, i.e. the normalized area
        !% divided by the normalized topwidth
        !% HACK -- this could be more efficient if the llgeo_... functions
        !% for the non-tabular are provided as normalized functions for
        !% area, topwidth using inputs of normalized depth. Below we must
        !% normalize the area and topwidth before we can use them to get
        !% Omega
        !%------------------------------------------------------------------
        !% Declarations
            logical, intent(in) :: isTabular
            integer, intent(in) :: geoType, eIn
            real(8), intent(in) :: Psi, Atable(:), Ttable(:)
            real(8), pointer    :: sedimentDepth, FullDepth, MaxTopwidth
            real(8), pointer    :: totalPipeArea
            real(8)             :: Atemp, Ttemp
        !%------------------------------------------------------------------
        !% Aliases
            if (.not. isTabular) then
                sedimentDepth => elemR(eIn,er_SedimentDepth)
                FullDepth     => elemR(eIn,er_FullDepth)
                MaxTopwidth   => elemR(eIn,er_BreadthMax)
                totalPipeArea => elemSGR(eIn,esgr_Filled_Circular_TotalPipeArea)
            end if
        !%------------------------------------------------------------------
        if (isTabular) then 
            !% --- simple table lookups
            Omega = xsect_table_lookup_singular(Psi,Atable) &
                  / xsect_table_lookup_singular(Psi,Ttable)
        else 
            select case (geoType)
            case (filled_circular)
                !% --- table lookups for circular pipe are with sediment removed
                !%     but the area must be renormalized to flow area (fulldepth)
                Omega  = ( (                                                          &
                            (totalPipeArea * xsect_table_lookup_singular(Psi,Atable)) &
                            - sedimentDepth                                           &
                           ) / FullDepth                                              &
                         ) / (xsect_table_lookup_singular(Psi,Ttable))

            case (mod_basket)
                Atemp = (llgeo_mod_basket_area_from_depth_singular (eIn, Psi*FullDepth, setting%ZeroValue%Area)) &
                        / FullDepth
                Ttemp = (llgeo_mod_basket_topwidth_from_depth_singular(eIn, Psi*FullDepth, setting%ZeroValue%Topwidth)) &
                        / MaxTopwidth
                Omega = Atemp / Ttemp

            case (rectangular_closed)
                !% --- rectangular case devolves to omega = psi
                Omega = Psi

            case (rect_round)
                Atemp = (llgeo_rect_round_area_from_depth_singular (eIn, Psi*FullDepth, setting%ZeroValue%Area)) &
                        / FullDepth
                Ttemp = (llgeo_rect_round_topwidth_from_depth_singular(eIn, Psi*FullDepth, setting%ZeroValue%Topwidth)) &
                        / MaxTopwidth
                Omega = Atemp / Ttemp

            case (rect_triang)
                Atemp = (llgeo_rectangular_triangular_area_from_depth_singular (eIn,Psi*FullDepth, setting%ZeroValue%Area)) &
                        / FullDepth
                Ttemp = (llgeo_rectangular_triangular_topwidth_from_depth_singular(eIn,Psi*FullDepth, setting%ZeroValue%Topwidth)) &
                        / MaxTopwidth
                Omega = Atemp / Ttemp

            case default
                print *, 'CODE ERROR: Unexpected case default'
                call util_crashpoint(55098723)
            end select
        end if

    end function culvert_omega
!%
!%==========================================================================  
!%==========================================================================  
!%   
    pure real(8) function culvert_residual_form1 &
            (Psi, Omega, halfM, Dhat, Khat, Bhat) result (outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% residual equation of form 1 culvert from SWWM5+ documentation
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(in) :: Dhat, Khat, Omega, halfM, Bhat, Psi
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        
        outvalue = Dhat - Khat * (Omega**halfM) - Bhat * Omega - Psi

    end function culvert_residual_form1
!%
!%==========================================================================  
!%==========================================================================  
!%
    pure real(8) function culvert_Bhat (eIn) result (outvalue)
        !%------------------------------------------------------------------
        !% Descriptions:
        !% computes the \hat{B} constant defined in the SWMM5+ documentation
        !% for equation form 1 culverts
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIn
        !%------------------------------------------------------------------

        outvalue = onehalfR * elemR(eIn,er_FullArea)                    &
                / ( elemR(eIn,er_BreadthMax) * elemR(eIn,er_FullDepth) )

    end function culvert_Bhat
!%
!%==========================================================================  
!%==========================================================================  
!%
    pure real(8) function culvert_Khat (eIn, Bhat) result (outvalue)
        !%------------------------------------------------------------------
        !% Descriptions:
        !% computes the \hat{K} constant defined in the SWMM5+ documentation
        !% for equation form 1 culverts
        !%
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIn
            real(8), intent(in) :: Bhat
        !%------------------------------------------------------------------

        outvalue = elemSR(eIn,esr_Conduit_Culvert_K)                          &
                    * ( (kUnit * sqrt(                                        &
                                      twoR * setting%Constant%gravity * Bhat  &
                                     )**elemSR(eIn,esr_Conduit_Culvert_M)     &
                        ) )
    
    end function culvert_Khat
!%
!%========================================================================== 
!%==========================================================================  
!%
    pure real(8) function culvert_Dhat (eIn, Dinlet, isreversed) result (outvalue)
        !%------------------------------------------------------------------
        !% Descriptions:
        !% computes the \hat{D} constant defined in the SWMM5+ documentation
        !% for equation form 1 culverts
        !% REQUIRES upstream head minus z bottom for culvert inlet to
        !% be stored in elemR(:,er_Temp01)
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIn
            real(8), intent(in) :: Dinlet
            logical, intent(in) :: isreversed
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        if (.not. isreversed) then 
            outvalue = (Dinlet / elemR(eIn,er_FullDepth) )  &
                     - elemSR(eIn,esr_Conduit_Culvert_SCF)                &
                     * elemR(eIn,er_BottomSlope)
        else
            outvalue = Dinlet / elemR(eIn,er_FullDepth)  
        end if

    end function culvert_Dhat
    !%
!%==========================================================================      
!%==========================================================================  
!%
    real(8) function culvert_QIC_unsubmerged_form2 (eIn, Dinlet) result (outvalue)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the inlet controlled flowrate for unsubmerged form 2
        !% equations
        !% eIn is the index of the element with true inflow to the culvert
        !% NOTE this does NOT provide the negative for reversed flow.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIn  
            real(8), intent(in) :: Dinlet
            real(8), pointer    :: Afull, Dfull, Kcoef, Mcoef
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        !% Aliases
            Afull  =>  elemR(eIn,er_FullArea)
            Dfull  =>  elemR(eIn,er_FullDepth)
            Kcoef  => elemSR(eIn,esr_Conduit_Culvert_K)
            Mcoef  => elemSR(eIn,esr_Conduit_Culvert_M)
        !%------------------------------------------------------------------
        
        if (Dinlet > zeroR) then
            outvalue = ( Afull * sqrt(Dfull) / kUnit )                    &
                * (                                                       &
                    ( Dinlet / (Dfull * Kcoef)) ** (oneR / Mcoef)    &
                ) 
        else 
            !% --- degenerate case, QIC is a large number so that
            !%     it is not a constraint
            outvalue = huge(nullvalueR)
        end if

    end function culvert_QIC_unsubmerged_form2 
!%
!%==========================================================================  
!% END OF MODULE
!%==========================================================================
end module culvert_elements