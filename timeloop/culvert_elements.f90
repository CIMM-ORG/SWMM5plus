module culvert_elements

    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use utility, only: util_CLprint
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
            integer, pointer :: thisCol, npack, eIn, eOut, EqForm
            integer, pointer :: fUp(:), fDn(:), thisE(:)
            real(8), pointer :: fHeadD(:), fHeadU(:), Dinlet(:), Flowrate(:)
            real(8), pointer :: Zbtm(:)
            integer          :: ii

            real(8) :: QIC

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

            Dinlet   => elemR(:,er_Temp01) !% storage for effective inlet depth
        !%------------------------------------------------------------------

        !% --- cycle through the culverts
        do ii=1,npack
            !% --- inlet and outlet
            eIn    => thisE(ii)
            eOut   => elemSI(eIn,esi_Conduit_Culvert_OutletID)
            EqForm => elemSI(eIn,esi_Conduit_Culvert_EquationForm)

            !% --- check flow direction
            if ((Flowrate(eIn) > zeroR) .and. (Flowrate(eOut) > zeroR)) then 
                !% --- downstream flow

                !% --- effective inlet flow depth
                Dinlet(eIn) = fHeadD(fUp(eIn)) - Zbtm(eIn)

                !% --- no action if negative inlet depth
                if (Dinlet(eIn) .le. zeroR) cycle
                
                !% --- submergence
                if (issubmerged(eIn,fUp(eIn),.false.)) then
                    QIC =  culvert_submerged (eIn, fUp(eIn), .false.)
                else

                    !% --- compute QIC
                    select case (EqForm)
                    case (1)
                    case (2)
                        QIC = culvert_unsubmerged_form2 (eIn,fUp(eIn))
                    case default
                        print *, 'CODE ERROR: unexpected case default'
                        call util_crashpoint(692873)
                    end select
                end if


            elseif ((Flowrate(eIn) < zeroR) .and. (Flowrate(eOut) < zeroR)) then 
                !% --- upstream (reversed) flow

                !% --- effective inlet flow depth
                Dinlet(eOut) = fHeadU(fDn(eOut)) - Zbtm(eOut)

                !% --- no action if negative inlet depth
                if (Dinlet(eOut) .le. zeroR) cycle
                
                !% --- submergence
                if (issubmerged(eOut,fDn(eOut),.true.)) then
                    QIC =  culvert_submerged (eOut, fDn(eOut), .true.)
                else
                    !% --- compute QIC
                    select case (EqForm)
                    case (1)
                    case (2)
                        QIC = - culvert_unsubmerged_form2 (eOut,fDn(eOut))
                    case default
                        print *, 'CODE ERROR: unexpected case default'
                        call util_crashpoint(6928735)
                    end select
                end if
            else 
                !% --- inconsistent flow directions
                !%     skip culvert computations entirely
            end if
            

        end do


    end subroutine culvert_toplevel
!%
!%==========================================================================  
!%==========================================================================  
!%   
    logical function issubmerged (eIn, fIn, isreversed) result (outvalue)
        !%------------------------------------------------------------------
        !% Description
        !% determines whether the culvert meets the submerged flow condition
        !% of inlet head greater than crown of culvert or equation 7-48
        !% of the SWMM5 Reference Manual Vol II Hydraulics.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIn, fIn
            logical, intent(in) :: isreversed
            real(8), pointer :: Dfull, Hup, Zbtm, Zcrown, slope
            real(8), pointer :: Ccoef, Ycoef, SCF
        !%------------------------------------------------------------------
        !% Aliases
            Dfull  =>  elemR(eIn,er_FullDepth)
            if (.not. isreversed) then
                Hup    =>  faceR(fIn,fr_Head_d)
            else
                Hup    =>  faceR(fIn,fr_Head_u)
            end if
            Zbtm   =>  elemR(eIn,er_Zbottom)
            Zcrown =>  elemR(eIn,er_Zcrown)
            slope  =>  elemR(eIn,er_BottomSlope)
            Ccoef  => elemSR(eIn,esr_Conduit_Culvert_C)
            Ycoef  => elemSR(eIn,esr_Conduit_Culvert_Y)
            SCF    => elemSR(eIn,esr_Conduit_Culvert_SCF)
        !%------------------------------------------------------------------


        if (.not. isreversed) then
            if ( ( Hup > Zcrown)                                              &
                .or.                                                          &
                 (Hup > (ZBtm + Dfull * (16.d0 * Ccoef + Ycoef + SCF*slope))) &
                ) then 

                outvalue = .true.

            else 
                outvalue = .false.
            end if
        else 
            !% --- reverse flow neglects slope term
            if ( (Hup > Zcrown)                                                &
                .or.                                                           &
                (Hup > (ZBtm + Dfull * (16.d0 * Ccoef + Ycoef)))               &
                ) then 

                outvalue = .true.

            else 
                outvalue = .false.
            end if
            
        end if

    end function issubmerged
!%
!%==========================================================================  
!%==========================================================================  
!%
    subroutine culvert_unsubmerged_form1 ()


        ! Dhat = culvert_Dhat (eIn, isreversed)

        ! Dhat = Di
        
    end subroutine culvert_unsubmerged_form1 
!%
!%==========================================================================  
!%==========================================================================  
!%
    pure real(8) function culvert_Dhat (eIn, isreversed) result (outvalue)
        !%------------------------------------------------------------------
        !% Descriptions:
        !% computes the \hat{D} constant defined in the SWMM5+ documentation
        !% for equation form 1 culverts
        !% REQUIRES upstream head minus z bottom for culvert inlet to
        !% be stored in elemR(:,er_Temp01)
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIn
            logical, intent(in) :: isreversed
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        if (.not. isreversed) then 
            outvalue = (elemR(eIn,er_Temp01) / elemR(eIn,er_FullDepth) )  &
                     - elemSI(eIn,esr_Conduit_Culvert_SCF)                &
                     * elemR(eIn,er_BottomSlope)
        else
            outvalue = elemR(eIn,er_Temp01) / elemR(eIn,er_FullDepth)  
        end if

    end function culvert_Dhat
    !%
!%==========================================================================  
!%==========================================================================  
!%
    ! pure real(8) function culvert_Bhat (eIn) result (outvalue)
    !     !%------------------------------------------------------------------
    !     !% Descriptions:
    !     !% computes the \hat{B} constant defined in the SWMM5+ documentation
    !     !% for equation form 1 culverts
    !     !%------------------------------------------------------------------

    !     outvalue = onehalfR * elemR(eIn,er_FullArea)
    !             / (elemR(eIn,er_T))
    ! end function culvert_Bhat
!%
!%==========================================================================  
!%==========================================================================  
!%
    pure real(8) function culvert_Khat (eIn) result (outvalue)
        !%------------------------------------------------------------------
        !% Descriptions:
        !% computes the \hat{K} constant defined in the SWMM5+ documentation
        !% for equation form 1 culverts
        !% REQUIRES the B term from SWMM5+ derivation to be stored in
        !% elemR(:,erTemp02)
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIn
        !%------------------------------------------------------------------


        outvalue = elemSR(eIn,esr_Conduit_Culvert_K)                  &
                    * (                                               &
                        (                                             &
                        kUnit * sqrt(twoR * setting%Constant%gravity  &
                                      * elemR(eIn,er_Temp02))            &
                        )**elemSR(eIn,esr_Conduit_Culvert_M)          &
                    )
    
    end function culvert_Khat
!%
!%==========================================================================  
!%==========================================================================  
!%
    real(8) function culvert_unsubmerged_form2 (eIn,fIn) result (outvalue)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the inlet controlled flowrate for unsubmerged form 2
        !% equations
        !% eIn is the index of the inlet to the culvert and fUp is the
        !% face upstream of the inlet. Note that for reversed flow
        !% the eIn stores the conventional "outlet" and the fUp is a
        !% the nominal downstream face that is feeding the culvert.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIn, fIn   
            real(8), pointer    :: Afull, Dfull, Dinlet, Kcoef, Mcoef
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        !% Aliases
            Afull  =>  elemR(eIn,er_FullArea)
            Dfull  =>  elemR(eIn,er_FullDepth)
            Dinlet =>  elemR(eIn,er_Temp01)
            Kcoef  => elemSR(eIn,esr_Conduit_Culvert_K)
            Mcoef  => elemSR(eIn,esr_Conduit_Culvert_M)
        !%------------------------------------------------------------------
        
        if (Dinlet > zeroR) then
            outvalue = ( Afull * sqrt(Dfull) / kUnit )                    &
                * (                                                       &
                    ( Dinlet / (Dfull * Kcoef)) ** (oneR / Mcoef)    &
                ) 
        else 
            !% --- degenerate case, QIC is a large number
            outvalue = huge(nullvalueR)
        end if


    end function culvert_unsubmerged_form2 
!%
!%==========================================================================  
!%==========================================================================  
!%
    real(8) function culvert_submerged (eIn, fUp, isreversed) result (outvalue)
        !%------------------------------------------------------------------
        !% Description
        !% Computes the inlet submerged form of the culvert equations
        !% Output is the inlet-controlled flowrate
        !% eIn is the index of the inlet to the culvert and fUp is the
        !% face upstream of the inlet. Note that for reversed flow
        !% the eIn stores the conventional "outlet" and the fUp is a
        !% the nominal downstream face that is feeding the culvert.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIn, fUp
            logical, intent(in) :: isreversed
            real(8), pointer :: Afull, Dfull, Dinlet, slope
            real(8), pointer :: Kcoef, Ccoef, Ycoef, SCF
        !%------------------------------------------------------------------
            Afull  =>  elemR(eIn,er_FullArea)
            Dfull  =>  elemR(eIn,er_FullDepth)
            Dinlet =>  elemR(eIn,er_Temp01)
            slope  =>  elemR(eIn,er_BottomSlope)
            Kcoef  => elemSR(eIn,esr_Conduit_Culvert_K)
            Ccoef  => elemSR(eIn,esr_Conduit_Culvert_C)
            Ycoef  => elemSR(eIn,esr_Conduit_Culvert_Y)
            SCF    => elemSR(eIn,esr_Conduit_Culvert_SCF)
        !%------------------------------------------------------------------

        if (.not. isreversed) then
            !% --- conventional downstream flow
            outvalue = ( Afull * sqrt(Dfull) / kUnit )                  &
                    * sqrt(                                             &
                            ( ( Dinlet / Dfull) - Ycoef - SCF * slope ) &
                            / Ccoef                                     &
                        ) 
        else
            !% --- reversed flow in culvert (neglect slope term)
            outvalue = ( Afull * sqrt(Dfull) / kUnit )                  &
                    * sqrt(                                             &
                            ( ( Dinlet / Dfull) - Ycoef )               &
                            / Ccoef                                     &
                        ) 
        end if


    end function culvert_submerged
!%
!%==========================================================================   
!% END OF MODULE
!%==========================================================================
end module culvert_elements