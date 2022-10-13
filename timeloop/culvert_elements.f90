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
            !% These are (in order) FORM, K, M, C, Y
            !%------------------------------------------------------------------
            !%------------------------------------------------------------------
            !%------------------------------------------------------------------

            !% Circular concrete
            culvertValue(1,:) = (/ 1.0d0, 0.0098d0, 2.00d0, 0.0398d0, 0.67d0 /) !% Square edge w/headwall
            culvertValue(2,:) = (/1.0d0, 0.0018d0, 2.00d0, 0.0292d0, 0.74d0/)  !% Groove end w/headwall
            culvertValue(3,:) = (/1.0d0, 0.0045d0, 2.00d0, 0.0317d0, 0.69d0/)  !% Groove end projecting
        
            !% Circular Corrugated Metal Pipe
            culvertValue(4,:) = (/1.0d0, 0.0078d0, 2.00d0, 0.0379d0, 0.69d0/)  !% Headwall
            culvertValue(5,:) = (/1.0d0, 0.0210d0, 1.33d0, 0.0463d0, 0.75d0/) !% Mitered to slope
            culvertValue(6,:) = (/1.0d0, 0.0340d0, 1.50d0, 0.0553d0, 0.54d0/)  !% Projecting
        
            !% Circular Pipe, Beveled Ring Entrance
            culvertValue(7,:) = (/1.0d0, 0.0018d0, 2.50d0, 0.0300d0, 0.74d0/)  !% Beveled ringd0, 45 deg bevels
            culvertValue(8,:) = (/1.0d0, 0.0018d0, 2.50d0, 0.0243d0, 0.83d0/)  !% Beveled ringd0, 33.7 deg bevels
        
            !% Rectangular Box with Flared Wingwalls
            culvertValue(9,:) = (/1.0d0, 0.026d0, 1.0d0,   0.0347d0, 0.81d0/)  !% 30-75 deg. wingwall flares
            culvertValue(10,:) = (/1.0d0, 0.061d0, 0.75d0,  0.0400d0, 0.80d0/)  !% 90 or 15 deg. wingwall flares
            culvertValue(11,:) = (/1.0d0, 0.061d0, 0.75d0,  0.0423d0, 0.82d0/)  !% 0 deg. wingwall flares (striaght sides)
        
            !% Rectanglar Box with Flared Wingwalls & Top Edge Bevel
            culvertValue(12,:) = (/2.0d0, 0.510d0, 0.667d0, 0.0309d0, 0.80d0/)  !% 45 deg. flare; 0.43D top edge bevel
            culvertValue(13,:) = (/2.0d0, 0.486d0, 0.667d0, 0.0249d0, 0.83d0/)  !% 18-33.7 deg flare; 0.083D top edge bevel
        
            !% Rectangular Box; 90-deg Headwall; Chamfered or Beveled Inlet Edges
            culvertValue(14,:) = (/2.0d0, 0.515d0, 0.667d0, 0.0375d0, 0.79d0/)  !% chamfered 3/4-in
            culvertValue(15,:) = (/2.0d0, 0.495d0, 0.667d0, 0.0314d0, 0.82d0/)  !% beveled 1/2-in/ft at 45 deg (1:1)
            culvertValue(16,:) = (/2.0d0, 0.486d0, 0.667d0, 0.0252d0, 0.865d0/) !% beveled 1-in/ft at 33.7 deg (1:1.5)
        
            !% Rectangular Box; Skewed Headwall; Chamfered or Beveled Inlet Edges
            culvertValue(17,:) = (/2.0d0, 0.545d0, 0.667d0, 0.04505d0,0.73d0/)  !% 3/4" chamfered edged0, 45 deg skewed headwall
            culvertValue(18,:) = (/2.0d0, 0.533d0, 0.667d0, 0.0425d0, 0.705d0/) !% 3/4" chamfered edged0, 30 deg skewed headwall
            culvertValue(19,:) = (/2.0d0, 0.522d0, 0.667d0, 0.0402d0, 0.68d0/)  !% 3/4" chamfered edged0, 15 deg skewed headwall
            culvertValue(20,:) = (/2.0d0, 0.498d0, 0.667d0, 0.0327d0, 0.75d0/)  !% 45 deg beveled edged0, 10-45 deg skewed headwall
        
            !% Rectangular box, Non-offset Flared Wingwalls; 3/4" Chamfer at Top of Inlet
            culvertValue(21,:) = (/2.0d0, 0.497d0, 0.667d0, 0.0339d0, 0.803d0/) !% 45 deg (1:1) wingwall flare
            culvertValue(22,:) = (/2.0d0, 0.493d0, 0.667d0, 0.0361d0, 0.806d0/) !% 18.4 deg (3:1) wingwall flare
            culvertValue(23,:) = (/2.0d0, 0.495d0, 0.667d0, 0.0386d0, 0.71d0/)  !% 18.4 deg (3:1) wingwall flared0, 30 deg inlet skew
        
            !% Rectangular box, Offset Flared Wingwallsd0, Beveled Edge at Inlet Top
            culvertValue(24,:) = (/2.0d0, 0.497d0, 0.667d0, 0.0302d0, 0.835d0/)  !% 45 deg (1:1) flared0, 0.042D top edge bevel
            culvertValue(25,:) = (/2.0d0, 0.495d0, 0.667d0, 0.0252d0, 0.881d0/)  !% 33.7 deg (1.5:1) flared0, 0.083D top edge bevel
            culvertValue(26,:) = (/2.0d0, 0.493d0, 0.667d0, 0.0227d0, 0.887d0/)  !% 18.4 deg (3:1) flared0, 0.083D top edge bevel
        
            !%  Corrugated Metal Box
            culvertValue(27,:) = (/1.0d0, 0.0083d0, 2.00d0, 0.0379d0, 0.69d0/)  !% 90 deg headwall
            culvertValue(28,:) = (/1.0d0, 0.0145d0, 1.75d0, 0.0419d0, 0.64d0/)  !% Thick wall projecting
            culvertValue(29,:) = (/1.0d0, 0.0340d0, 1.50d0, 0.0496d0, 0.57d0/)  !% Thin wall projecting
        
            !%  Horizontal Ellipse Concrete
            culvertValue(30,:) = (/1.0d0, 0.0100d0, 2.00d0, 0.0398d0, 0.67d0/)  !% Square edge w/headwall
            culvertValue(31,:) = (/1.0d0, 0.0018d0, 2.50d0, 0.0292d0, 0.74d0/)  !% Grooved end w/headwall
            culvertValue(32,:) = (/1.0d0, 0.0045d0, 2.00d0, 0.0317d0, 0.69d0/)  !% Grooved end projecting
        
            !%  Vertical Ellipse Concrete
            culvertValue(33,:) = (/1.0d0, 0.0100d0, 2.00d0, 0.0398d0, 0.67d0/)  !% Square edge w/headwall
            culvertValue(34,:) = (/1.0d0, 0.0018d0, 2.50d0, 0.0292d0, 0.74d0/)  !% Grooved end w/headwall
            culvertValue(35,:) = (/1.0d0, 0.0095d0, 2.00d0, 0.0317d0, 0.69d0/)  !% Grooved end projecting
        
            !%  Pipe Arch, 18" Corner Radiusd0, Corrugated Metal
            culvertValue(36,:) = (/1.0d0, 0.0083d0, 2.00d0, 0.0379d0, 0.69d0/)  !% 90 deg headwall
            culvertValue(37,:) = (/1.0d0, 0.0300d0, 1.00d0, 0.0463d0, 0.75d0/)  !% Mitered to slope
            culvertValue(38,:) = (/1.0d0, 0.0340d0, 1.50d0, 0.0496d0, 0.57d0/)  !% Projecting
        
            !%  Pipe Arch, 18" Corner Radiusd0, Corrugated Metal
            culvertValue(39,:) = (/1.0d0, 0.0300d0, 1.50d0, 0.0496d0, 0.57d0/)  !% Projecting
            culvertValue(40,:) = (/1.0d0, 0.0088d0, 2.00d0, 0.0368d0, 0.68d0/)  !% No bevels
            culvertValue(41,:) = (/1.0d0, 0.0030d0, 2.00d0, 0.0269d0, 0.77d0/)  !% 33.7 deg bevels
        
            !%  Pipe Arch, 31" Corner Radiusd0, Corrugated Metal
            culvertValue(42,:) = (/1.0d0, 0.0300d0, 1.50d0, 0.0496d0, 0.57d0/)  !% Projecting
            culvertValue(43,:) = (/1.0d0, 0.0088d0, 2.00d0, 0.0368d0, 0.68d0/)  !% No bevels
            culvertValue(44,:) = (/1.0d0, 0.0030d0, 2.00d0, 0.0269d0, 0.77d0/)  !% 33.7 deg. bevels
        
            !%  Arch, Corrugated Metal
            culvertValue(45,:) = (/1.0d0, 0.0083d0, 2.00d0, 0.0379d0, 0.69d0/)  !% 90 deg headwall
            culvertValue(46,:) = (/1.0d0, 0.0300d0, 1.00d0, 0.0473d0, 0.75d0/)  !% Mitered to slope                     !% (5.1.013)
            culvertValue(47,:) = (/1.0d0, 0.0340d0, 1.50d0, 0.0496d0, 0.57d0/)  !% Thin wall projecting
        
            !%  Circular Culvert
            culvertValue(48,:) = (/2.0d0, 0.534d0, 0.555d0, 0.0196d0, 0.90d0/)  !% Smooth tapered inlet throat
            culvertValue(49,:) = (/2.0d0, 0.519d0, 0.640d0, 0.0210d0, 0.90d0/)  !% Rough tapered inlet throat
        
            !%  Elliptical Inlet Face
            culvertValue(50,:) = (/2.0d0, 0.536d0, 0.622d0, 0.0368d0, 0.83d0/)  !% Tapered inletd0, beveled edges
            culvertValue(51,:) = (/2.0d0, 0.5035d0,0.719d0, 0.0478d0, 0.80d0/)  !% Tapered inletd0, square edges
            culvertValue(52,:) = (/2.0d0, 0.547d0, 0.800d0, 0.0598d0, 0.75d0/)  !% Tapered inletd0, thin edge projecting
        
            !%  Rectangular
            culvertValue(53,:) = (/2.0d0, 0.475d0, 0.667d0, 0.0179d0, 0.97d0/)  !% Tapered inlet throat
        
            !%  Rectangular Concrete
            culvertValue(54,:) = (/2.0d0, 0.560d0, 0.667d0, 0.0446d0, 0.85d0/)  !% Side taperedd0, less favorable edges
            culvertValue(55,:) = (/2.0d0, 0.560d0, 0.667d0, 0.0378d0, 0.87d0/)  !% Side taperedd0, more favorable edges
        
            !%  Rectangular Concrete
            culvertValue(56,:) = (/2.0d0, 0.500d0, 0.667d0, 0.0446d0, 0.65d0/) !% Slope taperedd0, less favorable edges
            culvertValue(57,:) = (/2.0d0, 0.500d0, 0.667d0, 0.0378d0, 0.71d0/)  !% Slope taperedd0, more favorable edges
        

        end subroutine culvert_parameter_values
    !%
    !%==========================================================================
    !% END OF MODULE
    !%=========================================================================
end module culvert_elements