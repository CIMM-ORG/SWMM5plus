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
            culvertValue(1,:) = (/ 1.0, 0.0098, 2.00, 0.0398, 0.67 /) !% Square edge w/headwall
            culvertValue(2,:) = (/1.0, 0.0018, 2.00, 0.0292, 0.74/)  !% Groove end w/headwall
            culvertValue(3,:) = (/1.0, 0.0045, 2.00, 0.0317, 0.69/)  !% Groove end projecting
        
            !% Circular Corrugated Metal Pipe
            culvertValue(4,:) = (/1.0, 0.0078, 2.00, 0.0379, 0.69/)  !% Headwall
            culvertValue(5,:) = (/1.0, 0.0210, 1.33, 0.0463, 0.75/) !% Mitered to slope
            culvertValue(6,:) = (/1.0, 0.0340, 1.50, 0.0553, 0.54/)  !% Projecting
        
            !% Circular Pipe, Beveled Ring Entrance
            culvertValue(7,:) = (/1.0, 0.0018, 2.50, 0.0300, 0.74/)  !% Beveled ring, 45 deg bevels
            culvertValue(8,:) = (/1.0, 0.0018, 2.50, 0.0243, 0.83/)  !% Beveled ring, 33.7 deg bevels
        
            !% Rectangular Box with Flared Wingwalls
            culvertValue(9,:) = (/1.0, 0.026, 1.0,   0.0347, 0.81/)  !% 30-75 deg. wingwall flares
            culvertValue(10,:) = (/1.0, 0.061, 0.75,  0.0400, 0.80/)  !% 90 or 15 deg. wingwall flares
            culvertValue(11,:) = (/1.0, 0.061, 0.75,  0.0423, 0.82/)  !% 0 deg. wingwall flares (striaght sides)
        
            !% Rectanglar Box with Flared Wingwalls & Top Edge Bevel
            culvertValue(12,:) = (/2.0, 0.510, 0.667, 0.0309, 0.80/)  !% 45 deg. flare; 0.43D top edge bevel
            culvertValue(13,:) = (/2.0, 0.486, 0.667, 0.0249, 0.83/)  !% 18-33.7 deg flare; 0.083D top edge bevel
        
            !% Rectangular Box; 90-deg Headwall; Chamfered or Beveled Inlet Edges
            culvertValue(14,:) = (/2.0, 0.515, 0.667, 0.0375, 0.79/)  !% chamfered 3/4-in
            culvertValue(15,:) = (/2.0, 0.495, 0.667, 0.0314, 0.82/)  !% beveled 1/2-in/ft at 45 deg (1:1)
            culvertValue(16,:) = (/2.0, 0.486, 0.667, 0.0252, 0.865/) !% beveled 1-in/ft at 33.7 deg (1:1.5)
        
            !% Rectangular Box; Skewed Headwall; Chamfered or Beveled Inlet Edges
            culvertValue(17,:) = (/2.0, 0.545, 0.667, 0.04505,0.73/)  !% 3/4" chamfered edge, 45 deg skewed headwall
            culvertValue(18,:) = (/2.0, 0.533, 0.667, 0.0425, 0.705/) !% 3/4" chamfered edge, 30 deg skewed headwall
            culvertValue(19,:) = (/2.0, 0.522, 0.667, 0.0402, 0.68/)  !% 3/4" chamfered edge, 15 deg skewed headwall
            culvertValue(20,:) = (/2.0, 0.498, 0.667, 0.0327, 0.75/)  !% 45 deg beveled edge, 10-45 deg skewed headwall
        
            !% Rectangular box, Non-offset Flared Wingwalls; 3/4" Chamfer at Top of Inlet
            culvertValue(21,:) = (/2.0, 0.497, 0.667, 0.0339, 0.803/) !% 45 deg (1:1) wingwall flare
            culvertValue(22,:) = (/2.0, 0.493, 0.667, 0.0361, 0.806/) !% 18.4 deg (3:1) wingwall flare
            culvertValue(23,:) = (/2.0, 0.495, 0.667, 0.0386, 0.71/)  !% 18.4 deg (3:1) wingwall flare, 30 deg inlet skew
        
            !% Rectangular box, Offset Flared Wingwalls, Beveled Edge at Inlet Top
            culvertValue(24,:) = (/2.0, 0.497, 0.667, 0.0302, 0.835/)  !% 45 deg (1:1) flare, 0.042D top edge bevel
            culvertValue(25,:) = (/2.0, 0.495, 0.667, 0.0252, 0.881/)  !% 33.7 deg (1.5:1) flare, 0.083D top edge bevel
            culvertValue(26,:) = (/2.0, 0.493, 0.667, 0.0227, 0.887/)  !% 18.4 deg (3:1) flare, 0.083D top edge bevel
        
            !%  Corrugated Metal Box
            culvertValue(27,:) = (/1.0, 0.0083, 2.00, 0.0379, 0.69/)  !% 90 deg headwall
            culvertValue(28,:) = (/1.0, 0.0145, 1.75, 0.0419, 0.64/)  !% Thick wall projecting
            culvertValue(29,:) = (/1.0, 0.0340, 1.50, 0.0496, 0.57/)  !% Thin wall projecting
        
            !%  Horizontal Ellipse Concrete
            culvertValue(30,:) = (/1.0, 0.0100, 2.00, 0.0398, 0.67/)  !% Square edge w/headwall
            culvertValue(31,:) = (/1.0, 0.0018, 2.50, 0.0292, 0.74/)  !% Grooved end w/headwall
            culvertValue(32,:) = (/1.0, 0.0045, 2.00, 0.0317, 0.69/)  !% Grooved end projecting
        
            !%  Vertical Ellipse Concrete
            culvertValue(33,:) = (/1.0, 0.0100, 2.00, 0.0398, 0.67/)  !% Square edge w/headwall
            culvertValue(34,:) = (/1.0, 0.0018, 2.50, 0.0292, 0.74/)  !% Grooved end w/headwall
            culvertValue(35,:) = (/1.0, 0.0095, 2.00, 0.0317, 0.69/)  !% Grooved end projecting
        
            !%  Pipe Arch, 18" Corner Radius, Corrugated Metal
            culvertValue(36,:) = (/1.0, 0.0083, 2.00, 0.0379, 0.69/)  !% 90 deg headwall
            culvertValue(37,:) = (/1.0, 0.0300, 1.00, 0.0463, 0.75/)  !% Mitered to slope
            culvertValue(38,:) = (/1.0, 0.0340, 1.50, 0.0496, 0.57/)  !% Projecting
        
            !%  Pipe Arch, 18" Corner Radius, Corrugated Metal
            culvertValue(39,:) = (/1.0, 0.0300, 1.50, 0.0496, 0.57/)  !% Projecting
            culvertValue(40,:) = (/1.0, 0.0088, 2.00, 0.0368, 0.68/)  !% No bevels
            culvertValue(41,:) = (/1.0, 0.0030, 2.00, 0.0269, 0.77/)  !% 33.7 deg bevels
        
            !%  Pipe Arch, 31" Corner Radius, Corrugated Metal
            culvertValue(42,:) = (/1.0, 0.0300, 1.50, 0.0496, 0.57/)  !% Projecting
            culvertValue(43,:) = (/1.0, 0.0088, 2.00, 0.0368, 0.68/)  !% No bevels
            culvertValue(44,:) = (/1.0, 0.0030, 2.00, 0.0269, 0.77/)  !% 33.7 deg. bevels
        
            !%  Arch, Corrugated Metal
            culvertValue(45,:) = (/1.0, 0.0083, 2.00, 0.0379, 0.69/)  !% 90 deg headwall
            culvertValue(46,:) = (/1.0, 0.0300, 1.00, 0.0473, 0.75/)  !% Mitered to slope                     !% (5.1.013)
            culvertValue(47,:) = (/1.0, 0.0340, 1.50, 0.0496, 0.57/)  !% Thin wall projecting
        
            !%  Circular Culvert
            culvertValue(48,:) = (/2.0, 0.534, 0.555, 0.0196, 0.90/)  !% Smooth tapered inlet throat
            culvertValue(49,:) = (/2.0, 0.519, 0.640, 0.0210, 0.90/)  !% Rough tapered inlet throat
        
            !%  Elliptical Inlet Face
            culvertValue(50,:) = (/2.0, 0.536, 0.622, 0.0368, 0.83/)  !% Tapered inlet, beveled edges
            culvertValue(51,:) = (/2.0, 0.5035,0.719, 0.0478, 0.80/)  !% Tapered inlet, square edges
            culvertValue(52,:) = (/2.0, 0.547, 0.800, 0.0598, 0.75/)  !% Tapered inlet, thin edge projecting
        
            !%  Rectangular
            culvertValue(53,:) = (/2.0, 0.475, 0.667, 0.0179, 0.97/)  !% Tapered inlet throat
        
            !%  Rectangular Concrete
            culvertValue(54,:) = (/2.0, 0.560, 0.667, 0.0446, 0.85/)  !% Side tapered, less favorable edges
            culvertValue(55,:) = (/2.0, 0.560, 0.667, 0.0378, 0.87/)  !% Side tapered, more favorable edges
        
            !%  Rectangular Concrete
            culvertValue(56,:) = (/2.0, 0.500, 0.667, 0.0446, 0.65/) !% Slope tapered, less favorable edges
            culvertValue(57,:) = (/2.0, 0.500, 0.667, 0.0378, 0.71/)  !% Slope tapered, more favorable edges
        

        end subroutine culvert_parameter_values
    !%
    !%==========================================================================
    !% END OF MODULE
    !%=========================================================================
end module culvert_elements