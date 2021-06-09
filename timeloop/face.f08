module face

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Provides computation of face values for timeloop of hydraulics
    !%
    !% METHOD:
    !% 
    !%

    private

    public :: face_interpolation_byMask
    public :: face_interpolation_byPack
    public :: face_interp_across_images
    public :: face_interp_interior 

    contains
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    !%
    subroutine face_interpolation_byMask (faceMaskCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Interpolates all faces using a mask -- assumes single processor
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: faceMaskCol !% Column in face array containing mask for all valid faces
        !%-----------------------------------------------------------------------------
        !%  
    end subroutine face_interpolation_byMask
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine face_interpolation_byPack (thisCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Interpolates all faces using a mask -- assumes single processor
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisCol !% Column in faceP array for needed pack
        integer, intent(in) :: Npack   !% expected number of packed rows in faceP.
        !%-----------------------------------------------------------------------------
        !%  
    end subroutine face_interpolation_byPack
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%   
    subroutine face_interp_across_images ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%  
    end subroutine face_interp_across_images   
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%
    subroutine face_interp_interior ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%  
    end subroutine face_interp_interior
    !%    
    !%==========================================================================
    !% PRIVATE
    !%==========================================================================   
    !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%  
    !%
    !%==========================================================================
    !% END OF MODULE
    !%+=========================================================================
end module face