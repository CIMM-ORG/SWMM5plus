module face

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use jump

    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Provides computation of face values for timeloop of hydraulics
    !%
    !% METHOD:
    !% 
    !%

    private

    public :: face_interpolation

    contains
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    !%
    subroutine face_interpolation (facecol, isMask)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Interpolates faces
        !%-----------------------------------------------------------------------------
        integer, intent(in)  :: faceCol
        logical, intent(in)  :: isMask
        integer, pointer :: Npack
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'face_interpolation'
        if (setting%Debug%File%face) print *, '*** enter ', subroutine_name
        !%-----------------------------------------------------------------------------
        !% 
        if (num_images() == 1) then
            if (isMask) then
                call face_interpolation_byMask(faceCol)
            else
                Npack => npack_faceP(faceCol)
                if (Npack > 0) then
                    call face_interpolation_byPack (faceCol, Npack)
                endif
            endif
        else
            call face_interp_across_images()
            call face_interp_interior()
        end if 
        
        if (setting%Debug%File%face)  print *, '*** leave ', subroutine_name
    end subroutine face_interpolation
    !%    
    !%==========================================================================
    !% PRIVATE
    !%==========================================================================  
    !%   
    subroutine face_interpolation_byMask (faceMaskCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Interpolates all faces using a mask -- assumes single processor
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: faceMaskCol !% Column in face array containing mask for all valid faces
        integer :: fGeoSetU(3), fGeoSetD(3), eGeoSet(3)
        integer :: fHeadSetU(1), fHeadSetD(1), eHeadSet(1)
        integer :: fFlowSet(1), eFlowSet(1)
        !%-----------------------------------------------------------------------------
        !% Face values are needed for
        !% Area_u, Area_d, Head_u, Head_d, Flowrate, 
        
        !% not sure if we need
        !% Topwidth_u, Topwidth_d, HydDepth_u, HydDepth_d
        !% Velocity_u, Velocity_d
        
        !% General approach
        !% interpolate to ..._u
        !% identify hydraulic jumps
        !% set .._u and ..d based on jumps
        
        !% set the matching sets
        !% THESE SHOULD BE DONE IN A GLOBAL -- MAYBE SETTINGS
        !% Note these can be expanded for other terms to be interpolated.
        fGeoSetU = [fr_Area_u, fr_Topwidth_u, fr_HydDepth_u]
        fGeoSetD = [fr_Area_d, fr_Topwidth_d, fr_HydDepth_d]
        eGeoSet  = [er_Area,   er_Topwidth,   er_HydDepth]
        
        fHeadSetU = [fr_Head_u]
        fHeadSetD = [fr_Head_d]
        eHeadSet = [er_Head]
        
        fFlowSet = [fr_Flowrate]
        eFlowSet = [er_Flowrate]
        
        !% two-sided interpolation to using the upstream face set
        call face_interp_set_byMask &
            (fGeoSetU, eGeoSet, er_InterpWeight_dG, er_InterpWeight_uG, faceMaskCol)
        call face_interp_set_byMask &
            (fHeadSetU, eHeadSet, er_InterpWeight_dH, er_InterpWeight_uH, faceMaskCol)
        call face_interp_set_byMask &
            (fFlowSet, eFlowSet, er_InterpWeight_dQ, er_InterpWeight_uQ, faceMaskCol)
        
        !% copy upstream to downstream storage at a face
        !% (only for Head and Geometry types as flow has a single value)
        !% note that these might be reset by hydraulic jump
        call face_copy_upstream_to_downstream_byMask (fGeoSetD,  fGeoSetU,  faceMaskCol)
        call face_copy_upstream_to_downstream_byMask (fHeadSetD, fHeadSetU, faceMaskCol)
        
        !% reset all the hydraulic jump faces
        call jump_compute
        
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
    !%==========================================================================    
    !%  
    subroutine face_interp_set_byMask &
        (fset, eset, eWdn, eWup, faceMaskCol) 
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Interpolates to a face for a set of variables using a mask
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: fset(:), eset(:), eWdn, eWup, faceMaskCol
        integer, pointer :: eup(:), edn(:)
        integer :: ii
        !%-----------------------------------------------------------------------------
        eup => faceI(:,fi_Melem_uL)
        edn => faceI(:,fi_Melem_dL) 
        !%-----------------------------------------------------------------------------

        !% cycle through each element in the set.
        do ii=1,size(fset)
            where (faceM(:,faceMaskCol))
                faceR(:,fset(ii)) = &
                    (+elemR(eup(:),eset(ii)) * elemR(edn(:),eWup) &
                     +elemR(edn(:),eset(ii)) * elemR(eup(:),eWdn) &
                    ) / &
                    ( elemR(edn(:),eWup) + elemR(eup(:),eWdn))
            endwhere    
        enddo
                
    end subroutine face_interp_set_byMask
    !%
    !%========================================================================== 
    !%==========================================================================    
    !%  
    subroutine face_copy_upstream_to_downstream_byMask &
        (downstreamSet, upstreamSet, faceMaskCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Copies the interpolated value on the upstrea side to the downstream side
        !% These values are later adjusted for hydraulic jumps
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: faceMaskCol, downstreamSet(:), upstreamSet(:)
        integer :: ii
        !%-----------------------------------------------------------------------------

        do ii=1,size(downstreamset)
            where (faceM(:,faceMaskCol))
                faceR(:,downstreamSet(ii)) = faceR(:,upstreamSet(ii))
            endwhere
        enddo
    
    end subroutine face_copy_upstream_to_downstream_byMask
    !%
    !%========================================================================== 
    !%==========================================================================    
    !%  
    !%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step 
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%  
    !%
    !%========================================================================== 
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
