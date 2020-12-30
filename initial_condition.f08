!
! module initial_condition
!
! Provides setup of initial conditions for geometry and dynamics
!
!==========================================================================
!
module initial_condition
    !
    ! initial conditions on the elements and faces
    !
    use array_index
    use bc
    use data_keys
    use element_geometry
    use element_dynamics
    use face_values
    use globals
    use junction
    use setting_definition
    use storage
    use utility
    use orifice
    use weir
    use xsect_tables

    implicit none

    private

    public  :: initial_condition_setup

    integer :: debuglevel = 0

contains
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine initial_condition_setup &
        (elem2R, elem2I, elem2YN, elemMR, elemMI, elemMYN, faceR, faceI, faceYN, &
        linkR, linkI, nodeR, nodeI, bcdataDn, bcdataUp, thisTime, ID,           &
        numberPairs, ManningsN, Length, zBottom, xDistance, Breadth,            &
        widthDepthData, cellType)

        character(64) :: subroutine_name = 'initial_condition_setup'

        real,      intent(in out)  :: elem2R(:,:),  elemMR(:,:), faceR(:,:)
        integer,   intent(in out)  :: elem2I(:,:),  elemMI(:,:), faceI(:,:)
        logical,   intent(in out)  :: elem2YN(:,:), elemMYN(:,:), faceYN(:,:)

        real,                intent(in)      :: linkR(:,:), nodeR(:,:)
        integer,   target,   intent(in)      :: linkI(:,:), nodeI(:,:)
        real,                intent(in)      :: thisTime

        type(bcType),        intent(in out)      :: bcdataDn(:), bcdataUp(:)

        integer, intent(inout)    :: ID(:)
        integer, intent(inout)    :: numberPairs(:)
        real,    intent(inout)    :: ManningsN(:)
        real,    intent(inout)    :: Length(:)
        real,    intent(inout)    :: zBottom(:)
        real,    intent(inout)    :: xDistance(:)
        real,    intent(inout)    :: Breadth(:)
        real,    intent(inout)    :: widthDepthData(:,:,:)

        type(string), intent(in out)   :: cellType(:)

        integer :: idx,ii

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% get data that can be extracted from links
        call initial_conditions_from_linkdata &
            (elem2R, elem2I, elemMR, elemMI, elem2YN, elemMYN, linkR, linkI)

        !% custom initial condition setup for special test cases
        call custom_initial_condition &
            (elem2R, elem2I, elemMR, elemMI, elem2YN, elemMYN, bcdataDn)

        call initial_junction_conditions &
           (faceR, faceI, elem2R, elem2I, elemMR, elemMYN, elemMI, nodeR, nodeI)

        call initial_storage_conditions &
            (faceR, faceI, elem2R, elem2I, elemMR, elemMI, nodeR, nodeI)

        !% set the bc elements (outside of face) to null values
        call bc_nullify_ghost_elem (elem2R, bcdataDn)
        call bc_nullify_ghost_elem (elem2R, bcdataUp)

        !% update the geometry
        call element_geometry_update &
            (elem2R, elem2I, elem2YN, e2r_Volume, &
            elemMR, elemMI, elemMYN, eMr_Volume, &
            faceR, faceI, bcdataDn, bcdataUp, thisTime, 0, &
            ID, numberPairs, ManningsN, Length, zBottom, xDistance, &
            Breadth, widthDepthData, cellType)
                    
        call meta_element_assign &
            (elem2I, e2i_elem_type, e2i_meta_elem_type)

        call meta_element_assign &
            (elemMI, eMi_elem_type, eMi_meta_elem_type)

        call element_dynamics_update &
            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
            bcdataDn, bcdataUp, e2r_Velocity, eMr_Velocity, &
            e2r_Volume, eMr_Volume, e2r_Flowrate, eMr_Flowrate, thisTime)

        call face_meta_element_assign (faceI, elem2I, N_face)

        call face_update &
            (elem2R, elem2I, elemMR, faceR, faceI, faceYN, &
            bcdataDn, bcdataUp, e2r_Velocity, eMr_Velocity,  &
            e2r_Volume, eMr_Volume, thisTime, 0)

        !% set the initial condition for Qonly element
        call QonlyElem_initial_condition_setup &
            (elem2R, elemMR, faceR, elem2I, elemMI, faceI, elem2YN, elemMYN, faceYN)

        !% select the solver based on Area/AreaFull
        call initial_solver_select &
            (elem2R, elemMR, elem2I, elemMI)

        !% set the element-specific smallvolume value
        !% HACK - THIS IS ONLY FOR RECTANGULAR ELEMENTS
        !% HACK - OTHER GEOMETRY TYPE NEEDED
        if (setting%SmallVolume%UseSmallVolumes) then
            elem2R(:,e2r_SmallVolume) = 0.01
            elemMR(:,eMr_SmallVolume) = 0.01
            where (elem2I(:,e2i_geometry) == eRectangular)
                elem2R(:,e2r_SmallVolume) = setting%SmallVolume%DepthCutoff * elem2R(:,e2r_BreadthScale) * &
                    elem2R(:,e2r_Length)
                !% HACK - Current version only allows junction to be rectangular
            elsewhere (elemMI(:,eMi_geometry) == eRectangular)
                elemMR(:,eMr_SmallVolume) = setting%SmallVolume%DepthCutoff * elemMR(:,eMr_BreadthScale) * &
                    elemMR(:,eMr_Length)
            elsewhere (elem2I(:,e2i_geometry) == eTriangular)
                elem2R(:,e2r_SmallVolume) = onehalfR * setting%SmallVolume%DepthCutoff * elem2R(:,e2r_BreadthScale) * &
                    elem2R(:,e2r_Length)
            elsewhere (elem2I(:,e2i_geometry) == eTrapezoidal)
                elem2R(:,e2r_SmallVolume) = (elem2R(:,e2r_BreadthScale) + onehalfR * (elem2R(:,e2r_LeftSlope) + &
                    elem2R(:,e2r_RightSlope)) * setting%SmallVolume%DepthCutoff) * &
                    setting%SmallVolume%DepthCutoff
            endwhere
        else
            elem2R(:,e2r_SmallVolume) = zeroR
            elemMR(:,eMr_SmallVolume) = zeroR
        endif

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine initial_condition_setup
    !
    !==========================================================================
    ! PRIVATE BELOW HERE
    !==========================================================================
    !
    subroutine initial_conditions_from_linkdata &
        (elem2R, elem2I, elemMR, elemMI, elem2YN, elemMYN, linkR, linkI)
        !
        ! The link data structure can store a variety of geometric data.
        ! This will be expanded in the future
        !
        character(64) :: subroutine_name = 'initial_conditions_from_linkdata'

        real,      intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        integer,   intent(in out)  :: elem2I(:,:),  elemMI(:,:)
        logical,   intent(in out)  :: elem2YN(:,:), elemMYN(:,:)

        real,      target,   intent(in)      :: linkR(:,:)
        integer,   target,   intent(in)      :: linkI(:,:)

        real               :: kappa
        real,      pointer :: dup, ddn
        integer,   pointer :: Lindx, LdepthType
        integer :: ii, ei_max, mm, nn

        real :: trapz_tanTheta, CC, BB, AoverAfull, YoverYfull

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        do ii=1,N_link
            Lindx      => linkI(ii,li_idx)
            LdepthType => linkI(ii,li_InitialDepthType)

            !% up and downstream depths on this link
            dup => linkR(ii,lr_InitialUpstreamDepth)
            ddn => linkR(ii,lr_InitialDnstreamDepth)


            select case (LdepthType)

              case (1)

                !%  Initial depth --------------------------------------------------
                if (linkR(ii,lr_InitialDepth) /= nullvalueR) then
                    !%  if the link has a uniform depth as an initial condition
                    where (elem2I(:,e2i_link_ID) == Lindx)
                        elem2R(:,e2r_Depth) = linkR(ii,lr_InitialDepth)
                    endwhere
                else
                    where (elem2I(:,e2i_link_ID) == Lindx)
                        elem2R(:,e2r_Depth) = 0.5*(dup + ddn)
                    endwhere
                endif
              case (2)
                !% if the link has linearly-varying depth
                !% depth at the downstream element (link position =1)
                where ( (elem2I(:,e2i_link_Pos) == 1) .and. (elem2I(:,e2i_link_ID) == Lindx) )
                    elem2R(:,e2r_Depth) = ddn
                endwhere

                !%  using a linear distribution over the links
                ei_max = maxval(elem2I(:,e2i_link_Pos),1,elem2I(:,e2i_link_ID) == Lindx)
                do mm=2,ei_max
                    !% find the element that is at the mm position in the link
                    where ( (elem2I(:,e2i_link_Pos) == mm) .and. (elem2I(:,e2i_link_ID) == Lindx) )
                        ! use a linear interp
                        elem2R(:,e2r_Depth) = ddn + (dup - ddn) * real(mm-1) / real(ei_max-1)
                    endwhere
                end do

              case (3)
                ! HACK - this needs work to make the exponent driven by dup and ddn
                ! to ensure decay to dup over length

                !% exponential decay


                !% if the link has linearly-varying depth
                !% depth at the downstream element (link position =1)
                where ( (elem2I(:,e2i_link_Pos) == 1) .and. (elem2I(:,e2i_link_ID) == Lindx) )
                    elem2R(:,e2r_Depth) = ddn
                endwhere

                !%  using a linear distribution over the links
                ei_max = maxval(elem2I(:,e2i_link_Pos),1,elem2I(:,e2i_link_ID) == Lindx)

                do mm=2,ei_max
                    kappa = real(ei_max-1)
                    !% find the element that is at the mm position in the link
                    if (ddn - dup > zeroR) then
                        !%  depth decreases exponentially going upstream
                        where ( (elem2I(:,e2i_link_Pos) == mm) .and. &
                            (elem2I(:,e2i_link_ID) == Lindx)        )
                            elem2R(:,e2r_Depth) = (ddn-dup) * exp(-real(mm-1)) + dup
                        endwhere
                    elseif (ddn - dup < zeroR) then
                        !%  depth increases exponentially going upstream
                        where ( (elem2I(:,e2i_link_Pos) == mm) .and. &
                            (elem2I(:,e2i_link_ID) == Lindx)        )
                            elem2R(:,e2r_Depth) = dup - (dup-ddn) * exp(-real(mm-1))
                        endwhere
                    else
                        !%  uniform depth
                        where ( (elem2I(:,e2i_link_Pos) == mm) .and. &
                            (elem2I(:,e2i_link_ID) == Lindx)        )
                            ! use a linear interp
                            elem2R(:,e2r_Depth) = ddn
                        endwhere

                    endif

                end do
            end select

            !%  handle all the initial conditions that don't depend on geometry type
            !%
            where (elem2I(:,e2i_link_ID) == Lindx)
                elem2I(:,e2i_roughness_type) = linkI(ii,li_roughness_type)
                elem2R(:,e2r_Roughness)      = linkR(ii,lr_Roughness)
                elem2R(:,e2r_Flowrate)       = linkR(ii,lr_InitialFlowrate)
                ! elem2R(:,e2r_Flowrate_N0)    = elem2R(:,e2r_Flowrate)
                ! elem2R(:,e2r_Flowrate_N1)    = elem2R(:,e2r_Flowrate) 
                elem2R(:,e2r_LeftSlope)      = linkR(ii,lr_LeftSlope)
                elem2R(:,e2r_RightSlope)     = linkR(ii,lr_RightSlope)
                elem2R(:,e2r_ParabolaValue)  = linkR(ii,lr_ParabolaValue)
            endwhere

            if (linkI(ii,li_geometry) == lRectangular ) then
                !% handle rectangular elements
                where (elem2I(:,e2i_link_ID) == Lindx)
                    elem2I(:,e2i_geometry)  = eRectangular
                    elem2R(:,e2r_HydDepth)  = elem2R(:,e2r_Depth)
                    elem2R(:,e2r_FullDepth) = linkR(ii,lr_FUllDepth)
                    elem2R(:,e2r_BreadthScale)   = linkR(ii,lr_BreadthScale)
                    elem2R(:,e2r_Topwidth)  = linkR(ii,lr_BreadthScale)
                    elem2R(:,e2r_Eta)       = elem2R(:,e2r_Zbottom)  + elem2R(:,e2r_HydDepth)
                    elem2R(:,e2r_Area)      = elem2R(:,e2r_HydDepth) * elem2R(:,e2r_BreadthScale)
                    elem2R(:,e2r_FullArea)  = elem2R(:,e2r_FullDepth) * elem2R(:,e2r_BreadthScale)
                    elem2R(:,e2r_Area_N0)   = elem2R(:,e2r_Area)
                    elem2R(:,e2r_Area_N1)   = elem2R(:,e2r_Area)
                    elem2R(:,e2r_Volume)    = elem2R(:,e2r_Area)     * elem2R(:,e2r_Length)
                    elem2R(:,e2r_Perimeter) = elem2R(:,e2r_BreadthScale) + twoR * elem2R(:,e2r_HydDepth)
                    elem2R(:,e2r_Zcrown)    = elem2R(:,e2r_Zbottom) + elem2R(:,e2r_FullDepth)
                endwhere

            elseif (linkI(ii,li_geometry) == lParabolic ) then
                !% handle parabolic elements
                ! Input Topwidth, InitialDepth
                where (elem2I(:,e2i_link_ID) == Lindx)
                    elem2I(:,e2i_geometry)  = eParabolic
                    elem2R(:,e2r_FullDepth) = linkR(ii,lr_FUllDepth)
                    ! calculate hyd depth from the depth
                    elem2R(:,e2r_HydDepth) = twothirdR * elem2R(:,e2r_Depth)

                    elem2R(:,e2r_BreadthScale) = zeroR

                    elem2R(:,e2r_Topwidth)  = twoR &
                        * sqrt(elem2R(:,e2r_Depth)/elem2R(:,e2r_ParabolaValue))

                    elem2R(:,e2r_Area)      = twothirdR * elem2R(:,e2r_Depth) &
                        * elem2R(:,e2r_Topwidth)

                    elem2R(:,e2r_FullArea)  = twothirdR * elem2R(:,e2r_FullDepth) &
                        * twoR * sqrt(elem2R(:,e2r_FullDepth)/elem2R(:,e2r_ParabolaValue))

                    elem2R(:,e2r_Area_N0)   = elem2R(:,e2r_Area)

                    elem2R(:,e2r_Area_N1)   = elem2R(:,e2r_Area)

                    elem2R(:,e2r_Perimeter) = onehalfR * elem2R(:,e2r_Topwidth) &
                        *( &
                        sqrt &
                        ( &
                        oneR  &
                        + (fourR  &
                        * elem2R(:,e2r_Depth)/elem2R(:,e2r_Topwidth))**twoR &
                        )  &
                        + (elem2R(:,e2r_Topwidth)/fourR * elem2R(:,e2r_Depth)) &
                        *log &
                        ( &
                        fourR * elem2R(:,e2r_Depth)/elem2R(:,e2r_Topwidth)  &
                        + sqrt &
                        ( &
                        oneR  &
                        + (fourR  &
                        * elem2R(:,e2r_Depth)/elem2R(:,e2r_Topwidth))**twoR &
                        ) &
                        )  &
                        )

                    elem2R(:,e2r_Eta)       = elem2R(:,e2r_Zbottom)                &
                        + elem2R(:,e2r_HydDepth)

                    elem2R(:,e2r_Volume)    = elem2R(:,e2r_Area) &
                        * elem2R(:,e2r_Length)

                    elem2R(:,e2r_Zcrown)    = elem2R(:,e2r_Zbottom) + elem2R(:,e2r_FullDepth)

                endwhere

            elseif (linkI(ii,li_geometry) == lTrapezoidal ) then
                !% handle trapezoidal elements
                ! Input: Left Slope, Right Slope, Bottom Width, InitialDepth
                where (elem2I(:,e2i_link_ID) == Lindx)
                    elem2I(:,e2i_geometry)  = eTrapezoidal

                    elem2R(:,e2r_BreadthScale) = linkR(ii,lr_BreadthScale)
                    elem2R(:,e2r_FullDepth)    = linkR(ii,lr_FUllDepth)

                    ! (Bottom width + averageSlope * hydraulicDepth)*hydraulicDepth
                    elem2R(:,e2r_Area)      = (elem2R(:,e2r_BreadthScale)           &
                        + onehalfR &
                        * (elem2R(:,e2r_LeftSlope) + elem2R(:,e2r_RightSlope)) &
                        * elem2R(:,e2r_Depth)) * elem2R(:,e2r_Depth)

                    elem2R(:,e2r_FullArea)  = (elem2R(:,e2r_BreadthScale)           &
                        + onehalfR &
                        * (elem2R(:,e2r_LeftSlope) + elem2R(:,e2r_RightSlope)) &
                        * elem2R(:,e2r_FullDepth)) * elem2R(:,e2r_FullDepth)

                    elem2R(:,e2r_Area_N0)   = elem2R(:,e2r_Area)

                    elem2R(:,e2r_Area_N1)   = elem2R(:,e2r_Area)

                    ! Bottom width + (lslope + rslope) * hydraulicDepth
                    elem2R(:,e2r_Topwidth)  = elem2R(:,e2r_BreadthScale)            &
                        + elem2R(:,e2r_Depth)                                   &
                        * (elem2R(:,e2r_LeftSlope) + elem2R(:,e2r_RightSlope))

                    elem2R(:,e2r_HydDepth) = elem2R(:,e2r_Area) / elem2R(:,e2r_Topwidth)

                    elem2R(:,e2r_Eta)       = elem2R(:,e2r_Zbottom)                &
                        + elem2R(:,e2r_HydDepth)

                    elem2R(:,e2r_Volume)    = elem2R(:,e2r_Area) &
                        * elem2R(:,e2r_Length)

                    ! Bottom width + hydraulicDepth*lengthSidewall
                    elem2R(:,e2r_Perimeter) = elem2R(:,e2r_BreadthScale) &
                        + elem2R(:,e2r_Depth) &
                        * (sqrt(oneR + elem2R(:,e2r_LeftSlope)**twoR) &
                        + sqrt(oneR + elem2R(:,e2r_RightSlope)**twoR))

                    elem2R(:,e2r_Zcrown)    = elem2R(:,e2r_Zbottom) + elem2R(:,e2r_FullDepth)
                endwhere
            elseif (linkI(ii,li_geometry) == lCircular ) then
                !% handle circular elements
    
                do nn=1, size(elem2I(:,e2i_link_ID),1)
                    if (elem2I(nn,e2i_link_ID) == Lindx) then

                        elem2I(nn,e2i_geometry)     = eCircular
                        
                        elem2R(nn,e2r_BreadthScale) = linkR(ii,lr_FUllDepth)

                        elem2R(nn,e2r_FullDepth) = linkR(ii,lr_FUllDepth)
                        
                        elem2R(nn,e2r_Radius)   = linkR(ii,lr_FUllDepth) / twoR
                        
                        elem2R(nn,e2r_Zcrown)   = elem2R(nn,e2r_Zbottom) + elem2R(nn,e2r_FullDepth)

                        elem2R(nn,e2r_FullArea) = onefourthR * pi * linkR(ii,lr_FUllDepth)  ** twoR

                        elem2R(nn,e2r_Eta)      = elem2R(nn,e2r_Zbottom) + elem2R(nn,e2r_Depth)

                        ! Set surcharge condition if the pipe is full
                        if (elem2R(nn,e2r_Eta) .GE. elem2R(nn,e2r_Zcrown)) then
                            elem2YN(nn,e2YN_IsSurcharged) = .true.
                            elem2R(nn,e2r_Area) = elem2R(nn,e2r_FullArea)
                            elem2R(nn,e2r_HydDepth) = zeroR  !this is the modified hydralic depth for pipe
                            elem2R(nn,e2r_Topwidth) = zeroR
                            elem2R(nn,e2r_HydRadius) = onefourthR * elem2R(nn,e2r_FullDepth)
                            YoverYfull = oneR
                            AoverAfull = oneR
                        else
                            elem2YN(nn,e2YN_IsSurcharged) = .false.
                            YoverYfull = elem2R(nn,e2r_Depth) / elem2R(nn,e2r_FullDepth)
                            elem2R(nn,e2r_Area)     = elem2R(nn,e2r_FullArea) * table_lookup(YoverYfull, ACirc, NACirc)
                            elem2R(nn,e2r_Topwidth) = elem2R(nn,e2r_FullDepth)* table_lookup(YoverYfull, WCirc, NWCirc)
                            elem2R(nn,e2r_HydRadius) = onefourthR * elem2R(nn,e2r_FullDepth) * &
                                                        table_lookup(YoverYfull, RCirc, NRCirc)
                            AoverAfull = elem2R(nn,e2r_Area) / elem2R(nn,e2r_FullArea)
                            if (AoverAfull .GT. onehalfR) then
                                elem2R(nn,e2r_HydDepth)   = elem2R(nn,e2r_Eta) - elem2R(nn,e2r_Zbottom) + &
                                                                elem2R(nn,e2r_Radius) * (onefourthR * pi - oneR)
                            elseif (AoverAfull .LE. onehalfR) then
                                elem2R(nn,e2r_HydDepth)   = max(elem2R(nn,e2r_Area)/ elem2R(nn,e2r_Topwidth), zeroR)
                            endif

                        endif

                        elem2R(nn,e2r_Area_N0)   = elem2R(nn,e2r_Area)
                        elem2R(nn,e2r_Area_N1)   = elem2R(nn,e2r_Area)
                        elem2R(nn,e2r_Volume)    = elem2R(nn,e2r_Area) * elem2R(nn,e2r_Length)
                        elem2R(nn,e2r_Perimeter) = elem2R(nn,e2r_Area) / elem2R(nn,e2r_HydRadius)
                    endif
                enddo

            elseif (linkI(ii,li_geometry) == lTriangular ) then
                !% handle triangle elements
                ! Input: Left Slope, Right Slope, InitialDepth
                where (elem2I(:,e2i_link_ID) == Lindx)
                    elem2I(:,e2i_geometry)  = eTriangular

                    elem2R(:,e2r_HydDepth)  = onehalfR * elem2R(:,e2r_Depth)
                    elem2R(:,e2r_FullDepth) = linkR(ii,lr_FUllDepth)
                    elem2R(:,e2r_BreadthScale) = elem2R(:,e2r_FullDepth)              &
                        * (elem2R(:,e2r_LeftSlope) + elem2R(:,e2r_RightSlope))
                    ! (averageSlope * hydraulicDepth)*hydraulicDepth
                    elem2R(:,e2r_Area) = onehalfR &
                        * (elem2R(:,e2r_LeftSlope) + elem2R(:,e2r_RightSlope)) &
                        * elem2R(:,e2r_Depth) * elem2R(:,e2r_Depth)

                    elem2R(:,e2r_FullArea) = onehalfR &
                        * (elem2R(:,e2r_LeftSlope) + elem2R(:,e2r_RightSlope)) &
                        * elem2R(:,e2r_FullDepth) * elem2R(:,e2r_FullDepth)

                    elem2R(:,e2r_Area_N0)   = elem2R(:,e2r_Area)

                    elem2R(:,e2r_Area_N1)   = elem2R(:,e2r_Area)

                    ! (lslope + rslope) * hydraulicDepth
                    elem2R(:,e2r_Topwidth) = elem2R(:,e2r_Depth)               &
                        * (elem2R(:,e2r_LeftSlope) + elem2R(:,e2r_RightSlope))

                    elem2R(:,e2r_Eta) = elem2R(:,e2r_Zbottom)                &
                        + elem2R(:,e2r_HydDepth)

                    elem2R(:,e2r_Volume) = elem2R(:,e2r_Area) &
                        * elem2R(:,e2r_Length)

                    ! hydraulicDepth*lengthSidewall
                    elem2R(:,e2r_Perimeter) = elem2R(:,e2r_Depth) &
                        * (sqrt(oneR + elem2R(:,e2r_LeftSlope)**twoR) &
                        + sqrt(oneR + elem2R(:,e2r_RightSlope)**twoR))
                    elem2R(:,e2r_Zcrown)    = elem2R(:,e2r_Zbottom) + elem2R(:,e2r_FullDepth)
                endwhere

            elseif (linkI(ii,li_geometry) == lWidthDepth ) then
                !% handle width-depth elements

                where (elem2I(:,e2i_link_ID) == Lindx)

                    elem2I(:,e2i_geometry)  = eWidthDepth
                    elem2R(:,e2r_HydDepth)  = elem2R(:,e2r_Depth)
                    elem2R(:,e2r_FullDepth) = linkR(ii,lr_FUllDepth)
                    elem2R(:,e2r_BreadthScale)   = linkR(ii,lr_BreadthScale)
                    elem2R(:,e2r_Topwidth)  = linkR(ii,lr_TopWidth)
                    elem2R(:,e2r_Eta)       = elem2R(:,e2r_Zbottom)  + elem2R(:,e2r_HydDepth)
                    elem2R(:,e2r_Area)      = elem2R(:,e2r_Topwidth) * elem2R(:,e2r_HydDepth)
                    elem2R(:,e2r_FullArea)  = elem2R(:,e2r_Topwidth) * elem2R(:,e2r_FullDepth)
                    elem2R(:,e2r_Area_N0)   = elem2R(:,e2r_Area)
                    elem2R(:,e2r_Area_N1)   = elem2R(:,e2r_Area)
                    elem2R(:,e2r_Volume)    = elem2R(:,e2r_Area) * elem2R(:,e2r_Length)
                    elem2R(:,e2r_Perimeter) = onehalfR * elem2R(:,e2r_Area) / elem2R(:,e2r_HydDepth)
                    elem2R(:,e2r_Zcrown)    = elem2R(:,e2r_Zbottom) + elem2R(:,e2r_FullDepth)
                endwhere
            else
                !% handle elements of other geometry types
                print *, 'error: initialization for non-defined elements needed in ',subroutine_name
                stop
            end if

            !%  Find if an element is surcharged. Set topwidth to zero for surcharged condition
            where ( (elem2I(:,e2i_elem_type) == ePipe) .and. &
                    (elem2R(:,e2r_Eta) .GE. elem2R(:,e2r_Zcrown)) )
                elem2YN(:,e2YN_IsSurcharged) = .true.
                elem2R(:,e2r_Topwidth) = zeroR 
            elsewhere ((elem2I(:,e2i_elem_type) == ePipe) .and. &
                       (elem2R(:,e2r_Eta) .LE. elem2R(:,e2r_Zcrown)) )
                elem2YN(:,e2YN_IsSurcharged) = .false.
            endwhere

            !%  Setting provisional geometry for weir and orifice element to correctly interpolate to faces
            !%  Also, setting up 
            where ( elem2I(:,e2i_meta_elem_type) == eQonly )
                elem2R(:,e2r_HydDepth)      = 1.0e-7
                elem2R(:,e2r_Topwidth)      = 1.0e-7
                elem2R(:,e2r_Eta)           = 1.0e-7
                elem2R(:,e2r_Area)          = 1.0e-7
                elem2R(:,e2r_Volume)        = 1.0e-7
                elem2R(:,e2r_Perimeter)     = 1.0e-7
                elem2R(:,e2r_Depth)         = 1.0e-7
                elem2R(:,e2r_Area_N0)       = elem2R(:,e2r_Area)
                elem2R(:,e2r_Area_N1)       = elem2R(:,e2r_Area)
            endwhere

            if (setting%BCondition%InflowRampup) then
                elem2R(:,e2r_Flowrate) = setting%BCondition%flowrateIC
            endif

            !%  Update velocity
            where (  (elem2I(:,e2i_link_ID) == Lindx) .and. (elem2R(:,e2r_Area) > zeroR) )
                elem2R(:,e2r_Velocity)  = elem2R(:,e2r_Flowrate) / elem2R(:,e2r_Area)
            endwhere

        enddo

        

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine initial_conditions_from_linkdata
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine custom_initial_condition &
        (elem2R, elem2I, elemMR, elemMI, elem2YN, elemMYN, bcdataDn)
    !    
    ! set up custom initial conditions for special cases like flow over a bump
    ! hard coded custon initial condition setup
    !
    character(64) :: subroutine_name = 'custom_initial_condition'

    real,           intent(in out)  :: elem2R(:,:),  elemMR(:,:)
    integer,        intent(in out)  :: elem2I(:,:),  elemMI(:,:)
    logical,        intent(in out)  :: elem2YN(:,:), elemMYN(:,:)

    type(bcType),   intent(in)      :: bcdataDn(:)
    real    :: thisVal
    integer :: ii

    !--------------------------------------------------------------------------
    if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        if (setting%CustomIC%UseCustomInitialCondition) then

            select case (setting%TestCase%TestName)

                case ('swashes_007')
                    !% there are problem with zbottom calculation for elements and reproducing
                    !% the same zbottom as SvePy.
                    !% the zbottom is hardcoded to make it consistant with the SvePy
                    !% find the x value for elements. Later move to a subroutine
                    !% hard coded for swashes test case
                    do ii = 2,N_elem2 -1 
                        elem2R(ii,e2r_X) = sum(elem2R(2:N_elem2 -1,e2r_Length)) - sum(elem2R(2:ii,e2r_Length)) &
                                              + elem2R(ii,e2r_Length)/2.0
                                              elem2R(ii,e2r_Zbottom) = (0.2 - 0.05 * ((elem2R(ii,e2r_X) - 10.0) ** 2.0)) 
                        if (elem2R(ii,e2r_X) < 8.0 ) then
                            elem2R(ii,e2r_Zbottom) = 0.0
                        elseif (elem2R(ii,e2r_X) > 12.0 ) then
                            elem2R(ii,e2r_Zbottom) = 0.0
                        endif
                    enddo

                    ! get the boundary eta and set that for every other elements
                    thisVal = elem2R(1,e2r_eta)
                    elem2R(:,e2r_Eta)       = thisVal
                    elem2R(:,e2r_Depth)     = elem2R(:,e2r_Eta) - elem2R(:,e2r_Zbottom)
                    elem2R(:,e2r_HydDepth)  = elem2R(:,e2r_Depth)
                    elem2R(:,e2r_Area)      = elem2R(:,e2r_HydDepth) * elem2R(:,e2r_BreadthScale)
                    elem2R(:,e2r_FullArea)  = elem2R(:,e2r_FullDepth) * elem2R(:,e2r_BreadthScale)
                    elem2R(:,e2r_Area_N0)   = elem2R(:,e2r_Area)
                    elem2R(:,e2r_Area_N1)   = elem2R(:,e2r_Area)
                    elem2R(:,e2r_Volume)    = elem2R(:,e2r_Area) * elem2R(:,e2r_Length)
                    elem2R(:,e2r_Perimeter) = elem2R(:,e2r_BreadthScale) + twoR * elem2R(:,e2r_HydDepth)
                    elem2R(:,e2r_Zcrown)    = elem2R(:,e2r_Zbottom) + elem2R(:,e2r_FullDepth)
                    elem2R(:,e2r_Velocity)  = elem2r(:,e2r_Flowrate) / elem2R(:,e2r_Area) 

                case ('simple_pipe_006')
                !% This custom initial contion propagates a wave through the pipe
                    elem2R(90:101,e2r_Eta) = elem2R(90:101,e2r_Eta) + 5.0
                    elem2R(:,e2r_Depth)     = elem2R(:,e2r_Eta) - elem2R(:,e2r_Zbottom)
                    elem2R(:,e2r_HydDepth)  = elem2R(:,e2r_Depth)
                    elem2R(:,e2r_Area)      = elem2R(:,e2r_HydDepth) * elem2R(:,e2r_BreadthScale)
                    elem2R(:,e2r_FullArea)  = elem2R(:,e2r_FullDepth) * elem2R(:,e2r_BreadthScale)
                    elem2R(:,e2r_Area_N0)   = elem2R(:,e2r_Area)
                    elem2R(:,e2r_Area_N1)   = elem2R(:,e2r_Area)
                    elem2R(:,e2r_Volume)    = elem2R(:,e2r_Area) * elem2R(:,e2r_Length)
                    elem2R(:,e2r_Perimeter) = elem2R(:,e2r_BreadthScale) + twoR * elem2R(:,e2r_HydDepth)
                    elem2R(:,e2r_Zcrown)    = elem2R(:,e2r_Zbottom) + elem2R(:,e2r_FullDepth)
                    elem2R(:,e2r_Velocity)  = elem2r(:,e2r_Flowrate) / elem2R(:,e2r_Area)
                case default
                print *, 'Warning '
                print *, setting%TestCase%TestName, ' does not need a custom initial condition'
            end select
        endif

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine custom_initial_condition
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine initial_junction_conditions &
        (faceR, faceI, elem2R, elem2I, elemMR, elemMYN, elemMI, nodeR, nodeI)

        character(64) :: subroutine_name = 'initial_junction_conditions'

        real,              intent(in out)  :: elemMR(:,:)
        logical,           intent(in out)  :: elemMYN(:,:)
        real,      target, intent(in)      :: elem2R(:,:), nodeR(:,:), faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:), elemMI(:,:), nodeI(:,:), faceI(:,:)

        integer,   pointer :: tface, telem

        real   :: upvalue(upstream_face_per_elemM), dnvalue(dnstream_face_per_elemM)

        integer :: ii, mm
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        if (N_elemM > 0) then
            !% initialize the free surface from the average of the adjacent elements
            call junction_adjacent_element_average &
                (elem2R, elemMR, elemMI, faceI, e2r_Eta, eMr_Eta)

            !% initialize the branch areas to the values of the adjacent elements
            call junction_adjacent_element_values_to_branches &
                (elem2R, elemMR, elemMI, faceI, e2r_Area, eMr_AreaUp, eMr_AreaDn)

            !% initialize the branch flowrates to the values of the adjacent elements
            call junction_adjacent_element_values_to_branches &
                (elem2R, elemMR, elemMI, faceI, e2r_Flowrate, eMr_FlowrateUp, eMr_FlowrateDn)

            !% initialize element momentum to the average of the net upstream and downstrea
            !% fluxes
            call junction_branch_average_of_inflows_and_outflows (elemMR, elemMI)

            !% here we assume the branch and junction topwidths are already initialized
            !% in a prior call to junction_geometry_setup
            !print *, elemMR(:,eMr_Topwidth)
            !print *, elemMR(:,eMr_TopwidthAll)

            where (elemMI(:,eMi_elem_type) == eJunctionChannel)
                elemMR(:,eMr_HydDepth) = elemMR(:,eMr_Eta) - elemMR(:,eMr_Zbottom)
            endwhere

            ! HACK -- need other geometry types

            where ((elemMI(:,eMi_geometry) == eRectangular) .and. &
                (elemMI(:,eMi_elem_type) == eJunctionChannel))
                elemMR(:,eMr_Area)      = elemMR(:,eMr_HydDepth) * elemMR(:,eMr_Topwidth)
                elemMR(:,eMr_Volume)    = elemMR(:,eMr_Area)     * elemMR(:,eMr_Length)
                elemMR(:,eMr_Perimeter) = elemMR(:,eMr_Breadthscale) + twoR * elemMR(:,eMr_HydDepth)
                elemMR(:,eMr_HydRadius) = elemMR(:,eMr_Area) / elemMR(:,eMr_Perimeter)
            endwhere

            where (elemMR(:,eMr_Eta) .GE. elem2R(:,eMr_Zcrown))
                elemMYN(:,eMYN_IsSurcharged) = .true.
                elemMR(:,eMr_Topwidth) = zeroR 
            elsewhere (elemMR(:,eMr_Eta) .LE. elem2R(:,eMr_Zcrown))
                elemMYN(:,eMYN_IsSurcharged) = .false.
            endwhere

            !% velocities
            call junction_branch_velocities (elemMR, elemMI)

            where ((elemMI(:,eMi_elem_type) == eJunctionChannel) .and. &
                (elemMr(:,eMr_Area) > zeroR))
                elemMR(:,eMr_Velocity) = elemMR(:,eMr_Flowrate) / elemMR(:,eMr_Area)
            endwhere

        end if


        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine initial_junction_conditions
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine initial_storage_conditions &
        (faceR, faceI, elem2R, elem2I, elemMR, elemMI, nodeR, nodeI)

        character(64) :: subroutine_name = 'initial_storage_conditions'

        real,              intent(in out)  :: elemMR(:,:)
        real,      target, intent(in)      :: elem2R(:,:), nodeR(:,:), faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:), elemMI(:,:), nodeI(:,:), faceI(:,:)

        integer,   pointer :: tface, telem

        real   :: upvalue(upstream_face_per_elemM), dnvalue(dnstream_face_per_elemM)

        integer,   dimension(4)    :: e2rset, eMrset

        integer :: ii, mm
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        if (N_elemM > 0) then

            e2rset = (/e2r_Eta, e2r_Topwidth,   e2r_Area,  e2r_Perimeter/)
            eMrset = (/eMr_Eta, eMr_Topwidth,   eMr_Area,  eMr_Perimeter/)

            do mm=1,size(e2rset)

                !% initialize the the average of the adjacent elements
                call storage_adjacent_element_average &
                    (elem2R, elemMR, elemMI, faceI, e2rset(mm), eMrset(mm))

            end do

            call storage_initialize_depth_volume (elemMR, elemMI)

            where (elemMI(:,eMi_elem_type) == eStorage)
                !% setting the flowrate in Storafe unit as zero to initialize.
                !% the flowrate in a storage element does not make any difference.
                !% because this flowrate will not be interpolated to the faces.
                !% Talk with Dr. Hodges about this matter
                elemMR(:,eMr_Flowrate) = 1.0E-7
                elemMR(:,eMr_Velocity) = 1.0E-7
            endwhere

        end if

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine initial_storage_conditions
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine meta_element_assign (elemI, ei_elem_type, ei_meta_elem_type)
        !
        ! Assign meta element type to elements
        !

        character(64) :: subroutine_name = 'meta_element_assign'

        integer,   target,     intent(inout)    :: elemI(:,:)
        integer,               intent(in)       :: ei_elem_type

        integer,               intent(in)       :: ei_meta_elem_type


        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        where ( (elemI(:,ei_elem_type) == eChannel)  .or. &
                (elemI(:,ei_elem_type) == ePipe) )

            elemI(:,ei_meta_elem_type) = eHQ2

        elsewhere ( (elemI(:,ei_elem_type) == eJunctionChannel) .or. &
                    (elemI(:,ei_elem_type) == eJunctionPipe)  )

            elemI(:,ei_meta_elem_type) = eHQM

        elsewhere ( (elemI(:,ei_elem_type) == eWeir)    .or. &
                    (elemI(:,ei_elem_type) == eorifice) .or. &
                    (elemI(:,ei_elem_type) == ePump)  )

            elemI(:,ei_meta_elem_type) = eQonly

        elsewhere ( (elemI(:,ei_elem_type) == eStorage) )

            elemI(:,ei_meta_elem_type) = eHonly

        elsewhere ( (elemI(:,ei_elem_type) == eBCup) .or. &
            (elemI(:,ei_elem_type) == eBCdn)  )
            ! Assigning nonHQ meta elem type to boundary conditions. Confirm this!
            elemI(:,ei_meta_elem_type) = eNonHQ
        end where

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine meta_element_assign
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine face_meta_element_assign (faceI, elemI, N_face)

        character(64) :: subroutine_name = 'face_meta_element_assign'

        integer,      target,     intent(in out)  :: faceI(:,:), elemI(:,:)
        integer,                  intent(in)      :: N_face

        integer :: ii

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        do ii=1, N_face
            if ( (faceI(ii,fi_etype_u) == eChannel) .or. &
                (faceI(ii,fi_etype_u) == ePipe) ) then

                faceI(ii,fi_meta_etype_u) = eHQ2

            elseif ( (faceI(ii,fi_etype_u) == eJunctionChannel) .or. &
                (faceI(ii,fi_etype_u) == eJunctionPipe) ) then

                faceI(ii,fi_meta_etype_u) = eHQm

            elseif ( (faceI(ii,fi_etype_u) == eWeir)    .or. &
                (faceI(ii,fi_etype_u) == eorifice) .or. &
                (faceI(ii,fi_etype_u) == ePump) ) then

                faceI(ii,fi_meta_etype_u) = eQonly

            elseif ( (faceI(ii,fi_etype_u) == eStorage) ) then

                faceI(ii,fi_meta_etype_u) = eHonly

            elseif ( (faceI(ii,fi_etype_u) == eBCdn)    .or. &
                (faceI(ii,fi_etype_u) == eBCup) ) then

                faceI(ii,fi_meta_etype_u) = eNonHQ

            else
                print*, 'undefined element type upstream of face', ii
                stop
            endif
        end do

        do ii=1, N_face
            if ( (faceI(ii,fi_etype_d) == eChannel) .or. &
                (faceI(ii,fi_etype_d) == ePipe) ) then

                faceI(ii,fi_meta_etype_d) = eHQ2

            elseif ( (faceI(ii,fi_etype_d) == eJunctionChannel) .or. &
                (faceI(ii,fi_etype_d) == eJunctionPipe) ) then

                faceI(ii,fi_meta_etype_d) = eHQm

            elseif ( (faceI(ii,fi_etype_d) == eWeir)    .or. &
                (faceI(ii,fi_etype_d) == eorifice) .or. &
                (faceI(ii,fi_etype_d) == ePump) ) then

                faceI(ii,fi_meta_etype_d) = eQonly

            elseif ( (faceI(ii,fi_etype_d) == eStorage) ) then

                faceI(ii,fi_meta_etype_d) = eHonly

            elseif ( (faceI(ii,fi_etype_d) == eBCdn)    .or. &
                (faceI(ii,fi_etype_d) == eBCup) ) then

                faceI(ii,fi_meta_etype_d) = eNonHQ

            else
                print*, 'undefined element type downstream of face', ii
                stop
            endif
        end do

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine face_meta_element_assign
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine QonlyElem_initial_condition_setup &
        (elem2R, elemMR, faceR, elem2I, elemMI, faceI, elem2YN, elemMYN, faceYN)
        ! this subroutine sets the wier and orifice initial condition.
        character(64) :: subroutine_name = 'QonlyElem_initial_condition_setup'


        real,      target, intent(in out)  :: elem2R(:,:),  elemMR(:,:)
        integer,           intent(in out)  :: faceI(:,:)
        real,      target, intent(in out)  :: faceR(:,:)
        integer,   target, intent(in)      :: elem2I(:,:),  elemMI(:,:)
        logical,   target, intent(in)      :: elem2YN(:,:), elemMYN(:,:)
        logical,   target, intent(in out)  :: faceYN(:,:)

        integer :: e2r_Volume_dummy, e2r_Velocity_dummy, eMr_Volume_dummy, eMr_Velocity_dummy

        real,      pointer  :: valueUp(:), valueDn(:)
        real,      pointer  :: weightUpQ(:), weightDnQ(:)
        real,      pointer  :: faceQ(:)
        logical,   pointer  :: facemask(:)


        integer :: mm
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% dummy pointers to call the wier/orifice step to find the flow through them.
        e2r_Volume_dummy = e2r_Temp(next_e2r_temparray)
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        e2r_Velocity_dummy = e2r_Temp(next_e2r_temparray)
        next_e2r_temparray = utility_advance_temp_array (next_e2r_temparray,e2r_n_temp)

        eMr_Volume_dummy = eMr_Temp(next_eMr_temparray)
        next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

        eMr_Velocity_dummy = eMr_Temp(next_eMr_temparray)
        next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

        valueUp => faceR(:,fr_Temp(next_fr_temparray))
        next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

        valueDn => faceR(:,fr_Temp(next_fr_temparray))
        next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

        weightUpQ => faceR(:,fr_Temp(next_fr_temparray))
        next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

        weightDnQ => faceR(:,fr_Temp(next_fr_temparray))
        next_fr_temparray = utility_advance_temp_array (next_fr_temparray, fr_n_temp)

        facemask   => faceYN(:,fYN_Temp(next_fYN_temparray))
        next_fYN_temparray = utility_advance_temp_array (next_fYN_temparray,fYN_n_temp)

        faceQ => faceR(:,fr_Flowrate)

        !% call weir and orifice step to initialize their geometry and flow
        call weir_step &
            (e2r_Volume_dummy, e2r_Velocity_dummy, eMr_Volume_dummy, eMr_Velocity_dummy, e2r_Volume, &
            e2r_Velocity, eMr_Volume, eMr_Velocity, elem2R, elemMR, &
            faceI, faceR, faceYN, elem2I, elemMI, elem2YN, elemMYN, 1.0)

        call orifice_step &
            (e2r_Volume_dummy, e2r_Velocity_dummy, eMr_Volume_dummy, eMr_Velocity_dummy, e2r_Volume, &
            e2r_Velocity, eMr_Volume, eMr_Velocity, elem2R, elemMR, &
            faceI, faceR, faceYN, elem2I, elemMI, elem2YN, elemMYN, 1.0)

        !% face reconstruction
        !% update the flow to their faces
        facemask = ( (faceI(:,fi_meta_etype_u) == eQonly) .or. &
            (faceI(:,fi_meta_etype_d) == eQonly) )

        weightUpQ = setting%Limiter%Timescale%Maximum
        weightDnQ = setting%Limiter%Timescale%Maximum

        where (facemask)
            weightUpQ = elem2R(faceI(:,fi_Melem_u),e2r_Timescale_Q_d)
            weightDnQ = elem2R(faceI(:,fi_Melem_d),e2r_Timescale_Q_u)
            valueUp  = elem2R(faceI(:,fi_Melem_u),e2r_Flowrate)
            valueDn  = elem2R(faceI(:,fi_Melem_d),e2r_Flowrate)
            !% linear interpolation
            faceQ = (weightUpQ * valueDn + weightDnQ * valueUp) /(weightUpQ + weightDnQ)
        endwhere

        valueUp    = nullvalueR
        valueDn    = nullvalueR
        weightUpQ  = nullvalueR
        weightDnQ  = nullvalueR

        nullify(valueUp, valueDn, weightUpQ, weightDnQ, facemask)

        next_e2r_temparray = next_e2r_temparray - 2
        next_eMr_temparray = next_eMr_temparray - 2
        next_fr_temparray  = next_fr_temparray  - 4
        next_fYN_temparray = next_fYN_temparray - 1

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

    end subroutine QonlyElem_initial_condition_setup
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine initial_solver_select &
        (elem2R, elemMR, elem2I, elemMI)

        character(64) :: subroutine_name = 'initial_solver_select'
        
        real,      intent(in)     :: elem2R(:,:), elemMR(:,:)
        integer,   intent(inout)  :: elem2I(:,:), elemMI(:,:)

        integer :: mm
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        where ( (elem2I(:,e2i_elem_type) == eChannel) .or. &
                (elem2I(:,e2i_elem_type) == ePipe   ) )
            ! elem2I(:,e2i_solver) = SVE
            ! elem2I(:,e2i_solver) = AC

            elem2I(:,e2i_solver) = AC
            ! elem2I(50:101,e2i_solver) = SVE
        endwhere
        
        where ( (elemMI(:,eMi_elem_type) == eJunctionChannel) .or. &
                (elemMI(:,eMi_elem_type) == eJunctionPipe   ) )
            ! elemMI(:,eMi_solver) = SVE
            elemMI(:,eMi_solver) = AC
        endwhere

        ! where ( (elem2I(:,e2i_elem_type) == ePipe)           .and. &
        !         (elem2R(:,e2r_Area)/elem2R(:,e2r_FullArea) .GE.  &
        !          setting%DefaultAC%Switch%Area) )
        !     elem2I(:,e2i_solver) = AC
        ! endwhere

        ! where ( (elem2I(:,eMi_elem_type) == eJunctionPipe)   .and. &
        !         (elemMR(:,eMr_Area)/elemMR(:,eMr_FullArea) .GE.  &
        !          setting%DefaultAC%Switch%Area) )
        !     elemMI(:,eMi_solver) = AC
        ! endwhere

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name
    end subroutine initial_solver_select
    !
    !==========================================================================
    !==========================================================================
    !
    !
    ! subroutine initial_condition_setupOLD &
    !    (elem2R, elemMR, elem2I, elemMI, elem2YN, elemMYN, &
    !     faceR, faceI, faceYN, bcdataDn, bcdataUp)
    !
    ! character(64) :: subroutine_name = 'initial_condition_setup'
    !
    ! real,      intent(in out)  :: elem2R(:,:), elemMR(:,:), faceR(:,:)
    !
    ! integer,   intent(in out)  :: elem2I(:,:), elemMI(:,:)
    !
    ! logical,   intent(in out)  :: elem2YN(:,:), elemMYN(:,:), faceYN(:,:)
    !
    ! integer,   intent(in out)      :: faceI(:,:)
    !
    ! type(bcType),  intent(in)  :: bcdataDn(:), bcdataUp(:)
    !
    ! real  ::  uniform_water_depth, uniform_bottom_roughness, uniform_flowrate
    !
    ! integer :: mm
    !!--------------------------------------------------------------------------
    ! if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name
    !
    ! select case (casename)
    !    case ('1link_network')
    !        uniform_water_depth = 0.5
    !        uniform_bottom_roughness = 0.03
    !        uniform_flowrate = 0.822225
    !
    !        print *, 'in ',subroutine_name,'---------------------------'
    !        print *, 'setting initial uniform water depth of ',uniform_water_depth
    !
    !!        call initial_condition_for_uniform_rectangular_channel &
    !!            (uniform_water_depth, uniform_bottom_roughness, &
    !!             elem2R, elemMR, elem2I, elemMI, elem2YN, elemMYN, &
    !!             faceI, bcdataDn, bcdataUp)
    !!
    !        call bc_nullify_ghost_elem (elem2R, bcdataDn)
    !        call bc_nullify_ghost_elem (elem2R, bcdataUp)
    !
    !        call element_geometry_update &
    !            (elem2R, elem2I, elem2YN, e2r_Volume, &
    !             elemMR, elemMI, elemMYN, eMr_Volume, &
    !             faceI, bcdataDn, bcdataUp)
    !
    !        ! HACK hard code setup of flowrate
    !
    !        elem2R(:,e2r_Flowrate) = uniform_flowrate
    !        elem2R(:,e2r_Velocity) = elem2R(:,e2r_Flowrate) / elem2R(:,e2r_Area)
    !        print *, 'setting initial flowrate of ',uniform_flowrate
    !
    !        call element_dynamics_update &
    !            (elem2R, elemMR, faceR, elem2I, elemMI, elem2YN, elemMYN, &
    !             bcdataDn, bcdataUp, e2r_Velocity, eMr_Velocity)
    !
    !        call face_update &
    !            (elem2R, elemMR, faceR, faceI, faceYN, &
    !             bcdataDn, bcdataUp, e2r_Velocity, eMr_Velocity)
    !
    !    case ('6link_1_line_network')
    !        print *, 'error initial condition not defined for ',casename,' in ',subroutine_name
    !        stop
    !
    !    case ('3link_Y_network')
    !        print *, 'error initial condition not defined for ',casename,' in ',subroutine_name
    !        stop
    !
    !    case ('6link_Y_network')
    !        print *, 'error initial condition not defined for ',casename,' in ',subroutine_name
    !        stop
    !
    !    case default
    !        print *
    !        print *, 'casename = ',casename
    !        print *, 'error: valid casename not selected in ',subroutine_name
    !        stop
    ! end select
    !
    ! ! remainder of geometry based on branches
    ! call junction_geometry_from_branches (elemMR, elemMI)
    !
    !
    ! if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    ! end subroutine initial_condition_setupOLD
    !!
    !==========================================================================
    !==========================================================================
    !
    ! subroutine initial_condition_for_uniform_rectangular_channelOLD &
    !    (uniform_water_depth, uniform_bottom_roughness, &
    !     elem2R, elemMR, elem2I, elemMI, elem2YN, elemMYN, &
    !     faceI, bcdataDn, bcdataUp)
    !!
    !! simple rectangular channels with initial uniform depth
    !!
    ! character(64) :: subroutine_name = 'initial_condition_for_uniform_rectangular_channel'
    !
    ! real,      intent(in out)  :: elem2R(:,:), elemMR(:,:)
    !
    ! integer,   intent(in out)  :: elem2I(:,:), elemMI(:,:)
    !
    ! logical,   intent(in out)  :: elem2YN(:,:), elemMYN(:,:)
    !
    ! integer,   intent(in)      :: faceI(:,:)
    !
    ! type(bcType),  target,  intent(in)  :: bcdataDn(:), bcdataUp(:)
    !
    ! real,  intent(in)  :: uniform_water_depth, uniform_bottom_roughness
    !
    ! integer,   pointer :: eID
    !
    ! integer :: ii
    !
    !!--------------------------------------------------------------------------
    ! if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name
    !
    ! ! 2-face elements - rectangular
    ! ! note that topwidth, zbottom, length were extracted from link data in network processing
    ! elem2R(:,e2r_Area)         = uniform_water_depth * elem2R(:,e2r_Topwidth)
    !
    ! elem2R(:,e2r_Volume)       = elem2R(:,e2r_Area) * elem2R(:,e2r_Length)
    !
    ! elem2R(:,e2r_Perimeter)    = elem2R(:,e2r_Topwidth) + twoR * uniform_water_depth
    !
    ! elem2R(:,e2r_Eta)          = uniform_water_depth + elem2R(:,e2r_Zbottom)
    !
    ! elem2R(:,e2r_Roughness)    = uniform_bottom_roughness
    !
    ! elem2R(:,e2r_Flowrate)     = zeroR
    ! elem2R(:,e2r_Velocity)     = zeroR
    !
    ! ! multi-face elements - rectangular
    ! ! note that topwidth, zbottom, length were extracted from link data in network processing
    ! elemMR(:,eMr_AreaAll)      = uniform_water_depth * elemMR(:,eMr_TopwidthAll)
    !
    ! elemMR(:,eMr_Eta)          = uniform_water_depth + elemMR(:,eMr_Zbottom)
    !
    ! elemMR(:,eMr_Perimeter)    = elemMR(:,eMr_Topwidth) + twoR * uniform_water_depth
    !
    ! elem2R(:,eMr_Roughness)    = uniform_bottom_roughness
    !
    ! elemMR(:,eMr_Flowrate)     = zeroR
    ! elemMR(:,eMr_Velocity)     = zeroR
    !
    !
    ! if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    ! end subroutine initial_condition_for_uniform_rectangular_channelOLD
    !!
    !==========================================================================
    ! END OF MODULE initial_condition
    !==========================================================================
end module initial_condition
