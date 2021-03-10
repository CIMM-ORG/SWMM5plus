! This module calculates the depth in storage nodes
!
!==========================================================================
!
module storage


    use adjustments
    use array_index
    use bc
    use data_keys
    use diagnostic
    use face_values
    use globals
    use setting_definition
    use utility

    implicit none

    private

    public :: storage_adjacent_element_average
    public :: storage_initialize_depth_volume
    public :: storage_step

    integer :: debuglevel = 0

contains
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine storage_adjacent_element_average &
        (elem2R, elemMR, elemMI, faceI, e2r_data, eMr_out)
        !
        ! this computes the average of values for all the elements upstream and
        ! downstream. Note that this should ONLY be use
        ! in setup routines (i.e. when initializing storage). This violates
        ! the "no-neighbor" rule and is time-consuming because it requires double
        ! mapping
        !
        character(64) :: subroutine_name = 'storage_adjacent_element_average'

        real(8),      target,     intent(in out)  :: elemMR(:,:)
        real(8),                  intent(in)      :: elem2R(:,:)
        integer,               intent(in)      :: elemMI(:,:), faceI(:,:)
        integer,               intent(in)      :: e2r_data, eMr_out

        real(8),      pointer :: Uvalue(:), Dvalue(:)
        integer :: eMr_tUp, eMr_tDn

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        eMr_tUp = eMr_Temp(next_eMr_temparray)
        Uvalue  => elemMR(:,eMr_tUp)
        next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

        eMr_tDn = eMr_Temp(next_eMr_temparray)
        Dvalue  => elemMR(:,eMr_tDn)
        next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

        Uvalue = zeroR
        call storage_summation_from_adjacent_elements_one_direction &
            (eMr_tUp, elem2R, elemMR, elemMI, faceI, &
            upstream_face_per_elemM, eMi_nfaces_u, eMi_MfaceUp, fi_Melem_u, e2r_data)

        Dvalue = zeroR
        call storage_summation_from_adjacent_elements_one_direction &
            (eMr_tDn, elem2R, elemMR, elemMI, faceI, &
            dnstream_face_per_elemM, eMi_nfaces_d, eMi_MfaceDn, fi_Melem_d, e2r_data)

        where (elemMI(:,eMi_elem_type) == eStorage)
            elemMR(:,eMr_out) = (Uvalue + Dvalue) / real( elemMI(:,eMi_nfaces_u) + elemMI(:,eMi_nfaces_d) )
        endwhere

        Dvalue = nullvalueR
        Uvalue = nullvalueR
        nullify(Dvalue,Uvalue)
        next_eMr_temparray = next_eMr_temparray-2

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine storage_adjacent_element_average
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine storage_initialize_depth_volume (elemMR, elemMI)

        character(64) :: subroutine_name = 'storage_initialize_depth_volume'
        !
        ! computes the volume and surface area in a storage unit by functional or tabular relationship
        ! this subroutine is only used in initialization of a storage unit
        !
        real(8),      target, intent(in out)  :: elemMR(:,:)
        integer,   target, intent(in)      :: elemMI(:,:)

        integer,   pointer :: cindx
        real(8),      pointer :: volume, depth, fullvolume, fulldepth
        real(8),      pointer :: surfarea, aConst, aCoeff, aExpon

        integer :: mm, ii

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        do ii=1,N_elemM
            cindx       => elemMI(ii,eMi_curve_type)
            depth       => elemMR(ii,eMr_HydDepth)
            fulldepth   => elemMR(ii,eMr_FullDepth)
            volume      => elemMR(ii,eMr_Volume)
            fullvolume  => elemMR(ii,eMr_FullVolume)
            surfarea    => elemMR(ii,eMr_SurfArea)
            ! functional pointers
            aConst   => elemMR(ii,eMr_StorageConstant)
            aCoeff   => elemMR(ii,eMr_StorageCoeff)
            aExpon   => elemMR(ii,eMr_StorageExponent)

            ! calculate storage depth and full volume
            if (elemMI(ii,eMi_elem_type) == eStorage) then
                !% depth will always be >= zeroR
                depth = max((elemMR(ii,eMr_Eta) - elemMR(ii,eMr_Zbottom)), zeroR)
                !% if calculated depth is greater than fulldepth, new depth is full depth
                depth = min(depth,fulldepth)
            endif

            select case (cindx)
              case (1)
                !% case where the relationship between depth and surface area is Functional
                if (elemMI(ii,eMi_elem_type) == eStorage) then

                    volume = aConst * depth + aCoeff / ((aExpon + oneR) &
                        * depth ** (aExpon + oneR))
                    surfarea   = aConst + aCoeff * depth ** aExpon

                    fullvolume = aConst * fulldepth + aCoeff / ((aExpon + oneR) &
                        * fulldepth ** (aExpon + oneR))
                endif

              case(2)
                !% case where the relationship between depth and surface area is Tabular
                if (elemMI(ii,eMi_elem_type) == eStorage) then

                    print*, 'error: tabular storage curve type in not handeled at this moment'
                    stop
                endif
            end select
        enddo
        ! print*, '-------------storage initialization----------'
        ! print*, depth,' <= initial depth'
        ! print*, volume, ' <= initial volume'
        ! print*, fullvolume, ' <= full volume'
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine storage_initialize_depth_volume
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine storage_step &
        (eMr_Volume_old, eMr_Velocity_old, eMr_Volume_new, eMr_Velocity_new,  &
        elemMR, faceR, elemMI, elemMYN, thiscoef)

        character(64) :: subroutine_name = 'storage_step'
        !
        ! This subroutine takes a time step for storage units
        !
        ! indexes for old/new volume and velocity storage
        integer,   intent(in) :: eMr_Volume_old, eMr_Volume_new
        integer,   intent(in) :: eMr_Velocity_old, eMr_Velocity_new

        real(8),      target, intent(in out)  :: elemMR(:,:)
        real(8),      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elemMI(:,:)
        logical,           intent(in out)  :: elemMYN(:,:)
        real(8),              intent(in)      :: thiscoef

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        call storge_net_inflow_outflow &
            (eMr_Volume_old, eMr_Volume_new, elemMR, faceR, elemMI, elemMYN, &
            thiscoef)
        call storage_get_depth (eMr_Volume_new, elemMR, faceR, elemMI, elemMYN)

        !% setting dummy zero velocity for storage elements
        where (elemMI(:,eMi_elem_type) == eStorage)
            elemMR(:,eMr_Velocity_new) = 1.0E-7
        endwhere

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine storage_step
    !
    !==========================================================================
    ! PRIVATE BELOW
    !==========================================================================
    !
    subroutine storage_summation_from_adjacent_elements_one_direction &
        (eMr_sumvalue, elem2R, elemMR, elemMI, faceI, &
        dir_face_per_elemM, eMi_nfaces_dir, eMi_MfaceDir, fi_Melem_dir, e2r_data)

        character(64) :: subroutine_name = 'storage_summation_from_adjacent_elements_one_direction'
        !
        ! computes the sum of all the elements adjacent to a junctio in either
        ! the upstream or downstream direction
        !
        ! THIS SHOULD ONLY BE USED IN SETUP AND INITIAL CONDITION ROUTINES
        !
        integer,           intent(in)      :: eMr_sumvalue
        real(8),      target, intent(in out)  :: elemMR(:,:)
        real(8),              intent(in)      :: elem2R(:,:)
        integer,   target, intent(in)      :: elemMI(:,:), faceI(:,:)

        integer,   intent(in)  :: dir_face_per_elemM, eMi_nfaces_dir, eMi_MfaceDir(:)
        integer,   intent(in)  :: fi_Melem_dir, e2r_data

        integer,   pointer :: tface, telem
        real(8)   :: thisvalue(dir_face_per_elemM)

        integer :: mm, ii

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


        do ii=1,N_elemM
            thisvalue = zeroR
            do mm=1,dir_face_per_elemM
                if ((elemMI(ii,eMi_nfaces_dir) >= mm) .and. &
                    (elemMI(ii,eMi_elem_type) == eStorage)) then
                    !% the face on a branch in the direction specified by dir
                    tface => elemMI(ii,eMi_MfaceDir(mm))
                    !% the element upstream of the face
                    telem => faceI(tface,fi_Melem_dir)
                    !% the value at the element
                    thisvalue(mm) = elem2R(telem,e2r_data)
                endif
            enddo
            elemMR(ii,eMr_sumvalue) = sum(thisvalue)
        enddo

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine storage_summation_from_adjacent_elements_one_direction
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine storge_net_inflow_outflow &
        (eMr_Volume_old, eMr_Volume_new, elemMR, faceR, elemMI, elemMYN, &
        thiscoef)

        character(64) :: subroutine_name = 'storge_net_inflow_outflow'
        !
        ! This subroutine calculates the new volume in a storage element
        ! from the continuity equation. Called in storage step
        !
        ! indexes for old/new volume and velocity storage
        integer,   intent(in) :: eMr_Volume_old, eMr_Volume_new

        real(8),      target, intent(in out)  :: elemMR(:,:)
        real(8),      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elemMI(:,:)
        logical,           intent(in out)  :: elemMYN(:,:)
        real(8),              intent(in)      :: thiscoef


        real(8),      pointer     :: fQ(:), dV(:), volumeMold(:), volumeMnew(:)
        integer,   pointer     :: iup(:), idn(:)

        integer        :: mm

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !%  pointers for volume storage (updating)
        volumeMold   => elemMR(:,eMr_Volume_old)
        volumeMnew   => elemMR(:,eMr_Volume_new)

        !%  pointers for convenience in notation
        fQ => faceR(:,fr_Flowrate)

        !%  temporary space for storage
        dV => elemMR(:,eMr_Temp(next_eMr_temparray))
        next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)
        !%  zero temporary arrays
        dV = zeroR

        !%  Storage (upstream faces)
        do mm=1,upstream_face_per_elemM
            iup   => elemMI(:,eMi_MfaceUp(mm))
            where ( (elemMI(:,eMi_elem_type) == eStorage).and. &
                (elemMI(:,eMi_nfaces_u) >= mm) )
                dV = dV + dt * fQ(iup)
            endwhere
        enddo

        !%  Storage (downstream faces)
        do mm=1,dnstream_face_per_elemM
            idn   => elemMI(:,eMi_MfaceDn(mm))
            where ( (elemMI(:,eMi_elem_type) == eStorage) .and. &
                (elemMI(:,eMi_nfaces_d) >= mm) )
                dV = dV - dt * fQ(idn)
            endwhere
        enddo

        where (elemMI(:,eMi_elem_type) == eStorage)
            volumeMnew   = min((volumeMold  + thiscoef * dV), elemMR(:,eMr_FullVolume))
        endwhere

        !%  remove negative volumes to prevent problems in velocity computation
        call adjust_negative_volume_reset (volumeMnew)

        !% release temporary arrays
        dV = nullvalueR
        nullify(dV)
        next_eMr_temparray = next_eMr_temparray - 1

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine storge_net_inflow_outflow
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine storage_get_depth &
        (eMr_Volume_new, elemMR, faceR, elemMI, elemMYN)

        character(64) :: subroutine_name = 'storage_get_depth'
        !
        ! This subroutine calculates the new depth in the storage unit
        !
        ! indexes for old/new volume and velocity storage
        integer,   intent(in) :: eMr_Volume_new

        real(8),      target, intent(in out)  :: elemMR(:,:)
        real(8),      target, intent(in)      :: faceR(:,:)
        integer,   target, intent(in)      :: elemMI(:,:)
        logical,           intent(in out)  :: elemMYN(:,:)

        real(8),      pointer :: volumeMnew, depth, eta, fulldepth, fullvolume, zbottom
        real(8),      pointer :: aConst, aCoeff, aExpon
        integer,   pointer :: cindx

        integer        :: ii, mm

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        do ii=1,N_elemM
            cindx       => elemMI(ii,eMi_curve_type)
            volumeMnew  => elemMR(ii,eMr_Volume_new)
            depth       => elemMR(ii,eMr_HydDepth)
            eta         => elemMR(ii,eMr_Eta)
            fulldepth   => elemMR(ii,eMr_FullDepth)
            fullvolume  => elemMR(ii,eMr_FullVolume)
            zbottom     => elemMR(ii,eMr_Zbottom)
            ! functional pointers
            aConst   => elemMR(ii,eMr_StorageConstant)
            aCoeff   => elemMR(ii,eMr_StorageCoeff)
            aExpon   => elemMR(ii,eMr_StorageExponent)

            if ((elemMI(ii,eMi_elem_type) == eStorage) .and. (volumeMnew .GE. fullvolume)) then

                depth = fulldepth
                eta = depth + zbottom
                !% need to calculate the overflow here

            elseif (elemMI(ii,eMi_elem_type) == eStorage)  then
                select case (cindx)
                  case (1)
                    !% case where the relationship between depth and surface area is Functional
                    if (aExpon == zeroR) then
                        depth = volumeMnew / (aConst + aCoeff)

                    elseif (aConst == zeroR) then

                        depth = (volumeMnew / aCoeff * (oneR / (aExpon + oneR))) ** (oneR / (aExpon + oneR))
                    else
                        print*, 'finding root using Newton Raphson method still in dev'
                        stop
                    endif

                  case(2)
                    !% case where the relationship between depth and surface area is Tabular
                    if (elemMI(ii,eMi_elem_type) == eStorage) then
                        print*, 'error: tabular storage curve type in not handeled at this moment'
                        stop
                    endif
                end select
                eta = depth + zbottom
            endif
        end do

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine storage_get_depth
    !
    !==========================================================================
    !==========================================================================
    !
end module storage
