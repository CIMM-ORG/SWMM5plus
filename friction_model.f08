!==========================================================================
!
module friction_model
    !
    ! Compute the friction term used in the time advance
    ! Note that this is done separately so that it can be limited to prevent
    ! overestimation of explicit friction when flow is reversing.
    !
    use array_index
    use data_keys
    use globals
    use setting_definition

    implicit none

    private

    public  :: friction_on_element

    integer :: debuglevel = 0

contains
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine friction_on_element &
        (elemR, elemI, er_Friction, er_Velocity, er_Volume, &
        er_Roughness, er_HydRadius, ei_elem_type, ei_roughness_type,ThisElemType)
        !
        ! friction term in momentum on one type of element
        ! The general form is grav * volume * FrictionSlope
        ! This model uses Manning's n form for FrictionSlope
        !
        character(64) :: subroutine_name = 'friction_on_element'

        real(4),      target, intent(in out)  ::  elemR(:,:)

        integer,           intent(in)      ::  elemI(:,:)

        integer,           intent(in)      ::  er_Friction, er_Velocity, er_Volume
        integer,           intent(in)      ::  er_Roughness, er_HydRadius
        integer,           intent(in)      ::  ei_elem_type, ei_roughness_type
        integer,           intent(in)      ::  ThisElemType

        real(4),  pointer :: friction(:), velocity(:), volume(:), manningsn(:), rh(:)

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        friction =>  elemR(:,er_Friction)
        velocity =>  elemR(:,er_Velocity)
        volume   =>  elemR(:,er_Volume)
        manningsn=>  elemR(:,er_Roughness)
        rh       =>  elemR(:,er_HydRadius)

        where ( ( elemI(:,ei_elem_type) == ThisElemType ) .and. &
            ( elemI(:,ei_roughness_type) == eManningsN ) )
            friction = sign(grav * (manningsn**2) * (velocity**2) * volume / (rh**(4.0/3.0)), velocity)
        endwhere

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine friction_on_element
    !
    !==========================================================================
    ! END OF MODULE friction_model
    !==========================================================================
end module friction_model
