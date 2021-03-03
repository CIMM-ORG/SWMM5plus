! This is a HACK to put control on gates and validate AC solver
!
!==========================================================================
!
module control

    use adjustments
    use array_index
    use data_keys
    use globals
    use setting_definition
    use utility
    use type_definitions
    use xsect_tables

    implicit none

    private

    public :: control_allocate
    public :: control_evaluate
    public :: control_assign

    integer :: debuglevel = 0
contains
    !
    !==========================================================================
    !==========================================================================
    !
    Subroutine control_allocate (gateSetting, N_Gates, nChanges)
    !
    ! allocate controls 
    !
    character(64) :: subroutine_name = 'control_allocate'

        type(controlType), dimension(:),   allocatable,    intent(out) :: gateSetting

        integer,   intent(in)  :: N_Gates, nChanges

        integer    :: ii

        integer            :: allocation_status
        character(len=99)  :: emsg

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        allocate( gateSetting(N_Gates), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation (allocation_status, emsg)

        do ii = 1,N_Gates
            gateSetting(ii)%Idx         = ii
            gateSetting(ii)%LinkId      = nullvalueI
            gateSetting(ii)%ElemId      = nullvalueI 

            !% these arrays are for furute development
            allocate( gateSetting(ii)%TimeArray(nChanges), stat=allocation_status, errmsg=emsg)
            call utility_check_allocation (allocation_status, emsg)

            allocate( gateSetting(ii)%HeightArray(nChanges), stat=allocation_status, errmsg=emsg)
            call utility_check_allocation (allocation_status, emsg)

            allocate( gateSetting(ii)%AreaArray(nChanges), stat=allocation_status, errmsg=emsg)
            call utility_check_allocation (allocation_status, emsg)
            gateSetting(ii)%HeightStart         = nullvalueR
            gateSetting(ii)%HeightNow           = nullvalueR
            gateSetting(ii)%AreaNow             = nullvalueR
            gateSetting(ii)%AreaPrior           = nullvalueR
            !% HACK code
            gateSetting(ii)%GateTimeChange1     = nullvalueR
            gateSetting(ii)%GateTimeChange2     = nullvalueR
            gateSetting(ii)%GateHeightChange1   = nullvalueR
            gateSetting(ii)%GateHeightChange2   = nullvalueR
            gateSetting(ii)%HeightMinimum       = nullvalueR
            gateSetting(ii)%GateSpeed           = nullvalueR
            gateSetting(ii)%CanMove             = nullvalueL
            gateSetting(ii)%MovedThisStep       = nullvalueL
        end do
     
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name    
    end Subroutine control_allocate
    !
    !==========================================================================
    !==========================================================================
    !
    Subroutine control_assign (elem2I, gateSetting, N_Gates)
    !
    ! assign element to control
    !
    character(64) :: subroutine_name = 'control_evaluate'

        integer,                   intent(in)     :: elem2I(:,:)
        type(controlType), target, intent(inout)  :: gateSetting(:)

        integer,                    intent(in)    :: N_Gates
        integer, pointer                          :: eID

        integer    :: ii, mm
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

        !%  cycle through gates
        do ii=1,N_Gates
            eID => gateSetting(ii)%ElemId
            do mm=1,N_elem2
                if (elem2I(mm,e2i_link_ID) == gateSetting(ii)%LinkId) then
                    eID = elem2I(mm,e2i_idx)
                endif
            enddo
        enddo

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name    
    end Subroutine control_assign
    !
    !==========================================================================
    !==========================================================================
    !
    Subroutine control_evaluate (elem2I, elem2R, gateSetting, N_Gates, StepTime)
    !
    ! HACK: this code is adapted from pipeAC. the sole purpose of this code
    ! is to test and compare AC-solver by utilizing trajkovic test cases 
    !
    character(64) :: subroutine_name = 'control_evaluate'

        integer,                   intent(in)     :: elem2I(:,:)
        real,                      intent(inout)  :: elem2R(:,:)
        type(controlType), target, intent(inout)  :: gateSetting(:)

        integer,             intent(in)    :: N_Gates
        real,                intent(in)    :: StepTime
        integer, pointer                   :: eID
        real                               :: gateInterval, gateMove, gateChange
        real                               :: newFullDepth, YoverYfull, newFullArea

        integer    :: ii, mm
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name

        !%  cycle through gates
        do ii=1,N_Gates

            eID => gateSetting(ii)%ElemId
            !%  for fixed gates
            if (.not. gateSetting(ii)%CanMove) then

                gateSetting(ii)%HeightNow = gateSetting(ii)%HeightStart 
            !%  for moving gates
            elseif (gateSetting(ii)%CanMove) then

                !%  if still in the first height regime
                if (StepTime < gateSetting(ii)%GateTimeChange1) then
                    gateSetting(ii)%HeightNow = gateSetting(ii)%HeightStart 

                !%  if in the second height regime 
                elseif ( (StepTime .GE. gateSetting(ii)%GateTimeChange1) .and. &
                         (StepTime .LT. gateSetting(ii)%GateTimeChange2)  ) then

                    gateInterval = StepTime - gateSetting(ii)%GateTimeChange1 
                    gateMove     = gateInterval * gateSetting(ii)%GateSpeed
                    gateChange   = gateSetting(ii)%HeightStart - gateSetting(ii)%GateHeightChange1

                    if (gateChange > zeroR) then
                        !%  gate is moving downward
                        gateSetting(ii)%HeightNow = gateSetting(ii)%HeightStart - gateMove

                        if (gateSetting(ii)%HeightNow < gateSetting(ii)%GateHeightChange1) then
                            gateSetting(ii)%HeightNow = gateSetting(ii)%GateHeightChange1
                        endif

                        gateSetting(ii)%MovedThisStep = .true.

                        ! if (gateSetting(ii)%MovedThisStep) then
                        !     print*, 'gate is moving downward ', gateSetting(ii)%HeightNow
                        ! endif

                    else
                        gateSetting(ii)%HeightNow = gateSetting(ii)%HeightStart + gateMove

                        if (gateSetting(ii)%HeightNow > gateSetting(ii)%GateHeightChange1) then
                            gateSetting(ii)%HeightNow = gateSetting(ii)%GateHeightChange1
                        endif
                    
                        gateSetting(ii)%MovedThisStep = .true.

                        ! if (gateSetting(ii)%MovedThisStep) then
                        !     print*, 'gate is moving upward ', gateSetting(ii)%HeightNow, ' at time ', StepTime
                        ! endif

                    endif

                !%  if in the third height regime
                elseif ( (StepTime .GE. gateSetting(ii)%GateTimeChange2                   ) .and. &
                         (gateSetting(ii)%HeightNow .NE. gateSetting(ii)%GateHeightChange2)  ) then

                    gateInterval = StepTime - gateSetting(ii)%GateTimeChange2 
                    gateMove     = gateInterval * gateSetting(ii)%GateSpeed
                    gateChange   = gateSetting(ii)%HeightStart - gateSetting(ii)%GateHeightChange2

                    if (gateChange > zeroR) then
                        !%  gate is moving downward
                        gateSetting(ii)%HeightNow = gateSetting(ii)%GateHeightChange1 - gateMove

                        if (gateSetting(ii)%HeightNow < gateSetting(ii)%GateHeightChange2) then
                            gateSetting(ii)%HeightNow = gateSetting(ii)%GateHeightChange2
                        endif

                        gateSetting(ii)%MovedThisStep = .true.

                        ! if (gateSetting(ii)%MovedThisStep) then
                        !     print*, 'gate is moving downward ', gateSetting(ii)%HeightNow, ' at time ', StepTime
                        ! endif

                    else
                        gateSetting(ii)%HeightNow = gateSetting(ii)%GateHeightChange1 + gateMove

                        if (gateSetting(ii)%HeightNow > gateSetting(ii)%GateHeightChange2) then
                            gateSetting(ii)%HeightNow = gateSetting(ii)%GateHeightChange2
                        endif

                        gateSetting(ii)%MovedThisStep = .true.

                        ! if (gateSetting(ii)%MovedThisStep) then
                        !     print*, 'gate is moving upward ', gateSetting(ii)%HeightNow, ' at time ', StepTime
                        ! endif

                    endif
                endif

            else
                print*, 'error -- unexpected else'
                print*, gateSetting(ii)%CanMove
                stop

            endif

            gateSetting(ii)%HeightNow = max(gateSetting(ii)%HeightNow, gateSetting(ii)%HeightMinimum)

            !%  calculate new full depth and area of the element
            !%  HACK: only works for orifice elements
            if ( (elem2I(eID,e2i_elem_type) == eOrifice) .and. &
                 (elem2I(eID,e2i_geometry) == eCircular)  ) then

                YoverYfull                = gateSetting(ii)%HeightNow / elem2r(eID,e2r_FullDepth)
                elem2R(eID,e2r_FullArea)  = elem2R(eID,e2r_FullArea) * table_lookup(YoverYfull, ACirc, NACirc)
                elem2R(eID,e2r_FullDepth) = gateSetting(ii)%HeightNow
                elem2R(eID,e2r_Zcrown)    = elem2R(eID,e2r_Zbottom) + elem2R(eID,e2r_FullDepth)

            elseif ( (elem2I(eID,e2i_elem_type) == eOrifice)   .and. &
                     (elem2I(eID,e2i_geometry) == eRectangular) ) then

                elem2R(eID,e2r_FullDepth) = gateSetting(ii)%HeightNow
                elem2R(eID,e2r_FullArea)  = elem2R(eID,e2r_FullDepth) * elem2R(eID,e2r_BreadthScale)
                elem2R(eID,e2r_Zcrown)    = elem2R(eID,e2r_Zbottom) + elem2R(eID,e2r_FullDepth)
            endif
        enddo


        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name    
    end Subroutine control_evaluate
end module control
