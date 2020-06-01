!==========================================================================
!
module output
        !
        use array_index
        use bc
        use data_keys
        use globals
        use setting_definition
        use utility


        implicit none

        private

        public  :: output_threaded_by_link_initialize
        public  :: output_all_threaded_data_by_link
        public  :: output_one_threaded_data_by_link

        integer :: debuglevel = 0

contains
        !
        !==========================================================================
        !==========================================================================
        !
        subroutine output_threaded_by_link_initialize (threadedfile)

                character(64) :: subroutine_name = 'output_threaded_by_link_initialize'

                type(threadedfileType), dimension(:), allocatable, intent(out) :: threadedfile

                logical    :: usethisdata(5)

                integer    :: ii, mm, nfile

                integer            :: allocation_status
                character(len=99)  :: emsg

                !--------------------------------------------------------------------------
                if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

                if (.not. setting%OutputThreadedLink%UseThisOutput) return

                if (setting%OutputThreadedLink%SuppressAllFiles) return

                usethisdata = .false.

                if (setting%OutputThreadedLink%area)       usethisdata(1) = .true.
                if (setting%OutputThreadedLink%flowrate)   usethisdata(2) = .true.
                if (setting%OutputThreadedLink%velocity)   usethisdata(3) = .true.
                if (setting%OutputThreadedLink%eta)        usethisdata(4) = .true.
                if (setting%OutputThreadedLink%depth)      usethisdata(5) = .true.

                nfile = count(usethisdata)

                if (nfile == 0) return

                allocate( threadedfile(nfile), stat=allocation_status, errmsg=emsg)
                call utility_check_allocation (allocation_status, emsg)

                mm = 1
                do ii=1,5
                        if (usethisdata(ii)) then
                                select case (ii)
                                    case (1)
                                        threadedfile(mm)%DataName = 'area'
                                        threadedfile(mm)%FileInfo%FileName = trim(setting%OutputThreadedLink%FileName)//'_area'
                                        mm=mm+1
                                    case (2)
                                        threadedfile(mm)%DataName = 'flowrate'
                                        threadedfile(mm)%FileInfo%FileName = trim(setting%OutputThreadedLink%FileName)//'_flowrate'
                                        mm=mm+1
                                    case (3)
                                        threadedfile(mm)%DataName = 'velocity'
                                        threadedfile(mm)%FileInfo%FileName = trim(setting%OutputThreadedLink%FileName)//'_velocity'
                                        mm=mm+1
                                    case (4)
                                        threadedfile(mm)%DataName = 'eta'
                                        threadedfile(mm)%FileInfo%FileName = trim(setting%OutputThreadedLink%FileName)//'_eta'
                                        mm=mm+1
                                    case (5)
                                        threadedfile(mm)%DataName = 'depth'
                                        threadedfile(mm)%FileInfo%FileName = trim(setting%OutputThreadedLink%FileName)//'_depth_'
                                        mm=mm+1
                                    case default
                                        print *, 'error: unexpected value'
                                        stop
                                end select
                        endif
                enddo

                do mm=1,size(threadedfile)
                        threadedfile(mm)%FileInfo%FolderName = trim(setting%OutputThreadedLink%FolderName)
                        threadedfile(mm)%FileInfo%FolderPath = trim(setting%OutputThreadedLink%FolderPath)
                        threadedfile(mm)%FileInfo%FileStatus = 'new'
                        threadedfile(mm)%FileInfo%IsOpen    = .false.

                        threadedfile(mm)%FileInfo%WriteName = trim(threadedfile(mm)%FileInfo%FolderPath)// &
                                trim(threadedfile(mm)%FileInfo%FolderName)//'/'// &
                                trim(threadedfile(mm)%FileInfo%FileName)// &
                                trim(setting%Time%DateTimeStamp) //&
                                '.txt'

                        call output_singlethreadedfile_open (threadedfile(mm))
                enddo

                if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
        end subroutine output_threaded_by_link_initialize
        !
        !==========================================================================
        !==========================================================================
        !
        subroutine output_all_threaded_data_by_link &
                (threadedfile, elem2R, elem2I, elemMR, elemMI, faceR, faceI, linkI,&
                bcdataUp, bcdataDn, thisstep)

                character(64) :: subroutine_name = 'output_all_threaded_data_by_link'

                type(threadedfileType), intent(in) :: threadedfile(:)

                real,      target, intent(in)  :: elem2R(:,:), elemMR(:,:), faceR(:,:)
                integer,   target, intent(in)  :: elem2I(:,:), elemMI(:,:), faceI(:,:)
                integer,   target, intent(in)  :: linkI(:,:)
                type(bcType),      intent(in)  :: bcdataUp(:), bcdataDn(:)
                integer,           intent(in)  :: thisstep

                logical    :: usethisdata(5)

                character(len=32)  :: outdataName
                integer            :: thisUnit, ii, mm, nfile
                !--------------------------------------------------------------------------
                if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

                if (setting%OutputThreadedLink%SuppressAllFiles) return

                if (setting%OutputThreadedLink%DisplayInterval <= 0) return

                if (mod(thisstep,setting%OutputThreadedLink%DisplayInterval) /= 0) return

                usethisdata = .false.



                if (setting%OutputThreadedLink%area)       usethisdata(1) = .true.
                if (setting%OutputThreadedLink%flowrate)   usethisdata(2) = .true.
                if (setting%OutputThreadedLink%velocity)   usethisdata(3) = .true.
                if (setting%OutputThreadedLink%eta)        usethisdata(4) = .true.
                if (setting%OutputThreadedLink%depth)      usethisdata(5) = .true.

                nfile = count(usethisdata)

                if (nfile == 0) return

                mm=1
                do ii=1,5
                        if (usethisdata(ii)) then
                                select case (ii)
                                    case (1)
                                        outdataName = 'area'
                                        thisUnit = threadedfile(mm)%FileInfo%Unitnumber
                                        mm=mm+1
                                    case (2)
                                        outdataName = 'flowrate'
                                        thisUnit = threadedfile(mm)%FileInfo%Unitnumber
                                        mm=mm+1
                                    case (3)
                                        outdataName = 'velocity'
                                        thisUnit = threadedfile(mm)%FileInfo%Unitnumber
                                        mm=mm+1
                                    case (4)
                                        outdataName = 'eta'
                                        thisUnit = threadedfile(mm)%FileInfo%Unitnumber
                                        mm=mm+1
                                    case (5)
                                        outdataName = 'depth'
                                        thisUnit = threadedfile(mm)%FileInfo%Unitnumber
                                        mm=mm+1
                                    case default
                                        print *, 'error: unexpected case value in ',trim(subroutine_name)
                                        stop
                                end select

                                call output_one_threaded_data_by_link &
                                        (outdataName, thisUnit, thisstep, &
                                        elem2R, elem2I, elemMR, elemMI, faceR, faceI, linkI, bcdataUp, bcdataDn, ii)

                                ! if (ii==2) then
                                !    print *, trim(subroutine_name)
                                !    stop
                                ! endif

                        end if
                end do

                if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
        end subroutine output_all_threaded_data_by_link
        !
        !==========================================================================
        !==========================================================================
        !
        subroutine output_one_threaded_data_by_link &
                (outdataName, thisUnit, thisstep, &
                elem2R, elem2I, elemMR, elemMI, faceR, faceI, linkI, bcdataUp, bcdataDn, itemp)
                !
                ! provides an output that threads the downstream BC, element or junction
                ! along with faces and links into a single array for visualization
                !
                character(64) :: subroutine_name = 'output_one_threaded_data_by_link'

                character(len=*), intent(in)  :: outdataName
                integer,           intent(in)  :: thisUnit

                real,      target, intent(in)  :: elem2R(:,:), elemMR(:,:), faceR(:,:)
                integer,   target, intent(in)  :: elem2I(:,:), elemMI(:,:), faceI(:,:)
                integer,   target, intent(in)  :: linkI(:,:)
                type(bcType),      intent(in)  :: bcdataUp(:), bcdataDn(:)
                integer,           intent(in)  :: thisstep, itemp

                real, dimension(:), allocatable, save  :: thisdata, thisx

                integer :: ElemCol, FaceCol_u, FaceCol_d, JunctionCol
                integer :: BranchColUp(upstream_face_per_elemM), BranchColDn(dnstream_face_per_elemM)

                integer, pointer :: fdn, ftypDn, edn, bcidxDn, bdn, Lidx
                integer, pointer :: fup, ftypUp, eup, bcidxUp, bup

                integer :: maxelements_in_link, nelem, eStart, eEnd, eLast

                integer :: en, ii, mm, ndata

                integer            :: allocation_status
                character(len=99)  :: emsg

                !--------------------------------------------------------------------------
                if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

                if (setting%OutputThreadedLink%SuppressAllFiles) return

                !% allocate space for combining data on elements, face up, face dn, and the upstream
                !% and downstream adjacent element or BC
                if (.not. allocated(thisdata)) then
                        maxelements_in_link = maxval(linkI(:,li_N_element),1,(linkI(:,li_N_element) /= nullvalueI))

                        allocate( thisdata(3*maxelements_in_link + 7), stat=allocation_status, errmsg=emsg)
                        call utility_check_allocation (allocation_status, emsg)

                        allocate( thisX(3*maxelements_in_link + 7), stat=allocation_status, errmsg=emsg)
                        call utility_check_allocation (allocation_status, emsg)
                endif
                thisdata = nullvalueR
                thisX    = nullvalueR

                select case (trim(outdataName))
                    case ('area')
                        ElemCol     = e2r_Area
                        FaceCol_u   = fr_Area_u
                        FaceCol_d   = fr_Area_d
                        JunctionCol = eMr_Area
                        BranchColUp = eMr_AreaUp
                        BranchColDn = eMr_AreaDn
                    case ('flowrate')
                        ElemCol     = e2r_Flowrate
                        FaceCol_u   = fr_Flowrate
                        FaceCol_d   = fr_Flowrate
                        JunctionCol = eMr_Flowrate
                        BranchColUp = eMr_FlowrateUp
                        BranchColDn = eMr_FlowrateDn
                    case ('velocity')
                        ElemCol     = e2r_Velocity
                        FaceCol_u   = fr_Velocity_u
                        FaceCol_d   = fr_Velocity_d
                        JunctionCol = eMr_Velocity
                        BranchColUp = eMr_VelocityUp
                        BranchColDn = eMr_VelocityDn
                    case ('eta')
                        ElemCol     = e2r_Eta
                        FaceCol_u   = fr_Eta_u
                        FaceCol_d   = fr_Eta_d
                        JunctionCol = eMr_Eta
                        BranchColUp = eMr_EtaUp
                        BranchColDn = eMr_EtaDn
                    case ('depth')
                        ElemCol     = e2r_HydDepth
                        FaceCol_u   = fr_HydDepth_u
                        FaceCol_d   = fr_HydDepth_d
                        JunctionCol = eMr_HydDepth
                        BranchColUp = eMr_HydDepthUp
                        BranchColDn = eMr_HydDepthDn
                    case default
                        print *, trim(subroutine_name)
                        print *, outdataName
                        print *, 'error: unexpected outdataName of ',outdataName,' in ',trim(subroutine_name)
                        stop
                end select

                do ii=1,N_link

                        !%  this link
                        Lidx => linkI(ii,li_idx)
                        !%  start and end elements on this link
                        eStart = minloc(elem2I(:,e2i_link_Pos),1,(  (elem2I(:,e2i_link_ID) == Lidx) &
                                .and. (elem2I(:,e2i_link_Pos) /= nullvalueI) ) )
                        eEnd   = maxloc(elem2I(:,e2i_link_Pos),1,(  (elem2I(:,e2i_link_ID) == Lidx) &
                                .and. (elem2I(:,e2i_link_Pos) /= nullvalueI) ) )
                        nelem = eEnd - eStart + 1

                        !%  the downstream element beyond this link
                        fdn     => elem2I(eStart,e2i_Mface_d)
                        ftypDn  => faceI(fdn,fi_type)
                        edn     => faceI(fdn,fi_Melem_d)
                        bcidxDn => faceI(fdn,fi_BC_ID)

                        thisdata = nullvalueR
                        thisX    = nullvalueR
                        ! print*, ii
                        ! get the data of the element (or BC) downstream of link
                        select case (ftypDn)
                            case (fChannel)
                                ! downstream last element from previous link - upstream face
                                thisdata(1) = faceR(elem2I(edn,e2i_Mface_d),FaceCol_u)
                                thisX(1)    = -elem2R(edn,e2r_Length)
                                ! downstream last element from previous linke - center
                                thisdata(2) =  elem2R(edn,elemCol)
                                thisX(2)    = -onehalfR * elem2R(edn,e2r_Length)
                                ! face downstream at link connection (node)
                                thisdata(3) = faceR(fdn,FaceCol_d)
                                thisX(3)    = -oneeighthR * elem2R(edn,e2r_Length)
                                ! face upstream at link connection (node)
                                thisdata(4) = faceR(fdn,FaceCol_u)
                                thisX(4)    = zeroR
                            case (fPipe)
                                print *, 'error: pipe output not handled'
                                stop
                            case (fWeir)
                                ! downstream last element from previous link - upstream face
                                thisdata(1) = faceR(elem2I(edn,e2i_Mface_d),FaceCol_u)
                                thisX(1)    = -elem2R(edn,e2r_Length)
                                ! downstream last element from previous linke - center
                                thisdata(2) =  elem2R(edn,elemCol)
                                thisX(2)    = -onehalfR * elem2R(edn,e2r_Length)
                                ! face downstream at link connection (node)
                                thisdata(3) = faceR(fdn,FaceCol_d)
                                thisX(3)    = -oneeighthR * elem2R(edn,e2r_Length)
                                ! face upstream at link connection (node)
                                thisdata(4) = faceR(fdn,FaceCol_u)
                                thisX(4)    = zeroR
                            case (fOrifice)
                                ! downstream last element from previous link - upstream face
                                thisdata(1) = faceR(elem2I(edn,e2i_Mface_d),FaceCol_u)
                                thisX(1)    = -elem2R(edn,e2r_Length)
                                ! downstream last element from previous linke - center
                                thisdata(2) =  elem2R(edn,elemCol)
                                thisX(2)    = -onehalfR * elem2R(edn,e2r_Length)
                                ! face downstream at link connection (node)
                                thisdata(3) = faceR(fdn,FaceCol_d)
                                thisX(3)    = -oneeighthR * elem2R(edn,e2r_Length)
                                ! face upstream at link connection (node)
                                thisdata(4) = faceR(fdn,FaceCol_u)
                                thisX(4)    = zeroR
                            case (fMultiple)
                                bdn => faceI(fdn,fi_branch_d)
                                ! downstram junction center
                                thisdata(1) = elemMR(edn,JunctionCol)
                                thisX(1)    = - onehalfR * elemMR(edn,eMr_Length)
                                ! downstrream branch
                                thisdata(2) = elemMR(edn,BranchColUp(bdn))
                                thisX(2)    = - onehalfR * elemMR(edn,eMr_LengthUp(bdn))
                                ! face downstream
                                thisdata(3) = faceR(fdn,FaceCol_d)
                                thisX(3)    = -oneeighthR * elemMR(edn,eMr_LengthUp(bdn))
                                ! face upstream
                                thisdata(4) = faceR(fdn,FaceCol_u)
                                thisX(4)    = zeroR
                            case (fBCup)
                                print *, 'error: fBCup not expected on a downstream face'
                                stop
                            case (fBCdn)

                                ! ghost element
                                thisdata(2) = elem2R(edn,ElemCol)
                                thisX(2)     = - onehalfR * elem2R(eStart,e2r_Length)
                                ! bc data
                                if (trim(outdataName) == 'eta') then
                                        thisdata(1) = bcdataDn(bcidxDn)%ThisValue
                                else
                                        thisdata(1) = thisdata(2)
                                endif
                                thisX(1)    = - elem2R(eStart,e2r_Length)

                                ! Face downstream
                                thisdata(3) = faceR(fdn,FaceCol_d)
                                thisX(3)    = -elem2R(eStart,e2r_Length) * 0.01
                                ! face upstream
                                thisdata(4) = faceR(fdn,FaceCol_u)
                                thisX(4)    = zeroR
                            case default
                                print *, trim(subroutine_name)
                                print *, ftypDn
                                print *, 'error: unexpected face type Dn of ',ftypDn,' in ',trim(subroutine_name)
                                stop
                        end select

                        !%  thread the data as element, upstream face, downstream face, element...
                        thisdata(5) = elem2R(eStart,ElemCol)
                        thisX(5)    = onehalfR * elem2R(eStart,e2r_Length)
                        en = eStart
                        do mm=0,nelem-1
                                fup => elem2I(en,e2i_Mface_u)
                                ! data on center
                                thisdata(3*mm+5) = elem2R(en,ElemCol)
                                thisX(3*mm+5)    = thisX(3*mm+4) + 0.5 * elem2R(en,e2r_Length)
                                ! downstream data on upstream face
                                thisdata(3*mm+6) = faceR(fup,FaceCol_d)
                                thisX(3*mm+6)    = thisX(3*mm+5) + 0.49 * elem2R(en,e2r_Length)
                                ! upstream data on upstream face
                                thisdata(3*mm+7) = faceR(fup,FaceCol_u)
                                thisX(3*mm+7)    = thisX(3*mm+5) + 0.5 * elem2R(en,e2r_Length)
                                en = en + 1
                        enddo
                        eLast = 3*nelem + 4

                        !% maps to the upstream face and element beyond this link
                        fup     => elem2I(eEnd,e2i_Mface_u)
                        ftypUp  => faceI(fup,fi_type)
                        eup     => faceI(fup,fi_Melem_u)
                        bcidxUp => faceI(fup,fi_BC_ID)


                        !    if (itemp == 2) then
                        !    if (ii == 2) then
                        !    print *, trim(subroutine_name)
                        !    print *, ftypUp
                        !    print *, fChannel, fPipe, fMultiple, fBCup, fBCdn
                        !    print *, eLast,'= eLast'
                        !    print *, eEnd,'=eEnd'
                        !    print *, eup,'=eup'
                        !    print *, bcidxUp,'=bcidxUp'
                        !    print *, size(thisdata,1),'=size thisdata'
                        !    print *, size(thisX,1),'=size thisx'
                        !    print *, size(elem2R,1),'=size elem2R'
                        !
                        !
                        !print *, thisdata(eLast+1)
                        !print *, elem2R(eup,ElemCol)
                        !print *, thisX(eLast+1)
                        !print *, thisX(eLast) + onehalfR * elem2R(eEnd,e2r_Length)
                        !!            ! BC
                        !            if (trim(outdataName) == 'flowrate') then
                        !print *, thisdata(eLast+2)
                        !print *, bcidxUp
                        !print *, size(bcdataUp)
                        !print *, bcdataUp(bcidxUp)%ThisValue
                        !            else
                        !!print *, thisdata(eLast+2)
                        !            endif
                        !!print *, thisX(eLast+2)
                        !!print *, thisX(eLast) +  elem2R(eEnd,e2r_Length)
                        !    stop
                        !    endif
                        !    endif
                        !
                        select case (ftypUp)
                            case (fChannel)
                                !%  First element of the next upstream link
                                thisdata(eLast+1) = elem2R(eup,ElemCol)
                                thisX(eLast+1)    = thisX(eLast) + onehalfR * elem2R(eup,e2r_Length)
                                thisdata(eLast+2) = faceR(elem2I(eup,e2i_Mface_u),FaceCol_d)
                                thisX(eLast+2)    = thisX(eLast) + elem2R(eup,e2r_Length)
                            case (fPipe)
                                print *, 'error: pipe not handled yet in ',trim(subroutine_name)
                                stop
                            case (fWeir)
                                !%  First element of the next upstream link
                                thisdata(eLast+1) = elem2R(eup,ElemCol)
                                thisX(eLast+1)    = thisX(eLast) + onehalfR * elem2R(eup,e2r_Length)
                                thisdata(eLast+2) = faceR(elem2I(eup,e2i_Mface_u),FaceCol_d)
                                thisX(eLast+2)    = thisX(eLast) + elem2R(eup,e2r_Length)
                            case (fOrifice)
                                !%  First element of the next upstream link
                                thisdata(eLast+1) = elem2R(eup,ElemCol)
                                thisX(eLast+1)    = thisX(eLast) + onehalfR * elem2R(eup,e2r_Length)
                                thisdata(eLast+2) = faceR(elem2I(eup,e2i_Mface_u),FaceCol_d)
                                thisX(eLast+2)    = thisX(eLast) + elem2R(eup,e2r_Length)
                            case (fMultiple)
                                !% an upstream junction
                                bup => faceI(fup,fi_branch_u)
                                ! branch values
                                thisdata(eLast+1) = elemMR(eup,BranchColDn(bup))
                                thisX(eLast+1)    = thisX(eLast) + onehalfR * elemMR(eup,eMr_LengthDn(bup))
                                ! junction center values
                                thisdata(eLast+2) = elemMR(eup,JunctionCol)
                                thisX(eLast+2)    = thisX(eLast) + onehalfR * elemMR(eup,eMr_Length)

                            case (fBCup)
                                ! an upstream BC with ghost element
                                thisdata(eLast+1) = elem2R(eup,ElemCol)
                                thisX(eLast+1)    = thisX(eLast) + onehalfR * elem2R(eEnd,e2r_Length)
                                ! BC
                                if (trim(outdataName) == 'flowrate') then
                                        thisdata(eLast+2) = bcdataUp(bcidxUp)%ThisValue
                                else
                                        thisdata(eLast+2) = thisdata(eLast+1)
                                endif
                                thisX(eLast+2)    = thisX(eLast) +  elem2R(eEnd,e2r_Length)
                            case (fBCdn)
                                print *, 'error: fBCdn not expected on an upstream face in ',trim(subroutine_name)
                                stop
                            case default
                                print *, trim(subroutine_name)
                                print *, ftypUp
                                print *, 'error: unexpected face type Up of ',ftypUp,' in ',trim(subroutine_name)
                                stop
                        end select
                        ndata = eLast + 2

                        !    do mm=1,ndata
                        !        print *, mm, thisX(mm), thisdata(mm)
                        !    end do

                        write(thisUnit,*)  thisstep,'=this step'
                        write(thisUnit,*)  Lidx,'=this_link_index'
                        write(thisUnit,*)  ndata,'=items_this_link'
                        write(thisUnit,*)  2,'=rows_this_link_X_data'
                        write(thisUnit,*)  thisX(1:ndata)
                        write(thisUnit,*)  thisdata(1:ndata)


                end do

                if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
        end subroutine output_one_threaded_data_by_link
        !==========================================================================
        !
        ! PRIVATE BELOW HERE
        !
        !==========================================================================
        !
        subroutine output_singlethreadedfile_open (threadedfile)
                !
                ! Opens a single file for output writing
                !
                character(64) :: subroutine_name = 'output_singlethreadedfile_open'

                type(threadedfileType), intent(in out) :: threadedfile

                integer                :: open_status
                character(len=512)     :: emsg

                !--------------------------------------------------------------------------
                if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

                open_status = 0

                threadedfile%FileInfo%UnitNumber = outputfile_next_unitnumber
                outputfile_next_unitnumber       = outputfile_next_unitnumber+1

                open(unit=threadedfile%FileInfo%Unitnumber, &
                        file=trim(threadedfile%FileInfo%WriteName), &
                        status = trim(threadedfile%FileInfo%FileStatus), &
                        access = 'sequential', &
                        form   = 'formatted', &
                        action = 'write', &
                        iostat = open_status)

                emsg = 'file exists or path/folder does not exist: file open failed in '//trim(subroutine_name) &
                        // '; filename = '//trim(threadedfile%FileInfo%WriteName)
                call utility_check_fileopen (open_status, emsg)

                if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
        end subroutine output_singlethreadedfile_open
        !
        !==========================================================================
        ! END OF MODULE output
        !==========================================================================
end module output
