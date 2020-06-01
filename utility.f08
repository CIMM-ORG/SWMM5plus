! module utility
!
! utility routines that may be called in a number of places
!
!==========================================================================
!
module utility

        use array_index
        use data_keys
        use setting_definition
        use globals

        implicit none

        private

        public  :: utility_advance_temp_array
        public  :: utility_average
        public  :: utility_check_allocation
        public  :: utility_check_fileopen
        public  :: utility_linear_interpolate_within_indexlist
        public  :: utility_get_datetime_stamp
        public  :: utility_print_values_by_link
        public  :: utility_round_to_significant_digits
        public  :: utility_scale_of_number
        public  :: utility_sign_with_ones

        integer :: debuglevel = 0

contains
        !
        !==========================================================================
        !==========================================================================
        !
        subroutine utility_print_values_by_link &
                (elem2R, elem2I, elemMR, elemMI, faceR, faceI, ilink, &
                fcolDn, fcolUp, e2col, eMcolC, eMcolDn, eMcolUp)
                !
                ! provides a printout of data along a link.
                ! fcolDn, fcolUp are (for example) fr_Area_d and fr_Area_u. If only
                ! one value (e.g. fr_Flowrate) then use for fcolDn and use 0 for fcolUp
                ! e2col is (for example) e2r_Area.
                ! eMcolC is (for example) eMr_Area.
                ! eMcolDn and eMcolUp are (for example) eMr_AreaUp and eMr_AreaDn.
                ! if eM up/dn do not exist (e.g. for eMr_Eta) then store a dummyarray of 0.
                ! note that the input dummy array must be an array rather than a scalar.
                !
                character(64) :: subroutine_name = 'utility_print_values_by_link'

                real,      intent(in)      :: elem2R(:,:), elemMR(:,:), faceR(:,:)
                integer,   intent(in)      :: elem2I(:,:), elemMI(:,:), faceI(:,:)
                integer,   intent(in)      :: ilink, fcolDn, fcolUp, e2col, eMcolC, eMcolDn(:), eMcolUp(:)

                integer :: eStart, eEnd, ii, fd, fu, eb, eu, ed, edtype, eutype
                !--------------------------------------------------------------------------
                if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

                eStart = minloc(elem2I(:,e2i_link_Pos),1,(  (elem2I(:,e2i_link_ID) == ilink) &
                        .and. (elem2I(:,e2i_link_Pos) /= nullvalueI) ) )
                eEnd   = maxloc(elem2I(:,e2i_link_Pos),1,(  (elem2I(:,e2i_link_ID) == ilink) &
                        .and. (elem2I(:,e2i_link_Pos) /= nullvalueI) ) )

                print *, ilink,'=link'

                !% info about the element downstream of the link
                fd = elem2I(eStart,e2i_Mface_d)
                edtype = faceI(fd,fi_etype_d)
                print *, edtype,'=type of downstream element'
                if (edtype == eJunctionChannel) then
                        ed = faceI(fd,fi_Melem_d)
                        eb = faceI(fd,fi_branch_d)
                        if (eMcolC /= 0) then
                                print *, elemMR(ed,eMcolC),'= elemM C ',ed
                        endif
                        if (eMcolUp(1) /= 0) then
                                print *, elemMR(ed,eMcolUp(eb)),'= branch up,',eb
                        endif
                endif

                !% print information about faces and elements of the link
                do ii = eStart,eEnd
                        fd = elem2I(ii,e2i_Mface_d)
                        print *, faceR(fd,fcolDn),'= faceDn ',fd
                        if (fcolUp /= 0) then
                                print *, faceR(fd,fcolUp),'= faceUp ',fd
                        endif
                        print *, elem2R(ii,e2col),'= elem   ',ii
                end do
                !% print info on the last upstream face
                fu = elem2I(eEnd,e2i_Mface_u)
                print *, faceR(fu,fcolDn),'= faceDn ',fu
                if (fcolUp /= 0) then
                        print *, faceR(fu,fcolUp),'= faceUp ',fu
                endif

                !% print info of the next element upstream
                eutype = faceI(fu,fi_etype_u)
                print *, eutype,'=type of upstream element'
                if (eutype == eJunctionChannel) then
                        eu = faceI(fu,fi_Melem_u)
                        eb = faceI(fu,fi_branch_u)
                        if (eMcolDn(1) /= 0) then
                                print *, elemMR(eu,eMcolDn(eb)),'= branch dn,',eb
                        endif
                        if (eMcolC /= 0) then
                                print *, elemMR(eu,eMcolC),'= elemM C ',eu
                        endif
                endif


                if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
        end subroutine utility_print_values_by_link
        !
        !==========================================================================
        !==========================================================================
        !
        integer function utility_advance_temp_array (next_temparray, temparraysize)
                !
                ! Advances a temp_array counter and checks to see if the array size is
                ! adequate. These arrays are columns of the 2D data structure that can be
                ! overwritten in any subroutine.
                !
                ! The temp_array (e.g. e2i_Temp) have sizes that are hard-coded in the
                ! array_index module and then are initialized in the initialization module.
                ! Care must be taken when increasing the temp_array index sizes for consistency.
                !
                integer, intent(in) :: next_temparray, temparraysize
                character(64) :: subroutine_name = 'utility_advance_temp_array'

                !--------------------------------------------------------------------------
                if (next_temparray > temparraysize) then
                        print *, 'code error: temparray is too small in ',subroutine_name
                        stop
                endif

                utility_advance_temp_array = next_temparray+1

        end function utility_advance_temp_array
        !
        !==========================================================================
        !==========================================================================
        !
        pure real function utility_average &
                (inarray) result(outscalar)
                !
                ! computes a simple average of an array
                !
                real,  intent(in)  :: inarray(:)

                !--------------------------------------------------------------------------

                if (size(inarray) > 0 ) then
                        outscalar = sum(inarray) / real(size(inarray))
                endif


        end function utility_average
        !
        !==========================================================================
        !==========================================================================
        !
        subroutine utility_check_fileopen &
                (open_status, emsg)
                !
                ! Checks to see if a file is already open. Stops as an error if file is open.
                !
                character(64) :: subroutine_name = 'utility_check_fileopen'

                integer,           intent(in)   :: open_status
                character(len=*),  intent(in)   :: emsg

                !--------------------------------------------------------------------------
                if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

                if (open_status > 0) then
                        print *, trim(emsg)
                        STOP
                endif

                if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
        end subroutine utility_check_fileopen
        !
        !==========================================================================
        !==========================================================================
        !
        subroutine utility_check_allocation &
                (allocation_status, emsg)
                !
                ! checks allocation status and stops if there is an error
                !
                character(64) :: subroutine_name = 'utility_check_allocation'

                integer,           intent(in)   :: allocation_status
                character(len=*),  intent(in)   :: emsg

                !--------------------------------------------------------------------------
                if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

                if (allocation_status > 0) then
                        print *, trim(emsg)
                        STOP
                endif

                if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
        end subroutine utility_check_allocation
        !
        !==========================================================================
        !==========================================================================
        !
        function utility_linear_interpolate_within_indexlist &
                (thisIndex, IndexArray,  ValueArray) result(thisValue)
                !
                ! IndexArray are ordered, non-repeating real numbers
                ! ValueArray are non-ordered values corresponding to IndexArray
                ! thisIndex is an index value within limits of IndexArray
                ! thisValue is a linearly-interpolated value
                !
                ! This function is used in boundary condition interpolation, but is more
                ! general and can be used elsewhere
                !
                character(64) :: subroutine_name = 'utility_linear_interpolate_within_indexlist'

                real,      intent(in)  :: IndexArray(:), ValueArray(:)
                real,      intent(in)  :: thisIndex

                real ::  thisValue

                integer :: closeloc

                !--------------------------------------------------------------------------
                if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

                closeloc = minloc(abs(thisIndex - IndexArray),1)

                !% get the value if the location is exactly on an index
                if (thisIndex - IndexArray(closeloc) == 0) then
                        thisValue = ValueArray(closeloc)

                        !% if the close location index is below the target
                elseif ((thisIndex - IndexArray(closeloc) > 0) .and. &
                        (closeloc < size(IndexArray,1)) ) then
                        thisValue = ( ValueArray(closeloc)   * ( IndexArray(closeloc + 1) - thisIndex)    &
                                +ValueArray(closeloc+1) * ( thisIndex - IndexArray(closeloc)    ) )  &
                                / ( IndexArray(closeloc+1) - IndexArray(closeloc) )

                        !% if the close location is above the target
                elseif ((thisIndex - IndexArray(closeloc) < 0) .and. &
                        (closeloc > 1) ) then
                        thisValue = ( ValueArray(closeloc)   * ( IndexArray(closeloc) - thisIndex   )    &
                                +ValueArray(closeloc-1) * ( thisIndex - IndexArray(closeloc-1) )  ) &
                                / ( IndexArray(closeloc) - IndexArray(closeloc - 1) )

                else
                        !% error condition - the index appears to be outside of the index array
                        print *, 'error: interpolation index (thisIndex) is outside the array bounds (IndexArray) in ' &
                                , subroutine_name
                        stop
                endif

                if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
        end function utility_linear_interpolate_within_indexlist
        !
        !==========================================================================
        !==========================================================================
        !
        subroutine utility_get_datetime_stamp &
                (datetimestamp)
                !
                !   Gets a character(14) datetimestamp in the form _yyyymmdd_hhmm .
                !   Note the leading _ is used so this is convenient for adding to filenames.
                !
                character(64) :: subroutine_name = 'utility_get_datetime_stamp'

                character(14), intent(out) :: datetimestamp ! in form _yyyymmdd_hhmm

                character(8) :: thisdate
                character(10):: thistime
                character(5) :: thiszone

                character(4) :: cyear
                character(2) :: cmonth
                character(2) :: cday
                character(2) :: chour
                character(2) :: cmin

                integer      :: thisvalues(8)
                integer      :: tyear, tmonth, tday, thour, tmin

                !--------------------------------------------------------------------------
                if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

                call date_and_time(thisdate, thistime, thiszone, thisvalues)

                tyear  = thisvalues(1)
                tmonth = thisvalues(2)
                tday   = thisvalues(3)
                thour  = thisvalues(5)
                tmin   = thisvalues(6)

                write(cyear,"(i4)") tyear

                if (tmonth >= 10) then
                        write(cmonth,"(i2)") tmonth
                else
                        write(cmonth,"(i1,i1)") 0, tmonth
                end if

                if (tday >= 10) then
                        write(cday,"(i2)") tday
                else
                        write(cday,"(i1,i1)") 0, tday
                end if

                if (thour >= 10) then
                        write(chour,"(i2)") thour
                else
                        write(chour,"(i1,i1)") 0, thour
                end if

                if (tmin >= 10) then
                        write(cmin,"(i2)") tmin
                else
                        write(cmin,"(i1,i1)") 0, tmin
                end if

                write(datetimestamp,"('_',a4,a2,a2,'_',a2,a2)") &
                        cyear, cmonth, cday, &
                        chour, cmin

                if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
        end subroutine utility_get_datetime_stamp
        !
        !==========================================================================
        !==========================================================================
        !
        pure elemental real function utility_round_to_significant_digits &
                (inarray,idigits) result (outarray)
                !
                ! returns the exponent scale of the largest n where 10^n < number
                !
                real,      intent(in)    :: inarray
                integer,   intent(in)    :: idigits

                integer :: inscale

                !--------------------------------------------------------------------------

                inscale = utility_scale_of_number(inarray)

                outarray = real(NINT(inarray / (tenR**(inscale - (idigits-1))))) &
                        * tenR**real(inscale - (idigits-1))


        end function utility_round_to_significant_digits
        !
        !==========================================================================
        !==========================================================================
        !
        pure elemental integer function utility_scale_of_number &
                (inarray) result (outarray)
                !
                ! returns the exponent scale of the largest n where 10^n < number
                !
                real,      intent(in)    :: inarray

                !--------------------------------------------------------------------------

                outarray = floor(log10(abs(inarray)))

        end function utility_scale_of_number
        !
        !==========================================================================
        !==========================================================================
        !
        pure elemental real function utility_sign_with_ones &
                (inarray) result (outarray)
                !
                ! returns is an array of real ones with the sign of the inarray argument
                !
                real,      intent(in)    :: inarray

                !--------------------------------------------------------------------------

                outarray = oneR
                outarray = sign(outarray,inarray)

        end function utility_sign_with_ones
        !
        !==========================================================================
        ! END OF MODULE utility
        !==========================================================================
end module utility
