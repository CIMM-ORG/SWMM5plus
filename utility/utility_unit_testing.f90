module utility_unit_testing
   !%==========================================================================
   !% SWMM5+ release, version 1.0.0
   !% 20230608
   !% Hydraulics engine that links with EPA SWMM-C
   !% June 8, 2023
   !%
   !% Description:
   !% The public subroutines are NOT normally called within SWMM5+
   !% Instead, these are designed to be called by programmers at almost any
   !% location by including USE utility_unit_testing in the preamble and
   !% then calling one of the public subroutine.
   !%
   !% The util_utest_CLprint() is a blank routine so that programmers can
   !% write their own command-line print statements to examine code behavior
   !%==========================================================================

   use interface_
   use utility_allocate
   use discretization
   use define_indexes
   use define_keys
   use define_globals
   use define_settings
   use storage_geometry, only : storage_implied_depth_from_volume

   use utility_crash, only : util_crashpoint

   implicit none

   private

   public :: util_utest_CLprint     !% custom command line printing for debugging
   public :: util_utest_checkIsNan  !% checks if selected data columns have NaN
   public :: util_utest_local_global !% checks that indexes are unique
   public :: util_utest_pack_arrays !% checks that pack arrays are unique
   public :: util_utest_node_link_image !% checks all images have nodes
   public :: util_utest_slope_checking  !% finds negative slopes
   public :: util_utest_global_index_check !% check global indexes

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine util_utest_CLprint (inputstring)
      !%------------------------------------------------------------------
      !% Description:
      !% Command-Line printing
      !% Blank subroutine for use for command-line write during debugging
      !%------------------------------------------------------------------
      !% Declarations:
         character (len=*), intent(in) :: inputstring

         integer, pointer :: fup(:), fdn(:), eup(:), edn(:)
         real(8), pointer :: dt, oneVec(:), grav
      !%------------------------------------------------------------------
      !% Preliminaries:
      !%------------------------------------------------------------------
      !% Aliases:
         fup => elemI(:,ei_Mface_uL)
         fdn => elemI(:,ei_Mface_dL)
         eup => faceI(:,fi_Melem_uL)
         edn => faceI(:,fi_Melem_dL)
         dt  => setting%Time%Hydraulics%Dt
         grav => setting%Constant%gravity
         oneVec   => elemR(:,er_ones)
      !%------------------------------------------------------------------

      !% --- USEFUL HEADER -----------------------------------------------
         print *, ' '
         write(*,"(A,A, e12.5)") ' ',trim(inputstring)     
         write(*,"(A,i7,A, f12.5, A, f12.5, A)") '        step = ',setting%Time%Step ,&
         '; dt = ',setting%Time%Hydraulics%Dt,&
         '; time = ',setting%Time%Now / 60.d0 , 'min'    
         
         return 

   end subroutine util_utest_CLprint   
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_utest_checkIsNan ()
      !%--------------------------------------------------------------------
      !% Description: 
      !% tests for NaN in arrays defined by eIsNan_... indexes
      !% Can be called anywhere to check whether arrays have NaN
      !%--------------------------------------------------------------------
      !% Declarations
         integer :: thisCol, ii, mm
         logical :: isThisElemNan(Ncol_elemIsNan) 
         logical :: isThisFaceNan(Ncol_faceIsNan)
         logical :: foundNanElem = .false.
         logical :: foundNanFace = .false.
         character (64)  :: colIndexName = ' '
      !%--------------------------------------------------------------------

      elemIsNan(:,:) = .false.
      faceIsNan(:,:) = .false.

      !% --- handle elements
      do ii = 1,Ncol_elemIsNan
         !% --- select the column in elemR(:,:)
         !%     this matches the eIsNan_... index to the er_... index
         thisCol = util_utest_get_elemR_col(ii)

         !% --- get the data name for the column in elemR
         ! call util_utest_get_elemR_indexName (ii, colIndexName)

         !% --- check for any NaN in the array
         isThisElemNan(ii) = util_utest_isThisCol_Nan(thisCol,.true.)

         if (isThisElemNan(ii)) then 
            !% ---- store locations of NaN
            elemIsNan(:,ii) = isnan(elemR(:,thisCol))
            foundNanElem = .true.
         else
            elemIsNan(:,ii) = .false.
         end if
      end do

      !% --- handle faces
      do ii = 1,Ncol_faceIsNan
         !% --- select the column in faceR(:,:)
         !%     this matches the fIsNan_... index to the fr_... index
         thisCol = util_utest_get_faceR_col(ii)

         !% --- get the data name for the column in faceR
         ! call util_utest_get_faceR_indexName (ii, colIndexName)

         !% --- check for any NaN in the array
         isThisFaceNan(ii) = util_utest_isThisCol_Nan(thisCol, .false.)

         if (isThisFaceNan(ii)) then 
            !% ---- store locations of NaN
            faceIsNan(:,ii) = isnan(faceR(:,thisCol))
            foundNanFace = .true.
         else
            faceIsNan(:,ii) = .false.
         end if
      end do

      !% --- report Nan
      if (foundNanElem) then 
         !% --- cycle through the checked columns
         do ii = 1,Ncol_elemIsNan
            if (isThisElemNan(ii)) then 
               !% --- get the column data name
               colIndexName = ' '
               ! call util_utest_get_elemR_indexName (ii, colIndexName)
               
               print *, ' '
               print *, 'CODE STOPPING DUE TO NaN in ',trim(colIndexName)
               print *, 'On processor image = ',this_image()
               print *, 'The following element indexes are involved '
               do mm=1,N_elem(this_image())
                  if (elemIsNan(mm,ii)) then 
                     !% --- list the element numbers with NaN values
                     print *, 'elem = ', mm
                  end if
                  print *, ' '
               end do
            end if
         end do
      end if

      if (foundNanFace) then 
         !% --- cycle through the checked columns
         do ii = 1,Ncol_faceIsNan
            if (isThisFaceNan(ii)) then 
               !% --- get the column data name
               colIndexName = ' '
               ! call util_utest_get_faceR_indexName (ii, colIndexName)
               
               print *, ' '
               print *, 'CODE STOPPING DUE TO NaN in ',trim(colIndexName)
               print *, 'On processor image = ',this_image()
               print *, 'The following face indexes are involved '
               do mm=1,N_face(this_image())
                  if (faceIsNan(mm,ii)) then 
                     !% --- list the face numbers with NaN values
                     print *, 'face = ', mm
                  end if
               end do
               print *, ' '
            end if
         end do
      end if     

      if ((foundNanFace) .or. (foundNanElem)) then
         call util_crashpoint(729873)
      end if
      
    end subroutine util_utest_checkIsNan
!%
!%==========================================================================
!%==========================================================================
!%
   subroutine util_utest_local_global
      !%------------------------------------------------------------------
      !% Description
      !% Checking the the local and global indexs 
      !% of the link,node,elem and face arrays are unique
      !% Can be called anywhere for debuggin5
      !%------------------------------------------------------------------
      !% Declarations
         integer ii, jj, kk, min_val, max_val
         logical dup_found
         character(64) :: subroutine_name = 'local_global_unique'
      !%------------------------------------------------------------------
      !% Preliminaries
         if (setting%Debug%File%initialization) &
                  write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
      !%------------------------------------------------------------------
      !% --- Looping through the link array and finding all of the unqiue values
      min_val = minval(link%I(:,li_idx)) - 1
      max_val = maxval(link%I(:,li_idx))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(link%I(:,li_idx),mask=link%I(:,li_idx)>min_val)
      end do

      if (ii /= size(link%I(:, li_idx))) then
         print *, "ERROR:::: link%I(:,li_idx) is not unique. This_image ::", this_image()
      else
         print *, "link%I(:,li_idx) is unique. This_image ::", this_image()
      end if

      !% --- checking node%I(:,:) indexes
      min_val = minval(node%I(:,ni_idx)) - 1
      max_val = maxval(node%I(:,ni_idx))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(node%I(:,ni_idx),mask=node%I(:,ni_idx)>min_val)
      end do

      if (ii /= size(node%I(:, ni_idx))) then
         print *, "ERROR:::: node%I(:,ni_idx) is not unique. This_image ::", this_image()
      else
         print *, "node%I(:,ni_idx) is unique. This_image ::", this_image()
      end if


      !% --- checking elemI(:,:) local indexes
      min_val = minval(elemI(1:N_elem(this_image()),ei_Lidx)) - 1
      max_val = maxval(elemI(1:N_elem(This_image()),ei_Lidx))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemI(1:N_elem(This_image()),ei_Lidx),mask=elemI(1:N_elem(This_image()),ei_Lidx)>min_val)
      end do

      if (ii /= n_elem(this_image())) then
         print *, "ERROR:::: elemI(:,ei_Lidx) is not unique. This_image ::", this_image()
      else
         print *, "elemI(:,ei_Lidx) is unique. This_image ::", this_image()
      end if

      !% --- checking faceI(:,:) indexes
      min_val = minval(faceI(1:N_face(this_image()),fi_Lidx)) - 1
      max_val = maxval(faceI(1:N_face(this_image()),fi_Lidx))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(faceI(1:N_face(this_image()),fi_Lidx),mask=faceI(1:N_face(this_image()),fi_Lidx)>min_val)
      end do

      if (ii /= n_face(this_image())) then
         print *, "ERROR:::: faceI(:,fi_Lidx) is not unique. This_image ::", this_image()
      else
         print *, "faceI(:,fi_Lidx) is unique. This_image ::", this_image()
      end if

      !% --- checking faceI(:,:) global indexes
      min_val = minval(faceI(1:N_face(this_image()),fi_Gidx)) - 1
      max_val = maxval(faceI(1:N_face(this_image()),fi_Gidx))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(faceI(1:N_face(this_image()),fi_Gidx),mask=faceI(1:N_face(this_image()),fi_Gidx)>min_val)
      end do

      if (ii /= n_face(this_image())) then
         print *, "ERROR:::: faceI(:,fi_Gidx) is not unique. This_image ::", this_image()
      else
         print *, "faceI(:,fi_Gidx) is unique. This_image ::", this_image()
      end if


      !% --- checking elemI(:,:) global indexes
      min_val = minval(elemI(1:N_elem(this_image()),ei_Gidx)) - 1
      max_val = maxval(elemI(1:N_elem(this_image()),ei_Gidx))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemI(1:N_elem(this_image()),ei_Gidx),mask=elemI(1:N_elem(this_image()),ei_Gidx)>min_val)
      end do

      if (ii /= N_elem(this_image())) then
         print *, "ERROR:::: elemI(:,ei_Gidx) is not unique. This_image ::", this_image()
      else
         print *, "elemI(:,ei_Gidx) is unique. This_image ::", this_image()
      end if

      !%------------------------------------------------------------------
      !% Closing
         if (setting%Debug%File%initialization)  &
                  write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

  end subroutine util_utest_local_global
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_utest_pack_arrays
      !%------------------------------------------------------------------
      !% Description
      !% Going through all of the pack arrays and making sure they are unique
      !%------------------------------------------------------------------
      !% Declarations
         integer ii, jj, kk, min_val, max_val
         logical dup_found
         character(64) :: subroutine_name = 'pack_arrays_unique'
      !%------------------------------------------------------------------
      !% Preliminaries
         if (setting%Debug%File%initialization) &
               write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
      !%------------------------------------------------------------------

      kk = 1
      dup_found = .false.

      !% --- checking elemP(:,ep_ALLtm) indexes
      min_val = minval(elemP(:,ep_CCJM)) - 1
      max_val = maxval(elemP(:,ep_CCJM))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_CCJM),mask=elemP(:,ep_CCJM)>min_val)
      end do

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if

      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CCJM) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemP(ep_CCJM)) then
         print *, "ERROR:::: elemP(:,ep_CCJM) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CCJM) is unique. This_image ::", this_image()
      end if

      !% --- checking elemP(:,ep_CC) indexes

      min_val = minval(elemP(:,ep_CC)) - 1
      max_val = maxval(elemP(:,ep_CC))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_CC),mask=elemP(:,ep_CC)>min_val)
      end do

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if

      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CC) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemP(ep_CC)) then
         print *, "ERROR:::: elemP(:,ep_CC) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CC) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_CC_H) indexes

      min_val = minval(elemP(:,ep_CC_H)) - 1
      max_val = maxval(elemP(:,ep_CC_H))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_CC_H),mask=elemP(:,ep_CC_H)>min_val)
      end do

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if

      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CC_H) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemP(ep_CC_H)) then
         print *, "ERROR:::: elemP(:,ep_CC_H) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CC_H) is unique. This_image ::", this_image()
      end if

      !% ---checking elemP(:,ep_CC_Q) indexes
      min_val = minval(elemP(:,ep_CC_Q)) - 1
      max_val = maxval(elemP(:,ep_CC_Q))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_CC_Q),mask=elemP(:,ep_CC_Q)>min_val)
      end do

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if

      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CC_Q) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemP(ep_CC_Q)) then
         print *, "ERROR:::: elemP(:,ep_CC_Q) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CC_Q is unique. This_image ::", this_image()
      end if

      !% --- checking elemP(:,ep_Diag) indexes
      min_val = minval(elemP(:,ep_Diag)) - 1
      max_val = maxval(elemP(:,ep_Diag))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_Diag),mask=elemP(:,ep_Diag)>min_val)
      end do

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if

      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_Diag) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemP(ep_Diag)) then
         print *, "ERROR:::: elemP(:,ep_Diag) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_Diag) is unique. This_image ::", this_image()
      end if

   
      !% --- checking elemPGetm(:,epg_CC_rectangular) indexes
      min_val = minval(elemPGetm(:,epg_CC_rectangular)) - 1
      max_val = maxval(elemPGetm(:,epg_CC_rectangular))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemPGetm(:,epg_CC_rectangular),&
              mask=elemPGetm(:,epg_CC_rectangular)>min_val)
      end do

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if

      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemPGetm(:,epg_CC_rectangular) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemPGetm(epg_CC_rectangular)) then
         print *, "ERROR:::: elemPGetm(:,epg_CC_rectangular) is not unique. This_image ::", this_image()

      else
         print *, "elemPGetm(:,epg_CC_rectangular) is unique. This_image ::", this_image()
      end if

      !% --- checking elemPGetm(:,epg_CC_trapezoidal) indexes
      min_val = minval(elemPGetm(:,epg_CC_trapezoidal)) - 1
      max_val = maxval(elemPGetm(:,epg_CC_trapezoidal))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemPGetm(:,epg_CC_trapezoidal), &
              mask=elemPGetm(:,epg_CC_trapezoidal)>min_val)
      end do

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if

      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemPGetm(:,epg_CC_trapezoidal) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemPGetm(epg_CC_trapezoidal)) then
         print *, "ERROR:::: elemPGetm(:,epg_CC_trapezoidal) is not unique. This_image ::", this_image()

      else
         print *, "elemPGetm(:,epg_CC_trapezoidal) is unique. This_image ::", this_image()
      end if

      !% --- checking faceP(:,fp_noBC_IorS) indexes
      min_val = minval(faceP(:,fp_noBC_IorS)) - 1
      max_val = maxval(faceP(:,fp_noBC_IorS))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(faceP(:,fp_noBC_IorS),mask=faceP(:,fp_noBC_IorS)>min_val)
      end do

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if

      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "faceP(:,fp_noBC_IorS) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_faceP(fp_noBC_IorS)) then
         print *, "ERROR:::: faceP(:,fp_noBC_IorS) is not unique. This_image ::", this_image()

      else
         print *, "faceP(:,fp_noBC_IorS) is unique. This_image ::", this_image()
      end if

      !% --- checking faceP(:,fp_Diag_IorS) indexes
      min_val = minval(faceP(:,fp_Diag_IorS)) - 1
      max_val = maxval(faceP(:,fp_Diag_IorS))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(faceP(:,fp_Diag_IorS),mask=faceP(:,fp_Diag_IorS)>min_val)
      end do

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if

      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "faceP(:,fp_Diag_IorS) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_faceP(fp_Diag_IorS)) then
         print *, "ERROR:::: faceP(:,fp_Diag_IorS) is not unique. This_image ::", this_image()

      else
         print *, "faceP(:,fp_Diag_IorS) is unique. This_image ::", this_image()
      end if

      !% --- checking faceP(:,fp_JumpDn_IorS) indexes
      min_val = minval(faceP(:,fp_JumpDn_IorS)) - 1
      max_val = maxval(faceP(:,fp_JumpDn_IorS))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(faceP(:,fp_JumpDn_IorS),mask=faceP(:,fp_JumpDn_IorS)>min_val)
      end do

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if

      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "faceP(:,fp_JumpDn_IorS) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_faceP(fp_JumpDn_IorS)) then
         print *, "ERROR:::: faceP(:,fp_JumpDn_IorS) is not unique. This_image ::", this_image()

      else
         print *, "faceP(:,fp_JumpDn_IorS) is unique. This_image ::", this_image()
      end if

      !% --- checking faceP(:,fp_JumpUp_IorS) indexes
      min_val = minval(faceP(:,fp_JumpUp_IorS)) - 1
      max_val = maxval(faceP(:,fp_JumpUp_IorS))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(faceP(:,fp_JumpUp_IorS),mask=faceP(:,fp_JumpUp_IorS)>min_val)
      end do

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if

      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "faceP(:,fp_JumpUp_IorS) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_faceP(fp_JumpUp_IorS)) then
         print *, "ERROR:::: faceP(:,fp_JumpUp_IorS) is not unique. This_image ::", this_image()

      else
         print *, "faceP(:,fp_JumpUp_IorS) is unique. This_image ::", this_image()
      end if


      if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_utest_pack_arrays
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_utest_node_link_image
      !%------------------------------------------------------------------
      !% Description
      !% Check whether the every image has at least one node and link 
      !% assigned to it.
      !%------------------------------------------------------------------
      !% Declarations
         integer :: ii, jj, kk, counter
         character(64) :: subroutine_name = 'init_face_check'
      !%------------------------------------------------------------------
      !% Preliminaries:
         if (setting%Debug%File%initialization) &
               write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
      !%------------------------------------------------------------------

      kk = 1
      counter = 0

      !% --- We can find if there is a node on each image by counting the amount 
      !%     of unique values there is  inside of node%I(:,ni_P_image)
      !%     So we use the code below to loop through node%I(:,ni_P_image) and 
      !%     find all the unique values

      do ii = 1, size(node%I(:,ni_P_image))
         do jj = 1, ii
            if ((node%I(ii, ni_P_image)) == node%I(jj, ni_P_image)) then
               exit
            end if
         end do
            if (ii == jj) then
               counter = counter + 1
            end if
      end do

      !% --- After that we compare the number of images we are using to what 
      !%     we counted. If it is correct they should be the same value, 
      !%     otherwise there was an error when paritioning the links and nodes

      if (num_images() /= counter) then
         print *, "error in NodeI images. This_image :: ", this_image()
      else
         print *, "correct number in node%I images. This_image :: ", this_image()
      end if

      !% --- We reset the counter and do the same process for the Links
      counter = 0

      do ii = 1, size(link%I(:,li_P_image))
         do jj = 1, ii
            if ((link%I(ii, li_P_image)) == link%I(jj, li_P_image)) then
               exit
            end if
         end do
            if (ii == jj) then
               counter = counter + 1

            end if
      end do

      if (num_images() /= counter) then
         print *, "error in link%I images. This_image :: ", this_image()
         print *, "counter", counter
      else
         print *, "correct number in link%I images.  This_image :: ", this_image()
      end if

      !%------------------------------------------------------------------
      !% Closing
         if (setting%Debug%File%initialization)  &
               write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_utest_node_link_image
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_utest_slope_checking
      !%------------------------------------------------------------------
      !% Description
      !% Checking that all of the slopes are postive.
      !% To do this and loop through and if we find a negative slope 
      !% then we exit the loop and report it.
      !% Note that SWMM5+ handles negative slopes without any changes
      !% but sometimes it is useful to know where the negative slopes
      !% occur because they cause EPA-SWMM-C to reverse the slope, which
      !% alters the comparison between models.
      !%------------------------------------------------------------------
      !% Description
         integer :: ii, jj
         logical :: invalid_slope
         character(64) :: subroutine_name = 'slope_checking'
      !%------------------------------------------------------------------
      !% Preliminaries
         if (setting%Debug%File%initialization) &
               write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
      !%------------------------------------------------------------------

      invalid_slope = .false.
      do ii = 1, size(link%R(:,lr_Slope))
         if (link%R(ii, lr_slope) < 0) then
            invalid_slope = .true.
            exit
         end if
      end do

      if (invalid_slope .eqv. .true.) then
         print *, "error found in link%R(:,lr_slope) slope is negative. This_image :: ", this_image()
      else
         print *, "all slopes are postive.  This_image :: ", this_image()
      end if

      !%------------------------------------------------------------------
      !% Closing
         if (setting%Debug%File%initialization)  &
               write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_utest_slope_checking
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_utest_global_index_check
      !%------------------------------------------------------------------
      !% Description
      !% Checking that the global indexs are correct by adding the first 
      !% valid global index of an image with the local index
      !%------------------------------------------------------------------
      !% Declarations
         integer :: ii, current_length, counter
         character(64) :: subroutine_name = 'global_index_checking'
      !%------------------------------------------------------------------
      !% Preliminaries
         if (setting%Debug%File%initialization) &
               write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
      !%------------------------------------------------------------------
     
      !% --- find the current length of the global index by looking at the 
      !%     first value of elemI on that image and subtracting one.

      current_length = elemI(1, ei_Gidx) - 1

      !% --- loop through elemI(:,ei_Gidx) and report and stop looping 
      !%     if there is an error in the global indexes
      do ii = 1, size(elemI(:,ei_Gidx))
         if (elemI(ii,ei_Gidx) /= current_length+elemI(ii,ei_Lidx) .and. elemI(ii,ei_Gidx)/= nullvalueI) then
            print *, "error in elem global indexes. Processor :: ", this_image()
            !%print *, "elemI(ii,ei_Gidx)", elemI(ii,ei_Gidx)
            !%print *, "elemI(ii,ei_Lidx)", elemI(ii,ei_Lidx)
            exit
         end if
      end do

      !% --- Now for faces we do something similar but we have to change it abit 
      !%     because there are certain faces that are shared among images, which 
      !%     means their global index could be different.
      !%     This means we can't do the same thing to find the current length, 
      !%     instead we need to find the first face that is not a shared face and 
      !%     then calculate the current length based of that.
      ii = 1

      do while(faceYN(ii,fYN_isSharedFace))
         ii = ii + 1
      end do

      current_length = faceI(ii, ei_Gidx) - ii

      !% --- Now that we have the correct current length we do the same thing as 
      !%     before with the elems, except we have to check if it a shared_face or not.
      !%     So we write an extra if statement before checking to skip those faces.

      do ii = 1, size(faceI(:,ei_Gidx))
         if (faceYN(ii,fYN_isSharedFace)) then
            cycle

         else if (faceI(ii,ei_Gidx) /= current_length+faceI(ii,ei_Lidx) .and. faceI(ii,ei_Gidx)/= nullvalueI) then
            print *, "error in face global indexes. Processor :: ", this_image()
            !% print *, "faceI(ii,ei_Gidx)", faceI(ii,ei_Gidx)
            !% print *, "faceI(ii,ei_Lidx)", faceI(ii,ei_Lidx)
            exit
         end if
      end do

      !%------------------------------------------------------------------
      !% Closing
      if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_utest_global_index_check
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%   
   integer function util_utest_get_elemR_col (eIsNanCol) result(elemRCol)
      !%--------------------------------------------------------------------
      !% Description:
      !% gets the column in the elemR array corresponding to the column in
      !% the elemIsNan array
      !%--------------------------------------------------------------------
      !% Declarations
         integer, intent(in) :: eIsNanCol
      !%--------------------------------------------------------------------

      select case (eIsNanCol)
         case (    eIsNan_Area)
            elemRCol = er_Area
         case (    eIsNan_Depth)
            elemRCol = er_Depth
         case (    eIsNan_EllDepth)
            elemRCol = er_EllDepth
         case (    eIsNan_Flowrate)
            elemRCol = er_Flowrate
         case (    eIsNan_FlowrateLateral)
            elemRCol = er_FlowrateLateral
         case (    eIsNan_FroudeNumber)
            elemRCol = er_FroudeNumber
         case (    eIsNan_Head)
            elemRCol = er_Head
         case (    eIsNan_HydRadius)
            elemRCol = er_HydRadius
         case (    eIsNan_InterpWeight_uG)
            elemRCol = er_InterpWeight_uG
         case (    eIsNan_InterpWeight_dG)
            elemRCol = er_InterpWeight_dG
         case (    eIsNan_InterpWeight_uH)
            elemRCol = er_InterpWeight_uH
         case (    eIsNan_InterpWeight_dH)
            elemRCol = er_InterpWeight_dH
         case (    eIsNan_InterpWeight_uQ)
            elemRCol = er_InterpWeight_uQ
         case (    eIsNan_InterpWeight_dQ)
            elemRCol = er_InterpWeight_dQ
         case (    eIsNan_InterpWeight_uP)
            elemRCol = er_InterpWeight_uP
         case (    eIsNan_InterpWeight_dP)
            elemRCol = er_InterpWeight_dP
         case (    eIsNan_Perimeter)
            elemRCol = er_Perimeter
         case (    eIsNan_SlotDepth)
            elemRCol = er_SlotDepth
         case (    eIsNan_SlotArea)
            elemRCol = er_SlotArea
         case (    eIsNan_SlotVolume)
            elemRCol = er_SlotVolume
         case (    eIsNan_SourceContinuity)
            elemRCol = er_SourceContinuity
         case (    eIsNan_SourceMomentum)
            elemRCol = er_SourceMomentum
         case (    eIsNan_Velocity)
            elemRCol = er_Velocity
         case (    eIsNan_Volume)
            elemRCol = er_Volume
         case (    eIsNan_WaveSpeed)
            elemRCol = er_WaveSpeed
         case default
            print *, 'CODE ERROR: unexpected case value'
            call util_crashpoint(7288734)
      end select

   end function util_utest_get_elemR_col 
!%
!%==========================================================================
!%==========================================================================
!%   
   subroutine util_utest_get_elemR_indexName (eIsNanCol, colIndexName)
      !%--------------------------------------------------------------------
      !% Description:
      !% gets the name of the data in the column in the elemR array corresponding 
      !% to the column in the elemIsNan array
      !%--------------------------------------------------------------------
      !% Declarations
         integer, intent(in)           :: eIsNanCol
         character (64), intent(inout) :: colIndexName
      !%--------------------------------------------------------------------

      select case (eIsNanCol)
         case (      eIsNan_Area)
            colIndexName = 'Area'
         case (      eIsNan_Depth)
            colIndexName = 'Depth'
         case (      eIsNan_EllDepth)
            colIndexName = 'EllDepth'
         case (      eIsNan_Flowrate)
            colIndexName = 'Flowrate'
         case (      eIsNan_FlowrateLateral)
            colIndexName = 'FlowrateLateral'
         case (      eIsNan_FroudeNumber)
            colIndexName = 'FroudeNumber'
         case (      eIsNan_Head)
            colIndexName = 'Head'
         case (      eIsNan_HydRadius)
            colIndexName = 'HydRadius'
         case (      eIsNan_InterpWeight_uG)
            colIndexName = 'InterpWeight_uG'
         case (      eIsNan_InterpWeight_dG)
            colIndexName = 'InterpWeight_dG'
         case (      eIsNan_InterpWeight_uH)
            colIndexName = 'InterpWeight_uH'
         case (      eIsNan_InterpWeight_dH)
            colIndexName = 'InterpWeight_dH'
         case (      eIsNan_InterpWeight_uQ)
            colIndexName = 'InterpWeight_uQ'
         case (      eIsNan_InterpWeight_dQ)
            colIndexName = 'InterpWeight_dQ'
         case (      eIsNan_InterpWeight_uP)
            colIndexName = 'InterpWeight_uP'
         case (      eIsNan_InterpWeight_dP)
            colIndexName = 'InterpWeight_dP'
         case (      eIsNan_Perimeter)
            colIndexName = 'Perimeter'
         case (      eIsNan_SlotDepth)
            colIndexName = 'SlotDepth'
         case (      eIsNan_SlotArea)
            colIndexName = 'SlotArea'
         case (      eIsNan_SlotVolume)
            colIndexName = 'SlotVolume'
         case (      eIsNan_SourceContinuity)
            colIndexName = 'SourceContinuity'
         case (      eIsNan_SourceMomentum)
            colIndexName = 'SourceMomentum'
         case (      eIsNan_Velocity)
            colIndexName = 'Velocity'
         case (      eIsNan_Volume)
            colIndexName = 'Volume'
         case (      eIsNan_WaveSpeed)
            colIndexName = 'WaveSpeed'
         case default
            print *, 'CODE ERROR: unexpected case value'
            call util_crashpoint(7288734)
      end select

end subroutine util_utest_get_elemR_indexName 
!%
!%==========================================================================
!%==========================================================================
!%   
   integer function util_utest_get_faceR_col (fIsNanCol) result(faceRCol)
      !%--------------------------------------------------------------------
      !% Description:
      !% gets the column in the elemR array corresponding to the column in
      !% the elemIsNan array
      !%--------------------------------------------------------------------
      !% Declarations
         integer, intent(in) :: fIsNanCol
      !%--------------------------------------------------------------------

      select case (fIsNanCol)
         case (    fIsNan_Area_d)
            faceRCol = fr_Area_d
         case (    fIsNan_Area_u)
            faceRCol = fr_Area_u
         case (    fIsNan_Depth_d)
            faceRCol = fr_Depth_d
         case (    fIsNan_Depth_u)
            faceRCol = fr_Depth_u
         case (    fIsNan_Flowrate)
            faceRCol = fr_Flowrate
         case (    fIsNan_Flowrate_Conservative)
            faceRCol = fr_Flowrate_Conservative
         case (    fIsNan_Head_u)
            faceRCol = fr_Head_u
         case (    fIsNan_Head_d)
            faceRCol = fr_Head_d
         case (    fIsNan_Velocity_d)
            faceRCol = fr_Velocity_d
         case (    fIsNan_Velocity_u)
            faceRCol = fr_Velocity_u
         case (    fIsNan_Preissmann_Number)
            faceRCol = fr_Preissmann_Number
         case default
            print *, 'CODE ERROR: unexpected case value'
            call util_crashpoint(94023)
      end select

   end function util_utest_get_faceR_col 
!%
!%==========================================================================
!%==========================================================================
!%   
   subroutine util_utest_get_faceR_indexName (eIsNanCol, colIndexName)
      !%--------------------------------------------------------------------
      !% Description:
      !% gets the name of the data in the column in the elemR array corresponding 
      !% to the column in the elemIsNan array
      !%--------------------------------------------------------------------
      !% Declarations
         integer, intent(in)           :: eIsNanCol
         character (64), intent(inout) :: colIndexName
      !%--------------------------------------------------------------------

      select case (eIsNanCol)
         case (      fIsNan_Area_d)
            colIndexName = 'Area_d'
         case (      fIsNan_Area_u)
            colIndexName = 'Area_u'
         case (      fIsNan_Depth_d)
            colIndexName = 'Depth_d'
         case (      fIsNan_Depth_u)
            colIndexName = 'Depth_u'
         case (      fIsNan_Flowrate)
            colIndexName = 'Flowrate'
         case (      fIsNan_Flowrate_Conservative)
            colIndexName = 'Flowrate_Conservative'
         case (      fIsNan_Head_d)
            colIndexName = 'Head_d'
         case (      fIsNan_Head_u)
            colIndexName = 'Head_'
         case (      fIsNan_Velocity_d)
            colIndexName = 'Velocity_d'
         case (      fIsNan_Velocity_u)
            colIndexName = 'Velocity_u'
         case (      fIsNan_Preissmann_Number)
            colIndexName = 'Preissmann_Number'
         case default
            print *, 'CODE ERROR: unexpected case value'
            call util_crashpoint(6111837)
      end select

   end subroutine util_utest_get_faceR_indexName 
!%
!%==========================================================================  
!%==========================================================================
!%
   logical function util_utest_isThisCol_Nan (thisCol, isElem) result(outvalue)
      !%--------------------------------------------------------------------
      !% Description:
      !% tests if there are any NaN in elemR(:,thisCol) vector
      !%--------------------------------------------------------------------
      !% Declarations
         integer, intent(in) :: thisCol
         logical, intent(in) :: isElem ! .true. if element, .false. if face
      !%--------------------------------------------------------------------

      if (isElem) then 
         if (any(isnan(elemR(:,thisCol)))) then 
            outvalue = .true.
         else 
            outvalue = .false.
         end if
      else 
         if (any(isnan(faceR(:,thisCol)))) then 
            outvalue = .true.
         else 
            outvalue = .false.
         end if
      end if

   end function util_utest_isThisCol_Nan   
!%
!%==========================================================================
!% END MODULE    
!%==========================================================================
!%
  end module utility_unit_testing
