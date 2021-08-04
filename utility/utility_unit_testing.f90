Module utility_unit_testing
  use interface
  use utility_allocate
  use discretization
  use define_indexes
  use define_keys
  use define_globals
  use define_settings
  use network_define

  implicit none

  private

  public :: util_utest_local_global
  public :: util_utest_pack_arrays
  public :: util_utest_node_link_image
  public :: util_utest_slope_checking
  public :: util_utest_global_index_check

contains

  subroutine util_utest_local_global

    !% In this subroutine we are checking the the local and global indexs of the link,node,elem and face arrays are unique

    integer ii, jj, kk, min_val, max_val
    logical dup_found
    character(64) :: subroutine_name = 'local_global_unique'
    if (setting%Debug%File%initialization) print *, '*** enter ', subroutine_name

    !% Looping through the array and finding all of the unqiue values
    min_val = minval(link%I(:,li_idx)) - 1
    max_val = maxval(link%I(:,li_idx))
    ii = 0

    do while(min_val<max_val)
       ii = ii + 1
       min_val = minval(link%I(:,li_idx),mask=link%I(:,li_idx)>min_val)
    end do

    if(ii /= size(link%I(:, li_idx))) then
       print *, "ERROR:::: link%I(:,li_idx) is not unique. This_image ::", this_image()
    else
       print *, "link%I(:,li_idx) is unique. This_image ::", this_image()
    end if

    !% checking node%I(:,:) indexes

    min_val = minval(node%I(:,ni_idx)) - 1
    max_val = maxval(node%I(:,ni_idx))
    ii = 0

    do while(min_val<max_val)
       ii = ii + 1
       min_val = minval(node%I(:,ni_idx),mask=node%I(:,ni_idx)>min_val)
    end do



    if(ii /= size(node%I(:, ni_idx))) then
       print *, "ERROR:::: node%I(:,ni_idx) is not unique. This_image ::", this_image()
    else
       print *, "node%I(:,ni_idx) is unique. This_image ::", this_image()
    end if


    !% checking elemI(:,:) local indexes

    min_val = minval(elemI(1:N_elem(this_image()),ei_Lidx)) - 1
    max_val = maxval(elemI(1:N_elem(This_image()),ei_Lidx))
    ii = 0

    do while(min_val<max_val)
       ii = ii + 1
       min_val = minval(elemI(1:N_elem(This_image()),ei_Lidx),mask=elemI(1:N_elem(This_image()),ei_Lidx)>min_val)
    end do



    if(ii /= n_elem(this_image())) then
       print *, "ERROR:::: elemI(:,ei_Lidx) is not unique. This_image ::", this_image()
    else
       print *, "elemI(:,ei_Lidx) is unique. This_image ::", this_image()
    end if



    !% checking faceI(:,:) indexes

    min_val = minval(faceI(1:N_face(this_image()),fi_Lidx)) - 1
    max_val = maxval(faceI(1:N_face(this_image()),fi_Lidx))
    ii = 0

    do while(min_val<max_val)
       ii = ii + 1
       min_val = minval(faceI(1:N_face(this_image()),fi_Lidx),mask=faceI(1:N_face(this_image()),fi_Lidx)>min_val)
    end do



    if(ii /= n_face(this_image())) then
       print *, "ERROR:::: faceI(:,fi_Lidx) is not unique. This_image ::", this_image()
    else
       print *, "faceI(:,fi_Lidx) is unique. This_image ::", this_image()
    end if

    !% checking faceI(:,:) global indexes

    min_val = minval(faceI(1:N_face(this_image()),fi_Gidx)) - 1
    max_val = maxval(faceI(1:N_face(this_image()),fi_Gidx))
    ii = 0

    do while(min_val<max_val)
       ii = ii + 1
       min_val = minval(faceI(1:N_face(this_image()),fi_Gidx),mask=faceI(1:N_face(this_image()),fi_Gidx)>min_val)
    end do



    if(ii /= n_face(this_Image())) then
       print *, "ERROR:::: faceI(:,fi_Gidx) is not unique. This_image ::", this_image()
    else
       print *, "faceI(:,fi_Gidx) is unique. This_image ::", this_image()
    end if


    !% checking elemI(:,:) global indexes

    min_val = minval(elemI(1:N_elem(this_image()),ei_Gidx)) - 1
    max_val = maxval(elemI(1:N_elem(this_image()),ei_Gidx))
    ii = 0

    do while(min_val<max_val)
       ii = ii + 1
       min_val = minval(elemI(1:N_elem(this_image()),ei_Gidx),mask=elemI(1:N_elem(this_image()),ei_Gidx)>min_val)
    end do


    if(ii /= N_elem(this_image())) then
       print *, "ERROR:::: elemI(:,ei_Gidx) is not unique. This_image ::", this_image()
    else
       print *, "elemI(:,ei_Gidx) is unique. This_image ::", this_image()
    end if


    if (setting%Debug%File%initialization)  print *, '*** leave ', subroutine_name

  end subroutine util_utest_local_global

    subroutine util_utest_pack_arrays

      !% Going through all of the pack arrays and making sure they are unique, this follows the same process as the subroutine above, with a slight change.
      integer ii, jj, kk, min_val, max_val
      logical dup_found
      character(64) :: subroutine_name = 'pack_arrays_unique'
      if (setting%Debug%File%initialization) print *, '*** enter ', subroutine_name



      kk = 1
      dup_found = .false.

      !% checking elemP(:,ep_AC) indexes
      min_val = minval(elemP(:,ep_AC)) - 1
      max_val = maxval(elemP(:,ep_AC))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_AC),mask=elemP(:,ep_AC)>min_val)
      end do

      !% Here we check if nullvalueI is part packed array and if so we subtract one.
      !% We do this because
      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_AC) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_ac)) then
         print *, "ERROR:::: elemP(:,ep_AC) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_AC) is unique. This_image ::", this_image()
      end if


      !% checking elemP(:,ep_ALLtm) indexes

      min_val = minval(elemP(:,ep_ALLtm)) - 1
      max_val = maxval(elemP(:,ep_ALLtm))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_ALLtm),mask=elemP(:,ep_ALLtm)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if

      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_ALLtm) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_ALLtm)) then
         print *, "ERROR:::: elemP(:,ep_ALLtm) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_ALLtm) is unique. This_image ::", this_image()
      end if


      !% checking elemP(:,ep_CC_AC) indexes

      min_val = minval(elemP(:,ep_CC_AC)) - 1
      max_val = maxval(elemP(:,ep_CC_AC))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_CC_AC),mask=elemP(:,ep_CC_AC)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if

      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CC_AC) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_CC_AC)) then
         print *, "ERROR:::: elemP(:,ep_CC_AC) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CC_AC) is unique. This_image ::", this_image()
      end if


      !% checking elemP(:,ep_CC_ALLtm) indexes

      min_val = minval(elemP(:,ep_CC_ALLtm)) - 1
      max_val = maxval(elemP(:,ep_CC_ALLtm))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_CC_ALLtm),mask=elemP(:,ep_CC_ALLtm)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if




      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CC_ALLtm) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_CC_Alltm)) then
         print *, "ERROR:::: elemP(:,ep_CC_ALLtm) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CC_ALLtm) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_CC_ETM) indexes

      min_val = minval(elemP(:,ep_CC_ETM)) - 1
      max_val = maxval(elemP(:,ep_CC_ETM))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_CC_ETM),mask=elemP(:,ep_CC_ETM)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CC_ETM) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_CC_ETM)) then
         print *, "ERROR:::: elemP(:,ep_CC_ETM) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CC_ETM) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_CC_H_ETM) indexes


      min_val = minval(elemP(:,ep_CC_H_ETM)) - 1
      max_val = maxval(elemP(:,ep_CC_H_ETM))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_CC_H_ETM),mask=elemP(:,ep_CC_H_ETM)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CC_H_ETM) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_CC_H_ETM)) then
         print *, "ERROR:::: elemP(:,ep_CC_H_ETM) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CC_H_ETM) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_CC_Q_AC) indexes


      min_val = minval(elemP(:,ep_CC_Q_AC)) - 1
      max_val = maxval(elemP(:,ep_CC_Q_AC))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_CC_Q_AC),mask=elemP(:,ep_CC_Q_AC)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CC_Q_AC) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_CC_Q_AC)) then
         print *, "ERROR:::: elemP(:,ep_CC_Q_AC) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CC_Q_AC) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_CC_Q_ETM) indexes

      min_val = minval(elemP(:,ep_CC_Q_ETM)) - 1
      max_val = maxval(elemP(:,ep_CC_Q_ETM))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_CC_Q_ETM),mask=elemP(:,ep_CC_Q_ETM)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CC_Q_ETM) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_CC_Q_ETM)) then
         print *, "ERROR:::: elemP(:,ep_CC_Q_ETM) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CC_Q_ETM) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_CCJB_AC_surcharged) indexes


      min_val = minval(elemP(:,ep_CCJB_AC_surcharged)) - 1
      max_val = maxval(elemP(:,ep_CCJB_AC_surcharged))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_CCJB_AC_surcharged),mask=elemP(:,ep_CCJB_AC_surcharged)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CCJB_AC_surcharged) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_CCJB_AC_surcharged)) then
         print *, "ERROR:::: elemP(:,ep_CCJB_AC_surcharged) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CCJB_AC_surcharged) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_CCJB_ALLtm) indexes


      min_val = minval(elemP(:,ep_CCJB_ALLtm)) - 1
      max_val = maxval(elemP(:,ep_CCJB_ALLtm))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_CCJB_ALLtm),mask=elemP(:,ep_CCJB_ALLtm)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CCJB_ALLtm) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_CCJB_ALLtm)) then
         print *, "ERROR:::: elemP(:,ep_CCJB_ALLtm) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CCJB_ALLtm) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_CCJB_AC) indexes


      min_val = minval(elemP(:,ep_CCJB_AC)) - 1
      max_val = maxval(elemP(:,ep_CCJB_AC))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_CCJB_AC),mask=elemP(:,ep_CCJB_AC)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CCJB_AC) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_CCJB_AC)) then
         print *, "ERROR:::: elemP(:,ep_CCJB_AC) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CCJB_AC) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_CCJB_ALLtm_surcharged) indexes

      min_val = minval(elemP(:,ep_CCJB_ALLtm_surcharged)) - 1
      max_val = maxval(elemP(:,ep_CCJB_ALLtm_surcharged))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_CCJB_ALLtm_surcharged),mask=elemP(:,ep_CCJB_ALLtm_surcharged)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CCJB_ALLtm_surcharged) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_CCJB_ALLtm_surcharged)) then
         print *, "ERROR:::: elemP(:,ep_CCJB_ALLtm_surcharged) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CCJB_ALLtm_surcharged) is unique. This_image ::", this_image()
      end if



      !% checking elemP(:,ep_CCJB_eETM_i_fAC) indexes


      min_val = minval(elemP(:,ep_CCJB_eETM_i_fAC)) - 1
      max_val = maxval(elemP(:,ep_CCJB_eETM_i_fAC))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_CCJB_eETM_i_fAC),mask=elemP(:,ep_CCJB_eETM_i_fAC)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CCJB_eETM_i_fAC) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= N_elem(this_image())) then
         print *, "ERROR:::: elemP(:,ep_CCJB_eETM_i_fAC) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CCJB_eETM_i_fAC) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_CCJB_ETM) indexes


      min_val = minval(elemP(:,ep_CCJB_ETM)) - 1
      max_val = maxval(elemP(:,ep_CCJB_ETM))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_CCJB_ETM),mask=elemP(:,ep_CCJB_ETM)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CCJB_ETM) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_CCJB_ETM)) then
         print *, "ERROR:::: elemP(:,ep_CCJB_ETM) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CCJB_ETM) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_Diag) indexes


      min_val = minval(elemP(:,ep_Diag)) - 1
      max_val = maxval(elemP(:,ep_Diag))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_Diag),mask=elemP(:,ep_Diag)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_Diag) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_Diag)) then
         print *, "ERROR:::: elemP(:,ep_Diag) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_Diag) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_ETM) indexes

      min_val = minval(elemP(:,ep_ETM)) - 1
      max_val = maxval(elemP(:,ep_ETM))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_ETM),mask=elemP(:,ep_ETM)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_ETM) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_ETM)) then
         print *, "ERROR:::: elemP(:,ep_ETM) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_ETM) is unique. This_image ::", this_image()
      end if


      !% checking elemP(:,ep_JM_AC) indexes

      min_val = minval(elemP(:,ep_JM_AC)) - 1
      max_val = maxval(elemP(:,ep_JM_AC))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_JM_AC),mask=elemP(:,ep_JM_AC)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_JM_AC) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_JM_AC)) then
         print *, "ERROR:::: elemP(:,ep_JM_AC) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_JM_AC) is unique. This_image ::", this_image()
      end if


      !% checking elemP(:,ep_JM_ALLtm) indexes


      min_val = minval(elemP(:,ep_JM_ALLtm)) - 1
      max_val = maxval(elemP(:,ep_JM_ALLtm))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_JM_ALLtm),mask=elemP(:,ep_JM_ALLtm)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_JM_ALLtm) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_JM_ALLtm)) then
         print *, "ERROR:::: elemP(:,ep_JM_ALLtm) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_JM_ALLtm) is unique. This_image ::", this_image()
      end if


      !% checking elemP(:,ep_JB_ETM) indexes

      min_val = minval(elemP(:,ep_JB_ETM)) - 1
      max_val = maxval(elemP(:,ep_JB_ETM))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_JB_ETM),mask=elemP(:,ep_JB_ETM)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_JB_ETM) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(EP_JB_ETM)) then
         print *, "ERROR:::: elemP(:,ep_JB_ETM) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_JB_ETM) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_NonSurcharged_AC) indexes


      min_val = minval(elemP(:,ep_NonSurcharged_AC)) - 1
      max_val = maxval(elemP(:,ep_NonSurcharged_AC))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_NonSurcharged_AC),mask=elemP(:,ep_NonSurcharged_AC)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_NonSurcharged_AC) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_NonSurcharged_AC)) then
         print *, "ERROR:::: elemP(:,ep_NonSurcharged_AC) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_NonSurcharged_AC) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_NonSurcharged_ALLtm) indexes

      min_val = minval(elemP(:,ep_NonSurcharged_ALLtm)) - 1
      max_val = maxval(elemP(:,ep_NonSurcharged_ALLtm))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_NonSurcharged_ALLtm),mask=elemP(:,ep_NonSurcharged_ALLtm)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_NonSurcharged_ALLtm) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_NonSurcharged_Alltm)) then
         print *, "ERROR:::: elemP(:,ep_NonSurcharged_ALLtm) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_NonSurcharged_ALLtm) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_NonSurcharged_ETM) indexes


      min_val = minval(elemP(:,ep_NonSurcharged_ETM)) - 1
      max_val = maxval(elemP(:,ep_NonSurcharged_ETM))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_NonSurcharged_ETM),mask=elemP(:,ep_NonSurcharged_ETM)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_NonSurcharged_ETM) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_NonSurcharged_ETM)) then
         print *, "ERROR:::: elemP(:,ep_NonSurcharged_ETM) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_NonSurcharged_ETM) is unique. This_image ::", this_image()
      end if


      !% checking elemP(:,ep_smallvolume_AC) indexes


      min_val = minval(elemP(:,ep_smallvolume_AC)) - 1
      max_val = maxval(elemP(:,ep_smallvolume_AC))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_smallvolume_AC),mask=elemP(:,ep_smallvolume_AC)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_smallvolume_AC) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_smallvolume_AC)) then
         print *, "ERROR:::: elemP(:,ep_smallvolume_AC) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_smallvolume_AC) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_smallvolume_ALLtm) indexes


      min_val = minval(elemP(:,ep_smallvolume_ALLtm)) - 1
      max_val = maxval(elemP(:,ep_smallvolume_ALLtm))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_smallvolume_ALLtm),mask=elemP(:,ep_smallvolume_ALLtm)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_smallvolume_ALLtm) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_smallvolume_Alltm)) then
         print *, "ERROR:::: elemP(:,ep_smallvolume_ALLtm) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_smallvolume_ALLtm) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_smallvolume_ETM) indexes

      min_val = minval(elemP(:,ep_smallvolume_ETM)) - 1
      max_val = maxval(elemP(:,ep_smallvolume_ETM))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_smallvolume_ETM),mask=elemP(:,ep_smallvolume_ETM)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_smallvolume_ETM) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_smallvolume_ETM)) then
         print *, "ERROR:::: elemP(:,ep_smallvolume_ETM) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_smallvolume_ETM) is unique. This_image ::", this_image()
      end if


      !% checking elemP(:,ep_Surcharged_AC) indexes


      min_val = minval(elemP(:,ep_Surcharged_AC)) - 1
      max_val = maxval(elemP(:,ep_Surcharged_AC))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_Surcharged_AC),mask=elemP(:,ep_Surcharged_AC)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_Surcharged_AC) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_Surcharged_AC)) then
         print *, "ERROR:::: elemP(:,ep_Surcharged_AC) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_Surcharged_AC) is unique. This_image ::", this_image()
      end if


      !% checking elemP(:,ep_Surcharged_ALLtm) indexes

      min_val = minval(elemP(:,ep_Surcharged_ALLtm)) - 1
      max_val = maxval(elemP(:,ep_Surcharged_ALLtm))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_Surcharged_ALLtm),mask=elemP(:,ep_Surcharged_ALLtm)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_Surcharged_ALLtm) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= N_elem(this_image())) then
         print *, "ERROR:::: elemP(:,ep_Surcharged_ALLtm) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_Surcharged_ALLtm) is unique. This_image ::", this_image()
      end if


      !% checking elemP(:,ep_Surcharged_ETM) indexes


      min_val = minval(elemP(:,ep_Surcharged_ETM)) - 1
      max_val = maxval(elemP(:,ep_Surcharged_ETM))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_Surcharged_ETM),mask=elemP(:,ep_Surcharged_ETM)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_Surcharged_ETM) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_Surcharged_ETM)) then
         print *, "ERROR:::: elemP(:,ep_Surcharged_ETM) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_Surcharged_ETM) is unique. This_image ::", this_image()
      end if


      !% checking elemP(:,ep_CCJM_H_AC_surcharged) indexes

      min_val = minval(elemP(:,ep_CCJM_H_AC_surcharged)) - 1
      max_val = maxval(elemP(:,ep_CCJM_H_AC_surcharged))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_CCJM_H_AC_surcharged),mask=elemP(:,ep_CCJM_H_AC_surcharged)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if

      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CCJM_H_AC_surcharged) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_CCJM_H_AC_surcharged)) then
         print *, "ERROR:::: elemP(:,ep_CCJM_H_AC_surcharged) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CCJM_H_AC_surcharged) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_CCJM_H_AC) indexes

      min_val = minval(elemP(:,ep_CCJM_H_AC)) - 1
      max_val = maxval(elemP(:,ep_CCJM_H_AC))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_CCJM_H_AC),mask=elemP(:,ep_CCJM_H_AC)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CCJM_H_AC) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemP(ep_CCJM_H_AC)) then
         print *, "ERROR:::: elemP(:,ep_CCJM_H_AC) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CCJM_H_AC) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_CCJB_eAC_i_fETM) indexes

      min_val = minval(elemP(:,ep_CCJB_eAC_i_fETM)) - 1
      max_val = maxval(elemP(:,ep_CCJB_eAC_i_fETM))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_CCJB_eAC_i_fETM),mask=elemP(:,ep_CCJB_eAC_i_fETM)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CCJB_eAC_i_fETM) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= N_elem(this_image())) then
         print *, "ERROR:::: elemP(:,ep_CCJB_eAC_i_fETM) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CCJB_eAC_i_fETM) is unique. This_image ::", this_image()
      end if

      !% checking elemPGalltm(:,epg_CCJM_rectangular_nonsurcharged) indexes


      min_val = minval(elemPGalltm(:,epg_CCJM_rectangular_nonsurcharged)) - 1
      max_val = maxval(elemPGalltm(:,epg_CCJM_rectangular_nonsurcharged))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemPGalltm(:,epg_CCJM_rectangular_nonsurcharged),&
              mask=elemPGalltm(:,epg_CCJM_rectangular_nonsurcharged)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemPGalltm(:,epg_CCJM_rectangular_nonsurcharged) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemPGalltm(epg_CCJM_rectangular_nonsurcharged)) then
         print *, "ERROR:::: elemPGalltm(:,epg_CCJM_rectangular_nonsurcharged) is not unique. This_image ::", this_image()

      else
         print *, "elemPGalltm(:,epg_CCJM_rectangular_nonsurcharged) is unique. This_image ::", this_image()
      end if

      !% checking elemPGalltm(:,epg_CCJM_trapezoidal_nonsurcharged) indexes

      min_val = minval(elemPGalltm(:,epg_CCJM_trapezoidal_nonsurcharged)) - 1
      max_val = maxval(elemPGalltm(:,epg_CCJM_trapezoidal_nonsurcharged))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemPGalltm(:,epg_CCJM_trapezoidal_nonsurcharged),&
              mask=elemPGalltm(:,epg_CCJM_trapezoidal_nonsurcharged)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemPGalltm(:,epg_CCJM_trapezoidal_nonsurcharged) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemPGalltm(epg_CCJM_trapezoidal_nonsurcharged)) then
         print *, "ERROR:::: elemPGalltm(:,epg_CCJM_trapezoidal_nonsurcharged) is not unique. This_image ::", this_image()

      else
         print *, "elemPGalltm(:,epg_CCJM_trapezoidal_nonsurcharged) is unique. This_image ::", this_image()
      end if


      !% checking elemPGalltm(:,epg_JB_rectangular) indexes


      min_val = minval(elemPGalltm(:,epg_JB_rectangular)) - 1
      max_val = maxval(elemPGalltm(:,epg_JB_rectangular))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemPGalltm(:,epg_JB_rectangular),mask=elemPGalltm(:,epg_JB_rectangular)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemPGalltm(:,epg_JB_rectangular) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemPGalltm(epg_JB_rectangular)) then
         print *, "ERROR:::: elemPGalltm(:,epg_JB_rectangular) is not unique. This_image ::", this_image()

      else
         print *, "elemPGalltm(:,epg_JB_rectangular) is unique. This_image ::", this_image()
      end if


      !% checking elemPGalltm(:,epg_JB_trapezoidal) indexes


      min_val = minval(elemPGalltm(:,epg_JB_trapezoidal)) - 1
      max_val = maxval(elemPGalltm(:,epg_JB_trapezoidal))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemPGalltm(:,epg_JB_trapezoidal),mask=elemPGalltm(:,epg_JB_trapezoidal)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemPGalltm(:,epg_JB_trapezoidal) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemPGalltm(epg_JB_trapezoidal)) then
         print *, "ERROR:::: elemPGalltm(:,epg_JB_trapezoidal) is not unique. This_image ::", this_image()

      else
         print *, "elemPGalltm(:,epg_JB_trapezoidal) is unique. This_image ::", this_image()
      end if


      !% --------------------------------------------------------------------------------------------------------------

      !% checking elemPGetm(:,epg_CCJM_rectangular_nonsurcharged) indexes


      min_val = minval(elemPGetm(:,epg_CCJM_rectangular_nonsurcharged)) - 1
      max_val = maxval(elemPGetm(:,epg_CCJM_rectangular_nonsurcharged))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemPGetm(:,epg_CCJM_rectangular_nonsurcharged),&
              mask=elemPGetm(:,epg_CCJM_rectangular_nonsurcharged)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if

      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemPGetm(:,epg_CCJM_rectangular_nonsurcharged) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemPGetm(epg_CCJM_rectangular_nonsurcharged)) then
         print *, "ERROR:::: elemPGetm(:,epg_CCJM_rectangular_nonsurcharged) is not unique. This_image ::", this_image()

      else
         print *, "elemPGetm(:,epg_CCJM_rectangular_nonsurcharged) is unique. This_image ::", this_image()
      end if

      !% checking elemPGetm(:,epg_CCJM_trapezoidal_nonsurcharged) indexes

      min_val = minval(elemPGetm(:,epg_CCJM_trapezoidal_nonsurcharged)) - 1
      max_val = maxval(elemPGetm(:,epg_CCJM_trapezoidal_nonsurcharged))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemPGetm(:,epg_CCJM_trapezoidal_nonsurcharged), &
              mask=elemPGetm(:,epg_CCJM_trapezoidal_nonsurcharged)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if

      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemPGetm(:,epg_CCJM_trapezoidal_nonsurcharged) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemPGetm(epg_CCJM_trapezoidal_nonsurcharged)) then
         print *, "ERROR:::: elemPGetm(:,epg_CCJM_trapezoidal_nonsurcharged) is not unique. This_image ::", this_image()

      else
         print *, "elemPGetm(:,epg_CCJM_trapezoidal_nonsurcharged) is unique. This_image ::", this_image()
      end if


      !% checking elemPGetm(:,epg_JB_rectangular) indexes


      min_val = minval(elemPGetm(:,epg_JB_rectangular)) - 1
      max_val = maxval(elemPGetm(:,epg_JB_rectangular))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemPGetm(:,epg_JB_rectangular),mask=elemPGetm(:,epg_JB_rectangular)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if

      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemPGetm(:,epg_JB_rectangular) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemPGetm(epg_JB_rectangular)) then
         print *, "ERROR:::: elemPGetm(:,epg_JB_rectangular) is not unique. This_image ::", this_image()

      else
         print *, "elemPGetm(:,epg_JB_rectangular) is unique. This_image ::", this_image()
      end if


      !% checking elemPGetm(:,epg_JB_trapezoidal) indexes


      min_val = minval(elemPGetm(:,epg_JB_trapezoidal)) - 1
      max_val = maxval(elemPGetm(:,epg_JB_trapezoidal))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemPGetm(:,epg_JB_trapezoidal),mask=elemPGetm(:,epg_JB_trapezoidal)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemPGetm(:,epg_JB_trapezoidal) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemPGetm(epg_JB_trapezoidal)) then
         print *, "ERROR:::: elemPGetm(:,epg_JB_trapezoidal) is not unique. This_image ::", this_image()

      else
         print *, "elemPGetm(:,epg_JB_trapezoidal) is unique. This_image ::", this_image()
      end if


      !% --------------------------------------------------------------------------------------------------------------

      !% checking elemPGac(:,epg_CCJM_rectangular_nonsurcharged) indexes


      min_val = minval(elemPGac(:,epg_CCJM_rectangular_nonsurcharged)) - 1
      max_val = maxval(elemPGac(:,epg_CCJM_rectangular_nonsurcharged))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemPGac(:,epg_CCJM_rectangular_nonsurcharged),&
              mask=elemPGac(:,epg_CCJM_rectangular_nonsurcharged)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemPGac(:,epg_CCJM_rectangular_nonsurcharged) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemPGac(epg_CCJM_rectangular_nonsurcharged)) then
         print *, "ERROR:::: elemPGac(:,epg_CCJM_rectangular_nonsurcharged) is not unique. This_image ::", this_image()

      else
         print *, "elemPGac(:,epg_CCJM_rectangular_nonsurcharged) is unique. This_image ::", this_image()
      end if

      !% checking elemPGac(:,epg_CCJM_trapezoidal_nonsurcharged) indexes

      min_val = minval(elemPGac(:,epg_CCJM_trapezoidal_nonsurcharged)) - 1
      max_val = maxval(elemPGac(:,epg_CCJM_trapezoidal_nonsurcharged))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemPGac(:,epg_CCJM_trapezoidal_nonsurcharged), &
              mask=elemPGac(:,epg_CCJM_trapezoidal_nonsurcharged)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemPGac(:,epg_CCJM_trapezoidal_nonsurcharged) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemPGac(epg_CCJM_trapezoidal_nonsurcharged)) then
         print *, "ERROR:::: elemPGac(:,epg_CCJM_trapezoidal_nonsurcharged) is not unique. This_image ::", this_image()

      else
         print *, "elemPGac(:,epg_CCJM_trapezoidal_nonsurcharged) is unique. This_image ::", this_image()
      end if


      !% checking elemPGac(:,epg_JB_rectangular) indexes


      min_val = minval(elemPGac(:,epg_JB_rectangular)) - 1
      max_val = maxval(elemPGac(:,epg_JB_rectangular))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemPGac(:,epg_JB_rectangular),mask=elemPGac(:,epg_JB_rectangular)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemPGac(:,epg_JB_rectangular) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemPGac(epg_JB_rectangular)) then
         print *, "ERROR:::: elemPGac(:,epg_JB_rectangular) is not unique. This_image ::", this_image()

      else
         print *, "elemPGac(:,epg_JB_rectangular) is unique. This_image ::", this_image()
      end if


      !% checking elemPGac(:,epg_JB_trapezoidal) indexes


      min_val = minval(elemPGac(:,epg_JB_trapezoidal)) - 1
      max_val = maxval(elemPGac(:,epg_JB_trapezoidal))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemPGac(:,epg_JB_trapezoidal),mask=elemPGac(:,epg_JB_trapezoidal)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "elemPGac(:,epg_JB_trapezoidal) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_elemPGac(epg_JB_trapezoidal)) then
         print *, "ERROR:::: elemPGac(:,epg_JB_trapezoidal) is not unique. This_image ::", this_image()

      else
         print *, "elemPGac(:,epg_JB_trapezoidal) is unique. This_image ::", this_image()
      end if

      !% checking faceP(:,fp_all) indexes


      min_val = minval(faceP(:,fp_all)) - 1
      max_val = maxval(faceP(:,fp_all))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(faceP(:,fp_all),mask=faceP(:,fp_all)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "faceP(:,fp_all) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_faceP(fp_all)) then
         print *, "ERROR:::: faceP(:,fp_all) is not unique. This_image ::", this_image()

      else
         print *, "faceP(:,fp_all) is unique. This_image ::", this_image()
      end if

      !% checking faceP(:,fp_AC) indexes


      min_val = minval(faceP(:,fp_AC)) - 1
      max_val = maxval(faceP(:,fp_AC))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(faceP(:,fp_AC),mask=faceP(:,fp_AC)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "faceP(:,fp_AC) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_faceP(fp_AC)) then
         print *, "ERROR:::: faceP(:,fp_AC) is not unique. This_image ::", this_image()

      else
         print *, "faceP(:,fp_AC) is unique. This_image ::", this_image()
      end if

      !% checking faceP(:,fp_Diag) indexes


      min_val = minval(faceP(:,fp_Diag)) - 1
      max_val = maxval(faceP(:,fp_Diag))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(faceP(:,fp_Diag),mask=faceP(:,fp_Diag)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "faceP(:,fp_Diag) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_faceP(fp_Diag)) then
         print *, "ERROR:::: faceP(:,fp_Diag) is not unique. This_image ::", this_image()

      else
         print *, "faceP(:,fp_Diag) is unique. This_image ::", this_image()
      end if

      !% checking faceP(:,fp_JumpDn) indexes


      min_val = minval(faceP(:,fp_JumpDn)) - 1
      max_val = maxval(faceP(:,fp_JumpDn))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(faceP(:,fp_JumpDn),mask=faceP(:,fp_JumpDn)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "faceP(:,fp_JumpDn) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_faceP(fp_JumpDn)) then
         print *, "ERROR:::: faceP(:,fp_JumpDn) is not unique. This_image ::", this_image()

      else
         print *, "faceP(:,fp_JumpDn) is unique. This_image ::", this_image()
      end if

      !% checking faceP(:,fp_JumpUp) indexes


      min_val = minval(faceP(:,fp_JumpUp)) - 1
      max_val = maxval(faceP(:,fp_JumpUp))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(faceP(:,fp_JumpUp),mask=faceP(:,fp_JumpUp)>min_val)
      end do

      if(min_val == nullvalueI) then
         ii = ii - 1
      end if


      if(ii == 0 .and. min_val == nullvalueI) then
         print *, "faceP(:,fp_JumpUp) is only filled with nullvalueI. This_image ::", this_image()

      else if(ii /= npack_faceP(fp_JumpUp)) then
         print *, "ERROR:::: faceP(:,fp_JumpUp) is not unique. This_image ::", this_image()

      else
         print *, "faceP(:,fp_JumpUp) is unique. This_image ::", this_image()
      end if


      if (setting%Debug%File%initialization)  print *, '*** leave ', subroutine_name

    end subroutine util_utest_pack_arrays

    subroutine util_utest_node_link_image

      !% In this subroutine we are checking whether the every image has atleast one node and link assigned to it.
      integer :: ii, jj, kk, counter
      character(64) :: subroutine_name = 'init_face_check'
      if (setting%Debug%File%initialization) print *, '*** enter ', subroutine_name

      kk = 1
      counter = 0

      !% We can find if there is a node on each image by counting the amount of unique values there is inside of node%I(:,ni_P_image)
      !% So we use the code below to loop through node%I(:,ni_P_image) and find all the unique values

      do ii = 1, size(node%I(:,ni_P_image))
         do jj = 1, ii

            if((node%I(ii, ni_P_image)) == node%I(jj, ni_P_image)) then
               exit
            end if
         end do

            if(ii == jj) then
               counter = counter + 1

            end if

      end do

      !% After that we compare the number of images we are using to what we counted.
      !% If it is correct they should be the same value, otherwise there was an error when paritioning the links and nodes

      if(num_images() /= counter) then
         print *, "error in NodeI images. This_image :: ", this_image()
      else
         print *, "correct number in nodeI images. This_image :: ", this_image()
      end if

      !% We reset the counter and do the same process for the Links
      counter = 0

      do ii = 1, size(link%I(:,li_P_image))
         do jj = 1, ii

            if((link%I(ii, li_P_image)) == link%I(jj, li_P_image)) then
               exit
            end if
         end do

            if(ii == jj) then
               counter = counter + 1

            end if

      end do


      if(num_images() /= counter) then
         print *, "error in linkI images. This_image :: ", this_image()
         print *, "counter", counter
      else
         print *, "correct number in linkI images.  This_image :: ", this_image()
      end if
      if (setting%Debug%File%initialization)  print *, '*** leave ', subroutine_name

    end subroutine util_utest_node_link_image

    subroutine util_utest_slope_checking
      !% In this subroutine we are checking that all of the slopes are postive.
      !% To do this and loop through and if we find a negative slope then we exit the loop and report it.

      integer :: ii, jj
      logical :: invalid_slope
      character(64) :: subroutine_name = 'slope_checking'
      if (setting%Debug%File%initialization) print *, '*** enter ', subroutine_name

      invalid_slope = .false.
      do ii = 1, size(link%R(:,lr_Slope))

         if(link%R(ii, lr_slope) < 0) then
            invalid_slope = .true.
            exit
         end if
      end do

      if(invalid_slope .eqv. .true.) then
         print *, "error found in link%R(:,lr_slope) slope is negative. This_image :: ", this_image()
      else
         print *, "all slopes are postive.  This_image :: ", this_image()
      end if

      if (setting%Debug%File%initialization)  print *, '*** leave ', subroutine_name

    end subroutine util_utest_slope_checking


    subroutine util_utest_global_index_check
      !% In this subroutine we are checking that the global indexs are correct by adding the first valid global index of an image with the local index

      integer :: ii, current_length, counter
      character(64) :: subroutine_name = 'global_index_checking'

      if (setting%Debug%File%initialization) print *, '*** enter ', subroutine_name

      !% here we find the current length of the global index by looking at the first value of elemI on that image and subtracting one.

      current_length = elemI(1, ei_Gidx) - 1

      !% Then we loop through elemI(:,ei_Gidx) and report and stop looping if there is an error in the global indexes
      do ii = 1, size(elemI(:,ei_Gidx))

         if(elemI(ii,ei_Gidx) /= current_length+elemI(ii,ei_Lidx) .and. elemI(ii,ei_Gidx)/= nullvalueI) then
            print *, "error in elem global indexes. Image :: ", this_Image()
            !%print *, "elemI(ii,ei_Gidx)", elemI(ii,ei_Gidx)
            !%print *, "elemI(ii,ei_Lidx)", elemI(ii,ei_Lidx)
            exit
         end if
      end do


      !% Now for faces we do something similar but we have to change it abit because there are certain faces that are shared among images which means their global index could be different.
      !% This means we can't do the same thing to find the current length, instead we need to find the first face that is not a shared face and then calculate the current length based of that.
      !% That is what is happening the lines below.
      ii = 1

      do while(faceYN(ii,fYN_isSharedFace))

         ii = ii + 1

      end do

      current_length = faceI(ii, ei_Gidx) - ii

      !% Now that we have the correct current length we do the same thing as before with the elems, except we have to check if it a shared_face or not.
      !% So we write an extra if statement before checking to skip those faces.

      do ii = 1, size(faceI(:,ei_Gidx))

         if(faceYN(ii,fYN_isSharedFace)) then
            cycle

         else if(faceI(ii,ei_Gidx) /= current_length+faceI(ii,ei_Lidx) .and. faceI(ii,ei_Gidx)/= nullvalueI) then
            print *, "error in face global indexes. Image :: ", this_Image()
            !% print *, "faceI(ii,ei_Gidx)", faceI(ii,ei_Gidx)
            !% print *, "faceI(ii,ei_Lidx)", faceI(ii,ei_Lidx)
            exit
         end if
      end do

      if (setting%Debug%File%initialization)  print *, '*** leave ', subroutine_name

    end subroutine util_utest_global_index_check

    !subroutine geometry_checking

      !integer ii, jj

      !do ii = 1, size(link%I(:,linktype))

         !if(link%I(ii, linktype) .eq. lchannel) then

          !  select case(link%I(ii,li_geometry)

     !          case(lRectangular)

    !end subroutine geometry_checking



  end Module utility_unit_testing
