module utility_unit_testing
  use interface_
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
  public :: util_utest_CLprint

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine util_utest_CLprint (inputstring)
        !%------------------------------------------------------------------
        !% Description:
        !% Use for command-line write during debugging
        !%------------------------------------------------------------------
        !% Declarations:
            character (len=*), intent(in) :: inputstring
            integer :: ii, jj, kk, nn, mm, fD, eD, thisCol, faceCol
            integer, pointer :: fup(:), fdn(:), eup(:), edn(:)
            integer, pointer :: thisP(:), Npack
            real(8), pointer :: dt, oneVec(:), grav
            real(8) :: hr, aa, bb

            ! integer   :: iet(4) =    (/49, 50, 51, 52/)
            ! integer   :: ift(3) = (/46, 47, 48/)
    
           integer :: iet(8) =  (/ 49,    51, 50, 52,     61,     63, 62,  64     /)
           integer :: ift(4)  = (/    42 ,            43,     52     ,         53    /)

            

            real(8) :: Qextra(4)
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
      
            return
            
        print *, ' '
        write(*,"(A,A, e12.5)") ' ',trim(inputstring)     
        write(*,"(A,i7,A, f12.5, A, f12.5, A)") '        step = ',setting%Time%Step ,&
        '; dt = ',setting%Time%Hydraulics%Dt,&
        '; time = ',setting%Time%Now, ' sec'

      !   do ii=1,N_elem(1)
      !     if ((elemI(ii,ei_elementType) == CC) .or. &
      !         (elemI(ii,ei_elementType) == JM) .or. &
      !         ((elemI(ii,ei_elementType) == JB) .and. (elemSI(ii,esi_JunctionBranch_Exists) == oneI))) then

      !          print *, ii, elemR(ii,er_Head), elemR(ii,er_Flowrate), trim(reverseKey(elemI(ii,ei_elementType)))

      !          !    write(*,"(i4,L,8f12.4)") ii, elemYN(ii,eYN_isPSsurcharged), elemR(ii,er_Head), elemR(ii,er_Depth), elemR(ii,er_SlotDepth), elemR(ii,er_SlotVolume), elemR(ii,er_Volume)
      !         end if
   
      !   end do

      ! do ii=1,N_elem(1)
      !    if ((elemI(ii,ei_elementType) == CC) .or. &
      !       ((elemI(ii,ei_elementType) == JB) .and. (elemSI(ii,esi_JunctionBranch_Exists) == oneI))) then
      !          !print *, ii, elemI(ii,ei_Mface_uL), elemI(ii,ei_Mface_dL)
      !          if (elemI(ii,ei_Mface_uL) .ne. nullvalueI) then 
      !           !  print *, ii, faceR(elemI(ii,ei_Mface_uL),fr_Head_u), faceR(elemI(ii,ei_Mface_uL),fr_Head_d),' for F ', elemI(ii,ei_Mface_uL)
      !             print *, ii, faceR(elemI(ii,ei_Mface_uL),fr_Flowrate), ' for F ', elemI(ii,ei_Mface_uL)
      !          end if
      !          !print *, ii, elemR(ii,er_Head), 'for C'
      !          print *, ii, elemR(ii,er_Flowrate), elemR(ii,er_Depth), 'for C'
      !          if (elemI(ii,ei_Mface_dL) .ne. nullvalueI) then 
      !            ! print *, ii, faceR(elemI(ii,ei_Mface_dL),fr_Head_u), faceR(elemI(ii,ei_Mface_dL),fr_Head_d),' for F ', elemI(ii,ei_Mface_dL)
      !             print *, ii, faceR(elemI(ii,ei_Mface_dL),fr_Flowrate), ' for F ', elemI(ii,ei_Mface_dL)
      !          end if
      !    end if
      ! end do


      ! do ii=1,N_elem(1)
      !    if ((elemI(ii,ei_elementType) == CC) .or. &
      !       ((elemI(ii,ei_elementType) == JB) .and. (elemSI(ii,esi_JunctionBranch_Exists) == oneI))) then
      !          if (elemI(ii,ei_Mface_uL) .ne. nullvalueI) then 
      !            print *, ii, faceR(elemI(ii,ei_Mface_uL),fr_Head_u), faceR(elemI(ii,ei_Mface_uL),fr_Head_d),' for F ', elemI(ii,ei_Mface_uL)
   
      !          end if
      !          print *, ii, elemR(ii,er_Head), 'for C'
              
      !          if (elemI(ii,ei_Mface_dL) .ne. nullvalueI) then 
      !            print *, ii, faceR(elemI(ii,ei_Mface_dL),fr_Head_u), faceR(elemI(ii,ei_Mface_dL),fr_Head_d),' for F ', elemI(ii,ei_Mface_dL)
      !          end if
      !    end if
      ! end do

      do ii=1,N_elem(1)
         if (elemI(ii,ei_elementType) == CC) then

            if (elemI(ii,ei_Mface_uL) .ne. nullvalueI) then 
               write(*,"(i3,i3, 4(f12.4))") ii, &
               faceI(elemI(ii,ei_Mface_uL),fi_zeroDepth), &
               faceR(elemI(ii,ei_Mface_uL),fr_Flowrate),  &
               onehalfR*(faceR(elemI(ii,ei_Mface_uL),fr_Depth_u) + faceR(elemI(ii,ei_Mface_uL),fr_Depth_d)), &
               faceR(elemI(ii,ei_Mface_uL),fr_Head_u), faceR(elemI(ii,ei_Mface_uL),fr_Head_d) 
            end if

            write(*,"(i3,A,L, 3(f12.4),A,A)") ii, ' ',&
               elemYN(ii,eYN_isZeroDepth), &
               elemR(ii,er_Flowrate), &
               elemR(ii,er_Depth), elemR(ii,er_Head) ,&
               ' ',  trim(reverseKey(elemI(ii,ei_elementType)))

            if (elemI(ii,ei_Mface_dL) .ne. nullvalueI) then 
               write(*,"(i3,i3, 4(f12.4))") ii, &
                  faceI(elemI(ii,ei_Mface_uL),fi_zeroDepth), &
                  faceR(elemI(ii,ei_Mface_dL),fr_Flowrate), &
                  onehalfR*(faceR(elemI(ii,ei_Mface_dL),fr_Depth_u) + faceR(elemI(ii,ei_Mface_dL),fr_Depth_d)), &
                  faceR(elemI(ii,ei_Mface_dL),fr_Head_u), faceR(elemI(ii,ei_Mface_dL),fr_Head_d)
            end if

         elseif (elemI(ii,ei_elementType) == JM) then 
            !% --- upstream JB
            write(*,"(i3,i3,4(f12.4))") ii+1,  &
               faceI(elemI(ii+1,ei_Mface_uL),fi_zeroDepth), &
               faceR(elemI(ii+1,ei_Mface_uL),fr_Flowrate), &
               onehalfR*(faceR(elemI(ii+1,ei_Mface_uL),fr_Depth_u) + faceR(elemI(ii+1,ei_Mface_uL),fr_Depth_d)), &
               faceR(elemI(ii+1,ei_Mface_uL),fr_Head_u), faceR(elemI(ii+1,ei_Mface_uL),fr_Head_d) 

            write(*,"(i3,A,L,3(f12.4),A,A)") ii+1, ' ', &
               elemYN(ii+1,eYN_isZeroDepth), &
               elemR(ii+1,er_Flowrate),&
               elemR(ii+1,er_Depth), elemR(ii+1,er_Head) ,&
               ' ',  trim(reverseKey(elemI(ii+1,ei_elementType)))

            write(*,"(i3,A, L,3(f12.4),A,A)") ii  , ' ',&
               elemYN(ii,eYN_isZeroDepth), &
               elemR(ii,er_Flowrate),&
               elemR(ii,er_Depth), elemR(ii  ,er_Head) , &
               ' ',  trim(reverseKey(elemI(ii  ,ei_elementType)))

            write(*,"(i3,A, L,3(f12.4),A,A)") ii+2, ' ',&
               elemYN(ii+2,eYN_isZeroDepth), &
               elemR (ii+2,er_Flowrate), &
               elemR(ii+2,er_Depth), elemR(ii+2,er_Head),   &
               ' ',  trim(reverseKey(elemI(ii+2,ei_elementType)))

            !% --- downstream JB
            write(*,"(i3,i3, 4(f12.4))") ii+2, &
               faceI(elemI(ii+1,ei_Mface_uL),fi_zeroDepth), & 
               faceR(elemI(ii+2,ei_Mface_dL),fr_Flowrate), &
               onehalfR*(faceR(elemI(ii+2,ei_Mface_dL),fr_Depth_u) + faceR(elemI(ii+2,ei_Mface_dL),fr_Depth_d)), &
               faceR(elemI(ii+2,ei_Mface_dL),fr_Head_u), faceR(elemI(ii+2,ei_Mface_dL),fr_Head_d)
         endif
      end do


      !    if ((elemI(ii,ei_elementType) == CC) .or. &
      !       ((elemI(ii,ei_elementType) == JB) .and. (elemSI(ii,esi_JunctionBranch_Exists) == oneI))) then

      !       if (elemI(ii,ei_Mface_uL) .ne. nullvalueI) then 
      !          write(*,"(i4,3(f12.4))") ii, faceR(elemI(ii,ei_Mface_uL),fr_Flowrate), &
      !          faceR(elemI(ii,ei_Mface_uL),fr_Head_u), faceR(elemI(ii,ei_Mface_uL),fr_Head_d) 
      !       end if

      !       write(*,"(i4,2(f12.4),A,A)") ii, elemR(ii,er_Flowrate), elemR(ii,er_Head) ,' ',  trim(reverseKey(elemI(ii,ei_elementType)))

      !       if (elemI(ii,ei_Mface_dL) .ne. nullvalueI) then 
      !          write(*,"(i4,3(f12.4))") ii, faceR(elemI(ii,ei_Mface_dL),fr_Flowrate), &
      !           faceR(elemI(ii,ei_Mface_dL),fr_Head_u), faceR(elemI(ii,ei_Mface_dL),fr_Head_d)
      !       end if
      !    elseif (elemI(ii,ei_elementType) == JM) then 
      !       write(*,"(i4,2(f12.4),A,A)") ii, elemR(ii,er_Flowrate), elemR(ii,er_Head) ,' ',  trim(reverseKey(elemI(ii,ei_elementType)))
      !    end if
      ! end do

      !   print *, ' '
      !   print *, elemR(5,er_Head), faceR(elemI(5,ei_Mface_dL),fr_Head_u),  faceR(elemI(5,ei_Mface_dL),fr_Head_d)
      !   print *, ' '

      !   print *, ' '
      !   print *, 'flowrate ',elemR(1,er_Flowrate) !%, elemR(9,er_Flowrate)

      !   do ii=1,N_elem(1)
      !       !print *, ii, trim(reverseKey(elemI(ii,ei_elementType))), elemR(ii,er_Depth)
      !       print *, ii, elemR(ii,er_Head), elemR(ii,er_Zbottom)
      !   end do


        return

      !   print *, ' '
      !   print *, 'face downstream zero '
      !   print *, faceP(:,fp_elem_downstream_is_zero)
      !   print *, ' '
      !   print *, 'face upstream zero '
      !   print *, faceP(:,fp_elem_upstream_is_zero)
      !   print *, ' '
      !   print *, 'face both zero'
      !   print *, faceP(:,fp_elem_bothsides_are_zero)
      !   print *, ' '

      !   print *, ' '
      !   print *, 'V, Q  ', elemR(49,er_Velocity), elemR(49,er_Flowrate)
      !   print *, 'FaceU ', faceR(42,fr_Velocity_d),faceR(42,fr_Velocity_u)
      !   print *, 'FaceU ', faceR(33,fr_Velocity_d),faceR(33,fr_Velocity_u)
      !   print *, 'FaceQ ', faceR(42,fr_Flowrate)
      !   print *, 'FaceQ ', faceR(33,fr_Flowrate)
      !   print *, 'FaceA ', faceR(42,fr_Area_d),faceR(42,fr_Area_u)
      !   print *, 'FaceA ', faceR(33,fr_Area_d),faceR(33,fr_Area_u)
      !   print *, ' '

      !   print *, ' '
      !   print *, 'Weights '
      !   !print *, 'H ',elemR(51,er_InterpWeight_uH),elemR(51,er_InterpWeight_dH)
      !   print *, 'G ',elemR(51,er_InterpWeight_uG),elemR(51,er_InterpWeight_dG)
      !   !print *, 'P ',elemR(51,er_InterpWeight_uP),elemR(51,er_InterpWeight_dP)
      !   !print *, 'Q ',elemR(51,er_InterpWeight_uQ),elemR(51,er_InterpWeight_dQ)
      !   print *, ' '
      ! !   print *, '?? ',elemR(51,35), elemR(51,36)
      !   print *, ' '
      !   print *, 'H' ,er_InterpWeight_uH,er_InterpWeight_dH
      !   print *, 'G' ,er_InterpWeight_uG,er_InterpWeight_dG
      !   print *, 'P' ,er_InterpWeight_uP,er_InterpWeight_dP
      !   print *, 'Q' ,er_InterpWeight_uQ,er_InterpWeight_dQ


   
      !   do ii=1,size(iet)
      !       print *, ' '
      !       print *, ii, '-------------------------------------'
      !       print *, 'elemI ',iet(ii)
      !       print *, ' node, link: ', elemI(iet(ii),ei_node_Gidx_SWMM), elemI(iet(ii),ei_link_Gidx_BIPquick)
      !       if (elemI(iet(ii),ei_link_Gidx_BIPquick) .ne. nullvalueI) then
      !           print *, 'link name    ', trim(link%Names(elemI(iet(ii),ei_link_Gidx_SWMM))%str)
      !       elseif (elemI(iet(ii),ei_node_Gidx_SWMM) .ne. nullvalueI) then
      !           print *, 'node name    ', trim(node%Names(elemI(iet(ii),ei_node_Gidx_SWMM))%str)
      !           if ( elemI(iet(ii),ei_elementType) == JB) then 
      !               print *, 'JM is ',elemSI(iet(ii),esi_JunctionBranch_Main_Index)
      !           end if
      !       end if
      !       print *, 'length       ', elemR(iet(ii),er_Length)
      !       print *, 'element type ', trim(reverseKey(elemI(iet(ii),ei_elementType))), elemI(iet(ii),ei_elementType)
      !   end do
      !   print *, ' '
      !   print *, 'up, face:    ',elemI(iet,ei_Mface_uL)
      !   print *, 'dn  face:    ',elemI(iet,ei_Mface_dL)
      !   print *, ' '
      !   print *, 'up elem      ',faceI(ift,fi_Melem_uL)
      !   print *, 'dn elem      ',faceI(ift,fi_Melem_dL)
      !   print *, ' '



      ! !  stop 4098734


      !   print *, 'JM nodes'
      !   do ii = 1,N_elem(1)
      !       !print *, ii, elemI(ii,ei_node_Gidx_SWMM)
      !       if (elemI(ii,ei_node_Gidx_SWMM) .ne. nullvalueI) then 
      !           if (elemI(ii,ei_elementType) .eq. JM) then 
      !               print *, ii, elemI(ii,ei_node_Gidx_SWMM), trim(node%Names(elemI(ii,ei_node_Gidx_SWMM))%str)
      !           end if
      !       end if
      !   end do

      !   print *, ' '
      !   print *, 'CC elements'
      !   do ii = 1,N_elem(1)
      !       !print *, ii, elemI(ii,ei_link_Gidx_SWMM)
      !       if (elemI(ii,ei_link_Gidx_SWMM) .ne. nullvalueI) then 
      !           !print *, elemI(ii,ei_elementType), trim(reverseKey(elemI(ii,ei_elementType))), CC
      !           if (elemI(ii,ei_elementType) .eq. CC) then 
      !               print *, ii, elemI(ii,ei_link_Gidx_SWMM), trim(link%Names(elemI(ii,ei_link_Gidx_SWMM))%str), elemR(ii,er_Length)
      !           end if
      !       end if
      !   end do
      !  stop 5587623


   

      !print *, 'volume ', elemR(50,er_Volume)

      
      ! print *, ' '
      ! print *, faceR(43,fr_Head_u), faceR(43,fr_Head_u), elemR(61, er_Head)
      ! print *, ' '
      

      write(*,"(24A)") '                 ',&
      trim(reverseKey(elemI(iet(1),ei_elementType))),&
      '               ', &
      '               ', &
      trim(reverseKey(elemI(iet(2),ei_elementType))),&
      '         ', &
      trim(reverseKey(elemI(iet(3),ei_elementType))),&
      '         ', &
      trim(reverseKey(elemI(iet(4),ei_elementType))),&
      '               ', &
      '               ', &
      trim(reverseKey(elemI(iet(5),ei_elementType))),&
      '               ', &
      '               ', &
      trim(reverseKey(elemI(iet(6),ei_elementType))),&
      '         ', &
      trim(reverseKey(elemI(iet(7),ei_elementType))),&
      '         ', &
      trim(reverseKey(elemI(iet(8),ei_elementType))),&
      '               ', &
      '               '


      write(*,"(A,L3, A, A,L3, A, L3,A, L3, A,A, L3, A,A, L3, A,   L3, A,  L3, A,A, L3, A,A, L3)"), '               ',&
      elemYN(iet(1),eYN_isZeroDepth),&
      '               ', &
      '               ', &
      elemYN(iet(2),eYN_isZeroDepth),&
      '       ', &
      elemYN(iet(3),eYN_isZeroDepth),&
      '       ', &
      elemYN(iet(4),eYN_isZeroDepth),&
      '               ', &
      '               ', &
      elemYN(iet(5),eYN_isZeroDepth),&
      '               ', &
      '               ', &
      elemYN(iet(6),eYN_isZeroDepth),&
      '       ', &
      elemYN(iet(7),eYN_isZeroDepth),&
      '       ', &
      elemYN(iet(8),eYN_isZeroDepth),&
      '               ', &
      '               '


      write(*,"(A, L3, A, A, L3, A,L3,A, L3, A,A, L3, A,A, L3, A,   L3, A,  L3, A,A, L3, A,A, L3)"), '               ',&
      elemYN(iet(1),eYN_isSmallDepth),&
      '               ', &
      '               ', &
      elemYN(iet(2),eYN_isSmallDepth),&
      '       ', &
      elemYN(iet(3),eYN_isSmallDepth),&
      '       ', &
      elemYN(iet(4),eYN_isSmallDepth),&
      '               ', &
      '               ', &
      elemYN(iet(5),eYN_isSmallDepth),&
      '               ', &
      '               ', &
      elemYN(iet(6),eYN_isSmallDepth),&
      '       ', &
      elemYN(iet(7),eYN_isSmallDepth),&
      '       ', &
      elemYN(iet(8),eYN_isSmallDepth),&
      '               ', &
      '               ', &
      elemYN(iet(9),eYN_isSmallDepth), &
      '               ', &
      '               '

      write(*,"(45A)") '           ',&
      trim(link%Names(elemI(iet(1),ei_link_Gidx_SWMM))%str), &
      '           ', &
      '             ', &
      '             ', &
      trim(node%Names(elemI(iet(3),ei_node_Gidx_SWMM))%str), &
      '           ', &
      '             ', &
      '             ', &
      trim(link%Names(elemI(iet(5),ei_link_Gidx_SWMM))%str), &
      '           ', &
      '             ', &
      '             ', &
      trim(node%Names(elemI(iet(7),ei_node_Gidx_SWMM))%str), &
      '           ', &
      '             ', &
      '             '

      write(*,"(35A)") '            ', &
      '    elem   ',&
      '    face   ',&
      '    face   ',&
      '    elem   ',&
      '    elem   ',&
      '    elem   ',&
      '    face   ',&
      '    face   ',&
      '    elem   ',&
      '    face   ',&
      '    face   ',&
      '    elem   ',&
      '    elem   ',&
      '    elem   ',&
      '    face   ',&
      '    face   ',&
      '    elem   '
                                  

      write(*,"(A,16f11.3)") '     H    '  ,          &
        (elemR(iet(1),er_Head) -198.12d0)* 1000.d0,        &
      (faceR(ift(1),fr_Head_u) -198.12d0)* 1000.d0,        &
      (faceR(ift(1),fr_Head_d)-198.12d0)* 1000.d0,        &
        (elemR(iet(2),er_Head) -198.12d0)* 1000.d0,        &
        (elemR(iet(3),er_Head) -198.12d0)* 1000.d0,        &
        (elemR(iet(4),er_Head) -198.12d0)* 1000.d0,        &
      (faceR(ift(2),fr_Head_u) -198.12d0)* 1000.d0,        &
      (faceR(ift(2),fr_Head_d)-198.12d0)* 1000.d0,        &
        (elemR(iet(5),er_Head)-198.12d0)* 1000.d0,        &
      (faceR(ift(3),fr_Head_u)-198.12d0)* 1000.d0,        &
      (faceR(ift(3),fr_Head_d)-198.12d0)* 1000.d0,        &
        (elemR(iet(6),er_Head)-198.12d0)* 1000.d0,        &
        (elemR(iet(7),er_Head)-198.12d0)* 1000.d0,        &
        (elemR(iet(8),er_Head)-198.12d0)* 1000.d0,        &
      (faceR(ift(4),fr_Head_u)-198.12d0)* 1000.d0,        &
      (faceR(ift(4),fr_Head_u)-198.12d0)* 1000.d0

      write(*,"(A,16e11.3)") '     V    '  ,          &
      elemR(iet(1),er_Volume)  ,        &
        0.0,        &
      0.0,        &
      0.0,        &
      elemR(iet(3),er_Volume)  ,        &
      0.0,        &
      0.0,        &
      0.0,        &
        elemR(iet(5),er_Volume) ,        &
        0.0,        &
      0.0,        &
        0.0,        &
        elemR(iet(7),er_Volume) ,        &
        0.0,        &
      0.0,        &
      0.0

      write(*,"(A,16f11.3)") '     Z    '  ,          &
        elemR(iet(1),er_Zbottom),        &
      faceR(ift(1),fr_Zbottom),        &
      faceR(ift(1),fr_Zbottom),         &
        elemR(iet(2),er_Zbottom),        &
        elemR(iet(3),er_Zbottom),        &
        elemR(iet(4),er_Zbottom),        &
      faceR(ift(2),fr_Zbottom),        &
      faceR(ift(2),fr_Zbottom),        &
        elemR(iet(5),er_Zbottom),        &
      faceR(ift(3),fr_Zbottom),        &
      faceR(ift(3),fr_Zbottom),        &
        elemR(iet(6),er_Zbottom),        &
        elemR(iet(7),er_Zbottom),        &
        elemR(iet(8),er_Zbottom),        &
      faceR(ift(4),fr_Zbottom),        &
      faceR(ift(4),fr_Zbottom)


      write(*,"(A,16e11.3)") '   H-Z    '  ,          &
        elemR(iet(1),er_Head) -  elemR(iet(1),er_Zbottom),        &
      faceR(ift(1),fr_Head_u) -faceR(ift(1),fr_Zbottom),        &
      faceR(ift(1),fr_Head_d) -faceR(ift(1),fr_Zbottom),        &
        elemR(iet(2),er_Head) -  elemR(iet(2),er_Zbottom),        &
        elemR(iet(3),er_Head) -  elemR(iet(3),er_Zbottom),        &
        elemR(iet(4),er_Head) -  elemR(iet(4),er_Zbottom),        &
      faceR(ift(2),fr_Head_u) -faceR(ift(2),fr_Zbottom),        &
      faceR(ift(2),fr_Head_d) -faceR(ift(2),fr_Zbottom),        &
        elemR(iet(5),er_Head) -  elemR(iet(5),er_Zbottom),        &
      faceR(ift(3),fr_Head_u) -faceR(ift(3),fr_Zbottom),        &
      faceR(ift(3),fr_Head_d) -faceR(ift(3),fr_Zbottom),        &
        elemR(iet(6),er_Head) -  elemR(iet(6),er_Zbottom),        &
        elemR(iet(7),er_Head) -  elemR(iet(7),er_Zbottom),        &
        elemR(iet(8),er_Head) -  elemR(iet(8),er_Zbottom),        &
      faceR(ift(4),fr_Head_u) -faceR(ift(4),fr_Zbottom),        &
      faceR(ift(4),fr_Head_u) -faceR(ift(4),fr_Zbottom)

      write(*,"(A,16e11.3)") '     A    '  ,          &
      elemR(iet(1),er_Area),        &
    faceR(ift(1),fr_Area_u),        &
    faceR(ift(1),fr_Area_d),        &
      elemR(iet(2),er_Area),        &
      elemR(iet(3),er_Area),        &
      elemR(iet(4),er_Area),        &
    faceR(ift(2),fr_Area_u),        &
    faceR(ift(2),fr_Area_d),        &
      elemR(iet(5),er_Area),        &
    faceR(ift(3),fr_Area_u),        &
    faceR(ift(3),fr_Area_d),        &
      elemR(iet(6),er_Area),        &
      elemR(iet(7),er_Area),        &
      elemR(iet(8),er_Area),        &
    faceR(ift(4),fr_Area_u),        &
    faceR(ift(4),fr_Area_u)

      write(*,"(A,16e11.3)") '     D    '  ,          &
        elemR(iet(1),er_Depth),        &
      faceR(ift(1),fr_Depth_u),        &
      faceR(ift(1),fr_Depth_d),        &
        elemR(iet(2),er_Depth),        &
        elemR(iet(3),er_Depth),        &
        elemR(iet(4),er_Depth),        &
      faceR(ift(2),fr_Depth_u),        &
      faceR(ift(2),fr_Depth_d),        &
        elemR(iet(5),er_Depth),        &
      faceR(ift(3),fr_Depth_u),        &
      faceR(ift(3),fr_Depth_d),        &
        elemR(iet(6),er_Depth),        &
        elemR(iet(7),er_Depth),        &
        elemR(iet(8),er_Depth),        &
      faceR(ift(4),fr_Depth_u),        &
      faceR(ift(4),fr_Depth_u)

    write(*,"(A,16e11.3)") '     Q    '  ,          &
      elemR(iet(1),er_Flowrate),        &
    faceR(ift(1),fr_Flowrate),        &
    faceR(ift(1),fr_Flowrate),        &
      elemR(iet(2),er_Flowrate),        &
      elemR(iet(3),er_Flowrate),        &
      elemR(iet(4),er_Flowrate),        &
    faceR(ift(2),fr_Flowrate),        &
    faceR(ift(2),fr_Flowrate),        &
      elemR(iet(5),er_Flowrate),        &
    faceR(ift(3),fr_Flowrate),        &
    faceR(ift(3),fr_Flowrate),        &
      elemR(iet(6),er_Flowrate),        &
      elemR(iet(7),er_Flowrate),        &
      elemR(iet(8),er_Flowrate),        &
    faceR(ift(4),fr_Flowrate),        &
    faceR(ift(4),fr_Flowrate)


    write(*,"(A,16e11.3)") '   Vel    '  ,          &
    elemR(iet(1),er_Velocity),        &
  faceR(ift(1),fr_Velocity_u),        &
  faceR(ift(1),fr_Velocity_d),        &
    elemR(iet(2),er_Velocity),        &
    elemR(iet(3),er_Velocity),        &
    elemR(iet(4),er_Velocity),        &
  faceR(ift(2),fr_Velocity_u),        &
  faceR(ift(2),fr_Velocity_d),        &
    elemR(iet(5),er_Velocity),        &
  faceR(ift(3),fr_Velocity_u),        &
  faceR(ift(3),fr_Velocity_d),        &
    elemR(iet(6),er_Velocity),        &
    elemR(iet(7),er_Velocity),        &
    elemR(iet(8),er_Velocity),        &
  faceR(ift(4),fr_Velocity_u),        &
  faceR(ift(4),fr_Velocity_d)

   !    print *, ' '

   ! write(*,"(A,16e11.3)") ' Qcons    '  ,          &
   !   elemR(iet(1),er_Flowrate),        &
   ! faceR(ift(1),fr_Flowrate_Conservative),        &
   ! faceR(ift(1),fr_Flowrate_Conservative),        &
   !    elemR(iet(2),er_Flowrate),        &
   !    elemR(iet(3),er_Flowrate),        &
   !    elemR(iet(4),er_Flowrate),        &
   !  faceR(ift(2),fr_Flowrate_Conservative),        &
   !  faceR(ift(2),fr_Flowrate_Conservative),        &
   !    elemR(iet(5),er_Flowrate),        &
   !  faceR(ift(3),fr_Flowrate_Conservative),        &
   !  faceR(ift(3),fr_Flowrate_Conservative),        &
   !    elemR(iet(6),er_Flowrate),        &
   !    elemR(iet(7),er_Flowrate),        &
   !    elemR(iet(8),er_Flowrate), &
   !  faceR(ift(4),fr_Flowrate_Conservative),        &
   !  faceR(ift(4),fr_Flowrate_Conservative)

   ! ! print *, ' '


        ! stop 2098734




    end subroutine util_utest_CLprint    
!%
!%==========================================================================
!%
  subroutine util_utest_local_global

    !% In this subroutine we are checking the the local and global indexs of the link,node,elem and face arrays are unique

    integer ii, jj, kk, min_val, max_val
    logical dup_found
    character(64) :: subroutine_name = 'local_global_unique'
    if (setting%Debug%File%initialization) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !% Looping through the array and finding all of the unqiue values
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

    !% checking node%I(:,:) indexes

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


    !% checking elemI(:,:) local indexes

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



    !% checking faceI(:,:) indexes

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

    !% checking faceI(:,:) global indexes

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


    !% checking elemI(:,:) global indexes

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


    if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

  end subroutine util_utest_local_global
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_utest_pack_arrays

      !% Going through all of the pack arrays and making sure they are unique, this follows the same process as the subroutine above, with a slight change.
      integer ii, jj, kk, min_val, max_val
      logical dup_found
      character(64) :: subroutine_name = 'pack_arrays_unique'
      if (setting%Debug%File%initialization) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"



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
      if (min_val == nullvalueI) then
         ii = ii - 1
      end if


      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_AC) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemP(ep_ac)) then
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

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if

      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_ALLtm) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemP(ep_ALLtm)) then
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

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if

      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CC_AC) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemP(ep_CC_AC)) then
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

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if




      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CC_ALLtm) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemP(ep_CC_Alltm)) then
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

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if


      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CC_ETM) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemP(ep_CC_ETM)) then
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

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if


      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CC_H_ETM) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemP(ep_CC_H_ETM)) then
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

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if


      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CC_Q_AC) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemP(ep_CC_Q_AC)) then
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

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if


      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_CC_Q_ETM) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemP(ep_CC_Q_ETM)) then
         print *, "ERROR:::: elemP(:,ep_CC_Q_ETM) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_CC_Q_ETM) is unique. This_image ::", this_image()
      end if

      !% checking elemP(:,ep_CCJB_ACsurcharged) indexes


      ! min_val = minval(elemP(:,ep_CCJB_ACsurcharged)) - 1
      ! max_val = maxval(elemP(:,ep_CCJB_ACsurcharged))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_CCJB_ACsurcharged),mask=elemP(:,ep_CCJB_ACsurcharged)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_CCJB_ACsurcharged) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_CCJB_ACsurcharged)) then
      !    print *, "ERROR:::: elemP(:,ep_CCJB_ACsurcharged) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_CCJB_ACsurcharged) is unique. This_image ::", this_image()
      ! end if

      !% checking elemP(:,ep_CCJB_ALLtm) indexes


      ! min_val = minval(elemP(:,ep_CCJB_ALLtm)) - 1
      ! max_val = maxval(elemP(:,ep_CCJB_ALLtm))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_CCJB_ALLtm),mask=elemP(:,ep_CCJB_ALLtm)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_CCJB_ALLtm) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_CCJB_ALLtm)) then
      !    print *, "ERROR:::: elemP(:,ep_CCJB_ALLtm) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_CCJB_ALLtm) is unique. This_image ::", this_image()
      ! end if

      !% checking elemP(:,ep_CCJB_AC) indexes


      ! min_val = minval(elemP(:,ep_CCJB_AC)) - 1
      ! max_val = maxval(elemP(:,ep_CCJB_AC))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_CCJB_AC),mask=elemP(:,ep_CCJB_AC)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_CCJB_AC) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_CCJB_AC)) then
      !    print *, "ERROR:::: elemP(:,ep_CCJB_AC) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_CCJB_AC) is unique. This_image ::", this_image()
      ! end if

      ! !% checking elemP(:,ep_CCJB_ALLtm_ACsurcharged) indexes

      ! min_val = minval(elemP(:,ep_CCJB_ALLtm_ACsurcharged)) - 1
      ! max_val = maxval(elemP(:,ep_CCJB_ALLtm_ACsurcharged))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_CCJB_ALLtm_ACsurcharged),mask=elemP(:,ep_CCJB_ALLtm_ACsurcharged)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_CCJB_ALLtm_ACsurcharged) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_CCJB_ALLtm_ACsurcharged)) then
      !    print *, "ERROR:::: elemP(:,ep_CCJB_ALLtm_ACsurcharged) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_CCJB_ALLtm_ACsurcharged) is unique. This_image ::", this_image()
      ! end if



      !% checking elemP(:,ep_CCJB_eETM_i_fAC) indexes

      !% HACK -- 20211212brh problems with CCJB_eETM_i_fAC array -- commented out
      !% DO NOT DELETE -- NEEDS TO BE RESTORED AND DEBUGGED
      ! min_val = minval(elemP(:,ep_CCJB_eETM_i_fAC)) - 1
      ! max_val = maxval(elemP(:,ep_CCJB_eETM_i_fAC))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_CCJB_eETM_i_fAC),mask=elemP(:,ep_CCJB_eETM_i_fAC)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if

      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_CCJB_eETM_i_fAC) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= N_elem(this_image())) then
      !    print *, "ERROR:::: elemP(:,ep_CCJB_eETM_i_fAC) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_CCJB_eETM_i_fAC) is unique. This_image ::", this_image()
      ! end if

      
      !% checking elemP(:,ep_CCJB_ETM) indexes
      ! min_val = minval(elemP(:,ep_CCJB_ETM)) - 1
      ! max_val = maxval(elemP(:,ep_CCJB_ETM))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_CCJB_ETM),mask=elemP(:,ep_CCJB_ETM)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_CCJB_ETM) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_CCJB_ETM)) then
      !    print *, "ERROR:::: elemP(:,ep_CCJB_ETM) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_CCJB_ETM) is unique. This_image ::", this_image()
      ! end if

      !% checking elemP(:,ep_Diag) indexes


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

      !% checking elemP(:,ep_ETM) indexes

      min_val = minval(elemP(:,ep_ETM)) - 1
      max_val = maxval(elemP(:,ep_ETM))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemP(:,ep_ETM),mask=elemP(:,ep_ETM)>min_val)
      end do

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if


      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_ETM) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemP(ep_ETM)) then
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

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if


      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_JM_AC) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemP(ep_JM_AC)) then
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

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if


      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemP(:,ep_JM_ALLtm) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemP(ep_JM_ALLtm)) then
         print *, "ERROR:::: elemP(:,ep_JM_ALLtm) is not unique. This_image ::", this_image()

      else
         print *, "elemP(:,ep_JM_ALLtm) is unique. This_image ::", this_image()
      end if


      !% checking elemP(:,ep_JB_ETM) indexes

      ! min_val = minval(elemP(:,ep_JB_ETM)) - 1
      ! max_val = maxval(elemP(:,ep_JB_ETM))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_JB_ETM),mask=elemP(:,ep_JB_ETM)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_JB_ETM) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(EP_JB_ETM)) then
      !    print *, "ERROR:::: elemP(:,ep_JB_ETM) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_JB_ETM) is unique. This_image ::", this_image()
      ! end if

      !% checking elemP(:,ep_AC_ACnonSurcharged) indexes


      ! min_val = minval(elemP(:,ep_AC_ACnonSurcharged)) - 1
      ! max_val = maxval(elemP(:,ep_AC_ACnonSurcharged))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_AC_ACnonSurcharged),mask=elemP(:,ep_AC_ACnonSurcharged)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_AC_ACnonSurcharged) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_AC_ACnonSurcharged)) then
      !    print *, "ERROR:::: elemP(:,ep_AC_ACnonSurcharged) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_AC_ACnonSurcharged) is unique. This_image ::", this_image()
      ! end if

      !% checking elemP(:,ep_ALLtm_NonSurcharged) indexes

      ! min_val = minval(elemP(:,ep_ALLtm_NonSurcharged)) - 1
      ! max_val = maxval(elemP(:,ep_ALLtm_NonSurcharged))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_ALLtm_NonSurcharged),mask=elemP(:,ep_ALLtm_NonSurcharged)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_ALLtm_NonSurcharged) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_ALLtm_NonSurcharged)) then
      !    print *, "ERROR:::: elemP(:,ep_ALLtm_NonSurcharged) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_ALLtm_NonSurcharged) is unique. This_image ::", this_image()
      ! end if

      !% checking elemP(:,ep_ETM_PSnonSurcharged) indexes


      ! min_val = minval(elemP(:,ep_ETM_PSnonSurcharged)) - 1
      ! max_val = maxval(elemP(:,ep_ETM_PSnonSurcharged))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_ETM_PSnonSurcharged),mask=elemP(:,ep_ETM_PSnonSurcharged)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_ETM_PSnonSurcharged) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_ETM_PSnonSurcharged)) then
      !    print *, "ERROR:::: elemP(:,ep_ETM_PSnonSurcharged) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_ETM_PSnonSurcharged) is unique. This_image ::", this_image()
      ! end if


      !% checking elemP(:,ep_SmallDepth_CC_ALLtm_posSlope) indexes


      ! min_val = minval(elemP(:,ep_SmallDepth_CC_ALLtm_posSlope)) - 1
      ! max_val = maxval(elemP(:,ep_SmallDepth_CC_ALLtm_posSlope))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_SmallDepth_CC_ALLtm_posSlope),mask=elemP(:,ep_SmallDepth_CC_ALLtm_posSlope)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_SmallDepth_CC_ALLtm_posSlope) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_SmallDepth_CC_ALLtm_posSLope)) then
      !    print *, "ERROR:::: elemP(:,ep_SmallDepth_CC_ALLtm_posSlope) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_SmallDepth_CC_ALLtm_posSlope) is unique. This_image ::", this_image()
      ! end if

      ! !% checking elemP(:,ep_SmallDepth_CC_ALLtm_negSlope) indexes


      ! min_val = minval(elemP(:,ep_SmallDepth_CC_ALLtm_negSlope)) - 1
      ! max_val = maxval(elemP(:,ep_SmallDepth_CC_ALLtm_negSLope))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_SmallDepth_CC_ALLtm_negSlope),mask=elemP(:,ep_SmallDepth_CC_ALLtm_negSlope)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_SmallDepth_CC_ALLtm_negSlope) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_SmallDepth_CC_ALLtm_negSlope)) then
      !    print *, "ERROR:::: elemP(:,ep_SmallDepth_CC_ALLtm_negSlope) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_SmallDepth_CC_ALLtm_negSlope) is unique. This_image ::", this_image()
      ! end if

      ! !% checking elemP(:,ep_SmallDepth_CC_ALLtm) indexes

      ! min_val = minval(elemP(:,ep_SmallDepth_CC_ALLtm)) - 1
      ! max_val = maxval(elemP(:,ep_SmallDepth_CC_ALLtm))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_SmallDepth_CC_ALLtm),mask=elemP(:,ep_SmallDepth_CC_ALLtm)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_SmallDepth_CC_ALLtm) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_SmallDepth_CC_ALLtm)) then
      !    print *, "ERROR:::: elemP(:,ep_SmallDepth_CC_ALLtm) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_SmallDepth_CC_ALLtm) is unique. This_image ::", this_image()
      ! end if


      !% checking elemP(:,ep_ACsurcharged) indexes


      ! min_val = minval(elemP(:,ep_ACsurcharged)) - 1
      ! max_val = maxval(elemP(:,ep_ACsurcharged))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_ACsurcharged),mask=elemP(:,ep_ACsurcharged)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_ACsurcharged) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_ACsurcharged)) then
      !    print *, "ERROR:::: elemP(:,ep_ACsurcharged) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_ACsurcharged) is unique. This_image ::", this_image()
      ! end if


      !% checking elemP(:,ep_ALLtmSurcharged) indexes

      ! min_val = minval(elemP(:,ep_ALLtmSurcharged)) - 1
      ! max_val = maxval(elemP(:,ep_ALLtmSurcharged))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_ALLtmSurcharged),mask=elemP(:,ep_ALLtmSurcharged)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_ALLtmSurcharged) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= N_elem(this_image())) then
      !    print *, "ERROR:::: elemP(:,ep_ALLtmSurcharged) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_ALLtmSurcharged) is unique. This_image ::", this_image()
      ! end if


      !% checking elemP(:,ep_PSsurcharged) indexes


      ! min_val = minval(elemP(:,ep_PSsurcharged)) - 1
      ! max_val = maxval(elemP(:,ep_PSsurcharged))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_PSsurcharged),mask=elemP(:,ep_PSsurcharged)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_PSsurcharged) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_PSsurcharged)) then
      !    print *, "ERROR:::: elemP(:,ep_PSsurcharged) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_PSsurcharged) is unique. This_image ::", this_image()
      ! end if


      !% checking elemP(:,ep_CCJM_H_ACsurcharged) indexes

      ! min_val = minval(elemP(:,ep_CCJM_H_ACsurcharged)) - 1
      ! max_val = maxval(elemP(:,ep_CCJM_H_ACsurcharged))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_CCJM_H_ACsurcharged),mask=elemP(:,ep_CCJM_H_ACsurcharged)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if

      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_CCJM_H_ACsurcharged) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_CCJM_H_ACsurcharged)) then
      !    print *, "ERROR:::: elemP(:,ep_CCJM_H_ACsurcharged) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_CCJM_H_ACsurcharged) is unique. This_image ::", this_image()
      ! end if

      !% checking elemP(:,ep_CCJM_H_AC) indexes

      ! min_val = minval(elemP(:,ep_CCJM_H_AC)) - 1
      ! max_val = maxval(elemP(:,ep_CCJM_H_AC))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_CCJM_H_AC),mask=elemP(:,ep_CCJM_H_AC)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_CCJM_H_AC) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_CCJM_H_AC)) then
      !    print *, "ERROR:::: elemP(:,ep_CCJM_H_AC) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_CCJM_H_AC) is unique. This_image ::", this_image()
      ! end if

      !% checking elemP(:,ep_CCJB_eAC_i_fETM) indexes

      ! min_val = minval(elemP(:,ep_CCJB_eAC_i_fETM)) - 1
      ! max_val = maxval(elemP(:,ep_CCJB_eAC_i_fETM))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_CCJB_eAC_i_fETM),mask=elemP(:,ep_CCJB_eAC_i_fETM)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_CCJB_eAC_i_fETM) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= N_elem(this_image())) then
      !    print *, "ERROR:::: elemP(:,ep_CCJB_eAC_i_fETM) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_CCJB_eAC_i_fETM) is unique. This_image ::", this_image()
      ! end if

      !% checking elemPGalltm(:,epg_CC_rectangular) indexes


      ! min_val = minval(elemPGalltm(:,epg_CC_rectangular)) - 1
      ! max_val = maxval(elemPGalltm(:,epg_CC_rectangular))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemPGalltm(:,epg_CC_rectangular),&
      !         mask=elemPGalltm(:,epg_CC_rectangular)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemPGalltm(:,epg_CC_rectangular) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemPGalltm(epg_CC_rectangular)) then
      !    print *, "ERROR:::: elemPGalltm(:,epg_CC_rectangular) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemPGalltm(:,epg_CC_rectangular) is unique. This_image ::", this_image()
      ! end if

      !% checking elemPGalltm(:,epg_CC_trapezoidal) indexes

      ! min_val = minval(elemPGalltm(:,epg_CC_trapezoidal)) - 1
      ! max_val = maxval(elemPGalltm(:,epg_CC_trapezoidal))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemPGalltm(:,epg_CC_trapezoidal),&
      !         mask=elemPGalltm(:,epg_CC_trapezoidal)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemPGalltm(:,epg_CC_trapezoidal) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemPGalltm(epg_CC_trapezoidal)) then
      !    print *, "ERROR:::: elemPGalltm(:,epg_CC_trapezoidal) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemPGalltm(:,epg_CC_trapezoidal) is unique. This_image ::", this_image()
      ! end if


      !% checking elemPGalltm(:,epg_JB_rectangular) indexes


      ! min_val = minval(elemPGalltm(:,epg_JB_rectangular)) - 1
      ! max_val = maxval(elemPGalltm(:,epg_JB_rectangular))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemPGalltm(:,epg_JB_rectangular),mask=elemPGalltm(:,epg_JB_rectangular)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemPGalltm(:,epg_JB_rectangular) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemPGalltm(epg_JB_rectangular)) then
      !    print *, "ERROR:::: elemPGalltm(:,epg_JB_rectangular) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemPGalltm(:,epg_JB_rectangular) is unique. This_image ::", this_image()
      ! end if


      !% checking elemPGalltm(:,epg_JB_trapezoidal) indexes


      ! min_val = minval(elemPGalltm(:,epg_JB_trapezoidal)) - 1
      ! max_val = maxval(elemPGalltm(:,epg_JB_trapezoidal))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemPGalltm(:,epg_JB_trapezoidal),mask=elemPGalltm(:,epg_JB_trapezoidal)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemPGalltm(:,epg_JB_trapezoidal) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemPGalltm(epg_JB_trapezoidal)) then
      !    print *, "ERROR:::: elemPGalltm(:,epg_JB_trapezoidal) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemPGalltm(:,epg_JB_trapezoidal) is unique. This_image ::", this_image()
      ! end if


      !% --------------------------------------------------------------------------------------------------------------

      !% checking elemPGetm(:,epg_CC_rectangular) indexes


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

      !% checking elemPGetm(:,epg_CC_trapezoidal) indexes

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


      !% checking elemPGetm(:,epg_JB_rectangular) indexes


      ! min_val = minval(elemPGetm(:,epg_JB_rectangular)) - 1
      ! max_val = maxval(elemPGetm(:,epg_JB_rectangular))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemPGetm(:,epg_JB_rectangular),mask=elemPGetm(:,epg_JB_rectangular)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if

      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemPGetm(:,epg_JB_rectangular) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemPGetm(epg_JB_rectangular)) then
      !    print *, "ERROR:::: elemPGetm(:,epg_JB_rectangular) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemPGetm(:,epg_JB_rectangular) is unique. This_image ::", this_image()
      ! end if


      !% checking elemPGetm(:,epg_JB_trapezoidal) indexes


      ! min_val = minval(elemPGetm(:,epg_JB_trapezoidal)) - 1
      ! max_val = maxval(elemPGetm(:,epg_JB_trapezoidal))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemPGetm(:,epg_JB_trapezoidal),mask=elemPGetm(:,epg_JB_trapezoidal)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemPGetm(:,epg_JB_trapezoidal) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemPGetm(epg_JB_trapezoidal)) then
      !    print *, "ERROR:::: elemPGetm(:,epg_JB_trapezoidal) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemPGetm(:,epg_JB_trapezoidal) is unique. This_image ::", this_image()
      ! end if


      !% --------------------------------------------------------------------------------------------------------------

      !% checking elemPGac(:,epg_CCJM_rectangular_nonsurcharged) indexes


      min_val = minval(elemPGac(:,epg_CC_rectangular)) - 1
      max_val = maxval(elemPGac(:,epg_CC_rectangular))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemPGac(:,epg_CC_rectangular),&
              mask=elemPGac(:,epg_CC_rectangular)>min_val)
      end do

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if


      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemPGac(:,epg_CC_rectangular) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemPGac(epg_CC_rectangular)) then
         print *, "ERROR:::: elemPGac(:,epg_CC_rectangular) is not unique. This_image ::", this_image()

      else
         print *, "elemPGac(:,epg_CC_rectangular) is unique. This_image ::", this_image()
      end if

      !% checking elemPGac(:,epg_CC_trapezoidal) indexes

      min_val = minval(elemPGac(:,epg_CC_trapezoidal)) - 1
      max_val = maxval(elemPGac(:,epg_CC_trapezoidal))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(elemPGac(:,epg_CC_trapezoidal), &
              mask=elemPGac(:,epg_CC_trapezoidal)>min_val)
      end do

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if


      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "elemPGac(:,epg_CC_trapezoidal) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_elemPGac(epg_CC_trapezoidal)) then
         print *, "ERROR:::: elemPGac(:,epg_CC_trapezoidal) is not unique. This_image ::", this_image()

      else
         print *, "elemPGac(:,epg_CC_trapezoidal) is unique. This_image ::", this_image()
      end if


      !% checking elemPGac(:,epg_JB_rectangular) indexes


      ! min_val = minval(elemPGac(:,epg_JB_rectangular)) - 1
      ! max_val = maxval(elemPGac(:,epg_JB_rectangular))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemPGac(:,epg_JB_rectangular),mask=elemPGac(:,epg_JB_rectangular)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemPGac(:,epg_JB_rectangular) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemPGac(epg_JB_rectangular)) then
      !    print *, "ERROR:::: elemPGac(:,epg_JB_rectangular) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemPGac(:,epg_JB_rectangular) is unique. This_image ::", this_image()
      ! end if


      !% checking elemPGac(:,epg_JB_trapezoidal) indexes


      ! min_val = minval(elemPGac(:,epg_JB_trapezoidal)) - 1
      ! max_val = maxval(elemPGac(:,epg_JB_trapezoidal))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemPGac(:,epg_JB_trapezoidal),mask=elemPGac(:,epg_JB_trapezoidal)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemPGac(:,epg_JB_trapezoidal) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemPGac(epg_JB_trapezoidal)) then
      !    print *, "ERROR:::: elemPGac(:,epg_JB_trapezoidal) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemPGac(:,epg_JB_trapezoidal) is unique. This_image ::", this_image()
      ! end if

      !% checking faceP(:,fp_all) indexes


      min_val = minval(faceP(:,fp_all)) - 1
      max_val = maxval(faceP(:,fp_all))
      ii = 0

      do while(min_val<max_val)
         ii = ii + 1
         min_val = minval(faceP(:,fp_all),mask=faceP(:,fp_all)>min_val)
      end do

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if


      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "faceP(:,fp_all) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_faceP(fp_all)) then
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

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if


      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "faceP(:,fp_AC) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_faceP(fp_AC)) then
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

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if


      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "faceP(:,fp_Diag) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_faceP(fp_Diag)) then
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

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if


      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "faceP(:,fp_JumpDn) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_faceP(fp_JumpDn)) then
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

      if (min_val == nullvalueI) then
         ii = ii - 1
      end if


      if (ii == 0 .and. min_val == nullvalueI) then
         print *, "faceP(:,fp_JumpUp) is only filled with nullvalueI. This_image ::", this_image()

      else if (ii /= npack_faceP(fp_JumpUp)) then
         print *, "ERROR:::: faceP(:,fp_JumpUp) is not unique. This_image ::", this_image()

      else
         print *, "faceP(:,fp_JumpUp) is unique. This_image ::", this_image()
      end if


      if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_utest_pack_arrays
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_utest_node_link_image

      !% In this subroutine we are checking whether the every image has atleast one node and link assigned to it.
      integer :: ii, jj, kk, counter
      character(64) :: subroutine_name = 'init_face_check'
      if (setting%Debug%File%initialization) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

      kk = 1
      counter = 0

      !% We can find if there is a node on each image by counting the amount of unique values there is inside of node%I(:,ni_P_image)
      !% So we use the code below to loop through node%I(:,ni_P_image) and find all the unique values

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

      !% After that we compare the number of images we are using to what we counted.
      !% If it is correct they should be the same value, otherwise there was an error when paritioning the links and nodes

      if (num_images() /= counter) then
         print *, "error in NodeI images. This_image :: ", this_image()
      else
         print *, "correct number in node%I images. This_image :: ", this_image()
      end if

      !% We reset the counter and do the same process for the Links
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
      if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_utest_node_link_image
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_utest_slope_checking
      !% In this subroutine we are checking that all of the slopes are postive.
      !% To do this and loop through and if we find a negative slope then we exit the loop and report it.

      integer :: ii, jj
      logical :: invalid_slope
      character(64) :: subroutine_name = 'slope_checking'
      if (setting%Debug%File%initialization) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

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

      if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_utest_slope_checking
!%
!%==========================================================================
!%==========================================================================
!%

    subroutine util_utest_global_index_check
      !% In this subroutine we are checking that the global indexs are correct by adding the first valid global index of an image with the local index

      integer :: ii, current_length, counter
      character(64) :: subroutine_name = 'global_index_checking'

      if (setting%Debug%File%initialization) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

      !% here we find the current length of the global index by looking at the first value of elemI on that image and subtracting one.

      current_length = elemI(1, ei_Gidx) - 1

      !% Then we loop through elemI(:,ei_Gidx) and report and stop looping if there is an error in the global indexes
      do ii = 1, size(elemI(:,ei_Gidx))

         if (elemI(ii,ei_Gidx) /= current_length+elemI(ii,ei_Lidx) .and. elemI(ii,ei_Gidx)/= nullvalueI) then
            print *, "error in elem global indexes. Processor :: ", this_image()
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

         if (faceYN(ii,fYN_isSharedFace)) then
            cycle

         else if (faceI(ii,ei_Gidx) /= current_length+faceI(ii,ei_Lidx) .and. faceI(ii,ei_Gidx)/= nullvalueI) then
            print *, "error in face global indexes. Processor :: ", this_image()
            !% print *, "faceI(ii,ei_Gidx)", faceI(ii,ei_Gidx)
            !% print *, "faceI(ii,ei_Lidx)", faceI(ii,ei_Lidx)
            exit
         end if
      end do

      if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_utest_global_index_check
!%
!%==========================================================================
!%==========================================================================
!%
    !subroutine geometry_checking

      !integer ii, jj

      !do ii = 1, size(link%I(:,linktype))

         !if (link%I(ii, linktype) == lchannel) then

          !  select case(link%I(ii,li_geometry)

     !          case(lRectangular)

    !end subroutine geometry_checking

!%
!%==========================================================================
!% END MODULE    
!%==========================================================================
!%
  end module utility_unit_testing
