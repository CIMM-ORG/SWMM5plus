module utility_unit_testing
  use interface_
  use utility_allocate
  use discretization
  use define_indexes
  use define_keys
  use define_globals
  use define_settings
  use network_define
  use utility_crash
  !use junction_elements, only: junction_conservation_residual

  implicit none

  private

  public :: util_utest_CLprint     !% custom command line printing for debugging
  public :: util_utest_checkIsNan  !% checks if selected data columns have NaN

  public :: util_utest_local_global
  public :: util_utest_pack_arrays
  public :: util_utest_node_link_image
  public :: util_utest_slope_checking
  public :: util_utest_global_index_check
 

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

            integer  :: tE, tFu, tFd

   
            ! !% YJ_Geom_Test_mod01
            ! integer :: iet(7) = (/  1,  2,   4,3,5, 14, 15/)
            ! integer :: ift(6)  = (/1,  2,  3,      4,  13, 14/)

            !% extran1_final_mod.inp
            ! integer ::  iet(7) = (/   29,   30,    32, 31, 33,    189,   190/)
            ! integer ::  ift(6)  = (/27, 28,     29,            30,    177,   178 /)

            !% extran1_final_mod2.inp
            ! integer :: iet(7) = (/   29, 30,   32,31,33,    42,  43/)
            ! integer :: ift(6) = (/27, 28, 29,           30,   39,   40/)

            !% T001 -- 100 or smaller
            ! integer :: iet(7) = (/ 4 , 5 , 6 , 7 , 8 , 9 , 10 /)
            ! integer :: ift(8) = (/4, 5 , 6 , 7 , 8 , 9 , 10 , 11 /)

            !% T002 dx 100
            ! integer :: iet(8) = (/3,    4,   5,   7, 6, 8,   17,  18/)
            ! integer :: ift(5) = (/   4 ,  5 ,   6,         7,   16/)

            !% T004 dx 10 implies storage
            !integer :: iet(8) = (/  8,  9,   10,    12,11,13,   22,  23/)
            !integer :: ift(7) = (/8, 9 , 10 ,   11,          12,   21, 22/)

            !% Example1.inp for node 10 (JM 13)
            ! integer :: iet(8) = (/  10,  11,   12,    14,13,15,   58,  59/) !% goes to 69
            ! integer :: ift(7) = (/10,  11, 12 ,   13,          14,   54, 55/)!% goes to 46 JM 47b

            !% Exampl.inp between node 13 and 47b
            ! integer :: iet(8) = (/15, 58,   59, 60, 61, 62, 63, 64/)
            ! integer :: ift(7) = (/  14, 54,   55, 56, 57, 58, 59 /)

             !% Example1.inp for node  (JM 47a)
            ! integer :: iet(8) = (/    44,  45,    46,    48,47,49,      70,   71/) !% goes to 78
            ! integer :: ift(7) = (/ 41,   42,   43 ,   44,            45,    65,  66/) !% goes to 73 JM 79


            !  !% Example1.inp for node  (JM 47b) extends for e58 in JM13
            ! integer :: iet(8) = (/    67,  68,    69,    50,47,49,      70,   71/) !% goes to 78
            ! integer :: ift(7) = (/ 62,   63,   64 ,   46,            45,    65,  66/) !% goes to 73 JM 79
 
            !% Example1.inp for node  (JM 79)
            ! integer :: iet(8) = (/    76,  77,    78,    80,79,81,      90,   91/)
            ! integer :: ift(7) = (/ 70,   71,   72 ,   73,            74,    83,  84/)
            
            !% Example1.inp for node  (JM 168)
            ! integer :: iet(8) = (/   165,  166,   167,    169,168,170,      179,   180/)
            ! integer :: ift(7) = (/153,  154,  155 ,   156,              157,    166,  167/)
    

             !% Example1.inp for node  (JM 30)
            ! integer :: iet(8) = (/    27,  28,    29,    31,30,32,      70,   71/)
            ! integer :: ift(7) = (/ 26,   27,   28 ,   29,            45,    65,  66/)

            !% T005_CC_Free_dx0050_ForceStorage

            ! integer :: iet(8) = (/   1,   1,   2,    4,3,5,   14,     15/)
            ! integer :: ift(7) = (/ 1,   1  , 2,   3,        4,   13,    14/)

            !% T005W_RO_Free_dx0010_ForceJM
            !integer :: iet(11) = (/ 93,   94,    96,95,97,   211,    107,106,108,  117,  118/)
            !integer :: ift(6) = (/     86,    87,         88,     97,            98,  107 /)

            ! !% T005Wr_RO_Free-dx0010_weirJM
            ! integer :: iet(11) = (/  49,   50,    52,51,53,   123,    63,62,64,   73,  74/)
            ! integer :: ift(6)  = (/     50,   51,         52,     61,         62,   71 /)


             !% T005Wr_RO_Free-dx0010_faceNJ2 with only upstream node
            ! integer :: iet(8) = (/  49,   50,    52,51,53,   112,    62,  63 /)
            ! integer :: ift(5)  = (/     50,   51,         52,     61,    62 /)


            !  !% T005Wr_RO_Free-dx0010_faceNJ2
            ! integer :: iet(5) = (/     49,   50,     101,      51,   52 /)
            ! integer :: ift(6)  = (/ 49,    50,   51,       52,   53,   54 /)


            !% T007
            ! integer :: iet(8) = (/  136,   136,    138,137,139,   148,      149,   150 /)
            ! integer :: ift(5)  = (/     124,   125,            126,     135,    136 /)

           

            ! !% T009
            ! integer :: iet(8) = (/  3,  4,   5,    7,6,8,   17,    18 /)
            ! integer :: ift(5)  = (/   4,  5,   6,      7,    16  /)

            !% T007 NJ2 
            ! integer :: iet(5) = (/ 1, 2, 3 , 4, 5/)
            ! integer :: ift(6) = (/1 ,2 ,3 ,4 , 5 , 6/)

            !% T009Wr_RO_Free_dx0100_NJ2
            ! integer :: iet(5) = (/  1,    2,    3,   4,    5/)
            ! integer :: ift(6) = (/1,    2,   3,    4,   5,    6/)


            ! integer :: iet(5) = (/   3,   4,    5,  7,  8/)
            ! integer :: ift(6) = (/3,    4,   5,    6, 7,  8/)

            ! integer :: iet(5) = (/  6,    7,    8,   9,    10/)
            ! integer :: ift(6) = (/6,    7,   8,    9,  10,   11/)

            ! integer :: iet(5) = (/  11,    12,    13,   14,    15/)
            ! integer :: ift(6) = (/11,    12,   13,    14,  15,   16/)

            ! integer :: iet(5) = (/  16,    17,    18,   19,    20/)
            ! integer :: ift(6) = (/16,    17,   18,    19,  20,   21/)

            ! integer :: iet(5) = (/  21,    22,    23,   24,    25/)
            ! integer :: ift(6) = (/21,   22,   23,    24,  25,    26/)

            ! integer :: iet(5) = (/  24,    25,    51,   26,    27/)
            ! integer :: ift(6) = (/24,    25,   26,    27,  28,   29/)

            ! integer :: iet(5) = (/  28,    29,    30,   31,    32/)
            ! integer :: ift(6) = (/29,   30,   31,    32,  33,    34/)

            ! integer :: iet(5) = (/  33,    34,    35,   36,    37/)
            ! integer :: ift(6) = (/34,   35,   36,    37,  38,    39/)

            ! integer :: iet(5) = (/  38,    39,    40,   41,    42/)
            ! integer :: ift(6) = (/39,   40,   41,    42,  43,    44/)

            ! integer :: iet(5) = (/  43,    44,    45,   46,    47/)
            ! integer :: ift(6) = (/44,   45,   46,    47,  48,    49/)

            ! integer :: iet(5) = (/  46,    47,    48,   49,    50/)
            ! integer :: ift(6) = (/47,   48,   49,    50,  51,     52/)

            ! !% T00_Wr_RO_Free_dx0100
            ! integer :: iet(11) = (/  68,   69,    71,70,72,   161,    82,81,83,   92,  93/)
            ! integer :: ift(6)  = (/     61,   62,          63,     72,          73,   82 /)

            !% extran1_final
            !integer :: iet(8) = (/127,    128,     129,        131,130,132,    760,   751 /)
            !integer :: ift(5) = (/     126,    127,      128,              129,    750  /)

            ! integer :: iet(7) = (/ 750,         717,    713,      715,714,716,    751  /)
            ! integer :: ift(3) = (/         710,              708,             709  /)

            !% T007a dx10
            ! integer :: iet(8) = (/100,       102,101,103,    112,   113,    114,    115 /)
            ! integer :: ift(5) = (/     101,              102,    111,   112,   113/)

            ! % T007a dx10  8 x 5 with JB as 4th elem
            ! integer :: iet(8) = (/98,  99,   100,   102,101,103,    112,   113 /)
            ! integer :: ift(5) = (/   99, 100,  101,              102,    111/)

                     !% T007a
            ! integer :: iet(8) = (/211,       213,212,214,    223,   224,    225,    226 /)
            ! integer :: ift(5) = (/     210,              211,    220,   221,   222/)

            ! !% T007 dx333
            ! integer :: iet(7) = (/ 102, 104,   113, 114,   115,  117, 116/)
            ! integer :: ift(4) = (/           89,   98,  99,    100 /)


            !% T007a dx250 1 junction
            ! integer :: iet(8) = (/    2,    3,    4,    6,5,7,   16,      17 /)
            ! integer :: ift(7)  = (/2,    3,     4,   5,        6,     15,    16 /)

            !% T007 dx250 9th junction
            ! integer :: iet(8) = (/    122,    123,    124,     126,125,127,     136,      137 /)
            ! integer :: ift(7)  = (/106,    107,     108,   109,            110,     119,    120 /)

            !% T007 dx250 4th junction
            ! integer :: iet(8) = (/    47,    48,    49,     51,50,52,     61,      62 /)
            ! integer :: ift(7)  = (/41,    42,     43,   44,            45,     54,    55 /)

            !% T007 7 junction
            ! integer :: iet(8) = (/31,       33,32,34,    43,   44,    45,    46 /)
            ! integer :: ift(5) = (/     30,            31,    40,   41,   42/)

            ! real(8) :: tempV, tempVold, tempFin, tempFout

            integer :: iet(8) = (/  1,  2,  3,     5,4,6,  15,  16/)
            integer :: ift(7) = (/1 , 2,  3,   4,        5,   14,  15/)

            ! integer :: iet(8) = (/  43,  44,  45,     47,46,48,    57,   58/)
            ! integer :: ift(7) = (/29 , 38,  39,   40,           41,   50,   51/)

            !% T009WR dx10  weir=611   11x6
            ! integer :: iet(11) = (/ 293,      294,     296,295,297,     611,     307,306,308,     317,     318 /)
            ! integer :: ift(6) = (/       286,      287,             288,      297,             298,    307/)

            !% user5 reduced04 dx150 from upstream
            ! integer :: iet(8) = (/ 1, 2, 3, 4,    6,5,7,  16/)
            ! integer :: ift(7) = (/ 1 ,2, 3, 4, 5,        6,  15/)

            !% user5 reduced04 dx150 from downstream
            ! integer :: iet(8) = (/   4,    6,5,7,  16, 17, 18, 19/)
            ! integer :: ift(7) = (/  4, 5,        6,  15, 16,  17, 18/)

            !% User5 reduced02 dx50 from downstream
            ! integer :: iet(8) = (/   72,   73,    75,74,76,     85, 86, 87/)
            ! integer :: ift(7) = (/62,   63,    64,          65,    74, 75, 76/)


             !% User5 reduced02 dx50 near downstream (2 junctions)
            ! integer :: iet(10) = (/   61,60,62,   71,    72,   73,     75,74,76,     85 /)
            ! integer :: ift(7) = (/52,           53,   62,    63,    64,          65,    74/)

            ! integer :: iet(6) = (/ 2049 ,       2050,         2022, 2019, 2021,    2023         /)
            ! integer :: ift(4) = (/      1920         ,1892,                   1891,    1893 /)
        !%------------------------------------------------------------------
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

            ! print *, ' '
            ! print *, '==========================================================='
            ! print *, 'elemR(1,er_Topwidth)', elemR(1,er_TopWidth)
            ! print *, ' '

           ! print *, 'plan area ',elemSR(4,esr_Storage_Plan_Area)
         ! print *, 'dt = ',dt
       ! return

            ! print *, reverseKey(elemI(307,ei_elementType))
            ! print *, reverseKey(elemI(320,ei_elementType))
            ! print *, elemI(307,ei_Mface_uL), elemI(307,ei_Mface_dL)
            ! print *, elemI(320,ei_Mface_uL), elemI(320,ei_Mface_dL)


            ! print *, faceI(elem)

            ! stop 2398743

       
          ! if (setting%Time%Step < 20733) return



            !   if (setting%Time%Step >1) then

            !     stop  6098723
            !   end if

   !% --- USEFUL HEADER -----------------------------------------------
            print *, ' '
            write(*,"(A,A, e12.5)") ' ',trim(inputstring)     
            write(*,"(A,i7,A, f12.5, A, f12.5, A)") '        step = ',setting%Time%Step ,&
            '; dt = ',setting%Time%Hydraulics%Dt,&
            '; time = ',setting%Time%Now / 60.d0, ' min'        

           
            

            ! if (elemR(46,er_Head) > 24) then 
            !    stop 5509872
            ! end if

            ! print *, elemR(102,er_Flowrate), elemR(101,er_Flowrate), elemR(103,er_Flowrate)
            ! print *, elemR(102,er_Depth), elemR(101,er_Depth), elemR(103,er_Depth)
            ! print *, elemR(91:110,er_Flowrate)
            ! return 

   !% --- SETUP  -----------------------------------------------
   !% branches and elements   

         ! print *, 'branches and elements'
         !       do ii= 1,N_elem(1)
         !          if (elemI(ii,ei_elementType) == CC) then
         !             if (faceI(fup(ii),fi_BCType) == BCup) then 
         !                write(*,"(A, i5, A)"), 'fbcU', fup(ii),       '             '//trim(reverseKey(faceI(fup(ii),fi_BCType))) 
         !             end if
         !                write(*,"(A, i5, A)"), 'fup ', fup(ii),       '             '//trim(reverseKey(faceI(fup(ii),fi_BCType)))         
         !                write(*,"(A, i5,A, f10.2)"), 'e   ', ii,            '             '//trim(reverseKey(elemI(    ii ,ei_elementType)))//'      '//trim(link%Names(elemI(ii,ei_link_Gidx_SWMM))%str), elemR(ii,er_Length)
         !                write(*,"(A, i5,A)"), 'fdn ', fdn(ii),       '             '//trim(reverseKey(faceI(fdn(ii),fi_BCType)))
         !             if (faceI(fup(ii),fi_BCType) == BCdn) then
         !                write(*,"(A, i5, A)"), 'fbcD', fdn(ii),       '             '//trim(reverseKey(faceI(fdn(ii),fi_BCType)))
         !             end if 
         !          elseif (elemI(ii,ei_elementType) == JM) then
         !             !% upstream branches
         !             do jj=1,max_branch_per_node,2
         !                if  (elemSI(ii+jj,esi_JunctionBranch_Exists) == oneI) then
         !                   write(*,"(A, i5,  A)"),     'fup ', fup(ii+jj),  '             '//trim(reverseKey(faceI(fup(ii+jj),fi_BCType)))
         !                   write(*,"(A, i5, A, i5)"), 'eJB ', ii+jj,      '             JB ',   elemSI(ii+jj,esi_JunctionBranch_IsUpstream)
         !                end if
         !             end do
         !             !%  JM
         !                   write(*,"(A, i5, A)"),      'eJM ', ii,  ' '//trim(reverseKey(elemI(ii,ei_elementType)))//'      '//trim(node%Names(elemI(ii,ei_node_Gidx_SWMM))%str)
         !             !% downstream branches
         !             do jj=2,max_branch_per_node,2
         !                if  (elemSI(ii+jj,esi_JunctionBranch_Exists) == oneI) then
         !                   write(*,"(A, i5,  A, i5)"), 'eJB ', ii+jj     ,  '            JB ', elemSI(ii+jj,esi_JunctionBranch_IsUpstream)
         !                   write(*,"(A, i5, A)"),     'fdn ', fdn(ii+jj),  '               '//trim(reverseKey(faceI(fdn(ii+jj),fi_BCType)))
         !                end if
         !             end do
         !          elseif ((elemI(ii,ei_elementType) == JB) ) then
         !             ! skip
         !          elseif ((elemI(ii,ei_elementType) == weir)) then 
         !             write(*,"(A, i5, A)"), 'fup ', fup(ii),       '             '//trim(reverseKey(faceI(fup(ii),fi_BCType)))         
         !             write(*,"(A, i5,A)"), 'e   ', ii,            '             '//trim(reverseKey(elemI(    ii ,ei_elementType)))//'      '//trim(link%Names(elemI(ii,ei_link_Gidx_SWMM))%str)
         !             write(*,"(A, i5,A)"), 'fdn ', fdn(ii),       '             '//trim(reverseKey(faceI(fdn(ii),fi_BCType)))
         !          else 
         !             print *, 'unexpected element type '
         !             print *, trim(reverseKey(elemI(ii,ei_elementType)))
         !             stop 77873
         !          end if
         !        end do
         !            stop 239874
         !           return

               ! print *, ' '
               ! print *, 'link,node names'
               ! do ii=1,N_elem(1)
               !    if (elemI(ii,ei_link_Gidx_SWMM) < nullvalueI) then 
               !       print *, ii, elemI(ii,ei_link_Gidx_SWMM), trim(link%Names(elemI(ii,ei_link_Gidx_SWMM))%str), ' islink'
               !    elseif (elemI(ii,ei_node_Gidx_SWMM) < nullvalueI) then
               !       print *, ii, elemI(ii,ei_node_Gidx_SWMM), trim(node%Names(elemI(ii,ei_node_Gidx_SWMM))%str), ' isnode'
               !    end if
               ! end do

               ! stop 239874
!%
!%  6 x 4 junction with two downstream
!%
            ! write(*,"(A9,21A11)"),' ','elem','faceU','faceD','elem','faceU','faceD',' jb ',' JM ',' jb ','faceU','faceD',' jb ','faceU','faceD'
   
         !  write(*,"(A,31e11.4)") '  Head    '  ,          &
         !    elemR(iet(1),er_Head), &
         !       faceR(ift(1),fr_Head_u), &
         !       faceR(ift(1),fr_Head_d), &
         !    elemR(iet(2),er_Head), &
         !       faceR(ift(2),fr_Head_u), &
         !       faceR(ift(2),fr_Head_d), &
         !    elemR(iet(3),er_Head), &
         !    elemR(iet(4),er_Head), &
         !    elemR(iet(5),er_Head), &
         !       faceR(ift(3),fr_Head_u), &
         !       faceR(ift(3),fr_Head_d), &
         !    elemR(iet(6),er_Head),&
         !       faceR(ift(4),fr_Head_u), &
         !       faceR(ift(4),fr_Head_d)

            ! write(*,"(A,31e11.4)") '  Depth   '  ,          &
            ! elemR(iet(1),er_Depth), &
            !    faceR(ift(1),fr_Depth_u), &
            !    faceR(ift(1),fr_Depth_d), &
            ! elemR(iet(2),er_Depth), &
            !    faceR(ift(2),fr_Depth_u), &
            !    faceR(ift(2),fr_Depth_d), &
            ! elemR(iet(3),er_Depth), &
            ! elemR(iet(4),er_Depth), &
            ! elemR(iet(5),er_Depth), &
            !    faceR(ift(3),fr_Depth_u), &
            !    faceR(ift(3),fr_Depth_d), &
            ! elemR(iet(6),er_Depth),&
            !    faceR(ift(4),fr_Depth_u), &
            !    faceR(ift(4),fr_Depth_d)

               ! write(*,"(A,31e11.4)") '  Area    '  ,          &
               ! elemR(iet(1),er_Area), &
               !    faceR(ift(1),fr_Area_u), &
               !    faceR(ift(1),fr_Area_d), &
               ! elemR(iet(2),er_Area), &
               !    faceR(ift(2),fr_Area_u), &
               !    faceR(ift(2),fr_Area_d), &
               ! elemR(iet(3),er_Area), &
               ! elemR(iet(4),er_Area), &
               ! elemR(iet(5),er_Area), &
               !    faceR(ift(3),fr_Area_u), &
               !    faceR(ift(3),fr_Area_d), &
               ! elemR(iet(6),er_Area),&
               !    faceR(ift(4),fr_Area_u), &
               !    faceR(ift(4),fr_Area_d)

            ! write(*,"(A,31e11.4)") '  Flow    '  ,          &
            !    elemR(iet(1),er_Flowrate), &
            !       faceR(ift(1),fr_Flowrate), &
            !       faceR(ift(1),fr_Flowrate), &
            !    elemR(iet(2),er_Flowrate), &
            !       faceR(ift(2),fr_Flowrate), &
            !       faceR(ift(2),fr_Flowrate), &
            !    elemR(iet(3),er_Flowrate), &
            !    elemR(iet(4),er_Flowrate), &
            !    elemR(iet(5),er_Flowrate), &
            !       faceR(ift(3),fr_Flowrate), &
            !       faceR(ift(3),fr_Flowrate), &
            !    elemR(iet(6),er_Flowrate),&
            !       faceR(ift(4),fr_Flowrate), &
            !       faceR(ift(4),fr_Flowrate)

            !       write(*,"(A,31e11.4)") '  Vel     '  ,          &
            !       elemR(iet(1),er_Velocity), &
            !          faceR(ift(1),fr_Velocity_u), &
            !          faceR(ift(1),fr_Velocity_d), &
            !       elemR(iet(2),er_Velocity), &
            !          faceR(ift(2),fr_Velocity_u), &
            !          faceR(ift(2),fr_Velocity_d), &
            !       elemR(iet(3),er_Velocity), &
            !       elemR(iet(4),er_Velocity), &
            !       elemR(iet(5),er_Velocity), &
            !          faceR(ift(3),fr_Velocity_u), &
            !          faceR(ift(3),fr_Velocity_d), &
            !       elemR(iet(6),er_Velocity),&
            !          faceR(ift(4),fr_Velocity_u), &
            !          faceR(ift(4),fr_Velocity_d)

            ! return

!%
!%  10 x 7 starting with 1 face, jbJM/jb, 3 elem, jb/jm/jb, 1 elem, 1 face
!%

         !  write(*,"(A9,25A11)"),' ','faceD',' jb ',' JM ',' jb ','faceU','faceD','elem','faceU','faceD','elem','faceU','faceD','elem','faceU','faceD',' jb ',' JM ',' jb ','faceU','faceD','elem','faceU'

         ! write(*,"(A,31f11.4)") '  Head    '  ,          &
         !             faceR(ift(1),fr_Head_d), &
         !          elemR(iet(1),er_Head), &
         !          elemR(iet(2),er_Head), &
         !          elemR(iet(3),er_Head), &
         !             faceR(ift(2),fr_Head_u), &
         !             faceR(ift(2),fr_Head_d), &
         !          elemR(iet(4),er_Head), &
         !             faceR(ift(3),fr_Head_u), &
         !             faceR(ift(3),fr_Head_d), &   
         !          elemR(iet(5),er_Head), &
         !             faceR(ift(4),fr_Head_u), &
         !             faceR(ift(4),fr_Head_d), &
         !          elemR(iet(6),er_Head), &
         !             faceR(ift(5),fr_Head_u), &
         !             faceR(ift(5),fr_Head_d), &
         !          elemR(iet(7),er_Head), &
         !          elemR(iet(8),er_Head), &
         !          elemR(iet(9),er_Head), &
         !             faceR(ift(6),fr_Head_u), &
         !             faceR(ift(6),fr_Head_d), &
         !          elemR(iet(10),er_Head), &
         !             faceR(ift(7),fr_Head_u)
         
            ! write(*,"(A,31f11.4)") '  Zbtm    '  ,          &
            !          faceR(ift(1),fr_Zbottom), &
            !       elemR(iet(1),er_Zbottom), &
            !       elemR(iet(2),er_Zbottom), &
            !       elemR(iet(3),er_Zbottom), &
            !          faceR(ift(2),fr_Zbottom), &
            !          faceR(ift(2),fr_Zbottom), &
            !       elemR(iet(4),er_Zbottom), &
            !          faceR(ift(3),fr_Zbottom), &
            !          faceR(ift(3),fr_Zbottom), &   
            !       elemR(iet(5),er_Zbottom), &
            !          faceR(ift(4),fr_Zbottom), &
            !          faceR(ift(4),fr_Zbottom), &
            !       elemR(iet(6),er_Zbottom), &
            !          faceR(ift(5),fr_Zbottom), &
            !          faceR(ift(5),fr_Zbottom), &
            !       elemR(iet(7),er_Zbottom), &
            !       elemR(iet(8),er_Zbottom), &
            !       elemR(iet(9),er_Zbottom), &
            !          faceR(ift(6),fr_Zbottom), &
            !          faceR(ift(6),fr_Zbottom), &
            !       elemR(iet(10),er_Zbottom), &
            !          faceR(ift(7),fr_Zbottom)


               !  write(*,"(A,31f11.4)") '  Depth   '  ,          &
               !       faceR(ift(1),fr_Depth_d), &
               !    elemR(iet(1),er_Depth), &
               !    elemR(iet(2),er_Depth), &
               !    elemR(iet(3),er_Depth), &
               !       faceR(ift(2),fr_Depth_u), &
               !       faceR(ift(2),fr_Depth_d), &
               !    elemR(iet(4),er_Depth), &
               !       faceR(ift(3),fr_Depth_u), &
               !       faceR(ift(3),fr_Depth_d), &   
               !    elemR(iet(5),er_Depth), &
               !       faceR(ift(4),fr_Depth_u), &
               !       faceR(ift(4),fr_Depth_d), &
               !    elemR(iet(6),er_Depth), &
               !       faceR(ift(5),fr_Depth_u), &
               !       faceR(ift(5),fr_Depth_d), &
               !    elemR(iet(7),er_Depth), &
               !    elemR(iet(8),er_Depth), &
               !    elemR(iet(9),er_Depth), &
               !       faceR(ift(6),fr_Depth_u), &
               !       faceR(ift(6),fr_Depth_d), &
               !    elemR(iet(10),er_Depth), &
               !       faceR(ift(7),fr_Depth_u)

               !       print *, ' '

                  ! write(*,"(A,31e11.3)") '  Flow    '  ,          &
                  !    faceR(ift(1),fr_Flowrate), &
                  ! elemR(iet(1),er_Flowrate), &
                  ! elemR(iet(2),er_Flowrate), &
                  ! elemR(iet(3),er_Flowrate), &
                  !    faceR(ift(2),fr_Flowrate), &
                  !    faceR(ift(2),fr_Flowrate), &
                  ! elemR(iet(4),er_Flowrate), &
                  !    faceR(ift(3),fr_Flowrate), &
                  !    faceR(ift(3),fr_Flowrate), &   
                  ! elemR(iet(5),er_Flowrate), &
                  !    faceR(ift(4),fr_Flowrate), &
                  !    faceR(ift(4),fr_Flowrate), &
                  ! elemR(iet(6),er_Flowrate), &
                  !    faceR(ift(5),fr_Flowrate), &
                  !    faceR(ift(5),fr_Flowrate), &
                  ! elemR(iet(7),er_Flowrate), &
                  ! elemR(iet(8),er_Flowrate), &
                  ! elemR(iet(9),er_Flowrate), &
                  !    faceR(ift(6),fr_Flowrate), &
                  !    faceR(ift(6),fr_Flowrate), &
                  ! elemR(iet(10),er_Flowrate), &
                  !    faceR(ift(7),fr_Flowrate)

                  !    print *, ' '

                  !    write(*,"(A,31f11.4)") '  Vol     '  ,          &
                  !    0.d0, &
                  ! elemR(iet(1),er_Volume), &
                  ! elemR(iet(2),er_Volume), &
                  ! elemR(iet(3),er_Volume), &
                  !    0.d0, &
                  !    0.d0, &
                  ! elemR(iet(4),er_Volume), &
                  !    0.d0, &
                  !    0.d0, &   
                  ! elemR(iet(5),er_Volume), &
                  !    0.d0, &
                  !    0.d0, &
                  ! elemR(iet(6),er_Volume), &
                  !    0.d0, &
                  !    0.d0, &
                  ! elemR(iet(7),er_Volume), &
                  ! elemR(iet(8),er_Volume), &
                  ! elemR(iet(9),er_Volume), &
                  !    0.d0, &
                  !    0.d0, &
                  ! elemR(iet(10),er_Volume), &
                  !    0.d0

                  !    print *, ' '
                  !    return

!%
!%  8 x 7 starting with 1 face, 1 elem, J, 4 elem
!% 
      !  write(*,"(A9,21A11)"),' ','faceD','elem','faceU','faceD',' jb ',' JM ',' jb ','faceU','faceD','elem','faceU','faceD','elem','faceU','faceD','elem','faceU','faceD','elem','faceU','faceD'


               ! write(*,"(A,31f11.4)") '  Head    '  ,          &
               !       faceR(ift(1),fr_Head_d), &
               !    elemR(iet(1),er_Head), &
               !       faceR(ift(2),fr_Head_u), &
               !       faceR(ift(2),fr_Head_d), &
               !    elemR(iet(2),er_Head), &
               !    elemR(iet(3),er_Head), &
               !    elemR(iet(4),er_Head), &
               !       faceR(ift(3),fr_Head_u), &
               !       faceR(ift(3),fr_Head_d), &
               !    elemR(iet(5),er_Head), &
               !       faceR(ift(4),fr_Head_u), &
               !       faceR(ift(4),fr_Head_d), &
               !    elemR(iet(6),er_Head), &
               !       faceR(ift(5),fr_Head_u), &
               !       faceR(ift(5),fr_Head_d), &
               !    elemR(iet(7),er_Head), &
               !       faceR(ift(6),fr_Head_u), &
               !       faceR(ift(6),fr_Head_d), &
               !    elemR(iet(8),er_Head), &
                     ! faceR(ift(7),fr_Head_u)

               ! write(*,"(A,31f11.4)") '  Zbtm    '  ,          &
               !       faceR(ift(1),fr_Zbottom), &
               !    elemR(iet(1),er_Zbottom), &
               !       faceR(ift(2),fr_Zbottom), &
               !       faceR(ift(2),fr_Zbottom), &
               !    elemR(iet(2),er_Zbottom), &
               !    elemR(iet(3),er_Zbottom), &
               !    elemR(iet(4),er_Zbottom), &
               !       faceR(ift(3),fr_Zbottom), &
               !       faceR(ift(3),fr_Zbottom), &
               !    elemR(iet(5),er_Zbottom), &
               !       faceR(ift(4),fr_Zbottom), &
               !       faceR(ift(4),fr_Zbottom), &
               !    elemR(iet(6),er_Zbottom), &
               !       faceR(ift(5),fr_Zbottom), &
               !       faceR(ift(5),fr_Zbottom), &
               !    elemR(iet(7),er_Zbottom), &
               !       faceR(ift(6),fr_Zbottom), &
               !       faceR(ift(6),fr_Zbottom), &
               !    elemR(iet(8),er_Zbottom), &
               !       faceR(ift(7),fr_Zbottom)

                  !    write(*,"(A,31f11.4)") '  Depth   '  ,          &
                  !    faceR(ift(1),fr_Depth_d), &
                  ! elemR(iet(1),er_Depth), &
                  !    faceR(ift(2),fr_Depth_u), &
                  !    faceR(ift(2),fr_Depth_d), &
                  ! elemR(iet(2),er_Depth), &
                  ! elemR(iet(3),er_Depth), &
                  ! elemR(iet(4),er_Depth), &
                  !    faceR(ift(3),fr_Depth_u), &
                  !    faceR(ift(3),fr_Depth_d), &
                  ! elemR(iet(5),er_Depth), &
                  !    faceR(ift(4),fr_Depth_u), &
                  !    faceR(ift(4),fr_Depth_d), &
                  ! elemR(iet(6),er_Depth), &
                  !    faceR(ift(5),fr_Depth_u), &
                  !    faceR(ift(5),fr_Depth_d), &
                  ! elemR(iet(7),er_Depth), &
                  !    faceR(ift(6),fr_Depth_u), &
                  !    faceR(ift(6),fr_Depth_d), &
                  ! elemR(iet(8),er_Depth), &
                  !    faceR(ift(7),fr_Depth_u)

                  ! print *, ' '

                  !    write(*,"(A,31f11.4)") '  Flow    '  ,          &
                  !    faceR(ift(1),fr_Flowrate), &
                  ! elemR(iet(1),er_Flowrate), &
                  !    faceR(ift(2),fr_Flowrate), &
                  !    faceR(ift(2),fr_Flowrate), &
                  ! elemR(iet(2),er_Flowrate), &
                  ! elemR(iet(3),er_Flowrate), &
                  ! elemR(iet(4),er_Flowrate), &
                  !    faceR(ift(3),fr_Flowrate), &
                  !    faceR(ift(3),fr_Flowrate), &
                  ! elemR(iet(5),er_Flowrate), &
                  !    faceR(ift(4),fr_Flowrate), &
                  !    faceR(ift(4),fr_Flowrate), &
                  ! elemR(iet(6),er_Flowrate), &
                  !    faceR(ift(5),fr_Flowrate), &
                  !    faceR(ift(5),fr_Flowrate), &
                  ! elemR(iet(7),er_Flowrate), &
                  !    faceR(ift(6),fr_Flowrate), &
                  !    faceR(ift(6),fr_Flowrate), &
                  ! elemR(iet(8),er_Flowrate), &
                  !    faceR(ift(7),fr_Flowrate)

                  !    return
!%
  
!% 
!% 8 x 7 starting with 1 face, 2 elem, J , 3 elem
!% 
         ! write(*,"(A9,21A11)"),' ','faceD','elem','faceU','faceD','elem','faceU','faceD',' jb ',' JM ',' jb ','faceU','faceD','elem','faceU','faceD','elem','faceU','faceD','elem','faceU'


               ! write(*,"(A,31f11.4)") '  Head    '  ,          &
               !       faceR(ift(1),fr_Head_d), &
               !    elemR(iet(1),er_Head), &
               !       faceR(ift(2),fr_Head_u), &
               !       faceR(ift(2),fr_Head_d), &
               !    elemR(iet(2),er_Head), &
               !       faceR(ift(3),fr_Head_u), &
               !       faceR(ift(3),fr_Head_d), &   
               !    elemR(iet(3),er_Head), &
               !    elemR(iet(4),er_Head), &
               !    elemR(iet(5),er_Head), &
               !       faceR(ift(4),fr_Head_u), &
               !       faceR(ift(4),fr_Head_d), &
               !    elemR(iet(6),er_Head), &
               !       faceR(ift(5),fr_Head_u), &
               !       faceR(ift(5),fr_Head_d), &
               !    elemR(iet(7),er_Head), &
               !       faceR(ift(6),fr_Head_u), &
               !       faceR(ift(6),fr_Head_d), &
               !    elemR(iet(8),er_Head), &
               !       faceR(ift(7),fr_Head_u)

                  ! write(*,"(A,31f11.4)") '  Zbtm    '  ,          &
                  !    faceR(ift(1),fr_Zbottom), &
                  ! elemR(iet(1),er_Zbottom), &
                  !    faceR(ift(2),fr_Zbottom), &
                  !    faceR(ift(2),fr_Zbottom), &
                  ! elemR(iet(2),er_Zbottom), &
                  !    faceR(ift(3),fr_Zbottom), &
                  !    faceR(ift(3),fr_Zbottom), &   
                  ! elemR(iet(3),er_Zbottom), &
                  ! elemR(iet(4),er_Zbottom), &
                  ! elemR(iet(5),er_Zbottom), &
                  !    faceR(ift(4),fr_Zbottom), &
                  !    faceR(ift(4),fr_Zbottom), &
                  ! elemR(iet(6),er_Zbottom), &
                  !    faceR(ift(5),fr_Zbottom), &
                  !    faceR(ift(5),fr_Zbottom), &
                  ! elemR(iet(7),er_Zbottom), &
                  !    faceR(ift(6),fr_Zbottom), &
                  !    faceR(ift(6),fr_Zbottom), &
                  ! elemR(iet(8),er_Zbottom), &
                  !    faceR(ift(7),fr_Zbottom)

                  ! write(*,"(A,31f11.4)") ' Depth    '  ,          &
                  !    faceR(ift(1),fr_Depth_d), &
                  ! elemR(iet(1),er_Depth), &
                  !    faceR(ift(2),fr_Depth_u), &
                  !    faceR(ift(2),fr_Depth_d), &
                  ! elemR(iet(2),er_Depth), &
                  !    faceR(ift(3),fr_Depth_u), &
                  !    faceR(ift(3),fr_Depth_d), &   
                  ! elemR(iet(3),er_Depth), &
                  ! elemR(iet(4),er_Depth), &
                  ! elemR(iet(5),er_Depth), &
                  !    faceR(ift(4),fr_Depth_u), &
                  !    faceR(ift(4),fr_Depth_d), &
                  ! elemR(iet(6),er_Depth), &
                  !    faceR(ift(5),fr_Depth_u), &
                  !    faceR(ift(5),fr_Depth_d), &
                  ! elemR(iet(7),er_Depth), &
                  !    faceR(ift(6),fr_Depth_u), &
                  !    faceR(ift(6),fr_Depth_d), &
                  ! elemR(iet(8),er_Depth), &
                  !    faceR(ift(7),fr_Depth_u)

               !       print *, ' '
            
               ! write(*,"(A,31e11.3)") '  Flow    '  ,          &
               !       faceR(ift(1),fr_Flowrate), &
               !    elemR(iet(1),er_Flowrate), &
               !       faceR(ift(2),fr_Flowrate), &
               !       faceR(ift(2),fr_Flowrate), &
               !    elemR(iet(2),er_Flowrate), &
               !       faceR(ift(3),fr_Flowrate), &
               !       faceR(ift(3),fr_Flowrate), &   
               !    elemR(iet(3),er_Flowrate), &
               !    elemR(iet(4),er_Flowrate), &
               !    elemR(iet(5),er_Flowrate), &
               !       faceR(ift(4),fr_Flowrate), &
               !       faceR(ift(4),fr_Flowrate), &
               !    elemR(iet(6),er_Flowrate), &
               !       faceR(ift(5),fr_Flowrate), &
               !       faceR(ift(5),fr_Flowrate), &
               !    elemR(iet(7),er_Flowrate), &
               !       faceR(ift(6),fr_Flowrate), &
               !       faceR(ift(6),fr_Flowrate), &
               !    elemR(iet(8),er_Flowrate), &
               !       faceR(ift(7),fr_Flowrate)

               !       print *, ' '

               !       return
!%            
!%  8 x 7 starting with 1 face, 4 elem
!% 

      ! write(*,"(A9,21A11)"),' ','faceD','elem','faceU','faceD','elem','faceU','faceD','elem','faceU','faceD','elem','faceU','faceD',' jb ',' JM ',' jb ','faceU','faceD','elem','faceU'

         !  write(*,"(A,31f11.4)") '  Head    '  ,          &
         !          faceR(ift(1),fr_Head_d), &
         !       elemR(iet(1),er_Head), &
         !          faceR(ift(2),fr_Head_u), &
         !          faceR(ift(2),fr_Head_d), &
         !       elemR(iet(2),er_Head), &
         !          faceR(ift(3),fr_Head_u), &
         !          faceR(ift(3),fr_Head_d), &
         !       elemR(iet(3),er_Head), &
         !          faceR(ift(4),fr_Head_u), &
         !          faceR(ift(4),fr_Head_d), &
         !       elemR(iet(4),er_Head), &
         !          faceR(ift(5),fr_Head_u), &
         !          faceR(ift(5),fr_Head_d), &
         !       elemR(iet(5),er_Head), &
         !       elemR(iet(6),er_Head), &
         !       elemR(iet(7),er_Head), &
         !          faceR(ift(6),fr_Head_u), &
         !          faceR(ift(6),fr_Head_d), &
         !       elemR(iet(8),er_Head), &
         !          faceR(ift(7),fr_Head_u)

            ! write(*,"(A,31f11.4)") '  Zbtm    '  ,          &
            !       faceR(ift(1),fr_Zbottom), &
            !    elemR(iet(1),er_Zbottom), &
            !       faceR(ift(2),fr_Zbottom), &
            !       faceR(ift(2),fr_Zbottom), &
            !    elemR(iet(2),er_Zbottom), &
            !       faceR(ift(3),fr_Zbottom), &
            !       faceR(ift(3),fr_Zbottom), &
            !    elemR(iet(3),er_Head), &
            !       faceR(ift(4),fr_Zbottom), &
            !       faceR(ift(4),fr_Zbottom), &
            !    elemR(iet(4),er_Zbottom), &
            !       faceR(ift(5),fr_Zbottom), &
            !       faceR(ift(5),fr_Zbottom), &
            !    elemR(iet(5),er_Zbottom), &
            !    elemR(iet(6),er_Zbottom), &
            !    elemR(iet(7),er_Zbottom), &
            !       faceR(ift(6),fr_Zbottom), &
            !       faceR(ift(6),fr_Zbottom), &
            !    elemR(iet(8),er_Zbottom), &
            !       faceR(ift(7),fr_Zbottom)      

            ! write(*,"(A,31e11.3)") '  Depth   '  ,          &
            !       faceR(ift(1),fr_Depth_d), &
            !    elemR(iet(1),er_Depth), &
            !       faceR(ift(2),fr_Depth_u), &
            !       faceR(ift(2),fr_Depth_d), &
            !    elemR(iet(2),er_Depth), &
            !       faceR(ift(3),fr_Depth_u), &
            !       faceR(ift(3),fr_Depth_d), &
            !    elemR(iet(3),er_Depth), &
            !       faceR(ift(4),fr_Depth_u), &
            !       faceR(ift(4),fr_Depth_d), &
            !    elemR(iet(4),er_Depth), &
            !       faceR(ift(5),fr_Depth_u), &
            !       faceR(ift(5),fr_Depth_d), &
            !    elemR(iet(5),er_Depth), &
            !    elemR(iet(6),er_Depth), &
            !    elemR(iet(7),er_Depth), &
            !       faceR(ift(6),fr_Depth_u), &
            !       faceR(ift(6),fr_Depth_d), &
            !    elemR(iet(8),er_Depth), &
            !       faceR(ift(7),fr_Depth_u)

            !    print *, ' '

            ! write(*,"(A,31e11.3)") '   Flow   '  ,          &
            !       faceR(ift(1),fr_Flowrate), &
            !    elemR(iet(1),er_Flowrate), &
            !       faceR(ift(2),fr_Flowrate), &
            !       faceR(ift(2),fr_Flowrate), &
            !    elemR(iet(2),er_Flowrate), &
            !       faceR(ift(3),fr_Flowrate), &
            !       faceR(ift(3),fr_Flowrate), &
            !    elemR(iet(3),er_Flowrate), &
            !       faceR(ift(4),fr_Flowrate), &
            !       faceR(ift(4),fr_Flowrate), &
            !    elemR(iet(4),er_Flowrate), &
            !       faceR(ift(5),fr_Flowrate), &
            !       faceR(ift(5),fr_Flowrate), &
            !    elemR(iet(5),er_Flowrate), &
            !    elemR(iet(6),er_Flowrate), &
            !    elemR(iet(7),er_Flowrate), &
            !       faceR(ift(6),fr_Flowrate), &
            !       faceR(ift(6),fr_Flowrate), &
            !    elemR(iet(8),er_Flowrate), &
            !       faceR(ift(7),fr_Flowrate)

            !       print *, ' '

               ! write(*,"(A,31e11.3)") '  Veloc   '  ,          &
               !    faceR(ift(1),fr_Velocity_d), &
               ! elemR(iet(1),er_Velocity), &
               !    faceR(ift(2),fr_Velocity_u), &
               !    faceR(ift(2),fr_Velocity_d), &
               ! elemR(iet(2),er_Velocity), &
               !    faceR(ift(3),fr_Velocity_u), &
               !    faceR(ift(3),fr_Velocity_d), &
               ! elemR(iet(3),er_Velocity), &
               !    faceR(ift(4),fr_Velocity_u), &
               !    faceR(ift(4),fr_Velocity_d), &
               ! elemR(iet(4),er_Velocity), &
               !    faceR(ift(5),fr_Velocity_u), &
               !    faceR(ift(5),fr_Velocity_d), &
               ! elemR(iet(5),er_Velocity), &
               ! elemR(iet(6),er_Velocity), &
               ! elemR(iet(7),er_Velocity), &
               !    faceR(ift(6),fr_Velocity_u), &
               !    faceR(ift(6),fr_Velocity_d), &
               ! elemR(iet(8),er_Velocity), &
               !    faceR(ift(7),fr_Velocity_u)

               !    print *, ' '

               ! write(*,"(A,31e11.3)") '   Volu   '  ,          &
               ! 0.d0, &
               !    elemR(iet(1),er_Volume), &
               !    0.d0, &
               !       0.d0, &
               !    elemR(iet(2),er_Volume), &
               !    0.d0, &
               !       0.d0, &
               !    elemR(iet(3),er_Volume), &
               !    0.d0, &
               !       0.d0, &
               !    elemR(iet(4),er_Volume), &
               !    0.d0, &
               !       0.d0, &
               !    elemR(iet(5),er_Volume), &
               !    elemR(iet(6),er_Volume), &
               !    elemR(iet(7),er_Volume), &
               !    0.d0, &
               !       0.d0, &
               !    elemR(iet(8),er_Volume), &
               !       0.d0

               !       return
!%
!% junction with two upstream branches
!% 
   ! write(*,"(A9,20A11)"),' ','elem','faceU','faceD','jb  ','xxxx','elem','faceU','faceD','jb  ','xxxx','JM','xxxx','jb','faceU','faceD','elem'
           !% junction with two upstream branches
         !    write(*,"(A,4e11.3,A,4e11.3,A,e11.3,A,4e11.3)") '  Head    '  ,          &
         !    elemR(iet(1),er_Head), &
         !       faceR(ift(1),fr_Head_u), &
         !       faceR(ift(1),fr_Head_d), &
         !    elemR(iet(2),er_Head), &
         !    '    -----    ',&
         !    elemR(iet(3),er_Head), &
         !       faceR(ift(2),fr_Head_u), &
         !       faceR(ift(2),fr_Head_d), &
         !    elemR(iet(4),er_Head), &
         !    '    -----    ',&
         !    elemR(iet(5),er_Head), &
         !    '    -----    ',&
         !    elemR(iet(6),er_Head), &
         !       faceR(ift(3),fr_Head_u), &
         !       faceR(ift(3),fr_Head_d), &
         !    elemR(iet(7),er_Head)

         !    write(*,"(A,4f11.3,A,4f11.3,A,f11.3,A,4f11.3)") '  Zbott   '  ,          &
         !    elemR(iet(1),er_Zbottom), &
         !       faceR(ift(1),fr_Zbottom), &
         !       faceR(ift(1),fr_Zbottom), &
         !    elemR(iet(2),er_Zbottom), &
         !    '    -----    ',&
         !    elemR(iet(3),er_Zbottom), &
         !       faceR(ift(2),fr_Zbottom), &
         !       faceR(ift(2),fr_Zbottom), &
         !    elemR(iet(4),er_Zbottom), &
         !    '    -----    ',&
         !    elemR(iet(5),er_Zbottom), &
         !    '    -----    ',&
         !    elemR(iet(6),er_Zbottom), &
         !       faceR(ift(3),fr_Zbottom), &
         !       faceR(ift(3),fr_Zbottom), &
         !    elemR(iet(7),er_Zbottom)

         !    write(*,"(A,4f11.3,A,4f11.3,A,f11.3,A,4f11.3)") '  Depth   '  ,          &
         !    elemR(iet(1),er_Depth), &
         !       faceR(ift(1),fr_Depth_u), &
         !       faceR(ift(1),fr_Depth_d), &
         !    elemR(iet(2),er_Depth), &
         !    '    -----    ',&
         !    elemR(iet(3),er_Depth), &
         !       faceR(ift(2),fr_Depth_u), &
         !       faceR(ift(2),fr_Depth_d), &
         !    elemR(iet(4),er_Depth), &
         !    '    -----    ',&
         !    elemR(iet(5),er_Depth), &
         !    '    -----    ',&
         !    elemR(iet(6),er_Depth), &
         !       faceR(ift(3),fr_Depth_u), &
         !       faceR(ift(3),fr_Depth_d), &
         !    elemR(iet(7),er_Depth)

         ! return

!%
!%  8/5  3 elements, J, 2 elements
!%
   ! write(*,"(A9,20A11)"),' ','elem','faceU','faceD','elem','faceU','faceD','elem','faceU','faceD','jb','JM','jb','faceU','faceD','elem','faceU','faceD','elem'
           !NJ2 with 1 junction     8 elm, 5 face looking at upstream
     
         ! write(*,"(A,31f11.4)") '  Head    '  ,          &
         !    elemR(iet(1),er_Head), &
         !       faceR(ift(1),fr_Head_u), &
         !       faceR(ift(1),fr_Head_d), &
         !    elemR(iet(2),er_Head), &
         !       faceR(ift(2),fr_Head_u), &
         !       faceR(ift(2),fr_Head_d), &
         !    elemR(iet(3),er_Head), &
         !       faceR(ift(3),fr_Head_u), &
         !       faceR(ift(3),fr_Head_d), &
         !    elemR(iet(4),er_Head), &
         !    elemR(iet(5),er_Head), &
         !    elemR(iet(6),er_Head), &
         !       faceR(ift(4),fr_Head_u), &
         !       faceR(ift(4),fr_Head_d), &
         !    elemR(iet(7),er_Head),&
         !       faceR(ift(5),fr_Head_u), &
         !       faceR(ift(5),fr_Head_d), &
         !    elemR(iet(8),er_Head)

            ! write(*,"(A,31f11.4)") '  SlotD   '  ,          &
            ! elemR(iet(1),er_SlotDepth), &
            ! 0.d0, &
            !    0.d0, &
            ! elemR(iet(2),er_SlotDepth), &
            ! 0.d0, &
            !    0.d0, &
            ! elemR(iet(3),er_SlotDepth), &
            !    0.d0, &
            !    0.d0, &
            ! elemR(iet(4),er_SlotDepth), &
            ! elemR(iet(5),er_SlotDepth), &
            ! elemR(iet(6),er_SlotDepth), &
            ! 0.d0, &
            !    0.d0, &
            ! elemR(iet(7),er_SlotDepth),&
            ! 0.d0, &
            !    0.d0, &
            ! elemR(iet(8),er_SlotDepth)

            ! write(*,"(A,31f11.4)") '  SlotA   '  ,          &
            ! elemR(iet(1),er_SlotArea), &
            ! 0.d0, &
            !    0.d0, &
            ! elemR(iet(2),er_SlotArea), &
            ! 0.d0, &
            !    0.d0, &
            ! elemR(iet(3),er_SlotArea), &
            !    0.d0, &
            !    0.d0, &
            ! elemR(iet(4),er_SlotArea), &
            ! elemR(iet(5),er_SlotArea), &
            ! elemR(iet(6),er_SlotArea), &
            ! 0.d0, &
            !    0.d0, &
            ! elemR(iet(7),er_SlotArea),&
            ! 0.d0, &
            !    0.d0, &
            ! elemR(iet(8),er_SlotArea)

            ! write(*,"(A,31f11.4)") '   ZCrown '  ,          &
            ! elemR(iet(1),er_Zcrown), &
            !    faceR(ift(1),fr_Zcrown_u), &
            !    faceR(ift(1),fr_Zcrown_d), &
            ! elemR(iet(2),er_Zcrown), &
            !    faceR(ift(2),fr_Zcrown_u), &
            !    faceR(ift(2),fr_Zcrown_d), &
            ! elemR(iet(3),er_Zcrown), &
            !    faceR(ift(3),fr_Zcrown_u), &
            !    faceR(ift(3),fr_Zcrown_d), &
            ! elemR(iet(4),er_Zcrown), &
            ! elemR(iet(5),er_Zcrown), &
            ! elemR(iet(6),er_Zcrown), &
            !    faceR(ift(4),fr_Zcrown_u), &
            !    faceR(ift(4),fr_Zcrown_d), &
            ! elemR(iet(7),er_Zcrown),&
            !    faceR(ift(5),fr_Zcrown_u), &
            !    faceR(ift(5),fr_Zcrown_d), &
            ! elemR(iet(8),er_Zcrown)

         ! write(*,"(A,31f11.4)") '   Zbot   '  ,          &
         !    elemR(iet(1),er_Zbottom), &
         !       faceR(ift(1),fr_Zbottom), &
         !       faceR(ift(1),fr_Zbottom), &
         !    elemR(iet(2),er_Zbottom), &
         !       faceR(ift(2),fr_Zbottom), &
         !       faceR(ift(2),fr_Zbottom), &
         !    elemR(iet(3),er_Zbottom), &
         !       faceR(ift(3),fr_Zbottom), &
         !       faceR(ift(3),fr_Zbottom), &
         !    elemR(iet(4),er_Zbottom), &
         !    elemR(iet(5),er_Zbottom), &
         !    elemR(iet(6),er_Zbottom), &
         !       faceR(ift(4),fr_Zbottom), &
         !       faceR(ift(4),fr_Zbottom), &
         !    elemR(iet(7),er_Zbottom),&
         !       faceR(ift(5),fr_Zbottom), &
         !       faceR(ift(5),fr_Zbottom), &
         !    elemR(iet(8),er_Zbottom)

            ! write(*,"(A,31e11.3)") ' Depth    '  ,          &
            ! elemR(iet(1),er_Depth), &
            !    faceR(ift(1),fr_Depth_u), &
            !    faceR(ift(1),fr_Depth_d), &
            ! elemR(iet(2),er_Depth), &
            !    faceR(ift(2),fr_Depth_u), &
            !    faceR(ift(2),fr_Depth_d), &
            ! elemR(iet(3),er_Depth), &
            !    faceR(ift(3),fr_Depth_u), &
            !    faceR(ift(3),fr_Depth_d), &
            ! elemR(iet(4),er_Depth), &
            ! elemR(iet(5),er_Depth), &
            ! elemR(iet(6),er_Depth), &
            !    faceR(ift(4),fr_Depth_u), &
            !    faceR(ift(4),fr_Depth_d), &
            ! elemR(iet(7),er_Depth),&
            !    faceR(ift(5),fr_Depth_u), &
            !    faceR(ift(5),fr_Depth_d), &
            ! elemR(iet(8),er_Depth)

            ! write(*,"(A,31f11.4)") '  EllD    '  ,          &
            ! elemR(iet(1),er_EllDepth), &
            ! 0.d0, &
            !    0.d0, &
            ! elemR(iet(2),er_EllDepth), &
            ! 0.d0, &
            !    0.d0, &
            ! elemR(iet(3),er_EllDepth), &
            !    0.d0, &
            !    0.d0, &
            ! elemR(iet(4),er_EllDepth), &
            ! elemR(iet(5),er_EllDepth), &
            ! elemR(iet(6),er_EllDepth), &
            ! 0.d0, &
            !    0.d0, &
            ! elemR(iet(7),er_EllDepth),&
            ! 0.d0, &
            !    0.d0, &
            ! elemR(iet(8),er_EllDepth)

         ! !    write(*,"(A,31e11.3)") ' Ell      '  ,          &
         ! !    elemR(iet(1),er_EllDepth), &
         ! !       faceR(ift(1),fr_Depth_u), &
         ! !       faceR(ift(1),fr_Depth_d), &
         ! !    elemR(iet(2),er_EllDepth), &
         ! !       faceR(ift(2),fr_Depth_u), &
         ! !       faceR(ift(2),fr_Depth_d), &
         ! !    elemR(iet(3),er_EllDepth), &
         ! !       faceR(ift(3),fr_Depth_u), &
         ! !       faceR(ift(3),fr_Depth_d), &
         ! !    elemR(iet(4),er_EllDepth), &
         ! !    elemR(iet(5),er_EllDepth), &
         ! !    elemR(iet(6),er_EllDepth), &
         ! !       faceR(ift(4),fr_Depth_u), &
         ! !       faceR(ift(4),fr_Depth_d), &
         ! !    elemR(iet(7),er_EllDepth),&
         ! !       faceR(ift(5),fr_Depth_u), &
         ! !       faceR(ift(5),fr_Depth_d), &
         ! !    elemR(iet(8),er_EllDepth)

            ! print *, ' '

            ! write(*,"(A,31e11.3)") '   Flow   '  ,          &
            ! elemR(iet(1),er_Flowrate), &
            !    faceR(ift(1),fr_Flowrate), &
            !    faceR(ift(1),fr_Flowrate), &
            ! elemR(iet(2),er_Flowrate), &
            !    faceR(ift(2),fr_Flowrate), &
            !    faceR(ift(2),fr_Flowrate), &
            ! elemR(iet(3),er_Flowrate), &
            !    faceR(ift(3),fr_Flowrate), &
            !    faceR(ift(3),fr_Flowrate), &
            ! elemR(iet(4),er_Flowrate), &
            ! elemR(iet(5),er_Flowrate), &
            ! elemR(iet(6),er_Flowrate), &
            !    faceR(ift(4),fr_Flowrate), &
            !    faceR(ift(4),fr_Flowrate), &
            ! elemR(iet(7),er_Flowrate),&
            !    faceR(ift(5),fr_Flowrate), &
            !    faceR(ift(5),fr_Flowrate), &
            ! elemR(iet(8),er_Flowrate)

            ! print *, ' '

         !    write(*,"(A,31e11.3)") '   dQ     '  ,          &
         !    elemR(iet(1),er_DeltaQ), &
         !       faceR(ift(1),fr_DeltaQ), &
         !       faceR(ift(1),fr_DeltaQ), &
         !    elemR(iet(2),er_DeltaQ), &
         !       faceR(ift(2),fr_DeltaQ), &
         !       faceR(ift(2),fr_DeltaQ), &
         !    elemR(iet(3),er_DeltaQ), &
         !       faceR(ift(3),fr_DeltaQ), &
         !       faceR(ift(3),fr_DeltaQ), &
         !    elemR(iet(4),er_DeltaQ), &
         !    elemR(iet(5),er_DeltaQ), &
         !    elemR(iet(6),er_DeltaQ), &
         !       faceR(ift(4),fr_DeltaQ), &
         !       faceR(ift(4),fr_DeltaQ), &
         !    elemR(iet(7),er_DeltaQ),&
         !       faceR(ift(5),fr_DeltaQ), &
         !       faceR(ift(5),fr_DeltaQ), &
         !    elemR(iet(8),er_DeltaQ)

         !    print *, ' '

            ! write(*,"(A,31e11.3)") '   Qcons  '  ,          &
            ! elemR(iet(1),er_Flowrate), &
            !    faceR(ift(1),fr_Flowrate_Conservative), &
            !    faceR(ift(1),fr_Flowrate_Conservative), &
            ! elemR(iet(2),er_Flowrate), &
            !    faceR(ift(2),fr_Flowrate_Conservative), &
            !    faceR(ift(2),fr_Flowrate_Conservative), &
            ! elemR(iet(3),er_Flowrate), &
            !    faceR(ift(3),fr_Flowrate_Conservative), &
            !    faceR(ift(3),fr_Flowrate_Conservative), &
            ! elemR(iet(4),er_Flowrate), &
            ! elemR(iet(5),er_Flowrate), &
            ! elemR(iet(6),er_Flowrate), &
            !    faceR(ift(4),fr_Flowrate_Conservative), &
            !    faceR(ift(4),fr_Flowrate_Conservative), &
            ! elemR(iet(7),er_Flowrate),&
            !    faceR(ift(5),fr_Flowrate_Conservative), &
            !    faceR(ift(5),fr_Flowrate_Conservative), &
            ! elemR(iet(8),er_Flowrate)

            ! print *, ' '

         ! write(*,"(A,31e11.3)") '   Volu   '  ,          &
         !    elemR(iet(1),er_Volume), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(2),er_Volume), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(3),er_Volume), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(4),er_Volume), &
         !    elemR(iet(5),er_Volume), &
         !    elemR(iet(6),er_Volume), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(7),er_Volume),&
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(8),er_Volume)

            ! write(*,"(A,31e11.3)") '  Vol/area'  ,          &
            ! elemR(iet(1),er_Volume) / (elemR(iet(1),er_Length) * elemR(iet(1),er_Topwidth)), &
            !    0.d0, &
            !    0.d0, &
            ! elemR(iet(2),er_Volume)  / (elemR(iet(2),er_Length) * elemR(iet(2),er_Topwidth)), &
            !    0.d0, &
            !    0.d0, &
            ! elemR(iet(3),er_Volume) / (elemR(iet(3),er_Length) * elemR(iet(3),er_Topwidth)), &
            !    0.d0, &
            !    0.d0, &
            ! elemR(iet(4),er_Volume) / (elemR(iet(3),er_Length) * elemR(iet(3),er_Topwidth)), &
            ! elemR(iet(5),er_Volume) / elemSR(iet(5),esr_Storage_Plan_Area), &
            ! elemR(iet(6),er_Volume) / (elemR(iet(3),er_Length) * elemR(iet(3),er_Topwidth)), &
            !    0.d0, &
            !    0.d0, &
            ! elemR(iet(7),er_Volume) / (elemR(iet(7),er_Length) * elemR(iet(7),er_Topwidth)),&
            !    0.d0, &
            !    0.d0, &
            ! elemR(iet(8),er_Volume) / (elemR(iet(8),er_Length) * elemR(iet(8),er_Topwidth))

            ! print *, ' '
          



         ! write(*,"(A,31f11.4)") '   Vel    '  ,          &
         !    elemR(iet(1),er_Velocity), &
         !       faceR(ift(1),fr_Velocity_u), &
         !       faceR(ift(1),fr_Velocity_d), &
         !    elemR(iet(2),er_Velocity), &
         !       faceR(ift(2),fr_Velocity_u), &
         !       faceR(ift(2),fr_Velocity_d), &
         !    elemR(iet(3),er_Velocity), &
         !       faceR(ift(3),fr_Velocity_u), &
         !       faceR(ift(3),fr_Velocity_d), &   
         !    elemR(iet(4),er_Velocity), &
         !    elemR(iet(5),er_Velocity), &
         !    elemR(iet(6),er_Velocity), &
         !       faceR(ift(4),fr_Velocity_u), &
         !       faceR(ift(4),fr_Velocity_d), &
         !    elemR(iet(7),er_Velocity),&
         !       faceR(ift(5),fr_Velocity_u), &
         !       faceR(ift(5),fr_Velocity_d), &
         !    elemR(iet(8),er_Velocity)

         !    write(*,"(A,31e11.3)") '   wave   '  ,          &
         !    elemR(iet(1),er_WaveSpeed), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(2),er_WaveSpeed), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(3),er_WaveSpeed), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(4),er_WaveSpeed), &
         !    elemR(iet(5),er_WaveSpeed), &
         !    elemR(iet(6),er_WaveSpeed), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(7),er_WaveSpeed),&
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(8),er_WaveSpeed)

         ! write(*,"(A,31e11.3)") '   IWQ    '  ,          &
         !    0.d0  , &
         !       elemR(iet(1),er_InterpWeight_dQ), &
         !       elemR(iet(2),er_InterpWeight_uQ), &
         !    0.d0, &
         !       elemR(iet(2),er_InterpWeight_dQ), &
         !       elemR(iet(3),er_InterpWeight_uQ), &
         !    0.d0, &
         !       elemR(iet(3),er_InterpWeight_dQ), &
         !       elemR(iet(4),er_InterpWeight_uQ), & 
         !    0.d0, &
         !    0.d0, &
         !    0.d0, &
         !       elemR(iet(6),er_InterpWeight_dQ), &
         !       elemR(iet(7),er_InterpWeight_uQ), &
         !    0.d0, &
         !       elemR(iet(7),er_InterpWeight_dQ), &
         !       elemR(iet(8),er_InterpWeight_uQ), &
         !    0.d0

         !    write(*,"(A,31e11.3)") '   iwG    '  ,          &
         !    0.d0  , &
         !       elemR(iet(1),er_InterpWeight_dG), &
         !       elemR(iet(2),er_InterpWeight_uG), &
         !    0.d0, &
         !       elemR(iet(2),er_InterpWeight_dG), &
         !       elemR(iet(3),er_InterpWeight_uG), &
         !    0.d0, &
         !       elemR(iet(3),er_InterpWeight_dG), &
         !       elemR(iet(4),er_InterpWeight_uG), & 
         !    0.d0, &
         !    0.d0, &
         !    0.d0, &
         !       elemR(iet(6),er_InterpWeight_dG), &
         !       elemR(iet(7),er_InterpWeight_uG), &
         !    0.d0, &
         !       elemR(iet(7),er_InterpWeight_dG), &
         !       elemR(iet(8),er_InterpWeight_uG), &
         !    0.d0

         !    write(*,"(A,31e11.3)") '   iwH    '  ,          &
         !    0.d0  , &
         !       elemR(iet(1),er_InterpWeight_dH), &
         !       elemR(iet(2),er_InterpWeight_uH), &
         !    0.d0, &
         !       elemR(iet(2),er_InterpWeight_dH), &
         !       elemR(iet(3),er_InterpWeight_uH), &
         !    0.d0, &
         !       elemR(iet(3),er_InterpWeight_dH), &
         !       elemR(iet(4),er_InterpWeight_uH), & 
         !    0.d0, &
         !    0.d0, &
         !    0.d0, &
         !       elemR(iet(6),er_InterpWeight_dH), &
         !       elemR(iet(7),er_InterpWeight_uH), &
         !    0.d0, &
         !       elemR(iet(7),er_InterpWeight_dH), &
         !       elemR(iet(8),er_InterpWeight_uH), &
         !    0.d0

            ! write(*,"(A,31e11.3)") '   Fr     '  ,          &
            ! elemR(iet(1),er_FroudeNumber), &
            !    0.d0, &
            !    0.d0, &
            ! elemR(iet(2),er_FroudeNumber), &
            !    0.d0, &
            !    0.d0, &
            ! elemR(iet(3),er_FroudeNumber), &
            !    0.d0, &
            !    0.d0, &
            ! elemR(iet(4),er_FroudeNumber), &
            ! elemR(iet(5),er_FroudeNumber), &
            ! elemR(iet(6),er_FroudeNumber), &
            !    0.d0, &
            !    0.d0, &
            ! elemR(iet(7),er_FroudeNumber),&
            !    0.d0, &
            !    0.d0, &
            ! elemR(iet(8),er_FroudeNumber)


            ! return

!%  
!% 8/5  1 elem  J 4 elem
!%
   !  write(*,"(A9,20A11)"),' ','elem','faceU','faceD','jb','JM','jb','faceU','faceD','elem','faceU','faceD','elem','faceU','faceD','elem'
           ! 1 junction     8 elm, 5 face looking at downstream
     
         ! write(*,"(A,31f11.4)") '  Head    '  ,          &
         !    elemR(iet(1),er_Head), &
         !       faceR(ift(1),fr_Head_u), &
         !       faceR(ift(1),fr_Head_d), &
         !    elemR(iet(2),er_Head), &
         !    elemR(iet(3),er_Head), &
         !    elemR(iet(4),er_Head), &
         !       faceR(ift(2),fr_Head_u), &
         !       faceR(ift(2),fr_Head_d), &
         !    elemR(iet(5),er_Head), &
         !       faceR(ift(3),fr_Head_u), &
         !       faceR(ift(3),fr_Head_d), &
         !    elemR(iet(6),er_Head), &
         !       faceR(ift(4),fr_Head_u), &
         !       faceR(ift(4),fr_Head_d), &
         !    elemR(iet(7),er_Head),&
         !       faceR(ift(5),fr_Head_u), &
         !       faceR(ift(5),fr_Head_d), &
         !    elemR(iet(8),er_Head)

         ! write(*,"(A,31f11.4)") '  Depth   '  ,          &
         !    elemR(iet(1),er_Depth), &
         !       faceR(ift(1),fr_Depth_u), &
         !       faceR(ift(1),fr_Depth_d), &
         !    elemR(iet(2),er_Depth), &
         !    elemR(iet(3),er_Depth), &
         !    elemR(iet(4),er_Depth), &
         !       faceR(ift(2),fr_Depth_u), &
         !       faceR(ift(2),fr_Depth_d), &
         !    elemR(iet(5),er_Depth), &
         !       faceR(ift(3),fr_Depth_u), &
         !       faceR(ift(3),fr_Depth_d), &
         !    elemR(iet(6),er_Depth), &
         !       faceR(ift(4),fr_Depth_u), &
         !       faceR(ift(4),fr_Depth_d), &
         !    elemR(iet(7),er_Depth),&
         !       faceR(ift(5),fr_Depth_u), &
         !       faceR(ift(5),fr_Depth_d), &
         !    elemR(iet(8),er_Depth)

            ! write(*,"(A,31f11.4)") '  Area    '  ,          &
            ! elemR(iet(1),er_Area), &
            !    faceR(ift(1),fr_Area_u), &
            !    faceR(ift(1),fr_Area_d), &
            ! elemR(iet(2),er_Area), &
            ! elemR(iet(3),er_Area), &
            ! elemR(iet(4),er_Area), &
            !    faceR(ift(2),fr_Area_u), &
            !    faceR(ift(2),fr_Area_d), &
            ! elemR(iet(5),er_Area), &
            !    faceR(ift(3),fr_Area_u), &
            !    faceR(ift(3),fr_Area_d), &
            ! elemR(iet(6),er_Area), &
            !    faceR(ift(4),fr_Area_u), &
            !    faceR(ift(4),fr_Area_d), &
            ! elemR(iet(7),er_Area),&
            !    faceR(ift(5),fr_Area_u), &
            !    faceR(ift(5),fr_Area_d), &
            ! elemR(iet(8),er_Area)

            ! print *, ' '

            ! write(*,"(A,31e11.3)") '  Flow    '  ,          &
            ! elemR(iet(1),er_Flowrate), &
            !    faceR(ift(1),fr_Flowrate), &
            !    faceR(ift(1),fr_Flowrate), &
            ! elemR(iet(2),er_Flowrate), &
            ! elemR(iet(3),er_Flowrate), &
            ! elemR(iet(4),er_Flowrate), &
            !    faceR(ift(2),fr_Flowrate), &
            !    faceR(ift(2),fr_Flowrate), &
            ! elemR(iet(5),er_Flowrate), &
            !    faceR(ift(3),fr_Flowrate), &
            !    faceR(ift(3),fr_Flowrate), &
            ! elemR(iet(6),er_Flowrate), &
            !    faceR(ift(4),fr_Flowrate), &
            !    faceR(ift(4),fr_Flowrate), &
            ! elemR(iet(7),er_Flowrate),&
            !    faceR(ift(5),fr_Flowrate), &
            !    faceR(ift(5),fr_Flowrate), &
            ! elemR(iet(8),er_Flowrate)

         !    write(*,"(A,31e11.3)") '  Qdif    '  ,          &
         !    elemR(iet(1),er_Flowrate) - elemR(iet(5),er_Flowrate), &
         !       faceR(ift(1),fr_Flowrate) - elemR(iet(5),er_Flowrate), &
         !       faceR(ift(1),fr_Flowrate) - elemR(iet(5),er_Flowrate), &
         !    elemR(iet(2),er_Flowrate) - elemR(iet(5),er_Flowrate), &
         !    elemR(iet(3),er_Flowrate) - elemR(iet(5),er_Flowrate), &
         !    elemR(iet(4),er_Flowrate) - elemR(iet(5),er_Flowrate), &
         !       faceR(ift(2),fr_Flowrate) - elemR(iet(5),er_Flowrate), &
         !       faceR(ift(2),fr_Flowrate) - elemR(iet(5),er_Flowrate), &
         !    elemR(iet(5),er_Flowrate) - elemR(iet(5),er_Flowrate), &
         !       faceR(ift(3),fr_Flowrate)- elemR(iet(5),er_Flowrate), &
         !       faceR(ift(3),fr_Flowrate) - elemR(iet(5),er_Flowrate), &
         !    elemR(iet(6),er_Flowrate) - elemR(iet(5),er_Flowrate), &
         !       faceR(ift(4),fr_Flowrate) - elemR(iet(5),er_Flowrate), &
         !       faceR(ift(4),fr_Flowrate) - elemR(iet(5),er_Flowrate), &
         !    elemR(iet(7),er_Flowrate) - elemR(iet(5),er_Flowrate),&
         !       faceR(ift(5),fr_Flowrate) -elemR(iet(5),er_Flowrate), &
         !       faceR(ift(5),fr_Flowrate) - elemR(iet(5),er_Flowrate), &
         !    elemR(iet(8),er_Flowrate) -elemR(iet(5),er_Flowrate)


         !    print *, ' '
         !    write(*,"(A,31e11.4)") '  gammaC  '  ,          &
         !    elemR(iet(1),er_gammaC ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(2),er_gammaC ), &
         !    elemR(iet(3),er_gammaC ), &
         !    elemR(iet(4),er_gammaC ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(5),er_gammaC ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(6),er_gammaC ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(7),er_gammaC ),&
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(8),er_gammaC )

         !    write(*,"(A,31e11.4)") '  gammaM  '  ,          &
         !    elemR(iet(1),er_gammaM ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(2),er_gammaM ), &
         !    elemR(iet(3),er_gammaM ), &
         !    elemR(iet(4),er_gammaM ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(5),er_gammaM ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(6),er_gammaM ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(7),er_gammaM ),&
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(8),er_gammaM )

         !    write(*,"(A,31e11.4)") '  SrcC    '  ,          &
         !    elemR(iet(1),er_SourceContinuity ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(2),er_SourceContinuity ), &
         !    elemR(iet(3),er_SourceContinuity ), &
         !    elemR(iet(4),er_SourceContinuity ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(5),er_SourceContinuity ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(6),er_SourceContinuity ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(7),er_SourceContinuity ),&
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(8),er_SourceContinuity )
 
         !    write(*,"(A,31e11.4)") '  SrcM    '  ,          &
         !    elemR(iet(1),er_SourceMomentum ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(2),er_SourceMomentum ), &
         !    elemR(iet(3),er_SourceMomentum ), &
         !    elemR(iet(4),er_SourceMomentum ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(5),er_SourceMomentum ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(6),er_SourceMomentum ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(7),er_SourceMomentum ),&
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(8),er_SourceMomentum )

         !    write(*,"(A,31e11.4)") '  Ksrc    '  ,          &
         !    elemR(iet(1),er_KSource ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(2),er_KSource ), &
         !    elemR(iet(3),er_KSource ), &
         !    elemR(iet(4),er_KSource ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(5),er_KSource ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(6),er_KSource ), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(7),er_KSource ),&
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(8),er_KSource )

         !    print *, ' '

            ! write(*,"(A,31f11.4)") '  Fr      '  ,          &
            ! elemR(iet(1),er_FroudeNumber), &
            !    faceR(ift(1),fr_FroudeNumber_u), &
            !    faceR(ift(1),fr_FroudeNumber_d), &
            ! elemR(iet(2),er_FroudeNumber), &
            ! elemR(iet(3),er_FroudeNumber), &
            ! elemR(iet(4),er_FroudeNumber), &
            !    faceR(ift(2),fr_FroudeNumber_u), &
            !    faceR(ift(2),fr_FroudeNumber_d), &
            ! elemR(iet(5),er_FroudeNumber), &
            !    faceR(ift(3),fr_FroudeNumber_u), &
            !    faceR(ift(3),fr_FroudeNumber_d), &
            ! elemR(iet(6),er_FroudeNumber), &
            !    faceR(ift(4),fr_FroudeNumber_u), &
            !    faceR(ift(4),fr_FroudeNumber_d), &
            ! elemR(iet(7),er_FroudeNumber),&
            !    faceR(ift(5),fr_FroudeNumber_u), &
            !    faceR(ift(5),fr_FroudeNumber_d), &
            ! elemR(iet(8),er_FroudeNumber)

            ! write(*,"(A,31f11.4)") '  Velo    '  ,          &
            ! elemR(iet(1),er_Velocity), &
            !    faceR(ift(1),fr_Velocity_u), &
            !    faceR(ift(1),fr_Velocity_d), &
            ! elemR(iet(2),er_Velocity), &
            ! elemR(iet(3),er_Velocity), &
            ! elemR(iet(4),er_Velocity), &
            !    faceR(ift(2),fr_Velocity_u), &
            !    faceR(ift(2),fr_Velocity_d), &
            ! elemR(iet(5),er_Velocity), &
            !    faceR(ift(3),fr_Velocity_u), &
            !    faceR(ift(3),fr_Velocity_d), &
            ! elemR(iet(6),er_Velocity), &
            !    faceR(ift(4),fr_Velocity_u), &
            !    faceR(ift(4),fr_Velocity_d), &
            ! elemR(iet(7),er_Velocity),&
            !    faceR(ift(5),fr_Velocity_u), &
            !    faceR(ift(5),fr_Velocity_d), &
            ! elemR(iet(8),er_Velocity)

         !    print *, ' '
         !    write(*,"(A,31e11.3)") '   iwG    '  ,          &
         !    0.d0  , &
         !       elemR(iet(1),er_InterpWeight_dG), &
         !       elemR(iet(2),er_InterpWeight_uG), &
         !    0.d0, &
         !    0.d0, &
         !    0.d0, &
         !       elemR(iet(4),er_InterpWeight_dG), &
         !       elemR(iet(5),er_InterpWeight_uG), &
         !    0.d0, &
         !       elemR(iet(5),er_InterpWeight_dG), &
         !       elemR(iet(6),er_InterpWeight_uG), & 
         !    0.d0, &
         !       elemR(iet(6),er_InterpWeight_dG), &
         !       elemR(iet(7),er_InterpWeight_uG), &
         !    0.d0, &
         !       elemR(iet(7),er_InterpWeight_dG), &
         !       elemR(iet(8),er_InterpWeight_uG), &
         !    0.d0



            ! return

!%
!%  11/6   Weir JM with 2 junctions
!%
   ! write(*,"(A9,30A11)"),' ','elem','faceU','faceD','elem','faceU','faceD','jb','JM','jb','faceU','faceD','WEIR','faceU','faceD','jb','JM','jb','faceU','faceD','elem','faceU','faceD','elem'
      !% 11 elem, 6 face WEIR JM -- 2 junctions
         ! write(*,"(A,31f11.4)") '  Head    '  ,          &
         !    elemR(iet(1),er_Head), &
         !       faceR(ift(1),fr_Head_u), &
         !       faceR(ift(1),fr_Head_d), &
         !    elemR(iet(2),er_Head), &
         !       faceR(ift(2),fr_Head_u), &
         !       faceR(ift(2),fr_Head_d), &
         !    elemR(iet(3),er_Head), &
         !    elemR(iet(4),er_Head), &
         !    elemR(iet(5),er_Head), &
         !       faceR(ift(3),fr_Head_u), &
         !       faceR(ift(3),fr_Head_d), &
         !    elemR(iet(6),er_Head),&
         !       faceR(ift(4),fr_Head_u), &
         !       faceR(ift(4),fr_Head_d), &
         !    elemR(iet(7),er_Head), &
         !    elemR(iet(8),er_Head), &
         !    elemR(iet(9),er_Head), &
         !       faceR(ift(5),fr_Head_u), &
         !       faceR(ift(5),fr_Head_d), &
         !    elemR(iet(10),er_Head), &
         !       faceR(ift(6),fr_Head_u), &
         !       faceR(ift(6),fr_Head_d), &
         !    elemR(iet(11),er_Head)

         ! write(*,"(A,31f11.4)") '  ZBtm    '  ,          &
         !    elemR(iet(1),er_Zbottom), &
         !       faceR(ift(1),fr_Zbottom), &
         !       faceR(ift(1),fr_Zbottom), &
         !    elemR(iet(2),er_Zbottom), &
         !       faceR(ift(2),fr_Zbottom), &
         !       faceR(ift(2),fr_Zbottom), &
         !    elemR(iet(3),er_Zbottom), &
         !    elemR(iet(4),er_Zbottom), &
         !    elemR(iet(5),er_Zbottom), &
         !       faceR(ift(3),fr_Zbottom), &
         !       faceR(ift(3),fr_Zbottom), &
         !    elemR(iet(6),er_Zbottom),&
         !       faceR(ift(4),fr_Zbottom), &
         !       faceR(ift(4),fr_Zbottom), &
         !    elemR(iet(7),er_Zbottom), &
         !    elemR(iet(8),er_Zbottom), &
         !    elemR(iet(9),er_Zbottom), &
         !       faceR(ift(5),fr_Zbottom), &
         !       faceR(ift(5),fr_Zbottom), &
         !    elemR(iet(10),er_Zbottom), &
         !       faceR(ift(6),fr_Zbottom), &
         !       faceR(ift(6),fr_Zbottom), &
         !    elemR(iet(11),er_Zbottom)


         ! write(*,"(A,31f11.4)") '  Depth   '  ,          &
         !    elemR(iet(1),er_Depth), &
         !       faceR(ift(1),fr_Depth_u), &
         !       faceR(ift(1),fr_Depth_d), &
         !    elemR(iet(2),er_Depth), &
         !       faceR(ift(2),fr_Depth_u), &
         !       faceR(ift(2),fr_Depth_d), &
         !    elemR(iet(3),er_Depth), &
         !    elemR(iet(4),er_Depth), &
         !    elemR(iet(5),er_Depth), &
         !       faceR(ift(3),fr_Depth_u), &
         !       faceR(ift(3),fr_Depth_d), &
         !    elemR(iet(6),er_Depth),&
         !       faceR(ift(4),fr_Depth_u), &
         !       faceR(ift(4),fr_Depth_d), &
         !    elemR(iet(7),er_Depth), &
         !    elemR(iet(8),er_Depth), &
         !    elemR(iet(9),er_Depth), &
         !       faceR(ift(5),fr_Depth_u), &
         !       faceR(ift(5),fr_Depth_d), &
         !    elemR(iet(10),er_Depth), &
         !       faceR(ift(6),fr_Depth_u), &
         !       faceR(ift(6),fr_Depth_d), &
         !    elemR(iet(11),er_Depth)

         !    print *, ' '

         ! write(*,"(A,31e11.3)") '  Flow    '  ,          &
         !    elemR(iet(1),er_Flowrate), &
         !       faceR(ift(1),fr_Flowrate), &
         !       faceR(ift(1),fr_Flowrate), &
         !    elemR(iet(2),er_Flowrate), &
         !       faceR(ift(2),fr_Flowrate), &
         !       faceR(ift(2),fr_Flowrate), &
         !    elemR(iet(3),er_Flowrate), &
         !    elemR(iet(4),er_Flowrate), &
         !    elemR(iet(5),er_Flowrate), &
         !       faceR(ift(3),fr_Flowrate), &
         !       faceR(ift(3),fr_Flowrate), &
         !    elemR(iet(6),er_Flowrate),&
         !       faceR(ift(4),fr_Flowrate), &
         !       faceR(ift(4),fr_Flowrate), &
         !    elemR(iet(7),er_Flowrate), &
         !    elemR(iet(8),er_Flowrate), &
         !    elemR(iet(9),er_Flowrate), &
         !       faceR(ift(5),fr_Flowrate), &
         !       faceR(ift(5),fr_Flowrate), &
         !    elemR(iet(10),er_Flowrate), &
         !       faceR(ift(6),fr_Flowrate), &
         !       faceR(ift(6),fr_Flowrate), &
         !    elemR(iet(11),er_Flowrate)


         ! write(*,"(A,31e11.3)") '  Qcons   '  ,          &
         !    elemR(iet(1),er_Flowrate), &
         !       faceR(ift(1),fr_Flowrate_Conservative), &
         !       faceR(ift(1),fr_Flowrate_Conservative), &
         !    elemR(iet(2),er_Flowrate), &
         !       faceR(ift(2),fr_Flowrate_Conservative), &
         !       faceR(ift(2),fr_Flowrate_Conservative), &
         !    elemR(iet(3),er_Flowrate), &
         !    elemR(iet(4),er_Flowrate), &
         !    elemR(iet(5),er_Flowrate), &
         !       faceR(ift(3),fr_Flowrate_Conservative), &
         !       faceR(ift(3),fr_Flowrate_Conservative), &
         !    elemR(iet(6),er_Flowrate),&
         !       faceR(ift(4),fr_Flowrate_Conservative), &
         !       faceR(ift(4),fr_Flowrate_Conservative), &
         !    elemR(iet(7),er_Flowrate), &
         !    elemR(iet(8),er_Flowrate), &
         !    elemR(iet(9),er_Flowrate), &
         !       faceR(ift(5),fr_Flowrate_Conservative), &
         !       faceR(ift(5),fr_Flowrate_Conservative), &
         !    elemR(iet(10),er_Flowrate), &
         !       faceR(ift(6),fr_Flowrate_Conservative), &
         !       faceR(ift(6),fr_Flowrate_Conservative), &
         !    elemR(iet(11),er_Flowrate)


         ! write(*,"(A,31e11.4)") '   Vol     '  ,          &
         !       elemR(iet(1),er_Volume), &
         !          0.d0, &
         !          0.d0, &
         !       elemR(iet(2),er_Volume), &
         !       0.d0, &
         !          0.d0, &
         !       elemR(iet(3),er_Volume), &
         !       elemR(iet(4),er_Volume), &
         !       elemR(iet(5),er_Volume), &
         !          0.d0, &
         !          0.d0, &
         !       elemR(iet(6),er_Volume),&
         !       0.d0, &
         !          0.d0, &
         !       elemR(iet(7),er_Volume), &
         !       elemR(iet(8),er_Volume), &
         !       elemR(iet(9),er_Volume), &
         !       0.d0, &
         !          0.d0, &
         !       elemR(iet(10),er_Volume), &
         !          0.d0, &
         !          0.d0, &
         !       elemR(iet(11),er_Flowrate)


                  ! return
!%    
!%   is zero formats        
!%            
   ! write(*,"(A,A,L3,  A,A,A,L3,   A,A,A,L3,  A,A,A,L3,    A,A,L3,  A,A,L3,  A,A,A,L3,   A,A,A,L3)"),&
                  ! 'Zero            ',&
                  ! '          ', &
                  ! elemYN(iet(1),eYN_isZeroDepth),&
                  ! '          ', &
                  ! '          ', &
                  ! '          ', &
                  ! elemYN(iet(2),eYN_isZeroDepth),&
                  ! '          ', &
                  ! '          ', &
                  ! '          ', &
                  ! elemYN(iet(3),eYN_isZeroDepth),&
                  ! '          ', &
                  ! '          ', &
                  ! '          ', &
                  ! elemYN(iet(4),eYN_isZeroDepth),&
                  ! '    ', &
                  ! '    ', &
                  ! elemYN(iet(5),eYN_isZeroDepth),&
                  ! '    ', &
                  ! '    ', &
                  ! elemYN(iet(6),eYN_isZeroDepth),&
                  ! '          ', &
                  ! '          ', &
                  ! '          ', &
                  ! elemYN(iet(7),eYN_isZeroDepth),&
                  ! '          ', &
                  ! '          ', &
                  ! '          ', &
                  ! elemYN(iet(8),eYN_isZeroDepth )
     

            ! write(*,"(A, A,L3, A,A,A, L3, A,A,A, L3, A,A, L3, A,A, L3, A,A,A, L3, A,A, A,  L3, A,A,  L3, A,A, L3, A, A, L3)"),&
               ! 'Small           ',&
               ! '          ', &
               ! elemYN(iet(1),eYN_isSmallDepth),&
               ! '          ', &
               ! '          ', &
               ! '          ', &
               ! elemYN(iet(2),eYN_isSmallDepth),&
               ! '          ', &
               ! '          ', &
               ! '          ', &
               ! elemYN(iet(3),eYN_isSmallDepth),&
               ! '    ', &
               ! '    ', &
               ! elemYN(iet(4),eYN_isSmallDepth),&
               ! '    ', &
               ! '    ', &
               ! elemYN(iet(5),eYN_isSmallDepth),&
               ! '          ', &
               ! '          ', &
               ! '          ', &
               ! elemYN(iet(6),eYN_isSmallDepth),&
               ! '          ', &
               ! '          ', &
               ! '          ', &
               ! elemYN(iet(7),eYN_isSmallDepth )

   
!%
!% 8/7 starts with face, 3 elem J 2 elem
!% 
           ! print *, 'APlan 46',elemSR(46,esr_Storage_Plan_Area)
    write(*,"(A9,20A11)"),' ','faceD','elem','faceU','faceD','elem','faceU','faceD','elem','faceU','faceD','jb','JM','jb','faceU','faceD','elem','faceU','faceD','elem', 'faceU','faceD'
            !% ` junction
         write(*,"(A,31f11.4)") '  Head    '  ,          &
                faceR(ift(1),fr_Head_d), &
            elemR(iet(1),er_Head), &
                faceR(ift(2),fr_Head_u), &
                faceR(ift(2),fr_Head_d), &
            elemR(iet(2),er_Head), &
                faceR(ift(3),fr_Head_u), &
                faceR(ift(3),fr_Head_d), &
            elemR(iet(3),er_Head), &
                faceR(ift(4),fr_Head_u), &
                faceR(ift(4),fr_Head_d), &
            elemR(iet(4),er_Head),&
            elemR(iet(5),er_Head),&
            elemR(iet(6),er_Head), &
                faceR(ift(5),fr_Head_u), &
                faceR(ift(5),fr_Head_d), &
            elemR(iet(7),er_Head), &
                faceR(ift(6),fr_Head_u), &
                faceR(ift(6),fr_Head_d), &
            elemR(iet(8),er_Head), &
                faceR(ift(7),fr_Head_u), &
                faceR(ift(7),fr_Head_d)

         write(*,"(A,31f11.4)") '  Zbtm    '  ,          &
            faceR(ift(1),fr_Zbottom), &
            elemR(iet(1),er_Zbottom), &
            faceR(ift(2),fr_Zbottom), &
            faceR(ift(2),fr_Zbottom), &
            elemR(iet(2),er_Zbottom), &
            faceR(ift(3),fr_Zbottom), &
            faceR(ift(3),fr_Zbottom), &
            elemR(iet(3),er_Zbottom), &
            faceR(ift(4),fr_Zbottom), &
            faceR(ift(4),fr_Zbottom), &
            elemR(iet(4),er_Zbottom),&
            elemR(iet(5),er_Zbottom),&
            elemR(iet(6),er_Zbottom), &
            faceR(ift(5),fr_Zbottom), &
            faceR(ift(5),fr_Zbottom), &
            elemR(iet(7),er_Zbottom), &
            faceR(ift(6),fr_Zbottom), &
            faceR(ift(6),fr_Zbottom), &
            elemR(iet(8),er_Zbottom), &
            faceR(ift(7),fr_Zbottom), &
            faceR(ift(7),fr_Zbottom)

            print *, ' '

         write(*,"(A,31f11.4)") '  Depth   '  ,          &
               faceR(ift(1),fr_Depth_d), &
            elemR(iet(1),er_Depth), &
               faceR(ift(2),fr_Depth_u), &
               faceR(ift(2),fr_Depth_d), &
            elemR(iet(2),er_Depth), &
               faceR(ift(3),fr_Depth_u), &
               faceR(ift(3),fr_Depth_d), &
            elemR(iet(3),er_Depth), &
               faceR(ift(4),fr_Depth_u), &
               faceR(ift(4),fr_Depth_d), &
            elemR(iet(4),er_Depth),&
            elemR(iet(5),er_Depth),&
            elemR(iet(6),er_Depth), &
               faceR(ift(5),fr_Depth_u), &
               faceR(ift(5),fr_Depth_d), &
            elemR(iet(7),er_Depth), &
               faceR(ift(6),fr_Depth_u), &
               faceR(ift(6),fr_Depth_d), &
            elemR(iet(8),er_Depth), &
               faceR(ift(7),fr_Depth_u), &
               faceR(ift(7),fr_Depth_d)

            write(*,"(A,31f11.4)") 'ELLDepth  '  ,          &
               faceR(ift(1),fr_Depth_d), &
            elemR(iet(1),er_EllDepth), &
               faceR(ift(2),fr_Depth_u), &
               faceR(ift(2),fr_Depth_d), &
            elemR(iet(2),er_EllDepth), &
               faceR(ift(3),fr_Depth_u), &
               faceR(ift(3),fr_Depth_d), &
            elemR(iet(3),er_EllDepth), &
               faceR(ift(4),fr_Depth_u), &
               faceR(ift(4),fr_Depth_d), &
            elemR(iet(4),er_EllDepth),&
            elemR(iet(5),er_EllDepth),&
            elemR(iet(6),er_EllDepth), &
               faceR(ift(5),fr_Depth_u), &
               faceR(ift(5),fr_Depth_d), &
            elemR(iet(7),er_EllDepth), &
               faceR(ift(6),fr_Depth_u), &
               faceR(ift(6),fr_Depth_d), &
            elemR(iet(8),er_EllDepth), &
               faceR(ift(7),fr_Depth_u), &
               faceR(ift(7),fr_Depth_d)


               print *, ' '

         write(*,"(A,31e11.4)") '  Area    '  ,          &
                     faceR(ift(1),fr_Area_d), &
                  elemR(iet(1),er_Area), &
                     faceR(ift(2),fr_Area_u), &
                     faceR(ift(2),fr_Area_d), &
                  elemR(iet(2),er_Area), &
                     faceR(ift(3),fr_Area_u), &
                     faceR(ift(3),fr_Area_d), &
                  elemR(iet(3),er_Area), &
                     faceR(ift(4),fr_Area_u), &
                     faceR(ift(4),fr_Area_d), &
                  elemR(iet(4),er_Area),&
                  elemR(iet(5),er_Area),&
                  elemR(iet(6),er_Area), &
                     faceR(ift(5),fr_Area_u), &
                     faceR(ift(5),fr_Area_d), &
                  elemR(iet(7),er_Area), &
                     faceR(ift(6),fr_Area_u), &
                     faceR(ift(6),fr_Area_d), &
                  elemR(iet(8),er_Area), &
                  faceR(ift(7),fr_Area_u), &
                  faceR(ift(7),fr_Area_d)

            write(*,"(A,31e11.4)") '  topW    '  ,          &
                  0.d0, &
               elemR(iet(1),er_Topwidth), &
                  0.d0, &
                  0.d0, &
               elemR(iet(2),er_Topwidth), &
                  0.d0, &
                  0.d0, &
               elemR(iet(3),er_Topwidth), &
                  0.d0, &
                  0.d0, &
               elemR(iet(4),er_Topwidth),&
               elemR(iet(5),er_Topwidth),&
               elemR(iet(6),er_Topwidth), &
                  0.d0, &
                  0.d0, &
               elemR(iet(7),er_Topwidth), &
                  0.d0, &
                  0.d0, &
               elemR(iet(8),er_Topwidth), &
                  0.d0, &
                  0.d0

                  print *, ' '

         write(*,"(A,31e11.3)") ' Flowrate '  ,          &
               faceR(ift(1),fr_Flowrate), &
            elemR(iet(1),er_Flowrate), &
               faceR(ift(2),fr_Flowrate), &
               faceR(ift(2),fr_Flowrate), &
            elemR(iet(2),er_Flowrate), &
               faceR(ift(3),fr_Flowrate), &
               faceR(ift(3),fr_Flowrate), &
            elemR(iet(3),er_Flowrate), &
               faceR(ift(4),fr_Flowrate), &
               faceR(ift(4),fr_Flowrate), &
            elemR(iet(4),er_Flowrate),&
            elemR(iet(5),er_Flowrate),&
            elemR(iet(6),er_Flowrate), &
               faceR(ift(5),fr_Flowrate), &
               faceR(ift(5),fr_Flowrate), &
            elemR(iet(7),er_Flowrate), &
               faceR(ift(6),fr_Flowrate), &
               faceR(ift(6),fr_Flowrate), &
            elemR(iet(8),er_Flowrate), &
               faceR(ift(7),fr_Flowrate), &
               faceR(ift(7),fr_Flowrate)

               print *, ' '

          write(*,"(A,31f11.4)") '  Froude  '  ,          &
               faceR(ift(1),fr_FroudeNumber_d), &
           elemR(iet(1),er_FroudeNumber), &
               faceR(ift(2),fr_FroudeNumber_u), &
               faceR(ift(2),fr_FroudeNumber_d), &
           elemR(iet(2),er_FroudeNumber), &
               faceR(ift(3),fr_FroudeNumber_u), &
               faceR(ift(3),fr_FroudeNumber_d), &
           elemR(iet(3),er_FroudeNumber), &
               faceR(ift(4),fr_FroudeNumber_u), &
               faceR(ift(4),fr_FroudeNumber_d), &
           elemR(iet(4),er_FroudeNumber),&
           elemR(iet(5),er_FroudeNumber),&
           elemR(iet(6),er_FroudeNumber), &
               faceR(ift(5),fr_FroudeNumber_u), &
               faceR(ift(5),fr_FroudeNumber_d), &
           elemR(iet(7),er_FroudeNumber), &
               faceR(ift(6),fr_FroudeNumber_u), &
               faceR(ift(6),fr_FroudeNumber_d), &
           elemR(iet(8),er_FroudeNumber), &
               faceR(ift(7),fr_FroudeNumber_u), &
               faceR(ift(7),fr_FroudeNumber_d)      


         !  write(*,"(A,31e11.3)") ' CFL V    '  ,          &
         !       0.d0, &
         !    elemR(iet(1),er_Volume) / elemR(iet(1),er_Flowrate), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(2),er_Volume) / elemR(iet(2),er_Flowrate), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(3),er_Volume) / elemR(iet(3),er_Flowrate)   , &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(4),er_Volume) / elemR(iet(4),er_Flowrate)  ,&
         !    elemR(iet(5),er_Volume) / elemR(iet(6),er_Flowrate),  &
         !    elemR(iet(6),er_Volume) / elemR(iet(6),er_Flowrate), &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(7),er_Volume)  / elemR(iet(7),er_Flowrate) , &
         !       0.d0, &
         !       0.d0, &
         !    elemR(iet(8),er_Volume) / elemR(iet(8),er_Flowrate) , &
         !       0.d0, &
         !       0.d0

            ! write(*,"(A,31e11.3)") ' FlowCons '  ,          &
            !          faceR(ift(1),fr_Flowrate_Conservative), &
            !       elemR(iet(1),er_Flowrate), &
            !          faceR(ift(2),fr_Flowrate_Conservative), &
            !          faceR(ift(2),fr_Flowrate_Conservative), &
            !       elemR(iet(2),er_Flowrate), &
            !          faceR(ift(3),fr_Flowrate_Conservative), &
            !          faceR(ift(3),fr_Flowrate_Conservative), &
            !       elemR(iet(3),er_Flowrate), &
            !          faceR(ift(4),fr_Flowrate_Conservative), &
            !          faceR(ift(4),fr_Flowrate_Conservative), &
            !       elemR(iet(4),er_Flowrate),&
            !       elemR(iet(5),er_Flowrate),&
            !       elemR(iet(6),er_Flowrate), &
            !          faceR(ift(5),fr_Flowrate_Conservative), &
            !          faceR(ift(5),fr_Flowrate_Conservative), &
            !       elemR(iet(7),er_Flowrate), &
            !          faceR(ift(6),fr_Flowrate_Conservative), &
            !          faceR(ift(6),fr_Flowrate_Conservative), &
            !       elemR(iet(8),er_Flowrate), &
            !          faceR(ift(7),fr_Flowrate_Conservative), &
            !          faceR(ift(7),fr_Flowrate_Conservative)

            ! write(*,"(A,31e11.3)") ' DeltaQ   '  ,          &
            !          faceR(ift(1),fr_DeltaQ), &
            !       elemR(iet(1),er_DeltaQ), &
            !          faceR(ift(2),fr_DeltaQ), &
            !          faceR(ift(2),fr_DeltaQ), &
            !       elemR(iet(2),er_DeltaQ), &
            !          faceR(ift(3),fr_DeltaQ), &
            !          faceR(ift(3),fr_DeltaQ), &
            !       elemR(iet(3),er_DeltaQ), &
            !          faceR(ift(4),fr_DeltaQ), &
            !          faceR(ift(4),fr_DeltaQ), &
            !       elemR(iet(4),er_DeltaQ),&
            !       elemR(iet(5),er_DeltaQ),&
            !       elemR(iet(6),er_DeltaQ), &
            !          faceR(ift(5),fr_DeltaQ), &
            !          faceR(ift(5),fr_DeltaQ), &
            !       elemR(iet(7),er_DeltaQ), &
            !          faceR(ift(6),fr_DeltaQ), &
            !          faceR(ift(6),fr_DeltaQ), &
            !       elemR(iet(8),er_DeltaQ), &
            !          faceR(ift(7),fr_DeltaQ), &
            !          faceR(ift(7),fr_DeltaQ)

                     print *, ' '

            ! write(*,"(A,31e11.3)") ' Vol      '  ,          &
            !       0.d0, &
            !       elemR(iet(1),er_Volume), &
            !       0.d0, &
            !       0.d0, &
            !       elemR(iet(2),er_Volume), &
            !       0.d0, &
            !       0.d0, &
            !       elemR(iet(3),er_Volume), &
            !       0.d0, &
            !       0.d0, &
            !       elemR(iet(4),er_Volume),&
            !       elemR(iet(5),er_Volume),&
            !       elemR(iet(6),er_Volume), &
            !       0.d0, &
            !       0.d0, &
            !       elemR(iet(7),er_Volume), &
            !       0.d0, &
            !       0.d0, &
            !       elemR(iet(8),er_Volume), &
            !       0.d0, &
            !       0.d0

                  print *, ' '


         write(*,"(A,31f11.4)") '  Vel     '  ,          &
               faceR(ift(1),fr_Velocity_d), &
            elemR(iet(1),er_Velocity), &
               faceR(ift(2),fr_Velocity_u), &
               faceR(ift(2),fr_Velocity_d), &
            elemR(iet(2),er_Velocity), &
               faceR(ift(3),fr_Velocity_u), &
               faceR(ift(3),fr_Velocity_d), &
            elemR(iet(3),er_Velocity), &
               faceR(ift(4),fr_Velocity_u), &
               faceR(ift(4),fr_Velocity_d), &
            elemR(iet(4),er_Velocity),&
            elemR(iet(5),er_Velocity),&
            elemR(iet(6),er_Velocity), &
               faceR(ift(5),fr_Velocity_u), &
               faceR(ift(5),fr_Velocity_d), &
            elemR(iet(7),er_Velocity), &
               faceR(ift(6),fr_Velocity_u), &
               faceR(ift(6),fr_Velocity_d), &
            elemR(iet(8),er_Velocity), &
               faceR(ift(7),fr_Velocity_u), &
               faceR(ift(7),fr_Velocity_d)

         ! write(*,"(A,31e11.4)") ' Volume    '  ,          &
         !       0.d0, &
         !       elemR(iet(1),er_Volume), &
         !       0.d0, &
         !       0.d0, &
         !       elemR(iet(2),er_Volume), &
         !       0.d0, &
         !       0.d0, &
         !       elemR(iet(3),er_Volume), &
         !       0.d0, &
         !       0.d0, &
         !       elemR(iet(4),er_Volume),&
         !       elemR(iet(5),er_Volume),&
         !       elemR(iet(6),er_Volume), &
         !       0.d0, &
         !       0.d0, &
         !       elemR(iet(7),er_Volume), &
         !       0.d0, &
         !       0.d0, &
         !       elemR(iet(8),er_Volume), &
         !       0.d0, &
         !       0.d0

         !    ! write(*,"(A,31f11.4)") ' Froude  '  ,          &
         !    !    0.d0, &
         !    !    elemR(iet(1),er_FroudeNumber), &
         !    !    0.d0, &
         !    !    0.d0, &
         !    !    elemR(iet(2),er_FroudeNumber), &
         !    !    0.d0, &
         !    !    0.d0, &
         !    !    elemR(iet(3),er_FroudeNumber), &
         !    !    0.d0, &
         !    !    0.d0, &
         !    !    elemR(iet(4),er_FroudeNumber),&
         !    !    elemR(iet(5),er_FroudeNumber),&
         !    !    elemR(iet(6),er_FroudeNumber), &
         !    !    0.d0, &
         !    !    0.d0, &
         !    !    elemR(iet(7),er_FroudeNumber), &
         !    !    0.d0, &
         !    !    0.d0, &
         !    !    elemR(iet(8),er_FroudeNumber), &
         !    !    0.d0, &
         !    !    0.d0

               ! write(*,"(A,31f11.4)") 'wavespeed'  ,          &
               ! 0.d0, &
               ! elemR(iet(1),er_WaveSpeed), &
               ! 0.d0, &
               ! 0.d0, &
               ! elemR(iet(2),er_WaveSpeed), &
               ! 0.d0, &
               ! 0.d0, &
               ! elemR(iet(3),er_WaveSpeed), &
               ! 0.d0, &
               ! 0.d0, &
               ! elemR(iet(4),er_WaveSpeed),&
               ! elemR(iet(5),er_WaveSpeed),&
               ! elemR(iet(6),er_WaveSpeed), &
               ! 0.d0, &
               ! 0.d0, &
               ! elemR(iet(7),er_WaveSpeed), &
               ! 0.d0, &
               ! 0.d0, &
               ! elemR(iet(8),er_WaveSpeed), &
               ! 0.d0, &
               ! 0.d0

         !    !    write(*,"(A,31f11.4)") '  length '  ,          &
         !    !    0.d0, &
         !    !    elemR(iet(1),er_Length), &
         !    !    0.d0, &
         !    !    0.d0, &
         !    !    elemR(iet(2),er_Length), &
         !    !    0.d0, &
         !    !    0.d0, &
         !    !    elemR(iet(3),er_Length), &
         !    !    0.d0, &
         !    !    0.d0, &
         !    !    elemR(iet(4),er_Length),&
         !    !    elemR(iet(5),er_Length),&
         !    !    elemR(iet(6),er_Length), &
         !    !    0.d0, &
         !    !    0.d0, &
         !    !    elemR(iet(7),er_Length), &
         !    !    0.d0, &
         !    !    0.d0, &
         !    !    elemR(iet(8),er_Length), &
         !    !    0.d0, &
         !    !    0.d0

         !    !    write(*,"(A,31f11.4)") '  Qlat   '  ,          &
         !    !    0.d0, &
         !    !    elemR(iet(1),er_FlowrateLateral), &
         !    !    0.d0, &
         !    !    0.d0, &
         !    !    elemR(iet(2),er_FlowrateLateral), &
         !    !    0.d0, &
         !    !    0.d0, &
         !    !    elemR(iet(3),er_FlowrateLateral), &
         !    !    0.d0, &
         !    !    0.d0, &
         !    !    elemR(iet(4),er_FlowrateLateral),&
         !    !    elemR(iet(5),er_FlowrateLateral),&
         !    !    elemR(iet(6),er_FlowrateLateral), &
         !    !    0.d0, &
         !    !    0.d0, &
         !    !    elemR(iet(7),er_FlowrateLateral), &
         !    !    0.d0, &
         !    !    0.d0, &
         !    !    elemR(iet(8),er_FlowrateLateral), &
         !    !    0.d0, &
         !    !    0.d0

         !    write(*,"(A,31e11.4)") '  IWuQ     '  ,          &
         !       0.d0, &
         !       elemR(iet(1),er_InterpWeight_uQ), &
         !       0.d0, &
         !       0.d0, &
         !       elemR(iet(2),er_InterpWeight_uQ), &
         !       0.d0, &
         !       0.d0, &
         !       elemR(iet(3),er_InterpWeight_uQ), &
         !       0.d0, &
         !       0.d0, &
         !       elemR(iet(4),er_InterpWeight_uQ),&
         !       elemR(iet(5),er_InterpWeight_uQ),&
         !       elemR(iet(6),er_InterpWeight_uQ), &
         !       0.d0, &
         !       0.d0, &
         !       elemR(iet(7),er_InterpWeight_uQ), &
         !       0.d0, &
         !       0.d0, &
         !       elemR(iet(8),er_InterpWeight_uQ), &
         !       0.d0, &
         !       0.d0

         !    write(*,"(A,31e11.4)") '  IWdQ     '  ,          &
         !       0.d0, &
         !       elemR(iet(1),er_InterpWeight_dQ), &
         !       0.d0, &
         !       0.d0, &
         !       elemR(iet(2),er_InterpWeight_dQ), &
         !       0.d0, &
         !       0.d0, &
         !       elemR(iet(3),er_InterpWeight_dQ), &
         !       0.d0, &
         !       0.d0, &
         !       elemR(iet(4),er_InterpWeight_dQ),&
         !       elemR(iet(5),er_InterpWeight_dQ),&
         !       elemR(iet(6),er_InterpWeight_dQ), &
         !       0.d0, &
         !       0.d0, &
         !       elemR(iet(7),er_InterpWeight_dQ), &
         !       0.d0, &
         !       0.d0, &
         !       elemR(iet(8),er_InterpWeight_dQ), &
         !       0.d0, &
         !       0.d0

         !    !if (setting%Time%Step > 490) stop 550987

         !    ! write(*,"(A,31f11.4)") '  Head    '  ,          &
         !    !    faceR(ift(1),fr__d), &
         !    !    elemR(iet(1),er_), &
         !    !    faceR(ift(2),fr__u), &
         !    !    faceR(ift(2),fr__d), &
         !    !    elemR(iet(2),er_), &
         !    !    faceR(ift(3),fr__u), &
         !    !    faceR(ift(3),fr__d), &
         !    !    elemR(iet(3),er_), &
         !    !    faceR(ift(4),fr__u), &
         !    !    faceR(ift(4),fr__d), &
         !    !    elemR(iet(4),er_),&
         !    !    elemR(iet(5),er_),&
         !    !    elemR(iet(6),er_), &
         !    !    faceR(ift(5),fr__u), &
         !    !    faceR(ift(5),fr__d), &
         !    !    elemR(iet(7),er_), &
         !    !    faceR(ift(6),fr__u), &
         !    !    faceR(ift(6),fr__d), &
         !    !    elemR(iet(8),er_), &
         !    !    faceR(ift(7),fr__u), &
         !    !    faceR(ift(7),fr__d)

         

         !    ! write(*,"(A,31f11.4)") '  Head    '  ,          &
         !    ! faceR(ift(1),fr__d), &
         !    ! elemR(iet(1),er_), &
         !    ! faceR(ift(2),fr__u), &
         !    ! faceR(ift(2),fr__d), &
         !    ! elemR(iet(2),er_), &
         !    ! faceR(ift(3),fr__u), &
         !    ! faceR(ift(3),fr__d), &
         !    ! elemR(iet(3),er_), &
         !    ! faceR(ift(4),fr__u), &
         !    ! faceR(ift(4),fr__d), &
         !    ! elemR(iet(4),er_),&
         !    ! elemR(iet(5),er_),&
         !    ! elemR(iet(6),er_), &
         !    ! faceR(ift(5),fr__u), &
         !    ! faceR(ift(5),fr__d), &
         !    ! elemR(iet(7),er_), &
         !    ! faceR(ift(6),fr__u), &
         !    ! faceR(ift(6),fr__d), &
         !    ! elemR(iet(8),er_), &
         !    ! faceR(ift(7),fr__u), &
         !    ! faceR(ift(7),fr__d)

               print *, ' '


               if (setting%Time%Hydraulics%Dt < 0.13) then 
                  print *, ' '
                  print *, 'STOPPING HERE ON SMALL DT '
                  print *, ' '
                  stop 660987223
               end if
          return     

!%  
!% 5 x 6 nonjunction
!%  
   ! write(*,"(A9,30A11)"),'    ','faceU','faceD','elem','faceU','faceD','elem','faceU','faceD','elem','faceU','faceD','elem', 'faceU','faceD','elem','faceU','faceD','elem','faceU','faceD','elem'

   ! print *, ' here '
   ! print *, ift
   ! print *, iet
   ! print *, ' '
         !%  e5, f6 no junctions
         !  write(*,"(A,31f11.4)") '  Head    '  ,          &
            !       faceR(ift(1),fr_Head_u), &
            !       faceR(ift(1),fr_Head_d), &
            !    elemR(iet(1),er_Head), &
            !       faceR(ift(2),fr_Head_u), &
            !       faceR(ift(2),fr_Head_d), &
            !    elemR(iet(2),er_Head), &
            !       faceR(ift(3),fr_Head_u), &
            !       faceR(ift(3),fr_Head_d), &
            !    elemR(iet(3),er_Head),&
            !       faceR(ift(4),fr_Head_u), &
            !       faceR(ift(4),fr_Head_d), &
            !    elemR(iet(4),er_Head),&
            !       faceR(ift(5),fr_Head_u), &
            !       faceR(ift(5),fr_Head_d), &
            !    elemR(iet(5),er_Head), &
            !       faceR(ift(6),fr_Head_u), &
            !       faceR(ift(6),fr_Head_d)
   
            ! write(*,"(A,31f11.4)") '  Zbot    '  ,          &
            !       faceR(ift(1),fr_Zbottom), &
            !       faceR(ift(1),fr_Zbottom), &
            !    elemR(iet(1),er_Zbottom), &
            !       faceR(ift(2),fr_Zbottom), &
            !       faceR(ift(2),fr_Zbottom), &
            !    elemR(iet(2),er_Zbottom), &
            !       faceR(ift(3),fr_Zbottom), &
            !       faceR(ift(3),fr_Zbottom), &
            !    elemR(iet(3),er_Zbottom),&
            !       faceR(ift(4),fr_Zbottom), &
            !       faceR(ift(4),fr_Zbottom), &
            !    elemR(iet(4),er_Zbottom),&
            !       faceR(ift(5),fr_Zbottom), &
            !       faceR(ift(5),fr_Zbottom), &
            !    elemR(iet(5),er_Zbottom), &
            !       faceR(ift(6),fr_Zbottom), &
            !       faceR(ift(6),fr_Zbottom)

            ! write(*,"(A,31f11.4)") '  Depth   '  ,          &
            !       faceR(ift(1),fr_Depth_u), &
            !       faceR(ift(1),fr_Depth_d), &
            !    elemR(iet(1),er_Depth), &
            !       faceR(ift(2),fr_Depth_u), &
            !       faceR(ift(2),fr_Depth_d), &
            !    elemR(iet(2),er_Depth), &
            !       faceR(ift(3),fr_Depth_u), &
            !       faceR(ift(3),fr_Depth_d), &
            !    elemR(iet(3),er_Depth),&
            !       faceR(ift(4),fr_Depth_u), &
            !       faceR(ift(4),fr_Depth_d), &
            !    elemR(iet(4),er_Depth),&
            !       faceR(ift(5),fr_Depth_u), &
            !       faceR(ift(5),fr_Depth_d), &
            !    elemR(iet(5),er_Depth), &
            !       faceR(ift(6),fr_Depth_u), &
            !       faceR(ift(6),fr_Depth_d)

            !       print *, ' '
         
            ! write(*,"(A,31f11.4)") '  Flow    '  ,          &
            !       faceR(ift(1),fr_Flowrate), &
            !       faceR(ift(1),fr_Flowrate), &
            !    elemR(iet(1),er_Flowrate), &
            !       faceR(ift(2),fr_Flowrate), &
            !       faceR(ift(2),fr_Flowrate), &
            !    elemR(iet(2),er_Flowrate), &
            !       faceR(ift(3),fr_Flowrate), &
            !       faceR(ift(3),fr_Flowrate), &
            !    elemR(iet(3),er_Flowrate),&
            !       faceR(ift(4),fr_Flowrate), &
            !       faceR(ift(4),fr_Flowrate), &
            !    elemR(iet(4),er_Flowrate),&
            !       faceR(ift(5),fr_Flowrate), &
            !       faceR(ift(5),fr_Flowrate), &
            !    elemR(iet(5),er_Flowrate), &
            !       faceR(ift(6),fr_Flowrate), &
            !       faceR(ift(6),fr_Flowrate)


            !    write(*,"(A,31e11.4)") '  Vel     '  ,          &
            !       faceR(ift(1),fr_Velocity_u), &
            !       faceR(ift(1),fr_Velocity_d), &
            !    elemR(iet(1),er_Velocity), &
            !       faceR(ift(2),fr_Velocity_u), &
            !       faceR(ift(2),fr_Velocity_d), &
            !    elemR(iet(2),er_Velocity), &
            !       faceR(ift(3),fr_Velocity_u), &
            !       faceR(ift(3),fr_Velocity_d), &
            !    elemR(iet(3),er_Velocity),&
            !       faceR(ift(4),fr_Velocity_u), &
            !       faceR(ift(4),fr_Velocity_d), &
            !    elemR(iet(4),er_Velocity),&
            !       faceR(ift(5),fr_Velocity_u), &
            !       faceR(ift(5),fr_Velocity_d), &
            !    elemR(iet(5),er_Velocity), &
            !       faceR(ift(6),fr_Velocity_u), &
            !       faceR(ift(6),fr_Velocity_d)

         ! write(*,"(A,31f11.4)") '  Qcon    '  ,          &
         !          faceR(ift(1),fr_Flowrate_Conservative), &
         !          faceR(ift(1),fr_Flowrate_Conservative), &
         !       elemR(iet(1),er_Flowrate), &
         !          faceR(ift(2),fr_Flowrate_Conservative), &
         !          faceR(ift(2),fr_Flowrate_Conservative), &
         !       elemR(iet(2),er_Flowrate), &
         !          faceR(ift(3),fr_Flowrate_Conservative), &
         !          faceR(ift(3),fr_Flowrate_Conservative), &
         !       elemR(iet(3),er_Flowrate),&
         !          faceR(ift(4),fr_Flowrate_Conservative), &
         !          faceR(ift(4),fr_Flowrate_Conservative), &
         !       elemR(iet(4),er_Flowrate),&
         !          faceR(ift(5),fr_Flowrate_Conservative), &
         !          faceR(ift(5),fr_Flowrate_Conservative), &
         !       elemR(iet(5),er_Flowrate), &
         !          faceR(ift(6),fr_Flowrate_Conservative), &
         !          faceR(ift(6),fr_Flowrate_Conservative)  
                  
                  
         ! write(*,"(A,31e11.4)") '  iw G    '  ,          &
         !          0.d0, &
         !          elemR(iet(1),er_InterpWeight_dG), &
         !       0.d0, &
         !          elemR(iet(1),er_InterpWeight_dG), &
         !          elemR(iet(2),er_InterpWeight_uG), &
         !       0.d0, &
         !          elemR(iet(2),er_InterpWeight_dG), &
         !          elemR(iet(3),er_InterpWeight_uG), &
         !       0.d0,&
         !          elemR(iet(3),er_InterpWeight_dG), &
         !          elemR(iet(4),er_InterpWeight_uG), &
         !       0.d0,&
         !          elemR(iet(4),er_InterpWeight_dG), &
         !          elemR(iet(5),er_InterpWeight_uG), &
         !       0.d0, &
         !          elemR(iet(5),er_InterpWeight_uG), &
         !          0.d0   


            ! print *, ' '
            ! return
!%
!%  7 x 4 2 junctions
!%
      ! write(*,"(A9,30A11)"),'    ','JM','jb','faceU','faceD','elem','faceU','faceD','elem','faceU','faceD','elem', 'faceU','faceD','jb','JM','faceD'
                  !% 2 junctions on edges i=7, f = 4

            ! write(*,"(A,31f11.4)") '  Head    '  ,          &
            !    elemR(iet(1),er_Head), &
            !    elemR(iet(2),er_Head), &
            !       faceR(ift(1),fr_Head_u), &
            !       faceR(ift(1),fr_Head_d), &
            !    elemR(iet(3),er_Head), &
            !       faceR(ift(2),fr_Head_u), &
            !       faceR(ift(2),fr_Head_d), &
            !    elemR(iet(4),er_Head),&
            !       faceR(ift(3),fr_Head_u), &
            !       faceR(ift(3),fr_Head_d), &
            !    elemR(iet(5),er_Head),&
            !       faceR(ift(4),fr_Head_u), &
            !       faceR(ift(4),fr_Head_d), &
            !    elemR(iet(6),er_Head), &
            !    elemR(iet(7),er_Head)

            ! write(*,"(A,31e11.4)") '  Depth   '  ,          &
            !    elemR(iet(1),er_Depth), &
            !    elemR(iet(2),er_Depth), &
            !       faceR(ift(1),fr_Depth_u), &
            !       faceR(ift(1),fr_Depth_d), &
            !    elemR(iet(3),er_Depth), &
            !       faceR(ift(2),fr_Depth_u), &
            !       faceR(ift(2),fr_Depth_d), &
            !    elemR(iet(4),er_Depth),&
            !       faceR(ift(3),fr_Depth_u), &
            !       faceR(ift(3),fr_Depth_d), &
            !    elemR(iet(5),er_Depth),&
            !       faceR(ift(4),fr_Depth_u), &
            !       faceR(ift(4),fr_Depth_d), &
            !    elemR(iet(6),er_Depth), &
            !    elemR(iet(7),er_Depth)

            ! print *, ' '
            ! write(*,"(A,31f11.4)") '  Q       '  ,          &
            !    elemR(iet(1),er_Flowrate), &
            !    elemR(iet(2),er_Flowrate), &
            !       faceR(ift(1),fr_Flowrate), &
            !       faceR(ift(1),fr_Flowrate), &
            !    elemR(iet(3),er_Flowrate), &
            !       faceR(ift(2),fr_Flowrate), &
            !       faceR(ift(2),fr_Flowrate), &
            !    elemR(iet(4),er_Flowrate),&
            !       faceR(ift(3),fr_Flowrate), &
            !       faceR(ift(3),fr_Flowrate), &
            !    elemR(iet(5),er_Flowrate),&
            !       faceR(ift(4),fr_Flowrate), &
            !       faceR(ift(4),fr_Flowrate), &
            !    elemR(iet(6),er_Flowrate), &
            !    elemR(iet(7),er_Flowrate)

    end subroutine util_utest_CLprint    
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_utest_checkIsNan ()
      !%--------------------------------------------------------------------
      !% Description: 
      !% tests for NaN in arrays defined by eIsNan_... indexes
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
!%==========================================================================
!%
    subroutine util_utest_report_Nan (isThisNan)
         !%--------------------------------------------------------------------
         !% Description: 
         !% tests for NaN in arrays defined by eIsNan_... indexes
         !%--------------------------------------------------------------------
         !% Declarations
            logical, intent(in) :: isThisNan(Ncol_elemIsNan) 
            integer :: ii
         !%--------------------------------------------------------------------

         ! do ii = 1,Ncol_elemIsNan
         !    if (isThisNan(ii)) then 
         !       write(*,*) 'Array '
         !    end if
         ! end do

    end subroutine util_utest_report_Nan
!%
!%==========================================================================
!%==========================================================================
!%
  subroutine util_utest_local_global

   !% In this subroutine we are checking the the local and global indexs 
   !% of the link,node,elem and face arrays are unique

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

      ! !% checking elemP(:,ep_AC) indexes
      ! min_val = minval(elemP(:,ep_AC)) - 1
      ! max_val = maxval(elemP(:,ep_AC))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_AC),mask=elemP(:,ep_AC)>min_val)
      ! end do

      ! !% Here we check if nullvalueI is part packed array and if so we subtract one.
      ! !% We do this because
      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_AC) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_ac)) then
      !    print *, "ERROR:::: elemP(:,ep_AC) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_AC) is unique. This_image ::", this_image()
      ! end if


      !% checking elemP(:,ep_ALLtm) indexes

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


      ! !% checking elemP(:,ep_CC_AC) indexes

      ! min_val = minval(elemP(:,ep_CC_AC)) - 1
      ! max_val = maxval(elemP(:,ep_CC_AC))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_CC_AC),mask=elemP(:,ep_CC_AC)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if

      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_CC_AC) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_CC_AC)) then
      !    print *, "ERROR:::: elemP(:,ep_CC_AC) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_CC_AC) is unique. This_image ::", this_image()
      ! end if


      !% checking elemP(:,ep_CC) indexes

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

      ! !% checking elemP(:,ep_CC_ETM) indexes

      ! min_val = minval(elemP(:,ep_CC_ETM)) - 1
      ! max_val = maxval(elemP(:,ep_CC_ETM))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_CC_ETM),mask=elemP(:,ep_CC_ETM)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_CC_ETM) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_CC_ETM)) then
      !    print *, "ERROR:::: elemP(:,ep_CC_ETM) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_CC_ETM) is unique. This_image ::", this_image()
      ! end if

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

      !% checking elemP(:,ep_CC_Q_AC) indexes


      ! min_val = minval(elemP(:,ep_CC_Q_AC)) - 1
      ! max_val = maxval(elemP(:,ep_CC_Q_AC))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_CC_Q_AC),mask=elemP(:,ep_CC_Q_AC)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_CC_Q_AC) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_CC_Q_AC)) then
      !    print *, "ERROR:::: elemP(:,ep_CC_Q_AC) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_CC_Q_AC) is unique. This_image ::", this_image()
      ! end if

      !% checking elemP(:,ep_CC_Q) indexes

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

      ! !% checking elemP(:,ep_ETM) indexes

      ! min_val = minval(elemP(:,ep_ETM)) - 1
      ! max_val = maxval(elemP(:,ep_ETM))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_ETM),mask=elemP(:,ep_ETM)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_ETM) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_ETM)) then
      !    print *, "ERROR:::: elemP(:,ep_ETM) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_ETM) is unique. This_image ::", this_image()
      ! end if


      !% checking elemP(:,ep_JM_AC) indexes

      ! min_val = minval(elemP(:,ep_JM_AC)) - 1
      ! max_val = maxval(elemP(:,ep_JM_AC))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_JM_AC),mask=elemP(:,ep_JM_AC)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_JM_AC) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_JM_AC)) then
      !    print *, "ERROR:::: elemP(:,ep_JM_AC) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_JM_AC) is unique. This_image ::", this_image()
      ! end if


      !% checking elemP(:,ep_JM_ALLtm) indexes


      ! min_val = minval(elemP(:,ep_JM_ALLtm)) - 1
      ! max_val = maxval(elemP(:,ep_JM_ALLtm))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(elemP(:,ep_JM_ALLtm),mask=elemP(:,ep_JM_ALLtm)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "elemP(:,ep_JM_ALLtm) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_elemP(ep_JM_ALLtm)) then
      !    print *, "ERROR:::: elemP(:,ep_JM_ALLtm) is not unique. This_image ::", this_image()

      ! else
      !    print *, "elemP(:,ep_JM_ALLtm) is unique. This_image ::", this_image()
      ! end if


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

      !% checking faceP(:,fp_noBC_IorS) indexes


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

      !% checking faceP(:,fp_AC) indexes


      ! min_val = minval(faceP(:,fp_AC)) - 1
      ! max_val = maxval(faceP(:,fp_AC))
      ! ii = 0

      ! do while(min_val<max_val)
      !    ii = ii + 1
      !    min_val = minval(faceP(:,fp_AC),mask=faceP(:,fp_AC)>min_val)
      ! end do

      ! if (min_val == nullvalueI) then
      !    ii = ii - 1
      ! end if


      ! if (ii == 0 .and. min_val == nullvalueI) then
      !    print *, "faceP(:,fp_AC) is only filled with nullvalueI. This_image ::", this_image()

      ! else if (ii /= npack_faceP(fp_AC)) then
      !    print *, "ERROR:::: faceP(:,fp_AC) is not unique. This_image ::", this_image()

      ! else
      !    print *, "faceP(:,fp_AC) is unique. This_image ::", this_image()
      ! end if

      !% checking faceP(:,fp_Diag_IorS) indexes


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

      !% checking faceP(:,fp_JumpDn_IorS) indexes


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

      !% checking faceP(:,fp_JumpUp_IorS) indexes


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
