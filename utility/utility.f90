module utility

    use define_indexes
    use define_keys
    use define_globals
    use define_settings, only: setting
    use utility_crash
    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

!-----------------------------------------------------------------------------
!
! Description:
!   Utility routines that may be called in a number of places
!
!-----------------------------------------------------------------------------

    private

    public :: util_print_programheader
    public :: util_count_node_types
    public :: util_sign_with_ones
    public :: util_sign_with_ones_or_zero
    public :: util_print_warning
    public :: util_linspace
    
    public :: util_accumulate_volume_conservation
    public :: util_total_volume_conservation

    public :: util_find_elements_in_link
    public :: util_find_elements_in_junction_node
    public :: util_find_neighbors_of_CC_element
    public :: util_find_neighbors_of_JM_element

    public :: util_CLprint
    public :: util_syncwrite

    public :: util_unique_rank

    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine util_CLprint (inputstring)
        !%------------------------------------------------------------------
        !% Description:
        !% Use for command-line write during debugging
        !%------------------------------------------------------------------
        !% Declarations:
            character (len = *), intent(in) :: inputstring
            integer :: ii, nn, fD, eD, jj
            integer, pointer :: fup(:), fdn(:), eup(:), edn(:)
            real(8), pointer :: dt, oneVec(:), grav
            real(8) :: hr, aa, bb, cc

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

        print *, ' '
        write(*,"(A)") trim(inputstring)
        write(*,"(A,i7,A, f12.5, A, f12.5)") 'step = ',setting%Time%Step,'; dt = ',setting%Time%Hydraulics%Dt,'; time = ',setting%Time%Now
        print *, ' '

        !!% STUFF FOR LOOKING AT NETWORK LAYOUT
        ! do ii=1,N_elem(this_image())
        !     print *, 'ELEMENT NUMBER ',ii,' ========================================='
        !     if (elemI(ii,ei_elementType) .eq. JB)  then
        !         if (elemSI(ii,esi_JunctionBranch_Exists) .eq. 1) then
        !             print *, '    elemnt type    ',trim(reverseKey(elemI(ii,ei_elementType))), ' branch exists'
        !         else
        !             print *, '    elemnt type    ',trim(reverseKey(elemI(ii,ei_elementType))), ' dead branch'
        !         end if
        !     else
        !         print *, '    elemnt type    ',trim(reverseKey(elemI(ii,ei_elementType)))
        !         if (elemI(ii,ei_link_Gidx_SWMM) .ne. nullvalueI) then
        !             print *, '    link # ',elemI(ii,ei_link_Gidx_SWMM), ';  name = ',trim(link%Names(elemI(ii,ei_link_Gidx_SWMM))%str)
        !             print *, '    link type ',trim(reverseKey(link%I(elemI(ii,ei_link_Gidx_SWMM),li_link_type)))
        !             print *, '    link start/end elem  ',link%I(elemI(ii,ei_link_Gidx_SWMM),li_first_elem_idx), link%I(elemI(ii,ei_link_Gidx_SWMM),li_last_elem_idx)
        !         elseif ( elemI(ii,ei_node_Gidx_SWMM)) then 
        !             print *, '    node # ', elemI(ii,ei_node_Gidx_SWMM),'; name = ',trim(node%Names( elemI(ii,ei_node_Gidx_SWMM))%str)
        !             print *, '    node type ',trim(reverseKey(node%I(elemI(ii,ei_node_Gidx_SWMM),ni_node_type)))
        !         end if
        !     end if    
        ! end do

        ! stop 2098374

        ! print *, ' '

        ! do ii=1,N_elem(this_image())
        !     if ( ((elemI(ii,ei_elementType) .ne. JB) .and. (elemI(ii,ei_elementType) .ne. JM)) .or. &
        !          ((elemI(ii,ei_elementType) .eq. JB) .and. (elemSI(ii,esi_JunctionBranch_Exists) .eq. 1) )) then
        !         print *, 'ELEMENT NUMBER ',ii,' --------------------'
        !         print *, '    upstream face       ', elemI(ii,ei_MFace_uL)
        !         if ( elemI(ii,ei_MFace_uL) .ne. nullvalueI) then
        !             if (faceI(elemI(ii,ei_MFace_uL),fi_node_idx_SWMM) .ne. nullvalueI) then
        !                 print *, '    upstream node index ', faceI(elemI(ii,ei_MFace_uL),fi_node_idx_SWMM), &
        !                     '; name = ',trim(node%Names( faceI(elemI(ii,ei_MFace_uL),fi_node_idx_SWMM))%str)
        !             end if
        !         end if
        !         print *, '    dnstream face       ', elemI(ii,ei_MFace_dL)
        !         if ( elemI(ii,ei_MFace_dL) .ne. nullvalueI) then
        !             if (faceI(elemI(ii,ei_MFace_dL),fi_node_idx_SWMM) .ne. nullvalueI) then
        !                 print *, '    dnstream node index ', faceI(elemI(ii,ei_MFace_dL),fi_node_idx_SWMM),  &
        !                     '; name = ',trim(node%Names(faceI(elemI(ii,ei_MFace_dL),fi_node_idx_SWMM))%str)
        !             end if
        !         end if
        !     end if
        ! end do

        ! stop 5987341


        ! print *, 'small depth ',elemYN(iet,eYN_isSmallDepth)
        ! print *, 'zero depth  ',elemYN(iet,eYN_isZeroDepth)
        ! print *, 'is slot',elemYN(iet,eYN_isSLot)

        ! write(*,"(A,10f12.5)") 'H face                ',faceR(elemI(iet(1),ei_Mface_dL),fr_Head_d), &
        !                                                 faceR(elemI(iet(2),ei_Mface_dL),fr_Head_u), &
        !                                                 faceR(elemI(iet(3),ei_Mface_dL),fr_Head_u), &
        !                                                 faceR(elemI(iet(4),ei_Mface_dL),fr_Head_u), &
        !                                                 faceR(elemI(iet(5),ei_Mface_dL),fr_Head_u), &
        !                                                 faceR(elemI(iet(6),ei_Mface_dL),fr_Head_u)

        write(*,"(A,10e12.4)") 'n elem          ',elemR(iet(1),er_Roughness_Dynamic), &
                                                        elemR(iet(2),er_Roughness_Dynamic), &
                                                        elemR(iet(3),er_Roughness_Dynamic), &
                                                        elemR(iet(4),er_Roughness_Dynamic), &
                                                        elemR(iet(5),er_Roughness_Dynamic), &
                                                        elemR(iet(6),er_Roughness_Dynamic), &
                                                        elemR(iet(7),er_Roughness_Dynamic)                                                        

        ! write(*,"(A,10f12.5)") 'H elem          ',elemR(iet(1),er_Head), &
        !                                           elemR(iet(2),er_Head), &
        !                                           elemR(iet(3),er_Head), &
        !                                           elemR(iet(4),er_Head), &
        !                                           elemR(iet(5),er_Head), &
        !                                           elemR(iet(6),er_Head), &
        !                                           elemR(iet(7),er_Head)
     
        ! write(*,"(A,10f12.5)") 'D elem          ',elemR(iet(1),er_Depth), &
        !                                           elemR(iet(2),er_Depth), &
        !                                           elemR(iet(3),er_Depth), &
        !                                           elemR(iet(4),er_Depth), &
        !                                           elemR(iet(5),er_Depth), &
        !                                           elemR(iet(6),er_Depth), &
        !                                           elemR(iet(7),er_Depth)

        ! write(*,"(A,10f12.5)") 'slot depth      ',elemR(iet(1),er_SlotDepth), &
        !                                           elemR(iet(2),er_SlotDepth), &
        !                                           elemR(iet(3),er_SlotDepth), &
        !                                           elemR(iet(4),er_SlotDepth), &
        !                                           elemR(iet(5),er_SlotDepth), &
        !                                           elemR(iet(6),er_SlotDepth), &
        !                                           elemR(iet(7),er_SlotDepth)

        ! write(*,"(A,10f12.5)") 'slot area       ',elemR(iet(1),er_SlotArea), &
        !                                           elemR(iet(2),er_SlotArea), &
        !                                           elemR(iet(3),er_SlotArea), &
        !                                           elemR(iet(4),er_SlotArea)                                                  


        ! write(*,"(A,10f12.5)") 'V norm          ',elemR(iet(1),er_Volume) / elemR(iet(1),er_FullVolume), &
        !                                           elemR(iet(2),er_Volume) / elemR(iet(2),er_FullVolume), &
        !                                           elemR(iet(3),er_Volume) / elemR(iet(3),er_FullVolume), &
        !                                           elemR(iet(4),er_Volume) / elemR(iet(4),er_FullVolume), &
        !                                           elemR(iet(4),er_Volume) / elemR(iet(5),er_FullVolume), &
        !                                           elemR(iet(5),er_Volume) / elemR(iet(6),er_FullVolume), &
        !                                           elemR(iet(6),er_Volume) / elemR(iet(7),er_FullVolume)      
                                                  
        ! write(*,"(A,10e12.3)") 'Q elem           ',elemR(iet(1),er_Flowrate), &
        !                                           elemR(iet(2),er_Flowrate), &
        !                                           elemR(iet(3),er_Flowrate), &
        !                                           elemR(iet(4),er_Flowrate), &
        !                                           elemR(iet(5),er_Flowrate), &
        !                                           elemR(iet(6),er_Flowrate), &
        !                                           elemR(iet(7),er_Flowrate)

        ! write(*,"(A,10e12.3)") 'Q Net            ',faceR(elemI(iet(1),ei_Mface_dL),fr_Flowrate_Conservative) - faceR(elemI(iet(2),ei_Mface_dL),fr_Flowrate_Conservative), &
        !                                            faceR(elemI(iet(2),ei_Mface_dL),fr_Flowrate_Conservative) - faceR(elemI(iet(3),ei_Mface_dL),fr_Flowrate_Conservative), &
        !                                            faceR(elemI(iet(3),ei_Mface_dL),fr_Flowrate_Conservative) - faceR(elemI(iet(4),ei_Mface_dL),fr_Flowrate_Conservative), &
        !                                            faceR(elemI(iet(4),ei_Mface_dL),fr_Flowrate_Conservative) - faceR(elemI(iet(5),ei_Mface_dL),fr_Flowrate_Conservative), &
        !                                            faceR(elemI(iet(5),ei_Mface_dL),fr_Flowrate_Conservative) - faceR(elemI(iet(6),ei_Mface_dL),fr_Flowrate_Conservative), &
        !                                            faceR(elemI(iet(6),ei_Mface_dL),fr_Flowrate_Conservative) - faceR(elemI(iet(7),ei_Mface_dL),fr_Flowrate_Conservative)                                            

        ! write(*,"(A,10e12.3)") 'Qcons face            ',faceR(elemI(iet(1),ei_Mface_dL),fr_Flowrate_Conservative), &
        !                                                 faceR(elemI(iet(2),ei_Mface_dL),fr_Flowrate_Conservative), &
        !                                                 faceR(elemI(iet(3),ei_Mface_dL),fr_Flowrate_Conservative), &
        !                                                 faceR(elemI(iet(4),ei_Mface_dL),fr_Flowrate_Conservative), &
        !                                                 faceR(elemI(iet(5),ei_Mface_dL),fr_Flowrate_Conservative), &
        !                                                 faceR(elemI(iet(6),ei_Mface_dL),fr_Flowrate_Conservative)   

        ! write(*,"(A,10e12.3)") 'Q face                ',faceR(elemI(iet(1),ei_Mface_dL),fr_Flowrate), &
        !                                           faceR(elemI(iet(2),ei_Mface_dL),fr_Flowrate), &
        !                                           faceR(elemI(iet(3),ei_Mface_dL),fr_Flowrate), &
        !                                           faceR(elemI(iet(4),ei_Mface_dL),fr_Flowrate), &
        !                                           faceR(elemI(iet(5),ei_Mface_dL),fr_Flowrate), &
        !                                           faceR(elemI(iet(6),ei_Mface_dL),fr_Flowrate)                                          
!% Vasconcelos at 15 m
        ! write(*,"(A,10f12.5)") 'H elem          ',elemR(iet(1),er_Head), &
        !                                     elemR(iet(2),er_Head), &
        !                                     elemR(iet(3),er_Head), &
        !                                     elemR(iet(4),er_Head), &
        !                                     elemR(iet(5),er_Head), &
        !                                     elemR(iet(6),er_Head), &                                      
        !                                     elemR(iet(7),er_Head), &
        !                                     elemR(iet(8),er_Head)       
        
        ! write(*,"(A,3f12.5,2A,2f12.5)") 'Hface U    ',  &
        !                                     faceR(ift(1),fr_Head_u), &
        !                                     faceR(ift(2),fr_Head_u), &
        !                                     faceR(ift(3),fr_Head_u), &
        !                                     '            ', &
        !                                     '            ', &
        !                                     faceR(ift(4),fr_Head_u), &
        !                                     faceR(ift(5),fr_Head_u)

        ! print *, ' '
        ! write(*,"(A,10f12.5)") 'Q elem          ',elemR(iet(1),er_Flowrate), &
        !                                     elemR(iet(2),er_Flowrate), &
        !                                     elemR(iet(3),er_Flowrate), &
        !                                     elemR(iet(4),er_Flowrate), &
        !                                     elemR(iet(5),er_Flowrate), &
        !                                     elemR(iet(6),er_Flowrate), &   
        !                                     elemR(iet(7),er_Flowrate), &
        !                                     elemR(iet(8),er_Flowrate)        
                                            
        ! write(*,"(A,3f12.5,2A,2f12.5)") 'Qface      ',  &
        !                                     faceR(ift(1),fr_Flowrate), &
        !                                     faceR(ift(2),fr_Flowrate), &
        !                                     faceR(ift(3),fr_Flowrate), &
        !                                     '            ', &
        !                                     '            ', &
        !                                     faceR(ift(4),fr_Flowrate), &
        !                                     faceR(ift(5),fr_Flowrate)     

        ! print *, ' '                                    
        ! write(*,"(A,10f12.5)") 'Vol norm        ',elemR(iet(1),er_Volume) / elemR(iet(1),er_FullVolume), &
        !                                     elemR(iet(2),er_Volume)/ elemR(iet(2),er_FullVolume), &
        !                                     elemR(iet(3),er_Volume)/ elemR(iet(3),er_FullVolume), &
        !                                     elemR(iet(4),er_Volume)/ elemR(iet(4),er_FullVolume), &
        !                                     elemR(iet(5),er_Volume)/ elemR(iet(5),er_FullVolume), &
        !                                     elemR(iet(6),er_Volume)/ elemR(iet(6),er_FullVolume), &                                      
        !                                     elemR(iet(7),er_Volume)/ elemR(iet(7),er_FullVolume), &
        !                                     elemR(iet(8),er_Volume)/ elemR(iet(8),er_FullVolume)                                      


!% vasconcelos at 5 m

        ! write(*,"(A,10f12.5)") 'H elem          ',elemR(iet(1),er_Head), &
        !                                     elemR(iet(2),er_Head), &
        !                                     elemR(iet(3),er_Head), &
        !                                     elemR(iet(4),er_Head), &
        !                                     elemR(iet(5),er_Head), &
        !                                     elemR(iet(6),er_Head), &                                      
        !                                     elemR(iet(7),er_Head), &
        !                                     elemR(iet(8),er_Head), &                                      
        !                                     elemR(iet(9),er_Head), &
        !                                     elemR(iet(10),er_Head)       
        
        ! write(*,"(A,4f12.5,2A,3f12.5)") 'Hface U    ',  &
        !                                     faceR(ift(1),fr_Head_u), &
        !                                     faceR(ift(2),fr_Head_u), &
        !                                     faceR(ift(3),fr_Head_u), &
        !                                     faceR(ift(4),fr_Head_u), &
        !                                     '            ', &
        !                                     '            ', &
        !                                     faceR(ift(5),fr_Head_u), &
        !                                     faceR(ift(6),fr_Head_u), &
        !                                     faceR(ift(7),fr_Head_u)

        ! print *, ' '
        ! write(*,"(A,10f12.5)") 'Q elem          ',elemR(iet(1),er_Flowrate), &
        !                                     elemR(iet(2),er_Flowrate), &
        !                                     elemR(iet(3),er_Flowrate), &
        !                                     elemR(iet(4),er_Flowrate), &
        !                                     elemR(iet(5),er_Flowrate), &
        !                                     elemR(iet(6),er_Flowrate), &   
        !                                     elemR(iet(7),er_Flowrate), &
        !                                     elemR(iet(8),er_Flowrate), &   
        !                                     elemR(iet(9),er_Flowrate), &
        !                                     elemR(iet(10),er_Flowrate)        
                                            
        ! write(*,"(A,4f12.5,2A,3f12.5)") 'Qface      ',  &
        !                                     faceR(ift(1),fr_Flowrate), &
        !                                     faceR(ift(2),fr_Flowrate), &
        !                                     faceR(ift(3),fr_Flowrate), &
        !                                     faceR(ift(4),fr_Flowrate), &
        !                                     '            ', &
        !                                     '            ', &
        !                                     faceR(ift(5),fr_Flowrate), &
        !                                     faceR(ift(6),fr_Flowrate), &
        !                                     faceR(ift(7),fr_Flowrate)     

        ! print *, ' '                                    
        ! write(*,"(A,10f12.5)") 'Vol norm        ',elemR(iet(1),er_Volume) / elemR(iet(1),er_FullVolume), &
        !                                     elemR(iet(2),er_Volume)/ elemR(iet(2),er_FullVolume), &
        !                                     elemR(iet(3),er_Volume)/ elemR(iet(3),er_FullVolume), &
        !                                     elemR(iet(4),er_Volume)/ elemR(iet(4),er_FullVolume), &
        !                                     elemR(iet(5),er_Volume)/ elemR(iet(5),er_FullVolume), &
        !                                     elemR(iet(6),er_Volume)/ elemR(iet(6),er_FullVolume), &                                      
        !                                     elemR(iet(7),er_Volume)/ elemR(iet(7),er_FullVolume), &
        !                                     elemR(iet(8),er_Volume)/ elemR(iet(8),er_FullVolume), &                                      
        !                                     elemR(iet(9),er_Volume)/ elemR(iet(9),er_FullVolume), &
        !                                     elemR(iet(10),er_Volume)/ elemR(iet(10),er_FullVolume)    


!% vasconcelos at 1 m
        ! write(*,"(A,12f10.4)") 'H elem          ',elemR(iet(1),er_Head), &
        !                                     elemR(iet(2),er_Head), &
        !                                     elemR(iet(3),er_Head), &
        !                                     elemR(iet(4),er_Head), &
        !                                     elemR(iet(5),er_Head), &
        !                                     elemR(iet(6),er_Head), &                                      
        !                                     elemR(iet(7),er_Head), &
        !                                     elemR(iet(8),er_Head), &                                      
        !                                     elemR(iet(9),er_Head), &
        !                                     elemR(iet(10),er_Head), &                                      
        !                                     elemR(iet(11),er_Head), &
        !                                     elemR(iet(12),er_Head)       
        
        ! write(*,"(A,5f10.4,2A,4f10.4)") 'Hface U    ',  &
        !                                     faceR(ift(1),fr_Head_u), &
        !                                     faceR(ift(2),fr_Head_u), &
        !                                     faceR(ift(3),fr_Head_u), &
        !                                     faceR(ift(4),fr_Head_u), &
        !                                     faceR(ift(5),fr_Head_u), &
        !                                     '          ', &
        !                                     '          ', &
        !                                     faceR(ift(6),fr_Head_u), &
        !                                     faceR(ift(7),fr_Head_u), &
        !                                     faceR(ift(8),fr_Head_u), &
        !                                     faceR(ift(9),fr_Head_u)

        ! print *, ' '
        ! write(*,"(A,12f10.4)") 'Q elem          ',elemR(iet(1),er_Flowrate), &
        !                                     elemR(iet(2),er_Flowrate), &
        !                                     elemR(iet(3),er_Flowrate), &
        !                                     elemR(iet(4),er_Flowrate), &
        !                                     elemR(iet(5),er_Flowrate), &
        !                                     elemR(iet(6),er_Flowrate), &   
        !                                     elemR(iet(7),er_Flowrate), &
        !                                     elemR(iet(8),er_Flowrate), &   
        !                                     elemR(iet(9),er_Flowrate), &
        !                                     elemR(iet(10),er_Flowrate), &   
        !                                     elemR(iet(11),er_Flowrate), &
        !                                     elemR(iet(12),er_Flowrate)        
                                            
        ! write(*,"(A,5f10.4,2A,4f10.4)") 'Qface      ',  &
        !                                     faceR(ift(1),fr_Flowrate), &
        !                                     faceR(ift(2),fr_Flowrate), &
        !                                     faceR(ift(3),fr_Flowrate), &
        !                                     faceR(ift(4),fr_Flowrate), &
        !                                     faceR(ift(5),fr_Flowrate), &
        !                                     '          ', &
        !                                     '          ', &
        !                                     faceR(ift(6),fr_Flowrate), &
        !                                     faceR(ift(7),fr_Flowrate), &
        !                                     faceR(ift(8),fr_Flowrate), &
        !                                     faceR(ift(9),fr_Flowrate)     

        ! print *, ' '                                    
        ! write(*,"(A,12f10.4)") 'Vol norm        ',elemR(iet(1),er_Volume) / elemR(iet(1),er_FullVolume), &
        !                                     elemR(iet(2),er_Volume)/ elemR(iet(2),er_FullVolume), &
        !                                     elemR(iet(3),er_Volume)/ elemR(iet(3),er_FullVolume), &
        !                                     elemR(iet(4),er_Volume)/ elemR(iet(4),er_FullVolume), &
        !                                     elemR(iet(5),er_Volume)/ elemR(iet(5),er_FullVolume), &
        !                                     elemR(iet(6),er_Volume)/ elemR(iet(6),er_FullVolume), &                                      
        !                                     elemR(iet(7),er_Volume)/ elemR(iet(7),er_FullVolume), &
        !                                     elemR(iet(8),er_Volume)/ elemR(iet(8),er_FullVolume), &                                      
        !                                     elemR(iet(9),er_Volume)/ elemR(iet(9),er_FullVolume), &
        !                                     elemR(iet(10),er_Volume)/ elemR(iet(10),er_FullVolume), &                                      
        !                                     elemR(iet(11),er_Volume)/ elemR(iet(11),er_FullVolume), &
        !                                     elemR(iet(12),er_Volume)/ elemR(iet(12),er_FullVolume)    






















                                            ! write(*,"(A,3f12.5,2A,2f12.5)") 'Hface D    ',  &
        !                                     faceR(ift(1),fr_Head_d), &
        !                                     faceR(ift(2),fr_Head_d), &
        !                                     faceR(ift(3),fr_Head_d), &
        !                                     '            ', &
        !                                     '            ', &
        !                                     faceR(ift(4),fr_Head_d), &
        !                                     faceR(ift(5),fr_Head_d)                                    

        ! ! write(*,"(A,10f12.5)") 'Setting   ',elemR(iet(1),er_Setting), &
        ! !                                     elemR(iet(2),er_Setting), &
        ! !                                     elemR(iet(3),er_Setting), &
        ! !                                     elemR(iet(4),er_Setting), &
        ! !                                     elemR(iet(5),er_Setting), &
        ! !                                     elemR(iet(6),er_Setting), &                                      
        ! !                                     elemR(iet(7),er_Setting)     


        ! write(*,"(A,10f12.5)") 'Qface                       ',  &
        !                                 faceR(ift(1),fr_Flowrate), &
        !                                 faceR(ift(2),fr_Flowrate), &
        !                                 faceR(ift(3),fr_Flowrate), &
        !                                 faceR(ift(4),fr_Flowrate), &
        !                                 faceR(ift(5),fr_Flowrate), &
        !                                 faceR(ift(6),fr_Flowrate)

        ! write(*,"(A,10f12.5)") 'Vface                       ',  &
        !                                 faceR(ift(1),fr_Velocity_d), &
        !                                 faceR(ift(2),fr_Velocity_d), &
        !                                 faceR(ift(3),fr_Velocity_d)
        ! write(*,"(A,10f12.5)") 'Vface                       ',  &
        !                                 faceR(ift(1),fr_Velocity_u), &
        !                                 faceR(ift(2),fr_Velocity_u), &
        !                                 faceR(ift(3),fr_Velocity_u)

        ! write(*,"(A,10f12.5)") 'Aface                       ',  &
        !                                 faceR(ift(1),fr_Area_d), &
        !                                 faceR(ift(2),fr_Area_d), &
        !                                 faceR(ift(3),fr_Area_d)
        ! write(*,"(A,10f12.5)") 'Aface                       ',  &
        !                                 faceR(ift(1),fr_Area_u), &
        !                                 faceR(ift(2),fr_Area_u), &
        !                                 faceR(ift(3),fr_Area_u)

        ! write(*,"(A,10f12.5)") 'Hface                       ',  &
        !                                 faceR(ift(1),fr_Head_d), &
        !                                 faceR(ift(2),fr_Head_d), &
        !                                 faceR(ift(3),fr_Head_d)




        ! write(*,"(A,10f12.5)") 'Qface Conservative          ',  &
        !                                 faceR(ift(1),fr_Flowrate_Conservative), &
        !                                 faceR(ift(2),fr_Flowrate_Conservative), &
        !                                 faceR(ift(3),fr_Flowrate_Conservative), &
        !                                 faceR(ift(4),fr_Flowrate_Conservative), &
        !                                 faceR(ift(5),fr_Flowrate_Conservative), &
        !                                 faceR(ift(6),fr_Flowrate_Conservative)

        !print *, elemR(iet(3),er_Depth), elemR(iet(3),er_Flowrate)
        ! print *, '0 depth ',setting%ZeroValue%Depth
        ! print *, '0 area  ',setting%ZeroValue%Area
        ! print *, '0 volume',setting%ZeroValue%Volume

    end subroutine util_CLprint    
!%
!%==========================================================================    
!%==========================================================================
!%
    subroutine util_print_programheader ()
        if (this_image() .ne. 1) return
        write(*,*) " "
        write(*,*) "*********************************************************************"
        write(*,*) "*                            SWMM5+                                 *"
        write(*,*) "*                   gamma 0.2 unreleased                            *"
        write(*,*) "*  A public-domain, finite-volume hydraulics engine for EPA SWMM.   *"
        write(*,*) "*                        developed by NCIMM                         *"
        write(*,*) "*                                                                   *"
        write(*,*) "* NCIMM is the National Center for Infrastructure Modeling and      *"
        write(*,*) "* Management funded under US EPA Cooperative Agreement 83595001     *" 
        write(*,*) "* awarded to the University of Texas at Austin, 2017-23.            *"
        write(*,*) "* PI: Prof. Ben R. Hodges                                           *"
        write(*,*) "*                                                                   *"
        write(*,*) "* Code authors:                                                     *"
        write(*,*) "*    2016-2022 Dr. Ben R. Hodges                                    *"
        write(*,*) "*    2016-2022 Dr. Edward Tiernan                                   *"
        write(*,*) "*    2019-2022 Sazzad Sharior                                       *"
        write(*,*) "*    2019-2022 Gerardo Riano-Briceno                                *"
        write(*,*) "*    2019-2022 Eric Jenkins                                         *"
        write(*,*) "*    2020-2022 Dr. Cheng-Wei (Justin) Yu                            *"
        write(*,*) "*    2021-2022 Christopher Brashear                                 *"
        write(*,*) "*    2021-2022 Abdulmuttalib Lokhandwala                            *"
        write(*,*) "*    2022-2022 Cesar E. Davila Hernandez                            *"
        write(*,*) "*    2018-2020 Dr. Ehsan Madadi-Kandjani                            *"
        write(*,*) "*********************************************************************"
        write(*,*)
        write(*,"(A,i5,A)") "Simulation starts with ",num_images()," processors"
        write(*,*) ' '
        write(*,"(A)") 'Using the following files and folders:'
        write(*,"(A)") '  SWMM Input file   : '//trim(setting%File%inp_file)
        write(*,"(A)") '  SWMM Report file  : '//trim(setting%File%rpt_file)
        write(*,"(A)") '  SWMM Output file  : '//trim(setting%File%out_file)
        write(*,"(A)") '  Output folder     : '//trim(setting%File%output_folder)
        write(*,"(A)") '  Settings file     : '//trim(setting%File%setting_file)
        write(*,"(A)") '  Library folder    : '//trim(setting%File%library_folder)
        write(*,"(A)") '  Project folder    : '//trim(setting%File%project_folder)
    end subroutine util_print_programheader  
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_count_node_types(N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2, N_nJ1)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% This subroutine uses the vectorized count() function to search the array for
        !% numberof instances of each node type
        !%-----------------------------------------------------------------------------
        integer, intent(in out) :: N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2, N_nJ1
        integer :: ii
        !%-----------------------------------------------------------------------------
        N_nBCup = count(node%I(:, ni_node_type) == nBCup)
        N_nBCdn = count(node%I(:, ni_node_type) == nBCdn)
        N_nJm = count(node%I(:, ni_node_type) == nJM)
        N_nStorage = count(node%I(:, ni_node_type) == nStorage)
        N_nJ2 = count(node%I(:, ni_node_type) == nJ2)
        N_nj1 = count(node%I(:, ni_node_type) == nJ1)

    end subroutine util_count_node_types
!%
!%==========================================================================
!%==========================================================================
!%
    pure elemental real(8) function util_sign_with_ones &
        (inarray) result (outarray)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% returns is an array of real ones with the sign of the inarray argument
        !%-----------------------------------------------------------------------------
        real(8),      intent(in)    :: inarray
        !%-----------------------------------------------------------------------------
        outarray = oneR
        outarray = sign(outarray,inarray)

    end function util_sign_with_ones

!%
!%==========================================================================
!%==========================================================================
!%
    function util_sign_with_ones_or_zero &
        (inarray) result (outarray)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% returns is an array of real ones with the sign of the inarray argument
        !% for non-zero inarray. For zero inarray returns zero.
        !%-----------------------------------------------------------------------------
        real(8),intent(in) :: inarray(:)
        real(8)            :: outarray(size(inarray,1))
        !%-----------------------------------------------------------------------------
        outarray = oneR
        outarray = sign(outarray,inarray)

        !% brh 20211227 
        where(inarray == zeroR)
            outarray = zeroI
        end where

    end function util_sign_with_ones_or_zero

!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_print_warning(msg,async)
        !%------------------------------------------------------------------
        !% Used for opening up the warning files and writing to the file
        !%------------------------------------------------------------------
            character(len = *), intent(in) :: msg
            logical, optional, intent(in) :: async
            logical :: async_actual
        !%------------------------------------------------------------------
        if (present(async)) then
            async_actual = async
        else
            async_actual = .true.
        end if
        if (this_image() == 1) then
            print *, "Warning: "//trim(msg)
        else if (async_actual) then
            print *, "Warning: "//trim(msg)
        end if

    end subroutine util_print_warning
!%
!%==========================================================================
!%==========================================================================
!%
    function util_linspace(startPoint,endPoint,N) result(outArray)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% similar to python/matlab linspace
        !%-----------------------------------------------------------------------------
        real(8), intent(in)  :: startPoint 
        real(8), intent(in)  :: endPoint
        integer, intent(in)  :: N
        real(8)              :: delta
        real(8), allocatable :: outArray(:)
        integer :: ii
        !%-----------------------------------------------------------------------------
        !% calculate step size
        delta = (endPoint - startPoint)/real(N-1,8)

        !% allocate the outArry based on number of samples
        allocate(outArray(N))

        do ii = 1, N
            outArray(ii) = startPoint + (ii-1)*delta
        end do

    end function util_linspace
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine util_accumulate_volume_conservation ()
        !%------------------------------------------------------------------
        !% Description:
        !% Computes/stores the cumulative mass conservation for each CC and JM
        !% element that is NOT small volume or zero depth
        !% HACK -- JM not completed
        !% HACK -- need a separate volume cons for the small/zero losses
        !%------------------------------------------------------------------
        !% Declarations:
            real(8), pointer :: eCons(:), fQ(:), eQLat(:), VolNew(:), VolOld(:), dt
            real(8), pointer :: VolOver(:), VolSlot(:)
            integer, pointer :: thisColCC, thisColJM, npack, thisP(:)
            integer, pointer :: fdn(:), fup(:), BranchExists(:)
            integer :: ii,kk
            real(8) :: netQ
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:
            thisColCC => col_elemP(ep_CC_ALLtm)
            thisColJM => col_elemP(ep_JM_ALLtm)
            fQ      => faceR(:,fr_Flowrate_Conservative)
            eCons   => elemR(:,er_VolumeConservation)
            eQLat   => elemR(:,er_FlowrateLateral)
            VolNew  => elemR(:,er_Volume)  
            VolOld  => elemR(:,er_Volume_N0) 
            VolOver => elemR(:,er_VolumeOverFlow) 
            VolSlot => elemR(:,er_SlotVolume)
            fup     => elemI(:,ei_Mface_uL)
            fdn     => elemI(:,ei_Mface_dL)
            dt      => setting%Time%Hydraulics%Dt
            BranchExists => elemSI(:,esi_JunctionBranch_Exists)
        !%------------------------------------------------------------------
        !% --- for the CC elements
        npack   => npack_elemP(thisColCC)

        ! print *
        ! print *, 'in here at top'
        !print *, 'aaa ',elemR(48,er_VolumeOverFlow)

        if (npack > 0) then
            thisP => elemP(1:npack,thisColCC)

            !where (.not. elemYN(thisP,eYN_isZeroDepth))
            !% --- sum of the net inflow and lateral flow should be the change in volume from head
            !% --- for output, use an accumulator
            eCons(thisP) = eCons(thisP)                                            &
                         + dt * ( fQ(fup(thisP)) - fQ(fdn(thisP)) + eQlat(thisP) ) &
                         - (VolNew(thisP) - VolOld(thisP))  - VolOver(thisP)

            !endwhere             
            !% --- for debugging, switch to using non-cumulative            
            !eCons(thisP) = dt * (fQ(fup(thisP)) - fQ(fdn(thisP)) + eQlat(thisP)) &
            !              - (VolNew(thisP) - VolOld(thisP))  - VolOver(thisP)      
                          
            do ii = 1,size(thisP)
                if (abs(eCons(thisP(ii))) > 1.0e-4) then
                    print *, 'CONSERVATION ISSUE CC', ii, thisP(ii), this_image()
                    print *,  eCons(thisP(ii))
                    print *, fQ(fup(thisP(ii))), fQ(fdn(thisP(ii))), eQlat(thisP(ii))
                    print *, ' vol ',VolNew(thisP(ii)), VolOld(thisP(ii))
                    print *, VolNew(thisP(ii)) - VolOld(thisP(ii)) & 
                            - dt * fQ(fup(thisP(ii))) + dt * fQ(fdn(thisP(ii))) - dt * eQlat(thisP(ii)) 
                    print *, 'vol overflow ',elemR(thisP(ii),er_VolumeOverFlow)
                    call util_crashpoint(358783)
                end if
            end do  

        end if

        !% for the JM elements
        npack   => npack_elemP(thisColJM)
        if (npack > 0) then
            thisP => elemP(1:npack,thisColJM)

            ! print *, 'printing stuff'
            ! do ii=1,npack
            !     if (thisP(ii) == 47) then
            !         print *, thisP(ii), eCons(thisP(ii)), eQlat(thisP(ii))
            !         print *, VolNew(thisP(ii)), VolOld(thisP(ii))
            !         print *, VolOver(thisP(ii)), VolSlot(thisP(ii))
            !     end if
            ! end do
            ! print *, 'xxx   ',eCons(47)
            ! print *, 'Vlat  ',dt * eQlat(47)
            ! print *, 'Vnew  ',VolNew(47)
            ! print *, 'Vold  ',VolOld(47)
            ! print *, 'Vover ',VolOver(47)
            ! print *, 'Vslt  ',VolSlot(47)
            ! print *, 'inV1  ',dt * fQ(fup(48))
            ! print *, 'inV2  ',dt * fQ(fup(50))
            ! print *, 'ouV1  ',dt * fQ(fdn(49))
            ! print *, ' '
            ! print *, 'eCons ',eCons(47)
            ! print *, ' '

            eCons(thisP) = eCons(thisP) + dt * eQlat(thisP) &
                    - (VolNew(thisP) - VolOld(thisP)) - VolOver(thisP) + VolSlot(thisP)

            ! do ii=1,npack
            !     if (thisP(ii) == 47) then
            !         print *, 'eCons ',thisP(ii),eCons(thisP(ii))
            !     end if
            ! end do


            !% debug, compute just this time step conservation
            !eCons(thisP) = dt * eQlat(thisP) - (VolNew(thisP) - VolOld(thisP)) - VolOver(thisP)

            ! print *, 'yyy ',eCons(47)

            ! print *, 'Vol Over ', VolOver(48), VolOver(49)
            ! print *, dt * (fQ(fup(48)) + fQ(fup(50)) - fQ(fdn(49)) )

            do ii=1,max_branch_per_node,2
                eCons(thisP) = eCons(thisP)                                                    &
                             + dt * ( real(BranchExists(thisP+ii  ),8) * fQ(fup(thisP+ii))     &
                                    - real(BranchExists(thisP+ii+1),8) * fQ(fdn(thisP+ii+1)) ) &
                             - VolOver(thisP+ii  ) * real(BranchExists(thisP+ii  ),8)          &
                             - VolOver(thisP+ii+1) * real(BranchExists(thisP+ii+1),8)                              
            end do
            
            ! print *, 'zzz ', eCons(47)

            do ii = 1,size(thisP)
                if (abs(eCons(thisP(ii))) > 1.0e-4) then
                    print *, ' '
                    print *, 'CONSERVATION ISSUE JM', ii, thisP(ii), this_image()
                    print *, 'is zero =',elemYN(thisP(ii),eYN_isZeroDepth), ';   is smalldepth = ',elemYN(thisP(ii),eYN_isSmallDepth)
                    print *,  'volume overflow ', VolOver(thisP(ii))
                    print *,  'net cons ',eCons(thisP(ii))
                    do kk = 1,max_branch_per_node,2
                        print *, 'branch Q',kk,   fQ(fup(thisP(ii)+kk  )) * real(BranchExists(thisP(ii)+kk  ),8)   &
                                           ,kk+1, fQ(fdn(thisP(ii)+kk+1)) * real(BranchExists(thisP(ii)+kk+1),8) 
                    end do
                    do kk = 1,max_branch_per_node,2
                        print *, 'branch VolOver',thisP(ii)+kk,   VolOver(thisP(ii)+kk  ) * real(BranchExists(thisP(ii)+kk  ),8)   &
                                                 ,thisP(ii)+kk+1, VolOver(thisP(ii)+kk+1) * real(BranchExists(thisP(ii)+kk+1),8) 
                    end do
                    print *, ' vol ', VolNew(thisP(ii)), VolOld(thisP(ii)), VolSlot(thisP(ii))
                    print *, 'd vol' , (VolNew(thisP(ii)) - VolOld(thisP(ii))) + VolSlot(thisP(ii))
                    netQ = eQlat(thisP(ii))
                    do kk = 1,max_branch_per_node,2
                        netQ = netQ + fQ(fup(thisP(ii)+kk))  * real(BranchExists(thisP(ii)+kk  ),8)  &
                                     -fQ(fdn(thisP(ii)+kk+1))* real(BranchExists(thisP(ii)+kk+1),8)
                    end do
                    !print *, ' netQ ',netQ
                    print *, ' netVol ', netQ * dt

                    
                    call util_crashpoint(3587832)
                end if
            end do   

        end if

        !%------------------------------------------------------------------

    end subroutine util_accumulate_volume_conservation
!%==========================================================================
!%==========================================================================
!%
    subroutine util_total_volume_conservation (volume_nonconservation)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes total volume non-conservation of all time-marching elements
        !%------------------------------------------------------------------
        !% Declarations:
            real(8), intent(inout) :: volume_nonconservation
            integer, pointer       :: npack, thisP(:), thisCol
            real(8), save          :: vstore[*]
            integer :: ii
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------
    
        !% for CC elements of any TM
        vstore = zeroR
        thisCol => col_elemP(ep_CC_ALLtm)
        npack   => npack_elemP(thisCol)
        if (npack > 0) then
            thisP => elemP(1:npack,thisCol)
            vstore = vstore + sum(elemR(thisP,er_VolumeConservation))
            !do ii = 1,npack
            !    if (abs(elemR(thisP(ii),er_VolumeConservation)) > 1.0) then
            !        print *, thisP(ii), elemR(thisP(ii),er_VolumeConservation)
            !    end if
            !end do
        end if
       ! print *, 'in util total_volume conservation ',vstore, this_image()

        !% for JM ETM elements
        thisCol =>col_elemP(ep_JM_ALLtm) 
        npack   => npack_elemP(thisCol)
        if (npack > 0) then
            thisP => elemP(1:npack,thisCol)
            vstore = vstore + sum(elemR(thisP,er_VolumeConservation))
            !do ii = 1,npack
            !    if (abs(elemR(thisP(ii),er_VolumeConservation)) > 1.0) then
            !        print *, thisP(ii), elemR(thisP(ii),er_VolumeConservation)
            !    end if
            !end do
        end if

        sync all
        call co_sum(vstore, result_image=1)

        volume_nonconservation = vstore
    
        !%------------------------------------------------------------------
        !% Closing:
    end subroutine util_total_volume_conservation
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_find_elements_in_link &
        (thislinkname, thislink_idx, thislink_image, elemInLink, nElemInLink)
        !%-------------------------------------------------------------------
        !% Description: 
        !% Finds the element indexes for the SWMM link with the name "thislinkname"
        !% stores the element indexes in the array elemInLink and
        !% returns the link index in link%I(:) and the image that host the link.
        !% Not efficiently written, but only purpose is for use in debugging.
        !%-------------------------------------------------------------------
            character(*), intent(in)   :: thislinkname
            integer, intent(inout)     :: elemInLink(:)
            integer, intent(inout)     :: thislink_idx, thislink_image
            integer, intent(inout)     :: nElemInLink
            integer ::  ii, jj
        !%-------------------------------------------------------------------
        elemInLink = 0
        thislink_idx = 0
        thislink_image = 0

        !% find the link idx and image for the input "thislinkname"
        if (this_image() == 1) then
            do ii = 1,size(link%I,dim=1)
                if (link%Names(ii)%str == thislinkname) then
                    write(*,"(A,A,A,i8,A,i6)")'link name ', trim(link%Names(ii)%str), ';  linkIdx= ', ii, ' On image = ',link%I(ii,li_P_image)
                    thislink_idx = ii
                    thislink_image = link%I(ii,li_P_image)
                end if
            end do
        end if
        !% broadcast result to all images
        call co_broadcast(thislink_idx,  source_image=1)
        call co_broadcast(thislink_image,source_image=1)
        sync all

    
        !% for the image that hosts the link, cycle through to find the elements in the link
        if (this_image() == thislink_image) then
            nElemInLink = 0
            do ii = 1,N_elem(this_image())
                if (elemI(ii,ei_link_Gidx_SWMM) == thislink_idx) then
                    nElemInLink = nElemInLink + 1
                    elemInLink(nElemInLink) = ii
                end if
            end do
            write(*,"(A,100i8)") 'elemIdx =',elemInLink(1:nElemInLink)
        end if
        call co_broadcast(elemInLink,source_image=thislink_image)
        sync all

    end subroutine util_find_elements_in_link
!%   
!%==========================================================================    
!%==========================================================================
!%
    subroutine util_find_elements_in_junction_node &
        (thisnodename, thisnode_idx, thisnode_image, elemJM_idx)
        !%------------------------------------------------------------------
        !% Description:
        !% Given the node name, this finds the node index, the image on
        !% which it is an element, and the JM element index
        !%
        !%------------------------------------------------------------------
        !% Declarations:
            character(*), intent(in)   :: thisnodename
            integer, intent(inout)     :: thisnode_idx, thisnode_image
            integer, intent(inout)     :: elemJM_idx
            integer ::  ii, jj
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------

        thisnode_idx = 0
        thisnode_image = 0

        !% find the node idx and image for the input "thisnodename"
        if (this_image() == 1) then
            do ii = 1,size(node%I,dim=1)
                !print *, ii, trim(node%Names(ii)%str)
                if (node%Names(ii)%str == thisnodename) then
                    write(*,"(A,A,A,i8,A,i6)")'node name ', trim(node%Names(ii)%str), ';  nodeIdx= ', ii, '; On image = ',node%I(ii,ni_P_image)
                    thisnode_idx = ii
                    thisnode_image = node%I(ii,ni_P_image)
                    exit !% the first node found should be the JM
                end if
            end do
        end if
        !% broadcast result to all images
        call co_broadcast(thisnode_idx,  source_image=1)
        call co_broadcast(thisnode_image,source_image=1)
        sync all

        elemJM_idx =0
        !% for the image that hosts the node, cycle through to find the element
        if (this_image() == thisnode_image) then
            do ii = 1,N_elem(this_image())
                !print *, ii, elemI(ii,ei_node_Gidx_SWMM)
                if (elemI(ii,ei_node_Gidx_SWMM) == thisnode_idx) then
                    write (*,"(A,i8,A,A)") 'elemIdx= ',ii,'; type = ',reverseKey(elemI(ii,ei_elementType))
                    elemJM_idx = ii
                    exit !% the first should be the correct nodes
                end if
            end do
        end if
        call co_broadcast(elemJM_idx,source_image=thisnode_image)
        sync all

        if (elemJM_idx == 0) then
            write(*,"(A,A)") 'No corresponding JM element found for the node ',trim(node%Names(thisnode_idx)%str)
            write(*,"(A)") 'This implies the node is a nJ1 or nJ2 and is either a boundary or a face'
        end if

        !%------------------------------------------------------------------
        !% Closing:
    end subroutine util_find_elements_in_junction_node
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_find_neighbors_of_CC_element (eIdx, iset)
        !%------------------------------------------------------------------
        !% Description:
        !% For debugging, given an element ID that is CC, find the upstream faces
        !% and elements. Output is in the form
        !% (upstream element, upstream face, element, downstream face, downstream element)\
        !% Note that this will not work across shared faces
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: eIdx
            integer, intent(inout) :: iset(5)
            integer :: ifaceUp, ifaceDn, ielemUp, ielemDn
            character(64) :: subroutine_name = 'util_find_neighbors_of_CC_element'
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (.not. (elemI(eIdx,ei_elementType) == CC)) then
                write(*,"(A,A,A,i8,A,A)") 'in ',trim(subroutine_name), ': element ',eIdx, 'is of type ',reverseKey(elemI(eIdx,ei_elementType))
                write(*,"(A)") 'but this procedure requires CC. Exiting with no result.'
                return
            end if
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------

        ifaceUp = elemI(eIdx,ei_Mface_uL)
        ifaceDn = elemI(eIdx,ei_Mface_dL)
        ielemUp = faceI(ifaceUp,fi_Melem_uL)
        ielemDn = faceI(ifaceDn,fi_Melem_dL)

        iset(1) = ielemUp
        iset(2) = iFaceUp
        iset(3) = eIdx
        iset(4) = ifaceDn
        iset(5) = ielemDn
    
        !%------------------------------------------------------------------
        !% Closing:
    end subroutine util_find_neighbors_of_CC_element
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_find_neighbors_of_JM_element &
        (eIdx, iUpSet, iDnSet, nUpBranch, nDnBranch)
        !%------------------------------------------------------------------
        !% Description:
        !% For debugging, given an element ID of type JM this finds the
        !% faces and neighbor elements
        !% Note that this will not work across shared faces
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: eIdx
            integer, intent(inout) :: iUpSet(max_up_branch_per_node,4)
            integer, intent(inout) :: iDnSet(max_dn_branch_per_node,4)
            integer, intent(inout) :: nUpBranch, nDnBranch
            integer, dimension(max_up_branch_per_node) :: ifaceUp, ielemUp, eIdx_BranchUp
            integer, dimension(max_dn_branch_per_node) :: ifaceDn, ielemDn, eIdx_BranchDn
            integer :: ii
            character(64) :: subroutine_name = 'util_find_neighbors_of_JM_element'
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (.not. (elemI(eIdx,ei_elementType) == JM)) then
                write(*,"(A,A,A,i8,A,A)") 'in ',trim(subroutine_name), ': element ',eIdx, 'is of type ',reverseKey(elemI(eIdx,ei_elementType))
                write(*,"(A)") 'but this procedure requires JM. Exiting with no result.'
                return
            end if
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------

        iUpSet = 0
        iDnSet = 0

        nUpBranch = 0
        do ii=1,max_branch_per_node,2
            if (elemSI(eIdx+ii,esi_JunctionBranch_Exists) > 0) then
                nUpBranch = nUpBranch + 1
                eIdx_BranchUp(nUpBranch) = eIdx + ii
                ifaceUp(nUpBranch) = elemI(eIdx_BranchUp(nUpBranch),ei_Mface_uL)
                ielemUp(nUpBranch) = faceI(ifaceUp(nUpBranch),fi_Melem_uL)
            end if
        end do    

        nDnBranch = 0
        do ii=2,max_branch_per_node,2
            if (elemSI(eIdx+ii,esi_JunctionBranch_Exists) > 0) then
                nDnBranch = nDnBranch + 1
                eIdx_BranchDn(nDnBranch) = eIdx + ii
                ifaceDn(nDnBranch) = elemI(eIdx_BranchDn(nDnBranch),ei_Mface_dL)
                ielemDn(nDnBranch) = faceI(ifaceUp(nDnBranch),fi_Melem_dL)
            end if
        end do    

        do ii=1,nUpBranch
            iUpSet(ii,1) = ielemUp(ii)
            iUpSet(ii,2) = ifaceUp(ii)
            iUpset(ii,3) = eIdx_BranchUp(ii)
            iUpSet(ii,4) = eIdx
        end do

        do ii=1,nDnBranch
            iDnSet(ii,1) = eIdx
            iDnset(ii,2) = eIdx_BranchDn(ii)
            iDnSet(ii,3) = ifaceDn(ii)
            iDnSet(ii,4) = ielemDn(ii)
        end do
            
    
        !%------------------------------------------------------------------
        !% Closing:
    end subroutine util_find_neighbors_of_JM_element
!%
!%==========================================================================    
!%==========================================================================
!%
    subroutine util_syncwrite
        !%------------------------------------------------------------------
        !% Description:
        !% writes the coarray outstring in processor order
        !% this is useful for debugging when an image is hanging.
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: ii
            character (len = 256) :: tstring(num_images())
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------

        sync all
        if (this_image() == 1) then
            tstring(1) = trim(outstring)
            do ii=2,num_images()
                tstring(ii) = trim(outstring[ii])
                ! sync all
                ! flush(6)
                ! if (ii==this_image()) then
                !     write(6,*) trim(thisstring),this_image(), thisI
                !     flush(6)
                !     !wait(6)
                !     !call sleep(1)
                ! end if
                ! sync all
            end do
            do ii=1,num_images()
                write(6,*) trim(tstring(ii)),ii
            end do
        end if
        sync all

    end subroutine util_syncwrite      
!%
!%==========================================================================    
!%==========================================================================
!%
    subroutine util_unique_rank (xInput, xRanked, Nunique)
        !%------------------------------------------------------------------
        !% Description:
        !% input array (integer) xInput is ranked (small to large)
        !% and duplicates discarded.
        !% Modified from public domain code ORDERPACK 2.0 
        !% written by Michel Olagnon, IFREMER Brest, michel.olagnon@ifremer.fr
        !% accessed from www.fortran-2000.com/rank/ on June 21, 2022
        !% Subroutine I_UNIRNK from module UNIRNK extracted and modified below
        !% The approach uses Merge-sort ranking of an array, with removal of
        !   duplicate entries.
        !   The routine is similar to pure merge-sort ranking, but on
        !   the last pass, it discards indices that correspond to
        !   duplicate entries.
        !   For performance reasons, the first 2 passes are taken
        !   out of the standard loop, and use dedicated coding.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, dimension (:), intent (in)  :: xInput
            integer, dimension (:), intent (out) :: xRanked
            integer, intent (out) :: Nunique

            integer, Dimension (SIZE(xRanked)) :: jwRankT
            integer :: LmtnA, LmtnC, iRanking, iRanking1, iRanking2
            integer :: nval, iInd, iwRankD, iwRank, iwRankF, jIndA, iIndA, iIndB
            integer :: xTst, xValA, xValB
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:   
        !%------------------------------------------------------------------
        nval    = Min (SIZE(xInput), SIZE(xRanked))
        Nunique = nval
    
        select case (nval)
        case (:0)
            return
        case (1)
            xRanked (1) = 1
            return
        case default
            continue
        end select
        !%
        !%  Fill-in the index array, creating ordered couples
        !%
        do iInd = 2, nval, 2
            if (xInput(iInd-1) < xInput(iInd)) then
                xRanked (iInd-1) = iInd - 1
                xRanked (iInd) = iInd
            else
                xRanked (iInd-1) = iInd
                xRanked (iInd) = iInd - 1
            end if
        end Do
        if (Modulo(nval, 2) /= 0) then
            xRanked (nval) = nval
        end if
        !%
        !%  We will now have ordered subsets A - B - A - B - ...
        !%  and merge A and B couples into     C   -   C   - ...
        !%
        LmtnA = 2
        LmtnC = 4
        !%
        !%  First iteration. The length of the ordered subsets goes from 2 to 4
        !%
        do
            if (nval <= 4) exit
            !%
            !%   Loop on merges of A and B into C
            !%
            do iwRankD = 0, nval - 1, 4
                if ((iwRankD+4) > nval) then
                    if ((iwRankD+2) >= nval) Exit
                    !%
                    !%   1 2 3
                    !%
                    if (xInput(xRanked(iwRankD+2)) <= xInput(xRanked(iwRankD+3))) exit
                    !%
                    !%   1 3 2
                    !%
                    if (xInput(xRanked(iwRankD+1)) <= xInput(xRanked(iwRankD+3))) then
                        iRanking2 = xRanked (iwRankD+2)
                        xRanked (iwRankD+2) = xRanked (iwRankD+3)
                        xRanked (iwRankD+3) = iRanking2
                    !%
                    !%   3 1 2
                    !%
                    else
                        iRanking1 = xRanked (iwRankD+1)
                        xRanked (iwRankD+1) = xRanked (iwRankD+3)
                        xRanked (iwRankD+3) = xRanked (iwRankD+2)
                        xRanked (iwRankD+2) = iRanking1
                    end if
                    exit
                end if
                !%
                !%   1 2 3 4
                !%
                if (xInput(xRanked(iwRankD+2)) <= xInput(xRanked(iwRankD+3))) cycle
                !%
                !%   1 3 x x
                !%
                if (xInput(xRanked(iwRankD+1)) <= xInput(xRanked(iwRankD+3))) then
                    iRanking2 = xRanked (iwRankD+2)
                    xRanked (iwRankD+2) = xRanked (iwRankD+3)
                    if (xInput(iRanking2) <= xInput(xRanked(iwRankD+4))) then
                        !%   1 3 2 4
                        xRanked (iwRankD+3) = iRanking2
                    else
                        !%   1 3 4 2
                        xRanked (iwRankD+3) = xRanked (iwRankD+4)
                        xRanked (iwRankD+4) = iRanking2
                    end if
                !%
                !%   3 x x x
                !%
                else
                    iRanking1 = xRanked (iwRankD+1)
                    iRanking2 = xRanked (iwRankD+2)
                    xRanked (iwRankD+1) = xRanked (iwRankD+3)
                    if (xInput(iRanking1) <= xInput(xRanked(iwRankD+4))) then
                        xRanked (iwRankD+2) = iRanking1
                        if (xInput(iRanking2) <= xInput(xRanked(iwRankD+4))) then
                            !%   3 1 2 4
                            xRanked (iwRankD+3) = iRanking2
                        else
                            !%   3 1 4 2
                            xRanked (iwRankD+3) = xRanked (iwRankD+4)
                            xRanked (iwRankD+4) = iRanking2
                        end If
                    else
                        !%   3 4 1 2
                        xRanked (iwRankD+2) = xRanked (iwRankD+4)
                        xRanked (iwRankD+3) = iRanking1
                        xRanked (iwRankD+4) = iRanking2
                    end if
                end if
            end do
            !%
            !%  The Cs become As and Bs
            !%
            LmtnA = 4
            exit
        end do
        !%
        !%  Iteration loop. Each time, the length of the ordered subsets
        !%  is doubled.
        !%
        do
            if (2*LmtnA >= nval) exit
            iwRankF = 0
            LmtnC = 2 * LmtnC
            !%
            !%   Loop on merges of A and B into C
            !%
            do
                iwRank  = iwRankF
                iwRankD = iwRankF + 1
                jIndA   = iwRankF + LmtnA
                iwRankF = iwRankF + LmtnC
                if (iwRankF >= nval) then
                    if (jIndA >= nval) exit
                    iwRankF = nval
                end if
                iIndA = 1
                iIndB = jIndA + 1
                !%
                !%  One steps in the C subset, that we create in the final rank array
                !%
                !%  Make a copy of the rank array for the iteration
                !%
                jwRankT (1:LmtnA) = xRanked (iwRankD:jIndA)
                xValA = xInput (jwRankT(iIndA))
                xValB = xInput (xRanked(iIndB))
                !%
                do
                    iwRank = iwRank + 1
                    !%
                    !%  We still have unprocessed values in both A and B
                    !%
                    if (xValA > xValB) then
                        xRanked (iwRank) = xRanked (iIndB)
                        iIndB = iIndB + 1
                        if (iIndB > iwRankF) then
                            !%  Only A still with unprocessed values
                            xRanked (iwRank+1:iwRankF) = jwRankT (iIndA:LmtnA)
                            exit
                        end if
                        xValB = xInput (xRanked(iIndB))
                    else
                        xRanked (iwRank) = jwRankT (iIndA)
                        iIndA = iIndA + 1
                        if (iIndA > LmtnA) exit! Only B still with unprocessed values
                        xValA = xInput (jwRankT(iIndA))
                    end if
                end do
            end do
            !%
            !%  The Cs become As and Bs
            !%
            LmtnA = 2 * LmtnA
        end Do
        !%
        !%   Last merge of A and B into C, with removal of duplicates.
        !%
        iIndA = 1
        iIndB = LmtnA + 1
        Nunique = 0
        !%
        !%  One steps in the C subset, that we create in the final rank array
        !%
        jwRankT (1:LmtnA) = xRanked (1:LmtnA)
        if (iIndB <= nval) then
            xTst = Min(xInput(jwRankT(1)), xInput(xRanked(iIndB))) - 1
        else
            xTst = xInput(jwRankT(1)) - 1
        end if
        do iwRank = 1, nval
            !%
            !%  We still have unprocessed values in both A and B
            !%
            if (iIndA <= LmtnA) then
                if (iIndB <= nval) then
                    if (xInput(jwRankT(iIndA)) > xInput(xRanked(iIndB))) then
                        iRanking = xRanked (iIndB)
                        iIndB = iIndB + 1
                    else
                        iRanking = jwRankT (iIndA)
                        iIndA = iIndA + 1
                    end if
                else
                !%
                !%  Only A still with unprocessed values
                !%
                    iRanking = jwRankT (iIndA)
                    iIndA = iIndA + 1
                end if
            else
                !%
                !%  Only B still with unprocessed values
                !%
                iRanking = xRanked (iwRank)
            end If
            if (xInput(iRanking) > xTst) then
                xTst = xInput (iRanking)
                Nunique = Nunique + 1
                xRanked (Nunique) = iRanking
            end if
        end do
            
        !%------------------------------------------------------------------
        !% Closing:
    end subroutine util_unique_rank        
!%
!%==========================================================================          
!%==========================================================================
!%
        !%------------------------------------------------------------------
        !% Description:
        !%
        !%------------------------------------------------------------------
        !% Declarations:
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------
    
    
        !%------------------------------------------------------------------
        !% Closing:
!%
!%==========================================================================
!% end OF MODULE
!%==========================================================================
end module utility
