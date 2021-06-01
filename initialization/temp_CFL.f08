! This is just a temporary modeul for claculating the number of element in each link
!
! For now, we use a fixed element length to get the linkI(:,li_N_element)
! However, an offical CPL module is needed for fine calculate the exact N_element


module temp_CFL

    use globals
    use array_index
    
    implicit none 
    
    public 

    subroutine pupolate_N_of_element()
        integer :: ii

        do ii= 1,N_link
            linkI(ii,li_N_element) = ceiling(linkR(ii,lr_Length)/element_length)
        end do

    end subroutine pupolate_N_of_element
