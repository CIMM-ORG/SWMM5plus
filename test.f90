program test
    integer :: ii(1)
    ii(1) = this_image()
    sync all
    print *, "My image id is: ", ii
    call co_min(ii)
    print *, "The smallest id is: ", ii
end program test