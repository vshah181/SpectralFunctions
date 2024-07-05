subroutine floquet_expansion(max_order, nbands, n_r_pts, r_ham_list,           &
    new_r_ham_list)
implicit none
    integer, intent(in) :: max_order, nbands, n_r_pts
    complex*16, intent(in) :: r_ham_list(n_r_pts, nbands, nbands)
    integer :: degree=(1+(2*max_order))*nbands
    integer :: m, irow, icol, i
    complex*16, intent(out) :: new_r_ham_list(n_r_pts, degree, degree)


    do m=-max_order, max_order
        irow=(num_bands*(abs(m)-m)/2)+1  ! 1 if +ve m, else |m|*num_bands
        icol=(num_bands*(abs(m)+m)/2)+1  ! 1 if -ve m, else |m|*num_bands
        do i=1, 1+(abs(i)-abs(order))
    end do
    
end subroutine floquet_expansion
