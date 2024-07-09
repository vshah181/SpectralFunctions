subroutine floquet_expansion(max_order, nbands, n_r_pts, r_ham_list,           &
    new_r_ham_list, r_list)
use, intrinsic :: iso_fortran_env, only: real64
implicit none
    integer, intent(in) :: max_order, nbands, n_r_pts, r_list(n_r_pts, 3)
    integer :: degree
    complex (real64), intent(in) :: r_ham_list(n_r_pts, nbands, nbands)
    complex (real64), intent(out) :: new_r_ham_list(n_r_pts,                   &
        (1+(max_order*2))*nbands, (1+(max_order*2))*nbands)
    complex (real64) :: new_r_ham((1+(max_order*2))*nbands,                    &
        (1+(max_order*2))*nbands), ham_m(nbands, nbands)
    integer :: m, irow, icol, i, ir

    degree=(1+(max_order*2))*nbands
    do ir=1, n_r_pts
        do m=-max_order, max_order
            irow=(nbands*(abs(m)-m)/2)+1 !1 if +ve m, else (|m|*nbands)+1
            icol=(nbands*(abs(m)+m)/2)+1 !1 if -ve m, else (|m|*nbands)+1
            call fourier_coefficient(nbands, r_ham_list(ir, :, :), m, 0d0,     &
                r_list(ir, :), ham_m)
            print*, ham_m
            stop
            do i=1, (max_order-abs(m))+1
                new_r_ham(irow:irow+nbands-1, icol:icol+nbands-1)=ham_m
                icol=icol+nbands
                irow=irow+nbands
            end do
        end do
    new_r_ham_list(ir, :, :)=new_r_ham
    end do
end subroutine floquet_expansion
