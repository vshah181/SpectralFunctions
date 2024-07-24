subroutine floquet_expansion(max_order, nbands, n_r_pts, r_ham_list,           &
    new_r_ham_list, r_list)
use, intrinsic :: iso_fortran_env, only: real64
use file_parsing, only: do_floquet, omega
use constants, only: reduced_planck_constant_ev
implicit none
    integer, intent(in) :: max_order, nbands, n_r_pts, r_list(n_r_pts, 3)
    complex (real64), intent(in) :: r_ham_list(n_r_pts, nbands, nbands)
    complex (real64), intent(out) :: new_r_ham_list(n_r_pts,                   &
        (1+max_order)*nbands, (1+max_order)*nbands)
    complex (real64) :: new_r_ham((1+max_order)*nbands, (1+max_order)*nbands), &
        ham_m(nbands, nbands), frequency_matrix((1+max_order)*nbands,          &
            (1+max_order)*nbands)
    real (real64) :: hw
    integer :: m, irow, icol, i, ir, im

    hw=omega*reduced_planck_constant_ev
    if(do_floquet) then
        frequency_matrix= cmplx(0d0, 0d0, kind=real64)
        im=max_order-1
        do i=1, nbands*(1+max_order), nbands
            frequency_matrix(i:i+(nbands-1), i:i+(nbands-1))=cmplx(im*hw, 0d0, &
                kind=real64)
            im=im-1
        enddo

        do ir=1, n_r_pts
            do m=-max_order, max_order
                irow=(nbands*(abs(m)-m)/2)+1  ! 1 if +ve m, else (|m|*nbands)+1
                icol=(nbands*(abs(m)+m)/2)+1  ! 1 if -ve m, else (|m|*nbands)+1
                call fourier_coefficient(nbands, r_ham_list(ir, :, :), m, 0d0, &
                    r_list(ir, :), ham_m)
                do i=1, (max_order-abs(m))+1
                    new_r_ham(irow:irow+nbands-1, icol:icol+nbands-1)=ham_m
                    icol=icol+nbands
                    irow=irow+nbands
                enddo
            enddo
            new_r_ham=new_r_ham+frequency_matrix
            new_r_ham_list(ir, :, :)=new_r_ham
        enddo
    else
        print*, 'Hello!'
        new_r_ham_list=r_ham_list
    endif
end subroutine floquet_expansion
