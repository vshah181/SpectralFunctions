subroutine ft_ham_r(num_bands, k, k_ham, r_list, weights, r_ham_list,          &
    num_r_pts, nlayers)
use constants, only : tau
use file_parsing, only : id => direction, bulk
use, intrinsic :: iso_fortran_env, only: real64 
implicit none
    integer, intent(in) :: num_r_pts, num_bands, nlayers
    complex (real64), intent(out) :: k_ham(num_bands*nlayers,                  &
        num_bands*nlayers)
    complex (real64), intent(in) :: r_ham_list(num_r_pts, num_bands, num_bands)
    real (real64) :: k(3), phase, r(3)
    integer, intent(in) :: r_list(num_r_pts, 3), weights(num_r_pts)
    integer :: ir, il, irow, icol

    k_ham=0d0
    if(.not. bulk) then
        do ir=1, num_r_pts
            r=r_list(ir, :)
            phase=dot_product(k*tau, r)
            irow=int((num_bands*(abs(r(id))-r(id))/2)+1)
            !irow = 1 if +ve z, else irow=|z|*num_bands
            icol=int((num_bands*(abs(r(id))+r(id))/2)+1)
            !icol=1 if -ve z, else icol=|z|*num_bands
            do il=1, nlayers-int(abs(r(id)))
                k_ham(irow:irow+num_bands-1, icol:icol+num_bands-1)=           &
                k_ham(irow:irow+num_bands-1, icol:icol+num_bands-1)+           &
                (r_ham_list(ir, :, :)*cmplx(cos(phase),                        &
                 -sin(phase), kind=real64))/weights(ir)
                irow=irow+num_bands
                icol=icol+num_bands
            enddo
        enddo
    else
        do ir=1, num_r_pts
            phase=dot_product(k*tau, r_list(ir, :))
            k_ham = k_ham + (r_ham_list(ir, :, :) * cmplx(cos(phase),          &
                 -sin(phase), kind=real64) / weights(ir))
        enddo
    endif
end subroutine ft_ham_r
