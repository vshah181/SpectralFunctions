subroutine ft_ham_r(num_bands, k, k_ham, r_list, weights, r_ham_list, num_r_pts,&
    nlayers)
use constants, only : tau
use file_parsing, only : id => direction
implicit none
    complex*16, intent(out) :: k_ham(num_bands*nlayers, num_bands*nlayers)
    complex*16, intent(in) :: r_ham_list(num_r_pts, num_bands, num_bands)
    complex*16 :: r_ham(num_bands, num_bands)
    real*8 :: k(3), phase, r(3)
    integer, intent(in) :: r_list(num_r_pts, 3), num_r_pts, num_bands, nlayers,&
        weights(num_r_pts)
    integer :: ir, w, il, irow, icol

    k_ham=0d0
    do ir=1, num_r_pts
        r=r_list(ir, :)
        phase=dot_product(k*tau, r)
        irow=(num_bands*(abs(r(id))-r(id))/2)+1 !1 if +ve z, else |z|*num_bands
        icol=(num_bands*(abs(r(id))+r(id))/2)+1 !1 if -ve z, else |z|*num_bands
        do il=1, nlayers-int(abs(r(id)))
            k_ham(irow:irow+num_bands-1, icol:icol+num_bands-1)=               &
            k_ham(irow:irow+num_bands-1, icol:icol+num_bands-1)+               &
            (r_ham_list(ir, :, :)*dcmplx(cos(phase),                           &
            sin(phase)))/weights(ir)
            irow=irow+num_bands
            icol=icol+num_bands
        end do
    end do
end subroutine ft_ham_r
