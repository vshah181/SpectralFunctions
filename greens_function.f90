subroutine greens_function(nene, nlayers, num_bands, omegas, eigkets, eigvals, &
    eta, green_func)
use, intrinsic :: iso_fortran_env, only: real64 
implicit none
    integer, intent(in) :: nlayers, num_bands, nene
    real (real64), intent(in) :: omegas(nene), eigvals(nlayers*num_bands), eta
    complex (real64), intent(in) :: eigkets(nlayers*num_bands,                 &
        nlayers*num_bands)
    complex (real64), intent(inout) :: green_func(nene*nlayers)
    integer :: ien, ib, il, layer_num, ig
    real (real64) :: wgt_ket(nlayers*num_bands), eigval, numerator
    complex (real64) :: eigket(nlayers*num_bands), denominator
    do ien=1, nene
        do ib=1, nlayers*num_bands
            eigket=eigkets(:, ib)
            wgt_ket=real(eigket*conjg(eigket))
            eigval=eigvals(ib)
            layer_num=1
            do il=1, nlayers*num_bands, num_bands
                ig=layer_num+(nlayers*(ien-1))
                numerator=sum(wgt_ket(il:il+(num_bands-1)))
                denominator=eigval-omegas(ien)-cmplx(0d0, eta, kind=real64)
                green_func(ig)=green_func(ig)+(numerator/denominator)
                layer_num=layer_num+1
            enddo
        enddo
    enddo
end subroutine greens_function
