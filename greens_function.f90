subroutine greens_function(nene, nlayers, num_bands, omegas, eigkets, eigvals, &
    eta, green_func)
    integer, intent(in) :: nlayers, num_bands, nene
    real*8, intent(in) :: omegas(nene), eigvals(nlayers*num_bands), eta
    complex*16, intent(in) :: eigkets(nlayers*num_bands, nlayers*num_bands)
    complex*16, intent(inout) :: green_func(nlayers, nene)
    integer :: ien, ib, il, layer_num
    real*8 :: wgt_ket(nlayers*num_bands), eigval, numerator
    complex*16 :: eigket(nlayers*num_bands), denominator
    do ien=1, nene
        do ib=1, nlayers*num_bands
            eigket=eigkets(:, ib)
            wgt_ket=eigket*dconjg(eigket)
            eigval=eigvals(ib)
            layer_num=1
            do il=1, nlayers*num_bands, num_bands
               numerator=sum(wgt_ket(il:il+(num_bands-1)))
               denominator=eigval-omegas(ien)-dcmplx(0.0, eta)
               green_func(layer_num, ien)=green_func(layer_num, ien)+numerator/denominator
               layer_num=layer_num+1
            end do
        end do
    end do
end subroutine greens_function
