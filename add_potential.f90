subroutine add_potential(k_ham, nlayers, num_bands)
use file_parsing, only : potential
implicit none
    integer, intent(in) :: nlayers, num_bands
    complex*16, intent(inout) :: k_ham(nlayers*num_bands, nlayers*num_bands)
    integer :: i, j, iblock
    real*8 :: layer_potential
    
    do i=1, nlayers
        layer_potential=potential(i)
        iblock=1+(i-1)*num_bands
        do j=iblock, iblock+(num_bands-1)
            k_ham(j, j)=k_ham(j, j)+layer_potential
        end do
    end do
end subroutine add_potential
