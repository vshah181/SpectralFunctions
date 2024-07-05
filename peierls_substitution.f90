subroutine peierls_substitution(num_r_pts, r_list, num_bands, r_ham_list)
use file_parsing, only : vector_potential
implicit none
    integer, intent(in) :: num_r_pts, r_list(num_r_pts, 3), num_bands
    complex*16, intent(inout) :: r_ham_list(num_r_pts, num_bands, num_bands)
    complex*16 :: r_ham(num_bands, num_bands), element
    integer ::  ir, ib, jb
    real*8 :: r(3), phase

    do ir=1, num_r_pts
        r = r_list(ir, :)
        r_ham = r_ham_list(ir, :, :)
        phase = dot_product(vector_potential, r)
        do ib=1, num_bands
            do jb=1, num_bands
                element = r_ham(ib, jb)
                element = element * dcmplx(dcos(phase), dsin(phase))
                r_ham(ib, jb) = element
            end do
        end do
        r_ham_list(ir, :, :) = r_ham
    end do
end subroutine peierls_substitution
