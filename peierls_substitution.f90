function t_ham(num_bands, r, time, r_ham) result(ham_t)
use file_parsing, only : a_0, omega, phase_shift
implicit none
    integer, intent(in) :: num_bands
    real*8, intent(in) :: r(3), time
    complex*16, intent(in) :: r_ham(num_bands, num_bands)
    complex*16 :: ham_t(num_bands, num_bands)
    complex*16 :: element
    integer ::  ib, jb
    real*8 :: phase, vector_potential(3)

    vector_potential(1)=a_0*dcos(time*omega)
    vector_potential(2)=a_0*dsin(phase_shift+(time*omega))
    vector_potential(3)=0d0
    phase = dot_product(vector_potential, r)
    do ib=1, num_bands
        do jb=1, num_bands
            element = r_ham(ib, jb)
            element = element * dcmplx(dcos(phase), dsin(phase))
            ham_t(ib, jb) = element
        end do
    end do
return
end function t_ham
