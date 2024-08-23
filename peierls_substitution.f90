module peierls_substitution
use, intrinsic :: iso_fortran_env, only: real64 
implicit none
contains
    function t_ham(nbands, r, time, r_ham) result(ham_t)
    use file_parsing, only : a_0, omega, phase_shift, avec
        integer, intent(in) :: nbands, r(3)
        real (real64), intent(in) :: time
        complex (real64), intent(in) :: r_ham(nbands, nbands)
        complex (real64) :: ham_t(nbands, nbands)
        integer ::  i
        real (real64) :: phase, vector_potential(3), r_real(3)
    
        r_real=0d0
        vector_potential(1)=-1d0*a_0*sin(time*omega)
        vector_potential(2)=a_0*cos(phase_shift+(time*omega))
        vector_potential(3)=0d0
        do i=1, 3
            r_real=r_real+(r(i)*avec(i, :))
        enddo
        phase=dot_product(vector_potential, r_real)
        ! Igore hbar*c
        ham_t=r_ham*cmplx(cos(phase), sin(phase), kind=real64)
    end function t_ham
end module peierls_substitution

