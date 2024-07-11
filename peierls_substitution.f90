module peierls_substitution
use, intrinsic :: iso_fortran_env, only: real64 
implicit none
contains
    function t_ham(nbands, r, time, r_ham) result(ham_t)
    use file_parsing, only : a_0, omega, phase_shift
    implicit none
        integer, intent(in) :: nbands, r(3)
        real (real64), intent(in) :: time
        complex (real64), intent(in) :: r_ham(nbands, nbands)
        complex (real64) :: ham_t(nbands, nbands)
        complex (real64) :: element
        integer ::  ib, jb
        real (real64) :: phase, vector_potential(3)
    
        vector_potential(1)=a_0*cos(time*omega)
        vector_potential(2)=a_0*sin(phase_shift+(time*omega))
        vector_potential(3)=0d0
        phase = dot_product(vector_potential, r)
        do ib=1, nbands
            do jb=1, nbands
                element = r_ham(ib, jb)
                element = element*cmplx(cos(phase), sin(phase), kind=real64)
                ham_t(ib, jb) = element
            end do
        end do
    end function t_ham
end module peierls_substitution

