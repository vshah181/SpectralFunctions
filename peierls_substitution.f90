module peierls_substitution
use, intrinsic :: iso_fortran_env, only: real64 
use file_parsing, only: wannier_centres
implicit none
contains
    function t_ham(nbands, r, time, r_ham) result(ham_t)
    use file_parsing, only : a_0, omega, phase_shift, avec
        integer, intent(in) :: nbands, r(3)
        real (real64), intent(in) :: time
        complex (real64), intent(in) :: r_ham(nbands, nbands)
        complex (real64) :: ham_t(nbands, nbands)
        integer ::  i, j, k
        real (real64) :: phase, vector_potential(3), r_real(3)
    
        r_real=0d0
        vector_potential(1)=-1d0*a_0*sin(time*omega)
        vector_potential(2)=a_0*cos(phase_shift+(time*omega))
        vector_potential(3)=0d0
        ! Natural units: hbar = c = -e = 1
        ! ham_t=r_ham*cmplx(cos(phase), sin(phase), kind=real64)
        ham_t=r_ham
        do i=1, nbands
            do j=1, nbands
                do k=1, 3
                    r_real=r_real+(r(k)*avec(k, :))
                enddo
                r_real=r_real+(wannier_centres(j, :)-wannier_centres(i, :))
                phase=dot_product(vector_potential, r_real)
                ham_t(i, j)=ham_t(i, j)*cmplx(cos(phase), sin(phase),          &
                    kind=real64)
            enddo
        enddo
    end function t_ham
end module peierls_substitution

