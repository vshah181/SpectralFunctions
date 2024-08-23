subroutine fourier_coefficient(nbands, r_ham, m, t0, r, ham_m)
use constants, only: tau
use file_parsing, only: omega
use peierls_substitution, only: t_ham
use, intrinsic :: iso_fortran_env, only: real64 
implicit none
    integer, parameter :: n=500
    integer, intent(in) :: nbands, m, r(3)
    real (real64), intent(in) :: t0
    complex (real64), intent(in) :: r_ham(nbands, nbands)
    integer :: i
    real (real64) :: period
    complex (real64), intent(out) :: ham_m(nbands, nbands)
    
    ham_m=0d0
    period=tau/omega
    ! Integrate using rectangles
    do i=0, n-1
        ham_m=ham_m+t_ham(nbands, r, t0+i*(period/n), r_ham)*exp(cmplx(0d0,    &
            m*omega*(t0+i*(period/n)), kind=real64))
    end do
    ham_m=ham_m*(1d0/n)
end subroutine fourier_coefficient

