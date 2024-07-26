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
    integer :: i, j
    real (real64) :: period
    complex (real64), intent(out) :: ham_m(nbands, nbands)
    
    period=tau/omega
    ! Integrate using the trapezium rule between t0 and t0+(2pi/omega)
    ham_m=0.5d0*(exp(cmplx(0d0, -m*omega*t0, kind=real64))*t_ham(nbands, r, t0,&
                 r_ham)+exp(cmplx(0d0, -m*omega*(t0+period), kind=real64))     &
                 *t_ham(nbands, r, t0+period, r_ham))
    do i=1, n-1
        ham_m=ham_m+(exp(cmplx(0, -m*omega*(t0+i*(period/n)), kind=real64))    &
        *t_ham(nbands, r, (t0+i*(period/n)), r_ham))
    enddo
    ham_m=ham_m/n
end subroutine fourier_coefficient

