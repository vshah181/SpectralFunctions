subroutine fourier_coefficient(nbands, r_ham, m, t0, r, ham_m)
use constants, only: tau
use file_parsing, only: omega
implicit none
    complex*16, external :: t_ham
    integer, parameter :: n=500
    integer, intent(in) :: nbands, m
    real*8, intent(in) :: t0, r(3)
    complex*16, intent(in) :: r_ham(nbands, nbands)
    integer :: i
    real*8 :: period
    complex*16, intent(out) :: ham_m(nbands, nbands)
    
    period=tau/omega
! Integrate using the trapezium rule between t0 and t0+(2pi/omega)
    ham_m=0.5d0*(zexp(dcmplx(0d0, -m*omega*t0))*t_ham(nbands, r, t0, r_ham) &
               +zexp(dcmplx(0d0, -m*omega*(t0+period)))*t_ham(nbands, r,    &
                                                             t0+period, r_ham))
    do i=1, n-1
        ham_m=ham_m+(zexp(dcmplx(0, -m*omega*(t0+i*(period/n))))               &
        *t_ham(nbands, r, (t0+i*(period/n)), r_ham))
    end do
    ham_m=ham_m/n

end subroutine fourier_coefficient

function t_ham(nbands, r, time, r_ham) result(ham_t)
use file_parsing, only : a_0, omega, phase_shift
implicit none
    integer, intent(in) :: nbands
    real*8, intent(in) :: r(3), time
    complex*16, intent(in) :: r_ham(nbands, nbands)
    complex*16 :: ham_t(nbands, nbands)
    complex*16 :: element
    integer ::  ib, jb
    real*8 :: phase, vector_potential(3)

    vector_potential(1)=a_0*dcos(time*omega)
    vector_potential(2)=a_0*dsin(phase_shift+(time*omega))
    vector_potential(3)=0d0
    phase = dot_product(vector_potential, r)
    do ib=1, nbands
        do jb=1, nbands
            element = r_ham(ib, jb)
            element = element * dcmplx(dcos(phase), dsin(phase))
            ham_t(ib, jb) = element
        end do
    end do
end function t_ham
