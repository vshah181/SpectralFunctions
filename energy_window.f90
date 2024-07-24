subroutine make_ene_window(nene, emin, de, omegas)
use, intrinsic :: iso_fortran_env, only: real64 
implicit none
    integer, intent(in) :: nene
    real (real64), intent(in) :: emin, de
    real (real64), intent(out) :: omegas(nene)
    integer :: i
    
    omegas(1) = emin
    do i=2, nene
        omegas(i)=omegas(i-1)+de
    enddo
end subroutine make_ene_window

