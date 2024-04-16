subroutine make_ene_window(nene, emin, emax, de, omegas)
implicit none
    integer, intent(in) :: nene
    real*8, intent(in) :: emin, emax, de
    real*8, intent(out) :: omegas(nene)
    integer :: i
    
    omegas(1) = emin
    do i=2, nene
        omegas(i)=omegas(i-1)+de
    end do
end subroutine make_ene_window

