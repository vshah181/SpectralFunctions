module file_parsing
use, intrinsic :: iso_fortran_env, only : iostat_end, real64
implicit none
private
    character(len=99) :: hr_file, seedname, basis, nnkp_file
    character(len=22), parameter :: kpt_file="kpoints", input_file='INPUT'
    complex (real64), allocatable :: r_ham_list(:, :, :)
    real (real64), allocatable :: high_sym_pts(:, :), potential(:)
    real (real64) :: bvec(3, 3)
    real (real64) :: e_fermi, emin, emax, de, eta, phase_shift, omega, a_0
    integer, allocatable :: r_list(:, :), weights(:)
    character, allocatable :: high_sym_pt_symbols(:)
    integer :: num_bands, num_r_pts, nkpt_per_path, nkpath, nlayers, direction,&
        max_order
    logical :: e_fermi_present, do_floquet
    public num_bands, num_r_pts, weights, r_list, r_ham_list, read_hr, eta, de,&
           write_spec_func, high_sym_pts, nkpath, nkpt_per_path, read_kpoints, &
           read_potential, nlayers, basis, bvec, emin, emax, do_floquet, a_0,  &
           potential, direction, read_vector_potential, omega, max_order, &
           phase_shift
contains
    subroutine read_input
        character(len=99) :: label, ival, line, temp_line
        integer :: i, eof, floquet_switch
        open(110, file=input_file)
        e_fermi_present = .false.
        do_floquet = .false.
        e_fermi = 0d0
        do while(eof .ne. iostat_end)
            read(110, '(a)', iostat=eof) line
            temp_line=adjustl(line)
            i=index(temp_line, ' ')
            label=temp_line(1:i)
            ival=temp_line(1+i:99)
            if(trim(adjustl(label)) .eq. 'seedname') then
                seedname=trim(adjustl(ival))
                write(hr_file, fmt='(2a)') trim(adjustl(seedname)), '_hr.dat'
                write(nnkp_file, fmt='(2a)') trim(adjustl(seedname)), '.nnkp'
            else if(trim(adjustl(label)) .eq. 'nlayers') then
                read(ival, *) nlayers
                allocate(potential(nlayers))
            else if(trim(adjustl(label)) .eq. 'basis') then
                read(ival, *) basis
            else if(trim(adjustl(label)) .eq. 'energy_range') then
                read(ival, *) emin, emax
            else if(trim(adjustl(label)) .eq. 'energy_step') then
                read(ival, *) de
            else if(trim(adjustl(label)) .eq. 'broadening_factor') then
                read(ival, *) eta
            else if(trim(adjustl(label)) .eq. 'direction') then
                read(ival, *) direction 
            else if(trim(adjustl(label)) .eq. 'e_fermi') then
                read(ival, *) e_fermi
                e_fermi_present = .true.
            else if(trim(adjustl(label)) .eq. 'floquet') then
                read(ival, *) floquet_switch
                if(floquet_switch .eq. 1) do_floquet=.true.
            end if
        end do
        close(110)
    end subroutine read_input

    subroutine read_hr
        integer :: hi_row, hi_col, ir, o_i, o_j
        real (real64) :: rp, ip

        call read_input
        call read_nnkp

        open(111, file=hr_file)
        read(111, *)
        read(111, *)num_bands, num_r_pts  ! bands and number of real-points
        allocate(weights(num_r_pts), r_list(num_r_pts, 3),                     &
                 r_ham_list(num_r_pts, num_bands, num_bands))
        read(111, *)weights ! degeneracy of each Wigner-Seitz grid point
        do ir=1, num_r_pts
            do o_i=1, num_bands
                do o_j=1, num_bands
                    read(111, *)r_list(ir, 1), r_list(ir, 2), r_list(ir, 3),   &
                              hi_row, hi_col, rp, ip
                    r_ham_list(ir, hi_row, hi_col)=cmplx(rp, ip, kind=real64)  
                end do
            end do
        end do
        close(111)
    end subroutine read_hr

    subroutine read_kpoints
        integer :: i
        open(112, file=kpt_file)
        read(112, *) nkpt_per_path
        read(112, *) nkpath
        allocate(high_sym_pts(nkpath, 3), high_sym_pt_symbols(nkpath))
        do i=1, nkpath
            read(112, *) high_sym_pts(i, :), high_sym_pt_symbols(i)
        end do
        close(112)
        nkpath=nkpath-1
    end subroutine read_kpoints

    subroutine read_nnkp
        character(len=99) :: line
        integer :: eof, i

        open (113, file=nnkp_file)

        do while(eof .ne. iostat_end)
            read(113, '(a)', iostat=eof) line
            if (trim(adjustl(line)) .eq. 'begin recip_lattice') then
                do i=1, 3
                    read(113, *, iostat=eof) bvec(i, :)
                end do
            end if
        end do
        close(113)
    end subroutine read_nnkp

    subroutine read_potential
        integer :: i
        open (114, file='potential.dat')
        do i=1, nlayers
            read(114, *) potential(i)
        end do
        close(114)
    end subroutine read_potential

    subroutine read_vector_potential
        use constants, only : pi, reduced_planck_constant_ev
        integer :: eof, i
        character(len=99) :: label, ival, line, temp_line
        open(115, file='vector_potential.dat')
        do while(eof .ne. iostat_end)
            read(115, '(a)', iostat=eof) line
            temp_line=adjustl(line)
            i=index(temp_line, ' ')
            label=temp_line(1:i)
            ival=temp_line(1+i:99)
            if(trim(adjustl(label)) .eq. 'phase_shift/pi') then
                read(ival, *) phase_shift
            else if(trim(adjustl(label)) .eq. 'A_0') then
                read(ival, *) a_0
            else if(trim(adjustl(label)) .eq. 'hbar*omega') then
                read(ival, *) omega
            else if(trim(adjustl(label)) .eq. 'max_order') then
                read(ival, *) max_order
            end if
        end do
        close(115)
        omega = omega / reduced_planck_constant_ev ! to get it in hertz
        phase_shift = phase_shift * pi ! to get it in radians
    end subroutine read_vector_potential


    subroutine write_spec_func(spec_func, nkp, nene)
        integer, intent(in) :: nene, nkp
        real (real64), intent(in) :: spec_func(nkp, nlayers, nene)
        integer :: ie, il
        real (real64) :: tot_spec_func(nkp, nene)
        character(len=99) :: ofname, fmt_string, nkp_string

        tot_spec_func=sum(spec_func, dim=2)
        write(ofname, '(2a)') trim(adjustl(seedname)), '_spec_func.dat'
        write(nkp_string, '(i10)') nkp
        write(fmt_string, '(3a)') '(', trim(adjustl(nkp_string)), '(1x,es19.9))'
        open (201, file=trim(adjustl(ofname)))
        do il=1, nlayers
            write(201, fmt='(a,i6)') 'layer= ', il
            do ie=1, nene
                write(201, fmt_string) spec_func(:, il, ie)
            end do
            write(201, fmt='(a)') ' '
        end do
        write(201, fmt='(a)') 'layer= all'
        do ie=1, nene
            write(201, fmt_string) tot_spec_func(:, ie)
        end do
        close(201)
    end subroutine write_spec_func

end module file_parsing
