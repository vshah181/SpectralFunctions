module file_parsing
use, intrinsic :: iso_fortran_env, only : iostat_end, real64
implicit none
private
    character(len=99) :: hr_file, seedname, basis, nnkp_file, wout_file
    character(len=22), parameter :: kpt_file="kpoints", input_file='INPUT'
    complex (real64), allocatable :: r_ham_list(:, :, :)
    real (real64), allocatable :: high_sym_pts(:, :), potential(:),            &
        projection_centres(:, :)
    real (real64) :: bvec(3, 3), avec(3, 3)
    real (real64) :: e_fermi, emin, emax, de, eta, phase_shift, omega, a_0
    integer, allocatable :: r_list(:, :), weights(:)
    integer :: num_bands, num_r_pts, nkpt_per_path, nkpath, nlayers, direction,&
        max_order
    logical :: e_fermi_present, do_floquet, bulk, bandplot, gfplot, use_soc
    public num_bands, num_r_pts, weights, r_list, r_ham_list, read_hr, eta, de,&
           write_spec_func, high_sym_pts, nkpath, nkpt_per_path, read_kpoints, &
           read_potential, nlayers, basis, bvec, avec, emin, emax, do_floquet, &
           a_0, potential, direction, read_vector_potential, omega, max_order, &
           phase_shift, bulk, gfplot, bandplot, write_energies,                &
           projection_centres 
contains
    subroutine read_input
        character(len=99) :: label, ival, line, temp_line
        integer :: i, eof, floquet_switch, bulk_switch, gf_switch, band_switch
        open(110, file=input_file, iostat=eof)
        e_fermi_present = .false.
        bulk=.false.
        use_soc = .false.
        do_floquet = .false.
        gfplot = .false.
        bandplot = .false.
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
                write(wout_file, fmt='(2a)') trim(adjustl(seedname)), '.wout'
            else if(trim(adjustl(label)) .eq. 'nlayers') then
                read(ival, *) nlayers
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
            else if(trim(adjustl(label)) .eq. 'bulk') then
                read(ival, *) bulk_switch
                if(bulk_switch .eq. 1) bulk=.true.
            else if(trim(adjustl(label)) .eq. 'bands_plot') then
                read(ival, *) band_switch
                if(band_switch .eq. 1) bandplot=.true.
            else if(trim(adjustl(label)) .eq. 'spectra_plot') then
                read(ival, *) gf_switch
                if(gf_switch .eq. 1) gfplot=.true.
            else if(trim(adjustl(label)) .eq. 'soc') then
                read(ival, *) use_soc
            endif
        enddo
        close(110)
        if(.not.(gfplot .or. bandplot)) then
            print*, "You want neither eigenvalues nor spectra? Doing nothing..."
            stop
        endif
        if(bulk) nlayers=1
        allocate(potential(nlayers))
    end subroutine read_input

    subroutine read_hr
        integer :: hi_row, hi_col, ir, o_i, o_j
        real (real64) :: rp, ip

        call read_input

        open(111, file=hr_file)
        read(111, *)
        read(111, *)num_bands, num_r_pts  ! bands and number of real-points
        allocate(weights(num_r_pts), r_list(num_r_pts, 3),                     &
                 r_ham_list(num_r_pts, num_bands, num_bands),                  &
                 projection_centres(num_bands, 3))
        read(111, *)weights ! degeneracy of each Wigner-Seitz grid point
        do ir=1, num_r_pts
            do o_i=1, num_bands
                do o_j=1, num_bands
                    read(111, *)r_list(ir, 1), r_list(ir, 2), r_list(ir, 3),   &
                              hi_row, hi_col, rp, ip
                    r_ham_list(ir, hi_row, hi_col)=cmplx(rp, ip, kind=real64)  
                enddo
            enddo
        enddo
        close(111)
        call read_nnkp
    end subroutine read_hr

    subroutine read_kpoints
        integer :: i
        open(112, file=kpt_file)
        read(112, *) nkpt_per_path
        read(112, *) nkpath
        allocate(high_sym_pts(nkpath, 3))
        do i=1, nkpath
            read(112, *) high_sym_pts(i, :)
        enddo
        close(112)
        nkpath=nkpath-1
    end subroutine read_kpoints

    subroutine read_nnkp
        character(len=99) :: line, temp_line
        integer :: eof, i, num_proj 

        open (113, file=nnkp_file, iostat=eof)

        do while(eof .ne. iostat_end)
            read(113, '(a)', iostat=eof) line
            if (trim(adjustl(line)) .eq. 'begin real_lattice') then
                do i=1, 3
                    read(113, *, iostat=eof) avec(i, :)
                enddo
            else if (trim(adjustl(line)) .eq. 'begin recip_lattice') then
                do i=1, 3
                    read(113, *, iostat=eof) bvec(i, :)
                enddo
            else if (trim(adjustl(line)) .eq. 'begin projections') then
                read(113, *) num_proj
                do i=1, num_proj
                    read(113, *) projection_centres(i, :)
                    read(113, *) temp_line
                enddo
                if ((num_proj*2 .eq. num_bands) .and. (basis.eq.'uudd')) then
                    do i=1, num_proj
                        projection_centres(i+num_proj, :)                      &
                            =projection_centres(i, :) 
                    enddo
                endif
            else if (trim(adjustl(line)) .eq. "begin spinor_projections") then
                read(113, *) num_proj
                do i=1, num_proj
                    read(113, *) projection_centres(i,:)
                    read(113, *) temp_line
                    read(113, *) temp_line
                enddo
            endif
        enddo
        close(113)
    end subroutine read_nnkp

    subroutine read_potential
        integer :: i
        open (115, file='potential.dat')
        do i=1, nlayers
            read(115, *) potential(i)
        enddo
        close(115)
    end subroutine read_potential

    subroutine read_vector_potential
        use constants, only : pi, reduced_planck_constant_ev
        integer :: eof, i
        real (real64) :: s
        character(len=99) :: label, ival, line, temp_line
        open(116, file='vector_potential.dat', iostat=eof)
        do while(eof .ne. iostat_end)
            read(116, '(a)', iostat=eof) line
            temp_line=adjustl(line)
            i=index(temp_line, ' ')
            label=temp_line(1:i)
            ival=temp_line(1+i:99)
            if(trim(adjustl(label)) .eq. 'phase_shift/pi') then
                read(ival, *) phase_shift
            else if(trim(adjustl(label)) .eq. 's') then
                read(ival, *) s
            else if(trim(adjustl(label)) .eq. 'hbar*omega') then
                read(ival, *) omega
            else if(trim(adjustl(label)) .eq. 'max_order') then
                read(ival, *) max_order
            endif
        enddo
        close(116)
        omega = omega / reduced_planck_constant_ev ! to get it in hertz
        phase_shift = phase_shift * pi ! to get it in radians
        a_0 = (s)!/norm2(avec(1, :))
        ! Ignore hbar, c and elementary charge.
    end subroutine read_vector_potential


    subroutine write_spec_func(spec_func, nkp, nene)
        integer, intent(in) :: nene, nkp
        real (real64), intent(in) :: spec_func(nene*nkp*nlayers)
        integer :: i
        character(len=99) :: ofname

        write(ofname, '(2a)') trim(adjustl(seedname)), '_spec_func.dat'
        print*, 'Writing output to: ', trim(adjustl(ofname))
        open (201, file=trim(adjustl(ofname)))
        do i=1, nene*nlayers*nkp
            write(201, *) spec_func(i)
        enddo
        close(201)
    end subroutine write_spec_func

    subroutine write_energies(energies, nkp, tot_bands)
        integer, intent(in) :: tot_bands, nkp
        real (real64), intent(in) :: energies(tot_bands*nkp)
        integer :: i
        character(len=99) :: ofname

        write(ofname, '(2a)') trim(adjustl(seedname)), '_eigenval.dat'
        print*, 'Writing output to: ', trim(adjustl(ofname))
        open(202, file=trim(adjustl(ofname)))
            do i=1, nkp*tot_bands
                write(202, *) energies(i)
            enddo
        close(202)
    end subroutine write_energies

end module file_parsing
