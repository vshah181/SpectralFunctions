program spectral_function
use mpi
use file_parsing
use, intrinsic :: iso_fortran_env, only: real64, int32
implicit none
    integer, parameter :: max_order=2
    real (real64), allocatable :: kp(:,:), kdists(:), hsym_kdists(:)
    integer :: nkp, ik, tot_bands, nene, nf_bands
    integer :: ibeg, iend, pid, ncpus, ierr, extra  ! for mpi
    integer (int32):: info, lwork
    real (real64), allocatable :: rwork(:), energies(:, :), omegas(:)
    complex (real64), allocatable :: work(:), kham(:, :), green_func(:, :, :),&
        green_func_glob(:, :, :), floquet_ham_list(:, :, :)

    ! Initialise MPI here
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ncpus, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)

    call read_kpoints ! nkpath, high_sym_pts, nkpt_per_path
    call read_hr ! r_list, r_ham_list, weights, num_r_pts, num_bands, nlayers
    call read_potential ! potential
    call read_vector_potential ! vector_potential

    nene=int((emax-emin)/de)
    tot_bands=num_bands*nlayers*(1+(2*max_order))
    nf_bands=num_bands*(1+(2*max_order))
    nkp=1+(nkpt_per_path*nkpath)
    allocate(kp(nkp, 3), kdists(nkp), hsym_kdists(1+nkpath), omegas(nene))
    call make_kpath(nkpath, high_sym_pts, nkpt_per_path, nkp, kp, kdists,      &
        hsym_kdists)
    call make_ene_window(nene, emin, de, omegas)

    allocate(energies(nkp, tot_bands), kham(tot_bands, tot_bands),             &
        green_func(nkp, nlayers, nene), green_func_glob(nkp, nlayers, nene),   &
        floquet_ham_list(num_r_pts, tot_bands, tot_bands))
    call floquet_expansion(max_order, num_bands, num_r_pts, r_ham_list,        &
        floquet_ham_list, r_list)

    lwork=max(1, 2*tot_bands-1)
    allocate(work(max(1, lwork)), rwork(max(1, 3*tot_bands-2)))

    ! Split the kpoints among CPUs
    extra=mod(nkp, ncpus)

    ibeg=1+(pid*(nkp/ncpus))
    if (pid<extra) then
        ibeg=pid+ibeg
        iend=ibeg+(nkp/ncpus)
    else
        ibeg=extra+ibeg
        iend=ibeg+(nkp/ncpus)-1
    end if
    do ik=ibeg, iend
        call ft_ham_r(nf_bands, kp(ik, :), kham, r_list, weights,              &
            floquet_ham_list, num_r_pts, nlayers)
        ! call add_potential(kham, nlayers, num_bands)
        call zheev('V', 'L', tot_bands, kham, tot_bands, energies(ik, :), work,&
            lwork, rwork, info)
        call greens_function(nene, nlayers, nf_bands, omegas, kham,            &
            energies(ik, :), eta, green_func(ik, :, :))
    end do

    call MPI_REDUCE(green_func, green_func_glob, nkp*nlayers*nene,             &
        MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (pid .eq. 0) call write_spec_func(dimag(green_func_glob), nkp, nene)
    call MPI_FINALIZE(ierr)
end program spectral_function
