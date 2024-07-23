program spectral_function
use mpi
use file_parsing
use, intrinsic :: iso_fortran_env, only: real64, int32
implicit none
    real (real64), allocatable :: kp(:,:), kdists(:), hsym_kdists(:)
    integer :: nkp, ik, tot_bands, nene, nf_bands, nkpar, ip, ikpar, ikgfb,    &
        ikgfe, ikeneb, ikenee
    integer :: ibeg, iend, pid, ncpus, ierr, extra  ! for mpi
    integer, allocatable :: sendcounts(:), displs(:)! for mpi
    integer, allocatable :: sendcounts_ene(:), displs_ene(:)! for mpi
    integer (int32):: info, lwork
    real (real64), allocatable :: rwork(:), energies(:), omegas(:),         &
        energies_glob(:)
    complex (real64), allocatable :: work(:), kham(:, :), green_func(:),       &
        green_func_glob(:), floquet_ham_list(:, :, :)

    ! Initialise MPI here
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ncpus, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
    ! Only the root process (0) needs the displs and sendcounts arrays
    if (pid .eq. 0) allocate(displs(ncpus), sendcounts(ncpus),                 &
                    sendcounts_ene(ncpus), displs_ene(ncpus))

    call read_kpoints ! nkpath, high_sym_pts, nkpt_per_path
    call read_hr ! r_list, r_ham_list, weights, num_r_pts, num_bands, nlayers
                 ! do_floquet, max_order
    call read_potential ! potential
    max_order=0
    if(do_floquet) call read_vector_potential ! vector_potential

    nene=int((emax-emin)/de)
    tot_bands=num_bands*nlayers*(1+max_order)
    nf_bands=num_bands*(1+max_order)
    nkp=1+(nkpt_per_path*nkpath)

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
    nkpar=(iend-(ibeg-1))
    if(pid .eq. 0) write(*, '(i6,x,a)') nkpar,                                 &
        'kpoints to be calculated on node 0'

    allocate(kp(nkp, 3), kdists(nkp), hsym_kdists(1+nkpath), omegas(nene))
    call make_kpath(nkpath, high_sym_pts, nkpt_per_path, nkp, kp, kdists,      &
        hsym_kdists)
    call make_ene_window(nene, emin, de, omegas)

    ! Only the root process needs the recvbuf array
    if (pid .eq. 0 .and. gfplot) allocate(green_func_glob(nene*nkp*nlayers))
    if (pid .eq. 0 .and. bandplot) allocate(energies_glob(tot_bands*nkp))

    ! Create the arrays required by MPI_GATHERV
    if (pid .eq. 0) then
        displs(1)=0
        displs_ene(1)=0
        do ip=1, ncpus
            if (ip .gt. extra) then
                if(gfplot) sendcounts(ip)=(nkp/ncpus)*nene*nlayers
                if(bandplot) sendcounts_ene(ip)=(nkp/ncpus)*tot_bands
            else
                if(gfplot) sendcounts(ip)=(1+(nkp/ncpus))*nene*nlayers
                if(bandplot) sendcounts_ene(ip)=(1+(nkp/ncpus))*tot_bands
            end if
            if(ip .lt. ncpus) then
                displs(1+ip)=displs(ip)+sendcounts(ip)
                displs_ene(1+ip)=displs_ene(ip)+sendcounts_ene(ip)
            end if
        end do
    end if

    allocate(energies(tot_bands*nkpar), kham(tot_bands, tot_bands),            &
        green_func(nene*nkpar*nlayers), floquet_ham_list(num_r_pts, tot_bands, &
        tot_bands))
    call floquet_expansion(max_order, num_bands, num_r_pts, r_ham_list,        &
        floquet_ham_list, r_list)

    lwork=max(1, 2*tot_bands-1)
    allocate(work(max(1, lwork)), rwork(max(1, 3*tot_bands-2)))

    if(gfplot) green_func=cmplx(0d0, 0d0, kind=real64)
    do ik=ibeg, iend
        if(pid .eq. 0) print*, 'Calculating k-point', ik, 'of', nkpar, '...' 
        ikpar=ik-(ibeg-1)
        ikeneb=((ikpar-1)*tot_bands)+1
        ikenee=ikpar*tot_bands
        ikgfb=((ikpar-1)*nlayers*nene)+1
        ikgfe=ikpar*nene*nlayers
        ! ikpar always starts at 1
        call ft_ham_r(nf_bands, kp(ik, :), kham, r_list, weights,              &
            floquet_ham_list, num_r_pts, nlayers)
        ! call add_potential(kham, nlayers, num_bands)
        call zheev('V', 'L', tot_bands, kham, tot_bands,                       &
            energies(ikeneb:ikenee), work, lwork, rwork, info)
        if(gfplot) call greens_function(nene, nlayers, nf_bands, omegas, kham, &
            energies(ikeneb:ikenee), eta, green_func(ikgfb:ikgfe), num_bands)
        if(pid .eq. 0) print*, 'done'
    end do

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if(pid .eq. 0) print*, 'Gathering data...'
    if(gfplot) call MPI_GATHERV(green_func, nkpar*nlayers*nene,                &
        MPI_DOUBLE_COMPLEX, green_func_glob, sendcounts, displs,               &
        MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
    if(bandplot) call MPI_GATHERV(energies, nkpar*tot_bands,                   &
        MPI_DOUBLE_PRECISION, energies_glob, sendcounts_ene, displs_ene,       &
        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    if (pid .eq. 0 .and. gfplot) call write_spec_func(dimag(green_func_glob),  &
        nkp, nene)
    if (pid .eq. 0 .and. bandplot) call write_energies(energies_glob, nkp,     &
        tot_bands)
    call MPI_FINALIZE(ierr)
end program spectral_function
