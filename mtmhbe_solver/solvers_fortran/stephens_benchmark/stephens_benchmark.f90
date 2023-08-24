
    
!  timedependentpulse.f90 
program stephens_benchmark

    use mdat,   only: N_gfrac, N
    use solver_qss, only: qss_solution
    use solver_td,  only: td_solution
    use solver_common, only: Xall, time, EN, Nt, Xeedf, &
                                mean_energy, close_solver, init_solver, &
                                initialize_eedf
    use timers, only: reset_timers, close_timers, set_timer_index, print_timers

    implicit none
    real(8) :: T, P, W, Texc, P_case(2), tmax(2), param(500)
    real(8) :: EN_array(3)
    real(8) :: rtol, atol, Te_0
    integer :: i, j, k, jac_iter, max_jac_iter, Nsim
    integer :: Ncase, N_EN, Nt_pw, Nwaves
    character(len=100) :: fname, dname, case_str
    real(8), allocatable :: ebar(:), Xstore(:)
    
    !call omp_set_num_threads(1)
    !call mkl_set_num_threads(4)
    
    ! Parameters
    param(:) = 0.0d0
    T = 0.0d0
    Texc = T    
    W = 6.911503837897545e+11
    Te_0 = 0.0001d0 
    rtol = 1.0d-8
    atol = 1.0d-15
    jac_iter = 20
    max_jac_iter = 500
    
    ! Reset parameters
    param(1:11) = (/ T, 101325.0d0, W, 0.0d0, Texc, Te_0, rtol, atol, &
                        0.0d0, dble(jac_iter), dble(max_jac_iter) /)
    param(21:22) = (/ 0.22d0, 0.78d0 /) ! Mass fraction
            
    ! Get folder path from command line arguments
    dname = '../data/stephens_benchmark/'
    
    ! Cases
    Ncase = 2
    P_case  = (/ 101325.0d0, 13332.2d0 /)
    tmax = (/ 10.0d-12, 25.0d-12 /)
    Nt_pw = 350 !256 ! 256
    
    ! Fields
    N_EN = 3
    EN_array = (/ 400.0d0, 700.0d0, 1000.0d0 /)
    EN_array = EN_array * sqrt(2.0d0)
    
    ! Set timer
    Nsim = 0
    call reset_timers(Ncase * N_EN * 2)

    ! Main Loop
    do k = 1, Ncase
        
        ! Set Pressure
        param(2) =  P_case(k)
        
        ! Case String
        write(case_str,'(I)') k
        fname = trim(adjustl(dname))//'stephens_P'//&
                trim(adjustl(case_str))//'.mat'
        write(*,*) 'FILE: ', trim(adjustl(fname))
        
        ! Reduced Field Sweep
        do j = 1, N_EN
            
            
            
            call init_solver(trim(adjustl(fname)), .true.)
            call initialize_eedf(Te_0)
            Nsim = Nsim + 1
            call set_timer_index(Nsim)
            param(4) = 0.01d0
            call qss_solution(param, Xeedf)
            if (allocated(Xstore)) deallocate(Xstore)
            allocate(Xstore(N))
            Xstore = Xeedf
            call close_solver()
            
            
            
            call init_solver(trim(adjustl(fname)), .false.)
            
            ! Set Current Timer Index
            Nsim = Nsim + 1
            call set_timer_index(Nsim)
            
            ! Initalize EEDF 
            !call initialize_eedf(Te_0)
            Xeedf = Xstore
            
            ! Time Dependent
            
            ! Set Time Steps
            Nwaves = nint(tmax(k) * W / 6.283185307179586d0)
            Nt = Nwaves * Nt_pw ! Nepsmber of steps to recover
            
            ! Reallocate Xall
            if (allocated(Xall)) deallocate(Xall)
            allocate(Xall(N,Nt))
            if (allocated(ebar)) deallocate(ebar)
            allocate(ebar(Nt))
            if (allocated(time)) deallocate(time)
            allocate(time(Nt))
            
            ! Zero Arrays
            Xall = 0.0d0
            ebar = 0.0d0
            time = 0.0d0
            
            ! Timers
            do i = 2, Nt
                time(i) = time(i-1) + (tmax(k)/dble(Nt-1))
            end do
            
            ! Set TD Field
            param(4) = EN_array(j)
            
            ! Call td_solution solution for Param
            call td_solution(param, Xeedf)
            
            ! Store Energies
            ebar(1) = mean_energy
            
            call close_solver()
            
        end do
    end do
    
    ! Print And close Timers
    call print_timers()
    call close_timers()
    if (allocated(Xstore)) deallocate(Xstore)
    

end program stephens_benchmark

