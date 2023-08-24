
    
!  timedependentpulse.f90 
program time_dependent_pulse

    use mdat,   only: N_gfrac, N
    use solver_td,  only: td_solution
    use solver_common, only: Xall, time, EN, Nt, Xeedf, &
                                mean_energy, close_solver, init_solver, &
                                initialize_eedf
    use timers, only: reset_timers, close_timers, set_timer_index, print_timers
    
    implicit none
    real(8) :: T, P, W, Texc, EN0, dt, param(500)
    real(8) :: rtol, atol, Te_0
    integer :: i, k, jac_iter, max_jac_iter, Nsim
    integer :: Neps_array(10)
    character(len=100) :: fname, dname
    real(8), allocatable :: ebar(:)
    
    !call omp_set_num_threads(1)
    !call mkl_set_num_threads(4)
    
    ! Parameters
    param(:) = 0.0d0
    T = 300.0d0
    P = 101325.0d0
    W = 1.539380400258999e+10
    Texc = T
    Te_0 = 2.585199978643554e-02
    rtol = 1.0d-8
    atol = 1.0d-15
    jac_iter = 5
    max_jac_iter = 500

    ! Collect parameters
    param(1:11) = (/ T, P, W, EN0, Texc, Te_0, rtol, atol, &
                     0.0d0, dble(jac_iter), dble(max_jac_iter) /)
    param(21) = 1.0 ! Mass fraction
    
    ! Get folder path from command line arguments
    if (iargc()>0) then
        call getarg(1, dname)
        write(*,*) 'Using Folder: ', dname
    else
        dname = '../data/loki_b_pulse/'
    end if
        
    ! Set timer
    Neps_array = (/ 100, 150, 200, 250, 300, 350, 400, 450, 500, 5000 /)
    !Neps_array = (/ 100, 100, 100, 100, 100, 100, 100, 100, 100, 500/)
    Nsim = 10
    call reset_timers(Nsim)
    
    ! Energy
    allocate(ebar(Nt))
    ebar(:) = 0.0d0
   
    !! ----------------------------------
    do k = 1,Nsim

        ! Initalize QSS Solver
        write(fname,'(I)') Neps_array(k)
        call init_solver( trim(adjustl(dname))//'loki_b_pulse_N'//&
                          trim(adjustl(fname))//'_NL2_NK8.mat', .false. )
        
        ! Initalize EEDF 
        call initialize_eedf(Te_0)
        
        ! Set Current Timer Index
        call set_timer_index(k)
        
        ! Timers
        Nt = 13000+1
        if (allocated(Xall)) deallocate(Xall)
        if (allocated(time)) deallocate(time)
        if (allocated(EN)) deallocate(EN)
        allocate(Xall(N,Nt))
        allocate(time(Nt), EN(Nt))
        dt = (log10(1.0e-05) - log10(1.0e-12))/dble(Nt-2)
        time(1) = 0.0e0
        time(2) = 1.0e-12
        do i = 3, Nt+1
            time(i) = 10.0d0**(log10(time(i-1))+dt)
        end do
        
        ! Set Field Profile
        EN = 50.0d0 * sqrt(2.0) * sqrt(time/1.0e-6) * dexp(-time/1.0e-6)
        
        ! Call td_solution solution for Param
        call td_solution(param, Xeedf)
        
        ! Store Energies
        ebar(1) = mean_energy
        
        ! Close Solver
        call close_solver()
        
        !write(*,*) k, Neps_array(k)

    end do

    ! Print And close Timers
    call print_timers()
    call close_timers()

end program time_dependent_pulse

