  
subroutine N2_bfactor(i, Texc, new_frac)
    use mdat, only: const_KB_QE
    integer, intent(in) :: i
    real(8), intent(in) :: Texc
    real(8), intent(out) :: new_frac(i)
    real(8) :: zsum
    real(8) :: de(59) = (/ 0.000000d0, 2.886000d-1, 5.737000d-1, 8.553000d-1, &
                           1.133500d0, 1.408200d0, 1.679400d0, 1.947000d0, &
                           2.211100d0, 2.471700d0, 2.728700d0, 2.982100d0, &
                           3.232000d0, 3.478200d0, 3.720800d0, 3.720800d0, &
                           4.195100d0, 4.426800d0, 4.654700d0, 4.878900d0, & 
                           5.099300d0, 5.315900d0, 5.528600d0, 5.737500d0, &
                           5.942300d0, 6.143200d0, 6.340000d0, 6.532700d0, &
                           6.721100d0, 6.905200d0, 7.085000d0, 7.260200d0, &
                           7.430900d0, 7.596800d0, 7.757900d0, 7.914000d0, &
                           8.065000d0, 8.210800d0, 8.351100d0, 8.485700d0, &
                           8.614600d0, 8.737400d0, 8.853900d0, 8.964000d0, &
                           9.067500d0, 9.163900d0, 9.253300d0, 9.335300d0, &
                           9.409800d0, 9.476700d0, 9.535900d0, 9.587500d0, &
                           9.631400d0, 9.667700d0, 9.696500d0, 9.718100d0, &
                           9.733300d0, 9.743300d0, 9.750100d0 /)
    
    new_frac = dexp(-de(1:i)/(Texc * const_KB_QE))
    zsum = sum(new_frac(:))
    new_frac = new_frac / zsum

end subroutine N2_bfactor
    
    
program laporta_grid_study

    use mdat,   only: N_gfrac
    use solver_qss,  only: initialize_eedf, qss_solution, Xeedf, &
                      init_solver, close_solver, mean_energy
    use timers, only: reset_timers, close_timers, set_timer_index, print_timers

    implicit none
    real(8) :: T, P, W, Texc, EN, param(500)
    real(8) :: rtol, atol, Te_0, Texc0, Texc1, dTexc
    integer :: j, k, jac_iter, max_jac_iter, NTexc
    integer :: grid_array(21), Nsim
    character(len=100) :: fname, dname
    real(8), allocatable :: ebar(:,:)
    
    !call omp_set_num_threads(1)
    !call mkl_set_num_threads(4)
    !call mkl_set_num_threads_local(4)
    !call mkl_domain_set_num_threads(4)
    
    ! Parameters
    param(:) = 0.0d0
    T = 300.0d0
    P = 101325.0d0
    W = 1.0d11
    EN = 10.0d0 * sqrt(2.0d0)
    Texc = T
    Te_0 = 0.5d0
    rtol = 1.0d-8
    atol = 1.0d-20
    jac_iter = 5
    max_jac_iter = 200

    ! Collect parameters
    param(1:11) = (/ T, P, W, EN, Texc, Te_0, rtol, atol, &
                     0.0d0, dble(jac_iter), dble(max_jac_iter) /)
    
    ! Texc array
    Texc0 = 300.0d0
    Texc1 = 3000.0d0
    dTexc = 50.0d0
    NTexc = int((Texc1-Texc0)/dTexc) + 1
    
    grid_array = (/ 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, &
                   700, 750, 800, 850, 900, 950, 1000, 4000 /)
    Nsim = size(grid_array, 1)

    ! Get folder path from command line arguments
    dname = '../data/grid_study_laporta_files/'
      
    ! Set timer
    call reset_timers(Nsim)
    
    ! Energy
    allocate(ebar(NTexc, Nsim))
    ebar(:,:) = 0.0d0
    
    !! ----------------------------------
    do k = 1,Nsim
        
        ! Initalize QSS Solver
        write(fname,'(I)') grid_array(k)
        call init_solver( trim(adjustl(dname))//'laporta_data_'//trim(adjustl(fname))//'.mat', .true. )
        
        ! Initalize EEDF 
        call initialize_eedf(Te_0)
        
        ! Set Current Timer Index
        call set_timer_index(k)
    
        ! Loop NTexc Cases
        do j = 1,NTexc 

            !write(*,*) j, '--------------------------------------'
            ! Set Texc 
            param(5) = Texc0 + dble(j-1)*dTexc ! Texc
            param(1) = param(5) ! T
            param(4) = EN * param(5) / Texc0
            call N2_bfactor(N_gfrac, param(5), param(21:20 + N_gfrac))
            
            ! Call qss_solution solution for Param
            if (j==1) then
              call qss_solution(param)
            else 
              call qss_solution(param, Xeedf)
            end if
            
            ! Store Energies
            ebar(j,k) = mean_energy
                        
        enddo
        
        ! Close Solver
        call close_solver()
        
        ! Print Status
        write(*,*) 'solution: ', real(k)/real(Nsim)
        
    end do
    
    ! Print Mean Energies
    write(*,*) 'Mean Energies'
    do k = 1, Nsim
        write(*,"(I','1X)", advance='no') k
        do j = 1, NTexc
            write(*,"(F9.6','1X)", advance='no') ebar(j,k)
        end do
        write(*,*)
    end do
    
    ! Print And close Timers
    call print_timers()
    call close_timers()
    
end program laporta_grid_study

