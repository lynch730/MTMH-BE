
! Module for holding matlab boltzmann solution
module solver_td

    ! Load Modules
    use timers, only: timer_start_all, timer_stop_all
    use mdat, only: N, Neps, zin_N, N_gfrac, grid_EC, &
                    grid_E_INT_MASS, grid_E_INT_ENERGY, &
                    grid_L_is_0, Nrates, rates_zid, &
                    rates_integral, NL0, zin_field, zin_mass_int
    use mbsol, only: init_mbsol, reset_A_factorization, close_mbsol
    use lsolve, only: is_not_factorized
    use bz_solution, only: new_bz, reset_bz, close_bz_solution
    use solver_common
    
    ! Custom Types
    implicit none
    
    contains
    
    ! td solver param indices
    ! 1 -> Tgas (K)
    ! 2 -> Pressure (Pa)
    ! 3 -> omega (rad/s)
    ! 4 -> EN (townsend)
    ! 5 -> Texc (excitation temp)
    ! 6 -> Te0 -> Inital Electron Temp (eV)
    ! 7 -> tol_rel
    ! 8 -> tol_abs
    ! 9 -> tol_eps
    ! 10 -> jac_iter
    ! 11 -> max_jac_iter
    ! 12 -> logical, whether to use old jacobian
    ! 13 -> Max time
    ! 14 -> Integer Nt
    ! 13:20 -> Reserved!
    ! 21:21+Nfrac -> Base species gas fractions
    
    ! Converg EEDF 
    subroutine td_solution(param, Xin, Xk)
        
        real(8), intent(inout) :: param(:)
        real(8), intent(in), optional :: Xin(:), Xk(:)
        real(8) :: ebar, bz(zin_N), dt
        real(8) ::  X0(N), B(N), X02(N), Xtmp(N)
        logical :: is_converged
        integer(4) :: cnt, iter, iter_total, i
        integer(4) :: max_total_iter, max_jac_iter
        
        ! Extract counters
        max_jac_iter = int(param(10))
        max_total_iter = int(param(11)) 
        
        ! Initialize EEDF, if not given
        if (not(present(Xin))) then
            call initialize_eedf(param(6))
        else
            Xeedf(:) = Xin(:)
        end if
        
        ! Set Jacobian to be regenerated
        reset_bz = .true.
        bz(:) = 0.0d0
        
        ! Pre-store 
        Xall(:, 1) = Xeedf
        
        ! Main Time Loop
        do i = 2, Nt

            ! Start Timer
            call timer_start_all()
            
            ! Time Step
            dt = (2.0 / 3.0) * (time(i)-time(i-1))
            
            ! Field Strength function (Manual toggle)
            if (allocated(EN)) then
                param(4) = EN(i)
            end if
            
            ! if requested, regenerate factorized Jacobian
            if (abs(param(12)) > 3.0d0*epsilon(1.0d0) .or. is_not_factorized()) then
                if (present(Xk)) then
                    call new_bz(bz, Xeedf, param(1:5), param(21:20+N_gfrac), Xk)
                else
                    call new_bz(bz, Xeedf, param(1:5), param(21:20+N_gfrac))
                end if
                bz(zin_mass_int) = 0.0d0
                call reset_A_factorization(bz, dt)
                if (meedf_verbose>0) write(*,*) 'New Jacobian'
            end if
            
            ! Store previous
            X02 = Xall(:, i-1)
            if (i > 2) X02 = X02*(4.0d0/3.0d0) - Xall(:, i-2)/3.0d0 
            
            ! Main loop
            is_converged = .false.
            iter = 0
            iter_total = 0
            if (meedf_verbose>0) write(*,*) 'New Solution'
            do while (not(is_converged) .and. iter_total<max_total_iter) 
                
                ! Advance counters
                iter = iter + 1
                iter_total = iter_total + 1
                
                ! Reset X0
                X0 = Xeedf
                
                ! New Bz
                if (present(Xk)) then
                   call new_bz(bz, X0, param(1:5), param(21:20+N_gfrac), Xk)
                else
                   call new_bz(bz, X0, param(1:5), param(21:20+N_gfrac))
                end if
                bz(zin_mass_int) = 0
            
                ! call iteration_step(bz, dt, X0, Xeedf)
                call compute_dxdt(bz, X0, B)
                call timer_start()
                B = X02 + dt*B - X0 ! X02 + dt.*B - X;
                call solve_Ax_B(B, Xtmp)
                Xeedf = Xeedf + Xtmp
                call timer_stop_sol()
                
                ! Test Convergence
                is_converged = test_convergence(X0, Xeedf, param(7:9))
                
                ! Test negative mean energy
                ebar = energy_moment()
                
                ! Reset Jacboian
                if (iter >= max_jac_iter) then
                    iter = 0
                    if (present(Xk)) then
                       call new_bz(bz, X0, param(1:5), param(21:20+N_gfrac), Xk)
                    else
                       call new_bz(bz, X0, param(1:5), param(21:20+N_gfrac))
                    end if
                    bz(zin_mass_int) = 0
                    call reset_A_factorization(bz, dt)
                    if (meedf_verbose>0) write(*,*) 'New Jacobian'
                end if
             
                ! Print Status
                if (meedf_verbose>0) write(*,*) 'iter:',iter, 'Total iter:', iter_total, 'ebar=', ebar, is_converged
                !write(*,*) 'iter:',iter, 'Total iter:', iter_total, 'ebar=', ebar, is_converged
                 
            end do
            

            cnt = cnt + 1
            if (cnt>=1000) then
              cnt = 0
              write(*,*) time(i), param(4), ebar!, iter, iter_total
            end if

            ! Store final value of Xeedf
            Xall(:, i) = Xeedf
            
            ! Finish Timer 
            call timer_stop_all()
        
        end do
    
        ! Call moments
        call eedf_moments(bz)
        
        return
        
    end subroutine td_solution

    
end module solver_td
    
