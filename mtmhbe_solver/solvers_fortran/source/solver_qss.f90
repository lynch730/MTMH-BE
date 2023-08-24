
! Module for holding matlab boltzmann solution
module solver_qss

    ! Load Modules
    use timers, only: timer_start_all, timer_stop_all
    use mdat, only: N, Neps, zin_N, N_gfrac, grid_EC, &
                    grid_E_INT_MASS, grid_E_INT_ENERGY, &
                    grid_L_is_0, Nrates, rates_zid, &
                    rates_integral, NL0, zin_field
    use mbsol, only: init_mbsol, compute_dxdt, &
                     reset_A_factorization, close_mbsol
    use lsolve, only: is_not_factorized
    use bz_solution, only: new_bz, reset_bz, &
                             close_bz_solution
    use solver_common

    ! Custom Types
    implicit none
    
    contains
    
    
    ! Converge EEDF 
    subroutine qss_solution(param, Xin, Xk)
        real(8), intent(in) :: param(:)
        real(8), intent(in), optional :: Xin(:), Xk(:)
        real(8) :: ebar, X0(N), Xin2(N), B(N)
        real(8) :: bz(zin_N)

        logical :: is_converged
        integer(4) :: iter, iter_total
        integer(4) :: max_total_iter, max_jac_iter
        
        ! Start Timer
        call timer_start_all()
        
        ! Initialize EEDF, if not given
        if (not(present(Xin))) then
            call initialize_eedf(param(6))
        else
            Xeedf(:) = Xin(:)
        end if

        ! if requested, regenerate factorized Jacobian
        reset_bz = .true.
        bz(:) = 0.0d0
        if (abs(param(12)) > 3.0d0*epsilon(1.0d0) .or. is_not_factorized()) then
            if (present(Xk)) then
              call new_bz(bz, Xeedf, param(1:5), param(21:20+N_gfrac), Xk)
            else
              call new_bz(bz, Xeedf, param(1:5), param(21:20+N_gfrac))
            end if
            bz_0 = bz
            call reset_A_factorization(bz)
            if (meedf_verbose>0) write(*,*) 'New Jacobian'
        end if
        
        ! Extract counters
        max_jac_iter = int(param(10))
        max_total_iter = int(param(11)) 
        
        ! Main loop
        is_converged = .false.
        iter = 0
        iter_total = 0
        if (meedf_verbose>0) write(*,*) 'New Solution'
        do while (not(is_converged) .and. iter_total<max_total_iter) 
              
            ! Advance counters
            iter = iter + 1
            iter_total = iter_total + 1
            
            ! Perform steps
            X0 = Xeedf
            if (present(Xk)) then
               call new_bz(bz, X0, param(1:5), param(21:20+N_gfrac), Xk)
            else
               call new_bz(bz, X0, param(1:5), param(21:20+N_gfrac))
            end if

            ! Compute dxdt
            call compute_dxdt(-(bz - bz_0), X0, B)
            
            ! Correct bz for problem type
            call timer_start()
            B(1) = 1.0d0
            call solve_Ax_B(B, Xeedf)
            !if (Xold_frac > 0.0d0) Xeedf = (1.0d0-Xold_frac)*Xeedf + Xold_frac*X0;
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
                bz_0 = bz
                call reset_A_factorization(bz)
                if (meedf_verbose>0) write(*,*) 'New Jacobian'
            end if
             
            ! Print Status
            if (meedf_verbose>0) write(*,*) 'iter:',iter, 'Total iter:', iter_total, 'ebar=', ebar, is_converged
            
        end do
        
        ! Call moments
        call eedf_moments(bz)
        
        ! Finish Timer 
        call timer_stop_all()
        
        return
    end subroutine qss_solution
    
end module solver_qss
    
