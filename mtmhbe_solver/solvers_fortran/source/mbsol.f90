
! Module for holding matlab boltzmann solution
module mbsol

    ! Load Modules
    use timers, only: timer_start, timer_stop_rhs, timer_stop_sol, &
                      timer_stop_jac, timer_stop_lum
    use mdat, only: zin_N, N, NY, C_zw, C_pw, C_Vw, &
                    load_vars, destroy_vars, ind_bc_zero
    use lsolve, only: Amn_Vp, V_J, init_Amn, update_A, mult_Ax, &
                      factorize_A, close_A, ind_diag
    use mkl_spblas, only: matrix_descr, mkl_sparse_d_create_csr, &
                    sparse_index_base_one, sparse_memory_aggressive, &
                    sparse_operation_non_transpose, sparse_matrix_t, &
                    mkl_sparse_set_mv_hint, mkl_sparse_set_memory_hint, &
                    mkl_sparse_d_update_values, mkl_sparse_d_mv, &
                    mkl_sparse_destroy, sparse_matrix_type_general, &
                    mkl_sparse_optimize
    
    implicit none
    
    ! Public Variables
    integer(4) :: mbsol_verbose=0
    
    ! Private Variables
    type(sparse_matrix_t), private :: Y
    !real(8), allocatable, private :: bz_0(:)
    
    contains
    
    ! Initialize solution
    subroutine init_mbsol(fname, is_qss)        
        character(len=*), intent(in) :: fname
        logical, intent(in) :: is_qss
        
        ! Load variable from file
        call load_vars(fname)

        ! Initialize A and Pardiso
        call init_Amn(is_qss)

        ! Create CSR sparse matrix for Y matrix
        call build_Y(is_qss)

    end subroutine init_mbsol
    
    
    ! Perform 
    subroutine compute_dxdt(bz, x, xdot)

        ! Input/output vars
        real(8), intent(in) :: bz(zin_N), x(N)
        real(8), intent(out) :: xdot(N)
        
        ! Update Jacobian
        !call Ybz2A( -(bz - bz_0), Amn_Vp)
        call Ybz2A(bz, Amn_Vp)
        
        ! Create A Matrix using Iu/Ju with new b vector
        call update_A()
        
        ! Get RHS
        call timer_start()
        call mult_Ax(x, xdot)
        call timer_stop_rhs()
        
        return 
        
    end subroutine compute_dxdt
    
    
    ! Create Sparse Matrix from CSR
    subroutine build_Y(is_qss)
        logical, intent(in) :: is_qss
        integer(4) :: info
        type(matrix_descr) :: descrA
        
        ! Clear First Row Y values
        if (is_qss) then
            C_Vw(ind_bc_zero) = 0.0d0
        end if

        ! Create A Matrix using Iu/Ju with new b vector
        info = mkl_sparse_d_create_csr(Y, sparse_index_base_one, &
                            nY, zin_N, C_pw(1:), C_pw(2:), C_zw, C_Vw)
        
        descrA%type = sparse_matrix_type_general
        info = mkl_sparse_set_mv_hint (Y, sparse_operation_non_transpose, &
                                      descrA, 1000)
        info = mkl_sparse_set_memory_hint(Y, sparse_memory_aggressive)
        info = mkl_sparse_optimize(Y)
        
        return
        
    end subroutine build_Y
    
    
    ! Reset A factorization
    ! ---------------------------------
    subroutine reset_A_factorization(bz, dt)
        real(8), intent(in) :: bz(:)
        real(8), intent(in), optional :: dt
        
        ! Start Timer
        call Ybz2A(bz, V_J)
        if (mbsol_verbose>0) write(*,*) 'RESET_A_FACT: Ybz2A()'
        
        ! A = 1-dt*J
        if (present(dt)) then 
            V_J = -dt * V_J
            V_J(ind_diag) = V_J(ind_diag) + 1.0d0
        end if

        ! Factorize Matrix
        call timer_start()
        call factorize_A()
        call timer_stop_lum()
        if (mbsol_verbose>0) write(*,*) 'RESET_A_FACT: factorize_A()'

        return
               
    end subroutine reset_A_factorization
    
    
    ! Compute (Y*bz+Cee) and form into A matrix
    ! Need to call Ybz2A(bz, I_A, J_A, V_A, A) for internal matrix
    subroutine Ybz2A(bz, V)
        real(8), intent(in) :: bz(:)
        real(8), intent(out) :: V(:)
        type(matrix_descr) :: descrA
        integer(4) :: info
        
        ! Timer
        call timer_start()

        descrA%type = sparse_matrix_type_general
        V(:) = 0.0d0

        ! Solve Mv
        info = mkl_sparse_d_mv(sparse_operation_non_transpose, &
                             1.0d0, Y, descrA, bz, 0.0d0, V)
        
        ! Finish Timer
        call timer_stop_jac()
        
        return
        
    end subroutine Ybz2A
    
    ! Close out all data
    !------------------------------
    subroutine close_mbsol()
        integer(4) :: info
        
        ! Destroy Y
        info = mkl_sparse_destroy(Y)

        ! Destroy A
        call close_A()
        
        ! Clear out mdat data
        call destroy_vars()
        
        return
        
    end subroutine close_mbsol
    
    
end module mbsol
