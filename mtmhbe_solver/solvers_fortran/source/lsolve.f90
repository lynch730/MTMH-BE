include 'mkl_pardiso.f90'

! Module for holding matlab boltzmann solution

! Summary:
!  init_lsove:      set N, init parameters
!  build_sparse_Amn:   
!  pointers_A:        c pointers to A matrix

module lsolve
    
    ! Load Modules
    use mdat, only: N, A_mp, A_np
    use mkl_pardiso, only: mkl_pardiso_handle
    use mkl_spblas, only: matrix_descr, mkl_sparse_d_create_csr, &
                    sparse_index_base_one, sparse_memory_aggressive, &
                    sparse_operation_non_transpose, sparse_matrix_t, &
                    mkl_sparse_set_mv_hint, mkl_sparse_set_memory_hint, &
                    mkl_sparse_d_update_values, mkl_sparse_d_mv, &
                    mkl_sparse_destroy, sparse_matrix_type_general, &
                    mkl_sparse_optimize, sparse_layout_column_major, &
                    mkl_sparse_convert_bsr
    implicit none
    
    ! Public Variables
    integer(4) :: lsolve_verbose = 0
    real(8), allocatable :: Amn_Vp(:), V_J(:)

    ! Private Variables
    type(sparse_matrix_t), private :: Amn
    type(mkl_pardiso_handle), private :: pt(64)
    integer(4), parameter, private :: nrhs=1, mtype=11, mnum=1, maxfct=1
    integer(4), allocatable, private :: perm(:)
    integer(4), private :: iparm(64)
    integer(4), allocatable :: ind_diag(:)
    logical, private :: not_factorized
    
    contains

    ! Return Fortran Pointers for Amn sparse matrix
    subroutine init_Amn(is_qss)
        logical, intent(in) :: is_qss
        
        integer(4) :: info, nY
        type(matrix_descr) :: descrA
        integer :: i, j
        
        not_factorized = .true.
        
        if (allocated(Amn_Vp)) deallocate(Amn_Vp)
        if (allocated(V_J)) deallocate(V_J)

        nY = size(A_np, 1)
        allocate(Amn_Vp(nY), V_J(nY))
        Amn_Vp(:) = 1.0d0
        V_J(:) = 1.0d0
        
        ! Get Index of diagonals
        if (allocated(ind_diag)) deallocate(ind_diag)
        allocate(ind_diag(N))
        ind_diag = -1
        do i = 1,N
            j = A_mp(i)
            do while (j < A_mp(i+1))
                if (i==A_np(j)) then
                    ind_diag(i) = j
                    j = A_mp(i+1)
                end if
                j = j + 1
            end do
            if (ind_diag(i)==-1) then
                write(*,*) 'ERROR, No diagonal index for i=', i
                stop
            end if
        end do
        
        ! Create Amn Matrix using Iu/Ju with new b vector
        info = mkl_sparse_d_create_csr(Amn, sparse_index_base_one, &
                                 N, N, A_mp(1:), A_mp(2:), A_np, Amn_Vp)
        
        info = mkl_sparse_convert_bsr(Amn, 1, sparse_layout_column_major, &
                                      sparse_operation_non_transpose, Amn)
        descrA%type = sparse_matrix_type_general
        info = mkl_sparse_set_mv_hint (Amn, sparse_operation_non_transpose, &
                                      descrA, 100)
        info = mkl_sparse_set_memory_hint(Amn, sparse_memory_aggressive)
        info = mkl_sparse_optimize(Amn)
        
        
        
        return
        
    end subroutine init_Amn
    
    
    ! Initialize Amn Matrix
    subroutine sym_factor_A()        
        integer(4), parameter :: phase = 11 ! factorization
        integer(4) :: i, error=0
        real(8) :: rdum
        
        ! Allocate permutation vectors
        if (allocated(perm))  deallocate(perm)
        allocate(perm(N))
        perm(:) = 0

        ! Initialize the internal solver memory pointer. This is only
        ! necessary for the FIRST call of the PARDISO solver.
        do i = 1, 64
           pt(i)%DUMMY =  0 
        end do
        !call pardisoinit(pt, mtype, iparm)
                
        ! Set up PARDISO control parameter
        iparm(:) = 0 
        iparm(1) = 1 ! no solver default
        iparm(2) = 2 ! fill-in reordering from METIS
        iparm(4) = 0 ! no iterative-direct algorithm
        iparm(5) = 0 ! no user fill-in reducing permutation
        iparm(6) = 0 ! =0 solution on the first n components of x
        iparm(8) = 200 ! max numbers of iterative refinement steps, alt 10
        iparm(10) = 13 ! perturb the pivot elements with 1E-13
        iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
        iparm(13) = 1 ! maximum weighted matching algorithm is on, alt 0
        iparm(14) = 0 ! Output: number of perturbed pivots
        iparm(18) = -1 ! Output: number of nonzeros in the factor LU
        iparm(19) = 0 ! Output: Mflops for LU factorization
        iparm(20) = 0 ! Output: Nepsmbers of CG Iterations
        iparm(24) = 0 ! 10 turns on advanced MKL solving routine, alt are 0 and 1
        iparm(25) = 0 ! 0, 1, 2 alternative parallelization
        iparm(27) = min(max(lsolve_verbose, 1), 0) ! checks input matrix is correct
    
        ! Call Reordering and Symbolic Factorization, This step also allocates
        ! all memory that is necessary for the factorization
        call pardiso(pt, maxfct, mnum, mtype, phase, N, V_J, A_mp, A_np, &
                     perm, nrhs, iparm, lsolve_verbose, rdum, rdum, error)
        
        ! Optional Output
        if (lsolve_verbose > 0) then
            write(*,*) 'Initialization complete'
            if (error /= 0) then
               write(*,*) 'The following error was detected: ', error
               stop
            end if
        end if
        
        
    end subroutine sym_factor_A
    
    
    subroutine update_A()
        integer :: info
        integer, pointer :: cnull(:) => null()
        info = mkl_sparse_d_update_values(Amn, 0, cnull, cnull, Amn_Vp)
        return 
    end subroutine 
    
    subroutine mult_Ax(x, b)
        real(8), intent(in) :: x(:)
        real(8), intent(out) :: b(N)
        integer(4) :: info
        type(matrix_descr) :: descrA
        descrA%type = sparse_matrix_type_general
        b(:) = 0.0d0
        info = mkl_sparse_d_mv(sparse_operation_non_transpose, &
                             1.0d0, Amn, descrA, x, 0.0d0, b)
        return
    end subroutine 
    
    ! Nepsmerical Factorization
    !------------------------------
    subroutine factorize_A()
    
        integer(4), parameter :: phase = 22 ! factorization
        integer(4) :: error = 0
        real(8) :: rdum = 0.0d0

        if (not_factorized) then
            call sym_factor_A()
            not_factorized = .false.
        end if
        
        ! Run factorization
        call pardiso(pt, maxfct, mnum, mtype, phase, N, V_J, A_mp, A_np, &
                      perm, nrhs, iparm, lsolve_verbose, rdum, rdum, error)
        
        ! Optional Output
        if (lsolve_verbose>0) then
            write(*,*) 'PARDISO: Factorization complete'
            if (error /= 0) then
                write(*,*) 'PARDISO: | The following error was detected: ', error
                stop
            endif
            write(*,*) 'PARDISO: | Factorization MFLOPS = ', iparm(19)
        endif
        
        return
        
    end subroutine factorize_A
    
    
    ! Solve RHS
    subroutine solve_Ax_B(b, x)
    
        ! Declare intent in
        real(8), intent(in) :: b(N)
        real(8), intent(inout) :: x(N)
        
        ! Private Variables
        integer(4), parameter :: phase = 33 ! pardiso solution step
        integer(4) :: error = 0
        
        ! Set initial values
        x(:) = 0.0d0
        
        ! Back substitution and iterative refinement
        call pardiso(pt, maxfct, mnum, mtype, phase, N, V_J, A_mp, A_np, &
                     perm, nrhs, iparm, lsolve_verbose, b, x, error)
    
        ! Optional Output
        if (lsolve_verbose>0) then
            write(*,*) 'PARDISO: Solution complete, iter:  ', iparm(7)
            if (error /= 0) then
                write(*,*) 'PARDISO: | The following error was detected: ', error
                stop
            endif
        endif
        
        return
    
    end subroutine solve_Ax_B
    

    ! Close out all data
    subroutine close_A()
    
        integer(4) :: idum, error, phase, info
        real(8) :: rdum
        
        ! Close Amn
        info = mkl_sparse_destroy(Amn)

        ! Clear Pardiso    
        error = 0
        phase = -1 ! release internal memory
        call pardiso(pt, maxfct, mnum, mtype, phase, N, rdum, idum, idum, &
                     perm, nrhs, iparm, lsolve_verbose, rdum, rdum, error)
        
        ! Deallocate mbsol data
        if (allocated(perm))  deallocate(perm)
        if (allocated(V_J))   deallocate(V_J)
        if (allocated(Amn_Vp))   deallocate(Amn_Vp)

    end subroutine close_A
    
    logical function is_not_factorized()
        is_not_factorized = not_factorized
    end function is_not_factorized
    
end module lsolve
    
