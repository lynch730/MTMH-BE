include 'mkl_blas.f90'

! Module for holding matlab boltzmann solution
module solver_common

    ! Load Modules
    use timers, only: timer_start, timer_start_all, timer_stop_all, timer_stop_sol
    use mdat, only: N, Neps, zin_N, N_gfrac, grid_EC, &
                    grid_E_INT_MASS, grid_E_INT_ENERGY, &
                    grid_L_is_0, Nrates, rates_zid, &
                    rates_integral, NL0, zin_field
    use mbsol, only: init_mbsol, compute_dxdt, close_mbsol
    use lsolve, only: is_not_factorized, solve_Ax_B
    use bz_solution, only: new_bz, reset_bz, &
                             close_bz_solution

    ! Custom Types
    implicit none
    logical :: is_qss
    real(8) :: mean_energy, mass_residual
    integer(4) :: meedf_verbose = 0 ! 0=no silent, 1=module, 2= all submodules
    
    ! Public Variables
    real(8), allocatable :: Xeedf(:), rates_red(:, :), rates_raw(:, :), bz_0(:) ! nu_i is number of L=0 terms
    real(8), allocatable :: Xall(:,:), time(:), EN(:)
    integer :: Nt
    
    contains
    
    ! qss solver param indices
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
    ! 13:20 -> Reserved!
    ! 21:21+Nfrac -> Base species gas fractions
    
    ! Initalize qss solver, load data, etc. 
    subroutine init_solver(fname, input_is_qss)
        character(len=*), intent(in) :: fname
        logical, intent(in) :: input_is_qss
        
        ! Initialize for particular use-case
        is_qss = input_is_qss
        call init_mbsol(fname, is_qss)
        
        ! Reallocate 
        if (allocated(Xeedf)) deallocate(Xeedf)
        if (allocated(rates_red)) deallocate(rates_red)
        if (allocated(rates_raw)) deallocate(rates_raw)
        if (allocated(bz_0)) deallocate(bz_0)
        allocate(Xeedf(N)) 
        allocate(rates_red(Nrates, NL0))
        allocate(rates_raw(Nrates, NL0))
        allocate(bz_0(zin_N))
        
        ! Set to initial State
        reset_bz = .true.
        Xeedf(:) = 0.0d0
        mean_energy = 0.0d0
        
        return
        
    end subroutine init_solver
    
    
    ! Compute moments, store key data in module vars
    subroutine eedf_moments(bz)
        real(8), intent(in) :: bz(:)
        real(8) :: FL0(Neps, NL0)
        
        ! Reshape
        FL0 = reshape(Xeedf(1:(NL0*Neps)), (/Neps, NL0/))

        ! All process rates
        call dgemm('n', 'n', Nrates, NL0, Neps, 1.0d0, rates_integral, &
                    Nrates, FL0, Neps, 0.0d0, rates_red, Nrates)

        ! Scale by number densities
        rates_raw = spread( bz(rates_zid), 2, NL0) * rates_red
        
        ! Mass and Energy Moments
        mass_residual= abs(mass_moment()-1.0d0)
        mean_energy = energy_moment()
        
        return
    end subroutine eedf_moments
    
    
    ! Test Convergence of iteration
    logical function test_convergence(a, b, tol)
        real(8), intent(in) :: a(:), b(:), tol(3)
        real(8) :: abmax
        integer :: i 
        test_convergence = not(any( abs(a-b) > max(tol(1)*abs(a), tol(2)) ))
        !test_convergence = not(any( abs( log10(abs(a))-log10(abs(b)) ) > max(tol(1)*abs(log10(abs(a))), log10(tol(2)) ) ))
        
        !if (not(test_convergence)) then
        !  do i = 1, size(a, 1)
        !    abmax = max(tol(1)*abs(a(i)), tol(2))
        !    if (abs(a(i)-b(i)) > abmax) then
        !       write(*,*) i, abs(a(i)-b(i)), tol(1)*abs(a(i))
        !    end if
        !  end do  
        !end if
        
    end function test_convergence

    
    ! Initalize EEDF
    subroutine initialize_eedf(Te_eV)
        real(8), intent(in) :: Te_eV
        real(8) :: mmoment
        real(8), parameter :: sqrt_pi_2 = 1.128379167095513d0 ! 2/sqrt(pi)
        Xeedf(:) = 0.0d0
        Xeedf(1:Neps) = sqrt_pi_2 * (Te_eV**-1.5d0) * exp(-grid_EC/Te_eV);
        mmoment = mass_moment()
        Xeedf(1:Neps) = Xeedf(1:Neps) / mmoment
        return
    end subroutine initialize_eedf
    
    
    ! Integral mass of X vector
     double precision function mass_moment()
        integer :: i
        mass_moment = dot_product(Xeedf(1:Neps), grid_E_INT_MASS)
     end function mass_moment
     
     
    ! Mean energy 
    double precision function energy_moment()
        energy_moment = dot_product(Xeedf(1:Neps), grid_E_INT_ENERGY)
    end function energy_moment
    
    
    ! Close all data
    subroutine close_solver()
        if (allocated(Xeedf)) deallocate(Xeedf)
        if (allocated(rates_red)) deallocate(rates_red)
        if (allocated(rates_raw)) deallocate(rates_raw)
        if (allocated(Xall)) deallocate(Xall)
        if (allocated(time)) deallocate(time)
        if (allocated(EN)) deallocate(EN)
        call close_mbsol()
        call close_bz_solution()
    end subroutine close_solver
    
    
end module solver_common
    
