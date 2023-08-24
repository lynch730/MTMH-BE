! Module for holding matlab boltzmann solution

! Summary:
!  init_lsove:      set N, init parameters
!  build_sparse_A:   
!  pointers_A:        c pointers to A matrix

module timers
    
    ! Load Modules
    use ifport, only: dclock
    implicit none
    
    ! Private
    real(8), allocatable, private :: t_rhs(:), t_sol(:), t_jac(:), t_lum(:), t_all(:)
    integer, allocatable, private :: n_rhs(:), n_sol(:), n_jac(:), n_lum(:), n_all(:)
    integer, private :: tdx, Ntimers ! time Index
    real(8), private :: tstart, tstart_all
    
    contains

    
    ! Initialize Timers
    subroutine reset_timers(N) 
        integer, intent(in) :: N
        call close_timers()
        allocate(t_rhs(N), t_sol(N), t_jac(N), t_lum(N), t_all(N))
        allocate(n_rhs(N), n_sol(N), n_jac(N), n_lum(N), n_all(N))
        Ntimers = N
        t_rhs(:) = 0.0d0
        t_sol(:) = 0.0d0
        t_jac(:) = 0.0d0
        t_lum(:) = 0.0d0
        t_all(:) = 0.0d0
        n_rhs(:) = 0.0d0
        n_sol(:) = 0.0d0
        n_jac(:) = 0.0d0
        n_lum(:) = 0.0d0
        n_all(:) = 0.0d0
        tdx = 1 ! Inital Time index
        return
    end subroutine reset_timers
    
    
    ! Set Index
    subroutine set_timer_index(i)
        integer, intent(in) :: i
        tdx = i
        return
    end subroutine set_timer_index    
    
    
    ! Timer Start
    subroutine timer_start()
        tstart = dclock()
    end subroutine timer_start
    
    
    ! Timer Start
    subroutine timer_start_all()
        tstart_all = dclock()
    end subroutine timer_start_all
    
    ! Generation of RHS, 
    subroutine timer_stop_rhs()
        t_rhs(tdx) = t_rhs(tdx) + (dclock()-tstart)
        n_rhs(tdx) = n_rhs(tdx) + 1
        return
    end subroutine timer_stop_rhs
        
    
    ! A\b Solution
    subroutine timer_stop_sol()
        t_sol(tdx) = t_sol(tdx) + (dclock()-tstart)
        n_sol(tdx) = n_sol(tdx) + 1
        return
    end subroutine timer_stop_sol
    
    
    ! Jacobian Generation
    subroutine timer_stop_jac()
        t_jac(tdx) = t_jac(tdx) + (dclock()-tstart)
        n_jac(tdx) = n_jac(tdx) + 1
        return
    end subroutine timer_stop_jac
    
    
    ! LU Decomposition
    subroutine timer_stop_lum()
        t_lum(tdx) = t_lum(tdx) + (dclock()-tstart)
        n_lum(tdx) = n_lum(tdx) + 1
        return
    end subroutine timer_stop_lum
    
    
    ! Total Time
    subroutine timer_stop_all()
        t_all(tdx) = t_all(tdx) + (dclock()-tstart_all)
        n_all(tdx) = n_all(tdx) + 1
        return
    end subroutine timer_stop_all
    
    
    ! Print Timers
    subroutine print_timers()
        integer :: i
        real(8) :: tres, dn
        character(len=12) :: fmt1 = "(F8.3','2X)"
        character(len=12) :: fmt2 = "(F8.6','2X)"
        character(len=8 ) :: fmti = "(I8',')"
        write(*,*) 'Timing Results: Ntimers=', Ntimers
        write(*,*) '------------------------------------------'
        do i = 1, Ntimers
            dn = dble(n_all(i)) / 1000.0d0 ! Convert to millisecond
            tres = t_all(i) - t_lum(i) - t_jac(i) - t_sol(i) - t_rhs(i)
            write(*, '(I3", ")', advance='no') i
            write(*, fmt1, advance='no') t_all(i)/dn
            write(*, fmt1, advance='no') t_lum(i)/dn
            write(*, fmt1, advance='no') t_jac(i)/dn
            write(*, fmt1, advance='no') t_sol(i)/dn
            write(*, fmt1, advance='no') t_rhs(i)/dn
            write(*, fmt1, advance='no') tres/dn
            write(*, fmt2, advance='no') t_lum(i)/dn/dble(n_lum(i))
            write(*, fmt2, advance='no') t_jac(i)/dn/dble(n_jac(i))
            write(*, fmt2, advance='no') t_sol(i)/dn/dble(n_sol(i))
            write(*, fmt2, advance='no') t_rhs(i)/dn/dble(n_rhs(i))
            write(*, fmti, advance='no') n_lum(i)
            write(*, fmti, advance='no') n_jac(i)
            write(*, fmti, advance='no') n_sol(i)
            write(*, fmti, advance='no') n_rhs(i)
            write(*, fmti, advance='no') n_all(i)
            write(*,*)
        end do
        write(*,*) '------------------------------------------'
        return
    end subroutine print_timers
    
    
    ! Reset Timers
    subroutine close_timers()
        if (allocated(t_rhs)) deallocate(t_rhs)
        if (allocated(t_sol)) deallocate(t_sol)
        if (allocated(t_jac)) deallocate(t_jac)
        if (allocated(t_lum)) deallocate(t_lum)
        if (allocated(t_all)) deallocate(t_all)
        if (allocated(n_rhs)) deallocate(n_rhs)
        if (allocated(n_sol)) deallocate(n_sol)
        if (allocated(n_jac)) deallocate(n_jac)
        if (allocated(n_lum)) deallocate(n_lum)
        if (allocated(n_all)) deallocate(n_all)
        return
    end subroutine close_timers
    
    
end module timers
    
