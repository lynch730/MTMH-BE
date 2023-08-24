
! Module for holding matlab boltzmann solution
module bz_solution

    ! Load Modules
    use mdat

    ! Public Variables
    implicit none
    logical :: reset_bz

    ! Private Variables.
    real(8), private :: T, Texc, P, W, EN
    real(8), private :: N0 !
    real(8), allocatable, private :: Xk(:), bz_store(:), gfrac(:)
    
    contains
    
    
    ! Close all data
    subroutine new_bz(bz, X, params, gfrac2, Xk2)
        real(8), intent(in) :: X(:), params(5), gfrac2(:)
        real(8), intent(in), optional :: Xk2(:)
        real(8), intent(out) :: bz(zin_N)
        real(8) :: T2, P2, W2, EN2, Texc2
        real(8) :: Tref, FL0(Neps, NL0)
        real(8) :: nu_bar(1, NL0), net_nu(1, Neps)
        logical :: new_N0, new_T, new_Xk, new_gfrac
        
        ! Copy
        T2 = params(1)
        P2 = params(2)
        W2 = params(3)
        EN2 = params(4) 
        Texc2 = params(5)
        
        ! Reset bz, true on first call
        if (reset_bz) then
            call close_bz_solution()
            allocate(Xk(nspec), bz_store(zin_N), gfrac(N_gfrac))
            bz_store(:) = 0.0d0
            gfrac(:) = -1.0d0
            Xk(:) = -1.0d0
            T = -1.0d0
            Texc = -1.0d0
            P = -1.0d0
            W = -1.0d0 
            EN = -1.0d0
        end if
        bz = bz_store

        ! Sore Nk reac
        new_T = is_not_eq(T, T2)
        new_N0 = (new_T .or. is_not_eq(P, P2) )
        if (new_N0) then
            Tref = T2;
            if (Tref <= 0.0d0) then
                Tref = 300.0d0;
            end if
            N0 = P2 / (const_KB * Tref);
        end if
        
        ! Determine if custom Xk or Texc
        if (present(Xk2)) then ! either Nk passed or Nbase_spec
            Xk = Xk2
            new_gfrac = .true.
            new_Xk = .true.
        else
            new_Xk = is_not_eq(Texc, Texc2) .or. array_is_not_eq(gfrac, gfrac2)
            if (new_Xk) then
                call species_fractions(Texc2, gfrac2) ! Update Xk from Texc
            end if
        end if
    
        ! Species Nepsmber Densities
        if (new_Xk .or. new_N0) then
            bz(1:nspec) = Xk * N0
            bz(zin_elastic_th) = bz(zin_elastic_th) * T2 * const_KB_QE
        elseif (new_T) then
            bz(zin_elastic_th) = bz(zin_elastic_th) * T2 / T
        end if
    
        ! Field Characteristics
        if (is_not_eq(EN, EN2) .or. new_N0) then
            bz(zin_field) = EN2 * N0 * const_VMM_PER_TD
        end if
        if (is_not_eq(W, W2)) then
            bz(zin_omega) = W2
        end if
    
        ! Initialize Moments for Neps_i
        FL0 = reshape(X(1:(NL0*Neps)), (/Neps, NL0/))
    
        ! Attachment
        nu_bar(1, :) = 0.0d0
        if (Natt>0) then
            net_nu(1, :) = sum( spread(Xk(zin_attach), 2, Neps) * integral_attachment, 1)
            call dgemm('n', 'n', 1, NL0, Neps, 1.0d0, net_nu, 1, FL0, Neps, 0.0d0, nu_bar, 1)
        end if
        
        ! Ionization
        if (Nion>0) then
            net_nu(1, :) = sum( spread(Xk(zin_ionization), 2, Neps) * integral_ionization, 1)
            call dgemm('n', 'n', 1, NL0, Neps, 1.0d0, net_nu, 1, FL0, Neps, -1.0d0, nu_bar, 1)            
        end if
        
        ! Scale by N0
        nu_bar = nu_bar * N0 ! Scale to number density
        
        ! Apply zin_nubar
        bz(zin_nubar) = nu_bar(1,:)
        !bz(zin_nubar_omega) = nu_bar(1,:) / W2
        
        ! If qss, apply mass int on first row
        bz(zin_mass_int) = 1.0d0
        
        ! Copy new data struct
        reset_bz = .false.
        T = T2
        Texc = Texc2
        P = P2
        W = W2
        EN = EN2
        gfrac(:) = gfrac2(:)
        bz_store(:) = bz(:)
        
    end subroutine new_bz
    
    
    ! Species Fractions, sets Xk
    subroutine species_fractions(Texc_new, gfrac_new)
        real(8), intent(in) :: Texc_new, gfrac_new(:)
        real(8) :: T_exc_qe, gg(nspec)!, Nj(nspec)
        integer :: i
        
        ! Excitation Energy
        T_exc_qe = Texc_new * const_KB_QE
        
        ! Get base fractions
        Xk = gfrac_new( zin_base_specid(zin_new_sid) )
        
        ! Partition Function (size of total x-fractions)
        if (zin_ensemble_type == 0) then

            gg = zin_grat * exp(-zin_de / T_exc_qe)
            do i = 1, nspec
                if (isnan(gg(i))) then
                    gg(i) = 0.0
                end if
            enddo
            Xk(zin_inelastic) = Xk(zin_inelastic) / (1.0 + gg(zin_inelastic))
            Xk(zin_superelastic) = Xk(zin_superelastic) * gg(zin_superelastic) &
                                    / (1.0+gg(zin_superelastic))
        
        elseif (zin_ensemble_type == 1) then ! Standard approach for Ensemble
            
            write(*,*) 'Ensemble Type 1 not prepared!'
            stop
            ! Nepsmber density of each species, weighted by boltzmann stats
            !Nj = zin_grat .* exp(-zin_de / T_exc_qe) ! get relative population
            !Nj(isnan(Nj) & zin_de>0.0) = 0.0
            !Nj(isnan(Nj) & zin_de==0.0) = zin_grat(isnan(Nj) & zin_de==0.0)

            ! Get partition functions for each Nj and normalize
            !Z = reshape(accumarray(zin_Zmap, Nj(:), [numel(X_total), 1]), 1, [])
            !Nj = Nj / Z(zin.species.base_specid)
            !Nj = Nj * Ngas *  X_total( zin.species.base_specid )

            ! Apply to inelastic and superelastic terms only
            ! Attachment, elastic, and ionization assume excited states
            ! all have the same corss-sections for thest as ground state
            ! only significant at very high Texc.
            !ind = zin.species.low_type == 0 | zin.species.low_type == -2
            !Nk(ind) = Nj(zin.species.new_sid(ind))

        elseif (zin_ensemble_type == 2) then ! Standard approach for Ensemble
        
        end if
    
        return
    end subroutine species_fractions
    
    
    ! Deallocate Arrays
    subroutine close_bz_solution()
        if (allocated(Xk)) deallocate(Xk)
        if (allocated(bz_store)) deallocate(bz_store)
        if (allocated(gfrac)) deallocate(gfrac)
        return
    end subroutine close_bz_solution

    
    ! Compare if real scalars are equal
    logical function is_not_eq(a,b)
        real(8), intent(in) :: a, b
        real(8), parameter :: thresh = 4.0d0 ! within 4 epsilon are equal
        is_not_eq = (abs(a-b) >= thresh*epsilon(a))
    end function is_not_eq
        
    ! Compare if real arrays are equal
    logical function array_is_not_eq(a,b)
        real(8), intent(in) :: a(:), b(:)
        real(8), parameter :: thresh = 4.0d0 ! within 4 epsilon are equal
        array_is_not_eq = (maxval(abs(a-b)) >= thresh*epsilon(minval(abs(a))))
    end function array_is_not_eq
    
end module bz_solution
    