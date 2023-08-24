#include "fintrf.h"
include 'mkl_spblas.f90'

! Module for holding static data
module mdat_alt

    implicit none

    ! Scalars
    integer(4) :: N, Neps, NL0, xsec_Nproc, xsec_Nproc_all, zin_N, nY, nspec
    integer(4) :: Natt, Nion, Nrates, N_gfrac, zin_ensemble_type
        
	real(8), parameter :: const_KB = 1.380649d-23
    real(8), parameter :: const_VMM_PER_TD = 1.0d-21
    real(8), parameter :: const_QE = 1.602176634d-19
    real(8), parameter :: const_KB_QE = 8.617333262145179d-05
    real(8), parameter :: const_gamma = 5.930969584768013d5
    
    ! Integer Arrays
    integer(4), allocatable, dimension(:) :: jac_Iu, jac_Ju, Y_I, Y_J, zin_Zmap, &
                                             zin_base_specid, zin_low_type, zin_new_sid, zin_omega, &
                                             zin_field, zin_ionization, zin_attach, &
                                             zin_inelastic, zin_superelastic, &
                                             zin_nubar_omega, zin_elastic_th, zin_nubar, &
                                             grid_L_is_0, xsec_is_rev, xsec_is_excitation, &
                                             xsec_is_ionization, xsec_is_attachment, &
                                             xsec_is_elastic, rates_zid, zin_mass_int, &
                                             ind_bc_zero

    ! Real Arrays
    real(8), allocatable, dimension(:) :: grid_EC, grid_OC, Y_V, zin_grat, zin_de, &
                                          grid_E_INT_ENERGY, grid_E_INT_MASS                           

    ! Real Matrices
    real(8), allocatable, dimension(:, :) :: integral_attachment, integral_ionization, &
                                             rates_integral

contains
    
    ! Load All variable from mat file
    subroutine load_vars(fname)
        character(len=*), intent(in) :: fname
        mwPointer :: matOpen, matClose, mp
         
        ! Open file and read full arrays
        mp = matOpen(fname, 'r')
        if (mp .eq. 0) then
            write(*,*) 'Cant open file'
            stop
        end if
        
        ! Scalar Integers
        call mat2int_scalar(mp, 'N', N)
        call mat2int_scalar(mp, 'Neps', Neps)
        call mat2int_scalar(mp, 'NL0', NL0)
        call mat2int_scalar(mp, 'zin_N', zin_N)
        call mat2int_scalar(mp, 'xsec_Nproc', xsec_Nproc)
        call mat2int_scalar(mp, 'zin_ensemble_type', zin_ensemble_type)

        ! Integer Arrays
        call mat2int(mp, 'rates_zid', rates_zid)
        call mat2int(mp, 'xsec_is_elastic', xsec_is_elastic)
        call mat2int(mp, 'xsec_is_attachment', xsec_is_attachment)
        call mat2int(mp, 'xsec_is_ionization', xsec_is_ionization)
        call mat2int(mp, 'xsec_is_excitation', xsec_is_excitation)
        call mat2int(mp, 'xsec_is_rev', xsec_is_rev) 
        call mat2int(mp, 'grid_L_is_0', grid_L_is_0)
        call mat2int(mp, 'zin_nubar', zin_nubar)
        call mat2int(mp, 'zin_elastic_th', zin_elastic_th)
        call mat2int(mp, 'zin_nubar_omega', zin_nubar_omega)
        call mat2int(mp, 'zin_attach', zin_attach)
        call mat2int(mp, 'zin_ionization', zin_ionization)
        call mat2int(mp, 'zin_inelastic', zin_inelastic)
        call mat2int(mp, 'zin_superelastic', zin_superelastic)
        call mat2int(mp, 'zin_field', zin_field)
        call mat2int(mp, 'zin_omega', zin_omega)
        call mat2int(mp, 'zin_new_sid', zin_new_sid)
        call mat2int(mp, 'zin_low_type', zin_low_type)
        call mat2int(mp, 'zin_base_specid', zin_base_specid)
        call mat2int(mp, 'zin_Zmap', zin_Zmap)
        call mat2int(mp, 'zin_mass_int', zin_mass_int)
        call mat2int(mp, 'jac_Iu', jac_Iu)
        call mat2int(mp, 'jac_Ju', jac_Ju)
        call mat2int(mp, 'Y_I', Y_I)
        call mat2int(mp, 'Y_J', Y_J)
        call mat2int(mp, 'ind_bc_zero', ind_bc_zero)

        ! Load all reals
        call mat2real(mp, 'grid_EC', grid_EC)
        call mat2real(mp, 'grid_OC', grid_OC)
        call mat2real(mp, 'grid_E_INT_ENERGY', grid_E_INT_ENERGY)
        call mat2real(mp, 'grid_E_INT_MASS', grid_E_INT_MASS)
        call mat2real(mp, 'Y_V', Y_V)
        call mat2real(mp, 'zin_grat', zin_grat)
        call mat2real(mp, 'zin_de', zin_de)
        
        ! Store rows of Y
        nY = size(jac_Ju, 1)
        nspec = size(zin_new_sid, 1)
        N_gfrac = maxval(zin_base_specid)
        
        ! Copy Allocataed arrays
        call mat2matrix(mp, 'rates_integral', rates_integral)
        Nrates = size(rates_integral, 1)
        if (size(rates_integral, 2) /= Neps) Nrates = 0

        call mat2matrix(mp, 'integral_attachment', integral_attachment)
        Natt = size(integral_attachment, 1)
        if (size(integral_attachment, 2) /= Neps) Natt = 0
        
        call mat2matrix(mp, 'integral_ionization', integral_ionization)
        Nion = size(integral_ionization, 1)
        if (size(integral_ionization, 2) /= Neps) Nion = 0
        
        ! Close File
        if (matClose(mp) .ne. 0) then
            write(*,*) 'Error closing smat file'
            stop
        end if
        
        return 
        
    end subroutine load_vars
    
    ! Load All variable from mat file
    subroutine destroy_vars()
    
        ! Deallocate integers
        if (allocated(rates_zid)) deallocate(rates_zid)
        if (allocated(xsec_is_elastic)) deallocate(xsec_is_elastic)
        if (allocated(xsec_is_attachment)) deallocate(xsec_is_attachment)
        if (allocated(xsec_is_ionization)) deallocate(xsec_is_ionization)
        if (allocated(xsec_is_excitation)) deallocate(xsec_is_excitation)
        if (allocated(xsec_is_rev)) deallocate(xsec_is_rev)
        if (allocated(grid_L_is_0)) deallocate(grid_L_is_0)
        if (allocated(zin_nubar)) deallocate(zin_nubar)
        if (allocated(zin_elastic_th)) deallocate(zin_elastic_th)
        if (allocated(zin_nubar_omega)) deallocate(zin_nubar_omega)
        if (allocated(zin_attach)) deallocate(zin_attach)
        if (allocated(zin_ionization)) deallocate(zin_ionization)
        if (allocated(zin_inelastic)) deallocate(zin_inelastic)
        if (allocated(zin_superelastic)) deallocate(zin_superelastic)
        if (allocated(zin_field)) deallocate(zin_field)
        if (allocated(zin_omega)) deallocate(zin_omega)
        if (allocated(zin_new_sid)) deallocate(zin_new_sid)
        if (allocated(zin_low_type)) deallocate(zin_low_type)
        if (allocated(zin_base_specid)) deallocate(zin_base_specid)
        if (allocated(zin_Zmap)) deallocate(zin_Zmap)
        if (allocated(jac_Iu)) deallocate(jac_Iu)
        if (allocated(jac_Ju)) deallocate(jac_Ju)
        if (allocated(Y_I)) deallocate(Y_I)
        if (allocated(Y_J)) deallocate(Y_J)
        
        ! Deallocate Floats
        if (allocated(grid_EC)) deallocate(grid_EC)   
        if (allocated(grid_OC)) deallocate(grid_OC)   
        if (allocated(grid_E_INT_ENERGY)) deallocate(grid_E_INT_ENERGY)   
        if (allocated(grid_E_INT_MASS)) deallocate(grid_E_INT_MASS)   
        if (allocated(Y_V)) deallocate(Y_V)   
        if (allocated(integral_attachment)) deallocate(integral_attachment)   
        if (allocated(integral_ionization)) deallocate(integral_ionization)   
        if (allocated(zin_grat)) deallocate(zin_grat)   
        if (allocated(zin_de)) deallocate(zin_de)   

        return
        
    end subroutine destroy_vars
    
    ! Load an individual variable
    subroutine mat2int(mp, matvar, fvar)
        mwPointer :: pa, matGetVariable, mxGetM, mxGetPr
        mwSize :: mrows
        mwPointer, intent(in) :: mp
        character(len=*), intent(in) :: matvar
        integer(4), allocatable, intent(out) :: fvar(:)
        pa = matGetVariable(mp, matvar)
        mrows = mxGetM(pa)
        allocate(fvar(mrows))
        call mxCopyPtrToInteger4(mxGetPr(pa),fvar,mrows)
        call mxDestroyArray(pa)
    end subroutine mat2int
    
    ! Load an individual variable
    subroutine mat2int_scalar(mp, matvar, fvar)
        mwPointer, intent(in) :: mp
        character(len=*), intent(in) :: matvar
        integer(4), intent(out) :: fvar
        integer(4), allocatable:: fvara(:)
        mwPointer :: pa, matGetVariable, mxGetM, mxGetPr
        mwSize :: mrows
        pa = matGetVariable(mp, matvar)
        mrows = mxGetM(pa)
        allocate(fvara(mrows))
        call mxCopyPtrToInteger4(mxGetPr(pa),fvara,mrows)
        fvar = fvara(1)
        call mxDestroyArray(pa)
        deallocate(fvara)
    end subroutine mat2int_scalar
    
    ! Load an individual variable
    subroutine mat2real(mp, matvar, fvar)
        mwPointer :: pa, matGetVariable, mxGetM, mxGetPr
        mwSize :: mrows
        mwPointer, intent(in) :: mp
        character(len=*), intent(in) :: matvar
        real(8), allocatable, intent(out) :: fvar(:)
        pa = matGetVariable(mp, matvar)
        mrows = mxGetM(pa)
        allocate(fvar(mrows))
        call mxCopyPtrToReal8(mxGetPr(pa),fvar,mrows)
        call mxDestroyArray(pa)
    end subroutine mat2real
    
    ! Load an individual matrix
    subroutine mat2matrix(mp, matvar, fvar)
        mwPointer, intent(in) :: mp
        character(len=*), intent(in) :: matvar
        real(8), allocatable, intent(out) :: fvar(:,:)
        mwPointer :: pa, matGetVariable, mxGetM, mxGetN, mxGetPr
        mwSize :: m, n
        pa = matGetVariable(mp, matvar)
        m = mxGetM(pa)
        n = mxGetN(pa)
        allocate(fvar(m, n))
        call mxCopyPtrToReal8(mxGetPr(pa),fvar,m*n)
        call mxDestroyArray(pa)
    end subroutine mat2matrix
    
end module mdat_alt
