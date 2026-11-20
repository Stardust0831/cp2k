! SPDX-License-Identifier: LGPL-3.0-or-later
module tblite_cp2k_scf
   use iso_c_binding
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use tblite_container, only: container_cache
   use tblite_coulomb_cache, only: coulomb_cache
   use tblite_integral_type, only: integral_type
   use tblite_scf_potential, only: potential_type
   use tblite_wavefunction, only: wavefunction_type
   use tblite_xtb_calculator, only: xtb_calculator
   implicit none
   private
   public :: cp2k_save_scf_solver, run_cp2k_scf

   ! Matches src/gxtb_scf.h. All buffers remain owned by this force evaluation.
   type, bind(C) :: problem
      integer(c_int) :: nao, nsh, nat, nspin, max_iter
      real(c_double) :: nel(2), energy_tol, residual_tol
      type(c_ptr) :: ao_shell, shell_atom, s, h0, dip_int, quad_int, neutral_shell
      type(c_ptr) :: p, c, eps, occ, qsh, qat, dip, quad
      type(c_ptr) :: vat, vsh, vao, vdip, vquad, exchange, energy, weighted_p
      type(c_ptr) :: z, pbc, xyz, cell, geometry_gradient, geometry_strain
      type(c_ptr) :: cn, gamma, tau, sd, dd, sq, dq, qq
   end type
   interface
      function cp2k_save_scf_solver() bind(C) result(f)
         import :: c_funptr
         type(c_funptr) :: f
      end function
   end interface
   abstract interface
      function solve(data, bytes, msg, len) bind(C) result(cycles)
         import :: problem, c_int, c_char, c_size_t
         type(problem), intent(inout) :: data
         integer(c_size_t), value :: bytes
         character(kind=c_char), intent(out) :: msg(*)
         integer(c_size_t), value :: len
         integer(c_int) :: cycles
      end function
   end interface
contains
   subroutine run_cp2k_scf(mol,calc,wfn,ints,pot,ccache,energy,weighted,cycles,error,gradient,sigma)
      type(structure_type), target, intent(in) :: mol
      type(xtb_calculator), target, intent(in) :: calc
      type(wavefunction_type), target, intent(inout) :: wfn
      type(integral_type), target, intent(in) :: ints
      type(potential_type), target, intent(inout) :: pot
      type(container_cache), target, intent(in) :: ccache
      real(wp), contiguous, target, intent(inout) :: energy(:)
      real(wp), allocatable, target, intent(out) :: weighted(:,:,:)
      integer, intent(out) :: cycles
      type(error_type), allocatable, intent(out) :: error
      real(wp), optional, intent(inout) :: gradient(:,:), sigma(:,:)
      integer(c_int), target :: numbers(mol%nat), pbc(3)
      real(wp), target :: cell(3,3), geom_grad(3,mol%nat), geom_sigma(3,3)
      type(coulomb_cache), pointer :: es
      type(problem) :: data
      procedure(solve), pointer :: driver
      integer :: i
      character(kind=c_char) :: msg(256)
      select type(raw => ccache%raw)
      type is(coulomb_cache)
         es => raw
      class default
         call fatal_error(error,"Missing g-xTB electrostatic integrals")
         return
      end select
      if (allocated(calc%interactions)) then
         call fatal_error(error,"Extra interaction containers are unsupported in the CP2K g-xTB solver")
         return
      end if
      allocate(weighted(calc%bas%nao,calc%bas%nao,wfn%nspin))
      numbers = mol%num(mol%id)
      pbc = 0
      pbc(:size(mol%periodic)) = merge(1,0,mol%periodic)
      cell = 0.0_wp
      if (any(pbc/=0)) cell = mol%lattice
      data=problem(calc%bas%nao,calc%bas%nsh,mol%nat,wfn%nspin,calc%iterator%max_iter, &
         wfn%nel(1:2),min(calc%iterator%econv,1.e-11_wp),min(calc%iterator%pconv,1.e-9_wp), &
         c_loc(calc%bas%ao2sh),c_loc(calc%bas%sh2at),c_loc(ints%overlap),c_loc(ints%hamiltonian), &
         c_loc(ints%dipole),c_loc(ints%quadrupole),c_loc(wfn%n0sh),c_loc(wfn%density),c_loc(wfn%coeff), &
         c_loc(wfn%emo),c_loc(wfn%focc),c_loc(wfn%qsh),c_loc(wfn%qat),c_loc(wfn%dpat),c_loc(wfn%qpat), &
         c_loc(pot%vat),c_loc(pot%vsh),c_loc(pot%vao),c_loc(pot%vdp),c_loc(pot%vqp),c_loc(pot%kao), &
         c_loc(energy),c_loc(weighted),c_loc(numbers),c_loc(pbc),c_loc(mol%xyz),c_loc(cell),c_loc(geom_grad),c_loc(geom_sigma), &
         c_loc(es%cn),c_loc(es%amat),c_loc(es%taumat),c_loc(es%amat_sd),c_loc(es%amat_dd), &
         c_loc(es%amat_sq),c_loc(es%amat_dq),c_loc(es%amat_qq))
      call c_f_procpointer(cp2k_save_scf_solver(),driver)
      cycles=driver(data,c_sizeof(data),msg,size(msg,kind=c_size_t))
      if (cycles==0) then
         do i=1,size(msg)
            if(msg(i)==c_null_char) exit
         end do
         call fatal_error(error,transfer(msg(:i-1),repeat(' ',i-1)))
      else
         if (present(gradient)) gradient = gradient + geom_grad
         if (present(sigma)) sigma = sigma + geom_sigma
      end if
   end subroutine

end module
