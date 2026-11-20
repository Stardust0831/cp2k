! SPDX-License-Identifier: LGPL-3.0-or-later
module tblite_disp_revd4
   use iso_c_binding, only : c_ptr, c_null_ptr, c_null_char, c_associated, c_int, c_double, c_size_t, c_char
   use iso_fortran_env, only : error_unit
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_container_cache, only : container_cache
   use tblite_disp_type, only : dispersion_type
   use tblite_scf_info, only : scf_info, atom_resolved
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private
   public :: revd4_dispersion, new_revd4

   type, extends(dispersion_type) :: revd4_dispersion
   contains
      procedure :: update, variable_info, get_engrad, get_energy, get_potential, get_gradient
   end type
   type :: geometry
      type(c_ptr) :: handle = c_null_ptr
      real(wp), allocatable :: atm(:), grad(:, :)
      real(wp) :: sigma(3,3)
   contains
      final :: release
   end type
   interface
      function prepare(n,z,xyz,cell,pbc,msg,len) bind(C,name="cp2k_save_d4_new") result(p)
         import :: c_int, c_double, c_ptr, c_char, c_size_t
         integer(c_int), value :: n
         integer(c_int), intent(in) :: z(*), pbc(3)
         real(c_double), intent(in) :: xyz(3,*), cell(3,3)
         character(kind=c_char), intent(out) :: msg(*)
         integer(c_size_t), value :: len
         type(c_ptr) :: p
      end function
      subroutine free(p) bind(C,name="cp2k_save_d4_free")
         import :: c_ptr
         type(c_ptr), value :: p
      end subroutine
      subroutine scf(p,q,e,v) bind(C,name="cp2k_save_d4_scf")
         import :: c_ptr, c_double
         type(c_ptr), value :: p
         real(c_double), intent(in) :: q(*)
         real(c_double), intent(out) :: e(*), v(*)
      end subroutine
      subroutine grad2(p,q,g,s) bind(C,name="cp2k_save_d4_derivatives")
         import :: c_ptr, c_double
         type(c_ptr), value :: p
         real(c_double), intent(in) :: q(*)
         real(c_double), intent(out) :: g(3,*), s(3,3)
      end subroutine
      subroutine atm(p,e,g,s) bind(C,name="cp2k_save_d4_atm")
         import :: c_ptr, c_double
         type(c_ptr), value :: p
         real(c_double), intent(out) :: e(*), g(3,*), s(3,3)
      end subroutine
   end interface
contains
   subroutine new_revd4(self, mol, error)
      type(revd4_dispersion), intent(out) :: self
      type(structure_type), intent(in) :: mol
      type(error_type), allocatable, intent(out) :: error
      self%label = "self-consistent g-xTB 2.0.1 D4Srev (CP2K)"
      if (any(mol%num < 1) .or. any(mol%num > 103)) call fatal_error(error, "D4Srev requires Z=1..103")
   end subroutine
   subroutine release(g)
      type(geometry), intent(inout) :: g
      if (c_associated(g%handle)) call free(g%handle)
      g%handle = c_null_ptr
   end subroutine
   pure function variable_info(self) result(info)
      class(revd4_dispersion), intent(in) :: self
      type(scf_info) :: info
      info = scf_info(charge=atom_resolved)
   end function
   subroutine update(self,mol,cache)
      class(revd4_dispersion), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(container_cache), intent(inout) :: cache
      character(kind=c_char) :: msg(256)
      integer :: i
      real(wp) :: cell(3,3)
      integer(c_int) :: pbc(3)
      if (allocated(cache%raw)) deallocate(cache%raw)
      allocate(geometry :: cache%raw)
      select type(g => cache%raw)
      type is(geometry)
         cell = 0.0_wp
         pbc = 0
         pbc(:size(mol%periodic)) = merge(1,0,mol%periodic)
         if (any(pbc/=0)) cell = mol%lattice
         g%handle = prepare(int(mol%nat,c_int),int(mol%num(mol%id),c_int),mol%xyz,cell,pbc,msg,size(msg,kind=c_size_t))
         allocate(g%atm(mol%nat),g%grad(3,mol%nat))
         if (c_associated(g%handle)) then
            call atm(g%handle,g%atm,g%grad,g%sigma)
         else
            do i = 1, size(msg)
               if (msg(i) == c_null_char) exit
            end do
            write(error_unit, '(a,*(a))') 'D4Srev geometry: ', msg(:i-1)
            g%atm = ieee_value(0.0_wp,ieee_quiet_nan)
            g%grad = g%atm(1)
            g%sigma = g%atm(1)
         end if
      end select
   end subroutine
   subroutine get_engrad(self,mol,cache,energies,gradient,sigma)
      class(revd4_dispersion), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(container_cache), intent(inout) :: cache
      real(wp), intent(inout) :: energies(:)
      real(wp), contiguous, intent(inout), optional :: gradient(:,:), sigma(:,:)
      select type(g => cache%raw)
      type is(geometry)
         energies = energies + g%atm
         if (present(gradient)) gradient = gradient + g%grad
         if (present(sigma)) sigma = sigma + g%sigma
      end select
   end subroutine
   subroutine get_energy(self,mol,cache,wfn,energies)
      class(revd4_dispersion), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(container_cache), intent(inout) :: cache
      type(wavefunction_type), intent(in) :: wfn
      real(wp), intent(inout) :: energies(:)
      real(wp) :: e(mol%nat),v(mol%nat)
      select type(g => cache%raw)
      type is(geometry)
         e = ieee_value(0.0_wp,ieee_quiet_nan)
         if (c_associated(g%handle)) call scf(g%handle,wfn%qat(:,1),e,v)
         energies = energies + e
      end select
   end subroutine
   subroutine get_potential(self,mol,cache,wfn,pot)
      class(revd4_dispersion), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(container_cache), intent(inout) :: cache
      type(wavefunction_type), intent(in) :: wfn
      type(potential_type), intent(inout) :: pot
      real(wp) :: e(mol%nat),v(mol%nat)
      select type(g => cache%raw)
      type is(geometry)
         v = ieee_value(0.0_wp,ieee_quiet_nan)
         if (c_associated(g%handle)) call scf(g%handle,wfn%qat(:,1),e,v)
         pot%vat(:,1) = pot%vat(:,1) + v
      end select
   end subroutine
   subroutine get_gradient(self,mol,cache,wfn,gradient,sigma)
      class(revd4_dispersion), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(container_cache), intent(inout) :: cache
      type(wavefunction_type), intent(in) :: wfn
      real(wp), contiguous, intent(inout) :: gradient(:,:), sigma(:,:)
      real(wp) :: d(3,mol%nat),s(3,3)
      select type(g => cache%raw)
      type is(geometry)
         if (c_associated(g%handle)) then
            call grad2(g%handle,wfn%qat(:,1),d,s)
            gradient = gradient + d
            sigma = sigma + s
         else
            gradient = ieee_value(0.0_wp,ieee_quiet_nan)
            sigma = ieee_value(0.0_wp,ieee_quiet_nan)
         end if
      end select
   end subroutine
end module
