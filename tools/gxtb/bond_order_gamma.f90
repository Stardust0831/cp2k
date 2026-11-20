! SPDX-License-Identifier: GPL-2.0-or-later
! Closed-shell reference contraction on a real, explicit Gamma supercell.
! pair_kernel contains symmetric atom-pair weights in hartree, including
! the desired image convention. owner maps each AO to its 1-based atom.
module gxtb_bond_order_gamma
   use, intrinsic :: iso_c_binding, only: c_double, c_int
   implicit none
   private
   public :: gxtb_gamma_bond_order
contains
   subroutine gxtb_gamma_bond_order(nao, nat, owner, pair_kernel, density, overlap, &
                                    energy, fock, overlap_response) bind(C)
      integer(c_int), value :: nao, nat
      integer(c_int), intent(in) :: owner(nao)
      real(c_double), intent(in) :: pair_kernel(nat, nat), density(nao, nao), overlap(nao, nao)
      real(c_double), intent(out) :: energy, fock(nao, nao), overlap_response(nao, nao)
      real(c_double), allocatable :: ps(:, :), response(:, :)
      integer :: iao, jao

      ! E = (1/2) sum_ij g_ij (PS)_ij (PS)_ji, using spin-summed P.
      ps = matmul(density, overlap)
      allocate (response(nao, nao))
      do jao = 1, nao
         do iao = 1, nao
            response(iao, jao) = pair_kernel(owner(iao), owner(jao))*ps(jao, iao)
         end do
      end do
      energy = 0.5_c_double*sum(ps*response)

      ! Only symmetric variations of the real AO matrices are independent.
      fock = matmul(response, overlap)
      fock = 0.5_c_double*(fock + transpose(fock))
      overlap_response = matmul(density, response)
      overlap_response = 0.5_c_double*(overlap_response + transpose(overlap_response))
   end subroutine gxtb_gamma_bond_order
end module gxtb_bond_order_gamma
