! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical primitives shared with the standard Quickstep tblite build.
MODULE gxtb_primitives
   USE, INTRINSIC :: iso_c_binding, ONLY: c_double, c_int
#if defined(__TBLITE)
   USE tblite_basis_type, ONLY: cgto_type
#endif
   IMPLICIT NONE
   PRIVATE
#if defined(__TBLITE)
   TYPE, BIND(C) :: shell
      INTEGER(c_int) :: l, n
      REAL(c_double) :: exponent(12), coefficient(12)
   END TYPE
   PUBLIC :: gxtb_overlap, gxtb_moments, gxtb_diatomic, gxtb_ewald_alpha
CONTAINS
   PURE FUNCTION gaussian(input) RESULT(output)
      TYPE(shell), INTENT(IN) :: input
      TYPE(cgto_type) :: output
      output%ang = input%l
      output%nprim = input%n
      output%alpha = input%exponent
      output%coeff = input%coefficient
   END FUNCTION

   SUBROUTINE gxtb_overlap(bra, ket, r, s) BIND(C)
      USE tblite_integral_native_integrals, ONLY: overlap_cgto
      TYPE(shell), INTENT(IN) :: bra, ket
      REAL(c_double), INTENT(IN) :: r(3)
      REAL(c_double), INTENT(OUT) :: s(*)
      CALL overlap_cgto(gaussian(bra), gaussian(ket), DOT_PRODUCT(r,r), r, 25.0_c_double, s)
   END SUBROUTINE

   SUBROUTINE gxtb_moments(bra, ket, r, s, d, q) BIND(C)
      USE tblite_integral_native_integrals, ONLY: multipole_cgto
      TYPE(shell), INTENT(IN) :: bra, ket
      REAL(c_double), INTENT(IN) :: r(3)
      REAL(c_double), INTENT(OUT) :: s(*), d(*), q(*)
      CALL multipole_cgto(gaussian(bra), gaussian(ket), DOT_PRODUCT(r,r), r, 25.0_c_double, s, d, q)
   END SUBROUTINE

   SUBROUTINE gxtb_diatomic(lb, la, r, factor, block) BIND(C)
      USE tblite_integral_diat_trafo, ONLY: diat_trafo, diat_trafo_cache, setup_diat_trafo
      INTEGER(c_int), VALUE :: lb, la
      REAL(c_double), INTENT(IN) :: r(3), factor(3)
      REAL(c_double), INTENT(INOUT) :: block((lb+1)**2,(la+1)**2)
      TYPE(diat_trafo_cache) :: frame
      CALL setup_diat_trafo(frame, r, INT(lb), INT(la), .FALSE.)
      CALL diat_trafo(frame, factor(1), factor(2), factor(3), block)
   END SUBROUTINE

   SUBROUTINE gxtb_ewald_alpha(cell, alpha) BIND(C)
      USE tblite_coulomb_ewald, ONLY: get_alpha
      REAL(c_double), INTENT(IN) :: cell(3,3)
      REAL(c_double), INTENT(OUT) :: alpha
      CALL get_alpha(cell, alpha, .FALSE.)
   END SUBROUTINE

#endif
END MODULE gxtb_primitives
