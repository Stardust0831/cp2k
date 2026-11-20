! SPDX-License-Identifier: GPL-2.0-or-later
MODULE gxtb_scf_linalg
   USE, INTRINSIC :: iso_c_binding, ONLY: c_double, c_int
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: gxtb_cholesky, gxtb_orbitals, gxtb_density, gxtb_product, gxtb_matvec
   PUBLIC :: gxtb_linear_solve
CONTAINS
   FUNCTION gxtb_linear_solve(n, matrix, rhs) RESULT(info) BIND(C)
      INTEGER(c_int), VALUE :: n
      REAL(c_double), INTENT(INOUT) :: matrix(n,n), rhs(n)
      INTEGER(c_int) :: info
      INTEGER :: status, pivots(n)
      CALL DGESV(INT(n), 1, matrix, INT(n), pivots, rhs, INT(n), status)
      info = INT(status, c_int)
   END FUNCTION

   SUBROUTINE gxtb_matvec(rows, cols, trans, scale, matrix, x, y) BIND(C)
      INTEGER(c_int), VALUE :: rows, cols, trans
      REAL(c_double), VALUE :: scale
      REAL(c_double), INTENT(IN) :: matrix(rows,cols), x(*)
      REAL(c_double), INTENT(INOUT) :: y(*)
      CALL DGEMV(MERGE("T", "N", trans /= 0), INT(rows), INT(cols), scale, matrix, INT(rows), x, 1, 1.0_c_double, y, 1)
   END SUBROUTINE

   SUBROUTINE gxtb_product(n, transa, transb, beta, a, b, c) BIND(C)
      INTEGER(c_int), VALUE :: n, transa, transb
      REAL(c_double), VALUE :: beta
      REAL(c_double), INTENT(IN) :: a(n,n), b(n,n)
      REAL(c_double), INTENT(INOUT) :: c(n,n)
      CALL DGEMM(MERGE("T", "N", transa /= 0), MERGE("T", "N", transb /= 0), INT(n), INT(n), INT(n), &
                 1.0_c_double, a, INT(n), b, INT(n), beta, c, INT(n))
   END SUBROUTINE

   FUNCTION gxtb_cholesky(n, overlap) RESULT(info) BIND(C)
      INTEGER(c_int), VALUE :: n
      REAL(c_double), INTENT(INOUT) :: overlap(n, n)
      INTEGER(c_int) :: info
      INTEGER :: status
      CALL DPOTRF("L", INT(n), overlap, INT(n), status)
      info = INT(status, c_int)
   END FUNCTION

   FUNCTION gxtb_orbitals(n, chol, fock, energies) RESULT(info) BIND(C)
      INTEGER(c_int), VALUE :: n
      REAL(c_double), INTENT(IN) :: chol(n, n)
      REAL(c_double), INTENT(INOUT) :: fock(n, n)
      REAL(c_double), INTENT(OUT) :: energies(n)
      INTEGER(c_int) :: info
      REAL(c_double) :: query(1)
      REAL(c_double), ALLOCATABLE :: work(:)
      INTEGER :: status, order

      order = INT(n)
      CALL DSYGST(1, "L", order, fock, order, chol, order, status)
      IF (status == 0) CALL DSYEV("V", "L", order, fock, order, energies, query, -1, status)
      IF (status == 0) THEN
         ALLOCATE (work(CEILING(query(1))))
         CALL DSYEV("V", "L", order, fock, order, energies, work, SIZE(work), status)
      END IF
      IF (status == 0) CALL DTRSM("L", "L", "T", "N", order, order, 1.0_c_double, chol, order, fock, order)
      info = INT(status, c_int)
   END FUNCTION

   SUBROUTINE gxtb_density(n, coeff, weights, p) BIND(C)
      INTEGER(c_int), VALUE :: n
      REAL(c_double), INTENT(IN) :: coeff(n, n), weights(n)
      REAL(c_double), INTENT(OUT) :: p(n, n)
      REAL(c_double), ALLOCATABLE :: weighted(:, :)
      INTEGER :: j
      ALLOCATE (weighted(n, n))
      DO j = 1, n
         weighted(:, j) = weights(j)*coeff(:, j)
      END DO
      CALL DGEMM("N", "T", INT(n), INT(n), INT(n), 1.0_c_double, weighted, INT(n), coeff, INT(n), 0.0_c_double, p, INT(n))
   END SUBROUTINE
END MODULE
