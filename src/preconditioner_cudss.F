!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Experimental sparse direct OT preconditioner. The factorization is retained until rebuild.
! **************************************************************************************************
MODULE preconditioner_cudss
   USE ISO_C_BINDING,                   ONLY: C_ASSOCIATED,&
                                              C_DOUBLE,&
                                              C_INT,&
                                              C_NULL_PTR,&
                                              C_PTR
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_desymmetrize, dbcsr_get_info, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_readonly_start, dbcsr_iterator_start, dbcsr_iterator_stop, &
        dbcsr_iterator_type, dbcsr_release, dbcsr_type
   USE cp_dbcsr_contrib,                ONLY: dbcsr_reserve_all_blocks
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_release,&
                                              mp_para_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner_cudss'
   TYPE, PUBLIC :: cudss_preconditioner_type
      TYPE(C_PTR) :: context = C_NULL_PTR
      TYPE(mp_para_env_type), POINTER :: para_env => NULL()
      INTEGER :: n = 0
   END TYPE cudss_preconditioner_type
   PUBLIC :: cudss_preconditioner_create, cudss_preconditioner_release, cudss_preconditioner_apply
   INTERFACE cudss_preconditioner_apply
      MODULE PROCEDURE apply_cudss_fm, apply_cudss_dbcsr
   END INTERFACE
#if defined(__CUDSS)
   INTERFACE
      SUBROUTINE cp_cudss_create(context, comm, n, nnz, batch_size, rows, cols, values) BIND(C)
         IMPORT C_PTR, C_INT, C_DOUBLE
         TYPE(C_PTR), INTENT(OUT) :: context
         INTEGER(C_INT), VALUE :: comm, n, nnz, batch_size
         INTEGER(C_INT), INTENT(IN) :: rows(*), cols(*)
         REAL(C_DOUBLE), INTENT(IN) :: values(*)
      END SUBROUTINE cp_cudss_create
      SUBROUTINE cp_cudss_solve(context, nrhs, values) BIND(C)
         IMPORT C_PTR, C_INT, C_DOUBLE
         TYPE(C_PTR), VALUE :: context
         INTEGER(C_INT), VALUE :: nrhs
         REAL(C_DOUBLE), INTENT(INOUT) :: values(*)
      END SUBROUTINE cp_cudss_solve
      SUBROUTINE cp_cudss_destroy(context) BIND(C)
         IMPORT C_PTR
         TYPE(C_PTR), INTENT(INOUT) :: context
      END SUBROUTINE cp_cudss_destroy
   END INTERFACE
#endif
CONTAINS

! **************************************************************************************************
!> \brief Use at most 1024 columns, targeting 128 MiB per RHS buffer.
!> \param n number of rows
!> \return maximum number of RHS columns per batch
! **************************************************************************************************
   INTEGER FUNCTION rhs_batch_size(n) RESULT(batch_size)
      INTEGER, INTENT(IN)                                :: n

      batch_size = MAX(1, MIN(1024, (16*1024*1024)/n))
   END FUNCTION rhs_batch_size

! **************************************************************************************************
!> \brief Collect the stored upper triangle without constructing a dense operator.
!> \param prec solver state
!> \param matrix SPD preconditioning operator
!> \param para_env communicator owning the operator
! **************************************************************************************************
   SUBROUTINE cudss_preconditioner_create(prec, matrix, para_env)
      TYPE(cudss_preconditioner_type), INTENT(INOUT) :: prec
      TYPE(dbcsr_type), INTENT(IN) :: matrix
      TYPE(mp_para_env_type), POINTER :: para_env
#if defined(__CUDSS)
      TYPE(dbcsr_type) :: full
      TYPE(dbcsr_iterator_type) :: iter
      REAL(dp), POINTER :: block(:, :)
      INTEGER(C_INT), ALLOCATABLE :: rows(:), cols(:)
      REAL(dp), ALLOCATABLE :: values(:)
      INTEGER :: n, ncol, nnz, pass, i, j, row, col, nr, nc, handle
      LOGICAL :: transposed
      REAL(dp) :: v
      CHARACTER(LEN=*), PARAMETER :: routineN = 'cudss_preconditioner_create'

      CALL timeset(routineN, handle)
      CALL cudss_preconditioner_release(prec)
      CALL dbcsr_get_info(matrix, nfullrows_total=n, nfullcols_total=ncol)
      CPASSERT(n > 0 .AND. n == ncol)
      CALL dbcsr_desymmetrize(matrix, full)
      DO pass = 1, 2
         nnz = 0
         CALL dbcsr_iterator_start(iter, full)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, block=block, transposed=transposed, &
                                           row_offset=row, col_offset=col, row_size=nr, col_size=nc)
            DO j = 1, nc
               DO i = 1, nr
                  IF (row + i > col + j) CYCLE
                  IF (transposed) THEN
                     v = block(j, i)
                  ELSE
                     v = block(i, j)
                  END IF
                  IF (v == 0.0_dp .AND. row + i /= col + j) CYCLE
                  nnz = nnz + 1
                  IF (pass == 2) THEN
                     rows(nnz) = INT(row + i - 2, C_INT)
                     cols(nnz) = INT(col + j - 2, C_INT)
                     values(nnz) = v
                  END IF
               END DO
            END DO
         END DO
         CALL dbcsr_iterator_stop(iter)
         IF (pass == 1) ALLOCATE (rows(nnz), cols(nnz), values(nnz))
      END DO
      CALL dbcsr_release(full)
      CALL cp_cudss_create(prec%context, INT(para_env%get_handle(), C_INT), INT(n, C_INT), &
                           INT(nnz, C_INT), INT(rhs_batch_size(n), C_INT), rows, cols, values)
      DEALLOCATE (rows, cols, values)
      prec%para_env => para_env
      CALL prec%para_env%retain()
      prec%n = n
      CALL timestop(handle)
#else
      MARK_USED(prec)
      MARK_USED(matrix)
      MARK_USED(para_env)
      CPABORT('PRECOND_SOLVER CUDSS requires a build with CP2K_USE_CUDSS and MPI')
#endif
   END SUBROUTINE cudss_preconditioner_create

! **************************************************************************************************
!> \brief Release GPU factors before their MPI environment.
!> \param prec solver state
! **************************************************************************************************
   SUBROUTINE cudss_preconditioner_release(prec)
      TYPE(cudss_preconditioner_type), INTENT(INOUT) :: prec
#if defined(__CUDSS)
      IF (C_ASSOCIATED(prec%context)) CALL cp_cudss_destroy(prec%context)
#endif
      IF (ASSOCIATED(prec%para_env)) CALL mp_para_env_release(prec%para_env)
      prec%n = 0
   END SUBROUTINE cudss_preconditioner_release

! **************************************************************************************************
!> \brief Apply retained factors to a batch of dense right hand sides.
!> \param prec solver state
!> \param values replicated batch; overwritten by the solution
! **************************************************************************************************
   SUBROUTINE solve_batch(prec, values)
      TYPE(cudss_preconditioner_type), INTENT(IN) :: prec
      REAL(dp), INTENT(INOUT) :: values(:, :)
#if defined(__CUDSS)
      CPASSERT(C_ASSOCIATED(prec%context))
      CALL cp_cudss_solve(prec%context, INT(SIZE(values, 2), C_INT), values)
#else
      MARK_USED(prec)
      MARK_USED(values)
      CPABORT('cuDSS is not enabled')
#endif
   END SUBROUTINE solve_batch

! **************************************************************************************************
!> \brief Apply sparse factors to a distributed full matrix in bounded batches.
!> \param prec solver state
!> \param matrix_in residual
!> \param matrix_out preconditioned residual
! **************************************************************************************************
   SUBROUTINE apply_cudss_fm(prec, matrix_in, matrix_out)
      TYPE(cudss_preconditioner_type), INTENT(IN)        :: prec
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_in, matrix_out

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'apply_cudss_fm'

      INTEGER                                            :: count, first, handle, k, n
      REAL(dp), ALLOCATABLE                              :: batch(:, :)

      CALL timeset(routineN, handle)
      CALL cp_fm_get_info(matrix_in, nrow_global=n, ncol_global=k)
      CPASSERT(n == prec%n)
      DO first = 1, k, rhs_batch_size(n)
         count = MIN(rhs_batch_size(n), k - first + 1)
         ALLOCATE (batch(n, count))
         CALL cp_fm_get_submatrix(matrix_in, batch, start_col=first)
         CALL solve_batch(prec, batch)
         CALL cp_fm_set_submatrix(matrix_out, batch, start_col=first)
         DEALLOCATE (batch)
      END DO
      CALL timestop(handle)
   END SUBROUTINE apply_cudss_fm

! **************************************************************************************************
!> \brief Apply sparse factors directly to DBCSR orbital gradients.
!> \param prec solver state
!> \param matrix_in residual
!> \param matrix_out preconditioned residual
! **************************************************************************************************
   SUBROUTINE apply_cudss_dbcsr(prec, matrix_in, matrix_out)
      TYPE(cudss_preconditioner_type), INTENT(IN)        :: prec
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_in
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_out

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'apply_cudss_dbcsr'

      INTEGER                                            :: col, count, first, handle, hi, j, k, lo, &
                                                            n, nc, nr, row
      LOGICAL                                            :: transposed
      REAL(dp), ALLOCATABLE                              :: batch(:, :)
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)
      CALL dbcsr_get_info(matrix_in, nfullrows_total=n, nfullcols_total=k)
      CPASSERT(n == prec%n)
      CALL dbcsr_reserve_all_blocks(matrix_out)
      DO first = 1, k, rhs_batch_size(n)
         count = MIN(rhs_batch_size(n), k - first + 1)
         ALLOCATE (batch(n, count))
         batch = 0.0_dp
         CALL dbcsr_iterator_readonly_start(iter, matrix_in)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, block=block, transposed=transposed, &
                                           row_offset=row, col_offset=col, row_size=nr, col_size=nc)
            lo = MAX(col, first); hi = MIN(col + nc - 1, first + count - 1)
            DO j = lo, hi
               IF (transposed) THEN
                  batch(row:row + nr - 1, j - first + 1) = block(j - col + 1, :)
               ELSE
                  batch(row:row + nr - 1, j - first + 1) = block(:, j - col + 1)
               END IF
            END DO
         END DO
         CALL dbcsr_iterator_stop(iter)
         CALL prec%para_env%sum(batch)
         CALL solve_batch(prec, batch)
         CALL dbcsr_iterator_start(iter, matrix_out)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, block=block, transposed=transposed, &
                                           row_offset=row, col_offset=col, row_size=nr, col_size=nc)
            lo = MAX(col, first); hi = MIN(col + nc - 1, first + count - 1)
            DO j = lo, hi
               IF (transposed) THEN
                  block(j - col + 1, :) = batch(row:row + nr - 1, j - first + 1)
               ELSE
                  block(:, j - col + 1) = batch(row:row + nr - 1, j - first + 1)
               END IF
            END DO
         END DO
         CALL dbcsr_iterator_stop(iter)
         DEALLOCATE (batch)
      END DO
      CALL timestop(handle)
   END SUBROUTINE apply_cudss_dbcsr
END MODULE preconditioner_cudss
