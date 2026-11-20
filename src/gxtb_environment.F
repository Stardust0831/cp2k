! SPDX-License-Identifier: GPL-2.0-or-later
MODULE gxtb_environment
   USE atomic_kind_types, ONLY: get_atomic_kind
   USE cp_log_handling, ONLY: cp_logger_get_default_io_unit
   USE cp_subsys_methods, ONLY: cp_subsys_create
   USE cp_subsys_types, ONLY: cp_subsys_set, cp_subsys_type
   USE distribution_1d_types, ONLY: distribution_1d_type, distribution_1d_release
   USE distribution_methods, ONLY: distribute_molecules_1d
   USE gxtb_types, ONLY: gxtb_env_type
   USE input_section_types, ONLY: section_vals_get, section_vals_get_subs_vals, section_vals_type, section_vals_val_get
   USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
   USE, INTRINSIC :: iso_c_binding, ONLY: c_associated, c_char, c_double, c_int, c_null_char, c_ptr, c_size_t
   USE kinds, ONLY: default_path_length, dp
   USE message_passing, ONLY: mp_para_env_type
#include "./base/base_uses.f90"
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: gxtb_env_type, gxtb_init, gxtb_evaluate

   INTERFACE
      FUNCTION cp_gxtb_open(path, msg, capacity) RESULT(model) BIND(C)
         IMPORT :: c_char, c_ptr, c_size_t
         CHARACTER(KIND=c_char), INTENT(IN) :: path(*)
         CHARACTER(KIND=c_char), INTENT(OUT) :: msg(*)
         INTEGER(c_size_t), VALUE :: capacity
         TYPE(c_ptr) :: model
      END FUNCTION
      FUNCTION cp_gxtb_eval(model, n, z, pos, cell, radii, forces, virial, energy, grad, strain, msg, capacity) RESULT(stat) BIND(C)
         IMPORT :: c_char, c_double, c_int, c_ptr, c_size_t
         TYPE(c_ptr), VALUE :: model
         INTEGER(c_int), VALUE :: n, forces, virial
         INTEGER(c_int), INTENT(IN) :: z(*)
         REAL(c_double), INTENT(IN) :: pos(3,*), cell(3,3), radii(2)
         REAL(c_double), INTENT(OUT) :: energy, grad(3,*), strain(3,3)
         CHARACTER(KIND=c_char), INTENT(OUT) :: msg(*)
         INTEGER(c_size_t), VALUE :: capacity
         INTEGER(c_int) :: stat
      END FUNCTION
   END INTERFACE
CONTAINS

   SUBROUTINE gxtb_init(env, root, fe, sub, para, motion)
      TYPE(gxtb_env_type), INTENT(INOUT) :: env
      TYPE(section_vals_type), POINTER :: root, fe, sub, opts, dft
      TYPE(mp_para_env_type), POINTER :: para
      LOGICAL, INTENT(IN) :: motion
      TYPE(distribution_1d_type), POINTER :: particles, molecules
      INTEGER, DIMENSION(:), POINTER :: repeats
      REAL(KIND=dp), DIMENSION(:), POINTER :: radii
      CHARACTER(LEN=default_path_length) :: library
      CHARACTER(LEN=512, KIND=c_char) :: msg
      INTEGER :: axis, n, term, unit_nr
      LOGICAL :: explicit

      dft => section_vals_get_subs_vals(fe, "DFT")
      CALL section_vals_get(dft, explicit=explicit)
      IF (explicit) CPABORT("Use &GXTB to configure METHOD GXTB. &DFT settings are unsupported.")
      opts => section_vals_get_subs_vals(fe, "GXTB")
      CALL section_vals_val_get(opts, "LIBRARY", c_val=library)
      CALL section_vals_val_get(opts, "REPEATS", i_vals=repeats)
      CALL section_vals_val_get(opts, "EXCHANGE_RADII", r_vals=radii)
      IF (ANY(repeats < 1)) CPABORT("GXTB REPEATS must be positive")
      env%repeats = repeats
      env%radii = radii
      CALL cp_subsys_create(env%subsys, para, root, fe, sub, motion)
      IF (ANY(env%subsys%cell%perd /= 1)) CPABORT("METHOD GXTB requires PERIODIC XYZ")
      n = env%subsys%particles%n_els
      IF (n < 1) CPABORT("METHOD GXTB requires atoms")
      DO axis = 1, 3
         IF (n > HUGE(0_c_int)/repeats(axis)) THEN
            CPABORT("GXTB supercell exceeds the provider atom count range")
         END IF
         n = n*repeats(axis)
      END DO
      env%copies = n/env%subsys%particles%n_els
      NULLIFY (particles, molecules)
      CALL distribute_molecules_1d(env%subsys%atomic_kinds%els, env%subsys%particles%els, particles, &
                                   env%subsys%molecule_kinds%els, env%subsys%molecules%els, molecules, fe)
      CALL cp_subsys_set(env%subsys, local_particles=particles, local_molecules=molecules)
      CALL distribution_1d_release(particles)
      CALL distribution_1d_release(molecules)
      IF (para%is_source()) THEN
         env%model = cp_gxtb_open(TRIM(library)//c_null_char, msg, LEN(msg, KIND=c_size_t))
         term = INDEX(msg, c_null_char)
         IF (term > 0) msg(term:) = " "
         IF (.NOT. c_associated(env%model)) CPABORT("Cannot open GXTB provider: "//TRIM(msg))
         unit_nr = cp_logger_get_default_io_unit()
         IF (unit_nr > 0) WRITE (unit_nr, '(1X,A,A,A,3(1X,I0))') "GXTB| ", TRIM(msg), "   Supercell:", repeats
      END IF
   END SUBROUTINE

   SUBROUTINE gxtb_evaluate(env, forces, virial)
      TYPE(gxtb_env_type), INTENT(INOUT) :: env
      LOGICAL, INTENT(IN) :: forces, virial
      REAL(KIND=dp), ALLOCATABLE :: pos(:,:), grad(:,:), force(:,:)
      INTEGER(c_int), ALLOCATABLE :: z(:)
      INTEGER :: n, total, a, b, c, i, j, last, stat, term
      REAL(KIND=dp) :: h(3,3), strain(3,3), shift(3)
      CHARACTER(LEN=512, KIND=c_char) :: msg
      TYPE(cp_subsys_type), POINTER :: sys

      sys => env%subsys
      n = sys%particles%n_els
      ALLOCATE (force(3,n))
      IF (sys%para_env%is_source()) THEN
         total = n*env%copies
         ALLOCATE (pos(3,total), grad(3,total), z(total))
         h = sys%cell%hmat
         DO i = 1, n
            CALL get_atomic_kind(sys%particles%els(i)%atomic_kind, z=z(i))
            IF (z(i) < 1 .OR. z(i) > 103) CPABORT("Unsupported atomic number in GXTB")
         END DO
         ! Each consecutive block contains the full basis in input atom order.
         j = 0
         DO a = 0, env%repeats(1)-1
            DO b = 0, env%repeats(2)-1
               DO c = 0, env%repeats(3)-1
                  shift = a*h(:,1) + b*h(:,2) + c*h(:,3)
                  DO i = 1, n
                     pos(:,j+i) = sys%particles%els(i)%r + shift
                  END DO
                  IF (j > 0) z(j+1:j+n) = z(:n)
                  j = j + n
               END DO
            END DO
         END DO
         DO i = 1, 3
            h(:,i) = h(:,i)*env%repeats(i)
         END DO
         stat = cp_gxtb_eval(env%model, INT(total,c_int), z, pos, h, env%radii, MERGE(1_c_int,0_c_int,forces), &
                             MERGE(1_c_int,0_c_int,virial), env%energy, grad, strain, &
                             msg, LEN(msg, KIND=c_size_t))
         IF (stat /= 0) THEN
            term = INDEX(msg, c_null_char)
            IF (term > 0) msg(term:) = " "
            CPABORT("GXTB evaluation failed: "//TRIM(msg))
         END IF
         IF (.NOT. (ieee_is_finite(env%energy) .AND. ALL(ieee_is_finite(grad)) .AND. ALL(ieee_is_finite(strain)))) THEN
            CPABORT("GXTB returned non-finite energy or derivatives")
         END IF
         ! Moving a primitive atom moves every image. Differentiate the energy
         ! per primitive cell; CP2K stores pressure-positive virials.
         force = 0.0_dp
         DO last = n, total, n
            force = force - grad(:,last-n+1:last)
         END DO
         force = force/env%copies
         env%energy = env%energy/env%copies
         strain = -strain/env%copies
      END IF
      CALL sys%para_env%bcast(env%energy)
      CALL sys%para_env%bcast(force)
      CALL sys%para_env%bcast(strain)
      DO i = 1, n
         sys%particles%els(i)%f = force(:,i)
      END DO
      IF (sys%virial%pv_availability .AND. .NOT. sys%virial%pv_numer) sys%virial%pv_virial = strain
   END SUBROUTINE
END MODULE
