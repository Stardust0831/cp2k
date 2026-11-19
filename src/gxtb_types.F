! SPDX-License-Identifier: GPL-2.0-or-later
MODULE gxtb_types
   USE cp_subsys_types, ONLY: cp_subsys_type, cp_subsys_release
   USE, INTRINSIC :: iso_c_binding, ONLY: c_associated, c_null_ptr, c_ptr
   USE kinds, ONLY: dp
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: gxtb_env_type, gxtb_release
   TYPE gxtb_env_type
      TYPE(cp_subsys_type), POINTER :: subsys => NULL()
      TYPE(c_ptr) :: model = c_null_ptr
      INTEGER :: repeats(3) = 1, copies = 1
      REAL(KIND=dp) :: energy = 0.0_dp, radii(2) = 0.0_dp
   END TYPE
CONTAINS
   SUBROUTINE gxtb_release(env)
      TYPE(gxtb_env_type), INTENT(INOUT) :: env
      INTERFACE
         SUBROUTINE cp_gxtb_close(model) BIND(C)
            IMPORT :: c_ptr
            TYPE(c_ptr), VALUE :: model
         END SUBROUTINE
      END INTERFACE
      IF (c_associated(env%model)) CALL cp_gxtb_close(env%model)
      CALL cp_subsys_release(env%subsys)
   END SUBROUTINE
END MODULE
