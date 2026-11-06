#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Find the optional cuDSS sparse direct solver.
find_path(
  CUDSS_INCLUDE_DIR cudss.h
  HINTS ${CP2K_CUDSS_ROOT} ENV CUDSS_ROOT
  PATH_SUFFIXES include)
find_library(
  CUDSS_LIBRARY
  NAMES cudss libcudss.so.1
  HINTS ${CP2K_CUDSS_ROOT} ENV CUDSS_ROOT
  PATH_SUFFIXES lib lib64)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CuDSS REQUIRED_VARS CUDSS_LIBRARY
                                                      CUDSS_INCLUDE_DIR)
if(CuDSS_FOUND AND NOT TARGET cp2k::CuDSS)
  add_library(cp2k::CuDSS UNKNOWN IMPORTED)
  set_target_properties(
    cp2k::CuDSS PROPERTIES IMPORTED_LOCATION "${CUDSS_LIBRARY}"
                           INTERFACE_INCLUDE_DIRECTORIES "${CUDSS_INCLUDE_DIR}")
endif()
mark_as_advanced(CUDSS_INCLUDE_DIR CUDSS_LIBRARY)
