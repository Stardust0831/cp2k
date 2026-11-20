cmake_minimum_required(VERSION 3.22)
project(cp2k_open_gxtb LANGUAGES C CXX Fortran)

# Link one set of Fortran modules privately into the adapter. CP2K can
# already have its own tblite and multicharge with a different Fortran ABI.
set(BUILD_SHARED_LIBS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(WITH_TESTS OFF CACHE BOOL "")
set(BUILD_TESTING OFF CACHE BOOL "")
set(WITH_JSON OFF CACHE BOOL "")
set(WITH_DDX OFF CACHE BOOL "")
set(WITH_TREXIO OFF CACHE BOOL "")
set(WITH_API ON CACHE BOOL "")
set(WITH_OpenMP ON CACHE BOOL "")
set(SDFTD3_WITH_TESTS OFF CACHE BOOL "")
set(TBLITE_WITH_DDX OFF CACHE BOOL "")
find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
function(add_dependency name location)
  add_subdirectory("${location}" "${name}")
  if(NOT TARGET "${name}::${name}")
    add_library("${name}::${name}" ALIAS "${name}")
  endif()
endfunction()
add_dependency(toml-f "${STOCK_SOURCE}/subprojects/toml-f")
add_dependency(mctc-lib "${STOCK_SOURCE}/subprojects/mctc-lib")
add_dependency(mstore "${STOCK_SOURCE}/subprojects/mstore")
add_dependency(multicharge "${MULTICHARGE_SOURCE}")
add_dependency(dftd4 "${STOCK_SOURCE}/subprojects/dftd4")
add_dependency(s-dftd3 "${STOCK_SOURCE}/subprojects/s-dftd3")
add_subdirectory("${PROVIDER_SOURCE}" save_tblite)
target_compile_options(tblite-lib PRIVATE
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffree-line-length-none>)

target_sources(tblite-exe PRIVATE host_d4.c ../../../src/gxtb_revd4.cpp)
target_compile_features(tblite-exe PRIVATE cxx_std_17)
target_compile_definitions(tblite-exe PRIVATE GXTB_STANDALONE)
add_library(gxtb_adapter OBJECT ../save_tblite.c host_d4.c)
target_link_libraries(gxtb_adapter PRIVATE tblite-lib)
add_library(gxtb_save SHARED $<TARGET_OBJECTS:gxtb_adapter>)
target_link_libraries(gxtb_save PRIVATE tblite-lib)
target_link_options(gxtb_save PRIVATE "LINKER:--exclude-libs,ALL"
  "LINKER:--version-script,${CMAKE_CURRENT_SOURCE_DIR}/provider.exports")
set_property(TARGET gxtb_save APPEND PROPERTY LINK_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/provider.exports")
set_target_properties(gxtb_save PROPERTIES INSTALL_RPATH "$ORIGIN")
set_target_properties(tblite-exe PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
# CP2K can also link stock tblite. Bind the private operator-module symbols
# within one relocatable object before putting it in the static archive.
# BLAS, OpenMP and compiler-runtime symbols resolve against the CP2K build.
set(kernel_archives tblite-lib toml-f-lib dftd4-lib multicharge-lib s-dftd3-lib mctc-lib-lib)
set(kernel_files)
foreach(archive IN LISTS kernel_archives)
  list(APPEND kernel_files "$<TARGET_FILE:${archive}>")
endforeach()
set(kernel_object "${CMAKE_CURRENT_BINARY_DIR}/gxtb_kernels.o")
add_custom_command(OUTPUT "${kernel_object}"
  COMMAND "${CMAKE_LINKER}" -r -o gxtb_kernels.raw.o $<TARGET_OBJECTS:gxtb_adapter>
          --start-group ${kernel_files} --end-group
  COMMAND "${CMAKE_OBJCOPY}"
          --keep-global-symbol=cp2k_gxtb_name_v1 --keep-global-symbol=cp2k_gxtb_eval_v2
          --keep-global-symbol=cp2k_gxtb_eval_v3 gxtb_kernels.raw.o "${kernel_object}"
  DEPENDS gxtb_adapter $<TARGET_OBJECTS:gxtb_adapter> ${kernel_archives}
  BYPRODUCTS gxtb_kernels.raw.o COMMAND_EXPAND_LISTS VERBATIM)
add_library(gxtb_kernels STATIC "${kernel_object}")
set_target_properties(gxtb_kernels PROPERTIES LINKER_LANGUAGE Fortran)
install(TARGETS gxtb_save gxtb_kernels tblite-exe
  ARCHIVE DESTINATION lib COMPONENT gxtb_runtime
  LIBRARY DESTINATION lib COMPONENT gxtb_runtime
  RUNTIME DESTINATION bin COMPONENT gxtb_runtime)
