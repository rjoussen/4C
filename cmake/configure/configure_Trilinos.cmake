# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# Kokkos is typically pulled in via Trilinos. If no location has been given,
# try the same location as Trilinos. If no Trilinos location exists, users
# will get an error to provide that one first.
set(Kokkos_FIND_QUIETLY TRUE)
if(Trilinos_ROOT AND NOT Kokkos_ROOT)
  set(Kokkos_ROOT
      ${Trilinos_ROOT}
      CACHE PATH "Path to Kokkos installation"
      )
endif()

# We only support Trilinos versions that provide a config file.
find_package(Trilinos REQUIRED)

message(STATUS "Trilinos version: ${Trilinos_VERSION}")
message(STATUS "Trilinos packages: ${Trilinos_PACKAGE_LIST}")

# Figure out the version.
if(EXISTS "${Trilinos_DIR}/../../../TrilinosRepoVersion.txt")
  file(STRINGS "${Trilinos_DIR}/../../../TrilinosRepoVersion.txt" TrilinosRepoVersionFile)
  # The hash is the first token on the second line.
  list(GET TrilinosRepoVersionFile 1 TrilinosRepoVersionFileLine2)
  separate_arguments(TrilinosRepoVersionFileLine2)
  list(GET TrilinosRepoVersionFileLine2 0 _sha)

  set(FOUR_C_Trilinos_GIT_HASH ${_sha})
else()
  set(FOUR_C_Trilinos_GIT_HASH "unknown")
endif()

target_link_libraries(
  four_c_all_enabled_external_dependencies INTERFACE Trilinos::all_selected_libs
  )

configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/templates/Trilinos.cmake.in
  ${PROJECT_BINARY_DIR}/cmake/templates/Trilinos.cmake
  @ONLY
  )
