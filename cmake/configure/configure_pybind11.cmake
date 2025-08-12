# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

if(NOT FOUR_C_WITH_PYTHON)
  message(FATAL_ERROR "FOUR_C_WITH_PYBIND11 requires FOUR_C_WITH_PYTHON to be enabled.")
endif()

message(STATUS "Fetch content for pybind11")
fetchcontent_declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG ed5057ded698e305210269dafa57574ecf964483 # version 3.0.0
  )
set(PYBIND11_FINDPYTHON NEW)
set(PYBIND11_INSTALL ON)
fetchcontent_makeavailable(pybind11)

set(FOUR_C_PYBIND11_ROOT "${CMAKE_INSTALL_PREFIX}")

four_c_add_external_dependency(
  four_c_all_enabled_external_dependencies pybind11::module pybind11::embed
  )

configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/templates/pybind11.cmake.in
  ${PROJECT_BINARY_DIR}/cmake/templates/pybind11.cmake
  @ONLY
  )
