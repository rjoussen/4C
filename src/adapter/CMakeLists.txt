# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    ale
    beam3
    beaminteraction
    config
    constraint
    contact
    core
    coupling
    fluid
    fluid_ele
    fluid_turbulence
    fluid_xfluid
    fpsi
    fsi
    global_data
    inpar
    levelset
    lubrication
    mat
    membrane
    mortar
    pasi
    poroelast
    porofluid_pressure_based
    scatra
    shell7p
    solid_3D_ele
    solver_nonlin_nox
    ssi
    structure
    structure_new
    thermo
    truss3
    w1
    xfem
    )

four_c_add_internal_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})

four_c_unity_build_compile_separately(${AUTO_DEFINED_MODULE_NAME} 4C_adapter_str_fsiwrapper.cpp)
four_c_unity_build_compile_separately(
  ${AUTO_DEFINED_MODULE_NAME} 4C_adapter_str_timint_adaptive.cpp
  )
four_c_unity_build_compile_separately(
  ${AUTO_DEFINED_MODULE_NAME} 4C_adapter_structure_scatra_ele.cpp
  )
