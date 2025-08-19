# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    ale
    config
    constraint
    core
    fbi
    fluid
    fluid_xfluid
    fsi_xfem
    global_data
    inpar
    mat
    mortar
    poroelast
    solver_nonlin_nox
    structure
    structure_new
    timestepping
    )

four_c_add_internal_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
