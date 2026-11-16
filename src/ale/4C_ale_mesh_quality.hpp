// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ALE_MESH_QUALITY_HPP
#define FOUR_C_ALE_MESH_QUALITY_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <array>
#include <span>

FOUR_C_NAMESPACE_OPEN

namespace ALE
{
  /*!
   * @brief Parameters of the mesh quality potential
   */
  struct MeshQualityParameters
  {
    //! weight of the size term relative to the shape term
    double size_weight = 1.0;

    //! weight of the motion term relative to the shape term
    double motion_weight = 0.0;

    //! whether the size term penalizes only volume ratios above one (growing elements)
    bool size_growth_only = false;

    //! exponent p >= 1 of the shape term (mu^p - 1) / p of the mean ratio measure mu >= 1, which
    //! penalizes strongly distorted elements more than the plain shape term mu - 1 (p = 1)
    double shape_exponent = 1.0;
  };

  /*!
   * @brief Mean ratio quality of an element
   *
   * The quality is evaluated at all corners of the element based on the Jacobian of the mapping
   * from parameter space, compared to the Jacobian of an ideal element (cube or regular
   * tetrahedron). It is 1 for an ideal element shape, tends to 0 for degenerated elements and is
   * negative for inverted elements. The minimum over the corners is returned.
   *
   * @param celltype Cell type (hex8 or tet4)
   * @param coordinates Nodal coordinates (x, y, z of each node)
   */
  double mean_ratio_quality(Core::FE::CellType celltype, std::span<const double> coordinates);

  /*!
   * @brief Evaluate the mesh quality potential of an element
   *
   * The potential penalizes the distortion of the element shape and changes of the element size.
   * At each corner of the element with Jacobian J, it contains the mean ratio distortion
   * mu - 1 of mu = |J W^-1|^2 / (3 det(J W^-1)^(2/3)) of the Jacobian relative to a target Jacobian
   * W, or (mu^p - 1) / p with an exponent p > 1 (see MeshQualityParameters), and the
   * size term v + 1/v - 2 with the volume ratio v = det(J) / det(W) (or relative to the target
   * volume, see MeshQualityParameters), or, for growing elements only, (v - 1)^3 / v^2 for v > 1
   * and zero otherwise, whose Hessian is continuous at v = 1. Both terms are zero in the optimum;
   * the shape term and the two-sided size term act as barriers against inversion. The motion term
   * is the mean over the nodes of |U W^-1 u|^2, the squared displacement u of a node relative to
   * the element size, with the Jacobian W at its corner (of the reference configuration if there is
   * no target) and the Jacobian U of the ideal element with unit edges; for an ideal element of
   * size h, it is (|u| / h)^2. Since the displacement is measured from the reference configuration
   * of the mesh movement, the term keeps the nodes close to it, e.g., close to the material if the
   * reference is the current configuration of a solid.
   *
   * With target coordinates, W is the Jacobian of the target configuration, such that the target
   * configuration is optimal (e.g., the initial mesh with its grading and anisotropy). Without, the
   * shape of W is that of an ideal element (cube or regular tetrahedron) and its size the one in
   * the reference configuration.
   *
   * @param celltype Cell type (hex8 or tet4)
   * @param reference_coordinates Nodal coordinates of the reference configuration of the mesh
   * movement
   * @param displacements Nodal displacements of the mesh movement
   * @param parameters Parameters of the potential
   * @param residual Gradient of the potential with respect to the displacements (if not null)
   * @param stiffness Hessian of the potential with respect to the displacements (if not null)
   * @param target_coordinates Nodal coordinates of the target configuration (optional)
   * @return The value of the potential, or NaN if the element or its target is inverted
   */
  double evaluate_mesh_quality_potential(Core::FE::CellType celltype,
      std::span<const double> reference_coordinates, std::span<const double> displacements,
      const MeshQualityParameters& parameters, Core::LinAlg::SerialDenseVector* residual,
      Core::LinAlg::SerialDenseMatrix* stiffness, std::span<const double> target_coordinates = {});

  /*!
   * @brief Size of an element
   *
   * The mean determinant of the Jacobians at the corners times the volume of the parameter space,
   * which is the volume of elements with a constant Jacobian.
   *
   * @param celltype Cell type (hex8 or tet4)
   * @param coordinates Nodal coordinates (x, y, z of each node)
   */
  double element_size(Core::FE::CellType celltype, std::span<const double> coordinates);
}  // namespace ALE

FOUR_C_NAMESPACE_CLOSE

#endif
