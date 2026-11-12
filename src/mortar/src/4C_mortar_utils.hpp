// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MORTAR_UTILS_HPP
#define FOUR_C_MORTAR_UTILS_HPP

#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_mortar_coupling3d_classes.hpp"
#include "4C_mortar_element.hpp"

#include <array>
#include <memory>

FOUR_C_NAMESPACE_OPEN


// forward declarations
namespace Core::LinAlg
{
  class SparseMatrix;
  class BlockSparseMatrixBase;
}  // namespace Core::LinAlg

namespace Mortar
{
  namespace Utils
  {
    /*!
    \brief Orthonormal coordinate system in a plane

    The origin is the first vertex of a polygon and the x-axis points to its second vertex.
    */
    class PlaneCoordinateSystem
    {
     public:
      PlaneCoordinateSystem(const std::vector<Vertex>& polygon, const double* normal);

      //! coordinates of the vertices of a polygon in the plane
      std::vector<std::array<double, 2>> to_plane(const std::vector<Vertex>& polygon) const;

      //! coordinates of a point of the plane in space
      std::vector<double> to_space(const std::array<double, 2>& point) const;

     private:
      std::array<double, 3> origin_;
      std::array<double, 3> x_axis_;
      std::array<double, 3> y_axis_;
    };

    /*!
    \brief Clip two convex polygons, which lie in a plane

    \param source_polygon (in): vertices of the source polygon
    \param target_polygon (in): vertices of the target polygon
    \param normal (in): normal of the plane
    \param tolerance (in): distance below which vertices are treated as identical
    \return vertices of the clip polygon (empty if one of the polygons is not convex or if the
    clip polygon has no area)
    */
    std::vector<Vertex> clip_polygons(const std::vector<Vertex>& source_polygon,
        const std::vector<Vertex>& target_polygon, const double* normal, double tolerance);

    /*!
    \brief copy the ghosting of dis_src to all discretizations with names in
           vector voldis. Material pointers can be added according to
           link_materials
    */
    void create_volume_ghosting(const Core::FE::Discretization& dis_src,
        const std::vector<std::shared_ptr<Core::FE::Discretization>>& voldis,
        std::vector<std::pair<int, int>> material_links, bool check_on_in = true,
        bool check_on_exit = true);


    /*!
    \brief Prepare mortar element for nurbs case


    store knot vector, zerosized information and normal factor
    */
    void prepare_nurbs_element(Core::FE::Discretization& discret,
        std::shared_ptr<Core::Elements::Element> ele, Mortar::Element& cele, int dim);

    /*!
    \brief Prepare mortar node for nurbs case

    store control point weight

    */
    void prepare_nurbs_node(Core::Nodes::Node* node, Mortar::Node& mnode);

    void mortar_matrix_condensation(std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase>& k,
        const std::vector<std::shared_ptr<Core::LinAlg::SparseMatrix>>& p);

    /*! \brief Perform static condensation of Jacobian with mortar matrix \f$D^{-1}M\f$
     *
     * @param[in/out] k Matrix to be condensed
     * @param[in] p_row Mortar projection operator for condensation of rows
     * @param[in] p_col Mortar projection operator for condenstaion of columns
     */
    void mortar_matrix_condensation(std::shared_ptr<Core::LinAlg::SparseMatrix>& k,
        const std::shared_ptr<const Core::LinAlg::SparseMatrix>& p_row,
        const std::shared_ptr<const Core::LinAlg::SparseMatrix>& p_col);

    void mortar_rhs_condensation(Core::LinAlg::Vector<double>& rhs, Core::LinAlg::SparseMatrix& p);

    void mortar_rhs_condensation(Core::LinAlg::Vector<double>& rhs,
        const std::vector<std::shared_ptr<Core::LinAlg::SparseMatrix>>& p);

    void mortar_recover(Core::LinAlg::Vector<double>& inc, Core::LinAlg::SparseMatrix& p);

    void mortar_recover(Core::LinAlg::Vector<double>& inc,
        const std::vector<std::shared_ptr<Core::LinAlg::SparseMatrix>>& p);
  }  // namespace Utils
}  // namespace Mortar

FOUR_C_NAMESPACE_CLOSE

#endif
