// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MORTAR_UTILS_HPP
#define FOUR_C_MORTAR_UTILS_HPP

#include "4C_config.hpp"

#include "4C_linalg_map.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_mortar_coupling3d_classes.hpp"

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

  /*!
  \brief Sort vector in ascending order

  This routine is taken from Trilinos MOERTEL package.

  \param dlist (in): vector to be sorted (unsorted on input, sorted on output)
  \param N (in):     length of vector to be sorted
  \param list2 (in): another vector which is sorted accordingly
  */
  void sort(double* dlist, int N, int* list2);

  /*!
  \brief Transform the row map of a matrix (only GIDs)

  This method changes the row map of an input matrix to a new row map
  with different GID numbering. However, the parallel distribution of
  this new row map is exactly the same as in the old row map. Thus, this
  is simply a processor-local 1:1 matching of old and new GIDs.

  @param inmat Matrix on which the row and column maps will be transformed
  @param newrowmap New row map used for the given input matrix

  \post Output matrix will be fill_complete()
  */
  std::shared_ptr<Core::LinAlg::SparseMatrix> matrix_row_transform_gids(
      const Core::LinAlg::SparseMatrix& inmat, const Core::LinAlg::Map& newrowmap);

  /*!
  \brief Transform the column map of a matrix (only GIDs)

  This method changes the column map of an input matrix to a new column
  map with different GID numbering (and the domain map, accordingly).
  However, the parallel distribution of the new domain map is exactly
  the same as in the old domain map. Thus, this is simply a processor-local
  1:1 matching of old and new GIDs.

  @param inmat Matrix on which the row and column maps will be transformed
  @param newdomainmap New domain map used for the given input matrix, which will indirectly
  transform the column map of the given input matrix

  \post Output matrix will be fill_complete()
  */
  std::shared_ptr<Core::LinAlg::SparseMatrix> matrix_col_transform_gids(
      const Core::LinAlg::SparseMatrix& inmat, const Core::LinAlg::Map& newdomainmap);


  /*!
  \brief Transform the row and column maps of a matrix (only GIDs)

  This method changes the row and column maps of an input matrix to new
  row and column maps with different GID numbering (and the domain map,
  accordingly). However, the parallel distribution of the new row and
  domain maps is exactly the same as in the old ones. Thus, this is simply
  a processor-local 1:1 matching of old and new GIDs.

  @param inmat Matrix on which the row and column maps will be transformed
  @param newrowmap New row map used for the given input matrix
  @param newdomainmap New domain map used for the given input matrix, which will indirectly
  transform the column map of the given input matrix

  \post Output matrix will be fill_complete()
  */
  std::shared_ptr<Core::LinAlg::SparseMatrix> matrix_row_col_transform_gids(
      const Core::LinAlg::SparseMatrix& inmat, const Core::LinAlg::Map& newrowmap,
      const Core::LinAlg::Map& newdomainmap);

  /*!
  \brief Convex hull points are sorted in order to obtain final clip polygon

  \param out (in): bool to switch output on/off
  \param transformed (in): coordinates of vertex objects transformed into auxiliary plane
  \param collconvexhull (in): vector of vertex objects to be sorted
  \param respoly (out): vector of vertex objects for result polygon
  \param tol (in): clipping tolerance for close vertices detection
  \return number of removed points from collconvexhull

  */
  int sort_convex_hull_points(bool out, Core::LinAlg::SerialDenseMatrix& transformed,
      std::vector<Vertex>& collconvexhull, std::vector<Vertex>& respoly, double& tol);

  namespace Utils
  {
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
