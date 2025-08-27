// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_GENERAL_UTILS_GAUSS_POINT_POSTPROCESS_HPP
#define FOUR_C_FEM_GENERAL_UTILS_GAUSS_POINT_POSTPROCESS_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_vector.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Elements
{
  class Element;
}

namespace Core::FE
{
  /*!
   * \brief Assemble nodal element count
   *
   * \param global_count Add a 1 to all nodes belonging to this element
   * \param ele element
   */
  void assemble_nodal_element_count(
      Core::LinAlg::Vector<int>& global_count, const Core::Elements::Element& ele);


  /*!
   * \brief Assemble Gauss point data into an array of global cell data
   *
   * \param global_data array of global cell data (length at least number of gauss points)
   * \param gp_data (numgp x size) matrix of the Gauss point data
   * \param ele element
   */
  void assemble_gauss_point_values(
      std::vector<std::shared_ptr<Core::LinAlg::MultiVector<double>>>& global_data,
      const Core::LinAlg::SerialDenseMatrix& gp_data, const Core::Elements::Element& ele);


  /*!
   * @brief Extrapolation of Gauss point quantities given in @data to the nodes of the Element @ele
   * using the shape functions of the element and assembly to the global nodal data @nodal_data.
   *
   * @note On shared nodes, the values of all participating elements will be averaged
   *
   * @param ele (in) : Reference to finite element
   * @param data (in) : Gauss point data in a Matrix (numgp x numdim of vector)
   * @param dis (in) : Reference to the discretization
   * @param nodal_data (out) : Assembled data
   */
  void extrapolate_gauss_point_quantity_to_nodes(Core::Elements::Element& ele,
      const Core::LinAlg::SerialDenseMatrix& data, const Core::FE::Discretization& dis,
      Core::LinAlg::MultiVector<double>& nodal_data);

  /*!
   * @brief Averaging of all Gauss point quantities in @data within the element @ele and assembly to
   * the element vector @element_data
   *
   * @param ele (in) : Reference to finite element
   * @param data (in) : Gauss point data in a Matrix (numgp x numdim of vector)
   * @param element_data (out) : Assembled data
   */
  void evaluate_gauss_point_quantity_at_element_center(Core::Elements::Element& ele,
      const Core::LinAlg::SerialDenseMatrix& data, Core::LinAlg::MultiVector<double>& element_data);

  /*!
   * \brief Assemble averaged data. The data at the Gauss points are averaged within the element.
   *
   * \tparam T Type of the data, either SerialDenseMatrix or Core::LinAlg::Matrix
   * \param global_data Global cell data
   * \param gp_data (numgp x size) matrix of the Gauss point data
   * \param ele element
   */
  template <class T>
  void assemble_averaged_element_values(Core::LinAlg::MultiVector<double>& global_data,
      const T& gp_data, const Core::Elements::Element& ele);


  // --- template and inline functions --- //
  template <class T>
  void assemble_averaged_element_values(Core::LinAlg::MultiVector<double>& global_data,
      const T& gp_data, const Core::Elements::Element& ele)
  {
    const Core::LinAlg::Map& elemap = global_data.get_map();

    int lid = elemap.lid(ele.id());
    if (lid != -1)
    {
      for (decltype(gp_data.numCols()) i = 0; i < gp_data.numCols(); ++i)
      {
        double& s = (global_data(i)).get_values()[lid];  // resolve pointer for faster access
        s = 0.;
        for (decltype(gp_data.numRows()) j = 0; j < gp_data.numRows(); ++j)
        {
          s += gp_data(j, i);
        }
        s *= 1.0 / gp_data.numRows();
      }
    }
  }

}  // namespace Core::FE

FOUR_C_NAMESPACE_CLOSE

#endif