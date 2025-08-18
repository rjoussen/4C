// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MORTAR_MULTIFIELD_COUPLING_HPP
#define FOUR_C_MORTAR_MULTIFIELD_COUPLING_HPP

#include "4C_config.hpp"

#include "4C_fem_general_shape_function_type.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_function_manager.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Core::LinAlg
{
  class SparseMatrix;
  class BlockSparseMatrixBase;
}  // namespace Core::LinAlg

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Mortar
{
  class MultiFieldCoupling
  {
   public:
    /// c-tor
    MultiFieldCoupling() {};


    /// add a new discretization to perform coupling on
    void push_back_coupling(
        const std::shared_ptr<Core::FE::Discretization>& dis,  ///< discretization
        const int nodeset,                                     ///< nodeset to couple
        const std::vector<int>& dofs_to_couple,                ///< dofs to couple
        const Teuchos::ParameterList& mortar_params, const Teuchos::ParameterList& contact_params,
        const Teuchos::ParameterList& binning_params,
        const std::map<std::string, std::shared_ptr<Core::FE::Discretization>>& discretization_map,
        const Core::Utils::FunctionManager& function_manager,
        std::shared_ptr<Core::IO::OutputControl> output_control,
        Core::FE::ShapeFunctionType shape_function_type, int ndim);

    /// Perform condensation in all blocks of the matrix
    void condense_matrix(std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase>& mat);

    /// Perform condensation in the right-hand side
    void condense_rhs(Core::LinAlg::Vector<double>& rhs);

    /// recover condensed primal slave-sided dofs
    void recover_incr(Core::LinAlg::Vector<double>& incr);

   private:
    std::vector<std::shared_ptr<Core::LinAlg::SparseMatrix>> p_;
  };
}  // namespace Mortar



FOUR_C_NAMESPACE_CLOSE

#endif
