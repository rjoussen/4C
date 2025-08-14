// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_NLN_LINEARSYSTEM_SCALING_HPP
#define FOUR_C_STRUCTURE_NEW_NLN_LINEARSYSTEM_SCALING_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"
#include "4C_linalg_map.hpp"
#include "4C_solver_nonlin_nox_scaling.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Solid
{
  namespace TimeInt
  {
    class BaseDataSDyn;
    class BaseDataGlobalState;
  }  // namespace TimeInt
}  // namespace Solid

namespace Core::LinAlg
{
  class SparseMatrix;
}

namespace NOX::Nln
{
  class LinearProblem;
}

namespace Solid
{
  namespace Nln
  {
    namespace LinSystem
    {
      class StcScaling : public NOX::Nln::Scaling
      {
       public:
        //! Constructor.
        StcScaling(const Solid::TimeInt::BaseDataSDyn& DataSDyn,
            Solid::TimeInt::BaseDataGlobalState& GState);

        //! Compute scaling.
        void compute_scaling(const NOX::Nln::LinearProblem& problem) override;

        //! Scales the linear system.
        void scale_linear_system(NOX::Nln::LinearProblem& problem) override;

        //! Remove the scaling from the linear system.
        void unscale_linear_system(NOX::Nln::LinearProblem& problem) override;

       private:
        //! stiffness matrix after scaling
        std::shared_ptr<Core::LinAlg::SparseMatrix> stiff_scaled_;

        //! scale thickness of shells
        const enum Inpar::Solid::StcScale stcscale_;

        //! number of layers for multilayered case
        const int stclayer_;

        //! scaling matrix for STC
        std::shared_ptr<Core::LinAlg::SparseMatrix> stcmat_;
      };
    }  // namespace LinSystem
  }  // namespace Nln
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif