// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_bele_bele3.hpp"
#include "4C_fem_condition.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_densematrix_multiply.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN

int Discret::Elements::Bele3Line::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  return 0;
}

int Discret::Elements::Bele3Line::evaluate_neumann(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, const Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  return 0;
}

FOUR_C_NAMESPACE_CLOSE
