// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLVER_NONLIN_NOX_AUX_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_AUX_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_enum_lists.hpp"
#include "4C_solver_nonlin_nox_forward_decl.hpp"
#include "4C_solver_nonlin_nox_statustest_factory.hpp"

#include <mpi.h>
#include <NOX_Abstract_Vector.H>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Core::LinAlg
{
  template <typename T>
  class Vector;
  class Solver;
  class SparseOperator;
}  // namespace Core::LinAlg

namespace NOX
{
  namespace Nln
  {
    namespace Aux
    {
      /*! Set printing parameters
       *
       *  Note: The Yes/No tuples are translated to booleans! */
      void set_printing_parameters(Teuchos::ParameterList& p_nox, MPI_Comm comm);

      /*! \brief Returns the type of operator that is passed in.
       *
       *   Uses dynamic casting to identify the underlying object type. */
      NOX::Nln::LinSystem::OperatorType get_operator_type(const Core::LinAlg::SparseOperator& op);

      /// return linear system type
      NOX::Nln::LinSystem::LinearSystemType get_linear_system_type(
          const std::map<enum NOX::Nln::SolutionType, Teuchos::RCP<Core::LinAlg::Solver>>&
              linsolvers);

      /*! \brief Calculate the root mean square for the NOX status test
       *  \f[
       *    \delta_{rms} = \sqrt{\frac{1}{N} \sum\limits_{i=1}^{N} \left( \frac{x_{i}^{k} -
       * x_{i}^{k-1}}{\mathcal{RTOL} | x_{i}^{k-1} | + \mathcal{ATOL}} \right)} \f]
       *
       *  \param atol  : absolute tolerance
       *  \param rtol  : relative tolerance
       *  \param xnew  : new / current iterate $x_{i}^{k}$
       *  \param xincr : current step increment $x_{i}^{k} - x_{i}^{k-1}$
       */
      double root_mean_square_norm(const double& atol, const double& rtol,
          const Core::LinAlg::Vector<double>& xnew, const Core::LinAlg::Vector<double>& xincr,
          const bool& disable_implicit_weighting = false);

      /*! \brief Do a recursive search for a NOX::Nln::StatusTest::NormWRMS object in the StatusTest
       * object list and return the class variable value of the desired quantity.
       *
       * \param test              : StatusTest object which will be scanned.
       * \param qType             : Quantity type of the NormWRMS test which we are looking for.
       * \param classVariableName : Name of the class variable which will be returned. (Type:
       * double) */
      double get_norm_wrms_class_variable(const ::NOX::StatusTest::Generic& test,
          const NOX::Nln::StatusTest::QuantityType& qType, const std::string& classVariableName);

      /*! \brief Do a recursive search for a NOX::Nln::StatusTest::NormF object in the StatusTest
       * object list and return the class variable value of the desired quantity.
       *
       * \param test              : StatusTest object which will be scanned.
       * \param qType             : Quantity type of the NormF test which we are looking for.
       * \param classVariableName : Name of the class variable which will be returned. (Type:
       * double) */
      double get_norm_f_class_variable(const ::NOX::StatusTest::Generic& test,
          const NOX::Nln::StatusTest::QuantityType& qType, const std::string& classVariableName);

      /*! Do a recursive search for a <T> status test and the given quantity
       *
       *  True is returned as soon as a status test of type <T> is found, which
       *  holds the given quantity. */
      template <class T>
      bool is_quantity(
          const ::NOX::StatusTest::Generic& test, const NOX::Nln::StatusTest::QuantityType& qtype);

      /*! \brief Do a recursive search for a <T> status test class and return the NormType of the
       * given quantity.
       *
       *  If there are more than one status tests of the type <T> which hold the given quantity, the
       * normtype of the first we can find, will be returned! */
      template <class T>
      int get_norm_type(
          const ::NOX::StatusTest::Generic& test, const NOX::Nln::StatusTest::QuantityType& qtype);

      /// \brief Do a recursive search for a <T> status test class.
      template <class T>
      ::NOX::StatusTest::Generic* get_outer_status_test(::NOX::StatusTest::Generic& full_otest);

      /*! \brief Do a recursive search for a <T> status test class and return its status.
       *
       * If more than one of the given status test objects is combined in a combination list,
       * the AND combination of the different status is returned. I.e. if one of the status
       * is unconverged, the return status is unconverged.
       * If we cannot find the given status test class, a default value of -100 is returned.
       *
       * \param test (in) : StatusTest object which will be scanned.
       */
      template <class T>
      int get_outer_status(const ::NOX::StatusTest::Generic& test);

      /*! \brief Convert the quantity type to a solution type
       *
       * \param qtype : Quantity type which has to be converted.
       */
      enum NOX::Nln::SolutionType convert_quantity_type_to_solution_type(
          const enum NOX::Nln::StatusTest::QuantityType& qtype);

      /*! \brief Map norm type stl_string to norm type enum
       *
       * \param name : Name of the vector norm type.
       */
      enum ::NOX::Abstract::Vector::NormType string_to_norm_type(const std::string& name);

      /// add pre/post operator to pre/post operator vector
      void add_to_pre_post_op_vector(
          Teuchos::ParameterList& p_nox_opt, const Teuchos::RCP<::NOX::Observer>& ppo_ptr);

      /// return the name of the parameter list corresponding to the set direction method
      std::string get_direction_method_list_name(const Teuchos::ParameterList& p);

      /// Wrapper to compute vector norms
      template <typename T>
      T calc_vector_norm(const Core::LinAlg::Vector<T>& vector,
          const ::NOX::Abstract::Vector::NormType type, const bool is_scaled = false);
    }  // namespace Aux
  }  // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
