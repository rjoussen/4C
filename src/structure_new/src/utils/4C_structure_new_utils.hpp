// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_UTILS_HPP
#define FOUR_C_STRUCTURE_NEW_UTILS_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"
#include "4C_solver_nonlin_nox_enum_lists.hpp"

#include <NOX_Abstract_Vector.H>

#include <set>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class Solver;
}  // namespace Core::LinAlg

namespace NOX
{
  namespace Nln
  {
    namespace CONSTRAINT
    {
      namespace Interface
      {
        class Required;
        class Preconditioner;
      }  // namespace Interface
    }  // namespace CONSTRAINT

    class Scaling;
  }  // namespace Nln
}  // namespace NOX

namespace Solid
{
  class Integrator;

  namespace TimeInt
  {
    class BaseDataSDyn;
    class BaseDataGlobalState;
  }  // namespace TimeInt

  namespace Nln
  {
    //! Convert the structural vector types to a corresponding nox norm type
    enum ::NOX::Abstract::Vector::NormType convert2_nox_norm_type(
        const enum Inpar::Solid::VectorNorm& normtype);

    /*! Convert the structural model type enums to nox nln solution
     *  type enums
     *
     *  Convert the model type enums to the nox internal solution type
     *  enums. This is necessary, because the nox framework is not
     *  supposed to be restricted to structural problems only.
     */
    void convert_model_type2_sol_type(std::vector<enum NOX::Nln::SolutionType>& soltypes,
        std::map<enum NOX::Nln::SolutionType, Teuchos::RCP<Core::LinAlg::Solver>>& slinsolvers,
        const std::set<enum Inpar::Solid::ModelType>& modeltypes,
        const std::map<enum Inpar::Solid::ModelType, std::shared_ptr<Core::LinAlg::Solver>>&
            mlinsolvers);

    /*! \brief Convert the structural model type enumerator to a nox nln solution
     *  type enumerator
     *
     *  \param modeltype (in) : Model type enumerator which has to be converted
     *  \param do_check  (in) : Check if a corresponding solution type exists */
    enum NOX::Nln::SolutionType convert_model_type2_sol_type(
        const enum Inpar::Solid::ModelType& modeltype, const bool& do_check);

    /*! \brief Convert the structural model type enumerator to a nox nln solution
     *  type enumerator and check if the conversion was successful
     *
     *  \param modeltype (in) : Model type enumerator which has to be converted. */
    inline enum NOX::Nln::SolutionType convert_model_type2_sol_type(
        const enum Inpar::Solid::ModelType& modeltype)
    {
      return convert_model_type2_sol_type(modeltype, true);
    }

    /*! \brief Convert the nox nln solution type enumerator to a structural model
     *  type enumerator
     *
     *  \param soltype (in)   : Solution type enumerator which has to be converted.
     *  \param do_check  (in) : Check if a corresponding model type exists. */
    enum Inpar::Solid::ModelType convert_sol_type2_model_type(
        const enum NOX::Nln::SolutionType& soltype, const bool& do_check);

    /*! \brief Convert the structural model type enumerator to a nox nln solution
     *  type enumerator and check if the conversion was successful
     *
     *  \param soltype (in) : Solution type enumerator which has to be converted. */
    inline enum Inpar::Solid::ModelType convert_sol_type2_model_type(
        const enum NOX::Nln::SolutionType& soltype)
    {
      return convert_sol_type2_model_type(soltype, true);
    }

    /*! \brief Convert the nox nln statustest quantity type enumerator to a structural model
     *  type enumerator
     *
     *  \param qtype (in)    : Quantity type enumerator which has to be converted.
     *  \param do_check (in) : Check if a corresponding model type exists. */
    enum Inpar::Solid::ModelType convert_quantity_type2_model_type(
        const enum NOX::Nln::StatusTest::QuantityType& qtype, const bool& do_check);

    /*! \brief Convert the nox nln statustest quantity type enumerator to a structural model
     *  type enumerator and check if the conversion was successful
     *
     *  \param qtype (in) : Quantity type enumerator which has to be converted. */
    inline enum Inpar::Solid::ModelType convert_quantity_type2_model_type(
        const enum NOX::Nln::StatusTest::QuantityType& qtype)
    {
      return convert_quantity_type2_model_type(qtype, true);
    }

    /*! \brief Convert the nox nln statustest quantity type enumerator to a structural model
     *  type enumerator
     *
     *  \param qtype (in)    : Quantity type enumerator which has to be converted. */
    enum Inpar::Solid::EleTech convert_quantity_type2_ele_tech(
        const enum NOX::Nln::StatusTest::QuantityType& qtype);

    //! Returns the optimization type of the underlying structural problem
    enum NOX::Nln::OptimizationProblemType optimization_type(
        const std::vector<enum NOX::Nln::SolutionType>& soltypes);

    /// convert structure condition number type to a nox condition number type
    NOX::Nln::LinSystem::ConditionNumber convert2_nox_condition_number_type(
        const Inpar::Solid::ConditionNumber stype);

    //! Set the constraint interfaces
    void create_constraint_interfaces(
        std::map<enum NOX::Nln::SolutionType,
            Teuchos::RCP<NOX::Nln::CONSTRAINT::Interface::Required>>& iconstr,
        Solid::Integrator& integrator, const std::vector<enum NOX::Nln::SolutionType>& soltypes);

    //! Set the constraint preconditioner interfaces
    void create_constraint_preconditioner(
        std::map<NOX::Nln::SolutionType,
            Teuchos::RCP<NOX::Nln::CONSTRAINT::Interface::Preconditioner>>& iconstr_prec,
        Solid::Integrator& integrator, const std::vector<enum NOX::Nln::SolutionType>& soltypes);

    //! Create object to scale linear system
    void create_scaling(std::shared_ptr<NOX::Nln::Scaling>& iscale,
        const Solid::TimeInt::BaseDataSDyn& DataSDyn, Solid::TimeInt::BaseDataGlobalState& GState);
  }  // namespace Nln
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
