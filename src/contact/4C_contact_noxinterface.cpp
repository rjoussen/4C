// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_contact_noxinterface.hpp"

#include "4C_contact_abstract_strategy.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_solver_nonlin_nox_aux.hpp"
#include "4C_utils_enum.hpp"

#include <NOX_Epetra_Vector.H>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::NoxInterface::NoxInterface()
    : isinit_(false),
      issetup_(false),
      strategy_ptr_(nullptr),
      cycling_maps_(std::vector<std::shared_ptr<Core::LinAlg::Map>>(0))
{
  // should stay empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::NoxInterface::init(const std::shared_ptr<CONTACT::AbstractStrategy>& strategy_ptr)
{
  issetup_ = false;

  strategy_ptr_ = strategy_ptr;

  // set flag at the end
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::NoxInterface::setup()
{
  check_init();

  // set flag at the end
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::NoxInterface::get_constraint_rhs_norms(const Core::LinAlg::Vector<double>& F,
    NOX::Nln::StatusTest::QuantityType checkQuantity, ::NOX::Abstract::Vector::NormType type,
    bool isScaled) const
{
  if (checkQuantity != NOX::Nln::StatusTest::quantity_contact_normal and
      checkQuantity != NOX::Nln::StatusTest::quantity_contact_friction)
    return -1.0;

  std::shared_ptr<const Core::LinAlg::Vector<double>> constrRhs =
      strategy().get_rhs_block_ptr_for_norm_check(CONTACT::VecBlockType::constraint);

  // no contact contributions present
  if (!constrRhs) return 0.0;

  // export the vector to the current redistributed map
  std::shared_ptr<Core::LinAlg::Vector<double>> constrRhs_red = nullptr;
  // Note: PointSameAs is faster than SameAs and should do the job right here,
  // since we replace the map afterwards anyway.               hiermeier 08/17
  if (not constrRhs->get_map().point_same_as(strategy().lm_dof_row_map(true)))
  {
    constrRhs_red = std::make_shared<Core::LinAlg::Vector<double>>(strategy().lm_dof_row_map(true));
    Core::LinAlg::export_to(*constrRhs, *constrRhs_red);
  }
  else
    constrRhs_red = std::make_shared<Core::LinAlg::Vector<double>>(*constrRhs);

  // replace the map
  constrRhs_red->replace_map(strategy().slave_dof_row_map(true));

  double constrNorm = -1.0;
  Teuchos::RCP<const ::NOX::Epetra::Vector> constrRhs_nox = Teuchos::null;
  switch (checkQuantity)
  {
    case NOX::Nln::StatusTest::quantity_contact_normal:
    {
      // create vector with redistributed slave dof row map in normal direction
      std::shared_ptr<Core::LinAlg::Vector<double>> nConstrRhs =
          Core::LinAlg::extract_my_vector(*constrRhs_red, strategy().slave_n_dof_row_map(true));

      constrNorm = NOX::Nln::Aux::calc_vector_norm(*nConstrRhs, type, isScaled);
      break;
    }
    case NOX::Nln::StatusTest::quantity_contact_friction:
    {
      // create vector with redistributed slave dof row map in tangential directions
      std::shared_ptr<Core::LinAlg::Vector<double>> tConstrRhs =
          Core::LinAlg::extract_my_vector(*constrRhs_red, strategy().slave_t_dof_row_map(true));

      constrNorm = NOX::Nln::Aux::calc_vector_norm(*tConstrRhs, type, isScaled);
      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported quantity type!");
      break;
    }
  }

  return constrNorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::NoxInterface::get_lagrange_multiplier_update_rms(
    const Core::LinAlg::Vector<double>& xNew, const Core::LinAlg::Vector<double>& xOld, double aTol,
    double rTol, NOX::Nln::StatusTest::QuantityType checkQuantity,
    bool disable_implicit_weighting) const
{
  if (checkQuantity != NOX::Nln::StatusTest::quantity_contact_normal and
      checkQuantity != NOX::Nln::StatusTest::quantity_contact_friction)
    return -1.0;

  double rms = -1.0;
  std::shared_ptr<Core::LinAlg::Vector<double>> z_ptr = nullptr;
  std::shared_ptr<Core::LinAlg::Vector<double>> zincr_ptr = nullptr;
  switch (checkQuantity)
  {
    case NOX::Nln::StatusTest::quantity_contact_normal:
    {
      // extract vectors with redistributed slave dof row map in normal direction
      z_ptr = Core::LinAlg::extract_my_vector(
          *strategy().lagrange_multiplier_np(true), strategy().slave_n_dof_row_map(true));
      zincr_ptr = Core::LinAlg::extract_my_vector(
          *strategy().lagrange_multiplier_increment(), strategy().slave_n_dof_row_map(true));

      break;
    }
    case NOX::Nln::StatusTest::quantity_contact_friction:
    {
      // extract vectors with redistributed slave dof row map in tangential directions
      z_ptr = Core::LinAlg::extract_my_vector(
          *strategy().lagrange_multiplier_np(true), strategy().slave_t_dof_row_map(true));
      zincr_ptr = Core::LinAlg::extract_my_vector(
          *strategy().lagrange_multiplier_increment(), strategy().slave_t_dof_row_map(true));

      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported quantity type!");
      break;
    }
  }

  rms = NOX::Nln::Aux::root_mean_square_norm(aTol, rTol, *strategy().lagrange_multiplier_np(true),
      *strategy().lagrange_multiplier_increment(), disable_implicit_weighting);

  return rms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::NoxInterface::get_lagrange_multiplier_update_norms(
    const Core::LinAlg::Vector<double>& xNew, const Core::LinAlg::Vector<double>& xOld,
    NOX::Nln::StatusTest::QuantityType checkQuantity, ::NOX::Abstract::Vector::NormType type,
    bool isScaled) const
{
  if (checkQuantity != NOX::Nln::StatusTest::quantity_contact_normal and
      checkQuantity != NOX::Nln::StatusTest::quantity_contact_friction)
    return -1.0;

  if (strategy().lagrange_multiplier_np(true) == nullptr) return 0.;

  std::shared_ptr<Core::LinAlg::Vector<double>> zincr_ptr = nullptr;
  switch (checkQuantity)
  {
    case NOX::Nln::StatusTest::quantity_contact_normal:
    {
      // extract vector with redistributed slave dof row map in normal direction
      zincr_ptr = Core::LinAlg::extract_my_vector(
          *strategy().lagrange_multiplier_increment(), strategy().slave_n_dof_row_map(true));
      break;
    }
    case NOX::Nln::StatusTest::quantity_contact_friction:
    {
      // extract vector with redistributed slave dof row map in tangential directions
      zincr_ptr = Core::LinAlg::extract_my_vector(
          *strategy().lagrange_multiplier_increment(), strategy().slave_t_dof_row_map(true));

      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported quantity type!");
      break;
    }
  }

  return NOX::Nln::Aux::calc_vector_norm(*zincr_ptr, type, isScaled);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::NoxInterface::get_previous_lagrange_multiplier_norms(
    const Core::LinAlg::Vector<double>& xOld, NOX::Nln::StatusTest::QuantityType checkQuantity,
    ::NOX::Abstract::Vector::NormType type, bool isScaled) const
{
  if (checkQuantity != NOX::Nln::StatusTest::quantity_contact_normal and
      checkQuantity != NOX::Nln::StatusTest::quantity_contact_friction)
    return -1.0;

  double zoldnorm = -1.0;

  if (strategy().lagrange_multiplier_np(true) == nullptr) return 0;

  /* lagrange multiplier of the previous Newton step
   * (NOT equal to zOld_, which is stored in the Strategy object!!!) */
  std::shared_ptr<Core::LinAlg::Vector<double>> zold_ptr =
      std::make_shared<Core::LinAlg::Vector<double>>(*strategy().lagrange_multiplier_np(true));
  zold_ptr->update(-1.0, *strategy().lagrange_multiplier_increment(), 1.0);
  std::shared_ptr<::NOX::Epetra::Vector> zold_nox_ptr = nullptr;
  switch (checkQuantity)
  {
    case NOX::Nln::StatusTest::quantity_contact_normal:
    {
      std::shared_ptr<Core::LinAlg::Vector<double>> znold_ptr =
          Core::LinAlg::extract_my_vector(*zold_ptr, strategy().slave_n_dof_row_map(true));

      zoldnorm = NOX::Nln::Aux::calc_vector_norm(*znold_ptr, type, isScaled);
      break;
    }
    case NOX::Nln::StatusTest::quantity_contact_friction:
    {
      std::shared_ptr<Core::LinAlg::Vector<double>> ztold_ptr =
          Core::LinAlg::extract_my_vector(*zold_ptr, strategy().slave_t_dof_row_map(true));

      zoldnorm = NOX::Nln::Aux::calc_vector_norm(*ztold_ptr, type, isScaled);
      break;
    }
    default:
    {
      FOUR_C_THROW("The given quantity type is unsupported!");
      break;
    }
  }

  // avoid very small norm values for the pure inactive case
  if (not strategy().is_in_contact()) zoldnorm = std::max(zoldnorm, 1.0);

  return zoldnorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum ::NOX::StatusTest::StatusType CONTACT::NoxInterface::get_active_set_info(
    NOX::Nln::StatusTest::QuantityType checkQuantity, int& activesetsize) const
{
  const bool semismooth = strategy().params().get<bool>("SEMI_SMOOTH_NEWTON");
  if (not semismooth) FOUR_C_THROW("Currently we support only the semi-smooth Newton case!");
  // ---------------------------------------------------------------------------
  // get the number of active nodes for the given active set type
  // ---------------------------------------------------------------------------
  switch (checkQuantity)
  {
    case NOX::Nln::StatusTest::quantity_contact_normal:
    {
      activesetsize = strategy().number_of_active_nodes();
      break;
    }
    case NOX::Nln::StatusTest::quantity_contact_friction:
    {
      activesetsize = strategy().number_of_slip_nodes();
      break;
    }
    default:
    {
      FOUR_C_THROW("The given quantity type is unsupported!");
      break;
    }
  }
  // ---------------------------------------------------------------------------
  // translate the active set semi-smooth Newton convergence flag
  // ---------------------------------------------------------------------------
  if (strategy().active_set_converged())
    return ::NOX::StatusTest::Converged;
  else
    return ::NOX::StatusTest::Unconverged;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Core::LinAlg::Map> CONTACT::NoxInterface::get_current_active_set_map(
    enum NOX::Nln::StatusTest::QuantityType checkQuantity) const
{
  switch (checkQuantity)
  {
    case NOX::Nln::StatusTest::quantity_contact_normal:
    {
      return Teuchos::rcp(strategy().active_row_nodes());
      break;
    }
    case NOX::Nln::StatusTest::quantity_contact_friction:
    {
      return Teuchos::rcp(strategy().slip_row_nodes());
      break;
    }
    default:
    {
      FOUR_C_THROW("The given active set type is unsupported!");
      break;
    }
  }  // switch (active_set_type)

  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Core::LinAlg::Map> CONTACT::NoxInterface::get_old_active_set_map(
    enum NOX::Nln::StatusTest::QuantityType checkQuantity) const
{
  switch (checkQuantity)
  {
    case NOX::Nln::StatusTest::quantity_contact_normal:
    {
      return Teuchos::rcp(strategy().get_old_active_row_nodes());
      break;
    }
    case NOX::Nln::StatusTest::quantity_contact_friction:
    {
      return Teuchos::rcp(strategy().get_old_slip_row_nodes());
      break;
    }
    default:
    {
      FOUR_C_THROW("The given active set type is unsupported!");
      break;
    }
  }  // switch (active_set_type)

  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::NoxInterface::get_model_value(NOX::Nln::MeritFunction::MeritFctName name) const
{
  switch (name)
  {
    case NOX::Nln::MeritFunction::mrtfct_infeasibility_two_norm:
    case NOX::Nln::MeritFunction::mrtfct_infeasibility_two_norm_active:
    {
      double val = strategy().get_potential_value(name);
      val = std::sqrt(val);

      return val;
    }
    case NOX::Nln::MeritFunction::mrtfct_energy:
    {
      // The energy of the primary field is considered, no contact contribution.
      return 0.0;
    }
    default:
      FOUR_C_THROW("Unsupported Merit function name! (enum = {})", name);
  }

  FOUR_C_THROW("Impossible to reach this point.");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::NoxInterface::get_linearized_model_terms(const Core::LinAlg::Vector<double>& dir,
    const enum NOX::Nln::MeritFunction::MeritFctName name,
    const enum NOX::Nln::MeritFunction::LinOrder linorder,
    const enum NOX::Nln::MeritFunction::LinType lintype) const
{
  switch (name)
  {
    case NOX::Nln::MeritFunction::mrtfct_infeasibility_two_norm:
    case NOX::Nln::MeritFunction::mrtfct_infeasibility_two_norm_active:
    {
      double lin_val =
          strategy().get_linearized_potential_value_terms(dir, name, linorder, lintype);
      const double modelvalue = get_model_value(name);
      if (modelvalue != 0.0) lin_val /= modelvalue;

      return lin_val;
    }
    default:
      FOUR_C_THROW("Unsupported Merit function name! (enum = {})", name);
  }

  FOUR_C_THROW("Impossible to reach this point.");
}

FOUR_C_NAMESPACE_CLOSE
