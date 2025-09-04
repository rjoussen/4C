// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_linalg_utils_densematrix_multiply.hpp"
#include "4C_mat_elchmat.hpp"
#include "4C_mat_elchphase.hpp"
#include "4C_mat_newman_multiscale.hpp"
#include "4C_scatra_ele_calc_elch_diffcond_multiscale.hpp"
#include "4C_scatra_ele_parameter_std.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<distype,
    probdim>::calculate_electrode_soc_and_c_rate(const Core::Elements::Element* const& ele,
    const Core::FE::Discretization& discretization, Core::Elements::LocationArray& la,
    Core::LinAlg::SerialDenseVector& scalars)
{
  // safety check
  if (my::numscal_ != 1)
    FOUR_C_THROW("Electrode state of charge can only be computed for one transported scalar!");

  // extract multi-scale material
  auto elchmat = std::dynamic_pointer_cast<const Mat::ElchMat>(ele->material());
  auto elchphase =
      std::dynamic_pointer_cast<const Mat::ElchPhase>(elchmat->phase_by_id(elchmat->phase_id(0)));
  auto newmanmultiscale =
      std::dynamic_pointer_cast<Mat::NewmanMultiScale>(elchphase->mat_by_id(elchphase->mat_id(0)));

  // initialize variables for integrals of concentration, its time derivative, and domain
  double intconcentration(0.0);
  double intconcentrationtimederiv(0.0);
  double intdomain(0.0);

  // integration points and weights
  const Core::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      ScaTra::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::eval_shape_func_and_derivs_at_int_point(intpoints, iquad);

    // calculate integral of concentration
    intconcentration += newmanmultiscale->evaluate_mean_concentration(iquad) * fac;

    // calculate integral of time derivative of concentration
    intconcentrationtimederiv +=
        newmanmultiscale->evaluate_mean_concentration_time_derivative(iquad) * fac;

    // calculate integral of domain
    intdomain += fac;
  }  // loop over integration points

  // safety check
  if (scalars.length() != 3 and scalars.length() != 6)
    FOUR_C_THROW("Result vector for electrode state of charge computation has invalid length!");

  // write results for concentration and domain integrals into result vector
  scalars(0) = intconcentration;
  scalars(1) = intconcentrationtimederiv;
  scalars(2) = intdomain;

  // set ale quantities to zero
  if (scalars.length() == 6) scalars(3) = scalars(4) = scalars(5) = 0.0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<distype,
    probdim>::calculate_mean_electrode_concentration(const Core::Elements::Element* const& ele,
    const Core::FE::Discretization& discretization, Core::Elements::LocationArray& la,
    Core::LinAlg::SerialDenseVector& conc)
{
  // safety check
  if (my::numscal_ != 1)
    FOUR_C_THROW("Electrode state of charge can only be computed for one transported scalar!");

  // extract multi-scale material
  auto elchmat = std::dynamic_pointer_cast<const Mat::ElchMat>(ele->material());
  auto elchphase =
      std::dynamic_pointer_cast<const Mat::ElchPhase>(elchmat->phase_by_id(elchmat->phase_id(0)));
  auto newmanmultiscale =
      std::dynamic_pointer_cast<Mat::NewmanMultiScale>(elchphase->mat_by_id(elchphase->mat_id(0)));

  // integration points and weights
  const Core::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      ScaTra::DisTypeToOptGaussRule<distype>::rule);

  if (intpoints.ip().nquad != nen_)
    FOUR_C_THROW(
        "number of element nodes must equal number of Gauss points for reasonable projection");

  // matrix of shape functions evaluated at Gauss points
  Core::LinAlg::SerialDenseMatrix N(nen_, nen_);

  // Gauss point concentration of electrode
  Core::LinAlg::SerialDenseVector conc_gp(nen_);

  // loop over integration points
  for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    my::eval_shape_func_and_derivs_at_int_point(intpoints, iquad);

    // calculate mean concentration at Gauss point and store in vector
    const double concentration_gp = newmanmultiscale->evaluate_mean_concentration(iquad);
    conc_gp(iquad) = concentration_gp;

    // build matrix of shape functions
    for (int node = 0; node < static_cast<int>(nen_); ++node) N(iquad, node) = my::funct_(node, 0);
  }

  // conc_gp = N * conc --> conc = N^-1 * conc_gp
  using ordinalType = Core::LinAlg::SerialDenseMatrix::ordinalType;
  using scalarType = Core::LinAlg::SerialDenseMatrix::scalarType;
  Teuchos::SerialDenseSolver<ordinalType, scalarType> invert;
  invert.setMatrix(Teuchos::rcpFromRef(N.base()));
  invert.invert();
  Core::LinAlg::multiply(conc, N, conc_gp);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<distype, probdim>::calculate_scalars(
    const Core::Elements::Element* ele, Core::LinAlg::SerialDenseVector& scalars,
    const bool inverting, const bool calc_grad_phi)
{
  my::calculate_scalars(ele, scalars, inverting, calc_grad_phi);

  // extract multi-scale material
  auto elchmat = std::dynamic_pointer_cast<const Mat::ElchMat>(ele->material());
  auto elchphase =
      std::dynamic_pointer_cast<const Mat::ElchPhase>(elchmat->phase_by_id(elchmat->phase_id(0)));
  auto newmanmultiscale =
      std::dynamic_pointer_cast<Mat::NewmanMultiScale>(elchphase->mat_by_id(elchphase->mat_id(0)));

  // initialize variables for integrals of concentration, its time derivative, and domain
  double intconcentration(0.0);

  // integration points and weights
  const Core::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      ScaTra::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::eval_shape_func_and_derivs_at_int_point(intpoints, iquad);

    // calculate integral of concentration
    intconcentration += newmanmultiscale->evaluate_mean_concentration(iquad) * fac;
  }

  scalars(scalars.length() - 1) = intconcentration;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
int Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<distype, probdim>::evaluate_action(
    Core::Elements::Element* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, const ScaTra::Action& action,
    Core::Elements::LocationArray& la, Core::LinAlg::SerialDenseMatrix& elemat1,
    Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseVector& elevec2, Core::LinAlg::SerialDenseVector& elevec3)
{
  // extract multi-scale material
  auto elchmat = std::dynamic_pointer_cast<const Mat::ElchMat>(ele->material());
  auto elchphase =
      std::dynamic_pointer_cast<const Mat::ElchPhase>(elchmat->phase_by_id(elchmat->phase_id(0)));
  auto newmanmultiscale =
      std::dynamic_pointer_cast<Mat::NewmanMultiScale>(elchphase->mat_by_id(elchphase->mat_id(0)));

  // determine and evaluate action
  switch (action)
  {
    case ScaTra::Action::micro_scale_initialize:
    {
      const Core::FE::IntPointsAndWeights<nsd_ele_> intpoints(
          ScaTra::DisTypeToOptGaussRule<distype>::rule);

      // loop over all Gauss points
      for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
        // initialize micro scale in multi-scale simulations
        newmanmultiscale->initialize(ele->id(), iquad, my::scatrapara_->is_ale());

      break;
    }

    case ScaTra::Action::micro_scale_prepare_time_step:
    case ScaTra::Action::micro_scale_solve:
    {
      // extract state variables at element nodes
      Core::FE::extract_my_values<Core::LinAlg::Matrix<nen_, 1>>(
          *discretization.get_state("phinp"), my::ephinp_, la[0].lm_);

      const Core::FE::IntPointsAndWeights<nsd_ele_> intpoints(
          ScaTra::DisTypeToOptGaussRule<distype>::rule);

      // loop over all Gauss points
      for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
      {
        // evaluate shape functions at Gauss point
        this->eval_shape_func_and_derivs_at_int_point(intpoints, iquad);

        // evaluate state variables at Gauss point
        this->set_internal_variables_for_mat_and_rhs();

        // initialize vector with macro-scale state variables
        std::vector<double> phinp(3, 0.);
        phinp[0] = my::scatravarmanager_->phinp(0);
        phinp[1] = my::funct_.dot(my::ephinp_[1]);
        phinp[2] = my::funct_.dot(my::ephinp_[2]);

        if (action == ScaTra::Action::micro_scale_prepare_time_step)
        {
          // prepare time step on micro scale
          newmanmultiscale->prepare_time_step(iquad, phinp);
        }
        else
        {
          // solve micro scale
          std::vector<double> dummy(3, 0.0);
          const double detF = my::eval_det_f_at_int_point(ele, intpoints, iquad);
          newmanmultiscale->evaluate(iquad, phinp, dummy[0], dummy, detF);
        }
      }

      break;
    }

    case ScaTra::Action::micro_scale_update:
    {
      const Core::FE::IntPointsAndWeights<nsd_ele_> intpoints(
          ScaTra::DisTypeToOptGaussRule<distype>::rule);

      // loop over all Gauss points
      for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
        // update multi-scale scalar transport material
        newmanmultiscale->update(iquad);

      break;
    }

    case ScaTra::Action::micro_scale_output:
    {
      const Core::FE::IntPointsAndWeights<nsd_ele_> intpoints(
          ScaTra::DisTypeToOptGaussRule<distype>::rule);

      // loop over all Gauss points
      for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
        // create output on micro scale
        newmanmultiscale->output(iquad);

      break;
    }

    case ScaTra::Action::collect_micro_scale_output:
    {
      const Core::FE::IntPointsAndWeights<nsd_ele_> intpoints(
          ScaTra::DisTypeToOptGaussRule<distype>::rule);

      // loop over all Gauss points
      for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
      {
        // create output on micro scale
        newmanmultiscale->collect_output_data(iquad);
      }

      break;
    }

    case ScaTra::Action::micro_scale_read_restart:
    {
      const Core::FE::IntPointsAndWeights<nsd_ele_> intpoints(
          ScaTra::DisTypeToOptGaussRule<distype>::rule);

      // loop over all Gauss points
      for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
        // read restart on micro scale
        newmanmultiscale->read_restart(iquad);

      break;
    }
    case ScaTra::Action::micro_scale_set_time:
    {
      const Core::FE::IntPointsAndWeights<nsd_ele_> intpoints(
          ScaTra::DisTypeToOptGaussRule<distype>::rule);
      for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
        newmanmultiscale->set_time_stepping(
            iquad, params.get<double>("dt"), params.get<double>("time"), params.get<int>("step"));
      [[fallthrough]];
    }

    default:
    {
      mydiffcond::evaluate_action(
          ele, params, discretization, action, la, elemat1, elemat2, elevec1, elevec2, elevec3);

      break;
    }
  }  // switch(action)

  return -1;
}

// template classes
// 1D elements
template class Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::line2, 1>;
template class Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::line2, 2>;
template class Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::line2, 3>;
template class Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::line3, 1>;

// 2D elements
template class Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::tri3, 2>;
template class Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::tri3, 3>;
template class Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::tri6, 2>;
template class Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::quad4, 2>;
template class Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::quad4, 3>;
template class Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::quad9, 2>;
template class Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::nurbs9,
    2>;

// 3D elements
template class Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::hex8, 3>;
template class Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::hex27, 3>;
template class Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::tet4, 3>;
template class Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::tet10, 3>;
template class Discret::Elements::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::pyramid5,
    3>;

FOUR_C_NAMESPACE_CLOSE
