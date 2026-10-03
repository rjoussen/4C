// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_thermo_ele_impl.hpp"

#include "4C_fem_condition_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_general_utils_nurbs_shapefunctions.hpp"
#include "4C_fem_geometry_position_array.hpp"
#include "4C_fem_nurbs_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_fixedsizematrix_solver.hpp"
#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor_conversion.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_mat_trait_thermo.hpp"
#include "4C_structure_new_input.hpp"
#include "4C_thermo_ele_action.hpp"
#include "4C_thermo_element.hpp"  // only for visualization of element data
#include "4C_thermo_input.hpp"
#include "4C_utils_function.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

Discret::Elements::TemperImplInterface* Discret::Elements::TemperImplInterface::impl(
    const Core::Elements::Element* ele)
{
  switch (ele->shape())
  {
    case Core::FE::CellType::hex8:
    {
      return TemperImpl<Core::FE::CellType::hex8>::instance();
    }
    case Core::FE::CellType::hex20:
    {
      return TemperImpl<Core::FE::CellType::hex20>::instance();
    }
    case Core::FE::CellType::hex27:
    {
      return TemperImpl<Core::FE::CellType::hex27>::instance();
    }
    case Core::FE::CellType::tet4:
    {
      return TemperImpl<Core::FE::CellType::tet4>::instance();
    }
    case Core::FE::CellType::tet10:
    {
      return TemperImpl<Core::FE::CellType::tet10>::instance();
    }
    case Core::FE::CellType::wedge6:
    {
      return TemperImpl<Core::FE::CellType::wedge6>::instance();
    }
    case Core::FE::CellType::pyramid5:
    {
      return TemperImpl<Core::FE::CellType::pyramid5>::instance();
    }
    case Core::FE::CellType::quad4:
    {
      return TemperImpl<Core::FE::CellType::quad4>::instance();
    }
    case Core::FE::CellType::quad8:
    {
      return TemperImpl<Core::FE::CellType::quad8>::instance();
    }
    case Core::FE::CellType::quad9:
    {
      return TemperImpl<Core::FE::CellType::quad9>::instance();
    }
    case Core::FE::CellType::tri3:
    {
      return TemperImpl<Core::FE::CellType::tri3>::instance();
    }
    case Core::FE::CellType::line2:
    {
      return TemperImpl<Core::FE::CellType::line2>::instance();
    }
    case Core::FE::CellType::nurbs27:
    {
      return TemperImpl<Core::FE::CellType::nurbs27>::instance();
    }
    default:
      FOUR_C_THROW("Element shape {} ({} nodes) not activated. Just do it.",
          Core::FE::cell_type_to_string(ele->shape()), ele->num_node());
      break;
  }
  return nullptr;

}  // TemperImperInterface::Impl()

template <Core::FE::CellType distype>
Discret::Elements::TemperImpl<distype>* Discret::Elements::TemperImpl<distype>::instance(
    Core::Utils::SingletonAction action)
{
  static auto singleton_owner = Core::Utils::make_singleton_owner(
      []()
      {
        return std::unique_ptr<Discret::Elements::TemperImpl<distype>>(
            new Discret::Elements::TemperImpl<distype>());
      });

  return singleton_owner.instance(action);
}

template <Core::FE::CellType distype>
Discret::Elements::TemperImpl<distype>::TemperImpl()
    : etempn_(Core::LinAlg::Initialization::uninitialized),
      xyze_(Core::LinAlg::Initialization::zero),
      radiation_(Core::LinAlg::Initialization::uninitialized),
      xsi_(Core::LinAlg::Initialization::zero),
      funct_(Core::LinAlg::Initialization::zero),
      deriv_(Core::LinAlg::Initialization::zero),
      xjm_(Core::LinAlg::Initialization::zero),
      xij_(Core::LinAlg::Initialization::zero),
      derxy_(Core::LinAlg::Initialization::zero),
      fac_(0.0),
      gradtemp_(Core::LinAlg::Initialization::zero),
      heatflux_(Core::LinAlg::Initialization::uninitialized),
      cmat_(Core::LinAlg::Initialization::uninitialized),
      dercmat_(Core::LinAlg::Initialization::zero),
      capacoeff_(0.0),
      dercapa_(0.0)

{
}

template <Core::FE::CellType distype>
int Discret::Elements::TemperImpl<distype>::evaluate(
    const Core::Elements::Element* ele, Teuchos::ParameterList& params,
    const Core::FE::Discretization& discretization, const Core::Elements::LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1,  // Tangent ("stiffness")
    Core::LinAlg::SerialDenseMatrix& elemat2,  // Capacity ("mass")
    Core::LinAlg::SerialDenseVector& elevec1,  // internal force vector
    Core::LinAlg::SerialDenseVector& elevec2,  // external force vector
    Core::LinAlg::SerialDenseVector& elevec3   // capacity vector
)
{
  prepare_nurbs_eval(ele, discretization);

  const auto action = Teuchos::getIntegralValue<Thermo::Action>(params, "action");

  // check length
  if (la[0].size() != nen_ * numdofpernode_) FOUR_C_THROW("Location vector length does not match!");

  // disassemble temperature
  if (discretization.has_state(0, "temperature"))
  {
    std::vector<double> mytempnp((la[0].lm_).size());
    std::shared_ptr<const Core::LinAlg::Vector<double>> tempnp =
        discretization.get_state(0, "temperature");
    if (tempnp == nullptr) FOUR_C_THROW("Cannot get state vector 'tempnp'");
    mytempnp = Core::FE::extract_values(*tempnp, la[0].lm_);
    // build the element temperature
    Core::LinAlg::Matrix<nen_ * numdofpernode_, 1> etempn(mytempnp.data(), true);  // view only!
    etempn_.update(etempn);                                                        // copy
  }

  if (discretization.has_state(0, "last temperature"))
  {
    std::vector<double> mytempn((la[0].lm_).size());
    std::shared_ptr<const Core::LinAlg::Vector<double>> tempn =
        discretization.get_state(0, "last temperature");
    if (tempn == nullptr) FOUR_C_THROW("Cannot get state vector 'tempn'");
    mytempn = Core::FE::extract_values(*tempn, la[0].lm_);
    // build the element temperature
    Core::LinAlg::Matrix<nen_ * numdofpernode_, 1> etemp(mytempn.data(), true);  // view only!
    etemp_.update(etemp);                                                        // copy
  }

  double time = 0.0;

  if (action != Thermo::calc_thermo_energy)
  {
    // extract time
    time = params.get<double>("total time");
  }


  //============================================================================
  // calculate tangent K and internal force F_int = K * Theta
  // --> for static case
  if (action == Thermo::calc_thermo_fintcond)
  {
    // set views
    Core::LinAlg::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_> etang(
        elemat1.values(), true);                                                   // view only!
    Core::LinAlg::Matrix<nen_ * numdofpernode_, 1> efint(elevec1.values(), true);  // view only!
    // ecapa, efext, efcap not needed for this action
    // econd: conductivity matrix
    // etang: tangent of thermal problem.
    // --> If dynamic analysis, i.e. T' != 0 --> etang consists of econd AND ecapa

    evaluate_tang_capa_fint(
        ele, time, discretization, la, &etang, nullptr, nullptr, &efint, params);
  }
  //============================================================================
  // calculate only the internal force F_int, needed for restart
  else if (action == Thermo::calc_thermo_fint)
  {
    // set views
    Core::LinAlg::Matrix<nen_ * numdofpernode_, 1> efint(elevec1.values(), true);  // view only!
    // etang, ecapa, efext, efcap not needed for this action

    evaluate_tang_capa_fint(
        ele, time, discretization, la, nullptr, nullptr, nullptr, &efint, params);
  }

  //============================================================================
  // calculate the capacity matrix and the internal force F_int
  // --> for dynamic case, called only once in determine_capa_consist_temp_rate()
  else if (action == Thermo::calc_thermo_fintcapa)
  {
    // set views
    Core::LinAlg::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_> ecapa(
        elemat2.values(), true);                                                   // view only!
    Core::LinAlg::Matrix<nen_ * numdofpernode_, 1> efint(elevec1.values(), true);  // view only!
    // etang, efext, efcap not needed for this action

    evaluate_tang_capa_fint(
        ele, time, discretization, la, nullptr, &ecapa, nullptr, &efint, params);

    // lumping
    if (params.get<bool>("lump capa matrix", false))
    {
      const auto timint =
          params.get<Thermo::DynamicType>("time integrator", Thermo::DynamicType::Undefined);
      switch (timint)
      {
        case Thermo::DynamicType::OneStepTheta:
        {
          calculate_lump_matrix(&ecapa);

          break;
        }
        case Thermo::DynamicType::GenAlpha:
        case Thermo::DynamicType::Statics:
        {
          FOUR_C_THROW("Lumped capacity matrix has not yet been tested");
          break;
        }
        case Thermo::DynamicType::Undefined:
        default:
        {
          FOUR_C_THROW("Undefined time integration scheme for thermal problem!");
          break;
        }
      }
    }
  }

  //============================================================================
  // called from overloaded function apply_force_tang_internal(), exclusively for
  // dynamic-timint (as OST, GenAlpha)
  // calculate effective dynamic tangent matrix K_{T, effdyn},
  // i.e. sum consistent capacity matrix C + its linearization and scaled conductivity matrix
  // --> for dynamic case
  else if (action == Thermo::calc_thermo_finttang)
  {
    // set views
    Core::LinAlg::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_> etang(
        elemat1.values(), true);  // view only!
    Core::LinAlg::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_> ecapa(
        Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<nen_ * numdofpernode_, 1> efint(elevec1.values(), true);  // view only!
    Core::LinAlg::Matrix<nen_ * numdofpernode_, 1> efcap(elevec3.values(), true);  // view only!

    // etang: effective dynamic tangent of thermal problem
    // --> etang == k_{T,effdyn}^{(e)} = timefac_capa ecapa + timefac_cond econd
    // econd: conductivity matrix
    // ecapa: capacity matrix
    // --> If dynamic analysis, i.e. T' != 0 --> etang consists of econd AND ecapa

    // helper matrix to store partial dC/dT*(T_{n+1} - T_n) linearization of capacity
    Core::LinAlg::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_> ecapalin(
        Core::LinAlg::Initialization::zero);

    evaluate_tang_capa_fint(
        ele, time, discretization, la, &etang, &ecapa, &ecapalin, &efint, params);


    if (params.get<bool>("lump capa matrix", false))
    {
      calculate_lump_matrix(&ecapa);
    }

    // explicitly insert capacity matrix into corresponding matrix if existing
    if (elemat2.values() != nullptr)
    {
      Core::LinAlg::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_> ecapa_export(
          elemat2.values(), true);  // view only!
      ecapa_export.update(ecapa);
    }

    // BUILD EFFECTIVE TANGENT AND RESIDUAL ACCORDING TO TIME INTEGRATOR
    // combine capacity and conductivity matrix to one global tangent matrix
    // check the time integrator
    // K_T = fac_capa . C + fac_cond . K
    const auto timint =
        params.get<Thermo::DynamicType>("time integrator", Thermo::DynamicType::Undefined);
    switch (timint)
    {
      case Thermo::DynamicType::Statics:
      {
        // continue
        break;
      }
      case Thermo::DynamicType::OneStepTheta:
      {
        // extract time values from parameter list
        const double theta = params.get<double>("theta");
        const double stepsize = params.get<double>("delta time");

        // ---------------------------------------------------------- etang
        // combine capacity and conductivity matrix to one global tangent matrix
        // etang = 1/Dt . ecapa + theta . econd
        // fac_capa = 1/Dt
        // fac_cond = theta
        etang.update(1.0 / stepsize, ecapa, theta);
        // add additional linearization term from variable capacity
        // + 1/Dt. ecapalin
        etang.update(1.0 / stepsize, ecapalin, 1.0);

        // ---------------------------------------------------------- efcap
        // fcapn = ecapa(T_{n+1}) .  (T_{n+1} -T_n) /Dt
        efcap.multiply(ecapa, etempn_);
        efcap.multiply(-1.0, ecapa, etemp_, 1.0);
        efcap.scale(1.0 / stepsize);
        break;
      }  // ost

      case Thermo::DynamicType::GenAlpha:
      {
        // extract time values from parameter list
        const double alphaf = params.get<double>("alphaf");
        const double alpham = params.get<double>("alpham");
        const double gamma = params.get<double>("gamma");
        const double stepsize = params.get<double>("delta time");

        // ---------------------------------------------------------- etang
        // combined tangent and conductivity matrix to one global matrix
        // etang = alpham/(gamma . Dt) . ecapa + alphaf . econd
        // fac_capa = alpham/(gamma . Dt)
        // fac_cond = alphaf
        double fac_capa = alpham / (gamma * stepsize);
        etang.update(fac_capa, ecapa, alphaf);

        // ---------------------------------------------------------- efcap
        // efcap = ecapa . R_{n+alpham}
        if (discretization.has_state(0, "mid-temprate"))
        {
          std::shared_ptr<const Core::LinAlg::Vector<double>> ratem =
              discretization.get_state(0, "mid-temprate");
          if (ratem == nullptr) FOUR_C_THROW("Cannot get mid-temprate state vector for fcap");
          std::vector<double> myratem((la[0].lm_).size());
          // fill the vector myratem with the global values of ratem
          myratem = Core::FE::extract_values(*ratem, la[0].lm_);
          // build the element mid-temperature rates
          Core::LinAlg::Matrix<nen_ * numdofpernode_, 1> eratem(
              myratem.data(), true);  // view only!
          efcap.multiply(ecapa, eratem);
        }  // ratem != nullptr
        break;

      }  // genalpha
      case Thermo::DynamicType::Undefined:
      default:
      {
        FOUR_C_THROW("Don't know what to do...");
        break;
      }
    }  // end of switch(timint)
  }  // action == Thermo::calc_thermo_finttang

  //============================================================================
  // Calculate/ evaluate heatflux q and temperature gradients gradtemp at
  // gauss points
  else if (action == Thermo::calc_thermo_heatflux)
  {
    // set views
    // efext, efcap not needed for this action, elemat1+2,elevec1-3 are not used anyway

    // working arrays
    Core::LinAlg::Matrix<nquad_, nsd_> eheatflux(Core::LinAlg::Initialization::uninitialized);
    Core::LinAlg::Matrix<nquad_, nsd_> etempgrad(Core::LinAlg::Initialization::uninitialized);

    // if ele is a thermo element --> the Thermo element method KinType() exists
    const auto* therm = dynamic_cast<const Thermo::Element*>(ele);
    const Solid::KinemType kintype = therm->kin_type();
    // thermal problem or geometrically linear TSI problem
    if (kintype == Solid::KinemType::linear)
    {
      linear_heatflux_tempgrad(ele, &eheatflux, &etempgrad);
    }  // TSI: (kintype_ == Solid::KinemType::linear)

    // geometrically nonlinear TSI problem
    if (kintype == Solid::KinemType::nonlinearTotLag)
    {
      // if it's a TSI problem and there are current displacements/velocities
      if (la.size() > 1)
      {
        if ((discretization.has_state(1, "displacement")) and
            (discretization.has_state(1, "velocity")))
        {
          std::vector<double> mydisp(((la[0].lm_).size()) * nsd_, 0.0);
          std::vector<double> myvel(((la[0].lm_).size()) * nsd_, 0.0);

          extract_disp_vel(discretization, la, mydisp, myvel);

          nonlinear_heatflux_tempgrad(ele, mydisp, myvel, &eheatflux, &etempgrad, params);
        }
      }
    }

    // Fill element-center averaged Multivectors for runtime VTK output
    auto heatflux = params.get<std::shared_ptr<Core::LinAlg::MultiVector<double>>>("heatflux");
    auto tempgrad = params.get<std::shared_ptr<Core::LinAlg::MultiVector<double>>>("tempgrad");
    int lid = ele->lid();
    if (lid != -1)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        {
          double s = 0.0;
          for (int jquad = 0; jquad < nquad_; ++jquad) s += eheatflux(jquad, idim);
          s /= nquad_;
          heatflux->get_vector(idim).get_values()[lid] = s;
        }
        {
          double s = 0.0;
          for (int jquad = 0; jquad < nquad_; ++jquad) s += etempgrad(jquad, idim);
          s /= nquad_;
          tempgrad->get_vector(idim).get_values()[lid] = s;
        }
      }
    }
  }  // action == Thermo::calc_thermo_heatflux

  //============================================================================
  else if (action == Thermo::integrate_shape_functions)
  {
    // calculate integral of shape functions
    const auto dofids = params.get<std::shared_ptr<Core::LinAlg::IntSerialDenseVector>>("dofids");
    integrate_shape_functions(ele, elevec1, *dofids);
  }

  //============================================================================
  else if (action == Thermo::calc_thermo_update_istep)
  {
    // call material specific update
    std::shared_ptr<Core::Mat::Material> material = ele->material();
    // we have to have a thermo-capable material here -> throw error if not
    std::shared_ptr<Mat::Trait::Thermo> thermoMat =
        std::dynamic_pointer_cast<Mat::Trait::Thermo>(material);

    Core::FE::IntPointsAndWeights<nsd_> intpoints(Thermo::DisTypeToOptGaussRule<distype>::rule);
    if (intpoints.ip().nquad != nquad_) FOUR_C_THROW("Trouble with number of Gauss points");
  }

  //==================================================================================
  // allowing the predictor TangTemp in input file --> can be decisive in compressible case!
  else if (action == Thermo::calc_thermo_reset_istep)
  {
    // we have to have a thermo-capable material here -> throw error if not
    std::shared_ptr<Mat::Trait::Thermo> thermoMat =
        std::dynamic_pointer_cast<Mat::Trait::Thermo>(ele->material());
    thermoMat->reset_current_state();
  }

  //============================================================================
  // evaluation of internal thermal energy
  else if (action == Thermo::calc_thermo_energy)
  {
    // check length of elevec1
    if (elevec1.length() < 1) FOUR_C_THROW("The given result vector is too short.");

    // get node coordinates
    Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, nen_>>(
        ele, xyze_);

    // declaration of internal variables
    double intenergy = 0.0;

    // ----------------------------- integration loop for one element

    // integrations points and weights
    Core::FE::IntPointsAndWeights<nsd_> intpoints(Thermo::DisTypeToOptGaussRule<distype>::rule);
    if (intpoints.ip().nquad != nquad_) FOUR_C_THROW("Trouble with number of Gauss points");

    // --------------------------------------- loop over Gauss Points
    for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
    {
      eval_shape_func_and_derivs_at_int_point(intpoints, iquad, ele->id());

      // call material law => sets capacoeff_
      materialize(ele, iquad);

      Core::LinAlg::Matrix<1, 1> temp(Core::LinAlg::Initialization::uninitialized);
      temp.multiply_tn(funct_, etempn_);

      // internal energy
      intenergy += capacoeff_ * fac_ * temp(0, 0);

    }  // -------------------------------- end loop over Gauss Points

    elevec1(0) = intenergy;

  }  // evaluation of internal energy

  //============================================================================
  // add linearistion of velocity for dynamic time integration to the stiffness term
  // calculate thermal mechanical tangent matrix K_Td
  else if (action == Thermo::calc_thermo_coupltang)
  {
    Core::LinAlg::Matrix<nen_ * numdofpernode_, nen_ * nsd_ * numdofpernode_> etangcoupl(
        elemat1.values(), true);

    // if it's a TSI problem and there are the current displacements/velocities
    evaluate_coupled_tang(ele, discretization, la, &etangcoupl, params);

  }  // action == "calc_thermo_coupltang"
  //============================================================================
  else if (action == Thermo::calc_thermo_error)
  {
    Core::LinAlg::Matrix<nen_ * numdofpernode_, 1> evector(elevec1.values(), true);  // view only!

    compute_error(ele, evector, params);
  }
  //============================================================================
  else
  {
    FOUR_C_THROW("Unknown type of action for Temperature Implementation: {}", action);
  }


  return 0;
}

template <Core::FE::CellType distype>
int Discret::Elements::TemperImpl<distype>::evaluate_neumann(const Core::Elements::Element* ele,
    const Teuchos::ParameterList& params, const Core::FE::Discretization& discretization,
    const std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  // prepare nurbs
  prepare_nurbs_eval(ele, discretization);

  // check length
  if (lm.size() != nen_ * numdofpernode_) FOUR_C_THROW("Location vector length does not match!");
  // set views
  Core::LinAlg::Matrix<nen_ * numdofpernode_, 1> efext(elevec1, true);  // view only!
  // disassemble temperature
  if (discretization.has_state(0, "temperature"))
  {
    std::vector<double> mytempnp(lm.size());
    std::shared_ptr<const Core::LinAlg::Vector<double>> tempnp =
        discretization.get_state("temperature");
    if (tempnp == nullptr) FOUR_C_THROW("Cannot get state vector 'tempnp'");
    mytempnp = Core::FE::extract_values(*tempnp, lm);
    Core::LinAlg::Matrix<nen_ * numdofpernode_, 1> etemp(mytempnp.data(), true);  // view only!
    etempn_.update(etemp);                                                        // copy
  }
  // check for the action parameter
  const auto action = Teuchos::getIntegralValue<Thermo::Action>(params, "action");
  // extract time
  const double time = params.get<double>("total time");

  // perform actions
  if (action == Thermo::calc_thermo_fext)
  {
    // so far we assume deformation INdependent external loads, i.e. NO
    // difference between geometrically (non)linear TSI

    // we prescribe a scalar value on the volume, constant for (non)linear analysis
    evaluate_fext(ele, time, efext);
  }
  else
  {
    FOUR_C_THROW("Unknown type of action for Temperature Implementation: {}", action);
  }

  return 0;
}

template <Core::FE::CellType distype>
void Discret::Elements::TemperImpl<distype>::evaluate_tang_capa_fint(
    const Core::Elements::Element* ele, const double time,
    const Core::FE::Discretization& discretization, const Core::Elements::LocationArray& la,
    Core::LinAlg::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* etang,
    Core::LinAlg::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* ecapa,
    Core::LinAlg::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* ecapalin,
    Core::LinAlg::Matrix<nen_ * numdofpernode_, 1>* efint, Teuchos::ParameterList& params)
{
  const auto* therm = dynamic_cast<const Thermo::Element*>(ele);
  const Solid::KinemType kintype = therm->kin_type();

  // initialise the vectors
  // evaluate() is called the first time in Thermo::BaseAlgorithm: at this stage
  // the coupling field is not yet known. Pass coupling vectors filled with zeros
  // the size of the vectors is the length of the location vector*nsd_
  std::vector<double> mydisp(((la[0].lm_).size()) * nsd_, 0.0);
  std::vector<double> myvel(((la[0].lm_).size()) * nsd_, 0.0);

  // if it's a TSI problem with displacementcoupling_ --> go on here!
  if (la.size() > 1)
  {
    extract_disp_vel(discretization, la, mydisp, myvel);
  }  // la.Size>1

  // geometrically linear TSI problem
  if ((kintype == Solid::KinemType::linear))
  {
    // purely thermal contributions (the mechanical heat source is evaluated by the structural
    // elements)
    linear_thermo_contribution(ele, time, etang,
        ecapa,     // capa matric
        ecapalin,  // capa linearization
        efint);
  }  // TSI: (kintype_ == Solid::KinemType::linear)

  // geometrically nonlinear TSI problem
  else if (kintype == Solid::KinemType::nonlinearTotLag)
  {
    nonlinear_thermo_disp_contribution(
        ele, time, mydisp, myvel, etang, ecapa, ecapalin, efint, params);
  }  // TSI: (kintype_ == Solid::KinemType::nonlinearTotLag)
}

template <Core::FE::CellType distype>
void Discret::Elements::TemperImpl<distype>::evaluate_coupled_tang(
    const Core::Elements::Element* ele, const Core::FE::Discretization& discretization,
    const Core::Elements::LocationArray& la,
    Core::LinAlg::Matrix<nen_ * numdofpernode_, nen_ * nsd_ * numdofpernode_>* etangcoupl,
    Teuchos::ParameterList& params)
{
  const auto* therm = dynamic_cast<const Thermo::Element*>(ele);
  const Solid::KinemType kintype = therm->kin_type();

  if (la.size() > 1)
  {
    std::vector<double> mydisp(((la[0].lm_).size()) * nsd_, 0.0);
    std::vector<double> myvel(((la[0].lm_).size()) * nsd_, 0.0);

    extract_disp_vel(discretization, la, mydisp, myvel);

    // if there is a strucutural vector available go on here
    // --> calculate coupling stiffness term in case of monolithic TSI
    // For geometrically linear problems, the thermal terms do not depend on the displacements.
    // The derivative of the mechanical heat source is evaluated by the structural elements.

    // geometrically nonlinear TSI problem
    if (kintype == Solid::KinemType::nonlinearTotLag)
    {
      nonlinear_coupled_tang(ele, mydisp, myvel, etangcoupl, params);
    }
  }
}

template <Core::FE::CellType distype>
void Discret::Elements::TemperImpl<distype>::evaluate_fext(
    const Core::Elements::Element* ele,                    // the element whose matrix is calculated
    const double time,                                     // current time
    Core::LinAlg::Matrix<nen_ * numdofpernode_, 1>& efext  // external force
)
{
  // get node coordinates
  Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, nen_>>(
      ele, xyze_);

  // ------------------------------- integration loop for one element

  // integrations points and weights
  Core::FE::IntPointsAndWeights<nsd_> intpoints(Thermo::DisTypeToOptGaussRule<distype>::rule);
  if (intpoints.ip().nquad != nquad_) FOUR_C_THROW("Trouble with number of Gauss points");

  // ----------------------------------------- loop over Gauss Points
  for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
  {
    eval_shape_func_and_derivs_at_int_point(intpoints, iquad, ele->id());

    // ---------------------------------------------------------------------
    // call routine for calculation of radiation in element nodes
    // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
    // ---------------------------------------------------------------------
    radiation(ele, time);
    // fext = fext + N . r. detJ . w(gp)
    // with funct_: shape functions, fac_:detJ . w(gp)
    efext.multiply_nn(fac_, funct_, radiation_, 1.0);
  }
}


/*----------------------------------------------------------------------*
 | calculate system matrix and rhs r_T(T), k_TT(T) (public) g.bau 08/08 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::Elements::TemperImpl<distype>::linear_thermo_contribution(
    const Core::Elements::Element* ele,  // the element whose matrix is calculated
    const double time,                   // current time
    Core::LinAlg::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
        econd,  // conductivity matrix
    Core::LinAlg::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* ecapa,  // capacity matrix
    Core::LinAlg::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
        ecapalin,                                          // linearization contribution of capacity
    Core::LinAlg::Matrix<nen_ * numdofpernode_, 1>* efint  // internal force
)
{
  // get node coordinates
  Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, nen_>>(
      ele, xyze_);

  // ------------------------------- integration loop for one element

  // integrations points and weights
  Core::FE::IntPointsAndWeights<nsd_> intpoints(Thermo::DisTypeToOptGaussRule<distype>::rule);
  if (intpoints.ip().nquad != nquad_) FOUR_C_THROW("Trouble with number of Gauss points");

  // ----------------------------------------- loop over Gauss Points
  for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
  {
    eval_shape_func_and_derivs_at_int_point(intpoints, iquad, ele->id());

    // gradient of current temperature value
    // grad T = d T_j / d x_i = L . N . T = B_ij T_j
    gradtemp_.multiply_nn(derxy_, etempn_);

    // call material law => cmat_,heatflux_
    // negative q is used for balance equation: -q = -(-k gradtemp)= k * gradtemp
    materialize(ele, iquad);


    // internal force vector
    if (efint != nullptr)
    {
      // fint = fint + B^T . q . detJ . w(gp)
      efint->multiply_tn(fac_, derxy_, heatflux_, 1.0);
    }

    // conductivity matrix
    if (econd != nullptr)
    {
      // ke = ke + ( B^T . C_mat . B ) * detJ * w(gp)  with C_mat = k * I
      Core::LinAlg::Matrix<nsd_, nen_> aop(Core::LinAlg::Initialization::uninitialized);  // (3x8)
      // -q = C * B
      aop.multiply_nn(cmat_, derxy_);              //(nsd_xnsd_)(nsd_xnen_)
      econd->multiply_tn(fac_, derxy_, aop, 1.0);  //(nen_xnen_)=(nen_xnsd_)(nsd_xnen_)

      // linearization of non-constant conductivity
      Core::LinAlg::Matrix<nen_, 1> dNgradT(Core::LinAlg::Initialization::uninitialized);
      dNgradT.multiply_tn(derxy_, gradtemp_);
      // TODO only valid for isotropic case
      econd->multiply_nt(dercmat_(0, 0) * fac_, dNgradT, funct_, 1.0);
    }

    // capacity matrix (equates the mass matrix in the structural field)
    if (ecapa != nullptr)
    {
      // ce = ce + ( N^T .  (rho * C_V) . N ) * detJ * w(gp)
      // (8x8)      (8x1)               (1x8)
      // caution: funct_ implemented as (8,1)--> use transposed in code for
      // theoretic part
      ecapa->multiply_nt((fac_ * capacoeff_), funct_, funct_, 1.0);
    }

    if (ecapalin != nullptr)
    {
      // calculate additional linearization d(C(T))/dT (3-tensor!)
      // multiply with temperatures to obtain 2-tensor
      //
      // ecapalin = dC/dT*(T_{n+1} -T_{n})
      //          = fac . dercapa . (T_{n+1} -T_{n}) . (N . N^T . T)^T
      Core::LinAlg::Matrix<1, 1> Netemp(Core::LinAlg::Initialization::uninitialized);
      Core::LinAlg::Matrix<numdofpernode_ * nen_, 1> difftemp(
          Core::LinAlg::Initialization::uninitialized);
      Core::LinAlg::Matrix<numdofpernode_ * nen_, 1> NNetemp(
          Core::LinAlg::Initialization::uninitialized);
      // T_{n+1} - T_{n}
      difftemp.update(1.0, etempn_, -1.0, etemp_);
      Netemp.multiply_tn(funct_, difftemp);
      NNetemp.multiply_nn(funct_, Netemp);
      ecapalin->multiply_nt((fac_ * dercapa_), NNetemp, funct_, 1.0);
    }

  }  // --------------------------------- end loop over Gauss Points

}  // linear_thermo_contribution


/*----------------------------------------------------------------------*
 | calculate coupled fraction for the system matrix          dano 05/10 |
 | and rhs: r_T(d), k_TT(d) (public)                                    |
 *----------------------------------------------------------------------*/

template <Core::FE::CellType distype>
void Discret::Elements::TemperImpl<distype>::nonlinear_thermo_disp_contribution(
    const Core::Elements::Element* ele,  // the element whose matrix is calculated
    const double time,                   // current time
    const std::vector<double>& disp,     // current displacements
    const std::vector<double>& vel,      // current velocities
    Core::LinAlg::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
        econd,  // conductivity matrix
    Core::LinAlg::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* ecapa,  // capacity matrix
    Core::LinAlg::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
        ecapalin,  //!< partial linearization dC/dT of capacity matrix
    Core::LinAlg::Matrix<nen_ * numdofpernode_, 1>* efint,  // internal force
    Teuchos::ParameterList& params)
{
  // The mechanical heat source is evaluated by the structural elements. Here, only the
  // conduction on the deformed configuration and the capacity are evaluated.

  // update element geometry
  Core::LinAlg::Matrix<nen_, nsd_> xcurr;      // current  coord. of element
  Core::LinAlg::Matrix<nen_, nsd_> xcurrrate;  // current  coord. of element
  initial_and_current_nodal_position_velocity(ele, disp, vel, xcurr, xcurrrate);

  // build the deformation gradient w.r.t. material configuration
  Core::LinAlg::Matrix<nsd_, nsd_> defgrd(Core::LinAlg::Initialization::uninitialized);
  // inverse of deformation gradient
  Core::LinAlg::Matrix<nsd_, nsd_> invdefgrd(Core::LinAlg::Initialization::uninitialized);

  // ----------------------------------- integration loop for one element

  // integrations points and weights
  Core::FE::IntPointsAndWeights<nsd_> intpoints(Thermo::DisTypeToOptGaussRule<distype>::rule);
  if (intpoints.ip().nquad != nquad_) FOUR_C_THROW("Trouble with number of Gauss points");

  // --------------------------------------------- loop over Gauss Points
  for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
  {
    // compute inverse Jacobian matrix and derivatives at GP w.r.t. material
    // coordinates
    eval_shape_func_and_derivs_at_int_point(intpoints, iquad, ele->id());

    // ------------------------------------------------- thermal gradient
    // gradient of current temperature value
    // Grad T = d T_j / d x_i = L . N . T = B_ij T_j
    gradtemp_.multiply_nn(derxy_, etempn_);

    // ---------------------------------------- call thermal material law
    // call material law => cmat_,heatflux_ and dercmat_
    // negative q is used for balance equation:
    // heatflux_ = k_0 . Grad T
    materialize(ele, iquad);
    // heatflux_ := qintermediate = k_0 . Grad T

    // -------------------------------------------- coupling to mechanics
    // (material) deformation gradient F
    // F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.multiply_tt(xcurr, derxy_);
    // inverse of deformation gradient
    invdefgrd.invert(defgrd);

    // Inverse right Cauchy-Green tensor C^{-1} = F^{-1} F^{-T}.
    Core::LinAlg::Matrix<nsd_, nsd_> Cinv(Core::LinAlg::Initialization::uninitialized);
    Cinv.multiply_nt(invdefgrd, invdefgrd);

    // initial heatflux Q = C^{-1} . qintermediate = k_0 . C^{-1} . B_T . T
    // the current heatflux q = detF . F^{-1} . q
    // store heatflux
    // (3x1)  (3x3) . (3x1)
    Core::LinAlg::Matrix<nsd_, 1> initialheatflux(Core::LinAlg::Initialization::uninitialized);
    initialheatflux.multiply(Cinv, heatflux_);
    // put the initial, material heatflux onto heatflux_
    heatflux_.update(initialheatflux);
    // from here on heatflux_ == -Q

    // ------------------------------ integrate internal force vector r_T
    // add the displacement-dependent terms to fint
    // fint = fint + fint_{Td}
    if (efint != nullptr)
    {
      // fint += B_T^T . Q . detJ * w(gp)
      //      += B_T^T . (k_0) . C^{-1} . B_T . T . detJ . w(gp)
      // (8x1)   (8x3) (3x1)
      efint->multiply_tn(fac_, derxy_, heatflux_, 1.0);
    }  // (efint != nullptr)

    // ------------------------------- integrate conductivity matrix k_TT
    // update conductivity matrix k_TT (with displacement dependent term)
    if (econd != nullptr)
    {
      // k^e_TT += ( B_T^T . C^{-1} . C_mat . B_T ) . detJ . w(gp)
      // 3D:        (8x3)    (3x3)    (3x3)   (3x8)
      // with C_mat = k_0 . I
      // -q = C_mat . C^{-1} . B
      Core::LinAlg::Matrix<nsd_, nen_> aop(Core::LinAlg::Initialization::uninitialized);  // (3x8)
      aop.multiply_nn(cmat_, derxy_);  // (nsd_xnsd_)(nsd_xnen_)
      Core::LinAlg::Matrix<nsd_, nen_> aop1(Core::LinAlg::Initialization::uninitialized);  // (3x8)
      aop1.multiply_nn(Cinv, aop);  // (nsd_xnsd_)(nsd_xnen_)

      // k^e_TT += ( B_T^T . C^{-1} . C_mat . B_T ) . detJ . w(gp)
      econd->multiply_tn(fac_, derxy_, aop1, 1.0);  //(8x8)=(8x3)(3x8)

      // linearization of non-constant conductivity
      // k^e_TT += ( B_T^T . C^{-1} . dC_mat . B_T . T . N) . detJ . w(gp)
      Core::LinAlg::Matrix<nsd_, 1> dCmatGradT(Core::LinAlg::Initialization::uninitialized);
      dCmatGradT.multiply_nn(dercmat_, gradtemp_);
      Core::LinAlg::Matrix<nsd_, 1> CinvdCmatGradT(Core::LinAlg::Initialization::uninitialized);
      CinvdCmatGradT.multiply_nn(Cinv, dCmatGradT);
      Core::LinAlg::Matrix<nsd_, nen_> CinvdCmatGradTN(Core::LinAlg::Initialization::uninitialized);
      CinvdCmatGradTN.multiply_nt(CinvdCmatGradT, funct_);
      econd->multiply_tn(fac_, derxy_, CinvdCmatGradTN, 1.0);  //(8x8)=(8x3)(3x8)
    }  // (econd != nullptr)

    // --------------------------------------- capacity matrix m_capa
    // capacity matrix is independent of deformation
    // m_capa corresponds to the mass matrix of the structural field
    if (ecapa != nullptr)
    {
      // m_capa = m_capa + ( N_T^T .  (rho_0 . C_V) . N_T ) . detJ . w(gp)
      //           (8x8)     (8x1)                 (1x8)
      // caution: funct_ implemented as (8,1)--> use transposed in code for
      // theoretic part
      ecapa->multiply_nt((fac_ * capacoeff_), funct_, funct_, 1.0);
    }  // (ecapa != nullptr)
    if (ecapalin != nullptr)
    {
      // calculate additional linearization d(C(T))/dT (3-tensor!)
      // multiply with temperatures to obtain 2-tensor
      //
      // ecapalin = dC/dT*(T_{n+1} -T_{n})
      //          = fac . dercapa . (T_{n+1} -T_{n}) . (N . N^T . T)^T
      Core::LinAlg::Matrix<1, 1> Netemp(Core::LinAlg::Initialization::uninitialized);
      Core::LinAlg::Matrix<numdofpernode_ * nen_, 1> difftemp(
          Core::LinAlg::Initialization::uninitialized);
      Core::LinAlg::Matrix<numdofpernode_ * nen_, 1> NNetemp(
          Core::LinAlg::Initialization::uninitialized);
      // T_{n+1} - T_{n}
      difftemp.update(1.0, etempn_, -1.0, etemp_);
      Netemp.multiply_tn(funct_, difftemp);
      NNetemp.multiply_nn(funct_, Netemp);
      ecapalin->multiply_nt((fac_ * dercapa_), NNetemp, funct_, 1.0);
    }

  }  // ---------------------------------- end loop over Gauss Points
}

template <Core::FE::CellType distype>
void Discret::Elements::TemperImpl<distype>::nonlinear_coupled_tang(
    const Core::Elements::Element* ele,  // the element whose matrix is calculated
    const std::vector<double>& disp,     // current displacements
    const std::vector<double>& vel,      // current velocities
    Core::LinAlg::Matrix<nen_ * numdofpernode_, nsd_ * nen_ * numdofpernode_>* etangcoupl,
    Teuchos::ParameterList& params  // parameter list
)
{
  // Derivative of the conduction on the deformed configuration w.r.t. the displacements. The
  // derivative of the mechanical heat source is evaluated by the structural elements.
  if (etangcoupl == nullptr) return;

  // update element geometry
  Core::LinAlg::Matrix<nen_, nsd_> xcurr(
      Core::LinAlg::Initialization::uninitialized);  // current  coord. of element
  Core::LinAlg::Matrix<nen_, nsd_> xcurrrate(
      Core::LinAlg::Initialization::uninitialized);  // current  velocity of element
  initial_and_current_nodal_position_velocity(ele, disp, vel, xcurr, xcurrrate);

  // the thermal internal force is weighted by the time integrator
  const double timefac = std::invoke(
      [&]()
      {
        const auto timint =
            params.get<Thermo::DynamicType>("time integrator", Thermo::DynamicType::Undefined);
        switch (timint)
        {
          case Thermo::DynamicType::Statics:
            return 1.0;
          case Thermo::DynamicType::OneStepTheta:
            return params.get<double>("theta");
          case Thermo::DynamicType::GenAlpha:
            return params.get<double>("alphaf");
          case Thermo::DynamicType::Undefined:
          default:
            FOUR_C_THROW("Add correct temporal coefficient here!");
        }
      });

  Core::LinAlg::Matrix<nsd_, nsd_> defgrd(Core::LinAlg::Initialization::uninitialized);
  Core::LinAlg::Matrix<nsd_, nsd_> invdefgrd(Core::LinAlg::Initialization::uninitialized);

  // integrations points and weights
  Core::FE::IntPointsAndWeights<nsd_> intpoints(Thermo::DisTypeToOptGaussRule<distype>::rule);
  if (intpoints.ip().nquad != nquad_) FOUR_C_THROW("Trouble with number of Gauss points");

  for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
  {
    eval_shape_func_and_derivs_at_int_point(intpoints, iquad, ele->id());

    // Grad T and the conductivity at the current state
    gradtemp_.multiply_nn(derxy_, etempn_);
    materialize(ele, iquad);

    defgrd.multiply_tt(xcurr, derxy_);
    invdefgrd.invert(defgrd);
    Core::LinAlg::Matrix<nsd_, nsd_> Cinv(Core::LinAlg::Initialization::uninitialized);
    Cinv.multiply_nt(invdefgrd, invdefgrd);

    // conduction: r_i = dN_i/dX . C^{-1} . G with G = C_mat . Grad T
    // dC^{-1}/dd_{nk} = - F^{-1} . (e_k (x) dN_n/dX) . C^{-1} - C^{-1} . (dN_n/dX (x) e_k) . F^{-T}
    // dr_i/dd_{nk} = - (dN_i/dX . F^{-1} . e_k) (dN_n/dX . C^{-1} . G)
    //                - (dN_i/dX . C^{-1} . dN_n/dX) (G . F^{-1} . e_k)
    Core::LinAlg::Matrix<nsd_, 1> G(Core::LinAlg::Initialization::uninitialized);
    G.multiply(cmat_, gradtemp_);
    Core::LinAlg::Matrix<nsd_, 1> CinvG(Core::LinAlg::Initialization::uninitialized);
    CinvG.multiply(Cinv, G);
    Core::LinAlg::Matrix<nen_, nsd_> dN_Finv(Core::LinAlg::Initialization::uninitialized);
    dN_Finv.multiply_tn(derxy_, invdefgrd);  // (dN_i/dX . F^{-1})_k
    Core::LinAlg::Matrix<nen_, 1> dN_CinvG(Core::LinAlg::Initialization::uninitialized);
    dN_CinvG.multiply_tn(derxy_, CinvG);  // dN_n/dX . C^{-1} . G
    Core::LinAlg::Matrix<nsd_, nen_> Cinv_dN(Core::LinAlg::Initialization::uninitialized);
    Cinv_dN.multiply(Cinv, derxy_);
    Core::LinAlg::Matrix<nen_, nen_> dN_Cinv_dN(Core::LinAlg::Initialization::uninitialized);
    dN_Cinv_dN.multiply_tn(derxy_, Cinv_dN);  // dN_i/dX . C^{-1} . dN_n/dX
    Core::LinAlg::Matrix<nsd_, 1> Finv_T_G(Core::LinAlg::Initialization::uninitialized);
    Finv_T_G.multiply_tn(invdefgrd, G);  // (G . F^{-1})_k

    for (int i = 0; i < nen_; ++i)
    {
      for (int n = 0; n < nen_; ++n)
      {
        for (int k = 0; k < nsd_; ++k)
        {
          (*etangcoupl)(i, n* nsd_ + k) -=
              fac_ * (dN_Finv(i, k) * dN_CinvG(n) + dN_Cinv_dN(i, n) * Finv_T_G(k));
        }
      }
    }
  }

  // scale total tangent with timefac
  etangcoupl->scale(timefac);
}

template <Core::FE::CellType distype>
void Discret::Elements::TemperImpl<distype>::linear_heatflux_tempgrad(
    const Core::Elements::Element* ele,
    Core::LinAlg::Matrix<nquad_, nsd_>* eheatflux,  // heat fluxes at Gauss points
    Core::LinAlg::Matrix<nquad_, nsd_>* etempgrad   // temperature gradients at Gauss points
)
{
  Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, nen_>>(
      ele, xyze_);

  // integrations points and weights
  Core::FE::IntPointsAndWeights<nsd_> intpoints(Thermo::DisTypeToOptGaussRule<distype>::rule);
  if (intpoints.ip().nquad != nquad_) FOUR_C_THROW("Trouble with number of Gauss points");

  // ----------------------------------------- loop over Gauss Points
  for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
  {
    eval_shape_func_and_derivs_at_int_point(intpoints, iquad, ele->id());

    // gradient of current temperature value
    // grad T = d T_j / d x_i = L . N . T = B_ij T_j
    gradtemp_.multiply_nn(derxy_, etempn_);

    // store the temperature gradient for postprocessing
    if (etempgrad != nullptr)
      for (int idim = 0; idim < nsd_; ++idim)
        // (8x3)                    (3x1)
        (*etempgrad)(iquad, idim) = gradtemp_(idim);

    // call material law => cmat_,heatflux_
    // negative q is used for balance equation: -q = -(-k gradtemp)= k * gradtemp
    materialize(ele, iquad);

    // store the heat flux for postprocessing
    if (eheatflux != nullptr)
      // negative sign for heat flux introduced here
      for (int idim = 0; idim < nsd_; ++idim) (*eheatflux)(iquad, idim) = -heatflux_(idim);
  }
}

template <Core::FE::CellType distype>
void Discret::Elements::TemperImpl<distype>::nonlinear_heatflux_tempgrad(
    const Core::Elements::Element* ele,             // the element whose matrix is calculated
    const std::vector<double>& disp,                // current displacements
    const std::vector<double>& vel,                 // current velocities
    Core::LinAlg::Matrix<nquad_, nsd_>* eheatflux,  // heat fluxes at Gauss points
    Core::LinAlg::Matrix<nquad_, nsd_>* etempgrad,  // temperature gradients at Gauss points
    Teuchos::ParameterList& params)
{
  // specific choice of heat flux / temperature gradient
  const auto ioheatflux =
      params.get<Thermo::HeatFluxType>("ioheatflux", Thermo::HeatFluxType::None);
  const auto iotempgrad =
      params.get<Thermo::TempGradType>("iotempgrad", Thermo::TempGradType::None);

  // update element geometry
  Core::LinAlg::Matrix<nen_, nsd_> xcurr;      // current  coord. of element
  Core::LinAlg::Matrix<nen_, nsd_> xcurrrate;  // current  coord. of element
  initial_and_current_nodal_position_velocity(ele, disp, vel, xcurr, xcurrrate);

  // build the deformation gradient w.r.t. material configuration
  Core::LinAlg::Matrix<nsd_, nsd_> defgrd(Core::LinAlg::Initialization::uninitialized);
  // inverse of deformation gradient
  Core::LinAlg::Matrix<nsd_, nsd_> invdefgrd(Core::LinAlg::Initialization::uninitialized);

  // ----------------------------------- integration loop for one element
  Core::FE::IntPointsAndWeights<nsd_> intpoints(Thermo::DisTypeToOptGaussRule<distype>::rule);
  if (intpoints.ip().nquad != nquad_) FOUR_C_THROW("Trouble with number of Gauss points");

  // --------------------------------------------- loop over Gauss Points
  for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
  {
    // compute inverse Jacobian matrix and derivatives at GP w.r.t. material
    // coordinates
    eval_shape_func_and_derivs_at_int_point(intpoints, iquad, ele->id());

    gradtemp_.multiply_nn(derxy_, etempn_);

    // ---------------------------------------- call thermal material law
    // call material law => cmat_,heatflux_ and dercmat_
    // negative q is used for balance equation:
    // heatflux_ = k_0 . Grad T
    materialize(ele, iquad);
    // heatflux_ := qintermediate = k_0 . Grad T

    // -------------------------------------------- coupling to mechanics
    // (material) deformation gradient F
    // F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.multiply_tt(xcurr, derxy_);
    // inverse of deformation gradient
    invdefgrd.invert(defgrd);

    Core::LinAlg::Matrix<nsd_, nsd_> Cinv(Core::LinAlg::Initialization::uninitialized);
    // build the inverse of the right Cauchy-Green deformation gradient C^{-1}
    // C^{-1} = F^{-1} . F^{-T}
    Cinv.multiply_nt(invdefgrd, invdefgrd);

    switch (iotempgrad)
    {
      case Thermo::TempGradType::Initial:
      {
        if (etempgrad == nullptr) FOUR_C_THROW("tempgrad data not available");
        // etempgrad = Grad T
        for (int idim = 0; idim < nsd_; ++idim) (*etempgrad)(iquad, idim) = gradtemp_(idim);
        break;
      }
      case Thermo::TempGradType::Current:
      {
        if (etempgrad == nullptr) FOUR_C_THROW("tempgrad data not available");
        // etempgrad = grad T = Grad T . F^{-1} =  F^{-T} . Grad T
        // (8x3)        (3x1)   (3x1)    (3x3)     (3x3)    (3x1)
        // spatial temperature gradient
        Core::LinAlg::Matrix<nsd_, 1> currentgradT(Core::LinAlg::Initialization::uninitialized);
        currentgradT.multiply_tn(invdefgrd, gradtemp_);
        for (int idim = 0; idim < nsd_; ++idim) (*etempgrad)(iquad, idim) = currentgradT(idim);
        break;
      }
      case Thermo::TempGradType::None:
      {
        // no postprocessing of temperature gradients
        break;
      }
      default:
        FOUR_C_THROW("requested tempgrad type not available");
        break;
    }  // iotempgrad

    switch (ioheatflux)
    {
      case Thermo::HeatFluxType::Initial:
      {
        if (eheatflux == nullptr) FOUR_C_THROW("heat flux data not available");
        Core::LinAlg::Matrix<nsd_, 1> initialheatflux(Core::LinAlg::Initialization::uninitialized);
        // eheatflux := Q = -k_0 . Cinv . Grad T
        initialheatflux.multiply(Cinv, heatflux_);
        for (int idim = 0; idim < nsd_; ++idim) (*eheatflux)(iquad, idim) = -initialheatflux(idim);
        break;
      }
      case Thermo::HeatFluxType::Current:
      {
        if (eheatflux == nullptr) FOUR_C_THROW("heat flux data not available");
        // eheatflux := q = - k_0 . 1/(detF) . F^{-T} . Grad T
        // (8x3)     (3x1)            (3x3)  (3x1)
        const double detF = defgrd.determinant();
        Core::LinAlg::Matrix<nsd_, 1> spatialq;
        spatialq.multiply_tn((1.0 / detF), invdefgrd, heatflux_);
        for (int idim = 0; idim < nsd_; ++idim) (*eheatflux)(iquad, idim) = -spatialq(idim);
        break;
      }
      case Thermo::HeatFluxType::None:
      {
        // no postprocessing of heat fluxes, continue!
        break;
      }
      default:
        FOUR_C_THROW("requested heat flux type not available");
        break;
    }  // ioheatflux
  }
}


template <Core::FE::CellType distype>
void Discret::Elements::TemperImpl<distype>::extract_disp_vel(
    const Core::FE::Discretization& discretization, const Core::Elements::LocationArray& la,
    std::vector<double>& mydisp, std::vector<double>& myvel) const
{
  if ((discretization.has_state(1, "displacement")) and (discretization.has_state(1, "velocity")))
  {
    // get the displacements
    std::shared_ptr<const Core::LinAlg::Vector<double>> disp =
        discretization.get_state(1, "displacement");
    if (disp == nullptr) FOUR_C_THROW("Cannot get state vectors 'displacement'");
    // extract the displacements
    mydisp = Core::FE::extract_values(*disp, la[1].lm_);

    // get the velocities
    std::shared_ptr<const Core::LinAlg::Vector<double>> vel =
        discretization.get_state(1, "velocity");
    if (vel == nullptr) FOUR_C_THROW("Cannot get state vectors 'velocity'");
    // extract the displacements
    myvel = Core::FE::extract_values(*vel, la[1].lm_);
  }
}

template <Core::FE::CellType distype>
void Discret::Elements::TemperImpl<distype>::calculate_lump_matrix(
    Core::LinAlg::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* ecapa) const
{
  // lump capacity matrix
  if (ecapa != nullptr)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned int c = 0; c < (*ecapa).n(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned int r = 0; r < (*ecapa).m(); ++r)  // parse rows
      {
        d += (*ecapa)(r, c);  // accumulate row entries
        (*ecapa)(r, c) = 0.0;
      }
      (*ecapa)(c, c) = d;  // apply sum of row entries on diagonal
    }
  }
}

template <Core::FE::CellType distype>
void Discret::Elements::TemperImpl<distype>::radiation(
    const Core::Elements::Element* ele, const double time)
{
  std::vector<const Core::Conditions::Condition*> myneumcond;

  // check whether all nodes have a unique VolumeNeumann condition
  switch (nsd_)
  {
    case 3:
      Core::Conditions::find_element_conditions(ele, "VolumeNeumann", myneumcond);
      break;
    case 2:
      Core::Conditions::find_element_conditions(ele, "SurfaceNeumann", myneumcond);
      break;
    case 1:
      Core::Conditions::find_element_conditions(ele, "LineNeumann", myneumcond);
      break;
    default:
      FOUR_C_THROW("Illegal number of space dimensions: {}", nsd_);
      break;
  }

  if (myneumcond.size() > 1) FOUR_C_THROW("more than one VolumeNeumann cond on one node");

  if (myneumcond.size() == 1)
  {
    // get node coordinates
    Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, nen_>>(
        ele, xyze_);

    // update element geometry
    Core::LinAlg::Matrix<nen_, nsd_> xrefe;  // material coord. of element
    auto nodes = ele->nodes();
    for (int i = 0; i < nen_; ++i)
    {
      const auto& x = nodes[i]->x();
      // (8x3) = (nen_xnsd_)
      for (int j = 0; j < nsd_; j++) xrefe(i, j) = x[j];
    }


    // integrations points and weights
    Core::FE::IntPointsAndWeights<nsd_> intpoints(Thermo::DisTypeToOptGaussRule<distype>::rule);
    if (intpoints.ip().nquad != nquad_) FOUR_C_THROW("Trouble with number of Gauss points");

    radiation_.clear();

    // compute the Jacobian matrix
    Core::LinAlg::Matrix<nsd_, nsd_> jac;
    jac.multiply(derxy_, xrefe);

    // compute determinant of Jacobian
    const double detJ = jac.determinant();
    if (detJ == 0.0)
      FOUR_C_THROW("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0)
      FOUR_C_THROW("NEGATIVE JACOBIAN DETERMINANT");

    const auto funct = myneumcond[0]->parameters().get<std::vector<std::optional<int>>>("FUNCT");

    Core::LinAlg::Matrix<nsd_, 1> xrefegp(Core::LinAlg::Initialization::uninitialized);
    // material/reference co-ordinates of Gauss point
    for (int dim = 0; dim < nsd_; dim++)
    {
      xrefegp(dim) = 0.0;
      for (int nodid = 0; nodid < nen_; ++nodid) xrefegp(dim) += funct_(nodid) * xrefe(nodid, dim);
    }

    // function evaluation
    FOUR_C_ASSERT(funct.size() == 1, "Need exactly one function.");

    double functfac = 1.0;
    if (funct[0].has_value() && funct[0].value() > 0)
      // evaluate function at current gauss point (3D position vector required!)
      functfac = Global::Problem::instance()
                     ->function_by_id<Core::Utils::FunctionOfSpaceTime>(funct[0].value())
                     .evaluate(xrefegp.as_span(), time, 0);

    // get values and switches from the condition
    const auto onoff = myneumcond[0]->parameters().get<std::vector<int>>("ONOFF");
    const auto val = myneumcond[0]->parameters().get<std::vector<double>>("VAL");

    // set this condition to the radiation array
    for (int idof = 0; idof < numdofpernode_; idof++)
    {
      radiation_(idof) = onoff[idof] * val[idof] * functfac;
    }
  }
  else
  {
    radiation_.clear();
  }
}


template <Core::FE::CellType distype>
void Discret::Elements::TemperImpl<distype>::materialize(
    const Core::Elements::Element* ele, const int gp)
{
  auto material = ele->material();

  // calculate the current temperature at the integration point
  Core::LinAlg::Matrix<1, 1> temp;
  temp.multiply_tn(1.0, funct_, etempn_, 0.0);

  auto thermoMaterial = std::dynamic_pointer_cast<Mat::Trait::Thermo>(material);
  thermoMaterial->reinit(temp(0), gp);
  thermoMaterial->evaluate(gradtemp_, cmat_, heatflux_, ele->id());
  capacoeff_ = thermoMaterial->capacity();
  thermoMaterial->conductivity_deriv_t(dercmat_);
  dercapa_ = thermoMaterial->capacity_deriv_t();
}

template <Core::FE::CellType distype>
void Discret::Elements::TemperImpl<distype>::eval_shape_func_and_derivs_at_int_point(
    const Core::FE::IntPointsAndWeights<nsd_>& intpoints,  // integration points
    const int iquad,                                       // id of current Gauss point
    const int eleid                                        // the element id
)
{
  // coordinates of the current (Gauss) integration point (xsi_)
  const double* gpcoord = (intpoints.ip().qxg)[iquad];
  for (int idim = 0; idim < nsd_; idim++)
  {
    xsi_(idim) = gpcoord[idim];
  }

  // shape functions (funct_) and their first derivatives (deriv_)
  // N, N_{,xsi}
  if (myknots_.size() == 0)
  {
    Core::FE::shape_function<distype>(xsi_, funct_);
    Core::FE::shape_function_deriv1<distype>(xsi_, deriv_);
  }
  else
    Core::FE::Nurbs::nurbs_get_3d_funct_deriv(funct_, deriv_, xsi_, myknots_, weights_, distype);

  // compute Jacobian matrix and determinant (as presented in FE lecture notes)
  // actually compute its transpose (compared to J in NiliFEM lecture notes)
  // J = dN/dxsi . x^{-}
  /*
   *   J-NiliFEM               J-FE
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
   */

  // derivatives at gp w.r.t. material coordinates (N_XYZ in solid)
  xjm_.multiply_nt(deriv_, xyze_);
  // xij_ = J^{-T}
  // det = J^{-T} *
  // J = (N_rst * X)^T (6.24 NiliFEM)
  const double det = xij_.invert(xjm_);

  if (det < 1e-16)
    FOUR_C_THROW("GLOBAL ELEMENT NO.{}\nZERO OR NEGATIVE JACOBIAN DETERMINANT: {}", eleid, det);

  // set integration factor: fac = Gauss weight * det(J)
  fac_ = intpoints.ip().qwgt[iquad] * det;

  // compute global derivatives
  derxy_.multiply(xij_, deriv_);
}

template <Core::FE::CellType distype>
void Discret::Elements::TemperImpl<distype>::initial_and_current_nodal_position_velocity(
    const Core::Elements::Element* ele, const std::vector<double>& disp,
    const std::vector<double>& vel, Core::LinAlg::Matrix<nen_, nsd_>& xcurr,
    Core::LinAlg::Matrix<nen_, nsd_>& xcurrrate)
{
  Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, nen_>>(
      ele, xyze_);
  for (int i = 0; i < nen_; ++i)
  {
    for (int j = 0; j < nsd_; ++j)
    {
      xcurr(i, j) = xyze_(j, i) + disp[i * nsd_ + j];
      xcurrrate(i, j) = vel[i * nsd_ + j];
    }
  }
}

template <Core::FE::CellType distype>
void Discret::Elements::TemperImpl<distype>::prepare_nurbs_eval(
    const Core::Elements::Element* ele,             // the element whose matrix is calculated
    const Core::FE::Discretization& discretization  // current discretisation
)
{
  if (ele->shape() != Core::FE::CellType::nurbs27)
  {
    myknots_.resize(0);
    return;
  }

  myknots_.resize(3);  // fixme: dimension
                       // get nurbs specific infos
  // cast to nurbs discretization
  const auto* nurbsdis =
      dynamic_cast<const Core::FE::Nurbs::NurbsDiscretization*>(&(discretization));
  if (nurbsdis == nullptr) FOUR_C_THROW("So_nurbs27 appeared in non-nurbs discretisation\n");

  // zero-sized element
  if ((*((*nurbsdis).get_knot_vector())).get_ele_knots(myknots_, ele->id())) return;

  // get weights from cp's
  for (int inode = 0; inode < nen_; inode++)
    weights_(inode) = dynamic_cast<const Core::FE::Nurbs::ControlPoint*>(ele->nodes()[inode])->w();
}

template <Core::FE::CellType distype>
void Discret::Elements::TemperImpl<distype>::integrate_shape_functions(
    const Core::Elements::Element* ele, Core::LinAlg::SerialDenseVector& elevec1,
    const Core::LinAlg::IntSerialDenseVector& dofids)
{
  // get node coordinates
  Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, nen_>>(
      ele, xyze_);

  // integrations points and weights
  Core::FE::IntPointsAndWeights<nsd_> intpoints(Thermo::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.ip().nquad; gpid++)
  {
    eval_shape_func_and_derivs_at_int_point(intpoints, gpid, ele->id());

    // compute integral of shape functions (only for dofid)
    for (int k = 0; k < numdofpernode_; k++)
    {
      if (dofids[k] >= 0)
      {
        for (int node = 0; node < nen_; node++)
        {
          elevec1[node * numdofpernode_ + k] += funct_(node) * fac_;
        }
      }
    }
  }  // loop over integration points

}  // TemperImpl<distype>::integrate_shape_function


template <Core::FE::CellType distype>
void Discret::Elements::TemperImpl<distype>::extrapolate_from_gauss_points_to_nodes(
    const Core::Elements::Element* ele,  // the element whose matrix is calculated
    const Core::LinAlg::Matrix<nquad_, nsd_>& gpheatflux,
    Core::LinAlg::Matrix<nen_ * numdofpernode_, 1>& efluxx,
    Core::LinAlg::Matrix<nen_ * numdofpernode_, 1>& efluxy,
    Core::LinAlg::Matrix<nen_ * numdofpernode_, 1>& efluxz)
{
  // this quick'n'dirty hack functions only for elements which has the same
  // number of gauss points AND same number of nodes
  if (not((distype == Core::FE::CellType::hex8) or (distype == Core::FE::CellType::hex27) or
          (distype == Core::FE::CellType::tet4) or (distype == Core::FE::CellType::quad4) or
          (distype == Core::FE::CellType::line2)))
    FOUR_C_THROW("Sorry, not implemented for element shape");

  // another check
  if (nen_ * numdofpernode_ != nquad_)
    FOUR_C_THROW("Works only if number of gauss points and nodes match");

  // integrations points and weights
  Core::FE::IntPointsAndWeights<nsd_> intpoints(Thermo::DisTypeToOptGaussRule<distype>::rule);
  if (intpoints.ip().nquad != nquad_) FOUR_C_THROW("Trouble with number of Gauss points");

  // build matrix of shape functions at Gauss points
  Core::LinAlg::Matrix<nquad_, nquad_> shpfctatgps;
  for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
  {
    // coordinates of the current integration point
    const double* gpcoord = (intpoints.ip().qxg)[iquad];
    for (int idim = 0; idim < nsd_; idim++) xsi_(idim) = gpcoord[idim];

    // shape functions and their first derivatives
    Core::FE::shape_function<distype>(xsi_, funct_);

    for (int inode = 0; inode < nen_; ++inode) shpfctatgps(iquad, inode) = funct_(inode);
  }

  // extrapolation
  Core::LinAlg::Matrix<nquad_, nsd_> ndheatflux;  //  objective nodal heatflux
  Core::LinAlg::Matrix<nquad_, nsd_> gpheatflux2(
      gpheatflux);  // copy the heatflux at the Gauss point
  {
    Core::LinAlg::FixedSizeSerialDenseSolver<nquad_, nquad_, nsd_> solver;  // must be quadratic
    solver.set_matrix(shpfctatgps);
    solver.set_vectors(ndheatflux, gpheatflux2);
    solver.solve();
  }

  // copy into component vectors
  for (int idof = 0; idof < nen_ * numdofpernode_; ++idof)
  {
    efluxx(idof) = ndheatflux(idof, 0);
    if (nsd_ > 1) efluxy(idof) = ndheatflux(idof, 1);
    if (nsd_ > 2) efluxz(idof) = ndheatflux(idof, 2);
  }
}

template <Core::FE::CellType distype>
double Discret::Elements::TemperImpl<distype>::calculate_char_ele_length() const
{
  // volume of the element (2D: element surface area; 1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double vol = fac_;

  // as shown in calc_char_ele_length() in ScaTraImpl
  // c) cubic/square root of element volume/area or element length (3-/2-/1-D)
  // cast dimension to a double variable -> pow()

  // get characteristic element length as cubic root of element volume
  // (2D: square root of element area, 1D: element length)
  // h = vol^(1/dim)
  double h = std::pow(vol, (1.0 / nsd_));

  return h;
}


template <Core::FE::CellType distype>
void Discret::Elements::TemperImpl<distype>::compute_error(
    const Core::Elements::Element* ele,  // the element whose matrix is calculated
    Core::LinAlg::Matrix<nen_ * numdofpernode_, 1>& elevec1,
    Teuchos::ParameterList& params  // parameter list
)
{
  // get node coordinates
  Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, nen_>>(
      ele, xyze_);

  // get scalar-valued element temperature
  // build the product of the shapefunctions and element temperatures T = N . T
  Core::LinAlg::Matrix<1, 1> NT(Core::LinAlg::Initialization::uninitialized);

  // analytical solution
  Core::LinAlg::Matrix<1, 1> T_analytical(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<1, 1> deltaT(Core::LinAlg::Initialization::zero);
  // ------------------------------- integration loop for one element

  // integrations points and weights
  Core::FE::IntPointsAndWeights<nsd_> intpoints(Thermo::DisTypeToOptGaussRule<distype>::rule);
  //  if (intpoints.ip().nquad != nquad_)
  //    FOUR_C_THROW("Trouble with number of Gauss points");

  const auto calcerr = Teuchos::getIntegralValue<Thermo::CalcError>(params, "calculate error");
  const int errorfunctno = params.get<int>("error function number");
  const double t = params.get<double>("total time");

  // ----------------------------------------- loop over Gauss Points
  for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
  {
    // compute inverse Jacobian matrix and derivatives
    eval_shape_func_and_derivs_at_int_point(intpoints, iquad, ele->id());

    // ------------------------------------------------ thermal terms

    // gradient of current temperature value
    // grad T = d T_j / d x_i = L . N . T = B_ij T_j
    gradtemp_.multiply_nn(derxy_, etempn_);

    // current element temperatures
    // N_T . T (funct_ defined as <nen,1>)
    NT.multiply_tn(funct_, etempn_);  // (1x8)(8x1)

    // H1 -error norm
    // compute first derivative of the displacement
    Core::LinAlg::Matrix<nsd_, 1> derT(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<nsd_, 1> deltaderT(Core::LinAlg::Initialization::zero);

    // Compute analytical solution
    switch (calcerr)
    {
      case Thermo::calcerror_byfunct:
      {
        // get coordinates at integration point
        // gp reference coordinates
        Core::LinAlg::Matrix<nsd_, 1> xyzint(Core::LinAlg::Initialization::zero);
        xyzint.multiply(xyze_, funct_);

        // function evaluation requires a 3D position vector!!
        double position[3] = {0.0, 0.0, 0.0};

        for (int dim = 0; dim < nsd_; ++dim) position[dim] = xyzint(dim);

        const double T_exact = Global::Problem::instance()
                                   ->function_by_id<Core::Utils::FunctionOfSpaceTime>(errorfunctno)
                                   .evaluate(position, t, 0);

        T_analytical(0, 0) = T_exact;

        std::vector<double> Tder_exact =
            Global::Problem::instance()
                ->function_by_id<Core::Utils::FunctionOfSpaceTime>(errorfunctno)
                .evaluate_spatial_derivative(position, t, 0);

        if (Tder_exact.size())
        {
          for (int dim = 0; dim < nsd_; ++dim) derT(dim) = Tder_exact[dim];
        }
      }
      break;
      default:
        FOUR_C_THROW("analytical solution is not defined");
        break;
    }

    // compute difference between analytical solution and numerical solution
    deltaT.update(1.0, NT, -1.0, T_analytical);

    // H1 -error norm
    // compute error for first velocity derivative
    deltaderT.update(1.0, gradtemp_, -1.0, derT);

    // 0: delta temperature for L2-error norm
    // 1: delta temperature for H1-error norm
    // 2: analytical temperature for L2 norm
    // 3: analytical temperature for H1 norm

    // the error for the L2 and H1 norms are evaluated at the Gauss point

    // integrate delta velocity for L2-error norm
    elevec1(0) += deltaT(0, 0) * deltaT(0, 0) * fac_;
    // integrate delta velocity for H1-error norm
    elevec1(1) += deltaT(0, 0) * deltaT(0, 0) * fac_;
    // integrate analytical velocity for L2 norm
    elevec1(2) += T_analytical(0, 0) * T_analytical(0, 0) * fac_;
    // integrate analytical velocity for H1 norm
    elevec1(3) += T_analytical(0, 0) * T_analytical(0, 0) * fac_;

    // integrate delta velocity derivative for H1-error norm
    elevec1(1) += deltaderT.dot(deltaderT) * fac_;
    // integrate analytical velocity for H1 norm
    elevec1(3) += derT.dot(derT) * fac_;
  }
}

FOUR_C_NAMESPACE_CLOSE
