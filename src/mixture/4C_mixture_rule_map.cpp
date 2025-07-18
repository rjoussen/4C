// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mixture_rule_map.hpp"

#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_mixture_constituent.hpp"
#include "4C_utils_exceptions.hpp"

#include <algorithm>
#include <iosfwd>
#include <memory>
#include <string>
#include <unordered_map>


FOUR_C_NAMESPACE_OPEN

namespace
{
  std::vector<double> get_validate_mass_fractions(
      const Core::IO::InputField<std::vector<double>>& mass_fractions_field, const int ele_id_key,
      const std::size_t num_constituents)
  {
    const std::vector<double> massfracs = mass_fractions_field.at(ele_id_key);
    if (massfracs.size() != num_constituents)
    {
      FOUR_C_THROW(
          "Number of mass fractions for element id {} does not match the number of constituents "
          "{}.",
          ele_id_key, num_constituents);
    }

    // check, whether the mass frac sums up to 1
    const double sum = std::accumulate(massfracs.begin(), massfracs.end(), 0.0);
    if (std::abs(1.0 - sum) > 1e-8)
      FOUR_C_THROW(
          "Mass fractions for element id {} don't sum up to 1, which is unphysical.", ele_id_key);

    return massfracs;
  }
}  // namespace

Mixture::PAR::MapMixtureRule::MapMixtureRule(const Core::Mat::PAR::Parameter::Data& matdata)
    : MixtureRule(matdata),
      initial_reference_density_(matdata.parameters.get<double>("DENS")),
      num_constituents_(matdata.parameters.get<int>("NUMCONST")),
      mass_fractions_field_(matdata.parameters.get<Core::IO::InputField<std::vector<double>>>(
          "MASSFRACMAPFILE_CONTENT")) {};

std::unique_ptr<Mixture::MixtureRule> Mixture::PAR::MapMixtureRule::create_rule()
{
  return std::make_unique<Mixture::MapMixtureRule>(this);
}

Mixture::MapMixtureRule::MapMixtureRule(Mixture::PAR::MapMixtureRule* params)
    : MixtureRule(params), params_(params)
{
}

void Mixture::MapMixtureRule::setup(const Teuchos::ParameterList& params, const int eleGID)
{
  MixtureRule::setup(params, eleGID);
}

void Mixture::MapMixtureRule::unpack_mixture_rule(Core::Communication::UnpackBuffer& buffer)
{
  Mixture::MixtureRule::unpack_mixture_rule(buffer);
}

void Mixture::MapMixtureRule::evaluate(const Core::LinAlg::Tensor<double, 3, 3>& F,
    const Core::LinAlg::SymmetricTensor<double, 3, 3>& E_strain,
    const Teuchos::ParameterList& params, Core::LinAlg::SymmetricTensor<double, 3, 3>& S_stress,
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat, const int gp, const int eleGID)
{
  // define temporary matrices
  Core::LinAlg::SymmetricTensor<double, 3, 3> cstress;
  Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> ccmat;

  // evaluate the mass fractions at the given element id
  auto massfracs =
      get_validate_mass_fractions(params_->mass_fractions_field_, eleGID, constituents().size());

  // Iterate over all constituents and add all stress/cmat contributions
  for (std::size_t i = 0; i < constituents().size(); ++i)
  {
    MixtureConstituent& constituent = *constituents()[i];
    cstress = {};
    ccmat = {};
    constituent.evaluate(F, E_strain, params, cstress, ccmat, gp, eleGID);

    // add stress contribution to global stress
    double constituent_density = params_->initial_reference_density_ * massfracs[i];
    S_stress += constituent_density * cstress;
    cmat += constituent_density * ccmat;
  }
}


FOUR_C_NAMESPACE_CLOSE
