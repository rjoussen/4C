// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fsi_nox_mpe.hpp"

#include "4C_global_data.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_vector.hpp"

#include <NOX_Abstract_Group.H>
#include <NOX_Epetra_Group.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_GlobalData.H>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

NOX::FSI::MinimalPolynomial::MinimalPolynomial(
    const Teuchos::RCP<::NOX::Utils>& utils, Teuchos::ParameterList& params)
    : utils_(utils)
{
  Teuchos::ParameterList& mpeparams = params.sublist("Extrapolation");
  kmax_ = mpeparams.get("kmax", 10);
  omega_ = mpeparams.get("omega", 0.01);
  eps_ = mpeparams.get("Tolerance", 1e-1);
  mpe_ = mpeparams.get("Method", "RRE") == "MPE";
}



bool NOX::FSI::MinimalPolynomial::reset(
    const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& params)
{
  utils_ = gd->getUtils();
  return true;
}


bool NOX::FSI::MinimalPolynomial::compute(
    ::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& grp, const ::NOX::Solver::Generic& solver)
{
  // We work in a local copy of the group so that we do not spoil the
  // current state.
  Teuchos::RCP<::NOX::Abstract::Group> copy_group_ptr = grp.clone();
  ::NOX::Abstract::Group& group = *copy_group_ptr;

  const ::NOX::Abstract::Vector& x = group.getX();

  std::vector<Teuchos::RCP<::NOX::Epetra::Vector>> q;
  Core::LinAlg::SerialDenseMatrix r(kmax_ + 1, kmax_ + 1, true);
  Core::LinAlg::SerialDenseVector c(kmax_ + 1, true);
  Core::LinAlg::SerialDenseVector gamma(kmax_ + 1, true);

  int k;
  for (k = 0; k < kmax_; ++k)
  {
    Teuchos::Time t("", true);
    ::NOX::Abstract::Group::ReturnType status;

    // Compute F at current solution
    status = group.computeF();
    if (status != ::NOX::Abstract::Group::Ok) throw_error("compute", "Unable to compute F");

    // get f = d(k+1) - d(k)
    const ::NOX::Epetra::Vector& f = dynamic_cast<const ::NOX::Epetra::Vector&>(group.getF());

    // We have to work on the scaled residual here.
    Teuchos::RCP<::NOX::Epetra::Vector> y = Teuchos::make_rcp<::NOX::Epetra::Vector>(f);
    y->scale(omega_);

    // modified Gram-Schmidt
    for (int j = 0; j < k; ++j)
    {
      r(j, k) = y->innerProduct(*q[j]);
      y->update(-r(j, k), *q[j], 1.);
    }
    r(k, k) = sqrt(y->innerProduct(*y));

    // store new direction
    if (r(k, k) > 1e-32 * r(0, 0) and k < kmax_)
    {
      y->scale(1. / r(k, k));
      q.push_back(y);
    }
    else if (r(k, k) <= 1e-32 * r(0, 0))
    {
      if (utils_->isPrintType(::NOX::Utils::Error))
        utils_->err() << "r(" << k << "," << k << ") <= " << 1e-32 << "*r(0,0)\n";
      break;
    }

    double res = 0;

    if (mpe_)
    {
      // MPE gamma calculation
      for (int i = k - 1; i >= 0; --i)
      {
        double ci = -r(i, k);
        for (int j = i + 1; j < k; ++j)
        {
          ci -= r(i, j) * c(j);
        }
        c(i) = ci / r(i, i);
      }
      c(k) = 1.;

      // sum over all entries of c
      double sc = 0.0;
      for (int l = 0; l < c.length(); l++) sc += c(l);

      if (fabs(sc) < 1e-16)
      {
        throw_error("compute", "sum(c) equals zero");
      }

      Core::LinAlg::update(1 / sc, c, 0.0, gamma);
      res = r(k, k) * fabs(gamma(k));

      if (utils_->isPrintType(::NOX::Utils::InnerIteration))
      {
        utils_->out() << "MPE:  k=" << k << "  res=" << res << "  eps*r(0,0)=" << eps_ * r(0, 0)
                      << "  r(k,k)=" << r(k, k) << std::endl;
      }
    }
    else
    {
      // RRE gamma calculation
      for (int i = 0; i <= k; ++i)
      {
        double ci = 1.;
        for (int j = 0; j < i; ++j)
        {
          ci -= r(j, i) * c(j);
        }
        c(i) = ci / r(i, i);
      }
      for (int i = k; i >= 0; --i)
      {
        double ci = c(i);
        for (int j = i + 1; j <= k; ++j)
        {
          ci -= r(i, j) * gamma(j);
        }
        gamma(i) = ci / r(i, i);
      }

      // sum over all entries of gamma
      double sc = 0.0;
      for (int l = 0; l < gamma.length(); l++) sc += gamma(l);

      gamma.scale(1 / sc);
      res = 1. / sqrt(fabs(sc));

      if (utils_->isPrintType(::NOX::Utils::InnerIteration))
      {
        utils_->out() << "RRE:  k=" << std::setw(2) << k << "  res=" << std::scientific << res
                      << "  eps*r(0,0)=" << std::scientific << eps_ * r(0, 0)
                      << "  r(k,k)=" << std::scientific << r(k, k) << "  time=" << std::scientific
                      << t.totalElapsedTime(true) << std::endl;
      }
    }

    // leave if we are close enough to the solution
    if (res <= eps_ * r(0, 0) or r(k, k) <= 1e-32 * r(0, 0))
    {
      k += 1;
      break;
    }

    // Update the group to go another round
    // Note: We do not use any extrapolated vector here but simply go
    // on with the series of vectors. The fixed relaxation is needed
    // to keep the iteration from diverging.
    group.computeX(group, f, omega_);
  }

  // calc extrapolated vector
  Core::LinAlg::SerialDenseVector xi(kmax_, true);
  xi(0) = 1. - gamma(0);
  for (int j = 1; j < k; ++j)
  {
    xi(j) = xi(j - 1) - gamma(j);
  }

  ::NOX::Epetra::Vector s(dynamic_cast<const ::NOX::Epetra::Vector&>(x));
  for (int j = 0; j < k; ++j)
  {
    double hp = 0.;
    for (int i = j; i < k; ++i)
    {
      hp += r(j, i) * xi(i);
    }
    s.update(hp, *q[j], 1.);
  }

  // set direction from original position
  dir.update(1., s, -1., x, 0.);

  return true;
}


bool NOX::FSI::MinimalPolynomial::compute(::NOX::Abstract::Vector& dir,
    ::NOX::Abstract::Group& group, const ::NOX::Solver::LineSearchBased& solver)
{
  return ::NOX::Direction::Generic::compute(dir, group, solver);
}


void NOX::FSI::MinimalPolynomial::throw_error(
    const std::string& functionName, const std::string& errorMsg)
{
  if (utils_->isPrintType(::NOX::Utils::Error))
    utils_->err() << "MinimalPolynomial::" << functionName << " - " << errorMsg << std::endl;
  throw "NOX Error";
}

FOUR_C_NAMESPACE_CLOSE
