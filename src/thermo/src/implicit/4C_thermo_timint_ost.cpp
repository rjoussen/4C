// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_thermo_timint_ost.hpp"

#include "4C_thermo_ele_action.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                               dano 06/13 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntOneStepTheta::verify_coeff()
{
  // theta
  if ((theta_ <= 0.0) or (theta_ > 1.0)) FOUR_C_THROW("theta out of range (0.0,1.0]");

  // done
  return;

}  // VerifyCoeff()


/*----------------------------------------------------------------------*
 | constructor                                               dano 08/09 |
 *----------------------------------------------------------------------*/
Thermo::TimIntOneStepTheta::TimIntOneStepTheta(const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& tdynparams, const Teuchos::ParameterList& xparams,
    std::shared_ptr<Core::FE::Discretization> actdis, std::shared_ptr<Core::LinAlg::Solver> solver,
    std::shared_ptr<Core::IO::DiscretizationWriter> output)
    : TimIntImpl(ioparams, tdynparams, xparams, actdis, solver, output),
      theta_(tdynparams.sublist("ONESTEPTHETA").get<double>("THETA")),
      tempt_(nullptr),
      fint_(nullptr),
      fintn_(nullptr),
      fcap_(nullptr),
      fcapn_(nullptr),
      fext_(nullptr),
      fextn_(nullptr)
{
  // info to user
  if (Core::Communication::my_mpi_rank(discret_->get_comm()) == 0)
  {
    // check if coefficient has admissible value
    verify_coeff();

    // print values of time integration parameters to screen
    std::cout << "with one-step-theta" << std::endl
              << "   theta = " << theta_ << std::endl
              << std::endl;
  }

  // determine capacity
  determine_capa_consist_temp_rate();

  // create state vectors
  // mid-temperatures
  tempt_ = Core::LinAlg::create_vector(*discret_->dof_row_map(), true);

  // create force vectors
  // internal force vector F_{int;n} at last time
  fint_ = Core::LinAlg::create_vector(*discret_->dof_row_map(), true);
  // internal force vector F_{int;n+1} at new time
  fintn_ = Core::LinAlg::create_vector(*discret_->dof_row_map(), true);
  // stored force vector F_{transient;n} at last time
  fcap_ = Core::LinAlg::create_vector(*discret_->dof_row_map(), true);
  // stored force vector F_{transient;n+1} at new time
  fcapn_ = Core::LinAlg::create_vector(*discret_->dof_row_map(), true);
  // set initial internal force vector
  apply_force_tang_internal(time_[0], dt_[0], temp_(0), zeros_, fcap_, fint_, tang_);

  // external force vector F_ext at last times
  fext_ = Core::LinAlg::create_vector(*discret_->dof_row_map(), true);
  // external force vector F_{n+1} at new time
  fextn_ = Core::LinAlg::create_vector(*discret_->dof_row_map(), true);
  // set initial external force vector
  apply_force_external(time_[0], temp_(0), *fext_);
  // set initial external force vector of convective heat transfer boundary
  // conditions
  apply_force_external_conv(time_[0], temp_(0), temp_(0), fext_, tang_);

  // have a nice day
  return;

}  // TimIntOneStepTheta()


/*----------------------------------------------------------------------*
 | consistent predictor with constant temperatures           dano 08/09 |
 | and consistent temperature rates and temperatures                    |
 *----------------------------------------------------------------------*/
void Thermo::TimIntOneStepTheta::predict_const_temp_consist_rate()
{
  // time step size
  const double dt = dt_[0];

  // constant predictor : temperature in domain
  tempn_->update(1.0, *temp_(0), 0.0);

  // new end-point temperature rates
  // R_{n+1}^{i+1} = -(1 - theta)/theta . R_n + 1/(theta . dt) . (T_{n+1}^{i+1} - T_n)
  raten_->update(1.0, *tempn_, -1.0, *temp_(0), 0.0);
  raten_->update(-(1.0 - theta_) / theta_, *rate_(0), 1.0 / (theta_ * dt));

  // watch out
  return;

}  // predict_const_temp_consist_rate()


/*----------------------------------------------------------------------*
 | evaluate residual force and its tangent, ie derivative    dano 08/09 |
 | with respect to end-point temperatures \f$T_{n+1}\f$                 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntOneStepTheta::evaluate_rhs_tang_residual()
{
  // theta-interpolate state vectors
  evaluate_mid_state();

  // build new external forces
  fextn_->put_scalar(0.0);

  // initialise tangent matrix to zero
  tang_->zero();

  // set initial external force vector of convective heat transfer boundary
  // conditions

  // if the boundary condition shall be dependent on the current temperature
  // solution T_n+1 --> linearisation must be uncommented
  // --> use tempn_

  // if the old temperature T_n  is sufficient --> no linearisation needed!
  // --> use temp_(0)
  apply_force_external_conv(timen_, temp_(0), tempn_, fextn_, tang_);

  apply_force_external(timen_, temp_(0), *fextn_);

  // initialise internal forces
  fintn_->put_scalar(0.0);
  fcapn_->put_scalar(0.0);

  // ordinary internal force and tangent
  apply_force_tang_internal(timen_, dt_[0], tempn_, tempi_, fcapn_, fintn_, tang_);

  // build residual  Res = R_{n+theta}
  //                     + F_{int;n+theta}
  //                     - F_{ext;n+theta}
  // with R_{n+theta}     = M_cap . ( T_{n+1} - T_n ) / dt = fcapn_ - fcap_
  //      F_{int;n+theta} = theta * F_{int;n+1} + (1 - theta) * F_{int;n}
  //      F_{ext;n+theta} = - theta * F_{ext;n+1} - (1 - theta) * F_{ext;n}

  // here the time derivative is introduced needed for fcap depending on T'!
  fres_->scale(1.0, *fcapn_);  // fcap already contains full R_{n+theta}
  fres_->update(theta_, *fintn_, (1.0 - theta_), *fint_, 1.0);
  // here is the negative sign for the external forces (heatfluxes)
  fres_->update(-theta_, *fextn_, -(1.0 - theta_), *fext_, 1.0);

  // no further modification on tang_ required
  // tang_ is already effective dynamic tangent matrix
  tang_->complete();  // close tangent matrix

  // hallelujah
  return;

}  // evaluate_rhs_tang_residual()


/*----------------------------------------------------------------------*
 | evaluate theta-state vectors by averaging                 dano 08/09 |
 | end-point vector                                                     |
 *----------------------------------------------------------------------*/
void Thermo::TimIntOneStepTheta::evaluate_mid_state()
{
  // mid-temperatures T_{n+1-alpha_f} (tempm)
  //    T_{n+theta} := theta * T_{n+1} + (1-theta) * T_{n}
  tempt_->update(theta_, *tempn_, 1.0 - theta_, *temp_(0), 0.0);

  // jump
  return;

}  // EvaluateMidState()


/*----------------------------------------------------------------------*
 | calculate characteristic/reference norms for              dano 08/09 |
 | temperatures originally by lw                                        |
 *----------------------------------------------------------------------*/
double Thermo::TimIntOneStepTheta::calc_ref_norm_temperature()
{
  // The reference norms are used to scale the calculated iterative
  // temperature norm and/or the residual force norm. For this
  // purpose we only need the right order of magnitude, so we don't
  // mind evaluating the corresponding norms at possibly different
  // points within the timestep (end point, generalized midpoint).

  double charnormtemp = 0.0;
  charnormtemp = Thermo::Aux::calculate_vector_norm(iternorm_, *temp_(0));

  // rise your hat
  return charnormtemp;

}  // calc_ref_norm_temperature()


/*----------------------------------------------------------------------*
 | calculate characteristic/reference norms for forces       dano 08/09 |
 | originally by lw                                                     |
 *----------------------------------------------------------------------*/
double Thermo::TimIntOneStepTheta::calc_ref_norm_force()
{
  // The reference norms are used to scale the calculated iterative
  // temperature norm and/or the residual force norm. For this
  // purpose we only need the right order of magnitude, so we don't
  // mind evaluating the corresponding norms at possibly different
  // points within the timestep (end point, generalized midpoint).

  // norm of the internal forces
  double fintnorm = 0.0;
  fintnorm = Thermo::Aux::calculate_vector_norm(iternorm_, *fintn_);

  // norm of the external forces
  double fextnorm = 0.0;
  fextnorm = Thermo::Aux::calculate_vector_norm(iternorm_, *fextn_);

  // norm of reaction forces
  double freactnorm = 0.0;
  freactnorm = Thermo::Aux::calculate_vector_norm(iternorm_, *freact_);

  // norm of stored forces
  double fcapnorm = 0.0;
  fcapnorm = Thermo::Aux::calculate_vector_norm(iternorm_, *fcap_);

  // return char norm
  return std::max(fcapnorm, std::max(fintnorm, std::max(fextnorm, freactnorm)));

}  // calc_ref_norm_force()


/*----------------------------------------------------------------------*
 | incremental iteration update of state                     dano 08/09 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntOneStepTheta::update_iter_incrementally()
{
  // Auxiliary vector holding new temperature rates
  // by extrapolation/scheme on __all__ DOFs. This includes
  // the Dirichlet DOFs as well. Thus we need to protect those
  // DOFs of overwriting; they already hold the
  // correctly 'predicted', final values.
  std::shared_ptr<Core::LinAlg::Vector<double>> aux =
      Core::LinAlg::create_vector(*discret_->dof_row_map(), false);

  // new end-point temperatures
  // T_{n+1}^{i+1} := T_{n+1}^{<k>} + IncT_{n+1}^{i}
  tempn_->update(1.0, *tempi_, 1.0);

  // new end-point temperature rates
  // aux = - (1-theta)/theta R_n + 1/(theta . dt) (T_{n+1}^{i+1} - T_{n+1}^i)
  aux->update(1.0, *tempn_, -1.0, *temp_(0), 0.0);
  aux->update(-(1.0 - theta_) / theta_, *rate_(0), 1.0 / (theta_ * dt_[0]));
  // put only to free/non-DBC DOFs
  dbcmaps_->insert_other_vector(*dbcmaps_->extract_other_vector(*aux), *raten_);

  // bye
  return;

}  // update_iter_incrementally()


/*----------------------------------------------------------------------*
 | iterative iteration update of state                       dano 08/09 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntOneStepTheta::update_iter_iteratively()
{
  // new end-point temperatures
  // T_{n+1}^{<k+1>} := T_{n+1}^{<k>} + IncT_{n+1}^{<k>}
  tempn_->update(1.0, *tempi_, 1.0);

  // new end-point temperature rates
  // R_{n+1}^{<k+1>} := R_{n+1}^{<k>} + 1/(theta . dt)IncT_{n+1}^{<k>}
  raten_->update(1.0 / (theta_ * dt_[0]), *tempi_, 1.0);

  // bye
  return;

}  // update_iter_iteratively()


/*----------------------------------------------------------------------*
 | update after time step                                    dano 08/09 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntOneStepTheta::update_step_state()
{
  // update state
  // new temperatures at t_{n+1} -> t_n
  //    T_{n} := T_{n+1}
  temp_.update_steps(*tempn_);
  // new temperature rates at t_{n+1} -> t_n
  //    R_{n} := R_{n+1}
  rate_.update_steps(*raten_);

  // update new external force
  //    F_{ext;n} := F_{ext;n+1}
  fext_->update(1.0, *fextn_, 0.0);

  // update new internal force
  //    F_{int;n} := F_{int;n+1}
  fint_->update(1.0, *fintn_, 0.0);

  // update new stored transient force
  //    F_{cap;n} := F_{cap;n+1}
  fcap_->update(1.0, *fcapn_, 0.0);


  // look out
  return;

}  // update_step_state()


/*----------------------------------------------------------------------*
 | update after time step after output on element level      dano 05/13 |
 | update anything that needs to be updated at the element level        |
 *----------------------------------------------------------------------*/
void Thermo::TimIntOneStepTheta::update_step_element()
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // other parameters that might be needed by the elements
  p.set("total time", timen_);
  p.set("delta time", dt_[0]);
  // action for elements
  p.set<Thermo::Action>("action", Thermo::calc_thermo_update_istep);
  // go to elements
  discret_->set_state(0, "temperature", *tempn_);
  discret_->evaluate(p, nullptr, nullptr, nullptr, nullptr, nullptr);

}  // update_step_element()


/*----------------------------------------------------------------------*
 | read restart forces                                       dano 08/09 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntOneStepTheta::read_restart_force()
{
  Core::IO::DiscretizationReader reader(
      discret_, Global::Problem::instance()->input_control_file(), step_);
  reader.read_vector(fext_, "fexternal");
  reader.read_vector(fint_, "fint");
  reader.read_vector(fcap_, "fcap");

  return;

}  // ReadRestartForce()


/*----------------------------------------------------------------------*
 | write internal and external forces for restart            dano 07/13 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntOneStepTheta::write_restart_force(
    std::shared_ptr<Core::IO::DiscretizationWriter> output)
{
  output->write_vector("fexternal", fext_);
  output->write_vector("fint", fint_);
  output->write_vector("fcap", fcap_);

  return;

}  // WriteRestartForce()


/*----------------------------------------------------------------------*
 | evaluate the internal force and the tangent               dano 08/09 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntOneStepTheta::apply_force_tang_internal(const double time,  //!< evaluation time
    const double dt,                                                           //!< step size
    const std::shared_ptr<Core::LinAlg::Vector<double>> temp,   //!< temperature state
    const std::shared_ptr<Core::LinAlg::Vector<double>> tempi,  //!< residual temperatures
    std::shared_ptr<Core::LinAlg::Vector<double>> fcap,         //!< capacity force
    std::shared_ptr<Core::LinAlg::Vector<double>> fint,         //!< internal force
    std::shared_ptr<Core::LinAlg::SparseMatrix> tang            //!< tangent matrix
)
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // set parameters
  p.set<double>("theta", theta_);
  p.set<double>("timefac", theta_);
  p.set<bool>("lump capa matrix", lumpcapa_);
  // call the base function
  TimInt::apply_force_tang_internal(p, time, dt, temp, tempi, fcap, fint, tang);
  // finish
  return;

}  // apply_force_tang_internal()


/*----------------------------------------------------------------------*
 | evaluate the internal force                               dano 08/09 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntOneStepTheta::apply_force_internal(const double time,  //!< evaluation time
    const double dt,                                                      //!< step size
    const std::shared_ptr<Core::LinAlg::Vector<double>> temp,             //!< temperature state
    const std::shared_ptr<Core::LinAlg::Vector<double>> tempi,  //!< incremental temperatures
    std::shared_ptr<Core::LinAlg::Vector<double>> fint          //!< internal force
)
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // set parameters
  p.set("theta", theta_);
  // call the base function
  TimInt::apply_force_internal(p, time, dt, temp, tempi, fint);
  // finish
  return;

}  // apply_force_tang_internal()


/*----------------------------------------------------------------------*
 | evaluate the convective boundary condition                dano 12/10 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntOneStepTheta::apply_force_external_conv(const double time,  //!< evaluation time
    const std::shared_ptr<Core::LinAlg::Vector<double>> tempn,  //!< old temperature state T_n
    const std::shared_ptr<Core::LinAlg::Vector<double>> temp,   //!< temperature state T_n+1
    std::shared_ptr<Core::LinAlg::Vector<double>> fext,         //!< external force
    std::shared_ptr<Core::LinAlg::SparseMatrix> tang            //!< tangent matrix
)
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // set parameters
  p.set<double>("theta", theta_);
  // call the base function
  TimInt::apply_force_external_conv(p, time, tempn, temp, fext, tang);
  // finish
  return;

}  // apply_force_external_conv()


/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
