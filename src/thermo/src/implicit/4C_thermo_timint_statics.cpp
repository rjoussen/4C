// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_thermo_timint_statics.hpp"

#include "4C_thermo_ele_action.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                               dano 08/09 |
 *----------------------------------------------------------------------*/
Thermo::TimIntStatics::TimIntStatics(const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& tdynparams, const Teuchos::ParameterList& xparams,
    std::shared_ptr<Core::FE::Discretization> actdis, std::shared_ptr<Core::LinAlg::Solver> solver,
    std::shared_ptr<Core::IO::DiscretizationWriter> output)
    : TimIntImpl(ioparams, tdynparams, xparams, actdis, solver, output),
      fint_(nullptr),
      fintn_(nullptr),
      fext_(nullptr),
      fextn_(nullptr)
{
  // info to user
  if (Core::Communication::my_mpi_rank(discret_->get_comm()) == 0)
  {
    std::cout << "with statics" << std::endl;
  }

  //! create force vectors

  //! internal force vector F_{int;n} at last time
  fint_ = Core::LinAlg::create_vector(*discret_->dof_row_map(), true);
  //! internal force vector F_{int;n+1} at new time
  fintn_ = Core::LinAlg::create_vector(*discret_->dof_row_map(), true);
  //! set initial internal force vector
  apply_force_tang_internal(time_[0], dt_[0], temp_(0), zeros_, fint_, tang_);
  //! external force vector F_ext at last times
  fext_ = Core::LinAlg::create_vector(*discret_->dof_row_map(), true);
  //! external force vector F_{n+1} at new time
  fextn_ = Core::LinAlg::create_vector(*discret_->dof_row_map(), true);

  // set initial external force vector of convective heat transfer boundary
  // conditions
  apply_force_external_conv(time_[0], temp_(0), temp_(0), fext_, tang_);

  //! set initial external force vector
  apply_force_external(time_[0], temp_(0), *fext_);

  //! have a nice day
  return;
}


/*----------------------------------------------------------------------*
 | consistent predictor with constant temperatures           dano 08/09 |
 | and consistent temperature rates and temperatures                    |
 *----------------------------------------------------------------------*/
void Thermo::TimIntStatics::predict_const_temp_consist_rate()
{
  //! constant predictor : temperature in domain
  // T_n+1,p = T_n
  tempn_->update(1.0, *temp_(0), 0.0);

  //! new end-point temperature rates, these stay zero in static calculation
  raten_->put_scalar(0.0);

  //! watch out
  return;
}  // predict_const_temp_consist_rate()


/*----------------------------------------------------------------------*
 | evaluate residual force and its tangent, ie derivative    dano 08/09 |
 | with respect to end-point temperatures \f$T_{n+1}\f$                 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntStatics::evaluate_rhs_tang_residual()
{
  //! build new external forces
  fextn_->put_scalar(0.0);

  //! initialise tangent matrix to zero
  tang_->zero();

  // set initial external force vector of convective heat transfer boundary
  // conditions
  // Warning: do not use convection boundary condition with T_n and statics
  // --> always use T_n+1 for statics!
  apply_force_external_conv(timen_, temp_(0), tempn_, fextn_, tang_);

  apply_force_external(timen_, temp_(0), *fextn_);

  //! initialise internal forces
  fintn_->put_scalar(0.0);

  //! ordinary internal force and tangent
  apply_force_tang_internal(timen_, dt_[0], tempn_, tempi_, fintn_, tang_);

  //! build residual  Res = F_{int;n+1} - F_{ext;n+1}
  fres_->update(-1.0, *fextn_, 0.0);
  fres_->update(1.0, *fintn_, 1.0);

  //! build tangent matrix : effective dynamic tangent matrix
  //!    K_{Teffdyn} = K_{T}
  //! i.e. do nothing here

  tang_->complete();  // close tangent matrix

  //! hallelujah
  return;
}  // evaluate_rhs_tang_residual()


/*----------------------------------------------------------------------*
 | calculate characteristic/reference norms for              dano 08/09 |
 | temperatures originally by lw                                        |
 *----------------------------------------------------------------------*/
double Thermo::TimIntStatics::calc_ref_norm_temperature()
{
  //! The reference norms are used to scale the calculated iterative
  //! temperature norm and/or the residual force norm. For this
  //! purpose we only need the right order of magnitude, so we don't
  //! mind evaluating the corresponding norms at possibly different
  //! points within the timestep (end point, generalized midpoint).

  double charnormtemp = 0.0;
  charnormtemp = Thermo::Aux::calculate_vector_norm(iternorm_, *temp_(0));

  //! rise your hat
  return charnormtemp;
}  // calc_ref_norm_temperature()


/*----------------------------------------------------------------------*
 | calculate characteristic/reference norms for forces       dano 08/09 |
 | originally by lw                                                     |
 *----------------------------------------------------------------------*/
double Thermo::TimIntStatics::calc_ref_norm_force()
{
  //! The reference norms are used to scale the calculated iterative
  //! temperature norm and/or the residual force norm. For this
  //! purpose we only need the right order of magnitude, so we don't
  //! mind evaluating the corresponding norms at possibly different
  //! points within the timestep (end point, generalized midpoint).

  //! norm of the internal forces
  double fintnorm = 0.0;
  fintnorm = Thermo::Aux::calculate_vector_norm(iternorm_, *fintn_);

  //! norm of the external forces
  double fextnorm = 0.0;
  fextnorm = Thermo::Aux::calculate_vector_norm(iternorm_, *fextn_);

  //! norm of reaction forces
  double freactnorm = 0.0;
  freactnorm = Thermo::Aux::calculate_vector_norm(iternorm_, *freact_);

  //! return char norm
  return std::max(fintnorm, std::max(fextnorm, freactnorm));
}  // calc_ref_norm_force()


/*----------------------------------------------------------------------*
 | incremental iteration update of state                     dano 08/09 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntStatics::update_iter_incrementally()
{
  //! new end-point temperatures
  //! T_{n+1}^{<k+1>} := T_{n+1}^{<k>} + IncT_{n+1}^{<k>}
  tempn_->update(1.0, *tempi_, 1.0);

  //! bye
  return;
}  // update_iter_incrementally()


/*----------------------------------------------------------------------*
 | iterative iteration update of state                       dano 08/09 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntStatics::update_iter_iteratively()
{
  //! new end-point temperatures
  //! T_{n+1}^{<k+1>} := T_{n+1}^{<k>} + IncT_{n+1}^{<k>}
  tempn_->update(1.0, *tempi_, 1.0);

  //! bye
  return;
}  // update_iter_iteratively()


/*----------------------------------------------------------------------*
 | update after time step                                    dano 08/09 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntStatics::update_step_state()
{
  //! update state
  //! new temperatures at t_{n+1} -> t_n
  //!    T_{n} := T_{n+1}
  temp_.update_steps(*tempn_);
  //! new temperature rates at t_{n+1} -> t_n
  //!    T'_{n} := T'_{n+1}
  rate_.update_steps(*raten_);  // this simply copies zero vectors

  //! update new external force
  //!    F_{ext;n} := F_{ext;n+1}
  fext_->update(1.0, *fextn_, 0.0);

  //! update new internal force
  //!    F_{int;n} := F_{int;n+1}
  fint_->update(1.0, *fintn_, 0.0);

  //! look out
  return;
}  // update_step_state()


/*----------------------------------------------------------------------*
 | update after time step after output on element level      dano 05/13 |
 | update anything that needs to be updated at the element level        |
 *----------------------------------------------------------------------*/
void Thermo::TimIntStatics::update_step_element()
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // other parameters that might be needed by the elements
  p.set("total time", timen_);
  p.set("delta time", dt_[0]);
  // action for elements
  p.set<Thermo::Action>("action", Thermo::calc_thermo_update_istep);
  // go to elements
  discret_->evaluate(p, nullptr, nullptr, nullptr, nullptr, nullptr);

}  // update_step_element()


/*----------------------------------------------------------------------*
 | read restart forces                                       dano 08/09 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntStatics::read_restart_force()
{
  // do nothing
  return;

}  // ReadRestartForce()


/*----------------------------------------------------------------------*
 | write internal and external forces for restart            dano 07/13 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntStatics::write_restart_force(
    std::shared_ptr<Core::IO::DiscretizationWriter> output)
{
  // do nothing
  return;

}  // WriteRestartForce()


/*----------------------------------------------------------------------*
 | evaluate the internal force and the tangent               dano 08/09 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntStatics::apply_force_tang_internal(const double time,  //!< evaluation time
    const double dt,                                                      //!< step size
    const std::shared_ptr<Core::LinAlg::Vector<double>> temp,             //!< temperature state
    const std::shared_ptr<Core::LinAlg::Vector<double>> tempi,            //!< residual temperatures
    std::shared_ptr<Core::LinAlg::Vector<double>> fint,                   //!< internal force
    std::shared_ptr<Core::LinAlg::SparseMatrix> tang                      //!< tangent matrix
)
{
  //! create the parameters for the discretization
  Teuchos::ParameterList p;
  //! set parameters
  p.set<double>("timefac", 1.);

  //! call the base function
  TimInt::apply_force_tang_internal(p, time, dt, temp, tempi, fint, tang);
  //! finish
  return;

}  // apply_force_tang_internal()


/*----------------------------------------------------------------------*
 | evaluate the internal force                               dano 08/09 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntStatics::apply_force_internal(const double time,  //!< evaluation time
    const double dt,                                                 //!< step size
    const std::shared_ptr<Core::LinAlg::Vector<double>> temp,        //!< temperature state
    const std::shared_ptr<Core::LinAlg::Vector<double>> tempi,       //!< incremental temperatures
    std::shared_ptr<Core::LinAlg::Vector<double>> fint               //!< internal force
)
{
  //! create the parameters for the discretization
  Teuchos::ParameterList p;
  //! set parameters
  // ...
  //! call the base function
  TimInt::apply_force_internal(p, time, dt, temp, tempi, fint);
  //! finish
  return;

}  // apply_force_tang_internal()


/*----------------------------------------------------------------------*
 | evaluate the convective boundary condition                dano 01/11 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntStatics::apply_force_external_conv(const double time,  //!< evaluation time
    const std::shared_ptr<Core::LinAlg::Vector<double>> tempn,  //!< old temperature state T_n
    const std::shared_ptr<Core::LinAlg::Vector<double>> temp,   //!< temperature state T_n+1
    std::shared_ptr<Core::LinAlg::Vector<double>> fext,         //!< external force
    std::shared_ptr<Core::LinAlg::SparseMatrix> tang            //!< tangent matrix
)
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // set parameters
  // ...
  // call the base function
  TimInt::apply_force_external_conv(p, time, tempn, temp, fext, tang);
  // finish
  return;

}  // apply_force_external_conv()


/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
