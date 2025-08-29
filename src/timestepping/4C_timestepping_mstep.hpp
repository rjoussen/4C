// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_TIMESTEPPING_MSTEP_HPP
#define FOUR_C_TIMESTEPPING_MSTEP_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_linalg_map.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_shared_ptr_from_ref.hpp"

#include <memory>
#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/* A collection of time integration utilities */
namespace TimeStepping
{
  /*====================================================================*/
  /*!
   * \brief This is the base object for holding multi-step solution quantities,
   *        e.g. displacements \f$D_n\f$, \f$D_{n-1}\f$, \f$D_{n-2}\f$,
   *        \f$\ldots\f$, or time points \f$t_n\f$, \f$t_{n-1}\f$,
   *        \f$t_{n-2}\f$, \f$\ldots\f$, etc.
   *
   * Multi-step quantities occur --surprisingly indeed-- in multi-step time
   * integrators. However, they are needed in single-step integrators
   * as well, if the auxiliary scheme is multi-step. The idea is to be
   * able to resize the multi-step quantities dynamically and have
   * intrinsic update mechanism such that a single-step integrator
   * can provide the data on which a multi-step auxiliary integrator
   * can work.
   *

   */
  template <typename STATE>
  class TimIntMStepBase
  {
   public:
    //! @name Access functions
    //@{

    //! Access state object by time step index
    STATE& operator[](const int step  //!< inquiry step
    )
    {
      if (not step_exists(step)) FOUR_C_THROW("Step {} is not admissible", step);
      return state_[index_by_step(step)];
    }

    //! Access state object by time step index (read-only)
    const STATE& operator[](const int step  //!< inquiry step
    ) const
    {
      if (not step_exists(step)) FOUR_C_THROW("Step {} is not admissible", step);
      return state_[index_by_step(step)];
    }

    //! Access state object by time step index and wrap into std::shared_ptr
    std::shared_ptr<STATE> operator()(const int step  //!< inquiry step
    )
    {
      if (not step_exists(step)) FOUR_C_THROW("Step {} is not admissible", step);
      return Core::Utils::shared_ptr_from_ref(state_[index_by_step(step)]);
    }

    std::shared_ptr<const STATE> operator()(const int step  //!< inquiry step
    ) const
    {
      if (not step_exists(step)) FOUR_C_THROW("Step {} is not admissible", step);
      return Core::Utils::shared_ptr_from_ref(state_[index_by_step(step)]);
    }

    //@}

    //! @name Query functions
    //@{

    //! step index of most past and future steps
    std::pair<int, int> get_steps() { return std::make_pair(steppast_, stepfuture_); }

    //! Number of steps
    int num_steps() { return steps_; }
    //@}

   protected:
    //! @name Allocation
    //@{

    //! Dummy constructor
    TimIntMStepBase() : steppast_(0), stepfuture_(0), steps_(0), state_() { return; }

    /*! Another dummy constructor which sets vector limits and its size
     *  but does not allocate anything
     */
    TimIntMStepBase(const int steppast,  //!< lower index bound
        const int stepfuture             //!< higher index bound, >= lower bound
        )
        : steppast_(steppast),
          stepfuture_(stepfuture),
          steps_(stepfuture >= steppast ? stepfuture - steppast + 1 : 0),
          state_()
    {
      // verify a positive #steps_
      if (steps_ <= 0) FOUR_C_THROW("Past step must be lower or equal to future step");

      return;
    }

    //@}

    //! @name Verification
    //@{

    //! Check sanity prior resize
    bool resize_sane(const int steppast,  //!< lower index bound
        const int stepfuture              //!< higher index bound, >= lower bound
    )
    {
      bool sane = true;

      sane = sane and (steppast <= stepfuture);
      if (not sane) FOUR_C_THROW("Past step must be lower than future step");

      sane = sane and (stepfuture == stepfuture_);
      if (not sane) FOUR_C_THROW("Future step cannot be changed");

      return sane;
    }

    //! Determine whether step lies in given bounds
    bool step_exists(const int step  //!< inquired step
    ) const
    {
      return ((step >= steppast_) and (step <= stepfuture_));
    }

    //! determine whether index lies in given bounds
    bool index_exists(const int index  //!< inquired index
    ) const
    {
      return ((index >= 0) and (index < (int)state_.size()));
    }

    //! Map step number to vector index
    int step_by_index(const int index  //!< vector index
    )
    {
      return index + steppast_;
    }

    //! Map vector index to step
    int index_by_step(const int step  //!< step number
    ) const
    {
      int index = step - steppast_;
      FOUR_C_ASSERT(index_exists(index), "step is not permissible!");
      return (unsigned)index;
    }

    //@}

    //! @name Steps
    //@{

    int steppast_;    ///< lowest number
    int stepfuture_;  ///< highest number
    int steps_;       ///< total

    //@}

    //! Multi-step quantity, it's stored vectorially
    std::vector<STATE> state_;
  };

  /*====================================================================*/
  /*!
   * \brief General derived class for multi-step quantities of
   *        simple type.
   *        This is going to be used with \c double and \c int
   */
  template <typename STATE>
  class TimIntMStep : public TimIntMStepBase<STATE>
  {
   protected:
    //! template base class
    using MStepBase = TimIntMStepBase<STATE>;

   public:
    //! @name Life
    //@{

    TimIntMStep() = default;

    //! Constructor
    TimIntMStep(const int steppast,  //!< lower index bound
        const int stepfuture,        //!< higher index bound, >= lower bound
        const STATE init             //!< initialise to init
        )
        : MStepBase(steppast, stepfuture)
    {
      // allocate the vectors themselves
      for (int index = 0; index < MStepBase::steps_; ++index)
      {
        MStepBase::state_.push_back(STATE(init));
      }

      return;
    }

    //! Resize
    void resize(const int steppast,  //!< lower index bound
        const int stepfuture,        //!< higher index bound, >= lower bound
        const STATE init             //!< initialise to init
    )
    {
      // check this
      {
        bool sane = MStepBase::resize_sane(steppast, stepfuture);
        if (not sane) FOUR_C_THROW("Sanity check not passed.");
      }

      // add states for steps in past
      if (steppast < MStepBase::steppast_)
      {
        for (int past = MStepBase::steppast_; past > steppast; --past)
        {
          MStepBase::state_.insert(MStepBase::state_.begin(), STATE(init));
        }
        MStepBase::steppast_ = steppast;
        MStepBase::steps_ = MStepBase::stepfuture_ - MStepBase::steppast_ + 1;
      }

      return;
    }

    //@}

    //! @name Actions
    //@{

    //! Set entry at #step to #value
    void set_step(const int step, const STATE value)
    {
      MStepBase::state_.at(MStepBase::index_by_step(step)) = value;
    }

    /*! \brief Update multi-step state
     *
     *  i.e. state_{n} := state_{n+1},
     *       state_{n-1} := state_{n},
     *       state_{n} := state_{n+1},
     *       etc.
     */
    void update_steps(const STATE staten  //!< state_{n+1}
    )
    {
      for (int ind = 0; ind < MStepBase::steps_ - 1; ++ind)
      {
        MStepBase::state_[ind] = MStepBase::state_[ind + 1];
      }
      MStepBase::state_[MStepBase::steps_ - 1] = staten;

      return;
    }

    //@}
  };

  /*====================================================================*/
  /*! \brief Specialize general #TimIntMStepBase object
   *         for \c Core::LinAlg::Vector<double> as needed for state vectors
   *         like displacements, velocities and accelerations
   */
  template <>
  class TimIntMStep<Core::LinAlg::Vector<double>>
      : public TimIntMStepBase<Core::LinAlg::Vector<double>>
  {
   protected:
    //! base template class
    using MStepBase = TimIntMStepBase<Core::LinAlg::Vector<double>>;

   public:
    //! @name Life
    //@{

    //! Dummy constructor
    TimIntMStep() : MStepBase::TimIntMStepBase() { ; }

    //! Constructor
    TimIntMStep(const int steppast,          //!< lower index bound
        const int stepfuture,                //!< higher index bound, >= lower bound
        const Core::LinAlg::Map* dofrowmap,  //!< vector layout from discretization
        const bool inittozero                //!< initialise to zero, if true
        )
        : MStepBase(steppast, stepfuture)
    {
      // allocate the vectors themselves
      for (int index = 0; index < steps_; ++index)
      {
        state_.push_back(Core::LinAlg::Vector<double>(*dofrowmap, inittozero));
      }

      return;
    }

    /*! \brief Resize
     *
     *  State vectors are added and placed according to their
     *  indices #steppast to #stepfuture
     */
    void resize(const int steppast,          //!< lower index bound
        const int stepfuture,                //!< higher index bound, >= lower bound
        const Core::LinAlg::Map* dofrowmap,  //!< vector layout from discretization
        const bool inittozero                //!< initialise to zero, if true
    )
    {
      // check this
      {
        bool sane = MStepBase::resize_sane(steppast, stepfuture);
        if (not sane) FOUR_C_THROW("Sanity check not passed.");
      }

      // add states for steps in past
      if (steppast < steppast_)
      {
        for (int past = MStepBase::steppast_; past > steppast; --past)
        {
          MStepBase::state_.insert(
              MStepBase::state_.begin(), Core::LinAlg::Vector<double>(*dofrowmap, inittozero));
        }
        MStepBase::steppast_ = steppast;
        MStepBase::steps_ = MStepBase::stepfuture_ - MStepBase::steppast_ + 1;
      }

      return;
    }

    //@}

    //! @name Actions
    //@{

    /*! \brief Update multi-step state
     *
     *  i.e. state_{n} := state_{n+1},
     *       state_{n-1} := state_{n},
     *       state_{n} := state_{n+1},
     *       etc.
     */
    void update_steps(const Core::LinAlg::Vector<double>& staten  //!< state_{n+1}
    )
    {
      for (int ind = 0; ind < MStepBase::steps_ - 1; ++ind)
      {
        (MStepBase::state_[ind]).update(1.0, (MStepBase::state_[ind + 1]), 0.0);
      }
      (MStepBase::state_[steps_ - 1]).update(1.0, staten, 0.0);

      return;
    }

    //@}

  };  // class TimIntMStep<Core::LinAlg::Vector<double>>
}  // namespace TimeStepping

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
