// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_TIMINT_EXPLICIT_HPP
#define FOUR_C_STRUCTURE_NEW_TIMINT_EXPLICIT_HPP


#include "4C_config.hpp"

#include "4C_structure_new_expl_generic.hpp"
#include "4C_structure_new_nln_solver_generic.hpp"
#include "4C_structure_new_timint_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Solid
{
  namespace TimeInt
  {
    /** \brief Explicit time integration strategy
     *
     * */
    class Explicit : public Base
    {
     public:
      //! constructor
      Explicit();

      void setup() override;

      void integrate() override;

      void integrate_step() override;

      void prepare_time_step() override;

      void update_state_incrementally(
          std::shared_ptr<const Core::LinAlg::Vector<double>> disiterinc) override;

      void determine_stress_strain() override;

      void evaluate() override;

      void evaluate(std::shared_ptr<const Core::LinAlg::Vector<double>> disiterinc) override;

      void set_state(const std::shared_ptr<Core::LinAlg::Vector<double>>& x) override;

      void reset_step() override;

      [[nodiscard]] Solid::StepStatus solve() override;

      void prepare_partition_step() override;

      void update(double endtime) override;

      void print_step() override;

      Solid::StcScale get_stc_algo() override;

      std::shared_ptr<Core::LinAlg::SparseMatrix> get_stc_mat() override;

      std::shared_ptr<const Core::LinAlg::Vector<double>> initial_guess() override;

      std::shared_ptr<const Core::LinAlg::Vector<double>> get_f() const override;

      std::shared_ptr<Core::LinAlg::Vector<double>> freact() override;

      std::shared_ptr<Core::LinAlg::SparseMatrix> system_matrix() override;

      std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> block_system_matrix() override;

      void use_block_matrix(std::shared_ptr<const Core::LinAlg::MultiMapExtractor> domainmaps,
          std::shared_ptr<const Core::LinAlg::MultiMapExtractor> rangemaps) override;

      ///@}

      //! @name Attribute access functions
      //@{

      Solid::DynamicType method_name() const override;

      bool is_implicit() const override { return false; }

      bool is_explicit() const override { return true; }

      int method_steps() const override;

      int method_order_of_accuracy_dis() const override;

      int method_order_of_accuracy_vel() const override;

      double method_lin_err_coeff_dis() const override;

      double method_lin_err_coeff_vel() const override;

      //@}

     protected:
      Solid::EXPLICIT::Generic& expl_int()
      {
        check_init_setup();
        return *explint_ptr_;
      };

      Solid::Nln::SOLVER::Generic& nln_solver()
      {
        check_init_setup();
        return *nlnsolver_ptr_;
      };

     private:
      //! ptr to the explicit time integrator object
      std::shared_ptr<Solid::EXPLICIT::Generic> explint_ptr_;

      //! ptr to the non-linear solver object
      std::shared_ptr<Solid::Nln::SOLVER::Generic> nlnsolver_ptr_;
    };
  }  // namespace TimeInt
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
