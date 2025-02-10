// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef FOUR_C_MAT_INELASTIC_DEFGRAD_FACTORS_SERVICE_HPP
#define FOUR_C_MAT_INELASTIC_DEFGRAD_FACTORS_SERVICE_HPP


#include "4C_config.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_comm_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io_runtime_csv_writer.hpp"
#include "4C_utils_exceptions.hpp"

#include <map>
#include <string>
#include <vector>


FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  /// enum class for error types in InelasticDefgradTransvIsotropElastViscoplast, used for
  /// triggering different procedures (e.g. substepping) during the
  /// Local Newton Loop
  enum class ViscoplastErrorType
  {
    NoErrors,
    NegativePlasticStrain,  // negative plastic strain which does not allow for evaluations
                            // inside the viscoplasticity laws
    OverflowError,  // overflow error of the term \f$ \Delta t \dot{\varepsilon}^{\text{p}} \f$ (and
                    // \f$ \mathsymbol{E}^{\text{p}}  = \exp(- \Delta t \dot{\varepsilon}^{\text{p}}
                    // \mathsymbol{N}^{\text{p}}) \f$) checked in the standard substepping procedure
    NoPlasticIncompressibility,  // no plastic incompressibility, meaning that our determinant
                                 // of the inelastic defgrad is far from
                                 // 1
    FailedSolLinSystLNL,  // solution of the linear system in the Local Newton-Raphson Loop failed
    FailedDetermLineSearchParam,  // the computation of a suitable line search parameter failed
    NoConvergenceLNL,  // the Local Newton Loop did not converge for the given loop settings
    SingularJacobian,  // singular Jacobian after converged LNL, which does not enable our
                       // analytical evaluation of the linearization
    FailedSolAnalytLinearization,  // solution of the linear system in the analytical linearization
                                   // failed
    FailedComputationFlowResistance,  // failed in the computation of the flow resistance via time
                                      // integration of the hardening-rate equation (e.g., Anand
                                      // model)
    FailedComputationFlowResistanceDerivs,  // failed in the computation of the flow resistance
                                            // derivatives (e.g., Anand model)
    FailedLogEval,  // failed evaluation of the matrix logarithm or its derivative
    FailedExpEval,  // failed evaluation of the matrix exponential or its derivative
  };

  /// enum class for error management actions in InelasticDefgradTransvIsotropElastViscoplast
  enum class ViscoplastErrorActions
  {
    Continue,             // continue without any errors (NoErrors)
    ReturnSolWithErrors,  // return the current solution with errors (if the maximum substepping
                          // settings have been reached)
    NextIter,             // go to next iteration after performing certain reset steps
  };

  /// to_string: error types to error messages in InelasticDefgradTransvIsotropElastViscoplast
  inline std::string to_string(Mat::ViscoplastErrorType err_type)
  {
    switch (err_type)
    {
      case Mat::ViscoplastErrorType::NegativePlasticStrain:
        return "Error in InelasticDefgradTransvIsotropElastViscoplast: negative plastic strain!";
      case Mat::ViscoplastErrorType::OverflowError:
        return "Error in InelasticDefgradTransvIsotropElastViscoplast: overflow error related to "
               "the evaluation of the plastic strain increment!";
      case Mat::ViscoplastErrorType::NoPlasticIncompressibility:
        return "Error in InelasticDefgradTransvIsotropElastViscoplast: plastic incompressibility "
               "not satisfied!";
      case Mat::ViscoplastErrorType::FailedSolLinSystLNL:
        return "Error in InelasticDefgradTransvIsotropElastViscoplast: solution of the linear "
               "system in the Local Newton Loop failed!";
      case Mat::ViscoplastErrorType::FailedDetermLineSearchParam:
        return "Error in InelasticDefgradTransvIsotropElastViscoplast: could not determine a "
               "suitable line search parameter!";
      case Mat::ViscoplastErrorType::NoConvergenceLNL:
        return "Error in InelasticDefgradTransvIsotropElastViscoplast: Local Newton Loop did not "
               "converge for the given loop settings!";
      case Mat::ViscoplastErrorType::SingularJacobian:
        return "Error in InelasticDefgradTransvIsotropElastViscoplast: singular Jacobian after "
               "converged Local Newton Loop, which does not allow for the analytical evaluation "
               "of "
               "the linearization!";
      case Mat::ViscoplastErrorType::FailedSolAnalytLinearization:
        return "Error in InelasticDefgradTransvIsotropElastViscoplast: solution of the linear "
               "system "
               "in the analytical linearization failed";
        break;
      case Mat::ViscoplastErrorType::FailedComputationFlowResistance:
        return "Error in InelasticDefgradTransvIsotropElastViscoplast: Failed while computing the "
               "flow resistance of the viscoplasticity law";
        break;
      case Mat::ViscoplastErrorType::FailedComputationFlowResistanceDerivs:
        return "Error in InelasticDefgradTransvIsotropElastViscoplast: Failed while computing the "
               "derivatives of the flow resistance of the viscoplasticity law";
        break;
      case Mat::ViscoplastErrorType::FailedLogEval:
        return "Error in InelasticDefgradTransvIsotropElastViscoplast: Failed in evaluating the "
               "matrix logarithm or its derivative with respect to the argument";
        break;
      case Mat::ViscoplastErrorType::FailedExpEval:
        return "Error in InelasticDefgradTransvIsotropElastViscoplast: Failed in evaluating the "
               "matrix exponential or its derivative with respect to the argument";
        break;
      default:
        FOUR_C_THROW("to_string(Mat::ViscoplastErrorType): You should not be here!");
    }
  }


  /// enum class for time integration types (integration of internal
  /// variables in the Local Newton Loop of InelasticDefgradTransvIsotropElastViscoplast)
  enum class ViscoplastTimIntType
  {
    Standard,     // standard time integration,
    Logarithmic,  // time integration with logarithmically transformed residual equation for the
                  // evolution of the plastic deformation gradient
  };

  /// enum class for material linearization types (linearization of
  /// InelasticDefgradTransvIsotropElastViscoplast)
  enum class ViscoplastLinearizationType
  {
    Analytic,  // analytical linearization involving the solution of a linear system of equations,
    PerturbBased,  // linearization based on perturbing the current state
  };


  /// get the time integration type (Local Newton Loop of
  /// InelasticDefgradTransvIsotropElastViscoplast) from the
  /// user-specified string in the input file
  ViscoplastTimIntType get_time_integration_type(const std::string& timint_string);

  /// get the material linearization type (InelasticDefgradTransvIsotropElastViscoplast) from the
  /// user-specified string in the input file
  ViscoplastLinearizationType get_linearization_type(const std::string& linearization_string);

  // names of the various error types
  inline std::map<ViscoplastErrorType, std::string> ViscoplastErrorNames = {
      {ViscoplastErrorType::NegativePlasticStrain, "NegativePlasticStrain"},
      {ViscoplastErrorType::OverflowError, "OverflowError"},
      {ViscoplastErrorType::NoPlasticIncompressibility, "NoPlasticIncompressibility"},
      {ViscoplastErrorType::FailedSolLinSystLNL, "FailedSolLinSystLNL"},
      {ViscoplastErrorType::FailedDetermLineSearchParam, "FailedDetermLineSearchParam"},
      {ViscoplastErrorType::NoConvergenceLNL, "NoConvergenceLNL"},
      {ViscoplastErrorType::SingularJacobian, "SingularJacobian"},
      {ViscoplastErrorType::FailedSolAnalytLinearization, "FailedSolAnalytLinearization"},
      {ViscoplastErrorType::FailedLogEval, "FailedLogEval"},
      {ViscoplastErrorType::FailedExpEval, "FailedExpEval"},
  };

  //! class containing utilities for analyzing the material time integration:
  //! error types, number of line searches, ... (InelastDefgradTransvIsotropElastViscoplast)
  class TimIntAnalysisUtils
  {
   public:
    //! number of substeps for the current evaluation (LNL)
    unsigned int eval_num_of_substeps_ = 0;

    //! total number of substeps
    unsigned int total_num_of_substeps_ = 0;

    //! number of iterations for the current evaluation (LNL)
    unsigned int eval_num_of_iters_ = 0;

    //! total number of LNL iterations
    unsigned int total_num_of_iters_ = 0;

    //! number of repredictorizations for the current evaluation (LNL)
    unsigned int eval_num_of_repredict_ = 0;

    //! total number of LNL repredictorizations
    unsigned int total_num_of_repredict_ = 0;

    //! number of line searches for the current evaluation (LNL)
    unsigned int eval_num_of_line_search_ = 0;

    //! line search: the number of times the step size \f$ \alpha \f$ of
    //! the last iteration (Local Newton Loop) deviates from 1.0
    //! (currently evaluated time step)
    unsigned int eval_num_of_alpha_neq_1_last_iter = 0;

    //! line search: the number of times the step size \f$ \alpha \f$
    //! deviates from 1.0 in all iterations of the Local Newton Loop
    //! (currently evaluated time step)
    unsigned int eval_num_of_alpha_neq_1 = 0;

    //! total number of LNL line searches
    unsigned int total_num_of_line_search_ = 0;

    //! line search: total number of times the step size \f$ \alpha \f$ of
    //! the last iteration (Local Newton Loop) deviates from 1.0
    unsigned int total_num_of_alpha_neq_1_last_iter = 0;

    //! line search: the number of times the step size \f$ \alpha \f$
    //! deviates from 1.0 in all iterations of the Local Newton Loop
    //! (currently evaluated time step)
    unsigned int total_num_of_alpha_neq_1 = 0;

    //! number of times the LNL convergences directly in its first
    //! iteration (due to a good predictor!) for the current evaluation
    unsigned int eval_num_of_first_iter_convergences = 0;

    //! total number of times the LNL convergences directly in its first
    //! iteration (due to a good predictor!)
    unsigned int total_num_of_first_iter_convergences = 0;

    //! timer for the current evaluation, from the start of preevaluate to the end of update
    Teuchos::Time eval_teuchos_timer_{
        "InelasticDefgradTransvIsotropElastViscoplast::from_preevaluate_to_update"};

    //! evaluation time
    double eval_time_;

    //! total time
    double total_time_;

    //! error map of the current evaluation, from the first preevaluate of this time step to the
    //! first preevaluate of the next
    std::map<Mat::ViscoplastErrorType, unsigned int> eval_error_map_ = {
        {Mat::ViscoplastErrorType::NegativePlasticStrain, 0},
        {Mat::ViscoplastErrorType::OverflowError, 0},
        {Mat::ViscoplastErrorType::NoPlasticIncompressibility, 0},
        {Mat::ViscoplastErrorType::FailedSolLinSystLNL, 0},
        {Mat::ViscoplastErrorType::FailedDetermLineSearchParam, 0},
        {Mat::ViscoplastErrorType::NoConvergenceLNL, 0},
        {Mat::ViscoplastErrorType::SingularJacobian, 0},
        {Mat::ViscoplastErrorType::FailedSolAnalytLinearization, 0},
        {Mat::ViscoplastErrorType::FailedLogEval, 0},
        {Mat::ViscoplastErrorType::FailedExpEval, 0},
    };

    //! error map of the total evaluation
    std::map<Mat::ViscoplastErrorType, unsigned int> total_error_map_ = {
        {Mat::ViscoplastErrorType::NegativePlasticStrain, 0},
        {Mat::ViscoplastErrorType::OverflowError, 0},
        {Mat::ViscoplastErrorType::NoPlasticIncompressibility, 0},
        {Mat::ViscoplastErrorType::FailedSolLinSystLNL, 0},
        {Mat::ViscoplastErrorType::FailedDetermLineSearchParam, 0},
        {Mat::ViscoplastErrorType::NoConvergenceLNL, 0},
        {Mat::ViscoplastErrorType::SingularJacobian, 0},
        {Mat::ViscoplastErrorType::FailedSolAnalytLinearization, 0},
        {Mat::ViscoplastErrorType::FailedLogEval, 0},
        {Mat::ViscoplastErrorType::FailedExpEval, 0},
    };

    //! runtime csv writer
    std::optional<Core::IO::RuntimeCsvWriter> csv_writer_;

    //! simulation time instant and time step
    double sim_time_ = 0.0;
    int sim_timestep_ = 0.0;

    //! was the pre_evaluate method of the first element called?
    bool pre_eval_called_ = false;

    //! how often was the update method called? (max. num_of_global_elements if one processor is
    //! considered)
    int num_update_calls_ = 0;

    //! reset method
    void reset()
    {
      eval_num_of_substeps_ = 0;
      eval_num_of_iters_ = 0;
      eval_num_of_repredict_ = 0;
      eval_num_of_line_search_ = 0;
      eval_teuchos_timer_.reset();
      eval_time_ = 0;
      eval_error_map_ = {
          {Mat::ViscoplastErrorType::NegativePlasticStrain, 0},
          {Mat::ViscoplastErrorType::OverflowError, 0},
          {Mat::ViscoplastErrorType::NoPlasticIncompressibility, 0},
          {Mat::ViscoplastErrorType::FailedSolLinSystLNL, 0},
          {Mat::ViscoplastErrorType::FailedDetermLineSearchParam, 0},
          {Mat::ViscoplastErrorType::NoConvergenceLNL, 0},
          {Mat::ViscoplastErrorType::SingularJacobian, 0},
          {Mat::ViscoplastErrorType::FailedSolAnalytLinearization, 0},
          {Mat::ViscoplastErrorType::FailedLogEval, 0},
          {Mat::ViscoplastErrorType::FailedExpEval, 0},
      };
      eval_num_of_alpha_neq_1 = 0;
      eval_num_of_alpha_neq_1_last_iter = 0;
      eval_num_of_first_iter_convergences = 0;
    }

    //! initialize csv_writer
    void init_csv_writer()
    {
      // get structure discretization
      std::shared_ptr<Core::FE::Discretization> structure_dis =
          Global::Problem::instance()->get_dis("structure");

      // check whether we are using a single processor! (no implementation for multiple processors
      // yet, and also not really required)
      int my_rank = Core::Communication::my_mpi_rank(structure_dis->get_comm());
      FOUR_C_ASSERT_ALWAYS(my_rank == 0,
          "InelasticDefgradTransvIsotropElastViscoplast: No implementation of time integration "
          "output "
          "for multiple processors");

      // create csv_writer and register its columns
      csv_writer_.emplace(
          my_rank, *Global::Problem::instance()->output_control_file(), "timint_output");
      csv_writer_->register_data_vector("Eval. substeps", 1, 16);
      csv_writer_->register_data_vector("Eval. iterations", 1, 16);
      csv_writer_->register_data_vector("Eval. repredictorizations", 1, 16);
      csv_writer_->register_data_vector("Eval. line searches", 1, 16);
      csv_writer_->register_data_vector("Eval. # of times: alpha neq 1 (all LNL iters)", 1, 16);
      csv_writer_->register_data_vector("Eval. # of times: alpha neq 1 (last LNL iter)", 1, 16);
      csv_writer_->register_data_vector("Eval. # of first LNL iter. convergences", 1, 16);
      csv_writer_->register_data_vector("Eval. time", 1, 16);
      csv_writer_->register_data_vector("Total substeps", 1, 16);
      csv_writer_->register_data_vector("Total iterations", 1, 16);
      csv_writer_->register_data_vector("Total repredictorizations", 1, 16);
      csv_writer_->register_data_vector("Total line searches", 1, 16);
      csv_writer_->register_data_vector("Total # of times: alpha neq 1 (all LNL iters)", 1, 16);
      csv_writer_->register_data_vector("Total # of times: alpha neq 1 (last LNL iter)", 1, 16);
      csv_writer_->register_data_vector("Total # of first LNL iter. convergences", 1, 16);
      csv_writer_->register_data_vector("Total time", 1, 16);
      for (const auto& [key, value] : Mat::ViscoplastErrorNames)
      {
        csv_writer_->register_data_vector(
            "Eval. Error " + std::to_string(static_cast<int>(key)) + ": " + value, 1, 16);
        csv_writer_->register_data_vector(
            "Total Error " + std::to_string(static_cast<int>(key)) + ": " + value, 1, 16);
      }
    }

    //! update total values
    void update_total()
    {
      total_num_of_substeps_ += eval_num_of_substeps_;
      total_num_of_iters_ += eval_num_of_iters_;
      total_num_of_repredict_ += eval_num_of_repredict_;
      total_num_of_line_search_ += eval_num_of_line_search_;
      total_num_of_alpha_neq_1 += eval_num_of_alpha_neq_1;
      total_num_of_alpha_neq_1_last_iter += eval_num_of_alpha_neq_1_last_iter;
      total_time_ += eval_time_;
      for (const auto& [error_type, error_count] : eval_error_map_)
      {
        total_error_map_[error_type] += error_count;
      }
    }

    //! write to csv after each timestep
    void write_to_csv()
    {
      // output data
      std::map<std::string, std::vector<double>> output_data;
      output_data["Eval. substeps"] = {static_cast<double>(eval_num_of_substeps_)};
      output_data["Total substeps"] = {static_cast<double>(total_num_of_substeps_)};
      output_data["Eval. iterations"] = {static_cast<double>(eval_num_of_iters_)};
      output_data["Total iterations"] = {static_cast<double>(total_num_of_iters_)};
      output_data["Eval. repredictorizations"] = {static_cast<double>(eval_num_of_repredict_)};
      output_data["Total repredictorizations"] = {static_cast<double>(total_num_of_repredict_)};
      output_data["Eval. line searches"] = {static_cast<double>(eval_num_of_line_search_)};
      output_data["Total line searches"] = {static_cast<double>(total_num_of_line_search_)};
      output_data["Eval. time"] = {static_cast<double>(eval_time_)};
      output_data["Total time"] = {static_cast<double>(total_time_)};
      output_data["Eval. # of times: alpha neq 1 (all LNL iters)"] = {
          static_cast<double>(eval_num_of_alpha_neq_1)};
      output_data["Eval. # of times: alpha neq 1 (last LNL iter)"] = {
          static_cast<double>(eval_num_of_alpha_neq_1_last_iter)};
      output_data["Eval. # of first LNL iter. convergences"] = {
          static_cast<double>(eval_num_of_first_iter_convergences)};
      output_data["Total # of times: alpha neq 1 (all LNL iters)"] = {
          static_cast<double>(total_num_of_alpha_neq_1)};
      output_data["Total # of times: alpha neq 1 (last LNL iter)"] = {
          static_cast<double>(total_num_of_alpha_neq_1_last_iter)};
      output_data["Total # of first LNL iter. convergences"] = {
          static_cast<double>(total_num_of_first_iter_convergences)};


      for (const auto& [key, value] : Mat::ViscoplastErrorNames)
      {
        output_data["Eval. Error " + std::to_string(static_cast<int>(key)) + ": " + value] = {
            static_cast<double>(eval_error_map_[key])};
        output_data["Total Error " + std::to_string(static_cast<int>(key)) + ": " + value] = {
            static_cast<double>(total_error_map_[key])};
      }

      // write output data to csv
      csv_writer_->write_data_to_file(sim_time_, sim_timestep_, output_data);
    }
  };

  /// enum class for state quantity evaluations in
  /// InelasticDefgradTransvIsotropElastViscoplast: what is the aim of
  /// the evaluation? (full evaluation, or only partial, e.g. only the
  /// plastic strain rate,...)
  enum class ViscoplastStateQuantityEvalType
  {
    FullEval,               // full evaluation (full call of the evaluate_state_quantities method)
    PlasticStrainRateOnly,  // return in evaluate_state_quantities once the plastic strain rate has
                            // been evaluated
  };

  /// enum class for evaluations of the state quantity derivatives in
  /// InelasticDefgradTransvIsotropElastViscoplast: what is the aim of
  /// the evaluation? (full evaluation, or only partial, e.g. only the
  /// derivatives of the plastic strain rate,...)
  enum class ViscoplastStateQuantityDerivEvalType
  {
    FullEval,  // full evaluation (full call of the evaluate_state_quantity_derivatives method)
    PlasticStrainRateDerivsOnly,  // return in evaluate_state_quantity_derivatives once the
                                  // derivatives of the plastic strain rate have been evaluated
  };

  /// defines
  // flag for debug output (viscoplastic material) related to
  // time integration
  // #define DEBUGVPLAST_TIMINT

  // flag for debug output (viscoplastic material) related to the
  // inverse inelastic defgrad computation; less detailed than
  // DEBUGVPLAST_TIMINT
  // #define DEBUGVPLAST_INELDEFGRAD

  // flag for debug output (viscoplastic material) related to the
  // stiffness contribution
  // #define DEBUGVPLAST_LINEARIZATION

}  // namespace Mat


FOUR_C_NAMESPACE_CLOSE

#endif