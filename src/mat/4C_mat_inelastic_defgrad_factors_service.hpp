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

  /// get the time integration type (Local Newton Loop of
  /// InelasticDefgradTransvIsotropElastViscoplast) from the
  /// user-specified string in the input file
  ViscoplastTimIntType get_time_integration_type(const std::string& timint_string);

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
  };

  //! class containing utilities for analyzing the material time integration:
  //! error types, number of line searches, ... (InelastDefgradTransvIsotropElastViscoplast)
  class TimIntAnalysisUtils
  {
   public:
    // number of substeps for the current evaluation (LNL)
    unsigned int eval_num_of_substeps_ = 0;

    // total number of substeps
    unsigned int total_num_of_substeps_ = 0;

    // number of iterations for the current evaluation (LNL)
    unsigned int eval_num_of_iters_ = 0;

    // total number of LNL iterations
    unsigned int total_num_of_iters_ = 0;

    // number of repredictorizations for the current evaluation (LNL)
    unsigned int eval_num_of_repredict_ = 0;

    // total number of LNL repredictorizations
    unsigned int total_num_of_repredict_ = 0;

    // number of line searches for the current evaluation (LNL)
    unsigned int eval_num_of_line_search_ = 0;

    // total number of LNL line searches
    unsigned int total_num_of_line_search_ = 0;

    // timer for the current evaluation, from the start of preevaluate to the end of update
    Teuchos::Time eval_teuchos_timer_{
        "InelasticDefgradTransvIsotropElastViscoplast::from_preevaluate_to_update"};

    // evaluation time
    double eval_time_;

    // total time
    double total_time_;

    // error map of the current evaluation, from the first preevaluate of this time step to the
    // first preevaluate of the next
    std::map<Mat::ViscoplastErrorType, unsigned int> eval_error_map_ = {
        {Mat::ViscoplastErrorType::NegativePlasticStrain, 0},
        {Mat::ViscoplastErrorType::OverflowError, 0},
        {Mat::ViscoplastErrorType::NoPlasticIncompressibility, 0},
        {Mat::ViscoplastErrorType::FailedSolLinSystLNL, 0},
        {Mat::ViscoplastErrorType::FailedDetermLineSearchParam, 0},
        {Mat::ViscoplastErrorType::NoConvergenceLNL, 0},
        {Mat::ViscoplastErrorType::SingularJacobian, 0},
        {Mat::ViscoplastErrorType::FailedSolAnalytLinearization, 0}};

    // error map of the total evaluation
    std::map<Mat::ViscoplastErrorType, unsigned int> total_error_map_ = {
        {Mat::ViscoplastErrorType::NegativePlasticStrain, 0},
        {Mat::ViscoplastErrorType::OverflowError, 0},
        {Mat::ViscoplastErrorType::NoPlasticIncompressibility, 0},
        {Mat::ViscoplastErrorType::FailedSolLinSystLNL, 0},
        {Mat::ViscoplastErrorType::FailedDetermLineSearchParam, 0},
        {Mat::ViscoplastErrorType::NoConvergenceLNL, 0},
        {Mat::ViscoplastErrorType::SingularJacobian, 0},
        {Mat::ViscoplastErrorType::FailedSolAnalytLinearization, 0}};


    // runtime csv writer
    std::optional<Core::IO::RuntimeCsvWriter> csv_writer_;

    // simulation time instant and time step
    double sim_time_ = 0.0;
    int sim_timestep_ = 0.0;

    // was the pre_evaluate method of the first element called?
    bool pre_eval_called_ = false;

    // how often was the update method called? (max. num_of_global_elements if one processor is
    // considered)
    int num_update_calls_ = 0;

    // reset method
    void reset()
    {
      eval_num_of_substeps_ = 0;
      eval_num_of_iters_ = 0;
      eval_num_of_repredict_ = 0;
      eval_num_of_line_search_ = 0;
      eval_teuchos_timer_.reset();
      eval_time_ = 0;
      eval_error_map_ = {{Mat::ViscoplastErrorType::NegativePlasticStrain, 0},
          {Mat::ViscoplastErrorType::OverflowError, 0},
          {Mat::ViscoplastErrorType::NoPlasticIncompressibility, 0},
          {Mat::ViscoplastErrorType::FailedSolLinSystLNL, 0},
          {Mat::ViscoplastErrorType::FailedDetermLineSearchParam, 0},
          {Mat::ViscoplastErrorType::NoConvergenceLNL, 0},
          {Mat::ViscoplastErrorType::SingularJacobian, 0},
          {Mat::ViscoplastErrorType::FailedSolAnalytLinearization, 0}};
    }

    // initialize csv_writer only if required
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
      csv_writer_->register_data_vector("Evaluated substeps", 1, 16);
      csv_writer_->register_data_vector("Evaluated iterations", 1, 16);
      csv_writer_->register_data_vector("Evaluated repredictorizations", 1, 16);
      csv_writer_->register_data_vector("Evaluated line searches", 1, 16);
      csv_writer_->register_data_vector("Evaluation time", 1, 16);
      csv_writer_->register_data_vector("Total substeps", 1, 16);
      csv_writer_->register_data_vector("Total iterations", 1, 16);
      csv_writer_->register_data_vector("Total repredictorizations", 1, 16);
      csv_writer_->register_data_vector("Total line searches", 1, 16);
      csv_writer_->register_data_vector("Total time", 1, 16);
      for (const auto& [key, value] : Mat::ViscoplastErrorNames)
      {
        csv_writer_->register_data_vector(
            "Evaluated Error " + std::to_string(static_cast<int>(key)) + ": " + value, 1, 16);
        csv_writer_->register_data_vector(
            "Total Error " + std::to_string(static_cast<int>(key)) + ": " + value, 1, 16);
      }
    }

    // update total values
    void update_total()
    {
      total_num_of_substeps_ += eval_num_of_substeps_;
      total_num_of_iters_ += eval_num_of_iters_;
      total_num_of_repredict_ += eval_num_of_repredict_;
      total_num_of_line_search_ += eval_num_of_line_search_;
      total_time_ += eval_time_;
      for (const auto& [error_type, error_count] : eval_error_map_)
      {
        total_error_map_[error_type] += error_count;
      }
    }

    // write to csv after each timestep
    void write_to_csv()
    {
      // output data
      std::map<std::string, std::vector<double>> output_data;
      output_data["Evaluated substeps"] = {static_cast<double>(eval_num_of_substeps_)};
      output_data["Total substeps"] = {static_cast<double>(total_num_of_substeps_)};
      output_data["Evaluated iterations"] = {static_cast<double>(eval_num_of_iters_)};
      output_data["Total iterations"] = {static_cast<double>(total_num_of_iters_)};
      output_data["Evaluated repredictorizations"] = {static_cast<double>(eval_num_of_repredict_)};
      output_data["Total repredictorizations"] = {static_cast<double>(total_num_of_repredict_)};
      output_data["Evaluated line searches"] = {static_cast<double>(eval_num_of_line_search_)};
      output_data["Total line searches"] = {static_cast<double>(total_num_of_line_search_)};
      output_data["Evaluation time"] = {static_cast<double>(eval_time_)};
      output_data["Total time"] = {static_cast<double>(total_time_)};

      for (const auto& [key, value] : Mat::ViscoplastErrorNames)
      {
        output_data["Evaluated Error " + std::to_string(static_cast<int>(key)) + ": " + value] = {
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



}  // namespace Mat


FOUR_C_NAMESPACE_CLOSE

#endif