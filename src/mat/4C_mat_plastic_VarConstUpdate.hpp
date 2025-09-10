// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_PLASTIC_VARCONSTUPDATE_HPP
#define FOUR_C_MAT_PLASTIC_VARCONSTUPDATE_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_plasticelasthyper.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration due to avoid header definition
namespace Mat
{
  namespace Elastic
  {
    class Summand;
  }

  // forward declaration
  class PlasticElastHyperVCU;

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// Collection of hyperelastic materials
    ///
    /// Storage map of hyperelastic summands.
    class PlasticElastHyperVCU : public Mat::PAR::PlasticElastHyper
    {
      friend class Mat::PlasticElastHyperVCU;

     public:
      /// standard constructor
      ///
      /// This constructor recursively calls the constructors of the
      /// parameter sets of the hyperelastic summands.
      PlasticElastHyperVCU(const Core::Mat::PAR::Parameter::Data& matdata);

      /// @name material parameters
      //@{

      /// provide access to material/summand by its ID
      std::shared_ptr<const Mat::Elastic::Summand> material_by_id(
          const int id  ///< ID to look for in collection of summands
      ) const;

      /// create material instance of matching type with my parameters
      std::shared_ptr<Core::Mat::Material> create_material() override;
      //@}

    };  // class PlasticElastHyper

  }  // namespace PAR

  class PlasticElastHyperVCUType : public Core::Communication::ParObjectType
  {
   public:
    std::string name() const override { return "PlasticElastHyperVCUType"; }

    static PlasticElastHyperVCUType& instance() { return instance_; };

    Core::Communication::ParObject* create(Core::Communication::UnpackBuffer& buffer) override;

   private:
    static PlasticElastHyperVCUType instance_;
  };


  /*----------------------------------------------------------------------*/
  /// Collection of hyperelastic materials
  class PlasticElastHyperVCU : public Mat::PlasticElastHyper
  {
   public:
    /// construct empty material object
    PlasticElastHyperVCU();

    /// construct the material object given material parameters
    explicit PlasticElastHyperVCU(Mat::PAR::PlasticElastHyperVCU* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of parobject.H (this file) and should return it in this method.
    int unique_par_object_id() const override
    {
      return PlasticElastHyperVCUType::instance().unique_par_object_id();
    }

    /// \brief Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by unique_par_object_id() which will then
    /// identify the exact class on the receiving processor.
    ///
    /// \param data (in/out): char vector to store class information
    void pack(Core::Communication::PackBuffer& data) const override;

    /// \brief Unpack data from a char vector into this class
    ///
    /// The vector data contains all information to rebuild the
    /// exact copy of an instance of a class on a different processor.
    /// The first entry in data has to be an integer which is the unique
    /// parobject id defined at the top of this file and delivered by
    /// unique_par_object_id().
    ///
    /// \param data (in) : vector storing all data to be unpacked into this
    ///                    instance.
    void unpack(Core::Communication::UnpackBuffer& buffer) override;

    //@}

    /// material type
    Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::m_plelasthyperVCU;
    }

    /// return copy of this material object
    std::shared_ptr<Core::Mat::Material> clone() const override
    {
      return std::make_shared<PlasticElastHyperVCU>(*this);
    }


    /// hyperelastic stress response plus elasticity tensor
    /// (pure virtual in material base class. Not allowed here)
    void evaluate(const Core::LinAlg::Tensor<double, 3, 3>* defgrad,
        const Core::LinAlg::SymmetricTensor<double, 3, 3>& glstrain,
        const Teuchos::ParameterList& params, Core::LinAlg::SymmetricTensor<double, 3, 3>& stress,
        Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat, int gp,
        int eleGID) override;  ///< Element GID

    /// setup material data
    void setup(int numgp, const Discret::Elements::Fibers& fibers,
        const std::optional<Discret::Elements::CoordinateSystem>& coord_system) override;

    /// update sumands
    void update() override;

    //! return names of visualization data
    void vis_names(std::map<std::string, int>& names) const override;

    //! return visualization data
    bool vis_data(
        const std::string& name, std::vector<double>& data, int numgp, int eleID) const override;

    // ********************************************************* //
    // ** here are some functions we don't need variationally ** //
    // ** consistent constitutive updates. They all give an   ** //
    // ** error when called. They are only here because of the * //
    // ** inheritance                                         ** //
    // ********************************************************* //

    /// evaluate quantities for elastic stiffness matrix
    /// in consideration of plastic history/deformation
    virtual void evaluate_elast(const Core::LinAlg::Matrix<3, 3>* defgrd,
        const Core::LinAlg::Matrix<3, 3>* deltaLp, const Teuchos::ParameterList& params,
        Core::LinAlg::Matrix<6, 1>* pk2, Core::LinAlg::Matrix<6, 6>* cmat, const int gp,
        const int eleGID)
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
    }

    /// evaluate stresses and stiffness contribution
    /// due to thermal expansion
    virtual void evaluate_thermal_stress(const Core::LinAlg::Matrix<3, 3>* defgrd,
        const double temp, const Teuchos::ParameterList& params, Core::LinAlg::Matrix<6, 1>* pk2,
        Core::LinAlg::Matrix<6, 6>* cmat, const int gp, const int eleGID)
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
    }

    /// evaluate stresses and stiffness contribution
    /// due to thermal expansion
    virtual void evaluate_c_tvol(const Core::LinAlg::Matrix<3, 3>* defgrd,
        const Teuchos::ParameterList& params, Core::LinAlg::Matrix<6, 1>* cTvol,
        Core::LinAlg::Matrix<6, 6>* dCTvoldE, const int gp, const int eleGID)
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
    }

    /// evaluate everything needed for the condensation of the plastic deformation
    /// at element level. (with zero plastic spin)
    virtual void evaluate_plast(const Core::LinAlg::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const Core::LinAlg::Matrix<3, 3>* deltaDp,  ///< symmetric part of plastic flow increment
        const double temp,                          ///< current temperature
        const Teuchos::ParameterList& params,       ///< Container for additional information
        Core::LinAlg::Matrix<6, 6>* dPK2dDp,        ///< derivative of PK2 w.r.t. F_p^{-1}
        Core::LinAlg::Matrix<6, 1>* NCP,            ///< NCP function
        Core::LinAlg::Matrix<6, 6>* dNCPdC,         ///< derivative of NCP function w.r.t. RCG
        Core::LinAlg::Matrix<6, 6>* dNCPdDp,        ///< derivative of NCP function w.r.t. deltaLp
        bool* active,                               ///< gauss point is active
        bool* elast,         ///< gauss point needs condensation if it is not elast
        bool* as_converged,  ///< convergence of active set (false, if as has changed)
        const int gp,        ///< gauss point
        Core::LinAlg::Matrix<6, 1>*
            dNCPdT,  ///< derivative of NCP function w.r.t. temperature (only in TSI case)
        Core::LinAlg::Matrix<6, 1>* dHdC,  ///< derivative of Heating w.r.t. RCG (only in TSI case)
        Core::LinAlg::Matrix<6, 1>*
            dHdDp,        ///< derivative of Heating w.r.t. deltaLp (only in TSI case)
        const double dt,  ///< time step size
        const int eleGID  ///< global ID of element
    )
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
    }

    /// evaluate everything needed for the condensation of the plastic deformation
    /// at element level. (with plastic spin)
    virtual void evaluate_plast(const Core::LinAlg::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const Core::LinAlg::Matrix<3, 3>*
            deltaLp,                           ///< plastic deformation gradient (non-symmetric)
        const double temp,                     ///< current temperature
        const Teuchos::ParameterList& params,  ///< Container for additional information
        Core::LinAlg::Matrix<6, 9>* dPK2dLp,   ///< derivative of PK2 w.r.t. F_p^{-1}
        Core::LinAlg::Matrix<9, 1>* NCP,       ///< NCP function
        Core::LinAlg::Matrix<9, 6>* dNCPdC,    ///< derivative of NCP function w.r.t. RCG
        Core::LinAlg::Matrix<9, 9>* dNCPdLp,   ///< derivative of NCP function w.r.t. deltaLp
        bool* active,                          ///< gauss point is active
        bool* elast,         ///< gauss point needs condensation if it is not elast
        bool* as_converged,  ///< convergence of active set (false, if as has changed)
        const int gp,        ///< gauss point
        Core::LinAlg::Matrix<9, 1>*
            dNCPdT,  ///< derivative of NCP function w.r.t. temperature (only in TSI case)
        Core::LinAlg::Matrix<6, 1>* dHdC,  ///< derivative of Heating w.r.t. RCG (only in TSI case)
        Core::LinAlg::Matrix<9, 1>*
            dHdLp,        ///< derivative of Heating w.r.t. deltaLp (only in TSI case)
        const double dt,  ///< time step size
        const int eleGID  ///< global ID of element
    )
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
    }

    /// update plastic history variables
    void update_gp(const int gp, const Core::LinAlg::Matrix<3, 3>* deltaDp) override
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
    }

    /// get plastic algorithm parameters
    void get_params(double s, double cpl) override
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
    }

    /// is this GP active
    virtual bool active(int gp)
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
      return false;
    }

    /// heating at this gp
    double& hep_diss(int gp) override
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
      static double a = 0.;
      return a;
    }

    /// derivative of heating at this gp
    Core::LinAlg::SerialDenseVector& d_hep_diss_dd(int gp) override
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
      static Core::LinAlg::SerialDenseVector tmp(0);
      return tmp;
    }

    // derivative of heating w.r.t. temperature
    double& d_hep_dt(int gp) override
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
      static double a = 0.;
      return a;
    }

    // derivative of heating at each gp w.r.t. nodal temperature vector
    // (only EAS contribution)
    std::shared_ptr<std::vector<Core::LinAlg::SerialDenseVector>> d_hep_d_teas() override
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
      return nullptr;
    }

   protected:
    virtual void evaluate_ncp(const Core::LinAlg::Matrix<3, 3>* mStr,
        const Core::LinAlg::Matrix<6, 6>* dMdC, const Core::LinAlg::Matrix<6, 9>* dMdFpinv,
        const Core::LinAlg::Matrix<6, 9>* dPK2dFpinv, const Core::LinAlg::Matrix<3, 3>* deltaDp,
        const int gp, const double temp, Core::LinAlg::Matrix<6, 1>* NCP,
        Core::LinAlg::Matrix<6, 6>* dNCPdC, Core::LinAlg::Matrix<6, 6>* dNCPdDp,
        Core::LinAlg::Matrix<6, 1>* dNCPdT, Core::LinAlg::Matrix<6, 6>* dPK2dDp, bool* active,
        bool* elast, bool* as_converged, Core::LinAlg::Matrix<6, 1>* dHdC,
        Core::LinAlg::Matrix<6, 1>* dHdDp, const Teuchos::ParameterList& params, const double dt)
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
    }


    virtual void evaluate_nc_pand_spin(const Core::LinAlg::Matrix<3, 3>* mStr,
        const Core::LinAlg::Matrix<6, 6>* dMdC, const Core::LinAlg::Matrix<6, 9>* dMdFpinv,
        const Core::LinAlg::Matrix<6, 9>* dPK2dFpinv, const Core::LinAlg::Matrix<3, 3>* deltaLp,
        const int gp, Core::LinAlg::Matrix<9, 1>* NCP, Core::LinAlg::Matrix<9, 6>* dNCPdC,
        Core::LinAlg::Matrix<9, 9>* dNCPdLp, Core::LinAlg::Matrix<6, 9>* dPK2dLp, bool* active,
        bool* elast, bool* as_converged, const double dt)
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
    }



    virtual void eval_dce_dlp(const Core::LinAlg::Matrix<3, 3> fpi,
        const Core::LinAlg::Matrix<3, 3>* defgrd, const Core::LinAlg::Matrix<6, 6> Dexp,
        const Core::LinAlg::Matrix<3, 3> cetrial, const Core::LinAlg::Matrix<3, 3> explp,
        Core::LinAlg::Matrix<6, 6>& dceDdeltalp, Core::LinAlg::Matrix<9, 6>& dFpiDdeltaDp);



    virtual void evaluate_rhs(const int gp, const Core::LinAlg::Matrix<3, 3> dLp,
        const Core::LinAlg::Matrix<3, 3> defgrd, Core::LinAlg::SymmetricTensor<double, 3, 3>& eeOut,
        Core::LinAlg::Matrix<5, 1>& rhs, Core::LinAlg::Matrix<5, 1>& rhsElast,
        Core::LinAlg::Matrix<6, 6>& dcedlp, Core::LinAlg::Matrix<9, 6>& dFpiDdeltaDp,
        const Teuchos::ParameterList& params, const int eleGID);


    virtual void yield_function(const double last_ai, const double norm_dLp,
        const Core::LinAlg::Matrix<3, 3> ExpEqui, const Core::LinAlg::Matrix<3, 3> cetr,
        const Core::LinAlg::Matrix<6, 1> str, double* yieldFunc,
        Core::LinAlg::Matrix<3, 3>& devMandelStr, Core::LinAlg::Matrix<3, 3>& MandelStr);

    virtual void comp_elast_quant(const Core::LinAlg::Matrix<3, 3>* defgrd,
        const Core::LinAlg::Matrix<3, 3> fpi, const Core::LinAlg::Matrix<3, 3> MatExp,
        Core::LinAlg::Matrix<3, 3>* cetrial, Core::LinAlg::SymmetricTensor<double, 3, 3>& Ee);

    virtual void evaluate_plast(Core::LinAlg::Matrix<6, 9>& dPK2dFpinvIsoprinc,
        const Core::LinAlg::Matrix<3, 1>& gamma, const Core::LinAlg::Matrix<8, 1>& delta,
        const Core::LinAlg::Matrix<3, 3>& id2, const Core::LinAlg::Matrix<6, 1>& Cpi,
        const Core::LinAlg::Matrix<3, 3>& Fpi, const Core::LinAlg::Matrix<3, 3>& CpiC,
        const Core::LinAlg::Matrix<9, 1>& CFpi, const Core::LinAlg::Matrix<9, 1>& CFpiCei,
        const Core::LinAlg::Matrix<6, 1>& ircg, const Core::LinAlg::Matrix<3, 3>& FpiCe,
        const Core::LinAlg::Matrix<9, 1>& CFpiCe, const Core::LinAlg::Matrix<6, 1>& CpiCCpi);

    virtual void evaluate_kin_quant_plast(int gp, int eleGID,
        const Core::LinAlg::Matrix<3, 3>* defgrd, const Core::LinAlg::Matrix<3, 3>* fpi,
        Core::LinAlg::Matrix<3, 1>& gamma, Core::LinAlg::Matrix<8, 1>& delta,
        Core::LinAlg::Matrix<3, 3>& id2, Core::LinAlg::Matrix<6, 1>& Cpi,
        Core::LinAlg::Matrix<3, 3>& CpiC, Core::LinAlg::Matrix<9, 1>& CFpi,
        Core::LinAlg::Matrix<9, 1>& CFpiCei, Core::LinAlg::Matrix<6, 1>& ircg,
        Core::LinAlg::Matrix<3, 3>& FpiCe, Core::LinAlg::Matrix<9, 1>& CFpiCe,
        Core::LinAlg::Matrix<6, 1>& CpiCCpi);

    virtual void dpk2d_fpi(int gp, int eleGID, const Core::LinAlg::Matrix<3, 3>* defgrd,
        const Core::LinAlg::Matrix<3, 3>* fpi, Core::LinAlg::Matrix<6, 9>& dPK2dFpinvIsoprinc);

    virtual void dpsiplast_dalphaiso(const double norm_dLp, const double last_alphaiso,
        const double isoHardMod, const double initYield, const double infYield,
        const double expIsoHard, double* dpsiplastdalphaiso);

    virtual void ce2nd_deriv(const Core::LinAlg::Matrix<3, 3>* defgrd,
        const Core::LinAlg::Matrix<3, 3> fpi, const Core::LinAlg::Matrix<3, 3> dLp,
        Core::LinAlg::Matrix<6, 6>* DDceDdLpDdLpVoigt);


    /// Access to material params
    Mat::PAR::PlasticElastHyperVCU* mat_params() const override { return params_; }

    /// get dissipation mode
    TSI::DissipationMode dis_mode() const override { return TSI::pl_flow; }

    /// inverse plastic deformation gradient for each Gauss point at current state
    std::vector<Core::LinAlg::Matrix<3, 3>> plastic_defgrd_inverse_;

    /// my material parameters
    Mat::PAR::PlasticElastHyperVCU* params_;
  };

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
