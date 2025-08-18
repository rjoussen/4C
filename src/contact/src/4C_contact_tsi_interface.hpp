// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CONTACT_TSI_INTERFACE_HPP
#define FOUR_C_CONTACT_TSI_INTERFACE_HPP

#include "4C_config.hpp"

#include "4C_contact_input.hpp"
#include "4C_contact_interface.hpp"

FOUR_C_NAMESPACE_OPEN


namespace CONTACT
{
  class TSIInterface : public Interface
  {
   public:
    /*!
    \brief Constructor

    */
    TSIInterface(const std::shared_ptr<Mortar::InterfaceDataContainer>& interfaceData_ptr,
        const int id, MPI_Comm comm, const int dim, const Teuchos::ParameterList& icontact,
        bool selfcontact);

    enum LinDmXMode
    {
      LinDM_Diss,
      LinDM_ThermoLM,
      linDM_ContactLMnormal
    };

    /*!
    \brief Assemble matrix LinStick containing linearizations

    This method builds an algebraic form of the FULL linearization
    of the tangential stick condition.

    */
    virtual void assemble_lin_stick(Core::LinAlg::SparseMatrix& linstickLMglobal,
        Core::LinAlg::SparseMatrix& linstickDISglobal,
        Core::LinAlg::SparseMatrix& linstickTEMPglobal,
        Core::LinAlg::Vector<double>& linstickRHSglobal);

    /*!
    \brief Assemble matrix LinSlip containing linearizations

    This method builds an algebraic form of the FULL linearization
    of the tangential slip condition. Concretely, this
    includes assembling the linearizations of the slave side
    nodal tangents and of the Mortar matrices D  and M.

    */
    virtual void assemble_lin_slip(Core::LinAlg::SparseMatrix& linslipLMglobal,
        Core::LinAlg::SparseMatrix& linslipDISglobal, Core::LinAlg::SparseMatrix& linslipTEMPglobal,
        Core::LinAlg::Vector<double>& linslipRHSglobal);


    /*!
    \brief Assemble contact heat conduction

    This method assembles the Thermo-Lagrange-Multiplier line
    in the global system of equations, containing linearizations
    wrt displacements, contact forces, temperatures
    */
    virtual void assemble_lin_conduct(Core::LinAlg::SparseMatrix& linConductDISglobal,
        Core::LinAlg::SparseMatrix& linConductTEMPglobal,
        Core::LinAlg::SparseMatrix& linConductThermoLMglobal,
        Core::LinAlg::SparseMatrix& linConductContactLMglobal);

    /*!
    \brief Assemble lumped mass matrix of the dual basis
           This is actually the D-Matrix. However, to keep
           the info where it comes from, we stick to that name

    */
    virtual void assemble_dual_mass_lumped(
        Core::LinAlg::SparseMatrix& dualMassGlobal,    /// dual mass matrix
        Core::LinAlg::SparseMatrix& linDualMassGlobal  // derivative of dual mass matrix wrt
                                                       // displacements multiplied with thermo-LM
    );

    /*!
      \brief Assemble the linearization of D and M times some nodal value

      This assembles  D_{jk,c) X_j and
                     -M_{jl,c} X_j (mind the minus sign!)
      where X_j is some (scalar) nodal value determined by "mode"
      */
    virtual void assemble_lin_dm_x(Core::LinAlg::SparseMatrix* linD_X,
        Core::LinAlg::SparseMatrix* linM_X, const double fac, const LinDmXMode mode,
        const std::shared_ptr<Core::LinAlg::Map> node_rowmap);

    /*!
      \brief Assemble D and M times the linearization of the nodal frictional dissipation

      This assembles  D_{jk) Diss_{j,c} and
                      M_{jl} Diss_{j,c}
                     where c is some discrete nodal DISPLACEMENT or LM dof.
      The dissipation is calculated by the tangential slip times the contact Lagrange multiplier
      Diss = \lambda^{contact}  \cdot  (1 - n \otimes n) \cdot jump
      */
    virtual void assemble_dm_lin_diss(Core::LinAlg::SparseMatrix* d_LinDissDISP,
        Core::LinAlg::SparseMatrix* m_LinDissDISP, Core::LinAlg::SparseMatrix* d_LinDissContactLM,
        Core::LinAlg::SparseMatrix* m_LinDissContactLM, const double fac);

    /*!
      \brief Assemble the linearization of D and M times the temperature

      This assembles the derivative wrt displacements
      linDM_Temp_{jl}= LMj_n* ( D_{jk,c) T_k
                               -M_{jm,c} T_m )
                       + LMj_d * n_{d,l} * ( D_{jk} T_k
                                            -M_{jm} T_m )
      and the derivative wrt the contact Lagrange multiplier
      lin_lm_{jk} = nj_{k} * (D_jk T_k - M_jl T_l)
                                 (mind the minus sign!)
      */
    virtual void assemble_lin_l_mn_dm_temp(
        const double fac, Core::LinAlg::SparseMatrix* lin_disp, Core::LinAlg::SparseMatrix* lin_lm);

    /*!
      \brief Assemble the D and M times the normal Lagrange multiplier

      This assembles lambda_n* ( D_{jk) LMn_j
                                -M_{jl} LMn_j ) (mind the minus sign! No sum over j!)
      */
    virtual void assemble_dm_l_mn(const double fac, Core::LinAlg::SparseMatrix* DM_LMn);

    /*!
      \brief Assemble inactive part of the thermal heat conduction equation

      This assembles lambda(thermo) = 0 for all inactive contact nodes
      */
    virtual void assemble_inactive(Core::LinAlg::SparseMatrix* linConductThermoLM);


    /*!
    \brief Initialize / reset interface for contact

    Derived version with some additional TSI related stuff!

    */
    void initialize() final;

   protected:
  };  // class


}  // namespace CONTACT

FOUR_C_NAMESPACE_CLOSE

#endif
