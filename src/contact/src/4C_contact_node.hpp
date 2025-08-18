// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CONTACT_NODE_HPP
#define FOUR_C_CONTACT_NODE_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mortar_node.hpp"

#include <unordered_map>

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  // forward declaration
  class Node;

  class NodeType : public Core::Communication::ParObjectType
  {
   public:
    std::string name() const final { return "CONTACT::NodeType"; }

    static NodeType& instance() { return instance_; };

    Core::Communication::ParObject* create(Core::Communication::UnpackBuffer& buffer) override;

   private:
    static NodeType instance_;
  };

  /*!
   \brief A class containing additional data from contact nodes

   This class contains additional information from contact nodes which are
   are not needed for contact search and therefore are only available on the
   node's processor (ColMap). The class NodeDataContainer must be declared
   before the Node itself.

   */
  class NodeDataContainer
  {
   public:
    //! @name Constructors and destructors and related methods

    /*!
     \brief Standard Constructor

     */
    NodeDataContainer();

    /*!
     \brief Destructor

     */
    virtual ~NodeDataContainer() = default;
    /*!
     \brief Pack this class so that it can be communicated

     This function packs the data container. This is only called
     when the class has been initialized and the pointer to this
     class exists.

     */
    virtual void pack(Core::Communication::PackBuffer& data) const;

    /*!
     \brief Unpack data from a vector into this class

     This function unpacks the data container. This is only called
     when the class has been initialized and the pointer to this
     class exists.

     */
    virtual void unpack(Core::Communication::UnpackBuffer& buffer);

    //@}

    //! @name Access methods

    /*!
     \brief Return current nodal tangent t_xi (only for slave side!) (length 3)
     */
    virtual double* txi() { return txi_; }

    /*!
     \brief Return current nodal tangent t_eta (only for slave side!) (length 3)
     */
    virtual double* teta() { return teta_; }

    /*!
     \brief Return the weighted gap (scalar) of this node
     */
    virtual double& getg() { return grow_; }

    /*!
     \brief Return the weighted gap (scalar) of this node for NTS
     */
    virtual double& getgnts() { return gnts_; }

    /*!
     \brief Return the weighted gap (scalar) of this node for LTS
     */
    virtual double& getglts() { return glts_; }

    /*!
    \brief Return the gap (vector) of this node for LTL
    */
    virtual double* getgltl() { return gltl_; }

    /*!
    \brief Return the gap (vector) of this node for LTL
    */
    virtual double* getjumpltl() { return jumpltl_; }

    /*!
    \brief Return contact status of last converged state n (active=true)

    */
    virtual inline bool& active_old() { return activeold_; }

    /*!
     \brief Return the 'DerivN' map (vector) of this node

     These maps contain the directional derivatives of the node's
     averaged normal vector with respect to the slave displacements.
     A vector is used because the normal itself is a vector (2 or 3 components).

     */
    virtual std::vector<Core::Gen::Pairedvector<int, double>>& get_deriv_n() { return derivn_; }
    const std::vector<Core::Gen::Pairedvector<int, double>>& get_deriv_n() const { return derivn_; }
    virtual std::vector<Core::Gen::Pairedvector<int, double>>& get_deriv_tangent()
    {
      return derivEdge_;
    }
    const std::vector<Core::Gen::Pairedvector<int, double>>& get_deriv_tangent() const
    {
      return derivEdge_;
    }

    /*!
     \brief Return the 'DerivTxi' map (vector) of this node

     These maps contain the directional derivatives of the node's
     unit tangent vector t_xi with respect to the slave displacements.
     A vector is used because the tangent itself is a vector (2 or 3 components).

     */
    virtual std::vector<Core::Gen::Pairedvector<int, double>>& get_deriv_txi() { return derivtxi_; }

    /*!
     \brief Return the 'DerivTeta' map (vector) of this node

     These maps contain the directional derivatives of the node's
     unit tangent vector t_eta with respect to the slave displacements.
     A vector is used because the tangent itself is a vector (2 or 3 components).

     */
    virtual std::vector<Core::Gen::Pairedvector<int, double>>& get_deriv_teta()
    {
      return derivteta_;
    }

    /*!
     \brief Return the 'DerivD' map of this node

     D-matrix entries with respect to the slave/master displacements.
     It is a "map of maps", with the outer map containing all master
     node indices k adjacent to this node i and the inner map containing
     all directional derivatives l existing for D_ik.

     Note: In earlier versions this was just a simple map accounting for the
     diagonality of D when using dual shape functions. Due to the extension to
     support arbitrary types of shape functions, this is not possible anymore.

     */
    virtual std::map<int, std::map<int, double>>& get_deriv_d() { return derivd_; }
    virtual std::map<int, std::map<int, double>>& get_deriv_dlts() { return derivdlts_; }
    virtual std::map<int, std::map<int, double>>& get_deriv_dltl() { return derivdltl_; }

    /*!
     \brief Return the 'DerivM' map of maps of this node

     This map contains the directional derivatives of the node's
     M-matrix entries with respect to the slave/master displacements.
     It is a "map of maps", with the outer map containing all master
     node indices k adjacent to this node i and the inner map containing
     all directional derivatives l existing for M_ik.

     */
    virtual std::map<int, std::map<int, double>>& get_deriv_m() { return derivm_; }
    virtual std::map<int, std::map<int, double>>& get_deriv_mnts() { return derivmnts_; }
    virtual std::map<int, std::map<int, double>>& get_deriv_mlts() { return derivmlts_; }
    virtual std::map<int, std::map<int, double>>& get_deriv_mltl() { return derivmltl_; }

    /*!
     \brief Return one specific 'DerivD' map of this node

     This method returns the map of directional derivatives of one
     specific D-matrix D_ik entry of this node i.

     */
    virtual std::map<int, double>& get_deriv_d(int& k)
    {
      using CI = std::map<int, std::map<int, double>>::const_iterator;
      CI p = derivd_.find(k);
      if (p == derivd_.end()) FOUR_C_THROW("GetDerivD: No map entry existing for given index");
      return derivd_[k];
    }

    /*!
     \brief Return one specific 'DerivM' map of this node

     This method returns the map of directional derivatives of one
     specific M-matrix M_ik entry of this node i.

     */
    virtual std::map<int, double>& get_deriv_m(int& k)
    {
      using CI = std::map<int, std::map<int, double>>::const_iterator;
      CI p = derivm_.find(k);
      if (p == derivm_.end()) FOUR_C_THROW("GetDerivM: No map entry existing for given index");
      return derivm_[k];
    }

    /*!
     \brief Return the 'DerivG' map of this node

     This map contains the directional derivatives of the node's
     weighted gap entry g~ with respect to the slave/master displacements.

     */
    virtual std::map<int, double>& get_deriv_g() { return derivg_; }
    virtual std::map<int, double>& get_deriv_gnts() { return derivgnts_; }
    virtual std::map<int, double>& get_deriv_glts() { return derivglts_; }
    virtual std::vector<std::map<int, double>>& get_deriv_gltl() { return derivgltl_; }
    virtual std::vector<std::map<int, double>>& get_deriv_jumpltl() { return derivjumpltl_; }

    /*!
     \brief Return the 'DerivGlm' map of this node

     */
    virtual std::map<int, double>& get_deriv_gw() { return derivgw_; }

    /*!
     \brief Return the 'DerivW' map of this node

     This map contains the directional derivatives of the node's
     weighted wear increment entry w~ with respect to the slave/master displacements.

     */
    virtual std::map<int, double>& get_deriv_w() { return derivw_; }

    /*!
     \brief Return the 'DerivW' map of this node

     This map contains the directional derivatives of the node's
     weighted wear increment entry w~ with respect to the lagr. mult.

     */
    virtual std::map<int, double>& get_deriv_wlm() { return derivw_lm_; }

    /*!
     \brief Return scaling factor for weighted gap

     Note: This is only calculated when performing a penalty strategy

     */
    virtual double& kappa() { return kappa_; }

    /*!
     \brief Return the 'DerivZ map of this node

     This map contains the directional derivatives of the node's
     lagrange multiplier entries with respect to the slave/master displacements.

     Note: This is only calculated when performing a penalty strategy

     */
    virtual std::vector<std::map<int, double>>& get_deriv_z() { return derivz_; }

    virtual Core::Gen::Pairedvector<int, double>& get_alpha() { return alpha_; };
    virtual double& get_alpha_n() { return nalpha_; };

    /*!
     \brief Return old normal vector (length 3)
     */
    virtual inline double* normal_old() { return n_old_; }

    //! @}

   protected:
    //! nodal tangent t_xi for contact methods at newton step n
    double txi_[3];

    //! nodal tangent t_eta for contact methods at newton step n
    double teta_[3];

    //! nodal entry of g vector
    double grow_;

    //! nodal entry of g vector for node-to-segment
    double gnts_;

    //! nodal entry of g vector for line-to-segment
    double glts_;

    //! nodal entry of g vector for line-to-segment
    double gltl_[3];

    //! nodal entry of jump vector for line-to-segment
    double jumpltl_[3];

    //! scaling factor for node-to-segment-mortar hybrid formulation
    double nalpha_;

    //! true if contact node was in contact (last converged state n)
    bool activeold_;

    //! directional derivative of nodal normal
    std::vector<Core::Gen::Pairedvector<int, double>> derivn_;

    //! directional derivative of nodal tangent
    std::vector<Core::Gen::Pairedvector<int, double>> derivEdge_;

    //! directional derivative of nodal tangent t_xi
    std::vector<Core::Gen::Pairedvector<int, double>> derivtxi_;

    //! directional derivative of nodal tangent t_eta
    std::vector<Core::Gen::Pairedvector<int, double>> derivteta_;

    //! directional derivative of nodal D-matrix value
    std::map<int, std::map<int, double>> derivd_;

    //! directional derivative of nodal D-matrix value line-to-segment
    std::map<int, std::map<int, double>> derivdlts_;

    //! directional derivative of nodal D-matrix value line-to-line
    std::map<int, std::map<int, double>> derivdltl_;

    //! directional derivative of nodal M-matrix values
    std::map<int, std::map<int, double>> derivm_;

    //! directional derivative of nodal M-matrix values node-to-segment
    std::map<int, std::map<int, double>> derivmnts_;

    //! directional derivative of nodal M-matrix values line-to-segment
    std::map<int, std::map<int, double>> derivmlts_;

    //! directional derivative of nodal M-matrix values line-to-line
    std::map<int, std::map<int, double>> derivmltl_;

    //! directional derivative of nodal weighted gap value
    std::map<int, double> derivg_;

    //! directional derivative of node-to-segment gap value
    std::map<int, double> derivgnts_;

    //! directional derivative of line-to-segment gap value
    std::map<int, double> derivglts_;

    //! directional derivative of line-to-line gap value
    std::vector<std::map<int, double>> derivgltl_;

    //! directional derivative of line-to-line gap value
    std::vector<std::map<int, double>> derivjumpltl_;

    //! directional derivative of nodal weighted gap value
    std::map<int, double> derivgw_;

    //! directional derivative of nodal weighted wear value
    std::map<int, double> derivw_;

    //! directional derivative of nodal weighted wear value
    std::map<int, double> derivw_lm_;

    //! lin. of scaling factor for hybrid formulation
    Core::Gen::Pairedvector<int, double> alpha_;

    //! @name Penalty-related quantities
    //! @{

    //! gap-scaling factor kappa
    double kappa_;

    //! @}

    //! direction derivative of nodal z-matrix value
    std::vector<std::map<int, double>> derivz_;

    double n_old_[3];  // normal of previous time step

   private:
    //! Operator= (disabled)
    NodeDataContainer operator=(const NodeDataContainer& old) = delete;

    //! Copy constructor (disabled)
    NodeDataContainer(const NodeDataContainer& old) = delete;
  };


  /*!
  \brief A class containing additional data from poro contact nodes

  This class contains additional information from poro contact nodes which are
  are not needed for normal contact nodes.
  NodePoroDataContainer must be declared before the Node itself.

  */
  class NodePoroDataContainer
  {
   public:
    //! @name Constructors and destructors and related methods

    /*!
    \brief Standard Constructor

    */
    NodePoroDataContainer();

    /*!
    \brief Destructor

    */
    virtual ~NodePoroDataContainer() = default;
    /*!
    \brief Pack this class so that it can be communicated

    This function packs the datacontainer. This is only called
    when the class has been initialized and the pointer to this
    class exists.

    */
    virtual void pack(Core::Communication::PackBuffer& data) const;

    /*!
    \brief Unpack data from a vector into this class

    This function unpacks the datacontainer. This is only called
    when the class has been initialized and the pointer to this
    class exists.

    */
    virtual void unpack(Core::Communication::UnpackBuffer& buffer);

    //! @name Access methods

    /*!
    \brief Return the normal coupling condition (scalar) of this node
    */
    virtual double& getn_coup() { return ncouprow_; }

    /*!
    \brief Return the 'DerivnCoup' map of this node

    This map contains the directional derivatives of the node's
    normal coupling condition with respect to the slave/master displacements.

    */
    virtual std::map<int, double>& get_derivn_coup() { return derivncoup_; }

    /*!
    \brief Return the 'VelDerivnCoup' map of this node

    This map contains the derivatives of the node's
    normal coupling condition with respect to the velocities. (for one sided contact just slave!)

    */
    virtual std::map<int, double>& get_vel_derivn_coup() { return velderivncoup_; }

    /*!
    \brief Return the 'PresDerivnCoup' map of this node

    This map contains the derivatives of the node's
    normal coupling condition with respect to the pressures. (for one sided contact just slave!)
    //h.Willmann
    */
    virtual std::map<int, double>& get_pres_derivn_coup() { return presderivncoup_; }

    /*!
    \brief Return current nodal fluid pressure (porous media!) (length 3)
    */
    virtual double* fpres() { return &fpres_; }

    /*!
    \brief Return current nodal fluid velocity (porous media!) (length 3)
    */
    virtual double* fvel() { return fvel_; }

    /*!
    \brief Return current nodal structural velocity (porous media!) (length 3)
    */
    virtual double* svel() { return svel_; }

    /*!
    \brief Return current nodal lagrangean multiplier (length 3)
    */
    virtual double* poro_lm() { return porolm_; }

   protected:
    // don't want = operator and cctor
    NodePoroDataContainer operator=(const NodePoroDataContainer& old) = delete;
    NodePoroDataContainer(const NodePoroDataContainer& old) = delete;


    double ncouprow_;  // nodal entry of n-coupling vector
    std::map<int, double>
        derivncoup_;  // directional derivative of nodal weighted n-coupling vector
    std::map<int, double>
        velderivncoup_;  // velocity derivative of nodal weighted n-coupling vector
    std::map<int, double>
        presderivncoup_;  // pressure derivative of nodal weighted n-coupling vector //h.Willmann

    double fvel_[3];    // fluid velocity for porous problem
    double fpres_;      // fluid pressure for porous problem
    double svel_[3];    // structural velocity for porous problem
    double porolm_[3];  // lagrange multiplier from poro no penetration condition!
  };  // class NodePoroDataContainer

  /*!
  \brief A class containing additional data from tsi contact nodes

  This class contains additional information from tsi contact nodes which are
  are not needed for normal contact nodes.

  */
  class NodeTSIDataContainer
  {
   public:
    //! @name Constructors and destructors and related methods

    /*!
    \brief Standard Constructor

    */
    NodeTSIDataContainer(double t_ref, double t_dam);

    /*!
    \brief empty Constructor: unpack data later

    */
    NodeTSIDataContainer() {};

    /*!
    \brief Destructor

    */
    virtual ~NodeTSIDataContainer() = default;
    /*!
    \brief Pack this class so that it can be communicated

    This function packs the datacontainer. This is only called
    when the class has been initialized and the pointer to this
    class exists.

    */
    virtual void pack(Core::Communication::PackBuffer& data) const;

    /*!
    \brief Unpack data from a vector into this class

    This function unpacks the datacontainer. This is only called
    when the class has been initialized and the pointer to this
    class exists.

    */
    virtual void unpack(Core::Communication::UnpackBuffer& buffer);


    //! @name Access methods

    /*!
    \brief Return max (Temp_slave , Temp_master)
    */
    double& temp_master() { return temp_master_; }
    double temp_master() const { return temp_master_; }

    /*!
    \brief Return temperature
    */
    double& temp() { return temp_; }
    double temp() const { return temp_; }

    /*!
    \brief Return reference temperature
    */
    double& temp_ref() { return t_ref_; }
    double temp_ref() const { return t_ref_; }

    /*!
    \brief Return temperature
    */
    double& temp_dam() { return t_dam_; }
    double temp_dam() const { return t_dam_; }

    /*!
    \brief Return thermo Lagrange multiplier
    */
    double& thermo_lm() { return thermo_lm_; }

    std::map<int, double>& deriv_temp_master_disp() { return derivTempMasterDisp_; }
    const std::map<int, double>& deriv_temp_master_disp() const { return derivTempMasterDisp_; }
    std::map<int, double>& deriv_temp_master_temp() { return derivTempMasterTemp_; }
    const std::map<int, double>& deriv_temp_master_temp() const { return derivTempMasterTemp_; }

    void clear();


   protected:
    // don't want = operator and cctor
    NodeTSIDataContainer operator=(const NodeTSIDataContainer& old) = delete;
    NodeTSIDataContainer(const NodeTSIDataContainer& old) = delete;

    double temp_;
    double t_ref_;
    double t_dam_;
    double thermo_lm_;

    double temp_master_;
    std::map<int, double> derivTempMasterDisp_;
    std::map<int, double> derivTempMasterTemp_;

  };  // class NodeTSIDataContainer

  /*!
  \brief A class containing additional data from EHL contact nodes

  This class contains additional information from EHL contact nodes which are
  are not needed for normal contact nodes.

  */
  class NodeEhlDataContainer
  {
   public:
    //! empty c-tor
    NodeEhlDataContainer() {}

    //! d-tor
    virtual ~NodeEhlDataContainer() = default;

    //! pack for parallel communication
    virtual void pack(Core::Communication::PackBuffer& data) const {
      /* no need to pack, since terms are re-evaluated after parallel communication */
    };

    //! unpack and re-init after parallel communication
    virtual void unpack(Core::Communication::UnpackBuffer& buffer) {
      /* no need to pack, since terms are re-evaluated after parallel communication */
    };

    //! clear all stored data
    void clear()
    {
      weighted_relTangVel_.clear();
      deriv_weighted_relTangVel_.clear();
      weighted_avTangVel_.clear();
      deriv_weighted_avTangVel_.clear();
      tang_grad_.clear();
      tang_grad_deriv_.clear();
      weighted_relTangVel_.clear();
      weighted_avTangVel_.clear();
    }

    Core::LinAlg::Matrix<3, 1>& get_weighted_rel_tang_vel() { return weighted_relTangVel_; }
    std::unordered_map<int, Core::LinAlg::Matrix<3, 1>>& get_weighted_rel_tang_vel_deriv()
    {
      return deriv_weighted_relTangVel_;
    }

    Core::LinAlg::Matrix<3, 1>& get_weighted_av_tang_vel() { return weighted_avTangVel_; }
    std::unordered_map<int, Core::LinAlg::Matrix<3, 1>>& get_weighted_av_tang_vel_deriv()
    {
      return deriv_weighted_avTangVel_;
    }

    std::unordered_map<int, Core::LinAlg::Matrix<3, 1>>& get_surf_grad() { return tang_grad_; }
    std::unordered_map<int, std::unordered_map<int, Core::LinAlg::Matrix<3, 1>>>&
    get_surf_grad_deriv()
    {
      return tang_grad_deriv_;
    }

   protected:
    // don't want = operator and cctor
    NodeEhlDataContainer operator=(const NodeEhlDataContainer& old) = delete;
    NodeEhlDataContainer(const NodeEhlDataContainer& old) = delete;

    // actual data
    Core::LinAlg::Matrix<3, 1> weighted_relTangVel_;
    std::unordered_map<int, Core::LinAlg::Matrix<3, 1>> deriv_weighted_relTangVel_;

    Core::LinAlg::Matrix<3, 1> weighted_avTangVel_;
    std::unordered_map<int, Core::LinAlg::Matrix<3, 1>> deriv_weighted_avTangVel_;

    std::unordered_map<int, Core::LinAlg::Matrix<3, 1>> tang_grad_;
    std::unordered_map<int, std::unordered_map<int, Core::LinAlg::Matrix<3, 1>>> tang_grad_deriv_;
  };

  /*!
   \brief A class for a contact node derived from Mortar::Node

   This class represents a finite element node capable of contact.

   */
  class Node : public Mortar::Node
  {
   public:
    //! @name Enums and Friends

    /*!
     \brief The discretization is a friend of Node
     */
    friend class Core::FE::Discretization;

    //@}

    //! @name Constructors and destructors and related methods

    /*!
     \brief Standard Constructor

     \param id     (in): A globally unique node id
     \param coords (in): vector of nodal coordinates, length 3
     \param owner  (in): Owner of this node.
     \param dofs   (in): list of global degrees of freedom
     \param isslave(in): flag indicating whether node is slave or master
     \param initactive (in): flag indicating whether initially set to active

     */
    Node(int id, const std::vector<double>& coords, const int owner, const std::vector<int>& dofs,
        const bool isslave, const bool initactive);

    /*!
     \brief Copy Constructor

     Makes a deep copy of a Node

     */
    Node(const CONTACT::Node& old);

    /*!
     \brief Deep copy the derived class and return pointer to it

     */
    CONTACT::Node* clone() const override;


    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of lib/parobject.H.

     */
    int unique_par_object_id() const override
    {
      return NodeType::instance().unique_par_object_id();
    }

    /*!
     \brief Pack this class so it can be communicated

     \ref pack and \ref unpack are used to communicate this node

     */
    void pack(Core::Communication::PackBuffer& data) const override;

    /*!
     \brief Unpack data from a char vector into this class

     \ref pack and \ref unpack are used to communicate this node

     */
    void unpack(Core::Communication::UnpackBuffer& buffer) override;

    //@}

    //! @name Access methods

    /*!
     \brief Print this contact node
     */
    void print(std::ostream& os) const override;

    /*!
     \brief Is Node initialized as active node (only slave nodes)
     */
    virtual bool is_init_active() const
    {
      if (!is_slave()) FOUR_C_THROW("InitActive requested for Master node");
      return initactive_;
    }

    /*!
     \brief Modify initial active status of slave node

     This belated modification is necessary to be able to use
     the binary search tree for contact initialization in the
     load-controlled quasi-static case (instead of input file
     information Active/Inactive)

     */
    virtual bool& set_init_active()
    {
      if (!is_slave()) FOUR_C_THROW("InitActive requested for Master node");
      return initactive_;
    }

    /*!
     \brief Return contact status of this node (active=true)
     */
    virtual bool& active() { return active_; }
    bool active() const { return active_; }

    virtual bool& involved_m() { return involvedm_; }

    /*!
     \brief Return data container of this node

     This method returns the data container of this node where additional
     contact specific quantities/information are stored.

     */
    inline CONTACT::NodeDataContainer& data() { return *codata_; }
    inline const CONTACT::NodeDataContainer& data() const { return *codata_; }

    inline CONTACT::NodePoroDataContainer& poro_data() { return *coporodata_; }
    inline const CONTACT::NodePoroDataContainer& poro_data() const { return *coporodata_; }
    inline CONTACT::NodeTSIDataContainer& tsi_data() { return *cTSIdata_; }
    inline const CONTACT::NodeTSIDataContainer& tsi_data() const { return *cTSIdata_; }
    inline CONTACT::NodeEhlDataContainer& ehl_data() { return *cEHLdata_; }
    inline const CONTACT::NodeEhlDataContainer& ehl_data() const { return *cEHLdata_; }

    //@}

    //! @name Evaluation methods

    /*!
     \brief Add a value to the weighted gap of this node

     This value is later assembled to the weighted gap vec.
     Note that grow_ here is a scalar.

     \param val : value to be added

     */
    void addg_value(const double val);

    /*!
     \brief Add a value to the point-wise gap of this node (NTS)

     \param val : value to be added

     */
    void addnts_gap_value(const double val);

    /*!
     \brief Add a value to the line-weighted gap of this node (LTS)

     \param val : value to be added

     */
    void addlts_gap_value(const double val);

    /*!
     \brief Add a value to the point-wise gap of this node (LTL)

     \param val : value to be added

     */
    void addltl_gap_value(const double* val);

    /*!
     \brief Add a value to the point-wise jump of this node (LTL)

     \param val : value to be added

     */
    void addltl_jump_value(const double* val);


    /*!
    \brief Add a value to the map of LM derivatives of this node

    The 'DerivZ' map is later assembled to the global DerivZ matrix.
    Note that derivz_ here is a vector.

    Note: This is only calculated when performing a penalty strategy

    \param row : local dof row id to add to (row-wise)
    \param val : value to be added
    \param col : global dof column id of the value added

    */
    void add_deriv_z_value(const int row, const int col, const double val);

    /*!
    \brief Add a value to the NCoup of this node

    \param val : value to be added

    */
    void add_ncoup_value(const double val);

    /*!
     \brief Build nodal normal
     */
    void build_averaged_normal() override;

    /*!
     \brief Build nodal edge tangent
     */
    void build_averaged_edge_tangent();

    /*!
     \brief Initializes the data container of the node

     With this function, the container with contact specific quantities/information
     is initialized.

     */
    void initialize_data_container() override;


    /*!
    \brief Initializes the poro data container of the node

    With this function, the container with contact specific quantities/information
    is initialized. --- Used to initialize PoroDataContainer for master nodes!

    */
    void initialize_poro_data_container() override;

    /*!
    \brief Initializes the ehl data container of the node

    With this function, the container with contact specific quantities/information
    is initialized.

    */
    void initialize_ehl_data_container() override;

    /*!
    \brief Initializes the TSI data container of the node
    */
    virtual void initialize_tsi_data_container(double t_ref, double t_dam);

    /*!
     \brief Resets the data container of the node

     With this function, the container with contact specific quantities/information
     is deleted / reset to nullptr pointer

     */
    void reset_data_container() override;

    /*!
     \brief Get number of linearization entries

     */
    virtual int& get_linsize() { return linsize_; };
    int get_linsize() const { return linsize_; };

    //! @}

    /*! @name Empty functions (friction only)
     *
     * All these functions only have functionality for friction nodes, thus they are
     * defined as empty here in the general mortar node. They can be called whenever you like.
     */

    virtual void add_s_node(int node) {}
    virtual void add_m_node(int node) {}

    //! @}

    virtual bool has_tsi_data() { return (cTSIdata_ != nullptr); }

    /*!
     \brief Write nodal normals to old ones
     */
    void store_old_normal();

   private:
    /*!
     \brief Build directional derivative of nodal normal + tangents

     This method will be called after having finished the method
     BuildAveragedNormal() for the computation of nodal normals. The result
     (directional derivatives) will be stored in the nodal maps derivn_ and
     will later be assembled to the global system of equations.
     Please note that we also compute the directional derivative of the
     nodal tangent here,and the results will be stored analogously in
     the nodal maps derivtxi_ (2D and 3D) and derivteta_ (only 3D).

     \param elens (in):  Matrix containing normals of adjacent elements
     \param length (in): Length of the nodal averaged normal
     \param ltxi (in):   Length of the nodal tangent txi

     */
    void deriv_averaged_normal(Core::LinAlg::SerialDenseMatrix& elens, double length, double ltxi);

   protected:
    //! true if contact node is in contact (active set strategy)
    bool active_;

    //! true if node is initialized as active node
    bool initactive_;

    //! cnode is an master node for which the mortar integration is performed (both-sided wear)
    bool involvedm_;

    //! number of lin entries per direction vector
    int linsize_;

    //! Additional information of proc's contact nodes
    std::shared_ptr<CONTACT::NodeDataContainer> codata_;

    //! Additional information of proc's poro contact nodes
    std::shared_ptr<CONTACT::NodePoroDataContainer> coporodata_;

    //! Additional information of TSI contact nodes
    std::shared_ptr<CONTACT::NodeTSIDataContainer> cTSIdata_;

    //! Additional information of EHL contact nodes
    std::shared_ptr<CONTACT::NodeEhlDataContainer> cEHLdata_;
  };
  // class Node
}  // namespace CONTACT

// << operator
std::ostream& operator<<(std::ostream& os, const CONTACT::Node& cnode);

FOUR_C_NAMESPACE_CLOSE

#endif
