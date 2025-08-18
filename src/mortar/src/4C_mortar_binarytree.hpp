// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MORTAR_BINARYTREE_HPP
#define FOUR_C_MORTAR_BINARYTREE_HPP

#include "4C_config.hpp"

#include "4C_inpar_mortar.hpp"
#include "4C_linalg_map.hpp"
#include "4C_mortar_base_binarytree.hpp"

#include <mpi.h>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Mortar
{
  // forward declarations

  //! @name Enums and Friends

  /// Type of binary tree node
  enum BinaryTreeNodeType
  {
    SLAVE_INNER,        ///< indicates a slave inner node (has children)
    SLAVE_LEAF,         ///< indicates a slave leaf node (no further children)
    MASTER_INNER,       ///< indicates a master inner node (has children)
    MASTER_LEAF,        ///< indicates a master leaf node (no further children)
    NOSLAVE_ELEMENTS,   ///< indicates that there are no slave elements on this (root) treenode
    NOMASTER_ELEMENTS,  ///< indicates that there are no master elements on this (root) treenode
    UNDEFINED           ///< indicates an undefined tree node
  };

  //@}

  /*!
  \brief A class representing one tree node of the binary search tree

  Refer also to the Semesterarbeit of Thomas Eberl, 2009

  */
  class BinaryTreeNode : public BaseBinaryTreeNode
  {
   public:
    /*!
    \brief constructor of a tree node

    \param type        type of BinaryTreeNode
    \param discret     interface discretization
    \param parent      points to parent tree node
    \param elelist     list of all elements in BinaryTreeNode
    \param dopnormals  reference to DOP normals
    \param kdop        reference to no. of vertices
    \param dim         dimension of problem
    \param useauxpos   bool whether auxiliary position is used when computing dops
    \param layer       current layer of tree node
    \param ...map      references to maps

    */
    BinaryTreeNode(BinaryTreeNodeType type, Core::FE::Discretization& discret,
        std::shared_ptr<BinaryTreeNode> parent, std::vector<int> elelist,
        const Core::LinAlg::SerialDenseMatrix& dopnormals, const int& kdop, const int& dim,
        const bool& useauxpos, const int layer,
        std::vector<std::vector<std::shared_ptr<BinaryTreeNode>>>& streenodesmap,
        std::vector<std::vector<std::shared_ptr<BinaryTreeNode>>>& mtreenodesmap,
        std::vector<std::vector<std::shared_ptr<BinaryTreeNode>>>& sleafsmap,
        std::vector<std::vector<std::shared_ptr<BinaryTreeNode>>>& mleafsmap);


    //! @name Evaluation methods

    /*!
    \brief Update slabs of current treenode in bottom up way

    */
    void update_slabs_bottom_up(double& enlarge) final;

    /*!
    \brief Initialize Tree

    */
    void initialize_tree(double& enlarge);

    /*!
    \brief Divide a TreeNode into two child nodes

    */
    void divide_tree_node();

    /*!
    \brief Print type of tree node to std::cout

    */
    void print_type() final;
    //@}

    //! @name Access and modification methods

    /*!
    \brief Get communicator

    */
    MPI_Comm get_comm() const;

    /*!
    \brief Return type of treenode

    */
    BinaryTreeNodeType type() const { return type_; }

    /*!
    \brief Return pointer to right child

    */
    std::shared_ptr<BinaryTreeNode> rightchild() const { return rightchild_; }

    /*!
    \brief Return pointer to left child

    */
    std::shared_ptr<BinaryTreeNode> leftchild() const { return leftchild_; }
    //@}

   private:
    // don't want = operator and cctor
    BinaryTreeNode operator=(const BinaryTreeNode& old);
    BinaryTreeNode(const BinaryTreeNode& old);

    //! type of BinaryTreeNode
    Mortar::BinaryTreeNodeType type_;

    // the pointers to the parent as well as to the left and right child are not moved to the
    // BaseBinaryTreeNode as this would require a lot of dynamic casting and thereby complicating
    // the readability of the code
    //! pointer to the parent BinaryTreeNode
    std::shared_ptr<BinaryTreeNode> parent_;

    //! pointer to the left child TreeNode
    std::shared_ptr<BinaryTreeNode> leftchild_;

    //! pointer to the right child TreeNode
    std::shared_ptr<BinaryTreeNode> rightchild_;

    //! reference to map of all slave treenodes, sorted by layer
    std::vector<std::vector<std::shared_ptr<BinaryTreeNode>>>& streenodesmap_;

    //! reference to map of all master treenodes, sorted by layer
    std::vector<std::vector<std::shared_ptr<BinaryTreeNode>>>& mtreenodesmap_;

    //! reference to map of all slave leaf treenodes
    std::vector<std::vector<std::shared_ptr<BinaryTreeNode>>>& sleafsmap_;

    //! reference to map of all master leaf treenodes
    std::vector<std::vector<std::shared_ptr<BinaryTreeNode>>>& mleafsmap_;

  };  // class BinaryTreeNode


  /*!
  \brief A class for performing search in 2D/3D based on binary trees

  Refer also to the Semesterarbeit of Thomas Eberl, 2009

  */
  class BinaryTree : public BaseBinaryTree
  {
   public:
    /*!
    \brief Standard constructor

    Constructs an instance of this class.<br>

    \param discret (in):     The interface discretization
    \param selements (in):   All slave elements (column map)
    \param melements (in):   All master elements (fully overlapping map)
    \param dim (in):         The problem dimension
    \param eps (in):         factor used to enlarge dops
    \param updatetype (in):  Defining type of binary tree update (top down, or bottom up)
    \param useauxpos (in):   flag indicating usage of auxiliary position for calculation of slabs

    */
    BinaryTree(Core::FE::Discretization& discret, std::shared_ptr<Core::LinAlg::Map> selements,
        std::shared_ptr<Core::LinAlg::Map> melements, int dim, double eps,
        Inpar::Mortar::BinaryTreeUpdateType updatetype, bool useauxpos);


    //! @name Query methods

    /*!
    \brief Evaluate search tree to get corresponding master elements for the slave elements

    */
    void evaluate_search() final;

    /*!
    \brief Initialize the binary tree

    */
    void init() final;

   private:
    /*!
    \brief clear found search elements

    */
    void init_search_elements();

    /*!
    \brief Print full tree

    */
    void print_tree(BinaryTreeNode& treenode);

    /*!
    \brief Print full tree out of map of treenodes

    */
    void print_tree_of_map(std::vector<std::vector<std::shared_ptr<BinaryTreeNode>>>& treenodesmap);

    //@}

    //! @name Access methods

    /*!
    \brief Get communicator

    */
    MPI_Comm get_comm() const;

    /*!
    \brief Return reference to slave treenodesmap

    */
    std::vector<std::vector<std::shared_ptr<BinaryTreeNode>>>& streenodesmap()
    {
      return streenodesmap_;
    }

    /*!
    \brief Return reference to master treenodesmap

    */
    std::vector<std::vector<std::shared_ptr<BinaryTreeNode>>>& mtreenodesmap()
    {
      return mtreenodesmap_;
    }

    /*!
    \brief Return reference to coupling treenodesmap

    */
    std::vector<std::vector<std::shared_ptr<BinaryTreeNode>>>& coupling_map()
    {
      return couplingmap_;
    }

    /*!
    \brief Return pointer to sroot-treenode

    */
    std::shared_ptr<BinaryTreeNode>& sroot() { return sroot_; }
    //@}

    //! @name Evaluation methods

    /*!
    \brief Initialize internal variables

     */
    void init_internal_variables() final;

    /*!
    \brief Calculate minimal element length / inflation factor "enlarge"

    */
    void set_enlarge() final;

    /*!
    \brief Update master and slave tree in a top down way

    */
    void update_tree_top_down()
    {
      evaluate_update_tree_top_down(*sroot_);
      evaluate_update_tree_top_down(*mroot_);
      return;
    }

    /*!
    \brief Evaluate update of master and slave tree in a top down way

    */
    void evaluate_update_tree_top_down(BinaryTreeNode& treenode);

    /*!
    \brief Updates master and slave tree in a bottom up way

    */
    void update_tree_bottom_up()
    {
      evaluate_update_tree_bottom_up(streenodesmap_);
      evaluate_update_tree_bottom_up(mtreenodesmap_);
      return;
    }

    /*!
    \brief Evaluate update of master and slave tree in a bottom up way

    */
    void evaluate_update_tree_bottom_up(
        std::vector<std::vector<std::shared_ptr<BinaryTreeNode>>>& treenodesmap);

    /*!
    \brief Evaluate binary search tree

    \note Search and update is carried out in a separate way. There has also been a combined
    approach, but this has been removed as part of GitLab Issue 181, as it is outperformed by the
    separate approach for large problems!

    */
    void evaluate_search(
        std::shared_ptr<BinaryTreeNode> streenode, std::shared_ptr<BinaryTreeNode> mtreenode);

    // don't want = operator and cctor
    BinaryTree operator=(const BinaryTree& old);
    BinaryTree(const BinaryTree& old);

    //! all slave elements on surface (column map)
    std::shared_ptr<Core::LinAlg::Map> selements_;
    //! all master elements on surface (full map)
    std::shared_ptr<Core::LinAlg::Map> melements_;
    //! map of all slave tree nodes, sorted by layers
    std::vector<std::vector<std::shared_ptr<BinaryTreeNode>>> streenodesmap_;
    //! map of all master tree nodes, sorted by layers
    std::vector<std::vector<std::shared_ptr<BinaryTreeNode>>> mtreenodesmap_;
    //! map of all tree nodes, that possibly couple, st/mt
    std::vector<std::vector<std::shared_ptr<BinaryTreeNode>>> couplingmap_;
    //! map of all slave leaf tree nodes, [0]=left child,[1]=right child
    std::vector<std::vector<std::shared_ptr<BinaryTreeNode>>> sleafsmap_;
    //! map of all master leaf tree nodes, [0]=left child,[1]=right child
    std::vector<std::vector<std::shared_ptr<BinaryTreeNode>>> mleafsmap_;
    //! slave root tree node
    std::shared_ptr<BinaryTreeNode> sroot_;
    //! master root tree node
    std::shared_ptr<BinaryTreeNode> mroot_;
    //! update type of binary tree
    const Inpar::Mortar::BinaryTreeUpdateType updatetype_;
    //! bool whether auxiliary position is used when computing dops
    bool useauxpos_;
  };  // class BinaryTree
}  // namespace Mortar

FOUR_C_NAMESPACE_CLOSE

#endif
