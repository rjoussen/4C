// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_dofset_base.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

#include <algorithm>
#include <iostream>

FOUR_C_NAMESPACE_OPEN


// list of all dof sets
std::list<Core::DOFSets::DofSetInterface*> Core::DOFSets::DofSetBase::static_dofsets_;

/*----------------------------------------------------------------------*
 |  ctor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
Core::DOFSets::DofSetBase::DofSetBase() : DofSetInterface() { return; }


/*----------------------------------------------------------------------*
 |  dtor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
Core::DOFSets::DofSetBase::~DofSetBase()
{
  // disconnect from proxies
  for (std::list<DofSetInterface*>::iterator i = registered_dofsets_.begin();
      i != registered_dofsets_.end(); ++i)
  {
    (*i)->disconnect(this);
  }
  // remove dofset from static list if necessary
  static_dofsets_.remove(this);
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetBase::add_dof_set_to_list()
{
  if (std::find(static_dofsets_.begin(), static_dofsets_.end(), this) == static_dofsets_.end())
  {
    static_dofsets_.push_back(this);
  }
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetBase::replace_in_static_dofsets(
    std::shared_ptr<DofSetInterface> olddofset)
{
  std::list<Core::DOFSets::DofSetInterface*>::iterator iterold =
      std::find(static_dofsets_.begin(), static_dofsets_.end(), &(*olddofset));
  if (iterold == static_dofsets_.end())
  {
    static_dofsets_.push_back(this);
  }
  else
  {
    static_dofsets_.insert(iterold, this);
    static_dofsets_.remove(&(*olddofset));
  }
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Core::DOFSets::DofSetBase::max_gid_in_list(MPI_Comm comm) const
{
  int count = -1;
  for (std::list<DofSetInterface*>::const_iterator i = static_dofsets_.begin();
      i != static_dofsets_.end(); ++i)
  {
    if (*i == this) break;

    // ignore empty (no yet initialized) dof row maps
    if ((*i)->initialized())
    {
      count = std::max((*i)->max_all_gid(), count);
    }
  }
  int max;
  max = Core::Communication::max_all(count, comm);
  return max;
}

void Core::DOFSets::DofSetBase::print_all_dofsets(MPI_Comm comm) const
{
  if (Core::Communication::my_mpi_rank(comm) == 0)
  {
    std::vector<int> min;
    std::vector<int> max;
    for (std::list<DofSetInterface*>::const_iterator i = static_dofsets_.begin();
        i != static_dofsets_.end(); ++i)
    {
      min.push_back((*i)->min_all_gid());
      max.push_back((*i)->max_all_gid());
    }
    if (min.size() < 1) return;

    std::vector<int>::const_iterator largest = max_element(max.begin(), max.end());
    int allmax = *largest;
    int availspace = 80;
    int arrowlen = availspace + 3;

    if (min[0] > 0) availspace -= 2;

    for (unsigned int i = 0; i < min.size(); ++i)
    {
      // left bar
      availspace--;
      // number
      availspace -= 12;
      // right bar
      if (i + 1 < min.size())
      {
        if (max[i] + 1 < min[i + 1]) availspace -= 2;
      }
      else
        availspace--;
    }

    if (availspace < 0)
    {
      arrowlen -= availspace;
      availspace = 0;
    }

    Core::IO::cout << "All registered dofsets:" << Core::IO::endl;

    if (min[0] > 0) Core::IO::cout << "| ";
    for (unsigned int i = 0; i < min.size(); ++i)
    {
      // left bar
      if (i > 0 and max[i - 1] > min[i])
        Core::IO::cout << "X";
      else
        Core::IO::cout << "|";

      // number
      Core::IO::cout << std::left << std::setw(6) << min[i];
      for (int j = 0; j < (int)(max[i] - min[i]) * 1.0 * availspace / allmax; ++j)
        Core::IO::cout << " ";
      Core::IO::cout << std::right << std::setw(6) << max[i];

      // right bar
      if (i + 1 < min.size())
      {
        if (max[i] + 1 < min[i + 1]) Core::IO::cout << "| ";
      }
      else
      {
        Core::IO::cout << "|";
      }
    }

    Core::IO::cout << Core::IO::endl << "+";
    for (int i = 0; i < arrowlen; ++i) Core::IO::cout << "-";
    Core::IO::cout << "> DofGID" << Core::IO::endl;
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetBase::register_proxy(DofSetInterface* dofset)
{
  registered_dofsets_.push_back(dofset);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetBase::unregister(DofSetInterface* dofset)
{
  registered_dofsets_.remove(dofset);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetBase::notify_assigned()
{
  for (std::list<DofSetInterface*>::iterator i = registered_dofsets_.begin();
      i != registered_dofsets_.end(); ++i)
    (*i)->notify_assigned();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetBase::notify_reset()
{
  for (std::list<DofSetInterface*>::iterator i = registered_dofsets_.begin();
      i != registered_dofsets_.end(); ++i)
    (*i)->notify_reset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetBase::print(std::ostream& os) const
{
  FOUR_C_THROW("print() is not implemented in base class. Override print() in subclass");
}

FOUR_C_NAMESPACE_CLOSE
