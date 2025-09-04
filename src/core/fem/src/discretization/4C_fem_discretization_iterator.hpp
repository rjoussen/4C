// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_DISCRETIZATION_ITERATOR_HPP
#define FOUR_C_FEM_DISCRETIZATION_ITERATOR_HPP

#include "4C_config.hpp"

#include <concepts>
#include <iterator>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;

  /**
   * Check whether T is a type that can be constructed from a Discretization pointer and an int.
   */
  template <typename T>
  concept DiscretizationReferenceType = std::constructible_from<T, Discretization*, int> ||
                                        std::constructible_from<T, const Discretization*, int>;

  /**
   * @brief A proxy object that behaves like a reference for iterator dereferencing.
   */
  template <typename RefType>
  class IteratorProxy
  {
   public:
    explicit IteratorProxy(RefType&& value) : value_(std::move(value)) {}

    const RefType* operator->() const { return &value_; }
    RefType* operator->() { return &value_; }

    const RefType& operator*() const { return value_; }
    RefType& operator*() { return value_; }

   private:
    RefType value_;
  };

  /**
   * @brief An iterator over entities in a Discretization.
   *
   * This iterator supports iteration over a contiguous index array. Dereferencing the iterator
   * returns an object of type RefType, referring to the entity at the current index.
   * The iterator supports iteration over const and non-const Discretization objects.
   */
  template <typename RefType>
  class DiscretizationIterator
  {
   public:
    //! The type of internal index that is used to construct RefType objects
    using IndexType = int;
    using DiscretizationType = typename RefType::DiscretizationType;

    using iterator_category = std::forward_iterator_tag;
    using value_type = RefType;
    using difference_type = std::ptrdiff_t;
    using reference = value_type;
    using pointer = IteratorProxy<value_type>;

    DiscretizationIterator(DiscretizationType* discretization, const IndexType* index_ptr)
        : discretization_(discretization), index_ptr_(index_ptr)
    {
    }

    reference operator*() const { return RefType(discretization_, *index_ptr_); }

    pointer operator->() const { return pointer(RefType(discretization_, *index_ptr_)); }

    DiscretizationIterator& operator++()
    {
      ++index_ptr_;
      return *this;
    }

    DiscretizationIterator operator++(int)
    {
      DiscretizationIterator tmp = *this;
      ++(*this);
      return tmp;
    }


    bool operator==(const DiscretizationIterator& other) const
    {
      return (discretization_ == other.discretization_) && (index_ptr_ == other.index_ptr_);
    }

    bool operator!=(const DiscretizationIterator& other) const { return !(*this == other); }

   private:
    //! Discretization containing the entities being iterated over
    DiscretizationType* discretization_;
    //! Pointer into the index array that we iterate over
    const IndexType* index_ptr_;
  };

  /**
   * A range of Iterator objects for use in range-based expressions.
   */
  template <typename Iterator>
  class IteratorRange
  {
   public:
    IteratorRange(Iterator begin, Iterator end) : begin_(begin), end_(end) {}

    Iterator begin() const { return begin_; }
    Iterator end() const { return end_; }

   private:
    Iterator begin_;
    Iterator end_;
  };
}  // namespace Core::FE

FOUR_C_NAMESPACE_CLOSE

#endif
