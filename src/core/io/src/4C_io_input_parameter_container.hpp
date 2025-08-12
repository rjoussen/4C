// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_INPUT_PARAMETER_CONTAINER_HPP
#define FOUR_C_IO_INPUT_PARAMETER_CONTAINER_HPP


#include "4C_config.hpp"

#include "4C_io_input_types.hpp"
#include "4C_utils_demangle.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <algorithm>
#include <any>
#include <functional>
#include <map>
#include <optional>
#include <ostream>
#include <ranges>
#include <string>
#include <typeindex>
#include <vector>


FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  /**
   * A container to store dynamic input parameters. The container can store arbitrary types of
   * parameters. Parameters can be grouped in sub-containers.
   *
   * This class is a core part of the input mechanism, as it contains the parsed data from the input
   * file and grants access to it.
   */
  class InputParameterContainer
  {
   public:
    /**
     * A type to store a list of InputParameterContainers. This type is used to represent what is
     * often called a *list*, *array*, or *sequence* of data in the input file.
     *
     * @note This type is called `List` to more clearly distinguish it from a simple `std::vector`
     * entry in the container. A `List` contains nested InputParameterContainers and encodes
     * rather complex data structures. Nevertheless, it is implemented as a `std::vector`.
     */
    using List = std::vector<InputParameterContainer>;

    using GroupStorage = std::map<std::string, InputParameterContainer>;

    /**
     * \brief Add @data to the container at the given key @name.
     *
     * If an entry with given @p name already exists, it will be overwritten. The type must be
     * one of the SupportedType of the input framework.
     */
    template <typename T>
    void add(const std::string& name, const T& data);

    /**
     * Access group @p name. If the group does not exist, it will be created.
     */
    InputParameterContainer& group(const std::string& name);

    /**
     * Access group @p name. This function throws an error if the group does not exist.
     */
    [[nodiscard]] const InputParameterContainer& group(const std::string& name) const;

    /**
     * Get a range view of all groups in the container.
     */
    [[nodiscard]] auto groups() const;

    /**
     * Check if a group with the given @p name exists.
     */
    [[nodiscard]] bool has_group(const std::string& name) const;

    /**
     * Ensure that exactly one group having one of the @p possible_group_names is present in the
     * container and return that group. Throws if there are no or multiple matching groups.
     */
    template <std::ranges::range R>
      requires std::convertible_to<std::ranges::range_value_t<R>, std::string>
    [[nodiscard]] std::pair<std::string_view, const InputParameterContainer&> exactly_one_group(
        const R& possible_group_names) const;

    /**
     * Add the list @p data at the given key @p name.
     *
     * @note This functions is a more obvious way to add a list to the container compared to
     * the add() function with a List template argument, although this is precisely what happens
     * internally.
     */
    void add_list(const std::string& name, List&& list);

    /**
     * Access the list @p name. This function throws an error if the list does not exist.
     *
     * @note This functions is a more obvious way to get a list from the container compared to
     * the get() function with a List template argument, although this is precisely
     * what happens internally.
     */
    [[nodiscard]] const List& get_list(const std::string& name) const;

    /**
     * Combine the data from another container with this one. Conflicting data will throw an
     * error.
     */
    void merge(const InputParameterContainer& other);

    /*!
     * Get a const reference to the data stored at the key @p name from the container. An error
     * is thrown in case no value of specified type is stored under @p name in the container.
     */
    template <typename T>
    const T& get(const std::string& name) const;

    /*!
     * Get the data stored at the key @p name from the container. Return the @p default_value if
     * no value of specified type is stored under @p name in the container.
     *
     * @note This function returns the value as a copy.
     */
    template <typename T>
    T get_or(const std::string& name, T default_value) const;

    /*!
     * Get the const pointer to the data stored at the key @p name from the container. Return a
     * `nullptr` if no value of specified type is stored under @p name in the container.
     */
    template <typename T>
    const T* get_if(const std::string& name) const;

    /**
     * Print the data in the container to the output stream @p os.
     */
    void print(std::ostream& os) const;

    /**
     * Clear the container.
     */
    void clear();

    /**
     * Convert the data in this container to a Teuchos::ParameterList. All groups are converted to
     * sublists.
     */
    void to_teuchos_parameter_list(Teuchos::ParameterList& list) const;

   private:
    //! Entry stored in the container.
    struct Entry
    {
      //! The actual data.
      std::any data;
    };

    //! Data stored in this container.
    std::map<std::string, Entry> entries_;

    //! Groups present in this container. Groups are InputParameterContainers themselves.
    std::map<std::string, InputParameterContainer> groups_;

    /**
     * Gather different actions that can be performed on a type.
     */
    struct TypeActions
    {
      //! The function to print the data.
      std::function<void(std::ostream&, const std::any&)> print;

      //! Function to write the data into a Teuchos::ParameterList.
      std::function<void(Teuchos::ParameterList&, const std::string& name, const std::any&)>
          write_to_pl;
    };

    /**
     * Add the type actions if not already present.
     */
    template <typename T>
    void ensure_type_action_registered();

    /**
     * Access the shared storage for all the type actions.
     */
    static std::map<std::type_index, TypeActions>& get_type_actions();
  };  // class InputParameterContainer

}  // namespace Core::IO


// --- template and inline functions ---//


namespace Core::IO::Internal::InputParameterContainerImplementation
{
  template <typename T>
  const T* try_get_any_data(const std::string& name, const std::any& data)
  {
    if (typeid(T) == data.type())
    {
      const T* any_ptr = std::any_cast<T>(&data);
      FOUR_C_ASSERT(any_ptr != nullptr, "Implementation error.");
      return any_ptr;
    }
    else
    {
      FOUR_C_THROW(
          "You tried to get the data named {} from the container as type '{}'.\n"
          "Actually, it has type '{}'.",
          name.c_str(), Core::Utils::get_type_name<T>().c_str(),
          Core::Utils::try_demangle(data.type().name()).c_str());
    }
  }
}  // namespace Core::IO::Internal::InputParameterContainerImplementation


inline auto Core::IO::InputParameterContainer::groups() const
{
  return std::ranges::subrange(groups_);
}


template <std::ranges::range R>
  requires std::convertible_to<std::ranges::range_value_t<R>, std::string>
std::pair<std::string_view, const Core::IO::InputParameterContainer&>
Core::IO::InputParameterContainer::exactly_one_group(const R& possible_group_names) const
{
  auto matching_group_names =
      possible_group_names |
      std::views::filter([this](const std::string& name) { return has_group(name); });

  FOUR_C_ASSERT_ALWAYS(std::ranges::distance(matching_group_names) == 1,
      "The data container must contain exactly one group that matches a physics name in the "
      "element definitions but found {} matching groups.",
      std::ranges::distance(matching_group_names));

  return {*matching_group_names.begin(), group(*matching_group_names.begin())};
}


template <typename T>
const T& Core::IO::InputParameterContainer::get(const std::string& name) const
{
  if (const T* p = get_if<T>(name))
    return *p;
  else
    FOUR_C_THROW("Key '{}' cannot be found in the container.", name);
}

template <typename T>
T Core::IO::InputParameterContainer::get_or(const std::string& name, T default_value) const
{
  if (const T* p = get_if<T>(name))
    return *p;
  else
    return default_value;
}

template <typename T>
const T* Core::IO::InputParameterContainer::get_if(const std::string& name) const
{
  const auto it = entries_.find(name);
  if (it != entries_.end())
  {
    return Internal::InputParameterContainerImplementation::try_get_any_data<T>(
        name, it->second.data);
  }
  else
  {
    return nullptr;
  }
}

FOUR_C_NAMESPACE_CLOSE

#endif
