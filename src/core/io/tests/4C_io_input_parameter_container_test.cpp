// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_io_input_parameter_container.templates.hpp"

#include <Teuchos_ParameterList.hpp>

#include <ostream>

namespace
{
  using namespace FourC;

  using namespace Core::IO;

  TEST(InputParameterContainerTest, ToTeuchosParameterList)
  {
    InputParameterContainer container;
    container.add<int>("a", 1);
    container.add<double>("b", 2.0);
    container.add<std::string>("c", "string");
    container.add<bool>("d", true);
    container.add("v", std::vector<int>{1, 2, 3});
    container.add("n", std::optional<int>{});
    container.group("group").add<std::string>("e", "group string");
    container.group("group").group("deeply").group("nested").add<int>("f", 42);

    Teuchos::ParameterList list("dummy");
    container.to_teuchos_parameter_list(list);
    EXPECT_EQ(list.name(), "dummy");
    EXPECT_EQ(list.get<int>("a"), 1);
    EXPECT_EQ(list.get<double>("b"), 2.0);
    EXPECT_EQ(list.get<std::string>("c"), "string");
    EXPECT_EQ(list.get<bool>("d"), true);
    EXPECT_EQ(list.get<std::vector<int>>("v"), (std::vector<int>{1, 2, 3}));
    EXPECT_EQ(list.get<std::optional<int>>("n"), std::nullopt);
    EXPECT_EQ(list.sublist("group").get<std::string>("e"), "group string");
    EXPECT_EQ(list.sublist("group").sublist("deeply").sublist("nested").get<int>("f"), 42);
  }

  TEST(InputParameterContainerTest, Lists)
  {
    InputParameterContainer container;
    InputParameterContainer::List list;

    list.emplace_back().add("a", 1);
    auto& list2 = list.emplace_back();
    list2.group("group").add("b", 2);

    container.add_list("list", std::move(list));

    Teuchos::ParameterList pl;
    container.to_teuchos_parameter_list(pl);

    EXPECT_EQ(
        pl.get<std::vector<Teuchos::ParameterList>>("list")[1].sublist("group").get<int>("b"), 2);
  }

  TEST(InputParameterContainerTest, PrintEnum)
  {
    enum class TestEnum
    {
      A,
      B,
      C
    };

    InputParameterContainer container;
    container.add("test_enum", TestEnum::A);

    // check print output
    std::ostringstream print_output{""};
    container.print(print_output);
    EXPECT_EQ(print_output.str(), "test_enum : A ");
  }

  TEST(InputParameterContainerTest, ExactlyOneGroup)
  {
    InputParameterContainer container;
    container.group("group1").add<int>("a", 1);
    container.group("group2").add<int>("b", 2);
    container.group("group3").add<int>("c", 3);

    {
      SCOPED_TRACE("Exactly one");
      std::array possible_names = {"group1", "not existing", ""};
      auto [name, group] = container.exactly_one_group(possible_names);
      EXPECT_EQ(name, "group1");
      EXPECT_EQ(group.get<int>("a"), 1);
    }

    {
      SCOPED_TRACE("Too many");
      std::vector<std::string> possible_names = {"group1", "group2", "not existing", ""};
      EXPECT_THROW(
          [[maybe_unused]] auto c = container.exactly_one_group(possible_names), Core::Exception);
    }

    {
      SCOPED_TRACE("Too few");
      std::array possible_names = {"not existing", ""};
      EXPECT_THROW(
          [[maybe_unused]] auto c = container.exactly_one_group(possible_names), Core::Exception);
    }
  }
}  // namespace