// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_io_input_spec_builders.hpp"
#include "4C_io_input_spec_builders.templates.hpp"
#include "4C_unittest_utils_assertions_test.hpp"


// Contains tests related to legacy dat-style input

namespace
{
  using namespace FourC;
  using namespace Core::IO;
  using namespace Core::IO::InputSpecBuilders;

  TEST(InputSpecLegacyTest, Simple)
  {
    auto spec = all_of({
        parameter<int>("a", {.description = "An integer", .default_value = 1}),
        parameter<double>("b"),
        parameter<bool>("d"),
    });
    InputParameterContainer container;
    std::string stream("b 2.0 d true // trailing comment");
    ValueParser parser(stream);
    spec.deprecated_parse(parser, container);
    EXPECT_EQ(container.get<int>("a"), 1);
    EXPECT_EQ(container.get<double>("b"), 2.0);
    EXPECT_EQ(container.get<bool>("d"), true);
  }

  TEST(InputSpecLegacyTest, OptionalLeftOut)
  {
    auto spec = all_of({
        parameter<int>("a"),
        parameter<double>("b"),
        parameter<std::string>("c", {.default_value = "default"}),
    });

    {
      InputParameterContainer container;
      std::string stream("a 1 b 2.0 // c 1");
      ValueParser parser(stream);
      spec.deprecated_parse(parser, container);
      EXPECT_EQ(container.get<int>("a"), 1);
      EXPECT_EQ(container.get<double>("b"), 2.0);
      EXPECT_EQ(container.get<std::string>("c"), "default");
    }
  }

  TEST(InputSpecLegacyTest, RequiredLeftOut)
  {
    auto spec = all_of({
        parameter<int>("a"),
        parameter<double>("b"),
        parameter<std::string>("c"),
    });

    {
      InputParameterContainer container;
      std::string stream("a 1 b 2.0");
      ValueParser parser(stream);
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(spec.deprecated_parse(parser, container), Core::Exception,
          "Required value 'c' not found in input line");
    }
  }

  TEST(InputSpecLegacyTest, EnumClassSelection)
  {
    enum class EnumClass
    {
      A,
      B,
      C,
    };

    auto spec = all_of({
        deprecated_selection<EnumClass>("enum",
            {
                {"A", EnumClass::A}, {"B", EnumClass::B},
                // Leave one out. Otherwise, we get an error to use a parameter<EnumClass> instead.
            }),
    });

    {
      InputParameterContainer container;
      std::string stream("enum A");
      ValueParser parser(stream);
      spec.deprecated_parse(parser, container);
      EXPECT_EQ(container.get<EnumClass>("enum"), EnumClass::A);
    }
  }

  TEST(InputSpecLegacyTest, MagicEnumParameter)
  {
    enum class EnumClass
    {
      A,
      B,
      C,
    };

    auto spec = parameter<EnumClass>(
        "enum", {.description = "An enum constant",
                    .enum_value_descriptions = {{EnumClass::A, "The option A"},
                        {EnumClass::B, "Other option"}, {EnumClass::C, "Other option"}}});

    {
      SCOPED_TRACE("Valid enum constant");
      InputParameterContainer container;
      std::string stream("enum A");
      ValueParser parser(stream);
      spec.deprecated_parse(parser, container);
      EXPECT_EQ(container.get<EnumClass>("enum"), EnumClass::A);
    }

    {
      SCOPED_TRACE("Invalid enum constant");
      InputParameterContainer container;
      std::string stream("enum XYZ");
      ValueParser parser(stream);
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(spec.deprecated_parse(parser, container), Core::Exception,
          "Could not parse value 'XYZ' as an enum constant of type 'EnumClass'");
    }

    {
      std::ostringstream out;
      ryml::Tree tree = init_yaml_tree_with_exceptions();
      ryml::NodeRef root = tree.rootref();
      YamlNodeRef yaml(root, "");
      spec.emit_metadata(yaml);
      out << tree;

      std::string expected = R"(name: enum
type: enum
description: "An enum constant"
required: true
choices:
  - name: A
    description: "The option A"
  - name: B
    description: "Other option"
  - name: C
    description: "Other option"
)";
      EXPECT_EQ(out.str(), expected);
    }
  }

  TEST(InputSpecLegacyTest, ParseSingleDefaultedEntryDat)
  {
    // This used to be a bug where a single default dat parameter was not accepted.
    auto spec = all_of({
        parameter<double>("a", {.default_value = 1.0}),
    });

    {
      InputParameterContainer container;
      std::string stream("");
      ValueParser parser(stream);
      spec.deprecated_parse(parser, container);
      EXPECT_EQ(container.get<double>("a"), 1.0);
    }
  }

  TEST(InputSpecLegacyTest, Vector)
  {
    auto spec = all_of({
        parameter<std::vector<std::vector<int>>>("a", {.size = {2, 2}}),
        parameter<std::vector<double>>("b", {.size = 3}),
    });

    {
      InputParameterContainer container;
      std::string stream("a 1 2 3 4 b 1.0 2.0 3.0");
      ValueParser parser(stream);
      spec.deprecated_parse(parser, container);
      const auto& a = container.get<std::vector<std::vector<int>>>("a");
      EXPECT_EQ(a.size(), 2);
      EXPECT_EQ(a[0].size(), 2);
      EXPECT_EQ(a[0][0], 1);
      EXPECT_EQ(a[0][1], 2);
      EXPECT_EQ(a[1].size(), 2);
      EXPECT_EQ(a[1][0], 3);
      EXPECT_EQ(a[1][1], 4);
      const auto& b = container.get<std::vector<double>>("b");
      EXPECT_EQ(b.size(), 3);
      EXPECT_EQ(b[0], 1.0);
      EXPECT_EQ(b[1], 2.0);
      EXPECT_EQ(b[2], 3.0);
    }

    {
      InputParameterContainer container;
      std::string stream("a 1 2 3 4 b 1.0 2.0 c");
      ValueParser parser(stream);
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(spec.deprecated_parse(parser, container), Core::Exception,
          "Could not parse 'c' as a double value");
    }
  }

  TEST(InputSpecLegacyTest, Optional)
  {
    auto spec = all_of({
        parameter<int>("size"),
        parameter<std::vector<std::optional<int>>>("vector_none",
            {
                .default_value = std::vector<std::optional<int>>{std::nullopt, 1},
                .size = from_parameter<int>("size"),
            }),
        parameter<std::optional<std::vector<int>>>(
            "none_vector", {.size = from_parameter<int>("size")}),
        parameter<std::optional<std::string>>("b", {.description = "b"}),
        parameter<std::optional<int>>("e"),
    });

    {
      SCOPED_TRACE("All values");
      InputParameterContainer container;
      std::string stream("size 3 vector_none 1 2 3 b none none_vector 1 2 3");
      ValueParser parser(stream);
      spec.deprecated_parse(parser, container);
      const auto& vector_none = container.get<std::vector<std::optional<int>>>("vector_none");
      EXPECT_EQ(vector_none.size(), 3);
      EXPECT_EQ(vector_none[0].has_value(), true);
      EXPECT_EQ(vector_none[0].value(), 1);
      EXPECT_EQ(vector_none[1].has_value(), true);
      EXPECT_EQ(vector_none[1].value(), 2);
      EXPECT_EQ(vector_none[2].has_value(), true);
      EXPECT_EQ(vector_none[2].value(), 3);

      EXPECT_TRUE(container.get<std::optional<std::vector<int>>>("none_vector").has_value());

      const auto& b = container.get<std::optional<std::string>>("b");
      EXPECT_EQ(b.has_value(), false);

      const auto& e = container.get<std::optional<int>>("e");
      EXPECT_EQ(e.has_value(), false);
    }

    {
      SCOPED_TRACE("None values");
      InputParameterContainer container;
      std::string stream("size 3 vector_none 1 none 3 b none e none none_vector none");
      ValueParser parser(stream);
      spec.deprecated_parse(parser, container);
      const auto& a = container.get<std::vector<std::optional<int>>>("vector_none");
      EXPECT_EQ(a.size(), 3);
      EXPECT_EQ(a[0].has_value(), true);
      EXPECT_EQ(a[0].value(), 1);
      EXPECT_EQ(a[1].has_value(), false);
      EXPECT_EQ(a[2].has_value(), true);
      EXPECT_EQ(a[2].value(), 3);

      const auto& b = container.get<std::optional<std::string>>("b");
      EXPECT_EQ(b.has_value(), false);

      const auto& e = container.get<std::optional<int>>("e");
      EXPECT_EQ(e.has_value(), false);
    }

    {
      SCOPED_TRACE("Defaults");
      InputParameterContainer container;
      std::string stream("size 3 b string e 42");
      ValueParser parser(stream);
      spec.deprecated_parse(parser, container);

      const auto& vector_none = container.get<std::vector<std::optional<int>>>("vector_none");
      EXPECT_EQ(vector_none.size(), 2);
      EXPECT_EQ(vector_none[0].has_value(), false);
      EXPECT_EQ(vector_none[1].has_value(), true);

      const auto& b = container.get<std::optional<std::string>>("b");
      EXPECT_EQ(b.has_value(), true);
      EXPECT_EQ(b.value(), "string");

      const auto& e = container.get<std::optional<int>>("e");
      EXPECT_EQ(e.has_value(), true);
      EXPECT_EQ(e.value(), 42);
    }
  }

  TEST(InputSpecLegacyTest, VectorWithParsedLength)
  {
    auto spec = all_of({
        parameter<int>("a"),
        parameter<std::vector<double>>("b", {.size = from_parameter<int>("a")}),
        parameter<std::string>("c"),
    });

    {
      InputParameterContainer container;
      std::string stream("a 3 b 1.0 2.0 3.0 c string");
      ValueParser parser(stream);
      spec.deprecated_parse(parser, container);
      EXPECT_EQ(container.get<int>("a"), 3);
      const auto& b = container.get<std::vector<double>>("b");
      EXPECT_EQ(b.size(), 3);
      EXPECT_EQ(b[0], 1.0);
      EXPECT_EQ(b[1], 2.0);
      EXPECT_EQ(b[2], 3.0);
      EXPECT_EQ(container.get<std::string>("c"), "string");
    }

    {
      InputParameterContainer container;
      std::string stream("a 3 b 1.0 2.0 c string");
      ValueParser parser(stream);
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(spec.deprecated_parse(parser, container), Core::Exception,
          "Could not parse 'c' as a double value");
    }
  }

  TEST(InputSpecLegacyTest, EntryWithCallback)
  {
    auto spec = all_of({
        parameter<int>("a"),
        parameter<double>("b"),
        parameter<std::string>("c",
            {
                .description = "A string",
                .default_value = "Not found",
                .on_parse_callback = [](InputParameterContainer& container)
                { container.add<int>("c_as_int", std::stoi(container.get<std::string>("c"))); },
            }),
        parameter<std::string>("s"),
    });

    {
      InputParameterContainer container;
      std::string stream("a 1 b 2.0 c 10 s hello");
      ValueParser parser(stream);
      spec.deprecated_parse(parser, container);
      EXPECT_EQ(container.get<int>("a"), 1);
      EXPECT_EQ(container.get<double>("b"), 2.0);
      EXPECT_EQ(container.get<std::string>("c"), "10");
      EXPECT_EQ(container.get<int>("c_as_int"), 10);
      EXPECT_EQ(container.get<std::string>("s"), "hello");
    }

    {
      InputParameterContainer container;
      std::string stream("a 1 b 2.0 c _ hello");
      ValueParser parser(stream);
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(
          spec.deprecated_parse(parser, container), std::invalid_argument, "stoi");
    }
  }

  TEST(InputSpecLegacyTest, Unparsed)
  {
    auto spec = all_of({
        parameter<int>("a"),
        parameter<int>("optional", {.default_value = 42}),
        parameter<double>("b"),
        parameter<std::string>("c"),
    });
    InputParameterContainer container;
    std::string stream("a 1 b 2.0 c string unparsed unparsed");
    ValueParser parser(stream);
    FOUR_C_EXPECT_THROW_WITH_MESSAGE(spec.deprecated_parse(parser, container), Core::Exception,
        "line still contains 'unparsed unparsed'");
  }


  TEST(InputSpecLegacyTest, Groups)
  {
    auto spec = all_of({
        parameter<int>("a"),
        group("group1",
            {
                parameter<double>("b"),
            }),
        group("group2",
            {
                parameter<double>("b", {.default_value = 3.0}),
                parameter<std::string>("c"),
            },
            {
                .required = false,
            }),
        group("group3",
            {
                parameter<std::string>("c", {.default_value = "default"}),
            },
            {
                .required = false,
            }),
        parameter<std::string>("c"),
    });

    {
      InputParameterContainer container;
      std::string stream("a 1 group1 b 2.0 c string");
      ValueParser parser(stream);
      spec.deprecated_parse(parser, container);
      const auto& const_container = container;
      EXPECT_EQ(const_container.get<int>("a"), 1);
      EXPECT_EQ(const_container.group("group1").get<double>("b"), 2.0);
      EXPECT_ANY_THROW([[maybe_unused]] const auto& c = const_container.group("group2"));
      // Group 3 only contains entries that have default values, so it implicitly has a default
      // value.
      EXPECT_EQ(const_container.group("group3").get<std::string>("c"), "default");
    }

    {
      InputParameterContainer container;
      std::string stream("a 1 group2 b 2.0 c string group1 b 4.0 c string");
      ValueParser parser(stream);
      spec.deprecated_parse(parser, container);
      const auto& const_container = container;
      EXPECT_EQ(const_container.get<int>("a"), 1);
      EXPECT_EQ(const_container.group("group2").get<double>("b"), 2.0);
      EXPECT_EQ(const_container.group("group1").get<double>("b"), 4.0);
      EXPECT_EQ(const_container.get<std::string>("c"), "string");
    }

    {
      InputParameterContainer container;
      std::string stream("a 1 group1 b 4.0 c string");
      ValueParser parser(stream);
      spec.deprecated_parse(parser, container);
      const auto& const_container = container;
      EXPECT_EQ(const_container.get<int>("a"), 1);
      EXPECT_ANY_THROW([[maybe_unused]] const auto& c = const_container.group("group2"));
      EXPECT_EQ(const_container.group("group1").get<double>("b"), 4.0);
      EXPECT_EQ(const_container.get<std::string>("c"), "string");
    }
  }

  TEST(InputSpecLegacyTest, NestedAllOf)
  {
    auto spec = all_of({
        parameter<int>("a"),
        all_of({
            all_of({
                parameter<double>("b"),
            }),
            // Not useful but might happen in practice, so ensure this can be handled.
            all_of({}),
        }),
    });

    {
      InputParameterContainer container;
      std::string stream("a 1 b 2.0");
      ValueParser parser(stream);
      spec.deprecated_parse(parser, container);
      const auto& const_container = container;
      EXPECT_EQ(const_container.get<int>("a"), 1);
      EXPECT_EQ(const_container.get<double>("b"), 2.0);
    }
  }

  TEST(InputSpecLegacyTest, OneOf)
  {
    auto spec = all_of({
        parameter<int>("a", {.default_value = 42}),
        one_of({
            parameter<double>("b"),
            group("group",
                {
                    parameter<std::string>("c"),
                    parameter<double>("d"),
                }),
        }),
    });

    {
      InputParameterContainer container;
      std::string stream("a 1 b 2");
      ValueParser parser(stream);
      spec.deprecated_parse(parser, container);
      EXPECT_EQ(container.get<int>("a"), 1);
      EXPECT_EQ(container.get<double>("b"), 2);
    }

    {
      InputParameterContainer container;
      std::string stream("group c string d 2.0");
      ValueParser parser(stream);
      spec.deprecated_parse(parser, container);
      EXPECT_EQ(container.get<int>("a"), 42);
      EXPECT_EQ(container.group("group").get<std::string>("c"), "string");
      EXPECT_EQ(container.group("group").get<double>("d"), 2);
    }

    {
      InputParameterContainer container;
      std::string stream("a 1 group c string d 2.0 b 3.0");
      ValueParser parser(stream);
      // More than one of the one_of entries is present. Refuse to parse any of them.
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(
          spec.deprecated_parse(parser, container), Core::Exception, "still contains 'b 3.0'");
    }

    {
      InputParameterContainer container;
      std::string stream("a 1 group c string");
      // Note: we start to parse the group, but the entries are not complete, so we backtrack.
      // The result is that the parts of the group remain unparsed.
      ValueParser parser(stream);
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(spec.deprecated_parse(parser, container), Core::Exception,
          "Required 'one_of' not found in input line");
    }

    {
      InputParameterContainer container;
      std::string stream("a 1");
      ValueParser parser(stream);
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(spec.deprecated_parse(parser, container), Core::Exception,
          "Required 'one_of' not found in input line");
    }
  }

  TEST(InputSpecLegacyTest, OneOfTopLevel)
  {
    auto spec = one_of(
        {
            all_of({
                parameter<int>("a"),
                parameter<double>("b"),
            }),
            all_of({
                parameter<std::string>("c"),
                parameter<double>("d"),
            }),
        },
        // Additionally store the index of the parsed group but map it to a different value.
        store_index_as<int>("index", /*reindex*/ {1, 10}));

    {
      InputParameterContainer container;
      std::string stream("a 1 b 2");
      ValueParser parser(stream);
      spec.deprecated_parse(parser, container);
      EXPECT_EQ(container.get<int>("a"), 1);
      EXPECT_EQ(container.get<double>("b"), 2);
      EXPECT_EQ(container.get<int>("index"), 1);
    }

    {
      InputParameterContainer container;
      std::string stream("c string d 2.0");
      ValueParser parser(stream);
      spec.deprecated_parse(parser, container);
      EXPECT_EQ(container.get<std::string>("c"), "string");
      EXPECT_EQ(container.get<double>("d"), 2);
      EXPECT_EQ(container.get<int>("index"), 10);
    }

    {
      InputParameterContainer container;
      std::string stream("a 1 b 2 c string d 2.0");
      ValueParser parser(stream);
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(
          spec.deprecated_parse(parser, container), Core::Exception, "Ambiguous input in one_of.");
    }

    {
      InputParameterContainer container;
      std::string stream("a 1 c string");
      ValueParser parser(stream);
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(spec.deprecated_parse(parser, container), Core::Exception,
          "None of the specs fit the input");
    }
  }

  TEST(InputSpecLegacyTest, DatToYaml)
  {
    auto spec = all_of({
        // awkward one_of where the first choice partially matches
        one_of({
            all_of({
                parameter<int>("a"),
                parameter<std::string>("b"),
            }),
            all_of({
                parameter<int>("a"),
                parameter<double>("d"),
            }),
        }),
        // group with all defaulted entries
        group("group",
            {
                parameter<double>("c", {.default_value = 1.0}),
            },
            {.required = false}),
        list("list",
            all_of({
                parameter<int>("l1"),
                parameter<double>("l2"),
            }),
            {.size = 2}),
        parameter<int>("i", {.default_value = 0}),
        parameter<std::vector<double>>("v", {.size = 3}),
    });

    std::string dat = "a 1 d 3.0 group c 1 i 42 v 1.0 2.0 3.0 list l1 1 l2 2.0 l1 3 l2 4.0";

    InputParameterContainer container;
    ValueParser parser(dat);
    spec.deprecated_parse(parser, container);

    {
      SCOPED_TRACE("Emit without default values");
      auto tree = init_yaml_tree_with_exceptions();
      YamlNodeRef yaml(tree.rootref(), "");
      spec.emit(yaml, container);

      std::ostringstream out;
      out << tree;
      std::string expected = R"(a: 1
d: 3
list:
  - l1: 1
    l2: 2
  - l1: 3
    l2: 4
i: 42
v: [1,2,3]
)";
      EXPECT_EQ(out.str(), expected);
    }

    {
      SCOPED_TRACE("Emit with defaulted values");
      auto tree = init_yaml_tree_with_exceptions();
      YamlNodeRef yaml(tree.rootref(), "");
      spec.emit(yaml, container, {.emit_defaulted_values = true});

      std::ostringstream out;
      out << tree;
      std::string expected = R"(a: 1
d: 3
group:
  c: 1
list:
  - l1: 1
    l2: 2
  - l1: 3
    l2: 4
i: 42
v: [1,2,3]
)";
      EXPECT_EQ(out.str(), expected);
    }
  }

}  // namespace