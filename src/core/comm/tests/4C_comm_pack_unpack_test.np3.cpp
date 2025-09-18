// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_comm_pack_helpers.hpp"

using namespace FourC;

namespace
{
  void fill_data(int& data) { data = 42; }

  void fill_data(double& data) { data = 3.14; }

  void fill_data(bool& data) { data = true; }

  void fill_data(char& data) { data = 'a'; }

  void fill_data(std::string& data) { data = "Hello, World!"; }

  template <typename T>
  void fill_data(std::vector<T>& data)
  {
    data.resize(2);
    fill_data(data[0]);
    fill_data(data[1]);
  }

  void fill_data(std::optional<int>& data) { data.reset(); }

  template <typename T>
  void fill_data(std::optional<T>& data)
  {
    T d;
    fill_data(d);
    data.emplace(d);
  }

  template <typename K, typename V>
  void fill_data(std::map<K, V>& data)
  {
    K k;
    fill_data(k);
    fill_data(data[k]);
  }

  template <typename First, typename Second>
  void fill_data(std::pair<First, Second>& data)
  {
    fill_data(data.first);
    fill_data(data.second);
  }

  struct TriviallyCopyable
  {
    int a;
    double b;
    char c;

    auto operator<=>(const TriviallyCopyable&) const = default;
  };
  static_assert(std::is_trivially_copyable_v<TriviallyCopyable>);

  void fill_data(TriviallyCopyable& data)
  {
    fill_data(data.a);
    fill_data(data.b);
    fill_data(data.c);
  }

  struct CustomPack
  {
    std::string data;

    void pack(Core::Communication::PackBuffer& buffer) const
    {
      Core::Communication::add_to_pack(buffer, data);
    }
    void unpack(Core::Communication::UnpackBuffer& buffer)
    {
      Core::Communication::extract_from_pack(buffer, data);
    }
    auto operator<=>(const CustomPack&) const = default;
  };

  void fill_data(CustomPack& data) { fill_data(data.data); }


  struct TriviallyCopyableAndCustomPack
  {
    int a;

    void pack(Core::Communication::PackBuffer& buffer) const
    {
      int a_plus_one = a + 1;
      Core::Communication::add_to_pack(buffer, a_plus_one);
    }

    void unpack(Core::Communication::UnpackBuffer& buffer)
    {
      int a_plus_one;
      Core::Communication::extract_from_pack(buffer, a_plus_one);
      a = a_plus_one - 1;
    }

    auto operator<=>(const TriviallyCopyableAndCustomPack&) const = default;
  };
  static_assert(std::is_trivially_copyable_v<TriviallyCopyableAndCustomPack>);

  void fill_data(TriviallyCopyableAndCustomPack& data) { fill_data(data.a); }


  /**
   * A test fixture to test that basic types can round-trip through the pack/unpack mechanism.
   */
  template <typename T>
  class PackUnpackRoundtripTest : public ::testing::Test
  {
  };

  using TestTypes = ::testing::Types<int, double, char, std::string, std::vector<int>,
      std::vector<std::vector<std::vector<int>>>, std::map<std::string, bool>, std::optional<int>,
      std::optional<std::vector<double>>,
      // custom types
      TriviallyCopyable, CustomPack, TriviallyCopyableAndCustomPack>;

  TYPED_TEST_SUITE(PackUnpackRoundtripTest, TestTypes);

  TYPED_TEST(PackUnpackRoundtripTest, Empty)
  {
    TypeParam data{};
    Core::Communication::PackBuffer pack_buffer;
    Core::Communication::add_to_pack(pack_buffer, data);

    TypeParam unpacked_data;
    Core::Communication::UnpackBuffer unpack_buffer(pack_buffer());
    Core::Communication::extract_from_pack(unpack_buffer, unpacked_data);

    EXPECT_EQ(data, unpacked_data);
  }

  TYPED_TEST(PackUnpackRoundtripTest, NonEmpty)
  {
    TypeParam data;
    fill_data(data);
    Core::Communication::PackBuffer pack_buffer;
    Core::Communication::add_to_pack(pack_buffer, data);

    TypeParam unpacked_data;
    Core::Communication::UnpackBuffer unpack_buffer(pack_buffer());
    Core::Communication::extract_from_pack(unpack_buffer, unpacked_data);

    EXPECT_EQ(data, unpacked_data);
  }

}  // namespace