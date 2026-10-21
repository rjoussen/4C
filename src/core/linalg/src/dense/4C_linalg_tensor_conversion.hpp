// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_TENSOR_CONVERSION_HPP
#define FOUR_C_LINALG_TENSOR_CONVERSION_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_linalg_tensor_internals.hpp"

#include <array>
#include <ranges>
#include <tuple>
#include <type_traits>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  /*!
   * @brief Create a 2-tensor given a fixed-size matrix.
   *
   * @param matrix
   * @return auto
   */
  template <typename T, unsigned int m, unsigned int n>
  auto make_tensor(const Matrix<m, n, T>& matrix)
  {
    using MatrixType = std::remove_cvref_t<decltype(matrix)>;
    Tensor<typename MatrixType::scalar_type, MatrixType::m(), MatrixType::n()> tensor;
    std::copy_n(matrix.data(), MatrixType::m() * MatrixType::n(), tensor.data());
    return tensor;
  }

  /*!
   * @brief Create a 2-tensor view onto a fixed-size matrix.
   *
   * @param matrix
   * @return auto
   */
  auto make_tensor_view(auto& matrix)
  {
    return make_tensor_view<std::remove_cvref_t<decltype(matrix)>::num_rows(),
        std::remove_cvref_t<decltype(matrix)>::num_cols()>(matrix.data());
  }

  /*!
   * @brief Creates a symmetric tensor from a stress-like Voigt notation matrix.
   */
  constexpr auto make_symmetric_tensor_from_stress_like_voigt_matrix(const auto& stress_like_voigt)
    requires(std::remove_cvref_t<decltype(stress_like_voigt)>::num_cols() == 1)
  {
    constexpr std::size_t voigt_size = std::remove_cvref_t<decltype(stress_like_voigt)>::num_rows();
    constexpr std::size_t symmetric_size = [](auto size) consteval
    {
      switch (size)
      {
        case 3:
          return 2;
        case 6:
          return 3;
        case 10:
          return 4;
        case 15:
          return 5;
        case 21:
          return 6;
        default:
          FOUR_C_THROW("Unsupported symmetric tensor size: " + std::to_string(size));
      }
    }(voigt_size);

    SymmetricTensor<typename std::remove_cvref_t<decltype(stress_like_voigt)>::scalar_type,
        symmetric_size, symmetric_size>
        symmetric_tensor;

    std::copy_n(stress_like_voigt.data(), voigt_size, symmetric_tensor.data());

    return symmetric_tensor;
  }


  /*!
   * @brief Creates a symmetric tensor from a strain-like Voigt notation matrix.
   */
  constexpr auto make_symmetric_tensor_from_strain_like_voigt_matrix(const auto& strain_like_voigt)
    requires(std::remove_cvref_t<decltype(strain_like_voigt)>::num_cols() == 1)
  {
    auto stress_like_voigt = strain_like_voigt;
    Voigt::Strains::to_stress_like(strain_like_voigt, stress_like_voigt);
    return make_symmetric_tensor_from_stress_like_voigt_matrix(stress_like_voigt);
  }

  /*!
   * @brief Create an arbitrarily-ranked tensor from a matrix.
   *
   * The number of elements of the tensor and the matrix must match. A common use case is to create
   * a rank-1 Tensor form a 3x1 Matrix.
   *
   * @tparam n Shape of the tensor
   */
  template <std::size_t... n>
  auto reinterpret_as_tensor(const auto& matrix)
    requires((n * ...) == std::remove_cvref_t<decltype(matrix)>::num_cols() *
                              std::remove_cvref_t<decltype(matrix)>::num_rows())
  {
    Tensor<typename std::remove_cvref_t<decltype(matrix)>::scalar_type, n...> tensor;
    std::copy_n(matrix.data(), (n * ...), tensor.data());
    return tensor;
  }

  /*!
   * @brief Create an arbitrarily-ranked tensor view from a matrix.
   *
   * The number of elements of the tensor view and the matrix must match. A common use case is to
   * create a rank-1 TensorView form a 3x1 Matrix-view.
   *
   * @tparam n
   */
  template <std::size_t... n>
  auto reinterpret_as_tensor_view(auto& matrix)
    requires((n * ...) == std::remove_cvref_t<decltype(matrix)>::num_cols() *
                              std::remove_cvref_t<decltype(matrix)>::num_rows())
  {
    return make_tensor_view<n...>(matrix.data());
  }

  /*!
   * @brief Creates a matrix that views a 2-tensor
   */
  auto make_matrix_view(auto& tensor)
    requires(std::remove_cvref_t<decltype(tensor)>::rank() == 2 &&
             !is_compressed_tensor<decltype(tensor)>)
  {
    using ValueType = std::remove_cvref_t<decltype(tensor)>::value_type;
    constexpr std::size_t n1 = std::remove_cvref_t<decltype(tensor)>::template extent<0>();
    constexpr std::size_t n2 = std::remove_cvref_t<decltype(tensor)>::template extent<1>();
    return Core::LinAlg::Matrix<n1, n2, ValueType>{tensor.data(), true};
  }

  /*!
   * @brief Creates a matrix of size n1xn2 that views a 1-tensor
   */
  template <std::size_t n1, std::size_t n2>
  auto make_matrix_view(auto& tensor)
    requires(std::remove_cvref_t<decltype(tensor)>::rank() == 1 &&
             n1 * n2 == std::remove_cvref_t<decltype(tensor)>::template extent<0>() &&
             !is_compressed_tensor<decltype(tensor)>)
  {
    using ValueType = std::remove_cvref_t<decltype(tensor)>::value_type;
    return Core::LinAlg::Matrix<n1, n2, ValueType>{tensor.data(), true};
  }

  /*!
   * @brief Creates a matrix view from a contiguous range of rank-1 tensors.
   *
   * The dimension of the tensors matches the number of rows and the number of tensors matches the
   * number of columns of the resulting matrix, e.g., given std::array<Tensor<T, dim>, n> the
   * resulting matrix will have shape (dim, n).
   */
  template <std::ranges::contiguous_range R>
  auto make_matrix_view(R& array_of_tensors)
    requires(Core::LinAlg::is_tensor<typename std::remove_cvref_t<R>::value_type> &&
             !Core::LinAlg::is_compressed_tensor<typename std::remove_cvref_t<R>::value_type> &&
             std::remove_cvref_t<R>::value_type::rank() == 1)
  {
    using T = typename std::remove_cvref_t<R>::value_type::value_type;
    constexpr std::size_t n1 = std::remove_cvref_t<R>::value_type::template extent<0>();
    constexpr std::size_t n2 = std::tuple_size_v<R>;
    return Core::LinAlg::Matrix<n1, n2, T>{array_of_tensors[0].data(), true};
  }

  /*!
   * @brief Creates a matrix from a 2-tensor
   */
  template <std::size_t n1, std::size_t n2, typename T>
  Core::LinAlg::Matrix<n1, n2, T> make_matrix(const Core::LinAlg::Tensor<T, n1, n2>& tensor)
  {
    return Core::LinAlg::Matrix<n1, n2, T>{tensor.data(), false};
  }

  /*!
   * @brief Creates a matrix of shape n1xn2 from a 1-tensor of shape n=n1*n2
   */
  template <std::size_t n1, std::size_t n2, typename T, std::size_t n>
    requires((n1 * n2) == n)
  Core::LinAlg::Matrix<n1, n2, T> make_matrix(const Core::LinAlg::Tensor<T, n>& tensor)
  {
    return Core::LinAlg::Matrix<n1, n2, T>{tensor.data(), false};
  }

  /*!
   * @brief Creates a stress-like Voigt notation matrix view from a symmetric tensor.
   *
   * If @p tensor is a symmetric 2-tensor of dimension 3, it will be converted to a 6x1 matrix.
   * If @p tensor is a symmetric 4-tensor of dimension 3, it will be converted to a 6x6 matrix.
   *
   */
  auto make_stress_like_voigt_view(auto& tensor)
    requires(is_symmetric_tensor<decltype(tensor)>)
  {
    using ValueType = std::remove_cvref_t<decltype(tensor)>::value_type;
    constexpr std::size_t rank = std::remove_cvref_t<decltype(tensor)>::rank();
    static_assert(rank == 2 || rank == 4,
        "Tensor must be a symmetric tensor of rank 2 or 4 for Voigt notation");

    if constexpr (rank == 2)
    {
      constexpr std::size_t compressed_size =
          std::remove_cvref_t<decltype(tensor)>::compressed_size;
      return Core::LinAlg::Matrix<compressed_size, 1, ValueType>(tensor.data(), true);
    }
    else if constexpr (rank == 4)
    {
      constexpr std::size_t size_left = std::remove_cvref_t<decltype(tensor)>::template extent<0>();
      constexpr std::size_t size_right =
          std::remove_cvref_t<decltype(tensor)>::template extent<2>();
      return Core::LinAlg::Matrix<size_left*(size_left + 1) / 2, size_right*(size_right + 1) / 2,
          ValueType>(tensor.data(), true);
    }
  }

  /*!
   * @brief Creates a 6x9 Voigt matrix from a 4th order tensor that has the minor symmetry in its
   * first index pair (C_ijkl = C_jikl).
   *
   * The first (symmetric) index pair is stored in stress-like 6-Voigt notation, the second
   * (non-symmetric) index pair as a 9-vector.
   */
  template <typename T>
  Core::LinAlg::Matrix<6, 9, T> make_6x9_voigt_matrix_from_tensor(
      const Core::LinAlg::Tensor<T, 3, 3, 3, 3>& tensor)
  {
    constexpr std::array<std::array<std::size_t, 2>, 6> row_index = {
        {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2}}};
    constexpr std::array<std::array<std::size_t, 2>, 9> col_index = {
        {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2}, {1, 0}, {2, 1}, {2, 0}}};

    Core::LinAlg::Matrix<6, 9, T> matrix_voigt;
    for (std::size_t r = 0; r < 6; ++r)
    {
      for (std::size_t c = 0; c < 9; ++c)
      {
        matrix_voigt(r, c) =
            0.5 * (tensor(row_index[r][0], row_index[r][1], col_index[c][0], col_index[c][1]) +
                      tensor(row_index[r][1], row_index[r][0], col_index[c][0], col_index[c][1]));
      }
    }
    return matrix_voigt;
  }

  /*!
   * @brief Creates a 9x6 Voigt matrix from a 4th order tensor that has the minor symmetry in its
   * second index pair (C_ijkl = C_ijlk).
   *
   * The first (non-symmetric) index pair is stored as a 9-vector, the second (symmetric) index
   * pair in stress-like 6-Voigt notation.
   */
  template <typename T>
  Core::LinAlg::Matrix<9, 6, T> make_9x6_voigt_matrix_from_tensor(
      const Core::LinAlg::Tensor<T, 3, 3, 3, 3>& tensor)
  {
    constexpr std::array<std::array<std::size_t, 2>, 9> row_index = {
        {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2}, {1, 0}, {2, 1}, {2, 0}}};
    constexpr std::array<std::array<std::size_t, 2>, 6> col_index = {
        {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2}}};

    Core::LinAlg::Matrix<9, 6, T> matrix_voigt;
    for (std::size_t r = 0; r < 9; ++r)
    {
      for (std::size_t c = 0; c < 6; ++c)
      {
        matrix_voigt(r, c) =
            0.5 * (tensor(row_index[r][0], row_index[r][1], col_index[c][0], col_index[c][1]) +
                      tensor(row_index[r][0], row_index[r][1], col_index[c][1], col_index[c][0]));
      }
    }
    return matrix_voigt;
  }

  /*!
   * @brief Creates a 6x6 stress-like Voigt matrix from a (minor-)symmetric 4th order tensor.
   *
   * The minor symmetries (IJ) and (KL) are encoded in the symmetric tensor, so its components map
   * directly to stress-like 6-Voigt notation. The result may still be a non-symmetric 6x6 matrix
   * if the tensor lacks major symmetry (C_IJKL != C_KLIJ).
   */
  template <typename T>
  Core::LinAlg::Matrix<6, 6, T> make_6x6_voigt_matrix_from_tensor(
      const Core::LinAlg::SymmetricTensor<T, 3, 3, 3, 3>& tensor)
  {
    // the symmetric tensor already stores its components in stress-like Voigt order
    return Core::LinAlg::Matrix<6, 6, T>(tensor.data());
  }

  /*!
   * @brief Creates a strain-like Voigt notation matrix from a symmetric tensor.
   *
   * If @p tensor is a symmetric 2-tensor of dimension 3, it will be converted to a 6x1 matrix in
   * strain-like Voigt notation.
   *
   */
  template <typename T, std::size_t size>
  Core::LinAlg::Matrix<size*(size + 1) / 2, 1, T> make_strain_like_voigt_matrix(
      const Core::LinAlg::SymmetricTensor<T, size, size>& tensor)
  {
    Core::LinAlg::Matrix<size*(size + 1) / 2, 1, T> matrix(tensor.data());
    Voigt::Stresses::to_strain_like(matrix, matrix);
    return matrix;
  }

  namespace Internal
  {
    template <typename T>
    struct NestedArrayHelper;

    template <typename T, std::size_t n1, std::size_t... n>
    struct NestedArrayHelper<Tensor<T, n1, n...>>
    {
      using type = std::array<typename NestedArrayHelper<Tensor<T, n...>>::type, n1>;
    };

    template <typename T, std::size_t n>
    struct NestedArrayHelper<Tensor<T, n>>
    {
      using type = std::array<T, n>;
    };
  }  // namespace Internal

  /*
   * Convert a nested std::array to a tensor.
   */
  template <typename T, std::size_t... n>
  Tensor<T, n...> make_tensor_from_nested_array(
      const typename Internal::NestedArrayHelper<Tensor<T, n...>>::type& nested_array)
  {
    constexpr std::array<std::size_t, (n * ...)> index_mapping =
        Internal::order_type_mapping<n...>();

    Tensor<T, n...> tensor;
    for (std::size_t i = 0; i < (n * ...); ++i)
    {
      auto view_to_first_element = [](const auto& nested_array)
      { return static_cast<const T*>(static_cast<const void*>(&nested_array)); };

      tensor.container()[index_mapping[i]] = *(view_to_first_element(nested_array) + i);
    }
    return tensor;
  }

  /*
   * Convert a Tensor to a nested std::array.
   */
  template <typename T, std::size_t... n>
  typename Internal::NestedArrayHelper<Tensor<T, n...>>::type make_nested_array_from_tensor(
      const Tensor<T, n...>& tensor)
  {
    constexpr std::array<std::size_t, (n * ...)> index_mapping =
        Internal::order_type_mapping<n...>();

    typename Internal::NestedArrayHelper<Tensor<T, n...>>::type nested_array;
    for (std::size_t i = 0; i < (n * ...); ++i)
    {
      auto view_to_first_element = [](auto& nested_array)
      { return static_cast<T*>(static_cast<void*>(&nested_array)); };

      *(view_to_first_element(nested_array) + i) = tensor.container()[index_mapping[i]];
    }
    return nested_array;
  }
}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif