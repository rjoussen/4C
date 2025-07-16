// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_COMM_PACK_BUFFER_HPP
#define FOUR_C_COMM_PACK_BUFFER_HPP

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"

#include <cstring>
#include <type_traits>
#include <typeinfo>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Communication
{
  class PotentiallyUnusedBufferScope;

  /**
   * @brief A class to pack data into a buffer.
   *
   * This class allows to add various types of data to a buffer, which is essentially a vector of
   * characters. The buffer can be used to send data over MPI or serialize the data to disk.
   *
   * @note It is important that data added to the PackBuffer is extracted in the same order from
   * a corresponding UnpackBuffer.
   */
  class PackBuffer
  {
   public:
    PackBuffer() = default;

    std::vector<char>& operator()() { return buf_; }

    const std::vector<char>& operator()() const { return buf_; }

    /// Add a trivially copyable object, i.e., an object of a type that can be copied with memcpy.
    template <typename T>
      requires std::is_trivially_copyable_v<T>
    void add_to_pack(const T& stuff)
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      // Write the type into the buffer
      const std::size_t hash = typeid(T).hash_code();
      buf_.resize(buf_.size() + sizeof(hash));
      std::memcpy(buf_.data() + buf_.size() - sizeof(hash), &hash, sizeof(hash));
#endif

      buf_.resize(buf_.size() + sizeof(T));
      std::memcpy(buf_.data() + buf_.size() - sizeof(T), &stuff, sizeof(T));
    }

    /// Add an array of trivially copyable objects, i.e., objects of a type that can be copied with
    /// memcpy.
    template <typename T>
      requires std::is_trivially_copyable_v<T>
    void add_to_pack(const T* stuff, std::size_t stuff_size)
    {
      FOUR_C_ASSERT(stuff_size % sizeof(T) == 0, "Size of stuff must be a multiple of sizeof(T).");

#ifdef FOUR_C_ENABLE_ASSERTIONS
      // Write the type into the buffer
      const std::size_t hash = typeid(T).hash_code();
      buf_.resize(buf_.size() + sizeof(hash));
      std::memcpy(buf_.data() + buf_.size() - sizeof(hash), &hash, sizeof(hash));
#endif

      buf_.resize(buf_.size() + stuff_size);
      if (stuff_size > 0) std::memcpy(buf_.data() + buf_.size() - stuff_size, stuff, stuff_size);
    }

   private:
    //! The actual buffer containing the packed data.
    std::vector<char> buf_;

    friend class PotentiallyUnusedBufferScope;
  };


  /**
   * The counterpart of the PackBuffer class. This class is used to unpack data from a buffer that
   * was packed using the PackBuffer class. Internally tracks the position in the buffer every time
   * an object is extracted.
   */
  class UnpackBuffer
  {
   public:
    explicit UnpackBuffer(const std::vector<char>& data, std::size_t start_position = 0)
        : data_(data), position_(start_position)
    {
    }

    /*!
     * \brief Extract stuff from the UnpackBuffer.
     *
     * This method is a template for all POD types like int, char, double, ... . It extracts the
     * next object into @p stuff and assumes that the raw bytes actually came from an object of the
     * type T.
     */
    template <typename T>
      requires std::is_trivially_copyable_v<T>
    void extract_from_pack(T& stuff)
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      // Check that the type matches the type that was packed.
      std::size_t hash;
      std::memcpy(&hash, data_.data() + position_, sizeof(hash));
      position_ += sizeof(hash);
      FOUR_C_ASSERT(hash == typeid(T).hash_code(),
          "Type mismatch during unpacking. Tried to extract type {}", typeid(T).name());
#endif

      // We know that stuff is trivially copyable, so we can safely use memcpy. We cast to void* to
      // circumvent a GCC warning
      std::memcpy(static_cast<void*>(&stuff), data_.data() + position_, sizeof(T));
      position_ += sizeof(T);
    }

    /*!
     * \brief Extract multiple entries into @p stuff.
     *
     * Same as the other method but extracts @p stuff_size entries into the array @p stuff.
     */
    template <typename T>
      requires std::is_trivially_copyable_v<T>
    void extract_from_pack(T* stuff, const std::size_t stuff_size)
    {
      FOUR_C_ASSERT(stuff_size % sizeof(T) == 0, "Size of stuff must be a multiple of sizeof(T).");

#ifdef FOUR_C_ENABLE_ASSERTIONS
      // Check that the type matches the type that was packed.
      std::size_t hash;
      std::memcpy(&hash, data_.data() + position_, sizeof(hash));
      position_ += sizeof(hash);
      FOUR_C_ASSERT(hash == typeid(T).hash_code(),
          "Type mismatch during unpacking. Tried to extract type {}", typeid(T).name());
#endif

      FOUR_C_ASSERT(position_ + stuff_size <= data_.size(),
          "UnpackBuffer: Trying to extract more data than available in the buffer.");
      // We know that stuff is trivially copyable, so we can safely use memcpy. We cast to void* to
      // circumvent a GCC warning
      if (stuff_size > 0)
        std::memcpy(static_cast<void*>(stuff), data_.data() + position_, stuff_size);
      position_ += stuff_size;
    }

    /**
     * Get @p stuff but also leave it in the buffer for the next extraction.
     */
    template <typename T>
      requires std::is_trivially_copyable_v<T>
    void peek(T& stuff) const
    {
      std::size_t position = position_;

#ifdef FOUR_C_ENABLE_ASSERTIONS
      // Check that the type matches the type that was packed.
      std::size_t hash;
      std::memcpy(&hash, data_.data() + position_, sizeof(hash));
      position += sizeof(hash);
      FOUR_C_ASSERT(hash == typeid(T).hash_code(), "Type mismatch during unpacking.");
#endif
      FOUR_C_ASSERT(position < data_.size(),
          "UnpackBuffer: Trying to extract more data than available in the buffer.");
      std::memcpy(&stuff, data_.data() + position, sizeof(T));
    }

    /**
     * Returns true if this buffer contains no more data to be extracted.
     */
    [[nodiscard]] bool at_end() const { return position_ >= data_.size(); }

   private:
    const std::vector<char>& data_;
    std::vector<char>::size_type position_{0};

    friend class PotentiallyUnusedBufferScope;
  };

  /**
   * Using this class marks all data added/extracted within a scope as potentially unused. An object
   * of this class should symmetrically appear in a pack() and unpack() method. In the pack() method
   * it will remember how much data is potentially unused. In the unpack() method it will
   * discard all the unused data from the buffer at the end of the scope.
   */
  class [[nodiscard]] PotentiallyUnusedBufferScope
  {
   public:
    PotentiallyUnusedBufferScope(PackBuffer& pack_buffer)
        : pack_buffer_(&pack_buffer), offset_(pack_buffer_->buf_.size())
    {
      // Add an entry which stores the amount of potentially unused data.
      // Note: use raw memory operations to not write metadata
      pack_buffer_->buf_.insert(pack_buffer_->buf_.end(), sizeof(std::size_t), 0);
    }

    PotentiallyUnusedBufferScope(UnpackBuffer& unpack_buffer) : unpack_buffer_(&unpack_buffer)
    {
      // Read the size of the potentially unused data.
      std::memcpy(&offset_, &unpack_buffer.data_[unpack_buffer.position_], sizeof(offset_));
      unpack_buffer.position_ += sizeof(offset_);

      // This is the position where the potentially unused data ends in the buffer.
      offset_ += unpack_buffer.position_;
    }

    ~PotentiallyUnusedBufferScope()
    {
      if (pack_buffer_ != nullptr)
      {
        std::size_t size_of_unused_data = pack_buffer_->buf_.size() - offset_ - sizeof(std::size_t);
        // Overwrite the special entry with the size of the unused data
        std::memcpy(
            &pack_buffer_->buf_[offset_], &size_of_unused_data, sizeof(size_of_unused_data));
      }
      else
      {
        FOUR_C_ASSERT(
            unpack_buffer_->position_ <= offset_, "Potentially unused data exceeds buffer size.");
        // Now simply skip the unused data.
        unpack_buffer_->position_ = offset_;
      }
    }

    PotentiallyUnusedBufferScope(const PotentiallyUnusedBufferScope&) = delete;
    PotentiallyUnusedBufferScope& operator=(const PotentiallyUnusedBufferScope&) = delete;
    PotentiallyUnusedBufferScope(PotentiallyUnusedBufferScope&&) = delete;
    PotentiallyUnusedBufferScope& operator=(PotentiallyUnusedBufferScope&&) = delete;

   private:
    PackBuffer* pack_buffer_{nullptr};
    UnpackBuffer* unpack_buffer_{nullptr};
    std::size_t offset_;
  };


  /**
   * A concept to check whether a type T supports the pack(PackBuffer&) method.
   */
  template <typename T>
  concept HasPack = requires(T t, PackBuffer& buffer) { t.pack(buffer); };

  /**
   * A concept to check whether a type T supports the unpack(UnpackBuffer&) method.
   */
  template <typename T>
  concept HasUnpack = requires(T t, UnpackBuffer& buffer) { t.unpack(buffer); };

}  // namespace Core::Communication

FOUR_C_NAMESPACE_CLOSE

#endif
