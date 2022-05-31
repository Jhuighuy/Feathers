/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person
/// obtaining a copy of this software and associated documentation
/// files (the "Software"), to deal in the Software without
/// restriction, including without limitation the rights  to use,
/// copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the
/// Software is furnished to do so, subject to the following
/// conditions:
///
/// The above copyright notice and this permission notice shall be
/// included in all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
/// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
/// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
/// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
/// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
/// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
/// OTHER DEALINGS IN THE SOFTWARE.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///

#pragma once

#include <concepts>

#include <stormBase.hxx>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Strict index.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Tag>
class Index {
private:
  size_t value_{};

public:

  /// @brief Default constructor.
  constexpr Index() noexcept = default;

  /// @brief Construct index with a value.
  constexpr explicit Index(size_t value) noexcept : value_(value) {}

  /// @brief Cast to the underlying value operator.
  template<std::integral Integer = size_t>
  constexpr explicit operator Integer() const noexcept {
    return value_;
  }

  /// @brief Cast to the other index operator.
  template<class OtherTag>
  constexpr explicit operator Index<OtherTag>() const noexcept {
    return Index<OtherTag>(value_);
  }

  /// @brief Index comparison operator.
  /// @{
#if 0 /// @todo Why the spaceship fails in ranges?
  constexpr auto operator<=>(Index other) const noexcept {
    return value_ <=> other.value_;
  }
  constexpr auto operator<=>(size_t value) const noexcept {
    return value_ <=> value;
  }
#else
#define Operator_(OP) \
  constexpr bool operator OP(Index other) const noexcept { \
    return value_ OP other.value_; \
  } \
  friend constexpr bool operator OP(Index first, size_t second) noexcept { \
    return first.value_ OP second; \
  } \
  friend constexpr bool operator OP(size_t first, Index second) noexcept { \
    return first OP second.value_; \
  }
  Operator_(==) Operator_(!=)
  Operator_(<) Operator_(<=) Operator_(>) Operator_(>=)
#undef Operator_
#endif
  /// @}

  /// @brief Increment operator.
  /// @{
  constexpr Index& operator++() noexcept {
    ++value_;
    return *this;
  }
  constexpr Index operator++(int) noexcept {
    Index index(*this);
    return ++(*this), index;
  }
  /// @}

  /// @brief Decrement operator.
  /// @{
  constexpr Index& operator--() noexcept {
    --value_;
    return *this;
  }
  constexpr Index operator--(int) noexcept {
    Index index(*this);
    return --(*this), index;
  }
  /// @}

  /// @brief Index addition operator.
  /// @{
  constexpr Index& operator+=(size_t value) noexcept {
    value_ += value;
    return *this;
  }
  friend constexpr Index operator+(Index first, size_t second) noexcept {
    return Index(first.value_ + second);
  }
  friend constexpr Index operator+(size_t first, Index second) noexcept {
    return Index(first + second.value_);
  }
  /// @}

  /// @brief Index subtraction operator.
  /// @{
  constexpr Index& operator-=(size_t value) noexcept {
    value_ -= value;
    return *this;
  }
  friend constexpr Index operator-(Index first, size_t second) noexcept {
    return Index(first.value_ - second);
  }
  friend constexpr Index operator-(size_t first, Index second) noexcept {
    return Index(first - second.value_);
  }
  /// @}

  /// @brief Index difference operator.
  constexpr ptrdiff_t operator-(Index other) const noexcept {
    return value_ - other.value_;
  }

  /// @brief Read @p index from the input @p stream.
  friend std::istream& operator>>(std::istream& stream, Index& index) {
    return stream >> index.value_;
  }

  /// @brief Write @p index to the output @p stream.
  friend std::ostream& operator<<(std::ostream& stream, Index index) {
    return stream << index.value_;
  }

}; // class Index<...>

/// @todo Move me to the better place!
template<class Value, class Index>
class Vector : public std::vector<Value> {
public:

  Vector() = default;
  Vector(std::initializer_list<Value> list) : std::vector<Value>(list) {}

  Value& operator[](Index index) noexcept {
    return std::vector<Value>::operator[]((size_t)index);
  }
  Value const& operator[](Index index) const noexcept {
    return std::vector<Value>::operator[]((size_t)index);
  }

};

} // namespace Storm
