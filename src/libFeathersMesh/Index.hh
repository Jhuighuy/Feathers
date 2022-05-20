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

namespace feathers {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Strict index.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Tag>
class Index {
private:
  size_t Value_{};

public:

  /// @brief Default constructor.
  constexpr Index() noexcept = default;

  /// @brief Construct index with a value.
  constexpr explicit Index(size_t value) noexcept : Value_(value) {}

  /// @brief Cast to the underlying value operator.
  constexpr explicit operator size_t() const noexcept {
    return Value_;
  }

  /// @brief Cast to the other index operator.
  template<class OtherTag>
  constexpr explicit operator Index<OtherTag>() const noexcept {
    return Index<OtherTag>(Value_);
  }

  /// @brief Index comparison operator.
  /// @{
#if 0 /// @todo Why the spaceship fails in ranges?
  constexpr auto operator<=>(Index other) const noexcept {
    return Value_ <=> other.Value_;
  }
  constexpr auto operator<=>(size_t value) const noexcept {
    return Value_ <=> value;
  }
#else
#define Operator_(OP) \
  constexpr bool operator OP(Index other) const noexcept { \
    return Value_ OP other.Value_; \
  } \
  friend constexpr bool operator OP(Index first, size_t second) noexcept { \
    return first.Value_ OP second; \
  } \
  friend constexpr bool operator OP(size_t first, Index second) noexcept { \
    return first OP second.Value_; \
  }
  Operator_(==) Operator_(!=)
  Operator_(<) Operator_(<=) Operator_(>) Operator_(>=)
#undef Operator_
#endif
  /// @}

  /// @brief Increment operator.
  /// @{
  constexpr Index& operator++() noexcept {
    ++Value_;
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
    --Value_;
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
    Value_ += value;
    return *this;
  }
  friend constexpr Index operator+(Index first, size_t second) noexcept {
    return Index(first.Value_ + second);
  }
  friend constexpr Index operator+(size_t first, Index second) noexcept {
    return Index(first + second.Value_);
  }
  /// @}

  /// @brief Index subtraction operator.
  /// @{
  constexpr Index& operator-=(size_t value) noexcept {
    Value_ -= value;
    return *this;
  }
  friend constexpr Index operator-(Index first, size_t second) noexcept {
    return Index(first.Value_ - second);
  }
  friend constexpr Index operator-(size_t first, Index second) noexcept {
    return Index(first - second.Value_);
  }
  /// @}

  /// @brief Index difference operator.
  constexpr ptrdiff_t operator-(Index other) const noexcept {
    return Value_ - other.Value_;
  }

}; // class Index<...>

/// @todo
template<class Value, class Index>
class IndexedVector : public std::vector<Value> {
public:

  IndexedVector() = default;
  IndexedVector(std::initializer_list<Value> list) : std::vector<Value>(list) {}

  Value& operator[](Index index) noexcept {
    return std::vector<Value>::operator[]((size_t)index);
  }
  Value const& operator[](Index index) const noexcept {
    return std::vector<Value>::operator[]((size_t)index);
  }

};

} // namespace feathers
