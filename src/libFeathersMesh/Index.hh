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

template<class Value, class Tag>
class Index {
private:
  Value Value_{};

public:

  constexpr Index() noexcept = default;

  constexpr explicit Index(Value value) noexcept : Value_(value) {}

  /// @brief Cast to the underlying value operator.
  constexpr explicit operator Value() const noexcept {
    return Value_;
  }

  /// @brief Cast to the other index operator.
  template<class OtherValue, class OtherTag>
  constexpr explicit operator Index<OtherValue, OtherTag>() const noexcept {
    return Index<OtherValue, OtherTag>(static_cast<OtherValue>(Value_));
  }

  /// @brief Difference operator.
  constexpr auto operator-(Index other) const noexcept {
    return Value_ - other.Value_;
  }

  /// @brief Equality operator.
  /// @{
  constexpr bool operator==(Value value) const noexcept {
    return Value_ == value;
  }
  constexpr bool operator==(Index other) const noexcept {
    return Value_ == other.Value_;
  }
  /// @}

  /// @brief Inequality operator.
  /// @{
  constexpr bool operator!=(Value value) const noexcept {
    return Value_ != value;
  }
  constexpr bool operator!=(Index other) const noexcept {
    return Value_ != other.Value_;
  }
  /// @}

  /// @brief Less than operator.
  /// @{
  constexpr bool operator<(Value value) const noexcept {
    return Value_ < value;
  }
  constexpr bool operator<(Index other) const noexcept {
    return Value_ < other.Value_;
  }
  /// @}

  /// @brief Less than or equal operator.
  /// @{
  constexpr bool operator<=(Value value) const noexcept {
    return Value_ <= value;
  }
  constexpr bool operator<=(Index other) const noexcept {
    return Value_ <= other.Value_;
  }
  /// @}

  /// @brief Greater than operator.
  /// @{
  constexpr bool operator>(Value value) const noexcept {
    return Value_ > value;
  }
  constexpr bool operator>(Index other) const noexcept {
    return Value_ > other.Value_;
  }
  /// @}

  /// @brief Greater than or equal operator.
  /// @{
  constexpr bool operator>=(Value value) const noexcept {
    return Value_ >= value;
  }
  constexpr bool operator>=(Index other) const noexcept {
    return Value_ >= other.Value_;
  }
  /// @}

  /// @brief Increment operator.
  /// @{
  constexpr Index& operator++() noexcept {
    ++Value_;
    return *this;
  }
  constexpr Index const operator++(int) noexcept {
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
  constexpr Index const operator--(int) noexcept {
    Index index(*this);
    return --(*this), index;
  }
  /// @}

  /// @brief Addition operator.
  /// @{
  constexpr Index& operator+=(uint_t offset) noexcept {
    Value_ += offset;
    return *this;
  }
  constexpr Index operator+(uint_t offset) const noexcept {
    return Index(*this) += offset;
  }
  /// @}

  /// @brief Subtraction operator.
  /// @{
  constexpr Index& operator-=(int_t offset) noexcept {
    Value_ -= offset;
    return *this;
  }
  constexpr Index operator-(int_t offset) const noexcept {
    return Index(*this) -= offset;
  }
  /// @}

}; // class Index<...>

/// @todo
template<class Value, class Index>
class IndexedVector : public std::vector<Value> {
public:

  IndexedVector() = default;
  IndexedVector(std::initializer_list<Value> list) : std::vector<Value>(list) {}

  Value& operator[](Index index) noexcept {
    return std::vector<Value>::operator[]((uint_t)index);
  }
  Value const& operator[](Index index) const noexcept {
    return std::vector<Value>::operator[]((uint_t)index);
  }

};

} // namespace feathers
