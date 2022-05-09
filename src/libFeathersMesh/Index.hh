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
  constexpr Index& operator+=(size_t value) noexcept {
    Value_ += value;
    return *this;
  }

  /// @brief Index subtraction operator.
  constexpr Index& operator-=(size_t value) noexcept {
    Value_ -= value;
    return *this;
  }

}; // class Index<...>

/// @brief Index difference operator.
template<class Tag>
constexpr ptrdiff_t operator-(Index<Tag> i1, Index<Tag> i2) noexcept {
  return size_t(i1) - size_t(i2);
}

/// @brief Index equality operator.
/// @{
template<class Tag>
constexpr bool operator==(Index<Tag> i1, size_t i2) noexcept {
  return size_t(i1) == i2;
}
template<class Tag>
constexpr bool operator==(Index<Tag> i1, Index<Tag> i2) noexcept {
  return size_t(i1) == size_t(i2);
}
/// @}

/// @brief Inequality operator.
/// @{
template<class Tag>
constexpr bool operator!=(Index<Tag> i1, size_t i2) noexcept {
  return size_t(i1) != i2;
}
template<class Tag>
constexpr bool operator!=(Index<Tag> i1, Index<Tag> i2) noexcept {
  return size_t(i1) != size_t(i2);
}
/// @}

/// @brief Index less than operator.
/// @{
template<class Tag>
constexpr bool operator<(Index<Tag> i1, size_t i2) noexcept {
  return size_t(i1) < i2;
}
template<class Tag>
constexpr bool operator<(Index<Tag> i1, Index<Tag> i2) noexcept {
  return size_t(i1) < size_t(i2);
}
/// @}

/// @brief Index less than or equal operator.
/// @{
template<class Tag>
constexpr bool operator<=(Index<Tag> i1, size_t i2) noexcept {
  return size_t(i1) <= i2;
}
template<class Tag>
constexpr bool operator<=(Index<Tag> i1, Index<Tag> i2) noexcept {
  return size_t(i1) <= size_t(i2);
}
/// @}

/// @brief Index greater than operator.
/// @{
template<class Tag>
constexpr bool operator>(Index<Tag> i1, size_t i2) noexcept {
  return size_t(i1) > i2;
}
template<class Tag>
constexpr bool operator>(Index<Tag> i1, Index<Tag> i2) noexcept {
  return size_t(i1) > size_t(i2);
}
/// @}

/// @brief Index greater than or equal operator.
/// @{
template<class Tag>
constexpr bool operator>=(Index<Tag> i1, size_t i2) noexcept {
  return size_t(i1) >= i2;
}
template<class Tag>
constexpr bool operator>=(Index<Tag> i1, Index<Tag> i2) noexcept {
  return size_t(i1) >= size_t(i2);
}
/// @}

/// @brief Index addition operator.
/// @{
template<class Tag>
constexpr Index<Tag> operator+(Index<Tag> i1, size_t i2) noexcept {
  return i1 += i2;
}
template<class Tag>
constexpr Index<Tag> operator+(size_t i1, Index<Tag> i2) noexcept {
  return Index<Tag>(i1) += size_t(i1);
}
/// @}

/// @brief Index subtraction operator.
/// @{
template<class Tag>
constexpr Index<Tag> operator-(Index<Tag> i1, size_t i2) noexcept {
  return i1 -= i2;
}
template<class Tag>
constexpr Index<Tag> operator-(size_t i1, Index<Tag> i2) noexcept {
  return Index<Tag>(i1) -= size_t(i1);
}
/// @}

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
