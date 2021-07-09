/**
 *    ______             __     __  _____ _____
 *   / __/ /____ _____  / /__  /  |/  / // / _ \
 *  _\ \/  '_/ // / _ \/  '_/ / /|_/ / _  / // /
 * /___/_/\_\\_,_/_//_/_/\_\ /_/  /_/_//_/____/
 *
 * Copyright (c) 2019 Oleg Butakov
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * iy the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included iy all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. iy NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER iy AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR iy CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS iy THE
 * SOFTWARE.
 */

#pragma once
#ifndef GEOM_MATRIX_HH_
#define GEOM_MATRIX_HH_

#include <SkunkBase.hh>

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

/** Fixed size matrix. */
template<typename scalar_t, uint_t nx, uint_t ny = nx>
class TMatrix;

/** Fixed size column vector. */
template<typename scalar_t, uint_t n>
using TVector = TMatrix<scalar_t, n, 1>;
/** Fixed size row vector. */
template<typename scalar_t, uint_t n>
using TRowVector = TMatrix<scalar_t, 1, n>;

/** Real 2D column vector. */
using vec2_t = TVector<real_t, 2>;
/** Real 3D column vector. */
using vec3_t = TVector<real_t, 3>;
/** Real 4D column vector. */
using vec4_t = TVector<real_t, 4>;
/** Real 5D column vector. */
using vec5_t = TVector<real_t, 5>;
/** Real 6D column vector. */
using vec6_t = TVector<real_t, 6>;
/** Real ND column vector. */
template<uint_t n>
using vecn_t = TVector<real_t, n>;

/** Integral 2D column vector. */
using ivec2_t = TVector<int_t, 2>;
/** Integral 3D column vector. */
using ivec3_t = TVector<int_t, 3>;
/** Integral 4D column vector. */
using ivec4_t = TVector<int_t, 4>;
/** Integral 5D column vector. */
using ivec5_t = TVector<int_t, 5>;
/** Integral 6D column vector. */
using ivec6_t = TVector<int_t, 6>;
/** Integral ND column vector. */
template<uint_t n>
using ivecn_t = TVector<int_t, n>;

/** Unsigned integral 2D column vector. */
using uvec2_t = TVector<uint_t, 2>;
/** Unsigned integral 3D column vector. */
using uvec3_t = TVector<uint_t, 3>;
/** Unsigned integral 4D column vector. */
using uvec4_t = TVector<uint_t, 4>;
/** Unsigned integral 5D column vector. */
using uvec5_t = TVector<uint_t, 5>;
/** Unsigned integral 6D column vector. */
using uvec6_t = TVector<uint_t, 6>;
/** Unsigned integral ND column vector. */
template<uint_t n>
using uvecn_t = TVector<uint_t, n>;

/** Real 2D column matrix. */
using mat2_t = TMatrix<real_t, 2>;
/** Real 3D matrix. */
using mat3_t = TMatrix<real_t, 3>;
/** Real 4D matrix. */
using mat4_t = TMatrix<real_t, 4>;
/** Real 5D matrix. */
using mat5_t = TMatrix<real_t, 5>;
/** Real 6D matrix. */
using mat6_t = TMatrix<real_t, 6>;
/** Real ND matrix. */
template<uint_t n>
using matn_t = TMatrix<real_t, n>;

/** Integral 2D column matrix. */
using imat2_t = TMatrix<int_t, 2>;
/** Integral 3D matrix. */
using imat3_t = TMatrix<int_t, 3>;
/** Integral 4D matrix. */
using imat4_t = TMatrix<int_t, 4>;
/** Integral 5D matrix. */
using imat5_t = TMatrix<int_t, 5>;
/** Integral 6D matrix. */
using imat6_t = TMatrix<int_t, 6>;
/** Integral ND matrix. */
template<uint_t n>
using imatn_t = TMatrix<int_t, n>;

/** Unsigned integral 2D column matrix. */
using umat2_t = TMatrix<uint_t, 2>;
/** Unsigned integral 3D matrix. */
using umat3_t = TMatrix<uint_t, 3>;
/** Unsigned integral 4D matrix. */
using umat4_t = TMatrix<uint_t, 4>;
/** Unsigned integral 5D matrix. */
using umat5_t = TMatrix<uint_t, 5>;
/** Unsigned integral 6D matrix. */
using umat6_t = TMatrix<uint_t, 6>;
/** Unsigned integral ND matrix. */
template<uint_t n>
using umatn_t = TMatrix<uint_t, n>;

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

/** Fixed size row major matrix. */
template<typename scalar_t, uint_t nx, uint_t ny>
class TMatrixBase {
public:
    /** Rows data. */
    std::array<std::array<scalar_t, ny>, nx> rows_data{};
};  // class TMatrixBase

template<typename scalar_t>
class TMatrixBase<scalar_t, 3, 3> {
public:
    std::array<std::array<scalar_t, 3>, 3> rows_data{};

    template<typename TMatrix>
    static scalar_t det(const TMatrix& m) {
        return m(0, 0)*(m(1, 1)*m(2, 2) - m(1, 2)*m(2, 1)) -
               m(0, 1)*(m(1, 0)*m(2, 2) - m(1, 2)*m(2, 0)) + m(0, 2)*(m(1, 0)*m(2, 1) - m(1, 1)*m(2, 0));
    }
    template<typename TMatrix>
    static TMatrix inv(const TMatrix& m) {
        TMatrix minv;
        scalar_t invdet{ 1.0/det(m) };
        minv(0, 0) = (m(1, 1)*m(2, 2) - m(2, 1)*m(1, 2))*invdet;
        minv(0, 1) = (m(0, 2)*m(2, 1) - m(0, 1)*m(2, 2))*invdet;
        minv(0, 2) = (m(0, 1)*m(1, 2) - m(0, 2)*m(1, 1))*invdet;
        minv(1, 0) = (m(1, 2)*m(2, 0) - m(1, 0)*m(2, 2))*invdet;
        minv(1, 1) = (m(0, 0)*m(2, 2) - m(0, 2)*m(2, 0))*invdet;
        minv(1, 2) = (m(1, 0)*m(0, 2) - m(0, 0)*m(1, 2))*invdet;
        minv(2, 0) = (m(1, 0)*m(2, 1) - m(2, 0)*m(1, 1))*invdet;
        minv(2, 1) = (m(2, 0)*m(0, 1) - m(0, 0)*m(2, 1))*invdet;
        minv(2, 2) = (m(0, 0)*m(1, 1) - m(1, 0)*m(0, 1))*invdet;
        return minv;
    }
};  // class TMatrixBase<1, 1>

/**************************************************************************/
/**************************************************************************/

/** Base 1D column vector. */
template<typename scalar_t>
class TMatrixBase<scalar_t, 1, 1> {
public:
    union {
        std::array<std::array<scalar_t, 1>, 1> rows_data{};
        struct { scalar_t x; };
    };
    /** Implicit cast to scalar.
     ** It is used to cast inner products into scalars. */
    operator scalar_t() const {
        return x;
    }
};  // class TMatrixBase<1, 1>

/** Base 2D column vector. */
template<typename scalar_t>
class TMatrixBase<scalar_t, 2, 1> {
public:
    union {
        std::array<std::array<scalar_t, 1>, 2> rows_data{};
        struct { scalar_t x, y; };
    };
#if 1
    /** Length of a 2D vector. */
    constexpr scalar_t len() const {
        return std::hypot(x, y);
    }
    /** Dot product of two 2D vectors. */
    template<typename TMatrix>
    constexpr scalar_t dot(const TMatrix& other) const {
        return x*other.x + y*other.y;
    }
    /** Determinant of the two 2D vectors. */
    template<typename TMatrix>
    constexpr scalar_t det(const TMatrix& other) const {
        return x*other.y - other.x*y;
    }
#endif
};  // class TMatrixBase<2, 1>

/** Base 3D column vector. */
template<typename scalar_t>
class TMatrixBase<scalar_t, 3, 1> {
public:
    union {
        /** Rows data. */
        std::array<std::array<scalar_t, 1>, 3> rows_data{};
        struct { scalar_t x, y, z; };
    };
#if 1
public:
    /** Length of a 3D vector. */
    constexpr scalar_t len() const {
#if SKUNK_HAS_CPP17
        return std::hypot(x, y, z);
#else
        return std::sqrt(x*x + y*y + z*z);
#endif
    }
    template<typename TMatrix>
    static scalar_t len(const TMatrix& other) {
        return other.len();
    }
    template<typename TMatrix>
    static TMatrix cross(const TMatrix& other, const TMatrix& other1) {
        return other.cross(other1);
    }
    /** Dot product of two 3D vectors.
     ** @todo Remove me! */
    template<typename TMatrix>
    constexpr scalar_t dot(const TMatrix& other) const {
        return x*other.x + y*other.y + z*other.z;
    }
    /** Cross product of two 3D vectors. */
    template<typename TMatrix>
    constexpr TMatrix cross(const TMatrix& other) const {
        return {y*other.z - other.y*z, z*other.x - x*other.z, x*other.y - other.x*y};
    }
    /** Determinant of the three 3D vectors. */
    template<typename TMatrix>
    constexpr scalar_t det(const TMatrix& other1, const TMatrix& other2) const {
        return dot(other1.cross(other2));
    }
#endif
};  // class TMatrixBase<3, 1>

/** Base 4D column vector. */
template<typename scalar_t>
class TMatrixBase<scalar_t, 4, 1> {
public:
    union {
        /** Rows data. */
        std::array<std::array<scalar_t, 1>, 4> rows_data{};
        struct { scalar_t x, y, z, u; };
    };
};  // class TMatrixBase<4, 1>

/** Base 5D column vector. */
template<typename scalar_t>
class TMatrixBase<scalar_t, 5, 1> {
public:
    union {
        /** Rows data. */
        std::array<std::array<scalar_t, 1>, 5> rows_data{};
        struct { scalar_t x, y, z, u, v; };
    };
};  // class TMatrixBase<5, 1>

/** Base 6D column vector. */
template<typename scalar_t>
class TMatrixBase<scalar_t, 6, 1> {
public:
    union {
        /** Rows data. */
        std::array<std::array<scalar_t, 1>, 6> rows_data{};
        struct { scalar_t x, y, z, u, v, w; };
    };
};  // class TMatrixBase<6, 1>

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

/** Fixed size row major matrix. */
template<typename scalar_t, uint_t nx, uint_t ny>
class TMatrix final : public TMatrixBase<scalar_t, nx, ny> {
public:
    using TMatrixBase<scalar_t, nx, ny>::rows_data;

    /**************************************************************************/
    /**************************************************************************/

    /** Construct a constant matrix. */
    constexpr TMatrix(const scalar_t& scalar = {}) {
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            at(ix, iy) = scalar;
        });
    }
    /** Construct a matrix with other matrix. */
    /** @{ */
    constexpr TMatrix(const TMatrix& mat) {
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            at(ix, iy) = mat(ix, iy);
        });
    }
    template<uint_t mx_t, uint_t my_t>
    constexpr explicit TMatrix(const TMatrix<scalar_t, mx_t, my_t>& mat,
                                             const scalar_t& scalar = {}) {
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            at(ix, iy) = scalar;
        });
        for_n(0u, std::min(nx, mx_t),
              0u, std::min(ny, my_t), [&](uint_t ix, uint_t iy) {
            at(ix, iy) = mat(ix, iy);
        });
    }
    /** @} */

    /** Construct a matrix with data. */
    template<typename iter_t>
    constexpr explicit TMatrix(iter_t data_iter) {
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            at(ix, iy) = *data_iter++;
        });
    }
    /** Construct a matrix with elements. */
    constexpr TMatrix(std::initializer_list<scalar_t> data) {
        SKUNK_ASSERT(data.size() == nx * ny);
        auto data_iter = data.begin();
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            at(ix, iy) = *data_iter++;
        });
    }

    /**************************************************************************/
    /**************************************************************************/

    /** Get element at row and column index. */
    /** @{ */
    constexpr scalar_t& at(uint_t ix, uint_t iy = 0) {
        SKUNK_ASSERT(ix < nx && iy < ny);
        return rows_data[ix][iy];
    }
    constexpr const scalar_t& at(uint_t ix, uint_t iy = 0) const {
        SKUNK_ASSERT(ix < nx && iy < ny);
        return rows_data[ix][iy];
    }
    constexpr scalar_t& operator()(uint_t ix, uint_t iy = 0) {
        return at(ix, iy);
    }
    constexpr const scalar_t& operator()(uint_t ix, uint_t iy = 0) const {
        return at(ix, iy);
    }
    /** @} */

    /** Get row at index. */
    constexpr auto row(int_t ix) const {
        SKUNK_ASSERT(ix < nx);
        return TRowVector<scalar_t, ny>(rows_data[ix].begin());
    }

    /**************************************************************************/
    /**************************************************************************/

    /** Copy contents into the iterator. */
    template<typename iter_t>
    void store(iter_t data_iter) const {
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            *data_iter++ = at(ix, iy);
        });
    }

    /**************************************************************************/
    /**************************************************************************/

    /** Equality operator. */
    constexpr bool operator==(const TMatrix& mat) const {
        return for_n_reduce(0u, nx, 0u, ny, true, std::logical_and<>(), [&](uint_t ix, uint_t iy) {
            return at(ix, iy) == mat(ix, iy);
        });
    }
    /** Inequality operator. */
    constexpr bool operator!=(const TMatrix& mat) const {
        return for_n_reduce(0u, nx, 0u, ny, false, std::logical_or<>(), [&](uint_t ix, uint_t iy) {
            return at(ix, iy) != mat(ix, iy);
        });
    }

    /** Element-wise less then operator. */
    constexpr bool operator<(const TMatrix& mat) const {
        return for_n_reduce(0u, nx, 0u, ny, true, std::logical_and<>(), [&](uint_t ix, uint_t iy) {
            return at(ix, iy) < mat(ix, iy);
        });
    }
    /** Element-wise less then or equal operator. */
    constexpr bool operator<=(const TMatrix& mat) const {
        return for_n_reduce(0u, nx, 0u, ny, true, std::logical_and<>(), [&](uint_t ix, uint_t iy) {
            return at(ix, iy) <= mat(ix, iy);
        });
    }
    /** Element-wise greater then operator. */
    constexpr bool operator>(const TMatrix& mat) const {
        return for_n_reduce(0u, nx, 0u, ny, true, std::logical_and<>(), [&](uint_t ix, uint_t iy) {
            return at(ix, iy) > mat(ix, iy);
        });
    }
    /** Element-wise greater then or equal operator. */
    constexpr bool operator>=(const TMatrix& mat) const {
        return for_n_reduce(0u, nx, 0u, ny, true, std::logical_and<>(), [&](uint_t ix, uint_t iy) {
            return at(ix, iy) >= mat(ix, iy);
        });
    }

    /**************************************************************************/
    /**************************************************************************/

    /** Unary plus operator. */
    constexpr TMatrix operator+() const {
        TMatrix out(*this);
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            out(ix, iy) = +out(ix, iy);
        });
        return out;
    }
    /** Addition operator. */
    constexpr TMatrix operator+(const TMatrix& mat) const {
        TMatrix out(*this);
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            out(ix, iy) += mat(ix, iy);
        });
        return out;
    }
    /** Addition assignment operator. */
    constexpr TMatrix& operator+=(const TMatrix& mat) {
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            at(ix, iy) += mat(ix, iy);
        });
        return *this;
    }

    /** Negation operator. */
    constexpr TMatrix operator-() const {
        TMatrix out(*this);
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            out(ix, iy) = -out(ix, iy);
        });
        return out;
    }
    /** Subtraction operator. */
    constexpr TMatrix operator-(const TMatrix& mat) const {
        TMatrix out(*this);
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            out(ix, iy) -= mat(ix, iy);
        });
        return out;
    }
    /** Subtraction assignment operator. */
    constexpr TMatrix& operator-=(const TMatrix& mat) {
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            at(ix, iy) -= mat(ix, iy);
        });
        return *this;
    }

    /** Multiplication operator. */
    constexpr TMatrix operator*(const scalar_t& mat) const {
        TMatrix out(*this);
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            out(ix, iy) *= mat;
        });
        return out;
    }
    /** Multiplication assignment operator. */
    constexpr TMatrix& operator*=(const scalar_t& mat) {
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            at(ix, iy) *= mat;
        });
        return *this;
    }
    /** Element-wise multiplication operator. */
    constexpr TMatrix operator*(const TMatrix& mat) const {
        TMatrix out(*this);
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            out(ix, iy) *= mat(ix, iy);
        });
        return out;
    }
    /** Element-wise multiplication assignment operator. */
    constexpr TMatrix& operator*=(const TMatrix& mat) {
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            at(ix, iy) *= mat(ix, iy);
        });
        return *this;
    }

    /** Division operator. */
    constexpr TMatrix operator/(const scalar_t& val) const {
        TMatrix out(*this);
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            out(ix, iy) /= val;
        });
        return out;
    }
    /** Division assignment operator. */
    constexpr TMatrix& operator/=(const scalar_t& val) {
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            at(ix, iy) /= val;
        });
        return *this;
    }
    /** Element-wise division operator. */
    constexpr TMatrix operator/(const TMatrix& mat) const {
        TMatrix out(*this);
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            out(ix, iy) /= mat(ix, iy);
        });
        return out;
    }
    /** Element-wise division assignment operator. */
    constexpr TMatrix& operator/=(const TMatrix& mat) {
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            at(ix, iy) /= mat(ix, iy);
        });
        return *this;
    }

    /** Per-element minimum. */
    constexpr TMatrix min(const TMatrix& mat) const {
        TMatrix out(*this);
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            out(ix, iy) = std::min(out(ix, iy), mat(ix, iy));
        });
        return out;
    }
    /** Per-element maximum. */
    constexpr TMatrix max(const TMatrix& mat) const {
        TMatrix out(*this);
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            out(ix, iy) = std::max(out(ix, iy), mat(ix, iy));
        });
        return out;
    }

    /**************************************************************************/
    /**************************************************************************/

    /** Transpose of a matrix. */
    constexpr auto transpose() const {
        TMatrix<scalar_t, ny, nx> out;
        for_n(0u, nx, 0u, ny, [&](uint_t ix, uint_t iy) {
            out(iy, ix) = at(ix, iy);
        });
        return out;
    }

    /** Matrix product. */
    template<uint_t nz>
    constexpr auto prod(const TMatrix<scalar_t, ny, nz>& mat) const {
        TMatrix<scalar_t, nx, nz> product;
        for_n(0u, nx, 0u, ny, 0u, nz, [&](uint_t ix, uint_t iy, uint_t iz) {
            product(ix, iz) += at(ix, iy)*mat(iy, iz);
        });
        return product;
    }
    /** Matrix inner product (@f$A\cdot\B:=A^TB@f$). */
    template<uint_t nz>
    constexpr auto inner_prod(const TMatrix<scalar_t, nx, nz>& mat) const {
        TMatrix<scalar_t, ny, nz> product;
        for_n(0u, nx, 0u, ny, 0u, nz, [&](uint_t ix, uint_t iy, uint_t iz) {
            product(iy, iz) += at(ix, iy)*mat(ix, iz);
        });
        return product;
    }
    /** Matrix outer product (@f$A\otimes\B:=AB^T@f$). */
    template<uint_t nz>
    constexpr auto outer_prod(const TMatrix<scalar_t, nz, ny>& mat) const {
        TMatrix<scalar_t, nx, nz> product;
        for_n(0u, nx, 0u, ny, 0u, nz, [&](uint_t ix, uint_t iy, uint_t iz) {
            product(ix, iz) += at(ix, iy)*mat(iz, iy);
        });
        return product;
    }

    /**************************************************************************/
    /**************************************************************************/

    /** Sum of matrix elements. */
    constexpr scalar_t sum_elements() const {
        return for_n_reduce(0u, nx, 0u, ny, scalar_t(0), std::plus<>(), *this);
    }
    /** Product of matrix elements. */
    constexpr scalar_t prod_elements() const {
        return for_n_reduce(0u, nx, 0u, ny, scalar_t(1), std::multiplies<>(), *this);
    }

    /** Minimum matrix element. */
    constexpr scalar_t min_element() const {
        return for_n_reduce(0u, nx, 0u, ny, std::numeric_limits<scalar_t>::max(),
                            [](scalar_t a, scalar_t b) { return std::min(a, b); }, *this);
    }
    /** Maximum matrix element. */
    constexpr scalar_t max_element() const {
        return for_n_reduce(0u, nx, 0u, ny, std::numeric_limits<scalar_t>::min(),
                            [](scalar_t a, scalar_t b) { return std::max(a, b); }, *this);
    }
};  // class TMatrix

/** Scalar-matrix multiplication operator. */
template<typename scalar_t, uint_t nx, uint_t ny>
TMatrix<scalar_t, nx, ny> operator*(const scalar_t& val,
                                                 const TMatrix<scalar_t, nx, ny>& mat) {
    return mat*val;
}   // operator*

/**************************************************************************/
/**************************************************************************/

/** Input a matrix. */
template<typename scalar_t, uint_t nx_t, uint_t ny_t>
std::istream& operator>>(std::istream& in,
                                      TMatrix<scalar_t, nx_t, ny_t>& mat) {
    for_n(0u, nx_t, 0u, ny_t, [&](uint_t ix, uint_t iy) {
        in >> mat(ix, iy);
    });
    return in;
}   // operator>>
/** Output a matrix. */
template<typename scalar_t, uint_t nx_t, uint_t ny_t>
std::ostream& operator<<(std::ostream& out,
                                      const TMatrix<scalar_t, nx_t, ny_t>& mat) {
    for_n(0u, nx_t, 0u, ny_t, [&](uint_t ix, uint_t iy) {
        out << mat(ix, iy);
        if (ix == nx_t - 1) {
            out << std::endl;
        } else {
            out << ' ';
        }
    });
    return out;
}   // operator<<

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif  // ifndef GEOM_MATRIX_HH_
