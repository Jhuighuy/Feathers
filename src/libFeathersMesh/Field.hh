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
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#pragma once
#ifndef FIELD_HH_
#define FIELD_HH_

#include "SkunkBase.hh"

namespace feathers {

template<typename data_t>
class tGenericSubField {
private:
    uint_t m_num_vars;
    data_t* m_elements;

public:
    tGenericSubField(uint_t num_vars, data_t* elements):
        m_num_vars(num_vars), m_elements(elements) {
    }

    template<typename data_u = data_t,
        typename = std::enable_if_t<!std::is_const_v<data_u>>>
    tGenericSubField& operator=(const std::initializer_list<data_u>& other) {
        if (other.size() == 0) {
            std::fill_n(m_elements, m_num_vars, data_u{});
        } else {
            std::copy(other.begin(), other.end(), m_elements);
        }
        return *this;
    }

    template<typename type_u = data_t>
    std::enable_if_t<!std::is_const_v<type_u>> fill(const type_u& value) {
        std::fill_n(m_elements, m_num_vars, value);
    }

    auto data() const {
        return m_elements;
    }

    auto& operator[](uint_t i) const {
        return m_elements[i];
    }
};

template<typename type_t, typename component_type_t = type_t,
    uint_t num_components = sizeof(type_t)/sizeof(component_type_t)>
class tGenericField {
private:
    uint_t m_num_vars;
    std::vector<component_type_t> m_elements;

public:
    tGenericField(uint_t num_vars, uint_t num_elements):
        m_num_vars(num_vars), m_elements(num_components*num_vars*num_elements) {
    }

    auto operator[](uint_t element_index) {
        return tGenericSubField<type_t>(m_num_vars,
            reinterpret_cast<type_t*>(&m_elements[num_components*m_num_vars*element_index]));
    }
    auto operator[](uint_t element_index) const {
        return tGenericSubField<const type_t>(m_num_vars,
            reinterpret_cast<const type_t*>(&m_elements[num_components*m_num_vars*element_index]));
    }

    void swap(tGenericField& other) {
        std::swap(m_num_vars, other.m_num_vars);
        std::swap(m_elements, other.m_elements);
    }
}; // class tGenericField

template<uint_t = 0>
using tScalarSubField = tGenericSubField<real_t>;
template<uint_t = 0>
using tVectorSubField = tGenericSubField<vec3_t>;
template<uint_t = 0>
using tMatrixSubField = tGenericSubField<mat3_t>;

template<uint_t = 0>
using tScalarConstSubField = tGenericSubField<const real_t>;
template<uint_t = 0>
using tVectorConstSubField = tGenericSubField<const vec3_t>;
template<uint_t = 0>
using tMatrixConstSubField = tGenericSubField<const mat3_t>;

template<uint_t = 0>
using tScalarField = tGenericField<real_t>;
template<uint_t = 0>
using tVectorField = tGenericField<vec3_t, real_t>;
template<uint_t = 0>
using tMatrixField = tGenericField<mat3_t, real_t>;

} // namespace feathers

#endif // FIELD_HH_
