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
#ifndef TABLE_HH_
#define TABLE_HH_

#include "SkunkBase.hh"

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/**
 * Dummy table class.
 */
class cDummyTable {
private:
    uint_t m_num_rows = 0;

public:

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Number of rows in the mesh. */
    uint_t num_rows() const {
        return m_num_rows;
    }

    /** Pointer to the beginning of the row. */
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_row, (uint_t row_index), {
        FEATHERS_ASSERT(row_index < num_rows());
        return nullptr;
    })

    /** Pointer to the end of the row. */
    FEATHERS_CONST_OVERLOAD(uint_t*, end_row, (uint_t row_index), {
        FEATHERS_ASSERT(row_index < num_rows());
        return nullptr;
    })

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    void reserve_rows(uint_t FEATHERS_NOT_USED(row_capacity)) {
    }

    /** Insert a row into the table. */
    /** @{ */
    void emplace_back_row(uint_t num_column_indices = 0, uint_t FEATHERS_NOT_USED(column_index) = npos) {
        FEATHERS_ASSERT(num_column_indices == 0);
        m_num_rows += 1;
    }
    template<typename tColumnIndexIter>
    void emplace_back_row(tColumnIndexIter first_index_iter, tColumnIndexIter last_index_iter) {
        FEATHERS_ASSERT(first_index_iter == last_index_iter);
        m_num_rows += 1;
    }
    /** @} */
};  // class cDummyTable

/**
 * Compressed sparse row (CSR) table class.
 */
class cCSRTable {
private:
    std::vector<uint_t> m_row_offsets{0};
    std::vector<uint_t> m_column_indices;

public:

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Number of rows in the mesh. */
    uint_t num_rows() const {
        return m_row_offsets.size() - 1;
    }

    /** Pointer to the beginning of the row. */
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_row, (uint_t row_index), {
        FEATHERS_ASSERT(row_index < num_rows());
        return &m_column_indices[m_row_offsets[row_index]];
    })

    /** Pointer to the end of the row. */
    FEATHERS_CONST_OVERLOAD(uint_t*, end_row, (uint_t row_index), {
        FEATHERS_ASSERT(row_index < num_rows());
        return &m_column_indices[m_row_offsets[row_index + 1]];
    })

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    void reserve_rows(uint_t row_capacity) {
        m_row_offsets.reserve(row_capacity + 1);
    }

    /** Insert a row into the table. */
    /** @{ */
    void emplace_back_row(uint_t num_column_indices = 0, uint_t column_index = npos) {
        m_column_indices.insert(
            m_column_indices.end(), num_column_indices, column_index);
        m_row_offsets.push_back(m_column_indices.size());
    }
    template<typename tColumnIndexIter>
    void emplace_back_row(tColumnIndexIter first_index_iter, tColumnIndexIter last_index_iter) {
        m_column_indices.insert(
            m_column_indices.end(), first_index_iter, last_index_iter);
        m_row_offsets.push_back(m_column_indices.size());
    }
    /** @} */
};  // class cCSRTable

/**
 * Table in some weird format.
 */
class cWCSRTable {
private:
    std::vector<uint_t> m_row_offsets;
    std::vector<uint_t> m_column_nums_and_indices;

public:

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Number of rows in the mesh. */
    uint_t num_rows() const {
        return m_row_offsets.size();
    }

    /** Pointer to the beginning of the row. */
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_row, (uint_t row_index), {
        FEATHERS_ASSERT(row_index < num_rows());
        const uint_t row_offset = m_row_offsets[row_index];
        return &m_column_nums_and_indices[row_offset + 1];
    })

    /** Pointer to the end of the row. */
    FEATHERS_CONST_OVERLOAD(uint_t*, end_row, (uint_t row_index), {
        FEATHERS_ASSERT(row_index < num_rows());
        const uint_t row_offset = m_row_offsets[row_index];
        const uint_t row_num_column_indices = m_column_nums_and_indices[row_offset];
        return &m_column_nums_and_indices[row_offset + 1 + row_num_column_indices];
    })

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    void reserve_rows(uint_t row_capacity) {
        m_row_offsets.reserve(row_capacity);
    }

    /** Insert a row into the table. */
    /** @{ */
    void emplace_back_row(uint_t num_column_indices = 0, uint_t column_index = npos) {
        m_row_offsets.push_back(m_column_nums_and_indices.size());
        m_column_nums_and_indices.emplace_back(num_column_indices);
        m_column_nums_and_indices.insert(
            m_column_nums_and_indices.end(), num_column_indices, column_index);
    }
    template<typename tIter>
    void emplace_back_row(tIter first_index, tIter last_index) {
        m_row_offsets.push_back(m_column_nums_and_indices.size());
        m_column_nums_and_indices.emplace_back(last_index - first_index);
        m_column_nums_and_indices.insert(
            m_column_nums_and_indices.end(), first_index, last_index);
    }
    /** @} */
};  // class cWCSRTable

using cTable = cCSRTable;

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

template<typename tIter, typename tTable>
void permute_rows(tIter first_permutation_iter,
                  tIter last_permutation_iter, tTable& table) {
    tTable permuted_table;
    for (tIter permutation_iter = first_permutation_iter;
         permutation_iter != last_permutation_iter; ++permutation_iter) {
        permuted_table.emplace_back_row(
            table.begin_row(*permutation_iter), table.end_row(*permutation_iter));
    }
    table = std::move(permuted_table);
}
template<typename tIter, typename tTable, typename... tTables>
void permute_rows(tIter first_permutation_iter,
                  tIter last_permutation_iter, tTable& table, tTables&... rest) {
    permute_rows(first_permutation_iter, last_permutation_iter, table);
    permute_rows(first_permutation_iter, last_permutation_iter, rest...);
}

} // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif // TABLE_HH_
