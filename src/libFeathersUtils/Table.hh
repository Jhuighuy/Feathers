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
 * Compressed sparse row table class.
 */
class cTable {
private:
    std::vector<uint_t> m_row_pointers{0};
    std::vector<uint_t> m_column_indices;

public:

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Number of rows in the mesh. */
    uint_t num_rows() const {
        return m_row_pointers.size() - 1;
    }

    /** Pointer to the beginning of the row. */
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_row, (uint_t row_index), {
        FEATHERS_ASSERT(row_index < num_rows());
        return m_column_indices.data() + m_row_pointers[row_index];
    })

    /** Pointer to the end of the row. */
    FEATHERS_CONST_OVERLOAD(uint_t*, end_row, (uint_t row_index), {
        FEATHERS_ASSERT(row_index < num_rows());
        return m_column_indices.data() + m_row_pointers[row_index + 1];
    })

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Insert a row into the table. */
    template<typename tIter>
    void emplace_back_row(tIter first_index, tIter last_index) {
        m_column_indices.insert(
            m_column_indices.end(), first_index, last_index);
        m_row_pointers.push_back(m_column_indices.size());
    }

    /** Insert a dummy row into the table. */
    void emplace_back_row(uint_t num_column_indices= 0, uint_t column_index= npos) {
        m_column_indices.insert(
            m_column_indices.end(), num_column_indices, column_index);
        m_row_pointers.push_back(m_column_indices.size());
    }
};  // class cTable

template<typename tIter>
void permute_tables(tIter first_permutation_iter,
                    tIter last_permutation_iter, cTable& table) {
    cTable permuted_table;
    for (tIter permutation_iter = first_permutation_iter;
         permutation_iter != last_permutation_iter; ++permutation_iter) {
        permuted_table.emplace_back_row(
            table.begin_row(*permutation_iter), table.end_row(*permutation_iter));
    }
    table = std::move(permuted_table);
}
template<typename tIter, typename... tTable>
void permute_tables(tIter first_permutation_iter,
                    tIter last_permutation_iter, cTable& table, tTable&... rest) {
    permute_tables(first_permutation_iter, last_permutation_iter, table);
    permute_tables(first_permutation_iter, last_permutation_iter, rest...);
}

} // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif // TABLE_HH_
