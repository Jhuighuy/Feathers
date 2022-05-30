/*
 *  ______  ______   ______   ______  __  __   ______   ______   ______
 * /\  ___\/\  ___\ /\  __ \ /\__  _\/\ \_\ \ /\  ___\ /\  __ \ /\  ___\
 * \ \  __\\ \  _\  \ \  __ \\/_/\ \/\ \  __ \\ \  __\ \ \  __/ \ \___  \
 *  \ \_\   \ \_____\\ \_\ \_\  \ \_\ \ \_\ \_\\ \_____\\ \_\ \_\\/\_____\
 *   \/_/    \/_____/ \/_/\/_/   \/_/  \/_/\/_/ \/_____/ \/_/ /_/ \/_____/
 *
 * Copyright (c) 2021 Oleg Butakov
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
#include "libFeathersMesh/Index.hxx"

namespace feathers {

/**
 * Compressed sparse row (CSR) table class.
 */
template<class RowIndex, class ColumnIndex>
class CsrTable {
private:
  IndexedVector<size_t, RowIndex> m_row_offsets{0};
  IndexedVector<ColumnIndex, size_t> m_column_indices;

public:

  // ---------------------------------------------------------------------- //
  // ---------------------------------------------------------------------- //

  /** Number of rows in the mesh. */
  size_t num_rows() const {
    return m_row_offsets.size() - 1;
  }

  void Clear() {
    m_row_offsets = {0};
    m_column_indices.clear();
  }

  void insert(RowIndex rowIndex, ColumnIndex columnIndex) {
    StormAssert(rowIndex < num_rows());
    m_column_indices.insert(
      m_column_indices.begin() + m_row_offsets[rowIndex + 1], columnIndex);
    std::for_each(m_row_offsets.begin() + (size_t)rowIndex + 1, m_row_offsets.end(),
      [](size_t& offset) { offset += 1; });
  }

  /** Pointer to the beginning of the row. */
  ConstOverload(ColumnIndex*, begin_row, (RowIndex row_index), {
    FEATHERS_ASSERT(row_index < num_rows());
    return &m_column_indices[m_row_offsets[row_index]];
  })
  ConstOverload(ColumnIndex*, begin_row, (size_t row_index), requires(!std::is_same_v<RowIndex, size_t>) {
    FEATHERS_ASSERT(row_index < num_rows());
    return &m_column_indices[m_row_offsets[RowIndex(row_index)]];
  })

  /** Pointer to the End of the row. */
  ConstOverload(ColumnIndex*, end_row, (RowIndex row_index), {
    FEATHERS_ASSERT(row_index < num_rows());
    return &m_column_indices[m_row_offsets[row_index + 1]];
  })
  ConstOverload(ColumnIndex*, end_row, (size_t row_index), requires(!std::is_same_v<RowIndex, size_t>) {
    FEATHERS_ASSERT(row_index < num_rows());
    return &m_column_indices[m_row_offsets[RowIndex(row_index) + 1]];
  })

  auto operator[](RowIndex row_index) noexcept {
    FEATHERS_ASSERT(row_index < num_rows());
    return ranges::subrange(
      &m_column_indices[m_row_offsets[row_index]],
      &m_column_indices[m_row_offsets[row_index + 1]]);
  }
  auto operator[](RowIndex row_index) const noexcept requires(!std::is_same_v<RowIndex, size_t>) {
    FEATHERS_ASSERT(row_index < num_rows());
    return ranges::subrange(
      &m_column_indices[m_row_offsets[row_index]],
      &m_column_indices[m_row_offsets[row_index + 1]]);
  }

  // ---------------------------------------------------------------------- //
  // ---------------------------------------------------------------------- //

  void clear() {
    m_row_offsets.assign(1, 0);
    m_column_indices.clear();
  }

  /** Insert a row into the table. */
  /** @{ */
  void emplaceRow(size_t num_column_indices = 0, ColumnIndex column_index = ColumnIndex{npos}) {
    m_column_indices.insert(
      m_column_indices.end(), num_column_indices, column_index);
    m_row_offsets.push_back(m_column_indices.size());
  }
  template<typename tColumnIndexIter>
  void emplaceRow(tColumnIndexIter first_index_iter, tColumnIndexIter last_index_iter) {
    m_column_indices.insert(
      m_column_indices.end(), first_index_iter, last_index_iter);
    m_row_offsets.push_back(m_column_indices.size());
  }
  void emplaceRow(ranges::range auto range) {
    emplaceRow(range.begin(), range.end());
  }
  /** @} */
};  // class CsrTable

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

template<typename tIter, typename tTable>
void permute_rows(tIter first_permutation_iter,
                  tIter last_permutation_iter, tTable& table) {
  tTable permuted_table;
  for (tIter permutation_iter = first_permutation_iter;
       permutation_iter != last_permutation_iter; ++permutation_iter) {
    permuted_table.emplaceRow(
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

#endif // TABLE_HH_
