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

#include "SkunkBase.hh"
#include "stormUtils/Index.hxx"

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Compressed sparse row (CSR) table class.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class RowIndex, class ColumnIndex>
class CsrTable {
private:

  Vector<size_t, RowIndex> m_row_offsets{0};
  Vector<ColumnIndex, size_t> m_column_indices;

public:

  size_t num_rows() const {
    return m_row_offsets.size() - 1;
  }

  auto operator[](RowIndex row_index) noexcept {
    STORM_ASSERT_(row_index < num_rows() && "Row index is out of range.");
    return ranges::subrange(&m_column_indices[m_row_offsets[row_index]],
                            &m_column_indices[m_row_offsets[row_index + 1]]);
  }
  auto operator[](RowIndex row_index) const noexcept {
    STORM_ASSERT_(row_index < num_rows() && "Row index is out of range.");
    return ranges::subrange(&m_column_indices[m_row_offsets[row_index]],
                            &m_column_indices[m_row_offsets[row_index + 1]]);
  }

  void insert(RowIndex rowIndex, ColumnIndex columnIndex) {
    STORM_ASSERT_(rowIndex < num_rows());
    m_column_indices.insert(
        m_column_indices.begin() + m_row_offsets[rowIndex + 1], columnIndex);
    std::for_each(m_row_offsets.begin() + (size_t) rowIndex + 1,
                  m_row_offsets.end(), [](size_t& offset) { offset += 1; });
  }

  void insert_row(size_t num_column_indices = 0,
                  ColumnIndex column_index = ColumnIndex{npos}) {
    m_column_indices.insert( //
        m_column_indices.end(), num_column_indices, column_index);
    m_row_offsets.push_back(m_column_indices.size());
  }
  template<class Iterator>
  void insert_row(Iterator column_indices_begin, Iterator column_indices_end) {
    m_column_indices.insert( //
        m_column_indices.end(), column_indices_begin, column_indices_end);
    m_row_offsets.push_back(m_column_indices.size());
  }
  void insert_row(ranges::range auto range) {
    insert_row(range.begin(), range.end());
  }

  /** Pointer to the beginning of the row. */
  ConstOverload(ColumnIndex*, begin_row, (RowIndex row_index),
                {
                  FEATHERS_ASSERT(row_index < num_rows());
                  return &m_column_indices[m_row_offsets[row_index]];
                })

      /** Pointer to the End of the row. */
      ConstOverload(ColumnIndex*, end_row, (RowIndex row_index), {
        FEATHERS_ASSERT(row_index < num_rows());
        return &m_column_indices[m_row_offsets[row_index + 1]];
      })

}; // class CsrTable

template<typename tIter, typename tTable>
void permute_rows(tIter first_permutation_iter, tIter last_permutation_iter,
                  tTable& table) {
  tTable permuted_table;
  for (tIter permutation_iter = first_permutation_iter;
       permutation_iter != last_permutation_iter; ++permutation_iter) {
    permuted_table.insert_row(table.begin_row(*permutation_iter),
                              table.end_row(*permutation_iter));
  }
  table = std::move(permuted_table);
}
template<typename tIter, typename tTable, typename... tTables>
void permute_rows(tIter first_permutation_iter, tIter last_permutation_iter,
                  tTable& table, tTables&... rest) {
  permute_rows(first_permutation_iter, last_permutation_iter, table);
  permute_rows(first_permutation_iter, last_permutation_iter, rest...);
}

} // namespace Storm
