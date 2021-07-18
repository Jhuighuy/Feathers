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
#ifndef SKUNK_MISC_REORDERING_HH
#define SKUNK_MISC_REORDERING_HH

#include <SkunkBase.hh>

#include <algorithm>

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/**
 * Inverse reordering. Complexity is linear.
 */
template<typename tIter, typename tInverseMutableIter>
void permutation_to_indices(tIter index_iter, tIter last_index_iter,
                            tInverseMutableIter inverse_index_iter) {
    size_t index = 0;
    for (; index_iter != last_index_iter; ++index_iter, ++index) {
        *(inverse_index_iter + *index_iter) = index;
    }
}   // permutation_to_indices

/**
 * Apply reordering.
 * Index iterators must be mutable. Complexity is linear.
 */
/** @{ */
template<typename tMutableIter, typename tSwapFunc>
void permute_swap(tMutableIter index_iter,
                  tMutableIter last_index_iter, tSwapFunc swap) {
    /* For implementation details see
     * https://devblogs.microsoft.com/oldnewthing/20170102-00/?p=95095 */
    size_t index = 0;
    for (; index_iter != last_index_iter; ++index_iter, ++index) {
        size_t current_index = index;
        tMutableIter current_index_iter = index_iter;
        while (*current_index_iter != index) {
            const size_t new_index = *current_index_iter;
            swap(current_index, new_index);
            *current_index_iter = current_index;
            current_index = new_index;
            current_index_iter = index_iter + (current_index - index);
        }
        *current_index_iter = current_index;
    }
}   // permute_swap
template<typename tMutableIter, typename tValueMutableIter>
void permute_inplace(tMutableIter first_index_iter,
                     tMutableIter last_index_iter, tValueMutableIter iter) {
    permute_swap(first_index_iter, last_index_iter, [&](auto index, auto new_index) {
        std::iter_swap(iter + index, iter + new_index);
    });
}   // permute_swap
template<typename tContainer, typename tValueIter>
void permute_inplace(tContainer indices, tValueIter iter) {
    permute_inplace(indices.begin(), indices.end(), iter);
}   // permute_inplace
/** @} */

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif  // ifndef SKUNK_MISC_REORDERING_HH
