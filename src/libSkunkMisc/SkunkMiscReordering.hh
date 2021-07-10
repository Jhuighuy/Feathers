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

/** Inverse reordering.
 ** Complexity is linear. */
template<typename TIter, typename TInvIter>
SKUNK_INLINE void inverse_reordering(TIter index_iter,
                                     TIter index_iter_end,
                                     TInvIter inv_index_iter) {
    size_t index = 0;
    for (; index_iter != index_iter_end; ++index_iter, ++index) {
        *(inv_index_iter + *index_iter) = index;
    }
}   // inverse_reordering

/** Apply reordering.
 ** Index iterators must be mutable. Complexity is linear. */
/** @{ */
template<typename TIter, typename TSwapFunc>
SKUNK_INLINE void reorder_swap(TIter index_iter,
                               TIter index_iter_end,
                               TSwapFunc&& swap) {
    /* For implementation details see
     * https://devblogs.microsoft.com/oldnewthing/20170102-00/?p=95095 */
    size_t index = 0;
    for (; index_iter != index_iter_end; ++index_iter, ++index) {
        size_t cur_index = index;
        TIter cur_index_iter = index_iter;
        while (*cur_index_iter != index) {
            const size_t new_index = *cur_index_iter;
            swap(cur_index, new_index);
            *cur_index_iter = cur_index;
            cur_index = new_index;
            cur_index_iter = index_iter + (cur_index - index);
        }
        *cur_index_iter = cur_index;
    }
}   // reorder_swap
template<typename TIter, typename TValIter>
SKUNK_INLINE void reorder(TIter ind_iter,
                          TIter ind_iter_end,
                          TValIter iter) {
    reorder_swap(ind_iter, ind_iter_end, [&](auto ind, auto new_ind) {
        std::iter_swap(iter + ind, iter + new_ind);
    });
}   // reorder
/** @} */

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif  // ifndef SKUNK_MISC_REORDERING_HH
