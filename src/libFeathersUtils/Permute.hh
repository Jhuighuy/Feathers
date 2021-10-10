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
#ifndef PERMUTATION_HH_
#define PERMUTATION_HH_

#include <SkunkBase.hh>

namespace feathers {

/**
 * Convert permutation to ordering. Complexity is linear.
 *
 * Beware that @c permutation is not the same as the @c ordering:
 * @verbatim
 *  new_values[i] = old_values[permutation[i]],
 *  new_values[ordering[i]] = old_values[i].
 * @endverbatim
 *
 * Consider the permutation is:
 * @verbatim
 *  {1, 3, 2, 0, 4, 5},
 * @endverbatim
 * For the considered permutation, ordering is:
 * @verbatim
 *  {3, 0, 2, 1, 4, 5}.
 * @endverbatim
 */
template<typename tPermutationMutableIter, typename tIndexMutableIter>
void convert_permutation_to_ordering(tPermutationMutableIter first_permutation_iter,
                                     tPermutationMutableIter last_permutation_iter,
                                     tIndexMutableIter first_index_iter) {
    using index_t = std::decay_t<decltype(*first_permutation_iter)>;
    index_t index = 0;
    for (tPermutationMutableIter permutation_iter = first_permutation_iter;
            permutation_iter != last_permutation_iter; ++permutation_iter, ++index) {
        FEATHERS_ASSERT(*permutation_iter != npos);
        *(first_index_iter + *permutation_iter) = index;
    }
}   // convert_permutation_to_ordering

/**
 * Permute iterators.
 * Index iterators must be mutable. Complexity is linear.
 *
 * Example usage:
 * @code
 * std::vector<T> values;
 * std::vector<size_t> permutation;
 * permute_inplace_swap(
 *     permutation.begin(), permutation.end(),
 *     [&](size_t index_1, size_t index_2) { std::swap(values[index_1], values[index_2]); });
 * @endcode
 */
template<typename tPermutationMutableIter, typename tSwapFunc>
void permute_inplace_swap(tPermutationMutableIter first_permutation_iter,
                          tPermutationMutableIter last_permutation_iter, tSwapFunc&& swap_func) {
    /* For implementation details see
     * https://devblogs.microsoft.com/oldnewthing/20170102-00/?p=95095 */
    using index_t = std::decay_t<decltype(*first_permutation_iter)>;
    index_t index = 0;
    for (tPermutationMutableIter permutation_iter = first_permutation_iter;
            permutation_iter != last_permutation_iter; ++permutation_iter, ++index) {
        FEATHERS_ASSERT(*permutation_iter != npos);
        index_t current_index = index;
        tPermutationMutableIter current_permutation_iter = permutation_iter;
        while (*current_permutation_iter != index) {
            FEATHERS_ASSERT(*current_permutation_iter != npos);
            const index_t new_index = *current_permutation_iter;
            swap_func(current_index, new_index);
            *current_permutation_iter = current_index;
            current_index = new_index;
            current_permutation_iter = permutation_iter + (current_index - index);
        }
        *current_permutation_iter = current_index;
    }
#ifndef NDEBUG
    std::fill(first_permutation_iter, last_permutation_iter, npos);
#endif
}   // permute_inplace_swap

/**
 * Permute iterators.
 * Index iterators must be mutable. Complexity is linear.
 *
 * Example usage:
 * @code
 * std::vector<T> values_1, values_2, values_3;
 * std::vector<size_t> permutation;
 * permute_inplace(
 *     permutation.begin(), permutation.end(),
 *     values_1.begin(), values_2.begin(), values_3.begin());
 * @endcode
 */
template<typename tPermutationMutableIter, typename... tValueMutableIter>
void permute_inplace(tPermutationMutableIter first_permutation_iter,
                     tPermutationMutableIter last_permutation_iter,
                     tValueMutableIter... first_value_iter) {
    permute_inplace_swap(first_permutation_iter,
                         last_permutation_iter, [&](auto index_1, auto index_2) {
        /* A dirty hack to perform the swap. */
        auto references_1 = std::tie(*(first_value_iter + index_1)...);
        auto references_2 = std::tie(*(first_value_iter + index_2)...);
        std::swap(references_1, references_2);
    });
}   // permute_inplace

}   // namespace feathers

#endif  // ifndef PERMUTATION_HH_
