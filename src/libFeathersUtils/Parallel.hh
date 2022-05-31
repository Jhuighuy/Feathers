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
#ifndef PARALLEL_HH_
#define PARALLEL_HH_

#include <SkunkBase.hh>

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/**
 * Get index of the current thread.
 */
static size_t get_thread_index();
/**
 * Get maximum number of threads.
 */
static size_t get_max_num_threads();
/**
 * Set maximum number of threads.
 */
static void set_max_num_threads(size_t num_threads);

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Run a range in parallel.
 */
/** @{ */
template<typename tIndex, typename tFunc>
void for_range(tIndex first, tIndex last, tFunc func);
template<typename tIndex, typename tFunc>
void for_range(tIndex first_1, tIndex last_1,
               tIndex first_2, tIndex last_2, tFunc func);
template<typename tIndex, typename tFunc>
void for_range(tIndex first_1, tIndex last_1,
               tIndex first_2, tIndex last_2,
               tIndex first_3, tIndex last_3, tFunc func);
template<class Range, class Func>
void for_each(Range&& range, Func func) {
  for_range(std::begin(range), std::end(range), [&](auto x) { func(*x); });
}
/** @} */

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Reduce a range in parallel.
 */
/** @{ */
template<typename tValue, typename tIndex, typename tFunc, typename tReduceFunc>
tValue for_range_reduce(tIndex first, tIndex last,
                        tValue init, tFunc func, tReduceFunc reduce_func);
template<typename tValue, typename tIndex, typename tFunc, typename tReduceFunc>
tValue for_range_reduce(tIndex first_1, tIndex last_1,
                        tIndex first_2, tIndex last_2,
                        tValue init, tFunc func, tReduceFunc reduce_func);
template<typename tValue, typename tIndex, typename tFunc, typename tReduceFunc>
tValue for_range_reduce(tIndex first_1, tIndex last_1,
                        tIndex first_2, tIndex last_2,
                        tIndex first_3, tIndex last_3,
                        tValue init, tFunc func, tReduceFunc reduce_func);
/** @} */

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Compute sum of the range values in parallel.
 */
/** @{ */
template<typename tValue, typename tIndex, typename tFunc>
tValue for_range_sum(tIndex first, tIndex last, tValue init, tFunc func);
template<typename tValue, typename tIndex, typename tFunc>
tValue for_range_sum(tIndex first_1, tIndex last_1,
                     tIndex first_2, tIndex last_2, tValue init, tFunc func);
template<typename tValue, typename tIndex, typename tFunc>
tValue for_range_sum(tIndex first_1, tIndex last_1,
                     tIndex first_2, tIndex last_2,
                     tIndex first_3, tIndex last_3, tValue init, tFunc func);
template<typename tValue, class Range, class Func>
auto ForEachSum(Range&& range, tValue init, Func func) {
  return for_range_sum(std::begin(range), std::end(range), init, [&](auto x) { return func(*x); });
}
/** @} */

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Compute minimum of the range values in parallel.
 */
/** @{ */
template<typename tValue, typename tIndex, typename tFunc>
tValue for_range_min(tIndex first, tIndex last, tValue init, tFunc func);
template<typename tValue, typename tIndex, typename tFunc>
tValue for_range_min(tIndex first_1, tIndex last_1,
                     tIndex first_2, tIndex last_2, tValue init, tFunc func);
template<typename tValue, typename tIndex, typename tFunc>
tValue for_range_min(tIndex first_1, tIndex last_1,
                     tIndex first_2, tIndex last_2,
                     tIndex first_3, tIndex last_3, tValue init, tFunc func);
/** @} */

/**
 * Compute maximum of the range values in parallel.
 */
/** @{ */
template<typename tValue, typename tIndex, typename tFunc>
tValue for_range_max(tIndex first, tIndex last, tValue init, tFunc func);
template<typename tValue, typename tIndex, typename tFunc>
tValue for_range_max(tIndex first_1, tIndex last_1,
                     tIndex first_2, tIndex last_2, tValue init, tFunc func);
template<typename tValue, typename tIndex, typename tFunc>
tValue for_range_max(tIndex first_1, tIndex last_1,
                     tIndex first_2, tIndex last_2,
                     tIndex first_3, tIndex last_3, tValue init, tFunc func);
/** @} */

/**
 * Compute minimum and maximum of the range values in parallel.
 */
/** @{ */
template<typename tValue, typename tIndex, typename tFunc>
std::pair<tValue, tValue> for_range_minmax(tIndex first, tIndex last,
                                           tValue min_init, tValue max_init, tFunc func);
template<typename tValue, typename tIndex, typename tFunc>
std::pair<tValue, tValue> for_range_minmax(tIndex first_1, tIndex last_1,
                                           tIndex first_2, tIndex last_2,
                                           tValue min_init, tValue max_init, tFunc func);
template<typename tValue, typename tIndex, typename tFunc>
std::pair<tValue, tValue> for_range_minmax(tIndex first_1, tIndex last_1,
                                           tIndex first_2, tIndex last_2,
                                           tIndex first_3, tIndex last_3,
                                           tValue min_init, tValue max_init, tFunc func);
template<typename tValue, class Range, class Func>
auto for_each_min_max(Range&& range, tValue min_init, tValue max_init, Func func) {
  return for_range_minmax(std::begin(range), std::end(range), min_init, max_init, [&](auto x) { return func(*x); });
}
/** @} */

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

/* Try to use some existing implementation, */
#if FEATHERS_HAS_TBB
#include "Parallel_TBB.inl"
#elif FEATHERS_HAS_OPENMP
#include "Parallel_OpenMP.inl"
#endif
/* or fallback to the generic implementation. */

namespace feathers {

#ifndef THREAD_IDS_DEFINED_
size_t get_thread_index() {
    return 0;
}
size_t get_max_num_threads() {
    return 1;
}
void set_max_num_threads(size_t FEATHERS_NOT_USED(num_threads)) {
}
#endif
#undef THREAD_IDS_DEFINED_

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

#ifndef FOR_RANGE_1_DEFINED_
template<typename tIndex, typename tFunc>
void for_range(tIndex first, tIndex last, tFunc func) {
    for (tIndex index = first; index < last; ++index) {
        func(index);
    }
}
#endif
#undef FOR_RANGE_1_DEFINED_
#ifndef FOR_RANGE_2_DEFINED_
template<typename tIndex, typename tFunc>
void for_range(tIndex first_1, tIndex last_1,
               tIndex first_2, tIndex last_2, tFunc func) {
    for_range(first_1, last_1, [&](tIndex index_1) {
        for (tIndex index_2 = first_2; index_2 < last_2; ++index_2) {
            func(index_1, index_2);
        }
    });
}
#endif
#undef FOR_RANGE_2_DEFINED_
#ifndef FOR_RANGE_3_DEFINED_
template<typename tIndex, typename tFunc>
void for_range(tIndex first_1, tIndex last_1,
               tIndex first_2, tIndex last_2,
               tIndex first_3, tIndex last_3, tFunc func) {
    for_range(first_1, last_1,
              first_2, last_2, [&](tIndex index_1, tIndex index_2) {
        for (tIndex index_3 = first_3; index_3 < last_3; ++index_3) {
            func(index_1, index_2, index_3);
        }
    });
}
#endif
#undef FOR_RANGE_3_DEFINED_

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

#ifndef FOR_RANGE_REDUCE_1_DEFINED_
template<typename tValue, typename tIndex, typename tFunc, typename tReduceFunc>
tValue for_range_reduce(tIndex first, tIndex last,
                        tValue init, tFunc func, tReduceFunc reduce_func) {
    std::vector<tValue> partial(get_max_num_threads(), init);
    for_range(first, last, [&](tIndex index) {
        partial[get_thread_index()] =
            reduce_func(partial[get_thread_index()], func(index));
    });
    init = std::accumulate(
        partial.begin(), partial.end() - 1, partial.back(), reduce_func);
    return init;
}
#endif
#undef FOR_RANGE_REDUCE_1_DEFINED_
#ifndef FOR_RANGE_REDUCE_2_DEFINED_
template<typename tValue, typename tIndex, typename tFunc, typename tReduceFunc>
tValue for_range_reduce(tIndex first_1, tIndex last_1,
                        tIndex first_2, tIndex last_2,
                        tValue init, tFunc func, tReduceFunc reduce_func) {
    std::vector<tValue> partial(get_max_num_threads(), init);
    for_range(first_1, last_1,
              first_2, last_2, [&](tIndex index_1, tIndex index_2) {
        partial[get_thread_index()] =
            reduce_func(partial[get_thread_index()], func(index_1, index_2));
    });
    init = std::accumulate(
        partial.begin(), partial.end() - 1, partial.back(), reduce_func);
    return init;
}
#endif
#undef FOR_RANGE_REDUCE_2_DEFINED_
#ifndef FOR_RANGE_REDUCE_3_DEFINED_
template<typename tValue, typename tIndex, typename tFunc, typename tReduceFunc>
tValue for_range_reduce(tIndex first_1, tIndex last_1,
                        tIndex first_2, tIndex last_2,
                        tIndex first_3, tIndex last_3,
                        tValue init, tFunc func, tReduceFunc reduce_func) {
    std::vector<tValue> partial(get_max_num_threads(), init);
    for_range(first_1, last_1, first_2, last_2,
              first_3, last_3, [&](tIndex index_1, tIndex index_2, tIndex index_3) {
        partial[get_thread_index()] =
            reduce_func(partial[get_thread_index()], func(index_1, index_2, index_3));
    });
    init = std::accumulate(
        partial.begin(), partial.end() - 1, partial.back(), reduce_func);
    return init;
}
#endif
#undef FOR_RANGE_REDUCE_3_DEFINED_

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

#ifndef FOR_RANGE_SUM_1_DEFINED_
template<typename tValue, typename tIndex, typename tFunc>
tValue for_range_sum(tIndex first, tIndex last, tValue init, tFunc func) {
    return for_range_reduce(first, last, init, func, std::plus<tValue>());
}
#endif
#undef FOR_RANGE_SUM_1_DEFINED_
#ifndef FOR_RANGE_SUM_2_DEFINED_
template<typename tValue, typename tIndex, typename tFunc>
tValue for_range_sum(tIndex first_1, tIndex last_1,
                     tIndex first_2, tIndex last_2, tValue init, tFunc func) {
    return for_range_reduce(first_1, last_1,
                            first_2, last_2, init, func, std::plus<tValue>());
}
#endif
#undef FOR_RANGE_SUM_2_DEFINED_
#ifndef FOR_RANGE_SUM_3_DEFINED_
template<typename tValue, typename tIndex, typename tFunc>
tValue for_range_sum(tIndex first_1, tIndex last_1,
                     tIndex first_2, tIndex last_2,
                     tIndex first_3, tIndex last_3, tValue init, tFunc func) {
    return for_range_reduce(first_1, last_1, first_2, last_2,
                            first_3, last_3, init, func, std::plus<tValue>());
}
#endif
#undef FOR_RANGE_SUM_3_DEFINED_

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

#ifndef FOR_RANGE_MIN_1_DEFINED_
template<typename tValue, typename tIndex, typename tFunc>
tValue for_range_min(tIndex first, tIndex last, tValue init, tFunc func) {
    return for_range_reduce(first, last, init, func, tMinFunc<tValue>());
}
#endif
#undef FOR_RANGE_MIN_1_DEFINED_
#ifndef FOR_RANGE_MIN_2_DEFINED_
template<typename tValue, typename tIndex, typename tFunc>
tValue for_range_min(tIndex first_1, tIndex last_1,
                     tIndex first_2, tIndex last_2, tValue init, tFunc func) {
    return for_range_reduce(first_1, last_1,
                            first_2, last_2, init, func, tMinFunc<tValue>());
}
#endif
#undef FOR_RANGE_MIN_2_DEFINED_
#ifndef FOR_RANGE_MIN_3_DEFINED_
template<typename tValue, typename tIndex, typename tFunc>
tValue for_range_min(tIndex first_1, tIndex last_1,
                     tIndex first_2, tIndex last_2,
                     tIndex first_3, tIndex last_3, tValue init, tFunc func) {
    return for_range_reduce(first_1, last_1, first_2, last_2,
                            first_3, last_3, init, func, tMinFunc<tValue>());
}
#endif
#undef FOR_RANGE_MIN_3_DEFINED_

#ifndef FOR_RANGE_MAX_1_DEFINED_
template<typename tValue, typename tIndex, typename tFunc>
tValue for_range_max(tIndex first, tIndex last, tValue init, tFunc func) {
    return for_range_reduce(first, last, init, func, tMaxFunc<tValue>());
}
#endif
#undef FOR_RANGE_MAX_1_DEFINED_
#ifndef FOR_RANGE_MAX_2_DEFINED_
template<typename tValue, typename tIndex, typename tFunc>
tValue for_range_max(tIndex first_1, tIndex last_1,
                     tIndex first_2, tIndex last_2, tValue init, tFunc func) {
    return for_range_reduce(first_1, last_1,
                            first_2, last_2, init, func, tMaxFunc<tValue>());
}
#endif
#undef FOR_RANGE_MAX_2_DEFINED_
#ifndef FOR_RANGE_MAX_3_DEFINED_
template<typename tValue, typename tIndex, typename tFunc>
tValue for_range_max(tIndex first_1, tIndex last_1,
                     tIndex first_2, tIndex last_2,
                     tIndex first_3, tIndex last_3, tValue init, tFunc func) {
    return for_range_reduce(first_1, last_1, first_2, last_2,
                            first_3, last_3, init, func, tMaxFunc<tValue>());
}
#endif
#undef FOR_RANGE_MAX_3_DEFINED_

#ifndef FOR_RANGE_MINMAX_1_DEFINED_
template<typename tValue, typename tIndex, typename tFunc>
std::pair<tValue, tValue> for_range_minmax(tIndex first, tIndex last,
                                           tValue min_init, tValue max_init, tFunc func) {
    return for_range_reduce(first, last,
                            std::make_pair(min_init, max_init), func, tMinMaxFunc<tValue>());
}
#endif
#undef FOR_RANGE_MINMAX_1_DEFINED_
#ifndef FOR_RANGE_MINMAX_2_DEFINED_
template<typename tValue, typename tIndex, typename tFunc>
std::pair<tValue, tValue> for_range_minmax(tIndex first_1, tIndex last_1,
                                           tIndex first_2, tIndex last_2,
                                           tValue min_init, tValue max_init, tFunc func) {
    return for_range_reduce(first_1, last_1, first_2, last_2,
                            std::make_pair(min_init, max_init), func, tMinMaxFunc<tValue>());
}
#endif
#undef FOR_RANGE_MINMAX_2_DEFINED_
#ifndef FOR_RANGE_MINMAX_3_DEFINED_
template<typename tValue, typename tIndex, typename tFunc>
std::pair<tValue, tValue> for_range_minmax(tIndex first_1, tIndex last_1,
                                           tIndex first_2, tIndex last_2,
                                           tIndex first_3, tIndex last_3,
                                           tValue min_init, tValue max_init, tFunc func) {
    return for_range_reduce(first_1, last_1,
                            first_2, last_2, first_3, last_3,
                            std::make_pair(min_init, max_init), func, tMinMaxFunc<tValue>());
}
#endif
#undef FOR_RANGE_MINMAX_3_DEFINED_

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif  // ifndef PARALLEL_HH_
