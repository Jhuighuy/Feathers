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
#ifndef SKUNK_MISC_PARALLEL_HH
#define SKUNK_MISC_PARALLEL_HH

#include <SkunkBase.hh>

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

/**************************************************************************/
/**************************************************************************/

/** Get maximum number of threads. */
SKUNK_INLINE uint_t get_num_threads();
/** Set maximum number of threads.
 ** @returns Number of threads set. */
SKUNK_INLINE uint_t set_num_threads(uint_t num_threads);

/** Run a range in parallel. */
/** @{ */
template<typename TIter, typename TFunc>
SKUNK_INLINE TFunc for_range(TIter iter, TIter iter_end, TFunc func);
template<typename TIter, typename func_t>
SKUNK_INLINE func_t for_range(TIter x_iter, TIter x_iter_end,
                              TIter y_iter, TIter y_iter_end, func_t func);
template<typename TIter, typename func_t>
SKUNK_INLINE func_t for_range(TIter x_iter, TIter x_iter_end,
                              TIter y_iter, TIter y_iter_end,
                              TIter z_iter, TIter z_iter_end, func_t func);
/** @} */

/**************************************************************************/
/**************************************************************************/

/** Compute sum of the range values in parallel. */
/** @{ */
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE type_t for_range_sum(iter_t iter_beg, iter_t iter_end, type_t init, func_t func);
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE type_t for_range_sum(iter_t x_iter_beg, iter_t x_iter_end,
                                  iter_t y_iter_beg, iter_t y_iter_end, type_t init, func_t func);
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE type_t for_range_sum(iter_t x_iter_beg, iter_t x_iter_end,
                                  iter_t y_iter_beg, iter_t y_iter_end,
                                  iter_t z_iter_beg, iter_t z_iter_end, type_t init, func_t func);
/** @} */

/** Compute minimum of the range values in parallel. */
/** @{ */
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE type_t for_range_min(iter_t iter_beg, iter_t iter_end, type_t init, func_t func);
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE type_t for_range_min(iter_t x_iter_beg, iter_t x_iter_end,
                                  iter_t y_iter_beg, iter_t y_iter_end, type_t init, func_t func);
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE type_t for_range_min(iter_t x_iter_beg, iter_t x_iter_end,
                                  iter_t y_iter_beg, iter_t y_iter_end,
                                  iter_t z_iter_beg, iter_t z_iter_end, type_t init, func_t func);
/** @} */

/** Compute maximum of the range values in parallel. */
/** @{ */
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE type_t for_range_max(iter_t iter_beg, iter_t iter_end, type_t init, func_t func);
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE type_t for_range_max(iter_t x_iter_beg, iter_t x_iter_end,
                                  iter_t y_iter_beg, iter_t y_iter_end, type_t init, func_t func);
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE type_t for_range_max(iter_t x_iter_beg, iter_t x_iter_end,
                                  iter_t y_iter_beg, iter_t y_iter_end,
                                  iter_t z_iter_beg, iter_t z_iter_end, type_t init, func_t func);
/** @} */

/** Compute minimum and maximum of the range values in parallel. */
/** @{ */
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE std::pair<type_t, type_t> for_range_minmax(iter_t iter_beg, iter_t iter_end,
                                                        type_t min_init, type_t max_init, func_t func);
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE std::pair<type_t, type_t> for_range_minmax(iter_t x_iter_beg, iter_t x_iter_end,
                                                        iter_t y_iter_beg, iter_t y_iter_end,
                                                        type_t min_init, type_t max_init, func_t func);
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE std::pair<type_t, type_t> for_range_minmax(iter_t x_iter_beg, iter_t x_iter_end,
                                                        iter_t y_iter_beg, iter_t y_iter_end,
                                                        iter_t z_iter_beg, iter_t z_iter_end,
                                                        type_t min_init, type_t max_init, func_t func);
/** @} */

/** Reduce a range in parallel. */
/** @{ */
template<typename iter_t, typename type_t, typename reduce_t, typename func_t>
SKUNK_INLINE type_t for_range_reduce(iter_t iter, iter_t iter_end, type_t init, reduce_t reduce, func_t func);
template<typename iter_t, typename type_t, typename reduce_t, typename func_t>
SKUNK_INLINE type_t for_range_reduce(iter_t x_iter, iter_t x_iter_end,
                                     iter_t y_iter, iter_t y_iter_end, type_t init, reduce_t reduce, func_t func);
template<typename iter_t, typename type_t, typename reduce_t, typename func_t>
SKUNK_INLINE type_t for_range_reduce(iter_t x_iter, iter_t x_iter_end,
                                     iter_t y_iter, iter_t y_iter_end,
                                     iter_t z_iter, iter_t z_iter_end, type_t init, reduce_t reduce, func_t func);
/** @} */

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#if SKUNK_HAS_OPENMP2
#include <omp.h>
#endif

namespace skunk {

#if SKUNK_HAS_OPENMP2

/** Get maximum number of threads. */
SKUNK_INLINE uint_t get_num_threads() {
    return static_cast<uint_t>(omp_get_max_threads());
}   // get_num_threads
/** Set maximum number of threads.
 ** @returns Number of threads set. */
SKUNK_INLINE uint_t set_num_threads(uint_t num_threads) {
    omp_set_num_threads(static_cast<int>(num_threads));
    return num_threads;
}   // set_num_threads

/**************************************************************************/
/**************************************************************************/

/** Run a range in parallel. */
/** @{ */
template<typename iter_t, typename func_t>
SKUNK_INLINE func_t for_range(iter_t iter_beg, iter_t iter_end, func_t func) {
#pragma omp parallel for default(none) shared(iter_beg, iter_end, func)
    for (iter_t iter = iter_beg; iter < iter_end; ++iter) {
        func(iter);
    }
    return func;
}   // for_range
#if SKUNK_HAS_OPENMP3
template<typename iter_t, typename func_t>
SKUNK_INLINE func_t for_range(iter_t x_iter_beg, iter_t x_iter_end,
                              iter_t y_iter_beg, iter_t y_iter_end, func_t func) {
#pragma omp parallel for collapse(2)
    for (iter_t x_iter = x_iter_beg; x_iter < x_iter_end; ++x_iter) {
        for (iter_t y_iter = y_iter_beg; y_iter < y_iter_end; ++y_iter) {
            func(x_iter, y_iter);
        }
    }
    return func;
}   // for_range
template<typename iter_t, typename func_t>
SKUNK_INLINE func_t for_range(iter_t x_iter_beg, iter_t x_iter_end,
                              iter_t y_iter_beg, iter_t y_iter_end,
                              iter_t z_iter_beg, iter_t z_iter_end, func_t func) {
#pragma omp parallel for collapse(3)
    for (iter_t x_iter = x_iter_beg; x_iter < x_iter_end; ++x_iter) {
        for (iter_t y_iter = y_iter_beg; y_iter < y_iter_end; ++y_iter) {
            for (iter_t z_iter = z_iter_beg; z_iter < z_iter_end; ++z_iter) {
                func(x_iter, y_iter, z_iter);
            }
        }
    }
    return func;
}   // for_range
#endif  // if SKUNK_HAS_OPENMP3
/** @} */

/**************************************************************************/
/**************************************************************************/

/** Compute sum of the range values in parallel. */
/** @{ */
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE type_t for_range_sum(iter_t iter_beg, iter_t iter_end, type_t init, func_t func) {
    type_t sum{};
#pragma omp parallel for reduction(+:sum)
    for (iter_t iter = iter_beg; iter < iter_end; ++iter) {
        sum += func(iter);
    }
    return init + sum;
}   // for_range_sum
#if SKUNK_HAS_OPENMP3
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE type_t for_range_sum(iter_t x_iter_beg, iter_t x_iter_end,
                                  iter_t y_iter_beg, iter_t y_iter_end, type_t init, func_t func) {
    type_t sum{};
#pragma omp parallel for collapse(2) reduction(+:sum)
    for (iter_t x_iter = x_iter_beg; x_iter < x_iter_end; ++x_iter) {
        for (iter_t y_iter = y_iter_beg; y_iter < y_iter_end; ++y_iter) {
            sum += func(x_iter, y_iter);
        }
    }
    return init + sum;
}   // for_range_sum
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE type_t for_range_sum(iter_t x_iter_beg, iter_t x_iter_end,
                                  iter_t y_iter_beg, iter_t y_iter_end,
                                  iter_t z_iter_beg, iter_t z_iter_end, type_t init, func_t func) {
    type_t sum{};
#pragma omp parallel for collapse(3) reduction(+:sum)
    for (iter_t x_iter = x_iter_beg; x_iter < x_iter_end; ++x_iter) {
        for (iter_t y_iter = y_iter_beg; y_iter < y_iter_end; ++y_iter) {
            for (iter_t z_iter = z_iter_beg; z_iter < z_iter_end; ++z_iter) {
                sum += func(x_iter, y_iter, z_iter);
            }
        }
    }
    return init + sum;
}   // for_range_sum
#endif  // if SKUNK_HAS_OPENMP3
/** @} */

#if SKUNK_HAS_OPENMP4
/** Compute minimum of the range values in parallel. */
/** @{ */
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE type_t for_range_min(iter_t iter_beg, iter_t iter_end, type_t init, func_t func) {
    type_t min_val{};
#pragma omp parallel for reduction(min:min_val)
    for (iter_t iter = iter_beg; iter < iter_end; ++iter) {
        min_val = std::min(min_val, func(iter));
    }
    return std::min(min_val, init);
}   // for_range_min
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE type_t for_range_min(iter_t x_iter_beg, iter_t x_iter_end,
                                  iter_t y_iter_beg, iter_t y_iter_end, type_t init, func_t func) {
    type_t min_val{};
#pragma omp parallel for collapse(2) reduction(min:min_val)
    for (iter_t x_iter = x_iter_beg; x_iter < x_iter_end; ++x_iter) {
        for (iter_t y_iter = y_iter_beg; y_iter < y_iter_end; ++y_iter) {
            min_val = std::min(min_val, func(x_iter, y_iter));
        }
    }
    return std::min(min_val, init);
}   // for_range_min
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE type_t for_range_min(iter_t x_iter_beg, iter_t x_iter_end,
                                  iter_t y_iter_beg, iter_t y_iter_end,
                                  iter_t z_iter_beg, iter_t z_iter_end, type_t init, func_t func) {
    type_t min_val{};
#pragma omp parallel for collapse(3) reduction(min:min_val)
    for (iter_t x_iter = x_iter_beg; x_iter < x_iter_end; ++x_iter) {
        for (iter_t y_iter = y_iter_beg; y_iter < y_iter_end; ++y_iter) {
            for (iter_t z_iter = z_iter_beg; z_iter < z_iter_end; ++z_iter) {
                min_val = std::min(min_val, func(x_iter, y_iter, z_iter));
            }
        }
    }
    return std::min(min_val, init);
}   // for_range_min
/** @} */
#endif  // if SKUNK_HAS_OPENMP4

#if SKUNK_HAS_OPENMP4
/** Compute maximum of the range values in parallel. */
/** @{ */
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE type_t for_range_max(iter_t iter_beg, iter_t iter_end, type_t init, func_t func) {
    type_t max_val{};
#pragma omp parallel for reduction(max:max_val)
    for (iter_t iter = iter_beg; iter < iter_end; ++iter) {
        max_val = std::max(max_val, func(iter));
    }
    return std::max(max_val, init);
}   // for_range_max
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE type_t for_range_max(iter_t x_iter_beg, iter_t x_iter_end,
                                  iter_t y_iter_beg, iter_t y_iter_end, type_t init, func_t func) {
    type_t max_val{};
#pragma omp parallel for collapse(2) reduction(max:max_val)
    for (iter_t x_iter = x_iter_beg; x_iter < x_iter_end; ++x_iter) {
        for (iter_t y_iter = y_iter_beg; y_iter < y_iter_end; ++y_iter) {
            max_val = std::max(max_val, func(x_iter, y_iter));
        }
    }
    return std::max(max_val, init);
}   // for_range_max
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE type_t for_range_max(iter_t x_iter_beg, iter_t x_iter_end,
                                  iter_t y_iter_beg, iter_t y_iter_end,
                                  iter_t z_iter_beg, iter_t z_iter_end, type_t init, func_t func) {
    type_t max_val{};
#pragma omp parallel for collapse(3) reduction(max:max_val)
    for (iter_t x_iter = x_iter_beg; x_iter < x_iter_end; ++x_iter) {
        for (iter_t y_iter = y_iter_beg; y_iter < y_iter_end; ++y_iter) {
            for (iter_t z_iter = z_iter_beg; z_iter < z_iter_end; ++z_iter) {
                max_val = std::max(max_val, func(x_iter, y_iter, z_iter));
            }
        }
    }
    return std::max(max_val, init);
}   // for_range_max
/** @} */
#endif  // if SKUNK_HAS_OPENMP4

#if SKUNK_HAS_OPENMP4
/** Compute minimum and maximum of the range values in parallel. */
/** @{ */
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE std::pair<type_t, type_t> for_range_minmax(iter_t iter_beg, iter_t iter_end,
                                                        type_t min_init, type_t max_init, func_t func) {
    type_t min_val{}, max_val{};
#pragma omp parallel for reduction(min:min_val) reduction(max:max_val)
    for (iter_t iter = iter_beg; iter < iter_end; ++iter) {
        min_val = std::min(min_val, func(iter));
        max_val = std::max(max_val, func(iter));
    }
    return {std::min(min_val, min_init), std::max(max_val, max_init)};
}   // for_range_minmax
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE std::pair<type_t, type_t> for_range_minmax(iter_t x_iter_beg, iter_t x_iter_end,
                                                        iter_t y_iter_beg, iter_t y_iter_end,
                                                        type_t min_init, type_t max_init, func_t func) {
    type_t min_val{}, max_val{};
#pragma omp parallel for collapse(2) reduction(min:min_val) reduction(max:max_val)
    for (iter_t x_iter = x_iter_beg; x_iter < x_iter_end; ++x_iter) {
        for (iter_t y_iter = y_iter_beg; y_iter < y_iter_end; ++y_iter) {
            min_val = std::min(min_val, func(x_iter, y_iter));
            max_val = std::max(max_val, func(x_iter, y_iter));
        }
    }
    return {std::min(min_val, min_init), std::max(max_val, max_init)};
}   // for_range_minmax
template<typename iter_t, typename type_t, typename func_t>
SKUNK_INLINE std::pair<type_t, type_t> for_range_minmax(iter_t x_iter_beg, iter_t x_iter_end,
                                                        iter_t y_iter_beg, iter_t y_iter_end,
                                                        iter_t z_iter_beg, iter_t z_iter_end,
                                                        type_t min_init, type_t max_init, func_t func) {
    type_t min_val{}, max_val{};
#pragma omp parallel for collapse(3) reduction(min:min_val) reduction(max:max_val)
    for (iter_t x_iter = x_iter_beg; x_iter < x_iter_end; ++x_iter) {
        for (iter_t y_iter = y_iter_beg; y_iter < y_iter_end; ++y_iter) {
            for (iter_t z_iter = z_iter_beg; z_iter < z_iter_end; ++z_iter) {
                min_val = std::min(min_val, func(x_iter, y_iter, z_iter));
                max_val = std::max(max_val, func(x_iter, y_iter, z_iter));
            }
        }
    }
    return {std::min(min_val, min_init), std::max(max_val, max_init)};
}   // for_range_minmax
/** @} */
#endif  // if SKUNK_HAS_OPENMP4

#endif  // if SKUNK_HAS_OPENMP2

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#if 0

#include <thread>
#include <functional>

namespace skunk {

/** Get maximum number of threads. */
SKUNK_INLINE size_t get_num_threads() {
    return std::max(1u, std::thread::hardware_concurrency());
}   // get_num_threads
/** Set maximum number of threads.
 ** @returns Number of threads set. */
SKUNK_INLINE size_t set_num_threads(size_t num_threads) {
    return get_num_threads();
}   // set_num_threads

/**************************************************************************/
/**************************************************************************/

/** Run a range in parallel. */
/** @{ */
template<typename iter_t, typename func_t>
SKUNK_INLINE func_t for_range(iter_t iter_beg, iter_t iter_end, func_t func) {
    const size_t num_threads = get_num_threads();
    if (num_threads <= 1) {
        for (auto iter = iter_beg; iter < iter_end; ++iter) {
            func(iter);
        }
        return func;
    }
    const size_t size = iter_end - iter_beg;
    const size_t batch_size = size/num_threads,
                 batch_remainder = size%num_threads;
    /* Launch threads. */
    std::vector<std::thread> threads(num_threads);
    for (size_t thread_ind = 0; thread_ind < num_threads; ++thread_ind) {
        const size_t beg = thread_ind*batch_size,
                     end = thread_ind*batch_size + batch_size;
        threads[thread_ind] = std::thread([&](size_t beg, size_t end) {
            for (auto iter = iter_beg + beg; iter < iter_beg + end; ++iter) {
                func(iter);
            }
        }, beg, end);
    }
    /* Compute the remainder. */
    const size_t beg = num_threads*batch_size,
                 end = num_threads*batch_size + batch_remainder;
    for (auto iter = iter_beg + beg; iter < iter_beg + end; ++iter) {
        func(iter);
    }
    /* Wait for threads. */
    std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
    return func;
}   // for_range
template<typename iter_t, typename func_t>
SKUNK_INLINE func_t for_range(iter_t x_iter_beg, iter_t x_iter_end,
                              iter_t y_iter_beg, iter_t y_iter_end, func_t func) {
    const size_t size = (x_iter_end - x_iter_beg)*(y_iter_end - y_iter_beg);
    for_range(0, size, [&](size_t batch_ind) {
        SKUNK_NOT_IMPLEMENTED();
    });
    return func;
}   // for_range
template<typename iter_t, typename func_t>
SKUNK_INLINE func_t for_range(iter_t x_iter_beg, iter_t x_iter_end,
                              iter_t y_iter_beg, iter_t y_iter_end,
                              iter_t z_iter_beg, iter_t z_iter_end, func_t func) {
    const size_t size = (x_iter_end - x_iter_beg)*(y_iter_end - y_iter_beg)*(z_iter_end - z_iter_beg);
    for_range(0, size, [&](size_t batch_ind) {
        SKUNK_NOT_IMPLEMENTED();
    });
    return func;
}   // for_range
/** @} */

}   // namespace skunk

#endif

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif  // ifndef SKUNK_MISC_PARALLEL_HH
