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

#ifndef FEATHERS_HAS_TBB
#error Thread Building Blocks should be enabled.
#endif

#include <tbb/tbb.h>

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

extern void set_max_num_threads_impl_(uint_t num_threads);

uint_t get_thread_index() {
    FEATHERS_NOT_IMPLEMENTED();
}
uint_t get_max_num_threads() {
    return tbb::global_control::active_value(
        tbb::global_control::max_allowed_parallelism);
}
static void set_max_num_threads(uint_t num_threads) {
    set_max_num_threads_impl_(num_threads);
}
#define THREAD_IDS_DEFINED_

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

template<typename tIndex, typename tFunc>
void for_range(tIndex first, tIndex last, tFunc func) {
    tbb::parallel_for(
        tbb::blocked_range<tIndex>(first, last),
        [&](const tbb::blocked_range<tIndex>& range) {
            for (tIndex index = range.begin(); index < range.end(); ++index) {
                func(index);
            }
        });
}
#define FOR_RANGE_1_DEFINED_
template<typename tIndex, typename tFunc>
void for_range(tIndex first_1, tIndex last_1,
               tIndex first_2, tIndex last_2, tFunc func) {
    tbb::parallel_for(
        tbb::blocked_range2d<tIndex>(first_1, last_1, first_2, last_2),
        [&](const tbb::blocked_range2d<tIndex>& range) {
            const auto& rows = range.rows(), & cols = range.cols();
            for (tIndex index_1 = rows.begin(); index_1 < rows.end(); ++index_1) {
                for (tIndex index_2 = cols.begin(); index_2 < cols.end(); ++index_2) {
                    func(index_1, index_2);
                }
            }
        });
}
#define FOR_RANGE_2_DEFINED_
template<typename tIndex, typename tFunc>
void for_range(tIndex first_1, tIndex last_1,
               tIndex first_2, tIndex last_2,
               tIndex first_3, tIndex last_3, tFunc func) {
    tbb::parallel_for(
        tbb::blocked_range3d<tIndex>(first_1, last_1, first_2, last_2, first_3, last_3),
        [&](const tbb::blocked_range3d<tIndex>& range) {
            const auto& pages = range.pages();
            const auto& rows = range.rows(), &cols = range.cols();
            for (tIndex index_1 = pages.begin(); index_1 < pages.end(); ++index_1) {
                for (tIndex index_2 = rows.begin(); index_2 < rows.end(); ++index_2) {
                    for (tIndex index_3 = cols.begin(); index_3 < cols.end(); ++index_3) {
                        func(index_1, index_2, index_3);
                    }
                }
            }
        });
}
#define FOR_RANGE_3_DEFINED_

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

template<typename tValue, typename tIndex, typename tFunc, typename tReduceFunc>
tValue for_range_reduce(tIndex first, tIndex last,
                        tValue init, tFunc func, tReduceFunc reduce_func) {
    return tbb::parallel_reduce(
        tbb::blocked_range<tIndex>(first, last),
        init,
        [&](const tbb::blocked_range<tIndex>& range, tValue current) {
            for (tIndex index = range.begin(); index < range.end(); ++index) {
                current = reduce_func(current, func(index));
            }
            return current;
        },
        reduce_func);
}
#define FOR_RANGE_REDUCE_1_DEFINED_
template<typename tValue, typename tIndex, typename tFunc, typename tReduceFunc>
tValue for_range_reduce(tIndex first_1, tIndex last_1,
                        tIndex first_2, tIndex last_2,
                        tValue init, tFunc func, tReduceFunc reduce_func) {
    return tbb::parallel_reduce(
        tbb::blocked_range2d<tIndex>(first_1, last_1, first_2, last_2),
        init,
        [&](const tbb::blocked_range2d<tIndex>& range, tValue current) {
            const auto& rows = range.rows(), &cols = range.cols();
            for (tIndex index_1 = rows.begin(); index_1 < rows.end(); ++index_1) {
                for (tIndex index_2 = cols.begin(); index_2 < cols.end(); ++index_2) {
                    current = reduce_func(current, func(index_1, index_2));
                }
            }
            return current;
        },
        reduce_func);
}
#define FOR_RANGE_REDUCE_2_DEFINED_
template<typename tValue, typename tIndex, typename tFunc, typename tReduceFunc>
tValue for_range_reduce(tIndex first_1, tIndex last_1,
                        tIndex first_2, tIndex last_2,
                        tIndex first_3, tIndex last_3,
                        tValue init, tFunc func, tReduceFunc reduce_func) {
    return tbb::parallel_reduce(
        tbb::blocked_range3d<tIndex>(first_1, last_1, first_2, last_2, first_3, last_3),
        init,
        [&](const tbb::blocked_range3d<tIndex>& range, tValue current) {
            const auto& pages = range.pages();
            const auto& rows = range.rows(), &cols = range.cols();
            for (tIndex index_1 = pages.begin(); index_1 < pages.end(); ++index_1) {
                for (tIndex index_2 = rows.begin(); index_2 < rows.end(); ++index_2) {
                    for (tIndex index_3 = cols.begin(); index_3 < cols.end(); ++index_3) {
                        current = reduce_func(current, func(index_1, index_2, index_3));
                    }
                }
            }
            return current;
        },
        reduce_func);
}
#define FOR_RANGE_REDUCE_3_DEFINED_

} // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
