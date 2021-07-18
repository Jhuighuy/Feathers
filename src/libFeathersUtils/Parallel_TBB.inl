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

extern void set_num_threads_tbb_(uint_t num_threads);

uint_t get_thread_index() {
    FEATHERS_NOT_IMPLEMENTED();
}
uint_t get_max_num_threads() {
    return tbb::global_control::active_value(
        tbb::global_control::max_allowed_parallelism);
}
static void set_max_num_threads(uint_t num_threads) {
    set_num_threads_tbb_(num_threads);
}
#define THREAD_IDS_DEFINED_

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

template<typename tIter, typename tFunc>
void for_range(tIter first, tIter last, tFunc func) {
    tbb::parallel_for(
        tbb::blocked_range<tIter>(first, last),
        [&](const tbb::blocked_range<tIter>& range) {
            for (tIter iter = range.begin(); iter < range.end(); ++iter) {
                func(iter);
            }
        });
}
#define FOR_RANGE_1_DEFINED_
template<typename tIter, typename tFunc>
void for_range(tIter first_1, tIter last_1,
               tIter first_2, tIter last_2, tFunc func) {
    tbb::parallel_for(
        tbb::blocked_range2d<tIter>(first_1, last_1, first_2, last_2),
        [&](const tbb::blocked_range2d<tIter>& range) {
            for (tIter iter_1 = range.rows().begin(); iter_1 < range.rows().end(); ++iter_1) {
                for (tIter iter_2 = range.cols().begin(); iter_2 < range.cols().end(); ++iter_2) {
                    func(iter_1, iter_2);
                }
            }
        });
}
#define FOR_RANGE_2_DEFINED_
template<typename tIter, typename tFunc>
void for_range(tIter first_1, tIter last_1,
               tIter first_2, tIter last_2,
               tIter first_3, tIter last_3, tFunc func) {
    tbb::parallel_for(
        tbb::blocked_range3d<tIter>(first_1, last_1, first_2, last_2, first_3, last_3),
        [&](const tbb::blocked_range3d<tIter>& range) {
            for (tIter iter_1 = range.pages().begin(); iter_1 < range.pages().end(); ++iter_1) {
                for (tIter iter_2 = range.rows().begin(); iter_2 < range.rows().end(); ++iter_2) {
                    for (tIter iter_3 = range.cols().begin(); iter_3 < range.cols().end(); ++iter_3) {
                        func(iter_1, iter_2, iter_3);
                    }
                }
            }
        });
}
#define FOR_RANGE_3_DEFINED_

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

template<typename tValue, typename tIter, typename tFunc, typename tReduceFunc>
tValue for_range_reduce(tIter first, tIter last,
                        tValue init, tFunc func, tReduceFunc reduce_func) {
    return tbb::parallel_reduce(
        tbb::blocked_range<tIter>(first, last),
        init,
        [&](const tbb::blocked_range<tIter>& range, tValue current) {
            for (tIter iter = range.begin(); iter < range.end(); ++iter) {
                current = reduce_func(current, func(iter));
            }
            return current;
        },
        reduce_func);
}
#define FOR_RANGE_REDUCE_1_DEFINED_
template<typename tValue, typename tIter, typename tFunc, typename tReduceFunc>
tValue for_range_reduce(tIter first_1, tIter last_1,
                        tIter first_2, tIter last_2,
                        tValue init, tFunc func, tReduceFunc reduce_func) {
    return tbb::parallel_reduce(
        tbb::blocked_range2d<tIter>(first_1, last_1, first_2, last_2),
        init,
        [&](const tbb::blocked_range2d<tIter>& range, tValue current) {
            for (tIter iter_1 = range.rows().begin(); iter_1 < range.rows().end(); ++iter_1) {
                for (tIter iter_2 = range.cols().begin(); iter_2 < range.cols().end(); ++iter_2) {
                    current = reduce_func(current, func(iter_1, iter_2));
                }
            }
            return current;
        },
        reduce_func);
}
#define FOR_RANGE_REDUCE_2_DEFINED_
template<typename tValue, typename tIter, typename tFunc, typename tReduceFunc>
tValue for_range_reduce(tIter first_1, tIter last_1,
                        tIter first_2, tIter last_2,
                        tIter first_3, tIter last_3,
                        tValue init, tFunc func, tReduceFunc reduce_func) {
    return tbb::parallel_reduce(
        tbb::blocked_range3d<tIter>(first_1, last_1, first_2, last_2, first_3, last_3),
        init,
        [&](const tbb::blocked_range3d<tIter>& range, tValue current) {
            for (tIter iter_1 = range.pages().begin(); iter_1 < range.pages().end(); ++iter_1) {
                for (tIter iter_2 = range.rows().begin(); iter_2 < range.rows().end(); ++iter_2) {
                    for (tIter iter_3 = range.cols().begin(); iter_3 < range.cols().end(); ++iter_3) {
                        current = reduce_func(current, func(iter_1, iter_2, iter_3));
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
