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

#ifndef FEATHERS_HAS_OPENMP
#error OpenMP should be enabled.
#endif

/** OpenMP 2.0 support. */
/** @{ */
#if (_OPENMP >= 200203)
#define FEATHERS_HAS_OPENMP_2_0 1
#else
#define FEATHERS_HAS_OPENMP_2_0 0
#endif
/** @} */

/** OpenMP 2.5 support. */
/** @{ */
#if (!FEATHERS_CONFIG_FORCE_OPENMP_2_0 && (_OPENMP >= 200505))
#define FEATHERS_HAS_OPENMP_2_5 1
#else
#define FEATHERS_HAS_OPENMP_2_5 0
#endif
/** @} */

/** OpenMP 3.0 support. */
/** @{ */
#if (!FEATHERS_CONFIG_FORCE_OPENMP_2_0 && (_OPENMP >= 200805))
#define FEATHERS_HAS_OPENMP_3_0 1
#else
#define FEATHERS_HAS_OPENMP_3_0 0
#endif
/** @} */

/** OpenMP 3.0 support. */
/** @{ */
#if (!FEATHERS_CONFIG_FORCE_OPENMP_2_0 && (_OPENMP >= 201107))
#define FEATHERS_HAS_OPENMP_3_1 1
#else
#define FEATHERS_HAS_OPENMP_3_1 0
#endif
/** @} */

/** OpenMP 4.0 support. */
/** @{ */
#if (!FEATHERS_CONFIG_FORCE_OPENMP_2_0 && (_OPENMP >= 201307))
#define FEATHERS_HAS_OPENMP_4_0 1
#else
#define FEATHERS_HAS_OPENMP_4_0 0
#endif
/** @} */

/** OpenMP 4.5 support. */
/** @{ */
#if (!FEATHERS_CONFIG_FORCE_OPENMP_2_0 && (_OPENMP >= 201511))
#define FEATHERS_HAS_OPENMP_4_5 1
#else
#define FEATHERS_HAS_OPENMP_4_5 0
#endif
/** @} */

/** OpenMP 5.0 support. */
/** @{ */
#if (!FEATHERS_CONFIG_FORCE_OPENMP_2_0 && (_OPENMP >= 201811))
#define FEATHERS_HAS_OPENMP_5_0 1
#else
#define FEATHERS_HAS_OPENMP_5_0 0
#endif
/** @} */

/** OpenMP 5.1 support. */
/** @{ */
#if (!FEATHERS_CONFIG_FORCE_OPENMP_2_0 && (_OPENMP >= 202011))
#define FEATHERS_HAS_OPENMP_5_1 1
#else
#define FEATHERS_HAS_OPENMP_5_1 0
#endif
/** @} */

#if !FEATHERS_HAS_OPENMP_2_0
#error Invalid OpenMP version.
#endif

#include <omp.h>

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

uint_t get_thread_index() {
    return static_cast<uint_t>(omp_get_thread_num());
}
uint_t get_max_num_threads() {
    return static_cast<uint_t>(omp_get_max_threads());
}
void set_max_num_threads(uint_t num_threads) {
    omp_set_num_threads(static_cast<int>(num_threads));
}
#define THREAD_IDS_DEFINED_

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

#if FEATHERS_HAS_OPENMP_3_0
template<typename tIter, typename tFunc>
void for_range(tIter first, tIter last, tFunc func) {
#pragma omp parallel for schedule(static)
    for (tIter iter = first; iter < last; ++iter) {
        func(iter);
    }
}
#define FOR_RANGE_1_DEFINED_
#elif FEATHERS_HAS_OPENMP_2_0
template<typename tIter, typename tFunc>
void for_range(tIter first, tIter last, tFunc func) {
    const ptrdiff_t count = last - first;
#pragma omp parallel for schedule(static)
    for (ptrdiff_t index = 0; index < count; ++index) {
        func(first + index);
    }
}
#define FOR_RANGE_1_DEFINED_
#endif
#if FEATHERS_HAS_OPENMP_3_0
template<typename tIter, typename tFunc>
void for_range(tIter first_1, tIter last_1,
               tIter first_2, tIter last_2, tFunc func) {
#pragma omp parallel for collapse(2) schedule(static)
    for (tIter iter_1 = first_1; iter_1 < last_1; ++iter_1) {
        for (tIter iter_2 = first_2; iter_2 < last_2; ++iter_2) {
            func(iter_1, iter_2);
        }
    }
}
#define FOR_RANGE_2_DEFINED_
template<typename tIter, typename tFunc>
void for_range(tIter first_1, tIter last_1,
               tIter first_2, tIter last_2,
               tIter first_3, tIter last_3, tFunc func) {
#pragma omp parallel for collapse(3) schedule(static)
    for (tIter iter_1 = first_1; iter_1 < last_1; ++iter_1) {
        for (tIter iter_2 = first_2; iter_2 < last_2; ++iter_2) {
            for (tIter iter_3 = first_3; iter_3 < last_3; ++iter_3) {
                func(iter_1, iter_2, iter_3);
            }
        }
    }
}
#define FOR_RANGE_3_DEFINED_
#endif // FEATHERS_HAS_OPENMP_3_0

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/* Generic reduction functionality is missing. */

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

#if FEATHERS_HAS_OPENMP_3_0
template<typename tValue, typename tIter, typename tFunc>
tValue for_range_sum(tIter first, tIter last, tValue init, tFunc func) {
#pragma omp parallel for reduction(+:init) schedule(static)
    for (tIter iter = first; iter < last; ++iter) {
        init += func(iter);
    }
    return init;
}
#define FOR_RANGE_SUM_1_DEFINED_
#elif FEATHERS_HAS_OPENMP_2_0
template<typename tValue, typename tIter, typename tFunc>
tValue for_range_sum(tIter first, tIter last, tValue init, tFunc func) {
    const ptrdiff_t count = last - first;
#pragma omp parallel for reduction(+:init) schedule(static)
    for (ptrdiff_t index = 0; index < count; ++index) {
        init += func(first + index);
    }
    return init;
}
#define FOR_RANGE_SUM_1_DEFINED_
#endif
#if FEATHERS_HAS_OPENMP_3_0
template<typename tValue, typename tIter, typename tFunc>
tValue for_range_sum(tIter first_1, tIter last_1,
                     tIter first_2, tIter last_2, tValue init, tFunc func) {
#pragma omp parallel for collapse(2) reduction(+:init) schedule(static)
    for (tIter iter_1 = first_1; iter_1 < last_1; ++iter_1) {
        for (tIter iter_2 = first_2; iter_2 < last_2; ++iter_2) {
            init += func(iter_1, iter_2);
        }
    }
    return init;
}
#define FOR_RANGE_SUM_2_DEFINED_
template<typename tValue, typename tIter, typename tFunc>
tValue for_range_sum(tIter first_1, tIter last_1,
                     tIter first_2, tIter last_2,
                     tIter first_3, tIter last_3, tValue init, tFunc func) {
#pragma omp parallel for collapse(3) reduction(+:init) schedule(static)
    for (tIter iter_1 = first_1; iter_1 < last_1; ++iter_1) {
        for (tIter iter_2 = first_2; iter_2 < last_2; ++iter_2) {
            for (tIter iter_3 = first_3; iter_3 < last_3; ++iter_3) {
                init += func(iter_1, iter_2, iter_3);
            }
        }
    }
    return init;
}
#define FOR_RANGE_SUM_3_DEFINED_
#endif // FEATHERS_HAS_OPENMP_3_0

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

#if FEATHERS_HAS_OPENMP_4_0
template<typename tValue, typename tIter, typename tFunc>
tValue for_range_min(tIter first, tIter last, tValue init, tFunc func) {
#pragma omp parallel for reduction(min:init) schedule(static)
    for (tIter iter = first; iter < last; ++iter) {
        init = std::min(init, func(iter));
    }
    return init;
}
#define FOR_RANGE_MIN_1_DEFINED_
template<typename tValue, typename tIter, typename tFunc>
tValue for_range_min(tIter first_1, tIter last_1,
                     tIter first_2, tIter last_2, tValue init, tFunc func) {
#pragma omp parallel for collapse(2) reduction(min:init) schedule(static)
    for (tIter iter_1 = first_1; iter_1 < last_1; ++iter_1) {
        for (tIter iter_2 = first_2; iter_2 < last_2; ++iter_2) {
            init = std::min(init, func(iter_1, iter_2));
        }
    }
    return init;
}
#define FOR_RANGE_MIN_2_DEFINED_
template<typename tValue, typename tIter, typename tFunc>
tValue for_range_min(tIter first_1, tIter last_1,
                     tIter first_2, tIter last_2,
                     tIter first_3, tIter last_3, tValue init, tFunc func) {
#pragma omp parallel for collapse(3) reduction(min:init) schedule(static)
    for (tIter iter_1 = first_1; iter_1 < last_1; ++iter_1) {
        for (tIter iter_2 = first_2; iter_2 < last_2; ++iter_2) {
            for (tIter iter_3 = first_3; iter_3 < last_3; ++iter_3) {
                init = std::min(init, func(iter_1, iter_2, iter_3));
            }
        }
    }
    return init;
}
#define FOR_RANGE_MIN_3_DEFINED_
#endif // FEATHERS_HAS_OPENMP_4_0

#if FEATHERS_HAS_OPENMP_4_0
template<typename tValue, typename tIter, typename tFunc>
tValue for_range_max(tIter first, tIter last, tValue init, tFunc func) {
#pragma omp parallel for reduction(max:init) schedule(static)
    for (tIter iter = first; iter < last; ++iter) {
        init = std::max(init, func(iter));
    }
    return init;
}
#define FOR_RANGE_MAX_1_DEFINED_
template<typename tValue, typename tIter, typename tFunc>
tValue for_range_max(tIter first_1, tIter last_1,
                     tIter first_2, tIter last_2, tValue init, tFunc func) {
#pragma omp parallel for collapse(2) reduction(max:init) schedule(static)
    for (tIter iter_1 = first_1; iter_1 < last_1; ++iter_1) {
        for (tIter iter_2 = first_2; iter_2 < last_2; ++iter_2) {
            init = std::max(init, func(iter_1, iter_2));
        }
    }
    return init;
}
#define FOR_RANGE_MAX_2_DEFINED_
template<typename tValue, typename tIter, typename tFunc>
tValue for_range_max(tIter first_1, tIter last_1,
                     tIter first_2, tIter last_2,
                     tIter first_3, tIter last_3, tValue init, tFunc func) {
#pragma omp parallel for collapse(3) reduction(max:init) schedule(static)
    for (tIter iter_1 = first_1; iter_1 < last_1; ++iter_1) {
        for (tIter iter_2 = first_2; iter_2 < last_2; ++iter_2) {
            for (tIter iter_3 = first_3; iter_3 < last_3; ++iter_3) {
                init = std::max(init, func(iter_1, iter_2, iter_3));
            }
        }
    }
    return init;
}
#define FOR_RANGE_MAX_3_DEFINED_
#endif // FEATHERS_HAS_OPENMP_4_0

#if FEATHERS_HAS_OPENMP_4_0
template<typename tValue, typename tIter, typename tFunc>
std::pair<tValue, tValue> for_range_minmax(tIter first, tIter last,
                                           tValue min_init, tValue max_init, tFunc func) {
#pragma omp parallel for \
        reduction(min:min_init) reduction(max:max_init) schedule(static)
    for (tIter iter = first; iter < last; ++iter) {
        const auto current = func(iter);
        min_init = std::min(min_init, current);
        max_init = std::max(max_init, current);
    }
    return std::make_pair(min_init, max_init);
}
#define FOR_RANGE_MINMAX_1_DEFINED_
template<typename tValue, typename tIter, typename tFunc>
std::pair<tValue, tValue> for_range_minmax(tIter first_1, tIter last_1,
                                           tIter first_2, tIter last_2,
                                           tValue min_init, tValue max_init, tFunc func) {
#pragma omp parallel for collapse(2) \
        reduction(min:min_init) reduction(max:max_init) schedule(static)
    for (tIter iter_1 = first_1; iter_1 < last_1; ++iter_1) {
        for (tIter iter_2 = first_2; iter_2 < last_2; ++iter_2) {
            const auto current = func(iter_1, iter_2);
            min_init = std::min(min_init, current);
            max_init = std::max(max_init, current);
        }
    }
    return std::make_pair(min_init, max_init);
}
#define FOR_RANGE_MINMAX_2_DEFINED_
template<typename tValue, typename tIter, typename tFunc>
std::pair<tValue, tValue> for_range_minmax(tIter first_1, tIter last_1,
                                           tIter first_2, tIter last_2,
                                           tIter first_3, tIter last_3,
                                           tValue min_init, tValue max_init, tFunc func) {
#pragma omp parallel for collapse(3) \
        reduction(min:min_init) reduction(max:max_init) schedule(static)
    for (tIter iter_1 = first_1; iter_1 < last_1; ++iter_1) {
        for (tIter iter_2 = first_2; iter_2 < last_2; ++iter_2) {
            for (tIter iter_3 = first_3; iter_3 < last_3; ++iter_3) {
                const auto current = func(iter_1, iter_2, iter_3);
                min_init = std::min(min_init, current);
                max_init = std::max(max_init, current);
            }
        }
    }
    return std::make_pair(min_init, max_init);
}
#define FOR_RANGE_MINMAX_3_DEFINED_
#endif // FEATHERS_HAS_OPENMP_4_0

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
