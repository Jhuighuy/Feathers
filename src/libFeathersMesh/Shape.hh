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
#ifndef SHAPES_HH_
#define SHAPES_HH_

#include "SkunkBase.hh"

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/**
 * Type of the shapes.
 */
enum class eShape : byte_t {
    null,
    node,
    segment_2,
    triangle_3,
    quadrangle_4,
    tetrahedron_4,
    pyramid_5,
    pentahedron_6,
    hexahedron_8,
};  // enum class eShape

/**
 * Shape part description.
 */
struct sShapePart {
    eShape part_shape;
    std::vector<uint_t> part_local_node_indices;

    template<typename... tIndex>
    explicit sShapePart(eShape part_shape, tIndex... local_node_index):
        part_shape(part_shape), part_local_node_indices{local_node_index...} {
    }
};  // sShapePart

/**
 * Array of shape parts.
 */
using tShapeParts = std::vector<sShapePart>;

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Shape pointer convenience wrapper.
 */
class iShapePtr final : public std::unique_ptr<class iShape> {
public:
    /** Construct the shape. */
    explicit iShapePtr(eShape shape_type);
};  // class iShapePtr

/**
 * Abstract shape class.
 */
class iShape {
protected:
    uint_t m_num_nodes = 0;
    const vec3_t* m_node_coords = nullptr;
    std::vector<uint_t> m_node_indices;

public:
    virtual ~iShape() = default;

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get shape node position. */
    const vec3_t& get_node_coords(uint_t node_local) const {
        FEATHERS_ASSERT(node_local < m_node_indices.size());
        const uint_t node_index = m_node_indices[node_local];
        FEATHERS_ASSERT(node_index < m_num_nodes);
        return m_node_coords[node_index];
    }
    /** Assign shape nodes. */
    template<typename tIndexIter>
    void assign_nodes(uint_t num_nodes, const vec3_t* node_coords,
                      tIndexIter first_shape_node_index,
                      tIndexIter last_shape_node_index) {
        m_num_nodes = num_nodes, m_node_coords = node_coords;
        FEATHERS_ASSERT(num_shape_nodes() ==
                        (last_shape_node_index - first_shape_node_index));
        m_node_indices.assign(
            first_shape_node_index, last_shape_node_index);
    }

    template<typename... tIndex>
    sShapePart get_part(eShape part_shape, tIndex... local_node_index) const {
        return sShapePart(part_shape, m_node_indices[local_node_index]...);
    }

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get shape diameter. */
    virtual real_t get_diameter() const {
        return qnan;
    }
    /** Get shape length/area/volume. */
    virtual real_t get_length_or_area_or_volume() const {
        return qnan;
    }
    /** Get normal to shape. */
    virtual vec3_t get_normal() const {
        return vec3_t(qnan);
    }
    /** Get shape direction. */
    virtual vec3_t get_direction() const {
        return vec3_t(qnan);
    }
    /** Get shape barycenter. */
    virtual vec3_t get_center_coords() const {
        return vec3_t(qnan);
    }

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Number of nodes in the shape. */
    virtual uint_t num_shape_nodes() const = 0;
    /** Get shape edges. */
    virtual tShapeParts get_shape_edges() const = 0;
    /** Get shape faces. */
    virtual tShapeParts get_shape_faces() const = 0;
};  // class iShape

/**
 * Abstract simplex splittable shape class.
 */
class iSplittableShape : public iShape {
public:
    real_t get_diameter() const final;
    real_t get_length_or_area_or_volume() const final;
    vec3_t get_normal() const final;
    vec3_t get_center_coords() const final;

    /** Split current shape into the set of simplex shapes. */
    virtual tShapeParts simplex_split() const = 0;

private:
    template<typename tFunc>
    void for_each_part_(tFunc func) const;
};  // class iSplittableShape

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Segment shape.
 * @verbatim
 *  n0 @ f0
 *      \         e0 = (n0,n1)
 *       v e0     f0 = (n0)
 *        \       f1 = (n1)
 *      n1 @ f1
 * @endverbatim
 */
class cSegmentShape final : public iShape {
public:
    real_t get_diameter() const final;
    real_t get_length_or_area_or_volume() const final;
    vec3_t get_normal() const final;
    vec3_t get_direction() const final;
    vec3_t get_center_coords() const final;

    uint_t num_shape_nodes() const final;
    tShapeParts get_shape_edges() const final;
    tShapeParts get_shape_faces() const final;
};  // class tSegmentShape

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Triangle shape.
 * @verbatim
 *          @ n2        e0 = f0 = (n0,n1)
 *         / \          e1 = f1 = (n1,n2)
 *  e2/f2 v   ^ e1/f1   e2 = f2 = (n2,n0)
 *       /     \
 *   n0 @--->---@ n1
 *        e0/f0
 * @endverbatim
 */
class cTriangleShape final : public iShape {
public:
    real_t get_diameter() const final;
    real_t get_length_or_area_or_volume() const final;
    vec3_t get_normal() const final;
    vec3_t get_center_coords() const final;

    uint_t num_shape_nodes() const final;
    tShapeParts get_shape_edges() const final;
    tShapeParts get_shape_faces() const final;
};  // class cTriangleShape

/**
 * Quadrangle shape.
 * @verbatim
 *            e2/f2
 *       n3 @-----<-----@ n2    e0 = f0 = (n0,n1)
 *         /           /        e1 = f2 = (n1,n2)
 *  e3/f3 v           ^ e1/f1   e2 = f2 = (n2,n3)
 *       /           /          e3 = f3 = (n3,n0)
 *   n0 @----->-----@ n1     split = ((n0,n1,n2),(n2,n3,n0))
 *        e0/f0
 * @endverbatim
 */
class cQuadrangleShape final : public iSplittableShape {
public:
    uint_t num_shape_nodes() const final;
    tShapeParts get_shape_edges() const final;
    tShapeParts get_shape_faces() const final;

    tShapeParts simplex_split() const final;
};  // class cQuadrangleShape

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Tetrahedron shape.
 * @verbatim
 *                    f3
 *               n3   ^
 *                @   |
 *         f1    /|\. |     f2         e0 = (n0,n1)
 *         ^    / | `\.     ^          e1 = (n1,n2)
 *          \  /  |f1 `\.  /           e2 = (n2,n0)
 *           \/   |   | `\/            e3 = (n0,n3)
 *           /\   |   o  /`\           e4 = (n1,n3)
 *       e3 ^  *  |     *   `^.e5      e5 = (n2,n3)
 *         /   e4 ^           `\       f0 = (n0,n2,n1)
 *     n0 @-------|---------<---@ n2   f1 = (n0,n1,n3)
 *         \      |  e2       ,/       f2 = (n1,n2,n3)
 *          \     |     o   ,/`        f3 = (n2,n0,n3)
 *           \    ^ e7  | ,/`
 *         e0 v   |     ,^ e1
 *             \  |   ,/`
 *              \ | ,/` |
 *               \|/`   |
 *                @ n1  v
 *                      f0
 * @endverbatim
 */
class cTetrahedronShape final : public iShape {
public:
    real_t get_length_or_area_or_volume() const final;
    vec3_t get_center_coords() const final;

    uint_t num_shape_nodes() const final;
    tShapeParts get_shape_edges() const final;
    tShapeParts get_shape_faces() const final;
};  // class cTetrahedronShape

/**
 * Pyramid shape.
 * @verbatim
 *                                n4
 *                  f3           ,@
 *                   ^        ,/`/|\     f1
 *                    \    ,/`  / | \    ^
 *                     \,/`    /  |  \  /
 *                e7 ,/`\     /   |   \/
 *                ,^`    o   /    |   /\
 *             ,/`          /     |  o  \
 *  f4 <------------*      /   e6 ^   *--\----> f2
 *       ,/`              /       |       \
 *   n3 @-----<----------/--------@. n2    ^ e5
 *       `\.  e2        /          `\.      \
 *          `v.        ^ e4           `\. e3 \
 *          e3 `\.    /       o          `^.  \
 *                `\./        |             `\.\
 *               n0 @-------------------->-----@ n1
 *                            |          e0
 *                            |
 *                            v
 *                            fo
 * @endverbatim
 *
 * @verbatim
 * e0 = (n0,n1)
 * e1 = (n1,n2)
 * e2 = (n2,n3)
 * e3 = (n3,n0)
 * e4 = (n0,n4)
 * e5 = (n1,n4)
 * e6 = (n2,n4)
 * e7 = (n3,n4)
 * f0 = (n0,n3,n2,n1)
 * f1 = (n0,n1,n4)
 * f2 = (n1,n2,n4)
 * f3 = (n2,n3,n4)
 * f4 = (n3,n0,n4)
 * @endverbatim
 */
class cPyramidShape final : public iSplittableShape {
public:
    uint_t num_shape_nodes() const final;
    tShapeParts get_shape_edges() const final;
    tShapeParts get_shape_faces() const final;

    tShapeParts simplex_split() const final;
};  // class cPyramidShape

/**
 * Pentahedron (triangular prism) shape.
 * @verbatim
 *                 f4
 *                 ^  f2
 *                 |  ^
 *             e8  |  |
 *      n3 @---<---|-------------@ n5        e0 = (n0,n1)
 *         |\      *  |        ,/|           e1 = (n1,n2)
 *         | \        o      ,/  |           e2 = (n2,n0)
 *         |  \         e7 ,^`   |           e3 = (n0,n3)
 *         |   v e6      ,/`     |           e4 = (n1,n4)
 *      e3 ^    \      ,/`       ^ e5        e5 = (n2,n5)
 *         |     \   ,/`         |           e6 = (n3,n4)
 *         |      \ /`        *-------> f1   e7 = (n4,n5)
 *  f0 <-------*   @ n4          |           e8 = (n5,n3)
 *         |       |             |           f0 = (n0,n1,n4,n3)
 *      n0 @-------|---------<---@ n2        f1 = (n1,n2,n5,n4)
 *          \      |        e2 ,/            f2 = (n2,n0,n3,n5)
 *           \     |     o   ,/`             f3 = (n0,n2,n1)
 *            \    ^ e4  | ,/`               f4 = (n3,n4,n5)
 *          e0 v   |     ,^ e1
 *              \  |   ,/|
 *               \ | ,/` |
 *                \|/`   |
 *                 @ n1  v
 *                       f3
 * @endverbatim
 */
class cPentahedronShape final : public iSplittableShape {
public:
    uint_t num_shape_nodes() const final;
    tShapeParts get_shape_edges() const final;
    tShapeParts get_shape_faces() const final;

    tShapeParts simplex_split() const final;
};  // class cPyramidShape

/**
 * Hexahedron shape.
 * @verbatim
 *                      f1
 *                      ^   f2
 *                      |   ^
 *                   e9 |   |
 *            n6 @---<--|----------@ n5         e0 = (n0,n1)
 *              /|      |   |     /|            e1 = (n1,n2)
 *             / |      |   o    / |            e2 = (n2,n3)
 *        e10 v  |      *    e8 ^  |            e3 = (n3,n0)
 *           /   ^ e6          /   ^ e5         e4 = (n0,n4)
 *          /    |      e11   /  *-------> f1   e5 = (n1,n5)
 *      n7 @------------->---@ n4  |            e6 = (n2,n6)
 *  f3 <---|--o  |           |     |            e7 = (n3,n7)
 *         |  n2 @---<-------|-----@ n1         e8 = (n4,n5)
 *         |    /    e1      |    /             e9 = (n5,n6)
 *      e7 ^   /          e4 ^   /             e10 = (n6,n7)
 *         |  v e2  *        |  ^ e0           e11 = (n7,n4)
 *         | /      |   o    | /                f0 = (n0,n3,n2,n1)
 *         |/       |   |    |/                 f1 = (n0,n1,n5,n4)
 *      n3 @--->----|--------@ n0               f2 = (n1,n2,n6,n5)
 *             e3   |   |                       f3 = (n2,n3,n7,n6)
 *                  |   v                       f4 = (n0,n4,n7,n3)
 *                  v   f0                      f5 = (n4,n5,n6,n7)
 *                  f4
 * @endverbatim
 */
class cHexahedronShape final : public iSplittableShape {
public:
    uint_t num_shape_nodes() const final;
    tShapeParts get_shape_edges() const final;
    tShapeParts get_shape_faces() const final;

    tShapeParts simplex_split() const final;
};  // class cHexahedronShape

} // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif // SHAPES_HH_
