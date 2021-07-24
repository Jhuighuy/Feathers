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
 * Element part description.
 */
struct sElementPart {
    eShape shape;
    std::vector<uint_t> node_indices;

    /* Construct an element part. */
    template<typename... tIndex>
    explicit sElementPart(eShape shape, tIndex... node_indices):
        shape(shape), node_indices{node_indices...} {
    }
};  // sElementPart

/**
 * Array of shape parts.
 */
using tElementPartList = std::vector<sElementPart>;

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Element pointer convenience wrapper.
 */
class iElementPtr final : public std::unique_ptr<class iElement> {
public:
    /** Construct the element. */
    explicit iElementPtr(eShape shape);
};  // class iElementPtr

/**
 * Abstract element class.
 */
class iElement {
protected:
    uint_t m_num_nodes = 0;
    const vec3_t* m_node_coords = nullptr;
    std::vector<uint_t> m_node_indices;

public:
    virtual ~iElement() = default;

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Assign element nodes. */
    template<typename tIndexIter>
    void assign_nodes(uint_t num_nodes, const vec3_t* node_coords,
                      tIndexIter first_shape_node_index,
                      tIndexIter last_shape_node_index) {
        m_num_nodes = num_nodes, m_node_coords = node_coords;
        FEATHERS_ASSERT(this->num_nodes() ==
                        (last_shape_node_index - first_shape_node_index));
        m_node_indices.assign(
            first_shape_node_index, last_shape_node_index);
    }

    /** Get node position. */
    const vec3_t& get_node_coords(uint_t node_local) const {
        FEATHERS_ASSERT(node_local < m_node_indices.size());
        const uint_t node_index = m_node_indices[node_local];
        FEATHERS_ASSERT(node_index < m_num_nodes);
        return m_node_coords[node_index];
    }

    /** Get element part. */
    template<typename... tIndex>
    sElementPart get_part(eShape part_shape, tIndex... node_locals) const {
        FEATHERS_ASSERT(
            iElementPtr(part_shape)->num_nodes() == sizeof...(node_locals));
        return sElementPart(part_shape, m_node_indices[node_locals]...);
    }

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get element diameter. */
    virtual real_t get_diameter() const {
        return qnan;
    }
    /** Get element length/area/volume. */
    virtual real_t get_length_or_area_or_volume() const {
        return qnan;
    }
    /** Get normal to element. */
    virtual vec3_t get_normal() const {
        return vec3_t(qnan);
    }
    /** Get element direction. */
    virtual vec3_t get_direction() const {
        return vec3_t(qnan);
    }
    /** Get element barycenter. */
    virtual vec3_t get_center_coords() const {
        return vec3_t(qnan);
    }

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Number of nodes in the element. */
    virtual uint_t num_nodes() const = 0;
    /** Get element edges. */
    virtual tElementPartList get_edges() const = 0;
    /** Get element faces. */
    virtual tElementPartList get_faces() const = 0;
};  // class iElement

/**
 * Abstract complex (not simplex) element class.
 */
class iComplexElement : public iElement {
public:
    real_t get_diameter() const final;
    real_t get_length_or_area_or_volume() const final;
    vec3_t get_normal() const final;
    vec3_t get_center_coords() const final;

    /** Get splitting into the simplex parts. */
    virtual tElementPartList get_simplicial_parts() const = 0;

private:
    template<typename tFunc>
    void for_each_part_(tFunc func) const;
};  // class iComplexElement

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Dummy nodal element.
 */
class cNode final : public iElement {
public:
    real_t get_diameter() const final;
    real_t get_length_or_area_or_volume() const final;
    vec3_t get_normal() const final;
    vec3_t get_center_coords() const final;

    uint_t num_nodes() const final;
    tElementPartList get_edges() const final;
    tElementPartList get_faces() const final;
};  // class cNode

/**
 * Segmental element.
 * @verbatim
 *  n0 @ f0
 *      \         e0 = (n0,n1)
 *       v e0     f0 = (n0)
 *        \       f1 = (n1)
 *      n1 @ f1
 * @endverbatim
 */
class cSegment final : public iElement {
public:
    real_t get_diameter() const final;
    real_t get_length_or_area_or_volume() const final;
    vec3_t get_normal() const final;
    vec3_t get_direction() const final;
    vec3_t get_center_coords() const final;

    uint_t num_nodes() const final;
    tElementPartList get_edges() const final;
    tElementPartList get_faces() const final;
};  // class tSegmentShape

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Triangular element.
 * @verbatim
 *          @ n2        e0 = f0 = (n0,n1)
 *         / \          e1 = f1 = (n1,n2)
 *  e2/f2 v   ^ e1/f1   e2 = f2 = (n2,n0)
 *       /     \
 *   n0 @--->---@ n1
 *        e0/f0
 * @endverbatim
 */
class cTriangle final : public iElement {
public:
    real_t get_diameter() const final;
    real_t get_length_or_area_or_volume() const final;
    vec3_t get_normal() const final;
    vec3_t get_center_coords() const final;

    uint_t num_nodes() const final;
    tElementPartList get_edges() const final;
    tElementPartList get_faces() const final;
};  // class cTriangle

/**
 * Quadrangular element.
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
class cQuadrangle final : public iComplexElement {
public:
    uint_t num_nodes() const final;
    tElementPartList get_edges() const final;
    tElementPartList get_faces() const final;

    tElementPartList get_simplicial_parts() const final;
};  // class cQuadrangle

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Tetrahedral element.
 * @verbatim
 *                    f3
 *               n3   ^
 *                @   |
 *         f1    /|\. |     f2         e0 = (n0,n1)
 *         ^    / | `\.     ^          e1 = (n1,n2)
 *          \  /  |f1 `\.  /           e2 = (n2,n0)
 *           \`   |   | `\/            e3 = (n0,n3)
 *           ,\   |   o  /`\           e4 = (n1,n3)
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
class cTetrahedron final : public iElement {
public:
    real_t get_length_or_area_or_volume() const final;
    vec3_t get_center_coords() const final;

    uint_t num_nodes() const final;
    tElementPartList get_edges() const final;
    tElementPartList get_faces() const final;
};  // class cTetrahedron

/**
 * Pyramidal element.
 * @verbatim
 *                                n4                      e0 = (n0,n1)
 *                  f3           ,@                       e1 = (n1,n2)
 *                   ^        ,/`/|\     f1               e2 = (n2,n3)
 *                    \    ,/`  / | \    ^                e3 = (n3,n0)
 *                     \,/`    /  |  \  /                 e4 = (n0,n4)
 *                e7 ,/`\     /   |   \/                  e5 = (n1,n4)
 *                ,^`    o   /    |   /\                  e6 = (n2,n4)
 *             ,/`          /     |  *  \                 e7 = (n3,n4)
 *  f4 <~~~~~~~~~~~~*      /   e6 ^   o~~\~~~~~~~~~> f2   f0 = (n0,n3,n2,n1)
 *       ,/`              /       |       \               f1 = (n0,n1,n4)
 *   n3 @-----<----------/--------@  n2    ^ e5           f2 = (n1,n2,n4)
 *       `\.  e2        /          `\.      \             f3 = (n2,n3,n4)
 *          `>.        ^ e4           `\. e1 \            f4 = (n3,n0,n4)
 *          e3 `\.    /       o          `<.  \        split = ((n0,n1,n2,n4),
 *                `\./        |             `\.\                (n2,n3,n0,n4))
 *               n0 @-------------------->-----@ n1
 *                            |          e0
 *                            |
 *                            v
 *                            fo
 * @endverbatim
 */
class cPyramid final : public iComplexElement {
public:
    uint_t num_nodes() const final;
    tElementPartList get_edges() const final;
    tElementPartList get_faces() const final;

    tElementPartList get_simplicial_parts() const final;
};  // class cPyramid

/**
 * Pentahedral element (triangular prism).
 * @verbatim
 *                 f4
 *                 ^  f2
 *                 |  ^
 *             e8  |  |
 *      n3 @---<---|-------------@ n5        e0 = (n0,n1)
 *         |\      *  |        ,/|           e1 = (n1,n2)
 *         | \        o      ,/` |           e2 = (n2,n0)
 *         |  \         e7 ,^`   |           e3 = (n0,n3)
 *         |   v e6      ,/`     |           e4 = (n1,n4)
 *      e3 ^    \      ,/`       ^ e5        e5 = (n2,n5)
 *         |     \   ,/`         |           e6 = (n3,n4)
 *         |      \ /`        *~~~~~~~> f1   e7 = (n4,n5)
 *  f0 <~~~~~~~*   @ n4          |           e8 = (n5,n3)
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
class cPentahedron final : public iComplexElement {
public:
    uint_t num_nodes() const final;
    tElementPartList get_edges() const final;
    tElementPartList get_faces() const final;

    tElementPartList get_simplicial_parts() const final;
};  // class cPyramid

/**
 * Hexahedral element.
 * @verbatim
 *                      f1
 *                      ^   f2
 *                      |   ^
 *                   e9 |   |
 *            n6 @---<--|----------@ n5         e0 = (n0,n1)
 *              /|      |   |     /|            e1 = (n1,n2)
 *             / |      |   o    / |            e2 = (n2,n3)
 *        e10 v  |      *    e8 ^  ^ e5         e3 = (n3,n0)
 *           /   ^ e6          /   |            e4 = (n0,n4)
 *          /    |      e11   /  *~~~~~~~> f1   e5 = (n1,n5)
 *      n7 @------------->---@ n4  |            e6 = (n2,n6)
 *  f3 <~~~|~~o  |           |     |            e7 = (n3,n7)
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
class cHexahedron final : public iComplexElement {
public:
    uint_t num_nodes() const final;
    tElementPartList get_edges() const final;
    tElementPartList get_faces() const final;

    tElementPartList get_simplicial_parts() const final;
};  // class cHexahedron

} // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif // SHAPES_HH_
