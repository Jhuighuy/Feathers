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
#ifndef SHAPES_HH_
#define SHAPES_HH_

#include "SkunkBase.hh"

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
 * Element description.
 */
struct sElementDesc {
    eShape shape;
    std::vector<size_t> node_indices;
};  // sElementDesc

/**
 * Array of shape descriptions.
 */
using tElementDescList = std::vector<sElementDesc>;

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Abstract element class.
 */
class iElement {
protected:
    size_t m_num_global_nodes = 0;
    const vec3_t* m_global_node_coords = nullptr;
    std::vector<size_t> m_node_indices;

public:

    /**
     * Construct a new element object.
     */
    static std::unique_ptr<iElement> make(sElementDesc&& desc,
                                          size_t num_global_nodes,
                                          const vec3_t* global_node_coords);

    virtual ~iElement() = default;

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get element node indices. */
    const std::vector<size_t>& get_nodes() const {
        return m_node_indices;
    }

    /** Get node position. */
    const vec3_t& get_node_coords(size_t node_local) const {
        FEATHERS_ASSERT(node_local < m_node_indices.size());
        return m_global_node_coords[m_node_indices[node_local]];
    }

    template<typename... tIndex>
    sElementDesc get_part(eShape part_shape, tIndex... node_locals) const {
        return { part_shape, std::vector<size_t>{ m_node_indices[node_locals]... } };
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

    /** Get element shape. */
    virtual eShape get_shape() const = 0;

    /** Number of nodes in the element. */
    virtual size_t num_nodes() const = 0;

    /** Number of edges in the element. */
    size_t num_edges() const {
        return get_edges_desc().size();
    }
    /** Get element edges description. */
    virtual tElementDescList get_edges_desc() const = 0;

    /** Number of faces in the element. */
    size_t num_faces() const {
        return get_faces_desc().size();
    }
    /** Get element faces description. */
    virtual tElementDescList get_faces_desc() const = 0;
};  // class iElement

/**
 * Abstract simplex element class.
 */
class iSimplexElement : public iElement {
};  // class iSimplexElement

/**
 * Abstract complex (not simplex) element class.
 */
class iComplexElement : public iElement {
public:
    real_t get_diameter() const final;
    real_t get_length_or_area_or_volume() const final;
    vec3_t get_normal() const final;
    vec3_t get_center_coords() const final;

    /**
     * Get splitting into the simplex parts.
     */
    virtual tElementDescList get_simplicial_parts(size_t partition_index) const = 0;

private:
    template<typename tFunc>
    void for_each_simplex_(tFunc func) const;
};  // class iComplexElement

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Dummy nodal element.
 */
class cNode final : public iSimplexElement {
public:
    real_t get_diameter() const final;
    real_t get_length_or_area_or_volume() const final;
    vec3_t get_normal() const final;
    vec3_t get_center_coords() const final;

    eShape get_shape() const final;
    size_t num_nodes() const final;
    tElementDescList get_edges_desc() const final;
    tElementDescList get_faces_desc() const final;
};  // class cNode

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Segmental element.
 * @verbatim
 *
 *  n0 O f0
 *      \
 *       \         e0 = (n0,n1)
 *        v e0     f0 = (n0)
 *         \       f1 = (n1)
 *          \
 *        n1 O f1
 *
 * @endverbatim
 */
class cSegment final : public iSimplexElement {
public:
    real_t get_diameter() const final;
    real_t get_length_or_area_or_volume() const final;
    vec3_t get_normal() const final;
    vec3_t get_direction() const final;
    vec3_t get_center_coords() const final;

    eShape get_shape() const final;
    size_t num_nodes() const final;
    tElementDescList get_edges_desc() const final;
    tElementDescList get_faces_desc() const final;
};  // class tSegmentShape

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Triangular element.
 * @verbatim
 *           n2
 *           O           e0 = f0 = (n0,n1)
 *          / \          e1 = f1 = (n1,n2)
 *         /   \         e2 = f2 = (n2,n0)
 *  e2/f2 v     ^ e1/f1
 *       /       \
 *      /         \
 *  n0 O----->-----O n1
 *        e0/f0
 * @endverbatim
 */
class cTriangle final : public iSimplexElement {
public:
    real_t get_diameter() const final;
    real_t get_length_or_area_or_volume() const final;
    vec3_t get_normal() const final;
    vec3_t get_center_coords() const final;

    eShape get_shape() const final;
    size_t num_nodes() const final;
    tElementDescList get_edges_desc() const final;
    tElementDescList get_faces_desc() const final;
};  // class cTriangle

/**
 * Quadrangular element.
 * @verbatim
 *               e2/f2
 *       n3 O-----<-----O n2    e0 = f0 = (n0,n1)
 *         /           /        e1 = f2 = (n1,n2)
 *  e3/f3 v           ^ e1/f1   e2 = f2 = (n2,n3)
 *       /           /          e3 = f3 = (n3,n0)
 *   n0 O----->-----O n1     split = ((n0,n1,n2),(n2,n3,n0))
 *          e0/f0
 * @endverbatim
 */
class cQuadrangle final : public iComplexElement {
public:
    eShape get_shape() const final;
    size_t num_nodes() const final;
    tElementDescList get_edges_desc() const final;
    tElementDescList get_faces_desc() const final;
    tElementDescList get_simplicial_parts(size_t partition_index) const final;
};  // class cQuadrangle

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Tetrahedral element.
 * @verbatim
 *                    f3
 *               n3   ^
 *                O   |
 *         f1    /|\. |     f2         e0 = (n0,n1)
 *         ^    / | `\.     ^          e1 = (n1,n2)
 *          \  /  |   `\.  /           e2 = (n2,n0)
 *           \`   |   | `\/            e3 = (n0,n3)
 *           ,\   |   o  /`\           e4 = (n1,n3)
 *       e3 ^  *  |     *   `^.e5      e5 = (n2,n3)
 *         /   e4 ^           `\       f0 = (n0,n2,n1)
 *     n0 O-------|--<----------O n2   f1 = (n0,n1,n3)
 *         \      |  e2       ,/       f2 = (n1,n2,n3)
 *          \     |     o   ,/`        f3 = (n2,n0,n3)
 *           \    ^ e4  | ,/`
 *         e0 v   |     ,^ e1
 *             \  |   ,/`
 *              \ | ,/` |
 *               \|/`   |
 *                O n1  v
 *                      f0
 * @endverbatim
 */
class cTetrahedron final : public iSimplexElement {
public:
    real_t get_diameter() const final;
    real_t get_length_or_area_or_volume() const final;
    vec3_t get_center_coords() const final;

    eShape get_shape() const final;
    size_t num_nodes() const final;
    tElementDescList get_edges_desc() const final;
    tElementDescList get_faces_desc() const final;
};  // class cTetrahedron

/**
 * Pyramidal element.
 * @verbatim
 *                                n4                      e0 = (n0,n1)
 *                  f3           ,O                       e1 = (n1,n2)
 *                   ^        ,/`/|\     f1               e2 = (n2,n3)
 *                    \    ,/`  / | \    ^                e3 = (n3,n0)
 *                     \,/`    /  |  \  /                 e4 = (n0,n4)
 *                e7 ,/`\     /   |   \/                  e5 = (n1,n4)
 *                ,^`    o   /    |   /\                  e6 = (n2,n4)
 *             ,/`          /     |  *  \                 e7 = (n3,n4)
 *  f4 <------------*      /   e6 ^   o--\---------> f2   f0 = (n0,n3,n2,n1)
 *       ,/`              /       |       \               f1 = (n0,n1,n4)
 *   n3 O-----<----------/--------O  n2    ^ e5           f2 = (n1,n2,n4)
 *       `\.  e2        /          `\.      \             f3 = (n2,n3,n4)
 *          `>.        ^ e4           `\. e1 \            f4 = (n3,n0,n4)
 *          e3 `\.    /       o          `<.  \        split = ((n0,n1,n2,n4),
 *                `\./        |             `\.\                (n2,n3,n0,n4))
 *               n0 O-------------------->-----O n1
 *                            |          e0
 *                            |
 *                            v
 *                            f0
 * @endverbatim
 */
class cPyramid final : public iComplexElement {
public:
    eShape get_shape() const final;
    size_t num_nodes() const final;
    tElementDescList get_edges_desc() const final;
    tElementDescList get_faces_desc() const final;
    tElementDescList get_simplicial_parts(size_t partition_index) const final;
};  // class cPyramid

/**
 * Pentahedral element (triangular prism).
 * @verbatim
 *                 f4
 *                 ^  f2
 *                 |  ^
 *             e8  |  |
 *      n3 O---<---|-------------O n5        e0 = (n0,n1)
 *         |\      *  |        ,/|           e1 = (n1,n2)
 *         | \        o      ,/` |           e2 = (n2,n0)
 *         |  \         e7 ,^`   |           e3 = (n0,n3)
 *         |   v e6      ,/`     |           e4 = (n1,n4)
 *      e3 ^    \      ,/`       ^ e5        e5 = (n2,n5)
 *         |     \   ,/`         |           e6 = (n3,n4)
 *         |      \ /`        *-------> f1   e7 = (n4,n5)
 *  f0 <-------*   @ n4          |           e8 = (n5,n3)
 *         |       |             |           f0 = (n0,n1,n4,n3)
 *      n0 O-------|---------<---O n2        f1 = (n1,n2,n5,n4)
 *          \      |        e2 ,/            f2 = (n2,n0,n3,n5)
 *           \     |     o   ,/`             f3 = (n0,n2,n1)
 *            \    ^ e4  | ,/`               f4 = (n3,n4,n5)
 *          e0 v   |     ,^ e1
 *              \  |   ,/|
 *               \ | ,/` |
 *                \|/`   |
 *                 O n1  v
 *                       f3
 * @endverbatim
 */
class cPentahedron final : public iComplexElement {
public:
    eShape get_shape() const final;
    size_t num_nodes() const final;
    tElementDescList get_edges_desc() const final;
    tElementDescList get_faces_desc() const final;
    tElementDescList get_simplicial_parts(size_t partition_index) const final;
};  // class cPyramid

/**
 * Hexahedral element.
 * @verbatim
 *                      f5
 *                      ^   f2
 *                      |   ^
 *                   e9 |   |
 *            n6 O---<--|----------O n5         e0 = (n0,n1)
 *              /|      |   |     /|            e1 = (n1,n2)
 *             / |      |   o    / |            e2 = (n2,n3)
 *        e10 v  |      *    e8 ^  ^ e5         e3 = (n3,n0)
 *           /   ^ e6          /   |            e4 = (n0,n4)
 *          /    |      e11   /  *-------> f1   e5 = (n1,n5)
 *      n7 O------------->---O n4  |            e6 = (n2,n6)
 *  f3 <---|--o  |           |     |            e7 = (n3,n7)
 *         |  n2 O---<-------|-----O n1         e8 = (n4,n5)
 *         |    /    e1      |    /             e9 = (n5,n6)
 *      e7 ^   /          e4 ^   /             e10 = (n6,n7)
 *         |  v e2  *        |  ^ e0           e11 = (n7,n4)
 *         | /      |   o    | /                f0 = (n0,n3,n2,n1)
 *         |/       |   |    |/                 f1 = (n0,n1,n5,n4)
 *      n3 O--->----|--------O n0               f2 = (n1,n2,n6,n5)
 *             e3   |   |                       f3 = (n2,n3,n7,n6)
 *                  |   v                       f4 = (n0,n4,n7,n3)
 *                  v   f0                      f5 = (n4,n5,n6,n7)
 *                  f4
 * @endverbatim
 */
class cHexahedron final : public iComplexElement {
public:
    eShape get_shape() const final;
    size_t num_nodes() const final;
    tElementDescList get_edges_desc() const final;
    tElementDescList get_faces_desc() const final;
    tElementDescList get_simplicial_parts(size_t partition_index) const final;
};  // class cHexahedron

} // namespace feathers

#endif // SHAPES_HH_
