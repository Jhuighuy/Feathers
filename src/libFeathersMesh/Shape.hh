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
 * Abstract shape class.
 */
class iShape {
protected:
    const vec3_t* m_node_coords = nullptr;
    std::vector<uint_t> m_local_indices;

public:
    virtual ~iShape() = default;

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get shape node position. */
    const vec3_t& get_node_coords(uint_t node_local) const {
        FEATHERS_ASSERT(node_local < m_local_indices.size());
        return m_node_coords[m_local_indices[node_local]];
    }
    /** Assign shape nodes. */
    template<class tIter>
    void assign_node_coords(const vec3_t* node_coords,
                            tIter begin_local_indices, tIter end_local_indices) {
        m_node_coords = node_coords;
        m_local_indices.assign(begin_local_indices, end_local_indices);
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
};  // class iShape

/**
 * Shape pointer convenience wrapper.
 */
class iShapePtr final : public std::unique_ptr<iShape> {
public:
    /** Construct the shape. */
    explicit iShapePtr(eShape shape_type);
};  // class iShapePtr

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Splitting scheme for simplex splittable shapes.
 */
struct tSplittingScheme {
    eShape part_shape;
    std::vector<std::vector<uint_t>> part_nodes;
};  // struct tSplittingScheme

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
    virtual tSplittingScheme simplex_split() const = 0;

private:
    template<typename tFunc>
    void for_each_part_(tFunc func) const;
};  // class iSplittableShape

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/**
 * Segment shape.
 */
class cSegmentShape final : public iShape {
public:
    real_t get_diameter() const final;
    real_t get_length_or_area_or_volume() const final;
    vec3_t get_normal() const final;
    vec3_t get_direction() const final;
    vec3_t get_center_coords() const final;
};  // class tSegmentShape

/**
 * Triangle shape.
 */
class cTriangleShape final : public iShape {
public:
    real_t get_diameter() const final;
    real_t get_length_or_area_or_volume() const final;
    vec3_t get_normal() const final;
    vec3_t get_center_coords() const final;
};  // class cTriangleShape

/**
 * Quadrangle shape.
 */
class cQuadrangleShape final : public iSplittableShape {
public:
    tSplittingScheme simplex_split() const final;
};  // class cQuadrangleShape

/**
 * Tetrahedron shape.
 */
class cTetrahedronShape final : public iShape {
public:
    real_t get_length_or_area_or_volume() const final;
    vec3_t get_center_coords() const final;
};  // class cTetrahedronShape

/**
 * Pyramid shape.
 */
class cPyramidShape final : public iSplittableShape {
public:
    tSplittingScheme simplex_split() const final;
};  // class cPyramidShape

/**
 * Pentahedron (triangular prism) shape.
 */
class cPentahedronShape final : public iSplittableShape {
public:
    tSplittingScheme simplex_split() const final;
};  // class cPyramidShape

/**
 * Hexahedron shape.
 */
class cHexahedronShape final : public iSplittableShape {
public:
    tSplittingScheme simplex_split() const final;
};  // class cHexahedronShape

} // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif // SHAPES_HH_
