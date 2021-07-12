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
    /**
     * Construct the shape.
     */
    static std::unique_ptr<iShape> construct(eShape shape_type);

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

/**
 * Edge (half-open interval) in 2D space.
 */
struct mhd_e2d final
{
    using mhd_vec2_t = feathers::vec2_t;
    mhd_vec2_t p1{};
    mhd_vec2_t p2{};

public:
    mhd_vec2_t vec() const
    {
        return p2 - p1;
    }

public:
    static real_t dot(const mhd_e2d& e1, const mhd_e2d& e2)
    {
        mhd_vec2_t p1{e1.vec()};
        mhd_vec2_t p2{e2.vec()};
        return glm::dot(p1, p2);
    }
    static real_t len(const mhd_e2d& e1)
    {
        mhd_vec2_t p1{e1.vec()};
        return glm::length(p1);
    }

    static real_t dist(const mhd_e2d& e, const mhd_vec2_t& p)
    {
        mhd_e2d e1{p, e.p1};
        mhd_e2d e2{p, e.p2};
        return det(e2, e1) / len(e);
    }

public:
    static real_t det(const mhd_e2d& e1, const mhd_e2d& e2)
    {
        mhd_vec2_t p1{e1.vec()};
        mhd_vec2_t p2{e2.vec()};
        return glm::determinant(feathers::mat2_t(p1, p2));
    }

public:
    static mhd_vec2_t normal(const mhd_e2d& e)
    {
        mhd_vec2_t p = e.vec();
        mhd_vec2_t n{p.y, -p.x};
        real_t l{glm::length(n)};
        if (l != 0.0) {
            return n / l;
        } else {
            std::cerr << "Warning: normal to singular edge." << std::endl;
            return {1.0, 0.0};
        }
    }
    static feathers::vec3_t normal3(const mhd_e2d& e) {
        const auto n{ normal(e) };
        return { n.x, n.y, 0.0 };
    }

};	// struct mhd_e2d

typedef mhd_e2d mhd_seg2_t;

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Edge (half-open interval) in 2D space.
 */
struct mhd_e3d final
{
    using mhd_vec3_t = feathers::vec3_t;
    mhd_vec3_t p1{};
    mhd_vec3_t p2{};

public:
    
    mhd_vec3_t vec() const
    {
        return p2 - p1;
    }

public:
    
    static real_t dot(const mhd_e3d& e1, const mhd_e3d& e2)
    {
        mhd_vec3_t p1{e1.vec()};
        mhd_vec3_t p2{e2.vec()};
        return glm::dot(p1, p2);
    }
    
    static real_t len(const mhd_e3d& e1)
    {
        mhd_vec3_t p1{e1.vec()};
        return glm::length(p1);
    }

public:
    
    static mhd_vec3_t cross(const mhd_e3d& e1, const mhd_e3d& e2)
    {
        mhd_vec3_t p1{e1.vec()};
        mhd_vec3_t p2{e2.vec()};
        return glm::cross(p1, p2);
    }

public:
    
    static mhd_vec3_t normal(const mhd_e3d& e1, const mhd_e3d& e2)
    {
        mhd_vec3_t n{cross(e1, e2)};
        real_t l{glm::length(n)};
        if (l != 0.0) {
            return n / l;
        } else {
            std::cerr << "Warning: normal to singular edges." << std::endl;
            return {1.0, 0.0, 0.0};
        }
    }
};  // struct mhd_e3d
