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

#include "Shape.hh"

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/**
 * Construct the shape.
 */
iShapePtr::iShapePtr(eShape shape_type) {
    switch (shape_type) {
    case eShape::segment_2:
        reset(new cSegmentShape());
        break;
    case eShape::triangle_3:
        reset(new cTriangleShape());
        break;
    case eShape::quadrangle_4:
        reset(new cQuadrangleShape());
        break;
    case eShape::tetrahedron_4:
        reset(new cTetrahedronShape());
        break;
    case eShape::pyramid_5:
        reset(new cPyramidShape());
        break;
    case eShape::pentahedron_6:
        reset(new cPentahedronShape());
        break;
    case eShape::hexahedron_8:
        reset(new cHexahedronShape());
        break;
    default:
        FEATHERS_NOT_REACHABLE();
    }
}   // iShape::construct

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

real_t cSegmentShape::get_diameter() const {
    const vec3_t delta =
        get_node_coords(1) - get_node_coords(0);
    return glm::length(delta);
}
real_t cSegmentShape::get_length_or_area_or_volume() const {
    return get_diameter();
}
vec3_t cSegmentShape::get_normal() const {
    const vec3_t delta =
        get_node_coords(1) - get_node_coords(0);
    static const vec3_t up(0.0, 0.0, 1.0);
    return glm::normalize(glm::cross(delta, up));
}
vec3_t cSegmentShape::get_direction() const {
    const vec3_t delta =
        get_node_coords(1) - get_node_coords(0);
    return glm::normalize(delta);
}
vec3_t cSegmentShape::get_center_coords() const {
    return 0.5*(get_node_coords(0) + get_node_coords(1));
}

real_t cTriangleShape::get_diameter() const {
    const vec3_t delta1 =
        get_node_coords(1) - get_node_coords(0);
    const vec3_t delta2 =
        get_node_coords(2) - get_node_coords(0);
    const vec3_t delta3 =
        get_node_coords(2) - get_node_coords(1);
    return std::max(glm::length(delta1),
                    std::max(glm::length(delta2), glm::length(delta3)));
}
real_t cTriangleShape::get_length_or_area_or_volume() const {
    const vec3_t delta1 =
        get_node_coords(1) - get_node_coords(0);
    const vec3_t delta2 =
        get_node_coords(2) - get_node_coords(0);
    return 0.5*glm::length(glm::cross(delta1, delta2));
}
vec3_t cTriangleShape::get_normal() const {
    const vec3_t delta1 =
        get_node_coords(1) - get_node_coords(0);
    const vec3_t delta2 =
        get_node_coords(2) - get_node_coords(0);
    return glm::normalize(glm::cross(delta1, delta2));
}
vec3_t cTriangleShape::get_center_coords() const {
    return (get_node_coords(0) +
            get_node_coords(1) + get_node_coords(2))/3.0;
}

real_t cTetrahedronShape::get_length_or_area_or_volume() const {
    const vec3_t delta1 =
        get_node_coords(1) - get_node_coords(0);
    const vec3_t delta2 =
        get_node_coords(2) - get_node_coords(0);
    const vec3_t delta3 =
        get_node_coords(3) - get_node_coords(0);
    return std::abs(glm::dot(delta1, glm::cross(delta2, delta3)))/6.0;
}
vec3_t cTetrahedronShape::get_center_coords() const {
    return 0.25*(get_node_coords(0) + get_node_coords(1) +
                 get_node_coords(2) + get_node_coords(3));
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

template<typename tFunc>
void iSplittableShape::for_each_part_(tFunc func) const {
    tSplittingScheme scheme(simplex_split());
    iShapePtr part_shape(scheme.part_shape);
    for (const auto& part : scheme.part_nodes) {
        part_shape->assign_node_coords(m_node_coords, part.begin(), part.end());
        func(*part_shape);
    }
}   // iSplittableShape::for_each_part_

real_t iSplittableShape::get_diameter() const {
    real_t diameter = 0.0;
    for_each_part_([&](iShape& shape) {
        diameter = std::max(diameter, shape.get_diameter());
    });
    return diameter;
}   // iSplittableShape::get_diameter

real_t iSplittableShape::get_length_or_area_or_volume() const {
    real_t length_or_area_or_volume = 0.0;
    for_each_part_([&](iShape& shape) {
        length_or_area_or_volume += shape.get_length_or_area_or_volume();
    });
    return length_or_area_or_volume;
}   // iSplittableShape::get_length_or_area_or_volume

vec3_t iSplittableShape::get_normal() const {
    vec3_t weighted_sum_of_normals(0.0);
    for_each_part_([&](iShape& shape) {
        real_t part_length_or_area_or_volume =
            shape.get_length_or_area_or_volume();
        weighted_sum_of_normals +=
            part_length_or_area_or_volume*shape.get_normal();
    });
    return glm::normalize(weighted_sum_of_normals);
}   // iSplittableShape::get_normal

vec3_t iSplittableShape::get_center_coords() const {
    vec3_t weighted_sum_of_center_coords(0.0);
    real_t length_or_area_or_volume = 0.0;
    for_each_part_([&](iShape& shape) {
        real_t part_length_or_area_or_volume =
            shape.get_length_or_area_or_volume();
        weighted_sum_of_center_coords +=
            part_length_or_area_or_volume*shape.get_center_coords();
        length_or_area_or_volume += part_length_or_area_or_volume;
    });
    return weighted_sum_of_center_coords/length_or_area_or_volume;
}   // iSplittableShape::get_center_coords

tSplittingScheme cQuadrangleShape::simplex_split() const {
    return tSplittingScheme{ eShape::triangle_3, { {0, 1, 2}, {2, 3, 1} } };
}
tSplittingScheme cPyramidShape::simplex_split() const {
    return tSplittingScheme{ eShape::tetrahedron_4, { {0, 1, 3, 4}, {1, 2, 3, 4} } };
}
tSplittingScheme cPentahedronShape::simplex_split() const {
    FEATHERS_NOT_IMPLEMENTED();
}
tSplittingScheme cHexahedronShape::simplex_split() const {
    return tSplittingScheme{ eShape::tetrahedron_4, { {0, 3, 1, 4}, {3, 2, 1, 6},
                                        {4, 5, 6, 1}, {4, 6, 7, 3}, {4, 3, 1, 6} } };
}

} // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
