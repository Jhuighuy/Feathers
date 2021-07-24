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

template<typename tFunc>
void iSplittableShape::for_each_part_(tFunc func) const {
    tShapeParts parts(simplex_split());
    for (const auto& part : parts) {
        iShapePtr part_shape(part.part_shape);
        part_shape->assign_nodes(
            m_num_nodes, m_node_coords,
            part.part_local_node_indices.begin(), part.part_local_node_indices.end());
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
uint_t cSegmentShape::num_shape_nodes() const {
    return 2;
}
tShapeParts cSegmentShape::get_shape_edges() const {
    return {
        get_part(eShape::segment_2, 0, 1),
    };
}
tShapeParts cSegmentShape::get_shape_faces() const {
    return {
        get_part(eShape::node, 0),
        get_part(eShape::node, 1),
    };
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

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
uint_t cTriangleShape::num_shape_nodes() const {
    return 3;
}
tShapeParts cTriangleShape::get_shape_edges() const {
    return {
        get_part(eShape::segment_2, 0, 1),
        get_part(eShape::segment_2, 1, 2),
        get_part(eShape::segment_2, 2, 0),
    };
}
tShapeParts cTriangleShape::get_shape_faces() const {
    return get_shape_edges();
}

uint_t cQuadrangleShape::num_shape_nodes() const {
    return 4;
}
tShapeParts cQuadrangleShape::get_shape_edges() const {
    return {
        get_part(eShape::segment_2, 0, 1),
        get_part(eShape::segment_2, 1, 2),
        get_part(eShape::segment_2, 2, 3),
        get_part(eShape::segment_2, 3, 0),
    };
}
tShapeParts cQuadrangleShape::get_shape_faces() const {
    return get_shape_edges();
}
tShapeParts cQuadrangleShape::simplex_split() const {
    return {
        get_part(eShape::triangle_3, 0, 1, 2),
        get_part(eShape::triangle_3, 2, 3, 0),
    };
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

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
uint_t cTetrahedronShape::num_shape_nodes() const {
    return 4;
}
tShapeParts cTetrahedronShape::get_shape_edges() const {
    return {
        get_part(eShape::segment_2, 0, 1),
        get_part(eShape::segment_2, 1, 2),
        get_part(eShape::segment_2, 2, 0),
        get_part(eShape::segment_2, 0, 3),
        get_part(eShape::segment_2, 1, 3),
        get_part(eShape::segment_2, 2, 3),
    };
}
tShapeParts cTetrahedronShape::get_shape_faces() const {
    return {
        get_part(eShape::triangle_3, 0, 2, 1),
        get_part(eShape::triangle_3, 0, 1, 3),
        get_part(eShape::triangle_3, 1, 2, 3),
        get_part(eShape::triangle_3, 2, 0, 3),
    };
}

uint_t cPyramidShape::num_shape_nodes() const {
    return 5;
}
tShapeParts cPyramidShape::get_shape_edges() const {
    FEATHERS_NOT_IMPLEMENTED();
}
tShapeParts cPyramidShape::get_shape_faces() const {
    return {
        get_part(eShape::quadrangle_4, 0, 3, 2, 1),
        get_part(eShape::triangle_3, 0, 1, 4),
        get_part(eShape::triangle_3, 1, 2, 4),
        get_part(eShape::triangle_3, 2, 3, 4),
        get_part(eShape::triangle_3, 3, 0, 4),
    };
}
tShapeParts cPyramidShape::simplex_split() const {
    return {
        get_part(eShape::tetrahedron_4, 0, 1, 2, 4),
        get_part(eShape::tetrahedron_4, 1, 2, 3, 4),
    };
}

uint_t cPentahedronShape::num_shape_nodes() const {
    return 6;
}
tShapeParts cPentahedronShape::get_shape_edges() const {
    return {
        get_part(eShape::segment_2, 0, 1),
        get_part(eShape::segment_2, 1, 2),
        get_part(eShape::segment_2, 2, 0),
        get_part(eShape::segment_2, 0, 3),
        get_part(eShape::segment_2, 1, 4),
        get_part(eShape::segment_2, 2, 5),
        get_part(eShape::segment_2, 3, 4),
        get_part(eShape::segment_2, 4, 5),
        get_part(eShape::segment_2, 5, 3),
    };
}
tShapeParts cPentahedronShape::get_shape_faces() const {
    return {
        get_part(eShape::quadrangle_4, 0, 1, 4, 3),
        get_part(eShape::quadrangle_4, 1, 2, 5, 4),
        get_part(eShape::quadrangle_4, 2, 0, 3, 5),
        get_part(eShape::triangle_3, 0, 2, 1),
        get_part(eShape::triangle_3, 3, 4, 5),
    };
}
tShapeParts cPentahedronShape::simplex_split() const {
    FEATHERS_NOT_IMPLEMENTED();
}

uint_t cHexahedronShape::num_shape_nodes() const {
    return 8;
}
tShapeParts cHexahedronShape::get_shape_edges() const {
    return {
        get_part(eShape::segment_2, 0, 1),
        get_part(eShape::segment_2, 1, 2),
        get_part(eShape::segment_2, 2, 3),
        get_part(eShape::segment_2, 3, 0),
        get_part(eShape::segment_2, 0, 4),
        get_part(eShape::segment_2, 1, 5),
        get_part(eShape::segment_2, 2, 6),
        get_part(eShape::segment_2, 3, 7),
        get_part(eShape::segment_2, 4, 5),
        get_part(eShape::segment_2, 5, 6),
        get_part(eShape::segment_2, 6, 7),
        get_part(eShape::segment_2, 7, 4),
    };
}
tShapeParts cHexahedronShape::get_shape_faces() const {
    return {
        get_part(eShape::quadrangle_4, 0, 3, 2, 1),
        get_part(eShape::quadrangle_4, 0, 1, 5, 4),
        get_part(eShape::quadrangle_4, 1, 2, 6, 5),
        get_part(eShape::quadrangle_4, 2, 3, 7, 6),
        get_part(eShape::quadrangle_4, 0, 4, 7, 3),
        get_part(eShape::quadrangle_4, 4, 5, 6, 7),
    };
}
tShapeParts cHexahedronShape::simplex_split() const {
    return {
        get_part(eShape::tetrahedron_4, 0, 3, 1, 4),
        get_part(eShape::tetrahedron_4, 3, 2, 1, 6),
        get_part(eShape::tetrahedron_4, 4, 5, 6, 1),
        get_part(eShape::tetrahedron_4, 4, 6, 7, 3),
        get_part(eShape::tetrahedron_4, 4, 3, 1, 6),
    };
}

} // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

/*
 * Alternative pyramid pictures.
 * @verbatim
 *                        @
 *                       /|\
 *                      /`|`\.
 *                    ,// | `|
 *                   / /  |  \.
 *                 ,/ /   |  `\
 *                /`  |   |   |
 *              ,/    |   |   \.
 *             /`    /    |   `\
 *           ,/     /     |    |
 *          /`     /      |    \.
 *        ,/      |       |    `|
 *       /`       |       |     |
 *     ,/        /        |     \.
 *   ,/`        /         |     `\
 *  @----------/----------@      |
 *   \.        |           \.    \.
 *    `\.      |            `\.  `\
 *      `\.   /               `\. |
 *        `\ /                  `\|
 *          @---------------------@
 *
 * @endverbatim
 *
 * @verbatim
 *                           f2
 *                     n4    ^
 *                     @.   /
 *                   ,/|`\./
 *              e7 ,^/ |  /\.
 *               ,/ /  | *  ^ e6
 *             ,/`./   ^ e5 `\.
 *           ,/` /`    |      `\
 *       n3 @---/------|---<---@ n2
 *         /   ^ e4    |  e2 ,/`
 *        /  ,/ *      |   ,/`
 *    e3 v /   / o     | ,^` e1
 *      /,/   /  |     |/`
 *  n0 @-----/----->---@ n1
 *          /    | e0
 *         /     |
 *        v      v
 *       f1      f0
 * @endverbatim
 */
