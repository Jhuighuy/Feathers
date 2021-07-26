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

#include "Element.hh"

#include <functional>

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/**
 * Construct the element.
 */
static std::unique_ptr<iElement> construct_element_(eShape shape) {
    switch (shape) {
    case eShape::segment_2:
        return std::make_unique<cSegment>();
    case eShape::triangle_3:
        return std::make_unique<cTriangle>();
    case eShape::quadrangle_4:
        return std::make_unique<cQuadrangle>();
    case eShape::tetrahedron_4:
        return std::make_unique<cTetrahedron>();
    case eShape::pyramid_5:
        return std::make_unique<cPyramid>();
    case eShape::pentahedron_6:
        return std::make_unique<cPentahedron>();
    case eShape::hexahedron_8:
        return std::make_unique<cHexahedron>();
    default:
        FEATHERS_NOT_REACHABLE();
    }
}   // construct_element_

/**
 * Construct a new element object.
 */
std::unique_ptr<iElement> iElement::make(sElementDesc&& desc,
                                         uint_t num_global_nodes,
                                         const vec3_t* global_node_coords) {
    std::unique_ptr<iElement> element = construct_element_(desc.shape);
    /* Assign a global node array. */
    element->m_num_global_nodes = num_global_nodes;
    element->m_global_node_coords = global_node_coords;
    /* Verify and assign node indices. */
    FEATHERS_ASSERT(
        desc.node_indices.size() <= num_global_nodes &&
        desc.node_indices.size() == element->num_nodes() &&
        std::all_of(desc.node_indices.begin(), desc.node_indices.end(),
                    [&](uint_t node_index) { return node_index < num_global_nodes; }));
    element->m_node_indices = std::move(desc.node_indices);
    return element;
}   // iElement::make

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

template<typename tFunc>
void iComplexElement::for_each_simplex_(tFunc func) const {
    tElementDescList simplices(get_simplicial_parts(0));
    for (sElementDesc& part : simplices) {
        std::unique_ptr<iElement> simplex =
            iElement::make(std::move(part), m_num_global_nodes, m_global_node_coords);
        func(*simplex);
    }
}   // iComplexElement::for_each_simplex_

real_t iComplexElement::get_diameter() const {
    real_t diameter = 0.0;
    for_each_simplex_([&](iElement& shape) {
        diameter = std::max(diameter, shape.get_diameter());
    });
    return diameter;
}   // iComplexElement::get_diameter

real_t iComplexElement::get_length_or_area_or_volume() const {
    real_t length_or_area_or_volume = 0.0;
    for_each_simplex_([&](iElement& shape) {
        length_or_area_or_volume += shape.get_length_or_area_or_volume();
    });
    return length_or_area_or_volume;
}   // iComplexElement::get_length_or_area_or_volume

vec3_t iComplexElement::get_normal() const {
    vec3_t weighted_sum_of_normals(0.0);
    for_each_simplex_([&](iElement& shape) {
        real_t part_length_or_area_or_volume =
            shape.get_length_or_area_or_volume();
        weighted_sum_of_normals +=
            part_length_or_area_or_volume * shape.get_normal();
    });
    return glm::normalize(weighted_sum_of_normals);
}   // iComplexElement::get_normal

vec3_t iComplexElement::get_center_coords() const {
    vec3_t weighted_sum_of_center_coords(0.0);
    real_t length_or_area_or_volume = 0.0;
    for_each_simplex_([&](iElement& shape) {
        real_t part_length_or_area_or_volume =
            shape.get_length_or_area_or_volume();
        weighted_sum_of_center_coords +=
            part_length_or_area_or_volume * shape.get_center_coords();
        length_or_area_or_volume += part_length_or_area_or_volume;
    });
    return weighted_sum_of_center_coords/length_or_area_or_volume;
}   // iComplexElement::get_center_coords

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

real_t cNode::get_diameter() const {
    return 0.0;
}
real_t cNode::get_length_or_area_or_volume() const {
    return 1.0;
}
vec3_t cNode::get_normal() const {
    static const vec3_t right(1.0, 0.0, 0.0);
    return right;
}
vec3_t cNode::get_center_coords() const {
    return get_node_coords(0);
}
eShape cNode::get_shape() const {
    return eShape::node;
}
uint_t cNode::num_nodes() const {
    return 1;
}
tElementDescList cNode::get_edges_desc() const {
    return {};
}
tElementDescList cNode::get_faces_desc() const {
    return {};
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

real_t cSegment::get_diameter() const {
    const vec3_t delta =
        get_node_coords(1) - get_node_coords(0);
    return glm::length(delta);
}
real_t cSegment::get_length_or_area_or_volume() const {
    return get_diameter();
}
vec3_t cSegment::get_normal() const {
    const vec3_t delta =
        get_node_coords(1) - get_node_coords(0);
    static const vec3_t up(0.0, 0.0, 1.0);
    return glm::normalize(glm::cross(delta, up));
}
vec3_t cSegment::get_direction() const {
    const vec3_t delta =
        get_node_coords(1) - get_node_coords(0);
    return glm::normalize(delta);
}
vec3_t cSegment::get_center_coords() const {
    return 0.5*(get_node_coords(0) + get_node_coords(1));
}
eShape cSegment::get_shape() const {
    return eShape::segment_2;
}
uint_t cSegment::num_nodes() const {
    return 2;
}
tElementDescList cSegment::get_edges_desc() const {
    return { get_part(eShape::segment_2, 0, 1) };
}
tElementDescList cSegment::get_faces_desc() const {
    return { get_part(eShape::node, 0), get_part(eShape::node, 1) };
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

real_t cTriangle::get_diameter() const {
    const vec3_t delta1 =
        get_node_coords(1) - get_node_coords(0);
    const vec3_t delta2 =
        get_node_coords(2) - get_node_coords(1);
    const vec3_t delta3 =
        get_node_coords(0) - get_node_coords(2);
    return std::max(glm::length(delta1),
                    std::max(glm::length(delta2), glm::length(delta3)));
}
real_t cTriangle::get_length_or_area_or_volume() const {
    const vec3_t delta1 =
        get_node_coords(1) - get_node_coords(0);
    const vec3_t delta2 =
        get_node_coords(2) - get_node_coords(0);
    return 0.5*glm::length(glm::cross(delta1, delta2));
}
vec3_t cTriangle::get_normal() const {
    const vec3_t delta1 =
        get_node_coords(1) - get_node_coords(0);
    const vec3_t delta2 =
        get_node_coords(2) - get_node_coords(0);
    return glm::normalize(glm::cross(delta1, delta2));
}
vec3_t cTriangle::get_center_coords() const {
    return (get_node_coords(0) +
            get_node_coords(1) + get_node_coords(2))/3.0;
}
eShape cTriangle::get_shape() const {
    return eShape::triangle_3;
}
uint_t cTriangle::num_nodes() const {
    return 3;
}
tElementDescList cTriangle::get_edges_desc() const {
    return { get_part(eShape::segment_2, 0, 1),
             get_part(eShape::segment_2, 1, 2),
             get_part(eShape::segment_2, 2, 0) };
}
tElementDescList cTriangle::get_faces_desc() const {
    return get_edges_desc();
}

uint_t cQuadrangle::num_nodes() const {
    return 4;
}
eShape cQuadrangle::get_shape() const {
    return eShape::quadrangle_4;
}
tElementDescList cQuadrangle::get_edges_desc() const {
    return { get_part(eShape::segment_2, 0, 1),
             get_part(eShape::segment_2, 1, 2),
             get_part(eShape::segment_2, 2, 3),
             get_part(eShape::segment_2, 3, 0) };
}
tElementDescList cQuadrangle::get_faces_desc() const {
    return get_edges_desc();
}
tElementDescList cQuadrangle::get_simplicial_parts(uint_t partition_index) const {
    switch (partition_index) {
    case 0:
        return { get_part(eShape::triangle_3, 0, 1, 2),
                 get_part(eShape::triangle_3, 2, 3, 0) };
    case 1:
        return { get_part(eShape::triangle_3, 0, 1, 3),
                 get_part(eShape::triangle_3, 1, 2, 3) };
    default:
        return {};
    }
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

real_t cTetrahedron::get_diameter() const {
    const vec3_t delta1 =
        get_node_coords(1) - get_node_coords(0);
    const vec3_t delta2 =
        get_node_coords(2) - get_node_coords(1);
    const vec3_t delta3 =
        get_node_coords(0) - get_node_coords(2);
    const vec3_t delta4 =
        get_node_coords(3) - get_node_coords(0);
    const vec3_t delta5 =
        get_node_coords(3) - get_node_coords(1);
    const vec3_t delta6 =
        get_node_coords(3) - get_node_coords(2);
    return std::max({ glm::length(delta1), glm::length(delta2),
                      glm::length(delta3), glm::length(delta4),
                      glm::length(delta5), glm::length(delta6) });
}
real_t cTetrahedron::get_length_or_area_or_volume() const {
    const vec3_t delta1 =
        get_node_coords(1) - get_node_coords(0);
    const vec3_t delta2 =
        get_node_coords(2) - get_node_coords(0);
    const vec3_t delta3 =
        get_node_coords(3) - get_node_coords(0);
    return std::abs(glm::dot(delta1,
                             glm::cross(delta2, delta3)))/6.0;
}
vec3_t cTetrahedron::get_center_coords() const {
    return 0.25*(get_node_coords(0) + get_node_coords(1) +
                 get_node_coords(2) + get_node_coords(3));
}
eShape cTetrahedron::get_shape() const {
    return eShape::tetrahedron_4;
}
uint_t cTetrahedron::num_nodes() const {
    return 4;
}
tElementDescList cTetrahedron::get_edges_desc() const {
    return { get_part(eShape::segment_2, 0, 1),
             get_part(eShape::segment_2, 1, 2),
             get_part(eShape::segment_2, 2, 0),
             get_part(eShape::segment_2, 0, 3),
             get_part(eShape::segment_2, 1, 3),
             get_part(eShape::segment_2, 2, 3) };
}
tElementDescList cTetrahedron::get_faces_desc() const {
    return { get_part(eShape::triangle_3, 0, 2, 1),
             get_part(eShape::triangle_3, 0, 1, 3),
             get_part(eShape::triangle_3, 1, 2, 3),
             get_part(eShape::triangle_3, 2, 0, 3) };
}

uint_t cPyramid::num_nodes() const {
    return 5;
}
eShape cPyramid::get_shape() const {
    return eShape::pyramid_5;
}
tElementDescList cPyramid::get_edges_desc() const {
    return { get_part(eShape::segment_2, 0, 1),
             get_part(eShape::segment_2, 1, 2),
             get_part(eShape::segment_2, 2, 3),
             get_part(eShape::segment_2, 3, 0),
             get_part(eShape::segment_2, 0, 4),
             get_part(eShape::segment_2, 1, 4),
             get_part(eShape::segment_2, 2, 4),
             get_part(eShape::segment_2, 3, 4) };
}
tElementDescList cPyramid::get_faces_desc() const {
    return { get_part(eShape::quadrangle_4, 0, 3, 2, 1),
             get_part(eShape::triangle_3, 0, 1, 4),
             get_part(eShape::triangle_3, 1, 2, 4),
             get_part(eShape::triangle_3, 2, 3, 4),
             get_part(eShape::triangle_3, 3, 0, 4) };
}
tElementDescList cPyramid::get_simplicial_parts(uint_t partition_index) const {
    switch (partition_index) {
    case 0:
        return { get_part(eShape::tetrahedron_4, 0, 1, 2, 4),
                 get_part(eShape::tetrahedron_4, 2, 3, 0, 4) };
    case 1:
        return { get_part(eShape::tetrahedron_4, 0, 1, 3, 4),
                 get_part(eShape::tetrahedron_4, 1, 2, 3, 4) };
    default:
        return {};
    }
}

uint_t cPentahedron::num_nodes() const {
    return 6;
}
eShape cPentahedron::get_shape() const {
    return eShape::pentahedron_6;
}
tElementDescList cPentahedron::get_edges_desc() const {
    return { get_part(eShape::segment_2, 0, 1),
             get_part(eShape::segment_2, 1, 2),
             get_part(eShape::segment_2, 2, 0),
             get_part(eShape::segment_2, 0, 3),
             get_part(eShape::segment_2, 1, 4),
             get_part(eShape::segment_2, 2, 5),
             get_part(eShape::segment_2, 3, 4),
             get_part(eShape::segment_2, 4, 5),
             get_part(eShape::segment_2, 5, 3) };
}
tElementDescList cPentahedron::get_faces_desc() const {
    return { get_part(eShape::quadrangle_4, 0, 1, 4, 3),
             get_part(eShape::quadrangle_4, 1, 2, 5, 4),
             get_part(eShape::quadrangle_4, 2, 0, 3, 5),
             get_part(eShape::triangle_3, 0, 2, 1),
             get_part(eShape::triangle_3, 3, 4, 5) };
}
tElementDescList cPentahedron::get_simplicial_parts(uint_t partition_index) const {
    FEATHERS_NOT_IMPLEMENTED();
}

uint_t cHexahedron::num_nodes() const {
    return 8;
}
eShape cHexahedron::get_shape() const {
    return eShape::hexahedron_8;
}
tElementDescList cHexahedron::get_edges_desc() const {
    return { get_part(eShape::segment_2, 0, 1),
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
             get_part(eShape::segment_2, 7, 4) };
}
tElementDescList cHexahedron::get_faces_desc() const {
    return { get_part(eShape::quadrangle_4, 0, 3, 2, 1),
             get_part(eShape::quadrangle_4, 0, 1, 5, 4),
             get_part(eShape::quadrangle_4, 1, 2, 6, 5),
             get_part(eShape::quadrangle_4, 2, 3, 7, 6),
             get_part(eShape::quadrangle_4, 0, 4, 7, 3),
             get_part(eShape::quadrangle_4, 4, 5, 6, 7) };
}
tElementDescList cHexahedron::get_simplicial_parts(uint_t partition_index) const {
    // TODO: alternative partitions!
    return { get_part(eShape::tetrahedron_4, 0, 3, 1, 4),
             get_part(eShape::tetrahedron_4, 3, 2, 1, 6),
             get_part(eShape::tetrahedron_4, 4, 5, 6, 1),
             get_part(eShape::tetrahedron_4, 4, 6, 7, 3),
             get_part(eShape::tetrahedron_4, 4, 3, 1, 6) };
}

} // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#if !FEATHERS_DOXYGEN
/**
 * Alternative pyramid pictures.
 * @verbatim
 *                            n4
 *                            @    f1
 *                           /|\   ^
 *                          /`|`\./
 *                        ,// | `/
 *                       /,/  | /\.
 *                     ,/ |`  |/ `\
 *          f3     e7 ^` ,|   /   ^ e5
 *           ^      ,/   /`  /|   \.
 *            \    /`   /   * |   `\
 *             \ ,/   ,/      |    |
 *              /`    |`   e6 ^    \.
 *            ,/ \    ^ e4    |  o-`|---------> f2
 *           /`   o  ,|       |     |
 *   <----------*    /`       |     \.
 *       ,/`        /         |     `\
 *   n3 @----<----,/----------@ n2   |
 *       \.  e2   |`           \.    \.
 *        `\.    ,|          e1 `^.  `\
 *      e3 `v.   /`     o         `\. |
 *           `\ /       |            `\|
 *         n0 @------------------->---@ n1
 *                      |        e0
 *                      |
 *                      v
 *                      f0
 * @endverbatim
 *
 * @verbatim
 *                           f2
 *                     n4    ^
 *                     @.   /
 *                   ,/|`\./
 *              e7 ,^/ |  /\.
 *               ,/ /  | *  ^ e6
 *             ,/`,/   ^ e5 `\.
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
#endif
