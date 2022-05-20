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

#include "Element.hh"

#include <functional>

namespace feathers {

/**
 * Construct the element.
 */
static std::unique_ptr<iElement> construct_element_(ShapeType shape) {
    switch (shape) {
        case ShapeType::Segment2:
            return std::make_unique<Segment>();
        case ShapeType::Triangle3:
            return std::make_unique<Triangle>();
        case ShapeType::Quadrangle4:
            return std::make_unique<Quadrangle>();
        case ShapeType::Tetrahedron4:
            return std::make_unique<Tetrahedron>();
        case ShapeType::Pyramid5:
            return std::make_unique<Pyramid>();
        case ShapeType::Pentahedron6:
            return std::make_unique<Pentahedron>();
        case ShapeType::Hexahedron8:
            return std::make_unique<Hexahedron>();
        default:
            FEATHERS_NOT_REACHABLE();
    }
}   // construct_element_

/**
 * Construct a new element object.
 */
std::unique_ptr<iElement> iElement::make(ElementDesc&& desc,
                                         size_t num_global_nodes,
                                         const vec3_t* global_node_coords) {
    std::unique_ptr<iElement> element = construct_element_(desc.Shape);
    /* Assign a global Node array. */
    element->NumGlobalNodes_ = num_global_nodes;
    element->GlobalNodeCoords_ = global_node_coords;
    /* Verify and assign Node indices. */
    FEATHERS_ASSERT(
      desc.NodeIndices.size() <= num_global_nodes &&
      desc.NodeIndices.size() == element->num_nodes() &&
      std::all_of(desc.NodeIndices.begin(), desc.NodeIndices.end(),
                  [&](size_t node_index) { return node_index < num_global_nodes; }));
    element->NodeIndices_ = std::move(desc.NodeIndices);
    return element;
}   // iElement::make

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

template<typename tFunc>
void iComplexElement::for_each_simplex_(tFunc func) const {
    tElementDescList simplices(get_simplicial_parts(0));
    for (ElementDesc& part : simplices) {
        std::unique_ptr<iElement> simplex =
            iElement::make(std::move(part), NumGlobalNodes_, GlobalNodeCoords_);
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

real_t iComplexElement::LenAreaOrVolume() const {
    real_t length_or_area_or_volume = 0.0;
    for_each_simplex_([&](iElement& shape) {
        length_or_area_or_volume += shape.LenAreaOrVolume();
    });
    return length_or_area_or_volume;
}   // iComplexElement::LenAreaOrVolume

vec3_t iComplexElement::Normal() const {
    vec3_t weighted_sum_of_normals(0.0);
    for_each_simplex_([&](iElement& shape) {
        real_t part_length_or_area_or_volume =
          shape.LenAreaOrVolume();
        weighted_sum_of_normals +=
            part_length_or_area_or_volume * shape.Normal();
    });
    return glm::normalize(weighted_sum_of_normals);
}   // iComplexElement::Normal

vec3_t iComplexElement::CenterPos() const {
    vec3_t weighted_sum_of_center_coords(0.0);
    real_t length_or_area_or_volume = 0.0;
    for_each_simplex_([&](iElement& shape) {
        real_t part_length_or_area_or_volume =
          shape.LenAreaOrVolume();
        weighted_sum_of_center_coords +=
            part_length_or_area_or_volume * shape.CenterPos();
        length_or_area_or_volume += part_length_or_area_or_volume;
    });
    return weighted_sum_of_center_coords/length_or_area_or_volume;
}   // iComplexElement::CenterPos

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

real_t cNode::get_diameter() const {
    return 0.0;
}
real_t cNode::LenAreaOrVolume() const {
    return 1.0;
}
vec3_t cNode::Normal() const {
    static const vec3_t right(1.0, 0.0, 0.0);
    return right;
}
vec3_t cNode::CenterPos() const {
    return get_node_coords(0);
}
ShapeType cNode::Shape() const {
    return ShapeType::Node;
}
size_t cNode::num_nodes() const {
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

real_t Segment::get_diameter() const {
    const vec3_t delta =
        get_node_coords(1) - get_node_coords(0);
    return glm::length(delta);
}
real_t Segment::LenAreaOrVolume() const {
    return get_diameter();
}
vec3_t Segment::Normal() const {
    const vec3_t delta =
        get_node_coords(1) - get_node_coords(0);
    static const vec3_t up(0.0, 0.0, 1.0);
    return glm::normalize(glm::cross(delta, up));
}
vec3_t Segment::Dir() const {
    const vec3_t delta =
        get_node_coords(1) - get_node_coords(0);
    return glm::normalize(delta);
}
vec3_t Segment::CenterPos() const {
    return 0.5*(get_node_coords(0) + get_node_coords(1));
}
ShapeType Segment::Shape() const {
    return ShapeType::Segment2;
}
size_t Segment::num_nodes() const {
    return 2;
}
tElementDescList Segment::get_edges_desc() const {
    return { get_part(ShapeType::Segment2, 0, 1) };
}
tElementDescList Segment::get_faces_desc() const {
    return { get_part(ShapeType::Node, 0), get_part(ShapeType::Node, 1) };
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

real_t Triangle::get_diameter() const {
    const vec3_t delta1 =
        get_node_coords(1) - get_node_coords(0);
    const vec3_t delta2 =
        get_node_coords(2) - get_node_coords(1);
    const vec3_t delta3 =
        get_node_coords(0) - get_node_coords(2);
    return std::max(glm::length(delta1),
                    std::max(glm::length(delta2), glm::length(delta3)));
}
real_t Triangle::LenAreaOrVolume() const {
    const vec3_t delta1 =
        get_node_coords(1) - get_node_coords(0);
    const vec3_t delta2 =
        get_node_coords(2) - get_node_coords(0);
    return 0.5*glm::length(glm::cross(delta1, delta2));
}
vec3_t Triangle::Normal() const {
    const vec3_t delta1 =
        get_node_coords(1) - get_node_coords(0);
    const vec3_t delta2 =
        get_node_coords(2) - get_node_coords(0);
    return glm::normalize(glm::cross(delta1, delta2));
}
vec3_t Triangle::CenterPos() const {
    return (get_node_coords(0) +
            get_node_coords(1) + get_node_coords(2))/3.0;
}
ShapeType Triangle::Shape() const {
    return ShapeType::Triangle3;
}
size_t Triangle::num_nodes() const {
    return 3;
}
tElementDescList Triangle::get_edges_desc() const {
    return { get_part(ShapeType::Segment2, 0, 1),
             get_part(ShapeType::Segment2, 1, 2),
             get_part(ShapeType::Segment2, 2, 0) };
}
tElementDescList Triangle::get_faces_desc() const {
    return get_edges_desc();
}

size_t Quadrangle::num_nodes() const {
    return 4;
}
ShapeType Quadrangle::Shape() const {
    return ShapeType::Quadrangle4;
}
tElementDescList Quadrangle::get_edges_desc() const {
    return { get_part(ShapeType::Segment2, 0, 1),
             get_part(ShapeType::Segment2, 1, 2),
             get_part(ShapeType::Segment2, 2, 3),
             get_part(ShapeType::Segment2, 3, 0) };
}
tElementDescList Quadrangle::get_faces_desc() const {
    return get_edges_desc();
}
tElementDescList Quadrangle::get_simplicial_parts(size_t partition_index) const {
    switch (partition_index) {
        case 0:
            return { get_part(ShapeType::Triangle3, 0, 1, 2),
                     get_part(ShapeType::Triangle3, 2, 3, 0) };
        case 1:
            return { get_part(ShapeType::Triangle3, 0, 1, 3),
                     get_part(ShapeType::Triangle3, 1, 2, 3) };
        default:
            return {};
    }
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

real_t Tetrahedron::get_diameter() const {
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
real_t Tetrahedron::LenAreaOrVolume() const {
    const vec3_t delta1 =
        get_node_coords(1) - get_node_coords(0);
    const vec3_t delta2 =
        get_node_coords(2) - get_node_coords(0);
    const vec3_t delta3 =
        get_node_coords(3) - get_node_coords(0);
    return std::abs(glm::dot(delta1,
                             glm::cross(delta2, delta3)))/6.0;
}
vec3_t Tetrahedron::CenterPos() const {
    return 0.25*(get_node_coords(0) + get_node_coords(1) +
                 get_node_coords(2) + get_node_coords(3));
}
ShapeType Tetrahedron::Shape() const {
    return ShapeType::Tetrahedron4;
}
size_t Tetrahedron::num_nodes() const {
    return 4;
}
tElementDescList Tetrahedron::get_edges_desc() const {
    return { get_part(ShapeType::Segment2, 0, 1),
             get_part(ShapeType::Segment2, 1, 2),
             get_part(ShapeType::Segment2, 2, 0),
             get_part(ShapeType::Segment2, 0, 3),
             get_part(ShapeType::Segment2, 1, 3),
             get_part(ShapeType::Segment2, 2, 3) };
}
tElementDescList Tetrahedron::get_faces_desc() const {
    return { get_part(ShapeType::Triangle3, 0, 2, 1),
             get_part(ShapeType::Triangle3, 0, 1, 3),
             get_part(ShapeType::Triangle3, 1, 2, 3),
             get_part(ShapeType::Triangle3, 2, 0, 3) };
}

size_t Pyramid::num_nodes() const {
    return 5;
}
ShapeType Pyramid::Shape() const {
    return ShapeType::Pyramid5;
}
tElementDescList Pyramid::get_edges_desc() const {
    return { get_part(ShapeType::Segment2, 0, 1),
             get_part(ShapeType::Segment2, 1, 2),
             get_part(ShapeType::Segment2, 2, 3),
             get_part(ShapeType::Segment2, 3, 0),
             get_part(ShapeType::Segment2, 0, 4),
             get_part(ShapeType::Segment2, 1, 4),
             get_part(ShapeType::Segment2, 2, 4),
             get_part(ShapeType::Segment2, 3, 4) };
}
tElementDescList Pyramid::get_faces_desc() const {
    return { get_part(ShapeType::Quadrangle4, 0, 3, 2, 1),
             get_part(ShapeType::Triangle3, 0, 1, 4),
             get_part(ShapeType::Triangle3, 1, 2, 4),
             get_part(ShapeType::Triangle3, 2, 3, 4),
             get_part(ShapeType::Triangle3, 3, 0, 4) };
}
tElementDescList Pyramid::get_simplicial_parts(size_t partition_index) const {
    switch (partition_index) {
    case 0:
        return { get_part(ShapeType::Tetrahedron4, 0, 1, 2, 4),
                 get_part(ShapeType::Tetrahedron4, 2, 3, 0, 4) };
    case 1:
        return { get_part(ShapeType::Tetrahedron4, 0, 1, 3, 4),
                 get_part(ShapeType::Tetrahedron4, 1, 2, 3, 4) };
    default:
        return {};
    }
}

size_t Pentahedron::num_nodes() const {
    return 6;
}
ShapeType Pentahedron::Shape() const {
    return ShapeType::Pentahedron6;
}
tElementDescList Pentahedron::get_edges_desc() const {
    return { get_part(ShapeType::Segment2, 0, 1),
             get_part(ShapeType::Segment2, 1, 2),
             get_part(ShapeType::Segment2, 2, 0),
             get_part(ShapeType::Segment2, 0, 3),
             get_part(ShapeType::Segment2, 1, 4),
             get_part(ShapeType::Segment2, 2, 5),
             get_part(ShapeType::Segment2, 3, 4),
             get_part(ShapeType::Segment2, 4, 5),
             get_part(ShapeType::Segment2, 5, 3) };
}
tElementDescList Pentahedron::get_faces_desc() const {
    return { get_part(ShapeType::Quadrangle4, 0, 1, 4, 3),
             get_part(ShapeType::Quadrangle4, 1, 2, 5, 4),
             get_part(ShapeType::Quadrangle4, 2, 0, 3, 5),
             get_part(ShapeType::Triangle3, 0, 2, 1),
             get_part(ShapeType::Triangle3, 3, 4, 5) };
}
tElementDescList Pentahedron::get_simplicial_parts(size_t partition_index) const {
    FEATHERS_NOT_IMPLEMENTED();
}

size_t Hexahedron::num_nodes() const {
    return 8;
}
ShapeType Hexahedron::Shape() const {
    return ShapeType::Hexahedron8;
}
tElementDescList Hexahedron::get_edges_desc() const {
    return { get_part(ShapeType::Segment2, 0, 1),
             get_part(ShapeType::Segment2, 1, 2),
             get_part(ShapeType::Segment2, 2, 3),
             get_part(ShapeType::Segment2, 3, 0),
             get_part(ShapeType::Segment2, 0, 4),
             get_part(ShapeType::Segment2, 1, 5),
             get_part(ShapeType::Segment2, 2, 6),
             get_part(ShapeType::Segment2, 3, 7),
             get_part(ShapeType::Segment2, 4, 5),
             get_part(ShapeType::Segment2, 5, 6),
             get_part(ShapeType::Segment2, 6, 7),
             get_part(ShapeType::Segment2, 7, 4) };
}
tElementDescList Hexahedron::get_faces_desc() const {
    return { get_part(ShapeType::Quadrangle4, 0, 3, 2, 1),
             get_part(ShapeType::Quadrangle4, 0, 1, 5, 4),
             get_part(ShapeType::Quadrangle4, 1, 2, 6, 5),
             get_part(ShapeType::Quadrangle4, 2, 3, 7, 6),
             get_part(ShapeType::Quadrangle4, 0, 4, 7, 3),
             get_part(ShapeType::Quadrangle4, 4, 5, 6, 7) };
}
tElementDescList Hexahedron::get_simplicial_parts(size_t partition_index) const {
    // TODO: alternative partitions!
    switch (partition_index) {
        case 0:
            return { get_part(ShapeType::Tetrahedron4, 0, 3, 1, 4),
                     get_part(ShapeType::Tetrahedron4, 3, 2, 1, 6),
                     get_part(ShapeType::Tetrahedron4, 4, 5, 6, 1),
                     get_part(ShapeType::Tetrahedron4, 4, 6, 7, 3),
                     get_part(ShapeType::Tetrahedron4, 4, 3, 1, 6) };
        default:
            return {};
    }
}

} // namespace feathers

#if !FEATHERS_DOXYGEN
/**
 * Alternative pyramid pictures.
 * @verbatim
 *                            n4
 *                            O    f1
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
 *   n3 O----<----,/----------O n2   |
 *       \.  e2   |`           \.    \.
 *        `\.    ,|          e1 `^.  `\
 *      e3 `v.   /`     o         `\. |
 *           `\ /       |            `\|
 *         n0 O------------------->---O n1
 *                      |        e0
 *                      |
 *                      v
 *                      f0
 * @endverbatim
 *
 * @verbatim
 *                           f2
 *                     n4    ^
 *                     O.   /
 *                   ,/|`\./
 *              e7 ,^/ |  /\.
 *               ,/ /  | *  ^ e6
 *             ,/`,/   ^ e5 `\.
 *           ,/` /`    |      `\
 *       n3 O---/------|---<---O n2
 *         /   ^ e4    |  e2 ,/`
 *        /  ,/ *      |   ,/`
 *    e3 v /   / o     | ,^` e1
 *      /,/   /  |     |/`
 *  n0 O-----/----->---O n1
 *          /    | e0
 *         /     |
 *        v      v
 *       f1      f0
 * @endverbatim
 */
#endif
