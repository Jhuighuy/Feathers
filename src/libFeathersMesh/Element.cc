/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person
/// obtaining a copy of this software and associated documentation
/// files (the "Software"), to deal in the Software without
/// restriction, including without limitation the rights  to use,
/// copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the
/// Software is furnished to do so, subject to the following
/// conditions:
///
/// The above copyright notice and this permission notice shall be
/// included in all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
/// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
/// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
/// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
/// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
/// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
/// OTHER DEALINGS IN THE SOFTWARE.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///

#include "Element.hh"

namespace feathers {

/**
 * Construct the element.
 */
static std::unique_ptr<Element> construct_element_(ShapeType shape) {

}   // construct_element_

/**
 * Construct a new element object.
 */
std::unique_ptr<Element> Element::make(ElementDesc&& desc,
                                       std::span<vec3_t const> nodePos) {

  // Construct an element by shape.
  auto element = [shape = desc.Shape]() -> std::unique_ptr<Element> {
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
  }();

  // Verify and assign nodes. */
  element->NodePos_ = nodePos;
  element->NodeIndices_ = std::move(desc.NodeIndices);
  StormAssert(
    desc.NodeIndices.size() == element->NumNodes() &&
    desc.NodeIndices.size() <= element->NodePos_.size() &&
    ranges::all_of(element->NodeIndices_, [&](size_t nodeIndex) {
      return nodeIndex < element->NodePos_.size();
    }));

  return element;

} // Element::make

template<class Func>
void ComplexElement::ForEachSimplex_(Func&& func) const {
  ElementDescArray simplexParts(MakeSimplicialPartsDesc());
  for (ElementDesc& part : simplexParts) {
    auto const simplex = Element::make(std::move(part), NodePos_);
    func(*simplex);
  }
} // ComplexElement::ForEachSimplex_

real_t ComplexElement::Diam() const {
  real_t diam = 0.0;
  ForEachSimplex_([&](Element& shape) {
    diam = std::max(diam, shape.Diam());
  });
  return diam;
}   // ComplexElement::Diam

real_t ComplexElement::LenAreaOrVolume() const {
  real_t lenAreaOrVolume = 0.0;
  ForEachSimplex_([&](Element& shape) {
    lenAreaOrVolume += shape.LenAreaOrVolume();
  });
  return lenAreaOrVolume;
} // ComplexElement::LenAreaOrVolume

vec3_t ComplexElement::Normal() const {
  vec3_t weightedSumOfNormals(0.0);
  ForEachSimplex_([&](Element& shape) {
    weightedSumOfNormals +=
      shape.LenAreaOrVolume() * shape.Normal();
  });
  return glm::normalize(weightedSumOfNormals);
} // ComplexElement::Normal

vec3_t ComplexElement::CenterPos() const {
  vec3_t weightedSumOfCenterPos(0.0);
  real_t lenAreaOrVolume{0.0};
  ForEachSimplex_([&](Element& shape) {
    real_t const partLenAreaOrVolume = shape.LenAreaOrVolume();
    weightedSumOfCenterPos += partLenAreaOrVolume*shape.CenterPos();
    lenAreaOrVolume += partLenAreaOrVolume;
  });
  return weightedSumOfCenterPos/lenAreaOrVolume;
} // ComplexElement::CenterPos

real_t Node::Diam() const {
  return 0.0;
}
real_t Node::LenAreaOrVolume() const {
  return 1.0;
}
vec3_t Node::Normal() const {
  static constexpr vec3_t right(1.0, 0.0, 0.0);
  return right;
}
vec3_t Node::CenterPos() const {
  return NodePos(0);
}
ShapeType Node::Shape() const noexcept {
  return ShapeType::Node;
}
size_t Node::NumNodes() const noexcept {
  return 1;
}
ElementDescArray Node::MakeEdgesDesc() const {
  return {};
}
ElementDescArray Node::MakeFacesDesc() const {
  return {};
}

real_t Segment::Diam() const {
  const vec3_t delta = NodePos(1) - NodePos(0);
  return glm::length(delta);
}
real_t Segment::LenAreaOrVolume() const {
  return Diam();
}
vec3_t Segment::Normal() const {
  const vec3_t delta = NodePos(1) - NodePos(0);
  static constexpr vec3_t up(0.0, 0.0, 1.0);
  return glm::normalize(glm::cross(delta, up));
}
vec3_t Segment::Dir() const {
  const vec3_t delta = NodePos(1) - NodePos(0);
  return glm::normalize(delta);
}
vec3_t Segment::CenterPos() const {
  return 0.5*(NodePos(0) + NodePos(1));
}
ShapeType Segment::Shape() const noexcept {
  return ShapeType::Segment2;
}
size_t Segment::NumNodes() const noexcept {
  return 2;
}
ElementDescArray Segment::MakeEdgesDesc() const {
  return {PartDesc_(ShapeType::Segment2, 0, 1) };
}
ElementDescArray Segment::MakeFacesDesc() const {
  return {PartDesc_(ShapeType::Node, 0), PartDesc_(ShapeType::Node, 1) };
}

real_t Triangle::Diam() const {
  const vec3_t delta1 = NodePos(1) - NodePos(0);
  const vec3_t delta2 = NodePos(2) - NodePos(1);
  const vec3_t delta3 = NodePos(0) - NodePos(2);
  return std::max(glm::length(delta1),
                  std::max(glm::length(delta2), glm::length(delta3)));
}
real_t Triangle::LenAreaOrVolume() const {
  const vec3_t delta1 = NodePos(1) - NodePos(0);
  const vec3_t delta2 = NodePos(2) - NodePos(0);
  return 0.5*glm::length(glm::cross(delta1, delta2));
}
vec3_t Triangle::Normal() const {
  const vec3_t delta1 = NodePos(1) - NodePos(0);
  const vec3_t delta2 = NodePos(2) - NodePos(0);
  return glm::normalize(glm::cross(delta1, delta2));
}
vec3_t Triangle::CenterPos() const {
  return (NodePos(0) + NodePos(1) + NodePos(2)) / 3.0;
}
ShapeType Triangle::Shape() const noexcept {
  return ShapeType::Triangle3;
}
size_t Triangle::NumNodes() const noexcept {
  return 3;
}
ElementDescArray Triangle::MakeEdgesDesc() const {
  return {
    PartDesc_(ShapeType::Segment2, 0, 1),
    PartDesc_(ShapeType::Segment2, 1, 2),
    PartDesc_(ShapeType::Segment2, 2, 0)};
}
ElementDescArray Triangle::MakeFacesDesc() const {
  return MakeEdgesDesc();
}

size_t Quadrangle::NumNodes() const noexcept {
  return 4;
}
ShapeType Quadrangle::Shape() const noexcept {
  return ShapeType::Quadrangle4;
}
ElementDescArray Quadrangle::MakeEdgesDesc() const {
  return {
    PartDesc_(ShapeType::Segment2, 0, 1),
    PartDesc_(ShapeType::Segment2, 1, 2),
    PartDesc_(ShapeType::Segment2, 2, 3),
    PartDesc_(ShapeType::Segment2, 3, 0)};
}
ElementDescArray Quadrangle::MakeFacesDesc() const {
  return MakeEdgesDesc();
}
ElementDescArray Quadrangle::MakeSimplicialPartsDesc() const {
  return {
    PartDesc_(ShapeType::Triangle3, 0, 1, 2),
    PartDesc_(ShapeType::Triangle3, 2, 3, 0)};
  //return {
  //  PartDesc_(ShapeType::Triangle3, 0, 1, 3),
  //  PartDesc_(ShapeType::Triangle3, 1, 2, 3)};
}

real_t Tetrahedron::Diam() const {
  const vec3_t delta1 = NodePos(1) - NodePos(0);
  const vec3_t delta2 = NodePos(2) - NodePos(1);
  const vec3_t delta3 = NodePos(0) - NodePos(2);
  const vec3_t delta4 = NodePos(3) - NodePos(0);
  const vec3_t delta5 = NodePos(3) - NodePos(1);
  const vec3_t delta6 = NodePos(3) - NodePos(2);
  return std::max({glm::length(delta1), glm::length(delta2),
                   glm::length(delta3), glm::length(delta4),
                   glm::length(delta5), glm::length(delta6)});
}
real_t Tetrahedron::LenAreaOrVolume() const {
  const vec3_t delta1 = NodePos(1) - NodePos(0);
  const vec3_t delta2 = NodePos(2) - NodePos(0);
  const vec3_t delta3 = NodePos(3) - NodePos(0);
  return std::abs(glm::dot(delta1, glm::cross(delta2, delta3)))/6.0;
}
vec3_t Tetrahedron::CenterPos() const {
  return 0.25*(NodePos(0) + NodePos(1) +
               NodePos(2) + NodePos(3));
}
ShapeType Tetrahedron::Shape() const noexcept {
  return ShapeType::Tetrahedron4;
}
size_t Tetrahedron::NumNodes() const noexcept {
  return 4;
}
ElementDescArray Tetrahedron::MakeEdgesDesc() const {
  return {
    PartDesc_(ShapeType::Segment2, 0, 1),
    PartDesc_(ShapeType::Segment2, 1, 2),
    PartDesc_(ShapeType::Segment2, 2, 0),
    PartDesc_(ShapeType::Segment2, 0, 3),
    PartDesc_(ShapeType::Segment2, 1, 3),
    PartDesc_(ShapeType::Segment2, 2, 3) };
}
ElementDescArray Tetrahedron::MakeFacesDesc() const {
  return {
    PartDesc_(ShapeType::Triangle3, 0, 2, 1),
    PartDesc_(ShapeType::Triangle3, 0, 1, 3),
    PartDesc_(ShapeType::Triangle3, 1, 2, 3),
    PartDesc_(ShapeType::Triangle3, 2, 0, 3) };
}

size_t Pyramid::NumNodes() const noexcept {
  return 5;
}
ShapeType Pyramid::Shape() const noexcept {
  return ShapeType::Pyramid5;
}
ElementDescArray Pyramid::MakeEdgesDesc() const {
  return {
    PartDesc_(ShapeType::Segment2, 0, 1),
    PartDesc_(ShapeType::Segment2, 1, 2),
    PartDesc_(ShapeType::Segment2, 2, 3),
    PartDesc_(ShapeType::Segment2, 3, 0),
    PartDesc_(ShapeType::Segment2, 0, 4),
    PartDesc_(ShapeType::Segment2, 1, 4),
    PartDesc_(ShapeType::Segment2, 2, 4),
    PartDesc_(ShapeType::Segment2, 3, 4)};
}
ElementDescArray Pyramid::MakeFacesDesc() const {
  return {
    PartDesc_(ShapeType::Quadrangle4, 0, 3, 2, 1),
    PartDesc_(ShapeType::Triangle3, 0, 1, 4),
    PartDesc_(ShapeType::Triangle3, 1, 2, 4),
    PartDesc_(ShapeType::Triangle3, 2, 3, 4),
    PartDesc_(ShapeType::Triangle3, 3, 0, 4)};
}
ElementDescArray Pyramid::MakeSimplicialPartsDesc() const {
  return {
    PartDesc_(ShapeType::Tetrahedron4, 0, 1, 2, 4),
    PartDesc_(ShapeType::Tetrahedron4, 2, 3, 0, 4)};
  //return {
  //  PartDesc_(ShapeType::Tetrahedron4, 0, 1, 3, 4),
  //  PartDesc_(ShapeType::Tetrahedron4, 1, 2, 3, 4)};
}

size_t Pentahedron::NumNodes() const noexcept {
  return 6;
}
ShapeType Pentahedron::Shape() const noexcept {
  return ShapeType::Pentahedron6;
}
ElementDescArray Pentahedron::MakeEdgesDesc() const {
  return {
    PartDesc_(ShapeType::Segment2, 0, 1),
    PartDesc_(ShapeType::Segment2, 1, 2),
    PartDesc_(ShapeType::Segment2, 2, 0),
    PartDesc_(ShapeType::Segment2, 0, 3),
    PartDesc_(ShapeType::Segment2, 1, 4),
    PartDesc_(ShapeType::Segment2, 2, 5),
    PartDesc_(ShapeType::Segment2, 3, 4),
    PartDesc_(ShapeType::Segment2, 4, 5),
    PartDesc_(ShapeType::Segment2, 5, 3)};
}
ElementDescArray Pentahedron::MakeFacesDesc() const {
  return {
    PartDesc_(ShapeType::Quadrangle4, 0, 1, 4, 3),
    PartDesc_(ShapeType::Quadrangle4, 1, 2, 5, 4),
    PartDesc_(ShapeType::Quadrangle4, 2, 0, 3, 5),
    PartDesc_(ShapeType::Triangle3, 0, 2, 1),
    PartDesc_(ShapeType::Triangle3, 3, 4, 5)};
}
ElementDescArray Pentahedron::MakeSimplicialPartsDesc() const {
  FEATHERS_NOT_IMPLEMENTED();
}

size_t Hexahedron::NumNodes() const noexcept {
  return 8;
}
ShapeType Hexahedron::Shape() const noexcept {
  return ShapeType::Hexahedron8;
}
ElementDescArray Hexahedron::MakeEdgesDesc() const {
  return {
    PartDesc_(ShapeType::Segment2, 0, 1),
    PartDesc_(ShapeType::Segment2, 1, 2),
    PartDesc_(ShapeType::Segment2, 2, 3),
    PartDesc_(ShapeType::Segment2, 3, 0),
    PartDesc_(ShapeType::Segment2, 0, 4),
    PartDesc_(ShapeType::Segment2, 1, 5),
    PartDesc_(ShapeType::Segment2, 2, 6),
    PartDesc_(ShapeType::Segment2, 3, 7),
    PartDesc_(ShapeType::Segment2, 4, 5),
    PartDesc_(ShapeType::Segment2, 5, 6),
    PartDesc_(ShapeType::Segment2, 6, 7),
    PartDesc_(ShapeType::Segment2, 7, 4)};
}
ElementDescArray Hexahedron::MakeFacesDesc() const {
  return {
    PartDesc_(ShapeType::Quadrangle4, 0, 3, 2, 1),
    PartDesc_(ShapeType::Quadrangle4, 0, 1, 5, 4),
    PartDesc_(ShapeType::Quadrangle4, 1, 2, 6, 5),
    PartDesc_(ShapeType::Quadrangle4, 2, 3, 7, 6),
    PartDesc_(ShapeType::Quadrangle4, 0, 4, 7, 3),
    PartDesc_(ShapeType::Quadrangle4, 4, 5, 6, 7)};
}
ElementDescArray Hexahedron::MakeSimplicialPartsDesc() const {
  return {
    PartDesc_(ShapeType::Tetrahedron4, 0, 3, 1, 4),
    PartDesc_(ShapeType::Tetrahedron4, 3, 2, 1, 6),
    PartDesc_(ShapeType::Tetrahedron4, 4, 5, 6, 1),
    PartDesc_(ShapeType::Tetrahedron4, 4, 6, 7, 3),
    PartDesc_(ShapeType::Tetrahedron4, 4, 3, 1, 6)};
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
