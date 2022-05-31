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

std::unique_ptr<Element> Element::Make(ShapeDesc&& desc,
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

  // Verify and assign nodes:
  element->NodePos_ = nodePos;
  element->NodeIndices_ = std::move(desc.NodeIndices);
  storm_assert(
    // we have the same number of node indices and nodes,
    element->NodeIndices_.size() == element->NumNodes() &&
    // nodes span has enough elements,
    element->NodeIndices_.size() <= element->NodePos_.size() &&
    // and all node indices are inside the range.
    ranges::all_of(element->NodeIndices_, [&](size_t nodeIndex) {
      return nodeIndex < element->NodePos_.size();
    }));

  return element;

} // Element::Make

real_t SimplexElement::Diam() const {
  if (NumNodes() == 1) {
    return 0.0;
  }
  real_t diam{glm::length(NodePos(0) - NodePos(1))};
  for (size_t i{2}; i < NumNodes(); ++i) {
    for (size_t j{0}; j < i; ++j) {
      diam = std::max(diam, glm::length(NodePos(i) - NodePos(j)));
    }
  }
  return diam;
} // SimplexElement::Diam

vec3_t SimplexElement::CenterPos() const {
  vec3_t sumOfNodePos{NodePos(0)};
  for (size_t i{1}; i < NumNodes(); ++i) {
    sumOfNodePos += NodePos(i);
  }
  return sumOfNodePos/real_t(NumNodes());
} // SimplexElement::centerPos

template<class Func>
void ComplexElement::ForEachSimplex_(Func&& func) const {
  ShapeDescArray simplicesDesc{MakeSimplicesDesc()};
  for (ShapeDesc& simplexDesc : simplicesDesc) {
    auto const simplex = Element::Make(std::move(simplexDesc), NodePos_);
    func(*simplex);
  }
} // ComplexElement::ForEachSimplex_

real_t ComplexElement::Diam() const {
  real_t diam{0.0};
  ForEachSimplex_([&](Element& shape) {
    diam = std::max(diam, shape.Diam());
  });
  return diam;
} // ComplexElement::Diam

real_t ComplexElement::Volume() const {
  real_t volume{0.0};
  ForEachSimplex_([&](Element& shape) {
    volume += shape.Volume();
  });
  return volume;
} // ComplexElement::volume

vec3_t ComplexElement::Normal() const {
  vec3_t weightedSumOfNormals(0.0);
  ForEachSimplex_([&](Element& shape) {
    weightedSumOfNormals += shape.Volume() * shape.Normal();
  });
  return glm::normalize(weightedSumOfNormals);
} // ComplexElement::normal

vec3_t ComplexElement::CenterPos() const {
  vec3_t weightedSumOfCenterPos(0.0);
  real_t volume{0.0};
  ForEachSimplex_([&](Element& shape) {
    real_t const partVolume{shape.Volume()};
    weightedSumOfCenterPos += partVolume * shape.CenterPos();
    volume += partVolume;
  });
  return weightedSumOfCenterPos / volume;
} // ComplexElement::centerPos

real_t Node::Volume() const {
  return 1.0;
}
vec3_t Node::Normal() const {
  static constexpr vec3_t right(1.0, 0.0, 0.0);
  return right;
}
ShapeDescArray Node::MakeEdgesDesc() const {
  return {};
}
ShapeDescArray Node::MakeFacesDesc() const {
  return {};
}

real_t Segment::Volume() const {
  return Diam();
}
vec3_t Segment::Normal() const {
  vec3_t const delta = NodePos(1) - NodePos(0);
  static constexpr vec3_t up(0.0, 0.0, 1.0);
  return glm::normalize(glm::cross(delta, up));
}
vec3_t Segment::Dir() const {
  vec3_t const delta = NodePos(1) - NodePos(0);
  return glm::normalize(delta);
}
ShapeDescArray Segment::MakeEdgesDesc() const {
  return {PartDesc_(ShapeType::Segment2, 0, 1)};
}
ShapeDescArray Segment::MakeFacesDesc() const {
  return {PartDesc_(ShapeType::Node, 0), PartDesc_(ShapeType::Node, 1)};
}

real_t Triangle::Volume() const {
  vec3_t const delta1 = NodePos(1) - NodePos(0);
  vec3_t const delta2 = NodePos(2) - NodePos(0);
  return 0.5*glm::length(glm::cross(delta1, delta2));
}
vec3_t Triangle::Normal() const {
  vec3_t const delta1 = NodePos(1) - NodePos(0);
  vec3_t const delta2 = NodePos(2) - NodePos(0);
  return glm::normalize(glm::cross(delta1, delta2));
}
ShapeDescArray Triangle::MakeEdgesDesc() const {
  return {
    PartDesc_(ShapeType::Segment2, 0, 1),
    PartDesc_(ShapeType::Segment2, 1, 2),
    PartDesc_(ShapeType::Segment2, 2, 0)};
}
ShapeDescArray Triangle::MakeFacesDesc() const {
  return MakeEdgesDesc();
}

ShapeDescArray Quadrangle::MakeEdgesDesc() const {
  return {
    PartDesc_(ShapeType::Segment2, 0, 1),
    PartDesc_(ShapeType::Segment2, 1, 2),
    PartDesc_(ShapeType::Segment2, 2, 3),
    PartDesc_(ShapeType::Segment2, 3, 0)};
}
ShapeDescArray Quadrangle::MakeFacesDesc() const {
  return MakeEdgesDesc();
}
ShapeDescArray Quadrangle::MakeSimplicesDesc() const {
  return {
    PartDesc_(ShapeType::Triangle3, 0, 1, 2),
    PartDesc_(ShapeType::Triangle3, 2, 3, 0)};
//return {
//  PartDesc_(ShapeType::Triangle3, 0, 1, 3),
//  PartDesc_(ShapeType::Triangle3, 1, 2, 3)};
}

real_t Tetrahedron::Volume() const {
  vec3_t const delta1 = NodePos(1) - NodePos(0);
  vec3_t const delta2 = NodePos(2) - NodePos(0);
  vec3_t const delta3 = NodePos(3) - NodePos(0);
  return std::abs(glm::dot(delta1, glm::cross(delta2, delta3)))/6.0;
}
ShapeDescArray Tetrahedron::MakeEdgesDesc() const {
  return {
    PartDesc_(ShapeType::Segment2, 0, 1),
    PartDesc_(ShapeType::Segment2, 1, 2),
    PartDesc_(ShapeType::Segment2, 2, 0),
    PartDesc_(ShapeType::Segment2, 0, 3),
    PartDesc_(ShapeType::Segment2, 1, 3),
    PartDesc_(ShapeType::Segment2, 2, 3) };
}
ShapeDescArray Tetrahedron::MakeFacesDesc() const {
  return {
    PartDesc_(ShapeType::Triangle3, 0, 2, 1),
    PartDesc_(ShapeType::Triangle3, 0, 1, 3),
    PartDesc_(ShapeType::Triangle3, 1, 2, 3),
    PartDesc_(ShapeType::Triangle3, 2, 0, 3) };
}

ShapeDescArray Pyramid::MakeEdgesDesc() const {
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
ShapeDescArray Pyramid::MakeFacesDesc() const {
  return {
    PartDesc_(ShapeType::Quadrangle4, 0, 3, 2, 1),
    PartDesc_(ShapeType::Triangle3, 0, 1, 4),
    PartDesc_(ShapeType::Triangle3, 1, 2, 4),
    PartDesc_(ShapeType::Triangle3, 2, 3, 4),
    PartDesc_(ShapeType::Triangle3, 3, 0, 4)};
}
ShapeDescArray Pyramid::MakeSimplicesDesc() const {
  return {
    PartDesc_(ShapeType::Tetrahedron4, 0, 1, 2, 4),
    PartDesc_(ShapeType::Tetrahedron4, 2, 3, 0, 4)};
//return {
//  PartDesc_(ShapeType::Tetrahedron4, 0, 1, 3, 4),
//  PartDesc_(ShapeType::Tetrahedron4, 1, 2, 3, 4)};
}

ShapeDescArray Pentahedron::MakeEdgesDesc() const {
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
ShapeDescArray Pentahedron::MakeFacesDesc() const {
  return {
    PartDesc_(ShapeType::Quadrangle4, 0, 1, 4, 3),
    PartDesc_(ShapeType::Quadrangle4, 1, 2, 5, 4),
    PartDesc_(ShapeType::Quadrangle4, 2, 0, 3, 5),
    PartDesc_(ShapeType::Triangle3, 0, 2, 1),
    PartDesc_(ShapeType::Triangle3, 3, 4, 5)};
}
ShapeDescArray Pentahedron::MakeSimplicesDesc() const {
  return {
    PartDesc_(ShapeType::Tetrahedron4, 0, 1, 2, 4),
    PartDesc_(ShapeType::Tetrahedron4, 2, 0, 3, 4),
    PartDesc_(ShapeType::Tetrahedron4, 3, 5, 2, 4)};
}

ShapeDescArray Hexahedron::MakeEdgesDesc() const {
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
ShapeDescArray Hexahedron::MakeFacesDesc() const {
  return {
    PartDesc_(ShapeType::Quadrangle4, 0, 3, 2, 1),
    PartDesc_(ShapeType::Quadrangle4, 0, 1, 5, 4),
    PartDesc_(ShapeType::Quadrangle4, 1, 2, 6, 5),
    PartDesc_(ShapeType::Quadrangle4, 2, 3, 7, 6),
    PartDesc_(ShapeType::Quadrangle4, 0, 4, 7, 3),
    PartDesc_(ShapeType::Quadrangle4, 4, 5, 6, 7)};
}
ShapeDescArray Hexahedron::MakeSimplicesDesc() const {
  return {
    PartDesc_(ShapeType::Tetrahedron4, 0, 3, 1, 4),
    PartDesc_(ShapeType::Tetrahedron4, 3, 2, 1, 6),
    PartDesc_(ShapeType::Tetrahedron4, 4, 5, 6, 1),
    PartDesc_(ShapeType::Tetrahedron4, 4, 6, 7, 3),
    PartDesc_(ShapeType::Tetrahedron4, 4, 3, 1, 6)};
}

} // namespace feathers

#if !FEATHERS_DOXYGEN
///
/// Alternative pyramid pictures.
/// @verbatim
///                            n4
///                            O    f1
///                           /|\   ^
///                          /`|`\./
///                        ,// | `/
///                       /,/  | /\.
///                     ,/ |`  |/ `\
///          f3     e7 ^` ,|   /   ^ e5
///           ^      ,/   /`  /|   \.
///            \    /`   /   * |   `\
///             \ ,/   ,/      |    |
///              /`    |`   e6 ^    \.
///            ,/ \    ^ e4    |  o-`|---------> f2
///           /`   o  ,|       |     |
///   <----------*    /`       |     \.
///       ,/`        /         |     `\
///   n3 O----<----,/----------O n2   |
///       \.  e2   |`           \.    \.
///        `\.    ,|          e1 `^.  `\
///      e3 `v.   /`     o         `\. |
///           `\ /       |            `\|
///         n0 O------------------->---O n1
///                      |        e0
///                      |
///                      v
///                      f0
/// @endverbatim
///
/// @verbatim
///                           f2
///                     n4    ^
///                     O.   /
///                   ,/|`\./
///              e7 ,^/ |  /\.
///               ,/ /  | *  ^ e6
///             ,/`,/   ^ e5 `\.
///           ,/` /`    |      `\
///       n3 O---/------|---<---O n2
///         /   ^ e4    |  e2 ,/`
///        /  ,/ *      |   ,/`
///    e3 v /   / o     | ,^` e1
///      /,/   /  |     |/`
///  n0 O-----/----->---O n1
///          /    | e0
///         /     |
///        v      v
///       f1      f0
/// @endverbatim
////
#endif
