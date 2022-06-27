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

namespace Storm {

std::unique_ptr<Shape> Shape::make(const ShapeDesc& shape_desc) {
  auto shape = [shape_type =
                    shape_desc.shape_type]() -> std::unique_ptr<Shape> {
    switch (shape_type) {
      case ShapeType::Node: return std::make_unique<Node>();
      case ShapeType::Segment: return std::make_unique<Segment>();
      case ShapeType::Triangle: return std::make_unique<Triangle>();
      case ShapeType::Quadrangle: return std::make_unique<Quadrangle>();
      case ShapeType::Tetrahedron: return std::make_unique<Tetrahedron>();
      case ShapeType::Pyramid: return std::make_unique<Pyramid>();
      case ShapeType::Pentahedron: return std::make_unique<Pentahedron>();
      case ShapeType::Hexahedron: return std::make_unique<Hexahedron>();
    }
    FEATHERS_NOT_REACHABLE();
  }();
  shape->node_indices_ = shape_desc.node_indices;
  return shape;
}

} // namespace Storm

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
