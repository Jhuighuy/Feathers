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

#pragma once

#include <span>

#include <stormMesh/Base.hxx>

namespace Storm {

/// @brief Shape type.
enum class ShapeType : byte_t {
  Null,
  Node,
  Segment,
  Triangle,
  Quadrangle,
  Tetrahedron,
  Pyramid,
  Pentahedron,
  Hexahedron,
}; // enum class ShapeType

/// @brief Shape description.
struct ShapeDesc {
  ShapeType Shape;
  std::vector<NodeIndex> NodeIndices;
}; // ShapeDesc

/// @brief Array of the shape descriptions.
using ShapeDescArray = std::vector<ShapeDesc>;

/// ----------------------------------------------------------------- ///
/// @brief Abstract element class.
/// ----------------------------------------------------------------- ///
class Element : public NonCopyable {
protected:
  std::span<vec3_t const> NodePos_;
  std::vector<NodeIndex> NodeIndices_;

  template<class... Indices>
  auto PartDesc_(ShapeType partShape, Indices... nodeLocals) const {
    return ShapeDesc{partShape, {NodeIndices_[nodeLocals]...}};
  }

  Element() = default;

public:
  /// @brief Virtual destructor.
  virtual ~Element() = default;

  /// @brief Construct a new element object \
  ///   with a description @p desc and a node position array @p NodeCoords.
  static std::unique_ptr<Element> Make(ShapeDesc&& desc,
                                       std::span<vec3_t const> nodePos);

  /** Get element node indices. */
  std::vector<NodeIndex> const& NodeIndices() const {
    return NodeIndices_;
  }

  /// @brief Get element shape.
  virtual ShapeType Shape() const noexcept = 0;

  /// @brief Get node @p position.
  vec3_t NodePos(size_t nodeLocal) const {
    StormAssert(nodeLocal < NodeIndices_.size());
    return NodePos_[(size_t) NodeIndices_[nodeLocal]];
  }

  /// @brief Compute the element diameter
  virtual real_t Diam() const {
    return qnan;
  }

  /// @brief Compute the element Volume (Area or length).
  virtual real_t Volume() const {
    return qnan;
  }

  /// @brief Compute the Normal to element.
  virtual vec3_t Normal() const {
    return vec3_t(qnan);
  }

  /// @brief Compute the element direction.
  virtual vec3_t Dir() const {
    return vec3_t(qnan);
  }

  /// @brief Compute the element center position.
  virtual vec3_t CenterPos() const {
    return vec3_t(qnan);
  }

  /// @brief Number of Nodes in the element.
  virtual size_t NumNodes() const noexcept = 0;

  /// @brief Number of Edges in the element.
  virtual size_t NumEdges() const {
    return make_edges_desc().size();
  }

  /// @brief Make element Edges description array.
  virtual ShapeDescArray make_edges_desc() const = 0;

  /// @brief Number of Faces in the element.
  size_t NumFaces() const {
    return make_faces_desc().size();
  }

  /// @brief Make element Faces description.
  virtual ShapeDescArray make_faces_desc() const = 0;

}; // class Element

using Shape = Element;

template<ShapeType Shape_, size_t NumNodes_, class Base>
class ElementHelper_ : public Base {
public:
  ShapeType Shape() const noexcept final {
    return Shape_;
  }
  size_t NumNodes() const noexcept final {
    return NumNodes_;
  }
}; // class ElementHelper_<...>

/// ----------------------------------------------------------------- ///
/// @brief Abstract simplex element class.
/// ----------------------------------------------------------------- ///
class SimplexElement : public Element {
public:
  real_t Diam() const final;
  vec3_t CenterPos() const final;
}; // SimplexElement

/// ----------------------------------------------------------------- ///
/// @brief Abstract complex (not simplex) element class.
/// ----------------------------------------------------------------- ///
class ComplexElement : public Element {
public:
  real_t Diam() const final;
  real_t Volume() const final;
  vec3_t Normal() const final;
  vec3_t CenterPos() const final;

  /// @brief Make splitting into the simplex parts.
  virtual ShapeDescArray MakeSimplicesDesc() const = 0;

private:
  template<class Func>
  void ForEachSimplex_(Func&& func) const;

}; // class ComplexElement

/// ----------------------------------------------------------------- ///
/// @brief Dummy nodal element.
/// ----------------------------------------------------------------- ///
class Node final : public ElementHelper_<ShapeType::Node, 1, SimplexElement> {
public:
  real_t Volume() const final;
  vec3_t Normal() const final;
  ShapeDescArray make_edges_desc() const final;
  ShapeDescArray make_faces_desc() const final;
}; // class Node

/// ----------------------------------------------------------------- ///
/// @brief Segmental element.
/// @verbatim
///
///  n0 O f0
///      \
///       \         e0 = (n0,n1)
///        v e0     f0 = (n0)
///         \       f1 = (n1)
///          \
///        n1 O f1
///
/// @endverbatim
/// ----------------------------------------------------------------- ///
class Segment final :
    public ElementHelper_<ShapeType::Segment, 2, SimplexElement> {
public:
  real_t Volume() const final;
  vec3_t Normal() const final;
  vec3_t Dir() const final;
  ShapeDescArray make_edges_desc() const final;
  ShapeDescArray make_faces_desc() const final;
}; // class tSegmentShape

/// ----------------------------------------------------------------- ///
/// Triangular element.
/// @verbatim
///           n2
///           O           e0 = f0 = (n0,n1)
///          / \          e1 = f1 = (n1,n2)
///         /   \         e2 = f2 = (n2,n0)
///  e2/f2 v     ^ e1/f1
///       /       \
///      /         \
///  n0 O----->-----O n1
///        e0/f0
/// @endverbatim
/// ----------------------------------------------------------------- ///
class Triangle final :
    public ElementHelper_<ShapeType::Triangle, 3, SimplexElement> {
public:
  real_t Volume() const final;
  vec3_t Normal() const final;
  ShapeDescArray make_edges_desc() const final;
  ShapeDescArray make_faces_desc() const final;
}; // class Triangle

/// ----------------------------------------------------------------- ///
/// @brief Quadrangular element.
/// @verbatim
///               e2/f2
///       n3 O-----<-----O n2    e0 = f0 = (n0,n1)
///         /           /        e1 = f2 = (n1,n2)
///  e3/f3 v           ^ e1/f1   e2 = f2 = (n2,n3)
///       /           /          e3 = f3 = (n3,n0)
///   n0 O----->-----O n1     split = ((n0,n1,n2),(n2,n3,n0))
///          e0/f0
/// @endverbatim
/// ----------------------------------------------------------------- ///
class Quadrangle final :
    public ElementHelper_<ShapeType::Quadrangle, 4, ComplexElement> {
public:
  ShapeDescArray make_edges_desc() const final;
  ShapeDescArray make_faces_desc() const final;
  ShapeDescArray MakeSimplicesDesc() const final;
}; // class Quadrangle

/// ----------------------------------------------------------------- ///
/// @brief Tetrahedral element.
/// @verbatim
///                    f3
///               n3   ^
///                O   |
///         f1    /|\. |     f2         e0 = (n0,n1)
///         ^    / | `\.     ^          e1 = (n1,n2)
///          \  /  |   `\.  /           e2 = (n2,n0)
///           \`   |   | `\/            e3 = (n0,n3)
///           ,\   |   o  /`\           e4 = (n1,n3)
///       e3 ^  *  |     *   `^.e5      e5 = (n2,n3)
///         /   e4 ^           `\       f0 = (n0,n2,n1)
///     n0 O-------|--<----------O n2   f1 = (n0,n1,n3)
///         \      |  e2       ,/       f2 = (n1,n2,n3)
///          \     |     o   ,/`        f3 = (n2,n0,n3)
///           \    ^ e4  | ,/`
///         e0 v   |     ,^ e1
///             \  |   ,/`
///              \ | ,/` |
///               \|/`   |
///                O n1  v
///                      f0
/// @endverbatim
/// ----------------------------------------------------------------- ///
class Tetrahedron final :
    public ElementHelper_<ShapeType::Tetrahedron, 4, SimplexElement> {
public:
  real_t Volume() const final;
  ShapeDescArray make_edges_desc() const final;
  ShapeDescArray make_faces_desc() const final;
}; // class Tetrahedron

/// ----------------------------------------------------------------- ///
/// @brief Pyramidal element.
/// @verbatim
///                                n4                      e0 = (n0,n1)
///                  f3           ,O                       e1 = (n1,n2)
///                   ^        ,/`/|\     f1               e2 = (n2,n3)
///                    \    ,/`  / | \    ^                e3 = (n3,n0)
///                     \,/`    /  |  \  /                 e4 = (n0,n4)
///                e7 ,/`\     /   |   \/                  e5 = (n1,n4)
///                ,^`    o   /    |   /\                  e6 = (n2,n4)
///             ,/`          /     |  *  \                 e7 = (n3,n4)
///  f4 <------------*      /   e6 ^   o--\---------> f2   f0 = (n0,n3,n2,n1)
///       ,/`              /       |       \               f1 = (n0,n1,n4)
///   n3 O-----<----------/--------O  n2    ^ e5           f2 = (n1,n2,n4)
///       `\.  e2        /          `\.      \             f3 = (n2,n3,n4)
///          `>.        ^ e4           `\. e1 \            f4 = (n3,n0,n4)
///          e3 `\.    /       o          `<.  \        split = ((n0,n1,n2,n4),
///                `\./        |             `\.\                (n2,n3,n0,n4))
///               n0 O-------------------->-----O n1
///                            |          e0
///                            |
///                            v
///                            f0
/// @endverbatim
/// ----------------------------------------------------------------- ///
class Pyramid final :
    public ElementHelper_<ShapeType::Pyramid, 5, ComplexElement> {
public:
  ShapeDescArray make_edges_desc() const final;
  ShapeDescArray make_faces_desc() const final;
  ShapeDescArray MakeSimplicesDesc() const final;
}; // class Pyramid

/// ----------------------------------------------------------------- ///
/// @brief Pentahedral element (triangular prism).
/// @verbatim
///                 f4
///                 ^  f2
///                 |  ^
///             e8  |  |
///      n3 O---<---|-------------O n5        e0 = (n0,n1)
///         |\      *  |        ,/|           e1 = (n1,n2)
///         | \        o      ,/` |           e2 = (n2,n0)
///         |  \         e7 ,^`   |           e3 = (n0,n3)
///         |   v e6      ,/`     |           e4 = (n1,n4)
///      e3 ^    \      ,/`       ^ e5        e5 = (n2,n5)
///         |     \   ,/`         |           e6 = (n3,n4)
///         |      \ /`        *-------> f1   e7 = (n4,n5)
///  f0 <-------*   @ n4          |           e8 = (n5,n3)
///         |       |             |           f0 = (n0,n1,n4,n3)
///      n0 O-------|---------<---O n2        f1 = (n1,n2,n5,n4)
///          \      |        e2 ,/            f2 = (n2,n0,n3,n5)
///           \     |     o   ,/`             f3 = (n0,n2,n1)
///            \    ^ e4  | ,/`               f4 = (n3,n4,n5)
///          e0 v   |     ,^ e1            split = ((n0,n1,n2,n4),
///              \  |   ,/|                         (n2,n0,n3,n4),
///               \ | ,/` |                         (n3,n5,n2,n4))
///                \|/`   |
///                 O n1  v
///                       f3
/// @endverbatim
/// ----------------------------------------------------------------- ///
class Pentahedron final :
    public ElementHelper_<ShapeType::Pentahedron, 6, ComplexElement> {
public:
  ShapeDescArray make_edges_desc() const final;
  ShapeDescArray make_faces_desc() const final;
  ShapeDescArray MakeSimplicesDesc() const final;
}; // class Pyramid

/// ----------------------------------------------------------------- ///
/// @brief Hexahedral element.
/// @verbatim
///                      f5
///                      ^       f2
///                      |       ^
///                   e9 |      /
///            n6 O---<--|----------O n5         e0 = (n0,n1)
///              /|      |    /    /|            e1 = (n1,n2)
///             / |      |   o    / |            e2 = (n2,n3)
///        e10 v  |      *    e8 ^  ^ e5         e3 = (n3,n0)
///           /   ^ e6          /   |            e4 = (n0,n4)
///          /    |      e11   /  *-------> f1   e5 = (n1,n5)
///      n7 O------------->---O n4  |            e6 = (n2,n6)
///  f3 <---|--o  |           |     |            e7 = (n3,n7)
///         |  n2 O---<-------|-----O n1         e8 = (n4,n5)
///         |    /    e1      |    /             e9 = (n5,n6)
///      e7 ^   /          e4 ^   /             e10 = (n6,n7)
///         |  v e2  *        |  ^ e0           e11 = (n7,n4)
///         | /     /    o    | /                f0 = (n0,n3,n2,n1)
///         |/     /     |    |/                 f1 = (n0,n1,n5,n4)
///      n3 O-----/-->--------O n0               f2 = (n1,n2,n6,n5)
///              /   e3  |                       f3 = (n2,n3,n7,n6)
///             /        v                       f4 = (n0,n4,n7,n3)
///            v         f0                      f5 = (n4,n5,n6,n7)
///            f4
/// @endverbatim
/// ----------------------------------------------------------------- ///
class Hexahedron final :
    public ElementHelper_<ShapeType::Hexahedron, 8, ComplexElement> {
public:
  ShapeDescArray make_edges_desc() const final;
  ShapeDescArray make_faces_desc() const final;
  ShapeDescArray MakeSimplicesDesc() const final;
}; // class Hexahedron

} // namespace Storm
