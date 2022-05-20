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

#include "SkunkBase.hh"

#include <span>

namespace feathers {

/// @brief Element Shape type.
enum class ShapeType : byte_t {
  Null,
  Node,
  Segment2,
  Triangle3,
  Quadrangle4,
  Tetrahedron4,
  Pyramid5,
  Pentahedron6,
  Hexahedron8,
}; // enum class ShapeType

/// @brief Element description.
struct ElementDesc {
  ShapeType Shape;
  std::vector<size_t> NodeIndices;
}; // ElementDesc

/// @brief Array of the element descriptions.
using tElementDescList = std::vector<ElementDesc>;

/// ----------------------------------------------------------------- ///
/// @brief Abstract element class.
/// ----------------------------------------------------------------- ///
class iElement {
protected:
  size_t NumGlobalNodes_ = 0;
  const vec3_t* GlobalNodeCoords_ = nullptr;
  std::vector<size_t> NodeIndices_;

public:

  /// @brief Construct a new element object \
  ///   with a @p description and a @p nodes_span.
  static std::unique_ptr<iElement> make(ElementDesc&& desc,
                                        size_t num_global_nodes,
                                        const vec3_t* global_node_coords);

  virtual ~iElement() = default;

  /** Get element node indices. */
  std::vector<size_t> const& get_nodes() const {
    return NodeIndices_;
  }

  /** Get node position. */
  const vec3_t& get_node_coords(size_t node_local) const {
      StormAssert(node_local < NodeIndices_.size());
    return GlobalNodeCoords_[NodeIndices_[node_local]];
  }

  template<typename... tIndex>
  ElementDesc get_part(ShapeType part_shape, tIndex... node_locals) const {
    return { part_shape, std::vector<size_t>{ NodeIndices_[node_locals]... } };
  }

  // ---------------------------------------------------------------------- //
  // ---------------------------------------------------------------------- //

  /** Get element diameter. */
  virtual real_t get_diameter() const {
    return qnan;
  }
  /** Get element length/area/volume. */
  virtual real_t LenAreaOrVolume() const {
    return qnan;
  }
  /** Get normal to element. */
  virtual vec3_t Normal() const {
    return vec3_t(qnan);
  }
  /** Get element direction. */
  virtual vec3_t Dir() const {
    return vec3_t(qnan);
  }
  /** Get element barycenter. */
  virtual vec3_t CenterPos() const {
    return vec3_t(qnan);
  }

  // ---------------------------------------------------------------------- //
  // ---------------------------------------------------------------------- //

  /** Get element Shape. */
  virtual ShapeType Shape() const = 0;

  /** Number of nodes in the element. */
  virtual size_t num_nodes() const = 0;

  /** Number of edges in the element. */
  size_t num_edges() const {
    return get_edges_desc().size();
  }
  /** Get element edges description. */
  virtual tElementDescList get_edges_desc() const = 0;

  /** Number of faces in the element. */
  size_t num_faces() const {
    return get_faces_desc().size();
  }
  /** Get element faces description. */
  virtual tElementDescList get_faces_desc() const = 0;

};  // class iElement

/**
 * Abstract simplex element class.
 */
class iSimplexElement : public iElement {
};  // class iSimplexElement

/**
 * Abstract complex (not simplex) element class.
 */
class iComplexElement : public iElement {
public:
  real_t get_diameter() const final;
  real_t LenAreaOrVolume() const final;
  vec3_t Normal() const final;
  vec3_t CenterPos() const final;

  /**
   * Get splitting into the simplex parts.
   */
  virtual tElementDescList get_simplicial_parts(size_t partition_index) const = 0;

private:
  template<typename tFunc>
  void for_each_simplex_(tFunc func) const;
};  // class iComplexElement


/**
 * Dummy nodal element.
 */
class cNode final : public iSimplexElement {
public:
  real_t get_diameter() const final;
  real_t LenAreaOrVolume() const final;
  vec3_t Normal() const final;
  vec3_t CenterPos() const final;

  ShapeType Shape() const final;
  size_t num_nodes() const final;
  tElementDescList get_edges_desc() const final;
  tElementDescList get_faces_desc() const final;
};  // class cNode

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
class Segment final : public iSimplexElement {
public:
  real_t get_diameter() const final;
  real_t LenAreaOrVolume() const final;
  vec3_t Normal() const final;
  vec3_t Dir() const final;
  vec3_t CenterPos() const final;

  ShapeType Shape() const final;
  size_t num_nodes() const final;
  tElementDescList get_edges_desc() const final;
  tElementDescList get_faces_desc() const final;
};  // class tSegmentShape

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
class Triangle final : public iSimplexElement {
public:
  real_t get_diameter() const final;
  real_t LenAreaOrVolume() const final;
  vec3_t Normal() const final;
  vec3_t CenterPos() const final;

  ShapeType Shape() const final;
  size_t num_nodes() const final;
  tElementDescList get_edges_desc() const final;
  tElementDescList get_faces_desc() const final;
};  // class Triangle

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
class Quadrangle final : public iComplexElement {
public:
  ShapeType Shape() const final;
  size_t num_nodes() const final;
  tElementDescList get_edges_desc() const final;
  tElementDescList get_faces_desc() const final;
  tElementDescList get_simplicial_parts(size_t partition_index) const final;
};  // class Quadrangle

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
class Tetrahedron final : public iSimplexElement {
public:
  real_t get_diameter() const final;
  real_t LenAreaOrVolume() const final;
  vec3_t CenterPos() const final;

  ShapeType Shape() const final;
  size_t num_nodes() const final;
  tElementDescList get_edges_desc() const final;
  tElementDescList get_faces_desc() const final;
};  // class Tetrahedron

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
class Pyramid final : public iComplexElement {
public:
  ShapeType Shape() const final;
  size_t num_nodes() const final;
  tElementDescList get_edges_desc() const final;
  tElementDescList get_faces_desc() const final;
  tElementDescList get_simplicial_parts(size_t partition_index) const final;
};  // class Pyramid

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
///          e0 v   |     ,^ e1
///              \  |   ,/|
///               \ | ,/` |
///                \|/`   |
///                 O n1  v
///                       f3
/// @endverbatim
/// ----------------------------------------------------------------- ///
class Pentahedron final : public iComplexElement {
public:
  ShapeType Shape() const final;
  size_t num_nodes() const final;
  tElementDescList get_edges_desc() const final;
  tElementDescList get_faces_desc() const final;
  tElementDescList get_simplicial_parts(size_t partition_index) const final;
}; // class Pyramid

/// ----------------------------------------------------------------- ///
/// @brief Hexahedral element.
/// @verbatim
///                      f5
///                      ^   f2
///                      |   ^
///                   e9 |   |
///            n6 O---<--|----------O n5         e0 = (n0,n1)
///              /|      |   |     /|            e1 = (n1,n2)
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
///         | /      |   o    | /                f0 = (n0,n3,n2,n1)
///         |/       |   |    |/                 f1 = (n0,n1,n5,n4)
///      n3 O--->----|--------O n0               f2 = (n1,n2,n6,n5)
///             e3   |   |                       f3 = (n2,n3,n7,n6)
///                  |   v                       f4 = (n0,n4,n7,n3)
///                  v   f0                      f5 = (n4,n5,n6,n7)
///                  f4
/// @endverbatim
/// ----------------------------------------------------------------- ///
class Hexahedron final : public iComplexElement {
public:
  ShapeType Shape() const final;
  size_t num_nodes() const final;
  tElementDescList get_edges_desc() const final;
  tElementDescList get_faces_desc() const final;
  tElementDescList get_simplicial_parts(size_t partition_index) const final;
}; // class Hexahedron

} // namespace feathers
