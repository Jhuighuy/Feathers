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

/// @brief shape_type type.
enum class ShapeType : byte_t {
  Node,
  Segment,
  Triangle,
  Quadrangle,
  Tetrahedron,
  Pyramid,
  Pentahedron,
  Hexahedron,
}; // enum class ShapeType

/// @brief shape_type description.
struct ShapeDesc {
  ShapeType shape_type;
  std::vector<NodeIndex> node_indices;
}; // struct ShapeDesc

/// ----------------------------------------------------------------- ///
/// @brief Abstract shape.
/// ----------------------------------------------------------------- ///
class Shape : public NonCopyable {
protected:

  std::vector<NodeIndex> node_indices_;

  template<class... Indices>
  /*constexpr*/ auto part_desc_(ShapeType shape_type,
                                Indices... node_locals) const {
    return ShapeDesc{shape_type, {node_indices_[node_locals]...}};
  }

public:

  /// ---------------------------------------------------------------- ///
  /// @name Construction.
  /// ---------------------------------------------------------------- ///
  /// @{

  /// @brief Construct a shape.
  Shape() = default;

  /// @brief Destruct a shape.
  virtual ~Shape() = default;

  /// @brief Construct a new shape
  ///   with a description @p desc and a node position array @p node_coords.
  static std::unique_ptr<Shape> make(const ShapeDesc& shape_desc);

  /// @} // Construction.

  /// ---------------------------------------------------------------- ///
  /// @name Geometry.
  /// ---------------------------------------------------------------- ///
  /// @{

  /// @brief Compute the shape diameter
  virtual real_t
  diam([[maybe_unused]] const NodeCoordsVector& node_coords) const {
    return qnan;
  }

  /// @brief Compute the shape volume (area or length).
  virtual real_t
  volume([[maybe_unused]] const NodeCoordsVector& node_coords) const {
    return qnan;
  }

  /// @brief Compute the shape to element.
  virtual vec3_t
  normal([[maybe_unused]] const NodeCoordsVector& node_coords) const {
    return vec3_t(qnan);
  }

  /// @brief Compute the element direction.
  virtual vec3_t
  dir([[maybe_unused]] const NodeCoordsVector& node_coords) const {
    return vec3_t(qnan);
  }

  /// @brief Compute the shape center position.
  virtual vec3_t
  center_pos([[maybe_unused]] const NodeCoordsVector& node_coords) const {
    return vec3_t(qnan);
  }

  /// @} // Geometry.

  /// ---------------------------------------------------------------- ///
  /// @name Topology.
  /// ---------------------------------------------------------------- ///
  /// @{

  /// @brief Get the shape type.
  /*constexpr*/ virtual ShapeType shape_type() const noexcept = 0;

  /** Get element node indices. */
  /*constexpr*/ const std::vector<NodeIndex>& node_indices() const {
    return node_indices_;
  }

  /// @brief Number of nodes in the shape.
  /*constexpr*/ virtual size_t num_nodes() const noexcept = 0;

  /// @brief Make shape edges description array.
  /*constexpr*/ virtual std::vector<ShapeDesc> make_edges_desc() const = 0;

  /// @brief Make shape faces description array.
  /*constexpr*/ virtual std::vector<ShapeDesc> make_faces_desc() const = 0;

  /// @} // Topology.

}; // class Shape

template<ShapeType ShapeType_, size_t NumNodes_, class Base>
class ShapeHelper_ : public Base {
public:

  /*constexpr*/ ShapeType shape_type() const noexcept final {
    return ShapeType_;
  }

  /*constexpr*/ size_t num_nodes() const noexcept final {
    return NumNodes_;
  }

}; // class ShapeHelper_<...>

/// ----------------------------------------------------------------- ///
/// @brief Abstract simplex element class.
/// ----------------------------------------------------------------- ///
class SimplexShape : public Shape {
public:

  real_t diam(const NodeCoordsVector& node_coords) const final {
    if (num_nodes() == 1) { return 0.0; }
    real_t diam{glm::length(node_coords[node_indices()[0]] -
                            node_coords[node_indices()[1]])};
    for (size_t i{2}; i < num_nodes(); ++i) {
      for (size_t j{0}; j < i; ++j) {
        diam = std::max(diam, glm::length(node_coords[node_indices()[i]] -
                                          node_coords[node_indices()[j]]));
      }
    }
    return diam;
  }

  vec3_t center_pos(const NodeCoordsVector& node_coords) const final {
    vec3_t sum_of_node_pos{node_coords[node_indices()[0]]};
    for (size_t i{1}; i < num_nodes(); ++i) {
      sum_of_node_pos += node_coords[node_indices()[i]];
    }
    return sum_of_node_pos / real_t(num_nodes());
  }

}; // SimplexElement

/// ----------------------------------------------------------------- ///
/// @brief Abstract complex (not simplex) shape.
/// ----------------------------------------------------------------- ///
class ComplexShape : public Shape {
public:

  real_t diam(const NodeCoordsVector& node_coords) const final {
    real_t diam{0.0};
    for_each_simplex_([&](const Shape& shape) {
      diam = std::max(diam, shape.diam(node_coords));
    });
    return diam;
  }

  real_t volume(const NodeCoordsVector& node_coords) const final {
    real_t volume{0.0};
    for_each_simplex_([&](const Shape& shape) { //
      volume += shape.volume(node_coords);
    });
    return volume;
  }

  vec3_t normal(const NodeCoordsVector& node_coords) const final {
    vec3_t weighted_sum_of_normals(0.0);
    for_each_simplex_([&](const Shape& shape) {
      weighted_sum_of_normals +=
          shape.volume(node_coords) * shape.normal(node_coords);
    });
    return glm::normalize(weighted_sum_of_normals);
  }

  vec3_t center_pos(const NodeCoordsVector& node_coords) const final {
    vec3_t weighted_sum_of_center_pos(0.0);
    real_t volume{0.0};
    for_each_simplex_([&](const Shape& shape) {
      real_t const part_volume{shape.volume(node_coords)};
      weighted_sum_of_center_pos += part_volume * shape.center_pos(node_coords);
      volume += part_volume;
    });
    return weighted_sum_of_center_pos / volume;
  }

  /// @brief Make splitting into the simplex parts.
  /*constexpr*/ virtual std::vector<ShapeDesc> make_simplices_desc() const = 0;

private:

  void for_each_simplex_(auto&& func) const {
    std::vector<ShapeDesc> simplices_desc{make_simplices_desc()};
    for (const ShapeDesc& simplex_desc : simplices_desc) {
      const auto simplex = Shape::make(simplex_desc);
      func(*simplex);
    }
  }

}; // class ComplexElement

/// ----------------------------------------------------------------- ///
/// @brief Nodal shape.
/// ----------------------------------------------------------------- ///
class Node final : public ShapeHelper_<ShapeType::Node, 1, SimplexShape> {
public:

  /*constexpr*/ real_t
  volume([[maybe_unused]] const NodeCoordsVector& node_coords) const final {
    return 1.0;
  }

  /*constexpr*/ vec3_t
  normal([[maybe_unused]] const NodeCoordsVector& node_coords) const final {
    /*constexpr*/ vec3_t right(1.0, 0.0, 0.0);
    return right;
  }

  /*constexpr*/ std::vector<ShapeDesc> make_edges_desc() const final {
    return {};
  }

  /*constexpr*/ std::vector<ShapeDesc> make_faces_desc() const final {
    return {};
  }

}; // class Node

/// ----------------------------------------------------------------- ///
/// @brief Segmental shape.
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
class Segment final : public ShapeHelper_<ShapeType::Segment, 2, SimplexShape> {
public:

  real_t volume(const NodeCoordsVector& node_coords) const final {
    return diam(node_coords);
  }

  vec3_t normal(const NodeCoordsVector& node_coords) const final {
    const vec3_t delta{node_coords[node_indices()[1]] -
                       node_coords[node_indices()[0]]};
    static /*constexpr*/ vec3_t up(0.0, 0.0, 1.0);
    return glm::normalize(glm::cross(delta, up));
  }

  vec3_t dir(const NodeCoordsVector& node_coords) const final {
    const vec3_t delta{node_coords[node_indices()[1]] -
                       node_coords[node_indices()[0]]};
    return glm::normalize(delta);
  }

  /*constexpr*/ std::vector<ShapeDesc> make_edges_desc() const final {
    return {part_desc_(ShapeType::Segment, 0, 1)};
  }

  /*constexpr*/ std::vector<ShapeDesc> make_faces_desc() const final {
    return {part_desc_(ShapeType::Node, 0), part_desc_(ShapeType::Node, 1)};
  }

}; // class Segment

/// ----------------------------------------------------------------- ///
/// Triangular shape.
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
    public ShapeHelper_<ShapeType::Triangle, 3, SimplexShape> {
public:

  real_t volume(const NodeCoordsVector& node_coords) const final {
    const vec3_t delta1{node_coords[node_indices()[1]] -
                        node_coords[node_indices()[0]]};
    const vec3_t delta2{node_coords[node_indices()[2]] -
                        node_coords[node_indices()[0]]};
    return 0.5 * glm::length(glm::cross(delta1, delta2));
  }

  vec3_t normal(const NodeCoordsVector& node_coords) const final {
    const vec3_t delta1{node_coords[node_indices()[1]] -
                        node_coords[node_indices()[0]]};
    const vec3_t delta2{node_coords[node_indices()[2]] -
                        node_coords[node_indices()[0]]};
    return glm::normalize(glm::cross(delta1, delta2));
  }

  /*constexpr*/ std::vector<ShapeDesc> make_edges_desc() const final {
    return {part_desc_(ShapeType::Segment, 0, 1),
            part_desc_(ShapeType::Segment, 1, 2),
            part_desc_(ShapeType::Segment, 2, 0)};
  }
  /*constexpr*/ std::vector<ShapeDesc> make_faces_desc() const final {
    return make_edges_desc();
  }

}; // class Triangle

/// ----------------------------------------------------------------- ///
/// @brief Quadrangular shape.
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
    public ShapeHelper_<ShapeType::Quadrangle, 4, ComplexShape> {
public:

  /*constexpr*/ std::vector<ShapeDesc> make_edges_desc() const final {
    return {part_desc_(ShapeType::Segment, 0, 1),
            part_desc_(ShapeType::Segment, 1, 2),
            part_desc_(ShapeType::Segment, 2, 3),
            part_desc_(ShapeType::Segment, 3, 0)};
  }

  /*constexpr*/ std::vector<ShapeDesc> make_faces_desc() const final {
    return make_edges_desc();
  }

  /*constexpr*/ std::vector<ShapeDesc> make_simplices_desc() const final {
    return {part_desc_(ShapeType::Triangle, 0, 1, 2),
            part_desc_(ShapeType::Triangle, 2, 3, 0)};
    // return {
    //   part_desc_(ShapeType::Triangle, 0, 1, 3),
    //   part_desc_(ShapeType::Triangle, 1, 2, 3)};
  }

}; // class Quadrangle

/// ----------------------------------------------------------------- ///
/// @brief Tetrahedral shape.
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
    public ShapeHelper_<ShapeType::Tetrahedron, 4, SimplexShape> {
public:

  real_t volume(const NodeCoordsVector& node_coords) const final {
    const vec3_t delta1{node_coords[node_indices()[1]] -
                        node_coords[node_indices()[0]]};
    const vec3_t delta2{node_coords[node_indices()[2]] -
                        node_coords[node_indices()[0]]};
    const vec3_t delta3{node_coords[node_indices()[3]] -
                        node_coords[node_indices()[0]]};
    return std::abs(glm::dot(delta1, glm::cross(delta2, delta3))) / 6.0;
  }

  /*constexpr*/ std::vector<ShapeDesc> make_edges_desc() const final {
    return {part_desc_(ShapeType::Segment, 0, 1),
            part_desc_(ShapeType::Segment, 1, 2),
            part_desc_(ShapeType::Segment, 2, 0),
            part_desc_(ShapeType::Segment, 0, 3),
            part_desc_(ShapeType::Segment, 1, 3),
            part_desc_(ShapeType::Segment, 2, 3)};
  }

  /*constexpr*/ std::vector<ShapeDesc> make_faces_desc() const final {
    return {part_desc_(ShapeType::Triangle, 0, 2, 1),
            part_desc_(ShapeType::Triangle, 0, 1, 3),
            part_desc_(ShapeType::Triangle, 1, 2, 3),
            part_desc_(ShapeType::Triangle, 2, 0, 3)};
  }

}; // class Tetrahedron

/// ----------------------------------------------------------------- ///
/// @brief Pyramidal shape.
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
class Pyramid final : public ShapeHelper_<ShapeType::Pyramid, 5, ComplexShape> {
public:

  /*constexpr*/ std::vector<ShapeDesc> make_edges_desc() const final {
    return {part_desc_(ShapeType::Segment, 0, 1),
            part_desc_(ShapeType::Segment, 1, 2),
            part_desc_(ShapeType::Segment, 2, 3),
            part_desc_(ShapeType::Segment, 3, 0),
            part_desc_(ShapeType::Segment, 0, 4),
            part_desc_(ShapeType::Segment, 1, 4),
            part_desc_(ShapeType::Segment, 2, 4),
            part_desc_(ShapeType::Segment, 3, 4)};
  }

  /*constexpr*/ std::vector<ShapeDesc> make_faces_desc() const final {
    return {part_desc_(ShapeType::Quadrangle, 0, 3, 2, 1),
            part_desc_(ShapeType::Triangle, 0, 1, 4),
            part_desc_(ShapeType::Triangle, 1, 2, 4),
            part_desc_(ShapeType::Triangle, 2, 3, 4),
            part_desc_(ShapeType::Triangle, 3, 0, 4)};
  }

  /*constexpr*/ std::vector<ShapeDesc> make_simplices_desc() const final {
    return {part_desc_(ShapeType::Tetrahedron, 0, 1, 2, 4),
            part_desc_(ShapeType::Tetrahedron, 2, 3, 0, 4)};
    // return {
    //   part_desc_(ShapeType::Tetrahedron, 0, 1, 3, 4),
    //   part_desc_(ShapeType::Tetrahedron, 1, 2, 3, 4)};
  }

}; // class Pyramid

/// ----------------------------------------------------------------- ///
/// @brief Pentahedral shape (triangular prism).
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
    public ShapeHelper_<ShapeType::Pentahedron, 6, ComplexShape> {
public:

  /*constexpr*/ std::vector<ShapeDesc> make_edges_desc() const final {
    return {part_desc_(ShapeType::Segment, 0, 1),
            part_desc_(ShapeType::Segment, 1, 2),
            part_desc_(ShapeType::Segment, 2, 0),
            part_desc_(ShapeType::Segment, 0, 3),
            part_desc_(ShapeType::Segment, 1, 4),
            part_desc_(ShapeType::Segment, 2, 5),
            part_desc_(ShapeType::Segment, 3, 4),
            part_desc_(ShapeType::Segment, 4, 5),
            part_desc_(ShapeType::Segment, 5, 3)};
  }

  /*constexpr*/ std::vector<ShapeDesc> make_faces_desc() const final {
    return {part_desc_(ShapeType::Quadrangle, 0, 1, 4, 3),
            part_desc_(ShapeType::Quadrangle, 1, 2, 5, 4),
            part_desc_(ShapeType::Quadrangle, 2, 0, 3, 5),
            part_desc_(ShapeType::Triangle, 0, 2, 1),
            part_desc_(ShapeType::Triangle, 3, 4, 5)};
  }

  /*constexpr*/ std::vector<ShapeDesc> make_simplices_desc() const final {
    return {part_desc_(ShapeType::Tetrahedron, 0, 1, 2, 4),
            part_desc_(ShapeType::Tetrahedron, 2, 0, 3, 4),
            part_desc_(ShapeType::Tetrahedron, 3, 5, 2, 4)};
  }

}; // class Pyramid

/// ----------------------------------------------------------------- ///
/// @brief Hexahedral shape.
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
    public ShapeHelper_<ShapeType::Hexahedron, 8, ComplexShape> {
public:

  /*constexpr*/ std::vector<ShapeDesc> make_edges_desc() const final {
    return {part_desc_(ShapeType::Segment, 0, 1),
            part_desc_(ShapeType::Segment, 1, 2),
            part_desc_(ShapeType::Segment, 2, 3),
            part_desc_(ShapeType::Segment, 3, 0),
            part_desc_(ShapeType::Segment, 0, 4),
            part_desc_(ShapeType::Segment, 1, 5),
            part_desc_(ShapeType::Segment, 2, 6),
            part_desc_(ShapeType::Segment, 3, 7),
            part_desc_(ShapeType::Segment, 4, 5),
            part_desc_(ShapeType::Segment, 5, 6),
            part_desc_(ShapeType::Segment, 6, 7),
            part_desc_(ShapeType::Segment, 7, 4)};
  }

  /*constexpr*/ std::vector<ShapeDesc> make_faces_desc() const final {
    return {part_desc_(ShapeType::Quadrangle, 0, 3, 2, 1),
            part_desc_(ShapeType::Quadrangle, 0, 1, 5, 4),
            part_desc_(ShapeType::Quadrangle, 1, 2, 6, 5),
            part_desc_(ShapeType::Quadrangle, 2, 3, 7, 6),
            part_desc_(ShapeType::Quadrangle, 0, 4, 7, 3),
            part_desc_(ShapeType::Quadrangle, 4, 5, 6, 7)};
  }

  /*constexpr*/ std::vector<ShapeDesc> make_simplices_desc() const final {
    return {part_desc_(ShapeType::Tetrahedron, 0, 3, 1, 4),
            part_desc_(ShapeType::Tetrahedron, 3, 2, 1, 6),
            part_desc_(ShapeType::Tetrahedron, 4, 5, 6, 1),
            part_desc_(ShapeType::Tetrahedron, 4, 6, 7, 3),
            part_desc_(ShapeType::Tetrahedron, 4, 3, 1, 6)};
  }

}; // class Hexahedron

} // namespace Storm
