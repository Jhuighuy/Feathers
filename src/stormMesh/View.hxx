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

#include <range/v3/all.hpp>

#include <stormMesh/Base.hxx>
#include <stormUtils/Parallel.hh>

namespace Storm {

template<class> class BaseNodeView;
template<class> class BaseEdgeView;
template<class> class BaseFaceView;
template<class> class BaseCellView;

template<class Mesh>
BaseNodeView(Mesh&, NodeIndex) -> BaseNodeView<Mesh>;
template<class Mesh>
BaseEdgeView(Mesh&, EdgeIndex) -> BaseEdgeView<Mesh>;
template<class Mesh>
BaseFaceView(Mesh&, FaceIndex) -> BaseFaceView<Mesh>;
template<class Mesh>
BaseCellView(Mesh&, CellIndex) -> BaseCellView<Mesh>;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base element view.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Mesh, class Tag>
class BaseElementView {
protected:
  Mesh* mesh_;
  Index<Tag> index_;

  template<class, class>
  friend class BaseElementView;

  // "NOLINT(...)" should not be here, this is due to a bug in clangd.
  BaseElementView( // NOLINT(cppcoreguidelines-pro-type-member-init)
    Mesh& mesh, Index<Tag> index) noexcept :
    mesh_{&mesh}, index_{index} {
    storm_assert(index_ != npos);
  }

  template<class OtherMesh>
  BaseElementView( // NOLINT(google-explicit-constructor,cppcoreguidelines-pro-type-member-init)
    BaseElementView<OtherMesh, Tag> const& other) noexcept :
      mesh_{other.mesh_}, index_{other.index_} {
    storm_assert(index_ != npos);
  }

public:

  /// @brief Cast to index operator.
  /// @{
  operator Index<Tag>() const noexcept {
    return index_;
  }
  FEATHERS_DEPRECATED operator size_t() const noexcept {
    return static_cast<size_t>(index_);
  }
  /// @}

  /// @brief Comparison operator.
  auto operator<=>(BaseElementView const& other) const noexcept {
    storm_assert(mesh_ == other.mesh_);
    return index_ <=> other.index_;
  }

  /// @brief Get mark. 
  Index<MarkTag<Tag>> mark() const noexcept {
    return mesh_->mark(index_);
  }

  /// @brief Get shape type.
  ShapeType shapeType() const noexcept {
    return mesh_->shapeType(index_);
  }

  /// @brief Get shape.
  std::unique_ptr<Element> shape() const {
    return mesh_->shape(index_);
  }

  /// @brief Ranges of the adjacent nodes.
  auto adjacent_nodes() const noexcept {
    return mesh_->adjacent_nodes(index_) |
      views::transform([&mesh = *mesh_](NodeIndex node_index) {
        return BaseNodeView(mesh, node_index);
      });
  }

  /// @brief Ranges of the adjacent edges.
  auto adjacent_edges() const noexcept {
    return mesh_->adjacent_edges(index_) |
      views::transform([&mesh = *mesh_](EdgeIndex edge_index) {
        return BaseEdgeView(mesh, edge_index);
      });
  }

  /// @brief Ranges of the adjacent faces.
  auto adjacent_faces() const noexcept {
    return mesh_->adjacent_faces(index_) |
      views::transform([&mesh = *mesh_](FaceIndex face_index) {
        return BaseFaceView(mesh, face_index);
      });
  }

  /// @brief Ranges of the adjacent cells.
  auto adjacent_cells() const noexcept {
    return mesh_->adjacent_cells(index_) |
      views::transform([&mesh = *mesh_](CellIndex cell_index) {
        return BaseCellView(mesh, cell_index);
      });
  }

  /// @brief Sequentially iterate through all the adjacent nodes.
  void for_each_node(auto&& func) const noexcept {
    ranges::for_each(adjacent_nodes(), func);
  }

  /// @brief Sequentially iterate through all the adjacent edges.
  void for_each_edge(auto&& func) const noexcept {
    ranges::for_each(adjacent_edges(), func);
  }

  /// @brief Sequentially iterate through all the adjacent faces.
  /// @{
  void for_each_face(auto&& func) const noexcept {
    ranges::for_each(adjacent_faces(), func);
  }
  void for_each_face_cells(auto&& func) const noexcept {
    ranges::for_each(adjacent_faces(), [&](BaseFaceView<Mesh> face) {
      func(face.inner_cell(), face.outer_cell());
    });
  }
  /// @}

  /// @brief Sequentially iterate through all the adjacent cells.
  void for_each_cell(auto&& func) const noexcept {
    ranges::for_each(adjacent_cells(), func);
  }

}; // class BaseElementView

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base node view.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Mesh>
class BaseNodeView final : public BaseElementView<Mesh, NodeTag> {
public:

  /// @brief Construct base node view.
  BaseNodeView(Mesh& mesh, NodeIndex index) noexcept :
    BaseElementView<Mesh, NodeTag>(mesh, index) {
  }

  /// @brief Copy constructor.
  BaseNodeView( // NOLINT(google-explicit-constructor)
      BaseNodeView<std::remove_const_t<Mesh>> const& other) noexcept :
    BaseElementView<Mesh, NodeTag>(other) {
  }

  /// @brief Get node position.
  vec3_t coords() const noexcept {
    return this->mesh_->node_coords(this->index_);
  }

  /// @brief Set node position @p coords.
  void set_coords(vec3_t const& coords) const noexcept requires (!std::is_const_v<Mesh>) {
    this->mesh_->set_node_coords(this->index_, coords);
  }

}; // class BaseNodeView<...>

/// @brief Mesh node view.
/// @{
using NodeView = BaseNodeView<Mesh const>;
using MutableNodeView = BaseNodeView<Mesh>;
/// @}

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base edge view.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Mesh>
class BaseEdgeView final : public BaseElementView<Mesh, EdgeTag> {
public:

  /// @brief Construct base edge view.
  BaseEdgeView(Mesh& mesh, EdgeIndex index) noexcept :
    BaseElementView<Mesh, EdgeTag>(mesh, index) {
  }

  /// @brief Copy constructor.
  BaseEdgeView( // NOLINT(google-explicit-constructor)
      BaseEdgeView<std::remove_const_t<Mesh>> const& other) noexcept :
    BaseElementView<Mesh, EdgeTag>(other) {
  }

  /// @brief Get edge length. 
  real_t len() const noexcept {
    return this->mesh_->edge_len(this->index_);
  }

  /// @brief Get edge direction. 
  vec3_t dir() const noexcept {
    return this->mesh_->edge_dir(this->index_);
  }

}; // class BaseEdgeView<...>

/// @brief Mesh edge view.
/// @{
using EdgeView = BaseEdgeView<Mesh const>;
using MutableEdgeView = BaseEdgeView<Mesh>;
/// @}

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base face view.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Mesh>
class BaseFaceView final : public BaseElementView<Mesh, FaceTag> {
public:

  /// @brief Construct base face view.
  BaseFaceView(Mesh& mesh, FaceIndex index) noexcept :
    BaseElementView<Mesh, FaceTag>(mesh, index) {
  }

  /// @brief Copy constructor.
  BaseFaceView( // NOLINT(google-explicit-constructor)
      BaseFaceView<std::remove_const_t<Mesh>> const& other) noexcept :
    BaseElementView<Mesh, FaceTag>(other) {
  }

  /// @brief Get connected inner cell. 
  auto inner_cell() const noexcept {
    storm_assert(this->adjacent_cells().size() == 2);
    return this->adjacent_cells()[face_inner_cell_];
  }

  /// @brief Get connected outer cell. 
  auto outer_cell() const noexcept {
    storm_assert(this->adjacent_cells().size() == 2);
    return this->adjacent_cells()[face_outer_cell_];
  }

  /// @brief Get face area/length. 
  real_t area() const noexcept {
    return this->mesh_->face_area(this->index_);
  }

  /// @brief Get face normal. 
  vec3_t normal() const noexcept {
    return this->mesh_->face_normal(this->index_);
  }

  /// @brief Get face center position.
  vec3_t barycenter() const noexcept {
    return this->mesh_->face_barycenter(this->index_);
  }

}; // class BaseFaceView<...>

/// @brief Mesh face view.
/// @{
using FaceView = BaseFaceView<Mesh const>;
using MutableFaceView = BaseFaceView<Mesh>;
/// @}

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base cell view.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Mesh>
class BaseCellView final : public BaseElementView<Mesh, CellTag> {
public:

  /// @brief Construct base cell view.
  BaseCellView(Mesh& mesh, CellIndex index) noexcept :
    BaseElementView<Mesh, CellTag>(mesh, index) {
  }

  /// @brief Copy constructor.
  BaseCellView( // NOLINT(google-explicit-constructor)
      BaseCellView<std::remove_const_t<Mesh>> const& other) noexcept :
    BaseElementView<Mesh, CellTag>(other) {
  }

  /// @brief Get cell volume/area/length.
  real_t volume() const noexcept {
    return this->mesh_->cell_volume(this->index_);
  }

  /// @brief Get cell center position.
  vec3_t barycenter() const noexcept {
    return this->mesh_->cell_barycenter(this->index_);
  }

}; // class BaseCellView<...>

/// @brief Mesh cell view.
/// @{
using CellView = BaseCellView<Mesh const>;
using MutableCellView = BaseCellView<Mesh>;
/// @}

/// @brief Range of the @p mesh nodes 
///   (or nodes with a @p node_mark, if present).
auto node_views(auto& mesh, std::same_as<NodeMark> auto... node_mark) noexcept {
  static_assert(sizeof...(node_mark) <= 1);
  return mesh.nodes(node_mark...) |
    views::transform([&mesh](NodeIndex node_index) {
      return BaseNodeView(mesh, node_index);
    });
}

/// @brief Range of the @p mesh edges 
///   (or edges with an @p edge_mark, if present).
auto edge_views(auto& mesh, std::same_as<EdgeMark> auto... edge_mark) noexcept {
  static_assert(sizeof...(edge_mark) <= 1);
  return mesh.edges(edge_mark...) |
    views::transform([&mesh](EdgeIndex edge_index) {
      return BaseEdgeView(mesh, edge_index);
    });
}

/// @brief Range of the @p mesh faces 
///   (or faces with a @p face_mark, if present).
auto face_views(auto& mesh, std::same_as<FaceMark> auto... face_mark) noexcept {
  static_assert(sizeof...(face_mark) <= 1);
  return mesh.faces(face_mark...) |
    views::transform([&mesh](FaceIndex face_index) {
      return BaseFaceView(mesh, face_index);
    });
}

/// @brief Range of the @p mesh nodes 
///   (or nodes with a @p cell_mark, if present).
auto cell_views(auto& mesh, std::same_as<CellMark> auto... cell_mark) noexcept {
  static_assert(sizeof...(cell_mark) <= 1);
  return mesh.cells(cell_mark...) |
    views::transform([&mesh](CellIndex cell_index) {
      return BaseCellView(mesh, cell_index);
    });
}

/// @brief Range of the interior @p mesh nodes.
auto intr_node_views(auto& mesh) noexcept {
  return node_views(mesh, NodeMark{0});
}

/// @brief Range of the interior @p mesh edges.
auto intr_edge_views(auto& mesh) noexcept {
  return edge_views(mesh, EdgeMark{0});
}

/// @brief Range of the interior @p mesh faces.
auto intr_face_views(auto& mesh) noexcept {
  return face_views(mesh, FaceMark{0});
}

/// @brief Range of the interior @p mesh cells.
auto intr_cell_views(auto& mesh) noexcept {
  return cell_views(mesh, CellMark{0});
}

/// @brief Range of the boundary @p mesh nodes.
auto bndr_node_views(auto& mesh) noexcept {
  return node_views(mesh) | views::drop(mesh.nodes(NodeMark{0}).size());
}

/// @brief Range of the boundary @p mesh edges.
auto bndr_edge_views(auto& mesh) noexcept {
  return edge_views(mesh) | views::drop(mesh.edges(EdgeMark{0}).size());
}

/// @brief Range of the boundary @p mesh faces.
auto bndr_face_views(auto& mesh) noexcept {
  return face_views(mesh) | views::drop(mesh.faces(FaceMark{0}).size());
}

auto bndr_face_cell_views(auto& mesh) noexcept {
  return face_views(mesh) |
         views::drop(mesh.faces(FaceMark{0}).size()) |
         views::transform([](BaseFaceView<Mesh> face) {
      return std::pair(face.inner_cell(), face.outer_cell());
    });
}

template<class Mesh>
void for_each_bndr_face_cells(Mesh& mesh, auto&& func) noexcept {
  for_each(bndr_face_views(mesh), [&](BaseFaceView<Mesh> face) {
    func(face.inner_cell(), face.outer_cell());
  });
}

/// @brief Range of the boundary @p mesh cells.
auto bndr_cell_views(auto& mesh) noexcept {
  return cell_views(mesh) | views::drop(mesh.cells(CellMark{0}).size());
}

} // namespace feathers
