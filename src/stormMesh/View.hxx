/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to
/// deal in the Software without restriction, including without limitation the
/// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
/// IN THE SOFTWARE.

#pragma once

#include <stormMesh/Base.hxx>
#include <stormUtils/Parallel.hh>

namespace Storm {

template<class>
class BaseNodeView;
template<class>
class BaseEdgeView;
template<class>
class BaseFaceView;
template<class>
class BaseCellView;

template<class Mesh>
BaseNodeView(Mesh&, NodeIndex) -> BaseNodeView<Mesh>;
template<class Mesh>
BaseEdgeView(Mesh&, EdgeIndex) -> BaseEdgeView<Mesh>;
template<class Mesh>
BaseFaceView(Mesh&, FaceIndex) -> BaseFaceView<Mesh>;
template<class Mesh>
BaseCellView(Mesh&, CellIndex) -> BaseCellView<Mesh>;

/// @brief Base element view.
// clang-format off
template<class Mesh, class Tag>
  requires std::is_object_v<Mesh>
class BaseElementView {
  // clang-format on
protected:

  Mesh* mesh_;
  Index<Tag> index_;

  template<class, class>
  friend class BaseElementView;

  /// @brief Initialize an element view with a @p mesh and @p index.
  constexpr BaseElementView(Mesh& mesh, Index<Tag> index) noexcept
      : mesh_{&mesh}, index_{index} {
    STORM_ASSERT_(index_ != npos);
  }

  /// @brief Copy construct an element view.
  template<class OtherMesh>
  constexpr BaseElementView( // NOLINT(google-explicit-constructor)
      const BaseElementView<OtherMesh, Tag>& other) noexcept
      : mesh_{other.mesh_}, index_{other.index_} {
    STORM_ASSERT_(index_ != npos);
  }

public:

  /// @brief Cast to index operator.
  /// @{
  [[nodiscard]] constexpr operator Index<Tag>() const noexcept {
    return index_;
  }
  FEATHERS_DEPRECATED constexpr operator size_t() const noexcept {
    return static_cast<size_t>(index_);
  }
  /// @}

  /// @brief Comparison operator.
  [[nodiscard]] constexpr auto
  operator<=>(const BaseElementView& other) const noexcept {
    STORM_ASSERT_(
        mesh_ == other.mesh_,
        "Can not compare element corresponding to the different meshes!");
    return index_ <=> other.index_;
  }

  /// @brief Get mark.
  [[nodiscard]] constexpr auto mark() const noexcept {
    return mesh_->mark(index_);
  }

  /// @brief Get shape type.
  [[nodiscard]] constexpr auto shape_type() const noexcept {
    return mesh_->shape_type(index_);
  }

  /// @brief Get shape.
  [[nodiscard]] constexpr auto shape() const {
    return mesh_->shape(index_);
  }

  /// @brief Ranges of the adjacent nodes.
  [[nodiscard]] constexpr auto adjacent_nodes() const noexcept {
    return mesh_->adjacent_nodes(index_) |
           views::transform([&mesh = *mesh_](NodeIndex node_index) {
             return BaseNodeView(mesh, node_index);
           });
  }

  /// @brief Ranges of the adjacent edges.
  [[nodiscard]] constexpr auto adjacent_edges() const noexcept {
    return mesh_->adjacent_edges(index_) |
           views::transform([&mesh = *mesh_](EdgeIndex edge_index) {
             return BaseEdgeView(mesh, edge_index);
           });
  }

  /// @brief Ranges of the adjacent faces.
  [[nodiscard]] constexpr auto adjacent_faces() const noexcept {
    return mesh_->adjacent_faces(index_) |
           views::transform([&mesh = *mesh_](FaceIndex face_index) {
             return BaseFaceView(mesh, face_index);
           });
  }

  /// @brief Ranges of the adjacent cells.
  [[nodiscard]] constexpr auto adjacent_cells() const noexcept {
    return mesh_->adjacent_cells(index_) |
           views::transform([&mesh = *mesh_](CellIndex cell_index) {
             return BaseCellView(mesh, cell_index);
           });
  }

  /// @brief Sequentially iterate through all the adjacent nodes.
  constexpr void for_each_node(auto func) const noexcept {
    ranges::for_each(adjacent_nodes(), func);
  }

  /// @brief Sequentially iterate through all the adjacent edges.
  constexpr void for_each_edge(auto func) const noexcept {
    ranges::for_each(adjacent_edges(), func);
  }

  /// @brief Sequentially iterate through all the adjacent faces.
  /// @{
  constexpr void for_each_face(auto func) const noexcept {
    ranges::for_each(adjacent_faces(), func);
  }
  constexpr void for_each_face_cells(auto func) const noexcept {
    ranges::for_each(adjacent_faces(), [&](BaseFaceView<Mesh> face) {
      func(face.inner_cell(), face.outer_cell());
    });
  }
  /// @}

  /// @brief Sequentially iterate through all the adjacent cells.
  constexpr void for_each_cell(auto func) const noexcept {
    ranges::for_each(adjacent_cells(), func);
  }

}; // class BaseElementView

/// @brief Base node view.
template<class Mesh>
class BaseNodeView final : public BaseElementView<Mesh, NodeTag> {
public:

  /// @brief Construct base node view.
  constexpr BaseNodeView(Mesh& mesh, NodeIndex index) noexcept
      : BaseElementView<Mesh, NodeTag>(mesh, index) {}

  /// @brief Copy construct a node view.
  constexpr BaseNodeView( // NOLINT(google-explicit-constructor)
      const BaseNodeView<std::remove_const_t<Mesh>>& other) noexcept
      : BaseElementView<Mesh, NodeTag>(other) {}

  /// @brief Get the node coordinates.
  [[nodiscard]] constexpr auto coords() const noexcept {
    return this->mesh_->node_coords(this->index_);
  }

  /// @brief Set the node coordinates @p coords.
  // clang-format off
  constexpr void set_coords(const auto& coords) const noexcept
      requires(!std::is_const_v<Mesh>) {
    // clang-format on
    this->mesh_->set_node_coords(this->index_, coords);
  }

}; // class BaseNodeView

/// @brief Mesh node view.
/// @{
using NodeView = BaseNodeView<const Mesh>;
using MutableNodeView = BaseNodeView<Mesh>;
/// @}

/// @brief Base edge view.
template<class Mesh>
class BaseEdgeView final : public BaseElementView<Mesh, EdgeTag> {
public:

  /// @brief Construct base edge view.
  constexpr BaseEdgeView(Mesh& mesh, EdgeIndex index) noexcept
      : BaseElementView<Mesh, EdgeTag>(mesh, index) {}

  /// @brief Copy construct an edge view.
  constexpr BaseEdgeView( // NOLINT(google-explicit-constructor)
      const BaseEdgeView<std::remove_const_t<Mesh>>& other) noexcept
      : BaseElementView<Mesh, EdgeTag>(other) {}

  /// @brief Get the edge length.
  [[nodiscard]] constexpr auto len() const noexcept {
    return this->mesh_->edge_len(this->index_);
  }

  /// @brief Get the edge direction.
  [[nodiscard]] constexpr auto dir() const noexcept {
    return this->mesh_->edge_dir(this->index_);
  }

}; // class BaseEdgeView

/// @brief Mesh edge view.
/// @{
using EdgeView = BaseEdgeView<const Mesh>;
using MutableEdgeView = BaseEdgeView<Mesh>;
/// @}

/// @brief Base face view.
template<class Mesh>
class BaseFaceView final : public BaseElementView<Mesh, FaceTag> {
public:

  /// @brief Construct base face view.
  constexpr BaseFaceView(Mesh& mesh, FaceIndex index) noexcept
      : BaseElementView<Mesh, FaceTag>(mesh, index) {}

  /// @brief Copy construct a face view.
  constexpr BaseFaceView( // NOLINT(google-explicit-constructor)
      const BaseFaceView<std::remove_const_t<Mesh>>& other) noexcept
      : BaseElementView<Mesh, FaceTag>(other) {}

  /// @brief Get the adjacent inner cell.
  [[nodiscard]] constexpr auto inner_cell() const noexcept {
    STORM_ASSERT_(detail_::face_inner_cell_ < this->adjacent_cells().size(),
                  "The face does not have an adjacent inner cell!");
    return this->adjacent_cells()[detail_::face_inner_cell_];
  }

  /// @brief Get the adjacent outer cell.
  [[nodiscard]] constexpr auto outer_cell() const noexcept {
    STORM_ASSERT_(detail_::face_outer_cell_ < this->adjacent_cells().size(),
                  "The face does not have an adjacent outer cell!");
    return this->adjacent_cells()[detail_::face_outer_cell_];
  }

  /// @brief Get the face area (or length in 2D).
  [[nodiscard]] constexpr auto area() const noexcept {
    return this->mesh_->face_area(this->index_);
  }

  /// @brief Get the face normal.
  [[nodiscard]] constexpr auto normal() const noexcept {
    return this->mesh_->face_normal(this->index_);
  }

  /// @brief Get the face center position.
  [[nodiscard]] constexpr auto center() const noexcept {
    return this->mesh_->face_center(this->index_);
  }

}; // class BaseFaceView

/// @brief Mesh face view.
/// @{
using FaceView = BaseFaceView<const Mesh>;
using MutableFaceView = BaseFaceView<Mesh>;
/// @}

/// @brief Base cell view.
template<class Mesh>
class BaseCellView final : public BaseElementView<Mesh, CellTag> {
public:

  /// @brief Construct base cell view.
  constexpr BaseCellView(Mesh& mesh, CellIndex index) noexcept
      : BaseElementView<Mesh, CellTag>(mesh, index) {}

  /// @brief Copy construct a cell view.
  constexpr BaseCellView( // NOLINT(google-explicit-constructor)
      const BaseCellView<std::remove_const_t<Mesh>>& other) noexcept
      : BaseElementView<Mesh, CellTag>(other) {}

  /// @brief Get the cell volume (or area in 2D).
  [[nodiscard]] constexpr auto volume() const noexcept {
    return this->mesh_->cell_volume(this->index_);
  }

  /// @brief Get cell center position.
  [[nodiscard]] constexpr auto center() const noexcept {
    return this->mesh_->cell_center(this->index_);
  }

}; // class BaseCellView

/// @brief Mesh cell view.
/// @{
using CellView = BaseCellView<const Mesh>;
using MutableCellView = BaseCellView<Mesh>;
/// @}

/// @brief Range of the @p mesh nodes
/// (or nodes with a @p node_mark, if present).
[[nodiscard]] constexpr auto
node_views(auto& mesh, std::same_as<NodeMark> auto... node_mark) noexcept {
  return mesh.nodes(node_mark...) |
         views::transform([&mesh](NodeIndex node_index) {
           return BaseNodeView(mesh, node_index);
         });
}

/// @brief Range of the @p mesh edges
/// (or edges with an @p edge_mark, if present).
[[nodiscard]] constexpr auto
edge_views(auto& mesh, std::same_as<EdgeMark> auto... edge_mark) noexcept {
  return mesh.edges(edge_mark...) |
         views::transform([&mesh](EdgeIndex edge_index) {
           return BaseEdgeView(mesh, edge_index);
         });
}

/// @brief Range of the @p mesh faces
/// (or faces with a @p face_mark, if present).
[[nodiscard]] constexpr auto
face_views(auto& mesh, std::same_as<FaceMark> auto... face_mark) noexcept {
  return mesh.faces(face_mark...) |
         views::transform([&mesh](FaceIndex face_index) {
           return BaseFaceView(mesh, face_index);
         });
}

/// @brief Range of the @p mesh cells
/// (or cells with a @p cell_mark, if present).
[[nodiscard]] constexpr auto
cell_views(auto& mesh, std::same_as<CellMark> auto... cell_mark) noexcept {
  return mesh.cells(cell_mark...) |
         views::transform([&mesh](CellIndex cell_index) {
           return BaseCellView(mesh, cell_index);
         });
}

/// @brief Range of the interior @p mesh nodes.
[[nodiscard]] constexpr auto int_node_views(auto& mesh) noexcept {
  return node_views(mesh, NodeMark{0});
}

/// @brief Range of the interior @p mesh edges.
[[nodiscard]] constexpr auto int_edge_views(auto& mesh) noexcept {
  return edge_views(mesh, EdgeMark{0});
}

/// @brief Range of the interior @p mesh faces.
[[nodiscard]] constexpr auto int_face_views(auto& mesh) noexcept {
  return face_views(mesh, FaceMark{0});
}

/// @brief Range of the interior @p mesh cells.
[[nodiscard]] constexpr auto int_cell_views(auto& mesh) noexcept {
  return cell_views(mesh, CellMark{0});
}

/// @brief Range of the boundary @p mesh nodes.
[[nodiscard]] constexpr auto bnd_node_views(auto& mesh) noexcept {
  return node_views(mesh) | views::drop(mesh.nodes(NodeMark{0}).size());
}

/// @brief Range of the boundary @p mesh edges.
[[nodiscard]] constexpr auto bnd_edge_views(auto& mesh) noexcept {
  return edge_views(mesh) | views::drop(mesh.edges(EdgeMark{0}).size());
}

/// @brief Range of the boundary @p mesh faces.
[[nodiscard]] constexpr auto bnd_face_views(auto& mesh) noexcept {
  return face_views(mesh) | views::drop(mesh.faces(FaceMark{0}).size());
}

template<class Mesh>
[[nodiscard]] constexpr auto bnd_face_cell_views(Mesh& mesh) noexcept {
  return bnd_face_views(mesh) | views::transform([](BaseFaceView<Mesh> face) {
           return std::pair(face.inner_cell(), face.outer_cell());
         });
}

template<class Mesh>
void for_each_bnd_face_cells(Mesh& mesh, auto&& func) noexcept {
  ForEach(bnd_face_views(mesh), [&](BaseFaceView<Mesh> face) {
    func(face.inner_cell(), face.outer_cell());
  });
}

/// @brief Range of the boundary @p mesh cells.
[[nodiscard]] constexpr auto bnd_cell_views(auto& mesh) noexcept {
  return cell_views(mesh) | views::drop(mesh.cells(CellMark{0}).size());
}

} // namespace Storm
