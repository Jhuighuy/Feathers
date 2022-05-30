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

#include <SkunkBase.hh>
#include <libFeathersUtils/Parallel.hh>

namespace feathers {

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
  Mesh* Mesh_;
  Index<Tag> Index_;

  template<class, class>
  friend class BaseElementView;

  // "NOLINT(...)" should not be here, this is due to a bug in clangd.
  BaseElementView( // NOLINT(cppcoreguidelines-pro-type-member-init)
    Mesh& mesh, Index<Tag> index) noexcept :
      Mesh_{&mesh}, Index_{index} {
    StormAssert(Index_ != npos);
  }

  template<class OtherMesh>
  BaseElementView( // NOLINT(google-explicit-constructor,cppcoreguidelines-pro-type-member-init)
    BaseElementView<OtherMesh, Tag> const& other) noexcept :
      Mesh_{other.Mesh_}, Index_{other.Index_} {
    StormAssert(Index_ != npos);
  }

public:

  /// @brief Cast to index operator.
  /// @{
  operator Index<Tag>() const noexcept {
    return Index_;
  }
  FEATHERS_DEPRECATED operator size_t() const noexcept {
    return static_cast<size_t>(Index_);
  }
  /// @}

  /// @brief Comparison operator.
  auto operator<=>(BaseElementView const& other) const noexcept {
    StormAssert(Mesh_ == other.Mesh_);
    return Index_ <=> other.Index_;
  }

  /// @brief Get mark. 
  Index<MarkTag<Tag>> mark() const noexcept {
    return Mesh_->mark(Index_);
  }

  /// @brief Get shape type.
  ShapeType shapeType() const noexcept {
    return Mesh_->shapeType(Index_);
  }

  /// @brief Get shape.
  std::unique_ptr<Element> shape() const {
    return Mesh_->shape(Index_);
  }

  /// @brief Ranges of the adjacent nodes.
  auto adjNodes() const noexcept {
    return Mesh_->adjNodes(Index_) |
      views::transform([&mesh = *Mesh_](NodeIndex nodeIndex) {
        return BaseNodeView(mesh, nodeIndex);
      });
  }

  /// @brief Ranges of the adjacent edges.
  auto adjEdges() const noexcept {
    return Mesh_->adjEdges(Index_) |
      views::transform([&mesh = *Mesh_](EdgeIndex edgeIndex) {
        return BaseEdgeView(mesh, edgeIndex);
      });
  }

  /// @brief Ranges of the adjacent faces.
  auto adjFaces() const noexcept {
    return Mesh_->adjFaces(Index_) |
      views::transform([&mesh = *Mesh_](FaceIndex faceIndex) {
        return BaseFaceView(mesh, faceIndex);
      });
  }

  /// @brief Ranges of the adjacent cells.
  auto adjCells() const noexcept {
    return Mesh_->adjCells(Index_) |
      views::transform([&mesh = *Mesh_](CellIndex cellIndex) {
        return BaseCellView(mesh, cellIndex);
      });
  }

  /// @brief Sequentially iterate through all the adjacent nodes.
  template<class Func>
  void forEachNode(Func&& func) const noexcept {
    ranges::for_each(adjNodes(), func);
  }

  /// @brief Sequentially iterate through all the adjacent edges.
  template<class Func>
  void forEachEdge(Func&& func) const noexcept {
    ranges::for_each(adjEdges(), func);
  }

  /// @brief Sequentially iterate through all the adjacent faces.
  /// @{
  template<class Func>
  void forEachFace(Func&& func) const noexcept {
    ranges::for_each(adjFaces(), func);
  }
  template<class Func>
  void forEachFaceCells(Func&& func) const noexcept {
    ranges::for_each(adjFaces(), [&](BaseFaceView<Mesh> face) {
      func(face.innerCell(), face.outerCell());
    });
  }
  /// @}

  /// @brief Sequentially iterate through all the adjacent cells.
  template<class Func>
  void forEachCell(Func&& func) const noexcept {
    ranges::for_each(adjCells(), func);
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
  vec3_t pos() const noexcept {
    return this->Mesh_->nodePos(this->Index_);
  }

  /// @brief Set node position @p pos.
  void setPos(vec3_t const& pos) const noexcept requires (!std::is_const_v<Mesh>) {
    this->Mesh_->setNodePos(this->Index_, pos);
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
    return this->Mesh_->edgeLen(this->Index_);
  }

  /// @brief Get edge direction. 
  vec3_t dir() const noexcept {
    return this->Mesh_->edgeDir(this->Index_);
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
  auto innerCell() const noexcept {
    StormAssert(this->adjCells().size() == 2);
    return this->adjCells()[FaceInnerCell_];
  }

  /// @brief Get connected outer cell. 
  auto outerCell() const noexcept {
    StormAssert(this->adjCells().size() == 2);
    return this->adjCells()[FaceOuterCell_];
  }

  /// @brief Get face area/length. 
  real_t area() const noexcept {
    return this->Mesh_->faceArea(this->Index_);
  }

  /// @brief Get face normal. 
  vec3_t normal() const noexcept {
    return this->Mesh_->faceNormal(this->Index_);
  }

  /// @brief Get face center position.
  vec3_t centerPos() const noexcept {
    return this->Mesh_->faceCenterPos(this->Index_);
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
    return this->Mesh_->cellVolume(this->Index_);
  }

  /// @brief Get cell center position.
  vec3_t centerPos() const noexcept {
    return this->Mesh_->cellCenterPos(this->Index_);
  }

}; // class BaseCellView<...>

/// @brief Mesh cell view.
/// @{
using CellView = BaseCellView<Mesh const>;
using MutableCellView = BaseCellView<Mesh>;
/// @}

/// @brief Range of the @p mesh nodes 
///   (or nodes with a @p nodeMark, if present).
auto nodeViews(auto& mesh, std::same_as<NodeMark> auto... nodeMark) noexcept {
  static_assert(sizeof...(nodeMark) <= 1);
  return mesh.nodes(nodeMark...) |
    views::transform([&mesh](NodeIndex nodeIndex) {
      return BaseNodeView(mesh, nodeIndex);
    });
}

/// @brief Range of the @p mesh edges 
///   (or edges with an @p edgeMark, if present).
auto edgeViews(auto& mesh, std::same_as<EdgeMark> auto... edgeMark) noexcept {
  static_assert(sizeof...(edgeMark) <= 1);
  return mesh.edges(edgeMark...) |
    views::transform([&mesh](EdgeIndex edgeIndex) {
      return BaseEdgeView(mesh, edgeIndex);
    });
}

/// @brief Range of the @p mesh faces 
///   (or faces with a @p faceMark, if present).
auto faceViews(auto& mesh, std::same_as<FaceMark> auto... faceMark) noexcept {
  static_assert(sizeof...(faceMark) <= 1);
  return mesh.faces(faceMark...) |
    views::transform([&mesh](FaceIndex faceIndex) {
      return BaseFaceView(mesh, faceIndex);
    });
}

/// @brief Range of the @p mesh nodes 
///   (or nodes with a @p cellMark, if present).
auto cellViews(auto& mesh, std::same_as<CellMark> auto... cellMark) noexcept {
  static_assert(sizeof...(cellMark) <= 1);
  return mesh.cells(cellMark...) |
    views::transform([&mesh](CellIndex cellIndex) {
      return BaseCellView(mesh, cellIndex);
    });
}

/// @brief Range of the interior @p mesh nodes.
auto intNodeViews(auto& mesh) noexcept {
  return nodeViews(mesh, NodeMark{0});
}

/// @brief Range of the interior @p mesh edges.
auto intEdgeViews(auto& mesh) noexcept {
  return edgeViews(mesh, EdgeMark{0});
}

/// @brief Range of the interior @p mesh faces.
auto intFaceViews(auto& mesh) noexcept {
  return faceViews(mesh, FaceMark{0});
}

/// @brief Range of the interior @p mesh cells.
auto intCellViews(auto& mesh) noexcept {
  return cellViews(mesh, CellMark{0});
}

/// @brief Range of the boundary @p mesh nodes.
auto bndNodeViews(auto& mesh) noexcept {
  return nodeViews(mesh) | views::drop(mesh.nodes(NodeMark{0}).size());
}

/// @brief Range of the boundary @p mesh edges.
auto bndEdgeViews(auto& mesh) noexcept {
  return edgeViews(mesh) | views::drop(mesh.edges(EdgeMark{0}).size());
}

/// @brief Range of the boundary @p mesh faces.
auto bndFaceViews(auto& mesh) noexcept {
  return faceViews(mesh) | views::drop(mesh.faces(FaceMark{0}).size());
}

template<class Mesh>
auto bndFaceCellViews(Mesh& mesh) noexcept {
  return faceViews(mesh) |
    views::drop(mesh.faceIndices(FaceMark{0}).size()) |
    views::transform([](BaseFaceView<Mesh> face) {
      return std::pair(face.innerCell(), face.outerCell());
    });
}

template<class Mesh>
void forEachBndFaceCells(Mesh& mesh, auto&& func) noexcept {
  ForEach(bndFaceViews(mesh), [&](BaseFaceView<Mesh> face) {
    func(face.innerCell(), face.outerCell());
  });
}

/// @brief Range of the boundary @p mesh cells.
auto bndCellViews(auto& mesh) noexcept {
  return cellViews(mesh) | views::drop(mesh.cells(CellMark{0}).size());
}

} // namespace feathers
