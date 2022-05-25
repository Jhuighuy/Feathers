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
    return Mesh_->Mark(Index_);
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
    return Mesh_->adjNodeIndices(Index_) |
      views::transform([&mesh = *Mesh_](NodeIndex nodeIndex) {
        return BaseNodeView(mesh, nodeIndex);
      });
  }

  /// @brief Ranges of the adjacent edges.
  auto adjEdges() const noexcept {
    return Mesh_->adjEdgeIndices(Index_) |
      views::transform([&mesh = *Mesh_](EdgeIndex edgeIndex) {
        return BaseEdgeView(mesh, edgeIndex);
      });
  }

  /// @brief Ranges of the adjacent faces.
  auto adjFaces() const noexcept {
    return Mesh_->adjFaceIndices(Index_) |
      views::transform([&mesh = *Mesh_](FaceIndex faceIndex) {
        return BaseFaceView(mesh, faceIndex);
      });
  }

  /// @brief Ranges of the adjacent cells.
  auto adjCells() const noexcept {
    return Mesh_->adjCellIndices(Index_) |
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
    return this->adjCells()[FaceInnerCell_];
  }

  /// @brief Get connected outer cell. 
  auto outerCell() const noexcept {
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

/// @brief Range of the @p mesh nodes \
///   (or nodes with a @p nodeMark, if present).
template<class Mesh, std::same_as<NodeMark>... NodeMark_>
  requires (sizeof...(NodeMark_) <= 1)
auto nodeViews(Mesh& mesh, NodeMark_... nodeMark) noexcept {
  return mesh.nodeIndices(nodeMark...) |
    views::transform([&mesh](NodeIndex nodeIndex) {
      return BaseNodeView(mesh, nodeIndex);
    });
}

/// @brief Range of the @p mesh edges \
///   (or edges with an @p edgeMark, if present).
template<class Mesh, std::same_as<EdgeMark>... EdgeMark_>
  requires (sizeof...(EdgeMark_) <= 1)
auto edgeViews(Mesh& mesh, EdgeMark_... edgeMark) noexcept {
  return mesh.edgeIndices(edgeMark...) |
    views::transform([&mesh](EdgeIndex edgeIndex) {
      return BaseEdgeView(mesh, edgeIndex);
    });
}

/// @brief Range of the @p mesh faces \
///   (or faces with a @p faceMark, if present).
template<class Mesh, std::same_as<FaceMark>... FaceMark_>
  requires (sizeof...(FaceMark_) <= 1)
auto faceViews(Mesh& mesh, FaceMark_... faceMark) noexcept {
  return mesh.faceIndices(faceMark...) |
    views::transform([&mesh](FaceIndex faceIndex) {
      return BaseFaceView(mesh, faceIndex);
    });
}

/// @brief Range of the @p mesh nodes \
///   (or nodes with a @p cellMark, if present).
template<class Mesh, std::same_as<CellMark>... CellMark_>
  requires (sizeof...(CellMark_) <= 1)
auto cellViews(Mesh& mesh, CellMark_... cellMark) noexcept {
  return mesh.cellIndices(cellMark...) |
    views::transform([&mesh](CellIndex cellIndex) {
      return BaseCellView(mesh, cellIndex);
    });
}

/// @brief Range of the interior @p mesh nodes.
template<class Mesh>
auto intNodeViews(Mesh& mesh) noexcept {
  return nodeViews(mesh, NodeMark{0});
}

/// @brief Range of the interior @p mesh edges.
template<class Mesh>
auto intEdgeViews(Mesh& mesh) noexcept {
  return edgeViews(mesh, EdgeMark{0});
}

/// @brief Range of the interior @p mesh faces.
template<class Mesh>
auto intFaceViews(Mesh& mesh) noexcept {
  return faceViews(mesh, FaceMark{0});
}

/// @brief Range of the interior @p mesh cells.
template<class Mesh>
auto intCellViews(Mesh& mesh) noexcept {
  return cellViews(mesh, CellMark{0});
}

/// @brief Range of the boundary @p mesh nodes.
template<class Mesh>
auto bndNodeViews(Mesh& mesh) noexcept {
  return nodeViews(mesh) | views::drop(mesh.nodeIndices(NodeMark{0}).size());
}

/// @brief Range of the boundary @p mesh edges.
template<class Mesh>
auto bndEdgeViews(Mesh& mesh) noexcept {
  return edgeViews(mesh) | views::drop(mesh.edgeIndices(EdgeMark{0}).size());
}

/// @brief Range of the boundary @p mesh faces.
template<class Mesh>
auto bndFaceViews(Mesh& mesh) noexcept {
  return faceViews(mesh) | views::drop(mesh.faceIndices(FaceMark{0}).size());
}

template<class Mesh>
auto bndFaceCellViews(Mesh& mesh) noexcept {
  return faceViews(mesh) |
    views::drop(mesh.faceIndices(FaceMark{0}).size()) |
    views::transform([](BaseFaceView<Mesh> face) {
      return std::make_pair(face.innerCell(), face.outerCell());
    });
}

template<class Mesh, class Func>
void forEachBndFaceCells(Mesh& mesh, Func func) noexcept {
  ForEach(bndFaceViews(mesh), [&](BaseFaceView<Mesh> face) {
    func(face.innerCell(), face.outerCell());
  });
}

/// @brief Range of the boundary @p mesh cells.
template<class Mesh>
auto bndCellViews(Mesh& mesh) noexcept {
  return cellViews(mesh) | views::drop(mesh.cellIndices(CellMark{0}).size());
}

} // namespace feathers
