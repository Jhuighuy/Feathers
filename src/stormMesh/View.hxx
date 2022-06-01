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
      Mesh& mesh, Index<Tag> index) noexcept
      : Mesh_{&mesh}, Index_{index} {
    StormAssert(Index_ != npos);
  }

  template<class OtherMesh>
  BaseElementView( // NOLINT(google-explicit-constructor,cppcoreguidelines-pro-type-member-init)
      BaseElementView<OtherMesh, Tag> const& other) noexcept
      : Mesh_{other.Mesh_}, Index_{other.Index_} {
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
    StormAssert(Mesh_ == other.Mesh_ &&
                "can not compare element views on the different meshes");
    return Index_ <=> other.Index_;
  }

  /// @brief Get Mark.
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

  /// @brief Ranges of the adjacent Nodes.
  auto AdjacentNodes() const noexcept {
    return Mesh_->AdjacentNodes(Index_) |
           views::transform([&mesh = *Mesh_](NodeIndex nodeIndex) {
             return BaseNodeView(mesh, nodeIndex);
           });
  }

  /// @brief Ranges of the adjacent Edges.
  auto AdjacentEdges() const noexcept {
    return Mesh_->AdjacentEdges(Index_) |
           views::transform([&mesh = *Mesh_](EdgeIndex edgeIndex) {
             return BaseEdgeView(mesh, edgeIndex);
           });
  }

  /// @brief Ranges of the adjacent Faces.
  auto AdjacentFaces() const noexcept {
    return Mesh_->AdjacentFaces(Index_) |
           views::transform([&mesh = *Mesh_](FaceIndex faceIndex) {
             return BaseFaceView(mesh, faceIndex);
           });
  }

  /// @brief Ranges of the adjacent Cells.
  auto AdjacentCells() const noexcept {
    return Mesh_->AdjacentCells(Index_) |
           views::transform([&mesh = *Mesh_](CellIndex cellIndex) {
             return BaseCellView(mesh, cellIndex);
           });
  }

  /// @brief Sequentially iterate through all the adjacent nodes.
  void ForEachNode(auto&& func) const noexcept {
    ranges::for_each(AdjacentNodes(), func);
  }

  /// @brief Sequentially iterate through all the adjacent edges.
  void ForEachEdge(auto&& func) const noexcept {
    ranges::for_each(AdjacentEdges(), func);
  }

  /// @brief Sequentially iterate through all the adjacent faces.
  /// @{
  void ForEachFace(auto&& func) const noexcept {
    ranges::for_each(AdjacentFaces(), func);
  }
  void ForEachFaceCells(auto&& func) const noexcept {
    ranges::for_each(AdjacentFaces(), [&](BaseFaceView<Mesh> face) {
      func(face.InnerCell(), face.OuterCell());
    });
  }
  /// @}

  /// @brief Sequentially iterate through all the adjacent cells.
  void ForEachCell(auto&& func) const noexcept {
    ranges::for_each(AdjacentCells(), func);
  }

}; // class BaseElementView

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base node view.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Mesh>
class BaseNodeView final : public BaseElementView<Mesh, NodeTag> {
public:

  /// @brief Construct base node view.
  BaseNodeView(Mesh& mesh, NodeIndex index) noexcept
      : BaseElementView<Mesh, NodeTag>(mesh, index) {}

  /// @brief Copy constructor.
  BaseNodeView( // NOLINT(google-explicit-constructor)
      BaseNodeView<std::remove_const_t<Mesh>> const& other) noexcept
      : BaseElementView<Mesh, NodeTag>(other) {}

  /// @brief Get node position.
  auto Coords() const noexcept {
    return this->Mesh_->NodeCoords(this->Index_);
  }

  /// @brief Set node position @p Coords.
  void SetCoords(auto const& coords) const noexcept
    requires(!std::is_const_v<Mesh>)
  {
    this->Mesh_->SetNodeCoords(this->Index_, coords);
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
  BaseEdgeView(Mesh& mesh, EdgeIndex index) noexcept
      : BaseElementView<Mesh, EdgeTag>(mesh, index) {}

  /// @brief Copy constructor.
  BaseEdgeView( // NOLINT(google-explicit-constructor)
      BaseEdgeView<std::remove_const_t<Mesh>> const& other) noexcept
      : BaseElementView<Mesh, EdgeTag>(other) {}

  /// @brief Get edge length.
  auto Len() const noexcept {
    return this->Mesh_->EdgeLen(this->Index_);
  }

  /// @brief Get edge direction.
  auto Dir() const noexcept {
    return this->Mesh_->EdgeDir(this->Index_);
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
  BaseFaceView(Mesh& mesh, FaceIndex index) noexcept
      : BaseElementView<Mesh, FaceTag>(mesh, index) {}

  /// @brief Copy constructor.
  BaseFaceView( // NOLINT(google-explicit-constructor)
      BaseFaceView<std::remove_const_t<Mesh>> const& other) noexcept
      : BaseElementView<Mesh, FaceTag>(other) {}

  /// @brief Get the connected inner cell.
  auto InnerCell() const noexcept {
    StormAssert(FaceInnerCell_ < this->AdjacentCells().size() &&
                "the face does not have an adjacent inner cell");
    return this->AdjacentCells()[FaceInnerCell_];
  }

  /// @brief Get the connected outer cell.
  auto OuterCell() const noexcept {
    StormAssert(FaceOuterCell_ < this->AdjacentCells().size() &&
                "the face does not have an adjacent outer cell");
    return this->AdjacentCells()[FaceOuterCell_];
  }

  /// @brief Get face area/length.
  auto Area() const noexcept {
    return this->Mesh_->FaceArea(this->Index_);
  }

  /// @brief Get face Normal.
  auto Normal() const noexcept {
    return this->Mesh_->FaceNormal(this->Index_);
  }

  /// @brief Get face center position.
  auto Center() const noexcept {
    return this->Mesh_->FaceCenter(this->Index_);
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
  BaseCellView(Mesh& mesh, CellIndex index) noexcept
      : BaseElementView<Mesh, CellTag>(mesh, index) {}

  /// @brief Copy constructor.
  BaseCellView( // NOLINT(google-explicit-constructor)
      BaseCellView<std::remove_const_t<Mesh>> const& other) noexcept
      : BaseElementView<Mesh, CellTag>(other) {}

  /// @brief Get cell volume/area/length.
  auto Volume() const noexcept {
    return this->Mesh_->CellVolume(this->Index_);
  }

  /// @brief Get cell center position.
  auto Center() const noexcept {
    return this->Mesh_->CellCenter(this->Index_);
  }

}; // class BaseCellView<...>

/// @brief Mesh cell view.
/// @{
using CellView = BaseCellView<Mesh const>;
using MutableCellView = BaseCellView<Mesh>;
/// @}

/// @brief Range of the @p mesh nodes
///   (or nodes with a @p nodeMark, if present).
auto NodeViews(auto& mesh, std::same_as<NodeMark> auto... nodeMark) noexcept {
  static_assert(sizeof...(nodeMark) <= 1);
  return mesh.Nodes(nodeMark...) |
         views::transform([&mesh](NodeIndex nodeIndex) {
           return BaseNodeView(mesh, nodeIndex);
         });
}

/// @brief Range of the @p mesh edges
///   (or edges with an @p edgeMark, if present).
auto EdgeViews(auto& mesh, std::same_as<EdgeMark> auto... edgeMark) noexcept {
  static_assert(sizeof...(edgeMark) <= 1);
  return mesh.Edges(edgeMark...) |
         views::transform([&mesh](EdgeIndex edgeIndex) {
           return BaseEdgeView(mesh, edgeIndex);
         });
}

/// @brief Range of the @p mesh faces
///   (or faces with a @p faceMark, if present).
auto FaceViews(auto& mesh, std::same_as<FaceMark> auto... faceMark) noexcept {
  static_assert(sizeof...(faceMark) <= 1);
  return mesh.Faces(faceMark...) |
         views::transform([&mesh](FaceIndex faceIndex) {
           return BaseFaceView(mesh, faceIndex);
         });
}

/// @brief Range of the @p mesh cells
///   (or cells with a @p cellMark, if present).
auto CellViews(auto& mesh, std::same_as<CellMark> auto... cellMark) noexcept {
  static_assert(sizeof...(cellMark) <= 1);
  return mesh.Cells(cellMark...) |
         views::transform([&mesh](CellIndex cellIndex) {
           return BaseCellView(mesh, cellIndex);
         });
}

/// @brief Range of the interior @p mesh nodes.
auto IntNodeViews(auto& mesh) noexcept {
  return NodeViews(mesh, NodeMark{0});
}

/// @brief Range of the interior @p mesh edges.
auto IntEdgeViews(auto& mesh) noexcept {
  return EdgeViews(mesh, EdgeMark{0});
}

/// @brief Range of the interior @p mesh faces.
auto IntFaceViews(auto& mesh) noexcept {
  return FaceViews(mesh, FaceMark{0});
}

/// @brief Range of the interior @p mesh cells.
auto IntCellViews(auto& mesh) noexcept {
  return CellViews(mesh, CellMark{0});
}

/// @brief Range of the boundary @p mesh nodes.
auto BndNodeViews(auto& mesh) noexcept {
  return NodeViews(mesh) | views::drop(mesh.Nodes(NodeMark{0}).size());
}

/// @brief Range of the boundary @p mesh edges.
auto BndEdgeViews(auto& mesh) noexcept {
  return EdgeViews(mesh) | views::drop(mesh.Edges(EdgeMark{0}).size());
}

/// @brief Range of the boundary @p mesh faces.
auto BndFaceViews(auto& mesh) noexcept {
  return FaceViews(mesh) | views::drop(mesh.Faces(FaceMark{0}).size());
}

template<class Mesh>
auto BndFaceCellViews(auto& mesh) noexcept {
  return BndFaceViews(mesh) | views::transform([](BaseFaceView<Mesh> face) {
           return std::pair(face.InnerCell(), face.OuterCell());
         });
}

template<class Mesh>
void ForEachBndFaceCells(Mesh& mesh, auto&& func) noexcept {
  ForEach(BndFaceViews(mesh), [&](BaseFaceView<Mesh> face) {
    func(face.InnerCell(), face.OuterCell());
  });
}

/// @brief Range of the boundary @p mesh Cells.
auto BndCellViews(auto& mesh) noexcept {
  return CellViews(mesh) | views::drop(mesh.Cells(CellMark{0}).size());
}

} // namespace Storm
