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

  // "NOLINT" is a bug here.
  BaseElementView( // NOLINT(cppcoreguidelines-pro-type-member-init)
    Mesh& mesh, Index<Tag> index) noexcept :
      Mesh_{&mesh}, Index_{index} {
    StormAssert(Index_ != npos);
  }

  // "NOLINT" for member initialization is a bug here.
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
  Index<MarkTag<Tag>> Mark() const noexcept {
    return Mesh_->Mark(Index_);
  }

  /// @brief Get Shape.
  ShapeType Shape() const noexcept {
    return Mesh_->Shape(Index_);
  }

  /// @brief Get element object. 
  std::unique_ptr<const Element> get_element_object() const {
    return Mesh_->get_object(Index_);
  }

  /// @brief Ranges of the adjacent nodes.
  auto Nodes() const noexcept {
    return Mesh_->AdjacentNodes(Index_) |
      views::transform([&mesh = *Mesh_](NodeIndex nodeIndex) {
        return BaseNodeView(mesh, nodeIndex);
      });
  }

  /// @brief Ranges of the adjacent edges.
  auto Edges() const noexcept {
    return Mesh_->AdjacentEdges(Index_) |
      views::transform([&mesh = *Mesh_](EdgeIndex edgeIndex) {
        return BaseEdgeView(mesh, edgeIndex);
      });
  }

  /// @brief Ranges of the adjacent faces.
  auto Faces() const noexcept {
    return Mesh_->AdjacentFaces(Index_) |
      views::transform([&mesh = *Mesh_](FaceIndex faceIndex) {
        return BaseFaceView(mesh, faceIndex);
      });
  }

  /// @brief Ranges of the adjacent cells.
  auto Cells() const noexcept {
    return Mesh_->AdjacentCells(Index_) |
      views::transform([&mesh = *Mesh_](CellIndex cellIndex) {
        return BaseCellView(mesh, cellIndex);
      });
  }

  /// @brief Iterate through all connected nodes. 
  template<class Func>
  void ForEachNode(Func&& func) const noexcept {
    ranges::for_each(Nodes(), func);
  }

  /// @brief Iterate through all connected edges. 
  template<class Func>
  void ForEachEdge(Func&& func) const noexcept {
    ranges::for_each(Edges(), func);
  }

  /// @brief Iterate through all connected faces. 
  /// @{
  template<class Func>
  void ForEachFace(Func&& func) const noexcept {
    ranges::for_each(Faces(), func);
  }
  template<class Func>
  void ForEachFaceCells(Func&& func) const noexcept {
    ranges::for_each(Faces(), [&](BaseFaceView<Mesh> face) {
      func(face.InnerCell(), face.OuterCell());
    });
  }
  /// @}

  /// @brief Iterate through all connected cells. 
  template<class Func>
  void ForEachCell(Func&& func) const noexcept {
    ranges::for_each(Cells(), func);
  }

}; // class BaseElementView

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base Node view.
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
  vec3_t Pos() const noexcept {
    return this->Mesh_->NodePos(this->Index_);
  }

  /// @brief Set node position @p pos.
  void SetPos(vec3_t const& pos) const noexcept 
                                 requires (!std::is_const_v<Mesh>) {
    this->Mesh_->SetNodePos(this->Index_, pos);
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
  real_t Len() const noexcept {
    return this->Mesh_->EdgeLen(this->Index_);
  }

  /// @brief Get edge direction. 
  vec3_t Dir() const noexcept {
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
  BaseFaceView(Mesh& mesh, FaceIndex index) noexcept :
    BaseElementView<Mesh, FaceTag>(mesh, index) {
  }

  /// @brief Copy constructor.
  BaseFaceView( // NOLINT(google-explicit-constructor)
      BaseFaceView<std::remove_const_t<Mesh>> const& other) noexcept :
    BaseElementView<Mesh, FaceTag>(other) {
  }

  /// @brief Get connected inner cell. 
  auto InnerCell() const noexcept {
    return this->Cells()[FaceInnerCell_];
  }

  /// @brief Get connected outer cell. 
  auto OuterCell() const noexcept {
    return this->Cells()[FaceOuterCell_];
  }

  /// @brief Get face area/length. 
  real_t Area() const noexcept {
    return this->Mesh_->FaceArea(this->Index_);
  }

  /// @brief Get face normal. 
  vec3_t Normal() const noexcept {
    return this->Mesh_->FaceNormal(this->Index_);
  }

  /// @brief Get face center position.
  vec3_t CenterPos() const noexcept {
    return this->Mesh_->FaceCenterPos(this->Index_);
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
  real_t Volume() const noexcept {
    return this->Mesh_->CellVolume(this->Index_);
  }

  /// @brief Get cell center position.
  vec3_t CenterPos() const noexcept {
    return this->Mesh_->CellCenterPos(this->Index_);
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
auto NodeViews(Mesh& mesh, NodeMark_... nodeMark) noexcept {
  return mesh.Nodes(nodeMark...) |
    views::transform([&mesh](NodeIndex nodeIndex) {
      return BaseNodeView(mesh, nodeIndex);
    });
}

/// @brief Range of the @p mesh edges \
///   (or edges with an @p edgeMark, if present).
template<class Mesh, std::same_as<EdgeMark>... EdgeMark_>
  requires (sizeof...(EdgeMark_) <= 1)
auto EdgeViews(Mesh& mesh, EdgeMark_... edgeMark) noexcept {
  return mesh.Edges(edgeMark...) |
    views::transform([&mesh](EdgeIndex edgeIndex) {
      return BaseEdgeView(mesh, edgeIndex);
    });
}

/// @brief Range of the @p mesh faces \
///   (or faces with a @p faceMark, if present).
template<class Mesh, std::same_as<FaceMark>... FaceMark_>
  requires (sizeof...(FaceMark_) <= 1)
auto FaceViews(Mesh& mesh, FaceMark_... faceMark) noexcept {
  return mesh.Faces(faceMark...) |
    views::transform([&mesh](FaceIndex faceIndex) {
      return BaseFaceView(mesh, faceIndex);
    });
}

/// @brief Range of the @p mesh nodes \
///   (or nodes with a @p cellMark, if present).
template<class Mesh, std::same_as<CellMark>... CellMark_>
  requires (sizeof...(CellMark_) <= 1)
auto CellViews(Mesh& mesh, CellMark_... cellMark) noexcept {
  return mesh.Cells(cellMark...) |
    views::transform([&mesh](CellIndex cellIndex) {
      return BaseCellView(mesh, cellIndex);
    });
}

/// @brief Range of the interior @p mesh nodes.
template<class Mesh>
auto InteriorNodeViews(Mesh& mesh) noexcept {
  return NodeViews(mesh, NodeMark{0});
}

/// @brief Range of the interior @p mesh edges.
template<class Mesh>
auto InteriorEdgeViews(Mesh& mesh) noexcept {
  return EdgeViews(mesh, EdgeMark{0});
}

/// @brief Range of the interior @p mesh nodes.
template<class Mesh>
auto InteriorFaceViews(Mesh& mesh) noexcept {
  return FaceViews(mesh, FaceMark{0});
}

/// @brief Range of the interior @p mesh nodes.
template<class Mesh>
auto InteriorCellViews(Mesh& mesh) noexcept {
  return CellViews(mesh, CellMark{0});
}

/// @brief Range of the boundary @p mesh nodes.
template<class Mesh>
auto BoundaryNodeViews(Mesh& mesh) noexcept {
  return NodeViews(mesh) | views::drop(mesh.Nodes(NodeMark{0}).size());
}

/// @brief Range of the boundary @p mesh edges.
template<class Mesh>
auto BoundaryEdgeViews(Mesh& mesh) noexcept {
  return EdgeViews(mesh) | views::drop(mesh.Edges(EdgeMark{0}).size());
}

/// @brief Range of the boundary @p mesh faces.
template<class Mesh>
auto BoundaryFaceViews(Mesh& mesh) noexcept {
  return FaceViews(mesh) | views::drop(mesh.Faces(FaceMark{0}).size());
}

template<class Mesh>
auto BoundaryFaceCellViews(Mesh& mesh) noexcept {
  return FaceViews(mesh) |
    views::drop(mesh.Faces(FaceMark{0}).size()) |
    views::transform([](BaseFaceView<Mesh> face) {
      return std::make_pair(face.InnerCell(), face.OuterCell());
    });
}

/// @brief Iterate through all boundary @p mesh faces.
template<class Mesh, class Func>
FEATHERS_DEPRECATED void ForEachBoundaryFaceCells(Mesh& mesh, Func func) noexcept {
  ForEach(BoundaryFaceViews(mesh), [&](BaseFaceView<Mesh> face) {
    func(face.InnerCell(), face.OuterCell());
  });
}

/// @brief Range of the boundary @p mesh cells.
template<class Mesh>
auto BoundaryCellViews(Mesh& mesh) noexcept {
  return CellViews(mesh) | views::drop(mesh.Cells(CellMark{0}).size());
}

} // namespace feathers
