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

#include <SkunkBase.hh>
#include <libFeathersUtils/Parallel.hh>

#include <range/v3/all.hpp>

namespace feathers {

template<class> class BaseNodeReference;
template<class> class BaseEdgeReference;
template<class> class BaseFaceReference;
template<class> class BaseCellReference;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base element reference.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Iterator, class Mesh, class Tag>
class BaseElementReference {
protected:
  Mesh* Mesh_;
  Index<Tag> Index_;

  template<class, class, class>
  friend class BaseElementReference;

  BaseElementReference( // NOLINT(cppcoreguidelines-pro-type-member-init)
    Mesh& mesh, Index<Tag> index) noexcept :
      Mesh_{&mesh}, Index_{index} {
    StormAssert(Index_ != npos);
  }

  template<class OtherIterator, class OtherMesh>
  BaseElementReference( // NOLINT(google-explicit-constructor,cppcoreguidelines-pro-type-member-init)
    BaseElementReference<OtherIterator, OtherMesh, Tag> const& other) noexcept :
      Mesh_{other.Mesh_}, Index_{other.Index_} {
    StormAssert(Index_ != npos);
  }

public:

  /// @brief Cast to index operator. 
  operator Index<Tag>() const noexcept {
    return Index_;
  }
  FEATHERS_DEPRECATED operator size_t() const noexcept {
    return static_cast<size_t>(Index_);
  }

  /// @brief Comparison operator.
  auto operator<=>(Iterator const& other) const noexcept {
    StormAssert(Mesh_ == other.Mesh_);
    return Index_ <=> other.Index_;
  }

  /// @brief Get mark. 
  Index<MarkTag_<Tag>> Mark() const noexcept {
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

  /// @brief Number of nodes in the element. 
  size_t NumNodes() const noexcept {
    return std::size(Mesh_->AdjacentNodes(Index_));
  }

  /// @brief Number of edges in the element. 
  size_t NumEdges() const noexcept {
    return std::size(Mesh_->AdjacentEdges(Index_));
  }

  /// @brief Number of faces in the element. 
  size_t NumFaces() const noexcept {
    return std::size(Mesh_->AdjacentFaces(Index_));
  }

  /// @brief Number of cells in the element. 
  size_t NumCells() const noexcept {
    return std::size(Mesh_->AdjacentCells(Index_));
  }

  /// @brief Ranges of the adjacent nodes.
  auto Nodes() const noexcept {
    return Mesh_->AdjacentNodes(Index_) |
      views::transform([&mesh = *Mesh_](NodeIndex nodeIndex) {
        return BaseNodeReference<Mesh>(mesh, nodeIndex);
      });
  }

  /// @brief Ranges of the adjacent edges.
  auto Edges() const noexcept {
    return Mesh_->AdjacentEdges(Index_) |
      views::transform([&mesh = *Mesh_](EdgeIndex edgeIndex) {
        return BaseEdgeReference<Mesh>(mesh, edgeIndex);
      });
  }

  /// @brief Ranges of the adjacent faces.
  auto Faces() const noexcept {
    return Mesh_->AdjacentFaces(Index_) |
      views::transform([&mesh = *Mesh_](FaceIndex faceIndex) {
        return BaseFaceReference<Mesh>(mesh, faceIndex);
      });
  }

  /// @brief Ranges of the adjacent cells.
  auto Cells() const noexcept {
    return Mesh_->AdjacentCells(Index_) |
      views::transform([&mesh = *Mesh_](CellIndex cellIndex) {
        return BaseCellReference<Mesh>(mesh, cellIndex);
      });
  }

  /// @brief Get the adjacent node at @p nodeLocal.
  auto Node(size_t nodeLocal) const noexcept {
    StormAssert(nodeLocal < NumNodes());
    return BaseNodeReference<Mesh>(
      *Mesh_, Mesh_->AdjacentNodes(Index_)[nodeLocal]);
  }

  /// @brief Get the adjacent edge at @p edgeLocal.
  auto Edge(size_t edgeLocal) const noexcept {
    StormAssert(edgeLocal < NumEdges());
    return BaseEdgeReference<Mesh>(
      *Mesh_, Mesh_->AdjacentEdges(Index_)[edgeLocal]);
  }

  /// @brief Get the adjacent face at @p faceLocal.
  auto Face(size_t faceLocal) const noexcept {
    StormAssert(faceLocal < NumFaces());
    return BaseFaceReference<Mesh>(
      *Mesh_, Mesh_->AdjacentFaces(Index_)[faceLocal]);
  }

  /// @brief Get the adjacent cell at @p cellLocal.
  auto Cell(size_t cellLocal) const noexcept {
    StormAssert(cellLocal < NumCells());
    return BaseCellReference<Mesh>(
      *Mesh_, Mesh_->AdjacentCells(Index_)[cellLocal]);
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
    ranges::for_each(Faces(), [&](BaseFaceReference<Mesh> face) {
      func(face.InnerCell(), face.OuterCell());
    });
  }
  /// @}

  /// @brief Iterate through all connected cells. 
  template<class Func>
  void ForEachCell(Func&& func) const noexcept {
    ranges::for_each(Cells(), func);
  }

}; // class BaseElementReference

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base Node reference.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Mesh>
class BaseNodeReference final :
  public BaseElementReference<BaseNodeReference<Mesh>, Mesh, NodeTag_> {
public:

  /// @brief Construct base Node reference.
  explicit BaseNodeReference(Mesh& mesh, NodeIndex index) noexcept :
    BaseElementReference<BaseNodeReference<Mesh>, Mesh, NodeTag_>(mesh, index) {
  }

  /// @brief Copy (make const) constructor. 
  BaseNodeReference( // NOLINT(google-explicit-constructor)
      BaseNodeReference<std::remove_const_t<Mesh>> const& other) noexcept :
    BaseElementReference<BaseNodeReference<Mesh>, Mesh, NodeTag_>(other) {
  }

  /// @brief Get node position.
  vec3_t Pos() const noexcept {
    return this->Mesh_->NodePos(this->Index_);
  }

  /// @brief Set node position @p pos.
  void SetPos(vec3_t const& pos) const noexcept requires(!std::is_const_v<Mesh>) {
    this->Mesh_->SetNodePos(this->Index_, pos);
  }

}; // class BaseNodeReference<...>

/// @brief Mesh nodes random-access reference.
/// @{
using NodeRef = BaseNodeReference<Mesh const>;
using NodeMutableRef = BaseNodeReference<Mesh>;
/// @}

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base edge reference.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Mesh>
class BaseEdgeReference final :
  public BaseElementReference<BaseEdgeReference<Mesh>, Mesh, EdgeTag_> {
public:

  /// @brief Construct base edge reference.
  explicit BaseEdgeReference(Mesh& mesh, EdgeIndex index) noexcept :
    BaseElementReference<BaseEdgeReference<Mesh>, Mesh, EdgeTag_>(mesh, index) {
  }

  /// @brief Copy (make const) constructor. 
  BaseEdgeReference( // NOLINT(google-explicit-constructor)
      BaseEdgeReference<std::remove_const_t<Mesh>> const& other) noexcept :
    BaseElementReference<BaseEdgeReference<Mesh>, Mesh, EdgeTag_>(other) {
  }

  /// @brief Get edge length. 
  real_t Len() const noexcept {
    return this->Mesh_->EdgeLen(this->Index_);
  }

  /// @brief Get edge direction. 
  vec3_t Dir() const noexcept {
    return this->Mesh_->EdgeDir(this->Index_);
  }

}; // class BaseEdgeReference<...>

/// @brief Mesh edges random-access reference.
/// @{
using EdgeRef = BaseEdgeReference<Mesh const>;
using EdgeMutableRef = BaseEdgeReference<Mesh>;
/// @}

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base face reference.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Mesh>
class BaseFaceReference final :
  public BaseElementReference<BaseFaceReference<Mesh>, Mesh, FaceTag_> {
public:

  /// @brief Construct base face reference.
  explicit BaseFaceReference(Mesh& mesh, FaceIndex index) noexcept :
    BaseElementReference<BaseFaceReference<Mesh>, Mesh, FaceTag_>(mesh, index) {
  }

  /// @brief Copy (make const) constructor. 
  BaseFaceReference( // NOLINT(google-explicit-constructor)
      BaseFaceReference<std::remove_const_t<Mesh>> const& other) noexcept :
    BaseElementReference<BaseFaceReference<Mesh>, Mesh, FaceTag_>(other) {
  }

  /// @brief Get connected inner cell. 
  auto InnerCell() const noexcept {
    return this->Cell(FaceInnerCell_);
  }

  /// @brief Get connected outer cell. 
  auto OuterCell() const noexcept {
    return this->Cell(FaceOuterCell_);
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

}; // class BaseFaceReference<...>

/// @brief Mesh faces random-access reference.
/// @{
using FaceRef = BaseFaceReference<Mesh const>;
using FaceMutableRef = BaseFaceReference<Mesh>;
/// @}

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base cell reference.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Mesh>
class BaseCellReference final :
  public BaseElementReference<BaseCellReference<Mesh>, Mesh, CellTag_> {
public:

  /// @brief Construct base cell reference.
  explicit BaseCellReference(Mesh& mesh, CellIndex index) noexcept :
    BaseElementReference<BaseCellReference<Mesh>, Mesh, CellTag_>(mesh, index) {
  }

  /// @brief Copy (make const) constructor. 
  BaseCellReference( // NOLINT(google-explicit-constructor)
      BaseCellReference<std::remove_const_t<Mesh>> const& other) noexcept :
    BaseElementReference<BaseCellReference<Mesh>, Mesh, CellTag_>(other) {
  }

  /// @brief Get cell volume/area/length.
  real_t Volume() const noexcept {
    return this->Mesh_->CellVolume(this->Index_);
  }

  /// @brief Get cell center position.
  vec3_t CenterPos() const noexcept {
    return this->Mesh_->CellCenterPos(this->Index_);
  }

}; // class BaseCellReference<...>

/// @brief Mesh cells random-access reference.
/// @{
using CellRef = BaseCellReference<Mesh const>;
using CellMutableRef = BaseCellReference<Mesh>;
/// @}

/// @brief Range of the @p mesh nodes.
template<class Mesh>
auto NodeRefs(Mesh& mesh) noexcept {
  return mesh.Nodes() |
    views::transform([&mesh](NodeIndex nodeIndex) {
      return BaseNodeReference<Mesh>(mesh, nodeIndex);
    });
}

/// @brief Range of the @p mesh edges.
template<class Mesh>
auto EdgeRefs(Mesh& mesh) noexcept {
  return mesh.Edges() |
    views::transform([&mesh](EdgeIndex edgeIndex) {
      return BaseEdgeReference<Mesh>(mesh, edgeIndex);
    });
}

/// @brief Range of the @p mesh nodes.
template<class Mesh>
auto FaceRefs(Mesh& mesh) noexcept {
  return mesh.Faces() |
    views::transform([&mesh](FaceIndex faceIndex) {
      return BaseFaceReference<Mesh>(mesh, faceIndex);
    });
}

/// @brief Range of the @p mesh nodes.
template<class Mesh>
auto CellRefs(Mesh& mesh) noexcept {
  return mesh.Cells() |
    views::transform([&mesh](CellIndex cellIndex) {
      return BaseCellReference<Mesh>(mesh, cellIndex);
    });
}

/// @brief Range of the @p mesh nodes with a @p nodeMark.
template<class Mesh>
auto NodeRefs(Mesh& mesh, NodeMark nodeMark) noexcept {
  return mesh.Nodes(nodeMark) |
    views::transform([&mesh](NodeIndex nodeIndex) {
      return BaseNodeReference<Mesh>(mesh, nodeIndex);
    });
}

/// @brief Range of the @p mesh edges with a @p edgeMark.
template<class Mesh>
auto EdgeRefs(Mesh& mesh, EdgeMark edgeMark) noexcept {
  return mesh.Edges(edgeMark) |
    views::transform([&mesh](EdgeIndex edgeIndex) {
      return BaseEdgeReference<Mesh>(mesh, edgeIndex);
    });
}

/// @brief Range of the @p mesh nodes with a @p faceMark.
template<class Mesh>
auto FaceRefs(Mesh& mesh, FaceMark faceMark) noexcept {
  return mesh.Faces(faceMark) |
    views::transform([&mesh](FaceIndex faceIndex) {
      return BaseFaceReference<Mesh>(mesh, faceIndex);
    });
}

/// @brief Range of the @p mesh nodes with a @p cellMark.
template<class Mesh>
auto CellRefs(Mesh& mesh, CellMark cellMark) noexcept {
  return mesh.Cells(cellMark) |
    views::transform([&mesh](CellIndex cellIndex) {
      return BaseCellReference<Mesh>(mesh, cellIndex);
    });
}

/// @brief Range of the interior @p mesh nodes.
template<class Mesh>
auto InteriorNodeRefs(Mesh& mesh) noexcept {
  return NodeRefs(mesh, NodeMark{0});
}

/// @brief Range of the interior @p mesh edges.
template<class Mesh>
auto InteriorEdgeRefs(Mesh& mesh) noexcept {
  return EdgeRefs(mesh, EdgeMark{0});
}

/// @brief Range of the interior @p mesh nodes.
template<class Mesh>
auto InteriorFaceRefs(Mesh& mesh) noexcept {
  return FaceRefs(mesh, FaceMark{0});
}

/// @brief Range of the interior @p mesh nodes.
template<class Mesh>
auto InteriorCellRefs(Mesh& mesh) noexcept {
  return CellRefs(mesh, CellMark{0});
}

#define FaceCellFunc_ \
  ([&](BaseFaceReference<Mesh> face) { \
     func(face.InnerCell(), face.OuterCell()); \
   })

/// @brief Range of the boundary @p mesh nodes.
template<class Mesh>
auto BoundaryNodeRefs(Mesh& mesh) noexcept {
  return NodeRefs(mesh) | views::drop(mesh.NumNodes({}));
}

/// @brief Range of the boundary @p mesh edges.
template<class Mesh>
auto BoundaryEdgeRefs(Mesh& mesh) noexcept {
  return EdgeRefs(mesh) | views::drop(mesh.NumEdges({}));
}

/// @brief Range of the boundary @p mesh nodes.
template<class Mesh>
auto BoundaryFaceRefs(Mesh& mesh) noexcept {
  return FaceRefs(mesh) | views::drop(mesh.NumFaces({}));
}

template<class Mesh>
auto BoundaryFaceCellRefs(Mesh& mesh) noexcept; /*{
  return FaceRefs(mesh) | views::drop(mesh.NumFaces({})) |
    views::transform([](BaseFaceReference<Mesh> face) {
      func(face.InnerCell(), face.OuterCell());
    });
}*/

/// @brief Range of the boundary @p mesh nodes.
template<class Mesh>
auto BoundaryCellRefs(Mesh& mesh) noexcept {
  return CellRefs(mesh) | views::drop(mesh.NumCells({}));
}

/// @brief Iterate through all boundary @p mesh faces.
template<class Mesh, class Func>
FEATHERS_DEPRECATED void ForEachBoundaryFaceCells(Mesh& mesh, Func func) noexcept {
  ForEach(BoundaryFaceRefs(mesh), FaceCellFunc_);
}

#undef FaceCellFunc_

} // namespace feathers
