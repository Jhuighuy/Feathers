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

#include <ranges>

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

  BaseElementReference(Mesh& mesh, Index<Tag> index) noexcept :
      Mesh_(&mesh), Index_(index) {
    StormAssert(Index_ != npos);
  }

  template<class OtherIterator, class OtherMesh>
  BaseElementReference(
    BaseElementReference<OtherIterator, OtherMesh, Tag> const& other) noexcept :
      Mesh_(other.Mesh_), Index_(other.Index_) {
    StormAssert(Index_ != npos);
  }

public:

  /// @brief Identity dereference operator (used for standard algorithms).
  /// @{
  Iterator& operator*() noexcept {
    return static_cast<Iterator&>(*this);
  }
  Iterator const& operator*() const noexcept {
    return static_cast<Iterator const&>(*this);
  }
  /// @}

  /// @brief Cast to index operator. 
  /*explicit*/ operator Index<Tag>() const noexcept {
    return Index_;
  }
  /*explicit*/ operator size_t() const noexcept {
    return static_cast<size_t>(Index_);
  }

  /// @brief Difference operator. 
  ptrdiff_t operator-(Iterator const& other) const noexcept {
    StormAssert(Mesh_ == other.Mesh_);
    return Index_ - other.Index_;
  }

  /// @brief Equality operator. 
  bool operator==(Iterator const& other) const noexcept {
    StormAssert(Mesh_ == other.Mesh_);
    return Index_ == other.Index_;
  }

  /// @brief Inequality operator. 
  bool operator!=(Iterator const& other) const noexcept {
    StormAssert(Mesh_ == other.Mesh_);
    return Index_ != other.Index_;
  }

  /// @brief Less than operator. 
  bool operator<(Iterator const& other) const noexcept {
    StormAssert(Mesh_ == other.Mesh_);
    return Index_ < other.Index_;
  }

  /// @brief Less than or equal operator. 
  bool operator<=(Iterator const& other) const noexcept {
    StormAssert(Mesh_ == other.Mesh_);
    return Index_ <= other.Index_;
  }

  /// @brief Greater than operator. 
  bool operator>(Iterator const& other) const noexcept {
    StormAssert(Mesh_ == other.Mesh_);
    return Index_ > other.Index_;
  }

  /// @brief Greater than or equal operator. 
  bool operator>=(Iterator const& other) const noexcept {
    StormAssert(Mesh_ == other.Mesh_);
    return Index_ >= other.Index_;
  }

  /// @brief Increment operator. 
  /// @{
  Iterator& operator++() noexcept {
    ++Index_;
    return static_cast<Iterator&>(*this);
  }
  Iterator const operator++(int) noexcept {
    Iterator iter(*this);
    return ++*this, iter;
  }
  /// @}

  /// @brief Decrement operator. 
  /// @{
  Iterator& operator--() noexcept {
    --Index_;
    return static_cast<Iterator&>(*this);
  }
  Iterator const operator--(int) noexcept {
    Iterator iter(static_cast<Iterator&>(*this));
    return --*this, iter;
  }
  /// @}

  /// @brief Addition operator. 
  /// @{
  Iterator& operator+=(ptrdiff_t offset) noexcept {
    Index_ += offset;
    return static_cast<Iterator&>(*this);
  }
  Iterator operator+(ptrdiff_t offset) const noexcept {
    return Iterator(static_cast<Iterator const&>(*this)) += offset;
  }
  /// @}

  /// @brief Subtraction operator. 
  /// @{
  Iterator& operator-=(ptrdiff_t offset) noexcept {
    Index_ -= offset;
    return static_cast<Iterator&>(*this);
  }
  Iterator operator-(ptrdiff_t offset) const noexcept {
    return Iterator(static_cast<Iterator const&>(*this)) -= offset;
  }
  /// @}

  /// @brief Get mark. 
  Index<MarkTag_<Tag>> Mark() const noexcept {
    return Mesh_->Mark(Index_);
  }

  /// @brief Get shape. 
  eShape Shape() const noexcept {
    return Mesh_->Shape(Index_);
  }

  /// @brief Get element object. 
  std::unique_ptr<const iElement> get_element_object() const {
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
      std::views::transform([&mesh = *Mesh_](NodeIndex nodeIndex) {
        return BaseNodeReference<Mesh>(mesh, nodeIndex);
      });
  }

  /// @brief Ranges of the adjacent edges.
  auto Edges() const noexcept {
    return Mesh_->AdjacentEdges(Index_) |
      std::views::transform([&mesh = *Mesh_](EdgeIndex edgeIndex) {
        return BaseEdgeReference<Mesh>(mesh, edgeIndex);
      });
  }

  /// @brief Ranges of the adjacent faces.
  auto Faces() const noexcept {
    return Mesh_->AdjacentFaces(Index_) |
      std::views::transform([&mesh = *Mesh_](FaceIndex faceIndex) {
        return BaseFaceReference<Mesh>(mesh, faceIndex);
      });
  }

  /// @brief Ranges of the adjacent cells.
  auto Cells() const noexcept {
    return Mesh_->AdjacentCells(Index_) |
      std::views::transform([&mesh = *Mesh_](CellIndex cellIndex) {
        return BaseCellReference<Mesh>(mesh, cellIndex);
      });
  }

  /// @brief Get adjacent node. 
  auto Node(size_t nodeLocal) const noexcept {
    StormAssert(nodeLocal < NumNodes());
    return BaseNodeReference<Mesh>(
      *Mesh_, Mesh_->AdjacentNodes(Index_)[nodeLocal]);
  }

  /// @brief Get adjacent edge. 
  auto Edge(size_t edgeLocal) const noexcept {
    StormAssert(edgeLocal < NumEdges());
    return BaseEdgeReference<Mesh>(
      *Mesh_, Mesh_->AdjacentEdges(Index_)[edgeLocal]);
  }

  /// @brief Get adjacent face. 
  auto Face(size_t faceLocal) const noexcept {
    StormAssert(faceLocal < NumFaces());
    return BaseFaceReference<Mesh>(
      *Mesh_, Mesh_->AdjacentFaces(Index_)[faceLocal]);
  }

  /// @brief Get adjacent cell. 
  auto Cell(size_t cellLocal) const noexcept {
    StormAssert(cellLocal < NumCells());
    return BaseCellReference<Mesh>(
      *Mesh_, Mesh_->AdjacentCells(Index_)[cellLocal]);
  }

  /// @brief Iterate through all connected nodes. 
  template<class Func>
  void ForEachNode(Func func) const noexcept {
    std::ranges::for_each(Nodes(), func);
  }

  /// @brief Iterate through all connected edges. 
  template<class Func>
  void ForEachEdge(Func func) const noexcept {
    std::ranges::for_each(Edges(), func);
  }

  /// @brief Iterate through all connected faces. 
  /// @{
  template<class Func>
  void ForEachFace(Func func) const noexcept {
    std::ranges::for_each(Faces(), func);
  }
  template<class Func>
  void ForEachFaceCells(Func func) const noexcept {
    std::ranges::for_each(Faces(), [&](BaseFaceReference<Mesh> face) {
      func(face.InnerCell(), face.OuterCell());
    });
  }
  /// @}

  /// @brief Iterate through all connected cells. 
  template<class Func>
  void ForEachCell(Func func) const noexcept {
    std::ranges::for_each(Cells(), func);
  }

}; // class BaseElementReference

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base node reference.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Mesh>
class BaseNodeReference final :
  public BaseElementReference<BaseNodeReference<Mesh>, Mesh, NodeTag_> {
public:

  /// @brief Construct base node reference.
  explicit BaseNodeReference(Mesh& mesh, NodeIndex index = {}) noexcept :
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

  /// @brief Set node position. 
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
  explicit BaseEdgeReference(Mesh& mesh, EdgeIndex index = {}) noexcept :
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

  /// @brief Set edge length. 
  void SetLen(real_t len) const noexcept requires(!std::is_const_v<Mesh>) {
    this->Mesh_->SetEdgeLen(this->Index_, len);
  }

  /// @brief Set edge direction. 
  void SetDir(vec3_t const& dir) const noexcept requires(!std::is_const_v<Mesh>) {
    this->Mesh_->SetEdgeDir(this->Index_, dir);
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
  explicit BaseFaceReference(Mesh& mesh, FaceIndex index = {}) noexcept :
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

  /// @brief Get face barycenter. 
  vec3_t CenterPos() const noexcept {
    return this->Mesh_->FaceCenterPos(this->Index_);
  }

  /// @brief Set face area/length. 
  void SetArea(real_t area) const noexcept requires(!std::is_const_v<Mesh>) {
    this->Mesh_->SetFaceArea(this->Index_, area);
  }

  /// @brief Set face normal. 
  void SetNormal(vec3_t const& normal) const noexcept requires(!std::is_const_v<Mesh>) {
    this->Mesh_->SetFaceNormal(this->Index_, normal);
  }

  /// @brief Set face barycenter. 
  void SetCenterPos(vec3_t const& centerPos) const noexcept requires(!std::is_const_v<Mesh>) {
    this->Mesh_->SetFaceCenterPos(this->Index_, centerPos);
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
  explicit BaseCellReference(Mesh& mesh, CellIndex index = {}) noexcept :
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

  /// @brief Get cell barycenter. 
  vec3_t CenterPos() const noexcept {
    return this->Mesh_->CellCenterPos(this->Index_);
  }

  /// @brief Set cell volume/area/length. 
  void SetVolume(real_t volume) const noexcept requires(!std::is_const_v<Mesh>) {
    this->Mesh_->SetCellVolume(this->Index_, volume);
  }

  /// @brief Set cell barycenter. 
  void SetCenterPos(vec3_t const& centerPos) const noexcept requires(!std::is_const_v<Mesh>) {
    this->Mesh_->SetCellCenterPos(this->Index_, centerPos);
  }

}; // class BaseCellReference<...>

/// @brief Mesh cells random-access reference.
/// @{
using CellRef = BaseCellReference<Mesh const>;
using CellMutableRef = BaseCellReference<Mesh>;
/// @}

template<class Mesh>
auto NodeIters(Mesh& mesh) noexcept {
  return mesh.Nodes()  |
    std::views::transform([&mesh](NodeIndex nodeIndex) {
      return BaseNodeReference<Mesh>(mesh, nodeIndex);
    });
}

template<class Mesh>
auto NodeIters(Mesh& mesh, NodeMark nodeMark) noexcept {
  return mesh.Nodes(nodeMark)  |
    std::views::transform([&mesh](NodeIndex nodeIndex) {
      return BaseNodeReference<Mesh>(mesh, nodeIndex);
    });
}

#define FaceCellFunc_ \
  ([&](BaseFaceReference<Mesh> face) { \
     func(face.InnerCell(), face.OuterCell()); \
   })

/// @brief Iterator pointing to the first node with a given mark or the first mark. 
/// @{
template<class Mesh>
auto BeginNode(Mesh& mesh) noexcept {
  return BaseNodeReference<Mesh>(mesh);
}
template<class Mesh>
auto BeginNode(Mesh& mesh, NodeMark nodeMark) noexcept {
  return BaseNodeReference<Mesh>(mesh, mesh.BeginNode(nodeMark));
}
/// @}

/// @brief Iterator pointing to the first edge with a given mark or the first mark. 
/// @{
template<class Mesh>
auto BeginEdge(Mesh& mesh) noexcept {
  return BaseEdgeReference<Mesh>(mesh);
}
template<class Mesh>
auto BeginEdge(Mesh& mesh, EdgeMark edgeMark) noexcept {
  return BaseEdgeReference<Mesh>(mesh, mesh.BeginEdge(edgeMark));
}
/// @}

/// @brief Iterator pointing to the first face with a given mark or the first mark. 
/// @{
template<class Mesh>
auto BeginFace(Mesh& mesh) noexcept {
  return BaseFaceReference<Mesh>(mesh);
}
template<class Mesh>
auto BeginFace(Mesh& mesh, FaceMark faceMark) noexcept {
  return BaseFaceReference<Mesh>(mesh, mesh.BeginFace(faceMark));
}
/// @}

/// @brief Iterator pointing to the first cell with a given mark or the first mark. 
/// @{
template<class Mesh>
auto BeginCell(Mesh& mesh) noexcept {
  return BaseCellReference<Mesh>(mesh);
}
template<class Mesh>
auto BeginCell(Mesh& mesh, CellMark cellMark) noexcept {
  return BaseCellReference<Mesh>(mesh, mesh.BeginCell(cellMark));
}
/// @}

/// @brief Iterator pointing to a node after the last node with a given mark or the last mark. 
/// @{
template<class Mesh>
auto EndNode(Mesh& mesh) noexcept {
  return BaseNodeReference<Mesh>(mesh, NodeIndex(mesh.NumNodes()));
}
template<class Mesh>
auto EndNode(Mesh& mesh, NodeMark nodeMark) noexcept {
  return BaseNodeReference<Mesh>(mesh, mesh.EndNode(nodeMark));
}
/// @}

/// @brief Iterator pointing to an edge after the last edge with a given mark or the last mark. 
/// @{
template<class Mesh>
auto EndEdge(Mesh& mesh) noexcept {
  return BaseEdgeReference<Mesh>(mesh, EdgeIndex(mesh.NumEdges()));
}
template<class Mesh>
auto EndEdge(Mesh& mesh, EdgeMark edgeMark) noexcept {
  return BaseEdgeReference<Mesh>(mesh, mesh.EndEdge(edgeMark));
}
/// @}

/// @brief Iterator pointing to a face after the last face with a given mark or the last mark. 
/// @{
template<class Mesh>
auto EndFace(Mesh& mesh) noexcept {
  return BaseFaceReference<Mesh>(mesh, FaceIndex(mesh.NumFaces()));
}
template<class Mesh>
auto EndFace(Mesh& mesh, FaceMark faceMark) noexcept {
  return BaseFaceReference<Mesh>(mesh, mesh.EndFace(faceMark));
}
/// @}

/// @brief Iterator pointing to a cell after the last cell with a given mark or the last mark. 
/// @{
template<class Mesh>
auto EndCell(Mesh& mesh) noexcept {
  return BaseCellReference<Mesh>(mesh, CellIndex(mesh.NumCells()));
}
template<class Mesh>
auto EndCell(Mesh& mesh, CellMark cellMark) noexcept {
  return BaseCellReference<Mesh>(mesh, mesh.EndCell(cellMark));
}
/// @}

/// @brief Iterate through all nodes with a given mark or all marks. 
/// @{
template<class Mesh, class Func>
void ForEachNode(Mesh& mesh, Func func) noexcept {
  for_range(BeginNode(mesh), EndNode(mesh), func);
}
template<class Mesh, class Func>
void ForEachNode(Mesh& mesh, NodeMark nodeMark, Func func) noexcept {
  for_range(BeginNode(mesh, nodeMark), EndNode(mesh, nodeMark), func);
}
/// @}

/// @brief Iterate through all edges with a given mark or all marks. 
/// @{
template<class Mesh, class Func>
void ForEachEdge(Mesh& mesh, Func func) noexcept {
  for_range(BeginEdge(mesh), EndEdge(mesh), func);
}
template<class Mesh, class Func>
void ForEachEdge(Mesh& mesh, EdgeMark edgeMark, Func func) noexcept {
  for_range(BeginEdge(mesh, edgeMark), EndEdge(mesh, edgeMark), func);
}
/// @}

/// @brief Iterate through all faces with a given mark or all marks. 
/// @{
template<class Mesh, class Func>
void ForEachFace(Mesh& mesh, Func func) noexcept {
  for_range(BeginFace(mesh), EndFace(mesh), func);
}
template<class Mesh, class Func>
void ForEachFace(Mesh& mesh, FaceMark faceMark, Func func) noexcept {
  for_range(BeginFace(mesh, faceMark), EndFace(mesh, faceMark), func);
}
template<class Mesh, class Func>
void ForEachFaceCells(Mesh& mesh, Func func) noexcept {
  for_range(BeginFace(mesh), EndFace(mesh), FaceCellFunc_);
}
template<class Mesh, class Func>
void ForEachFaceCells(Mesh& mesh, FaceMark faceMark, Func func) noexcept {
  for_range(BeginFace(mesh, faceMark), EndFace(mesh, faceMark), FaceCellFunc_);
}
/// @}

/// @brief Iterate through all cells with a given mark or all marks. 
/// @{
template<class Mesh, class Func>
void ForEachCell(Mesh& mesh, Func func) noexcept {
  for_range(BeginCell(mesh), EndCell(mesh), func);
}
template<class Mesh, class Func>
void ForEachCell(Mesh& mesh, CellMark cellMark, Func func) noexcept {
  for_range(BeginCell(mesh, cellMark), EndCell(mesh, cellMark), func);
}
/// @}

/// @brief Iterator pointing to the first interior node.
template<class Mesh>
auto BeginInteriorNode(Mesh& mesh) noexcept {
  return BeginNode(mesh, {});
}

/// @brief Iterator pointing to the first interior edge. 
template<class Mesh>
auto BeginInteriorEdge(Mesh& mesh) noexcept {
  return BeginEdge(mesh, {});
}

/// @brief Iterator pointing to the first interior face. 
template<class Mesh>
auto BeginInteriorFace(Mesh& mesh) noexcept {
  return BeginFace(mesh, {});
}

/// @brief Iterator pointing to the first interior cell. 
template<class Mesh>
auto BeginInteriorCell(Mesh& mesh) noexcept {
  return BeginCell(mesh, {});
}

/// @brief Iterator pointing to a node after the last node. 
template<class Mesh>
auto EndInteriorNode(Mesh& mesh) noexcept {
  return EndNode(mesh, {});
}

/// @brief Iterator pointing to an edge after the last node. 
template<class Mesh>
auto EndInteriorEdge(Mesh& mesh) noexcept {
  return EndEdge(mesh, {});
}

/// @brief Iterator pointing to a face after the last node. 
template<class Mesh>
auto EndInteriorFace(Mesh& mesh) noexcept {
  return EndFace(mesh, {});
}

/// @brief Iterator pointing to a cell after the last node. 
template<class Mesh>
auto EndInteriorCell(Mesh& mesh) noexcept {
  return EndCell(mesh, {});
}

/// @brief Iterate through all interior nodes. 
template<class Mesh, class Func>
void ForEachInteriorNode(Mesh& mesh, Func func) noexcept {
  for_range(BeginInteriorNode(mesh), EndInteriorNode(mesh), func);
}

/// @brief Iterate through all interior edges. 
template<class Mesh, class Func>
void ForEachInteriorEdge(Mesh& mesh, Func func) noexcept {
  for_range(BeginInteriorEdge(mesh), EndInteriorEdge(mesh), func);
}

/// @brief Iterate through all interior faces. 
/// @{
template<class Mesh, class Func>
void ForEachInteriorCace(Mesh& mesh, Func func) noexcept {
  for_range(BeginInteriorFace(mesh), EndInteriorFace(mesh), func);
}
template<class Mesh, class Func>
void ForEachInteriorFaceCells(Mesh& mesh, Func func) noexcept {
  for_range(BeginInteriorFace(mesh), EndInteriorFace(mesh), FaceCellFunc_);
}
/// @}

/// @brief Iterate through all interior cells. 
template<class Mesh, class Func>
void ForEachInteriorCell(Mesh& mesh, Func func) noexcept {
  for_range(BeginInteriorCell(mesh), EndInteriorCell(mesh), func);
}

/// @brief Iterator pointing to the boundary first node with a given mark or the first boundary mark.
template<class Mesh>
auto BeginBoundaryNode(Mesh& mesh, NodeMark nodeMark = NodeMark{1}) {
  StormAssert(nodeMark >= 1);
  return BeginNode(mesh, nodeMark);
}

/// @brief Iterator pointing to the boundary first edge with a given mark or the first boundary mark.
template<class Mesh>
auto BeginBoundaryEdge(Mesh& mesh, EdgeMark edgeMark = EdgeMark{1}) {
  StormAssert(edgeMark >= 1);
  return BeginEdge(mesh, edgeMark);
}

/// @brief Iterator pointing to the boundary first face with a given mark or the first boundary mark.
template<class Mesh>
auto BeginBoundaryFace(Mesh& mesh, FaceMark faceMark = FaceMark{1}) {
  StormAssert(faceMark >= 1);
  return BeginFace(mesh, faceMark);
}

/// @brief Iterator pointing to the boundary first cell with a given mark or the first boundary mark.
template<class Mesh>
auto BeginBoundaryCell(Mesh& mesh, CellMark cellMark = CellMark{1}) {
  StormAssert(cellMark >= 1);
  return BeginCell(mesh, cellMark);
}

/// @brief Iterator pointing to a node after the last node with a given or the last boundary mark. 
/// @{
template<class Mesh>
auto EndBoundaryNode(Mesh& mesh) noexcept {
  return EndNode(mesh);
}
template<class Mesh>
auto EndBoundaryNode(Mesh& mesh, NodeMark nodeMark) noexcept {
  StormAssert(nodeMark >= 1);
  return EndNode(mesh, nodeMark);
}
/// @}

/// @brief Iterator pointing to an edge after the last edge with a given or the last boundary mark. 
/// @{
template<class Mesh>
auto EndBoundaryEdge(Mesh& mesh) noexcept {
  return EndEdge(mesh);
}
template<class Mesh>
auto EndBoundaryEdge(Mesh& mesh, EdgeMark edgeMark) noexcept {
  StormAssert(edgeMark >= 1);
  return EndEdge(mesh, edgeMark);
}
/// @}

/// @brief Iterator pointing to a face after the last face with a given or the last boundary mark. 
/// @{
template<class Mesh>
auto EndBoundaryFace(Mesh& mesh) noexcept {
  return EndFace(mesh);
}
template<class Mesh>
auto EndBoundaryFace(Mesh& mesh, FaceMark faceMark) noexcept {
  StormAssert(faceMark >= 1);
  return EndFace(mesh, faceMark);
}
/// @}

/// @brief Iterator pointing to a cell after the last cell with a given or the last boundary mark. 
/// @{
template<class Mesh>
auto EndBoundaryCell(Mesh& mesh) noexcept {
  return EndCell(mesh);
}
template<class Mesh>
auto EndBoundaryCell(Mesh& mesh, CellMark cellMark) noexcept {
  StormAssert(cellMark >= 1);
  return EndCell(mesh, cellMark);
}
/// @}

/// @brief Iterate through all boundary nodes with a given mark or all boundary marks. 
/// @{
template<class Mesh, class Func>
void ForEachBoundaryNode(Mesh& mesh, Func func) noexcept {
  for_range(BeginBoundaryNode(mesh), EndBoundaryNode(mesh), func);
}
template<class Mesh, class Func>
void ForEachBoundaryNode(Mesh& mesh, NodeMark nodeMark, Func func) noexcept {
  for_range(BeginBoundaryNode(mesh, nodeMark), EndBoundaryNode(mesh, nodeMark), func);
}
/// @}

/// @brief Iterate through all boundary edges with a given mark or all boundary marks. 
/// @{
template<class Mesh, class Func>
void ForEachBoundaryEdge(Mesh& mesh, Func func) noexcept {
  for_range(BeginBoundaryEdge(mesh), EndBoundaryEdge(mesh), func);
}
template<class Mesh, class Func>
void ForEachBoundaryEdge(Mesh& mesh, EdgeMark edgeMark, Func func) noexcept {
  for_range(BeginBoundaryEdge(mesh, edgeMark), EndBoundaryEdge(mesh, edgeMark), func);
}
/// @}

/// @brief Iterate through all boundary faces with a given mark or all boundary marks. 
/// @{
template<class Mesh, class Func>
void ForEachBoundaryFace(Mesh& mesh, Func func) noexcept {
  for_range(BeginBoundaryFace(mesh), EndBoundaryFace(mesh), func);
}
template<class Mesh, class Func>
void ForEachBoundaryFace(Mesh& mesh, FaceMark faceMark, Func func) noexcept {
  for_range(BeginBoundaryFace(mesh, faceMark), EndBoundaryFace(mesh, faceMark), func);
}
template<class Mesh, class Func>
void ForEachBoundaryFaceCells(Mesh& mesh, Func func) noexcept {
  for_range(BeginBoundaryFace(mesh), EndBoundaryFace(mesh), FaceCellFunc_);
}
template<class Mesh, class Func>
void ForEachBoundaryFaceCells(Mesh& mesh, FaceMark faceMark, Func func) noexcept {
  for_range(BeginBoundaryFace(mesh, faceMark), EndBoundaryFace(mesh, faceMark), FaceCellFunc_);
}
/// @}

/// @brief Iterate through all boundary cells with a given mark or all boundary marks. 
/// @{
template<class Mesh, class Func>
void ForEachBoundaryCell(Mesh& mesh, Func func) noexcept {
  for_range(BeginBoundaryCell(mesh), EndBoundaryCell(mesh), func);
}
template<class Mesh, class Func>
void ForEachBoundaryCell(Mesh& mesh, CellMark cellMark, Func func) noexcept {
  for_range(BeginBoundaryCell(mesh, cellMark), EndBoundaryCell(mesh, cellMark), func);
}
/// @}

#undef FaceCellFunc_

} // namespace feathers
