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

template<class Mesh>
class BaseNodeIterator;
template<class Mesh>
class BaseEdgeIterator;
template<class Mesh>
class BaseFaceIterator;
template<class Mesh>
class BaseCellIterator;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base element iterator.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Iterator, class Mesh, class Tag>
class BaseElementIterator {
private:
  Mesh* Mesh_;
  Index<uint_t, Tag> Index_;

protected:

  BaseElementIterator(Mesh& mesh, Index<uint_t, Tag> index) noexcept :
      Mesh_(&mesh), Index_(index) {
    StormAssert(Index_ != npos);
  }

  template<class OtherIterator, class OtherMesh>
  BaseElementIterator(
      BaseElementIterator<OtherIterator, OtherMesh, Tag> const& other) noexcept :
      Mesh_(&other.GetMesh()), Index_(other.GetIndex()) {
    StormAssert(Index_ != npos);
  }

public:

  /// @brief Get associated mesh. 
  Mesh& GetMesh() const noexcept {
    return *Mesh_;
  }

  /// @brief Get associated element index. 
  Index<uint_t, Tag> GetIndex() const noexcept {
    return Index_;
  }

  /// @brief Cast to index operator. 
  /*explicit*/ operator Index<uint_t, Tag>() const noexcept {
    return Index_;
  }
  /*explicit*/ operator uint_t() const noexcept {
    return static_cast<uint_t>(Index_);
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

  /// @brief Identity dereference operator (used for standard algorithms). 
  /// @{
  Iterator& operator*() noexcept {
    return static_cast<Iterator&>(*this);
  }
  Iterator const& operator*() const noexcept {
    return static_cast<Iterator const&>(*this);
  }
  /// @}

  /// @brief Get mark. 
  Index<uint_t, MarkTag_<Tag>> Mark() const noexcept {
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
    //return std::size(Mesh_->AdjacentNodes(Index_));
    return std::size(Nodes());
  }

  /// @brief Number of edges in the element. 
  size_t NumEdges() const noexcept {
    return std::size(Edges());
  }

  /// @brief Number of faces in the element. 
  size_t NumFaces() const noexcept {
    return std::size(Faces());
  }

  /// @brief Number of cells in the element. 
  size_t NumCells() const noexcept {
    return std::size(Cells());
  }

  /// @brief Ranges of the adjacent nodes.
  auto Nodes() const noexcept {
    return //Mesh_->AdjacentNodes(Index_) |
      std::ranges::subrange(Mesh_->BeginAdjacentNode(Index_), Mesh_->EndAdjacentNode(Index_)) |
      std::views::transform([this](NodeIndex nodeIndex) {
        return BaseNodeIterator<Mesh>(*Mesh_, nodeIndex);
      });
  }

  /// @brief Ranges of the adjacent edges.
  auto Edges() const noexcept {
    return //Mesh_->AdjacentEdges(Index_) |
      std::ranges::subrange(Mesh_->BeginAdjacentEdge(Index_), Mesh_->EndAdjacentEdge(Index_)) |
      std::views::transform([this](EdgeIndex edgeIndex) {
        return BaseEdgeIterator<Mesh>(*Mesh_, edgeIndex);
      });
  }

  /// @brief Ranges of the adjacent faces.
  auto Faces() const noexcept {
    return //Mesh_->AdjacentFaces(Index_) |
      std::ranges::subrange(Mesh_->BeginAdjacentFace(Index_), Mesh_->EndAdjacentFace(Index_)) |
      std::views::transform([this](FaceIndex faceIndex) {
        return BaseFaceIterator<Mesh>(*Mesh_, faceIndex);
      });
  }

  /// @brief Ranges of the adjacent cells.
  auto Cells() const noexcept {
    return //Mesh_->AdjacentCells(Index_) |
      std::ranges::subrange(Mesh_->BeginAdjacentCell(Index_), Mesh_->EndAdjacentCell(Index_)) |
      std::views::transform([this](CellIndex cellIndex) {
        return BaseCellIterator<Mesh>(*Mesh_, cellIndex);
      });
  }

  /// @brief Ranges of the adjacent elements.
  template<class OtherTag>
  auto Adjacent() const noexcept {
    if constexpr (std::is_same_v<OtherTag, NodeTag_>) {
      return Nodes();
    } else if constexpr (std::is_same_v<OtherTag, EdgeTag_>) {
      return Edges();
    } else if constexpr (std::is_same_v<OtherTag, FaceTag_>) {
      return Faces();
    } else if constexpr (std::is_same_v<OtherTag, CellTag_>) {
      return Cells();
    }
  }

  /// @brief Pointer to the beginning of the adjacent elements. 
  /// @{
  [[deprecated("")]] auto Begin(NodeTag_ const&) const noexcept {
    return Mesh_->BeginAdjacentNode(Index_);
  }
  [[deprecated("")]] auto Begin(EdgeTag_ const&) const noexcept {
    return Mesh_->BeginAdjacentEdge(Index_);
  }
  [[deprecated("")]] auto Begin(FaceTag_ const&) const noexcept {
    return Mesh_->BeginAdjacentFace(Index_);
  }
  [[deprecated("")]] auto Begin(CellTag_ const&) const noexcept {
    return Mesh_->BeginAdjacentCell(Index_);
  }
  /// @}

  /// @brief Pointer to the end of the adjacent elements. 
  /// @{
  [[deprecated("")]] auto End(NodeTag_ const&) const {
    return Mesh_->EndAdjacentNode(Index_);
  }
  [[deprecated("")]] auto End(EdgeTag_ const&) const {
    return Mesh_->EndAdjacentEdge(Index_);
  }
  [[deprecated("")]] auto End(FaceTag_ const&) const {
    return Mesh_->EndAdjacentFace(Index_);
  }
  [[deprecated("")]] auto End(CellTag_ const&) const {
    return Mesh_->EndAdjacentCell(Index_);
  }
  /// @}

  /// @brief Get adjacent node. 
  auto Node(uint_t nodeLocal) const noexcept {
    StormAssert(nodeLocal < NumNodes());
    return BaseNodeIterator<Mesh>(
      *Mesh_, Mesh_->BeginAdjacentNode(Index_)[nodeLocal]);
  }

  /// @brief Get adjacent edge. 
  auto Edge(uint_t edgeLocal) const noexcept {
    StormAssert(edgeLocal < NumEdges());
    return BaseEdgeIterator<Mesh>(
      *Mesh_, Mesh_->BeginAdjacentEdge(Index_)[edgeLocal]);
  }

  /// @brief Get adjacent face. 
  auto Face(uint_t faceLocal) const noexcept {
    StormAssert(faceLocal < NumFaces());
    return BaseFaceIterator<Mesh>(
      *Mesh_, Mesh_->BeginAdjacentFace(Index_)[faceLocal]);
  }

  /// @brief Get adjacent cell. 
  auto Cell(uint_t cellLocal) const noexcept {
    StormAssert(cellLocal < NumCells());
    return BaseCellIterator<Mesh>(
      *Mesh_, Mesh_->BeginAdjacentCell(Index_)[cellLocal]);
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
    std::ranges::for_each(Faces(), [&](BaseFaceIterator<Mesh> face) {
      func(face.InnerCell(), face.OuterCell());
    });
  }
  /// @}

  /// @brief Iterate through all connected cells. 
  template<class Func>
  void ForEachCell(Func func) const noexcept {
    std::ranges::for_each(Cells(), func);
  }

}; // class BaseElementIterator

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base node iterator.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Mesh>
class BaseNodeIterator final :
  public BaseElementIterator<BaseNodeIterator<Mesh>, Mesh, NodeTag_> {
public:

  /// @brief Construct base node iterator. 
  explicit BaseNodeIterator(Mesh& mesh, NodeIndex index = {}) noexcept :
    BaseElementIterator<BaseNodeIterator<Mesh>, Mesh, NodeTag_>(mesh, index) {
  }

  /// @brief Copy (make const) constructor. 
  BaseNodeIterator( // NOLINT(google-explicit-constructor)
      BaseNodeIterator<std::remove_const_t<Mesh>> const& other) noexcept :
    BaseElementIterator<BaseNodeIterator<Mesh>, Mesh, NodeTag_>(other) {
  }

  /// @brief Get node position. 
  vec3_t Pos() const noexcept {
    return this->GetMesh().NodePos(this->GetIndex());
  }

  /// @brief Set node position. 
  void SetPos(vec3_t const& pos) const noexcept requires(!std::is_const_v<Mesh>) {
    this->GetMesh().SetNodePos(this->GetIndex(), pos);
  }

}; // class BaseNodeIterator<...>

/// @brief Mesh nodes random-access iterator. 
/// @{
using NodeIter = BaseNodeIterator<Mesh const>;
using NodeMutableIter = BaseNodeIterator<Mesh>;
/// @}

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base edge iterator.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Mesh>
class BaseEdgeIterator final :
  public BaseElementIterator<BaseEdgeIterator<Mesh>, Mesh, EdgeTag_> {
public:

  /// @brief Construct base edge iterator. 
  explicit BaseEdgeIterator(Mesh& mesh, EdgeIndex index = {}) noexcept :
    BaseElementIterator<BaseEdgeIterator<Mesh>, Mesh, EdgeTag_>(mesh, index) {
  }

  /// @brief Copy (make const) constructor. 
  BaseEdgeIterator( // NOLINT(google-explicit-constructor)
      BaseEdgeIterator<std::remove_const_t<Mesh>> const& other) noexcept :
    BaseElementIterator<BaseEdgeIterator<Mesh>, Mesh, EdgeTag_>(other) {
  }

  /// @brief Get edge length. 
  real_t Len() const noexcept {
    return this->GetMesh().EdgeLen(this->GetIndex());
  }

  /// @brief Get edge direction. 
  vec3_t Dir() const noexcept {
    return this->GetMesh().EdgeDir(this->GetIndex());
  }

  /// @brief Set edge length. 
  void SetLen(real_t len) const noexcept requires(!std::is_const_v<Mesh>) {
    this->GetMesh().SetEdgeLen(this->GetIndex(), len);
  }

  /// @brief Set edge direction. 
  void SetDir(vec3_t const& dir) const noexcept requires(!std::is_const_v<Mesh>) {
    this->GetMesh().SetEdgeDir(this->GetIndex(), dir);
  }

}; // class BaseEdgeIterator<...>

/// @brief Mesh edges random-access iterator. 
/// @{
using EdgeIter = BaseEdgeIterator<Mesh const>;
using EdgeMutableIter = BaseEdgeIterator<Mesh>;
/// @}

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base face iterator.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Mesh>
class BaseFaceIterator final :
  public BaseElementIterator<BaseFaceIterator<Mesh>, Mesh, FaceTag_> {
public:

  /// @brief Construct base face iterator. 
  explicit BaseFaceIterator(Mesh& mesh, FaceIndex index = {}) noexcept :
    BaseElementIterator<BaseFaceIterator<Mesh>, Mesh, FaceTag_>(mesh, index) {
  }

  /// @brief Copy (make const) constructor. 
  BaseFaceIterator( // NOLINT(google-explicit-constructor)
      BaseFaceIterator<std::remove_const_t<Mesh>> const& other) noexcept :
    BaseElementIterator<BaseFaceIterator<Mesh>, Mesh, FaceTag_>(other) {
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
    return this->GetMesh().FaceArea(this->GetIndex());
  }

  /// @brief Get face normal. 
  vec3_t Normal() const noexcept {
    return this->GetMesh().FaceNormal(this->GetIndex());
  }

  /// @brief Get face barycenter. 
  vec3_t CenterPos() const noexcept {
    return this->GetMesh().FaceCenterPos(this->GetIndex());
  }

  /// @brief Set face area/length. 
  void SetArea(real_t area) const noexcept requires(!std::is_const_v<Mesh>) {
    this->GetMesh().SetFaceArea(this->GetIndex(), area);
  }

  /// @brief Set face normal. 
  void SetNormal(vec3_t const& normal) const noexcept requires(!std::is_const_v<Mesh>) {
    this->GetMesh().SetFaceNormal(this->GetIndex(), normal);
  }

  /// @brief Set face barycenter. 
  void SetCenterPos(vec3_t const& centerPos) const noexcept requires(!std::is_const_v<Mesh>) {
    this->GetMesh().SetFaceCenterPos(this->GetIndex(), centerPos);
  }

}; // class BaseFaceIterator<...>

/// @brief Mesh faces random-access iterator. 
/// @{
using FaceIter = BaseFaceIterator<Mesh const>;
using FaceMutableIter = BaseFaceIterator<Mesh>;
/// @}

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base cell iterator.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Mesh>
class BaseCellIterator final :
  public BaseElementIterator<BaseCellIterator<Mesh>, Mesh, CellTag_> {
public:

  /// @brief Construct base cell iterator. 
  explicit BaseCellIterator(Mesh& mesh, CellIndex index = {}) noexcept :
    BaseElementIterator<BaseCellIterator<Mesh>, Mesh, CellTag_>(mesh, index) {
  }

  /// @brief Copy (make const) constructor. 
  BaseCellIterator( // NOLINT(google-explicit-constructor)
      BaseCellIterator<std::remove_const_t<Mesh>> const& other) noexcept :
    BaseElementIterator<BaseCellIterator<Mesh>, Mesh, CellTag_>(other) {
  }

  /// @brief Get cell volume/area/length.
  real_t Volume() const noexcept {
    return this->GetMesh().CellVolume(this->GetIndex());
  }

  /// @brief Get cell barycenter. 
  vec3_t CenterPos() const noexcept {
    return this->GetMesh().CellCenterPos(this->GetIndex());
  }

  /// @brief Set cell volume/area/length. 
  void SetVolume(real_t volume) const noexcept requires(!std::is_const_v<Mesh>) {
    this->GetMesh().SetCellVolume(this->GetIndex(), volume);
  }

  /// @brief Set cell barycenter. 
  void SetCenterPos(vec3_t const& centerPos) const noexcept requires(!std::is_const_v<Mesh>) {
    this->GetMesh().SetCellCenterPos(this->GetIndex(), centerPos);
  }

}; // class BaseCellIterator<...>

/// @brief Mesh cells random-access iterator. 
/// @{
using CellIter = BaseCellIterator<Mesh const>;
using CellMutableIter = BaseCellIterator<Mesh>;
/// @}

#define FaceCellFunc_ \
  ([&](BaseFaceIterator<Mesh> face) { \
     func(face.InnerCell(), face.OuterCell()); \
   })

/// @brief Iterator pointing to the first node with a given mark or the first mark. 
/// @{
template<class Mesh>
auto BeginNode(Mesh& mesh) noexcept {
  return BaseNodeIterator<Mesh>(mesh);
}
template<class Mesh>
auto BeginNode(Mesh& mesh, NodeMark nodeMark) noexcept {
  return BaseNodeIterator<Mesh>(mesh, mesh.BeginNode(nodeMark));
}
/// @}

/// @brief Iterator pointing to the first edge with a given mark or the first mark. 
/// @{
template<class Mesh>
auto BeginEdge(Mesh& mesh) noexcept {
  return BaseEdgeIterator<Mesh>(mesh);
}
template<class Mesh>
auto BeginEdge(Mesh& mesh, EdgeMark edgeMark) noexcept {
  return BaseEdgeIterator<Mesh>(mesh, mesh.BeginEdge(edgeMark));
}
/// @}

/// @brief Iterator pointing to the first face with a given mark or the first mark. 
/// @{
template<class Mesh>
auto BeginFace(Mesh& mesh) noexcept {
  return BaseFaceIterator<Mesh>(mesh);
}
template<class Mesh>
auto BeginFace(Mesh& mesh, FaceMark faceMark) noexcept {
  return BaseFaceIterator<Mesh>(mesh, mesh.BeginFace(faceMark));
}
/// @}

/// @brief Iterator pointing to the first cell with a given mark or the first mark. 
/// @{
template<class Mesh>
auto BeginCell(Mesh& mesh) noexcept {
  return BaseCellIterator<Mesh>(mesh);
}
template<class Mesh>
auto BeginCell(Mesh& mesh, CellMark cellMark) noexcept {
  return BaseCellIterator<Mesh>(mesh, mesh.BeginCell(cellMark));
}
/// @}

/// @brief Iterator pointing to a node after the last node with a given mark or the last mark. 
/// @{
template<class Mesh>
auto EndNode(Mesh& mesh) noexcept {
  return BaseNodeIterator<Mesh>(mesh, NodeIndex(mesh.NumNodes()));
}
template<class Mesh>
auto EndNode(Mesh& mesh, NodeMark nodeMark) noexcept {
  return BaseNodeIterator<Mesh>(mesh, mesh.EndNode(nodeMark));
}
/// @}

/// @brief Iterator pointing to an edge after the last edge with a given mark or the last mark. 
/// @{
template<class Mesh>
auto EndEdge(Mesh& mesh) noexcept {
  return BaseEdgeIterator<Mesh>(mesh, EdgeIndex(mesh.NumEdges()));
}
template<class Mesh>
auto EndEdge(Mesh& mesh, EdgeMark edgeMark) noexcept {
  return BaseEdgeIterator<Mesh>(mesh, mesh.EndEdge(edgeMark));
}
/// @}

/// @brief Iterator pointing to a face after the last face with a given mark or the last mark. 
/// @{
template<class Mesh>
auto EndFace(Mesh& mesh) noexcept {
  return BaseFaceIterator<Mesh>(mesh, FaceIndex(mesh.NumFaces()));
}
template<class Mesh>
auto EndFace(Mesh& mesh, FaceMark faceMark) noexcept {
  return BaseFaceIterator<Mesh>(mesh, mesh.EndFace(faceMark));
}
/// @}

/// @brief Iterator pointing to a cell after the last cell with a given mark or the last mark. 
/// @{
template<class Mesh>
auto EndCell(Mesh& mesh) noexcept {
  return BaseCellIterator<Mesh>(mesh, CellIndex(mesh.NumCells()));
}
template<class Mesh>
auto EndCell(Mesh& mesh, CellMark cellMark) noexcept {
  return BaseCellIterator<Mesh>(mesh, mesh.EndCell(cellMark));
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
