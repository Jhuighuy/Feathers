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

#include <map>

#include "SkunkBase.hh"
#include "Index.hxx"
#include "libFeathersUtils/Table.hh"
#include "libFeathersUtils/Image.hh"
#include "Field.hh"
#include "Element.hh"

namespace feathers {
template<class>
inline constexpr bool AlwaysFalse = false;
template<class RowIndex, class ColumnIndex>
class CsrMapping {
public:
  std::span<ColumnIndex> operator[](RowIndex) const;
};
template<class RowIndex, class ColumnIndex>
class BiCsrMapping {
public:
  std::span<ColumnIndex> operator[](RowIndex) const;
  CsrMapping<ColumnIndex, RowIndex>& transposed() const;
};
template<class RowIndex, class ColumnIndex>
class TriCsrMapping {
public:
  std::span<ColumnIndex> operator[](RowIndex) const;
  CsrMapping<ColumnIndex, RowIndex>& transposed() const;
  CsrMapping<RowIndex, RowIndex>& symmetric() const;
};
}

#define StormThis_

#define StormDeducingThisEmulation_(name) \
  /** @copydoc name */ \
  /** @{ */ \
  template<std::same_as<void> = void, class... Args> \
  auto name(Args const&... args) noexcept { \
    using Self = std::decay_t<decltype(*this)>; \
    return name<Self>(*this, args...); \
  } \
  template<std::same_as<void> = void, class... Args> \
  auto name(Args const&... args) const noexcept { \
    using Self = std::decay_t<decltype(*this)>; \
    return const_cast<Self&>(*this).name<Self const>(*this, args...); \
  } \
  /** @} */

#define StormDeducingThisEmulationT_(name) \
  /** @copydoc name */ \
  /** @{ */ \
  template<class T, std::same_as<void> = void, class... Args> \
  auto name(Args const&... args) noexcept { \
    using Self = std::decay_t<decltype(*this)>; \
    return name<T, Self>(*this, args...); \
  } \
  template<class T, std::same_as<void> = void, class... Args> \
  auto name(Args const&... args) const noexcept { \
    using Self = std::decay_t<decltype(*this)>; \
    return const_cast<Self&>(*this).name<T, Self const>(*this, args...); \
  } \
  /** @} */

namespace feathers {

enum : size_t {
  FaceInnerCell_ = 0,
  FaceOuterCell_ = 1,
}; // enum

class NodeTag;
class EdgeTag;
class FaceTag;
class CellTag;
template<class> class MarkTag;

/// @brief Node index.
using NodeIndex = Index<NodeTag>;

/// @brief Edge index.
using EdgeIndex = Index<EdgeTag>;

/// @brief Face index.
using FaceIndex = Index<FaceTag>;

/// @brief Cell index.
using CellIndex = Index<CellTag>;

/// @brief Node mark index.
using NodeMark = Index<MarkTag<NodeTag>>;

/// @brief Edge mark index.
using EdgeMark = Index<MarkTag<EdgeTag>>;

/// @brief Face mark index.
using FaceMark = Index<MarkTag<FaceTag>>;

/// @brief Cell mark index.
using CellMark = Index<MarkTag<CellTag>>;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Hybrid unstructured multidimensional mesh.
/// @todo Do not use iterators in the mesh implementation.
/// @todo Switch to ranges.
/// @todo "Element" -> "shapeType".
/// @todo Make Dim a template parameter.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
//template<size_t SDim, size_t TDim, template<class, class, class> class Table>
class Mesh : public tObject<Mesh> {
private:
  size_t NumNodes_ = 0, NumEdges_ = 0;
  size_t NumFaces_ = 0, NumCells_ = 0;

  IndexedVector<NodeMark, NodeIndex> NodeMarks_;
  IndexedVector<EdgeMark, EdgeIndex> EdgeMarks_;
  IndexedVector<FaceMark, FaceIndex> FaceMarks_;
  IndexedVector<CellMark, CellIndex> CellMarks_;
  IndexedVector<NodeIndex, NodeMark> NodeRanges_;
  IndexedVector<EdgeIndex, EdgeMark> EdgeRanges_;
  IndexedVector<FaceIndex, FaceMark> FaceRanges_;
  IndexedVector<CellIndex, CellMark> CellRanges_;

  IndexedVector<ShapeType, EdgeIndex> EdgeShapeTypes_;
  IndexedVector<ShapeType, FaceIndex> FaceShapeTypes_;
  IndexedVector<ShapeType, CellIndex> CellShapeTypes_;
  // TODO: edge length + direction -> 4D oriented direction.
  // TODO: face area + normal -> 4D oriented area.
  IndexedVector<vec3_t, NodeIndex> NodePos_;
  IndexedVector<real_t, EdgeIndex> EdgeLens_;
  IndexedVector<vec3_t, EdgeIndex> EdgeDirs_;
  IndexedVector<real_t, FaceIndex> FaceAreas_;
  IndexedVector<vec3_t, FaceIndex> FaceNormals_;
  IndexedVector<vec3_t, FaceIndex> FaceCenterPos_;
  IndexedVector<real_t, CellIndex> CellVolumes_;
  IndexedVector<vec3_t, CellIndex> CellCenterPos_;
  real_t MaxEdgeLen_ = 0.0, MinEdgeLen_ = 0.0;
  real_t MaxFaceArea_ = 0.0, MinFaceArea_ = 0.0;
  real_t MaxCellVolume_ = 0.0, MinCellVolume_ = 0.0;

  TriCsrMapping<EdgeIndex, NodeIndex> EdgeNodes_;
  BiCsrMapping<FaceIndex, NodeIndex> FaceNodes_;
  BiCsrMapping<CellIndex, NodeIndex> CellNodes_;
  TriCsrMapping<FaceIndex, EdgeIndex> FaceEdges_;
  BiCsrMapping<CellIndex, EdgeIndex> CellEdges_;
  TriCsrMapping<CellIndex, FaceIndex> CellFaces_;
  CsrMapping<NodeIndex, NodeIndex> NodeNodes_;

#if 0
  CsrTable<EdgeIndex, NodeIndex> EdgeNodes_; CsrTable<NodeIndex, EdgeIndex> NodeEdges_; CsrTable<EdgeIndex, EdgeIndex> EdgeEdges_;
  CsrTable<FaceIndex, NodeIndex> FaceNodes_; CsrTable<NodeIndex, FaceIndex> NodeFaces_;
  CsrTable<CellIndex, NodeIndex> CellNodes_; CsrTable<NodeIndex, CellIndex> NodeCells_;
  
  CsrTable<FaceIndex, EdgeIndex> FaceEdges_; CsrTable<EdgeIndex, FaceIndex> EdgeFaces_; CsrTable<FaceIndex, FaceIndex> FaceFaces_;
  CsrTable<CellIndex, EdgeIndex> CellEdges_; CsrTable<EdgeIndex, CellIndex> EdgeCells_;
  
  CsrTable<CellIndex, FaceIndex> CellFaces_; CsrTable<FaceIndex, CellIndex> FaceCells_; CsrTable<CellIndex, CellIndex> CellCells_;

  CsrTable<NodeIndex, NodeIndex> NodeNodes_;
#endif

public:

  /// ---------------------------------------------------------------- ///
  /// @name Constructors.
  /// ---------------------------------------------------------------- ///
  /// @{

  /// @brief Initialize an empty mesh.
  Mesh() noexcept = default;

  bool ReadFromTriangle(std::string const& path);

  bool ReadFromTetgen(std::string const& path);

  bool ReadFromImage(const char* path,
                     const std::map<Pixel, size_t>& markColors,
                     Pixel fluidColor = eBlackPixel,
                     vec2_t pixelSize = vec2_t(1.0, 1.0));

  void save_vtk(const char* path,
                const std::vector<sFieldDesc>& fields) const;

  void finalize() {
    FinalizeFaces_();
    generate_boundary_cells();
    reorder_faces();
    UpdateElementsGeometry();
  }

  /// @}

  /// ---------------------------------------------------------------- ///
  /// @name Element ranges.
  /// ---------------------------------------------------------------- ///
  /// @{

  /// @brief Range of node indices.
  auto nodeIndices() const noexcept {
    return views::iota(NodeIndex{0}, NodeIndex{NumNodes_});
  }

  /// @brief Range of edge indices.
  auto edgeIndices() const noexcept {
    return views::iota(EdgeIndex{0}, EdgeIndex{NumEdges_});
  }

  /// @brief Range of face indices.
  auto faceIndices() const noexcept {
    return views::iota(FaceIndex{0}, FaceIndex{NumFaces_});
  }

  /// @brief Range of cell indices.
  auto cellIndices() const noexcept {
    return views::iota(CellIndex{0}, CellIndex{NumCells_});
  }

  /// @}

  /// ---------------------------------------------------------------- ///
  /// @name Marks.
  /// ---------------------------------------------------------------- ///
  /// @{

  /// @brief Number of node marks.
  size_t numNodeMarks() const noexcept {
    StormAssert(!NodeRanges_.empty());
    return NodeRanges_.size() - 1;
  }

  /// @brief Number of edge marks.
  size_t numEdgeMarks() const noexcept {
    StormAssert(!EdgeRanges_.empty());
    return EdgeRanges_.size() - 1;
  }

  /// @brief Number of face marks.
  size_t numFaceMarks() const noexcept {
    StormAssert(!FaceRanges_.empty());
    return FaceRanges_.size() - 1;
  }

  /// @brief Number of cell marks.
  size_t numCellMarks() const noexcept {
    StormAssert(!CellRanges_.empty());
    return CellRanges_.size() - 1;
  }

  /// @brief Range of node indices with a @p nodeMark.
  auto nodeIndices(NodeMark nodeMark) const noexcept {
    StormAssert(nodeMark < numNodeMarks());
    return views::iota(NodeRanges_[nodeMark], NodeRanges_[nodeMark + 1]);
  }

  /// @brief Range of edge indices with a @p edgeMark.
  auto edgeIndices(EdgeMark edgeMark) const noexcept {
    StormAssert(edgeMark < numEdgeMarks());
    return views::iota(EdgeRanges_[edgeMark], EdgeRanges_[edgeMark + 1]);
  }

  /// @brief Range of face indices with a @p faceMark.
  auto faceIndices(FaceMark faceMark) const noexcept {
    StormAssert(faceMark < numFaceMarks());
    return views::iota(FaceRanges_[faceMark], FaceRanges_[faceMark + 1]);
  }

  /// @brief Range of cell indices with a @p cellMark.
  auto cellIndices(CellMark cellMark) const noexcept {
    StormAssert(cellMark < numCellMarks());
    return views::iota(CellRanges_[cellMark], CellRanges_[cellMark + 1]);
  }

  /// @brief Get node @p nodeIndex mark.
  NodeMark mark(NodeIndex nodeIndex) const noexcept {
    StormAssert(nodeIndex < NumNodes_);
    return NodeMarks_[nodeIndex];
  }

  /// @brief Get edge @p edgeIndex mark.
  EdgeMark mark(EdgeIndex edgeIndex) const noexcept {
    StormAssert(edgeIndex < NumEdges_);
    return EdgeMarks_[edgeIndex];
  }

  /// @brief Get face @p faceIndex mark.
  FaceMark mark(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces_);
    return FaceMarks_[faceIndex];
  }

  /// @brief Get cell @p cellIndex mark.
  CellMark mark(CellIndex cellIndex) const noexcept {
    StormAssert(cellIndex < NumCells_);
    return CellMarks_[cellIndex];
  }

  /// @}

  /// ---------------------------------------------------------------- ///
  /// @name Shapes.
  /// ---------------------------------------------------------------- ///
  /// @{

  /// @brief Get edge @p edgeIndex shape type.
  ShapeType shapeType(EdgeIndex edgeIndex) const noexcept {
    StormAssert(edgeIndex < NumEdges_);
    return EdgeShapeTypes_[edgeIndex];
  }

  /// @brief Get face @p faceIndex shape type.
  ShapeType shapeType(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces_);
    return FaceShapeTypes_[faceIndex];
  }

  /// @brief Get cell @p cellIndex shape type.
  ShapeType shapeType(CellIndex cellIndex) const noexcept {
    StormAssert(cellIndex < NumCells_);
    return CellShapeTypes_[cellIndex];
  }

private:

  auto makeShape_(ShapeDesc&& desc) const {
    return Element::Make(std::forward<ShapeDesc>(desc), NodePos_);
  }

public:

  /// @brief Get element object.
  template<class Tag>
  auto shape(Index<Tag> index) const {
    auto const nodeIndices = adjNodeIndices(index) |
      views::transform([](NodeIndex nodeIndex) {
        return static_cast<size_t>(nodeIndex);
      });
    return makeShape_({shapeType(index),
                        std::vector(nodeIndices.begin(), nodeIndices.end())});
  }

  /// @brief Get node @p nodeIndex position.
  vec3_t nodePos(NodeIndex nodeIndex) const noexcept {
    StormAssert(nodeIndex < NumNodes_);
    return NodePos_[nodeIndex];
  }

  /// @brief Set node @p nodeIndex position @p pos.
  void setNodePos(NodeIndex nodeIndex, vec3_t const& pos) noexcept {
    StormAssert(nodeIndex < NumNodes_);
    NodePos_[nodeIndex] = pos;
  }

  /// @brief Get edge @p edgeIndex length.
  real_t edgeLen(EdgeIndex edgeIndex) const noexcept {
    StormAssert(edgeIndex < NumEdges_);
    return EdgeLens_[edgeIndex];
  }

  /// @brief Get edge @p edgeIndex direction.
  vec3_t edgeDir(EdgeIndex edgeIndex) const noexcept {
    StormAssert(edgeIndex < NumEdges_);
    return EdgeDirs_[edgeIndex];
  }

  /// @brief Get face @p faceIndex area/length.
  real_t faceArea(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces_);
    return FaceAreas_[faceIndex];
  }

  /// @brief Get face @p faceIndex normal.
  vec3_t faceNormal(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces_);
    return FaceNormals_[faceIndex];
  }

  /// @brief Get face @p faceIndex barycenter.
  vec3_t faceCenterPos(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces_);
    return FaceCenterPos_[faceIndex];
  }

  /// @brief Get cell @p cellIndex volume/area/length.
  real_t cellVolume(CellIndex cellIndex) const noexcept {
    StormAssert(cellIndex < NumCells_);
    return CellVolumes_[cellIndex];
  }

  /// @brief Get cell @p cellIndex center position.
  vec3_t cellCenterPos(CellIndex cellIndex) const noexcept {
    StormAssert(cellIndex < NumCells_);
    return CellCenterPos_[cellIndex];
  }

  /// @brief Get maximal edge length.
  real_t maxEdgeLen() const noexcept {
    return MaxEdgeLen_;
  }

  /// @brief Get minimal edge length.
  real_t minEdgeLen() const noexcept {
    return MinEdgeLen_;
  }

  /// @brief Get maximal face area.
  real_t maxFaceArea() const noexcept {
    return MaxFaceArea_;
  }

  /// @brief Get minimal face area.
  real_t minFaceArea() const noexcept {
    return MinFaceArea_;
  }

  /// @brief Get maximal cell volume.
  real_t maxCellVolume() const noexcept {
    return MaxCellVolume_;
  }

  /// @brief Get minimal cell volume.
  real_t minCellVolume() const noexcept {
    return MinCellVolume_;
  }

  /// @brief Compute all elements geometry properties: \
  ///   edge lengths and directions; face areas, normals and \
  ///   center positions; cell volumes and center positions.
  ///
  /// This function should be called \
  ///   after each node position modification.
  void UpdateElementsGeometry();

  /// @}

  /// ---------------------------------------------------------------- ///
  /// @name Adjacency.
  /// ---------------------------------------------------------------- ///
  /// @{

  /// @name Primary adjacency: adjacency, \
  ///   that can be extracted from the shape only.
  /// @{

  /// @brief Range of the edge @p edgeIndex adjacent node indices.
  /// Denote a node to be adjacent to an edge if it one its nodes.
  template<class Self>
  auto adjNodeIndices(StormThis_ Self& self, EdgeIndex edgeIndex) noexcept {
    StormAssert(edgeIndex < self.NumEdges_);
    return self.EdgeNodes_[edgeIndex];
  }

  /// @brief Range of the face @p faceIndex adjacent node indices.
  /// Denote a node to be adjacent to a face if it one its nodes.
  template<class Self>
  auto adjNodeIndices(StormThis_ Self& self, FaceIndex faceIndex) noexcept {
    StormAssert(faceIndex < self.NumFaces_);
    return self.FaceNodes_[faceIndex];
  }

  /// @brief Range of the cell @p cellIndex adjacent node indices.
  /// Denote a node to be adjacent to a cell if it one its nodes.
  template<class Self>
  auto adjNodeIndices(StormThis_ Self& self, CellIndex cellIndex) noexcept {
    StormAssert(cellIndex < self.NumCells_);
    return self.CellNodes_[cellIndex];
  }

  /// @brief Range of the face @p faceIndex adjacent edge indices.
  /// Denote an edge to be adjacent to a face if it one its edges.
  template<class Self>
  auto adjEdgeIndices(StormThis_ Self& self, FaceIndex faceIndex) noexcept {
    StormAssert(faceIndex < self.NumFaces_);
    return self.FaceEdges_[faceIndex];
  }

  /// @brief Range of the cell @p cellIndex adjacent edge indices.
  /// Denote an edge to be adjacent to a cell if it one its edges.
  template<class Self>
  auto adjEdgeIndices(StormThis_ Self& self, CellIndex cellIndex) noexcept {
    StormAssert(cellIndex < self.NumCells_);
    return self.CellEdges_[cellIndex];
  }

  /// @brief Range of the cell @p cellIndex adjacent face indices.
  /// Denote a face to be adjacent to a cell if it one its faces.
  template<class Self>
  auto adjFaceIndices(StormThis_ Self& self, CellIndex cellIndex) noexcept {
    StormAssert(cellIndex < self.NumCells_);
    return self.CellFaces_[cellIndex];
  }

  /// @}

  /// @name Secondary adjacency: adjacency, \
  ///   that is a transpose of the primary adjacency.
  /// @{

  /// @brief Range of the node @p nodeIndex adjacent edge indices.
  template<class Self>
  auto adjEdgeIndices(StormThis_ Self& self, NodeIndex nodeIndex) noexcept {
    StormAssert(nodeIndex < self.NumNodes_);
    return self.EdgeNodes_.transposed()[nodeIndex];
  }

  /// @brief Range of the node @p nodeIndex adjacent faces indices.
  template<class Self>
  auto adjFaceIndices(StormThis_ Self& self, NodeIndex nodeIndex) noexcept {
    StormAssert(nodeIndex < self.NumNodes_);
    return self.FaceNodes_.transposed()[nodeIndex];
  }

  /// @brief Range of the node @p nodeIndex adjacent cell indices.
  template<class Self>
  auto adjCellIndices(StormThis_ Self& self, NodeIndex nodeIndex) noexcept {
    StormAssert(nodeIndex < self.NumNodes_);
    return self.CellNodes_.transposed()[nodeIndex];
  }

  /// @brief Range of the edge @p edgeIndex adjacent face indices.
  template<class Self>
  auto adjFaceIndices(StormThis_ Self& self, EdgeIndex edgeIndex) noexcept {
    StormAssert(edgeIndex < self.NumEdges_);
    return self.FaceEdges_.transposed()[edgeIndex];
  }

  /// @brief Range of the edge @p edgeIndex adjacent cell indices.
  template<class Self>
  auto adjCellIndices(StormThis_ Self& self, EdgeIndex edgeIndex) noexcept {
    StormAssert(edgeIndex < self.NumEdges_);
    return self.CellEdges_.transposed()[edgeIndex];
  }

  /// @brief Range of the face @p faceIndex adjacent cell indices.
  template<class Self>
  auto adjCellIndices(StormThis_ Self& self, FaceIndex faceIndex) noexcept {
    StormAssert(faceIndex < self.NumFaces_);
    return self.CellFaces_.transposed()[faceIndex];
  }

  /// @}

  /// @name Symmetric adjacency: adjacency, that is a product \
  ///   of the corresponding primary and secondary adjacencies.
  /// @{

  /// @brief Range of the edge @p edgeIndex adjacent edge indices.
  /// Denote two edges as adjacent if they share a common node.
  template<class Self>
  auto adjEdgeIndices(StormThis_ Self& self, EdgeIndex edgeIndex) noexcept {
    StormAssert(edgeIndex < self.NumEdges_);
    return self.EdgeNodes_.square()[edgeIndex];
  }

  /// @brief Range of the face @p faceIndex adjacent face indices.
  /// Denote two faces as adjacent if they share a common edge.
  template<class Self>
  auto adjFaceIndices(StormThis_ Self& self, FaceIndex faceIndex) noexcept {
    StormAssert(faceIndex < self.NumFaces_);
    return self.FaceEdges_.square()[faceIndex];
  }

  /// @brief Range of the cell @p cellIndex adjacent cell indices.
  /// Denote two cells as adjacent if they share a common face.
  template<class Self>
  auto adjCellIndices(StormThis_ Self& self, CellIndex cellIndex) noexcept {
    StormAssert(cellIndex < self.NumCells_);
    return self.CellFaces_.square()[cellIndex];
  }

  /// @}

  /// @brief Range of the node @p nodeIndex adjacent node indices.
  /// Denote two nodes as adjacent if there is an edge connecting them.
  template<class Self>
  auto adjNodeIndices(StormThis_ Self& self, NodeIndex nodeIndex) noexcept {
    StormAssert(nodeIndex < self.NumNodes_);
    return self.NodeNodes_.square()[nodeIndex];
  }

  StormDeducingThisEmulation_(adjNodeIndices)
  StormDeducingThisEmulation_(adjEdgeIndices)
  StormDeducingThisEmulation_(adjFaceIndices)
  StormDeducingThisEmulation_(adjCellIndices)

private:

  template<class OtherTag, class Self, class Tag>
  auto adjElemIndices_(StormThis_ Self& self, Index<Tag> elementIndex) noexcept {
    if constexpr (std::is_same_v<OtherTag, NodeTag>) {
      return self.adjNodeIndices(elementIndex);
    } else if constexpr (std::is_same_v<OtherTag, EdgeTag>) {
      return self.adjEdgeIndices(elementIndex);
    } else if constexpr (std::is_same_v<OtherTag, FaceTag>) {
      return self.adjFaceIndices(elementIndex);
    } else if constexpr (std::is_same_v<OtherTag, CellTag>) {
      return self.adjCellIndices(elementIndex);
    } else {
      static_assert(AlwaysFalse<OtherTag>, "Invalid tag.");
    }
  }

  StormDeducingThisEmulationT_(adjElemIndices_)

public:

  /// @}

  /// ---------------------------------------------------------------- ///
  /// @name Insertions.
  /// ---------------------------------------------------------------- ///
  /// @{

  /// @brief Insert a new node with a position @p pos \
  ///   and a mark @p nodeMark into the mesh.
  /// @returns Index of the inserted node.
  NodeIndex insertNode(vec3_t const& nodePos,
                       NodeMark nodeMark = {});

  /// @brief Find or emplace a new edge with a shape @p edgeShape
  ///   and node indices @p nodeIndices.
  ///
  /// Insertion would update the edge-node, edge-edge and
  ///   node-node topologies.
  ///
  /// @param edgeMark Mark that would be assigned to an edge if
  ///   the insertion took place, otherwise ignored.
  /// @returns A pair of an index of the found or inserted edge
  ///   and a boolean value denoting whether the insertion took place.
  std::pair<EdgeIndex, bool>
    findOrInsertEdge(std::shared_ptr<Shape> edgeShape,
                     EdgeMark edgeMark = {});

  /// @brief Find or emplace a new face with a shape @p faceShape
  ///   and node indices @p nodeIndices.
  ///
  /// Insertion would implicitly insert the missing edges and
  ///   update the face-edge and face-face topologies. The
  ///   implicitly inserted edges would inherit the mark from
  ///   explicitly inserted face.
  ///
  /// @param faceMark Mark that would be assigned to a face if
  ///   the insertion took place, otherwise ignored.
  /// @returns A pair of an index of the found or inserted face
  ///   and a boolean value denoting whether the insertion took place.
  std::pair<FaceIndex, bool>
    findOrInsertFace(std::shared_ptr<Shape> faceShape,
                     FaceMark faceMark = {});

  /// @brief Find or emplace a new cell with a shape @p cellShape
  ///   and node indices @p nodeIndices.
  ///
  /// Insertion would implicitly insert the missing faces and
  ///   update the cell-face and cell-cell topologies. The
  ///   implicitly inserted face would inherit the mark from
  ///   explicitly inserted cell.
  ///
  /// @param cellMark Mark that would be assigned to a cell if
  ///   the insertion took place, otherwise ignored.
  /// @returns A pair of an index of the found or inserted cell
  ///   and a boolean value denoting whether the insertion took place.
  std::pair<CellIndex, bool>
    findOrInsertCell(std::shared_ptr<Shape> cellShape,
                     CellMark cellMark = {});

  /// @}

  /// ---------------------------------------------------------------- ///
  /// ---------------------------------------------------------------- ///

  /// @brief Emplace a new node into the mesh.
  /// @returns Index of the inserted node.
  NodeIndex EmplaceNode(vec3_t const& nodePos, NodeMark nodeMark = {});

  /// @brief Emplace a new edge into the mesh.
  /// @returns Index of the inserted edge.
  /// @{
  EdgeIndex EmplaceEdge(std::unique_ptr<Element>&& edge, EdgeMark edgeMark = {});
  EdgeIndex EmplaceEdge(ShapeDesc&& edgeDesc, EdgeMark edgeMark = {}) {
    return EmplaceEdge(makeShape_(std::forward<ShapeDesc>(edgeDesc)), edgeMark);
  }
  /// @}

  /// @brief Emplace a new face into the mesh.
  /// @returns Index of the inserted face.
  /// @{
  FaceIndex EmplaceFace(std::unique_ptr<Element>&& face, FaceMark faceMark = {});
  FaceIndex EmplaceFace(ShapeDesc&& faceDesc, FaceMark faceMark = {}) {
    return EmplaceFace(makeShape_(std::forward<ShapeDesc>(faceDesc)), faceMark);
  }
  /// @}

  /// @brief Emplace a new cell into the mesh.
  /// @returns Index of the inserted cell.
  /// @{
  CellIndex EmplaceCell(std::unique_ptr<Element>&& cell, CellMark cellMark = {});
  CellIndex EmplaceCell(ShapeDesc&& cellDesc, CellMark cellMark = {}) {
    return EmplaceCell(makeShape_(std::forward<ShapeDesc>(cellDesc)), cellMark);
  }
  /// @}

private:

  template<class Tag>
  void FixPermutationAndAdjacency_(std::vector<size_t>& cellIndex);

public:

  /// @brief Change order of all nodes.
  void PermuteNodes(std::vector<size_t>&& nodePermutation = {});

  /// @brief Change order of all edges.
  void PermuteEdges(std::vector<size_t>&& edgePermutation = {});

  /// @brief Change order of all faces.
  void PermuteFaces(std::vector<size_t>&& facePermutation = {});

  /// @brief Change order of all cells.
  void PermuteCells(std::vector<size_t>&& cellPermutation = {});

protected:

  void reorder_faces();

  // ---------------------------------------------------------------------- //
  // ---------------------------------------------------------------------- //

  /// @brief Generate edges using the face to Node connectivity.
  /// @warning This function may be slow and memory-consuming.
  void FinalizeEdges_();

  /// @brief Generate faces using the cell to Node connectivity.
  /// @warning This function may be slow and memory-consuming.
  void FinalizeFaces_();

  /// @brief Generate boundary cells to complete face connectivity.
  void generate_boundary_cells();

}; // class Mesh

} // namespace feathers

// Include iterators.
#include "View.hxx"

/// @brief @todo Remove me.
using cMesh = feathers::Mesh;
