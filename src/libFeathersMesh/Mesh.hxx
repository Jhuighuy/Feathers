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

  CsrTable<NodeIndex, NodeIndex> NodeNodes_;
  CsrTable<EdgeIndex, NodeIndex> EdgeNodes_;
  CsrTable<FaceIndex, NodeIndex> FaceNodes_;
  CsrTable<CellIndex, NodeIndex> CellNodes_;

  CsrTable<NodeIndex, EdgeIndex> NodeEdges_;
  CsrTable<EdgeIndex, EdgeIndex> EdgeEdges_;
  CsrTable<FaceIndex, EdgeIndex> FaceEdges_;
  CsrTable<CellIndex, EdgeIndex> CellEdges_;

  CsrTable<NodeIndex, FaceIndex> NodeFaces_;
  CsrTable<EdgeIndex, FaceIndex> EdgeFaces_;
  CsrTable<FaceIndex, FaceIndex> FaceFaces_;
  CsrTable<CellIndex, FaceIndex> CellFaces_;

  CsrTable<NodeIndex, CellIndex> NodeCells_;
  CsrTable<EdgeIndex, CellIndex> EdgeCells_;
  CsrTable<FaceIndex, CellIndex> FaceCells_;
  CsrTable<CellIndex, CellIndex> CellCells_;

public:

  // ---------------------------------------------------------------- //
  // Constructors.
  // ---------------------------------------------------------------- //

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

  // ---------------------------------------------------------------- //
  // Element ranges.
  // ---------------------------------------------------------------- //

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

  // ---------------------------------------------------------------- //
  // Marks.
  // ---------------------------------------------------------------- //

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
  NodeMark Mark(NodeIndex nodeIndex) const noexcept {
    StormAssert(nodeIndex < NumNodes_);
    return NodeMarks_[nodeIndex];
  }

  /// @brief Get edge @p edgeIndex mark.
  EdgeMark Mark(EdgeIndex edgeIndex) const noexcept {
    StormAssert(edgeIndex < NumEdges_);
    return EdgeMarks_[edgeIndex];
  }

  /// @brief Get face @p faceIndex mark.
  FaceMark Mark(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces_);
    return FaceMarks_[faceIndex];
  }

  /// @brief Get cell @p cellIndex mark.
  CellMark Mark(CellIndex cellIndex) const noexcept {
    StormAssert(cellIndex < NumCells_);
    return CellMarks_[cellIndex];
  }

  // ---------------------------------------------------------------- //
  // Shapes.
  // ---------------------------------------------------------------- //

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

  auto MakeElement_(ElementDesc&& desc) const {
    return Element::Make(std::forward<ElementDesc>(desc), NodePos_);
  }

public:

  /// @brief Get element object.
  template<class Tag>
  auto get_object(Index<Tag> index) const {
    auto const nodeIndices =
      adjNodeIndices(index) | views::transform([](NodeIndex nodeIndex) {
        return static_cast<size_t>(nodeIndex);
      });
    return MakeElement_({shapeType(index),
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

  // ---------------------------------------------------------------- //
  // Adjacency.
  // ---------------------------------------------------------------- //

  /// @brief Range of the node @p nodeIndex adjacent nodes.
  StormAutoConstOverload_(adjNodeIndices, (NodeIndex nodeIndex), noexcept {
    StormAssert(nodeIndex < NumNodes_);
    return NodeNodes_[nodeIndex];
  })

  /// @brief Range of the edge @p edgeIndex adjacent nodes.
  StormAutoConstOverload_(adjNodeIndices, (EdgeIndex edgeIndex), noexcept {
    StormAssert(edgeIndex < NumEdges_);
    return EdgeNodes_[edgeIndex];
  })

  /// @brief Range of the face @p faceIndex adjacent nodes.
  StormAutoConstOverload_(adjNodeIndices, (FaceIndex faceIndex), noexcept {
    StormAssert(faceIndex < NumFaces_);
    return FaceNodes_[faceIndex];
  })

  /// @brief Range of the cell @p cellIndex adjacent nodes.
  StormAutoConstOverload_(adjNodeIndices, (CellIndex cellIndex), noexcept {
    StormAssert(cellIndex < NumCells_);
    return CellNodes_[cellIndex];
  })

  /// @brief Range of the node @p nodeIndex adjacent edges.
  StormAutoConstOverload_(adjEdgeIndices, (NodeIndex nodeIndex), noexcept {
    StormAssert(nodeIndex < NumNodes_);
    return NodeEdges_[nodeIndex];
  })

  /// @brief Range of the edge @p edgeIndex adjacent edges.
  StormAutoConstOverload_(adjEdgeIndices, (EdgeIndex edgeIndex), noexcept {
    StormAssert(edgeIndex < NumEdges_);
    return EdgeEdges_[edgeIndex];
  })

  /// @brief Range of the face @p faceIndex adjacent edges.
  StormAutoConstOverload_(adjEdgeIndices, (FaceIndex faceIndex), noexcept {
    StormAssert(faceIndex < NumFaces_);
    return FaceEdges_[faceIndex];
  })

  /// @brief Range of the cell @p cellIndex adjacent edges.
  StormAutoConstOverload_(adjEdgeIndices, (CellIndex cellIndex), noexcept {
    StormAssert(cellIndex < NumCells_);
    return CellEdges_[cellIndex];
  })

  /// @brief Range of the node @p nodeIndex adjacent faces.
  StormAutoConstOverload_(adjFaceIndices, (NodeIndex nodeIndex), noexcept {
    StormAssert(nodeIndex < NumNodes_);
    return NodeFaces_[nodeIndex];
  })

  /// @brief Range of the edge @p edgeIndex adjacent faces.
  StormAutoConstOverload_(adjFaceIndices, (EdgeIndex edgeIndex), noexcept {
    StormAssert(edgeIndex < NumEdges_);
    return EdgeFaces_[edgeIndex];
  })

  /// @brief Range of the face @p faceIndex adjacent faces.
  StormAutoConstOverload_(adjFaceIndices, (FaceIndex faceIndex), noexcept {
    StormAssert(faceIndex < NumFaces_);
    return FaceFaces_[faceIndex];
  })

  /// @brief Range of the cell @p cellIndex adjacent faces.
  StormAutoConstOverload_(adjFaceIndices, (CellIndex cellIndex), noexcept {
    StormAssert(cellIndex < NumCells_);
    return CellFaces_[cellIndex];
  })

  /// @brief Range of the node @p nodeIndex adjacent cells.
  StormAutoConstOverload_(adjCellIndices, (NodeIndex nodeIndex), noexcept {
    StormAssert(nodeIndex < NumNodes_);
    return NodeCells_[nodeIndex];
  })

  /// @brief Range of the edge @p edgeIndex adjacent cells.
  StormAutoConstOverload_(adjCellIndices, (EdgeIndex edgeIndex), noexcept {
    StormAssert(edgeIndex < NumEdges_);
    return EdgeCells_[edgeIndex];
  })

  /// @brief Range of the face @p faceIndex adjacent cells.
  StormAutoConstOverload_(adjCellIndices, (FaceIndex faceIndex), noexcept {
    StormAssert(faceIndex < NumFaces_);
    return FaceCells_[faceIndex];
  })

  /// @brief Range of the cell @p cellIndex adjacent cells.
  StormAutoConstOverload_(adjCellIndices, (CellIndex cellIndex), noexcept {
    StormAssert(cellIndex < NumCells_);
    return CellCells_[cellIndex];
  })

private:

  StormAutoConstOverloadT_(StormPass_(template<class OtherTag, class Tag>),
      AdjacentElements_, (Index<Tag> elementIndex), noexcept {
    if constexpr (std::is_same_v<OtherTag, NodeTag>) {
      return adjNodeIndices(elementIndex);
    } else if constexpr (std::is_same_v<OtherTag, EdgeTag>) {
      return adjEdgeIndices(elementIndex);
    } else if constexpr (std::is_same_v<OtherTag, FaceTag>) {
      return adjFaceIndices(elementIndex);
    } else if constexpr (std::is_same_v<OtherTag, CellTag>) {
      return adjCellIndices(elementIndex);
    }
  })

public:

  // ---------------------------------------------------------------- //
  // ---------------------------------------------------------------- //

  /// @brief Emplace a new node into the mesh.
  /// @returns Index of the inserted node.
  NodeIndex EmplaceNode(vec3_t const& nodePos, NodeMark nodeMark = {});

  /// @brief Emplace a new edge into the mesh.
  /// @returns Index of the inserted edge.
  /// @{
  EdgeIndex EmplaceEdge(std::unique_ptr<Element>&& edge, EdgeMark edgeMark = {});
  EdgeIndex EmplaceEdge(ElementDesc&& edgeDesc, EdgeMark edgeMark = {}) {
    return EmplaceEdge(MakeElement_(std::forward<ElementDesc>(edgeDesc)), edgeMark);
  }
  /// @}

  /// @brief Emplace a new face into the mesh.
  /// @returns Index of the inserted face.
  /// @{
  FaceIndex EmplaceFace(std::unique_ptr<Element>&& face, FaceMark faceMark = {});
  FaceIndex EmplaceFace(ElementDesc&& faceDesc, FaceMark faceMark = {}) {
    return EmplaceFace(MakeElement_(std::forward<ElementDesc>(faceDesc)), faceMark);
  }
  /// @}

  /// @brief Emplace a new cell into the mesh.
  /// @returns Index of the inserted cell.
  /// @{
  CellIndex EmplaceCell(std::unique_ptr<Element>&& cell, CellMark cellMark = {});
  CellIndex EmplaceCell(ElementDesc&& cellDesc, CellMark cellMark = {}) {
    return EmplaceCell(MakeElement_(std::forward<ElementDesc>(cellDesc)), cellMark);
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
