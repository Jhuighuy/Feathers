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

#include "SkunkBase.hh"
#include "Index.hh"
#include "libFeathersUtils/Table.hh"
#include "libFeathersUtils/Image.hh"
#include "Field.hh"
#include "Element.hh"

#include <ranges>
#include <map>

namespace feathers {

enum : size_t {
  FaceInnerCell_ = 0,
  FaceOuterCell_ = 1,
}; // enum

/// @todo Remove "{}".
class NodeTag_ {};
class EdgeTag_ {};
class FaceTag_ {};
class CellTag_ {};
template<class> class MarkTag_;

/// @brief Node index.
using NodeIndex = Index<NodeTag_>;

/// @brief Edge index.
using EdgeIndex = Index<EdgeTag_>;

/// @brief Face index.
using FaceIndex = Index<FaceTag_>;

/// @brief Cell index.
using CellIndex = Index<CellTag_>;

/// @brief Node mark index.
using NodeMark = Index<MarkTag_<NodeTag_>>;

/// @brief Edge mark index.
using EdgeMark = Index<MarkTag_<EdgeTag_>>;

/// @brief Face mark index.
using FaceMark = Index<MarkTag_<FaceTag_>>;

/// @brief Cell mark index.
using CellMark = Index<MarkTag_<CellTag_>>;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Hybrid unstructured multidimensional mesh.
/// @todo Use strict indices instead of pure tags.
/// @todo Do not use iterators in the mesh implementation.
/// @todo Swith to ranges.
/// @todo Make Dim a template parameter.
/// @todo Switch to size_t.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
class Mesh : public tObject<Mesh> {
public:
  size_t m_dim;

  size_t NumNodes_ = 0, NumEdges_ = 0;
  size_t NumFaces_ = 0, NumCells_ = 0;

  IndexedVector<NodeMark, NodeIndex> NodeMarks_;
  IndexedVector<EdgeMark, EdgeIndex> EdgeMarks_;
  IndexedVector<FaceMark, FaceIndex> FaceMarks_;
  IndexedVector<CellMark, CellIndex> CellMarks_;
  IndexedVector<NodeIndex, NodeMark> NodeRanges_{NodeIndex{0}};
  IndexedVector<EdgeIndex, EdgeMark> EdgeRanges_{EdgeIndex{0}};
  IndexedVector<FaceIndex, FaceMark> FaceRanges_{FaceIndex{0}};
  IndexedVector<CellIndex, CellMark> CellRanges_{CellIndex{0}};

  IndexedVector<eShape, EdgeIndex> EdgeShapes_;
  IndexedVector<eShape, FaceIndex> FaceShapes_;
  IndexedVector<eShape, CellIndex> CellShapes_;
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

  /// @brief Initialize an empty mesh.
  explicit Mesh(size_t dim = 0) : m_dim(dim) {
    StormAssert(0 <= m_dim && m_dim <= 3);
  }

  bool read_triangle(const char* path);
  bool read_tetgen(const char* path);
  bool read_image2D(const char* path,
                    const std::map<sPixel, size_t>& mark_colors,
                    sPixel fluid_color = eBlackPixel,
                    vec2_t pixel_size = vec2_t(1.0, 1.0));

  void save_vtk(const char* path,
                const std::vector<sFieldDesc>& fields) const;

  void finalize() {
    FinalizeFaces_();
    generate_boundary_cells();
    reorder_faces();
    ComputeShapeProperties();
  }

  /// @brief Total number of nodes in the mesh.
  size_t NumNodes() const noexcept {
    return NumNodes_;
  }

  /// @brief Total number of edges in the mesh.
  size_t NumEdges() const noexcept {
    return NumEdges_;
  }

  /// @brief Total number of faces in the mesh.
  size_t NumFaces() const noexcept {
    return NumFaces_;
  }

  /// @brief Total number of cells in the mesh.
  size_t NumCells() const noexcept {
    return NumCells_;
  }

  /// @brief Range of nodes.
  auto Nodes() const noexcept {
    return std::views::iota(NodeIndex{0}, NodeIndex{NumNodes()});
  }

  /// @brief Range of edges.
  auto Edges() const noexcept {
    return std::views::iota(EdgeIndex{0}, EdgeIndex{NumEdges()});
  }

  /// @brief Range of faces.
  auto Faces() const noexcept {
    return std::views::iota(FaceIndex{0}, FaceIndex{NumFaces()});
  }

  /// @brief Range of cells.
  auto Cells() const noexcept {
    return std::views::iota(CellIndex{0}, CellIndex{NumCells()});
  }

  /// @brief Number of node marks.
  size_t NumNodeMarks() const noexcept {
    StormAssert(!NodeRanges_.empty());
    return NodeRanges_.size() - 1;
  }

  /// @brief Number of edge marks. 
  size_t NumEdgeMarks() const noexcept {
    StormAssert(!EdgeRanges_.empty());
    return EdgeRanges_.size() - 1;
  }

  /// @brief Number of face marks. 
  size_t NumFaceMarks() const noexcept {
    StormAssert(!FaceRanges_.empty());
    return FaceRanges_.size() - 1;
  }

  /// @brief Number of cell marks. 
  size_t NumCellMarks() const noexcept {
    StormAssert(!CellRanges_.empty());
    return CellRanges_.size() - 1;
  }

  /// @brief Range of nodes with a given mark.
  auto Nodes(NodeMark nodeMark) const noexcept {
    StormAssert(nodeMark < NumNodeMarks());
    return std::views::iota(
      NodeRanges_[nodeMark], NodeRanges_[nodeMark + 1]);
  }

  /// @brief Range of edges with a given mark.
  auto Edges(EdgeMark edgeMark) const noexcept {
    StormAssert(edgeMark < NumEdgeMarks());
    return std::views::iota(
      EdgeRanges_[edgeMark], EdgeRanges_[edgeMark + 1]);
  }

  /// @brief Range of faces with a given mark.
  auto Faces(FaceMark faceMark) const noexcept {
    StormAssert(faceMark < NumFaceMarks());
    return std::views::iota(
      FaceRanges_[faceMark], FaceRanges_[faceMark + 1]);
  }

  /// @brief Range of cells with a given mark.
  auto Cells(CellMark cellMark) const noexcept {
    StormAssert(cellMark < NumCellMarks());
    return std::views::iota(
      CellRanges_[cellMark], CellRanges_[cellMark + 1]);
  }

  /// @brief Number of marked cells in the mesh.
  size_t NumCells(CellMark cellMark) const noexcept {
    return std::size(Cells(cellMark));
  }

#if 1 // -------------------------------------------------------------------

  /// @brief Index of the first node with a given mark.
  FEATHERS_DEPRECATED NodeIndex BeginNode(NodeMark nodeMark) const noexcept {
    StormAssert(nodeMark < NumNodeMarks());
    return NodeRanges_[nodeMark];
  }

  /// @brief Index of the first edge with a given mark.
  FEATHERS_DEPRECATED EdgeIndex BeginEdge(EdgeMark edgeMark) const noexcept {
    StormAssert(edgeMark < NumEdgeMarks());
    return EdgeRanges_[edgeMark];
  }

  /// @brief Index of the first face with a given mark.
  FEATHERS_DEPRECATED FaceIndex BeginFace(FaceMark faceMark) const noexcept {
    StormAssert(faceMark < NumFaceMarks());
    return FaceRanges_[faceMark];
  }

  /// @brief Index of the first cell with a given mark.
  FEATHERS_DEPRECATED CellIndex BeginCell(CellMark cellMark) const noexcept {
    StormAssert(cellMark < NumCellMarks());
    return CellRanges_[cellMark];
  }

  /// @brief Index of a node after the last node with a given mark.
  FEATHERS_DEPRECATED NodeIndex EndNode(NodeMark nodeMark) const noexcept {
    StormAssert(nodeMark < NumNodeMarks());
    return NodeRanges_[nodeMark + 1];
  }

  /// @brief Index of an edge after the last edge with a given mark.
  FEATHERS_DEPRECATED EdgeIndex EndEdge(EdgeMark edgeMark) const noexcept {
    StormAssert(edgeMark < NumEdgeMarks());
    return EdgeRanges_[edgeMark + 1];
  }

  /// @brief Index of a face after the last face with a given mark.
  FEATHERS_DEPRECATED FaceIndex EndFace(FaceMark faceMark) const noexcept {
    StormAssert(faceMark < NumFaceMarks());
    return FaceRanges_[faceMark + 1];
  }

  /// @brief Index a cell after of the last cell with a given cellMark.
  FEATHERS_DEPRECATED CellIndex EndCell(CellMark cellMark) const noexcept {
    StormAssert(cellMark < NumCellMarks());
    return CellRanges_[cellMark + 1];
  }

#endif // ------------------------------------------------------------------

  /// @brief Get element mark.
  /// @{
  NodeMark Mark(NodeIndex nodeIndex) const noexcept {
    StormAssert(nodeIndex < NumNodes());
    return NodeMarks_[nodeIndex];
  }
  EdgeMark Mark(EdgeIndex edgeIndex) const noexcept {
    StormAssert(edgeIndex < NumEdges());
    return EdgeMarks_[edgeIndex];
  }
  FaceMark Mark(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces());
    return FaceMarks_[faceIndex];
  }
  CellMark Mark(CellIndex cellIndex) const noexcept {
    StormAssert(cellIndex < NumCells());
    return CellMarks_[cellIndex];
  }
  /// @}

  /// @brief Get element shape.
  /// @{
  eShape Shape(EdgeIndex edgeIndex) const noexcept {
    StormAssert(edgeIndex < NumEdges());
    return EdgeShapes_[edgeIndex];
  }
  eShape Shape(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces());
    return FaceShapes_[faceIndex];
  }
  eShape Shape(CellIndex cellIndex) const noexcept {
    StormAssert(cellIndex < NumCells());
    return CellShapes_[cellIndex];
  }
  /// @}

  /// @brief Make a mesh element from description.
  FEATHERS_CONST_OVERLOAD_R(
    std::unique_ptr<iElement>,
    std::unique_ptr<iElement const>, make_element, (sElementDesc&& desc), {
      return iElement::make(std::forward<sElementDesc>(desc), NumNodes_, NodePos_.data());
  })

  /// @brief Get element object.
  FEATHERS_CONST_OVERLOAD_R_T(
    template<class Tag>,
    std::unique_ptr<iElement>,
    std::unique_ptr<const iElement>, get_object, (Index<Tag> index), {
      return make_element(
        {Shape(index), std::vector<size_t>(std::begin(AdjacentNodes(index)),
                                           std::end(AdjacentNodes(index)))});
  })

  /// @brief Compute edge shape properties.
  void ComputeEdgeShapeProperties();

  /// @brief Compute face shape properties.
  void ComputeFaceShapeProperties();

  /// @brief Compute cell shape properties.
  void ComputeCellShapeProperties();

  /// @brief Compute all elements shape properties.
  void ComputeShapeProperties() {
    ComputeEdgeShapeProperties();
    ComputeFaceShapeProperties();
    ComputeCellShapeProperties();
  }

  /// @brief Get node position.
  vec3_t const& NodePos(NodeIndex nodeIndex) const noexcept {
    StormAssert(nodeIndex < NumNodes());
    return NodePos_[nodeIndex];
  }

  /// @brief Set node position.
  void SetNodePos(NodeIndex nodeIndex, vec3_t const& pos) noexcept {
    StormAssert(nodeIndex < NumNodes());
    NodePos_[nodeIndex] = pos;
  }

  /// @brief Get edge length.
  real_t EdgeLen(EdgeIndex edgeIndex) const noexcept {
    StormAssert(edgeIndex < NumEdges());
    return EdgeLens_[edgeIndex];
  }

  /// @brief Set edge length.
  /// @todo Remove me!
  FEATHERS_DEPRECATED void SetEdgeLen(EdgeIndex edgeIndex, real_t len) noexcept {
    StormAssert(edgeIndex < NumEdges());
    EdgeLens_[edgeIndex] = len;
  }

  /// @brief Get edge direction.
  vec3_t EdgeDir(EdgeIndex edgeIndex) const noexcept {
    StormAssert(edgeIndex < NumEdges());
    return EdgeDirs_[edgeIndex];
  }

  /// @brief Set edge direction.
  /// @todo Remove me!
  FEATHERS_DEPRECATED void SetEdgeDir(EdgeIndex edgeIndex, vec3_t const& dir) noexcept {
    StormAssert(edgeIndex < NumEdges());
    EdgeDirs_[edgeIndex] = dir;
  }

  /// @brief Get face area/length.
  real_t FaceArea(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces());
    return FaceAreas_[faceIndex];
  }

  /// @brief Set face area/length.
  FEATHERS_DEPRECATED void SetFaceArea(FaceIndex faceIndex, real_t area) noexcept {
    StormAssert(faceIndex < NumFaces());
    FaceAreas_[faceIndex] = area;
  }

  /// @brief Get face normal.
  vec3_t FaceNormal(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces());
    return FaceNormals_[faceIndex];
  }

  /// @brief Set face normal.
  /// @todo Remove me!
  FEATHERS_DEPRECATED void SetFaceNormal(FaceIndex faceIndex, vec3_t const& normal) noexcept {
    StormAssert(faceIndex < NumFaces());
    FaceNormals_[faceIndex] = normal;
  }

  /// @brief Get face barycenter.
  vec3_t FaceCenterPos(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces());
    return FaceCenterPos_[faceIndex];
  }

  /// @brief Set face barycenter.
  /// @todo Remove me!
  FEATHERS_DEPRECATED void SetFaceCenterPos(FaceIndex faceIndex, vec3_t const& centerPos) noexcept {
    StormAssert(faceIndex < NumFaces());
    FaceCenterPos_[faceIndex] = centerPos;
  }

  /// @brief Get cell volume/area/length.
  real_t CellVolume(CellIndex cellIndex) const noexcept {
    StormAssert(cellIndex < NumCells());
    return CellVolumes_[cellIndex];
  }

  /// @brief Set cell volume/area/length.
  /// @todo Remove me!
  FEATHERS_DEPRECATED void SetCellVolume(CellIndex cellIndex, real_t volume) noexcept {
    StormAssert(cellIndex < NumCells());
    CellVolumes_[cellIndex] = volume;
  }

  /// @brief Get cell barycenter.
  vec3_t CellCenterPos(CellIndex cellIndex) const noexcept {
    StormAssert(cellIndex < NumCells());
    return CellCenterPos_[cellIndex];
  }

  /// @brief Get cell barycenter.
  /// @todo Remove me!
  FEATHERS_DEPRECATED void SetCellCenterPos(CellIndex cellIndex, vec3_t const& centerPos) noexcept {
    StormAssert(cellIndex < NumCells());
    CellCenterPos_[cellIndex] = centerPos;
  }

  /// @brief Get maximal edge length.
  real_t MaxEdgeLen() const noexcept {
    return MaxEdgeLen_;
  }

  /// @brief Get minimal edge length.
  real_t MinEdgeLen() const noexcept {
    return MinEdgeLen_;
  }

  /// @brief Get maximal face area.
  real_t MaxFaceArea() const noexcept {
    return MaxFaceArea_;
  }

  /// @brief Get minimal face area.
  real_t MinFaceArea() const noexcept {
    return MinFaceArea_;
  }

  /// @brief Get maximal cell volume.
  real_t MaxCellVolume() const noexcept {
    return MaxCellVolume_;
  }

  /// @brief Get minimal cell volume.
  real_t MinCellVolume() const noexcept {
    return MinCellVolume_;
  }

  /// @brief Range of the element adjacent nodes.
  /// @{
  ConstOverload(auto, AdjacentNodes, (NodeIndex nodeIndex), noexcept {
    StormAssert(nodeIndex < NumNodes());
    return NodeNodes_[nodeIndex];
  })
  ConstOverload(auto, AdjacentNodes, (EdgeIndex edgeIndex), noexcept {
    StormAssert(edgeIndex < NumEdges());
    return EdgeNodes_[edgeIndex];
  })
  ConstOverload(auto, AdjacentNodes, (FaceIndex faceIndex), noexcept {
    StormAssert(faceIndex < NumFaces());
    return FaceNodes_[faceIndex];
  })
  ConstOverload(auto, AdjacentNodes, (CellIndex cellIndex), noexcept {
    StormAssert(cellIndex < NumCells());
    return CellNodes_[cellIndex];
  })
  /// @}

  /// @brief Range of the element adjacent edges.
  /// @{
  ConstOverload(auto, AdjacentEdges, (NodeIndex nodeIndex), noexcept {
    StormAssert(nodeIndex < NumNodes());
    return NodeEdges_[nodeIndex];
  })
  ConstOverload(auto, AdjacentEdges, (EdgeIndex edgeIndex), noexcept {
    StormAssert(edgeIndex < NumEdges());
    return EdgeEdges_[edgeIndex];
  })
  ConstOverload(auto, AdjacentEdges, (FaceIndex faceIndex), noexcept {
    StormAssert(faceIndex < NumFaces());
    return FaceEdges_[faceIndex];
  })
  ConstOverload(auto, AdjacentEdges, (CellIndex cellIndex), noexcept {
    StormAssert(cellIndex < NumCells());
    return CellEdges_[cellIndex];
  })
  /// @}

  /// @brief Range of the element adjacent faces.
  /// @{
  ConstOverload(auto, AdjacentFaces, (NodeIndex nodeIndex), noexcept {
    StormAssert(nodeIndex < NumNodes());
    return NodeFaces_[nodeIndex];
  })
  ConstOverload(auto, AdjacentFaces, (EdgeIndex edgeIndex), noexcept {
    StormAssert(edgeIndex < NumEdges());
    return EdgeFaces_[edgeIndex];
  })
  ConstOverload(auto, AdjacentFaces, (FaceIndex faceIndex), noexcept {
    StormAssert(faceIndex < NumFaces());
    return FaceFaces_[faceIndex];
  })
  ConstOverload(auto, AdjacentFaces, (CellIndex cellIndex), noexcept {
    StormAssert(cellIndex < NumCells());
    return CellFaces_[cellIndex];
  })
  /// @}

  /// @brief Range of the element adjacent cells.
  /// @{
  ConstOverload(auto, AdjacentCells, (NodeIndex nodeIndex), noexcept {
    StormAssert(nodeIndex < NumNodes());
    return NodeCells_[nodeIndex];
  })
  ConstOverload(auto, AdjacentCells, (EdgeIndex edgeIndex), noexcept {
    StormAssert(edgeIndex < NumEdges());
    return EdgeCells_[edgeIndex];
  })
  ConstOverload(auto, AdjacentCells, (FaceIndex faceIndex), noexcept {
    StormAssert(faceIndex < NumFaces());
    return FaceCells_[faceIndex];
  })
  ConstOverload(auto, AdjacentCells, (CellIndex cellIndex), noexcept {
    StormAssert(cellIndex < NumCells());
    return CellCells_[cellIndex];
  })

private:

  /// @brief Ranges of the adjacent elements.
  template<class OtherTag, class Tag>
  auto AdjacentElements_(Index<Tag> elementIndex) noexcept {
    if constexpr (std::is_same_v<OtherTag, NodeTag_>) {
      return AdjacentNodes(elementIndex);
    } else if constexpr (std::is_same_v<OtherTag, EdgeTag_>) {
      return AdjacentEdges(elementIndex);
    } else if constexpr (std::is_same_v<OtherTag, FaceTag_>) {
      return AdjacentFaces(elementIndex);
    } else if constexpr (std::is_same_v<OtherTag, CellTag_>) {
      return AdjacentCells(elementIndex);
    }
  }

public:

  /// @brief Emplace a new node into the mesh.
  /// @returns Index of the inserted node.
  NodeIndex EmplaceNode(vec3_t const& nodePos, NodeMark nodeMark = {});

  /// @brief Emplace a new edge into the mesh.
  /// @returns Index of the inserted edge.
  /// @{
  EdgeIndex EmplaceEdge(std::unique_ptr<iElement>&& edge, EdgeMark edgeMark = {});
  EdgeIndex EmplaceEdge(sElementDesc&& edgeDesc, EdgeMark edgeMark = {}) {
    return EmplaceEdge(make_element(std::forward<sElementDesc>(edgeDesc)), edgeMark);
  }
  /// @}

  /// @brief Emplace a new face into the mesh.
  /// @returns Index of the inserted face.
  /// @{
  FaceIndex EmplaceFace(std::unique_ptr<iElement>&& face, FaceMark faceMark = {});
  FaceIndex EmplaceFace(sElementDesc&& faceDesc, FaceMark faceMark = {}) {
    return EmplaceFace(make_element(std::forward<sElementDesc>(faceDesc)), faceMark);
  }
  /// @}

  /// @brief Emplace a new cell into the mesh.
  /// @returns Index of the inserted cell.
  /// @{
  CellIndex EmplaceCell(std::unique_ptr<iElement>&& cell, CellMark cellMark = {});
  CellIndex EmplaceCell(sElementDesc&& cellDesc, CellMark cellMark = {}) {
    return EmplaceCell(make_element(std::forward<sElementDesc>(cellDesc)), cellMark);
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

  /// @brief Generate edges using the face to node connectivity.
  /// @warning This function may be slow and memory-consuming.
  void FinalizeEdges_();

  /// @brief Generate faces using the cell to node connectivity.
  /// @warning This function may be slow and memory-consuming.
  void FinalizeFaces_();

  /// @brief Generate boundary cells to complete face connectivity.
  void generate_boundary_cells();

}; // class Mesh

} // namespace feathers

// Include iterators.
#include "MeshReferences.hh"

/// @brief @todo Remove me. 
using cMesh = feathers::Mesh;
#include "Field.hh"
