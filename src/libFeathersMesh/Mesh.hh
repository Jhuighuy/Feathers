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

#include <map>

namespace feathers {

enum : uint_t {
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
using NodeIndex = Index<uint_t, NodeTag_>;

/// @brief Edge index.
using EdgeIndex = Index<uint_t, EdgeTag_>;

/// @brief Face index.
using FaceIndex = Index<uint_t, FaceTag_>;

/// @brief Cell index.
using CellIndex = Index<uint_t, CellTag_>;

/// @brief Node mark index.
using NodeMark = Index<uint_t, MarkTag_<NodeTag_>>;

/// @brief Edge mark index.
using EdgeMark = Index<uint_t, MarkTag_<EdgeTag_>>;

/// @brief Face mark index.
using FaceMark = Index<uint_t, MarkTag_<FaceTag_>>;

/// @brief Cell mark index.
using CellMark = Index<uint_t, MarkTag_<CellTag_>>;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Hybrid unstructured multidimensional mesh.
/// @todo Use strict indices instead of pure tags.
/// @todo Make Dim a template parameter.
/// @todo Switch to size_t.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
class Mesh : public tObject<Mesh> {
public:
  uint_t m_dim;

  size_t NumNodes_ = 0, NumEdges_ = 0;
  size_t NumFaces_ = 0, NumCells_ = 0;

  IndexedVector<NodeMark, NodeIndex> NodeMarks_;
  IndexedVector<EdgeMark, EdgeIndex> EdgeMarks_;
  IndexedVector<FaceMark, FaceIndex> FaceMarks_;
  IndexedVector<CellMark, CellIndex> CellMarks_;
  IndexedVector<NodeIndex, NodeMark> NodeRanges_{NodeIndex(0)};
  IndexedVector<EdgeIndex, EdgeMark> EdgeRanges_{EdgeIndex(0)};
  IndexedVector<FaceIndex, FaceMark> FaceRanges_{FaceIndex(0)};
  IndexedVector<CellIndex, CellMark> CellRanges_{CellIndex(0)};

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
  explicit Mesh(uint_t dim = 0) : m_dim(dim) {
    StormAssert(0 <= m_dim && m_dim <= 3);
  }

  bool read_triangle(const char* path);
  bool read_tetgen(const char* path);
  bool read_image2D(const char* path,
                    const std::map<sPixel, uint_t>& mark_colors,
                    sPixel fluid_color = eBlackPixel,
                    vec2_t pixel_size = vec2_t(1.0, 1.0));

  void save_vtk(const char* path,
                const std::vector<sFieldDesc>& fields) const;

  void finalize() {
    generate_faces();
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

  /// @brief Index of the first node with a given mark.
  NodeIndex BeginNode(NodeMark nodeMark) const noexcept {
    StormAssert(nodeMark < NumNodeMarks());
    return NodeRanges_[nodeMark];
  }

  /// @brief Index of the first edge with a given mark.
  EdgeIndex BeginEdge(EdgeMark edgeMark) const noexcept {
    StormAssert(edgeMark < NumEdgeMarks());
    return EdgeRanges_[edgeMark];
  }

  /// @brief Index of the first face with a given mark.
  FaceIndex BeginFace(FaceMark faceMark) const noexcept {
    StormAssert(faceMark < NumFaceMarks());
    return FaceRanges_[faceMark];
  }

  /// @brief Index of the first cell with a given mark.
  CellIndex BeginCell(CellMark cellMark) const noexcept {
    StormAssert(cellMark < NumCellMarks());
    return CellRanges_[cellMark];
  }

  /// @brief Index of a node after the last node with a given mark.
  NodeIndex EndNode(NodeMark nodeMark) const noexcept {
    StormAssert(nodeMark < NumNodeMarks());
    return NodeRanges_[nodeMark + 1];
  }

  /// @brief Index of an edge after the last edge with a given mark.
  EdgeIndex EndEdge(EdgeMark edgeMark) const noexcept {
    StormAssert(edgeMark < NumEdgeMarks());
    return EdgeRanges_[edgeMark + 1];
  }

  /// @brief Index of a face after the last face with a given mark.
  FaceIndex EndFace(FaceMark faceMark) const noexcept {
    StormAssert(faceMark < NumFaceMarks());
    return FaceRanges_[faceMark + 1];
  }

  /// @brief Index a cell after of the last cell with a given cellMark.
  CellIndex EndCell(CellMark cellMark) const noexcept {
    StormAssert(cellMark < NumCellMarks());
    return CellRanges_[cellMark + 1];
  }

  /// @brief Number of marked cells in the mesh. 
  uint_t NumCells(CellMark mark) const noexcept {
    return EndCell(mark) - BeginCell(mark);
  }

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
    std::unique_ptr<const iElement>, get_object, (Index<uint_t, Tag> index), {
      return make_element(
        {Shape(index), std::vector<uint_t>(BeginAdjacentNode(index),
                                             EndAdjacentNode(index))});
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
  void SetEdgeLen(EdgeIndex edgeIndex, real_t len) noexcept {
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
  void SetEdgeDir(EdgeIndex edgeIndex, vec3_t const& dir) noexcept {
    StormAssert(edgeIndex < NumEdges());
    EdgeDirs_[edgeIndex] = dir;
  }

  /// @brief Get face area/length.
  real_t FaceArea(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces());
    return FaceAreas_[faceIndex];
  }

  /// @brief Set face area/length.
  void SetFaceArea(FaceIndex faceIndex, real_t area) noexcept {
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
  void SetFaceNormal(FaceIndex faceIndex, vec3_t const& normal) noexcept {
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
  void SetFaceCenterPos(FaceIndex faceIndex, vec3_t const& centerPos) noexcept {
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
  void SetCellVolume(CellIndex cellIndex, real_t volume) noexcept {
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
  void SetCellCenterPos(CellIndex cellIndex, vec3_t const& centerPos) noexcept {
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

  /// @brief Pointer to the beginning of the element adjacent nodes.
  /// @{
  ConstOverload(NodeIndex*, BeginAdjacentNode, (NodeIndex nodeIndex), {
    StormAssert(nodeIndex < NumNodes());
    return NodeNodes_.begin_row(nodeIndex);
  })
  ConstOverload(NodeIndex*, BeginAdjacentNode, (EdgeIndex edgeIndex), {
    StormAssert(edgeIndex < NumEdges());
    return EdgeNodes_.begin_row(edgeIndex);
  })
  ConstOverload(NodeIndex*, BeginAdjacentNode, (FaceIndex faceIndex), {
    StormAssert(faceIndex < NumFaces());
    return FaceNodes_.begin_row(faceIndex);
  })
  ConstOverload(NodeIndex*, BeginAdjacentNode, (CellIndex cellIndex), {
    StormAssert(cellIndex < NumCells());
    return CellNodes_.begin_row(cellIndex);
  })
  /// @}

  /// @brief Pointer to the beginning of the element adjacent edges.
  /// @{
  ConstOverload(EdgeIndex*, BeginAdjacentEdge, (NodeIndex nodeIndex), {
    StormAssert(nodeIndex < NumNodes());
    return NodeEdges_.begin_row(nodeIndex);
  })
  ConstOverload(EdgeIndex*, BeginAdjacentEdge, (EdgeIndex edgeIndex), {
    StormAssert(edgeIndex < NumEdges());
    return EdgeEdges_.begin_row(edgeIndex);
  })
  ConstOverload(EdgeIndex*, BeginAdjacentEdge, (FaceIndex faceIndex), {
    StormAssert(faceIndex < NumFaces());
    return FaceEdges_.begin_row(faceIndex);
  })
  ConstOverload(EdgeIndex*, BeginAdjacentEdge, (CellIndex cellIndex), {
    StormAssert(cellIndex < NumCells());
    return CellEdges_.begin_row(cellIndex);
  })
  /// @}

  /// @brief Pointer to the beginning of the element adjacent faces.
  /// @{
  ConstOverload(FaceIndex*, BeginAdjacentFace, (NodeIndex nodeIndex), {
    StormAssert(nodeIndex < NumNodes());
    return NodeFaces_.begin_row(nodeIndex);
  })
  ConstOverload(FaceIndex*, BeginAdjacentFace, (EdgeIndex edgeIndex), {
    StormAssert(edgeIndex < NumEdges());
    return EdgeFaces_.begin_row(edgeIndex);
  })
  ConstOverload(FaceIndex*, BeginAdjacentFace, (FaceIndex faceIndex), {
    StormAssert(faceIndex < NumFaces());
    return FaceFaces_.begin_row(faceIndex);
  })
  ConstOverload(FaceIndex*, BeginAdjacentFace, (CellIndex cellIndex), {
    StormAssert(cellIndex < NumCells());
    return CellFaces_.begin_row(cellIndex);
  })
  /// @}

  /// @brief Pointer to the beginning of the element adjacent cells.
  /// @{
  ConstOverload(CellIndex*, BeginAdjacentCell, (NodeIndex nodeIndex), {
    StormAssert(nodeIndex < NumNodes());
    return NodeCells_.begin_row(nodeIndex);
  })
  ConstOverload(CellIndex*, BeginAdjacentCell, (EdgeIndex edgeIndex), {
    StormAssert(edgeIndex < NumEdges());
    return EdgeCells_.begin_row(edgeIndex);
  })
  ConstOverload(CellIndex*, BeginAdjacentCell, (FaceIndex faceIndex), {
    StormAssert(faceIndex < NumFaces());
    return FaceCells_.begin_row(faceIndex);
  })
  ConstOverload(CellIndex*, BeginAdjacentCell, (CellIndex cellIndex), {
    StormAssert(cellIndex < NumCells());
    return CellCells_.begin_row(cellIndex);
  })

  /// @brief Pointer to the End of the element adjacent nodes.
  /// @{
  ConstOverload(NodeIndex*, EndAdjacentNode, (NodeIndex nodeIndex), {
    StormAssert(nodeIndex < NumNodes());
    return NodeNodes_.end_row(nodeIndex);
  })
  ConstOverload(NodeIndex*, EndAdjacentNode, (EdgeIndex edgeIndex), {
    StormAssert(edgeIndex < NumEdges());
    return EdgeNodes_.end_row(edgeIndex);
  })
  ConstOverload(NodeIndex*, EndAdjacentNode, (FaceIndex faceIndex), {
    StormAssert(faceIndex < NumFaces());
    return FaceNodes_.end_row(faceIndex);
  })
  ConstOverload(NodeIndex*, EndAdjacentNode, (CellIndex cellIndex), {
    StormAssert(cellIndex < NumCells());
    return CellNodes_.end_row(cellIndex);
  })
  /// @}

  /// @brief Pointer to the End of the element adjacent edges.
  /// @{
  ConstOverload(EdgeIndex*, EndAdjacentEdge, (NodeIndex nodeIndex), {
    StormAssert(nodeIndex < NumNodes());
    return NodeEdges_.end_row(nodeIndex);
  })
  ConstOverload(EdgeIndex*, EndAdjacentEdge, (EdgeIndex edgeIndex), {
    StormAssert(edgeIndex < NumEdges());
    return EdgeEdges_.end_row(edgeIndex);
  })
  ConstOverload(EdgeIndex*, EndAdjacentEdge, (FaceIndex faceIndex), {
    StormAssert(faceIndex < NumFaces());
    return FaceEdges_.end_row(faceIndex);
  })
  ConstOverload(EdgeIndex*, EndAdjacentEdge, (CellIndex cellIndex), {
    StormAssert(cellIndex < NumCells());
    return CellEdges_.end_row(cellIndex);
  })
  /// @}

  /// @brief Pointer to the End of the element adjacent faces.
  /// @{
  ConstOverload(FaceIndex*, EndAdjacentFace, (NodeIndex nodeIndex), {
    StormAssert(nodeIndex < NumNodes());
    return NodeFaces_.end_row(nodeIndex);
  })
  ConstOverload(FaceIndex*, EndAdjacentFace, (EdgeIndex edgeIndex), {
    StormAssert(edgeIndex < NumEdges());
    return EdgeFaces_.end_row(edgeIndex);
  })
  ConstOverload(FaceIndex*, EndAdjacentFace, (FaceIndex faceIndex), {
    StormAssert(faceIndex < NumFaces());
    return FaceFaces_.end_row(faceIndex);
  })
  ConstOverload(FaceIndex*, EndAdjacentFace, (CellIndex cellIndex), {
    StormAssert(cellIndex < NumCells());
    return CellFaces_.end_row(cellIndex);
  })
  /// @}

  /// @brief Pointer to the End of the element adjacent cells.
  /// @{
  ConstOverload(CellIndex*, EndAdjacentCell, (NodeIndex nodeIndex), {
    StormAssert(nodeIndex < NumNodes());
    return NodeCells_.end_row(nodeIndex);
  })
  ConstOverload(CellIndex*, EndAdjacentCell, (EdgeIndex edgeIndex), {
    StormAssert(edgeIndex < NumEdges());
    return EdgeCells_.end_row(edgeIndex);
  })
  ConstOverload(CellIndex*, EndAdjacentCell, (FaceIndex faceIndex), {
    StormAssert(faceIndex < NumFaces());
    return FaceCells_.end_row(faceIndex);
  })
  ConstOverload(CellIndex*, EndAdjacentCell, (CellIndex cellIndex), {
    StormAssert(cellIndex < NumCells());
    return CellCells_.end_row(cellIndex);
  })
  /// @}

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
  void FixPermutationAndAdjacency_(std::vector<uint_t>& index1);

public:

  /// @brief Change order of all nodes.
  void PermuteNodes(std::vector<uint_t>&& nodePermutation = {});

  /// @brief Change order of all edges.
  void PermuteEdges(std::vector<uint_t>&& edgePermutation = {});

  /// @brief Change order of all faces.
  void PermuteFaces(std::vector<uint_t>&& facePermutation = {});

  /// @brief Change order of all cells.
  void PermuteCells(std::vector<uint_t>&& cellPermutation = {});

protected:

  void reorder_faces();

  // ---------------------------------------------------------------------- //
  // ---------------------------------------------------------------------- //

  /// @brief Generate edges using the face to node connectivity.
  /// @warning This function may be slow and memory-consuming.
  void generate_edges();

  /// @brief Generate faces using the cell to node connectivity.
  /// @warning This function may be slow and memory-consuming.
  void generate_faces();

  /// @brief Generate boundary cells to complete face connectivity.
  void generate_boundary_cells();

}; // class Mesh

} // namespace feathers

// Include iterators.
#include "MeshIterators.hh"

/// @brief @todo Remove me. 
using cMesh = feathers::Mesh;
#include "Field.hh"
