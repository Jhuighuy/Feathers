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
#include <set>

#include "Element.hh"
#include "Field.hh"
#include <stormMesh/Base.hxx>
#include <stormUtils/Image.hh>
#include <stormUtils/Table.hh>

namespace Storm {

inline constexpr size_t FaceInnerCell_{0};
inline constexpr size_t FaceOuterCell_{1};

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Hybrid unstructured multidimensional mesh.
/// @todo On our way to stateless mesh:
/// @todo 1. Polish the insertion behaviour: marks assignments and others.
/// @todo 2. Implement the symmetric topoligies insertion.
/// @todo 3. Switch to the faster lookup with set algorithms.
/// @todo 3. ShapeType + ShapeDesc vs. Shape. Ghosts?
/// @todo Do not use iterators in the mesh implementation.
/// @todo "Element" -> "Shape".
/// @todo Make Dim a template parameter.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
// template<size_t SDim, size_t TDim, template<class, class, class> class Table>
class Mesh : public tObject<Mesh> {
private:
  size_t NumNodes_{0}, NumEdges_{0};
  size_t NumFaces_{0}, NumCells_{0};

  Vector<NodeMark, NodeIndex> NodeMarks_;
  Vector<EdgeMark, EdgeIndex> EdgeMarks_;
  Vector<FaceMark, FaceIndex> FaceMarks_;
  Vector<CellMark, CellIndex> CellMarks_;
  Vector<NodeIndex, NodeMark> NodeRanges_;
  Vector<EdgeIndex, EdgeMark> EdgeRanges_;
  Vector<FaceIndex, FaceMark> FaceRanges_;
  Vector<CellIndex, CellMark> CellRanges_;

  Vector<ShapeType, EdgeIndex> EdgeShapes_;
  Vector<ShapeType, FaceIndex> FaceShapes_;
  Vector<ShapeType, CellIndex> CellShapes_;
  Vector<vec3_t, NodeIndex> NodeCoords_;
  Vector<real_t, EdgeIndex> EdgeLens_;
  Vector<vec3_t, EdgeIndex> EdgeDirs_;
  Vector<real_t, FaceIndex> FaceAreas_;
  Vector<vec3_t, FaceIndex> FaceNormals_;
  Vector<vec3_t, FaceIndex> FaceCenters_;
  Vector<real_t, CellIndex> CellVolumes_;
  Vector<vec3_t, CellIndex> CellCenters_;
  real_t MinEdgeLen_{qnan}, MaxEdgeLen_{qnan};
  real_t MinFaceArea_{qnan}, MaxFaceArea_{qnan};
  real_t MinCellVolume_{qnan}, MaxCellVolume_{qnan};

  CsrTable<EdgeIndex, NodeIndex> EdgeNodes_;
  CsrTable<FaceIndex, NodeIndex> FaceNodes_;
  CsrTable<CellIndex, NodeIndex> CellNodes_;
  CsrTable<FaceIndex, EdgeIndex> FaceEdges_;
  CsrTable<CellIndex, EdgeIndex> CellEdges_;
  CsrTable<CellIndex, FaceIndex> CellFaces_;

  CsrTable<NodeIndex, EdgeIndex> NodeEdges_;
  CsrTable<NodeIndex, FaceIndex> NodeFaces_;
  CsrTable<NodeIndex, CellIndex> NodeCells_;
  CsrTable<EdgeIndex, FaceIndex> EdgeFaces_;
  CsrTable<EdgeIndex, CellIndex> EdgeCells_;
  CsrTable<FaceIndex, CellIndex> FaceCells_;

  CsrTable<NodeIndex, NodeIndex> NodeNodes_;
  CsrTable<EdgeIndex, EdgeIndex> EdgeEdges_;
  CsrTable<FaceIndex, FaceIndex> FaceFaces_;
  CsrTable<CellIndex, CellIndex> CellCells_;

  std::map<std::set<NodeIndex>, EdgeIndex> EdgeLookup_;
  std::map<std::set<NodeIndex>, FaceIndex> FaceLookup_;
  std::map<std::set<NodeIndex>, CellIndex> CellLookup_;

public:
  /// ---------------------------------------------------------------- ///
  /// @name Constructors.
  /// ---------------------------------------------------------------- ///
  /// @{

  /// @brief Initialize an empty mesh.
  Mesh() noexcept = default;

  bool read_from_triangle(std::string const& path);

  bool read_from_tetgen(std::string const& path);

  bool read_from_image(const char* path,
                       const std::map<Pixel, size_t>& mark_colors,
                       Pixel fluid_color = eBlackPixel,
                       vec2_t pixel_size = vec2_t(1.0, 1.0));

  void save_vtk(const char* path, const std::vector<sFieldDesc>& fields) const;

  void finalize() {
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
  auto Nodes() const noexcept {
    return views::iota(NodeIndex{0}, NodeIndex{NumNodes_});
  }

  /// @brief Range of edge indices.
  auto Edges() const noexcept {
    return views::iota(EdgeIndex{0}, EdgeIndex{NumEdges_});
  }

  /// @brief Range of face indices.
  auto Faces() const noexcept {
    return views::iota(FaceIndex{0}, FaceIndex{NumFaces_});
  }

  /// @brief Range of cell indices.
  auto Cells() const noexcept {
    return views::iota(CellIndex{0}, CellIndex{NumCells_});
  }

  /// @}

  /// ---------------------------------------------------------------- ///
  /// @name Marks.
  /// ---------------------------------------------------------------- ///
  /// @{

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

  /// @brief Range of node indices with a @p nodeMark.
  auto Nodes(NodeMark nodeMark) const noexcept {
    StormAssert(nodeMark < NumNodeMarks());
    return views::iota(NodeRanges_[nodeMark], NodeRanges_[nodeMark + 1]);
  }

  /// @brief Range of edge indices with a @p edgeMark.
  auto Edges(EdgeMark edgeMark) const noexcept {
    StormAssert(edgeMark < NumEdgeMarks());
    return views::iota(EdgeRanges_[edgeMark], EdgeRanges_[edgeMark + 1]);
  }

  /// @brief Range of face indices with a @p faceMark.
  auto Faces(FaceMark faceMark) const noexcept {
    StormAssert(faceMark < NumFaceMarks());
    return views::iota(FaceRanges_[faceMark], FaceRanges_[faceMark + 1]);
  }

  /// @brief Range of cell indices with a @p cellMark.
  auto Cells(CellMark cellMark) const noexcept {
    StormAssert(cellMark < NumCellMarks());
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

  /// @}

  /// ---------------------------------------------------------------- ///
  /// @name Shapes.
  /// ---------------------------------------------------------------- ///
  /// @{

  /// @brief Get edge @p edgeIndex shape type.
  ShapeType shapeType(EdgeIndex edgeIndex) const noexcept {
    StormAssert(edgeIndex < NumEdges_);
    return EdgeShapes_[edgeIndex];
  }

  /// @brief Get face @p faceIndex shape type.
  ShapeType shapeType(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces_);
    return FaceShapes_[faceIndex];
  }

  /// @brief Get cell @p cellIndex shape type.
  ShapeType shapeType(CellIndex cellIndex) const noexcept {
    StormAssert(cellIndex < NumCells_);
    return CellShapes_[cellIndex];
  }

private:
  auto makeShape_(ShapeDesc&& desc) const {
    return Element::Make(std::forward<ShapeDesc>(desc), NodeCoords_);
  }

public:
  /// @brief Get element object.
  template<class Tag>
  auto shape(Index<Tag> index) const {
    auto const node_indices = AdjacentNodes(index);
    return makeShape_({shapeType(index),
                       std::vector(node_indices.begin(), node_indices.end())});
  }

  /// @brief Get node @p nodeIndex coordinates.
  vec3_t NodeCoords(NodeIndex nodeIndex) const noexcept {
    StormAssert(nodeIndex < NumNodes_ && "nodeIndex is out of range");
    return NodeCoords_[nodeIndex];
  }

  /// @brief Set node @p nodeIndex position @p Coords.
  void SetNodeCoords(NodeIndex nodeIndex, vec3_t const& coords) noexcept {
    StormAssert(nodeIndex < NumNodes_ && "nodeIndex is out of range");
    NodeCoords_[nodeIndex] = coords;
  }

  /// @brief Get edge @p edgeIndex length.
  real_t EdgeLen(EdgeIndex edgeIndex) const noexcept {
    StormAssert(edgeIndex < NumEdges_ && "edgeIndex is out of range");
    return EdgeLens_[edgeIndex];
  }

  /// @brief Get edge @p edgeIndex direction.
  vec3_t EdgeDir(EdgeIndex edgeIndex) const noexcept {
    StormAssert(edgeIndex < NumEdges_ && "edgeIndex is out of range");
    return EdgeDirs_[edgeIndex];
  }

  /// @brief Get face @p faceIndex Area/length.
  real_t FaceArea(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces_ && "faceIndex is out of range");
    return FaceAreas_[faceIndex];
  }

  /// @brief Get face @p faceIndex Normal.
  vec3_t FaceNormal(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces_ && "faceIndex is out of range");
    return FaceNormals_[faceIndex];
  }

  /// @brief Get face @p faceIndex Center.
  vec3_t FaceCenter(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces_ && "faceIndex is out of range");
    return FaceCenters_[faceIndex];
  }

  /// @brief Get cell @p cellIndex Volume/Area/length.
  real_t CellVolume(CellIndex cellIndex) const noexcept {
    StormAssert(cellIndex < NumCells_ && "cellIndex is out of range");
    return CellVolumes_[cellIndex];
  }

  /// @brief Get cell @p cellIndex Center.
  vec3_t CellCenter(CellIndex cellIndex) const noexcept {
    StormAssert(cellIndex < NumCells_ && "cellIndex is out of range");
    return CellCenters_[cellIndex];
  }

  /// @brief Get the minimal edge length.
  real_t MinEdgeLen() const noexcept {
    StormAssert(std::isfinite(MinEdgeLen_));
    return MinEdgeLen_;
  }

  /// @brief Get the maximal edge length.
  real_t MaxEdgeLen() const noexcept {
    StormAssert(std::isfinite(MaxEdgeLen_));
    return MaxEdgeLen_;
  }

  /// @brief Get the minimal face Area.
  real_t MinFaceArea() const noexcept {
    StormAssert(std::isfinite(MinFaceArea_));
    return MinFaceArea_;
  }

  /// @brief Get the maximal face Area.
  real_t MaxFaceArea() const noexcept {
    StormAssert(std::isfinite(MaxFaceArea_));
    return MaxFaceArea_;
  }

  /// @brief Get the minimal cell Volume.
  real_t MinCellVolume() const noexcept {
    StormAssert(std::isfinite(MinCellVolume_));
    return MinCellVolume_;
  }

  /// @brief Get the maximal cell Volume.
  real_t MaxCellVolume() const noexcept {
    StormAssert(std::isfinite(MaxCellVolume_));
    return MaxCellVolume_;
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
  /// @name Adjacency. @todo remove "const_cast<Mesh*>(this)->"
  /// ---------------------------------------------------------------- ///
  /// @{

  /// @name Primary adjacency: adjacency,
  ///   that can be extracted from the shape only.
  /// @{

  /// @brief Range of the edge @p edgeIndex adjacent node indices.
  /// Denote a node to be adjacent to an edge if it one its nodes.
  auto AdjacentNodes(EdgeIndex edgeIndex) const noexcept {
    StormAssert(edgeIndex < NumEdges_ && "edgeIndex is out of range");
    return const_cast<Mesh*>(this)->EdgeNodes_[edgeIndex];
  }

  /// @brief Range of the face @p faceIndex adjacent node indices.
  /// Denote a node to be adjacent to a face if it one its nodes.
  auto AdjacentNodes(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces_ && "faceIndex is out of range");
    return const_cast<Mesh*>(this)->FaceNodes_[faceIndex];
  }

  /// @brief Range of the cell @p cellIndex adjacent node indices.
  /// Denote a node to be adjacent to a cell if it one its nodes.
  auto AdjacentNodes(CellIndex cellIndex) const noexcept {
    StormAssert(cellIndex < NumCells_ && "cellIndex is out of range");
    return const_cast<Mesh*>(this)->CellNodes_[cellIndex];
  }

  /// @brief Range of the face @p faceIndex adjacent edge indices.
  /// Denote an edge to be adjacent to a face if it one its edges.
  auto AdjacentEdges(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces_ && "faceIndex is out of range");
    return const_cast<Mesh*>(this)->FaceEdges_[faceIndex];
  }

  /// @brief Range of the cell @p cellIndex adjacent edge indices.
  /// Denote an edge to be adjacent to a cell if it one its edges.
  auto AdjacentEdges(CellIndex cellIndex) const noexcept {
    StormAssert(cellIndex < NumCells_ && "cellIndex is out of range");
    return const_cast<Mesh*>(this)->CellEdges_[cellIndex];
  }

  /// @brief Range of the cell @p cellIndex adjacent face indices.
  /// Denote a face to be adjacent to a cell if it one its faces.
  auto AdjacentFaces(CellIndex cellIndex) const noexcept {
    StormAssert(cellIndex < NumCells_ && "cellIndex is out of range");
    return const_cast<Mesh*>(this)->CellFaces_[cellIndex];
  }

  /// @}

  /// @name Secondary adjacency: adjacency,
  ///   that is a transpose of the primary adjacency.
  /// @{

  /// @brief Range of the node @p nodeIndex adjacent edge indices.
  auto AdjacentEdges(NodeIndex nodeIndex) const noexcept {
    StormAssert(nodeIndex < NumNodes_ && "nodeIndex is out of range");
    return const_cast<Mesh*>(this)->NodeEdges_[nodeIndex];
  }

  /// @brief Range of the node @p nodeIndex adjacent face indices.
  auto AdjacentFaces(NodeIndex nodeIndex) const noexcept {
    StormAssert(nodeIndex < NumNodes_ && "nodeIndex is out of range");
    return const_cast<Mesh*>(this)->NodeFaces_[nodeIndex];
  }

  /// @brief Range of the node @p nodeIndex adjacent cell indices.
  auto AdjacentCells(NodeIndex nodeIndex) const noexcept {
    StormAssert(nodeIndex < NumNodes_ && "nodeIndex is out of range");
    return const_cast<Mesh*>(this)->NodeCells_[nodeIndex];
  }

  /// @brief Range of the edge @p edgeIndex adjacent face indices.
  auto AdjacentFaces(EdgeIndex edgeIndex) const noexcept {
    StormAssert(edgeIndex < NumEdges_ && "edgeIndex is out of range");
    return const_cast<Mesh*>(this)->EdgeFaces_[edgeIndex];
  }

  /// @brief Range of the edge @p edgeIndex adjacent cell indices.
  auto AdjacentCells(EdgeIndex edgeIndex) const noexcept {
    StormAssert(edgeIndex < NumEdges_ && "edgeIndex is out of range");
    return const_cast<Mesh*>(this)->EdgeCells_[edgeIndex];
  }

  /// @brief Range of the face @p faceIndex adjacent cell indices.
  auto AdjacentCells(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces_ && "faceIndex is out of range");
    return const_cast<Mesh*>(this)->FaceCells_[faceIndex];
  }

  /// @}

  /// @name Symmetric adjacency: adjacency, that is a product
  ///   of the corresponding primary and secondary adjacencies.
  /// @{

  /// @brief Range of the node @p nodeIndex adjacent node indices.
  /// Denote two nodes as adjacent if there is an edge connecting them.
  auto AdjacentNodes(NodeIndex nodeIndex) const noexcept {
    StormAssert(nodeIndex < NumNodes_ && "nodeIndex is out of range");
    return const_cast<Mesh*>(this)->NodeNodes_[nodeIndex];
  }

  /// @brief Range of the edge @p edgeIndex adjacent edge indices.
  /// Denote two edges as adjacent if they share a common node.
  auto AdjacentEdges(EdgeIndex edgeIndex) const noexcept {
    StormAssert(edgeIndex < NumEdges_ && "edgeIndex is out of range");
    return const_cast<Mesh*>(this)->EdgeEdges_[edgeIndex];
  }

  /// @brief Range of the face @p faceIndex adjacent face indices.
  /// Denote two faces as adjacent if they share a common edge.
  auto AdjacentFaces(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumFaces_ && "faceIndex is out of range");
    return const_cast<Mesh*>(this)->FaceFaces_[faceIndex];
  }

  /// @brief Range of the cell @p cellIndex adjacent cell indices.
  /// Denote two cells as adjacent if they share a common face.
  auto AdjacentCells(CellIndex cellIndex) const noexcept {
    StormAssert(cellIndex < NumCells_ && "cellIndex is out of range");
    return const_cast<Mesh*>(this)->CellCells_[cellIndex];
  }

  /// @}

private:
  template<class Tag>
  auto AdjacentElements_(auto elementIndex) const noexcept {
    if constexpr (std::same_as<Tag, NodeTag>) {
      return AdjacentNodes(elementIndex);
    } else if constexpr (std::same_as<Tag, EdgeTag>) {
      return AdjacentEdges(elementIndex);
    } else if constexpr (std::same_as<Tag, FaceTag>) {
      return AdjacentFaces(elementIndex);
    } else if constexpr (std::same_as<Tag, CellTag>) {
      return AdjacentCells(elementIndex);
    } else {
      static_assert(always_false<Tag>, "Invalid tag.");
    }
  }

public:
  /// @}

  /// ---------------------------------------------------------------- ///
  /// @name Insertions.
  /// ---------------------------------------------------------------- ///
  /// @{

  /// @brief Insert a new node with a position @p Coords
  ///   and a Mark @p nodeMark into the mesh.
  /// @returns Index of the inserted node.
  NodeIndex InsertNode(vec3_t const& coords, NodeMark nodeMark = {});

  /// @brief Find or emplace a new edge_shape with a shape @p edgeShape
  ///   and node indices @p Nodes.
  ///
  /// Insertion would update the edge_shape-node, edge_shape-edge_shape and
  ///   node-node topologies.
  ///
  /// @param edgeMark Mark that would be assigned to an edge_shape if
  ///   the insertion took place, otherwise ignored.
  /// @returns A pair of an index of the found or inserted edge_shape
  ///   and a boolean value denoting whether the insertion took place.
  // std::pair<EdgeIndex, bool>
  //   findOrInsertEdge(std::shared_ptr<Shape> edgeShape,
  //                    EdgeMark edgeMark = {});
  EdgeIndex InsertEdge(std::unique_ptr<Element>&& edgeShape,
                       EdgeMark edgeMark = {});
  EdgeIndex InsertEdge(ShapeDesc&& edgeDesc, EdgeMark edgeMark = {}) {
    return InsertEdge(makeShape_(std::forward<ShapeDesc>(edgeDesc)), edgeMark);
  }

  /// @brief Find or emplace a new face_shape with a shape @p faceShape
  ///   and node indices @p Nodes.
  ///
  /// Insertion would implicitly insert the missing Edges and
  ///   update the face_shape-edge and face_shape-face_shape topologies. The
  ///   implicitly inserted Edges would inherit the Mark from
  ///   explicitly inserted face_shape.
  ///
  /// @param faceMark Mark that would be assigned to a face_shape if
  ///   the insertion took place, otherwise ignored.
  /// @returns A pair of an index of the found or inserted face_shape
  ///   and a boolean value denoting whether the insertion took place.
  // std::pair<FaceIndex, bool>
  //   findOrInsertFace(std::shared_ptr<Shape> faceShape,
  //                    FaceMark faceMark = {});
  FaceIndex InsertFace(std::unique_ptr<Element>&& faceShape,
                       FaceMark faceMark = {});
  FaceIndex InsertFace(ShapeDesc&& faceDesc, FaceMark faceMark = {}) {
    return InsertFace(makeShape_(std::forward<ShapeDesc>(faceDesc)), faceMark);
  }

  /// @brief Find or emplace a new cell with a shape @p cellShape
  ///   and node indices @p Nodes.
  ///
  /// Insertion would implicitly insert the missing Faces and
  ///   update the cell-face and cell-cell topologies. The
  ///   implicitly inserted face would inherit the Mark from
  ///   explicitly inserted cell.
  ///
  /// @param cellMark Mark that would be assigned to a cell if
  ///   the insertion took place, otherwise ignored.
  /// @returns A pair of an index of the found or inserted cell
  ///   and a boolean value denoting whether the insertion took place.
  // std::pair<CellIndex, bool>
  //   findOrInsertCell(std::shared_ptr<Shape> cellShape,
  //                    CellMark cellMark = {});
  CellIndex insert_cell(std::unique_ptr<Element>&& cell, CellMark cellMark = {},
                        bool ghost = false);
  CellIndex insert_cell(ShapeDesc&& cellDesc, CellMark cellMark = {},
                        bool ghost = false) {
    return insert_cell(makeShape_(std::forward<ShapeDesc>(cellDesc)), cellMark,
                       ghost);
  }

  /// @brief Generate boundary Cells to complete face connectivity.
  void generate_boundary_cells(FaceIndex ff = FaceIndex{npos});

  /// @}

  /// ---------------------------------------------------------------- ///
  /// @name Permutations.
  /// ---------------------------------------------------------------- ///
  /// @{

  /// @brief Flip the face @p faceIndex.
  /// A face can be flipped only if it is not adjacent to any Cells
  /// or it is adjacent to the exactly two interior Cells.
  void flip_face(FaceIndex faceIndex) noexcept;

  /// @brief Change order of all Nodes.
  void PermuteNodes(std::vector<size_t>&& nodePermutation = {});

  /// @brief Change order of all Edges.
  void PermuteEdges(std::vector<size_t>&& edgePermutation = {});

  /// @brief Change order of all Faces.
  void PermuteFaces(std::vector<size_t>&& facePermutation = {});

  /// @brief Change order of all Cells.
  void PermuteCells(std::vector<size_t>&& cellPermutation = {});

private:
  template<class Tag>
  void FixPermutationAndAdjacency_(std::vector<size_t>& cellIndex);

protected:
  void reorder_faces();

  /// @}

}; // class Mesh

} // namespace Storm

// Include iterators.
#include "View.hxx"

/// @brief @todo Remove me.
using cMesh = Storm::Mesh;
