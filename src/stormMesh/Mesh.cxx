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

#include <boost/container/small_vector.hpp>

#include <stormMesh/Element.hh>
#include <stormMesh/Mesh.hxx>
#include <stormUtils/Parallel.hh>
#include <stormUtils/Permute.hh>

namespace Storm {

template<class T, size_t N>
using small_vector = boost::container::small_vector<T, N>;

void Mesh::UpdateElementsGeometry() {
  // Compute edge lengths and directions.
  std::tie(MinEdgeLen_, MaxEdgeLen_) =
      for_each_min_max(Edges(), +huge, -huge, [this](EdgeIndex edgeIndex) {
        std::unique_ptr<Shape> const edge_shape = shape(edgeIndex);
        EdgeLens_[edgeIndex] = edge_shape->Volume();
        EdgeDirs_[edgeIndex] = edge_shape->Dir();
        return EdgeLens_[edgeIndex];
      });

  // Compute face areas, normals and center positions.
  std::tie(MinFaceArea_, MaxFaceArea_) =
      for_each_min_max(Faces(), +huge, -huge, [&](FaceIndex faceIndex) {
        std::unique_ptr<Shape> const face_shape = shape(faceIndex);
        FaceAreas_[faceIndex] = face_shape->Volume();
        FaceNormals_[faceIndex] = face_shape->Normal();
        FaceCenters_[faceIndex] = face_shape->CenterPos();
        return FaceAreas_[faceIndex];
      });

  // Compute cell volumes and center positions.
  std::tie(MinCellVolume_, MaxCellVolume_) =
      for_each_min_max(Cells(), +huge, -huge, [&](CellIndex cellIndex) {
        std::unique_ptr<Shape> const cell_shape = shape(cellIndex);
        CellVolumes_[cellIndex] = cell_shape->Volume();
        CellCenters_[cellIndex] = cell_shape->CenterPos();
        return CellVolumes_[cellIndex];
      });

} // Mesh::UpdateElementsGeometry

NodeIndex Mesh::InsertNode(vec3_t const& coords, NodeMark nodeMark) {
  NodeIndex const nodeIndex{NumNodes_++};

  // Emplace node properties.
  NodeMarks_.emplace_back(nodeMark);
  NodeCoords_.emplace_back(coords);

  // Emplace empty rows that would be filled later on.
  NodeNodes_.insert_row();
  NodeEdges_.insert_row();
  NodeFaces_.insert_row();
  NodeCells_.insert_row();

  return nodeIndex;

} // Mesh::insert_node

EdgeIndex Mesh::InsertEdge(std::unique_ptr<Element>&& edgeShape,
                           EdgeMark edgeMark) {

  // Try to find an edge_shape first.
  std::set const edgeKey(edgeShape->NodeIndices().begin(),
                         edgeShape->NodeIndices().end());
  if (auto const it = EdgeLookup_.find(edgeKey); it != EdgeLookup_.end()) {
    return it->second;
  }

  EdgeIndex const edgeIndex{NumEdges_++};
  EdgeLookup_[edgeKey] = edgeIndex;

  // Emplace the edge properties.
  EdgeMarks_.emplace_back(edgeMark);
  EdgeShapes_.emplace_back(edgeShape->Shape());
  EdgeLens_.emplace_back(edgeShape->Volume());
  EdgeDirs_.emplace_back(edgeShape->Dir());

  // Fill the edge-nodes and node-edges, edge-edges and node-nodes.
  EdgeNodes_.insert_row(edgeShape->NodeIndices());
  for (NodeIndex nodeIndex : EdgeNodes_[edgeIndex]) {
    NodeEdges_.insert(nodeIndex, edgeIndex);
  }

  // Emplace empty rows that would be filled later on.
  EdgeEdges_.insert_row();
  EdgeFaces_.insert_row();
  EdgeCells_.insert_row();

  return edgeIndex;

} // Mesh::insert_edge

FaceIndex Mesh::InsertFace(std::unique_ptr<Element>&& faceShape,
                           FaceMark faceMark) {

#if 0
  auto const this_face_node_faces =
    face_shape->NodeIndices() |
    views::transform([this](NodeIndex nodeIndex) {
      small_vector<FaceIndex, 8> node_faces(
        adjacent_faces(nodeIndex).begin(), adjacent_faces(nodeIndex).end());
      ranges::sort(node_faces);
      return node_faces;
    });
  auto const found_faces = ranges::accumulate(
    this_face_node_faces | views::drop(1),
    this_face_node_faces.front(), [](auto&& first, auto&& second) {
      small_vector<FaceIndex, 8> out;
      ranges::set_intersection(first, second, std::back_inserter(out));
      return out;
    });
  storm_ensure(found_faces.size() <= 1 && "");
#endif

  // Try to find an edge first.
  std::set const faceKey(faceShape->NodeIndices().begin(),
                         faceShape->NodeIndices().end());
  if (auto const it = FaceLookup_.find(faceKey); it != FaceLookup_.end()) {
    return it->second;
  }

  FaceIndex const faceIndex{NumFaces_++};
  FaceLookup_[faceKey] = faceIndex;

  // Emplace the face properties.
  FaceMarks_.emplace_back(faceMark);
  FaceShapes_.emplace_back(faceShape->Shape());
  FaceAreas_.emplace_back(faceShape->Volume());
  FaceNormals_.emplace_back(faceShape->Normal());
  FaceCenters_.emplace_back(faceShape->CenterPos());

  // Fill the face-nodes and node-faces.
  FaceNodes_.insert_row(faceShape->NodeIndices());
  for (NodeIndex nodeIndex : FaceNodes_[faceIndex]) {
    NodeFaces_.insert(nodeIndex, faceIndex);
  }

  // Fill the face-edges, edge-faces and face-faces.
  FaceEdges_.insert_row([&]() {
    boost::container::small_vector<EdgeIndex, 4> thisFaceEdges{};
    for (ShapeDesc& edgeDesc : faceShape->make_edges_desc()) {
      EdgeIndex const edgeIndex{InsertEdge(std::move(edgeDesc))};
      thisFaceEdges.emplace_back(edgeIndex);
      EdgeFaces_.insert(edgeIndex, faceIndex);
    }
    return thisFaceEdges;
  }());

  // Emplace empty rows that would be filled later on.
  FaceFaces_.insert_row();
  FaceCells_.insert_row();

  return faceIndex;

} // Mesh::insert_face

CellIndex Mesh::insert_cell(std::unique_ptr<Element>&& cell, CellMark cellMark,
                            bool ghost) {
  CellIndex const cellIndex{NumCells_++};

  // Emplace cell properties.
  CellMarks_.emplace_back(cellMark);
  CellShapes_.emplace_back(cell->Shape());
  CellVolumes_.emplace_back(cell->Volume());
  CellCenters_.emplace_back(cell->CenterPos());

  // Fill the cell Nodes and node Cells.
  CellNodes_.insert_row(cell->NodeIndices());
  for (NodeIndex nodeIndex : CellNodes_[cellIndex]) {
    NodeCells_.insert(nodeIndex, cellIndex);
  }

  if (ghost) {
    CellEdges_.insert_row();
    CellFaces_.insert_row();
    CellCells_.insert_row();
    return cellIndex;
  }

  // Fill the cell Edges and edge Cells.
  CellEdges_.insert_row([&]() {
    boost::container::small_vector<EdgeIndex, 12> thisCellEdges{};
    for (ShapeDesc& edgeDesc : cell->make_edges_desc()) {
      EdgeIndex const edgeIndex{InsertEdge(std::move(edgeDesc))};
      thisCellEdges.emplace_back(edgeIndex);
      EdgeCells_.insert(edgeIndex, cellIndex);
    }
    return thisCellEdges;
  }());

  // Fill the cell Faces and face Cells.
  CellFaces_.insert_row([&]() {
    boost::container::small_vector<FaceIndex, 6> thisCellFaces{};
    for (ShapeDesc& faceDesc : cell->make_faces_desc()) {
      FaceIndex const faceIndex{InsertFace(ShapeDesc(faceDesc))};
      thisCellFaces.emplace_back(faceIndex);
      if (bool const faceOrientedAsInner =
              ranges::equal(AdjacentNodes(faceIndex), faceDesc.NodeIndices);
          FaceCells_[faceIndex].empty() && !faceOrientedAsInner) {
        flip_face(faceIndex);
      } else {
        StormAssert(faceOrientedAsInner);
      }
      FaceCells_.insert(faceIndex, cellIndex);
    }
    return thisCellFaces;
  }());

  /// @todo Fill me!
  CellCells_.insert_row();

  return cellIndex;

} // Mesh::insert_cell

// ------------------------------------------------------------------------------------
// //
// ------------------------------------------------------------------------------------
// //

void Mesh::flip_face(FaceIndex faceIndex) noexcept {
  StormAssert(faceIndex < NumFaces_ && "faceIndex is out of range");
  StormAssert([&]() {
    auto face_cells{FaceCells_[faceIndex]};
    return face_cells.size() != 1 &&
           ranges::all_of(face_cells, [this](CellIndex cellIndex) {
             return Mark(cellIndex) == CellMark{0};
           });
  }() && "face at faceIndex can not be flipped");

  ranges::reverse(FaceNodes_[faceIndex]);
  ranges::reverse(FaceEdges_[faceIndex]);
  ranges::reverse(FaceCells_[faceIndex]);
  FaceNormals_[faceIndex] = -FaceNormals_[faceIndex];

} // Mesh::flip_face

template<class Tag>
void Mesh::FixPermutationAndAdjacency_(std::vector<size_t>& permutation) {
  /* Fix permutations by resorting it by marks.
   * Stable sort is used here to preserve the
   * permutation in the best possible way. */
  std::stable_sort(permutation.begin(), permutation.end(),
                   [&](auto index1, auto index2) {
                     return Mark(Index<Tag>(index1)) < Mark(Index<Tag>(index2));
                   });

  /* Fix adjacency tables. */
  std::vector<size_t> inversePermutation(permutation.size(), npos);
  InversePermutation(permutation.begin(), permutation.end(),
                     inversePermutation.begin());

  auto const reorderFunc = [&inversePermutation](Index<Tag> index) {
    return Index<Tag>(index != npos ? inversePermutation[size_t(index)] : npos);
  };

  ranges::for_each(Nodes(), [&](NodeIndex nodeIndex) {
    ranges::transform(AdjacentElements_<Tag>(nodeIndex),
                      std::begin(AdjacentElements_<Tag>(nodeIndex)),
                      reorderFunc);
  });

  ranges::for_each(Edges(), [&](EdgeIndex edgeIndex) {
    ranges::transform(AdjacentElements_<Tag>(edgeIndex),
                      std::begin(AdjacentElements_<Tag>(edgeIndex)),
                      reorderFunc);
  });

  ranges::for_each(Faces(), [&](FaceIndex faceIndex) {
    ranges::transform(AdjacentElements_<Tag>(faceIndex),
                      std::begin(AdjacentElements_<Tag>(faceIndex)),
                      reorderFunc);
  });

  ranges::for_each(Cells(), [&](CellIndex cellIndex) {
    ranges::transform(AdjacentElements_<Tag>(cellIndex),
                      std::begin(AdjacentElements_<Tag>(cellIndex)),
                      reorderFunc);
  });

} // Mesh::FixPermutationAndAdjacency_

void Mesh::PermuteNodes(std::vector<size_t>&& nodePermutation) {
  /* Permute Node properties and fix the adjacency tables. */
  FixPermutationAndAdjacency_<NodeTag>(nodePermutation);
  permute_rows(nodePermutation.begin(), nodePermutation.end(), NodeNodes_,
               NodeEdges_, NodeFaces_, NodeCells_);
  permute_inplace(nodePermutation.begin(), nodePermutation.end(),
                  NodeMarks_.begin(), NodeCoords_.begin());

  /* Generate the Node ranges. */
  NodeRanges_.clear();
  ranges::for_each(NodeViews(*this), [&](NodeView node) {
    NodeRanges_.resize((size_t) node.mark() + 2);
    NodeRanges_[node.mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&) NodeRanges_;
  std::partial_sum(r.begin(), r.end(), r.begin());

} // Mesh::PermuteNodes

void Mesh::PermuteEdges(std::vector<size_t>&& edgePermutation) {
  /* Permute edge properties and fix the adjacency tables. */
  FixPermutationAndAdjacency_<EdgeTag>(edgePermutation);
  permute_rows(edgePermutation.begin(), edgePermutation.end(), EdgeNodes_,
               EdgeEdges_, EdgeFaces_, EdgeCells_);
  permute_inplace(edgePermutation.begin(), edgePermutation.end(),
                  EdgeMarks_.begin(), EdgeShapes_.begin(), EdgeLens_.begin(),
                  EdgeDirs_.begin());

  /* Generate the edge ranges. */
  EdgeRanges_.clear();
  ranges::for_each(EdgeViews(*this), [&](EdgeView edge) {
    EdgeRanges_.resize((size_t) edge.mark() + 2);
    EdgeRanges_[edge.mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&) EdgeRanges_;
  std::partial_sum(r.begin(), r.end(), r.begin());

} // Mesh::PermuteEdges

void Mesh::PermuteFaces(std::vector<size_t>&& facePermutation) {
  /* Permute data. */
  FixPermutationAndAdjacency_<FaceTag>(facePermutation);
  permute_rows(facePermutation.begin(), facePermutation.end(), FaceNodes_,
               FaceEdges_, FaceFaces_, FaceCells_);
  permute_inplace(facePermutation.begin(), facePermutation.end(),
                  FaceMarks_.begin(), FaceShapes_.begin(), FaceAreas_.begin(),
                  FaceNormals_.begin(), FaceCenters_.begin());

  /* Generate Mark ranges. */
  FaceRanges_.clear();
  ranges::for_each(FaceViews(*this), [&](FaceView face) {
    FaceRanges_.resize((size_t) face.mark() + 2);
    FaceRanges_[face.mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&) FaceRanges_;
  std::partial_sum(r.begin(), r.end(), r.begin());

} // Mesh::PermuteFaces

void Mesh::PermuteCells(std::vector<size_t>&& cellPermutation) {
  /* Permute data. */
  FixPermutationAndAdjacency_<CellTag>(cellPermutation);
  permute_rows(cellPermutation.begin(), cellPermutation.end(), CellNodes_,
               CellEdges_, CellFaces_, CellCells_);
  permute_inplace(cellPermutation.begin(), cellPermutation.end(),
                  CellMarks_.begin(), CellShapes_.begin(), CellVolumes_.begin(),
                  CellCenters_.begin());

  /* Generate Mark ranges. */
  CellRanges_.clear();
  ranges::for_each(CellViews(*this), [&](CellView cell) {
    CellRanges_.resize((size_t) cell.mark() + 2);
    CellRanges_[cell.mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&) CellRanges_;
  std::partial_sum(r.begin(), r.end(), r.begin());

} // Mesh::PermuteCells

// ------------------------------------------------------------------------------------
// //
// ------------------------------------------------------------------------------------
// //

void Mesh::reorder_faces() {
  /** @todo Refactor me! */
  {
    std::vector<size_t> node_reordering(Nodes().size());
    std::iota(node_reordering.begin(), node_reordering.end(), 0);
    PermuteNodes(std::move(node_reordering));
  }
  {
    std::vector<size_t> edge_reordering(Edges().size());
    std::iota(edge_reordering.begin(), edge_reordering.end(), 0);
    PermuteEdges(std::move(edge_reordering));
  }
  {
    std::vector<size_t> face_reordering(Faces().size());
    std::iota(face_reordering.begin(), face_reordering.end(), 0);
    PermuteFaces(std::move(face_reordering));
  }
  {
    std::vector<size_t> cell_reordering(Cells().size());
    std::iota(cell_reordering.begin(), cell_reordering.end(), 0);
    PermuteCells(std::move(cell_reordering));
  }

} // Mesh::PermuteFaces

// ------------------------------------------------------------------------------------
// //
// ------------------------------------------------------------------------------------
// //

/* A Node and edge flip table for various face types. */
static const std::map<ShapeType,
                      std::pair<std::vector<size_t>, std::vector<size_t>>>
    g_face_shape_to_nodes_and_edges_flip{
        /* 1D Faces. */
        {ShapeType::Node, {{0}, {0}}},
        /* 2D Faces. */
        {ShapeType::Segment, {{1, 0}, {1, 0}}},
        /* 3D Faces. */
        {ShapeType::Triangle, {{0, 2, 1}, {0, 2, 1}}},
        {ShapeType::Quadrangle, {{0, 3, 2, 1}, {0, 3, 2, 1}}},
    };

/**
 * Generate boundary Cells to complete face connectivity.
 */
void Mesh::generate_boundary_cells(FaceIndex ff) {
  auto const f = [&](MutableFaceView face) {
    if (face.mark() == 0) { return; }

    CellView cell = face.InnerCell();

    /* Generate the boundary cell: reflect a connected interior cell. */
    std::vector<NodeIndex> ghost_cell_nodes;
    cell.ForEachNode([&](NodeView node) {
      if (std::find(face.AdjacentNodes().begin(), face.AdjacentNodes().end(),
                    node) == face.AdjacentNodes().end()) {
        /* Reflect an interior cell Node. */
        // TODO: face normals are not computed here!
        // TODO:
        // https://glm.g-truc.net/0.9.5/api/a00157.html#gab63646fc36b81cf69d3ce123a72f76f2
        vec3_t node_coords = node.Coords();
        const vec3_t delta = node_coords - face.Center();
        node_coords -= 2.0 * glm::dot(delta, face.Normal()) * face.Normal();
        ghost_cell_nodes.push_back(
            InsertNode(node_coords, (NodeMark) face.mark()));
      } else {
        /* Insert a boundary face Node. */
        ghost_cell_nodes.push_back(node);
      }
    });

    /* Insert the boundary cell3. */
    const CellIndex bndr_cellIndex = insert_cell(
        {cell.shapeType(), ghost_cell_nodes}, (CellMark) face.mark(), true);
    FaceCells_.insert(face, bndr_cellIndex);
  };

  if (ff != npos) {
    f(MutableFaceView(*this, ff));
    return;
  }
  ranges::for_each(FaceViews(*this), f);

} // Mesh::generate_boundary_cells

} // namespace Storm
