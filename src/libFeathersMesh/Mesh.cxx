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

#include <libFeathersMesh/Mesh.hxx>
#include <libFeathersMesh/Element.hh>
#include <libFeathersUtils/Parallel.hh>
#include <libFeathersUtils/Permute.hh>

#include <set>
#include <map>
#include <ranges>

namespace feathers {

void Mesh::UpdateElementsGeometry() {

  // Compute edge lengths and directions.
  std::tie(MinEdgeLen_, MaxEdgeLen_) =
    ForEachMinMax(edges(), +huge, -huge, [this](EdgeIndex edgeIndex) {
      std::unique_ptr<Shape> const edgeShape = shape(edgeIndex);
      EdgeLens_[edgeIndex] = edgeShape->Volume();
      EdgeDirs_[edgeIndex] = edgeShape->Dir();
      return EdgeLens_[edgeIndex];
    });

  // Compute face areas, normals and center positions.
  std::tie(MinFaceArea_, MaxFaceArea_) =
    ForEachMinMax(faces(), +huge, -huge, [&](FaceIndex faceIndex) {
      std::unique_ptr<Shape> faceShape = shape(faceIndex);
      FaceAreas_[faceIndex] = faceShape->Volume();
      FaceNormals_[faceIndex] = faceShape->Normal();
      FaceCenterPos_[faceIndex] = faceShape->CenterPos();
      return FaceAreas_[faceIndex];
    });

  // Compute cell volumes and center positions.
  std::tie(MinFaceArea_, MaxFaceArea_) =
    ForEachMinMax(cells(), +huge, -huge, [&](CellIndex cellIndex) {
      std::unique_ptr<Shape> const cellShape = shape(cellIndex);
      CellVolumes_[cellIndex] = cellShape->Volume();
      CellCenterPos_[cellIndex] = cellShape->CenterPos();
      return CellVolumes_[cellIndex];
    });

} // Mesh::UpdateElementsGeometry

NodeIndex Mesh::insertNode(vec3_t const& nodePos, NodeMark nodeMark) {

  NodeIndex const nodeIndex(NumNodes_++);

  // Emplace node properties.
  NodeMarks_.emplace_back(nodeMark);
  NodePos_.emplace_back(nodePos);

  // These should be filled later.
  NodeNodes_.emplaceRow();
  NodeEdges_.emplaceRow();
  NodeFaces_.emplaceRow();
  NodeCells_.emplaceRow();

  return nodeIndex;

} // Mesh::insertNode

EdgeIndex Mesh::insertEdge(std::unique_ptr<Element>&& edge, EdgeMark edgeMark) {

  // Try to find an edge first.
  std::set<size_t> const edgeKey(
    edge->NodeIndices().begin(), edge->NodeIndices().end());
  if (EdgeLookup_.contains(edgeKey)) {
    return EdgeLookup_.at(edgeKey);
  }

  EdgeIndex const edgeIndex(NumEdges_++);
  EdgeLookup_[edgeKey] = edgeIndex;

  // Emplace edge properties.
  EdgeMarks_.emplace_back(edgeMark);
  EdgeShapeTypes_.emplace_back(edge->Shape());
  EdgeLens_.emplace_back(edge->Volume());
  EdgeDirs_.emplace_back(edge->Dir());

  // Fill the edge nodes and node edges.
  EdgeNodes_.emplaceRow(edge->NodeIndices() |
    views::transform([&](size_t nodeIndex_) {
      NodeIndex const nodeIndex{nodeIndex_};
      NodeEdges_.insert(nodeIndex, edgeIndex);
      return nodeIndex;
    }));

  // These should be filled later.
  EdgeEdges_.emplaceRow();
  EdgeFaces_.emplaceRow();
  EdgeCells_.emplaceRow();

  /// @todo Fill me!
  EdgeEdges_.emplaceRow();
  NodeNodes_.emplaceRow();

  return edgeIndex;

} // Mesh::insertEdge

FaceIndex Mesh::insertFace(std::unique_ptr<Element>&& face, FaceMark faceMark) {

  // Try to find an edge first.
  std::set<size_t> const faceKey(
    face->NodeIndices().begin(), face->NodeIndices().end());
  if (FaceLookup_.contains(faceKey)) {
    return FaceLookup_.at(faceKey);
  }

  FaceIndex const faceIndex(NumFaces_++);
  FaceLookup_[faceKey] = faceIndex;

  // Emplace face properties.
  FaceMarks_.emplace_back(faceMark);
  FaceShapeTypes_.emplace_back(face->Shape());
  FaceAreas_.emplace_back(face->Volume());
  FaceNormals_.emplace_back(face->Normal());
  FaceCenterPos_.emplace_back(face->CenterPos());

  // Fill the face nodes and node faces.
  FaceNodes_.emplaceRow(face->NodeIndices() |
    views::transform([&](size_t nodeIndex_) {
      NodeIndex const nodeIndex{nodeIndex_};
      NodeFaces_.insert(nodeIndex, faceIndex);
      return nodeIndex;
    }));

  // Fill the face edges and edge faces.
  //std::vector<EdgeIndex> faceEdges;
  //ranges::transform(face->MakeEdgesDesc(),
  //  std::back_inserter(faceEdges), [&](ShapeDesc edgeDesc) {
  //    EdgeIndex const edgeIndex = insertEdge(std::move(edgeDesc));
  //    EdgeFaces_.insert(edgeIndex, faceIndex);
  //    return edgeIndex;
  //  });
  //FaceEdges_.emplaceRow(faceEdges.begin(), faceEdges.end());
  FaceEdges_.emplaceRow();

  // This should be filled later.
  FaceCells_.emplaceRow(2); // @todo here should be no 2!

  /// @todo Fill me!
  FaceFaces_.emplaceRow();

  return faceIndex;

} // Mesh::insertFace

CellIndex Mesh::insertCell(std::unique_ptr<Element>&& cell, CellMark cellMark, bool ghost) {

  CellIndex const cellIndex(NumCells_++);

  // Emplace cell properties.
  CellMarks_.emplace_back(cellMark);
  CellShapeTypes_.emplace_back(cell->Shape());
  CellVolumes_.emplace_back(cell->Volume());
  CellCenterPos_.emplace_back(cell->CenterPos());

  // Fill the cell nodes and node cells.
  CellNodes_.emplaceRow(cell->NodeIndices() |
    views::transform([&](size_t nodeIndex_) {
      NodeIndex const nodeIndex{nodeIndex_};
      NodeCells_.insert(nodeIndex, cellIndex);
      return nodeIndex;
    }));

  if (ghost) {
    CellEdges_.emplaceRow();
    CellFaces_.emplaceRow();
    CellEdges_.emplaceRow();
    return cellIndex;
  }

  // Fill the cell edges and edge cells.
  //std::vector<EdgeIndex> cellEdges;
  //ranges::transform(cell->MakeEdgesDesc(),
  //  std::back_inserter(cellEdges), [&](ShapeDesc edgeDesc) {
  //    EdgeIndex const edgeIndex = insertEdge(std::move(edgeDesc));
  //    EdgeCells_.insert(edgeIndex, cellIndex);
  //    return edgeIndex;
  //  });
  //CellEdges_.emplaceRow(cellEdges.begin(), cellEdges.end());
  CellEdges_.emplaceRow();

  // Fill the cell faces and face cells.
  std::vector<FaceIndex> cellFaces;
  ranges::transform(cell->MakeFacesDesc(),
    std::back_inserter(cellFaces), [&](ShapeDesc faceDesc) {
      FaceIndex const faceIndex = EmplaceFace(std::move(ShapeDesc(faceDesc)));
      if (std::equal(adjNodes(faceIndex).begin(), adjNodes(faceIndex).end(), faceDesc.NodeIndices.begin())) {
        FaceCells_[faceIndex][FaceInnerCell_] = cellIndex;
      } else {
        FaceCells_[faceIndex][FaceOuterCell_] = cellIndex;
      }
      return faceIndex;
    });
  CellFaces_.emplaceRow(cellFaces.begin(), cellFaces.end());

  /// @todo Fill me!
  CellCells_.emplaceRow();

  return cellIndex;

} // Mesh::insertCell

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

template<class Tag>
void Mesh::FixPermutationAndAdjacency_(std::vector<size_t>& permutation) {

  /* Fix permutations by resorting it by marks.
   * Stable sort is used here to preserve the
   * permutation in the best possible way. */
  std::stable_sort(permutation.begin(), permutation.end(),
    [&](auto index1, auto index2) {
      return mark(Index<Tag>(index1)) < mark(Index<Tag>(index2));
    });

  /* Fix adjacency tables. */
  std::vector<size_t> inversePermutation(permutation.size(), npos);
  InversePermutation(permutation.begin(), permutation.end(), inversePermutation.begin());

  auto const reorderFunc = [&inversePermutation](Index<Tag> index) {
    return Index<Tag>(index != npos ? inversePermutation[size_t(index)] : npos);
  };

  ranges::for_each(nodes(), [&](NodeIndex nodeIndex) {
    ranges::transform(
      AdjacentElements_<Tag>(nodeIndex),
      std::begin(AdjacentElements_<Tag>(nodeIndex)), reorderFunc);
  });

  ranges::for_each(edges(), [&](EdgeIndex edgeIndex) {
    ranges::transform(
      AdjacentElements_<Tag>(edgeIndex),
      std::begin(AdjacentElements_<Tag>(edgeIndex)), reorderFunc);
  });

  ranges::for_each(faces(), [&](FaceIndex faceIndex) {
    ranges::transform(
      AdjacentElements_<Tag>(faceIndex),
      std::begin(AdjacentElements_<Tag>(faceIndex)), reorderFunc);
  });

  ranges::for_each(cells(), [&](CellIndex cellIndex) {
    ranges::transform(
      AdjacentElements_<Tag>(cellIndex),
      std::begin(AdjacentElements_<Tag>(cellIndex)), reorderFunc);
  });

} // Mesh::FixPermutationAndAdjacency_

void Mesh::PermuteNodes(std::vector<size_t>&& nodePermutation) {

  /* Permute Node properties and fix the adjacency tables. */
  FixPermutationAndAdjacency_<NodeTag>(nodePermutation);
  permute_rows(nodePermutation.begin(), nodePermutation.end(),
               NodeNodes_, NodeEdges_, NodeFaces_, NodeCells_);
  permute_inplace(nodePermutation.begin(), nodePermutation.end(),
                  NodeMarks_.begin(), NodePos_.begin());

  /* Generate the Node ranges. */
  NodeRanges_.clear();
  ranges::for_each(nodeViews(*this), [&](NodeView node) {
    NodeRanges_.resize((size_t) node.mark() + 2);
    NodeRanges_[node.mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&)NodeRanges_;
  std::partial_sum(r.begin(), r.end(), r.begin());

} // Mesh::PermuteNodes

void Mesh::PermuteEdges(std::vector<size_t>&& edgePermutation) {

  /* Permute edge properties and fix the adjacency tables. */
  FixPermutationAndAdjacency_<EdgeTag>(edgePermutation);
  permute_rows(edgePermutation.begin(), edgePermutation.end(),
               EdgeNodes_, EdgeEdges_, EdgeFaces_, EdgeCells_);
  permute_inplace(edgePermutation.begin(), edgePermutation.end(),
                  EdgeMarks_.begin(), EdgeShapeTypes_.begin(),
                  EdgeLens_.begin(), EdgeDirs_.begin());

  /* Generate the edge ranges. */
  EdgeRanges_.clear();
  ranges::for_each(edgeViews(*this), [&](EdgeView edge) {
    EdgeRanges_.resize((size_t) edge.mark() + 2);
    EdgeRanges_[edge.mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&)EdgeRanges_;
  std::partial_sum(r.begin(), r.end(), r.begin());

} // Mesh::PermuteEdges

void Mesh::PermuteFaces(std::vector<size_t>&& facePermutation) {

  /* Permute data. */
  FixPermutationAndAdjacency_<FaceTag>(facePermutation);
  permute_rows(
    facePermutation.begin(), facePermutation.end(),
    FaceNodes_, FaceEdges_, FaceFaces_, FaceCells_);
  permute_inplace(
    facePermutation.begin(), facePermutation.end(),
    FaceMarks_.begin(), FaceShapeTypes_.begin(),
    FaceAreas_.begin(), FaceNormals_.begin(), FaceCenterPos_.begin());

  /* Generate mark ranges. */
  FaceRanges_.clear();
  ranges::for_each(faceViews(*this), [&](FaceView face) {
    FaceRanges_.resize((size_t) face.mark() + 2);
    FaceRanges_[face.mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&)FaceRanges_;
  std::partial_sum(r.begin(), r.end(), r.begin());

} // Mesh::PermuteFaces

void Mesh::PermuteCells(std::vector<size_t>&& cellPermutation) {

  /* Permute data. */
  FixPermutationAndAdjacency_<CellTag>(cellPermutation);
  permute_rows(
    cellPermutation.begin(), cellPermutation.end(),
    CellNodes_, CellEdges_, CellFaces_, CellCells_);
  permute_inplace(
    cellPermutation.begin(), cellPermutation.end(),
    CellMarks_.begin(), CellShapeTypes_.begin(),
    CellVolumes_.begin(), CellCenterPos_.begin());

  /* Generate mark ranges. */
  CellRanges_.clear();
  ranges::for_each(cellViews(*this), [&](CellView cell) {
    CellRanges_.resize((size_t) cell.mark() + 2);
    CellRanges_[cell.mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&)CellRanges_;
  std::partial_sum(r.begin(), r.end(), r.begin());

} // Mesh::PermuteCells

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

void Mesh::reorder_faces() {

  /** @todo Refactor me! */
  {
    std::vector<size_t> node_reordering(nodes().size());
    std::iota(node_reordering.begin(), node_reordering.end(), 0);
    PermuteNodes(std::move(node_reordering));
  }
  {
    std::vector<size_t> edge_reordering(edges().size());
    std::iota(edge_reordering.begin(), edge_reordering.end(), 0);
    PermuteEdges(std::move(edge_reordering));
  }
  {
    std::vector<size_t> face_reordering(faces().size());
    std::iota(face_reordering.begin(), face_reordering.end(), 0);
    PermuteFaces(std::move(face_reordering));
  }
  {
    std::vector<size_t> cell_reordering(cells().size());
    std::iota(cell_reordering.begin(), cell_reordering.end(), 0);
    PermuteCells(std::move(cell_reordering));
  }

} // Mesh::PermuteFaces

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/* A Node and edge flip table for various face types. */
static const std::map<ShapeType, std::pair<std::vector<size_t>, std::vector<size_t>>>
  g_face_shape_to_nodes_and_edges_flip {
  /* 1D faces. */
  { ShapeType::Node, { {0}, {0} } },
  /* 2D faces. */
  { ShapeType::Segment2, { {1, 0}, {1, 0} } },
  /* 3D faces. */
  { ShapeType::Triangle3, { {0, 2, 1}, {0, 2, 1} } },
  { ShapeType::Quadrangle4, { {0, 3, 2, 1}, {0, 3, 2, 1} } },
};

/**
 * Generate boundary cells to complete face connectivity.
 */
void Mesh::generate_boundary_cells() {

  ranges::for_each(faceViews(*this), [&](MutableFaceView face) {
    if (face.mark() == 0) {
      return;
    }

    /* Boundary faces should be oriented outwards from the mesh. */
    if (face.innerCell() == npos) {
      /* Flip normal and cell connectivity. */
      FaceNormals_[face] = -FaceNormals_[face];
      std::swap(adjCells(face)[FaceInnerCell_],
                adjCells(face)[FaceOuterCell_]);
      /* Flip Node and edge connectivity. */
      std::vector<size_t> node_permutation;
      std::vector<size_t> edge_permutation;
      std::tie(node_permutation, edge_permutation) =
        g_face_shape_to_nodes_and_edges_flip.at(face.shapeType());
      permute_inplace(
        node_permutation.begin(), node_permutation.end(), std::begin(adjNodes(face)));
      permute_inplace(
        edge_permutation.begin(), edge_permutation.end(), std::begin(adjNodes(face)));
    }
    CellView cell = face.innerCell();

    /* Generate the boundary cell: reflect a connected interior cell. */
    std::vector<size_t> ghost_cell_nodes;
#if 0
    ghost_cell_nodes.assign(face.begin_node(), face.end_node());
#endif
#if 1
    cell.forEachNode([&](NodeView node) {
      if (ranges::find(adjNodes(face), (NodeIndex) node) == std::end(adjNodes(face))) {
        /* Reflect an interior cell Node. */
        // TODO: face normals are not computed here!
        // TODO: https://glm.g-truc.net/0.9.5/api/a00157.html#gab63646fc36b81cf69d3ce123a72f76f2
        vec3_t node_coords = node.pos();
        const vec3_t delta = node_coords - face.centerPos();
        node_coords -= 2.0 * glm::dot(delta, face.normal()) * face.normal();
        ghost_cell_nodes.push_back(
          (size_t) insertNode(node_coords, (NodeMark) face.mark()));
      } else {
        /* Insert a boundary face Node. */
        ghost_cell_nodes.push_back((size_t) node);
      }
    });
#endif
    /* Insert the boundary cell. */
    // TODO:
    const CellIndex boundaryCellIndex =
      EmplaceCell({cell.shapeType(), ghost_cell_nodes}, (CellMark) face.mark(), true);
    adjFaces(boundaryCellIndex)[0] = face;
#if 0
    Cell& boundary_cell = get_cell(boundaryCellIndex);
        while (boundary_cell._num_faces() != 1) {
            boundary_cell.erase_face(0);
        }
        boundary_cell.begin_face()[0] = face;
#endif

    StormAssert(adjCells(face)[FaceOuterCell_] == npos);
    adjCells(face)[FaceOuterCell_] = boundaryCellIndex;

  });
} // Mesh::generate_boundary_cells

} // namespace feathers
