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
  std::tie(min_edge_len_, max_edge_len_) =
    for_each_min_max(edges(), +huge, -huge, [this](EdgeIndex edge_index) {
      std::unique_ptr<Shape> const edgeShape = shape(edge_index);
      edge_lens_[edge_index] = edgeShape->Volume();
      edge_dirs_[edge_index] = edgeShape->Dir();
      return edge_lens_[edge_index];
    });

  // Compute face areas, normals and center positions.
  std::tie(min_face_area_, max_face_area_) =
    for_each_min_max(faces(), +huge, -huge, [&](FaceIndex face_index) {
      std::unique_ptr<Shape> faceShape = shape(face_index);
      face_areas_[face_index] = faceShape->Volume();
      face_normals_[face_index] = faceShape->Normal();
      FaceCenterPos_[face_index] = faceShape->CenterPos();
      return face_areas_[face_index];
    });

  // Compute cell volumes and center positions.
  std::tie(min_cell_volume_, max_cell_volume_) =
    for_each_min_max(cells(), +huge, -huge, [&](CellIndex cell_index) {
      std::unique_ptr<Shape> const cellShape = shape(cell_index);
      cell_volumes_[cell_index] = cellShape->Volume();
      CellCenterPos_[cell_index] = cellShape->CenterPos();
      return cell_volumes_[cell_index];
    });

} // Mesh::UpdateElementsGeometry

NodeIndex Mesh::insertNode(vec3_t const& nodePos, NodeMark nodeMark) {

  NodeIndex const node_index(num_nodes_++);

  // Emplace node properties.
  node_marks_.emplace_back(nodeMark);
  NodePos_.emplace_back(nodePos);

  // These should be filled later.
  node_nodes_.emplaceRow();
  node_edges_.emplaceRow();
  node_faces_.emplaceRow();
  node_cells_.emplaceRow();

  return node_index;

} // Mesh::insertNode

EdgeIndex Mesh::insertEdge(std::unique_ptr<Element>&& edge, EdgeMark edgeMark) {

  // Try to find an edge first.
  std::set<size_t> const edgeKey(
    edge->NodeIndices().begin(), edge->NodeIndices().end());
  if (EdgeLookup_.contains(edgeKey)) {
    return EdgeLookup_.at(edgeKey);
  }

  EdgeIndex const edge_index(num_edges_++);
  EdgeLookup_[edgeKey] = edge_index;

  // Emplace edge properties.
  edge_marks_.emplace_back(edgeMark);
  EdgeShapeTypes_.emplace_back(edge->Shape());
  edge_lens_.emplace_back(edge->Volume());
  edge_dirs_.emplace_back(edge->Dir());

  // Fill the edge nodes and node edges.
  edge_nodes_.emplaceRow(edge->NodeIndices() |
    views::transform([&](size_t node_index_) {
      NodeIndex const node_index{node_index_};
      node_edges_.insert(node_index, edge_index);
      return node_index;
    }));

  // These should be filled later.
  edge_edges_.emplaceRow();
  edge_faces_.emplaceRow();
  edge_cells_.emplaceRow();

  /// @todo Fill me!
  edge_edges_.emplaceRow();
  node_nodes_.emplaceRow();

  return edge_index;

} // Mesh::insertEdge

FaceIndex Mesh::insertFace(std::unique_ptr<Element>&& face, FaceMark faceMark) {

  // Try to find an edge first.
  std::set<size_t> const faceKey(
    face->NodeIndices().begin(), face->NodeIndices().end());
  if (FaceLookup_.contains(faceKey)) {
    return FaceLookup_.at(faceKey);
  }

  FaceIndex const face_index(num_faces_++);
  FaceLookup_[faceKey] = face_index;

  // Emplace face properties.
  face_marks_.emplace_back(faceMark);
  FaceShapeTypes_.emplace_back(face->Shape());
  face_areas_.emplace_back(face->Volume());
  face_normals_.emplace_back(face->Normal());
  FaceCenterPos_.emplace_back(face->CenterPos());

  // Fill the face nodes and node faces.
  face_nodes_.emplaceRow(face->NodeIndices() |
    views::transform([&](size_t node_index_) {
      NodeIndex const node_index{node_index_};
      node_faces_.insert(node_index, face_index);
      return node_index;
    }));

  // Fill the face edges and edge faces.
  //std::vector<EdgeIndex> faceEdges;
  //ranges::transform(face->MakeEdgesDesc(),
  //  std::back_inserter(faceEdges), [&](ShapeDesc edgeDesc) {
  //    EdgeIndex const edge_index = insertEdge(std::move(edgeDesc));
  //    edge_faces_.insert(edge_index, face_index);
  //    return edge_index;
  //  });
  //face_edges_.emplaceRow(faceEdges.begin(), faceEdges.end());
  face_edges_.emplaceRow();

  // This should be filled later.
  face_cells_.emplaceRow(2); // @todo here should be no 2!

  /// @todo Fill me!
  face_faces_.emplaceRow();

  return face_index;

} // Mesh::insertFace

CellIndex Mesh::insertCell(std::unique_ptr<Element>&& cell, CellMark cellMark, bool ghost) {

  CellIndex const cell_index(num_cells_++);

  // Emplace cell properties.
  cell_marks_.emplace_back(cellMark);
  CellShapeTypes_.emplace_back(cell->Shape());
  cell_volumes_.emplace_back(cell->Volume());
  CellCenterPos_.emplace_back(cell->CenterPos());

  // Fill the cell nodes and node cells.
  cell_nodes_.emplaceRow(cell->NodeIndices() |
    views::transform([&](size_t node_index_) {
      NodeIndex const node_index{node_index_};
      node_cells_.insert(node_index, cell_index);
      return node_index;
    }));

  if (ghost) {
    cell_edges_.emplaceRow();
    cell_faces_.emplaceRow();
    cell_edges_.emplaceRow();
    return cell_index;
  }

  // Fill the cell edges and edge cells.
  //std::vector<EdgeIndex> cellEdges;
  //ranges::transform(cell->MakeEdgesDesc(),
  //  std::back_inserter(cellEdges), [&](ShapeDesc edgeDesc) {
  //    EdgeIndex const edge_index = insertEdge(std::move(edgeDesc));
  //    edge_cells_.insert(edge_index, cell_index);
  //    return edge_index;
  //  });
  //cell_edges_.emplaceRow(cellEdges.begin(), cellEdges.end());
  cell_edges_.emplaceRow();

  // Fill the cell faces and face cells.
  std::vector<FaceIndex> cellFaces;
  ranges::transform(cell->MakeFacesDesc(),
    std::back_inserter(cellFaces), [&](ShapeDesc faceDesc) {
      FaceIndex const face_index = EmplaceFace(std::move(ShapeDesc(faceDesc)));
      if (std::equal(adjacent_nodes(face_index).begin(), adjacent_nodes(face_index).end(), faceDesc.NodeIndices.begin())) {
        face_cells_[face_index][face_inner_cell_] = cell_index;
      } else {
        face_cells_[face_index][face_outer_cell_] = cell_index;
      }
      return face_index;
    });
  cell_faces_.emplaceRow(cellFaces.begin(), cellFaces.end());

  /// @todo Fill me!
  cell_cells_.emplaceRow();

  return cell_index;

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

  ranges::for_each(nodes(), [&](NodeIndex node_index) {
    ranges::transform(
      adjacent_elements<Tag>(node_index),
      std::begin(adjacent_elements<Tag>(node_index)), reorderFunc);
  });

  ranges::for_each(edges(), [&](EdgeIndex edge_index) {
    ranges::transform(
      adjacent_elements<Tag>(edge_index),
      std::begin(adjacent_elements<Tag>(edge_index)), reorderFunc);
  });

  ranges::for_each(faces(), [&](FaceIndex face_index) {
    ranges::transform(
      adjacent_elements<Tag>(face_index),
      std::begin(adjacent_elements<Tag>(face_index)), reorderFunc);
  });

  ranges::for_each(cells(), [&](CellIndex cell_index) {
    ranges::transform(
      adjacent_elements<Tag>(cell_index),
      std::begin(adjacent_elements<Tag>(cell_index)), reorderFunc);
  });

} // Mesh::FixPermutationAndAdjacency_

void Mesh::PermuteNodes(std::vector<size_t>&& nodePermutation) {

  /* Permute Node properties and fix the adjacency tables. */
  FixPermutationAndAdjacency_<NodeTag>(nodePermutation);
  permute_rows(nodePermutation.begin(), nodePermutation.end(),
               node_nodes_, node_edges_, node_faces_, node_cells_);
  permute_inplace(nodePermutation.begin(), nodePermutation.end(),
                  node_marks_.begin(), NodePos_.begin());

  /* Generate the Node ranges. */
  node_ranges_.clear();
  ranges::for_each(node_views(*this), [&](NodeView node) {
    node_ranges_.resize((size_t) node.mark() + 2);
    node_ranges_[node.mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&)node_ranges_;
  std::partial_sum(r.begin(), r.end(), r.begin());

} // Mesh::PermuteNodes

void Mesh::PermuteEdges(std::vector<size_t>&& edgePermutation) {

  /* Permute edge properties and fix the adjacency tables. */
  FixPermutationAndAdjacency_<EdgeTag>(edgePermutation);
  permute_rows(edgePermutation.begin(), edgePermutation.end(),
               edge_nodes_, edge_edges_, edge_faces_, edge_cells_);
  permute_inplace(edgePermutation.begin(), edgePermutation.end(),
                  edge_marks_.begin(), EdgeShapeTypes_.begin(),
                  edge_lens_.begin(), edge_dirs_.begin());

  /* Generate the edge ranges. */
  edge_ranges_.clear();
  ranges::for_each(edge_views(*this), [&](EdgeView edge) {
    edge_ranges_.resize((size_t) edge.mark() + 2);
    edge_ranges_[edge.mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&)edge_ranges_;
  std::partial_sum(r.begin(), r.end(), r.begin());

} // Mesh::PermuteEdges

void Mesh::PermuteFaces(std::vector<size_t>&& facePermutation) {

  /* Permute data. */
  FixPermutationAndAdjacency_<FaceTag>(facePermutation);
  permute_rows(
    facePermutation.begin(), facePermutation.end(),
    face_nodes_, face_edges_, face_faces_, face_cells_);
  permute_inplace(
    facePermutation.begin(), facePermutation.end(),
    face_marks_.begin(), FaceShapeTypes_.begin(),
    face_areas_.begin(), face_normals_.begin(), FaceCenterPos_.begin());

  /* Generate mark ranges. */
  face_ranges_.clear();
  ranges::for_each(face_views(*this), [&](FaceView face) {
    face_ranges_.resize((size_t) face.mark() + 2);
    face_ranges_[face.mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&)face_ranges_;
  std::partial_sum(r.begin(), r.end(), r.begin());

} // Mesh::PermuteFaces

void Mesh::PermuteCells(std::vector<size_t>&& cellPermutation) {

  /* Permute data. */
  FixPermutationAndAdjacency_<CellTag>(cellPermutation);
  permute_rows(
    cellPermutation.begin(), cellPermutation.end(),
    cell_nodes_, cell_edges_, cell_faces_, cell_cells_);
  permute_inplace(
    cellPermutation.begin(), cellPermutation.end(),
    cell_marks_.begin(), CellShapeTypes_.begin(),
    cell_volumes_.begin(), CellCenterPos_.begin());

  /* Generate mark ranges. */
  cell_ranges_.clear();
  ranges::for_each(cell_views(*this), [&](CellView cell) {
    cell_ranges_.resize((size_t) cell.mark() + 2);
    cell_ranges_[cell.mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&)cell_ranges_;
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

  ranges::for_each(face_views(*this), [&](MutableFaceView face) {
    if (face.mark() == 0) {
      return;
    }

    /* Boundary faces should be oriented outwards from the mesh. */
    if (face.inner_cell() == npos) {
      /* Flip normal and cell connectivity. */
      face_normals_[face] = -face_normals_[face];
      std::swap(adjacent_cells(face)[face_inner_cell_],
                adjacent_cells(face)[face_outer_cell_]);
      /* Flip Node and edge connectivity. */
      std::vector<size_t> node_permutation;
      std::vector<size_t> edge_permutation;
      std::tie(node_permutation, edge_permutation) =
        g_face_shape_to_nodes_and_edges_flip.at(face.shapeType());
      permute_inplace(
        node_permutation.begin(), node_permutation.end(), std::begin(adjacent_nodes(face)));
      permute_inplace(
        edge_permutation.begin(), edge_permutation.end(), std::begin(adjacent_nodes(face)));
    }
    CellView cell = face.inner_cell();

    /* Generate the boundary cell: reflect a connected interior cell. */
    std::vector<size_t> ghost_cell_nodes;
#if 0
    ghost_cell_nodes.assign(face.begin_node(), face.end_node());
#endif
#if 1
    cell.for_each_node([&](NodeView node) {
      if (ranges::find(adjacent_nodes(face), (NodeIndex) node) == std::end(adjacent_nodes(face))) {
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
    adjacent_faces(boundaryCellIndex)[0] = face;
#if 0
    Cell& boundary_cell = get_cell(boundaryCellIndex);
        while (boundary_cell._num_faces() != 1) {
            boundary_cell.erase_face(0);
        }
        boundary_cell.begin_face()[0] = face;
#endif

    storm_assert(adjacent_cells(face)[face_outer_cell_] == npos);
    adjacent_cells(face)[face_outer_cell_] = boundaryCellIndex;

  });
} // Mesh::generate_boundary_cells

} // namespace feathers
