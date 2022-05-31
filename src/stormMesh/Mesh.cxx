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

#include <stormMesh/Mesh.hxx>
#include <stormMesh/Element.hh>
#include <stormUtils/Parallel.hh>
#include <stormUtils/Permute.hh>

namespace Storm {

void Mesh::UpdateElementsGeometry() {

  // Compute edge lengths and directions.
  std::tie(min_edge_len_, max_edge_len_) =
    for_each_min_max(edges(), +huge, -huge, [this](EdgeIndex edge_index) {
      std::unique_ptr<Shape> const edge_shape = shape(edge_index);
      edge_lens_[edge_index] = edge_shape->Volume();
      edge_dirs_[edge_index] = edge_shape->Dir();
      return edge_lens_[edge_index];
    });

  // Compute face areas, normals and center positions.
  std::tie(min_face_area_, max_face_area_) =
    for_each_min_max(faces(), +huge, -huge, [&](FaceIndex face_index) {
      std::unique_ptr<Shape> const face_shape = shape(face_index);
      face_areas_[face_index] = face_shape->Volume();
      face_normals_[face_index] = face_shape->Normal();
      face_barycenters_[face_index] = face_shape->CenterPos();
      return face_areas_[face_index];
    });

  // Compute cell volumes and center positions.
  std::tie(min_cell_volume_, max_cell_volume_) =
    for_each_min_max(cells(), +huge, -huge, [&](CellIndex cell_index) {
      std::unique_ptr<Shape> const cell_shape = shape(cell_index);
      cell_volumes_[cell_index] = cell_shape->Volume();
      cell_barycenters_[cell_index] = cell_shape->CenterPos();
      return cell_volumes_[cell_index];
    });

} // Mesh::UpdateElementsGeometry

NodeIndex Mesh::insert_node(vec3_t const& coords, NodeMark node_mark) {

  NodeIndex const node_index{num_nodes_++};

  // Emplace node properties.
  node_marks_.emplace_back(node_mark);
  node_coords_.emplace_back(coords);

  // Emplace empty rows that would be filled later on.
  node_nodes_.insert_row();
  node_edges_.insert_row();
  node_faces_.insert_row();
  node_cells_.insert_row();

  return node_index;

} // Mesh::insert_node

EdgeIndex Mesh::insert_edge(std::unique_ptr<Element>&& edge_shape,
                            EdgeMark edge_mark) {

  // Try to find an edge_shape first.
  std::set<NodeIndex> const edge_key(
    edge_shape->NodeIndices().begin(), edge_shape->NodeIndices().end());
  if (edge_lookup_.contains(edge_key)) {
    return edge_lookup_.at(edge_key);
  }

  EdgeIndex const edge_index{num_edges_++};
  edge_lookup_[edge_key] = edge_index;

  // Emplace edge_shape properties.
  edge_marks_.emplace_back(edge_mark);
  edge_shapes_.emplace_back(edge_shape->Shape());
  edge_lens_.emplace_back(edge_shape->Volume());
  edge_dirs_.emplace_back(edge_shape->Dir());

  // Fill the edge_shape-nodes and node-edges, edge-edges and node-nodes.
  edge_nodes_.insert_row(edge_shape->NodeIndices());
  for (NodeIndex node_index : edge_nodes_[edge_index]) {
    node_edges_.insert(node_index, edge_index);
  }

  // Emplace empty rows that would be filled later on.
  edge_faces_.insert_row();
  edge_cells_.insert_row();

  /// @todo Fill me!
  edge_edges_.insert_row();
  node_nodes_.insert_row();

  return edge_index;

} // Mesh::insert_edge

FaceIndex Mesh::insert_face(std::unique_ptr<Element>&& face_shape,
                            FaceMark face_mark) {

  // Try to find an edge first.
  std::set<NodeIndex> const face_key(
    face_shape->NodeIndices().begin(), face_shape->NodeIndices().end());
  if (face_lookup_.contains(face_key)) {
    return face_lookup_.at(face_key);
  }

  FaceIndex const face_index{num_faces_++};
  face_lookup_[face_key] = face_index;

  // Emplace the face properties.
  face_marks_.emplace_back(face_mark);
  face_shapes_.emplace_back(face_shape->Shape());
  face_areas_.emplace_back(face_shape->Volume());
  face_normals_.emplace_back(face_shape->Normal());
  face_barycenters_.emplace_back(face_shape->CenterPos());

  // Fill the face-nodes and node-faces.
  face_nodes_.insert_row(face_shape->NodeIndices());
  for (NodeIndex node_index : face_nodes_[face_index]) {
    node_faces_.insert(node_index, face_index);
  }

  // Fill the face-edges, edge-faces and face-faces.
  face_edges_.insert_row([&]() {
    boost::container::small_vector<EdgeIndex, 4> this_face_edges{};
    for (ShapeDesc& edge_desc : face_shape->make_edges_desc()) {
      EdgeIndex const edge_index{insert_edge(std::move(edge_desc))};
      this_face_edges.emplace_back(edge_index);
      edge_faces_.insert(edge_index, face_index);
    }
    return this_face_edges;
  }());

  // Emplace empty rows that would be filled later on.
  face_faces_.insert_row();
  face_cells_.insert_row();

  return face_index;

} // Mesh::insert_face

CellIndex Mesh::insert_cell(std::unique_ptr<Element>&& cell, CellMark cell_mark, bool ghost) {

  CellIndex const cell_index{num_cells_++};

  // Emplace cell properties.
  cell_marks_.emplace_back(cell_mark);
  cell_shapes_.emplace_back(cell->Shape());
  cell_volumes_.emplace_back(cell->Volume());
  cell_barycenters_.emplace_back(cell->CenterPos());

  // Fill the cell nodes and node cells.
  cell_nodes_.insert_row(cell->NodeIndices());
  for (NodeIndex node_index : cell_nodes_[cell_index]) {
    node_cells_.insert(node_index, cell_index);
  }

  if (ghost) {
    cell_edges_.insert_row();
    cell_faces_.insert_row();
    cell_cells_.insert_row();
    return cell_index;
  }

  // Fill the cell edges and edge cells.
  cell_edges_.insert_row([&]() {
    boost::container::small_vector<EdgeIndex, 12> this_cell_edges{};
    for (ShapeDesc& edge_desc : cell->make_edges_desc()) {
      EdgeIndex const edge_index{insert_edge(std::move(edge_desc))};
      this_cell_edges.emplace_back(edge_index);
      edge_cells_.insert(edge_index, cell_index);
    }
    return this_cell_edges;
  }());

  // Fill the cell faces and face cells.
  cell_faces_.insert_row([&]() {
    boost::container::small_vector<FaceIndex, 6> this_cell_faces{};
    for (ShapeDesc& face_desc : cell->make_faces_desc()) {
      FaceIndex const face_index{insert_face(ShapeDesc(face_desc))};
      this_cell_faces.emplace_back(face_index);
      if (face_cells_[face_index].empty() &&
          !ranges::equal(adjacent_nodes(face_index), face_desc.NodeIndices)) {
        flip_face(face_index);
      }
      face_cells_.insert(face_index, cell_index);
      if (mark(face_index) != FaceMark{0}) {
        face_cells_.insert(face_index, CellIndex{npos});
      }
    }
    return this_cell_faces;
  }());

  /// @todo Fill me!
  cell_cells_.insert_row();

  return cell_index;

} // Mesh::insert_cell

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

void Mesh::flip_face(FaceIndex face_index) noexcept {
  storm_assert(face_index < num_faces_ && "face_index is out of range");
  storm_assert([&]() {
    auto face_cells{face_cells_[face_index]};
    return face_cells.size() != 1 &&
      ranges::all_of(face_cells, [this](CellIndex cell_index) {
        return mark(cell_index) == CellMark{0};
      });
    }() && "face at face_index can not be flipped");

  ranges::reverse(face_nodes_[face_index]);
  ranges::reverse(face_edges_[face_index]);
  ranges::reverse(face_cells_[face_index]);
  face_normals_[face_index] = -face_normals_[face_index];

} // Mesh::flip_face

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
                  node_marks_.begin(), node_coords_.begin());

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
                  edge_marks_.begin(), edge_shapes_.begin(),
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
    face_marks_.begin(), face_shapes_.begin(),
    face_areas_.begin(), face_normals_.begin(), face_barycenters_.begin());

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
    cell_marks_.begin(), cell_shapes_.begin(),
    cell_volumes_.begin(), cell_barycenters_.begin());

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
  { ShapeType::Segment, { {1, 0}, {1, 0} } },
  /* 3D faces. */
  { ShapeType::Triangle, { {0, 2, 1}, {0, 2, 1} } },
  { ShapeType::Quadrangle, { {0, 3, 2, 1}, {0, 3, 2, 1} } },
};

/**
 * Generate boundary cells to complete face connectivity.
 */
void Mesh::generate_boundary_cells(FaceIndex ff) {

  auto const f = [&](MutableFaceView face) {
    if (face.mark() == 0) {
      return;
    }

    /* Boundary faces should be oriented outwards from the mesh. */
    if (ff == npos && face.inner_cell() == npos) {
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
    std::vector<NodeIndex> ghost_cell_nodes;
#if 0
    ghost_cell_nodes.assign(face.begin_node(), face.end_node());
#endif
#if 1
    cell.for_each_node([&](NodeView node) {
      if (std::find(face.adjacent_nodes().begin(), face.adjacent_nodes().end(), node) == face.adjacent_nodes().end()) {
        /* Reflect an interior cell Node. */
        // TODO: face normals are not computed here!
        // TODO: https://glm.g-truc.net/0.9.5/api/a00157.html#gab63646fc36b81cf69d3ce123a72f76f2
        vec3_t node_coords = node.coords();
        const vec3_t delta = node_coords - face.barycenter();
        node_coords -= 2.0 * glm::dot(delta, face.normal()) * face.normal();
        ghost_cell_nodes.push_back(
          insert_node(node_coords, (NodeMark) face.mark()));
      } else {
        /* Insert a boundary face Node. */
        ghost_cell_nodes.push_back(node);
      }
    });
#endif
    /* Insert the boundary cell3. */
    // TODO:
    const CellIndex boundaryCellIndex =
      insert_cell({cell.shapeType(), ghost_cell_nodes}, (CellMark) face.mark(), true);
    if (ff == npos) {
      //adjacent_faces(boundaryCellIndex)[0] = face;
      storm_assert(adjacent_cells(face)[face_outer_cell_] == npos);
      adjacent_cells(face)[face_outer_cell_] = boundaryCellIndex;
    } else {
      face_cells_.insert(face, boundaryCellIndex);
    }

  };

  if (ff != npos) {
    f(MutableFaceView(*this, ff));
    return;
  }
  ranges::for_each(face_views(*this), f);

} // Mesh::generate_boundary_cells

} // namespace feathers
