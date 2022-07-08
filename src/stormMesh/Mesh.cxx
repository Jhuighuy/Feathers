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

#include <ranges>

#include <stormMesh/Element.hh>
#include <stormMesh/Mesh.hxx>
#include <stormUtils/Parallel.hh>
#include <stormUtils/Permute.hh>

namespace Storm {

template<class T, size_t N>
using small_vector = boost::container::small_vector<T, N>;

void Mesh::update_elements_geometry() {
  // Compute edge lengths and directions.
  std::tie(min_edge_len_, max_edge_len_) =
      for_each_min_max(edges(), +huge, -huge, [this](EdgeIndex edge_index) {
        const std::unique_ptr<Shape> edge_shape = shape(edge_index);
        edge_lens_[edge_index] = edge_shape->volume(node_coords_);
        edge_dirs_[edge_index] = edge_shape->dir(node_coords_);
        return edge_lens_[edge_index];
      });

  // Compute face areas, normals and center positions.
  std::tie(min_face_area_, max_face_area_) =
      for_each_min_max(faces(), +huge, -huge, [&](FaceIndex face_index) {
        const std::unique_ptr<Shape> face_shape = shape(face_index);
        face_areas_[face_index] = face_shape->volume(node_coords_);
        face_normals_[face_index] = face_shape->normal(node_coords_);
        face_centers_[face_index] = face_shape->center_pos(node_coords_);
        return face_areas_[face_index];
      });

  // Compute cell volumes and center positions.
  std::tie(min_cell_volume_, max_cell_volume_) =
      for_each_min_max(cells(), +huge, -huge, [&](CellIndex cell_index) {
        const std::unique_ptr<Shape> cell_shape = shape(cell_index);
        cell_volumes_[cell_index] = cell_shape->volume(node_coords_);
        cell_centers_[cell_index] = cell_shape->center_pos(node_coords_);
        return cell_volumes_[cell_index];
      });

} // Mesh::UpdateElementsGeometry

NodeIndex Mesh::insert_node(const vec3_t& coords, NodeMark node_mark) {
  const NodeIndex node_index{num_nodes_};

  // Emplace node properties.
  num_nodes_ += 1;
  node_marks_.emplace_back(node_mark);
  node_coords_.emplace_back(coords);

  // Emplace empty rows that would be filled later on.
  node_nodes_.insert_row();
  node_edges_.insert_row();
  node_faces_.insert_row();
  node_cells_.insert_row();

  return node_index;

} // Mesh::insert_node

EdgeIndex Mesh::insert_edge(std::unique_ptr<Shape>&& edge_shape,
                            std::optional<EdgeMark> edge_mark) {
  // Try to find an existing edge.
  const auto [edge_iter, edge_inserted] = //
      edges_lookup_table_.try_emplace(
          ranges::to<std::set>(edge_shape->node_indices()), num_edges_);
  const EdgeIndex edge_index = edge_iter->second;
  if (!edge_inserted) {
    if (edge_mark.has_value()) { edge_marks_[edge_index] = *edge_mark; }
    return edge_index;
  }

  // Emplace the edge properties.
  num_edges_ += 1;
  edge_marks_.emplace_back(edge_mark.value_or(EdgeMark{0}));
  edge_shapes_.emplace_back(edge_shape->shape_type());
  edge_lens_.emplace_back(edge_shape->volume(node_coords_));
  edge_dirs_.emplace_back(edge_shape->dir(node_coords_));

  // Fill the edge-nodes and node-edges, edge-edges and node-nodes.
  edge_nodes_.insert_row(edge_shape->node_indices());
  for (NodeIndex node_index : edge_nodes_[edge_index]) {
    node_edges_.insert(node_index, edge_index);
  }
  edge_edges_.insert_row([&]() {
    std::vector<EdgeIndex> this_edge_edges{};
    for (NodeIndex node_index :
         {edge_nodes_[edge_index].front(), edge_nodes_[edge_index].back()}) {
      for (EdgeIndex other_edge_index : node_edges_[node_index]) {
        this_edge_edges.emplace_back(other_edge_index);
        if (edge_index != other_edge_index) {
          edge_edges_.insert(other_edge_index, edge_index);
        }
      }
    }
    return this_edge_edges;
  }());
  node_nodes_.insert(edge_nodes_[edge_index].front(),
                     edge_nodes_[edge_index].back());
  node_nodes_.insert(edge_nodes_[edge_index].back(),
                     edge_nodes_[edge_index].front());

  // Emplace empty rows that would be filled later on.
  edge_edges_.insert_row();
  edge_faces_.insert_row();
  edge_cells_.insert_row();

  return edge_index;

} // Mesh::insert_edge

FaceIndex Mesh::insert_face(std::unique_ptr<Shape>&& face_shape,
                            std::optional<FaceMark> face_mark) {
  // Try to find an existing face.
  const auto [face_iter, face_inserted] = //
      faces_lookup_table_.try_emplace(
          ranges::to<std::set>(face_shape->node_indices()), num_faces_);
  const FaceIndex face_index = face_iter->second;
  if (!face_inserted) {
    if (face_mark.has_value()) { face_marks_[face_index] = *face_mark; }
    return face_index;
  }

  // Emplace the face properties.
  num_faces_ += 1;
  face_marks_.emplace_back(face_mark.value_or(FaceMark{0}));
  face_shapes_.emplace_back(face_shape->shape_type());
  face_areas_.emplace_back(face_shape->volume(node_coords_));
  face_normals_.emplace_back(face_shape->normal(node_coords_));
  face_centers_.emplace_back(face_shape->center_pos(node_coords_));

  // Fill the face-nodes and node-faces.
  face_nodes_.insert_row(face_shape->node_indices());
  for (NodeIndex node_index : face_nodes_[face_index]) {
    node_faces_.insert(node_index, face_index);
  }

  // Fill the face-edges, edge-faces and face-faces.
  face_edges_.insert_row([&]() {
    small_vector<EdgeIndex, 4> this_face_edges{};
    for (const ShapeDesc& edges_desc : face_shape->make_edges_desc()) {
      const EdgeIndex edge_index{insert_edge(edges_desc)};
      this_face_edges.emplace_back(edge_index);
      edge_faces_.insert(edge_index, face_index);
    }
    return this_face_edges;
  }());
  face_faces_.insert_row([&]() {
    std::vector<FaceIndex> this_face_faces{};
    for (EdgeIndex edge_index : face_edges_[face_index]) {
      for (FaceIndex other_face_index : edge_faces_[edge_index]) {
        this_face_faces.emplace_back(other_face_index);
        if (face_index != other_face_index) {
          face_faces_.insert(other_face_index, face_index);
        }
      }
    }
    return this_face_faces;
  }());

  // Emplace empty rows that would be filled later on.
  face_cells_.insert_row();

  return face_index;

} // Mesh::insert_face

CellIndex Mesh::insert_cell(std::unique_ptr<Shape>&& cell_shape,
                            CellMark cell_mark, bool ghost) {
  // Allocate the cell index.
  const CellIndex cell_index{num_cells_};

  // Emplace cell properties.
  num_cells_ += 1;
  cell_marks_.emplace_back(cell_mark);
  cell_shapes_.emplace_back(cell_shape->shape_type());
  cell_volumes_.emplace_back(cell_shape->volume(node_coords_));
  cell_centers_.emplace_back(cell_shape->center_pos(node_coords_));

  // Fill the cell nodes and node cells.
  cell_nodes_.insert_row(cell_shape->node_indices());
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
    small_vector<EdgeIndex, 12> this_cell_edges{};
    for (const ShapeDesc& edges_desc : cell_shape->make_edges_desc()) {
      const EdgeIndex edge_index{insert_edge(edges_desc)};
      this_cell_edges.emplace_back(edge_index);
      edge_cells_.insert(edge_index, cell_index);
    }
    return this_cell_edges;
  }());

  // Fill the cell faces and face cells.
  cell_faces_.insert_row([&]() {
    small_vector<FaceIndex, 6> this_cell_faces{};
    for (const ShapeDesc& face_desc : cell_shape->make_faces_desc()) {
      const FaceIndex face_index{insert_face(face_desc)};
      this_cell_faces.emplace_back(face_index);
      if (bool face_oriented_as_inner =
              ranges::equal(adjacent_nodes(face_index), face_desc.node_indices);
          face_cells_[face_index].empty() && !face_oriented_as_inner) {
        flip_face(face_index);
      } else {
        STORM_ASSERT_(face_oriented_as_inner);
      }
      face_cells_.insert(face_index, cell_index);
    }
    return this_cell_faces;
  }());

  /// @todo Fill me!
  cell_cells_.insert_row();

  return cell_index;

} // Mesh::insert_cell

// ------------------------------------------------------------------------------------
// //
// ------------------------------------------------------------------------------------
// //

void Mesh::flip_face(FaceIndex face_index) noexcept {
  STORM_ASSERT_(face_index < num_faces_ && "face_index is out of range");
  STORM_ASSERT_([&]() {
    auto const this_face_cells{face_cells_[face_index]};
    return this_face_cells.size() != 1 &&
           ranges::all_of(this_face_cells, [this](CellIndex cell_index) {
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
  InversePermutation(permutation.begin(), permutation.end(),
                     inversePermutation.begin());

  auto const reorderFunc = [&inversePermutation](Index<Tag> index) {
    return Index<Tag>(index != npos ? inversePermutation[size_t(index)] : npos);
  };

  ranges::for_each(nodes(), [&](NodeIndex node_index) {
    ranges::transform(adjacent_elements_<Tag>(node_index),
                      std::begin(adjacent_elements_<Tag>(node_index)),
                      reorderFunc);
  });

  ranges::for_each(edges(), [&](EdgeIndex edge_index) {
    ranges::transform(adjacent_elements_<Tag>(edge_index),
                      std::begin(adjacent_elements_<Tag>(edge_index)),
                      reorderFunc);
  });

  ranges::for_each(faces(), [&](FaceIndex face_index) {
    ranges::transform(adjacent_elements_<Tag>(face_index),
                      std::begin(adjacent_elements_<Tag>(face_index)),
                      reorderFunc);
  });

  ranges::for_each(cells(), [&](CellIndex cell_index) {
    ranges::transform(adjacent_elements_<Tag>(cell_index),
                      std::begin(adjacent_elements_<Tag>(cell_index)),
                      reorderFunc);
  });

} // Mesh::FixPermutationAndAdjacency_

void Mesh::PermuteNodes(std::vector<size_t>&& nodePermutation) {
  /* Permute Node properties and fix the adjacency tables. */
  FixPermutationAndAdjacency_<NodeTag>(nodePermutation);
  permute_rows(((std::vector<NodeIndex>&) nodePermutation).begin(),
               ((std::vector<NodeIndex>&) nodePermutation).end(), node_nodes_,
               node_edges_, node_faces_, node_cells_);
  permute_inplace(nodePermutation.begin(), nodePermutation.end(),
                  node_marks_.begin(), node_coords_.begin());

  /* Generate the Node ranges. */
  node_ranges_.clear();
  ranges::for_each(node_views(*this), [&](NodeView node) {
    node_ranges_.resize((size_t) node.mark() + 2);
    node_ranges_[node.mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&) node_ranges_;
  std::partial_sum(r.begin(), r.end(), r.begin());

} // Mesh::PermuteNodes

void Mesh::PermuteEdges(std::vector<size_t>&& edgePermutation) {
  /* Permute edge properties and fix the adjacency tables. */
  FixPermutationAndAdjacency_<EdgeTag>(edgePermutation);
  permute_rows(((std::vector<EdgeIndex>&) edgePermutation).begin(),
               ((std::vector<EdgeIndex>&) edgePermutation).end(), edge_nodes_,
               edge_edges_, edge_faces_, edge_cells_);
  permute_inplace(edgePermutation.begin(), edgePermutation.end(),
                  edge_marks_.begin(), edge_shapes_.begin(), edge_lens_.begin(),
                  edge_dirs_.begin());

  /* Generate the edge ranges. */
  edge_ranges_.clear();
  ranges::for_each(edge_views(*this), [&](EdgeView edge) {
    edge_ranges_.resize((size_t) edge.mark() + 2);
    edge_ranges_[edge.mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&) edge_ranges_;
  std::partial_sum(r.begin(), r.end(), r.begin());

} // Mesh::PermuteEdges

void Mesh::PermuteFaces(std::vector<size_t>&& facePermutation) {
  /* Permute data. */
  FixPermutationAndAdjacency_<FaceTag>(facePermutation);
  permute_rows(((std::vector<FaceIndex>&) facePermutation).begin(),
               ((std::vector<FaceIndex>&) facePermutation).end(), face_nodes_,
               face_edges_, face_faces_, face_cells_);
  permute_inplace(facePermutation.begin(), facePermutation.end(),
                  face_marks_.begin(), face_shapes_.begin(),
                  face_areas_.begin(), face_normals_.begin(),
                  face_centers_.begin());

  /* Generate Mark ranges. */
  face_ranges_.clear();
  ranges::for_each(face_views(*this), [&](FaceView face) {
    face_ranges_.resize((size_t) face.mark() + 2);
    face_ranges_[face.mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&) face_ranges_;
  std::partial_sum(r.begin(), r.end(), r.begin());

} // Mesh::PermuteFaces

void Mesh::PermuteCells(std::vector<size_t>&& cellPermutation) {
  /* Permute data. */
  FixPermutationAndAdjacency_<CellTag>(cellPermutation);
  permute_rows(((std::vector<CellIndex>&) cellPermutation).begin(),
               ((std::vector<CellIndex>&) cellPermutation).end(), cell_nodes_,
               cell_edges_, cell_faces_, cell_cells_);
  permute_inplace(cellPermutation.begin(), cellPermutation.end(),
                  cell_marks_.begin(), cell_shapes_.begin(),
                  cell_volumes_.begin(), cell_centers_.begin());

  /* Generate Mark ranges. */
  cell_ranges_.clear();
  ranges::for_each(cell_views(*this), [&](CellView cell) {
    cell_ranges_.resize((size_t) cell.mark() + 2);
    cell_ranges_[cell.mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&) cell_ranges_;
  std::partial_sum(r.begin(), r.end(), r.begin());

} // Mesh::PermuteCells

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

/// ---------------------------------------------------------------- ///
/// ---------------------------------------------------------------- ///

/**
 * Generate boundary cells to complete face connectivity.
 */
void Mesh::generate_boundary_cells(FaceIndex ff) {
  auto const f = [&](MutableFaceView face) {
    if (face.mark() == 0) { return; }

    CellView cell = face.inner_cell();

    /* Generate the boundary cell: reflect a connected interior cell. */
    std::vector<NodeIndex> ghost_cell_nodes;
    cell.for_each_node([&](NodeView node) {
      if (std::find(face.adjacent_nodes().begin(), face.adjacent_nodes().end(),
                    node) == face.adjacent_nodes().end()) {
        /* Reflect an interior cell Node. */
        // TODO: face normals are not computed here!
        // TODO:
        // https://glm.g-truc.net/0.9.5/api/a00157.html#gab63646fc36b81cf69d3ce123a72f76f2
        vec3_t node_coords = node.coords();
        const vec3_t delta = node_coords - face.center();
        node_coords -= 2.0 * glm::dot(delta, face.normal()) * face.normal();
        ghost_cell_nodes.push_back(
            insert_node(node_coords, (NodeMark) face.mark()));
      } else {
        /* Insert a boundary face Node. */
        ghost_cell_nodes.push_back(node);
      }
    });

    /* Insert the boundary cell3. */
    const CellIndex bndr_cell_index = insert_cell(
        {cell.shape_type(), ghost_cell_nodes}, (CellMark) face.mark(), true);
    face_cells_.insert(face, bndr_cell_index);
  };

  if (ff != npos) {
    f(MutableFaceView(*this, ff));
    return;
  }
  ranges::for_each(face_views(*this), f);

} // Mesh::generate_boundary_cells

} // namespace Storm
