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

#include <set>
#include <map>

#include <stormMesh/Base.hxx>
#include <stormUtils/Table.hh>
#include <stormUtils/Image.hh>
#include "Field.hh"
#include "Element.hh"

namespace Storm {

inline constexpr size_t face_inner_cell_{0};
inline constexpr size_t face_outer_cell_{1};

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
//template<size_t SDim, size_t TDim, template<class, class, class> class Table>
class Mesh : public tObject<Mesh> {
private:
  size_t num_nodes_{0}, num_edges_{0};
  size_t num_faces_{0}, num_cells_{0};

  Vector<NodeMark, NodeIndex> node_marks_;
  Vector<EdgeMark, EdgeIndex> edge_marks_;
  Vector<FaceMark, FaceIndex> face_marks_;
  Vector<CellMark, CellIndex> cell_marks_;
  Vector<NodeIndex, NodeMark> node_ranges_;
  Vector<EdgeIndex, EdgeMark> edge_ranges_;
  Vector<FaceIndex, FaceMark> face_ranges_;
  Vector<CellIndex, CellMark> cell_ranges_;

  Vector<ShapeType, EdgeIndex> edge_shapes_;
  Vector<ShapeType, FaceIndex> face_shapes_;
  Vector<ShapeType, CellIndex> cell_shapes_;
  // TODO: edge length + direction -> 4D oriented direction.
  // TODO: face area + normal -> 4D oriented area.
  Vector<vec3_t, NodeIndex> node_coords_;
  Vector<real_t, EdgeIndex> edge_lens_;
  Vector<vec3_t, EdgeIndex> edge_dirs_;
  Vector<real_t, FaceIndex> face_areas_;
  Vector<vec3_t, FaceIndex> face_normals_;
  Vector<vec3_t, FaceIndex> face_barycenters_;
  Vector<real_t, CellIndex> cell_volumes_;
  Vector<vec3_t, CellIndex> cell_barycenters_;
  real_t min_edge_len_{qnan}, max_edge_len_{qnan};
  real_t min_face_area_{qnan}, max_face_area_{qnan};
  real_t min_cell_volume_{qnan}, max_cell_volume_{qnan};

  CsrTable<EdgeIndex, NodeIndex> edge_nodes_;
  CsrTable<FaceIndex, NodeIndex> face_nodes_;
  CsrTable<CellIndex, NodeIndex> cell_nodes_;
  CsrTable<FaceIndex, EdgeIndex> face_edges_;
  CsrTable<CellIndex, EdgeIndex> cell_edges_;
  CsrTable<CellIndex, FaceIndex> cell_faces_;

  CsrTable<NodeIndex, EdgeIndex> node_edges_;
  CsrTable<NodeIndex, FaceIndex> node_faces_;
  CsrTable<NodeIndex, CellIndex> node_cells_;
  CsrTable<EdgeIndex, FaceIndex> edge_faces_;
  CsrTable<EdgeIndex, CellIndex> edge_cells_;
  CsrTable<FaceIndex, CellIndex> face_cells_;

  CsrTable<NodeIndex, NodeIndex> node_nodes_;
  CsrTable<EdgeIndex, EdgeIndex> edge_edges_;
  CsrTable<FaceIndex, FaceIndex> face_faces_;
  CsrTable<CellIndex, CellIndex> cell_cells_;

  std::map<std::set<NodeIndex>, EdgeIndex> edge_lookup_;
  std::map<std::set<NodeIndex>, FaceIndex> face_lookup_;
  std::map<std::set<NodeIndex>, CellIndex> cell_lookup_;

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

  void save_vtk(const char* path,
                const std::vector<sFieldDesc>& fields) const;

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
  auto nodes() const noexcept {
    return views::iota(NodeIndex{0}, NodeIndex{num_nodes_});
  }

  /// @brief Range of edge indices.
  auto edges() const noexcept {
    return views::iota(EdgeIndex{0}, EdgeIndex{num_edges_});
  }

  /// @brief Range of face indices.
  auto faces() const noexcept {
    return views::iota(FaceIndex{0}, FaceIndex{num_faces_});
  }

  /// @brief Range of cell indices.
  auto cells() const noexcept {
    return views::iota(CellIndex{0}, CellIndex{num_cells_});
  }

  /// @}

  /// ---------------------------------------------------------------- ///
  /// @name Marks.
  /// ---------------------------------------------------------------- ///
  /// @{

  /// @brief Number of node marks.
  size_t num_node_marks() const noexcept {
    storm_assert(!node_ranges_.empty());
    return node_ranges_.size() - 1;
  }

  /// @brief Number of edge marks.
  size_t num_edge_marks() const noexcept {
    storm_assert(!edge_ranges_.empty());
    return edge_ranges_.size() - 1;
  }

  /// @brief Number of face marks.
  size_t num_face_marks() const noexcept {
    storm_assert(!face_ranges_.empty());
    return face_ranges_.size() - 1;
  }

  /// @brief Number of cell marks.
  size_t num_cell_marks() const noexcept {
    storm_assert(!cell_ranges_.empty());
    return cell_ranges_.size() - 1;
  }

  /// @brief Range of node indices with a @p node_mark.
  auto nodes(NodeMark node_mark) const noexcept {
    storm_assert(node_mark < num_node_marks());
    return views::iota(node_ranges_[node_mark], node_ranges_[node_mark + 1]);
  }

  /// @brief Range of edge indices with a @p edge_mark.
  auto edges(EdgeMark edge_mark) const noexcept {
    storm_assert(edge_mark < num_edge_marks());
    return views::iota(edge_ranges_[edge_mark], edge_ranges_[edge_mark + 1]);
  }

  /// @brief Range of face indices with a @p face_mark.
  auto faces(FaceMark face_mark) const noexcept {
    storm_assert(face_mark < num_face_marks());
    return views::iota(face_ranges_[face_mark], face_ranges_[face_mark + 1]);
  }

  /// @brief Range of cell indices with a @p cell_mark.
  auto cells(CellMark cell_mark) const noexcept {
    storm_assert(cell_mark < num_cell_marks());
    return views::iota(cell_ranges_[cell_mark], cell_ranges_[cell_mark + 1]);
  }

  /// @brief Get node @p node_index mark.
  NodeMark mark(NodeIndex node_index) const noexcept {
    storm_assert(node_index < num_nodes_);
    return node_marks_[node_index];
  }

  /// @brief Get edge @p edge_index mark.
  EdgeMark mark(EdgeIndex edge_index) const noexcept {
    storm_assert(edge_index < num_edges_);
    return edge_marks_[edge_index];
  }

  /// @brief Get face @p face_index mark.
  FaceMark mark(FaceIndex face_index) const noexcept {
    storm_assert(face_index < num_faces_);
    return face_marks_[face_index];
  }

  /// @brief Get cell @p cell_index mark.
  CellMark mark(CellIndex cell_index) const noexcept {
    storm_assert(cell_index < num_cells_);
    return cell_marks_[cell_index];
  }

  /// @}

  /// ---------------------------------------------------------------- ///
  /// @name Shapes.
  /// ---------------------------------------------------------------- ///
  /// @{

  /// @brief Get edge @p edge_index shape type.
  ShapeType shapeType(EdgeIndex edge_index) const noexcept {
    storm_assert(edge_index < num_edges_);
    return edge_shapes_[edge_index];
  }

  /// @brief Get face @p face_index shape type.
  ShapeType shapeType(FaceIndex face_index) const noexcept {
    storm_assert(face_index < num_faces_);
    return face_shapes_[face_index];
  }

  /// @brief Get cell @p cell_index shape type.
  ShapeType shapeType(CellIndex cell_index) const noexcept {
    storm_assert(cell_index < num_cells_);
    return cell_shapes_[cell_index];
  }

private:

  auto makeShape_(ShapeDesc&& desc) const {
    return Element::Make(std::forward<ShapeDesc>(desc), node_coords_);
  }

public:

  /// @brief Get element object.
  template<class Tag>
  auto shape(Index<Tag> index) const {
    auto const node_indices = adjacent_nodes(index);
    return makeShape_({shapeType(index), std::vector(node_indices.begin(), node_indices.end())});
  }

  /// @brief Get node @p node_index coordinates.
  vec3_t node_coords(NodeIndex node_index) const noexcept {
    storm_assert(node_index < num_nodes_ && "node_index is out of range");
    return node_coords_[node_index];
  }

  /// @brief Set node @p node_index position @p coords.
  void set_node_coords(NodeIndex node_index, vec3_t const& coords) noexcept {
    storm_assert(node_index < num_nodes_ && "node_index is out of range");
    node_coords_[node_index] = coords;
  }

  /// @brief Get edge @p edge_index length.
  real_t edge_len(EdgeIndex edge_index) const noexcept {
    storm_assert(edge_index < num_edges_ && "edge_index is out of range");
    return edge_lens_[edge_index];
  }

  /// @brief Get edge @p edge_index direction.
  vec3_t edge_dir(EdgeIndex edge_index) const noexcept {
    storm_assert(edge_index < num_edges_ && "edge_index is out of range");
    return edge_dirs_[edge_index];
  }

  /// @brief Get face @p face_index area/length.
  real_t face_area(FaceIndex face_index) const noexcept {
    storm_assert(face_index < num_faces_ && "face_index is out of range");
    return face_areas_[face_index];
  }

  /// @brief Get face @p face_index normal.
  vec3_t face_normal(FaceIndex face_index) const noexcept {
    storm_assert(face_index < num_faces_ && "face_index is out of range");
    return face_normals_[face_index];
  }

  /// @brief Get face @p face_index barycenter.
  vec3_t face_barycenter(FaceIndex face_index) const noexcept {
    storm_assert(face_index < num_faces_ && "face_index is out of range");
    return face_barycenters_[face_index];
  }

  /// @brief Get cell @p cell_index volume/area/length.
  real_t cell_volume(CellIndex cell_index) const noexcept {
    storm_assert(cell_index < num_cells_ && "cell_index is out of range");
    return cell_volumes_[cell_index];
  }

  /// @brief Get cell @p cell_index barycenter.
  vec3_t cell_barycenter(CellIndex cell_index) const noexcept {
    storm_assert(cell_index < num_cells_ && "cell_index is out of range");
    return cell_barycenters_[cell_index];
  }

  /// @brief Get the minimal edge length.
  real_t min_edge_len() const noexcept {
    storm_assert(std::isfinite(min_cell_volume_));
    return min_edge_len_;
  }

  /// @brief Get the maximal edge length.
  real_t max_edge_len() const noexcept {
    storm_assert(std::isfinite(min_cell_volume_));
    return max_edge_len_;
  }

  /// @brief Get the minimal face area.
  real_t min_face_area() const noexcept {
    storm_assert(std::isfinite(min_cell_volume_));
    return min_face_area_;
  }

  /// @brief Get the maximal face area.
  real_t max_face_area() const noexcept {
    storm_assert(std::isfinite(min_cell_volume_));
    return max_face_area_;
  }

  /// @brief Get the minimal cell volume.
  real_t min_cell_volume() const noexcept {
    storm_assert(std::isfinite(min_cell_volume_));
    return min_cell_volume_;
  }

  /// @brief Get the maximal cell volume.
  real_t max_cell_volume() const noexcept {
    storm_assert(std::isfinite(max_cell_volume_));
    return max_cell_volume_;
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

  /// @brief Range of the edge @p edge_index adjacent node indices.
  /// Denote a node to be adjacent to an edge if it one its nodes.
  auto adjacent_nodes(EdgeIndex edge_index) const noexcept {
    storm_assert(edge_index < num_edges_ && "edge_index is out of range");
    return const_cast<Mesh*>(this)->edge_nodes_[edge_index];
  }

  /// @brief Range of the face @p face_index adjacent node indices.
  /// Denote a node to be adjacent to a face if it one its nodes.
  auto adjacent_nodes(FaceIndex face_index) const noexcept {
    storm_assert(face_index < num_faces_ && "face_index is out of range");
    return const_cast<Mesh*>(this)->face_nodes_[face_index];
  }

  /// @brief Range of the cell @p cell_index adjacent node indices.
  /// Denote a node to be adjacent to a cell if it one its nodes.
  auto adjacent_nodes(CellIndex cell_index) const noexcept {
    storm_assert(cell_index < num_cells_ && "cell_index is out of range");
    return const_cast<Mesh*>(this)->cell_nodes_[cell_index];
  }

  /// @brief Range of the face @p face_index adjacent edge indices.
  /// Denote an edge to be adjacent to a face if it one its edges.
  auto adjacent_edges(FaceIndex face_index) const noexcept {
    storm_assert(face_index < num_faces_ && "face_index is out of range");
    return const_cast<Mesh*>(this)->face_edges_[face_index];
  }

  /// @brief Range of the cell @p cell_index adjacent edge indices.
  /// Denote an edge to be adjacent to a cell if it one its edges.
  auto adjacent_edges(CellIndex cell_index) const noexcept {
    storm_assert(cell_index < num_cells_ && "cell_index is out of range");
    return const_cast<Mesh*>(this)->cell_edges_[cell_index];
  }

  /// @brief Range of the cell @p cell_index adjacent face indices.
  /// Denote a face to be adjacent to a cell if it one its faces.
  auto adjacent_faces(CellIndex cell_index) const noexcept {
    storm_assert(cell_index < num_cells_ && "cell_index is out of range");
    return const_cast<Mesh*>(this)->cell_faces_[cell_index];
  }

  /// @}

  /// @name Secondary adjacency: adjacency,
  ///   that is a transpose of the primary adjacency.
  /// @{

  /// @brief Range of the node @p node_index adjacent edge indices.
  auto adjacent_edges(NodeIndex node_index) const noexcept {
    storm_assert(node_index < num_nodes_ && "node_index is out of range");
    return const_cast<Mesh*>(this)->node_edges_[node_index];
  }

  /// @brief Range of the node @p node_index adjacent faces indices.
  auto adjacent_faces(NodeIndex node_index) const noexcept {
    storm_assert(node_index < num_nodes_ && "node_index is out of range");
    return const_cast<Mesh*>(this)->node_faces_[node_index];
  }

  /// @brief Range of the node @p node_index adjacent cell indices.
  auto adjacent_cells(NodeIndex node_index) const noexcept {
    storm_assert(node_index < num_nodes_ && "node_index is out of range");
    return const_cast<Mesh*>(this)->node_cells_[node_index];
  }

  /// @brief Range of the edge @p edge_index adjacent face indices.
  auto adjacent_faces(EdgeIndex edge_index) const noexcept {
    storm_assert(edge_index < num_edges_ && "edge_index is out of range");
    return const_cast<Mesh*>(this)->edge_faces_[edge_index];
  }

  /// @brief Range of the edge @p edge_index adjacent cell indices.
  auto adjacent_cells(EdgeIndex edge_index) const noexcept {
    storm_assert(edge_index < num_edges_ && "edge_index is out of range");
    return const_cast<Mesh*>(this)->edge_cells_[edge_index];
  }

  /// @brief Range of the face @p face_index adjacent cell indices.
  auto adjacent_cells(FaceIndex face_index) const noexcept {
    storm_assert(face_index < num_faces_ && "face_index is out of range");
    return const_cast<Mesh*>(this)->face_cells_[face_index];
  }

  /// @}

  /// @name Symmetric adjacency: adjacency, that is a product
  ///   of the corresponding primary and secondary adjacencies.
  /// @{

  /// @brief Range of the node @p node_index adjacent node indices.
  /// Denote two nodes as adjacent if there is an edge connecting them.
  auto adjacent_nodes(NodeIndex node_index) const noexcept {
    storm_assert(node_index < num_nodes_ && "node_index is out of range");
    return const_cast<Mesh*>(this)->node_nodes_[node_index];
  }

  /// @brief Range of the edge @p edge_index adjacent edge indices.
  /// Denote two edges as adjacent if they share a common node.
  auto adjacent_edges(EdgeIndex edge_index) const noexcept {
    storm_assert(edge_index < num_edges_ && "edge_index is out of range");
    return const_cast<Mesh*>(this)->edge_edges_[edge_index];
  }

  /// @brief Range of the face @p face_index adjacent face indices.
  /// Denote two faces as adjacent if they share a common edge.
  auto adjacent_faces(FaceIndex face_index) const noexcept {
    storm_assert(face_index < num_faces_ && "face_index is out of range");
    return const_cast<Mesh*>(this)->face_faces_[face_index];
  }

  /// @brief Range of the cell @p cell_index adjacent cell indices.
  /// Denote two cells as adjacent if they share a common face.
  auto adjacent_cells(CellIndex cell_index) const noexcept {
    storm_assert(cell_index < num_cells_ && "cell_index is out of range");
    return const_cast<Mesh*>(this)->cell_cells_[cell_index];
  }

  /// @}

private:

  template<class OutTag, class Tag>
  auto adjacent_elements(Index<Tag> element_index) const noexcept {
    if constexpr (std::same_as<OutTag, NodeTag>) {
      return adjacent_nodes(element_index);
    } else if constexpr (std::same_as<OutTag, EdgeTag>) {
      return adjacent_edges(element_index);
    } else if constexpr (std::same_as<OutTag, FaceTag>) {
      return adjacent_faces(element_index);
    } else if constexpr (std::same_as<OutTag, CellTag>) {
      return adjacent_cells(element_index);
    } else {
      static_assert(always_false<OutTag>, "Invalid tag.");
    }
  }

public:

  /// @}

  /// ---------------------------------------------------------------- ///
  /// @name Insertions.
  /// ---------------------------------------------------------------- ///
  /// @{

  /// @brief Insert a new node with a position @p coords
  ///   and a mark @p node_mark into the mesh.
  /// @returns Index of the inserted node.
  NodeIndex insert_node(vec3_t const& coords,
                        NodeMark node_mark = {});

  /// @brief Find or emplace a new edge_shape with a shape @p edgeShape
  ///   and node indices @p nodes.
  ///
  /// Insertion would update the edge_shape-node, edge_shape-edge_shape and
  ///   node-node topologies.
  ///
  /// @param edge_mark Mark that would be assigned to an edge_shape if
  ///   the insertion took place, otherwise ignored.
  /// @returns A pair of an index of the found or inserted edge_shape
  ///   and a boolean value denoting whether the insertion took place.
  //std::pair<EdgeIndex, bool>
  //  findOrInsertEdge(std::shared_ptr<Shape> edgeShape,
  //                   EdgeMark edge_mark = {});
  EdgeIndex insert_edge(std::unique_ptr<Element>&& edge_shape, EdgeMark edge_mark = {});
  EdgeIndex insert_edge(ShapeDesc&& edgeDesc, EdgeMark edge_mark = {}) {
    return insert_edge(makeShape_(std::forward<ShapeDesc>(edgeDesc)), edge_mark);
  }

  /// @brief Find or emplace a new face_shape with a shape @p faceShape
  ///   and node indices @p nodes.
  ///
  /// Insertion would implicitly insert the missing edges and
  ///   update the face_shape-edge and face_shape-face_shape topologies. The
  ///   implicitly inserted edges would inherit the mark from
  ///   explicitly inserted face_shape.
  ///
  /// @param face_mark Mark that would be assigned to a face_shape if
  ///   the insertion took place, otherwise ignored.
  /// @returns A pair of an index of the found or inserted face_shape
  ///   and a boolean value denoting whether the insertion took place.
  //std::pair<FaceIndex, bool>
  //  findOrInsertFace(std::shared_ptr<Shape> faceShape,
  //                   FaceMark face_mark = {});
  FaceIndex insert_face(std::unique_ptr<Element>&& face_shape, FaceMark face_mark = {});
  FaceIndex insert_face(ShapeDesc&& faceDesc, FaceMark faceMark = {}) {
    return insert_face(makeShape_(std::forward<ShapeDesc>(faceDesc)), faceMark);
  }

  /// @brief Find or emplace a new cell with a shape @p cellShape
  ///   and node indices @p nodes.
  ///
  /// Insertion would implicitly insert the missing faces and
  ///   update the cell-face and cell-cell topologies. The
  ///   implicitly inserted face would inherit the mark from
  ///   explicitly inserted cell.
  ///
  /// @param cell_mark Mark that would be assigned to a cell if
  ///   the insertion took place, otherwise ignored.
  /// @returns A pair of an index of the found or inserted cell
  ///   and a boolean value denoting whether the insertion took place.
  //std::pair<CellIndex, bool>
  //  findOrInsertCell(std::shared_ptr<Shape> cellShape,
  //                   CellMark cell_mark = {});
  CellIndex insert_cell(std::unique_ptr<Element>&& cell, CellMark cell_mark = {}, bool ghost = false);
  CellIndex insert_cell(ShapeDesc&& cellDesc, CellMark cellMark = {}, bool ghost = false) {
    return insert_cell(makeShape_(std::forward<ShapeDesc>(cellDesc)), cellMark, ghost);
  }

  /// @brief Generate boundary cells to complete face connectivity.
  void generate_boundary_cells(FaceIndex ff = FaceIndex{npos});

  /// @}

  /// ---------------------------------------------------------------- ///
  /// @name Permutations.
  /// ---------------------------------------------------------------- ///
  /// @{

  /// @brief Flip the face @p face_index.
  /// A face can be flipped only if it is not adjacent to any cells
  /// or it is adjacent to the exactly two interior cells.
  void flip_face(FaceIndex face_index) noexcept;

  /// @brief Change order of all nodes.
  void PermuteNodes(std::vector<size_t>&& nodePermutation = {});

  /// @brief Change order of all edges.
  void PermuteEdges(std::vector<size_t>&& edgePermutation = {});

  /// @brief Change order of all faces.
  void PermuteFaces(std::vector<size_t>&& facePermutation = {});

  /// @brief Change order of all cells.
  void PermuteCells(std::vector<size_t>&& cellPermutation = {});

private:

  template<class Tag>
  void FixPermutationAndAdjacency_(std::vector<size_t>& cell_index);

protected:

  void reorder_faces();

  /// @}

}; // class Mesh

} // namespace feathers

// Include iterators.
#include "View.hxx"

/// @brief @todo Remove me.
using cMesh = Storm::Mesh;
