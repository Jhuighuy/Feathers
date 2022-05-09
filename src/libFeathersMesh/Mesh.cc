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

#include <libFeathersMesh/Mesh.hh>
#include <libFeathersMesh/Element.hh>
#include <libFeathersUtils/Parallel.hh>
#include <libFeathersUtils/Permute.hh>

#include <set>
#include <map>

namespace feathers {

void Mesh::ComputeEdgeShapeProperties() {

  std::tie(MinEdgeLen_, MaxEdgeLen_) =
    for_range_minmax(feathers::BeginEdge(*this), feathers::EndEdge(*this),
                     +huge, -huge, [&](EdgeMutableRef edge) {

        std::unique_ptr<const iElement> edgeObj = edge.get_element_object();

        edge.SetLen(edgeObj->get_length_or_area_or_volume());
        edge.SetDir(edgeObj->get_direction());

        return edge.Len();

      });

} // Mesh<...>::ComputeEdgeShapeProperties

void Mesh::ComputeFaceShapeProperties() {

  std::tie(MinFaceArea_, MaxFaceArea_) =
    for_range_minmax(feathers::BeginFace(*this), feathers::EndFace(*this),
                     +huge, -huge, [&](FaceMutableRef face) {

        std::unique_ptr<const iElement> face_element = face.get_element_object();

        face.SetArea(face_element->get_length_or_area_or_volume());
        face.SetNormal(face_element->get_normal());
        face.SetCenterPos(face_element->get_center_coords());

        return face.Area();

      });

} // Mesh<...>::ComputeFaceShapeProperties

void Mesh::ComputeCellShapeProperties() {

  std::tie(MinFaceArea_, MaxFaceArea_) =
    for_range_minmax(feathers::BeginCell(*this), feathers::EndCell(*this),
                     +huge, -huge, [&](CellMutableRef cell) {

        std::unique_ptr<const iElement> cell_element = cell.get_element_object();

        cell.SetVolume(cell_element->get_length_or_area_or_volume());
        cell.SetCenterPos(cell_element->get_center_coords());
        return cell.Volume();

      });

} // Mesh<...>::ComputeCellShapeProperties

NodeIndex Mesh::EmplaceNode(vec3_t const& nodePos, NodeMark nodeMark) {

  NodeIndex const nodeIndex(NumNodes_++);

  // Emplace the node properties.
  NodeMarks_.emplace_back(nodeMark);

  NodePos_.emplace_back(nodePos);

  // Emplace empty rows into the all node
  // adjacency tables to keep the mesh consistent.
  NodeNodes_.emplace_back_row(/* dynamic value */);
  NodeEdges_.emplace_back_row(/* dynamic value */);
  NodeFaces_.emplace_back_row(/* dynamic value */);
  NodeCells_.emplace_back_row(/* dynamic value */);

  return nodeIndex;

} // Mesh<...>::EmplaceNode

EdgeIndex Mesh::EmplaceEdge(std::unique_ptr<iElement>&& edge, EdgeMark edgeMark) {

  EdgeIndex const edgeIndex(NumEdges_++);

  // Emplace the edge properties.
  EdgeMarks_.emplace_back(edgeMark);

  EdgeShapes_.emplace_back(edge->get_shape());
  EdgeLens_.emplace_back(edge->get_length_or_area_or_volume());
  EdgeDirs_.emplace_back(edge->get_direction());

  // Fill the edge nodes.
  ((CsrTable<size_t, size_t>&) EdgeNodes_).emplace_back_row(edge->get_nodes().begin(), edge->get_nodes().end());

  // Emplace empty rows into the remaining edge
  // adjacency tables to keep the mesh consistent. */
  EdgeEdges_.emplace_back_row(/* dynamic value */);
  EdgeFaces_.emplace_back_row(/* dynamic value */);
  EdgeCells_.emplace_back_row(/* dynamic value */);

  return edgeIndex;

} // Mesh<...>::EmplaceEdge

FaceIndex Mesh::EmplaceFace(std::unique_ptr<iElement>&& face, FaceMark faceMark) {

  FaceIndex const faceIndex(NumFaces_++);

  // Emplace the face properties.
  FaceMarks_.emplace_back(faceMark);

  FaceShapes_.emplace_back(face->get_shape());
  FaceAreas_.emplace_back(face->get_length_or_area_or_volume());
  FaceNormals_.emplace_back(face->get_normal());
  FaceCenterPos_.emplace_back(face->get_center_coords());

  // Fill the face nodes.
  ((CsrTable<size_t, size_t>&) FaceNodes_).emplace_back_row(face->get_nodes().begin(), face->get_nodes().end());

  // Preallocate/emplace empty rows into the remaining face
  // adjacency tables to keep the mesh consistent.
  FaceEdges_.emplace_back_row(face->num_edges());
  FaceFaces_.emplace_back_row(/* dynamic value */);
  FaceCells_.emplace_back_row(2);

  return faceIndex;

} // Mesh<...>::EmplaceFace

CellIndex Mesh::EmplaceCell(std::unique_ptr<iElement>&& cell, CellMark cellMark) {

  CellIndex const cellIndex(NumCells_++);

  // Emplace the cell properties.
  CellMarks_.emplace_back(cellMark);

  CellShapes_.emplace_back(cell->get_shape());
  CellVolumes_.emplace_back(cell->get_length_or_area_or_volume());
  CellCenterPos_.emplace_back(cell->get_center_coords());

  // Fill the cell nodes.
  ((CsrTable<size_t, size_t>&) CellNodes_).emplace_back_row(cell->get_nodes().begin(), cell->get_nodes().end());

  // Preallocate rows into the remaining cell
  // adjacency tables to keep the mesh consistent.
  CellEdges_.emplace_back_row(cell->num_edges());
  CellFaces_.emplace_back_row(cell->num_faces());
  CellCells_.emplace_back_row(cell->num_faces());

  return cellIndex;

} // Mesh<...>::EmplaceCell

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

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
  InversePermutation(permutation.begin(), permutation.end(), inversePermutation.begin());

  auto const reorderFunc = [&inversePermutation](Index<Tag> index) {
    return Index<Tag>(index != npos ? inversePermutation[size_t(index)] : npos);
  };

  std::ranges::for_each(Nodes(), [&](NodeIndex nodeIndex) {
    std::ranges::transform(
      AdjacentElements_<Tag>(nodeIndex),
      std::begin(AdjacentElements_<Tag>(nodeIndex)), reorderFunc);
  });

  std::ranges::for_each(Edges(), [&](EdgeIndex edgeIndex) {
    std::ranges::transform(
      AdjacentElements_<Tag>(edgeIndex),
      std::begin(AdjacentElements_<Tag>(edgeIndex)), reorderFunc);
  });

  std::ranges::for_each(Faces(), [&](FaceIndex faceIndex) {
    std::ranges::transform(
      AdjacentElements_<Tag>(faceIndex),
      std::begin(AdjacentElements_<Tag>(faceIndex)), reorderFunc);
  });

  std::ranges::for_each(Cells(), [&](CellIndex cellIndex) {
    std::ranges::transform(
      AdjacentElements_<Tag>(cellIndex),
      std::begin(AdjacentElements_<Tag>(cellIndex)), reorderFunc);
  });

} // Mesh::FixPermutationAndAdjacency_

void Mesh::PermuteNodes(std::vector<size_t>&& nodePermutation) {

  /* Permute node properties and fix the adjacency tables. */
  FixPermutationAndAdjacency_<NodeTag_>(nodePermutation);
  permute_rows(nodePermutation.begin(), nodePermutation.end(),
               NodeNodes_, NodeEdges_, NodeFaces_, NodeCells_);
  permute_inplace(nodePermutation.begin(), nodePermutation.end(),
                  NodeMarks_.begin(), NodePos_.begin());

  /* Generate the node ranges. */
  NodeRanges_.clear();
  std::for_each(feathers::BeginNode(*this), feathers::EndNode(*this), [&](NodeRef node) {
    NodeRanges_.resize((size_t)node.Mark() + 2);
    NodeRanges_[node.Mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&)NodeRanges_;
  std::partial_sum(r.begin(), r.end(), r.begin());

} // Mesh<...>::PermuteNodes

void Mesh::PermuteEdges(std::vector<size_t>&& edgePermutation) {

  /* Permute edge properties and fix the adjacency tables. */
  FixPermutationAndAdjacency_<EdgeTag_>(edgePermutation);
  permute_rows(edgePermutation.begin(), edgePermutation.end(),
               EdgeNodes_, EdgeEdges_, EdgeFaces_, EdgeCells_);
  permute_inplace(edgePermutation.begin(), edgePermutation.end(),
                  EdgeMarks_.begin(), EdgeShapes_.begin(),
                  EdgeLens_.begin(), EdgeDirs_.begin());

  /* Generate the edge ranges. */
  EdgeRanges_.clear();
  std::for_each(feathers::BeginEdge(*this), feathers::EndEdge(*this), [&](EdgeRef edge) {
    EdgeRanges_.resize((size_t)edge.Mark() + 2);
    EdgeRanges_[edge.Mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&)EdgeRanges_;
  std::partial_sum(r.begin(), r.end(), r.begin());

} // Mesh<...>::PermuteEdges

void Mesh::PermuteFaces(std::vector<size_t>&& facePermutation) {

  /* Permute data. */
  FixPermutationAndAdjacency_<FaceTag_>(facePermutation);
  permute_rows(
    facePermutation.begin(), facePermutation.end(),
    FaceNodes_, FaceEdges_, FaceFaces_, FaceCells_);
  permute_inplace(
    facePermutation.begin(), facePermutation.end(),
    FaceMarks_.begin(), FaceShapes_.begin(),
    FaceAreas_.begin(), FaceNormals_.begin(), FaceCenterPos_.begin());

  /* Generate mark ranges. */
  FaceRanges_.clear();
  std::for_each(feathers::BeginFace(*this), feathers::EndFace(*this), [&](FaceRef face) {
    FaceRanges_.resize((size_t)face.Mark() + 2);
    FaceRanges_[face.Mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&)FaceRanges_;
  std::partial_sum(r.begin(), r.end(), r.begin());

} // Mesh<...>::PermuteFaces

void Mesh::PermuteCells(std::vector<size_t>&& cellPermutation) {

  /* Permute data. */
  FixPermutationAndAdjacency_<CellTag_>(cellPermutation);
  permute_rows(
    cellPermutation.begin(), cellPermutation.end(),
    CellNodes_, CellEdges_, CellFaces_, CellCells_);
  permute_inplace(
    cellPermutation.begin(), cellPermutation.end(),
    CellMarks_.begin(), CellShapes_.begin(),
    CellVolumes_.begin(), CellCenterPos_.begin());

  /* Generate mark ranges. */
  CellRanges_.clear();
  std::for_each(feathers::BeginCell(*this), feathers::EndCell(*this), [&](CellRef cell) {
    CellRanges_.resize((size_t)cell.Mark() + 2);
    CellRanges_[cell.Mark() + 1] += 1;
  });

  auto& r = (std::vector<size_t>&)CellRanges_;
  std::partial_sum(r.begin(), r.end(), r.begin());

} // Mesh<...>::PermuteCells

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

void Mesh::reorder_faces() {

  /** @todo Refactor me! */
  {
    std::vector<size_t> node_reordering(NumNodes());
    std::iota(node_reordering.begin(), node_reordering.end(), 0);
    PermuteNodes(std::move(node_reordering));
  }
  {
    std::vector<size_t> edge_reordering(NumEdges());
    std::iota(edge_reordering.begin(), edge_reordering.end(), 0);
    PermuteEdges(std::move(edge_reordering));
  }
  {
    std::vector<size_t> face_reordering(NumFaces());
    std::iota(face_reordering.begin(), face_reordering.end(), 0);
    PermuteFaces(std::move(face_reordering));
  }
  {
    std::vector<size_t> cell_reordering(NumCells());
    std::iota(cell_reordering.begin(), cell_reordering.end(), 0);
    PermuteCells(std::move(cell_reordering));
  }

} // Mesh<...>::PermuteFaces

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

void Mesh::FinalizeEdges_() {

  // ----------------------
  // An edge lookup table.
  // We assume that the edge can be uniquely identified by the set of its nodes.
  // ----------------------
  std::map<std::set<NodeIndex>, EdgeIndex> edgeLookupTable;

  // ----------------------
  // Add the existing edges to the lookup table.
  // ----------------------
  std::ranges::for_each(Edges(), [&](EdgeIndex edgeIndex) {
    std::set<NodeIndex> edgeLookupKey(
      std::begin(AdjacentNodes(edgeIndex)), std::end(AdjacentNodes(edgeIndex)));
    edgeLookupTable.emplace(std::move(edgeLookupKey), edgeIndex);
  });

  // ----------------------
  // For each face-to-edge adjacency table entry:
  // find the edge in the lookup table or emplace the new edge.
  // ----------------------
  std::for_each(feathers::BeginFace(*this), feathers::EndFace(*this), [&](FaceMutableRef face) {
    tElementDescList edgesDesc = face.get_element_object()->get_edges_desc();
    for (size_t edgeLocal = 0; edgeLocal < face.NumEdges(); ++edgeLocal) {
      EdgeIndex& edgeIndex = AdjacentEdges(face)[edgeLocal];
      if (edgeIndex != npos) {
        continue;
      }

      // Create the face or add current cell to the adjacency list.
      sElementDesc& edgeDesc = edgesDesc[edgeLocal];
      std::set<NodeIndex> edgeLookupKey(
        edgeDesc.node_indices.begin(), edgeDesc.node_indices.end());
      if (edgeLookupTable.count(edgeLookupKey) == 0) {
        // Create a brand-new edge.
        edgeIndex = EmplaceEdge(std::move(edgeDesc), (EdgeMark)face.Mark());
        edgeLookupTable.emplace(edgeLookupKey, edgeIndex);
      } else {
        // Edge exists.
        edgeIndex = edgeLookupTable[edgeLookupKey];
      }
    }
  });

  // ----------------------
  // Check faces:
  // each face should be connected to all edges.
  // ----------------------
  std::ranges::for_each(Faces(), [&](FaceIndex faceIndex) {
    StormEnsure("Face-edge connectivity is broken" &&
      std::ranges::all_of(AdjacentEdges(faceIndex), is_not_npos<EdgeIndex>));
  });

} // Mesh<...>::FinalizeEdges_

void Mesh::FinalizeFaces_() {

  // ----------------------
  // A face lookup table.
  // We assume that the face can be uniquely identified by the set of its nodes.
  // ----------------------
  std::map<std::set<NodeIndex>, FaceIndex> faceLookupTable;

  // ----------------------
  // Add the existing faces to the lookup table.
  // ----------------------
  std::ranges::for_each(Faces(), [&](FaceIndex faceIndex) {
    std::set<NodeIndex> faceLookupKey(
      std::begin(AdjacentNodes(faceIndex)), std::end(AdjacentNodes(faceIndex)));
    faceLookupTable.emplace(std::move(faceLookupKey), faceIndex);
  });

  // ----------------------
  // For each cell-to-face adjacency table entry:
  // find the face in the lookup table or emplace the new face.
  // Also fill the face-to-cell adjacency table.
  // ----------------------
  std::for_each(feathers::BeginCell(*this), feathers::EndCell(*this), [&](CellMutableRef cell) {
    tElementDescList facesDesc = cell.get_element_object()->get_faces_desc();
    for (size_t faceLocal = 0; faceLocal < cell.NumFaces(); ++faceLocal) {
      FaceIndex& faceIndex = AdjacentFaces(cell)[faceLocal];
      if (faceIndex != npos) {
        continue;
      }

      /* Create the face or add current cell to the adjacency list. */
      sElementDesc& faceDesc = facesDesc[faceLocal];
      std::set<NodeIndex> faceLookupKey(
        faceDesc.node_indices.begin(), faceDesc.node_indices.end());
      if (faceLookupTable.count(faceLookupKey) == 0) {
        /* Create a brand-new face.
         * Assign the current cell as the inner one. */
        faceIndex = EmplaceFace(std::move(faceDesc));
        faceLookupTable.emplace(std::move(faceLookupKey), faceIndex);
        AdjacentCells(faceIndex)[FaceInnerCell_] = cell;
      } else {
        /* Face exists.
         * Determine orientation and link it. */
        faceIndex = faceLookupTable[faceLookupKey];
        FaceMutableRef face(*this, faceIndex);
        if (std::equal(
            std::begin(AdjacentNodes(face)), std::end(AdjacentNodes(face)), faceDesc.node_indices.begin())) {
          /* Face node order matches the order
           * in the face description: face is inner. */
          FEATHERS_ASSERT(AdjacentCells(face)[FaceInnerCell_] == npos);
          AdjacentCells(face)[FaceInnerCell_] = cell;
        } else {
          /* Otherwise, the face is outer. */
          FEATHERS_ASSERT(AdjacentCells(face)[FaceOuterCell_] == npos);
          AdjacentCells(face)[FaceOuterCell_] = cell;
        }
      }
    }
  });

  // ----------------------
  // Check cells:
  // each cell should be connected to all faces.
  // ----------------------
  std::ranges::for_each(Cells(), [&](CellIndex cellIndex) {
    StormEnsure("Cell-face connectivity is broken" &&
      std::ranges::all_of(AdjacentFaces(cellIndex), is_not_npos<FaceIndex>));
  });

  // ----------------------
  // Check faces:
  // each internal face should be connected to two cells;
  // each boundary face should be connected to at least one cell.
  // ----------------------
  std::ranges::for_each(Faces(), [&](FaceIndex faceIndex) {
    if (Mark(faceIndex) == 0) {
      StormEnsure("Interior face-cell connectivity is broken" &&
        std::ranges::all_of(AdjacentCells(faceIndex), is_not_npos<CellIndex>));
    } else {
      StormEnsure("Boundary face-cell connectivity is broken" &&
        std::ranges::count_if(AdjacentCells(faceIndex), is_not_npos<CellIndex>) == 1);
    }
  });

} // Mesh<...>::FinalizeFaces_

/* A node and edge flip table for various face types. */
static const std::map<eShape, std::pair<std::vector<size_t>, std::vector<size_t>>>
  g_face_shape_to_nodes_and_edges_flip {
  /* 1D faces. */
  { eShape::node, { {0}, {0} } },
  /* 2D faces. */
  { eShape::segment_2, { {1, 0}, {1, 0} } },
  /* 3D faces. */
  { eShape::triangle_3, { {0, 2, 1}, {0, 2, 1} } },
  { eShape::quadrangle_4, { {0, 3, 2, 1}, {0, 3, 2, 1} } },
};

/**
 * Generate boundary cells to complete face connectivity.
 */
void Mesh::generate_boundary_cells() {

  std::for_each(feathers::BeginFace(*this), feathers::EndFace(*this), [&](FaceMutableRef face) {
    if (face.Mark() == 0) {
      return;
    }

    /* Boundary faces should be oriented outwards from the mesh. */
    if (face.InnerCell() == npos) {
      /* Flip normal and cell connectivity. */
      face.SetNormal(-face.Normal());
      std::swap(AdjacentCells(face)[FaceInnerCell_],
                AdjacentCells(face)[FaceOuterCell_]);
      /* Flip node and edge connectivity. */
      std::vector<size_t> node_permutation, edge_permutation;
      std::tie(node_permutation, edge_permutation) =
        g_face_shape_to_nodes_and_edges_flip.at(face.Shape());
      permute_inplace(
        node_permutation.begin(), node_permutation.end(), std::begin(AdjacentNodes(face)));
      permute_inplace(
        edge_permutation.begin(), edge_permutation.end(), std::begin(AdjacentNodes(face)));
    }
    CellRef cell = face.InnerCell();

    /* Generate the boundary cell: reflect a connected interior cell. */
    std::vector<size_t> ghost_cell_nodes;
#if 0
    ghost_cell_nodes.assign(face.begin_node(), face.end_node());
#endif
#if 1
    cell.ForEachNode([&](NodeRef node) {
      if (std::ranges::find(AdjacentNodes(face), (NodeIndex)node) == std::end(AdjacentNodes(face))) {
        /* Reflect an interior cell node. */
        // TODO: face normals are not computed here!
        // TODO: https://glm.g-truc.net/0.9.5/api/a00157.html#gab63646fc36b81cf69d3ce123a72f76f2
        vec3_t node_coords = node.Pos();
        const vec3_t delta = node_coords - face.CenterPos();
        node_coords -= 2.0 * glm::dot(delta, face.Normal()) * face.Normal();
        ghost_cell_nodes.push_back(
          (size_t)EmplaceNode(node_coords, (NodeMark)face.Mark()));
      } else {
        /* Insert a boundary face node. */
        ghost_cell_nodes.push_back((size_t)node);
      }
    });
#endif
    /* Insert the boundary cell. */
    // TODO:
    const CellIndex boundaryCellIndex =
      EmplaceCell({cell.Shape(), ghost_cell_nodes}, (CellMark)face.Mark());
    AdjacentFaces(boundaryCellIndex)[0] = face;
#if 0
    Cell& boundary_cell = get_cell(boundaryCellIndex);
        while (boundary_cell._num_faces() != 1) {
            boundary_cell.erase_face(0);
        }
        boundary_cell.begin_face()[0] = face;
#endif

    StormAssert(AdjacentCells(face)[FaceOuterCell_] == npos);
    AdjacentCells(face)[FaceOuterCell_] = boundaryCellIndex;

  });
} // Mesh<...>::generate_boundary_cells

} // namespace feathers
