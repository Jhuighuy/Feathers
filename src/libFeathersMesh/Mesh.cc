/**
 *    ______             __     __  _____ _____
 *   / __/ /____ _____  / /__  /  |/  / // / _ \
 *  _\ \/  '_/ // / _ \/  '_/ / /|_/ / _  / // /
 * /___/_/\_\\_,_/_//_/_/\_\ /_/  /_/_//_/____/
 *
 * Copyright (c) 2019 Oleg Butakov
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "Mesh.hh"
#include <libFeathersMesh/Shape.hh>
#include <libSkunkMisc/SkunkMiscParallel.hh>
#include <libSkunkMisc/SkunkMiscReordering.hh>

#include <set>
#include <map>
#include <fstream>

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

/** @todo Refactor me! */
#define assert1(x) do { if(!(x)) {throw std::runtime_error("hui");} } while(false)

bool cMesh::read_triangle(const char *path) {
    std::string line;

    std::ifstream node_file(path + std::string("node"));
    FEATHERS_ENSURE(node_file.is_open());
    uint_t num_nodes = 0;
    node_file >> num_nodes >> m_dim;
    std::getline(node_file, line);
    for (uint_t i = 0; i < num_nodes; ++i) {
        uint_t node_index = 0;
        vec3_t node_pos;
        node_file >> node_index >> node_pos.x >> node_pos.y;
        std::getline(node_file, line);
        FEATHERS_ENSURE(node_index == insert_node(mesh_node1_t(), node_pos));
    }

    std::ifstream face_file(path + std::string("edge"));
    assert1(face_file.is_open());
    uint_t num_faces = 0;
    face_file >> num_faces;
    std::getline(face_file, line);
    for (uint_t i = 0; i < num_faces; ++i) {
        uint_t face_index = 0;
        std::array<uint_t, 2> face_nodes{npos};
        uint_t mark = 0;
        face_file >> face_index >> face_nodes[0] >> face_nodes[1] >> mark;
        FEATHERS_ENSURE(
            face_index == insert_face(mesh_face_segment2_t(face_nodes.begin(), face_nodes.end()), mark));
        std::getline(face_file, line);
    }

    std::ifstream cell_file(path + std::string("ele"));
    assert1(cell_file.is_open());
    uint_t num_cells = 0;
    cell_file >> num_cells;
    std::getline(cell_file, line);
    for (uint_t i = 0; i < num_cells; ++i) {
        uint_t cell_index = 0;
        std::array<uint_t, 3> cell_nodes{npos};
        cell_file >> cell_index >> cell_nodes[0] >> cell_nodes[1] >> cell_nodes[2];
        FEATHERS_ENSURE(
            cell_index == insert_cell(mesh_cell_triangle3_t(cell_nodes.begin(), cell_nodes.end())));
        std::getline(cell_file, line);
    }

    generate_faces();
    generate_boundary_cells();
    reorder_faces();
    compute_all_shape_properties();
    return true;
}   // cMesh::read_triangle

bool cMesh::read_tetgen(const char *path) {
    std::string line;

    std::ifstream node_file(path + std::string("node"));
    FEATHERS_ENSURE(node_file.is_open());
    uint_t num_nodes = 0;
    node_file >> num_nodes >> m_dim;
    std::getline(node_file, line);
    for (uint_t i = 0; i < num_nodes; ++i) {
        uint_t node_index = 0;
        vec3_t node_pos;
        node_file >> node_index >> node_pos.x >> node_pos.y >> node_pos.z;
        FEATHERS_ENSURE(
            node_index == insert_node(mesh_node1_t(), node_pos));
        std::getline(node_file, line);
    }

    std::ifstream face_file(path + std::string("face"));
    assert1(face_file.is_open());
    uint_t num_faces = 0;
    face_file >> num_faces;
    std::getline(face_file, line);
    for (uint_t i = 0; i < num_faces; ++i) {
        uint_t face_index = 0;
        std::array<uint_t, 3> face_nodes{npos};
        uint_t mark = 0;
        face_file >> face_index >> face_nodes[0] >> face_nodes[1] >> face_nodes[2] >> mark;
        FEATHERS_ENSURE(
            face_index == insert_face(mesh_face_triangle3_t(face_nodes.begin(), face_nodes.end()), mark));
        std::getline(face_file, line);
    }

    std::ifstream cell_file(path + std::string("ele"));
    FEATHERS_ENSURE(cell_file.is_open());
    uint_t num_cells = 0;
    cell_file >> num_cells;
    std::getline(cell_file, line);
    for (uint_t i = 0; i < num_cells; ++i) {
        uint_t cell_index = 0;
        std::array<uint_t, 4> cell_nodes{npos};
        cell_file >> cell_index >> cell_nodes[0] >> cell_nodes[1] >> cell_nodes[2] >> cell_nodes[3];
        FEATHERS_ENSURE(
            cell_index == insert_cell(mesh_cell_tetrahedron4_t(cell_nodes.begin(), cell_nodes.end())));
        std::getline(cell_file, line);
    }

    generate_faces();
    generate_boundary_cells();
    reorder_faces();
    compute_all_shape_properties();
    return true;
}   // cMesh::read_tetgen

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/** Compute edge shape properties. */
void cMesh::compute_edge_shape_properties() {
    for_each_edge(*this, [&](tEdgeMutableIter edge) {
        iShapePtr edge_shape = edge.get_shape_ptr();
        edge.set_length(edge_shape->get_length_or_area_or_volume());
        edge.set_direction(edge_shape->get_direction());
    });
}   // cMesh::compute_edge_shape_properties

/** Compute face shape properties. */
void cMesh::compute_face_shape_properties() {
    for_each_face(*this, [&](tFaceMutableIter face) {
        iShapePtr face_shape = face.get_shape_ptr();
        face.set_area(face_shape->get_length_or_area_or_volume());
        face.set_normal(face_shape->get_normal());
        face.set_center_coords(face_shape->get_center_coords());
    });
}   // cMesh::compute_face_shape_properties

/** Compute cell shape properties. */
void cMesh::compute_cell_shape_properties() {
    for_each_cell(*this, [&](tCellMutableIter cell) {
        iShapePtr cell_shape = cell.get_shape_ptr();
        cell.set_volume(cell_shape->get_length_or_area_or_volume());
        cell.set_center_coords(cell_shape->get_center_coords());
    });
}   // cMesh::compute_cell_shape_properties

/**
 * Get minimal edge length.
 */
real_t cMesh::get_min_edge_length() const {
    // TODO: refactor with mesh iterators!
    return
        m_min_edge_length > 0.0 ? m_min_edge_length :
        (m_min_edge_length = for_range_min(0u, num_edges(), huge, [&](uint_t edge_index) {
            return get_edge_length(edge_index);
        }));
}   // cMesh::get_min_edge_length
/**
 * Get maximal edge length.
 */
real_t cMesh::get_max_edge_length() const {
    // TODO: refactor with mesh iterators!
    return
        m_max_edge_length > 0.0 ? m_max_edge_length :
        (m_max_edge_length = for_range_max(0u, num_edges(), 0.0, [&](uint_t edge_index) {
            return get_edge_length(edge_index);
        }));
}   // cMesh::get_max_edge_length

/**
 * Get minimal face area.
 */
real_t cMesh::get_min_face_area() const {
    // TODO: refactor with mesh iterators!
    return
        m_min_face_area > 0.0 ? m_min_face_area :
        (m_min_face_area = for_range_min(0u, num_faces(), huge, [&](uint_t face_index) {
            return get_face_area(face_index);
        }));
}   // cMesh::get_min_face_area
/**
 * Get maximal face area.
 */
real_t cMesh::get_max_face_area() const {
    // TODO: refactor with mesh iterators!
    return
        m_max_face_area > 0.0 ? m_max_face_area :
        (m_max_face_area = for_range_max(0u, num_faces(), 0.0, [&](uint_t face_index) {
            return get_face_area(face_index);
        }));
}   // cMesh::get_max_face_area

/**
 * Get minimal cell volume.
 */
real_t cMesh::get_min_cell_volume() const {
    // TODO: refactor with mesh iterators!
    return
        m_min_cell_volume > 0.0 ? m_min_cell_volume :
        (m_min_cell_volume = for_range_min(0u, num_cells(), huge, [&](uint_t cell_index) {
            return get_cell_volume(cell_index);
        }));
}   // cMesh::get_min_cell_volume
/**
 * Get maximal cell volume.
 */
real_t cMesh::get_max_cell_volume() const {
    // TODO: refactor with mesh iterators!
    return
        m_max_cell_volume > 0.0 ? m_max_cell_volume :
        (m_max_cell_volume = for_range_max(0u, num_cells(), 0.0, [&](uint_t cell_index) {
            return get_cell_volume(cell_index);
        }));
}   // cMesh::get_max_cell_volume

/**
 * Compute mesh orthogonality.
 * Mesh orthogonality is defined as a...
 */
real_t cMesh::get_orthogonality() const {
    // TODO: refactor with mesh iterators!
    return for_range_min(begin_face(0), end_face(0), huge, [&](uint_t face_index) {
        tFaceIter face(this, face_index);
        const vec3_t delta =
            face.get_outer_cell().get_center_coords() -
            face.get_inner_cell().get_center_coords();
        const real_t orth = glm::dot(face.get_normal(), glm::normalize(delta));
        return orth;
    });
    return huge;
}   // cMesh::get_orthogonality

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/** Insert a new node into the mesh.
 ** @param node_ Edge nodes.
 ** @returns Index of the inserted node. */
uint_t cMesh::insert_node(const mesh_node1_t& node_, vec3_t p, uint_t mark) {
    /** @todo Refactor me! */
    const uint_t node_index = allocate_node_(node_);
    set_mark(eNodeTag, node_index, mark);
    set_node_coords(node_index, p);
    const vec3_t& node_geometry = p;
    if (m_dim <= 2) {
        FEATHERS_ASSERT(node_geometry.z == 0.0);
        if (m_dim == 1) {
            FEATHERS_ASSERT(node_geometry.y == 0.0);
        }
    }
    return node_index;
}   // cMesh::insert_node

/** Insert a new edge into the mesh.
 ** @param edge_ Edge nodes.
 ** @returns Index of the inserted edge. */
/** @{ */
uint_t cMesh::insert_edge(const mesh_edge_node1_t& edge_, uint_t mark) {
    const uint_t edge_index = allocate_edge_(edge_);
    set_mark(eEdgeTag, edge_index, mark);
    set_shape(eEdgeTag, edge_index, eShape::node);
    return edge_index;
}   // cMesh::insert_edge
uint_t cMesh::insert_edge(const mesh_edge_segment2_t& edge_, uint_t mark) {
    const uint_t edge_index = allocate_edge_(edge_);
    set_mark(eEdgeTag, edge_index, mark);
    set_shape(eEdgeTag, edge_index, eShape::segment_2);
    return edge_index;
}   // cMesh::insert_edge
/** @} */

/** Insert a new edge into the mesh.
 ** Type of the edge is detected by the set of nodes and dimension.
 ** @param edge_nodes Set of the edge nodes.
 ** @param mark Boundary mark.
 ** @param dim Overridden dimension.
 ** @returns Index of the inserted edge. */
uint_t cMesh::insert_edge(const std::vector<uint_t>& edge_nodes, uint_t mark, uint_t dim) {
    if (dim == 0) {
        dim = m_dim;
    }
    FEATHERS_ASSERT(1 <= dim && dim <= 3);
    if (dim <= 2) {
        switch (edge_nodes.size()) {
        /* Dummy 1D/2D edge. */
        case 1:
            return insert_edge(mesh_edge_node1_t(edge_nodes.cbegin(), edge_nodes.cend()), mark);
        default:;
        }
    } else if (dim == 3) {
        switch (edge_nodes.size()) {
        /* Segment edge. */
        case 2:
            return insert_edge(mesh_edge_segment2_t(edge_nodes.cbegin(), edge_nodes.cend()), mark);
        default:;
        }
    }
    FEATHERS_ASSERT_FALSE("Invalid edge nodes: cell type cannot be uniquely matched.");
}   // cMesh::insert_edge

/** Insert a new face into the mesh.
 ** @param face_ Face nodes.
 ** @returns Index of the inserted face. */
/** @{ */
uint_t cMesh::insert_face(const mesh_face_node1_t& face_, uint_t mark) {
    const uint_t face_index = allocate_face_(face_);
    set_mark(eFaceTag, face_index, mark);
    set_shape(eFaceTag, face_index, eShape::node);
    return face_index;
}   // cMesh::insert_face
uint_t cMesh::insert_face(const mesh_face_segment2_t& face_, uint_t mark) {
    const uint_t face_index = allocate_face_(face_);
    set_mark(eFaceTag, face_index, mark);
    set_shape(eFaceTag, face_index, eShape::segment_2);
    return face_index;
}   // cMesh::insert_face
uint_t cMesh::insert_face(const mesh_face_triangle3_t& face_, uint_t mark) {
    const uint_t face_index = allocate_face_(face_);
    set_mark(eFaceTag, face_index, mark);
    set_shape(eFaceTag, face_index, eShape::triangle_3);
    return face_index;
}   // cMesh::insert_face
uint_t cMesh::insert_face(const mesh_face_quadrangle4_t& face_, uint_t mark) {
    const uint_t face_index = allocate_face_(face_);
    set_mark(eFaceTag, face_index, mark);
    set_shape(eFaceTag, face_index, eShape::quadrangle_4);
    return face_index;
}   // cMesh::insert_face
/** @} */

/** Insert a new face into the mesh.
 ** Type of the face is detected by the set of nodes and dimension.
 ** @param face_nodes Set of the face nodes.
 ** @param mark Boundary mark.
 ** @param dim Overridden dimension.
 ** @returns Index of the inserted face. */
uint_t cMesh::insert_face(const std::vector<uint_t>& face_nodes, uint_t mark, uint_t dim) {
    if (dim == 0) {
        dim = m_dim;
    }
    FEATHERS_ASSERT(1 <= dim && dim <= 3);
    if (dim == 1) {
        switch (face_nodes.size()) {
        /* Dummy 1D face. */
        case 1:
            return insert_face(mesh_face_node1_t(face_nodes.cbegin(), face_nodes.cend()), mark);
        default:;
        }
    } else if (dim == 2) {
        switch (face_nodes.size()) {
        /* Segment face. */
        case 2:
            return insert_face(mesh_face_segment2_t(face_nodes.cbegin(), face_nodes.cend()), mark);
        default:;
        }
    } else if (dim == 3) {
        switch (face_nodes.size()) {
        /* Triangular face. */
        case 3:
            return insert_face(mesh_face_triangle3_t(face_nodes.cbegin(), face_nodes.cend()), mark);
        /* Quadrangular face. */
        case 4:
            return insert_face(mesh_face_quadrangle4_t(face_nodes.cbegin(), face_nodes.cend()), mark);
        default:;
        }
    }
    FEATHERS_ASSERT_FALSE("Invalid face nodes: face type cannot be uniquely matched.");
}   // cMesh::insert_face

/** Insert a new cell into the mesh.
 ** @param cell_ Cell nodes.
 ** @returns Index of the inserted face. */
/** @{ */
uint_t cMesh::insert_cell(const mesh_cell_segment2_t& cell_, uint_t mark) {
    const uint_t cell_index = allocate_cell_(cell_);
    set_mark(eCellTag, cell_index, mark);
    set_shape(eCellTag, cell_index, eShape::segment_2);
    return cell_index;
}   // cMesh::insert_cell
uint_t cMesh::insert_cell(const mesh_cell_triangle3_t& cell_, uint_t mark) {
    const uint_t cell_index = allocate_cell_(cell_);
    set_mark(eCellTag, cell_index, mark);
    set_shape(eCellTag, cell_index, eShape::triangle_3);
    return cell_index;
}   // cMesh::insert_cell
uint_t cMesh::insert_cell(const mesh_cell_quadrangle4_t& cell_, uint_t mark) {
    const uint_t cell_index = allocate_cell_(cell_);
    set_mark(eCellTag, cell_index, mark);
    set_shape(eCellTag, cell_index, eShape::quadrangle_4);
    return cell_index;
}   // cMesh::insert_cell
uint_t cMesh::insert_cell(const mesh_cell_tetrahedron4_t& cell_, uint_t mark) {
    const uint_t cell_index = allocate_cell_(cell_);
    set_mark(eCellTag, cell_index, mark);
    set_shape(eCellTag, cell_index, eShape::tetrahedron_4);
    return cell_index;
}   // cMesh::insert_face
uint_t cMesh::insert_cell(const mesh_cell_pyramid5_t& cell_, uint_t mark) {
    const uint_t cell_index = allocate_cell_(cell_);
    set_mark(eCellTag, cell_index, mark);
    set_shape(eCellTag, cell_index, eShape::pyramid_5);
    return cell_index;
}
uint_t cMesh::insert_cell(const mesh_cell_pentahedron6_t& cell_, uint_t mark) {
    const uint_t cell_index = allocate_cell_(cell_);
    set_mark(eCellTag, cell_index, mark);
    set_shape(eCellTag, cell_index, eShape::pentahedron_6);
    return cell_index;
}   // cMesh::insert_face
uint_t cMesh::insert_cell(const feathers::mesh_cell_hexahedron8_t& cell_, uint_t mark) {
    const uint_t cell_index = allocate_cell_(cell_);
    set_mark(eCellTag, cell_index, mark);
    set_shape(eCellTag, cell_index, eShape::hexahedron_8);
    return cell_index;
}
/** @} */

/** Insert a new cell into the mesh.
 ** Type of the cell is detected by the set of nodes and dimension.
 ** @param cell_nodes Set of the cell nodes.
 ** @param mark Boundary mark.
 ** @param dim Overridden dimension.
 ** @returns Index of the inserted cell. */
uint_t cMesh::insert_cell(const std::vector<uint_t>& cell_nodes, uint_t mark, uint_t dim) {
    if (dim == 0) {
        dim = m_dim;
    }
    FEATHERS_ASSERT(1 <= dim && dim <= 3);
    if (dim == 1) {
        switch (cell_nodes.size()) {
        /* Segment cell. */
        case 2:
            return insert_cell(mesh_cell_segment2_t(cell_nodes.cbegin(), cell_nodes.cend()), mark);
        default:;
        }
    } else if (dim == 2) {
        switch (cell_nodes.size()) {
        /* Triangular cell. */
        case 3:
            return insert_cell(mesh_cell_triangle3_t(cell_nodes.cbegin(), cell_nodes.cend()), mark);
        /* Quadrangular cell. */
        case 4:
            return insert_cell(mesh_cell_quadrangle4_t(cell_nodes.cbegin(), cell_nodes.cend()), mark);
        default:;
        }
    } else if (dim == 3) {
        switch (cell_nodes.size()) {
        /* Tetrahedral cell. */
        case 4:
            return insert_cell(mesh_cell_tetrahedron4_t(cell_nodes.cbegin(), cell_nodes.cend()), mark);
        /* Pyramidal cell. */
        case 5:
            return insert_cell(mesh_cell_pyramid5_t(cell_nodes.cbegin(), cell_nodes.cend()), mark);
        /* Pentahedral cell. */
        case 6:
            return insert_cell(mesh_cell_pentahedron6_t(cell_nodes.cbegin(), cell_nodes.cend()), mark);
        /* Hexahedral cell. */
        case 8:
            return insert_cell(mesh_cell_hexahedron8_t(cell_nodes.cbegin(), cell_nodes.cend()), mark);
        default:;
        }
    }
    FEATHERS_ASSERT(!"Invalid cell nodes: cell type cannot be uniquely matched.");
}   // cMesh::insert_cell

/** @internal
 ** Allocate and insert the new element into the storage. */
template<typename mesh_elem_struct_t>
SKUNK_INLINE uint_t allocate_element_(const mesh_elem_struct_t& elem, size_t alignment,
                                      std::vector<byte_t>& elem_storage, std::vector<size_t>& elem_offsets) {
    FEATHERS_ASSERT(alignment > 0);
    /* Align size of the element. */
    size_t delta_elem = sizeof(elem);
    delta_elem += (alignment - delta_elem%alignment)%alignment;
    /* Allocate and initialize a new element. */
    elem_offsets.push_back(elem_storage.size());
    elem_storage.resize(elem_storage.size() + delta_elem);
    new (elem_storage.data() + elem_offsets.back()) mesh_elem_struct_t(elem);
    const auto elem_index = static_cast<uint_t>(elem_offsets.size()) - 1;
    return elem_index;
}   // allocate_element_
/** @internal
 ** Allocate and insert a new node into the mesh. */
template<typename mesh_node_struct_t>
SKUNK_INLINE uint_t cMesh::allocate_node_(const mesh_node_struct_t& node, size_t alignment) {
    m_node_marks.emplace_back();
    m_node_coords.emplace_back();
    return allocate_element_(node, alignment, m_node_storage, m_node_offsets);
}   // cMesh::allocate_node_
/** @internal
 ** Allocate and insert a new edge into the mesh. */
template<typename mesh_edge_struct_t>
SKUNK_INLINE uint_t cMesh::allocate_edge_(const mesh_edge_struct_t& edge, size_t alignment) {
    m_edge_marks.emplace_back();
    m_edge_shapes.emplace_back();
    m_edge_lengths.emplace_back();
    m_edge_directions.emplace_back();
    return allocate_element_(edge, alignment, m_edge_storage, m_edge_offsets);
}   // cMesh::allocate_edge_
/** @internal
 ** Allocate and insert a new face into the mesh. */
template<typename mesh_face_struct_t>
SKUNK_INLINE uint_t cMesh::allocate_face_(const mesh_face_struct_t& face, size_t alignment) {
    m_face_marks.emplace_back();
    m_face_shapes.emplace_back();
    m_face_areas.emplace_back();
    m_face_normals.emplace_back();
    m_face_center_coords.emplace_back();
    return allocate_element_(face, alignment, m_face_storage, m_face_offsets);
}   // cMesh::allocate_face_
/** @internal
 ** Allocate and insert a new cell into the mesh. */
template<typename mesh_cell_struct_t>
SKUNK_INLINE uint_t cMesh::allocate_cell_(const mesh_cell_struct_t& cell, size_t alignment) {
    m_cell_marks.emplace_back();
    m_cell_shapes.emplace_back();
    m_cell_volumes.emplace_back();
    m_cell_center_coords.emplace_back();
    return allocate_element_(cell, alignment, m_cell_storage, m_cell_offsets);
}   // cMesh::allocate_cell_

/**************************************************************************/
/**************************************************************************/

/** @internal
 ** Swap bytes of two elements. */
SKUNK_INLINE static void swap_elem_bytes_(uint_t first_elem_index, uint_t second_elem_index,
                                          std::vector<byte_t>& elem_storage, std::vector<size_t>& elem_offsets) {
    /* Find ranges for element storage. */
    FEATHERS_ASSERT(first_elem_index <= elem_offsets.size());
    const size_t first_elem_beg = elem_offsets[first_elem_index];
    const size_t first_elem_end = elem_offsets.size() > first_elem_index + 1 ?
                                  elem_offsets[first_elem_index + 1] : elem_storage.size();

    /* Find ranges for new element storage. */
    FEATHERS_ASSERT(second_elem_index <= elem_offsets.size());
    const size_t second_elem_beg = elem_offsets[second_elem_index];
    const size_t second_elem_end = elem_offsets.size() > second_elem_index + 1 ?
                                   elem_offsets[second_elem_index + 1] : elem_storage.size();

    const size_t first_elem_range = first_elem_end - first_elem_beg;
    const size_t second_elem_range = second_elem_end - second_elem_beg;
    if (first_elem_range == second_elem_range) {
        /* Identical sizes -- fast inplace swap. */
        std::swap_ranges(elem_storage.begin() + first_elem_beg,
                         elem_storage.begin() + first_elem_end, elem_storage.begin() + second_elem_beg);
    } else {
        /* Identical sizes -- slow inplace rotation. */
        const auto iter = std::rotate(elem_storage.begin() + first_elem_beg,
                                      elem_storage.begin() + second_elem_beg, elem_storage.begin() + second_elem_end);
        std::rotate(iter, iter + first_elem_range, elem_storage.begin() + second_elem_end);
        const ptrdiff_t range = second_elem_range - first_elem_range;
        for (; first_elem_index <= second_elem_index; ++first_elem_index) {
            FEATHERS_ASSERT(elem_offsets[first_elem_index] >= first_elem_beg);
            elem_offsets[first_elem_index] += range;
        }
    }
}   // swap_elem_bytes_

class tReorderFunc {
private:
    const std::vector<uint_t>& m_indices;
public:
    explicit tReorderFunc(const std::vector<uint_t>& indices)
        : m_indices(indices) {
    }
    void operator()(uint_t& index) const {
        if (index != npos) {
            index = m_indices.at(static_cast<size_t>(index));
        }
    }
};  // class tReorderFunc

/** Change order of all nodes. */
void cMesh::reorder_nodes(const std::vector<uint_t>& node_permutation) {
    {
        /* Reconnect edge, face and cells with nodes. */
        std::vector<uint_t> node_indices(node_permutation.size());
        convert_permutation_to_indices(
            node_permutation.begin(), node_permutation.end(), node_indices.begin());
        for_range(0u, num_edges(), [&](uint_t edge_index) {
            std::for_each(begin_adjacent_node(eEdgeTag, edge_index),
                          end_adjacent_node(eEdgeTag, edge_index), tReorderFunc(node_indices));
        });
        for_range(0u, num_faces(), [&](uint_t face_index) {
            std::for_each(begin_adjacent_node(eFaceTag, face_index),
                          end_adjacent_node(eFaceTag, face_index), tReorderFunc(node_indices));
        });
        for_range(0u, num_cells(), [&](uint_t cell_index) {
            std::for_each(begin_adjacent_node(eCellTag, cell_index),
                          end_adjacent_node(eCellTag, cell_index), tReorderFunc(node_indices));
        });
    }
    /* Permute node attributes. */
    reorder_nodes_bytes_(node_permutation);
    permute_inplace(node_permutation, m_node_marks.begin());
    permute_inplace(node_permutation, m_node_coords.begin());
}   // cMesh::reorder_nodes

/** Change order of all edges. */
void cMesh::reorder_edges(const std::vector<uint_t>& edge_permutation) {
    {
        /* Reconnect faces with edges. */
        std::vector<uint_t> edge_indices(edge_permutation.size());
        convert_permutation_to_indices(
            edge_permutation.begin(), edge_permutation.end(), edge_indices.begin());
        for_range(0u, num_faces(), [&](uint_t face_index) {
            std::for_each(begin_adjacent_edge(eFaceTag, face_index),
                          end_adjacent_edge(eFaceTag, face_index), tReorderFunc(edge_indices));
        });
    }
    /* Permute edge attributes. */
    reorder_edges_bytes_(edge_permutation);
    permute_inplace(edge_permutation, m_edge_marks.begin());
    permute_inplace(edge_permutation, m_edge_shapes.begin());
    permute_inplace(edge_permutation, m_edge_lengths.begin());
    permute_inplace(edge_permutation, m_edge_directions.begin());
}   // cMesh::reorder_edges
/** Change order of all faces. */
void cMesh::reorder_faces(const std::vector<uint_t>& face_permutation) {
    {
        /* Reconnect cells with faces. */
        std::vector<uint_t> face_indices(face_permutation.size());
        convert_permutation_to_indices(
            face_permutation.begin(), face_permutation.end(), face_indices.begin());
        for_range(0u, num_cells(), [&](uint_t cell_index) {
            std::for_each(begin_adjacent_face(eCellTag, cell_index),
                          end_adjacent_face(eCellTag, cell_index), tReorderFunc(face_indices));
        });
    }
    /* Permute face attributes. */
    reorder_faces_bytes_(face_permutation);
    permute_inplace(face_permutation, m_face_marks.begin());
    permute_inplace(face_permutation, m_face_shapes.begin());
    permute_inplace(face_permutation, m_face_areas.begin());
    permute_inplace(face_permutation, m_face_normals.begin());
    permute_inplace(face_permutation, m_face_center_coords.begin());
}   // cMesh::reorder_faces
/** Change order of all cells. */
void cMesh::reorder_cells(const std::vector<uint_t>& cell_permutation) {
    {
        /* Reconnect faces with cells. */
        std::vector<uint_t> cell_indices(cell_permutation.size());
        convert_permutation_to_indices(
            cell_permutation.begin(), cell_permutation.end(), cell_indices.begin());
        for_range(0u, num_faces(), [&](uint_t face_index) {
            std::for_each(begin_adjacent_node(eFaceTag, face_index),
                          end_adjacent_node(eFaceTag, face_index), tReorderFunc(cell_permutation));
        });
    }
    /* Permute cell attributes. */
    reorder_cells_bytes_(cell_permutation);
    permute_inplace(cell_permutation, m_cell_marks.begin());
    permute_inplace(cell_permutation, m_cell_shapes.begin());
    permute_inplace(cell_permutation, m_cell_volumes.begin());
    permute_inplace(cell_permutation, m_cell_center_coords.begin());
}   // cMesh::reorder_cells

/** @internal
 ** Swap bytes of two elements. */
SKUNK_INLINE static void reorder_elems_bytes_(std::vector<uint_t> elem_reordering,
                                              std::vector<byte_t>& elem_storage,
                                              std::vector<size_t>& elem_offsets) {
    reorder_swap(elem_reordering.begin(), elem_reordering.end(),
                 [&](uint_t first_elem_index, uint_t second_elem_index) {
        swap_elem_bytes_(first_elem_index, second_elem_index, elem_storage, elem_offsets);
    });
}   // reorder_elems_bytes_
/** @internal
 ** Change order of bytes of all nodes. */
SKUNK_INLINE void cMesh::reorder_nodes_bytes_(const std::vector<uint_t>& node_reordering) {
    reorder_elems_bytes_(node_reordering, m_node_storage, m_node_offsets);
}   // cMesh::reorder_nodes_bytes_
/** @internal
 ** Change order of bytes of all edges. */
SKUNK_INLINE void cMesh::reorder_edges_bytes_(const std::vector<uint_t>& edge_reordering) {
    reorder_elems_bytes_(edge_reordering, m_edge_storage, m_edge_offsets);
}   // cMesh::reorder_edges_bytes_
/** @internal
 ** Change order of bytes of all faces. */
SKUNK_INLINE void cMesh::reorder_faces_bytes_(const std::vector<uint_t>& face_reordering) {
    reorder_elems_bytes_(face_reordering, m_face_storage, m_face_offsets);
}   // cMesh::reorder_faces_bytes_
/** @internal
 ** Change order of bytes of all cells. */
SKUNK_INLINE void cMesh::reorder_cells_bytes_(const std::vector<uint_t>& cell_reordering) {
    reorder_elems_bytes_(cell_reordering, m_cell_storage, m_cell_offsets);
}   // cMesh::reorder_cells_bytes_

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

void cMesh::reorder_faces() {
    /** @todo Refactor me! */
    {
        std::vector<uint_t> node_reordering(num_nodes());
        std::iota(node_reordering.begin(), node_reordering.end(), 0);
        std::stable_sort(node_reordering.begin(), node_reordering.end(),
                         [&](uint_t node_ind_1, uint_t node_ind_2) {
            return get_mark(eNodeTag, node_ind_1) < get_mark(eNodeTag, node_ind_2);
        });
        reorder_nodes(node_reordering);
        for (uint_t node_index = 0; node_index < num_nodes(); ++node_index) {
            m_marked_node_ranges.resize(get_mark(eNodeTag, node_index) + 2);
            m_marked_node_ranges[get_mark(eNodeTag, node_index) + 1] += 1;
        }
        std::partial_sum(
            m_marked_node_ranges.begin(), m_marked_node_ranges.end(), m_marked_node_ranges.begin());
    }
    {
        std::vector<uint_t> edge_reordering(num_edges());
        std::iota(edge_reordering.begin(), edge_reordering.end(), 0);
        std::stable_sort(edge_reordering.begin(), edge_reordering.end(),
                         [&](uint_t edge_ind_1, uint_t edge_ind_2) {
            return get_mark(eEdgeTag, edge_ind_1) < get_mark(eEdgeTag, edge_ind_2);
        });
        reorder_edges(edge_reordering);
        for (uint_t edge_index = 0; edge_index < num_edges(); ++edge_index) {
            m_marked_edge_ranges.resize(get_mark(eEdgeTag, edge_index) + 2);
            m_marked_edge_ranges[get_mark(eEdgeTag, edge_index) + 1] += 1;
        }
        std::partial_sum(
            m_marked_edge_ranges.begin(), m_marked_edge_ranges.end(), m_marked_edge_ranges.begin());
    }
    {
        std::vector<uint_t> face_reordering(num_faces());
        std::iota(face_reordering.begin(), face_reordering.end(), 0);
        std::stable_sort(face_reordering.begin(), face_reordering.end(),
                         [&](uint_t face_ind_1, uint_t face_ind_2) {
            return get_mark(eFaceTag, face_ind_1) < get_mark(eFaceTag, face_ind_2);
        });
        reorder_faces(face_reordering);
        for (uint_t face_index = 0; face_index < num_faces(); ++face_index) {
            m_marked_face_ranges.resize(get_mark(eFaceTag, face_index) + 2);
            m_marked_face_ranges[get_mark(eFaceTag, face_index) + 1] += 1;
        }
        std::partial_sum(
            m_marked_face_ranges.begin(), m_marked_face_ranges.end(), m_marked_face_ranges.begin());
    }
    {
        std::vector<uint_t> cell_reordering(num_cells());
        std::iota(cell_reordering.begin(), cell_reordering.end(), 0);
        std::stable_sort(cell_reordering.begin(), cell_reordering.end(),
                         [&](uint_t cell_ind_1, uint_t cell_ind_2) {
            return get_mark(eCellTag, cell_ind_1) < get_mark(eCellTag, cell_ind_2);
        });
        reorder_cells(cell_reordering);
        for (uint_t cell_index = 0; cell_index < num_cells(); ++cell_index) {
            m_marked_cell_ranges.resize(get_mark(eCellTag, cell_index) + 2);
            m_marked_cell_ranges[get_mark(eCellTag, cell_index) + 1] += 1;
        }
        std::partial_sum(
            m_marked_cell_ranges.begin(), m_marked_cell_ranges.end(), m_marked_cell_ranges.begin());
    }
}   // cMesh::reorder_faces

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Generate edges using the face to node connectivity.
 * @warning This function may be slow and memory-consuming.
 */
void cMesh::generate_edges() {
    /* Face edge node index tables for various face types. */
    static const std::map<eShape, std::vector<std::vector<uint_t>>> face_edge_nodes = {
        /* 1D faces. */
        {eShape::node,         {{0} } },
        /* 2D faces. */
        {eShape::segment_2,    {{0},    {1} } },
        /* 3D faces. */
        {eShape::triangle_3,   {{1, 2}, {2, 0}, {0, 1} } },
        {eShape::quadrangle_4, {{0, 1}, {1, 2}, {2, 3}, {3, 0} } },
    };

    /* Build a map of the existing edges.
     * ( An edge is uniquely defined by a set of it's nodes. ) */
    std::map<std::set<uint_t>, uint_t> edge_cache;
    for (uint_t edge_index = 0; edge_index < num_edges(); ++edge_index) {
        const Edge& edge = get_edge(edge_index);
        edge_cache.emplace(std::set<uint_t>(edge.begin_node(), edge.end_node()), edge_index);
    }

    /* Generate missing edges. */
    for (uint_t face_index = 0; face_index < num_faces(); ++face_index) {
        Face& face = get_face(face_index);
        for (uint_t edge_local = 0; edge_local < face._num_edges(); ++edge_local) {
            uint_t& edge_index = face.begin_edge()[edge_local];
            if (edge_index != npos) {
                continue;
            }
            /* Collect edge nodes using the table. */
            std::vector<uint_t> edge_nodes;
            for (uint_t node_local : face_edge_nodes.at(face.get_type()).at(edge_local)) {
                FEATHERS_ASSERT(node_local < face._num_nodes());
                const uint_t node_index = face.begin_node()[node_local];
                edge_nodes.push_back(node_index);
            }
            /* Create the face or add current cell to the adjacency list. */
            std::set<uint_t> edge_nodes_set(edge_nodes.begin(), edge_nodes.end());
            const auto edge_iter = edge_cache.find(edge_nodes_set);
            if (edge_iter != edge_cache.cend()) {
                /* Locate the existing edge. */
                edge_index = edge_iter->second;
            } else {
                /* Create a brand-new edge. */
                edge_index = insert_edge(edge_nodes);
                edge_cache.emplace(edge_nodes_set, edge_index);
                set_mark(eEdgeTag, edge_index, get_mark(eFaceTag, face_index));
            }
        }
    }

    /* Check faces:
     * each face should be connected to all edges. */
    for (uint_t face_index = 0; face_index < num_faces(); ++face_index) {
        const Face& face = get_face(face_index);
        FEATHERS_ASSERT(std::all_of(face.begin_edge(), face.end_edge(), is_not_npos));
    }
}   // cMesh::generate_edges

/* Cell face node index tables for various cell types. */
static const std::map<
    eShape, std::vector<std::vector<uint_t>>> cell_face_nodes = {
    /* 1D cells. */
    {eShape::segment_2,
        {{0}, {1} } },
    /* 2D cells. */
    {eShape::triangle_3,
        {{1, 2},       {2, 0},       {0, 1} } },
    {eShape::quadrangle_4,
        {{0, 1},       {1, 2},       {2, 3},       {3, 0} } },
    /* 3D cells.
     * TODO: Add other cell types! */
    {eShape::tetrahedron_4,
        {{0, 2, 1},    {0, 1, 3},    {1, 2, 3},    {2, 0, 3} } },
    {eShape::hexahedron_8,
        {{0, 4, 7, 3}, {1, 2, 6, 5}, {3, 7, 6, 2}, {0, 1, 5, 4}, {0, 3, 2, 1}, {4, 5, 6, 7} } },
};

/**
 * Generate faces using the cell to node connectivity.
 * @warning This function may be slow and memory-consuming.
 */
void cMesh::generate_faces() {
    /* Build a map of the existing faces.
     * ( A face is uniquely defined by a set of it's nodes. ) */
    std::map<std::set<uint_t>, uint_t> face_cache;
    for (uint_t face_index = 0; face_index < num_faces(); ++face_index) {
        const Face& face = get_face(face_index);
        face_cache.emplace(std::set<uint_t>(face.begin_node(), face.end_node()), face_index);
    }

    /* Generate missing faces. */
    for (uint_t cell_index = 0; cell_index < num_cells(); ++cell_index) {
        Cell& cell = get_cell(cell_index);
        for (uint_t face_local = 0; face_local < cell._num_faces(); ++face_local) {
            uint_t& face_index = cell.begin_face()[face_local];
            if (face_index != npos) {
                continue;
            }
            /* Collect face nodes using the table. */
            std::vector<uint_t> face_nodes;
            for (uint_t node_local : cell_face_nodes.at(cell.get_type()).at(face_local)) {
                FEATHERS_ASSERT(node_local < cell._num_nodes());
                const uint_t node_index = cell.begin_node()[node_local];
                face_nodes.push_back(node_index);
            }
            /* Create the face or add current cell to the adjacency list. */
            const std::set<uint_t> face_nodes_set(face_nodes.begin(), face_nodes.end());
            const auto face_iter = face_cache.find(face_nodes_set);
            if (face_iter != face_cache.cend()) {
                /* Locate the existing face.
                 * ( We should be careful with face orientation here. ) */
                face_index = face_iter->second;
                Face& face = get_face(face_index);
                if (std::equal(face.begin_node(), face.end_node(), face_nodes.cbegin())) {
                    FEATHERS_ASSERT(face.get_inner_cell() == npos);
                    face.get_inner_cell() = cell_index;
                } else {
                    FEATHERS_ASSERT(face.get_outer_cell() == npos);
                    face.get_outer_cell() = cell_index;
                }
            } else {
                /* Create a brand-new face.
                 * Assign the current cell as the inner one. */
                face_index = insert_face(face_nodes);
                face_cache.emplace(face_nodes_set, face_index);
                Face& face = get_face(face_index);
                face.get_inner_cell() = cell_index;
            }
        }
    }

    /* Check cells:
     * each cell should be connected to all faces. */
    for (uint_t cell_index = 0; cell_index < num_cells(); ++cell_index) {
        const Cell& cell = get_cell(cell_index);
        FEATHERS_ASSERT(std::all_of(cell.begin_face(), cell.end_face(), is_not_npos));
    }
    /* Check faces:
     * each internal face should be connected to two cells;
     * each boundary face should be connected to at least one cell. */
    for (uint_t face_index = 0; face_index < num_faces(); ++face_index) {
        const Face& face = get_face(face_index);
        if (get_mark(eFaceTag, face_index) == 0) {
            FEATHERS_ASSERT(std::all_of(face.begin_cell(), face.end_cell(), is_not_npos));
        } else {
            FEATHERS_ASSERT(std::any_of(face.begin_cell(), face.end_cell(), is_not_npos));
        }
    }
}   // cMesh::generate_faces

/**
 * Generate boundary cells to complete face connectivity.
 */
void cMesh::generate_boundary_cells() {
    /* A node and edge flip table for various face types. */
    static const std::map<
        eShape, std::pair<std::vector<uint_t>, std::vector<uint_t>>> face_nodes_edges_flip_table {
        /* 1D faces. */
        {eShape::node,         {{0},          {0} } },
        /* 2D faces. */
        {eShape::segment_2,    {{1, 0},       {1, 0} } },
        /* 3D faces. */
        {eShape::triangle_3,   {{0, 2, 1},    {0, 2, 1} } },
        {eShape::quadrangle_4, {{0, 3, 2, 1}, {0, 3, 2, 1} } },
    };

    /* Generate boundary cells and fix boundary faces orientation. */
    for (uint_t face_index = 0; face_index < num_faces(); ++face_index) {
        Face& face = get_face(face_index);
        if (get_mark(eFaceTag, face_index) == 0) {
            continue;
        }
        /* Boundary faces should be oriented outwards from the mesh. */
        if (face.get_inner_cell() == npos) {
            /* Flip normal and cell connectivity. */
            set_face_normal(face_index, -get_face_normal(face_index));
            std::swap(face.get_inner_cell(), face.get_outer_cell());
            /* Flip node and edge connectivity. */
            std::vector<uint_t> node_reordering, edge_reordering;
            std::tie(node_reordering, edge_reordering) = face_nodes_edges_flip_table.at(face.get_type());
            reorder(node_reordering.begin(), node_reordering.end(), face.begin_node());
            reorder(edge_reordering.begin(), edge_reordering.end(), face.begin_edge());
        }
        /* Generate the boundary cell: reflect a connected interior cell. */
        const Cell& cell = get_cell(face.get_inner_cell());
        uint_t boundary_cell_dim;
        std::vector<uint_t> boundary_cell_nodes;
#if 0
        boundary_cell_dim = m_dim - 1;
        boundary_cell_nodes.assign(face.begin_node(), face.end_node());
#endif
#if 1
        boundary_cell_dim = m_dim;
        std::for_each(cell.begin_node(), cell.end_node(), [&](uint_t node_index) {
            const auto node_iter = std::find(face.begin_node(), face.end_node(), node_index);
            if (node_iter == face.end_node()) {
                /* Reflect an interior cell node. */
                vec3_t node_pos = get_node_coords(node_index);
                const vec3_t delta = node_pos - get_face_center_coords(face_index);
                node_pos -= 2.0*glm::dot(delta, get_face_normal(face_index))*get_face_normal(face_index);
                const mesh_node1_t node;
                boundary_cell_nodes.push_back(
                    insert_node(node, node_pos, get_mark(eFaceTag, face_index)));
            } else {
                /* Insert a boundary face node. */
                boundary_cell_nodes.push_back(node_index);
            }
        });
#endif
        /* Insert the boundary cell. */
        const uint_t boundary_cell_index =
            insert_cell(boundary_cell_nodes, get_mark(eFaceTag, face_index), boundary_cell_dim);
        Cell& boundary_cell = get_cell(boundary_cell_index);
        while (boundary_cell._num_faces() != 1) {
            boundary_cell.erase_face(0);
        }
        boundary_cell.begin_face()[0] = face_index;
        face.get_outer_cell() = boundary_cell_index;
    }
}   // cMesh::generate_boundary_cells

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
