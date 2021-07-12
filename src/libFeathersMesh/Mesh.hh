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

#pragma once
#ifndef MESH_HH_
#define MESH_HH_

#include <SkunkBase.hh>
#include "Shape.hh"

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/** Mesh element connectivity structure. */
template<eShape type_t,
         uint_t num_nodes_t, uint_t num_edges_t, uint_t num_faces_t, uint_t num_cells_t>
class mesh_elem_struct_t {
private:
    /** @internal
     ** @warning Order of the fields should not be changed! */
    /** @{ */
    eShape m_type : 8;
    uint_t m_num_nodes : 6, m_num_edges : 6, m_num_faces : 6, m_num_cells : 6;
    uint_t m_elem_conn[std::max(1u, num_nodes_t + num_edges_t + num_faces_t + num_cells_t)];
    /** @} */

    /**************************************************************************/
    /**************************************************************************/

protected:

    /** Construct mesh element storage. */
    mesh_elem_struct_t()
        : m_type(type_t),
          m_num_nodes(num_nodes_t), m_num_edges(num_edges_t),
          m_num_faces(num_faces_t), m_num_cells(num_cells_t), m_elem_conn() {
        std::fill(std::begin(m_elem_conn), std::end(m_elem_conn), npos);
    }
    /** Construct mesh element storage with nodes connectivity. */
    /** @{ */
    mesh_elem_struct_t(std::initializer_list<uint_t> nodes)
        : m_type(type_t),
          m_num_nodes(num_nodes_t), m_num_edges(num_edges_t),
          m_num_faces(num_faces_t), m_num_cells(num_cells_t), m_elem_conn() {
        std::fill(std::begin(m_elem_conn), std::end(m_elem_conn), npos);
        std::copy(nodes.begin(), nodes.end(), begin_node());
    }
    template<typename node_iter_t>
    mesh_elem_struct_t(node_iter_t begin_nodes, node_iter_t end_nodes)
        : m_type(type_t),
          m_num_nodes(num_nodes_t), m_num_edges(num_edges_t),
          m_num_faces(num_faces_t), m_num_cells(num_cells_t), m_elem_conn() {
        std::fill(std::begin(m_elem_conn), std::end(m_elem_conn), npos);
        std::copy(begin_nodes, end_nodes, begin_node());
    }
    /** @} */

public://private:
    friend class cMesh;

    /**************************************************************************/
    /**************************************************************************/

    /** Type of the element. */
    eShape get_type() const {
        return m_type;
    }

    /**************************************************************************/
    /**************************************************************************/

    /** Number of nodes in the element. */
    uint_t _num_nodes() const {
        return m_num_nodes;
    }
    /** Number of edges in the element. */
    uint_t _num_edges() const {
        return m_num_edges;
    }
    /** Number of faces in the element. */
    uint_t _num_faces() const {
        return m_num_faces;
    }
    /** Number of cells in the element. */
    uint_t _num_cells() const {
        return m_num_cells;
    }

    /** Pointer to the beginning of the connected nodes. */
    /** @{ */
    uint_t* begin_node() {
        return m_elem_conn;
    }
    const uint_t* begin_node() const {
        return m_elem_conn;
    }
    /** @} */
    /** Pointer to the beginning of the connected edges. */
    /** @{ */
    uint_t* begin_edge() {
        return m_elem_conn + m_num_nodes;
    }
    const uint_t* begin_edge() const {
        return m_elem_conn + m_num_nodes;
    }
    /** @} */
    /** Pointer to the beginning of the connected faces. */
    /** @{ */
    uint_t* begin_face() {
        return m_elem_conn + m_num_nodes + m_num_edges;
    }
    const uint_t* begin_face() const {
        return m_elem_conn + m_num_nodes + m_num_edges;
    }
    /** @} */
    /** Pointer to the beginning of the connected cells. */
    /** @{ */
    uint_t* begin_cell() {
        return m_elem_conn + m_num_nodes + m_num_edges + m_num_faces;
    }
    const uint_t* begin_cell() const {
        return m_elem_conn + m_num_nodes + m_num_edges + m_num_faces;
    }
    /** @} */

    /** Pointer to the end of the connected nodes. */
    /** @{ */
    uint_t* end_node() {
        return m_elem_conn + m_num_nodes;
    }
    const uint_t* end_node() const {
        return m_elem_conn + m_num_nodes;
    }
    /** @} */
    /** Pointer to the end of the connected edges. */
    /** @{ */
    uint_t* end_edge() {
        return m_elem_conn + m_num_nodes + m_num_edges;
    }
    const uint_t* end_edge() const {
        return m_elem_conn + m_num_nodes + m_num_edges;
    }
    /** @} */
    /** Pointer to the end of the connected faces. */
    /** @{ */
    uint_t* end_face() {
        return m_elem_conn + m_num_nodes + m_num_edges + m_num_faces;
    }
    const uint_t* end_face() const {
        return m_elem_conn + m_num_nodes + m_num_edges + m_num_faces;
    }
    /** @} */
    /** Pointer to the end of the connected cells. */
    /** @{ */
    uint_t* end_cell() {
        return m_elem_conn + m_num_nodes + m_num_edges + m_num_faces + m_num_cells;
    }
    const uint_t* end_cell() const {
        return m_elem_conn + m_num_nodes + m_num_edges + m_num_faces + m_num_cells;
    }
    /** @} */

    /**************************************************************************/
    /**************************************************************************/

    /** Erase the face from element.
     ** @warning Element type would be set to unknown. */
    void erase_face(uint_t face_loc) {
        FEATHERS_ASSERT(face_loc < _num_faces());
        std::rotate(begin_face() + face_loc, begin_face() + face_loc + 1, end_cell());
        m_type = eShape::null;
        m_num_faces -= 1;
    }
};  // class mesh_elem_struct_t

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/** Node element. */
template<eShape type_t>
class mesh_node_struct_t final : public mesh_elem_struct_t<type_t, 0, 0, 0, 0> {
};  // class mesh_node_struct_t

/** Generic node element. */
using Node = mesh_node_struct_t<eShape::null>;

/** Node element. */
using mesh_node1_t = mesh_node_struct_t<eShape::node>;

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/** Edge element. */
template<eShape type_t, uint_t num_nodes_t>
class mesh_edge_struct_t : public mesh_elem_struct_t<type_t, num_nodes_t, 0, 0, 0> {
public:
    /** Construct an edge element with nodes connectivity. */
    /** @{ */
    mesh_edge_struct_t(std::initializer_list<uint_t> nodes)
        : mesh_elem_struct_t<type_t, num_nodes_t, 0, 0, 0>(nodes) {
    }
    template<typename node_iter_t>
    mesh_edge_struct_t(node_iter_t begin_nodes, node_iter_t end_nodes)
        : mesh_elem_struct_t<type_t, num_nodes_t, 0, 0, 0>(begin_nodes, end_nodes) {
    }
    /** @} */
};  // class mesh_edge_struct_t

/** Generic edge element. */
using Edge = mesh_edge_struct_t<eShape::null, 0>;

/** Dummy "node" edge element.
 ** Acts as a placeholder fo 1D/2D meshes. */
using mesh_edge_node1_t = mesh_edge_struct_t<eShape::node, 1>;

/** Linear segment edge element. */
using mesh_edge_segment2_t = mesh_edge_struct_t<eShape::segment_2, 2>;

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/** Face element. */
template<eShape type_t, uint_t num_nodes_t, uint_t num_edges_t>
class mesh_face_struct_t : public mesh_elem_struct_t<type_t, num_nodes_t, num_edges_t, 0, 2> {
public:
    /** Construct a face element with nodes connectivity. */
    /** @{ */
    mesh_face_struct_t(std::initializer_list<uint_t> nodes)
        : mesh_elem_struct_t<type_t, num_nodes_t, num_edges_t, 0, 2>(nodes) {
    }
    template<typename node_iter_t>
    mesh_face_struct_t(node_iter_t begin_nodes, node_iter_t end_nodes)
        : mesh_elem_struct_t<type_t, num_nodes_t, num_edges_t, 0, 2>(begin_nodes, end_nodes) {
    }
    /** @} */

    /** Index of the inner cell. */
    /** @{ */
    uint_t& get_inner_cell() {
        return this->begin_cell()[0];
    }
    const uint_t& get_inner_cell() const {
        return this->begin_cell()[0];
    }
    /** @} */
    /** Index of the outer cell. */
    /** @{ */
    uint_t& get_outer_cell() {
        return this->begin_cell()[1];
    }
    const uint_t& get_outer_cell() const {
        return this->begin_cell()[1];
    }
    /** @} */
};  // class mesh_face_struct_t

/** Generic face element. */
using Face = mesh_face_struct_t<eShape::null, 0, 0>;

/** Dummy "node" face element.
 ** Acts as a placeholder for 1D meshes. */
using mesh_face_node1_t = mesh_face_struct_t<eShape::node, 1, 1>;

/** Linear segment face element for 2D meshes. */
using mesh_face_segment2_t = mesh_face_struct_t<eShape::segment_2, 2, 2>;

/** Linear triangle face element. */
using mesh_face_triangle3_t = mesh_face_struct_t<eShape::triangle_3, 3, 3>;

/** Linear quadrangle face element. */
using mesh_face_quadrangle4_t = mesh_face_struct_t<eShape::quadrangle_4, 4, 4>;

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/** Cell element. */
template<eShape type_t, uint_t num_nodes_t, uint_t num_faces_t>
class mesh_cell_struct_t : public mesh_elem_struct_t<type_t, num_nodes_t, 0, num_faces_t, 0> {
public:
    /** Construct a cell element with nodes connectivity. */
    /** @{ */
    mesh_cell_struct_t(std::initializer_list<uint_t> nodes)
        : mesh_elem_struct_t<type_t, num_nodes_t, 0, num_faces_t, 0>(nodes) {
    }
    template<typename node_iter_t>
    mesh_cell_struct_t(node_iter_t begin_nodes, node_iter_t end_nodes)
        : mesh_elem_struct_t<type_t, num_nodes_t, 0, num_faces_t, 0>(begin_nodes, end_nodes) {
    }
    /** @} */
};  // class mesh_cell_struct_t

/** Generic cell element. */
using Cell = mesh_cell_struct_t<eShape::null, 0, 0>;

/** Linear segment cell element for 1D meshes. */
using mesh_cell_segment2_t = mesh_cell_struct_t<eShape::segment_2, 2, 2>;

/** Linear triangle cell element for 2D meshes. */
using mesh_cell_triangle3_t = mesh_cell_struct_t<eShape::triangle_3, 3, 3>;

/** Linear quadrangle cell element for 2D meshes. */
using mesh_cell_quadrangle4_t = mesh_cell_struct_t<eShape::quadrangle_4, 4, 4>;

/** Linear tetrahedron cell element. */
using mesh_cell_tetrahedron4_t = mesh_cell_struct_t<eShape::tetrahedron_4, 4, 4>;

/** Linear pyramid cell element. */
using mesh_cell_pyramid5_t = mesh_cell_struct_t<eShape::pyramid_5, 5, 5>;

/** Linear pentahedron cell element. */
using mesh_cell_pentahedron6_t = mesh_cell_struct_t<eShape::pentahedron_6, 6, 5>;

/** Linear hexahedron cell element. */
using mesh_cell_hexahedron8_t = mesh_cell_struct_t<eShape::hexahedron_8, 8, 6>;

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/** Node element tag. */
enum tNodeTag { eNodeTag = 1 };
/** Edge element tag. */
enum tEdgeTag { eEdgeTag = 2 };
/** Face element tag. */
enum tFaceTag { eFaceTag = 3 };
/** Cell element tag. */
enum tCellTag { eCellTag = 4 };

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/** Hybrid unstructured finite element mesh. */
class cMesh : public tObject<cMesh> {
private:
    std::vector<eShape> m_edge_types, m_face_types, m_cell_types;

private:
    uint_t m_dim;

    std::vector<byte_t> m_node_storage;
    std::vector<byte_t> m_edge_storage;
    std::vector<byte_t> m_face_storage;
    std::vector<byte_t> m_cell_storage;
    std::vector<size_t> m_node_offsets;
    std::vector<size_t> m_edge_offsets;
    std::vector<size_t> m_face_offsets;
    std::vector<size_t> m_cell_offsets;

    std::vector<uint_t> m_node_marks;
    std::vector<uint_t> m_edge_marks;
    std::vector<uint_t> m_face_marks;
    std::vector<uint_t> m_cell_marks;
    std::vector<uint_t> m_marked_node_ranges{0};
    std::vector<uint_t> m_marked_edge_ranges{0};
    std::vector<uint_t> m_marked_face_ranges{0};
    std::vector<uint_t> m_marked_cell_ranges{0};

    std::vector<vec3_t> m_node_coords;
    std::vector<real_t> m_edge_lengths;
    std::vector<vec3_t> m_edge_directions;
    std::vector<real_t> m_face_areas;
    std::vector<vec3_t> m_face_normals;
    std::vector<vec3_t> m_face_center_coords;
    std::vector<real_t> m_cell_volumes;
    std::vector<vec3_t> m_cell_center_coords;
    mutable real_t m_min_edge_length = 0.0, m_max_edge_length = 0.0;
    mutable real_t m_min_face_area = 0.0, m_max_face_area = 0.0;
    mutable real_t m_min_cell_volume = 0.0, m_max_cell_volume = 0.0;

public:

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Initialize an empty mesh. */
    explicit cMesh(uint_t dim = 0) : m_dim(dim) {
        FEATHERS_ASSERT(0 <= m_dim && m_dim <= 3);
    }

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    bool read_triangle(const char* path);
    bool read_tetgen(const char* path);

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Number of the mesh dimensions. */
    uint_t num_dim() const {
        return m_dim;
    }

    /** Total number of nodes in the mesh. */
    uint_t num_nodes() const {
        return m_node_offsets.size();
    }
    /** Total number of edges in the mesh. */
    uint_t num_edges() const {
        return m_edge_offsets.size();
    }
    /** Total number of faces in the mesh. */
    uint_t num_faces() const {
        return m_face_offsets.size();
    }
    /** Total number of cells in the mesh. */
    uint_t num_cells() const {
        return m_cell_offsets.size();
    }

    /** Get node element at global index. */
    /** @{ */
    Node& get_node(uint_t node_index) {
        FEATHERS_ASSERT(node_index < num_nodes());
        return reinterpret_cast<Node&>(m_node_storage[m_node_offsets[node_index]]);
    }
    const Node& get_node(uint_t node_index) const {
        FEATHERS_ASSERT(node_index < num_nodes());
        return reinterpret_cast<const Node&>(m_node_storage[m_node_offsets[node_index]]);
    }
    /** @} */
    /** Get edge element at global index. */
    /** @{ */
    Edge& get_edge(uint_t edge_index) {
        FEATHERS_ASSERT(edge_index < num_edges());
        return reinterpret_cast<Edge&>(m_edge_storage[m_edge_offsets[edge_index]]);
    }
    const Edge& get_edge(uint_t edge_index) const {
        FEATHERS_ASSERT(edge_index < num_edges());
        return reinterpret_cast<const Edge&>(m_edge_storage[m_edge_offsets[edge_index]]);
    }
    /** @} */
    /** Get face element at global index. */
    /** @{ */
    Face& get_face(uint_t face_index) {
        FEATHERS_ASSERT(face_index < num_faces());
        return reinterpret_cast<Face&>(m_face_storage[m_face_offsets[face_index]]);
    }
    const Face& get_face(uint_t face_index) const {
        FEATHERS_ASSERT(face_index < num_faces());
        return reinterpret_cast<const Face&>(m_face_storage[m_face_offsets[face_index]]);
    }
    /** @} */
    /** Get cell element at global index. */
    /** @{ */
    Cell& get_cell(uint_t cell_index) {
        FEATHERS_ASSERT(cell_index < num_cells());
        return reinterpret_cast<Cell&>(m_cell_storage[m_cell_offsets[cell_index]]);
    }
    const Cell& get_cell(uint_t cell_index) const {
        FEATHERS_ASSERT(cell_index < num_cells());
        return reinterpret_cast<const Cell&>(m_cell_storage[m_cell_offsets[cell_index]]);
    }
    /** @} */

    /** Pointer to the beginning of the element adjacent nodes. */
    /** @{ */
#if FEATHERS_DOXYGEN
    FEATHERS_CONST_OVERLOAD_T(template<typename tTag>,
        uint_t*, begin_adjacent_node, (tTag, uint_t index), ;)
#else
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_node, (tNodeTag, uint_t node_index), {
        FEATHERS_ASSERT(node_index < num_nodes());
        return get_node(node_index).begin_node();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_node, (tEdgeTag, uint_t edge_index), {
        FEATHERS_ASSERT(edge_index < num_edges());
        return get_edge(edge_index).begin_node();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_node, (tFaceTag, uint_t face_index), {
        FEATHERS_ASSERT(face_index < num_faces());
        return get_face(face_index).begin_node();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_node, (tCellTag, uint_t cell_index), {
        FEATHERS_ASSERT(cell_index < num_cells());
        return get_cell(cell_index).begin_node();
    })
#endif
    /** @} */
    /** Pointer to the beginning of the element adjacent edges. */
    /** @{ */
#if FEATHERS_DOXYGEN
    FEATHERS_CONST_OVERLOAD_T(template<typename tTag>,
        uint_t*, begin_adjacent_edge, (tTag, uint_t index), ;)
#else
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_edge, (tNodeTag, uint_t node_index), {
        FEATHERS_ASSERT(node_index < num_nodes());
        return get_node(node_index).begin_edge();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_edge, (tEdgeTag, uint_t edge_index), {
        FEATHERS_ASSERT(edge_index < num_edges());
        return get_edge(edge_index).begin_edge();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_edge, (tFaceTag, uint_t face_index), {
        FEATHERS_ASSERT(face_index < num_faces());
        return get_face(face_index).begin_edge();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_edge, (tCellTag, uint_t cell_index), {
        FEATHERS_ASSERT(cell_index < num_cells());
        return get_cell(cell_index).begin_edge();
    })
#endif
    /** @} */
    /** Pointer to the beginning of the element adjacent faces. */
    /** @{ */
#if FEATHERS_DOXYGEN
    FEATHERS_CONST_OVERLOAD_T(template<typename tTag>,
        uint_t*, begin_adjacent_face, (tTag, uint_t index), ;)
#else
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_face, (tNodeTag, uint_t node_index), {
        FEATHERS_ASSERT(node_index < num_nodes());
        return get_node(node_index).begin_face();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_face, (tEdgeTag, uint_t edge_index), {
        FEATHERS_ASSERT(edge_index < num_edges());
        return get_edge(edge_index).begin_face();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_face, (tFaceTag, uint_t face_index), {
        FEATHERS_ASSERT(face_index < num_faces());
        return get_face(face_index).begin_face();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_face, (tCellTag, uint_t cell_index), {
        FEATHERS_ASSERT(cell_index < num_cells());
        return get_cell(cell_index).begin_face();
    })
#endif
    /** @} */
    /** Pointer to the beginning of the element adjacent cells. */
    /** @{ */
#if FEATHERS_DOXYGEN
    FEATHERS_CONST_OVERLOAD_T(template<typename tTag>,
        uint_t*, begin_adjacent_cell, (tTag, uint_t index), ;)
#else
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_cell, (tNodeTag, uint_t node_index), {
        FEATHERS_ASSERT(node_index < num_nodes());
        return get_node(node_index).begin_cell();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_cell, (tEdgeTag, uint_t edge_index), {
        FEATHERS_ASSERT(edge_index < num_edges());
        return get_edge(edge_index).begin_cell();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_cell, (tFaceTag, uint_t face_index), {
        FEATHERS_ASSERT(face_index < num_faces());
        return get_face(face_index).begin_cell();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_cell, (tCellTag, uint_t cell_index), {
        FEATHERS_ASSERT(cell_index < num_cells());
        return get_cell(cell_index).begin_cell();
    })
#endif
    /** @} */

    /** Pointer to the end of the element adjacent nodes. */
    /** @{ */
#if FEATHERS_DOXYGEN
    FEATHERS_CONST_OVERLOAD_T(template<typename tTag>,
        uint_t*, end_adjacent_node, (tTag, uint_t index), ;)
#else
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_node, (tNodeTag, uint_t node_index), {
        FEATHERS_ASSERT(node_index < num_nodes());
        return get_node(node_index).end_node();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_node, (tEdgeTag, uint_t edge_index), {
        FEATHERS_ASSERT(edge_index < num_edges());
        return get_edge(edge_index).end_node();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_node, (tFaceTag, uint_t face_index), {
        FEATHERS_ASSERT(face_index < num_faces());
        return get_face(face_index).end_node();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_node, (tCellTag, uint_t cell_index), {
        FEATHERS_ASSERT(cell_index < num_cells());
        return get_cell(cell_index).end_node();
    })
#endif
    /** @} */
    /** Pointer to the end of the element adjacent edges. */
    /** @{ */
#if FEATHERS_DOXYGEN
    FEATHERS_CONST_OVERLOAD_T(template<typename tTag>,
        uint_t*, end_adjacent_edge, (tTag, uint_t index), ;)
#else
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_edge, (tNodeTag, uint_t node_index), {
        FEATHERS_ASSERT(node_index < num_nodes());
        return get_node(node_index).end_edge();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_edge, (tEdgeTag, uint_t edge_index), {
        FEATHERS_ASSERT(edge_index < num_edges());
        return get_edge(edge_index).end_edge();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_edge, (tFaceTag, uint_t face_index), {
        FEATHERS_ASSERT(face_index < num_faces());
        return get_face(face_index).end_edge();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_edge, (tCellTag, uint_t cell_index), {
        FEATHERS_ASSERT(cell_index < num_cells());
        return get_cell(cell_index).end_edge();
    })
#endif
    /** @} */
    /** Pointer to the end of the element adjacent faces. */
    /** @{ */
#if FEATHERS_DOXYGEN
    FEATHERS_CONST_OVERLOAD_T(template<typename tTag>,
        uint_t*, end_adjacent_face, (tTag, uint_t index), ;)
#else
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_face, (tNodeTag, uint_t node_index), {
        FEATHERS_ASSERT(node_index < num_nodes());
        return get_node(node_index).end_face();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_face, (tEdgeTag, uint_t edge_index), {
        FEATHERS_ASSERT(edge_index < num_edges());
        return get_edge(edge_index).end_face();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_face, (tFaceTag, uint_t face_index), {
        FEATHERS_ASSERT(face_index < num_faces());
        return get_face(face_index).end_face();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_face, (tCellTag, uint_t cell_index), {
        FEATHERS_ASSERT(cell_index < num_cells());
        return get_cell(cell_index).end_face();
    })
#endif
    /** @} */
    /** Pointer to the end of the element adjacent cells. */
    /** @{ */
#if FEATHERS_DOXYGEN
    FEATHERS_CONST_OVERLOAD_T(template<typename tTag>,
        uint_t*, end_adjacent_cell, (tTag, uint_t index), ;)
#else
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_cell, (tNodeTag, uint_t node_index), {
        FEATHERS_ASSERT(node_index < num_nodes());
        return get_node(node_index).end_cell();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_cell, (tEdgeTag, uint_t edge_index), {
        FEATHERS_ASSERT(edge_index < num_edges());
        return get_edge(edge_index).end_cell();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_cell, (tFaceTag, uint_t face_index), {
        FEATHERS_ASSERT(face_index < num_faces());
        return get_face(face_index).end_cell();
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_cell, (tCellTag, uint_t cell_index), {
        FEATHERS_ASSERT(cell_index < num_cells());
        return get_cell(cell_index).end_cell();
    })
#endif
    /** @} */

    /** Number of nodes in the element. */
    template<typename tTag>
    uint_t num_adjacent_nodes(tTag tag, uint_t index) const {
        return end_adjacent_node(tag, index) - begin_adjacent_node(tag, index);
    }
    /** Number of edges in the element. */
    template<typename tTag>
    uint_t num_adjacent_edges(tTag tag, uint_t index) const {
        return end_adjacent_edge(tag, index) - begin_adjacent_edge(tag, index);
    }
    /** Number of faces in the element. */
    template<typename tTag>
    uint_t num_adjacent_faces(tTag tag, uint_t index) const {
        return end_adjacent_face(tag, index) - begin_adjacent_face(tag, index);
    }
    /** Number of cells in the element. */
    template<typename tTag>
    uint_t num_adjacent_cells(tTag tag, uint_t index) const {
        return end_adjacent_cell(tag, index) - begin_adjacent_cell(tag, index);
    }
    
    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Number of node marks. */
    uint_t num_node_marks() const {
        return m_marked_node_ranges.size() - 1;
    }
    /** Number of edge marks. */
    uint_t num_edge_marks() const {
        return m_marked_edge_ranges.size() - 1;
    }
    /** Number of face marks. */
    uint_t num_face_marks() const {
        return m_marked_face_ranges.size() - 1;
    }
    /** Number of cell marks. */
    uint_t num_cell_marks() const {
        return m_marked_cell_ranges.size() - 1;
    }

    /** Get element mark. */
    /** @{ */
#if FEATHERS_DOXYGEN
    template<typename tTag>
    uint_t get_mark(tTag, uint_t index) const;
#else
    uint_t get_mark(tNodeTag, uint_t node_index) const {
        FEATHERS_ASSERT(node_index < num_nodes());
        return m_node_marks[node_index];
    }
    uint_t get_mark(tEdgeTag, uint_t edge_index) const {
        FEATHERS_ASSERT(edge_index < num_edges());
        return m_edge_marks[edge_index];
    }
    uint_t get_mark(tFaceTag, uint_t face_index) const {
        FEATHERS_ASSERT(face_index < num_faces());
        return m_face_marks[face_index];
    }
    uint_t get_mark(tCellTag, uint_t cell_index) const {
        FEATHERS_ASSERT(cell_index < num_cells());
        return m_cell_marks[cell_index];
    }
#endif
    /** @} */
    /** Set element mark. */
    /** @{ */
#if FEATHERS_DOXYGEN
    template<typename tTag>
    void set_mark(tTag, uint_t index, uint_t mark);
#else
    void set_mark(tNodeTag, uint_t node_index, uint_t mark) {
        FEATHERS_ASSERT(node_index < num_nodes());
        m_node_marks[node_index] = mark;
    }
    void set_mark(tEdgeTag, uint_t edge_index, uint_t mark) {
        FEATHERS_ASSERT(edge_index < num_edges());
        m_edge_marks[edge_index] = mark;
    }
    void set_mark(tFaceTag, uint_t face_index, uint_t mark) {
        FEATHERS_ASSERT(face_index < num_faces());
        m_face_marks[face_index] = mark;
    }
    void set_mark(tCellTag, uint_t cell_index, uint_t mark) {
        FEATHERS_ASSERT(cell_index < num_cells());
        m_cell_marks[cell_index] = mark;
    }
#endif
    /** @} */

    /** Index of the first node with a given mark. */
    uint_t begin_node(uint_t mark) const {
        FEATHERS_ASSERT(mark < num_node_marks());
        return m_marked_node_ranges[mark];
    }
    /** Index of the first edge with a given mark. */
    uint_t begin_edge(uint_t mark) const {
        FEATHERS_ASSERT(mark < num_edge_marks());
        return m_marked_edge_ranges[mark];
    }
    /** Index of the first face with a given mark. */
    uint_t begin_face(uint_t mark) const {
        FEATHERS_ASSERT(mark < num_face_marks());
        return m_marked_face_ranges[mark];
    }
    /** Index of the first cell with a given mark. */
    uint_t begin_cell(uint_t mark) const {
        FEATHERS_ASSERT(mark < num_cell_marks());
        return m_marked_cell_ranges[mark];
    }

    /** Index of a node after the last node with a given mark. */
    uint_t end_node(uint_t mark) const {
        FEATHERS_ASSERT(mark < num_node_marks());
        return m_marked_node_ranges[mark + 1];
    }
    /** Index of an edge after the last edge with a given mark. */
    uint_t end_edge(uint_t mark) const {
        FEATHERS_ASSERT(mark < num_edge_marks());
        return m_marked_edge_ranges[mark + 1];
    }
    /** Index of a face after the last face with a given mark. */
    uint_t end_face(uint_t mark) const {
        FEATHERS_ASSERT(mark < num_face_marks());
        return m_marked_face_ranges[mark + 1];
    }
    /** Index a cell after of the last cell with a given mark. */
    uint_t end_cell(uint_t mark) const {
        FEATHERS_ASSERT(mark < num_cell_marks());
        return m_marked_cell_ranges[mark + 1];
    }

    /** Number of marked cells in the mesh. */
    uint_t num_marked_cells(uint_t mark) const {
        return end_cell(mark) - begin_cell(mark);
    }

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get node position. */
    const vec3_t& get_node_coords(uint_t node_index) const {
        FEATHERS_ASSERT(node_index < num_nodes());
        return m_node_coords[node_index];
    }
    /** Set node position. */
    void set_node_coords(uint_t node_index, const vec3_t& node_coords) {
        FEATHERS_ASSERT(node_index < num_nodes());
        m_node_coords[node_index] = node_coords;
    }

    /** Get minimal edge length. */
    real_t get_min_edge_length() const;
    /** Get maximal edge length. */
    real_t get_max_edge_length() const;

    /** Get edge length. */
    real_t get_edge_length(uint_t edge_index) const {
        FEATHERS_ASSERT(edge_index < num_edges());
        return m_edge_lengths[edge_index];
    }
    /** Set edge length. */
    void set_edge_length(uint_t edge_index, real_t edge_length) {
        FEATHERS_ASSERT(edge_index < num_edges());
        m_edge_lengths[edge_index] = edge_length;
    }
    /** Get edge direction. */
    const vec3_t& get_edge_direction(uint_t edge_index) const {
        FEATHERS_ASSERT(edge_index < num_edges());
        return m_edge_directions[edge_index];
    }
    /** Set edge direction. */
    void set_edge_direction(uint_t edge_index, const vec3_t& edge_direction) {
        FEATHERS_ASSERT(edge_index < num_edges());
        m_edge_directions[edge_index] = edge_direction;
    }

    /** Get minimal face area. */
    real_t get_min_face_area() const;
    /** Get maximal face area. */
    real_t get_max_face_area() const;

    /** Get face area/length. */
    real_t get_face_area(uint_t face_index) const {
        FEATHERS_ASSERT(face_index < num_faces());
        return m_face_areas[face_index];
    }
    /** Set face area/length. */
    void set_face_area(uint_t face_index, real_t face_area) {
        FEATHERS_ASSERT(face_index < num_faces());
        m_face_areas[face_index] = face_area;
    }
    /** Get face normal. */
    const vec3_t& get_face_normal(uint_t face_index) const {
        FEATHERS_ASSERT(face_index < num_faces());
        return m_face_normals[face_index];
    }
    /** Set face normal. */
    void set_face_normal(uint_t face_index, const vec3_t& face_normal) {
        FEATHERS_ASSERT(face_index < num_faces());
        m_face_normals[face_index] = face_normal;
    }
    /** Get face barycenter. */
    const vec3_t& get_face_center_coords(uint_t face_index) const {
        FEATHERS_ASSERT(face_index < num_faces());
        return m_face_center_coords[face_index];
    }
    /** Set face barycenter. */
    void set_face_center_coords(uint_t face_index, const vec3_t& face_center_coords) {
        FEATHERS_ASSERT(face_index < num_faces());
        m_face_center_coords[face_index] = face_center_coords;
    }

    /** Get minimal cell volume. */
    real_t get_min_cell_volume() const;
    /** Get maximal cell volume. */
    real_t get_max_cell_volume() const;

    /** Get cell volume/area/length. */
    real_t get_cell_volume(uint_t cell_index) const {
        FEATHERS_ASSERT(cell_index < num_cells());
        return m_cell_volumes[cell_index];
    }
    /** Set cell volume/area/length. */
    void set_cell_volume(uint_t cell_index, real_t cell_volume) {
        FEATHERS_ASSERT(cell_index < num_cells());
        m_cell_volumes[cell_index] = cell_volume;
    }
    /** Get cell barycenter. */
    const vec3_t& get_cell_center_coords(uint_t cell_index) const {
        FEATHERS_ASSERT(cell_index < num_cells());
        return m_cell_center_coords[cell_index];
    }
    /** Get cell barycenter. */
    void set_cell_center_coords(uint_t cell_index, const vec3_t& cell_center_coords) {
        FEATHERS_ASSERT(cell_index < num_cells());
        m_cell_center_coords[cell_index] = cell_center_coords;
     }

    /** Compute mesh orthogonality.
     ** Mesh orthogonality is defined as a ... */
    real_t get_orthogonality() const;

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

protected:

    /** Insert a new node into the mesh.
     ** @param node Edge nodes.
     ** @returns Index of the inserted node. */
    uint_t insert_node(const mesh_node1_t& node, vec3_t p, uint_t mark = 0);

    /** Insert a new edge into the mesh.
     ** @param edge Edge nodes.
     ** @returns Index of the inserted edge. */
    /** @{ */
    uint_t insert_edge(const mesh_edge_node1_t& edge, uint_t mark = 0);
    uint_t insert_edge(const mesh_edge_segment2_t& edge, uint_t mark = 0);
    /** @} */

    /** Insert a new edge into the mesh.
     ** Type of the edge is detected by the set of nodes and dimension.
     ** @param edge_nodes Set of the edge nodes.
     ** @param mark Boundary mark.
     ** @param dim Overridden dimension.
     ** @returns Index of the inserted edge. */
    uint_t insert_edge(const std::vector<uint_t>& edge_nodes, uint_t mark = 0, uint_t dim = 0);

    /** Insert a new face into the mesh.
     ** @param face Face nodes.
     ** @returns Index of the inserted face. */
    /** @{ */
    uint_t insert_face(const mesh_face_node1_t& face, uint_t mark = 0);
    uint_t insert_face(const mesh_face_segment2_t& face, uint_t mark = 0);
    uint_t insert_face(const mesh_face_triangle3_t& face, uint_t mark = 0);
    uint_t insert_face(const mesh_face_quadrangle4_t& face, uint_t mark = 0);
    /** @} */

    /** Insert a new face into the mesh.
     ** Type of the face is detected by the set of nodes and dimension.
     ** @param face_nodes Set of the face nodes.
     ** @param mark Boundary mark.
     ** @param dim Overridden dimension.
     ** @returns Index of the inserted face. */
    uint_t insert_face(const std::vector<uint_t>& face_nodes, uint_t mark = 0, uint_t dim = 0);

    /** Insert a new cell into the mesh.
     ** @param cell Cell nodes.
     ** @returns Index of the inserted face. */
    /** @{ */
    uint_t insert_cell(const mesh_cell_segment2_t& cell, uint_t mark = 0);
    uint_t insert_cell(const mesh_cell_triangle3_t& cell, uint_t mark = 0);
    uint_t insert_cell(const mesh_cell_quadrangle4_t& cell, uint_t mark = 0);
    uint_t insert_cell(const mesh_cell_tetrahedron4_t& cell, uint_t mark = 0);
    uint_t insert_cell(const mesh_cell_pyramid5_t& cell, uint_t mark = 0);
    uint_t insert_cell(const mesh_cell_pentahedron6_t& cell, uint_t mark = 0);
    uint_t insert_cell(const mesh_cell_hexahedron8_t& cell, uint_t mark = 0);
    /** @} */

    /** Insert a new cell into the mesh.
     ** Type of the cell is detected by the set of nodes and dimension.
     ** @param cell_nodes Set of the cell nodes.
     ** @param mark Boundary mark.
     ** @param dim Overridden dimension.
     ** @returns Index of the inserted cell. */
    uint_t insert_cell(const std::vector<uint_t>& cell_nodes, uint_t mark = 0, uint_t dim = 0);

private:

    /** @internal
     ** Default alignment for elements storage. */
    static constexpr size_t s_default_alignment = sizeof(size_t);

    /** @internal
     ** Allocate and insert a new node into the mesh. */
    template<typename mesh_node_struct_t>
    uint_t allocate_node_(const mesh_node_struct_t& node, size_t alignment = s_default_alignment);
    /** @internal
     ** Allocate and insert a new edge into the mesh. */
    template<typename mesh_edge_struct_t>
    uint_t allocate_edge_(const mesh_edge_struct_t& edge, size_t alignment = s_default_alignment);
    /** @internal
     ** Allocate and insert a new face into the mesh. */
    template<typename mesh_face_struct_t>
    uint_t allocate_face_(const mesh_face_struct_t& face, size_t alignment = s_default_alignment);
    /** @internal
     ** Allocate and insert a new cell into the mesh. */
    template<typename mesh_cell_struct_t>
    uint_t allocate_cell_(const mesh_cell_struct_t& cell, size_t alignment = s_default_alignment);

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

protected:

    /** Change order of all nodes. */
    void reorder_nodes(const std::vector<uint_t>& node_reordering);
    /** Change order of all edges. */
    void reorder_edges(const std::vector<uint_t>& edge_reordering);
    /** Change order of all faces. */
    void reorder_faces(const std::vector<uint_t>& face_reordering);
    /** Change order of all cells. */
    void reorder_cells(const std::vector<uint_t>& cell_reordering);

private:

    /** @internal
     ** Change order of bytes of all nodes. */
    void reorder_nodes_bytes_(const std::vector<uint_t>& node_reordering);
    /** @internal
     ** Change order of bytes of all edges. */
    void reorder_edges_bytes_(const std::vector<uint_t>& edge_reordering);
    /** @internal
     ** Change order of bytes of all faces. */
    void reorder_faces_bytes_(const std::vector<uint_t>& face_reordering);
    /** @internal
     ** Change order of bytes of all cells. */
    void reorder_cells_bytes_(const std::vector<uint_t>& cell_reordering);

protected:

    void reorder_faces();

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /**
     * Generate edges using the face to node connectivity.
     * @warning This function may be slow and memory-consuming.
     */
    void generate_edges();

    /**
     * Generate faces using the cell to node connectivity.
     * @warning This function may be slow and memory-consuming.
     */
    void generate_faces();

    /** Generate boundary cells to complete face connectivity. */
    void generate_boundary_cells();
};  // class cMesh

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

/* Include iterators. */
#include "MeshIterators.hh"

/** @todo Remove me. */
using cMesh = feathers::cMesh;
#include "libSkunkSparse/SkunkSparseField.hh"

#endif  // ifndef MESH_HH_
