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
#ifndef SKUNK_MESH_HH
#define SKUNK_MESH_HH

#include <SkunkBase.hh>
#include <libSkunkGeom/SkunkGeomMatrix.hh>

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

/** Type of the mesh elements. */
enum class mesh_elem_type_t : byte_t {
    null,
    node,
    segment_2, segment_3,
    triangle_3, triangle_6,
    quadrangle_4, quadrangle_8, quadrangle_9,
    tetrahedron_4, tetrahedron_10,
    pyramid_5, pyramid_14,
    pentahedron_6, pentahedron_15, pentahedron_18,
    hexahedron_8, hexahedron_20, hexahedron_27,
};  // enum class mesh_elem_type_t

/** Mesh element connectivity structure. */
template<mesh_elem_type_t type_t,
         uint_t num_nodes_t, uint_t num_edges_t, uint_t num_faces_t, uint_t num_cells_t>
class mesh_elem_struct_t {
private:
    /** @internal
     ** @warning Order of the fields should not be changed! */
    /** @{ */
    mesh_elem_type_t m_type : 8;
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

public:

    /**************************************************************************/
    /**************************************************************************/

    /** Type of the element. */
    mesh_elem_type_t get_type() const {
        return m_type;
    }

    /**************************************************************************/
    /**************************************************************************/

    /** Number of nodes in the element. */
    uint_t num_nodes() const {
        return m_num_nodes;
    }
    /** Number of edges in the element. */
    uint_t num_edges() const {
        return m_num_edges;
    }
    /** Number of faces in the element. */
    uint_t num_faces() const {
        return m_num_faces;
    }
    /** Number of cells in the element. */
    uint_t num_cells() const {
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

    /** Assign node indices. */
    void set_nodes(std::initializer_list<uint_t> nodes) {
        std::copy(nodes.begin(), nodes.end(), begin_node());
    }
    /** Assign edge indices. */
    void set_edges(std::initializer_list<uint_t> edges) {
        std::copy(edges.begin(), edges.end(), begin_edge());
    }
    /** Assign face indices. */
    void set_faces(std::initializer_list<uint_t> faces) {
        std::copy(faces.begin(), faces.end(), begin_face());
    }
    /** Assign cell indices. */
    void set_cells(std::initializer_list<uint_t> cells) {
        std::copy(cells.begin(), cells.end(), begin_cell());
    }

    /**************************************************************************/
    /**************************************************************************/

    /** Erase the node from element.
     ** @warning Element type would be set to unknown. */
    void erase_node(uint_t node_loc) {
        SKUNK_ASSERT(node_loc < num_nodes());
        std::rotate(begin_face() + node_loc, begin_face() + node_loc + 1, end_cell());
        m_type = mesh_elem_type_t::null;
        m_num_nodes -= 1;
    }
    /** Erase the edge from element.
     ** @warning Element type would be set to unknown. */
    void erase_edge(uint_t edge_loc) {
        SKUNK_ASSERT(edge_loc < num_edges());
        std::rotate(begin_face() + edge_loc, begin_face() + edge_loc + 1, end_cell());
        m_type = mesh_elem_type_t::null;
        m_num_edges -= 1;
    }
    /** Erase the face from element.
     ** @warning Element type would be set to unknown. */
    void erase_face(uint_t face_loc) {
        SKUNK_ASSERT(face_loc < num_faces());
        std::rotate(begin_face() + face_loc, begin_face() + face_loc + 1, end_cell());
        m_type = mesh_elem_type_t::null;
        m_num_faces -= 1;
    }
    /** Erase the cell from element.
     ** @warning Element type would be set to unknown. */
    void erase_cell(uint_t cell_loc) {
        SKUNK_ASSERT(cell_loc < num_cells());
        std::rotate(begin_cell() + cell_loc, begin_cell() + cell_loc + 1, end_cell());
        m_type = mesh_elem_type_t::null;
        m_num_cells -= 1;
    }
};  // class mesh_elem_struct_t

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

template<mesh_elem_type_t type_t>
class mesh_node_struct_t;

/** Generic node element. */
using mesh_node_t
    = mesh_node_struct_t<mesh_elem_type_t::null>;

/** Node element. */
using mesh_node1_t
    = mesh_node_struct_t<mesh_elem_type_t::node>;

/** Base node element. */
class mesh_node_struct_base_t {
private:
    uint_t m_mark;
    vec3_t m_position;

protected:
    /** Construct base node element. */
    explicit mesh_node_struct_base_t(uint_t mark = 0)
        : m_mark(mark) {
    }

public:
    /** Get node mark. */
    uint_t get_mark() const {
        return m_mark;
    }
    /** Set node mark. */
    void set_mark(uint_t mark) {
        m_mark = mark;
    }

    /** Get node position. */
    const vec3_t& get_position() const {
        return m_position;
    }
    /** Set node position. */
    void set_position(const vec3_t& position) {
        m_position = position;
    }
};  // class mesh_node_struct_base_t

/** Node element. */
template<mesh_elem_type_t type_t>
class mesh_node_struct_t final : public mesh_node_struct_base_t,
                                 public mesh_elem_struct_t<type_t, 0, 0, 0, 0> {
public:
    /** Construct a node element. */
    explicit mesh_node_struct_t(uint_t mark = 0)
        : mesh_node_struct_base_t(mark) {
    }
    /** Construct a node element with position. */
    explicit mesh_node_struct_t(const vec3_t& position, uint_t mark = 0)
        : mesh_node_struct_base_t(mark) {
        this->set_position(position);
    }
};  // class mesh_node_struct_t

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

template<mesh_elem_type_t type_t, uint_t num_nodes_t>
class mesh_edge_struct_t;

/** Generic edge element. */
using mesh_edge_t
    = mesh_edge_struct_t<mesh_elem_type_t::null, 0>;

/** Dummy "node" edge element.
 ** Acts as a placeholder fo 1D/2D meshes. */
using mesh_edge_node1_t
    = mesh_edge_struct_t<mesh_elem_type_t::node, 1>;

/** Linear segment edge element. */
using mesh_edge_segment2_t
    = mesh_edge_struct_t<mesh_elem_type_t::segment_2, 2>;
/** Quadratic segment edge element. */
using mesh_edge_segment3_t
    = mesh_edge_struct_t<mesh_elem_type_t::segment_3, 3>;

/** Base edge element. */
class mesh_edge_struct_base_t {
private:
    uint_t m_mark;
    real_t m_length;
    vec3_t m_direction;

protected:
    /** Construct base edge element. */
    explicit mesh_edge_struct_base_t(uint_t mark = 0)
        : m_mark(mark), m_length() {
    }

public:
    /** Get edge mark. */
    uint_t get_mark() const {
        return m_mark;
    }
    /** Set edge mark. */
    void set_mark(uint_t mark) {
        m_mark = mark;
    }

    /** Get edge length. */
    real_t get_length() const {
        return m_length;
    }
    /** Set edge length. */
    void set_length(real_t length) {
        m_length = length;
    }

    /** Get edge direction. */
    const vec3_t& get_direction() const {
        return m_direction;
    }
    /** Set edge direction. */
    void set_direction(const vec3_t& direction) {
        m_direction = direction;
    }
};  // class mesh_edge_struct_base_t

/** Edge element. */
template<mesh_elem_type_t type_t, uint_t num_nodes_t>
class mesh_edge_struct_t : public mesh_edge_struct_base_t,
                           public mesh_elem_struct_t<type_t, num_nodes_t, 0, 0, 0> {
public:
    /** Construct an edge element. */
    explicit mesh_edge_struct_t(uint_t mark = 0)
        : mesh_edge_struct_base_t(mark) {
    }
    /** Construct an edge element with nodes connectivity. */
    /** @{ */
    mesh_edge_struct_t(std::initializer_list<uint_t> nodes, uint_t mark = 0)
        : mesh_edge_struct_base_t(mark),
          mesh_elem_struct_t<type_t, num_nodes_t, 0, 0, 0>(nodes) {
    }
    template<typename node_iter_t>
    mesh_edge_struct_t(node_iter_t begin_nodes, node_iter_t end_nodes, uint_t mark = 0)
        : mesh_edge_struct_base_t(mark),
          mesh_elem_struct_t<type_t, num_nodes_t, 0, 0, 0>(begin_nodes, end_nodes) {
    }
    /** @} */
};  // class mesh_edge_struct_t

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

template<mesh_elem_type_t type_t, uint_t num_nodes_t, uint_t num_edges_t>
class mesh_face_struct_t;

/** Generic face element. */
using mesh_face_t
    = mesh_face_struct_t<mesh_elem_type_t::null, 0, 0>;

/** Dummy "node" face element.
 ** Acts as a placeholder for 1D meshes. */
using mesh_face_node1_t
    = mesh_face_struct_t<mesh_elem_type_t::node, 1, 1>;

/** Linear segment face element for 2D meshes. */
using mesh_face_segment2_t
    = mesh_face_struct_t<mesh_elem_type_t::segment_2, 2, 2>;
/** Quadratic segment face element for 2D meshes. */
using mesh_face_segment3_t
    = mesh_face_struct_t<mesh_elem_type_t::segment_3, 3, 2>;

/** Linear triangle face element. */
using mesh_face_triangle3_t
    = mesh_face_struct_t<mesh_elem_type_t::triangle_3, 3, 3>;
/** Quadratic triangle face element. */
using mesh_face_triangle6_t
    = mesh_face_struct_t<mesh_elem_type_t::triangle_6, 6, 3>;

/** Linear quadrangle face element. */
using mesh_face_quadrangle4_t
    = mesh_face_struct_t<mesh_elem_type_t::quadrangle_4, 4, 4>;
/** Quadratic quadrangle face element. */
/** @{ */
using mesh_face_quadrangle8_t
    = mesh_face_struct_t<mesh_elem_type_t::quadrangle_8, 8, 4>;
using mesh_face_quadrangle9_t
    = mesh_face_struct_t<mesh_elem_type_t::quadrangle_9, 9, 4>;
/** @} */

/** Base face element. */
class mesh_face_struct_base_t {
private:
    uint_t m_mark;
    real_t m_area;
    vec3_t m_normal, m_center;

protected:
    /** Construct base face element. */
    explicit mesh_face_struct_base_t(uint_t mark = 0)
        : m_mark(mark), m_area() {
    }

public:
    /** Get face mark. */
    uint_t get_mark() const {
        return m_mark;
    }
    /** Set face mark. */
    void set_mark(uint_t mark) {
        m_mark = mark;
    }

    /** Get face area/length. */
    real_t get_area() const {
        return m_area;
    }
    /** Set face area/length. */
    void set_area(real_t area) {
        m_area = area;
    }

    /** Get face normal. */
    const vec3_t& get_normal() const {
        return m_normal;
    }
    /** Set face normal. */
    void set_normal(const vec3_t& normal) {
        m_normal = normal;
    }

    /** Get face barycenter. */
    const vec3_t& get_center_position() const {
        return m_center;
    }
    /** Set face barycenter. */
    void set_center_position(const vec3_t& center) {
        m_center = center;
    }
};  // class mesh_face_struct_base_t

/** Face element. */
template<mesh_elem_type_t type_t, uint_t num_nodes_t, uint_t num_edges_t>
class mesh_face_struct_t : public mesh_face_struct_base_t,
                           public mesh_elem_struct_t<type_t, num_nodes_t, num_edges_t, 0, 2> {
public:
    /** Construct a face element. */
    explicit mesh_face_struct_t(uint_t mark = 0)
        : mesh_face_struct_base_t(mark) {
    }
    /** Construct a face element with nodes connectivity. */
    /** @{ */
    mesh_face_struct_t(std::initializer_list<uint_t> nodes, uint_t mark = 0)
        : mesh_face_struct_base_t(mark),
          mesh_elem_struct_t<type_t, num_nodes_t, num_edges_t, 0, 2>(nodes) {
    }
    template<typename node_iter_t>
    mesh_face_struct_t(node_iter_t begin_nodes, node_iter_t end_nodes, uint_t mark = 0)
        : mesh_face_struct_base_t(mark),
          mesh_elem_struct_t<type_t, num_nodes_t, num_edges_t, 0, 2>(begin_nodes, end_nodes) {
    }
    /** @} */

    /**************************************************************************/
    /**************************************************************************/

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

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

template<mesh_elem_type_t type_t, uint_t num_nodes_t, uint_t num_faces_t>
class mesh_cell_struct_t;

/** Generic cell element. */
using mesh_cell_t
    = mesh_cell_struct_t<mesh_elem_type_t::null, 0, 0>;

/** Linear segment cell element for 1D meshes. */
using mesh_cell_segment2_t
    = mesh_cell_struct_t<mesh_elem_type_t::segment_2, 2, 2>;
/** Quadratic segment cell element for 1D meshes. */
using mesh_cell_segment3_t
    = mesh_cell_struct_t<mesh_elem_type_t::segment_3, 3, 2>;

/** Linear triangle cell element for 2D meshes. */
using mesh_cell_triangle3_t
    = mesh_cell_struct_t<mesh_elem_type_t::triangle_3, 3, 3>;
/** Quadratic triangle cell element for 2D meshes. */
using mesh_cell_triangle6_t
    = mesh_cell_struct_t<mesh_elem_type_t::triangle_6, 6, 3>;

/** Linear quadrangle cell element for 2D meshes. */
using mesh_cell_quadrangle4_t
    = mesh_cell_struct_t<mesh_elem_type_t::quadrangle_4, 4, 4>;
/** Quadratic quadrangle cell element for 2D meshes. */
/** @{ */
using mesh_cell_quadrangle8_t
    = mesh_cell_struct_t<mesh_elem_type_t::quadrangle_8, 8, 4>;
using mesh_cell_quadrangle9_t
    = mesh_cell_struct_t<mesh_elem_type_t::quadrangle_9, 9, 4>;
/** @} */

/** Linear tetrahedron cell element. */
using mesh_cell_tetrahedron4_t
    = mesh_cell_struct_t<mesh_elem_type_t::tetrahedron_4, 4, 4>;
/** Quadratic tetrahedron cell element. */
using mesh_cell_tetrahedron10_t
    = mesh_cell_struct_t<mesh_elem_type_t::tetrahedron_10, 10, 4>;

/** Linear pyramid cell element. */
using mesh_cell_pyramid5_t
    = mesh_cell_struct_t<mesh_elem_type_t::pyramid_5, 5, 5>;
/** Quadratic pyramid cell element. */
using mesh_cell_pyramid14_t
    = mesh_cell_struct_t<mesh_elem_type_t::pyramid_14, 14, 5>;

/** Linear pentahedron cell element. */
using mesh_cell_pentahedron6_t
    = mesh_cell_struct_t<mesh_elem_type_t::pentahedron_6, 6, 5>;
/** Quadratic pentahedron cell element. */
/** @{ */
using mesh_cell_pentahedron15_t
    = mesh_cell_struct_t<mesh_elem_type_t::pentahedron_15, 15, 5>;
using mesh_cell_pentahedron18_t
    = mesh_cell_struct_t<mesh_elem_type_t::pentahedron_18, 18, 5>;
/** @} */

/** Linear hexahedron cell element. */
using mesh_cell_hexahedron8_t
    = mesh_cell_struct_t<mesh_elem_type_t::hexahedron_8, 8, 6>;
/** Quadratic hexahedron cell element. */
/** @{ */
using mesh_cell_hexahedron20_t
    = mesh_cell_struct_t<mesh_elem_type_t::hexahedron_20, 20, 8>;
using mesh_cell_hexahedron27_t
    = mesh_cell_struct_t<mesh_elem_type_t::hexahedron_27, 27, 8>;
/** @} */

/** Base cell element. */
class mesh_cell_struct_base_t {
private:
    uint_t m_mark;
    real_t m_volume;
    vec3_t m_center;

protected:
    /** Construct base cell element. */
    explicit mesh_cell_struct_base_t(uint_t mark = 0)
        : m_mark(mark), m_volume() {
    }

public:
    /** Get cell mark. */
    int_t get_mark() const {
        return m_mark;
    }
    /** Set cell mark. */
    void set_mark(int_t mark) {
        m_mark = mark;
    }

    /** Get cell volume/area/length. */
    real_t get_volume() const {
        return m_volume;
    }
    /** Set cell volume/area/length. */
    void set_volume(real_t volume) {
        m_volume = volume;
    }

    /** Get cell barycenter. */
    const vec3_t& get_center_position() const {
        return m_center;
    }
    /** Set cell barycenter. */
    void set_center_position(const vec3_t& center) {
        m_center = center;
    }
};  // class mesh_cell_struct_base_t

/** Cell element. */
template<mesh_elem_type_t type_t, uint_t num_nodes_t, uint_t num_faces_t>
class mesh_cell_struct_t : public mesh_cell_struct_base_t,
                           public mesh_elem_struct_t<type_t, num_nodes_t, 0, num_faces_t, 0> {
public:
    /** Construct a cell element. */
    explicit mesh_cell_struct_t(uint_t mark = 0)
        : mesh_cell_struct_base_t(mark) {
    }
    /** Construct a cell element with nodes connectivity. */
    /** @{ */
    mesh_cell_struct_t(std::initializer_list<uint_t> nodes, uint_t mark = 0)
        : mesh_cell_struct_base_t(mark),
          mesh_elem_struct_t<type_t, num_nodes_t, 0, num_faces_t, 0>(nodes) {
    }
    template<typename node_iter_t>
    mesh_cell_struct_t(node_iter_t begin_nodes, node_iter_t end_nodes, uint_t mark = 0)
        : mesh_cell_struct_base_t(mark),
          mesh_elem_struct_t<type_t, num_nodes_t, 0, num_faces_t, 0>(begin_nodes, end_nodes) {
    }
    /** @} */
};  // class mesh_cell_struct_t

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

/** Hybrid unstructured finite element mesh. */
class mesh_t : public TObject<mesh_t> {
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
    std::vector<uint_t> m_marked_node_ranges{0};
    std::vector<uint_t> m_marked_edge_ranges{0};
    std::vector<uint_t> m_marked_face_ranges{0};
    std::vector<uint_t> m_marked_cell_ranges{0};

public:

    /**************************************************************************/
    /**************************************************************************/

    /** Initialize an empty mesh. */
    explicit mesh_t(uint_t dim = 0)
        : m_dim(dim) {
        SKUNK_ASSERT(0 <= m_dim && m_dim <= 3);
    }

    /**************************************************************************/
    /**************************************************************************/

    bool read_su2(const char* path);
    bool read_tetgen(const char* path);
    bool read_triangle(const char* path);

    /**************************************************************************/
    /**************************************************************************/

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
    mesh_node_t& get_node(uint_t node_ind) {
        SKUNK_ASSERT(node_ind < num_nodes());
        return reinterpret_cast<mesh_node_t&>(m_node_storage[m_node_offsets[node_ind]]);
    }
    const mesh_node_t& get_node(uint_t node_ind) const {
        SKUNK_ASSERT(node_ind < num_nodes());
        return reinterpret_cast<const mesh_node_t&>(m_node_storage[m_node_offsets[node_ind]]);
    }
    /** @} */
    /** Get edge element at global index. */
    /** @{ */
    mesh_edge_t& get_edge(uint_t edge_ind) {
        SKUNK_ASSERT(edge_ind < num_edges());
        return reinterpret_cast<mesh_edge_t&>(m_edge_storage[m_edge_offsets[edge_ind]]);
    }
    const mesh_edge_t& get_edge(uint_t edge_ind) const {
        SKUNK_ASSERT(edge_ind < num_edges());
        return reinterpret_cast<const mesh_edge_t&>(m_edge_storage[m_edge_offsets[edge_ind]]);
    }
    /** @} */
    /** Get face element at global index. */
    /** @{ */
    mesh_face_t& get_face(uint_t face_ind) {
        SKUNK_ASSERT(face_ind < num_faces());
        return reinterpret_cast<mesh_face_t&>(m_face_storage[m_face_offsets[face_ind]]);
    }
    const mesh_face_t& get_face(uint_t face_ind) const {
        SKUNK_ASSERT(face_ind < num_faces());
        return reinterpret_cast<const mesh_face_t&>(m_face_storage[m_face_offsets[face_ind]]);
    }
    /** @} */
    /** Get cell element at global index. */
    /** @{ */
    mesh_cell_t& get_cell(uint_t cell_ind) {
        SKUNK_ASSERT(cell_ind < num_cells());
        return reinterpret_cast<mesh_cell_t&>(m_cell_storage[m_cell_offsets[cell_ind]]);
    }
    const mesh_cell_t& get_cell(uint_t cell_ind) const {
        SKUNK_ASSERT(cell_ind < num_cells());
        return reinterpret_cast<const mesh_cell_t&>(m_cell_storage[m_cell_offsets[cell_ind]]);
    }
    /** @} */

    /**************************************************************************/
    /**************************************************************************/

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

    /** Index of the first node with a given mark. */
    uint_t begin_node(uint_t mark) const {
        SKUNK_ASSERT(mark < num_node_marks());
        return m_marked_node_ranges[mark];
    }
    /** Index of the first edge with a given mark. */
    uint_t begin_edge(uint_t mark) const {
        SKUNK_ASSERT(mark < num_edge_marks());
        return m_marked_edge_ranges[mark];
    }
    /** Index of the first face with a given mark. */
    uint_t begin_face(uint_t mark) const {
        SKUNK_ASSERT(mark < num_face_marks());
        return m_marked_face_ranges[mark];
    }
    /** Index of the first cell with a given mark. */
    uint_t begin_cell(uint_t mark) const {
        SKUNK_ASSERT(mark < num_cell_marks());
        return m_marked_cell_ranges[mark];
    }

    /** Index of a node after the last node with a given mark. */
    uint_t end_node(uint_t mark) const {
        SKUNK_ASSERT(mark < num_node_marks());
        return m_marked_node_ranges[mark + 1];
    }
    /** Index of an edge after the last edge with a given mark. */
    uint_t end_edge(uint_t mark) const {
        SKUNK_ASSERT(mark < num_edge_marks());
        return m_marked_edge_ranges[mark + 1];
    }
    /** Index of a face after the last face with a given mark. */
    uint_t end_face(uint_t mark) const {
        SKUNK_ASSERT(mark < num_face_marks());
        return m_marked_face_ranges[mark + 1];
    }
    /** Index a cell after of the last cell with a given mark. */
    uint_t end_cell(uint_t mark) const {
        SKUNK_ASSERT(mark < num_cell_marks());
        return m_marked_cell_ranges[mark + 1];
    }

    /** Get face element index at marker index. */
    uint_t get_marked_face_index(uint_t face_mark_ind, uint_t mark) const {
        return begin_face(mark) + face_mark_ind;
    }
    /** Get cell element index at marker index. */
    uint_t get_marked_cell_index(uint_t cell_mark_ind, uint_t mark) const {
        return begin_cell(mark) + cell_mark_ind;
    }

    /** Number of marked faces in the mesh. */
    uint_t num_marked_faces(uint_t mark) const {
        return end_face(mark) - begin_face(mark);
    }
    /** Number of marked cells in the mesh. */
    uint_t num_marked_cells(uint_t mark) const {
        return end_cell(mark) - begin_cell(mark);
    }

    /**************************************************************************/
    /**************************************************************************/

    /** Get node position. */
    const vec3_t& get_node_position(uint_t node_ind) const {
        return get_node(node_ind).get_position();
    }

    /** Get edge length. */
    real_t get_edge_length(uint_t edge_ind) const {
        return get_edge(edge_ind).get_length();
    }
    /** Get edge direction. */
    const vec3_t& get_edge_direction(uint_t edge_ind) const {
        return get_edge(edge_ind).get_direction();
    }

    /** Get face area/length. */
    real_t get_face_area(uint_t face_ind) const {
        return get_face(face_ind).get_area();
    }
    /** Get face normal. */
    const vec3_t& get_face_normal(uint_t face_ind) const {
        return get_face(face_ind).get_normal();
    }
    /** Get face barycenter. */
    const vec3_t& get_face_center_position(uint_t face_ind) const {
        return get_face(face_ind).get_center_position();
    }

    /** Get cell volume/area/length. */
    real_t get_cell_volume(uint_t cell_ind) const {
        return get_cell(cell_ind).get_volume();
    }
    /** Get cell barycenter. */
    const vec3_t& get_cell_center_position(uint_t cell_ind) const {
        return get_cell(cell_ind).get_center_position();
    }

    /** Get minimal edge length. */
    SKUNK_EXTERN real_t get_min_edge_length() const;
    /** Get maximal edge length. */
    SKUNK_EXTERN real_t get_max_edge_length() const;

    /** Get minimal face area. */
    SKUNK_EXTERN real_t get_min_face_area() const;
    /** Get maximal face area. */
    SKUNK_EXTERN real_t get_max_face_area() const;

    /** Get minimal cell volume. */
    SKUNK_EXTERN real_t get_min_cell_volume() const;
    /** Get maximal cell volume. */
    SKUNK_EXTERN real_t get_max_cell_volume() const;

    /** Compute mesh orthogonality.
     ** Mesh orthogonality is defined as a. */
    SKUNK_EXTERN real_t get_orthogonality() const;

    /**************************************************************************/
    /**************************************************************************/

    /** Rescale the mesh. */
    SKUNK_EXTERN void scale(real_t scale);

    /**************************************************************************/
    /**************************************************************************/

protected:

    /** Insert a new node into the mesh.
     ** @param node Edge nodes.
     ** @returns Index of the inserted node. */
    SKUNK_EXTERN virtual uint_t insert_node(const mesh_node1_t& node);

    /** Insert a new edge into the mesh.
     ** @param edge Edge nodes.
     ** @returns Index of the inserted edge. */
    /** @{ */
    SKUNK_EXTERN virtual uint_t insert_edge(const mesh_edge_node1_t& edge);
    SKUNK_EXTERN virtual uint_t insert_edge(const mesh_edge_segment2_t& edge);
    SKUNK_EXTERN virtual uint_t insert_edge(const mesh_edge_segment3_t& edge);
    /** @} */

    /** Insert a new edge into the mesh.
     ** Type of the edge is detected by the set of nodes and dimension.
     ** @param edge_nodes Set of the edge nodes.
     ** @param mark Boundary mark.
     ** @param dim Overridden dimension.
     ** @returns Index of the inserted edge. */
    SKUNK_EXTERN virtual uint_t insert_edge(const std::vector<uint_t>& edge_nodes, uint_t mark = 0, uint_t dim = 0);

    /** Insert a new face into the mesh.
     ** @param face Face nodes.
     ** @returns Index of the inserted face. */
    /** @{ */
    SKUNK_EXTERN virtual uint_t insert_face(const mesh_face_node1_t& face);
    SKUNK_EXTERN virtual uint_t insert_face(const mesh_face_segment2_t& face);
    SKUNK_EXTERN virtual uint_t insert_face(const mesh_face_segment3_t& face);
    SKUNK_EXTERN virtual uint_t insert_face(const mesh_face_triangle3_t& face);
    SKUNK_EXTERN virtual uint_t insert_face(const mesh_face_triangle6_t& face);
    SKUNK_EXTERN virtual uint_t insert_face(const mesh_face_quadrangle4_t& face);
    SKUNK_EXTERN virtual uint_t insert_face(const mesh_face_quadrangle8_t& face);
    SKUNK_EXTERN virtual uint_t insert_face(const mesh_face_quadrangle9_t& face);
    /** @} */

    /** Insert a new face into the mesh.
     ** Type of the face is detected by the set of nodes and dimension.
     ** @param face_nodes Set of the face nodes.
     ** @param mark Boundary mark.
     ** @param dim Overridden dimension.
     ** @returns Index of the inserted face. */
    SKUNK_EXTERN virtual uint_t insert_face(const std::vector<uint_t>& face_nodes, uint_t mark = 0, uint_t dim = 0);

    /** Insert a new cell into the mesh.
     ** @param cell Cell nodes.
     ** @returns Index of the inserted face. */
    /** @{ */
    SKUNK_EXTERN virtual uint_t insert_cell(const mesh_cell_segment2_t& cell);
    SKUNK_EXTERN virtual uint_t insert_cell(const mesh_cell_segment3_t& cell);
    SKUNK_EXTERN virtual uint_t insert_cell(const mesh_cell_triangle3_t& cell);
    SKUNK_EXTERN virtual uint_t insert_cell(const mesh_cell_triangle6_t& cell);
    SKUNK_EXTERN virtual uint_t insert_cell(const mesh_cell_quadrangle4_t& cell);
    SKUNK_EXTERN virtual uint_t insert_cell(const mesh_cell_quadrangle8_t& cell);
    SKUNK_EXTERN virtual uint_t insert_cell(const mesh_cell_quadrangle9_t& cell);
    SKUNK_EXTERN virtual uint_t insert_cell(const mesh_cell_tetrahedron4_t&  cell);
    SKUNK_EXTERN virtual uint_t insert_cell(const mesh_cell_tetrahedron10_t& cell);
    SKUNK_EXTERN virtual uint_t insert_cell(const mesh_cell_pyramid5_t&  cell);
    SKUNK_EXTERN virtual uint_t insert_cell(const mesh_cell_pyramid14_t& cell);
    SKUNK_EXTERN virtual uint_t insert_cell(const mesh_cell_pentahedron6_t&  cell);
    SKUNK_EXTERN virtual uint_t insert_cell(const mesh_cell_pentahedron15_t& cell);
    SKUNK_EXTERN virtual uint_t insert_cell(const mesh_cell_pentahedron18_t& cell);
    SKUNK_EXTERN virtual uint_t insert_cell(const mesh_cell_hexahedron8_t&  cell);
    SKUNK_EXTERN virtual uint_t insert_cell(const mesh_cell_hexahedron20_t& cell);
    SKUNK_EXTERN virtual uint_t insert_cell(const mesh_cell_hexahedron27_t& cell);
    /** @} */

    /** Insert a new cell into the mesh.
     ** Type of the cell is detected by the set of nodes and dimension.
     ** @param cell_nodes Set of the cell nodes.
     ** @param mark Boundary mark.
     ** @param dim Overridden dimension.
     ** @returns Index of the inserted cell. */
    SKUNK_EXTERN virtual uint_t insert_cell(const std::vector<uint_t>& cell_nodes, uint_t mark = 0, uint_t dim = 0);

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

    /**************************************************************************/
    /**************************************************************************/

protected:

    /** Erase the node.
     ** @param replacement_node_ind Node to connect links of adjacent elements. */
    SKUNK_EXTERN void erase_node(uint_t node_ind, uint_t replacement_node_ind = npos);
    /** Erase the edge.
     ** @param replacement_edge_ind Edge to connect links of adjacent elements. */
    SKUNK_EXTERN void erase_edge(uint_t edge_ind, uint_t replacement_edge_ind = npos);
    /** Erase the face.
     ** @param replacement_face_ind Face to connect links of adjacent elements. */
    SKUNK_EXTERN void erase_face(uint_t face_ind, uint_t replacement_face_ind = npos);
    /** Erase the cell.
     ** @param replacement_cell_ind Cell to connect links of adjacent elements. */
    SKUNK_EXTERN void erase_cell(uint_t cell_ind, uint_t replacement_cell_ind = npos);

private:

    /** @internal
     ** Deallocate and erase the node. */
    void deallocate_node_(uint_t node_ind);
    /** @internal
     ** Deallocate and erase the edge. */
    void deallocate_edge_(uint_t edge_ind);
    /** @internal
     ** Deallocate and erase the face. */
    void deallocate_face_(uint_t face_ind);
    /** @internal
     ** Deallocate and erase the cell. */
    void deallocate_cell_(uint_t cell_ind);

protected:

    /** Erase edges that's length is less than specified tolerance. */
    SKUNK_EXTERN void erase_degenerate_edges(real_t min_edge_length);
    /** Erase faces that's area is less than specified tolerance. */
    SKUNK_EXTERN void erase_degenerate_faces(real_t min_face_area);

    /**************************************************************************/
    /**************************************************************************/

protected:

    /** Change order of two nodes. */
    SKUNK_EXTERN void reorder_nodes(uint_t first_node_ind, uint_t second_node_ind);
    /** Change order of two edges. */
    SKUNK_EXTERN void reorder_edges(uint_t first_edge_ind, uint_t second_edge_ind);
    /** Change order of two faces. */
    SKUNK_EXTERN void reorder_faces(uint_t first_face_ind, uint_t second_face_ind);
    /** Change order of two cells. */
    SKUNK_EXTERN void reorder_cells(uint_t first_cell_ind, uint_t second_cell_ind);

private:

    /** @internal
     ** Swap bytes of two nodes. */
    void swap_node_bytes_(uint_t first_node_ind, uint_t second_node_ind);
    /** @internal
     ** Swap bytes of two nodes. */
    void swap_edge_bytes_(uint_t first_edge_ind, uint_t second_edge_ind);
    /** @internal
     ** Swap bytes of two nodes. */
    void swap_face_bytes_(uint_t first_face_ind, uint_t second_face_ind);
    /** @internal
     ** Swap bytes of two nodes. */
    void swap_cell_bytes_(uint_t first_cell_ind, uint_t second_cell_ind);

protected:

    /** Change order of all nodes. */
    SKUNK_EXTERN void reorder_nodes(const std::vector<uint_t>& node_reordering);
    /** Change order of all edges. */
    SKUNK_EXTERN void reorder_edges(const std::vector<uint_t>& edge_reordering);
    /** Change order of all faces. */
    SKUNK_EXTERN void reorder_faces(const std::vector<uint_t>& face_reordering);
    /** Change order of all cells. */
    SKUNK_EXTERN void reorder_cells(const std::vector<uint_t>& cell_reordering);

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

    SKUNK_EXTERN void reorder_faces();

    /**************************************************************************/
    /**************************************************************************/

    /** Generate edges using the face to node connectivity.
     ** @warning This function may be slow and memory-consuming. */
    SKUNK_EXTERN void generate_edges();
    /** Generate faces using the cell to node connectivity.
     ** @warning This function may be slow and memory-consuming. */
    SKUNK_EXTERN void generate_faces();

    /** Generate boundary cells to complete face connectivity. */
    SKUNK_EXTERN void generate_boundary_cells();
};  // class mesh_t

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

/* Include iterators. */
#include <libSkunkMesh/SkunkMeshIterator.hh>

#endif  // ifndef SKUNK_MESH_HH

/** @todo Remove me. */
using UMesh = skunk::mesh_t;
#include "libSkunkSparse/SkunkSparseField.hh"
