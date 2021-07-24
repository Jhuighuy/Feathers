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

#include "SkunkBase.hh"
#include "libFeathersUtils/Table.hh"
#include "Element.hh"

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

enum : uint_t {
    /** Local index of the face inner cell. */
    eFaceInnerCell = 0,
    /** Local index of the face outer cell. */
    eFaceOuterCell = 1
};  // enum

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Hybrid unstructured multidimensional mesh.
 */
class cMesh : public tObject<cMesh> {
private:
    uint_t m_dim;

    uint_t m_num_nodes = 0, m_num_edges = 0;
    uint_t m_num_faces = 0, m_num_cells = 0;

    std::vector<uint_t> m_node_marks;
    std::vector<uint_t> m_edge_marks;
    std::vector<uint_t> m_face_marks;
    std::vector<uint_t> m_cell_marks;
    std::vector<uint_t> m_marked_node_ranges{0};
    std::vector<uint_t> m_marked_edge_ranges{0};
    std::vector<uint_t> m_marked_face_ranges{0};
    std::vector<uint_t> m_marked_cell_ranges{0};

    std::vector<eShape> m_edge_shapes;
    std::vector<eShape> m_face_shapes;
    std::vector<eShape> m_cell_shapes;
    // TODO: edge length + direction -> 4D oriented direction.
    // TODO: face area + normal -> 4D oriented area.
    std::vector<vec3_t> m_node_coords;
    std::vector<real_t> m_edge_lengths;
    std::vector<vec3_t> m_edge_directions;
    std::vector<real_t> m_face_areas;
    std::vector<vec3_t> m_face_normals;
    std::vector<vec3_t> m_face_center_coords;
    std::vector<real_t> m_cell_volumes;
    std::vector<vec3_t> m_cell_center_coords;
    real_t m_min_edge_length = 0.0, m_max_edge_length = 0.0;
    real_t m_min_face_area = 0.0, m_max_face_area = 0.0;
    real_t m_min_cell_volume = 0.0, m_max_cell_volume = 0.0;

    cTable m_node_nodes, m_edge_nodes, m_face_nodes, m_cell_nodes;
    cTable m_node_edges, m_edge_edges, m_face_edges, m_cell_edges;
    cTable m_node_faces, m_edge_faces, m_face_faces, m_cell_faces;
    cTable m_node_cells, m_edge_cells, m_face_cells, m_cell_cells;

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
        return m_num_nodes;
    }
    /** Total number of edges in the mesh. */
    uint_t num_edges() const {
        return m_num_edges;
    }
    /** Total number of faces in the mesh. */
    uint_t num_faces() const {
        return m_num_faces;
    }
    /** Total number of cells in the mesh. */
    uint_t num_cells() const {
        return m_num_cells;
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
    uint_t begin_marked_node(uint_t mark) const {
        FEATHERS_ASSERT(mark < num_node_marks());
        return m_marked_node_ranges[mark];
    }
    /** Index of the first edge with a given mark. */
    uint_t begin_marked_edge(uint_t mark) const {
        FEATHERS_ASSERT(mark < num_edge_marks());
        return m_marked_edge_ranges[mark];
    }
    /** Index of the first face with a given mark. */
    uint_t begin_marked_face(uint_t mark) const {
        FEATHERS_ASSERT(mark < num_face_marks());
        return m_marked_face_ranges[mark];
    }
    /** Index of the first cell with a given mark. */
    uint_t begin_marked_cell(uint_t mark) const {
        FEATHERS_ASSERT(mark < num_cell_marks());
        return m_marked_cell_ranges[mark];
    }

    /** Index of a node after the last node with a given mark. */
    uint_t end_marked_node(uint_t mark) const {
        FEATHERS_ASSERT(mark < num_node_marks());
        return m_marked_node_ranges[mark + 1];
    }
    /** Index of an edge after the last edge with a given mark. */
    uint_t end_marked_edge(uint_t mark) const {
        FEATHERS_ASSERT(mark < num_edge_marks());
        return m_marked_edge_ranges[mark + 1];
    }
    /** Index of a face after the last face with a given mark. */
    uint_t end_marked_face(uint_t mark) const {
        FEATHERS_ASSERT(mark < num_face_marks());
        return m_marked_face_ranges[mark + 1];
    }
    /** Index a cell after of the last cell with a given mark. */
    uint_t end_marked_cell(uint_t mark) const {
        FEATHERS_ASSERT(mark < num_cell_marks());
        return m_marked_cell_ranges[mark + 1];
    }

    /** Number of marked cells in the mesh. */
    uint_t num_marked_cells(uint_t mark) const {
        return end_marked_cell(mark) - begin_marked_cell(mark);
    }

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get element shape. */
    /** @{ */
#if FEATHERS_DOXYGEN
    template<typename tTag>
    uint_t get_shape(tTag, uint_t index) const;
#else
    eShape get_shape(tEdgeTag, uint_t edge_index) const {
        FEATHERS_ASSERT(edge_index < num_edges());
        return m_edge_shapes[edge_index];
    }
    eShape get_shape(tFaceTag, uint_t face_index) const {
        FEATHERS_ASSERT(face_index < num_faces());
        return m_face_shapes[face_index];
    }
    eShape get_shape(tCellTag, uint_t cell_index) const {
        FEATHERS_ASSERT(cell_index < num_cells());
        return m_cell_shapes[cell_index];
    }
#endif
    /** @} */
    /** Set element shape. */
    /** @{ */
#if FEATHERS_DOXYGEN
    template<typename tTag>
    uint_t set_shape(tTag, uint_t index, eShape shape) const;
#else
    void set_shape(tEdgeTag, uint_t edge_index, eShape edge_shape) {
        FEATHERS_ASSERT(edge_index < num_edges());
        m_edge_shapes[edge_index] = edge_shape;
    }
    void set_shape(tFaceTag, uint_t face_index, eShape face_shape) {
        FEATHERS_ASSERT(face_index < num_faces());
        m_face_shapes[face_index] = face_shape;
    }
    void set_shape(tCellTag, uint_t cell_index, eShape cell_shape) {
        FEATHERS_ASSERT(cell_index < num_cells());
        m_cell_shapes[cell_index] = cell_shape;
    }
#endif
    /** @} */

    /** Get element object. */
    template<typename tTag>
    std::unique_ptr<const iElement> get_element_object(tTag tag, uint_t index) const {
        // TODO: element factory.
        iElementPtr element(get_shape(tag, index));
        element->assign_nodes(m_num_nodes, m_node_coords.data(),
                              begin_adjacent_node(tag, index),
                              end_adjacent_node(tag, index));
        return std::unique_ptr<const iElement>(element.release());
    }

    /** Compute edge shape properties. */
    void compute_edge_shape_properties();
    /** Compute face shape properties. */
    void compute_face_shape_properties();
    /** Compute cell shape properties. */
    void compute_cell_shape_properties();
    /** Compute all elements shape properties. */
    void compute_all_shape_properties() {
        // TODO:
        //compute_edge_shape_properties();
        compute_face_shape_properties();
        compute_cell_shape_properties();
    }

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
    real_t get_min_edge_length() const {
        return m_min_edge_length;
    }
    /** Get maximal edge length. */
    real_t get_max_edge_length() const {
        return m_max_edge_length;
    }
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
    real_t get_min_face_area() const{
        return m_min_face_area;
    }
    /** Get maximal face area. */
    real_t get_max_face_area() const {
        return m_max_face_area;
    }
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
    real_t get_min_cell_volume() const {
        return m_min_cell_volume;
    }
    /** Get maximal cell volume. */
    real_t get_max_cell_volume() const {
        return m_max_cell_volume;
    }
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

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Pointer to the beginning of the element adjacent nodes. */
    /** @{ */
#if FEATHERS_DOXYGEN
    FEATHERS_CONST_OVERLOAD_T(template<typename tTag>,
        uint_t*, begin_adjacent_node, (tTag, uint_t index), ;)
#else
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_node, (tNodeTag, uint_t node_index), {
        FEATHERS_ASSERT(node_index < num_nodes());
        return m_node_nodes.begin_row(node_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_node, (tEdgeTag, uint_t edge_index), {
        FEATHERS_ASSERT(edge_index < num_edges());
        return m_edge_nodes.begin_row(edge_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_node, (tFaceTag, uint_t face_index), {
        FEATHERS_ASSERT(face_index < num_faces());
        return m_face_nodes.begin_row(face_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_node, (tCellTag, uint_t cell_index), {
        FEATHERS_ASSERT(cell_index < num_cells());
        return m_cell_nodes.begin_row(cell_index);
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
        return m_node_edges.begin_row(node_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_edge, (tEdgeTag, uint_t edge_index), {
        FEATHERS_ASSERT(edge_index < num_edges());
        return m_edge_edges.begin_row(edge_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_edge, (tFaceTag, uint_t face_index), {
        FEATHERS_ASSERT(face_index < num_faces());
        return m_face_edges.begin_row(face_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_edge, (tCellTag, uint_t cell_index), {
        FEATHERS_ASSERT(cell_index < num_cells());
        return m_cell_edges.begin_row(cell_index);
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
        return m_node_faces.begin_row(node_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_face, (tEdgeTag, uint_t edge_index), {
        FEATHERS_ASSERT(edge_index < num_edges());
        return m_edge_faces.begin_row(edge_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_face, (tFaceTag, uint_t face_index), {
        FEATHERS_ASSERT(face_index < num_faces());
        return m_face_faces.begin_row(face_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_face, (tCellTag, uint_t cell_index), {
        FEATHERS_ASSERT(cell_index < num_cells());
        return m_cell_faces.begin_row(cell_index);
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
        return m_node_cells.begin_row(node_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_cell, (tEdgeTag, uint_t edge_index), {
        FEATHERS_ASSERT(edge_index < num_edges());
        return m_edge_cells.begin_row(edge_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_cell, (tFaceTag, uint_t face_index), {
        FEATHERS_ASSERT(face_index < num_faces());
        return m_face_cells.begin_row(face_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, begin_adjacent_cell, (tCellTag, uint_t cell_index), {
        FEATHERS_ASSERT(cell_index < num_cells());
        return m_cell_cells.begin_row(cell_index);
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
        return m_node_nodes.end_row(node_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_node, (tEdgeTag, uint_t edge_index), {
        FEATHERS_ASSERT(edge_index < num_edges());
        return m_edge_nodes.end_row(edge_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_node, (tFaceTag, uint_t face_index), {
        FEATHERS_ASSERT(face_index < num_faces());
        return m_face_nodes.end_row(face_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_node, (tCellTag, uint_t cell_index), {
        FEATHERS_ASSERT(cell_index < num_cells());
        return m_cell_nodes.end_row(cell_index);
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
        return m_node_edges.end_row(node_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_edge, (tEdgeTag, uint_t edge_index), {
        FEATHERS_ASSERT(edge_index < num_edges());
        return m_edge_edges.end_row(edge_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_edge, (tFaceTag, uint_t face_index), {
        FEATHERS_ASSERT(face_index < num_faces());
        return m_face_edges.end_row(face_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_edge, (tCellTag, uint_t cell_index), {
        FEATHERS_ASSERT(cell_index < num_cells());
        return m_cell_edges.end_row(cell_index);
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
        return m_node_faces.end_row(node_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_face, (tEdgeTag, uint_t edge_index), {
        FEATHERS_ASSERT(edge_index < num_edges());
        return m_edge_faces.end_row(edge_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_face, (tFaceTag, uint_t face_index), {
        FEATHERS_ASSERT(face_index < num_faces());
        return m_face_faces.end_row(face_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_face, (tCellTag, uint_t cell_index), {
        FEATHERS_ASSERT(cell_index < num_cells());
        return m_cell_faces.end_row(cell_index);
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
        return m_node_cells.end_row(node_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_cell, (tEdgeTag, uint_t edge_index), {
        FEATHERS_ASSERT(edge_index < num_edges());
        return m_edge_cells.end_row(edge_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_cell, (tFaceTag, uint_t face_index), {
        FEATHERS_ASSERT(face_index < num_faces());
        return m_face_cells.end_row(face_index);
    })
    FEATHERS_CONST_OVERLOAD(uint_t*, end_adjacent_cell, (tCellTag, uint_t cell_index), {
        FEATHERS_ASSERT(cell_index < num_cells());
        return m_cell_cells.end_row(cell_index);
    })
#endif
    /** @} */

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /**
     * Insert a new node into the mesh.
     * @returns Index of the inserted node.
     */
    uint_t insert_node(const vec3_t& node_coords, uint_t mark = 0);

    /**
     * Insert a new edge into the mesh.
     * @returns Index of the inserted edge.
     */
    uint_t insert_edge(eShape edge_shape, const std::vector<uint_t>& edge_nodes, uint_t mark=0);

    /**
     * Insert a new face into the mesh.
     * @returns Index of the inserted face.
     */
    uint_t insert_face(eShape face_shape, const std::vector<uint_t>& face_nodes, uint_t mark=0);

    /**
     * Insert a new cell into the mesh.
     * @returns Index of the inserted cell.
     */
    uint_t insert_cell(eShape cell_shape, const std::vector<uint_t>& cell_nodes, uint_t mark=0);

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

private:
    template<typename tTag>
    void fix_permutation_and_adjacency_(tTag tag, std::vector<uint_t>& permutation);

public:

    /** Change order of all nodes. */
    void permute_nodes(std::vector<uint_t>&& node_permutation);
    /** Change order of all edges. */
    void permute_edges(std::vector<uint_t>&& edge_permutation);
    /** Change order of all faces. */
    void permute_faces(std::vector<uint_t>&& face_permutation);
    /** Change order of all cells. */
    void permute_cells(std::vector<uint_t>&& cell_permutation);

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
};  // class tMesh

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
