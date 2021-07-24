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
#ifndef MESH_ITERATORS_HH_
#define MESH_ITERATORS_HH_

#include <SkunkBase.hh>
#include "Mesh.hh"
#include <libFeathersUtils/Parallel.hh>

namespace feathers {

template<typename tMesh>
class tNodeIterBase;
template<typename tMesh>
class tEdgeIterBase;
template<typename tMesh>
class tFaceIterBase;
template<typename tMesh>
class tCellIterBase;

/**
 * Base element iterator.
 */
template<class tIter, typename tMesh, typename tTag, tTag tag>
class tElementIterBase {
private:
    tMesh* m_mesh;
    uint_t m_index;

protected:
    tElementIterBase(tMesh& mesh, uint_t index):
        m_mesh(&mesh), m_index(index) {
        FEATHERS_ASSERT(m_index != npos);
    }

    template<typename uIter, typename uMesh>
    tElementIterBase( // NOLINT(google-explicit-constructor)
            const tElementIterBase<uIter, uMesh, tTag, tag>& other):
        m_mesh(&other.get_mesh()), m_index(other.get_index()) {
        FEATHERS_ASSERT(m_index != npos);
    }

public:

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get associated mesh. */
    tMesh& get_mesh() const {
        return *m_mesh;
    }
    /** Get associated element index. */
    uint_t get_index() const {
        return m_index;
    }

    /** Cast to index operator. */
    operator uint_t() const { // NOLINT(google-explicit-constructor)
        return m_index;
    }

public:

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Difference operator. */
    int_t operator-(const tIter& other) const {
        FEATHERS_ASSERT(m_mesh == other.m_mesh);
        return m_index - other.m_index;
    }

    /** Equality operator. */
    bool operator==(const tIter& other) const {
        FEATHERS_ASSERT(m_mesh == other.m_mesh);
        return m_index == other.m_index;
    }
    /** Inequality operator. */
    bool operator!=(const tIter& other) const {
        FEATHERS_ASSERT(m_mesh == other.m_mesh);
        return m_index != other.m_index;
    }

    /** Lexicographical less than operator. */
    bool operator<(const tIter& other) const {
        FEATHERS_ASSERT(m_mesh == other.m_mesh);
        return m_index < other.m_index;
    }
    /** Lexicographical less than or equal operator. */
    bool operator<=(const tIter& other) const {
        FEATHERS_ASSERT(m_mesh == other.m_mesh);
        return m_index <= other.m_index;
    }

    /** Lexicographical greater than operator. */
    bool operator>(const tIter& other) const {
        FEATHERS_ASSERT(m_mesh == other.m_mesh);
        return m_index > other.m_index;
    }
    /** Lexicographical greater than or equal operator. */
    bool operator>=(const tIter& other) const {
        FEATHERS_ASSERT(m_mesh == other.m_mesh);
        return m_index >= other.m_index;
    }

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Increment operator. */
    /** @{ */
    tIter& operator++() {
        ++m_index;
        return static_cast<tIter&>(*this);
    }
    const tIter operator++(int) {
        tIter iter(*this);
        return ++*this, iter;
    }
    /** @} */

    /** Decrement operator. */
    /** @{ */
    tIter& operator--() {
        --m_index;
        return static_cast<tIter&>(*this);
    }
    const tIter operator--(int) {
        tIter iter(static_cast<tIter&>(*this));
        return --*this, iter;
    }
    /** @} */

    /** Addition operator. */
    /** @{ */
    tIter& operator+=(int_t offset) {
        m_index += offset;
        return static_cast<tIter&>(*this);
    }
    tIter& operator+=(uint_t offset) {
        m_index += offset;
        return static_cast<tIter&>(*this);
    }
    tIter operator+(int_t offset) const {
        return tIter(static_cast<const tIter&>(*this)) += offset;
    }
    tIter operator+(uint_t offset) const {
        return tIter(static_cast<const tIter&>(*this)) += offset;
    }
    /** @} */

    /** Subtraction operator. */
    /** @{ */
    tIter& operator-=(int_t offset) {
        m_index -= offset;
        return static_cast<tIter&>(*this);
    }
    tIter operator-(int_t offset) const {
        return tIter(static_cast<const tIter&>(*this)) -= offset;
    }
    /** @} */

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Identity dereference operator (used for standard algorithms). */
    /** @{ */
    tIter& operator*() {
        return static_cast<tIter&>(*this);
    }
    const tIter& operator*() const {
        return static_cast<const tIter&>(*this);
    }
    /** @} */

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get mark. */
    uint_t get_mark() const {
        return m_mesh->get_mark(tag, m_index);
    }

    /** Get shape. */
    eShape get_shape() const {
        return m_mesh->get_shape(tag, m_index);
    }
    iShapePtr get_shape_ptr() const {
        return m_mesh->get_shape_ptr(tag, m_index);
    }

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Number of nodes in the element. */
    uint_t num_nodes() const {
        return end(eNodeTag) - begin(eNodeTag);
    }
    /** Number of edges in the element. */
    uint_t num_edges() const {
        return end(eEdgeTag) - begin(eEdgeTag);
    }
    /** Number of faces in the element. */
    uint_t num_faces() const {
        return end(eFaceTag) - begin(eFaceTag);
    }
    /** Number of cells in the element. */
    uint_t num_cells() const {
        return end(eCellTag) - begin(eCellTag);
    }

    /** Pointer to the beginning of the adjacent elements. */
    /** @{ */
#if FEATHERS_DOXYGEN
    template<typename uTag>
    auto begin(uTag);
#else
    auto begin(tNodeTag) const {
        return m_mesh->begin_adjacent_node(tag, m_index);
    }
    auto begin(tEdgeTag) const {
        return m_mesh->begin_adjacent_edge(tag, m_index);
    }
    auto begin(tFaceTag) const {
        return m_mesh->begin_adjacent_face(tag, m_index);
    }
    auto begin(tCellTag) const {
        return m_mesh->begin_adjacent_cell(tag, m_index);
    }
#endif

    /** Pointer to the end of the adjacent elements. */
#if FEATHERS_DOXYGEN
    template<typename uTag>
    auto end(uTag);
#else
    auto end(tNodeTag) const {
        return m_mesh->end_adjacent_node(tag, m_index);
    }
    auto end(tEdgeTag) const {
        return m_mesh->end_adjacent_edge(tag, m_index);
    }
    auto end(tFaceTag) const {
        return m_mesh->end_adjacent_face(tag, m_index);
    }
    auto end(tCellTag) const {
        return m_mesh->end_adjacent_cell(tag, m_index);
    }
#endif

    /** Get adjacent node. */
    auto get_node(uint_t node_local) const {
        FEATHERS_ASSERT(node_local < num_nodes());
        return tNodeIterBase<tMesh>(*m_mesh, begin(eNodeTag)[node_local]);
    }
    /** Get adjacent edge. */
    auto get_edge(uint_t edge_local) const {
        FEATHERS_ASSERT(edge_local < num_edges());
        return tEdgeIterBase<tMesh>(*m_mesh, begin(eEdgeTag)[edge_local]);
    }
    /** Get adjacent face. */
    auto get_face(uint_t face_local) const {
        FEATHERS_ASSERT(face_local < num_faces());
        return tFaceIterBase<tMesh>(*m_mesh, begin(eFaceTag)[face_local]);
    }
    /** Get adjacent cell. */
    auto get_cell(uint_t cell_local) const {
        FEATHERS_ASSERT(cell_local < num_cells());
        return tCellIterBase<tMesh>(*m_mesh, begin(eCellTag)[cell_local]);
    }

    /** Iterate through all connected nodes. */
    template<typename tFunc>
    void for_each_node(tFunc func) const {
        std::for_each(begin(eNodeTag), end(eNodeTag), [&](uint_t node_index) {
            func(tNodeIterBase<tMesh>(*m_mesh, node_index));
        });
    }
    /** Iterate through all connected edges. */
    template<typename tFunc>
    void for_each_edge(tFunc func) const {
        std::for_each(begin(eEdgeTag), end(eEdgeTag), [&](uint_t edge_index) {
            func(tEdgeIterBase<tMesh>(*m_mesh, edge_index));
        });
    }
    /** Iterate through all connected faces. */
    /** @{ */
    template<typename tFunc>
    void for_each_face(tFunc func) const {
        std::for_each(begin(eFaceTag), end(eFaceTag), [&](uint_t face_index) {
            func(tFaceIterBase<tMesh>(*m_mesh, face_index));
        });
    }
    template<typename tFunc>
    void for_each_face_cells(tFunc func) const {
        std::for_each(begin(eFaceTag), end(eFaceTag), [&](uint_t face_index) {
            tFaceIterBase<tMesh> face(*m_mesh, face_index);
            func(face.get_inner_cell(), face.get_outer_cell());
        });
    }
    /** @} */
    /** Iterate through all connected cells. */
    template<typename tFunc>
    void for_each_cell(tFunc func) const {
        std::for_each(begin(eCellTag), end(eCellTag), [&](uint_t cell_index) {
            func(tCellIterBase<tMesh>(*m_mesh, cell_index));
        });
    }
};  // class tElementIterBase

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Base node iterator.
 */
template<typename tMesh>
class tNodeIterBase final :
    public tElementIterBase<tNodeIterBase<tMesh>, tMesh, tNodeTag, eNodeTag> {
public:

    /** Construct base node iterator. */
    explicit tNodeIterBase(tMesh& mesh, uint_t node_index = 0):
        tElementIterBase<tNodeIterBase<tMesh>, tMesh, tNodeTag, eNodeTag>(mesh, node_index) {
    }
    /** Copy (make const) constructor. */
    tNodeIterBase( // NOLINT(google-explicit-constructor)
            const tNodeIterBase<std::remove_const_t<tMesh>>& other):
        tElementIterBase<tNodeIterBase<tMesh>, tMesh, tNodeTag, eNodeTag>(other) {
    }

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get node position. */
    const vec3_t& get_coords() const {
        return this->get_mesh().get_node_coords(this->get_index());
    }
    /** Set node position. */
    template<typename uMesh = tMesh>
    std::enable_if_t<!std::is_const_v<uMesh>> set_coords(const vec3_t& node_coords) {
        this->get_mesh().set_node_coords(this->get_index(), node_coords);
    }
};  // class tNodeIterBase

/**
 * Mesh nodes random-access iterator.
 */
/** @{ */
using tNodeIter = tNodeIterBase<const cMesh>;
using tNodeMutableIter = tNodeIterBase<cMesh>;
/** @} */

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Base edge iterator.
 */
template<typename tMesh>
class tEdgeIterBase final :
    public tElementIterBase<tEdgeIterBase<tMesh>, tMesh, tEdgeTag, eEdgeTag> {
public:

    /** Construct base edge iterator. */
    explicit tEdgeIterBase(tMesh& mesh, uint_t edge_index = 0):
        tElementIterBase<tEdgeIterBase<tMesh>, tMesh, tEdgeTag, eEdgeTag>(mesh, edge_index) {
    }
    /** Copy (make const) constructor. */
    tEdgeIterBase( // NOLINT(google-explicit-constructor)
            const tEdgeIterBase<std::remove_const_t<tMesh>>& other):
        tElementIterBase<tEdgeIterBase<tMesh>, tMesh, tEdgeTag, eEdgeTag>(other) {
    }

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get edge length. */
    real_t get_length() const {
        return this->get_mesh().get_edge_length(this->get_index());
    }
    /** Set edge length. */
    template<typename uMesh = tMesh>
    std::enable_if_t<!std::is_const_v<uMesh>> set_length(real_t edge_length) const {
        this->get_mesh().set_edge_length(this->get_index(), edge_length);
    }

    /** Get edge direction. */
    const vec3_t& get_direction() const {
        return this->get_mesh().get_edge_direction(this->get_index());
    }
    /** Set edge direction. */
    template<typename uMesh = tMesh>
    std::enable_if_t<!std::is_const_v<uMesh>> set_direction(const vec3_t& edge_direction) const {
        this->get_mesh().set_edge_direction(this->get_index(), edge_direction);
    }
};  // class tEdgeIterBase

/**
 * Mesh edges random-access iterator.
 */
/** @{ */
using tEdgeIter = tEdgeIterBase<const cMesh>;
using tEdgeMutableIter = tEdgeIterBase<cMesh>;
/** @} */

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Base face iterator.
 */
template<typename tMesh>
class tFaceIterBase final :
    public tElementIterBase<tFaceIterBase<tMesh>, tMesh, tFaceTag, eFaceTag> {
public:

    /** Construct base face iterator. */
    explicit tFaceIterBase(tMesh& mesh, uint_t face_index = 0):
        tElementIterBase<tFaceIterBase<tMesh>, tMesh, tFaceTag, eFaceTag>(mesh, face_index) {
    }
    /** Copy (make const) constructor. */
    tFaceIterBase( // NOLINT(google-explicit-constructor)
            const tFaceIterBase<std::remove_const_t<tMesh>>& other):
        tElementIterBase<tFaceIterBase<tMesh>, tMesh, tFaceTag, eFaceTag>(other) {
    }

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get connected inner cell. */
    auto get_inner_cell() const {
        return this->get_cell(eFaceInnerCell);
    }
    /** Get connected outer cell. */
    auto get_outer_cell() const {
        return this->get_cell(eFaceOuterCell);
    }

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get face area/length. */
    real_t get_area() const {
        return this->get_mesh().get_face_area(this->get_index());
    }
    /** Set face area/length. */
    template<typename uMesh = tMesh>
    std::enable_if_t<!std::is_const_v<uMesh>> set_area(real_t face_area) const {
        this->get_mesh().set_face_area(this->get_index(), face_area);
    }

    /** Get face normal. */
    const vec3_t& get_normal() const {
        return this->get_mesh().get_face_normal(this->get_index());
    }
    /** Set face normal. */
    template<typename uMesh = tMesh>
    std::enable_if_t<!std::is_const_v<uMesh>> set_normal(const vec3_t& face_normal) const {
        this->get_mesh().set_face_normal(this->get_index(), face_normal);
    }

    /** Get face barycenter. */
    const vec3_t& get_center_coords() const {
        return this->get_mesh().get_face_center_coords(this->get_index());
    }
    /** Set face barycenter. */
    template<typename uMesh = tMesh>
    std::enable_if_t<!std::is_const_v<uMesh>> set_center_coords(const vec3_t& face_center_coords) const {
        this->get_mesh().set_face_center_coords(this->get_index(), face_center_coords);
    }
};  // class tFaceIterBase

/**
 * Mesh faces random-access iterator.
 */
/** @{ */
using tFaceIter = tFaceIterBase<const cMesh>;
using tFaceMutableIter = tFaceIterBase<cMesh>;
/** @} */

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Base cell iterator.
 */
template<typename tMesh>
class tCellIterBase final :
    public tElementIterBase<tCellIterBase<tMesh>, tMesh, tCellTag, eCellTag> {
public:

    /** Construct base cell iterator. */
    explicit tCellIterBase(tMesh& mesh, uint_t cell_index = 0):
        tElementIterBase<tCellIterBase<tMesh>, tMesh, tCellTag, eCellTag>(mesh, cell_index) {
    }
    /** Copy (make const) constructor. */
    tCellIterBase( // NOLINT(google-explicit-constructor)
            const tCellIterBase<std::remove_const_t<tMesh>>& other):
        tElementIterBase<tCellIterBase<tMesh>, tMesh, tCellTag, eCellTag>(other) {
    }

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get cell volume/area/length. */
    real_t get_volume() const {
        return this->get_mesh().get_cell_volume(this->get_index());
    }
    /** Set cell volume/area/length. */
    template<typename uMesh = tMesh>
    std::enable_if_t<!std::is_const_v<uMesh>> set_volume(real_t cell_volume) const {
        this->get_mesh().set_cell_volume(this->get_index(), cell_volume);
    }

    /** Get cell barycenter. */
    const vec3_t& get_center_coords() const {
        return this->get_mesh().get_cell_center_coords(this->get_index());
    }
    /** Set cell barycenter. */
    template<typename uMesh = tMesh>
    std::enable_if_t<!std::is_const_v<uMesh>> set_center_coords(const vec3_t& cell_center_coords) const {
        this->get_mesh().set_cell_center_coords(this->get_index(), cell_center_coords);
    }
};  // class tCellIterBase

/**
 * Mesh cells random-access iterator.
 */
/** @{ */
using tCellIter = tCellIterBase<const cMesh>;
using tCellMutableIter = tCellIterBase<cMesh>;
/** @} */

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#define FACE_CELL_FUNC_ \
    ([&](tFaceIterBase<tMesh> face) { \
         func(face.get_inner_cell(), face.get_outer_cell()); \
     })

/** Iterator pointing to the first node with a given mark or the first mark. */
/** @{ */
template<typename tMesh>
auto begin_node(tMesh& mesh) {
    return tNodeIterBase<tMesh>(mesh);
}
template<typename tMesh>
auto begin_node(tMesh& mesh, uint_t mark) {
    return tNodeIterBase<tMesh>(mesh, mesh.begin_marked_node(mark));
}
/** @} */
/** Iterator pointing to the first edge with a given mark or the first mark. */
/** @{ */
template<typename tMesh>
auto begin_edge(tMesh& mesh) {
    return tEdgeIterBase<tMesh>(mesh);
}
template<typename tMesh>
auto begin_edge(tMesh& mesh, uint_t mark) {
    return tEdgeIterBase<tMesh>(mesh, mesh.begin_marked_edge(mark));
}
/** @} */
/** Iterator pointing to the first face with a given mark or the first mark. */
/** @{ */
template<typename tMesh>
auto begin_face(tMesh& mesh) {
    return tFaceIterBase<tMesh>(mesh);
}
template<typename tMesh>
auto begin_face(tMesh& mesh, uint_t mark) {
    return tFaceIterBase<tMesh>(mesh, mesh.begin_marked_face(mark));
}
/** @} */
/** Iterator pointing to the first cell with a given mark or the first mark. */
/** @{ */
template<typename tMesh>
auto begin_cell(tMesh& mesh) {
    return tCellIterBase<tMesh>(mesh);
}
template<typename tMesh>
auto begin_cell(tMesh& mesh, uint_t mark) {
    return tCellIterBase<tMesh>(mesh, mesh.begin_marked_cell(mark));
}
/** @} */

/** Iterator pointing to a node after the last node with a given mark or the last mark. */
/** @{ */
template<typename tMesh>
auto end_node(tMesh& mesh) {
    return tNodeIterBase<tMesh>(mesh, mesh.num_nodes());
}
template<typename tMesh>
auto end_node(tMesh& mesh, uint_t mark) {
    return tNodeIterBase<tMesh>(mesh, mesh.end_marked_node(mark));
}
/** @} */
/** Iterator pointing to an edge after the last edge with a given mark or the last mark. */
/** @{ */
template<typename tMesh>
auto end_edge(tMesh& mesh) {
    return tEdgeIterBase<tMesh>(mesh, mesh.num_edges());
}
template<typename tMesh>
auto end_edge(tMesh& mesh, uint_t mark) {
    return tEdgeIterBase<tMesh>(mesh, mesh.end_marked_edge(mark));
}
/** @} */
/** Iterator pointing to a face after the last face with a given mark or the last mark. */
/** @{ */
template<typename tMesh>
auto end_face(tMesh& mesh) {
    return tFaceIterBase<tMesh>(mesh, mesh.num_faces());
}
template<typename tMesh>
auto end_face(tMesh& mesh, uint_t mark) {
    return tFaceIterBase<tMesh>(mesh, mesh.end_marked_face(mark));
}
/** @} */
/** Iterator pointing to a cell after the last cell with a given mark or the last mark. */
/** @{ */
template<typename tMesh>
auto end_cell(tMesh& mesh) {
    return tCellIterBase<tMesh>(mesh, mesh.num_cells());
}
template<typename tMesh>
auto end_cell(tMesh& mesh, uint_t mark) {
    return tCellIterBase<tMesh>(mesh, mesh.end_marked_cell(mark));
}
/** @} */

/** Iterate through all nodes with a given mark or all marks. */
/** @{ */
template<typename tMesh, typename tFunc>
void for_each_node(tMesh& mesh, tFunc func) {
    for_range(begin_node(mesh), end_node(mesh), func);
}
template<typename tMesh, typename tFunc>
void for_each_node(tMesh& mesh, uint_t mark, tFunc func) {
    for_range(begin_node(mesh, mark), end_node(mesh, mark), func);
}
/** @} */
/** Iterate through all edges with a given mark or all marks. */
/** @{ */
template<typename tMesh, typename tFunc>
void for_each_edge(tMesh& mesh, tFunc func) {
    for_range(begin_edge(mesh), end_edge(mesh), func);
}
template<typename tMesh, typename tFunc>
void for_each_edge(tMesh& mesh, uint_t mark, tFunc func) {
    for_range(begin_edge(mesh, mark), end_edge(mesh, mark), func);
}
/** @} */
/** Iterate through all faces with a given mark or all marks. */
/** @{ */
template<typename tMesh, typename tFunc>
void for_each_face(tMesh& mesh, uint_t mark, tFunc func) {
    for_range(begin_face(mesh, mark), end_face(mesh, mark), func);
}
template<typename tMesh, typename tFunc>
void for_each_face(tMesh& mesh, tFunc func) {
    for_range(begin_face(mesh), end_face(mesh), func);
}
template<typename tMesh, typename tFunc>
void for_each_face_cells(tMesh& mesh, tFunc func) {
    for_range(begin_face(mesh), end_face(mesh), FACE_CELL_FUNC_);
}
template<typename tMesh, typename tFunc>
void for_each_face_cells(tMesh& mesh, uint_t mark, tFunc func) {
    for_range(begin_face(mesh, mark), end_face(mesh, mark), FACE_CELL_FUNC_);
}
/** @} */
/** Iterate through all cells with a given mark or all marks. */
/** @{ */
template<typename tMesh, typename tFunc>
void for_each_cell(tMesh& mesh, tFunc func) {
    for_range(begin_cell(mesh), end_cell(mesh), func);
}
template<typename tMesh, typename tFunc>
void for_each_cell(tMesh& mesh, uint_t mark, tFunc func) {
    for_range(begin_cell(mesh, mark), end_cell(mesh, mark), func);
}
/** @} */

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/** Iterator pointing to the first interior node. */
template<typename tMesh>
auto begin_interior_node(tMesh& mesh) {
    return begin_node(mesh, 0);
}
/** Iterator pointing to the first interior edge. */
template<typename tMesh>
auto begin_interior_edge(tMesh& mesh) {
    return begin_edge(mesh, 0);
}
/** Iterator pointing to the first interior face. */
template<typename tMesh>
auto begin_interior_face(tMesh& mesh) {
    return begin_face(mesh, 0);
}
/** Iterator pointing to the first interior cell. */
template<typename tMesh>
auto begin_interior_cell(tMesh& mesh) {
    return begin_cell(mesh, 0);
}

/** Iterator pointing to a node after the last node. */
template<typename tMesh>
auto end_interior_node(tMesh& mesh) {
    return end_node(mesh, 0);
}
/** Iterator pointing to an edge after the last node. */
template<typename tMesh>
auto end_interior_edge(tMesh& mesh) {
    return end_edge(mesh, 0);
}
/** Iterator pointing to a face after the last node. */
template<typename tMesh>
auto end_interior_face(tMesh& mesh) {
    return end_face(mesh, 0);
}
/** Iterator pointing to a cell after the last node. */
template<typename tMesh>
auto end_interior_cell(tMesh& mesh) {
    return end_cell(mesh, 0);
}

/** Iterate through all interior nodes. */
template<typename tMesh, typename tFunc>
void for_each_interior_node(tMesh& mesh, tFunc func) {
    for_range(begin_interior_node(mesh), end_interior_node(mesh), func);
}
/** Iterate through all interior edges. */
template<typename tMesh, typename tFunc>
void for_each_interior_edge(tMesh& mesh, tFunc func) {
    for_range(begin_interior_edge(mesh), end_interior_edge(mesh), func);
}
/** Iterate through all interior faces. */
/** @{ */
template<typename tMesh, typename tFunc>
void for_each_interior_face(tMesh& mesh, tFunc func) {
    for_range(begin_interior_face(mesh), end_interior_face(mesh), func);
}
template<typename tMesh, typename tFunc>
void for_each_interior_face_cells(tMesh& mesh, tFunc func) {
    for_range(begin_interior_face(mesh), end_interior_face(mesh), FACE_CELL_FUNC_);
}
/** @} */
/** Iterate through all interior cells. */
template<typename tMesh, typename tFunc>
void for_each_interior_cell(tMesh& mesh, tFunc func) {
    for_range(begin_interior_cell(mesh), end_interior_cell(mesh), func);
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/** Iterator pointing to the boundary first node with a given mark or the first boundary mark. */
template<typename tMesh>
auto begin_boundary_node(tMesh& mesh, uint_t mark = 1) {
    FEATHERS_ASSERT(mark >= 1);
    return begin_node(mesh, mark);
}
/** Iterator pointing to the boundary first edge with a given mark or the first boundary mark. */
template<typename tMesh>
auto begin_boundary_edge(tMesh& mesh, uint_t mark = 1) {
    FEATHERS_ASSERT(mark >= 1);
    return begin_edge(mesh, mark);
}
/** Iterator pointing to the boundary first face with a given mark or the first boundary mark. */
template<typename tMesh>
auto begin_boundary_face(tMesh& mesh, uint_t mark = 1) {
    FEATHERS_ASSERT(mark >= 1);
    return begin_face(mesh, mark);
}
/** Iterator pointing to the boundary first cell with a given mark or the first boundary mark. */
template<typename tMesh>
auto begin_boundary_cell(tMesh& mesh, uint_t mark = 1) {
    FEATHERS_ASSERT(mark >= 1);
    return begin_cell(mesh, mark);
}

/** Iterator pointing to a node after the last node with a given or the last boundary mark. */
/** @{ */
template<typename tMesh>
auto end_boundary_node(tMesh& mesh) {
    return end_node(mesh);
}
template<typename tMesh>
auto end_boundary_node(tMesh& mesh, uint_t mark) {
    FEATHERS_ASSERT(mark >= 1);
    return end_node(mesh, mark);
}
/** @} */
/** Iterator pointing to an edge after the last edge with a given or the last boundary mark. */
/** @{ */
template<typename tMesh>
auto end_boundary_edge(tMesh& mesh) {
    return end_edge(mesh);
}
template<typename tMesh>
auto end_boundary_edge(tMesh& mesh, uint_t mark) {
    FEATHERS_ASSERT(mark >= 1);
    return end_edge(mesh, mark);
}
/** @} */
/** Iterator pointing to a face after the last face with a given or the last boundary mark. */
/** @{ */
template<typename tMesh>
auto end_boundary_face(tMesh& mesh) {
    return end_face(mesh);
}
template<typename tMesh>
auto end_boundary_face(tMesh& mesh, uint_t mark) {
    FEATHERS_ASSERT(mark >= 1);
    return end_face(mesh, mark);
}
/** @} */
/** Iterator pointing to a cell after the last cell with a given or the last boundary mark. */
/** @{ */
template<typename tMesh>
auto end_boundary_cell(tMesh& mesh) {
    return end_cell(mesh);
}
template<typename tMesh>
auto end_boundary_cell(tMesh& mesh, uint_t mark) {
    FEATHERS_ASSERT(mark >= 1);
    return end_cell(mesh, mark);
}
/** @} */

/** Iterate through all boundary nodes with a given mark or all boundary marks. */
/** @{ */
template<typename tMesh, typename tFunc>
void for_each_boundary_node(tMesh& mesh, tFunc func) {
    for_range(begin_boundary_node(mesh), end_boundary_node(mesh), func);
}
template<typename tMesh, typename tFunc>
void for_each_boundary_node(tMesh& mesh, uint_t mark, tFunc func) {
    for_range(begin_boundary_node(mesh, mark), end_boundary_node(mesh, mark), func);
}
/** @} */
/** Iterate through all boundary edges with a given mark or all boundary marks. */
/** @{ */
template<typename tMesh, typename tFunc>
void for_each_boundary_edge(tMesh& mesh, tFunc func) {
    for_range(begin_boundary_edge(mesh), end_boundary_edge(mesh), func);
}
template<typename tMesh, typename tFunc>
void for_each_boundary_edge(tMesh& mesh, uint_t mark, tFunc func) {
    for_range(begin_boundary_edge(mesh, mark), end_boundary_edge(mesh, mark), func);
}
/** @} */
/** Iterate through all boundary faces with a given mark or all boundary marks. */
/** @{ */
template<typename tMesh, typename tFunc>
void for_each_boundary_face(tMesh& mesh, tFunc func) {
    for_range(begin_boundary_face(mesh), end_boundary_face(mesh), func);
}
template<typename tMesh, typename tFunc>
void for_each_boundary_face(tMesh& mesh, uint_t mark, tFunc func) {
    for_range(begin_boundary_face(mesh, mark), end_boundary_face(mesh, mark), func);
}
template<typename tMesh, typename tFunc>
void for_each_boundary_face_cells(tMesh& mesh, tFunc func) {
    for_range(begin_boundary_face(mesh), end_boundary_face(mesh), FACE_CELL_FUNC_);
}
template<typename tMesh, typename tFunc>
void for_each_boundary_face_cells(tMesh& mesh, uint_t mark, tFunc func) {
    for_range(begin_boundary_face(mesh, mark), end_boundary_face(mesh, mark), FACE_CELL_FUNC_);
}
/** @} */
/** Iterate through all boundary cells with a given mark or all boundary marks. */
/** @{ */
template<typename tMesh, typename tFunc>
void for_each_boundary_cell(tMesh& mesh, tFunc func) {
    for_range(begin_boundary_cell(mesh), end_boundary_cell(mesh), func);
}
template<typename tMesh, typename tFunc>
void for_each_boundary_cell(tMesh& mesh, uint_t mark, tFunc func) {
    for_range(begin_boundary_cell(mesh, mark), end_boundary_cell(mesh, mark), func);
}
/** @} */

#undef FACE_CELL_FUNC_

} // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif // MESH_ITERATORS_HH_
