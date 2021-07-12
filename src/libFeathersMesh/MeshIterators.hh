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
#include <libSkunkMisc/SkunkMiscParallel.hh>

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

template<class cMesh>
class tNodeIterBase;
template<class cMesh>
class tEdgeIterBase;
template<class cMesh>
class tFaceIterBase;
template<class cMesh>
class tCellIterBase;

/**
 * Base element iterator.
 */
template<class tIter, class cMesh, typename tTag, tTag eTag>
class tElementIterBase {
private:
    cMesh* m_mesh;
    uint_t m_index;

protected:
    template<class cMeshPtr>
    tElementIterBase(const cMeshPtr& mesh_ptr, uint_t index):
        m_mesh(&*mesh_ptr), m_index(index) {
    }

public:

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get associated mesh. */
    cMesh* get_mesh() const {
        return m_mesh;
    }
    /** Get associated element index. */
    uint_t get_index() const {
        return m_index;
    }

    /** Cast to index operator. */
    operator uint_t() const {
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
    const tIter operator++(int_t) {
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
    const tIter operator--(int_t) {
        tIter iter(*this);
        return --*this, iter;
    }
    /** @} */

    /** Addition operator. */
    /** @{ */
    tIter& operator+=(int_t offset) {
        m_index += offset;
        return static_cast<tIter&>(*this);
    }
    tIter operator+(int_t offset) const {
        return tIter(*this) += offset;
    }
    /** @} */

    /** Subtraction operator. */
    /** @{ */
    tIter& operator-=(int_t offset) {
        m_index -= offset;
        return static_cast<tIter&>(*this);
    }
    tIter operator-(int_t offset) const {
        return tIter(*this) -= offset;
    }
    /** @} */

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get mark. */
    uint_t get_mark() const {
        return m_mesh->get_mark(eTag, m_index);
    }

    /** Get shape. */
    eShape get_shape() const {
        return m_mesh->get_shape(eTag, m_index);
    }
    iShapePtr get_shape_ptr() const {
        return m_mesh->get_shape_ptr(eTag, m_index);
    }

    /** Get connected node. */
    auto get_node(uint_t node_local) const {
        FEATHERS_ASSERT(
            node_local < m_mesh->num_adjacent_nodes(eTag, m_index));
        return tNodeIterBase<cMesh>(
            m_mesh, m_mesh->begin_adjacent_node(eTag, m_index)[node_local]);
    }
    /** Get connected edge. */
    auto get_edge(uint_t edge_local) const {
        FEATHERS_ASSERT(
            edge_local < m_mesh->num_adjacent_edges(eTag, m_index));
        return tEdgeIterBase<cMesh>(
            m_mesh, m_mesh->begin_adjacent_edge(eTag, m_index)[edge_local]);
    }
    /** Get connected face. */
    auto get_face(uint_t face_local) const {
        FEATHERS_ASSERT(
            face_local < m_mesh->num_adjacent_faces(eTag, m_index));
        return tFaceIterBase<cMesh>(
            m_mesh, m_mesh->begin_adjacent_face(eTag, m_index)[face_local]);
    }
    /** Get connected cell. */
    auto get_cell(uint_t cell_local) const {
        FEATHERS_ASSERT(
            cell_local < m_mesh->num_adjacent_cells(eTag, m_index));
        return tCellIterBase<cMesh>(
            m_mesh, m_mesh->begin_adjacent_cell(eTag, m_index)[cell_local]);
    }

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Iterate through all connected nodes. */
    template<typename tFunc>
    tFunc for_each_node(tFunc func) const {
        std::for_each(m_mesh->begin_adjacent_node(eTag, m_index),
                      m_mesh->end_adjacent_node(eTag, m_index), [&](uint_t node_index) {
            func(tNodeIterBase<cMesh>(m_mesh, node_index));
        });
        return func;
    }
    /** Iterate through all connected edges. */
    template<typename tFunc>
    tFunc for_each_edge(tFunc func) const {
        std::for_each(m_mesh->begin_adjacent_edge(eTag, m_index),
                      m_mesh->end_adjacent_edge(eTag, m_index), [&](uint_t edge_index) {
            func(tEdgeIterBase<cMesh>(m_mesh, edge_index));
        });
        return func;
    }
    /** Iterate through all connected faces. */
    /** @{ */
    template<typename tFunc>
    tFunc for_each_face(tFunc func) const {
        std::for_each(m_mesh->begin_adjacent_face(eTag, m_index),
                      m_mesh->end_adjacent_face(eTag, m_index), [&](uint_t face_index) {
            func(tFaceIterBase<cMesh>(m_mesh, face_index));
        });
        return func;
    }
    template<typename tFunc>
    tFunc for_each_face_cells(tFunc func) const {
        std::for_each(m_mesh->begin_adjacent_face(eTag, m_index),
                      m_mesh->end_adjacent_face(eTag, m_index), [&](uint_t face_index) {
            tFaceIterBase<cMesh> face(m_mesh, face_index);
            func(face.get_inner_cell(), face.get_outer_cell());
        });
        return func;
    }
    /** @} */
    /** Iterate through all connected cells. */
    template<typename tFunc>
    tFunc for_each_cell(tFunc func) const {
        std::for_each(m_mesh->begin_adjacent_cell(eTag, m_index),
                      m_mesh->end_adjacent_cell(eTag, m_index), [&](uint_t cell_index) {
            func(tCellIterBase<cMesh>(m_mesh, cell_index));
        });
        return func;
    }
};  // class tElementIterBase

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Base node iterator.
 */
template<class cMesh>
class tNodeIterBase final :
    public tElementIterBase<tNodeIterBase<cMesh>, cMesh, tNodeTag, eNodeTag> {
public:

    /** Construct base node iterator. */
    template<class cMeshPtr>
    explicit tNodeIterBase(const cMeshPtr& mesh_ptr, uint_t node_index = 0):
        tElementIterBase<tNodeIterBase<cMesh>, cMesh, tNodeTag, eNodeTag>(mesh_ptr, node_index) {
    }

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get node position. */
    const vec3_t& get_coords() const {
        return this->get_mesh()->get_node_coords(this->get_index());
    }
    /** Set node position. */
    template<class tMesh = cMesh>
    std::enable_if_t<!std::is_const_v<tMesh>> set_coords(const vec3_t& node_coords) {
        this->get_mesh()->set_node_coords(this->get_index(), node_coords);
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
template<class cMesh>
class tEdgeIterBase final :
    public tElementIterBase<tEdgeIterBase<cMesh>, cMesh, tEdgeTag, eEdgeTag> {
public:

    /** Construct base edge iterator. */
    template<class cMeshPtr>
    explicit tEdgeIterBase(const cMeshPtr& mesh_ptr, uint_t edge_index = 0):
        tElementIterBase<tEdgeIterBase<cMesh>, cMesh, tEdgeTag, eEdgeTag>(mesh_ptr, edge_index) {
    }

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get edge length. */
    real_t get_length() const {
        return this->get_mesh()->get_edge_length(this->get_index());
    }
    /** Set edge length. */
    template<class tMesh = cMesh>
    std::enable_if_t<!std::is_const_v<tMesh>> set_length(real_t edge_length) {
        this->get_mesh()->set_edge_length(this->get_index(), edge_length);
    }

    /** Get edge direction. */
    const vec3_t& get_direction() const {
        return this->get_mesh()->get_edge_direction(this->get_index());
    }
    /** Set edge direction. */
    template<class tMesh = cMesh>
    std::enable_if_t<!std::is_const_v<tMesh>> set_direction(const vec3_t& edge_direction) {
        this->get_mesh()->set_edge_direction(this->get_index(), edge_direction);
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
template<class cMesh>
class tFaceIterBase final :
    public tElementIterBase<tFaceIterBase<cMesh>, cMesh, tFaceTag, eFaceTag> {
public:

    /** Construct base face iterator. */
    template<class cMeshPtr>
    explicit tFaceIterBase(const cMeshPtr& mesh_ptr, uint_t face_index = 0):
        tElementIterBase<tFaceIterBase<cMesh>, cMesh, tFaceTag, eFaceTag>(mesh_ptr, face_index) {
    }

    /** Get connected inner cell. */
    auto get_inner_cell() const {
        // TODO: encapsulate indices!
        return this->get_cell(0);
    }
    /** Get connected outer cell. */
    auto get_outer_cell() const {
        // TODO: encapsulate indices!
        return this->get_cell(1);
    }

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get face area/length. */
    real_t get_area() const {
        return this->get_mesh()->get_face_area(this->get_index());
    }
    /** Set face area/length. */
    template<class tMesh = cMesh>
    std::enable_if_t<!std::is_const_v<tMesh>> set_area(real_t face_area) {
        this->get_mesh()->set_face_area(this->get_index(), face_area);
    }

    /** Get face normal. */
    const vec3_t& get_normal() const {
        return this->get_mesh()->get_face_normal(this->get_index());
    }
    /** Set face normal. */
    template<class tMesh = cMesh>
    std::enable_if_t<!std::is_const_v<tMesh>> set_normal(const vec3_t& face_normal) {
        this->get_mesh()->set_face_normal(this->get_index(), face_normal);
    }

    /** Get face barycenter. */
    const vec3_t& get_center_coords() const {
        return this->get_mesh()->get_face_center_coords(this->get_index());
    }
    /** Set face barycenter. */
    template<class tMesh = cMesh>
    std::enable_if_t<!std::is_const_v<tMesh>> set_center_coords(const vec3_t& face_center_coords) {
        this->get_mesh()->set_face_center_coords(this->get_index(), face_center_coords);
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
template<class cMesh>
class tCellIterBase final :
    public tElementIterBase<tCellIterBase<cMesh>, cMesh, tCellTag, eCellTag> {
public:

    /** Construct base cell iterator. */
    template<class cMeshPtr>
    explicit tCellIterBase(const cMeshPtr& mesh, uint_t cell_index = 0):
        tElementIterBase<tCellIterBase<cMesh>, cMesh, tCellTag, eCellTag>(mesh, cell_index) {
    }

    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    /** Get cell volume/area/length. */
    real_t get_volume() const {
        return this->get_mesh()->get_cell_volume(this->get_index());
    }
    /** Set cell volume/area/length. */
    template<class tMesh = cMesh>
    std::enable_if_t<!std::is_const_v<tMesh>> set_volume(real_t cell_volume) {
        this->get_mesh()->set_cell_volume(this->get_index(), cell_volume);
    }

    /** Get cell barycenter. */
    const vec3_t& get_center_coords() const {
        return this->get_mesh()->get_cell_center_coords(this->get_index());
    }
    /** Set cell barycenter. */
    template<class tMesh = cMesh>
    std::enable_if_t<!std::is_const_v<tMesh>> set_center_coords(const vec3_t& cell_center_coords) {
        this->get_mesh()->set_cell_center_coords(this->get_index(), cell_center_coords);
    }
};  // class tCellIterBase

/**
 * Mesh cells random-access iterator.
 */
/** @{ */
using tCellIter = tCellIterBase<const cMesh>;
using tCellMutableIter = tCellIterBase<cMesh>;
/** @} */

} // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

#define FACE_CELL_FUNC_ \
    ([&](tFaceIterBase<cMesh> face) { func(face.get_inner_cell(), face.get_outer_cell()); })

/** Iterator pointing to the first node with a given mark or the first mark. */
template<class cMesh>
auto begin_node(cMesh& mesh, uint_t mark = 0) {
    return tNodeIterBase<cMesh>(&mesh, mesh.begin_node(mark));
}
/** Iterator pointing to the first edge with a given mark or the first mark. */
template<class cMesh>
auto begin_edge(cMesh& mesh, uint_t mark = 0) {
    return tEdgeIterBase<cMesh>(&mesh, mesh.begin_edge(mark));
}
/** Iterator pointing to the first face with a given mark or the first mark. */
template<class cMesh>
auto begin_face(cMesh& mesh, uint_t mark = 0) {
    return tFaceIterBase<cMesh>(&mesh, mesh.begin_face(mark));
}
/** Iterator pointing to the first cell with a given mark or the first mark. */
template<class cMesh>
auto begin_cell(cMesh& mesh, uint_t mark = 0) {
    return tCellIterBase<cMesh>(&mesh, mesh.begin_cell(mark));
}

/** Iterator pointing to a node after the last node with a given mark or the last mark. */
/** @{ */
template<class cMesh>
auto end_node(cMesh& mesh) {
    const uint_t mark = mesh.num_node_marks() - 1;
    return tNodeIterBase<cMesh>(&mesh, mesh.end_node(mark));
}
template<class cMesh>
auto end_node(cMesh& mesh, uint_t mark) {
    return tNodeIterBase<cMesh>(&mesh, mesh.end_node(mark));
}
/** @} */
/** Iterator pointing to an edge after the last edge with a given mark or the last mark. */
/** @{ */
template<class cMesh>
auto end_edge(cMesh& mesh) {
    const uint_t mark = mesh.num_edge_marks() - 1;
    return tEdgeIterBase<cMesh>(&mesh, mesh.end_edge(mark));
}
template<class cMesh>
auto end_edge(cMesh& mesh, uint_t mark) {
    return tEdgeIterBase<cMesh>(&mesh, mesh.end_edge(mark));
}
/** @} */
/** Iterator pointing to a face after the last face with a given mark or the last mark. */
/** @{ */
template<class cMesh>
auto end_face(cMesh& mesh) {
    const uint_t mark = mesh.num_face_marks() - 1;
    return tFaceIterBase<cMesh>(&mesh, mesh.end_face(mark));
}
template<class cMesh>
auto end_face(cMesh& mesh, uint_t mark) {
    return tFaceIterBase<cMesh>(&mesh, mesh.end_face(mark));
}
/** @} */
/** Iterator pointing to a cell after the last cell with a given mark or the last mark. */
/** @{ */
template<class cMesh>
auto end_cell(cMesh& mesh) {
    const uint_t mark = mesh.num_cell_marks() - 1;
    return tCellIterBase<cMesh>(&mesh, mesh.end_cell(mark));
}
template<class cMesh>
auto end_cell(cMesh& mesh, uint_t mark) {
    return tCellIterBase<cMesh>(&mesh, mesh.end_cell(mark));
}
/** @} */

/** Iterate through all nodes with a given mark or all marks. */
/** @{ */
template<class cMesh, typename tFunc>
tFunc for_each_node(cMesh& mesh, tFunc func) {
    return for_range(begin_node(mesh), end_node(mesh), func);
}
template<class cMesh, typename tFunc>
tFunc for_each_node(cMesh& mesh, uint_t mark, tFunc func) {
    return for_range(begin_node(mesh, mark), end_node(mesh, mark), func);
}
/** @} */
/** Iterate through all edges with a given mark or all marks. */
/** @{ */
template<class cMesh, typename tFunc>
tFunc for_each_edge(cMesh& mesh, tFunc func) {
    return for_range(begin_edge(mesh), end_edge(mesh), func);
}
template<class cMesh, typename tFunc>
tFunc for_each_edge(cMesh& mesh, uint_t mark, tFunc func) {
    return for_range(begin_edge(mesh, mark), end_edge(mesh, mark), func);
}
/** @} */
/** Iterate through all faces with a given mark or all marks. */
/** @{ */
template<class cMesh, typename tFunc>
tFunc for_each_face(cMesh& mesh, uint_t mark, tFunc func) {
    return for_range(begin_face(mesh, mark), end_face(mesh, mark), func);
}
template<class cMesh, typename tFunc>
tFunc for_each_face(cMesh& mesh, tFunc func) {
    return for_range(begin_face(mesh), end_face(mesh), func);
}
template<class cMesh, typename tFunc>
tFunc for_each_face_cells(cMesh& mesh, tFunc func) {
    return for_range(begin_face(mesh), end_face(mesh), FACE_CELL_FUNC_), func;
}
template<class cMesh, typename tFunc>
tFunc for_each_face_cells(cMesh& mesh, uint_t mark, tFunc func) {
    return for_range(begin_face(mesh, mark), end_face(mesh, mark), FACE_CELL_FUNC_), func;
}
/** @} */
/** Iterate through all cells with a given mark or all marks. */
/** @{ */
template<class cMesh, typename tFunc>
tFunc for_each_cell(cMesh& mesh, tFunc func) {
    return for_range(begin_cell(mesh), end_cell(mesh), func);
}
template<class cMesh, typename tFunc>
tFunc for_each_cell(cMesh& mesh, uint_t mark, tFunc func) {
    return for_range(begin_cell(mesh, mark), end_cell(mesh, mark), func);
}
/** @} */

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/** Iterator pointing to the first interior node. */
template<class cMesh>
auto begin_interior_node(cMesh& mesh) {
    return begin_node(mesh, 0);
}
/** Iterator pointing to the first interior edge. */
template<class cMesh>
auto begin_interior_edge(cMesh& mesh) {
    return begin_edge(mesh, 0);
}
/** Iterator pointing to the first interior face. */
template<class cMesh>
auto begin_interior_face(cMesh& mesh) {
    return begin_face(mesh, 0);
}
/** Iterator pointing to the first interior cell. */
template<class cMesh>
auto begin_interior_cell(cMesh& mesh) {
    return begin_cell(mesh, 0);
}

/** Iterator pointing to a node after the last node. */
template<class cMesh>
auto end_interior_node(cMesh& mesh) {
    return end_node(mesh, 0);
}
/** Iterator pointing to an edge after the last node. */
template<class cMesh>
auto end_interior_edge(cMesh& mesh) {
    return end_edge(mesh, 0);
}
/** Iterator pointing to a face after the last node. */
template<class cMesh>
auto end_interior_face(cMesh& mesh) {
    return end_face(mesh, 0);
}
/** Iterator pointing to a cell after the last node. */
template<class cMesh>
auto end_interior_cell(cMesh& mesh) {
    return end_cell(mesh, 0);
}

/** Iterate through all interior nodes. */
template<class cMesh, typename tFunc>
tFunc for_each_interior_node(cMesh& mesh, tFunc func) {
    return for_range(begin_interior_node(mesh), end_interior_node(mesh), func);
}
/** Iterate through all interior edges. */
template<class cMesh, typename tFunc>
tFunc for_each_interior_edge(cMesh& mesh, tFunc func) {
    return for_range(begin_interior_edge(mesh), end_interior_edge(mesh), func);
}
/** Iterate through all interior faces. */
/** @{ */
template<class cMesh, typename tFunc>
tFunc for_each_interior_face(cMesh& mesh, tFunc func) {
    return for_range(begin_interior_face(mesh), end_interior_face(mesh), func);
}
template<class cMesh, typename tFunc>
tFunc for_each_interior_face_cells(cMesh& mesh, tFunc func) {
    return for_range(begin_interior_face(mesh), end_interior_face(mesh), FACE_CELL_FUNC_), func;
}
/** @} */
/** Iterate through all interior cells. */
template<class cMesh, typename tFunc>
tFunc for_each_interior_cell(cMesh& mesh, tFunc func) {
    return for_range(begin_interior_cell(mesh), end_interior_cell(mesh), func);
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/** Iterator pointing to the boundary first node with a given mark or the first boundary mark. */
template<class cMesh>
auto begin_boundary_node(cMesh& mesh, uint_t mark = 1) {
    FEATHERS_ASSERT(mark >= 1);
    return begin_node(mesh, mark);
}
/** Iterator pointing to the boundary first edge with a given mark or the first boundary mark. */
template<class cMesh>
auto begin_boundary_edge(cMesh& mesh, uint_t mark = 1) {
    FEATHERS_ASSERT(mark >= 1);
    return begin_edge(mesh, mark);
}
/** Iterator pointing to the boundary first face with a given mark or the first boundary mark. */
template<class cMesh>
auto begin_boundary_face(cMesh& mesh, uint_t mark = 1) {
    FEATHERS_ASSERT(mark >= 1);
    return begin_face(mesh, mark);
}
/** Iterator pointing to the boundary first cell with a given mark or the first boundary mark. */
template<class cMesh>
auto begin_boundary_cell(cMesh& mesh, uint_t mark = 1) {
    FEATHERS_ASSERT(mark >= 1);
    return begin_cell(mesh, mark);
}

/** Iterator pointing to a node after the last node with a given or the last boundary mark. */
/** @{ */
template<class cMesh>
auto end_boundary_node(cMesh& mesh) {
    return end_node(mesh);
}
template<class cMesh>
auto end_boundary_node(cMesh& mesh, uint_t mark) {
    FEATHERS_ASSERT(mark >= 1);
    return end_node(mesh, mark);
}
/** @} */
/** Iterator pointing to an edge after the last edge with a given or the last boundary mark. */
/** @{ */
template<class cMesh>
auto end_boundary_edge(cMesh& mesh) {
    return end_edge(mesh);
}
template<class cMesh>
auto end_boundary_edge(cMesh& mesh, uint_t mark) {
    FEATHERS_ASSERT(mark >= 1);
    return end_edge(mesh, mark);
}
/** @} */
/** Iterator pointing to a face after the last face with a given or the last boundary mark. */
/** @{ */
template<class cMesh>
auto end_boundary_face(cMesh& mesh) {
    return end_face(mesh);
}
template<class cMesh>
auto end_boundary_face(cMesh& mesh, uint_t mark) {
    FEATHERS_ASSERT(mark >= 1);
    return end_face(mesh, mark);
}
/** @} */
/** Iterator pointing to a cell after the last cell with a given or the last boundary mark. */
/** @{ */
template<class cMesh>
auto end_boundary_cell(cMesh& mesh) {
    return end_cell(mesh);
}
template<class cMesh>
auto end_boundary_cell(cMesh& mesh, uint_t mark) {
    FEATHERS_ASSERT(mark >= 1);
    return end_cell(mesh, mark);
}
/** @} */

/** Iterate through all boundary nodes with a given mark or all boundary marks. */
/** @{ */
template<class cMesh, typename tFunc>
tFunc for_each_boundary_node(cMesh& mesh, tFunc func) {
    return for_range(begin_boundary_node(mesh), end_boundary_node(mesh), func);
}
template<class cMesh, typename tFunc>
tFunc for_each_boundary_node(cMesh& mesh, uint_t mark, tFunc func) {
    return for_range(begin_boundary_node(mesh, mark), end_boundary_node(mesh, mark), func);
}
/** @} */
/** Iterate through all boundary edges with a given mark or all boundary marks. */
/** @{ */
template<class cMesh, typename tFunc>
tFunc for_each_boundary_edge(cMesh& mesh, tFunc func) {
    return for_range(begin_boundary_edge(mesh), end_boundary_edge(mesh), func);
}
template<class cMesh, typename tFunc>
tFunc for_each_boundary_edge(cMesh& mesh, uint_t mark, tFunc func) {
    return for_range(begin_boundary_edge(mesh, mark), end_boundary_edge(mesh, mark), func);
}
/** @} */
/** Iterate through all boundary faces with a given mark or all boundary marks. */
/** @{ */
template<class cMesh, typename tFunc>
tFunc for_each_boundary_face(cMesh& mesh, tFunc func) {
    return for_range(begin_boundary_face(mesh), end_boundary_face(mesh), func);
}
template<class cMesh, typename tFunc>
tFunc for_each_boundary_face(cMesh& mesh, uint_t mark, tFunc func) {
    return for_range(begin_boundary_face(mesh, mark), end_boundary_face(mesh, mark), func);
}
template<class cMesh, typename tFunc>
tFunc for_each_boundary_face_cells(cMesh& mesh, tFunc func) {
    return for_range(begin_boundary_face(mesh), end_boundary_face(mesh), FACE_CELL_FUNC_), func;
}
template<class cMesh, typename tFunc>
tFunc for_each_boundary_face_cells(cMesh& mesh, uint_t mark, tFunc func) {
    return for_range(begin_boundary_face(mesh, mark), end_boundary_face(mesh, mark), FACE_CELL_FUNC_), func;
}
/** @} */
/** Iterate through all boundary cells with a given mark or all boundary marks. */
/** @{ */
template<class cMesh, typename tFunc>
tFunc for_each_boundary_cell(cMesh& mesh, tFunc func) {
    return for_range(begin_boundary_cell(mesh), end_boundary_cell(mesh), func);
}
template<class cMesh, typename tFunc>
tFunc for_each_boundary_cell(cMesh& mesh, uint_t mark, tFunc func) {
    return for_range(begin_boundary_cell(mesh, mark), end_boundary_cell(mesh, mark), func);
}
/** @} */

#undef FACE_CELL_FUNC_

} // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif // MESH_ITERATORS_HH_
