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
#ifndef SKUNK_MESH_ITER_HH
#define SKUNK_MESH_ITER_HH

#include <SkunkBase.hh>
#include "Mesh.hh"
#include <libSkunkMisc/SkunkMiscParallel.hh>

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

template<typename mesh_t>
class mesh_node_iter_struct_t;
template<typename mesh_t>
class mesh_edge_iter_struct_t;
template<typename mesh_t>
class mesh_face_iter_struct_t;
template<typename mesh_t>
class mesh_cell_iter_struct_t;

/** Base element iterator.  */
template<typename iter_t, typename mesh_t, typename base_elem_t>
class mesh_elem_iter_struct_t {
private:
    using elem_t = std::conditional_t<std::is_const<mesh_t>::value, std::add_const_t<base_elem_t>, base_elem_t>;
    mesh_t* m_mesh;
    elem_t* m_elem;
    uint_t  m_elem_ind;

    /**************************************************************************/
    /**************************************************************************/

protected:

    /** @internal
     ** Construct base element iterator. */
    template<typename mesh_ptr_t>
    SKUNK_INLINE mesh_elem_iter_struct_t(const mesh_ptr_t& mesh, uint_t elem_ind)
        : m_mesh(&*mesh),
          m_elem(iter_t::get_elem_(m_mesh, elem_ind)), m_elem_ind(elem_ind) {
    }

    /**************************************************************************/
    /**************************************************************************/

public:

    /** Get associated mesh. */
    SKUNK_INLINE mesh_t* get_mesh() const {
        return m_mesh;
    }

    /** Cast to index operator. */
    SKUNK_INLINE operator uint_t() const {
        return m_elem_ind;
    }
    /** Dereference operator. */
    SKUNK_INLINE elem_t& operator*() const {
        return *m_elem;
    }
    /** Dereference operator. */
    SKUNK_INLINE elem_t* operator->() const {
        return m_elem;
    }

    /**************************************************************************/
    /**************************************************************************/

    /** Difference operator. */
    SKUNK_INLINE int_t operator-(const iter_t& other) const {
        SKUNK_ASSERT(m_mesh == other.m_mesh);
        return m_elem_ind - other.m_elem_ind;
    }

    /** Equality operator. */
    SKUNK_INLINE bool operator==(const iter_t& other) const {
        SKUNK_ASSERT(m_mesh == other.m_mesh);
        return m_elem_ind == other.m_elem_ind;
    }
    /** Inequality operator. */
    SKUNK_INLINE bool operator!=(const iter_t& other) const {
        SKUNK_ASSERT(m_mesh == other.m_mesh);
        return m_elem_ind != other.m_elem_ind;
    }

    /** Lexicographical less than operator. */
    SKUNK_INLINE bool operator<(const iter_t& other) const {
        SKUNK_ASSERT(m_mesh == other.m_mesh);
        return m_elem_ind < other.m_elem_ind;
    }
    /** Lexicographical less than or equal operator. */
    SKUNK_INLINE bool operator<=(const iter_t& other) const {
        SKUNK_ASSERT(m_mesh == other.m_mesh);
        return m_elem_ind <= other.m_elem_ind;
    }

    /** Lexicographical greater than operator. */
    SKUNK_INLINE bool operator>(const iter_t& other) const {
        SKUNK_ASSERT(m_mesh == other.m_mesh);
        return m_elem_ind > other.m_elem_ind;
    }
    /** Lexicographical greater than or equal operator. */
    SKUNK_INLINE bool operator>=(const iter_t& other) const {
        SKUNK_ASSERT(m_mesh == other.m_mesh);
        return m_elem_ind >= other.m_elem_ind;
    }

    /**************************************************************************/
    /**************************************************************************/

    /** Increment operator. */
    /** @{ */
    SKUNK_INLINE iter_t& operator++() {
        m_elem = iter_t::get_elem_(m_mesh, ++m_elem_ind);
        return static_cast<iter_t&>(*this);
    }
    SKUNK_INLINE const iter_t operator++(int_t) {
        iter_t iter(*this);
        return ++*this, iter;
    }
    /** @} */

    /** Decrement operator. */
    /** @{ */
    SKUNK_INLINE iter_t& operator--() {
        m_elem = iter_t::get_elem_(m_mesh, --m_elem_ind);
        return static_cast<iter_t&>(*this);
    }
    SKUNK_INLINE const iter_t operator--(int_t) {
        iter_t iter(*this);
        return --*this, iter;
    }
    /** @} */

    /** Addition operator. */
    /** @{ */
    SKUNK_INLINE iter_t& operator+=(int_t offset) {
        m_elem = iter_t::get_elem_(m_mesh, m_elem_ind += offset);
        return static_cast<iter_t&>(*this);
    }
    SKUNK_INLINE iter_t operator+(int_t offset) const {
        return iter_t(*this) += offset;
    }
    /** @} */

    /** Subtraction operator. */
    /** @{ */
    SKUNK_INLINE iter_t& operator-=(int_t offset) {
        m_elem = iter_t::get_elem_(m_mesh, m_elem_ind -= offset);
        return static_cast<iter_t&>(*this);
    }
    SKUNK_INLINE iter_t operator-(int_t offset) const {
        return iter_t(*this) -= offset;
    }
    /** @} */

    /**************************************************************************/
    /**************************************************************************/

    /** Get connected node. */
    SKUNK_INLINE auto get_node(uint_t node_loc) const {
        SKUNK_ASSERT(node_loc < m_elem->num_nodes());
        return mesh_node_iter_struct_t<mesh_t>(m_mesh, m_elem->begin_node()[node_loc]);
    }
    /** Get connected edge. */
    SKUNK_INLINE auto get_edge(uint_t edge_loc) const {
        SKUNK_ASSERT(edge_loc < m_elem->num_edges());
        return mesh_edge_iter_struct_t<mesh_t>(m_mesh, m_elem->begin_edge()[edge_loc]);
    }
    /** Get connected face. */
    SKUNK_INLINE auto get_face(uint_t face_loc) const {
        SKUNK_ASSERT(face_loc < m_elem->num_faces());
        return mesh_face_iter_struct_t<mesh_t>(m_mesh, m_elem->begin_face()[face_loc]);
    }
    /** Get connected cell. */
    SKUNK_INLINE auto get_cell(uint_t cell_loc) const {
        SKUNK_ASSERT(cell_loc < m_elem->num_cells());
        return mesh_cell_iter_struct_t<mesh_t>(m_mesh, m_elem->begin_cell()[cell_loc]);
    }

    /**************************************************************************/
    /**************************************************************************/

    /** Iterate through all connected nodes. */
    template<typename func_t>
    SKUNK_INLINE func_t for_each_node(func_t func) const {
        elem_t& elem = **this;
        std::for_each(elem.begin_node(), elem.end_node(), [&](int_t node_ind) {
            func(mesh_node_iter_struct_t<mesh_t>(m_mesh, node_ind));
        });
        return func;
    }
    /** Iterate through all connected edges. */
    template<typename func_t>
    SKUNK_INLINE func_t for_each_edge(func_t func) const {
        elem_t& elem = **this;
        std::for_each(elem.begin_edge(), elem.end_edge(), [&](int_t edge_ind) {
            func(mesh_edge_iter_struct_t<mesh_t>(m_mesh, edge_ind));
        });
        return func;
    }
    /** Iterate through all connected faces. */
    /** @{ */
    template<typename func_t>
    SKUNK_INLINE func_t for_each_face(func_t func) const {
        elem_t& elem = **this;
        std::for_each(elem.begin_face(), elem.end_face(), [&](int_t face_ind) {
            func(mesh_face_iter_struct_t<mesh_t>(m_mesh, face_ind));
        });
        return func;
    }
    template<typename func_t>
    SKUNK_INLINE func_t for_each_face_cells(func_t func) const {
        elem_t& elem = **this;
        std::for_each(elem.begin_face(), elem.end_face(), [&](int_t face_ind) {
            mesh_face_iter_struct_t<mesh_t> face(m_mesh, face_ind);
            func(face.get_inner_cell(), face.get_outer_cell());
        });
        return func;
    }
    /** @} */
    /** Iterate through all connected cells. */
    template<typename func_t>
    SKUNK_INLINE func_t for_each_cell(func_t func) const {
        elem_t& elem = **this;
        std::for_each(elem.begin_cell(), elem.end_cell(), [&](int_t cell_ind) {
            func(mesh_cell_iter_struct_t<mesh_t>(m_mesh, cell_ind));
        });
        return func;
    }
};  // class mesh_elem_iter_struct_t

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/** Mesh nodes random-access iterator. */
/** @{ */
using NodeIter = mesh_node_iter_struct_t<UMesh>;
using NodeConstIter = mesh_node_iter_struct_t<const UMesh>;
/** @} */

/** Mesh edges random-access iterator. */
/** @{ */
using EdgeIter = mesh_edge_iter_struct_t<UMesh>;
using EdgeConstIter = mesh_edge_iter_struct_t<const UMesh>;
/** @} */

/** Mesh faces random-access iterator. */
/** @{ */
using FaceIter = mesh_face_iter_struct_t<UMesh>;
using FaceConstIter = mesh_face_iter_struct_t<const UMesh>;
/** @} */

/** Mesh cells random-access iterator. */
/** @{ */
using CellIter = mesh_cell_iter_struct_t<UMesh>;
using CellConstIter = mesh_cell_iter_struct_t<const UMesh>;
/** @} */

/**************************************************************************/
/**************************************************************************/

/** Iterator pointing to the first node with a given mark or the first mark. */
template<typename mesh_t>
SKUNK_INLINE auto begin_node(mesh_t& mesh, uint_t mark = 0) {
    return mesh_node_iter_struct_t<mesh_t>(&mesh, mesh.begin_node(mark));
}   // begin_node
/** Iterator pointing to the first edge with a given mark or the first mark. */
template<typename mesh_t>
SKUNK_INLINE auto begin_edge(mesh_t& mesh, uint_t mark = 0) {
    return mesh_edge_iter_struct_t<mesh_t>(&mesh, mesh.begin_edge(mark));
}   // begin_edge
/** Iterator pointing to the first face with a given mark or the first mark. */
template<typename mesh_t>
SKUNK_INLINE auto begin_face(mesh_t& mesh, uint_t mark = 0) {
    return mesh_face_iter_struct_t<mesh_t>(&mesh, mesh.begin_face(mark));
}   // begin_face
/** Iterator pointing to the first cell with a given mark or the first mark. */
template<typename mesh_t>
SKUNK_INLINE auto begin_cell(mesh_t& mesh, uint_t mark = 0) {
    return mesh_cell_iter_struct_t<mesh_t>(&mesh, mesh.begin_cell(mark));
}   // begin_cell

/** Iterator pointing to a node after the last node with a given mark or the last mark. */
/** @{ */
template<typename mesh_t>
SKUNK_INLINE auto end_node(mesh_t& mesh, uint_t mark) {
    return mesh_node_iter_struct_t<mesh_t>(&mesh, mesh.end_node(mark));
}   // end_node
template<typename mesh_t>
SKUNK_INLINE auto end_node(mesh_t& mesh) {
    const uint_t mark = mesh.num_node_marks() - 1;
    return mesh_node_iter_struct_t<mesh_t>(&mesh, mesh.end_node(mark));
}   // end_node
/** @} */
/** Iterator pointing to an edge after the last edge with a given mark or the last mark. */
/** @{ */
template<typename mesh_t>
SKUNK_INLINE auto end_edge(mesh_t& mesh, uint_t mark) {
    return mesh_edge_iter_struct_t<mesh_t>(&mesh, mesh.end_edge(mark));
}   // end_edge
template<typename mesh_t>
SKUNK_INLINE auto end_edge(mesh_t& mesh) {
    const uint_t mark = mesh.num_edge_marks() - 1;
    return mesh_edge_iter_struct_t<mesh_t>(&mesh, mesh.end_edge(mark));
}   // end_edge
/** @} */
/** Iterator pointing to a face after the last face with a given mark or the last mark. */
/** @{ */
template<typename mesh_t>
SKUNK_INLINE auto end_face(mesh_t& mesh, uint_t mark) {
    return mesh_face_iter_struct_t<mesh_t>(&mesh, mesh.end_face(mark));
}   // end_face
template<typename mesh_t>
SKUNK_INLINE auto end_face(mesh_t& mesh) {
    const uint_t mark = mesh.num_face_marks() - 1;
    return mesh_face_iter_struct_t<mesh_t>(&mesh, mesh.end_face(mark));
}   // end_face
/** @} */
/** Iterator pointing to a cell after the last cell with a given mark or the last mark. */
/** @{ */
template<typename mesh_t>
SKUNK_INLINE auto end_cell(mesh_t& mesh, uint_t mark) {
    return mesh_cell_iter_struct_t<mesh_t>(&mesh, mesh.end_cell(mark));
}   // end_cell
template<typename mesh_t>
SKUNK_INLINE auto end_cell(mesh_t& mesh) {
    const uint_t mark = mesh.num_cell_marks() - 1;
    return mesh_cell_iter_struct_t<mesh_t>(&mesh, mesh.end_cell(mark));
}   // end_cell
/** @} */

/** Iterate through all nodes with a given mark or all marks. */
/** @{ */
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_node(mesh_t& mesh, uint_t mark, func_t func) {
    return for_range(begin_node(mesh, mark), end_node(mesh, mark), func);
}   // for_each_node
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_node(mesh_t& mesh, func_t func) {
    return for_range(begin_node(mesh), end_node(mesh), func);
}   // for_each_node
/** @} */
/** Iterate through all edges with a given mark or all marks. */
/** @{ */
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_edge(mesh_t& mesh, uint_t mark, func_t func) {
    return for_range(begin_edge(mesh, mark), end_edge(mesh, mark), func);
}   // for_each_edge
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_edge(mesh_t& mesh, func_t func) {
    return for_range(begin_edge(mesh), end_edge(mesh), func);
}   // for_each_edge
/** @} */
/** Iterate through all faces with a given mark or all marks. */
/** @{ */
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_face(mesh_t& mesh, uint_t mark, func_t func) {
    return for_range(begin_face(mesh, mark), end_face(mesh, mark), func);
}   // for_each_face
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_face(mesh_t& mesh, func_t func) {
    return for_range(begin_face(mesh), end_face(mesh), func);
}   // for_each_face
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_face_cells(mesh_t& mesh, uint_t mark, func_t func) {
    for_range(begin_face(mesh, mark),
                end_face(mesh, mark), [&](mesh_face_iter_struct_t<mesh_t> face) {
        func(face.get_inner_cell(), face.get_outer_cell());
    });
    return func;
}   // for_each_face
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_face_cells(mesh_t& mesh, func_t func) {
    for_range(begin_face(mesh),
                end_face(mesh), [&](mesh_face_iter_struct_t<mesh_t> face) {
        func(face.get_inner_cell(), face.get_outer_cell());
    });
    return func;
}   // for_each_face
/** @} */
/** Iterate through all cells with a given mark or all marks. */
/** @{ */
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_cell(mesh_t& mesh, uint_t mark, func_t func) {
    return for_range(begin_cell(mesh, mark), end_cell(mesh, mark), func);
}   // for_each_cell
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_cell(mesh_t& mesh, func_t func) {
    return for_range(begin_cell(mesh), end_cell(mesh), func);
}   // for_each_cell
/** @} */

/**************************************************************************/
/**************************************************************************/

/** Iterator pointing to the first interior node. */
template<typename mesh_t>
SKUNK_INLINE auto begin_interior_node(mesh_t& mesh) {
    return begin_node(mesh, 0);
}   // begin_interior_node
/** Iterator pointing to the first interior edge. */
template<typename mesh_t>
SKUNK_INLINE auto begin_interior_edge(mesh_t& mesh) {
    return begin_edge(mesh, 0);
}   // begin_interior_edge
/** Iterator pointing to the first interior face. */
template<typename mesh_t>
SKUNK_INLINE auto begin_interior_face(mesh_t& mesh) {
    return begin_face(mesh, 0);
}   // begin_interior_face
/** Iterator pointing to the first interior cell. */
template<typename mesh_t>
SKUNK_INLINE auto begin_interior_cell(mesh_t& mesh) {
    return begin_cell(mesh, 0);
}   // begin_interior_cell

/** Iterator pointing to a node after the last node. */
template<typename mesh_t>
SKUNK_INLINE auto end_interior_node(mesh_t& mesh) {
    return end_node(mesh, 0);
}   // end_interior_node
/** Iterator pointing to an edge after the last node. */
template<typename mesh_t>
SKUNK_INLINE auto end_interior_edge(mesh_t& mesh) {
    return end_edge(mesh, 0);
}   // end_interior_edge
/** Iterator pointing to a face after the last node. */
template<typename mesh_t>
SKUNK_INLINE auto end_interior_face(mesh_t& mesh) {
    return end_face(mesh, 0);
}   // end_interior_face
/** Iterator pointing to a cell after the last node. */
template<typename mesh_t>
SKUNK_INLINE auto end_interior_cell(mesh_t& mesh) {
    return end_cell(mesh, 0);
}   // end_interior_cell

/** Iterate through all interior nodes. */
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_interior_node(mesh_t& mesh, func_t func) {
    return for_range(begin_interior_node(mesh), end_interior_node(mesh), func);
}   // for_each_interior_node
/** Iterate through all interior edges. */
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_interior_edge(mesh_t& mesh, func_t func) {
    return for_range(begin_interior_edge(mesh), end_interior_edge(mesh), func);
}   // for_each_interior_edge
/** Iterate through all interior faces. */
/** @{ */
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_interior_face(mesh_t& mesh, func_t func) {
    return for_range(begin_interior_face(mesh), end_interior_face(mesh), func);
}   // for_each_interior_face
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_interior_face_cells(mesh_t& mesh, func_t func) {
    for_range(begin_interior_face(mesh),
                end_interior_face(mesh), [&](mesh_face_iter_struct_t<mesh_t> face) {
        func(face.get_inner_cell(), face.get_outer_cell());
    });
    return func;
}   // for_each_interior_face_cells
/** @} */
/** Iterate through all interior cells. */
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_interior_cell(mesh_t& mesh, func_t func) {
    return for_range(begin_interior_cell(mesh), end_interior_cell(mesh), func);
}   // for_each_interior_cell

/**************************************************************************/
/**************************************************************************/

/** Iterator pointing to the boundary first node with a given mark or the first boundary mark. */
template<typename mesh_t>
SKUNK_INLINE auto begin_boundary_node(mesh_t& mesh, uint_t mark = 1) {
    SKUNK_ASSERT(1 <= mark);
    return begin_node(mesh, mark);
}   // begin_boundary_node
/** Iterator pointing to the boundary first edge with a given mark or the first boundary mark. */
template<typename mesh_t>
SKUNK_INLINE auto begin_boundary_edge(mesh_t& mesh, uint_t mark = 1) {
    SKUNK_ASSERT(1 <= mark);
    return begin_edge(mesh, mark);
}   // begin_boundary_edge
/** Iterator pointing to the boundary first face with a given mark or the first boundary mark. */
template<typename mesh_t>
SKUNK_INLINE auto begin_boundary_face(mesh_t& mesh, uint_t mark = 1) {
    SKUNK_ASSERT(1 <= mark);
    return begin_face(mesh, mark);
}   // begin_boundary_face
/** Iterator pointing to the boundary first cell with a given mark or the first boundary mark. */
template<typename mesh_t>
SKUNK_INLINE auto begin_boundary_cell(mesh_t& mesh, uint_t mark = 1) {
    SKUNK_ASSERT(1 <= mark);
    return begin_cell(mesh, mark);
}   // begin_boundary_cell

/** Iterator pointing to a node after the last node with a given or the last boundary mark. */
/** @{ */
template<typename mesh_t>
SKUNK_INLINE auto end_boundary_node(mesh_t& mesh, uint_t mark) {
    SKUNK_ASSERT(1 <= mark);
    return end_node(mesh, mark);
}   // end_boundary_node
template<typename mesh_t>
SKUNK_INLINE auto end_boundary_node(mesh_t& mesh) {
    return end_node(mesh);
}   // end_boundary_node
/** @} */
/** Iterator pointing to an edge after the last edge with a given or the last boundary mark. */
/** @{ */
template<typename mesh_t>
SKUNK_INLINE auto end_boundary_edge(mesh_t& mesh, uint_t mark) {
    SKUNK_ASSERT(1 <= mark);
    return end_edge(mesh, mark);
}   // end_boundary_edge
template<typename mesh_t>
SKUNK_INLINE auto end_boundary_edge(mesh_t& mesh) {
    return end_edge(mesh);
}   // end_boundary_edge
/** @} */
/** Iterator pointing to a face after the last face with a given or the last boundary mark. */
/** @{ */
template<typename mesh_t>
SKUNK_INLINE auto end_boundary_face(mesh_t& mesh, uint_t mark) {
    SKUNK_ASSERT(1 <= mark);
    return end_face(mesh, mark);
}   // end_boundary_face
template<typename mesh_t>
SKUNK_INLINE auto end_boundary_face(mesh_t& mesh) {
    return end_face(mesh);
}   // end_boundary_face
/** @} */
/** Iterator pointing to a cell after the last cell with a given or the last boundary mark. */
/** @{ */
template<typename mesh_t>
SKUNK_INLINE auto end_boundary_cell(mesh_t& mesh, uint_t mark) {
    SKUNK_ASSERT(1 <= mark);
    return end_cell(mesh, mark);
}   // end_boundary_cell
template<typename mesh_t>
SKUNK_INLINE auto end_boundary_cell(mesh_t& mesh) {
    return end_cell(mesh);
}   // end_boundary_cell
/** @} */

/** Iterate through all boundary nodes with a given mark or all boundary marks. */
/** @{ */
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_boundary_node(mesh_t& mesh, uint_t mark, func_t func) {
    return for_range(begin_boundary_node(mesh, mark), end_boundary_node(mesh, mark), func);
}   // for_each_boundary_node
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_boundary_node(mesh_t& mesh, func_t func) {
    return for_range(begin_boundary_node(mesh), end_boundary_node(mesh), func);
}   // for_each_boundary_node
/** @} */
/** Iterate through all boundary edges with a given mark or all boundary marks. */
/** @{ */
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_boundary_edge(mesh_t& mesh, uint_t mark, func_t func) {
    return for_range(begin_boundary_edge(mesh, mark), end_boundary_edge(mesh, mark), func);
}   // for_each_boundary_edge
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_boundary_edge(mesh_t& mesh, func_t func) {
    return for_range(begin_boundary_edge(mesh), end_boundary_edge(mesh), func);
}   // for_each_boundary_edge
/** @} */
/** Iterate through all boundary faces with a given mark or all boundary marks. */
/** @{ */
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_boundary_face(mesh_t& mesh, uint_t mark, func_t func) {
    return for_range(begin_boundary_face(mesh, mark), end_boundary_face(mesh, mark), func);
}   // for_each_boundary_face
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_boundary_face(mesh_t& mesh, func_t func) {
    return for_range(begin_boundary_face(mesh), end_boundary_face(mesh), func);
}   // for_each_boundary_face
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_boundary_face_cells(mesh_t& mesh, uint_t mark, func_t func) {
    for_range(begin_boundary_face(mesh, mark),
                end_boundary_face(mesh, mark), [&](mesh_face_iter_struct_t<mesh_t> face) {
        func(face.get_inner_cell(), face.get_outer_cell());
    });
    return func;
}   // for_each_boundary_face_cells
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_boundary_face_cells(mesh_t& mesh, func_t func) {
    for_range(begin_boundary_face(mesh),
                end_boundary_face(mesh), [&](mesh_face_iter_struct_t<mesh_t> face) {
        func(face.get_inner_cell(), face.get_outer_cell());
    });
    return func;
}   // for_each_boundary_face_cells
/** @} */
/** Iterate through all boundary cells with a given mark or all boundary marks. */
/** @{ */
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_boundary_cell(mesh_t& mesh, uint_t mark, func_t func) {
    return for_range(begin_boundary_cell(mesh, mark), end_boundary_cell(mesh, mark), func);
}   // for_each_boundary_cell
template<typename mesh_t, typename func_t>
SKUNK_INLINE func_t for_each_boundary_cell(mesh_t& mesh, func_t func) {
    return for_range(begin_boundary_cell(mesh), end_boundary_cell(mesh), func);
}   // for_each_boundary_cell
/** @} */

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/** Base node iterator. */
template<typename mesh_t>
class mesh_node_iter_struct_t final :
    public mesh_elem_iter_struct_t<mesh_node_iter_struct_t<mesh_t>, mesh_t, mesh_node_t> {
private:
    friend class mesh_elem_iter_struct_t<mesh_node_iter_struct_t<mesh_t>, mesh_t, mesh_node_t>;

    /** @internal
     ** Get mesh element at index. */
    SKUNK_INLINE static auto get_elem_(mesh_t* mesh, uint_t node_ind) {
        SKUNK_ASSERT(mesh != nullptr);
        return node_ind < mesh->num_nodes() ? &mesh->get_node(node_ind) : nullptr;
    }

public:

    /** Construct base node iterator. */
    template<typename mesh_ptr_t>
    SKUNK_INLINE explicit mesh_node_iter_struct_t(const mesh_ptr_t& mesh, uint_t node_ind = 0)
        : mesh_elem_iter_struct_t<mesh_node_iter_struct_t<mesh_t>, mesh_t, mesh_node_t>(mesh, node_ind) {
    }
};  // class NodeIter

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/** Base edge iterator. */
template<typename mesh_t>
class mesh_edge_iter_struct_t final :
    public mesh_elem_iter_struct_t<mesh_edge_iter_struct_t<mesh_t>, mesh_t, Edge> {
private:
    friend class mesh_elem_iter_struct_t<mesh_edge_iter_struct_t<mesh_t>, mesh_t, Edge>;

    /** @internal
     ** Get mesh element at index. */
    SKUNK_INLINE static auto get_elem_(mesh_t* mesh, uint_t edge_ind) {
        SKUNK_ASSERT(mesh != nullptr);
        return edge_ind < mesh->num_edges() ? &mesh->get_edge(edge_ind) : nullptr;
    }

public:

    /** Construct base edge iterator. */
    template<typename mesh_ptr_t>
    SKUNK_INLINE explicit mesh_edge_iter_struct_t(const mesh_ptr_t& mesh, uint_t edge_ind = 0)
        : mesh_elem_iter_struct_t<mesh_edge_iter_struct_t<mesh_t>, mesh_t, Edge>(mesh, edge_ind) {
    }
};  // class EdgeIter

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/** Base face iterator. */
template<typename mesh_t>
class mesh_face_iter_struct_t final :
    public mesh_elem_iter_struct_t<mesh_face_iter_struct_t<mesh_t>, mesh_t, Face> {
private:
    friend class mesh_elem_iter_struct_t<mesh_face_iter_struct_t<mesh_t>, mesh_t, Face>;

    /** @internal
     ** Get mesh element at index. */
    SKUNK_INLINE static auto get_elem_(mesh_t* mesh, uint_t face_ind) {
        SKUNK_ASSERT(mesh != nullptr);
        return face_ind < mesh->num_faces() ? &mesh->get_face(face_ind) : nullptr;
    }

public:

    /** Construct base face iterator. */
    template<typename mesh_ptr_t>
    SKUNK_INLINE explicit mesh_face_iter_struct_t(const mesh_ptr_t& mesh, uint_t face_ind = 0)
        : mesh_elem_iter_struct_t<mesh_face_iter_struct_t<mesh_t>, mesh_t, Face>(mesh, face_ind) {
    }

    /**************************************************************************/
    /**************************************************************************/

    /** Get connected inner cell. */
    SKUNK_INLINE auto get_inner_cell() const {
        return mesh_cell_iter_struct_t<mesh_t>(this->get_mesh(), (*this)->get_inner_cell());
    }
    /** Get connected outer cell. */
    SKUNK_INLINE auto get_outer_cell() const {
        return mesh_cell_iter_struct_t<mesh_t>(this->get_mesh(), (*this)->get_outer_cell());
    }
};  // class FaceIter

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

/** Base cell iterator. */
template<typename mesh_t>
class mesh_cell_iter_struct_t final :
    public mesh_elem_iter_struct_t<mesh_cell_iter_struct_t<mesh_t>, mesh_t, Cell> {
private:
    friend class mesh_elem_iter_struct_t<mesh_cell_iter_struct_t<mesh_t>, mesh_t, Cell>;

    /** @internal
     ** Get mesh element at index. */
    SKUNK_INLINE static auto get_elem_(mesh_t* mesh, uint_t cell_ind) {
        SKUNK_ASSERT(mesh != nullptr);
        return cell_ind < mesh->num_cells() ? &mesh->get_cell(cell_ind) : nullptr;
    }

public:

    /** Construct base cell iterator. */
    template<typename mesh_ptr_t>
    SKUNK_INLINE explicit mesh_cell_iter_struct_t(const mesh_ptr_t& mesh, uint_t cell_ind = 0)
        : mesh_elem_iter_struct_t<mesh_cell_iter_struct_t<mesh_t>, mesh_t, Cell>(mesh, cell_ind) {
    }
};  // class CellIter

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#endif  // ifndef MHD_MESH_ITER_HH
