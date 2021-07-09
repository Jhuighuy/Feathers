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

#include <libSkunkMesh/SkunkMesh.hh>
#include <libSkunkGeom/SkunkGeomEdge.hh>
#include <libSkunkGeom/SkunkGeomTriangle.hh>
#include <libSkunkMisc/SkunkMiscParallel.hh>
#include <libSkunkMisc/SkunkMiscReordering.hh>

#include <new>
#include <set>
#include <map>
#include <fstream>

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

#define assert1(x) do { if(!(x)) {throw std::runtime_error("hui");} } while(false)

/** @todo Refactor me! */
bool UMesh::read_tetgen(const char *path) {
    std::string line;

    std::ifstream node_file(path + std::string("node"));
    assert1(node_file.is_open());
    uint_t num_nodes = 0;
    node_file >> num_nodes;
    node_file >> m_dim;
    std::getline(node_file, line);
    for (uint_t i = 0; i < num_nodes; ++i) {
        uint_t node_ind = 0;
        vec3_t node_pos;
        node_file >> node_ind >> node_pos.x >> node_pos.y >> node_pos.z;
        std::getline(node_file, line);
        assert1(insert_node(mesh_node1_t(node_pos)) == node_ind);
    }

    std::ifstream face_file(path + std::string("face"));
    assert1(face_file.is_open());
    uint_t num_faces = 0;
    face_file >> num_faces;
    std::getline(face_file, line);
    for (uint_t i = 0; i < num_faces; ++i) {
        uint_t face_ind = 0;
        std::vector<uint_t> face_nodes(3);
        uint_t mark = 0;
        face_file >> face_ind >> face_nodes[0] >> face_nodes[1] >> face_nodes[2] >> mark;
        std::getline(face_file, line);
        assert1(insert_face(mesh_face_triangle3_t(face_nodes.begin(), face_nodes.end(), mark)) == face_ind);
    }

    std::ifstream cell_file(path + std::string("ele"));
    assert1(cell_file.is_open());
    uint_t num_cells = 0;
    cell_file >> num_cells;
    std::getline(cell_file, line);
    for (uint_t i = 0; i < num_cells; ++i) {
        uint_t cell_ind = 0;
        std::vector<uint_t> cell_nodes(4);
        cell_file >> cell_ind >> cell_nodes[0] >> cell_nodes[1] >> cell_nodes[2] >> cell_nodes[3];
        std::getline(cell_file, line);
        assert1(insert_cell(mesh_cell_tetrahedron4_t(cell_nodes.begin(), cell_nodes.end())) == cell_ind);
    }

    generate_faces();
    generate_boundary_cells();
    reorder_faces();
    return true;
}
bool UMesh::read_triangle(const char *path) {
    std::string line;

    std::ifstream node_file(path + std::string("node"));
    assert1(node_file.is_open());
    uint_t num_nodes = 0;
    node_file >> num_nodes;
    node_file >> m_dim;
    std::getline(node_file, line);
    for (uint_t i = 0; i < num_nodes; ++i) {
        uint_t node_ind = 0;
        vec3_t node_pos;
        node_file >> node_ind >> node_pos.x >> node_pos.y;
        std::getline(node_file, line);
        assert1(insert_node(mesh_node1_t(node_pos)) == node_ind);
    }

    std::ifstream face_file(path + std::string("edge"));
    assert1(face_file.is_open());
    uint_t num_faces = 0;
    face_file >> num_faces;
    std::getline(face_file, line);
    for (uint_t i = 0; i < num_faces; ++i) {
        uint_t face_ind = 0;
        std::vector<uint_t> face_nodes(2);
        uint_t mark = 0;
        face_file >> face_ind >> face_nodes[0] >> face_nodes[1] >> mark;
        std::getline(face_file, line);
        assert1(insert_face(mesh_face_segment2_t(face_nodes.begin(), face_nodes.end(), mark)) == face_ind);
    }

    std::ifstream cell_file(path + std::string("ele"));
    assert1(cell_file.is_open());
    uint_t num_cells = 0;
    cell_file >> num_cells;
    std::getline(cell_file, line);
    for (uint_t i = 0; i < num_cells; ++i) {
        uint_t cell_ind = 0;
        std::vector<uint_t> cell_nodes(3);
        cell_file >> cell_ind >> cell_nodes[0] >> cell_nodes[1] >> cell_nodes[2];
        std::getline(cell_file, line);
        assert1(insert_cell(mesh_cell_triangle3_t(cell_nodes.begin(), cell_nodes.end())) == cell_ind);
    }

    generate_faces();
    generate_boundary_cells();
    //generate_ghost_cells_and_fix_boundaries_();
    reorder_faces();
    return true;
}

/** @todo Refactor me! */
bool UMesh::read_su2(const char* path) {
    std::ifstream file(path);
    std::string line;

    /* Parse dimension. */
    std::string ndime;
    std::getline(file, ndime, '=');
    if (ndime != "NDIME") {
        assert(false); return false;
    }

    std::string ndime_str;
    std::getline(file, ndime_str);
    m_dim = std::atoi(ndime_str.c_str());
    if (!(1 <= m_dim && m_dim <= 3)) {
        assert(false); return false;
    }
    /*
     * Parse nodes.
     */
    std::string npoin;
    std::getline(file, npoin, '=');
    if (npoin != "NPOIN") {
        assert(false); return false;
    }
    std::string npoin_str;
    std::getline(file, npoin_str);
    const uint_t num_nodes = std::atoi(npoin_str.c_str());
    for (uint_t node_ind = 0;
        node_ind < num_nodes; ++node_ind) {
        vec3_t node_coord;
        file >> node_coord.x;
        if (m_dim >= 2) {
            file >> node_coord.y;
            if (m_dim == 3) {
                file >> node_coord.z;
            }
        }
        insert_node(skunk::mesh_node1_t(node_coord));
        std::getline(file, line);
    }
    /*
     * Parse cells.
     */
    std::string nelem;
    std::getline(file, nelem, '=');
    if (nelem != "NELEM") {
        assert(false); return false;
    }
    std::string nelem_str;
    std::getline(file, nelem_str);
    const uint_t num_cells = std::atoi(nelem_str.c_str());
    for (uint_t cell_ind = 0;
                    cell_ind < num_cells; ++cell_ind) {
        uint_t nelem_typ;
        file >> nelem_typ;
        switch (nelem_typ) {
            /* Triangle cell. */
            case 5: {
                std::array<uint_t, 3> cell_nodes{};
                file >> cell_nodes[0];
                file >> cell_nodes[1];
                file >> cell_nodes[2];
                insert_cell(skunk::mesh_cell_triangle3_t(std::begin(cell_nodes), std::end(cell_nodes)));
                std::getline(file, line);
                break;
            }
            default: {
                assert(false); return false;
            }
        }
    }
    /*
     * Read boundary faces.
     */
    std::string nmark;
    std::getline(file, nmark, '=');
    if (nmark != "NMARK") {
        assert(false); return false;
    }
    std::string nmark_str;
    std::getline(file, nmark_str);
    const uint_t num_marks = std::atoi(nmark_str.c_str());
    for (uint_t mark_ind = 0; mark_ind < num_marks; ++mark_ind) {
        std::string marker_tag;
        std::getline(file, marker_tag, '=');
        if (marker_tag != "MARKER_TAG") {
            assert(false); return false;
        }
        std::string marker_tag_str;
        std::getline(file, marker_tag_str);
        const uint_t marker_val = std::atoi(marker_tag_str.c_str());
        std::string marker_elems;
        std::getline(file, marker_elems, '=');
        if (marker_elems != "MARKER_ELEMS") {
            assert(false); return false;
        }
        std::string marker_elems_str;
        std::getline(file, marker_elems_str);
        const uint_t num_faces = std::atoi(marker_elems_str.c_str());
        for (uint_t face_loc = 0; face_loc < num_faces; ++face_loc) {
            uint_t nelem_typ;
            file >> nelem_typ;
            switch (nelem_typ) {
                case 3: {
                    std::array<uint_t, 2> face_nodes{};
                    file >> face_nodes[0];
                    file >> face_nodes[1];
                    insert_face(skunk::mesh_face_segment2_t(face_nodes.cbegin(), face_nodes.cend(), marker_val));
                    std::getline(file, line);
                    break;
                }
                default: {
                    assert(false);
                    return false;
                }
            }
        }
    }
    generate_faces();
    generate_edges();
    //generate_ghost_cells_and_fix_boundaries_();
    generate_boundary_cells();
    reorder_faces();
    return true;
}

namespace skunk {

/** Get minimal edge length. */
SKUNK_EXTERN real_t mesh_t::get_min_edge_length() const {
    return for_range_min(0u, num_edges(), 1e+6, [&](uint_t edge_ind) {
        return get_edge_length(edge_ind);
    });
}   // mesh_t::get_min_edge_length
/** Get maximal edge length. */
SKUNK_EXTERN real_t mesh_t::get_max_edge_length() const {
    return for_range_max(0u, num_edges(), 0e+0, [&](uint_t edge_ind) {
        return get_edge_length(edge_ind);
    });
}   // mesh_t::get_max_edge_length

/** Get minimal face area. */
SKUNK_EXTERN real_t mesh_t::get_min_face_area() const {
    return for_range_min(0u, num_faces(), 1e+6, [&](uint_t face_ind) {
        return get_face_area(face_ind);
    });
}   // mesh_t::get_min_face_area
/** Get maximal face area. */
SKUNK_EXTERN real_t mesh_t::get_max_face_area() const {
    return for_range_max(0u, num_faces(), 0e+0, [&](uint_t face_ind) {
        return get_face_area(face_ind);
    });
}   // mesh_t::get_max_face_area

/** Get minimal cell volume. */
SKUNK_EXTERN real_t mesh_t::get_min_cell_volume() const {
    return for_range_min(0u, num_cells(), 1e+6, [&](uint_t cell_ind) {
        return get_cell_volume(cell_ind);
    });
}   // mesh_t::get_min_cell_volume
/** Get maximal cell volume. */
SKUNK_EXTERN real_t mesh_t::get_max_cell_volume() const {
    return for_range_max(0u, num_cells(), 0e+0, [&](uint_t cell_ind) {
        return get_cell_volume(cell_ind);
    });
}   // mesh_t::get_max_cell_volume

/** Compute mesh orthogonality.
 ** Mesh orthogonality is defined as a. */
SKUNK_EXTERN real_t mesh_t::get_orthogonality() const {
    return for_range_min(begin_face(0), end_face(0), 1e+6, [&](uint_t face_ind) {
        const mesh_face_t& face = get_face(face_ind);
        const vec3_t direction = get_cell_center_position(face.get_outer_cell()) -
                                 get_cell_center_position(face.get_inner_cell());
        const real_t orth = face.get_normal().inner_prod(direction)/direction.len();
        return orth;
    });
}   // mesh_t::get_orthogonality

/**************************************************************************/
/**************************************************************************/

/** Rescale the mesh. */
SKUNK_EXTERN void mesh_t::scale(real_t scale) {
    /* Compute the scaling factors for face areas and cell volumes. */
    real_t scale_area, scale_volume;
    switch (m_dim) {
    case 1:
        scale_area   = scale;
        scale_volume = scale;
        break;
    case 2:
        scale_area   = scale;
        scale_volume = std::pow(scale, 2);
        break;
    case 3:
        scale_area   = std::pow(scale, 2);
        scale_volume = std::pow(scale, 3);
        break;
    default:
        SKUNK_ASSERT_FALSE("Invalid amount of dimensions.");
    }
    /* Rescale nodes. */
    for_range(0u, num_nodes(), [&](uint_t node_ind) {
        mesh_node_t& node = get_node(node_ind);
        node.set_position(scale*node.get_position());
    });
    /* Rescale edges. */
    for_range(0u, num_edges(), [&](uint_t edge_ind) {
        mesh_edge_t& edge = get_edge(edge_ind);
        edge.set_length(scale*edge.get_length());
    });
    /* Rescale faces. */
    for_range(0u, num_faces(), [&](uint_t face_ind) {
        mesh_face_t& face = get_face(face_ind);
        face.set_area(scale_area*face.get_area());
        face.set_center_position(scale*face.get_center_position());
    });
    /* Rescale cells. */
    for_range(0u, num_cells(), [&](uint_t cell_ind) {
        mesh_cell_t& cell = get_cell(cell_ind);
        cell.set_volume(scale_volume*cell.get_volume());
        cell.set_center_position(scale*cell.get_center_position());
    });
}   // mesh_t::scale

/**************************************************************************/
/**************************************************************************/

/** Insert a new node into the mesh.
 ** @param node_ Edge nodes.
 ** @returns Index of the inserted node. */
SKUNK_EXTERN uint_t mesh_t::insert_node(const mesh_node1_t& node_) {
    /** @todo Refactor me! */
    const uint_t node_ind = allocate_node_(node_);
    const mesh_node_t& node = get_node(node_ind);
    const vec3_t& node_geometry = node.get_position();
    if (m_dim <= 2) {
        SKUNK_ASSERT(node_geometry.z == 0.0);
        if (m_dim == 1) {
            SKUNK_ASSERT(node_geometry.y == 0.0);
        }
    }
    return node_ind;
}   // mesh_t::insert_node

/** Insert a new edge into the mesh.
 ** @param edge_ Edge nodes.
 ** @returns Index of the inserted edge. */
/** @{ */
/*
 * Dummy 1D/2D edge.
 */
SKUNK_EXTERN uint_t mesh_t::insert_edge(const mesh_edge_node1_t& edge_) {
    const uint_t edge_ind = allocate_edge_(edge_);
    /** @todo Implement me! */
    return edge_ind;
}   // mesh_t::insert_edge
/*
 * Segment edge.
 */
SKUNK_EXTERN uint_t mesh_t::insert_edge(const mesh_edge_segment2_t& edge_) {
    const uint_t edge_ind = allocate_edge_(edge_);
    /** @todo Implement me! */
    return edge_ind;
}   // mesh_t::insert_edge
SKUNK_EXTERN uint_t mesh_t::insert_edge(const mesh_edge_segment3_t& edge_) {
    const uint_t edge_ind = allocate_edge_(edge_);
    /** @todo Implement me! */
    return edge_ind;
}   // mesh_t::insert_edge
/** @} */

/** Insert a new edge into the mesh.
 ** Type of the edge is detected by the set of nodes and dimension.
 ** @param edge_nodes Set of the edge nodes.
 ** @param mark Boundary mark.
 ** @param dim Overridden dimension.
 ** @returns Index of the inserted edge. */
SKUNK_EXTERN uint_t mesh_t::insert_edge(const std::vector<uint_t>& edge_nodes, uint_t mark, uint_t dim) {
    if (dim == 0) {
        dim = m_dim;
    }
    SKUNK_ASSERT(1 <= dim && dim <= 3);
    if (dim <= 2) {
        switch (edge_nodes.size()) {
        /* Dummy 1D/2D edge. */
        case 1:
            return insert_edge(mesh_edge_node1_t(edge_nodes.cbegin(), edge_nodes.cend(), mark));
        default:;
        }
    } else if (dim == 3) {
        switch (edge_nodes.size()) {
        /* Segment edge. */
        case 2:
            return insert_edge(mesh_edge_segment2_t(edge_nodes.cbegin(), edge_nodes.cend(), mark));
        case 3:
            return insert_edge(mesh_edge_segment3_t(edge_nodes.cbegin(), edge_nodes.cend(), mark));
        default:;
        }
    }
    SKUNK_ASSERT_FALSE("Invalid edge nodes: cell type cannot be uniquely matched.");
}   // mesh_t::insert_edge

/** Insert a new face into the mesh.
 ** @param face_ Face nodes.
 ** @returns Index of the inserted face. */
/** @{ */
/*
 * Dummy 1D face.
 */
SKUNK_EXTERN uint_t mesh_t::insert_face(const mesh_face_node1_t& face_) {
    /** @todo Refactor me! */
    /* Insert the face. */
    const uint_t face_ind = allocate_face_(face_);
    /* Set up the face geometry. */
    mesh_face_t& face = get_face(face_ind);
    face.set_area(1.0);
    face.set_center_position(get_node_position(face.begin_node()[0]));
    face.set_normal({1.0, 0.0, 0.0});
    return face_ind;
}   // mesh_t::insert_face
/*
 * Segment face.
 */
SKUNK_EXTERN uint_t mesh_t::insert_face(const mesh_face_segment2_t& face_) {
    /** @todo Refactor me! */
    /* Insert the face. */
    const uint_t face_ind = allocate_face_(face_);
    /* Set up the face geometry. */
    mesh_face_t& face = get_face(face_ind);
    const mhd_seg2_t face_geometry = {
        vec2_t(get_node_position(face.begin_node()[0])),
        vec2_t(get_node_position(face.begin_node()[1])),
    };
    face.set_area(mhd_seg2_t::len(face_geometry));
    face.set_center_position(vec3_t(0.5*(face_geometry.p1 + face_geometry.p2)));
    face.set_normal(mhd_seg2_t::normal3(face_geometry));
    return face_ind;
}   // mesh_t::insert_face
SKUNK_EXTERN uint_t mesh_t::insert_face(const mesh_face_segment3_t& face_) {
    SKUNK_NOT_IMPLEMENTED();
}   // mesh_t::insert_face
/*
 * Triangular face.
 */
SKUNK_EXTERN uint_t mesh_t::insert_face(const mesh_face_triangle3_t& face_) {
    /** @todo Refactor me! */
    /* Insert the face. */
    const uint_t face_ind = allocate_face_(face_);
    /* Set up the face geometry. */
    mesh_face_t& face = get_face(face_ind);
    const mhd_tri3_t face_geometry = {
        get_node_position(face.begin_node()[0]),
        get_node_position(face.begin_node()[1]),
        get_node_position(face.begin_node()[2]),
    };
    face.set_area(mhd_tri3_t::area(face_geometry));
    face.set_center_position(mhd_tri3_t::barycenter(face_geometry));
    face.set_normal(mhd_tri3_t::normal(face_geometry));
    return face_ind;
}   // mesh_t::insert_face
SKUNK_EXTERN uint_t mesh_t::insert_face(const mesh_face_triangle6_t& face_) {
    SKUNK_NOT_IMPLEMENTED();
}   // mesh_t::insert_face
/*
 * Quadrangular face.
 */
SKUNK_EXTERN uint_t mesh_t::insert_face(const mesh_face_quadrangle4_t& face_) {
    /** @todo Refactor me! */
    /* Insert the face. */
    const uint_t face_ind = allocate_face_(face_);
    /* Set up the face geometry. */
    mesh_face_t& face = get_face(face_ind);
    const mhd_tri3_t face_geometry1 = {
        get_node_position(face.begin_node()[0]),
        get_node_position(face.begin_node()[1]),
        get_node_position(face.begin_node()[2]),
    };
    const mhd_tri3_t face_geometry2 = {
        get_node_position(face.begin_node()[2]),
        get_node_position(face.begin_node()[3]),
        get_node_position(face.begin_node()[0]),
    };
    const auto a1 = mhd_tri3_t::area(face_geometry1);
    const auto a2 = mhd_tri3_t::area(face_geometry2);
    //SKUNK_ASSERT(a1 > 0.0);
    //SKUNK_ASSERT(a2 > 0.0);
    //SKUNK_ASSERT((a1 + a2) > 0.0);
    const auto a = a1 + a2;
    if (a != 0.0) {
        face.set_area(a1 + a2);
        face.set_center_position(a1/(a1 + a2)*mhd_tri3_t::barycenter(face_geometry1) +
                                 a2/(a1 + a2)*mhd_tri3_t::barycenter(face_geometry2));
    } else {
        face.set_area(0.0);
        face.set_center_position(0.5*mhd_tri3_t::barycenter(face_geometry1) +
                                 0.5*mhd_tri3_t::barycenter(face_geometry2));
    }
    auto nn = a1*mhd_tri3_t::normal(face_geometry1) + a2*mhd_tri3_t::normal(face_geometry2);
    face.set_normal(nn*safe_inv(nn.len()));
    return face_ind;
}   // mesh_t::insert_face
SKUNK_EXTERN uint_t mesh_t::insert_face(const mesh_face_quadrangle8_t& face_) {
    SKUNK_NOT_IMPLEMENTED();
}   // mesh_t::insert_face
SKUNK_EXTERN uint_t mesh_t::insert_face(const mesh_face_quadrangle9_t& face_) {
    SKUNK_NOT_IMPLEMENTED();
}   // mesh_t::insert_face
/** @} */

/** Insert a new face into the mesh.
 ** Type of the face is detected by the set of nodes and dimension.
 ** @param face_nodes Set of the face nodes.
 ** @param mark Boundary mark.
 ** @param dim Overridden dimension.
 ** @returns Index of the inserted face. */
SKUNK_EXTERN uint_t mesh_t::insert_face(const std::vector<uint_t>& face_nodes, uint_t mark, uint_t dim) {
    if (dim == 0) {
        dim = m_dim;
    }
    SKUNK_ASSERT(1 <= dim && dim <= 3);
    if (dim == 1) {
        switch (face_nodes.size()) {
        /* Dummy 1D face. */
        case 1:
            return insert_face(mesh_face_node1_t(face_nodes.cbegin(), face_nodes.cend(), mark));
        default:;
        }
    } else if (dim == 2) {
        switch (face_nodes.size()) {
        /* Segment face. */
        case 2:
            return insert_face(mesh_face_segment2_t(face_nodes.cbegin(), face_nodes.cend(), mark));
        case 3:
            return insert_face(mesh_face_segment3_t(face_nodes.cbegin(), face_nodes.cend(), mark));
        default:;
        }
    } else if (dim == 3) {
        switch (face_nodes.size()) {
        /* Triangular face. */
        case 3:
            return insert_face(mesh_face_triangle3_t(face_nodes.cbegin(), face_nodes.cend(), mark));
        case 6:
            return insert_face(mesh_face_triangle6_t(face_nodes.cbegin(), face_nodes.cend(), mark));
        /* Quadrangular face. */
        case 4:
            return insert_face(mesh_face_quadrangle4_t(face_nodes.cbegin(), face_nodes.cend(), mark));
        case 8:
            return insert_face(mesh_face_quadrangle8_t(face_nodes.cbegin(), face_nodes.cend(), mark));
        case 9:
            return insert_face(mesh_face_quadrangle9_t(face_nodes.cbegin(), face_nodes.cend(), mark));
        default:;
        }
    }
    SKUNK_ASSERT_FALSE("Invalid face nodes: face type cannot be uniquely matched.");
}   // mesh_t::insert_face

/** Initialize cell geometry. */
/** @{ */
template<typename geom_t, typename mesh_cell_t>
SKUNK_INLINE void init_cell_geometry_2d_(const mesh_t& mesh, mesh_cell_t& cell) {
    geom_t cell_geometry{};
    for (uint_t node_loc = 0; node_loc < cell.num_nodes(); ++node_loc) {
        const uint_t node_ind = cell.begin_node()[node_loc];
        cell_geometry.node(node_loc) = vec2_t(mesh.get_node_position(node_ind));
    }
    cell.set_center_position(cell_geometry.barycenter());
    cell.set_volume(cell_geometry.area());
}   // init_cell_geometry_
template<typename geom_t, typename mesh_cell_t>
SKUNK_INLINE void init_cell_geometry_23d_(const mesh_t& mesh, mesh_cell_t& cell) {
    geom_t cell_geometry{};
    for (uint_t node_loc = 0; node_loc < cell.num_nodes(); ++node_loc) {
        const uint_t node_ind = cell.begin_node()[node_loc];
        cell_geometry.node(node_loc) = (mesh.get_node_position(node_ind));
    }
    cell.set_center_position(cell_geometry.barycenter(cell_geometry));
    cell.set_volume(cell_geometry.area(cell_geometry));
}   // init_cell_geometry_
template<typename geom_t, typename mesh_cell_t>
SKUNK_INLINE void init_cell_geometry_3d_(const mesh_t& mesh, mesh_cell_t& cell) {
    geom_t cell_geometry{};
    for (uint_t node_loc = 0; node_loc < cell.num_nodes(); ++node_loc) {
        const uint_t node_ind = cell.begin_node()[node_loc];
        cell_geometry.node(node_loc) = mesh.get_node_position(node_ind);
    }
    cell.set_center_position(cell_geometry.barycenter());
    cell.set_volume(cell_geometry.volume());
}   // init_cell_geometry_
/** @} */

/** Insert a new cell into the mesh.
 ** @param cell_ Cell nodes.
 ** @returns Index of the inserted face. */
/** @{ */
/*
 * Segment cell.
 */
SKUNK_EXTERN uint_t mesh_t::insert_cell(const mesh_cell_segment2_t& cell_) {
    const uint_t cell_ind = allocate_cell_(cell_);
    mesh_cell_t& cell = get_cell(cell_ind);
    /*
     * Calculate cell area and verify orientation.
     */
    const mhd_seg2_t cell_geometry = {
            vec2_t(get_node_position(cell.begin_node()[0])),
            vec2_t(get_node_position(cell.begin_node()[1])),
    };
    cell.set_center_position(vec3_t(0.5*(cell_geometry.p1 + cell_geometry.p2)));
    cell.set_volume(mhd_seg2_t::len(cell_geometry));
    return cell_ind;
}   // mesh_t::insert_cell
SKUNK_EXTERN uint_t mesh_t::insert_cell(const mesh_cell_segment3_t& cell_) {
    SKUNK_NOT_IMPLEMENTED();
}   // mesh_t::insert_face
/*
 * Triangular cell.
 */
SKUNK_EXTERN uint_t mesh_t::insert_cell(const mesh_cell_triangle3_t& cell_) {
    const uint_t cell_ind = allocate_cell_(cell_);
    mesh_cell_t& cell = get_cell(cell_ind);
    /** @todo A 3D triangle? */
    init_cell_geometry_23d_<mhd_tri3d>(*this, cell);
    return cell_ind;
}   // UMesh::insert_cell
SKUNK_EXTERN uint_t mesh_t::insert_cell(const mesh_cell_triangle6_t& cell_) {
    SKUNK_NOT_IMPLEMENTED();
}   // mesh_t::insert_face
/*
 * Quadrangular cell.
 */
SKUNK_EXTERN uint_t mesh_t::insert_cell(const mesh_cell_quadrangle4_t& cell_) {
    const uint_t cell_ind = allocate_cell_(cell_);
    mesh_cell_t& cell = get_cell(cell_ind);
    /** @todo A 3D quadrangle? */
    init_cell_geometry_2d_<geom_quad2_t>(*this, cell);
    return cell_ind;
}   // mesh_t::insert_cell
SKUNK_EXTERN uint_t mesh_t::insert_cell(const mesh_cell_quadrangle8_t& cell_) {
    SKUNK_NOT_IMPLEMENTED();
}   // mesh_t::insert_face
SKUNK_EXTERN uint_t mesh_t::insert_cell(const mesh_cell_quadrangle9_t& cell_) {
    SKUNK_NOT_IMPLEMENTED();
}   // mesh_t::insert_face
/*
 * Tetrahedral cell.
 */
SKUNK_EXTERN uint_t mesh_t::insert_cell(const mesh_cell_tetrahedron4_t& cell_) {
    const uint_t cell_ind = allocate_cell_(cell_);
    mesh_cell_t& cell = get_cell(cell_ind);
    /** @todo A 3D triangle? */
    init_cell_geometry_3d_<mhd_tetr3d_t>(*this, cell);
    return cell_ind;
}   // mesh_t::insert_face
SKUNK_EXTERN uint_t mesh_t::insert_cell(const mesh_cell_tetrahedron10_t& cell_) {
    SKUNK_NOT_IMPLEMENTED();
}   // mesh_t::insert_face
/*
 * Pyramidal cell.
 */
SKUNK_EXTERN uint_t mesh_t::insert_cell(const mesh_cell_pyramid5_t& cell_) {
    const uint_t cell_ind = allocate_cell_(cell_);
    mesh_cell_t& cell = get_cell(cell_ind);
    mhd_pyra3d_t cell_shape{
            get_node_position(cell.begin_node()[0]),
            get_node_position(cell.begin_node()[1]),
            get_node_position(cell.begin_node()[2]),
            get_node_position(cell.begin_node()[3]),
            get_node_position(cell.begin_node()[4]),
    };
    cell.set_volume(cell_shape.volume());
    cell.set_center_position(cell_shape.barycenter());
    return cell_ind;
}
SKUNK_EXTERN uint_t mesh_t::insert_cell(const mesh_cell_pyramid14_t& cell_) {
    SKUNK_NOT_IMPLEMENTED();
}   // mesh_t::insert_face
/*
 * Pentahedral cell.
 */
SKUNK_EXTERN uint_t mesh_t::insert_cell(const mesh_cell_pentahedron6_t& cell_) {
    SKUNK_NOT_IMPLEMENTED();
}   // mesh_t::insert_face
SKUNK_EXTERN uint_t mesh_t::insert_cell(const mesh_cell_pentahedron15_t& cell_) {
    SKUNK_NOT_IMPLEMENTED();
}   // mesh_t::insert_face
SKUNK_EXTERN uint_t mesh_t::insert_cell(const mesh_cell_pentahedron18_t& cell_) {
    SKUNK_NOT_IMPLEMENTED();
}   // mesh_t::insert_face
/*
 * Hexahedral cell.
 */
SKUNK_EXTERN uint_t mesh_t::insert_cell(const skunk::mesh_cell_hexahedron8_t& cell_) {
    const uint_t cell_ind = allocate_cell_(cell_);
    mesh_cell_t& cell = get_cell(cell_ind);
    mhd_hexa3d_t cell_shape{
            get_node_position(cell.begin_node()[0]),
            get_node_position(cell.begin_node()[1]),
            get_node_position(cell.begin_node()[2]),
            get_node_position(cell.begin_node()[3]),
            get_node_position(cell.begin_node()[4]),
            get_node_position(cell.begin_node()[5]),
            get_node_position(cell.begin_node()[6]),
            get_node_position(cell.begin_node()[7]),
    };
    cell.set_volume(cell_shape.volume());
    cell.set_center_position(cell_shape.barycenter());
    return cell_ind;
}
SKUNK_EXTERN uint_t mesh_t::insert_cell(const mesh_cell_hexahedron20_t& cell_) {
    SKUNK_NOT_IMPLEMENTED();
}   // mesh_t::insert_face
SKUNK_EXTERN uint_t mesh_t::insert_cell(const mesh_cell_hexahedron27_t& cell_) {
    SKUNK_NOT_IMPLEMENTED();
}   // mesh_t::insert_face
/** @} */

/** Insert a new cell into the mesh.
 ** Type of the cell is detected by the set of nodes and dimension.
 ** @param cell_nodes Set of the cell nodes.
 ** @param mark Boundary mark.
 ** @param dim Overridden dimension.
 ** @returns Index of the inserted cell. */
SKUNK_EXTERN uint_t mesh_t::insert_cell(const std::vector<uint_t>& cell_nodes, uint_t mark, uint_t dim) {
    if (dim == 0) {
        dim = m_dim;
    }
    SKUNK_ASSERT(1 <= dim && dim <= 3);
    if (dim == 1) {
        switch (cell_nodes.size()) {
        /* Segment cell. */
        case 2:
            return insert_cell(mesh_cell_segment2_t(cell_nodes.cbegin(), cell_nodes.cend(), mark));
        case 3:
            return insert_cell(mesh_cell_segment3_t(cell_nodes.cbegin(), cell_nodes.cend(), mark));
        default:;
        }
    } else if (dim == 2) {
        switch (cell_nodes.size()) {
        /* Triangular cell. */
        case 3:
            return insert_cell(mesh_cell_triangle3_t(cell_nodes.cbegin(), cell_nodes.cend(), mark));
        case 6:
            return insert_cell(mesh_cell_triangle6_t(cell_nodes.cbegin(), cell_nodes.cend(), mark));
        /* Quadrangular cell. */
        case 4:
            return insert_cell(mesh_cell_quadrangle4_t(cell_nodes.cbegin(), cell_nodes.cend(), mark));
        case 8:
            return insert_cell(mesh_cell_quadrangle8_t(cell_nodes.cbegin(), cell_nodes.cend(), mark));
        case 9:
            return insert_cell(mesh_cell_quadrangle9_t(cell_nodes.cbegin(), cell_nodes.cend(), mark));
        default:;
        }
    } else if (dim == 3) {
        switch (cell_nodes.size()) {
        /* Tetrahedral cell. */
        case 4:
            return insert_cell(mesh_cell_tetrahedron4_t(cell_nodes.cbegin(), cell_nodes.cend(), mark));
        case 10:
            return insert_cell(mesh_cell_tetrahedron10_t(cell_nodes.cbegin(), cell_nodes.cend(), mark));
        /* Pyramidal cell. */
        case 5:
            return insert_cell(mesh_cell_pyramid5_t(cell_nodes.cbegin(), cell_nodes.cend(), mark));
        case 14:
            return insert_cell(mesh_cell_pyramid14_t(cell_nodes.cbegin(), cell_nodes.cend(), mark));
        /* Pentahedral cell. */
        case 6:
            return insert_cell(mesh_cell_pentahedron6_t(cell_nodes.cbegin(), cell_nodes.cend(), mark));
        case 15:
            return insert_cell(mesh_cell_pentahedron15_t(cell_nodes.cbegin(), cell_nodes.cend(), mark));
        case 18:
            return insert_cell(mesh_cell_pentahedron18_t(cell_nodes.cbegin(), cell_nodes.cend(), mark));
        /* Hexahedral cell. */
        case 8:
            return insert_cell(mesh_cell_hexahedron8_t(cell_nodes.cbegin(), cell_nodes.cend(), mark));
        case 20:
            return insert_cell(mesh_cell_hexahedron20_t(cell_nodes.cbegin(), cell_nodes.cend(), mark));
        case 27:
            return insert_cell(mesh_cell_hexahedron27_t(cell_nodes.cbegin(), cell_nodes.cend(), mark));
        default:;
        }
    }
    SKUNK_ASSERT(!"Invalid cell nodes: cell type cannot be uniquely matched.");
}   // mesh_t::insert_cell

/** @internal
 ** Allocate and insert the new element into the storage. */
template<typename mesh_elem_struct_t>
SKUNK_INLINE uint_t allocate_element_(const mesh_elem_struct_t& elem, size_t alignment,
                                      std::vector<byte_t>& elem_storage, std::vector<size_t>& elem_offsets) {
    SKUNK_ASSERT(alignment > 0);
    /* Align size of the element. */
    size_t delta_elem = sizeof(elem);
    delta_elem += (alignment - delta_elem%alignment)%alignment;
    /* Allocate and initialize a new element. */
    elem_offsets.push_back(elem_storage.size());
    elem_storage.resize(elem_storage.size() + delta_elem);
    new (elem_storage.data() + elem_offsets.back()) mesh_elem_struct_t(elem);
    const auto elem_ind = static_cast<uint_t>(elem_offsets.size()) - 1;
    return elem_ind;
}   // allocate_element_
/** @internal
 ** Allocate and insert a new node into the mesh. */
template<typename mesh_node_struct_t>
SKUNK_INLINE uint_t mesh_t::allocate_node_(const mesh_node_struct_t& node, size_t alignment) {
    return allocate_element_(node, alignment, m_node_storage, m_node_offsets);
}   // mesh_t::allocate_node_
/** @internal
 ** Allocate and insert a new edge into the mesh. */
template<typename mesh_edge_struct_t>
SKUNK_INLINE uint_t mesh_t::allocate_edge_(const mesh_edge_struct_t& edge, size_t alignment) {
    return allocate_element_(edge, alignment, m_edge_storage, m_edge_offsets);
}   // mesh_t::allocate_edge_
/** @internal
 ** Allocate and insert a new face into the mesh. */
template<typename mesh_face_struct_t>
SKUNK_INLINE uint_t mesh_t::allocate_face_(const mesh_face_struct_t& face, size_t alignment) {
    return allocate_element_(face, alignment, m_face_storage, m_face_offsets);
}   // mesh_t::allocate_face_
/** @internal
 ** Allocate and insert a new cell into the mesh. */
template<typename mesh_cell_struct_t>
SKUNK_INLINE uint_t mesh_t::allocate_cell_(const mesh_cell_struct_t& cell, size_t alignment) {
    return allocate_element_(cell, alignment, m_cell_storage, m_cell_offsets);
}   // mesh_t::allocate_cell_

/**************************************************************************/
/**************************************************************************/

namespace {
/** @internal
 ** Functor that replaces connected element indices when an element is erased. */
class mesh_erase_func_t {
private:
    uint_t m_elem_ind;
    uint_t m_elem_replacement_ind;
public:
    SKUNK_INLINE explicit mesh_erase_func_t(uint_t elem_ind, uint_t elem_replacement_ind = -1)
        : m_elem_ind(elem_ind),
          m_elem_replacement_ind(elem_replacement_ind) {
        if (m_elem_replacement_ind > m_elem_ind) {
            m_elem_replacement_ind -= 1;
        }
    }
    SKUNK_INLINE void operator()(uint_t& elem_ind) const {
        if (elem_ind == m_elem_ind) {
            elem_ind = m_elem_replacement_ind;
        } else if (elem_ind > m_elem_ind) {
            elem_ind -= 1;
        }
    }
};  // class mesh_erase_func_t
}   // namespace

/** Erase the node.
 ** All connections from other mesh elements would be set to -1. */
SKUNK_EXTERN void mesh_t::erase_node(uint_t node_ind, uint_t replacement_node_ind) {
    SKUNK_ASSERT(node_ind < num_nodes());
    SKUNK_ASSERT(replacement_node_ind < num_nodes() || replacement_node_ind == npos);
    /* Fix connections with edges. */
    for_range(0u, num_edges(), [&](uint_t edge_ind) {
        mesh_edge_t& edge = get_edge(edge_ind);
        std::for_each(edge.begin_node(), edge.end_node(), mesh_erase_func_t(node_ind, replacement_node_ind));
    });
    /* Fix connections with faces. */
    for_range(0u, num_faces(), [&](uint_t face_ind) {
        mesh_face_t& face = get_face(face_ind);
        std::for_each(face.begin_node(), face.end_node(), mesh_erase_func_t(node_ind, replacement_node_ind));
    });
    /* Fix connections with cells. */
    for_range(0u, num_cells(), [&](uint_t cell_ind) {
        mesh_cell_t& cell = get_cell(cell_ind);
        std::for_each(cell.begin_node(), cell.end_node(), mesh_erase_func_t(node_ind, replacement_node_ind));
    });
    /* Deallocate the node. */
    deallocate_node_(node_ind);
}   // mesh_t::erase_node
/** Erase the edge.
 ** @param replacement_edge_ind Edge to connect links of adjacent elements. */
SKUNK_EXTERN void mesh_t::erase_edge(uint_t edge_ind, uint_t replacement_edge_ind) {
    SKUNK_ASSERT(edge_ind < num_edges());
    SKUNK_ASSERT(replacement_edge_ind < num_edges() || replacement_edge_ind == npos);
    /* Fix connections with faces. */
    for_range(0u, num_faces(), [&](uint_t face_ind) {
        mesh_face_t& face = get_face(face_ind);
        std::for_each(face.begin_edge(), face.end_edge(), mesh_erase_func_t(edge_ind, replacement_edge_ind));
    });
    /* Deallocate the edge. */
    deallocate_edge_(edge_ind);
}   // mesh_t::erase_edge
/** Erase the face.
 ** @param replacement_face_ind Face to connect links of adjacent elements. */
SKUNK_EXTERN void mesh_t::erase_face(uint_t face_ind, uint_t replacement_face_ind) {
    SKUNK_ASSERT(face_ind < num_faces());
    SKUNK_ASSERT(replacement_face_ind < num_faces() || replacement_face_ind == npos);
    /* Fix connections with cells. */
    for_range(0u, num_cells(), [&](uint_t cell_ind) {
        mesh_cell_t& cell = get_cell(cell_ind);
        std::for_each(cell.begin_face(), cell.end_face(), mesh_erase_func_t(face_ind, replacement_face_ind));
    });
    /* Deallocate the face. */
    deallocate_face_(face_ind);
}   // mesh_t::erase_face
/** Erase the cell.
 ** @param replacement_cell_ind Cell to connect links of adjacent elements. */
SKUNK_EXTERN void mesh_t::erase_cell(uint_t cell_ind, uint_t replacement_cell_ind) {
    SKUNK_ASSERT(cell_ind < num_cells());
    SKUNK_ASSERT(replacement_cell_ind < num_cells() || replacement_cell_ind == npos);
    /* Fix connections with faces. */
    for_range(0u, num_faces(), [&](uint_t face_ind) {
        mesh_face_t& face = get_face(face_ind);
        std::for_each(face.begin_cell(), face.end_cell(), mesh_erase_func_t(cell_ind, replacement_cell_ind));
    });
    /* Deallocate the cell. */
    deallocate_cell_(cell_ind);
}   // mesh_t::erase_face

/** @internal
 ** Deallocate and erase element. */
template<typename mesh_elem_struct_t>
SKUNK_INLINE static void deallocate_element_(const mesh_elem_struct_t& elem, uint_t elem_ind,
                                             std::vector<byte_t>& elem_storage, std::vector<size_t>& elem_offsets) {
    /* Find ranges for element storage. */
    SKUNK_ASSERT(elem_ind <= elem_offsets.size());
    const size_t elem_beg = elem_offsets[elem_ind];
    const size_t elem_end = elem_offsets.size() > elem_ind + 1 ?
                            elem_offsets[elem_ind + 1] : elem_storage.size();

    /* Erase element offset and memory storage. */
    elem_offsets.erase(elem_offsets.begin() + elem_ind);
    elem_storage.erase(elem_storage.begin() + elem_beg, elem_storage.begin() + elem_end);
    /* Fix offsets for the following elements. */
    const size_t delta_elem = elem_end - elem_beg;
    for (; elem_ind < elem_offsets.size(); ++elem_ind) {
        SKUNK_ASSERT(elem_offsets[elem_ind] >= elem_beg);
        elem_offsets[elem_ind] -= delta_elem;
    }
}   // deallocate_element_
/** @internal
 ** Deallocate and erase the node. */
SKUNK_INLINE void mesh_t::deallocate_node_(uint_t node_ind) {
    deallocate_element_(get_node(node_ind), node_ind, m_node_storage, m_node_offsets);
}   // mesh_t::deallocate_node_
/** @internal
 ** Deallocate and erase the edge. */
SKUNK_INLINE void mesh_t::deallocate_edge_(uint_t edge_ind) {
    deallocate_element_(get_edge(edge_ind), edge_ind, m_edge_storage, m_edge_offsets);
}   // mesh_t::deallocate_edge_
/** @internal
 ** Deallocate and erase the face. */
SKUNK_INLINE void mesh_t::deallocate_face_(uint_t face_ind) {
    deallocate_element_(get_face(face_ind), face_ind, m_face_storage, m_face_offsets);
}   // mesh_t::deallocate_face_
/** @internal
 ** Deallocate and erase the cell. */
SKUNK_INLINE void mesh_t::deallocate_cell_(uint_t cell_ind) {
    deallocate_element_(get_cell(cell_ind), cell_ind, m_cell_storage, m_cell_offsets);
}   // mesh_t::deallocate_cell_

/** Erase edges that's length is less than specified tolerance. */
SKUNK_EXTERN void mesh_t::erase_degenerate_edges(real_t min_edge_length) {
    /* Build a list of the edges to erase
     * and remove the edges from faces connectivity. */
    std::set<uint_t> edges_to_erase;
    for (uint_t face_ind = 0; face_ind < num_faces(); ++face_ind) {
        mesh_face_t& face = get_face(face_ind);
        for (uint_t edge_loc = 0; edge_loc < face.num_edges(); ++edge_loc) {
            const uint_t edge_ind = face.begin_edge()[edge_loc];
            if (get_edge_length(edge_ind) <= min_edge_length) {
                edges_to_erase.insert(edge_ind);
                face.erase_edge(edge_loc);
            }
        }
    }
    /* Erase the edges in reverse order. */
    std::for_each(edges_to_erase.rbegin(), edges_to_erase.rend(), [&](uint_t edge_ind) {
        deallocate_edge_(edge_ind);
    });
}   // mesh_t::erase_degenerate_edges
/** Erase faces that's area is less than specified tolerance. */
SKUNK_EXTERN void mesh_t::erase_degenerate_faces(real_t min_face_area) {
    /* Build a list of the faces to erase
     * and remove the faces from cells connectivity. */
    std::set<uint_t> faces_to_erase;
    for (uint_t cell_ind = 0; cell_ind < num_cells(); ++cell_ind) {
        mesh_cell_t& cell = get_cell(cell_ind);
        for (uint_t face_loc = 0; face_loc < cell.num_faces(); ++face_loc) {
            const uint_t face_ind = cell.begin_face()[face_loc];
            const real_t face_area = get_face_area(face_ind);
            if (face_area <= min_face_area) {
                faces_to_erase.insert(face_ind);
                cell.erase_face(face_loc);
            }
        }
    }
    /* Erase the faces in reverse order. */
    std::for_each(faces_to_erase.rbegin(), faces_to_erase.rend(), [&](uint_t face_ind) {
        erase_face(face_ind);
    });
}   // mesh_t::erase_degenerate_faces

/**************************************************************************/
/**************************************************************************/

namespace {
/** @internal
 ** Functor that swaps connected element indices when two elements are reordered. */
class mesh_swap_func_t {
private:
    uint_t m_first_elem_ind, m_second_elem_ind;
public:
    SKUNK_INLINE mesh_swap_func_t(uint_t first_elem_ind, uint_t second_elem_ind)
        : m_first_elem_ind(first_elem_ind), m_second_elem_ind(second_elem_ind) {
    }
    SKUNK_INLINE void operator()(uint_t& elem_ind) const {
        if (elem_ind == m_first_elem_ind) {
            elem_ind = m_second_elem_ind;
        } else if (elem_ind == m_second_elem_ind) {
            elem_ind = m_first_elem_ind;
        }
    }
};  // class mesh_swap_func_t
}   // namespace

/** Change order of two nodes. */
SKUNK_EXTERN void mesh_t::reorder_nodes(uint_t first_node_ind, uint_t second_node_ind) {
    if (first_node_ind == second_node_ind) {
        return;
    }
    /* Replace connections with edges. */
    for_range(0u, num_edges(), [&](uint_t edge_ind) {
        mesh_edge_t& edge = get_edge(edge_ind);
        std::for_each(edge.begin_node(), edge.end_node(), mesh_swap_func_t(first_node_ind, second_node_ind));
    });
    /* Replace connections with faces. */
    for_range(0u, num_faces(), [&](uint_t face_ind) {
        mesh_face_t& face = get_face(face_ind);
        std::for_each(face.begin_node(), face.end_node(), mesh_swap_func_t(first_node_ind, second_node_ind));
    });
    /* Replace connections with cells. */
    for_range(0u, num_cells(), [&](uint_t cell_ind) {
        mesh_cell_t& cell = get_cell(cell_ind);
        std::for_each(cell.begin_node(), cell.end_node(), mesh_swap_func_t(first_node_ind, second_node_ind));
    });
    /* Swap node bytes. */
    swap_node_bytes_(first_node_ind, second_node_ind);
}   // mesh_t::reorder_nodes
/** Change order of two edges. */
SKUNK_EXTERN void mesh_t::reorder_edges(uint_t first_edge_ind, uint_t second_edge_ind) {
    if (first_edge_ind == second_edge_ind) {
        return;
    }
    /* Replace connections with faces. */
    for_range(0u, num_faces(), [&](uint_t face_ind) {
        mesh_face_t& face = get_face(face_ind);
        std::for_each(face.begin_edge(), face.end_edge(), 
                      mesh_swap_func_t(first_edge_ind, second_edge_ind));
    });
    /* Swap edge bytes. */
    swap_edge_bytes_(first_edge_ind, second_edge_ind);
}   // mesh_t::reorder_edges
/** Change order of two faces. */
SKUNK_EXTERN void mesh_t::reorder_faces(uint_t first_face_ind, uint_t second_face_ind) {
    if (first_face_ind == second_face_ind) {
        return;
    }
    /* Replace connections with cells. */
    for_range(0u, num_cells(), [&](uint_t cell_ind) {
        mesh_cell_t& cell = get_cell(cell_ind);
        std::for_each(cell.begin_face(), cell.end_face(), mesh_swap_func_t(first_face_ind, second_face_ind));
    });
    /* Swap face bytes. */
    swap_face_bytes_(first_face_ind, second_face_ind);
}   // mesh_t::reorder_faces
/** Change order of two cells. */
SKUNK_EXTERN void mesh_t::reorder_cells(uint_t first_cell_ind, uint_t second_cell_ind) {
    if (first_cell_ind == second_cell_ind) {
        return;
    }
    /* Replace connections with cells. */
    for_range(0u, num_faces(), [&](uint_t face_ind) {
        mesh_face_t& face = get_face(face_ind);
        std::for_each(face.begin_cell(), face.end_cell(), mesh_swap_func_t(first_cell_ind, second_cell_ind));
    });
    /* Swap cell bytes. */
    swap_cell_bytes_(first_cell_ind, second_cell_ind);
}   // mesh_t::reorder_cells

/** @internal
 ** Swap bytes of two elements. */
SKUNK_INLINE static void swap_elem_bytes_(uint_t first_elem_ind, uint_t second_elem_ind,
                                          std::vector<byte_t>& elem_storage, std::vector<size_t>& elem_offsets) {
    /* Find ranges for element storage. */
    SKUNK_ASSERT(first_elem_ind <= elem_offsets.size());
    const size_t first_elem_beg = elem_offsets[first_elem_ind];
    const size_t first_elem_end = elem_offsets.size() > first_elem_ind + 1 ?
                                  elem_offsets[first_elem_ind + 1] : elem_storage.size();

    /* Find ranges for new element storage. */
    SKUNK_ASSERT(second_elem_ind <= elem_offsets.size());
    const size_t second_elem_beg = elem_offsets[second_elem_ind];
    const size_t second_elem_end = elem_offsets.size() > second_elem_ind + 1 ?
                                   elem_offsets[second_elem_ind + 1] : elem_storage.size();

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
        for (; first_elem_ind <= second_elem_ind; ++first_elem_ind) {
            SKUNK_ASSERT(elem_offsets[first_elem_ind] >= first_elem_beg);
            elem_offsets[first_elem_ind] += range;
        }
    }
}   // swap_elem_bytes_
/** @internal
 ** Swap bytes of two nodes. */
SKUNK_INLINE void mesh_t::swap_node_bytes_(uint_t first_node_ind, uint_t second_node_ind) {
    swap_elem_bytes_(first_node_ind, second_node_ind, m_node_storage, m_node_offsets);
}   // mesh_t::swap_node_bytes_
/** @internal
 ** Swap bytes of two nodes. */
SKUNK_INLINE void mesh_t::swap_edge_bytes_(uint_t first_edge_ind, uint_t second_edge_ind) {
    swap_elem_bytes_(first_edge_ind, second_edge_ind, m_edge_storage, m_edge_offsets);
}   // mesh_t::swap_edge_bytes_
/** @internal
 ** Swap bytes of two nodes. */
SKUNK_INLINE void mesh_t::swap_face_bytes_(uint_t first_face_ind, uint_t second_face_ind) {
    swap_elem_bytes_(first_face_ind, second_face_ind, m_face_storage, m_face_offsets);
}   // mesh_t::swap_face_bytes_
/** @internal
 ** Swap bytes of two nodes. */
SKUNK_INLINE void mesh_t::swap_cell_bytes_(uint_t first_cell_ind, uint_t second_cell_ind) {
    swap_elem_bytes_(first_cell_ind, second_cell_ind, m_cell_storage, m_cell_offsets);
}   // mesh_t::swap_cell_bytes_

namespace {
/** @internal
 ** Functor that replaces connected element indices when elements are reordered. */
class mesh_reorder_func_t {
private:
    const std::vector<uint_t>& m_elem_inv_reordering;
public:
    SKUNK_INLINE explicit mesh_reorder_func_t(const std::vector<uint_t>& elem_reordering)
        : m_elem_inv_reordering(elem_reordering) {
    }
    SKUNK_INLINE void operator()(uint_t& elem_ind) const {
        if (elem_ind != npos) {
            elem_ind = m_elem_inv_reordering.at(static_cast<size_t>(elem_ind));
        }
    }
};  // class mesh_reorder_func_t
}   // namespace

/** Change order of all nodes. */
SKUNK_EXTERN void mesh_t::reorder_nodes(const std::vector<uint_t>& node_reordering) {
    /* Inverse reordering. */
    std::vector<uint_t> node_inv_reordering(node_reordering.size());
    inverse_reordering(node_reordering.begin(), node_reordering.end(), node_inv_reordering.begin());
    /* Replace connections with edges. */
    for_range(0u, num_edges(), [&](uint_t edge_ind) {
        mesh_edge_t& edge = get_edge(edge_ind);
        std::for_each(edge.begin_node(), edge.end_node(), mesh_reorder_func_t(node_inv_reordering));
    });
    /* Replace connections with faces. */
    for_range(0u, num_faces(), [&](uint_t face_ind) {
        mesh_face_t& face = get_face(face_ind);
        std::for_each(face.begin_node(), face.end_node(), mesh_reorder_func_t(node_inv_reordering));
    });
    /* Replace connections with cells. */
    for_range(0u, num_cells(), [&](uint_t cell_ind) {
        mesh_cell_t& cell = get_cell(cell_ind);
        std::for_each(cell.begin_node(), cell.end_node(), mesh_reorder_func_t(node_inv_reordering));
    });
    /* Reorder nodes bytes. */
    reorder_nodes_bytes_(node_reordering);
}   // mesh_t::reorder_nodes
/** Change order of all edges. */
SKUNK_EXTERN void mesh_t::reorder_edges(const std::vector<uint_t>& edge_reordering) {
    /* Inverse reordering. */
    std::vector<uint_t> edge_inv_reordering(edge_reordering.size());
    inverse_reordering(edge_reordering.begin(), edge_reordering.end(), edge_inv_reordering.begin());
    /* Replace connections with faces. */
    for_range(0u, num_faces(), [&](uint_t face_ind) {
        mesh_face_t& face = get_face(face_ind);
        std::for_each(face.begin_edge(), face.end_edge(), mesh_reorder_func_t(edge_inv_reordering));
    });
    /* Reorder edges bytes. */
    reorder_edges_bytes_(edge_reordering);
}   // mesh_t::reorder_edges
/** Change order of all faces. */
SKUNK_EXTERN void mesh_t::reorder_faces(const std::vector<uint_t>& face_reordering) {
    /* Inverse reordering. */
    std::vector<uint_t> face_inv_reordering(face_reordering.size());
    inverse_reordering(face_reordering.begin(), face_reordering.end(), face_inv_reordering.begin());
    /* Replace connections with cells. */
    for_range(0u, num_cells(), [&](uint_t cell_ind) {
        mesh_cell_t& cell = get_cell(cell_ind);
        std::for_each(cell.begin_face(), cell.end_face(), mesh_reorder_func_t(face_inv_reordering));
    });
    /* Reorder faces bytes. */
    reorder_faces_bytes_(face_reordering);
}   // mesh_t::reorder_faces
/** Change order of all cells. */
SKUNK_EXTERN void mesh_t::reorder_cells(const std::vector<uint_t>& cell_reordering) {
    /* Inverse reordering. */
    std::vector<uint_t> cell_inv_reordering(cell_reordering.size());
    inverse_reordering(cell_reordering.begin(), cell_reordering.end(), cell_inv_reordering.begin());
    /* Replace connections with faces. */
    for_range(0u, num_faces(), [&](uint_t face_ind) {
        mesh_face_t& face = get_face(face_ind);
        std::for_each(face.begin_cell(), face.end_cell(), mesh_reorder_func_t(cell_inv_reordering));
    });
    /* Reorder cells bytes. */
    reorder_cells_bytes_(cell_reordering);
}   // mesh_t::reorder_cells

/** @internal
 ** Swap bytes of two elements. */
SKUNK_INLINE static void reorder_elems_bytes_(std::vector<uint_t> elem_reordering,
                                              std::vector<byte_t>& elem_storage, std::vector<size_t>& elem_offsets) {
    reorder_swap(elem_reordering.begin(), elem_reordering.end(), [&](uint_t first_elem_ind, uint_t second_elem_ind) {
        swap_elem_bytes_(first_elem_ind, second_elem_ind, elem_storage, elem_offsets);
    });
}   // reorder_elems_bytes_
/** @internal
 ** Change order of bytes of all nodes. */
SKUNK_INLINE void mesh_t::reorder_nodes_bytes_(const std::vector<uint_t>& node_reordering) {
    reorder_elems_bytes_(node_reordering, m_node_storage, m_node_offsets);
}   // mesh_t::reorder_nodes_bytes_
/** @internal
 ** Change order of bytes of all edges. */
SKUNK_INLINE void mesh_t::reorder_edges_bytes_(const std::vector<uint_t>& edge_reordering) {
    reorder_elems_bytes_(edge_reordering, m_edge_storage, m_edge_offsets);
}   // mesh_t::reorder_edges_bytes_
/** @internal
 ** Change order of bytes of all faces. */
SKUNK_INLINE void mesh_t::reorder_faces_bytes_(const std::vector<uint_t>& face_reordering) {
    reorder_elems_bytes_(face_reordering, m_face_storage, m_face_offsets);
}   // mesh_t::reorder_faces_bytes_
/** @internal
 ** Change order of bytes of all cells. */
SKUNK_INLINE void mesh_t::reorder_cells_bytes_(const std::vector<uint_t>& cell_reordering) {
    reorder_elems_bytes_(cell_reordering, m_cell_storage, m_cell_offsets);
}   // mesh_t::reorder_cells_bytes_

/**************************************************************************/
/**************************************************************************/

SKUNK_EXTERN void mesh_t::reorder_faces() {
    /** @todo Refactor me! */
    {
        std::vector<uint_t> node_reordering(num_nodes());
        std::iota(node_reordering.begin(), node_reordering.end(), 0);
        std::stable_sort(node_reordering.begin(), node_reordering.end(), [&](uint_t node_ind_1, uint_t node_ind_2) {
            return get_node(node_ind_1).get_mark() < get_node(node_ind_2).get_mark();
        });
        reorder_nodes(node_reordering);
        for (uint_t node_ind = 0; node_ind < num_nodes(); ++node_ind) {
            mesh_node_t& node = get_node(node_ind);
            m_marked_node_ranges.resize(node.get_mark() + 2);
            m_marked_node_ranges[node.get_mark() + 1] += 1;
        }
        std::partial_sum(m_marked_node_ranges.begin(), m_marked_node_ranges.end(), m_marked_node_ranges.begin());
    }
    {
        std::vector<uint_t> edge_reordering(num_edges());
        std::iota(edge_reordering.begin(), edge_reordering.end(), 0);
        std::stable_sort(edge_reordering.begin(), edge_reordering.end(), [&](uint_t edge_ind_1, uint_t edge_ind_2) {
            return get_edge(edge_ind_1).get_mark() < get_edge(edge_ind_2).get_mark();
        });
        reorder_edges(edge_reordering);
        for (uint_t edge_ind = 0; edge_ind < num_edges(); ++edge_ind) {
            mesh_edge_t& edge = get_edge(edge_ind);
            m_marked_edge_ranges.resize(edge.get_mark() + 2);
            m_marked_edge_ranges[edge.get_mark() + 1] += 1;
        }
        std::partial_sum(m_marked_edge_ranges.begin(), m_marked_edge_ranges.end(), m_marked_edge_ranges.begin());
    }
    {
        std::vector<uint_t> face_reordering(num_faces());
        std::iota(face_reordering.begin(), face_reordering.end(), 0);
        std::stable_sort(face_reordering.begin(), face_reordering.end(), [&](uint_t face_ind_1, uint_t face_ind_2) {
            return get_face(face_ind_1).get_mark() < get_face(face_ind_2).get_mark();
        });
        reorder_faces(face_reordering);
        for (uint_t face_ind = 0; face_ind < num_faces(); ++face_ind) {
            mesh_face_t& face = get_face(face_ind);
            m_marked_face_ranges.resize(face.get_mark() + 2);
            m_marked_face_ranges[face.get_mark() + 1] += 1;
        }
        std::partial_sum(m_marked_face_ranges.begin(), m_marked_face_ranges.end(), m_marked_face_ranges.begin());
    }
    {
        std::vector<uint_t> cell_reordering(num_cells());
        std::iota(cell_reordering.begin(), cell_reordering.end(), 0);
        std::stable_sort(cell_reordering.begin(), cell_reordering.end(), [&](uint_t cell_ind_1, uint_t cell_ind_2) {
            return get_cell(cell_ind_1).get_mark() < get_cell(cell_ind_2).get_mark();
        });
        reorder_cells(cell_reordering);
        for (uint_t cell_ind = 0; cell_ind < num_cells(); ++cell_ind) {
            mesh_cell_t& cell = get_cell(cell_ind);
            m_marked_cell_ranges.resize(cell.get_mark() + 2);
            m_marked_cell_ranges[cell.get_mark() + 1] += 1;
        }
        std::partial_sum(m_marked_cell_ranges.begin(), m_marked_cell_ranges.end(), m_marked_cell_ranges.begin());
    }
}   // mesh_t::reorder_faces

/**************************************************************************/
/**************************************************************************/

/** Generate edges using the face to node connectivity.
 ** @warning This function may be slow and memory-consuming. */
SKUNK_EXTERN void mesh_t::generate_edges() {
    /* Face edge node index tables for various face types. */
    static const std::map<mesh_elem_type_t, std::vector<std::vector<uint_t>>> face_edge_nodes = {
        /* 1D faces. */
        { mesh_elem_type_t::node, {
            {0}
        } },
        /* 2D faces. */
        { mesh_elem_type_t::segment_2, {
            {0}, {1}
        } },
        { mesh_elem_type_t::segment_3, {
            {0}, {1}
        } },
        /* 3D faces. */
        { mesh_elem_type_t::triangle_3, {
            {1, 2},
            {2, 0}, {0, 1}
        } },
        { mesh_elem_type_t::triangle_6, {
            {1, 4, 2},
            {2, 5, 0}, {0, 3, 1}
        } },
        { mesh_elem_type_t::quadrangle_4, {
            {0, 1}, {1, 2},
            {2, 3}, {3, 0}
        } },
        { mesh_elem_type_t::quadrangle_8, {
            {0, 4, 1}, {1, 5, 2},
            {2, 6, 3}, {3, 7, 0}
        } },
        { mesh_elem_type_t::quadrangle_9, {
            {0, 4, 1}, {1, 5, 2},
            {2, 6, 3}, {3, 7, 0}
        } },
    };
    /* Build a map of the existing edges.
     * ( An edge is uniquely defined by a set of it's nodes. ) */
    std::map<std::set<uint_t>, uint_t> edge_cache;
    for (uint_t edge_ind = 0; edge_ind < num_edges(); ++edge_ind) {
        const mesh_edge_t& edge = get_edge(edge_ind);
        edge_cache.emplace(std::set<uint_t>(edge.begin_node(), edge.end_node()), edge_ind);
    }
    /* Generate missing edges. */
    for (uint_t face_ind = 0; face_ind < num_faces(); ++face_ind) {
        mesh_face_t& face = get_face(face_ind);
        for (uint_t edge_loc = 0; edge_loc < face.num_edges(); ++edge_loc) {
            uint_t& edge_ind = face.begin_edge()[edge_loc];
            if (edge_ind != npos) {
                continue;
            }
            /* Collect edge nodes using the table. */
            std::vector<uint_t> edge_nodes;
            for (uint_t node_loc : face_edge_nodes.at(face.get_type()).at(edge_loc)) {
                SKUNK_ASSERT(node_loc < face.num_nodes());
                const uint_t node_ind = face.begin_node()[node_loc];
                edge_nodes.push_back(node_ind);
            }
            /* Create the face or add current cell to the adjacency list. */
            std::set<uint_t> edge_nodes_set(edge_nodes.begin(), edge_nodes.end());
            const auto edge_iter = edge_cache.find(edge_nodes_set);
            if (edge_iter != edge_cache.cend()) {
                /* Locate the existing edge. */
                edge_ind = edge_iter->second;
            } else {
                /* Create a brand-new edge. */
                edge_ind = insert_edge(edge_nodes);
                edge_cache.emplace(edge_nodes_set, edge_ind);
                mesh_edge_t& edge = get_edge(edge_ind);
                edge.set_mark(face.get_mark());
            }
        }
    }
    /* Check faces:
     * each face should be connected to all edges. */
    for (uint_t face_ind = 0; face_ind < num_faces(); ++face_ind) {
        const mesh_face_t& face = get_face(face_ind);
        SKUNK_ASSERT(std::all_of(face.begin_edge(), face.end_edge(), is_not_npos));
    }
}   // mesh_t::generate_edges
/** Generate faces using the cell to node connectivity.
 ** @warning This function may be slow and memory-consuming. */
SKUNK_EXTERN void mesh_t::generate_faces() {
    /* Cell face node index tables for various cell types. */
    static const std::map<mesh_elem_type_t, std::vector<std::vector<uint_t>>> cell_face_nodes = {
        /* 1D cells. */
        { mesh_elem_type_t::segment_2, {
            {0}, {1},
        } },
        { mesh_elem_type_t::segment_3, {
            {0}, {1},
        } },
        /* 2D cells. */
        { mesh_elem_type_t::triangle_3, {
            {1, 2},
            {2, 0}, {0, 1},
        } },
        { mesh_elem_type_t::triangle_6, {
            {1, 4, 2},
            {2, 5, 0}, {0, 3, 1},
        } },
        { mesh_elem_type_t::quadrangle_4, {
            {0, 1}, {1, 2},
            {2, 3}, {3, 0},
        } },
        { mesh_elem_type_t::quadrangle_8, {
            {0, 4, 1}, {1, 5, 2},
            {2, 6, 3}, {3, 7, 0},
        } },
        { mesh_elem_type_t::quadrangle_9, {
            {0, 4, 1}, {1, 5, 2},
            {2, 6, 3}, {3, 7, 0},
        } },
        /* 3D cells.
         * @todo Add other cell types! */
        { mesh_elem_type_t::tetrahedron_4, {
            {0, 2, 1}, {0, 1, 3},
            {1, 2, 3}, {2, 0, 3}
        } },
        { mesh_elem_type_t::hexahedron_8, {
            {0, 4, 7, 3}, {1, 2, 6, 5},
            {3, 7, 6, 2}, {0, 1, 5, 4},
            {0, 3, 2, 1}, {4, 5, 6, 7},
        } },
        { mesh_elem_type_t::hexahedron_20, {
            {0, 4, 7, 3, 12, 19, 15, 11}, {1, 2, 6, 5,  9, 14, 17, 13},
            {3, 7, 6, 2, 15, 18, 14, 10}, {0, 1, 5, 4,  8, 13, 16, 12},
            {0, 3, 2, 1, 11, 10,  9,  8}, {4, 5, 6, 7, 16, 17, 18, 19},
        } },
        { mesh_elem_type_t::hexahedron_27, {
            {0, 4, 7, 3, 12, 19, 15, 11, 23}, {1, 2, 6, 5,  9, 14, 17, 13, 21},
            {3, 7, 6, 2, 15, 18, 14, 10, 22}, {0, 1, 5, 4,  8, 13, 16, 12, 24},
            {0, 3, 2, 1, 11, 10,  9,  8, 20}, {4, 5, 6, 7, 16, 17, 18, 19, 25},
        } },
    };
    /* Build a map of the existing faces.
     * ( A face is uniquely defined by a set of it's nodes. ) */
    std::map<std::set<uint_t>, uint_t> face_cache;
    for (uint_t face_ind = 0; face_ind < num_faces(); ++face_ind) {
        const mesh_face_t& face = get_face(face_ind);
        face_cache.emplace(std::set<uint_t>(face.begin_node(), face.end_node()), face_ind);
    }
    /* Generate missing faces. */
    for (uint_t cell_ind = 0; cell_ind < num_cells(); ++cell_ind) {
        mesh_cell_t& cell = get_cell(cell_ind);
        for (uint_t face_loc = 0; face_loc < cell.num_faces(); ++face_loc) {
            uint_t& face_ind = cell.begin_face()[face_loc];
            if (face_ind != npos) {
                continue;
            }
            /* Collect face nodes using the table. */
            std::vector<uint_t> face_nodes;
            for (uint_t node_loc : cell_face_nodes.at(cell.get_type()).at(face_loc)) {
                SKUNK_ASSERT(node_loc < cell.num_nodes());
                const uint_t node_ind = cell.begin_node()[node_loc];
                face_nodes.push_back(node_ind);
            }
            /* Create the face or add current cell to the adjacency list. */
            const std::set<uint_t> face_nodes_set(face_nodes.begin(), face_nodes.end());
            const auto face_iter = face_cache.find(face_nodes_set);
            if (face_iter != face_cache.cend()) {
                /* Locate the existing face.
                 * ( We should be careful with face orientation here. ) */
                face_ind = face_iter->second;
                mesh_face_t& face = get_face(face_ind);
                if (std::equal(face.begin_node(), face.end_node(), face_nodes.cbegin())) {
                    SKUNK_ASSERT(face.get_inner_cell() == npos);
                    face.get_inner_cell() = cell_ind;
                } else {
                    SKUNK_ASSERT(face.get_outer_cell() == npos);
                    face.get_outer_cell() = cell_ind;
                }
            } else {
                /* Create a brand-new face.
                 * Assign the current cell as the inner one. */
                face_ind = insert_face(face_nodes);
                face_cache.emplace(face_nodes_set, face_ind);
                mesh_face_t& face = get_face(face_ind);
                face.get_inner_cell() = cell_ind;
            }
        }
    }
    /* Check cells:
     * each cell should be connected to all faces. */
    for (uint_t cell_ind = 0; cell_ind < num_cells(); ++cell_ind) {
        const mesh_cell_t& cell = get_cell(cell_ind);
        SKUNK_ASSERT(std::all_of(cell.begin_face(), cell.end_face(), is_not_npos));
    }
    /* Check faces:
     * each internal face should be connected to two cells;
     * each boundary face should be connected to at least one cell. */
    for (uint_t face_ind = 0; face_ind < num_faces(); ++face_ind) {
        const mesh_face_t& face = get_face(face_ind);
        if (face.get_mark() == 0) {
            SKUNK_ASSERT(std::all_of(face.begin_cell(), face.end_cell(), is_not_npos));
        } else {
            SKUNK_ASSERT(std::any_of(face.begin_cell(), face.end_cell(), is_not_npos));
        }
    }
}   // mesh_t::generate_faces

/** Generate boundary cells to complete face connectivity. */
SKUNK_EXTERN void mesh_t::generate_boundary_cells() {
    /* A node and edge flip table for various face types. */
    static const std::map<mesh_elem_type_t, std::pair<std::vector<uint_t>,
                                                      std::vector<uint_t>>> face_nodes_edges_flip_table {
        /* 1D faces. */
        { mesh_elem_type_t::node, {
          {0}, {0}
        } },
        /* 2D faces. */
        { mesh_elem_type_t::segment_2, {
          {1, 0}, {1, 0}
        } },
        { mesh_elem_type_t::segment_3, {
          {1, 0, 2}, {1, 0}
        } },
        /* 3D faces. */
        { mesh_elem_type_t::triangle_3, {
          {0, 2, 1}, {0, 2, 1}
        } },
        { mesh_elem_type_t::triangle_6, {
          {0, 2, 1, 5, 4, 3}, {0, 2, 1}
        } },
        { mesh_elem_type_t::quadrangle_4, {
          {0, 3, 2, 1}, {0, 3, 2, 1}
        } },
        { mesh_elem_type_t::quadrangle_8, {
          {0, 3, 2, 1, 7, 6, 5, 4}, {0, 3, 2, 1}
        } },
        { mesh_elem_type_t::quadrangle_9, {
          {0, 3, 2, 1, 7, 6, 5, 4, 9}, {0, 3, 2, 1}
        } },
    };
    /* Generate boundary cells and fix boundary faces orientation. */
    for (uint_t face_ind = 0; face_ind < num_faces(); ++face_ind) {
        mesh_face_t& face = get_face(face_ind);
        if (face.get_mark() == 0) {
            continue;
        }
        /* Boundary faces should be oriented outwards from the mesh. */
        if (face.get_inner_cell() == npos) {
            /* Flip normal and cell connectivity. */
            face.set_normal(-face.get_normal());
            std::swap(face.get_inner_cell(), face.get_outer_cell());
            /* Flip node and edge connectivity. */
            std::vector<uint_t> node_reordering, edge_reordering;
            std::tie(node_reordering, edge_reordering) = face_nodes_edges_flip_table.at(face.get_type());
            reorder(node_reordering.begin(), node_reordering.end(), face.begin_node());
            reorder(edge_reordering.begin(), edge_reordering.end(), face.begin_edge());
        }
        /* Generate the boundary cell: reflect a connected interior cell. */
        const mesh_cell_t& cell = get_cell(face.get_inner_cell());
        uint_t boundary_cell_dim;
        std::vector<uint_t> boundary_cell_nodes;
#if 0
        boundary_cell_dim = m_dim - 1;
        boundary_cell_nodes.assign(face.begin_node(), face.end_node());
#endif
#if 1
        boundary_cell_dim = m_dim;
        std::for_each(cell.begin_node(), cell.end_node(), [&](uint_t node_ind) {
            const auto node_iter = std::find(face.begin_node(), face.end_node(), node_ind);
            if (node_iter == face.end_node()) {
                /* Reflect an interior cell node. */
                vec3_t node_pos = get_node_position(node_ind);
                const vec3_t delta = node_pos - face.get_center_position();
                node_pos -= 2.0*(real_t)delta.inner_prod(face.get_normal())*face.get_normal();
                const mesh_node1_t node(node_pos, face.get_mark());
                boundary_cell_nodes.push_back(insert_node(node));
            } else {
                /* Insert a boundary face node. */
                boundary_cell_nodes.push_back(node_ind);
            }
        });
#endif
        /* Insert the boundary cell. */
        const uint_t boundary_cell_ind = insert_cell(boundary_cell_nodes, face.get_mark(), boundary_cell_dim);
        mesh_cell_t& boundary_cell = get_cell(boundary_cell_ind);
        while (boundary_cell.num_faces() != 1) {
            boundary_cell.erase_face(0);
        }
        boundary_cell.begin_face()[0] = face_ind;
        face.get_outer_cell() = boundary_cell_ind;
    }
}   // mesh_t::generate_boundary_cells

}   // namespace skunk

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
