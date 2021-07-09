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

#include <libSkunkMesh/SkunkMeshStructured.hh>

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace skunk {

inline
vec3_t cyl(const vec3_t& sph) {
    const auto r = sph.x, phi = sph.y, z = sph.z;
    return { r*std::cos(phi), r*std::sin(phi), z };
}

inline
vec3_t sph(const vec3_t& sph) {
    const auto r = sph.x, phi = sph.y, psi = sph.z;
    vec3_t s = { r*std::cos(phi)*std::cos(psi), r*std::sin(phi)*std::cos(psi), r*std::sin(psi) };
    auto e = 1e-10;
    if (s.x <= 1e-3 || s.y <= 1e-3) {
        s.x = std::copysign(e*std::round(std::abs(s.x)/e), s.x);
        s.y = std::copysign(e*std::round(std::abs(s.y)/e), s.y);
        s.z = std::copysign(e*std::round(std::abs(s.z)/e), s.z);
    }
    return s;
}

inline
vec3_t permute(const vec3_t& r) {
    //real_t s = 1000.0*r.x*r.y + 10.0*r.z;
    real_t s = c_pi/RAND_MAX*std::rand();
    if (r.x == 0.0 || r.x == 10.0) return r;
    if (r.y == 0.0 || r.y == 2.50) return r;
    const real_t e = 1e-2;
    return {
        r.x + e*std::cos(s),
        r.y + e*std::sin(s),
        r.z + 0*e*std::sin(s)
    };
}

/** 3D index. */
class index3_t {
private:
    int_t m_x_begin, m_x_end,
          m_y_begin, m_y_end,
          m_z_begin, m_z_end, m_offset;
public:
    /** Construct a 3D index. */
    index3_t(int_t x_begin, int_t x_end,
             int_t y_begin, int_t y_end,
             int_t z_begin, int_t z_end, int_t offset = 0)
        : m_x_begin(x_begin), m_x_end(x_end),
          m_y_begin(y_begin), m_y_end(y_end),
          m_z_begin(z_begin), m_z_end(z_end), m_offset(offset) {
    }

    /** Total size of elements */
    size_t size() const {
        return (m_x_end-m_x_begin)*(m_y_end-m_y_begin)*(m_z_end-m_z_begin);
    }
    uint_t operator()(int_t ix, int_t iy, int_t iz) const {
        SKUNK_ASSERT(m_x_begin <= ix && ix < m_x_end);
        SKUNK_ASSERT(m_y_begin <= iy && iy < m_y_end);
        SKUNK_ASSERT(m_z_begin <= iz && iz < m_z_end);
        return ix + (m_x_end-m_x_begin)*(iy + iz*(m_y_end-m_y_begin)) + m_offset;
    }
};  // class MhdIndex3D

/** Construct a uniform structured 2D segment mesh. */
SKUNK_EXTERN structured_mesh_t::structured_mesh_t(real_t x0, real_t x1, uint_t nx,
                                                  uint_t x0_mark, uint_t x1_mark)
    : mesh_t(1) {
}   // structured_mesh_t::structured_mesh_t

/** Construct a uniform structured 2D rectangular mesh. */
SKUNK_EXTERN structured_mesh_t::structured_mesh_t(real_t x0, real_t x1, uint_t nx,
                                                  real_t y0, real_t y1, uint_t ny,
                                                  uint_t x0_mark, uint_t x1_mark,
                                                  uint_t y0_mark, uint_t y1_mark,
                                                  real_t skew_angle, bool triangulate)
    : mesh_t(2) {
}   // structured_mesh_t::structured_mesh_t

/** Construct a uniform structured 3D hexahedral mesh. */
SKUNK_EXTERN structured_mesh_t::structured_mesh_t(real_t x0, real_t x1, uint_t nx,
                                                  real_t y0, real_t y1, uint_t ny,
                                                  real_t z0, real_t z1, uint_t nz,
                                                  uint_t x0_mark, uint_t x1_mark,
                                                  uint_t y0_mark, uint_t y1_mark,
                                                  uint_t z0_mark, uint_t z1_mark)
    : mesh_t(3) {
    /* Uniform steps. */
    const real_t hx = (x1 - x0)/nx;
    const real_t hy = (y1 - y0)/ny;
    const real_t hz = (z1 - z0)/nz;

    /* Construct nodes. */
    const index3_t node_xyz(0, nx+1, 0, ny+1, 0, nz+1);
    for (uint_t iz = 0; iz <= nz; ++iz) {
        for (uint_t iy = 0; iy <= ny; ++iy) {
            for (uint_t ix = 0; ix <= nx; ++ix) {
                /* Create a node. */
                mesh_node1_t node({x0 + ix*hx, y0 + iy*hy, z0 + iz*hz });
                /* Insert the node. */
                const uint_t node_ind = insert_node(node);
                SKUNK_ASSERT(node_ind == node_xyz(ix, iy, iz));
            }
        }
    }

    /* Construct edges. */
    const index3_t edge_x_xyz(0, nx, 0, ny+1, 0, nz+1);
    for (uint_t iz = 0; iz <= nz; ++iz) {
        for (uint_t iy = 0; iy <= ny; ++iy) {
            for (uint_t ix = 0; ix < nx; ++ix) {
                /* Create an edge. */
                mesh_edge_segment2_t edge;
                edge.set_nodes({ node_xyz(ix  , iy  , iz  ),
                                 node_xyz(ix+1, iy  , iz  ) });
                /* Set boundary marks. */
                if (iy == 0) {
                    edge.set_mark(y0_mark);
                } else if (iy == ny) {
                    edge.set_mark(y1_mark);
                } else if (iz == 0) {
                    edge.set_mark(z0_mark);
                } else if (iz == nz) {
                    edge.set_mark(z1_mark);
                }
                /* Insert the edge. */
                const uint_t edge_ind = insert_edge(edge);
                SKUNK_ASSERT(edge_ind == edge_x_xyz(ix, iy, iz));
            }
        }
    }
    const index3_t edge_y_xyz(0, nx+1, 0, ny, 0, nz+1, edge_x_xyz.size());
    for (uint_t iz = 0; iz <= nz; ++iz) {
        for (uint_t iy = 0; iy < ny; ++iy) {
            for (uint_t ix = 0; ix <= nx; ++ix) {
                /* Create an edge. */
                mesh_edge_segment2_t edge;
                edge.set_nodes({ node_xyz(ix  , iy  , iz  ),
                                 node_xyz(ix  , iy+1, iz  ) });
                /* Set boundary marks. */
                if (ix == 0) {
                    edge.set_mark(x0_mark);
                } else if (ix == nx) {
                    edge.set_mark(x1_mark);
                } else if (iz == 0) {
                    edge.set_mark(z0_mark);
                } else if (iz == nz) {
                    edge.set_mark(z1_mark);
                }
                /* Insert the edge. */
                const uint_t edge_ind = insert_edge(edge);
                SKUNK_ASSERT(edge_ind == edge_y_xyz(ix, iy, iz));
            }
        }
    }
    const index3_t edge_z_xyz(0, nx+1, 0, ny+1, 0, nz, edge_x_xyz.size() + edge_y_xyz.size());
    for (uint_t iz = 0; iz < nz; ++iz) {
        for (uint_t iy = 0; iy <= ny; ++iy) {
            for (uint_t ix = 0; ix <= nx; ++ix) {
                /* Create an edge. */
                mesh_edge_segment2_t edge;
                edge.set_nodes({ node_xyz(ix  , iy  , iz  ),
                                 node_xyz(ix  , iy  , iz+1) });
                /* Set boundary marks. */
                if (ix == 0) {
                    edge.set_mark(x0_mark);
                } else if (ix == nx) {
                    edge.set_mark(x1_mark);
                } else if (iy == 0) {
                    edge.set_mark(y0_mark);
                } else if (iy == ny) {
                    edge.set_mark(y1_mark);
                }
                /* Insert the edge. */
                const uint_t edge_ind = insert_edge(edge);
                SKUNK_ASSERT(edge_ind == edge_z_xyz(ix, iy, iz));
            }
        }
    }

    /* Construct faces. */
    const index3_t face_yz_xyz(0, nx+1, 0, ny, 0, nz);
    for (uint_t iz = 0; iz < nz; ++iz) {
        for (uint_t iy = 0; iy < ny; ++iy) {
            for (uint_t ix = 0; ix <= nx; ++ix) {
                /* Create a face. */
                mesh_face_quadrangle4_t face;
                face.set_nodes({ node_xyz(ix  , iy  , iz  ),
                                 node_xyz(ix  , iy  , iz+1),
                                 node_xyz(ix  , iy+1, iz+1),
                                 node_xyz(ix  , iy+1, iz  ), });
                face.set_edges({ edge_z_xyz(ix  , iy  , iz  ),
                                 edge_y_xyz(ix  , iy  , iz+1),
                                 edge_z_xyz(ix  , iy+1, iz  ),
                                 edge_y_xyz(ix  , iy  , iz  ), });
                /* Set marks. */
                if (ix == 0) {
                    face.set_mark(x0_mark);
                } else if (ix == nx) {
                    face.set_mark(x1_mark);
                }
                /* Insert the face. */
                const uint_t face_ind = insert_face(face);
                SKUNK_ASSERT(face_ind == face_yz_xyz(ix, iy, iz));
            }
        }
    }
    const index3_t face_xz_xyz(0, nx, 0, ny+1, 0, nz, face_yz_xyz.size());
    for (uint_t iz = 0; iz < nz; ++iz) {
        for (uint_t iy = 0; iy <= ny; ++iy) {
            for (uint_t ix = 0; ix < nx; ++ix) {
                /* Create a face. */
                mesh_face_quadrangle4_t face;
                face.set_nodes({ node_xyz(ix  , iy  , iz  ),
                                 node_xyz(ix+1, iy  , iz  ),
                                 node_xyz(ix+1, iy  , iz+1),
                                 node_xyz(ix  , iy  , iz+1), });
                face.set_edges({ edge_x_xyz(ix  , iy  , iz  ),
                                 edge_z_xyz(ix+1, iy  , iz  ),
                                 edge_x_xyz(ix  , iy  , iz+1),
                                 edge_z_xyz(ix  , iy  , iz  ), });
                /* Set marks. */
                if (iy == 0) {
                    face.set_mark(y0_mark);
                } else if (iy == ny) {
                    face.set_mark(y1_mark);
                }
                /* Insert the face. */
                const uint_t face_ind = insert_face(face);
                SKUNK_ASSERT(face_ind == face_xz_xyz(ix, iy, iz));
            }
        }
    }
    const index3_t face_xy_xyz(0, nx, 0, ny, 0, nz+1, face_yz_xyz.size() + face_xz_xyz.size());
    for (uint_t iz = 0; iz <= nz; ++iz) {
        for (uint_t iy = 0; iy < ny; ++iy) {
            for (uint_t ix = 0; ix < nx; ++ix) {
                /* Create a face. */
                mesh_face_quadrangle4_t face;
                face.set_nodes({ node_xyz(ix  , iy  , iz  ),
                                 node_xyz(ix  , iy+1, iz  ),
                                 node_xyz(ix+1, iy+1, iz  ),
                                 node_xyz(ix+1, iy  , iz  ), });
                face.set_edges({ edge_x_xyz(ix  , iy  , iz  ),
                                 edge_y_xyz(ix  , iy  , iz  ),
                                 edge_x_xyz(ix  , iy+1, iz  ),
                                 edge_y_xyz(ix+1, iy  , iz  ), });
                /* Set marks. */
                if (iz == 0) {
                    face.set_mark(z0_mark);
                } else if (iz == nz) {
                    face.set_mark(z1_mark);
                }
                /* Insert the face. */
                const uint_t face_ind = insert_face(face);
                SKUNK_ASSERT(face_ind == face_xy_xyz(ix, iy, iz));
            }
        }
    }

    /* Construct cells. */
    const index3_t cell_xyz(0, nx, 0, ny, 0, nz);
    for (uint_t iz = 0; iz < nz; ++iz) {
        for (uint_t iy = 0; iy < ny; ++iy) {
            for (uint_t ix = 0; ix < nx; ++ix) {
                /* Create a cell. */
                mesh_cell_hexahedron8_t cell;
                cell.set_nodes({ node_xyz(ix  , iy  , iz  ),
                                 node_xyz(ix+1, iy  , iz  ),
                                 node_xyz(ix+1, iy+1, iz  ),
                                 node_xyz(ix  , iy+1, iz  ),
                                 node_xyz(ix  , iy  , iz+1),
                                 node_xyz(ix+1, iy  , iz+1),
                                 node_xyz(ix+1, iy+1, iz+1),
                                 node_xyz(ix  , iy+1, iz+1), });
                cell.set_faces({ face_yz_xyz(ix  , iy  , iz  ),
                                 face_yz_xyz(ix+1, iy  , iz  ),
                                 face_xz_xyz(ix  , iy  , iz  ),
                                 face_xz_xyz(ix  , iy+1, iz  ),
                                 face_xy_xyz(ix  , iy  , iz  ),
                                 face_xy_xyz(ix  , iy  , iz+1), });
                /* Insert the cell. */
                const uint_t cell_ind = insert_cell(cell);
                SKUNK_ASSERT(cell_ind == cell_xyz(ix, iy, iz));
                /* Connect faces with a cell. */
                for (uint_t face_loc = 0; face_loc < cell.num_faces(); face_loc += 2) {
                    get_face(cell.begin_face()[face_loc]).get_inner_cell() = cell_ind;
                }
                for (uint_t face_loc = 1; face_loc < cell.num_faces(); face_loc += 2) {
                    get_face(cell.begin_face()[face_loc]).get_outer_cell() = cell_ind;
                }
            }
        }
    }

    //const real_t mrl = get_min_edge_length();
    //const real_t mfa = get_min_face_area();
    //erase_degenerate_edges((1.0+1.0e-2)*mrl);
    //erase_degenerate_faces(1e-10);
    generate_boundary_cells();
    //generate_ghost_cells_and_fix_boundaries_();
    reorder_faces();
}   // structured_mesh_t::structured_mesh_t

/**************************************************************************/
/**************************************************************************/

}   // namespace skunk

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //