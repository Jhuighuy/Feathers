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
#include <libFeathersMesh/Element.hh>
#include <libFeathersUtils/Parallel.hh>
#include <libSkunkMisc/SkunkMiscReordering.hh>

#include <set>
#include <map>
#include <fstream>


// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace feathers {

bool cMesh::read_triangle(const char *path) {
    std::string line;

    std::ifstream node_file(path + std::string("node"));
    FEATHERS_ENSURE(node_file.is_open());
    uint_t num_nodes = 0;
    node_file >> num_nodes >> m_dim;
    std::getline(node_file, line);
    for (uint_t i = 0; i < num_nodes; ++i) {
        uint_t node_index = 0;
        vec3_t node_coords(0.0);
        node_file >> node_index >> node_coords.x >> node_coords.y;
        std::getline(node_file, line);
        FEATHERS_ENSURE(node_index == insert_node(node_coords));
    }

    std::ifstream face_file(path + std::string("edge"));
    FEATHERS_ENSURE(face_file.is_open());
    uint_t num_faces = 0;
    face_file >> num_faces;
    std::getline(face_file, line);
    for (uint_t i = 0; i < num_faces; ++i) {
        uint_t face_index = 0;
        std::vector<uint_t> face_nodes(2, npos);
        uint_t mark = 0;
        face_file >> face_index >> face_nodes[0] >> face_nodes[1] >> mark;
        FEATHERS_ENSURE(
            face_index == insert_face(eShape::segment_2, face_nodes, mark));
        std::getline(face_file, line);
    }

    std::ifstream cell_file(path + std::string("ele"));
    FEATHERS_ENSURE(cell_file.is_open());
    uint_t num_cells = 0;
    cell_file >> num_cells;
    std::getline(cell_file, line);
    for (uint_t i = 0; i < num_cells; ++i) {
        uint_t cell_index = 0;
        std::vector<uint_t> cell_nodes(3, npos);
        cell_file >> cell_index >> cell_nodes[0] >> cell_nodes[1] >> cell_nodes[2];
        FEATHERS_ENSURE(
            cell_index == insert_cell(eShape::triangle_3, cell_nodes));
        std::getline(cell_file, line);
    }

    generate_faces();
    generate_boundary_cells();
    reorder_faces();
    compute_all_shape_properties();
    return true;
}   // сMesh::read_triangle

bool cMesh::read_tetgen(const char *path) {
    std::string line;

    std::ifstream node_file(path + std::string("node"));
    FEATHERS_ENSURE(node_file.is_open());
    uint_t num_nodes = 0;
    node_file >> num_nodes >> m_dim;
    std::getline(node_file, line);
    for (uint_t i = 0; i < num_nodes; ++i) {
        uint_t node_index = 0;
        vec3_t node_coords(0.0);
        node_file >> node_index >> node_coords.x >> node_coords.y >> node_coords.z;
        FEATHERS_ENSURE(
            node_index == insert_node(node_coords));
        std::getline(node_file, line);
    }

    std::ifstream face_file(path + std::string("face"));
    FEATHERS_ENSURE(face_file.is_open());
    uint_t num_faces = 0;
    face_file >> num_faces;
    std::getline(face_file, line);
    for (uint_t i = 0; i < num_faces; ++i) {
        uint_t face_index = 0;
        std::vector<uint_t> face_nodes(3, npos);
        uint_t mark = 0;
        face_file >> face_index >> face_nodes[0] >> face_nodes[1] >> face_nodes[2] >> mark;
        FEATHERS_ENSURE(
            face_index == insert_face(eShape::triangle_3, face_nodes, mark));
        std::getline(face_file, line);
    }

    std::ifstream cell_file(path + std::string("ele"));
    FEATHERS_ENSURE(cell_file.is_open());
    uint_t num_cells = 0;
    cell_file >> num_cells;
    std::getline(cell_file, line);
    for (uint_t i = 0; i < num_cells; ++i) {
        uint_t cell_index = 0;
        std::vector<uint_t> cell_nodes(4, npos);
        cell_file >> cell_index >> cell_nodes[0] >> cell_nodes[1] >> cell_nodes[2] >> cell_nodes[3];
        FEATHERS_ENSURE(
            cell_index == insert_cell(eShape::tetrahedron_4, cell_nodes));
        std::getline(cell_file, line);
    }

    generate_faces();
    generate_boundary_cells();
    reorder_faces();
    compute_all_shape_properties();
    return true;
}   // сMesh::read_tetgen

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Compute edge shape properties.
 */
void cMesh::compute_edge_shape_properties() {
    std::tie(m_min_edge_length, m_max_edge_length) =
        for_range_minmax(begin_edge(*this), end_edge(*this),
                         +huge, -huge, [&](tEdgeMutableIter edge) {
            std::unique_ptr<const iElement> edge_element = edge.get_element_object();
            edge.set_length(edge_element->get_length_or_area_or_volume());
            edge.set_direction(edge_element->get_direction());
            return edge.get_length();
        });
}   // сMesh::compute_edge_shape_properties

/**
 * Compute face shape properties.
 */
void cMesh::compute_face_shape_properties() {
    std::tie(m_min_face_area, m_max_face_area) =
        for_range_minmax(begin_face(*this), end_face(*this),
                         +huge, -huge, [&](tFaceMutableIter face) {
            std::unique_ptr<const iElement> face_element = face.get_element_object();
            face.set_area(face_element->get_length_or_area_or_volume());
            face.set_normal(face_element->get_normal());
            face.set_center_coords(face_element->get_center_coords());
            return face.get_area();
        });
}   // сMesh::compute_face_shape_properties

/**
 * Compute cell shape properties.
 */
void cMesh::compute_cell_shape_properties() {
    std::tie(m_min_face_area, m_max_face_area) =
        for_range_minmax(begin_cell(*this), end_cell(*this),
                         +huge, -huge, [&](tCellMutableIter cell) {
            std::unique_ptr<const iElement> cell_element = cell.get_element_object();
            cell.set_volume(cell_element->get_length_or_area_or_volume());
            cell.set_center_coords(cell_element->get_center_coords());
            return cell.get_volume();
        });
}   // сMesh::compute_cell_shape_properties

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Insert a new node into the mesh.
 * @returns Index of the inserted node.
 */
uint_t cMesh::insert_node(const vec3_t& node_coords, uint_t mark) {
    const uint_t node_index = m_num_nodes;
    m_num_nodes += 1;
    m_node_marks.emplace_back(mark);
    m_node_coords.emplace_back(node_coords);
    m_node_nodes.emplace_back_row();
    m_node_edges.emplace_back_row();
    m_node_faces.emplace_back_row();
    m_node_cells.emplace_back_row();
    return node_index;
}   // сMesh::insert_node

/**
 * Insert a new edge into the mesh.
 * @returns Index of the inserted edge.
 */
uint_t cMesh::insert_edge(const std::unique_ptr<const iElement>& edge, uint_t mark) {
    const uint_t edge_index = m_num_edges;
    m_num_edges += 1;
    m_edge_marks.emplace_back(mark);
    m_edge_shapes.emplace_back(edge->get_shape());
    m_edge_lengths.emplace_back(edge->get_length_or_area_or_volume());
    m_edge_directions.emplace_back(edge->get_direction());
    m_edge_nodes.emplace_back_row(edge->get_nodes().begin(), edge->get_nodes().end());
    m_edge_edges.emplace_back_row();
    m_edge_faces.emplace_back_row();
    m_edge_cells.emplace_back_row();
    FEATHERS_NOT_IMPLEMENTED();
    return edge_index;
}   // сMesh::insert_edge

/**
 * Insert a new face into the mesh.
 * @returns Index of the inserted face.
 */
uint_t cMesh::insert_face(eShape face_shape, const std::vector<uint_t>& face_nodes, uint_t mark) {
    const uint_t face_index = m_num_faces;
    m_num_faces += 1;
    m_face_marks.emplace_back(mark);
    m_face_shapes.emplace_back(face_shape);
    m_face_areas.emplace_back();
    m_face_normals.emplace_back();
    m_face_center_coords.emplace_back();
    m_face_nodes.emplace_back_row(face_nodes.begin(), face_nodes.end());
    m_face_edges.emplace_back_row();
    m_face_faces.emplace_back_row();
    m_face_cells.emplace_back_row(2); // TODO
    return face_index;
}   // сMesh::insert_face

/**
 * Insert a new cell into the mesh.
 * @returns Index of the inserted cell.
 */
uint_t cMesh::insert_cell(eShape cell_shape, const std::vector<uint_t>& cell_nodes, uint_t mark) {
    const uint_t cell_index = m_num_cells;
    m_num_cells += 1;
    m_cell_marks.emplace_back(mark);
    m_cell_shapes.emplace_back(cell_shape);
    m_cell_volumes.emplace_back();
    m_cell_center_coords.emplace_back();
    m_cell_nodes.emplace_back_row(cell_nodes.begin(), cell_nodes.end());
    m_cell_edges.emplace_back_row();
    m_cell_faces.emplace_back_row(3); // TODO
    m_cell_cells.emplace_back_row();
    return cell_index;
}   // сMesh::insert_cell

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

class tReorderFunc {
private:
    const std::vector<uint_t>& m_indices;
public:
    explicit tReorderFunc(const std::vector<uint_t>& indices)
        : m_indices(indices) {
    }
    uint_t operator()(uint_t index) const {
        return index != npos ? m_indices.at(index) : npos;
    }
};  // class tReorderFunc

template<typename tTag>
void cMesh::fix_permutation_and_adjacency_(tTag tag, std::vector<uint_t>& permutation) {
    /* Fix permutations by resorting it by marks.
     * Stable sort is used here to preserve the permutation in the best possible way. */
    std::stable_sort(permutation.begin(), permutation.end(), [&](uint_t index_1, uint_t index_2) {
        return get_mark(tag, index_1) < get_mark(tag, index_2);
    });
    /* Fix adjacency lists. */
    std::vector<uint_t> ordering(permutation.size());
    convert_permutation_to_ordering(
        permutation.begin(), permutation.end(), ordering.begin());
    for_each_node(*this, [&](tNodeMutableIter node) {
        std::transform(node.begin(tag), node.end(tag),
                       node.begin(tag), tReorderFunc(ordering));
    });
    for_each_edge(*this, [&](tEdgeMutableIter edge) {
        std::transform(edge.begin(tag), edge.end(tag),
                       edge.begin(tag), tReorderFunc(ordering));
    });
    for_each_face(*this, [&](tFaceMutableIter face) {
        std::transform(face.begin(tag), face.end(tag),
                       face.begin(tag), tReorderFunc(ordering));
    });
    for_each_cell(*this, [&](tCellMutableIter cell) {
        std::transform(cell.begin(tag), cell.end(tag),
                       cell.begin(tag), tReorderFunc(ordering));
    });
}   // cMesh::fix_permutation_and_adjacency_

/**
 * Change order of all nodes.
 */
void cMesh::permute_nodes(std::vector<uint_t>&& node_permutation) {
    fix_permutation_and_adjacency_(eNodeTag, node_permutation);
    permute_rows(
        node_permutation.begin(), node_permutation.end(),
        m_node_nodes, m_node_edges, m_node_faces, m_node_cells);
    permute_inplace(
        node_permutation.begin(), node_permutation.end(),
        m_node_marks.begin(), m_node_coords.begin());
}   // сMesh::permute_nodes

/**
 * Change order of all edges.
 */
void cMesh::permute_edges(std::vector<uint_t>&& edge_permutation) {
    fix_permutation_and_adjacency_(eEdgeTag, edge_permutation);
    permute_rows(
        edge_permutation.begin(), edge_permutation.end(),
        m_edge_nodes, m_edge_edges, m_edge_faces, m_edge_cells);
    permute_inplace(
        edge_permutation.begin(), edge_permutation.end(),
        m_edge_marks.begin(), m_edge_shapes.begin(),
        m_edge_lengths.begin(), m_edge_directions.begin());
}   // сMesh::permute_edges

/**
 * Change order of all faces.
 */
void cMesh::permute_faces(std::vector<uint_t>&& face_permutation) {
    fix_permutation_and_adjacency_(eFaceTag, face_permutation);
    permute_rows(
        face_permutation.begin(), face_permutation.end(),
        m_face_nodes, m_face_edges, m_face_faces, m_face_cells);
    permute_inplace(
        face_permutation.begin(), face_permutation.end(),
        m_face_marks.begin(), m_face_shapes.begin(),
        m_face_areas.begin(), m_face_normals.begin(), m_face_center_coords.begin());
}   // сMesh::permute_faces

/**
 * Change order of all cells.
 */
void cMesh::permute_cells(std::vector<uint_t>&& cell_permutation) {
    fix_permutation_and_adjacency_(eCellTag, cell_permutation);
    permute_rows(
        cell_permutation.begin(), cell_permutation.end(),
        m_cell_nodes, m_cell_edges, m_cell_faces, m_cell_cells);
    permute_inplace(
        cell_permutation.begin(), cell_permutation.end(),
        m_cell_marks.begin(), m_cell_shapes.begin(),
        m_cell_volumes.begin(), m_cell_center_coords.begin());
}   // сMesh::permute_cells

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

void cMesh::reorder_faces() {
    /** @todo Refactor me! */
    {
        std::vector<uint_t> node_reordering(num_nodes());
        std::iota(node_reordering.begin(), node_reordering.end(), 0);
        permute_nodes(std::move(node_reordering));
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
        permute_edges(std::move(edge_reordering));
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
        permute_faces(std::move(face_reordering));
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
        permute_cells(std::move(cell_reordering));
        for (uint_t cell_index = 0; cell_index < num_cells(); ++cell_index) {
            m_marked_cell_ranges.resize(get_mark(eCellTag, cell_index) + 2);
            m_marked_cell_ranges[get_mark(eCellTag, cell_index) + 1] += 1;
        }
        std::partial_sum(
            m_marked_cell_ranges.begin(), m_marked_cell_ranges.end(), m_marked_cell_ranges.begin());
    }
}   // сMesh::permute_faces

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/* Face edge node index tables for various face types. */
static const std::map<eShape, std::vector<std::vector<uint_t>>> g_face_to_edge_index = {
    /* 1D faces. */
    { eShape::node, { {0} } },
    /* 2D faces. */
    { eShape::segment_2, { {0}, {1} } },
    /* 3D faces. */
    { eShape::triangle_3, { {1, 2}, {2, 0}, {0, 1} } },
    { eShape::quadrangle_4, { {0, 1}, {1, 2}, {2, 3}, {3, 0} } },
};

/**
 * Generate edges using the face to node connectivity.
 * @warning This function may be slow and memory-consuming.
 */
void cMesh::generate_edges() {
    /* Build a map of the existing edges.
     * ( An edge is uniquely defined by a set of it's nodes. ) */
    std::map<std::set<uint_t>, uint_t> edge_cache;
    std::for_each(begin_edge(*this), end_edge(*this), [&](tEdgeIter edge) {
        edge_cache.emplace(
            std::set<uint_t>(edge.begin(eNodeTag), edge.end(eNodeTag)), edge);
    });

    /* Generate missing edges. */
    std::for_each(begin_face(*this), end_face(*this), [&](tFaceMutableIter face) {
        for (uint_t edge_local = 0; edge_local < face.num_edges(); ++edge_local) {
            uint_t& edge_index = face.begin(eEdgeTag)[edge_local];
            if (edge_index != npos) {
                continue;
            }

            /* Collect edge nodes using the table. */
            std::set<uint_t> edge_nodes;
            for (uint_t node_local : g_face_to_edge_index.at(face.get_shape()).at(edge_local)) {
                FEATHERS_ASSERT(node_local < face.num_nodes());
                const uint_t node_index = face.begin(eNodeTag)[node_local];
                edge_nodes.insert(node_index);
            }

            /* Create the face or add current cell to the adjacency list. */
            const auto edge_cache_iter = edge_cache.find(edge_nodes);
            if (edge_cache_iter != edge_cache.cend()) {
                /* Locate the existing edge. */
                edge_index = edge_cache_iter->second;
            } else {
                /* Create a brand-new edge. */
                edge_index = insert_edge( // TODO:
                    { eShape::segment_2, std::vector<uint_t>(edge_nodes.begin(), edge_nodes.end()) });
                edge_cache.emplace(edge_nodes, edge_index);
                set_mark(eEdgeTag, edge_index, face.get_mark());
            }
        }
    });

    /* Check faces:
     * each face should be connected to all edges. */
    for_each_face(*this, [](tFaceIter face) {
        FEATHERS_ENSURE(
            std::all_of(face.begin(eEdgeTag), face.end(eEdgeTag), is_not_npos));
    });
}   // сMesh::generate_edges

/* Cell face node index tables for various cell types. */
static const std::map<eShape, std::vector<std::vector<uint_t>>> g_shape_to_face_nodes = {
    /* 1D cells. */
    { eShape::segment_2,
        { {0}, {1} } },
    /* 2D cells. */
    { eShape::triangle_3,
        { {1, 2}, {2, 0}, {0, 1} } },
    { eShape::quadrangle_4,
        { {0, 1}, {1, 2}, {2, 3}, {3, 0} } },
    /* 3D cells.
     * TODO: Add other cell types! */
    { eShape::tetrahedron_4,
        { {0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3} } },
    { eShape::hexahedron_8,
        { {0, 4, 7, 3}, {1, 2, 6, 5}, {3, 7, 6, 2}, {0, 1, 5, 4}, {0, 3, 2, 1}, {4, 5, 6, 7} } },
};

/**
 * Generate faces using the cell to node connectivity.
 * @warning This function may be slow and memory-consuming.
 */
void cMesh::generate_faces() {
    /* Build a map of the existing faces.
     * ( A face can be uniquely defined by the set of it's nodes. ) */
    std::map<std::set<uint_t>, uint_t> face_cache;
    std::for_each(begin_face(*this), end_face(*this), [&](tFaceIter face) {
        face_cache.emplace(
            std::set<uint_t>(face.begin(eNodeTag), face.end(eNodeTag)), face);
    });

    /* Generate missing faces. */
    std::for_each(begin_cell(*this), end_cell(*this), [&](tCellMutableIter cell) {
        for (uint_t face_local = 0; face_local < cell.num_faces(); ++face_local) {
            uint_t& face_index = cell.begin(eFaceTag)[face_local];
            if (face_index != npos) {
                continue;
            }

            /* Collect face nodes using the table. */
            std::vector<uint_t> face_nodes;
            for (uint_t node_local : g_shape_to_face_nodes.at(cell.get_shape()).at(face_local)) {
                FEATHERS_ASSERT(node_local < cell.num_nodes());
                const uint_t node_index = cell.begin(eNodeTag)[node_local];
                face_nodes.push_back(node_index);
            }

            /* Create the face or add current cell to the adjacency list. */
            std::set<uint_t> face_nodes_set(face_nodes.begin(), face_nodes.end());
            const auto face_cache_iter = face_cache.find(face_nodes_set);
            if (face_cache_iter != face_cache.cend()) {
                /* Locate the existing face.
                 * ( We should be careful with face orientation here. ) */
                face_index = face_cache_iter->second;
                tFaceMutableIter face(*this, face_index);
                if (std::equal(face.begin(eNodeTag), face.end(eNodeTag), face_nodes.cbegin())) {
                    FEATHERS_ASSERT(face.begin(eCellTag)[eFaceInnerCell] == npos);
                    face.begin(eCellTag)[eFaceInnerCell] = cell;
                } else {
                    FEATHERS_ASSERT(face.begin(eCellTag)[eFaceOuterCell] == npos);
                    face.begin(eCellTag)[eFaceOuterCell] = cell;
                }
            } else {
                /* Create a brand-new face.
                 * Assign the current cell as the inner one. */
                // TODO: correct face shape!
                face_index = insert_face(eShape::segment_2, face_nodes);
                face_cache.emplace(face_nodes_set, face_index);
                tFaceMutableIter face(*this, face_index);
                face.begin(eCellTag)[eFaceInnerCell] = cell;
            }
        }
    });

    /* Check cells:
     * each cell should be connected to all faces. */
    for_each_cell(*this, [](tCellIter cell) {
        FEATHERS_ENSURE(std::all_of(cell.begin(eFaceTag), cell.end(eFaceTag), is_not_npos));
    });
    /* Check faces:
     * each internal face should be connected to two cells;
     * each boundary face should be connected to at least one cell. */
    for_each_face(*this, [](tFaceIter face) {
        if (face.get_mark() == 0) {
            FEATHERS_ENSURE(
                std::all_of(face.begin(eCellTag), face.end(eCellTag), is_not_npos));
        } else {
            FEATHERS_ENSURE(
                std::any_of(face.begin(eCellTag), face.end(eCellTag), is_not_npos));
        }
    });
}   // сMesh::generate_faces

/* A node and edge flip table for various face types. */
static const std::map<eShape, std::pair<std::vector<uint_t>, std::vector<uint_t>>>
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
void cMesh::generate_boundary_cells() {
    std::for_each(begin_face(*this), end_face(*this), [&](tFaceMutableIter face) {
        if (face.get_mark() == 0) {
            return;
        }

        /* Boundary faces should be oriented outwards from the mesh. */
        if (face.get_inner_cell() == npos) {
            /* Flip normal and cell connectivity. */
            face.set_normal(-face.get_normal());
            std::swap(face.begin(eCellTag)[eFaceInnerCell],
                      face.begin(eCellTag)[eFaceOuterCell]);
            /* Flip node and edge connectivity. */
            std::vector<uint_t> node_permutation, edge_permutation;
            std::tie(node_permutation, edge_permutation) =
                g_face_shape_to_nodes_and_edges_flip.at(face.get_shape());
            permute_inplace(
                node_permutation.begin(), node_permutation.end(), face.begin(eNodeTag));
            permute_inplace(
                edge_permutation.begin(), edge_permutation.end(), face.begin(eNodeTag));
        }
        tCellIter cell = face.get_inner_cell();

        /* Generate the boundary cell: reflect a connected interior cell. */
        uint_t ghost_cell_dim;
        std::vector<uint_t> ghost_cell_nodes;
#if 0
        ghost_cell_dim = m_dim - 1;
        ghost_cell_nodes.assign(face.begin_node(), face.end_node());
#endif
#if 1
        ghost_cell_dim = m_dim;
        cell.for_each_node([&](tNodeIter node) {
            if (std::find(face.begin(eNodeTag), face.end(eNodeTag), node) == face.end(eNodeTag)) {
                /* Reflect an interior cell node. */
                // TODO: face normals are not computed here!
                // TODO: https://glm.g-truc.net/0.9.5/api/a00157.html#gab63646fc36b81cf69d3ce123a72f76f2
                vec3_t node_coords = node.get_coords();
                const vec3_t delta = node_coords - face.get_center_coords();
                node_coords -= 2.0*glm::dot(delta, face.get_normal())*face.get_normal();
                ghost_cell_nodes.push_back(
                    insert_node(node_coords, face.get_mark()));
            } else {
                /* Insert a boundary face node. */
                ghost_cell_nodes.push_back(node);
            }
        });
#endif
        /* Insert the boundary cell. */
        // TODO:
        const uint_t boundary_cell_index =
            insert_cell(eShape::triangle_3, ghost_cell_nodes, face.get_mark());
#if 0
        Cell& boundary_cell = get_cell(boundary_cell_index);
        while (boundary_cell._num_faces() != 1) {
            boundary_cell.erase_face(0);
        }
        boundary_cell.begin_face()[0] = face;
#endif
        face.begin(eCellTag)[eFaceOuterCell] = boundary_cell_index;
    });
}   // сMesh::generate_boundary_cells

}   // namespace feathers

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //
