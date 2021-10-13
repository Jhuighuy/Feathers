/*
 *  ______  ______   ______   ______  __  __   ______   ______   ______
 * /\  ___\/\  ___\ /\  __ \ /\__  _\/\ \_\ \ /\  ___\ /\  __ \ /\  ___\
 * \ \  __\\ \  _\  \ \  __ \\/_/\ \/\ \  __ \\ \  __\ \ \  __/ \ \___  \
 *  \ \_\   \ \_____\\ \_\ \_\  \ \_\ \ \_\ \_\\ \_____\\ \_\ \_\\/\_____\
 *   \/_/    \/_____/ \/_/\/_/   \/_/  \/_/\/_/ \/_____/ \/_/ /_/ \/_____/
 *
 * Copyright (c) 2021 Oleg Butakov
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
#include <libFeathersUtils/Permute.hh>

#include <set>
#include <map>
#include <fstream>
#include <iomanip>

namespace feathers {

bool cMesh::read_triangle(const char *path) {
    std::string line;

    std::ifstream node_file(path + std::string("node"));
    FEATHERS_ENSURE(node_file.is_open());
    uint_t num_nodes = 0;
    node_file >> num_nodes >> m_dim;
    std::getline(node_file, line);
    reserve_nodes(num_nodes);
    for (uint_t i = 0; i < num_nodes; ++i) {
        uint_t node_index = 0;
        vec3_t node_coords(0.0);
        node_file >> node_index >> node_coords.x >> node_coords.y;
        std::getline(node_file, line);
        FEATHERS_ENSURE(node_index == emplace_back_node(node_coords));
    }

    std::ifstream face_file(path + std::string("edge"));
    FEATHERS_ENSURE(face_file.is_open());
    uint_t num_faces = 0;
    face_file >> num_faces;
    std::getline(face_file, line);
    reserve_faces(num_faces);
    for (uint_t i = 0; i < num_faces; ++i) {
        uint_t face_index = 0;
        std::vector<uint_t> face_nodes(2, npos);
        uint_t mark = 0;
        face_file >> face_index >> face_nodes[0] >> face_nodes[1] >> mark;
        FEATHERS_ENSURE(
            face_index == emplace_back_face({eShape::segment_2, face_nodes}, mark));
        std::getline(face_file, line);
    }

    std::ifstream cell_file(path + std::string("ele"));
    FEATHERS_ENSURE(cell_file.is_open());
    uint_t num_cells = 0;
    cell_file >> num_cells;
    std::getline(cell_file, line);
    reserve_cells(num_cells);
    for (uint_t i = 0; i < num_cells; ++i) {
        uint_t cell_index = 0;
        std::vector<uint_t> cell_nodes(3, npos);
        cell_file >> cell_index >> cell_nodes[0] >> cell_nodes[1] >> cell_nodes[2];
        FEATHERS_ENSURE(
            cell_index == emplace_back_cell({eShape::triangle_3, cell_nodes}));
        std::getline(cell_file, line);
    }

    generate_faces();
    generate_boundary_cells();
    reorder_faces();
    compute_all_shape_properties();
    return true;
} // сMesh::read_triangle

bool cMesh::read_tetgen(const char *path) {
    std::string line;

    std::ifstream node_file(path + std::string("node"));
    FEATHERS_ENSURE(node_file.is_open());
    uint_t num_nodes = 0;
    node_file >> num_nodes >> m_dim;
    std::getline(node_file, line);
    reserve_nodes(num_nodes);
    for (uint_t i = 0; i < num_nodes; ++i) {
        uint_t node_index = 0;
        vec3_t node_coords(0.0);
        node_file >> node_index >> node_coords.x >> node_coords.y >> node_coords.z;
        FEATHERS_ENSURE(
            node_index == emplace_back_node(node_coords));
        std::getline(node_file, line);
    }

    std::ifstream face_file(path + std::string("face"));
    FEATHERS_ENSURE(face_file.is_open());
    uint_t num_faces = 0;
    face_file >> num_faces;
    std::getline(face_file, line);
    reserve_faces(num_faces);
    for (uint_t i = 0; i < num_faces; ++i) {
        uint_t face_index = 0;
        std::vector<uint_t> face_nodes(3, npos);
        uint_t mark = 0;
        face_file >> face_index >> face_nodes[0] >> face_nodes[1] >> face_nodes[2] >> mark;
        FEATHERS_ENSURE(
            face_index == emplace_back_face({eShape::triangle_3, face_nodes}, mark));
        std::getline(face_file, line);
    }

    std::ifstream cell_file(path + std::string("ele"));
    FEATHERS_ENSURE(cell_file.is_open());
    uint_t num_cells = 0;
    cell_file >> num_cells;
    std::getline(cell_file, line);
    reserve_cells(num_cells);
    for (uint_t i = 0; i < num_cells; ++i) {
        uint_t cell_index = 0;
        std::vector<uint_t> cell_nodes(4, npos);
        cell_file >> cell_index >> cell_nodes[0] >> cell_nodes[1] >> cell_nodes[2] >> cell_nodes[3];
        FEATHERS_ENSURE(
            cell_index == emplace_back_cell({eShape::tetrahedron_4, cell_nodes}));
        std::getline(cell_file, line);
    }

    generate_faces();
    generate_boundary_cells();
    reorder_faces();
    compute_all_shape_properties();
    return true;
} // сMesh::read_tetgen

bool cMesh::read_image(const char *path,
                       const std::map<sPixel, uint_t>& mark_colors,
                       sPixel fluid_color,
                       vec2_t pixel_size) {
    cImage image;
    image.load(path);

    cImage nodes_image;
    nodes_image.init(image.width() + 1, image.height() + 1, sPixel(0, 0, 0, 0));

    uint_t node_index = 0;
    for (uint_t y = 1; y < image.height() - 1; ++y) {
        for (uint_t x = 1; x < image.width() - 1; ++x) {
            if (image(x, y).rgba != fluid_color.rgba) {
                continue;
            }

            const vec2_t cell_center_coords = pixel_size * vec2_t(x - 0.5, y - 0.5);

            /* Insert or query the cell nodes. */
            uint_t& node_index_ll = nodes_image(x + 0, y + 0).rgba;
            if (node_index_ll == 0) {
                node_index_ll = node_index++;
                const vec3_t node_coords(
                    cell_center_coords + pixel_size*vec2_t(-0.5, -0.5), 0.0);
                FEATHERS_ENSURE(node_index_ll == emplace_back_node(node_coords));
            }
            uint_t& node_index_lr = nodes_image(x + 1, y + 0).rgba;
            if (node_index_lr == 0) {
                node_index_lr = node_index++;
                const vec3_t node_coords(
                    cell_center_coords + pixel_size*vec2_t(+0.5, -0.5), 0.0);
                FEATHERS_ENSURE(node_index_lr == emplace_back_node(node_coords));
            }
            uint_t& node_index_ur = nodes_image(x + 1, y + 1).rgba;
            if (node_index_ur == 0) {
                node_index_ur = node_index++;
                const vec3_t node_coords(
                    cell_center_coords + pixel_size*vec2_t(+0.5, +0.5), 0.0);
                FEATHERS_ENSURE(node_index_ur == emplace_back_node(node_coords));
            }
            uint_t& node_index_ul = nodes_image(x + 0, y + 1).rgba;
            if (node_index_ul == 0) {
                node_index_ul = node_index++;
                const vec3_t node_coords(
                    cell_center_coords + pixel_size*vec2_t(-0.5, +0.5), 0.0);
                FEATHERS_ENSURE(node_index_ul == emplace_back_node(node_coords));
            }

            /* Insert the cell. */
            emplace_back_cell(
                {eShape::quadrangle_4, {node_index_ll, node_index_lr, node_index_ur, node_index_ul}});

            /* Insert the boundary faces. */
            if (const sPixel lower_pixel = image(x, y - 1); lower_pixel.rgba != fluid_color.rgba) {
                const uint_t mark = mark_colors.at(lower_pixel);
                emplace_back_face({eShape::segment_2, {node_index_ll, node_index_lr}}, mark);
            }
            if (const sPixel right_pixel = image(x + 1, y); right_pixel.rgba != fluid_color.rgba) {
                const uint_t mark = mark_colors.at(right_pixel);
                emplace_back_face({eShape::segment_2, {node_index_lr, node_index_ur}}, mark);
            }
            if (const sPixel upper_pixel = image(x, y + 1); upper_pixel.rgba != fluid_color.rgba) {
                const uint_t mark = mark_colors.at(upper_pixel);
                emplace_back_face({eShape::segment_2, {node_index_ur, node_index_ul}}, mark);
            }
            if (const sPixel left_pixel = image(x - 1, y); left_pixel.rgba != fluid_color.rgba) {
                const uint_t mark = mark_colors.at(left_pixel);
                emplace_back_face({eShape::segment_2, {node_index_ul, node_index_ll}}, mark);
            }
        }
    }
    
    generate_faces();
    generate_boundary_cells();
    reorder_faces();
    compute_all_shape_properties();
    return true;
} // cMesh::read_image

void cMesh::save_vtk(const char* path,
                     const std::vector<sFieldDesc>& fields) const {
    std::ofstream file(path);
    file << std::setprecision(std::numeric_limits<real_t>::digits10 + 1);
    file << "# vtk DataFile Version 2.0" << std::endl;
    file << "kek" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;

    file << "POINTS " << num_nodes() << " double" << std::endl;
    std::for_each(begin_node(*this), end_node(*this), [&](tNodeIter node) {
        const vec3_t& node_coords = node.get_coords();
        file << node_coords.x << " " << node_coords.y << " " << node_coords.z << std::endl;
    });
    file << std::endl;

    const size_t total_num_cell_nodes = for_range_sum(
        begin_interior_cell(*this), end_interior_cell(*this), size_t(0), [](tCellIter cell) {
            return cell.num_nodes() + 1;
        });
    file << "CELLS " << num_marked_cells(0) << " " << total_num_cell_nodes << std::endl;
    std::for_each(begin_interior_cell(*this), end_interior_cell(*this), [&](tCellIter cell) {
        file << cell.num_nodes() << " ";
        cell.for_each_node([&](uint_t node_index) {
            file << node_index << " ";
        });
        file << std::endl;
    });
    file << std::endl;

    file << "CELL_TYPES " << num_marked_cells(0) << std::endl;
    std::for_each(begin_interior_cell(*this), end_interior_cell(*this), [&](tCellIter cell) {
        static const std::map<eShape, const char*> shapes = {
            { eShape::node, "1" }, { eShape::segment_2, "2" },
            { eShape::triangle_3, "5" }, { eShape::quadrangle_4, "9" },
            { eShape::tetrahedron_4, "10" }, { eShape::pyramid_5, "14" },
            { eShape::pentahedron_6, "13" }, { eShape::hexahedron_8, "12" }
        };
        file << shapes.at(cell.get_shape()) << std::endl;
    });
    file << std::endl;

    file << "CELL_DATA " << num_marked_cells(0) << std::endl;
    for (const sFieldDesc& field : fields) {
        file << "SCALARS " << field.name << " double 1" << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;
        std::for_each(begin_interior_cell(*this), end_interior_cell(*this), [&](tCellIter cell) {
            file << (*field.scalar)[cell][field.var_index] << std::endl;
        });
    }
} // cMesh::save_vtk

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
} // сMesh::compute_edge_shape_properties

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
} // сMesh::compute_face_shape_properties

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
} // сMesh::compute_cell_shape_properties

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Preallocate the node storage.
 */
void cMesh::reserve_nodes(uint_t node_capacity) {
    m_node_marks.reserve(node_capacity);

    m_node_coords.reserve(node_capacity);

    m_node_nodes.reserve_rows(node_capacity);
    m_node_edges.reserve_rows(node_capacity);
    m_node_faces.reserve_rows(node_capacity);
    m_node_cells.reserve_rows(node_capacity);
} // cMesh::reserve_nodes

/**
 * Preallocate the edge storage.
 */
void cMesh::reserve_edges(uint_t edge_capacity) {
    m_edge_marks.reserve(edge_capacity);

    m_edge_shapes.reserve(edge_capacity);
    m_edge_lengths.reserve(edge_capacity);
    m_edge_directions.reserve(edge_capacity);

    m_edge_nodes.reserve_rows(edge_capacity);
    m_edge_edges.reserve_rows(edge_capacity);
    m_edge_faces.reserve_rows(edge_capacity);
    m_edge_cells.reserve_rows(edge_capacity);
} // cMesh::reserve_edges

/**
 * Preallocate the face storage.
 */
void cMesh::reserve_faces(uint_t face_capacity) {
    m_face_marks.reserve(face_capacity);

    m_face_shapes.reserve(face_capacity);
    m_face_areas.reserve(face_capacity);
    m_face_normals.reserve(face_capacity);
    m_face_center_coords.reserve(face_capacity);

    m_face_nodes.reserve_rows(face_capacity);
    m_face_edges.reserve_rows(face_capacity);
    m_face_faces.reserve_rows(face_capacity);
    m_face_cells.reserve_rows(face_capacity);
} // cMesh::reserve_faces

/**
 * Preallocate the cell storage.
 */
void cMesh::reserve_cells(uint_t cell_capacity) {
    m_cell_marks.reserve(cell_capacity);

    m_cell_shapes.reserve(cell_capacity);
    m_cell_volumes.reserve(cell_capacity);
    m_cell_center_coords.reserve(cell_capacity);

    m_cell_nodes.reserve_rows(cell_capacity);
    m_cell_edges.reserve_rows(cell_capacity);
    m_cell_faces.reserve_rows(cell_capacity);
    m_cell_cells.reserve_rows(cell_capacity);
} // cMesh::reserve_cells

/**
 * Emplace a new node into the mesh.
 * @returns Index of the inserted node.
 */
uint_t cMesh::emplace_back_node(const vec3_t& node_coords, uint_t mark) {
    const uint_t node_index = m_num_nodes++;

    m_node_marks.emplace_back(mark);

    m_node_coords.emplace_back(node_coords);

    m_node_nodes.emplace_back_row();
    m_node_edges.emplace_back_row();
    m_node_faces.emplace_back_row();
    m_node_cells.emplace_back_row();

    return node_index;
} // сMesh::emplace_back_node

/**
 * Emplace a new edge into the mesh.
 * @returns Index of the inserted edge.
 */
uint_t cMesh::emplace_back_edge(std::unique_ptr<iElement>&& edge, uint_t mark) {
    const uint_t edge_index = m_num_edges++;

    m_edge_marks.emplace_back(mark);

    m_edge_shapes.emplace_back(edge->get_shape());
    m_edge_lengths.emplace_back(edge->get_length_or_area_or_volume());
    m_edge_directions.emplace_back(edge->get_direction());

    m_edge_nodes.emplace_back_row(edge->get_nodes().begin(), edge->get_nodes().end());
    m_edge_edges.emplace_back_row();
    m_edge_faces.emplace_back_row();
    m_edge_cells.emplace_back_row();

    return edge_index;
} // сMesh::emplace_back_edge

/**
 * Emplace a new face into the mesh.
 * @returns Index of the inserted face.
 */
uint_t cMesh::emplace_back_face(std::unique_ptr<iElement>&& face, uint_t mark) {
    const uint_t face_index = m_num_faces++;

    m_face_marks.emplace_back(mark);

    m_face_shapes.emplace_back(face->get_shape());
    m_face_areas.emplace_back(face->get_length_or_area_or_volume());
    m_face_normals.emplace_back(face->get_normal());
    m_face_center_coords.emplace_back(face->get_center_coords());

    m_face_nodes.emplace_back_row(face->get_nodes().begin(), face->get_nodes().end());
    m_face_edges.emplace_back_row(face->num_edges());
    m_face_faces.emplace_back_row();
    m_face_cells.emplace_back_row(2);

    return face_index;
} // сMesh::emplace_back_face

/**
 * Insert a new cell into the mesh.
 * @returns Index of the inserted cell.
 */
uint_t cMesh::emplace_back_cell(std::unique_ptr<iElement>&& cell, uint_t mark) {
    const uint_t cell_index = m_num_cells++;

    m_cell_marks.emplace_back(mark);

    m_cell_shapes.emplace_back(cell->get_shape());
    m_cell_volumes.emplace_back(cell->get_length_or_area_or_volume());
    m_cell_center_coords.emplace_back(cell->get_center_coords());

    m_cell_nodes.emplace_back_row(cell->get_nodes().begin(), cell->get_nodes().end());
    m_cell_edges.emplace_back_row(cell->num_edges());
    m_cell_faces.emplace_back_row(cell->num_faces());
    m_cell_cells.emplace_back_row();

    return cell_index;
} // сMesh::emplace_back_cell

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

template<typename tTag>
void cMesh::fix_permutation_and_adjacency_(tTag tag, std::vector<uint_t>& permutation) {
    //if (permutation.empty()) {
    //    /* Generate identity permutation. */
    //    permutation.reserve();
    //}
    /* Fix permutations by resorting it by marks.
     * Stable sort is used here to preserve the permutation in the best possible way. */
    std::stable_sort(permutation.begin(), permutation.end(),
                     [&](uint_t index_1, uint_t index_2) {
        return get_mark(tag, index_1) < get_mark(tag, index_2);
    });
    /* Fix adjacency lists. */
    std::vector<uint_t> ordering(permutation.size(), npos);
    convert_permutation_to_ordering(
        permutation.begin(), permutation.end(), ordering.begin());
    auto reorder_func = [&ordering](uint_t index) {
        return index != npos ? ordering[index] : npos;
    };
    for_each_node(*this, [&](tNodeMutableIter node) {
        std::transform(node.begin(tag), node.end(tag), node.begin(tag), reorder_func);
    });
    for_each_edge(*this, [&](tEdgeMutableIter edge) {
        std::transform(edge.begin(tag), edge.end(tag), edge.begin(tag), reorder_func);
    });
    for_each_face(*this, [&](tFaceMutableIter face) {
        std::transform(face.begin(tag), face.end(tag), face.begin(tag), reorder_func);
    });
    for_each_cell(*this, [&](tCellMutableIter cell) {
        std::transform(cell.begin(tag), cell.end(tag), cell.begin(tag), reorder_func);
    });
} // cMesh::fix_permutation_and_adjacency_

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
} // сMesh::permute_nodes

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
} // сMesh::permute_edges

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
} // сMesh::permute_faces

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
} // сMesh::permute_cells

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
} // сMesh::permute_faces

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Generate edges using the face to node connectivity.
 * @warning This function may be slow and memory-consuming.
 */
void cMesh::generate_edges() {
    /* An edge lookup table.
     * We assume that the edge can be uniquely identified by the set of its nodes. */
    std::map<std::set<uint_t>, uint_t> edge_lookup_table;

    std::for_each(begin_edge(*this), end_edge(*this), [&](tEdgeIter edge) {
        std::set<uint_t> edge_lookup_key(edge.begin(eNodeTag), edge.end(eNodeTag));
        edge_lookup_table.emplace(std::move(edge_lookup_key), edge);
    });

    /* For each cell-to-edge adjacency table entry:
     * find the face in the lookup table or emplace the new face. */
    std::for_each(begin_face(*this), end_face(*this), [&](tFaceMutableIter face) {
        tElementDescList edges_desc =
            face.get_element_object()->get_edges_desc();
        for (uint_t edge_local = 0; edge_local < face.num_edges(); ++edge_local) {
            uint_t& edge_index = face.begin(eEdgeTag)[edge_local];
            if (edge_index != npos) {
                continue;
            }

            /* Create the face or add current cell to the adjacency list. */
            sElementDesc& edge_desc = edges_desc[edge_local];
            std::set<uint_t> edge_lookup_key(
                edge_desc.node_indices.begin(), edge_desc.node_indices.end());
            if (edge_lookup_table.count(edge_lookup_key) == 0) {
                /* Create a brand-new edge. */
                edge_index = emplace_back_edge(std::move(edge_desc));
                edge_lookup_table.emplace(edge_lookup_key, edge_index);
                set_mark(eEdgeTag, edge_index, face.get_mark());
            } else {
                /* Edge exists. */
                edge_index = edge_lookup_table[edge_lookup_key];
            }
        }
    });

    /* Check faces:
     * each face should be connected to all edges. */
    for_each_face(*this, [](tFaceIter face) {
        FEATHERS_ENSURE(std::all_of(
            face.begin(eEdgeTag), face.end(eEdgeTag), is_not_npos));
    });
} // сMesh::generate_edges

/**
 * Generate faces using the cell to node connectivity.
 * @warning This function may be slow and memory-consuming.
 */
void cMesh::generate_faces() {
    /* A face lookup table.
     * We assume that the face can be uniquely identified by the set of its nodes. */
    std::map<std::set<uint_t>, uint_t> face_lookup_table;

    /* Add the existing faces to the lookup table. */
    std::for_each(begin_face(*this), end_face(*this), [&](tFaceIter face) {
        std::set<uint_t> face_lookup_key(face.begin(eNodeTag), face.end(eNodeTag));
        face_lookup_table.emplace(std::move(face_lookup_key), face);
    });

    /* For each cell-to-face adjacency table entry:
     * find the face in the lookup table or emplace the new face.
     * Also fill the face-to-cell adjacency table. */
    std::for_each(begin_cell(*this), end_cell(*this), [&](tCellMutableIter cell) {
        tElementDescList faces_desc =
            cell.get_element_object()->get_faces_desc();
        for (uint_t face_local = 0; face_local < cell.num_faces(); ++face_local) {
            uint_t& face_index = cell.begin(eFaceTag)[face_local];
            if (face_index != npos) {
                continue;
            }

            /* Create the face or add current cell to the adjacency list. */
            sElementDesc& face_desc = faces_desc[face_local];
            std::set<uint_t> face_lookup_key(face_desc.node_indices.begin(), face_desc.node_indices.end());
            if (face_lookup_table.count(face_lookup_key) == 0) {
                /* Create a brand-new face.
                 * Assign the current cell as the inner one. */
                face_index = emplace_back_face(std::move(face_desc));
                face_lookup_table.emplace(std::move(face_lookup_key), face_index);
                begin_adjacent_cell(eFaceTag, face_index)[eFaceInnerCell] = cell;
            } else {
                /* Face exists.
                 * Determine orientation and link it. */
                face_index = face_lookup_table[face_lookup_key];
                tFaceMutableIter face(*this, face_index);
                if (std::equal(face.begin(eNodeTag), face.end(eNodeTag), face_desc.node_indices.begin())) {
                    /* Face node order matches the order
                     * in the face description: face is inner. */
                    FEATHERS_ASSERT(face.begin(eCellTag)[eFaceInnerCell] == npos);
                    face.begin(eCellTag)[eFaceInnerCell] = cell;
                } else {
                    /* Otherwise, the face is outer. */
                    FEATHERS_ASSERT(face.begin(eCellTag)[eFaceOuterCell] == npos);
                    face.begin(eCellTag)[eFaceOuterCell] = cell;
                }
            }
        }
    });

    /* Check cells:
     * each cell should be connected to all faces. */
    for_each_cell(*this, [](tCellIter cell) {
        FEATHERS_ENSURE(std::all_of(
            cell.begin(eFaceTag), cell.end(eFaceTag), is_not_npos));
    });
    /* Check faces:
     * each internal face should be connected to two cells;
     * each boundary face should be connected to at least one cell. */
    for_each_face(*this, [](tFaceIter face) {
        if (face.get_mark() == 0) {
            std::vector c(face.begin(eCellTag), face.end(eCellTag));
            FEATHERS_ENSURE(std::all_of(
                face.begin(eCellTag), face.end(eCellTag), is_not_npos));
        } else {
            FEATHERS_ENSURE(std::any_of(
                face.begin(eCellTag), face.end(eCellTag), is_not_npos));
        }
    });
} // сMesh::generate_faces

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
                    emplace_back_node(node_coords, face.get_mark()));
            } else {
                /* Insert a boundary face node. */
                ghost_cell_nodes.push_back(node);
            }
        });
#endif
        /* Insert the boundary cell. */
        // TODO:
        const uint_t boundary_cell_index =
            emplace_back_cell({cell.get_shape(), ghost_cell_nodes}, face.get_mark());
#if 0
        Cell& boundary_cell = get_cell(boundary_cell_index);
        while (boundary_cell._num_faces() != 1) {
            boundary_cell.erase_face(0);
        }
        boundary_cell.begin_face()[0] = face;
#endif
        face.begin(eCellTag)[eFaceOuterCell] = boundary_cell_index;
    });
} // сMesh::generate_boundary_cells

} // namespace feathers
