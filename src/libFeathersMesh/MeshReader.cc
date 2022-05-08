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

#include <libFeathersMesh/Mesh.hh>

#include <fstream>
#include <iomanip>

namespace feathers {

bool Mesh::read_triangle(const char *path) {
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
    FEATHERS_ENSURE(node_index == EmplaceNode(node_coords));
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
      face_index == EmplaceFace({eShape::segment_2, face_nodes}, FaceMark(mark)));
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
      cell_index == EmplaceCell({eShape::triangle_3, cell_nodes}));
    std::getline(cell_file, line);
  }

  finalize();
  return true;
} // сMesh::read_triangle

#if 0
bool Mesh::read_tetgen(const char *path) {
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
      node_index == EmplaceNode(node_coords));
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
      face_index == EmplaceFace({eShape::triangle_3, face_nodes}, mark));
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
      cell_index == EmplaceCell({eShape::tetrahedron_4, cell_nodes}));
    std::getline(cell_file, line);
  }

  finalize();
  return true;
} // сMesh::read_tetgen

bool Mesh::read_image2D(const char *path,
                        const std::map<sPixel, uint_t>& mark_colors,
                        sPixel fluid_color,
                        vec2_t pixel_size) {
  cImage2D image;
  image.load(path);

  cImage2D nodes_image;
  nodes_image.init(image.width() + 1, image.height() + 1, sPixel(0, 0, 0, 0));

  uint_t node_index = 0;
  for (uint_t y = 1; y < image.height() - 1; ++y) {
    for (uint_t x = 1; x < image.width() - 1; ++x) {
      if (image(x, y).rgba != fluid_color.rgba) {
        continue;
      }

      const vec2_t cell_center_coords = pixel_size*vec2_t(x - 0.5, y - 0.5);

      /* Insert or query the cell nodes. */
      uint_t& lw_lt_node_index = nodes_image(x + 0, y + 0).rgba;
      if (lw_lt_node_index == 0) {
        lw_lt_node_index = node_index++;
        const vec3_t node_coords(cell_center_coords + pixel_size*vec2_t(-0.5, -0.5), 0.0);
        FEATHERS_ENSURE(lw_lt_node_index == EmplaceNode(node_coords));
      }
      uint_t& lw_rt_node_index = nodes_image(x + 1, y + 0).rgba;
      if (lw_rt_node_index == 0) {
        lw_rt_node_index = node_index++;
        const vec3_t node_coords(cell_center_coords + pixel_size*vec2_t(+0.5, -0.5), 0.0);
        FEATHERS_ENSURE(lw_rt_node_index == EmplaceNode(node_coords));
      }
      uint_t& up_rt_node_index = nodes_image(x + 1, y + 1).rgba;
      if (up_rt_node_index == 0) {
        up_rt_node_index = node_index++;
        const vec3_t node_coords(cell_center_coords + pixel_size*vec2_t(+0.5, +0.5), 0.0);
        FEATHERS_ENSURE(up_rt_node_index == EmplaceNode(node_coords));
      }
      uint_t& up_lt_node_index = nodes_image(x + 0, y + 1).rgba;
      if (up_lt_node_index == 0) {
        up_lt_node_index = node_index++;
        const vec3_t node_coords(cell_center_coords + pixel_size*vec2_t(-0.5, +0.5), 0.0);
        FEATHERS_ENSURE(up_lt_node_index == EmplaceNode(node_coords));
      }

      /* Insert the cell. */
      EmplaceCell({eShape::quadrangle_4,
                    {lw_lt_node_index, lw_rt_node_index, up_rt_node_index, up_lt_node_index}});

      /* Insert the boundary faces. */
      if (const sPixel loPixel = image(x, y - 1); loPixel.rgba != fluid_color.rgba) {
        const uint_t mark = mark_colors.at(loPixel);
        EmplaceFace({eShape::segment_2, {lw_lt_node_index, lw_rt_node_index}}, mark);
      }
      if (const sPixel rtPixel = image(x + 1, y); rtPixel.rgba != fluid_color.rgba) {
        const uint_t mark = mark_colors.at(rtPixel);
        EmplaceFace({eShape::segment_2, {lw_rt_node_index, up_rt_node_index}}, mark);
      }
      if (const sPixel upPixel = image(x, y + 1); upPixel.rgba != fluid_color.rgba) {
        const uint_t mark = mark_colors.at(upPixel);
        EmplaceFace({eShape::segment_2, {up_rt_node_index, up_lt_node_index}}, mark);
      }
      if (const sPixel ltPixel = image(x - 1, y); ltPixel.rgba != fluid_color.rgba) {
        const uint_t mark = mark_colors.at(ltPixel);
        EmplaceFace({eShape::segment_2, {up_lt_node_index, lw_lt_node_index}}, mark);
      }
    }
  }

  finalize();
  return true;
} // Mesh::read_image2D
#endif

void Mesh::save_vtk(const char* path,
                    const std::vector<sFieldDesc>& fields) const {
  std::ofstream file(path);
  file << std::setprecision(std::numeric_limits<real_t>::digits10 + 1);
  file << "# vtk DataFile Version 2.0" << std::endl;
  file << "# Generated by Feathers/StormRuler/Mesh2VTK" << std::endl;
  file << "ASCII" << std::endl;
  file << "DATASET UNSTRUCTURED_GRID" << std::endl;

  file << "POINTS " << NumNodes() << " double" << std::endl;
  std::for_each(feathers::BeginNode(*this), feathers::EndNode(*this), [&](NodeIter node) {
    const vec3_t& node_coords = node.Pos();
    file << node_coords.x << " " << node_coords.y << " " << node_coords.z << std::endl;
  });
  file << std::endl;

  size_t const total_num_cell_nodes = for_range_sum(
    BeginInteriorCell(*this), EndInteriorCell(*this), size_t(0), [](CellIter cell) {
      return cell.NumNodes() + 1;
    });
  file << "CELLS " << NumCells({}) << " " << total_num_cell_nodes << std::endl;
  std::for_each(BeginInteriorCell(*this), EndInteriorCell(*this), [&](CellIter cell) {
    file << cell.NumNodes() << " ";
    cell.ForEachNode([&](uint_t node_index) {
      file << node_index << " ";
    });
    file << std::endl;
  });
  file << std::endl;

  file << "CELL_TYPES " << NumCells({}) << std::endl;
  std::for_each(BeginInteriorCell(*this), EndInteriorCell(*this), [&](CellIter cell) {
    static const std::map<eShape, const char*> shapes = {
      { eShape::node, "1" }, { eShape::segment_2, "2" },
      { eShape::triangle_3, "5" }, { eShape::quadrangle_4, "9" },
      { eShape::tetrahedron_4, "10" }, { eShape::pyramid_5, "14" },
      { eShape::pentahedron_6, "13" }, { eShape::hexahedron_8, "12" }
    };
    file << shapes.at(cell.Shape()) << std::endl;
  });
  file << std::endl;

  file << "CELL_DATA " << NumCells({}) << std::endl;
  for (const sFieldDesc& field : fields) {
    file << "SCALARS " << field.name << " double 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    std::for_each(BeginInteriorCell(*this), EndInteriorCell(*this), [&](CellIter cell) {
      file << (*field.scalar)[cell][field.var_index] << std::endl;
    });
  }
  file << std::endl;
} // Mesh::save_vtk

} // namespace feathers
