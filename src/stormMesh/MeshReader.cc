/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person
/// obtaining a copy of this software and associated documentation
/// files (the "Software"), to deal in the Software without
/// restriction, including without limitation the rights  to use,
/// copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the
/// Software is furnished to do so, subject to the following
/// conditions:
///
/// The above copyright notice and this permission notice shall be
/// included in all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
/// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
/// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
/// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
/// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
/// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
/// OTHER DEALINGS IN THE SOFTWARE.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///

#include <fstream>
#include <iomanip>

#include <stormMesh/Mesh.hxx>

namespace Storm {

bool Mesh::read_from_triangle(std::string const& path) {
  std::string line;

  std::ifstream nodeStream(path + std::string("node"));
  StormEnsure(nodeStream.is_open());
  size_t numNodes{0}, dim{0};
  nodeStream >> numNodes >> dim;
  std::getline(nodeStream, line);
  for (size_t i{0}; i < numNodes; ++i) {
    NodeIndex nodeIndex{0};
    vec3_t nodeCoords(0.0);
    nodeStream >> nodeIndex >> nodeCoords.x >> nodeCoords.y;
    std::getline(nodeStream, line);
    StormEnsure(nodeIndex == insert_node(nodeCoords));
  }

  std::ifstream faceStream(path + std::string("edge"));
  StormEnsure(faceStream.is_open());
  size_t numFaces{0};
  faceStream >> numFaces;
  std::getline(faceStream, line);
  for (size_t i{0}; i < numFaces; ++i) {
    FaceIndex faceIndex{0};
    std::vector<NodeIndex> faceNodes(2);
    size_t faceMark{0};
    faceStream >> faceIndex >> faceNodes[0] >> faceNodes[1] >> faceMark;
    StormEnsure(faceIndex == insert_face({ShapeType::Segment, faceNodes},
                                         FaceMark(faceMark)));
    std::getline(faceStream, line);
  }

  std::ifstream cellStream(path + std::string("ele"));
  StormEnsure(cellStream.is_open());
  size_t numCells{0};
  cellStream >> numCells;
  std::getline(cellStream, line);
  for (size_t i{0}; i < numCells; ++i) {
    CellIndex cellIndex{0};
    std::vector<NodeIndex> cellNodes(3);
    cellStream >> cellIndex >> cellNodes[0] >> cellNodes[1] >> cellNodes[2];
    StormEnsure(cellIndex == insert_cell({ShapeType::Triangle, cellNodes}));
    std::getline(cellStream, line);
  }

  finalize();
  return true;

} // сMesh::read_from_triangle

#if 0
bool Mesh::read_from_tetgen(std::string const& path) {

  std::string line;

  std::ifstream nodeFile(path + std::string("node"));
  storm_ensure(nodeFile.is_open());
  size_t numNodes{0}, dim{0};
  nodeFile >> numNodes >> dim;
  storm_ensure(dim == 2);
  std::getline(nodeFile, line);
  for (size_t i{0}; i < numNodes; ++i) {
    size_t nodeIndex{0};
    vec3_t nodePos(0.0);
    nodeFile >> nodeIndex >> nodePos.x >> nodePos.y >> nodePos.z;
    storm_ensure(
      nodeIndex == insertNode(nodePos));
    std::getline(nodeFile, line);
  }

  std::ifstream faceFile(path + std::string("face"));
  storm_ensure(faceFile.is_open());
  size_t numFaces{0};
  faceFile >> numFaces;
  std::getline(faceFile, line);
  for (size_t i{0}; i < numFaces; ++i) {
    size_t faceIndex{0};
    std::vector<size_t> faceNodes(3, npos);
    size_t faceMark{0};
    faceFile >> faceIndex >> faceNodes[0] >> faceNodes[1] >> faceNodes[2] >> faceMark;
    storm_ensure(
      faceIndex == EmplaceFace({ShapeType::Triangle, faceNodes}, FaceMark{faceMark}));
    std::getline(faceFile, line);
  }

  std::ifstream cellFile(path + std::string("ele"));
  storm_ensure(cellFile.is_open());
  size_t numCells{0};
  cellFile >> numCells;
  std::getline(cellFile, line);
  for (size_t i{0}; i < numCells; ++i) {
    size_t cellIndex{0};
    std::vector<size_t> cellNodes(4, npos);
    cellFile >> cellIndex >> cellNodes[0] >> cellNodes[1] >> cellNodes[2] >> cellNodes[3];
    storm_ensure(
      cellIndex == EmplaceCell({ShapeType::Tetrahedron, cellNodes}));
    std::getline(cellFile, line);
  }

  finalize();
  return true;

} // сMesh::read_from_tetgen

bool Mesh::read_from_image(const char *path,
                           const std::map<Pixel, size_t>& mark_colors,
                           Pixel fluid_color,
                           vec2_t pixel_size) {

  Image2D image;
  image.load(path);

  Image2D nodesImage;
  nodesImage.init(image.width() + 1, image.height() + 1, Pixel(0, 0, 0, 0));

  size_t nodeIndex = 0;
  for (size_t y = 1; y < image.height() - 1; ++y) {
    for (size_t x = 1; x < image.width() - 1; ++x) {
      if (image(x, y).rgba != fluid_color.rgba) {
        continue;
      }

      vec2_t const cellCenterPos = pixel_size * vec2_t(real_t(x) - 0.5, real_t(y) - 0.5);

      // Insert or query the cell Nodes.
      size_t& swNodeIndex = nodesImage(x + 0, y + 0).rgba;
      if (swNodeIndex == 0) {
        swNodeIndex = nodeIndex++;
        vec3_t const nodePos(cellCenterPos + pixel_size * vec2_t(-0.5, -0.5), 0.0);
        storm_ensure(swNodeIndex == insertNode(nodePos));
      }
      size_t& seNodeIndex = nodesImage(x + 1, y + 0).rgba;
      if (seNodeIndex == 0) {
        seNodeIndex = nodeIndex++;
        vec3_t const nodePos(cellCenterPos + pixel_size * vec2_t(+0.5, -0.5), 0.0);
        storm_ensure(seNodeIndex == insertNode(nodePos));
      }
      size_t& neNodeIndex = nodesImage(x + 1, y + 1).rgba;
      if (neNodeIndex == 0) {
        neNodeIndex = nodeIndex++;
        vec3_t const nodePos(cellCenterPos + pixel_size * vec2_t(+0.5, +0.5), 0.0);
        storm_ensure(neNodeIndex == insertNode(nodePos));
      }
      size_t& nwNodeIndex = nodesImage(x + 0, y + 1).rgba;
      if (nwNodeIndex == 0) {
        nwNodeIndex = nodeIndex++;
        vec3_t const nodePos(cellCenterPos + pixel_size * vec2_t(-0.5, +0.5), 0.0);
        storm_ensure(nwNodeIndex == insertNode(nodePos));
      }

      // Insert the cell.
      EmplaceCell({ShapeType::Quadrangle,
                    {swNodeIndex, seNodeIndex, neNodeIndex, nwNodeIndex}});

      // Insert the boundary Faces.
      if (Pixel const sPixel = image(x, y - 1); sPixel.rgba != fluid_color.rgba) {
        FaceMark const faceMark{mark_colors.at(sPixel)};
        EmplaceFace({ShapeType::Segment, {swNodeIndex, seNodeIndex}}, faceMark);
      }
      if (Pixel const ePixel = image(x + 1, y); ePixel.rgba != fluid_color.rgba) {
        FaceMark const faceMark{mark_colors.at(ePixel)};
        EmplaceFace({ShapeType::Segment, {seNodeIndex, neNodeIndex}}, faceMark);
      }
      if (Pixel const nPixel = image(x, y + 1); nPixel.rgba != fluid_color.rgba) {
        FaceMark const faceMark{mark_colors.at(nPixel)};
        EmplaceFace({ShapeType::Segment, {neNodeIndex, nwNodeIndex}}, faceMark);
      }
      if (Pixel const wPixel = image(x - 1, y); wPixel.rgba != fluid_color.rgba) {
        FaceMark const faceMark{mark_colors.at(wPixel)};
        EmplaceFace({ShapeType::Segment, {nwNodeIndex, swNodeIndex}}, faceMark);
      }
    }
  }

  finalize();
  return true;

} // Mesh::read_from_image
#endif

void Mesh::save_vtk(const char* path,
                    const std::vector<sFieldDesc>& fields) const {
  std::ofstream file(path);
  file << std::setprecision(std::numeric_limits<real_t>::digits10 + 1);
  file << "# vtk DataFile Version 2.0" << std::endl;
  file << "# Generated by Feathers/StormRuler/Mesh2VTK" << std::endl;
  file << "ASCII" << std::endl;
  file << "DATASET UNSTRUCTURED_GRID" << std::endl;

  file << "POINTS " << nodes().size() << " double" << std::endl;
  ranges::for_each(node_views(*this), [&](NodeView node) {
    const vec3_t& pos = node.coords();
    file << pos.x << " " << pos.y << " " << pos.z << std::endl;
  });
  file << std::endl;

  size_t const sumNumCellAdjNodes =
      ForEachSum(int_cell_views(*this), size_t(0), [](CellView cell) {
        return cell.adjacent_nodes().size() + 1;
      });
  file << "CELLS " << cells({}).size() << " " << sumNumCellAdjNodes
       << std::endl;
  ranges::for_each(int_cell_views(*this), [&](CellView cell) {
    file << cell.adjacent_nodes().size() << " ";
    cell.for_each_node([&](size_t node_index) { file << node_index << " "; });
    file << std::endl;
  });
  file << std::endl;

  file << "CELL_TYPES " << cells({}).size() << std::endl;
  ranges::for_each(int_cell_views(*this), [&](CellView cell) {
    static const std::map<ShapeType, const char*> shapes = {
        {ShapeType::Node, "1"},         {ShapeType::Segment, "2"},
        {ShapeType::Triangle, "5"},     {ShapeType::Quadrangle, "9"},
        {ShapeType::Tetrahedron, "10"}, {ShapeType::Pyramid, "14"},
        {ShapeType::Pentahedron, "13"}, {ShapeType::Hexahedron, "12"}};
    file << shapes.at(cell.shape_type()) << std::endl;
  });
  file << std::endl;

  file << "CELL_DATA " << cells({}).size() << std::endl;
  for (const sFieldDesc& field : fields) {
    file << "SCALARS " << field.name << " double 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    ranges::for_each(int_cell_views(*this), [&](CellView cell) {
      file << (*field.scalar)[cell][field.var_index] << std::endl;
    });
  }
  file << std::endl;
} // Mesh::save_vtk

} // namespace Storm
