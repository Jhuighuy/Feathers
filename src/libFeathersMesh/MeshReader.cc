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

bool Mesh::ReadFromTriangle(std::string const& path) {

  std::string line;

  std::ifstream nodeFile(path + std::string("node"));
  StormEnsure(nodeFile.is_open());
  size_t numNodes{0}, dim{0};
  nodeFile >> numNodes >> dim;
  std::getline(nodeFile, line);
  for (size_t i{0}; i < numNodes; ++i) {
    size_t nodeIndex{0};
    vec3_t nodePos(0.0);
    nodeFile >> nodeIndex >> nodePos.x >> nodePos.y;
    std::getline(nodeFile, line);
    StormEnsure(nodeIndex == EmplaceNode(nodePos));
  }

  std::ifstream faceFile(path + std::string("edge"));
  StormEnsure(faceFile.is_open());
  size_t numFaces{0};
  faceFile >> numFaces;
  std::getline(faceFile, line);
  for (size_t i{0}; i < numFaces; ++i) {
    size_t faceIndex{0};
    std::vector<size_t> faceNodes(2, npos);
    size_t faceMark{0};
    faceFile >> faceIndex >> faceNodes[0] >> faceNodes[1] >> faceMark;
    StormEnsure(
      faceIndex == EmplaceFace({ShapeType::Segment2, faceNodes}, FaceMark(faceMark)));
    std::getline(faceFile, line);
  }

  std::ifstream cellFile(path + std::string("ele"));
  FEATHERS_ENSURE(cellFile.is_open());
  size_t numCells{0};
  cellFile >> numCells;
  std::getline(cellFile, line);
  for (size_t i{0}; i < numCells; ++i) {
    size_t cellIndex{0};
    std::vector<size_t> cellNodes(3, npos);
    cellFile >> cellIndex >> cellNodes[0] >> cellNodes[1] >> cellNodes[2];
    StormEnsure(
      cellIndex == EmplaceCell({ShapeType::Triangle3, cellNodes}));
    std::getline(cellFile, line);
  }

  finalize();
  return true;

} // сMesh::ReadFromTriangle

bool Mesh::ReadFromTetgen(std::string const& path) {

  std::string line;

  std::ifstream nodeFile(path + std::string("node"));
  StormEnsure(nodeFile.is_open());
  size_t numNodes{0}, dim{0};
  nodeFile >> numNodes >> dim;
  StormEnsure(dim == 2);
  std::getline(nodeFile, line);
  for (size_t i{0}; i < numNodes; ++i) {
    size_t nodeIndex{0};
    vec3_t nodePos(0.0);
    nodeFile >> nodeIndex >> nodePos.x >> nodePos.y >> nodePos.z;
    StormEnsure(
      nodeIndex == EmplaceNode(nodePos));
    std::getline(nodeFile, line);
  }

  std::ifstream faceFile(path + std::string("face"));
  StormEnsure(faceFile.is_open());
  size_t numFaces{0};
  faceFile >> numFaces;
  std::getline(faceFile, line);
  for (size_t i{0}; i < numFaces; ++i) {
    size_t faceIndex{0};
    std::vector<size_t> faceNodes(3, npos);
    size_t faceMark{0};
    faceFile >> faceIndex >> faceNodes[0] >> faceNodes[1] >> faceNodes[2] >> faceMark;
    StormEnsure(
      faceIndex == EmplaceFace({ShapeType::Triangle3, faceNodes}, FaceMark{faceMark}));
    std::getline(faceFile, line);
  }

  std::ifstream cellFile(path + std::string("ele"));
  StormEnsure(cellFile.is_open());
  size_t numCells{0};
  cellFile >> numCells;
  std::getline(cellFile, line);
  for (size_t i{0}; i < numCells; ++i) {
    size_t cellIndex{0};
    std::vector<size_t> cellNodes(4, npos);
    cellFile >> cellIndex >> cellNodes[0] >> cellNodes[1] >> cellNodes[2] >> cellNodes[3];
    StormEnsure(
      cellIndex == EmplaceCell({ShapeType::Tetrahedron4, cellNodes}));
    std::getline(cellFile, line);
  }

  finalize();
  return true;

} // сMesh::ReadFromTetgen

bool Mesh::ReadFromImage(const char *path,
                         const std::map<Pixel, size_t>& markColors,
                         Pixel fluidColor,
                         vec2_t pixelSize) {

  Image2D image;
  image.load(path);

  Image2D nodesImage;
  nodesImage.init(image.width() + 1, image.height() + 1, Pixel(0, 0, 0, 0));

  size_t nodeIndex = 0;
  for (size_t y = 1; y < image.height() - 1; ++y) {
    for (size_t x = 1; x < image.width() - 1; ++x) {
      if (image(x, y).rgba != fluidColor.rgba) {
        continue;
      }

      vec2_t const cellCenterPos = pixelSize * vec2_t(real_t(x) - 0.5, real_t(y) - 0.5);

      // Insert or query the cell nodes.
      size_t& swNodeIndex = nodesImage(x + 0, y + 0).rgba;
      if (swNodeIndex == 0) {
        swNodeIndex = nodeIndex++;
        vec3_t const nodePos(cellCenterPos + pixelSize * vec2_t(-0.5, -0.5), 0.0);
        StormEnsure(swNodeIndex == EmplaceNode(nodePos));
      }
      size_t& seNodeIndex = nodesImage(x + 1, y + 0).rgba;
      if (seNodeIndex == 0) {
        seNodeIndex = nodeIndex++;
        vec3_t const nodePos(cellCenterPos + pixelSize * vec2_t(+0.5, -0.5), 0.0);
        StormEnsure(seNodeIndex == EmplaceNode(nodePos));
      }
      size_t& neNodeIndex = nodesImage(x + 1, y + 1).rgba;
      if (neNodeIndex == 0) {
        neNodeIndex = nodeIndex++;
        vec3_t const nodePos(cellCenterPos + pixelSize * vec2_t(+0.5, +0.5), 0.0);
        StormEnsure(neNodeIndex == EmplaceNode(nodePos));
      }
      size_t& nwNodeIndex = nodesImage(x + 0, y + 1).rgba;
      if (nwNodeIndex == 0) {
        nwNodeIndex = nodeIndex++;
        vec3_t const nodePos(cellCenterPos + pixelSize * vec2_t(-0.5, +0.5), 0.0);
        StormEnsure(nwNodeIndex == EmplaceNode(nodePos));
      }

      // Insert the cell.
      EmplaceCell({ShapeType::Quadrangle4,
                    {swNodeIndex, seNodeIndex, neNodeIndex, nwNodeIndex}});

      // Insert the boundary faces.
      if (Pixel const sPixel = image(x, y - 1); sPixel.rgba != fluidColor.rgba) {
        FaceMark const faceMark{markColors.at(sPixel)};
        EmplaceFace({ShapeType::Segment2, {swNodeIndex, seNodeIndex}}, faceMark);
      }
      if (Pixel const ePixel = image(x + 1, y); ePixel.rgba != fluidColor.rgba) {
        FaceMark const faceMark{markColors.at(ePixel)};
        EmplaceFace({ShapeType::Segment2, {seNodeIndex, neNodeIndex}}, faceMark);
      }
      if (Pixel const nPixel = image(x, y + 1); nPixel.rgba != fluidColor.rgba) {
        FaceMark const faceMark{markColors.at(nPixel)};
        EmplaceFace({ShapeType::Segment2, {neNodeIndex, nwNodeIndex}}, faceMark);
      }
      if (Pixel const wPixel = image(x - 1, y); wPixel.rgba != fluidColor.rgba) {
        FaceMark const faceMark{markColors.at(wPixel)};
        EmplaceFace({ShapeType::Segment2, {nwNodeIndex, swNodeIndex}}, faceMark);
      }
    }
  }

  finalize();
  return true;

} // Mesh::ReadFromImage

void Mesh::save_vtk(const char* path,
                    const std::vector<sFieldDesc>& fields) const {
  std::ofstream file(path);
  file << std::setprecision(std::numeric_limits<real_t>::digits10 + 1);
  file << "# vtk DataFile Version 2.0" << std::endl;
  file << "# Generated by Feathers/StormRuler/Mesh2VTK" << std::endl;
  file << "ASCII" << std::endl;
  file << "DATASET UNSTRUCTURED_GRID" << std::endl;

  file << "POINTS " << NumNodes() << " double" << std::endl;
  ranges::for_each(NodeRefs(*this), [&](NodeRef node) {
    const vec3_t& pos = node.Pos();
    file << pos.x << " " << pos.y << " " << pos.z << std::endl;
  });
  file << std::endl;

  size_t const sumNumCellAdjNodes =
    ForEachSum(InteriorCellRefs(*this), size_t(0), [](CellRef cell) {
      return cell.NumNodes() + 1;
    });
  file << "CELLS " << Cells({}).size() << " " << sumNumCellAdjNodes << std::endl;
  ranges::for_each(InteriorCellRefs(*this), [&](CellRef cell) {
    file << cell.NumNodes() << " ";
    cell.ForEachNode([&](size_t node_index) {
      file << node_index << " ";
    });
    file << std::endl;
  });
  file << std::endl;

  file << "CELL_TYPES " << Cells({}).size() << std::endl;
  ranges::for_each(InteriorCellRefs(*this), [&](CellRef cell) {
    static const std::map<ShapeType, const char*> shapes = {
      { ShapeType::Node, "1" }, { ShapeType::Segment2, "2" },
      { ShapeType::Triangle3, "5" }, { ShapeType::Quadrangle4, "9" },
      { ShapeType::Tetrahedron4, "10" }, { ShapeType::Pyramid5, "14" },
      { ShapeType::Pentahedron6, "13" }, { ShapeType::Hexahedron8, "12" }
    };
    file << shapes.at(cell.Shape()) << std::endl;
  });
  file << std::endl;

  file << "CELL_DATA " << Cells({}).size() << std::endl;
  for (const sFieldDesc& field : fields) {
    file << "SCALARS " << field.name << " double 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    ranges::for_each(InteriorCellRefs(*this), [&](CellRef cell) {
      file << (*field.scalar)[cell][field.var_index] << std::endl;
    });
  }
  file << std::endl;
} // Mesh::save_vtk

} // namespace feathers
