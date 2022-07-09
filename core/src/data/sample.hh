#pragma once

#include <deal.II/base/function.h>
#include <deal.II/grid/tria.h>

#include <vector>

namespace plksim {
namespace data {

using namespace dealii;

template <int dim>
struct Sample {
  using BoundaryValues = Functions::ZeroFunction<dim>;

  class RightHandSide : public Functions::ConstantFunction<dim> {
  public:
    RightHandSide() : Functions::ConstantFunction<dim>(1.) {
    }
  };

  static void create_coarse_grid(Triangulation<dim>& coarse_grid);
};

template <>
void Sample<2>::create_coarse_grid(Triangulation<2>& coarse_grid) {
  const unsigned int dim = 2;

  const std::vector<Point<2>> vertices = {
      {-1.0, -1.0}, {-0.5, -1.0}, {+0.0, -1.0}, {+0.5, -1.0}, {+1.0, -1.0}, {-1.0, -0.5},
      {-0.5, -0.5}, {+0.0, -0.5}, {+0.5, -0.5}, {+1.0, -0.5}, {-1.0, +0.0}, {-0.5, +0.0},
      {+0.5, +0.0}, {+1.0, +0.0}, {-1.0, +0.5}, {-0.5, +0.5}, {+0.0, +0.5}, {+0.5, +0.5},
      {+1.0, +0.5}, {-1.0, +1.0}, {-0.5, +1.0}, {+0.0, +1.0}, {+0.5, +1.0}, {+1.0, +1.0}};

  const std::vector<std::array<int, GeometryInfo<dim>::vertices_per_cell>> cell_vertices = {
      {{0, 1, 5, 6}},     {{1, 2, 6, 7}},     {{2, 3, 7, 8}},     {{3, 4, 8, 9}},
      {{5, 6, 10, 11}},   {{8, 9, 12, 13}},   {{10, 11, 14, 15}}, {{12, 13, 17, 18}},
      {{14, 15, 19, 20}}, {{15, 16, 20, 21}}, {{16, 17, 21, 22}}, {{17, 18, 22, 23}}};

  const unsigned int n_cells = cell_vertices.size();

  std::vector<CellData<dim>> cells(n_cells, CellData<dim>());
  for (unsigned int i = 0; i < n_cells; ++i) {
    for (unsigned int j = 0; j < cell_vertices[i].size(); ++j) {
      cells[i].vertices[j] = cell_vertices[i][j];
    }

    cells[i].material_id = 0;
  }

  coarse_grid.create_triangulation(vertices, cells, SubCellData());
  coarse_grid.refine_global(1);
}

} // namespace data
} // namespace plksim