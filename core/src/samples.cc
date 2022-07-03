#include "samples.hh"

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <sstream>

using namespace dealii;

namespace plksim {

std::string sample_mesh_svg() {
  Triangulation<2> triangulation;

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(4);

  std::ostringstream out;
  GridOut grid_out;
  grid_out.write_svg(triangulation, out);

  return out.str();
};

} // namespace plksim