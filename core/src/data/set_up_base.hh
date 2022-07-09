#pragma once

#include <deal.II/base/function.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/grid/tria.h>

namespace plksim {
namespace data {

using namespace dealii;

template <int dim>
struct SetUpBase : public Subscriptor {
  virtual const Function<dim>& get_boundary_values() const = 0;
  virtual const Function<dim>& get_right_hand_side() const = 0;

  virtual void create_coarse_grid(Triangulation<dim>& coarse_grid) const = 0;
};

} // namespace data
} // namespace plksim