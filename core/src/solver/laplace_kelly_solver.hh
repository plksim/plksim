#pragma once

#include "laplace_solver.hh"

#include <deal.II/grid/grid_refinement.h>

namespace plksim {
namespace solver {

template <int dim>
class LaplaceKellySolver : public LaplaceSolver<dim> {
public:
  LaplaceKellySolver(Triangulation<dim>& coarse_grid, const FiniteElement<dim>& fe,
                     const Quadrature<dim>& quadrature, const Quadrature<dim - 1>& face_quadrature,
                     const Function<dim>& rhs_function, const Function<dim>& boundary_values);

  virtual void refine_grid() override;
};

template <int dim>
LaplaceKellySolver<dim>::LaplaceKellySolver(Triangulation<dim>& coarse_grid,
                                            const FiniteElement<dim>& fe,
                                            const Quadrature<dim>& quadrature,
                                            const Quadrature<dim - 1>& face_quadrature,
                                            const Function<dim>& rhs_function,
                                            const Function<dim>& boundary_values)
    : SolverBase<dim>(coarse_grid), LaplaceSolver<dim>(coarse_grid, fe, quadrature, face_quadrature,
                                                       rhs_function, boundary_values) {
}

template <int dim>
void LaplaceKellySolver<dim>::refine_grid() {
  Vector<float> estimated_error_per_cell(this->triangulation->n_active_cells());
  KellyErrorEstimator<dim>::estimate(this->dof_handler, QGauss<dim - 1>(this->fe->degree + 1),
                                     std::map<types::boundary_id, const Function<dim>*>(),
                                     this->solution, estimated_error_per_cell);
  GridRefinement::refine_and_coarsen_fixed_number(*this->triangulation, estimated_error_per_cell,
                                                  0.3, 0.03);
  this->triangulation->execute_coarsening_and_refinement();
}

} // namespace solver
} // namespace plksim