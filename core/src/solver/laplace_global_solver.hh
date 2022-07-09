#pragma once

#include "laplace_solver.hh"

namespace plksim {
namespace solver {

template <int dim>
class LaplaceGlobalSolver : public LaplaceSolver<dim> {
public:
  LaplaceGlobalSolver(Triangulation<dim>& coarse_grid, const FiniteElement<dim>& fe,
                      const Quadrature<dim>& quadrature, const Quadrature<dim - 1>& face_quadrature,
                      const Function<dim>& rhs_function, const Function<dim>& boundary_values);

  virtual void refine_grid() override;
};

template <int dim>
LaplaceGlobalSolver<dim>::LaplaceGlobalSolver(Triangulation<dim>& coarse_grid,
                                              const FiniteElement<dim>& fe,
                                              const Quadrature<dim>& quadrature,
                                              const Quadrature<dim - 1>& face_quadrature,
                                              const Function<dim>& rhs_function,
                                              const Function<dim>& boundary_values)
    : SolverBase<dim>(coarse_grid), LaplaceSolver<dim>(coarse_grid, fe, quadrature, face_quadrature,
                                                       rhs_function, boundary_values) {
}

template <int dim>
void LaplaceGlobalSolver<dim>::refine_grid() {
  this->triangulation->refine_global(1);
}

} // namespace solver
} // namespace plksim