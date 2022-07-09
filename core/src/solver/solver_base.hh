#pragma once

#include "evaluation/evaluation_base.hh"

#include "io/output_format.hh"

#include <deal.II/base/smartpointer.h>
#include <deal.II/base/types.h>
#include <deal.II/grid/tria.h>

#include <iostream>

namespace plksim {
namespace solver {

using namespace dealii;

template <int dim>
class SolverBase {
public:
  SolverBase(Triangulation<dim>& coarse_grid);
  virtual ~SolverBase() = default;

  virtual void solve() = 0;
  virtual void postprocess(const evaluation::EvaluationBase<dim>& postprocessor) const = 0;
  virtual void refine_grid() = 0;
  virtual unsigned int n_dofs() const = 0;
  virtual void output_solution(std::ostream& out_stream,
                               const io::OutputFormat out_format) const = 0;

  virtual void set_refinement_cycle(const unsigned int v);

protected:
  const SmartPointer<Triangulation<dim>> triangulation;

  unsigned int refinement_cycle;
};

template <int dim>
SolverBase<dim>::SolverBase(Triangulation<dim>& coarse_grid)
    : triangulation(&coarse_grid), refinement_cycle(numbers::invalid_unsigned_int) {
}

template <int dim>
void SolverBase<dim>::set_refinement_cycle(const unsigned int v) {
  refinement_cycle = v;
}

} // namespace solver
} // namespace plksim