#pragma once

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/vector.h>

namespace plksim {
namespace evaluation {

using namespace dealii;

template <int dim>
class EvaluationBase {
public:
  virtual ~EvaluationBase() = default;

  virtual void operator()(const DoFHandler<dim>& dof_handler,
                          const Vector<double>& solution) const = 0;

  void set_refinement_cycle(const unsigned int v);

protected:
  unsigned int refinement_cycle;
};

template <int dim>
void EvaluationBase<dim>::set_refinement_cycle(const unsigned int v) {
  refinement_cycle = v;
}

} // namespace evaluation
} // namespace plksim