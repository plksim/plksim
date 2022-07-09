#pragma once

#include "evaluation_base.hh"

#include <deal.II/base/point.h>

#include <exception>
#include <sstream>

namespace plksim {
namespace evaluation {

template <int dim>
class PointValueEvaluation : public EvaluationBase<dim> {
public:
  PointValueEvaluation(const Point<dim>& evaluation_point);

  virtual void operator()(const DoFHandler<dim>& dof_handler,
                          const Vector<double>& solution) const override;

private:
  const Point<dim> evaluation_point;
};

template <int dim>
PointValueEvaluation<dim>::PointValueEvaluation(const Point<dim>& evaluation_point)
    : evaluation_point(evaluation_point) {
}

template <int dim>
void PointValueEvaluation<dim>::operator()(const DoFHandler<dim>& dof_handler,
                                           const Vector<double>& /*solution*/) const {
  bool evaluation_point_found = false;

  for (const auto& cell : dof_handler.active_cell_iterators()) {
    if (!evaluation_point_found) {
      for (const auto vertex : cell->vertex_indices()) {
        if (cell->vertex(vertex).distance(evaluation_point) < cell->diameter() * 1e-8) {
          evaluation_point_found = true;
          break;
        }
      }
    }
  }

  if (!evaluation_point_found) {
    std::ostringstream ss;
    ss << "NOT_FOUND: evaluation point " << evaluation_point << " is not found among the vertices";
    throw std::runtime_error(ss.str());
  }
}

} // namespace evaluation
} // namespace plksim