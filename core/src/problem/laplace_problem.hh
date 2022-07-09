#pragma once

#include "data/sample.hh"
#include "data/set_up.hh"
#include "evaluation/point_value_evaluation.hh"
#include "problem_base.hh"
#include "problem_log.hh"
#include "solver/laplace_global_solver.hh"
#include "solver/laplace_kelly_solver.hh"
#include "util/output_format_util.hh"

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/tria.h>

#include <iostream>
#include <memory>

namespace plksim {
namespace problem {

using namespace dealii;

template <int dim>
class LaplaceProblem : public ProblemBase {
public:
  LaplaceProblem() = default;
  void compute() override;
};

template <int dim>
void LaplaceProblem<dim>::compute() {
  write_log("STARTED", {}, log_stream);

  unsigned int primal_fe_degree = 1;
  unsigned int dual_fe_degree = 2;
  unsigned int max_degrees_of_freedom = 20000;

  auto data = std::make_unique<data::SetUp<data::Sample<dim>, dim>>();

  Triangulation<dim> triangulation(Triangulation<dim>::smoothing_on_refinement);
  data->create_coarse_grid(triangulation);

  const Point<dim> evaluation_point(0.75, 0.75);
  evaluation::PointValueEvaluation<dim> postprocessor_eval(evaluation_point);

  std::list<evaluation::EvaluationBase<dim>*> evaluator_list;
  evaluator_list.push_back(&postprocessor_eval);

  const FE_Q<dim> primal_fe(primal_fe_degree);
  const FE_Q<dim> dual_fe(dual_fe_degree);
  const QGauss<dim> quadrature(dual_fe_degree + 1);
  const QGauss<dim - 1> face_quadrature(dual_fe_degree + 1);

  std::unique_ptr<solver::SolverBase<dim>> solver;
  if (refine_criterion == "kelly") {
    solver = std::make_unique<solver::LaplaceKellySolver<dim>>( //
        triangulation,                                          //
        primal_fe,                                              //
        quadrature,                                             //
        face_quadrature,                                        //
        data->get_right_hand_side(),                            //
        data->get_boundary_values());
  } else {
    solver = std::make_unique<solver::LaplaceGlobalSolver<dim>>( //
        triangulation,                                           //
        primal_fe,                                               //
        quadrature,                                              //
        face_quadrature,                                         //
        data->get_right_hand_side(),                             //
        data->get_boundary_values());
  }

  for (unsigned int step = 0; true; ++step) {
    solver->set_refinement_cycle(step);
    solver->solve();

    for (const auto& evaluator : evaluator_list) {
      evaluator->set_refinement_cycle(step);
      solver->postprocess(*evaluator);
    }

    if (solver->n_dofs() < max_degrees_of_freedom) {
      solver->refine_grid();
    } else {
      break;
    }
  }

  if (output_solution_stream != nullptr) {
    solver->output_solution(*output_solution_stream, output_solution_format);
  }

  write_log("FINISHED", {}, log_stream);
}

} // namespace problem
} // namespace plksim