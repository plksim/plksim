#pragma once

#include "problem_base.hh"
#include "problem_log.hh"
#include "utils/output_format_utils.hh"

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

namespace plksim {
namespace problems {

using namespace dealii;

template <int dim>
class LaplaceProblem : public ProblemBase {
public:
  LaplaceProblem();
  ~LaplaceProblem();
  void compute(std::ostream& out, const io::OutputFormat out_format, std::ostream& log) override;

private:
  void make_grid();
  void setup_system();
  void assemble_system();
  void solve();
  void refine_grid();
  void output(std::ostream& out, const io::OutputFormat out_format) const;

  double coefficient(const Point<dim>& p) const;

  Triangulation<dim> triangulation;
  FE_Q<dim> fe;
  DoFHandler<dim> dof_handler;

  AffineConstraints<double> constraints;

  SparseMatrix<double> system_matrix;
  SparsityPattern sparsity_pattern;

  Vector<double> solution;
  Vector<double> system_rhs;
};

template <int dim>
LaplaceProblem<dim>::LaplaceProblem() : fe(2), dof_handler(triangulation){};

template <int dim>
LaplaceProblem<dim>::~LaplaceProblem(){};

template <int dim>
void LaplaceProblem<dim>::compute(std::ostream& out_stream, const io::OutputFormat out_format,
                                  std::ostream& log_stream) {
  write_log("STARTED", {}, log_stream);

  make_grid();

  for (int cycle = 0; cycle < 8; ++cycle) {
    if (cycle > 0) {
      refine_grid();
    }

    setup_system();
    assemble_system();
    solve();
  }

  output(out_stream, out_format);

  write_log("FINISHED", {}, log_stream);
};

template <int dim>
void LaplaceProblem<dim>::make_grid() {
  GridGenerator::hyper_ball(triangulation);
  triangulation.refine_global(1);
};

template <int dim>
void LaplaceProblem<dim>::setup_system() {
  dof_handler.distribute_dofs(fe);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  VectorTools::interpolate_boundary_values(dof_handler, 0, Functions::ZeroFunction<dim>(), constraints);

  constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);

  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);
};

template <int dim>
void LaplaceProblem<dim>::assemble_system() {
  QGauss<dim> quadrature_formula(fe.degree + 1);

  FEValues<dim> fe_values(fe, quadrature_formula,
                          update_values | update_gradients | update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto& cell : dof_handler.active_cell_iterators()) {
    fe_values.reinit(cell);

    cell_matrix = 0;
    cell_rhs = 0;

    for (const unsigned int q_index : fe_values.quadrature_point_indices()) {
      const double current_coefficient = coefficient(fe_values.quadrature_point(q_index));
      for (const unsigned int i : fe_values.dof_indices()) {
        for (const unsigned int j : fe_values.dof_indices())
          cell_matrix(i, j) += (current_coefficient *              // a(xQ)
                                fe_values.shape_grad(i, q_index) * // grad phi_i(xQ)
                                fe_values.shape_grad(j, q_index) * // grad phi_j(xQ)
                                fe_values.JxW(q_index));           // dx

        cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(xQ)
                        1.0 *                               // f(x)
                        fe_values.JxW(q_index));            // dx
      }
    }

    cell->get_dof_indices(local_dof_indices);
    constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
  }
};

template <int dim>
void LaplaceProblem<dim>::solve() {
  SolverControl sc(1000, 1e-12);
  SolverCG<Vector<double>> solver(sc);

  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  solver.solve(system_matrix, solution, system_rhs, preconditioner);

  constraints.distribute(solution);
};

template <int dim>
void LaplaceProblem<dim>::refine_grid() {
  Vector<float> estimated_error_perCell(triangulation.n_active_cells());

  KellyErrorEstimator<dim>::estimate(dof_handler, QGauss<dim - 1>(fe.degree + 1), {}, solution,
                                     estimated_error_perCell);

  GridRefinement::refine_and_coarsen_fixed_number(triangulation, estimated_error_perCell, 0.3, 0.03);

  triangulation.execute_coarsening_and_refinement();
};

template <int dim>
void LaplaceProblem<dim>::output(std::ostream& out_stream, const io::OutputFormat out_format) const {
  DataOut<dim> dataOut;

  dataOut.attach_dof_handler(dof_handler);
  dataOut.add_data_vector(solution, "result");
  dataOut.build_patches();

  auto format = utils::output_format_to_dealii(out_format);
  dataOut.write(out_stream, format);
};

template <int dim>
double LaplaceProblem<dim>::coefficient(const Point<dim>& p) const {
  if (p.square() < 0.5 * 0.5)
    return 10;
  else
    return 1;
};

} // namespace problems
} // namespace plksim