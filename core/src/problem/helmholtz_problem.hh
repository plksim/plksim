#pragma once

#include "problem_base.hh"
#include "problem_log.hh"
#include "util/output_format_util.hh"

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
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
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <array>
#include <sstream>

namespace plksim {
namespace problem {

using namespace dealii;

template <int dim>
class SolutionBase {
protected:
  static const std::array<Point<dim>, 3> source_centers;
  static const double width;
};

template <>
const std::array<Point<1>, 3> SolutionBase<1>::source_centers = {
    {Point<1>(-1.0 / 3.0), Point<1>(0.0), Point<1>(+1.0 / 3.0)}};

template <>
const std::array<Point<2>, 3> SolutionBase<2>::source_centers = {
    {Point<2>(-0.5, +0.5), Point<2>(-0.5, -0.5), Point<2>(+0.5, -0.5)}};

template <int dim>
const double SolutionBase<dim>::width = 0.125;

template <int dim>
class Solution : public Function<dim>, protected SolutionBase<dim> {
public:
  virtual double value(const Point<dim>& p, const unsigned int component = 0) const override;

  virtual Tensor<1, dim> gradient(const Point<dim>& p,
                                  const unsigned int component = 0) const override;
};

template <int dim>
double Solution<dim>::value(const Point<dim>& p, const unsigned int) const {
  double return_value = 0;
  for (const auto& center : this->source_centers) {
    const Tensor<1, dim> x_minus_xi = p - center;
    return_value += std::exp(-x_minus_xi.norm_square() / (this->width * this->width));
  }

  return return_value;
}

template <int dim>
Tensor<1, dim> Solution<dim>::gradient(const Point<dim>& p, const unsigned int) const {
  Tensor<1, dim> return_value;

  for (const auto& center : this->source_centers) {
    const Tensor<1, dim> x_minus_xi = p - center;

    return_value +=
        (-2. / (this->width * this->width) *
         std::exp(-x_minus_xi.norm_square() / (this->width * this->width)) * x_minus_xi);
  }

  return return_value;
}

template <int dim>
class HelmholtzRhs : public Function<dim>, protected SolutionBase<dim> {
public:
  virtual double value(const Point<dim>& p, const unsigned int component = 0) const override;
};

template <int dim>
double HelmholtzRhs<dim>::value(const Point<dim>& p, const unsigned int) const {
  double return_value = 0;
  for (const auto& center : this->source_centers) {
    const Tensor<1, dim> x_minus_xi = p - center;

    return_value += ((2. * dim - 4. * x_minus_xi.norm_square() / (this->width * this->width)) /
                     (this->width * this->width) *
                     std::exp(-x_minus_xi.norm_square() / (this->width * this->width)));
    return_value += std::exp(-x_minus_xi.norm_square() / (this->width * this->width));
  }

  return return_value;
}

template <int dim>
class HelmholtzProblem : public ProblemBase {
public:
  HelmholtzProblem();
  void compute() override;

private:
  void make_grid();
  void setup_system();
  void assemble_system();
  void solve();
  void refine_grid();
  void process_solution(const unsigned int cycle);
  void output() const;

  Triangulation<dim> triangulation;
  FE_Q<dim> fe;
  DoFHandler<dim> dof_handler;

  AffineConstraints<double> constraints;

  SparseMatrix<double> system_matrix;
  SparsityPattern sparsity_pattern;

  Vector<double> solution;
  Vector<double> system_rhs;

  ConvergenceTable convergence_table;
};

template <int dim>
HelmholtzProblem<dim>::HelmholtzProblem() : fe(2), dof_handler(triangulation){};

template <int dim>
void HelmholtzProblem<dim>::compute() {
  write_log("STARTED", {}, log_stream);

  make_grid();

  for (int cycle = 0; cycle < 8; ++cycle) {
    if (cycle > 0) {
      refine_grid();
    }

    setup_system();
    assemble_system();
    solve();

    process_solution(cycle);
  }

  output();

  write_log("FINISHED", {}, log_stream);
};

template <int dim>
void HelmholtzProblem<dim>::make_grid() {
  GridGenerator::hyper_cube(triangulation, -1.0, 1.0);
  triangulation.refine_global(3);

  for (const auto& cell : triangulation.cell_iterators())
    for (const auto& face : cell->face_iterators()) {
      const auto center = face->center();
      if ((std::fabs(center(0) - (-1.0)) < 1e-12) || (std::fabs(center(1) - (-1.0)) < 1e-12))
        face->set_boundary_id(1);
    }
};

template <int dim>
void HelmholtzProblem<dim>::setup_system() {
  dof_handler.distribute_dofs(fe);
  DoFRenumbering::Cuthill_McKee(dof_handler);

  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  constraints.condense(dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
};

template <int dim>
void HelmholtzProblem<dim>::assemble_system() {
  QGauss<dim> quadrature_formula(fe.degree + 1);
  QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);

  const unsigned int n_q_points = quadrature_formula.size();
  const unsigned int n_face_q_points = face_quadrature_formula.size();

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  FEValues<dim> fe_values(fe, quadrature_formula,
                          update_values | update_gradients | update_quadrature_points |
                              update_JxW_values);

  FEFaceValues<dim> fe_face_values(fe, face_quadrature_formula,
                                   update_values | update_normal_vectors |
                                       update_quadrature_points | update_JxW_values);

  const HelmholtzRhs<dim> right_hand_side;
  std::vector<double> rhs_values(n_q_points);

  Solution<dim> exact_solution;

  for (const auto& cell : dof_handler.active_cell_iterators()) {
    cell_matrix = 0;
    cell_rhs = 0;

    fe_values.reinit(cell);

    right_hand_side.value_list(fe_values.get_quadrature_points(), rhs_values);

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      for (unsigned int i = 0; i < dofs_per_cell; ++i) {
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          cell_matrix(i, j) += ((fe_values.shape_grad(i, q_point) *         // grad phi_i(x_q)
                                     fe_values.shape_grad(j, q_point)       // grad phi_j(x_q)
                                 + fe_values.shape_value(i, q_point) *      // phi_i(x_q)
                                       fe_values.shape_value(j, q_point)) * // phi_j(x_q)
                                fe_values.JxW(q_point));                    // dx

        cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q)
                        rhs_values[q_point] *               // f(x_q)
                        fe_values.JxW(q_point));            // dx
      }

    for (const auto& face : cell->face_iterators())
      if (face->at_boundary() && (face->boundary_id() == 1)) {
        fe_face_values.reinit(cell, face);

        for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point) {
          const double neumann_value =
              (exact_solution.gradient(fe_face_values.quadrature_point(q_point)) *
               fe_face_values.normal_vector(q_point));

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            cell_rhs(i) += (fe_face_values.shape_value(i, q_point) * // phi_i(x_q)
                            neumann_value *                          // g(x_q)
                            fe_face_values.JxW(q_point));            // dx
        }
      }

    cell->get_dof_indices(local_dof_indices);

    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
        system_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));

      system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
  }

  constraints.condense(system_matrix);
  constraints.condense(system_rhs);

  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler, 0, Solution<dim>(), boundary_values);
  MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution, system_rhs);
};

template <int dim>
void HelmholtzProblem<dim>::solve() {
  SolverControl sc(1000, 1e-12);
  SolverCG<Vector<double>> solver(sc);

  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  solver.solve(system_matrix, solution, system_rhs, preconditioner);

  constraints.distribute(solution);
};

template <int dim>
void HelmholtzProblem<dim>::refine_grid() {
  Vector<float> estimatedErrorPerCell(triangulation.n_active_cells());

  KellyErrorEstimator<dim>::estimate(dof_handler, QGauss<dim - 1>(fe.degree + 1),
                                     std::map<types::boundary_id, const Function<dim>*>(), solution,
                                     estimatedErrorPerCell);

  GridRefinement::refine_and_coarsen_fixed_number(triangulation, estimatedErrorPerCell, 0.3, 0.03);

  triangulation.execute_coarsening_and_refinement();
};

template <int dim>
void HelmholtzProblem<dim>::process_solution(const unsigned int cycle) {
  Vector<float> difference_per_cell(triangulation.n_active_cells());

  VectorTools::integrate_difference(dof_handler, solution, Solution<dim>(), difference_per_cell,
                                    QGauss<dim>(fe.degree + 1), VectorTools::L2_norm);
  const double L2_error =
      VectorTools::compute_global_error(triangulation, difference_per_cell, VectorTools::L2_norm);

  VectorTools::integrate_difference(dof_handler, solution, Solution<dim>(), difference_per_cell,
                                    QGauss<dim>(fe.degree + 1), VectorTools::H1_seminorm);
  const double H1_error = VectorTools::compute_global_error(triangulation, difference_per_cell,
                                                            VectorTools::H1_seminorm);

  const QTrapezoid<1> q_trapez;
  const QIterated<dim> q_iterated(q_trapez, fe.degree * 2 + 1);
  VectorTools::integrate_difference(dof_handler, solution, Solution<dim>(), difference_per_cell,
                                    q_iterated, VectorTools::Linfty_norm);
  const double Linfty_error = VectorTools::compute_global_error(triangulation, difference_per_cell,
                                                                VectorTools::Linfty_norm);

  const unsigned int n_active_cells = triangulation.n_active_cells();
  const unsigned int n_dofs = dof_handler.n_dofs();

  convergence_table.add_value("cycle", cycle);
  convergence_table.add_value("cells", n_active_cells);
  convergence_table.add_value("dofs", n_dofs);
  convergence_table.add_value("L2", L2_error);
  convergence_table.add_value("H1", H1_error);
  convergence_table.add_value("Linfty", Linfty_error);

  convergence_table.set_precision("L2", 3);
  convergence_table.set_precision("H1", 3);
  convergence_table.set_precision("Linfty", 3);

  convergence_table.set_scientific("L2", true);
  convergence_table.set_scientific("H1", true);
  convergence_table.set_scientific("Linfty", true);

  convergence_table.set_tex_caption("cells", "\\# cells");
  convergence_table.set_tex_caption("dofs", "\\# dofs");
  convergence_table.set_tex_caption("L2", "@f$L^2@f$-error");
  convergence_table.set_tex_caption("H1", "@f$H^1@f$-error");
  convergence_table.set_tex_caption("Linfty", "@f$L^\\infty@f$-error");

  convergence_table.set_tex_format("cells", "r");
  convergence_table.set_tex_format("dofs", "r");
};

template <int dim>
void HelmholtzProblem<dim>::output() const {
  if (output_solution_stream != nullptr) {
    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "result");
    data_out.build_patches(fe.degree);

    auto format = util::output_format_to_dealii(output_solution_format);
    data_out.write(*output_solution_stream, format);
  }

  std::ostringstream table_stream;
  convergence_table.write_text(table_stream);
  write_log("CONVERGENCE", {{"table", table_stream.str()}}, log_stream);
};

} // namespace problem
} // namespace plksim