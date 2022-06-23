#pragma once

#include "problem_base.hh"
#include "utils/output_format_utils.hh"

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

namespace plksim {
namespace problems {

using namespace dealii;

template <int dim>
class LaplaceProblem : public ProblemBase {
public:
  LaplaceProblem();
  ~LaplaceProblem();
  void compute(std::ostream& out, const io::OutputFormat outFormat);

private:
  void makeGrid();
  void setupSystem();
  void assembleSystem();
  void solve();
  void output(std::ostream& out, const io::OutputFormat outFormat);

  Triangulation<dim> triangulation;
  FE_Q<dim> fe;
  DoFHandler<dim> dofHandler;

  SparsityPattern sparsityPattern;
  SparseMatrix<double> systemMatrix;

  Vector<double> solution;
  Vector<double> systemRhs;
};

template <int dim>
LaplaceProblem<dim>::LaplaceProblem() : fe(1), dofHandler(triangulation){};

template <int dim>
LaplaceProblem<dim>::~LaplaceProblem(){};

template <int dim>
void LaplaceProblem<dim>::compute(std::ostream& out, const io::OutputFormat outFormat) {
  makeGrid();
  setupSystem();
  assembleSystem();
  solve();
  output(out, outFormat);
};

template <int dim>
void LaplaceProblem<dim>::makeGrid() {
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(5);
};

template <int dim>
void LaplaceProblem<dim>::setupSystem() {
  dofHandler.distribute_dofs(fe);

  DynamicSparsityPattern dsp(dofHandler.n_dofs());
  DoFTools::make_sparsity_pattern(dofHandler, dsp);
  sparsityPattern.copy_from(dsp);

  systemMatrix.reinit(sparsityPattern);

  solution.reinit(dofHandler.n_dofs());
  systemRhs.reinit(dofHandler.n_dofs());
};

template <int dim>
void LaplaceProblem<dim>::assembleSystem() {
  QGauss<dim> quadratureFormula(fe.degree + 1);
  FEValues<dim> feValues(fe, quadratureFormula, update_values | update_gradients | update_JxW_values);

  const unsigned int dofsPerCell = fe.n_dofs_per_cell();

  FullMatrix<double> cellMatrix(dofsPerCell, dofsPerCell);
  Vector<double> cellRhs(dofsPerCell);

  std::vector<types::global_dof_index> localDofIndices(dofsPerCell);

  for (const auto& cell : dofHandler.active_cell_iterators()) {
    feValues.reinit(cell);

    cellMatrix = 0;
    cellRhs = 0;

    for (const unsigned int idx : feValues.quadrature_point_indices()) {
      for (const unsigned int i : feValues.dof_indices())
        for (const unsigned int j : feValues.dof_indices())
          cellMatrix(i, j) += (feValues.shape_grad(i, idx) * // grad phi_i(x_q)
                               feValues.shape_grad(j, idx) * // grad phi_j(x_q)
                               feValues.JxW(idx));           // dx

      for (const unsigned int i : feValues.dof_indices())
        cellRhs(i) += (feValues.shape_value(i, idx) * // phi_i(x_q)
                       1. *                           // f(x_q)
                       feValues.JxW(idx));            // dx
    }
    cell->get_dof_indices(localDofIndices);

    for (const unsigned int i : feValues.dof_indices())
      for (const unsigned int j : feValues.dof_indices())
        systemMatrix.add(localDofIndices[i], localDofIndices[j], cellMatrix(i, j));

    for (const unsigned int i : feValues.dof_indices())
      systemRhs(localDofIndices[i]) += cellRhs(i);
  }

  std::map<types::global_dof_index, double> boundaryValues;
  VectorTools::interpolate_boundary_values(dofHandler, 0, Functions::ZeroFunction<dim>(), boundaryValues);
  MatrixTools::apply_boundary_values(boundaryValues, systemMatrix, solution, systemRhs);
};

template <int dim>
void LaplaceProblem<dim>::solve() {
  SolverControl sc(1000, 1e-12);
  SolverCG<Vector<double>> solver(sc);

  solver.solve(systemMatrix, solution, systemRhs, PreconditionIdentity());
};

template <int dim>
void LaplaceProblem<dim>::output(std::ostream& out, const io::OutputFormat outFormat) {
  DataOut<dim> dataOut;

  dataOut.attach_dof_handler(dofHandler);
  dataOut.add_data_vector(solution, "result");
  dataOut.build_patches();

  auto format = utils::outputFormatToDealii(outFormat);
  dataOut.write(out, format);
};

} // namespace problems
} // namespace plksim