#include "laplace_problem.hh"

namespace plksim {
namespace problems {

void LaplaceProblem::compute(std::ostream& out, const io::OutputFormat outFormat) {
  makeGrid();
  setupSystem();
  assembleSystem();
  solve();
  output(out, outFormat);
};

void LaplaceProblem::makeGrid() {
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(5);
};

void LaplaceProblem::setupSystem() {
  dofHandler.distribute_dofs(fe);

  DynamicSparsityPattern dsp(dofHandler.n_dofs());
  DoFTools::make_sparsity_pattern(dofHandler, dsp);
  sparsityPattern.copy_from(dsp);

  systemMatrix.reinit(sparsityPattern);

  solution.reinit(dofHandler.n_dofs());
  systemRhs.reinit(dofHandler.n_dofs());
};

void LaplaceProblem::assembleSystem() {
  QGauss<2> quadratureFormula(fe.degree + 1);
  FEValues<2> feValues(fe, quadratureFormula, update_values | update_gradients | update_JxW_values);

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
  VectorTools::interpolate_boundary_values(dofHandler, 0, Functions::ZeroFunction<2>(), boundaryValues);
  MatrixTools::apply_boundary_values(boundaryValues, systemMatrix, solution, systemRhs);
};

void LaplaceProblem::solve() {
  SolverControl sc(1000, 1e-12);
  SolverCG<Vector<double>> solver(sc);

  solver.solve(systemMatrix, solution, systemRhs, PreconditionIdentity());
};

void LaplaceProblem::output(std::ostream& out, const io::OutputFormat outFormat) {
  DataOut<2> dataOut;
  dataOut.attach_dof_handler(dofHandler);
  dataOut.add_data_vector(solution, "result");
  dataOut.build_patches();

  switch (outFormat) {
  case io::OutputFormat::dx:
    dataOut.write(out, DataOutBase::dx);
    break;
  case io::OutputFormat::eps:
    dataOut.write(out, DataOutBase::eps);
    break;
  case io::OutputFormat::gmv:
    dataOut.write(out, DataOutBase::gmv);
    break;
  case io::OutputFormat::gnuplot:
    dataOut.write(out, DataOutBase::gnuplot);
    break;
  case io::OutputFormat::hdf5:
    dataOut.write(out, DataOutBase::hdf5);
    break;
  case io::OutputFormat::povray:
    dataOut.write(out, DataOutBase::povray);
    break;
  case io::OutputFormat::svg:
    dataOut.write(out, DataOutBase::svg);
    break;
  case io::OutputFormat::tecplot:
    dataOut.write(out, DataOutBase::tecplot);
    break;
  case io::OutputFormat::ucd:
    dataOut.write(out, DataOutBase::ucd);
    break;
  case io::OutputFormat::vtk:
    dataOut.write(out, DataOutBase::vtk);
    break;
  case io::OutputFormat::vtu:
    dataOut.write(out, DataOutBase::vtu);
    break;
  }
};

} // namespace problems
} // namespace plksim