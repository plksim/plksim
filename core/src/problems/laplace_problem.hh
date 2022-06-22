#pragma once

#include "problem_base.hh"

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

class LaplaceProblem : public ProblemBase {
public:
  LaplaceProblem() : fe(1), dofHandler(triangulation){};
  ~LaplaceProblem(){};

  void compute(std::ostream& out, const io::OutputFormat outFormat);

private:
  void makeGrid();
  void setupSystem();
  void assembleSystem();
  void solve();
  void output(std::ostream& out, const io::OutputFormat outFormat);

  Triangulation<2> triangulation;
  FE_Q<2> fe;
  DoFHandler<2> dofHandler;

  SparsityPattern sparsityPattern;
  SparseMatrix<double> systemMatrix;

  Vector<double> solution;
  Vector<double> systemRhs;
};

} // namespace problems
} // namespace plksim