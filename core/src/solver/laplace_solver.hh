#pragma once

#include "solver.hh"

namespace plksim {
namespace solver {

template <int dim>
class LaplaceSolver : public Solver<dim> {
public:
  LaplaceSolver(Triangulation<dim>& triangulation, const FiniteElement<dim>& fe,
                const Quadrature<dim>& quadrature, const Quadrature<dim - 1>& face_quadrature,
                const Function<dim>& rhs_function, const Function<dim>& boundary_values);

  virtual void output_solution(std::ostream& out_stream,
                               const io::OutputFormat out_format) const override;

protected:
  const SmartPointer<const Function<dim>> rhs_function;
  virtual void assemble_rhs(Vector<double>& rhs) const override;
};

template <int dim>
LaplaceSolver<dim>::LaplaceSolver(Triangulation<dim>& triangulation, const FiniteElement<dim>& fe,
                                  const Quadrature<dim>& quadrature,
                                  const Quadrature<dim - 1>& face_quadrature,
                                  const Function<dim>& rhs_function,
                                  const Function<dim>& boundary_values)
    : SolverBase<dim>(triangulation), Solver<dim>(triangulation, fe, quadrature, face_quadrature,
                                                  boundary_values),
      rhs_function(&rhs_function) {
}

template <int dim>
void LaplaceSolver<dim>::output_solution(std::ostream& out_stream,
                                         const io::OutputFormat out_format) const {
  DataOut<dim> data_out;
  data_out.attach_dof_handler(this->dof_handler);
  data_out.add_data_vector(this->solution, "result");
  data_out.build_patches();

  auto format = util::output_format_to_dealii(out_format);
  data_out.write(out_stream, format);
}

template <int dim>
void LaplaceSolver<dim>::assemble_rhs(Vector<double>& rhs) const {
  FEValues<dim> fe_values(*this->fe, *this->quadrature,
                          update_values | update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = this->fe->n_dofs_per_cell();
  const unsigned int n_q_points = this->quadrature->size();

  Vector<double> cell_rhs(dofs_per_cell);
  std::vector<double> rhs_values(n_q_points);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto& cell : this->dof_handler.active_cell_iterators()) {
    cell_rhs = 0;

    fe_values.reinit(cell);

    rhs_function->value_list(fe_values.get_quadrature_points(), rhs_values);

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
      for (unsigned int i = 0; i < dofs_per_cell; ++i) {
        cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q)
                        rhs_values[q_point] *               // f((x_q)
                        fe_values.JxW(q_point));            // dx
      }
    }

    cell->get_dof_indices(local_dof_indices);
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      rhs(local_dof_indices[i]) += cell_rhs(i);
    }
  }
}

} // namespace solver
} // namespace plksim