#pragma once

#include "set_up_base.hh"

namespace plksim {
namespace data {

template <class T, int dim>
struct SetUp : public SetUpBase<dim> {
  virtual const Function<dim>& get_boundary_values() const override;
  virtual const Function<dim>& get_right_hand_side() const override;

  virtual void create_coarse_grid(Triangulation<dim>& coarse_grid) const override;

private:
  static const typename T::BoundaryValues boundary_values;
  static const typename T::RightHandSide right_hand_side;
};

template <class T, int dim>
const typename T::BoundaryValues SetUp<T, dim>::boundary_values;

template <class T, int dim>
const typename T::RightHandSide SetUp<T, dim>::right_hand_side;

template <class T, int dim>
const Function<dim>& SetUp<T, dim>::get_boundary_values() const {
  return boundary_values;
}

template <class T, int dim>
const Function<dim>& SetUp<T, dim>::get_right_hand_side() const {
  return right_hand_side;
}

template <class T, int dim>
void SetUp<T, dim>::create_coarse_grid(Triangulation<dim>& coarse_grid) const {
  T::create_coarse_grid(coarse_grid);
}

} // namespace data
} // namespace plksim