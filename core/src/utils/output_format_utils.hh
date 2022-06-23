#pragma once

#include "io/output_format.hh"

#include <deal.II/numerics/data_out.h>

namespace plksim {
namespace utils {

dealii::DataOutBase::OutputFormat outputFormatToDealii(const io::OutputFormat format) {
  switch (format) {
  case io::OutputFormat::dx:
    return dealii::DataOutBase::dx;
  case io::OutputFormat::eps:
    return dealii::DataOutBase::eps;
  case io::OutputFormat::gmv:
    return dealii::DataOutBase::gmv;
  case io::OutputFormat::gnuplot:
    return dealii::DataOutBase::gnuplot;
  case io::OutputFormat::hdf5:
    return dealii::DataOutBase::hdf5;
  case io::OutputFormat::povray:
    return dealii::DataOutBase::povray;
  case io::OutputFormat::svg:
    return dealii::DataOutBase::svg;
  case io::OutputFormat::tecplot:
    return dealii::DataOutBase::tecplot;
  case io::OutputFormat::ucd:
    return dealii::DataOutBase::ucd;
  case io::OutputFormat::vtk:
    return dealii::DataOutBase::vtk;
  case io::OutputFormat::vtu:
    return dealii::DataOutBase::vtu;
  default:
    return dealii::DataOutBase::none;
  }
};

} // namespace utils
} // namespace plksim