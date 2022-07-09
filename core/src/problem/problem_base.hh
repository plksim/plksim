#pragma once

#include "io/output_format.hh"

#include <iostream>

namespace plksim {
namespace problem {

class ProblemBase {
public:
  virtual ~ProblemBase() = default;

  virtual void compute() = 0;

  void set_refine_criterion(const std::string& v);

  void set_log_stream(std::ostream& log_stream);
  void set_output_grid_stream(std::ostream& out_stream, const io::OutputFormat out_format);
  void set_output_solution_stream(std::ostream& out_stream, const io::OutputFormat out_format);

protected:
  std::string refine_criterion;

  std::ostream* log_stream;

  std::ostream* output_grid_stream;
  io::OutputFormat output_grid_format = io::OutputFormat::vtu;

  std::ostream* output_solution_stream;
  io::OutputFormat output_solution_format = io::OutputFormat::vtu;
};

void ProblemBase::set_refine_criterion(const std::string& v) {
  refine_criterion = v;
}

void ProblemBase::set_log_stream(std::ostream& stream) {
  log_stream = &stream;
};

void ProblemBase::set_output_grid_stream(std::ostream& stream, const io::OutputFormat format) {
  output_grid_stream = &stream;
  output_grid_format = format;
};

void ProblemBase::set_output_solution_stream(std::ostream& stream, const io::OutputFormat format) {
  output_solution_stream = &stream;
  output_solution_format = format;
};

} // namespace problem
} // namespace plksim