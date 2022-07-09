#pragma once

#include "problem.hh"

#include <string>

namespace plksim {

class Problem::Impl final {
public:
  Impl();
  ~Impl();

  void load(const std::string& json_str);

  void set_log_stream(std::ostream& stream);
  void set_output_grid_stream(std::ostream& stream, const io::OutputFormat format);
  void set_output_solution_stream(std::ostream& stream, const io::OutputFormat format);

  void compute();

  std::string get_type() const;

private:
  std::ostream* log_stream;

  std::ostream* output_grid_stream;
  io::OutputFormat output_grid_format = io::OutputFormat::vtu;

  std::ostream* output_solution_stream;
  io::OutputFormat output_solution_format = io::OutputFormat::vtu;

  std::string type;
};

} // namespace plksim