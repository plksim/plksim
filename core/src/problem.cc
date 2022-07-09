#include "problem_impl.hh"

#include <memory>

namespace plksim {

Problem::Problem() : impl(std::make_unique<Impl>()){};

Problem::~Problem() = default;

void Problem::load(const std::string& json_str) {
  return impl->load(json_str);
};

void Problem::set_log_stream(std::ostream& stream) {
  return impl->set_log_stream(stream);
};

void Problem::set_output_grid_stream(std::ostream& stream, const io::OutputFormat format) {
  return impl->set_output_grid_stream(stream, format);
};

void Problem::set_output_solution_stream(std::ostream& stream, const io::OutputFormat format) {
  return impl->set_output_solution_stream(stream, format);
};

void Problem::compute() {
  return impl->compute();
};

std::string Problem::get_type() const {
  return impl->get_type();
};

} // namespace plksim