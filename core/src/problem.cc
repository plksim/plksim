#include "problem_impl.hh"

#include <memory>

namespace plksim {

Problem::Problem() : impl(std::make_unique<Impl>()){};

Problem::~Problem() = default;

void Problem::load(const std::string& json_str) {
  return impl->load(json_str);
};

void Problem::compute(std::ostream& out_stream, const io::OutputFormat out_format, std::ostream& log_stream) {
  return impl->compute(out_stream, out_format, log_stream);
};

std::string Problem::get_type() const {
  return impl->get_type();
};

} // namespace plksim