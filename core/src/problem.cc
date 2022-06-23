#include "problem_impl.hh"

#include <memory>

namespace plksim {

Problem::Problem() : impl(std::make_unique<Impl>()){};

Problem::~Problem() = default;

void Problem::load(const std::string& jsonStr) {
  return impl->load(jsonStr);
};

void Problem::compute(std::ostream& out, const io::OutputFormat outFormat) {
  return impl->compute(out, outFormat);
};

std::string Problem::getType() const {
  return impl->getType();
};

} // namespace plksim