#include "problem_impl.hh"

#include <memory>

namespace plksim {

Problem::Problem() : mImpl(std::make_unique<Impl>()) {
}

Problem::~Problem() = default;

void Problem::load(const std::string& jsonStr) {
  mImpl->load(jsonStr);
};

std::string Problem::getType() {
  return mImpl->getType();
};

} // namespace plksim