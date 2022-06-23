#include "problem_impl.hh"

#include "problems/laplace_problem.hh"

#include <nlohmann/json.hpp>

#include <exception>
#include <sstream>

using json = nlohmann::json;

namespace plksim {

Problem::Impl::Impl(){};
Problem::Impl::~Impl(){};

void Problem::Impl::load(const std::string& jsonStr) {
  auto obj = json::parse(jsonStr);

  if (!obj["type"].is_string()) {
    throw std::runtime_error("wrong format: type is not defined");
  } else {
    type = obj["type"].get<std::string>();
  }
};

void Problem::Impl::compute(std::ostream& out, const io::OutputFormat outFormat) {
  if (type == "laplace") {
    problems::LaplaceProblem<2> prob;
    return prob.compute(out, outFormat);
  }
};

std::string Problem::Impl::getType() const {
  return type;
};

} // namespace plksim