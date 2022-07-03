#include "problem_impl.hh"

#include "problems/helmholtz_problem.hh"
#include "problems/laplace_problem.hh"

#include <nlohmann/json.hpp>

#include <exception>
#include <sstream>

using json = nlohmann::json;

namespace plksim {

Problem::Impl::Impl(){};
Problem::Impl::~Impl(){};

void Problem::Impl::load(const std::string& json_str) {
  auto obj = json::parse(json_str);

  if (!obj["type"].is_string()) {
    throw std::runtime_error("wrong format: type is not defined");
  } else {
    type = obj["type"].get<std::string>();
  }
};

void Problem::Impl::compute(std::ostream& out_stream, const io::OutputFormat out_format, std::ostream& log_stream) {
  if (type == "laplace") {
    problems::LaplaceProblem<2> prob;
    prob.compute(out_stream, out_format, log_stream);
  } else if (type == "helmholtz") {
    problems::HelmholtzProblem<2> prob;
    prob.compute(out_stream, out_format, log_stream);
  }
};

std::string Problem::Impl::get_type() const {
  return type;
};

} // namespace plksim