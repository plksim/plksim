#include "problem_impl.hh"

#include "problem/advection_problem.hh"
#include "problem/elastic_problem.hh"
#include "problem/helmholtz_problem.hh"
#include "problem/laplace_problem.hh"

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
    throw std::runtime_error("INVALID_FORMAT: type is not defined");
  } else {
    type = obj["type"].get<std::string>();
  }
};

void Problem::Impl::set_log_stream(std::ostream& stream) {
  log_stream = &stream;
};

void Problem::Impl::set_output_grid_stream(std::ostream& stream, const io::OutputFormat format) {
  output_grid_stream = &stream;
  output_grid_format = format;
};

void Problem::Impl::set_output_solution_stream(std::ostream& stream,
                                               const io::OutputFormat format) {
  output_solution_stream = &stream;
  output_solution_format = format;
};

void Problem::Impl::compute() {
  std::unique_ptr<problem::ProblemBase> prob;
  if (type == "laplace") {
    prob = std::make_unique<problem::LaplaceProblem<2>>();
    // prob.set_refine_criterion("kelly");
  } else if (type == "helmholtz") {
    prob = std::make_unique<problem::HelmholtzProblem<2>>();
  } else if (type == "elastic") {
    prob = std::make_unique<problem::ElasticProblem<2>>();
  } else if (type == "advection") {
    prob = std::make_unique<problem::AdvectionProblem<2>>();
  } else {
    throw std::runtime_error("ILLEGAL_TYPE: type is not legal");
  }

  if (prob != nullptr) {
    if (log_stream != nullptr) {
      prob->set_log_stream(*log_stream);
    }

    if (output_solution_stream != nullptr) {
      prob->set_output_solution_stream(*output_solution_stream, output_solution_format);
    }

    if (output_grid_stream != nullptr) {
      prob->set_output_grid_stream(*output_grid_stream, output_grid_format);
    }

    prob->compute();
  }
};

std::string Problem::Impl::get_type() const {
  return type;
};

} // namespace plksim