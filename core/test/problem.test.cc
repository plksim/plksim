#include "problem.hh"
#include "io/output_format.hh"

#include <fstream>

#include <gtest/gtest.h>

namespace plksim_test {

TEST(problem, load) {
  plksim::Problem prob;
  prob.load(R"(
    {
      "type":"test"
    }
  )");

  ASSERT_EQ(prob.get_type(), "test");
}

TEST(problem, compute_laplace) {
  plksim::Problem prob;
  prob.load(R"(
    {
      "type":"laplace"
    }
  )");

  prob.set_log_stream(std::cout);

  std::ofstream out("laplace.vtu");
  prob.set_output_solution_stream(out);

  prob.compute();
}

TEST(problem, compute_helmholtz) {
  plksim::Problem prob;
  prob.load(R"(
    {
      "type":"helmholtz"
    }
  )");

  prob.set_log_stream(std::cout);

  std::ofstream out("helmholtz.vtu");
  prob.set_output_solution_stream(out);

  prob.compute();
}

TEST(problem, compute_elastic) {
  plksim::Problem prob;
  prob.load(R"(
    {
      "type":"elastic"
    }
  )");

  prob.set_log_stream(std::cout);

  std::ofstream out("elastic.vtu");
  prob.set_output_solution_stream(out);

  prob.compute();
}

TEST(problem, compute_advection) {
  plksim::Problem prob;
  prob.load(R"(
    {
      "type":"advection"
    }
  )");

  prob.set_log_stream(std::cout);

  std::ofstream out("advection.vtu");
  prob.set_output_solution_stream(out);

  prob.compute();
}

} // namespace plksim_test