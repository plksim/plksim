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

  std::ofstream out("laplace.vtk");
  prob.compute(out);
}

TEST(problem, compute_helmholtz) {
  plksim::Problem prob;
  prob.load(R"(
    {
      "type":"helmholtz"
    }
  )");

  std::ofstream out("helmholtz.vtk");
  std::ofstream log("helmholtz.log");
  prob.compute(out, plksim::io::OutputFormat::vtk, log);
}

} // namespace plksim_test