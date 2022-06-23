#include "problem.hh"

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

  ASSERT_EQ(prob.getType(), "test");
}

TEST(problem, compute) {
  plksim::Problem prob;
  prob.load(R"(
    {
      "type":"laplace"
    }
  )");

  std::ofstream out("test.vtk");
  prob.compute(out);
}

} // namespace plksim_test