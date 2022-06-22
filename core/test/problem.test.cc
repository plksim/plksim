#include "problem.hh"

#include <fstream>

#include <gtest/gtest.h>

namespace plksim_test {

TEST(problem, load) {
  plksim::Problem pb;
  pb.load(R"(
    {
      "type":"test"
    }
  )");

  ASSERT_EQ(pb.getType(), "test");
}

TEST(problem, compute) {
  plksim::Problem pb;
  pb.load(R"(
    {
      "type":"laplace"
    }
  )");

  std::ofstream out("test.vtk");
  pb.compute(out);
}

} // namespace plksim_test