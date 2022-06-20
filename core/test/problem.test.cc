#include "problem.hh"

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

} // namespace plksim_test