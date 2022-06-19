#include "samples.hh"

#include <gtest/gtest.h>

namespace plksim_test {

TEST(samples, sampleMeshSvg) {
  ASSERT_NO_THROW({ auto svg = plksim::sampleMeshSvg(); });
}

} // namespace plksim_test