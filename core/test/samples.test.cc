#include "../include/plksim/samples.hh"

#include <gtest/gtest.h>

namespace plksim {

TEST(samples, sampleMeshSvg) {
  ASSERT_NO_THROW({ auto svg = sampleMeshSvg(); });
}

} // namespace plksim