#include "samples.hh"

#include <gtest/gtest.h>

namespace plksim_test {

TEST(samples, sample_mesh_svg) {
  ASSERT_NO_THROW({ auto svg = plksim::sample_mesh_svg(); });
}

} // namespace plksim_test