#include "version.hh"

#include <gtest/gtest.h>

namespace plksim_test {

TEST(version, version_nums) {
  ASSERT_GE(plksim::version_major, 0);
  ASSERT_GE(plksim::version_minor, 0);
  ASSERT_GE(plksim::version_patch, 0);
}

TEST(version, version_string) {
  std::cout << "[**********] version = " << plksim::version << std::endl;
  std::ostringstream buffer;
  buffer << plksim::version_major << "." << plksim::version_minor << "." << plksim::version_patch;
  ASSERT_EQ(plksim::version, buffer.str());
}

} // namespace plksim_test