#include "../include/plksim/version.hh"

#include <gtest/gtest.h>

namespace plksim {

TEST(version, versionCode) {
  ASSERT_GE(MAJOR_VERSION, 0);
  ASSERT_GE(MINOR_VERSION, 0);
  ASSERT_GE(PATCH_VERSION, 0);
}

TEST(version, versionString) {
  std::cout << "[**********] VERSION_STRING = " << VERSION_STRING << std::endl;
}

} // namespace plksim