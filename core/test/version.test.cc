#include "version.hh"

#include <gtest/gtest.h>

namespace plksim_test {

TEST(version, versionCode) {
  ASSERT_GE(plksim::versionMajor, 0);
  ASSERT_GE(plksim::versionMinor, 0);
  ASSERT_GE(plksim::versionPatch, 0);
}

TEST(version, versionString) {
  std::cout << "[**********] version = " << plksim::version << std::endl;
  std::ostringstream buffer;
  buffer << plksim::versionMajor << "." << plksim::versionMinor << "." << plksim::versionPatch;
  ASSERT_EQ(plksim::version, buffer.str());
}

} // namespace plksim_test