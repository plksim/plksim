#include "../include/plksim/version.hh"

#include <gmock/gmock.h>

using ::testing::MatchesRegex;

namespace plksim {

TEST(version, versionCode) {
  ASSERT_GE(MAJOR_VERSION, 0);
  ASSERT_GE(MINOR_VERSION, 0);
  ASSERT_GE(PATCH_VERSION, 0);
}

TEST(version, versionString) {
  std::cout << "[**********] VERSION = " << VERSION << std::endl;
  ASSERT_THAT(VERSION, MatchesRegex("^[0-9]*\\.[0-9]*\\.[0-9]*$"));
}

} // namespace plksim