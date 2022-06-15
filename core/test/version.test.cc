#include "../include/plksim/version.hh"

#include <gmock/gmock.h>

using ::testing::MatchesRegex;

namespace plksim {

TEST(version, versionCode) {
  ASSERT_GE(versionMajor, 0);
  ASSERT_GE(versionMinor, 0);
  ASSERT_GE(versionPatch, 0);
}

TEST(version, versionString) {
  std::cout << "[**********] version = " << version << std::endl;
  ASSERT_THAT(version, MatchesRegex("^[0-9]*\\.[0-9]*\\.[0-9]*$"));
}

} // namespace plksim