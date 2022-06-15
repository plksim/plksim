#include "../include/plksim/version.hh"

#include <gmock/gmock.h>

using ::testing::MatchesRegex;

namespace plksim {

TEST(version, versionCode) {
  ASSERT_GE(VERSION_MAJOR, 0);
  ASSERT_GE(VERSION_MINOR, 0);
  ASSERT_GE(VERSION_PATCH, 0);
}

TEST(version, versionString) {
  std::cout << "[**********] VERSION = " << VERSION << std::endl;
  ASSERT_THAT(VERSION, MatchesRegex("^[0-9]*\\.[0-9]*\\.[0-9]*$"));
}

} // namespace plksim