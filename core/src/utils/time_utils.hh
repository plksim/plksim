#pragma once

#include <ctime>
#include <iomanip>
#include <sstream>

namespace plksim {
namespace utils {

std::string now_utc_time_iso() {
  auto now = std::chrono::system_clock::now();
  auto itt = std::chrono::system_clock::to_time_t(now);

  auto sec = std::chrono::duration_cast<std::chrono::seconds>(now.time_since_epoch());
  auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch());
  auto ms = msec - sec;

  std::ostringstream ss;
  ss << std::put_time(gmtime(&itt), "%FT%T") << "." << std::setfill('0') << std::setw(3) << ms.count() << "Z";

  return ss.str();
};

} // namespace utils
} // namespace plksim