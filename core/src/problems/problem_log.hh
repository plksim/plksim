#pragma once

#include "utils/time_utils.hh"

#include <nlohmann/json.hpp>

#include <iostream>

namespace plksim {
namespace problems {

void write_log(const std::string tag, const nlohmann::json data, std::ostream& log_stream) {
  nlohmann::json j;
  j["time"] = utils::now_utc_time_iso();
  j["tag"] = tag;
  j["data"] = data;

  log_stream << j.dump() << std::endl;
};

} // namespace problems
} // namespace plksim