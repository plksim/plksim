#pragma once

#include "util/time_util.hh"

#include <nlohmann/json.hpp>

#include <iostream>
#include <memory>

namespace plksim {
namespace problem {

void write_log(const std::string tag, const nlohmann::json data, std::ostream* log_stream) {
  if (log_stream != nullptr) {
    nlohmann::json j;

    j["time"] = util::now_utc_time_iso();
    j["tag"] = tag;
    j["data"] = data;

    *log_stream << j.dump() << std::endl;
  }
};

} // namespace problem
} // namespace plksim