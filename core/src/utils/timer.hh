#pragma once

#include <functional>
#include <map>
#include <shared_mutex>

namespace plksim {
namespace utils {

uint64_t set_timeout(std::function<void()> func, uint32_t delay_ms);
void clear_timeout(uint64_t timer_id);

static std::shared_mutex timer_mutex;
static std::map<uint64_t, bool> timer_map;

} // namespace utils
} // namespace plksim
