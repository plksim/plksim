#pragma once

#include <functional>
#include <map>
#include <shared_mutex>

namespace plksim {
namespace utils {

uint64_t setTimeout(std::function<void()> func, uint32_t delayMs);
void clearTimeout(uint64_t timerId);

static std::shared_mutex timerMutex;
static std::map<uint64_t, bool> timerMap;

} // namespace utils
} // namespace plksim
