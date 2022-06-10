#include "timer.hh"

#include <atomic>
#include <chrono>
#include <thread>

namespace plksim {
namespace utils {

static std::atomic_ulong timerIdSeq = 0;

uint64_t setTimeout(std::function<void()> func, uint32_t delayMs) {
  uint64_t timerId = ++timerIdSeq;

  timerMutex.lock();
  timerMap.insert({timerId, false});
  timerMutex.unlock();

  std::thread t([timerId, func, delayMs]() {
    std::this_thread::sleep_for(std::chrono::milliseconds(delayMs));

    bool cancelled = false;

    timerMutex.lock();
    auto itr = timerMap.find(timerId);
    if (itr != timerMap.end()) {
      cancelled = itr->second;

      timerMap.erase(itr);
    }
    timerMutex.unlock();

    if (cancelled) {
      return;
    }

    func();
  });
  t.detach();

  return timerId;
};

void clearTimeout(uint64_t timerId) {
  timerMutex.lock();

  auto itr = timerMap.find(timerId);
  if (itr != timerMap.end()) {
    itr->second = true;
  }

  timerMutex.unlock();
};

} // namespace utils
} // namespace plksim