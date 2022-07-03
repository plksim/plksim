#include "timer.hh"

#include <atomic>
#include <chrono>
#include <thread>

namespace plksim {
namespace utils {

static std::atomic_ulong timer_id_seq = 0;

uint64_t set_timeout(std::function<void()> func, uint32_t delay_ms) {
  uint64_t timer_id = ++timer_id_seq;

  timer_mutex.lock();
  timer_map.insert({timer_id, false});
  timer_mutex.unlock();

  std::thread t([timer_id, func, delay_ms]() {
    std::this_thread::sleep_for(std::chrono::milliseconds(delay_ms));

    bool cancelled = false;

    timer_mutex.lock();
    auto itr = timer_map.find(timer_id);
    if (itr != timer_map.end()) {
      cancelled = itr->second;

      timer_map.erase(itr);
    }
    timer_mutex.unlock();

    if (cancelled) {
      return;
    }

    func();
  });
  t.detach();

  return timer_id;
};

void clear_timeout(uint64_t timer_id) {
  timer_mutex.lock();

  auto itr = timer_map.find(timer_id);
  if (itr != timer_map.end()) {
    itr->second = true;
  }

  timer_mutex.unlock();
};

} // namespace utils
} // namespace plksim