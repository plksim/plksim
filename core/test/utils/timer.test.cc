#include "utils/timer.hh"

#include <gtest/gtest.h>

#include <atomic>
#include <chrono>
#include <iostream>
#include <thread>

namespace plksim_test {

TEST(timer, set_timeout) {
  std::atomic<bool> tick(false);

  auto t0 = std::chrono::steady_clock::now();

  plksim::utils::set_timeout([&]() { tick.store(true); }, 500);
  std::this_thread::sleep_for(std::chrono::milliseconds(600));

  auto t1 = std::chrono::steady_clock::now();
  auto ms1 = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  ASSERT_GE(ms1, 600);

  ASSERT_EQ(tick.load(), true);
}

TEST(timer, clear_timeout) {
  std::atomic<bool> tick(false);

  auto t0 = std::chrono::steady_clock::now();

  auto timer = plksim::utils::set_timeout([&]() { tick.store(true); }, 500);

  std::this_thread::sleep_for(std::chrono::milliseconds(100));
  auto t1 = std::chrono::steady_clock::now();
  auto ms1 = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  plksim::utils::clear_timeout(timer);
  ASSERT_GE(ms1, 100);
  ASSERT_LT(ms1, 120);

  std::this_thread::sleep_for(std::chrono::milliseconds(500));
  auto t2 = std::chrono::steady_clock::now();
  auto ms2 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t0).count();
  ASSERT_GE(ms2, 600);

  ASSERT_EQ(tick.load(), false);
}

} // namespace plksim_test