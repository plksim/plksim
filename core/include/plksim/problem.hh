#pragma once

#include <memory>
#include <string>

namespace plksim {

class Problem {
public:
  Problem();
  ~Problem();

  /**
   * Load problem definition from JSON.
   *
   * @param jsonStr    JSON string.
   */
  void load(const std::string& jsonStr);

  /**
   * Get the type.
   */
  std::string getType();

private:
  class Impl;
  const std::unique_ptr<Impl> mImpl;
};

} // namespace plksim