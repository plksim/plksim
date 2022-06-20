#pragma once

#include "problem.hh"

#include <string>

namespace plksim {

class Problem::Impl final {
public:
  Impl();
  ~Impl();

  void load(const std::string& jsonStr);

  std::string getType();

private:
  std::string mType;
};

} // namespace plksim