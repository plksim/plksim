#pragma once

#include "problem.hh"

#include <string>

namespace plksim {

class Problem::Impl final {
public:
  Impl();
  ~Impl();

  void load(const std::string& jsonStr);
  void compute(std::ostream& out, const io::OutputFormat outFormat);

  std::string getType() const;

private:
  std::string type;
};

} // namespace plksim