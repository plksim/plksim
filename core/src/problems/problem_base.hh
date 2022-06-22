#pragma once

#include "io/output_format.hh"

#include <iostream>

namespace plksim {
namespace problems {

class ProblemBase {
public:
  virtual void compute(std::ostream& out, const io::OutputFormat outFormat) = 0;
};
} // namespace problems
} // namespace plksim