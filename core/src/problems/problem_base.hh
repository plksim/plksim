#pragma once

#include "io/output_format.hh"

#include <iostream>

namespace plksim {
namespace problems {

class ProblemBase {
public:
  virtual void compute(std::ostream& out_stream, const io::OutputFormat out_format, std::ostream& log_stream) = 0;
};

} // namespace problems
} // namespace plksim