#pragma once

#include "problem.hh"

#include <string>

namespace plksim {

class Problem::Impl final {
public:
  Impl();
  ~Impl();

  void load(const std::string& json_str);
  void compute(std::ostream& out_stream, const io::OutputFormat out_format, std::ostream& log_stream);

  std::string get_type() const;

private:
  std::string type;
};

} // namespace plksim