#pragma once

#include "io/output_format.hh"

#include <iostream>
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
   * @param json_str    JSON string.
   */
  void load(const std::string& json_str);

  /**
   * Compute the problem, write result into out stream, and write processing log into log stream.
   */
  void compute(std::ostream& out_stream, const io::OutputFormat out_format = io::OutputFormat::vtk,
               std::ostream& log_stream = std::cout);

  /**
   * Get the type.
   */
  std::string get_type() const;

private:
  class Impl;
  const std::unique_ptr<Impl> impl;
};

} // namespace plksim