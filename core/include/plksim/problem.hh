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
   * Set the stream to output log.
   */
  void set_log_stream(std::ostream& stream);

  /**
   * Set the stream to output grid.
   */
  void set_output_grid_stream(std::ostream& stream,
                              const io::OutputFormat format = io::OutputFormat::vtu);

  /**
   * Set the stream to output solution.
   */
  void set_output_solution_stream(std::ostream& stream,
                                  const io::OutputFormat format = io::OutputFormat::vtu);

  /**
   * Compute the problem, and write result and log into streams.
   */
  void compute();

  /**
   * Get the type.
   */
  std::string get_type() const;

private:
  class Impl;
  const std::unique_ptr<Impl> impl;
};

} // namespace plksim