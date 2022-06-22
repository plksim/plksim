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
   * @param jsonStr    JSON string.
   */
  void load(const std::string& jsonStr);

  /**
   * Compute the problem, write result into out stream.
   */
  void compute(std::ostream& out, const io::OutputFormat outFormat = io::OutputFormat::vtk);

  /**
   * Get the type.
   */
  std::string getType() const;

private:
  class Impl;
  const std::unique_ptr<Impl> impl;
};

} // namespace plksim