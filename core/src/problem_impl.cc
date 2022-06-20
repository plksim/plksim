#include "problem_impl.hh"

#include <nlohmann/json.hpp>

#include <exception>
#include <iostream>

using json = nlohmann::json;

namespace plksim {

Problem::Impl::Impl(){};
Problem::Impl::~Impl(){};

void Problem::Impl::load(const std::string& jsonStr) {
  auto obj = json::parse(jsonStr);

  if (!obj["type"].is_string()) {
    throw std::runtime_error("wrong format: type is not defined");
  } else {
    mType = obj["type"].get<std::string>();
  }
};

std::string Problem::Impl::getType() {
  return mType;
};

} // namespace plksim