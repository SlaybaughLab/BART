#ifndef BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_VALIDATOR_I_HPP_
#define BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_VALIDATOR_I_HPP_

#include "problem/parameters_i.h"
#include "instrumentation/port.hpp"
#include "framework/framework_parameters.hpp"
#include "problem/parameters_i.h"

namespace bart::framework::builder {

enum class FrameworkPart {
  ScatteringSourceUpdate = 0,
  FissionSourceUpdate = 1,
  AngularSolutionStorage = 2
};

class FrameworkValidatorI {
 public:
  virtual auto AddPart(const FrameworkPart to_add) -> FrameworkValidatorI& = 0;
  virtual auto Parse(const framework::FrameworkParameters) -> void = 0;
  virtual auto Parse(const problem::ParametersI& to_parse) -> void = 0;
  virtual auto ReportValidation() -> void = 0;

  virtual auto HasNeededParts() const -> bool = 0;
  virtual auto HasUnneededParts() const -> bool = 0;
  virtual auto NeededParts() const -> std::set<FrameworkPart> = 0;
  virtual auto Parts() const -> std::set<FrameworkPart> = 0;
  virtual auto UnneededParts() const -> std::set<FrameworkPart> = 0;
};

} // namespace bart::framework::builder

#endif //BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_VALIDATOR_I_HPP_
