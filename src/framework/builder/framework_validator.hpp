#ifndef BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_VALIDATOR_HPP_
#define BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_VALIDATOR_HPP_

#include "problem/parameters_i.h"
#include "instrumentation/port.hpp"
#include "framework/framework_parameters.hpp"
#include "instrumentation/port.h"
#include "problem/parameters_i.h"
#include "utility/colors.hpp"
#include "utility/uncopyable.h"

#include <set>

namespace bart::framework::builder {

enum class FrameworkPart {
  ScatteringSourceUpdate = 0,
  FissionSourceUpdate = 1,
  AngularSolutionStorage = 2
};

namespace data_port {
struct ValidatorStatus;
using ValidatorStatusPort = instrumentation::Port<std::pair<std::string, utility::Color>, ValidatorStatus>;
} // namespace data_port

class FrameworkValidator : public data_port::ValidatorStatusPort, public utility::Uncopyable {
 public:
  auto AddPart(const FrameworkPart to_add) -> FrameworkValidator&;
  auto Parse(const framework::FrameworkParameters) -> void;
  auto Parse(const problem::ParametersI& to_parse) -> void;
  auto ReportValidation() -> void;

  [[nodiscard]] auto HasNeededParts() const -> bool { return needed_parts_.size() > 0; }
  [[nodiscard]] auto HasUnneededParts() const -> bool { return parts_.size() > needed_parts_.size(); }
  [[nodiscard]] auto NeededParts() const -> std::set<FrameworkPart> { return needed_parts_; }
  [[nodiscard]] auto Parts() const -> std::set<FrameworkPart> { return parts_; }
  [[nodiscard]] auto UnneededParts() const -> std::set<FrameworkPart>;
 private:
  std::set<FrameworkPart> needed_parts_{};
  std::set<FrameworkPart> parts_{};

  std::map<FrameworkPart, std::string> framework_part_descriptions_{
      {FrameworkPart::ScatteringSourceUpdate, "scattering source updater"},
      {FrameworkPart::FissionSourceUpdate, "fission source updater"},
      {FrameworkPart::AngularSolutionStorage, "angular solution storage"}
  };
};

} // namespace bart::framework::builder

#endif //BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_VALIDATOR_HPP_
