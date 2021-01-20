#ifndef BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_VALIDATOR_HPP_
#define BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_VALIDATOR_HPP_

#include "framework/builder/framework_validator_i.hpp"
#include "utility/colors.hpp"
#include "utility/uncopyable.h"

#include <set>

namespace bart::framework::builder {

namespace data_port {
struct ValidatorStatus;
using ValidatorStatusPort = instrumentation::Port<std::pair<std::string, utility::Color>, ValidatorStatus>;
} // namespace data_port

class FrameworkValidator : public FrameworkValidatorI,
                           public data_port::ValidatorStatusPort,
                           public utility::Uncopyable {
 public:
  auto AddPart(const FrameworkPart to_add) -> FrameworkValidator& override;
  auto Parse(const framework::FrameworkParameters) -> void override;
  [[deprecated]] auto Parse(const problem::ParametersI& to_parse) -> void override;
  auto ReportValidation() -> void override;

  [[nodiscard]] auto HasNeededParts() const -> bool  override { return needed_parts_.size() > 0; }
  [[nodiscard]] auto HasUnneededParts() const -> bool  override { return parts_.size() > needed_parts_.size(); }
  [[nodiscard]] auto NeededParts() const -> std::set<FrameworkPart>  override { return needed_parts_; }
  [[nodiscard]] auto Parts() const -> std::set<FrameworkPart>  override { return parts_; }
  [[nodiscard]] auto UnneededParts() const -> std::set<FrameworkPart> override;
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
