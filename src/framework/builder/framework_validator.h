#ifndef BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_VALIDATOR_H_
#define BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_VALIDATOR_H_

#include "problem/parameters_i.h"
#include "instrumentation/port.h"
#include "utility/colors.hpp"

#include <set>

namespace bart {

namespace framework {

namespace builder {

enum class FrameworkPart {
  ScatteringSourceUpdate = 0,
  FissionSourceUpdate = 1,
  AngularSolutionStorage = 2
};

namespace data_port {
struct ValidatorStatus;
using ValidatorStatusPort = instrumentation::Port<std::pair<std::string, utility::Color>, ValidatorStatus>;
} // namespace data_port

class FrameworkValidator : public data_port::ValidatorStatusPort {
 public:
  FrameworkValidator& AddPart(const FrameworkPart to_add);

  bool HasNeededParts() const {
    return needed_parts_.size() > 0; }

  bool HasUnneededParts() const {
    return parts_.size() > needed_parts_.size(); }
  std::set<FrameworkPart> NeededParts() const {
    return needed_parts_; }
  void Parse(const problem::ParametersI& to_parse);
  std::set<FrameworkPart> Parts() const {
    return parts_; }
  void ReportValidation();
  std::set<FrameworkPart> UnneededParts() const;
 private:

  std::set<FrameworkPart> needed_parts_{};
  std::set<FrameworkPart> parts_{};

  std::map<FrameworkPart, std::string> framework_part_descriptions_{
      {FrameworkPart::ScatteringSourceUpdate, "scattering source updater"},
      {FrameworkPart::FissionSourceUpdate, "fission source updater"},
      {FrameworkPart::AngularSolutionStorage, "angular solution storage"}
  };

};

} // namespace builder

} // namespace framework

} // namespace bart

#endif //BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_VALIDATOR_H_
