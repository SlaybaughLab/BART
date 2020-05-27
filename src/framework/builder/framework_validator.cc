#include "framework/builder/framework_validator.h"

namespace bart {

namespace framework {

namespace builder {

FrameworkValidator &FrameworkValidator::AddPart(FrameworkPart to_add) {
  parts_.insert(to_add);
  return *this;
}

void FrameworkValidator::Parse(const problem::ParametersI& to_parse) {
  needed_parts_ = {FrameworkPart::ScatteringSourceUpdate};
  if (to_parse.IsEigenvalueProblem())
    needed_parts_.insert(FrameworkPart::FissionSourceUpdate);
  if (to_parse.HaveReflectiveBC() == true &&
      to_parse.TransportModel() == problem::EquationType::kSelfAdjointAngularFlux) {
    needed_parts_.insert(FrameworkPart::AngularSolutionStorage);
  }
}

std::set<FrameworkPart> FrameworkValidator::UnneededParts() const {
  std::set<FrameworkPart> return_set;
  if (HasUnneededParts()) {
    for (const auto part : parts_) {
      if (needed_parts_.count(part) == 0)
        return_set.insert(part);
    }
  }
  return return_set; }

} // namespace builder

} // namespace framework

} // namespace bart