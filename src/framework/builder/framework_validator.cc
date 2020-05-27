#include "framework/builder/framework_validator.h"

namespace bart {

namespace framework {

namespace builder {

void FrameworkValidator::Parse(const problem::ParametersI& to_parse) {
  needed_parts_are_present_ = {};
  needed_parts_are_present_.insert({FrameworkPart::ScatteringSourceUpdate,
                                    false});
  if (to_parse.IsEigenvalueProblem())
    needed_parts_are_present_.insert({FrameworkPart::FissionSourceUpdate,
                                      false});
  if (to_parse.HaveReflectiveBC() == true &&
      to_parse.TransportModel() == problem::EquationType::kSelfAdjointAngularFlux) {
    needed_parts_are_present_.insert({FrameworkPart::AngularSolutionStorage,
                                      false});
  }
}
} // namespace builder

} // namespace framework

} // namespace bart