#include "framework/builder/framework_validator.h"

namespace bart {

namespace framework {

namespace builder {

void FrameworkValidator::Parse(const problem::ParametersI& to_parse) {
  needed_parts_ = {FrameworkPart::ScatteringSourceUpdate};
  if (to_parse.IsEigenvalueProblem())
    needed_parts_.insert(FrameworkPart::FissionSourceUpdate);
  if (to_parse.HaveReflectiveBC() == true &&
      to_parse.TransportModel() == problem::EquationType::kSelfAdjointAngularFlux) {
    needed_parts_.insert(FrameworkPart::AngularSolutionStorage);
  }
}
} // namespace builder

} // namespace framework

} // namespace bart