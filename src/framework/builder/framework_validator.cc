#include "framework/builder/framework_validator.h"

namespace bart {

namespace framework {

namespace builder {

void FrameworkValidator::Parse(const problem::ParametersI& to_parse) {
  needed_parts_are_present_.insert({FrameworkPart::ScatteringSourceUpdate,
                                    false});
  if (to_parse.IsEigenvalueProblem())
    needed_parts_are_present_.insert({FrameworkPart::FissionSourceUpdate, false});
}
} // namespace builder

} // namespace framework

} // namespace bart