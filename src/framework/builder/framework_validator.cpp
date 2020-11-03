#include "framework/builder/framework_validator.hpp"

#include <algorithm>

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

  const auto reflective_boundary = to_parse.ReflectiveBoundary();
  bool has_reflective = std::any_of(
      reflective_boundary.begin(),
      reflective_boundary.end(),
      [](std::pair<problem::Boundary, bool> pair){ return pair.second; });

  if (has_reflective &&
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

void FrameworkValidator::ReportValidation() {
  using Color = utility::Color;
  bool issue = false;
  data_port::ValidatorStatusPort::Expose({"Validating framework components:\n",
                                          Color::kReset});
  for (const auto part : parts_) {
    std::string description = framework_part_descriptions_.at(part);
    description[0] = std::toupper(description[0]);
    Expose({"\t", Color::kReset});
    if (needed_parts_.count(part) > 0) {
      Expose({description + "\n", Color::kGreen});
    } else {
      Expose({description + " (Unneeded)\n", Color::kYellow});
      issue = true;
    }
  }

  for (const auto part : needed_parts_) {
    if (parts_.count(part) == 0) {
      std::string description = framework_part_descriptions_.at(part);
      description[0] = std::toupper(description[0]);
      Expose({"\t" + description + "(Missing)\n", Color::kRed});
      issue = true;
    }
  }

  if (issue) {
    Expose({"Warning: one or more issues identified during "
            "framework validation\n", Color::kYellow});
  }
}

} // namespace builder

} // namespace framework

} // namespace bart