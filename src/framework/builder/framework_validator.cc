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

void FrameworkValidator::ReportValidation(
    utility::reporter::BasicReporterI& to_report) const {
  using Color = utility::reporter::Color;
  bool issue = false;
  to_report << "Validating framework components:\n";

  for (const auto part : parts_) {
    std::string description = framework_part_descriptions_.at(part);
    description[0] = std::toupper(description[0]);
    to_report << "\t";
    if (needed_parts_.count(part) > 0) {
      to_report << Color::Green << description << "\n";
    } else {
      to_report << Color::Yellow << description << " (Unneeded)\n";
      issue = true;
    }
  }

  for (const auto part : needed_parts_) {
    if (parts_.count(part) == 0) {
      std::string description = framework_part_descriptions_.at(part);
      description[0] = std::toupper(description[0]);
      to_report << "\t" << Color::Red << description << " (Missing)\n";
      issue = true;
    }
  }

  if (issue) {
    to_report.Report("Warning: one or more issues identified during "
                     "framework validation\n", Color::Yellow);
  }
}

} // namespace builder

} // namespace framework

} // namespace bart