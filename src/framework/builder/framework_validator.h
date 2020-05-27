#ifndef BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_VALIDATOR_H_
#define BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_VALIDATOR_H_

#include "problem/parameters_i.h"

#include <map>

namespace bart {

namespace framework {

namespace builder {

enum class FrameworkPart {
  ScatteringSourceUpdate = 0,
  FissionSourceUpdate = 1,
  AngularSolutionStorage = 2
};

class FrameworkValidator {
 public:
  void Parse(const problem::ParametersI& to_parse);

  bool HasNeededParts() const {
    return needed_parts_are_present_.size() > 0; }

  std::vector<FrameworkPart> NeededParts() const {
    std::vector<FrameworkPart> return_vector;
    for (const auto part_pair : needed_parts_are_present_)
      return_vector.push_back(part_pair.first);
    return return_vector; }
 private:

  std::map<FrameworkPart, bool> needed_parts_are_present_{};

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
