#ifndef BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_VALIDATOR_H_
#define BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_VALIDATOR_H_

#include "problem/parameters_i.h"

#include <set>

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
    return needed_parts_.size() > 0; }

  bool HasUnneededParts() const {
    return parts_.size() > needed_parts_.size(); }

  std::set<FrameworkPart> NeededParts() const {
    return needed_parts_; }

  std::set<FrameworkPart> UnneededParts() const {
    std::set<FrameworkPart> return_set;
    if (HasUnneededParts()) {
      for (const auto part : needed_parts_) {
        if (parts_.count(part) > 0)
          return_set.insert(part);
      }
    }
    return return_set; }
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
