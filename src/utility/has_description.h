#ifndef BART_SRC_UTILITY_IMPLEMENTATION_HAS_DESCRIPTION_H_
#define BART_SRC_UTILITY_IMPLEMENTATION_HAS_DESCRIPTION_H_

#include <string>
#include "utility/named_type.h"

namespace bart {

namespace utility {

class HasDescription {
 public:
  ~HasDescription() = default;
  void set_description(const std::string to_set) { description_ = to_set; };
  void is_default_implementation(const bool to_set) {
    is_default_implementation_ = to_set; }
  std::string description() const {
    if (is_default_implementation_)
      return "(default) " + description_;
    return description_; }
 private:
  bool is_default_implementation_ = false;
  std::string description_ = "";
};

} // namespace utility

} // namespace bart

#endif //BART_SRC_UTILITY_IMPLEMENTATION_HAS_DESCRIPTION_H_
