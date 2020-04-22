#ifndef BART_SRC_UTILITY_IMPLEMENTATION_HAS_DESCRIPTION_H_
#define BART_SRC_UTILITY_IMPLEMENTATION_HAS_DESCRIPTION_H_

#include <string>
#include "utility/named_type.h"

namespace bart {

namespace utility {

using DefaultImplementation = NamedType<bool, struct DefaultImplementationStruct>;

class HasDescription {
 public:
  ~HasDescription() = default;

  inline HasDescription& set_description(const std::string to_set,
                                         const DefaultImplementation is_default) {
    description_ = to_set;
    is_default_implementation_ = is_default.get();
    return *this;
  }

  inline HasDescription& set_description(const std::string to_set) {
    return set_description(to_set,
                           DefaultImplementation(is_default_implementation_));
  }

  inline HasDescription& set_is_default_implementation(const bool to_set) {
    is_default_implementation_ = to_set;
    return *this;
  }

  inline std::string description() const {
    if (is_default_implementation_)
      return "(default) " + description_;
    return description_;
  }

 private:
  bool is_default_implementation_ = false;
  std::string description_ = "";
};

} // namespace utility

} // namespace bart

#endif //BART_SRC_UTILITY_IMPLEMENTATION_HAS_DESCRIPTION_H_
