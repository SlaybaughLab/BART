#ifndef BART_SRC_UTILITY_HAS_DEPENDENCIES_H_
#define BART_SRC_UTILITY_HAS_DEPENDENCIES_H_

#include <deal.II/base/exceptions.h>
#include <string_view>

namespace bart {

namespace utility {

class HasDependencies {
 public:
  template <typename T>
  void AssertPointerNotNull(T* pointer,
                            const std::string& dependency_name,
                            const std::string& calling_location) {
    std::string error_string{"Error in " + calling_location +
        ", dependency pointer: " + dependency_name + " is null."};
    AssertThrow(pointer != nullptr, dealii::ExcMessage(error_string));
  }
};

} // namespace utility

} // namespace bart

#endif //BART_SRC_UTILITY_HAS_DEPENDENCIES_H_
