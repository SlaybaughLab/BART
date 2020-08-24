#ifndef BART_SRC_UTILITY_RUNTIME_RUNTIME_HELPER_H_
#define BART_SRC_UTILITY_RUNTIME_RUNTIME_HELPER_H_

#include <string>

namespace bart {

namespace utility {

namespace runtime {

class RuntimeHelper {
 public:
  RuntimeHelper(std::string version)
      : version_(version) {};

  /// \brief Returns the header for the BART program including version number
  std::string ProgramHeader() const;

  /// \brief Returns a string containing the version number
  std::string version() const { return version_; }
 private:
  const std::string version_;
};

} // namespace runtime

} // namespace utility

} // namespace bart

#endif //BART_SRC_UTILITY_RUNTIME_RUNTIME_HELPER_H_
