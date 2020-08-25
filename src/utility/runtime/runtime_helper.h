#ifndef BART_SRC_UTILITY_RUNTIME_RUNTIME_HELPER_H_
#define BART_SRC_UTILITY_RUNTIME_RUNTIME_HELPER_H_

#include <string>


namespace bart {

namespace utility {

namespace runtime {

class RuntimeHelper {
 public:
  RuntimeHelper(std::string version);

  /// \brief Parses program arguments
  void ParseArguments(int argc, char** argv);

  /// \brief Returns the help message for the BART program
  std::string HelpMessage() const;

  /// \brief Returns the header for the BART program including version number
  std::string ProgramHeader() const;

  /// \brief Returns a string containing the version number
  std::string version() const { return version_; }

  /// \brief Returns bool indicating if the program should immediately terminate
  bool show_help() const { return show_help_; }

  /// \brief Returns bool indicating if a pause is required before running
  bool do_pause() const {return do_pause_; };

  /// \brief Returns a string containing the filename passed to bart
  std::string filename() const { return filename_; };

 private:
  bool show_help_{false};
  bool do_pause_{false};
  std::string filename_{};
  const std::string version_;
};

} // namespace runtime

} // namespace utility

} // namespace bart

#endif //BART_SRC_UTILITY_RUNTIME_RUNTIME_HELPER_H_
