#ifndef BART_SRC_UTILITY_REPORTER_BASIC_REPORTER_I_H_
#define BART_SRC_UTILITY_REPORTER_BASIC_REPORTER_I_H_

#include <string>

#include "utility/reporter/colors.h"

namespace bart {

namespace utility {

namespace reporter {

/*! \brief Reports convergence status and other text to stdout with MPI.
 *
 * Uses the dealii::ConditionalOStream object to safely output information when
 * using multiple processors.
 *
 */
class BasicReporterI {
 public:
  virtual ~BasicReporterI() = default;
  /*! \brief Report a given string.
   *
   * \param to_report string to output.
   */
  virtual void Report(const std::string &to_report) = 0;
  virtual void Report(const std::string &to_report, Color) = 0;
  virtual BasicReporterI& operator<<(const std::string&) = 0;
  virtual BasicReporterI& operator<<(const Color) = 0;
};

} // namespace reporter

} // namespace utility

} // namespace bart

#endif //BART_SRC_UTILITY_REPORTER_BASIC_REPORTER_I_H_