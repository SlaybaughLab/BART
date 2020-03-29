#ifndef BART_SRC_UTILITY_REPORTER_MPI_I_H_
#define BART_SRC_UTILITY_REPORTER_MPI_I_H_

#include <string>

#include "convergence/status.h"

namespace bart {

namespace utility {

namespace reporter {

/*! \brief Reports convergence status and other text to stdout with MPI.
 *
 * Uses the dealii::ConditionalOStream object to safely output information when
 * using multiple processors.
 *
 */
class MpiI {
 public:
  virtual ~MpiI() = default;
  /*! \brief Report a given string.
   *
   * \param to_report string to output.
   */
  virtual void Report(const std::string &to_report) = 0;
};

} // namespace reporter

} // namespace utility

} // namespace bart

#endif // BART_SRC_UTILITY_REPORTER_MPI_I_H_