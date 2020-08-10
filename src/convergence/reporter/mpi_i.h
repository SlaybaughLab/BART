#ifndef BART_SRC_CONVERGENCE_REPORTER_MPI_I_H_
#define BART_SRC_CONVERGENCE_REPORTER_MPI_I_H_

#include <string>

#include "convergence/status.h"

namespace bart {

namespace convergence {

namespace reporter {

/*! \brief Reports convergence status and other text to stdout with MPI.
 *
 * Uses the dealii::ConditionalOStream object to safely output information when
 * using multiple processors.
 *
 */
class [[deprecated("Convergence reporters have been replaced by new "
                   "convergence status instrumentation")]] MpiI {
 public:
  virtual ~MpiI() = default;
  /*! \brief Report status of convergence.
   *
   * \param to_report convergence::Status holding the status of convergence.
   */
  virtual void Report(const bart::convergence::Status &to_report) = 0;

  /*! \brief Report a given string.
   *
   * \param to_report string to output.
   */
  virtual void Report(const std::string &to_report) = 0;
};

} // namespace reporter

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_REPORTER_MPI_I_H_