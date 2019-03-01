#ifndef BART_SRC_POST_REPORTER_MPI_NOISY_H_
#define BART_SRC_POST_REPORTER_MPI_NOISY_H_

#include "post/reporter/mpi_i.h"

namespace bart {

namespace post {

namespace reporter {

/*! \brief Noisy reporter, will report each time Report is called for
 * convergence, or a passed string.
 */

class MpiNoisy : public MpiI {
 public:
  MpiNoisy() = default;
  ~MpiNoisy() = default;

};

} // namespace reporter

} // namespace post

} // namespace bart

#endif // BART_SRC_POST_REPORTER_MPI_NOISY_H_