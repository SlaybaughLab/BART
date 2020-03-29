#ifndef BART_SRC_UTILITY_REPORTER_MPI_H_
#define BART_SRC_UTILITY_REPORTER_MPI_H_

#include <memory>

#include <deal.II/base/conditional_ostream.h>
#include <ostream>

#include "utility/reporter/mpi_i.h"
#include "utility/uncopyable.h"

namespace bart {

namespace utility {

namespace reporter {

/*! \brief Noisy reporter, will report each time Report is called for
 * convergence, or a passed string.
 */

class Mpi : public MpiI, private Uncopyable {
 public:
  Mpi(std::unique_ptr<dealii::ConditionalOStream> pout_ptr)
      : pout_ptr_(std::move(pout_ptr)) {};
  ~Mpi() = default;

  void Report(const std::string &to_report) override {
    *pout_ptr_ << to_report;
  }

 private:
  std::unique_ptr<dealii::ConditionalOStream> pout_ptr_;


};
} // namespace reporter

} // namespace utility

} // namespace bart

#endif // BART_SRC_UTILITY_REPORTER_MPI_H_