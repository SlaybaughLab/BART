#include "mpi_noisy.h"

namespace bart {

namespace post {

namespace reporter {

void MpiNoisy::Report(const bart::convergence::Status &to_report) {
  std::ostringstream status_report;
  status_report << "Iteration: " << to_report.iteration_number << "/"
                << to_report.max_iterations << "\tdelta: " << to_report.delta.value()
                << "\tidx: " << to_report.failed_index.value() << std::endl;
  *pout_ptr_ << status_report.str();
}

} // namespace reporter

} // namespace post

} // namespace bart