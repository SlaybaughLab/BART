#include "convergence/reporter/mpi_noisy.h"

namespace bart {

namespace convergence {

namespace reporter {

void MpiNoisy::Report(const bart::convergence::Status &to_report) {
  std::ostringstream status_report;

  status_report << "Iteration: " << to_report.iteration_number << "/"
                << to_report.max_iterations << "\tdelta: ";

  if (to_report.delta.has_value()) {
    status_report << to_report.delta.value();
  } else {
    status_report << "N/A";
  }

  status_report << "\tidx: ";

  if (to_report.failed_index.has_value()) {
    status_report << to_report.failed_index.value();
  } else {
    status_report << "N/A";
  }

  status_report << std::endl;

  *pout_ptr_ << status_report.str();
}

} // namespace reporter

} // namespace convergence

} // namespace bart