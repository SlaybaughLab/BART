#ifndef BART_SRC_CONVERGENCE_REPORTER_TESTS_MPI_MOCK_H_
#define BART_SRC_CONVERGENCE_REPORTER_TESTS_MPI_MOCK_H_

#include "convergence/reporter/mpi_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace convergence {

namespace reporter {

class MpiMock : public MpiI {
 public:
  MOCK_METHOD1(Report, void(const bart::convergence::Status &to_report));
  MOCK_METHOD1(Report, void(const std::string &to_report));
};

} // namespace reporter

} // namespace convergence

} // namespace bart

#endif //BART_SRC_CONVERGENCE_REPORTER_TESTS_MPI_MOCK_H_
