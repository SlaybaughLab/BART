#ifndef BART_SRC_CONVERGENCE_REPORTER_TESTS_MPI_MOCK_H_
#define BART_SRC_CONVERGENCE_REPORTER_TESTS_MPI_MOCK_H_

#include "utility/reporter/basic_reporter_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace convergence {

namespace reporter {

class MpiMock : public BasicReporterI {
 public:
  MOCK_METHOD(void, Report, (const std::string &to_report), (override));
  MOCK_METHOD(void, Report, (const std::string &to_report, Color), (override));
};

} // namespace reporter

} // namespace convergence

} // namespace bart

#endif //BART_SRC_CONVERGENCE_REPORTER_TESTS_MPI_MOCK_H_
