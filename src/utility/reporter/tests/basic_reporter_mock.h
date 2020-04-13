#ifndef BART_SRC_UTILITY_REPORTER_TESTS_MPI_MOCK_H_
#define BART_SRC_UTILITY_REPORTER_TESTS_MPI_MOCK_H_

#include "utility/reporter/basic_reporter_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace utility {

namespace reporter {

class BasicReporterMock : public BasicReporterI {
 public:
  MOCK_METHOD(void, Report, (const std::string &to_report), (override));
  MOCK_METHOD(void, Report, (const std::string &to_report, Color), (override));
  MOCK_METHOD(BasicReporterMock&, Instream, (const std::string&));
  MOCK_METHOD(BasicReporterMock&, Instream, (const Color));
  BasicReporterMock& operator<<(const std::string& to_report) override {
    return Instream(to_report); }
  BasicReporterMock& operator<<(const Color color) override {
    return Instream(color); }
};

} // namespace reporter

} // namespace utility

} // namespace bart

#endif //BART_SRC_UTILITY_REPORTER_TESTS_MPI_MOCK_H_
